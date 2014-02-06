//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_TRA.C,v $
//
// Purpose        : Implement lossless transmission line
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Component Information and Models
//
// Creation Date  : 06/14/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.147.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_TRA.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_BreakPoint.h>

#include <N_UTL_Functors.h>

namespace Xyce {
namespace Device {
template<>
ParametricData<TRA::Instance>::ParametricData()
{
    // Set up configuration constants:
    setNumNodes(4);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(0);

    // Set up double precision variables:
    addPar ("Z0", 0.0, false, ParameterType::NO_DEP,
      &TRA::Instance::Z0,
      NULL,U_OHM,CAT_NONE,"Characteristic Impedance");

    addPar ("ZO", 0.0, false, ParameterType::NO_DEP,
      &TRA::Instance::ZO,
      NULL,U_OHM,CAT_NONE,"Characteristic Impedance");

    addPar ("TD", 0.0, false, ParameterType::NO_DEP,
      &TRA::Instance::td,
      NULL,U_SECOND,CAT_NONE,"Time delay");

    addPar ("F", 0.0, false, ParameterType::NO_DEP,
      &TRA::Instance::freq,
      NULL,U_HZ,CAT_NONE,"Frequency");

    addPar ("NL", 0.0, false, ParameterType::NO_DEP,
      &TRA::Instance::NL,
      NULL,U_NONE,CAT_NONE,"Length in wavelengths");
}

template<>
ParametricData<TRA::Model>::ParametricData()
{
}

namespace TRA {

vector< vector<int> > Instance::jacStamp;



ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 6/15/01
//-----------------------------------------------------------------------------

Instance::Instance(InstanceBlock & IB,
                                     Model & Miter,
                                     MatrixLoadData & mlData1,
                                     SolverState &ss1,
                                     ExternData  &ed1,
                                     DeviceOptions & do1)

  : DeviceInstance(IB,mlData1,ss1,ed1,do1),
    model_(Miter),
    Z0(0.0),
    G0(0.0),
    td(0.0),
    freq(0.0),
    NL(0.25),
    newtonIterOld(0),
    timeOld(-1.0),
    li_Pos1(-1),
    li_Neg1(-1),
    li_Int1(-1),
    li_Ibr1(-1),
    li_Pos2(-1),
    li_Neg2(-1),
    li_Int2(-1),
    li_Ibr2(-1),
    li_store_dev_i1(-1),
    li_store_dev_i2(-1),
    APos1EquPos1NodeOffset(-1),
    APos1EquInt1NodeOffset(-1),
    AInt1EquPos1NodeOffset(-1),
    AInt1EquInt1NodeOffset(-1),
    AInt1EquIbr1NodeOffset(-1),
    ANeg1EquIbr1NodeOffset(-1),
    AIbr1EquInt1NodeOffset(-1),
    AIbr1EquNeg1NodeOffset(-1),
    APos2EquPos2NodeOffset(-1),
    APos2EquInt2NodeOffset(-1),
    AInt2EquPos2NodeOffset(-1),
    AInt2EquInt2NodeOffset(-1),
    AInt2EquIbr2NodeOffset(-1),
    ANeg2EquIbr2NodeOffset(-1),
    AIbr2EquInt2NodeOffset(-1),
    AIbr2EquNeg2NodeOffset(-1),
    AIbr1EquPos2NodeOffset(-1),
    AIbr1EquNeg2NodeOffset(-1),
    AIbr1EquIbr2NodeOffset(-1),
    AIbr2EquPos1NodeOffset(-1),
    AIbr2EquNeg1NodeOffset(-1),
    AIbr2EquIbr1NodeOffset(-1),
    first_BP_call_done(false),
    last_t(0.0),
    v1(0.0),
    v2(0.0),
    newBreakPoint(false),
    newBreakPointTime(0.0)
{
  numIntVars   = 4;
  numExtVars   = 4;
  numStateVars = 0;
  numLeadCurrentStoreVars = 2;  // lead currents i1 and i2

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 1;
  devConMap[2] = 2;
  devConMap[3] = 2;

  if( jacStamp.empty() )
  {
    jacStamp.resize(8);
    jacStamp[0].resize(2);
    jacStamp[0][0]=0;
    jacStamp[0][1]=4;
    jacStamp[1].resize(1);
    jacStamp[1][0]=5;
    jacStamp[2].resize(2);
    jacStamp[2][0]=2;
    jacStamp[2][1]=6;
    jacStamp[3].resize(1);
    jacStamp[3][0]=7;
    jacStamp[4].resize(3);
    jacStamp[4][0]=0;
    jacStamp[4][1]=4;
    jacStamp[4][2]=5;
    jacStamp[5].resize(5);
    jacStamp[5][0]=1;
    jacStamp[5][1]=2;
    jacStamp[5][2]=3;
    jacStamp[5][3]=4;
    jacStamp[5][4]=7;
    jacStamp[6].resize(3);
    jacStamp[6][0]=2;
    jacStamp[6][1]=6;
    jacStamp[6][2]=7;
    jacStamp[7].resize(5);
    jacStamp[7][0]=0;
    jacStamp[7][1]=1;
    jacStamp[7][2]=3;
    jacStamp[7][3]=5;
    jacStamp[7][4]=6;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  if (!given("Z0"))
  {
    if (given("ZO"))
      Z0 = ZO;
    else
    {
      string msg = "Z0 not given.";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr:: report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
  }
  if (Z0>0)
  {
    G0 = 1.0/Z0;
  }
  else
  {
     string msg = "Invalid (zero or negative) impedance given.";
     std::ostringstream oss;
     oss << "Error in " << netlistLocation() << "\n" << msg;
     N_ERH_ErrorMgr:: report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  // Must give either TD or F.
  if (!given("TD") && !given("F"))
  {
    string msg = "Neither time delay (TD) nor frequency (F) given.";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr:: report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }
  if (given("TD") && given("F"))
  {
    string msg = "Both time delay (TD) and frequency (F) given.  Pick one.";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr:: report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  if (!given("TD") && freq > 0)
  {
    td = NL/freq;
  }
  else if (!given("TD"))
  {
    string msg = "Zero or negative frequency given.";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr:: report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }
  if (td == 0)
  {
    string msg = "Zero time delay.";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr:: report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }


  processParams();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    std::cout << " Z0 = " << Z0 << std::endl;
    std::cout << " td = " << td << std::endl;
    std::cout << " freq = " << freq << std::endl;
    std::cout << " NL = " << NL << std::endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 7/02/04
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       : function for registering, and setting up, local ID's.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                                      const vector<int> & extLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    std::cout << std::endl << dashedline << std::endl;
    std::cout << "In Instance::registerLIDs\n\n";
    std::cout << "name             = " << getName() << std::endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    std::cout << "number of internal variables: " << numInt << std::endl;
    std::cout << "number of external variables: " << numExt << std::endl;
  }
#endif

  if (numInt != numIntVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.

  li_Pos1 = extLIDVec[0];
  li_Neg1 = extLIDVec[1];
  li_Pos2 = extLIDVec[2];
  li_Neg2 = extLIDVec[3];

  // Now do the internal variables

  li_Int1 = intLIDVec[0];
  li_Ibr1 = intLIDVec[1];
  li_Int2 = intLIDVec[2];
  li_Ibr2 = intLIDVec[3];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    std::cout << " VARIABLE Indicies " << std::endl;
    std::cout << "li_Pos1 = " << li_Pos1 << std::endl;
    std::cout << "li_Neg1 = " << li_Neg1 << std::endl;
    std::cout << "li_Int1 = " << li_Int1 << std::endl;
    std::cout << "li_Ibr1 = " << li_Ibr1 << std::endl;
    std::cout << "li_Pos2 = " << li_Pos2 << std::endl;
    std::cout << "li_Neg2 = " << li_Neg2 << std::endl;
    std::cout << "li_Int2 = " << li_Int2 << std::endl;
    std::cout << "li_Ibr2 = " << li_Ibr2 << std::endl;

    std::cout << "\nEnd of Instance::register LIDs\n";
    std::cout << dashedline << std::endl;
  }
#endif
}
//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {

    // set up internal name map
    string tmpstr;
    tmpstr = getName()+"_int1"; spiceInternalName (tmpstr);
    intNameMap[ li_Int1 ] = tmpstr;

    tmpstr = getName()+"_int2"; spiceInternalName (tmpstr);
    intNameMap[ li_Int2 ] = tmpstr;

    tmpstr = getName()+"_i1"; spiceInternalName (tmpstr);
    intNameMap[ li_Ibr1 ] = tmpstr;

    tmpstr = getName()+"_i2"; spiceInternalName (tmpstr);
    intNameMap[ li_Ibr2 ] = tmpstr;
  }

  return intNameMap;
}



//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef)
{
  string msg;

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

}

//-----------------------------------------------------------------------------
// Function      : N_DEV_TRAInstance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/05/2013
//-----------------------------------------------------------------------------
void N_DEV_TRAInstance::registerStoreLIDs(const vector<int> & stoLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSto = stoLIDVecRef.size();

  if (numSto != getNumStoreVars())
  {
    msg = "N_DEV_TRAInstance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
  if( loadLeadCurrent )
  {
    li_store_dev_i1 = stoLIDVecRef[0];
    li_store_dev_i2 = stoLIDVecRef[1];
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_TRAInstance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/05/2013
//-----------------------------------------------------------------------------
map<int,string> & N_DEV_TRAInstance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if( loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    string modName(getName());
    spiceInternalName(modName);
    string tmpstr;
    tmpstr = modName+":DEV_I1";
    storeNameMap[ li_store_dev_i1 ] = tmpstr;
    tmpstr = modName+":DEV_I2";
    storeNameMap[ li_store_dev_i2 ] = tmpstr;
  }

  return storeNameMap;
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_TRAInstance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  APos1EquPos1NodeOffset = jacLIDVec[0][0];
  APos1EquInt1NodeOffset = jacLIDVec[0][1];

  ANeg1EquIbr1NodeOffset = jacLIDVec[1][0];

  APos2EquPos2NodeOffset = jacLIDVec[2][0];
  APos2EquInt2NodeOffset = jacLIDVec[2][1];

  ANeg2EquIbr2NodeOffset = jacLIDVec[3][0];

  AInt1EquPos1NodeOffset = jacLIDVec[4][0];
  AInt1EquInt1NodeOffset = jacLIDVec[4][1];
  AInt1EquIbr1NodeOffset = jacLIDVec[4][2];

  AIbr1EquNeg1NodeOffset = jacLIDVec[5][0];
  AIbr1EquPos2NodeOffset = jacLIDVec[5][1];
  AIbr1EquNeg2NodeOffset = jacLIDVec[5][2];
  AIbr1EquInt1NodeOffset = jacLIDVec[5][3];
  AIbr1EquIbr2NodeOffset = jacLIDVec[5][4];

  AInt2EquPos2NodeOffset = jacLIDVec[6][0];
  AInt2EquInt2NodeOffset = jacLIDVec[6][1];
  AInt2EquIbr2NodeOffset = jacLIDVec[6][2];

  AIbr2EquPos1NodeOffset = jacLIDVec[7][0];
  AIbr2EquNeg1NodeOffset = jacLIDVec[7][1];
  AIbr2EquNeg2NodeOffset = jacLIDVec[7][2];
  AIbr2EquIbr1NodeOffset = jacLIDVec[7][3];
  AIbr2EquInt2NodeOffset = jacLIDVec[7][4];

}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 TRA instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;

  double coef_pos1;
  double coef_neg1;
  double coef_int1;
  double coef_ibr1;
  double coef_pos2;
  double coef_neg2;
  double coef_int2;
  double coef_ibr2;

  double * fVec = extData.daeFVectorRawPtr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << dashedline2 << std::endl;
    std::cout << "  Instance::loadDAEFVector" << std::endl;
    std::cout << "  name = " << getName() <<std::endl;
  }
#endif

  // Most of the work has already been done by uIVB.
  coef_pos1 = (Vpos1-Vint1)*G0;
  coef_neg1 = -Ibr1;
  coef_int1 = -(Vpos1-Vint1)*G0+Ibr1;
  coef_ibr1 = ((Vint1-Vneg1)-v1);
  coef_pos2 = (Vpos2-Vint2)*G0;
  coef_neg2 = -Ibr2;
  coef_int2 = -(Vpos2-Vint2)*G0+Ibr2;
  coef_ibr2 = ((Vint2-Vneg2)-v2);


  fVec[li_Pos1] += coef_pos1;
  fVec[li_Neg1] += coef_neg1;
  fVec[li_Int1] += coef_int1;
  fVec[li_Ibr1] += coef_ibr1;
  fVec[li_Pos2] += coef_pos2;
  fVec[li_Neg2] += coef_neg2;
  fVec[li_Int2] += coef_int2;
  fVec[li_Ibr2] += coef_ibr2;

  if( loadLeadCurrent )
  {
    double * stoVec = extData.nextStoVectorRawPtr;
    stoVec[li_store_dev_i1] = Ibr1;
    stoVec[li_store_dev_i2] = -Ibr2;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single device instance.
//
// Special Notes : The F-vector is an algebraic constaint.
//
//                 The special notes below are those that were taken from
//                 the old loadAnalyticJacobian header.  The matrix
//                 it describes is the full jacobian matrix:
//----------------------------------------------------------------------
// Special Notes : This is based on there being two two-node ports
//                 and a model of the following sort:
//
//  Pos1  o-----+          +------o Pos2
//              |          |
//              \          \
//              /          /
//              \Z0        \ Z0
//              /          /
//              |          |
//              oInt1      o Int2
//              |          |
//           ++++++      ++++++
//           | V1 |      | V2 |
//           ------      ------
//              |          |
//   Neg1 o-----+          +------o Neg2
//
// There are also two branch currents, Ibr1 and Ibr2 for left and right
// sides as well.
//
//  The matrix for this ends up being:
//            V_Pos1  V_Neg1  V_Int1  Ibr1  V_Pos2  V_Neg2  V_Int2  V_Ibr2
//            ------------------------------------------------------------
//  KCL Pos1    a               b
//  KCL Neg1                           c
//  KCL Int1    d               e      f
//  KCL Ibr1            g       h             i       j               k
//  KCL Pos2                                  l                m
//  KCL Neg2                                                          n
//  KCL Int2                                  o                p      q
//  KCL Ibr2    r       s              t              u        v
//
//  When doing time integration, i,j,k,r,s and t are zero, those dependences
//  are time-delayed, i.e. the equations for the output depend on time-delayed
//  of the input.  For DC calculations, i,j,k,r,s and t are non-zero.
//
// The right hand sides are:
// Pos1:  (V_int1-V_Pos1)*G0            Pos2: (V_Int2-V_Pos2)*G0
// Neg1:  -Ibr1                         Neg2: -Ibr2
// Int1:  (V_Pos1-V_Int1)*G0+Ibr1       Int2: (V_Pos2-V_Int2)*G0+Ibr2
// Ibr1:  (V_Int1-V_Neg1)-V1            Ibr2: (V_Int2-V_Neg2)-V2
//
// For transient operation, v1 and v2 depend on values of voltage and
//  current at delayed time at the opposite port:
//   V1 = DeltaV2(t-td)+Z0*Ibr2(t-td)
//   V2 = DeltaV1(t-td)+Z0*Ibr1(t-td)
//
// For DC operation V1=VPos2-Vneg2+Ibr2*Z0, V2 = Vpos1-Vneg1+Ibr1*Z0
//--------------------------------------------------------------------
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    std::cout << dashedline2 << std::endl;
    std::cout << "  name             = " << getName() << std::endl;
  }
#endif


  dFdx[li_Pos1][APos1EquPos1NodeOffset] += G0;

  dFdx[li_Pos1][APos1EquInt1NodeOffset] -= G0;


  dFdx[li_Int1][AInt1EquPos1NodeOffset] -= G0;

  dFdx[li_Int1][AInt1EquInt1NodeOffset] += G0;

  dFdx[li_Int1][AInt1EquIbr1NodeOffset] += 1.0;


  dFdx[li_Neg1][ANeg1EquIbr1NodeOffset] -= 1.0;


  dFdx[li_Ibr1][AIbr1EquInt1NodeOffset] += 1.0;

  dFdx[li_Ibr1][AIbr1EquNeg1NodeOffset] -= 1.0;
  if( DCMODE )
  {

    dFdx[li_Ibr1][AIbr1EquPos2NodeOffset] -= 1.0;

    dFdx[li_Ibr1][AIbr1EquNeg2NodeOffset] += 1.0;

    dFdx[li_Ibr1][AIbr1EquIbr2NodeOffset] -= Z0;
  }


  dFdx[li_Pos2][APos2EquPos2NodeOffset] += G0;

  dFdx[li_Pos2][APos2EquInt2NodeOffset] -= G0;


  dFdx[li_Int2][AInt2EquPos2NodeOffset] -= G0;

  dFdx[li_Int2][AInt2EquInt2NodeOffset] += G0;

  dFdx[li_Int2][AInt2EquIbr2NodeOffset] += 1.0;


  dFdx[li_Neg2][ANeg2EquIbr2NodeOffset] -= 1.0;


  dFdx[li_Ibr2][AIbr2EquInt2NodeOffset] += 1.0;

  dFdx[li_Ibr2][AIbr2EquNeg2NodeOffset] -= 1.0;
  if( DCMODE )
  {

    dFdx[li_Ibr2][AIbr2EquPos1NodeOffset] -= 1.0;

    dFdx[li_Ibr2][AIbr2EquNeg1NodeOffset] += 1.0;

    dFdx[li_Ibr2][AIbr2EquIbr1NodeOffset] -= Z0;
  }


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    std::cout << dashedline2 << std::endl;
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       : update primary state for one TRA instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << std::endl << dashedline2 << std::endl;
    std::cout << "In TRA::updatePrimaryState\n";
    std::cout << " last_t is " << last_t << std::endl;
    std::cout << " v1 is " << v1 << std::endl;
    std::cout << " v2 is " << v2 << std::endl;
  }
#endif

  return updateIntermediateVars ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one TRA instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;    // the current guess

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << std::endl << dashedline2 << std::endl;
    std::cout << "  In ::updateIntermediateVars\n\n";
  }
#endif

  Vpos1 = Vpos2 = Vneg1 = Vneg2 = Vint1 = Vint2 = 0.0;

  Vpos1 = solVec[li_Pos1];
  Vneg1 = solVec[li_Neg1];
  Vint1 = solVec[li_Int1];
  Ibr1  = solVec[li_Ibr1];
  Vpos2 = solVec[li_Pos2];
  Vneg2 = solVec[li_Neg2];
  Vint2 = solVec[li_Int2];
  Ibr2  = solVec[li_Ibr2];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << " Vpos1 = " << Vpos1 << std::endl;
    std::cout << " Vneg1 = " << Vneg1 << std::endl;
    std::cout << " Vint1 = " << Vint1 << std::endl;
    std::cout << " Ibr1 = " << Ibr1 << std::endl;
    std::cout << " Vpos2 = " << Vpos2 << std::endl;
    std::cout << " Vneg2 = " << Vneg2 << std::endl;
    std::cout << " Vint2 = " << Vint2 << std::endl;
    std::cout << " Ibr2 = " << Ibr2 << std::endl;
  }
#endif

  // Test if we're doing DC or Transient
  if ((getSolverState().dcopFlag))
  {
    // DC operation
    DCMODE=true;
    v1 = (Vpos2-Vneg2)+Z0*Ibr2;
    v2 = (Vpos1-Vneg1)+Z0*Ibr1;
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      std::cout << "DC Mode, V1 = " << v1 <<  ", V2 = " << v2  << std::endl;
#endif
  }
  else
  {
    double currentTime = getSolverState().currTime;
    DCMODE=false;
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      std::cout << "Not DC, newtonIter = " << getSolverState().newtonIter;
      std::cout << " Time is " << currentTime << std::endl;
      std::cout << "        newtonIterOld = " << newtonIterOld << " timeOld is " << timeOld << std::endl;
    }
#endif
    // Transient operation
    // Now determine if we're on the first newton step of an iteration
    if (getSolverState().newtonIter == 0 && (currentTime != timeOld))
    {
      timeOld = currentTime;
      // we are, so need to manipulate history and calculate v1,v2.
      // If we're the first time step, we need to initialize it
      if (getSolverState().initTranFlag)
      {
        last_t = currentTime;
        v1 = (Vpos2-Vneg2)+Z0*Ibr2;
        v2 = (Vpos1-Vneg1)+Z0*Ibr1;

        history.clear();
        history.push_back(History(-2*td,v1,v2));
        history.push_back(History(-td,v1,v2));
        history.push_back(History(0,v1,v2));
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          std::cout << "Transient, first time, T = ";
          std::cout << currentTime;
          std::cout << ", V1 = " << v1;
          std::cout << ", V2 = " << v2  << std::endl;
        }
#endif
      }
      else
      {
        double delayedTime = currentTime-td;

        // now get the values of v1 and v2 from the delayed-time
        // information
        InterpV1V2FromHistory(delayedTime, &v1, &v2);
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          std::cout << "  Done with interpolation to delayedTime=";
          std::cout << delayedTime;
          std::cout << ", have v1="<<v1 << " and v2=" << v2 << std::endl;
          std::cout << " INTERP " << delayedTime << " " << v1 << " " << v2 << std::endl;
          std::cout << " Set last_t to " << currentTime << std::endl;
        }
#endif
        // now save the current time so we can have it next time
        // we get to this block (i.e. on the next time step)
        last_t = currentTime;
      }
    }
    else
    {
      // we're on the second iteration or later of the second time
      // step or later.  Re-use the values of v1 and v2 from the
      // first iteration  of this time step.  We don't care  what time
      // it is
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        std::cout << "second or later iteration, t is " << currentTime;
        std::cout << " have last_t = " << last_t << " v1="<<v1 << " and v2=" << v2 << std::endl;
      }
#endif
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pruneHistory
// Purpose       : sift through the transmission line state history and
//                 delete records that are so old they'll never be used again
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 6/15/2001
//-----------------------------------------------------------------------------

void Instance::pruneHistory(double t1)
{

  // The input t is the oldest time for which we'll ever interpolate again.
  // That means we only need two times in the history that are older than t1,
  // so this routine drops everything off the head but the most recent 2 that
  // are older than t.

  vector<History>::iterator first = history.begin();
  vector<History>::iterator it1;
  vector<History>::iterator last = history.end();
  int i;

  last--; // point to last stored item, not end of the list!
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << "----------------------------------" << std::endl;
    std::cout << "Pruning for time t1="<<t1 << std::endl;
    std::cout << " Oldest in list is t="<<first->t<<" v1 = "<<first->v1 <<
      " v2="<<first->v2 << std::endl;
    std::cout << " latest in list is t="<<last->t<<" v1 = "<<last->v1
         << " v2="<<last->v2 << std::endl;
  }
#endif
  // First find the first element for which the stored time is greater than t
  for (it1 = first, i = 0; it1->t < t1 && it1 != last; ++it1, ++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      std::cout << "i = " << i << " t = " << it1->t;
      std::cout << " v1 = " << it1->v1;
      std::cout << " v2 = " << it1->v2 << std::endl;
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << "   i ="<<i << std::endl;
  }
#endif

  //   Now it1 points to the first element with t>t1
  if (i > 2)
    {
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        std::cout << "Need to prune.  Keeping " << it1->t << std::endl;
      }
#endif
      // if i>2 we have  too many old ones.
      // back up 2
      it1--;
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        std::cout << "                Keeping " << it1->t << std::endl;
      }
#endif

      it1--;
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        std::cout << "                Keeping " << it1->t << std::endl;
      }
#endif
      // delete everything from the first to it1, not counting it1
      history.erase(first,it1);
    }
  // otherwise we don't need to do anything.
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << "----------------------------------" << std::endl;
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : Instance::InterpV1V2FromHistory
// Purpose       : Use 3-point lagrange interpolation to determine
//                 v1(t) and v2(t) at a specified time in the past
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 6/15/2001
//-----------------------------------------------------------------------------
void Instance::InterpV1V2FromHistory(double t, double * v1p,
                                              double *v2p)
{
  vector<History>::iterator first = history.begin();
  vector<History>::iterator it1;
  vector<History>::iterator last = history.end();
  double t1,t2,t3;
  double dt1,dt2,dt3;
  double v11,v21,v12,v22,v13,v23;
  double dt12,dt13,dt23;
  double f1,f2,f3;    // interpolating functions

  if (history.size() <= 0)
  {
    string msg;
    msg="Instance::InterpV1V2FromHistory called but history list is"
      " empty.  Might be due to trying to restart this netlist.\n"
      "Restarts of netlists with transmission lines does not work yet.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  last--; // point to the last stored item, not the tail of the list!
  // sanity clause (you canna foola me, I know they're ain'ta no sanity
  // clause!)
  //  if (t < first->t || t > last->t)
  if (t - first->t < -N_UTL_MachineDependentParams::MachinePrecision()
      || t - last->t > N_UTL_MachineDependentParams::MachinePrecision() )
  {
    string msg;
    char msg2[256];
    sprintf(msg2, "Instance::InterpV1V2FromHistory: "
            "Asked to interpolate to a time (%20.17lg) prior to oldest(%20.17g) or "
            "after newest(%20.17lg) in history\n",t,first->t,last->t);
    msg = msg2;
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << " interpolating for t = " << t  << std::endl;
  }
#endif

  // If we are within roundoff of the endpoints of the history, just use
  // the endpoints, otherwise interpolate to get it.
  if ( fabs(t-first->t)<N_UTL_MachineDependentParams::MachinePrecision())
  {
    *v1p = first->v1;
    *v2p = first->v2;
  }
  else if ( fabs(t-last->t)<N_UTL_MachineDependentParams::MachinePrecision())
  {
    *v1p = last->v1;
    *v2p = last->v2;
  }
  else
  {

    LessThan<History,double> lessFunct;
    it1 = lower_bound(history.begin(),history.end(),t,lessFunct);

    // Now it1 points to the first element with time > t
    t3 = it1->t;
    v13 = it1->v1;
    v23 = it1->v2;
    it1--;
    t2 = it1->t;
    v12 = it1->v1;
    v22 = it1->v2;
    it1--;
    t1 = it1->t;
    v11 = it1->v1;
    v21 = it1->v2;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      std::cout << "Using time t3="<<t3<<" v1(t3)="<<v13<<" v2(t3)="<<v23  << std::endl;
      std::cout << "Using time t2="<<t2<<" v1(t2)="<<v12<<" v2(t2)="<<v22  << std::endl;
      std::cout << "Using time t1="<<t1<<" v1(t1)="<<v11<<" v2(t1)="<<v21  << std::endl;
    }
#endif

    //  now we have three values of each function to be interpolated, and three
    // times.  t3 is after the desired time, t1 and t2 are before (t2 might be
    // equal to the desired time)
    // Set up the differences for lagrange interpolation:
    dt12 = t1-t2;
    dt13 = t1-t3;
    dt23 = t2-t3;
    dt1 = t-t1;
    dt2 = t-t2;
    dt3 = t-t3;
    // now we set up the lagrange interpolating functions
    // e.g. f1 = (t-t2)*(t-t3)/((t1-t2)*(t1-t3))
    // so that fi is 1 at ti and 0 at the other times.
    f1 = dt2*dt3;
    f2 = dt1*dt3;
    f3 = dt1*dt2;
    if (dt12 != 0)
    {
      f1 /= dt12;
      f2 /= -dt12;
    }
    else
    {
      f1 = f2 = 0.0;
    }
    if (dt13 != 0)
    {
      f1 /= dt13;
      f3 /= -dt13;
    }
    else
    {
      f1 = f2 = 0.0;
    }
    if (dt23 != 0)
    {
      f2 /= dt23;
      f3 /= -dt23;
    }
    else
    {
      f2 = f3 = 0.0;
    }
    // that's it, we have the interpolation functions evaluated at the time t,
    // and the values of v1 and v2 at the points, perform  the interpolation
    *v1p = f1*v11+f2*v12+f3*v13;
    *v2p = f1*v21+f2*v22+f3*v23;
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
//                 It does not bother to check them in any way, or put them
//                 in order.  It only adds them in.
//
// Special Notes : The guts of this has been moved to acceptStep, which
//                 actually computes the breakpoints if needed.  We only add
//                 them to the list here if necessary.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/08/01
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints ( vector<N_UTL_BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  double currentTime = getSolverState().currTime;
  int timeStep = getSolverState().timeStepNumber;

  //  We're called once prior to any newton iterations, not even the
  // DC Op point.  Never do anything if first_BP_call_done is false.

  if (timeStep != 0 && first_BP_call_done)
  {
    if (newBreakPoint)
    {
      breakPointTimes.push_back(newBreakPointTime);
      newBreakPoint = false;
    }
  }
  else
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        std::cout << " In Instance::getBreakPoints "<<std::endl;
        std::cout << " First time step, I don't get to set breakpoints.  Time is ";
        std::cout << currentTime << std::endl;
      }
#endif
  }

  first_BP_call_done=true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       : This function saves the values of v1 and v2 along with
//                 the current time.  It is to be called ONLY at the point
//                 when the time integrator has determined we've got a
//                 converged, acceptable solution and is accepting it,
//                 but before it's updated its times and rotated vectors.
//
// Special Notes : In SPICE this same stuff was done in the "TRAaccept" function.
//
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 01/23/07
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  if (!getSolverState().dcopFlag)
  {
    double currentTime = getSolverState().currTime;

    double d11, d21, d12, d22;
    N_LAS_Vector *theSolVectorPtr = extData.nextSolVectorPtr;// the accepted
    // values from this
    // step

    vector<History>::iterator last = history.end();

    last--;  // point to last item, not past last item.

    //  We're called once prior to any newton iterations, not even the
    // DC Op point.  Never do anything if first_BP_call_done is false.
    double oVp1,oVp2,oVn1,oVn2,oI1,oI2;
    double ov1,ov2;
    double tmp_v1,tmp_v2, tmp_t;

#ifdef Xyce_DEBUG_DEVICE
    double oVi1,oVi2;
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      std::cout << " In Instance::acceptStep "<<std::endl;
      std::cout << "I want breakpoints.  Time is " << currentTime << std::endl;
      std::cout << "   timeOld is " << timeOld << std::endl;
    }
#endif

    // we're the end of a time step, the solution has been accepted.
    // clean up the history by deleting records of times so far
    // back that they'll never be used for interpolation again
    // never try to prune history for anything but times that have
    // been accepted already.
    // TVR: The goal of this was to prune the early history so we don't
    // get unbounded growth of the history vector, with the intent of making
    // the interpolation method faster.  Turns out that deleting these vector
    // elements is very expensive, much more expensive than using "lower_bound"
    // to find a value in the long list.  So I'm commenting this out.
    // double delayedTime;
    //  if (timeOld != -1)
    //  {
    //    delayedTime = timeOld-td;
    //    pruneHistory(delayedTime);
    //  }

    oVp1 = (*theSolVectorPtr)[li_Pos1];
    oVn1 = (*theSolVectorPtr)[li_Neg1];
    oI1  = (*theSolVectorPtr)[li_Ibr1];
    oVp2 = (*theSolVectorPtr)[li_Pos2];
    oVn2 = (*theSolVectorPtr)[li_Neg2];
    oI2  = (*theSolVectorPtr)[li_Ibr2];

    // Having the old values means we can calculate what v1 and v2
    // were for that time.
    ov1=(oVp2-oVn2)+Z0*oI2;
    ov2=(oVp1-oVn1)+Z0*oI1;

#ifdef Xyce_DEBUG_DEVICE
    oVi1 = (*theSolVectorPtr)[li_Int1];
    oVi2 = (*theSolVectorPtr)[li_Int2];

    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      std::cout << " ----- New time step -----" << std::endl;
      std::cout << " Last solution : " << std::endl;
      std::cout << " vpos1 = " << oVp1 << std::endl;
      std::cout << " vneg1 = " << oVn1 << std::endl;
      std::cout << " vint1 = " << oVi1 << std::endl;
      std::cout << " ibr1 = " << oI1 << std::endl;
      std::cout << " vpos2 = " << oVp2 << std::endl;
      std::cout << " vneg2 = " << oVn2 << std::endl;
      std::cout << " vint2 = " << oVi2 << std::endl;
      std::cout << " ibr2 = " << oI2 << std::endl;
      std::cout << "in set breakpoints, saving for time=" << currentTime << ", V1 = " << ov1 <<  ", V2 = " << ov2  << std::endl;
      std::cout << " V1V2DBG " << currentTime << " " << ov1 << " " << ov2 << std::endl;
    }
#endif

    history.push_back(History(currentTime,ov1,ov2));

    last = history.end();
    last--; // point to last item, not past last item.
    // Now calculate derivatives based on history
    tmp_v1 = last->v1; tmp_v2 = last->v2; tmp_t = last->t;
    last--;
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      std::cout << "tmp_t="  << tmp_t  << " last->t =" << last->t << std::endl;
      std::cout << "tmp_v1=" << tmp_v1 << " last->v1=" << last->v1 << std::endl;
      std::cout << "tmp_v2=" << tmp_v2 << " last->v2=" << last->v2 << std::endl;
    }
#endif
    d11 = (tmp_v1-last->v1)/(tmp_t-last->t);
    d12 = (tmp_v2-last->v2)/(tmp_t-last->t);
    tmp_v1 = last->v1; tmp_v2 = last->v2; tmp_t = last->t;
    last--;
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      std::cout << "tmp_t="  << tmp_t  << " last->t =" << last->t << std::endl;
      std::cout << "tmp_v1=" << tmp_v1 << " last->v1=" << last->v1 << std::endl;
      std::cout << "tmp_v2=" << tmp_v2 << " last->v2=" << last->v2 << std::endl;
    }
#endif
    d21 = (tmp_v1-last->v1)/(tmp_t-last->t);
    d22 = (tmp_v2-last->v2)/(tmp_t-last->t);
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      std::cout << "Derivs are " << d11 << " " << d21 << std::endl;
      std::cout << "   and " << d12 << " " <<d22 << std::endl;
      std::cout << " fabs(d11-d21) = " << fabs(d11-d21) << std::endl;
      std::cout << " fabs(d12-d22) = " << fabs(d12-d22) << std::endl;
      std::cout << "D1D2DBG " << currentTime << " " << d11 << " " << d12 << std::endl;
    }
#endif

    if ((fabs(d11-d21) >= .99*Xycemax(fabs(d11),fabs(d21))+1) ||
        (fabs(d12-d22) >= .99*Xycemax(fabs(d12),fabs(d22))+1))
    {
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        std::cout << "Derivative is changing enough, I want to set a break point ";
        std::cout << td << " ahead of discontinuity, which is ";
        std::cout << tmp_t+td<<std::endl;
      }
#endif
      newBreakPointTime = (tmp_t+td);
      newBreakPoint = true;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInternalState
// Purpose       : Generates an DeviceState object and populates
//                 it with the contents of the history vector for use by
//                 restarts
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 09/03/04
//-----------------------------------------------------------------------------

DeviceState * Instance::getInternalState()
{
  int hsize,i,j;
  // allocate object to return
  DeviceState * myState = new DeviceState;


  myState->ID=getName();
  // We'll pack our history data into the single vector of doubles
  myState->data.resize(history.size()*3);
  hsize=history.size();
  for (i=0;i<hsize;++i)
  {
    j=i*3;
    myState->data[j]=history[i].t;
    myState->data[j+1]=history[i].v1;
    myState->data[j+2]=history[i].v2;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << "---------------------------------------------------------------- "
         << std::endl;
    std::cout << " In Instance::getInternalState " << std::endl;
    std::cout << "   name=" << getName() << std::endl;
    std::cout << "   history size = " << hsize << std::endl;
    std::cout << "   history  data: " << std::endl;
    for (i = 0 ; i < hsize ; ++i)
    {
      std::cout << "   (" << history[i].t << ", " << history[i].v1 << ", "
           << history[i].v2 << ")"<< std::endl;
    }

    std::cout << "   DeviceState ID = " << myState->ID << std::endl;
    std::cout << "   DeviceState data size " << myState->data.size() << std::endl;
    std::cout << "   Device State data: " << std::endl;
    for (i = 0 ; i < myState->data.size() ; ++i)
    {
      std::cout << "    " << myState->data[i] << std::endl;
    }
    std::cout << "---------------------------------------------------------------- "
         << std::endl;
  }
#endif

  return myState;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setInternalState
// Purpose       : Reload history data from restart
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 09/03/04
//-----------------------------------------------------------------------------
bool Instance::setInternalState(const DeviceState &state)
{
  int dsize=state.data.size();
  int hsize,i,j;
  if ( state.ID != getName())
  {
    string msg;
    msg = "Instance::setInternalState:  ID ("+state.ID+")";
    msg += "from restart does not match my name ("+getName()+")!\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  if (dsize%3 != 0)
  {
    string msg;
    char msg2[256];
    sprintf(msg2, "Instance::setInternalState: "
            "Data size from restart (%d) not a multiple of 3!",
            dsize);
    msg=msg2;
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  hsize=dsize/3;
  history.clear();
  history.resize(hsize);
  for ( i=0; i<hsize; ++i)
  {
    j=i*3;
    history[i].t=state.data[j];
    history[i].v1=state.data[j+1];
    history[i].v2=state.data[j+2];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << "---------------------------------------------------------------- "
         << std::endl;
    std::cout << " In Instance::setInternalState " << std::endl;
    std::cout << "   name=" << getName() << std::endl;
    std::cout << "   history size = " << hsize << std::endl;
    std::cout << "   history  data: " << std::endl;
    for (i = 0 ; i < hsize ; ++i)
    {
      std::cout << "   (" << history[i].t << ", " << history[i].v1 << ", "
           << history[i].v2 << ")"<< std::endl;
    }

    std::cout << "   DeviceState ID = " << state.ID << std::endl;
    std::cout << "   DeviceState data size " << state.data.size() << std::endl;
    std::cout << "   Device State data: " << std::endl;
    for (i = 0 ; i < state.data.size() ; ++i)
    {
      std::cout << "    " << state.data[i] << std::endl;
    }
    std::cout << "---------------------------------------------------------------- "
         << std::endl;
  }
#endif
  return true;
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 9/25/02
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{
  // there are no model parameters to process.
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams(string param)
{

  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                      SolverState & ss1,
                                      DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1)
{
}


//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::~Model ()
{
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     getModelName()  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << (*iter)->getModelName();

    os << std::endl;
    os << "Z0 = " << (*iter)->Z0 << std::endl;
    os << "G0 = " << (*iter)->G0 << std::endl;
    os << "TD = " << (*iter)->td << std::endl;
    os << "FREQ = " << (*iter)->freq << std::endl;
    os << "NL = " << (*iter)->NL << std::endl;

    os << std::endl;
  }
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 8/01/01
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize ()
{
  return td;
}

// Additional Declarations

// History member (trivial) functions

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : default constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/14/01
//-----------------------------------------------------------------------------
History::History()
  : t(0),v1(0),v2(0)
{
}

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/14/01
//-----------------------------------------------------------------------------
History::~History()
{
}
//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/14/01
//-----------------------------------------------------------------------------
History::History(const History &right)
  : t(right.t),v1(right.v1),v2(right.v2)
{
}

//-----------------------------------------------------------------------------
// Function      : History::History
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/14/01
//-----------------------------------------------------------------------------
History::History(double a, double b, double c)
  : t(a),v1(b),v2(c)
{
}

} // namespace TRA
} // namespace Device
} // namespace Xyce
