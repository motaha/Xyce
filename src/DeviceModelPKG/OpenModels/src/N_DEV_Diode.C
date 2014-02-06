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
// Filename       : $RCSfile: N_DEV_Diode.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.268.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_Diode.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Diode::Instance>::ParametricData()
{
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(1);
  addModelType("D");

  // Set up double precision variables:
  addPar ("AREA",  1.0, false,   ParameterType::NO_DEP,
          &Diode::Instance::Area,
          NULL,
          U_NONE, CAT_GEOMETRY, "Area scaling value (scales IS, ISR, IKF, RS, CJ0, and IBV)");

  addPar ("IC",    0.0, false,   ParameterType::NO_DEP,
          &Diode::Instance::InitCond,
          &Diode::Instance::InitCondGiven,
          STANDARD, CAT_NONE, "");

  addPar ("TEMP",  0.0, false, ParameterType::TIME_DEP,
          &Diode::Instance::Temp,
          NULL,
          STANDARD, CAT_NONE, "");

  // Set up non-double precision variables:
  addPar ("OFF", 0, false, ParameterType::NO_DEP,
          &Diode::Instance::off, NULL,
          U_LOGIC, CAT_CONTROL, "Initial voltage drop across device set to zero");
  addPar ("LAMBERTW", 0, false, ParameterType::NO_DEP,
          &Diode::Instance::lambertWFlag, NULL,
          U_LOGIC, CAT_CONTROL, "Option to solve diode equations with the Lambert-W function");
}

template<>
ParametricData<Diode::Model>::ParametricData()
{
  addPar ("IS", 1.0e-14, false, ParameterType::NO_DEP,
          &Diode::Model::IS,
          NULL,
          U_AMP, CAT_CURRENT, "Saturation current");

  // synonym for IS
  addPar ("JS", 1.0e-14, false, ParameterType::NO_DEP,
          &Diode::Model::IS,
          NULL,
          U_AMP, CAT_CURRENT, "Saturation current");

  addPar ("RS",  0.0, false, ParameterType::MIN_RES,
          &Diode::Model::RS,
          NULL,
          U_OHM, CAT_RES, "Parasitic resistance");

  addPar ("N",  1.0, true, ParameterType::NO_DEP,
          &Diode::Model::N,
          NULL,
          U_NONE, CAT_PROCESS, "Emission coefficient");

  addPar ("ISR", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::ISR,
          NULL,
          U_AMP, CAT_CURRENT, "Recombination current parameter (level 2)");

  addPar ("NR", 2.0, false, ParameterType::NO_DEP,
          &Diode::Model::NR,
          NULL,
          U_NONE, CAT_NONE, "Emission coefficient for ISR (level 2)");

  addPar ("IKF", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::IKF,
          NULL,
          U_AMP, CAT_CURRENT, "High-injection \"knee\" current (level 2)");

  addPar ("TT", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::TT,
          NULL,
          U_SECOND, CAT_PROCESS, "Transit time");

  addPar ("CJO", 0.0, false, ParameterType::MIN_CAP,
          &Diode::Model::CJO,
          NULL,
          U_FARAD, CAT_CAP, "Zero-bias p-n depletion capacitance");

  // synonyms for CJO
  addPar ("CJ", 0.0, false, ParameterType::MIN_CAP,
          &Diode::Model::CJO,
          NULL,
          U_FARAD, CAT_CAP, "Zero-bias p-n depletion capacitance");
  addPar ("CJ0", 0.0, false, ParameterType::MIN_CAP,
          &Diode::Model::CJO,
          NULL,
          U_FARAD, CAT_CAP, "Zero-bias p-n depletion capacitance");

  addPar ("VJ", 1.0, false, ParameterType::NO_DEP,
          &Diode::Model::VJ,
          NULL,
          U_VOLT, CAT_VOLT, "Potential for p-n junction");

  addPar ("M", 0.5, false, ParameterType::NO_DEP,
          &Diode::Model::M,
          NULL,
          U_NONE, CAT_PROCESS, "Grading parameter for p-n junction");

  addPar ("EG", 1.11, false, ParameterType::NO_DEP,
          &Diode::Model::EG,
          NULL,
          U_EV, CAT_PROCESS, "Bandgap voltage (barrier height)");

  addPar ("XTI", 3.0, false, ParameterType::NO_DEP,
          &Diode::Model::XTI,
          NULL,
          U_NONE, CAT_TEMP, "IS temperature exponent");

  addPar ("TIKF", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::TIKF,
          NULL,
          U_DEGCM1, CAT_TEMP, "IKF temperature coefficient (linear) (level 2)");

  addPar ("TBV1", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::TBV1,
          NULL,
          U_DEGCM1, CAT_TEMP, "BV temperature coefficient (linear) (level 2)");

  addPar ("TBV2", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::TBV2,
          NULL,
          U_DEGCM2, CAT_TEMP, "BV temperature coefficient (quadratic) (level 2)");

  addPar ("TRS1", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::TRS1,
          NULL,
          U_DEGCM1, CAT_TEMP, "RS temperature coefficient (linear) (level 2)");

  addPar ("TRS2", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::TRS2,
          NULL,
          U_DEGCM2, CAT_TEMP, "RS temperature coefficient (quadratic) (level 2)");
  addPar ("FC", 0.5, false, ParameterType::NO_DEP,
          &Diode::Model::FC,
          NULL,
          U_NONE, CAT_CAP, "Forward-bias depletion capacitance coefficient");

  addPar ("BV", 1E99, false, ParameterType::NO_DEP,
          &Diode::Model::BV,
          &Diode::Model::BVGiven,
          U_VOLT, CAT_VOLT, "Reverse breakdown \"knee\" voltage");

  // synonym for BV
  addPar ("VB", 1E99, false, ParameterType::NO_DEP,
          &Diode::Model::BV,
          &Diode::Model::BVGiven,
          U_VOLT, CAT_VOLT, "Reverse breakdown \"knee\" voltage");

  addPar ("IBV", 1.0e-3, false, ParameterType::NO_DEP,
          &Diode::Model::IBV,
          NULL,
          U_AMP, CAT_CURRENT, "Reverse breakdown \"knee\" current");

  addPar ("IRF", 1.0, false, ParameterType::NO_DEP,
          &Diode::Model::IRF,
          NULL,
          U_NONE, CAT_CURRENT, "Reverse current fitting factor");

  addPar ("NBV", 1.0, false, ParameterType::NO_DEP,
          &Diode::Model::NBV,
          NULL,
          U_NONE, CAT_PROCESS, "Reverse breakdown ideality factor (level 2)");

  addPar ("IBVL", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::IBVL,
          NULL,
          U_AMP, CAT_CURRENT, "Low-level reverse breakdown \"knee\" current (level 2)");

  addPar ("NBVL", 1.0, false, ParameterType::NO_DEP,
          &Diode::Model::NBVL,
          NULL,
          U_NONE, CAT_PROCESS, "Low-level reverse breakdown ideality factor (level 2)");

  addPar ("TNOM", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::TNOM,
          NULL,
          STANDARD, CAT_NONE, "");

  addPar ("KF", 0.0, false, ParameterType::NO_DEP,
          &Diode::Model::KF,
          NULL,
          U_NONE, CAT_FLICKER, "Flicker noise coefficient");

  addPar ("AF", 1.0, false, ParameterType::NO_DEP,
          &Diode::Model::AF,
          NULL,
          U_NONE, CAT_FLICKER, "Flicker noise exponent");
}

  namespace Diode {

vector< vector<int> > Instance::jacStamp_RS;
vector< vector<int> > Instance::jacStamp;

vector<int> Instance::jacMap_RS;
vector<int> Instance::jacMap;

vector< vector<int> > Instance::jacMap2_RS;
vector< vector<int> > Instance::jacMap2;



ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams(string param)
{
  updateTemperature( Temp );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
                   Model & Miter,
                   MatrixLoadData & mlData1,
                   SolverState &ss1,
                   ExternData  &ed1,
                   DeviceOptions & do1)

  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(Miter),
    off(0),
    Area(1.0),
    InitCond(0.0),
    Temp(getDeviceOptions().temp.dVal()),
    lambertWFlag(0),
    InitCondGiven(false),
    tJctPot(0.0),
    tJctCap(0.0),
    tDepCap(0.0),
    tSatCur(0.0),
    tVcrit(0.0),
    tF1(0.0),
    tBrkdwnV(0.0),
    tSatCurR(0.0),
    tIKF(0.0),
    tRS(0.0),
    tIRF(1.0),
    Id(0.0),
    Gd(0.0),
    Cd(0.0),
    Gcd(0.0),
    Icd(0.0),
    Qd(0.0),
    Gspr(0.0),
    Vpp(0.0),
    Vp(0.0),
    Vn(0.0),
    Vc(0.0),
    Vd(0.0),
    Vd_old(0.0),
    Vd_orig(0.0),
    newtonIterOld(0),
    li_Pos(-1),
    li_Neg(-1),
    li_Pri(-1),
    APosEquPosNodeOffset(-1),
    APosEquPriNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    ANegEquPriNodeOffset(-1),
    APriEquPosNodeOffset(-1),
    APriEquNegNodeOffset(-1),
    APriEquPriNodeOffset(-1),
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    fPosEquPosNodePtr(0),
    fPosEquPriNodePtr(0),
    fNegEquNegNodePtr(0),
    fNegEquPriNodePtr(0),
    fPriEquPosNodePtr(0),
    fPriEquNegNodePtr(0),
    fPriEquPriNodePtr(0),
    qPosEquPosNodePtr(0),
    qPosEquPriNodePtr(0),
    qNegEquNegNodePtr(0),
    qNegEquPriNodePtr(0),
    qPriEquPosNodePtr(0),
    qPriEquNegNodePtr(0),
    qPriEquPriNodePtr(0),
#endif
    li_storevd(-1),
    li_store_dev_i(-1)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;
  setNumStoreVars(1);
  numLeadCurrentStoreVars = 1; // lead current DEV_I

  if( jacStamp.empty() )
  {
    jacStamp_RS.resize(3);
    jacStamp_RS[0].resize(2);
    jacStamp_RS[0][0]=0;
    jacStamp_RS[0][1]=2;
    jacStamp_RS[1].resize(2);
    jacStamp_RS[1][0]=1;
    jacStamp_RS[1][1]=2;
    jacStamp_RS[2].resize(3);
    jacStamp_RS[2][0]=0;
    jacStamp_RS[2][1]=1;
    jacStamp_RS[2][2]=2;

    jacMap_RS.clear();
    jacStampMap(jacStamp_RS, jacMap_RS, jacMap2_RS,
                jacStamp,    jacMap, jacMap2, 2, 0, 3);

  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("LAMBERTW"))
    lambertWFlag = getDeviceOptions().lambertWFlag;
  if (!given("TEMP"))
    Temp = getDeviceOptions().temp.dVal();
  if ( (model_.RS == 0.0) || (lambertWFlag == 1) )
    numIntVars = 0;
  if ( model_.CJO == 0.0 )
    numStateVars = 1;
  if ( (model_.RS == 0.0) && (lambertWFlag == 1) )
    model_.RS =1.0e-12;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

// Additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       : function for registering, and setting up, local ID's.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                             const vector<int> & extLIDVecRef)
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << dashedline << endl;
    cout << "In Instance::register LIDs\n\n";
    cout << "name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "number of internal variables: " << numInt << endl;
    cout << "number of external variables: " << numExt << endl;
  }
#endif

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  //Test for Ohmic Resistance.  If RS=0, Pri=Pos Node
  //-------------------------------------------------
  if ( numInt != numIntVars )
  {
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  //Setup of Pri node indices
  //If RS=0, Pri=Pos Node
  //--------------------------
  if( model_.RS && (lambertWFlag != 1) )
    li_Pri = intLIDVec[0];
  else
    li_Pri = li_Pos;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "\nSolution and RHS variables:\n";
    cout << "\nli_Pos = ";
    cout.width(4);
    cout << li_Pos << endl;

    cout << "\nli_Neg = ";
    cout.width(4);
    cout << li_Neg << endl;

    cout << "\nli_Pri = ";
    cout.width(4);
    cout << li_Pri << endl;

  }
#endif



#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "\nEnd of Instance::register LIDs\n";
    cout << dashedline << endl;
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
    string tmpstr;
    if ( li_Pos != li_Pri )
    {
      tmpstr = getName()+"_internal";
      spiceInternalName (tmpstr);
      intNameMap[ li_Pri ] = tmpstr;
    }
  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
  string msg;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl;
    cout << dashedline << endl;
    cout << "  In Instance::registerStateLIDs\n\n";
    cout << "  name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists equals the proper # of variables.
  int numSta = staLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  Number of State LIDs: " << numSta << endl;
  }
#endif
  //-------------------------------------------------------------
  if( numSta != numStateVars )
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  i=0;
  // copy over the global ID lists:
  staLIDVec = staLIDVecRef;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(
  const vector<int> & stoLIDVecRef )
{

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl;
    cout << dashedline << endl;
    cout << "  In Instance::registerStoreLIDs\n\n";
    cout << "  name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists equals the proper # of variables.
  int numSto = stoLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  Number of Store LIDs: " << numSto << endl;
  }
#endif

  //-------------------------------------------------------------
  if( numSto != getNumStoreVars() )
  {
    string msg = "Instance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;
  li_storevd = stoLIDVec[0];
  if( loadLeadCurrent )
  {
    li_store_dev_i = stoLIDVec[1];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "li_storevd = " << li_storevd;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_DEV_DiodeInstance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/27/2013
//-----------------------------------------------------------------------------
std::map<int, std::string> & N_DEV_DiodeInstance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if( storeNameMap.empty () )
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    string modName(getName());
    spiceInternalName(modName);
    string tmpstr;
    tmpstr = modName+":vd";
    storeNameMap[ li_storevd ] = tmpstr;
    if( loadLeadCurrent )
    {
      tmpstr = modName+":DEV_I";
      storeNameMap[ li_store_dev_i ] = tmpstr;
    }
  }
  return storeNameMap;
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_DiodeInstance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  if( model_.RS && (lambertWFlag != 1) )
    return jacStamp_RS;
  else
    return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  vector<int> map;
  vector< vector<int> > map2;

  if( model_.RS && (lambertWFlag != 1) )
  {
    map = jacMap_RS;
    map2 = jacMap2_RS;
  }
  else
  {
    map = jacMap;
    map2 = jacMap2;
  }

  APosEquPosNodeOffset = jacLIDVec[map[0]][map2[0][0]];
  APosEquPriNodeOffset = jacLIDVec[map[0]][map2[0][1]];

  ANegEquNegNodeOffset = jacLIDVec[map[1]][map2[1][0]];
  ANegEquPriNodeOffset = jacLIDVec[map[1]][map2[1][1]];

  APriEquPosNodeOffset = jacLIDVec[map[2]][map2[2][0]];
  APriEquNegNodeOffset = jacLIDVec[map[2]][map2[2][1]];
  APriEquPriNodeOffset = jacLIDVec[map[2]][map2[2][2]];

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  fPosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  fPosEquPriNodePtr = &(dFdx[li_Pos][APosEquPriNodeOffset]);
  fNegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
  fNegEquPriNodePtr = &(dFdx[li_Neg][ANegEquPriNodeOffset]);
  fPriEquPosNodePtr = &(dFdx[li_Pri][APriEquPosNodeOffset]);
  fPriEquNegNodePtr = &(dFdx[li_Pri][APriEquNegNodeOffset]);
  fPriEquPriNodePtr = &(dFdx[li_Pri][APriEquPriNodeOffset]);

  qPosEquPosNodePtr = &(dQdx[li_Pos][APosEquPosNodeOffset]);
  qPosEquPriNodePtr = &(dQdx[li_Pos][APosEquPriNodeOffset]);
  qNegEquNegNodePtr = &(dQdx[li_Neg][ANegEquNegNodeOffset]);
  qNegEquPriNodePtr = &(dQdx[li_Neg][ANegEquPriNodeOffset]);
  qPriEquPosNodePtr = &(dQdx[li_Pri][APriEquPosNodeOffset]);
  qPriEquNegNodePtr = &(dQdx[li_Pri][APriEquNegNodeOffset]);
  qPriEquPriNodePtr = &(dQdx[li_Pri][APriEquPriNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 diode instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  // note: the only capacitor goes from the negative to the
  // positive-prime node, so there is not contribution in this
  // function to the positive node.

  // load in the KCL for the negative node:
  double coef = Qd;
  (extData.daeQVectorRawPtr)[li_Neg] -= coef;

  // load in the KCL for the positive prime node:
  coef *= -1.0; // -Qd
  (extData.daeQVectorRawPtr)[li_Pri] -= coef;

  // load the voltage limiter vector.
  if( getDeviceOptions().voltageLimiterFlag )
  {
    double Vd_diff = Vd - Vd_orig;
    double Cd_Jdxp = 0.0;
    Cd_Jdxp = -( Cd ) * Vd_diff;

    // Load the dQdxdVp vector
    (extData.dQdxdVpVectorRawPtr)[li_Neg] += Cd_Jdxp;
    (extData.dQdxdVpVectorRawPtr)[li_Pri] -= Cd_Jdxp;
  }

  if( loadLeadCurrent && (model_.CJO != 0.0))
  {
    (extData.storeLeadCurrQCompRawPtr)[li_store_dev_i] = Qd;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 diode instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  // 3f5 compatible currents
  // Including derivation of Vd_diff and Limiting Correction
  //---------------------------------------------------------
  double Ir = Gspr * (Vp - Vpp);

  // load in the KCL for the positive node:
  double coef = -Ir;
  (extData.daeFVectorRawPtr)[li_Pos] -= coef;

  // load in the KCL for the negative node:
  coef = Id;
  (extData.daeFVectorRawPtr)[li_Neg] -= coef;

  // load in the KCL for the positive prime node:
  coef *= -1;
  coef += Ir;
  (extData.daeFVectorRawPtr)[li_Pri] -= coef;

  // load the voltage limiter vector.
  if( getDeviceOptions().voltageLimiterFlag )
  {
    double Vd_diff = Vd - Vd_orig;
    double Gd_Jdxp = 0.0;
    Gd_Jdxp = -( Gd ) * Vd_diff;

    // Load the dFdxdVp vector
    (extData.dFdxdVpVectorRawPtr)[li_Neg] += Gd_Jdxp;
    (extData.dFdxdVpVectorRawPtr)[li_Pri] -= Gd_Jdxp;
  }

  if( loadLeadCurrent )
  {
    (extData.nextStoVectorRawPtr)[li_store_dev_i] = Id;
  }


  if( loadLeadCurrent )
  {
    (extData.nextStoVectorRawPtr)[li_store_dev_i] = Id;
  }


  if( loadLeadCurrent )
  {
    (extData.nextStoVectorRawPtr)[li_store_dev_i] = Id;
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 diode instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

// New Spice3f5 matched conductivities
// Only major difference in loads is support for Pos=Pri node when RS=0
// For this DAE dQdx load, the capacitance contribution (Gcd) is the
// only part used.
//---------------------------------------------------------------------

  dQdx[li_Neg][ANegEquNegNodeOffset] += Cd;
  dQdx[li_Neg][ANegEquPriNodeOffset] -= Cd;

  dQdx[li_Pri][APriEquNegNodeOffset] -= Cd;
  dQdx[li_Pri][APriEquPriNodeOffset] += Cd;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 diode instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

// New Spice3f5 matched conductivities
// Only major difference in loads is support for Pos=Pri node when RS=0
// For this DAE dFdx load, the capacitance contribution (Gcd) is removed.
//---------------------------------------------------------------------

  dFdx[li_Pos][APosEquPosNodeOffset] += Gspr;
  dFdx[li_Pos][APosEquPriNodeOffset] -= Gspr;

  dFdx[li_Neg][ANegEquNegNodeOffset] += Gd;
  dFdx[li_Neg][ANegEquPriNodeOffset] -= Gd;

  dFdx[li_Pri][APriEquPosNodeOffset] -= Gspr;
  dFdx[li_Pri][APriEquNegNodeOffset] -= Gd;
  dFdx[li_Pri][APriEquPriNodeOffset] += Gspr + Gd;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       : update primary state for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  double * stoVec = extData.nextStoVectorRawPtr;
  stoVec[li_storevd] = Vd;

  //Qd - capacitor charge generated in updateIntermediateVars
  bool bsuccess = updateIntermediateVars ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  double * solVec = extData.nextSolVectorRawPtr;

  //diode parameters
  double M;       // grading parameter
  double BV;      // breakdown voltage
  double IBV;     // reverse breakdown current
  double NBV;     // reverse breakdown ideality factor
  double IBVL;    // low-level reverse breakdown current
  double NBVL;    // low-level reverse breakdown ideality factor
  double N;       // non-ideality factor.
  double NR;      // emission coeff. for ISR.
  double TT;      // transit time.
  double F2;      // capacitive polynomial factor
  double F3;      // capacitive polynomial factor
  double Inorm;   // normal diffusion current
  double Irec;    // recombination current
  double Kgen;    // generation factor
  double Khi;     // high-injection factor
  double Gd1, DKgen;
  double Gd2, DKhi;

  //Model Diode parameters
  M    = model_.M;
  BV   = model_.BV;
  IBV  = model_.IBV;
  NBV  = model_.NBV;
  IBVL = model_.IBVL;
  NBVL = model_.NBVL;
  N    = model_.N;
  NR   = model_.NR;
  TT   = model_.TT;
  F2   = model_.F2;
  F3   = model_.F3;

  // obtain voltage drop accross the capacitor:
  Vp  = Vn = Vpp = 0.0;
  Vpp = solVec[li_Pri];
  Vn  = solVec[li_Neg];
  Vp  = solVec[li_Pos];

  // Junction Voltage
  Vd = Vpp - Vn;

  double Isat  = tSatCur * Area;
  double IsatR = tSatCurR * Area;
  double Vt    = CONSTKoverQ * Temp;
  double Vte   = N * Vt;
  double VteR  = NR * Vt;

  Gspr = tCOND * Area;

  Vd_orig = Vd;
  origFlag = true;

  // Setup initial junction conditions if UIC enabled
  //------------------------------------------------
  if (getSolverState().newtonIter == 0)
  {
    newtonIterOld = 0;
    //Vd_old = Vd;
    if (getSolverState().initJctFlag && getDeviceOptions().voltageLimiterFlag)
    {
      if (InitCondGiven)
      {
        Vd = InitCond;
        origFlag = false;
      }
      else if (off)
      {
        Vd = 0.0;
        origFlag = false;
      }
      else
      {
        if (getSolverState().inputOPFlag)
        {
          N_LAS_Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
          if ((*flagSolVectorPtr)[li_Pos] == 0 ||
              (*flagSolVectorPtr)[li_Neg] == 0 ||
              (*flagSolVectorPtr)[li_Pri] == 0)
          {
            Vd=tVcrit;
            origFlag = false;
          }
        }
        else
        {
          Vd=tVcrit;
          origFlag = false;
        }
      }
    }

    // assume there is no history -- then check if the
    // state vector can overwrite this
    Vd_old = Vd;

    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      Vd_old = (extData.currStoVectorRawPtr)[li_storevd];
    }
  }
  else  // just do this whenever it isn't the first iteration:
  {
    Vd_old = (extData.nextStoVectorRawPtr)[li_storevd];
  }

  // Voltage limiting based on mode of diode
  //---------------------------------------
  if (getDeviceOptions().voltageLimiterFlag)
  {
    int ichk = 0;

    if (getSolverState().newtonIter >= 0)
    {
      //Test if breakdown voltage given or not
      if (model_.BVGiven && (Vd < Xycemin(0.0, -BV + 10.0 * Vte)))
      {
        double Vdtmp = -( BV + Vd );
        Vdtmp = devSupport.pnjlim(Vdtmp, -(Vd_old+BV), Vte, tVcrit, &ichk);
        Vd    = -(Vdtmp + BV);
      }
      else
        Vd = devSupport.pnjlim(Vd, Vd_old, Vte, tVcrit, &ichk);

      if (ichk) origFlag = false;
    }
  }

  // update the "old" newton iteration number.
  if (getSolverState().newtonIter != 0 && getSolverState().newtonIter != newtonIterOld)
  {
    newtonIterOld = getSolverState().newtonIter;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
  {
    cout << "----------------------------------------------" << endl;
    cout << "Instance::updateIntermediateVars " << getName()<<endl;
  }
#endif

  // Current and Conductivity
  //----------------------------------------------

  int level = model_.getLevel();
  if(level == 1)
  {
    // Using LambertW function
    if (lambertWFlag == 1)
    {
      double RS = model_.RS;
      if (Vd >= -3.0 * Vte)
      {
        lambertWCurrent(Isat, Vte, RS);
      }
      // linear reverse bias
      else if ( !tBrkdwnV || (Vd >= -tBrkdwnV) )
      {
        lambertWLinearReverseBias(Isat, Vte, RS);
      }
      // reverse breakdown
      else
      {
        lambertWBreakdownCurrent(Isat, Vte, RS);
      }

      double Vrs;
      Vrs = (Id + Icd)*RS;
      Vc = Vd - Vrs;
    }
    else if (lambertWFlag == 2)
    {
      if (Vd >= -3.0 * Vte)
      {
        lambertWCurrent(Isat, Vte, 1.0e-15);
      }
      // linear reverse bias
      else if ( !tBrkdwnV || (Vd >= -tBrkdwnV) )
      {
        lambertWLinearReverseBias(Isat, Vte, 1.0e-15);
      }
      // reverse breakdown
      else
      {
        lambertWBreakdownCurrent(Isat, Vte, 1.0e-15);
      }
      Vc = Vd;
    }

    // Normal exponential
    else
    {
      // Adjustment for linear portion of reverse current
      double IRfactor;
      //if(Vd >= 0) IRfactor = 1.0;  //error?  shouldn't this be Vd>=-3*Vte????
      if(Vd >= -3.0 * Vte) IRfactor = 1.0;  //changing this to be consistent
      else IRfactor = tIRF;                 //with model in reference guide.
      Isat *= IRfactor;

      if (Vd >= -3.0 * Vte)
      {
        double arg1 = Vd / Vte;
        arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
        double evd = exp(arg1);

        Id = Isat * (evd - 1.0) + getDeviceOptions().gmin * Vd;
        Gd = Isat * evd / Vte + getDeviceOptions().gmin;
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
        {
          cout  << "Normal exponential regime." << endl;
          cout  << " Vd  = " << Vd << endl;
          cout  << " Vte = " << Vte << endl;
          cout  << " Id  = " << Id  << endl;
          cout  << " Gd  = " << Gd  << endl;
        }
#endif

      }

      // Linear reverse bias
      else if(!tBrkdwnV || (Vd >= -tBrkdwnV))
      {
        double arg = 3.0 * Vte / (Vd * CONSTe);
        arg = arg * arg * arg;
        Id = -Isat * (1.0 + arg) + getDeviceOptions().gmin * Vd;
        Gd = Isat * 3.0 * arg / Vd + getDeviceOptions().gmin;
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
        {
          cout  << "Linear reverse bias regime." << endl;
          cout  << " Vd       = " << Vd << endl;
          cout  << " tBrkdwnV = " << tBrkdwnV << endl;
          cout  << " Id       = " << Id  << endl;
          cout  << " Gd       = " << Gd  << endl;
        }
#endif
      }

      // Reverse breakdown
      else
      {
        double arg1 = -(tBrkdwnV + Vd) / Vte;
        arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
        double evrev = exp(arg1);

#ifdef Xyce_BREAKDOWN_ORIGINAL
        Id = -Isat * evrev + getDeviceOptions().gmin * Vd;
        Gd = Isat * evrev / Vte + getDeviceOptions().gmin;
#else
        //added by K. Santarelli 9/18/07 to account for change in tBrkdwnV
        //calculation.
        double arg2=3.0*Vte/(CONSTe*tBrkdwnV);
        arg2=arg2*arg2*arg2;
        double Isat_tBrkdwnV=Isat*(1-arg2);
        Id = -Isat_tBrkdwnV * evrev + getDeviceOptions().gmin * Vd;
        Gd = Isat_tBrkdwnV * evrev / Vte + getDeviceOptions().gmin;
#endif

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
        {
          cout  << "Reverse breakdown regime." << endl;
          cout  << " Vd       = " << Vd << endl;
          cout  << " tBrkdwnV = " << tBrkdwnV << endl;
          cout  << " Id       = " << Id  << endl;
          cout  << " Gd       = " << Gd  << endl;
        }
#endif
      }
      Vc = Vd;
    }

  }
  else if(level == 2)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
    {
      cout << " Level 2 diode code " << endl;
    }
#endif
    // Adjustment for linear portion of reverse current
    double IRfactor;
    //if(Vd >= 0) IRfactor = 1.0;    //error?  Shouldn't this be Vd >= -3*Vte?
    if(Vd >= -3.0*Vte) IRfactor=1.0; //changing it to be consistent with model
    else IRfactor = tIRF;            //in the reference manual.
    Isat *= IRfactor;

    if (Vd >= -3.0 * Vte)
    {
      double arg1 = Vd / Vte;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double evd = exp(arg1);
      Inorm = Isat * (evd - 1.0) + getDeviceOptions().gmin * Vd;
      Gd1 = Isat*evd/Vte + getDeviceOptions().gmin;

      arg1 = Vd / VteR;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      evd  = exp(arg1);
      Irec = IsatR * (evd - 1.0);
      Gd2  = IsatR*evd/VteR;

      Khi = 1;
      DKhi = 0;
      if(tIKF > 0)
      {
        Khi = sqrt(tIKF/(tIKF+Inorm));
        DKhi = 0.5*Khi*Gd1/(tIKF+Inorm);
      }
      Kgen = 0;
      DKgen = 0;
      if(Irec != 0)
      {
        Kgen = sqrt( pow(((1-Vd/tJctPot)*(1-Vd/tJctPot) + 0.005),M) );
        DKgen = -M*(1-Vd/tJctPot)/(tJctPot*Kgen);
      }

      Id = Inorm*Khi + Irec*Kgen;
      Gd = Gd1*Khi + Inorm*DKhi + Gd2*Kgen + Irec*DKgen;
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
      {
        cout  << "L2 Normal exponential regime." << endl;
        cout  << " Vd  = " << Vd << endl;
        cout  << " Vte = " << Vte << endl;
        cout  << " Id  = " << Id  << endl;
        cout  << " Irec= " << Irec  << endl;
        cout  << " Gd  = " << Gd  << endl;
        cout  << " Gd1 = " << Gd1 << endl;
        cout  << " Gd2 = " << Gd2 << endl;
        cout  << " Khi = " << Khi << endl;
        cout  << " Kgen=" << Kgen << endl;
        cout  << "DKgen=" <<DKgen << endl;
      }
#endif
    }

    // Linear reverse bias
    else if(!tBrkdwnV || (Vd >= -tBrkdwnV))
    {
      double arg = 3.0 * Vte / (Vd * CONSTe);
      arg = arg * arg * arg;
      Id = -Isat * (1.0 + arg) + getDeviceOptions().gmin * Vd;
      Gd = Isat * 3.0 * arg / Vd + getDeviceOptions().gmin;
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
      {
        cout  << "L2 Linear reverse bias regime." << endl;
        cout  << " Vd       = " << Vd << endl;
        cout  << " tBrkdwnV = " << tBrkdwnV << endl;
        cout  << " Id       = " << Id  << endl;
        cout  << " Gd       = " << Gd  << endl;
      }
#endif
    }

    // Reverse breakdown
    else
    {
      double arg1 = -(tBrkdwnV + Vd) / Vte;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double evrev = exp(arg1);

#ifdef Xyce_BREAKDOWN_ORIGINAL
      Id = -Isat * evrev + getDeviceOptions().gmin * Vd;
      Gd = Isat * evrev / Vte + getDeviceOptions().gmin;
#else
      //added 9/18/07 by K. Santarelli to account for change in tBrkdwnV
      //calculation.
      double arg2=3.0*Vte/(CONSTe*tBrkdwnV);
      arg2=arg2*arg2*arg2;
      double Isat_tBrkdwnV=Isat*(1-arg2);
      Id = -Isat_tBrkdwnV * evrev + getDeviceOptions().gmin * Vd;
      Gd = Isat_tBrkdwnV * evrev / Vte + getDeviceOptions().gmin;
#endif
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
      {
        cout  << "L2 Reverse breakdown regime." << endl;
        cout  << " Vd       = " << Vd << endl;
        cout  << " tBrkdwnV = " << tBrkdwnV << endl;
        cout  << " Id       = " << Id  << endl;
        cout  << " Gd       = " << Gd  << endl;
      }
#endif
    }
    Vc = Vd;

  }  // level

  // Only compute if Capacitance is non-zero
  //---------------------------------------
  if (tJctCap != 0.0)
  {
    // Charge Storage
    double Czero = tJctCap * Area;
    if (Vc < tDepCap)
    {
      //double arg = 1.0 - Vd/VJ;
      double arg = 1.0 - Vc / tJctPot;
      double arg1 = -M * log(arg);
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double sarg = exp(arg1);

      //Qd = TT*Id + VJ*Czero*(1.0-arg*sarg)/(1.0-M);
      Qd = TT * Id + tJctPot * Czero * (1.0 - arg * sarg) / (1.0 - M);
      Cd = TT * Gd + Czero * sarg;
    }
    else
    {
      double Czof2 = Czero / F2;
      double MotJctPot = M / tJctPot;

      //Qd = TT*Id + Czero*tF1 + Czof2*(F3*(Vd-tDepCap)+(M/(VJ+VJ))
      Qd = TT * Id + Czero * tF1 +
           Czof2 * (F3 * (Vc - tDepCap) + (0.5 * MotJctPot) *
                    (Vc * Vc - tDepCap * tDepCap));
      //Cd = TT*Gd + Czof2*(F3+M*Vd/VJ);
      Cd = TT * Gd + Czof2 * (F3 + MotJctPot * Vc);
    }
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
    {
      cout << "Qd = " << Qd << endl;
      cout << "Cd = " << Qd << endl;
    }
#endif
  }
  else
  {
    Qd = 0.0;
    Cd = 0.0;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
  {
    cout << "----------------------------------------------" << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature( const double & temp )
{

  double vtnom = CONSTKoverQ * model_.TNOM;

  double xfc = log( 1.0 - model_.FC );

  if( temp != -999.0 ) Temp = temp;
  double TNOM = model_.TNOM;

  double vt = CONSTKoverQ * Temp;
  double fact2 = Temp / CONSTREFTEMP;
  double egfet = CONSTEg0 - (CONSTalphaEg*Temp*Temp)/(Temp+CONSTbetaEg);
  double arg = -egfet/(2.0*CONSTboltz*Temp) +
               CONSTEg300/(CONSTboltz*(CONSTREFTEMP+CONSTREFTEMP));
  double pbfact = -2.0*vt*(1.5*log(fact2)+CONSTQ*arg);
  double egfet1 = CONSTEg0 - (CONSTalphaEg*model_.TNOM*model_.TNOM)/
                  (model_.TNOM+CONSTbetaEg);
  double arg1 = -egfet1/(2.0*CONSTboltz*model_.TNOM) +
		CONSTEg300/(2.0*CONSTboltz*CONSTREFTEMP);
  double fact1 = model_.TNOM/CONSTREFTEMP;
  double pbfact1 = -2.0*vtnom*(1.5*log(fact1)+CONSTQ*arg1);

  double pbo = (model_.VJ-pbfact1)/fact1;
  double gmaold = (model_.VJ-pbo)/pbo;

  tJctCap = model_.CJO/
            (1.0+model_.M*(4.0e-4*(model_.TNOM-CONSTREFTEMP) -gmaold));

  tJctPot = pbfact+fact2*pbo;

  double gmanew = (tJctPot-pbo)/pbo;

  tJctCap *= 1.0 + model_.M*(4.0e-4*(Temp-CONSTREFTEMP)-gmanew);

  tSatCur = model_.IS*exp(((Temp/model_.TNOM)-1.0)*
                          model_.EG/(model_.N*vt)+
                          model_.XTI/model_.N*log(Temp/model_.TNOM));

  tF1 = tJctPot*(1.0-exp((1.0-model_.M)*xfc))/(1.0-model_.M);

  tDepCap = model_.FC*tJctPot;

  double vte = model_.N*vt;
  double tempBV;

  tVcrit = vte*log(vte/(CONSTroot2*tSatCur));
  tRS   = model_.RS;
  tCOND = model_.COND;
  tIRF  = model_.IRF*pow(fact2,1.6);

  int level = model_.getLevel();
  if(level == 2)   // this section is PSPICE compatible
  {
    tSatCurR = model_.ISR*exp((Temp/TNOM - 1.0)*
                              model_.EG/(model_.NR*vt)+
                              model_.XTI/model_.NR*log(Temp/TNOM));

    tIKF = model_.IKF*(1 + model_.TIKF*(Temp-TNOM));

    tempBV = model_.BV*(1 + (Temp-TNOM)*
                        ( model_.TBV1 + model_.TBV2*(Temp-TNOM) ));

    tRS = model_.RS*(1 + (Temp-TNOM)*
                     ( model_.TRS1 + model_.TRS2*(Temp-TNOM) ));

    tCOND = 0.0;
    if(tRS != 0.0) tCOND = 1.0/tRS;

    tJctPot = (model_.VJ - egfet1)*fact2 - 3*vt*log(fact2) + egfet;

    tJctCap = model_.CJO/(1.0 +
                          model_.M*(4.0e-4*(Temp-TNOM) + (1-tJctPot/model_.VJ)));
  }
  else
  {
    tempBV=model_.BV;
  }

#ifdef Xyce_BREAKDOWN_ORIGINAL
  //Changed 9/18/07, K. R. Santarelli.  This is the original (broken) version
  //of the breakdown voltage computation.  It has two known issues:
  //
  //1.  It assumes N (the emission coefficient) is always 1 (division by
  //    vt instead of vte).
  //2.  While the for loop terminates due to the fact that it's implemented via
  //    a counter, the value of xbv does not converge in many cases, so the
  //    value of xbv upon terminating is often garbage.
  //
  //For consistency with old, existing simulations, this version of the
  //breakdown voltage calculation (along with the appropriate code for
  //computing the breakdown current in the level 1 and level 2 models---see
  //the appropriate sections of
  //Instance::updateIntermediateVars) can be invoked by
  //specifiying the CPPFLAG "-DXyce_BREAKDOWN_ORIGINAL" when creating a Xyce
  //Build.  The default, however, is the code which follows #else.

  double reltol = 1.0e-3;
  if( model_.BVGiven )
  {
    double IRfactor = tIRF;
    double cbv = model_.IBV;
    double xbv, xcbv;
    if( cbv < IRfactor*tSatCur*tempBV/vt )
    {
      cbv = IRfactor*tSatCur*tempBV/vt;
      xbv = tempBV;
    }
    else
    {
      double tol = reltol*cbv;
      xbv = tempBV-vt*log(1.0+cbv/(IRfactor*tSatCur));
      for( int i = 0; i < 25; ++i )
      {
        xbv = tempBV-vt*log(cbv/(IRfactor*tSatCur)+1.0-xbv/vt);
        xcbv = IRfactor*tSatCur*(exp((tempBV-xbv)/vt)-1.0+xbv/vt);
        if(fabs(xcbv-cbv)<=tol) break;
      }
    }
    tBrkdwnV = xbv;
  }
#else
  if( model_.BVGiven)
  {
    double IRFactor=tIRF;
    double cbv = model_.IBV;
    double xbv;
    double cthreshlow; //lower threshold for IBV
    double cthreshhigh; //high threshold for IBV
    int iter_count;
    const int ITER_COUNT_MAX=8;
    double arg2;

    double arg1=3.0*vte/(CONSTe*tempBV);
    arg1=arg1*arg1*arg1;
    cthreshlow=tSatCur*IRFactor*(1-arg1);
    cthreshhigh=tSatCur*IRFactor*(1-1.0/(CONSTe*CONSTe*CONSTe)) *
                exp(-1.0*(3.0*vte-tempBV)/vte);

    if(cbv >= cthreshhigh)
    {                     //if IBV is too high, tBrkdwnV will go below 3NVt.
      tBrkdwnV=3.0*vte;   //Clip tBrkdwnV to 3*N*Vt in this case (and hence
    }                     //clip IBV to cthreshhigh).


    else if(cbv <= cthreshlow)
    {                          //if IBV is too low, tBrkdwnV will go above
      tBrkdwnV=tempBV;  //BV.  Clip tBrkdwnV to BV in this case (and
    }                          //hence clip IBV to cthreshlow).


    //If IBV is in an acceptable range, perform a Picard iteration to find
    //tBrkdwnV, starting with an initial guess of tBrkdwnV=tempBV, and
    //running through the iteration ITER_COUNT_MAX times.

    else
    {
      xbv=tempBV;
      for(iter_count=0; iter_count < ITER_COUNT_MAX; iter_count++)
      {
        arg2=3.0*vte/(CONSTe*xbv);
        arg2=arg2*arg2*arg2;
        xbv=tempBV-vte*log(cbv/(tSatCur*IRFactor))+vte *
            log(1-arg2);
      }
      tBrkdwnV=xbv;
    }

//Note that, not only is tBrkdwnV adjusted using this method, but the effective
//value of Is (tSatCur) is adjusted as well.  The code used to use Is before,
//but now it will use Is*IRFactor*(1-(3*N*Vt/(e*tBrkdwnV))^3) instead.  This is
//reflected in changes made to Instance::updateIntermediateVars
//(for "normal" level 1 and 2 diode models) for the reverse breakdown region.
//Changes have not been made to the associated LambertW functions as of yet
//since there seem to be other issues involved with those functions independent
// of this change.
  }
#endif


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
  {
    cout << "----------------------------------------------" << endl;
    cout << "Instance::UpdateTemperature" << getName() <<endl;
    cout << " IS = " << model_.IS << endl;
    cout << " vtnom   = " << vtnom << endl;
    cout << " xfc     = " << xfc << endl;
    cout << " TNOM    = " << TNOM << endl;
    cout << " vt      = " << vt << endl;
    cout << " fact2   = " << fact2 << endl;
    cout << " egfet   = " << egfet << endl;
    cout << " arg     = " << arg   << endl;
    cout << " pbfact  = " << pbfact   << endl;
    cout << " egfet1  = " << egfet1   << endl;
    cout << " arg1    = " << arg1     << endl;
    cout << " fact1   = " << fact1    << endl;
    cout << " pbfact1 = " << pbfact1  << endl;
    cout << " pbo     = " << pbo      << endl;
    cout << " gmaold  = " << gmaold   << endl;
    cout << " gmanew  = " << gmanew   << endl;
    cout << " tJctCap = " << tJctCap  << endl;
    cout << " tJctPot = " << tJctPot  << endl;
    cout << " tSatCur = " << tSatCur  << endl;
    cout << " tF1     = " << tF1      << endl;
    cout << " tDepCap = " << tDepCap  << endl;
    cout << " vte     = " << vte      << endl;
    cout << " tempBV  = " << tempBV   << endl;
    cout << " tVcrit  = " << tVcrit   << endl;
    cout << " tRS     = " << tRS      << endl;
    cout << " tCOND   = " << tCOND    << endl;
    cout << " tIRF    = " << tIRF     << endl;
    cout << " tBrkdwnV= " << tBrkdwnV << endl;
  }
#endif

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::lambertWCurrent
// Purpose       : Determine the diode current using Lambert-W function
// Special Notes :
// Scope         : public
// Creator       : Nick Johnson, Summer Intern
// Creation Date : 7/5/02
//-----------------------------------------------------------------------------
bool Instance::lambertWCurrent(double Isat, double Vte, double RS)
{
  // with capacitor current and using LambertW accros whole model
  /*  if (Cd != 0.0 && Icd != 0 && (lambertWFlag == 1))
      {
      double AA = Vte*(1.0 + RS*getSolverState().pdt*Cd);
      double XX = Icd - getSolverState().pdt*Qd;
      double arg1 = (Vd + RS*Isat - RS*XX)/AA;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double evd = exp(arg1);
      double lambWArg = Isat*RS/AA * evd;
      double lambWReturn;
      int ierr;
      double lambWError;
      devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

      #ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
      if (ierr == 0)
      cout << "Safe LambertW return" << endl;
      else if (ierr == 1)
      cout << "LambertW argument not in domain" << endl;
      else
      cout << "Arithmetic problems with LambertW" << endl;
      }
      #endif

      double prefac = AA/RS;
      Id = -Isat + prefac*lambWReturn;
      Gd = lambWReturn / ((1 + lambWReturn)*RS);

      #ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
      cout << "lambWArg   = " << lambWArg << endl;
      cout << "lambWError = " << lambWError << endl;
      cout << "Id         = " << Id << endl;
      cout << "Gd         = " << Gd << endl;
      cout << "Icd        = " << Icd << endl;
      cout << "Cd         = " << Cd  << endl;
      cout << "Qd         = " << Qd  << endl;
      cout << "AA         = " << AA << endl;
      cout << "XX         = " << XX << endl;
      cout << "Using lambertwCurrent w/ capacitance" << endl;
      }
      #endif
      }

      // when capacitor current=0, there is no capacitance, or using LambertW
      // only across the diode element
      else
      { */
  double arg1 = (Vd + Isat*RS)/Vte;
  arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
  double evd = exp(arg1);
  double lambWArg = Isat*RS*evd/Vte;
  double lambWReturn;
  int ierr;
  double lambWError;
  devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

  Id = -Isat+Vte*(lambWReturn)/RS;
  Gd = lambWReturn / ((1 + lambWReturn)*RS);

  // }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::lambertWLinearReverseBias
// Purpose       : Function to maintain continuity between forward bias and
//                 breakdown voltages of diode
// Special Notes :
// Scope         : public
// Creator       : Nick Johnson, Summer Intern
// Creation Date : 7/30/02
//-----------------------------------------------------------------------------
bool Instance::lambertWLinearReverseBias(double Isat, double Vte, double RS)
{
  double FF1 = (Vd + tBrkdwnV)/(-3.0 * Vte + tBrkdwnV);
  double FF2 = (1/2 - FF1)*(-2);
  double arg = Vte/RS;
  double arg1 = (Isat/arg - 3);
  double arg2 = FF1*arg1;
  arg1 = Xycemin(CONSTMAX_EXP_ARG, arg2);
  double evd = exp(arg2);
  double lambWArg = Isat*evd/arg;
  double lambWReturn;
  int ierr;
  double lambWError;
  devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    if (ierr == 0)
      cout << "Safe LambertW return" << endl;
    else if (ierr == 1)
      cout << "LambertW argument not in domain" << endl;
    else
      cout << "Arithmetic problems with LambertW" << endl;
  }
#endif

  Id = -Isat*FF1 + FF2*lambWReturn*arg;

  double GdFF1 = 1/(-3.0*Vte + tBrkdwnV);
  double GdFF2 = 2 * GdFF1;
  double GdW = arg*lambWReturn*GdFF1*arg1/(1 + lambWReturn);
  Gd = -Isat*GdFF1 + GdFF2*arg*lambWReturn + FF2*GdW;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "lambWArg   = " << lambWArg << endl;
    cout << "lambWError = " << lambWError << endl;
    cout << "lambWReturn= " << lambWReturn << endl;
    cout << "Id         = " << Id << endl;
    cout << "Gd         = " << Gd << endl;
    cout << "Using lambertwReverseBias" << endl;
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::lambertWBreakdownCurrent
// Purpose       : Determine the diode breakdown current using Lambert-W function
// Special Notes :
// Scope         : public
// Creator       : Nick Johnson, Summer Intern
// Creation Date : 7/11/02
//-----------------------------------------------------------------------------
bool Instance::lambertWBreakdownCurrent(double Isat, double Vte, double RS)
{
  // with capacitor current and applying LambertW accros whole diode model
  /*  if (Cd != 0.0 && Icd != 0 && (lambertWFlag == 1))
      {
      double AA = Vte*(1 + RS*getSolverState().pdt*Cd);
      double XX = Icd - getSolverState().pdt*Qd;
      double arg1 = (RS*XX - Vd)/AA - tBrkdwnV/Vte;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double evd = exp(arg1);
      double lambWArg = Isat*RS*evd/AA;
      double lambWReturn;
      int ierr;
      double lambWError;
      devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

      #ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
      if (ierr == 0)
      cout << "Safe LambertW return" << endl;
      else if (ierr == 1)
      cout << "LambertW argument not in domain" << endl;
      else
      cout << "Arithmetic problems with LambertW" << endl;
      }
      #endif

      Id = -AA*lambWReturn/RS;
      Gd = (lambWReturn / (1 + lambWReturn)) * (1/RS);

      #ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
      cout << "lambWArg   = " << lambWArg << endl;
      cout << "lambWError = " << lambWError << endl;
      cout << "Id         = " << Id << endl;
      cout << "Gd         = " << Gd << endl;
      cout << "Icd        = " << Icd << endl;
      cout << "Cd         = " << Cd  << endl;
      cout << "Qd         = " << Qd  << endl;
      cout << "AA         = " << AA << endl;
      cout << "XX         = " << XX << endl;
      cout << "Using lambertwBreakdown w/ capacitance" << endl;
      }
      #endif
      }

      // without capacitor current or applying lambertW just accros diode element
      else
      {*/
  double arg1 = (-Vd - tBrkdwnV)/ Vte;
  arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
  double evd = exp(arg1);
  double lambWArg = Isat*RS*evd/Vte;
  double lambWReturn;
  int ierr;
  double lambWError;
  devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

  Id = -Vte*lambWReturn/RS;
  Gd = lambWReturn / ((1 + lambWReturn)*RS);

  // }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{
  //limit grading coeff
  if( M > 0.9 ) M = 0.9;

  //limit activation energy
  if( EG < 0.1 ) EG = 0.1;

  //limit depl cap coeff
  if( FC > 0.95 ) FC = 0.95;

  if( RS==0.0 )
    COND = 0.0;
  else
    COND = 1.0/RS;

  double xfc = log(1.0-FC);
  F2 = exp((1.0+M)*xfc);
  F3 = 1.0-FC*(1.0+M);

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
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::Model(const ModelBlock & MB,
             SolverState & ss1,
             DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
    IS(1.0e-14),
    RS(0.0),
    COND(0.0),
    N(1.0),
    ISR(0.0),
    NR(2.0),
    IKF(0.0),
    TT(0.0),
    CJO(0.0),
    VJ(1.0),
    M(0.5),
    EG(1.11),
    XTI(3.0),
    TIKF(0.0),
    TBV1(0.0),
    TBV2(0.0),
    TRS1(0.0),
    TRS2(0.0),
    FC(0.5),
    BV(1e99),
    IBV(1.0e-10),
    IRF(1.0),
    NBV(1.0),
    IBVL(0.0),
    NBVL(1.0),
    F2(0.0),
    F3(0.0),
    TNOM(27),
    KF(0.0),
    AF(1.0),
    BVGiven(false)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    TNOM = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  // Note: Level 2 only params are: ISR, NR, IKF, NBV, IBVL, NBVL, TIKF, TBV1, TBV2, TRS1, TRS2
  if (getLevel() == 1)
  {
    string msg = "Illegal parameter(s) given for level 1 diode:";
    int okSize = msg.size();
    if (given("ISR"))
      msg += " ISR";
    if (given("NR"))
      msg += " NR";
    if (given("IKF"))
      msg += " IKF";
    if (given("NBV"))
      msg += " NBV";
    if (given("IBVL"))
      msg += " IBVL";
    if (given("NBVL"))
      msg += " NBVL";
    if (given("TIKF"))
      msg += " TIKF";
    if (given("TBV1"))
      msg += " TBV1";
    if (given("TBV2"))
      msg += " TBV2";
    if (given("TRS1"))
      msg += " TRS1";
    if (given("TRS2"))
      msg += " TRS2";
    if (msg.size() > okSize)
    {
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
  }

  processParams ();
}


//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::~Model()
{
  vector<Instance*>::iterator iterI;
  vector<Instance*>::iterator firstI = instanceContainer.begin ();
  vector<Instance*>::iterator lastI  = instanceContainer.end ();

  // loop over instances:
  for (iterI = firstI; iterI != lastI; ++iterI)
  {
    delete (*iterI);
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
  os << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << (*iter)->getModelName();

    os << endl;
    os << "AREA   = " << (*iter)->Area << endl;
    os << "  IC   = " << (*iter)->InitCond   << endl;
    os << "TEMP   = " << (*iter)->Temp << endl;
    os << " off   = " << (*iter)->off  << endl;

    os << endl;
  }

  os << endl;

  return os;
}

//-----------------------------------------------------------------------------
// Diode Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);

    // save voltage drops
    double * stoVec = di.extData.nextStoVectorRawPtr;
    stoVec[di.li_storevd] = di.Vd;

    bool btmp = di.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);
    // load F:
    double Ir = di.Gspr * (di.Vp - di.Vpp);

    fVec[di.li_Pos] -= -Ir;
    fVec[di.li_Neg] -= di.Id;
    fVec[di.li_Pri] -= (-di.Id + Ir);

    // load Q:
    qVec[di.li_Neg] -= di.Qd;
    qVec[di.li_Pri] -= -di.Qd;

    // voltage limiter vectors.
    if( getDeviceOptions().voltageLimiterFlag )
    {
      double Vd_diff = di.Vd - di.Vd_orig;
      double Cd_Jdxp = -( di.Cd ) * Vd_diff;
      double Gd_Jdxp = -( di.Gd ) * Vd_diff;

      double * dFdxdVp = di.extData.dFdxdVpVectorRawPtr;
      // dFdxdVp vector
      dFdxdVp[di.li_Neg] += Gd_Jdxp;
      dFdxdVp[di.li_Pri] -= Gd_Jdxp;

      double * dQdxdVp = di.extData.dQdxdVpVectorRawPtr;
      // dQdxdVp vector
      dQdxdVp[di.li_Neg] += Cd_Jdxp;
      dQdxdVp[di.li_Pri] -= Cd_Jdxp;
    }

    if( di.loadLeadCurrent )
    {
      storeLeadF[di.li_store_dev_i] = di.Id;
      if (di.model_.CJO != 0.0)
      {
        storeLeadQ[di.li_store_dev_i] = di.Qd;
      }
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // F-matrix:
    *di.fPosEquPosNodePtr += di.Gspr;
    *di.fPosEquPriNodePtr -= di.Gspr;
    *di.fNegEquNegNodePtr += di.Gd;
    *di.fNegEquPriNodePtr -= di.Gd;
    *di.fPriEquPosNodePtr -= di.Gspr;
    *di.fPriEquNegNodePtr -= di.Gd;
    *di.fPriEquPriNodePtr += di.Gspr + di.Gd;

    // Q-matrix:
    *di.qNegEquNegNodePtr += di.Cd;
    *di.qNegEquPriNodePtr -= di.Cd;
    *di.qPriEquNegNodePtr -= di.Cd;
    *di.qPriEquPriNodePtr += di.Cd;
#else
    // F-matrix:
    dFdx[di.li_Pos][di.APosEquPosNodeOffset] += di.Gspr;
    dFdx[di.li_Pos][di.APosEquPriNodeOffset] -= di.Gspr;
    dFdx[di.li_Neg][di.ANegEquNegNodeOffset] += di.Gd;
    dFdx[di.li_Neg][di.ANegEquPriNodeOffset] -= di.Gd;
    dFdx[di.li_Pri][di.APriEquPosNodeOffset] -= di.Gspr;
    dFdx[di.li_Pri][di.APriEquNegNodeOffset] -= di.Gd;
    dFdx[di.li_Pri][di.APriEquPriNodeOffset] += di.Gspr + di.Gd;

    // Q-matrix:
    dQdx[di.li_Neg][di.ANegEquNegNodeOffset] += di.Cd;
    dQdx[di.li_Neg][di.ANegEquPriNodeOffset] -= di.Cd;
    dQdx[di.li_Pri][di.APriEquNegNodeOffset] -= di.Cd;
    dQdx[di.li_Pri][di.APriEquPriNodeOffset] += di.Cd;
#endif
  }
  return true;
}

  } // namespace Diode
} // namespace Device
} // namespace Xyce
