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
// Filename       : $RCSfile: N_DEV_VCCS.C,v $
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
// Revision Number: $Revision: 1.110.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:39 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_VCCS.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {
template<>
ParametricData<VCCS::Instance>::ParametricData()
{
    setNumNodes(4);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(0);
    setPrimaryParameter("T");

    // Set up double precision variables:
    addPar ("T", 0.0, false, ParameterType::NO_DEP,
            &VCCS::Instance::Transconductance,
        NULL,U_NONE,CAT_NONE,"Transconductance");
}

template<>
ParametricData<VCCS::Model>::ParametricData()
{
}

namespace VCCS {

// static member components
vector< vector<int> > Instance::jacStamp;




ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : ::
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 12/20/00
//-----------------------------------------------------------------------------
Instance::Instance (
  InstanceBlock & IBref,
    Model & Viter,
    MatrixLoadData & mlData1,
    SolverState &ss1,
    ExternData  &ed1,
    DeviceOptions & do1)

  : DeviceInstance(IBref, mlData1, ss1, ed1, do1),
    IB(IBref),
    model_(Viter),
    Transconductance(1.0),
    li_Pos(-1),
    li_Neg(-1),
    li_ContPos(-1),
    li_ContNeg(-1),
    li_store_dev_i(-1),
    APosEquContPosVarOffset(-1),
    APosEquContNegVarOffset(-1),
    ANegEquContPosVarOffset(-1),
    ANegEquContNegVarOffset(-1),
    f_PosEquContPosVarPtr(0),
    f_PosEquContNegVarPtr(0),
    f_NegEquContPosVarPtr(0),
    f_NegEquContNegVarPtr(0)
{
  numIntVars   = 0;
  numExtVars   = 4;
  numStateVars = 0;
  numLeadCurrentStoreVars = 1; // lead current

  if( jacStamp.empty() )
  {
    jacStamp.resize(4);
    jacStamp[0].resize(2);
    jacStamp[0][0]=2;
    jacStamp[0][1]=3;
    jacStamp[1].resize(2);
    jacStamp[1][0]=2;
    jacStamp[1][1]=3;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("T"))
  {
    string msg = "Could not find Transconductance parameter in instance.";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

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

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Component Information and Models
// Creation Date : 6/21/01
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
	                               const vector<int> & extLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
      "-----------------------------------------------------------------------------";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << dashedline << endl;
    cout << "  VCCSInstance::registerLIDs" << endl;
    cout << "  name = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_ContPos = extLIDVec[2];
  li_ContNeg = extLIDVec[3];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  li_Pos = " << li_Pos << endl;
    cout << "  li_Neg = " << li_Neg << endl;
    cout << "  li_ContPos = " << li_ContPos << endl;
    cout << "  li_ContNeg = " << li_ContNeg << endl;

    cout << dashedline << endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoektra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef)
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSTa != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 4/23/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const vector<int> & stoLIDVecRef )
{
  int numSto = stoLIDVecRef.size();
  if (numSto != getNumStoreVars())
  {
    string msg = "Instance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;
  if( loadLeadCurrent )
  {
    li_store_dev_i = stoLIDVec[0];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 4/23/2013
//-----------------------------------------------------------------------------
map<int,string> & Instance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if( loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    string modName(getName());
    spiceInternalName(modName);
    string tmpstr;
    tmpstr = modName+":DEV_I";
    storeNameMap[ li_store_dev_i ] = tmpstr;
  }

  return storeNameMap;
}
//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoektra, SNL, Parallel Computational Sciences
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
// Creator       : Robert Hoektra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquContPosVarOffset = jacLIDVec[0][0];
  APosEquContNegVarOffset = jacLIDVec[0][1];
  ANegEquContPosVarOffset = jacLIDVec[1][0];
  ANegEquContNegVarOffset = jacLIDVec[1][1];
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  f_PosEquContPosVarPtr = &(dFdx[li_Pos][APosEquContPosVarOffset]);
  f_PosEquContNegVarPtr = &(dFdx[li_Pos][APosEquContNegVarOffset]);
  f_NegEquContPosVarPtr = &(dFdx[li_Neg][ANegEquContPosVarOffset]);
  f_NegEquContNegVarPtr = &(dFdx[li_Neg][ANegEquContNegVarOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 VCCS instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;

  double v_cont_pos = solVec[li_ContPos];
  double v_cont_neg = solVec[li_ContNeg];

  fVec[li_Pos] += Transconductance * ( v_cont_pos - v_cont_neg );
  fVec[li_Neg] += -Transconductance * ( v_cont_pos - v_cont_neg );

  if( loadLeadCurrent )
  {
    double * storeLeadF = extData.nextStoVectorRawPtr;
    storeLeadF[li_store_dev_i] = Transconductance * ( v_cont_pos - v_cont_neg );
  }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 VCCS instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquContPosVarOffset] += Transconductance;
  dFdx[li_Pos][APosEquContNegVarOffset] -= Transconductance;
  dFdx[li_Neg][ANegEquContPosVarOffset] -= Transconductance;
  dFdx[li_Neg][ANegEquContNegVarOffset] += Transconductance;

  return true;
}

// Class Model
//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 12/21/00
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
  os << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;

  return os;
}

// VCCS Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & vi = *(*it);
    double v_cont_pos = solVec[vi.li_ContPos];
    double v_cont_neg = solVec[vi.li_ContNeg];


    fVec[vi.li_Pos] += vi.Transconductance * ( v_cont_pos - v_cont_neg );
    fVec[vi.li_Neg] += -vi.Transconductance * ( v_cont_pos - v_cont_neg );
    if( vi.loadLeadCurrent )
    {
      storeLeadF[vi.li_store_dev_i] = vi.Transconductance * ( v_cont_pos - v_cont_neg );
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
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & vi = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *vi.f_PosEquContPosVarPtr += vi.Transconductance;
    *vi.f_PosEquContNegVarPtr -= vi.Transconductance;
    *vi.f_NegEquContPosVarPtr -= vi.Transconductance;
    *vi.f_NegEquContNegVarPtr += vi.Transconductance;
#else

    dFdx[vi.li_Pos][vi.APosEquContPosVarOffset] += vi.Transconductance;
    dFdx[vi.li_Pos][vi.APosEquContNegVarOffset] -= vi.Transconductance;
    dFdx[vi.li_Neg][vi.ANegEquContPosVarOffset] -= vi.Transconductance;
    dFdx[vi.li_Neg][vi.ANegEquContNegVarOffset] += vi.Transconductance;
#endif
  }

  return true;
}

} // namespace VCCS
} // namespace Device
} // namespace Xyce
