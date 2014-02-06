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

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_ACC.C,v $
//
// Purpose        : Provide a solution to the simple second order
//                  initial value problem:
//                  d^2x/dt^2 = a(t); x(0) = x0, dx/dt(0) = v0
//
// Special Notes  : Intended for use in FY07 CoilGun LDRD work
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 10/25/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_ACC.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<ACC::Instance>::ParametricData()
{
  // Set up configuration constants:
  setNumNodes(3);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(0);
  // no primary parameter, no setPrimaryParameter call


  // Set up double precision variables:
  addPar ("V0", 0.0, false, ParameterType::NO_DEP,
          &ACC::Instance::v0,
          NULL, U_MSM1, CAT_NONE, "Initial Velocity");

  addPar ("X0", 0.0, false, ParameterType::NO_DEP,
          &ACC::Instance::x0,
          NULL, U_METER, CAT_NONE, "Initial Position");
}

template<>
ParametricData<ACC::Model>::ParametricData()
{}

namespace ACC {

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
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
Instance::Instance(
  InstanceBlock &               IB,
  Model & Miter,
  MatrixLoadData &              mlData1,
  SolverState &                 ss1,
  ExternData  &                 ed1,
  DeviceOptions &               do1)
  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(Miter),
    v0(0.0),
    x0(0.0),
    li_Acc(-1),
    li_Velocity(-1),
    li_Position(-1),
    li_state_vel(-1),
    li_state_pos(-1),
    AVelEquAccNodeOffset(-1),
    AVelEquVelNodeOffset(-1),
    APosEquVelNodeOffset(-1),
    APosEquPosNodeOffset(-1)
{
  numIntVars   = 0;
  numExtVars   = 3;
  numStateVars = 2; // position and velocity saved in state so we can get
                    //derivs

  devConMap.resize(3);
  devConMap[0]=1;
  devConMap[1]=1;
  devConMap[2]=1;

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    //jacStamp[0].resize(0); // the row for acceleration has nothing.
    jacStamp[1].resize(2); // velocity row
    jacStamp[1][0]=0;      //velocity-acceleration
    jacStamp[1][1]=1;      //velocity-velocity
    jacStamp[2].resize(2);
    jacStamp[2][0]=1;      // position-velocity
    jacStamp[2][1]=2;      // position-position
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:

  setParams (IB.params);

  // Set any non-constant parameter defaults:


}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                             const vector<int> & extLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "-------------------------------------------------------------------------"
    "----";
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << dashedline << endl;
    cout << "  ACCInstance::registerLIDs" << endl;
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

  li_Acc = extLIDVec[0];
  li_Velocity = extLIDVec[1];
  li_Position = extLIDVec[2];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << "  li_Acc = " << li_Acc << endl;
    cout << "  li_Velocity = " << li_Velocity << endl;
    cout << "  li_Position = " << li_Position << endl;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
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
  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Set up offsets so we can store quantities that need to be
//                 differentiated.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const vector<int> & staLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  li_state_vel = staLIDVec[0];
  li_state_pos = staLIDVec[1];

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/01
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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  AVelEquAccNodeOffset = jacLIDVec[1][0];
  AVelEquVelNodeOffset = jacLIDVec[1][1];
  APosEquVelNodeOffset = jacLIDVec[2][0];
  APosEquPosNodeOffset = jacLIDVec[2][1];
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  double * solVec = extData.nextSolVectorRawPtr;

  position = solVec[li_Position];
  velocity = solVec[li_Velocity];
  acceleration = solVec[li_Acc];

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::updatePrimaryState" <<endl;
  }
#endif
  bool bsuccess = updateIntermediateVars ();
  double * staVec = extData.nextStaVectorRawPtr;

  // We need to save the position and velocity so we can differentiate them
  // to give what we *think* the velocity and acceleration are.
  staVec[li_state_pos] = position;
  staVec[li_state_vel] = velocity;

  // if this is the first newton step of the first time step
  // of the transient simulation, we need to enforce that the
  // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
  // compatibility.  ERK.

  if (!(getSolverState().dcopFlag) && (getSolverState().initTranFlag) && getSolverState().newtonIter==0)
  {
    double * currStaVec = extData.currStaVectorRawPtr;
    currStaVec[li_state_pos] = position;
    currStaVec[li_state_vel] = velocity;
  }

  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  double * staDerivVec = extData.nextStaDerivVectorRawPtr;

  xdot = staDerivVec[li_state_pos];
  vdot = staDerivVec[li_state_vel];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 ACC  instance.
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;

  // Load DAE F-vector

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Instance::loadDAEFVector" << endl;
    cout << "  name = " << getName() <<endl;
    cout << "  velocity = " << velocity << endl;
    cout << "  position = " << position << endl;
  }
#endif

  if (getSolverState().dcopFlag)
  {

    fVec[li_Velocity] += velocity-v0;

    fVec[li_Position] += position-x0;
  }
  else
  {

    fVec[li_Velocity] += -acceleration;

    fVec[li_Position] += -velocity;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 ACC  instance.
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;


  qVec[li_Velocity] += velocity;

  qVec[li_Position] += position;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 ACC  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if (getSolverState().dcopFlag)
  {

    dFdx[li_Velocity][AVelEquVelNodeOffset] += 1.0;

    dFdx[li_Position][APosEquPosNodeOffset] += 1.0;
  }
  else
  {

    dFdx[li_Velocity][AVelEquAccNodeOffset] += -1.0;

    dFdx[li_Position][APosEquVelNodeOffset] += -1.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx ()
//
// Purpose       : Loads the Q-vector contributions for a single
//                 ACC  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  if (!getSolverState().dcopFlag)
  {

    dQdx[li_Velocity][AVelEquVelNodeOffset] += 1.0;

    dQdx[li_Position][APosEquPosNodeOffset] += 1.0;
  }

  return true;
}

// These are all just placeholder functions, as there is nothing to the model
// class.
// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
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
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
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
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/25/07
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

} // namespace ACC
} // namespace Device
} // namespace Xyce
