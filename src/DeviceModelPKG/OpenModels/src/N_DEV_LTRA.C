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
// Filename       : $RCSfile: N_DEV_LTRA.C,v $
//
// Purpose        : lossy transmission line
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 06/16/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.35.2.5 $
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
#include <N_DEV_LTRA.h>
#include <N_DEV_LTRA_Faddeeva.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_BreakPoint.h>

#include <N_UTL_Functors.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<LTRA::Instance>::ParametricData()
{
    setNumNodes(4);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(1);
    addModelType("LTRA");

    // Set up double precision variables:
    addPar ("V1", 0.0, false, ParameterType::NO_DEP,
        &LTRA::Instance::initVolt1,
        &LTRA::Instance::initVolt1Given,
        U_VOLT,CAT_NONE, "Initial voltage at end 1");

    addPar ("V2", 0.0, false, ParameterType::NO_DEP,
        &LTRA::Instance::initVolt2,
        &LTRA::Instance::initVolt2Given,
        U_VOLT,CAT_NONE, "Initial voltage at end 2");

    addPar ("I1", 0.0, false, ParameterType::NO_DEP,
        &LTRA::Instance::initCur1,
        &LTRA::Instance::initCur1Given,
        U_AMP,CAT_NONE, "Initial current at end 1");

    addPar ("I2", 0.0, false, ParameterType::NO_DEP,
        &LTRA::Instance::initCur2,
        &LTRA::Instance::initCur2Given,
        U_AMP,CAT_NONE, "Initial current at end 2");
}

template<>
ParametricData<LTRA::Model>::ParametricData()
{
    addPar ("R", 0.0, false, ParameterType::NO_DEP,
            &LTRA::Model::resist,
            &LTRA::Model::resistGiven, U_OHMMM1, CAT_NONE,
            "Resistance per unit length");

    addPar ("L", 0.0, false, ParameterType::NO_DEP,
            &LTRA::Model::induct,
            &LTRA::Model::inductGiven, U_HMM1, CAT_NONE,
            "Inductance per unit length");

    addPar ("G", 0.0, false, ParameterType::NO_DEP,
            &LTRA::Model::conduct,
            &LTRA::Model::conductGiven, U_OHMM1MM1, CAT_NONE,
            "Conductance per unit length");

    addPar ("C", 0.0, false, ParameterType::NO_DEP,
            &LTRA::Model::capac,
            &LTRA::Model::capacGiven, U_FARADMM1, CAT_NONE,
            "Capacitance per unit length");

    addPar ("LEN", 0.0, false, ParameterType::NO_DEP,
            &LTRA::Model::length,
            &LTRA::Model::lengthGiven, U_METER, CAT_NONE,
            "length of line");

    addPar ("REL", 1.0, false, ParameterType::NO_DEP,
            &LTRA::Model::reltol,
            &LTRA::Model::reltolGiven, U_NONE, CAT_NONE,
            "Rel. rate of change of deriv. for bkpt");

    addPar ("ABS", 1.0, false, ParameterType::NO_DEP,
            &LTRA::Model::abstol,
            &LTRA::Model::abstolGiven, U_NONE, CAT_NONE,
            "Abs. rate of change of deriv. for bkpt");

    addPar ("STEPLIMIT", true, false, ParameterType::NO_DEP,
            &LTRA::Model::stepLimit,
            &LTRA::Model::stepLimitGiven, U_LOGIC, CAT_NONE,
            "limit timestep size based on the time constant of the line");

    addPar ("NOSTEPLIMIT", false, false, ParameterType::NO_DEP,
            &LTRA::Model::noStepLimit,
            &LTRA::Model::noStepLimitGiven, U_LOGIC, CAT_NONE,
            "don't limit timestep size based on the time constant of the line");

    addPar ("COMPLEXSTEPCONTROL", false, false, ParameterType::NO_DEP,
            &LTRA::Model::lteTimeStepControl,
            &LTRA::Model::lteTimeStepControlGiven, U_LOGIC, CAT_NONE,
            "do complex time step control using local truncation error estimation");

    addPar ("LININTERP", false, false, ParameterType::NO_DEP,
            &LTRA::Model::linInterp,
            &LTRA::Model::linInterpGiven, U_LOGIC, CAT_NONE,
            "use linear interpolation");

    addPar ("QUADINTERP", true, false, ParameterType::NO_DEP,
            &LTRA::Model::quadInterp,
            &LTRA::Model::quadInterpGiven, U_LOGIC, CAT_NONE,
            "use quadratic interpolation");

    addPar ("MIXEDINTERP", false, false, ParameterType::NO_DEP,
            &LTRA::Model::mixedInterp,
            &LTRA::Model::mixedInterpGiven, U_LOGIC, CAT_NONE,
            "use linear interpolation if quadratic results look unacceptable");

    addPar ("COMPACTREL", 1.0e-3, false, ParameterType::NO_DEP,
            &LTRA::Model::stLineReltol,
            &LTRA::Model::stLineReltolGiven, U_NONE, CAT_NONE,
            "special reltol for straight line checking");

    addPar ("COMPACTABS", 1.0e-12, false, ParameterType::NO_DEP,
            &LTRA::Model::stLineAbstol,
            &LTRA::Model::stLineAbstolGiven, U_NONE, CAT_NONE,
            "special abstol for straight line checking");

    addPar ("TRUNCNR", false, false, ParameterType::NO_DEP,
            &LTRA::Model::truncNR,
            &LTRA::Model::truncNRGiven,
            U_LOGIC, CAT_NONE, "use N-R iterations for step calculation in LTRAtrunc");

    addPar ("TRUNCDONTCUT", false, false, ParameterType::NO_DEP,
            &LTRA::Model::truncDontCut,
            &LTRA::Model::truncDontCutGiven,
            U_LOGIC, CAT_NONE, "don't limit timestep to keep impulse response calculation errors low");
}

namespace LTRA {

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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------

Instance::Instance(InstanceBlock& IB,
                                     Model & Miter,
                                     MatrixLoadData& mlData1,
                                     SolverState& ss1,
                                     ExternData& ed1,
                                     DeviceOptions& do1)

  : DeviceInstance(IB,mlData1,ss1,ed1,do1),
    model_(Miter),

    input1(0.0),
    input2(0.0),
    listSize(0),
    initVolt1(0.0),
    initVolt2(0.0),
    initCur1(0.0),
    initCur2(0.0),

    initVolt1Given(false),
    initVolt2Given(false),
    initCur1Given(false),
    initCur2Given(false),

    li_Pos1(-1),
    li_Neg1(-1),

    li_Pos2(-1),
    li_Neg2(-1),

    li_Ibr1(-1),
    li_Ibr2(-1),

    APos1EquPos1NodeOffset(-1),
    APos1EquIbr1NodeOffset(-1),

    ANeg1EquNeg1NodeOffset(-1),
    ANeg1EquIbr1NodeOffset(-1),

    APos2EquPos2NodeOffset(-1),
    APos2EquIbr2NodeOffset(-1),

    ANeg2EquNeg2NodeOffset(-1),
    ANeg2EquIbr2NodeOffset(-1),

    AIbr1EquPos1NodeOffset(-1),
    AIbr1EquNeg1NodeOffset(-1),
    AIbr1EquPos2NodeOffset(-1),
    AIbr1EquNeg2NodeOffset(-1),
    AIbr1EquIbr1NodeOffset(-1),
    AIbr1EquIbr2NodeOffset(-1),

    AIbr2EquPos1NodeOffset(-1),
    AIbr2EquNeg1NodeOffset(-1),
    AIbr2EquPos2NodeOffset(-1),
    AIbr2EquNeg2NodeOffset(-1),
    AIbr2EquIbr1NodeOffset(-1),
    AIbr2EquIbr2NodeOffset(-1),

    pos1Pos1Ptr(0),
    pos1Ibr1Ptr(0),

    neg1Neg1Ptr(0),
    neg1Ibr1Ptr(0),

    pos2Pos2Ptr(0),
    pos2Ibr2Ptr(0),

    neg2Neg2Ptr(0),
    neg2Ibr2Ptr(0),

    ibr1Pos1Ptr(0),
    ibr1Neg1Ptr(0),
    ibr1Pos2Ptr(0),
    ibr1Neg2Ptr(0),
    ibr1Ibr1Ptr(0),
    ibr1Ibr2Ptr(0),

    ibr2Pos1Ptr(0),
    ibr2Neg1Ptr(0),
    ibr2Pos2Ptr(0),
    ibr2Neg2Ptr(0),
    ibr2Ibr1Ptr(0),
    ibr2Ibr2Ptr(0),

    first_BP_call_done(false),
    newBreakPoint(false),
    newBreakPointTime(0.0)
{
  numIntVars   = 2;
  numExtVars   = 4;
  numStateVars = 0;

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 1;
  devConMap[2] = 2;
  devConMap[3] = 2;

  if( jacStamp.empty() )
  {
    jacStamp.resize(6);

    jacStamp[0].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 4;

    jacStamp[1].resize(2);
    jacStamp[1][0] = 1;
    jacStamp[1][1] = 4;

    jacStamp[2].resize(2);
    jacStamp[2][0] = 2;
    jacStamp[2][1] = 5;

    jacStamp[3].resize(2);
    jacStamp[3][0] = 3;
    jacStamp[3][1] = 5;

    jacStamp[4].resize(6);
    jacStamp[4][0] = 0;
    jacStamp[4][1] = 1;
    jacStamp[4][2] = 2;
    jacStamp[4][3] = 3;
    jacStamp[4][4] = 4;
    jacStamp[4][5] = 5;

    jacStamp[5].resize(6);
    jacStamp[5][0] = 0;
    jacStamp[5][1] = 1;
    jacStamp[5][2] = 2;
    jacStamp[5][3] = 3;
    jacStamp[5][4] = 4;
    jacStamp[5][5] = 5;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();

}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/16/10
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
// Creation Date : 06/16/10
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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int>& intLIDVecRef,
                                       const vector<int>& extLIDVecRef )
{
  ostringstream msg;

#if defined(Xyce_DEBUG_DEVICE)
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << dashedline << endl;
    cout << "[LTRA-DBG-DEV] In Instance::registerLIDs\n\n";
    cout << "name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

#if defined(Xyce_DEBUG_DEVICE)
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "[LTRA-DBG-DEV] number of internal variables: " << numInt << endl;
    cout << "[LTRA-DBG-DEV] number of external variables: " << numExt << endl;
  }
#endif

  if (numInt != numIntVars)
  {
    msg.clear();
    msg << "Instance::registerLIDs: ";
    msg << "numInt != numIntVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg.str());
  }

  if (numExt != numExtVars)
  {
    msg.clear();
    msg << "Instance::registerLIDs: ";
    msg << "numExt != numExtVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg.str());
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

  li_Ibr1 = intLIDVec[0];
  li_Ibr2 = intLIDVec[1];

#if defined(Xyce_DEBUG_DEVICE)
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "[LTRA-DBG-DEV]  VARIABLE Indicies " << endl;
    cout << "li_Pos1 = " << li_Pos1 << endl;
    cout << "li_Neg1 = " << li_Neg1 << endl;
    cout << "li_Pos2 = " << li_Pos2 << endl;
    cout << "li_Neg2 = " << li_Neg2 << endl;
    cout << "li_Ibr1 = " << li_Ibr1 << endl;
    cout << "li_Ibr1 = " << li_Ibr1 << endl;
    cout << "li_Ibr2 = " << li_Ibr2 << endl;
    cout << "li_Ibr2 = " << li_Ibr2 << endl;

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
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
map<int,string>& Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {

    // set up internal name map
    string tmpstr;
    tmpstr = getName()+"_branch1"; spiceInternalName (tmpstr);
    intNameMap[ li_Ibr1 ] = tmpstr;

    tmpstr = getName()+"_branch2"; spiceInternalName (tmpstr);
    intNameMap[ li_Ibr2 ] = tmpstr;

  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int>& staLIDVecRef)
{

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    ostringstream msg;
    msg << "Instance::registerStateLIDs: ";
    msg << "numSta != numStateVars";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg.str());
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
const vector< vector<int> >& Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> >& jacLIDVec )
{

  DeviceInstance::registerJacLIDs( jacLIDVec );

  APos1EquPos1NodeOffset = jacLIDVec[0][0];
  APos1EquIbr1NodeOffset = jacLIDVec[0][1];

  ANeg1EquNeg1NodeOffset = jacLIDVec[1][0];
  ANeg1EquIbr1NodeOffset = jacLIDVec[1][1];

  APos2EquPos2NodeOffset = jacLIDVec[2][0];
  APos2EquIbr2NodeOffset = jacLIDVec[2][1];

  ANeg2EquNeg2NodeOffset = jacLIDVec[3][0];
  ANeg2EquIbr2NodeOffset = jacLIDVec[3][1];

  AIbr1EquPos1NodeOffset = jacLIDVec[4][0];
  AIbr1EquNeg1NodeOffset = jacLIDVec[4][1];
  AIbr1EquPos2NodeOffset = jacLIDVec[4][2];
  AIbr1EquNeg2NodeOffset = jacLIDVec[4][3];
  AIbr1EquIbr1NodeOffset = jacLIDVec[4][4];
  AIbr1EquIbr2NodeOffset = jacLIDVec[4][5];

  AIbr2EquPos1NodeOffset = jacLIDVec[5][0];
  AIbr2EquNeg1NodeOffset = jacLIDVec[5][1];
  AIbr2EquPos2NodeOffset = jacLIDVec[5][2];
  AIbr2EquNeg2NodeOffset = jacLIDVec[5][3];
  AIbr2EquIbr1NodeOffset = jacLIDVec[5][4];
  AIbr2EquIbr2NodeOffset = jacLIDVec[5][5];

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
//
// Purpose       : Sets up pointers!?
//
// Special Notes :
//
// Scope         : public
// Creator       : Gary Hennigan, SNL
// Creation Date : 10/11/2012
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix& dFdx = *(extData.dFdxMatrixPtr);

  pos1Ibr1Ptr = &(dFdx[li_Pos1][APos1EquIbr1NodeOffset]);
  pos1Pos1Ptr = &(dFdx[li_Pos1][APos1EquPos1NodeOffset]);

  neg1Ibr1Ptr = &(dFdx[li_Neg1][ANeg1EquIbr1NodeOffset]);
  neg1Neg1Ptr = &(dFdx[li_Neg1][ANeg1EquNeg1NodeOffset]);

  pos2Ibr2Ptr = &(dFdx[li_Pos2][APos2EquIbr2NodeOffset]);
  pos2Pos2Ptr = &(dFdx[li_Pos2][APos2EquPos2NodeOffset]);

  neg2Ibr2Ptr = &(dFdx[li_Neg2][ANeg2EquIbr2NodeOffset]);
  neg2Neg2Ptr = &(dFdx[li_Neg2][ANeg2EquNeg2NodeOffset]);

  ibr1Ibr1Ptr = &(dFdx[li_Ibr1][AIbr1EquIbr1NodeOffset]);
  ibr1Ibr2Ptr = &(dFdx[li_Ibr1][AIbr1EquIbr2NodeOffset]);
  ibr1Pos1Ptr = &(dFdx[li_Ibr1][AIbr1EquPos1NodeOffset]);
  ibr1Neg1Ptr = &(dFdx[li_Ibr1][AIbr1EquNeg1NodeOffset]);
  ibr1Pos2Ptr = &(dFdx[li_Ibr1][AIbr1EquPos2NodeOffset]);
  ibr1Neg2Ptr = &(dFdx[li_Ibr1][AIbr1EquNeg2NodeOffset]);

  ibr2Ibr1Ptr = &(dFdx[li_Ibr2][AIbr2EquIbr1NodeOffset]);
  ibr2Ibr2Ptr = &(dFdx[li_Ibr2][AIbr2EquIbr2NodeOffset]);
  ibr2Pos1Ptr = &(dFdx[li_Ibr2][AIbr2EquPos1NodeOffset]);
  ibr2Neg1Ptr = &(dFdx[li_Ibr2][AIbr2EquNeg1NodeOffset]);
  ibr2Pos2Ptr = &(dFdx[li_Ibr2][AIbr2EquPos2NodeOffset]);
  ibr2Neg2Ptr = &(dFdx[li_Ibr2][AIbr2EquNeg2NodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 LTRA instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       : update primary state for one LTRA instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------

bool Instance::updatePrimaryState ()
{
  return updateIntermediateVars ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       : update secondary state for one LTRA instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------

bool Instance::updateSecondaryState()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one LTRA instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double* solVec = extData.nextSolVectorRawPtr;

  return true;
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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints ( vector<N_UTL_BreakPoint>& breakPointTimes )
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
#if defined(Xyce_DEBUG_DEVICE)
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        cout << "[LTRA-DBG-DEV]  In Instance::getBreakPoints "<<endl;
        cout << " First time step, I don't get to set breakpoints.  Time is ";
        cout << currentTime << endl;
      }
#endif
  }

  first_BP_call_done=true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       :
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{

  // This stores the voltage and current time history at the ports. Note
  // that both the dc-op and first time step have timeStepNumber 0 so we
  // have to distinguish between them. For the purposes of these
  // histories the dc-op is stored in index 0 and the first time step is
  // at index 1. This is consistent with NG-Spice.
  if (getSolverState().dcopFlag && listSize == 0)
  {
    listSize = 10;

    v1.resize(listSize);
    v2.resize(listSize);
    i1.resize(listSize);
    i2.resize(listSize);
  }
  else if (getSolverState().ltraTimeIndex >= listSize)
  {
    listSize += 10;

    v1.resize(listSize);
    v2.resize(listSize);
    i1.resize(listSize);
    i2.resize(listSize);
  }

  // because of DCOP at index 0
  v1[getSolverState().ltraTimeIndex] = vpos1 - vneg1;
  v2[getSolverState().ltraTimeIndex] = vpos2 - vneg2;

  i1[getSolverState().ltraTimeIndex] = currp1;
  i2[getSolverState().ltraTimeIndex] = currp2;

  // Allocate storage for time history entities
  if (getSolverState().initTranFlag && !getSolverState().dcopFlag)
  {
    model_.listSize = 10;

    model_.h1dashCoeffs.resize(model_.listSize);
    model_.h2Coeffs.resize(model_.listSize);
    model_.h3dashCoeffs.resize(model_.listSize);
  }
  else if (!getSolverState().dcopFlag && getSolverState().ltraTimeIndex >= model_.listSize)
  {
    model_.listSize += 10;

    model_.h1dashCoeffs.resize(model_.listSize);
    model_.h2Coeffs.resize(model_.listSize);
    model_.h3dashCoeffs.resize(model_.listSize);
  }

  bool compact = false;
  if (getDeviceOptions().tryToCompact && getSolverState().ltraTimeIndex >= 2) {

    // Figure out if the last 3 points line on a straight line for all
    // the termainal variables
    double t1 = getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex - 2];
    double t2 = getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex - 1];
    double t3 = getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex];

    compact = model_.straightLineCheck_(t1, v1[getSolverState().ltraTimeIndex-2],
                                            t2, v1[getSolverState().ltraTimeIndex-1],
                                            t3, v1[getSolverState().ltraTimeIndex],
                                            model_.stLineReltol,
                                            model_.stLineAbstol);
    if (compact) {
      compact = model_.straightLineCheck_(t1, v2[getSolverState().ltraTimeIndex-2],
                                              t2, v2[getSolverState().ltraTimeIndex-1],
                                              t3, v2[getSolverState().ltraTimeIndex],
                                              model_.stLineReltol,
                                              model_.stLineAbstol);
    }
    if (compact) {
      compact = model_.straightLineCheck_(t1, i1[getSolverState().ltraTimeIndex-2],
                                              t2, i1[getSolverState().ltraTimeIndex-1],
                                              t3, i1[getSolverState().ltraTimeIndex],
                                              model_.stLineReltol,
                                              model_.stLineAbstol);
    }
    if (compact) {
      compact = model_.straightLineCheck_(t1, i2[getSolverState().ltraTimeIndex-2],
                                              t2, i2[getSolverState().ltraTimeIndex-1],
                                              t3, i2[getSolverState().ltraTimeIndex],
                                              model_.stLineReltol,
                                              model_.stLineAbstol);
    }
  }

  if (getSolverState().ltraTimeIndex > 0)
  {
    double v1_ = (v1[getSolverState().ltraTimeIndex] + i1[getSolverState().ltraTimeIndex] *
                  model_.imped) * model_.attenuation;

    double v2_ = (v1[getSolverState().ltraTimeIndex - 1] +
                  i1[getSolverState().ltraTimeIndex - 1] *
                  model_.imped) * model_.attenuation;

    double v3_ = getSolverState().ltraTimeIndex < 2 ? v2_ : (v1[getSolverState().ltraTimeIndex-2] +
                                                      i1[getSolverState().ltraTimeIndex-2] *
                                                      model_.imped) * model_.attenuation;
    double v4_ = (v2[getSolverState().ltraTimeIndex] +
                  i2[getSolverState().ltraTimeIndex] *
                  model_.imped) * model_.attenuation;

    double v5_ = (v2[getSolverState().ltraTimeIndex-1] +
                  i2[getSolverState().ltraTimeIndex-1] *
                  model_.imped) * model_.attenuation;

    double v6_ = getSolverState().ltraTimeIndex < 2 ? v5_ : (v2[getSolverState().ltraTimeIndex-2] +
                                                      i2[getSolverState().ltraTimeIndex-2] *
                                                      model_.imped) * model_.attenuation;

    double d1_ = (v1_ - v2_) / (getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex]
                                - getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1]);

    double d2_ = getSolverState().ltraTimeIndex < 2 ? d1_ :
      (v2_ - v3_) / (getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1]
                     - getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-2]);

    double d3_ = (v4_ - v5_) / (getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex]
                                - getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1]);

    double d4_ = getSolverState().ltraTimeIndex < 2 ? d3_ :
      (v5_ - v6_) / (getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1]
                     - getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-2]);

#define CHECK(a,b,c) (Xycemax(Xycemax(a,b),c)-Xycemin(Xycemin(a,b),c) >= \
                      fabs(50.0*(getDeviceOptions().reltol/3.0*(a+b+c) + \
                                 getDeviceOptions().abstol)))

    bool tmp_test = (fabs(d1_ - d2_) > model_.reltol * Xycemax(fabs(d1_), fabs(d2_)) +
                     model_.abstol) && CHECK(v1_,v2_,v3_);

    if (tmp_test || ((fabs(d3_ - d4_)
                      >= model_.reltol * Xycemax(fabs(d3_), fabs(d4_)) +
                      model_.abstol) && CHECK(v4_,v5_,v6_)))
    {
      // Set breakpoint here
      newBreakPoint = true;
      newBreakPointTime = getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1] + model_.td;

#ifdef Xyce_DEBUG_DEVICE
      std::cout << "[LTRA-DBG-DEV]: At simulation time " << getSolverState().currTime
                << " adding a breakpoint at time " << newBreakPointTime
                << std::endl;
#endif
    }
  }

  if (getDeviceOptions().tryToCompact && compact && getSolverState().ltraTimeIndex >= 2)
  {
    v1[getSolverState().ltraTimeIndex-1] = v1[getSolverState().ltraTimeIndex];
    v2[getSolverState().ltraTimeIndex-1] = v2[getSolverState().ltraTimeIndex];
    i1[getSolverState().ltraTimeIndex-1] = i1[getSolverState().ltraTimeIndex];
    i2[getSolverState().ltraTimeIndex-1] = i2[getSolverState().ltraTimeIndex];

    getSolverState().ltraDoCompact = true;
  }

  calculateMaxTimeStep_();

#ifdef Xyce_DEBUG_DEVICE
  std::cout << "[LTRA-DBG-DEV]: At time: " << getSolverState().currTime
            << " max time step set to: " << model_.maxTimeStep
            << std::endl;
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::calculateMaxTimeStep_
// Purpose       : Calculates a maximum safe time step to avoid excessive
//                 errors
//
// Special Notes :
//
// Scope         : private
// Creator       : Gary Hennigan, SNL
// Creation Date : 12/04/2012
//-----------------------------------------------------------------------------
void Instance::calculateMaxTimeStep_()
{

  Model& model = model_;
  model.maxTimeStep = 1.0e99;

  if (getSolverState().ltraTimeIndex < 2)
  {
    model.maxTimeStep = Xycemin(model.td, model.maxSafeStep);
    return;
  }

  switch (model.specialCase)
  {
    case LTRA_MOD_LC:
    case LTRA_MOD_RLC:

      if (model.stepLimitType == LTRA_MOD_STEPLIMIT)
      {
        model.maxTimeStep = model.td;
      }
      else
      {
        size_t ti = getSolverState().ltraTimeIndex;

        // Approximate derivative to detect changing slope and adjust
        // time step accordingly
        double i1_ = (v2[ti] * model.admit + i2[ti]) * model.attenuation;
        double i2_ = (v2[ti-1] * model.admit + i2[ti-1]) * model.attenuation;
        double i3_ = (v2[ti-2] * model.admit + i2[ti-2]) * model.attenuation;

        double i4_ = (v1[ti] * model.admit + i1[ti]) * model.attenuation;
        double i5_ = (v1[ti-1] * model.admit + i1[ti-1]) * model.attenuation;
        double i6_ = (v1[ti-2] * model.admit + i1[ti-2]) * model.attenuation;

        double d1_ = (i1_ - i2_) /
          (getSolverState().ltraTimePoints[ti] - getSolverState().ltraTimePoints[ti-1]);

        double d2_ = (i2_ - i3_) /
          (getSolverState().ltraTimePoints[ti-1] - getSolverState().ltraTimePoints[ti-2]);

        double d3_ = (i4_ - i5_) /
          (getSolverState().ltraTimePoints[ti] - getSolverState().ltraTimePoints[ti-1]);

        double d4_ = (i5_ - i6_) /
          (getSolverState().ltraTimePoints[ti-1] - getSolverState().ltraTimePoints[ti-2]);

        if ((fabs(d1_-d2_) >= model.reltol * Xycemax(fabs(d1_), fabs(d2_)) + model.abstol) ||
            (fabs(d3_-d4_) >= model.reltol * Xycemax(fabs(d3_), fabs(d4_)) + model.abstol))
        {
          model.maxTimeStep = Xycemin(model.maxTimeStep, model.td);
        }

      }
      break;

    case LTRA_MOD_RC:
    case LTRA_MOD_RG:
      break;

    default:
        ostringstream msg;
        msg << "**********" << std::endl;
        msg << ": Error. Case not handled in calculateMaxTimeStep_()";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg.str());
        return;

  }

  /*
   * the above was for the parts of the equations that resemble the
   * lossless equations. Now we need to estimate the local truncation
   * error in each of the three convolution equations, and if possible
   * adjust the timestep so that all of them remain within some bound.
   * Unfortunately, the expression for the LTE in a convolution
   * operation is complicated and costly to evaluate; in addition, no
   * explicit inverse exists.
   *
   * So what we do here (for the moment) is check to see the current
   * error is acceptable. If so, the timestep is not changed. If not,
   * then an estimate is made for the new timestep using a few
   * iterations of the newton-raphson method.
   *
   * modification: we change the timestep to half its previous value
   */
  if ((model.specialCase == LTRA_MOD_RLC) && !(model.truncDontCut))
  {
    model.maxTimeStep = Xycemin(model.maxTimeStep, model.maxSafeStep);
  }

  // NOTE-GLH: None of the following code has been tested. As far as I
  // can tell there is no user option in Spice3 to turn this bit of code
  // on and as such there is no Xyce option to turn it on. Just to
  // reiterate it has NOT been tested or even run.
  if (model.lteTimeStepControl) {

    double current_lte;
    double tolerance;
    switch (model.specialCase) {

      case LTRA_MOD_RLC:
      case LTRA_MOD_RC:
        tolerance = 7.0 *
                    (getDeviceOptions().reltol * (fabs(input1) + fabs(input2)) + getDeviceOptions().abstol);

        current_lte = model.lteCalculate_(*this, getSolverState().currTime);

        if (current_lte >= tolerance) {
          if (model.truncNR) {

            double ti = getSolverState().ltraTimeIndex;
            double x = getSolverState().ltraTimePoints[ti];
            double y = current_lte;
            for (;;)
            {
              double deriv_delta = 0.01 * (x - getSolverState().ltraTimePoints[ti-1]);

#ifdef Xyce_DEBUG_DEVICE
              if (deriv_delta <= 0.0)
                std::cout << "LTRAtrunc: error: timestep is now less than zero" << std::endl;
#endif
              double deriv = model.lteCalculate_(*this, x + deriv_delta) - y;

              deriv /= deriv_delta;
              double change = (tolerance - y) / deriv;
              x += change;

              int maxiter=2;
              int iterations=0;
              if (maxiter == 0) {
                if (fabs(change) <= fabs(deriv_delta))
                  break;
              } else {
                iterations++;
                if (iterations >= maxiter)
                  break;
              }
              y = model.lteCalculate_(*this, x);
            }

            double tmp = x - getSolverState().ltraTimePoints[ti-1];
            model.maxTimeStep = Xycemin(model.maxTimeStep, tmp);
          }
          else
            model.maxTimeStep *= 0.5;
        }
        break;

      case LTRA_MOD_RG:
      case LTRA_MOD_LC:
        break;

      default:
        ostringstream msg;
        msg << "**********" << std::endl;
        msg << ": Error. Case not handled in calculateMaxTimeStep_() [2]";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg.str());
        return;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInternalState
// Purpose       : Generates an DeviceState object and populates
//                 it with the contents of the history vector for use by
//                 restarts.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
DeviceState * Instance::getInternalState()
{
  int i,j;
  // allocate obiect to return
  DeviceState * myState = new DeviceState;

  myState->ID=getName();

  // stuff owned by the instance:
  myState->dataSizeT.resize(2);
  myState->dataSizeT[0] = listSize;

  int origSize = myState->data.size();
  myState->data.resize(origSize + 4*listSize + 6);

  myState->data[origSize  ]=input1;
  myState->data[origSize+1]=input2;
  myState->data[origSize+2]= initVolt1 ;
  myState->data[origSize+3]= initVolt2 ;
  myState->data[origSize+4]= initCur1 ;
  myState->data[origSize+5]= initCur2 ;

#ifdef Xyce_DEBUG_RESTART
  cout.width(18); cout.precision(10); cout.setf(ios::scientific);
  std::cout << "LTRA::getInternalState:  input1="<<input1<<" input2="<<input2<<std::endl;
  std::cout << "LTRA::getInternalState:  initVolt11="<<initVolt1<<" initVolt2="<<initVolt2<<std::endl;
  std::cout << "LTRA::getInternalState:  initCur11="<<initCur1<<" initCur2="<<initCur2<<std::endl;
#endif

  for (i=0;i<listSize;++i)
  {
    j=(origSize+6)+i*4;
    myState->data[j  ]=v1[i];
    myState->data[j+1]=v2[i];
    myState->data[j+2]=i1[i];
    myState->data[j+3]=i2[i];
#ifdef Xyce_DEBUG_RESTART
    cout.width(18); cout.precision(10); cout.setf(ios::scientific);
    std::cout <<
      "LTRA::getInternalState:  v1["<<i<<"]="<<v1[i]
      <<" v2["<<i<<"]="<<v2[i]
      <<" i1["<<i<<"]="<<i1[i]
      <<" i2["<<i<<"]="<<i2[i]
      <<std::endl;
#endif
  }

  // stuff owned by the model:
  //if (!(model_.restartStoredFlag))
  //{
    myState->dataSizeT[1] = model_.listSize;

    origSize = myState->data.size();
    myState->data.resize(origSize+model_.listSize*3);
    for (i=0;i<model_.listSize;++i)
    {
      j=origSize+i*3;
      myState->data[j]=model_.h1dashCoeffs[i];
      myState->data[j+1]=model_.h2Coeffs[i];
      myState->data[j+2]=model_.h3dashCoeffs[i];
#ifdef Xyce_DEBUG_RESTART
    cout.width(18); cout.precision(10); cout.setf(ios::scientific);
    std::cout <<
      "LTRA::getInternalState:  h1dashCoeffs["<<i<<"] =" << model_.h1dashCoeffs[i]
          <<" h2Coeffs["<<i<<"] =" << model_.h2Coeffs[i]
          <<" h3dashCoeffs["<<i<<"] =" << model_.h3dashCoeffs[i]<<std::endl;
#endif
    }

    //model_.restartStoredFlag=true;
  //}

  return myState;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setInternalState
// Purpose       : Reload history data from restart
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Instance::setInternalState(const DeviceState &state)
{
  int i,j;
  if ( state.ID != getName())
  {
    string msg;
    msg = "Instance::setInternalState:  ID ("+state.ID+")";
    msg += "from restart does not match my name ("+getName()+")!\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // stuff owned by the instance:
  listSize=state.dataSizeT[0];
  v1.clear(); v2.clear(); i1.clear(); i2.clear();
  v1.resize(listSize); v2.resize(listSize); i1.resize(listSize); i2.resize(listSize);

  input1=state.data[0];
  input2=state.data[1];
  initVolt1=state.data[2];
  initVolt2=state.data[3];
  initCur1=state.data[4];
  initCur2=state.data[5];

#ifdef Xyce_DEBUG_RESTART
  cout.width(18); cout.precision(10); cout.setf(ios::scientific);
  std::cout << "LTRA::setInternalState:  input1="<<input1<<" input2="<<input2<<std::endl;
  std::cout << "LTRA::setInternalState:  initVolt11="<<initVolt1<<" initVolt2="<<initVolt2<<std::endl;
  std::cout << "LTRA::setInternalState:  initCur11="<<initCur1<<" initCur2="<<initCur2<<std::endl;
#endif

  for ( i=0; i<listSize; ++i)
  {
    j=6+i*4;
    v1[i]= state.data[j  ];
    v2[i]= state.data[j+1];
    i1[i]= state.data[j+2];
    i2[i]= state.data[j+3];
#ifdef Xyce_DEBUG_RESTART
    cout.width(18); cout.precision(10); cout.setf(ios::scientific);
    std::cout <<
      "LTRA::setInternalState:  v1["<<i<<"]="<<v1[i]
      <<" v2["<<i<<"]="<<v2[i]
      <<" i1["<<i<<"]="<<i1[i]
      <<" i2["<<i<<"]="<<i2[i]
      <<std::endl;
#endif
  }

  // stuff owned by the model:
  model_.listSize=state.dataSizeT[1];

  model_.h1dashCoeffs.clear();
  model_.h2Coeffs.clear();
  model_.h3dashCoeffs.clear();

  model_.h1dashCoeffs.resize(model_.listSize);
  model_.h2Coeffs.resize(model_.listSize);
  model_.h3dashCoeffs.resize(model_.listSize);

  for ( i=0; i<model_.listSize; ++i)
  {
    j=(listSize*4+6)+i*3;
    model_.h1dashCoeffs[i]= state.data[j];
    model_.h2Coeffs[i]= state.data[j+1];
    model_.h3dashCoeffs[i]= state.data[j+2];
#ifdef Xyce_DEBUG_RESTART
    cout.width(18); cout.precision(10); cout.setf(ios::scientific);
    std::cout <<
      "LTRA::setInternalState:  h1dashCoeffs["<<i<<"] =" << model_.h1dashCoeffs[i]
          <<" h2Coeffs["<<i<<"] =" << model_.h2Coeffs[i]
          <<" h3dashCoeffs["<<i<<"] =" << model_.h3dashCoeffs[i]<<std::endl;
#endif
  }

  return true;
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{

  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock& MB,
                                      SolverState& ss1,
                                      DeviceOptions& do1)
  : DeviceModel(MB,ss1,do1),

    h1dashFirstVal(0.0),
    h2FirstVal(0.0),
    h3dashFirstVal(0.0),
    h1dashFirstCoeff(0.0),
    h2FirstCoeff(0.0),
    h3dashFirstCoeff(0.0),
    listSize(0),
    resist(0.0),
    induct(0.0),
    conduct(0.0),
    capac(0.0),
    length(0.0),
    reltol(0.0),
    abstol(0.0),

    noStepLimit(false),
    stepLimit(true),
    stepLimitType(LTRA_MOD_STEPLIMIT),

    linInterp(false),
    quadInterp(true),
    mixedInterp(false),

    stLineReltol(0.0),
    stLineAbstol(0.0),

    truncNR(false),
    truncDontCut(false),

    resistGiven(false),
    inductGiven(false),
    conductGiven(false),
    capacGiven(false),
    lengthGiven(false),
    reltolGiven(false),
    abstolGiven(false),
    noStepLimitGiven(false),
    stepLimitGiven(false),
    linInterpGiven(false),
    quadInterpGiven(false),
    mixedInterpGiven(false),
    stLineReltolGiven(false),
    stLineAbstolGiven(false),
    truncNRGiven(false),
    truncDontCutGiven(false),

    td(0.0),
    imped(0.0),
    admit(0.0),
    alpha(0.0),
    beta(0.0),
    attenuation(0.0),
    cByR(0.0),
    rclsqr(0.0),
    intH1dash(0.0),
    intH2(0.0),
    intH3dash(0.0),

    coshlrootGR(0.0),
    rRsLrGRorG(0.0),
    rGsLrGRorR(0.0),
    auxIndex(0),
    chopReltol(0.0),
    chopAbstol(0.0),

    maxSafeStep(1.0e99),
    maxTimeStep(1.0e99),
    lteTimeStepControl(false),
    howToInterp(LTRA_MOD_QUADINTERP),
    printFlag(false),
    specialCase(0),
    tdover(false),
    restartStoredFlag(false)

{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  processParams ();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
bool Instance::setIC()
{

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
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
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
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

    os << endl;
  }
  os << endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::modelCalculations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Model::modelCalculations_(int& isaved,
                                         double& qf1, double& qf2, double& qf3,
                                         double& lf2, double& lf3)
{
  double t1(0.0),t2(0.0),t3(0.0);
  double v1d(0.0), v2d(0.0), i1d(0.0), i2d(0.0);
  double dummy1(0.0), dummy2(0.0);

  // Initialize index
  isaved = 0;

  if( getSolverState().dcopFlag)
  {
    switch (specialCase)
    {
      case LTRA_MOD_RG:
        dummy1 = length*sqrt(resist*conduct);
        dummy2 = exp(-dummy1);
        dummy1 = exp(dummy1); // warning: may overflow!
        coshlrootGR = 0.5*(dummy1 + dummy2);

        if (conduct <= 1.0e-10)
        { // hack!
          rRsLrGRorG = length*resist;
        }
        else
        {
          rRsLrGRorG = 0.5*(dummy1 - dummy2)*sqrt(resist/conduct);
        }

        if (resist <= 1.0e-10)
        { // hack!
          rGsLrGRorR = length*conduct;
        }
        else
        {
          rGsLrGRorR = 0.5*(dummy1 - dummy2)*sqrt(conduct/resist);
        }
        break;

      case LTRA_MOD_RC:
      case LTRA_MOD_LC:
      case LTRA_MOD_RLC:
        // simple resistor-like behaviour nothing to set up
        break;

      default:
        ostringstream msg;
        msg << "**********" << std::endl;
        msg << ": Error. Case not handled in modelCalculations_()";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg.str());
        return false;
    }
  }
  else
  {
    switch (specialCase)
    {
      case LTRA_MOD_RLC:
      case LTRA_MOD_LC:

        if (getSolverState().currTime > td)
          tdover = true;
        else
          tdover = false;

      default:
        break;
    }

    switch (specialCase)
    {
      case LTRA_MOD_RLC:

        // set up lists of values of the functions at the
        // necessary timepoints.

        // set up coefficient lists LTRAh1dashCoeffs,
        // LTRAh2Coeffs, LTRAh3dashCoeffs for current
        // timepoint

        // NOTE: h1, h2 and h3 here actually refer to h1tilde,
        // h2tilde, h3tilde in the paper

        // Note: many function evaluations are saved by doing
        // the following all together in one procedure

        (void) rlcCoeffsSetup_( h1dashFirstCoeff, h2FirstCoeff, h3dashFirstCoeff,
                                h1dashCoeffs, h2Coeffs, h3dashCoeffs,
                                listSize,
                                td, alpha, beta,
                                getSolverState().currTime,
                                getSolverState().ltraTimePoints,
                                getSolverState().ltraTimeIndex,
                                chopReltol,
                                &(auxIndex));

      case LTRA_MOD_LC:
        // setting up the coefficients for interpolation
        if (tdover)
        { // serious hack -fix!
          int i = 0;
          for (i = getSolverState().ltraTimeIndex; i >= 0; i--)
          {
            if (getSolverState().ltraTimePoints[i] < getSolverState().currTime - td)
              break;

          }
#ifdef Xyce_DEBUG_DEVICE
          if (i == getSolverState().ltraTimeIndex)
          {
            std::cout << "[LTRA-DBG-DEV] LTRAload: Warning: timestep larger than delay of line" << std::endl;
            std::cout << "\tTime now: " << getSolverState().currTime << std::endl << std::endl;
          }
#endif

          if (i == getSolverState().ltraTimeIndex)
            i--;

          if (i == -1)
          {
#ifdef Xyce_DEBUG_DEVICE
            std::cout << "[LTRA-DBG-DEV] LTRAload: mistake: cannot find delayed timepoint" << std::endl;
#endif
            ostringstream msg;
            msg << "************" << std::endl;
            msg << ": Error. Delayed time point not found.";
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg.str());
            return false;
          }

          isaved = i;

          t2 = getSolverState().ltraTimePoints[i];
          t3 = getSolverState().ltraTimePoints[i+1];

          // quadratic interpolation
          if ((i != 0) && ((howToInterp == LTRA_MOD_QUADINTERP)|| (howToInterp == LTRA_MOD_MIXEDINTERP)))
          {
            t1 = getSolverState().ltraTimePoints[i-1];
            quadInterp_(getSolverState().currTime-td, t1, t2, t3, qf1, qf2, qf3);
          }

          // linear interpolation
          if ( (i == 0) || (howToInterp == LTRA_MOD_MIXEDINTERP) || (howToInterp ==  LTRA_MOD_LININTERP))
          {
            linInterp_(getSolverState().currTime-td, t2, t3, lf2, lf3);
          }
        }

        // interpolation coefficients set-up
        break;

      case LTRA_MOD_RC:

        /*
         * set up lists of values of the coefficients at the
         * necessary timepoints.
         */

        /* set up coefficient lists LTRAh1dashCoeffs,
           LTRAh2Coeffs, LTRAh3dashCoeffs for current
           timepoint*/

        /* Note: many function evaluations are saved by doing the
         * following all together in one procedure
         */

        (void)
          rcCoeffsSetup_(h1dashFirstCoeff, h2FirstCoeff, h3dashFirstCoeff,
                         h1dashCoeffs, h2Coeffs, h3dashCoeffs,
                         listSize,
                         cByR, rclsqr,
                         getSolverState().currTime,
                         getSolverState().ltraTimePoints,
                         getSolverState().ltraTimeIndex,
                         chopReltol);

        break;

      case LTRA_MOD_RG:
        break;

      default:
        return false;
        //return(E_BADPARM);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::quadInterp_
// Purpose       :
//
// quadratic interpolation function
// t = timepoint where value wanted
// t1, t2, t3 are three timepoints where the value is known
// c1, c2, c3 are set to the proper coefficients by the function
// the interpolated value is c1*v1 + c2*v2 + c3*v3; this should be
// done in the calling program; (v1,v2,v3 are the known values at
// t1,t2,t3)
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
int Model::quadInterp_ (double t, double t1, double t2, double t3, double& c1, double& c2, double& c3)
{
  double f1, f2, f3;

  if (t == t1)
  {
    c1 = 1.0;
    c2 = 0.0;
    c3 = 0.0;
    return(0);
  }
  if (t == t2)
  {
    c1 = 0.0;
    c2 = 1.0;
    c3 = 0.0;
    return(0);
  }
  if (t == t3)
  {
    c1 = 0.0;
    c2 = 0.0;
    c3 = 1.0;
    return(0);
  }
  if( (t2-t1)==0  || (t3-t2) == 0 || (t1 - t3) ==0) return(1);

  f1 = (t - t2) * (t - t3) ;
  f2 = (t - t1) * (t - t3) ;
  f3 = (t - t1) * (t - t2) ;
  if((t2-t1)==0)
  { /* should never happen, but don't want
     * to divide by zero, EVER... */
    f1=0;
    f2=0;
  }
  else
  {
    f1 /= (t1-t2);
    f2 /= (t2-t1);
  }
  if((t3-t2)==0)
  { /* should never happen, but don't want
     * to divide by zero, EVER... */
    f2=0;
    f3=0;
  }
  else
  {
    f2 /= (t2-t3);
    f3 /= (t2-t3);
  }
  if((t3-t1)==0)
  { /* should never happen, but don't want
     * to divide by zero, EVER... */
    f1=0;
    f2=0;
  }
  else
  {
    f1 /= (t1-t3);
    f3 /= (t1-t3);
  }
  c1 = f1;
  c2 = f2;
  c3 = f3;
  return(0);
}

//-----------------------------------------------------------------------------
// Function      : Model::linInterp_
// Purpose       : linear interpolation
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
int Model::linInterp_  (double t, double t1, double t2, double& c1, double& c2)
{
  double temp;

  if (t1 == t2) return(1);

  if (t==t1)
  {
    c1 = 1.0;
    c2 = 0.0;
    return(0);
  }

  if (t==t2)
  {
    c1 = 0.0;
    c2 = 1.0;
    return(0);
  }

  temp = (t-t1)/(t2-t1);
  c2 = temp;
  c1 = 1-temp;

  return(0);
}

//-----------------------------------------------------------------------------
// Function      : Model::intlinfunc_
// Purpose       :
//
// intlinfunc returns \int_lolimit^hilimit h(\tau) d \tau, where
// h(\tau) is assumed to be linear, with values lovalue and hivalue
// \tau = t1 and t2 respectively
// this is used only locally
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::intlinfunc_ (double lolimit, double hilimit,
                                     double lovalue, double hivalue,
                                     double t1, double t2)
{
  double width, m;

  width = t2 - t1;
  if (width == 0.0) return(0.0);
  m = (hivalue - lovalue)/width;

  return( (hilimit-lolimit)*lovalue + 0.5*m*((hilimit-t1)*(hilimit-t1)
    - (lolimit - t1)*(lolimit - t1)));
}

//-----------------------------------------------------------------------------
// Function      : Model::twiceintlinfunc_
// Purpose       :
//
// twiceintlinfunc returns \int_lolimit^hilimit \int_otherlolimit^\tau
// h(\tau') d \tau' d \tau , where
// h(\tau') is assumed to be linear, with values lovalue and hivalue
// \tau = t1 and t2 respectively
// this is used only locally
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::twiceintlinfunc_(double lolimit, double hilimit,
                                         double otherlolimit, double lovalue,
                                         double hivalue, double t1, double t2)
{
  double width, m, dummy;
  double temp1, temp2, temp3;

  width = t2 - t1;
  if (width == 0.0) return(0.0);
  m = (hivalue - lovalue)/width;

  temp1 = hilimit - t1;
  temp2 = lolimit - t1;
  temp3 = otherlolimit - t1;
  dummy = lovalue*((hilimit - otherlolimit)*(hilimit - otherlolimit) -
    (lolimit - otherlolimit)*(lolimit - otherlolimit));
  dummy += m*((temp1*temp1*temp1 - temp2*temp2*temp2)/3.0 -
    temp3*temp3*(hilimit - lolimit));
  return(dummy*0.5);
}


//-----------------------------------------------------------------------------
// Function      : Model::thriceintlinfunc_
// Purpose       :
//
// thriceintlinfunc returns \int_lolimit^hilimit \int_secondlolimit^\tau
// \int_thirdlolimit^\tau' h(\tau'') d \tau'' d \tau' d \tau , where
// h(\tau'') is assumed to be linear, with values lovalue and hivalue
// \tau = t1 and t2 respectively
// this is used only locally
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::thriceintlinfunc_(double lolimit, double hilimit,
                                          double secondlolimit, double thirdlolimit,
                                          double lovalue, double  hivalue, double t1, double t2)
{
  double width, m, dummy;
  double temp1, temp2, temp3, temp4;
  double temp5, temp6, temp7, temp8, temp9, temp10;


  width = t2 - t1;
  if (width == 0.0) return(0.0);
  m = (hivalue - lovalue)/width;

  temp1 = hilimit - t1;
  temp2 = lolimit - t1;
  temp3 = secondlolimit - t1;
  temp4 = thirdlolimit - t1;
  temp5 = hilimit - thirdlolimit;
  temp6 = lolimit - thirdlolimit;
  temp7 = secondlolimit - thirdlolimit;
  temp8 = hilimit - lolimit;
  temp9 = hilimit - secondlolimit;
  temp10 = lolimit - secondlolimit;
  dummy = lovalue*((temp5*temp5*temp5 - temp6*temp6*temp6)/3 -
    temp7*temp5*temp8);
  dummy += m*(((temp1*temp1*temp1*temp1 - temp2*temp2*temp2*temp2)*0.25 -
    temp3*temp3*temp3*temp8)/3 - temp4*temp4*0.5*(temp9*temp9 -
    temp10*temp10));
  return(dummy*0.5);
}

// from numerical recipies in C:
//-----------------------------------------------------------------------------
// Function      : Model::bessI0_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::bessI0_(double x)
{
  double ax,ans;
  double y;

  if ((ax=fabs(x)) < 3.75)
  {
    y=x/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
                                         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  }
  else
  {
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))*
      (0.39894228+y*(0.1328592e-1+y*(0.225319e-2+
                                     y*(-0.157565e-2+y*(0.916281e-2+
                                       y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
                                         +y*0.392377e-2))))))));
  }
  return(ans);
}

//-----------------------------------------------------------------------------
// Function      : Model::bessI1_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::bessI1_(double x)
{
  double ax,ans;
  double y;

  if ((ax=fabs(x)) < 3.75)
  {
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
                                               +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  }
  else
  {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
                                         -y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
                                       +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return(x < 0.0 ? -ans : ans);
}

//-----------------------------------------------------------------------------
// Function      : Model::bessI1xOverX_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::bessI1xOverX_(double x)
{
  double ax,ans;
  double y;

  if ((ax=fabs(x)) < 3.75) {
    y=x/3.75;
    y*=y;
    ans=0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
                                           +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))));
  }
  else
  {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
                                         -y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
                                       +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/(ax*sqrt(ax)));
  }
  return(ans);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH1dashFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH1dashFunc_(double time, double T, double alpha, double beta)
{
  double besselarg, exparg, returnval;
  /* T is not used in this function */

  /* result = alpha * e^{- beta*time} * {I_1(alpha*time) -
  * I_0(alpha*time)}
  */

  if (alpha == 0.0) return(0.0);

  exparg = - beta * time;
  besselarg = alpha*time;

  returnval = (bessI1_(besselarg)-bessI0_(besselarg))* alpha * exp(exparg);
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH2Func_
// Purpose       : first impulse response function
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH2Func_(double time, double T, double alpha, double beta)
{
  double besselarg, exparg, returnval;

  /*
  * result = 0, time < T
  *      = (alpha*T*e^{-beta*time})/sqrt(t^2 - T^2) *
  *        I_1(alpha*sqrt(t^2 - T^2)), time >= T
  */

  if (alpha == 0.0) return(0.0);
  if (time < T) return(0.0);

  if (time != T) {
    besselarg = alpha*sqrt(time*time - T*T);
  } else {
    besselarg = 0.0;
  }
  exparg = -beta*time;

  returnval = alpha*alpha*T*exp(exparg)*bessI1xOverX_(besselarg);
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH3dashFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH3dashFunc_(double time, double T, double alpha, double beta)
{
  double exparg,besselarg,returnval;

  /*
  * result = 0, time < T
  *      = alpha*e^{-beta*time}*(t/sqrt(t^2-T^2)*
  *      I_1(alpha*sqrt(t^2-T^2)) - I_0(alpha*sqrt(t^2-T^2)))
  */

  if (alpha == 0.0) return(0.0);
  if (time < T) return(0.0);

  exparg = - beta*time;
  if (time != T) {
    besselarg = alpha*sqrt(time*time - T*T);
  } else {
    besselarg = 0.0;
  }

  returnval = alpha*time*bessI1xOverX_(besselarg) - bessI0_(besselarg);
  returnval *= alpha*exp(exparg);
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH1dashTwiceIntFunc_
//
// Purpose       : Twice repeated integral of h1dash for the
//                 special case of G = 0
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH1dashTwiceIntFunc_(double time, double beta)
{
  double arg, returnval;

  /* result = time * e^{- beta*time} * {I_0(beta*time) +
  * I_1(beta*time)} - time
  */

  if (beta == 0.0) return(time);
  arg = beta*time;
  if (arg == 0.0) return(0.0);

  returnval = (bessI1_(arg)+bessI0_(arg))* time * exp(-arg) - time;
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcH3dashIntFunc_
//
// Purpose       : twice repeated integral of h1dash for the
//                 special case of G = 0
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rlcH3dashIntFunc_(double time, double T, double beta)
{
  double exparg, besselarg;
  double returnval;

  if (time <= T) return(0.0);
  if (beta == 0.0) return(0.0);
  exparg = -beta*time;
  besselarg = beta*sqrt(time*time-T*T);
  returnval = exp(exparg)* bessI0_(besselarg) - exp(-beta*T);
  return(returnval);
}

//-----------------------------------------------------------------------------
// Function      : Model::rcH1dashTwiceIntFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rcH1dashTwiceIntFunc_(double time, double cbyr)
{
  return(sqrt(4*cbyr*time/M_PI));
}

//-----------------------------------------------------------------------------
// Function      : Model::rcH2TwiceIntFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rcH2TwiceIntFunc_(double time, double rclsqr)
{
  double temp(0.0);
  if (time != 0.0)
  {
    temp = rclsqr/(4*time);

    // FIXME: The Intel compiler's location of erf/erfc requires an
    // include of mathimf.h, which is incompatible with math.h.  This is
    // completely screwed, because math.h is used EVERYWHERE in Xyce.
    // Must either refactor all of Xyce not to depend everywhere on
    // math.h, or must refactor this device to provide an alternative
    // erf/erfc function on Windows. Right now I'm incorporating a
    // third-party implementation of erfc() to avoid this problem

    double erfc_res = Faddeeva::erfc(sqrt(temp));

    return((time + rclsqr*0.5)*erfc_res - sqrt(time*rclsqr/M_PI)*exp(- temp));
  }
  else
  {
    return(0.0);
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::rcH3dashTwiceIntFunc_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::rcH3dashTwiceIntFunc_(double time, double cbyr, double rclsqr)
{
  double temp;
  if (time != 0.0)
  {
    temp =  rclsqr/(4*time);

    // see note in rcH2TwiceIntFunc_ about intel compilers
    double erfc_res = Faddeeva::erfc(sqrt(temp));

    temp = 2*sqrt(time/M_PI)*exp(-temp) - sqrt(rclsqr)*erfc_res;
    return(sqrt(cbyr)*temp);
  }
  else
  {
    return(0.0);
  }
}

// coefficient setups:
//-----------------------------------------------------------------------------
// Function      : Model::rcCoeffsSetup_
// Purpose       :
//
// Sets up the all coefficient lists for the special case where L=G=0
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
void Model::rcCoeffsSetup_(
    double& h1dashfirstcoeff,
    double& h2firstcoeff,
    double& h3dashfirstcoeff,
    vector<double>& h1dashcoeffs,
    vector<double>& h2coeffs,
    vector<double>& h3dashcoeffs,
    size_t listsize, double cbyr, double rclsqr, double curtime,
    vector<double>& timelist, int timeindex, double reltol)
{
  double delta1, delta2;
  double h1dummy1, h1dummy2;
  double h2dummy1, h2dummy2;
  double h3dummy1, h3dummy2;
  double lolimit1,lolimit2, hilimit1, hilimit2;
  double h1lovalue1,h1lovalue2,h1hivalue1,h1hivalue2;
  double h2lovalue1,h2lovalue2,h2hivalue1,h2hivalue2;
  double h3lovalue1,h3lovalue2,h3hivalue1,h3hivalue2;
  double temp, temp2, temp3, temp4, temp5;
  double h1relval, h2relval, h3relval;
  int doh1=1, doh2=1, doh3=1;
  int i,auxindex;

  /* coefflists should already have been allocated to the necessary size */

#ifdef Xyce_DEBUG_DEVICE
  if (listsize < timeindex) {
    std::cout << "[LTRA-DBG-DEV]: LTRAcoeffSetup: not enough space in coefflist" << std::endl;
  }
#endif

  auxindex = timeindex;

  /* the first coefficients */

  delta1 = curtime - timelist[auxindex];
  lolimit1 = 0.0;
  hilimit1 = delta1;

  h1lovalue1 = 0.0;
  h1hivalue1 = /* LTRArcH1dashTwiceIntFunc(hilimit1,cbyr); */
        sqrt(4*cbyr*hilimit1/M_PI);
  h1dummy1 = h1hivalue1/delta1;
  h1dashfirstcoeff = h1dummy1;
  h1relval = fabs(h1dummy1*reltol);

  temp = rclsqr/(4*hilimit1);

  // see note in :rcH2TwiceIntFunc_ re intel compiler
  temp2 = (temp >= 100.0 ? 0.0 : Faddeeva::erfc(sqrt(temp)));
  temp3 = exp(-temp);
  temp4 = sqrt(rclsqr);
  temp5 = sqrt(cbyr);

  h2lovalue1 = 0.0;
  h2hivalue1 = /* LTRArcH2TwiceIntFunc(hilimit1,rclsqr); */
  (hilimit1 != 0.0?  (hilimit1 + rclsqr*0.5)*temp2 - sqrt(hilimit1*rclsqr/M_PI)*temp3 : 0.0);


  h2dummy1 = h2hivalue1/delta1;
  h2firstcoeff = h2dummy1;
  h2relval = fabs(h2dummy1*reltol);

  h3lovalue1 = 0.0;
  h3hivalue1 = /* LTRArcH3dashTwiceIntFunc(hilimit1,cbyr,rclsqr); */
  (hilimit1 != 0.0? temp = 2*sqrt(hilimit1/M_PI)*temp3 - temp4*temp2, (temp5*temp): 0.0);


  h3dummy1 = h3hivalue1/delta1;
  h3dashfirstcoeff = h3dummy1;
  h3relval = fabs(h3dummy1*reltol);

  /* the coefficients for the rest of the timepoints */

  for (i=auxindex; i>0; i--)
  {
    delta2 = delta1; /* previous delta1 */
    lolimit2 = lolimit1; /* previous lolimit1 */
    hilimit2 = hilimit1; /*previous hilimit1 */

    delta1 = timelist[i] - timelist[i - 1];
    lolimit1 = hilimit2;
    hilimit1 = curtime - timelist[i - 1];

    if (doh1)
    {
      h1lovalue2 = h1lovalue1; /* previous lovalue1 */
      h1hivalue2 = h1hivalue1; /*previous hivalue1 */
      h1dummy2 = h1dummy1; /* previous dummy1 */

      h1lovalue1 = h1hivalue2;
      h1hivalue1 = /* LTRArcH1dashTwiceIntFunc(hilimit1,cbyr); */
            sqrt(4*cbyr*hilimit1/M_PI);
      h1dummy1 = (h1hivalue1 - h1lovalue1)/delta1;
      h1dashcoeffs[i] = h1dummy1 - h1dummy2;
      if (fabs(h1dashcoeffs[i]) < h1relval) doh1=0;
    }
    else
    {
      h1dashcoeffs[i] = 0.0;
    }

    if (doh2 || doh3) {
    temp = rclsqr/(4*hilimit1);
    // see note in :rcH2TwiceIntFunc_ re intel compiler
    temp2 = (temp >= 100.0 ? 0.0 : Faddeeva::erfc(sqrt(temp)));
    temp3 = exp(-temp);
    }

    if (doh2)
    {
      h2lovalue2 = h2lovalue1; /* previous lovalue1 */
      h2hivalue2 = h2hivalue1; /*previous hivalue1 */
      h2dummy2 = h2dummy1; /* previous dummy1 */

      h2lovalue1 = h2hivalue2;
      h2hivalue1 = /* LTRArcH2TwiceIntFunc(hilimit1,rclsqr); */
          (hilimit1 != 0.0?  (hilimit1 + rclsqr*0.5)*temp2 - sqrt(hilimit1*rclsqr/M_PI)*temp3 : 0.0);
      h2dummy1 = (h2hivalue1 - h2lovalue1)/delta1;
      h2coeffs[i] = h2dummy1 - h2dummy2;
      if (fabs(h2coeffs[i]) < h2relval) doh2=0;
    }
    else
    {
      h2coeffs[i] = 0.0;
    }

    if (doh3)
    {
      h3lovalue2 = h3lovalue1; /* previous lovalue1 */
      h3hivalue2 = h3hivalue1; /*previous hivalue1 */
      h3dummy2 = h3dummy1; /* previous dummy1 */

      h3lovalue1 = h3hivalue2;
      h3hivalue1 = /*LTRArcH3dashTwiceIntFunc(hilimit1,cbyr,rclsqr);*/
          (hilimit1 != 0.0? temp = 2*sqrt(hilimit1/M_PI)*temp3 - temp4*temp2, (temp5*temp): 0.0);
      h3dummy1 = (h3hivalue1 - h3lovalue1)/delta1;
      h3dashcoeffs[i] = h3dummy1 - h3dummy2;
      if (fabs(h3dashcoeffs[i]) < h3relval) doh3=0;
    }
    else
    {
      h3dashcoeffs[i] = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::rlcCoeffsSetup_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
void Model::rlcCoeffsSetup_(
    double& h1dashfirstcoeff,
    double& h2firstcoeff,
    double& h3dashfirstcoeff,
    vector<double>& h1dashcoeffs,
    vector<double>& h2coeffs,
    vector<double>& h3dashcoeffs,
    size_t listsize,
    double T, double alpha, double beta, double curtime,
    vector<double>& timelist, int timeindex, double reltol, int* auxindexptr)
{
  unsigned exact;
  double lolimit1,lolimit2,hilimit1,hilimit2;
  double delta1, delta2;

  double h1dummy1, h1dummy2;
  double h1lovalue1,h1lovalue2,h1hivalue1,h1hivalue2;

  double h2dummy1, h2dummy2;
  double h2lovalue1,h2lovalue2,h2hivalue1,h2hivalue2;

  double h3dummy1, h3dummy2;
  double h3lovalue1,h3lovalue2,h3hivalue1,h3hivalue2;

  double exparg, besselarg, expterm, bessi1overxterm, bessi0term;
  double expbetaTterm, alphasqTterm;
  double h1relval, h2relval, h3relval;
  int doh1=1, doh2=1, doh3=1;

  int i,auxindex;

  /* coefflists should already have been allocated to the necessary size */

#ifdef Xyce_DEBUG_DEVICE
  if (listsize < timeindex) {
    std::cout << "[LTRA-DBG-DEV]: LTRArlcCoeffsSetup_: not enough space in coefflist" << std::endl;
  }
#endif


  /*
  * we assume a piecewise linear function, and we calculate the
  * coefficients using this assumption in the integration of the
  * function
  */

  if (T == 0.0) {
    auxindex = timeindex;
  } else {

    if (curtime - T <= 0.0) {
      auxindex = 0;
    } else {
      exact = 0;
      for (i = timeindex; i>= 0; i--) {
        if (curtime - timelist[i] ==  T) {
          exact =1;
          break;
        }
        if (curtime - timelist[i] > T) break;
      }

#ifdef Xyce_DEBUG_DEVICE
      if ((i < 0) || ((i==0) && (exact==1)))
        std::cout << "[LTRA-DBG-DEV]: LTRAcoeffSetup: i <= 0: some mistake!" << std::endl;
#endif

      if (exact == 1) {
        auxindex = i-1;
      } else {
        auxindex = i;
      }
    }
  }
  /* the first coefficient */

  if (auxindex != 0)
  {
    lolimit1 = T;
    hilimit1 = curtime - timelist[auxindex];
    delta1 = hilimit1 - lolimit1;

    h2lovalue1 = rlcH2Func_(T,T,alpha,beta);
    besselarg = (hilimit1 > T) ? alpha*sqrt(hilimit1*hilimit1-T*T):0.0;
    exparg = -beta*hilimit1;
    expterm = exp(exparg);
    bessi1overxterm = bessI1xOverX_(besselarg);
    alphasqTterm = alpha*alpha*T;
    h2hivalue1 = /* LTRArlcH2Func(hilimit1,T,alpha,beta); */
    ((alpha == 0.0) || (hilimit1 < T)) ? 0.0: alphasqTterm*expterm*bessi1overxterm;

    h2dummy1 = twiceintlinfunc_(lolimit1,hilimit1,lolimit1,h2lovalue1,
      h2hivalue1,lolimit1,hilimit1)/delta1;
    h2firstcoeff = h2dummy1;
    h2relval = fabs(reltol*h2dummy1);

    h3lovalue1 = 0.0; /* E3dash should be consistent with this */
    bessi0term = bessI0_(besselarg);
    expbetaTterm = exp(-beta*T);
    h3hivalue1 = /*LTRArlcH3dashIntFunc(hilimit1,T,beta);*/
    ((hilimit1 <= T) || (beta == 0.0)) ? 0.0: expterm* bessi0term-expbetaTterm;
    h3dummy1 = intlinfunc_(lolimit1,hilimit1,h3lovalue1,
      h3hivalue1,lolimit1,hilimit1)/delta1;
    h3dashfirstcoeff = h3dummy1;
    h3relval = fabs(h3dummy1*reltol);
  }
  else
  {
    h2firstcoeff = h3dashfirstcoeff = 0.0;
  }

  lolimit1 = 0.0;
  hilimit1 = curtime - timelist[timeindex];
  delta1 = hilimit1 - lolimit1;
  exparg = -beta*hilimit1;
  expterm = exp(exparg);

  h1lovalue1 = 0.0;
  h1hivalue1 = /*LTRArlcH1dashTwiceIntFunc(hilimit1,beta);*/
    (beta == 0.0) ? hilimit1 : ((hilimit1 == 0.0) ? 0.0 :
                               (bessI1_(-exparg)+bessI0_(-exparg))* hilimit1 * expterm - hilimit1);
  h1dummy1 = h1hivalue1/delta1;
  h1dashfirstcoeff = h1dummy1;
  h1relval = fabs(h1dummy1*reltol);


  /* the coefficients for the rest of the timepoints */

  for (i=timeindex; i>0; i--)
  {
    if (doh1 || doh2 || doh3)
    {
      lolimit2 = lolimit1; /* previous lolimit1 */
      hilimit2 = hilimit1; /*previous hilimit1 */
      delta2 = delta1; /* previous delta1 */

      lolimit1 = hilimit2;
      hilimit1 = curtime - timelist[i - 1];
      delta1 = timelist[i] - timelist[i - 1];

      exparg = -beta*hilimit1;
      expterm = exp(exparg);
    }

    if (doh1)
    {
      h1lovalue2 = h1lovalue1; /* previous lovalue1 */
      h1hivalue2 = h1hivalue1; /*previous hivalue1 */
      h1dummy2 = h1dummy1; /* previous dummy1 */

      h1lovalue1 = h1hivalue2;
      h1hivalue1 = /* LTRArlcH1dashTwiceIntFunc(hilimit1,beta);*/
        (beta == 0.0) ? hilimit1 : ((hilimit1 == 0.0) ? 0.0 :
                                    (bessI1_(-exparg)+bessI0_(-exparg))* hilimit1 * expterm - hilimit1);
      h1dummy1 = (h1hivalue1 - h1lovalue1)/delta1;

      h1dashcoeffs[i] = h1dummy1 - h1dummy2;
      if (fabs(h1dashcoeffs[i]) <= h1relval) doh1 = 0;
    }
    else
    {
      h1dashcoeffs[i] = 0.0;
    }

    if (i <= auxindex)
    {
      /*
      if (i == auxindex) {
      lolimit2 = T;
      delta2 = hilimit2 - lolimit2;
      }
      */

      if (doh2 || doh3)
      {
        besselarg = (hilimit1 > T) ? alpha*sqrt(hilimit1*hilimit1-T*T):0.0;
      }

      if (doh2)
      {
        h2lovalue2 = h2lovalue1; /* previous lovalue1 */
        h2hivalue2 = h2hivalue1; /*previous hivalue1 */
        h2dummy2 = h2dummy1; /* previous dummy1 */

        h2lovalue1 = h2hivalue2;
        bessi1overxterm = bessI1xOverX_(besselarg);
        h2hivalue1 = /*rlcH2Func(hilimit1,T,alpha,beta);*/
        ((alpha == 0.0) || (hilimit1 < T)) ? 0.0: alphasqTterm*expterm*bessi1overxterm;
        h2dummy1 = twiceintlinfunc_(lolimit1,hilimit1,lolimit1,
        h2lovalue1,h2hivalue1,lolimit1,hilimit1)/delta1;

        h2coeffs[i] = h2dummy1 - h2dummy2 + intlinfunc_(lolimit2,hilimit2,
        h2lovalue2,h2hivalue2,lolimit2,hilimit2);
        if (fabs(h2coeffs[i]) <= h2relval) doh2 = 0;
      }
      else
      {
        h2coeffs[i] = 0.0;
      }

      if (doh3)
      {
        h3lovalue2 = h3lovalue1; /* previous lovalue1 */
        h3hivalue2 = h3hivalue1; /*previous hivalue1 */
        h3dummy2 = h3dummy1; /* previous dummy1 */

        h3lovalue1 = h3hivalue2;
        bessi0term = bessI0_(besselarg);
        h3hivalue1 = /*LTRArlcH3dashIntFunc(hilimit1,T,beta);*/
        ((hilimit1 <= T) || (beta == 0.0)) ? 0.0: expterm* bessi0term-expbetaTterm;
        h3dummy1 = intlinfunc_(lolimit1,hilimit1,h3lovalue1,h3hivalue1,lolimit1,hilimit1)/delta1;

        h3dashcoeffs[i] = h3dummy1 - h3dummy2;
        if (fabs(h3dashcoeffs[i]) <= h3relval) doh3 = 0;
      }
      else
      {
        h3dashcoeffs[i] = 0.0;
      }
    }
  }
  *auxindexptr = auxindex;
}

//-----------------------------------------------------------------------------
// Function      : Model::straightLineCheck_
// Purpose       :
//
// takes the co-ordinates of three points,
// finds the area of the triangle enclosed by these points and
// compares this area with the area of the quadrilateral formed by
// the line between the first point and the third point, the
// perpendiculars from the first and third points to the x-axis, and
// the x-axis. If within reltol, then it returns 1, else 0. The
// purpose of this function is to determine if three points lie
// acceptably close to a straight line. This area criterion is used
// because it is related to integrals and convolution
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
bool Model::straightLineCheck_(double x1, double y1,
                                         double x2, double y2,
                                         double x3, double y3,
                                         double reltol, double abstol)
{
  /*
  double asqr, bsqr, csqr, c, c1sqr;
  double htsqr;
  */
  double TRarea, QUADarea1,QUADarea2,QUADarea3, area;

  /*
  asqr = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
  bsqr = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2);
  csqr = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1);
  c = sqrt(csqr);
  c1sqr = (asqr - bsqr + csqr)/(2*c);
  c1sqr *= c1sqr;
  htsqr = asqr - c1sqr;
  TRarea = c*sqrt(htsqr)*0.5;
  */
  /* this should work if y1,y2,y3 all have the same sign and x1,x2,x3
  are in increasing order*/

  QUADarea1 = (fabs(y2)+fabs(y1))*0.5*fabs(x2-x1);
  QUADarea2 = (fabs(y3)+fabs(y2))*0.5*fabs(x3-x2);
  QUADarea3 = (fabs(y3)+fabs(y1))*0.5*fabs(x3-x1);
  TRarea = fabs( QUADarea3 - QUADarea1 - QUADarea2);
  area = QUADarea1 + QUADarea2;
  if (area*reltol + abstol > TRarea)
    return(true);
  else
    return(false);
}

// i is the index of the latest value,
// a,b,c values correspond to values at t_{i-2}, t{i-1} and t_i
//
// ERK: Note: check curtime.
#define SECONDDERIV(i,a,b,c) \
  (oof = (i==getSolverState().ltraTimeIndex?getSolverState().currTime:                  \
          (getSolverState().ltraTimePoints[i])),                                \
   (( c - b )/(oof-(getSolverState().ltraTimePoints[i-1])) -                    \
    ( b - a )/((getSolverState().ltraTimePoints[i-1])-                          \
               (getSolverState().ltraTimePoints[i-2])))/(oof -                  \
                                                 (getSolverState().ltraTimePoints[i-2])))


//-----------------------------------------------------------------------------
// Function      : Model::SECONDDERIV_
// Purpose       :
// Special Notes : see macro, above, modified from spice3
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::SECONDDERIV_(int i, double a, double b, double c)
{
  double oof=0.0;
  return SECONDDERIV(i,a,b,c);
}

//-----------------------------------------------------------------------------
// Function      : Model::lteCalculate_
// Purpose       :
//
// returns sum of the absolute values of the total
// local truncation error of the 2 equations for the LTRAline
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/05/10
//-----------------------------------------------------------------------------
double Model::lteCalculate_ (
  Instance& instance,
  double curtime
  )
{
  double h1dashTfirstCoeff;
  double h2TfirstCoeff;
  double h3dashTfirstCoeff;
  double dashdash;
  double hilimit1, lolimit1, hivalue1, lovalue1, f1i, g1i;
  double eq1LTE=0.0, eq2LTE=0.0;
  int auxindex, tdover, i, exact;

  switch(specialCase)
  {
    case LTRA_MOD_LC:
    case LTRA_MOD_RG:
      return(0.0);
      break;

    case LTRA_MOD_RLC:

      if (curtime > td)
      {
        tdover = 1;

        exact = 0;

        for (i=(getSolverState().ltraTimeIndex-1); i >= 0; i--)
        {
          if (curtime - getSolverState().ltraTimePoints[i] ==  td)
          {
            exact = 1;
            break;
          }
          if (curtime - getSolverState().ltraTimePoints[i] > td)
          {
            break;
          }
        }

#ifdef Xyce_DEBUG_DEVICE
        if ((i < 0) || ((i==0) && (exact==1)))
          std::cout << "[LTRA-DBG-DEV]: lteCalculate_: i <= 0: some mistake!" << std::endl;
#endif

        if (exact == 1)
        {
          auxindex = i-1;
        }
        else
        {
          auxindex = i;
        }
      }
      else
      {
        tdover = 0;
      }

      hilimit1 = curtime - getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1];
      lolimit1 = 0.0;
      hivalue1 = rlcH1dashTwiceIntFunc_(hilimit1,beta);
      lovalue1 = 0.0;

      f1i = hivalue1;
      g1i = intlinfunc_(lolimit1,hilimit1,lovalue1,hivalue1,
                        lolimit1,hilimit1);
      h1dashTfirstCoeff = 0.5 * f1i *
        (curtime - getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1]) - g1i;

      if (tdover)
      {
        hilimit1 = curtime - getSolverState().ltraTimePoints[auxindex];
        lolimit1 = getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1] - getSolverState().ltraTimePoints[auxindex];
        lolimit1 = Xycemax(td,lolimit1);

        // are the following really doing the operations in the write-up?
        hivalue1 = rlcH2Func_(hilimit1,td,alpha,beta);
        lovalue1 = rlcH2Func_(lolimit1,td,alpha,beta);
        f1i = twiceintlinfunc_(lolimit1,hilimit1,lolimit1,lovalue1,hivalue1,lolimit1,
                               hilimit1);
        g1i = thriceintlinfunc_(lolimit1,hilimit1,lolimit1,lolimit1,lovalue1,
                                hivalue1,lolimit1,hilimit1);

        h2TfirstCoeff = 0.5*f1i*(curtime-td-getSolverState().ltraTimePoints[auxindex]) - g1i;

        hivalue1 = rlcH3dashIntFunc_(hilimit1,td,beta);
        lovalue1 = rlcH3dashIntFunc_(lolimit1,td,beta);
        f1i = intlinfunc_(lolimit1,hilimit1,lovalue1,hivalue1,lolimit1,
                          hilimit1);
        g1i = twiceintlinfunc_(lolimit1,hilimit1,lolimit1,lovalue1,
                               hivalue1,lolimit1,hilimit1);
        h3dashTfirstCoeff = 0.5*f1i*(curtime-td-getSolverState().ltraTimePoints[auxindex]) - g1i;
      }


      /* LTEs for convolution with v1 */
      /* get divided differences for v1 (2nd derivative estimates) */

      /*
        * no need to subtract operating point values because
        * taking differences anyway
        */


      dashdash = SECONDDERIV_(getSolverState().ltraTimeIndex,
                              instance.v1[getSolverState().ltraTimeIndex-2],
                              instance.v1[getSolverState().ltraTimeIndex-1],
                              instance.v1[getSolverState().ltraTimeIndex]);
      eq1LTE += admit*fabs(dashdash * h1dashTfirstCoeff);

      // not bothering to interpolate since everything is approximate
      // anyway
      if (tdover)
      {
        dashdash = SECONDDERIV_(auxindex+1,
                                instance.v1[auxindex - 1],
                                instance.v1[auxindex],
                                instance.v1[auxindex + 1]) ;

        eq2LTE += admit*fabs(dashdash * h3dashTfirstCoeff);
      }
      /* end LTEs for convolution with v1 */

      /* LTEs for convolution with v2 */
      /* get divided differences for v2 (2nd derivative estimates) */

      dashdash = SECONDDERIV_(getSolverState().ltraTimeIndex,
                              instance.v2[getSolverState().ltraTimeIndex-2],
                              instance.v2[getSolverState().ltraTimeIndex-1],
                              instance.v2[getSolverState().ltraTimeIndex]);

      eq2LTE += admit*fabs(dashdash * h1dashTfirstCoeff);

      if (tdover)
      {
        dashdash = SECONDDERIV_(auxindex+1,
                                instance.v2[auxindex - 1],
                                instance.v2[auxindex],
                                instance.v2[auxindex + 1]);

        eq1LTE += admit*fabs(dashdash * h3dashTfirstCoeff);
      }

      /* end LTEs for convolution with v2 */

      /* LTE for convolution with i1 */
      /* get divided differences for i1 (2nd derivative estimates) */

      if (tdover)
      {
        dashdash = SECONDDERIV_(auxindex+1,
                                instance.i1[auxindex - 1],
                                instance.i1[auxindex],
                                instance.i1[auxindex + 1]) ;

        eq2LTE += fabs(dashdash * h2TfirstCoeff);
      }
      /* end LTE for convolution with i1 */

      /* LTE for convolution with i2 */
      /* get divided differences for i2 (2nd derivative estimates) */

      if (tdover)
      {
        dashdash = SECONDDERIV_(auxindex+1,
                                instance.i2[auxindex - 1],
                                instance.i2[auxindex],
                                instance.i2[auxindex + 1]) ;

        eq1LTE += fabs(dashdash * h2TfirstCoeff);
      }

      /* end LTE for convolution with i1 */

      break;

    case LTRA_MOD_RC:

      hilimit1 = curtime - getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1];
      lolimit1 = 0.0;

      hivalue1 = rcH1dashTwiceIntFunc_(hilimit1,cByR);
      lovalue1 = 0.0;

      f1i = hivalue1;
      g1i = intlinfunc_(lolimit1,hilimit1,lovalue1,hivalue1,lolimit1,hilimit1);

      h1dashTfirstCoeff = 0.5*f1i*(curtime-getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1]) - g1i;

      hivalue1 = rcH2TwiceIntFunc_(hilimit1,rclsqr);
      lovalue1 = 0.0;

      f1i = hivalue1;
      g1i = intlinfunc_(lolimit1,hilimit1,lovalue1,hivalue1,lolimit1,hilimit1);
      h1dashTfirstCoeff = 0.5*f1i*(curtime-getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1]) - g1i;

      hivalue1 = rcH2TwiceIntFunc_(hilimit1,rclsqr);
      lovalue1 = 0.0;

      f1i = hivalue1;
      g1i = intlinfunc_(lolimit1,hilimit1,lovalue1,
                        hivalue1,lolimit1,hilimit1);
      h1dashTfirstCoeff = 0.5*f1i*(curtime-getSolverState().ltraTimePoints[getSolverState().ltraTimeIndex-1]) - g1i;

      /* LTEs for convolution with v1 */
      /* get divided differences for v1 (2nd derivative estimates) */

      /*
        * no need to subtract operating point values because
        * taking differences anyway
        */

      dashdash = SECONDDERIV_( getSolverState().ltraTimeIndex,
                               instance.v1[getSolverState().ltraTimeIndex-2],
                               instance.v1[getSolverState().ltraTimeIndex-1],
                               instance.v1[getSolverState().ltraTimeIndex] );

      eq1LTE += fabs(dashdash * h1dashTfirstCoeff);
      eq2LTE += fabs(dashdash * h3dashTfirstCoeff);

      /* end LTEs for convolution with v1 */

      /* LTEs for convolution with v2 */
      /* get divided differences for v2 (2nd derivative estimates) */

      dashdash = SECONDDERIV_( getSolverState().ltraTimeIndex,
                               instance.v2[getSolverState().ltraTimeIndex-2],
                               instance.v2[getSolverState().ltraTimeIndex-1],
                               instance.v2[getSolverState().ltraTimeIndex] );

      eq2LTE += fabs(dashdash * h1dashTfirstCoeff);
      eq1LTE += fabs(dashdash * h3dashTfirstCoeff);

      /* end LTEs for convolution with v2 */

      /* LTE for convolution with i1 */
      /* get divided differences for i1 (2nd derivative estimates) */

      dashdash = SECONDDERIV_( getSolverState().ltraTimeIndex,
                               instance.i1[getSolverState().ltraTimeIndex-2],
                               instance.i1[getSolverState().ltraTimeIndex-1],
                               instance.i1[getSolverState().ltraTimeIndex] );

      eq2LTE += fabs(dashdash * h2TfirstCoeff);

      /* end LTE for convolution with i1 */

      /* LTE for convolution with i2 */
      /* get divided differences for i2 (2nd derivative estimates) */

      dashdash = SECONDDERIV_( getSolverState().ltraTimeIndex,
                               instance.i2[getSolverState().ltraTimeIndex-2],
                               instance.i2[getSolverState().ltraTimeIndex-1],
                               instance.i2[getSolverState().ltraTimeIndex] );

      eq1LTE += fabs(dashdash * h2TfirstCoeff);

      /* end LTE for convolution with i1 */

      break;

    default:
      return(1/*error*/);
  }

#ifdef Xyce_DEBUG_DEVICE
  std::cout << "[LTRA-DBG-DEV] " << instance.getName() << ": LTE/input for Eq1 at time "
            << curtime << " is: " << eq1LTE/instance.input1 << std::endl;

  std::cout << "[LTRA-DBG-DEV] " << instance.getName() << ": LTE/input for Eq2 at time "
            << curtime << " is: " << eq2LTE/instance.input1 << std::endl;
#endif

  return(fabs(eq1LTE) + fabs(eq2LTE));
}

//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize ()
{
  return model_.maxTimeStep;
}

// LTRA Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Master::updateState (double* solVec, double* staVec, double* stoVec)
{

  // Compute some quantities derived from input parameters and do some
  // error checking.
  if (!vars_initialized)
  {
    initialize_vars_();

    // Set flag to avoid multiple invocations
    vars_initialized = true;
  }

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance& di = *(*it);

    // Update current state
    di.vpos1 = solVec[di.li_Pos1];
    di.vneg1 = solVec[di.li_Neg1];

    di.vpos2 = solVec[di.li_Pos2];
    di.vneg2 = solVec[di.li_Neg2];

    di.currp1 = solVec[di.li_Ibr1];
    di.currp2 = solVec[di.li_Ibr2];

    // Initial state, generally the end result of the DC-OP calculation
    if (getSolverState().dcopFlag)
    {
      di.initVolt1 = di.vpos1 - di.vneg1;
      di.initVolt2 = di.vpos2 - di.vneg2;

      di.initCur1 = di.currp1;
      di.initCur2 = di.currp2;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : Error check and set some variables derived from input
//                 parameters.
// Special Notes : This should only be invoked once after the parameters from
//                 the input file have been initialized. This should be
//                 invoked from a function only once after the paramater
//                 initialization and before the solve starts.
// Scope         : private
// Creator       : Gary Hennigan
// Creation Date : 10/25/2012
//-----------------------------------------------------------------------------
void Master::initialize_vars_(void)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {

    Model& m = (*it)->getModel();

    // If tolerances for the line interpolation aren't given set them to
    // the same as the tolerances for the device options
    if (m.stLineReltol == 0.0)
      m.stLineReltol = getDeviceOptions().reltol;
    if (m.stLineAbstol == 0.0)
      m.stLineAbstol = getDeviceOptions().abstol;

    // Initialize which case this is based on the nonzero user-specified
    // parameters.
    if ((m.resist == 0) && (m.conduct == 0) && (m.capac != 0) && (m.induct != 0))
      m.specialCase = LTRA_MOD_LC;

    else if ((m.resist != 0) && (m.conduct == 0) && (m.capac != 0) && (m.induct != 0))
      m.specialCase = LTRA_MOD_RLC;

    else if ((m.resist != 0) && (m.conduct == 0) && (m.capac != 0) && (m.induct == 0))
      m.specialCase = LTRA_MOD_RC;

    else if ((m.resist != 0) && (m.conduct == 0) && (m.capac == 0) && (m.induct != 0))
      m.specialCase = LTRA_MOD_RL;

    else if ((m.resist != 0) && (m.conduct != 0) && (m.capac == 0) && (m.induct == 0))
      m.specialCase = LTRA_MOD_RG;

    else if ((m.conduct != 0) && ((m.capac != 0) || (m.induct != 0)))
    {
      ostringstream msg;
      msg << "**********" << std::endl;
      msg << "LTRA: Error. RL line not supported. "
          << "Modes supported: RC, RG, LC, RLC";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, msg.str());
      m.specialCase = LTRA_MOD_LTRA;
    }
    else if ((m.resist != 0) && (m.conduct == 0) && (m.capac == 0) && (m.induct != 0))
    {
      ostringstream msg;
      msg << "**********" << std::endl;
      msg << "LTRA: Error. Nonzero G (except RG) transmission line not supported. "
          << "Modes supported: RC, RG, LC, RLC";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, msg.str());
      m.specialCase = LTRA_MOD_LTRA;
    }

    if ((m.resist == 0.0 ? 0 : 1) + (m.conduct == 0.0 ? 0 : 1) +
        (m.induct == 0.0 ? 0 : 1) + (m.capac == 0.0 ? 0 : 1) <= 1)
    {
      ostringstream msg;
      msg << "************" << std::endl;
      msg << "LTRA: Error. Invalid specification. Specify at least "
          << "two of R, L, G, or C with nonzero values. "
          << "Modes supported: RC, RG, LC, RLC";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, msg.str());
      m.specialCase = LTRA_MOD_LTRA;
    }

    // Override the interpolation, either default or user specified, if
    // the TRYTOCOMPACT option is specified.
    if (getDeviceOptions().tryToCompact) {
      m.howToInterp = LTRA_MOD_LININTERP;
    }

    if (m.stepLimit && m.noStepLimit)
    {
      ostringstream msg;
      msg << "************" << std::endl
          << "LTRA: Warning. Conflicting options STEPLIMIT and NOSTEPLIMIT "
          << "given. Using STEPLIMIT";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, msg.str());
      m.stepLimitType = LTRA_MOD_STEPLIMIT;
    }
    else if (m.stepLimit || !m.noStepLimit)
    {
      m.stepLimitType = LTRA_MOD_STEPLIMIT;
    }
    else if (m.noStepLimit || !m.stepLimit)
    {
      m.stepLimitType = LTRA_MOD_NOSTEPLIMIT;
    }
    else
    { // default
      m.stepLimitType = LTRA_MOD_STEPLIMIT;
    }

    // Calculate some derived parameters
    switch (m.specialCase)
    {

      case LTRA_MOD_LC:
        m.imped = sqrt(m.induct / m.capac);
        m.admit = 1.0 / m.imped;
        m.td = sqrt(m.induct*m.capac) * m.length;
        m.attenuation = 1.0;
        break;

      case LTRA_MOD_RLC:
        m.imped = sqrt(m.induct / m.capac);
        m.admit = 1.0 / m.imped;
        m.td = sqrt(m.induct * m.capac) * m.length;
        m.alpha = 0.5 * (m.resist / m.induct);
        m.beta = m.alpha;
        m.attenuation = exp(-m.beta * m.td);

        if (m.alpha > 0.0)
        {
          m.intH1dash = -1.0;
          m.intH2 = 1.0 - m.attenuation;
          m.intH3dash = -m.attenuation;
        }
        else
        {
          m.intH1dash = m.intH2 = m.intH3dash = 0.0;
        }

        // Sanity check
        if (m.alpha < 0.0) {
          ostringstream msg;
          msg << "**********" << std::endl;
          msg << ": Error. Resistance and inductance ";
          msg << "must be greater than zero";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg.str());
          return;
        }

        // Calculate the time step size limit in order to keep
        // impulse-response errors low
        if (!m.truncDontCut) {
          double xbig, xsmall, xmid, y1big, y1small, y1mid;
          double y2big, y2small, y2mid;
          int done = 0, maxiter = 50, iters = 0;

          xbig = 10.0 * m.td;
          xsmall = m.td;
          xmid = 0.5 * (xbig + xsmall);
          y1small = m.rlcH2Func_(xsmall, m.td, m.alpha, m.beta);
          y2small = m.rlcH3dashFunc_(xsmall, m.td, m.beta, m.beta);
          iters = 0;
          for (;;) {

            iters++;
            y1big = m.rlcH2Func_(xbig, m.td, m.alpha, m.beta);
            y1mid = m.rlcH2Func_(xmid, m.td, m.alpha, m.beta);
            y2big = m.rlcH3dashFunc_(xbig, m.td, m.beta, m.beta);
            y2mid = m.rlcH3dashFunc_(xmid, m.td, m.beta, m.beta);
            done = m.straightLineCheck_(xbig, y1big, xmid, y1mid, xsmall,
                                        y1small, m.stLineReltol,
                                        m.stLineAbstol) +
	      m.straightLineCheck_(xbig, y1big, xmid, y1mid, xsmall,
                                   y1small, m.stLineReltol,
                                   m.stLineAbstol);
            if ((done == 2) || (iters > maxiter))
              break; // out of "for (;;)"
            xbig = xmid;
            xmid = 0.5 * (xbig + xsmall);
          }
          m.maxSafeStep = xbig - m.td;
        }
        break;

      case LTRA_MOD_RC:
        m.cByR = m.capac / m.resist;
        m.rclsqr = m.resist * m.capac * m.length * m.length;
        m.intH1dash = 0.0;
        m.intH2 = 1.0;
        m.intH3dash = 0.0;
        break;

      case LTRA_MOD_RG:
        break;

      default:
      {
        ostringstream msg;
        msg << "************" << std::endl;
        msg << ": Error. Unhandled LTRA special case encountered.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg.str());
        return;
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  double max(0.0),min(0.0);
  double v1d(0.0), v2d(0.0), i1d(0.0), i2d(0.0);
  double dummy1(0.0), dummy2(0.0);
  ostringstream msg;

  int sizeInstances = getInstanceVector().size();

  for (int i=0; i<sizeInstances; ++i)
  {
    Instance& di = *getInstanceVector()[i];

    if( getSolverState().dcopFlag || (di.getModel().specialCase == LTRA_MOD_RG))
    {
      switch (di.getModel().specialCase)
      {
        case LTRA_MOD_RG:
          dummy1 = di.getModel().length * std::sqrt(di.getModel().resist *
                                                    di.getModel().conduct);
          dummy2 = exp(-dummy1);
          dummy1 = exp(dummy1);  // May overflow
          di.getModel().coshlrootGR = 0.5 * (dummy1 + dummy2);

          if (di.getModel().conduct <= 1.0e-10)
          {	// Spice3 hack!
            di.getModel().rRsLrGRorG = di.getModel().length *
              di.getModel().resist;
          }
          else
          {
            di.getModel().rRsLrGRorG =
              0.5 * (dummy1 - dummy2) * sqrt(di.getModel().resist /
                                             di.getModel().conduct);
          }

          if (di.getModel().resist <= 1.0e-10)
          {	// Spice3 hack!
            di.getModel().rGsLrGRorR = di.getModel().length *
	      di.getModel().conduct;
          }
          else
          {
            di.getModel().rGsLrGRorR =
              0.5 * (dummy1 - dummy2) * sqrt(di.getModel().conduct /
                                             di.getModel().resist);
          }

          fVec[di.li_Ibr1] += (di.vpos1 -
                               di.vneg1 -
                               di.getModel().coshlrootGR * di.vpos2 +
                               di.getModel().coshlrootGR * di.vneg2 +
                               (1.0 + getDeviceOptions().gmin) * di.getModel().rRsLrGRorG * di.currp2);

          fVec[di.li_Ibr2] += (di.getModel().coshlrootGR * di.currp2 -
                               (1.0 + getDeviceOptions().gmin) * di.getModel().rGsLrGRorR * di.vpos2 +
                               (1.0 + getDeviceOptions().gmin) * di.getModel().rGsLrGRorR * di.vneg2 +
                               di.currp1);

          break;

        // load a simple resistor (DC case, C=open, L=short). In the
        // lossless case (R=0.0) the port voltages are equal.
        case LTRA_MOD_LC:
        case LTRA_MOD_RC:
        case LTRA_MOD_RLC:

          // i_1 + i_2 = 0
          fVec[di.li_Ibr1] += (di.currp1 + di.currp2);

          // v_(n1+) - v_(n2+) - R_ltra * i_1 = 0
          fVec[di.li_Ibr2] += (di.vpos1 - di.vpos2 -
                               di.currp1 *
                               di.getModel().resist *
                               di.getModel().length);

          break;

        default:

          msg.clear();
          msg << "**********" << std::endl;
          msg << ": Error. Unknown LTRA configuration, ";
          msg << di.getModel().specialCase;
          msg << ". Must be one of RG, LC, RC, or RLC.";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg.str());

          return false;
      }

      // These are common for all DC cases. They are just the residuals
      // that enforce that the current out of the positive terminal of
      // the TL is equal to the current in to the negative terminal at
      // the same end of the TL.
      fVec[di.li_Pos1] += di.currp1;  // i_(n1+) = i_1
      fVec[di.li_Neg1] += -di.currp1; // i_(n1-) = -i_1

      fVec[di.li_Pos2] += di.currp2;  // i_(n2+) = i_2
      fVec[di.li_Neg2] += -di.currp2; // i_(n2-) = -i_2

    }
    else
    {
      // all cases other than DC or the RG case

      int isaved = 0;
      double qf1, qf2, qf3;
      double lf2, lf3;

      qf1 = qf2 = qf3 = 0.0;
      lf2 = lf3 = 0.0;

      di.getModel().modelCalculations_(isaved, qf1, qf2, qf3, lf2, lf3);

      di.input1 = di.input2 = 0.0;

      switch (di.getModel().specialCase)
      {
        case LTRA_MOD_LC:
        case LTRA_MOD_RLC:

          if (di.getModel().tdover)
          {
            /* have to interpolate values */
            if ((isaved != 0) &&
                ((di.getModel().howToInterp == LTRA_MOD_QUADINTERP) || (di.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)))
            {
              v1d = di.v1[isaved-1] * qf1
                + di.v1[isaved]   * qf2
                + di.v1[isaved+1] * qf3;

              max = Xycemax(di.v1[isaved-1], di.v1[isaved]);
              max = Xycemax(max,di.v1[isaved+1]);
              min = Xycemin(di.v1[isaved-1], di.v1[isaved]);
              min = Xycemin(min,di.v1[isaved+1]);
            }

            if ((di.getModel().howToInterp == LTRA_MOD_LININTERP) || (isaved == 0) ||
                ((isaved != 0) && ((di.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                                   (di.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)) && ((v1d > max) || (v1d < min))))
            {
              if ((isaved != 0) && (di.getModel().howToInterp == LTRA_MOD_QUADINTERP))
              {
#ifdef Xyce_DEBUG_DEVICE
                std::cout << "[LTRA-DBG-DEV] load: warning: interpolated v1 is out of range after timepoint "
                          << getSolverState().ltraTimeIndex << std::endl;
                std::cout << "         values: "
                          << di.v1[isaved-1] << "  "
                          << di.v1[isaved] << "  "
                          << di.v1[isaved+1] << "; interpolated: "
                          << v1d << std::endl;
                std::cout << "        timepoints are: "
                          << getSolverState().currTime - di.getModel().td << std::endl;
#endif
              }
              else
              {
                v1d = di.v1[isaved] * lf2 + di.v1[isaved+1] * lf3;
              }
            }

            if ((isaved != 0) &&
                ((di.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                 (di.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)))
            {
              i1d = di.i1[isaved-1] * qf1
                + di.i1[isaved] * qf2
                + di.i1[isaved+1] * qf3;

              max = Xycemax(di.i1[isaved-1], di.i1[isaved]);
              max = Xycemax(max,di.i1[isaved+1]);
              min = Xycemin(di.i1[isaved-1], di.i1[isaved]);
              min = Xycemin(min,di.i1[isaved+1]);
            }

            if ((di.getModel().howToInterp == LTRA_MOD_LININTERP) || (isaved == 0) ||
                ((isaved != 0) && ((di.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                                   (di.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)) &&
                 ((i1d > max) || (i1d < min))))
            {

              if ((isaved != 0) && (di.getModel().howToInterp == LTRA_MOD_QUADINTERP))
              {
#ifdef Xyce_DEBUG_DEVICE
                std::cout << "[LTRA-DBG-DEV] load: warning: interpolated i1 is out of range after timepoint "
                          << getSolverState().ltraTimeIndex << std::endl;
                std::cout << "         values: "
                          << di.i1[isaved-1] << "  "
                          << di.i1[isaved] << "  "
                          << di.i1[isaved+1] << "; interpolated: "
                          << i1d << std::endl;
                std::cout << "        timepoints are: "
                          << getSolverState().currTime - di.getModel().td << std::endl;
#endif
              }
              else
              {
                i1d = di.i1[isaved] * lf2 + di.i1[isaved+1] * lf3;
              }
            }

            if ((isaved != 0) &&
                ((di.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                 (di.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)))
            {
              v2d = di.v2[isaved-1] * qf1
                + di.v2[isaved] * qf2
                + di.v2[isaved+1] * qf3;

              max = Xycemax(di.v2[isaved-1], di.v2[isaved]);
              max = Xycemax(max,di.v2[isaved+1]);
              min = Xycemin(di.v2[isaved-1], di.v2[isaved]);
              min = Xycemin(min,di.v2[isaved+1]);
            }

            if ((di.getModel().howToInterp ==
                 LTRA_MOD_LININTERP) || (isaved == 0) ||
                ((isaved != 0) &&
                 ((di.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                  (di.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)) &&
                 ((v2d > max) || (v2d < min))))
            {

              if ((isaved != 0) &&
                  (di.getModel().howToInterp == LTRA_MOD_QUADINTERP))
              {
#ifdef Xyce_DEBUG_DEVICE
                std::cout << "[LTRA-DBG-DEV] load: warning: interpolated v2 is out of range after timepoint "
                          << getSolverState().ltraTimeIndex << std::endl;
                std::cout << "         values: "
                          << di.v2[isaved-1] << "  "
                          << di.v2[isaved] << "  "
                          << di.v2[isaved+1] << "; interpolated: "
                          << v2d << std::endl;
                std::cout << "        timepoints are: "
                          << getSolverState().currTime - di.getModel().td << std::endl;
#endif
              }
              else
              {
                v2d = di.v2[isaved] * lf2
                  + di.v2[isaved+1] *
                  lf3;
              }
            }

            if ((isaved != 0) &&
                ((di.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                 (di.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)))
            {
              i2d = di.i2[isaved-1] * qf1
                + di.i2[isaved] * qf2
                + di.i2[isaved+1] * qf3;

              max = Xycemax(di.i2[isaved-1], di.i2[isaved]);
              max = Xycemax(max,di.i2[isaved+1]); min = Xycemin(di.i2[isaved-1], di.i2[isaved]);
              min = Xycemin(min,di.i2[isaved+1]);
            }

            if ((di.getModel().howToInterp == LTRA_MOD_LININTERP) || (isaved == 0) ||
                ((isaved != 0) && ((di.getModel().howToInterp == LTRA_MOD_QUADINTERP) ||
                                   (di.getModel().howToInterp == LTRA_MOD_MIXEDINTERP)) &&
                 ((i2d > max) || (i2d < min))))
            {
              if ((isaved != 0) && (di.getModel().howToInterp == LTRA_MOD_QUADINTERP))
              {
#ifdef Xyce_DEBUG_DEVICE
                std::cout << "[LTRA-DBG-DEV] load: warning: interpolated i2 is out of range after timepoint "
                          << getSolverState().ltraTimeIndex << std::endl;
                std::cout << "         values: "
                          << di.i2[isaved-1] << "  "
                          << di.i2[isaved] << "  "
                          << di.i2[isaved+1] << "; interpolated: "
                          << i2d << std::endl;
                std::cout << "        timepoints are: "
                          << getSolverState().currTime - di.getModel().td << std::endl;
#endif
              }
              else
              {
                i2d = di.i2[isaved] * lf2 + di.i2[isaved+1] * lf3;
              }
            }
          }

          /* interpolation done */
          break;

      case LTRA_MOD_RC:
        break;

      default:
        return false;
        // return(E_BADPARM);
      }

      switch (di.getModel().specialCase)
      {
        case LTRA_MOD_RLC:

          /* begin convolution parts */

          /* convolution of h1dash with v1 and v2 */
          /* the matrix has already been loaded above */

          dummy1 = dummy2 = 0.0;
          for (int j = getSolverState().ltraTimeIndex; j > 0; j--)
          {
            if (di.getModel().h1dashCoeffs[j] != 0.0)
            {
              dummy1 += di.getModel().h1dashCoeffs[j] * (di.v1[j] - di.initVolt1);
              dummy2 += di.getModel().h1dashCoeffs[j] * (di.v2[j] - di.initVolt2);
            }
          }

          dummy1 += di.initVolt1 * di.getModel().intH1dash;
          dummy2 += di.initVolt2 * di.getModel().intH1dash;

          dummy1 -= di.initVolt1 * di.getModel().h1dashFirstCoeff;
          dummy2 -= di.initVolt2 * di.getModel().h1dashFirstCoeff;

          di.input1 -= dummy1 * di.getModel().admit;
          di.input2 -= dummy2 * di.getModel().admit;

          /* end convolution of h1dash with v1 and v2 */

          /* convolution of h2 with i2 and i1 */

          dummy1 = dummy2 = 0.0;
          if (di.getModel().tdover)
          {
            /* the term for ckt->CKTtime - di.getModel().td */
            dummy1 = (i2d - di.initCur2)* di.getModel().h2FirstCoeff;
            dummy2 = (i1d - di.initCur1)* di.getModel().h2FirstCoeff;

            /* the rest of the convolution */

            for (int j= /*di.getModel().h2Index*/di.getModel().auxIndex; j > 0; j--)
            {

              if (di.getModel().h2Coeffs[j] != 0.0)
              {
                dummy1 += di.getModel().h2Coeffs[j] * (di.i2[j] - di.initCur2);
                dummy2 += di.getModel().h2Coeffs[j] * (di.i1[j] - di.initCur1);
              }
            }
          }

          /* the initial-condition terms */

          dummy1 += di.initCur2 * di.getModel().intH2;
          dummy2 += di.initCur1 * di.getModel().intH2;

          di.input1 += dummy1;
          di.input2 += dummy2;

          /* end convolution of h2 with i2 and i1 */
          /* convolution of h3dash with v2 and v1 */
          /* the term for ckt->CKTtime - di.getModel().td */

          dummy1 = dummy2 = 0.0;
          if (di.getModel().tdover)
          {
            dummy1 = (v2d - di.initVolt2)* di.getModel().h3dashFirstCoeff;
            dummy2 = (v1d - di.initVolt1)* di.getModel().h3dashFirstCoeff;

            /* the rest of the convolution */

            for (int j= /*di.getModel().h3dashIndex*/di.getModel().auxIndex; j > 0; j--)
            {
              if (di.getModel().h3dashCoeffs[j] != 0.0)
              {
                dummy1 += di.getModel().h3dashCoeffs[j] * (di.v2[j] - di.initVolt2);
                dummy2 += di.getModel().h3dashCoeffs[j] * (di.v1[j] - di.initVolt1);
              }
            }
          }

          /* the initial-condition terms */

          dummy1 += di.initVolt2 * di.getModel().intH3dash;
          dummy2 += di.initVolt1 * di.getModel().intH3dash;

          di.input1 += di.getModel().admit*dummy1;
          di.input2 += di.getModel().admit*dummy2;

          /* end convolution of h3dash with v2 and v1 */

          /* NOTE: this switch passes through to following case */

        case LTRA_MOD_LC:
          /* begin lossless-like parts */

          if (!di.getModel().tdover)
          {
            di.input1 += di.getModel().attenuation * (di.initVolt2*di.getModel().admit + di.initCur2);
            di.input2 += di.getModel().attenuation * (di.initVolt1*di.getModel().admit + di.initCur1);
          }
          else
          {
            di.input1 += di.getModel().attenuation * (v2d*di.getModel().admit + i2d);
            di.input2 += di.getModel().attenuation * (v1d*di.getModel().admit + i1d);
          }

          // Residuals for the internal equations. These are for both
          // the RLC and LC case.
          fVec[di.li_Ibr1] += ((di.getModel().admit * (di.getModel().h1dashFirstCoeff + 1.0)) *
                               (di.vpos1-di.vneg1) - di.currp1) - di.input1;

          fVec[di.li_Ibr2] += ((di.getModel().admit * (di.getModel().h1dashFirstCoeff + 1.0)) *
                               (di.vpos2-di.vneg2) - di.currp2) - di.input2;

          /* end lossless-like parts */
          break;

        case LTRA_MOD_RC:

          /* begin convolution parts */

          /* convolution of h1dash with v1 and v2 */
          /* the matrix has already been loaded above */

          dummy1 = 0.0;
          dummy2 = 0.0;
          for (int j = getSolverState().ltraTimeIndex; j > 0; j--)
          {
            if (di.getModel().h1dashCoeffs[j]!= 0.0)
            {
              dummy1 += di.getModel().h1dashCoeffs[j] * (di.v1[j] - di.initVolt1);
              dummy2 += di.getModel().h1dashCoeffs[j] * (di.v2[j] - di.initVolt2);
            }
          }

          /* the initial condition terms */

          dummy1 += di.initVolt1 * di.getModel().intH1dash; dummy2 += di.initVolt2 * di.getModel().intH1dash;

          /* the constant contributed by the init
           * condition and the latest timepoint */

          dummy1 -= di.initVolt1* di.getModel().h1dashFirstCoeff;
          dummy2 -= di.initVolt2* di.getModel().h1dashFirstCoeff;

          di.input1 -= dummy1;
          di.input2 -= dummy2;

          /* end convolution of h1dash with v1 and v2 */
          /* convolution of h2 with i2 and i1 */

          dummy1=dummy2=0.0;

          for (int j=getSolverState().ltraTimeIndex; j > 0; j--)
          {
            if (di.getModel().h2Coeffs[i] != 0.0)
            {
              dummy1 += di.getModel().h2Coeffs[j] * (di.i2[j] - di.initCur2);
              dummy2 += di.getModel().h2Coeffs[j] * (di.i1[j] - di.initCur1);
            }
          }

          /* the initial-condition terms */
          dummy1 += di.initCur2 * di.getModel().intH2;
          dummy2 += di.initCur1 * di.getModel().intH2;

          dummy1 -= di.initCur2* di.getModel().h2FirstCoeff;
          dummy2 -= di.initCur1* di.getModel().h2FirstCoeff;

          di.input1 += dummy1;
          di.input2 += dummy2;

          /* end convolution of h2 with i2 and i1 */
          /* convolution of h3dash with v2 and v1 */

          dummy1 = dummy2 = 0.0;

          for (int j=getSolverState().ltraTimeIndex; j > 0; j--)
          {
            if (di.getModel().h3dashCoeffs[j] != 0.0)
            {
              dummy1 += di.getModel().h3dashCoeffs[j] * (di.v2[j] - di.initVolt2);
              dummy2 += di.getModel().h3dashCoeffs[j] * (di.v1[j] - di.initVolt1);
            }
          }

          /* the initial-condition terms */

          dummy1 += di.initVolt2 * di.getModel().intH3dash;
          dummy2 += di.initVolt1 * di.getModel().intH3dash;

          dummy1 -= di.initVolt2* di.getModel().h3dashFirstCoeff;
          dummy2 -= di.initVolt1* di.getModel().h3dashFirstCoeff;

          di.input1 += dummy1;
          di.input2 += dummy2;

          // Residuales for the internal equations.
          fVec[di.li_Ibr1] += ((di.getModel().h1dashFirstCoeff * (di.vpos1-di.vneg1) -
                                di.getModel().h3dashFirstCoeff * (di.vpos2-di.vneg2) -
                                di.getModel().h2FirstCoeff * di.currp2 -
                                di.currp1) -
                               di.input1);

          fVec[di.li_Ibr2] += ((di.getModel().h1dashFirstCoeff * (di.vpos2-di.vneg2) -
                                di.getModel().h3dashFirstCoeff * (di.vpos1-di.vneg1) -
                                di.getModel().h2FirstCoeff * di.currp1 -
                                 di.currp2) -
                               di.input2);

          /* end convolution of h3dash with v2 and v1 */

          break;

        default:
          return false;
            //return(E_BADPARM);
      }

      // Residuals (KCL) for "normal" nodes and common between all cases
      fVec[di.li_Pos1] += di.currp1;
      fVec[di.li_Neg1] -= di.currp1;

      fVec[di.li_Pos2] += di.currp2;
      fVec[di.li_Neg2] -= di.currp2;

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
// Creation Date : 06/16/10
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix& dFdx, N_LAS_Matrix& dQdx)
{
  double dummy1(0.0), dummy2(0.0);
  ostringstream msg;

  // this is commented out for now because the loop contains return statements which
  // break the OMP threading on the loop -- RLS 8/20/2010
  // #ifdef _OMP
  // #pragma omp parallel for
  // #endif
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance& di = *(*it);

    if( getSolverState().dcopFlag || (di.getModel().specialCase == LTRA_MOD_RG))
    {
      switch (di.getModel().specialCase)
      {
        case LTRA_MOD_RG:
          *(di.ibr1Pos1Ptr) +=  1.0;
          *(di.ibr1Neg1Ptr) += -1.0;
          *(di.ibr1Pos2Ptr) += -di.getModel().coshlrootGR;
          *(di.ibr1Neg2Ptr) +=  di.getModel().coshlrootGR;
          *(di.ibr1Ibr2Ptr) +=  (1.0 + getDeviceOptions().gmin) * di.getModel().rRsLrGRorG;

          *(di.ibr2Ibr2Ptr) +=  di.getModel().coshlrootGR;
          *(di.ibr2Pos2Ptr) += -(1.0 + getDeviceOptions().gmin) * di.getModel().rGsLrGRorR;
          *(di.ibr2Neg2Ptr) +=  (1.0 + getDeviceOptions().gmin) * di.getModel().rGsLrGRorR;
          *(di.ibr2Ibr1Ptr) +=  1.0;

          *(di.pos1Ibr1Ptr) +=  1.0;
          *(di.neg1Ibr1Ptr) += -1.0;
          *(di.pos2Ibr2Ptr) +=  1.0;
          *(di.neg2Ibr2Ptr) += -1.0;

          break;

        case LTRA_MOD_LC:
        case LTRA_MOD_RLC:
        case LTRA_MOD_RC: // load a simple resistor

          *(di.pos1Ibr1Ptr) +=  1.0;
          *(di.neg1Ibr1Ptr) += -1.0;

          *(di.pos2Ibr2Ptr) +=  1.0;
          *(di.neg2Ibr2Ptr) += -1.0;

          *(di.ibr1Ibr1Ptr) +=  1.0;
          *(di.ibr1Ibr2Ptr) +=  1.0;

          *(di.ibr2Pos1Ptr) +=  1.0;
          *(di.ibr2Pos2Ptr) += -1.0;
          *(di.ibr2Ibr1Ptr) += -di.getModel().resist*di.getModel().length;

          break;

        default:
          msg.clear();
          msg << "***************" << std::endl;
          msg << ": Error. Unknown LTRA configuration, ";
          msg << di.getModel().specialCase;
          msg << ". Must be one of RG, LC, RLC or RC.";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR, msg.str());

          return false;
      }

    }
    else
    {
      /* all cases other than DC or the RG case */

      /* matrix loading - done every time load is called */
      switch (di.getModel().specialCase)
      {
        case LTRA_MOD_RLC:
          /* loading for convolution parts' first terms */

          dummy1 = di.getModel().admit * di.getModel().h1dashFirstCoeff;

          *(di.ibr1Pos1Ptr) += dummy1;
          *(di.ibr1Neg1Ptr) -= dummy1;

          *(di.ibr2Pos2Ptr) += dummy1;
          *(di.ibr2Neg2Ptr) -= dummy1;
          /* end loading for convolution parts' first terms */

          /* NOTE: This case intentionally falls through to the next case */

        case LTRA_MOD_LC:
          /* this section loads for the parts of the equations that
             resemble the lossless equations */

          *(di.ibr1Pos1Ptr) += di.getModel().admit;
          *(di.ibr1Neg1Ptr) -= di.getModel().admit;

          *(di.ibr1Ibr1Ptr) -= 1.0;

          *(di.pos1Ibr1Ptr) += 1.0;
          *(di.neg1Ibr1Ptr) -= 1.0;

          *(di.ibr2Pos2Ptr) += di.getModel().admit;
          *(di.ibr2Neg2Ptr) -= di.getModel().admit;

          *(di.ibr2Ibr2Ptr) -= 1.0;

          *(di.pos2Ibr2Ptr) += 1.0;
          *(di.neg2Ibr2Ptr) -= 1.0;

          /* loading for lossless-like parts over */
          break;

        case LTRA_MOD_RC:

          /* this section loads for the parts of the equations that
             have no convolution */

          *(di.ibr1Ibr1Ptr) -= 1.0;

          *(di.pos1Ibr1Ptr) += 1.0;
          *(di.neg1Ibr1Ptr) -= 1.0;

          *(di.ibr2Ibr2Ptr) -= 1.0;

          *(di.pos2Ibr2Ptr) += 1.0;
          *(di.neg2Ibr2Ptr) -= 1.0;

          /* loading for non-convolution parts over */
          /* loading for convolution parts' first terms */

          dummy1 = di.getModel().h1dashFirstCoeff;

          *(di.ibr1Pos1Ptr) += dummy1;
          *(di.ibr1Neg1Ptr) -= dummy1;

          *(di.ibr2Pos2Ptr) += dummy1;
          *(di.ibr2Neg2Ptr) -= dummy1;

          dummy1 = di.getModel().h2FirstCoeff;

          *(di.ibr1Ibr2Ptr) -= dummy1;
          *(di.ibr2Ibr1Ptr) -= dummy1;

          dummy1 = di.getModel().h3dashFirstCoeff;

          *(di.ibr1Pos2Ptr) -= dummy1;
          *(di.ibr1Neg2Ptr) += dummy1;

          *(di.ibr2Pos1Ptr) -= dummy1;
          *(di.ibr2Neg1Ptr) += dummy1;

          /* end loading for convolution parts' first terms */

          break;

        default:
          return false;
      }
    }
  }

  return true;
}

} // namespace LTRA
} // namespace Device
} // namespace Xyce
