//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_Synapse4.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Christina Warrender, SNL, Cognitive Modeling
//
// Creation Date  : 10/12/11
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.29.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:34 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#include <Xyce_config.h>
//#define Xyce_FullSynapseJac 1

// ---------- Standard Includes ----------
// used to get time in seconds to seed random number generator.
#include<time.h>

// ----------   Xyce Includes   ----------
//
#include <N_DEV_Const.h>
#include <N_DEV_Synapse4.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Synapse4::Instance>::ParametricData()
{
  // Set up configuration constants:
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(0);
  addModelType("SYNAPSE");

  // Set up map for double precision variables:
  addPar ("GMAX", 0.01, false, ParameterType::NO_DEP,
          &Synapse4::Instance::gMax,
          &Synapse4::Instance::gMaxGiven,
          U_OHMM1, CAT_NONE, "Maximal Synaptic Conductance");
  
}

template<>
ParametricData<Synapse4::Model>::ParametricData()
{
  // Set up map for double precision variables:
  addPar ("VTHRESH", 0.01, false, ParameterType::NO_DEP,
          &Synapse4::Model::vThresh, NULL,
          U_VOLT, CAT_NONE, "Presynaptic voltage spike threhsold");

  addPar ("DELAY", 0.001, false, ParameterType::NO_DEP,
          &Synapse4::Model::delay, NULL,
          U_SECOND, CAT_NONE, "Time delay between presynaptic signal and postsynaptic response");

  addPar ("GMAX", 0.01, false, ParameterType::NO_DEP,
          &Synapse4::Model::gMax, NULL,
          U_OHMM1, CAT_NONE, "Maximal Synaptic Conductance");

  addPar ("EREV", 0.0, false, ParameterType::NO_DEP,
          &Synapse4::Model::eRev, NULL,
          U_VOLT, CAT_NONE, "Postsynaptic Reversal Potential");

  addPar ("TAU1", 0.0001, false, ParameterType::NO_DEP,
          &Synapse4::Model::tau1, NULL,
          U_SECM1, CAT_NONE, "Rise time constant");

  addPar ("TAU2", 0.01, false, ParameterType::NO_DEP,
          &Synapse4::Model::tau2, NULL,
          U_SECM1, CAT_NONE, "Decay time constant");
}

namespace Synapse4 {

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
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{
  // initialization
  respondTime = numeric_limits<double>::max( );
  ready = true;
  //active = false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
Instance::Instance(
  InstanceBlock & IB,
  Model & Riter,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)
  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(Riter),
    li_Prev(-1),
    li_Post(-1),
    li_A0_store(-1),
    li_B0_store(-1),
    li_t0_store(-1),
    li_store_dev_i(-1),
#ifdef Xyce_FullSynapseJac
    APostEquPostNodeOffset(-1),
    f_PostEquPostNodePtr(0),
#endif
    ipost(0),
    didVpost(0),
    randInitialized(false)
{
  numIntVars   = 0;   // A and B   2
  numExtVars   = 2;   // presynaptic V and postsynaptic V
  setNumStoreVars(3);   // A0, B0, t0
  numLeadCurrentStoreVars = 1;

  setName(IB.getName());
  setModelName(model_.getName());

  if( jacStamp.empty() )
  {
    jacStamp.resize(2);
    jacStamp[0].resize(0);    // presynaptic V not changed
#ifdef Xyce_FullSynapseJac
    jacStamp[1].resize(1);    // postsynaptic V depends on itself
#else
    jacStamp[1].resize(0);    // postsynaptic V depends on itself
#endif
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  if (!gMaxGiven )
  {
    gMax = model_.gMax;
  }

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
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
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
    cout << "  SynapseInstance::registerLIDs" << endl;
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

  li_Prev = extLIDVec[0];
  li_Post = extLIDVec[1];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << "  li_Prev = " << li_Prev << endl;
    cout << "  li_Post = " << li_Post << endl;
  }
#endif

  /*
    li_AVar = intLIDVec[0];
    li_BVar = intLIDVec[1];
  */

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
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 10/17/11
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  /*
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
  string tmpstr;

  tmpstr = name + "_A";
  spiceInternalName (tmpstr);
  intNameMap[ li_AVar ] = tmpstr;
  tmpstr = name + "_B";
  spiceInternalName (tmpstr);
  intNameMap[ li_BVar ] = tmpstr;
  }
  */
  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const vector<int> & stoLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSto = stoLIDVecRef.size();

  if (numSto != getNumStoreVars())
  {
    msg = "Instance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
  // copy over the global ID lists.
  stoLIDVec = stoLIDVecRef;

  li_A0_store = stoLIDVec[0];
  li_B0_store = stoLIDVec[1];
  li_t0_store = stoLIDVec[2];
  if( loadLeadCurrent )
  {
    li_store_dev_i = stoLIDVec[3];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 08/01/2012
//-----------------------------------------------------------------------------
map<int,string> & Instance::getStoreNameMap()
{
  // set up the internal name map, if it hasn't been already.
  if (stateNameMap.empty ())
  {
    string baseString(getName() + "_");
    string tempString;

    tempString = baseString + "A0";
    storeNameMap[ li_A0_store ] = tempString;
    tempString = baseString + "B0";
    storeNameMap[ li_B0_store ] = tempString;
    tempString = baseString + "T0";
    storeNameMap[ li_t0_store ] = tempString;
    if( loadLeadCurrent )
    {
      tempString = getName() + ":DEV_I";
      storeNameMap[ li_store_dev_i ] = tempString;
    }
  }
  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
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
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
#ifdef Xyce_FullSynapseJac
  APostEquPostNodeOffset = jacLIDVec[1][0];
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
#ifdef Xyce_FullSynapseJac
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  f_PostEquPostNodePtr = &(dFdx[li_Post][APostEquPostNodeOffset]);
#endif
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one synapse instance
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;

  double * lastSolVecPtr = extData.lastSolVectorRawPtr;
  double * stoVector = extData.nextStoVectorRawPtr;
  double * lastStoVec = extData.lastStoVectorRawPtr;

  // initialized random number generator if needed
  if( !randInitialized )
  {
    if( getDeviceOptions().randomSeed != 0 )
    {
      devSupport.SetSeed( getDeviceOptions().randomSeed );
    }
    else
    {
      unsigned int aSeed = static_cast<unsigned int>(time( NULL ) ) ;
      devSupport.SetSeed( aSeed );
    }
    randInitialized=true;
  }

  // This check need to
  // Check for presynaptic spike start, set time to respond
  double vPre  = lastSolVecPtr[li_Prev];
  double vPost = lastSolVecPtr[li_Post];
  double time = getSolverState().currTime;
  double vThresh = model_.vThresh;
  double delay = model_.delay;
  if (ready)
  {
    if (vPre > vThresh)
    {
      ready=false;
      respondTime = time + delay;
    }
  }
  else  // already had spike start, looking for end
  {
    if (vPre < vThresh)
    {
      ready=true;
    }
  }

  // only need to update state variables and synaptic current when synapse active
  // disabling testing for this because determining the point at which the synaptic current
  // becomes negligible is somewhat problematic, and there don't appear to be significant savings
  // from avoiding loads when synapse is inactive
  //if (active)
  {
    // handle decay of A and B
    double tau1 = model_.tau1;
    double tau2 = model_.tau2;

    double t0 = stoVector[li_t0_store];
    double Anow = stoVector[li_A0_store] * exp( - ( time - t0 ) / tau1 );
    double Bnow =  stoVector[li_B0_store] * exp( - ( time - t0 ) / tau2 );

    // set up variables for load methods
    // current equation is the same whether we're responding to spike or not,
    // assuming current A and B values set appropriately above
    // ipost = (B-A)*(V-Erev)

    double eRev = model_.eRev;
    double cond = Bnow-Anow;
    ipost = cond*(vPost-eRev);
    didVpost = cond;

#ifdef Xyce_DEBUG_DEVICE
    const string dashedline =
      "-------------------------------------------------------------------------"
      "----";
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << endl << dashedline << endl;
      cout << "  SynapseInstance::updateIntermediateVars" << endl;
      cout << "Anow:  " << Anow << endl;
      cout << "Bnow:  " << Bnow << endl;
      cout << "vPost:  " << vPost << endl;
      cout << "eRev:  " << eRev << endl;
      cout << "ipost:  " << ipost << endl;
      cout << "didVpost:  " << didVpost << endl;
    }
#endif

  }	// end if synapse active

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = updateIntermediateVars();
  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState()
{
  return  true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Synapse4
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
//                 The "Q" vector contains charges and fluxes, mostly.
//                 However, it is ordered like the solution vector, and as
//                 it is part of the KCL formulation, the terms in Q will
//                 actually be *sums* of charges, rather than single
//                 distinct charges.
//
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  return true;
}
//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Synapse4  instance.
//
// Special Notes : This is an algebraic constaint, and as such the Synapse4
//                 does make a contribution to it.
//
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  N_LAS_Vector *  fVecPtr = extData.daeFVectorPtr;
  (*fVecPtr)[li_Prev] += 0.0;
  (*fVecPtr)[li_Post] += ipost;
  if( loadLeadCurrent )
  {
    double * stoVec = extData.nextStoVectorRawPtr;
    stoVec[li_store_dev_i] = ipost;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 Synapse4 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
//                 The "Q" vector contains charges and fluxes, mostly.
//                 However, it is ordered like the solution vector, and as
//                 it is part of the KCL formulation, the terms in Q will
//                 actually be *sums* of charges, rather than single
//                 distinct charges.
//
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Synapse4  instance.
//
// Special Notes : This is an algebraic constaint, and as such the Synapse4
//                 does make a contribution to it.
//
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
#ifdef Xyce_FullSynapseJac
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Post][APostEquPostNodeOffset] += didVpost;
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : InstanceInstance::outputPlotFiles
// Purpose       : If requested by the user output all the variables
//                 associated with the population
// Special Notes : We're actually using this method not for output, but because
//                 it's called after the system has converged.  In this case,
//                 that let us mark the end of handling a presynaptic event.
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 10/25/2011
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles ()
{
  bool bsuccess = true;

  // cew 11/3/11:  changing this back to just storing A0 and B0 for updateIntermediateVars
  // to use in calculating A, B, and ipost
  // But when incrementing A0 (B0), current value of A (B) must be used.

  double time = getSolverState().currTime;
  //N_LAS_Vector * staVector = extData.nextStaVectorPtr;
  double * stoVector = extData.nextStoVectorRawPtr;
  if (time >= respondTime)
  {
    //std::cout << "Instance::outputPlotFiles() adjusting A0, B0 and t0 = " << getSolverState().currTime << std::endl;
    // succesfully processed a step, so adjust the next
    // respondTime to distant future
    respondTime = numeric_limits<double>::max( );
    //active = true;

    // now we also update the A0, B0 and t0 in the state vector
    double factor = model_.factor;
    double deltaAB = factor*gMax;
    double tau1 = model_.tau1;
    double tau2 = model_.tau2;
    double t0 = stoVector[li_t0_store];
    double Anow = stoVector[li_A0_store] * exp( - ( time - t0 ) / tau1 );
    double Bnow =  stoVector[li_B0_store] * exp( - ( time - t0 ) / tau2 );
    stoVector[li_A0_store] = Anow + deltaAB;
    stoVector[li_B0_store] = Bnow + deltaAB;
    stoVector[li_t0_store] = getSolverState().currTime;
    //std::cout << "A:  " << (*staVector)[li_A0_state] << " B: " << (*staVector)[li_B0_state] << endl;
  }

/* disabling this because determining the point at which the synaptic current becomes negligible
   is somewhat problematic, and there don't appear to be significant savings from avoiding loads
   when synapse is inactive
   if (active)
   {
   // determine when we can stop updating synaptic current, when it's negligible
   // can't use magnitude of current, because it starts off small initially
   // instead, determine whether time since last spike is greater than reasonable multiple of
   // larger time constant
   if ( (time-(*staVector)[li_t0_state])>3.0*model_.maxtau )
   {
   active = false;
   ipost = 0.0;
   didVpost = 0.0;
   }
   }
*/

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{
  // initialize variables needed to calculate synaptic dynamics
  if (tau1/tau2 > .9999) {
    tau1 = .9999*tau2;
  }
  tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1);
  factor = -exp(-tp/tau1) + exp(-tp/tau2);
  factor = 1/factor;
  maxtau = (tau1>tau2)?tau1:tau2;

  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 10/12/11
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
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
              SolverState & ss1,
              DeviceOptions & do1)
  : DeviceModel(MB, ss1,do1),
    vThresh(0.0),
    delay(0.0),
    gMax(0.0),
    eRev(0.0),
    tau1(0.0),
    tau2(0.0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
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
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;
  isize = instanceContainer.size();
  os << endl;
  os << "Number of Synapse4 Instances: " << isize << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 07/18/12
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  if( getSolverState().dcopFlag )
  {
    // no synaptic activity during dcop
    return true;
  }

  bool bsuccess = true;

#ifdef _OMP
#pragma omp parallel for reduction (&&:bsuccess)
#endif

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    bool tmpBool = (*it)->updateIntermediateVars();	// skipping call to updatePrimaryState
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 7/18/12
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState (double * staDerivVec, double *stoVec)
{
  // not used for this device
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 07/16/12
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors(double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ)
{
  if( getSolverState().dcopFlag )
  {
    // no synaptic activity during dcop
    return true;
  }

  bool bsuccess = true;

  // We don't need Q loads for this synapse device at all
  // Only need to do F loads if there's a nonzero synaptic current
  // disabling testing for the latter because determining the point at which the synaptic current
  // becomes negligible is somewhat problematic, and there don't appear to be significant savings
  // from avoiding loads when synapse is inactive

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    bool tmpBool = (*it)->loadDAEFVector();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 07/16/12
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  if( getSolverState().dcopFlag )
  {
    // no synaptic activity during dcop
    return true;
  }

  bool bsuccess = true;

  // We don't need Q loads for this synapse device at all
  // Only need to do F loads if there's a nonzero synaptic current
  // disabling testing for the latter because determining the point at which the synaptic current
  // becomes negligible is somewhat problematic, and there don't appear to be significant savings
  // from avoiding loads when synapse is inactive

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    bool tmpBool = (*it)->loadDAEdFdx();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

} // namespace Synapse4
} // namespace Device
} // namespace Xyce
