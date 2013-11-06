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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Neuron7.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical and Microsytem Modeling
//
// Creation Date  : 06/10/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_Neuron7.h>
#include <N_DEV_Neuron_CommonEquations.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Neuron7::Instance>::ParametricData()
{
  setNumNodes(1);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(1);
  setPrimaryParameter("");
  addModelType("NEURON");

  // Set up map for normal (double) param variables:
  addPar ("MEMC", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::memCap,
          &Neuron7::Instance::memCapGiven,
          U_FARAD, CAT_NONE, "Membrane capacitance");

  addPar ("VT", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::Vt,
          &Neuron7::Instance::VtGiven,
          U_VOLT, CAT_NONE, "Instantaneous threshold voltage");

  addPar ("VR", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::Vr,
          &Neuron7::Instance::VrGiven,
          U_VOLT, CAT_NONE, "Resting membrane potential");

  addPar ("VP", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::Vpeak,
          &Neuron7::Instance::VpeakGiven,
          U_VOLT, CAT_NONE, "Peak voltage");

  addPar ("K", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::k,
          &Neuron7::Instance::kGiven,
          U_NONE, CAT_NONE, "instanceing parameter");

  addPar ("A", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::a,
          &Neuron7::Instance::aGiven,
          U_NONE, CAT_NONE, "instanceing parameter");

  addPar ("B", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::b,
          &Neuron7::Instance::bGiven,
          U_NONE, CAT_NONE, "instanceing parameter");

  addPar ("C", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::c,
          &Neuron7::Instance::cGiven,
          U_NONE, CAT_NONE, "instanceing parameter");

  addPar ("D", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Instance::d,
          &Neuron7::Instance::dGiven,
          U_NONE, CAT_NONE, "instanceing parameter");

  addPar ("USCALE", 1.0e-9, false, ParameterType::NO_DEP,
          &Neuron7::Instance::uscale,
          &Neuron7::Instance::uscaleGiven,
          U_NONE, CAT_NONE, "scaling for u variable");

  addPar ("FALLRATE", 1.0e3, false, ParameterType::NO_DEP,
          &Neuron7::Instance::fallRate,
          &Neuron7::Instance::fallRateGiven,
          U_NONE, CAT_NONE, "recovery rate");
}

template<>
ParametricData<Neuron7::Model>::ParametricData()
{
  addPar ("MEMC", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::memCap,
          &Neuron7::Model::memCapGiven,
          U_FARAD, CAT_NONE, "Membrane capacitance");

  addPar ("VT", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::Vt,
          &Neuron7::Model::VtGiven,
          U_VOLT, CAT_NONE, "Instantaneous threshold voltage");

  addPar ("VR", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::Vr,
          &Neuron7::Model::VrGiven,
          U_VOLT, CAT_NONE, "Resting membrane potential");

  addPar ("VP", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::Vpeak,
          &Neuron7::Model::VpeakGiven,
          U_VOLT, CAT_NONE, "Peak voltage");

  addPar ("K", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::k,
          &Neuron7::Model::kGiven,
          U_NONE, CAT_NONE, "Neuron7::Modeling parameter");

  addPar ("A", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::a,
          &Neuron7::Model::aGiven,
          U_NONE, CAT_NONE, "Neuron7::Modeling parameter");

  addPar ("B", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::b,
          &Neuron7::Model::bGiven,
          U_NONE, CAT_NONE, "Neuron7::Modeling parameter");

  addPar ("C", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::c,
          &Neuron7::Model::cGiven,
          U_NONE, CAT_NONE, "Neuron7::Modeling parameter");

  addPar ("D", 0.0, false, ParameterType::NO_DEP,
          &Neuron7::Model::d,
          &Neuron7::Model::dGiven,
          U_NONE, CAT_NONE, "Neuron7::Modeling parameter");

  addPar ("USCALE", 1.0e-9, false, ParameterType::NO_DEP,
          &Neuron7::Model::uscale,
          &Neuron7::Model::uscaleGiven,
          U_NONE, CAT_NONE, "scaling for u variable");

  addPar ("FALLRATE", 1.0e3, false, ParameterType::NO_DEP,
          &Neuron7::Model::fallRate,
          &Neuron7::Model::fallRateGiven,
          U_NONE, CAT_NONE, "recovery rate");
}

namespace Neuron7 {

//
// static class member inits
//
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
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
                   Model & Miter,
                   MatrixLoadData & mlData1,
                   SolverState &ss1,
                   ExternData  &ed1,
                   DeviceOptions & do1)
  : DeviceInstance (IB, mlData1, ss1, ed1, do1),
    model_(Miter),
    memCap(0.0),
    Vt(0.0),
    Vr(0.0),
    Vpeak(0.0),
    k(0.0),
    a(0.0),
    b(0.0),
    c(0.0),
    d(0.0),
    uscale(1.0e-9),
    fallRate(1.0e3),
    memCapGiven(false),
    VtGiven(false),
    VrGiven(false),
    VpeakGiven(false),
    kGiven(false),
    aGiven(false),
    bGiven(false),
    cGiven(false),
    dGiven(false),
    uscaleGiven(false),
    fallRateGiven(false),
    vEquFvalue(0.0),
    vEquQvalue(0.0),
    vEqudFdv(0.0),
    vEqudFdu(0.0),
    vEqudQdv(0.0),
    uEquFvalue(0.0),
    uEquQvalue(0.0),
    uEqudFdv(0.0),
    uEqudFdu(0.0),
    uEqudQdu(0.0),
    resetting(false),
    uPeak(0.0),
    li_V(-1),
    li_U(-1),
    vEquVOffset(-1),
    vEquUOffset(-1),
    uEquVOffset(-1),
    uEquUOffset(-1)
{
  setName(IB.getName());
  setModelName(model_.getName());

  numExtVars = 1;  // membrane voltage

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  //if (!given("TEMP"))
  //  temp = getDeviceOptions().temp.dVal();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // pull unspecified params out of the model if they weren't specified here
  if( !memCapGiven && model_.memCapGiven )
  {
    memCap = model_.memCap;
    memCapGiven = true;
  }
  if( !VtGiven && model_.VtGiven )
  {
    Vt = model_.Vt;
    VtGiven = true;
  }
  if( !VrGiven && model_.VrGiven )
  {
    Vr = model_.Vr;
    VrGiven = true;
  }
  if( !VpeakGiven && model_.VpeakGiven )
  {
    Vpeak = model_.Vpeak;
    VpeakGiven = true;
  }
  if( !kGiven && model_.kGiven )
  {
    k = model_.k;
    kGiven = true;
  }
  if( !aGiven && model_.aGiven )
  {
    a = model_.a;
    aGiven = true;
  }
  if( !bGiven && model_.bGiven )
  {
    b = model_.b;
    bGiven = true;
  }
  if( !cGiven && model_.cGiven )
  {
    c = model_.c;
    cGiven = true;
  }
  if( !dGiven && model_.dGiven )
  {
    d = model_.d;
    dGiven = true;
  }
  if( !uscaleGiven && model_.uscaleGiven )
  {
    uscale = model_.uscale;
    uscaleGiven = true;
  }
  if( !fallRateGiven && model_.fallRateGiven )
  {
    fallRate = model_.fallRate;
    fallRateGiven = true;
  }

  /*
    std::cout << "Instance::Instance" << std::endl
    << "memC  = " << memCap << std::endl
    << "Vt    = " << Vt << std::endl
    << "Vr    = " << Vr << std::endl
    << "Vp    = " << Vpeak << std::endl
    << "k     = " << k << std::endl
    << "a     = " << a << std::endl
    << "b     = " << b << std::endl
    << "c     = " << c << std::endl
    << "d     = " << d << std::endl
    << "uscale = " << uscale << std::endl
    << "fallRate = " << fallRate << std::endl;
  */

  numIntVars = 1;  // u
  numStateVars = 0;

  // total up number of vars.
  int numVars = numExtVars + numIntVars;


  //
  // equations are:
  // v equation
  //  k * (v-Vr) * (v-Vt) - u + I - C dv/dt = 0
  // u euqtion
  //  a * ( b * (v-Vr) - u ) - du/dt = 0

  //  df/dx
  //            V                u
  // Vequ  2k V - Vt - Vr       -1
  // Uequ      a b              -a
  //
  // dq/dx
  //            V                u
  // Vequ    -memCap             0
  // Uequ       0               -1
  //
  // so jacStamp is dense.


  // set up jacStamp
  if( jacStamp.empty() )
  {
    jacStamp.resize(2);
    jacStamp[0].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1].resize(2);
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
  }

  /*
  // print out jacStamp
  std::cout << "jacStamp for Neuron6" << std::endl;
  int numRows = jacStamp.size();
  for( int i=0; i< numRows; i++ )
  {
  int numCol = jacStamp[i].size();
  std::cout << "jacStamp[ " << i << " ] = { ";
  for(int j=0; j<numCol; j++)
  {
  std::cout << jacStamp[i][j] << "  ";
  }
  std::cout << " } " <<  std::endl;
  }
  std::cout << std::endl;
  */

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::processParams(string param)
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  //updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const vector<int> & intLIDVecRef,
                            const vector<int> & extLIDVecRef)
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "-------------------------------------------------------------------------"
    "----";
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << dashedline << endl;
    cout << "  Instance::registerLIDs" << endl;
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

  li_V  = extLIDVec[0];
  li_U  = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << "  li_V = " << li_V << endl
         << "  li_U = " << li_U << endl;

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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    string tmpstr;
    tmpstr = getName() + "_" + "U" ;
    spiceInternalName (tmpstr);
    intNameMap[ li_U ] = tmpstr;
  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
  // no state vars so this is a no op
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDeviceMask
//
// Purpose       : Loads the zero elements of the device mask
//
// Special Notes : elements of the error vector associated with zero
//                 elements of the mask will not be included in weighted
//                 norms by the time integrator.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDeviceMask ()
{
  // no masking so return false
  return false;

  //bool returnVal=false;
  //N_LAS_Vector * maskVectorPtr = extData.deviceMaskVectorPtr;
  //
  // return (returnVal);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  vEquVOffset = jacLIDVec[0][0];
  vEquUOffset = jacLIDVec[0][1];
  uEquVOffset = jacLIDVec[1][0];
  uEquUOffset = jacLIDVec[1][1];
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;

  // here we take the current solutions for V1, V2, n, m and h
  // and use those to calculate all the terms needed for the next
  // load cycle (F, Q, dFdX, dQdX)
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  double vVal = (*solVectorPtr)[li_V];
  double uVal = (*solVectorPtr)[li_U];

  // not sure if I need to do this as it should just
  // fall out of the solution when d/dt = 0.  Possibly
  // it doesn't because there are two solutions when d/dt = 0,
  // either vVal = Vr or vVal = Vt we want vVal=Vr.
  if( getSolverState().dcopFlag  )
  {
    // vVal - Vr = 0 and uVal=0
    vEquFvalue = vVal - Vr;
    vEquQvalue = - memCap * vVal;
    vEqudFdv   = 1.0;
    vEqudFdu   = 0.0;
    vEqudQdv   = - memCap;
    uEquFvalue = uVal*uscale;
    uEquQvalue = -uVal*uscale;
    uEqudFdv   = 0.0;
    uEqudFdu   = 1.0*uscale;
    uEqudQdu   = uscale;
  }
  else
  {
    // need to cacth case where v> Vpeak
    if( !resetting && (vVal >= Vpeak) )
    {
      resetting = true;
      uPeak = uVal;
      // ask the device manager to call a breakpoint at this step (discontinuity step)
      //extData.devMgrPtr->declareCurrentStepAsBreakpoint();
    }
    else
    {
      // only break out of the resetting phase at the start of a newton solve
      if( (getSolverState().newtonIter == 0 ) && (vVal <= (c+0.0001) ) && (uVal >= (uPeak + (d-1.0e-13)/uscale) ))
      {
        // reset is done, go back to old form
        resetting = false;
      }
    }

    if ( resetting )
    {
      //extData.devMgrPtr->declareCurrentStepAsBreakpoint();
#if 0
      std::cout << "In resetting section uPeak + d = " << (uPeak+d*uscale) << " u - (uP +d) = " << (uVal - (uPeak + d*uscale)) << std::endl;
      // in this case vVal - c = 0 and uVal = uVal + d
      vEquFvalue = vVal - c;
      vEquQvalue = 0.0;
      vEqudFdv   = 1.0;
      vEqudFdu   = 0.0;
      vEqudQdv   = 0.0;
      uEquFvalue = uscale * (uVal - (uPeak + d*uscale));
      uEquQvalue = 0.0;
      uEqudFdv   = 0.0;
      uEqudFdu   = uscale;
      uEqudQdu   = 0.0;
      if( vVal <= c )
      {
        // reset is done, go back to old form
        resetting = false;
      }
#endif
//#if 0

      //std::cout << "In resetting section. (uPeak*uscale + d) = " << (uPeak + d/uscale) << std::endl;
      vEquFvalue = -fallRate*(vVal - c) - uVal*uscale;
      vEquQvalue = -memCap* vVal;
      vEqudFdv   = -fallRate;
      vEqudFdu   = -uscale;
      vEqudQdv   = -memCap;

      uEquFvalue = -fallRate*(uVal - (uPeak + d/uscale) );
      uEquQvalue = -uVal;
      uEqudFdv   =  0.0;
      uEqudFdu   = -fallRate;
      uEqudQdu   = -1.0;

//#endif
    }
    else
    {
      //std::cout << "In normal section." << std::endl;
      vEquFvalue = k * (vVal - Vr) * (vVal - Vt) - uVal*uscale;
      vEquQvalue = - memCap * vVal;
      vEqudFdv   = k * ( 2*vVal - Vt - Vr);
      vEqudFdu   = -uscale;
      vEqudQdv   = -memCap;

      uEquFvalue = a * ( b * (vVal - Vr)/uscale - uVal );
      uEquQvalue = -uVal;
      uEqudFdv   =  a * b / uscale;
      uEqudFdu   = -a;
      uEqudQdu   = -1.0;
    }


  }

  /*
    std::cout << "Instance::updateIntermediateVars()" << std::endl
    << "vEquFvalue = " <<  vEquFvalue << std::endl
    << "vEquQvalue = " << vEquQvalue << std::endl
    << "vEqudFdv   = " << vEqudFdv << std::endl
    << "vEqudFdu   = " << vEqudFdu << std::endl
    << "vEqudQdv   = " << vEqudQdv << std::endl
    << "uEquFvalue = " << uEquFvalue << std::endl
    << "uEquQvalue = " << uEquQvalue << std::endl
    << "uEqudFdv   = " << uEqudFdv << std::endl
    << "uEqudFdu   = " << uEqudFdu << std::endl
    << "uEqudQdu   = " << uEqudQdu << std::endl
    << "vVal       = " << vVal << std::endl
    << "uVal       = " << uVal << std::endl;
  */
  return bsuccess;
}

/*
  this approach can be problematic if many of these devices are out of
  phase and issueing breakpoints all the time.  Thus, not very workable -- RLS
//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints(
vector<N_UTL_BreakPoint> &breakPointTimes)
{
breakPointTimes.push_back(breakPoint);
return true;
}
*/

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

  updateIntermediateVars();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 7 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  N_LAS_Vector * daeQVecPtr = extData.daeQVectorPtr;
  (*daeQVecPtr)[li_V] += vEquQvalue;
  (*daeQVecPtr)[li_U] += uEquQvalue;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 7 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  N_LAS_Vector * daeFVecPtr = extData.daeFVectorPtr;

  (*daeFVecPtr)[li_V]  += vEquFvalue;
  (*daeFVecPtr)[li_U]  += uEquFvalue;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 7 instance.
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  N_LAS_Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  (*dQdxMatPtr)[li_V][vEquVOffset] += vEqudQdv;
  (*dQdxMatPtr)[li_U][uEquUOffset] += uEqudQdu;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 7 instance.
//
// Special Notes : This is an algebraic constaint.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[li_V][vEquVOffset] += vEqudFdv;
  (*dFdxMatPtr)[li_V][vEquUOffset] += vEqudFdu;

  (*dFdxMatPtr)[li_U][uEquVOffset] += uEqudFdv;
  (*dFdxMatPtr)[li_U][uEquUOffset] += uEqudFdu;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::varTypes( vector<char> & varTypeVec )
{
  //varTypeVec.resize(1);
  //varTypeVec[0] = 'I';
}


//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
              SolverState & ss1,
              DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
    memCap(0.0),
    Vt(0.0),
    Vr(0.0),
    Vpeak(0.0),
    k(0.0),
    a(0.0),
    b(0.0),
    c(0.0),
    d(0.0),
    uscale(1.0e-9),
    fallRate(1.0e3),
    memCapGiven(false),
    VtGiven(false),
    VrGiven(false),
    VpeakGiven(false),
    kGiven(false),
    aGiven(false),
    bGiven(false),
    cGiven(false),
    dGiven(false),
    uscaleGiven(false),
    fallRateGiven(false)
{

  /*
    std::cout << "Model::Model" << std::endl
    << "memC  = " << memCap << std::endl
    << "Vt    = " << Vt << std::endl
    << "Vr    = " << Vr << std::endl
    << "Vp    = " << Vpeak << std::endl
    << "k     = " << k << std::endl
    << "a     = " << a << std::endl
    << "b     = " << b << std::endl
    << "c     = " << c << std::endl
    << "d     = " << d << std::endl
    << "uscale = " << uscale << std::endl
    << "fallRate = " << fallRate << std::endl;
  */

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  //if (!given("TNOM"))
  //  tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  processParams ();
}


//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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

// additional Declarations
//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << endl;
  os << "Number of Neuron instances: " << isize << endl;
  os << "    name=\t\tmodelName\tParameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;
}

} // namespace Neuron2
} // namespace Device
} // namespace Xyce
