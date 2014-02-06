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
// Filename       : $RCSfile: N_DEV_ADC.C,v $
//
// Purpose        : This file implements the ADC digital to analog conversion
//                  device used in the integration of Xyce withe SAVANT VHDL
//                  simulator.
//
// Special Notes  : The ADC looks to the analog simulation like a resistor.
//
// Creator        : Lon Waters
//
// Creation Date  : 07/26/2002
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revsion$
//
// Revsion Date   : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#include <N_UTL_Misc.h>

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <algorithm>

// ----------   Xyce Includes   ----------
#include <N_DEV_ADC.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<ADC::Instance>::ParametricData()
{
  // Set up configuration constants:
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(0);
  addModelType("ADC");

  // Set up double precision variables:
  addPar ("R", 1.e+12, false, ParameterType::NO_DEP, &ADC::Instance::R,
          NULL, U_OHM, CAT_NONE, "internal Resistance");
}

template<>
ParametricData<ADC::Model>::ParametricData()
{
  addPar ("LOWERVOLTAGELIMIT", 0.0, false, ParameterType::NO_DEP,
          &ADC::Model::lowerVoltageLimit_,
          NULL,U_VOLT,CAT_NONE,"Lower limit of ADC voltage range");

  addPar ("UPPERVOLTAGELIMIT", 5.0, false, ParameterType::NO_DEP,
          &ADC::Model::upperVoltageLimit_,
          NULL,U_VOLT,CAT_NONE,"Upper limit of ADC voltage range");

  addPar ("SETTLINGTIME", 1.0e-8, false, ParameterType::NO_DEP,
          &ADC::Model::settlingTime_,
          NULL,U_SECOND,CAT_NONE,"Settling time");
}

namespace ADC {

vector< vector<int> > Instance::jacStamp;


ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

//----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Instance::processParams (string param)
{

  return true;
}

//----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
Instance::Instance(
  InstanceBlock & IB,
  Model & ADCiter,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)
  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    R(1.0e12),
    G(0.0),
    v_pos(0.0),
    v_neg(0.0),
    i0(0.0),
    outputBitVectorWidth_(1),
    nQuantLevels_(0),
    lastOutputLevel_(0),
    model_(ADCiter),
    li_Pos(-1),
    li_Neg(-1),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1)
{

  setName(IB.getName());
  setModelName(model_.getName());

  numIntVars   = 0;
  numExtVars   = 2;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    jacStamp.resize(2);
    jacStamp[0].resize(2);
    jacStamp[1].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  if (R != 0.0) G = 1.0/R;
  else          G = 0.0;

  processParams ();
}

//----------------------------------------------------------------------------
// Function       : Instance::~Instance
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
Instance::~Instance()
{
}

//----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                             const vector<int> & extLIDVecRef)
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "-----------------------------------------------------------------------------";
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << endl << dashedline << endl;
    cout << "  ADCInstance::registerLIDs" << endl;
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

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
    cout << "  li_Pos = " << li_Pos << endl;
#endif

  li_Neg = extLIDVec[1];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
    cout << "  li_Neg = " << li_Neg << endl;
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
    cout << dashedline << endl;
#endif
}

//----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef)
{
  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    string msg("Instance::registerStateLIDs:");
    msg += "numSTa != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
}

//----------------------------------------------------------------------------
// Function       : jacobianStamp
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 11/12/2002
//----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//----------------------------------------------------------------------------
// Function       : registerJacLIDs
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 11/12/2002
//----------------------------------------------------------------------------
void Instance::registerJacLIDs(const vector< vector<int> >& jacLIDVec)
{
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
}

//----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;
  double * solVector = extData.nextSolVectorRawPtr;

  v_pos = solVector[li_Pos];
  v_neg = solVector[li_Neg];
  i0 = (v_pos-v_neg)*G;

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  return updateIntermediateVars ();
}

//----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  // There is no secondary state
  return true;
}

//----------------------------------------------------------------------------
// Function       : Instance::getTVVec
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/28/2003
//----------------------------------------------------------------------------
void Instance::getTVVEC(vector< pair<double, double> > & TVVEC_Out)
{
  TVVEC_Out.clear();

#ifdef HAVE_FLEXIBLE_INSERT
  TVVEC_Out.insert(TVVEC_Out.end(), TVVEC.begin(), TVVEC.end());
#else
  vector< pair<double, double> >::iterator iterTV = TVVEC.begin();
  vector< pair<double, double> >::iterator iterEnd = TVVEC.end();
  for(; iterTV != iterEnd; ++iterTV)
    TVVEC_Out.push_back( *iterTV );
#endif

  // Now that the digital simulator knows about these, let's forget them
  // and not worry about maintaining the list anymore
  TVVEC.clear();
}


//----------------------------------------------------------------------------
// Function      : Instance::trimTVVEC
// Purpose       : clear out old time-voltage pairs
// Special Notes : ASSUMES the vector of t-v pairs is sorted by time!
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/10/2004
//----------------------------------------------------------------------------
void Instance::trimTVVEC(double earliestTime)
{
  vector< pair<double,double> >::iterator itVec;

  // get reference pointing to first element that exceeds earliestTime
  itVec = lower_bound(TVVEC.begin(),TVVEC.end(),pair<double,double>(earliestTime,0.0));
  // delete everything prior to that.
  TVVEC.erase(TVVEC.begin(),itVec);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for this device which
//                 effectively is a single instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;

  double * daeFVec = extData.daeFVectorRawPtr;

  // Load DAE F-vector

  daeFVec[li_Pos] += i0;

  daeFVec[li_Neg] += -i0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for this device which
//                 effectively is a single instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline(
    "-----------------------------------------------------------------------------");
  const string dashedline2("---------------------");

  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 <<endl;
    cout << "  name = " << getName() << endl;
    cout << "  G = " << G << endl;
  }
#endif


  (*dFdxMatPtr)[li_Pos][APosEquPosNodeOffset] += G;

  (*dFdxMatPtr)[li_Pos][APosEquNegNodeOffset] -= G;

  (*dFdxMatPtr)[li_Neg][ANegEquPosNodeOffset] -= G;

  (*dFdxMatPtr)[li_Neg][ANegEquNegNodeOffset] += G;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/04
//-----------------------------------------------------------------------------

bool Instance::getInstanceBreakPoints ( vector<N_UTL_BreakPoint> & breakPointTimes)
{
  bool bsuccess = true;

  // getBreakPoints is called after a timesteps solution is converged and
  // the step is finished.  So we can count on our voltage being the real one
  // that is valid for this time.

  // In this routine we need to check our current voltage value and see if it
  // has changed enough since our last time step to change a bit in the digital
  // output.  If so, we need to set a time/voltage pair of the current value
  // and (the current time + this device's conversion time), then set a
  // breakpoint of type PAUSE_BREAKPOINT for the that time

  double vPos(0.0);
  double vNeg(0.0);
  double deltaV(0.0), vFrac(0.0);
  double currentTime = getSolverState().currTime;
  int newState;

  // Get the pointer to the vector of accepted solution values
  double * solVector = extData.nextSolVectorRawPtr;

  vPos = solVector[li_Pos];
  vNeg = solVector[li_Neg];

  deltaV = vPos-vNeg;

  // This "upper voltage limit" and "lower voltage limit" crap is not right,
  // and needs to be replaced with a "Vref+"  node against which
  // Vin is compared, with a common negative reference (e.g. ground)
  // For now, let's always just document this failing, and tell the users
  // to wire vNeg to ground, and use 0.0 as the
  // lower limit.  Worry about doing it right later.
  vFrac = deltaV/(model_.upperVoltageLimit_
                  - model_.lowerVoltageLimit_);

#if 1
  // I believe this is how it SHOULD be done
  if (vFrac < (1.0)/(nQuantLevels_) )
  {
    newState = 0;
  }
  else if (vFrac >= (nQuantLevels_-1.0)/(nQuantLevels_))
  {
    newState = nQuantLevels_ -1;
  }
  else
  {
    newState =  int(vFrac*nQuantLevels_);
  }
#else
  // but this is how simbus is doing it:
  newState = int(vFrac*(nQuantLevels_-1));
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << "In Instance::getInstanceBreakPoints.  deltaV = " << deltaV
              << " and output state is " << newState
              << " vFrac = " << vFrac << " nQuantLevels = " << nQuantLevels_ << std::endl;
  }
#endif

  if (newState != lastOutputLevel_)
  {
    // get time rounded to nearest femptosecond
    long long int timeInFS =static_cast<long long int>(
      (currentTime+model_.settlingTime_+6e-16)/1e-15);

    // we need to pause
#ifdef Xyce_OLD_PRE_ROLLBACK
    breakPointTimes.push_back( N_UTL_BreakPoint(timeInFS*1e-15,PAUSE_BREAKPOINT));
#else

    // See bug # 1720 on charleston bugzilla.  This breakpoint pushback is being
    // commented out to fix this bug.
    //
    // However, this should not be a long-term fix!

#endif

    // since time is always advancing, we can simply push_back the pair,
    // and always be sure
    TVVEC.push_back(pair<double,double>(timeInFS*1e-15,deltaV));
    lastOutputLevel_ = newState;
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "ADC ----------------------------------------" << endl;
      cout << "ADC Debug output.  name = " << getName() << endl;
      cout << "ADC setting pause breakpoint for " << currentTime << endl;
      double approxTime = timeInFS*1e-15;
      cout << "ADC Approximated time, to nearest femptosecond: " << approxTime << endl;
      cout << "ADC Time value pairs: " << endl;

      vector< pair<double, double> >::iterator beg = TVVEC.begin();
      vector< pair<double, double> >::iterator end = TVVEC.end();
      vector< pair<double, double> >::iterator itTV = beg;

      for (; itTV!=end; ++itTV)
      {
        pair<double, double>  & tvPair = *itTV;
        cout << "ADC time: " << tvPair.first << "  value: " << tvPair.second << endl;
      }

      cout << "ADC ----------------------------------------" << endl;
    }
#endif
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/05/08
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  double vPos(0.0), vNeg(0.0), deltaV(0.0), vFrac(0.0);
  double currentTime = getSolverState().currTime;

  if (getSolverState().dcopFlag)
  {
    double * solVector = extData.nextSolVectorRawPtr;

    vPos = solVector[li_Pos];
    vNeg = solVector[li_Neg];
    deltaV = vPos-vNeg;

    TVVEC.push_back(pair<double,double>(0.0,deltaV));
  }

}

//----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Model::processParams(string param)
{
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

//----------------------------------------------------------------------------
// Function       : Instance::getInstanceParamsMap
// Purpose        : Return a map of name,parameter for this instance
// Special Notes  : used by API for mixed signal
//                  in zero order version, these parameters happen to be
//                  model parameters, but that is not necessarily where they'll
//                  stay (the voltage limits, for example, should actually be
//                  taken from the nodal voltages of the reference nodes of the
//                  device instance)
// Scope          : public
// Creator        : Tom Russo, SNL, Component Information and Models
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Instance::getInstanceParamsMap(map<string,double>& paramsMap)
{
  paramsMap.clear();

  paramsMap["lowerVoltageLimit"] = model_.lowerVoltageLimit_;
  paramsMap["upperVoltageLimit"] = model_.upperVoltageLimit_;
  paramsMap["settlingTime"] = model_.settlingTime_;

  return true;
}

//----------------------------------------------------------------------------
// Function       : Instance::setBitVectorWidth
// Purpose        : set the number of bits in the output of this ADC
// Special Notes  :
// Scope          : public
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Instance::setBitVectorWidth(int width)
{
  outputBitVectorWidth_ = width;
  nQuantLevels_ = 1;
  for (int i=0; i<width;++i)
    nQuantLevels_ *= 2;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << "Instance::setBitVectorWidth name = "
              << getName() << " width = " << width
              << " nQuantLevels_ = " << nQuantLevels_ << std::endl;
  }
#endif

  return true;
}
//----------------------------------------------------------------------------
// Function       : Model::Model
// Purpose        : constructor
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
Model::Model(const ModelBlock& MB,
             SolverState& ss1,
             DeviceOptions& do1)
  : DeviceModel(MB, ss1, do1),
    lowerVoltageLimit_(0.0),
    upperVoltageLimit_(5.0),
    settlingTime_(1e-8)

{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//----------------------------------------------------------------------------
// Function       : Model::Model
// Purpose        : destructor
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
Model::~Model()
{
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

//----------------------------------------------------------------------------
// Function       : printOutInstances
// Purpose        : debugging tool
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << endl;
  os << "    name\t\tmodelName\tParameters" << endl;

  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;

  return os;
}

// ADC Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);

    di.v_pos = solVec[di.li_Pos];
    di.v_neg = solVec[di.li_Neg];
    di.i0 = (di.v_pos-di.v_neg)*di.G;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);

    fVec[di.li_Pos] += di.i0;

    fVec[di.li_Neg] += -di.i0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);

    dFdx[di.li_Pos][di.APosEquPosNodeOffset] += di.G;

    dFdx[di.li_Pos][di.APosEquNegNodeOffset] -= di.G;

    dFdx[di.li_Neg][di.ANegEquPosNodeOffset] -= di.G;

    dFdx[di.li_Neg][di.ANegEquNegNodeOffset] += di.G;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ::getADCMap
// Purpose       : Get a list of names of all ADCs in the circuit
// Special Notes : We assume the map we're given is cleared already
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::getADCMap( map<string, map<string,double> > & ADCMap )
{
  bool bSuccess = true;
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);
    std::string adcName(di.getName());

    map<string,double> TmpMap;
    di.getInstanceParamsMap(TmpMap);
    ADCMap[adcName] = TmpMap;
  }
  return bSuccess;
}

//-----------------------------------------------------------------------------
// Function      : ::getTimeVoltagePairs
// Purpose       : Get time/voltage pairs for all ADCs in the circuit
// Special Notes : We assume the map we're given is cleared already
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::getTimeVoltagePairs( map<string, vector< pair<double,double> > >& TimeVoltageMap )
{
  bool bSuccess = true;
  double vPos, vNeg;

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);
    std::string adcName(di.getName());

    vector<pair <double,double> > TmpVec;
    di.getTVVEC(TmpVec);
    double currentTime = getSolverState().currTime;
    double * solVector = di.extData.nextSolVectorRawPtr;

    vPos = solVector[di.li_Pos];
    vNeg = solVector[di.li_Neg];
    TmpVec.push_back(pair<double,double>(currentTime, vPos-vNeg));
    TimeVoltageMap[adcName] = TmpVec;
  }
  return bSuccess;
}

//-----------------------------------------------------------------------------
// Function      : ::getBreakPoints
// Purpose       : getBreakPoints for all instances
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::getBreakPoints ( vector<N_UTL_BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);
    bool tmpBool = di.getInstanceBreakPoints(breakPointTimes);
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;

}

} // namespace ADC
} // namespace Device
} // namespace Xyce
