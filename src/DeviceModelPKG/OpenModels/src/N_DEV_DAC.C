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
// Filename       : $RCSfile: N_DEV_DAC.C,v $
//
// Purpose        : This file implements the DAC digital to analog conversion
//                  device used in the integration of Xyce withe SAVANT VHDL
//                  simulator.
//
// Special Notes  :
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

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#include <algorithm>

// ----------   Xyce Includes   ----------
#include <N_DEV_DAC.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<DAC::Instance>::ParametricData()
{
  // Set up configuration constants:
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(0);
  addModelType("DAC");

  // Set up double precision variables:

  // Set up exceptions (ie variables that are not doubles):

  addPar ("FILE", string(""), false, ParameterType::NO_DEP,
          &DAC::Instance::file, NULL);
}

template<>
ParametricData<DAC::Model>::ParametricData()
{
  // Set up double precision variables:
  addPar ("TR", 1.e-9, false,   ParameterType::NO_DEP,
          &DAC::Model::riseTime,
          NULL,U_SECOND,CAT_NONE,"Rise Time");

  addPar ("TF", 1.e-9, false,   ParameterType::NO_DEP,
          &DAC::Model::fallTime,
          NULL,U_SECOND,CAT_NONE,"Fall Time");

  addPar ("R",   0.01, false,   ParameterType::NO_DEP,
          &DAC::Model::R,
          NULL,U_OHM,CAT_NONE,"Resitance");

  addPar ("L",  1.e-5, false,   ParameterType::NO_DEP,
          &DAC::Model::L,
          NULL,U_HENRY,CAT_NONE,"Inductance");

  addPar ("C",    0.0, false,   ParameterType::NO_DEP,
          &DAC::Model::C,
          NULL,U_FARAD,CAT_NONE,"Capacitance");

  // Set up non-double precision variables:
  addPar ("TRANBP", true, false, ParameterType::NO_DEP,
          &DAC::Model::includeTransitionBP_, NULL,
          U_NONE, CAT_CONTROL, "Flag for transitional breakpoints");
}

namespace DAC {


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
  Model & DACiter,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)
  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(DACiter),
    v_pos(0),
    v_neg(0),
    i_bra(0),
    vDrop(0),
    voltage_(0),
    file(""),
    loc_(0),
    numTVpairs_(0),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1)
{

  setName(IB.getName());
  setModelName(model_.getName());

  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    jacStamp[0].resize(1);
    jacStamp[0][0] = 2;
    jacStamp[1].resize(1);
    jacStamp[1][0] = 2;
    jacStamp[2].resize(2);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

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
    cout << "  DACInstance::registerLIDs" << endl;
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

  li_Bra = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
    cout << "  li_Bra = " << li_Bra << endl;
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
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSi
// Creation Date : 11/20/09
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    string tmpstr;
    tmpstr = getName()+"_branch";
    spiceInternalName (tmpstr);
    intNameMap[li_Bra] = tmpstr;
  }

  return intNameMap;
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
  APosEquBraVarOffset = jacLIDVec[0][0];
  ANegEquBraVarOffset = jacLIDVec[1][0];
  ABraEquPosNodeOffset = jacLIDVec[2][0];
  ABraEquNegNodeOffset = jacLIDVec[2][1];
}

//----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;
  double * solVector = extData.nextSolVectorRawPtr;

  // Get the value for the source.
  updateVoltage(getSolverState().acceptedTime);

  // get the value for v_pos, v_neg, i_bra
  v_pos = solVector[li_Pos];
  v_neg = solVector[li_Neg];
  i_bra  = solVector[li_Bra];
  vDrop = (v_pos-v_neg-voltage_);

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
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;
  bsuccess = updateIntermediateVars ();
  return bsuccess;
}

//----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Instance::updateSecondaryState()
{
  bool bsuccess = true;

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : Instance::updateTVVEC
// Purpose        : Append the contents of newPairs to TVVEC
// Special Notes  :
// Scope          : public
// Creator        : Lisa Maynes & Lon Waters
// Creation Date  : 06/10/2003
//----------------------------------------------------------------------------
bool Instance::updateTVVEC (
  vector< pair<double, double> > const & newPairsIn )
{
  int i, last, newStart;
  double transitionTime;
  bool bsuccess = true;
  vector < pair<double,double> >::iterator itVec, itVec_end;
  vector< pair<double, double> > newPairs;
  vector < pair<double,double> >::const_iterator const_itVec;
  map <double, double> tv;
  map <double, double>::iterator tv_i, tv_end;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "In device " << getName() << endl;
    cout << "At time = " <<  getSolverState().acceptedTime <<
      ", at beginning of Instance::updateTVVEC():\n"
         << "   TVVEC size = " << numTVpairs_ << "\n"
         << "   TVVEC loc  = " << loc_ << "\n"
         << "   TVVEC contents:\n" ;
    itVec = TVVEC.begin();
    itVec_end = TVVEC.end();
    for( ; itVec != itVec_end; ++itVec )
    {
      cout << "   " << (*itVec).first
           << "s, " << (*itVec).second
           << "V\n";
    }
    cout << newPairsIn.size() << " New pairs:" << endl;
    for ( const_itVec = newPairsIn.begin() ; const_itVec != newPairsIn.end() ; ++const_itVec )
      cout << (*const_itVec).first << "  " << (*const_itVec).second << endl;
  }
#endif

  updateVoltage(getSolverState().acceptedTime);

  itVec = TVVEC.begin();
  itVec_end = TVVEC.end();
  for( ; itVec != itVec_end; ++itVec )
    tv[(*itVec).first] = (*itVec).second;

  if (!newPairsIn.empty())
  {
    if (getSolverState().acceptedTime == 0)
    {
      if (TVVEC.size() > 0 && TVVEC[0].first == 0) {
        TVVEC.resize(1);
      }
      else
      {
        TVVEC.push_back(pair<double,double>(0,(*(newPairsIn.end()-1)).second));
      }
      if ((*(newPairsIn.end()-1)).first > TVVEC[0].first)
        TVVEC.push_back(*(newPairsIn.end()-1));
      numTVpairs_ = TVVEC.size();
      updateVoltage(getSolverState().acceptedTime);
      return bsuccess;
    }
    vector< pair<double, double> >::const_iterator n_i, n_end;
    vector< pair<double, double> >::iterator t_i, t_end;
    n_i = newPairsIn.begin();
    n_end = newPairsIn.end();
    for ( ; n_i != n_end ; ++n_i)
    {
      if ((*n_i).first < 0)
      {
        double d = -(*n_i).first;
        tv_i = lower_bound(tv.begin(), tv.end(), pair<const double, double>(d,0));
        tv.erase(tv_i,tv.end());
      }
    }
    n_i = newPairsIn.begin();
    for ( ; n_i != n_end ; ++n_i)
    {
      if ((*n_i).first >= getSolverState().acceptedTime)
      {
        transitionTime = model_.riseTime;
        if (transitionTime > 0) {
          tv_i = lower_bound(tv.begin(), tv.end(), pair<const double, double>((*n_i).first,0));
          if (tv_i != tv.begin())
          {
            --tv_i;
            if ((*n_i).second < (*tv_i).second)
              transitionTime = model_.fallTime;
          }
          tv[(*n_i).first] = (*tv_i).second;
        }
        tv[(*n_i).first + transitionTime] = (*n_i).second;
      }
    }
    tv[getSolverState().acceptedTime] = voltage_;
  }
  tv_i = lower_bound(tv.begin(), tv.end(), pair<const double, double>(getSolverState().acceptedTime,0));
  if (tv_i != tv.begin())
    --tv_i;
  tv.erase(tv.begin(), tv_i);
  double lastTimeEntry = -1;
  TVVEC.clear();
  tv_i = tv.begin();
  tv_end = tv.end();
  for ( ; tv_i != tv_end ; ++tv_i)
  {
    if ((*tv_i).first - lastTimeEntry > 1e-15)
    {
      TVVEC.push_back(pair<double, double>((*tv_i).first,(*tv_i).second));
      lastTimeEntry = (*tv_i).first;
    }
  }
  numTVpairs_ = TVVEC.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::updateTVVEC():\n"
         << "   TVVEC size = " << numTVpairs_ << "\n"
         << "   TVVEC loc  = " << loc_ << "\n"
         << "   TVVEC contents:\n" ;
    vector< pair<double, double> >::iterator tv_i = TVVEC.begin();
    for( ; tv_i != TVVEC.end(); ++tv_i )
    {
      cout << "   " << (*tv_i).first
           << "s, " << (*tv_i).second
           << "V\n";
    }
  }
#endif

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : Instance::updateVoltage
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 11/12/2002
//----------------------------------------------------------------------------
bool Instance::updateVoltage(double time)
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << endl;
    cout << "  DACInstance::updateVoltage\n";
    cout << "-------------------------------" << endl;
    cout << "  Time = " << time << endl;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag &&  numTVpairs_ > 0)
  {
    cout << "    TVVEC[numTVpairs_-1].first = "
         << TVVEC[numTVpairs_-1].first << endl;
  }
  for (int i=0 ; i<numTVpairs_ ; ++i) {
    cout << TVVEC[i].first << " :: " << TVVEC[i].second << endl;
  }
#endif

  if( numTVpairs_ > 0 && time >= TVVEC[0].first )
  {
    if( time < TVVEC[numTVpairs_-1].first )
    {
      for( int i = 0; i < numTVpairs_ - 1; ++i )
      {
        if( time >= TVVEC[i].first && time <= TVVEC[i+1].first)
        {
          double time1 = TVVEC[i].first;
          double voltage1 = TVVEC[i].second;

          double time2 = TVVEC[i+1].first;
          double voltage2 = TVVEC[i+1].second;

          voltage_ = voltage1 + (voltage2 - voltage1) * (time - time1) / (time2 - time1);
          break;
        }
      }
    }
    else
    {
      voltage_ = TVVEC[numTVpairs_-1].second;
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "  voltage_ = " << voltage_ << endl;
    cout << "-------------------------------" << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for this device.
//
// Special Notes : This is an algebraic constaint, and as such the resistor
//                 does make a contribution to it.
//
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;

  double * daeFVec = extData.daeFVectorRawPtr;


  daeFVec[li_Pos] += i_bra;

  daeFVec[li_Neg] += -i_bra;

  daeFVec[li_Bra] += vDrop;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for this device.
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


  (*dFdxMatPtr)[li_Pos][APosEquBraVarOffset] += 1.0;

  (*dFdxMatPtr)[li_Neg][ANegEquBraVarOffset] -= 1.0;

  (*dFdxMatPtr)[li_Bra][ABraEquPosNodeOffset] += 1.0;

  (*dFdxMatPtr)[li_Bra][ABraEquNegNodeOffset] -= 1.0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 2/16/04
//-----------------------------------------------------------------------------

bool Instance::getInstanceBreakPoints ( vector<N_UTL_BreakPoint> & breakPointTimes)
{
  bool bsuccess = true;
  double currentTime = getSolverState().currTime;
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "In ::getInstanceBreakPoints " << endl;
    cout << " I want breakpoints.  Current time is " << currentTime << endl;
  }
#endif

  for (int i = 0; i < numTVpairs_ ; ++i)
  {
    // ERK:  the 1e-15 is a tolerance.  Fix for bug 1766.  Possibly use bpTol instead?
    // DNS:  No, bpTol tends to be ridiculously small and this might cause
    //       excessively small time steps.  I would agree if bpTol had a floor value,
    //       such as 1e-15.  Steps smaller than this are unjustified and numerically
    //       problematic.  For reference, light travels 0.3 micron in 1e-15 seconds.
    if (TVVEC[i].first >= currentTime - 1e-15 && model_.riseTime != 0 && model_.fallTime != 0)
    {
      breakPointTimes.push_back(N_UTL_BreakPoint(TVVEC[i].first,
                                                 SIMPLE_BREAKPOINT));
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "DAC ----------------------------------------" << endl;
    cout << "DAC getInstanceBreakPoints " << endl;
    cout << "DAC Debug output.  name = " << getName() << endl;
    cout << "DAC setting breakpoints at currentTime = " << currentTime << endl;
    cout << "DAC breakpoints: " << endl;

    vector< N_UTL_BreakPoint  >::iterator beg = breakPointTimes.begin();
    vector< N_UTL_BreakPoint  >::iterator end = breakPointTimes.end();
    vector< N_UTL_BreakPoint  >::iterator itBP = beg;
    for (;itBP!=end;++itBP)
    {
      N_UTL_BreakPoint & bp = *itBP;
      cout << "DAC breakpoint: " << bp.value() << endl;
    }

    cout << "DAC ----------------------------------------" << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInternalState
// Purpose       : Generates an DeviceState object and populates
//                 it with the contents of the TVVEC vector for use by
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
  int vsize,i,j;
  // allocate object to return
  DeviceState * myState = new DeviceState;

  myState->ID=getName();
  vsize=TVVEC.size();
  // pack the pairs into the single vector of doubles.
  myState->data.resize(vsize*2);
  for (i=0;i<vsize;++i)
  {
    j=i*2;
    myState->data[j]=TVVEC[i].first;
    myState->data[j+1]=TVVEC[i].second;
  }

  return myState;
}
//-----------------------------------------------------------------------------
// Function      : Instance::setInternalState
// Purpose       : Reload TVVEC data from restart
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
  int vsize,i,j;
  if ( state.ID != getName())
  {
    string msg;
    msg = "Instance::setInternalState:  ID ("+state.ID+")";
    msg += "from restart does not match my name ("+getName()+")!\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  if (dsize%2 != 0)
  {
    string msg;
    char msg2[256];
    sprintf(msg2, "Instance::setInternalState: "
            "Data size from restart (%d) not a multiple of 2!",
            dsize);
    msg=msg2;
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  vsize=dsize/2;
  TVVEC.clear();
  TVVEC.resize(vsize);
  for (i=0;i<vsize;++i)
  {
    j=i*2;
    TVVEC[i].first=state.data[j];
    TVVEC[i].second=state.data[j+1];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 05/12/09
//-----------------------------------------------------------------------------
void Instance::varTypes( vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
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
    riseTime(1.0e-9),
    fallTime(1.0e-9),
    R(.01),
    C(0.0),
    L(0.0),
    includeTransitionBP_(true)
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
    os << "\t\tfile = " << (*iter)->file;
    os << endl;
  }

  os << endl;

  return os;
}


// DAC Master functions:

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

    // here we update the voltage
    di.updateVoltage(getSolverState().currTime);

    // get the value for v_pos, v_neg, i_bra
    di.v_pos = solVec[di.li_Pos];
    di.v_neg = solVec[di.li_Neg];
    di.i_bra  = solVec[di.li_Bra];
    di.vDrop = (di.v_pos-di.v_neg-di.voltage_);
#ifdef Xyce_DEBUG_DEVICE
    cout << "DAC ----------------------------------------" << endl;
    cout << "DAC Debug output.  name = " << di.getName() << endl;
    cout << "DAC v_pos = " << di.v_pos <<endl;
    cout << "DAC v_neg = " << di.v_neg <<endl;
    cout << "DAC i_bra = " << di.i_bra <<endl;
    cout << "DAC vDrop = " << di.vDrop <<endl;
    cout << "DAC voltage_ = " << di.voltage_ <<endl;
    cout << "DAC ----------------------------------------" << endl;
#endif
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
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);

    fVec[di.li_Pos] += di.i_bra;

    fVec[di.li_Neg] += -di.i_bra;

    fVec[di.li_Bra] += di.vDrop;
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

    dFdx[di.li_Pos][di.APosEquBraVarOffset] += 1.0;

    dFdx[di.li_Neg][di.ANegEquBraVarOffset] -= 1.0;

    dFdx[di.li_Bra][di.ABraEquPosNodeOffset] += 1.0;

    dFdx[di.li_Bra][di.ABraEquNegNodeOffset] -= 1.0;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::getDACDeviceNames
// Purpose       : Get a list of names of all DACs in the circuit
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and MOdels
// Creation Date : 05/05/04
//-----------------------------------------------------------------------------
bool Master::getDACDeviceNames( vector<string> & dacNames)
{
  bool bSuccess = true;

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & di = *(*it);
    dacNames.push_back(di.getName());
  }
  return bSuccess;
}

} // namespace DAC
} // namespace Device
} // namespace Xyce
