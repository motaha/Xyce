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
// Filename       : $RCSfile: N_DEV_ExternDevice.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/09/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.84.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:35 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Const.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_DEV_XyceInterface.h>
#include <N_DEV_CharonInterface.h>
#include <N_DEV_1D_Interface.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_PDS_Comm.h>

#include <N_TIA_TimeIntInfo.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<VoltageNode>::ParametricData()
{
  // Set up map for double precision variables:
  addPar("INITVAL", 0.0, false, &VoltageNode::initVal);
  addPar("LIMITHIGH", 0.0, false, &VoltageNode::limValHigh);
  addPar("LIMITLOW", 0.0, false, &VoltageNode::limValLow);

  // Set up exceptions (ie variables that are not doubles):
  addPar("NAME", "VCONNECT0000", false, &VoltageNode::vsrcName);
}

ParametricData<VoltageNode> &VoltageNode::getParametricData() {
  static ParametricData<VoltageNode> parMap;

  return parMap;
}


//-----------------------------------------------------------------------------
// Function      : VoltageNode::VoltageNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/06
//-----------------------------------------------------------------------------
VoltageNode::VoltageNode ()
  : CompositeParam (),
    vsrcName(""),
    initVal(0.0),
    limValHigh(0.15),
    limValLow(0.10)
{

  // Set up mapping from param names to class variables:

  // if (parMap.empty())
  // {
  //   // Set up map for double precision variables:
  //   addPar ("INITVAL", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&VoltageNode::initVal),
  //           NULL);

  //   addPar ("LIMITHIGH", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&VoltageNode::limValHigh),
  //           NULL);

  //   addPar ("LIMITLOW", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&VoltageNode::limValLow),
  //           NULL);

  //   // Set up exceptions (ie variables that are not doubles):
  //   addPar ("NAME", "VCONNECT0000", false,
  //           static_cast <string CompositeParam:: *> (&VoltageNode::vsrcName), NULL);
  // }
}

//-----------------------------------------------------------------------------
// Function      : VoltageNode::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/06
//-----------------------------------------------------------------------------
void VoltageNode::processParams ()
{

}

template<>
ParametricData<ExternDevice::Instance>::ParametricData()
{
  // Set up configuration constants:
  setNumNodes(2);
  setNumOptionalNodes(1000);
  setNumFillNodes(0);
  setModelRequired(0);

  // Set up map for double precision variables:

  // Set up map for non-double precision variables:
  addPar("EXTERNCODE", string("xyce"), false, ParameterType::NO_DEP, &ExternDevice::Instance::externCode_);
  addPar("NETLIST", string("input.cir"), false, ParameterType::NO_DEP, &ExternDevice::Instance::netlistFileName_);
  addPar("VOLTLIM", false, false, ParameterType::NO_DEP, &ExternDevice::Instance:: voltageLimiterFlag);

  addComposite("NODE", VoltageNode::getParametricData(), &ExternDevice::Instance::nodeMap);
}

template<>
ParametricData<ExternDevice::Model>::ParametricData()
{}


namespace ExternDevice {



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
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::processParams(string param)
{

  string msg;

  if ( externCode_ == "xyce" )
  {
    if (extCodePtr_==0)
    {
      extCodePtr_ = new XyceInterface ( devOptions, getDeviceOptions().commandLine, netlistFileName_ );
    }
  }
  else if ( externCode_ == "1d")
  {
    msg = "Instance::processParams.  1D interface not suppported.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }
  else if ( externCode_ == "charon")
  {
    if (extCodePtr_==0)
    {
      extCodePtr_ = new CharonInterface ( devOptions, netlistFileName_,
						solState);
    }
  }
  else
  {
    msg = "Instance::processParams. " + externCode_;
    msg += " is not a recognized external device\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------

Instance::Instance(
    InstanceBlock & instance_block,
    Model &               model,
    MatrixLoadData &      matrix_load_data,
    SolverState &         solver_state,
    ExternData  &         extern_data,
    DeviceOptions &       device_options)
  : DeviceInstance(instance_block, matrix_load_data, solver_state, extern_data, device_options),
    model_(model),
    externCode_("xyce"),
    netlistFileName_("input.cir"),
    extCodePtr_(0),
    initializeFlag_(false),
    innerSolveStatus_(false),
    locallyConnected_(true),
    owningProc_(-1),
    comm_(0),
    newtonIterOld(0),
    voltageLimiterFlag(false),
    initJctGiven_ (false),
    nodeGiven_(false)
{

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << dashedline << endl;
    cout << "In Instance::Instance() \n\n";
    cout << "name             = " << getName() << endl;
  }
#endif

  // allow for a variable numbers of nodes
  numExtVars = instance_block.numExtVars;
  numStateVars = instance_block.numExtVars;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "numExtVars = " << numExtVars << endl;

    vector<Param>::const_iterator iter_t;
    vector<Param>::const_iterator begin_t = instance_block.params.begin();
    vector<Param>::const_iterator end_t   = instance_block.params.end();

    for (iter_t = begin_t; iter_t !=  end_t;  ++iter_t)
    {
      cout << "Tag: " << iter_t->tag();
      cout << "  Value = " << iter_t->sVal() << "  Given: " << iter_t->given() << endl;
    }
  }
#endif

  int size = numExtVars;
  if( jacStamp.empty() )
  {
    jacStamp.resize(size);
    for( int i=0; i<size ; i++ )
    {
      jacStamp[i].resize(size);
      for( int j=0; j<size; j++ )
      {
        jacStamp[i][j] = j;
      }
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    for( int i=0; i<size; i++)
    {
      for( int j=0; j<size; j++)
      {
        cout << "jacStamp[ " << i << ", " << j << "] = "
             << jacStamp[i][j] << endl;
      }
    }
    cout << endl;
  }
#endif

  currentOutputVector_.resize(size,0.0);
  conductanceJacobian_.resize(size);
  for (int i=0;i<size;++i)
  {
    conductanceJacobian_[i].resize(size,0.0);

    // This name MUST be all caps!
    char tmp[16];
    sprintf(tmp,"VCONNECT%04d",i);
    string vname(tmp);
    voltageInputMap_[vname] = 0.0;
  }

  // voltlim vectors:
  voltageOld_.resize      (size,0.0);
  voltageOrig_.resize     (size,0.0);
  voltageLastCall_.resize (size,0.0);
  voltageDiff_.resize     (size,0.0);
  voltageFactor_.resize   (size,0.0);
  voltageStateID_.resize  (size,0);


  // Set params to constant default values:

  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:

  setParams (instance_block.params);

  // Set any non-constant parameter defaults:

  if (!given("EXTERNCODE"))
    externCode_ = "xyce";
  if (!given("NETLIST"))
    netlistFileName_ = "input.cir";

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
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
Instance::~Instance()
{
  if (extCodePtr_!=0)
  {
    delete extCodePtr_ ;
    extCodePtr_ = 0;
  }

  // delete the pointers that are in the voltage node map.
  vector<VoltageNode *>::iterator iterNode = voltageNodeVec.begin();
  vector<VoltageNode *>::iterator endNode  = voltageNodeVec.end ();
  for (;iterNode!=endNode;++iterNode)
  {
    delete (*iterNode);
    (*iterNode) = 0;
  }
}

// Additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       : function for registering, and setting up, local ID's.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter , SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
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

  if ( numInt != numIntVars )
  {
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "\nEnd of Instance::register LIDs\n";
    cout << dashedline << endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter , SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
  string msg;

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
  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  Number of State LIDs: " << numSta << endl;
  }
#endif


  if( numSta != numStateVars )
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  int i=0;
  for ( ; i<numStateVars;++i)
  {
    voltageStateID_[i] = staLIDVecRef[i];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter , SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
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
// Creator       : Eric R. Keiter , SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  jacLIDs = jacLIDVec;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 external device instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 external device instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;
  N_LAS_Vector * daeFVecPtr = extData.daeFVectorPtr;
  N_LAS_Vector * dFdxdVpPtr = extData.dFdxdVpVectorPtr;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "In ::loadDAEFVector " << getName() << endl;
  }
#endif

  if( locallyConnected_ )
  {
    // sum the current vector into the RHS:
    for (int i=0;i<currentOutputVector_.size();++i)
    {
      int iRow = extLIDVec[i];

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        cout <<"Instance: iRow = " << iRow << "  current = ";
        cout << currentOutputVector_[i] <<endl;
      }
#endif
      (*daeFVecPtr)[iRow] -= currentOutputVector_[i];

      if (voltageLimiterFlag && getDeviceOptions().voltageLimiterFlag && !origFlag)
      {
        (*dFdxdVpPtr)[iRow] += voltageFactor_[i];
      }
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Done with Instance::loadDAEFVector "<< endl;
  }
#endif


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 external device instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 external device instance.
//
// Special Notes : This is an algebraic constaint, and as such the external device
//                 does make a contribution to it.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;
  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  if( locallyConnected_ )
  {
    // Load in the conductance Jacobian:
    int iRowLoc,iColLoc;
    int size = conductanceJacobian_.size();
    for (iRowLoc=0;iRowLoc < size; ++iRowLoc)
    {
      for (iColLoc=0;iColLoc < size; ++iColLoc)
      {
        int iRow = extLIDVec[iRowLoc];
        int iCol = jacLIDs[iRowLoc][iColLoc];
        double val = conductanceJacobian_[iRowLoc][iColLoc];
        (*dFdxMatPtr)[iRow][iCol] += val;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       : update primary state for one external device
//                 instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess=true;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Instance::updatePrimaryState" << endl;
    cout << "  name   = " << getName() << endl;
  }
#endif

//  bsuccess = bsuccess && updateIntermediateVars();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       : update secondary state for one external device
//                 instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------

bool Instance::updateSecondaryState()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/06/04
//-----------------------------------------------------------------------------

bool Instance::getBreakPoints ( vector<N_UTL_BreakPoint> & breakPointTimes)
{
  // if necessary, initialize
  initialize ();
  bool bsuccess = extCodePtr_->getBreakPoints (breakPointTimes);
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::runExternalDevice
// Purpose       : actually run the solve of the external device
//
// Special Notes : updateIntermediateVars is called here, instead of
//                 in updatePrimaryState, b/c these devices need to do
//                 synchronized loads in parallel.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Elec. & Micro Systems
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Instance::runExternalDevice()
{
  bool bsuccess = true;
  bsuccess = bsuccess && updateIntermediateVars();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupVoltageLimiting_
// Purpose       :
// Special Notes : The limits are hardwired for now, but should be made to
//                 be more flexible later.  Specifically, the limits should
//                 take advantage of gradient information.
//
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 03/24/06
//-----------------------------------------------------------------------------
bool Instance::setupVoltageLimiting_ ()
{
  bool bsuccess = true;
  int iBC(0),i(0);

#ifdef Xyce_DEBUG_DEVICE
  cout << "In setupVoltageLimiting" << endl;
#endif

  origFlag=true;

  map<string,double>::iterator iterM = voltageInputMap_.begin();
  map<string,double>::iterator endM  = voltageInputMap_.end ();

  int voltSize = voltageOrig_.size();
  // save the orig, unlimited voltages:
  for(i=0;i<voltSize;++i,++iterM)
  {
    voltageOrig_[i] = iterM->second;
#ifdef Xyce_DEBUG_DEVICE
    cout << "voltageOrig_["<<i<<"] = " << voltageOrig_[i] << endl;
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  cout << getName() << "  Before Limiting:"<<endl;
  iterM = voltageInputMap_.begin();
  for( ; iterM != endM; ++iterM)
  {
    cout << "Blim: ";
    cout << iterM->first;
    cout << "\t";
    cout << iterM->second;
    cout << endl;
  }
#endif

  // grab the "old" voltages from either the "last call" vector, or the
  // state vector.
  if (getSolverState().newtonIter == 0)
  {
    newtonIterOld=0;
#ifdef Xyce_DEBUG_DEVICE
    cout << "Inside the newtonIter==0... " << endl;
#endif

    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
#ifdef Xyce_DEBUG_DEVICE
      cout << "Inside !dcopFlag || (locaEnabledFlago && dcopFlag) " << endl;
#endif
      for(i=0;i<voltSize;++i)
      {
        int li_state = voltageStateID_[i];
        //voltageOld_[i] = (*(*extData.currStaVectorPtrPtr))[li_state];
        voltageOld_[i] = (*extData.currStaVectorPtr)[li_state];
      }
    }
    else  // no history, just set old voltages to the current ones.
    {
#ifdef Xyce_DEBUG_DEVICE
      cout << "Inside no-history" << endl;
#endif
      iterM = voltageInputMap_.begin();
      for(i=0;i<voltSize;++i,++iterM)
      {
        voltageOld_[i] = iterM->second;
      }
    }
  }
  else if (getSolverState().newtonIter != 0 && getSolverState().newtonIter != newtonIterOld)
  {
#ifdef Xyce_DEBUG_DEVICE
    cout << "Inside lastCall. " << endl;
#endif
    for(i=0;i<voltSize;++i)
    {
      voltageOld_[i] = voltageLastCall_[i];
    }
  }

  bool setInitJct = false;
  if (getSolverState().initJctFlag && getDeviceOptions().voltageLimiterFlag)
  {
#ifdef Xyce_DEBUG_DEVICE
    cout << "Inside the initJctFlag==true." << endl;
#endif
    if (initJctGiven_)
    {
#ifdef Xyce_DEBUG_DEVICE
      cout << "Inside the local initJctGiven_==true." << endl;
#endif
      origFlag=false;
      setInitJct = true;

      iterM = voltageInputMap_.begin();
      for(i=0;i<voltSize;++i,++iterM)
      {
        iterM->second = voltageNodeVec[i]->initVal;
#ifdef Xyce_DEBUG_DEVICE
        cout << "initVal : " << iterM->first << " = ";
        cout << voltageNodeVec[i]->initVal << endl;
#endif
      }
    }
  }

  if (!setInitJct)
  {
    iterM = voltageInputMap_.begin();
    endM = voltageInputMap_.end ();
    for(i=0; iterM != endM; ++iterM, ++i)
    {
      double v1     = iterM->second;
      double v1_old = voltageOld_[i];
      double delV1 = v1 - v1_old;

      // Note: For now we are using hardwired values.
      // In the future, this should
      // be a more sophisticated limiter(s).
      //double limValHigh = voltageNodeVec[i].limValHigh;
      double limValHigh = 0.15;
      if ( delV1 > +limValHigh)
      {
        v1 = v1_old + limValHigh;
        origFlag=false;
      }
      //double limValLow = voltageNodeVec[i].limValLow;
      double limValLow = 0.10;
      if ( delV1 < -limValLow)
      {
        v1 = v1_old - limValLow;
        origFlag=false;
      }
      iterM->second = v1;
    }
  }

  // update the "old" variables:
  if (getSolverState().newtonIter != 0 && getSolverState().newtonIter != newtonIterOld)
  {
    newtonIterOld = getSolverState().newtonIter;
  }

  // save the junction voltages in case the newton iteration is about
  // to be declared as successful.  If so, these "lastcall" variables will
  // get copied into their "old" versions.
  iterM = voltageInputMap_.begin();
  for(i=0;i<voltSize;++i,++iterM)
  {
    voltageLastCall_[i] = iterM->second;
  }

#ifdef Xyce_DEBUG_DEVICE
  cout << getName() << "  After Limiting:"<<endl ;
  iterM = voltageInputMap_.begin();
  for( ; iterM != endM; ++iterM)
  {
    cout << "Alim: ";
    cout << iterM->first;
    cout << "\t";
    cout << iterM->second;
    cout << endl;
  }
  if (origFlag)
    cout << "origFlag is TRUE" << endl;
  else
    cout << "origFlag is FALSE" << endl;
#endif

  iterM = voltageInputMap_.begin();
  for(i=0;i<voltSize;++i,++iterM)
  {
    voltageDiff_[i] = iterM->second - voltageOrig_[i];
  }

  // Finally, save voltages in the state vector:
  N_LAS_Vector * staVectorPtr;
  //staVectorPtr = *(extData.nextStaVectorPtrPtr);
  staVectorPtr = extData.nextStaVectorPtr;
  for(i=0;i<voltSize;++i)
  {
    int li_state = voltageStateID_[i];
    (*staVectorPtr)[li_state] = voltageOld_[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcVoltLimFactors_
//
// Purpose       : This function calculates the factors that will need to
//                 be loaded into the RHS vector to support voltage limiting.
//
// Special Notes : This function does a very small matvec, to produce
//                 jdxp terms.
//
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 03/24/06
//-----------------------------------------------------------------------------
bool Instance::calcVoltLimFactors_ ()
{
  bool bsuccess = true;

  int i=0;
  int j=0;
  int voltSize = voltageDiff_.size();
  for (;i<voltSize;++i)
  {
    double fac = 0.0;
    for (j=0;j<voltSize;++j)
    {
      fac += conductanceJacobian_[i][j] * voltageDiff_[j];
    }
    voltageFactor_[i] = fac;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::initialize
// Purpose       : If this device hasn't been initialized yet, do it here.
//                 If it has already been initialized, then this function is
//                 a no-op.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Elec. & Micro Systems
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool Instance::initialize ()
{
  // If necessary, initialize the external device/code.
  if (!initializeFlag_)
  {
    if( locallyConnected_ )
    {
      // Do some initializations related to the voltageNodeVec.
      // I'm not sure where else to put this, so put it here.
      // bool nodeGiven = given("NODE");
      // vector composites don't seem to use the given function correctly...
      if (nodeGiven_)
      {
        // Check that the voltageNodeVec is the right size.
        int sizeVoltMap = voltageNodeVec.size();
        if (sizeVoltMap != numExtVars)
        {
          string msg = "Instance::initialize: \n";
          msg += "Number of specified nodes != numExtVars";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        }

        // The name parameter was required, so force the voltageInputMap
        // to have the same keys as the voltageNodeVec.
        voltageInputMap_.clear();
        int i=0; int vnSize=voltageNodeVec.size();
        for (i=0;i<vnSize;++i)
        {
          voltageInputMap_[voltageNodeVec[i]->vsrcName] = 0.0;
        }

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          map <string,double>::const_iterator iterVI = voltageInputMap_.begin();
          map <string,double>::const_iterator endVI = voltageInputMap_.end ();
          for (;iterVI!=endVI;++iterVI)
          {
            cout << "voltageInputMap["<<iterVI->first<<"] = ";
            cout << iterVI->second << endl;
          }
        }
#endif

        // Now check if initial values were specified.
        bool allInitSet = true;
        bool oneOrMoreInitSet = false;
        for (i=0;i<vnSize;++i)
        {
          VoltageNode & vn = *(voltageNodeVec[i]);
          bool tmpBool = vn.given("INITVAL");
          allInitSet = allInitSet && tmpBool;
          oneOrMoreInitSet = oneOrMoreInitSet || tmpBool;
        }

        if (oneOrMoreInitSet && !allInitSet)
        {
          string msg = "Instance::initialize: \n";
          msg += "You have set an initial value on at least one node,\nbut not all of them.  You must either set all or none.\n";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        }

        initJctGiven_ = allInitSet;
      }
      else // if the "node" vector-composite was *not* specified, we need to use the
           // defaults, and make sure that they are set up.
      {
        int sizeVoltMap = voltageNodeVec.size();
        if (sizeVoltMap != numExtVars)
        {
          int i=0;
          voltageNodeVec.resize(numExtVars);
          for (i=0;i<numExtVars;++i)
          {
            voltageNodeVec[i] = new VoltageNode();
            char tmp[16];
            sprintf(tmp,"VCONNECT%04d",i);
            voltageNodeVec[i]->vsrcName = string(tmp);
          }
        }
      }
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "Initializing internal Xyce" << endl;
    }
#endif
    extCodePtr_->initialize(
#ifdef Xyce_PARALLEL_MPI
      comm_
#endif
                            );

    startTimeStep ();
    initializeFlag_ = true;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::constructComposite
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/03/06
//-----------------------------------------------------------------------------
CompositeParam *Instance::constructComposite(string & cName, string & pName)
{
  CompositeParam * cpPtr = NULL;

  if (cName == "NODE")
  {
    nodeGiven_ = true;
    VoltageNode *vNodePtr = new VoltageNode ();
    voltageNodeVec.push_back(vNodePtr);
    cpPtr = (static_cast<CompositeParam *> (vNodePtr));
  }
  else
  {
    string msg =
      "Instance::ConstructComposite: unrecognized composite name: ";
    msg += cName;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  return cpPtr;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one external device
//                 instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  //N_LAS_Vector * solPtr = *(extData.nextSolVectorPtrPtr);
  N_LAS_Vector * solPtr = extData.nextSolVectorPtr;
  int i=0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::updateIntermediateVars:" << endl;

    int vnSize = voltageNodeVec.size();

    for (i=0;i<vnSize;++i)
    {
      cout << voltageNodeVec[i];
    }
  }
#endif

  // If necessary, initialize the external device/code.
  initialize ();

  innerSolveStatus_ = false;

#ifdef Xyce_PARALLEL_MPI
  if( comm_->procID() != owningProc_ ) locallyConnected_ = false;
#endif

  if( locallyConnected_ )
  {
    // Determine the input variables, etc., and place them in
    // the voltageInputMap
    for (i=0;i<extLIDVec.size();++i)
    {
      voltageInputMap_[voltageNodeVec[i]->vsrcName] = (*solPtr)[extLIDVec[i]];
    }

    if (voltageLimiterFlag && getDeviceOptions().voltageLimiterFlag)
    {
      setupVoltageLimiting_ ();
    }
  }

#ifdef Xyce_PARALLEL_MPI
  assert( owningProc_ != -1 );
  assert( comm_ != 0 );

  if( comm_->procID() != owningProc_ ) locallyConnected_ = false;

  map<string,double>::iterator iterM = voltageInputMap_.begin();
  map<string,double>::iterator endM  = voltageInputMap_.end ();

  for( ; iterM != endM; ++iterM )
  {
    double val = iterM->second;
    double result = 0.0;
    if( !locallyConnected_ ) val = 0.0;
    comm_->sumAll( &val, &result, 1 );
    iterM->second = result;
  }
#endif

  // Now, have the code perform a short simulation, and return currents
  // and conductances.
  //
  // had to change the first two arguments from
  //   currtime-currTimeStep, currTime
  // to
  //   currTime, currTime + currTimeStep
  // to avoid negative intial times and cases where this is
  // no time difference between the two arguments.
  //

  innerSolveStatus_ = extCodePtr_->simulateStep (
    solState,
    voltageInputMap_,
    currentOutputVector_,
    conductanceJacobian_,
    tlError_
                                                 );

  if( locallyConnected_ )
  {
    // Do some voltlim stuff here, if necessary.
    if (voltageLimiterFlag && getDeviceOptions().voltageLimiterFlag)
    {
      calcVoltLimFactors_ ();
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       : update intermediate variables for one external
//                 device instance
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Instance::updateTemperature( const double & temp )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void Instance::homotopyStepSuccess
(const vector<string> & paramNames,
 const vector<double> & paramVals)
{
  extCodePtr_->homotopyStepSuccess ( paramNames, paramVals);
  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void Instance::homotopyStepFailure ()
{
  extCodePtr_->homotopyStepFailure ();
  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void Instance::stepSuccess (int analysis)
{
  extCodePtr_->stepSuccess(analysis);
  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void Instance::stepFailure (int analysis)
{
  extCodePtr_->stepFailure(analysis);
  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  extCodePtr_->getInitialQnorm (tlError_);
  tle = tlError_;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInnerLoopErrorSum
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::getInnerLoopErrorSum (N_TIA_TwoLevelError & tle)
{
  tle = tlError_;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateStateArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::updateStateArrays ()
{
  bool bsuccess = true;
  if (extCodePtr_)
  {
    bool bs1 = extCodePtr_->updateStateArrays ();
    bsuccess = bsuccess && bs1;
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::startTimeStep ()
{
  bool bsuccess = true;
  if (extCodePtr_)
  {
    bool bs1 = extCodePtr_->startTimeStep (getSolverState().tiInfo);
    bsuccess = bsuccess && bs1;
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setInternalParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool Instance::setInternalParam(string & name1, double val)
{
  bool bsuccess = true;
  initialize ();
  if (extCodePtr_)
  {
    bool bs1 = extCodePtr_->setInternalParam(name1,val);
    bsuccess = bsuccess && bs1;
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
Model::Model(const ModelBlock & MB,
             SolverState & ss1,
             DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1)
{

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
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
// Creation Date : 04/09/05
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
}

} // namespace ExternDevice
} // namespace Device
} // namespace Xyce
