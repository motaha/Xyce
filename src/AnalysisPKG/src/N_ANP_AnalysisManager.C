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

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_ANP_AnalysisManager.C,v $
// Purpose       : This file contains the functions which define the time
//                 domain & integration algorithm classes.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.87.2.5 $
// Revision Date  : $Date: 2013/10/03 17:23:31 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#include <iostream>
#include <ctime>

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

#include <sstream>

// ----------   Xyce Includes   ----------

#include <N_TIA_Assembler.h>
#include <N_ANP_AnalysisManager.h>

#include <N_ANP_Transient.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_Step.h>
#include <N_ANP_Dakota.h>
#include <N_ANP_MPDE.h>
#include <N_ANP_HB.h>
#include <N_ANP_AC.h>
#include <N_ANP_MOR.h>

#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TIAParams.h>

#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TwoLevelError.h>

#include <N_TIA_MPDEInterface.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_System.h>
#include <N_LAS_LAFactory.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_NLS_Manager.h>

#include <N_LOA_Loader.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>
#include <N_UTL_BreakPoint.h>

#include <N_IO_OutputMgr.h>
#include <N_ANP_OutputMgrAdapter.h>

#include <N_IO_RestartMgr.h>

#include <N_IO_CmdParse.h>

#include <N_MPDE_Manager.h>

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::N_ANP_AnalysisManager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------

N_ANP_AnalysisManager::N_ANP_AnalysisManager(N_IO_CmdParse & cp, N_ANP_AnalysisInterface * anaIntPtr_tmp)
  :
  calledBeforeTwoLevelTran_(false),
  breakPointRestartStep(0),
  tiaParams(cp),
  anaIntPtr(anaIntPtr_tmp,false),
  analysis(ANP_MODE_TRANSIENT),
  analysisParamsRegistered(false),
  currentMode_(0),
  firstTime(true),
  oldPercentComplete(0.0),
  startSimTime(-1.0),
  switchIntegrator_(false),
  initializeAllFlag_(false),
  startTRANtime(0.0),
  stepLoopFlag_(false),
  stepLoopInitialized_(false),
  dcLoopInitialized_(false),
  gui_(false),
  daeStateDerivFlag_(true),
  initializeSolvers_mixedSignal_(false),
  dcLoopSize_(0),
  sweepSourceResetFlag_(true),
  progressFlag_(true),
  solverStartTime_(0.0),
  dakotaRunFlag_(false),
  dakotaIterationNumber_(0),
  saveTime_(0.0),
  saveTimeGiven_(false),
  saveFlag_(false),
  savedAlready_(false),
  dcopRestartFlag_(false),
  dotOpSpecified_(false),
  initialOutputInterval_(0.0),
  nextOutputTime_(0.0),
  initialRestartInterval_(0.0),
  nextRestartSaveTime_(0.0),
  blockAnalysisFlag_(false),
  hbFlag_(false),
  mpdeFlag_(false),
  sensFlag_(false),
  commandLine_(cp)
{
  gui_ = commandLine_.argExists("-gui");

  // check for maximum order on the command line
  if ( commandLine_.argExists ("-maxord") )
  {
    tiaParams.maxOrder =
      atoi( commandLine_.getArgumentValue( "-maxord" ).c_str() );

    if ( tiaParams.maxOrder < 1) tiaParams.maxOrder = 1;
    if ( tiaParams.maxOrder > 5) tiaParams.maxOrder = 5;
  }

  // Create the MPDEIface object and register it with tiaControl
  tiaMPDEIfacePtr_ = rcp(new N_TIA_MPDEInterface(tiaParams));
  tiaMPDEIfacePtr_->registerTIAControl(this);

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::~N_ANP_AnalysisManager
//
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
N_ANP_AnalysisManager::~N_ANP_AnalysisManager()

{
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::resetAll
//
// Purpose       : just like a destructor without the death
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::resetAll()
{
  // dsPtr_ is created in initializeAll
  dsPtr_ = Teuchos::null;
  // secPtr_ is created in initializeAll
  secPtr_ = Teuchos::null;
  // wimPtr is created in initializeAll
  wimPtr = Teuchos::null;

  // xyceTranTimerPtr_ is created in run
  xyceTranTimerPtr_ = Teuchos::null;

  // tiaMPDEIfacePtr_'s copy to dsPtr_ and secPtr_ are reset in intializeAll, we don't need to delete it.

  // assemblerPtr is created in initializeAll
  assemblerPtr = Teuchos::null;

  // Reset step statistics to zero.
  primaryAnalysisObject_->resetAll ();

  initializeAllFlag_ = false;

}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getBlockAnalysisFlag
// Purpose       : "get" function for MPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getBlockAnalysisFlag () const
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return blockAnalysisFlag_;
  }
  return mpdeMgrPtr_->blockAnalysisFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getMPDEFlag
// Purpose       : "get" function for MPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getMPDEFlag ()
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return mpdeFlag_;
  }
  return mpdeMgrPtr_->getMPDEFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getMPDEStartupFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getMPDEStartupFlag ()
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return false;
  }
  return mpdeMgrPtr_->getMPDEStartupFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getMPDEIcFlag
// Purpose       : get function for MPDE initial condition flag (true if MPDE & IC
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getMPDEIcFlag()
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return false;
  }
  return mpdeMgrPtr_->getMPDEIcFlag();
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getWaMPDEFlag ()
// Purpose       : "get" function for WaMPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getWaMPDEFlag ()
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return false;
  }
  return mpdeMgrPtr_->getWaMPDEFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::isPaused
//
// Purpose       : return the true if the simulation is currently paused
//
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical and Microsystems Simulation
// Creation Date : 02/18/2008
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::isPaused()
{
  return secPtr_->isPauseTime();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::initializeAll
// Purpose       : This function performs the final initializations.  Mainly,
//                 it initializes the two data container classes, DataStore and
//                 StepErrorControl.  It also registers the neccessary vectors
//                 with the LAS system class.
// Special Notes : This function should *only* be called after all the
//                 registrations have all been performed.  In particular, the
//                 N_LAS_System class *must* be registered before this function
//                 is called.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::initializeAll()
{
  string msg;

  if (Teuchos::is_null(lasSysPtr))
  {
    msg = "N_ANP_AnalysisManager::initializeAll ";
    msg += " You need to register the LAS system first.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  tiaParams.solutionSize = lasSysPtr->getSolutionSize();
  tiaParams.stateSize    = lasSysPtr->getStateSize();

  // allocate data store class, which will allocate all the vectors.
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    std::cout << "N_ANP_AnalysisManager::initializeAll.  " ;
    std::cout << "  maximum order = " << tiaParams.maxOrder << std::endl;
  }
#endif

  dsPtr_ = rcp(new N_TIA_DataStore(&tiaParams, &*lasSysPtr));
  secPtr_ = rcp(new N_TIA_StepErrorControl(commandLine_,*this,tiaParams));

  // Now that data store has been created, we can also create the
  // working integration method object.
  wimPtr = rcp(new N_TIA_WorkingIntegrationMethod (tiaParams,*secPtr_,*dsPtr_));

  secPtr_->registerWIMPtr(&*wimPtr);

  assemblerPtr = rcp(new N_TIA_DAE_Assembler(*dsPtr_, *loaderPtr, *wimPtr, *pdsMgrPtr, daeStateDerivFlag_));

  tiaMPDEIfacePtr_->registerTIADataStore(&*dsPtr_);
  tiaMPDEIfacePtr_->registerTIAStepErrorControl(&*secPtr_);
  // The final time has to be forced to be a "paused" breakpoint,
  // after any registerTIAParams call.
  // RLS: May need to conditionally call this only when this is an MPDE run
  // not sure if we know that when this routine is called.
  //setPauseTime(tiaParams.finalTime);

  registerOutputIntervals();
  registerRestartIntervals();

  // register current:
  lasSysPtr->registerCurrStaVector(&(dsPtr_->currStatePtr));
  lasSysPtr->registerCurrSolVector(&(dsPtr_->currSolutionPtr));

  // register next:
  lasSysPtr->registerNextStaVector(&(dsPtr_->nextStatePtr));
  lasSysPtr->registerNextSolVector(&(dsPtr_->nextSolutionPtr));

  // register last:
  lasSysPtr->registerLastStaVector(&(dsPtr_->lastStatePtr));
  lasSysPtr->registerLastSolVector(&(dsPtr_->lastSolutionPtr));

  lasSysPtr->registerFlagSolVector(&(dsPtr_->flagSolutionPtr));

  // register temporaries:
  lasSysPtr->registerTmpSolVector(&(dsPtr_->tmpSolVectorPtr));
  lasSysPtr->registerTmpStaVector(&(dsPtr_->tmpStaVectorPtr));
  lasSysPtr->registerTmpStaDerivVector(&(dsPtr_->tmpStaDerivPtr));
  lasSysPtr->registerTmpStaDivDiffVector(&(dsPtr_->tmpStaDivDiffPtr));

  // register next derivatives:
  lasSysPtr->registerNextSolDerivVector(&(dsPtr_->nextSolutionDerivPtr));
  lasSysPtr->registerNextStaDerivVector(&(dsPtr_->nextStateDerivPtr));

  // register the device mask
  lasSysPtr->registerDeviceMaskVector(&(dsPtr_->deviceMaskPtr));

  // Get the RHS and the Jacobian
  dsPtr_->JMatrixPtr    = lasSysPtr->getJacobianMatrix();
  dsPtr_->RHSVectorPtr  = lasSysPtr->getRHSVector();

  //RefCountPtr<N_TIA_DataStore> dsDaePtr = Teuchos::rcp_dynamic_cast<N_TIA_DataStore>(dsPtr_);
  // DAE formulation vectors
  lasSysPtr->registerDAEQVector     ( dsPtr_->daeQVectorPtr );
  lasSysPtr->registerDAEFVector     ( dsPtr_->daeFVectorPtr );

  // DAE formulation matrices
  lasSysPtr->registerDAEdQdxMatrix  ( dsPtr_->dQdxMatrixPtr );
  lasSysPtr->registerDAEdFdxMatrix  ( dsPtr_->dFdxMatrixPtr );

  // Get the RHS and the Jacobian
  //dsPtr_->JMatrixPtr    = lasSysPtr->getJacobianMatrix();
  //dsPtr_->RHSVectorPtr  = lasSysPtr->getRHSVector();

  dsPtr_->dFdxdVpVectorPtr = lasSysPtr->getdFdxdVpVector ();
  dsPtr_->dQdxdVpVectorPtr = lasSysPtr->getdQdxdVpVector ();

  dsPtr_->limiterFlag = loaderPtr->getLimiterFlag ();

  // This should probably be moved elsewhere later.  If the user has
  // specified that steps should only be accepted when the nonlinear solver
  // truly converges, then the "nearConvergence" return code needs to be
  // negative.  ERK. 7/02/03
  if ( !(tiaParams.nlNearConvFlag) )
  {
    N_NLS_ReturnCodes retCodes;
    retCodes.nearConvergence = -3;
    nlsMgrPtr->setReturnCodes(retCodes);
  }

  // same for the "small update" case.
  if ( !(tiaParams.nlSmallUpdateFlag) )
  {
    N_NLS_ReturnCodes retCodes;
    retCodes.smallUpdate = -4;
    nlsMgrPtr->setReturnCodes(retCodes);
  }

  // check if analysis was specified.  If not, but .OP was specified,
  // then set up a DC calculation.
  if ( !analysisParamsRegistered )
  {
    if ( dotOpSpecified_ )
    {
      analysis = ANP_MODE_DC_SWEEP;
      analysisParamsRegistered = true;
    }
    else // flag an error.
    {
      string msg =  "No analysis statement in the netlist\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
    }
  }

  // Allocate analysis objects, and also set up params.
  allocateAnalysisObject_ ();

  initializeAllFlag_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::run
// Purpose       : Execute the control loop for the set analysis type.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::run()
{
  string msg;

  if (initializeAllFlag_ == false)
  {
    msg = "N_ANP_AnalysisManager::run - "
      "you need to call the initializeAll function first\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  if (analysisParamsRegistered == false)
  {
    msg = "N_ANP_AnalysisManager::run: "
      "no analysis statement in the netlist\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
  }

  // now that the problem is set up, we can tell the datastore what the var type are
  // this info can be used by the time integrator to control error in a way that is
  // appropriate for different variable types (like V's and I's).
  vector<char> varTypes;
  topoMgrPtr->returnVarTypeVec( varTypes );
  dsPtr_->varTypeVec = varTypes;

  // This checks to make sure that all quantities in the .print are valid before
  // we continue.
  outputMgrAdapterRCPtr_->check_output( analysis ); //, *dsPtr_->currSolutionPtr, *dsPtr_->currStatePtr, *dsPtr_->currStorePtr );

#ifdef Xyce_VERBOSE_TIME
  tiaParams.printParams(anpAnalysisModeToNLS(analysis));
#endif

#ifndef Xyce_NO_MASKED_WRMS_NORMS
  if (loaderPtr->loadDeviceMask())
  {
     //std::cout << " Nontrivial mask! " << std::endl;
  }
  else
  {
    // std::cout << " Trivial mask! " << std::endl;
  }
#endif
//

  // 6/24/2010: ERK Note: this needs a refactor!
  // For HB and MPDE, it is necessary to reallocate different analysis types as the simulation
  // progresses from the initial condition phase to the block analysis phase.  For example,
  // in MPDE, it is common to solve a DC, then a series of transients, and then the full MPDE
  // simulation.
  //
  // I recently moved some stuff from the analysis manager down into specific analysis type
  // classes, which is where they really belong.  However, doing this required that the
  // order of setup change somewhat.  DC and sweep parameters now cannot be processed until
  // the N_ANP_DCSweep or N_ANP_Step classes are allocated.  They are now primarily allocated
  // in the intializeAll function.  They need to be allocated prior to this ::run function,
  // because they need to happen before restart files are read in.
  //
  // However, to preserve MPDE and HB, there has to be an option of reallocating analysis
  // classes here.  So, that is still the case, but this should be refactored to make
  // it cleaner.
  if ( getBlockAnalysisFlag () )
  {
    allocateAnalysisObject_ ();
  }

  solverStartTime_ = elapsedTimerPtr_->elapsedTime();

  // Start the solvers timers.
  xyceTranTimerPtr_ = rcp(new N_UTL_Timer(*(pdsMgrPtr->getPDSComm())));

  bool runStatus = analysisObject_->run();

  if (tiaParams.condTestFlag)
  {
    conductanceTest ();
  }

  return runStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::allocateAnalysisObject_
// Purpose       : Allocate analysis objects, and also setup params.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/10
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::allocateAnalysisObject_ ()
{
  string msg("");

  if( !tiaParams.resume )
  {
    if (analysis == ANP_MODE_TRANSIENT )
    {
      analysisObject_ = Teuchos::rcp(new N_ANP_Transient(this));
      analysisObject_->setAnalysisParams(tranParamsBlock);
      secPtr_->resetAll();
    }
    else if (analysis == ANP_MODE_DC_SWEEP)
    {
      analysisObject_ = Teuchos::rcp(new N_ANP_DCSweep(this));
      for (int i=0;i<(int)dcParamsBlockVec.size();++i)
      {
        analysisObject_->setAnalysisParams(dcParamsBlockVec[i]);
      }
    }
    else if (analysis == ANP_MODE_MPDE)
    {
      analysisObject_ = Teuchos::rcp(new N_ANP_MPDE(this));
      analysisObject_->setAnalysisParams(mpdeParamsBlock);
    }
    else if (analysis == ANP_MODE_HB)
    {
      analysisObject_ = Teuchos::rcp(new N_ANP_HB(this));
      analysisObject_->setAnalysisParams(hbParamsBlock);
      Teuchos::rcp_dynamic_cast<N_ANP_HB>(analysisObject_)->setHBOptions(hbOptionsBlock);
      Teuchos::rcp_dynamic_cast<N_ANP_HB>(analysisObject_)->setHBLinSol(hbLinSolBlock);
      Teuchos::rcp_dynamic_cast<N_ANP_HB>(analysisObject_)->setLinSol(linSolBlock);
      setHBFlag( true );
    }
    else if ( analysis == ANP_MODE_AC)
    {
      analysisObject_ = Teuchos::rcp(new N_ANP_AC(this));
      analysisObject_->setAnalysisParams(acParamsBlock);
    }
    else if ( analysis == ANP_MODE_MOR)
    {
      analysisObject_ = Teuchos::rcp(new N_ANP_MOR(this));
      analysisObject_->setAnalysisParams(morParamsBlock);
    }
    else
    {
      msg = "N_ANP_AnalysisManager::initializeAll : unknown type of analysis\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }

    if( !is_null(analysisObject_) )
    {
      analysisObject_->setParamsWithOutputMgrAdapter ( outputMgrAdapterRCPtr_ );
    }

    // need to throw an error here as this really shouldn't happen.
    if( is_null( analysisObject_ ) )
    {
      msg = "N_ANP_AnalysisManager::initializeAll : unable to allocate analysis type.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }

    // ultimately sensitivity calculation should get its own analysis object.
    if (sensFlag_)
    {
      analysisObject_->setSensFlag();
    }

    // The primary analysis object is stored in case there is an outer
    // loop analysis such as .STEP or .DAKOTA.   Many accessors such
    // as "getStepNumber" need to access this object rather than the
    // step object.
    primaryAnalysisObject_ = analysisObject_;

    if( stepLoopFlag_ )
    {
      stepAnalysisTarget_ = analysisObject_;
      analysisObject_ = Teuchos::rcp(new N_ANP_Step(this, stepAnalysisTarget_.get()));
      for (int i=0;i<(int)stepParamsBlockVec.size();++i)
      {
        analysisObject_->setAnalysisParams(stepParamsBlockVec[i]);
      }
      analysisObject_->setParamsWithOutputMgrAdapter ( outputMgrAdapterRCPtr_ );
    }

    if (dakotaRunFlag_ && is_null(dakotaAnalysisTarget_) )
    {
      dakotaAnalysisTarget_ = analysisObject_;
      analysisObject_ = Teuchos::rcp(new N_ANP_Dakota(this, dakotaAnalysisTarget_.get()));
      analysisObject_->setAnalysisParams(dakotaParamsBlock);
      analysisObject_->setParamsWithOutputMgrAdapter ( outputMgrAdapterRCPtr_ );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::printLoopInfo
// Purpose       : Prints out time loop information.
// Special Notes : Prints stats from save point start to save point finish.
//                 Special case 0,0 is entire run to this point
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::printLoopInfo(int start, int finish)
{
  return primaryAnalysisObject_->printLoopInfo (start,finish);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::testDCOPOutputTime_
// Purpose       : Similar to testOutputTime, except that this is for
//                 DCOP restart files.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::testDCOPOutputTime_()
{
  bool flag(true);

  if( !dcopRestartFlag_ )
  {
    flag = false;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::testSaveOutputTime_
// Purpose       : Similar to testOutputTime, except that this is for
//                 .SAVE files.
// Special Notes : Only outputs 1x.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::testSaveOutputTime_()
{
  bool flag(true);

  if( !saveFlag_ )
  {
    flag = false;
  }
  else if (secPtr_->currentTime < saveTime_)
  {
    flag = false;
  }
  else if (savedAlready_)
  {
    flag = false;
  }

  if (flag==true)
  {
    savedAlready_ = true;
    std::cout <<"Calling SAVE outputs!" <<std::endl;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::testRestartSaveTime_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel ComputationalSciences.
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::testRestartSaveTime_()
{
  bool flag;

#ifdef Xyce_DEBUG_RESTART
  std::cout << "TESTING FOR RESTART SAVE\n";
  std::cout << "------------------------\n";
  std::cout << "secPtr_->currentTime: " << secPtr_->currentTime << std::endl;
  std::cout << "nextSaveTime: " << nextRestartSaveTime_ << std::endl;
  std::cout << "initialRestartInterval_: " << initialRestartInterval_ << std::endl;
  if (!(restartIntervals_.empty()))
  {
    std::cout << "First restart interval: " << restartIntervals_[0].first << std::endl;
  }
  else
  {
    std::cout << "restartIntervals_ is empty" << std::endl;
  }
#endif

  if (initialRestartInterval_ == 0.0)
  {
    flag = false;
  }
  else if (secPtr_->currentTime < nextRestartSaveTime_)
  {
    flag = false;
  }
  else if (restartIntervals_.empty())
  {
    while (nextRestartSaveTime_ <= secPtr_->currentTime)
    {
      nextRestartSaveTime_ += initialRestartInterval_;
    }
    flag = true;
  }
  else if (secPtr_->currentTime < restartIntervals_[0].first)
  {
    while (nextRestartSaveTime_ <= secPtr_->currentTime)
    {
      nextRestartSaveTime_ += initialRestartInterval_;
    }
    if (nextRestartSaveTime_ > restartIntervals_[0].first)
    {
      nextRestartSaveTime_ = restartIntervals_[0].first;
    }
    flag = true;
  }
  else
  {
    pair<double, double> currInterval, nextInterval;
    int size = restartIntervals_.size();
    for (int i = 0; i < size; ++i)
    {
      if (restartIntervals_[i].first <= secPtr_->currentTime)
      {
        currInterval = restartIntervals_[i];
        if ((i+1) < (int)restartIntervals_.size())
        {
          nextInterval = restartIntervals_[i+1];
        }
      }
    }
    int step = static_cast <int> ((secPtr_->currentTime-currInterval.first) /
                                  currInterval.second);
    nextRestartSaveTime_ = currInterval.first + (step+1)*currInterval.second;

    if (nextInterval.first && (nextInterval.first!=currInterval.first)
        && (nextRestartSaveTime_>=nextInterval.first))
    {
      nextRestartSaveTime_ = nextInterval.first;
    }
    flag = true;
  }

#ifdef Xyce_DEBUG_RESTART
  std::cout << "new nextSaveTime: " << nextRestartSaveTime_ << std::endl;
  std::cout << "restart flag: " << flag << std::endl;
  std::cout << "-----------------------\n";
#endif

  return flag;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::partialTimeDerivative
// Purpose       : Returns the current partial time derivative for either the
//                 solution or state vector.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::partialTimeDerivative()
{
  double partDT = wimPtr->partialTimeDeriv();

  // Add a check here to try to prevent capacitive spiral of death.
  // This is still a "research" option, so is not on by default.
  if (tiaParams.jacLimitFlag)
  {
#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      std::cout << "N_ANP_AnalysisManager::partialTimeDerivative.";
      std::cout << "   Using jac limit = " << tiaParams.jacLimit << std::endl;
    }
#endif
    if (partDT > tiaParams.jacLimit)
    {
      partDT = tiaParams.jacLimit;
    }
  }

  return partDT;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getBreakpointTol
// Purpose       : Returns the breakpoint tolerance.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/7/02
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getBreakpointTol()
{
  return N_UTL_BreakPoint::getBPTol();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setBreakpointTol
// Purpose       : Sets the breakpoint tolerance.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/7/02
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setBreakpointTol(double bptol)
{
  N_UTL_BreakPoint::setBPTol(bptol);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerRestartIntervals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerRestartIntervals()
{
  double initInt;
  vector< pair<double,double> > intPairs;
  restartPtr->getRestartIntervals(initInt, intPairs);
  initialRestartInterval_ = initInt;
  nextRestartSaveTime_    = secPtr_->initialTime;
  restartIntervals_       = intPairs;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerOutputIntevals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerOutputIntervals()
{
  double initInt;
  vector< pair<double,double> > intPairs;
  outputMgrAdapterRCPtr_->getOutputIntervals( initInt, & intPairs );
  initialOutputInterval_ = initInt;
  nextOutputTime_ = secPtr_->initialTime;
  outputIntervals_ = intPairs;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setTranAnalysisParams
// Purpose       : Sets transient analysis parameters (from .TRAN)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setTranAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysis = ANP_MODE_TRANSIENT;
  tranParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setDCAnalysisParams
// Purpose       : Sets the DC sweep calculation parameters (from .DC)
//
// Special Notes : This function will be called multiple times if there is more
//                 than one sweep variable.  The parser separates each variable
//                 into separate option blocks prior to calling this function.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setDCAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysis = ANP_MODE_DC_SWEEP;
  // before saving a copy of paramsBlock, check and see if it already
  // is in the dcParamsBlockVec.  If so then replace the original copy.
  // This replacement is necessary if the first copy had an expression
  // element that was resolved in later parsing.

  bool foundMatch = false;
  vector<N_UTL_OptionBlock>::iterator paramsBlockVecItr = dcParamsBlockVec.begin();
  vector<N_UTL_OptionBlock>::iterator paramsBlockVecEnd = dcParamsBlockVec.end();
  while( paramsBlockVecItr != paramsBlockVecEnd )
  {
    if( paramsBlockVecItr->compareParamLists( paramsBlock ) )
    {
      // these are the same
      foundMatch = true;
      break;
    }
    paramsBlockVecItr++;
  }

  if( foundMatch )
  {
    // replace the existing one with the new one
    *paramsBlockVecItr = paramsBlock;
  }
  else
  {
    // save the new one.
    dcParamsBlockVec.push_back (paramsBlock); // save a copy for later.
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setOPAnalysisParams
// Purpose       : Handle OP statement. (.OP)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setOPAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  dotOpSpecified_ = true;
  opParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setSTEPAnalysisParams
// Purpose       : Sets the STEP calculation parameters. (from .STEP)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/03
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setSTEPAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  stepLoopFlag_ = true;
  // before saving a copy of paramsBlock, check and see if it already
  // is in the stepParamsBlockVec.  If so then replace the original copy.
  // This replacement is necessary if the first copy had an expression
  // element that was resolved in later parsing.
  bool foundMatch = false;
  vector<N_UTL_OptionBlock>::iterator paramsBlockVecItr = stepParamsBlockVec.begin();
  vector<N_UTL_OptionBlock>::iterator paramsBlockVecEnd = stepParamsBlockVec.end();
  while( paramsBlockVecItr != paramsBlockVecEnd )
  {
    if( paramsBlockVecItr->compareParamLists( paramsBlock ) )
    {
      // these are the same
      foundMatch = true;
      break;
    }
    paramsBlockVecItr++;
  }

  if( foundMatch )
  {
    // replace the existing one with the new one
    *paramsBlockVecItr = paramsBlock;
  }
  else
  {
    // save the new one.
    stepParamsBlockVec.push_back (paramsBlock); // save a copy for later.
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setSaveOptions
// Purpose       : Sets the Save parameters.
// Special Notes : Most of these parameters are handled in the output manager,
//                 rather than here.  So, most params are a noop here, except
//                 for "TIME".
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setSaveOptions(
  const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    std::cout << "In N_ANP_AnalysisManager::setSaveOptions" << std::endl;
  }
#endif

  saveFlag_ = true;

  list<N_UTL_Param>::const_iterator iterPL = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
    std::cout << "iterPL->tag = " << iterPL->tag() << std::endl;
#endif
    if (iterPL->tag() == "TYPE")
    {
      // noop
    }
    else if (iterPL->tag() == "FILE")
    {
      // noop
    }
    else if (iterPL->tag() == "TIME")
    {
      saveTime_ = iterPL->dVal();
      saveTimeGiven_ = true;
    }
    else if (iterPL->tag() == "LEVEL")
    {
      // noop
    }
    else
    {
      // noop.  If it gets here there is an error on the .SAVE line.
      // However, the IO manager has a trap for this, so do nothing
      // here.
    }

    ++iterPL;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setACAnalysisParams
// Purpose       : Sets the AC sweep calculation parameters (from .AC)
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setACAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysis = ANP_MODE_AC;
  acParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setMORAnalysisParams
// Purpose       : Sets the MOR calculation parameters (from .MOR)
//
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 05/29/12
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setMORAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysis = ANP_MODE_MOR;
  morParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setMOROptions
// Purpose       :
// Special Notes : These are from '.options mor'
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 05/31/12
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setMOROptions(const N_UTL_OptionBlock & OB)
{
  list<N_UTL_Param>::const_iterator it_tpL;
  list<N_UTL_Param>::const_iterator first = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag()=="METHOD")
    {
      ExtendedString stringVal ( it_tpL->sVal() );
      stringVal.toUpper();
      tiaParams.morMethod = stringVal;
    }
    else if (it_tpL->uTag()=="SAVEREDSYS")
    {
      tiaParams.morSaveRedSys = static_cast<bool>(it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="COMPORIGTF")
    {
      tiaParams.morCompOrigTF = static_cast<bool>(it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="COMPREDTF")
    {
      tiaParams.morCompRedTF = static_cast<bool>(it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="COMPTYPE")
    {
      ExtendedString stringVal ( it_tpL->sVal() );
      stringVal.toUpper();
      tiaParams.morCompType = stringVal;
    }
    else if (it_tpL->uTag()=="COMPNP")
    {
      tiaParams.morCompNP = it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="COMPFSTART")
    {
      tiaParams.morCompFStart = it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="COMPFSTOP")
    {
      tiaParams.morCompFStop = it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="EXPPOINT")
    {
      tiaParams.morExpPoint = it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="SCALETYPE")
    {
      tiaParams.morScaleType = it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="SCALEFACTOR")
    {
      tiaParams.morScaleFactor = it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="SCALEFACTOR1")
    {
      tiaParams.morScaleFactor1 = it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="SPARSIFICATIONTYPE")
    {
      tiaParams.morSparsificationType = it_tpL->iVal();
    }
    else
    {
      string tmp = " ***** ERROR: " + it_tpL->uTag() +
        " is not a recognized model-order reduction option.\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, tmp);
    }
  }

  // If we are computing the transfer function, make sure the frequency range is valid.
  if (tiaParams.morCompOrigTF || tiaParams.morCompRedTF)
  {
    if (tiaParams.morCompFStop < tiaParams.morCompFStart)
    {
       ostringstream err("");
       err << " ***** ERROR: .options mor COMPFSTART = " << tiaParams.morCompFStart << " > " << tiaParams.morCompFStop << " = COMPFSTOP!" << std::endl;
       N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, err.str());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setDCOPRestartParams
// Purpose       : Sets the DCOP restart parameters.
// Special Notes : Most of the dcop restart parameters are used and handled by
//                 the N_IO_OutputMgr class, so this function here doens't need
//                 to do very much.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setDCOPRestartParams(const N_UTL_OptionBlock & OB)
{
  dcopRestartFlag_ = true;

  list<N_UTL_Param>::const_iterator iterPL = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
    std::cout << "iterPL->tag = " << iterPL->tag() << std::endl;
#endif

    if (iterPL->tag() == "INPUT")
    {
      // do nothing, this is handled in the output manager.
    }
    else if (iterPL->tag() == "OUTPUT")
    {
      // do nothing, this is handled in the output manager.
    }
    else if (iterPL->tag() == "TIME")
    {
      saveTime_ = iterPL->dVal();
      saveTimeGiven_ = true;
    }
    else
    {
      // noop.  If it gets here there is an error on the .SAVE line.
      // However, the IO manager has a trap for this, so do nothing
      // here.
    }

    ++iterPL;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setTranOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setTranOptions(const N_UTL_OptionBlock & OB)
{
  list<N_UTL_Param>::const_iterator it_tpL;
  list<N_UTL_Param>::const_iterator first = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag()=="METHOD")
    {
//      tiaParams.integrationMethod=it_tpL->iVal();

      if (it_tpL->isInteger())
        tiaParams.integrationMethod=it_tpL->iVal();
      else
      {

        ExtendedString stringVal ( it_tpL->sVal() );
        stringVal.toUpper();

        if (stringVal == "TRAP" || stringVal == "TRAPEZOIDAL")
          tiaParams.integrationMethod = 7;
        else if (stringVal == "BDF")
          tiaParams.integrationMethod = 6;
        else if (stringVal == "GEAR")
          tiaParams.integrationMethod = 8;
        else
        {
          string msg=
        "N_ANP_AnalysisManager::setTranOptions: unsupported transient method type";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
        }
      }

    }
#ifdef Xyce_DEBUG_ANALYSIS
    else if (it_tpL->uTag()=="CONSTSTEP")
    {
      tiaParams.constantStepSize
        =  static_cast<bool> (it_tpL->iVal());
    }
#endif
    else if (it_tpL->uTag()=="USEDEVICEMAX")
    {
      tiaParams.useDeviceTimeStepMax =  static_cast<bool> (it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="RELTOL")
    {
      tiaParams.relErrorTol=it_tpL->dVal();
      tiaParams.relErrorTolGiven = true;
    }
    else if (it_tpL->uTag()=="ABSTOL")
    {
      tiaParams.absErrorTol=it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="DOUBLEDCOPSTEP")
    {
      tiaParams.doubleDCOPStep=it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="FIRSTDCOPSTEP")
    {
      tiaParams.firstDCOPStep=it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="LASTDCOPSTEP")
    {
      tiaParams.lastDCOPStep=it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="BPENABLE" )
    {
      tiaParams.bpEnable =  static_cast<bool> (it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="RESTARTSTEPSCALE" )
    {
      tiaParams.restartTimeStepScale=it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="EXITTIME" )
    {
      tiaParams.exitTime=it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="EXITSTEP" )
    {
      tiaParams.exitStep=it_tpL->iVal();
    }
    else if (it_tpL->uTag() == "MINTIMESTEPSBP")
    {
      tiaParams.minTimeStepsBP = it_tpL->iVal();
      tiaParams.minTimeStepsBPGiven = true;
    }
    else if (it_tpL->uTag()=="ERROPTION" )
    {
      tiaParams.errorAnalysisOption=it_tpL->iVal();
      // tscoffe/tmei 12/7/07:  If error option = 1 (NO LTE) then make sure minTimeStepBP is enabled.
      if (tiaParams.errorAnalysisOption == 1)
      {
        if (!tiaParams.minTimeStepsBPGiven)
        {
          tiaParams.minTimeStepsBPGiven = true;
        }
      }
    }
    else if (it_tpL->uTag()=="NLNEARCONV" )
    {
      tiaParams.nlNearConvFlag = static_cast<bool> (it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="NLSMALLUPDATE" )
    {
      tiaParams.nlSmallUpdateFlag = static_cast<bool> (it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="JACLIMITFLAG" )
    {
      tiaParams.jacLimitFlag= static_cast<bool> (it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="JACLIMIT" )
    {
      tiaParams.jacLimit = it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="DAESTATEDERIV" )
    {
      daeStateDerivFlag_ = static_cast<bool> (it_tpL->iVal());
    }
    else if (it_tpL->uTag()=="MAXORD" )
    {
/*      if (tiaParams.integrationMethod == 7)
      {
        string msg = "***** ERROR: maxord is not a valid option for variable order Trapezoid method!";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      } */

      tiaParams.maxOrder = it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="MINORD" )
    {
/*      if (tiaParams.integrationMethod == 7)
      {
        string msg = "***** ERROR: minord is not a valid option for variable order Trapezoid method!";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      } */
      tiaParams.minOrder = it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="TIMESTEPSREVERSAL" )
    {
      tiaParams.timestepsReversal =it_tpL->bVal();
    }
    else if (it_tpL->uTag()=="TESTFIRSTSTEP" )
    {
      tiaParams.testFirstStep = static_cast<bool> (it_tpL->iVal ());
    }
    else if (it_tpL->uTag()=="DELMAX" )
    {
      tiaParams.delmax =it_tpL->dVal();
      tiaParams.delmaxGiven = true;
//      tiaParams.delmax =it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="NLMIN" )
    {
      tiaParams.NLmin=it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="NLMAX" )
    {
      tiaParams.NLmax=it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="OUTPUTINTERPMPDE")
    {
      tiaParams.outputInterpMPDE = static_cast<bool> (it_tpL->iVal ());
    }
    else if (it_tpL->uTag()=="NEWLTE")
    {
      tiaParams.newLte =it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="NEWBPSTEPPING")
    {
      tiaParams.newBPStepping = it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="INTERPOUTPUT")
    {
      tiaParams.interpOutputFlag = static_cast<bool> (it_tpL->iVal ());
    }
    else if (it_tpL->uTag()=="CONDTEST")
    {
      tiaParams.condTestFlag = static_cast<bool> (it_tpL->iVal ());
    }
    else if (it_tpL->uTag()=="CONDTESTDEVICENAME")
    {
      tiaParams.condTestDeviceNames.push_back(it_tpL->sVal () );
    }
    else if (it_tpL->uTag() == "DTMIN")
    {
      tiaParams.userSpecMinTimeStep = it_tpL->dVal();
      tiaParams.userSpecMinTimeStepGiven = true;
    }
    else if (it_tpL->uTag() == "PASSNLSTALL")
    {
      tiaParams.passNLStall = it_tpL->bVal();
    }
    else if (it_tpL->uTag() == "FASTTESTS")
    {
      tiaParams.fastTests = it_tpL->bVal();
    }
    else if (it_tpL->uTag() == "MINTIMESTEPRECOVERY")
    {
      tiaParams.minTimeStepRecoveryCounter = it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="VOLTZEROTOL" )
    {
      tiaParams.voltZeroTol=it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="CURRZEROTOL" )
    {
      tiaParams.currZeroTol=it_tpL->dVal();
    }
    else if (it_tpL->uTag()=="DEBUGLEVEL" )
    {
#ifdef Xyce_DEBUG_TIME
      // set (or override) debug levels based on command line options
      if ( commandLine_.argExists( "-tdl" ) )
      {
        tiaParams.debugLevel = atoi( commandLine_.getArgumentValue( "-tdl" ).c_str() );
      }

      else
      {
        tiaParams.debugLevel = it_tpL->iVal();
      }
#endif
    }
    else if (it_tpL->uTag()=="HISTORYTRACKINGDEPTH" )
    {
      tiaParams.historyTrackingDepth = it_tpL->iVal();
    }
    else
    {
      string tmp = " ***** ERROR: " + it_tpL->uTag() +
        " is not a recognized time integration option.\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, tmp);
    }
  }

  if (tiaParams.NLmin > tiaParams.NLmax)
  {
     ostringstream err("");
      err << " ***** ERROR: .options timeint NLMIN = " << tiaParams.NLmin << " > " << tiaParams.NLmax << " = NLMAX!" << std::endl;
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, err.str());
  }

  if (tiaParams.firstDCOPStep < 0) tiaParams.firstDCOPStep = 0;
  if (tiaParams.firstDCOPStep > 1) tiaParams.firstDCOPStep = 1;
  if (tiaParams.lastDCOPStep < 0)  tiaParams.lastDCOPStep  = 0;
  if (tiaParams.lastDCOPStep > 1)  tiaParams.lastDCOPStep  = 1;

  // check for maximum order on the command line
  if ( commandLine_.argExists ("-maxord") )
  {
/*    if (tiaParams.integrationMethod == 7)
    {
      string msg = "***** ERROR: maxord is not a valid option for variable order Trapezoid method!";
       N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    } */

    tiaParams.maxOrder =
      atoi( commandLine_.getArgumentValue( "-maxord" ).c_str() );
  }

  if ( tiaParams.maxOrder < 1) tiaParams.maxOrder = 1;
  if ( tiaParams.maxOrder > 5) tiaParams.maxOrder = 5;

  if (tiaParams.newLte == true)
  {
    if (tiaParams.relErrorTolGiven != true)
      tiaParams.relErrorTol = 1.0e-3;
  }

//  if (tiaParams.integrationMethod == 7)
//  {
//     string msg = "***** ERROR: maxord is not a valid option for variable order Trapezoid method!"
//     N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
//  }

/*  if (tiaParams.errorAnalysisOption == 1 && tiaParams.maxTimeStepGiven == false)
  {
    tiaParams.maxTimeStep = 0.01 * tiaParams.finalTime;
  }
*/
  return true;

}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setMPDEAnalysisParams
// Purpose       : Sets the MPDE sweep calculation parameters (from .MPDE)
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setMPDEAnalysisParams(const N_UTL_OptionBlock & OB)
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    setupMPDEMgr_();
  }
  bool bsuccess = mpdeMgrPtr_->setMPDEAnalysisParams(OB);
  analysis = ANP_MODE_MPDE;
  analysisParamsRegistered = true;
  mpdeParamsBlock = OB;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setMPDEOptions
// Purpose       :
// Special Notes : from '.options mpdeint'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setMPDEOptions(const N_UTL_OptionBlock & OB)
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    setupMPDEMgr_();
  }
  mpdeMgrPtr_->setMPDEOptions(OB);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setHBAnalysisParams
// Purpose       : Sets the HB sweep calculation parameters (from .HB)
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/30/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setHBAnalysisParams(const N_UTL_OptionBlock & OB)
{
  list<N_UTL_Param>::const_iterator it_tpL;
  list<N_UTL_Param>::const_iterator first = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag() == "FREQ")
    {
      tiaParams.freq= it_tpL->dVal();
      tiaParams.freqGiven = true;
    }
  }

  if (tiaParams.freq <= 0.0 )
  {
    ostringstream ost;
    ost << " N_ANP_AnalysisManager::setHBAnalysisParams: " << std::endl;
    ost << " Frequency of oscillation " << tiaParams.freq
        << " is less than or equal to zero.  "
        << " Check netlist for invalid .HB specification " << std::endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, ost.str());
  }

#ifdef Xyce_DEBUG_ANALYSIS
  string dashedline = "-------------------------------------------------"
                       "----------------------------\n";

  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  HB transient simulation parameters");

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  HB frequency = ",
            tiaParams.freq);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       dashedline);
  }
#endif

  setBlockAnalysisFlag( true );
  devInterfacePtr->setBlockAnalysisFlag( true );
  analysis = ANP_MODE_HB;
  analysisParamsRegistered = true;
  hbParamsBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setHBOptions
// Purpose       :
// Special Notes : from '.options hb'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setHBOptions(const N_UTL_OptionBlock & OB)
{
  hbOptionsBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setLinSol(const N_UTL_OptionBlock & OB)
{
  linSolBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setHBLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setHBLinSol(const N_UTL_OptionBlock & OB)
{
  hbLinSolBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setTRANMPDEOptions
// Purpose       :
// Special Notes : from '.options timeint-mpde'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setTRANMPDEOptions(const N_UTL_OptionBlock & OB)
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    setupMPDEMgr_();
  }
  mpdeMgrPtr_->registerTranMPDEOptions(OB);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setSensOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 6/5/13
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setSensOptions(const N_UTL_OptionBlock & OB)
{
  sensFlag_=true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::completeOPStartStep
// Purpose       : Call to rotate next state to current state following a
//               : constrained DCOP solve when using a previous operating point
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 06/27/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::completeOPStartStep  ( )
{
  bool bsuccess = true;

  bsuccess = dsPtr_->updateStateDataArrays ();
  dsPtr_->setConstantHistory();
  dsPtr_->equateTmpVectors();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::completeHomotopyStep
// Purpose       : Call to do final cleanup after a single,
//                 successful homotopy step.
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::completeHomotopyStep
    ( const vector<string> & paramNames,
      const vector<double> & paramVals,
      N_LAS_Vector * solnVecPtr )
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  string netListFile = commandLine_.getArgumentValue("netlist");
  std::cout << "\n " << netListFile;
  std::cout << " N_ANP_AnalysisManager::completeHomotopyStep " << std::endl;
#endif

  // Rotate the data vectors:
  bool bs1 = dsPtr_->updateStateDataArrays ();    bsuccess = bsuccess && bs1;
  dsPtr_->setConstantHistory();
  dsPtr_->equateTmpVectors();

  // Pass info in to the lower level solver:
  loaderPtr->homotopyStepSuccess (paramNames,paramVals);

  // Call output
  outputMgrAdapterRCPtr_->outputHomotopy( paramNames, paramVals, *solnVecPtr );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::failHomotopyStep
// Purpose       : Call to do final cleanup after a single,
//                 failed homotopy step.
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::failHomotopyStep ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  string netListFile = commandLine_.getArgumentValue("netlist");
  std::cout << "\n " << netListFile;
  std::cout << " N_ANP_AnalysisManager::failHomotopyStep " << std::endl;
#endif

  // Pass info in to the lower level solver:
  loaderPtr->homotopyStepFailure ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setupMPDEMgr_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setupMPDEMgr_()
{

  // Allocate new MPDE Manager
  mpdeMgrPtr_ = rcp(new N_MPDE_Manager(commandLine_));

  // Register MPDE manager with everything
  bool bs1 = true;
  bool bsuccess = true;
  bs1 = outMgrPtr->registerMPDEManager(&*mpdeMgrPtr_);
  bsuccess = bsuccess && bs1;
  mpdeMgrPtr_->registerAnalysisInterface(anaIntPtr);
  mpdeMgrPtr_->registerNonlinearSolver(nlsMgrPtr);
  mpdeMgrPtr_->registerDeviceInterface(devInterfacePtr);
  mpdeMgrPtr_->registerParallelManager(pdsMgrPtr);
  mpdeMgrPtr_->registerTopology(topoMgrPtr);
  mpdeMgrPtr_->registerRestartManager(restartPtr);
  mpdeMgrPtr_->registerOutputManager(outMgrPtr);
  mpdeMgrPtr_->registerApplicationLoader(loaderPtr);
  mpdeMgrPtr_->registerNonlinearEquationLoader(nonlinearEquationLoaderPtr);
  mpdeMgrPtr_->registerApplicationBuilder(appBuilderPtr);
  mpdeMgrPtr_->registerLinearSystem(lasSysPtr);
  mpdeMgrPtr_->registerElapsedTimer(elapsedTimerPtr_);

  mpdeMgrPtr_->registerTIAMPDEInterface(tiaMPDEIfacePtr_);

  if (!bsuccess)
  {
    string msg("N_ANP_AnalysisManager::setupMPDEMgr_ - Error!  A registration function failed!");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::initializeTransientModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::initializeTransientModel()
{
  wimPtr->createTimeIntegMethod(TIAMethod_BACKWARD_DIFFERENTIATION_15);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::evalTransientModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::evalTransientModel(
    double t,
    N_LAS_Vector * SolVectorPtr,
    N_LAS_Vector * CurrSolVectorPtr,
    N_LAS_Vector * LastSolVectorPtr,
    N_LAS_Vector * StaVectorPtr,
    N_LAS_Vector * CurrStaVectorPtr,
    N_LAS_Vector * LastStaVectorPtr,
    N_LAS_Vector * StaDerivVectorPtr,
    N_LAS_Vector * StoVectorPtr,
    N_LAS_Vector * CurrStoVectorPtr,
    N_LAS_Vector * LastStoVectorPtr,
    N_LAS_Vector * stoLeadCurrQCompVectorPtr,
    N_LAS_Vector * QVectorPtr,
    N_LAS_Vector * FVectorPtr,
    N_LAS_Vector * dFdxdVpVectorPtr,
    N_LAS_Vector * dQdxdVpVectorPtr,
    N_LAS_Matrix * dQdxMatrixPtr,
    N_LAS_Matrix * dFdxMatrixPtr
    )
{
  // This is F,Q load:
  bool bsuccess = loaderPtr->loadDAEVectors(
        SolVectorPtr,
        CurrSolVectorPtr,
        LastSolVectorPtr,
        StaVectorPtr,
        CurrStaVectorPtr,
        LastStaVectorPtr,
        StaDerivVectorPtr,
        StoVectorPtr,
        CurrStoVectorPtr,
        LastStoVectorPtr,
        stoLeadCurrQCompVectorPtr,
        QVectorPtr,
        FVectorPtr,
        dFdxdVpVectorPtr,
        dQdxdVpVectorPtr
        );
  // This is dQdx, dFdx load:
  bsuccess = bsuccess && loaderPtr->loadDAEMatrices(
      SolVectorPtr,
      StaVectorPtr,
      StaDerivVectorPtr,
      StoVectorPtr,
      dQdxMatrixPtr,
      dFdxMatrixPtr);
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::evalTransientModelState
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::evalTransientModelState(
    double t,
    N_LAS_Vector * SolVectorPtr,
    N_LAS_Vector * StaVectorPtr,
    N_LAS_Vector * StoVectorPtr
    )
{
  // This is part of state vector load:
  secPtr_->currentTime = t;
  secPtr_->nextTime = t;
  loaderPtr->updateSources(); // updates source values wrt time t.
  bool bsuccess = loaderPtr->updateState(
      SolVectorPtr,
      SolVectorPtr,
      SolVectorPtr,
      StaVectorPtr,
      StaVectorPtr,
      StaVectorPtr,
      StoVectorPtr,
      StoVectorPtr,
      StoVectorPtr
      );
  return bsuccess;
}

// ***** Accessor methods *****
//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setBeginningIntegrationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setBeginningIntegrationFlag(bool bif)
{
  primaryAnalysisObject_->setBeginningIntegrationFlag(bif);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getBeginningIntegrationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getBeginningIntegrationFlag()
{
  return primaryAnalysisObject_->getBeginningIntegrationFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setIntegrationMethod(int im)
{
  primaryAnalysisObject_->setIntegrationMethod (im);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
unsigned int N_ANP_AnalysisManager::getIntegrationMethod ()
{
  return primaryAnalysisObject_->getIntegrationMethod ();
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTotalLinearSolutionTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getTotalLinearSolutionTime() const
{
  return primaryAnalysisObject_->getTotalLinearSolutionTime();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTotalResidualLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getTotalResidualLoadTime() const
{
  return primaryAnalysisObject_->getTotalResidualLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTotalJacobianLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getTotalJacobianLoadTime() const
{
  return primaryAnalysisObject_->getTotalJacobianLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getDoubleDCOPEnabled ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getDoubleDCOPEnabled ()
{
  return primaryAnalysisObject_->getDoubleDCOPEnabled ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::simulationComplete
// Purpose       : return boolean signifying whether simulation complete or
//                 not.
//
// Special Notes :THIS VERSION IS ONLY VALID FOR TRANSIENT RUNS, where
//                 completion of the simulation means integration to final
//                 time.
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::simulationComplete()
{
  if (analysis == ANP_MODE_TRANSIENT)
  {
    return (secPtr_->finished());
  }
  else
  {
    string msg = "N_ANP_AnalysisManager::simulationComplete called for non-transient run... not currently valid!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    // This return here solely to shut the SGI compiler up
    return(true);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getPauseTime
//
// Purpose       : return the time at which the simulation will pause
//
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getPauseTime()
{
  return(tiaParams.pauseTime);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getStepNumber ()
{
  int number=0;
  if( !is_null( primaryAnalysisObject_ ) )
  {
    number = primaryAnalysisObject_->getStepNumber();
  }
  return number;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTranStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getTranStepNumber ()
{
  int number=0;
  if (analysis == ANP_MODE_TRANSIENT && !is_null(primaryAnalysisObject_) )
  {
    number = primaryAnalysisObject_->getTranStepNumber();
  }
  return number;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setStepNumber (int step)
{
  if( !is_null( primaryAnalysisObject_ ) )
  {
    primaryAnalysisObject_->setStepNumber(step);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setTranStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setTranStepNumber (int step)
{
  if( !is_null( primaryAnalysisObject_ ) )
  {
    primaryAnalysisObject_->setTranStepNumber(step);
  }
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getInitTranFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getInitTranFlag()
{
  int stepNum =  getStepNumber();
  return (stepNum <= 0);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTime
// Purpose       : Gets the next time value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getTime() const
{
  return secPtr_->nextTime;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getCurrentTime
// Purpose       : Gets the current time value.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 8/21/2009
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getCurrentTime() const
{
  return secPtr_->currentTime;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getFinalTime
// Purpose       : Gets the final time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getFinalTime() const
{
  return secPtr_->finalTime;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getInitialTime
// Purpose       : Gets the initial time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getInitialTime() const
{
  return secPtr_->initialTime;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getStartingTimeStep
// Purpose       : Gets the starting time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getStartingTimeStep ()
{
  return secPtr_->startingTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTimeIntMode
// Purpose       : Gets the time-integration method.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getTimeIntMode()
{
  return primaryAnalysisObject_->getIntegrationMethod();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getSteadyStateFlag
// Purpose       : Gets a flag indicating we are in steady state.
//                  (steady=true)
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getSteadyStateFlag ()
{
  return ((primaryAnalysisObject_->getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getInputOPFlag
// Purpose       : Gets a flag indicating we are starting from a previous OP
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/13/06
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getInputOPFlag ()
{
  return primaryAnalysisObject_->getInputOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTranOPFlag
// Purpose       : Gets a flag indicating we are in a DCOP calculation that is
//                 the initialization for a transient simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getTranOPFlag ()
{
  return ((analysis == ANP_MODE_TRANSIENT || primaryAnalysisObject_->isAnalysis(ANP_MODE_TRANSIENT))
     && (primaryAnalysisObject_->getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getACOPFlag
// Purpose       : Gets a flag indicating we are in a DCOP calculation that is
//                 the initialization for a transient simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/16/12
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getACOPFlag ()
{
  return ((analysis == ANP_MODE_AC || primaryAnalysisObject_->isAnalysis(ANP_MODE_AC))
     && (primaryAnalysisObject_->getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getDCSweepFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getDCSweepFlag ()
{
  return ((analysis == ANP_MODE_DC_SWEEP || primaryAnalysisObject_->isAnalysis(ANP_MODE_DC_SWEEP)));
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTransientFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getTransientFlag () const
{
  return (((analysis == ANP_MODE_TRANSIENT || primaryAnalysisObject_->isAnalysis(ANP_MODE_TRANSIENT))
     && (primaryAnalysisObject_->getIntegrationMethod()) !=TIAMethod_NONE));
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getDoubleDCOPStep
// Purpose       : Gets the double DC Operating Point step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/25/01
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getDoubleDCOPStep()
{
  return primaryAnalysisObject_->getDoubleDCOPStep();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getCurrentStepSize
// Purpose       : Returns the "current" time step size.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/15/01
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getCurrentStepSize()
{
  return secPtr_->currentTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getLastStepSize
// Purpose       : Returns the "last" time step size.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/15/01
//-----------------------------------------------------------------------------
double N_ANP_AnalysisManager::getLastStepSize()
{
  return secPtr_->lastTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::updateDerivs
// Purpose       : Calls the time  int. method to update the corrector
//                 derivatives.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::updateDerivs()
{
  wimPtr->obtainCorrectorDeriv();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::updateDivDiffs
// Purpose       : Updates the divided difference values.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::updateDivDiffs()
{
  dsPtr_->computeDividedDifferences();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::updateDerivsBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::updateDerivsBlock(
  const list<index_pair> & solGIDList,
  const list<index_pair> & staGIDList)
{
  dsPtr_->computeDivDiffsBlock(solGIDList,staGIDList);
  wimPtr->updateDerivsBlock (solGIDList, staGIDList);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::equateTmpVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::equateTmpVectors()
{
  return dsPtr_->equateTmpVectors();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerElapsedTimer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/24/06
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerElapsedTimer(N_UTL_Timer * et)
{
  elapsedTimerPtr_ = rcp(et,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::restartDataSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::restartDataSize( bool pack )
{
  return secPtr_->restartDataSize( pack );
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::dumpRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::dumpRestartData
  (char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  return secPtr_->dumpRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::restoreRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::restoreRestartData
  (char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  return secPtr_->restoreRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getSolnVarData( const int & gid,
                                             vector<double> & varData )
{
  return dsPtr_->getSolnVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getStateVarData( const int & gid,
                                              vector<double> & varData )
{
  return dsPtr_->getStateVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getStoreVarData( const int & gid,
                                              vector<double> & varData )
{
  return dsPtr_->getStoreVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getOrder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 02/15/07
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getOrder ()
{
  return wimPtr->getOrder();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getNumberOfSteps
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getNumberOfSteps ()
{
  return wimPtr->getNumberOfSteps();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getUsedOrder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getUsedOrder ()
{
  return wimPtr->getUsedOrder();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getNscsco
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getNscsco ()
{
  return wimPtr->getNscsco();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setSolnVarData( const int & gid,
                                             const vector<double> & varData )
{
  return dsPtr_->setSolnVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setStateVarData( const int & gid,
                                              const vector<double> & varData )
{
  return dsPtr_->setStateVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setStoreVarData( const int & gid,
                                              const vector<double> & varData )
{
  return dsPtr_->setStoreVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 09/06/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp)
{
  tiaParams = tiaParams_tmp;

  if (!Teuchos::is_null(secPtr_))
  {
    secPtr_->setTIAParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerMPDEInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/26/04
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerMPDEInterface( N_TIA_MPDEInterface * tiaMPDEIface_tmp)
{
  tiaMPDEIfacePtr_ = rcp(tiaMPDEIface_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerLinearSystem
// Purpose       : Added to make it easier for the lasSysPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerLinearSystem(N_LAS_System * lasSysPtr_tmp)
{
  lasSysPtr = rcp(lasSysPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerNLSManager
// Purpose       : Added to make it easier for the nlsSysPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerNLSManager(N_NLS_Manager * nlsMgrPtr_tmp)
{
  nlsMgrPtr = rcp(nlsMgrPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerLoader
// Purpose       : Added to make it easier for the loaderPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerLoader(N_LOA_Loader * loaderPtr_tmp)
{
  loaderPtr = rcp(loaderPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerOutputMgr
// Purpose       : Added to make it easier for the outputPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerOutputMgr(N_IO_OutputMgr * outputPtr_tmp)
{
  outMgrPtr = rcp(outputPtr_tmp,false);
  outputMgrAdapterRCPtr_ = rcp( new N_ANP_OutputMgrAdapter() );
  outputMgrAdapterRCPtr_->registerOutputMgr( outputPtr_tmp );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerRestartMgr
// Purpose       : Added to make it easier for the restartPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerRestartMgr(
  N_IO_RestartMgr * restartPtr_tmp)
{
  restartPtr = rcp(restartPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerNonlinearEquationLoader( N_LOA_NonlinearEquationLoader * nonlinearEquationLoaderPtr_tmp )
{
  nonlinearEquationLoaderPtr = rcp(nonlinearEquationLoaderPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerDeviceInterface ( N_DEV_DeviceInterface * devInterfacePtr_tmp )
{
  devInterfacePtr = rcp(devInterfacePtr_tmp,false);
  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerTopology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerTopology( N_TOP_Topology * topoMgrPtr_tmp )
{
  topoMgrPtr = rcp(topoMgrPtr_tmp,false);
  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerApplicationBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerApplicationBuilder( N_LAS_Builder * appBuilderPtr_tmp )
{
  appBuilderPtr = rcp(appBuilderPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::registerParallelServices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/19/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::registerParallelServices(N_PDS_Manager * pds_tmp)
{
  pdsMgrPtr = rcp(pds_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setNextSolVectorPtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 04/22/2003
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::setNextSolVectorPtr (N_LAS_Vector * solVecPtr)
{
  return dsPtr_->setNextSolVectorPtr (solVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setPauseTime
//
// Purpose       : Set the time at which to pause the simulation
//
// Special Notes : This is a temporary implementation that I'm using to
//                 begin the process of hiding time integrator internals from
//                 the "N_CIR_Xyce::simulateUntil" method, so that I can
//                 ultimately change those internals without the simulateUntil
//                 method
//                 In the zero order version, just sets tiaParams.pauseTime
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setPauseTime(double pauseTime)
{
  secPtr_->setBreakPoint(N_UTL_BreakPoint(pauseTime,PAUSE_BREAKPOINT));
}



//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::resumeSimulation
//
// Purpose       : set flag to signify that simulation is continuation of
//                 previously paused simulation.
//
// Special Notes : This is a temporary implementation that I'm using to
//                 begin the process of hiding time integrator internals from
//                 the "N_CIR_Xyce::simulateUntil" method, so that I can
//                 ultimately change those internals without the simulateUntil
//                 method
//                 In the zero order version, just sets tiaParams.resume=true
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::resumeSimulation()
{
  tiaParams.resume=true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::unset_resumeSimulation
// Purpose       : UN-set flag to signify that simulation is continuation of
//                 previously paused simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 06/12/2013
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::unset_resumeSimulation()
{
  tiaParams.resume=false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::outputIntervalSpecified_
// Purpose       : Return true if user has specified output control options
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems modeling
// Creation Date : 2/12/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::outputIntervalSpecified_()
{
  // This is necessary and sufficient ---
  return (initialOutputInterval_ > 0.0);
}

// routines to get/set Dakota run flags and actually run a Dakota iteration
//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getDakotaRunFlag()
{
  return dakotaRunFlag_;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setDakotaRunFlag( bool flag )
{
  dakotaRunFlag_ = flag;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
int N_ANP_AnalysisManager::getDakotaIteration()
{
  return dakotaIterationNumber_;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::setDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::setDakotaIteration( int iterNumber )
{
    dakotaIterationNumber_ = iterNumber;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getTimeIntInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::getTimeIntInfo (N_TIA_TimeIntInfo & tiInfo)
{
  tiInfo.currentOrder         = getOrder ();
  tiInfo.numberOfSteps        = getNumberOfSteps ();
  tiInfo.usedOrder            = getUsedOrder ();
  tiInfo.nscsco               = getNscsco ();

  tiInfo.pdt                  = partialTimeDerivative(); // alpha/DT
  tiInfo.nextTimeStep         = getCurrentStepSize();
  tiInfo.currTimeStep         = getLastStepSize();
  tiInfo.currentTime          = secPtr_->currentTime;
  tiInfo.nextTime             = getTime();
  tiInfo.finalTime            = getFinalTime();
  tiInfo.startingTimeStep     = getStartingTimeStep();
  tiInfo.bpTol                = getBreakpointTol();

  tiInfo.dcopFlag             = getSteadyStateFlag ();
  tiInfo.inputOPFlag          = getInputOPFlag ();
  tiInfo.tranopFlag           = getTranOPFlag ();
  tiInfo.acopFlag             = getACOPFlag ();
  tiInfo.transientFlag        = getTransientFlag ();
  tiInfo.dcsweepFlag          = getDCSweepFlag ();
  tiInfo.sweepSourceResetFlag = getSweepSourceResetFlag ();

  tiInfo.timeStepNumber       = getStepNumber();
  tiInfo.initTranFlag         = getInitTranFlag ();
  tiInfo.beginIntegrationFlag = getBeginningIntegrationFlag();

  tiInfo.doubleDCOPStep       = getDoubleDCOPStep();
  tiInfo.doubleDCOPEnabled    = getDoubleDCOPEnabled ();

  tiInfo.stepLoopIter         = 0;
  if( stepLoopFlag_ )
  {
    analysisObject_->getStepIter();
  }

  tiInfo.timeIntMode          = getTimeIntMode();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::silenceProgress
// Purpose       : Shut up "Percent Complete" noises
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::silenceProgress()
{
  progressFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::enableProgress
// Purpose       : stop shutting up "Percent Complete" noises
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::enableProgress()
{
  progressFlag_ = true;
}

