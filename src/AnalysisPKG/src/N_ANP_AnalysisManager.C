//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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
// Revision Number: $Revision: 1.113.2.1 $
// Revision Date  : $Date: 2014/03/04 23:50:53 $
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

#include <N_UTL_fwd.h>
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
#include <N_ANP_Report.h>

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

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::AnalysisManager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------

AnalysisManager::AnalysisManager(N_IO_CmdParse & cp, AnalysisInterface * anaIntPtr_tmp)
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
// Function      : AnalysisManager::~AnalysisManager
//
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
AnalysisManager::~AnalysisManager()
{
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::resetAll
//
// Purpose       : just like a destructor without the death
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
void AnalysisManager::resetAll()
{
  // tiaDataStore_ is created in initializeAll
  tiaDataStore_ = Teuchos::null;
  // secPtr_ is created in initializeAll
  secPtr_ = Teuchos::null;
  // wimPtr is created in initializeAll
  wimPtr = Teuchos::null;

  // xyceTranTimerPtr_ is created in run
  xyceTranTimerPtr_ = Teuchos::null;

  // tiaMPDEIfacePtr_'s copy to tiaDataStore_ and secPtr_ are reset in intializeAll, we don't need to delete it.

  // assemblerPtr is created in initializeAll
  assemblerPtr = Teuchos::null;

  // Reset step statistics to zero.
  primaryAnalysisObject_->resetAll ();

  initializeAllFlag_ = false;

}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBlockAnalysisFlag
// Purpose       : "get" function for MPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getBlockAnalysisFlag () const
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return blockAnalysisFlag_;
  }
  return mpdeMgrPtr_->blockAnalysisFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getMPDEFlag
// Purpose       : "get" function for MPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getMPDEFlag ()
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return mpdeFlag_;
  }
  return mpdeMgrPtr_->getMPDEFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getMPDEStartupFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getMPDEStartupFlag ()
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return false;
  }
  return mpdeMgrPtr_->getMPDEStartupFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getMPDEIcFlag
// Purpose       : get function for MPDE initial condition flag (true if MPDE & IC
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getMPDEIcFlag()
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return false;
  }
  return mpdeMgrPtr_->getMPDEIcFlag();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getWaMPDEFlag ()
// Purpose       : "get" function for WaMPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getWaMPDEFlag ()
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    return false;
  }
  return mpdeMgrPtr_->getWaMPDEFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::isPaused
//
// Purpose       : return the true if the simulation is currently paused
//
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical and Microsystems Simulation
// Creation Date : 02/18/2008
//-----------------------------------------------------------------------------
bool AnalysisManager::isPaused()
{
  return secPtr_->isPauseTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::initializeAll
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
bool AnalysisManager::initializeAll()
{
  std::string msg;

  if (Teuchos::is_null(lasSysPtr))
  {
    Report::DevelFatal0().in("AnalysisManager::initializeAll")
      << "Register LAS system first.";
  }

  tiaParams.solutionSize = lasSysPtr->getSolutionSize();
  tiaParams.stateSize    = lasSysPtr->getStateSize();

  // allocate data store class, which will allocate all the vectors.
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    Xyce::dout() << "AnalysisManager::initializeAll.  " ;
    Xyce::dout() << "  maximum order = " << tiaParams.maxOrder << std::endl;
  }
#endif

  tiaDataStore_ = rcp(new N_TIA_DataStore(&tiaParams, &*lasSysPtr));
  secPtr_ = rcp(new N_TIA_StepErrorControl(commandLine_,*this,tiaParams));

  // Now that data store has been created, we can also create the
  // working integration method object.
  wimPtr = rcp(new N_TIA_WorkingIntegrationMethod (tiaParams,*secPtr_,*tiaDataStore_));

  secPtr_->registerWIMPtr(&*wimPtr);

  assemblerPtr = rcp(new N_TIA_DAE_Assembler(*tiaDataStore_, *loaderPtr, *wimPtr, *pdsMgrPtr, daeStateDerivFlag_));

  tiaMPDEIfacePtr_->registerTIADataStore(&*tiaDataStore_);
  tiaMPDEIfacePtr_->registerTIAStepErrorControl(&*secPtr_);
  // The final time has to be forced to be a "paused" breakpoint,
  // after any registerTIAParams call.
  // RLS: May need to conditionally call this only when this is an MPDE run
  // not sure if we know that when this routine is called.
  //setPauseTime(tiaParams.finalTime);

  registerOutputIntervals();
  registerRestartIntervals();

  // register current:
  lasSysPtr->registerCurrStaVector(&(tiaDataStore_->currStatePtr));
  lasSysPtr->registerCurrSolVector(&(tiaDataStore_->currSolutionPtr));

  // register next:
  lasSysPtr->registerNextStaVector(&(tiaDataStore_->nextStatePtr));
  lasSysPtr->registerNextSolVector(&(tiaDataStore_->nextSolutionPtr));

  // register last:
  lasSysPtr->registerLastStaVector(&(tiaDataStore_->lastStatePtr));
  lasSysPtr->registerLastSolVector(&(tiaDataStore_->lastSolutionPtr));

  lasSysPtr->registerFlagSolVector(&(tiaDataStore_->flagSolutionPtr));

  // register temporaries:
  lasSysPtr->registerTmpSolVector(&(tiaDataStore_->tmpSolVectorPtr));
  lasSysPtr->registerTmpStaVector(&(tiaDataStore_->tmpStaVectorPtr));
  lasSysPtr->registerTmpStaDerivVector(&(tiaDataStore_->tmpStaDerivPtr));
  lasSysPtr->registerTmpStaDivDiffVector(&(tiaDataStore_->tmpStaDivDiffPtr));

  // register next derivatives:
  lasSysPtr->registerNextSolDerivVector(&(tiaDataStore_->nextSolutionDerivPtr));
  lasSysPtr->registerNextStaDerivVector(&(tiaDataStore_->nextStateDerivPtr));

  // register the device mask
  lasSysPtr->registerDeviceMaskVector(&(tiaDataStore_->deviceMaskPtr));

  // Get the RHS and the Jacobian
  tiaDataStore_->JMatrixPtr    = lasSysPtr->getJacobianMatrix();
  tiaDataStore_->RHSVectorPtr  = lasSysPtr->getRHSVector();

  // DAE formulation vectors
  lasSysPtr->registerDAEQVector     ( tiaDataStore_->daeQVectorPtr );
  lasSysPtr->registerDAEFVector     ( tiaDataStore_->daeFVectorPtr );

  // DAE formulation matrices
  lasSysPtr->registerDAEdQdxMatrix  ( tiaDataStore_->dQdxMatrixPtr );
  lasSysPtr->registerDAEdFdxMatrix  ( tiaDataStore_->dFdxMatrixPtr );

  // Get the RHS and the Jacobian
  //tiaDataStore_->JMatrixPtr    = lasSysPtr->getJacobianMatrix();
  //tiaDataStore_->RHSVectorPtr  = lasSysPtr->getRHSVector();

  tiaDataStore_->dFdxdVpVectorPtr = lasSysPtr->getdFdxdVpVector ();
  tiaDataStore_->dQdxdVpVectorPtr = lasSysPtr->getdQdxdVpVector ();

  tiaDataStore_->limiterFlag = loaderPtr->getLimiterFlag ();

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
      Report::UserError0() << "No analysis statement in the netlist";
      return false;
    }
  }

  // Allocate analysis objects, and also set up params.
  allocateAnalysisObject_ ();

  initializeAllFlag_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::run
// Purpose       : Execute the control loop for the set analysis type.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
bool AnalysisManager::run()
{
  std::string msg;

  if (initializeAllFlag_ == false)
  {
    Report::DevelFatal0().in("AnalysisManager::run")
      << "Call the initializeAll function first";
  }

  if (!analysisParamsRegistered)
  {
    Report::UserError0() << "No analysis statement in the netlist";
    return false;
  }

  // now that the problem is set up, we can tell the datastore what the var type are
  // this info can be used by the time integrator to control error in a way that is
  // appropriate for different variable types (like V's and I's).
  std::vector<char> varTypes;
  topoMgrPtr->returnVarTypeVec( varTypes );
  tiaDataStore_->varTypeVec = varTypes;

  // This checks to make sure that all quantities in the .print are valid before
  // we continue.
  outputMgrAdapterRCPtr_->check_output( analysis );

#ifdef Xyce_VERBOSE_TIME
  tiaParams.printParams(Xyce::lout(), anpAnalysisModeToNLS(analysis));
#endif

#ifndef Xyce_NO_MASKED_WRMS_NORMS
  if (loaderPtr->loadDeviceMask())
  {
     //Xyce::dout() << " Nontrivial mask! " << std::endl;
  }
  else
  {
    // Xyce::dout() << " Trivial mask! " << std::endl;
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
  // the DCSweep or Step classes are allocated.  They are now primarily allocated
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

  Report::safeBarrier(pdsMgrPtr->getPDSComm()->comm());
  
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
// Function      : AnalysisManager::allocateAnalysisObject_
// Purpose       : Allocate analysis objects, and also setup params.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/10
//-----------------------------------------------------------------------------
void AnalysisManager::allocateAnalysisObject_ ()
{
  std::string msg("");

  if( !tiaParams.resume )
  {
    if (analysis == ANP_MODE_TRANSIENT )
    {
      analysisObject_ = Teuchos::rcp(new Transient(this));
      analysisObject_->setAnalysisParams(tranParamsBlock);
      secPtr_->resetAll();
    }
    else if (analysis == ANP_MODE_DC_SWEEP)
    {
      analysisObject_ = Teuchos::rcp(new DCSweep(this));
      for (int i=0;i<(int)dcParamsBlockVec.size();++i)
      {
        analysisObject_->setAnalysisParams(dcParamsBlockVec[i]);
      }
    }
    else if (analysis == ANP_MODE_MPDE)
    {
      analysisObject_ = Teuchos::rcp(new MPDE(this));
      analysisObject_->setAnalysisParams(mpdeParamsBlock);
    }
    else if (analysis == ANP_MODE_HB)
    {
      analysisObject_ = Teuchos::rcp(new HB(this));
      analysisObject_->setAnalysisParams(hbParamsBlock);
      Teuchos::rcp_dynamic_cast<HB>(analysisObject_)->setHBOptions(hbOptionsBlock);
      Teuchos::rcp_dynamic_cast<HB>(analysisObject_)->setHBLinSol(hbLinSolBlock);
      Teuchos::rcp_dynamic_cast<HB>(analysisObject_)->setLinSol(linSolBlock);
      setHBFlag( true );
    }
    else if ( analysis == ANP_MODE_AC)
    {
      analysisObject_ = Teuchos::rcp(new AC(this));
      analysisObject_->setAnalysisParams(acParamsBlock);
    }
    else if ( analysis == ANP_MODE_MOR)
    {
      analysisObject_ = Teuchos::rcp(new MOR(this));
      analysisObject_->setAnalysisParams(morParamsBlock);
    }
    else
    {
      Report::UserError0() <<  "Unknown type of analysis";
      return;
    }

    if( !is_null(analysisObject_) )
    {
      analysisObject_->setParamsWithOutputMgrAdapter ( outputMgrAdapterRCPtr_ );
    }

    // need to throw an error here as this really shouldn't happen.
    if( is_null( analysisObject_ ) )
    {
      Report::DevelFatal0().in("AnalysisManager::initializeAll")
        << "Unable to allocate analysis type";
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
      analysisObject_ = Teuchos::rcp(new Step(this, stepAnalysisTarget_.get()));
      for (int i=0;i<(int)stepParamsBlockVec.size();++i)
      {
        analysisObject_->setAnalysisParams(stepParamsBlockVec[i]);
      }
      analysisObject_->setParamsWithOutputMgrAdapter ( outputMgrAdapterRCPtr_ );
    }

    if (dakotaRunFlag_ && is_null(dakotaAnalysisTarget_) )
    {
      dakotaAnalysisTarget_ = analysisObject_;
      analysisObject_ = Teuchos::rcp(new Dakota(this, dakotaAnalysisTarget_.get()));
      analysisObject_->setAnalysisParams(dakotaParamsBlock);
      analysisObject_->setParamsWithOutputMgrAdapter ( outputMgrAdapterRCPtr_ );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::printLoopInfo
// Purpose       : Prints out time loop information.
// Special Notes : Prints stats from save point start to save point finish.
//                 Special case 0,0 is entire run to this point
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
bool AnalysisManager::printLoopInfo(int start, int finish)
{
  return primaryAnalysisObject_->printLoopInfo (start,finish);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::testDCOPOutputTime_
// Purpose       : Similar to testOutputTime, except that this is for
//                 DCOP restart files.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisManager::testDCOPOutputTime_()
{
  bool flag(true);

  if( !dcopRestartFlag_ )
  {
    flag = false;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::testSaveOutputTime_
// Purpose       : Similar to testOutputTime, except that this is for
//                 .SAVE files.
// Special Notes : Only outputs 1x.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisManager::testSaveOutputTime_()
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
    Xyce::dout() <<"Calling SAVE outputs!" <<std::endl;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::testRestartSaveTime_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel ComputationalSciences.
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::testRestartSaveTime_()
{
  bool flag;

#ifdef Xyce_DEBUG_RESTART
  Xyce::dout() << "TESTING FOR RESTART SAVE" << std::endl
               << Xyce::subsection_divider << std::endl
               << "secPtr_->currentTime: " << secPtr_->currentTime << std::endl
               << "nextSaveTime: " << nextRestartSaveTime_ << std::endl
               << "initialRestartInterval_: " << initialRestartInterval_ << std::endl;
  if (!(restartIntervals_.empty()))
  {
    Xyce::dout() << "First restart interval: " << restartIntervals_[0].first << std::endl;
  }
  else
  {
    Xyce::dout() << "restartIntervals_ is empty" << std::endl;
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
    std::pair<double, double> currInterval, nextInterval;
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
  Xyce::dout() << "new nextSaveTime: " << nextRestartSaveTime_ << std::endl
               << "restart flag: " << flag << std::endl
               << Xyce::subsection_divider << std::endl;
#endif

  return flag;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::partialTimeDerivative
// Purpose       : Returns the current partial time derivative for either the
//                 solution or state vector.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
double AnalysisManager::partialTimeDerivative()
{
  double partDT = wimPtr->partialTimeDeriv();

  // Add a check here to try to prevent capacitive spiral of death.
  // This is still a "research" option, so is not on by default.
  if (tiaParams.jacLimitFlag)
  {
#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      Xyce::dout() << "AnalysisManager::partialTimeDerivative.";
      Xyce::dout() << "   Using jac limit = " << tiaParams.jacLimit << std::endl;
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
// Function      : AnalysisManager::getBreakpointTol
// Purpose       : Returns the breakpoint tolerance.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/7/02
//-----------------------------------------------------------------------------
double AnalysisManager::getBreakpointTol()
{
  return N_UTL_BreakPoint::getBPTol();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setBreakpointTol
// Purpose       : Sets the breakpoint tolerance.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/7/02
//-----------------------------------------------------------------------------
void AnalysisManager::setBreakpointTol(double bptol)
{
  N_UTL_BreakPoint::setBPTol(bptol);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerRestartIntervals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::registerRestartIntervals()
{
  double initInt;
  std::vector< std::pair<double,double> > intPairs;
  restartPtr->getRestartIntervals(initInt, intPairs);
  initialRestartInterval_ = initInt;
  nextRestartSaveTime_    = secPtr_->initialTime;
  restartIntervals_       = intPairs;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerOutputIntevals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::registerOutputIntervals()
{
  double initInt;
  std::vector< std::pair<double,double> > intPairs;
  outputMgrAdapterRCPtr_->getOutputIntervals( initInt, & intPairs );
  initialOutputInterval_ = initInt;
  nextOutputTime_ = secPtr_->initialTime;
  outputIntervals_ = intPairs;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTranAnalysisParams
// Purpose       : Sets transient analysis parameters (from .TRAN)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool AnalysisManager::setTranAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysis = ANP_MODE_TRANSIENT;
  tranParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDCAnalysisParams
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
bool AnalysisManager::setDCAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysis = ANP_MODE_DC_SWEEP;
  // before saving a copy of paramsBlock, check and see if it already
  // is in the dcParamsBlockVec.  If so then replace the original copy.
  // This replacement is necessary if the first copy had an expression
  // element that was resolved in later parsing.

  bool foundMatch = false;
  std::vector<N_UTL_OptionBlock>::iterator paramsBlockVecItr = dcParamsBlockVec.begin();
  std::vector<N_UTL_OptionBlock>::iterator paramsBlockVecEnd = dcParamsBlockVec.end();
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
// Function      : AnalysisManager::setOPAnalysisParams
// Purpose       : Handle OP statement. (.OP)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
bool AnalysisManager::setOPAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  dotOpSpecified_ = true;
  opParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSTEPAnalysisParams
// Purpose       : Sets the STEP calculation parameters. (from .STEP)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/03
//-----------------------------------------------------------------------------
bool AnalysisManager::setSTEPAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  stepLoopFlag_ = true;
  // before saving a copy of paramsBlock, check and see if it already
  // is in the stepParamsBlockVec.  If so then replace the original copy.
  // This replacement is necessary if the first copy had an expression
  // element that was resolved in later parsing.
  bool foundMatch = false;
  std::vector<N_UTL_OptionBlock>::iterator paramsBlockVecItr = stepParamsBlockVec.begin();
  std::vector<N_UTL_OptionBlock>::iterator paramsBlockVecEnd = stepParamsBlockVec.end();
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
// Function      : AnalysisManager::setSaveOptions
// Purpose       : Sets the Save parameters.
// Special Notes : Most of these parameters are handled in the output manager,
//                 rather than here.  So, most params are a noop here, except
//                 for "TIME".
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisManager::setSaveOptions(
  const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    Xyce::dout() << "In AnalysisManager::setSaveOptions" << std::endl;
  }
#endif

  saveFlag_ = true;

  std::list<N_UTL_Param>::const_iterator iterPL = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
    Xyce::dout() << "iterPL->tag = " << iterPL->tag() << std::endl;
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
      saveTime_ = iterPL->getImmutableValue<double>();
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
// Function      : AnalysisManager::setACAnalysisParams
// Purpose       : Sets the AC sweep calculation parameters (from .AC)
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::setACAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysis = ANP_MODE_AC;
  acParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setMORAnalysisParams
// Purpose       : Sets the MOR calculation parameters (from .MOR)
//
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 05/29/12
//-----------------------------------------------------------------------------
bool AnalysisManager::setMORAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysis = ANP_MODE_MOR;
  morParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setMOROptions
// Purpose       :
// Special Notes : These are from '.options mor'
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 05/31/12
//-----------------------------------------------------------------------------
bool AnalysisManager::setMOROptions(const N_UTL_OptionBlock & OB)
{
  std::list<N_UTL_Param>::const_iterator it_tpL;
  std::list<N_UTL_Param>::const_iterator first = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag()=="METHOD")
    {
      ExtendedString stringVal ( it_tpL->stringValue() );
      stringVal.toUpper();
      tiaParams.morMethod = stringVal;
    }
    else if (it_tpL->uTag()=="SAVEREDSYS")
    {
      tiaParams.morSaveRedSys = static_cast<bool>(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="COMPORIGTF")
    {
      tiaParams.morCompOrigTF = static_cast<bool>(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="COMPREDTF")
    {
      tiaParams.morCompRedTF = static_cast<bool>(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="COMPTYPE")
    {
      ExtendedString stringVal ( it_tpL->stringValue() );
      stringVal.toUpper();
      tiaParams.morCompType = stringVal;
    }
    else if (it_tpL->uTag()=="COMPNP")
    {
      tiaParams.morCompNP = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="COMPFSTART")
    {
      tiaParams.morCompFStart = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="COMPFSTOP")
    {
      tiaParams.morCompFStop = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="EXPPOINT")
    {
      tiaParams.morExpPoint = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="SCALETYPE")
    {
      tiaParams.morScaleType = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="SCALEFACTOR")
    {
      tiaParams.morScaleFactor = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="SCALEFACTOR1")
    {
      tiaParams.morScaleFactor1 = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="SPARSIFICATIONTYPE")
    {
      tiaParams.morSparsificationType = it_tpL->getImmutableValue<int>();
    }
    else
    {
      Report::UserError0() << it_tpL->uTag() << " is not a recognized model-order reduction option.";
    }
  }

  // If we are computing the transfer function, make sure the frequency range is valid.
  if (tiaParams.morCompOrigTF || tiaParams.morCompRedTF)
  {
    if (tiaParams.morCompFStop < tiaParams.morCompFStart)
    {
      Report::UserError() << ".options mor COMPFSTART = " << tiaParams.morCompFStart << " > " << tiaParams.morCompFStop << " = COMPFSTOP!";
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDCOPRestartParams
// Purpose       : Sets the DCOP restart parameters.
// Special Notes : Most of the dcop restart parameters are used and handled by
//                 the N_IO_OutputMgr class, so this function here doens't need
//                 to do very much.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisManager::setDCOPRestartParams(const N_UTL_OptionBlock & OB)
{
  dcopRestartFlag_ = true;

  std::list<N_UTL_Param>::const_iterator iterPL = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
    Xyce::dout() << "iterPL->tag = " << iterPL->tag() << std::endl;
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
      saveTime_ = iterPL->getImmutableValue<double>();
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
// Function      : AnalysisManager::setTranOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::setTranOptions(const N_UTL_OptionBlock & OB)
{
  std::list<N_UTL_Param>::const_iterator it_tpL;
  std::list<N_UTL_Param>::const_iterator first = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag()=="METHOD")
    {
//      tiaParams.integrationMethod=it_tpL->iVal();

      if (it_tpL->isInteger())
        tiaParams.integrationMethod=it_tpL->getImmutableValue<int>();
      else
      {

        ExtendedString stringVal ( it_tpL->stringValue() );
        stringVal.toUpper();

        if (stringVal == "TRAP" || stringVal == "TRAPEZOIDAL")
          tiaParams.integrationMethod = 7;
        else if (stringVal == "BDF")
          tiaParams.integrationMethod = 6;
        else if (stringVal == "GEAR")
          tiaParams.integrationMethod = 8;
        else
        {
          Report::UserError0() << "Unsupported transient method type";
        }
      }

    }
#ifdef Xyce_DEBUG_ANALYSIS
    else if (it_tpL->uTag()=="CONSTSTEP")
    {
      tiaParams.constantStepSize
        =  static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
#endif
    else if (it_tpL->uTag()=="USEDEVICEMAX")
    {
      tiaParams.useDeviceTimeStepMax =  static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="RELTOL")
    {
      tiaParams.relErrorTol=it_tpL->getImmutableValue<double>();
      tiaParams.relErrorTolGiven = true;
    }
    else if (it_tpL->uTag()=="ABSTOL")
    {
      tiaParams.absErrorTol=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="DOUBLEDCOPSTEP")
    {
      tiaParams.doubleDCOPStep=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="FIRSTDCOPSTEP")
    {
      tiaParams.firstDCOPStep=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="LASTDCOPSTEP")
    {
      tiaParams.lastDCOPStep=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="BPENABLE" )
    {
      tiaParams.bpEnable =  static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="RESTARTSTEPSCALE" )
    {
      tiaParams.restartTimeStepScale=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="EXITTIME" )
    {
      tiaParams.exitTime=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="EXITSTEP" )
    {
      tiaParams.exitStep=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag() == "MINTIMESTEPSBP")
    {
      tiaParams.minTimeStepsBP = it_tpL->getImmutableValue<int>();
      tiaParams.minTimeStepsBPGiven = true;
    }
    else if (it_tpL->uTag()=="ERROPTION" )
    {
      tiaParams.errorAnalysisOption=it_tpL->getImmutableValue<int>();
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
      tiaParams.nlNearConvFlag = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="NLSMALLUPDATE" )
    {
      tiaParams.nlSmallUpdateFlag = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="JACLIMITFLAG" )
    {
      tiaParams.jacLimitFlag= static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="JACLIMIT" )
    {
      tiaParams.jacLimit = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="DAESTATEDERIV" )
    {
      daeStateDerivFlag_ = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="MAXORD" )
    {
      tiaParams.maxOrder = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="MINORD" )
    {
      tiaParams.minOrder = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="TIMESTEPSREVERSAL" )
    {
      tiaParams.timestepsReversal =it_tpL->getImmutableValue<bool>();
    }
    else if (it_tpL->uTag()=="TESTFIRSTSTEP" )
    {
      tiaParams.testFirstStep = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="DELMAX" )
    {
      tiaParams.delmax =it_tpL->getImmutableValue<double>();
      tiaParams.delmaxGiven = true;
//      tiaParams.delmax =it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="NLMIN" )
    {
      tiaParams.NLmin=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="NLMAX" )
    {
      tiaParams.NLmax=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="OUTPUTINTERPMPDE")
    {
      tiaParams.outputInterpMPDE = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="NEWLTE")
    {
      tiaParams.newLte =it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="NEWBPSTEPPING")
    {
      tiaParams.newBPStepping = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="INTERPOUTPUT")
    {
      tiaParams.interpOutputFlag = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="CONDTEST")
    {
      tiaParams.condTestFlag = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="CONDTESTDEVICENAME")
    {
      tiaParams.condTestDeviceNames.push_back(it_tpL->stringValue() );
    }
    else if (it_tpL->uTag() == "DTMIN")
    {
      tiaParams.userSpecMinTimeStep = it_tpL->getImmutableValue<double>();
      tiaParams.userSpecMinTimeStepGiven = true;
    }
    else if (it_tpL->uTag() == "PASSNLSTALL")
    {
      tiaParams.passNLStall = it_tpL->getImmutableValue<bool>();
    }
    else if (it_tpL->uTag() == "FASTTESTS")
    {
      tiaParams.fastTests = it_tpL->getImmutableValue<bool>();
    }
    else if (it_tpL->uTag() == "MINTIMESTEPRECOVERY")
    {
      tiaParams.minTimeStepRecoveryCounter = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="VOLTZEROTOL" )
    {
      tiaParams.voltZeroTol=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="CURRZEROTOL" )
    {
      tiaParams.currZeroTol=it_tpL->getImmutableValue<double>();
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
        tiaParams.debugLevel = it_tpL->getImmutableValue<int>();
      }
#endif
    }
    else if (it_tpL->uTag()=="HISTORYTRACKINGDEPTH" )
    {
      tiaParams.historyTrackingDepth = it_tpL->getImmutableValue<int>();
    }
    else
    {
      Report::UserError() << it_tpL->uTag() << " is not a recognized time integration option";
    }
  }

  if (tiaParams.NLmin > tiaParams.NLmax)
  {
    Report::UserError() << ".options timeint NLMIN = " << tiaParams.NLmin << " > " << tiaParams.NLmax << " = NLMAX!";
  }

  if (tiaParams.firstDCOPStep < 0) tiaParams.firstDCOPStep = 0;
  if (tiaParams.firstDCOPStep > 1) tiaParams.firstDCOPStep = 1;
  if (tiaParams.lastDCOPStep < 0)  tiaParams.lastDCOPStep  = 0;
  if (tiaParams.lastDCOPStep > 1)  tiaParams.lastDCOPStep  = 1;

  // check for maximum order on the command line
  if ( commandLine_.argExists ("-maxord") )
  {

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

  return true;

}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setMPDEAnalysisParams
// Purpose       : Sets the MPDE sweep calculation parameters (from .MPDE)
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setMPDEAnalysisParams(const N_UTL_OptionBlock & OB)
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
// Function      : AnalysisManager::setMPDEOptions
// Purpose       :
// Special Notes : from '.options mpdeint'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setMPDEOptions(const N_UTL_OptionBlock & OB)
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    setupMPDEMgr_();
  }
  mpdeMgrPtr_->setMPDEOptions(OB);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setHBAnalysisParams
// Purpose       : Sets the HB sweep calculation parameters (from .HB)
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/30/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setHBAnalysisParams(const N_UTL_OptionBlock & OB)
{
  std::list<N_UTL_Param>::const_iterator it_tpL;
  std::list<N_UTL_Param>::const_iterator first = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag() == "FREQ")
    {

      tiaParams.freqs = it_tpL->getValue<std::vector<double> >();
      tiaParams.freqGiven = true;
    }
  }

  if (tiaParams.freqs[0] <= 0.0 )
  {
    Report::UserError() << "Frequency of oscillation " << tiaParams.freqs[0] << " is less than or equal to zero, invalid .HB specification";
  }

  if ((DEBUG_ANALYSIS & tiaParams.debugLevel) > 0)
  {
    dout() << section_divider << std::endl
           << "HB transient simulation parameters" 
           //<< Util::push << std::endl
           << std::endl
           << "HB frequency = " << tiaParams.freqs[0] << std::endl
           //<< Util::pop << std::endl;
           << std::endl;
  }

  setBlockAnalysisFlag( true );
  devInterfacePtr->setBlockAnalysisFlag( true );
  analysis = ANP_MODE_HB;
  analysisParamsRegistered = true;
  hbParamsBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setHBOptions
// Purpose       :
// Special Notes : from '.options hb'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setHBOptions(const N_UTL_OptionBlock & OB)
{
  hbOptionsBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------
bool AnalysisManager::setLinSol(const N_UTL_OptionBlock & OB)
{
  linSolBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setHBLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setHBLinSol(const N_UTL_OptionBlock & OB)
{
  hbLinSolBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTRANMPDEOptions
// Purpose       :
// Special Notes : from '.options timeint-mpde'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setTRANMPDEOptions(const N_UTL_OptionBlock & OB)
{
  if (Teuchos::is_null(mpdeMgrPtr_))
  {
    setupMPDEMgr_();
  }
  mpdeMgrPtr_->registerTranMPDEOptions(OB);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSensOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 6/5/13
//-----------------------------------------------------------------------------
bool AnalysisManager::setSensOptions(const N_UTL_OptionBlock & OB)
{
  sensFlag_=true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::completeOPStartStep
// Purpose       : Call to rotate next state to current state following a
//               : constrained DCOP solve when using a previous operating point
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 06/27/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::completeOPStartStep  ( )
{
  bool bsuccess = true;

  bsuccess = tiaDataStore_->updateStateDataArrays ();
  tiaDataStore_->setConstantHistory();
  tiaDataStore_->equateTmpVectors();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::completeHomotopyStep
// Purpose       : Call to do final cleanup after a single,
//                 successful homotopy step.
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::completeHomotopyStep
    ( const std::vector<std::string> & paramNames,
      const std::vector<double> & paramVals,
      N_LAS_Vector * solnVecPtr )
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  std::string netListFile = commandLine_.getArgumentValue("netlist");
  Xyce::dout() << "\n " << netListFile;
  Xyce::dout() << " AnalysisManager::completeHomotopyStep " << std::endl;
#endif

  // Rotate the data vectors:
  bool bs1 = tiaDataStore_->updateStateDataArrays ();    bsuccess = bsuccess && bs1;
  tiaDataStore_->setConstantHistory();
  tiaDataStore_->equateTmpVectors();

  // Pass info in to the lower level solver:
  loaderPtr->homotopyStepSuccess (paramNames,paramVals);

  // Call output
  outputMgrAdapterRCPtr_->outputHomotopy( paramNames, paramVals, *solnVecPtr );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::failHomotopyStep
// Purpose       : Call to do final cleanup after a single,
//                 failed homotopy step.
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::failHomotopyStep ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  std::string netListFile = commandLine_.getArgumentValue("netlist");
  Xyce::dout() << "\n " << netListFile;
  Xyce::dout() << " AnalysisManager::failHomotopyStep " << std::endl;
#endif

  // Pass info in to the lower level solver:
  loaderPtr->homotopyStepFailure ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setupMPDEMgr_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
void AnalysisManager::setupMPDEMgr_()
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
    Report::DevelFatal0().in("AnalysisManager::setupMPDEMgr_") << "Registration function failed";
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::initializeTransientModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
void AnalysisManager::initializeTransientModel()
{
  wimPtr->createTimeIntegMethod(TIAMethod_BACKWARD_DIFFERENTIATION_15);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::evalTransientModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//-----------------------------------------------------------------------------
bool AnalysisManager::evalTransientModel(
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
// Function      : AnalysisManager::evalTransientModelState
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
bool AnalysisManager::evalTransientModelState(
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
// Function      : AnalysisManager::setBeginningIntegrationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setBeginningIntegrationFlag(bool bif)
{
  primaryAnalysisObject_->setBeginningIntegrationFlag(bif);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBeginningIntegrationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getBeginningIntegrationFlag()
{
  return primaryAnalysisObject_->getBeginningIntegrationFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setIntegrationMethod(int im)
{
  primaryAnalysisObject_->setIntegrationMethod (im);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
unsigned int AnalysisManager::getIntegrationMethod ()
{
  return primaryAnalysisObject_->getIntegrationMethod ();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalLinearSolutionTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalLinearSolutionTime() const
{
  return primaryAnalysisObject_->getTotalLinearSolutionTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalResidualLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalResidualLoadTime() const
{
  return primaryAnalysisObject_->getTotalResidualLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalJacobianLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalJacobianLoadTime() const
{
  return primaryAnalysisObject_->getTotalJacobianLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDoubleDCOPEnabled ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getDoubleDCOPEnabled ()
{
  return primaryAnalysisObject_->getDoubleDCOPEnabled ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::simulationComplete
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
bool AnalysisManager::simulationComplete()
{
  if (analysis == ANP_MODE_TRANSIENT)
  {
    return (secPtr_->finished());
  }
  else
  {
    Report::DevelFatal0().in("AnalysisManager::simulationComplete") << "Called for non-transient run, not currently valid";
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getPauseTime
//
// Purpose       : return the time at which the simulation will pause
//
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
double AnalysisManager::getPauseTime()
{
  return(tiaParams.pauseTime);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
int AnalysisManager::getStepNumber ()
{
  int number=0;
  if( !is_null( primaryAnalysisObject_ ) )
  {
    number = primaryAnalysisObject_->getStepNumber();
  }
  return number;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTranStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
int AnalysisManager::getTranStepNumber ()
{
  int number=0;
  if (analysis == ANP_MODE_TRANSIENT && !is_null(primaryAnalysisObject_) )
  {
    number = primaryAnalysisObject_->getTranStepNumber();
  }
  return number;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
void AnalysisManager::setStepNumber (int step)
{
  if( !is_null( primaryAnalysisObject_ ) )
  {
    primaryAnalysisObject_->setStepNumber(step);
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTranStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
void AnalysisManager::setTranStepNumber (int step)
{
  if( !is_null( primaryAnalysisObject_ ) )
  {
    primaryAnalysisObject_->setTranStepNumber(step);
  }
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInitTranFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
bool AnalysisManager::getInitTranFlag()
{
  int stepNum =  getStepNumber();
  return (stepNum <= 0);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTime
// Purpose       : Gets the next time value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getTime() const
{
  return secPtr_->nextTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getCurrentTime
// Purpose       : Gets the current time value.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 8/21/2009
//-----------------------------------------------------------------------------
double AnalysisManager::getCurrentTime() const
{
  return secPtr_->currentTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getFinalTime
// Purpose       : Gets the final time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getFinalTime() const
{
  return secPtr_->finalTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInitialTime
// Purpose       : Gets the initial time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getInitialTime() const
{
  return secPtr_->initialTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStartingTimeStep
// Purpose       : Gets the starting time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getStartingTimeStep ()
{
  return secPtr_->startingTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTimeIntMode
// Purpose       : Gets the time-integration method.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
int AnalysisManager::getTimeIntMode()
{
  return primaryAnalysisObject_->getIntegrationMethod();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getSteadyStateFlag
// Purpose       : Gets a flag indicating we are in steady state.
//                  (steady=true)
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getSteadyStateFlag ()
{
  return ((primaryAnalysisObject_->getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInputOPFlag
// Purpose       : Gets a flag indicating we are starting from a previous OP
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/13/06
//-----------------------------------------------------------------------------
bool AnalysisManager::getInputOPFlag ()
{
  return primaryAnalysisObject_->getInputOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTranOPFlag
// Purpose       : Gets a flag indicating we are in a DCOP calculation that is
//                 the initialization for a transient simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getTranOPFlag ()
{
  return ((analysis == ANP_MODE_TRANSIENT || primaryAnalysisObject_->isAnalysis(ANP_MODE_TRANSIENT))
     && (primaryAnalysisObject_->getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getACOPFlag
// Purpose       : Gets a flag indicating we are in a DCOP calculation that is
//                 the initialization for a transient simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/16/12
//-----------------------------------------------------------------------------
bool AnalysisManager::getACOPFlag ()
{
  return ((analysis == ANP_MODE_AC || primaryAnalysisObject_->isAnalysis(ANP_MODE_AC))
     && (primaryAnalysisObject_->getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDCSweepFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getDCSweepFlag ()
{
  return ((analysis == ANP_MODE_DC_SWEEP || primaryAnalysisObject_->isAnalysis(ANP_MODE_DC_SWEEP)));
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTransientFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getTransientFlag () const
{
  return (((analysis == ANP_MODE_TRANSIENT || primaryAnalysisObject_->isAnalysis(ANP_MODE_TRANSIENT))
     && (primaryAnalysisObject_->getIntegrationMethod()) !=TIAMethod_NONE));
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDoubleDCOPStep
// Purpose       : Gets the double DC Operating Point step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/25/01
//-----------------------------------------------------------------------------
int AnalysisManager::getDoubleDCOPStep()
{
  return primaryAnalysisObject_->getDoubleDCOPStep();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getCurrentStepSize
// Purpose       : Returns the "current" time step size.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/15/01
//-----------------------------------------------------------------------------
double AnalysisManager::getCurrentStepSize()
{
  return secPtr_->currentTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getLastStepSize
// Purpose       : Returns the "last" time step size.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/15/01
//-----------------------------------------------------------------------------
double AnalysisManager::getLastStepSize()
{
  return secPtr_->lastTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::updateDerivs
// Purpose       : Calls the time  int. method to update the corrector
//                 derivatives.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool AnalysisManager::updateDerivs()
{
  wimPtr->obtainCorrectorDeriv();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::updateDivDiffs
// Purpose       : Updates the divided difference values.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
bool AnalysisManager::updateDivDiffs()
{
  tiaDataStore_->computeDividedDifferences();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::updateDerivsBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
bool AnalysisManager::updateDerivsBlock(
  const std::list<index_pair> & solGIDList,
  const std::list<index_pair> & staGIDList)
{
  tiaDataStore_->computeDivDiffsBlock(solGIDList,staGIDList);
  wimPtr->updateDerivsBlock (solGIDList, staGIDList);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::equateTmpVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool AnalysisManager::equateTmpVectors()
{
  return tiaDataStore_->equateTmpVectors();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerElapsedTimer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/24/06
//-----------------------------------------------------------------------------
bool AnalysisManager::registerElapsedTimer(N_UTL_Timer * et)
{
  elapsedTimerPtr_ = rcp(et,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::restartDataSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
int AnalysisManager::restartDataSize( bool pack )
{
  return secPtr_->restartDataSize( pack );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::dumpRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::dumpRestartData
  (char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  return secPtr_->dumpRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::restoreRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::restoreRestartData
  (char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  return secPtr_->restoreRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getSolnVarData( const int & gid,
                                             std::vector<double> & varData )
{
  return tiaDataStore_->getSolnVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getStateVarData( const int & gid,
                                              std::vector<double> & varData )
{
  return tiaDataStore_->getStateVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getStoreVarData( const int & gid,
                                              std::vector<double> & varData )
{
  return tiaDataStore_->getStoreVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getOrder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 02/15/07
//-----------------------------------------------------------------------------
int AnalysisManager::getOrder ()
{
  return wimPtr->getOrder();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNumberOfSteps
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int AnalysisManager::getNumberOfSteps ()
{
  return wimPtr->getNumberOfSteps();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getUsedOrder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int AnalysisManager::getUsedOrder ()
{
  return wimPtr->getUsedOrder();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNscsco
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int AnalysisManager::getNscsco ()
{
  return wimPtr->getNscsco();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::setSolnVarData( const int & gid,
                                             const std::vector<double> & varData )
{
  return tiaDataStore_->setSolnVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::setStateVarData( const int & gid,
                                              const std::vector<double> & varData )
{
  return tiaDataStore_->setStateVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::setStoreVarData( const int & gid,
                                              const std::vector<double> & varData )
{
  return tiaDataStore_->setStoreVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 09/06/01
//-----------------------------------------------------------------------------
bool AnalysisManager::registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp)
{
  tiaParams = tiaParams_tmp;

  if (!Teuchos::is_null(secPtr_))
  {
    secPtr_->setTIAParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerMPDEInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/26/04
//-----------------------------------------------------------------------------
bool AnalysisManager::registerMPDEInterface( N_TIA_MPDEInterface * tiaMPDEIface_tmp)
{
  tiaMPDEIfacePtr_ = rcp(tiaMPDEIface_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerLinearSystem
// Purpose       : Added to make it easier for the lasSysPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerLinearSystem(N_LAS_System * lasSysPtr_tmp)
{
  lasSysPtr = rcp(lasSysPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerNLSManager
// Purpose       : Added to make it easier for the nlsSysPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerNLSManager(N_NLS_Manager * nlsMgrPtr_tmp)
{
  nlsMgrPtr = rcp(nlsMgrPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerLoader
// Purpose       : Added to make it easier for the loaderPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerLoader(N_LOA_Loader * loaderPtr_tmp)
{
  loaderPtr = rcp(loaderPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerOutputMgr
// Purpose       : Added to make it easier for the outputPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerOutputMgr(N_IO_OutputMgr * outputPtr_tmp)
{
  outMgrPtr = rcp(outputPtr_tmp,false);
  outputMgrAdapterRCPtr_ = rcp( new OutputMgrAdapter() );
  outputMgrAdapterRCPtr_->registerOutputMgr( outputPtr_tmp );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerRestartMgr
// Purpose       : Added to make it easier for the restartPtr pointer to be
//                 non-static.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerRestartMgr(
  N_IO_RestartMgr * restartPtr_tmp)
{
  restartPtr = rcp(restartPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::registerNonlinearEquationLoader( N_LOA_NonlinearEquationLoader * nonlinearEquationLoaderPtr_tmp )
{
  nonlinearEquationLoaderPtr = rcp(nonlinearEquationLoaderPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::registerDeviceInterface ( N_DEV_DeviceInterface * devInterfacePtr_tmp )
{
  devInterfacePtr = rcp(devInterfacePtr_tmp,false);
  return true;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerTopology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::registerTopology( N_TOP_Topology * topoMgrPtr_tmp )
{
  topoMgrPtr = rcp(topoMgrPtr_tmp,false);
  return true;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerApplicationBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::registerApplicationBuilder( N_LAS_Builder * appBuilderPtr_tmp )
{
  appBuilderPtr = rcp(appBuilderPtr_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerParallelServices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/19/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerParallelServices(N_PDS_Manager * pds_tmp)
{
  pdsMgrPtr = rcp(pds_tmp,false);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setNextSolVectorPtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 04/22/2003
//-----------------------------------------------------------------------------
bool AnalysisManager::setNextSolVectorPtr (N_LAS_Vector * solVecPtr)
{
  return tiaDataStore_->setNextSolVectorPtr (solVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setPauseTime
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
void AnalysisManager::setPauseTime(double pauseTime)
{
  secPtr_->setBreakPoint(N_UTL_BreakPoint(pauseTime, Util::PAUSE_BREAKPOINT));
}



//-----------------------------------------------------------------------------
// Function      : AnalysisManager::resumeSimulation
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
void AnalysisManager::resumeSimulation()
{
  tiaParams.resume=true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::unset_resumeSimulation
// Purpose       : UN-set flag to signify that simulation is continuation of
//                 previously paused simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 06/12/2013
//-----------------------------------------------------------------------------
void AnalysisManager::unset_resumeSimulation()
{
  tiaParams.resume=false;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::outputIntervalSpecified_
// Purpose       : Return true if user has specified output control options
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems modeling
// Creation Date : 2/12/07
//-----------------------------------------------------------------------------
bool AnalysisManager::outputIntervalSpecified_()
{
  // This is necessary and sufficient ---
  return (initialOutputInterval_ > 0.0);
}

// routines to get/set Dakota run flags and actually run a Dakota iteration
//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getDakotaRunFlag()
{
  return dakotaRunFlag_;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setDakotaRunFlag( bool flag )
{
  dakotaRunFlag_ = flag;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
int AnalysisManager::getDakotaIteration()
{
  return dakotaIterationNumber_;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setDakotaIteration( int iterNumber )
{
    dakotaIterationNumber_ = iterNumber;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTimeIntInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::getTimeIntInfo (N_TIA_TimeIntInfo & tiInfo)
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
// Function      : AnalysisManager::silenceProgress
// Purpose       : Shut up "Percent Complete" noises
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void AnalysisManager::silenceProgress()
{
  progressFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::enableProgress
// Purpose       : stop shutting up "Percent Complete" noises
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void AnalysisManager::enableProgress()
{
  progressFlag_ = true;
}

} // namespace Analysis
} // namespace Xyce
