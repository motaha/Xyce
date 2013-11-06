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
// Filename      : $RCSfile: N_ANP_Transient.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.86.2.4 $
// Revision Date  : $Date: 2013/10/03 17:23:31 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <sstream>
#include <iomanip>


// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_TIA_Assembler.h>
#include <N_TIA_DataStore.h>
#include <N_LOA_Loader.h>
#include <N_LAS_System.h>
#include <N_MPDE_Manager.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_RestartMgr.h>
#include <N_TIA_StepErrorControl.h>
#include <N_UTL_Timer.h>
#include <N_UTL_NoCase.h>

#include <N_IO_CmdParse.h>

#include<N_ANP_Transient.h>

#include<N_UTL_ExpressionData.h>
#include<N_NLS_ReturnCodes.h>

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::N_ANP_Transient( N_ANP_AnalysisManager * )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
N_ANP_Transient::N_ANP_Transient( N_ANP_AnalysisManager * anaManagerPtr ) :
  N_ANP_AnalysisBase(anaManagerPtr),
  initialIntegrationMethod_(TIAMethod_ONESTEP),
  firstTranOutput_(true),
  isPaused(false),
  dcopFlag_(true),
  startDCOPtime(0.0),
  endTRANtime(0.0),
  gui_(false),
  historyTrackingOn_(true),
  maxTimeStepExpressionGiven_(false),
  maxTimeStepExpressionAsString_(""),
  firstTime(true),
  oldPercentComplete(0.0),
  startSimTime(-1.0),
  dcStats(0),
  tranStats(0)
{
  gui_ = commandLine_.argExists("-gui");
  restartMgrRCPtr_ = anaManagerPtr->restartPtr;

  // initialize this to an empty vector
  outputInterpolationTimes_.clear();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .TRAN statement.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/21/10
//-----------------------------------------------------------------------------
bool N_ANP_Transient::setAnalysisParams(const N_UTL_OptionBlock & paramsBlock)
{
  list<N_UTL_Param>::const_iterator it_tp;
  list<N_UTL_Param>::const_iterator first = paramsBlock.getParams().begin();
  list<N_UTL_Param>::const_iterator last = paramsBlock.getParams().end();
  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag()      == "TSTART")
    {
      tiaParams.tStart = it_tp->dVal();
      tiaParams.tStartGiven = true;
    }
    else if (it_tp->uTag() == "TSTOP")
    {
      tiaParams.finalTime   = it_tp->dVal();
    }
    else if (it_tp->uTag() == "TSTEP")
    {
      tiaParams.userSpecified_startingTimeStep = it_tp->dVal();
    }
    else if (it_tp->uTag() == "NOOP" ||
             it_tp->uTag() == "UIC")
    {
      tiaParams.NOOP = 1;
    }
    else if (it_tp->uTag() == "DTMAX")
    {
      tiaParams.maxTimeStep = it_tp->dVal();
      tiaParams.maxTimeStepGiven = true;
#ifdef Xyce_DEBUG_ANALYSIS
      if (tiaParams.debugLevel > 0)
      {
	      std::cout << "setting maxTimeStep = " << tiaParams.maxTimeStep << std::endl;
      }
#endif
    }
    else if (it_tp->uTag() == "MAXTIMEEXPRESSION")
    {
      // an expression was given which should specify what
      // max time step to use based on the current simulation
      // time.  The expected format of this is
      // {schedule( t0, max_del_t0, t1, max_del_t1, ...)}
      // we'll just store the expression as a string for now
      // as we need to make sure other classes are up and running
      // before we can create the actual expression.  That will
      // happen in N_ANP_Transient class as it will ultimately use this
      maxTimeStepExpressionGiven_ = true;
      maxTimeStepExpressionAsString_ = it_tp->sVal();
    }
  }

  if (tiaParams.finalTime <= tiaParams.tStart || tiaParams.finalTime <= 0 || tiaParams.tStart < 0)
  {
    ostringstream ost;
    ost << " In N_TIA_StepErrorControl::setTranAnalysisParams: " << std::endl;
    ost << "Final time of " << tiaParams.finalTime
        << " is earlier or same as start time of "
        << tiaParams.tStart << std::endl;
    ost << " Check netlist for invalid .TRAN specification " << std::endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, ost.str());
  }

  // if starting time steps is baloney, set then to a default value.
  // ERK.  This trap is redudant with an identical one in step error
  // control.
  if ( tiaParams.userSpecified_startingTimeStep <= 0.0 )
  {
    tiaParams.userSpecified_startingTimeStep = 1.0e-10;
  }

#ifdef Xyce_DEBUG_ANALYSIS
  string dashedline = "-------------------------------------------------"
                       "----------------------------\n";

  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  Transient simulation parameters");

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  initial time = ",
            tiaParams.initialTime);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  final time   = ",
            tiaParams.finalTime);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  starting time step = ",
            tiaParams.userSpecified_startingTimeStep);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  tStart (time of first output) = ",
            tiaParams.tStart);

    if (!(tiaParams.NOOP))
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
         "  NOOP/UIC is NOT set");
    }
    else
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  NOOP/UIC is set");
    }

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       dashedline);
  }
#endif

  if( maxTimeStepExpressionGiven_ )
  {
    // set-up expression object for getting a user specified, time dependent
    // max time step value
    maxTimeStepExpressionGiven_ = true;
	  maxTimeStepExpressionRCPtr_ = rcp( new N_UTL_ExpressionData( maxTimeStepExpressionAsString_,  outputMgrAdapterRCPtr_->getOutputMgrPtr() ) );
  }

  // queueSize_ should can be set from the .options timeint line.  Use default in tiaParams
  // as that will be overwritten if the user set it via .options.  Also, if it's less than
  // or equal to zero assume that the user wants this turned off.
  // We'll store the attempted time step, time step status,
  // estimated error over tolerance, non-linear solver status and non-linear solver norm
  // for each step, up to the last queueSize_ steps.

  if (tiaParams.historyTrackingDepth > 0 )
  {
    historyTrackingOn_ = true;
    queueSize_= tiaParams.historyTrackingDepth;
    timeQueue_.set_size( queueSize_ );
    timeStepQueue_.set_size( queueSize_ );
    stepStatusQueue_.set_size( queueSize_ );
    estErrorOverTolQueue_.set_size( queueSize_ );
    nonlinearSolverStatusQueue_.set_size( queueSize_ );
    nonlinearSolverNumIterationsQueue_.set_size( queueSize_ );
    nonlinearSolverMaxNormQueue_.set_size( queueSize_ );
    nonlinearSolverMaxNormIndexQueue_.set_size( queueSize_ );
  }
  else
  {
    historyTrackingOn_ = false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::run()
{
  bool bsuccess = true;
  isPaused = false;

  if (anaManagerRCPtr_->getHBFlag() && !anaManagerRCPtr_->getMPDEStartupFlag())
    resetForHB();  
    
  bsuccess = bsuccess & init();
  bsuccess = bsuccess & loopProcess();

  // if processing the loop failed,
  // then skip finish step
  if( bsuccess )
  {
    bsuccess = bsuccess & finish();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::init()
{
  bool bsuccess = true;

  initialIntegrationMethod_ = TIAMethod_ONESTEP;
  if (tiaParams.integrationMethod != TIAMethod_ONESTEP)
  {
    initialIntegrationMethod_ = tiaParams.integrationMethod;
  }

  if (!tiaParams.resume)
  {
    dcopFlag_ = !tiaParams.NOOP;
    if(dcopFlag_)
    {
      anaManagerRCPtr_->currentMode_ = 0;
    }
    else
    {
      anaManagerRCPtr_->currentMode_ = 1;
    }

    doubleDCOPFlag_ = loaderRCPtr_->getDoubleDCOPFlag();
    doubleDCOPStep_ = tiaParams.firstDCOPStep;

    if (restartMgrRCPtr_->isRestart())
    {
      loaderRCPtr_->getInitialQnorm(dsRCPtr_->innerErrorInfoVec);

//      secRCPtr_->initialTime = secRCPtr_->currentTime;
      wimRCPtr_->createTimeIntegMethod(integrationMethod_);
      wimRCPtr_->initialize();

      dcopFlag_ = false;
      anaManagerRCPtr_->currentMode_ = 1;

#ifdef Xyce_PARALLEL_MPI
      // Update vectors with off proc values.
      lasSystemRCPtr_->updateExternValsSolnVector(dsRCPtr_->nextSolutionPtr);
      lasSystemRCPtr_->updateExternValsSolnVector(dsRCPtr_->currSolutionPtr);
      lasSystemRCPtr_->updateExternValsStateVector(dsRCPtr_->nextStatePtr);
      lasSystemRCPtr_->updateExternValsStateVector(dsRCPtr_->currStatePtr);
      lasSystemRCPtr_->updateExternValsStoreVector(dsRCPtr_->nextStorePtr);
      lasSystemRCPtr_->updateExternValsStoreVector(dsRCPtr_->currStorePtr);
#endif

      // Set the nonlinear solver parameters to those appropriate for the
      // transient solution.  If starting from here, we're not doing a DCOP
      // calculation first - we're jumping right into the transient.
      nlsMgrRCPtr_->setAnalysisMode(anpAnalysisModeToNLS(ANP_MODE_TRANSIENT));
    }
    else
    {
      if (dcopFlag_)
      {
        // Get set to do the operating point.
        integrationMethod_ = TIAMethod_NONE;
        wimRCPtr_->createTimeIntegMethod(integrationMethod_);
      }
      else // otherwise NOOP/UIC
      {
        integrationMethod_ = initialIntegrationMethod_;
        wimRCPtr_->createTimeIntegMethod(integrationMethod_);
      }

      // This setInitialGuess call is to up an initial guess in the
      // devices that have them (usually PDE devices).  This is DIFFERENT
      // than an initial condition.
      loaderRCPtr_->setInitialGuess (dsRCPtr_->nextSolutionPtr);
     
      if (anaManagerRCPtr_->getBlockAnalysisFlag())
      {
        beginningIntegration = true;
      }

      // 12/13/07 tscoffe/tmei:  .ic should not be applied to MPDE currently.
      // A new mechanism must be set up for this, e.g. .mpde_ic
      if (!anaManagerRCPtr_->getMPDEFlag() && !anaManagerRCPtr_->getMPDEStartupFlag())
      {
        // If available, set initial solution.  This may be from a DCOP restart,
        // a .IC, or a .NODESET line.  These can be used in both the UIC/NOOP
        // case, as well as the DCOP case, so they need to be outside the if-statement.
        inputOPFlag_ = outputMgrAdapterRCPtr_->setupInitialConditions
          ( *(dsRCPtr_->nextSolutionPtr),*(dsRCPtr_->flagSolutionPtr));
      }

      if (!dcopFlag_ && !anaManagerRCPtr_->getMPDEFlag())
      {
        // this "initializeProblem" call is to set the IC's in devices that
        // have them.  This is done IN PLACE of the operating point.
        loaderRCPtr_->initializeProblem
              ((dsRCPtr_->nextSolutionPtr),
               (dsRCPtr_->currSolutionPtr),
               (dsRCPtr_->lastSolutionPtr),
               (dsRCPtr_->nextStatePtr),
               (dsRCPtr_->currStatePtr),
               (dsRCPtr_->lastStatePtr),
               (dsRCPtr_->nextStateDerivPtr),
               (dsRCPtr_->nextStorePtr),
               (dsRCPtr_->currStorePtr),
               (dsRCPtr_->lastStorePtr),
               (dsRCPtr_->daeQVectorPtr),
               (dsRCPtr_->daeFVectorPtr),
               (dsRCPtr_->dFdxdVpVectorPtr),
               (dsRCPtr_->dQdxdVpVectorPtr) );

        // Do this to populate the q-vector:
        // since we're also skipping the DC OP, this call will
        // also force the loader to propagate the state vector to the
        // device manager which updates state vector seen by the devices.
        assemblerRCPtr_->loadRHS();
      }

      // Set a constant history.
      dsRCPtr_->setConstantHistory();
      dsRCPtr_->computeDividedDifferences();
      wimRCPtr_->obtainCorrectorDeriv();

#ifdef Xyce_PARALLEL_MPI
      // Update vectors with off proc values.
      lasSystemRCPtr_->updateExternValsSolnVector(dsRCPtr_->nextSolutionPtr);
      lasSystemRCPtr_->updateExternValsSolnVector(dsRCPtr_->currSolutionPtr);
      lasSystemRCPtr_->updateExternValsStateVector(dsRCPtr_->nextStatePtr);
      lasSystemRCPtr_->updateExternValsStateVector(dsRCPtr_->currStatePtr);
      lasSystemRCPtr_->updateExternValsStoreVector(dsRCPtr_->nextStorePtr);
      lasSystemRCPtr_->updateExternValsStoreVector(dsRCPtr_->currStorePtr);
#endif

      if (!dcopFlag_ && !tiaParams.resume)
      {
        // Now do a NOOP output.  ERK: 08/30/2007  This isn't the best place to
        // put this, but it will have to do for now.  If this isn't here, then
        // NOOP/UIC simulations don't output at t=0.0 in the *prn file.
        noopOutputs ();
      }

      stepNumber            = 0;
      tranStepNumber        = 0;
      anaManagerRCPtr_->breakPointRestartStep = 0;
    }
  }
  else
  {
    if(anaManagerRCPtr_->currentMode_ == 1)
    {
      // we're resuming from a previous transient simulation
      // this is an ugly side effect of the analysis manager's deleting
      // the transient object after it retured in a paused state
      // need to fix this in Xyce 6.0 .  RLS 12/22/2009
      dcopFlag_=false;

      // we default firstTranOutput_ to true to force the first time point to be output
      // however, if we're resuming, then we need this to be false
      // this is really a defect of a paused simulation's N_ANP_Transient object
      // being deleted by the AnalysisManager and then a new one being created
      // when the simulation is resumed.  Need to fully fix this.  RLS 12/22/2009
      firstTranOutput_ = false;
    }

    // in a 2 level problem, we can be asked to resume
    // without first making the integration method.  So
    // catch that.
    if (wimRCPtr_->getIntegMethodPtr() == NULL)
    {
      integrationMethod_ = tiaParams.integrationMethod;
      wimRCPtr_->createTimeIntegMethod(integrationMethod_);
    }

#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
      std::cout << "  transient loop called with resume true " << std::endl;
#endif
  }

  anaManagerRCPtr_->switchIntegrator_ = false;

  if (!tiaParams.resume)
  {
    // if we're resuming, this was already done prior to pausing, and doing
    // it again screws us up
    double suggestedMaxTime=0.0;
    if( maxTimeStepExpressionGiven_ )
    {
      suggestedMaxTime = maxTimeStepExpressionRCPtr_->evaluate(
        dsRCPtr_->currSolutionPtr, dsRCPtr_->currStatePtr, dsRCPtr_->currStorePtr);
    }
    secRCPtr_->updateMaxTimeStep( suggestedMaxTime );
    secRCPtr_->updateMinTimeStep();
    secRCPtr_->updateBreakPoints();
  }

  // reset min error tracking variables
  // if the step number is less than zero we'll assume there is no valid
  // min. estimated error over tol or an associated time step.  This frees us
  // from putting Machine::Big here and trying to do something with that.
  stepNumberAtMinEstErrorOverTol = -1;
  minEstErrorOverTol = 0.0;
  timeStepAtMinEstErrorOverTol = 0.0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::loopProcess()
// Purpose       : Conduct the time stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::loopProcess()
{
  bool bsuccess = true;

  // Transient time stepping loop:
  while (!(secRCPtr_->finished()))
  {
#ifdef Xyce_VERBOSE_TIME
    this->printStepHeader();
#endif
    this->printProgress();

    // ------------------------------------------------------------------------
    // If the flag is set to switch integration methods, do that here.
    // For example, switch from operating point to transient backward euler.

    if (anaManagerRCPtr_->switchIntegrator_)
    {
      wimRCPtr_->createTimeIntegMethod(integrationMethod_);
    }

    // ------------------------------------------------------------------------
    // Set the step size, current time and next time.

    secRCPtr_->updateStopTime();

#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      std::cout << std::endl;
      std::cout << "N_ANP_Transient::loopProcess()" << std::endl;
      std::cout << "beginningIntegration = " << beginningIntegration << std::endl;
      std::cout << "secRCPtr_->stepAttemptStatus = " << secRCPtr_->stepAttemptStatus << std::endl;
    }
#endif

    if (beginningIntegration &&
        secRCPtr_->stepAttemptStatus)
    {
      // ------------------------------------------------------------------------
      // 07/29/04 TSC:  initial step-size selection is now done in the
      // integration method initialize call under the DAE formulation.  This
      // segregates codes changes better and makes sense to do as an
      // initialization step even if its changed later.
      loaderRCPtr_->getInitialQnorm(dsRCPtr_->innerErrorInfoVec);
      double suggestedMaxTime=0.0;
      if( maxTimeStepExpressionGiven_ )
      {
        suggestedMaxTime = maxTimeStepExpressionRCPtr_->evaluate(
          dsRCPtr_->currSolutionPtr, dsRCPtr_->currStatePtr, dsRCPtr_->currStorePtr);
      }
      secRCPtr_->updateMaxTimeStep( suggestedMaxTime );
      wimRCPtr_->initialize();
    }

    // ------------------------------------------------------------------------
    // If we've switched the integration method, we need to obtain the
    // corrector derivative only after we've updated the TimeInfo.
    if (anaManagerRCPtr_->switchIntegrator_)
    {
      anaManagerRCPtr_->switchIntegrator_ = false;
      wimRCPtr_->obtainCorrectorDeriv();
    }

#ifdef Xyce_VERBOSE_TIME
    if (!dcopFlag_)
      secRCPtr_->outputTimeInfo();
#endif

    // ------------------------------------------------------------------------
    // Set the nonlinear solver parameters to those appropriate for the
    // transient solution, if neccessary.
    if (!dcopFlag_)
    {
      nlsMgrRCPtr_->setAnalysisMode(anpAnalysisModeToNLS(ANP_MODE_TRANSIENT));
    }

    // Ask the method to update its coefficients
    wimRCPtr_->updateCoeffs();

    // ------------------------------------------------------------------------
    // Perform the time step:
    takeAnIntegrationStep_();

    // ------------------------------------------------------------------------

    if (dcopFlag_ && secRCPtr_->stepAttemptStatus)
    {
      processSuccessfulDCOP();
    }
    else if (dcopFlag_ && !secRCPtr_->stepAttemptStatus)
    {
      processFailedDCOP();
      bsuccess = false;
      break;
    }

    // Transient
    else
    {
      if (secRCPtr_->stepAttemptStatus)
      {
        processSuccessfulStep();
      }
      else if( (tiaParams.passNLStall == true) &&
               !(secRCPtr_->stepAttemptStatus)  &&
               (secRCPtr_->currentTimeStep < (4*secRCPtr_->minTimeStep)) )
      {
        // potentially a VERY dangerous options.
        // if the non-linear solver is stalling, and we're very close to a min
        // time step, then calls this failure a pass
        if( secRCPtr_->newtonConvergenceStatus == -3)
        {
          string msg = "N_ANP_Transient::loopProcess() Nonlinear solver stalled. Calling this a pass";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
          processSuccessfulStep();
        }
        // another VERY dangerous options.
        // if the non-linear solver is reporting too big of an update, and we're very close to a min
        // time step, then calls this failure a pass
        if( secRCPtr_->newtonConvergenceStatus == -2)
        {
          string msg = "N_ANP_Transient::loopProcess() Update too big. Calling this a pass";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
          processSuccessfulStep();
        }
        // another VERY dangerous options.
        // if the non-linear solver is not converging in the max number of steps,
        // and we're very close to a min time step, then calls this failure a pass
        /*
        if( secRCPtr_->newtonConvergenceStatus == -1)
        {
          string msg = "N_ANP_Transient::loopProcess() Too many steps. Calling this a pass";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
          processSuccessfulStep();
        }
        */
        else
        {
          // process this failed step as we would have by default.
          bool b1 = processFailedStep();
          if (!b1)
          {
            bsuccess = false;
            break;
          }
        }
      }
      else // stepAttemptStatus  (ie do this if the step FAILED)
      {
        bool b1 = processFailedStep();
        if (!b1)
        {
          bsuccess = false;
          break;
        }
      } // stepAttemptStatus

    } // transient

#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      std::cout << std::endl;
      std::cout << "   Here we are, just before checking whether to pause. " << std::endl;
      std::cout << "   minTimeStep = " << secRCPtr_->minTimeStep << std::endl;
      std::cout << "   final time = " << tiaParams.finalTime << std::endl;
      std::cout << "   pause time = " << anaManagerRCPtr_->getPauseTime() << std::endl;
      std::cout << "   initial time = " << secRCPtr_->initialTime << std::endl;
      std::cout << "   current time = " << secRCPtr_->currentTime << std::endl;
      if (anaManagerRCPtr_->getPauseTime() == secRCPtr_->currentTime)
      {
        std::cout << "    Pause time and current time equal " << std::endl;
      }
      else
      {
        std::cout << "     difference between current and pause times is " 
          << anaManagerRCPtr_->getPauseTime() - secRCPtr_->currentTime << std::endl;
      }
      if (anaManagerRCPtr_->getPauseTime() == secRCPtr_->initialTime)
      {
        std::cout << "    Pause time and initial time equal " << std::endl;
      }
    }
#endif

    if (secRCPtr_->isPauseTime())
    {
      // Failure at this point only indicates that the simulation
      // is paused and may be resumed.
#ifdef Xyce_DEBUG_ANALYSIS
      if (tiaParams.debugLevel > 0)
      {
        std::cout << "N_ANP_Transient::loopProcess():   pausing simulation " << std::endl;
      }
#endif
      secRCPtr_->simulationPaused();
      isPaused = true;
      bsuccess = true;
      break;
    }

    // If the exit time has been exceeded, exit.
    if (tiaParams.exitTime != 0.0 &&
        secRCPtr_->currentTime > tiaParams.exitTime)
    {
        string msg = "N_ANP_Transient::loopProcess(): "
          "  Exit Time exceeded.  Exiting transient loop\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);
        bsuccess = true;
      break;
    }

    if (tiaParams.exitStep != -1 &&
          static_cast<int>(stepNumber) == tiaParams.exitStep)
    {
        string msg = "N_ANP_Transient::loopProcess(): "
          "  Exit Step.  Exiting transient loop\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);
        bsuccess = true;
      break;
    }

  } // end of time loop
  endTRANtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
  tranStats = saveLoopInfo();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::mixedSignalStep
// Purpose       :
// Special Notes : Habanero API function
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
bool N_ANP_Transient::mixedSignalStep()
{
  takeAnIntegrationStep_();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::preStepDetails
// Purpose       :
// Special Notes : Habanero API function
// Scope         : private_
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
void N_ANP_Transient::preStepDetails (double maxTimeStepFromHabanero)
{
#ifdef Xyce_VERBOSE_TIME
  this->printStepHeader();
#endif
  this->printProgress();

  // ------------------------------------------------------------------------
  // If the flag is set to switch integration methods, do that here.
  // For example, switch from operating point to transient backward euler.

  if (anaManagerRCPtr_->switchIntegrator_)
  {
    wimRCPtr_->createTimeIntegMethod(integrationMethod_);
  }

  // ------------------------------------------------------------------------
  // Set the step size, current time and next time.
  secRCPtr_->updateStopTime();

  // ------------------------------------------------------------------------
  // If a max time is set from Habanero, impose it now that stopTime is updated.
  if (maxTimeStepFromHabanero > 0)
  {
    double currentTimeStep = Xycemin(maxTimeStepFromHabanero, secRCPtr_->currentTimeStep);
    secRCPtr_->setTimeStep(currentTimeStep);
  }

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    std::cout << std::endl;
    std::cout << "N_ANP_Transient::loopProcess()" << std::endl;
    std::cout << "beginningIntegration = " << beginningIntegration << std::endl;
    std::cout << "secRCPtr_->stepAttemptStatus = " << secRCPtr_->stepAttemptStatus << std::endl;
  }
#endif


  if (beginningIntegration &&
      secRCPtr_->stepAttemptStatus)
  {
    // ------------------------------------------------------------------------
    // 07/29/04 TSC:  initial step-size selection is now done in the
    // integration method initialize call under the DAE formulation.  This
    // segregates codes changes better and makes sense to do as an
    // initialization step even if its changed later.
    loaderRCPtr_->getInitialQnorm(dsRCPtr_->innerErrorInfoVec);
    double suggestedMaxTime=0.0;
    if( maxTimeStepExpressionGiven_ )
    {
      suggestedMaxTime = maxTimeStepExpressionRCPtr_->evaluate(
        dsRCPtr_->currSolutionPtr, dsRCPtr_->currStatePtr, dsRCPtr_->currStorePtr);
    }
    secRCPtr_->updateMaxTimeStep( suggestedMaxTime );
    wimRCPtr_->initialize();
  }

  // ------------------------------------------------------------------------
  // If we've switched the integration method, we need to obtain the
  // corrector derivative only after we've updated the TimeInfo.
  if (anaManagerRCPtr_->switchIntegrator_)
  {
    anaManagerRCPtr_->switchIntegrator_ = false;
    wimRCPtr_->obtainCorrectorDeriv();
  }

#ifdef Xyce_VERBOSE_TIME
  if (!dcopFlag_)
  {
    secRCPtr_->outputTimeInfo();
  }
#endif

  // ------------------------------------------------------------------------
  // Set the nonlinear solver parameters to those appropriate for the
  // transient solution, if neccessary.
  if (!dcopFlag_)
  {
    nlsMgrRCPtr_->setAnalysisMode(anpAnalysisModeToNLS(ANP_MODE_TRANSIENT));
  }

  // Ask the method to update its coefficients
  wimRCPtr_->updateCoeffs();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::finalizeStep
// Purpose       :
// Special Notes : Habanero API function
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
bool N_ANP_Transient::finalizeStep ()
{
  bool recoverableFailureFlag = true;

  if (dcopFlag_ && secRCPtr_->stepAttemptStatus)
  {
    processSuccessfulDCOP();
  }
  else if (dcopFlag_ && !secRCPtr_->stepAttemptStatus)
  {
    processFailedDCOP();
    recoverableFailureFlag = false;
  }
  // Transient
  else
  {
    if (secRCPtr_->stepAttemptStatus)
    {
      processSuccessfulStep();
    }
    else if( (tiaParams.passNLStall == true) &&
              !(secRCPtr_->stepAttemptStatus)  &&
              (secRCPtr_->currentTimeStep < (4*secRCPtr_->minTimeStep)) )
    {
      // potentially a VERY dangerous options.
      // if the non-linear solver is stalling, and we're very close to a min
      // time step, then calls this failure a pass

      if( secRCPtr_->newtonConvergenceStatus == -3)
      {
        string msg = "N_ANP_Transient::loopProcess() Nonlinear solver stalled. Calling this a pass";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
        processSuccessfulStep();
      }
      // another VERY dangerous options.
      // if the non-linear solver is reporting too big of an update, and we're very close to a min
      // time step, then calls this failure a pass
      if( secRCPtr_->newtonConvergenceStatus == -2)
      {
        string msg = "N_ANP_Transient::loopProcess() Update too big. Calling this a pass";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
        processSuccessfulStep();
      }
      // another VERY dangerous options.
      // if the non-linear solver is not converging in the max number of steps,
      // and we're very close to a min time step, then calls this failure a pass
      /*
      if( secRCPtr_->newtonConvergenceStatus == -1)
      {
        string msg = "N_ANP_Transient::loopProcess() Too many steps. Calling this a pass";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
        processSuccessfulStep();
      }
      */
      else
      {
        // process this failed step as we would have by default.
        recoverableFailureFlag = processFailedStep();
      }
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      recoverableFailureFlag = processFailedStep();
    }
  } // transient

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    std::cout << std::endl;
    std::cout << "   Here we are, just before checking whether to pause. " << std::endl;
    std::cout << "   minTimeStep = " << secRCPtr_->minTimeStep << std::endl;
    std::cout << "   final time = " << tiaParams.finalTime << std::endl;
    std::cout << "   pause time = " << anaManagerRCPtr_->getPauseTime() << std::endl;
    std::cout << "   initial time = " << secRCPtr_->initialTime << std::endl;
    std::cout << "   current time = " << secRCPtr_->currentTime << std::endl;
    if (anaManagerRCPtr_->getPauseTime() == secRCPtr_->currentTime)
    {
      std::cout << "    Pause time and current time equal " << std::endl;
    }
    else
    {
      std::cout << "     difference between current and pause times is " << anaManagerRCPtr_->getPauseTime() - secRCPtr_->currentTime << std::endl;
    }
    if (anaManagerRCPtr_->getPauseTime() == secRCPtr_->initialTime)
    {
      std::cout << "    Pause time and initial time equal " << std::endl;
    }
  }
#endif

  if (secRCPtr_->isPauseTime())
  {
    // Failure at this point only indicates that the simulation
    // is paused and may be resumed.
#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      std::cout << "N_ANP_Transient::loopProcess():   pausing simulation " << std::endl;
    }
#endif
    secRCPtr_->simulationPaused();
    isPaused = true;
    recoverableFailureFlag = false;
  }

  // If the exit time has been exceeded, exit.
  if (tiaParams.exitTime != 0.0 &&
      secRCPtr_->currentTime > tiaParams.exitTime)
  {
    string msg = "N_ANP_Transient::loopProcess(): "
      "  Exit Time exceeded.  Exiting transient loop\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);
    recoverableFailureFlag = false;
  }

  if (tiaParams.exitStep != -1 &&
        static_cast<int>(stepNumber) == tiaParams.exitStep)
  {
    string msg = "N_ANP_Transient::loopProcess(): "
      "  Exit Step.  Exiting transient loop\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);
    recoverableFailureFlag = false;
  }

  return recoverableFailureFlag;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::processSuccessfulDCOP()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::processSuccessfulDCOP()
{
  bool bsuccess = true;

  loaderRCPtr_->stepSuccess (anaManagerRCPtr_->currentMode_);

  totalNumberSuccessfulStepsTaken_ += 1;

  // Communicate down to the device level that the step has been accepted.
  // This must be done before the "times" get rotated, and certainly before
  // the vectors get rotated.  The primary target for this call is the
  // transmission line device, which needs to know when it has a valid
  // solution so it can save history.  (And now the ADC/DAC devices too)
  // needs to happen before dcopFlag_ is set to false.
  loaderRCPtr_->acceptStep();

  // Reset some settings (to switch from DCOP to transient, if not the
  // first step of a "double" DCOP.
  if ( firstDoubleDCOPStep_ () )  // (pde-only)
  {
    dcopFlag_          = true;
    anaManagerRCPtr_->currentMode_       = 0;
    integrationMethod_ = TIAMethod_NONE;
  }
  else
  {
    dcopFlag_                              = false;
    anaManagerRCPtr_->currentMode_         = 1;
    anaManagerRCPtr_->switchIntegrator_    = true;
    integrationMethod_   = initialIntegrationMethod_;
    beginningIntegration = true;
  }

  dsRCPtr_->setConstantHistory();
  dsRCPtr_->computeDividedDifferences();
  wimRCPtr_->obtainCorrectorDeriv();

  dsRCPtr_->updateSolDataArrays ();

  tranopOutputs ();

  // Now that output has been called, update the doubleDCOP step
  // if neccessary. (Only matters for pde problems)
  doubleDCOPStep_ = tiaParams.lastDCOPStep;

  //Test and save restart if necessary
  if( anaManagerRCPtr_->testRestartSaveTime_() )
  {
#ifdef Xyce_DEBUG_RESTART
    string netListFile = commandLine_.getArgumentValue("netlist");
    std::cout << "\n " << netListFile;
    std::cout << "  Calling dumpRestartData" << std::endl;
#endif
    restartMgrRCPtr_->dumpRestartData( secRCPtr_->currentTime );
#ifdef Xyce_DEBUG_RESTART
    std::cout << "  Done Calling dumpRestartData" << std::endl;
#endif
  }

  // This output call is for device-specific output, such as .OP,
  // or internal plot output from PDE(TCAD) devices.  
  loaderRCPtr_->output();

  nlsMgrRCPtr_->allocateTranSolver();
  secRCPtr_->previousCallStepSuccessful = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::processSuccessfulStep()
{
  bool bsuccess = true;
  loaderRCPtr_->stepSuccess (anaManagerRCPtr_->currentMode_);

  if( historyTrackingOn_ )
  {
    // store status of this step
    timeQueue_.push_back( secRCPtr_->currentTime );
    timeStepQueue_.push_back( secRCPtr_->currentTimeStep );
    stepStatusQueue_.push_back( 1 );
    estErrorOverTolQueue_.push_back( secRCPtr_->estOverTol_);
    nonlinearSolverStatusQueue_.push_back( secRCPtr_->newtonConvergenceStatus );
    nonlinearSolverNumIterationsQueue_.push_back( secRCPtr_->nIterations);
    nonlinearSolverMaxNormQueue_.push_back( nlsMgrRCPtr_->getMaxNormF() );
    nonlinearSolverMaxNormIndexQueue_.push_back( nlsMgrRCPtr_->getMaxNormFindex () );
  }

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  N_ANP_Transient::processSuccessfulStep()");
    std::cout << "Newton step succeeded:" << std::endl;
    std::cout << "nextSolutionPtr: " << std::endl;
    dsRCPtr_->nextSolutionPtr->printPetraObject();
    std::cout << std::endl;
  }
#endif
  // Set things up for the next time step, based on if this one was
  // successful.

  // This output call is for device-specific output (like from a PDE device,
  // outputting mesh-based tecplot files).  It will only work in parallel if on
  // a machine where all processors have I/O capability.
  // Note: this output needs to happen before the "times" get rotated.
  loaderRCPtr_->output();

  // Communicate down to the device level that the step has been accepted.
  // This must be done before the "times" get rotated, and certainly before
  // the vectors get rotated.  The primary target for this call is the
  // transmission line device, which needs to know when it has a valid
  // solution so it can save history.
  loaderRCPtr_->acceptStep();

  // current time will get updated in completeStep().  We'll save its value
  // for the moment so it can be saved if needed with the rest of the
  // solution if tiaParams.saveTimeStepsFlag is set.
  // This fixes an off by one bug in getting the right time value and
  // keeps the real solutions associated with that value too.
  double currentTime = secRCPtr_->currentTime;
  double suggestedMaxTime=0.0;
  if( maxTimeStepExpressionGiven_ )
  {
    suggestedMaxTime = maxTimeStepExpressionRCPtr_->evaluate(
      dsRCPtr_->currSolutionPtr, dsRCPtr_->currStatePtr, dsRCPtr_->currStorePtr);
  }
  secRCPtr_->updateMaxTimeStep( suggestedMaxTime );
  secRCPtr_->updateMinTimeStep();
  secRCPtr_->updateBreakPoints();

#ifdef Xyce_VERBOSE_TIME
  if (tiaParams.debugLevel > 0)
  {
    std::cout << "Transient Analysis:  accepting time step" << std::endl;
  }
#endif // Xyce_VERBOSE_TIME

  wimRCPtr_->completeStep();

#ifdef Xyce_VERBOSE_TIME
  if (tiaParams.errorAnalysisOption == 1)
  {
    std::cout.precision(15);
    std::cout << "ERROROPTION=1: TimeStepLimitedbyBP = " << tiaParams.TimeStepLimitedbyBP << "\n" << std::endl;
    std::cout << "ERROROPTION=1: NL Its =  " << secRCPtr_->nIterations << "\n" << std::endl;
    std::cout << "ERROROPTION=1: New DeltaT = " << secRCPtr_->currentTimeStep << "\n" << std::endl;
  }
#endif // Xyce_VERBOSE_TIME

  stepNumber     += 1;
  tranStepNumber += 1;
  totalNumberSuccessStepsThisParameter_ += 1;
  totalNumberSuccessfulStepsTaken_ += 1;

  secRCPtr_->numberSuccessiveFailures -= 1;
  if (secRCPtr_->numberSuccessiveFailures < 0)
    secRCPtr_->numberSuccessiveFailures = 0;

  // --------------------------------------------------------------------
  // Check to see if we have hit the end or a discontinuity point.  The
  // next discontinuity point is represented by the stopTime variable.

  // Note: make sure that when we check for discontinuity point, we take
  // into account the possibility of roundoff error.  Use the same bpTol
  // as is used in function updateBreakPoints.
  double bpTol = anaManagerRCPtr_->getBreakpointTol();

  if (tiaParams.bpEnable)
  {
    double timeDiff1 = secRCPtr_->currentTime - secRCPtr_->stopTime;
    double timeDiff2 = secRCPtr_->currentTime - tiaParams.finalTime;
    timeDiff1 = fabs(timeDiff1);
    timeDiff2 = fabs(timeDiff2);

#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      std::cout << " Checking whether to set breakpointrestartstep" << std::endl;
      std::cout << "   current - stop  = " << timeDiff1 << std::endl;
      std::cout << "   current - final = " << timeDiff2 << std::endl;
      std::cout << "   bpTol           = " << bpTol << std::endl;
      if (timeDiff1 <= bpTol && timeDiff2 > bpTol)
        std::cout << "    setting breakPointRestartStep to " << tranStepNumber;
    }
#endif
    if (timeDiff1 <= bpTol && timeDiff2 > bpTol)
    {
      anaManagerRCPtr_->breakPointRestartStep = tranStepNumber;
    }
  }

  if (anaManagerRCPtr_->breakPointRestartStep == tranStepNumber)
  {
    beginningIntegration = true;
  }
  else
  {
    beginningIntegration = false;
  }

  if (tiaParams.errorAnalysisOptionResetIter > 1)
  {
    tiaParams.errorAnalysisOptionResetIter--;
  }
  else if(tiaParams.errorAnalysisOptionResetIter == 1)
  {
    tiaParams.errorAnalysisOption = 0;
    tiaParams.errorAnalysisOptionResetIter = 0;
  }

  if (tiaParams.saveTimeStepsFlag)
  {
    dsRCPtr_->timeSteps.push_back(currentTime);
    dsRCPtr_->timeStepsBreakpointFlag.push_back(beginningIntegration);
    N_LAS_Vector * aVecPtr = new N_LAS_Vector( *(dsRCPtr_->currSolutionPtr) );
    dsRCPtr_->fastTimeSolutionVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(dsRCPtr_->currStatePtr) );
    dsRCPtr_->fastTimeStateVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(dsRCPtr_->daeQVectorPtr) );
    dsRCPtr_->fastTimeQVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(dsRCPtr_->currStorePtr) );
    dsRCPtr_->fastTimeStoreVec.push_back( aVecPtr );
  }

  // 03/16/04 tscoffe:  This is where the solution pointers are rotated.
  dsRCPtr_->updateSolDataArrays();

#ifdef Xyce_DEBUG_ANALYSIS
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
    dsRCPtr_->outputSolDataArrays();
#endif
#endif

  tranStepOutputs ();

  // Test and save restart if necessary
  if (anaManagerRCPtr_->testRestartSaveTime_())
    restartMgrRCPtr_->dumpRestartData(secRCPtr_->currentTime);

  secRCPtr_->previousCallStepSuccessful = true;

  // reset min error tracking variables
  // if the step number is less than zero we'll assume there is no valid
  // min. estimated error over tol or an associated time step.  This frees us
  // from putting Machine::Big here and trying to do something with that.
  stepNumberAtMinEstErrorOverTol = -1;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::processFailedStep()
{
  bool bsuccess = true;

  if( historyTrackingOn_ )
  {
    // store status of this step
    timeQueue_.push_back( secRCPtr_->currentTime );
    timeStepQueue_.push_back( secRCPtr_->currentTimeStep );
    stepStatusQueue_.push_back( 0 );
    estErrorOverTolQueue_.push_back( secRCPtr_->estOverTol_);
    nonlinearSolverStatusQueue_.push_back( secRCPtr_->newtonConvergenceStatus );
    nonlinearSolverNumIterationsQueue_.push_back( secRCPtr_->nIterations);
    nonlinearSolverMaxNormQueue_.push_back( nlsMgrRCPtr_->getMaxNormF() );
    nonlinearSolverMaxNormIndexQueue_.push_back( nlsMgrRCPtr_->getMaxNormFindex() );
  }

  // save some info about this step
  double estOverTol = secRCPtr_->getEstOverTol();
  if( (stepNumberAtMinEstErrorOverTol < 0) || ( estOverTol < minEstErrorOverTol ) )
  {
    // our first failed step, so automatically save info
    stepNumberAtMinEstErrorOverTol = stepNumber;
    minEstErrorOverTol = secRCPtr_->getEstOverTol();
    timeStepAtMinEstErrorOverTol = secRCPtr_->currentTimeStep;
  }

  loaderRCPtr_->stepFailure (anaManagerRCPtr_->currentMode_);

#ifdef Xyce_VERBOSE_TIME 
  // DO NOT REMOVE THIS OUTPUT LINE.  It is depended upon by the TIA/ERROROPTION
  // test case, which will fail if this output doesn't happen.
  std::cout << "Transient Analysis:  rejecting time step" << std::endl;
#endif // Xyce_VERBOSE_TIME 

  wimRCPtr_->rejectStep();

#ifdef Xyce_VERBOSE_TIME
  if (tiaParams.errorAnalysisOption == 1)
  {
    std::cout.precision(15);
    std::cout << "ERROROPTION=1: TimeStepLimitedbyBP = " << tiaParams.TimeStepLimitedbyBP << "\n" << std::endl;
    std::cout << "ERROROPTION=1: NL Its =  " << secRCPtr_->nIterations << "\n" << std::endl;
    std::cout << "ERROROPTION=1: New DeltaT = " << secRCPtr_->currentTimeStep << "\n" << std::endl;
  }
#endif // Xyce_VERBOSE_TIME

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  N_ANP_Transient::processFailedStep");
    std::cout << "Newton step failed:" << std::endl;
    std::cout << "nextSolutionPtr: " << std::endl;
    dsRCPtr_->nextSolutionPtr->printPetraObject();
    std::cout << std::endl;
  }
#endif
  totalNumberFailedStepsAttempted_  += 1;
  secRCPtr_->numberSuccessiveFailures += 1;

  if (secRCPtr_->currentTimeStep <= secRCPtr_->minTimeStep)
  {
#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      std::cout << "currTimeStep: " << secRCPtr_->currentTimeStep << std::endl;
      std::cout << "minTimeStep:  " << secRCPtr_->minTimeStep << std::endl;
    }
#endif

    // before we exit with a time step too small error, check if minTimeStepRecoveryCounter is greater than zero
    // if it is, then the user wants us to try to accept the step that had the minimum
    // estimated error over tol.
    if( tiaParams.minTimeStepRecoveryCounter > 0 )
    {
      ostringstream msgStream;
      msgStream
        << "N_ANP_Transient::processFailedStep:  Attempting to retake and accept step where estimated error over tolerance was: "
        << minEstErrorOverTol
        << " and time step was: "
        << timeStepAtMinEstErrorOverTol << std::endl;

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msgStream.str() );
      tiaParams.minTimeStepRecoveryCounter--;
      bsuccess = retakeAndAcceptTimeStep( timeStepAtMinEstErrorOverTol );
    }
    else
    {
      outputQueuedData();
      string msg, msg2;
      msg = "\nN_ANP_Transient::processFailedStep:"
        " Time step too small near step number: ";
      msg2  = "  Exiting transient loop.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg, stepNumber, msg2);

      bsuccess = false;
    }
  }

  if (tiaParams.constantStepSize)
  {
    outputQueuedData();
    string msg = "N_ANP_Transient::processFailedStep: ";
    msg += "  Newton solver failed in constant time step mode.  "
      "Exiting transient loop.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);

    bsuccess = false;
  }

  if (tiaParams.exitStep != -1 &&
          static_cast<int>(totalNumberSuccessStepsThisParameter_) == (tiaParams.exitStep-1))
  {
    outputQueuedData();
    string msg = "N_ANP_Transient::processFailedStep: ";
    msg += "  Exit Step.  Exiting transient loop\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);
    bsuccess = false;
  }
#ifdef Xyce_VERBOSE_TIME
  if (!bsuccess)
  {
    endTRANtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
    tranStats = saveLoopInfo();
    finalVerboseOutput();
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::processFailedDCOP
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::processFailedDCOP()
{
  bool bsuccess = true;

  loaderRCPtr_->stepFailure (anaManagerRCPtr_->currentMode_);

  // DC Operating Point failed.
  (totalNumberFailedStepsAttempted_)++;
  (secRCPtr_->numberSuccessiveFailures)++;

  string msg = "N_ANP_Transient::processFailedDCOP - ";
  msg += "DC Operating Point Failed.  Exiting transient loop.\n";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::finish()
{
  bool bsuccess = true;

  if (tiaParams.saveTimeStepsFlag && anaManagerRCPtr_->getHBFlag())
  {
    N_TIA_DataStore * dsPtr_ = dsRCPtr_.get();
    dsPtr_->timeSteps.push_back(secRCPtr_->currentTime);
    dsPtr_->timeStepsBreakpointFlag.push_back(beginningIntegration);
    N_LAS_Vector * aVecPtr = new N_LAS_Vector( *(dsPtr_->currSolutionPtr) );
    dsPtr_->fastTimeSolutionVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(dsPtr_->currStatePtr) );
    dsPtr_->fastTimeStateVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(dsPtr_->daeQVectorPtr) );
    dsPtr_->fastTimeQVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(dsPtr_->currStorePtr) );
    dsPtr_->fastTimeStoreVec.push_back( aVecPtr );
  }

  if (!isPaused)
  {
#ifdef Xyce_DEBUG_ANALYSIS
    std::cout << "Calling finishOutput" << std::endl;
#endif
    outputMgrAdapterRCPtr_->finishOutput();

    // This output call is for device-specific output (like from a PDE device,
    // outputting mesh-based tecplot files).  It will only work in parallel if on
    // a machine where all processors have I/O capability, as devices are
    // local to a processor.
    loaderRCPtr_->finishOutput();

    finalVerboseOutput();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::handlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool N_ANP_Transient::handlePredictor()
{
  dsRCPtr_->setErrorWtVector();
  wimRCPtr_->obtainPredictor();
  wimRCPtr_->obtainPredictorDeriv();

#ifdef Xyce_DEBUG_ANALYSIS
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
           "  N_ANP_Transient::handlePredictor");
    dsRCPtr_->outputPredictedSolution();
    dsRCPtr_->outputPredictedDerivative();
  }
#endif
#endif

  // Now, in case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  loaderRCPtr_->startTimeStep ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::resetForStepAnalysis()
// Purpose       : When doing a .STEP sweep, some data must be reset to its
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool N_ANP_Transient::resetForStepAnalysis()
{
  totalNumberSuccessStepsThisParameter_ = 0;
  stepNumber = 0;
  beginningIntegration = true;

  nlsMgrRCPtr_->resetAll(DC_OP);
  secRCPtr_->resetAll();
  dcopFlag_ = false;
  anaManagerRCPtr_->nextOutputTime_ = 0.0;

  if (historyTrackingOn_)
  {
    // these set_size calls will also "reset" the queue.
    timeQueue_.set_size( queueSize_ );
    timeStepQueue_.set_size( queueSize_ );
    stepStatusQueue_.set_size( queueSize_ );
    estErrorOverTolQueue_.set_size( queueSize_ );
    nonlinearSolverStatusQueue_.set_size( queueSize_ );
    nonlinearSolverNumIterationsQueue_.set_size( queueSize_ );
    nonlinearSolverMaxNormQueue_.set_size( queueSize_ );
    nonlinearSolverMaxNormIndexQueue_.set_size( queueSize_ );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::resetForHB()
// Purpose       : When doing initial transient run, some analyses require
//                 a reset function
//                 they can fill this in if needed.
// Special Notes :
// Scope         : public
// Creator       : T. Mei, SNL, Parallel Computational Sciences
// Creation Date : 2/26/09
//-----------------------------------------------------------------------------
bool N_ANP_Transient::resetForHB()
{
  nlsMgrRCPtr_->resetAll(DC_OP);
  secRCPtr_->resetAll();
  dsRCPtr_->resetAll();
  dcopFlag_ = false;
  anaManagerRCPtr_->nextOutputTime_ = 0.0;

  if (historyTrackingOn_)
  {
    // these set_size calls will also "reset" the queue.
    timeQueue_.set_size( queueSize_ );
    timeStepQueue_.set_size( queueSize_ );
    stepStatusQueue_.set_size( queueSize_ );
    estErrorOverTolQueue_.set_size( queueSize_ );
    nonlinearSolverStatusQueue_.set_size( queueSize_ );
    nonlinearSolverNumIterationsQueue_.set_size( queueSize_ );
    nonlinearSolverMaxNormQueue_.set_size( queueSize_ );
    nonlinearSolverMaxNormIndexQueue_.set_size( queueSize_ );
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::finalVerboseOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::finalVerboseOutput()
{
  bool bsuccess = true;

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      string(" ***** Problem read in and set up time: "),
               anaManagerRCPtr_->solverStartTime_, string(" seconds"));

  if (anaManagerRCPtr_->analysis == ANP_MODE_TRANSIENT)
  {
    double time;
    if (anaManagerRCPtr_->startTRANtime > startDCOPtime)
    {
      time = anaManagerRCPtr_->startTRANtime - startDCOPtime;
    }
    else
    {
      time = endTRANtime - startDCOPtime;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
        string(" ***** DCOP time: "),
           time, string(" seconds.  Breakdown follows:"));

    printLoopInfo(0, dcStats);
  }

  if (anaManagerRCPtr_->analysis == ANP_MODE_TRANSIENT &&
      endTRANtime >= anaManagerRCPtr_->startTRANtime)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
        string(" ***** Transient Stepping time: "),
           endTRANtime-anaManagerRCPtr_->startTRANtime,
           string(" seconds.  Breakdown follows:"));

    printLoopInfo(dcStats, tranStats);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::outputQueuedData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
void N_ANP_Transient::outputQueuedData()
{
  if( historyTrackingOn_ )
  {
#ifndef Xyce_PARALLEL_MPI
    // if necessary set up the names map.
    Xyce::NodeNamePairMap & nodeNames = outputMgrAdapterRCPtr_->getAllNodes();
    int nodeNameSize = nodeNames.size();
    nameVec_.resize(nodeNameSize+1,"gnd");
    Xyce::NodeNamePairMap::iterator mapI, mapEnd;
    mapEnd = nodeNames.end();
    mapI =  nodeNames.begin();
    for ( ; mapI != mapEnd ; ++mapI)
    {
      nameVec_[(*mapI).second.first] = (*mapI).first;
    }
#endif

    // get the current non-linear solver return codes
    N_NLS_ReturnCodes nlReturnCodes = nlsMgrRCPtr_->getReturnCodes();

    ostringstream outStringStream;
    outStringStream << " *** Transient failure history: " << std::endl;
    //                  123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    if (tiaParams.errorAnalysisOption == 1)
    {
      // truncation error is not used here so est err over tol is not useful in the output
      //                  123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
      outStringStream << "Time        Time      Step      Non-Linear Solver      node    node" << std::endl;
      outStringStream << "(sec)       Step     Status   Status   Iters   ||F||   index   name" << std::endl;
    }
    else
    {
      //                  123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
      outStringStream << "Time        Time      Step   EstErr      Non-Linear Solver      node     node" << std::endl;
      outStringStream << "(sec)       Step     Status  OverTol   Status  Iters   ||F||    index    name"
        << std::endl;
    }

    for( int i=0; i<queueSize_; i++ )
    {
      int fieldWidth=10;
      outStringStream << scientific << setprecision(fieldWidth-7) << setfill(' ') << right << setw( fieldWidth )
        << timeQueue_.at_from_tail(i) << "  "
        << timeStepQueue_.at_from_tail(i) << "  ";
      if( stepStatusQueue_.at_from_tail(i) == 1 )
      {
        outStringStream << "pass  ";
      }
      else
      {
        outStringStream << "fail  ";
      }
      if (tiaParams.errorAnalysisOption == 1)
      {
      }
      else
      {
        outStringStream << estErrorOverTolQueue_.at_from_tail(i) << "  ";
      }
      int nlStatus = nonlinearSolverStatusQueue_.at_from_tail(i);
      outStringStream << setw(7) << right;
      if( nlStatus == nlReturnCodes.normTooSmall )
      {
        outStringStream << "P:s nrm";
      }
      else if( nlStatus == nlReturnCodes.normalConvergence)
      {
        outStringStream << "pass   ";
      }
      else if( nlStatus == nlReturnCodes.nearConvergence )
      {
        outStringStream << "P:near ";
      }
      else if( nlStatus == nlReturnCodes.smallUpdate )
      {
        outStringStream << "P:s up ";
      }
      else if( nlStatus == nlReturnCodes.nanFail )
      {
        outStringStream << "F:NaN  ";
      }
      else if( nlStatus == nlReturnCodes.tooManySteps )
      {
        outStringStream << "F:max s";
      }
      else if( nlStatus == nlReturnCodes.tooManyTranSteps )
      {
        outStringStream << "F:max s";
      }
      else if( nlStatus == nlReturnCodes.updateTooBig )
      {
        outStringStream << "F:big u";
      }
      else if( nlStatus == nlReturnCodes.stalled )
      {
        outStringStream << "F:stall";
      }
      else if( nlStatus == nlReturnCodes.wrmsExactZero )
      {
        outStringStream << "F:n zro";
      }
      else if( nlStatus == nlReturnCodes.innerSolveFailed )
      {
        outStringStream << "F:in Fl";
      }
      else
      {
        outStringStream << "code=" <<
            nonlinearSolverStatusQueue_.at_from_tail(i) << "  ";
      }

      outStringStream << right << setw( 4 )
        << nonlinearSolverNumIterationsQueue_.at_from_tail(i) << "  "
        << nonlinearSolverMaxNormQueue_.at_from_tail(i) ;

      int outIndex = nonlinearSolverMaxNormIndexQueue_.at_from_tail(i) ;
      outStringStream << right << fixed << setw( 7 ) << outIndex;

      std::string outIndexName("");
      if (!(nameVec_.empty()))
      {
        int nsize = nameVec_.size();
        if (nsize >  outIndex && outIndex >=0)
        {
          outIndexName = nameVec_[outIndex];
        }
      }
      else
      {
        outIndexName = "N/A";
      }
      outStringStream << left << "    " << outIndexName;
      outStringStream << std::endl;
    }

    string message(outStringStream.str());
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, message );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::takeAnIntegrationStep_
// Purpose       : Take a transient integration step.
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void N_ANP_Transient::takeAnIntegrationStep_()
{
  handlePredictor();
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();

#ifdef Xyce_DEBUG_ANALYSIS
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
        "N_ANP_Transient::takeAnIntegrationStep_:  newtonConvergenceStatus = ", 
        secRCPtr_->newtonConvergenceStatus);
#endif

  dsRCPtr_->stepLinearCombo ();
  gatherStepStatistics_ ();
  secRCPtr_->nIterations = nlsMgrRCPtr_->getNumIterations();
  secRCPtr_->evaluateStepError ();

  return;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::twoLevelStep
//
// Purpose       : Take a transient integration step, for inner 2-level solve.
//
// Special Notes : Same as takeAnIntegrationStep, but does not do the
//                 prediction. (that is in handlePredictor_).
//
//                 The prediction is handled separately, as for a 2-level
//                 solve, you only want to do the prediction on the first
//                 solve of the attempted time step.
//
// Scope         : public
// Creator       :
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_Transient::twoLevelStep ()
{
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();
  dsRCPtr_->stepLinearCombo ();

  gatherStepStatistics_ ();
  secRCPtr_->evaluateStepError ();

  endTRANtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
  tranStats = saveLoopInfo();

  return (secRCPtr_->stepAttemptStatus);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::retakeAndAcceptTimeStep
// Purpose       : Do a requested time step and accept it
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 01/23/09
//-----------------------------------------------------------------------------
bool N_ANP_Transient::retakeAndAcceptTimeStep( double aTimeStep )
{
  bool bsuccess=true;
  // This function was put in place to handle the following situation.
  // If a series of time steps are rejected because the estimated error
  // over tol. is too high, Xyce may exit if the next time step falls under
  // the minTimeStep.  A user could adjust the tolerances of a simulation
  // to try and get around this, but in UQ studies where there are 1,000s of
  // simulations, it can be difficult to get all of them to run.  So,
  // if the user has set the netlist option timeint MINTIMESTEPRECOVERY=<int>
  // and we're about to exit because the time step has fallen under the minTimeStep
  // we will try and retake the time step that had the min. est. error over tol.
  // and accept that.
  //
  // At this point, N_ANP_Transient::processFailedStep() has already determined
  // that the preconditions outlined above are in place (i.e. the user requested
  // this and we're about to exit with a time step too small error), so lets
  // try the step that had the min. est error over tol.

  // set time step
  secRCPtr_->currentTimeStep = timeStepAtMinEstErrorOverTol;

  // take step
  takeAnIntegrationStep_();

  // can't accept step if the non-linear solver failed
  if(secRCPtr_->nIterations==0)
  {
    string msg, msg2;
    msg = "\nN_ANP_Transient::retakeAndAcceptTimeStep:"
        " Time step too small near step number: ";
    msg2  = "  Exiting transient loop.\n";

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg, stepNumber, msg2);

    bsuccess = false;
  }
  else
  {
    processSuccessfulStep();
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::printStepHeader()
// Purpose       : Prints out time step information.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void N_ANP_Transient::printStepHeader()
{
#ifdef Xyce_VERBOSE_TIME

#ifdef Xyce_DEBUG_ANALYSIS
  string netListFile = commandLine_.getArgumentValue("netlist");
  string banner = "***** " + netListFile + "  ";
#else
  string banner = "***** ";
#endif
  string crStr("\n");
  string tmpStr;

  //N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, crStr);

  if (dcopFlag_)
  {
    tmpStr = banner + "Start of DCOP STEP                        # ";
  }
  else
  {
    if (beginningIntegration)
    {
      if (secRCPtr_->currentTime == tiaParams.initialTime)
      {
        tmpStr = banner + "Start of Time Step (INITIAL STEP)         # ";
      }
      else
      {
        tmpStr = banner + "Start of Time Step (DISCONTINUITY STEP)   # ";
      }
    }
    else
    {
      if (!secRCPtr_->stepAttemptStatus)
      {
        tmpStr = banner + "Start of Time Step (AFTER FAILED STEP)    # ";
      }
      else
      {
        tmpStr = banner + "Start of Time Step (AFTER SUCCESS STEP)   # ";
      }
    }
  }

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmpStr,
                          totalNumberSuccessStepsThisParameter_ + 1);
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::printProgress()
// Purpose       : Outputs run completion percentage and estimated
//                 time-to-completion.
//
// Special Notes : This will need some fixing to work with .STEP.
//
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 06/07/2002
//-----------------------------------------------------------------------------
void N_ANP_Transient::printProgress()
{
  if (anaManagerRCPtr_->progressFlag_)
  {
    // Percent of the overall requested simulation time which is completed.
    double percentComplete;
    // Average CPU time per time step.
    double aveCPUTimePerStep;
    // Average clock time per time step.
    double aveSimTimePerStep;
    // Estimated CPU time to complete the simulation.
    double estCompletionTime = 0.0;

    string banner(" ***** ");
    string crStr("\n");
    string tmpStr;

    // Report the beginning of the DC OP calculation.  First call in OP.
    if (dcopFlag_)
    {
      startDCOPtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
      if (gui_)
      {
        tmpStr = "Xyce in DC Operating Point Calculation";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::GUI_PROGRESS, tmpStr);
      }
      tmpStr = banner + "Beginning DC Operating Point Calculation...\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmpStr);
    }
    else if (firstTime && totalNumberSuccessStepsThisParameter_ == 1)
    {
      anaManagerRCPtr_->startTRANtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
      dcStats = saveLoopInfo();
      firstTime = false;
      if (gui_)
      {
        tmpStr = "Xyce in Transient Calculation";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::GUI_PROGRESS, tmpStr);
      }
      tmpStr = banner + "Beginning Transient Calculation...\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmpStr);
    }
    if (anaManagerRCPtr_->analysis == ANP_MODE_TRANSIENT && totalNumberSuccessStepsThisParameter_  > 0)
    {
      if (startSimTime == -1.0)
        startSimTime = secRCPtr_->currentTime;

      double diff1 = fabs(secRCPtr_->currentTime - tiaParams.initialTime);
      double diff2 = fabs(tiaParams.finalTime - tiaParams.initialTime);

      // 02/05/08 tscoffe:  optimization differences after the TIA_REFACTOR
      // caused a floating point difference here and this resolves that
      // difference.  After the refactor is merged back onto the main trunk, this
      // stuff should be removed.
      percentComplete = 100.0 * diff1/diff2;

      if (fabs(percentComplete - oldPercentComplete) > 1.0)
      {
        oldPercentComplete = percentComplete;

        aveCPUTimePerStep = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime() /
          (totalNumberSuccessfulStepsTaken_ + 1);
        aveSimTimePerStep = (secRCPtr_->currentTime - startSimTime) /
          (totalNumberSuccessfulStepsTaken_ + 1);

        if (aveSimTimePerStep > N_UTL_MachineDependentParams::MachineEpsilon())
          estCompletionTime = aveCPUTimePerStep *
            fabs(tiaParams.finalTime - secRCPtr_->currentTime) /
            aveSimTimePerStep;

        if (!gui_)
        {
          tmpStr = banner + "Percent complete: ";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmpStr,
                                 percentComplete, " %");
        }

        if (estCompletionTime > N_UTL_MachineDependentParams::MachineEpsilon())
        {
          unsigned int days, hours, minutes, seconds;
          days    = static_cast<int> (estCompletionTime / 86400);
          hours   = static_cast<int> ((estCompletionTime - days * 86400) / 3600);
          minutes = static_cast<int> ((estCompletionTime - days * 86400 - hours * 3600) / 60);
          seconds = static_cast<int> (estCompletionTime - days * 86400 - hours * 3600 - minutes * 60);

          char timeStr[256];

#ifdef Xyce_PARALLEL_MPI

          // get current local system time
          time_t t = time( NULL );
          struct tm * now = localtime( &t );

          // format and display output
          if ( ( t != (time_t)-1 ) && ( strftime( timeStr, 255, "%c", now ) != 0 ) )
          {
            tmpStr = banner + "Current system time: " + timeStr;
          }

          else
          {
            tmpStr = banner + "Current system time could not be determined.";
          }

          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmpStr );

#endif

          if (days > 0)
            sprintf(timeStr, "%3d days, %2d hrs., %2d min., %2d sec.", days,
                    hours, minutes, seconds);
          else if (hours > 0)
            sprintf(timeStr, "%2d hrs., %2d min., %2d sec.", hours, minutes,
                    seconds);
          else if (minutes > 0)
            sprintf(timeStr, "%2d min., %2d sec.", minutes, seconds);
          else
            sprintf(timeStr, "%2d sec.", seconds);

          if (gui_)
          {
            ostringstream ost;
            ost << "Xyce transient ";
            if (percentComplete < 10)
              ost.precision(2);
            else
              ost.precision(3);
            ost << percentComplete;
            ost << "%% complete ... Estimated time to completion: ";
            ost << timeStr;
            ost << crStr;
            tmpStr = ost.str();
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::GUI_PROGRESS, tmpStr);
          }
          else
          {
            tmpStr = banner + "Estimated time to completion: ";
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmpStr + timeStr +
                                   crStr);
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::noopOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void N_ANP_Transient::noopOutputs ()
{
  if( testOutputTime_() )
  {
    if ( !firstDoubleDCOPStep_ () )
    {
      bool printIC = true;
//      if ( (!Teuchos::is_null(mpdeMgrPtr_)) && (mpdeMgrPtr_->getMPDEIcFlag()) )
//      {
        // don't print the dc-op in the MPDE IC calculation
//        printIC = false;
//     }
      if (printIC)
      {
        outputMgrAdapterRCPtr_->tranOutput(
            secRCPtr_->currentTime, *dsRCPtr_->currSolutionPtr, *dsRCPtr_->currStatePtr, *dsRCPtr_->currStorePtr);
      }
      if ( ( !Teuchos::is_null(anaManagerRCPtr_->mpdeMgrPtr_)  )  &&
           (  anaManagerRCPtr_->mpdeMgrPtr_->getMPDEFlag()     )  &&
           ( !(anaManagerRCPtr_->mpdeMgrPtr_->getMPDEIcFlag()) )
         )
      {
        // output the operating point if we can in MPDE mode
        outputMgrAdapterRCPtr_->outputMPDE(secRCPtr_->currentTime, *dsRCPtr_->nextSolutionPtr);
      }
    }
    updateOutputTime_(secRCPtr_->currentTime);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::tranopOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void N_ANP_Transient::tranopOutputs ()
{
  //Test and output if necessary
  if( testOutputTime_() )
  {
    // Make sure this isn't the NLP step of a PDE DCOP.
    if ( !firstDoubleDCOPStep_ () )
    {
#ifdef Xyce_DEBUG_ANALYSIS
      std::cout << "Calling conventional TRANOP outputs!" << std::endl;
#endif
      outputMgrAdapterRCPtr_->tranOutput(secRCPtr_->currentTime, *(dsRCPtr_->currSolutionPtr), *dsRCPtr_->currStatePtr, *dsRCPtr_->currStorePtr ) ;

      outputMgrAdapterRCPtr_->outputMPDE(secRCPtr_->currentTime, *(dsRCPtr_->nextSolutionPtr));
      // check this - may be a mistake.
    }
    updateOutputTime_(secRCPtr_->currentTime);
  }

  // SAVE and DCOP restart:
  if ( anaManagerRCPtr_->testDCOPOutputTime_() || anaManagerRCPtr_->testSaveOutputTime_() )
  {
    outputMgrAdapterRCPtr_->outputDCOP( *(dsRCPtr_->currSolutionPtr) );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::tranStepOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void N_ANP_Transient::tranStepOutputs ()
{
  // The printOutputSolution function will sometimes perform
  // interpolations between the current step and the previous step,
  // if the integrator is using a high order of integration during
  // a new-DAE simulation.

  // If the user has specified tstart, and this happens to be
  // the first output of the simulation (at t=tstart, instead of
  // t=0) the force the interpolator to *not* interpolate, even
  // if it is doing so normally.
  bool doNotInterpolate=
      (beginningIntegration && tiaParams.tStartGiven && firstTranOutput_);

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    if (doNotInterpolate)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
          "doNotInterpolate is TRUE!");
    }
    else
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
          "doNotInterpolate is FALSE!");
    }
  }
#endif

  if (testOutputTime_())
  {
    if ((!Teuchos::is_null(anaManagerRCPtr_->mpdeMgrPtr_)) && (anaManagerRCPtr_->mpdeMgrPtr_->blockAnalysisFlag()==true))
    {
#ifdef Xyce_DEBUG_ANALYSIS
      if (tiaParams.debugLevel > 0)
      {
        std::cout << "Calling MPDE outputs!" << std::endl;
      }
#endif
      outputMgrAdapterRCPtr_->outputMPDE(secRCPtr_->currentTime, *(dsRCPtr_->nextSolutionPtr));

      // If we are actually in the MPDE phase, rather than the initial
      // condition, then output here.
      if (anaManagerRCPtr_->mpdeMgrPtr_->getMPDEFlag()==true && tiaParams.outputInterpMPDE)
      {
        if (anaManagerRCPtr_->mpdeMgrPtr_->getWaMPDEFlag()==false)
        {
          std::vector<double> fastTimes = anaManagerRCPtr_->mpdeMgrPtr_->getFastTimePoints();
          wimRCPtr_->printMPDEOutputSolution(
                  outputMgrAdapterRCPtr_, secRCPtr_->currentTime,
                  dsRCPtr_->currSolutionPtr,
                  fastTimes );
        }
        else
        {
          std::vector<double> fastTimes = anaManagerRCPtr_->mpdeMgrPtr_->getFastTimePoints();
          int phiGID = anaManagerRCPtr_->mpdeMgrPtr_->getPhiGID();
          wimRCPtr_->printWaMPDEOutputSolution(
                  outputMgrAdapterRCPtr_, secRCPtr_->currentTime,
                  dsRCPtr_->currSolutionPtr,
                  fastTimes, phiGID );
        }
      }
      // uncomment the following else block to get mpde initial
      // condition calculations printed out
      else
      {
        // try normal new dae output
        computeOutputInterpolationTimes_(secRCPtr_->currentTime);
        wimRCPtr_->printOutputSolution(
              outputMgrAdapterRCPtr_, secRCPtr_->currentTime,
              dsRCPtr_->currSolutionPtr,
              doNotInterpolate,
              outputInterpolationTimes_,
              false) ;
      }
    }
    else
    {
      computeOutputInterpolationTimes_(secRCPtr_->currentTime);
      wimRCPtr_->printOutputSolution(
              outputMgrAdapterRCPtr_, secRCPtr_->currentTime,
              dsRCPtr_->currSolutionPtr,
              doNotInterpolate,
              outputInterpolationTimes_,
              false) ;
    }

    firstTranOutput_ = false;

    updateOutputTime_(secRCPtr_->currentTime);
  }
  else
  {
    // just call the output manager to udpate non print line elements
    computeOutputInterpolationTimes_(secRCPtr_->currentTime);
    wimRCPtr_->printOutputSolution(
              outputMgrAdapterRCPtr_, secRCPtr_->currentTime,
              dsRCPtr_->currSolutionPtr,
              doNotInterpolate,
              outputInterpolationTimes_,
              true) ;
  }

  // SAVE output:
  if ( anaManagerRCPtr_->testSaveOutputTime_() )
  {
    wimRCPtr_->saveOutputSolution(
                outputMgrAdapterRCPtr_,
                dsRCPtr_->currSolutionPtr,
                anaManagerRCPtr_->saveTime_, doNotInterpolate);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::testOutputTime_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel ComputationalSciences.
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_Transient::testOutputTime_()
{
  bool flag;

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    std::cout.width(21); std::cout.precision(13); std::cout.setf(ios::scientific);
    std::cout << "N_ANP_Transient::testOutputTime.  secPtr_->currentTime = " << secRCPtr_->currentTime << std::endl;
    std::cout << "N_ANP_Transient::testOutputTime.  tStart      = " << tiaParams.tStart << std::endl;
    std::cout << "N_ANP_Transient::testOutputTime.  nextOutputTime_ = " << anaManagerRCPtr_->nextOutputTime_ << std::endl;
    for (int i=0;i<anaManagerRCPtr_->outputIntervals_.size();++i)
    {
      std::cout << "outputIntervals_["<<i<<"].first = " << anaManagerRCPtr_->outputIntervals_[i].first;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
#endif

  if (secRCPtr_->currentTime < tiaParams.tStart)
  {
    flag = false;
  }
  else if (secRCPtr_->currentTime >= tiaParams.finalTime)
  {
    flag = true;
  }
  else if (anaManagerRCPtr_->initialOutputInterval_ == 0.0)
  {
    flag = true;
  }
  else if (secRCPtr_->currentTime < anaManagerRCPtr_->nextOutputTime_)
  {
    flag = false;
  }
  else
  {
    flag = true;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::updateOutputTime_
// Purpose       : Advance output time so it is the next one after currTime
// Special Notes : formerly part of testOutputTime_
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 02/14/07
//-----------------------------------------------------------------------------
void N_ANP_Transient::updateOutputTime_(double currTime)
{
  // We only need to bother with this if the user has specified output control
  // options
  if (anaManagerRCPtr_->outputIntervalSpecified_())
  {
    if (anaManagerRCPtr_->outputIntervals_.empty())
    {
      // Only an initial interval was specified
      while (anaManagerRCPtr_->nextOutputTime_ <= currTime)
        anaManagerRCPtr_->nextOutputTime_ += anaManagerRCPtr_->initialOutputInterval_;
    }
    else if (currTime < anaManagerRCPtr_->outputIntervals_[0].first)
    {
      // Multiple intervals, but we're still in the first one
      while (anaManagerRCPtr_->nextOutputTime_ <= currTime)
        anaManagerRCPtr_->nextOutputTime_ += anaManagerRCPtr_->initialOutputInterval_;
      if (anaManagerRCPtr_->nextOutputTime_ > anaManagerRCPtr_->outputIntervals_[0].first)
        anaManagerRCPtr_->nextOutputTime_ = anaManagerRCPtr_->outputIntervals_[0].first;
    }
    else
    {
      // Multiple intervals specified, we're past the first one
      pair<double, double> currInterval, nextInterval;
      int size = anaManagerRCPtr_->outputIntervals_.size();
      for (int i = 0; i < size; ++i)
        if (anaManagerRCPtr_->outputIntervals_[i].first <= currTime)
        {
          currInterval = anaManagerRCPtr_->outputIntervals_[i];
          if ((i+1) < static_cast<int>(anaManagerRCPtr_->outputIntervals_.size()))
            nextInterval = anaManagerRCPtr_->outputIntervals_[i+1];
        }
      int step = static_cast<int> ((currTime-currInterval.first) /
                                   currInterval.second);
      anaManagerRCPtr_->nextOutputTime_ = currInterval.first + (step+1)*currInterval.second;
      if (nextInterval.first && (nextInterval.first!=currInterval.first)
          && (anaManagerRCPtr_->nextOutputTime_>=nextInterval.first))
        anaManagerRCPtr_->nextOutputTime_ = nextInterval.first;
    }

    if (anaManagerRCPtr_->nextOutputTime_ >= tiaParams.finalTime)
      anaManagerRCPtr_->nextOutputTime_ = tiaParams.finalTime;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_Transient::computeOutputInterpolationTimes_
// Purpose       : When we pass "nextOutputTime_", we might have skipped
//                 over points where output was requested.  Make a list of
//                 those times so we can interpolate to them.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling.
// Creation Date : 02/14/07
//-----------------------------------------------------------------------------
void N_ANP_Transient::computeOutputInterpolationTimes_(double currTime)
{
  double t;
  outputInterpolationTimes_.clear();
    // if there are no output control options specified, we do nothing.
  if (anaManagerRCPtr_->outputIntervalSpecified_())
  {
    if (anaManagerRCPtr_->outputIntervals_.empty() || currTime <= anaManagerRCPtr_->outputIntervals_[0].first)
    {
      // we're either in the first interval, or there only *is* one interval
      t=anaManagerRCPtr_->nextOutputTime_;
      while (t < currTime)
      {
        outputInterpolationTimes_.push_back(t);
        t += anaManagerRCPtr_->initialOutputInterval_;
      }

      if( (t - currTime) <= 100 * N_UTL_MachineDependentParams::MachinePrecision() )
      {
        outputInterpolationTimes_.push_back(currTime);
      }
  
      if( (t - tiaParams.finalTime) > 100 * N_UTL_MachineDependentParams::MachinePrecision() )
      {
        // tack on finalTime or we'll miss it
        outputInterpolationTimes_.push_back( tiaParams.finalTime );
      }
    }
    else // we're out of the first interval, and there is more than one
    {
      int outInt,lastInt;

      lastInt=anaManagerRCPtr_->outputIntervals_.size()-1;
      // find which interval nextOutputTime_ is in
      for (outInt=0;
           outInt<lastInt&&anaManagerRCPtr_->outputIntervals_[outInt+1].first<=anaManagerRCPtr_->nextOutputTime_;
           ++outInt) ;

      t=anaManagerRCPtr_->nextOutputTime_;
      while (t <=  currTime)
      {
        outputInterpolationTimes_.push_back(t);
        t += anaManagerRCPtr_->outputIntervals_[outInt].second;
        if (outInt != lastInt && t >= anaManagerRCPtr_->outputIntervals_[outInt+1].first)
        {
          ++outInt;
          t = anaManagerRCPtr_->outputIntervals_[outInt].first;
        }
      }
    }
  }
}

