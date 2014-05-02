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
// Filename      : $RCSfile: N_ANP_Transient.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.112.2.1 $
// Revision Date  : $Date: 2014/03/04 23:50:53 $
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
#include <N_ERH_Progress.h>

#include <N_IO_CmdParse.h>

#include<N_ANP_Transient.h>

#include<N_UTL_ExpressionData.h>
#include<N_NLS_ReturnCodes.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : Transient::Transient( AnalysisManager * )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
Transient::Transient( AnalysisManager * anaManagerPtr ) :
  AnalysisBase(anaManagerPtr),
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
// Function      : Transient::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .TRAN statement.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/21/10
//-----------------------------------------------------------------------------
bool Transient::setAnalysisParams(const N_UTL_OptionBlock & paramsBlock)
{
  std::list<N_UTL_Param>::const_iterator it_tp;
  std::list<N_UTL_Param>::const_iterator first = paramsBlock.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last = paramsBlock.getParams().end();
  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag()      == "TSTART")
    {
      tiaParams.tStart = it_tp->getImmutableValue<double>();
      tiaParams.tStartGiven = true;
    }
    else if (it_tp->uTag() == "TSTOP")
    {
      tiaParams.finalTime   = it_tp->getImmutableValue<double>();
    }
    else if (it_tp->uTag() == "TSTEP")
    {
      tiaParams.userSpecified_startingTimeStep = it_tp->getImmutableValue<double>();
    }
    else if (it_tp->uTag() == "NOOP" ||
             it_tp->uTag() == "UIC")
    {
      tiaParams.NOOP = 1;
    }
    else if (it_tp->uTag() == "DTMAX")
    {
      tiaParams.maxTimeStep = it_tp->getImmutableValue<double>();
      tiaParams.maxTimeStepGiven = true;
#ifdef Xyce_DEBUG_ANALYSIS
      if (tiaParams.debugLevel > 0)
      {
        dout() << "setting maxTimeStep = " << tiaParams.maxTimeStep << std::endl;
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
      // happen in Transient class as it will ultimately use this
      maxTimeStepExpressionGiven_ = true;
      maxTimeStepExpressionAsString_ = it_tp->stringValue();
    }
  }

  if (tiaParams.finalTime <= tiaParams.tStart || tiaParams.finalTime <= 0 || tiaParams.tStart < 0)
  {
    std::ostringstream ost;
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
  if (tiaParams.debugLevel > 0)
  {
    dout() << std::endl
           << section_divider << std::endl
           << "  Transient simulation parameters" << std::endl
           << "  initial time = " << tiaParams.initialTime << std::endl
           << "  final time   = " << tiaParams.finalTime << std::endl
           << "  starting time step = " << tiaParams.userSpecified_startingTimeStep << std::endl
           << "  tStart (time of first output) = " << tiaParams.tStart << std::endl;
    
    if (!(tiaParams.NOOP))
    {
      dout() << "  NOOP/UIC is NOT set" << std::endl;
    }
    else
    {
      dout() << "  NOOP/UIC is set" << std::endl;
    }

    dout() << section_divider << std::endl;
  }
#endif

  if( maxTimeStepExpressionGiven_ )
  {
    // set-up expression object for getting a user specified, time dependent
    // max time step value
    maxTimeStepExpressionGiven_ = true;
    maxTimeStepExpressionRCPtr_ = rcp( new N_UTL_ExpressionData( maxTimeStepExpressionAsString_,  *outputMgrAdapterRCPtr_->getOutputMgrPtr() ) );
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
// Function      : Transient::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::run()
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
// Function      : Transient::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::init()
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
      loaderRCPtr_->getInitialQnorm(anaManagerRCPtr_->getTIADataStore()->innerErrorInfoVec);

//      secRCPtr_->initialTime = secRCPtr_->currentTime;
      wimRCPtr_->createTimeIntegMethod(integrationMethod_);
      wimRCPtr_->initialize();

      dcopFlag_ = false;
      anaManagerRCPtr_->currentMode_ = 1;

#ifdef Xyce_PARALLEL_MPI
      // Update vectors with off proc values.
      lasSystemRCPtr_->updateExternValsSolnVector(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr);
      lasSystemRCPtr_->updateExternValsSolnVector(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr);
      lasSystemRCPtr_->updateExternValsStateVector(anaManagerRCPtr_->getTIADataStore()->nextStatePtr);
      lasSystemRCPtr_->updateExternValsStateVector(anaManagerRCPtr_->getTIADataStore()->currStatePtr);
      lasSystemRCPtr_->updateExternValsStoreVector(anaManagerRCPtr_->getTIADataStore()->nextStorePtr);
      lasSystemRCPtr_->updateExternValsStoreVector(anaManagerRCPtr_->getTIADataStore()->currStorePtr);
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
      loaderRCPtr_->setInitialGuess (anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr);

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
          ( *(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr),*(anaManagerRCPtr_->getTIADataStore()->flagSolutionPtr));
      }

      if (!dcopFlag_ && !anaManagerRCPtr_->getMPDEFlag())
      {
        // this "initializeProblem" call is to set the IC's in devices that
        // have them.  This is done IN PLACE of the operating point.
        loaderRCPtr_->initializeProblem
              ((anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->currSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->lastSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->currStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->lastStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStateDerivPtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->currStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->lastStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->daeQVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->daeFVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->dFdxdVpVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->dQdxdVpVectorPtr) );

        // Do this to populate the q-vector:
        // since we're also skipping the DC OP, this call will
        // also force the loader to propagate the state vector to the
        // device manager which updates state vector seen by the devices.
        assemblerRCPtr_->loadRHS();
      }

      // Set a constant history.
      anaManagerRCPtr_->getTIADataStore()->setConstantHistory();
      anaManagerRCPtr_->getTIADataStore()->computeDividedDifferences();
      wimRCPtr_->obtainCorrectorDeriv();

#ifdef Xyce_PARALLEL_MPI
      // Update vectors with off proc values.
      lasSystemRCPtr_->updateExternValsSolnVector(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr);
      lasSystemRCPtr_->updateExternValsSolnVector(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr);
      lasSystemRCPtr_->updateExternValsStateVector(anaManagerRCPtr_->getTIADataStore()->nextStatePtr);
      lasSystemRCPtr_->updateExternValsStateVector(anaManagerRCPtr_->getTIADataStore()->currStatePtr);
      lasSystemRCPtr_->updateExternValsStoreVector(anaManagerRCPtr_->getTIADataStore()->nextStorePtr);
      lasSystemRCPtr_->updateExternValsStoreVector(anaManagerRCPtr_->getTIADataStore()->currStorePtr);
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
      // this is really a defect of a paused simulation's Transient object
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
      dout() << "  transient loop called with resume true " << std::endl;
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
        anaManagerRCPtr_->getTIADataStore()->currSolutionPtr, anaManagerRCPtr_->getTIADataStore()->currStatePtr, anaManagerRCPtr_->getTIADataStore()->currStorePtr);
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
// Function      : Transient::loopProcess()
// Purpose       : Conduct the time stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::loopProcess()
{
  bool bsuccess = true;

  // Transient time stepping loop:
  while (!(secRCPtr_->finished()))
  {
#ifdef Xyce_VERBOSE_TIME
    printStepHeader(Xyce::lout());
#endif
    printProgress(Xyce::lout());

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
      dout() << std::endl;
      dout() << "Transient::loopProcess()" << std::endl;
      dout() << "beginningIntegration = " << beginningIntegration << std::endl;
      dout() << "secRCPtr_->stepAttemptStatus = " << secRCPtr_->stepAttemptStatus << std::endl;
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
      loaderRCPtr_->getInitialQnorm(anaManagerRCPtr_->getTIADataStore()->innerErrorInfoVec);
      double suggestedMaxTime=0.0;
      if( maxTimeStepExpressionGiven_ )
      {
        suggestedMaxTime = maxTimeStepExpressionRCPtr_->evaluate(
          anaManagerRCPtr_->getTIADataStore()->currSolutionPtr, anaManagerRCPtr_->getTIADataStore()->currStatePtr, anaManagerRCPtr_->getTIADataStore()->currStorePtr);
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
      secRCPtr_->outputTimeInfo(lout());
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
          std::string msg = "Transient::loopProcess() Nonlinear solver stalled. Calling this a pass";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
          processSuccessfulStep();
        }
        // another VERY dangerous options.
        // if the non-linear solver is reporting too big of an update, and we're very close to a min
        // time step, then calls this failure a pass
        if( secRCPtr_->newtonConvergenceStatus == -2)
        {
          std::string msg = "Transient::loopProcess() Update too big. Calling this a pass";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
          processSuccessfulStep();
        }
        // another VERY dangerous options.
        // if the non-linear solver is not converging in the max number of steps,
        // and we're very close to a min time step, then calls this failure a pass
        //
        // if( secRCPtr_->newtonConvergenceStatus == -1)
        // {
        //   std::string msg = "Transient::loopProcess() Too many steps. Calling this a pass";
        //   N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
        //   processSuccessfulStep();
        // }
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
      dout() << std::endl;
      dout() << "   Here we are, just before checking whether to pause. " << std::endl;
      dout() << "   minTimeStep = " << secRCPtr_->minTimeStep << std::endl;
      dout() << "   final time = " << tiaParams.finalTime << std::endl;
      dout() << "   pause time = " << anaManagerRCPtr_->getPauseTime() << std::endl;
      dout() << "   initial time = " << secRCPtr_->initialTime << std::endl;
      dout() << "   current time = " << secRCPtr_->currentTime << std::endl;
      if (anaManagerRCPtr_->getPauseTime() == secRCPtr_->currentTime)
      {
        dout() << "    Pause time and current time equal " << std::endl;
      }
      else
      {
        dout() << "     difference between current and pause times is "
          << anaManagerRCPtr_->getPauseTime() - secRCPtr_->currentTime << std::endl;
      }
      if (anaManagerRCPtr_->getPauseTime() == secRCPtr_->initialTime)
      {
        dout() << "    Pause time and initial time equal " << std::endl;
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
        dout() << "Transient::loopProcess():   pausing simulation " << std::endl;
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
      lout() << "Exit time exceeded.  Exiting transient loop\n" << std::endl;
      bsuccess = true;
      break;
    }

    if (tiaParams.exitStep != -1 &&
          static_cast<int>(stepNumber) == tiaParams.exitStep)
    {
      lout() << "Exit step.  Exiting transient loop\n" << std::endl;
      bsuccess = true;
      break;
    }

  } // end of time loop
  endTRANtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
  tranStats = saveLoopInfo();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Transient::mixedSignalStep
// Purpose       :
// Special Notes : Habanero API function
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
bool Transient::mixedSignalStep()
{
  takeAnIntegrationStep_();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Transient::preStepDetails
// Purpose       :
// Special Notes : Habanero API function
// Scope         : private_
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
void Transient::preStepDetails (double maxTimeStepFromHabanero)
{
#ifdef Xyce_VERBOSE_TIME
  printStepHeader(Xyce::lout());
#endif
  printProgress(Xyce::lout());

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
    dout() << std::endl;
    dout() << "Transient::loopProcess()" << std::endl;
    dout() << "beginningIntegration = " << beginningIntegration << std::endl;
    dout() << "secRCPtr_->stepAttemptStatus = " << secRCPtr_->stepAttemptStatus << std::endl;
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
    loaderRCPtr_->getInitialQnorm(anaManagerRCPtr_->getTIADataStore()->innerErrorInfoVec);
    double suggestedMaxTime=0.0;
    if( maxTimeStepExpressionGiven_ )
    {
      suggestedMaxTime = maxTimeStepExpressionRCPtr_->evaluate(
        anaManagerRCPtr_->getTIADataStore()->currSolutionPtr, anaManagerRCPtr_->getTIADataStore()->currStatePtr, anaManagerRCPtr_->getTIADataStore()->currStorePtr);
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
    secRCPtr_->outputTimeInfo(lout());
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
// Function      : Transient::finalizeStep
// Purpose       :
// Special Notes : Habanero API function
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/04/09
//-----------------------------------------------------------------------------
bool Transient::finalizeStep ()
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
        std::string msg = "Transient::loopProcess() Nonlinear solver stalled. Calling this a pass";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
        processSuccessfulStep();
      }
      // another VERY dangerous options.
      // if the non-linear solver is reporting too big of an update, and we're very close to a min
      // time step, then calls this failure a pass
      if( secRCPtr_->newtonConvergenceStatus == -2)
      {
        std::string msg = "Transient::loopProcess() Update too big. Calling this a pass";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
        processSuccessfulStep();
      }
      // another VERY dangerous options.
      // if the non-linear solver is not converging in the max number of steps,
      // and we're very close to a min time step, then calls this failure a pass
      // 
      // if( secRCPtr_->newtonConvergenceStatus == -1)
      // {
      //   std::string msg = "Transient::loopProcess() Too many steps. Calling this a pass";
      //   N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
      //   processSuccessfulStep();
      // }
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
    dout() << std::endl;
    dout() << "   Here we are, just before checking whether to pause. " << std::endl;
    dout() << "   minTimeStep = " << secRCPtr_->minTimeStep << std::endl;
    dout() << "   final time = " << tiaParams.finalTime << std::endl;
    dout() << "   pause time = " << anaManagerRCPtr_->getPauseTime() << std::endl;
    dout() << "   initial time = " << secRCPtr_->initialTime << std::endl;
    dout() << "   current time = " << secRCPtr_->currentTime << std::endl;
    if (anaManagerRCPtr_->getPauseTime() == secRCPtr_->currentTime)
    {
      dout() << "    Pause time and current time equal " << std::endl;
    }
    else
    {
      dout() << "     difference between current and pause times is " << anaManagerRCPtr_->getPauseTime() - secRCPtr_->currentTime << std::endl;
    }
    if (anaManagerRCPtr_->getPauseTime() == secRCPtr_->initialTime)
    {
      dout() << "    Pause time and initial time equal " << std::endl;
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
      dout() << "Transient::loopProcess():   pausing simulation " << std::endl;
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
    
    lout() << "Exit time exceeded.  Exiting transient loop\n" << std::endl;
    recoverableFailureFlag = false;
  }

  if (tiaParams.exitStep != -1 &&
        static_cast<int>(stepNumber) == tiaParams.exitStep)
  {
    lout() <<"Exit step.  Exiting transient loop\n" << std::endl;
    recoverableFailureFlag = false;
  }

  return recoverableFailureFlag;
}

//-----------------------------------------------------------------------------
// Function      : Transient::processSuccessfulDCOP()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::processSuccessfulDCOP()
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

  anaManagerRCPtr_->getTIADataStore()->setConstantHistory();
  anaManagerRCPtr_->getTIADataStore()->computeDividedDifferences();
  wimRCPtr_->obtainCorrectorDeriv();

  anaManagerRCPtr_->getTIADataStore()->updateSolDataArrays ();

  tranopOutputs ();

  // Now that output has been called, update the doubleDCOP step
  // if neccessary. (Only matters for pde problems)
  doubleDCOPStep_ = tiaParams.lastDCOPStep;

  //Test and save restart if necessary
  if( anaManagerRCPtr_->testRestartSaveTime_() )
  {
    if (DEBUG_RESTART)
      dout() << "\n " << commandLine_.getArgumentValue("netlist")
                   << "  Calling dumpRestartData" << std::endl;

    restartMgrRCPtr_->dumpRestartData( secRCPtr_->currentTime );

    if (DEBUG_RESTART)
      dout() << "  Done Calling dumpRestartData" << std::endl;
  }

  // This output call is for device-specific output, such as .OP,
  // or internal plot output from PDE(TCAD) devices.
  loaderRCPtr_->output();

  nlsMgrRCPtr_->allocateTranSolver();
  secRCPtr_->previousCallStepSuccessful = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::processSuccessfulStep()
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
    dout() << "  Transient::processSuccessfulStep()" << std::endl
           << "Newton step succeeded:" << std::endl
           << "nextSolutionPtr: " << std::endl;
    
    anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr->printPetraObject(dout());
    dout() << std::endl;
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
      anaManagerRCPtr_->getTIADataStore()->currSolutionPtr, anaManagerRCPtr_->getTIADataStore()->currStatePtr, anaManagerRCPtr_->getTIADataStore()->currStorePtr);
  }
  secRCPtr_->updateMaxTimeStep( suggestedMaxTime );
  secRCPtr_->updateMinTimeStep();
  secRCPtr_->updateBreakPoints();

#ifdef Xyce_VERBOSE_TIME
  if (tiaParams.debugLevel > 0)
  {
    dout() << "Transient Analysis:  accepting time step" << std::endl;
  }
#endif // Xyce_VERBOSE_TIME

  wimRCPtr_->completeStep();

#ifdef Xyce_VERBOSE_TIME
  if (tiaParams.errorAnalysisOption == 1)
  {
    dout().precision(15);
    dout() << "ERROROPTION=1: TimeStepLimitedbyBP = " << tiaParams.TimeStepLimitedbyBP << "\n" << std::endl;
    dout() << "ERROROPTION=1: NL Its =  " << secRCPtr_->nIterations << "\n" << std::endl;
    dout() << "ERROROPTION=1: New DeltaT = " << secRCPtr_->currentTimeStep << "\n" << std::endl;
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
      dout() << " Checking whether to set breakpointrestartstep" << std::endl;
      dout() << "   current - stop  = " << timeDiff1 << std::endl;
      dout() << "   current - final = " << timeDiff2 << std::endl;
      dout() << "   bpTol           = " << bpTol << std::endl;
      if (timeDiff1 <= bpTol && timeDiff2 > bpTol)
        dout() << "    setting breakPointRestartStep to " << tranStepNumber;
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
    anaManagerRCPtr_->getTIADataStore()->timeSteps.push_back(currentTime);
    anaManagerRCPtr_->getTIADataStore()->timeStepsBreakpointFlag.push_back(beginningIntegration);
    N_LAS_Vector * aVecPtr = new N_LAS_Vector( *(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr) );
    anaManagerRCPtr_->getTIADataStore()->fastTimeSolutionVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(anaManagerRCPtr_->getTIADataStore()->currStatePtr) );
    anaManagerRCPtr_->getTIADataStore()->fastTimeStateVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(anaManagerRCPtr_->getTIADataStore()->daeQVectorPtr) );
    anaManagerRCPtr_->getTIADataStore()->fastTimeQVec.push_back( aVecPtr );
    aVecPtr = new N_LAS_Vector( *(anaManagerRCPtr_->getTIADataStore()->currStorePtr) );
    anaManagerRCPtr_->getTIADataStore()->fastTimeStoreVec.push_back( aVecPtr );
  }

  // 03/16/04 tscoffe:  This is where the solution pointers are rotated.
  anaManagerRCPtr_->getTIADataStore()->updateSolDataArrays();

#ifdef Xyce_DEBUG_ANALYSIS
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
    anaManagerRCPtr_->getTIADataStore()->outputSolDataArrays(Xyce::dout());
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
// Function      : Transient::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::processFailedStep()
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
  dout() << "Transient Analysis:  rejecting time step" << std::endl;
#endif // Xyce_VERBOSE_TIME

  wimRCPtr_->rejectStep();

#ifdef Xyce_VERBOSE_TIME
  if (tiaParams.errorAnalysisOption == 1)
  {
    dout().precision(15);
    dout() << "ERROROPTION=1: TimeStepLimitedbyBP = " << tiaParams.TimeStepLimitedbyBP << "\n" << std::endl;
    dout() << "ERROROPTION=1: NL Its =  " << secRCPtr_->nIterations << "\n" << std::endl;
    dout() << "ERROROPTION=1: New DeltaT = " << secRCPtr_->currentTimeStep << "\n" << std::endl;
  }
#endif // Xyce_VERBOSE_TIME

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    dout() << "  Transient::processFailedStep" << std::endl
           << "Newton step failed:" << std::endl
           << "nextSolutionPtr: " << std::endl;
    anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr->printPetraObject(dout());
    dout() << std::endl;
  }
#endif
  totalNumberFailedStepsAttempted_  += 1;
  secRCPtr_->numberSuccessiveFailures += 1;

  if (secRCPtr_->currentTimeStep <= secRCPtr_->minTimeStep)
  {
#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      dout() << "currTimeStep: " << secRCPtr_->currentTimeStep << std::endl;
      dout() << "minTimeStep:  " << secRCPtr_->minTimeStep << std::endl;
    }
#endif

    // before we exit with a time step too small error, check if minTimeStepRecoveryCounter is greater than zero
    // if it is, then the user wants us to try to accept the step that had the minimum
    // estimated error over tol.
    if( tiaParams.minTimeStepRecoveryCounter > 0 )
    {
      lout() << "Attempting to retake and accept step where estimated error over tolerance was: " << minEstErrorOverTol
             << " and time step was: " << timeStepAtMinEstErrorOverTol << std::endl;

      tiaParams.minTimeStepRecoveryCounter--;
      bsuccess = retakeAndAcceptTimeStep( timeStepAtMinEstErrorOverTol );
    }
    else
    {
      outputQueuedData();
      lout() << "Time step too small near step number: " <<  stepNumber << "  Exiting transient loop.\n" << std::endl;

      bsuccess = false;
    }
  }

  if (tiaParams.constantStepSize)
  {
    outputQueuedData();
    lout() << "Newton solver failed in constant time step mode.  Exiting transient loop.\n" << std::endl;

    bsuccess = false;
  }

  if (tiaParams.exitStep != -1 &&
          static_cast<int>(totalNumberSuccessStepsThisParameter_) == (tiaParams.exitStep-1))
  {
    outputQueuedData();
    lout() << "Exit Step.  Exiting transient loop\n" << std::endl;
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
// Function      : Transient::processFailedDCOP
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::processFailedDCOP()
{
  bool bsuccess = true;

  loaderRCPtr_->stepFailure (anaManagerRCPtr_->currentMode_);

  // DC Operating Point failed.
  (totalNumberFailedStepsAttempted_)++;
  (secRCPtr_->numberSuccessiveFailures)++;

  lout() << "DC Operating Point Failed.  Exiting transient loop" << std::endl;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::finish()
{
  bool bsuccess = true;

  if (tiaParams.saveTimeStepsFlag && anaManagerRCPtr_->getHBFlag())
  {
    RefCountPtr<N_TIA_DataStore> dsPtr_ = anaManagerRCPtr_->getTIADataStore();
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
    dout() << "Calling finishOutput" << std::endl;
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
// Function      : Transient::handlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool Transient::handlePredictor()
{
  anaManagerRCPtr_->getTIADataStore()->setErrorWtVector();
  wimRCPtr_->obtainPredictor();
  wimRCPtr_->obtainPredictorDeriv();

#ifdef Xyce_DEBUG_ANALYSIS
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
  {
    dout() << "  Transient::handlePredictor" << std::endl;
    anaManagerRCPtr_->getTIADataStore()->outputPredictedSolution(Xyce::dout());
    anaManagerRCPtr_->getTIADataStore()->outputPredictedDerivative(Xyce::dout());
  }
#endif
#endif

  // Now, in case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  loaderRCPtr_->startTimeStep ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Transient::resetForStepAnalysis()
// Purpose       : When doing a .STEP sweep, some data must be reset to its
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool Transient::resetForStepAnalysis()
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
// Function      : Transient::resetForHB()
// Purpose       : When doing initial transient run, some analyses require
//                 a reset function
//                 they can fill this in if needed.
// Special Notes :
// Scope         : public
// Creator       : T. Mei, SNL, Parallel Computational Sciences
// Creation Date : 2/26/09
//-----------------------------------------------------------------------------
bool Transient::resetForHB()
{
  nlsMgrRCPtr_->resetAll(DC_OP);
  secRCPtr_->resetAll();
  anaManagerRCPtr_->getTIADataStore()->resetAll();
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
// Function      : Transient::finalVerboseOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool Transient::finalVerboseOutput()
{
  bool bsuccess = true;

  lout() << "***** Problem read in and set up time: " << anaManagerRCPtr_->solverStartTime_ << " seconds" << std::endl;

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
    lout() << " ***** DCOP time: " << time << " seconds.  Breakdown follows:" << std::endl;

    printLoopInfo(0, dcStats);
  }

  if (anaManagerRCPtr_->analysis == ANP_MODE_TRANSIENT &&
      endTRANtime >= anaManagerRCPtr_->startTRANtime)
  {
    lout() << " ***** Transient Stepping time: " << endTRANtime-anaManagerRCPtr_->startTRANtime << " seconds.  Breakdown follows:" << std::endl;

    printLoopInfo(dcStats, tranStats);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::outputQueuedData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
void Transient::outputQueuedData()
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

    lout() << " *** Transient failure history: " << std::endl;
    if (tiaParams.errorAnalysisOption == 1)
    {
      // truncation error is not used here so est err over tol is not useful in the output
      lout() << "Time        Time      Step      Non-Linear Solver      node    node" << std::endl;
      lout() << "(sec)       Step     Status   Status   Iters   ||F||   index   name" << std::endl;
    }
    else
    {
      lout() << "Time        Time      Step   EstErr      Non-Linear Solver      node     node" << std::endl;
      lout() << "(sec)       Step     Status  OverTol   Status  Iters   ||F||    index    name" << std::endl;
    }

    for( int i=0; i<queueSize_; i++ )
    {
      int fieldWidth=10;
      lout() << std::scientific << std::setprecision(fieldWidth-7) << std::setfill(' ') << std::right << std::setw( fieldWidth )
        << timeQueue_.at_from_tail(i) << "  "
        << timeStepQueue_.at_from_tail(i) << "  ";
      if( stepStatusQueue_.at_from_tail(i) == 1 )
      {
        lout() << "pass  ";
      }
      else
      {
        lout() << "fail  ";
      }
      if (tiaParams.errorAnalysisOption == 1)
      {
      }
      else
      {
        lout() << estErrorOverTolQueue_.at_from_tail(i) << "  ";
      }
      int nlStatus = nonlinearSolverStatusQueue_.at_from_tail(i);
      lout() << std::setw(7) << std::right;
      if( nlStatus == nlReturnCodes.normTooSmall )
      {
        lout() << "P:s nrm";
      }
      else if( nlStatus == nlReturnCodes.normalConvergence)
      {
        lout() << "pass   ";
      }
      else if( nlStatus == nlReturnCodes.nearConvergence )
      {
        lout() << "P:near ";
      }
      else if( nlStatus == nlReturnCodes.smallUpdate )
      {
        lout() << "P:s up ";
      }
      else if( nlStatus == nlReturnCodes.nanFail )
      {
        lout() << "F:NaN  ";
      }
      else if( nlStatus == nlReturnCodes.tooManySteps )
      {
        lout() << "F:max s";
      }
      else if( nlStatus == nlReturnCodes.tooManyTranSteps )
      {
        lout() << "F:max s";
      }
      else if( nlStatus == nlReturnCodes.updateTooBig )
      {
        lout() << "F:big u";
      }
      else if( nlStatus == nlReturnCodes.stalled )
      {
        lout() << "F:stall";
      }
      else if( nlStatus == nlReturnCodes.wrmsExactZero )
      {
        lout() << "F:n zro";
      }
      else if( nlStatus == nlReturnCodes.innerSolveFailed )
      {
        lout() << "F:in Fl";
      }
      else
      {
        lout() << "code=" <<
            nonlinearSolverStatusQueue_.at_from_tail(i) << "  ";
      }

      lout() << std::right << std::setw( 4 )
        << nonlinearSolverNumIterationsQueue_.at_from_tail(i) << "  "
        << nonlinearSolverMaxNormQueue_.at_from_tail(i) ;

      int outIndex = nonlinearSolverMaxNormIndexQueue_.at_from_tail(i) ;
      lout() << std::right << std::fixed << std::setw( 7 ) << outIndex;

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
      lout() << std::left << "    " << outIndexName;
      lout() << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::takeAnIntegrationStep_
// Purpose       : Take a transient integration step.
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void Transient::takeAnIntegrationStep_()
{
  handlePredictor();
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();

#ifdef Xyce_DEBUG_ANALYSIS
  dout() << "Transient::takeAnIntegrationStep_:  newtonConvergenceStatus = " <<  secRCPtr_->newtonConvergenceStatus << std::endl;
#endif

  anaManagerRCPtr_->getTIADataStore()->stepLinearCombo ();
  gatherStepStatistics_ ();
  secRCPtr_->nIterations = nlsMgrRCPtr_->getNumIterations();
  secRCPtr_->evaluateStepError ();

  return;
}


//-----------------------------------------------------------------------------
// Function      : Transient::twoLevelStep
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
bool Transient::twoLevelStep ()
{
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();
  anaManagerRCPtr_->getTIADataStore()->stepLinearCombo ();

  gatherStepStatistics_ ();
  secRCPtr_->evaluateStepError ();

  endTRANtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
  tranStats = saveLoopInfo();

  return (secRCPtr_->stepAttemptStatus);
}

//-----------------------------------------------------------------------------
// Function      : Transient::retakeAndAcceptTimeStep
// Purpose       : Do a requested time step and accept it
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 01/23/09
//-----------------------------------------------------------------------------
bool Transient::retakeAndAcceptTimeStep( double aTimeStep )
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
  // At this point, Transient::processFailedStep() has already determined
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
    lout() << "Time step too small near step number: " <<  stepNumber << "  Exiting transient loop.\n" << std::endl;
    bsuccess = false;
  }
  else
  {
    processSuccessfulStep();
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::printStepHeader()
// Purpose       : Prints out time step information.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void Transient::printStepHeader(std::ostream &os)
{
  if (VERBOSE_TIME) {
    os << "***** " << (DEBUG_ANALYSIS ? commandLine_.getArgumentValue("netlist") : "") << "  ";

    if (dcopFlag_)
    {
      os << "Start of DCOP STEP                        # ";
    }
    else
    {
      if (beginningIntegration)
      {
        if (secRCPtr_->currentTime == tiaParams.initialTime)
        {
          os << "Start of Time Step (INITIAL STEP)         # ";
        }
        else
        {
          os << "Start of Time Step (DISCONTINUITY STEP)   # ";
        }
      }
      else
      {
        if (!secRCPtr_->stepAttemptStatus)
        {
          os << "Start of Time Step (AFTER FAILED STEP)    # ";
        }
        else
        {
          os <<  "Start of Time Step (AFTER SUCCESS STEP)   # ";
        }
      }
    }

    os << totalNumberSuccessStepsThisParameter_ + 1 << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::printProgress()
// Purpose       : Outputs run completion percentage and estimated
//                 time-to-completion.
//
// Special Notes : This will need some fixing to work with .STEP.
//
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 06/07/2002
//-----------------------------------------------------------------------------
void Transient::printProgress(std::ostream &os)
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

    // Report the beginning of the DC OP calculation.  First call in OP.
    if (dcopFlag_)
    {
      startDCOPtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
      if (gui_)
      {
        Report::signalProgress("Xyce in DC Operating Point Calculation");
      }
      os << "***** Beginning DC Operating Point Calculation...\n" << std::endl;
    }
    else if (firstTime && totalNumberSuccessStepsThisParameter_ == 1)
    {
      anaManagerRCPtr_->startTRANtime = anaManagerRCPtr_->xyceTranTimerPtr_->elapsedTime();
      dcStats = saveLoopInfo();
      firstTime = false;
      if (gui_)
      {
        Report::signalProgress("Xyce in Transient Calculation");
      }
      os << "***** Beginning Transient Calculation...\n" << std::endl;
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
          os << "***** Percent complete: " <<  percentComplete <<  " %" << std::endl;
        }

        if (estCompletionTime > N_UTL_MachineDependentParams::MachineEpsilon())
        {
          unsigned int days, hours, minutes, seconds;
          days    = static_cast<int> (estCompletionTime / 86400);
          hours   = static_cast<int> ((estCompletionTime - days * 86400) / 3600);
          minutes = static_cast<int> ((estCompletionTime - days * 86400 - hours * 3600) / 60);
          seconds = static_cast<int> (estCompletionTime - days * 86400 - hours * 3600 - minutes * 60);

          char timeStr[256];
          for (char *c = timeStr; c != timeStr + sizeof(timeStr); ++c)
            *c = 0;

#ifdef Xyce_PARALLEL_MPI
          // get current local system time
          time_t t = time( NULL );
          struct tm * now = localtime( &t );

          // format and display output
          if ( ( t != (time_t)-1 ) && ( strftime( timeStr, 255, "%c", now ) != 0 ) )
          {
            os << "***** Current system time: " << timeStr << std::endl;
          }

          else
          {
            os << "***** Current system time could not be determined." << std::endl;
          }
#endif

          if (days > 0)
            sprintf(timeStr, "%3d days, %2d hrs., %2d min., %2d sec.", days, hours, minutes, seconds);
          else if (hours > 0)
            sprintf(timeStr, "%2d hrs., %2d min., %2d sec.", hours, minutes, seconds);
          else if (minutes > 0)
            sprintf(timeStr, "%2d min., %2d sec.", minutes, seconds);
          else
            sprintf(timeStr, "%2d sec.", seconds);

          if (gui_)
          {
            std::ostringstream ost;
            ost << "Xyce transient ";
            if (percentComplete < 10)
              ost.precision(2);
            else
              ost.precision(3);
            ost << percentComplete
                << "%% complete ... Estimated time to completion: "
                << timeStr << std::endl;
            Report::signalProgress(ost.str());
          }
          else
          {
            os << "***** Estimated time to completion: " << timeStr << std::endl << std::endl;
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::noopOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void Transient::noopOutputs ()
{
  if( testOutputTime_() )
  {
    if ( !firstDoubleDCOPStep_ () )
    {
      bool printIC = true;
      if (printIC)
      {
        outputMgrAdapterRCPtr_->tranOutput(
            secRCPtr_->currentTime, 
            *(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr), 
            *(anaManagerRCPtr_->getTIADataStore()->currStatePtr), 
            *(anaManagerRCPtr_->getTIADataStore()->currStorePtr),
            objectiveVec_, 
            dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
      }
      if ( ( !Teuchos::is_null(anaManagerRCPtr_->mpdeMgrPtr_)  )  &&
           (  anaManagerRCPtr_->mpdeMgrPtr_->getMPDEFlag()     )  &&
           ( !(anaManagerRCPtr_->mpdeMgrPtr_->getMPDEIcFlag()) )
         )
      {
        // output the operating point if we can in MPDE mode
        outputMgrAdapterRCPtr_->outputMPDE(secRCPtr_->currentTime, 
            *(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr));
      }
    }
    updateOutputTime_(secRCPtr_->currentTime);
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::tranopOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void Transient::tranopOutputs ()
{
  //Test and output if necessary
  if( testOutputTime_() )
  {
    // Make sure this isn't the NLP step of a PDE DCOP.
    if ( !firstDoubleDCOPStep_ () )
    {
#ifdef Xyce_DEBUG_ANALYSIS
      dout() << "Calling conventional TRANOP outputs!" << std::endl;
#endif
      outputMgrAdapterRCPtr_->tranOutput(secRCPtr_->currentTime, 
          *(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr), 
          *(anaManagerRCPtr_->getTIADataStore()->currStatePtr), 
          *(anaManagerRCPtr_->getTIADataStore()->currStorePtr),
          objectiveVec_, 
            dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

      outputMgrAdapterRCPtr_->outputMPDE(secRCPtr_->currentTime, 
          *(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr));
      // check this - may be a mistake.
    }
    updateOutputTime_(secRCPtr_->currentTime);
  }

  // SAVE and DCOP restart:
  if ( anaManagerRCPtr_->testDCOPOutputTime_() || anaManagerRCPtr_->testSaveOutputTime_() )
  {
    outputMgrAdapterRCPtr_->outputDCOP( *(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr) );
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::tranStepOutputs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/31/2007
//-----------------------------------------------------------------------------
void Transient::tranStepOutputs ()
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
      dout() << "doNotInterpolate is TRUE!" << std::endl;
    }
    else
    {
      dout() << "doNotInterpolate is FALSE!" << std::endl;
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
        dout() << "Calling MPDE outputs!" << std::endl;
      }
#endif
      outputMgrAdapterRCPtr_->outputMPDE(secRCPtr_->currentTime, *(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr));

      // If we are actually in the MPDE phase, rather than the initial
      // condition, then output here.
      if (anaManagerRCPtr_->mpdeMgrPtr_->getMPDEFlag()==true && tiaParams.outputInterpMPDE)
      {
        if (anaManagerRCPtr_->mpdeMgrPtr_->getWaMPDEFlag()==false)
        {
          std::vector<double> fastTimes = anaManagerRCPtr_->mpdeMgrPtr_->getFastTimePoints();
          wimRCPtr_->printMPDEOutputSolution(
                  outputMgrAdapterRCPtr_, secRCPtr_->currentTime,
                  anaManagerRCPtr_->getTIADataStore()->currSolutionPtr,
                  fastTimes );
        }
        else
        {
          std::vector<double> fastTimes = anaManagerRCPtr_->mpdeMgrPtr_->getFastTimePoints();
          int phiGID = anaManagerRCPtr_->mpdeMgrPtr_->getPhiGID();
          wimRCPtr_->printWaMPDEOutputSolution(
                  outputMgrAdapterRCPtr_, secRCPtr_->currentTime,
                  anaManagerRCPtr_->getTIADataStore()->currSolutionPtr,
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
              anaManagerRCPtr_->getTIADataStore()->currSolutionPtr,
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
              anaManagerRCPtr_->getTIADataStore()->currSolutionPtr,
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
              anaManagerRCPtr_->getTIADataStore()->currSolutionPtr,
              doNotInterpolate,
              outputInterpolationTimes_,
              true) ;
  }

  // SAVE output:
  if ( anaManagerRCPtr_->testSaveOutputTime_() )
  {
    wimRCPtr_->saveOutputSolution(
                outputMgrAdapterRCPtr_,
                anaManagerRCPtr_->getTIADataStore()->currSolutionPtr,
                anaManagerRCPtr_->saveTime_, doNotInterpolate);
  }
}

//-----------------------------------------------------------------------------
// Function      : Transient::testOutputTime_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel ComputationalSciences.
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool Transient::testOutputTime_()
{
  bool flag;

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    dout().width(21); dout().precision(13); dout().setf(std::ios::scientific);
    dout() << "Transient::testOutputTime.  secPtr_->currentTime = " << secRCPtr_->currentTime << std::endl;
    dout() << "Transient::testOutputTime.  tStart      = " << tiaParams.tStart << std::endl;
    dout() << "Transient::testOutputTime.  nextOutputTime_ = " << anaManagerRCPtr_->nextOutputTime_ << std::endl;
    for (int i=0;i<anaManagerRCPtr_->outputIntervals_.size();++i)
    {
      dout() << "outputIntervals_["<<i<<"].first = " << anaManagerRCPtr_->outputIntervals_[i].first;
      dout() << std::endl;
    }
    dout() << std::endl;
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
// Function      : Transient::updateOutputTime_
// Purpose       : Advance output time so it is the next one after currTime
// Special Notes : formerly part of testOutputTime_
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 02/14/07
//-----------------------------------------------------------------------------
void Transient::updateOutputTime_(double currTime)
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
      std::pair<double, double> currInterval, nextInterval;
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
// Function      : Transient::computeOutputInterpolationTimes_
// Purpose       : When we pass "nextOutputTime_", we might have skipped
//                 over points where output was requested.  Make a list of
//                 those times so we can interpolate to them.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling.
// Creation Date : 02/14/07
//-----------------------------------------------------------------------------
void Transient::computeOutputInterpolationTimes_(double currTime)
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

} // namespace Analysis
} // namespace Xyce
