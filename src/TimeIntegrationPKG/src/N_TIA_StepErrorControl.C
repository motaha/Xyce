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
// Filename      : $RCSfile: N_TIA_StepErrorControl.C,v $
//
// Purpose       : This file contains the functions which define the
//		             time integration stepsize control algorithm.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.214.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <iostream>
#include <stdio.h>

#include <sstream>
 
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

#include <N_ANP_AnalysisManager.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntInfo.h>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_Functors.h>

#include <N_PDS_Comm.h>

#include <N_ERH_ErrorMgr.h>

#include <N_IO_CmdParse.h>

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::N_TIA_StepErrorControl
// Purpose       : Non-argument constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/06/01
//-----------------------------------------------------------------------------
N_TIA_StepErrorControl::N_TIA_StepErrorControl(
    N_IO_CmdParse & cp,
    N_ANP_AnalysisManager & anaManager,
    N_TIA_TIAParams & tiaP)
  :
    startingTimeStep(1.0e-10),
    currentTimeStep(1.0e-10),
    lastAttemptedTimeStep(1.0e-10),
    lastTimeStep(1.0e-10),
    minTimeStep(0.0),
    maxTimeStep(0.0),
    maxTimeStepUser(1.0e+99),
    maxTimeStepBP(0.0),
    savedTimeStep(1.0e-10),
    lastTime(0.0),
    currentTime(0.0),
    nextTime(0.0),
    stopTime(0.0),
    initialTime(0.0),
    finalTime(0.0),
    currentTimeStepRatio(0.0),
    currentTimeStepSum(0.0),
    lastTimeStepRatio(0.0),
    lastTimeStepSum(0.0),
    newtonConvergenceStatus(-1),
    nIterations(0),
    numberSuccessiveFailures(0),
    stepAttemptStatus(true),
    previousCallStepSuccessful(false),
    estOverTol_(0.0),
    initializeFlag_(false),
    minStepPrecisionFac_(10.0),
    newtonStepReduction_(0.25),
    restartTimeStepScale_(0.005),
    tolAimFac_(0.5),
    tiaParams_(tiaP),
    anaManager_(anaManager),
    wimPtr_(NULL),
    loader_( *(anaManager_.loaderPtr) ),
    commandLine_(cp),
    // define "heuristic" StepSize and Error Control parameters.

    // new-DAE variables:
    currentOrder_(1), // Current order of integration
    oldOrder_(1), // previous order of integration
    minOrder_(1),  // minimum order = max(1,user option minord) - see below.
    maxOrder_(5),  // maximum order = min(5,user option maxord) - see below.
    usedOrder_(1),  // order used in current step (used after currentOrder_ is updated)
    alphas_(-1.0),  // $\alpha_s$ fixed-leading coefficient of this BDF method
    alpha_(6,0.0),  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
                    // note:   $h_n$ = current step size, n = current time step
    alpha0_(0.0),   // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
    cj_ (0.0),      // $-\alpha_s/h_n$ coefficient used in local error test
    ck_ (0.0),      // local error coefficient gamma_[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
    sigma_(6,0.0),  // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
    gamma_(6,0.0),  // calculate time derivative of history array for predictor 
    beta_(6,0.0),   // coefficients used to evaluate predictor from history array
    psi_(6,0.0),    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
                    // compute $\beta_j(n)$
    numberOfSteps_(0),   // number of total time integration steps taken
    nef_(0),
    usedStep_(0.0),
    nscsco_(0),
    Ek_(0.0),
    Ekm1_(0.0),
    Ekm2_(0.0),
    Ekp1_(0.0),
    Est_(0.0),
    Tk_(0.0),
    Tkm1_(0.0),
    Tkm2_(0.0),
    Tkp1_(0.0),
    newOrder_(1),
    initialPhase_(true),
    h0_safety_(2.0),
    h0_max_factor_(0.005),  // New value, to match old-DAE.
    //h0_max_factor_(0.0001), // this is the new-DAE equivalent of restartTimeStepScale_(0.005)
    h_phase0_incr_(2.0),
    h_max_inv_(0.0),
    Tkm1_Tk_safety_(2.0),
    Tkp1_Tk_safety_(0.5),
#ifndef Xyce_USE_Q_NORM
    r_factor_(1.0),
#else
    r_factor_(0.9),
#endif
    r_safety_(2.0),
    r_fudge_(0.0001),
    //r_min_(0.125),
    r_min_(0.25),  // r_min_ is the same as the old-DAE variable, minFailStepFac_.
    r_max_(0.9),   // r_max_ is the same as the old-DAE variable, maxFailStepFac_.
    r_hincr_test_(2.0),
    r_hincr_(2.0),
    max_LET_fail_(15)
{
  // make sure we initialize the iterator currentPauseBP to a valid but 
  // bad value until it is calulcated
  currentPauseBP = breakPoints_.end();

  setTIAParams();
  setBreakPoint(N_UTL_BreakPoint(tiaParams_.finalTime,PAUSE_BREAKPOINT));

  h0_max_factor_ = tiaParams_.restartTimeStepScale;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::~N_TIA_StepErrorControl
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------

N_TIA_StepErrorControl::~N_TIA_StepErrorControl()
{
  return ;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::setTIAParams
// Purpose       : This function copies stuff out of the tiaParams object into
//                 variables that are local to the step error control class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 09/06/01
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::setTIAParams()
{
  startingTimeStep = tiaParams_.userSpecified_startingTimeStep;
  currentTimeStep  = tiaParams_.userSpecified_startingTimeStep;
  initialTime      = tiaParams_.initialTime;
  finalTime        = tiaParams_.finalTime;
  currentTime      = tiaParams_.initialTime;
  nextTime         = tiaParams_.initialTime;
  lastTime         = tiaParams_.initialTime;

  if (tiaParams_.maxTimeStepGiven)
  {
    maxTimeStepUser  = tiaParams_.maxTimeStep;
    maxTimeStep      = tiaParams_.maxTimeStep;
  }
  else
  {
    maxTimeStep  = 0.1*(tiaParams_.finalTime-tiaParams_.initialTime);
  }

  initializeFlag_  = true;

  restartTimeStepScale_ = tiaParams_.restartTimeStepScale;

  // if initial time steps are baloney, set then to a default value.
  if (startingTimeStep <= 0.0) startingTimeStep = 1.0e-10;
  if (currentTimeStep <= 0.0)  currentTimeStep  = 1.0e-10;

  initializeBreakPoints();

  // 12/15/06 tscoffe:  The MPDE Manger calls resetAll and then setTIAParams.
  // resetAll correctly sets up the currentPauseBP iterator, but setTIAParams
  // did not.  This little block of code should fix that.
  if( anaManager_.getBlockAnalysisFlag() == true )
  {
    set<N_UTL_BreakPoint>::iterator lastBP  = breakPoints_.end();
    lastBP--;
    updatePauseTime(*lastBP);
  }

  h0_max_factor_ = tiaParams_.restartTimeStepScale;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::getEstOverTol
// Purpose       : This function lets the controlling class (a transient
//                 analysis) get the estimated error over tol from the
//                 last step
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 01/23/09
//-----------------------------------------------------------------------------
double N_TIA_StepErrorControl::getEstOverTol() const
{
  return estOverTol_;
}



//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::getTranOPFlag
// Purpose       : Get flag from analysis manager for DC op part of a 
//                 transient run.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 01/23/09
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::getTranOPFlag() const
{
  return anaManager_.getTranOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::setTimeStep
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 08/09/09
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::setTimeStep(double newTimeStep)
{
  newTimeStep = Xycemax(newTimeStep, minTimeStep);
  newTimeStep = Xycemin(newTimeStep, maxTimeStep);

  double nextTimePt = currentTime + newTimeStep;

  if (nextTimePt > stopTime)
  {
    nextTimePt  = stopTime;
    newTimeStep = stopTime - currentTime;
    tiaParams_.TimeStepLimitedbyBP = true;
  }

  nextTime = nextTimePt;

  currentTimeStepRatio = newTimeStep/lastTimeStep;
  currentTimeStepSum   = newTimeStep + lastTimeStep;

  currentTimeStep = newTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::resetAll
//
// Purpose       : This function resets everything so that a transient loop
//                 can be started from the beginning.
//
// Special Notes : This function was needed for the .STEP capability.
// 
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 11/04/03
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::resetAll ()
{
  // reset a bunch of variables.  
  tiaParams_.pauseSetAtZero = false;
  tiaParams_.pauseTime = 0.0;

  startingTimeStep = tiaParams_.userSpecified_startingTimeStep;
  currentTimeStep  = tiaParams_.userSpecified_startingTimeStep;
  initialTime      = tiaParams_.initialTime;
  finalTime        = tiaParams_.finalTime;

  lastTimeStep     = tiaParams_.userSpecified_startingTimeStep;
  lastAttemptedTimeStep = tiaParams_.userSpecified_startingTimeStep;

  currentTime = initialTime;
  lastTime    = initialTime;
  nextTime    = initialTime;

  currentTimeStepRatio = 1.0;
  currentTimeStepSum   = 2.0*currentTimeStep;

  lastTimeStep      = currentTimeStep;
  lastTimeStepRatio = currentTimeStepRatio;
  lastTimeStepSum   = currentTimeStepSum;

  newtonConvergenceStatus = -1;
  numberSuccessiveFailures = 0;
  stepAttemptStatus        = true;

  minTimeStep = 0.0;
  estOverTol_ = 0.0;

  if (tiaParams_.maxTimeStepGiven)
  {
    maxTimeStepUser  = tiaParams_.maxTimeStep;
    maxTimeStep      = tiaParams_.maxTimeStep;
  }
  else
  {
    maxTimeStep  = 0.1* (tiaParams_.finalTime - tiaParams_.initialTime);
  }

  restartTimeStepScale_ = tiaParams_.restartTimeStepScale;

  // if initial time steps are baloney, set then to a default value.
  if (startingTimeStep <= 0.0) startingTimeStep = 1.0e-10;
  if (currentTimeStep <= 0.0)  currentTimeStep  = 1.0e-10;

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams_.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_StepErrorControl::resetAll");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  before initializeBreakPoints()");
    printBreakPoints(cout);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  currentPauseBP = ", currentPauseBP->value());

  }
#endif // Xyce_DEBUG_TIME

  initializeBreakPoints();

#ifdef Xyce_DEBUG_TIME
  if (tiaParams_.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  after initializeBreakPoints()");
    printBreakPoints(cout);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  currentPauseBP = ", currentPauseBP->value());
  }
#endif // Xyce_DEBUG_TIME
 
   // ERK: Note:  always do this, not just for block solves?
  //if( anaManager_.getBlockAnalysisFlag() == true ) // if MPDE or HB
  //{
    set<N_UTL_BreakPoint>::iterator lastBP  = breakPoints_.end();
    lastBP--;
    updatePauseTime(*lastBP);
  //}

  // need to set a pause breakpoint at the final time.
  setBreakPoint(N_UTL_BreakPoint(tiaParams_.finalTime ,PAUSE_BREAKPOINT));

#ifdef Xyce_DEBUG_TIME
  if (tiaParams_.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  after updatePauseTime() & setBreakPoint");
    printBreakPoints(cout);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  currentPauseBP = ", currentPauseBP->value());
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::initializeStepSizeVariables
// Purpose       : Remainder of initialization for the StepSize Process.
//                 This should be called every time an integration start-up is
//                 needed -- such as at the beginning of integration and after
//                 a point of discontinuity is encountered, but not upon
//                 continuation from a saved solution state (i.e. a restart
//                 file).
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::initializeStepSizeVariables()
{
  double time_to_stop = stopTime - currentTime;

#ifdef Xyce_VERBOSE_TIME
  const string crMsg         = "\n";
  const string startMsg      = ("* Initial Time Value:\t");
  const string stopMsg       = (" *  Ending Time Value:\t");
  const string secondsMsg    = (" secs");
  const string timeToStopMsg = (" * Time to Stop Value:\t");

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, crMsg + startMsg,
                         currentTime, secondsMsg);

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, stopMsg,
                         stopTime, secondsMsg);

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, timeToStopMsg,
                         time_to_stop, secondsMsg);

#endif

  // At some point, we might want to consider computing an on-scale starting
  // step-size based on local behavior of problem solution

  if (tiaParams_.constantStepSize)
  {
    currentTimeStep = 0.1 * time_to_stop;
  }
  else
  {
    // The original of this was just:
    currentTimeStep = restartTimeStepScale_ * time_to_stop;

    // This can cause problems in the cases of very steep pulsed sources,
    // e.g. those with nanosecond scale rise and fall times.  Then we get
    // order picosecond timesteps after breakpoints, and this is often too
    // bloody small.
    // I (TVR) tried this but it is not smart enough to cover all cases.
    //    currentTimeStep = Xycemax(0.005*time_to_stop,5e-11);
  }

  if (currentTime == initialTime || tiaParams_.constantStepSize)
  {
    currentTimeStep = Xycemin(startingTimeStep, currentTimeStep);
  }
  else
  {
    currentTimeStep = Xycemin(lastTimeStep, currentTimeStep);
  }

#ifdef Xyce_VERBOSE_TIME
  const string stepMsg      = (" * Initial Step Size: \t");

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, stepMsg,
                         currentTimeStep, secondsMsg + crMsg);
#endif // Xyce_VERBOSE_TIME

  currentTimeStepRatio = 1.0;
  currentTimeStepSum   = 2.0*currentTimeStep;

  lastTimeStep      = currentTimeStep;
  lastTimeStepRatio = currentTimeStepRatio;
  lastTimeStepSum   = currentTimeStepSum;

  numberSuccessiveFailures = 0;
  stepAttemptStatus        = true;

  nextTime = currentTime + currentTimeStep;
  
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::updateStopTime
// Purpose       : The "stop time" is either the next discontinuity point,
//                 a pause point, or the final time, whichever comes first.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::updateStopTime()
{

//#ifdef Xyce_MPDE
  // 10/04/05 tscoffe:  In MPDE mode, we need to set the final time, but
  // StepErrorControl only reads finalTime from TIAParams at construction.  So
  // here, we're forcing it to reload it if its in MPDE mode.
//  finalTime = tiaParams_.finalTime;
//#endif // Xyce_MPDE

  double oldStopTime = stopTime;
  double diffStopTime(0.0);

  if (tiaParams_.bpEnable)
  {
    set<N_UTL_BreakPoint>::iterator itBP;
    set<N_UTL_BreakPoint>::iterator firstBP = breakPoints_.begin();
    set<N_UTL_BreakPoint>::iterator lastBP  = breakPoints_.end();

    // Find the first breakpoint equal to or larger than the
    // current time.

    itBP = upper_bound(firstBP,lastBP,N_UTL_BreakPoint(currentTime));

    stopTime =  Xycemin(finalTime, itBP->value());

    // The breakpoint could be a pause breakpoint, in which case we might
    // need to update the pauseTime:
    if (itBP->bptype() == PAUSE_BREAKPOINT)
    {
      updatePauseTime(*itBP);
    }
    stopTime =  Xycemin(currentPauseBP->value(), stopTime);

    // if this is a breakpoint step, make sure the new stop
    // time is for the next breakpoint, not the current one.
    // This check is neccessary because of roundoff error.

    diffStopTime = fabs(stopTime-oldStopTime);
    if(diffStopTime < anaManager_.getBreakpointTol() &&
       anaManager_.getBeginningIntegrationFlag() &&
       stopTime != tiaParams_.pauseTime &&
       stopTime != finalTime )
    {
      ++itBP;
      stopTime = itBP->value();
    }

#ifdef  Xyce_PARALLEL_MPI
    double sT = stopTime;
    double mST;

    N_PDS_Comm * comm = anaManager_.pdsMgrPtr->getPDSComm();
  // NOTE:  This needs to be a minAll !
    comm->minAll ( &sT, &mST, 1);
    stopTime = mST;
#endif // Xyce_PARALLEL_MPI 

  }
  else
  {
    stopTime = Xycemin(tiaParams_.pauseTime, finalTime);
  }

  if ( anaManager_.getBeginningIntegrationFlag()) 
  {
    double time_to_stop = stopTime - currentTime;
    if (tiaParams_.minTimeStepsBPGiven && (tiaParams_.minTimeStepsBP > 0) ) 
    {
      maxTimeStepBP = time_to_stop/tiaParams_.minTimeStepsBP;
    }
  }


#ifdef Xyce_DEBUG_TIME
  if (tiaParams_.debugLevel >0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  stopTime    = ", stopTime);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  pauseTime    = ", tiaParams_.pauseTime);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  currentTime = ", currentTime);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  oldStopTime = ", oldStopTime);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  finalTime   = ", finalTime);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  maxTimeStepBP = ", maxTimeStepBP);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  beginningIntegration = ", anaManager_.getBeginningIntegrationFlag());
    const string dashedline = 
      "-------------------------------------------------"
                       "----------------------------\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME 

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::findNextStopTime
// Purpose       : Determine the next stop time, that comes immediately after
//                 the current one.
//
// Special Notes : Used by mixed-signal rollback.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/9/09
//-----------------------------------------------------------------------------
double N_TIA_StepErrorControl::findNextStopTime()
{
  double nextStop=stopTime;

  if (tiaParams_.bpEnable)
  {
    set<N_UTL_BreakPoint>::iterator itBP;
    set<N_UTL_BreakPoint>::iterator firstBP = breakPoints_.begin();
    set<N_UTL_BreakPoint>::iterator lastBP  = breakPoints_.end();

    // Find the first breakpoint equal to or larger than the
    // current time.

    itBP = upper_bound(firstBP,lastBP,N_UTL_BreakPoint(currentTime));
    itBP++; // go to the next one.

    nextTime =  Xycemin(finalTime, itBP->value());
    nextTime =  Xycemin(currentPauseBP->value(), nextTime);

#ifdef  Xyce_PARALLEL_MPI
    double nT = nextTime;
    double mNT;

    N_PDS_Comm * comm = anaManager_.pdsMgrPtr->getPDSComm();
  // NOTE:  This needs to be a minAll !
    comm->minAll ( &nT, &mNT, 1);
    nextTime = mNT;
#endif // Xyce_PARALLEL_MPI 

  }
  else
  {
    nextTime = Xycemin(tiaParams_.pauseTime, finalTime);
  }

  return nextStop;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::evaluateStepError
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 1/28/07
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::evaluateStepError ()
{
  bool step_attempt_status( newtonConvergenceStatus >= 0);
  bool sAStatus(false);
  bool errorOptionStatus(true);
  bool testTimeIntegrationError(false);

  // If we are running with constant step size, or are on the first pass
  // through the transient loop, only base success on the Newton loop.

//  if (anaManager_.getIntegrationMethod() == 7 && tiaParams_.newLte == true)
  if (tiaParams_.newBPStepping == true)
  {
    if (currentTime == tiaParams_.initialTime)
    {
      testTimeIntegrationError = (anaManager_.getStepNumber() >= 1 && !(anaManager_.getBeginningIntegrationFlag()));
    }
    else
    {
      testTimeIntegrationError = (anaManager_.getStepNumber() >= 1); 
    }
  }
  else
  {
    testTimeIntegrationError = (anaManager_.getStepNumber() >= 1 && !(anaManager_.getBeginningIntegrationFlag())); 
  }

  if (tiaParams_.testFirstStep)
  {
    testTimeIntegrationError = true;
  }

  // if the step status is already false, don't do any more work.
  if (!step_attempt_status)
  {
    testTimeIntegrationError = false;
  }

  
  if (testTimeIntegrationError)
  {
#ifdef Xyce_EXTDEV
    // Needed for 2-level Solves:
    loader_.getInnerLoopErrorSums (anaManager_.dsPtr_->innerErrorInfoVec);
#endif

    estOverTol_ = wimPtr_->computeErrorEstimate();

    if (estOverTol_ <= tiaParams_.errTolAcceptance)
    {
      sAStatus = true;
    }
    else
    {
      sAStatus = false;
    }
  
    if (tiaParams_.timestepsReversal == true)
    { 
      if (nIterations <= tiaParams_.NLmax)
        errorOptionStatus = true;
      else
        errorOptionStatus = false;
    }
#ifdef Xyce_VERBOSE_TIME
    if (tiaParams_.errorAnalysisOption == 1) 
    {
      cout << "ERROROPTION=1:  DOREJECTSTEP = ";
      if (tiaParams_.timestepsReversal == true) 
      {
        cout << "1" << endl;
      } 
      else 
      {
        cout << "0" << endl;
      }
    }
#endif // Xyce_VERBOSE_TIME
    
    if ( tiaParams_.userSpecMinTimeStepGiven && (currentTimeStep < tiaParams_.userSpecMinTimeStep) )
    {
      // This step has dropped under the user specified min time step, so only
      // test if the sovler converged to accept the step.
      // don't do the step_attempt_status && sAStatus;
#ifdef Xyce_DEBUG_TIME
      if (tiaParams_.debugLevel > 2)
      {
        std::cout << "Trying to skip time integrator error checks: " << currentTimeStep 
          << " newton status " << step_attempt_status << std::endl;
      }
#endif
    }
//    else if (!(tiaParams_.constantStepSize) && ((tiaParams_.errorAnalysisOption == 0)) )
//    {
//      step_attempt_status = step_attempt_status && sAStatus;
//    }
      else if (!(tiaParams_.constantStepSize))
      {
        if (tiaParams_.errorAnalysisOption == 1)
          step_attempt_status = step_attempt_status && errorOptionStatus;  
        else
          step_attempt_status = step_attempt_status && sAStatus;
      }
    }


#ifdef Xyce_VERBOSE_TIME
#ifdef Xyce_DEBUG_TIME
  if (tiaParams_.debugLevel > 0)
  {
    integrationStepReport_ (step_attempt_status, sAStatus, testTimeIntegrationError);
  }
  else
#endif
  {
    terseIntegrationStepReport_ 
      (step_attempt_status, sAStatus, testTimeIntegrationError);
  }
#endif

  // Now that the status has been completely determined, 
  // set the class variable for step attempt 
  stepAttemptStatus = step_attempt_status;

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::terseIntegrationStepReport_
// Purpose       : This gives a one-line description of the step accept/reject.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/27/07
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::terseIntegrationStepReport_
  (bool step_attempt_status, bool sAStatus, bool testedError)
{
    ostringstream output1;

#ifdef Xyce_DEBUG_TIME
  string netListFile = commandLine_.getArgumentValue("netlist");
  output1 << netListFile << "  STEP STATUS:";
#else
  output1 << "STEP STATUS:";
#endif

  if (step_attempt_status)
  {
    output1 << " success ";
  }
  else
  {
    output1 << " fail    ";
  }

  output1 << " Newton: " << newtonConvergenceStatus;

  if (testedError && !(tiaParams_.constantStepSize))
  {
    output1 << "   estOverTol: " << estOverTol_;
  }
  else
  {
    output1 << "   estOverTol: " << estOverTol_  << " (not used for this step)";
  }

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, output1.str());
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::integrationStepReport_
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::integrationStepReport_
  (bool step_attempt_status, bool sAStatus, bool testedError)
{
      ostringstream output1;

  string dashedline = "-------------------------------------------------"
                       "----------------------------\n";

  if (tiaParams_.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "\n estOverTol      = ", estOverTol_);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                              "  error tolerance = ",
                              tiaParams_.errTolAcceptance);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");


    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "\nSTEP ATTEMPT STATUS:");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "NOTE:");

    if (!(tiaParams_.constantStepSize) &&
          anaManager_.getStepNumber() >= 1             &&
        !(anaManager_.getBeginningIntegrationFlag()))
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  We are running in variable stepsize mode,");

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  and we have NOT just passed a breakpoint.  As such,");

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  for an integration step to succeed, the ");

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  nonlinear solver must succeed, AND the predictor");

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  and corrector need to be close within a tolerance.");

      if (tiaParams_.errorAnalysisOption == 1)
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                                "ADDENDUM:  This is with erroption=1, "
                                "so predictor-corrector is ignored for "
                                "step error control.");
      }
    }
    else
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  We are either running constant stepsize,");

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  or we just passed a breakpoint.  As such,");

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  the only criteria we use in accepting/rejecting");

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  an integration step is the nonlinear solver");

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  success/failure.");
    }

    if (step_attempt_status)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "\n  This has been a successful step:");
    }
    else
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "\n  This has NOT been a successful step:");
    }

    if ( newtonConvergenceStatus > 0)
    {
      output1 << "    - Newton solver succeded with return code " 
              << newtonConvergenceStatus;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              output1.str());
    }
    else
    {
      output1 << "    - Newton solver failed with return code " 
              << newtonConvergenceStatus;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              output1.str());
    }

    if (testedError)
    {
      if (!(tiaParams_.constantStepSize))
      {
        if (sAStatus)
        {
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                  "   - predictor vs. corrector analysis succeeded.");
        }
        else
        {
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                  "   - predictor vs. corrector analysis failed.");
        }

        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                "     (compare estOverTol with error tolerance, above.)");
      }
      else
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                "If we had been using it, ");

        if (sAStatus)
        {
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "   - predictor vs. corrector analysis would have succeeded.");
        }
        else
        {
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "   - predictor vs. corrector analysis would have failed.");
        }

        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "     (compare estOverTol with error tolerance, above.)");
      }
    }
    else
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  predictor vs. corrector was not tested");
    }

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

  } // end of debugLevel if statement
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::initializeBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 06/11/01
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::initializeBreakPoints()
{
  bool bsuccess = true;

  breakPoints_.clear ();
  // if currentPauseBP was set, it is now invalid so set it to a safe but
  // invalid value until it is calulcated
  currentPauseBP = breakPoints_.end();

  // first breakpoint is the start time, last one is the final time.
  setBreakPoint(initialTime);

  // if tStart is nonzero, then make it a breakpoint.
  if (tiaParams_.tStart > initialTime &&
      tiaParams_.tStart < finalTime)
  {
    setBreakPoint(tiaParams_.tStart);
  }

  // The final time needs to be the very last breakpoint.
  //setBreakPoint(finalTime);
  setBreakPoint(N_UTL_BreakPoint(tiaParams_.finalTime ,PAUSE_BREAKPOINT));

  // Now add in the user break points. NOT DONE YET.
  //vector<double> ;
  //tiaParams_.userSpecBreakPoints;


  // other breakpoints will later be added/removed dynamically over the
  // course of the transient simulation in calls to the loader.

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::updateBreakPoints
// Purpose       : Requests dynamic breakpoint information from the
//                 loader.  Adds, subtracts from the breakpoints array.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 06/11/01
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::updateBreakPoints ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
"----------------------------------------------------------------------------";
  if (tiaParams_.debugLevel >0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  N_TIA_StepErrorControl::updateBreakPoints.  time = ", currentTime);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
  }
#endif

  vector<N_UTL_BreakPoint> tmpBP;
  tmpBP.clear ();

  loader_.getBreakPoints(tmpBP);

  // debug outputs:
  vector<N_UTL_BreakPoint>::iterator iter;
  vector<N_UTL_BreakPoint>::iterator first = tmpBP.begin();
  vector<N_UTL_BreakPoint>::iterator last  = tmpBP.end();

  // add any new breakpoints to the vector of breakpoints:
  set<N_UTL_BreakPoint>::iterator itBP;
  set<N_UTL_BreakPoint>::iterator itBP_2;
  set<N_UTL_BreakPoint>::iterator firstBP = breakPoints_.begin();
  set<N_UTL_BreakPoint>::iterator lastBP  = breakPoints_.end();

  // Add new breakpoints to the set:
  for (iter=first; iter!=last; ++iter)
  {
    if (iter->value() < finalTime && iter->value() > lastTime) 
    {
      setBreakPoint(*iter);
    }
  }

#ifdef Xyce_DEBUG_TIME
  char tmp[128];
  int i;
  if (tiaParams_.debugLevel >0)
  {
    firstBP = breakPoints_.begin();

    string netListFile = commandLine_.getArgumentValue("netlist");
    string msg = netListFile + "  breakPoints_ vector container, before any removals:";

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,msg);

    for (i=0, itBP=firstBP;itBP!=lastBP;++i,++itBP)
    {
      if (i==0)
      {
        sprintf(tmp,"%4d %16.8e",i,itBP->value());
      }
      else
      {
        sprintf(tmp,"%4d %16.8e diff=%16.8e", i, itBP->value(),(itBP->value()-itBP_2->value()));
      }

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,tmp);
      itBP_2 = itBP;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,"");
  }
#endif

  // Remove breakpoints which are now obsolete (old):
  LessThan<N_UTL_BreakPoint,double> LessFunct;
  itBP_2 = lower_bound(firstBP,lastBP,lastTime,LessFunct);
  breakPoints_.erase(firstBP,itBP_2);

  // Remove breakpoints which are too close together.
  // ERK.  This is kind of ugly and can undoubtably be done
  // in a much more elegant way using STL.  Note that
  // the value of "bpTol" is pulled out of my ear.

  // TVR: only need to do this if the BP tolerance has changed from what it
  // was set to in the breakpoint class --- as things were inserted, near
  // values were rejected if within tolerance.  Only need to do this if the 
  // tolerance has gotten looser

  // TVR: Cannot use the STL "unique" function to accomplish this, because
  // STL sets are designed to prevent changing the set elements after they're
  // inserted.  Have to do this by removing duplicate entries explicitly.

  double bpTol = 2.0 * minTimeStep;
  if (bpTol != anaManager_.getBreakpointTol())
  {
    // set the class's tolerance to this new one
    anaManager_.setBreakpointTol(bpTol);
#ifdef Xyce_DEBUG_TIME
    if (tiaParams_.debugLevel >0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "  bpTol = ",bpTol);
      cout << " Must now eliminate new duplicates " << endl;
    }
#endif

    bool doneRemove = false;
    while (!doneRemove)
    {
      doneRemove = true;
      firstBP = breakPoints_.begin();
      lastBP  = breakPoints_.end();
      int icount;
      for (icount = 0, itBP=firstBP, itBP_2=firstBP;
           itBP!=lastBP;
           ++icount, ++itBP)
      {
        double diff = (itBP->value() - itBP_2->value());
        
        if (icount != 0)
        {
          if (fabs(diff) < bpTol)
          {
            // If both are simple, just toss the later one
            if (itBP->bptype() == SIMPLE_BREAKPOINT && 
                itBP_2->bptype() == SIMPLE_BREAKPOINT)
            {
              if (diff > 0.0)
                breakPoints_.erase(itBP);
              else
                breakPoints_.erase(itBP_2);
            }
            else
            {
              // one of these breakpoints is not simple!  Determine the 
              // overriding type, then set a breakpoint at the earliest time
              // with the overriding type:
              BreakpointType overridingType = itBP->bptype();
              double minTime=Xycemin(itBP->value(),itBP_2->value());
              
              // The following line will need to be changed if any other types
              // besides SIMPLE_BREAKPOINT and PAUSE_BREAKPOINT are ever
              // introduced and any more complex precedence is defined.
              if (itBP_2->bptype() != SIMPLE_BREAKPOINT)
              {
                overridingType = itBP_2->bptype();
              }
              breakPoints_.erase(itBP);
              breakPoints_.erase(itBP_2);
              N_UTL_BreakPoint tmpBP(minTime,overridingType);
              breakPoints_.insert(tmpBP);
              // and just to be very careful, make sure to update the pause
              // time, lest the currentPauseBP iterator be confused.
#ifdef Xyce_DEBUG_TIME
              cout << " Purging breakpoints, overriding with breakpoint of type " << tmpBP.bptype();
#endif
              updatePauseTime(tmpBP);
            }

            doneRemove = false;
            break;
          }
        }
        itBP_2 = itBP;
      }
    }
  }

#ifdef Xyce_DEBUG_TIME
  if (tiaParams_.debugLevel >0)
  {
    firstBP = breakPoints_.begin();

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
       "  breakPoints_ vector container, after:");

    for (i=0, itBP=firstBP;itBP!=lastBP;++i,++itBP)
    {
      if (i==0)
      {
        sprintf(tmp,"%4d %16.8e  type=%d",i,itBP->value(),itBP->bptype());
      }
      else
      {
        sprintf(tmp,"%4d %16.8e type=%d diff=%16.8e", i, itBP->value(),
                itBP->bptype(),(itBP->value()-itBP_2->value()));
      }

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,tmp);
      itBP_2 = itBP;
    }

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,"");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,dashedline);
  }
#endif

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::updateMaxTimeStep
// Purpose       : Requests dynamic time step information from the loader.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::updateMaxTimeStep(double suggestedMaxTimeStep)
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
"----------------------------------------------------------------------------";
  if (tiaParams_.debugLevel >0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  N_TIA_StepErrorControl::updateMaxTimeStep");
  }
#endif

  double maxDevStep = 1.0e+99;
  if (tiaParams_.useDeviceTimeStepMax)
  {
    maxDevStep = loader_.getMaxTimeStepSize();
  }

  if (tiaParams_.maxTimeStepGiven || tiaParams_.delmaxGiven)
  {
    maxTimeStep = Xycemin(tiaParams_.maxTimeStep, tiaParams_.delmax);
  }
  else
  {
    maxTimeStep  = 0.1*(tiaParams_.finalTime-tiaParams_.initialTime);
  }

  // if the default arg is not zero, then a suggested max time step was
  // passed in.  Test if it is feasible to use that at this time
  if( suggestedMaxTimeStep > 0.0 )
  {
    maxTimeStep = Xycemin( maxTimeStep, suggestedMaxTimeStep );
  }

  if ((maxTimeStepBP > 0.0) && (maxTimeStep > maxTimeStepBP)) 
  {
    maxTimeStep = maxTimeStepBP;
  }

  if (maxDevStep > 0.0)
  {
    maxTimeStep = Xycemin(maxTimeStep, maxDevStep);
  }

  if (tiaParams_.maxTimeStepGiven)
  {
    maxTimeStep = Xycemin(maxTimeStep, maxTimeStepUser);
  }

#ifdef  Xyce_PARALLEL_MPI
  double mTS = maxTimeStep;
  double mMTS;

  N_PDS_Comm * comm = anaManager_.pdsMgrPtr->getPDSComm();

  comm->minAll ( &mTS, &mMTS, 1);

  maxTimeStep = mMTS;
#endif

#ifdef Xyce_DEBUG_TIME
  if (tiaParams_.debugLevel >0)
  {
    if (!(tiaParams_.maxTimeStepGiven))
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  User did not specify a maximum time step.");
    }
    else
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  User specified a maximum time step. = ", maxTimeStepUser);
    }

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  maxDevStep  = ", maxDevStep);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
      "  maxTimeStep = ", maxTimeStep);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,dashedline);
  }
#endif

  if(maxTimeStep<=0.0)
  {
    const string msg = "Maximum Time step is invalid!\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::updateMinTimeStep
// Purpose       : Sets the minimum time step based on machine precision.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 08/02/01
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::updateMinTimeStep()
{
  bool bsuccess = true;

  minTimeStep = currentTime * minStepPrecisionFac_ *
    N_UTL_MachineDependentParams::MachinePrecision();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::setBreakPoint
// Purpose       : public method to set individual breakpoint
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::setBreakPoint(N_UTL_BreakPoint theBP)
{
  if (theBP.bptype() == SIMPLE_BREAKPOINT)
  {
    breakPoints_.insert(theBP);
  }
  else
  {
#ifdef Xyce_DEBUG_TIME
  ostringstream ost;
    ost << "In setBreakPoint, got non-simple breakpoint of type " 
         << theBP.bptype() << " at time " << theBP.value()  << endl;

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, ost.str());
#endif
    // We're a pause breakpoint, and must override any simple breakpoint
    // at the same time.
    // erase any breakpoint that is the "same" time as this one
    breakPoints_.erase(theBP);
    // and put the new one in instead
    breakPoints_.insert(theBP);
    // force code to recalculate pause time NOW.
    updatePauseTime(theBP);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::setBreakPoint
// Purpose       : public method to set individual SIMPLE breakpoint
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::setBreakPoint(double theBP)
{
    breakPoints_.insert(theBP);
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::updatePauseTime
// Purpose       : private method to recalculate pause time
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::updatePauseTime(N_UTL_BreakPoint BP)
{
  // gotta handle case where pauseTime is still at its initial value of
  // 0.0, or we already passed the last pause time!
  // But we mustn't reset it if it's equal to the current time, because
  // that means we need to stop NOW and would overwrite that.
  //   
  // If a pause break point is specifically set at 0, then 
  // we shouldn't ignore that here.  So set the pauseSetAtZero flag
  // here if needed.
  //
  // ERK: type>0 means breakpoint is "not simple", ie PAUSE breakpoint.
  if ((BP.bptype() > 0) && (BP.value() == 0.0)) 
  {
    tiaParams_.pauseSetAtZero = true;
  }
  
  if (tiaParams_.pauseTime < currentTime ||
     ((tiaParams_.pauseTime == tiaParams_.initialTime) && !tiaParams_.pauseSetAtZero))
  {
    tiaParams_.pauseTime = BP.value();
  }
  else
  {
    tiaParams_.pauseTime = Xycemin(tiaParams_.pauseTime, BP.value());
  }

  // If we used this BP for the pause time, save the iterator to this bp in 
  // the list so we can use it later.
  if (tiaParams_.pauseTime == BP.value())
  {
    currentPauseBP = breakPoints_.find(BP);
#ifdef Xyce_DEBUG_TIME
    ostringstream ost;
    string netListFile = commandLine_.getArgumentValue("netlist");
    ost << "\n" << netListFile;
    ost << "  UPDATING PAUSE TIME TO " << currentPauseBP->value() 
      << " encountered breakpoint " << BP.value() << " current time  is " 
      << currentTime << endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, ost.str());
#endif
  }

}


//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::simulationPaused
// Purpose       : public method to clear out breakpoint list and reset
//                 pause time when pause time is reached
// Special Notes : Need this method because the prior values get in the way 
//                 when we resume.
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/30/04
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::simulationPaused()
{
  breakPoints_.erase(currentTime);  // clear current breakpoint
  tiaParams_.pauseTime = tiaParams_.initialTime;          // unset this
  currentPauseBP = breakPoints_.end(); // make this invalid
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::printBreakPoints
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::printBreakPoints (ostream & os) const
{
  set<N_UTL_BreakPoint>::const_iterator itBP;
  set<N_UTL_BreakPoint>::const_iterator itBP2;
  set<N_UTL_BreakPoint>::const_iterator firstBP = breakPoints_.begin();
  set<N_UTL_BreakPoint>::const_iterator lastBP  = breakPoints_.end();

  char tmp[128];

  int i;
  for (i=0, itBP=firstBP;itBP!=lastBP;++i,++itBP)
  {
    if (i==0)
      sprintf(tmp,"%4d %16.8e  type=%d",i,itBP->value(),itBP->bptype());
    else
      sprintf(tmp,"%4d %16.8e type=%d diff=%16.8e", i, itBP->value(),
	      itBP->bptype(),(itBP->value()-itBP2->value()));

    os << string(tmp);
    itBP2 = itBP;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::restartDataSize
// Purpose       :
// Special Notes : This gives the *total* size:  both the base 
//                 N_TIA_StepErrorControl and the derived 
//                 N_TIA_StepErrorControlDAE, summed together.
//
//                 Don't sum them 2x!
//
//                 ERK:  6/20/2010:  
//                 There are no longer 2 distinct classes.  The DAE version is
//                 now part of the original, as one single class.  The two 
//                 blocks of code in this function correspond to the original
//                 class (first block), and the DAE class (second block).
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 07/27/06
//-----------------------------------------------------------------------------
int N_TIA_StepErrorControl::restartDataSize( bool pack )
{
  int totalSize = 0;

  // original set of vars:
  int numdoubles = 21;
  int numints = 9;
  int count = sizeof(double) * (numdoubles);
  count += sizeof(int) * numints;

  // Must include the bp type now, not just the value
  count += sizeof(N_UTL_BreakPoint)*breakPoints_.size();

  //overestimate buffer size for unpacked data
  // assumes there are fewer than 100 possible breakpoints types (i.e.
  // because the type can be represented by two ascii characters)
  if( !pack ) 
  {                  // 34
    count = 24*((numdoubles+numints) + (breakPoints_.size()))
                + 2*(breakPoints_.size());
  }

  int baseClassSize = count;

  // another set of vars: (newDAE)
  numdoubles = 57;
  numints = 10;

  totalSize = baseClassSize;
  totalSize += sizeof(double) * numdoubles;
  totalSize += sizeof(int) * numints;

  //overestimate buffer size for unpacked data
  if ( !pack )
  {
    totalSize = baseClassSize + 24*(numdoubles+numints);
  }

  return totalSize;
}
//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::dumpRestartData
// Purpose       : 
// Special Notes : ERK:  6/20/2010:  
//                 There are no longer 2 distinct classes.  The DAE version is
//                 now part of the original, as one single class.  The two 
//                 blocks of code in this function correspond to the original
//                 class (first block), and the DAE class (second block).
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 7/28/06
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::dumpRestartData
  (char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack )
{

  // Set this variable up for later.  Note that pos means different things
  // for packed vs. unpacked.  For unpacked, it is the current index into
  // the buf array.
  int newPos = pos + N_TIA_StepErrorControl::restartDataSize(false);

  // if unpacked, initialize the buf array before calling the base class
  // function.  At this point pos is zero, probably.  The derived DAE
  // dataSize will be the size of the entire buf array, including the
  // base class size.
  if ( !pack )
  {
    for( int i = pos; i < (newPos); ++i) buf[i] = ' ';
  }

#ifdef Xyce_DEBUG_RESTART
  string netListFile = commandLine_.getArgumentValue("netlist");
  cout << "TIA Restart Data DUMP!  " << netListFile << "\n";
  cout << "----------------------------------------\n";
  cout << "startingTimeStep: " << startingTimeStep << endl;
  cout << "currentTimeStep: " << currentTimeStep << endl;
  cout << "lastAttemptedTimeStep: " << lastAttemptedTimeStep << endl;
  cout << "lastTimeStep: " << lastTimeStep << endl;
  cout << "minTimeStep: " << minTimeStep << endl;
  cout << "maxTimeStep: " << maxTimeStep << endl;
  cout << "maxTimeStepUser: " << maxTimeStepUser << endl;
  cout << "lastTime: " << lastTime << endl;
  cout << "currentTime: " << currentTime << endl;
  cout << "nextTime: " << nextTime << endl;
  cout << "initialTime: " << initialTime << endl;
  cout << "estOverTol_: " << estOverTol_ << endl;
  cout << "breakpts: ";

  for (set<N_UTL_BreakPoint>::iterator iterSD = breakPoints_.begin();
       iterSD != breakPoints_.end(); ++iterSD )
    cout << iterSD->value() << " ";
  cout << endl;
  cout << "integMethod: " << anaManager_.getIntegrationMethod() << endl;
  cout << "stepNumber: " << anaManager_.getStepNumber() << endl;
  cout << "transStepNumber: " << anaManager_.getTranStepNumber() << endl;
  cout << "breakPointRestartNumber: " << anaManager_.breakPointRestartStep << endl;
  cout << "----------------------------------------\n\n";
#endif
  if( pack )
  {
    comm->pack( &startingTimeStep, 1, buf, bsize, pos );
    comm->pack( &currentTimeStep, 1, buf, bsize, pos );
    comm->pack( &lastAttemptedTimeStep, 1, buf, bsize, pos );
    comm->pack( &lastTimeStep, 1, buf, bsize, pos );
    comm->pack( &minTimeStep, 1, buf, bsize, pos );
    comm->pack( &maxTimeStep, 1, buf, bsize, pos );
    comm->pack( &maxTimeStepUser, 1, buf, bsize, pos );
    comm->pack( &lastTime, 1, buf, bsize, pos );
    comm->pack( &currentTime, 1, buf, bsize, pos );
    comm->pack( &nextTime, 1, buf, bsize, pos );
    comm->pack( &initialTime, 1, buf, bsize, pos );
    comm->pack( &currentTimeStepRatio, 1, buf, bsize, pos );
    comm->pack( &currentTimeStepSum, 1, buf, bsize, pos );
    comm->pack( &lastTimeStepRatio, 1, buf, bsize, pos );
    comm->pack( &lastTimeStepSum, 1, buf, bsize, pos );
    comm->pack( &newtonConvergenceStatus, 1, buf, bsize, pos );
    comm->pack( &numberSuccessiveFailures, 1, buf, bsize, pos );
    int flag = stepAttemptStatus;
    comm->pack( &flag, 1, buf, bsize, pos );
    comm->pack( &minStepPrecisionFac_, 1, buf, bsize, pos );
    comm->pack( &newtonStepReduction_, 1, buf, bsize, pos );
    comm->pack( &tolAimFac_, 1, buf, bsize, pos );
    comm->pack( &estOverTol_, 1, buf, bsize, pos );
    // Subtract one, because we won't write out the pause breakpoint at the
    // final time
    int size = breakPoints_.size() -1 ;
    set<N_UTL_BreakPoint>::iterator bpStart = breakPoints_.begin();
    set<N_UTL_BreakPoint>::iterator bpEnd = breakPoints_.end();
    comm->pack( &size, 1, buf, bsize, pos );

    double val;
    int bptype;
    for( set<N_UTL_BreakPoint>::iterator iterSD = bpStart;
       iterSD != bpEnd; ++iterSD)
    {
      val=iterSD->value();
      bptype=iterSD->bptype();
      if (!(bptype == PAUSE_BREAKPOINT && val == finalTime))
      {
        comm->pack( &(val), 1, buf, bsize, pos );
        comm->pack( &(bptype), 1, buf, bsize, pos );
      }
    }
    int im = anaManager_.getIntegrationMethod();
    comm->pack( &im, 1, buf, bsize, pos );
    int sN = anaManager_.getStepNumber();
    comm->pack( &sN, 1, buf, bsize, pos );
    int tSN = anaManager_.getTranStepNumber();
    comm->pack( &tSN, 1, buf, bsize, pos );
    int bPRS = anaManager_.breakPointRestartStep;
    comm->pack( &bPRS, 1, buf, bsize, pos );
    int beginFlag = (anaManager_.getBeginningIntegrationFlag())?1:0;
    comm->pack( &beginFlag, 1, buf, bsize, pos );
  }
  else
  {
    // count here will be the size for the base N_TIA_StepErrorControl 
    // class *only*.
    int count = N_TIA_StepErrorControl::restartDataSize( false );
    int startIndex = pos;

    // Clobber any data in buf lest we leave garbage
    for( int i = startIndex; i < (startIndex+count); ++i) buf[i] = ' ';

    ostringstream ost;
    ost.width(24);ost.precision(16);ost.setf(ios::scientific);
    ost << startingTimeStep << " ";
    ost << currentTimeStep << " ";
    ost << lastAttemptedTimeStep << " ";
    ost << lastTimeStep << " ";
    ost << minTimeStep << " ";
    ost << maxTimeStep << " ";
    ost << maxTimeStepUser << " ";
    ost << lastTime << " ";
    ost << currentTime << " ";
    ost << nextTime << " ";
    ost << initialTime << " ";
    ost << currentTimeStepRatio << " ";
    ost << currentTimeStepSum << " ";
    ost << lastTimeStepRatio << " ";
    ost << lastTimeStepSum << " ";
    ost << newtonConvergenceStatus << " ";
    ost << numberSuccessiveFailures << " ";
    int flag = (stepAttemptStatus)?1:0;
    ost << flag << " ";
    ost << minStepPrecisionFac_ << " ";
    ost << newtonStepReduction_ << " ";
    ost << tolAimFac_ << " ";
    ost << estOverTol_ << " ";
    // Subtract one because we won't write out the pause breakpoint at the
    // final time
    int size = breakPoints_.size() - 1;
    ost << size << " ";

    set<N_UTL_BreakPoint>::iterator bpStart = breakPoints_.begin();
    set<N_UTL_BreakPoint>::iterator bpEnd = breakPoints_.end();
    for( set<N_UTL_BreakPoint>::iterator iterSD = bpStart;
         iterSD != bpEnd; ++iterSD )
    {
      if (!(iterSD->bptype() == PAUSE_BREAKPOINT && 
            iterSD->value() == finalTime))
      {
        ost << iterSD->value() << " ";
        ost << iterSD->bptype() << " ";
      }
    }      
    int im = anaManager_.getIntegrationMethod();
    ost << im << " ";
    int sN = anaManager_.getStepNumber();
    ost << sN << " ";
    int tSN = anaManager_.getTranStepNumber();
    ost << tSN << " ";
    int bPRS = anaManager_.breakPointRestartStep;
    ost << bPRS << " ";
    int beginFlag = (anaManager_.getBeginningIntegrationFlag())?1:0;
    ost << beginFlag << " ";

    string data( ost.str() );
    for( unsigned int i = 0; i < data.length(); ++i ) buf[startIndex+i] = data[i];
    // The line above copies the characters of the data string into buf,
    // but doesn't null-terminate buf.
    // it is essential to terminate the buffer with a null, or attempts
    // to construct a string object from it will get memory access problems.
    buf[startIndex+data.length()] = '\0';
    pos += data.length();

#ifdef Xyce_DEBUG_RESTART
    string outputString(buf);

    cout << "StepErrorControl UNPACKED output buffer:" << endl;
    cout << outputString << endl;
#endif
  }

#ifdef Xyce_DEBUG_RESTART
  cout << "TIA Restart Data DUMP (DAE)!  " << netListFile << "\n";
  cout << "----------------------------------------\n"<<endl;
  cout << "alphas_ = " <<alphas_<<endl;
  cout << "alpha0_ = " <<alpha0_<<endl;
  cout << "cj_ = " <<cj_<<endl;
  cout << "ck_ = " <<ck_<<endl;
  cout << "usedStep_ = " <<usedStep_<<endl;
  cout << "Ek_ = " <<Ek_<<endl;
  cout << "Ekm1_ = " <<Ekm1_<<endl;
  cout << "Ekm2_ = " <<Ekm2_<<endl;
  cout << "Ekp1_ = " <<Ekp1_<<endl;
  cout << "Est_ = " <<Est_<<endl;
  cout << "Tk_ = " <<Tk_<<endl;
  cout << "Tkm1_ = " <<Tkm1_<<endl;
  cout << "Tkm2_ = " <<Tkm2_<<endl;
  cout << "Tkp1_ = " <<Tkp1_<<endl;
  cout << "h0_safety_ = " <<h0_safety_<<endl;
  cout << "h0_max_factor_ = " <<h0_max_factor_<<endl;
  cout << "h_phase0_incr_ = " <<h_phase0_incr_<<endl;
  cout << "h_max_inv_ = " <<h_max_inv_<<endl;
  cout << "Tkm1_Tk_safety_ = " <<Tkm1_Tk_safety_<<endl;
  cout << "Tkp1_Tk_safety_ = " <<Tkp1_Tk_safety_<<endl;
  cout << "r_factor_ = " <<r_factor_<<endl;
  cout << "r_safety_ = " <<r_safety_<<endl;
  cout << "r_fudge_ = " <<r_fudge_<<endl;
  cout << "r_min_ = " <<r_min_<<endl;
  cout << "r_max_ = " <<r_max_<<endl;
  cout << "r_hincr_test_ = " <<r_hincr_test_<<endl;
  cout << "r_hincr_ = " <<r_hincr_<<endl;

  for (int i=0;i<6;++i)
  {
    cout << "  alpha_["<<i<<"] = " <<  alpha_[i]<<endl;
    cout << "  sigma_["<<i<<"] = " <<  sigma_[i]<<endl;
    cout << "  gamma_["<<i<<"] = " <<  gamma_[i]<<endl;
    cout << "  beta_["<<i<<"] = " <<  beta_[i]<<endl;
    cout << "  psi_["<<i<<"] = " <<  psi_[i]<<endl;
  }

  cout << "----------------------------------------\n\n"<<endl;
#endif

  if( pack )
  {
    // doubles:
    comm->pack( &alphas_  , 1, buf, bsize, pos );
    comm->pack( &alpha0_  , 1, buf, bsize, pos );
    comm->pack( &cj_  , 1, buf, bsize, pos );
    comm->pack( &ck_  , 1, buf, bsize, pos );
    comm->pack( &usedStep_  , 1, buf, bsize, pos );
    comm->pack( &Ek_  , 1, buf, bsize, pos );
    comm->pack( &Ekm1_  , 1, buf, bsize, pos );
    comm->pack( &Ekm2_  , 1, buf, bsize, pos );
    comm->pack( &Ekp1_  , 1, buf, bsize, pos );
    comm->pack( &Est_ , 1, buf, bsize, pos );
    comm->pack( &Tk_  , 1, buf, bsize, pos );
    comm->pack( &Tkm1_  , 1, buf, bsize, pos );
    comm->pack( &Tkm2_  , 1, buf, bsize, pos );
    comm->pack( &Tkp1_  , 1, buf, bsize, pos );
    comm->pack( &h0_safety_ , 1, buf, bsize, pos );
    comm->pack( &h0_max_factor_ , 1, buf, bsize, pos );
    comm->pack( &h_phase0_incr_ , 1, buf, bsize, pos );
    comm->pack( &h_max_inv_ , 1, buf, bsize, pos );
    comm->pack( &Tkm1_Tk_safety_  , 1, buf, bsize, pos );
    comm->pack( &Tkp1_Tk_safety_  , 1, buf, bsize, pos );
    comm->pack( &r_factor_  , 1, buf, bsize, pos );
    comm->pack( &r_safety_  , 1, buf, bsize, pos );
    comm->pack( &r_fudge_ , 1, buf, bsize, pos );
    comm->pack( &r_min_ , 1, buf, bsize, pos );
    comm->pack( &r_max_ , 1, buf, bsize, pos );
    comm->pack( &r_hincr_test_  , 1, buf, bsize, pos );
    comm->pack( &r_hincr_ , 1, buf, bsize, pos );

    // vectors of doubles:
    for (int i=0;i<6;++i)
    {
      comm->pack( &alpha_[i], 1, buf, bsize, pos );
      comm->pack( &sigma_[i], 1, buf, bsize, pos );
      comm->pack( &gamma_[i], 1, buf, bsize, pos );
      comm->pack( &beta_[i], 1, buf, bsize, pos );
      comm->pack( &psi_[i], 1, buf, bsize, pos );
    }

    // ints:
    comm->pack ( &currentOrder_, 1, buf, bsize, pos);
    comm->pack ( &oldOrder_, 1, buf, bsize, pos);
    comm->pack ( &maxOrder_, 1, buf, bsize, pos);
    comm->pack ( &minOrder_, 1, buf, bsize, pos);
    comm->pack ( &usedOrder_, 1, buf, bsize, pos);
    comm->pack ( &numberOfSteps_, 1, buf, bsize, pos);
    comm->pack ( &nef_, 1, buf, bsize, pos);
    comm->pack ( &nscsco_, 1, buf, bsize, pos);
    comm->pack ( &newOrder_, 1, buf, bsize, pos);
    comm->pack ( &max_LET_fail_, 1, buf, bsize, pos);

    // bools:
    int iP = (initialPhase_)?1:0;
    comm->pack ( &iP, 1, buf, bsize, pos);
  }
  else
  {
    ostringstream ost;
    ost.width(24);ost.precision(16);ost.setf(ios::scientific);

    // doubles:
    ost << alphas_  << " ";
    ost << alpha0_  << " ";
    ost << cj_  << " ";
    ost << ck_  << " ";
    ost << usedStep_  << " ";
    ost << Ek_  << " ";
    ost << Ekm1_  << " ";
    ost << Ekm2_  << " ";
    ost << Ekp1_  << " ";
    ost << Est_ << " ";
    ost << Tk_  << " ";
    ost << Tkm1_  << " ";
    ost << Tkm2_  << " ";
    ost << Tkp1_  << " ";
    ost << h0_safety_ << " ";
    ost << h0_max_factor_ << " ";
    ost << h_phase0_incr_ << " ";
    ost << h_max_inv_ << " ";
    ost << Tkm1_Tk_safety_  << " ";
    ost << Tkp1_Tk_safety_  << " ";
    ost << r_factor_  << " ";
    ost << r_safety_  << " ";
    ost << r_fudge_ << " ";
    ost << r_min_ << " ";
    ost << r_max_ << " ";
    ost << r_hincr_test_  << " ";
    ost << r_hincr_ << " ";

    // vectors of doubles:
    for (int i=0;i<6;++i)
    {
      ost << alpha_[i]<< " " ;
      ost << sigma_[i]<< " " ;
      ost << gamma_[i]<< " " ;
      ost << beta_[i]<< " " ;
      ost << psi_[i]<< " " ;
    }

    // ints:
    ost << currentOrder_ << " ";
    ost << oldOrder_ << " ";
    ost << maxOrder_ << " ";
    ost << minOrder_ << " ";
    ost << usedOrder_ << " ";
    ost << numberOfSteps_ << " ";
    ost << nef_ << " ";
    ost << nscsco_ << " ";
    ost << newOrder_ << " ";
    ost << max_LET_fail_ << " ";

    // bools:
    int iP = (initialPhase_)?1:0;
    ost << iP << " ";

    string data( ost.str() );
    for( unsigned int i = 0; i < data.length(); ++i ) buf[pos+i] = data[i];
    // null terminate buf, essential if buf is used as a string anywhere,
    // including in the string constructor below:
    buf[pos+data.length()] = '\0';

#ifdef Xyce_DEBUG_RESTART
    string outputString(buf);

    cout << "StepErrorControlDAE  UNPACKED output buffer:" << endl;
    cout << outputString << endl;
#endif
    pos = newPos;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::restoreRestartData
// Purpose       : 
// Special Notes : ERK:  6/20/2010:  
//                 There are no longer 2 distinct classes.  The DAE version is
//                 now part of the original, as one single class.  The two 
//                 blocks of code in this function correspond to the original
//                 class (first block), and the DAE class (second block).
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 7/28/06
//-----------------------------------------------------------------------------
bool N_TIA_StepErrorControl::restoreRestartData
  (char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{

  // original class variables:
  if( pack )
  {
    comm->unpack(buf, bsize, pos, &startingTimeStep, 1);
    comm->unpack(buf, bsize, pos, &currentTimeStep, 1);
    comm->unpack(buf, bsize, pos, &lastAttemptedTimeStep, 1);
    comm->unpack(buf, bsize, pos, &lastTimeStep, 1);
    comm->unpack(buf, bsize, pos, &minTimeStep, 1);
    comm->unpack(buf, bsize, pos, &maxTimeStep, 1);
    comm->unpack(buf, bsize, pos, &maxTimeStepUser, 1);
    comm->unpack(buf, bsize, pos, &lastTime, 1);
    comm->unpack(buf, bsize, pos, &currentTime, 1);
    comm->unpack(buf, bsize, pos, &nextTime, 1);
    comm->unpack(buf, bsize, pos, &initialTime, 1);
    comm->unpack(buf, bsize, pos, &currentTimeStepRatio, 1);
    comm->unpack(buf, bsize, pos, &currentTimeStepSum, 1);
    comm->unpack(buf, bsize, pos, &lastTimeStepRatio, 1);
    comm->unpack(buf, bsize, pos, &lastTimeStepSum, 1);
    comm->unpack(buf, bsize, pos, &newtonConvergenceStatus, 1);
    comm->unpack(buf, bsize, pos, &numberSuccessiveFailures, 1);
    int flag;
    comm->unpack(buf, bsize, pos, &flag, 1);
    stepAttemptStatus = flag;
    comm->unpack(buf, bsize, pos, &minStepPrecisionFac_, 1);
    comm->unpack(buf, bsize, pos, &newtonStepReduction_, 1);
    comm->unpack(buf, bsize, pos, &tolAimFac_, 1);
    comm->unpack(buf, bsize, pos, &estOverTol_, 1);

    double bpTol = 2.0 * minTimeStep;
    anaManager_.setBreakpointTol(bpTol);
    tiaParams_.initialTime=initialTime;
    tiaParams_.pauseTime=initialTime;
    set<N_UTL_BreakPoint> tmpSet = breakPoints_;
    breakPoints_.clear();
    set<N_UTL_BreakPoint>::iterator iterSD;
    set<N_UTL_BreakPoint>::iterator firstSD = tmpSet.begin();
    set<N_UTL_BreakPoint>::iterator lastSD = tmpSet.end();
    for (iterSD = firstSD; iterSD != lastSD; ++iterSD)
    {
      if (iterSD->value() > currentTime) 
      {
        breakPoints_.insert(*iterSD);
        if (iterSD->bptype() == PAUSE_BREAKPOINT)
          updatePauseTime(*iterSD);
      }
    }

    int size;
    double val;
    int bptype;
    N_UTL_BreakPoint TmpBP;
    comm->unpack(buf, bsize, pos, &size, 1);
    
    for (int i = 0; i < size; ++i)
    {
      comm->unpack(buf, bsize, pos, &val, 1);
      comm->unpack(buf, bsize, pos, &bptype, 1);
      if (val > currentTime)
      {
        TmpBP.set(val,bptype);
        breakPoints_.insert(TmpBP);
        if (TmpBP.bptype() == PAUSE_BREAKPOINT)
        {
          updatePauseTime(TmpBP);
        }
      }
    }

    int im;
    comm->unpack(buf, bsize, pos, &im, 1);
    anaManager_.setIntegrationMethod(im);
    comm->unpack(buf, bsize, pos, &im, 1);
    anaManager_.setStepNumber(im);
    comm->unpack(buf, bsize, pos, &im, 1);
    anaManager_.setTranStepNumber(im);
    comm->unpack(buf, bsize, pos, &im, 1);
    anaManager_.breakPointRestartStep = im;
    comm->unpack(buf, bsize, pos, &im, 1);
    anaManager_.setBeginningIntegrationFlag(im==1);
  }
  else
  {
    string str1(buf);
    int length = str1.size() - pos;
    string str2(str1,pos,length);

    istringstream ist( str2 );

    // count here will be the size for the base N_TIA_StepErrorControl 
    // class *only*.

    // THIS IS INCORRECT!  restartDataSize returns the MAXIMUM the thing
    // is allowed to be, which allows for 24 characters for each integer
    // value in the output.  This is almost always an overestimate, and 
    // using it as a way of pointing to the next character after our
    // data is a sure-fire way to break downstream processing!
    // Instead, we'll use the tellg function from the stringstream to 
    // return the offset after we read everything out.
    //    int count = N_TIA_StepErrorControl::restartDataSize( false );
    //    pos += count;

    ist >> startingTimeStep;
    ist >> currentTimeStep;
    ist >> lastAttemptedTimeStep;
    ist >> lastTimeStep;
    ist >> minTimeStep;
    ist >> maxTimeStep;
    ist >> maxTimeStepUser;
    ist >> lastTime;
    ist >> currentTime;
    ist >> nextTime;
    ist >> initialTime;
    ist >> currentTimeStepRatio;
    ist >> currentTimeStepSum;
    ist >> lastTimeStepRatio;
    ist >> lastTimeStepSum;
    ist >> newtonConvergenceStatus;
    ist >> numberSuccessiveFailures;
    int flag;
    ist >> flag;
    stepAttemptStatus = flag;
    ist >> minStepPrecisionFac_;
    ist >> newtonStepReduction_;
    ist >> tolAimFac_;
    ist >> estOverTol_;

    double bpTol = 2.0 * minTimeStep;
    anaManager_.setBreakpointTol(bpTol);
    tiaParams_.initialTime=initialTime;
    tiaParams_.pauseTime=initialTime;

    set<N_UTL_BreakPoint> tmpSet = breakPoints_;
    breakPoints_.clear();
    set<N_UTL_BreakPoint>::iterator iterSD;
    set<N_UTL_BreakPoint>::iterator firstSD = tmpSet.begin();
    set<N_UTL_BreakPoint>::iterator lastSD = tmpSet.end();
    for (iterSD = firstSD; iterSD != lastSD; ++iterSD)
    {
      if (iterSD->value() > currentTime)
      { 
        breakPoints_.insert(*iterSD);
        if (iterSD->bptype() == PAUSE_BREAKPOINT)
          updatePauseTime(*iterSD);
      }
    }

    int size;
    double val;
    int bptype;
    N_UTL_BreakPoint TmpBP;
    ist >> size;

    for( int i = 0; i < size; ++i )
    {
      ist >> val;
      ist >> bptype;
      if (val > currentTime)
        {
          TmpBP.set(val,bptype);
          breakPoints_.insert( TmpBP );
          if (TmpBP.bptype() == PAUSE_BREAKPOINT)
          {
            updatePauseTime(TmpBP);
          }
        }
    }

    int tmpInt;
    ist >> tmpInt;
    anaManager_.setIntegrationMethod(tmpInt);
    ist >> tmpInt;
    anaManager_.setStepNumber(tmpInt);
    ist >> tmpInt;
    anaManager_.setTranStepNumber(tmpInt);
    ist >> tmpInt;
    anaManager_.breakPointRestartStep = tmpInt;
    ist >> tmpInt;
    anaManager_.setBeginningIntegrationFlag(tmpInt==1);

    pos += ist.tellg();
  }

#ifdef Xyce_DEBUG_RESTART
  string netListFile = commandLine_.getArgumentValue("netlist");
  cout << "TIA Restart Data RESTORE!  " << netListFile << "\n";
  cout << "----------------------------------------\n";
  cout << "startingTimeStep: " << startingTimeStep << endl;
  cout << "currentTimeStep: " << currentTimeStep << endl;
  cout << "lastAttemptedTimeStep: " << lastAttemptedTimeStep << endl;
  cout << "lastTimeStep: " << lastTimeStep << endl;
  cout << "minTimeStep: " << minTimeStep << endl;
  cout << "maxTimeStep: " << maxTimeStep << endl;
  cout << "maxTimeStepUser: " << maxTimeStepUser << endl;
  cout << "lastTime: " << lastTime << endl;
  cout << "currentTime: " << currentTime << endl;
  cout << "nextTime: " << nextTime << endl;
  cout << "initialTime: " << initialTime << endl;
  cout << "estOverTol_: " << estOverTol_ << endl;
  cout << "breakpts: ";
  for (set<N_UTL_BreakPoint>::iterator iterSD = breakPoints_.begin();
       iterSD != breakPoints_.end(); ++iterSD)
  {
    cout << iterSD->value() << " " << iterSD->bptype() << endl;
  }
  cout << endl;
  cout << "integMethod: " << anaManager_.getIntegrationMethod() << endl;
  cout << "stepNumber: " << anaManager_.getStepNumber() << endl;
  cout << "tranStepNumber: " << anaManager_.getTranStepNumber() << endl;
  cout << "breakPointRestartStep: " << anaManager_.breakPointRestartStep <<
    endl;
  cout << "----------------------------------------\n\n";
#endif

  // DAE class variables:
  if( pack )
  {
    // doubles:
    comm->unpack(buf, bsize, pos, &alphas_  , 1);
    comm->unpack(buf, bsize, pos, &alpha0_  , 1);
    comm->unpack(buf, bsize, pos, &cj_  , 1);
    comm->unpack(buf, bsize, pos, &ck_  , 1);
    comm->unpack(buf, bsize, pos, &usedStep_  , 1);
    comm->unpack(buf, bsize, pos, &Ek_  , 1);
    comm->unpack(buf, bsize, pos, &Ekm1_  , 1);
    comm->unpack(buf, bsize, pos, &Ekm2_  , 1);
    comm->unpack(buf, bsize, pos, &Ekp1_  , 1);
    comm->unpack(buf, bsize, pos, &Est_ , 1);
    comm->unpack(buf, bsize, pos, &Tk_  , 1);
    comm->unpack(buf, bsize, pos, &Tkm1_  , 1);
    comm->unpack(buf, bsize, pos, &Tkm2_  , 1);
    comm->unpack(buf, bsize, pos, &Tkp1_  , 1);
    comm->unpack(buf, bsize, pos, &h0_safety_ , 1);
    comm->unpack(buf, bsize, pos, &h0_max_factor_ , 1);
    comm->unpack(buf, bsize, pos, &h_phase0_incr_ , 1);
    comm->unpack(buf, bsize, pos, &h_max_inv_ , 1);
    comm->unpack(buf, bsize, pos, &Tkm1_Tk_safety_  , 1);
    comm->unpack(buf, bsize, pos, &Tkp1_Tk_safety_  , 1);
    comm->unpack(buf, bsize, pos, &r_factor_  , 1);
    comm->unpack(buf, bsize, pos, &r_safety_  , 1);
    comm->unpack(buf, bsize, pos, &r_fudge_ , 1);
    comm->unpack(buf, bsize, pos, &r_min_ , 1);
    comm->unpack(buf, bsize, pos, &r_max_ , 1);
    comm->unpack(buf, bsize, pos, &r_hincr_test_  , 1);
    comm->unpack(buf, bsize, pos, &r_hincr_ , 1);

    // vectors of doubles:
    for (int i=0;i<6;++i)
    {
      comm->unpack(buf, bsize, pos, &alpha_[i], 1);
      comm->unpack(buf, bsize, pos, &sigma_[i], 1);
      comm->unpack(buf, bsize, pos, &gamma_[i], 1);
      comm->unpack(buf, bsize, pos, &beta_[i], 1);
      comm->unpack(buf, bsize, pos, &psi_[i], 1);
    }

    // ints:
    comm->unpack(buf, bsize, pos, &currentOrder_, 1);
    comm->unpack(buf, bsize, pos, &oldOrder_, 1);
    comm->unpack(buf, bsize, pos, &maxOrder_, 1);
    comm->unpack(buf, bsize, pos, &minOrder_, 1);
    comm->unpack(buf, bsize, pos, &usedOrder_, 1);
    comm->unpack(buf, bsize, pos, &numberOfSteps_, 1);
    comm->unpack(buf, bsize, pos, &nef_, 1);
    comm->unpack(buf, bsize, pos, &nscsco_, 1);
    comm->unpack(buf, bsize, pos, &newOrder_, 1);
    comm->unpack(buf, bsize, pos, &max_LET_fail_, 1);

    // bools:
    int iP;
    comm->unpack (buf, bsize, pos, &iP, 1);
    if (iP == 0) initialPhase_ = false;
    else         initialPhase_ = true;
  }
  else
  {
    // want the string stream to only represent the new-DAE part of 
    // the buf array.
    string str1(buf);
    int length = str1.size() - pos;
    string str2(str1,pos,length);

    istringstream ist( str2 );

    int count = N_TIA_StepErrorControl::restartDataSize(false);
    pos = count;  // can update here, as pos isn't used after this.

    // doubles:
    ist >> alphas_;
    ist >> alpha0_;
    ist >> cj_;
    ist >> ck_;
    ist >> usedStep_;
    ist >> Ek_;
    ist >> Ekm1_;
    ist >> Ekm2_;
    ist >> Ekp1_;
    ist >> Est_;
    ist >> Tk_;
    ist >> Tkm1_;
    ist >> Tkm2_;
    ist >> Tkp1_;
    ist >> h0_safety_;
    ist >> h0_max_factor_;
    ist >> h_phase0_incr_;
    ist >> h_max_inv_;
    ist >> Tkm1_Tk_safety_;
    ist >> Tkp1_Tk_safety_;
    ist >> r_factor_;
    ist >> r_safety_;
    ist >> r_fudge_;
    ist >> r_min_;
    ist >> r_max_;
    ist >> r_hincr_test_;
    ist >> r_hincr_;

    // vectors of doubles:
    for (int i=0;i<6;++i)
    {
      ist >> alpha_[i];
      ist >> sigma_[i];
      ist >> gamma_[i];
      ist >> beta_[i];
      ist >> psi_[i];
    }

    // ints:
    ist >> currentOrder_;
    ist >> oldOrder_;
    ist >> maxOrder_;
    ist >> minOrder_;
    ist >> usedOrder_;
    ist >> numberOfSteps_;
    ist >> nef_;
    ist >> nscsco_;
    ist >> newOrder_;
    ist >> max_LET_fail_;

    // bools:
    int iP;
    ist >> iP;
    if (iP == 0) initialPhase_ = false;
    else         initialPhase_ = true;
  }

#ifdef Xyce_DEBUG_RESTART
  cout << "TIA Restart Data RESTORE (DAE)!  " << netListFile << "\n";
  cout << "----------------------------------------\n"<<endl;
  cout << "alphas_ = " <<alphas_<<endl;
  cout << "alpha0_ = " <<alpha0_<<endl;
  cout << "cj_ = " <<cj_<<endl;
  cout << "ck_ = " <<ck_<<endl;
  cout << "usedStep_ = " <<usedStep_<<endl;
  cout << "Ek_ = " <<Ek_<<endl;
  cout << "Ekm1_ = " <<Ekm1_<<endl;
  cout << "Ekm2_ = " <<Ekm2_<<endl;
  cout << "Ekp1_ = " <<Ekp1_<<endl;
  cout << "Est_ = " <<Est_<<endl;
  cout << "Tk_ = " <<Tk_<<endl;
  cout << "Tkm1_ = " <<Tkm1_<<endl;
  cout << "Tkm2_ = " <<Tkm2_<<endl;
  cout << "Tkp1_ = " <<Tkp1_<<endl;
  cout << "h0_safety_ = " <<h0_safety_<<endl;
  cout << "h0_max_factor_ = " <<h0_max_factor_<<endl;
  cout << "h_phase0_incr_ = " <<h_phase0_incr_<<endl;
  cout << "h_max_inv_ = " <<h_max_inv_<<endl;
  cout << "Tkm1_Tk_safety_ = " <<Tkm1_Tk_safety_<<endl;
  cout << "Tkp1_Tk_safety_ = " <<Tkp1_Tk_safety_<<endl;
  cout << "r_factor_ = " <<r_factor_<<endl;
  cout << "r_safety_ = " <<r_safety_<<endl;
  cout << "r_fudge_ = " <<r_fudge_<<endl;
  cout << "r_min_ = " <<r_min_<<endl;
  cout << "r_max_ = " <<r_max_<<endl;
  cout << "r_hincr_test_ = " <<r_hincr_test_<<endl;
  cout << "r_hincr_ = " <<r_hincr_<<endl;

  for (int i=0;i<6;++i)
  {
    cout << "  alpha_["<<i<<"] = " <<  alpha_[i]<<endl;
    cout << "  sigma_["<<i<<"] = " <<  sigma_[i]<<endl;
    cout << "  gamma_["<<i<<"] = " <<  gamma_[i]<<endl;
    cout << "  beta_["<<i<<"] = " <<  beta_[i]<<endl;
    cout << "  psi_["<<i<<"] = " <<  psi_[i]<<endl;
  }

  cout << "----------------------------------------\n\n"<<endl;
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::updateTwoLevelTimeInfo
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel ComputationalSciences.
// Creation Date : 3/14/06
//-----------------------------------------------------------------------------
void N_TIA_StepErrorControl::updateTwoLevelTimeInfo (const N_TIA_TimeIntInfo & tiInfo)
{
  updateStopTime();

  // Need to get this right.
  if (previousCallStepSuccessful == true)
  {
    lastTimeStep      = currentTimeStep;
    lastTimeStepRatio = currentTimeStepRatio;
    lastTimeStepSum   = currentTimeStepSum;
    previousCallStepSuccessful = false;
  }

  // If the step needs to be adjusted:
  lastAttemptedTimeStep = currentTimeStep;
  double newTimeStep = tiInfo.nextTimeStep;
  nextTime = tiInfo.nextTime;
  currentTimeStepRatio = newTimeStep/lastTimeStep;
  currentTimeStepSum   = newTimeStep + lastTimeStep;
  currentTimeStep = newTimeStep;
  currentOrder_ = tiInfo.currentOrder;
}

#ifdef Xyce_VERBOSE_TIME
//-----------------------------------------------------------------------------
// Function      : N_TIA_StepErrorControl::outputTimeInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 03/04/08
//-----------------------------------------------------------------------------

void N_TIA_StepErrorControl::outputTimeInfo()
{
  ostringstream ost;
  streamsize width=ost.width();
  streamsize precision=ost.precision();
  ios_base::fmtflags ff=ost.flags();

#ifdef Xyce_DEBUG_TIME
  string netListFile = commandLine_.getArgumentValue("netlist");
  ost << " " << netListFile  << " ";
#else
  ost << "     ";
#endif

  ost << "Current,Next,Step: ";

#ifdef Xyce_DEBUG_TIME
  ost.width(14);    ost.precision(7); ost.setf(ios::scientific); 
#else
  ost.width(16);    ost.precision(9); ost.setf(ios::scientific); 
#endif

  //ost << currentTime << ", " << nextTime << ", " << currentTimeStep << endl;
  ost << currentTime << ", " << nextTime << ", " << currentTimeStep;

  // Restore width,precision and flags
  ost.width(width);    ost.precision(precision); ost.flags(ff);

  ost << "    ("<<numberOfSteps_<<") Used, Next Order:  ";
  //ost << usedOrder_ << ", " << currentOrder_ << endl;
  ost << usedOrder_ << ", " << currentOrder_;

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, ost.str());
}
#endif
