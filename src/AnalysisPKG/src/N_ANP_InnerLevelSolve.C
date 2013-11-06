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
// Filename      : $RCSfile: N_ANP_InnerLevelSolve.C,v $
// Purpose       : This file contains functions from control-algorithm that are
//                 used by the inner solve of a 2-level xyce-to-xyce solve.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 1/23/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.39.2.3 $
// Revision Date  : $Date: 2013/10/03 17:23:31 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#define Xyce_MPDE_IC
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

// ----------   Xyce Includes   ----------

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_Transient.h>
#include <N_ANP_DCSweep.h>

#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TwoLevelError.h>

#include <N_ERH_ErrorMgr.h>

#include <N_NLS_Manager.h>
#include <N_LOA_Loader.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>
#include <N_UTL_BreakPoint.h>

#include <N_IO_CmdParse.h>

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::provisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/04/2009
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::provisionalStep (double maxTimeStep,  double &timeStep)
{
  bool bsuccess = true;
  bool b1 = true;
  string msg;
  bool dcopFlag = true;

  if (!initializeSolvers_mixedSignal_)
  {
    if (analysis == ANP_MODE_TRANSIENT)
    {
      mixedSignalAnalysisObject_ = Teuchos::rcp(new N_ANP_Transient(this));
      mixedSignalAnalysisObject_->setAnalysisParams(tranParamsBlock);
      secPtr_->resetAll();
    }
    else if (analysis == ANP_MODE_DC_SWEEP)
    {
      mixedSignalAnalysisObject_ = Teuchos::rcp(new N_ANP_DCSweep(this));
      for (int i=0;i<dcParamsBlockVec.size();++i)
      {
        mixedSignalAnalysisObject_->setAnalysisParams(dcParamsBlockVec[i]);
      }
    }
    else
    {
      msg = "N_ANP_AnalysisManager::provisionalStep - "
        "unknown type of analysis\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
    mixedSignalAnalysisObject_->init();

    // This checks to make sure that all quantities in the .print are
    // valid before we continue.  This only needs to be done once.
    outputMgrAdapterRCPtr_->check_output( analysis ); //, *dsPtr_->currSolutionPtr, *dsPtr_->currStatePtr, *dsPtr_->currStorePtr );

    // Start the solvers timer.
    xyceTranTimerPtr_ = rcp(new N_UTL_Timer( *(pdsMgrPtr->getPDSComm()) ));

    initializeSolvers_mixedSignal_=true;
  }

  primaryAnalysisObject_ = mixedSignalAnalysisObject_;

  RefCountPtr<N_ANP_Transient> mixedSignalTransientAnalysisObject =
    Teuchos::rcp_dynamic_cast<N_ANP_Transient>(mixedSignalAnalysisObject_);

  if (!Teuchos::is_null(mixedSignalTransientAnalysisObject))
  {
    dcopFlag = mixedSignalTransientAnalysisObject->getDCOPFlag();
  }

  // Now save time step info, in case this step gets rejected.

  if (!(secPtr_->finished()))
  {
    bool stepSuccess = false;

    if (dcopFlag) // if dcop step, make one attempt.
    {
      mixedSignalAnalysisObject_->preStepDetails (maxTimeStep);
      b1 = mixedSignalAnalysisObject_->mixedSignalStep();

      // only call finalize step here if we have failed.
      if (!secPtr_->stepAttemptStatus)
      {
        mixedSignalAnalysisObject_->finalizeStep();
      }
      stepSuccess = secPtr_->stepAttemptStatus;
    }
    else // else, if transient step, keep re-taking the step
         // until it succeeds, or gets an unrecoverable failure,
         // such as time-step-too-small.
    {
      bool recoverableFailureFlag=true;
      while (!stepSuccess && recoverableFailureFlag)
      {
        mixedSignalAnalysisObject_->preStepDetails (maxTimeStep);
        b1 = mixedSignalAnalysisObject_->mixedSignalStep();

        // Only call finalize step here if step has failed.
        // If we succeed, we want to give Habanero the opportunity
        // to reject the step, after this function (provisionalStep)
        // exits.
        if (!secPtr_->stepAttemptStatus)
        {
          recoverableFailureFlag = mixedSignalAnalysisObject_->finalizeStep ();
        }
        else
        {
          stepSuccess = true;
        }
      }
    }
    bsuccess = stepSuccess;
  }

  // get the step information.
  //

  if (dcopFlag)
  {
    timeStep = 0.0;
  }
  else
  {
    N_TIA_TimeIntInfo tiInfo;
    getTimeIntInfo(tiInfo);
    timeStep = tiInfo.nextTimeStep;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::acceptProvisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::acceptProvisionalStep ()
{
  mixedSignalAnalysisObject_->finalizeStep ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::rejectProvisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::rejectProvisionalStep ()
{
  secPtr_->stepAttemptStatus = false;
  secPtr_->updateBreakPoints();

  bool dcopFlag = false;
  RefCountPtr<N_ANP_Transient> mixedSignalTransientAnalysisObject =
    Teuchos::rcp_dynamic_cast<N_ANP_Transient>(mixedSignalAnalysisObject_);

  if (!Teuchos::is_null(mixedSignalTransientAnalysisObject))
  {
    dcopFlag = mixedSignalTransientAnalysisObject->getDCOPFlag();
  }

  if (dcopFlag)
  {
    mixedSignalAnalysisObject_->finalizeStep ();
  }
  // Transient
  else
  {
     //   bool b1 = processFailedStep();
    loaderPtr->stepFailure (currentMode_);
    wimPtr->rejectStepForHabanero();

    mixedSignalTransientAnalysisObject->totalNumberFailedStepsAttempted_  += 1;
    secPtr_->numberSuccessiveFailures += 1;

  } // transient

#if  0
  if (secPtr_->isPauseTime())
  {
    // Failure at this point only indicates that the simulation
    // is paused and may be resumed.
    secPtr_->simulationPaused();
    isPaused = true;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::runStep
//
// Purpose       : This function is similar to "run" except that only a single
//                 integration (or DC sweep) step will be executed.
//
// Special Notes : Used for multi-level Newton solves, for levels other
//                 than the top level.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::runStep
    (const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError)
{
  string msg;

  if (initializeAllFlag_ == false)
  {
    msg = "N_ANP_AnalysisManager::runStep - "
      "you need to call the initializeAll function first\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  if (analysisParamsRegistered == false)
  {
    msg = "N_ANP_AnalysisManager::runStep: "
      "no analysis statement in the netlist\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  bool integration_status = false;

#ifdef Xyce_VERBOSE_TIME
  tiaParams.printParams(anpAnalysisModeToNLS(analysis));
#endif

  solverStartTime_ = elapsedTimerPtr_->elapsedTime();

  if (stepLoopFlag_)
  {
    msg = "N_ANP_AnalysisManager::runStep - "
      "Not valid to use .STEP statements in an inner solve.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  else
  {
    if (analysis == ANP_MODE_TRANSIENT)
    {
#ifdef Xyce_VERBOSE_TIME
      std::cout << "N_ANP_AnalysisManager::runStep:" << std::endl;
      std::cout << "nextTime = " << tiInfo.nextTime << std::endl;
      std::cout << "stepSize = " << tiInfo.nextTimeStep << std::endl;
#endif
    }
    else if (analysis == ANP_MODE_DC_SWEEP)
    {
      // do nothing
    }
    else
    {
      msg = "N_ANP_AnalysisManager::runStep - "
        "unknown type of analysis\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
    integration_status  = twoLevelAnalysisObject_->twoLevelStep();
    wimPtr->setupTwoLevelError(tlError);
  }

  return integration_status;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::startTimeStep_
// Purpose       : used by 2-level solves.
// Special Notes : One of the primary purposes for this function is to impose
//                 a lot of upper level information from the top level time
//                 integrator on the inner level time integrator.  This
//                 information is contained in the N_TIA_TimeIntInfo
//                 class (tiInfo here).  This includes things like the
//                 time step size, the time integration order, etc.
//
//                 In general, in a 2-level solve, an inner solver doesn't
//                 have enough information to correctly determine the
//                 step size, order, etc.  This is in part because the
//                 inner solver, while it knows about the top level solver,
//                 it cannot know about any OTHER inner solvers.
//
//                 The top level solver, however, does have enough information.
//                 It gathers break point, error analysis, and other info
//                 from all the inner solves.  So, the top level solver
//                 makes all the decisions and imposes them on the inner
//                 solves.  This function is where it does that impose.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::startTimeStep (const N_TIA_TimeIntInfo & tiInfo)
{
  // Beginning Integration flag is the only piece of data in tiInfo that is
  // currently owned by the control algorithm class.  Everything else
  // is owned by the step error control, or the integration method.
  twoLevelAnalysisObject_->setBeginningIntegrationFlag(tiInfo.beginIntegrationFlag);

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    std::cout << "N_ANP_AnalysisManager::startTimeStep:" << std::endl;
  }
#endif

#ifdef Xyce_VERBOSE_TIME
  if( !is_null(twoLevelAnalysisObject_) )
  {
    // an MPDE run also traverses this area so catch case where this is null
    twoLevelAnalysisObject_->printStepHeader();
  }
#endif
  if (switchIntegrator_) wimPtr->createTimeIntegMethod( twoLevelAnalysisObject_->getIntegrationMethod() );

  // ------------------------------------------------------------------------
  // Set the step size, current time and next time.
#ifndef Xyce_CHARON
  if ( twoLevelAnalysisObject_->getIntegrationMethod() != TIAMethod_NONE)
#endif
  {
    secPtr_->updateTwoLevelTimeInfo(tiInfo);
  }


  if (twoLevelAnalysisObject_->getBeginningIntegrationFlag() && secPtr_->stepAttemptStatus)
  {
    wimPtr->setTwoLevelTimeInfo(tiInfo); // new-dae only
  }

  // ------------------------------------------------------------------------
  // If we've switched the integration method, we need to obtain the
  // corrector derivative only after we've updated the TimeInfo.
  if (switchIntegrator_)
  {
    switchIntegrator_ = false;
    wimPtr->obtainCorrectorDeriv();
  }

  bool dcopFlag = true;
  RefCountPtr<N_ANP_Transient> twoLevelTransientAnalysisObject = Teuchos::rcp_dynamic_cast<N_ANP_Transient>(twoLevelAnalysisObject_);
  if (!Teuchos::is_null(twoLevelTransientAnalysisObject)) {
    dcopFlag = twoLevelTransientAnalysisObject->getDCOPFlag();
  }
#ifdef Xyce_VERBOSE_TIME
  if (!dcopFlag) secPtr_->outputTimeInfo();
#endif

  // ------------------------------------------------------------------------
  // Set the nonlinear solver parameters to those appropriate for the
  // transient solution, if neccessary.
  if (!dcopFlag)
  {
    nlsMgrPtr->setAnalysisMode(anpAnalysisModeToNLS(ANP_MODE_TRANSIENT));
  }

  // Ask the method to update its coefficients
  wimPtr->updateCoeffs(); 
  twoLevelAnalysisObject_->handlePredictor();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::conductanceTest
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::conductanceTest ()
{
  map<string,double> inputMap;
  vector<double> outputVector;
  vector< vector<double> > jacobian;

  // load inputMap from tiaParam.condTestDeviceNames option
  std::list< std::string >::iterator currentDeviceName = tiaParams.condTestDeviceNames.begin();
  std::list< std::string >::iterator endDeviceName = tiaParams.condTestDeviceNames.end();
  while( currentDeviceName != endDeviceName )
  {
    inputMap[ *currentDeviceName ] = 0.0;
    ++currentDeviceName;
  }

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    std::cout << "N_ANP_AnalysisManager::conductanceTest()" << std::endl;
    currentDeviceName = tiaParams.condTestDeviceNames.begin();
    while( currentDeviceName != endDeviceName )
    {
      std::cout << "currentDeviceName = \"" << *currentDeviceName << "\" added to inputMap[ "
            << *currentDeviceName << " ] = " << inputMap[ *currentDeviceName ] << std::endl;
      ++currentDeviceName;
    }
  }
#endif

  int isize=inputMap.size();
  outputVector.resize(isize,0.0);
  jacobian.resize(isize);
  for (int i=0;i<isize;++i)
  {
    jacobian[i].resize(isize,0.0);
  }

  bool b1 = nlsMgrPtr->obtainConductances(
        inputMap,
        outputVector,
        jacobian
    );

  int iE1, iE2;
  int numElectrodes = isize;

  FILE *fp1;
  fp1 = fopen("conductance.txt","w");

  fprintf(fp1, "%s", "Conductance array: \n");
  fprintf(fp1,"%s", "              ");
  if (b1)
  {
    map<string,double>::iterator iterM = inputMap.begin();
    map<string,double>::iterator  endM = inputMap.end  ();
    for (iE2 = 0; iE2 < numElectrodes; ++iE2,++iterM)
    {
      string srcname = iterM->first;
      fprintf(fp1,"\t%14s",srcname.c_str());
    }
    fprintf(fp1,"%s", "\n");

    iterM = inputMap.begin();
    for (iE1 = 0; iE1 < numElectrodes; ++iE1, ++iterM)
    {
      string srcname = iterM->first;
      fprintf(fp1,"%14s",srcname.c_str());
      for (iE2 = 0; iE2 < numElectrodes; ++iE2)
      {
        fprintf(fp1,"\t%14.4e",jacobian[iE1][iE2]);
      }
      fprintf(fp1,"%s", "\n");
    }
    fprintf(fp1,"%s", "\n");
  }
  else
  {
    fprintf(fp1,"%s", "\nConductance calculation failed!\n");
  }

  fclose(fp1);

}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::startupSolvers
// Purpose       :
// Special Notes : Used only for 2-level solves.
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 3/10/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::startupSolvers ()
{
  bool bsuccess = true;
  string msg;

  if (analysis == ANP_MODE_TRANSIENT)
  {
    nlsMgrPtr->resetAll(DC_OP);
    twoLevelAnalysisObject_ = Teuchos::rcp(new N_ANP_Transient(this));
    twoLevelAnalysisObject_->setAnalysisParams(tranParamsBlock);
    secPtr_->resetAll();
  }
  else if (analysis == ANP_MODE_DC_SWEEP)
  {
    secPtr_->setTIAParams();
    twoLevelAnalysisObject_ = Teuchos::rcp(new N_ANP_DCSweep(this));
    for (int i=0;i<dcParamsBlockVec.size();++i)
    {
      twoLevelAnalysisObject_->setAnalysisParams(dcParamsBlockVec[i]);
    }
  }
  else
  {
    msg = "N_ANP_AnalysisManager::startupSolvers: Multi-Level Newton solves only supports DC and Transient analysis\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  primaryAnalysisObject_ = twoLevelAnalysisObject_;
  twoLevelAnalysisObject_->init();

  // This checks to make sure that all quantities in the .print are
  // valid before we continue.  This only needs to be done once.
  outputMgrAdapterRCPtr_->check_output( analysis ); //, *dsPtr_->currSolutionPtr, *dsPtr_->currStatePtr, *dsPtr_->currStorePtr );


  // Start the solvers timer.
  xyceTranTimerPtr_ = rcp(new N_UTL_Timer( *(pdsMgrPtr->getPDSComm()) ));

  // Hardwire the erroption parameter to 1, which will force the inner
  // solve (initiated by this function) to only use the Newton success/failure
  // as step criteria.  Predictor-corrector information is handled in the
  // upper level solver.
  tiaParams.errorAnalysisOption = 1;
  // tscoffe 03/07/08  With the new errorAnalysisOption = 1 features of nlmin
  // and nlmax, we should enforce the original mode where nlmin=nlmax=maxNLSteps.
  //int maxNLSteps = ???;
  //tiaParams.NLmax = maxNLSteps;
  //tiaParams.NLmin = maxNLSteps;

 return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::finishSolvers
// Purpose       :
// Special Notes : Used only for 2-level solves.
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/10/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::finishSolvers ()
{
  bool bsuccess = true;
  string msg;

  twoLevelAnalysisObject_->finish();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::homotopyStepSuccess
// Purpose       : Lower-level processing of a successful homotopy step,
//                 which was controlled from the upper level of a 2-level solve.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::homotopyStepSuccess
    ( const vector<string> & paramNames,
      const vector<double> & paramVals)
{
#ifdef Xyce_DEBUG_ANALYSIS
  string netListFile = commandLine_.getArgumentValue("netlist");
  std::cout << "\n " << netListFile;
  std::cout << " N_ANP_AnalysisManager::homotopyStepSuccess " << std::endl;
#endif
  // output:
  outputMgrAdapterRCPtr_->outputHomotopy( paramNames, paramVals, *dsPtr_->nextSolutionPtr );

  // update the data arrays:
  dsPtr_->updateSolDataArrays();

  // pass info to the next level down, if it exists.
  loaderPtr->homotopyStepSuccess (paramNames,paramVals);

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::homotopyStepFailure
//
// Purpose       : Lower-level processing of a failed homotopy step,
//                 which was controlled from the upper level of a
//                 2-level solve.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::homotopyStepFailure ()
{
#ifdef Xyce_DEBUG_ANALYSIS
  string netListFile = commandLine_.getArgumentValue("netlist");
  std::cout << "\n " << netListFile;
  std::cout << " N_ANP_AnalysisManager::homotopyStepFailure " << std::endl;
#endif

  // The solutions currently in place represent failure.  Get rid of them.
  dsPtr_->usePreviousSolAsPredictor ();

  // pass info to the next level down, if it exists.
  loaderPtr->homotopyStepFailure ();

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::stepSuccess (int analysisUpper)
{

#ifdef Xyce_DEBUG_ANALYSIS
  string netListFile = commandLine_.getArgumentValue("netlist");
  std::cout << "\n " << netListFile;
  std::cout << " N_ANP_AnalysisManager::stepSuccess " << std::endl;
#endif

  currentMode_ = analysisUpper;
  secPtr_->stepAttemptStatus = true;
  switch( analysisUpper )
  {
    case 0:
      {
        N_ANP_Transient * twoLevelTransientAnalysisObject = dynamic_cast<N_ANP_Transient*>(&*twoLevelAnalysisObject_);
        if (twoLevelTransientAnalysisObject)
        {
          twoLevelTransientAnalysisObject->processSuccessfulDCOP();
        }
        else
        {
          string msg = "N_ANP_AnalysisManager::stepSuccess - "
            "Failed dynamic_cast of twoLevelAnalysisObject to N_ANP_Transient.\n";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
        }
      }
      break;
    case 1:
      twoLevelAnalysisObject_->processSuccessfulStep();
      break;
    case 2:
      twoLevelAnalysisObject_->processSuccessfulStep();
      break;
    default:
      string msg = "N_ANP_AnalysisManager::stepSuccess - "
        "unknown type of analysis\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void N_ANP_AnalysisManager::stepFailure (int analysisUpper)
{
  currentMode_ = analysisUpper;
  secPtr_->stepAttemptStatus = false;
  switch( analysisUpper )
  {
    case 0:
      {
        N_ANP_Transient * twoLevelTransientAnalysisObject = dynamic_cast<N_ANP_Transient*>(&*twoLevelAnalysisObject_);
        if (twoLevelTransientAnalysisObject)
        {
          twoLevelTransientAnalysisObject->processFailedDCOP();
        }
        else
        {
          string msg = "N_ANP_AnalysisManager::stepSuccess - "
            "Failed dynamic_cast of twoLevelAnalysisObject to N_ANP_Transient.\n";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
        }
      }
      break;
    case 1:
      twoLevelAnalysisObject_->processFailedStep();
      break;
    case 2:
      twoLevelAnalysisObject_->processFailedStep();
      break;
    default:
      string msg = "N_ANP_AnalysisManager::stepFailure - "
        "unknown type of analysis\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getInitialQnorm
// Purpose       : Used for 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  bool bsuccess = true;
  wimPtr->getInitialQnorm (tle);
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisManager::getBreakPoints
// Purpose       : Used for 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::getBreakPoints
     (vector<N_UTL_BreakPoint> &breakPointTimes)
{
  bool bsuccess = true;
  loaderPtr->getBreakPoints(breakPointTimes);
  return bsuccess;
}

