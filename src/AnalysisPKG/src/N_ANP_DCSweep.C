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
// Filename      : $RCSfile: N_ANP_DCSweep.C,v $
// Purpose       : DC Sweep class analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.30.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:31 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_HB.h>

#include <N_TIA_DataStore.h>
#include <N_LOA_Loader.h>
#include <N_MPDE_Manager.h>
#include <N_IO_CmdParse.h>

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/22/10
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::setAnalysisParams(const N_UTL_OptionBlock & paramsBlock)
{
#ifdef Xyce_DEBUG_ANALYSIS
  string dashedline = "-------------------------------------------------"
                       "----------------------------\n";
  if (tiaParams.debugLevel > 0)
  {
    std::cout << std::endl << dashedline << std::endl;
    std::cout << "N_ANP_DCSweep::setDCAnalysisParams" << std::endl;
  }
#endif

  list<N_UTL_Param>::const_iterator it_tp;
  list<N_UTL_Param>::const_iterator it_param;
  list<N_UTL_Param>::const_iterator it_type;
  list<N_UTL_Param>::const_iterator first = paramsBlock.getParams().begin();
  list<N_UTL_Param>::const_iterator last = paramsBlock.getParams().end();

  string msg;

  N_ANP_SweepParam sp;
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    for (it_tp = first; it_tp != last; ++it_tp)
    {
      std::cout << it_tp->uTag() ;
      std::cout << "\t";
      if (it_tp->uTag() == "PARAM" || it_tp->uTag() == "TYPE")
      {
        std::cout << it_tp->sVal ();
      }
      else
      {
        std::cout << it_tp->dVal ();
      }
      std::cout << std::endl;
    }
  }
#endif

  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag() == "TYPE")
    {
      it_type = it_tp;
      sp.type = it_tp->sVal();
    }

    if (it_tp->uTag() == "PARAM")
    {
      it_param = it_tp;
      sp.name = it_tp->sVal ();
    }
  }

  it_tp = it_param;
  ++it_tp;
  if (sp.type == "LIN") // default
  {
    sp.startVal = it_tp->dVal (); ++it_tp;
    sp.stopVal  = it_tp->dVal (); ++it_tp;
    sp.stepVal  = it_tp->dVal (); ++it_tp;
  }
  else if (sp.type == "DEC" || sp.type == "OCT")
  {
    sp.startVal = it_tp->dVal (); ++it_tp;
    sp.stopVal  = it_tp->dVal (); ++it_tp;
    sp.numSteps = it_tp->iVal (); ++it_tp;
  }
  else if (sp.type == "LIST")
  {
    for (;it_tp!=last;++it_tp)
    {
      sp.valList.push_back(it_tp->dVal());
    }
  }
  else
  {
    msg =  "N_ANP_DCSweep::setDCParams: ";
    msg += " unsupported DC type\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  // need to do a bunch of stuff to initialize the DC loop.
  (*dcParamVec_).push_back(sp);

  // set analysis type
  //analysis = ANP_MODE_DC_SWEEP;

  // good idea to put this here, but this code isn't touched in
  // no .dc runs and check_outputs will expcect the param vec
  // to exist (but be empty) even if this isn't a .dc run
  //outputMgrAdapterRCPtr_->setDCParamVec( & dcParamVec_ );

  //analysisParamsRegistered = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::outputFailureStats
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/21/2010
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::outputFailureStats ()
{
  if (!(dcSweepFailures_.empty()))
  {
    string msg = ("\tFailed DC sweep steps:\t\t");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, msg);
    list<int>::iterator first = dcSweepFailures_.begin();
    list<int>::iterator last  = dcSweepFailures_.end();
    list<int>::iterator iter;

    for (iter = first; iter != last; ++iter)
    {
      int itmp = *iter;
      msg = "\t\tDC Step # ";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, msg, itmp);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::run()
// Purpose       : This is the main controlling loop for DC sweep analysis.
//                 This loop calls a series of operating point calculations
//                 (ie calculations in which there is no time integration,
//                 and the time integration method is always set to "none").
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::run()
{
  bool bsuccess = true;
  bsuccess = init();
  bsuccess &= loopProcess();
  bsuccess &= finish();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::init()
{
  bool bsuccess = true;

  if (sensFlag_)
  {
    nlsMgrRCPtr_->enableSensitivity();
  }

  // set up the various sweep variables
  // the following should only happen once, but must be isolated to
  // prevent a step loop outside of this dc sweep from causing to
  // occur more often
  if( !dcLoopInitialized_ )
  {
#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams.debugLevel > 0)
    {
      std::cout << std::endl << std::endl;
      std::cout << "-----------------" << std::endl;
      std::cout << "N_ANP_DCSweep::run()" << std::endl;
    }
#endif

    dcLoopSize_ = setupSweepLoop_( *dcParamVec_ );

    outputMgrAdapterRCPtr_->setDCAnalysisMaxSteps( dcLoopSize_ );

    dcLoopInitialized_ = true;
  }

  anaManagerRCPtr_->currentMode_ = 2;

  //setup for operating pt calculation
  integrationMethod_ = TIAMethod_NONE;
  wimRCPtr_->createTimeIntegMethod(integrationMethod_);

  stepNumber = 0;
  doubleDCOPFlag_ = loaderRCPtr_->getDoubleDCOPFlag ();
  doubleDCOPStep_ = tiaParams.firstDCOPStep;

  initializeSolution_();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::initializeSolution()
// Purpose       : Move solution-initialization steps to separate routine so
//                 the process may be repeated cleanly when necessary.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystem Modeling
// Creation Date : 2/17/2010
//-----------------------------------------------------------------------------
void N_ANP_DCSweep::initializeSolution_()
{
  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loaderRCPtr_->setInitialGuess (dsRCPtr_->nextSolutionPtr);

  // If available, set initial solution (.IC, .NODESET, etc).
  inputOPFlag_ =
    outputMgrAdapterRCPtr_->setupInitialConditions( *(dsRCPtr_->nextSolutionPtr),
                        *(dsRCPtr_->flagSolutionPtr));

  // Set a constant history for operating point calculation
  dsRCPtr_->setConstantHistory();
  dsRCPtr_->computeDividedDifferences();
  wimRCPtr_->obtainCorrectorDeriv();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::loopProcess()
{
  bool bsuccess = true;
  // DC Sweep loop:

  int currentStep = 0;
  int finalStep = dcLoopSize_;
  while (currentStep < finalStep)
  {
    outputMgrAdapterRCPtr_->setDCAnalysisStepNumber( currentStep );

#ifdef Xyce_VERBOSE_TIME
    this->printStepHeader();
    secRCPtr_->outputTimeInfo();
#endif

    updateSweepParams_(currentStep, *dcParamVec_);
    if (currentStep != 0 &&(anaManagerRCPtr_->getSweepSourceResetFlag()))
    {
      dsRCPtr_->setZeroHistory();
      initializeSolution_();
    }

    // Perform the step:
    takeStep_ ();

    // Set things up for the next time step, based on if this one was
    // successful.
    if (secRCPtr_->stepAttemptStatus)
    {
      processSuccessfulStep();
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      processFailedStep();
    }

    // we don't really control the loop counter here.  processSuccessfulStep
    // and processFailedStep do that work.  This needs to be cleaned up. RLS
    currentStep = stepNumber;
  } // end of sweep loop

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::processSuccessfulStep()
//
// Purpose       : Used by both function dcSweepLoop and 2-level solve
//                 function calls to process successful DC steps.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::processSuccessfulStep()
{
  bool bsuccess = true;
  loaderRCPtr_->stepSuccess(anaManagerRCPtr_->currentMode_);

  // This output call is for device-specific output, such as .OP,
  // or internal plot output from PDE(TCAD) devices.  
  loaderRCPtr_->output();

  if (sensFlag_ && !firstDoubleDCOPStep_() )
  {
    nlsMgrRCPtr_->calcSensitivity();
  }

  // Do some statistics, as long as this isn't the first "double"
  // DCOP step. (that one doesn't count)
  if ( !firstDoubleDCOPStep_() )
  {
    stepNumber += 1;
    totalNumberSuccessStepsThisParameter_ += 1;
    totalNumberSuccessfulStepsTaken_ += 1;
  }

  // update the data arrays, output:
  dsRCPtr_->updateSolDataArrays();

#ifdef Xyce_DEBUG_ANALYSIS
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
    dsRCPtr_->outputSolDataArrays();
#endif
#endif

  dcSweepOutput();

  // now that output has been called, update the doubleDCOP step
  // if neccessary. (pde-only)
  doubleDCOPStep_ = tiaParams.lastDCOPStep;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::processFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::processFailedStep()
{
  loaderRCPtr_->stepFailure (anaManagerRCPtr_->currentMode_);

  stepNumber += 1;
  dcSweepFailures_.push_back(stepNumber);
  totalNumberFailedStepsAttempted_  += 1;
  secRCPtr_->numberSuccessiveFailures += 1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::finish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::finish()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  std::cout << "Calling N_ANP_DCSweep::finish() outputs!" << std::endl;
#endif

  outputMgrAdapterRCPtr_->finishOutput();
  if (!(dcSweepFailures_.empty()))
  {
    bsuccess = false;
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::handlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::handlePredictor()
{
  dsRCPtr_->setErrorWtVector();
  wimRCPtr_->obtainPredictor();
  wimRCPtr_->obtainPredictorDeriv();

  // In case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  loaderRCPtr_->startTimeStep ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::takeStep_
// Purpose       : Take a DC Sweep integration step.
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void N_ANP_DCSweep::takeStep_()
{
  handlePredictor();
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();
  dsRCPtr_->stepLinearCombo ();
  gatherStepStatistics_ ();
  secRCPtr_->evaluateStepError ();

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::twoLevelStep
//
// Purpose       : Used by 2-level Newton solves to execute a single DC sweep
//                 step.
//
// Special Notes : This is mostly what happens on the inner loop of
//                 dcSweepLoop, except that DC parameters are not updated,
//                 and success/failure of the step is not determined.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::twoLevelStep()
{
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();
  dsRCPtr_->stepLinearCombo ();
  gatherStepStatistics_ ();
  secRCPtr_->evaluateStepError ();

  return secRCPtr_->stepAttemptStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::dcSweepOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/31/08
//-----------------------------------------------------------------------------
void N_ANP_DCSweep::dcSweepOutput()
{
  // if this is an MPDE run, don't output anything.
  bool mpdeFlag = anaManagerRCPtr_->getMPDEFlag();
  bool hbFlag = anaManagerRCPtr_->getHBFlag();

  // Make sure this isn't the NLP step of a PDE DCOP.
  if ( !mpdeFlag  )
  { 
    if (!firstDoubleDCOPStep_())
    {
      outputMgrAdapterRCPtr_->dcOutput( 
          stepNumber, 
          *dsRCPtr_->currSolutionPtr, 
          *dsRCPtr_->currStatePtr, 
          *dsRCPtr_->currStorePtr );

      outputMgrAdapterRCPtr_->outputDCOP( *(dsRCPtr_->currSolutionPtr) );
    }
  }
  else
  {
    if ( hbFlag )
    {
      std::vector<double> timePoints, freqPoints;
      Teuchos::RefCountPtr<N_LAS_BlockVector> timeDomainSolnVec;
      Teuchos::RefCountPtr<N_LAS_BlockVector> freqDomainSolnVecReal;
      Teuchos::RefCountPtr<N_LAS_BlockVector> freqDomainSolnVecImaginary;
  
      Teuchos::RCP<const N_ANP_AnalysisBase> analysisObject 
        = anaManagerRCPtr_->getAnalysisObject();
 
      Teuchos::rcp_dynamic_cast<const N_ANP_HB>( analysisObject )->prepareHBOutput(
          *(dsRCPtr_->currSolutionPtr), 
          timePoints, 
          freqPoints, 
          timeDomainSolnVec, 
          freqDomainSolnVecReal, 
          freqDomainSolnVecImaginary
        );

      outputMgrAdapterRCPtr_->outputHB(
          timePoints,
          freqPoints,
          *timeDomainSolnVec,
          *freqDomainSolnVecReal, 
          *freqDomainSolnVecImaginary
        );

    }
  } 
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::printStepHeader()
// Purpose       : Prints out time step information.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void N_ANP_DCSweep::printStepHeader()
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

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, crStr);

  tmpStr = banner + "Start of DCOP STEP                        # ";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmpStr,
                           stepNumber+1, crStr);
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_ANP_DCSweep::printLoopInfo
// Purpose       :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
bool N_ANP_DCSweep::printLoopInfo(int start, int finish)
{
  bool bsuccess = N_ANP_AnalysisBase::printLoopInfo(start, finish);
  if (start == 0 && finish == 0)
  {
    outputFailureStats ();
  }
  return bsuccess;
}

