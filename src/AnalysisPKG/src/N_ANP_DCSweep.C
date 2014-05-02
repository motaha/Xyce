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
// Filename      : $RCSfile: N_ANP_DCSweep.C,v $
// Purpose       : DC Sweep class analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.51 $
// Revision Date  : $Date: 2014/02/24 23:49:12 $
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

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : DCSweep::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/22/10
//-----------------------------------------------------------------------------
bool DCSweep::setAnalysisParams(const N_UTL_OptionBlock & paramsBlock)
{
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "DCSweep::setDCAnalysisParams" << std::endl;
  }
#endif

  std::list<N_UTL_Param>::const_iterator it_tp;
  std::list<N_UTL_Param>::const_iterator it_param;
  std::list<N_UTL_Param>::const_iterator it_type;
  std::list<N_UTL_Param>::const_iterator first = paramsBlock.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last = paramsBlock.getParams().end();

  std::string msg;

  SweepParam sp;
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    for (it_tp = first; it_tp != last; ++it_tp)
    {
      Xyce::dout() << it_tp->uTag() ;
      Xyce::dout() << "\t";
      if (it_tp->uTag() == "PARAM" || it_tp->uTag() == "TYPE")
      {
        Xyce::dout() << it_tp->stringValue();
      }
      else
      {
        Xyce::dout() << it_tp->getImmutableValue<double>();
      }
      Xyce::dout() << std::endl;
    }
  }
#endif

  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag() == "TYPE")
    {
      it_type = it_tp;
      sp.type = it_tp->stringValue();
    }

    if (it_tp->uTag() == "PARAM")
    {
      it_param = it_tp;
      sp.name = it_tp->stringValue();
    }
  }

  it_tp = it_param;
  ++it_tp;
  if (sp.type == "LIN") // default
  {
    sp.startVal = it_tp->getImmutableValue<double>(); ++it_tp;
    sp.stopVal  = it_tp->getImmutableValue<double>(); ++it_tp;
    sp.stepVal  = it_tp->getImmutableValue<double>(); ++it_tp;
  }
  else if (sp.type == "DEC" || sp.type == "OCT")
  {
    sp.startVal = it_tp->getImmutableValue<double>(); ++it_tp;
    sp.stopVal  = it_tp->getImmutableValue<double>(); ++it_tp;
    sp.numSteps = it_tp->getImmutableValue<int>(); ++it_tp;
  }
  else if (sp.type == "LIST")
  {
    for (;it_tp!=last;++it_tp)
    {
      sp.valList.push_back(it_tp->getImmutableValue<double>());
    }
  }
  else
  {
    msg =  "DCSweep::setDCParams: ";
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
// Function      : DCSweep::outputFailureStats
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/21/2010
//-----------------------------------------------------------------------------
bool DCSweep::outputFailureStats(std::ostream &os)
{
  if (!(dcSweepFailures_.empty()))
  {
    os << "\tFailed DC sweep steps:\t\t" << std::endl;
    
    for (std::list<int>::iterator iter = dcSweepFailures_.begin(); iter != dcSweepFailures_.end(); ++iter)
    {
      os << "\t\tDC Step # " << *iter << std::endl;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::run()
// Purpose       : This is the main controlling loop for DC sweep analysis.
//                 This loop calls a series of operating point calculations
//                 (ie calculations in which there is no time integration,
//                 and the time integration method is always set to "none").
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::run()
{
  bool bsuccess = true;
  bsuccess = init();
  bsuccess &= loopProcess();
  bsuccess &= finish();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::init()
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
      Xyce::dout() << std::endl << std::endl
                   << Xyce::subsection_divider << std::endl
                   << "DCSweep::run()" << std::endl;
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
// Function      : DCSweep::initializeSolution()
// Purpose       : Move solution-initialization steps to separate routine so
//                 the process may be repeated cleanly when necessary.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystem Modeling
// Creation Date : 2/17/2010
//-----------------------------------------------------------------------------
void DCSweep::initializeSolution_()
{
  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loaderRCPtr_->setInitialGuess (anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr);

  // If available, set initial solution (.IC, .NODESET, etc).
  inputOPFlag_ =
    outputMgrAdapterRCPtr_->setupInitialConditions( *(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr),
                        *(anaManagerRCPtr_->getTIADataStore()->flagSolutionPtr));

  // Set a constant history for operating point calculation
  anaManagerRCPtr_->getTIADataStore()->setConstantHistory();
  anaManagerRCPtr_->getTIADataStore()->computeDividedDifferences();
  wimRCPtr_->obtainCorrectorDeriv();
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::loopProcess()
{
  bool bsuccess = true;
  // DC Sweep loop:

  int currentStep = 0;
  int finalStep = dcLoopSize_;
  while (currentStep < finalStep)
  {
    outputMgrAdapterRCPtr_->setDCAnalysisStepNumber( currentStep );

#ifdef Xyce_VERBOSE_TIME
    printStepHeader(Xyce::lout());
    secRCPtr_->outputTimeInfo(Xyce::lout());
#endif

    updateSweepParams_(currentStep, *dcParamVec_);
    if (currentStep != 0 &&(anaManagerRCPtr_->getSweepSourceResetFlag()))
    {
      anaManagerRCPtr_->getTIADataStore()->setZeroHistory();
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
// Function      : DCSweep::processSuccessfulStep()
//
// Purpose       : Used by both function dcSweepLoop and 2-level solve
//                 function calls to process successful DC steps.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::processSuccessfulStep()
{
  bool bsuccess = true;
  loaderRCPtr_->stepSuccess(anaManagerRCPtr_->currentMode_);

  // This output call is for device-specific output, such as .OP,
  // or internal plot output from PDE(TCAD) devices.
  loaderRCPtr_->output();

  if (sensFlag_ && !firstDoubleDCOPStep_() )
  {
    nlsMgrRCPtr_->calcSensitivity(objectiveVec_, 
        dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
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
  anaManagerRCPtr_->getTIADataStore()->updateSolDataArrays();

#ifdef Xyce_DEBUG_ANALYSIS
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 1)
    anaManagerRCPtr_->getTIADataStore()->outputSolDataArrays(Xyce::dout());
#endif
#endif

  dcSweepOutput();

  // now that output has been called, update the doubleDCOP step
  // if neccessary. (pde-only)
  doubleDCOPStep_ = tiaParams.lastDCOPStep;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::processFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::processFailedStep()
{
  loaderRCPtr_->stepFailure (anaManagerRCPtr_->currentMode_);

  stepNumber += 1;
  dcSweepFailures_.push_back(stepNumber);
  totalNumberFailedStepsAttempted_  += 1;
  secRCPtr_->numberSuccessiveFailures += 1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::finish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/28/08
//-----------------------------------------------------------------------------
bool DCSweep::finish()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  Xyce::dout() << "Calling DCSweep::finish() outputs!" << std::endl;
#endif

  outputMgrAdapterRCPtr_->finishOutput();
  if (!(dcSweepFailures_.empty()))
  {
    bsuccess = false;
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DCSweep::handlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool DCSweep::handlePredictor()
{
  anaManagerRCPtr_->getTIADataStore()->setErrorWtVector();
  wimRCPtr_->obtainPredictor();
  wimRCPtr_->obtainPredictorDeriv();

  // In case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  loaderRCPtr_->startTimeStep ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::takeStep_
// Purpose       : Take a DC Sweep integration step.
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void DCSweep::takeStep_()
{
  handlePredictor();
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();
  anaManagerRCPtr_->getTIADataStore()->stepLinearCombo ();
  gatherStepStatistics_ ();
  secRCPtr_->evaluateStepError ();

  return;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::twoLevelStep
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
bool DCSweep::twoLevelStep()
{
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();
  anaManagerRCPtr_->getTIADataStore()->stepLinearCombo ();
  gatherStepStatistics_ ();
  secRCPtr_->evaluateStepError ();

  return secRCPtr_->stepAttemptStatus;
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::dcSweepOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 1/31/08
//-----------------------------------------------------------------------------
void DCSweep::dcSweepOutput()
{
  // if this is an MPDE run, don't output anything.
  bool mpdeFlag = anaManagerRCPtr_->getMPDEFlag();
  bool hbFlag = anaManagerRCPtr_->getHBFlag();

  // Make sure this isn't the NLP step of a PDE DCOP.
  if ( !mpdeFlag  )
  {
    if (!firstDoubleDCOPStep_())
    {
      // conventional .PRINT output
      outputMgrAdapterRCPtr_->dcOutput(
          stepNumber,
          *(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr),
          *(anaManagerRCPtr_->getTIADataStore()->currStatePtr),
          *(anaManagerRCPtr_->getTIADataStore()->currStorePtr),
          objectiveVec_,
          dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

      // output for DCOP restart
      outputMgrAdapterRCPtr_->outputDCOP( *(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr) );
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
      Teuchos::RefCountPtr<N_LAS_BlockVector> timeDomainStoreVec;
      Teuchos::RefCountPtr<N_LAS_BlockVector> freqDomainStoreVecReal;
      Teuchos::RefCountPtr<N_LAS_BlockVector> freqDomainStoreVecImaginary;  


      Teuchos::RCP<const AnalysisBase> analysisObject
        = anaManagerRCPtr_->getAnalysisObject();

      Teuchos::rcp_dynamic_cast<const HB>( analysisObject )->prepareHBOutput(
          *(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr),
          timePoints,
          freqPoints,
          timeDomainSolnVec,
          freqDomainSolnVecReal,
          freqDomainSolnVecImaginary,
          timeDomainStoreVec,
          freqDomainStoreVecReal,
          freqDomainStoreVecImaginary
        );

      outputMgrAdapterRCPtr_->outputHB(
          timePoints,
          freqPoints,
          *timeDomainSolnVec,
          *freqDomainSolnVecReal,
          *freqDomainSolnVecImaginary,
          *timeDomainStoreVec,
          *freqDomainStoreVecReal,
          *freqDomainStoreVecImaginary
        );

    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::printStepHeader()
// Purpose       : Prints out time step information.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void DCSweep::printStepHeader(std::ostream &os)
{
#ifdef Xyce_VERBOSE_TIME
  lout() << std::endl << std::endl
         << "***** "<< (DEBUG_ANALYSIS ? commandLine_.getArgumentValue("netlist") : "") << "  Start of DCOP STEP                        # " << stepNumber+1
         << std::endl << std::endl;
#endif
}

//-----------------------------------------------------------------------------
// Function      : DCSweep::printLoopInfo
// Purpose       :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
bool DCSweep::printLoopInfo(int start, int finish)
{
  bool bsuccess = AnalysisBase::printLoopInfo(start, finish);
  if (start == 0 && finish == 0)
  {
    outputFailureStats(Xyce::lout());
  }
  return bsuccess;
}

} // namespace Analysis
} // namespace Xyce
