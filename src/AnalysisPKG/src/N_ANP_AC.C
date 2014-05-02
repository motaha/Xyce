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
//-------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_ANP_AC.C,v $
// Purpose       : AC analysis functions.
// Special Notes :
// Creator       : Ting Mei
// Creation Date :  7/11
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.55 $
// Revision Date  : $Date: 2014/02/24 23:49:12 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iomanip>


// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_TIA_Assembler.h>
#include <N_TIA_DataStore.h>
#include <N_LOA_Loader.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_MPDE_Manager.h>
#include <N_IO_OutputMgr.h>
#include <N_TIA_StepErrorControl.h>
#include <N_UTL_Timer.h>
#include <N_UTL_LogStream.h>

#include <N_IO_CmdParse.h>

#include<N_ANP_AC.h>

#include <N_UTL_Misc.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockSystemHelpers.h>

#include <N_PDS_ParMap.h>
#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_ParComm.h>
#include <mpi.h>
#else
#include <N_PDS_SerialComm.h>
#endif

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Epetra_LinearProblem.h>
#include <Amesos.h>

#include<N_UTL_ExpressionData.h>
#include<N_NLS_ReturnCodes.h>

namespace Xyce {
namespace Analysis {
    
//-----------------------------------------------------------------------------
// Function      : AC::AC( AnalysisManager * )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
AC::AC( AnalysisManager * anaManagerPtr ) :
  AnalysisBase(anaManagerPtr),
//  initialIntegrationMethod_(TIAMethod_BACKWARD_EULER),
  dcopFlag_(true),
  bVecRealPtr(0),
  bVecImagPtr(0),
  acLoopSize_(0),
  stepMult_(0.0),
  fstep_(0.0),
  currentFreq_(0.0)
{
  bVecRealPtr = lasSystemRCPtr_->builder().createVector();
  bVecImagPtr = lasSystemRCPtr_->builder().createVector();

  bVecRealPtr->putScalar(0.0);
  bVecImagPtr->putScalar(0.0);

}

//-----------------------------------------------------------------------------
// Function      : AC::~AC()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
AC::~AC()
{

  if (bVecRealPtr) { delete bVecRealPtr; bVecRealPtr=0; }

  if (bVecImagPtr) { delete bVecImagPtr; bVecImagPtr=0; }

}

//-----------------------------------------------------------------------------
// Function      : AC::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .AC statement.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 6/11
//-----------------------------------------------------------------------------
bool AC::setAnalysisParams(const N_UTL_OptionBlock & paramsBlock)
{
  std::list<N_UTL_Param>::const_iterator it_tp;
  std::list<N_UTL_Param>::const_iterator first = paramsBlock.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last = paramsBlock.getParams().end();
  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag()      == "TYPE")
    {
      tiaParams.type = it_tp->stringValue();
//      tiaParams.tStartGiven = true;
    }
    else if (it_tp->uTag() == "NP")
    {
      tiaParams.np = it_tp->getImmutableValue<double>();
    }
    else if (it_tp->uTag() == "FSTART")
    {
      tiaParams.fStart = it_tp->getImmutableValue<double>();
    }
    else if (it_tp->uTag() == "FSTOP")
    {
      tiaParams.fStop = it_tp->getImmutableValue<double>();
    }
  }

//  tiaParams.debugLevel = 1;
  if (DEBUG_ANALYSIS && tiaParams.debugLevel > 0)
  {
    dout() << section_divider << std::endl
           << "AC simulation parameters" 
           //<< Util::push << std::endl
           << std::endl
           << "number of points  = " << tiaParams.np << std::endl
           << "starting frequency = " << tiaParams.fStart << std::endl
           << "stop frequency = " << tiaParams.fStop 
           //<< Util::pop << std::endl;
           << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::run()
{
  bool bsuccess = true;

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
// Function      : AC::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::init()
{
  bool bsuccess = true;

  acLoopSize_ = setupSweepParam_();

  if(dcopFlag_)
  {
    anaManagerRCPtr_->currentMode_ = 0;
  }

  // Get set to do the operating point.
  integrationMethod_ = TIAMethod_NONE;
  wimRCPtr_->createTimeIntegMethod(integrationMethod_);

  stepNumber            = 0;
  doubleDCOPFlag_ = loaderRCPtr_->getDoubleDCOPFlag ();
  doubleDCOPStep_ = tiaParams.firstDCOPStep;

  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loaderRCPtr_->setInitialGuess (anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr);

  // If available, set initial solution (.IC, .NODESET, etc).
  inputOPFlag_ = outputMgrAdapterRCPtr_->setupInitialConditions( *(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr), *(anaManagerRCPtr_->getTIADataStore()->flagSolutionPtr));

  // Set a constant history for operating point calculation
  anaManagerRCPtr_->getTIADataStore()->setConstantHistory();
  anaManagerRCPtr_->getTIADataStore()->computeDividedDifferences();
  wimRCPtr_->obtainCorrectorDeriv();

  // solving for DC op
  handlePredictor();
  loaderRCPtr_->updateSources();
  secRCPtr_->newtonConvergenceStatus = nlsMgrRCPtr_->solve();
  anaManagerRCPtr_->getTIADataStore()->stepLinearCombo ();
  gatherStepStatistics_ ();
  secRCPtr_->evaluateStepError ();

  if ( secRCPtr_->newtonConvergenceStatus <= 0)
  {
    std::string msg=
        "Solving for DC operating point failed! Cannot continue AC analysis.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  // only output DC op if the .op was specified in the netlist
  // or if this AC analysis is called from a .step loop
  if ( anaManagerRCPtr_->getStepFlag() || anaManagerRCPtr_->getDotOpFlag() )
  {
    outputMgrAdapterRCPtr_->dcOutput(
        stepNumber, 
        *anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr, 
        *anaManagerRCPtr_->getTIADataStore()->nextStatePtr, 
        *anaManagerRCPtr_->getTIADataStore()->nextStorePtr,
        objectiveVec_, 
          dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
    outputMgrAdapterRCPtr_->finishOutput();
  }

  // This for saving the data from the DC op.  different from the above where we are
  // concerned with generating normal output.
  outputMgrAdapterRCPtr_->outputDCOP( *(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr) );
  loaderRCPtr_->loadBVectorsforAC (bVecRealPtr, bVecImagPtr);

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : AC::loopProcess()
// Purpose       : Conduct the stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::loopProcess()
{
  bool bsuccess = true;

  anaManagerRCPtr_->getTIADataStore()->daeQVectorPtr->putScalar(0.0);
  anaManagerRCPtr_->getTIADataStore()->daeFVectorPtr->putScalar(0.0);

  anaManagerRCPtr_->getTIADataStore()->dFdxdVpVectorPtr->putScalar(0.0);
  anaManagerRCPtr_->getTIADataStore()->dQdxdVpVectorPtr->putScalar(0.0);
  anaManagerRCPtr_->getTIADataStore()->dQdxMatrixPtr->put(0.0);
  anaManagerRCPtr_->getTIADataStore()->dFdxMatrixPtr->put(0.0);

//  integrationMethod_ = 6;

  loaderRCPtr_->updateState
                ((anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->currSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->lastSolutionPtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->currStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->lastStatePtr),
               (anaManagerRCPtr_->getTIADataStore()->nextStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->currStorePtr),
               (anaManagerRCPtr_->getTIADataStore()->lastStorePtr)
               );

  loaderRCPtr_->loadDAEVectors
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
               (anaManagerRCPtr_->getTIADataStore()->nextStoreLeadCurrQCompPtr),
               (anaManagerRCPtr_->getTIADataStore()->daeQVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->daeFVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->dFdxdVpVectorPtr),
               (anaManagerRCPtr_->getTIADataStore()->dQdxdVpVectorPtr) );

//  loaderRCPtr_->loadDAEMatrices(anaManagerRCPtr_->getTIADataStore()->currSolutionPtr,
//      anaManagerRCPtr_->getTIADataStore()->currStatePtr, anaManagerRCPtr_->getTIADataStore()->currStateDerivPtr,
//      anaManagerRCPtr_->getTIADataStore()->dQdxMatrixPtr,  anaManagerRCPtr_->getTIADataStore()->dFdxMatrixPtr);

  loaderRCPtr_->loadDAEMatrices(anaManagerRCPtr_->getTIADataStore()->nextSolutionPtr,
      anaManagerRCPtr_->getTIADataStore()->nextStatePtr, anaManagerRCPtr_->getTIADataStore()->nextStateDerivPtr,
      anaManagerRCPtr_->getTIADataStore()->nextStorePtr,
      anaManagerRCPtr_->getTIADataStore()->dQdxMatrixPtr,  anaManagerRCPtr_->getTIADataStore()->dFdxMatrixPtr);

  CPtr_ = rcp(anaManagerRCPtr_->getTIADataStore()->dQdxMatrixPtr, false);
  GPtr_ = rcp(anaManagerRCPtr_->getTIADataStore()->dFdxMatrixPtr, false);

  if (tiaParams.debugLevel > 0)
  {
    Xyce::dout() << "dQdxMatrixPtr:" << std::endl;
    anaManagerRCPtr_->getTIADataStore()->dQdxMatrixPtr->printPetraObject( Xyce::dout() );

    Xyce::dout() << "dFdxMatrixPtr:" << std::endl;
    anaManagerRCPtr_->getTIADataStore()->dFdxMatrixPtr->printPetraObject( Xyce::dout() );

    Xyce::dout() << std::endl;
  }  

  createLinearSystem_();

  int currentStep = 0;
  int finalStep = acLoopSize_;

  bool stepAttemptStatus;

  while (currentStep < finalStep)
  {
    updateCurrentFreq_(currentStep);

    updateLinearSystemFreq_();

    stepAttemptStatus = solveLinearSystem_();

    currentStep++;

    if (stepAttemptStatus)
    {
      processSuccessfulStep();
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      processFailedStep();
    }

  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : AC::createLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------
bool AC::createLinearSystem_()
{

  bool bsuccess = true;

  RCP<N_PDS_Manager> pdsMgrPtr_;
  pdsMgrPtr_ = rcp(lasSystemRCPtr_->getPDSManager(), false);

  RCP<N_PDS_ParMap> baseMap;
  RCP<Epetra_CrsGraph> BaseFullGraph_;

  baseMap = rcp(pdsMgrPtr_->getParallelMap( "SOLUTION" ), false);
  BaseFullGraph_ = rcp( pdsMgrPtr_->getMatrixGraph("JACOBIAN"), false );

  int numBlocks = 2;

  RCP<N_PDS_ParMap> blockMap = createBlockParMap(numBlocks, *baseMap);

  // Create a block vector
  BPtr_ = rcp ( new N_LAS_BlockVector ( numBlocks, blockMap, baseMap ) );

  // -----------------------------------------------------
  // Now test block graphs.
  // -----------------------------------------------------

  std::vector<std::vector<int> > blockPattern(2);
  blockPattern[0].resize(2);
  blockPattern[0][0] = 0; blockPattern[0][1] = 1;
  blockPattern[1].resize(2);
  blockPattern[1][0] = 0; blockPattern[1][1] = 1;

  int offset = generateOffset( *baseMap );

  RCP<Epetra_CrsGraph> blockGraph = createBlockGraph( offset, blockPattern, *blockMap, *BaseFullGraph_);

  ACMatrixPtr_ = rcp ( new N_LAS_BlockMatrix( numBlocks, offset, blockPattern, *blockGraph, *BaseFullGraph_) );

  // First diagonal block
  ACMatrixPtr_->put( 0.0 ); // Zero out whole matrix
  ACMatrixPtr_->block( 0, 0 ).add(*GPtr_);
  // Second diagonal block
  ACMatrixPtr_->block( 1, 1 ).add(*GPtr_);

  BPtr_->putScalar( 0.0 );
  BPtr_->block( 0 ).addVec( 1.0, *bVecRealPtr);
  BPtr_->block( 1 ).addVec( 1.0, *bVecImagPtr);

  Amesos amesosFactory;

  XPtr_ = rcp ( new N_LAS_BlockVector (numBlocks, blockMap, baseMap) );
  XPtr_->putScalar( 0.0 );

  blockProblem = rcp(new Epetra_LinearProblem(&ACMatrixPtr_->epetraObj(), &XPtr_->epetraObj(), &BPtr_->epetraObj() ) );

  blockSolver = rcp( amesosFactory.Create( "Klu", *blockProblem ) );

  // Need to reindex the linear system because of the noncontiguous block map.
  Teuchos::ParameterList params;
  params.set( "Reindex", true );
  blockSolver->SetParameters( params );

  // Call symbolic factorization without syncronizing the values, since they are not necessary here. 
  int linearStatus = blockSolver->SymbolicFactorization();

  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;
    bsuccess = false;
  }

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : AC::updateLinearSystemFreq_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------

bool AC::updateLinearSystemFreq_()
{
  double omega;

  omega =  2.0 * M_PI * currentFreq_;

  ACMatrixPtr_->block( 0, 1).put( 0.0);
  ACMatrixPtr_->block( 0, 1).add(*CPtr_);
  ACMatrixPtr_->block( 0, 1).scale(-omega);

  ACMatrixPtr_->block(1, 0).put( 0.0);
  ACMatrixPtr_->block(1, 0).add(*CPtr_);
  ACMatrixPtr_->block(1, 0).scale(omega);

  // Copy the values loaded into the blocks into the global matrix for the solve.
  ACMatrixPtr_->assembleGlobalMatrix();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::solveLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :  Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/2011
//-----------------------------------------------------------------------------

bool AC::solveLinearSystem_()
{

  bool bsuccess = true;
 // Solve the block problem

  int linearStatus = blockSolver->NumericFactorization();

  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
    bsuccess = false;
  }

  linearStatus = blockSolver->Solve();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos solve exited with error: " << linearStatus;
    bsuccess = false;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Transient::resetForStepAnalysis()
// Purpose       : When doing a .STEP sweep, some data must be reset to its
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : T. Mei
// Creation Date : 04/16/12
//-----------------------------------------------------------------------------
bool AC::resetForStepAnalysis()
{
  totalNumberSuccessStepsThisParameter_ = 0;
  stepNumber = 0;
  beginningIntegration = true;

//  nlsMgrRCPtr_->resetAll(DC_OP);
//  secRCPtr_->resetAll();

  bVecRealPtr->putScalar(0.0);
  bVecImagPtr->putScalar(0.0);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::processSuccessfulStep()
{
  bool bsuccess = true;

// Output x.
  outputMgrAdapterRCPtr_->outputAC (currentFreq_, XPtr_->block(0), XPtr_-> block(1));

  if ( !firstDoubleDCOPStep_() )
  {
      stepNumber += 1;
      totalNumberSuccessStepsThisParameter_ += 1;
      totalNumberSuccessfulStepsTaken_ += 1;
  }

  // This output call is for device-specific output (like from a PDE device,
  // outputting mesh-based tecplot files).  It will only work in parallel if on
  // a machine where all processors have I/O capability.
  loaderRCPtr_->output();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : AC::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::processFailedStep()
{
  bool bsuccess = true;

  stepNumber += 1;
  acSweepFailures_.push_back(stepNumber);
  totalNumberFailedStepsAttempted_  += 1;
  secRCPtr_->numberSuccessiveFailures += 1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AC::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::finish()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  Xyce::dout() << "Calling AC::finish() outputs!" << std::endl;
#endif
  outputMgrAdapterRCPtr_->finishOutput();
  if (!(acSweepFailures_.empty()))
  {
    bsuccess = false;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : AC::handlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool AC::handlePredictor()
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
// Function      : AC::updateCurrentFreq_(int stepNumber )
// Purpose       :
// Special Notes : Used for AC analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::updateCurrentFreq_(int stepNumber)
{
  double fstart;
  fstart = tiaParams.fStart;

//  tiaParams.debugLevel = 1;

    if (tiaParams.type=="LIN")
    {

      currentFreq_  = fstart + static_cast<double>(stepNumber)*fstep_;
    }
    else if(tiaParams.type=="DEC" || tiaParams.type=="OCT")
    {

      currentFreq_ = fstart*pow(stepMult_, static_cast<double>(stepNumber) );
    }
    else
    {
      std::string msg=
        "AC::updateCurrentFreq_: unsupported STEP type";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);

    }

    if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
    {
      Xyce::dout() << "currentFreq_  = " << currentFreq_ << std::endl;
    }   

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::setupSweepParam_
// Purpose       : Processes sweep parameters.
// Special Notes : Used for AC analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
int AC::setupSweepParam_()
{
  double fstart, fstop;
  double fcount = 0.0;

  fstart = tiaParams.fStart;
  fstop = tiaParams.fStop;

#ifdef Xyce_DEBUG_ANALYSIS
  if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "AC::setupSweepParam_" << std::endl;
  }
#endif

    if (tiaParams.type=="LIN")
    {
      int np = static_cast<int> (tiaParams.np);

      if ( np == 1)
        fstep_ = 0;
      else
        fstep_  = (fstop - fstart)/(tiaParams.np - 1);

      fcount = tiaParams.np;
#ifdef Xyce_DEBUG_ANALYSIS
      if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
      {
        Xyce::dout() << "fstep   = " << fstep_  << std::endl;
      }
#endif
    }
    else if(tiaParams.type=="DEC")
    {
      stepMult_ = pow(10,(1/tiaParams.np));
      fcount   = floor(fabs(log10(fstart) - log10(fstop)) * tiaParams.np + 1);
#ifdef Xyce_DEBUG_ANALYSIS
      if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
#endif
    }
    else if(tiaParams.type=="OCT")
    {
      stepMult_ = pow(2,1/(tiaParams.np));

      // changed to remove dependence on "log2" function, which apparently
      // doesn't exist in the math libraries of FreeBSD or the mingw
      // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
      double ln2=log(2.0);
      fcount   = floor(fabs(log(fstart) - log(fstop))/ln2 * tiaParams.np + 1);
#ifdef Xyce_DEBUG_ANALYSIS
      if (anaManagerRCPtr_->tiaParams.debugLevel > 0)
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
#endif
    }
    else
    {
      std::string msg=
        "AC::setupSweepParam: unsupported type";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);

    }

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return static_cast<int> (fcount);
}

} // namespace Analysis
} // namespace Xyce
