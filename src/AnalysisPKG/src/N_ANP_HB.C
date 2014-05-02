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
// Filename      : $RCSfile: N_ANP_HB.C,v $
// Purpose       : HB analysis functions.
// Special Notes :
// Creator       : Todd Coffey, 1414, Ting Mei, 1437
// Creation Date : 07/23/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.60 $
// Revision Date  : $Date: 2014/02/24 23:49:12 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_HB.h>
#include <N_ANP_Transient.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_Report.h>
#include <N_MPDE_Manager.h>
#include <N_MPDE_Discretization.h>

#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_DataStore.h>

#include <N_LOA_HBLoader.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_HBPrecondFactory.h>
#include <N_LAS_PrecondFactory.h>
#include <N_LAS_System.h>
#include <N_LAS_BlockSystemHelpers.h>

#include <N_NLS_Manager.h>

#include <N_IO_OutputMgr.h>

#include <N_UTL_FFTInterface.hpp>

#include <Teuchos_Utils.hpp>

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : HB::HB( AnalysisManager * )
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
HB::HB( AnalysisManager * anaManagerPtr ) :
  AnalysisBase(anaManagerPtr),
  debugLevel(0),
  isPaused(false),
  startDCOPtime(0.0),
  endTRANtime(0.0),
  isTransient_(false),
  isDCSweep_(false),
  test_(false),
  size_(21),
  period_(1.0),
  startUpPeriods_(0),
  startUpPeriodsGiven_(false),
  startUpPeriodsFinished_(false),
  saveIcData_(false),
  tiaParams_( anaManagerPtr->tiaParams ),
  voltLimFlag_(1), 
  taHB_(1),
  fastTimeDisc_(0),
  fastTimeDiscOrder_(1),
  hbTotalNumberSuccessfulStepsTaken_(0),
  hbTotalNumberFailedStepsAttempted_(0),
  hbTotalNumberJacobiansEvaluated_(0),
  hbTotalNumberIterationMatrixFactorizations_(0),
  hbTotalNumberLinearSolves_(0),
  hbTotalNumberFailedLinearSolves_(0),
  hbTotalNumberLinearIters_(0),
  hbTotalNumberResidualEvaluations_(0),
  hbTotalNonlinearConvergenceFailures_(0),
  hbTotalResidualLoadTime_(0.0),
  hbTotalJacobianLoadTime_(0.0),
  hbTotalLinearSolutionTime_(0.0),
  resetForStepCalledBefore_(false)
{
  devInterfacePtr_ = anaManagerRCPtr_->devInterfacePtr;
  topoMgrPtr_ = anaManagerRCPtr_->topoMgrPtr;
  nonlinearEquationLoaderPtr_ = anaManagerRCPtr_->nonlinearEquationLoaderPtr;
  appBuilderPtr_ = anaManagerRCPtr_->appBuilderPtr;
  pdsMgrPtr_ = anaManagerRCPtr_->pdsMgrPtr;
}

//-----------------------------------------------------------------------------
// Function      : HB::getStepNumber()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
int HB::getStepNumber ()
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    return analysisObject_->getStepNumber();
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : HB::setStepNumber()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
void HB::setStepNumber (int step)
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    analysisObject_->setStepNumber( step );
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::setBeginningIntegrationFlag()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
void HB::setBeginningIntegrationFlag(bool bif)
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    analysisObject_->setBeginningIntegrationFlag( bif );
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::getBeginningIntegrationFlag()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
bool HB::getBeginningIntegrationFlag()
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    return analysisObject_->getBeginningIntegrationFlag();
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
void HB::setIntegrationMethod(int im)
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    analysisObject_->setIntegrationMethod( im );
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::getIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
unsigned int HB::getIntegrationMethod ()
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    return analysisObject_->getIntegrationMethod ();
  }
  return TIAMethod_NONE;
}

//-----------------------------------------------------------------------------
// Function      : HB::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::run()
{

  // initializeAll_();
  //   get TIAParams from AnalysisInterface
  //   create HBBuilder
  //     generateMaps
  //     generateStateMaps
  //   create vectors for IC
  //
  // if (test_) {
  //   runTests_();
  // } else {
  //
  //   computeInitialCondition_();
  //     run startup periods
  //     integrate one period
  //     interpolate to evenly spaced points
  //     set this as IC for HB nonlinear problem
  //
  //   setupHBProblem_();
  //
  //   runHBProblem_();
  //
  // }

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
// Function      : HB::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::init()
{
  bool returnValue=true;

  Xyce::lout() << " ***** Running HB initial conditions....\n" << std::endl;
#ifdef Xyce_DEBUG_HB
  Xyce::dout() << std::endl
         << section_divider << std::endl
         << "  HB::init()" << std::endl;
#endif // Xyce_DEBUG_HB

  //Store copy of transient TIAParams for HB run
  tiaParams_ = anaManagerRCPtr_->tiaParams;
  //
  // If it was requested, advance the solution a fixed number of startup periods
  //
  period_ =   1.0/tiaParams_.freqs[0];


  if (taHB_ == 1)
  {
  bool retTol1 = runTol_(); returnValue = returnValue && retTol1;

  // Start up periods need to be run before the initial condition is computed, otherwise
  // just used the solution from the tolerance calculation.
  if( startUpPeriodsGiven_ )
  {
    bool startupPeriodsSuccess = runStartupPeriods_();
    if (!startupPeriodsSuccess)
    {
      std::string msg = "HB::init().  Failed to calculate the startup periods.\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
    }
    returnValue = returnValue && startupPeriodsSuccess;

    bool icSuccess = runTransientIC_();
    if (!icSuccess)
    {
      std::string msg = "HB::init().  Initial HB Transient failed \n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
    }
    returnValue = returnValue && icSuccess;
  }

  interpolateIC_();


  }
  else
  {

    double TimeStep = period_/static_cast<double>(size_);
    timeSteps_.push_back( TimeStep );
    fastTimes_.resize(size_+1);

    goodTimePoints_.resize(size_+1);
    for( int i = 0; i <= size_; ++i )
    {
      fastTimes_[i] = tiaParams_.initialTime + static_cast<double>(i) * TimeStep;
    }
    goodTimePoints_ = fastTimes_;

  }

  freqPoints_.resize(size_);
  //   double fastStep = period_/static_cast<double>(size_);
  int i=0;
  for( i = 0; i < size_; ++i )
  {
    if (i < (size_-1)/2)
      freqPoints_[i] = -static_cast<double>((size_-1)/2 - i) * tiaParams_.freqs[0];
    else
      freqPoints_[i] = static_cast<double>(i - (size_-1)/2) * tiaParams_.freqs[0];
  }


  // now that we have size_, continue with the initialization of objects for HB
  mpdeDiscPtr_ = rcp(new N_MPDE_Discretization(
     static_cast<N_MPDE_Discretization::Type>(fastTimeDisc_), fastTimeDiscOrder_ ));

  hbBuilderPtr_ = rcp(new N_LAS_HBBuilder( size_, mpdeDiscPtr_ ));

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB::init():  Generate Maps\n";
#endif // Xyce_DEBUG_HB
  hbBuilderPtr_->generateMaps( rcp(pdsMgrPtr_->getParallelMap( "SOLUTION" ), false),
                               rcp(pdsMgrPtr_->getParallelMap( "SOLUTION_OVERLAP_GND" ), false) );
  hbBuilderPtr_->generateStateMaps( rcp(pdsMgrPtr_->getParallelMap( "STATE" ),false) );
  hbBuilderPtr_->generateStoreMaps( rcp(pdsMgrPtr_->getParallelMap( "STORE" ),false) );
  hbBuilderPtr_->generateGraphs( *(pdsMgrPtr_->getMatrixGraph( "JACOBIAN" )) );

  HBICVectorPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();
  HBICStateVectorPtr_ = hbBuilderPtr_->createTimeDomainStateBlockVector();
  HBICStoreVectorPtr_ = hbBuilderPtr_->createTimeDomainStoreBlockVector();
  HBICQVectorPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();

  HBICVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();

//  HBICStateVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeStateBlockVector();
//  HBICQVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();
//  HBICStoreVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeStoreBlockVector();

  tiaParams_.finalTime =  period_;
  anaManagerRCPtr_->registerTIAParams (tiaParams_);

  // set the fast source flag in the devices
  std::vector<std::string> srcVec;
  devInterfacePtr_->registerFastSources( srcVec );

  //anaManagerRCPtr_->anaIntPtr->getnlHBOptions(saved_nlHBOB_);
  nlsMgrRCPtr_->getHBOptions(saved_nlHBOB_);

  // Create HB Loader.
  hbLoaderPtr_ = rcp(new N_LOA_HBLoader( mpdeState_, mpdeDiscPtr_ ) );
  hbLoaderPtr_->registerHBBuilder(hbBuilderPtr_);
  hbLoaderPtr_->registerAppBuilder(appBuilderPtr_);

  if (taHB_==1)
  {
  // Pick up IC data from the initial transient.
  for (int i=0 ; i<size_ ; ++i)
  {
#ifdef Xyce_DEBUG_HB
    if( debugLevel > 0 )
    {
      Xyce::dout() << "HB::init():  Loading initial condition data from time: fastTimes_["
                << i << "] = " << fastTimes_[i] << std::endl;
    }
#endif // Xyce_DEBUG_HB

    HBICVectorPtr_->block(i) = *(goodSolutionVec_[i]);
    HBICStateVectorPtr_->block(i) = *(goodStateVec_[i]);
    HBICQVectorPtr_->block(i) = *(goodQVec_[i]);
    HBICStoreVectorPtr_->block(i) = *(goodStoreVec_[i]);
  }

  hbLoaderPtr_->permutedFFT(*HBICVectorPtr_, &*HBICVectorFreqPtr_);

  }
#ifdef Xyce_DEBUG_HB
  if ( debugLevel > 1 )
  {
    Xyce::dout() << "HB Initial Condition Solution!\n";
    HBICVectorPtr_->printPetraObject(std::cout);
    Xyce::dout() << "HB Initial Condition State Vector!\n";
    HBICStateVectorPtr_->printPetraObject(std::cout);
    Xyce::dout() << "HB Initial Condition Store Vector!\n";
    HBICStoreVectorPtr_->printPetraObject(std::cout);
  }
#endif // Xyce_DEBUG_HB

  //Destroy Solvers, etc. from IC phase and prepare for HB
  //-----------------------------------------

//  if ( voltLimFlag_ == 0 )
//    devInterfacePtr_->setVoltageLimiterFlag (false);
  devInterfacePtr_->setVoltageLimiterFlag (voltLimFlag_);

  devInterfacePtr_->setMPDEFlag( true );
  anaManagerRCPtr_->resetAll();

  //-----------------------------------------

  //Finish setup of HB Loader
  //-----------------------------------------
  hbLoaderPtr_->registerAppLoader( loaderRCPtr_ );
  hbLoaderPtr_->registerDeviceInterface( devInterfacePtr_ );
  goodTimePoints_.resize(size_+1);
  goodTimePoints_[size_] = period_;

  for( int i = 0; i < size_; ++i )
  {
    goodTimePoints_[i] = goodTimePoints_[i] - tiaParams_.initialTime;
  }

  hbLoaderPtr_->setFastTimes(goodTimePoints_);

  //-----------------------------------------
  //Construct Solvers, etc. for HB Phase
  //-----------------------------------------
  lasHBSysPtr_ = rcp(new N_LAS_System());
  //-----------------------------------------

  anaManagerRCPtr_->registerTIAParams( tiaParams_ );
  anaManagerRCPtr_->registerLinearSystem( &*lasHBSysPtr_ );
  anaManagerRCPtr_->registerLoader( &*hbLoaderPtr_ );

  //hack needed by TIA initialization currently
  hbBuilderPtr_->registerPDSManager( &*pdsMgrPtr_ );

  lasHBSysPtr_->registerANPInterface( (anaManagerRCPtr_->anaIntPtr).get() );
  lasHBSysPtr_->registerPDSManager( &*pdsMgrPtr_ );
  lasHBSysPtr_->registerBuilder( &*hbBuilderPtr_ );

  //need to cut out unnecessary stuff from this call for new dae
  lasHBSysPtr_->initializeSystem();

  // Give NLS Manager the same old nonlinearEquationLoader as it just calls the TIA loader in newDAE
  nlsMgrRCPtr_->registerLinearSystem( &*lasHBSysPtr_ );
  nlsMgrRCPtr_->setLinSolOptions( saved_lsHBOB_ );
  nlsMgrRCPtr_->setMatrixFreeFlag( true );

  // Let the HB loader know that the application of the operator is matrix free
  hbLoaderPtr_->setMatrixFreeFlag( true );

  if (Teuchos::is_null( precFactory_ ))
  {
    // Generate the HB preconditioner factory.
    precFactory_ = rcp( new N_LAS_HBPrecondFactory( saved_lsHBOB_ ) );
  }

  // Register application loader with preconditioner factory
  RCP<N_LAS_HBPrecondFactory> tmpPrecFactory
    = rcp_dynamic_cast<N_LAS_HBPrecondFactory>( precFactory_ );

  tmpPrecFactory->registerAppBuilder( appBuilderPtr_ );
  tmpPrecFactory->registerHBLoader( hbLoaderPtr_ );
  tmpPrecFactory->registerHBBuilder( hbBuilderPtr_ );
  tmpPrecFactory->setFastTimes( goodTimePoints_ );
  tmpPrecFactory->setTimeSteps( timeSteps_ );

  nlsMgrRCPtr_->registerPrecondFactory( precFactory_ );
  //-----------------------------------------

  //Initialization of Solvers, etc. for HB Phase
  //-----------------------------------------
  //Dummy call to setup time integrator for transient
  anaManagerRCPtr_->resumeSimulation();
  anaManagerRCPtr_->initializeAll();
  anaManagerRCPtr_->setMPDEFlag( true );

  nlsMgrRCPtr_->initializeAll();
  nlsMgrRCPtr_->setAnalysisMode(anpAnalysisModeToNLS(ANP_MODE_HB));

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << section_divider << std::endl;
#endif // Xyce_DEBUG_HB

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::loopProcess()
{
  bool returnValue = true;

  Xyce::lout() << " ***** Beginning full HB simulation....\n" << std::endl;

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << std::endl
         << section_divider << std::endl
         << "  HB::loopProcess" << std::endl;
#endif // Xyce_DEBUG_HB

  RCP<N_TIA_DataStore> dsPtr = anaManagerRCPtr_->getTIADataStore();
  *(dsPtr->nextSolutionPtr) = *(HBICVectorFreqPtr_.get());
  *(dsPtr->nextStatePtr) = *(HBICStateVectorPtr_.get());
  *(dsPtr->nextStorePtr) = *(HBICStoreVectorPtr_.get());

  // try to run the problem
  analysisObject_ = Teuchos::rcp(new DCSweep( anaManagerRCPtr_.get() ));
  returnValue = analysisObject_->run();

  // Add in simulation times
  accumulateStatistics_();

  // print out analysis info
  Xyce::lout() << " ***** Harmonic Balance Computation Summary *****" << std::endl;
  analysisObject_->printLoopInfo( 0, 0 );

#ifdef Xyce_DEBUG_HB
  dout() << section_divider << std::endl;
#endif // Xyce_DEBUG_HB

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::processSuccessfulDCOP()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processSuccessfulDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processSuccessfulStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processFailedStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processFailedDCOP
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processFailedDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::finish()
{
  // Move statistics into common variable
  totalNumberSuccessfulStepsTaken_ = hbTotalNumberSuccessfulStepsTaken_;
  totalNumberFailedStepsAttempted_ = hbTotalNumberFailedStepsAttempted_;
  totalNumberJacobiansEvaluated_ = hbTotalNumberJacobiansEvaluated_;
  totalNumberIterationMatrixFactorizations_ = hbTotalNumberIterationMatrixFactorizations_;
  totalNumberLinearSolves_ = hbTotalNumberLinearSolves_;
  totalNumberFailedLinearSolves_ = hbTotalNumberFailedLinearSolves_;
  totalNumberLinearIters_ = hbTotalNumberLinearIters_;
  totalNumberResidualEvaluations_ = hbTotalNumberResidualEvaluations_;
  totalNonlinearConvergenceFailures_ = hbTotalNonlinearConvergenceFailures_;
  totalResidualLoadTime_ = hbTotalResidualLoadTime_;
  totalJacobianLoadTime_ = hbTotalJacobianLoadTime_;
  totalLinearSolutionTime_ = hbTotalLinearSolutionTime_;

  return true;
}

bool HB::handlePredictor()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::resetForStepAnalysis()
// Purpose       : When doing a .STEP sweep, some data must be reset to its
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::resetForStepAnalysis()
{
  totalNumberSuccessStepsThisParameter_ = 0;

  if (resetForStepCalledBefore_)
  {
    goodSolutionVec_.clear();
    goodStateVec_.clear();
    goodQVec_.clear();
    goodStoreVec_.clear();

    secRCPtr_->resetAll();

    anaManagerRCPtr_->setNextOutputTime(0.0);
    anaManagerRCPtr_->registerLinearSystem( &*lasSystemRCPtr_ );
    anaManagerRCPtr_->registerLoader( loaderRCPtr_.get() );
    anaManagerRCPtr_->resumeSimulation();

    nlsMgrRCPtr_->resetAll(DC_OP);
    nlsMgrRCPtr_->setMatrixFreeFlag( false );
    nlsMgrRCPtr_->registerLinearSystem( &*lasSystemRCPtr_ );
    nlsMgrRCPtr_->registerLoader( nonlinearEquationLoaderPtr_.get() );
    nlsMgrRCPtr_->setLinSolOptions( saved_lsOB_ );

    devInterfacePtr_->registerLinearSystem( &*lasSystemRCPtr_ );

    // un-set the fast source flag in the devices
    std::vector<std::string> srcVec;
    devInterfacePtr_->deRegisterFastSources( srcVec );

    anaManagerRCPtr_->initializeAll();

    devInterfacePtr_->resetForStepAnalysis();
    devInterfacePtr_->initializeAll();
    devInterfacePtr_->setMPDEFlag( false );

    nlsMgrRCPtr_->initializeAll();

    anaManagerRCPtr_->unset_resumeSimulation();

    // this would normally be allocated in the analysis manager::run function.
    anaManagerRCPtr_->xyceTranTimerPtr_ = rcp(new N_UTL_Timer(*(pdsMgrPtr_->getPDSComm())));
  }

  resetForStepCalledBefore_=true;

  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::finalVerboseOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::finalVerboseOutput()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::setHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 05/13/13
//-----------------------------------------------------------------------------
bool HB::setHBOptions(const N_UTL_OptionBlock & OB)
{
  std::list<N_UTL_Param>::const_iterator iterPL = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator  endPL = OB.getParams().end();

  for( ; iterPL != endPL; ++iterPL )
  {
    ExtendedString tag = iterPL->tag();
    tag.toUpper();

    if ( tag == "NUMFREQ" )
    {
      size_    = iterPL->getImmutableValue<int>();
    }
    else if ( tag == "STARTUPPERIODS" )
    {
      startUpPeriods_ = iterPL->getImmutableValue<int>();

      if (startUpPeriods_ > 0)
        startUpPeriodsGiven_ = true;
    }
    else if( tag == "SAVEICDATA" )
    {
      saveIcData_ = true;
    }
    else if( tag == "TEST" )
    {
      test_     = static_cast<bool> (iterPL->getImmutableValue<int>());
    }
    else if (tag == "DEBUGLEVEL" )
    {
      debugLevel = iterPL->getImmutableValue<int>();
    }
    else if ( tag == "TAHB" )
    {
      taHB_ = iterPL->getImmutableValue<int>();
    }
    else if ( tag == "VOLTLIM" )
    {
      voltLimFlag_ = static_cast<bool> (iterPL->getImmutableValue<int>());
    }  
    else
    {
      UserWarning(*this) << "Unrecognized HBINT option " << tag;
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : HB::setLinSol
// Purpose       : this is needed for .STEP to work with HB
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/12/2013
//-----------------------------------------------------------------------------
bool HB::setLinSol(const N_UTL_OptionBlock & OB)
{
  // Save the non-HB linear solver option block
  saved_lsOB_ = OB;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setHBLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/13/13
//-----------------------------------------------------------------------------
bool HB::setHBLinSol(const N_UTL_OptionBlock & OB)
{
  // Save the HB linear solver option block
  saved_lsHBOB_ = OB;

  // Generate the HB preconditioner factory.
  precFactory_ = rcp( new N_LAS_HBPrecondFactory( OB ) );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::isAnalysis
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/13/13
// Notes         : Alternatively, we could try to cast the analysis object
//               : However, this method is called a lot.
//-----------------------------------------------------------------------------
bool HB::isAnalysis( int analysis_type )
{
  bool returnValue = false;

  if ( analysis_type == ANP_MODE_TRANSIENT )
  {
    returnValue = isTransient_;
  }
  if ( analysis_type == ANP_MODE_DC_SWEEP )
  {
    returnValue = isDCSweep_;
  }
  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::prepareHBOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 08/20/07
//-----------------------------------------------------------------------------
void HB::prepareHBOutput(
    N_LAS_Vector & solnVecPtr,
    std::vector<double> & timePoints,
    std::vector<double> & freqPoints,
    RCP<N_LAS_BlockVector> & timeDomainSolnVec,
    RCP<N_LAS_BlockVector> & freqDomainSolnVecReal,
    RCP<N_LAS_BlockVector> & freqDomainSolnVecImag,
    RCP<N_LAS_BlockVector> & timeDomainStoreVec,
    RCP<N_LAS_BlockVector> & freqDomainStoreVecReal,
    RCP<N_LAS_BlockVector> & freqDomainStoreVecImag
    ) const
{
  N_LAS_BlockVector & blockSolVecPtr = dynamic_cast<N_LAS_BlockVector &>(solnVecPtr);

  Teuchos::RCP<N_LAS_BlockVector> bStoreVecFreqPtr_ = hbLoaderPtr_->getStoreVecFreqPtr();

  timeDomainStoreVec = hbBuilderPtr_->createTimeDomainStoreBlockVector();

  if (bStoreVecFreqPtr_->blockCount() > 0 )
  {
    hbLoaderPtr_->permutedIFT(*bStoreVecFreqPtr_, &*timeDomainStoreVec);
//  bStoreVecFreqPtr_->printPetraObject(std::cout);
  } 

  //TD solution
  timeDomainSolnVec = hbBuilderPtr_->createTimeDomainBlockVector(); 
  int blockCount = timeDomainSolnVec->blockCount();
  int N = timeDomainSolnVec->block(0).globalLength(); 

  timePoints.resize(size_);

  for( int i = 0; i < size_; ++i )
  {
    timePoints[i] = fastTimes_[i] - tiaParams_.initialTime;
  }

  freqPoints = freqPoints_;

  // Create individual block vectors to store the real and imaginary parts separately.
  Teuchos::RCP<N_PDS_ParMap> baseMap = Teuchos::rcp_const_cast<N_PDS_ParMap>( hbBuilderPtr_->getBaseSolutionMap() );
  Teuchos::RCP<N_PDS_ParMap> globalMap = createBlockParMap( blockCount, *baseMap );
  freqDomainSolnVecReal = Teuchos::rcp( new N_LAS_BlockVector( blockCount, globalMap, baseMap ) );
  freqDomainSolnVecImag = Teuchos::rcp( new N_LAS_BlockVector( blockCount, globalMap, baseMap ) );

  hbLoaderPtr_->permutedIFT(blockSolVecPtr, &*timeDomainSolnVec); 

//  Xyce::dout() << "HB X Vector TD" << std::endl;
//  timeDomainSolnVec->printPetraObject(std::cout);

//  Xyce::dout() << "HB Store Vector TD" << std::endl;
//  timeDomainStoreVec->printPetraObject(std::cout);

  // Now copy over the frequency domain solution, real and imaginary parts separately, into the output vectors.
  for (int j=0; j<N; j++)
  {
    // See if this time-domain solution variable is owned by the local processor.
    // If so, this processor owns the entire j-th block of the blockSolVecPtr vector,
    // and the j-th entry of every block in the freqDomainSolnVec[Real/Imag] vector.
    int lid = baseMap->globalToLocalIndex( j );
    N_LAS_Vector& solBlock = blockSolVecPtr.block( j );

    N_LAS_Vector& realVecRef =  freqDomainSolnVecReal->block((blockCount-1)/2);
    N_LAS_Vector& imagVecRef =  freqDomainSolnVecImag->block((blockCount-1)/2);

    if (lid >= 0)
    { 
      realVecRef[lid] = solBlock[0];  
      imagVecRef[lid] = solBlock[1];
    }

    for (int i=1; i <= (blockCount-1)/2; ++i)
    {
      N_LAS_Vector& realVecRef_neg =  freqDomainSolnVecReal->block((blockCount-1)/2 - i);
      N_LAS_Vector& imagVecRef_neg =  freqDomainSolnVecImag->block((blockCount-1)/2 - i);
      N_LAS_Vector& realVecRef_pos =  freqDomainSolnVecReal->block((blockCount-1)/2 + i);
      N_LAS_Vector& imagVecRef_pos =  freqDomainSolnVecImag->block((blockCount-1)/2 + i);

      if (lid >= 0)
      {
        realVecRef_neg[lid] = solBlock[ 2*(blockCount-i) ];
        imagVecRef_neg[lid] = solBlock[ 2*(blockCount-i) + 1 ];
        realVecRef_pos[lid] = solBlock[ 2*i ];
        imagVecRef_pos[lid] = solBlock[ 2*i+1 ];
      }
    }
  } 

  // proceed to store variables
  Teuchos::RCP<N_PDS_ParMap> baseStoreMap = Teuchos::rcp_const_cast<N_PDS_ParMap>( hbBuilderPtr_->getBaseStoreMap() );
  Teuchos::RCP<N_PDS_ParMap> globalStoreMap = createBlockParMap( blockCount, *baseStoreMap );
  freqDomainStoreVecReal = Teuchos::rcp( new N_LAS_BlockVector( blockCount, globalStoreMap, baseStoreMap ) );
  freqDomainStoreVecImag = Teuchos::rcp( new N_LAS_BlockVector( blockCount, globalStoreMap, baseStoreMap ) );

//  hbLoaderPtr_->permutedIFT(blockSolVecPtr, &*timeDomainSolnVec);

  N = timeDomainStoreVec->block(0).globalLength(); 

  for (int j=0; j<N; j++)
  {
    // See if this time-domain solution variable is owned by the local processor.
    // If so, this processor owns the entire j-th block of the blockSolVecPtr vector,
    // and the j-th entry of every block in the freqDomainSolnVec[Real/Imag] vector.
    int lid = baseStoreMap->globalToLocalIndex( j );
    N_LAS_Vector& storeBlock =  bStoreVecFreqPtr_->block( j );

    N_LAS_Vector& realVecRef =  freqDomainStoreVecReal->block((blockCount-1)/2);
    N_LAS_Vector& imagVecRef =  freqDomainStoreVecImag->block((blockCount-1)/2);

    if (lid >= 0)
    { 
      realVecRef[lid] = storeBlock[0];  
      imagVecRef[lid] = storeBlock[1];
    }

    for (int i=1; i <= (blockCount-1)/2; ++i)
    {
      N_LAS_Vector& realVecRef_neg =  freqDomainStoreVecReal->block((blockCount-1)/2 - i);
      N_LAS_Vector& imagVecRef_neg =  freqDomainStoreVecImag->block((blockCount-1)/2 - i);
      N_LAS_Vector& realVecRef_pos =  freqDomainStoreVecReal->block((blockCount-1)/2 + i);
      N_LAS_Vector& imagVecRef_pos =  freqDomainStoreVecImag->block((blockCount-1)/2 + i);

      if (lid >= 0)
      {
        realVecRef_neg[lid] = storeBlock[ 2*(blockCount-i) ];
        imagVecRef_neg[lid] = storeBlock[ 2*(blockCount-i) + 1 ];
        realVecRef_pos[lid] = storeBlock[ 2*i ];
        imagVecRef_pos[lid] = storeBlock[ 2*i+1 ];
      }
    }
  } 

//  Xyce::dout() << "HB X Vector FD" << std::endl;
//  freqDomainSolnVecReal->printPetraObject(std::cout);
//  freqDomainSolnVecImag->printPetraObject(std::cout);

//  Xyce::dout() << "HB Store Vector FD" << std::endl;

//  freqDomainStoreVecReal->printPetraObject(std::cout);
//  freqDomainStoreVecImag->printPetraObject(std::cout); 

}


//-----------------------------------------------------------------------------
// Function      : HB::accumulateStatistics()
// Purpose       : Add in the statistics from the current analysis object
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, 1355, Electrical Models & Simulation
// Creation Date : 05/29/13
//-----------------------------------------------------------------------------
void HB::accumulateStatistics_()
{
  hbTotalNumberSuccessfulStepsTaken_ += analysisObject_->totalNumberSuccessfulStepsTaken_;
  hbTotalNumberFailedStepsAttempted_ += analysisObject_->totalNumberFailedStepsAttempted_;
  hbTotalNumberJacobiansEvaluated_ += analysisObject_->totalNumberJacobiansEvaluated_;
  hbTotalNumberIterationMatrixFactorizations_ += analysisObject_->totalNumberIterationMatrixFactorizations_;
  hbTotalNumberLinearSolves_ += analysisObject_->totalNumberLinearSolves_;
  hbTotalNumberFailedLinearSolves_ += analysisObject_->totalNumberFailedLinearSolves_;
  hbTotalNumberLinearIters_ += analysisObject_->totalNumberLinearIters_;
  hbTotalNumberResidualEvaluations_ += analysisObject_->totalNumberResidualEvaluations_;
  hbTotalNonlinearConvergenceFailures_ += analysisObject_->totalNonlinearConvergenceFailures_;
  hbTotalResidualLoadTime_ += analysisObject_->totalResidualLoadTime_;
  hbTotalJacobianLoadTime_ += analysisObject_->totalJacobianLoadTime_;
  hbTotalLinearSolutionTime_ += analysisObject_->totalLinearSolutionTime_;
}


//-----------------------------------------------------------------------------
// Function      : HB::runTol_
// Purpose       : Conducts transient run to determine right tolerance
//                 parameters for IC calculation
// Special Notes :
// Scope         : private
// Creator       : T. Mei, 1437, Electrical and Micro Modeling
// Creation Date : 02/23/09
//-----------------------------------------------------------------------------
bool HB::runTol_()
{
  bool returnValue = true;

  Xyce::lout() << " ***** Computing tolerance parameters for HB IC calculation....\n" << std::endl;

  N_TIA_TIAParams & tiaParams = anaManagerRCPtr_->tiaParams;
  N_TIA_TIAParams tiaParamsSave = tiaParams;

  // now try the real thing
  tiaParams.initialTime = 0;
  tiaParams.finalTime = period_;
  tiaParams.pauseTime = tiaParams.finalTime;
  tiaParams.resume = false;
  tiaParams.maxOrder = 1;

  // If start up periods are not being run then we can use this transient to compute HB ICs.
  if (!startUpPeriodsGiven_)
  {
    tiaParams.saveTimeStepsFlag = true;
  }

  // register new tiaParams with time integrator
  anaManagerRCPtr_->registerTIAParams (tiaParams);

  // Create a transient analysis object for this section.
  isTransient_ = true;
  analysisObject_ = Teuchos::rcp( new Transient( anaManagerRCPtr_.get() ) );
  analysisObject_->setAnalysisParams( N_UTL_OptionBlock() );
  //analysisObject_->registerLoader( loaderRCPtr_.get() );
  Teuchos::rcp_dynamic_cast<Transient>(analysisObject_)->resetForHB();
  returnValue = analysisObject_->run();

  if (!returnValue)
  {
    std::string msg = "Calculation of tolerance parameters failed for relErrorTol = "
                      + Teuchos::Utils::toString(tiaParams.relErrorTol) + ".\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0, msg );
  }


  int numPoints = anaManagerRCPtr_->getStepNumber();

  while((numPoints < (1.2*size_)) && (tiaParams.relErrorTol>= 1e-6))
  {
      std::string msg = "Tolerance parameters refined, re-running with relErrorTol = "
                        + Teuchos::Utils::toString(tiaParams.relErrorTol/10) + ".\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg );

      if (!startUpPeriodsGiven_)
      {
        // Clear the fast time data storage before performing the next transient
        RCP<N_TIA_DataStore> dsPtr = anaManagerRCPtr_->getTIADataStore();
        dsPtr->resetFastTimeData();
      }

      tiaParams.relErrorTol =  tiaParams.relErrorTol/10;
      anaManagerRCPtr_->registerTIAParams (tiaParams);

      // Create a transient analysis object for this section.
      analysisObject_ = Teuchos::rcp( new Transient( anaManagerRCPtr_.get() ) );
      analysisObject_->setAnalysisParams( N_UTL_OptionBlock() );
      //analysisObject_->registerLoader( loaderRCPtr_.get() );
      Teuchos::rcp_dynamic_cast<Transient>(analysisObject_)->resetForHB();
      bool retV = analysisObject_->run();

      if (!retV)
      {
        std::string msg = "Calculation of tolerance parameters failed for relErrorTol = "
                          + Teuchos::Utils::toString(tiaParams.relErrorTol) + ".\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0, msg );
      }
      returnValue = retV && returnValue;


    numPoints = anaManagerRCPtr_->getStepNumber();
      
  }

  // Add in simulation times
  accumulateStatistics_();

  // Reset parameters
  tiaParamsSave.relErrorTol  = tiaParams.relErrorTol;
  anaManagerRCPtr_->registerTIAParams (tiaParamsSave);

  // Reset transient flag before exiting.
  isTransient_ = false;

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::runStartupPeriods_()
// Purpose       : Runs normal transient problem through the requested
//                 number of startup periods
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool HB::runStartupPeriods_()
{
  bool returnValue = true;

  Xyce::lout() << "  ***** Computing " << startUpPeriods_ << " start up periods for HB IC calculation...." << std::endl;

  // need to advance time by startUpPeriods_ * period_
  N_TIA_TIAParams & tiaParams = anaManagerRCPtr_->tiaParams;
  N_TIA_TIAParams tiaParamsSave = tiaParams;

  // set DAE initial time = 0.0
  tiaParams.initialTime = 0.0;
  tiaParams.finalTime = startUpPeriods_ * period_;
  tiaParams.pauseTime = tiaParams.finalTime;

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB::runStartupPeriods_():  Advancing time through "
            << startUpPeriods_ << " startup periods"
            << " initialTime = " << tiaParams.initialTime
            << " finalTime = " << tiaParams.finalTime << std::endl;
#endif // Xyce_DEBUG_HB

  // tell the output manager to save this data as netlist.cir.startup.prn
  // anaManagerRCPtr_->outMgrPtr->setOutputFilenameSuffix(".startup");

  {
    Xyce::IO::OutputMgr::ActiveOutput x(*anaManagerRCPtr_->outMgrPtr);
    x.add(Xyce::IO::PrintType::HB_STARTUP);

    // register new tiaParams with time integrator
    anaManagerRCPtr_->registerTIAParams (tiaParams);

    // Create a transient analysis object for this section.
    isTransient_ = true;
    analysisObject_ = Teuchos::rcp( new Transient( anaManagerRCPtr_.get() ) );
    analysisObject_->setAnalysisParams( N_UTL_OptionBlock() );
    Teuchos::rcp_dynamic_cast<Transient>(analysisObject_)->resetForHB();
    returnValue = analysisObject_->run();
    isTransient_ = false;

    // Add in simulation times
    accumulateStatistics_();

    anaManagerRCPtr_->outMgrPtr->finishOutput();
  }

  // reset the output filename suffix
  // anaManagerRCPtr_->outMgrPtr->setOutputFilenameSuffix("");

  // put the dsPtr->currentSolutionPtr into dcOpSol and State Vec so that it
  // is used as our initial condition for the pending fast time scale runs
  RCP<N_TIA_DataStore> dsPtr = anaManagerRCPtr_->getTIADataStore();
  dcOpSolVecPtr_ = rcp( new N_LAS_Vector( *(dsPtr->currSolutionPtr) ));
  dcOpStateVecPtr_ = rcp( new N_LAS_Vector( *(dsPtr->currStatePtr) ));
  dcOpQVecPtr_ = rcp( new N_LAS_Vector( *(dsPtr->daeQVectorPtr) ));
  dcOpStoreVecPtr_ = rcp( new N_LAS_Vector( *(dsPtr->currStorePtr) ));

  // tell HB to start after this startup period
  tiaParamsSave.initialTime = startUpPeriods_ * period_;
  tiaParams_.initialTime = startUpPeriods_ * period_;
  anaManagerRCPtr_->registerTIAParams (tiaParamsSave);
  startUpPeriodsFinished_ = true;

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::runTransientIC_
// Purpose       : Conducts a regular transient run for HB initial conditions
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 10/03/2008
//-----------------------------------------------------------------------------
bool HB::runTransientIC_()
{
  bool returnValue = true;

  Xyce::lout() << " ***** Running transient to compute HB initial condition....\n" << std::endl;

  // this prevents extra DC op data from being printed.
  devInterfacePtr_->setMPDEFlag( true );

  if(saveIcData_)
  {
    // Keep the initial condition data
    // anaManagerRCPtr_->outMgrPtr->setOutputFilenameSuffix( ".hb_ic" );
  }

  // use an initial transient run to create a set of time points for the fast time scale
  N_TIA_TIAParams & tiaParams = anaManagerRCPtr_->tiaParams;
  N_TIA_TIAParams tiaParamsSave = tiaParams;

  //tiaParams.initialTime = 0; // should be tiaParams_.initialTime;
  tiaParams.initialTime = tiaParams_.initialTime;
  //tiaParams.finalTime = period_; // should be tiaParams_.initialTime + period_;
  tiaParams.finalTime = tiaParams_.initialTime + period_;
  tiaParams.pauseTime = tiaParams.finalTime;
  tiaParams.saveTimeStepsFlag = true;
  tiaParams.maxOrder = 1;

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB::runTransientIC_():  Advancing time from"
       << " initialTime = " << tiaParams.initialTime
       << " finalTime = " << tiaParams.finalTime << std::endl;
#endif // Xyce_DEBUG_HB

  // Initial conditions will be set if startup periods were run.
  if ( startUpPeriodsGiven_ )
  {
    tiaParams.NOOP = true;

    RCP<N_TIA_DataStore> dsPtr = anaManagerRCPtr_->getTIADataStore();
    *(dsPtr->nextSolutionPtr) = *(dcOpSolVecPtr_.get());
    *(dsPtr->nextStatePtr) = *(dcOpStateVecPtr_.get());
    *(dsPtr->daeQVectorPtr) = *(dcOpQVecPtr_.get());
    *(dsPtr->nextStorePtr) = *(dcOpStoreVecPtr_.get());
  }

  // register new tiaParams with time integrator
  anaManagerRCPtr_->registerTIAParams (tiaParams);

  // Create a transient analysis object for this section.
  isTransient_ = true;
  analysisObject_ = Teuchos::rcp( new Transient( anaManagerRCPtr_.get() ) );
  analysisObject_->setAnalysisParams( N_UTL_OptionBlock() );
  Teuchos::rcp_dynamic_cast<Transient>(analysisObject_)->resetForHB();
  returnValue = analysisObject_->run();
  isTransient_ = false;

  // Add in simulation times
  accumulateStatistics_();

  if(saveIcData_)
  {
    // reset suffix
    // anaManagerRCPtr_->outMgrPtr->setOutputFilenameSuffix( "" );
  }

  // restore the saved copy of time integration parameters.
  anaManagerRCPtr_->registerLoader( loaderRCPtr_.get() );

  tiaParamsSave.initialTime += period_;  // start HB problem after this transient init.
  anaManagerRCPtr_->registerTIAParams (tiaParamsSave);
  devInterfacePtr_->setMPDEFlag( false );

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::interpolateIC_()
// Purpose       : Tries to filter the fast time points from a transient run
//                 so that points are not too close together
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool HB::interpolateIC_()
{
  Xyce::lout() << " ***** Interpolating transient solution for IC calculation....\n" << std::endl;

  RCP<N_TIA_DataStore> dsPtr = anaManagerRCPtr_->getTIADataStore();
  int numPoints = dsPtr->timeSteps.size();

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB::interpolateIC_(): Initial transient run produced " << numPoints << " points." << std::endl;
#endif

  std::vector<int> goodIndicies;
  goodTimePoints_.resize(size_);

  double TimeStep = period_/static_cast<double>(size_);
  timeSteps_.push_back( TimeStep );
  for( int i = 0; i < size_; ++i )
  {
    goodTimePoints_[i] = tiaParams_.initialTime + static_cast<double>(i) * TimeStep;
  }
  fastTimes_ = goodTimePoints_;
  fastTimes_.resize(size_+1);
  fastTimes_[size_] = tiaParams_.initialTime + period_;

  int breakpoints = 0;          // need to keep track of how many breakpoints there are
  int startIndex = 0;

  // always keep first point
  goodIndicies.push_back(startIndex);
  int GoodTimePointIndex = startIndex + 1;

  for( int i=startIndex; i < numPoints - 1 ; i++ )
  {
    // count up breakpoints
    if( dsPtr->timeStepsBreakpointFlag[i] == true )
    {
      breakpoints++;
    }

#ifdef Xyce_DEBUG_HB
    if( debugLevel > 0 )
    {
      Xyce::dout() << "\t\t timeStep[ " << i << " ] = " << dsPtr->timeSteps[i];
      if( dsPtr->timeStepsBreakpointFlag[i] == true )
      {
        Xyce::dout() << "  Breakpoint";
      }
      Xyce::dout() << std::endl;
    }
#endif
    while( ( GoodTimePointIndex < size_ )  && (dsPtr->timeSteps[i] <= goodTimePoints_[GoodTimePointIndex]) && (goodTimePoints_[GoodTimePointIndex] < dsPtr->timeSteps[i+1]))
    {
      // found a good point so save the index
      goodIndicies.push_back( i );
      GoodTimePointIndex =  GoodTimePointIndex+1;
    }
  }

  for(int i=0; i<size_; i++ )
  {
    int currentIndex = goodIndicies[i];
    N_LAS_Vector * firstSolVecPtr = dsPtr->fastTimeSolutionVec[currentIndex];
    N_LAS_Vector * secondSolVecPtr = dsPtr->fastTimeSolutionVec[currentIndex+1];

    N_LAS_Vector * firstStateVecPtr = dsPtr->fastTimeStateVec[currentIndex];
    N_LAS_Vector * secondStateVecPtr = dsPtr->fastTimeStateVec[currentIndex+1];

    N_LAS_Vector * firstQVecPtr = dsPtr->fastTimeQVec[currentIndex];
    N_LAS_Vector * secondQVecPtr = dsPtr->fastTimeQVec[currentIndex+1];

    N_LAS_Vector * firstStoreVecPtr = dsPtr->fastTimeStoreVec[currentIndex];
    N_LAS_Vector * secondStoreVecPtr = dsPtr->fastTimeStoreVec[currentIndex+1];

    double fraction = (goodTimePoints_[i] -  dsPtr->timeSteps[currentIndex])/(dsPtr->timeSteps[currentIndex+1] -  dsPtr->timeSteps[currentIndex]);

    RCP<N_LAS_Vector> InterpICSolVecPtr = rcp( new N_LAS_Vector( *secondSolVecPtr ) );
    RCP<N_LAS_Vector> InterpICStateVecPtr = rcp( new N_LAS_Vector( *secondStateVecPtr ) );
    RCP<N_LAS_Vector> InterpICQVecPtr = rcp( new N_LAS_Vector( *secondQVecPtr ) );
    RCP<N_LAS_Vector> InterpICStoreVecPtr = rcp( new N_LAS_Vector( *secondStoreVecPtr ) );

    InterpICSolVecPtr->putScalar(0.0);
    InterpICStateVecPtr->putScalar(0.0);
    InterpICQVecPtr->putScalar(0.0);
    InterpICStoreVecPtr->putScalar(0.0);

    InterpICSolVecPtr->linearCombo(-1.0, *firstSolVecPtr, 1.0, *secondSolVecPtr );
    InterpICSolVecPtr->linearCombo(1.0, *firstSolVecPtr, fraction , *InterpICSolVecPtr);

    InterpICStateVecPtr->linearCombo(-1.0, *firstStateVecPtr, 1.0, *secondStateVecPtr );
    InterpICStateVecPtr->linearCombo(1.0, *firstStateVecPtr, fraction , *InterpICStateVecPtr);

    InterpICQVecPtr->linearCombo(-1.0, *firstQVecPtr, 1.0, *secondQVecPtr );
    InterpICQVecPtr->linearCombo(1.0, *firstQVecPtr, fraction , *InterpICQVecPtr);

    InterpICStoreVecPtr->linearCombo(-1.0, *firstStoreVecPtr, 1.0, *secondStoreVecPtr );
    InterpICStoreVecPtr->linearCombo(1.0, *firstStoreVecPtr, fraction , *InterpICStoreVecPtr);

    goodSolutionVec_.push_back(InterpICSolVecPtr);
    goodStateVec_.push_back(InterpICStateVecPtr);
    goodQVec_.push_back(InterpICQVecPtr);
    goodStoreVec_.push_back(InterpICStoreVecPtr);
  }

  // Clean up the fast time data since we are finished computing the initial condition.
  // The fast time data can take a considerable amount of memory for large problems.
  dsPtr->resetFastTimeData();

  return true;
}

} // namespace Analysis
} // namespace Xyce

