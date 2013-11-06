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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_ParameterSet.C,v $
//
// Purpose        : Interface to Xyce vectors for NOX.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.88.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:48 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>



// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include<list>

// ----------   Xyce Includes   ----------

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "N_NLS_NOX_Interface.h"
#include "N_NLS_NOX_XyceTests.h"
#include "N_NLS_NOX_FastTests.h"
#include "N_NLS_NOX_WeightedUpdateTest.h"
#include "N_NLS_NOX_UpdateTooSmallTest.h"
#include "N_NLS_NOX_NearConvergenceTest.h"
#include "N_NLS_NOX_StagnationTest.h"
#include "N_NLS_NOX_LargeUpdateTest.h"
#include "N_ERH_ErrorMgr.h"
#include "N_UTL_Param.h"
#include "N_UTL_OptionBlock.h"
#include "LOCA.H"
#include "N_NLS_ReturnCodes.h"
#include "N_LOA_Loader.h"

// The following are needed to build AugmentLinSys strategies
#include "N_NLS_NOX_AugmentLinSys.h"
#include "N_NLS_NOX_AugmentLinSys_PseudoTransient.h"
#include "N_NLS_NOX_AugmentLinSys_GStepping.h"
#include "N_NLS_NOX_AugmentLinSys_OPStart.h"
#include "N_NLS_NOX_AugmentLinSys_IC.h"
#include "N_NLS_NOX_AugmentLinSys_IC_Gmin.h"
#include "N_LAS_Builder.h"
#include "N_LAS_System.h"
#include "N_LAS_QueryUtil.h"
#include "Epetra_MapColoring.h"

// -----------  Forward declarations  -------
class N_PDS_Comm;

using namespace N_NLS_NOX;

//-----------------------------------------------------------------------------
// Function      : ParameterSet::ParameterSet
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
ParameterSet::ParameterSet(AnalysisMode mode) :
  allParams_(Teuchos::rcp(new Teuchos::ParameterList)),
  noxParams_(allParams_->sublist("NOX")),
  locaParams_(allParams_->sublist("LOCA")),
  debugParams_(allParams_->sublist("DEBUG")),
  comboPtr_(Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR))),
  isParamsSet_(false),
  isStatusTestsSet_(false),
  mode_(mode),
  noxSolver(0),
  voltageListType_(VLT_None),
  gstepping_minimum_conductance_(0.0),
  savedLocaOptions_(false),
#ifdef Xyce_DEBUG_NONLINEAR
  debugLevel_(0),
  debugMinTimeStep_(0),
  debugMaxTimeStep_(N_UTL_MachineDependentParams::IntMax()),
  debugMinTime_(0.0),
  debugMaxTime_(N_UTL_MachineDependentParams::DoubleMax()),
  screenOutputFlag_(false),
#endif
  voltageScaleFactor_(1.0)
{
  // Add the main status test to the list of tests to delete in dtor
  tests_.push_back(comboPtr_);

  // Default printing options
#ifdef Xyce_VERBOSE_NOX
  noxParams_.sublist("Printing")
    .set("Output Information",
		  NOX::Utils::Error +
		  NOX::Utils::Warning +
		  NOX::Utils::OuterIteration +
		  NOX::Utils::OuterIterationStatusTest +
		  NOX::Utils::InnerIteration +
		  NOX::Utils::Details +
		  NOX::Utils::StepperIteration +
		  NOX::Utils::StepperDetails +
		  NOX::Utils::Parameters
		  );

#else
#ifdef Xyce_VERBOSE_NONLINEAR
  noxParams_.sublist("Printing")
    .set("Output Information",
		  NOX::Utils::Error+
		  NOX::Utils::Warning +
		  NOX::Utils::OuterIteration +
		  NOX::Utils::StepperIteration
		  );
#else
  noxParams_.sublist("Printing")
    .set("Output Information", NOX::Utils::Error);
#endif
#endif

  // Defaults that are mode dependent
  switch (mode_)
  {
  case TRANSIENT:
    // These values correspond to N_NLS_DampedNewton.C; see the constructor
    statusTestParams_.set("ABSTOL", 1.0e-6);
    statusTestParams_.set("RELTOL", 1.0e-2);
    statusTestParams_.set("DELTAXTOL", 0.33);
    statusTestParams_.set("RHSTOL", 1.0e-2);
    statusTestParams_.set("MAXSTEP", 20);
    noxParams_.set("Nonlinear Solver", "Line Search Based");
    noxParams_.sublist("Line Search").set("Method", "Full Step");
    //noxParams_.sublist("Line Search").set("Method", "Polynomial");
    //noxParams_.sublist("Line Search").sublist("Polynomial")
    //  .set("Max Iters", 2);
    //noxParams_.sublist("Line Search").sublist("Polynomial")
    //  .set("Interpolation Type", "Quadratic");
    //noxParams_.sublist("Direction").sublist("Newton")
    //.sublist("Linear Solver").set("Tolerance", 1.0e-9);
    break;
  case HB:
    // These values correspond to N_NLS_DampedNewton.C; see the constructor
    statusTestParams_.set("ABSTOL", 1.0e-9);
    statusTestParams_.set("RELTOL", 1.0e-3);
    statusTestParams_.set("DELTAXTOL", 1.0);
    statusTestParams_.set("RHSTOL", 1.0e-4);
    statusTestParams_.set("MAXSTEP", 200);
    noxParams_.set("Nonlinear Solver", "Line Search Based");
    noxParams_.sublist("Line Search").set("Method", "Full Step");
    //noxParams_.sublist("Line Search").set("Method", "Polynomial");
    //noxParams_.sublist("Line Search").sublist("Polynomial")
    //  .set("Max Iters", 2);
    //noxParams_.sublist("Line Search").sublist("Polynomial")
    //  .set("Interpolation Type", "Quadratic");
    //noxParams_.sublist("Direction").sublist("Newton")
    //.sublist("Linear Solver").set("Tolerance", 1.0e-9);
    break;
  default:
    statusTestParams_.set("ABSTOL", 1.0e-12);
    statusTestParams_.set("RELTOL", 1.0e-3);
    statusTestParams_.set("DELTAXTOL", 1.0);
    statusTestParams_.set("RHSTOL", 1.0e-6);
    statusTestParams_.set("MAXSTEP", 200);
    noxParams_.set("Nonlinear Solver", "Line Search Based");
    noxParams_.sublist("Line Search").set("Method", "Full Step");
    noxParams_.sublist("Direction").sublist("Newton")
      .sublist("Linear Solver").set("Tolerance", 1.0e-12);
    break;
  }

  // Parameters that should always be set
  noxParams_.sublist("Line Search").sublist("Polynomial")
    .set("Recovery Step Type", "Last Computed Step");


  // Set default loca options in case this is a loca run.
  Teuchos::ParameterList& stepperList = locaParams_.sublist("Stepper");
  Teuchos::ParameterList& predictorList = locaParams_.sublist("Predictor");
  Teuchos::ParameterList& stepSizeList = locaParams_.sublist("Step Size");

  stepperList.set("Continuation Method", "Natural");
  stepperList.set("Skip df/dp", true);
  predictorList.set("Method", "Tangent");
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::~ParameterSet
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
ParameterSet::~ParameterSet()
{
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::setOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::setOptions(const N_UTL_OptionBlock& OB)
{
  // Parse the option block
  bool parseok = parseOptionBlock_(OB);
  if (!parseok)
  {
    return false;
  }

  isParamsSet_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::setOutputOptions
// Purpose       : Set output for parallel runs
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::setOutputOptions(int myPID, int outputProcess)
{
  noxParams_.sublist("Printing").set("MyPID", myPID);
  noxParams_.sublist("Printing").set("Output Processor",
					      outputProcess);
  locaParams_.sublist("Utilities").set("MyPID", myPID);
  locaParams_.sublist("Utilities").set("Output Processor",
						outputProcess);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::createStatusTests
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::createStatusTests(N_LAS_Vector** currSolVectorPtrPtr,
				     N_LOA_Loader& loader, vector<char> & varTypeVec
#ifdef Xyce_NLS_MASKED_WRMS_NORMS
             , bool nonTrivialDeviceMaskFlag,  N_LAS_Vector * maskVectorPtr)
#else
             )
#endif

{
  /*
  // Test 1 - Make sure the residual isn't too small, hardwired tolerances

  // Test 1a - Max Norm F is less than machine precision
  NOX::StatusTest::NormF* test1a =
    new NOX::StatusTest::NormF(N_UTL_MachineDependentParams::MachineEpsilon(),
			       NOX::Abstract::Vector::MaxNorm,
			       NOX::StatusTest::NormF::Unscaled);
  tests_.push_back(test1a);

  // Test 1b - Max Norm F is less than requested tolerance
  NOX::StatusTest::NormF* test1b =
    new NOX::StatusTest::NormF(statusTestParams_.get("RHSTOL", 1.0e-6),
			       NOX::Abstract::Vector::MaxNorm,
			       NOX::StatusTest::NormF::Unscaled);
  tests_.push_back(test1b);

  NOX::StatusTest::Combo* test1 =
    new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND, *test1a, *test1b);
  tests_.push_back(test1);
  comboPtr_->addStatusTest(*test1);


  // Test 2 - Normal convergence based on rhs residual (2a) and
  // update norm (2b).

  // Test 2a - Max Residual
  NOX::StatusTest::NormF* test2a =
    new NOX::StatusTest::NormF(statusTestParams_.get("RHSTOL", 1.0e-6),
			       NOX::Abstract::Vector::MaxNorm,
			       NOX::StatusTest::NormF::Unscaled);
  tests_.push_back(test2a);

  // Test 2b - weighted update
  bool isTransient = false;
  if (mode_ == TRANSIENT)
    isTransient = true;
  */
  /*
  N_NLS_NOX::WeightedUpdateTest* test2b =
    new N_NLS_NOX::WeightedUpdateTest(currSolVectorPtrPtr,
		      statusTestParams_.get("ABSTOL", 1.0e-12),
		      statusTestParams_.get("RELTOL", 1.0e-3),
		      statusTestParams_.get("DELTAXTOL", 1.0),
		      isTransient);
  tests_.push_back(test2b);
 */
  /*
  //RPP Hack for test2b to get Homotopy to work
  NOX::StatusTest::NormF* test2b =
    new NOX::StatusTest::NormF(statusTestParams_.get("RHSTOL", 1.0e-6),
			       NOX::Abstract::Vector::MaxNorm,
			       NOX::StatusTest::NormF::Unscaled);
  tests_.push_back(test2b);


  NOX::StatusTest::Combo* test2 =
    new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND, *test2a, *test2b);
  tests_.push_back(test2);
  comboPtr_->addStatusTest(*test2);

  // Test 3 - Near Convergence - Hit max iterations but residual
  // and convergence rate indicate we may be near a converged solution.
  // Therefore, let the time stepper decide whether or  not the step is ok.
  // Transient mode ONLY!
  NOX::StatusTest::Generic* test3 = 0;
  if (mode_ == TRANSIENT)
  {
    test3 = new N_NLS_NOX::NearConvergenceTest(statusTestParams_.get("MAXSTEP", 200), 1.0);
    tests_.push_back(test3);
    comboPtr_->addStatusTest(*test3);
  }
  */
  /*
  // Test 4 - Update is too small
  N_NLS_NOX::UpdateTooSmallTest* test4 =
    new N_NLS_NOX::UpdateTooSmallTest(*test2b, 1.0e-6);
  tests_.push_back(test4);
  comboPtr_->addStatusTest(*test4);
  */
  /*
  // Test 5 - Max nonlinear iterations (if transient, this will be checked
  // in the NearConvergence test (#3)
  if (mode_ != TRANSIENT)
  {
    NOX::StatusTest::MaxIters* test5 =
      new NOX::StatusTest::MaxIters(statusTestParams_.get("MAXSTEP", 200));
    tests_.push_back(test5);
    comboPtr_->addStatusTest(*test5);
  }

  // Test 6 - update is too big
  N_NLS_NOX::LargeUpdateTest* test6 =
    new N_NLS_NOX::LargeUpdateTest(0.5*N_UTL_MachineDependentParams::DoubleMax());
  tests_.push_back(test6);
  comboPtr_->addStatusTest(*test6);

  // Test 7 - Stall in the convergence rate. Transient mode ONLY!
  if (mode_ == TRANSIENT)
  {
    N_NLS_NOX::Stagnation* test7 = new N_NLS_NOX::Stagnation(50);
    tests_.push_back(test7);
    comboPtr_->addStatusTest(*test7);
  }

  */

  // RPP: moving inside allTests so we don't have to dynamic cast
  // NaN/Inf check on the residual vector 2-norm
  //NOX::StatusTest::FiniteValue* fvTest = new NOX::StatusTest::FiniteValue;
  //tests_.push_back(fvTest);
  //comboPtr_->addStatusTest(*fvTest);

  // Tests All - replaces tests 1-7
  bool isTransient = false;
  if (mode_ == TRANSIENT)
  {
    isTransient = true;
  }

  // make the correct testing object
  Teuchos::RefCountPtr<N_NLS_NOX::XyceTests> allTests;
  if( statusTestParams_.get("FASTTESTS",false) )
  {
    // make the FastTests object
    allTests = Teuchos::rcp(new N_NLS_NOX::FastTests(isTransient,
              statusTestParams_.get("RHSTOL", 1.0e-6),
              N_UTL_MachineDependentParams::MachineEpsilon(),
              currSolVectorPtrPtr,
              statusTestParams_.get("ABSTOL", 1.0e-12),
              statusTestParams_.get("RELTOL", 1.0e-3),
              statusTestParams_.get("DELTAXTOL", 1.0),
              statusTestParams_.get("MAXSTEP", 200),
              0.9,
              1.0,
              0.5*N_UTL_MachineDependentParams::DoubleMax(),
              1.0e-3,
              5,
              statusTestParams_.get("ENFORCEDEVICECONV", 1),
              statusTestParams_.get("SMALLUPDATETOL", 1.0e-6),
              &loader,
              varTypeVec,
              statusTestParams_.get("VOLTZEROTOL", 1.0e-6),
              statusTestParams_.get("CURRZEROTOL", 1.0e-6)
  #ifdef Xyce_NLS_MASKED_WRMS_NORMS
              , nonTrivialDeviceMaskFlag, maskVectorPtr
  #endif
              ));
	}
	else
	{
	  // here we make the default XyceTests object
     allTests = Teuchos::rcp(new N_NLS_NOX::XyceTests(isTransient,
              statusTestParams_.get("RHSTOL", 1.0e-6),
              N_UTL_MachineDependentParams::MachineEpsilon(),
              currSolVectorPtrPtr,
              statusTestParams_.get("ABSTOL", 1.0e-12),
              statusTestParams_.get("RELTOL", 1.0e-3),
              statusTestParams_.get("DELTAXTOL", 1.0),
              statusTestParams_.get("MAXSTEP", 200),
              0.9,
              1.0,
              0.5*N_UTL_MachineDependentParams::DoubleMax(),
              1.0e-3,
              5,
              statusTestParams_.get("ENFORCEDEVICECONV", 1),
              statusTestParams_.get("SMALLUPDATETOL", 1.0e-6),
              &loader
  #ifdef Xyce_NLS_MASKED_WRMS_NORMS
              , nonTrivialDeviceMaskFlag, maskVectorPtr
  #endif
              ));
	}
  tests_.push_back(allTests);
  comboPtr_->addStatusTest(allTests);
  isStatusTestsSet_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getStatusTests
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<NOX::StatusTest::Generic> ParameterSet::getStatusTests()
{
  if (!isStatusTestsSet_)
  {
    const string message = "Error: N_NLS::NOX::ParameterSet::getStatusTests() - Status tests are not set!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return comboPtr_;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getAllParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<Teuchos::ParameterList> ParameterSet::getAllParams()
{
  return allParams_;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getNoxParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<Teuchos::ParameterList> ParameterSet::getNoxParams()
{
  return Teuchos::rcp(&noxParams_, false);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getLocaParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<Teuchos::ParameterList> ParameterSet::getLocaParams()
{
  return Teuchos::rcp(&locaParams_, false);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getDebugParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<Teuchos::ParameterList> ParameterSet::getDebugParams()
{
  return Teuchos::rcp(&debugParams_, false);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::unsupportedOption_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void ParameterSet::unsupportedOption_(const string& tag)
{
  const string warning =
    "Tag \"" + tag + "\" is unsupported by the NOX interface at this time.\n";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, warning);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::parseOptionBlock_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::parseOptionBlock_(const N_UTL_OptionBlock& OB)
{
  /*
     RPP: Some parameters can't be immediately set in the nox
     parameter list because they are in a sublist dependent upon
     another parameter being set first.  We have to store these
     parameters until the entire option block is parsed and then set
     them in the correct sublist.  These parameters are listed below.
  */
  int maxSearchStep = 2;
  int in_Forcing = 0;
  double AZ_tol = 1.0e-12;
  int recoveryStepType = 0;
  double recoveryStep = 1.0;
  int memory = 400;


  // Loop over all parameters in the option block
  for (list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
       it_tpL != OB.getParams().end(); ++ it_tpL)
  {
    const string tag = it_tpL->uTag();


    // Parameters for nonlinear convergence tests
    if (tag == "ABSTOL")
    {
      statusTestParams_.set("ABSTOL", it_tpL->dVal());
    }
    else if (tag == "RELTOL")
    {
      statusTestParams_.set("RELTOL", it_tpL->dVal());
    }
    else if (tag == "DELTAXTOL")
    {
      statusTestParams_.set("DELTAXTOL", it_tpL->dVal());
    }
    else if (tag == "RHSTOL")
    {
      statusTestParams_.set("RHSTOL", it_tpL->dVal());
    }
    else if (tag == "MAXSTEP")
    {
      statusTestParams_.set("MAXSTEP", it_tpL->iVal());
    }
    else if (tag == "SMALLUPDATETOL")
    {
      statusTestParams_.set("SMALLUPDATETOL", it_tpL->dVal());
    }

    // Check devices for convergence (by calls to cktloader)
    else if (tag == "ENFORCEDEVICECONV")
    {
      statusTestParams_.set("ENFORCEDEVICECONV", it_tpL->iVal());
    }
	  else if (tag == "FASTTESTS")
    {
      statusTestParams_.set("FASTTESTS", it_tpL->bVal());
    }
    else if (tag == "VOLTZEROTOL")
    {
      statusTestParams_.set("VOLTZEROTOL", it_tpL->dVal());
    }
    else if (tag == "ABSZEROTOL")
    {
      statusTestParams_.set("ABSZEROTOL", it_tpL->dVal());
    }

    // Unsupported options
    else if (tag == "LINOPT")
    {
      unsupportedOption_(tag);
    }
    else if (tag == "CONSTRAINTBT")
    {
      unsupportedOption_(tag);
    }
    else if (tag == "CONSTRAINTMAX")
    {
      unsupportedOption_(tag);
    }
    else if (tag == "CONSTRAINTMIN")
    {
      unsupportedOption_(tag);
    }
    else if (tag == "CONSTRAINTCHANGE")
    {
      unsupportedOption_(tag);
    }
    else if (tag == "NORMLVL")
    {
      if (it_tpL->iVal() != 2)
        unsupportedOption_(tag);
    }
#ifdef Xyce_DEBUG_NONLINEAR
    else if (tag == "DEBUGLEVEL")
    {
      debugLevel_ = it_tpL->iVal();
    }
    else if (tag == "SCREENOUTPUT")
    {
      screenOutputFlag_ = it_tpL->bVal();
    }
    else if (tag == "DEBUGMINTIMESTEP")
    {
      debugMinTimeStep_ = it_tpL->iVal();
    }
    else if (tag == "DEBUGMAXTIMESTEP")
    {
      debugMaxTimeStep_ = it_tpL->iVal();
    }
    else if (tag == "DEBUGMINTIME")
    {
      debugMinTime_ = it_tpL->dVal();
    }
    else if (tag == "DEBUGMAXTIME")
    {
      debugMaxTime_ = it_tpL->dVal();
    }
#endif
    else if (tag == "DLSDEBUG")
      unsupportedOption_(tag);

    // Nonlinear Strategy
    else if (tag == "NLSTRATEGY")
    {
      int val = it_tpL->iVal();
      if (val == 0)		// Newton
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction").set("Method", "Newton");
      }
      else if (val == 1) 	// Steepest descent
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction")
          .set("Method", "Steepest Descent");
      }
      else if (val == 2)	// Trust Region
      {
        noxParams_.set("Nonlinear Solver", "Trust Region Based");
      }
      else if (val == 3) 	// Modified Newton
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction")
          .set("Method", "Modified-Newton");
      }
      else if (val == 4)	// BFGS
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction").set("Method", "Quasi-Newton");
      }
      else if (val == 5)	// Broyden
      {
        noxParams_.set("Nonlinear Solver", "Line Search Based");
        noxParams_.sublist("Direction").set("Method", "Broyden");
      }
      else if (val == 6)	// Tensor
      {
        noxParams_.set("Nonlinear Solver", "Tensor Based");
        noxParams_.sublist("Direction").set("Method", "Tensor");
        noxParams_.sublist("Direction").sublist("Tensor")
          .sublist("Linear Solver").set("Compute Step", "Newton");
        noxParams_.sublist("Direction").sublist("Tensor")
          .sublist("Linear Solver").set("Reorthogonalize", "Always");
        noxParams_.sublist("Line Search").set("Method", "Tensor");
        noxParams_.sublist("Line Search").sublist("Tensor")
          .set("Submethod", "Full Step");
      }
      else if (val == 7)	// Fast Newton Direction
      {
        // RPP: No longer supported
        const string warning =
          "NLStrategy = 7 is no longer supported.\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, warning);
      }
      else if (val == 8)	// Newton/Steepest Descent Combo Direction
      {
        // RPP: No longer supported
        const string warning =
          "NLStrategy = 8 is no longer supported.\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, warning);
      }
      else
      {
        const string warning =
          "NLStrategy is not found!\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, warning);
      }
    }

    // Line search method
    else if (tag == "SEARCHMETHOD")
    {
      int val = it_tpL->iVal();
      if (val == 0)
      {
        noxParams_.sublist("Line Search").set("Method", "Full Step");
      }
      else if (val == 1)
      {
        noxParams_.sublist("Line Search").set("Method", "Backtrack");
      }
      else if (val == 2)
      {
        noxParams_.sublist("Line Search").set("Method", "Polynomial");
        noxParams_.sublist("Line Search").sublist("Polynomial")
          .set("Interpolation Type", "Quadratic");
      }
      else if (val == 3)
      {
        noxParams_.sublist("Line Search").set("Method", "Polynomial");
        noxParams_.sublist("Line Search").sublist("Polynomial")
          .set("Interpolation Type", "Cubic");
      }
      else if (val == 4)
      {
        noxParams_.sublist("Line Search")
          .set("Method", "More'-Thuente");
      }
      else
      {
        std::ostringstream ost;
        ost << "N_NLS_NOX::ParameterSet::parseOptionBlock_ - "
            << "SEARCHMETHOD = " << val
            << " not supported by Xyce at this time." << endl;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, ost.str());
      }
    }

    // Trust Region auxiliary parameters
    else if (tag == "TRMINRADIUS")
    {
      noxParams_.sublist("Trust Region")
        .set("Minimum Trust Region Radius", it_tpL->dVal());
    }
    else if (tag == "TRMAXRADIUS")
    {
      noxParams_.sublist("Trust Region")
        .set("Maximum Trust Region Radius", it_tpL->dVal());
    }
    else if (tag == "TRMINIMPROVEMENTRATIO")
    {
      noxParams_.sublist("Trust Region")
        .set("Minimum Improvement Ratio", it_tpL->dVal());
    }
    else if (tag == "TRCONTRACTIONRATIO")
    {
      noxParams_.sublist("Trust Region")
        .set("Contraction Trigger Ratio", it_tpL->dVal());
    }
    else if (tag == "TRCONTRACTIONFACTOR")
    {
      noxParams_.sublist("Trust Region")
        .set("Contraction Factor", it_tpL->dVal());
    }
    else if (tag == "TREXPANSIONRATIO")
    {
      noxParams_.sublist("Trust Region")
        .set("Expansion Trigger Ratio", it_tpL->dVal());
    }
    else if (tag == "TREXPANSIONFACTOR")
    {
      noxParams_.sublist("Trust Region")
        .set("Expansion Factor", it_tpL->dVal());
    }
    else if (tag == "TRRECOVERYSTEP")
    {
      noxParams_.sublist("Trust Region")
        .set("Recovery Step", it_tpL->dVal());
    }

    // RPP: Why this is here???
    else if (tag == "NOX")
    {
      // do nothing (this option is handled in the manager)
    }

    // LOCA Continuation control
    // 0 = Nox solve (no continuation)
    // 1 = Natural Parameter Continuation
    // 2 = Mosfet Specific Dual Parameter Continuation
    // 3 = gmin stepping.
    // 33 = Artificial Parameter Continuation
    else if (tag == "CONTINUATION")
    {
      if (it_tpL->isNumeric())
      {
        noxSolver = it_tpL->iVal();
      }
      else
      {
        ExtendedString p(it_tpL->sVal());
        p.toUpper();
        if (p.substr(0,4) == "STAN")
        {
          noxSolver = 0;
        }
        else if (p.substr(0,3) == "NAT")
        {
          noxSolver = 1;
        }
        else if (p.substr(0,3) == "MOS")
        {
          noxSolver = 2;
        }
        else if (p.substr(0,4) == "GMIN")
        {
          noxSolver = 3;
        }
        else if (p.substr(0,6) == "NEWMOS")
        {
          noxSolver = 4;
        }
        else if (p.substr(0,9) == "BSIM3INV1")
        {
          noxSolver = 5;
        }
        else if (p.substr(0,9) == "BSIM3INV2")
        {
          noxSolver = 6;
        }
        else if (p.substr(0,9) == "BLOCKGAIN")
        {
          noxSolver = 7;
        }
        else if (p.substr(0,4) == "TEST")
        {
          noxSolver = 8;
        }
        else if (p.substr(0,6) == "PSEUDO")
        {
          noxSolver = 9;
        }
        else if (p.substr(0,5) == "POWER")
        {
          noxSolver = 10;
        }
        else if (p.substr(0,3) == "ART")
        {
          noxSolver = 33;
        }
        else
        {
          string message = "Unknown specification in .options for 'continuation': ";
          message += it_tpL->sVal();
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
        }

	// Bail out if option 33 (artif. homotopy) is chosen but
	// required support is not built in.
#ifndef Xyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT
        if (noxSolver == 33)
        {
          const string message = "Error: N_NLS_NOX::ParameterSet::parseOptionBlock() - option \"continuation=33\" requires Xyce to be built with -DXyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT in the configure script.  Please rebuild Xyce or choose a different continuation method.\n";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
        }
#endif
      }
    }

    // Parameters that can't be set in the list until all options
    // have been parsed
    else if (tag == "MAXSEARCHSTEP")
    {
      maxSearchStep = it_tpL->iVal();
    }
    else if (tag == "IN_FORCING")
    {
      in_Forcing = it_tpL->iVal();
    }
    else if (tag == "AZ_TOL")
    {
      AZ_tol = it_tpL->dVal();
    }
    else if (tag == "RECOVERYSTEPTYPE")
    {
      recoveryStepType = it_tpL->iVal();
    }
    else if (tag == "RECOVERYSTEP")
    {
      recoveryStep = it_tpL->dVal();
    }
    else if (tag == "MEMORY")
    {
      memory = it_tpL->iVal();
    }

    // Parameters that can't be set in the list until all options
    // have been parsed
    else if (tag == "MAXSEARCHSTEP")
    {
      maxSearchStep = it_tpL->iVal();
    }
    else if (tag == "IN_FORCING")
    {
      in_Forcing = it_tpL->iVal();
    }
    else if (tag == "AZ_TOL")
    {
      AZ_tol = it_tpL->dVal();
    }
    else if (tag == "RECOVERYSTEPTYPE")
    {
      recoveryStepType = it_tpL->iVal();
    }
    else if (tag == "RECOVERYSTEP")
    {
      recoveryStep = it_tpL->dVal();
    }
    else if (tag == "MEMORY")
    {
      memory = it_tpL->iVal();
    }

    // Warn user about unrecognized solver option
    else
    {
      string tmp = "  Warning: " + tag +
      " is not a recognized nonlinear solver option.\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0,  tmp);
    }

  } // end loop over all options in block


  /*
     Set parameters that are dependent upon other parameters

     RPP: Some parameters can't be immediately set in the nox
     parameter list because they are in a sublist dependent upon
     another parameter being set first.  We have to store these
     parameters until the entire option block is parsed and then set
     them in the correct sublist.
  */
  string directionMethod =
    noxParams_.sublist("Direction").get("Method", "Newton");

  string lineSearchMethod =
    noxParams_.sublist("Line Search").get("Method", "Full Step");

  // MAXSEARCHSTEP
  noxParams_.sublist("Line Search").sublist(lineSearchMethod)
    .set("Max Iters", maxSearchStep);

  // In_FORCING
  if (in_Forcing == 0)
  {
    noxParams_.sublist("Direction").sublist(directionMethod)
      .set("Forcing Term Method", "Constant");
  }
  else if (in_Forcing == 1)
  {
    noxParams_.sublist("Direction").sublist(directionMethod)
      .set("Forcing Term Method", "Type 1");
    noxParams_.sublist("Direction").sublist("Newton")
      .set("Forcing Term Minimum Tolerance", AZ_tol);
  }
  else
  {
    noxParams_.sublist("Direction").sublist(directionMethod)
      .set("Forcing Term Method", "Type 2");
    noxParams_.sublist("Direction").sublist("Newton")
      .set("Forcing Term Minimum Tolerance", AZ_tol);
  }

  // RECOVERYSTEPTYPE
  if (recoveryStepType == 1)
  {
    // Recovery the NOX way
    noxParams_.sublist("Line Search").sublist(lineSearchMethod)
      .set("Recovery Step Type", "Constant");
  }
  else
  {
    // Recovery the Xyce way
    noxParams_.sublist("Line Search").sublist(lineSearchMethod)
      .set("Recovery Step Type", "Last Computed Step");
  }

  // RECOVERYSTEP
  noxParams_.sublist("Line Search").sublist(lineSearchMethod)
    .set("Recovery Step", recoveryStep);

  // MEMORY
  if ((directionMethod == "Quasi-Newton") ||
      (directionMethod == "Broyden"))
    noxParams_.sublist("Direction").sublist(directionMethod)
      .set("Memory", memory);


  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool ParameterSet::setLocaOptions(const N_UTL_OptionBlock& OB, bool saveCopy)
{
  Teuchos::ParameterList& stepperList = locaParams_.sublist("Stepper");
  Teuchos::ParameterList& predictorList = locaParams_.sublist("Predictor");
  Teuchos::ParameterList& stepSizeList = locaParams_.sublist("Step Size");

  bool stepperGiven=false;
  bool predictorGiven=false;

  if (saveCopy)
  {
    savedLocaOptions_ = true;
    savedLocaOB_ = OB; // save a copy to re-assert defaults later, if needed.
  }

  for (list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
       it_tpL != OB.getParams().end(); ++ it_tpL)
  {
    const string tag = it_tpL->uTag();
    string baseTag=tag.substr(0,8);
    bool isVectorParam=false;
    if (baseTag == "CONPARAM" ||
        baseTag == "MINVALUE" ||
        baseTag == "MAXVALUE" ||
        baseTag == "INITIALV" ||
        baseTag == "INITIALS" ||
        baseTag == "MINSTEPS" ||
        (baseTag == "MAXSTEPS" && tag.substr(0,11)=="MAXSTEPSIZE")||
        baseTag == "AGGRESSI")
    {
      if (baseTag == "INITIALV")
      {
        baseTag = "INITIALVALUE";
      }
      else if (baseTag == "INITIALS")
      {
        baseTag = "INITIALSTEPSIZE";
      }
      else if (baseTag == "MINSTEPS")
      {
        baseTag = "MINSTEPSIZE";
      }
      else if (baseTag == "MAXSTEPS")
      {
        baseTag = "MAXSTEPSIZE";
      }
      else if (baseTag == "AGGRESSI")
      {
        baseTag = "AGGRESSIVENESS";
      }

      string num=tag.substr(baseTag.size(),tag.size()-baseTag.size());
      int index=ExtendedString(num).Ival();
      vectorParams[baseTag].push_back(*it_tpL);
      isVectorParam=true;
    }
    if (tag == "STEPPER")
    {
      stepperGiven=true;
      if (it_tpL->isNumeric())
      {
        int iType = it_tpL->iVal();
        if (iType == 1)
        {
          stepperList.set("Continuation Method", "Arc Length");
          stepperList.set("Skip df/dp", false);
        }
        else
        {
          stepperList.set("Continuation Method", "Natural");
          stepperList.set("Skip df/dp", true);
        }
      }
      else
      {
        ExtendedString p(it_tpL->sVal());
        p.toUpper();
        if (p.substr(0,3) == "ARC")
        {
          stepperList.set("Continuation Method", "Arc Length");
          stepperList.set("Skip df/dp", false);
        }
        else if (p.substr(0,3) == "NAT")
        {
          stepperList.set("Continuation Method", "Natural");
          stepperList.set("Skip df/dp", true);
        }
        else
        {
          string message = "Unknown specification in .options for 'stepper': ";
          message += it_tpL->sVal();
          message += ".  Legal choices are ARCLENGTH or NATURAL, which may be abbreviated to three characters.";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
        }
      }
    }
    else if (tag == "PREDICTOR")
    {
      predictorGiven=true;
      if (it_tpL->isNumeric())
      {
        int iType = it_tpL->iVal();
        if (iType == 1)
        {
          predictorList.set("Method", "Tangent");
        }
        else if (iType == 2)
        {
          predictorList.set("Method", "Secant");
        }
        else if (iType == 3)
        {
          predictorList.set("Method", "Random");
        }
        else
        {
         predictorList.set("Method", "Constant");
        }
      }
      else
      {
        ExtendedString p(it_tpL->sVal());
        p.toUpper();
        if (p.substr(0,3) == "TAN")
        {
          predictorList.set("Method", "Tangent");
        }
        else if (p.substr(0,3) == "SEC")
        {
          predictorList.set("Method", "Secant");
        }
        else if (p.substr(0,3) == "RAN")
        {
          predictorList.set("Method", "Random");
        }
        else if (p.substr(0,3) == "CON")
        {
         predictorList.set("Method", "Constant");
        }
        else
        {
          string message = "Unknown specification in .options for 'predictor': ";
          message += it_tpL->sVal();
          message += ".  Legal choices are TANGENT, SECANT, RANDOM, CONSTANT, which may be abbreviated to three characters.";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
        }
      }
    }
    // if we're the very first instance of one of these vector params,
    // then set our stepperList (this is the old behavior from the
    // "alternate" vector handling
    else if (tag == "CONPARAM1") // continuation parameter
    {
      stepperList.set("Continuation Parameter", it_tpL->sVal());
    }
    else if (tag == "INITIALVALUE1")
    {
      stepperList.set("Initial Value", it_tpL->dVal());
    }
    else if (tag == "MINVALUE1")
    {
      gstepping_min_value_ = it_tpL->dVal();
      stepperList.set("Min Value", it_tpL->dVal());
    }
    else if (tag == "RESIDUALCONDUCTANCE")
    {
      gstepping_minimum_conductance_ = it_tpL->dVal();
      if (gstepping_minimum_conductance_ > 0)
      {
        string message = "You have specified a non-zero value for the GMIN Stepping residual conductance (RESIDUALCONDUCTANCE= ";
        message += it_tpL->sVal();
        message += ").\nThis option should never be used unless absolutely necessary to obtain an initial condition for transient runs with ill-posed DC operating points.  The operating point obtained by GMIN Stepping will not be a valid steady state condition for the circuit as defined in the netlist, but might possibly produce a reasonable initial condition for transient runs.\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, message);
      }
    }
    else if (tag == "INITIALSTEPSIZE1")
    {
      stepSizeList.set("Initial Step Size", it_tpL->dVal());
    }
    else if (tag == "MINSTEPSIZE1")
    {
      stepSizeList.set("Min Step Size", it_tpL->dVal());
    }
    else if (tag == "MAXSTEPSIZE1")
    {
      stepSizeList.set("Max Step Size", it_tpL->dVal());
    }
    else if (tag == "AGGRESSIVENESS1")
    {
      stepSizeList.set("Aggressiveness", it_tpL->dVal());
    }
    else if (tag == "MAXVALUE1")
    {
      stepperList.set("Max Value", it_tpL->dVal());
    }
    else if (tag == "BIFPARAM") // bifurcation parameter
    {
      stepperList.set("Bifurcation Parameter", it_tpL->sVal());
    }
    else if (tag == "MAXSTEPS")
    {
      stepperList.set("Max Steps", it_tpL->iVal());
    }
    else if (tag == "MAXNLITERS")
    {
      stepperList.set("Max Nonlinear Iterations", it_tpL->iVal());
    }
    else if (tag == "STEPCONTROL")
    {
      if (it_tpL->isNumeric())
      {
        int iType = it_tpL->iVal();
        if (iType == 1)
          stepSizeList.set("Method", "Adaptive");
        else
          stepSizeList.set("Method", "Constant");
      }
      else
      {
        ExtendedString p(it_tpL->sVal());
        p.toUpper();
        if (p.substr(0,3) == "ADA")
        {
          stepSizeList.set("Method", "Adaptive");
        }
        else if (p.substr(0,3) == "CON")
        {
          stepSizeList.set("Method", "Constant");
        }
        else
        {
          string message = "Unknown specification in .options for 'stepcontrol': ";
          message += it_tpL->sVal();
          message += ".  Legal choices are ADAPTIVE or CONSTANT, which may be abbreviated to three characters.";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
        }
      }
    }
    else if (tag == "POWERNODE") // continuation parameter
    {
      stepperList.set("Power Node", it_tpL->sVal());
    }
    else if (tag == "VOLTAGELIST")
    {
      ExtendedString p(it_tpL->sVal());
      p.toUpper();
      if (p.substr(0,3) == "DOF") // DOFS - degrees of freedom
      {
        voltageListType_ = VLT_DOFS;
      }
      else if (p.substr(0,3) == "NOD") // NODES - voltage nodes
      {
        voltageListType_ = VLT_Node;
      }
      else
      {
        string message = "Unknown specification in .options loca for 'voltagelist': ";
        message += it_tpL->sVal();
        message += ".  Legal choices are DOFS or NODES, which may be abbreviated to three characters.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
      }
    }
    else if(tag == "VOLTAGESCALEFACTOR")
    {
      voltageScaleFactor_ = it_tpL->dVal();
    }
    // Start of the new parameter section
    else if (string(tag,0,10) == "PARAMLIST")
    {
      // don't know what to do yet.
      cout << "tag = " << tag << endl;
    }
    else
    {
      // if "isVectorParam" we've already handled this at the beginning
      // of the loop.  Otherwise it's an unrecognized parameter.
      if (!isVectorParam)
      {
        string tmp =
          "  Warning: " + tag +
          " is not a recognized loca option.\n";
        N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_WARNING_0,  tmp);
      }
    }
  }

  // Insure that the correct defaults are set.  Sometimes the LOCA
  // defaults change in the LOCA library from release to release.  However
  // Xyce's defaults should not necessarily change in that case.
  if (!stepperGiven)
  {
    stepperList.set("Continuation Method", "Natural");
    stepperList.set("Skip df/dp", true);
  }

  if (!predictorGiven)
  {
    predictorList.set("Method", "Tangent");
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getVectorParam
// Purpose       : Obtain a parameter specified in the option line in vector syntax
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/02/06
//-----------------------------------------------------------------------------
bool ParameterSet::getVectorParam (const string & tag, int index, double & value)
{
  if (vectorParams.find(tag) != vectorParams.end() &&
      vectorParams[tag].size() > index)
  {
    value = vectorParams[tag][index].dVal();
    return true;
  }
  else
  {
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getVectorParam
// Purpose       : Obtain a parameter specified in the option line in vector syntax
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/02/06
//-----------------------------------------------------------------------------
bool ParameterSet::getVectorParam (const string & tag, int index, string & value)
{
  if (vectorParams.find(tag) != vectorParams.end() &&
      vectorParams[tag].size() > index)
  {
    value = vectorParams[tag][index].sVal();
    return true;
  }
  else
  {
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getVectorParamSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int ParameterSet::getVectorParamSize(const string& tag)
{
  if (vectorParams.find(tag) != vectorParams.end())
  {
    return static_cast<int>(vectorParams[tag].size());
  }
  else
  {
    std::string msg = "N_NLS_NOX::ParameterSet::getVectorParam - ";
    msg += "the parameter \"";
    msg += tag;
    msg += "\" is required for parameter continuation!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  return -1;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getStatusTestReturnCode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int ParameterSet::getStatusTestReturnCode() const
{
  // Get the main Xyce Test
  Teuchos::RefCountPtr<N_NLS_NOX::XyceTests> testPtr =
    Teuchos::rcp_dynamic_cast<N_NLS_NOX::XyceTests>(tests_[1]);
  if (Teuchos::is_null(testPtr))
  {
    const string message = "N_NLS_NOX::getStatusTestReturnCode - Dynamic cast on Xyce Tests check failed.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return testPtr->getXyceReturnCode();
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::setStatusTestReturnCode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void ParameterSet::setStatusTestReturnCodes
  (const N_NLS_ReturnCodes & retCodesTmp)
{
  // Get the main Xyce Test
  Teuchos::RefCountPtr<N_NLS_NOX::XyceTests> testPtr =
    Teuchos::rcp_dynamic_cast<N_NLS_NOX::XyceTests>(tests_[1]);
  if (Teuchos::is_null(testPtr))
  {
    const string message = "N_NLS_NOX::setStatusTestReturnCode - Dynamic cast on Xyce Tests check failed.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return testPtr->setReturnCodes (retCodesTmp);
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getMaxNormF
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
double ParameterSet::getMaxNormF() const
{
  // Get the main Xyce Test
  Teuchos::RefCountPtr<N_NLS_NOX::XyceTests> testPtr =
    Teuchos::rcp_dynamic_cast<N_NLS_NOX::XyceTests>(tests_[1]);
  if (Teuchos::is_null(testPtr))
  {
    const string message = "N_NLS_NOX::getMaxNormF - Dynamic cast on Xyce Tests check failed.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return testPtr->getMaxNormF();
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getMaxNormFindex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int ParameterSet::getMaxNormFindex () const
{
  // Get the main Xyce Test
  Teuchos::RefCountPtr<N_NLS_NOX::XyceTests> testPtr =
    Teuchos::rcp_dynamic_cast<N_NLS_NOX::XyceTests>(tests_[1]);
  if (Teuchos::is_null(testPtr))
  {
    const string message = "N_NLS_NOX::getMaxNormFindex - Dynamic cast on Xyce Tests check failed.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return testPtr->getMaxNormFindex ();
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::getNoxSolverType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int ParameterSet::getNoxSolverType() const
{
  return noxSolver;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::createAugmentLinearSystem
// Purpose       : creates an AugmentLinSys strategy object.
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 1416
// Creation Date : 03/08/06
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>
ParameterSet::createAugmentLinearSystem(N_LAS_System* ls) const
{
  Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als;

  if (noxSolver == 9)
  {
    Teuchos::RefCountPtr<Epetra_MapColoring> color_map =
      Teuchos::rcp(ls->builder().createSolnColoring());

    if (voltageScaleFactor_ == 1.0)
    {
      als = Teuchos::rcp( new N_NLS_NOX::
			  AugmentLinSysPseudoTransient(color_map,
						       ls->getRHSVector()) );
    }
    else
    {
      als = Teuchos::rcp( new N_NLS_NOX::
			  AugmentLinSysPseudoTransient(color_map,
						       ls->getRHSVector(),
						       true,
						       voltageScaleFactor_) );
    }

  }
  else if(noxSolver == 1 || noxSolver == 3)
  {
    if (voltageListType_ == VLT_DOFS)
    {
      Teuchos::RefCountPtr<Epetra_MapColoring> color_map =
        Teuchos::rcp(ls->builder().createSolnColoring());

      als = Teuchos::rcp( new N_NLS_NOX::GStepping(color_map,
						   ls->getRHSVector(),
						   gstepping_min_value_,
               gstepping_minimum_conductance_) );
    }
    else
    {
      als = Teuchos::rcp( new N_NLS_NOX::
			  GStepping(ls->getQueryUtil()->vnodeGIDVec(),
				    ls->getRHSVector(),
				    gstepping_min_value_,
            gstepping_minimum_conductance_) );
    }
  }
  else
  {
    string message = "N_NLS_NOX::createAugmentLinearSystem - The \'continuation\' ";
    message += "parameter in the .options nox list must be set to PSEUDO or NATURAL for ";
    message += "this function to be called!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return als;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::createAugmentLinearSystem
// Purpose       : creates an AugmentLinSys strategy object.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/08/06
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>
ParameterSet::createAugmentLinearSystem(N_LAS_System* ls, Xyce::NodeNamePairMap & op,
#ifdef Xyce_PARALLEL_MPI
   Xyce::NodeNamePairMap & allNodes, N_PDS_Comm * pdsCommPtr
#else
   Xyce::NodeNamePairMap & allNodes
#endif
   ) const
{
  Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als;

#ifdef Xyce_PARALLEL_MPI
  als = Teuchos::rcp( new N_NLS_NOX::AugmentLinSysOPStart(op, allNodes, pdsCommPtr) );
#else
  als = Teuchos::rcp( new N_NLS_NOX::AugmentLinSysOPStart(op, allNodes) );
#endif

  return als;
}

//-----------------------------------------------------------------------------
// Function      : ParameterSet::createAugmentLinearSystem
// Purpose       : creates an AugmentLinSys strategy object for .IC with
//                 or without gmin stepping
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/08/07
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>
ParameterSet::createAugmentLinearSystem(N_LAS_System* ls, Xyce::NodeNamePairMap & op,
   bool gminStepping) const
{
  Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als;

  Teuchos::RefCountPtr<Epetra_MapColoring> ICcolor_map =
        Teuchos::rcp(ls->builder().createInitialConditionColoring());

  if (gminStepping==false)
  {
    als = Teuchos::rcp( new N_NLS_NOX::AugmentLinSysIC(op, ICcolor_map, ls->getRHSVector() ) );
  }
  else
  {
    if (voltageListType_ == VLT_DOFS)
    {
      Teuchos::RefCountPtr<Epetra_MapColoring> GMINcolor_map =
        Teuchos::rcp(ls->builder().createSolnColoring());

      als = Teuchos::rcp( new N_NLS_NOX::AugmentLinSysIC_Gmin(
                op,
                ICcolor_map,
                GMINcolor_map,
                ls->getRHSVector(),
						    gstepping_min_value_,
                gstepping_minimum_conductance_) );
    }
    else
    {
      als = Teuchos::rcp( new N_NLS_NOX::AugmentLinSysIC_Gmin(
                op,
                ICcolor_map,
                //AugmentLinSysIC_Gmin(ls->getQueryUtil()->vnodeGIDVec()),
                ls->getQueryUtil()->vnodeGIDVec(),
                ls->getRHSVector(),
						    gstepping_min_value_,
                gstepping_minimum_conductance_) );
    }
  }

  return als;
}

