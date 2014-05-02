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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_XyceTests.C,v $
//
// Purpose        : Status test.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 04/15/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.44 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include "N_NLS_NOX_XyceTests.h"
#include "N_NLS_NOX_Vector.h"
#include "N_LAS_Vector.h"
#include "N_LOA_Loader.h"
#include "NOX.H"
#include "NOX_Solver_LineSearchBased.H"

// ----------   Namespaces   ----------

using namespace N_NLS_NOX;

// ----------   Code   ----------

XyceTests::XyceTests(bool isTransient,
		     double normF,
		     double machPrec,
		     N_LAS_Vector** currSolVectorPtrPtr,
		     double epsilon_a, 
		     double epsilon_r, 
		     double tol,
		     int maxIters,
		     double convRate,
		     double relConvRate,
		     double maxConvRate,
		     double stagnationTol,
		     int maxBadSteps,
		     int checkDeviceConvergence,
		     double smallUpdateTol,
		     N_LOA_Loader* loader
#ifdef Xyce_NLS_MASKED_WRMS_NORMS
         , bool nonTrivialDeviceMaskFlag,
         N_LAS_Vector * maskVectorPtr
#endif		     
		     ) :
  
  status_(NOX::StatusTest::Unconverged),
  returnTest_(0),
  isTransient_(isTransient),
  niters_(-1),
  maxNormFindex_(-1),
  maxNormF_(0.0),
  requestedMaxNormF_(normF),
  requestedMachPrecTol_(machPrec),
  oldTimeStepVectorPtrPtr_(currSolVectorPtrPtr),
  weightsVectorPtr_(0),
  updateVectorPtr_(0),
  tmpVectorPtr_(0),
  epsilon_a_(epsilon_a),
  epsilon_r_(epsilon_r),
  tol_(tol),
  weightedUpdate_(0.0),
  maxIters_(maxIters),
  requestedConvRate_(convRate),
  currentConvRate_(1.0),
  requestedRelativeConvRate_(relConvRate),
  currentRelativeConvRate_(1.0),
  normResidualInit_(1.0),
  maxConvRate_(maxConvRate),
  lastIteration_(-1),
  badStepCount_(0),
  minConvRate_(1.0),
  stagnationTol_(stagnationTol),
  maxBadSteps_(maxBadSteps),
  xyceReturnCode_(0),
  smallUpdateTol_(smallUpdateTol),
  checkDeviceConvergence_(checkDeviceConvergence),
  loaderPtr_(loader),
  allDevicesConverged_(false),
  innerDevicesConverged_(false)
#ifdef Xyce_NLS_MASKED_WRMS_NORMS
  , deviceMaskFlag_( false ), 
  weightMaskVectorPtr_( maskVectorPtr)
#endif
{

}

XyceTests::~XyceTests()
{
  if( weightsVectorPtr_ != 0 )
  {
    delete weightsVectorPtr_;
    delete updateVectorPtr_;
    delete tmpVectorPtr_;
  }
}


NOX::StatusTest::StatusType XyceTests::
checkStatus(const NOX::Solver::Generic& problem,
	    NOX::StatusTest::CheckType checkType)
{
  status_ = NOX::StatusTest::Unconverged;
  xyceReturnCode_ = 0;
  niters_ = problem.getNumIterations();

  returnTest_ = 0;




  /*
  problem.getSolutionGroup();
  const NOX::Abstract::Vector& v = problem.getSolutionGroup().getX();
  cout << "ROGER - norm = " << v.norm() << endl;
  const N_NLS_NOX::Vector* testPtr = 0;
  testPtr = dynamic_cast<const N_NLS_NOX::Vector*>(&v);
  cout << "ROGER - ptr = " << testPtr << endl;
  
  if (testPtr == 0) {
    const N_NLS_NOX::Vector* testPtr = 0;
    testPtr = dynamic_cast<const N_NLS_NOX::Vector*>(&v);
  }
    //const N_NLS_NOX::Vector& v_nox = dynamic_cast<const N_NLS_NOX::Vector&>
    //(problem.getSolutionGroup().getX());

    //const N_LAS_Vector& x__ = (dynamic_cast<const N_NLS_NOX::Vector&>
    //(problem.getSolutionGroup().getX())).getNativeVectorRef();
  exit(0);
  */


  // Get the current and previous solutions
  const N_LAS_Vector& x = (dynamic_cast<const N_NLS_NOX::Vector&>
    (problem.getSolutionGroup().getX())).getNativeVectorRef();
  const N_LAS_Vector& oldX = (dynamic_cast<const N_NLS_NOX::Vector&>
    (problem.getPreviousSolutionGroup().getX())).getNativeVectorRef();

  // Test0 - NaN/Inf checker
  NOX::StatusTest::StatusType check = finiteTest_.checkStatus(problem,
							      checkType);
  if (check == NOX::StatusTest::Failed) 
  {
    status_ = check;
    returnTest_ = 0;
    xyceReturnCode_ = retCodes_.nanFail; // default: -6
    return status_;
  }

  // This test is for 2-level solves only.  
  // If the inner solve failed, then the whole thing needs to fail.
  // If 2-level is not being used, this test doesn't do anything.
  innerDevicesConverged_ = loaderPtr_->innerDevsConverged();
  if (!innerDevicesConverged_)
  {
    status_ = NOX::StatusTest::Failed;
    returnTest_ = 9;
    xyceReturnCode_ = retCodes_.innerSolveFailed; // default: -5
    return status_;
  }

  // Test 8 - Devices need to satisfy their own convergence criteria
  if (checkDeviceConvergence_) 
  {
    allDevicesConverged_ = loaderPtr_->allDevsConverged();
    if (!allDevicesConverged_ && (niters_ < maxIters_)) 
    {
      status_ = NOX::StatusTest::Unconverged;
      returnTest_ = 8;
      xyceReturnCode_ = 0;
      return status_;
    }
  }

  //Test 1 - Make sure the residual isn't too small, hardwired tolerances
  maxNormF_ = problem.getSolutionGroup().getF()
    .norm(NOX::Abstract::Vector::MaxNorm);

  const N_LAS_Vector& F = (dynamic_cast<const N_NLS_NOX::Vector&>
    (problem.getSolutionGroup().getF())).getNativeVectorRef();

  std::vector<int> index(1, -1);
  F.infNormIndex( &index[0] );
  maxNormFindex_ = index[0];

  if ((maxNormF_ < requestedMaxNormF_) && 
      (maxNormF_ < requestedMachPrecTol_)) 
  {
    status_ = NOX::StatusTest::Converged;
    returnTest_ = 1;
    xyceReturnCode_ = retCodes_.normTooSmall; // default: 1
    return status_;
  }

  // Test 2 - Normal convergence based on rhs residual (2a) and 
  // update norm (2b).

  // Copy into local reference
  N_LAS_Vector& oldTimeStepX = **oldTimeStepVectorPtrPtr_;

  // Allocate space if necessary
  if (weightsVectorPtr_ == 0) 
  {
    weightsVectorPtr_ = new N_LAS_Vector(x);
    // when we create weightsVectorPtr_ from the latest solution, there
    // is a chance that one of the values will be zero.  If this isn't 
    // the DC op step or the first iteration of a time step, then
    // we'll end up dividing by zero in the wMaxNorm function below.
    // So to be safe we'll just add epsilon_a_ on the when we create 
    // this vector.
    for (int i=0; i< x.localLength() ; ++i ) 
    {
      (*weightsVectorPtr_)[i] += epsilon_a_;
    }    
    updateVectorPtr_ = new N_LAS_Vector(x);
    tmpVectorPtr_ = new N_LAS_Vector(x);
  }

  // Local references
  N_LAS_Vector& weights = *weightsVectorPtr_;
  N_LAS_Vector& update = *updateVectorPtr_;
  N_LAS_Vector& tmp = *tmpVectorPtr_;

  // Compute local portion of weights vector
  // Weights are recomputed at each nonlinear iteration of a DC Op calc
  // but only at the beginning of a transient nonlinear solve.
  if ((!isTransient_) || (niters_ == 0)) 
  {
    int length = x.localLength();
    for (int i = 0; i < length; ++i ) 
    {
      //update[i] = x[i] - oldX[i];
      weights[i] =
       epsilon_r_ * Xycemax(fabs(x[i]), fabs(oldTimeStepX[i])) + epsilon_a_;
#ifdef Xyce_NLS_MASKED_WRMS_NORMS
      if( deviceMaskFlag_ && ((*weightMaskVectorPtr_)[i] == 0.0) )
      {
        weights[i] = N_UTL_MachineDependentParams::MachineBig();
      }
#endif    
    }
  }
  


  if (niters_ < 1) 
  {
    weightedUpdate_ = 1.0;
  }
  else 
  {
    // Next compute the update
    update.update(1.0, x, -1.0, oldX, 0.0);

    // Compute final result
#ifdef Xyce_SPICE_NORMS
    update.wMaxNorm(weights,tmp,&weightedUpdate_);
#else
    update.wRMSNorm(weights,&weightedUpdate_);
#endif
    
    // RPP: If a line search is being used, we must account for any 
    // damping of the step length.  Otherwise delta X could be small due 
    // the line search and not due to being close to a solution.
    const NOX::Solver::LineSearchBased* test = 0;
    test = dynamic_cast<const NOX::Solver::LineSearchBased*>(&problem);
    if (test != 0) 
    {
      weightedUpdate_ = weightedUpdate_/(test->getStepSize());
    }

    // RPP: 11/11/2003 - Fix for Bug 354 - Xyce fails in DC Op calc 
    // but proceeds to transient.
    // Check to see if WRMS is exactly zero.  If so the linear solver 
    // has failed.  Return a failure. This should be commented out for 
    // Broyden runs.
    // RPP: 02/18/2004 - Causing premature failures if deltaX is ~1e-16.
    // Adding a check to look at norm of dx
    // RPP: 08/09/2004 - Still goes into transient from failed DC Op on
    // freebsd platforms for Down_8-bit_03.cir (bug 354).  Changed return
    // code from +4 to -4.  This forces TIA to assume failure.
    if ((weightedUpdate_ == 0.0) && (niters_ > 0))
    {
      if (!(problem.getPreviousSolutionGroup().isNewton()))
      {
        NOX::Abstract::Group & tmpGrp = 
        (const_cast<NOX::Abstract::Group&>(problem.getPreviousSolutionGroup()));
        Teuchos::ParameterList tmpParams;
        tmpGrp.computeNewton (tmpParams);
      }

      const N_LAS_Vector& dx = (dynamic_cast<const N_NLS_NOX::Vector&>
	    (problem.getPreviousSolutionGroup().getNewton())).getNativeVectorRef();

      double tmp = 0.0;
      dx.lpNorm(2, &tmp );
      if (tmp == 0.0) 
      {
        status_ = NOX::StatusTest::Failed;
        returnTest_ = 4;
        xyceReturnCode_ = retCodes_.wrmsExactZero; // default: -4
        return status_;
      }
    }

    if ((weightedUpdate_ < tol_) &&(maxNormF_ < requestedMaxNormF_)) 
    {
      status_ = NOX::StatusTest::Converged;
      returnTest_ = 2;
      xyceReturnCode_ = retCodes_.normalConvergence; // default: 2
      return status_;
    }
  }
  
  // Test 3 - Near Convergence - Hit max iterations but residual 
  // and convergence rate indicate we may be near a converged solution.  
  // Therefore, let the time stepper decide whether or  not the step is ok.
  // Transient mode ONLY!
  // NOTE: Convergence rates are based on the 2-Norm, not max norm!
  if (niters_ > 0) 
  {
    // ||F(x_current)|| / ||F(x_previous)||
    currentConvRate_ = (problem.getSolutionGroup().getNormF()) / 
      (problem.getPreviousSolutionGroup().getNormF());

    //  ||F(x)|| / ||F(x_init)||
    currentRelativeConvRate_ = (problem.getSolutionGroup().getNormF()) / 
    normResidualInit_;
  }    
  else 
  {
    currentConvRate_ = 1.0;
    currentRelativeConvRate_ = 1.0;
  }

  if (isTransient_) 
  {
    if (niters_ == 0) 
    {
      normResidualInit_ = problem.getSolutionGroup().getNormF();
    }
 
    // Test only if we hit the max number of iterations
    if (niters_ >= maxIters_) 
    {
      if ((currentConvRate_ <= requestedConvRate_) && 
          (currentRelativeConvRate_ <= requestedRelativeConvRate_)) 
      {
        status_ = NOX::StatusTest::Converged;
        xyceReturnCode_ = retCodes_.nearConvergence; // default: 3

        if (xyceReturnCode_ < 0)
          status_ = NOX::StatusTest::Failed;
        else
          status_ = NOX::StatusTest::Converged;
      }
      else
      {
        status_ = NOX::StatusTest::Failed;
        xyceReturnCode_ = retCodes_.tooManySteps; // default: -1
      }
      
      returnTest_ = 3;
      return status_;
    }
  } // end test 3

  // Test 4 - Update is too small
  if ((niters_ > 0) && (weightedUpdate_ < smallUpdateTol_) && (niters_ < maxIters_)) 
  {
    if (isTransient_) // Let the time integrator determine convergence (+4)
    {
      xyceReturnCode_ = retCodes_.smallUpdate; // default: 4
      status_ = NOX::StatusTest::Failed;
    }
    else // Steady state should always return a hard failure -4 
    {
      xyceReturnCode_ = 0; // neither pass nor fail.
      status_ = NOX::StatusTest::Unconverged;

      //xyceReturnCode_ = retCodes_.wrmsExactZero; // default: -4
      //status_ = NOX::StatusTest::Failed;
    }

    returnTest_ = 4;

    return status_;
  }

  // Test 5 - Max nonlinear iterations (if transient, this will be checked 
  // in the NearConvergence test (#3)
  if ((!isTransient_) && (niters_ >= maxIters_)) 
  {
    status_ = NOX::StatusTest::Failed;
    returnTest_ = 5;
    xyceReturnCode_ = retCodes_.tooManySteps; // default: -1
    return status_;
  }

  // Test 6 - update is too big
  if (currentConvRate_ > maxConvRate_) 
  {
    status_ = NOX::StatusTest::Failed;
    returnTest_ = 6;
    xyceReturnCode_ = retCodes_.updateTooBig; // default: -2
    return status_;
  }

  // Test 7 - Stall in the convergence rate. Transient mode ONLY! 
  if (isTransient_) 
  {  
    // First time through we don't do anything but reset the counters
    if (niters_ == 0) 
    {
      badStepCount_ = 0;
      lastIteration_ = 0;
      //minConvRate = 1.0;  // Don't reset this.  Xyce solver never does.
    } 

    // Make sure we have not already counted the last nonlinear iteration.
    // This protects against multiple calls to checkStatus() in between 
    // nonlinear iterations.
    bool isCounted = false;
    if (niters_ == lastIteration_) 
    {
      isCounted = true;
    }
    else
    {
      lastIteration_ = niters_;
    }
    
    // Set counter appropriately
    if (!isCounted) 
    {
      if (fabs(currentConvRate_ - 1.0) <= stagnationTol_) 
      {
        if ((badStepCount_ == 0) || (currentConvRate_ < minConvRate_)) 
        {
          minConvRate_ = currentConvRate_;
        }
        ++badStepCount_ ;
      }
      else
      {
        badStepCount_ = 0;
      }
    }

    if (badStepCount_ >= maxBadSteps_) 
    {
      if ((currentRelativeConvRate_ <= 0.9) && (minConvRate_ <= 1.0)) 
      {
        status_ = NOX::StatusTest::Converged;
        returnTest_ = 7;
        xyceReturnCode_ = retCodes_.nearConvergence;   // default: 3
                // note - I'm not sure if this is 
               // really a near convergece test - but 3 is the code for it...

        if (xyceReturnCode_ < 0)
          status_ = NOX::StatusTest::Failed;
        else
          status_ = NOX::StatusTest::Converged;
      }
      else 
      {
        status_ = NOX::StatusTest::Failed;
        returnTest_ = 7;
        xyceReturnCode_ = retCodes_.stalled; // default: -3
      }
    }
  }

  return status_;
}

std::ostream& XyceTests::print(std::ostream& stream, int indent) const
{
  // precision
  int p = 5;

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << status_ << "by Test #" << returnTest_ << "\n";

  indent += 4;

  //for (int j = 0; j < indent; ++j )
  // stream << ' ';
  finiteTest_.print(stream, indent);

  if (checkDeviceConvergence_) {
    for (int j = 0; j < indent; ++j )
      stream << ' ';
    stream << "8. Devices are Converged: ";
    if (allDevicesConverged_)
      stream << "true" << "\n";
    else
      stream << "false" << "\n";
  }

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "1. Inf-Norm F too small" << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Machine Precision: " << NOX::Utils::sciformat(maxNormF_, p)
	 << " < " << NOX::Utils::sciformat(requestedMachPrecTol_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Requested Tolerance: " << NOX::Utils::sciformat(maxNormF_, p)
	 << " < " << NOX::Utils::sciformat(requestedMaxNormF_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "2. Normal Convergence" << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Inf-Norm F: " << NOX::Utils::sciformat(maxNormF_, p)
	 << " < " << NOX::Utils::sciformat(requestedMaxNormF_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Weighted Update: " << NOX::Utils::sciformat(weightedUpdate_, p)
	 << " < " << NOX::Utils::sciformat(tol_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "3. Near Convergence" << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Max Iters: " << niters_
	 << " < " << maxIters_ << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Convergence Rate: " 
	 << NOX::Utils::sciformat(currentConvRate_, p)
	 << " < " << NOX::Utils::sciformat(requestedConvRate_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Relative Convergence Rate: " 
	 << NOX::Utils::sciformat(currentRelativeConvRate_, p)
	 << " < " << NOX::Utils::sciformat(requestedRelativeConvRate_, p) 
	 << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "4. Small Weighted Update: " 
	 << NOX::Utils::sciformat(weightedUpdate_, p)
	 << " < " << NOX::Utils::sciformat(smallUpdateTol_, p) << "\n";

  if (!isTransient_) {

    for (int j = 0; j < indent; ++j )
      stream << ' ';
    stream << "5. Maximum Iterations: " 
	   << niters_
	   << " < " << maxIters_ << "\n";
    
    for (int j = 0; j < indent; ++j )
      stream << ' ';
    stream << "6. Large Conv Rate: " 
	   << NOX::Utils::sciformat(currentConvRate_, p)
	   << " < " << NOX::Utils::sciformat(maxConvRate_, p) << "\n";
  }

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "7. Stagnation " << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Bad Step Count: " 
	 << badStepCount_ << " < " << maxBadSteps_ << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Stagnation Tolerance: " 
	 << NOX::Utils::sciformat(fabs(currentConvRate_ - 1.0), p)
	 << " < " << NOX::Utils::sciformat(stagnationTol_, p) << "\n";

  stream << std::endl;
  return stream;
}

int XyceTests::getXyceReturnCode() const
{
  return xyceReturnCode_;
}

double XyceTests::getMaxNormF() const
{
  return maxNormF_;
}

int XyceTests::getMaxNormFindex() const
{
  return maxNormFindex_;
}

