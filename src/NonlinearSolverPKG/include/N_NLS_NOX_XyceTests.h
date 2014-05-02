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
// Filename       : $RCSfile: N_NLS_NOX_XyceTests.h,v $
//
// Purpose        : Particular Status Test Based on the Weighed Norm of the
//                  Update
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
// Revision Number: $Revision: 1.23 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#ifndef Xyce_N_NLS_NOX_XyceTests_h
#define Xyce_N_NLS_NOX_XyceTests_h

// ---------- Standard Includes ----------
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_NLS_ReturnCodes.h>

// ----------   NOX Includes   ----------

#include "NOX_StatusTest_Generic.H"	// base class
#include "NOX_StatusTest_FiniteValue.H"

// ---------- Forward Declarations ----------

class N_LAS_Vector;
class N_LAS_System;
class N_LOA_Loader;

// ---------- Namespace Declarations ----------

// N_NLS namespace is for the Xyce Nonlinear Solver Package
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Class         : N_NLS_NOX::XyceTests
//
// Purpose       :
//
//      NOX Vector Interface for Xyce vectors.
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 2/1/02
//-----------------------------------------------------------------------------

class XyceTests : public NOX::StatusTest::Generic {

public:

  //---------------------------------------------------------------------------
  // Function      : XyceTests (constructor)
  //
  // Purpose       : Constructs a NOX-compatiable status test based
  //                 the weighted norm of the update.
  //
  //---------------------------------------------------------------------------
  XyceTests(bool isTransient, 
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
	    );
  
  //---------------------------------------------------------------------------
  // Function      : Destructor
  //---------------------------------------------------------------------------
  ~XyceTests();

  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  NOX::StatusTest::StatusType 
    checkStatus(const NOX::Solver::Generic& problem,
		NOX::StatusTest::CheckType checkType);
  
  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  NOX::StatusTest::StatusType getStatus() const { return status_; };
  
  //---------------------------------------------------------------------------
  // Purpose       : Get the return code to send to the time stepper.
  //           
  //---------------------------------------------------------------------------
  int getXyceReturnCode() const;
  
  //---------------------------------------------------------------------------
  // Purpose       : Output formatted description of stopping test to
  //                 output stream.
  //---------------------------------------------------------------------------
  std::ostream& print(std::ostream& stream, int indent = 0) const;
  
  //---------------------------------------------------------------------------
  // Purpose       : Set a specific set of return codes to be used.
  //---------------------------------------------------------------------------
  void setReturnCodes (const N_NLS_ReturnCodes & retCodesTmp);

  double getMaxNormF() const;

  int getMaxNormFindex () const;
  
protected:
  
  NOX::StatusTest::StatusType status_;
  int returnTest_;
  bool isTransient_;
  int niters_;
  
  N_NLS_ReturnCodes retCodes_;
  
  //Test #0 
  //****************************************
  NOX::StatusTest::FiniteValue finiteTest_;
    
  //Test #1 
  //****************************************
  int maxNormFindex_;
  double maxNormF_;
  double requestedMaxNormF_;
  double requestedMachPrecTol_;

  // Test #2
  //****************************************
  N_LAS_Vector** oldTimeStepVectorPtrPtr_;
  N_LAS_Vector* weightsVectorPtr_;
  N_LAS_Vector* updateVectorPtr_;
  N_LAS_Vector* tmpVectorPtr_;
  const double epsilon_a_;
  const double epsilon_r_;
  const double tol_;
  double weightedUpdate_;

  // Test #3
  //****************************************
  // Maximum number of nonlinear iterations allowed
  int maxIters_;

  // ||F(x_current)|| / ||F(x_previous)||
  const double requestedConvRate_;
  double currentConvRate_;

  //  ||F(x)|| / ||F(x_init)||
  const double requestedRelativeConvRate_;
  double currentRelativeConvRate_;

  // Initial norm of the RHS used to calculate ratio
  double normResidualInit_;
  
  // test #4
  double smallUpdateTol_;

  // Test #6
  const double maxConvRate_;
    
  // Test #7
  int lastIteration_;
  int badStepCount_;
  const int maxBadSteps_;
  double minConvRate_;
  const double stagnationTol_;
    
  // Xyce return Code
  int xyceReturnCode_;

  // For Device Specific Convergence cirteria
  int checkDeviceConvergence_;
  N_LOA_Loader* loaderPtr_;
#ifdef Xyce_NLS_MASKED_WRMS_NORMS
  bool deviceMaskFlag_;
  N_LAS_Vector* weightMaskVectorPtr_;
#endif
  bool allDevicesConverged_;
  bool innerDevicesConverged_;

}; // class SharedSystem

inline void XyceTests::setReturnCodes 
  (const N_NLS_ReturnCodes & retCodesTmp)
{
  retCodes_ = retCodesTmp;
}


} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_SharedSystem_h

