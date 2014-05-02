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
// Filename       : $RCSfile: N_NLS_NOX_FastTests.h,v $
//
// Purpose        : Status tests based designed for faster, spice-like convergence
//
// Special Notes  :
//
// Creator        : Richard Schiek
//
// Creation Date  : 12/03/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#ifndef Xyce_N_NLS_NOX_FastTests_h
#define Xyce_N_NLS_NOX_FastTests_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_NLS_ReturnCodes.h>

// ----------   NOX Includes   ----------

#include <N_NLS_NOX_XyceTests.h> 	// base class

// ---------- Forward Declarations ----------


// ---------- Namespace Declarations ----------

// N_NLS namespace is for the Xyce Nonlinear Solver Package
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Class         : N_NLS_NOX::FastTests
//
// Purpose       :
//
//      NOX Vector Interface for Xyce vectors.
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 2/1/02
//-----------------------------------------------------------------------------

class FastTests : public XyceTests {

public:

  //---------------------------------------------------------------------------
  // Function      : FastTests (constructor)
  //
  // Purpose       : Constructs a NOX-compatiable status test based
  //                 the weighted norm of the update.
  //
  //---------------------------------------------------------------------------
  FastTests(bool isTransient, 
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
	    N_LOA_Loader* loader,
      std::vector<char> & varTypeVec,
      double voltZeroTol,
      double currZeroTol
#ifdef Xyce_NLS_MASKED_WRMS_NORMS
      , bool nonTrivialDeviceMaskFlag,
      N_LAS_Vector * maskVectorPtr
#endif
	    );
  
  //---------------------------------------------------------------------------
  // Function      : Destructor
  //---------------------------------------------------------------------------
  ~FastTests();

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
  
private:
    
  // vector of variable types
  std::vector<char>  varTypeVec_;
  
  // zero tolerances for voltage and current
  double voltZeroTol_;
  double currZeroTol_;

}; // class FastTests


} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_FastTests_h


