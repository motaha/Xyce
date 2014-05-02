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
// Filename       : $RCSfile: N_NLS_NOX_WeightedUpdateTest.h,v $
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
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_WeightedUpdateTest_h
#define Xyce_N_NLS_NOX_WeightedUpdateTest_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

// ----------   NOX Includes   ----------

#include "NOX_StatusTest_Generic.H"	// base class

// ---------- Forward Declarations ----------

class N_LAS_Vector;
class N_LAS_System;

// ---------- Namespace Declarations ----------

// N_NLS namespace is for the Xyce Nonlinear Solver Package
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Class         : N_NLS_NOX::WeightedUpdateTest
//
// Purpose       :
//
//      NOX Vector Interface for Xyce vectors.
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 2/1/02
//-----------------------------------------------------------------------------

class WeightedUpdateTest : public NOX::StatusTest::Generic {

public:

  //---------------------------------------------------------------------------
  // Function      : WeightedUpdateTest (constructor)
  //
  // Purpose       : Constructs a NOX-compatiable status test based
  //                 the weighted norm of the update.
  //
  //---------------------------------------------------------------------------
  WeightedUpdateTest(N_LAS_Vector** currSolVectorPtrPtr, double epsilon_a, 
		     double epsilon_r, double tol, bool isTransient);

  //---------------------------------------------------------------------------
  // Function      : Destructor
  //---------------------------------------------------------------------------
  ~WeightedUpdateTest();

  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  NOX::StatusTest::StatusType checkStatus(const NOX::Solver::Generic& problem);

  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  NOX::StatusTest::StatusType getStatus() const { return status_; };

  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  double getWeightedUpdate() const { return result_; };

  //---------------------------------------------------------------------------
  // Purpose       : Output formatted description of stopping test to
  //                 output stream.
  //---------------------------------------------------------------------------
  std::ostream& print(std::ostream& stream, int indent = 0) const;

private:

  double result_;

  NOX::StatusTest::StatusType status_;

  // Pointer-pointer to solution from last time step
  N_LAS_Vector** oldTimeStepVectorPtrPtr_;
  N_LAS_Vector* weightsVectorPtr_;
  N_LAS_Vector* updateVectorPtr_;
  N_LAS_Vector* tmpVectorPtr_;

  const double epsilon_a_;
  const double epsilon_r_;
  const double tol_;

  bool isTransient_;

}; // class SharedSystem
} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_SharedSystem_h

