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
// Filename       : $RCSfile: N_NLS_NOX_NearConvergenceTest.h,v $
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
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_NearConvergenceTest_h
#define Xyce_N_NLS_NOX_NearConvergenceTest_h

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
// Creator       : Roger Pawlowski, SNL 9233
//
// Creation Date : 4/2/2003
//-----------------------------------------------------------------------------

class NearConvergenceTest : public NOX::StatusTest::Generic {

public:

  //---------------------------------------------------------------------------
  // Function      : NearConvergenceTest (constructor)
  //
  // Purpose       : Constructs a NOX-compatiable status test based
  //                 the weighted norm of the update.
  //
  //---------------------------------------------------------------------------
  NearConvergenceTest(int maxIters, double convRate, 
		      double relativeConvRate = 1.0);

  //---------------------------------------------------------------------------
  // Function      : Destructor
  //---------------------------------------------------------------------------
  ~NearConvergenceTest();

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
  // Purpose       : Output formatted description of stopping test to
  //                 output stream.
  //---------------------------------------------------------------------------
  std::ostream& print(std::ostream& stream, int indent = 0) const;

private:
  
  // Maximum number of nonlinear iterations allowed
  int maxIters_;
  int nIters_;

  // ||F(x_current)|| / ||F(x_previous)||
  double requestedConvRate_;
  double currentConvRate_;

  //  ||F(x)|| / ||F(x_init)||
  double requestedRelativeConvRate_;
  double currentRelativeConvRate_;

  // Initial norm of the RHS used to calculate ratio
  double normResidualInit_;

  NOX::StatusTest::StatusType status_;

}; // class NearConvergenceTest
} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_NearConvergenceTest

