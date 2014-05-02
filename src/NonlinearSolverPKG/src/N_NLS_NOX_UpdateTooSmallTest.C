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
// Filename       : $RCSfile: N_NLS_NOX_UpdateTooSmallTest.C,v $
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
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX_UpdateTooSmallTest.h"
#include "N_NLS_NOX_WeightedUpdateTest.h"
#include "NOX.H"

// ----------   Namespaces   ----------

using namespace N_NLS_NOX;

// ----------   Code   ----------

UpdateTooSmallTest::UpdateTooSmallTest(const N_NLS_NOX::WeightedUpdateTest& t, 
				       double tol) :
  test_(t),
  tol_(tol)
{
  result_ = 0.0;
  status_ = NOX::StatusTest::Unconverged;
}

UpdateTooSmallTest::~UpdateTooSmallTest()
{
}


NOX::StatusTest::StatusType UpdateTooSmallTest::checkStatus(const NOX::Solver::Generic& problem)
{
  status_ = NOX::StatusTest::Unconverged;
  
  if (problem.getNumIterations() < 1)
    return status_;

  result_ = test_.getWeightedUpdate();

  if (result_ < tol_) 
    status_ = NOX::StatusTest::Converged;

  return status_;
}

std::ostream& UpdateTooSmallTest::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << status_;
  stream << "Weighted Update Too Small = " 
	 << NOX::Utils::sciformat(result_, 3) 
	 << " < " << NOX::Utils::sciformat(tol_, 3);
  stream << std::endl;

  return stream;
}

