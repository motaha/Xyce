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
// Filename       : $RCSfile: N_NLS_NOX_LargeUpdateTest.C,v $
//
// Purpose        : Status test.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 04/02/2003
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX_LargeUpdateTest.h"
#include "N_NLS_NOX_Vector.h"
#include "N_LAS_Vector.h"
#include "NOX.H"
#include "NOX_Solver_LineSearchBased.H"

// ----------   Namespaces   ----------

using namespace N_NLS_NOX;

// ----------   Code   ----------

LargeUpdateTest::LargeUpdateTest(double maxConvRate) :
  
  maxConvRate_(maxConvRate),
  currentConvRate_(1.0),
  status_(NOX::StatusTest::Unconverged)
{

}

LargeUpdateTest::~LargeUpdateTest()
{

}


NOX::StatusTest::StatusType LargeUpdateTest::checkStatus(const NOX::Solver::Generic& problem)
{
  status_ = NOX::StatusTest::Unconverged;

  if (problem.getNumIterations() < 1) {
    currentConvRate_ = 1.0;
    return status_;
  }

  // ||F(x_current)|| / ||F(x_previous)||
  currentConvRate_ = (problem.getSolutionGroup().getNormF()) / 
    (problem.getPreviousSolutionGroup().getNormF());

    if (currentConvRate_ > maxConvRate_) {
      status_ = NOX::StatusTest::Failed;
    }

  return status_;
}

std::ostream& LargeUpdateTest::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; ++j )
    stream << ' ';

  stream << status_;

  stream << "Update Too Big = " << currentConvRate_ << " > " << maxConvRate_
	 << std::endl;

  return stream;
}

