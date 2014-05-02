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
// Filename       : $RCSfile: N_NLS_ReturnCodes.C,v $
//
// Purpose        : This is a simple container class for different convergence
//                  status "return codes".
//
// Special Notes  : This point of having a class like this is to make the
//                  code more flexible in terms of what is considered to be
//                  a successful Newton solve.
//
//                  In the class N_NLS_DampedNewton, the function
//                  "converged_" performs a series of tests to determine
//                  convergence status, and depending on the results of
//                  these tests, it returns a value.  If the value is
//                  positive, then the solve is considered to be converged,
//                  and if it is <=0, it isn't.  
//
//                  This convention is used both inside the damped newton
//                  class, and outside.  It is used inside when the solver
//                  is trying to determine for itself whether or not to
//                  continue with the solve, or exit.  It is used outside
//                  by the time integrator, and also by continuation loops
//                  in the 2-level solver, in determining if the step was
//                  successful.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/27/03
//
// Revision Information:
// ----------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author $
//-------------------------------------------------------------------------

#include <iostream>

#include <N_NLS_ReturnCodes.h>
#include <N_ERH_ErrorMgr.h>

//-----------------------------------------------------------------------------
// Function      : N_NLS_ReturnCodes::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/02/03
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const N_NLS_ReturnCodes & rc)
{
  os << "\n\n-----------------------------------------" << std::endl
     << "\tNonlinear Solver Return Codes:\n"
     << "\t\tnormTooSmall      = " << rc.normTooSmall << "\n"
     << "\t\tnormalConvergence = " << rc.normalConvergence << "\n"
     << "\t\tnearConvergence   = " << rc.nearConvergence << "\n"
     << "\t\tsmallUpdate       = " << rc.smallUpdate << "\n"
     << "\t\tnanFail           = " << rc.nanFail << "\n"
     << "\t\ttooManySteps      = " << rc.tooManySteps << "\n"
     << "\t\ttooManyTranSteps  = " << rc.tooManyTranSteps << "\n"
     << "\t\tupdateTooBig      = " << rc.updateTooBig << "\n"
     << "\t\tstalled           = " << rc.stalled << "\n"
     << "\t\twrmsExactZero     = " << rc.wrmsExactZero << "\n"
     << Xyce::section_divider << std::endl
     << std::endl;

  return os;
}
