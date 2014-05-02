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
// Filename       : $RCSfile: N_NLS_ReturnCodes.h,v $
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
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_ReturnCodes_h
#define Xyce_N_NLS_ReturnCodes_h

#include <iosfwd>

#include <N_UTL_Xyce.h>

//-----------------------------------------------------------------------------
// Class         : N_NLS_ReturnCodes
// Purpose       : Container class for solver success/failure return codes.
//
// Special Notes : Any result that you want the code to perceive as a
//                 "failure" should be <= 0.  Any result you want to code
//                 to perceive as a "success" should be >0.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/27/03
//-----------------------------------------------------------------------------
class N_NLS_ReturnCodes
{
public:
  N_NLS_ReturnCodes ():
    normTooSmall      (1),
    normalConvergence (2),
    nearConvergence   (3),   // (near convergence, but assume success)
    smallUpdate       (4),
    nanFail           (-6),
    tooManySteps      (-1),
    tooManyTranSteps  (3),
    updateTooBig      (-2),
    stalled           (-3),  // (near convergence, but fail anyway)
    wrmsExactZero     (-4),
    innerSolveFailed  (-5)
  {};

public:
  int normTooSmall;        // default = 1
  int normalConvergence;   // default = 2
  int nearConvergence;     // default = 3
  int smallUpdate;         // default = 4
  int nanFail;             // default = -6
  int tooManySteps;        // default = -1 
  int tooManyTranSteps;    // default = 3
  int updateTooBig;        // default = -2
  int stalled;             // default = -3;
  int wrmsExactZero;       // default = -4;
  int innerSolveFailed;    // default = -5;
};

std::ostream & operator<<(std::ostream & os, const N_NLS_ReturnCodes & rc);

#endif // Xyce_N_NLS_ReturnCodes_h

