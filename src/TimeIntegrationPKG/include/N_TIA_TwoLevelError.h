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

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_TIA_TwoLevelError.h,v $
// Purpose       : 
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 1/29/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_TWO_LEVEL_ERROR_H
#define Xyce_N_TIA_TWO_LEVEL_ERROR_H

#include <iosfwd>

//-----------------------------------------------------------------------------
// Class         : N_TIA_TwoLevelError
//
// Purpose       : This class contains error information from an inner solve.
//                 Any wrms norm that is taken in the course of time integration
//                 needs to include information from the inner problems.
//                 This class contains that info for a single inner problem.
//                 There will generally be a class like this for each 
//                 inner problem.
//
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 1/29/07
//-----------------------------------------------------------------------------
class N_TIA_TwoLevelError
{
public:
  N_TIA_TwoLevelError():
    xErrorSum(0.0),
    qErrorSum(0.0),
    innerSize(0.0),
    xErrorSum_m1(0.0),
    xErrorSum_m2(0.0),
    xErrorSum_p1(0.0),
    q1HistorySum(0.0)
  { };

  virtual ~N_TIA_TwoLevelError() {};

  double xErrorSum;
  double qErrorSum;
  double innerSize;

  double xErrorSum_m1;
  double xErrorSum_m2;
  double xErrorSum_p1;

  double q1HistorySum;
};

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for two level error class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const N_TIA_TwoLevelError & tlerror);

#endif // Xyce_N_TIA_TWO_LEVEL_ERROR_H

