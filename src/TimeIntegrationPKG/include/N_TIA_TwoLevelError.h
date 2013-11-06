//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Revision Number: $Revision: 1.9.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_TWO_LEVEL_ERROR_H
#define Xyce_N_TIA_TWO_LEVEL_ERROR_H

// ---------- Standard Declarations ----------
#ifdef Xyce_DEBUG_TIME
#include <iostream>
#endif

// ---------- Forward Declarations ----------


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

#ifdef Xyce_DEBUG_TIME
//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for two level error class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const N_TIA_TwoLevelError & tlerror)
{
  os.width(20);os.precision(12);os.setf(std::ios::scientific);
  os << "\n-----------------------------------------\n";
  os << "\tTwoLevelError:\n";
  os << "\t    innerSize:\t" << tlerror.innerSize << std::endl;
  os << "\t    xErrorSum:\t" << tlerror.xErrorSum << std::endl;
  os << "\t    qErrorSum:\t" << tlerror.qErrorSum << std::endl;
  os << "\t xErrorSum_m1:\t" << tlerror.xErrorSum_m1 << std::endl;
  os << "\t xErrorSum_m2:\t" << tlerror.xErrorSum_m2 << std::endl;
  os << "\t xErrorSum_p1:\t" << tlerror.xErrorSum_p1 << std::endl;
  os << "\t q1HistorySum:\t" << tlerror.q1HistorySum << std::endl;
  os << "-----------------------------------------\n";
  os << std::endl;

  return os;
}
#endif // debug time

#endif // Xyce_N_TIA_TWO_LEVEL_ERROR_H

