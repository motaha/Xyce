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
// Filename      : $RCSfile: N_TIA_TwoLevelError.C,v $
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
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <iostream>

#include <N_TIA_TwoLevelError.h>
#include <N_ERH_ErrorMgr.h>

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for two level error class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const N_TIA_TwoLevelError & tlerror)
{
  os.width(20);os.precision(12);os.setf(std::ios::scientific);
  os << "\n-----------------------------------------" << std::endl;
  os << "\tTwoLevelError:\n";
  os << "\t    innerSize:\t" << tlerror.innerSize << std::endl;
  os << "\t    xErrorSum:\t" << tlerror.xErrorSum << std::endl;
  os << "\t    qErrorSum:\t" << tlerror.qErrorSum << std::endl;
  os << "\t xErrorSum_m1:\t" << tlerror.xErrorSum_m1 << std::endl;
  os << "\t xErrorSum_m2:\t" << tlerror.xErrorSum_m2 << std::endl;
  os << "\t xErrorSum_p1:\t" << tlerror.xErrorSum_p1 << std::endl;
  os << "\t q1HistorySum:\t" << tlerror.q1HistorySum << std::endl;
  os << Xyce::section_divider << std::endl;
  os << std::endl;

  return os;
}
