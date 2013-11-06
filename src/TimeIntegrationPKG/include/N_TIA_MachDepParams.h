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
// Filename      : $RCSfile: N_TIA_MachDepParams.h,v $
//
// Purpose       : This file defines the machine dependent parameters to be
//                 used in the Time Integration Algorithm.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_MACH_DEP_PARAMS_H_
#define Xyce_N_TIA_MACH_DEP_PARAMS_H_

#ifndef HAVE_NUMERIC_LIMITS
#include <N_TIA_NumericalLimits.h>
#else
#include <limits>
#endif

#include <N_UTL_Misc.h>

#ifndef HAVE_NUMERIC_LIMITS
//-----------------------------------------------------------------------------
// Class         : N_TIA_MachineDependentParams
// Purpose       :
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class N_TIA_MachineDependentParams
{
public:

  // 10 * Minimum Floating Point Value
  inline static double MachineZero()
  { return 10.0 * numeric_limitsDOUBLE::min(); }

  // SquareRoot(Maximum Floating Point Value)
  inline static double MachineBig()
  { return sqrt(numeric_limitsDOUBLE::max()); }

  // 4 * Machine Unit Roundoff
  inline static double MachinePrecision()
  { return 4.0 * numeric_limitsDOUBLE::epsilon(); }

  // Machine roundoff (epsilon)
  inline static double MachineEpsilon()
  { return numeric_limitsDOUBLE::epsilon(); }
};
#else
//-----------------------------------------------------------------------------
// Class         : N_TIA_MachineDependentParams
// Purpose       :
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class N_TIA_MachineDependentParams
{
public:

  // 10 * Minimum Floating Point Value
  inline static double MachineZero()
  { return 10.0 * numeric_limits<double>::min(); }

  // SquareRoot(Maximum Floating Point Value)
  inline static double MachineBig()
  { return sqrt(numeric_limits<double>::max()); }

  // N_TIA_NumericalLimits' numeric_limitsDOUBLE::epsilon does NOT return
  //  "machine epsilon" as in "DBL_EPSILON", it returns .5*10^(-DBL_DIG)
  inline static double MachinePrecision()
  { return 2.0 * pow(10.0,-(numeric_limits<double>::digits10)); }

  // Machine Unit Roundoff
  inline static double MachineEpsilon()
  { return numeric_limits<double>::epsilon(); }
};
#endif
#endif     // Xyce_N_TIA_MACH_DEP_PARAMS_H
