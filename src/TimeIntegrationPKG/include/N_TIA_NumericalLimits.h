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
// Filename      : $RCSfile: N_TIA_NumericalLimits.h,v $
//
// Purpose       : This file represents an interim fix for ultimate use of the
//                 C++ numeric_limits Template package. Parameters defined are
//                 intended for use in the Time Integration Algorithm.
//
// Special Notes :
//                 This class should be deprecated, as it is a relic of the days
//                 when "numeric_limits" was not universally available in
//                 C++ compilers.  Not only that, but it is also
//                 a relic of the days when *templating* was not well 
//                 supported in all C++ compilers, and it was necessary to
//                 provide a kludged, non-templated version.  In 2012, 
//                 numeric_limits is so universal that it should be possible
//                 to do away with this class altogether. 
//
// Creator       : Buddy Watts
//
// Creation Date : 1/19/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_NUMERICAL_LIMITS_H_
#define Xyce_N_TIA_NUMERICAL_LIMITS_H_


#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef HAVE_CLIMITS
#include <climits>
#else
#include <limits.h>
#endif

#ifdef HAVE_NUMERIC_LIMITS
#include <limits>
#endif

#ifndef DBL_MAX
#ifdef HAVE_VALUES_H
#ifdef VALUES_H_HAS_DBL_MAX
#include <values.h>
#endif
#endif

#ifdef HAVE_LIMITS_H
#ifdef LIMITS_H_HAS_DBL_MAX
#include <limits.h>
#endif
#endif
#ifdef HAVE_FLOAT_H
#ifdef FLOAT_H_HAS_DBL_MAX
#include <float.h>
#endif
#endif
#endif

#if 0
//-----------------------------------------------------------------------------
// Class         : template<class T> class numeric_limits
// Purpose       :
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 1/19/00
//-----------------------------------------------------------------------------


template < class T >
class numeric_limits
{
public:
  inline static T min();
  inline static T max();
  inline static T epsilon();
};

#endif

//-----------------------------------------------------------------------------
// Class         : numeric_limits<float>
// Purpose       :
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 1/19/00
//-----------------------------------------------------------------------------

//class numeric_limits<float> {
class numeric_limitsFLOAT
{
public:
  inline static float min() { return FLT_MIN; }

  inline static float max() { return FLT_MAX; }

  inline static float epsilon() { return 0.2 * pow(10.0, -FLT_DIG); }
};


//-----------------------------------------------------------------------------
// Class         : numeric_limits<double>
// Purpose       :
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 1/19/00
//-----------------------------------------------------------------------------
//class numeric_limits<double> {
class numeric_limitsDOUBLE
{
public:
  inline static double min() { return DBL_MIN; }

  inline static double max() { return DBL_MAX; }

  inline static double epsilon() { return 0.5 * pow(10.0, -DBL_DIG); }
};

#endif

