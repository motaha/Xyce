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
// Filename       : $RCSfile: N_UTL_BreakPoint.h,v $
//
// Purpose        : Provide a basic "breakpoint" class that can be used in
//                  STL sets, with appropriate operators so that the set will
//                  automatically eliminate duplicates on insert, with 
//                  "duplicate" meaning "within a certain tolerance"
//
// Special Notes  : There remains some lingering issue with the assignment
//                  operator, so the STL "unique" algorithm doesn't work yet
//                  
//
// Creator        : Tom Russo, SNL, Component Information and Models
//
// Creation Date  : 04/28/2004
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#ifndef Xyce_N_UTL_BreakPoint_h
#define Xyce_N_UTL_BreakPoint_h
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

namespace Xyce {
namespace Util {

enum BreakpointType 
  {
    SIMPLE_BREAKPOINT,
    PAUSE_BREAKPOINT
  };

class BreakPoint
{

 public:
  // constructors
  // default
  BreakPoint () : time(0.0), type(SIMPLE_BREAKPOINT) {}
  // from a double and type enum
  BreakPoint(double t, BreakpointType ty=SIMPLE_BREAKPOINT);
  // copy
  BreakPoint( const BreakPoint& right);
  // assignment operators
  BreakPoint& operator= ( const BreakPoint& b);
  BreakPoint& operator= ( const double& t);
  // logical overloaded operators
  inline  bool operator==(const BreakPoint b) const;
  inline  bool operator<(const BreakPoint b) const;
  inline  bool operator>(const BreakPoint b) const;
  inline  void set(double t, BreakpointType ty);
  void set(double t, int ty);

  inline double value() const;
  inline BreakpointType bptype() const; 
  inline static void setBPTol(double tol);
  inline static double getBPTol();

 private:
  double time;
  BreakpointType type;
  static double bptol_;
};

//-----------------------------------------------------------------------------
// Function      : BreakPoint::value
// Purpose       : return the double part of a breakpoint
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
inline double BreakPoint::value() const
{
  return time;
}

//-----------------------------------------------------------------------------
// Function      : BreakPoint::bptype
// Purpose       : return the type of a breakpoint
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
inline BreakpointType BreakPoint::bptype() const
{
  return type;
}

//-----------------------------------------------------------------------------
// Function      : BreakPoint::setBPTol
// Purpose       : Sets the static (class) variable bptol_, which 
//                 is used in comparisons between breakpoints
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
inline void BreakPoint::setBPTol(double tol)
{
  bptol_=tol;
}

//-----------------------------------------------------------------------------
// Function      : BreakPoint::getBPTol
// Purpose       : returns the static (class) variable bptol_, which 
//                 is used in comparisons between breakpoints
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
inline double BreakPoint::getBPTol()
{
  return bptol_;
}

//-----------------------------------------------------------------------------
// Function      : BreakPoint::operator==
// Function      : BreakPoint::operator<
// Function      : BreakPoint::operator>
// Purpose       : compares two breakpoints, using the "bptol_" as tolerance
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
// comparison (equality)
inline bool BreakPoint::operator== (const BreakPoint b) const
{ return (fabs(time-b.time) <= bptol_);}

// comparison (less than)
inline bool BreakPoint::operator< (const BreakPoint b) const
{ return (time<b.time && b.time-time>bptol_);}

// comparison (greater than)
inline bool BreakPoint::operator> (const BreakPoint b) const
{ return (time>b.time && time-b.time>bptol_);}

//-----------------------------------------------------------------------------
// Function      : BreakPoint::set
// Purpose       : Used to set the value and type of an existing object
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
inline void BreakPoint::set(double t, BreakpointType ty)
{
  time=t;
  type=ty;
}

} // namespace Util
} // namespace Xyce

typedef Xyce::Util::BreakPoint N_UTL_BreakPoint;

#endif
