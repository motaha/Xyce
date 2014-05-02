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
// Filename       : $RCSfile: N_UTL_BreakPoint.C,v $
//
// Purpose        : Essential methods for N_UTL_BreakPoint class
//
// Special Notes  : 
//
// Creator        : Tom Russo, SNL, Component Information and Models
//
// Creation Date  : 04/28/2004
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <iostream>
#include <string>

#include <N_UTL_Xyce.h>
#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Util {

// Set the thing for the first time, so the static (class) variable gets
// instantiated
double BreakPoint::bptol_=1.0e-20;

//-----------------------------------------------------------------------------
// Function      : BreakPoint::BreakPoint
// Purpose       : constructor
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
BreakPoint::BreakPoint(double t, BreakpointType ty):time(t),type(ty)
{}

//-----------------------------------------------------------------------------
// Function      : BreakPoint::BreakPoint
// Purpose       : copy constructor
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
BreakPoint::BreakPoint(const BreakPoint& right)
  :time(right.time),
   type(right.type)
{}

//-----------------------------------------------------------------------------
// Function      : BreakPoint::operator=
// Purpose       : copy assignment
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
BreakPoint&  BreakPoint::operator= ( const BreakPoint& b)
{ time=b.time; type=b.type; return *this;}

//-----------------------------------------------------------------------------
// Function      : BreakPoint::BreakPoint
// Purpose       : conversion assignment from double
// Special Notes : allows one to assign a double to a breakpoint, defaults
//                 the type to simple
//                 Use the "set" method instead to set both value and type
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
BreakPoint&  BreakPoint::operator= ( const double& t)
{ time=t; type=SIMPLE_BREAKPOINT; return *this;}


//-----------------------------------------------------------------------------
// Function      : BreakPoint::set
// Purpose       : Set function that lets one give an int instead of an enum
//                 for the type.  Does that by kludging.
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
void BreakPoint::set(double t, int ty)
{
  time=t;
  switch (ty)
  {
  case 0:
    type=SIMPLE_BREAKPOINT;
    break;
  case 1:
    type=PAUSE_BREAKPOINT;
    break;
  default:
    type=SIMPLE_BREAKPOINT;
  }
}

} // namespace Util
} // namespace Xyce
