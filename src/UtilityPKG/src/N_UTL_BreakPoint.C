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
// Revision Number: $Revision: 1.6.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// These are included because the SGI won't accept "using namespace std;" 
// in N_UTL_Xyce.h unless there's something that actually uses the std
// namespace.  Apparently the math header doesn't.
#include <iostream>
#include <string>

#include <N_UTL_Xyce.h>
#include <N_UTL_BreakPoint.h>

// Set the thing for the first time, so the static (class) variable gets
// instantiated
double N_UTL_BreakPoint::bptol_=1.0e-20;

//-----------------------------------------------------------------------------
// Function      : N_UTL_BreakPoint::N_UTL_BreakPoint
// Purpose       : constructor
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
N_UTL_BreakPoint::N_UTL_BreakPoint(double t, BreakpointType ty):time(t),type(ty)
{}

//-----------------------------------------------------------------------------
// Function      : N_UTL_BreakPoint::N_UTL_BreakPoint
// Purpose       : copy constructor
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
N_UTL_BreakPoint::N_UTL_BreakPoint(const N_UTL_BreakPoint& right)
  :time(right.time),
   type(right.type)
{}

//-----------------------------------------------------------------------------
// Function      : N_UTL_BreakPoint::operator=
// Purpose       : copy assignment
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
N_UTL_BreakPoint&  N_UTL_BreakPoint::operator= ( const N_UTL_BreakPoint& b)
{ time=b.time; type=b.type; return *this;}

//-----------------------------------------------------------------------------
// Function      : N_UTL_BreakPoint::N_UTL_BreakPoint
// Purpose       : conversion assignment from double
// Special Notes : allows one to assign a double to a breakpoint, defaults
//                 the type to simple
//                 Use the "set" method instead to set both value and type
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
N_UTL_BreakPoint&  N_UTL_BreakPoint::operator= ( const double& t)
{ time=t; type=SIMPLE_BREAKPOINT; return *this;}


//-----------------------------------------------------------------------------
// Function      : N_UTL_BreakPoint::set
// Purpose       : Set function that lets one give an int instead of an enum
//                 for the type.  Does that by kludging.
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/28/2004
//-----------------------------------------------------------------------------
void N_UTL_BreakPoint::set(double t, int ty)
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
