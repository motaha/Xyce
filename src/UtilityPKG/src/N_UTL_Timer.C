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
// Filename      : $RCSfile: N_UTL_Timer.C,v $
//
// Purpose       : This file contains the functions for the N_UTL_Timing
//                 class.
//
// Special Notes : Adapted from the timing routines in Petra.
//
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date : 9/18/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_Timer.h>

#include <N_PDS_Comm.h>

#include <Epetra_Time.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : Timer
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
Timer::Timer(N_PDS_Comm & comm )
{
  aTimer_ = new Epetra_Time( *(comm.petraComm()) );
}

//-----------------------------------------------------------------------------
// Function      : Timer
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/24/04
//-----------------------------------------------------------------------------
Timer::Timer( const Epetra_Comm & comm )
{
  aTimer_ = new Epetra_Time( comm );
}

//-----------------------------------------------------------------------------
// Function      : Timer
// Purpose       : Copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
Timer::Timer(const Timer &right)
{
  // Copy to a new Petra_Time object using Petra's copy constructor.
  aTimer_ = new Epetra_Time( *(right.aTimer_) );
}

//-----------------------------------------------------------------------------
// Function      : ~Timer
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
Timer::~Timer()
{
  delete aTimer_;
}

//-----------------------------------------------------------------------------
// Function      : Timer::wallTime
// Purpose       : Returns the wall-clock time in seconds. A code section can
//                 be timed by putting it between two calls to wallTime and
//                 taking the difference of the times.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
double Timer::wallTime() const
{
  return aTimer_->WallTime();
}

//-----------------------------------------------------------------------------
// Function      : Timer::resetStartTime
// Purpose       : Resets the start time for the timer object to the current
//                 time.  A code section can be timed by putting it between a
//                 call to resetStartTime and elapsedTime.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
void Timer::resetStartTime()
{
  aTimer_->ResetStartTime();
}

//-----------------------------------------------------------------------------
// Function      : Timer::elapsedTime
// Purpose       : Returns the elapsed time in seconds since the timer object
//                 was constructed, or since the resetStartTime function was
//                 called. A code section can be timed by putting it between
//                 the Timer constructor and a call to elapsedTime, or
//                 between a call to resetStartTime and elapsedTime.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
double Timer::elapsedTime() const
{
  return aTimer_->ElapsedTime();
}

} // namespace Util
} // namespace Xyce

