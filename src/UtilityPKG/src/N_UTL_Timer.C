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
// Revision Number: $Revision: 1.8.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Timer.h>

#include <N_PDS_Comm.h>

// ---------  Other Includes  -----------

#include <Epetra_Time.h>

// Class N_UTL_Timer

//-----------------------------------------------------------------------------
// Function      : N_UTL_Timer
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
N_UTL_Timer::N_UTL_Timer(N_PDS_Comm & comm )
{
  aTimer_ = new Epetra_Time( *(comm.petraComm()) );
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Timer
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/24/04
//-----------------------------------------------------------------------------
N_UTL_Timer::N_UTL_Timer( const Epetra_Comm & comm )
{
  aTimer_ = new Epetra_Time( comm );
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Timer
// Purpose       : Copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
N_UTL_Timer::N_UTL_Timer(const N_UTL_Timer &right)
{
  // Copy to a new Petra_Time object using Petra's copy constructor.
  aTimer_ = new Epetra_Time( *(right.aTimer_) );
}

//-----------------------------------------------------------------------------
// Function      : ~N_UTL_Timer
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
N_UTL_Timer::~N_UTL_Timer()
{
  delete aTimer_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Timer::wallTime
// Purpose       : Returns the wall-clock time in seconds. A code section can
//                 be timed by putting it between two calls to wallTime and
//                 taking the difference of the times.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
double N_UTL_Timer::wallTime() const
{
  return aTimer_->WallTime();
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Timer::resetStartTime
// Purpose       : Resets the start time for the timer object to the current
//                 time.  A code section can be timed by putting it between a
//                 call to resetStartTime and elapsedTime.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
void N_UTL_Timer::resetStartTime()
{
  aTimer_->ResetStartTime();
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Timer::elapsedTime
// Purpose       : Returns the elapsed time in seconds since the timer object
//                 was constructed, or since the resetStartTime function was
//                 called. A code section can be timed by putting it between
//                 the N_UTL_Timer constructor and a call to elapsedTime, or
//                 between a call to resetStartTime and elapsedTime.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/27/00
//-----------------------------------------------------------------------------
double N_UTL_Timer::elapsedTime() const
{
  return aTimer_->ElapsedTime();
}

