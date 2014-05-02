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
// Filename       : $RCSfile: N_UTL_Timer.h,v $
//
// Purpose        : Contains wrappers for timing classes.
//
// Special Notes  : Initial implementation is based upon the Petra timing
//                  classes.
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 9/18/00
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

#ifndef Xyce_UTL_Timing_H
#define Xyce_UTL_Timing_H

#include <N_PDS_fwd.h>

class Epetra_Comm;
class Epetra_Time;

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Class         : Timer
// Purpose       : Wraps Petra timing class.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
//-----------------------------------------------------------------------------
class Timer
{

public:

  // Constructors
  Timer( N_PDS_Comm & comm );

  Timer( const Epetra_Comm & comm );

  // Destructor
  ~Timer();

private:

  // Copy constructor (private).
  Timer(const Timer & right);

  // Assignment operator (private).
  Timer & operator = (const Timer & right);

  // Equality Operations (private).
  int operator == (const Timer & right) const;
  int operator != (const Timer & right) const;

public:

  // Wall-clock time function
  double wallTime() const;

  // Resets the start time for a timer object
  void resetStartTime();

  // Elapsed time function
  double elapsedTime() const;

private:

  // Pointer to timer class (Petra default)
  Epetra_Time * aTimer_;

};

} // namespace Util
} // namespace Xyce

typedef Xyce::Util::Timer N_UTL_Timer;

#endif // Xyce_UTL_Timing_H
