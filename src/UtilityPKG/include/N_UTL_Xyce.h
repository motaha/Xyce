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
// Filename       : $RCSfile: N_UTL_Xyce.h,v $
//
// Purpose        : This file is similar to Petra_Petra.h.  It  contains
//                  a number of Xyce-wide definitions and includes.
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.22 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Xyce_H
#define Xyce_N_UTL_Xyce_H

namespace Xyce {

#ifdef Xyce_DEBUG_DEVICE
static const int DEBUG_DEVICE = 1;
#else
static const int DEBUG_DEVICE = 0;
#endif

#ifdef Xyce_DEBUG_ANALYSIS
static const int DEBUG_ANALYSIS = 1;
#else
static const int DEBUG_ANALYSIS = 0;
#endif

#ifdef Xyce_DEBUG_IO
static const int DEBUG_IO = 1;
#else
static const int DEBUG_IO = 0;
#endif

#ifdef Xyce_DEBUG_EXPRESSION
static const int DEBUG_EXPRESSION = 1;
#else
static const int DEBUG_EXPRESSION = 0;
#endif

#ifdef Xyce_DEBUG_RESTART
static const int DEBUG_RESTART = 1;
#else
static const int DEBUG_RESTART = 0;
#endif

#ifdef Xyce_DEBUG_TIME
static const int DEBUG_TIME = 1;
#else
static const int DEBUG_TIME = 0;
#endif

#ifdef Xyce_DEBUG_CIRCUIT
static const int DEBUG_CIRCUIT = 1;
#else
static const int DEBUG_CIRCUIT = 0;
#endif

#ifdef Xyce_DEBUG_NONLINEAR
static const int DEBUG_NONLINEAR = 1;
#else
static const int DEBUG_NONLINEAR = 0;
#endif

#ifdef Xyce_DEBUG_DISTRIBUTION
static const int DEBUG_DISTRIBUTION = 1;
#else
static const int DEBUG_DISTRIBUTION = 0;
#endif

#ifdef Xyce_VERBOSE_TIME
static const int VERBOSE_TIME = 1;
#else
static const int VERBOSE_TIME = 0;
#endif

} // namespace Xyce

#endif // Xyce_N_UTL_Xyce_H
