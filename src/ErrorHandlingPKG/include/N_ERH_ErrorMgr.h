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
// Filename      : $RCSfile: N_ERH_ErrorMgr.h,v $
//
// Purpose       : Contains the class definition for the N_ERH_ErrorMgr
//                 class.  This class is the error handler for Xyce.
//
// Special Notes : 
//
// Creator       : Eric Keiter, SNL,  Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.47.2.2 $
//
// Revision Date  : $Date: 2014/02/27 00:52:16 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_ERH_ErrorMgr_h
#define Xyce_ERH_ErrorMgr_h

// ---------- Standard Includes ----------

#include <iosfwd>
#include <fstream>
#include <string>

#include <N_UTL_Xyce.h>
#include <N_UTL_fwd.h>
#include <N_ERH_fwd.h>
#include <N_PDS_fwd.h>

#include <N_ERH_Message.h>

namespace Xyce {
namespace Report {

void xyce_report_handler(const char *message, unsigned message_typed);
void abort();

// Output message
inline void report(MessageType message_type, const std::string & message) 
{
  Message(message_type) << message;
}

// Method which registers the comm object
void registerComm(Parallel::Machine comm);

void safeBarrier(Parallel::Machine comm);

} // namespace Report
} // namespace Xyce

#endif // Xyce_ERH_ErrorMgr_h
