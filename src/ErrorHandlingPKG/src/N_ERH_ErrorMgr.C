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
// Filename      : $RCSfile: N_ERH_ErrorMgr.C,v $
//
// Purpose       : This file contains the functions for the N_ERH_ErrorMgr
//                 class.
//
// Special Notes : 
//
// Creator       : Eric Keiter,  SNL, Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.84.2.1 $
//
// Revision Date  : $Date: 2014/02/27 00:52:17 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_Misc.h>

#include <iostream>
#include <sstream>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <N_UTL_Misc.h>
#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Messenger.h>
#include <N_PDS_ParallelMachine.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_ReportHandler.h>

namespace Xyce {
namespace Report {

Parallel::Machine comm_;

void trim(std::string &);

//-----------------------------------------------------------------------------
// Function      : abort
// Purpose       : This function handles the termination of execution.
// Special Notes : The MPI stuff has not been tested yet.
// Scope         :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

void abort()
{
  // Can't flush the parallel output buffer, so flush it to the standard error
  // stream.

  lout() << std::endl;
  std::cerr << std::endl << std::endl << "*** Xyce Abort ***" << std::endl;

  Xyce_exit(-1); // Must do this because we need our MPI_Finalize to happen
}

//-----------------------------------------------------------------------------
// Function      : registerComm
// Purpose       : This function registers the N_PDS_Comm object.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/03/00
//-----------------------------------------------------------------------------
void registerComm(Parallel::Machine comm)
{
  comm_ = comm;
}

//-----------------------------------------------------------------------------
// Function      : safeBarrier
// Purpose       : This barrier will exit cleanly in parallel if one PE exits
//                 with an error.  The region covered is that following a
//                 previous call to startSafeBarrier
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSi
// Creation Date : 07/18/05
//-----------------------------------------------------------------------------
void safeBarrier(Parallel::Machine comm)
{
  // Collect all pending message to the log file.
  pout(comm);

  unsigned count = get_message_count(MSG_FATAL) + get_message_count(MSG_ERROR);

  if (Parallel::is_parallel_run(comm))
    Parallel::AllReduce(comm, MPI_SUM, &count, 1);

  if (count > 0) {
    UserFatal0().die() << "Simulation aborted due to error";

    Xyce_exit(-1);
  }
}

//-----------------------------------------------------------------------------
// Function      : trim
// Purpose       : Trim processor number from front of an error message
// Special Notes :
// Scope         : Private
// Creator       : Dave Shirley, PSSi
// Creation Date : 12/15/05
//-----------------------------------------------------------------------------

void trim(std::string & s)
{
  if (s.size() < 3)
    return;
  if (s[0] == 'P')
  {
    int i = s.find_first_of(':');
    if (i == std::string::npos)
      return;
    s.erase(0,i+2);
  }
}

//-----------------------------------------------------------------------------
// Function      : outputMessage
// Purpose       : Outputs message parameters based on the report_mask and the
//                 message content.
// Special Notes :
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 9/26/00
//-----------------------------------------------------------------------------

void
xyce_report_handler(
  const char *  message,
  unsigned      report_mask)
{
  // if ( !comm_->isSerial() && !(report_mask & MSG_SYMMETRIC))
  //   os << "P" << comm_->procID() << " - ";

  bool fatal = report_mask & MSG_TERMINATE;

  std::ostringstream oss;

  if (Parallel::size(comm_) > 0)
    oss << "Processor " << Parallel::rank(comm_) << ":" << std::endl;

  Util::word_wrap(oss, message, 78, " ", "");

  // If symmetric then all processors are getting the same message, only write to p0 and not to ~p0 backlog.  If
  // asymetric then one processor is getting the message, write to per processor stream which writes to per processor
  // log file and to backlog.
  if (report_mask & MSG_SYMMETRIC)
    lout() << oss.str();
  else
    pout() << oss.str();

  // If fatal error also send the message to the standard error file:
  // Also save it for output on proc 0 if running in parallel
  if (fatal)
  {
    std::cerr << oss.str() << std::endl;
  }

}


} // namespace Report
} // namespace Xyce
