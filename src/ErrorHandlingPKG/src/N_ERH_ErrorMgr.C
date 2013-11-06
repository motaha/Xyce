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
// Filename      : $RCSfile: N_ERH_ErrorMgr.C,v $
//
// Purpose       : This file contains the functions for the N_ERH_ErrorMgr
//                 class.
//
// Special Notes : Adapted fromt the error handling  routines in SIERRA.
//
// Creator       : Eric Keiter,  SNL, Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.52.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:39 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <iostream>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <sstream>

#include <stdio.h>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_PDS_ParComm.h>
#include <N_PDS_MPIComm.h>

#include <N_ERH_ErrorMgr.h>

ofstream * N_ERH_ErrorMgr::output    = NULL;
bool      N_ERH_ErrorMgr::suppress = false;

RefCountPtr<N_PDS_Comm> N_ERH_ErrorMgr::comm_;
bool N_ERH_ErrorMgr::inSafeBarrier = false;
Teuchos::oblackholestream N_ERH_ErrorMgr::myBHS_;

#ifdef Xyce_PARALLEL_MPI
string heldMessage;
string lastMessage;
#endif

int charsOut = 0;
bool format = true;

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::report
// Purpose       : This function reports error messages to standard error.  It
//                 also will call the abort function, depending on whether or
//                 not the error is a fatal error.
// Special Notes : ReportLevel is an enumerated data type defined in the
//                 N_ERH_ErrorMgr class definition.  The postfix, "_0" is
//                 attached to half of the members of this type, and indicates
//                 that error messages of these types should only print on
//                 processor #0.
// Scope         :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::report(ReportLevel errorType, const string & message)
{
  bool task_0_only = false;
  Disposition disposition;

  string msgHeading;

  constructMessage(errorType, task_0_only, disposition, msgHeading);
  outputMessage(task_0_only, disposition, msgHeading, message);
}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::report
// Purpose       : Overload version of "report" method above which will take
//                 an unsigned integer argument at the end of the string.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 9/26/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::report(ReportLevel errorType, const string & message,
                            unsigned int iValue)

{

  bool task_0_only = false;
  Disposition disposition;

  string msgHeading, msg;

  ostringstream ost;

  ost << iValue;
  string tmp = ost.str();

  msg = message;
  msg += tmp;

  constructMessage(errorType, task_0_only, disposition, msgHeading);
  outputMessage(task_0_only, disposition, msgHeading, msg);

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::report
// Purpose       : Overload version of "report" method above which will take
//                 an unsigned integer argument at the end of the string and
//                 append another string to it.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 12/1/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::report(ReportLevel errorType,
                            const string & message1, unsigned int iValue,
                            const string & message2)

{

  bool task_0_only = false;
  Disposition disposition;

  string msgHeading, msg;

  ostringstream ost;

  ost << iValue;
  string tmp = ost.str();

  msg = message1;
  msg += tmp;
  msg += message2;

  constructMessage(errorType, task_0_only, disposition, msgHeading);
  outputMessage(task_0_only, disposition, msgHeading, msg);

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::report
// Purpose       : Overload version of "report" method above which will take
//                 an integer argument at the end of the string.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 12/1/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::report(ReportLevel errorType,
                            const string & message1, int iValue)

{

  bool task_0_only = false;
  Disposition disposition;

  string msgHeading, msg;

  ostringstream ost;

  ost << iValue;
  string tmp = ost.str();

  msg = message1;
  msg += tmp;

  constructMessage(errorType, task_0_only, disposition, msgHeading);
  outputMessage(task_0_only, disposition, msgHeading, msg);

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::report
// Purpose       : Overload version of "report" method above which will take
//                 an integer argument at the end of the string and appending
//                 another string at the end.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 12/1/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::report(ReportLevel errorType,
                            const string & message1, int iValue,
                            const string & message2)

{

  bool task_0_only = false;
  Disposition disposition;

  string msgHeading, msg;

  ostringstream ost;

  ost << iValue;
  string tmp = ost.str();

  msg = message1;
  msg += tmp;
  msg += message2;

  constructMessage(errorType, task_0_only, disposition, msgHeading);
  outputMessage(task_0_only, disposition, msgHeading, msg);

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::report
// Purpose       : Overload version of "report" method above which will take
//                 a double argument at the end of the string.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 9/29/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::report(ReportLevel errorType, const string & message,
                            double value)

{

  bool task_0_only = false;
  Disposition disposition;

  string msgHeading, msg;

  ostringstream ost;

  ost << value;
  string tmp = ost.str();

  msg = message;
  msg += tmp;

  constructMessage(errorType, task_0_only, disposition, msgHeading);
  outputMessage(task_0_only, disposition, msgHeading, msg);

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::report
// Purpose       : Overload version of "report" method above which will take
//                 a double argument at the end of the string and append
//                 another string to this end.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 10/2/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::report(ReportLevel errorType,
                            const string & message1, double value,
                            const string & message2)

{

  bool task_0_only = false;
  Disposition disposition;

  string msgHeading, msg;

  ostringstream ost;

  ost << value;
  string tmp = ost.str();

  msg = message1;
  msg += tmp;
  msg += message2;

  constructMessage(errorType, task_0_only, disposition, msgHeading);
  outputMessage(task_0_only, disposition, msgHeading, msg);

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::report
// Purpose       : Overload version of "report" method above which will
//                 prepend the message with the file name and line number
//                 on which the error occurred. Used for error reporting
//                 during netlist parsing.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 9/26/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::report(ReportLevel errorType, const string & message,
                            const string & fileName, const int lineNumber)

{

  bool task_0_only = false;
  Disposition disposition;

  string msgHeading, msg;

  if (lineNumber > 0)
  {
    ostringstream ost;

    ost << lineNumber;
    string tmp = ost.str();

    msg = "Error in file " + fileName + " at or near line " + tmp + "\n";
  }

  msg += message;

  constructMessage(errorType, task_0_only, disposition, msgHeading);
  outputMessage(task_0_only, disposition, msgHeading, msg);

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::constructMessage
// Purpose       : Construct message parameters based on the errorType and the
//                 message content.
// Special Notes :
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 9/26/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::constructMessage(ReportLevel errorType,
                                      bool & task_0_only, Disposition & disposition,
                                      string &msgHeading)

{
  const string userStr    = "User ";
  const string devStr     = "Dev ";
  const string fatalStr   = "Fatal: ";
  const string errorStr   = "Error: ";
  const string warningStr = "Warning: ";
  const string noteStr    = "Note: ";
  const string debugStr   = "Debug: ";

  string              labelProc  = "\n";

#ifdef Xyce_PARALLEL_MPI
  heldMessage = "";
#endif

  charsOut = 0;
  format = true;

  switch (errorType) {

    // check for user errors:
  case USR_FATAL_0:
    msgHeading  = userStr + fatalStr;
    task_0_only = true;
    disposition = ERR_TERMINATE;
    break;

  case USR_FATAL:
    msgHeading = userStr + fatalStr;
    disposition = ERR_TERMINATE;
    break;

  case USR_ERROR_0:
    msgHeading  = userStr + errorStr;
    task_0_only = true;
    disposition = ERR_NORMAL;
    break;

  case USR_ERROR:
    msgHeading = userStr + errorStr;
    disposition = ERR_NORMAL;
    break;

  case USR_WARNING_0:
    msgHeading  = userStr + warningStr;
    task_0_only = true;
    disposition = ERR_NORMAL;
    break;

  case USR_WARNING:
    msgHeading = userStr + warningStr;
    disposition = ERR_NORMAL;
    break;

  case USR_INFO_0:
    msgHeading  = labelProc;  // Special case - no header output.
    task_0_only = true;
    disposition = ERR_NORMAL;
    break;

  case USR_INFO:
    msgHeading = labelProc;  // Special case - only proc # in header output.
    disposition = ERR_NORMAL;
    break;

  case USR_OUT_0:
    msgHeading  = "";
    task_0_only = true;
    format = false;
    disposition = ERR_NORMAL;
    break;

  case USR_OUT:
    msgHeading  = "";
    format = false;
    disposition = ERR_NORMAL;
    break;

    // check for developer errors:
  case DEV_FATAL_0:
    msgHeading  = devStr + fatalStr;
    task_0_only = true;
    disposition = ERR_TERMINATE;
    break;

  case DEV_FATAL:
    msgHeading = devStr + fatalStr;
    disposition = ERR_TERMINATE;
    break;

  case DEV_ERROR_0:
    msgHeading  = devStr + errorStr;
    task_0_only = true;
    disposition = ERR_NORMAL;
    break;

  case DEV_ERROR:
    msgHeading = devStr + errorStr;
    disposition = ERR_NORMAL;
    break;

  case DEV_WARNING_0:
    msgHeading  = devStr + warningStr;
    task_0_only = true;
    disposition = ERR_NORMAL;
    break;

  case DEV_WARNING:
    msgHeading = devStr + warningStr;
    disposition = ERR_NORMAL;
    break;

  case DEV_DEBUG_0:
    msgHeading  = devStr + debugStr;
    task_0_only = true;
    disposition = ERR_NORMAL;
    break;

  case DEV_DEBUG:
    msgHeading = devStr + debugStr;
    disposition = ERR_NORMAL;
    break;

  case GUI_PROGRESS:
    msgHeading = "";
    task_0_only = true;
    disposition = ERR_PROGRESS;
    break;
  }

  if (suppress && disposition == ERR_NORMAL)
  {
    disposition = ERR_SUPPRESS;
    return;
  }

#ifdef Xyce_PARALLEL_MPI

  if ( !comm_->isSerial() )
  {
    ostringstream ost;

    ost << comm_->procID();
    string tmp = ost.str();

    labelProc = "P" + tmp + ": ";
  }

  if (!task_0_only)
    msgHeading = labelProc + msgHeading;
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::outputMessage
// Purpose       : Outputs message parameters based on the errorType and the
//                 message content.
// Special Notes :
// Scope         : Private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 9/26/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::outputMessage(const bool & task_0_only,
                                   Disposition disposition, string &msgHeading,
                                   const string &message)

{
  bool fatal;

  if (disposition == ERR_PROGRESS)
  {
    ofstream statusOut;
    statusOut.open("Xyce.sta");
    if ( !statusOut.is_open() )
    {
      string msg = "Unable to open status file";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
    statusOut << message << endl;;
    statusOut.close();
    return;
  }

#ifdef Xyce_PARALLEL_MPI
  heldMessage = "";
#endif

  charsOut = 0;
  if (disposition == ERR_SUPPRESS)
    return;
  if (disposition == ERR_TERMINATE)
    fatal = true;
  else
    fatal = false;

  // The only condition for which output will not happen is
  //     task_0_only = true AND parallel_rank != 0
  //
  // Output will happen for the following  combinations:
  //
  //     task_0_only = true,  parallel_rank  = 0
  //     task_0_only = false, parallel_rank  = 0
  //     task_0_only = false, parallel_rank != 0
  //
  // Recall that parallel_rank specifies which processor we are on.

  format_and_print_error_message(msgHeading, message);
  if (comm_->isSerial())
  {
    (Output() << msgHeading).flush();
    charsOut += msgHeading.size();
#ifdef Xyce_PARALLEL_MPI
    // we only use this in parallel, so I'm not sure why in serial it's saved
    lastMessage = msgHeading;
#endif
  }
  else
  {
#ifdef Xyce_PARALLEL_MPI
    if (!task_0_only || comm_->procID() == 0)
    {
      if( comm_->procID() == 0 )
        (Output() << msgHeading).flush();
      else
        (Output() << msgHeading).flush();
      charsOut += msgHeading.size();
      lastMessage = msgHeading;
    }
    else
    {
      heldMessage = msgHeading;
    }
    if (fatal)
    {
      heldMessage = msgHeading;
    }
#endif
  }
  // If fatal error also send the message to the standard error file:
  // Also save it for output on proc 0 if running in parallel
  if (charsOut > 0 && fatal)
  {
    (cerr << msgHeading).flush();
  }

  if (fatal) N_ERH_ErrorMgr::abort();

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::format_and_print_error_message
// Purpose       : Formats and prints an error message which includes the error
//                 type string ("out") and the message string ("message").
// Special Notes : The original version of this function comes from Sierra, and
//                 did not use STL.  It has been re-written to use STL, in
//                 particular the string class.
// Scope         :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::format_and_print_error_message(string & out,
                                                    const string &message)

{
  const string whitespace = "                                   "
    "     ";
  const int spaces = whitespace.length();
  const int chars_per_line = 78;

  int indent = out.length();

  int line_length      = chars_per_line - indent;
  int lines            = 0;
  int last_white_space = 0;

  const char* char_ptr = message.c_str();
  const char* out_ptr  = message.c_str();

  char b2[2]; b2[1] = 0;

  if (format)
  {
    // Loop over the message string.  Get the array index for the last
    // whitespace.
    while (*char_ptr)
    {
      if (*char_ptr == ' ' || *char_ptr == '\t')
        last_white_space = char_ptr - out_ptr;

      if (char_ptr - out_ptr >= line_length && last_white_space >= 0)
      {

        //if (lines) out += &whitespace[spaces-indent];
        if (lines) out += whitespace[spaces-indent];
        while (*out_ptr == ' ' || *out_ptr == '\t')
        {
          last_white_space--; ++out_ptr;
        }
        while (last_white_space-- && *out_ptr != '\n')
        {
          b2[0] = *out_ptr++;
          out += b2;
        }
        out += "\n";
        ++out_ptr;
        ++lines;
      }
      ++char_ptr;
    }

    if (*out_ptr)
    {
      if (lines)
        //out += &whitespace[spaces-indent];
        out += whitespace[spaces-indent];

      while (*out_ptr)
      {
        b2[0] = *out_ptr;
        out += b2;

         if (*out_ptr == '\n')
           //out += &whitespace[spaces-indent];
           out += whitespace[spaces-indent];
        ++out_ptr;
      }

      out += "\n";
    }
  }
  else
  {
    out = message + "\n";
  }

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::abort
// Purpose       : This function handles the termination of execution.
// Special Notes : The MPI stuff has not been tested yet.
// Scope         :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::abort()

{

  // Can't flush the parallel output buffer, so flush it to the standard error
  // stream.

#if 0
#ifdef Xyce_PARALLEL_MPI
  if (NULL != output)
  {
    cerr << endl
              << "*** Xyce Abort *** Contents of local output buffer to "
      "follow."
              << endl << endl;

    const char * buf     = NULL;
    size_t       buf_len = 0;

    output_buf.get_buffer(buf, buf_len);

    cerr.write(buf, buf_len);
  }
#endif
#endif

  (Output() << endl).flush();
  (cerr << endl << endl << "*** Xyce Abort ***" << endl).flush();

#ifdef Xyce_PARALLEL_MPI
  if (MPI_COMM_NULL !=
	(reinterpret_cast<N_PDS_ParComm*>(comm_.get()))->mpiComm()->comm())
  {
    // Don't do this:
    //    MPI_Finalize();
    //  It'll be done by Xyce_exit(-1) anyway.
    // that had replaced
    //    MPI_Abort();
    // which apparently caused voluminous stack trace output on SGI
    if (N_ERH_ErrorMgr::inSafeBarrier)
    {
      // sync up with processors that have not encountered errors, but indicate to
      // all that an error has occurred
      N_ERH_ErrorMgr::safeBarrier(-1);
    }
  }
  else
  {
    (cerr << endl << "MPI not available, calling exit(-1)" <<
     endl).flush();
  }
#endif

  Xyce_exit(-1); // Must do this because we need our MPI_Finalize to happen

}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::Output
// Purpose       : Output function...
// Special Notes :
// Scope         :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

ostream & N_ERH_ErrorMgr::Output(const bool & task_0_only)
{
  if (task_0_only) {
    if( comm_->procID() == 0 )
      return NULL != output ? (ostream &) *output : (ostream &) cout;
    else
      return myBHS_;
  }
  else
      return NULL != output ? (ostream &) *output : (ostream &) cout;
}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::registerComm
// Purpose       : This function registers the N_PDS_Comm object.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/03/00
//-----------------------------------------------------------------------------
void N_ERH_ErrorMgr::registerComm(N_PDS_Comm *comm)
{
  RefCountPtr<N_PDS_Comm> commRCPtr = rcp(comm,false);
  N_ERH_ErrorMgr::registerComm(commRCPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::registerComm
// Purpose       : This function registers the N_PDS_Comm object.
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/11/08
//-----------------------------------------------------------------------------
void N_ERH_ErrorMgr::registerComm(RefCountPtr<N_PDS_Comm> comm)
{
  N_ERH_ErrorMgr::comm_ = comm;

#ifdef Xyce_PARALLEL_MPI
  heldMessage = "";
#endif

  charsOut = 0;
  N_ERH_ErrorMgr::inSafeBarrier = false;
}

#ifdef Xyce_PARALLEL_MPI
//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::startSafeBarrier
// Purpose       : This is called before code that may error out on one PE so
//                 that the other PEs can exit cleanly from the safeBarrier()
//                 following the dangerous code segment
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSi
// Creation Date : 07/18/05
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::startSafeBarrier()
{
  heldMessage = "";
  charsOut = 0;
  N_ERH_ErrorMgr::inSafeBarrier = true;
}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::safeBarrier
// Purpose       : This barrier will exit cleanly in parallel if one PE exits
//                 with an error.  The region covered is that following a
//                 previous call to startSafeBarrier
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSi
// Creation Date : 07/18/05
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::safeBarrier(int stat)
{
  if( comm_->numProc() > 1 )
  {
		int sum;
		comm_->sumAll(&stat, &sum, 1);
		if (sum != 0)
		{
			// Death is imminent.  Output any held messages that may be present due to
			// errors limited to processor zero or other conditions.  Output is limited
			// to messages that are different than what was (or not) output on proc 0.
			int procCnt = comm_->numProc();
			int procID = comm_->procID();
			int i;
			int *heldSize = new int[procCnt];
			int *heldSize_all = new int[procCnt];

			for (i=0 ; i<procCnt ; ++i)
				heldSize[i] = 0;
			heldSize[procID] = heldMessage.size();
			comm_->sumAll(heldSize, heldSize_all, procCnt);
			if (procID == 0)
			{
				trim(heldMessage);
				trim(lastMessage);
				string msg;
				bool out=false;
				for (i=1 ; i<procCnt ; ++i)
				{
					if (heldSize_all[i] > 0)
					{
						char *buf = new char[heldSize_all[i]];
						comm_->recv (buf, heldSize_all[i], i);
						msg.assign(buf,heldSize_all[i]);
						trim(msg);
						if (msg != heldMessage && msg != lastMessage)
						{
							if (!out)
							{
								out = true;
								Output() << endl << "Messages from other procs follow:";
							}
							Output() << endl;
							Output() << "P" << i << ": ";
							Output() << msg;
						}
					}
				}
				Output() << endl;
				Output() << "Please note that duplicate error messages in parallel may result on" << endl;
				Output() << "some systems, and do not necessarily indicate duplicate errors.";
				Output() << endl;
				Output().flush();
			}
			else
			{
				if (heldSize_all[procID] > 0)
				{
					comm_->send (heldMessage.c_str(), heldSize_all[procID], 0);
				}
			}
			comm_->barrier();
			Xyce_exit(sum);
		}
    heldMessage = "";
		charsOut = 0;
		N_ERH_ErrorMgr::inSafeBarrier = false;
	}
}

//-----------------------------------------------------------------------------
// Function      : N_ERH_ErrorMgr::trim
// Purpose       : Trim processor number from front of an error message
// Special Notes :
// Scope         : Private
// Creator       : Dave Shirley, PSSi
// Creation Date : 12/15/05
//-----------------------------------------------------------------------------

void N_ERH_ErrorMgr::trim (string & s)
{
  if (s.size() < 3)
    return;
  if (s[0] == 'P')
  {
    int i = s.find_first_of(':');
    if (i == string::npos)
      return;
    s.erase(0,i+2);
  }
}
#endif
