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
// Filename      : $RCSfile: N_ERH_ErrorMgr.h,v $
//
// Purpose       : Contains the class definition for the N_ERH_ErrorMgr
//                 class.  This class is the error handler for Xyce.
//
// Special Notes : For the initial implementation, this stuff was taken
//                 directly from the SIERRA class, Fmwk_Env.
//
// Creator       : Eric Keiter, SNL,  Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.25.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:39 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_ERH_MGR_H
#define Xyce_ERH_MGR_H

// ---------- Standard Includes ----------

#include <iosfwd>
#include <fstream>
#include <string>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_oblackholestream.hpp>

// ---------- Forward Declarations ----------

class N_PDS_Comm;

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_ERH_fwd.h>

using Teuchos::RefCountPtr;
using Teuchos::rcp;

//-----------------------------------------------------------------------------
// Class         : N_ERH_ErrorMgr
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/15/00
//-----------------------------------------------------------------------------
class N_ERH_ErrorMgr
{

  // Attributes:

  // N_PDS_Comm object pointer.
  static RefCountPtr<N_PDS_Comm> comm_;
  static bool inSafeBarrier;
  static Teuchos::oblackholestream myBHS_;

public:

  static ofstream * output;
  static bool       suppress;

  // Reporting errors, warning, and general information.  A report level with
  // the "_0" ending will only generate output from process # 0.

  enum ReportLevel {

    // user error types:
    USR_FATAL       ,
    USR_ERROR       ,
    USR_WARNING     ,
    USR_INFO        ,
    USR_OUT         ,
    USR_FATAL_0     ,
    USR_ERROR_0     ,
    USR_WARNING_0   ,
    USR_INFO_0      ,
    USR_OUT_0       ,

    // developer error types:
    DEV_FATAL       ,
    DEV_ERROR       ,
    DEV_WARNING     ,
    DEV_DEBUG       ,
    DEV_FATAL_0     ,
    DEV_ERROR_0     ,
    DEV_WARNING_0   ,
    DEV_DEBUG_0     ,

    // Progress messages for GUI:
    GUI_PROGRESS
  };

  enum Disposition {
    ERR_TERMINATE ,
    ERR_NORMAL    ,
    ERR_SUPPRESS  ,
    ERR_PROGRESS
  };

  // Functions:
public:

  // Various overloaded versions of the "report" method used to output messages

  //
  static void report(ReportLevel errorType, std::ostringstream &oss);

  // Output single string message
  static void report(ReportLevel errorType, const string & message);

  // Output single string message and a following unsigned integer
  static void report(ReportLevel errorType, const string & message,
                     unsigned int iValue);

  // Output two messages with an unsigned integer between them
  static void report(ReportLevel errorType, const string & message1,
                     unsigned int iValue, const string & message2);

  // Output single string message and a following integer
  static void report(ReportLevel errorType, const string & message1,
                     int iValue);

  // Output two messages with an integer between them
  static void report(ReportLevel errorType, const string & message1,
                     int iValue, const string & message2);

  // Output single string message and a following double
  static void report(ReportLevel errorType, const string & message,
                     double value);

  // Output two messages with an double between them
  static void report(ReportLevel errorType, const string & message1,
                     double value, const string & message2);

  // Output a message with the file name and line number. Used for
  // error reporting during netlist parsing.
  static void report(ReportLevel errorType, const string & message,
                     const string & fileName, const int lineNumber);

  // A FATAL report level will always cause an abort.  The abort procedure is
  // sensitive to the sequential / parallel environment (e.g. MPI_Abort).

  // Make a "good effort" at aborting Xyce on all processors
  static void abort();

  // Takes the message string, formats and outputs the message
  static void format_and_print_error_message(string & out,
                                             const string & message);

  // Buffered Output to command line '-o' file,
  // By default, output on all processors.
  static ostream & Output(const bool & task_0_only = false);

  // Method which registers the comm object
  static void registerComm(N_PDS_Comm *);

  // Method which registers the comm object
  static void registerComm(RefCountPtr<N_PDS_Comm>);

  static void quiet()
  {
    suppress = true;
  }
  static void normal()
  {
    suppress = false;
  }

#ifdef Xyce_PARALLEL_MPI
  static void startSafeBarrier();
  static void safeBarrier(int);
#endif

private:
#ifdef Xyce_PARALLEL_MPI
  static void trim (string &);
#endif

  // Constructs message for output based on the report level, string, etc.
  static void constructMessage(ReportLevel errorType, bool & task_0_only,
                               Disposition & disposition, string &msgHeading);

  // Output constructed message including header and message
  static void outputMessage(const bool & task_0_only, Disposition disposition,
                            string &msgHeading,
                            const string &message);

private:

};

namespace Xyce {
namespace Report {

class Message
{
public:
  Message()
    : errorType(N_ERH_ErrorMgr::USR_FATAL),
      oss(new std::ostringstream)
  {}

  Message(N_ERH_ErrorMgr::ReportLevel error_type)
    : errorType(error_type),
      oss(new std::ostringstream)
  {}

  Message(const Message &right)
    : errorType(right.errorType),
      oss(right.oss)
  {}

  ~Message() {
    N_ERH_ErrorMgr::report(errorType, oss->str());
    delete oss;
  }

private:
  Message &operator=(const Message &);

public:
  std::ostringstream &errorStream() {
    return *oss;
  }

private:
  N_ERH_ErrorMgr::ReportLevel            errorType;
  std::ostringstream *   oss;
};

template<class T>
inline Message &operator<<(Message &eos, const T &t) {
  eos.errorStream() << t;
  return eos;
}

} // namespace Report
} // namespace Xyce

#endif
