//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_LogStream.C,v $
//
// Purpose        : Describe the purpose of the contents of the file. If the
//                  contains the header file of a class, provide a clear
//                  description of the nature of the class.
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : David Baur
//
// Creation Date  : 3/28/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.22.2.1 $
//
// Revision Date  : $Date: 2014/03/03 18:29:29 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

/**
 * @file   N_UTL_LogStream.C
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 06:38:04 2013
 *
 * @brief  Logging and diagnostic stream management
 *
 * There are three output stream to help manage logging for a multi-processor run.  These
 * streams provide basic logging, diagnostic logging and message forwarding.  The logging
 * stream provide a few features for massively parallel processing using Tee and
 * Indentation stream buffers.
 *
 * lout() - Regular logging output.
 *
 *   - Before initializeLogStream(rank, size), writes to std::cout.
 *
 *   - After initializeLogStream(rank, size), writes to a tee stream buffer with std::cout
 *     added on rank 0 and nothing added to other rank processors.
 *
 *   - If openLogFile(path, false), writes to path on rank 0 and nothing on the
 *     other rank processors.
 *
 *   - If outputLogFile(path, true), writes to path on rank 0 and to path.<rank>.<size>
 *     on other rank processors.
 *
 *
 * dout() - Diagnostic output
 *
 *   - Before initializeLogStream(rank, size), writes to std::cout.
 *
 *   - After initializeLogStream(rank, size), write to indent stream buffer which then
 *     writes to lout().
 *
 *
 * pout() - Log messages to be forwarded to the rank 0 processor.
 *
 *   - Before initializeLogStream(rank, size), writes to std::cout.
 *
 *   - After initializeLogStream(rank, size), writes to stream buffer with lout() added on
 *     all processors.
 *
 *   - When pout(comm) is called, all messages for each processor are forwarded to the
 *     rank 0 processor and written to the log file, clearing the backlog.
 *
 *
 * The end result is logging happens on all processors.  The rank 0 processor is the
 * regular log output.  If per processor logging is enabled, then each processor gets log
 * output.  Diagnostic messages are written to the log file for each processor.  Processor
 * unique messages are written to pout() and logged when communication is possible and
 * makes sense.
 *
 *
 * Threading support for lout():
 *
 * While Xyce is not implemented to for thread safety, it is possible to execute several
 * instances of Xyce each on a single thread.  To support this, the lout() stream has a
 * thread to stream map that selects a unique output stream for each thread.
 *
 *   - Call initializeLogStreamByThread() to activate the per-thread stream.
 *
 *   - Each thread should call addThreadStream(std::ostream *os) to set the lout() stream
 *     destination for that thread.
 *
 *   - Each thread should call removeThreadStream(std::ostream *os) to remove the stream
 *     as a destination.
 *
 * Note that if no stream is associated with the thread, the original lout() stream is
 * used.
 *
 */

// NEVER, leave this out!
#include <Xyce_config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <N_UTL_PThread.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_UTL_TeeStreamBuf.h>

#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

namespace Xyce {

const char *section_divider  = "------------------------------------------------------------";
const char *subsection_divider = "----------------------------------------";

namespace {

int                             s_rank = 0;                                     ///< Processor rank
int                             s_size = 1;                                     ///< Number of processors
Util::tee_streambuf             s_loutStreambuf;                                ///< Tee'd stream buffer for logging output
Util::indent_streambuf          s_doutStreambuf(std::cout.rdbuf());             ///< Indentation stream buffer for diagnostic output
Util::tee_streambuf             s_poutStreambuf;                                ///< Processor zero forwarder
std::ofstream *                 s_logFileStream;                                ///< Log file output stream

std::ostream                    s_lout(std::cout.rdbuf());                      ///< Logging output stream
std::ostream                    s_dout(std::cout.rdbuf());                      ///< Diagnostic output stream
std::ostream                    s_pout(std::cout.rdbuf());                      ///< Processor zero forwarded output stream
std::ostringstream              s_poutBacklog;                                  ///< Stream holding messages to be forwarded

const char *section_divider  = "------------------------------------------------------------";
const char *subsection_divider = "----------------------------------------";

typedef std::map<Util::xyce_pthread_t, std::ostream *> StreamMap;

Util::xyce_pthread_mutex_t *    s_threadLoutMutex;                              ///< 0 means not threaded, otherwise mutex for m_threadLoutMap
StreamMap                       s_threadLoutMap;                                ///< Destination output streams to write to

struct Lock
{
  Lock(Util::xyce_pthread_mutex_t *lock)
    : lock_(lock)
  {
    if (lock_)
      Util::xyce_pthread_mutex_lock(lock_);
  }

    ~Lock()
    {
      if (lock_)
        Util::xyce_pthread_mutex_unlock(lock_);
    }

    Util::xyce_pthread_mutex_t *lock_;
};

/**
 * Initializes output streams and stream buffers.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 06:37:28 2013
 */
struct INIT
{

  /**
   * Connects the output stream via tee and indentation stream buffers.
   *
   * lout() to std::cout, dout() to lout() via indentation stream buffer, and pout() to
   * lout().  pout() is not initialized to write to the backlog.
   * initializeLogStream(rank, size) handles that since that is the indication that
   * parallel processing is happening and that message forwarding will be done.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Thu Oct  3 07:05:15 2013
   */
  INIT()
  {
    s_loutStreambuf.add(&std::cout);
    s_doutStreambuf.redirect(&s_loutStreambuf);
    s_poutStreambuf.add(&s_lout);

    s_lout.rdbuf(&s_loutStreambuf);
    s_dout.rdbuf(&s_doutStreambuf);
    s_pout.rdbuf(&s_poutStreambuf);
  }

  /**
   * Disconnect the output streams from the tee and indentation stream buffers.
   *
   * lout(), dout() and pout() once again write to std::cout.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Thu Oct  3 07:07:49 2013
   */
  ~INIT()
  {
    s_dout.flush();
    s_pout.flush();
    s_lout.flush();

    s_loutStreambuf.clear();
    s_poutStreambuf.clear();

    delete s_logFileStream;

    s_lout.rdbuf(std::cout.rdbuf());
    s_dout.rdbuf(std::cout.rdbuf());
    s_pout.rdbuf(std::cout.rdbuf());

    s_logFileStream = 0;

    delete s_threadLoutMutex;

    s_threadLoutMutex = 0;
  }
};

/**
 * Let the C++ abi handle the initialization and cleanup.
 *
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:10:01 2013
 */
inline void init()
{
  static INIT s_init;
}

} // namespace <unnamed>


/**
 * Initialize the log streams based on processor rank and size.
 *
 * On rank 0 processor, lout() write to std::cout, nowhere on other ranked processors.
 * pout() writes to the forward message backlog.
 *
 * @param rank This processor's rank
 * @param size Number of processor's in parallel run
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:10:38 2013
 */
void initializeLogStream(int rank, int size)
{
  init();

  s_rank = rank;
  s_size = size;

  if (rank != 0)
  {
    s_poutStreambuf.add(&s_poutBacklog);
    s_loutStreambuf.remove(&std::cout);
  }
}

/**
 * Logging output stream.
 *
 * If threaded output has been activated, the per thread stream is returned if one has
 * been added to the per-thread stream map.
 *
 * @return reference to the logging output stream.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:15:07 2013
 */
std::ostream &lout()
{
  init();

  if (s_threadLoutMutex) {
    Util::xyce_pthread_t thread_id = Util::xyce_pthread_self();

    Lock lock(s_threadLoutMutex);

    StreamMap::iterator it = s_threadLoutMap.find(thread_id);

    if (it != s_threadLoutMap.end())
      return *(*it).second;
  }

  return s_lout;
}

/**
 * Diagnostic output stream.
 *
 * @return reference to the diagnostic output stream.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:15:07 2013
 */
std::ostream &dout()
{
  init();

  return s_dout;
}

/**
 * Processor zero forwarding output stream.
 *
 * @return reference to the processor zero forwarding output stream.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:15:07 2013
 */
std::ostream &pout()
{
  init();

  return s_pout;
}

/**
 * Open output file for logging.
 *
 * Open the specified path on processor rank 0 for the lout() stream.  If per_processor is
 * specified, then the other ranked processors lout() writes to the file path.r.s where r
 * is rank and s is size from initializeLogStream(rank, size).
 *
 * @param path Path of the logging file to open
 * @param per_processor true if per processor logging is to occur
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:16:28 2013
 */
bool openLogFile(const std::string &path, bool per_processor)
{
  // Allocate the output stream
  if (path == "cout")
  {
    s_loutStreambuf.remove(&std::cout);
    if (per_processor || s_rank == 0)
      s_loutStreambuf.add(&std::cout);

    return true;
  }

  std::ofstream *ofs = 0;

  if (per_processor && s_rank != 0)
  {
    std::ostringstream perprocess_path;
    perprocess_path << path << "." << s_rank << "." << s_size;

    ofs = new std::ofstream();
    ofs->open(perprocess_path.str().c_str());
  }
  else if (s_rank == 0)
  {
    ofs = new std::ofstream();
    ofs->open(path.c_str());
  }

  // If opening the stream was unsuccessful, silently ignore logging to the file.
  if (ofs && ofs->fail())
  {
    delete ofs;
    ofs = 0;
  }

  if (ofs)
  {
    s_logFileStream = ofs;
    s_loutStreambuf.add(s_logFileStream);
    s_loutStreambuf.remove(&std::cout);
  }

  return ofs != 0;
}

/**
 * Close output log file.
 *
 * Closes the output log file on all processors.  On processor rank 0, the output is
 * written to cout.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:23:25 2013
 */
void closeLogFile()
{
  pout().flush();
  dout().flush();
  lout().flush();

  s_loutStreambuf.remove(s_logFileStream);

  delete s_logFileStream;

  s_logFileStream = 0;

  if (s_rank != 0)
  {
    s_loutStreambuf.remove(&std::cout);
  }
}

/**
 * Initialize the log streams for thread based logging.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:10:38 2013
 */
void initializeLogStreamByThread()
{
  init();

  s_threadLoutMutex = new Util::xyce_pthread_mutex_t();
}

/**
 * Associate an output stream with the calling thread.
 *
 *
 * @param os    Output stream pointer to associate with the calling thread
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Wed Dec  4 10:13:42 2013
 */
void addThreadStream(std::ostream *os)
{
  if (!s_threadLoutMutex)
    throw std::runtime_error("Must initializeLogStreamByThread() first");

  Util::xyce_pthread_t thread_id = Util::xyce_pthread_self();

  Lock lock(s_threadLoutMutex);

  s_threadLoutMap.insert(StreamMap::value_type(thread_id, os));
}

/**
 * Dissociate the output stream from the calling thread.
 *
 *
 * @param os    Output stream pointer to dissociate with the calling thread
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Wed Dec  4 10:14:36 2013
 */
void removeThreadStream(std::ostream *os)
{
  if (!s_threadLoutMutex)
    throw std::runtime_error("Must initializeLogStreamByThread() first");

  Util::xyce_pthread_t thread_id = Util::xyce_pthread_self();

  Lock lock(s_threadLoutMutex);

  s_threadLoutMap.erase(thread_id);
}

/** 
 * Gather the stored per-processor output to processor rank zero.
 *
 * 
 * @param comm  Communicator to use to gather the messages
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Wed Dec  4 10:15:26 2013
 */
void pout(Parallel::Machine comm)
{
  Parallel::AllWriteString(comm, lout(), s_poutBacklog.str());

  // Clear the backlog
  s_poutBacklog.str("");
  s_poutBacklog.clear();
}

} // namespace Xyce
