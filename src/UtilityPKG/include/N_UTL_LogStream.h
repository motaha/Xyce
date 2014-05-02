//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_UTL_LogStream.h,v $
//
// Purpose        : Output stream the uses the indentation stream buffer
//
// Special Notes  : 
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13.2.1 $
//
// Revision Date  : $Date: 2014/03/03 18:29:29 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_LogStream_h
#define Xyce_N_UTL_LogStream_h

#include <iosfwd>

#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_UTL_IndentStreamBuf.h>

namespace Xyce {

std::ostream &lout();
std::ostream &dout();
std::ostream &pout();

void pout(Parallel::Machine comm);

void initializeLogStream(int rank, int size);

bool openLogFile(const std::string &path, bool per_processor);
void closeLogFile();

void initializeLogStreamByThread();
void addThreadStream(std::ostream *os);
void removeThreadStream(std::ostream *os);

} // namespace Xyce

#endif // Xyce_N_UTL_LogStream_h
