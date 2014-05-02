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
// Filename       : $RCSfile: N_PDS_ParallelMachine.h,v $
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
// Revision Number: $Revision: 1.6 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_PDS_ParallelMachine_h
#define Xyce_N_PDS_ParallelMachine_h

#include <N_PDS_fwd.h>

namespace Xyce {
namespace Parallel {

#ifdef Xyce_PARALLEL_MPI

int rank(MPI_Comm comm);
int size(MPI_Comm comm);

#else

inline int rank(Machine comm) {
  return 0;
}

inline int size(Machine comm) {
  return 1;
}

#endif

inline bool is_parallel_run(Machine comm) {
  return size(comm) > 1;
}

} // namespace Parallel
} // namespace Xyce

#endif /// Xyce_N_PDS_ParallelMachine_h
