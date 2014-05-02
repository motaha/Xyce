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
// Filename       : $RCSfile: N_PDS_CommFactory.h,v $
//
// Purpose        : Specification file for the comm abstract factory.
//
// Special Notes  : GoF Abstract Factory
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/26/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_CommFactory_h
#define Xyce_N_PDS_CommFactory_h

#include <list>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <N_PDS_fwd.h>

class Epetra_Comm;

//-----------------------------------------------------------------------------
// Class         : N_PDS_CommFactory
// Purpose       : Comm abstract factory (GoF).
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
class N_PDS_CommFactory
{
private:

  // Default constructor.
  N_PDS_CommFactory() { }

  // Copy constructor.
  N_PDS_CommFactory(const N_PDS_CommFactory & c) { }

  // Assignment operator.
  N_PDS_CommFactory & operator = (const N_PDS_CommFactory & c)
    { return * this; }

public:

  // Return a new ParMap
    static N_PDS_Comm * create(int iargs = 0, char * cargs[] = 0, Xyce::Parallel::Machine comm = 0);
};

#endif

