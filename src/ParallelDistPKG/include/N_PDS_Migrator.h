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
// Filename       : $RCSfile: N_PDS_Migrator.h,v $
//
// Purpose        : Migrator tool using Zoltan utilities
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/07/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:48 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Migrator_h
#define Xyce_N_PDS_Migrator_h

// ---------- Standard Includes ----------

#include <string>
#include <vector>

// ---------- Forward Declarations -------

class N_PDS_Comm;
class Packable;

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

// ----------   Other Includes   ----------


//-----------------------------------------------------------------------------
// Class         : N_PDS_Migrator
// Purpose       : Migrator tool
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/23/01
//-----------------------------------------------------------------------------
class N_PDS_Migrator
{

public:

  // Constructor
  N_PDS_Migrator(N_PDS_Comm * comm);

  // Destructor
  ~N_PDS_Migrator();

private:
  // No public copy construction, assignment, or equality operators.

  // Copy constructor (private).
  N_PDS_Migrator(const N_PDS_Migrator & right);

  // Assignment operator (private).
  N_PDS_Migrator & operator = (const N_PDS_Migrator & right);

  // Equality operation (private).
  bool operator == (const N_PDS_Migrator & right) const;

  // Non-equality operation (private).
  bool operator != (const N_PDS_Migrator & right) const;

public:

#ifdef Xyce_PARALLEL_MPI

  void migratePackable(const vector < int > & procVec,
                       const vector < Packable * > & exportVec,
                       vector < Packable * > & importVec);
  void reverseMigratePackable(const vector < int > & procVec,
                              const vector < Packable * > & exportVec,
                              vector < Packable * > & importVec);

  void migrateString(const vector < int > & procVec,
                     const vector < string > & exportVec,
                     vector < string > & importVec);

  void migrateInt(const vector < int > & procVec,
                  const vector < pair < string, int > > & exportVec,
                  vector < pair < string, int > > & importVec);
  void reverseMigrateInt(const vector < int > & procVec,
                         const vector < pair < string, int > > & exportVec,
                         vector < pair < string, int > > & importVec);

  void migrateIntVec(const vector < int > & procVec,
                  const vector < pair < string, vector < int > > > & exportVec,
                     vector < pair < string, vector < int > > > & importVec);
  void reverseMigrateIntVec(const vector < int > & procVec,
  	const vector < pair < string, vector < int > > > & exportVec,
  	vector < pair < string, vector < int > > > & importVec);

#endif

protected:

  // Pointer to the PDS_Comm object.
  N_PDS_Comm * pdsComm_;
  bool commOwned_;

};

#endif
