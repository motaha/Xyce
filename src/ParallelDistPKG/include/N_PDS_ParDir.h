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
// Filename       : $RCSfile: N_PDS_ParDir.h,v $
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
// Revision Number: $Revision: 1.10.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:48 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_ParDir_h
#define Xyce_N_PDS_ParDir_h

// ---------- Standard Includes ----------
#include <string>
#include <list>
#include <map>
#include <vector>

// ---------- Forward Declarations -------
class N_PDS_Comm;
class N_PDS_Migrator;

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_TOP_ParNode.h>

// ----------   Other Includes   ----------


//-----------------------------------------------------------------------------
// Class         : N_PDS_ParDir
// Purpose       : Distributed directory for object lookup
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
class N_PDS_ParDir
{

public:

  // Constructors
  N_PDS_ParDir(N_PDS_Comm * comm = NULL, N_PDS_Migrator * migrator = NULL);

  // Destructor
  ~N_PDS_ParDir();

private:
  // No public copy construction, assignment, or equality operators

  // Copy constructor (private).
  N_PDS_ParDir(const N_PDS_ParDir & right);

  // Assignment operation (private).
  N_PDS_ParDir & operator = (const N_PDS_ParDir & right);

  // Equality operation (private).
  bool operator == (const N_PDS_ParDir & right) const;

  // Non-equality operation (private).
  bool operator != (const N_PDS_ParDir & right) const;

public:

  // Registration method for the PDS_Comm pointer.
  void registerPDSComm(N_PDS_Comm * comm);

  // Add objects from directory.
  void addItems(const vector < N_TOP_ParNode * > & nodeVec);

  // Remove objects from directory.
  void deleteItems(const vector < NodeID > & idVec);

  // Access directory info.

  // Get the items in the directory.
  void getItems(const vector < NodeID > & idVec,
                vector < N_TOP_ParNode * > & nodeVec);

  // Get the GIDs in the directory.
  void getGIDs(const vector < NodeID > & idVec, vector < int > & gidVec);

  // Get the processors in the directory.
  void getProcs(const vector < NodeID > & idVec, vector < int > & procVec);

  bool debugDump(ostream & os) const;

#ifdef Xyce_PARALLEL_MPI

  // Change ghost listing.
  void addGhosts(const vector < pair < NodeID, int > > & idVec);
  void deleteGhosts(const vector < pair < NodeID, int > > & idVec);
  void clearGhosts(const vector < NodeID > & idVec);
  void clearGhosts();
  void getGhosts(const vector < NodeID > & idVec,
                 vector < vector < int > > & ghostVec);

protected:
  int parKey(const string & token);
#endif

protected:

  // Pointer to the PDS_Comm object.
  N_PDS_Comm * pdsComm_;
  bool commOwned_;

  // Pointer to the PDS_Migrator object.
  N_PDS_Migrator * pdsMigrator_;
  bool migratorOwned_;

  map < string, N_TOP_ParNode > nodeMap_;

};

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::N_PDS_ParDir
// Purpose       : copy constructor
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
inline N_PDS_ParDir::N_PDS_ParDir(const N_PDS_ParDir & right)
{
  // no_op
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::operator=
// Purpose       : assignment operator
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
inline N_PDS_ParDir & N_PDS_ParDir::operator = (const N_PDS_ParDir & right)
{
  return (* this);
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::operator==
// Purpose       : equality operator
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
inline bool N_PDS_ParDir::operator == (const N_PDS_ParDir & right) const
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::operator!=
// Purpose       : inequality operator
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
inline bool N_PDS_ParDir::operator != (const N_PDS_ParDir & right) const
{
  return true;
}

#endif
