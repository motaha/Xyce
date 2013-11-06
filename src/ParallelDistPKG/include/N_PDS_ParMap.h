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
// Filename       : $RCSfile: N_PDS_ParMap.h,v $
//
// Purpose        : Specification file for abstract base class for the parallel
//                  map data and functions.
//
// Special Notes  : Part of a GoF Abstract Factory.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:48 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_ParMap_h
#define Xyce_N_PDS_ParMap_h

// ---------- Standard Includes ----------

#include <vector>
#include <map>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>

// ----------  Other Includes   ----------

// ----------  Fwd Declarations  ---------

class N_PDS_Comm;

class Epetra_Map;
class Epetra_BlockMap;

//-----------------------------------------------------------------------------
// Class         : N_PDS_ParMap
// Purpose       : Abstract base class for the parallel map data and functions.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/08/00
//-----------------------------------------------------------------------------
class N_PDS_ParMap
{

public:
  // Constructors
  N_PDS_ParMap( int & numGlobalEntities,
                const int & numLocalEntities,
                const vector<int> & lbMap,
                const int index_base = 0,
                N_PDS_Comm * aComm = 0);

  // Destructor
  ~N_PDS_ParMap();

  // Constructor which takes a EPetra map (private).
  N_PDS_ParMap(Epetra_Map * pMap,
               N_PDS_Comm * aComm);

private:

  // Copy constructor (private).
  N_PDS_ParMap(const N_PDS_ParMap & right);

  // Assignment operator (private).
  N_PDS_ParMap & operator=(const N_PDS_ParMap & right);

  // Equality operator (private).
  bool operator==(const N_PDS_ParMap & right) const;

  // Non-equality operator (private).
  bool operator!=(const N_PDS_ParMap & right) const;

public:

  // Number of global "entities" represented as vertices in the graph. These
  // may be, for example, equations for the linear algebra quantities or
  // devices/nodes for the circuit graph.
  int numGlobalEntities() const;

  // Number of local (on processor) "entities".
  int numLocalEntities() const;

  // Indexing base (0 or 1) used for the maps.
  int indexBase() const;

  // Maximum globally-numbered identifier.
  // This is used predominately for block linear algebra.
  int maxGlobalEntity() const;
/*
  // List of globally-numbered "entities" owned on this processor.
  const int * parMap() const;
*/

  // Comm object
  N_PDS_Comm * pdsComm() { return comm_; }

  // Accessor functions (overridden in derived classes) for the pointer to the
  // library map object.
  Epetra_Map * petraMap() { return petraMap_; }

  Epetra_BlockMap * petraBlockMap();

  // dereference global index to get local index
  int globalToLocalIndex(const int & global_index);
  // dereference local index to get global index
  int localToGlobalIndex(const int & local_index);

protected:

  // Pointer to comm object.
  N_PDS_Comm * comm_;
  bool commOwned_;

  // Pointer to Petra map object.
  Epetra_Map * petraMap_;
  bool mapOwned_;

};

#endif
