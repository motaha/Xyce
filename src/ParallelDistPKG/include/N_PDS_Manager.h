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
// Filename       : $RCSfile: N_PDS_Manager.h,v $
//
// Purpose        : Manager for the parallel load-balance, distribution and
//                  migration tools.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.35 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Manager_h
#define Xyce_N_PDS_Manager_h

#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>

class N_PDS_GlobalAccessor;
class N_PDS_ParMap;
class N_PDS_Migrator;
class N_PDS_ParDir;

// ---------- Standard Includes ----------

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <vector>
#include <string>
#include <map>

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

// ----------   Fwd Declares    ----------

class Epetra_CrsGraph;

namespace EpetraExt {

class CrsGraph_View;

}

//-----------------------------------------------------------------------------
// Class         : N_PDS_Manager
// Purpose       : Manager for the parallel load-balance, distribution and
//                 migration tools.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 03/08/00
//-----------------------------------------------------------------------------
class N_PDS_Manager
{

public:
  //Constructors
    N_PDS_Manager(bool & isSerial, bool & procFlag, int iargs = 0, char **cargs = 0, Xyce::Parallel::Machine comm = 0);

  //Destructor
  ~N_PDS_Manager();

  bool registerTopology(N_TOP_Topology * topo);

  // Parallel Map manipulators
  // add map to stl::map container.
  bool addParallelMap(const std::string & name, N_PDS_ParMap * map);

  // Delete map from stl::map container.
  bool deleteParallelMap(const std::string & name);

  // Gets the parallel map object.
  N_PDS_ParMap * getParallelMap(const std::string & name);

  // Create a parallel map from a global id array:
  // if num_global = -1, num_global will be calculated and returned
  N_PDS_ParMap * createParallelMap(int & num_global,
                                   const int & num_local,
                                   const std::vector<int> & gid_map,
                                   const int index_base = 0 );
  
  // Add a global accessor object associated with a parallel map.
  bool addGlobalAccessor(const std::string & name);

  // Delete a global accessor object associated with a parallel map.
  bool deleteGlobalAccessor(const std::string & name);

  // Method which greats a global accessor object.
  N_PDS_GlobalAccessor * createGlobalAccessor(const std::string & name = "");

  // Gets a global accessor object.
  N_PDS_GlobalAccessor * getGlobalAccessor(const std::string & name);

  // Return ptr to common Comm object owned by Mgr.
  N_PDS_Comm * getPDSComm() const { return Comm_; }
  
  // Return ptr to Topology object
  N_TOP_Topology * getTopology() const { return Topo_; }

  N_PDS_Migrator * createMigrator() const;
  N_PDS_ParDir * createParDir() const;

  // Add a matrix graph object
  bool addMatrixGraph( const std::string & name,
                       Epetra_CrsGraph * graph,
                       EpetraExt::CrsGraph_View * trans = 0 );

  // Gets a matrix graph object.
  Epetra_CrsGraph * getMatrixGraph(const std::string & name);

  // Delete a matrix graph object
  bool deleteMatrixGraph(const std::string & name);

private:

  // Copy constructor - private for "singleton" class.
  N_PDS_Manager(const N_PDS_Manager & right);

  // Assignment operator - private for "singleton" class.
  N_PDS_Manager & operator = (const N_PDS_Manager & right);

  // Equality operation - private for "singleton" class.
  bool operator == (const N_PDS_Manager & right) const;

  // Non-equality operation - private for "singleton" class.
  bool operator != (const N_PDS_Manager & right) const;

private:

  // Pointer to the PDS_Comm object.
  N_PDS_Comm * Comm_;

  // Pointer to the topology object.
  N_TOP_Topology * Topo_;

  std::map<std::string, N_PDS_ParMap*> pm_Map_;
  std::map<std::string, N_PDS_GlobalAccessor*> ga_Map_;

  std::map< std::string, Epetra_CrsGraph*> mg_Map_;
  std::map< std::string, EpetraExt::CrsGraph_View* > mgvt_Map_;

};

#endif
