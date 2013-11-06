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
// Filename       : $RCSfile: N_TOP_TopoLSUtil.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.44.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_TopoLSUtil_h
#define N_TOP_TopoLSUtil_h 1

// ---------- Standard Includes ----------

#include <list>
#include <iosfwd>
#include <vector>
#include <set>
#include <map>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>

#include <N_LAS_QueryUtil.h>

#include <N_IO_PkgOptionsMgr.h>

// ---------- Forward Declarations ----------

class N_PDS_Manager;
class N_PDS_GlobalAccessor;

class N_TOP_Topology;

class N_IO_CmdParse;

//-----------------------------------------------------------------------------
// Class         : N_TOP_TopoLSUtil
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
class N_TOP_TopoLSUtil : public N_LAS_QueryUtil
{
public:

  // Constructor
  //N_TOP_TopoLSUtil(N_TOP_Topology * topo = 0, N_IO_CmdParse & cp);
  N_TOP_TopoLSUtil(N_TOP_Topology * topo, N_IO_CmdParse & cp);

  // Destructor
  ~N_TOP_TopoLSUtil() {}

  // Register the pointer to the topology object.
  bool registerTopology(N_TOP_Topology * topo) { return (topoPtr_ = topo); }

  // Register the pointer to the parallel services manager object.
  bool registerParallelMgr(N_PDS_Manager * pdsmgr) { return (pdsMgrPtr_ = pdsmgr); }

  // Method to register the package options manager
  bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );
  bool registerTimeOptions(const N_UTL_OptionBlock & OB);
  bool registerOptions(const N_UTL_OptionBlock & OB);

  // Generation of Linear system data:

  // Setup row/col data for linear solver including reorder.
  bool setupRowCol();

  // Generate Ordering and Var GIDs.
  bool setupNodeGIDs();

  // Generate Ordering and Var GIDs
  bool setupSolnAndStateGIDs();

  // Generate row/col data for linear solver.
  void generateRowColData();

  // Reorder GIDs using based on orderedNodeList.
  void reorderGIDs();

  // Register External GIDs and generate Migration Plans
  bool setupGlobalAccessors();
  
  // get access to the supernode flag
  bool supernodeFlag() {return supernode_;}

  // Accessor methods.

  const vector<int> & vnodeGIDVec() const { return vnodeGIDVector_; }
  const vector<int> & vsrcGIDVec() const { return vsrcGIDVector_; }
  //const vector<int> & noDCPathGIDVec() const { return noDCPathGIDVector_; }
  //const vector<int> & connToOneTermGIDVec() const { return connToOneTermGIDVector_; }

  int numGlobalNodes() const { return numGlobalNodes_; }
  int numLocalNodes() const { return numLocalNodes_; }

  int numGlobalRows() const { return numGlobalRows_; }
  int numLocalRows() const { return numLocalRows_; }
  int numExternRows() const { return numExternRows_; }
  int numGlobalExternRows() const { return numGlobalExternRows_; }
  int baseRowGID() const { return baseRowGID_; }
  const vector<int> & rowList_GID() const { return rowList_GID_; }
  const vector< pair<int,int> > & rowList_ExternGID() const { return rowList_ExternGID_; }

  int numGlobalStateVars() const { return numGlobalStateVars_; }
  int numLocalStateVars() const { return numLocalStateVars_; }
  int numExternStateVars() const { return numExternStateVars_; }
  int numGlobalExternStateVars() const { return numGlobalExternStateVars_; }
  int baseStateVarGID() const { return baseStateVarGID_; }
  const vector<int> & rowList_StateGID() const { return rowList_StateGID_; }
  const vector< pair<int,int> > & rowList_ExternStateGID() const { return rowList_ExternStateGID_; }

  int numGlobalStoreVars() const { return numGlobalStoreVars_; }
  int numLocalStoreVars() const { return numLocalStoreVars_; }
  int numExternStoreVars() const { return numExternStoreVars_; }
  int numGlobalExternStoreVars() const { return numGlobalExternStoreVars_; }
  int baseStoreVarGID() const { return baseStoreVarGID_; }
  const vector<int> & rowList_StoreGID() const { return rowList_StoreGID_; }
  const vector< pair<int,int> > & rowList_ExternStoreGID() const { return rowList_ExternStoreGID_; }

  int numGlobalNZs() const { return numGlobalNZs_; }
  int numLocalNZs() const { return numLocalNZs_; }
  const vector<int> & rowList_NumNZs() const { return rowList_NumNZs_; }
  const vector< list<int> > & rowList_ColList() const { return rowList_ColList_; }

  const vector<char> & rowList_VarType() const { return rowList_VarType_; }

  struct N_TOP_TopoLSUtil_OptionsReg : public N_IO_PkgOptionsReg
  {
    N_TOP_TopoLSUtil_OptionsReg( N_TOP_TopoLSUtil & tlsu )
    : topLSutil(tlsu)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return topLSutil.registerOptions( options ); }

    N_TOP_TopoLSUtil & topLSutil;
  };

  struct N_TOP_TopoLSUtil_TimeOptionsReg : public N_IO_PkgOptionsReg
  {
    N_TOP_TopoLSUtil_TimeOptionsReg( N_TOP_TopoLSUtil & tlsu )
    : topLSutil(tlsu)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return topLSutil.registerTimeOptions( options ); }

    N_TOP_TopoLSUtil & topLSutil;
  };

private:

  // Pointer to the topology object.
  N_TOP_Topology * topoPtr_;

  // command line object:
  N_IO_CmdParse & commandLine_;
  
  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;
  
  // Pointer to the parallel services manager object.
  N_PDS_Manager * pdsMgrPtr_;

  N_PDS_GlobalAccessor * nodeGlobalAccessorPtr_;
  N_PDS_GlobalAccessor * solnGlobalAccessorPtr_;
  N_PDS_GlobalAccessor * stateGlobalAccessorPtr_;
  N_PDS_GlobalAccessor * storeGlobalAccessorPtr_;

  // Number of global (across all processors) nodes in the topology.
  int numGlobalNodes_;
  // Number of local (on processor) nodes in the topology.
  int numLocalNodes_;
  int baseNodeGID_;
  vector<int> nodeList_GID_;
  vector< pair<int,int> > nodeList_ExternGID_;
  map<int,int> nodeGtoL_Map_;

  // Number of global (across all processors) rows in the linear system.
  int numGlobalRows_;
  // Number of local (on processor) rows in the linear system.
  int numLocalRows_;
  int numExternRows_;
  int numGlobalExternRows_;
  int baseRowGID_;
  vector<int> rowList_GID_;
  vector< pair<int,int> > rowList_ExternGID_;
  map<int,int> GtoL_Map_;

  // Number of global (across all processors) state-variables associated with
  // the linear system.
  int numGlobalStateVars_;
  // Number of local (on processor) state-variables associated with the linear
  // system.
  int numLocalStateVars_;
  int numExternStateVars_;
  int numGlobalExternStateVars_;
  int baseStateVarGID_;
  vector<int> rowList_StateGID_;
  vector< pair<int,int> > rowList_ExternStateGID_;

  // Number of global (across all processors) store-variables associated with
  // the linear system.
  int numGlobalStoreVars_;
  // Number of local (on processor) store-variables associated with the linear
  // system.
  int numLocalStoreVars_;
  int numExternStoreVars_;
  int numGlobalExternStoreVars_;
  int baseStoreVarGID_;
  vector<int> rowList_StoreGID_;
  vector< pair<int,int> > rowList_ExternStoreGID_;

  // Variable Type info used for scaling
  vector<char> rowList_VarType_;

  //Graph info for jacobian matrix
  int numGlobalNZs_;
  int numLocalNZs_;
  vector<int> rowList_NumNZs_;
  vector< list<int> > rowList_ColList_;

  // new DAE boolean.
  bool checkConnectivity_;
  bool supernode_;

  vector<int> vnodeGIDVector_;
  vector<int> vsrcGIDVector_;

  //Adding these lists to detect the IDs of nodes for which there is
  //no DC path to ground or for which the node is only connected to one device
  //terminal.

  vector<string> noDCPathIDVector_;
  vector<string> connToOneTermIDVector_;

private:

  //testing routine for problems with voltage node connectivity
  bool testVoltageNodeConnectivity_();
  void comm_boundaries (map<int, vector<int> > & gid_map,
                      vector<int> & actual_buf_in, vector<int> & actual_buf_out,
                      vector<int> & buf_len, vector<int> & buf_dest,
                      vector<int *> & buf_in, vector<int *> & buf_out, int mode);

  void outputTopoWarnings(vector<int> &, map<int,string> &, string);

  // Don't allow copy construction or assignment.
  // Copy constructor (private)
  N_TOP_TopoLSUtil(const N_TOP_TopoLSUtil & right);
  // Assignment operator (private).
  N_TOP_TopoLSUtil & operator = (const N_TOP_TopoLSUtil & right);

  friend ostream & operator << (ostream & os, const N_TOP_TopoLSUtil & tlsu);

};

#endif
