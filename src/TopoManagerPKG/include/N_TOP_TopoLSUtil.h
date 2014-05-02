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
// Revision Number: $Revision: 1.51.2.1 $
//
// Revision Date  : $Date: 2014/03/06 00:41:20 $
//
// Current Owner  : $Author: erkeite $
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

#include <N_IO_fwd.h>
#include <N_TOP_fwd.h>
#include <N_PDS_fwd.h>

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_LAS_QueryUtil.h>
#include <N_IO_PkgOptionsMgr.h>


namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : TopoLSUtil
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
class TopoLSUtil : public N_LAS_QueryUtil
{
public:

  // Constructor
  //TopoLSUtil(Topology * topo = 0, IO::CmdParse & cp);
  TopoLSUtil(Topology * topo, IO::CmdParse & cp);

  // Destructor
  ~TopoLSUtil() {}

  // Register the pointer to the topology object.
  bool registerTopology(Topology * topo) { return (topoPtr_ = topo); }

  // Register the pointer to the parallel services manager object.
  bool registerParallelMgr(N_PDS_Manager * pdsmgr) { return (pdsMgrPtr_ = pdsmgr); }

  // Method to register the package options manager
  bool registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr );
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

  // get access to the names file flag
  bool namesFileFlag() {return namesFile_;}

  // Accessor methods.

  const std::vector<int> & vnodeGIDVec() const { return vnodeGIDVector_; }
  const std::vector<int> & vsrcGIDVec() const { return vsrcGIDVector_; }
  //const std::vector<int> & noDCPathGIDVec() const { return noDCPathGIDVector_; }
  //const std::vector<int> & connToOneTermGIDVec() const { return connToOneTermGIDVector_; }

  int numGlobalNodes() const { return numGlobalNodes_; }
  int numLocalNodes() const { return numLocalNodes_; }

  int numGlobalRows() const { return numGlobalRows_; }
  int numLocalRows() const { return numLocalRows_; }
  int numExternRows() const { return numExternRows_; }
  int numGlobalExternRows() const { return numGlobalExternRows_; }
  int baseRowGID() const { return baseRowGID_; }
  const std::vector<int> & rowList_GID() const { return rowList_GID_; }
  const std::vector< std::pair<int,int> > & rowList_ExternGID() const { return rowList_ExternGID_; }

  int numGlobalStateVars() const { return numGlobalStateVars_; }
  int numLocalStateVars() const { return numLocalStateVars_; }
  int numExternStateVars() const { return numExternStateVars_; }
  int numGlobalExternStateVars() const { return numGlobalExternStateVars_; }
  int baseStateVarGID() const { return baseStateVarGID_; }
  const std::vector<int> & rowList_StateGID() const { return rowList_StateGID_; }
  const std::vector< std::pair<int,int> > & rowList_ExternStateGID() const { return rowList_ExternStateGID_; }

  int numGlobalStoreVars() const { return numGlobalStoreVars_; }
  int numLocalStoreVars() const { return numLocalStoreVars_; }
  int numExternStoreVars() const { return numExternStoreVars_; }
  int numGlobalExternStoreVars() const { return numGlobalExternStoreVars_; }
  int baseStoreVarGID() const { return baseStoreVarGID_; }
  const std::vector<int> & rowList_StoreGID() const { return rowList_StoreGID_; }
  const std::vector< std::pair<int,int> > & rowList_ExternStoreGID() const { return rowList_ExternStoreGID_; }

  int numGlobalNZs() const { return numGlobalNZs_; }
  int numLocalNZs() const { return numLocalNZs_; }
  const std::vector<int> & rowList_NumNZs() const { return rowList_NumNZs_; }
  const std::vector< std::list<int> > & rowList_ColList() const { return rowList_ColList_; }

  const std::vector<char> & rowList_VarType() const { return rowList_VarType_; }

  struct TopoLSUtil_OptionsReg : public IO::PkgOptionsReg
  {
    TopoLSUtil_OptionsReg( TopoLSUtil & tlsu )
    : topLSutil(tlsu)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return topLSutil.registerOptions( options ); }

    TopoLSUtil & topLSutil;
  };

  struct TopoLSUtil_TimeOptionsReg : public IO::PkgOptionsReg
  {
    TopoLSUtil_TimeOptionsReg( TopoLSUtil & tlsu )
    : topLSutil(tlsu)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return topLSutil.registerTimeOptions( options ); }

    TopoLSUtil & topLSutil;
  };

private:

  // Pointer to the topology object.
  Topology * topoPtr_;

  // command line object:
  IO::CmdParse & commandLine_;
  
  // package options manager
  IO::PkgOptionsMgr * pkgOptMgrPtr_;
  
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
  std::vector<int> nodeList_GID_;
  std::vector< std::pair<int,int> > nodeList_ExternGID_;
  std::map<int,int> nodeGtoL_Map_;

  // Number of global (across all processors) rows in the linear system.
  int numGlobalRows_;
  // Number of local (on processor) rows in the linear system.
  int numLocalRows_;
  int numExternRows_;
  int numGlobalExternRows_;
  int baseRowGID_;
  std::vector<int> rowList_GID_;
  std::vector< std::pair<int,int> > rowList_ExternGID_;
  std::map<int,int> GtoL_Map_;

  // Number of global (across all processors) state-variables associated with
  // the linear system.
  int numGlobalStateVars_;
  // Number of local (on processor) state-variables associated with the linear
  // system.
  int numLocalStateVars_;
  int numExternStateVars_;
  int numGlobalExternStateVars_;
  int baseStateVarGID_;
  std::vector<int> rowList_StateGID_;
  std::vector< std::pair<int,int> > rowList_ExternStateGID_;

  // Number of global (across all processors) store-variables associated with
  // the linear system.
  int numGlobalStoreVars_;
  // Number of local (on processor) store-variables associated with the linear
  // system.
  int numLocalStoreVars_;
  int numExternStoreVars_;
  int numGlobalExternStoreVars_;
  int baseStoreVarGID_;
  std::vector<int> rowList_StoreGID_;
  std::vector< std::pair<int,int> > rowList_ExternStoreGID_;

  // Variable Type info used for scaling
  std::vector<char> rowList_VarType_;

  //Graph info for jacobian matrix
  int numGlobalNZs_;
  int numLocalNZs_;
  std::vector<int> rowList_NumNZs_;
  std::vector< std::list<int> > rowList_ColList_;

  bool checkConnectivity_;
  bool supernode_;
  bool namesFile_;

  std::vector<int> vnodeGIDVector_;
  std::vector<int> vsrcGIDVector_;

  //Adding these lists to detect the IDs of nodes for which there is
  //no DC path to ground or for which the node is only connected to one device
  //terminal.

  std::vector<std::string> noDCPathIDVector_;
  std::vector<std::string> connToOneTermIDVector_;

private:

  //testing routine for problems with voltage node connectivity
  bool testVoltageNodeConnectivity_();
  void comm_boundaries (std::map<int, std::vector<int> > & gid_map,
                      std::vector<int> & actual_buf_in, std::vector<int> & actual_buf_out,
                      std::vector<int> & buf_len, std::vector<int> & buf_dest,
                      std::vector<int *> & buf_in, std::vector<int *> & buf_out, int mode);

  void outputTopoWarnings(std::vector<int> &, std::map<int,std::string> &, std::string);

  // Don't allow copy construction or assignment.
  // Copy constructor (private)
  TopoLSUtil(const TopoLSUtil & right);
  // Assignment operator (private).
  TopoLSUtil & operator = (const TopoLSUtil & right);

    friend std::ostream & operator << (std::ostream & os, const TopoLSUtil & tlsu);

};

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::TopoLSUtil N_TOP_TopoLSUtil;

#endif
