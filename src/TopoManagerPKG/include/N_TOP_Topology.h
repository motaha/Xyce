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
// Filename       : $RCSfile: N_TOP_Topology.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/20/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.64.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_Topology_h
#define N_TOP_Topology_h 1

// ---------- Standard Includes ----------
#include <string>
#include <list>
#include <map>
#include <set>
#include <iosfwd>
#include <vector>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_UTL_fwd.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_NoCase.h>
#include <N_TOP_Misc.h>

// ---------- Forward Declarations ----------

class N_TOP_CktGraph;
class N_TOP_CktGraphCreator;
class N_TOP_CktNodeCreator;
class N_TOP_CktNode;
class N_TOP_NodeDevBlock;
class N_TOP_NodeBlock;
class N_TOP_TopoLSUtil;
class N_TOP_Directory;

class N_ANP_AnalysisInterface;

class N_PDS_Manager;

class N_IO_RestartNode;
class N_IO_CmdParse;
class N_IO_PkgOptionsMgr;

//-----------------------------------------------------------------------------
// Class         : N_TOP_Topology
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class N_TOP_Topology
{

public:

  // Default constructor.
  N_TOP_Topology(N_IO_CmdParse & cp);

  // Constructor.
  N_TOP_Topology(N_DEV_DeviceInterface * devInt, N_IO_CmdParse & cp);

  // Constructor.
  N_TOP_Topology(
      N_DEV_DeviceInterface * devInt,
      const string & graphType,
      N_IO_CmdParse & cp);

  // Destructor
  ~N_TOP_Topology();

  // Registration functions:

  // Register the pointer to the device interface.
  bool registerDeviceInterface(N_DEV_DeviceInterface * devInt);
  // Register the pointer to the parallel services manager.
  bool registerParallelMgr(N_PDS_Manager * pdsmgr);
  // Register the pointer to the time-integration manager.
  bool registeranaInt(N_ANP_AnalysisInterface * anaInt);
  // Method to register the package options manager
  bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

  bool registerICs(const N_UTL_OptionBlock & ob);

  // Setup linear system information.
  bool setupGlobalIndices();
  bool setupGlobalAccessors();

  // Current add methods for voltage nodes tailored to work with device
  // ghosting.
  void addVoltageNode(const N_TOP_NodeBlock & nb);
  // Current add methods for devices tailored to work with device ghosting.
  void addDevice(const N_TOP_NodeBlock & nb, const Teuchos::RefCountPtr<N_DEV_InstanceBlock> ibPtr);

  // Main call to register global id's with device instances.
  void registerGIDswithDevs();

  // Main call to register state-variable global id's with device instances.
  void registerStateGIDswithDevs();

  // Main call to register store-variable global id's with device instances.
  void registerStoreGIDswithDevs();

  // Main call to register local id's with device instances.
  void registerLIDswithDevs();

  // Main call to register local id's with device instances.
  void registerJacLIDswithDevs();

  // Resolution of Secondary Dependencies for devices.
  void resolveDependentVars();

  // Output to screen, debug, traversals of graph.
  void OutputBFSGraphLists();
  void OutputDFSGraphLists();

  // Generate Ordered node list using BFS.
  void setOrderedNodeList() const;

#ifdef Xyce_PARALLEL_MPI
  // merge the off processor superNodeList and communicate the same list
  // to all procs so topology reduction is the same on all procs
  void mergeOffProcTaggedNodesAndDevices();
#endif

  // this functions builds up the supernode list.  Called at the end of
  // parsing from n_cir_xyce
  void verifyNodesAndDevices();

 // remove any nodes and devices that were tagged for removal during parsing
  void removeTaggedNodesAndDevices();

  // Used for delayed instantiation of devices
  void instantiateDevices();

  // Use ordered node list to generate ordered global id list.
  void returnNodeGIDVec(vector<int> & nodeGIDVec);

  // Generate ordered externnode global id vector using orderedNodeList.
  void returnExternNodeGIDVec(vector< pair<int,int> > & nodeGIDVec);

  // Generate ordered soln var global id vector using orderedNodeList.
  void returnSVarGIDVec(vector<int> & sVarGIDVec);

  // Generate ordered soln var global id vector using orderedNodeList for
  // external nodes.
  void returnExternSVarGIDVec(vector< pair<int,int> > & sVarGIDVec);

  // Generate ordered state var global id vector using orderedNodeList.
  void returnStateVarGIDVec(vector<int> & sVarGIDVec);

  // Generate ordered store var global id vector using orderedNodeList.
  void returnStoreVarGIDVec(vector<int> & sVarGIDVec);

  // Generate ordered state var global id vector using orderedNodeList for
  // external nodes.
  void returnExternStateVarGIDVec(vector< pair<int,int> > & sVarGIDVec);

  // Generate ordered store var global id vector using orderedNodeList for
  // external nodes.
  void returnExternStoreVarGIDVec(vector< pair<int,int> > & sVarGIDVec);

  // Generate ordered list of variable types
  void returnVarTypeVec(vector<char> & varTypeVec) const;

  // Generate list of voltage node sVar GIDs
  void returnSVarVNodeGIDVec(vector<int> & sVarVNodeGIDVec);

  // Generate list of vsrc sVar GIDs
  void returnSVarVsrcGIDVec(vector<int> & sVarVsrcGIDVec);

  // Generate list of node voltage IDs that have no DC path to ground
  // NOTE:  The IDs returned from this method are ONLY related to _VNODE.
  void returnSVarNoDCPathIDVec( vector<string> & SVarNoDCPathGIDVec );

  // Generate list of node voltage IDs that are connected to only one device
  // terminal.
  // NOTE:  The IDs returned from this method are ONLY related to _VNODE.
  void returnSVarConnToOneTermIDVec( vector<string> & SVarConnToOneTermGIDVec );

  // Return list of solution var indices for named node.  Returns false if node
  // not owned or not local.
  bool getNodeSVarGIDs( const NodeID& id,
                        list<int> & sVarGIDList,
                        list<int> & extSVarGIDList,
                        char & type );

  // Accessor to get utility class for linear solver.
  N_TOP_TopoLSUtil * get_LinSolvUtil();

  // Migrate Functionality:

  // Extracts nodes and devices to be migrated.
  void extractMigrateNodes(const int & num, const int * nodeGIDs,
                           const int * procs);

  // Inserts nodes and devices from migration Map.
  void insertMigrateNodes(const int & num, const int * nodeGIDs,
                          const int * procs);

  // Clears migration map.
  void clearMigrateNodeMap();

  // Get node from migration node map.
  map < int, N_TOP_NodeDevBlock * > & getMigrateNode(const int & id);

  // Add node to migration node map.
  void addMigrateNode(const int & id, map < int,
                      N_TOP_NodeDevBlock * > & ndbL);

  // Prune out unnecessary circuit nodes
  void pruneCkt(const int & num, const int * nodeGIDs);

  // Prune out unnecessary device nodes
  void pruneDevNode(const int & node);

  void regenerateGIDNodeMap();

  bool generateICLoader();

  // Restart capability.
  bool getRestartNodes(vector < N_IO_RestartNode * > & nodeVec);
  bool restoreRestartNodes(vector < N_IO_RestartNode * > & nodeVec);

  // Directory generation.
  bool generateDirectory();

#ifdef Xyce_TEST_SOLN_VAR_MAP
  // Solution variable name output.
  bool outputNameFile();
#endif

  // These calls generate and return node maps.
  // NOTE:  The node type is known, this map will have unique keys.
  void getNodeNames (map<string, pair<int, double>, Xyce::LessNoCase > &);
  void getStateNodeNames (map<string, pair<int, double>, Xyce::LessNoCase > &);
  void getStoreNodeNames (map<string, pair<int, double>, Xyce::LessNoCase > &);
  void getExternNodeNames (map<string, pair<int, double>, Xyce::LessNoCase > &);
  void getVsrcNodes (set<string> & vsrcSet);

  // get node ids, names, and types
  bool getRawData( map< int, string > & nRef, vector< char > & tRef );

  N_ANP_AnalysisInterface * getTIAManager() const { return anaIntPtr_; }

  // These functions are added to augment a copy of the netlist file to
  // include resistors which connect ground to "dangling" nodes.
  // NOTE: These IDs are assumed to be _VNODE, so no need to use NodeID.

  void addResistors(const vector<string> & inputVec, const string & resValue,
		    const string & netlistFile, bool oneTermNotNoDCPath);

  void appendEndStatement(const string & netlistFile);

private:

  // Don't allow copy construction or assignment.

  // Copy constructor (private)
  N_TOP_Topology(const N_TOP_Topology & right);

  // Assignment operator (private)
  N_TOP_Topology & operator = (const N_TOP_Topology & right);

private:

  // command line object
  N_IO_CmdParse & commandLine_;

  // maximum number of tries to find the graph center.
  int maxTries_;

  // Utility class to extract alloc info for linear solver.
  N_TOP_TopoLSUtil * lsUtilPtr_;
  N_TOP_Directory * dirPtr_;
  friend class N_TOP_TopoLSUtil;
  friend class N_TOP_Directory;

  list < N_TOP_CktGraph * > graphList_;

  N_TOP_CktGraph * mainGraphPtr_;

  N_TOP_CktGraphCreator * graphCreatorPtr_;

  N_TOP_CktNodeCreator * nodeCreatorPtr_;

  N_DEV_DeviceInterface * devIntPtr_;
  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;

  // Pointer to the parallel services manager object.
  N_PDS_Manager * pdsMgrPtr_;

  // Pointer to the time-integration manager object.
  N_ANP_AnalysisInterface * anaIntPtr_;

  mutable list < N_TOP_CktNode * > * orderedNodeListPtr_;
  list < N_TOP_CktNode * >  deviceNodeListPtr_;

  map < string, N_DEV_InstanceBlock * > devInstBlockMap_;
  map < string, int > devInstMap_;

  map < int, map < int, N_TOP_NodeDevBlock * > > migrateNodeMap_;

  N_UTL_OptionBlock * icSettings_;

  map < int, int > depSolnGIDMap_;
  map < int, int > depStateGIDMap_;
  map < int, int > depStoreGIDMap_;

  // this is a list of nodes to be supernoded
  // format is: nodeToBeReplaced, nodeToBeKept in each pair.
  list < pair<NodeID, NodeID> > superNodeList;
#ifdef Xyce_PARALLEL_MPI
  list < pair<NodeID, NodeID> > globalSuperNodeList;
#endif

  friend ostream & operator << (ostream & os, const N_TOP_Topology & topo);
};

//-----------------------------------------------------------------------------
// Function      : N_TOP_Topology::get_LinSolvUtil
// Purpose       : Accessor to get utility for linear solver data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
inline N_TOP_TopoLSUtil * N_TOP_Topology::get_LinSolvUtil()
{
  return lsUtilPtr_;
}

#endif
