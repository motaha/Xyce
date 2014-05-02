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
// Revision Number: $Revision: 1.72.2.2 $
//
// Revision Date  : $Date: 2014/03/06 17:23:42 $
//
// Current Owner  : $Author: erkeite $
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

#include <N_TOP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_NoCase.h>
#include <N_TOP_Misc.h>

class N_PDS_Manager;

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : Topology
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class Topology
{

public:

  // Default constructor.
  Topology(IO::CmdParse & cp);

  // Constructor.
  Topology(Device::DeviceInterface * devInt, IO::CmdParse & cp);

  // Constructor.
  Topology(
      Device::DeviceInterface * devInt,
      const std::string & graphType,
      IO::CmdParse & cp);

  // Destructor
  ~Topology();

  // Registration functions:

  // Register the pointer to the device interface.
  bool registerDeviceInterface(Device::DeviceInterface * devInt);
  // Register the pointer to the parallel services manager.
  bool registerParallelMgr(N_PDS_Manager * pdsmgr);
  // Register the pointer to the time-integration manager.
  bool registeranaInt(Xyce::Analysis::AnalysisInterface * anaInt);
  // Method to register the package options manager
  bool registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr );

  bool registerICs(const N_UTL_OptionBlock & ob);

  // Setup linear system information.
  bool setupGlobalIndices();
  bool setupGlobalAccessors();

  // Current add methods for voltage nodes tailored to work with device
  // ghosting.
  void addVoltageNode(const NodeBlock & nb);
  // Current add methods for devices tailored to work with device ghosting.
  void addDevice(const NodeBlock & nb, const Teuchos::RefCountPtr<Device::InstanceBlock> ibPtr);

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
  void returnNodeGIDVec(std::vector<int> & nodeGIDVec);

  // Generate ordered externnode global id vector using orderedNodeList.
    void returnExternNodeGIDVec(std::vector<std::pair<int,int> > & nodeGIDVec);

  // Generate ordered soln var global id vector using orderedNodeList.
  void returnSVarGIDVec(std::vector<int> & sVarGIDVec);

  // Generate ordered soln var global id vector using orderedNodeList for
  // external nodes.
  void returnExternSVarGIDVec(std::vector< std::pair<int,int> > & sVarGIDVec);

  // Generate ordered state var global id vector using orderedNodeList.
  void returnStateVarGIDVec(std::vector<int> & sVarGIDVec);

  // Generate ordered store var global id vector using orderedNodeList.
  void returnStoreVarGIDVec(std::vector<int> & sVarGIDVec);

  // Generate ordered state var global id vector using orderedNodeList for
  // external nodes.
  void returnExternStateVarGIDVec(std::vector< std::pair<int,int> > & sVarGIDVec);

  // Generate ordered store var global id vector using orderedNodeList for
  // external nodes.
  void returnExternStoreVarGIDVec(std::vector< std::pair<int,int> > & sVarGIDVec);

  // Generate ordered list of variable types
  void returnVarTypeVec(std::vector<char> & varTypeVec) const;

  // Generate list of voltage node sVar GIDs
  void returnSVarVNodeGIDVec(std::vector<int> & sVarVNodeGIDVec);

  // Generate list of vsrc sVar GIDs
  void returnSVarVsrcGIDVec(std::vector<int> & sVarVsrcGIDVec);

  // Generate list of node voltage IDs that have no DC path to ground
  // NOTE:  The IDs returned from this method are ONLY related to _VNODE.
  void returnSVarNoDCPathIDVec( std::vector<std::string> & SVarNoDCPathGIDVec );

  // Generate list of node voltage IDs that are connected to only one device
  // terminal.
  // NOTE:  The IDs returned from this method are ONLY related to _VNODE.
  void returnSVarConnToOneTermIDVec( std::vector<std::string> & SVarConnToOneTermGIDVec );

  // Return list of solution var indices for named node.  Returns false if node
  // not owned or not local.
  bool getNodeSVarGIDs( const NodeID& id,
                        std::list<int> & sVarGIDList,
                        std::list<int> & extSVarGIDList,
                        char & type );

  // Accessor to get utility class for linear solver.
  TopoLSUtil * get_LinSolvUtil();

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
  std::map< int, NodeDevBlock * > & getMigrateNode(const int & id);

  // Add node to migration node map.
  void addMigrateNode(const int & id, std::map< int,
                      NodeDevBlock * > & ndbL);

  // Prune out unnecessary circuit nodes
  void pruneCkt(const int & num, const int * nodeGIDs);

  // Prune out unnecessary device nodes
  void pruneDevNode(const int & node);

  void regenerateGIDNodeMap();

  bool generateICLoader();

  // Restart capability.
  bool getRestartNodes(std::vector< IO::RestartNode * > & nodeVec);
  bool restoreRestartNodes(std::vector< IO::RestartNode * > & nodeVec);

  // Directory generation.
  bool generateDirectory();

  // Solution variable name output.
  bool outputNameFile(bool overRideOutput=false);

  // These calls generate and return node maps.
  // NOTE:  The node type is known, this map will have unique keys.
  void getNodeNames (std::map<std::string, std::pair<int, double>, Xyce::LessNoCase > &);
  void getStateNodeNames (std::map<std::string, std::pair<int, double>, Xyce::LessNoCase > &);
  void getStoreNodeNames (std::map<std::string, std::pair<int, double>, Xyce::LessNoCase > &);
  void getExternNodeNames (std::map<std::string, std::pair<int, double>, Xyce::LessNoCase > &);
    void getVsrcNodes (std::set<std::string> & vsrcSet);

  // get node ids, names, and types
  bool getRawData( std::map< int, std::string > & nRef, std::vector< char > & tRef );

  Xyce::Analysis::AnalysisInterface * getTIAManager() const { return anaIntPtr_; }

  // These functions are added to augment a copy of the netlist file to
  // include resistors which connect ground to "dangling" nodes.
  // NOTE: These IDs are assumed to be _VNODE, so no need to use NodeID.

  void addResistors(const std::vector<std::string> & inputVec, const std::string & resValue,
		    const std::string & netlistFile, bool oneTermNotNoDCPath);

  void appendEndStatement(const std::string & netlistFile);

private:

  // Don't allow copy construction or assignment.

  // Copy constructor (private)
  Topology(const Topology & right);

  // Assignment operator (private)
  Topology & operator = (const Topology & right);

private:

  // command line object
  IO::CmdParse & commandLine_;

  // maximum number of tries to find the graph center.
  int maxTries_;

  // Utility class to extract alloc info for linear solver.
  TopoLSUtil * lsUtilPtr_;
  Directory * dirPtr_;
  friend class TopoLSUtil;
  friend class Directory;

  std::list< CktGraph * > graphList_;

  CktGraph * mainGraphPtr_;

  CktGraphCreator * graphCreatorPtr_;

  CktNodeCreator * nodeCreatorPtr_;

  Device::DeviceInterface * devIntPtr_;
  // package options manager
  IO::PkgOptionsMgr *pkgOptMgrPtr_;

  // Pointer to the parallel services manager object.
  N_PDS_Manager * pdsMgrPtr_;

  // Pointer to the time-integration manager object.
  Xyce::Analysis::AnalysisInterface * anaIntPtr_;

  mutable std::list< CktNode * > * orderedNodeListPtr_;
  std::list< CktNode * >  deviceNodeListPtr_;

  std::map< std::string, Device::InstanceBlock * > devInstBlockMap_;
  std::map< std::string, int > devInstMap_;

  std::map< int, std::map< int, NodeDevBlock * > > migrateNodeMap_;

  N_UTL_OptionBlock * icSettings_;

  std::map< int, int > depSolnGIDMap_;
  std::map< int, int > depStateGIDMap_;
  std::map< int, int > depStoreGIDMap_;

  // this is a list of nodes to be supernoded
  // format is: nodeToBeReplaced, nodeToBeKept in each pair.
  std::list< std::pair<NodeID, NodeID> > superNodeList;
#ifdef Xyce_PARALLEL_MPI
  std::list< std::pair<NodeID, NodeID> > globalSuperNodeList;
#endif

    friend std::ostream & operator << (std::ostream & os, const Topology & topo);
};

//-----------------------------------------------------------------------------
// Function      : Topology::get_LinSolvUtil
// Purpose       : Accessor to get utility for linear solver data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
inline TopoLSUtil * Topology::get_LinSolvUtil()
{
  return lsUtilPtr_;
}

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::Topology N_TOP_Topology;

#endif
