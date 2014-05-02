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
// Filename       : $RCSfile: N_TOP_CktGraph.h,v $
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
// Revision Number: $Revision: 1.33 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktGraph_h
#define N_TOP_CktGraph_h 1

#include <N_TOP_Misc.h>

#include <string>
#include <ostream>
#include <list>
#include <vector>
#include <map>

#include <N_TOP_fwd.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktGraph
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktGraph
{
public:

  // Default constructor.
  CktGraph() {}

  // Constructor
  CktGraph(const CktGraph& right) {}

  // Constructor
  CktGraph(const std::string& cgID)
  : id_(cgID)
  {}

  // Constructor
  CktGraph(const std::string& cgID, const std::list<NodeID>& nodelist)
  : id_(cgID), internNodeList_(nodelist)
  {}

  // Destructor
  virtual ~CktGraph() { }

  // Inserts graph node for ckt node if does not exist (abstract).
  virtual void InsertNode(CktNode* node,
                          const std::list<NodeID>& neighborList) = 0;

  // Destroys graph node and returns the ckt node (abstract).
  virtual CktNode* ExtractNode(const NodeID& cgID) = 0;
  // Destroys graph node and returns the ckt node (abstract).
  virtual CktNode* ExtractNode(const int& cgID) = 0;

  // Returns pointer to specified ckt node (abstract).
  virtual CktNode* FindCktNode(const NodeID& cnID) = 0;
  // Returns pointer to specified ckt node (abstract).
  virtual CktNode* FindCktNode(const int& cnID) = 0;

  const std::string& getParentGraph() const { return parentGraph_; }

  const std::list<NodeID>& getExternalNodeList() const { return externNodeList_; }
  const std::list<NodeID>& getInternalNodeList() const { return internNodeList_; }

  // Produces list of ckt nodes in breadth-first traversal order (abstract).
  virtual std::list<CktNode*>* getBFSNodeList() = 0;

  // Produces list of ckt nodes in depth-first traversal order (abstract).
  virtual std::list<CktNode*>* getDFSNodeList() = 0;

  // Returns the node list from the circuit graph without any specific ordering
  virtual const std::map<NodeID, CktNode *> getNodeList() = 0;

  // Loop over nodes and register int and ext global ids with each device
  // (abstract).
  virtual void registerGIDswithDevs() = 0;

  // Loop over nodes and register int and ext global ids with each device
  // (abstract).
  virtual void registerStateGIDswithDevs() = 0;
  virtual void registerStoreGIDswithDevs() = 0;

  // Loop over nodes and register int and ext global ids with each device
  virtual void registerLIDswithDevs( Indexor & indexor ) = 0;
  // Loop over nodes and register state local ids with each device
  virtual void registerStateLIDswithDevs( Indexor & indexor ) = 0;
  virtual void registerStoreLIDswithDevs( Indexor & indexor ) = 0;

  // Loop over nodes and register dep. global ids with each device
  virtual void registerDepLIDswithDevs( Indexor & indexor ) = 0;
  // Loop over nodes and register dep. state local ids with each device
  virtual void registerDepStateLIDswithDevs( Indexor & indexor ) = 0;
  virtual void registerDepStoreLIDswithDevs( Indexor & indexor ) = 0;

  // Registration of jacobian local offsets with devices
  virtual void registerJacLIDswithDevs( Indexor & indexor ) = 0;

  // Returns the number of nodes (abstract).
  virtual int returnNumNodes() = 0;

  // Loop through ordered node list and generate ordered list of global id's
  // (abstract).
  virtual void returnGIDList(std::list<int>& gidList) = 0;

  // Loop through node list and produce ordered list of global id's for soln
  // var's (abstract).
  virtual void returnSVGIDList(std::list<int>& svGIDList) = 0;

  // Loop through node list and produce ordered list or local id's (std::strings)
  // for nodes (abstract).
  virtual void returnLIDList(std::list<NodeID>& lidList) = 0;

  // Return number of connected edges for given node (abstract).
  virtual int returnNumEdges(const NodeID& id) = 0;

  // Return number of connected edges for given node (Global) (abstract).
  virtual int returnNumEdges(const int& globalID) = 0;

  // Returns adj IDs to the given ID
  virtual void returnAdjIDs( const NodeID& id, std::vector<NodeID>& adj_ids ) = 0;

  // Returns adj GIDs to the given GID
  virtual void returnAdjGIDs( int gid, std::vector<int>& adj_gids ) = 0;

  // Create and return a NodeDevBlock for the given GID
  virtual NodeDevBlock* returnNodeDevBlock( int gid ) = 0;

  // Loop over adjacent nodes creating ordered lists of neighboring global id's
  // and owning processor numbers (abstract).
  virtual void returnAdjNodes(const NodeID& id,
                              std::list<int>& gidList,
                              std::list<int>& svGIDList,
                              std::list<int>& procList,
                              std::list<NodeID>& idList) = 0;

  // Loop over adjacent nodes creating ordered lists of neighboring global id's
  // and owning processor numbers (abstract).
  virtual void returnAdjNodes(const int& globalID,
                              std::list<int>& gidList,
                              std::list<int>& svGIDList,
                              std::list<int>& procList,
                              std::list<NodeID>& idList) = 0;

  // Loop over adjacent nodes connected to ground creating ordered lists of
  // neighboring global id's and owning processor numbers (abstract).
  virtual void returnAdjNodesWithGround(const int& globalID,
                                        std::list<int>& gidList,
                                        std::list<int>& svGIDList,
                                        std::list<int>& procList,
                                        std::list<NodeID >& idList) = 0;

  // Redo global index node map (abstract).
  virtual void regenerateGIDNodeMap() = 0;
  
  // supernode given nodes
  virtual CktNode * replaceNode( const NodeID nodeToBeReplaced, const NodeID nodeToKeep ) = 0;
  virtual void removeRedundantDevices(std::vector< CktNode * > & removedDevices) = 0;

  const std::string& get_id() const { return id_; }

protected:

  void setParentGraph(const std::string & graphname) { parentGraph_ = graphname; }
  void setExternalNodeList(const std::list<NodeID>& nodelist) { externNodeList_ = nodelist; }
  void setInternalNodeList(const std::list<NodeID>& nodelist) { internNodeList_ = nodelist; }

protected:

  void set_id(std::string & value) { id_ = value; }

protected:

  std::string id_;

  std::string parentGraph_;
  std::string templateGraph_;

  std::list<NodeID > externNodeList_;
  std::list<NodeID > internNodeList_;

  virtual std::ostream & put(std::ostream & os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const CktGraph& cg);
};

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::CktNode N_TOP_CktNode;

#endif
