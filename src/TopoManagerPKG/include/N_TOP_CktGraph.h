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
// Revision Number: $Revision: 1.28.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:50 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktGraph_h
#define N_TOP_CktGraph_h 1

// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>
#include <N_TOP_Misc.h>

#include <string>
#ifdef HAVE_OSTREAM
#include <ostream>
#else
#ifdef HAVE_OSTREAM_H
#include <ostream.h>
#else
#error Neither ostream nor ostream.h detected
#endif
#endif

#include <list>
#include <vector>
#include <map>

// ------------ Xyce Includes ------------
#include <N_UTL_Xyce.h>

// ---------- Forward Declarations ----------
class N_TOP_CktNode;
class N_TOP_Topology;

class N_TOP_NodeDevBlock;

class N_TOP_Indexor;

//-----------------------------------------------------------------------------
// Class         : N_TOP_CktGraph
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class N_TOP_CktGraph
{

public:

  // Default constructor.
  N_TOP_CktGraph() {}

  // Constructor
  N_TOP_CktGraph(const N_TOP_CktGraph& right) {}

  // Constructor
  N_TOP_CktGraph(const string& cgID)
  : id_(cgID)
  {}

  // Constructor
  N_TOP_CktGraph(const string& cgID, const list<NodeID>& nodelist)
  : id_(cgID), internNodeList_(nodelist)
  {}

  // Destructor
  virtual ~N_TOP_CktGraph() { }

  // Inserts graph node for ckt node if does not exist (abstract).
  virtual void InsertNode(N_TOP_CktNode* node,
                          const list<NodeID>& neighborList) = 0;

  // Destroys graph node and returns the ckt node (abstract).
  virtual N_TOP_CktNode* ExtractNode(const NodeID& cgID) = 0;
  // Destroys graph node and returns the ckt node (abstract).
  virtual N_TOP_CktNode* ExtractNode(const int& cgID) = 0;

  // Returns pointer to specified ckt node (abstract).
  virtual N_TOP_CktNode* FindCktNode(const NodeID& cnID) = 0;
  // Returns pointer to specified ckt node (abstract).
  virtual N_TOP_CktNode* FindCktNode(const int& cnID) = 0;

  const string& getParentGraph() const { return parentGraph_; }

  const list<NodeID>& getExternalNodeList() const { return externNodeList_; }
  const list<NodeID>& getInternalNodeList() const { return internNodeList_; }

  // Produces list of ckt nodes in breadth-first traversal order (abstract).
  virtual list<N_TOP_CktNode*>* getBFSNodeList() = 0;

  // Produces list of ckt nodes in depth-first traversal order (abstract).
  virtual list<N_TOP_CktNode*>* getDFSNodeList() = 0;

  // Returns the node list from the circuit graph without any specific ordering
  virtual const map<NodeID, N_TOP_CktNode *> getNodeList() = 0;

  // Loop over nodes and register int and ext global ids with each device
  // (abstract).
  virtual void registerGIDswithDevs() = 0;

  // Loop over nodes and register int and ext global ids with each device
  // (abstract).
  virtual void registerStateGIDswithDevs() = 0;
  virtual void registerStoreGIDswithDevs() = 0;

  // Loop over nodes and register int and ext global ids with each device
  virtual void registerLIDswithDevs( N_TOP_Indexor & indexor ) = 0;
  // Loop over nodes and register state local ids with each device
  virtual void registerStateLIDswithDevs( N_TOP_Indexor & indexor ) = 0;
  virtual void registerStoreLIDswithDevs( N_TOP_Indexor & indexor ) = 0;

  // Loop over nodes and register dep. global ids with each device
  virtual void registerDepLIDswithDevs( N_TOP_Indexor & indexor ) = 0;
  // Loop over nodes and register dep. state local ids with each device
  virtual void registerDepStateLIDswithDevs( N_TOP_Indexor & indexor ) = 0;
  virtual void registerDepStoreLIDswithDevs( N_TOP_Indexor & indexor ) = 0;

  // Registration of jacobian local offsets with devices
  virtual void registerJacLIDswithDevs( N_TOP_Indexor & indexor ) = 0;

  // Returns the number of nodes (abstract).
  virtual int returnNumNodes() = 0;

  // Loop through ordered node list and generate ordered list of global id's
  // (abstract).
  virtual void returnGIDList(list<int>& gidList) = 0;

  // Loop through node list and produce ordered list of global id's for soln
  // var's (abstract).
  virtual void returnSVGIDList(list<int>& svGIDList) = 0;

  // Loop through node list and produce ordered list or local id's (strings)
  // for nodes (abstract).
  virtual void returnLIDList(list<NodeID>& lidList) = 0;

  // Return number of connected edges for given node (abstract).
  virtual int returnNumEdges(const NodeID& id) = 0;

  // Return number of connected edges for given node (Global) (abstract).
  virtual int returnNumEdges(const int& globalID) = 0;

  // Returns adj IDs to the given ID
  virtual void returnAdjIDs( const NodeID& id, vector<NodeID>& adj_ids ) = 0;

  // Returns adj GIDs to the given GID
  virtual void returnAdjGIDs( int gid, vector<int>& adj_gids ) = 0;

  // Create and return a N_TOP_NodeDevBlock for the given GID
  virtual N_TOP_NodeDevBlock* returnNodeDevBlock( int gid ) = 0;

  // Loop over adjacent nodes creating ordered lists of neighboring global id's
  // and owning processor numbers (abstract).
  virtual void returnAdjNodes(const NodeID& id,
                              list<int>& gidList,
                              list<int>& svGIDList,
                              list<int>& procList,
                              list<NodeID>& idList) = 0;

  // Loop over adjacent nodes creating ordered lists of neighboring global id's
  // and owning processor numbers (abstract).
  virtual void returnAdjNodes(const int& globalID,
                              list<int>& gidList,
                              list<int>& svGIDList,
                              list<int>& procList,
                              list<NodeID>& idList) = 0;

  // Loop over adjacent nodes connected to ground creating ordered lists of
  // neighboring global id's and owning processor numbers (abstract).
  virtual void returnAdjNodesWithGround(const int& globalID,
                                        list<int>& gidList,
                                        list<int>& svGIDList,
                                        list<int>& procList,
                                        list<NodeID >& idList) = 0;

  // Redo global index node map (abstract).
  virtual void regenerateGIDNodeMap() = 0;
  
  // supernode given nodes
  virtual N_TOP_CktNode * replaceNode( const NodeID nodeToBeReplaced, const NodeID nodeToKeep ) = 0;
  virtual void removeRedundantDevices(vector< N_TOP_CktNode * > & removedDevices) = 0;

  const string& get_id() const { return id_; }

protected:

  void setParentGraph(const string & graphname) { parentGraph_ = graphname; }
  void setExternalNodeList(const list<NodeID>& nodelist) { externNodeList_ = nodelist; }
  void setInternalNodeList(const list<NodeID>& nodelist) { internNodeList_ = nodelist; }

protected:

  void set_id(string & value) { id_ = value; }

protected:

  string id_;

  string parentGraph_;
  string templateGraph_;

  list<NodeID > externNodeList_;
  list<NodeID > internNodeList_;

  virtual ostream & put(ostream & os) const = 0;

  friend ostream& operator<<(ostream& os, const N_TOP_CktGraph& cg);

};

#endif
