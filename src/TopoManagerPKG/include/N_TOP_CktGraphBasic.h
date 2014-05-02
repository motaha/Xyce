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
// Filename       : $RCSfile: N_TOP_CktGraphBasic.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Electrical & Microsystems
//
// Creation Date  : 08/10/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktGraphBasic_h
#define N_TOP_CktGraphBasic_h 1

#include <iosfwd>
#include <string>
#include <map>

#include <N_TOP_fwd.h>
#include <N_TOP_CktGraph.h>
#include <N_UTL_Graph.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktGraphBasic
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Electrical & Microsystems
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
class CktGraphBasic : public CktGraph
{

public:

  // Default constructor.
  CktGraphBasic(const int maxTries);

  // Constructor
  CktGraphBasic(const std::string & cgID, const int maxTries);

  // Constructor
  CktGraphBasic(const std::string & cgID, const std::list<NodeID> & nodelist, const int maxTries);

  // Copy constructor.
  CktGraphBasic(const CktGraphBasic & right);

  // Destructor
  ~CktGraphBasic();

  // Inserts graph node for ckt node if does not exist
  void InsertNode(CktNode * cktnode,
                  const std::list<NodeID> & neighborList);

  // Destroys graph node and returns the ckt node
  CktNode * ExtractNode(const NodeID & cnID);
  // Destroys graph node and returns the ckt node
  CktNode * ExtractNode(const int & cnID);

  // Returns pointer to specified ckt node
  CktNode * FindCktNode(const NodeID & cnID);
  // Returns pointer to specified ckt node
  CktNode * FindCktNode(const int & cnID);

  // Breadth first, depth first, and basic traversals of graph
  //----------------------------------------------------------

  // Produces list of ckt nodes in breadth-first traversal order
  std::list<CktNode*>* getBFSNodeList();
  // Produces list of ckt nodes in depth-first traversal order
  std::list<CktNode*>* getDFSNodeList();
  
  // Returns the node list from the circuit graph without any specific ordering
  const std::map<NodeID, CktNode *> getNodeList() { return cktgph.getData1Map(); }

  // Returns the number of nodes
  int returnNumNodes();

  // production of global and local id lists for load balance and solver
  //----------------------------------------------------------

  // Loop through ordered node list and generate ordered list of global id's
  void returnGIDList(std::list<int>& gidList);

  // Loop through node list and produce ordered list of global id's for soln
  void returnSVGIDList(std::list<int>& svGIDList);

  // Loop through node list and produce ordered list or local id's (strings)
  // for nodes
  void returnLIDList(std::list<NodeID>& lidList);

  // assigns global ids to device int. and ext. vars
  //----------------------------------------------------------

  // Loop over nodes and register int and ext global ids with each device
  void registerGIDswithDevs();

  // Loop over nodes and register int and ext global ids with each device
  void registerStateGIDswithDevs();
  void registerStoreGIDswithDevs();

  // Loop over nodes and register int and ext local ids with each device
  void registerLIDswithDevs( Indexor & indexor );
  // Loop over nodes and register state local ids with each device
  void registerStateLIDswithDevs( Indexor & indexor );
  void registerStoreLIDswithDevs( Indexor & indexor );

  // Loop over nodes and register dep. local ids with each device
  void registerDepLIDswithDevs( Indexor & indexor );
  // Loop over nodes and register dep. state local ids with each device
  void registerDepStateLIDswithDevs( Indexor & indexor );
  void registerDepStoreLIDswithDevs( Indexor & indexor );

  //Loop over nodes and register jacobian offsets with each device
  void registerJacLIDswithDevs( Indexor & indexor );

  // methods used to build Zoltan Query Functions
  //----------------------------------------------------------

  // Return number of connected edges for given node
  int returnNumEdges(const NodeID & id);

  // Return number of connected edges for given node
  int returnNumEdges(const int & globalID);

  // Returns vector of adj ids
  void returnAdjIDs( const NodeID & id, std::vector<NodeID> & adj_ids );

  // Returns vector of adj gids
  void returnAdjGIDs( int gid, std::vector<int> & adj_gids );

  // Creates and returns a NodeDevBlock for the specified GID
  NodeDevBlock * returnNodeDevBlock( int gid );

  // Loop over adjacent nodes creating ordered lists of neighboring global id's
  // and owning processor numbers
  void returnAdjNodes(const NodeID & id, std::list<int> & gidList,
                      std::list<int> & svGIDList, std::list<int> & procList,
                      std::list<NodeID> & idList);

  // Loop over adjacent nodes creating ordered lists of neighboring global id's
  // and owning processor numbers
  void returnAdjNodes(const int & globalID, std::list<int> & gidList,
                      std::list<int> & svGIDList, std::list<int> & procList,
                      std::list<NodeID> & idList);

  // Loop over adjacent nodes connected to ground creating ordered lists of
  // neighboring global id's and owning processor numbers 
  void returnAdjNodesWithGround(const int & globalID, std::list<int> & gidList,
  	std::list<int> & svGIDList, std::list<int> & procList,
  	std::list<NodeID> & idList);

  // Redo global index node map
  void regenerateGIDNodeMap();
  
  // supernode given nodes
  CktNode * replaceNode( const NodeID nodeToBeReplaced, const NodeID nodeToKeep );
  void removeRedundantDevices(std::vector< CktNode * > & removedDevices);

private:

  // circuit graph.
  // pair = <id, node type>, int = gid, CktNode = data
    N_UTL_Graph<NodeID,int,CktNode*> cktgph;

  // List of ckt nodes in breadth-first traversal order.
  std::list<CktNode*> BFSNodeList_;
  // List of ckt nodes in depth-first traversal order.
  std::list<CktNode*> DFSNodeList_;

  // Don't update traversals if not modified - flag.
  bool isModified_;

  // Maximum number of attempts to compute the graph center.
  int maxTries_;

public:
  std::ostream & put(std::ostream & os) const;
};

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::CktGraphBasic N_TOP_CktGraphBasic;

#endif
