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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_TOP_CktGraphBasic.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/10/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.28 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

#include <N_UTL_Misc.h>
#include <N_TOP_Misc.h>

#include <N_TOP_CktGraphBasic.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_Dev.h>

#include <N_TOP_Indexor.h>

#include <N_UTL_Functors.h>

#include <N_TOP_NodeDevBlock.h>

#include <N_DEV_DeviceBlock.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::CktGraphBasic
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktGraphBasic::CktGraphBasic(const int maxTries)
 : isModified_(true),
   maxTries_(maxTries)
{
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::CktGraphBasic
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktGraphBasic::CktGraphBasic(const std::string &cgID, 
   const int maxTries)
 : CktGraph(cgID),
   isModified_(true),
   maxTries_(maxTries)
{
}


//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::CktGraphBasic
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktGraphBasic::CktGraphBasic(const std::string & cgID,
   const std::list<NodeID> & nodelist,
   const int maxTries)
 : CktGraph(cgID,nodelist),
   isModified_(true),
   maxTries_(maxTries)
{
}


//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::~CktGraphBasic()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktGraphBasic::~CktGraphBasic()
{
  std::list<CktNode*>::iterator it_nL = BFSNodeList_.begin();
  std::list<CktNode*>::iterator it_nL_end = BFSNodeList_.end();
  for( ; it_nL != it_nL_end; ++it_nL )
    if( *it_nL ) delete *it_nL;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::InsertNode
// Purpose       : Inserts graph node for ckt node if does not exist
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::InsertNode(CktNode* cktnode,
                                     const std::list<NodeID> &neighborList)
{
#ifdef Xyce_DEBUG_TOPOLOGY
 Xyce::dout() << "Inserting Node: " << cktnode->get_id() << ", type = " << cktnode->type() << std::endl;
 Xyce::dout() << *cktnode;
#endif

  std::vector<NodeID> adjVec( neighborList.begin(), neighborList.end() );

  bool inserted = cktgph.insertNode( NodeID(cktnode->get_id(),cktnode->type()), 
                                     cktnode->get_gID(), adjVec, cktnode );
  
  if( !inserted )
  {
    //------- If node already exists, check to see if it should be
    //------- changed to be owned.

    if( cktnode->get_IsOwned() )
    {
      cktgph.getData(NodeID(cktnode->get_id(),cktnode->type()))->set_IsOwned(true);
      cktgph.getData(NodeID(cktnode->get_id(),cktnode->type()))->set_ProcNum(cktnode->get_ProcNum());
      cktgph.getData(NodeID(cktnode->get_id(),cktnode->type()))->set_gID(cktnode->get_gID());
      cktgph.chgKey2(NodeID(cktnode->get_id(),cktnode->type()),cktnode->get_gID());
    }

    delete cktnode;
  }

  //------ Graph has been changed so traversals will change
  isModified_ = true;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::ExtractNode
// Purpose       : Destroys graph node and returns the ckt node
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktNode* CktGraphBasic::ExtractNode (const NodeID &cnID)
{
  // get CktNode*
  CktNode * nodeToBeRemovedCktNodePtr = FindCktNode( cnID );
   
  // Remove nodeToBeReplaced and return it's associated CktNode object
  cktgph.removeKey( cnID );
  isModified_=true;

  return nodeToBeRemovedCktNodePtr;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::ExtractNode
// Purpose       : Destroys graph node and returns the ckt node
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktNode* CktGraphBasic::ExtractNode(const int &cnID)
{
  // get CktNode*
  CktNode * nodeToBeRemovedCktNodePtr = FindCktNode( cnID );
   
  // Remove nodeToBeReplaced and return it's associated CktNode object
  cktgph.removeKey( cnID );
  isModified_=true;

  return nodeToBeRemovedCktNodePtr;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::FindCktNode
// Purpose       : Returns ptr to specified ckt node.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktNode* CktGraphBasic::FindCktNode (const NodeID &cnID)
{
  // If the node type cannot be determined before hand, try to see if there is
  // a unique node associated with the string in cnID.
  if (cnID.second != -1)
  {
    if( cktgph.checkKey(cnID) )
      return cktgph.getData(cnID);
  }
  else
  {
    // Loop through all the possible node types to see if this node exists.
    // Only return it if there is one unique node.  Otherwise, this operation
    // is not well defined.
    int numFound = 0;
    CktNode* foundNode = 0;
    for (int i=0; i<_NUM_NODE_TYPES; ++i)
    {
      NodeID tmpID( cnID.first, i );
      if (cktgph.checkKey(cnID) )
      {
        foundNode = cktgph.getData(cnID);
        numFound++;
      }
    }
    if (numFound == 1)
      return foundNode;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::FindCktNode
// Purpose       : Returns ptr to specified cktnode
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
CktNode* CktGraphBasic::FindCktNode(const int &cnID)
{
  if( cktgph.checkKey(cnID) )
    return cktgph.getData(cnID);

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::getBFSNodeList
// Purpose       : Produces list of ckt nodes in breadth first traversal order
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
std::list<CktNode*> * CktGraphBasic::getBFSNodeList()
{
  //------- Check if list not created or graph modified

  if( BFSNodeList_.empty() || isModified_ )
  {
    //------- Include disconnected parts for graph in traversal
    BFSNodeList_.clear();

    int numNodes = cktgph.numNodes();

    if( numNodes )
    {
      // Compute the graph center if the maximum number of tries to find a good
      // start node is > 1, else just generate a traversal with the first node.
      if (maxTries_ > 1)
      {
        NodeID start_node(cktgph.getCenter( 0.33, maxTries_ ));

        //------- Produce traversal as object in albBFS_
        //------- and extract to BFSNodeList_
        cktgph.generateBFT( start_node );
      }
      else
      {
        cktgph.generateBFT();
      }

      const std::vector<NodeID> & bfsVec = cktgph.getBFT();

      for( int i = 0; i < numNodes; ++i )
        BFSNodeList_.push_front( cktgph.getData(bfsVec[i]) );
    }

    isModified_ = false;
  }

  return &BFSNodeList_;
}


//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::getDFSNodeList
// Purpose       : Produces list of ckt nodes in depth first traversal order
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
std::list<CktNode*> * CktGraphBasic::getDFSNodeList()
{
  N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
                          std::string( "CktGraphBasic::getDFSNodeList() NOT IMPLEMENTED!" ) );
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnNumNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
int CktGraphBasic::returnNumNodes()
{
  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  int cnt = 0;
  std::list<CktNode*>::iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::iterator it_nL_end = tmpNodeList->end();
  for( ; it_nL != it_nL_end; ++it_nL )
    if( (*it_nL)->get_IsOwned() ) ++cnt;

  return cnt;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnGIDList
// Purpose       : Loop through ordered node list and generate ordered
//                 list of global id's
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnGIDList( std::list<int> & gidList )
{
  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::iterator it_nL_end = tmpNodeList->end();

  for( ; it_nL != it_nL_end; ++it_nL )
  {
    if( (*it_nL)->get_IsOwned() )
      gidList.push_back( (*it_nL)->get_gID() );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnSVGIDList
// Purpose       : Loop through node list and produce ordered list
//                 of global id's for soln var's
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnSVGIDList( std::list<int> & svGIDList )
{
  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  for( ; it_nL != it_nL_end; ++it_nL )
  {
    svGIDList.insert( svGIDList.end(), 
                      (*it_nL)->get_SolnVarGIDList().begin(),
                      (*it_nL)->get_SolnVarGIDList().end() );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnLIDList
// Purpose       : Loop through node list and produce ordered list
//                 or local id's (strings) for nodes
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnLIDList( std::list<NodeID> & lidList )
{
  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::iterator it_nL_end = tmpNodeList->end();
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    if( (*it_nL)->get_IsOwned() )
      lidList.push_back( NodeID((*it_nL)->get_id(),(*it_nL)->type()) );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnNumEdges
// Purpose       : Return number of connected edges for given node
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
int CktGraphBasic::returnNumEdges( const NodeID & id )
{
  int cnt=0;

  if( cktgph.checkKey(id) )
  {
    std::vector<NodeID> adjIDs( cktgph.getAdjacent(id) );

    int numAdj = adjIDs.size();
    for( int i = 0; i < numAdj; ++i )
      if( cktgph.getData(adjIDs[i])->get_gID() != -1 ) cnt++;
  }
  else
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
                          std::string( "CktGraphBasic::returnNumEdges(const std::string &) NODE NOT FOUND!" ) );

  return(cnt);
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnNumEdges
// Purpose       : Return number of connected edges for given node (Global)
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
int CktGraphBasic::returnNumEdges(const int & globalID)
{
  int cnt=0;

  if( cktgph.checkKey(globalID) )
  {
    NodeID id(cktgph.getKey1(globalID));
    std::vector<NodeID> adjIDs = cktgph.getAdjacent(id);

    int numAdj = adjIDs.size();
    for( int i = 0; i < numAdj; ++i )
      if( cktgph.getKey2(adjIDs[i]) != -1 ) cnt++;
  }
  else
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
                          std::string( "CktGraphBasic::returnNumEdges(const std::string &) NODE NOT FOUND!" ) );

  return(cnt);
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnAdjIDs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnAdjIDs( const NodeID & id,
                                        std::vector<NodeID> & adj_ids )
{
  adj_ids.clear();

  std::vector<NodeID> adjIDs = cktgph.getAdjacent(id);

  int numAdj = adjIDs.size();
  for( int i = 0; i < numAdj; ++i )
    if( adjIDs[i].first != "0" ) adj_ids.push_back( adjIDs[i] );
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnAdjGIDs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnAdjGIDs( int gid,
                                       std::vector<int> & adj_gids )
{
  adj_gids.clear();

  std::vector<int> adjGIDs = cktgph.getAdjacent(gid);

  int numAdj = adjGIDs.size();
  for( int i = 0; i < numAdj; ++i )
    if( adjGIDs[i] != -1 ) adj_gids.push_back( adjGIDs[i] );
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnNodeDevBlock
// Purpose       : create a NodeDevBlock for the node gid given
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
NodeDevBlock * CktGraphBasic::returnNodeDevBlock( int gid )
{
  if( cktgph.checkKey(gid) )
  {
    NodeBlock nb;

    CktNode * cn = cktgph.getData(gid);
    nb.set_id( cn->get_id() );
    nb.set_gID( cn->get_gID() );
    nb.set_ProcNum( cn->get_ProcNum() );
    nb.set_IsOwned( cn->get_IsOwned() );

    std::list<int> gids;
    std::list<int> sv_gids;
    std::list<int> pids;
    std::list<NodeID> ids;

    returnAdjNodesWithGround( gid, gids, sv_gids, pids, ids );

    std::list<tagged_param> nList;
    std::list<tagged_param> npList;
    std::list<int>::iterator itGL = gids.begin();
    std::list<int>::iterator enGL = gids.end();
    std::list<int>::iterator itPL = pids.begin();
    std::list<NodeID>::iterator itSL = ids.begin();
    for( ; itGL != enGL; ++itGL, ++itPL, ++itSL )
    {
      nList.push_back( tagged_param( (*itSL).first, *itGL ) );
      npList.push_back( tagged_param( (*itSL).first, *itPL ) );
    }

    nb.set_NodeList( nList );
    nb.set_NodeProcList( npList );

    Device::InstanceBlock * ib = 0;
    if( cn->type() ==_DNODE )
    {
      // this seems a little heavy, but NodeDevBlock() is designed to take copies
      // of a node block and instance block in its constructor.  So, mimic that behavior
      // without just passing in the underlying object of the ref counted pointer (if I 
      // did that, then I would end up breaking the ref count pointer).
      Teuchos::RCP< Device::InstanceBlock > ibRcp( (dynamic_cast<CktNode_Dev*>(cn))->devBlock() );
      ib = new Device::InstanceBlock( *ibRcp );
    }
    else
    {
      ib = new Device::InstanceBlock();
    }

    return new NodeDevBlock( nb, *ib );
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnAdjNodes
// Purpose       : Loop over adjacent nodes creating ordered lists
//                 of neighboring global id's and owning proc nums
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnAdjNodes( const NodeID & id,
                                          std::list<int> & gidList,
                                          std::list<int> & svGIDList,
                                          std::list<int> & procList,
                                          std::list<NodeID> & idList )
{
  if( cktgph.checkKey(id) )
  {
    std::vector<NodeID> adjIDs = cktgph.getAdjacent(id);

    int numAdj = adjIDs.size();
    for( int i = 0; i < numAdj; ++i )
    {
      if(cktgph.getKey2(adjIDs[i]) != -1)
      {
        CktNode * cktnode = cktgph.getData(adjIDs[i]);
        gidList.push_back( cktnode->get_gID() );
        idList.push_back( NodeID(cktnode->get_id(),cktnode->type()) );

        svGIDList.insert( svGIDList.end(),
        cktnode->get_SolnVarGIDList().begin(),
        cktnode->get_SolnVarGIDList().end() );

        procList.push_back( cktnode->get_ProcNum() );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnAdjNodes
// Purpose       : Loop over adjacent nodes creating ordered lists
//                 of neighboring global id's and owning proc nums
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnAdjNodes( const int & globalID,
  std::list<int> & gidList,
  std::list<int> & svGIDList,
  std::list<int> & procList,
  std::list<NodeID> & idList )
{
  if( (globalID!=-1) && cktgph.checkKey(globalID) )
  {
    returnAdjNodes( cktgph.getKey1(globalID), gidList, svGIDList, procList, idList );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::returnAdjNodesWithGround
// Purpose       : Loop over adjacent nodes creating ordered lists
//                 of neighboring global id's and owning proc nums
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::returnAdjNodesWithGround( const int & globalID,
  std::list<int> & gidList,
  std::list<int> & svGIDList,
  std::list<int> & procList,
  std::list<NodeID> & idList )
{
  if( cktgph.checkKey(globalID) )
  {
    NodeID id(cktgph.getKey1(globalID));

    std::vector<NodeID> adjIDs = cktgph.getAdjacent(id);

    int numAdj = adjIDs.size();
    for( int i = 0; i < numAdj; ++i )
    {
      CktNode * cktnode = cktgph.getData(adjIDs[i]);
      gidList.push_back( cktnode->get_gID() );
      idList.push_back( NodeID(cktnode->get_id(),cktnode->type()) );

      svGIDList.insert( svGIDList.end(),
        cktnode->get_SolnVarGIDList().begin(),
        cktnode->get_SolnVarGIDList().end() );

      procList.push_back( cktnode->get_ProcNum() );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerGIDswithDevs
// Purpose       : Loop over nodes and register int and ext global ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerGIDswithDevs()
{
  std::list<int> gidList, svGIDList, procList, ownedList;
  std::list<int>::const_iterator it_iL, end_iL, it_iL2, it_iL3;
  std::list<index_pair> tmpIPList1, tmpIPList2;
  std::list<NodeID> idList;

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    //------- Clear lists for each node
    gidList.clear();
    procList.clear();
    svGIDList.clear();
    ownedList.clear();
    idList.clear();

    //------  Generate global id list for external variables

    // initialize the local gidList, the svGIDList and the procList.
    returnAdjNodesWithGround( (*it_nL)->get_gID(), gidList, svGIDList, procList, idList );

    for( it_iL=gidList.begin(),end_iL=gidList.end(); it_iL!=end_iL; ++it_iL )
    {
      if( *it_iL != -1 ) ownedList.push_back(1);
      else               ownedList.push_back(0);
    }

    //------- Push list of int. and ext. global id's to ckt node

    // first arg. is internal variable list
    // 2nd   arg. is external variable list

    tmpIPList1.clear();
    tmpIPList2.clear();

    it_iL = (*it_nL)->get_SolnVarGIDList().begin();
    end_iL = (*it_nL)->get_SolnVarGIDList().end();
    if( (*it_nL)->get_IsOwned() )
      for( ; it_iL != end_iL; ++it_iL )
        tmpIPList1.push_back( index_pair( *it_iL, 1 ) );
    else
      for( ; it_iL != end_iL; ++it_iL )
        tmpIPList1.push_back( index_pair( *it_iL, 0 ) );

    int offset = (*it_nL)->get_Offset();

    it_iL = svGIDList.begin();
    end_iL = svGIDList.end();
    it_iL2 = ownedList.begin();

    std::list<int> extGIDs;

    for( ; it_iL != end_iL; ++it_iL, ++it_iL2 )
      if( offset )
        tmpIPList2.push_back( index_pair( *it_iL, *it_iL2, offset ) );
      else
        tmpIPList2.push_back( index_pair( *it_iL, *it_iL2 ) );

    (*it_nL)->registerGIDswithDev( tmpIPList1, tmpIPList2 );

    (*it_nL)->set_ExtSolnVarGIDList( svGIDList );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerStateGIDswithDevs
// Purpose       : Loop over nodes and register int and ext global ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerStateGIDswithDevs()
{
  std::list<int>::const_iterator it_iL, end_iL;
  std::list<index_pair> tmpIPList1;

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    tmpIPList1.clear();

    it_iL = (*it_nL)->get_StateVarGIDList().begin();
    end_iL = (*it_nL)->get_StateVarGIDList().end();
    if( (*it_nL)->get_IsOwned() )
    {
      for( ; it_iL != end_iL; ++it_iL )
        tmpIPList1.push_back( index_pair( *it_iL, 1 ) );
    }
    else
    {
      for( ; it_iL != end_iL; ++it_iL )
        tmpIPList1.push_back( index_pair( *it_iL, 0 ) );
    }

    //------- Push list of state global id's to ckt node
    (*it_nL)->registerStateGIDswithDev( tmpIPList1 );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerStoreGIDswithDevs
// Purpose       : Loop over nodes and register int and ext global ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktGraphBasic::registerStoreGIDswithDevs()
{
  std::list<int>::const_iterator it_iL, end_iL;
  std::list<index_pair> tmpIPList1;

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    tmpIPList1.clear();

    it_iL = (*it_nL)->get_StoreVarGIDList().begin();
    end_iL = (*it_nL)->get_StoreVarGIDList().end();
    if( (*it_nL)->get_IsOwned() )
    {
      for( ; it_iL != end_iL; ++it_iL )
        tmpIPList1.push_back( index_pair( *it_iL, 1 ) );
    }
    else
    {
      for( ; it_iL != end_iL; ++it_iL )
        tmpIPList1.push_back( index_pair( *it_iL, 0 ) );
    }

    //------- Push list of store global id's to ckt node
    (*it_nL)->registerStoreGIDswithDev( tmpIPList1 );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerLIDswithDevs
// Purpose       : Loop over nodes and register int and ext local ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerLIDswithDevs( Indexor & indexor )
{
  std::list<int> gidList, svGIDList, procList, ownedList;
  std::list<int>::const_iterator it_iL, end_iL, it_iL2, it_iL3;
  std::vector<int> intVec, extVec;
  std::list<NodeID> idList;

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    //------- Clear lists for each node
    gidList.clear();
    procList.clear();
    svGIDList.clear();
    ownedList.clear();
    idList.clear();

    //------  Generate global id list for external variables

    // initialize the local gidList, the svGIDList and the procList.
    returnAdjNodesWithGround( (*it_nL)->get_gID(), gidList, svGIDList, procList, idList );

    for( it_iL=gidList.begin(),end_iL=gidList.end(); it_iL!=end_iL; ++it_iL )
    {
      if( *it_iL != -1 ) ownedList.push_back(1);
      else               ownedList.push_back(0);
    }

    //------- Push list of int. and ext. local id's to ckt node

    // first arg. is internal variable list
    // 2nd   arg. is external variable list

    it_iL = (*it_nL)->get_SolnVarGIDList().begin();
    end_iL = (*it_nL)->get_SolnVarGIDList().end();

    std::string mapName("SOLUTION_OVERLAP_GND");

    intVec.assign( it_iL, end_iL );

    bool success = indexor.globalToLocal( mapName, intVec );

    it_iL = svGIDList.begin();
    end_iL = svGIDList.end();

    extVec.assign( it_iL, end_iL );

    success = success && indexor.globalToLocal( mapName, extVec );

    (*it_nL)->registerLIDswithDev( intVec, extVec );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerStateLIDswithDevs
// Purpose       : Loop over nodes and register int and ext global ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerStateLIDswithDevs( Indexor & indexor )
{
  std::list<int>::const_iterator it_iL, end_iL;
  std::vector<int> stateVec;

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    it_iL = (*it_nL)->get_StateVarGIDList().begin();
    end_iL = (*it_nL)->get_StateVarGIDList().end();

    std::string mapName("STATE_OVERLAP");

    stateVec.assign( it_iL, end_iL );

    indexor.globalToLocal( mapName, stateVec );

    //------- Push list of state global id's to ckt node
    (*it_nL)->registerStateLIDswithDev( stateVec );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerStoreLIDswithDevs
// Purpose       : Loop over nodes and register int and ext global ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktGraphBasic::registerStoreLIDswithDevs( Indexor & indexor )
{
  std::list<int>::const_iterator it_iL, end_iL;
  std::vector<int> storeVec;

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    it_iL = (*it_nL)->get_StoreVarGIDList().begin();
    end_iL = (*it_nL)->get_StoreVarGIDList().end();

    std::string mapName("STORE_OVERLAP");

    storeVec.assign( it_iL, end_iL );

    indexor.globalToLocal( mapName, storeVec );

    //------- Push list of store global id's to ckt node
    (*it_nL)->registerStoreLIDswithDev( storeVec );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerDepLIDswithDevs
// Purpose       : Loop over nodes and register dep local ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerDepLIDswithDevs( Indexor & indexor )
{
  std::vector< std::vector<int> > indexVec;
  std::string mapName("SOLUTION_OVERLAP_GND");

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    indexVec = (*it_nL)->get_DepSolnGIDVec();

    for( unsigned int i = 0; i < indexVec.size(); ++i )
      indexor.globalToLocal( mapName, indexVec[i] );

    //------- Push list of local id's to ckt node

    (*it_nL)->registerDepLIDswithDev( indexVec );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerDepStateLIDswithDevs
// Purpose       : Loop over nodes and register dep state local ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerDepStateLIDswithDevs( Indexor & indexor )
{
  std::vector< std::vector<int> > indexVec;
  std::string mapName("STATE_OVERLAP");

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    indexVec = (*it_nL)->get_DepStateGIDVec();

    for( unsigned int i = 0; i < indexVec.size(); ++i )
      indexor.globalToLocal( mapName, indexVec[i] );

    //------- Push list of state local id's to ckt node
    (*it_nL)->registerDepStateLIDswithDev( indexVec );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerDepStoreLIDswithDevs
// Purpose       : Loop over nodes and register dep store local ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktGraphBasic::registerDepStoreLIDswithDevs( Indexor & indexor )
{
  std::vector< std::vector<int> > indexVec;
  std::string mapName("STORE_OVERLAP");

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    indexVec = (*it_nL)->get_DepStoreGIDVec();

    for( unsigned int i = 0; i < indexVec.size(); ++i )
      indexor.globalToLocal( mapName, indexVec[i] );

    //------- Push list of store local id's to ckt node
    (*it_nL)->registerDepStoreLIDswithDev( indexVec );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::registerJacLIDswithDevs
// Purpose       : Loop over nodes and register jacobian local ids
//                 with each device
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::registerJacLIDswithDevs( Indexor & indexor )
{
  std::string graphName("JACOBIAN_OVERLAP_GND");
  std::string mapName("SOLUTION_OVERLAP_GND");
  std::vector< std::vector<int> > stampVec;

  indexor.setupAcceleratedMatrixIndexing( graphName );

  std::list<CktNode*> * tmpNodeList = getBFSNodeList();

  std::list<CktNode*>::const_iterator it_nL = tmpNodeList->begin();
  std::list<CktNode*>::const_iterator it_nL_end = tmpNodeList->end();

  // loop over nodes:
  for( ; it_nL != it_nL_end; ++it_nL )
  {
    CktNode & cn = **it_nL;
    if( cn.type() == _DNODE )
    {
      const std::list<int> & intGIDs = cn.get_SolnVarGIDList();
      const std::list<int> & extGIDs = cn.get_ExtSolnVarGIDList();
      const std::vector<int> & depGIDs = cn.get_DepSolnGIDJacVec();
      std::vector<int> gids( intGIDs.size() + extGIDs.size() + depGIDs.size() );
      copy( extGIDs.begin(), extGIDs.end(), gids.begin() );
      copy( intGIDs.begin(), intGIDs.end(), gids.begin() + extGIDs.size() );
      copy( depGIDs.begin(), depGIDs.end(), gids.begin() + extGIDs.size() + intGIDs.size() );

      stampVec = cn.jacobianStamp();

      int numRows = stampVec.size();
      for( int i = 0; i < numRows; ++i )
      {
        int numCols = stampVec[i].size();
        for( int j = 0; j < numCols; ++j )
          stampVec[i][j] = gids[ stampVec[i][j] ];
      }

      std::vector<int> counts(3);
      counts[0] = extGIDs.size();
      counts[1] = intGIDs.size();
      counts[2] = depGIDs.size();
      cn.registerGIDDataWithDev( counts, gids, stampVec );

      indexor.matrixGlobalToLocal( graphName, gids, stampVec );

      cn.registerJacLIDswithDev( stampVec );
    }
  }

  indexor.deleteAcceleratedMatrixIndexing();
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::regenerateGIDNodeMap
// Purpose       : Redo global index node map
// Special Notes : Should find a way to automate this
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
void CktGraphBasic::regenerateGIDNodeMap()
{
  std::list<CktNode*>::const_iterator it_cnL = BFSNodeList_.begin();
  std::list<CktNode*>::const_iterator it_cnL_end = BFSNodeList_.end();
  for( ; it_cnL != it_cnL_end ; ++it_cnL )
  {
    cktgph.chgKey2( NodeID((*it_cnL)->get_id(),(*it_cnL)->type()), (*it_cnL)->get_gID() );
  }
}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::replaceNode
// Purpose       : replace supernode two nodes
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/10/10
//-----------------------------------------------------------------------------
CktNode * CktGraphBasic::replaceNode( const NodeID nodeToBeReplaced, 
                                                  const NodeID nodeToKeep )
{
  // procedure is as follows:
  // 1. get the adjacency of both nodes
  // 2. combine to the two adjacencies
  // 3. set nodeToKeep's adjacency to what we found in step 2.
  // 4. tell the graph to replace any other adjacencies that refer to nodeToBeReplaced with nodeToKeep
  // 5. Remove nodeToBeReplaced and return it's associated CktNode object
  
  // look up key for these nodes in graph.
  CktNode * nodeToBeReplacedCktNodePtr = FindCktNode( nodeToBeReplaced );
  
  // 1. get the adjacency of both nodes
  std::vector<NodeID> adjNodeToBeReplaced, adjNodeToKeep;
  returnAdjIDs( nodeToBeReplaced, adjNodeToBeReplaced );
  returnAdjIDs( nodeToKeep, adjNodeToKeep );
  
  /*
  // 2. combine to the two adjacencies
  // loop over ajdNodeToBeReplaced and add any nodes not found in ajdNodeToKeep to ajdNodeToKeep
  std::vector<std::string>::iterator currNodeToBeReplacedAdjItr = adjNodeToBeReplaced.begin();
  std::vector<std::string>::iterator endNodeToBeReplacedAdjItr = adjNodeToBeReplaced.end();
  while( currNodeToBeReplacedAdjItr != endNodeToBeReplacedAdjItr )
  {
    std::vector<std::string>::iterator beginNodeToKeepAdjItr = ajdNodeToKeep.begin();
    std::vector<std::string>::iterator endNodeToKeepAdjItr = ajdNodeToKeep.end();
    std::vector<std::string>::iterator locationNodeToKeepAdjItr = find( beginNodeToKeepAdjItr, endNodeToKeepAdjItr, *currNodeToBeReplacedAdjItr);
    if( locationNodeToKeepAdjItr == endNodeToKeepAdjItr )
    {
      // this adjacency was not in the new list.  So add it in
      ajdNodeToKeep.push_back( *currNodeToBeReplacedAdjItr );
    }
    currNodeToBeReplacedAdjItr++;
  }
  
  
  // 3. set nodeToKeep's adjacency to what we found in step 2.
  cktgph.setAdjacent( nodeToBeReplaced, nodeToKeep, ajdNodeToKeep );
  */
  cktgph.addToAdjacent( nodeToBeReplaced, nodeToKeep, adjNodeToBeReplaced );
  
  // 4. tell the graph to replace any other adjacencies that refer to nodeToBeReplaced with nodeToKeep
  cktgph.replaceAdjacent( nodeToBeReplaced, nodeToKeep );
  
  // 5. Remove nodeToBeReplaced and return it's associated CktNode object
  cktgph.removeKey( nodeToBeReplaced );

  // changed ordering so set isModified_ flag
  isModified_ = true;
  
  return nodeToBeReplacedCktNodePtr;
}


//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::removeRedundantDevices
// Purpose       : remove any device nodes that are only connected to one ckt node
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/10/10
//-----------------------------------------------------------------------------
void CktGraphBasic::removeRedundantDevices(std::vector< CktNode * > & removedDevices)
{
  std::ostringstream outputStringStream;
  // Collect the deviceIDs for all the devices to be removed from the graph
  std::vector< NodeID> deviceIDs;
  const std::map< NodeID, CktNode* > dataMap = cktgph.getData1Map();

  std::map< NodeID, CktNode* >::const_iterator currentCktNodeItr = dataMap.begin();
  std::map< NodeID, CktNode* >::const_iterator endCktNodeItr = dataMap.end();
  while( currentCktNodeItr != endCktNodeItr )
  {
    if( ((*currentCktNodeItr).second)->type() == _DNODE )
    {
      CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>((*currentCktNodeItr).second);
      Teuchos::RCP<Device::InstanceBlock> deviceInstanceBlockPtr = cktNodeDevPtr->devBlock();
      if( Teuchos::nonnull( deviceInstanceBlockPtr ) )
      {
        // have a valid device.
        NodeID deviceID = (*currentCktNodeItr).first;
        int deviceIndex = cktgph.getIndex( deviceID );
        std::vector< int > adjacentIDs = cktgph.getAdjacentRow( deviceIndex );
        int numAdjIDs = adjacentIDs.size(); 
 
        bool removeDevice=false;
        if( numAdjIDs == 2 )
        {
          if( adjacentIDs[0] == adjacentIDs[1] )
            removeDevice=true;
        }
        else if( numAdjIDs == 3 )
        {
          if( ( adjacentIDs[0] == adjacentIDs[1] ) && (adjacentIDs[1] == adjacentIDs[2]) )
            removeDevice=true;
        }
        else if( numAdjIDs == 4 )
        {
          if( ( adjacentIDs[0] == adjacentIDs[1] ) && (adjacentIDs[1] == adjacentIDs[2])
            && (adjacentIDs[2] == adjacentIDs[3]) )
            removeDevice=true;
        }
        if( removeDevice )
        {
          deviceIDs.push_back(deviceID);
          removedDevices.push_back((*currentCktNodeItr).second);
        }
      }
      else
      {
        // issue fatal error as this case shouldn't occur
        std::string msg("Topology::removeRedundantDevices() null Device Instance Block pointer.");
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }
    currentCktNodeItr++;
  }

  // If there are any devices to remove, remove them from the graph and mark the object as modified.
  if( removedDevices.size() > 0 )
  {
    cktgph.removeKeys( deviceIDs );

    // After the devices are removed, check the graph to see if any singletons were created.
    // These are most likely ghost nodes that connected the removed devices to the rest of the circuit.
    // NOTE:  This can be a real issue in practice, so don't remove this search.
    std::vector< NodeID> singletonIDs = cktgph.getSingletons();
    if (singletonIDs.size() > 0)
      cktgph.removeKeys( singletonIDs );

    isModified_=true;
  }

}

//-----------------------------------------------------------------------------
// Function      : CktGraphBasic::put
// Purpose       : Allows virtual override of operator<<
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
std::ostream& CktGraphBasic::put(std::ostream& os) const
{
  int i = 0;

  std::list<CktNode*>::const_iterator it_cnL = BFSNodeList_.begin();
  std::list<CktNode*>::const_iterator it_cnL_end = BFSNodeList_.end();
  for( ; it_cnL != it_cnL_end ; ++it_cnL, ++i )
  {
    os << "[" << i << "]:" << **it_cnL << std::endl;
  }

  // Print out information about the circuit graph.
  cktgph.print( os );

  return os;
}

} // namespace Topo
} // namespace Xyce
