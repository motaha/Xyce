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
// Filename       : $RCSfile: N_TOP_TopoLSUtil.C,v $
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
// Revision Number: $Revision: 1.136.2.2 $
//
// Revision Date  : $Date: 2014/03/06 04:01:46 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <numeric>
#include <assert.h>

#include <N_TOP_TopoLSUtil.h>
#include <N_TOP_Topology.h>
#include <N_TOP_Misc.h>
#include <N_TOP_NodeBlock.h>
#include <N_TOP_CktGraphSupport.h>
#include <N_TOP_CktGraphCreator.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_CktNodeCreator.h>
#include <N_TOP_CktNode.h>

#include <N_DEV_DeviceInterface.h>

#include <N_PDS_Manager.h>
#include <N_PDS_GlobalAccessor.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>

#include <N_PDS_Directory.h>
#include <N_PDS_Migrate.h>
#include <N_PDS_Node.h>

#include <N_UTL_OptionBlock.h>

#include <N_IO_CmdParse.h>
#include <N_IO_PkgOptionsMgr.h>

#include <N_ERH_ErrorMgr.h>

#include <Epetra_Util.h>
#include <Teuchos_Utils.hpp>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::TopoLSUtil
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter  SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
TopoLSUtil::TopoLSUtil(Topology * topo, IO::CmdParse & cp)
  : topoPtr_(topo),
    commandLine_ (cp),
    pdsMgrPtr_(0),
    nodeGlobalAccessorPtr_(0),
    solnGlobalAccessorPtr_(0),
    stateGlobalAccessorPtr_(0),
    storeGlobalAccessorPtr_(0),
    numGlobalNodes_(0),
    numLocalNodes_(0),
    baseNodeGID_(0),
    numGlobalRows_(0),
    numLocalRows_(0),
    numExternRows_(0),
    numGlobalExternRows_(0),
    baseRowGID_(0),
    numGlobalStateVars_(0),
    numLocalStateVars_(0),
    numExternStateVars_(0),
    numGlobalExternStateVars_(0),
    baseStateVarGID_(0),
    numGlobalStoreVars_(0),
    numLocalStoreVars_(0),
    numExternStoreVars_(0),
    numGlobalExternStoreVars_(0),
    baseStoreVarGID_(0),
    numGlobalNZs_(0),
    numLocalNZs_(0),
    checkConnectivity_(true),
    supernode_(false),
#ifdef Xyce_TEST_SOLN_VAR_MAP
    namesFile_(true)
#else
    namesFile_(false)
#endif
{

}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::operator<<
// Purpose       : generate utility with reference to topology
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/15/00
//-----------------------------------------------------------------------------
std::ostream & operator<< (std::ostream & os, const TopoLSUtil & tlsu )
{
  os << "Topology LS Utility" << std::endl;
  os << "-------------------" << std::endl;
  os << "ProcID: " << tlsu.pdsMgrPtr_->getPDSComm()->procID() << std::endl;
  os << "Num Global Nodes: " << tlsu.numGlobalNodes_ << std::endl;
  os << "Num Local Nodes: " << tlsu.numLocalNodes_ << std::endl;
  os << "Num Global Rows: " << tlsu.numGlobalRows_ << std::endl;
  os << "Num Local Rows: " << tlsu.numLocalRows_ << std::endl;
  os << "Num Extern Rows: " << tlsu.numExternRows_ << std::endl;
  os << "Num Global Extern Rows: " << tlsu.numGlobalExternRows_ << std::endl;
  os << "Row GID Base: " << tlsu.baseRowGID_ << std::endl;

  os << "Num Global State Vars: " << tlsu.numGlobalStateVars_ << std::endl;
  os << "Num Local State Vars: " << tlsu.numLocalStateVars_ << std::endl;
  os << "Num Extern State Vars: " << tlsu.numExternStateVars_ << std::endl;
  os << "Num Global Extern State Vars: " << tlsu.numGlobalExternStateVars_ << std::endl;
  os << "State Var GID Base: " << tlsu.baseStateVarGID_ << std::endl;

  os << "Num Global Store Vars: " << tlsu.numGlobalStoreVars_ << std::endl;
  os << "Num Local Store Vars: " << tlsu.numLocalStoreVars_ << std::endl;
  os << "Num Extern Store Vars: " << tlsu.numExternStoreVars_ << std::endl;
  os << "Num Global Extern Store Vars: " << tlsu.numGlobalExternStoreVars_ << std::endl;
  os << "Store Var GID Base: " << tlsu.baseStoreVarGID_ << std::endl;

  os << "Num Global NZs: " << tlsu.numGlobalNZs_ << std::endl;
  os << "Num Local NZs: " << tlsu.numLocalNZs_ << std::endl;

  os << "Node Array: ";
  for( int i = 0; i < tlsu.numLocalNodes_; ++i )
    os << tlsu.nodeList_GID_[i] << " ";
  os << std::endl;

  os << "Extern Node Array: ";
  for( unsigned int i = 0; i < tlsu.nodeList_ExternGID_.size(); ++i )
    os << tlsu.nodeList_ExternGID_[i].first << " " <<
	tlsu.nodeList_ExternGID_[i].second << "   ";
  os << std::endl;

  os << "GID Array: ";
  for( int i = 0; i < tlsu.numLocalRows_; ++i )
    os << tlsu.rowList_GID_[i] << " ";
  os << std::endl;

  os << "Extern GID Array: ";
  for( unsigned int i = 0; i < tlsu.rowList_ExternGID_.size(); ++i )
    os << tlsu.rowList_ExternGID_[i].first << " " <<
	tlsu.rowList_ExternGID_[i].second << "   ";
  os << std::endl;

  os << "State GID Array: ";
  for( int i = 0; i < tlsu.numLocalStateVars_; ++i )
    os << tlsu.rowList_StateGID_[i] << " ";
  os << std::endl;

  os << "Extern State GID Array: ";
  for( unsigned int i = 0; i < tlsu.rowList_ExternStateGID_.size(); ++i )
    os << tlsu.rowList_ExternStateGID_[i].first << " " <<
	tlsu.rowList_ExternStateGID_[i].second << "   ";
  os << std::endl;

  os << "Store GID Array: ";
  for( int i = 0; i < tlsu.numLocalStoreVars_; ++i )
    os << tlsu.rowList_StoreGID_[i] << " ";
  os << std::endl;

  os << "Extern Store GID Array: ";
  for( unsigned int i = 0; i < tlsu.rowList_ExternStoreGID_.size(); ++i )
    os << tlsu.rowList_ExternStoreGID_[i].first << " " <<
	tlsu.rowList_ExternStoreGID_[i].second << "   ";
  os << std::endl;

  os << "NZ Array: ";
  for( int i = 0; i < tlsu.numLocalRows_; ++i )
    os << tlsu.rowList_NumNZs_[i] << " ";
  os << std::endl;

  os << "Col Index Array: " << std::endl;
  for( int i = 0; i < tlsu.numLocalRows_; ++i )
  {
    os << tlsu.rowList_GID_[i] << ": ";
    for( std::list<int>::const_iterator it_iL = tlsu.rowList_ColList_[i].begin();
	it_iL != tlsu.rowList_ColList_[i].end(); ++it_iL )
      os << (*it_iL) << " ";
    os << std::endl;
  }

  return os;
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::setupRowCol
// Purpose       : Setup row/col data for linear solver including reorder
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/1/01
//-----------------------------------------------------------------------------
bool TopoLSUtil::setupRowCol()
{
  topoPtr_->setOrderedNodeList();

  N_PDS_Comm & comm = *(pdsMgrPtr_->getPDSComm());
  int procCnt = comm.numProc();
  int procID = comm.procID();

#ifdef Xyce_DEBUG_PIO_TOPO
  Xyce::dout() << "DEBUGGING NEW TOPOLOGY SETUP FOR PARALLEL IO: " << procID << "\n";
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "<<<<<<<<<<<<<<<<<NODES>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
  std::list<CktNode*>::iterator itL = topoPtr_->orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator std::endl = topoPtr_->orderedNodeListPtr_->end();
  for( int i = 0 ; itL != std::endl; ++itL, ++i )
  {
    Xyce::dout() << "Proc: " << procID << "\t" << "Node: " << i << std::endl;
    Xyce::dout() << **itL;
  }
  Xyce::dout() << "<<<<<<<<<<<<<<NODES END>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
#endif

  //construct v-node and d-node directories
  typedef RCP< Xyce::Parallel::IndexNode > IndexNodePtr;

  typedef std::multimap< std::string, IndexNodePtr > VNodeContainer;
  typedef std::map< std::string, IndexNodePtr > DNodeContainer;

  typedef std::map< std::string, IndexNodePtr > INodeContainer;

  typedef Xyce::Parallel::Hash<std::string> StringHash;
  typedef Xyce::Parallel::Migrate<std::string,Xyce::Parallel::IndexNode> INMigrate;

  typedef Xyce::Parallel::Directory< std::string,
                             Xyce::Parallel::IndexNode,
                             StringHash,
                             VNodeContainer,
                             INMigrate >
          VNodeDir;
  typedef Xyce::Parallel::Directory< std::string,
                             Xyce::Parallel::IndexNode,
                             StringHash,
                             DNodeContainer,
                             INMigrate >
          DNodeDir;

  StringHash SHobj( procCnt );
  INMigrate Mobj( comm );

  VNodeDir VDir( Mobj, SHobj );
  DNodeDir DDir( Mobj, SHobj );

  //loop over node list and setup IndexNodes for Vs and Ds
  //register with VDir and DDir

  //data for directory registration
  INodeContainer VData;
  INodeContainer DData;

  //vectors for retrieval
  std::vector<std::string> VNames;
  std::vector<std::string> DNames;

  //find all voltage nodes connected to a voltage source in case we
  // need to force these nodes to be on the same processor
  std::set<std::string> Vsrc_Connected_Nodes_;
  std::list<CktNode*>::iterator it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator end_cnL = topoPtr_->orderedNodeListPtr_->end();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    std::vector<NodeID> adj_ids;

    const std::string & id = (*it_cnL)->get_id();
    std::string::size_type col = id.find_first_of(':');

#ifdef Xyce_DEBUG_TOPOLOGY
      Xyce::dout() << "Considering nodes for: " << id << std::endl;
#endif
    int type = (*it_cnL)->type();
    if  (type == _DNODE)
    {
      if ( id[col+1] == 'V' || id[col+1] == 'v'
                          || id.substr(col+1,col+7) == "y%iso2"
                          || id.substr(col+1,col+7) == "Y%ISO2"
                          || id.substr(col+1,col+6) == "y%ext"
                          || id.substr(col+1,col+6) == "Y%EXT" )
      {
#ifdef Xyce_DEBUG_TOPOLOGY
        Xyce::dout() << "Getting adjacent nodes for: " << id << std::endl;
#endif
        topoPtr_->mainGraphPtr_->returnAdjIDs( NodeID(id,type), adj_ids );
        int adjSize = adj_ids.size();
        for( int i = 0; i < adjSize; ++i )
        {
#ifdef Xyce_DEBUG_TOPOLOGY
          Xyce::dout() << "adj_ids["<<i<<"] = " << adj_ids[i] << std::endl;
#endif
          if( adj_ids[i].first != "0" )
          {
            Vsrc_Connected_Nodes_.insert(adj_ids[i].first);
          }
        }
      }
    }
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Vsrc_Connected_Nodes:"<<std::endl;
  std::set<std::string>::iterator its = Vsrc_Connected_Nodes_.begin();
  std::set<std::string>::iterator fts = Vsrc_Connected_Nodes_.end();
  for( ; its != fts; ++its ) Xyce::dout() << *its << std::endl;
#endif

  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    const std::string & id = (*it_cnL)->get_id();
    int type = (*it_cnL)->type();

    if( id != "0" )
    {
      IndexNodePtr inode( new Xyce::Parallel::IndexNode( -99, procID ) );

      if( type == _VNODE )
      {
        if( Vsrc_Connected_Nodes_.count(id) ) inode->gid = -98;

#ifdef Xyce_DEBUG_TOPOLOGY
        if( inode->gid == -98 ) Xyce::dout() << "Node: " << id << std::endl;
#endif

        VData[id] = inode;
        VNames.push_back( id );
      }
      else if( type == _DNODE )
      {
        DData[id] = inode;
        DNames.push_back( id );
      }
      else
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
          "TopoLSUtil::setupRowCol - Unrecoqnized Node Type: " + id + "\n" );
    }
  }

  VDir.addEntries( VData );
  DDir.addEntries( DData );

  //v-node use multimap, pick owner
  //locally loop over container in VDir, pick an owner
  //random choice of owner
  VNodeContainer & VNodes = VDir.container();
  VNodeContainer::iterator iterVN = VNodes.begin();
  VNodeContainer::iterator endVN = VNodes.end();
  VNodeContainer OwnedVNodes;
  std::string id("");
  while( iterVN != endVN )
  {
      int proc = -1;

      id = iterVN->first;
      if( iterVN->second->gid == -98 ) proc = iterVN->second->pid;

      IndexNodePtr inode( new Xyce::Parallel::IndexNode( -99, iterVN->second->pid ) );

      ++iterVN;
      if( (iterVN != endVN) && (id == iterVN->first) )
      {
        std::vector<int> intVec;
        intVec.push_back(inode->pid);
        while( (iterVN != endVN) && (iterVN->first == id) )
        {
          if( iterVN->second->gid == -98 ) proc = iterVN->second->pid;

          intVec.push_back( iterVN->second->pid );
          ++iterVN;
        }
        random_shuffle( intVec.begin(), intVec.end() );
#ifdef Xyce_EXTDEV
        if( proc != -1 )
        {
          inode->pid = proc;
#ifdef Xyce_DEBUG_TOPOLOGY
          Xyce::dout() << "pNode: " << id << " " << proc << std::endl;
#endif
        }
        else
#endif
          inode->pid = *(intVec.begin());
      }

      OwnedVNodes.insert( VNodeContainer::value_type( id, inode ) );
  }
  if( VNodes.size() > OwnedVNodes.size() )
  {
    VNodes.clear();
    VNodes = OwnedVNodes;
  }
  OwnedVNodes.clear();

  //gids for both
  //locally loop over containers in VDir and DDir and set GIDs
  int VSize = VNodes.size();

  double tmpVar1, tmpVar2;

  tmpVar1 = static_cast<double>(VSize);
  comm.scanSum( &tmpVar1, &tmpVar2, 1 );
  int baseVNodeGID = static_cast<int>(tmpVar2) - VSize;

  comm.sumAll( &tmpVar1, &tmpVar2, 1 );
  int globalVNodeCnt = static_cast<int>(tmpVar2);

  int currGID = baseVNodeGID;
  iterVN = VNodes.begin();
  for( ; iterVN != endVN; ++iterVN, ++currGID )
    iterVN->second->gid = currGID;

  //same for DNodes
  DNodeContainer & DNodes = DDir.container();
  DNodeContainer::iterator iterDN = DNodes.begin();
  DNodeContainer::iterator endDN = DNodes.end();

  int DSize = DNodes.size();

  tmpVar1 = static_cast<double>(DSize);
  comm.scanSum( &tmpVar1, &tmpVar2, 1 );
  int baseDNodeGID = static_cast<int>(tmpVar2) - DSize + globalVNodeCnt;

  currGID = baseDNodeGID;
  for( ; iterDN != endDN; ++iterDN, ++currGID )
    iterDN->second->gid = currGID;

  //push indices back
  //do gets on VDir and DDir to get ownership and GIDs
  VData.clear();
  DData.clear();

  VDir.getEntries( VNames, VData );
  DDir.getEntries( DNames, DData );

  IndexNodePtr inode;
  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    CktNode & cn = **it_cnL;
    const std::string & id = cn.get_id();
    int type = cn.type();

    if( id != "0" )
    {
      if( type == _VNODE )
        inode = VData[id];
      else if( type == _DNODE )
        inode = DData[id];

      cn.set_gID( inode->gid );
      cn.set_ProcNum( inode->pid );
      cn.set_IsOwned( procID == inode->pid );
    }
    else
    {
      cn.set_gID( -1 );
      cn.set_ProcNum( -1 );
      cn.set_IsOwned( false );
    }
  }

#ifdef Xyce_DEBUG_PIO_TOPO
  comm.barrier();

  for( int i = 0; i < procCnt; ++i )
  {
    if( i == procID )
    {
      Xyce::dout() << "<<<<<<<<<<<<<REINDEXED NODES>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
      std::list<CktNode*>::iterator itL2 = topoPtr_->orderedNodeListPtr_->begin();
      std::list<CktNode*>::iterator std::endl2 = topoPtr_->orderedNodeListPtr_->end();
      for( int i = 0 ; itL2 != std::endl2; ++itL2, ++i )
      {
        Xyce::dout() << "Proc: " << procID << "\t" << "Node: " << i << std::endl;
        Xyce::dout() << **itL2;
      }
      Xyce::dout() << "<<<<<<<<REINDEXED NODES END>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
    }
    comm.barrier();
  }

  comm.barrier();
#endif

  //Setup GID maps in CktGraph
  topoPtr_->regenerateGIDNodeMap();

  //Check Node GID setup
  //Check Soln/State GID setup

  //Register GIDs should be fine
  //Check Generation of RowCol Data

  //Directory Construction Replace?
  //Reset of RowCol Data?

#ifdef Xyce_PARALLEL_IO_GHOST
  typedef RCP<NodeDevBlock> NodeDevBlockPtr;
  typedef std::vector<NodeDevBlockPtr> GhostContainer;

  std::vector<int> pids;
  GhostContainer Ghosts;

  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  end_cnL = topoPtr_->orderedNodeListPtr_->end();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    const std::string & id = (*it_cnL)->get_id();
    int type = (*it_cnL)->type();
    bool owned = (*it_cnL)->get_IsOwned();
    int proc = (*it_cnL)->get_ProcNum();
    int gid = (*it_cnL)->get_gID();

    if( type == _VNODE && !owned  )
    {
      //find all connected devices to be ghosted
      std::vector<int> adj_gids;
      topoPtr_->mainGraphPtr_->returnAdjGIDs( gid, adj_gids );

      for( int i = 0; i < adj_gids.size(); ++i )
      {
        NodeDevBlockPtr ndb( topoPtr_->mainGraphPtr_->returnNodeDevBlock( adj_gids[i] ) );
        ndb->node.set_IsOwned( false );
        Ghosts.push_back( ndb );
        pids.push_back( proc );
      }
    }
  }

  typedef Xyce::Parallel::Migrate<NodeDevBlock> NDBMigrateType;
  NDBMigrateType NDBMigrate( comm );

  GhostContainer RecdGhosts;
  NDBMigrate( pids, Ghosts, RecdGhosts );

  //loop over rec'd ghosts and add to local ckt
#endif

#ifndef Xyce_STATIC_SET_GID

  setupNodeGIDs();

#ifdef Xyce_DEBUG_PIO_TOPO
  comm.barrier();

  for( int i = 0; i < procCnt; ++i )
  {
    if( i == procID )
    {
      Xyce::dout() << "<<<<<<<<<<<<<REINDEXED NODES>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
      std::list<CktNode*>::iterator itL2 = topoPtr_->orderedNodeListPtr_->begin();
      std::list<CktNode*>::iterator std::endl2 = topoPtr_->orderedNodeListPtr_->end();
      for( int i = 0 ; itL2 != std::endl2; ++itL2, ++i )
      {
        Xyce::dout() << "Proc: " << procID << "\t" << "Node: " << i << std::endl;
        Xyce::dout() << **itL2;
      }
      Xyce::dout() << "<<<<<<<<REINDEXED NODES END>>>>>>>>>>>>>>>>>>>>>> " << procID << "\n";
    }
    comm.barrier();
  }

  comm.barrier();
#endif

  setupSolnAndStateGIDs();

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << *topoPtr_ << std::endl;
#endif

  topoPtr_->registerGIDswithDevs();

  //  generateRowColData(); (HKT, 01/17/07: This was a duplicate call to generateRowColData())

#else

#ifdef Xyce_VERBOSE_SETUPTOPO
  Xyce::dout() << "Generating Row Col Data" << std::endl;
#endif
  //Get initial row/col data
  //------------------------
  generateRowColData();

#ifdef Xyce_VERBOSE_SETUPTOPO
  Xyce::dout() << "Generating initial Global Accessors" << std::endl;
#endif
  //Setup temp global accessors
  //---------------------------
  solnGlobalAccessorPtr_ = pdsMgrPtr_->createGlobalAccessor();
  stateGlobalAccessorPtr_ = pdsMgrPtr_->createGlobalAccessor();
  storeGlobalAccessorPtr_ = pdsMgrPtr_->createGlobalAccessor();

  solnGlobalAccessorPtr_->registerExternGIDVector( rowList_ExternGID_ );
  stateGlobalAccessorPtr_->registerExternGIDVector( rowList_ExternStateGID_ );
  storeGlobalAccessorPtr_->registerExternGIDVector( rowList_ExternStoreGID_ );

  solnGlobalAccessorPtr_->generateMigrationPlan();
  stateGlobalAccessorPtr_->generateMigrationPlan();
  storeGlobalAccessorPtr_->generateMigrationPlan();

#ifdef Xyce_VERBOSE_SETUPTOPO
  Xyce::dout() << "Reordering GIDs" << std::endl;
#endif
  //Reorder GIDs
  //------------
  reorderGIDs();

#ifdef Xyce_VERBOSE_SETUPTOPO
  Xyce::dout() << "Registering GIDs" << std::endl;
#endif
  //Register GIDs with Devices
  //--------------------------
  topoPtr_->registerGIDswithDevs();

#ifdef Xyce_VERBOSE_SETUPTOPO
  Xyce::dout() << "Generating Row Col Data" << std::endl;
#endif
  //Get final row/col data
  //------------------------
  generateRowColData();

#ifdef Xyce_VERBOSE_SETUPTOPO
  Xyce::dout() << "Deleting initial Global Accessors" << std::endl;
#endif
  //delete temporary global accessors
  //---------------------------------
  delete solnGlobalAccessorPtr_;
  delete stateGlobalAccessorPtr_;
  delete storeGlobalAccessorPtr_;

#endif

  if( checkConnectivity_ )
    testVoltageNodeConnectivity_();

  topoPtr_->generateDirectory();

  //resolve late dependencies for devices
  topoPtr_->resolveDependentVars();

  generateRowColData();

  topoPtr_->returnSVarVNodeGIDVec( vnodeGIDVector_ );

  topoPtr_->returnSVarVsrcGIDVec( vsrcGIDVector_ );

  topoPtr_->returnSVarNoDCPathIDVec( noDCPathIDVector_ );

  topoPtr_->returnSVarConnToOneTermIDVec( connToOneTermIDVector_ );


  //Added 1/3/08, KRS:  For correct parallel implementation, we need to
  //put all of the information from each proc in the noDCPathGIDVector_ and
  //connToOneTermGIDVector_ back onto proc 0.  If we don't, and the parallel
  //run is taking place on different machines (rather than one machine with
  //multiple procs), then each machine will open up a local copy of the netlist
  //copy file and write to this file *only* those resistors which correspond to
  //nodes owned by that machine.

#ifdef Xyce_PARALLEL_MPI
  //Need to compute and store the number of nodes in noDCPathGIDVector_ on
  //each processor

  std::vector<int> nodclength_local(procCnt,0);
  std::vector<int> nodclength_sum(procCnt,0);
  int j;

  for (j=0; j < procCnt ; j++)
  {
    if (j != procID)
      nodclength_local[j] = 0;
    else
      nodclength_local[j] = noDCPathIDVector_.size();
  }
  N_ERH_ErrorMgr::safeBarrier(pdsMgrPtr_->getPDSComm()->comm());
  //  for (j=0; j < procCnt; j++)
  //  comm.sumAll(&nodclength_local[j],&nodclength_sum[j],procCnt);

  comm.sumAll(&nodclength_local[0],&nodclength_sum[0],procCnt);

  //Broadcast the nodes to proc 0....
  if (procID != 0)
  {
    std::vector<std::string>::iterator nodcit = noDCPathIDVector_.begin();
    std::vector<std::string>::iterator nodcend = noDCPathIDVector_.end();
    int templength;

    while (nodcit != nodcend)
    {
      templength = (*nodcit).length();
      comm.send(&templength,1,0);
      comm.send(&(*nodcit)[0],templength,0);
      nodcit++;
    }
  }
  //....and let proc 0 receive.
  else
  {
    int templength;
    for (j=1; j < procCnt; j++)
    {
      for (int i=0; i < nodclength_sum[j]; i++)
      {
        comm.recv(&templength,1,j);
        std::string tempstring;
        tempstring.resize(templength);
        comm.recv(&tempstring[0],templength,j);
        noDCPathIDVector_.push_back(tempstring);
      }
    }
  }
  N_ERH_ErrorMgr::safeBarrier(pdsMgrPtr_->getPDSComm()->comm());

  //Need to compute and store the number of nodes in connToOneTermGIDVector_ on
  //each processor

  std::vector<int> connoneterm_local(procCnt,0);
  std::vector<int> connoneterm_sum(procCnt,0);

  for (j=0; j < procCnt ; j++)
    {
      if (j != procID)
        connoneterm_local[j] = 0;
      else
        connoneterm_local[j] = connToOneTermIDVector_.size();
    }
  N_ERH_ErrorMgr::safeBarrier(pdsMgrPtr_->getPDSComm()->comm());
  //for (j=0; j < procCnt; j++)
  //  comm.sumAll(&connoneterm_local[j],&connoneterm_sum[j],procCnt);

  comm.sumAll(&connoneterm_local[0],&connoneterm_sum[0],procCnt);

  //Broadcast the nodes to proc 0....
  if (procID != 0)
  {
    std::vector<std::string>::iterator connoneit = connToOneTermIDVector_.begin();
    std::vector<std::string>::iterator connoneend = connToOneTermIDVector_.end();
    int templength2;

    while (connoneit != connoneend)
    {
      templength2 = (*connoneit).length();
      comm.send(&templength2,1,0);
      comm.send(&(*connoneit)[0],templength2,0);
      connoneit++;
    }
  }
  //....and let proc 0 receive.
  else
  {
    int templength2;
    for (j=1; j < procCnt; j++)
    {
      for (int i=0; i < connoneterm_sum[j]; i++)
      {
	comm.recv(&templength2,1,j);
	std::string tempstring2;
	tempstring2.resize(templength2);
	comm.recv(&tempstring2[0],templength2,j);
	connToOneTermIDVector_.push_back(tempstring2);
      }
    }
  }
  N_ERH_ErrorMgr::safeBarrier(pdsMgrPtr_->getPDSComm()->comm());

#endif

  //Here's where we make the function calls to append resistors to nodes with
  //only one terminal connection and nodes which have no DC path to ground.
  //In parallel, we only need to do this on proc 0.

#ifdef Xyce_PARALLEL_MPI
  if (procID == 0)
  {

#endif
    bool netlistcopy = commandLine_.getHangingResistor().getNetlistCopy();
    std::string netlistFile("");
    if (commandLine_.getArgumentValue("netlist") != "")
    {
      netlistFile = commandLine_.getArgumentValue("netlist");
    }
    bool oneTermNotNoDCPath = true;
    //We use this boolean to print a different banner depending upon whether
    //resistors are being added because they are connected to only one device
    //terminal, or if they are being added because they have no DC path to
    //ground.


    //append resistors to nodes with only one terminal connection.

    if (netlistcopy && !connToOneTermIDVector_.empty())
    {
      std::string onetermres(commandLine_.getHangingResistor().getOneTermRes());
      topoPtr_->addResistors(connToOneTermIDVector_,onetermres,netlistFile,
                oneTermNotNoDCPath);
    }


    //append resistors to nodes with no dc path to ground.
    if (netlistcopy && !noDCPathIDVector_.empty())
    {
     std::string nodcpathres(commandLine_.getHangingResistor().getNoDCPathRes());
     topoPtr_->addResistors(noDCPathIDVector_,nodcpathres,netlistFile,
           !oneTermNotNoDCPath);

    }

    //if we've requested to produce a netlist file with resistors between
    //dangling nodes and ground, but it turns out that there aren't any
    //dangling nodes, we just need to add a ".END" to the end of the netlist
    //file copy that we produced in the I/O Interface Package (not critical).
    if (netlistcopy)
    {
      topoPtr_->appendEndStatement(netlistFile);
    }

#if defined(Xyce_PARALLEL_MPI)
  }
#endif

  //Don't remove this!!!!!!!!
  return true;
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::testVoltageNodeConnectivity
// Purpose       : testing of voltage node connectivity for problems
// Special Notes : initially, just warn if a voltage node has only 1 connection
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/31/03
//-----------------------------------------------------------------------------
bool TopoLSUtil::testVoltageNodeConnectivity_()
{
  std::map<int,std::string> cName;

  std::list<CktNode*>::iterator it_cnL =
    topoPtr_->orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator end_cnL =
    topoPtr_->orderedNodeListPtr_->end();

#if defined(Xyce_PARALLEL_MPI)
  N_PDS_Comm & comm = *(pdsMgrPtr_->getPDSComm());
  int procCnt = comm.numProc();
  int procID = comm.procID();
  int m, n;
  int proc;
#endif
  int i, j, k;
  int max_gid=0;
  int gid, num_gid;
  std::list<int> gidList;
  std::list<int> svGIDList;
  std::list<int> procList;
  std::list<NodeID> idList;
  std::vector<int> a,b;

#ifdef DEBUG_TOPO_DIAGS
#if defined(Xyce_PARALLEL_MPI)
  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  for (i=0 ; i<procCnt ; ++i)
  {
    if (i == procID)
    {
      Xyce::dout() << std::endl << "In TopoLSUtil::testVoltageNodeConnectivity_:"
       " Processing data on PE: " << i << std::endl;
#endif

      for( ; it_cnL != end_cnL; ++it_cnL )
      {
        const std::string & id = ((*it_cnL)->get_id()).first;
        int type = (*it_cnL)->type();
        int owned = (*it_cnL)->get_IsOwned();
        int proc = (*it_cnL)->get_ProcNum();
        gid = (*it_cnL)->get_gID();
        Xyce::dout() << id << " : type: " << type << " owned: " << owned
             << " proc: " << proc << " gid: " << gid << " connections:";
        gidList.clear();
        svGIDList.clear();
        procList.clear();
        idList.clear();
        topoPtr_->mainGraphPtr_->returnAdjNodes
          (gid, gidList, svGIDList, procList, idList);
        k = 0;
        std::list<int>::iterator gl = gidList.begin();
        std::list<int>::iterator end_gl = gidList.end();
        for ( ;gl != end_gl; ++gl)
        {
          if (k++ >= 20)
          {
            Xyce::dout() << " [connections terminated because > 20]";
            break;
          }
          Xyce::dout() << " " << *gl;
        }
        if ((*it_cnL)->type() == _DNODE)
        {
          const std::list<int> & GIDs = (*it_cnL)->get_ExtSolnVarGIDList();
          Xyce::dout() << " Ext:";
          std::list<int>::const_iterator ext_i=GIDs.begin();
          std::list<int>::const_iterator end_i= GIDs.end();
          for ( ; ext_i!=end_i; ++ext_i)
          {
            Xyce::dout() << " " << *ext_i;
          }
          Xyce::dout() << "  Connected leads:";
          std::vector<int> lead_conn;
          (*it_cnL)->leadConnect(lead_conn);
          for (j=0 ; j<(int)lead_conn.size() ; ++j)
          {
            Xyce::dout() << " " << lead_conn[j];
          }
        }
        Xyce::dout() << std::endl;
      }
#if defined(Xyce_PARALLEL_MPI)
    }
    comm.barrier();
  }
#endif

#endif

  std::vector<int> candidates;
  std::map<int,int> my_vnodes;
  std::map<int, std::vector<int> > gid_map;
  std::map<int, std::vector<int> >::iterator gm_i;
  std::map<int, std::vector<int> >::iterator gm_end;

#if defined(Xyce_PARALLEL_MPI)
  std::vector<int> comm_pe;
  std::map< int, std::set<int> > comm_gid;
  std::map< int, std::set<int> >::iterator cg_i;
  std::map< int, std::set<int> >::iterator cg_end;
  std::set<int>::iterator vm_i;
  std::set<int>::iterator vm_end;

  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      if (gid >= 0 && !(*it_cnL)->get_IsOwned())
      {
        proc = (*it_cnL)->get_ProcNum();
        comm_gid[proc].insert(gid);
      }
    }
  }

  int num_proc = comm_gid.size();
  int max_num_proc = 0;
  comm.maxAll (&num_proc, &max_num_proc, 1);

  int tmpSize = procCnt*max_num_proc;
  std::vector<int> pcomm    (tmpSize,0);
  std::vector<int> pcomm_all(tmpSize,0);
  std::vector<int> ncomm    (tmpSize,0);
  std::vector<int> ncomm_all(tmpSize,0);

  i = 0;
  cg_i   = comm_gid.begin();
  cg_end = comm_gid.end() ;
  for (; cg_i!=cg_end; ++cg_i)
  {
    pcomm[procID*max_num_proc+i] = (*cg_i).first+1;
    ncomm[procID*max_num_proc+i] = (*cg_i).second.size();
    i++;
  }
  comm.sumAll(&pcomm[0], &pcomm_all[0], tmpSize);
  comm.sumAll(&ncomm[0], &ncomm_all[0], tmpSize);

  std::vector<std::vector<int> > buf;
  std::vector<int> sendBuf;
  std::vector<int> src;
  int numMsg = 0;

  // Allocate receive buffers
  for (i=0 ; i<procCnt*max_num_proc ; ++i)
  {
    if (pcomm_all[i] == procID+1)
    {
      j = i/max_num_proc;
      k = ncomm_all[i];
      buf.resize(numMsg+1);
      buf[numMsg].resize(k);
      src.push_back(j);
      numMsg++;
    }
  }

  // Issue receives
  for (i=0 ; i<numMsg ; ++i)
  {
    comm.iRecv (&buf[i][0], buf[i].size(), src[i]);
  }

  // Issue sends
  cg_i=comm_gid.begin();
  for ( ; cg_i!=cg_end; ++cg_i)
  {
    i = (*cg_i).second.size();
    sendBuf.resize(i);
    j = 0;
    vm_i=(*cg_i).second.begin();
    vm_end=(*cg_i).second.end();
    for ( ; vm_i !=vm_end; ++vm_i)
    {
      sendBuf[j++] = *vm_i;
    }
    comm.send (&sendBuf[0], i, (*cg_i).first);
  }

  // Wait for messages to complete
  comm.waitAll();

  // Use received data
  for (i=0 ; i<numMsg ; ++i)
  {
    k = buf[i].size();
    for (m=0 ; m<k ; ++m)
      comm_gid[src[i]].insert(buf[i][m]);
  }

  int n_bufs = comm_gid.size();
  std::vector<int> buf_dest(n_bufs), buf_len(n_bufs);
  std::vector<int *> buf_in(n_bufs), buf_out(n_bufs);
  std::vector< std::vector<int> > buf_gid(n_bufs);

  int buf_tot = 0;
  cg_i   = comm_gid.begin();
  cg_end = comm_gid.end();
  for ( ; cg_i!=cg_end; ++cg_i)
  {
    buf_tot += (*cg_i).second.size();
  }

  std::vector<int> actual_buf_in(buf_tot);
  std::vector<int> actual_buf_out(buf_tot);

  i = 0;
  m = 0;
  n = 0;
  cg_i   = comm_gid.begin();
  cg_end = comm_gid.end();
  for ( ; cg_i!=cg_end; ++cg_i)
  {
    j = (*cg_i).second.size();
    buf_dest[i] = (*cg_i).first;
    buf_len[i] = j;
    buf_in[i] = &actual_buf_in[m];
    buf_out[i] = &actual_buf_out[m];
    buf_gid[i].resize(j);
    m += j;

    k = 0;
    vm_i   = (*cg_i).second.begin();
    vm_end = (*cg_i).second.end();
    for ( ; vm_i !=vm_end; ++vm_i)
    {
      buf_in[i][k] = 0;
      buf_out[i][k] = 0;
      buf_gid[i][k] = *vm_i;
      gid_map[*vm_i].push_back(n+k);
      ++k;
    }
    n += k;
    ++i;
  }

#ifdef DEBUG_TOPO_DIAGS
  for (i=0 ; i<procCnt ; ++i)
  {
    if (i == procID)
    {
      Xyce::dout() << "Connections for PE: " << i << std::endl;
      cg_i=comm_gid.begin();
      cg_end=comm_gid.end();
      for ( ; cg_i!=cg_end ; ++cg_i)
      {
        Xyce::dout() << "owner: " << (*cg_i).first << " :: ";
        k = 0;
        vm_i   = (*cg_i).second.begin();
        vm_end = (*cg_i).second.end();
        for ( ; vm_i !=vm_end; ++vm_i)
        {
          Xyce::dout() << "  " << *vm_i << "(" << gid_map[*vm_i][0];
          for (j=1 ; j<(int)gid_map[*vm_i].size() ; ++j)
          {
            Xyce::dout() << "," << gid_map[*vm_i][j];
          }
          Xyce::dout() << ")";
          if (++k%20 == 0)
            Xyce::dout() << std::endl;
        }
        Xyce::dout() << std::endl;
      }
    }
    comm.barrier();
  }
#endif

  for (i=0 ; i<buf_tot ; ++i)
  {
    actual_buf_out[i] = 0;
  }

  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      if (gid >= 0)
      {
        gidList.clear();
        svGIDList.clear();
        procList.clear();
        idList.clear();
        topoPtr_->mainGraphPtr_->returnAdjNodes
          (gid, gidList, svGIDList, procList, idList);
        num_gid = gidList.size();
        if (num_gid > 0 && !((*it_cnL)->get_IsOwned()))
        {
          if (gid_map.find(gid) != gid_map.end())
          {
            actual_buf_out[gid_map[gid][0]] = num_gid;
          }
        }
      }
    }
  }

  comm_boundaries (gid_map, actual_buf_in, actual_buf_out,
                   buf_len, buf_dest, buf_in, buf_out, 1);
#endif

  bool oneTerm = commandLine_.getHangingResistor().getOneTerm();

  candidates.clear();
  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      if (gid >= 0)
      {
        gidList.clear();
        svGIDList.clear();
        procList.clear();
        idList.clear();
        topoPtr_->mainGraphPtr_->returnAdjNodes
          (gid, gidList, svGIDList, procList, idList);
        num_gid = gidList.size();
#if defined(Xyce_PARALLEL_MPI)
        if (num_gid <= 1 && (*it_cnL)->get_IsOwned())
        {
          if (gid_map.find(gid) == gid_map.end() ||
             (gid_map.find(gid) != gid_map.end() &&
              num_gid + actual_buf_in[gid_map[gid][0]] <= 1))
          {
#else
        if (num_gid <= 1)
        {
          {                                                  // }}
#endif
            candidates.push_back(gid);
            cName[gid] = (*it_cnL)->get_id();

            if(oneTerm)
            {
              (*it_cnL)->setTrueConnToOneTermVar();
            }
          }
        }
      }
    }
  }

  outputTopoWarnings (candidates, cName,
  std::string("connected to only 1 device Terminal"));


  std::map<int, int> gid_pos;
  std::map<int, int>::iterator gp_i;
  std::map<int, int>::iterator gp_end;
  candidates.clear();

  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  int num_nodes = 0;
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      if (gid >= 0)
        gid_pos[gid] = num_nodes++;
    }
  }
  std::vector<int> node_val (num_nodes,0);
  std::vector<int> ext_gid  (gid_pos.size(),0);

  gp_i   = gid_pos.begin();
  gp_end = gid_pos.end();
  for ( ; gp_i != gp_end; ++gp_i)
  {
    node_val[(*gp_i).second] = (*gp_i).first;
  }
#if defined(Xyce_PARALLEL_MPI)
  i = 0;
  gm_i   = gid_map.begin();
  gm_end = gid_map.end();
  for ( ; gm_i != gm_end; ++gm_i)
  {
    if (gid_pos.find((*gm_i).first) != gid_pos.end())
    {
      ext_gid[i] = gid_pos[(*gm_i).first];
      ++i;
    }
    else
    {
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0, std::string(
      "TopoLSUtil::testVoltageNodeConnectivity_: External GID not found internally" ) );
    }
  }
#endif

  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  cName.clear();
#ifdef DEBUG_TOPO_DIAGS
  Xyce::dout() << std::endl;
#endif
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _DNODE)
    {
      int gid = (*it_cnL)->get_gID();
      gidList.clear();
      svGIDList.clear();
      procList.clear();
      idList.clear();
      topoPtr_->mainGraphPtr_->returnAdjNodes
        (gid, gidList, svGIDList, procList, idList);
      std::list<int>::iterator gl = gidList.begin();
      const std::list<int> & solnList = (*it_cnL)->get_ExtSolnVarGIDList();
      std::vector<int> GIDs;

      std::list<int>::const_iterator sol_i   = solnList.begin();
      std::list<int>::const_iterator sol_end = solnList.end();
      for ( ; sol_i!=sol_end; ++sol_i)
      {
        if (*sol_i == -1)
        {
          GIDs.push_back(-1);
        }
        else
        {
          GIDs.push_back(*(gl++));
        }
      }
      std::vector<int> lead_conn;
      (*it_cnL)->leadConnect(lead_conn);
      i = lead_conn.size();
      for (k=0 ; k<(int)lead_conn.size() ; ++k)
      {
        if (lead_conn[k] == 0)
        {
          if (GIDs[k] >= 0)
          {
            node_val[gid_pos[GIDs[k]]] = -1;
          }
          i--;
        }
      }
      for (j=1 ; j<10 ; ++j)
      {
        int last = -100;
        for (k=0 ; k<(int)lead_conn.size() ; ++k)
        {
          if (lead_conn[k] == j)
          {
            i--;
            if (last == -100)
            {
              last = GIDs[k];
            }
            else
            {
              if (last == -1)
              {
                if (GIDs[k] >= 0)
                  node_val[gid_pos[GIDs[k]]] = -1;
              }
              else
              {
                if (GIDs[k] >= 0)
                {
                  a.push_back(gid_pos[last]);
                  b.push_back(gid_pos[GIDs[k]]);
                }
                else
                  node_val[gid_pos[last]] = -1;
              }
              if (GIDs[k] < last)
                last = GIDs[k];
            }
          }
        }
        if (i == 0)
        {
          break;
        }
      }
      if (i != 0)
      {
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
        std::string( " Connectivity checker: lead index not found.  This checker is limited to devices with less than ten leads.  If you have gotten this error, and the circuit contains a device with more connections than this, you may need to run Xyce with this diagnostic turned off, via:\n  .options topology CHECK_CONNECTIVITY=0 " ) );
      }
    }
  }
  int n_pair = a.size();

  bool same_local, same = false;
  while (!same)
  {
    same = true;
    same_local = false;
    while (!same_local)
    {
      same_local = true;
      for (i=0 ; i<(int)a.size() ; ++i)
      {
        if (node_val[a[i]] < node_val[b[i]])
        {
          node_val[b[i]] = node_val[a[i]];
          same_local = false;
        }
        else if (node_val[b[i]] < node_val[a[i]])
        {
          node_val[a[i]] = node_val[b[i]];
          same_local = false;
        }
      }
      if (!same_local)
        same = false;
    }
#if defined(Xyce_PARALLEL_MPI)
    i = 0;
    gm_i   = gid_map.begin();
    gm_end = gid_map.end();
    for ( ; gm_i != gm_end; ++gm_i)
    {
      actual_buf_out[(*gm_i).second[0]] = node_val[ext_gid[i]];
      ++i;
    }
    comm_boundaries
      (gid_map, actual_buf_in, actual_buf_out,
       buf_len, buf_dest, buf_in, buf_out, 2);
    i = 0;
    gm_i   = gid_map.begin();
    gm_end = gid_map.end();
    for ( ; gm_i != gm_end; ++gm_i)
    {
      if (node_val[ext_gid[i]] > actual_buf_in[(*gm_i).second[0]])
      {
        node_val[ext_gid[i]] = actual_buf_in[(*gm_i).second[0]];
        same = false;
      }
      ++i;
    }
    if (same)
    {
      j = 1;
    }
    else
    {
      j = 0;
    }
    comm.sumAll(&j, &k, 1);
    if (k == procCnt)
    {
      same = true;
    }
    else
    {
      same = false;
    }
#endif
  }

  bool noDCPath = commandLine_.getHangingResistor().getNoDCPath();

  cName.clear();
  candidates.clear();
  for (i=0 ; i<(int)node_val.size() ; ++i)
  {
    if (node_val[i] != -1)
    {
      it_cnL = topoPtr_->orderedNodeListPtr_->begin();
      for( ; it_cnL != end_cnL; ++it_cnL )
      {
        if ((*it_cnL)->type() == _VNODE)
        {
          gid = (*it_cnL)->get_gID();
          if (gid >= 0 && node_val[gid_pos[gid]] >= 0)
          {
            candidates.push_back(gid);
            const std::string & id = (*it_cnL)->get_id();
            cName[gid] = id;
            if (!(*it_cnL)->getConnToOneTermVar() && noDCPath)
                                                //We let 'connected to one
            {                                   //terminal' take precedence
              (*it_cnL)->setTrueNoDCPathVar();  //over 'no DC path.'  (Don't
            }                                   //want to label a node as
          }                                     //*both* being connected to only one
        }                                       //terminal and having no DC path to ground
      }                                         //since this will cause the addition of two
      break;                                    //resistors instead of just one.
    }
  }

  outputTopoWarnings (candidates, cName,
      std::string("does not have a DC path to ground"));

  return true;


  // This is the old non-scalable method:
  candidates.clear();
  cName.clear();

  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  j = 0;
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      if (gid >= 0)
      {
        if (gid > j)
          j = gid;
        gidList.clear();
        svGIDList.clear();
        procList.clear();
        idList.clear();
        topoPtr_->mainGraphPtr_->returnAdjNodes(gid, gidList, svGIDList, procList, idList);
        num_gid = gidList.size();
        if ((*it_cnL)->get_IsOwned())
        {
          if (num_gid == 1)
          {
            candidates.push_back(gid);
            cName[gid] = (*it_cnL)->get_id();
          }
        }
        else
        {
          my_vnodes[gid] = num_gid;
        }
      }
    }
  }
  max_gid = j;

#if defined(Xyce_PARALLEL_MPI)
  int nc, nc_base, nc_all;
  nc = candidates.size();
  comm.scanSum(&nc, &nc_base, 1);
  comm.sumAll(&nc, &nc_all, 1);
  comm.maxAll (&j, &max_gid, 1);
  nc_base -= nc;

  std::vector<int> candidates_all(nc_all,0);
  std::vector<int> candidates_sum(nc_all,0);

  for (i=0 ; i<nc ; ++i)
    candidates_all[nc_base+i] = candidates[i];
  comm.sumAll(&candidates_all[0], &candidates_sum[0], nc_all);

  // Now candidates_sum has all of the candidate gids (based on a single node
  // connection on one PE)

  for (i=0 ; i<nc_all ; ++i)
  {
    candidates_all[i] = candidates_sum[i];
    if (i < nc_base || i >= nc_base+nc)
    {
      j = candidates_sum[i];
      if (my_vnodes.find(j) != my_vnodes.end())
      {
        candidates_all[i] = -1;
      }
    }
  }
  comm.minAll(&candidates_all[0], &candidates_sum[0], nc_all);

  // Now candidates_sum has all of the verified gids

  k = 0;
  for (i=0 ; i<nc ; ++i)
    if (candidates_sum[nc_base+i] >= 0)
      candidates[k++] = candidates[i];
  candidates.resize(k);
#endif

  outputTopoWarnings (candidates, cName, std::string("CONNECTED to only 1 device Terminal"));

  std::vector<int> gmin (max_gid+1,0);

  for (i=0 ; i<=max_gid ; ++i)
  {
    gmin[i] = i;
  }
  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  cName.clear();
#ifdef DEBUG_TOPO_DIAGS
  Xyce::dout() << std::endl;
#endif
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if ((*it_cnL)->type() == _DNODE)
    {
      int gid = (*it_cnL)->get_gID();
      gidList.clear();
      svGIDList.clear();
      procList.clear();
      idList.clear();
      topoPtr_->mainGraphPtr_->returnAdjNodes(gid, gidList, svGIDList, procList, idList);
      std::list<int>::iterator gl = gidList.begin();
      const std::list<int> & solnList = (*it_cnL)->get_ExtSolnVarGIDList();
      std::vector<int> GIDs;

      std::list<int>::const_iterator sol_i=solnList.begin() ;
      std::list<int>::const_iterator sol_end=solnList.end();
      for ( ; sol_i!=sol_end; ++sol_i)
      {
        if (*sol_i == -1)
          GIDs.push_back(-1);
        else
          GIDs.push_back(*(gl++));
      }
      std::vector<int> lead_conn;
      (*it_cnL)->leadConnect(lead_conn);
      i = lead_conn.size();
      for (k=0 ; k<(int)lead_conn.size() ; ++k)
      {
        if (lead_conn[k] == 0)
        {
          a.push_back(GIDs[k]);
          b.push_back(-1);
          i--;
        }
      }
      for (j=1 ; j<10 ; ++j)
      {
        int last = -100;
        for (k=0 ; k<(int)lead_conn.size() ; ++k)
        {
          if (lead_conn[k] == j)
          {
            i--;
            if (last == -100)
            {
              last = GIDs[k];
            }
            else
            {
              if (GIDs[k] > last)
              {
                a.push_back(last);
                b.push_back(GIDs[k]);
#ifdef DEBUG_TOPO_DIAGS
                Xyce::dout() << "Pair: " << a[a.size()-1] << " : " << b[b.size()-1] << std::endl;
#endif
              }
              else if (last > GIDs[k])
              {
                a.push_back(GIDs[k]);
                b.push_back(last);
#ifdef DEBUG_TOPO_DIAGS
                Xyce::dout() << "Pair: " << a[a.size()-1] << " : " << b[b.size()-1] << std::endl;
#endif
              }
              last = GIDs[k];
            }
          }
        }
        if (i == 0)
          break;
      }
      if (i != 0)
      {
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
        std::string( " Connectivity checker: lead index not found.  This checker is limited to devices with less than ten leads.  If you have gotten this error, and the circuit contains a device with more connections than this, you may need to run Xyce with this diagnostic turned off, via:\n  .options topology CHECK_CONNECTIVITY=0 " ) );
      }
    }
    else if ((*it_cnL)->type() == _VNODE)
    {
      gid = (*it_cnL)->get_gID();
      if (gid >= 0)
      {
        const std::string & id = (*it_cnL)->get_id();
        cName[gid] = id;
      }
    }
  }
  n_pair = a.size();

#if defined(Xyce_PARALLEL_MPI)
  int diff, diff_g;
  diff_g = 1;

  std::vector<int> gmin_g (max_gid+1,0);
  std::vector<int> gmin_t (max_gid+1,0);

  while (diff_g > 0)
  {
#endif

    bool changed = true;
    while (changed)
    {
      changed = false;
      for (i=0 ; i<n_pair ; ++i)
      {
        if (a[i] == -1 && gmin[b[i]] >= 0)
        {
          gmin[b[i]] = -1;
          changed = true;
        }
        else if (gmin[a[i]] != gmin[b[i]])
        {
          if (gmin[a[i]] > gmin[b[i]])
          {
            if (a[i] >= 0)
            {
              changed = true;
              gmin[a[i]] = gmin[b[i]];
            }
          }
          else
          {
            if (b[i] >= 0)
            {
              changed = true;
              gmin[b[i]] = gmin[a[i]];
            }
          }
        }
      }
    }

#if defined(Xyce_PARALLEL_MPI)
    comm.minAll(&gmin[0], &gmin_g[0], max_gid+1);
    diff = 0;
    for (i=0 ; i<=max_gid ; ++i)
    {
      if (gmin[i] != gmin_g[i])
        diff = 1;
    }
    comm.sumAll(&diff, &diff_g, 1);
    gmin_t = gmin;
    gmin = gmin_g;
    gmin_g = gmin_t;
  }
#endif

  candidates.resize(0);
  std::map<int,std::string>::iterator cn = cName.begin() ;
  std::map<int,std::string>::iterator cn_end = cName.end() ;
  for ( ; cn != cn_end; ++cn)
  {
    if (gmin[cn->first] >= 0)
      candidates.push_back(cn->first);
  }

  outputTopoWarnings (candidates, cName, std::string("DOES NOT have a DC path to ground"));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::comm_boundaries
// Purpose       : communicate boundary node data for topological checks
// Special Notes : mode=1 is sum, mode=2 is min
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/05
//-----------------------------------------------------------------------------
void TopoLSUtil::comm_boundaries (std::map<int, std::vector<int> > & gid_map,
                      std::vector<int> & actual_buf_in, std::vector<int> & actual_buf_out,
                      std::vector<int> & buf_len, std::vector<int> & buf_dest,
                      std::vector<int *> & buf_in, std::vector<int *> & buf_out, int mode)

{
  N_PDS_Comm & comm = *(pdsMgrPtr_->getPDSComm());
  unsigned int i;
  unsigned int n_bufs = buf_len.size();
  std::map< int, std::map<int, bool> >::iterator cg_i;
  std::map<int, std::vector<int> >::iterator g_i = gid_map.begin() ;
  std::map<int, std::vector<int> >::iterator g_end = gid_map.end();
  for ( ; g_i != g_end; ++g_i)
  {
    if ((*g_i).second.size() > 1)
    {
      for (i=1 ; i<(*g_i).second.size() ; ++i)
        actual_buf_out[(*g_i).second[i]] = actual_buf_out[(*g_i).second[0]];
    }
  }

  for (i = 0 ; i < n_bufs ; ++i)
  {
    comm.iRecv (buf_in[i], buf_len[i], buf_dest[i]);
  }
  for (i = 0 ; i < n_bufs ; ++i)
  {
    comm.send (buf_out[i], buf_len[i], buf_dest[i]);
  }
  comm.waitAll();

  g_i = gid_map.begin();
  g_end = gid_map.end();
  for ( ; g_i != g_end; ++g_i)
  {
    if ((*g_i).second.size() > 1)
    {
      for (i=1 ; i<(*g_i).second.size() ; ++i)
      {
        if (mode == 1)
          actual_buf_in[(*g_i).second[0]] += actual_buf_in[(*g_i).second[i]];
        else if (mode == 2)
        {
          if (actual_buf_in[(*g_i).second[i]] < actual_buf_in[(*g_i).second[0]])
            actual_buf_in[(*g_i).second[0]] = actual_buf_in[(*g_i).second[i]];
        }
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::outputTopoWarnings
// Purpose       : Output warnings from node connectivity checks
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/21/05
//-----------------------------------------------------------------------------
void TopoLSUtil::outputTopoWarnings(std::vector<int> & candidates,
                                          std::map<int,std::string> & cName,
                                          std::string errMsg)
{
  std::set<std::string> warnings;

#if defined(Xyce_PARALLEL_MPI)
  N_PDS_Comm & comm = *(pdsMgrPtr_->getPDSComm());
  int procCnt = comm.numProc();
  int procID = comm.procID();

  std::vector<int> candidates_all (procCnt,0);
  std::vector<int> candidates_sum (procCnt,0);

  int bs;

  for (int i=0 ; i<procCnt ; ++i)
    candidates_all[i] = 0;
  candidates_all[procID] = candidates.size();
  comm.sumAll(&candidates_all[0], &candidates_sum[0], procCnt);
#else
  int procID = 0;
#endif
  if (candidates.size() > 0)
  {
    for (unsigned i=0 ; i < candidates.size() ; ++i)
    {
      if (procID == 0)
      {
        std::string msg("Voltage Node (" + cName[candidates[i]] + ") " + errMsg);
        warnings.insert(msg);
      }
#if defined(Xyce_PARALLEL_MPI)
      else
      {
        bs = cName[candidates[i]].size();
        comm.send (&bs, 1, 0);
        comm.send (&cName[candidates[i]][0], cName[candidates[i]].size(), 0);
      }
#endif
    }
  }

#if defined(Xyce_PARALLEL_MPI)
  if (procID == 0)
  {
    std::string bad;
    std::string buf;

    for (int i=1 ; i<procCnt ; ++i)
    {
      for (int j=0 ; j<candidates_sum[i] ; ++j)
      {
        comm.recv (&bs, 1, i);
        buf.resize(bs);
        comm.recv (&buf[0], bs, i);
        bad = buf;
        std::string msg("Voltage Node (" + bad + ") " + errMsg);
        warnings.insert(msg);
      }
    }
  }
#endif
  if (procID == 0 && !warnings.empty())
  {
    std::set<std::string>::iterator warn = warnings.begin();
    for ( ; warn != warnings.end() ; ++warn)
    {
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING_0, *warn );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::setupNodeGIDs
// Purpose       : Generate Ordering and Var GIDs.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/22/01
//-----------------------------------------------------------------------------
bool TopoLSUtil::setupNodeGIDs()
{
  //get lists of owned and boundary/ghost node GIDs
  topoPtr_->returnNodeGIDVec( nodeList_GID_ );
  topoPtr_->returnExternNodeGIDVec( nodeList_ExternGID_ );
  numLocalNodes_ = nodeList_GID_.size();

  //calculate base GID for this processors nodes, lex. ordering
  double tmpVar1, tmpVar2;
  tmpVar1 = numLocalNodes_;
  pdsMgrPtr_->getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseNodeGID_ = static_cast<int>(tmpVar2 - numLocalNodes_);
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalNodes_ = static_cast<int>(tmpVar2);

  //Setup temp global accessors
  //Use to get new node GIDs for boundary/ghost nodes
  //---------------------------
  nodeGlobalAccessorPtr_ = pdsMgrPtr_->createGlobalAccessor();
  nodeGlobalAccessorPtr_->registerExternGIDVector( nodeList_ExternGID_ );
  nodeGlobalAccessorPtr_->generateMigrationPlan();

  std::map<int,int> nodeGIDMap;
  for( int iSV = 0; iSV < numLocalNodes_; ++iSV )
    nodeGIDMap[ nodeList_GID_[ iSV ] ] = iSV + baseNodeGID_;

/*
Xyce::dout() << "Before" << std::endl;
Xyce::dout() << "---------------------" << std::endl;
Xyce::dout() << pdsMgrPtr_->getPDSComm()->procID() << " Int Map ";
for( std::map<int,int>::iterator iter = nodeGIDMap.begin();
     iter != nodeGIDMap.end(); ++iter )
  Xyce::dout() << " " << iter->first << " " << iter->second << std::endl;
Xyce::dout() << "---------------------" << std::endl;
*/

  std::map<int,int> externGIDMap;
  nodeGlobalAccessorPtr_->migrateIntArray( nodeGIDMap, externGIDMap );

/*
Xyce::dout() << "After" << std::endl;
Xyce::dout() << "---------------------" << std::endl;
Xyce::dout() << pdsMgrPtr_->getPDSComm()->procID() << " Int Map ";
for( std::map<int,int>::iterator iter = nodeGIDMap.begin();
     iter != nodeGIDMap.end(); ++iter )
  Xyce::dout() << " " << iter->first << " " << iter->second << std::endl;
Xyce::dout() << pdsMgrPtr_->getPDSComm()->procID() << " Ext Map ";
for( std::map<int,int>::iterator iter = externGIDMap.begin();
     iter != externGIDMap.end(); ++iter )
  Xyce::dout() << " " << iter->first << " " << iter->second << std::endl;
Xyce::dout() << "---------------------" << std::endl;
*/

//Xyce::dout() << "Before Node GIDs" << std::endl << *topoPtr_;

  //loop over nodes and reset all GIDs to new lex. ordering
  std::list<CktNode*>::iterator it_cnL, end_cnL;
  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  end_cnL = topoPtr_->orderedNodeListPtr_->end();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if( nodeGIDMap.find( (*it_cnL)->get_gID() ) != nodeGIDMap.end() )
    {
      (*it_cnL)->set_gID( nodeGIDMap[ (*it_cnL)->get_gID() ] );

#ifdef Xyce_DEBUG_TOPOLOGY
Xyce::dout() << "Node: " << (*it_cnL)->get_id() << " " << (*it_cnL)->get_gID() << std::endl;
#endif

    }
    else if( externGIDMap.find( (*it_cnL)->get_gID() ) != externGIDMap.end() )
      (*it_cnL)->set_gID( externGIDMap[ (*it_cnL)->get_gID() ] );
    else if( (*it_cnL)->get_gID() != -1 )
    {
       std::string err = "P" +  Teuchos::Utils::toString(pdsMgrPtr_->getPDSComm()->procID())
          + ": Node: " + (*it_cnL)->get_id() + ", global index ("
          + Teuchos::Utils::toString( (*it_cnL)->get_gID() ) + ") is NOT found!";
       N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0, err );
    }
  }

//Xyce::dout() << "Finished Node GIDs" << std::endl << *topoPtr_;

  //Reset global accessor for new lex. indexing
  topoPtr_->returnNodeGIDVec( nodeList_GID_ );
  topoPtr_->returnExternNodeGIDVec( nodeList_ExternGID_ );

  //Register Node Parallel Map for later lookups
  N_PDS_ParMap * nodeMap = pdsMgrPtr_->createParallelMap( numGlobalNodes_,
	                                                  numLocalNodes_,
                                                          nodeList_GID_ );
  pdsMgrPtr_->addParallelMap( "NODE", nodeMap );

  pdsMgrPtr_->addGlobalAccessor( "NODE" );
  delete nodeGlobalAccessorPtr_;
  nodeGlobalAccessorPtr_ = pdsMgrPtr_->getGlobalAccessor( "NODE" );
  nodeGlobalAccessorPtr_->registerExternGIDVector( nodeList_ExternGID_ );
  nodeGlobalAccessorPtr_->generateMigrationPlan();

  //regenerate graph's Node GID lookup map since GIDs have changed
  topoPtr_->regenerateGIDNodeMap();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::setupSolnAndStateGIDs
// Purpose       : Generate Ordering and Var GIDs.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/2/01
//-----------------------------------------------------------------------------
bool TopoLSUtil::setupSolnAndStateGIDs()
{
  //get soln and state var counts from all owned v-nodes and d-nodes
  std::vector<int> rowCountVec(numLocalNodes_);
  std::vector<int> stateCountVec(numLocalNodes_);
  std::vector<int> storeCountVec(numLocalNodes_);
  std::list<CktNode*>::iterator it_cnL, end_cnL;
  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  end_cnL = topoPtr_->orderedNodeListPtr_->end();
  int Loc = 0;
  for( ; it_cnL != end_cnL; ++it_cnL )
    if( (*it_cnL)->get_IsOwned() && (*it_cnL)->get_gID() != -1 )
    {
      rowCountVec[Loc] = (*it_cnL)->solnVarCount();
      stateCountVec[Loc] = (*it_cnL)->stateVarCount();
      storeCountVec[Loc] = (*it_cnL)->storeVarCount();
      ++Loc;
    }

  //get global sum of counts
  numLocalRows_ = 0;
  numLocalStateVars_ = 0;
  numLocalRows_ = accumulate( rowCountVec.begin(), rowCountVec.end(), 0 );
  numLocalStateVars_ = accumulate( stateCountVec.begin(), stateCountVec.end(), 0 );
  numLocalStoreVars_ = 0;
  numLocalStoreVars_ = accumulate( storeCountVec.begin(), storeCountVec.end(), 0 );

  //calculate base soln and state GIDs for this processor
  double tmpVar1, tmpVar2;
  tmpVar1 = numLocalRows_;
  pdsMgrPtr_->getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseRowGID_ = static_cast<int>(tmpVar2 - numLocalRows_);
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalRows_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalStateVars_;
  pdsMgrPtr_->getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseStateVarGID_ = static_cast<int>(tmpVar2 - numLocalStateVars_);
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalStateVars_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalStoreVars_;
  pdsMgrPtr_->getPDSComm()->scanSum( &tmpVar1, &tmpVar2, 1 );
  baseStoreVarGID_ = static_cast<int>(tmpVar2 - numLocalStoreVars_);
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalStoreVars_ = static_cast<int>(tmpVar2);

  //Use Global Accessors to get GIDs for boundary/ghost nodes
  std::map< int,std::vector<int> > rowGIDMap, stateGIDMap;
  std::map< int,std::vector<int> > externRowGIDMap, externStateGIDMap;
  std::map< int,std::vector<int> > storeGIDMap;
  std::map< int,std::vector<int> > externStoreGIDMap;

  std::vector<int> tmpVec;
  int currRowLoc = baseRowGID_;
  int currStateLoc = baseStateVarGID_;
  int currStoreLoc = baseStoreVarGID_;
  for( int i = 0; i < numLocalNodes_; ++i )
  {
    tmpVec.resize(rowCountVec[i]);
    iota( tmpVec.begin(), tmpVec.end(), currRowLoc );
    rowGIDMap[ nodeList_GID_[i] ] = tmpVec;
    currRowLoc += rowCountVec[i];

    tmpVec.resize(stateCountVec[i]);
    iota( tmpVec.begin(), tmpVec.end(), currStateLoc );
    stateGIDMap[ nodeList_GID_[i] ] = tmpVec;
    currStateLoc += stateCountVec[i];

    tmpVec.resize(storeCountVec[i]);
    iota( tmpVec.begin(), tmpVec.end(), currStoreLoc );
    storeGIDMap[ nodeList_GID_[i] ] = tmpVec;
    currStoreLoc += storeCountVec[i];
  }

  nodeGlobalAccessorPtr_->migrateIntVecs( rowGIDMap, externRowGIDMap );
  nodeGlobalAccessorPtr_->migrateIntVecs( stateGIDMap, externStateGIDMap );
  nodeGlobalAccessorPtr_->migrateIntVecs( storeGIDMap, externStoreGIDMap );

  //calculate number of boundary/ghost soln and state vars
  std::map< int,std::vector<int> >::iterator iterIVM = externRowGIDMap.begin();
  std::map< int,std::vector<int> >::iterator endIVM = externRowGIDMap.end();
  numExternRows_ = 0;
  for( ; iterIVM != endIVM; ++iterIVM )
    numExternRows_ += iterIVM->second.size();
  tmpVar1 = numExternRows_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalExternRows_ = static_cast<int>(tmpVar2);

  iterIVM = externStateGIDMap.begin();
  endIVM = externStateGIDMap.end();
  numExternStateVars_ = 0;
  for( ; iterIVM != endIVM; ++iterIVM )
    numExternStateVars_ += iterIVM->second.size();
  tmpVar1 = numExternStateVars_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalExternStateVars_ = static_cast<int>(tmpVar2);

  iterIVM = externStoreGIDMap.begin();
  endIVM = externStoreGIDMap.end();
  numExternStoreVars_ = 0;
  for( ; iterIVM != endIVM; ++iterIVM )
    numExternStoreVars_ += iterIVM->second.size();
  tmpVar1 = numExternStoreVars_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalExternStoreVars_ = static_cast<int>(tmpVar2);

  //loop over nodes and assign soln and state GIDs
  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  end_cnL = topoPtr_->orderedNodeListPtr_->end();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    if( rowGIDMap.count( (*it_cnL)->get_gID() ) )
    {
      (*it_cnL)->set_SolnVarGIDVec( rowGIDMap[ (*it_cnL)->get_gID() ] );
      (*it_cnL)->set_StateVarGIDVec( stateGIDMap[ (*it_cnL)->get_gID() ] );
      (*it_cnL)->set_StoreVarGIDVec( storeGIDMap[ (*it_cnL)->get_gID() ] );
    }
    else if( externRowGIDMap.count( (*it_cnL)->get_gID() ) )
    {
      (*it_cnL)->set_SolnVarGIDVec( externRowGIDMap[ (*it_cnL)->get_gID() ] );
      (*it_cnL)->set_StateVarGIDVec( externStateGIDMap[ (*it_cnL)->get_gID() ] );
      (*it_cnL)->set_StoreVarGIDVec( externStoreGIDMap[ (*it_cnL)->get_gID() ] );
    }
    else if( (*it_cnL)->get_gID() == -1 )
      (*it_cnL)->set_SolnVarGIDVec( std::vector<int>(1,-1) );
    else
    {
       std::string err = "P" +  Teuchos::Utils::toString(pdsMgrPtr_->getPDSComm()->procID())
          + ": Node: " + (*it_cnL)->get_id() + ", global index ("
          + Teuchos::Utils::toString( (*it_cnL)->get_gID() ) + ") is NOT found!";
       N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0, err );
    }
  }

  return true;

}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::generateRowColData
// Purpose       : Generate row/col data for linear solver.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
void TopoLSUtil::generateRowColData()
{
  int procID = pdsMgrPtr_->getPDSComm()->procID();

  double tmpVar1, tmpVar2;

  std::list<index_pair> tmpIPList;

  //--- extract ordered list of GIDs from node list
  topoPtr_->returnSVarGIDVec( rowList_GID_ );

  //--- extract list of variable types
  topoPtr_->returnVarTypeVec( rowList_VarType_ );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "Ordered List of GIDs extracted" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif

  //--- extract ordered list of extern GIDs from node list
  topoPtr_->returnExternSVarGIDVec( rowList_ExternGID_ );

  //--- add in dep soln var stuff
  if( !topoPtr_->depSolnGIDMap_.empty() )
  {
    std::map<int,int> tmpMap;
    std::map<int,int>::iterator iterIIM = topoPtr_->depSolnGIDMap_.begin();
    std::map<int,int>::iterator endIIM = topoPtr_->depSolnGIDMap_.end();
    for( ; iterIIM != endIIM; ++iterIIM )
      if( iterIIM->second != procID ) tmpMap.insert( *iterIIM );

    for( unsigned int i = 0; i < rowList_ExternGID_.size(); ++i )
      if( rowList_ExternGID_[i].second  != procID )
        tmpMap.insert( rowList_ExternGID_[i] );

    rowList_ExternGID_.resize( tmpMap.size() );
    iterIIM = tmpMap.begin();
    for( unsigned int i = 0; i < tmpMap.size(); ++iterIIM, ++i )
      rowList_ExternGID_[i] = *iterIIM;
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "Ordered List of Extern GIDs extracted" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif

  //--- extract ordered list of State Var GIDs from node list
  topoPtr_->returnStateVarGIDVec( rowList_StateGID_ );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "Ordered List of State GIDs extracted" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif

  //--- extract ordered list of extern state GIDs from node list
  topoPtr_->returnExternStateVarGIDVec( rowList_ExternStateGID_ );

  //--- add in dep state var stuff
  if( !topoPtr_->depStateGIDMap_.empty() )
  {
    std::map<int,int> tmpMap;
    std::map<int,int>::iterator iterIIM = topoPtr_->depStateGIDMap_.begin();
    std::map<int,int>::iterator endIIM = topoPtr_->depStateGIDMap_.end();
    for( ; iterIIM != endIIM; ++iterIIM )
      if( iterIIM->second != procID ) tmpMap.insert( *iterIIM );

    for( unsigned int i = 0; i < rowList_ExternStateGID_.size(); ++i )
      if( rowList_ExternStateGID_[i].second  != procID )
        tmpMap.insert( rowList_ExternStateGID_[i] );

    rowList_ExternStateGID_.resize( tmpMap.size() );
    iterIIM = tmpMap.begin();
    endIIM = tmpMap.end();
    for( int i = 0; iterIIM != endIIM; ++iterIIM, ++i )
      rowList_ExternStateGID_[i] = *iterIIM;
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "Ordered List of Extern State GIDs extracted" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif


  //--- extract ordered list of Store Var GIDs from node list
  topoPtr_->returnStoreVarGIDVec( rowList_StoreGID_ );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "Ordered List of Store GIDs extracted" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif

  //--- extract ordered list of extern store GIDs from node list
  topoPtr_->returnExternStoreVarGIDVec( rowList_ExternStoreGID_ );

  //--- add in dep store var stuff
  if( !topoPtr_->depStoreGIDMap_.empty() )
  {
    std::map<int,int> tmpMap;
    std::map<int,int>::iterator iterIIM = topoPtr_->depStoreGIDMap_.begin();
    std::map<int,int>::iterator endIIM = topoPtr_->depStoreGIDMap_.end();
    for( ; iterIIM != endIIM; ++iterIIM )
      if( iterIIM->second != procID ) tmpMap.insert( *iterIIM );

    for( unsigned int i = 0; i < rowList_ExternStoreGID_.size(); ++i )
      if( rowList_ExternStoreGID_[i].second  != procID )
        tmpMap.insert( rowList_ExternStoreGID_[i] );

    rowList_ExternStoreGID_.resize( tmpMap.size() );
    iterIIM = tmpMap.begin();
    endIIM = tmpMap.end();
    for( int i = 0; iterIIM != endIIM; ++iterIIM, ++i )
      rowList_ExternStoreGID_[i] = *iterIIM;
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "Ordered List of Extern Store GIDs extracted" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif


  //--- set numLocalRows_, numLocalStateVars_, and resize lists
  numLocalRows_ = rowList_GID_.size();
  numExternRows_ = rowList_ExternGID_.size();
  numLocalStateVars_ = rowList_StateGID_.size();
  numExternStateVars_ = rowList_ExternStateGID_.size();
  numLocalStoreVars_ = rowList_StoreGID_.size();
  numExternStoreVars_ = rowList_ExternStoreGID_.size();

  tmpVar1 = numLocalRows_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalRows_ = static_cast<int>(tmpVar2);

  tmpVar1 = numExternRows_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalExternRows_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalStateVars_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalStateVars_ = static_cast<int>(tmpVar2);

  tmpVar1 = numExternStateVars_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalExternStateVars_ = static_cast<int>(tmpVar2);

  tmpVar1 = numLocalStoreVars_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalStoreVars_ = static_cast<int>(tmpVar2);

  tmpVar1 = numExternStoreVars_;
  pdsMgrPtr_->getPDSComm()->sumAll( &tmpVar1, &tmpVar2, 1 );
  numGlobalExternStoreVars_ = static_cast<int>(tmpVar2);

  //--- block extern gids by processor to support aztec requirements
  if( numExternRows_ )
  {
    std::vector<int> externGIDs( numExternRows_ );
    std::vector<int> PIDs( numExternRows_ );
    for( int i = 0; i < numExternRows_; ++i )
    {
      externGIDs[i] = rowList_ExternGID_[i].first;
      PIDs[i] = rowList_ExternGID_[i].second;
    }

    Epetra_Util Util;
    int ** listPtr = new int *[1];
    listPtr[0] = &externGIDs[0];
    Util.Sort( true, numExternRows_, &PIDs[0], 0, 0, 1, listPtr );
    delete [] listPtr;

    for( int i = 0; i < numExternRows_; ++i )
      rowList_ExternGID_[i] = std::pair<int,int>( externGIDs[i], PIDs[i] );
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "TopoUtil Vals: " << std::endl;
  Xyce::dout() << numGlobalRows_ << " " << numLocalRows_ << std::endl;
  Xyce::dout() << numGlobalStateVars_ << " " << numLocalStateVars_ << std::endl;
  Xyce::dout() << numGlobalStoreVars_ << " " << numLocalStoreVars_ << std::endl;
#endif

  int numRows = numLocalRows_ + numExternRows_ + 1;

  //Build global to local map for speed
  GtoL_Map_.clear();
  for( int i = 0; i < numLocalRows_; ++i )
    GtoL_Map_[ rowList_GID_[i] ] = i;
  for( int i = 0; i < numExternRows_; ++i )
    GtoL_Map_[ rowList_ExternGID_[i].first ] = i + numLocalRows_;
  GtoL_Map_[ -1 ] = numLocalRows_ + numExternRows_;

  rowList_ColList_.resize( numRows );
  rowList_NumNZs_.resize( numRows );

  int rcCnt = rowList_ColList_.size();
  for( int i = 0; i < rcCnt; ++i ) rowList_ColList_[i].clear();

  if( commandLine_.getArgumentValue( "-dma" ) == "off" )
  {
    //    to compile a list of col indices for each row
    std::list<CktNode*>::iterator it_cnL, end_cnL;
    it_cnL = topoPtr_->orderedNodeListPtr_->begin();
    end_cnL = topoPtr_->orderedNodeListPtr_->end();
    for( ; it_cnL != end_cnL; ++it_cnL )
      if( (*it_cnL)->type() == _DNODE )
      {
        tmpIPList.clear();
        (*it_cnL)->getRowColPairs( tmpIPList );

        if( tmpIPList.size() > 0 )
        {
          std::list<index_pair>::iterator it_ipL = tmpIPList.begin();
          std::list<index_pair>::iterator it_ipL_end = tmpIPList.end();
          for( ; it_ipL != it_ipL_end; ++it_ipL )
          {
            rowList_ColList_[ GtoL_Map_[it_ipL->row] ].push_back
              ( it_ipL->col );
          }
        }
      }
  }
  else
  {
    //    to compile a list of col indices for each row
    std::list<CktNode*>::iterator it_cnL = topoPtr_->orderedNodeListPtr_->begin();
    std::list<CktNode*>::iterator end_cnL = topoPtr_->orderedNodeListPtr_->end();
    for( ; it_cnL != end_cnL; ++it_cnL )
    {
      if( (*it_cnL)->type() == _DNODE )
      {
        const std::vector< std::vector<int> > & stamp = (*it_cnL)->jacobianStamp();
        const std::list<int> & intGIDs = (*it_cnL)->get_SolnVarGIDList();
        const std::list<int> & extGIDs = (*it_cnL)->get_ExtSolnVarGIDList();
        const std::vector<int> & depGIDs = (*it_cnL)->get_DepSolnGIDJacVec();
        std::vector<int> gids( intGIDs.size() + extGIDs.size() + depGIDs.size() );
        copy( extGIDs.begin(), extGIDs.end(), gids.begin() );
        copy( intGIDs.begin(), intGIDs.end(), gids.begin() + extGIDs.size() );
        copy( depGIDs.begin(), depGIDs.end(), gids.begin() + extGIDs.size() + intGIDs.size() );

        assert( ( extGIDs.size() + intGIDs.size() ) == stamp.size() );

        for( unsigned int i = 0; i < stamp.size(); ++i )
        {
          int length = stamp[i].size();
          for( int j = 0; j < length; ++j )
          {
            if (stamp[i][j] < (int)gids.size())
              rowList_ColList_[ GtoL_Map_[ gids[i] ] ].push_back( gids[ stamp[i][j] ] );
          }
        }
      }
    }

  }


#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "Row List of NZ Cols extracted" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif

  numLocalNZs_ = 0;

#ifdef Xyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT
    //add in diagonal for homotopy support
    for( int i = 0; i < numLocalRows_; ++i )
      rowList_ColList_[i].push_back( rowList_GID_[i] );
    if( commandLine_.getArgumentValue( "-dma" ) != "off" )
      for( int i = 0; i < numExternRows_; ++i )
        rowList_ColList_[i+numLocalRows_].push_back( rowList_ExternGID_[i].first );
#endif

  //--- sort list of col indices and get rid of redundancies
  //    generate num of NZs data
  for( int i = 0; i < numRows; ++i )
  {
    rowList_ColList_[i].sort();
    rowList_ColList_[i].unique();

    rowList_NumNZs_[i] = rowList_ColList_[i].size();

    numLocalNZs_ += rowList_NumNZs_[i];
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "Row List of NZ Cols cleaned up" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << *this;
#endif

}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::setupGlobalAccessors
// Purpose       : Register External GIDs and generate Migration Plans.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/00
//-----------------------------------------------------------------------------
bool TopoLSUtil::setupGlobalAccessors()
{
  solnGlobalAccessorPtr_ = pdsMgrPtr_->createGlobalAccessor( "SOLUTION" );
  stateGlobalAccessorPtr_ = pdsMgrPtr_->createGlobalAccessor( "STATE" );
  storeGlobalAccessorPtr_ = pdsMgrPtr_->createGlobalAccessor( "STORE" );

  solnGlobalAccessorPtr_->registerExternGIDVector( rowList_ExternGID_ );
  stateGlobalAccessorPtr_->registerExternGIDVector( rowList_ExternStateGID_ );
  storeGlobalAccessorPtr_->registerExternGIDVector( rowList_ExternStoreGID_ );

  solnGlobalAccessorPtr_->generateMigrationPlan();
  stateGlobalAccessorPtr_->generateMigrationPlan();
  storeGlobalAccessorPtr_->generateMigrationPlan();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::reorderGIDs
// Purpose       : Reorder GIDs using based on orderedNodeList.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/00
//-----------------------------------------------------------------------------
void TopoLSUtil::reorderGIDs()
{
  std::map<int,int> GIDMap, stateGIDMap;
  std::map<int,int> storeGIDMap;

  for( int iSV = 0; iSV < numLocalRows_; ++iSV )
  {
    GIDMap[ rowList_GID_[ iSV ] ] = iSV + baseRowGID_;
  }

  for( int iStV = 0; iStV < numLocalStateVars_; ++iStV )
  {
    stateGIDMap[ rowList_StateGID_[ iStV ] ] = iStV + baseStateVarGID_;
  }
  for( int iStV = 0; iStV < numLocalStoreVars_; ++iStV )
  {
    storeGIDMap[ rowList_StoreGID_[ iStV ] ] = iStV + baseStoreVarGID_;
  }

  std::map<int,int> externGIDMap, externStateGIDMap;
  solnGlobalAccessorPtr_->migrateIntArray( GIDMap, externGIDMap );
  stateGlobalAccessorPtr_->migrateIntArray( stateGIDMap, externStateGIDMap );

  std::map<int,int> externStoreGIDMap;
  storeGlobalAccessorPtr_->migrateIntArray( storeGIDMap, externStoreGIDMap );

  std::map<int,int>::iterator it_iiM;
  std::map<int,int>::iterator it_iiM_end;

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "GIDMap" << std::endl;
  it_iiM = GIDMap.begin();
  it_iiM_end = GIDMap.end();
  for( ; it_iiM != it_iiM_end; ++it_iiM )
    Xyce::dout() << it_iiM->first << " " << it_iiM->second << std::endl;
  Xyce::dout() << "externGIDMap" << std::endl;
  it_iiM = externGIDMap.begin();
  it_iiM_end = externGIDMap.end();
  for( ; it_iiM != it_iiM_end; ++it_iiM )
    Xyce::dout() << it_iiM->first << " " << it_iiM->second << std::endl;

  Xyce::dout() << "stateGIDMap" << std::endl;
  it_iiM = stateGIDMap.begin();
  it_iiM_end = stateGIDMap.end();
  for( ; it_iiM != it_iiM_end; ++it_iiM )
    Xyce::dout() << it_iiM->first << " " << it_iiM->second << std::endl;
  Xyce::dout() << "externStateGIDMap" << std::endl;
  it_iiM = externStateGIDMap.begin();
  it_iiM_end = externStateGIDMap.end();
  for( ; it_iiM != it_iiM_end; ++it_iiM )
    Xyce::dout() << it_iiM->first << " " << it_iiM->second << std::endl;
  Xyce::dout() << std::endl;

  Xyce::dout() << "storeGIDMap" << std::endl;
  it_iiM = storeGIDMap.begin();
  it_iiM_end = storeGIDMap.end();
  for( ; it_iiM != it_iiM_end; ++it_iiM )
    Xyce::dout() << it_iiM->first << " " << it_iiM->second << std::endl;
  Xyce::dout() << "externStoreGIDMap" << std::endl;
  it_iiM = externStoreGIDMap.begin();
  it_iiM_end = externStoreGIDMap.end();
  for( ; it_iiM != it_iiM_end; ++it_iiM )
    Xyce::dout() << it_iiM->first << " " << it_iiM->second << std::endl;
  Xyce::dout() << std::endl;

#endif

  std::list<int> tmpIList;

  std::list<CktNode*>::iterator it_cnL, end_cnL;
  it_cnL = topoPtr_->orderedNodeListPtr_->begin();
  end_cnL = topoPtr_->orderedNodeListPtr_->end();
  for( ; it_cnL != end_cnL; ++it_cnL )
  {
    tmpIList.clear();

    std::list<int>::const_iterator it_iL, end_iL;
    it_iL = (*it_cnL)->get_SolnVarGIDList().begin();
    end_iL = (*it_cnL)->get_SolnVarGIDList().end();
    for( ; it_iL != end_iL ; ++it_iL )
    {
      if( *it_iL != -1 )
      {
        if( GIDMap.find( *it_iL ) != GIDMap.end() )
        {
          tmpIList.push_back( GIDMap[ *it_iL ] );
        }
        else
        {
          if( externGIDMap.find( *it_iL ) != externGIDMap.end() )
          {
            tmpIList.push_back( externGIDMap[ *it_iL ] );
          }
          else
          {
            N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
              std::string( " Global index not found " ) );
          }
        }
      }
      else
      {
        tmpIList.push_back( -1 );
      }
    }

    (*it_cnL)->set_SolnVarGIDList( tmpIList );

    tmpIList.clear();

    it_iL = (*it_cnL)->get_StateVarGIDList().begin();
    end_iL = (*it_cnL)->get_StateVarGIDList().end();
    for( ; it_iL != end_iL ; ++it_iL )
    {
      if( *it_iL != -1 )
      {
        if( stateGIDMap.find( *it_iL ) != stateGIDMap.end() )
        {
          tmpIList.push_back( stateGIDMap[ *it_iL ] );
        }
        else
        {
          if( externStateGIDMap.find( *it_iL ) != externStateGIDMap.end() )
          {
            tmpIList.push_back( externStateGIDMap[ *it_iL ] );
          }
          else
          {
            N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
              std::string( " State Global index not found " ) );
          }
        }
      }
      else
      {
        tmpIList.push_back( -1 );
      }
    }

    (*it_cnL)->set_StateVarGIDList( tmpIList );

    tmpIList.clear();

    it_iL = (*it_cnL)->get_StoreVarGIDList().begin();
    end_iL = (*it_cnL)->get_StoreVarGIDList().end();
    for( ; it_iL != end_iL ; ++it_iL )
    {
      if( *it_iL != -1 )
      {
        if( storeGIDMap.find( *it_iL ) != storeGIDMap.end() )
        {
          tmpIList.push_back( storeGIDMap[ *it_iL ] );
        }
        else
        {
          if( externStoreGIDMap.find( *it_iL ) != externStoreGIDMap.end() )
          {
            tmpIList.push_back( externStoreGIDMap[ *it_iL ] );
          }
          else
          {
            N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
              std::string( " Store Global index not found " ) );
          }
        }
      }
      else
      {
        tmpIList.push_back( -1 );
      }
    }

    (*it_cnL)->set_StoreVarGIDList( tmpIList );
  }

}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::registerTimeOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/29/05
//-----------------------------------------------------------------------------
bool TopoLSUtil::registerTimeOptions(const N_UTL_OptionBlock & OB)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool TopoLSUtil::registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  std::string netListFile("");
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }

  pkgOptMgrPtr_->submitRegistration( "TIMEINT", netListFile,
                                  new TopoLSUtil_TimeOptionsReg( *this ) );

  pkgOptMgrPtr_->submitRegistration( "TOPOLOGY", netListFile,
                                  new TopoLSUtil_OptionsReg( *this ) );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : TopoLSUtil::registerOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/05
//-----------------------------------------------------------------------------
bool TopoLSUtil::registerOptions(const N_UTL_OptionBlock & OB)
{
  std::list<N_UTL_Param>::const_iterator it_tpL;
  std::list<N_UTL_Param>::const_iterator first = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag()=="CHECK_CONNECTIVITY")
    {
      checkConnectivity_ = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if(it_tpL->uTag()=="SUPERNODE")
    {
      supernode_ = static_cast<bool>(it_tpL->getImmutableValue<bool>());
    }
    else if(it_tpL->uTag()=="OUTPUTNAMESFILE")
    {
      namesFile_ = static_cast<bool>(it_tpL->getImmutableValue<bool>());
    }
  }

  return true;
}

} // namespace Topo
} // namespace Xyce
