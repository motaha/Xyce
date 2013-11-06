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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_TOP_Directory.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/02/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>
#include <N_TOP_Misc.h>

#include <iostream>

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif

// ---------- Xyce Includes --------------
#include <N_TOP_Directory.h>

#include <N_TOP_Topology.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_ParNode.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_Misc.h>

#include <N_UTL_Functors.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_PDS_Directory.h>
#include <N_PDS_Migrate.h>

//-----------------------------------------------------------------------------
// Class         : N_TOP_DirectoryData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/26/04
//-----------------------------------------------------------------------------
struct N_TOP_DirectoryData
{
  typedef RCP<N_TOP_ParNode> NodePtr;
  typedef map<NodeID,NodePtr> NodeContainer;
                                                                                           
  typedef Xyce::Parallel::Hash<NodeID> NodeIDHash;
  typedef Xyce::Parallel::Migrate<NodeID,N_TOP_ParNode> NodeMigrate;

  typedef Xyce::Parallel::Directory< NodeID,
                                     N_TOP_ParNode,
                                     NodeIDHash,
                                     NodeContainer,
                                     NodeMigrate >  NodeDir;

  N_TOP_DirectoryData( N_PDS_Comm & comm )
  : hash(comm.numProc()),
    migrate(comm),
    directory(migrate,hash)
  {}

  ~N_TOP_DirectoryData() {}

  NodeIDHash hash;
  NodeMigrate migrate;

  NodeDir directory;
};
  

//-----------------------------------------------------------------------------
// Function      : N_TOP_Directory::~N_TOP_Directory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/02/01
//-----------------------------------------------------------------------------
N_TOP_Directory::~N_TOP_Directory()
{
  if( data_ ) delete data_;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_Directory::generateDirectory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/02/01
//-----------------------------------------------------------------------------
bool N_TOP_Directory::generateDirectory()
{
  int procID = pdsMgr_->getPDSComm()->procID();

  if( data_ ) delete data_;
  data_ = new N_TOP_DirectoryData( *(pdsMgr_->getPDSComm()) );

  N_TOP_DirectoryData::NodeContainer nodes;

  list<N_TOP_CktNode*>::iterator iterCN = topMgr_->orderedNodeListPtr_->begin();
  list<N_TOP_CktNode*>::iterator endCN = topMgr_->orderedNodeListPtr_->end();
  
  for( ; iterCN != endCN; ++iterCN )
    if( (*iterCN)->get_IsOwned() && (*iterCN)->get_gID() != -1 )
    {
      N_TOP_DirectoryData::NodePtr new_node(
              new N_TOP_ParNode( NodeID( (*iterCN)->get_id(), (*iterCN)->get_gID() ),
                                 true, procID ) );
      nodes[NodeID( (*iterCN)->get_id(), (*iterCN)->type() )] = new_node;
    }

  data_->directory.addEntries( nodes );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_Directory::getGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool N_TOP_Directory::getGIDs( const vector<NodeID> & idVec,
                               vector<int> & gidVec )
{
  N_TOP_DirectoryData::NodeDir::DataMap nodes;

  vector<NodeID> ids( idVec );
  data_->directory.getEntries( ids, nodes );

  int size = idVec.size();
  gidVec.resize( size );

  for( int i = 0; i < size; ++i )
    gidVec[i] = nodes[idVec[i]]->GID();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_Directory::getProcs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool N_TOP_Directory::getProcs( const vector<NodeID> & idVec,
                                vector<int> & procVec )
{
  N_TOP_DirectoryData::NodeDir::DataMap nodes;

  vector<NodeID> ids( idVec );
  data_->directory.getEntries( ids, nodes );

  int size = idVec.size();
  procVec.resize( size );

  for( int i = 0; i < size; ++i )
    procVec[i] = nodes[idVec[i]]->proc();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_Directory::getSolnGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool N_TOP_Directory::getSolnGIDs( const vector<NodeID> & idVec,
                                   vector< vector<int> > & gidVec,
                                   vector<int> & procVec )
{
  N_TOP_CktNode * cnp;
  gidVec.resize( idVec.size() );

  getProcs( idVec, procVec );

#ifdef Xyce_PARALLEL_MPI
  typedef Xyce::Parallel::Migrate< NodeID, vector<int> > SGMigrate;
  SGMigrate migrator( *(pdsMgr_->getPDSComm()) );

  vector<int> sortedProcVec( procVec );
  vector<NodeID> ids( idVec );
  SortContainer2( sortedProcVec, ids );

  vector<NodeID> inIDs;
  migrator( sortedProcVec, ids, inIDs );

  SGMigrate::DataMap outGIDs;
  for( unsigned int i = 0; i < inIDs.size(); ++i )
  {
    RCP< vector<int> > gids( rcp( new vector<int>() ) );
    cnp = topMgr_->mainGraphPtr_->FindCktNode( inIDs[i] );
#ifdef HAVE_FLEXIBLE_INSERT
    gids->assign( cnp->get_SolnVarGIDList().begin(),
                    cnp->get_SolnVarGIDList().end() );
#else
    list<int>::const_iterator iterIL = cnp->get_SolnVarGIDList().begin();
    list<int>::const_iterator endIL = cnp->get_SolnVarGIDList().end();
    for( ; iterIL != endIL; ++iterIL )
      gids->push_back( *iterIL );
#endif
    outGIDs[inIDs[i]] = gids;
  }

  SGMigrate::DataMap inGIDs;
  migrator.rvs( sortedProcVec, inIDs, outGIDs, inGIDs );

  for( unsigned int i = 0; i < idVec.size(); ++i )
    gidVec[i] = *(inGIDs[idVec[i]]);

#else

  for( unsigned int i = 0; i < idVec.size(); ++i )
  {
    cnp = topMgr_->mainGraphPtr_->FindCktNode( idVec[i] );
    if( !cnp )
    {
      string msg("Directory node not found: " + idVec[i].first + "\n");
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg );
    }
#ifdef HAVE_FLEXIBLE_INSERT
    gidVec[i].assign( cnp->get_SolnVarGIDList().begin(),
                      cnp->get_SolnVarGIDList().end() );
#else
    gidVec[i].clear();
    list<int>::const_iterator iterIL = cnp->get_SolnVarGIDList().begin();
    list<int>::const_iterator endIL = cnp->get_SolnVarGIDList().end();
    for( ; iterIL != endIL; ++iterIL )
      gidVec[i].push_back( *iterIL );
#endif
  }

#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_Directory::getStateGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool N_TOP_Directory::getStateGIDs( const vector<NodeID> & idVec,
                                    vector< vector<int> > & gidVec,
                                    vector<int> & procVec )
{
  N_TOP_CktNode * cnp;
  gidVec.resize( idVec.size() );

  getProcs( idVec, procVec );

#ifdef Xyce_PARALLEL_MPI
  typedef Xyce::Parallel::Migrate< NodeID, vector<int> > SGMigrate;
  SGMigrate migrator( *(pdsMgr_->getPDSComm()) );

  vector<int> sortedProcVec( procVec );
  vector<NodeID> ids( idVec );
  SortContainer2( sortedProcVec, ids );

  vector<NodeID> inIDs;
  migrator( sortedProcVec, ids, inIDs );

  SGMigrate::DataMap outGIDs;
  for( unsigned int i = 0; i < inIDs.size(); ++i )
  {
    RCP< vector<int> > gids( rcp(new vector<int>()) );
    cnp = topMgr_->mainGraphPtr_->FindCktNode( inIDs[i] );
#ifdef HAVE_FLEXIBLE_INSERT
    gids->assign( cnp->get_StateVarGIDList().begin(),
                    cnp->get_StateVarGIDList().end() );
#else
    list<int>::const_iterator iterIL = cnp->get_StateVarGIDList().begin();
    list<int>::const_iterator endIL = cnp->get_StateVarGIDList().end();
    for( ; iterIL != endIL; ++iterIL )
      gids->push_back( *iterIL );
#endif
    outGIDs[inIDs[i]] = gids;
  }

  SGMigrate::DataMap inGIDs;
  migrator.rvs( sortedProcVec, inIDs, outGIDs, inGIDs );

  for( unsigned int i = 0; i < idVec.size(); ++i )
    gidVec[i] = *(inGIDs[idVec[i]]);

#else

  for( unsigned int i = 0; i < idVec.size(); ++i )
  {
    cnp = topMgr_->mainGraphPtr_->FindCktNode( idVec[i] );
#ifdef HAVE_FLEXIBLE_INSERT
    gidVec[i].assign( cnp->get_StateVarGIDList().begin(),
                      cnp->get_StateVarGIDList().end() );
#else
    gidVec[i].clear();
    list<int>::const_iterator iterIL = cnp->get_StateVarGIDList().begin();
    list<int>::const_iterator endIL = cnp->get_StateVarGIDList().end();
    for( ; iterIL != endIL; ++iterIL )
      gidVec[i].push_back( *iterIL );
#endif
  }

#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_Directory::getStoreGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_TOP_Directory::getStoreGIDs( const vector<NodeID> & idVec,
                                    vector< vector<int> > & gidVec,
                                    vector<int> & procVec )
{
  N_TOP_CktNode * cnp;
  gidVec.resize( idVec.size() );

  getProcs( idVec, procVec );

#ifdef Xyce_PARALLEL_MPI
  typedef Xyce::Parallel::Migrate< NodeID, vector<int> > SGMigrate;
  SGMigrate migrator( *(pdsMgr_->getPDSComm()) );

  vector<int> sortedProcVec( procVec );
  vector<NodeID> ids( idVec );
  SortContainer2( sortedProcVec, ids );

  vector<NodeID> inIDs;
  migrator( sortedProcVec, ids, inIDs );

  SGMigrate::DataMap outGIDs;
  for( unsigned int i = 0; i < inIDs.size(); ++i )
  {
    RCP< vector<int> > gids( rcp(new vector<int>()) );
    cnp = topMgr_->mainGraphPtr_->FindCktNode( inIDs[i] );
#ifdef HAVE_FLEXIBLE_INSERT
    gids->assign( cnp->get_StoreVarGIDList().begin(),
                    cnp->get_StoreVarGIDList().end() );
#else
    list<int>::const_iterator iterIL = cnp->get_StoreVarGIDList().begin();
    list<int>::const_iterator endIL = cnp->get_StoreVarGIDList().end();
    for( ; iterIL != endIL; ++iterIL )
      gids->push_back( *iterIL );
#endif
    outGIDs[inIDs[i]] = gids;
  }

  SGMigrate::DataMap inGIDs;
  migrator.rvs( sortedProcVec, inIDs, outGIDs, inGIDs );

  for( unsigned int i = 0; i < idVec.size(); ++i )
    gidVec[i] = *(inGIDs[idVec[i]]);

#else

  for( unsigned int i = 0; i < idVec.size(); ++i )
  {
    cnp = topMgr_->mainGraphPtr_->FindCktNode( idVec[i] );
#ifdef HAVE_FLEXIBLE_INSERT
    gidVec[i].assign( cnp->get_StoreVarGIDList().begin(),
                      cnp->get_StoreVarGIDList().end() );
#else
    gidVec[i].clear();
    list<int>::const_iterator iterIL = cnp->get_StoreVarGIDList().begin();
    list<int>::const_iterator endIL = cnp->get_StoreVarGIDList().end();
    for( ; iterIL != endIL; ++iterIL )
      gidVec[i].push_back( *iterIL );
#endif
  }

#endif

  return true;
}

