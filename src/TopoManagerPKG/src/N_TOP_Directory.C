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
// Revision Number: $Revision: 1.28 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_TOP_Misc.h>

#include <iostream>
#include <algorithm>

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

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : DirectoryData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/26/04
//-----------------------------------------------------------------------------
struct DirectoryData
{
  typedef RCP<ParNode> NodePtr;
  typedef std::map<NodeID,NodePtr> NodeContainer;
                                                                                           
  typedef Xyce::Parallel::Hash<NodeID> NodeIDHash;
  typedef Xyce::Parallel::Migrate<NodeID,ParNode> NodeMigrate;

  typedef Xyce::Parallel::Directory< NodeID,
                                     ParNode,
                                     NodeIDHash,
                                     NodeContainer,
                                     NodeMigrate >  NodeDir;

  DirectoryData( N_PDS_Comm & comm )
  : hash(comm.numProc()),
    migrate(comm),
    directory(migrate,hash)
  {}

  ~DirectoryData() {}

  NodeIDHash hash;
  NodeMigrate migrate;

  NodeDir directory;
};
  

//-----------------------------------------------------------------------------
// Function      : Directory::~Directory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/02/01
//-----------------------------------------------------------------------------
Directory::~Directory()
{
  if( data_ ) delete data_;
}

//-----------------------------------------------------------------------------
// Function      : Directory::generateDirectory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/02/01
//-----------------------------------------------------------------------------
bool Directory::generateDirectory()
{
  int procID = pdsMgr_->getPDSComm()->procID();

  if( data_ ) delete data_;
  data_ = new DirectoryData( *(pdsMgr_->getPDSComm()) );

  DirectoryData::NodeContainer nodes;

  std::list<CktNode*>::iterator iterCN = topMgr_->orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator endCN = topMgr_->orderedNodeListPtr_->end();
  
  for( ; iterCN != endCN; ++iterCN )
    if( (*iterCN)->get_IsOwned() && (*iterCN)->get_gID() != -1 )
    {
      DirectoryData::NodePtr new_node(
              new ParNode( NodeID( (*iterCN)->get_id(), (*iterCN)->get_gID() ),
                                 true, procID ) );
      nodes[NodeID( (*iterCN)->get_id(), (*iterCN)->type() )] = new_node;
    }

  data_->directory.addEntries( nodes );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Directory::getGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool Directory::getGIDs( const std::vector<NodeID> & idVec,
                               std::vector<int> & gidVec )
{
  DirectoryData::NodeDir::DataMap nodes;

  std::vector<NodeID> ids( idVec );
  data_->directory.getEntries( ids, nodes );

  int size = idVec.size();
  gidVec.resize( size );

  for( int i = 0; i < size; ++i )
    gidVec[i] = nodes[idVec[i]]->GID();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Directory::getProcs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool Directory::getProcs( const std::vector<NodeID> & idVec,
                                std::vector<int> & procVec )
{
  DirectoryData::NodeDir::DataMap nodes;

  std::vector<NodeID> ids( idVec );
  data_->directory.getEntries( ids, nodes );

  int size = idVec.size();
  procVec.resize( size );

  for( int i = 0; i < size; ++i )
    procVec[i] = nodes[idVec[i]]->proc();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Directory::getSolnGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool Directory::getSolnGIDs( const std::vector<NodeID> & idVec,
                                   std::vector< std::vector<int> > & gidVec,
                                   std::vector<int> & procVec )
{
  CktNode * cnp;
  gidVec.resize( idVec.size() );

  getProcs( idVec, procVec );

#ifdef Xyce_PARALLEL_MPI
  typedef Xyce::Parallel::Migrate< NodeID, std::vector<int> > SGMigrate;
  SGMigrate migrator( *(pdsMgr_->getPDSComm()) );

  std::vector<int> sortedProcVec( procVec );
  std::vector<NodeID> ids( idVec );
  SortContainer2( sortedProcVec, ids );

  std::vector<NodeID> inIDs;
  migrator( sortedProcVec, ids, inIDs );

  SGMigrate::DataMap outGIDs;
  for( unsigned int i = 0; i < inIDs.size(); ++i )
  {
    RCP< std::vector<int> > gids( rcp( new std::vector<int>() ) );
    cnp = topMgr_->mainGraphPtr_->FindCktNode( inIDs[i] );
    gids->assign( cnp->get_SolnVarGIDList().begin(),
                    cnp->get_SolnVarGIDList().end() );
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
      std::string msg("Directory node not found: " + idVec[i].first + "\n");
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, msg );
    }
    gidVec[i].assign( cnp->get_SolnVarGIDList().begin(),
                      cnp->get_SolnVarGIDList().end() );
  }

#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Directory::getStateGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/19/01
//-----------------------------------------------------------------------------
bool Directory::getStateGIDs( const std::vector<NodeID> & idVec,
                                    std::vector< std::vector<int> > & gidVec,
                                    std::vector<int> & procVec )
{
  CktNode * cnp;
  gidVec.resize( idVec.size() );

  getProcs( idVec, procVec );

#ifdef Xyce_PARALLEL_MPI
  typedef Xyce::Parallel::Migrate< NodeID, std::vector<int> > SGMigrate;
  SGMigrate migrator( *(pdsMgr_->getPDSComm()) );

  std::vector<int> sortedProcVec( procVec );
  std::vector<NodeID> ids( idVec );
  SortContainer2( sortedProcVec, ids );

  std::vector<NodeID> inIDs;
  migrator( sortedProcVec, ids, inIDs );

  SGMigrate::DataMap outGIDs;
  for( unsigned int i = 0; i < inIDs.size(); ++i )
  {
    RCP< std::vector<int> > gids( rcp(new std::vector<int>()) );
    cnp = topMgr_->mainGraphPtr_->FindCktNode( inIDs[i] );
    gids->assign( cnp->get_StateVarGIDList().begin(),
                    cnp->get_StateVarGIDList().end() );
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
    gidVec[i].assign( cnp->get_StateVarGIDList().begin(),
                      cnp->get_StateVarGIDList().end() );
  }

#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Directory::getStoreGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Directory::getStoreGIDs( const std::vector<NodeID> & idVec,
                                    std::vector< std::vector<int> > & gidVec,
                                    std::vector<int> & procVec )
{
  CktNode * cnp;
  gidVec.resize( idVec.size() );

  getProcs( idVec, procVec );

#ifdef Xyce_PARALLEL_MPI
  typedef Xyce::Parallel::Migrate< NodeID, std::vector<int> > SGMigrate;
  SGMigrate migrator( *(pdsMgr_->getPDSComm()) );

  std::vector<int> sortedProcVec( procVec );
  std::vector<NodeID> ids( idVec );
  SortContainer2( sortedProcVec, ids );

  std::vector<NodeID> inIDs;
  migrator( sortedProcVec, ids, inIDs );

  SGMigrate::DataMap outGIDs;
  for( unsigned int i = 0; i < inIDs.size(); ++i )
  {
    RCP< std::vector<int> > gids( rcp(new std::vector<int>()) );
    cnp = topMgr_->mainGraphPtr_->FindCktNode( inIDs[i] );
    gids->assign( cnp->get_StoreVarGIDList().begin(),
                    cnp->get_StoreVarGIDList().end() );
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
    gidVec[i].assign( cnp->get_StoreVarGIDList().begin(),
                      cnp->get_StoreVarGIDList().end() );
  }

#endif

  return true;
}

} // namespace Topo
} // namespace Xyce
