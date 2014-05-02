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
// Filename       : $RCSfile: N_PDS_ParDir.C,v $
//
// Purpose        : Distributed directory for circuit node info
//
// Special Notes  :
//
// Creator        : Rob Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.24 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <N_PDS_ParDir.h>
#include <N_PDS_Comm.h>
#include <N_PDS_CommFactory.h>
#include <N_PDS_Migrator.h>
#include <N_UTL_Functors.h>
#include <N_ERH_ErrorMgr.h>

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::N_PDS_ParDir
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
N_PDS_ParDir::N_PDS_ParDir( N_PDS_Comm * comm, N_PDS_Migrator * migrator )
 : pdsComm_(comm),
   commOwned_(false),
   pdsMigrator_(migrator),
   migratorOwned_(false)
{
  if( !pdsComm_ )
  {
    pdsComm_ = N_PDS_CommFactory::create();
    commOwned_ = true;
  }

  if( !pdsMigrator_ )
  {
    pdsMigrator_ = new N_PDS_Migrator( pdsComm_ );
    migratorOwned_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::registerPDSComm
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::registerPDSComm( N_PDS_Comm * comm )
{
  if( pdsComm_ && commOwned_ ) delete pdsComm_;

  pdsComm_ = comm;

  if( !pdsMigrator_ )
  {
    pdsMigrator_ = new N_PDS_Migrator( pdsComm_ );
    migratorOwned_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::~N_PDS_ParDir
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
N_PDS_ParDir::~N_PDS_ParDir()
{
  if( commOwned_ ) delete pdsComm_;
  if( migratorOwned_ ) delete pdsMigrator_;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::addItems
// Purpose       : Add node items to directory
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::addItems(const std::vector<N_TOP_ParNode *> & nodeVec)
{
  int size = nodeVec.size();

#ifdef Xyce_PARALLEL_MPI

  std::vector<int> assign(size);
  for (int i = 0; i < size; ++i)
    assign[i] = parKey(nodeVec[i]->ID());

  std::vector<Packable *> pNodeVec(size);
  for (int i = 0; i < size; ++i)
    pNodeVec[i] = nodeVec[i];

  std::vector<Packable *> pNewNodeVec;

  pdsMigrator_->migratePackable(assign, pNodeVec, pNewNodeVec);

  std::vector<N_TOP_ParNode *> newNodeVec(pNewNodeVec.size());
  int pSize = pNewNodeVec.size();
  for (int i = 0; i < pSize; ++i)
    newNodeVec[i] = reinterpret_cast<N_TOP_ParNode *> (pNewNodeVec[i]);

  int newSize = newNodeVec.size();
  for( int i = 0; i < newSize; ++i )
    if( !nodeMap_.count(newNodeVec[i]->ID()) )
      nodeMap_[newNodeVec[i]->ID()] = *(newNodeVec[i]);

  for (int i = 0; i < newSize; ++i)
    delete newNodeVec[i];

#else

  for (int i = 0; i < size; ++i)
    nodeMap_[nodeVec[i]->ID()] = *(nodeVec[i]);

#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::deleteItems
// Purpose       : Delete node items to directory
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::deleteItems(const std::vector<NodeID> & idVec)
{
  int size = idVec.size();

#ifdef Xyce_PARALLEL_MPI

  std::vector<int> assign(size);
  for (int i = 0; i < size; ++i)
    assign[i] = parKey(idVec[i].first);

  std::vector<std::string> oldIDVec;
  oldIDVec.reserve(size);
  transform(idVec.begin(), idVec.end(), oldIDVec.begin(),
		FirstOfPair<NodeID,std::string>());
  std::vector<std::string> newIDVec;
  pdsMigrator_->migrateString(assign, oldIDVec, newIDVec);

  int newSize = newIDVec.size();
  for (int i = 0; i < newSize; ++i)
    nodeMap_.erase(newIDVec[i]);

#else

  for (int i = 0; i < size; ++i)
    nodeMap_.erase(idVec[i].first);

#endif
}

#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::addGhosts
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::addGhosts( const std::vector< std::pair<NodeID,int> > & ghostVec )
{
  int size = ghostVec.size();

  std::vector<int> assign(size);
  for( int i = 0; i < size; ++i )
    assign[i] = parKey( ghostVec[i].first.first );

  std::vector< std::pair<std::string,int> > oldGhostVec;
  oldGhostVec.reserve(size);
  for( int i = 0; i < size; ++i )
    oldGhostVec[i] = std::pair<std::string,int>( ghostVec[i].first.first,
                                       ghostVec[i].second );
  std::vector< std::pair<std::string,int> > newGhostVec;
  pdsMigrator_->migrateInt( assign, oldGhostVec, newGhostVec );

  for( int i = 0; i < newGhostVec.size(); ++i )
    nodeMap_[ newGhostVec[i].first ].ghosts().insert( newGhostVec[i].second );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::deleteGhosts
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::deleteGhosts( const std::vector< std::pair<NodeID,int> > & ghostVec )
{
  int size = ghostVec.size();

  std::vector<int> assign(size);
  for( int i = 0; i < size; ++i )
    assign[i] = parKey( ghostVec[i].first.first );

  std::vector< std::pair<std::string,int> > oldGhostVec;
  oldGhostVec.reserve(size);
  for( int i = 0; i < size; ++i )
    oldGhostVec[i] = std::pair<std::string,int>( ghostVec[i].first.first,
                                       ghostVec[i].second );
  std::vector< std::pair<std::string,int> > newGhostVec;
  pdsMigrator_->migrateInt( assign, oldGhostVec, newGhostVec );

  for( int i = 0; i < newGhostVec.size(); ++i )
    nodeMap_[ newGhostVec[i].first ].ghosts().erase( newGhostVec[i].second );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::clearGhosts
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::clearGhosts( const std::vector<NodeID> & idVec )
{
  int size = idVec.size();

  std::vector<int> assign(size);
  for( int i = 0; i < size; ++i )
    assign[i] = parKey( idVec[i].first );

  std::vector<std::string> oldIDVec;
  oldIDVec.reserve(size);
  transform( idVec.begin(), idVec.end(), oldIDVec.begin(),
		FirstOfPair<NodeID,std::string>() );
  std::vector<std::string> newIDVec;
  pdsMigrator_->migrateString( assign, oldIDVec, newIDVec );

  for( int i = 0; i < newIDVec.size(); ++i )
    nodeMap_[ newIDVec[i] ].ghosts().clear();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::clearGhosts
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::clearGhosts()
{
  for( std::map<std::string,N_TOP_ParNode>::iterator it_M = nodeMap_.begin();
       it_M != nodeMap_.end(); ++it_M )
    it_M->second.ghosts().clear();
}

#endif

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::getItems
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::getItems( const std::vector<NodeID> & idVec,
                             std::vector<N_TOP_ParNode *> & nodeVec)
{
  int size = idVec.size();

#ifdef Xyce_PARALLEL_MPI

  std::vector<int> assign(size);
  for (int i = 0; i < size; ++i)
    assign[i] = parKey(idVec[i].first);

  std::vector<std::string> oldIDVec;
  oldIDVec.reserve(size);
  transform(idVec.begin(), idVec.end(), oldIDVec.begin(), FirstOfPair<NodeID,std::string>());
  std::vector<std::string> newIDVec;
  pdsMigrator_->migrateString(assign, oldIDVec, newIDVec);

  int newSize = newIDVec.size();
  std::vector<Packable *> newNodeVec(newSize);
  for (int i = 0; i < newSize; ++i)
  {
    if( !nodeMap_.count( newIDVec[i] ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Node not in directory: " + newIDVec[i] + "\n" );
    newNodeVec[i] = &(nodeMap_[newIDVec[i]]);
  }

  std::vector<Packable *> pNodeVec;
  pdsMigrator_->reverseMigratePackable(assign, newNodeVec, pNodeVec);

  int pSize = pNodeVec.size();
  nodeVec.resize(pSize);
  for (int i = 0; i < pSize; ++i)
    nodeVec[i] = PNP_CAST(pNodeVec[i]);

#else

  nodeVec.resize(size);
  for (int i = 0; i < size; ++i)
  {
    if( !nodeMap_.count( idVec[i].first ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Node not in directory: " + idVec[i].first + "\n" );

    nodeVec[i] = new N_TOP_ParNode(nodeMap_[ idVec[i].first ]);
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::getGIDs
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::getGIDs( const std::vector<NodeID> & idVec,
			    std::vector<int> & gidVec )
{
  int size = idVec.size();

#ifdef Xyce_PARALLEL_MPI

  std::vector<int> assign(size);
  for (int i = 0; i < size; ++i)
    assign[i] = parKey(idVec[i].first);

  std::vector<std::string> oldIDVec;
  oldIDVec.reserve(size);
  transform(idVec.begin(), idVec.end(), oldIDVec.begin(),
		FirstOfPair<NodeID,std::string>());
  std::vector<std::string> newIDVec;
  pdsMigrator_->migrateString(assign, oldIDVec, newIDVec);

  int newSize = newIDVec.size();
  std::vector< std::pair<std::string,int> > newGIDVec(newSize);
  for (int i = 0; i < newSize; ++i)
  {
    if( !nodeMap_.count( newIDVec[i] ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Node not in directory: " + newIDVec[i] + "\n" );

    newGIDVec[i].first  = newIDVec[i];
    newGIDVec[i].second = nodeMap_[newIDVec[i]].GID();
  }

  std::vector< std::pair<std::string,int> > tmpGIDVec;
  pdsMigrator_->reverseMigrateInt(assign, newGIDVec, tmpGIDVec);

  std::map< std::string, int > tmpMap;
  int tmpSize = tmpGIDVec.size();
  for (int i = 0; i < tmpSize; ++i)
    tmpMap[tmpGIDVec[i].first] = tmpGIDVec[i].second;

  int idSize = idVec.size();
  for (int i = 0; i < idSize; ++i)
    gidVec[i] = tmpMap[idVec[i].first];

#else

  gidVec.resize(size);
  for (int i = 0; i < size; ++i)
  {
    if( !nodeMap_.count( idVec[i].first ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Node not in directory: " + idVec[i].first + "\n" );

    gidVec[i] = nodeMap_[idVec[i].first].GID();
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::getProcs
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::getProcs( const std::vector<NodeID> & idVec,
			     std::vector<int> & procVec )
{
  int size = idVec.size();

#ifdef Xyce_PARALLEL_MPI

  std::vector<int> assign(size);
  for (int i = 0; i < size; ++i)
    assign[i] = parKey(idVec[i].first);

  std::vector<std::string> oldIDVec;
  oldIDVec.resize(size);
  transform(idVec.begin(), idVec.end(), oldIDVec.begin(),
		FirstOfPair<NodeID,std::string>());

  std::vector<std::string> newIDVec;
  pdsMigrator_->migrateString(assign, oldIDVec, newIDVec);

  int newSize = newIDVec.size();
  std::vector< std::pair<std::string,int> > newProcVec(newSize);
  for (int i = 0; i < newSize; ++i )
  {
    if( !nodeMap_.count( newIDVec[i] ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Node not in directory: " + newIDVec[i] + "\n");

    newProcVec[i].first  = newIDVec[i];
    newProcVec[i].second = nodeMap_[newIDVec[i]].proc();
  }

  std::vector< std::pair<std::string,int> > tmpProcVec;
  pdsMigrator_->reverseMigrateInt(assign, newProcVec, tmpProcVec);

  std::map< std::string, int > tmpMap;
  int tmpSize = tmpProcVec.size();
  for (int i = 0; i < tmpSize; ++i)
    tmpMap[tmpProcVec[i].first] = tmpProcVec[i].second;

  int idSize = idVec.size();
  procVec.resize(idSize);
  for (int i = 0; i < idSize; ++i)
    procVec[i] = tmpMap[idVec[i].first];

#else

  procVec.resize(size);
  for (int i = 0; i < size; ++i)
  {
    if( !nodeMap_.count( idVec[i].first ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Node not in directory: " + idVec[i].first + "\n" );

    procVec[i] = nodeMap_[idVec[i].first].proc();
  }

#endif
}

#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::getGhosts
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_ParDir::getGhosts( const std::vector<NodeID> & idVec,
			      std::vector< std::vector<int> > & ghostVec )
{
  int size = idVec.size();

  std::vector<int> assign(size);
  for( int i = 0; i < size; ++i )
    assign[i] = parKey( idVec[i].first );

  std::vector<std::string> oldIDVec;
  oldIDVec.reserve(size);
  transform( idVec.begin(), idVec.end(), oldIDVec.begin(),
		FirstOfPair<NodeID,std::string>() );
  std::vector<std::string> newIDVec;
  pdsMigrator_->migrateString( assign, oldIDVec, newIDVec );

  std::vector< std::pair< std::string,std::vector<int> > > newGhostVec( newIDVec.size() );
  for( int i = 0; i < newIDVec.size(); ++i )
  {
    if( !nodeMap_.count( newIDVec[i] ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Node not in directory: " + newIDVec[i] + "\n" );

    newGhostVec[i].first = newIDVec[i];
    newGhostVec[i].second.resize( nodeMap_[ newIDVec[i] ].ghosts().size() );
    std::set<int>::iterator it_M = nodeMap_[ newIDVec[i] ].ghosts().begin();
    for( int j = 0; j < newGhostVec[i].second.size(); ++j, ++it_M )
      newGhostVec[i].second[j] = *it_M;
  }

  std::vector< std::pair< std::string,std::vector<int> > > tmpGhostVec;
  pdsMigrator_->reverseMigrateIntVec( assign, newGhostVec, tmpGhostVec );

  std::map< std::string, std::vector<int> > tmpMap;
  for( int i = 0; i < tmpGhostVec.size(); ++i )
    tmpMap[ tmpGhostVec[i].first ] = tmpGhostVec[i].second;

  for( int i = 0; i < idVec.size(); ++i )
    ghostVec[i] = tmpMap[ idVec[i].first ];

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::parKey
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
int N_PDS_ParDir::parKey(const std::string & token)
{
  int sum = 0;

  for( int i = 0; i < token.length(); ++i )
    sum += token[i];

  float fsum = sum;
  float fproc = pdsComm_->numProc();

  return fmod(fsum, fproc);
}

#endif

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::debugDump
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
bool N_PDS_ParDir::debugDump( std::ostream & os ) const
{
#ifdef Xyce_PARALLEL_MPI
  pdsComm_->barrier();
#endif

  int numProc = pdsComm_->numProc();
  int myProc = pdsComm_->procID();

  for( int i = 0; i < numProc; ++i )
  {
    if( i == myProc )
    {
      os << "Directory for Proc: " << myProc << std::endl;
      os << Xyce::section_divider << std::endl;
      for( std::map<std::string,N_TOP_ParNode>::const_iterator iterSP =
	    nodeMap_.begin(); iterSP != nodeMap_.end(); ++iterSP )
        os << iterSP->second;
      os << Xyce::section_divider << std::endl;
    }

#ifdef Xyce_PARALLEL_MPI
    pdsComm_->barrier();
#endif
  }

  return true;
}

