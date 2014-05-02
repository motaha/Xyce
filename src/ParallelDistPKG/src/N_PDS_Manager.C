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
// Filename       : $RCSfile: N_PDS_Manager.C,v $
//
// Purpose        : Manager for the parallel load-balance, distribution and
//                  migration tools.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.44 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_CommFactory.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_ParMapFactory.h>
#include <N_PDS_GlobalAccessor.h>
#include <N_PDS_Migrator.h>
#include <N_PDS_ParDir.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_LAFactory.h>
#include <Epetra_CrsGraph.h>
#undef HAVE_LIBPARMETIS
#include <EpetraExt_View_CrsGraph.h>

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::N_PDS_Manager
// Purpose       : Default constructor.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/09/00
//-----------------------------------------------------------------------------
N_PDS_Manager::N_PDS_Manager(bool & isSerial, bool & procFlag, int iargs, char **cargs, Xyce::Parallel::Machine comm)
  : Comm_(0),
    Topo_(0)
{
  Comm_ = N_PDS_CommFactory::create(iargs, cargs, comm);

  N_ERH_ErrorMgr::registerComm(Comm_->comm());
  
  isSerial = Comm_->isSerial();
  procFlag = (Comm_->procID() == 0);
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::~N_PDS_Manager
// Purpose       : Destructor - closes up the load balance list and shuts down
//                 MPI.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/09/00
//-----------------------------------------------------------------------------
N_PDS_Manager::~N_PDS_Manager()
{
  std::map<std::string,N_PDS_ParMap *>::iterator it_spmM = pm_Map_.begin();
  std::map<std::string,N_PDS_ParMap *>::iterator end_spmM = pm_Map_.end();
  for( ; it_spmM != end_spmM; ++it_spmM )
  {
    deleteParallelMap( it_spmM->first );
  }

  std::map<std::string,Epetra_CrsGraph*>::iterator it_mgM = mg_Map_.begin();
  std::map<std::string,Epetra_CrsGraph*>::iterator end_mgM = mg_Map_.end();
  for( ; it_mgM != end_mgM; ++it_mgM )
  {
    deleteMatrixGraph( it_mgM->first );
  }

  if( Comm_ ) delete Comm_;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::registerTopology
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/09/00
//-----------------------------------------------------------------------------
bool N_PDS_Manager::registerTopology( N_TOP_Topology * topo )
{
  if( topo )
  {
    Topo_ = topo;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::createParallelMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/31/01
//-----------------------------------------------------------------------------
N_PDS_ParMap * N_PDS_Manager::createParallelMap( int & num_global,
			                         const int & num_local,
                                                 const std::vector<int> & gid_map,
                                                 const int index_base )
{
  return N_PDS_ParMapFactory::create( num_global, num_local, gid_map, index_base, Comm_ );
}

//-----------------------------------------------------------------------------
// Method        : N_PDS_Manager::addParallelMap
// Purpose       : add N_PDS_ParMap obj to container, key{string}
// Scope         : public
// Special Notes :
// Creator       : Robert J Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 01/31/01
//-----------------------------------------------------------------------------
bool N_PDS_Manager::addParallelMap( const std::string & name, N_PDS_ParMap * map )
{
  if( !pm_Map_.count( name ) )
  {
    pm_Map_[ name ] = map;
    return true;
  }
  else
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
      "Parallel Map already exists!\n" );
    return false;
  }
}

//-----------------------------------------------------------------------------
// Method        : N_PDS_Manager::getParallelMap
// Purpose       : get N_PDS_ParMap obj from container, key{string}
// Scope         : public
// Special Notes :
// Creator       : Robert J Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 01/31/01
//-----------------------------------------------------------------------------
N_PDS_ParMap * N_PDS_Manager::getParallelMap( const std::string & name )
{
  if( pm_Map_.count( name ) )
    return pm_Map_[ name ];
  else
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
      "Parallel Map not found!\n" );
    return 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::deleteParallelMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/31/01
//-----------------------------------------------------------------------------
bool N_PDS_Manager::deleteParallelMap( const std::string & name )
{
  if( pm_Map_.count( name ) )
  {
    deleteGlobalAccessor( name );

    if( pm_Map_[name] )
    {
      delete pm_Map_[name];
      pm_Map_[name] = 0;
      return true;
    }
    else
    {
      Xyce::Report::DevelFatal().in("N_PDS_Manager::deleteParallelMap") << "Parallel Map is NULL!";
      return false;
    }
  }
  else
  {
    Xyce::Report::UserWarning() << "Deleting Non-existent Parallel Map!";
    
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::addGlobalAccessor
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/31/01
//-----------------------------------------------------------------------------
bool N_PDS_Manager::addGlobalAccessor( const std::string & name )
{
  if( pm_Map_.count( name ) )
  {
    if( !ga_Map_.count( name ) )
    {
      ga_Map_[ name ] = new N_PDS_GlobalAccessor( pm_Map_[name]->pdsComm() );
      return true;
    }
    else
    {
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
        "Global Accessor already Exists\n" );
      return false;
    }
  }
  else
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
      "Creating Global Accessor without Parallel Map!\n" );
    return false;
  }
}

//-----------------------------------------------------------------------------
// Method        : N_PDS_Manager::getGlobalAccessor
// Purpose       : get N_PDS_GlobalAccessor obj from container, key{string}
// Scope         : public
// Special Notes :
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/31/01
//-----------------------------------------------------------------------------
N_PDS_GlobalAccessor * N_PDS_Manager::getGlobalAccessor( const std::string & name )
{
  if( ga_Map_.count(name) )
    return ga_Map_[ name ];
  else
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
      "Global Accessor not found!\n" );
    return 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::deleteGlobalAccessor
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/31/01
//-----------------------------------------------------------------------------
bool N_PDS_Manager::deleteGlobalAccessor( const std::string & name )
{
  if( ga_Map_.count( name ) )
  {
    if( ga_Map_[name] )
    {
      delete ga_Map_[name];
      ga_Map_.erase( name );
      return true;
    }
    else
    {
      ga_Map_.erase( name );
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
        "Deleting NULL Global Accessor!\n" );
      return false;
    }
  }
  else
    return false;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::createGlobalAccessor
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/1/01
//-----------------------------------------------------------------------------
N_PDS_GlobalAccessor * N_PDS_Manager::createGlobalAccessor
				( const std::string & name )
{
  if( name == "" )
    return new N_PDS_GlobalAccessor( Comm_ );
  else
  {
    addGlobalAccessor( name );
    return getGlobalAccessor( name );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::createMigrator
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
N_PDS_Migrator * N_PDS_Manager::createMigrator() const
{
  return new N_PDS_Migrator( Comm_ );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::createParDir
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
N_PDS_ParDir * N_PDS_Manager::createParDir() const
{
  return new N_PDS_ParDir( Comm_ );
}

//-----------------------------------------------------------------------------
// Method        : N_PDS_Manager::addMatrixGraph
// Purpose       :
// Scope         : public
// Special Notes :
// Creator       : Robert J Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 08/23/02
//-----------------------------------------------------------------------------
bool N_PDS_Manager::addMatrixGraph( const std::string & name,
                                    Epetra_CrsGraph * graph,
                                    EpetraExt::CrsGraph_View * trans )
{
  if( !mg_Map_.count( name ) )
  {
    mg_Map_[ name ] = graph;

    if( trans ) // store owning transform as well
      mgvt_Map_[ name ] = trans;

    return true;
  }
  else
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Matrix Graph already exists!\n" );
    return false;
  }
}

//-----------------------------------------------------------------------------
// Method        : N_PDS_Manager::getMatrixGraph
// Purpose       :
// Scope         : public
// Special Notes :
// Creator       : Robert J Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 08/23/02
//-----------------------------------------------------------------------------
Epetra_CrsGraph * N_PDS_Manager::getMatrixGraph( const std::string & name )
{
  if( mg_Map_.count( name ) )
    return mg_Map_[ name ];
  else
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Matrix Graph not found!\n" );
    return 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Manager::deleteMatrixGraph
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/31/03
//-----------------------------------------------------------------------------
bool N_PDS_Manager::deleteMatrixGraph( const std::string & name )
{
  if( mg_Map_.count( name ) )
  {
    if( mgvt_Map_.count( name ) )
    {
      assert( mgvt_Map_[name]!=0 );
      delete mgvt_Map_[name];
      mgvt_Map_[name] = 0;
      mg_Map_[name] = 0;
    }
    else
    {
      assert( mg_Map_[name]!= 0 );
      delete mg_Map_[name];
      mg_Map_[name] = 0;
    }
    return true;
  }
  else
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
      "Deleting Non-existent Parallel Map!\n" );
    return false;
  }
}

