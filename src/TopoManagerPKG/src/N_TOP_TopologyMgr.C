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
// Filename       : $RCSfile: N_TOP_TopologyMgr.C,v $
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
// Revision Number: $Revision: 1.11.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#ifdef HAVE_CASSERT
 #include <cassert>
#else
 #include <assert.h>
#endif

// ---------- Xyce Includes --------------

#include <N_TOP_TopologyMgr.h>

#include <N_TOP_Topology.h>
#include <N_TOP_DevInsertionTool.h>
#include <N_IO_CmdParse.h>
#include <N_IO_PkgOptionsMgr.h>

namespace Xyce {
namespace Topology {

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::~Manager
// Purpose       : Destructor
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
Manager::~Manager()
{
  if( topo_ ) delete topo_;
  if( currentDevInsertionTool != NULL )
  {
    delete currentDevInsertionTool;
    currentDevInsertionTool = NULL;
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::createSystem
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
/*
System & Manager::createSystem( const string & type, const string & name )
{
  //Check that a System of this 'name' has not been already created
  assert( systems_.count( name ) == 0 );

  System * newSystem = 0;
    
  if( type == "CIRCUIT" )
     newSys = dynamic_cast<System*>( new Circuit::System( name, pdsMgr_ ) );
  else
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
     "Xyce::Topology::Manager::instantiateSystem(...) : System Type Not Recognized>> " + type );

  systems_[name] = newSys;

  return *newSys;
}
*/

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::getSystem
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
/*
System & Manager::getSystem( const string & name )
{
  assert( systems_.count( name ) != 0 );

  return systems_[name];
}
*/

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::deleteSystem
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
/*
bool Manager::deleteSystem( const string & name )
{
  //Check that the System to be deleted is in the Container 'systems_'
  assert( systems_.count( name ) != 0 );

  delete systems_[name];
  systems_.erase( name );
}
*/
    
//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::getPartitionTool
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
/*
Parallel::PartitionTool * Manager::getPartitionTool( const string & sysName,
                                                     const string & typeName )
{
  return PartitionToolFactory::create( getSystem( sysName ), typeName );
}
*/

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::createTopology
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
N_TOP_Topology * Manager::createTopology(N_IO_CmdParse & cp)
{
  return (topo_ = new N_TOP_Topology(cp));
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::getInsertionTool
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
InsertionTool * Manager::getInsertionTool( const string & sysName )
{
//  return InsertionToolFactory::create( getSystem( sysName ) );
  if( currentDevInsertionTool == NULL )
  {
    currentDevInsertionTool = new DevInsertionTool( *topo_ );
  }
  return currentDevInsertionTool;
}

//-----------------------------------------------------------------------------
// Function      :  Xyce::Topology::Manager::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool  Xyce::Topology::Manager::registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  topo_->registerPkgOptionsMgr(pkgOptMgrPtr_);
  return true;
}

/*
//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::getRestartDataTool
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
InputOutput::RestartDataTool * Manager::getRestartDataTool( const string & sysName )
{
  return RestartDataToolFactory::create( getSystem( sysName ) );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::getNodeTool
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
NodeTool * Manager::getNodeTool( const string & sysName )
{
  return new NodeTool( getSystem( sysName ) );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::Manager::getLASQueryUtil
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
N_LAS_QueryUtil * Manager::getLASQueryUtil( const string & sysName )
{
  return getSystem( sysName ).getLASQueryUtil();
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes : friend
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
ostream & operator<<( ostream & os, const Manager & mgr );
{
  return os;
}
*/

} //namespace Topology
} //namespace Xyce
