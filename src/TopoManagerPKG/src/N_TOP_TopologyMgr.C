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
// Revision Number: $Revision: 1.18 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#ifdef HAVE_CASSERT
 #include <cassert>
#else
 #include <assert.h>
#endif

#include <N_TOP_TopologyMgr.h>

#include <N_TOP_Topology.h>
#include <N_TOP_DevInsertionTool.h>
#include <N_IO_CmdParse.h>
#include <N_IO_PkgOptionsMgr.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::~Manager
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
// Function      : Xyce::Topo::Manager::createSystem
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
/*
System & Manager::createSystem( const std::string & type, const std::string & name )
{
  //Check that a System of this 'name' has not been already created
  assert( systems_.count( name ) == 0 );

  System * newSystem = 0;
    
  if( type == "CIRCUIT" )
     newSys = dynamic_cast<System*>( new Circuit::System( name, pdsMgr_ ) );
  else
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
     "Xyce::Topo::Manager::instantiateSystem(...) : System Type Not Recognized>> " + type );

  systems_[name] = newSys;

  return *newSys;
}
*/

//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::getSystem
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
/*
System & Manager::getSystem( const std::string & name )
{
  assert( systems_.count( name ) != 0 );

  return systems_[name];
}
*/

//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::deleteSystem
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
/*
bool Manager::deleteSystem( const std::string & name )
{
  //Check that the System to be deleted is in the Container 'systems_'
  assert( systems_.count( name ) != 0 );

  delete systems_[name];
  systems_.erase( name );
}
*/
    
//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::getPartitionTool
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
/*
Parallel::PartitionTool * Manager::getPartitionTool( const std::string & sysName,
                                                     const std::string & typeName )
{
  return PartitionToolFactory::create( getSystem( sysName ), typeName );
}
*/

//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::createTopo
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
N_TOP_Topology * Manager::createTopology(IO::CmdParse & cp)
{
  return (topo_ = new N_TOP_Topology(cp));
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::getInsertionTool
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
InsertionTool * Manager::getInsertionTool( const std::string & sysName )
{
//  return InsertionToolFactory::create( getSystem( sysName ) );
  if( currentDevInsertionTool == NULL )
  {
    currentDevInsertionTool = new DevInsertionTool( *topo_ );
  }
  return currentDevInsertionTool;
}

//-----------------------------------------------------------------------------
// Function      :  Xyce::Topo::Manager::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool  Xyce::Topo::Manager::registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  topo_->registerPkgOptionsMgr(pkgOptMgrPtr_);
  return true;
}

/*
//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::getRestartDataTool
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
InputOutput::RestartDataTool * Manager::getRestartDataTool( const std::string & sysName )
{
  return RestartDataToolFactory::create( getSystem( sysName ) );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::getNodeTool
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
NodeTool * Manager::getNodeTool( const std::string & sysName )
{
  return new NodeTool( getSystem( sysName ) );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Topo::Manager::getLASQueryUtil
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/24/03
//-----------------------------------------------------------------------------
N_LAS_QueryUtil * Manager::getLASQueryUtil( const std::string & sysName )
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

} //namespace Topo
} //namespace Xyce
