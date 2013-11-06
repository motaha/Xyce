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
// Filename       : $RCSfile: N_TOP_CktNode_Dev.C,v $
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
// Revision Number: $Revision: 1.47.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ---------
#include <N_UTL_Misc.h>
#include <N_TOP_Misc.h>

#ifdef HAVE_CASSERT
 #include <cassert>
#else
 #include <assert.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_TOP_CktNode_Dev.h>

#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInterface.h>

#include <N_UTL_Interface_Enum_Types.h>
#include <N_UTL_Expression.h>

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::~N_TOP_CktNode_Dev
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/22/03
//-----------------------------------------------------------------------------
N_TOP_CktNode_Dev::~N_TOP_CktNode_Dev()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::instantiate
// Purpose       : instantiate device with dev pkg
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/22/03
//-----------------------------------------------------------------------------
bool N_TOP_CktNode_Dev::instantiate()
{
  if( devPtr_ )
    return false;
  else
  {
    instance_.assert_not_null();
    assert( iface_ != 0 );
    devPtr_ = iface_->addDeviceInstance( *instance_ );
    assert( devPtr_ != 0 );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::getDevState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
N_DEV_DeviceState * N_TOP_CktNode_Dev::getDevState()
{
  return devPtr_->getInternalState();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::setDevState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
bool N_TOP_CktNode_Dev::setDevState( const N_DEV_DeviceState & state )
{
  return devPtr_->setInternalState( state );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::solnVarCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/26/01
//-----------------------------------------------------------------------------
int N_TOP_CktNode_Dev::solnVarCount()
{
  return devPtr_->getNumIntVars();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::stateVarCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/26/01
//-----------------------------------------------------------------------------
int N_TOP_CktNode_Dev::stateVarCount()
{
  return devPtr_->getNumStateVars();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::storeVarCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
int N_TOP_CktNode_Dev::storeVarCount()
{
  return devPtr_->getNumStoreVars();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::leadConnect
// Purpose       : Find which leads have connections to each other
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 6/20/05
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::leadConnect(vector<int> & connect)
{
  devPtr_->getDevConMap(connect);
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerGIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/22/00
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerGIDswithDev(
		const list<index_pair> & intGIDList,
		const list<index_pair> & extGIDList )
{
  devPtr_->registerGIDs( intGIDList, extGIDList );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerStateGIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerStateGIDswithDev(
		const list<index_pair> & stateGIDList )
{
#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "reg state gids: " << id_ << endl;

  list<index_pair>::const_iterator it_ipL = stateGIDList.begin(); 
  list<index_pair>::const_iterator it_ipL_end = stateGIDList.end();
  for( ; it_ipL != it_ipL_end; ++it_ipL )
  {
    cout << "   " << it_ipL->row << " " << it_ipL->col;
  }

  cout << endl;
#endif

  devPtr_->registerStateGIDs( stateGIDList );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerStoreGIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerStoreGIDswithDev(
		const list<index_pair> & storeGIDList )
{
#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "reg store gids: " << id_ << endl;

  list<index_pair>::const_iterator it_ipL = storeGIDList.begin(); 
  list<index_pair>::const_iterator it_ipL_end = storeGIDList.end();
  for( ; it_ipL != it_ipL_end; ++it_ipL )
  {
    cout << "   " << it_ipL->row << " " << it_ipL->col;
  }

  cout << endl;
#endif

  devPtr_->registerStoreGIDs( storeGIDList );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/02
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerLIDswithDev( const vector<int> & intLIDVec,
                                             const vector<int> & extLIDVec )
{
  devPtr_->registerLIDs( intLIDVec, extLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerStateLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/02
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerStateLIDswithDev( const vector<int> & stateLIDVec )
{
  devPtr_->registerStateLIDs( stateLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerStoreLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerStoreLIDswithDev( const vector<int> & storeLIDVec )
{
  devPtr_->registerStoreLIDs( storeLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerDepLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerDepLIDswithDev( const vector< vector<int> > & depLIDVec )
{
  devPtr_->registerDepSolnLIDs( depLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerDepStateLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerDepStateLIDswithDev( const vector< vector<int> > & depStateLIDVec )
{
  devPtr_->registerDepStateLIDs( depStateLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerDepStoreLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerDepStoreLIDswithDev( const vector< vector<int> > & depStoreLIDVec )
{
  devPtr_->registerDepStoreLIDs( depStoreLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/05/01
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::getDepSolnVars( vector< NodeID >& dsVars )
{
  // Get dependent solution variables from the device and the sDepend structure to determine
  // which node type to look for in the graph.
  dsVars.clear();
  std::vector< std::string > tmpSolnVars = devPtr_->getDepSolnVars();

  std::vector< std::string >::iterator it_sV = tmpSolnVars.begin();
  std::vector< std::string >::iterator it_sV_end = tmpSolnVars.end();

  int type = 0;

  for( ; it_sV != it_sV_end; ++it_sV )
  {
    bool found = false;
    std::vector< sDepend > depParams = devPtr_->getDependentParams();
    std::vector< sDepend >::iterator it_sdV = depParams.begin();
    std::vector< sDepend >::iterator it_sdV_end = depParams.end();

    for( ; it_sdV != it_sdV_end; ++it_sdV )
    {
      type = (it_sdV->expr)->get_type( *it_sV );
      if (type == XEXP_NODE)
      {
        dsVars.push_back( NodeID( *it_sV, _VNODE ) );
        found = true;
        break;
      }
      if (type == XEXP_INSTANCE || type == XEXP_LEAD)
      {
        dsVars.push_back( NodeID( *it_sV, _DNODE ) );
        found = true;
        break;
      }
    }
    // Punt if the node has not been found, hope that the graph only 
    // has one node associated with this name.
    if (!found)
    {
      dsVars.push_back( NodeID( *it_sV, -1 ) );
    }
  } 
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerDepSolnGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerDepSolnGIDs(
		const vector< vector<int> > & dsGIDs )
{

#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "reg secondary gids: " << id_ << endl;

  for( unsigned int i = 0; i < dsGIDs.size(); ++i )
  {
    for( unsigned int j = 0; j < dsGIDs[i].size(); ++j )
      cout << " " << dsGIDs[i][j];
    cout << endl;
  }

  cout << endl;
#endif

  devPtr_->registerDepSolnGIDs( dsGIDs );

  devPtr_->getDepSolnGIDVec( depSolnGIDJacVec_ );

}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::getDepStateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/05/01
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::getDepStateVars( vector< NodeID >& dsVars )
{
  // Clear input vector before loading it with dependent state variables for this node.
  dsVars.clear();
  std::vector< std::string > tmpStateVars = devPtr_->getDepStateVars();

  std::vector< std::string >::iterator it_sV = tmpStateVars.begin();
  std::vector< std::string >::iterator it_sV_end = tmpStateVars.end();

  for( ; it_sV != it_sV_end; ++it_sV )
  {
    bool found = false;
    std::vector< sDepend > depParams = devPtr_->getDependentParams();
    std::vector< sDepend >::iterator it_sdV = depParams.begin();
    std::vector< sDepend >::iterator it_sdV_end = depParams.end();

    for( ; it_sdV != it_sdV_end; ++it_sdV )
    {
      int type = (it_sdV->expr)->get_type( *it_sV );
      if (type == XEXP_NODE)
      {
        dsVars.push_back( NodeID( *it_sV, _VNODE ) );
        found = true;
        break;
      }
      if (type == XEXP_INSTANCE)
      {
        dsVars.push_back( NodeID( *it_sV, _DNODE ) );
        found = true;
        break;
      }
    }
    // Punt if the node has not been found, hope that the graph only 
    // has one node associated with this name.
    if (!found)
    {
      dsVars.push_back( NodeID( *it_sV, -1 ) );
    }
  } 
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerDepStateGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerDepStateGIDs(
		const vector< vector<int> > & dsGIDs )
{
#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "reg secondary gids: " << id_ << endl;

  for( unsigned int i = 0; i < dsGIDs.size(); ++i )
  {
    for( unsigned int j = 0; j < dsGIDs[i].size(); ++j )
      cout << " " << dsGIDs[i][j];
    cout << endl;
  }

  cout << endl;
#endif

  devPtr_->registerDepStateGIDs( dsGIDs );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::getDepStoreVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::getDepStoreVars( vector< NodeID >& dsVars )
{
  // Clear input vector before loading it with dependent store variables for this node.
  dsVars.clear();
  std::vector< std::string > tmpStoreVars = devPtr_->getDepStoreVars();

  std::vector< std::string >::iterator it_sV = tmpStoreVars.begin();
  std::vector< std::string >::iterator it_sV_end = tmpStoreVars.end();

  for( ; it_sV != it_sV_end; ++it_sV )
  {
    bool found = false;
    std::vector< sDepend > depParams = devPtr_->getDependentParams();
    std::vector< sDepend >::iterator it_sdV = depParams.begin();
    std::vector< sDepend >::iterator it_sdV_end = depParams.end();

    for( ; it_sdV != it_sdV_end; ++it_sdV )
    {
      int type = (it_sdV->expr)->get_type( *it_sV );
      if (type == XEXP_NODE)
      {
        dsVars.push_back( NodeID( *it_sV, _VNODE ) );
        found = true;
        break;
      }
      if (type == XEXP_INSTANCE)
      {
        dsVars.push_back( NodeID( *it_sV, _DNODE ) );
        found = true;
        break;
      }
    }
    // Punt if the node has not been found, hope that the graph only 
    // has one node associated with this name.
    if (!found)
    {
      dsVars.push_back( NodeID( *it_sV, -1 ) );
    }
  } 
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerDepStoreGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerDepStoreGIDs(
		const vector< vector<int> > & dsGIDs )
{
#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "reg secondary gids: " << id_ << endl;

  for( unsigned int i = 0; i < dsGIDs.size(); ++i )
  {
    for( unsigned int j = 0; j < dsGIDs[i].size(); ++j )
      cout << " " << dsGIDs[i][j];
    cout << endl;
  }

  cout << endl;
#endif

  devPtr_->registerDepStoreGIDs( dsGIDs );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::getRowColPairs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::getRowColPairs( list<index_pair> & rcList )
{
  devPtr_->getIndexPairList( rcList );

#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "Index list for: " << id_ << endl;

  list<index_pair>::iterator it_ipL = rcList.begin();
  list<index_pair>::iterator it_ipL_end = rcList.end();
  for( ; it_ipL != it_ipL_end; ++it_ipL )
  {
    cout << "  " << it_ipL->row << " " << it_ipL->col;
  }

  cout << endl << endl;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::put
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
ostream& N_TOP_CktNode_Dev::put(ostream& os) const
{
  os << "CN_Dev: " << id_ << endl;
  os << "   GID: " << gID_ << "  Proc: " << procNum_ << endl;
  os << "   Owned: " << isOwned_ << endl;
  os << "   Offset: " << Offset_ << endl;
  os << "   Soln Var GID List: ";
  int count=0;
  list<int>::const_iterator it_iL = solnVarGIDList_.begin();
  list<int>::const_iterator it_iL_end = solnVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL )
  {
    os << *it_iL << "  ";
    if (count >= 12) {os << endl;count=0;}
    else ++count;
  }
  os << endl;

  os << "   State Var GID List: ";
  it_iL = stateVarGIDList_.begin();
  it_iL_end = stateVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL )
  {
    os << *it_iL << "  ";
  }

  os << "   Store Var GID List: ";
  it_iL = storeVarGIDList_.begin();
  it_iL_end = storeVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL )
  {
    os << *it_iL << "  ";
  }

  return os << endl;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/07/01
//-----------------------------------------------------------------------------
map<int,string> & N_TOP_CktNode_Dev::getIntNameMap ()
{
  return devPtr_->getIntNameMap ();
}


//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::getStateNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 08/01/2012
//-----------------------------------------------------------------------------
map<int,string> & N_TOP_CktNode_Dev::getStateNameMap ()
{
  return devPtr_->getStateNameMap ();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 08/01/2012
//-----------------------------------------------------------------------------
map<int,string> & N_TOP_CktNode_Dev::getStoreNameMap ()
{
  return devPtr_->getStoreNameMap ();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/20/02
//-----------------------------------------------------------------------------
const vector< vector<int> > & N_TOP_CktNode_Dev::jacobianStamp() const
{
  return devPtr_->jacobianStamp();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerJacLIDswithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/02
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerJacLIDswithDev( const vector< vector<int> > & jacLIDVec )
{
  devPtr_->registerJacLIDs( jacLIDVec );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::registerGIDDataWithDev
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 12/13/04
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::registerGIDDataWithDev(
                        const vector<int> & counts,
                        const vector<int> & GIDs,
                        const vector< vector<int> > & jacGIDs )
{
  devPtr_->registerGIDData( counts, GIDs, jacGIDs );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode_Dev::varTypeList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/02
//-----------------------------------------------------------------------------
void N_TOP_CktNode_Dev::varTypeList( vector<char> & varTypeVec )
{
  devPtr_->varTypes( varTypeVec );
}

