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
// Filename       : $RCSfile: N_TOP_Topology.C,v $
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
// Revision Number: $Revision: 1.121.2.2 $
//
// Revision Date  : $Date: 2014/03/06 17:23:43 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <sstream>

#include <fstream>

#include <N_TOP_Topology.h>

#include <N_TOP_TopoLSUtil.h>
#include <N_TOP_Directory.h>

#include <N_TOP_Indexor.h>

#include <N_TOP_CktNodeCreator.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_Dev.h>
#include <N_TOP_NodeBlock.h>
#include <N_TOP_NodeDevBlock.h>

#include <N_TOP_CktGraphSupport.h>
#include <N_TOP_CktGraphCreator.h>
#include <N_TOP_CktGraph.h>


#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>

#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#endif

#include <N_ANP_AnalysisInterface.h>

#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceState.h>

#include <N_UTL_OptionBlock.h>
#include <N_UTL_Functors.h>

#include <N_IO_RestartNode.h>
#include <N_IO_CmdParse.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : Topology::Topology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Topology::Topology(IO::CmdParse & cp)
  : commandLine_(cp),
    maxTries_(1),
    lsUtilPtr_(0),
    dirPtr_(0),
    mainGraphPtr_(0),
    graphCreatorPtr_(0),
    nodeCreatorPtr_(0),
    devIntPtr_(0),
    pkgOptMgrPtr_(0),
    pdsMgrPtr_(0),
    anaIntPtr_(0),
    icSettings_(0)
{
  lsUtilPtr_ = new TopoLSUtil( this, commandLine_ );

  // check for maximum tries to compute graph center on command line.
  if ( commandLine_.argExists ("-maxgraphcentertries") )
  {
    maxTries_ =
      atoi( commandLine_.getArgumentValue( "-maxgraphcentertries" ).c_str() );
  }

  // check for graph type command-line override:
  std::string graphType("Basic");
  graphCreatorPtr_ = CktGraphSupport::factory( graphType, maxTries_ );

  nodeCreatorPtr_ = CktNodeCreator::instance();

  mainGraphPtr_ = graphCreatorPtr_->create( std::string("") );

  graphList_.push_back(mainGraphPtr_);
}

//-----------------------------------------------------------------------------
// Function      : Topology::Topology
// Purpose       : generate topology with empty graph and put on list
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Topology::Topology(Device::DeviceInterface* devInt, IO::CmdParse & cp)
  : commandLine_(cp),
    maxTries_(1),
    lsUtilPtr_(0),
    dirPtr_(0),
    mainGraphPtr_(0),
    graphCreatorPtr_(0),
    nodeCreatorPtr_(0),
    devIntPtr_(devInt),
    pkgOptMgrPtr_(0),
    pdsMgrPtr_(0),
    anaIntPtr_(0),
    icSettings_(0)
{
  lsUtilPtr_ = new TopoLSUtil( this, commandLine_ );

  // check for maximum tries to compute graph center on command line.
  if ( commandLine_.argExists ("-maxgraphcentertries") )
  {
    maxTries_ =
      atoi( commandLine_.getArgumentValue( "-maxgraphcentertries" ).c_str() );
  }

  // check for graph type command-line override:
  std::string graphType("Basic");
  graphCreatorPtr_ = CktGraphSupport::factory( graphType, maxTries_ );

  mainGraphPtr_ = graphCreatorPtr_->create( std::string("") );

  graphList_.push_back(mainGraphPtr_);
}

//-----------------------------------------------------------------------------
// Function      : Topology::Topology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Topology::Topology(Device::DeviceInterface * devInt,
		const std::string & graphType, IO::CmdParse & cp)
  : commandLine_(cp),
    maxTries_(1),
    lsUtilPtr_(0),
    dirPtr_(0),
    mainGraphPtr_(0),
    graphCreatorPtr_(0),
    nodeCreatorPtr_(0),
    devIntPtr_(devInt),
    pkgOptMgrPtr_(0),
    pdsMgrPtr_(0),
    anaIntPtr_(0),
    icSettings_(0)
{
  lsUtilPtr_ = new TopoLSUtil( this, commandLine_ );

  // check for maximum tries to compute graph center on command line.
  if ( commandLine_.argExists ("-maxgraphcentertries") )
  {
    maxTries_ =
      atoi( commandLine_.getArgumentValue( "-maxgraphcentertries" ).c_str() );
  }

  graphCreatorPtr_ = CktGraphSupport::factory( graphType, maxTries_ );

  mainGraphPtr_ = graphCreatorPtr_->create( std::string("") );

  graphList_.push_back( mainGraphPtr_ );
}

//-----------------------------------------------------------------------------
// Function      : Topology::~Topology
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Topology::~Topology()
{
  delete lsUtilPtr_;

  if( icSettings_ ) delete icSettings_;

  std::list<CktGraph*>::iterator it_gL, end_gL;
  for ( it_gL = graphList_.begin(), end_gL = graphList_.end();
	it_gL != end_gL; ++it_gL )
    delete *it_gL;

  for( std::map<std::string,Device::InstanceBlock*>::iterator it_ibM =
	devInstBlockMap_.begin(); it_ibM != devInstBlockMap_.end(); ++it_ibM )
    if( it_ibM->second ) delete it_ibM->second;

  clearMigrateNodeMap();

  if( dirPtr_ ) delete dirPtr_;

  if (graphCreatorPtr_) delete graphCreatorPtr_;
  if (nodeCreatorPtr_) delete nodeCreatorPtr_;
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerDeviceInterface
// Purpose       : Basic device processor addes dev node and connected v-nodes.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/8/00
//-----------------------------------------------------------------------------
bool Topology::registerDeviceInterface ( Device::DeviceInterface * devInt)
{
  return (devIntPtr_ = devInt);
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerParallelMgr
// Purpose       : Registers parallel mgr.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/1/01
//-----------------------------------------------------------------------------
bool Topology::registerParallelMgr ( N_PDS_Manager * pdsmgr)
{
  if( (pdsMgrPtr_ = pdsmgr))
    return lsUtilPtr_->registerParallelMgr( pdsmgr );
  else
    return false;
}

//-----------------------------------------------------------------------------
// Function      : Topology::registeranaInt
// Purpose       : Basic device processor addes dev node and connected v-nodes
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/8/00
//-----------------------------------------------------------------------------
bool Topology::registeranaInt ( N_ANP_AnalysisInterface * anaInt)
{
  return (anaIntPtr_ = anaInt);
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerICs
// Purpose       : Registers parallel mgr
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/1/01
//-----------------------------------------------------------------------------
bool Topology::registerICs( const N_UTL_OptionBlock & ob )
{
  if( icSettings_ ) delete icSettings_;
  icSettings_ = new N_UTL_OptionBlock( ob );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool Topology::registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  lsUtilPtr_->registerPkgOptionsMgr( pkgOptMgrPtr_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Topology::setupGlobalIndices
// Purpose       : Lin system data setup
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/1/01
//-----------------------------------------------------------------------------
bool Topology::setupGlobalIndices()
{
  return lsUtilPtr_->setupRowCol();
}

//-----------------------------------------------------------------------------
// Function      : Topology::setupGlobalAccessors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/1/01
//-----------------------------------------------------------------------------
bool Topology::setupGlobalAccessors()
{
  return lsUtilPtr_->setupGlobalAccessors();
}

//-----------------------------------------------------------------------------
// Function      : Topology::addVoltageNode
// Purpose       : New v-node instantiator (planarized_ ckts only)
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/20/00
//-----------------------------------------------------------------------------
void Topology::addVoltageNode( const NodeBlock & nb )
{
  std::list<NodeID> emptyList;

  mainGraphPtr_->InsertNode( nodeCreatorPtr_->CreateVoltageNode( nb ),
                             emptyList);
}

//-----------------------------------------------------------------------------
// Function      : Topology::addDevice
// Purpose       : New dev-node instantiator (planarized ckts only)
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/20/00
//-----------------------------------------------------------------------------
void Topology::addDevice( const NodeBlock & nb, const Teuchos::RefCountPtr<Device::InstanceBlock> ibPtr )
{
  std::list<NodeID> emptyNLList, nlList;
  //------ Add connected voltage nodes

  std::list<tagged_param>::const_iterator it_nlL = nb.get_NodeList().begin();
  std::list<tagged_param>::const_iterator end_nlL = nb.get_NodeList().end();

  for( ; it_nlL != end_nlL ; ++it_nlL )
  {
    //----- insert each v-node: Unless the v-node already exists,
    //----- it will be instantiated but not owned and ProcNum_
    //----- identifies the owning processor
    mainGraphPtr_->InsertNode( nodeCreatorPtr_->CreateVoltageNode( it_nlL->tag ),
                               emptyNLList );

    nlList.push_back( NodeID( it_nlL->tag, _VNODE ) );
  }

  //------- Now instantiate the device node,  DistribMgr has
  //------- already set ownership

  devInstMap_[ibPtr->getName()] = 1;
  mainGraphPtr_->InsertNode( nodeCreatorPtr_->CreateDeviceNode( nb, ibPtr, *devIntPtr_ ),
                             nlList );
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerGIDswithDevs
// Purpose       : register the int. and ext. global ids stored by
//                 the cktnodes with their respective devices
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/22/00
//-----------------------------------------------------------------------------
void Topology::registerGIDswithDevs()
{
  mainGraphPtr_->registerGIDswithDevs();
  mainGraphPtr_->registerStateGIDswithDevs();
  mainGraphPtr_->registerStoreGIDswithDevs();
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerLIDswithDevs
// Purpose       : register the int. and ext. local ids stored by
//                 the cktnodes with their respective devices
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
void Topology::registerLIDswithDevs()
{
  if( commandLine_.getArgumentValue( "-dva" ) != "off" )
  {
    Indexor indexor( pdsMgrPtr_ );

    mainGraphPtr_->registerLIDswithDevs( indexor );
    mainGraphPtr_->registerStateLIDswithDevs( indexor );
    mainGraphPtr_->registerStoreLIDswithDevs( indexor );

    mainGraphPtr_->registerDepLIDswithDevs( indexor );
    mainGraphPtr_->registerDepStateLIDswithDevs( indexor );
    mainGraphPtr_->registerDepStoreLIDswithDevs( indexor );
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::registerJacLIDswithDevs
// Purpose       : register the jacobian local offsets with devices
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/02
//-----------------------------------------------------------------------------
void Topology::registerJacLIDswithDevs()
{
  if( commandLine_.getArgumentValue( "-dma" ) != "off" )
  {
    Indexor indexor( pdsMgrPtr_ );

    mainGraphPtr_->registerJacLIDswithDevs( indexor );
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::resolveDependentVars
// Purpose       : loop through devices and resolve their secondary
//                 dependencies
// Special Notes : Used to resolve Expressions and Current Dependencies
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/05/01
//-----------------------------------------------------------------------------
void Topology::resolveDependentVars()
{
  int count = 0;
  std::vector<int> locVec;
  std::vector<NodeID> nidVec;
  std::vector<NodeID> idVec;
  std::list<CktNode*>::iterator it_cnL;
  std::list<CktNode*>::iterator it_cnL_end;

  it_cnL = mainGraphPtr_->getBFSNodeList()->begin();
  it_cnL_end = mainGraphPtr_->getBFSNodeList()->end();

  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    (*it_cnL)->getDepSolnVars( idVec );
    if( !idVec.empty() )
    {
      locVec.push_back( count );
      for( unsigned int i = 0; i < idVec.size(); ++i )
        nidVec.push_back( idVec[i] );
      count += idVec.size();
    }
  }
  locVec.push_back( count );

  std::vector< std::vector<int> > gidVec, indexVec;
  std::vector<int> procVec;
  dirPtr_->getSolnGIDs( nidVec, gidVec, procVec );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Resolution of Dependent Solution Variables!\n";
  for( unsigned int i = 0; i < nidVec.size(); ++i )
  {
    Xyce::dout() << " Var: " << nidVec[i].first;
    Xyce::dout() << " Proc: " << procVec[i] << std::endl;
    for( unsigned int ii = 0; ii < gidVec[i].size(); ++ii )
      Xyce::dout() << " " << gidVec[i][ii];
    Xyce::dout() << std::endl;
  }
#endif

  count = 0;
  it_cnL = mainGraphPtr_->getBFSNodeList()->begin();
  it_cnL_end = mainGraphPtr_->getBFSNodeList()->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    (*it_cnL)->getDepSolnVars( idVec );
    if( !idVec.empty() )
    {
      indexVec.assign( gidVec.begin()+locVec[count],
                       gidVec.begin()+locVec[count+1] );
      (*it_cnL)->set_DepSolnGIDVec( indexVec );
      (*it_cnL)->registerDepSolnGIDs( indexVec );
      ++count;
    }
  }

  depSolnGIDMap_.clear();
  for( unsigned int i = 0; i < gidVec.size(); ++i )
    for( unsigned int ii = 0; ii < gidVec[i].size(); ++ii )
      depSolnGIDMap_[ gidVec[i][ii] ] = procVec[i];

  count = 0;
  locVec.clear();
  nidVec.clear();
  it_cnL = mainGraphPtr_->getBFSNodeList()->begin();
  it_cnL_end = mainGraphPtr_->getBFSNodeList()->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    (*it_cnL)->getDepStateVars( idVec );
    if( !idVec.empty() )
    {
      locVec.push_back( count );
      for( unsigned int i = 0; i < idVec.size(); ++i )
        nidVec.push_back( idVec[i] );
      count += idVec.size();
    }
  }
  locVec.push_back( count );

  gidVec.clear();
  procVec.clear();
  dirPtr_->getStateGIDs( nidVec, gidVec, procVec );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Resolution of Dependent State Variables!\n";
  for( unsigned int i = 0; i < nidVec.size(); ++i )
  {
    Xyce::dout() << " Var: " << nidVec[i].first << " Proc: " << procVec[i] << std::endl;
    for( unsigned int ii = 0; ii < gidVec[i].size(); ++ii )
      Xyce::dout() << " " << gidVec[i][ii];
    Xyce::dout() << std::endl;
  }
#endif

  count = 0;
  it_cnL = mainGraphPtr_->getBFSNodeList()->begin();
  it_cnL_end = mainGraphPtr_->getBFSNodeList()->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    (*it_cnL)->getDepStateVars( idVec );
    if( !idVec.empty() )
    {
      indexVec.assign( gidVec.begin()+locVec[count],
                       gidVec.begin()+locVec[count+1] );
      (*it_cnL)->set_DepStateGIDVec( indexVec );
      (*it_cnL)->registerDepStateGIDs( indexVec );
      ++count;
    }
  }

  depStateGIDMap_.clear();
  for( unsigned int i = 0; i < gidVec.size(); ++i )
  {
    for( unsigned int ii = 0; ii < gidVec[i].size(); ++ii )
    {
      depStateGIDMap_[ gidVec[i][ii] ] = procVec[i];
    }
  }

  count = 0;
  locVec.clear();
  nidVec.clear();
  it_cnL = mainGraphPtr_->getBFSNodeList()->begin();
  it_cnL_end = mainGraphPtr_->getBFSNodeList()->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    (*it_cnL)->getDepStoreVars( idVec );
    if( !idVec.empty() )
    {
      locVec.push_back( count );
      for( unsigned int i = 0; i < idVec.size(); ++i )
        nidVec.push_back( idVec[i] );
      count += idVec.size();
    }
  }
  locVec.push_back( count );

  gidVec.clear();
  procVec.clear();
  dirPtr_->getStoreGIDs( nidVec, gidVec, procVec );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Resolution of Dependent Store Variables!\n";
  for( unsigned int i = 0; i < nidVec.size(); ++i )
  {
    Xyce::dout() << " Var: " << nidVec[i].first << " Proc: " << procVec[i] << std::endl;
    for( unsigned int ii = 0; ii < gidVec[i].size(); ++ii )
      Xyce::dout() << " " << gidVec[i][ii];
    Xyce::dout() << std::endl;
  }
#endif

  count = 0;
  it_cnL = mainGraphPtr_->getBFSNodeList()->begin();
  it_cnL_end = mainGraphPtr_->getBFSNodeList()->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    (*it_cnL)->getDepStoreVars( idVec );
    if( !idVec.empty() )
    {
      indexVec.assign( gidVec.begin()+locVec[count],
                       gidVec.begin()+locVec[count+1] );
      (*it_cnL)->set_DepStoreGIDVec( indexVec );
      (*it_cnL)->registerDepStoreGIDs( indexVec );
      ++count;
    }
  }

  depStoreGIDMap_.clear();
  for( unsigned int i = 0; i < gidVec.size(); ++i )
  {
    for( unsigned int ii = 0; ii < gidVec[i].size(); ++ii )
    {
      depStoreGIDMap_[ gidVec[i][ii] ] = procVec[i];
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::OutputBFSGraphLists
// Purpose       : Output to Xyce::dout() BFS node list for debugging
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
void Topology::OutputBFSGraphLists()
{
  std::list<CktGraph*>::iterator it_cgL, end_cgL;
  std::list<CktNode*>::iterator it_cnL, end_cnL;

  Xyce::dout() << "BFS Node Listing for Graphs" << std::endl;

  //------ Loop over ckts in graphList_
  for( it_cgL = graphList_.begin(), end_cgL = graphList_.end();
	it_cgL != end_cgL; ++it_cgL )
  {
    std::list<CktNode*> * tmpCNL = (*it_cgL)->getBFSNodeList();
    for( it_cnL = tmpCNL->begin(), end_cnL = tmpCNL->end();
	it_cnL != end_cnL; ++it_cnL )
    {
      Xyce::dout() << *(*it_cnL) << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::OutputDFSGraphLists
// Purpose       : Output to Xyce::dout() DFS nodelist for debugging
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
void Topology::OutputDFSGraphLists()
{
  std::list<CktGraph*>::iterator it_cgL, end_cgL;
  std::list<CktNode*>::iterator it_cnL, end_cnL;

  Xyce::dout() << "DFS Node Listing for Graphs" << std::endl;

  for( it_cgL = graphList_.begin(), end_cgL = graphList_.end();
	it_cgL != end_cgL; ++it_cgL )
  {
    std::list<CktNode*> * tmpCNL = (*it_cgL)->getDFSNodeList();
    for( it_cnL = tmpCNL->begin(), end_cnL = tmpCNL->end();
         it_cnL != end_cnL; ++it_cnL )
    {
      Xyce::dout() << *(*it_cnL) << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::setOrderedNodeList
// Purpose       : Currently sets orderedNodeListPtr_ attribute to
//                 BFS traversal of main ckt.
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
void Topology::setOrderedNodeList() const
{
  orderedNodeListPtr_ = mainGraphPtr_->getBFSNodeList();
}

//-----------------------------------------------------------------------------
// Function      : verifyNodesAndDevices -- scans existing topology
//                 and asks device manager to verify each device.
//                 if the device manager returns false on any verification,
//                 then the nodes on that device are added to a list of nodes
//                 to be supernoded and the device will later be removed as redundant
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 2/18/2010
//-----------------------------------------------------------------------------
void Topology::verifyNodesAndDevices()
{
  if( lsUtilPtr_->supernodeFlag() )
  {
    int badDeviceCount=0;
    std::ostringstream outputStringStream;

    const std::map< NodeID, CktNode* > dataMap = mainGraphPtr_->getNodeList();
    std::map<NodeID, CktNode*>::const_iterator currentCktNodeItr = dataMap.begin();
    std::map<NodeID, CktNode*>::const_iterator endCktNodeItr = dataMap.end();
    while( currentCktNodeItr != endCktNodeItr )
    {
      if( ((*currentCktNodeItr).second)->type() == _DNODE )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>((*currentCktNodeItr).second);
        Teuchos::RCP<Device::InstanceBlock> deviceInstanceBlockPtr = cktNodeDevPtr->devBlock();
        if( Teuchos::nonnull( deviceInstanceBlockPtr ) )
        {
          bool deviceInstanceOk = devIntPtr_->verifyDeviceInstance( *deviceInstanceBlockPtr );
          if( !deviceInstanceOk )
          {
            badDeviceCount++;
            // place the nodes for this device in the superNodeList
            // collected nodes so that they can be combined (supernoded)
            // It's ok to let the device node get created and inserted into the
            // topology.  we can remove it later
            NodeID deviceID = (*currentCktNodeItr).first;

            std::vector< NodeID > adjacentIDs;
            mainGraphPtr_->returnAdjIDs( deviceID, adjacentIDs );

            std::vector< NodeID >::iterator currentID = adjacentIDs.begin();
            std::vector< NodeID >::iterator endID = adjacentIDs.end();
            std::vector< NodeID >::iterator nextID = currentID;
            nextID++;

            while( nextID != endID )
            {
              if ( (*currentID).first != (*nextID).first )
              {
                // these nodes are to be supernoded.  Take the lexically smaller one
                // as the node to keep (i.e. A < B )
                if( (*currentID).first < (*nextID).first )
                {
                  // format of pair is nodeToReplace, nodeToKeep
                  superNodeList.push_back( make_pair( *nextID, *currentID ) );
                }
                else
                {
                  // format of pair is nodeToReplace, nodeToKeep
                  superNodeList.push_back( make_pair( *currentID, *nextID ) );
                }
              }
              currentID++;
              nextID++;
            }
          }

        }
        else
        {
          // issue fatal error as this case shouldn't occur
          std::string msg("Topology::verifyNodesAndDevices() null Device Instance Block pointer.");
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        }
      }
      currentCktNodeItr++;
    }

    outputStringStream << "Device verification found " << badDeviceCount << " device(s) to remove";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, outputStringStream.str());
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::removeTaggedNodesAndDevices
// Purpose       : Remove devices and nodes that were tagged for removal
//                 during parsing.  Node removal is done through supernoding,
//                 .i.e. replacing one node globally with another.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 2/2/2010
//-----------------------------------------------------------------------------
void Topology::removeTaggedNodesAndDevices()
{
  if( lsUtilPtr_->supernodeFlag() )
  {

#ifdef Xyce_DEBUG_TOPOLOGY
    Xyce::dout() << "Topology::removeTaggedNodesAndDevices" << std::endl;
    Xyce::dout() << *this << std::endl;
#endif

    std::ostringstream outputStringStream;
    // storage for oldNode that we'll delete when done with this routine
    // use a set because we can get the same oldNode CktNode pointer
    // multiple times
    std::list<CktNode *> oldNodeList;

    // print out current state of supernode list
    std::list< std::pair<NodeID, NodeID> >::iterator currentNodePair = superNodeList.begin();
    std::list< std::pair<NodeID, NodeID> >::iterator endNodePair = superNodeList.end();
    std::set<NodeID> nodesReplaced;
    while ( currentNodePair != endNodePair )
    {
      NodeID nodeToBeReplaced( currentNodePair->first );
      NodeID replacementNode( currentNodePair->second );
      if( nodeToBeReplaced != replacementNode)
      {
        std::set<NodeID>::iterator nodesReplacedEnd = nodesReplaced.end();
        std::set<NodeID>::iterator nodeReplacedLoc = nodesReplaced.find( nodeToBeReplaced );
        if( nodeReplacedLoc == nodesReplacedEnd )
        {
          // this node hasn't been done before so replace it
          outputStringStream.str("");
          outputStringStream << "Replacing node \"" << nodeToBeReplaced << "\" with \"" << replacementNode << "\"";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, outputStringStream.str());

#ifdef Xyce_DEBUG_TOPOLOGY
          if( nodeToBeReplaced.first < replacementNode.first )
          {
            Xyce::dout() << "Ordering is wrong on nodes!" << std::endl;
          }
#endif
          CktNode * oldNode = mainGraphPtr_->replaceNode( nodeToBeReplaced, replacementNode );
          nodesReplaced.insert( nodeToBeReplaced );

          // If we delete this old node now, we'll make the orderedNodeListPtr_ untraversable.
          // this would be ok as we can regenerate it, but that takes time and we would only
          // invalidate it when we do our next delete.  So, store up the oldNode so we can
          // delete them when were done
          oldNodeList.push_back( oldNode );

          // now that we've replaced "nodeToBeReplaced" with "replacementNode" we need to
          // search the superNodeList from this point on also doing this same substitution
          //
          // for example if our super node list was:
          //
          //    B   A
          //    C   B
          //
          // If after the B->A substitution we didn't update our list, then we would next bring
          // back the B's with C->B.
          //
          std::list< std::pair<NodeID, NodeID> >::iterator nextNodePair = currentNodePair;
          nextNodePair++;
          while ( nextNodePair != endNodePair )
          {
            if( nextNodePair->first == nodeToBeReplaced )
            {
              if( replacementNode.first < (nextNodePair->second).first )
              {
                // need to swap on insert
                superNodeList.insert( nextNodePair, make_pair( nextNodePair->second, replacementNode ));
              }
              else
              {
                // just insert in order
                superNodeList.insert( nextNodePair, make_pair( replacementNode, nextNodePair->second ));
              }
              // store iterator to new pair before erasing old pair
              std::list< std::pair<NodeID, NodeID> >::iterator nextNodePairTmp = nextNodePair;
              nextNodePairTmp--;
              superNodeList.erase( nextNodePair );
              nextNodePair = nextNodePairTmp;
            }
            else if( nextNodePair->second == nodeToBeReplaced )
            {
              if( replacementNode.first < (nextNodePair->first).first )
              {
                // no swap needed
                superNodeList.insert( nextNodePair, make_pair( nextNodePair->first, replacementNode ));
              }
              else
              {
                // need to swap
                superNodeList.insert( nextNodePair, make_pair( replacementNode, nextNodePair->first ));
              }
              // store iterator to new pair before erasing old pair
              std::list< std::pair<NodeID, NodeID> >::iterator nextNodePairTmp = nextNodePair;
              nextNodePairTmp--;
              superNodeList.erase( nextNodePair );
              nextNodePair = nextNodePairTmp;
            }
            nextNodePair++;
          }
        }
      }
      currentNodePair++;
    }

    if( nodesReplaced.size() > 0 )
    {
      outputStringStream.str("");
      outputStringStream << "After combining equivalent nodes, " << nodesReplaced.size() << " nodes were removed." << std::endl;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, outputStringStream.str());
    }

    std::vector< CktNode * > removedDevices;
    mainGraphPtr_->removeRedundantDevices(removedDevices);

    // now it's safe to delete the old nodes and device nodes
    std::list< CktNode * >::iterator currentOldNodeItr = oldNodeList.begin();
    std::list< CktNode * >::iterator endOldNodeItr = oldNodeList.end();
    while( currentOldNodeItr != endOldNodeItr )
    {
      delete *currentOldNodeItr;
      currentOldNodeItr++;
    }

    if( removedDevices.size() > 0 )
    {
      outputStringStream.str("");
      outputStringStream <<  "After removing devices connected to one terminal, " << removedDevices.size() << " devices were removed." << std::endl;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, outputStringStream.str());

      // delete old devices that were removed
      std::vector< CktNode * >::iterator currentCktNodeItr = removedDevices.begin();
      std::vector< CktNode * >::iterator endCktNodeItr = removedDevices.end();
      while( currentCktNodeItr != endCktNodeItr )
      {
        CktNode_Dev * cktNodeDevPtr = dynamic_cast<CktNode_Dev*>(*currentCktNodeItr);
        Teuchos::RCP<Device::InstanceBlock> deviceInstanceBlockPtr = cktNodeDevPtr->devBlock();
        if( Teuchos::nonnull(deviceInstanceBlockPtr ))
        {
          // have a valid device.
#ifdef Xyce_DEBUG_TOPOLOGY
          std::string deviceID = (*currentCktNodeItr)->get_id();
          Xyce::dout() << "Device id: \"" << deviceID << "\"" << std::endl;
#endif
          // delete the device node
          delete cktNodeDevPtr;
        }
        currentCktNodeItr++;
      }
    }
  }
}


#ifdef Xyce_PARALLEL_MPI
//-----------------------------------------------------------------------------
// Function      : Topology::mergeOffProcTaggedNodesAndDevices
// Purpose       : Merge the off processor superNodeList and communicate
//                 the same list to all procs so topology reduction is the
//                 same on all procs
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 2/2/2010
//-----------------------------------------------------------------------------
void Topology::mergeOffProcTaggedNodesAndDevices()
{
  if( lsUtilPtr_->supernodeFlag() )
  {
    N_PDS_Comm * commPtr = pdsMgrPtr_->getPDSComm();
    int numProcs = commPtr->numProc();
    int thisProc = commPtr->procID();

    // Set up a pointer to tell us what processor the supernodes come from.
    std::vector<int> procPtr(numProcs+1);
    procPtr[0] = 0;

    int localNodes = superNodeList.size();

    // Count the bytes for packing these strings
    int byteCount = 0;

    // First count the size of the local superNodeList
    byteCount += sizeof(int);

    // Now count all the NodeIDs in the list
    for (std::list< std::pair<NodeID, NodeID> >::iterator nodePair = superNodeList.begin(); nodePair != superNodeList.end(); nodePair++) {
      byteCount += ((nodePair->first).first).length() + 2*sizeof(int);
      byteCount += ((nodePair->second).first).length() + 2*sizeof(int);
    }

    for( int p = 0; p < numProcs; ++p )
    {
      commPtr->barrier();

      // Broadcast the buffer size for this processor.
      int bsize=0;
      if (p==thisProc) { bsize = byteCount; }
      commPtr->bcast( &bsize, 1, p );

      // Create buffer.
      int pos = 0;
      char * superNodeBuffer = new char[bsize];

      if (p==thisProc) {

        // Pack number of supernodes on this processor and place in globalSuperNodeList
        commPtr->pack( &localNodes, 1, superNodeBuffer, bsize, pos );
        for (std::list< std::pair<NodeID, NodeID> >::iterator nodePair = superNodeList.begin(); nodePair != superNodeList.end(); nodePair++) {
          // Pack first entry of pair.
          int length = ((nodePair->first).first).size();
          commPtr->pack( &length, 1, superNodeBuffer, bsize, pos );
          commPtr->pack( ((nodePair->first).first).c_str(), length, superNodeBuffer, bsize, pos );
          commPtr->pack( &((nodePair->first).second), 1, superNodeBuffer, bsize, pos );

          // Pack second entry of pair.
          length = ((nodePair->second).first).size();
          commPtr->pack( &length, 1, superNodeBuffer, bsize, pos );
          commPtr->pack( ((nodePair->second).first).c_str(), length, superNodeBuffer, bsize, pos );
          commPtr->pack( &((nodePair->second).second), 1, superNodeBuffer, bsize, pos );
        }
        // Update the processor pointer.
        procPtr[p+1] = procPtr[p] + localNodes;

        // Broadcast packed buffer.
        commPtr->bcast( superNodeBuffer, bsize, p );

      }
      else {

        // Unpack buffer and place in globalSuperNodeList
        commPtr->bcast( superNodeBuffer, bsize, p );

        // Get number of supernodes from that processor
        int numSuperNodes = 0;
        commPtr->unpack( superNodeBuffer, bsize, pos, &numSuperNodes, 1 );

        // Update the processor pointer.
        procPtr[p+1] = procPtr[p] + numSuperNodes;

        // Extract supernode pairs and push to the back of the globalSuperNodeList
        int length=0;
        for (int i=0; i<numSuperNodes; ++i) {

          // Unpack first pair.
          commPtr->unpack( superNodeBuffer, bsize, pos, &length, 1 );
          std::string first_val( (superNodeBuffer+pos), length );
          pos += length;
          int first_type = 0;
          commPtr->unpack( superNodeBuffer, bsize, pos, &first_type, 1 );

          // Unpack second pair.
          commPtr->unpack( superNodeBuffer, bsize, pos, &length, 1 );
          std::string second_val( (superNodeBuffer+pos), length );
          pos += length;
          int second_type = 0;
          commPtr->unpack( superNodeBuffer, bsize, pos, &second_type, 1 );

          // Push to back of globalSuperNodeList
          globalSuperNodeList.push_back( make_pair( NodeID(first_val,first_type), NodeID(second_val,second_type) ) );
        }
      }
      // Clean up.
      delete [] superNodeBuffer;
    }

    // Now go through the local superNodeList and tack on any boundary cases, where adjacencies may effect device removal.
    for (std::list< std::pair<NodeID, NodeID> >::iterator nodePair = superNodeList.begin(); nodePair != superNodeList.end(); nodePair++)
    {
      NodeID nodeToBeRemoved = nodePair->first;
      NodeID replacementNode = nodePair->second;
      std::list< std::pair<NodeID, NodeID> >::iterator currentGlobalSN = globalSuperNodeList.begin();
      std::list< std::pair<NodeID, NodeID> >::iterator endGlobalSN = globalSuperNodeList.end();
      while ( currentGlobalSN != endGlobalSN )
      {
        // Add pair if replacement node is node to be removed by another processor.
        if (currentGlobalSN->first == replacementNode)
          superNodeList.push_back( *currentGlobalSN );
        // Add pair if node to be removed is node a replacement node on another processor.
        if (currentGlobalSN->second == nodeToBeRemoved)
          superNodeList.push_back( *currentGlobalSN );
        currentGlobalSN++;
      }
    }
  }
}

#endif

//-----------------------------------------------------------------------------
// Function      : Topology::instantiateDevices
// Purpose       : Delayed instantiation of devices
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/03
//-----------------------------------------------------------------------------
void Topology::instantiateDevices()
{
  setOrderedNodeList();

  std::list<CktNode*>::iterator iterCN = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator endCN = orderedNodeListPtr_->end();
  for( ; iterCN != endCN; ++iterCN )
  {
    CktNode_Dev * cnd = dynamic_cast<CktNode_Dev*>(*iterCN);
    if( cnd ) cnd->instantiate();
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnNodeGIDVec
// Purpose       : Generate ordered node global id vector using orderedNodeList
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/25/01
//-----------------------------------------------------------------------------
void Topology::returnNodeGIDVec( std::vector<int> & nodeGIDVec )
{
  setOrderedNodeList();

  nodeGIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  int i = 0;
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
    {
      ++i;
    }
  }

  nodeGIDVec.reserve( i );

  it_cnL = orderedNodeListPtr_->begin();
  it_cnL_end = orderedNodeListPtr_->end();
  for( ; it_cnL != it_cnL_end ; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
    {
      nodeGIDVec.push_back( (*it_cnL)->get_gID() );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnExternNodeGIDVec
// Purpose       : Generate ordered externnode global id vector using
//                 orderedNodeList.
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/23/01
//-----------------------------------------------------------------------------
void Topology::returnExternNodeGIDVec( std::vector< std::pair<int,int> >
							& nodeGIDVec )
{
  setOrderedNodeList();

  nodeGIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  int i = 0;
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( !( (*it_cnL)->get_IsOwned() ) && ( (*it_cnL)->get_gID() != -1 ) )
    {
      ++i;
    }
  }

  nodeGIDVec.reserve( i );

  it_cnL = orderedNodeListPtr_->begin();
  it_cnL_end = orderedNodeListPtr_->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( !( (*it_cnL)->get_IsOwned() ) && ( (*it_cnL)->get_gID() != -1 ) )
    {
      nodeGIDVec.push_back( std::pair<int,int>( (*it_cnL)->get_gID(),
					(*it_cnL)->get_ProcNum() ) );
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Topology::returnSVarGIDVec
// Purpose       : Generate ordered soln var global id vector using
//                 orderedNodeList.
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/25/01
//-----------------------------------------------------------------------------
void Topology::returnSVarGIDVec( std::vector<int> & sVarGIDVec )
{

  setOrderedNodeList();

  sVarGIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  int i = 0;
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
    {
      i += (*it_cnL)->get_SolnVarGIDList().size();
    }
  }

  sVarGIDVec.reserve( i );

  it_cnL = orderedNodeListPtr_->begin();
  it_cnL_end = orderedNodeListPtr_->end();

  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
      sVarGIDVec.insert( sVarGIDVec.end(),
        (*it_cnL)->get_SolnVarGIDList().begin(),
        (*it_cnL)->get_SolnVarGIDList().end() );
  }

}

//-----------------------------------------------------------------------------
// Function      : Topology::returnExternSVarGIDVec
// Purpose       : Generate ordered soln var global id vector using
//                 orderedNodeList for external nodes.
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/25/01
//-----------------------------------------------------------------------------
void Topology::returnExternSVarGIDVec(
		std::vector< std::pair<int,int> > & sVarGIDVec )
{
  setOrderedNodeList();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  sVarGIDVec.clear();

  int i = 0;
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( !( (*it_cnL)->get_IsOwned() ) && ( (*it_cnL)->get_gID() != -1 ) )
    {
      i += (*it_cnL)->get_SolnVarGIDList().size();
    }
  }

  sVarGIDVec.reserve( i );

  it_cnL = orderedNodeListPtr_->begin();
  it_cnL_end = orderedNodeListPtr_->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( !( (*it_cnL)->get_IsOwned() ) && ( (*it_cnL)->get_gID() != -1 ) )
    {
      std::list<int>::const_iterator it_svL =
       (*it_cnL)->get_SolnVarGIDList().begin();
      std::list<int>::const_iterator it_svL_end =
       (*it_cnL)->get_SolnVarGIDList().end();

      for(  ; it_svL != it_svL_end; ++it_svL )
      {
        sVarGIDVec.push_back( std::pair<int,int>( *it_svL,
                                  (*it_cnL)->get_ProcNum() ) );
      }
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Topology::returnStateVarGIDVec
// Purpose       : Generate ordered state var global id vector using
//                 orderedNodeList.
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/25/01
//-----------------------------------------------------------------------------
void Topology::returnStateVarGIDVec( std::vector<int> & sVarGIDVec )
{
  setOrderedNodeList();
  sVarGIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  int i = 0;
  for( ; it_cnL != it_cnL_end; ++it_cnL)
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
    {
      i += (*it_cnL)->get_StateVarGIDList().size();
    }
  }

  sVarGIDVec.reserve( i );

  it_cnL = orderedNodeListPtr_->begin();
  it_cnL_end = orderedNodeListPtr_->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
      sVarGIDVec.insert( sVarGIDVec.end(),
                        (*it_cnL)->get_StateVarGIDList().begin(),
                        (*it_cnL)->get_StateVarGIDList().end() );
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnStoreVarGIDVec
// Purpose       : Generate ordered store var global id vector using
//                 orderedNodeList.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Topology::returnStoreVarGIDVec( std::vector<int> & sVarGIDVec )
{
  setOrderedNodeList();
  sVarGIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  int i = 0;
  for( ; it_cnL != it_cnL_end; ++it_cnL)
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
    {
      i += (*it_cnL)->get_StoreVarGIDList().size();
    }
  }

  sVarGIDVec.reserve( i );

  it_cnL = orderedNodeListPtr_->begin();
  it_cnL_end = orderedNodeListPtr_->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
      sVarGIDVec.insert( sVarGIDVec.end(),
                        (*it_cnL)->get_StoreVarGIDList().begin(),
                        (*it_cnL)->get_StoreVarGIDList().end() );
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnExternStateVarGIDVec
// Purpose       : Generate ordered state var global id vector using
//                 orderedNodeList for external nodes.
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/25/01
//-----------------------------------------------------------------------------
void Topology::returnExternStateVarGIDVec(
		std::vector< std::pair<int,int> > & sVarGIDVec )
{
  setOrderedNodeList();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  sVarGIDVec.clear();

  int i = 0;
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( !( (*it_cnL)->get_IsOwned() ) && ( (*it_cnL)->get_gID() != -1 ) )
    {
      i += (*it_cnL)->get_StateVarGIDList().size();
    }
  }

  sVarGIDVec.reserve( i );

  it_cnL     = orderedNodeListPtr_->begin();
  it_cnL_end = orderedNodeListPtr_->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( !( (*it_cnL)->get_IsOwned() ) && ( (*it_cnL)->get_gID() != -1 ) )
    {
      std::list<int>::const_iterator it_svL =
       (*it_cnL)->get_StateVarGIDList().begin();
      std::list<int>::const_iterator it_svL_end =
       (*it_cnL)->get_StateVarGIDList().end();

      for( ; it_svL != it_svL_end; ++it_svL )
      {
        sVarGIDVec.push_back( std::pair<int,int>( *it_svL,
                                (*it_cnL)->get_ProcNum() ) );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnExternStoreVarGIDVec
// Purpose       : Generate ordered store var global id vector using
//                 orderedNodeList for external nodes.
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/25/01
//-----------------------------------------------------------------------------
void Topology::returnExternStoreVarGIDVec(
		std::vector< std::pair<int,int> > & sVarGIDVec )
{
  setOrderedNodeList();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  sVarGIDVec.clear();

  int i = 0;
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( !( (*it_cnL)->get_IsOwned() ) && ( (*it_cnL)->get_gID() != -1 ) )
    {
      i += (*it_cnL)->get_StoreVarGIDList().size();
    }
  }

  sVarGIDVec.reserve( i );

  it_cnL     = orderedNodeListPtr_->begin();
  it_cnL_end = orderedNodeListPtr_->end();
  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( !( (*it_cnL)->get_IsOwned() ) && ( (*it_cnL)->get_gID() != -1 ) )
    {
      std::list<int>::const_iterator it_svL =
       (*it_cnL)->get_StoreVarGIDList().begin();
      std::list<int>::const_iterator it_svL_end =
       (*it_cnL)->get_StoreVarGIDList().end();

      for( ; it_svL != it_svL_end; ++it_svL )
      {
        sVarGIDVec.push_back( std::pair<int,int>( *it_svL,
                                (*it_cnL)->get_ProcNum() ) );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnVarTypeVec
// Purpose       : Generate ordered list of variable types
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/25/01
//-----------------------------------------------------------------------------
void Topology::returnVarTypeVec( std::vector<char> & varTypeVec ) const
{
  setOrderedNodeList();

  varTypeVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();

  for( ; it_cnL != orderedNodeListPtr_->end(); ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 ) )
    {
      if( (*it_cnL)->type() == _VNODE )
        varTypeVec.push_back( 'V' );
      else if( (*it_cnL)->solnVarCount() )
      {
        std::vector<char> typeList;
        (*it_cnL)->varTypeList( typeList ); // types set by individual devices.

        if( typeList.empty() )
        {
          int cnt = (*it_cnL)->solnVarCount();
          for( int i = 0; i < cnt; ++i ) varTypeVec.push_back( 'V' );
        }
        else
        {
          int cnt = typeList.size();
          for( int i = 0; i < cnt; ++i ) varTypeVec.push_back( typeList[i] );
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnSVarVNodeGIDVec
// Purpose       : Generate ordered list of variable types
// Special Notes :
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/25/01
//-----------------------------------------------------------------------------
void Topology::returnSVarVNodeGIDVec( std::vector<int> & SVarVNodeGIDVec )
{
  setOrderedNodeList();

  SVarVNodeGIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 )
         && ( (*it_cnL)->type() == _VNODE ) )
    {
      SVarVNodeGIDVec.push_back( *((*it_cnL)->get_SolnVarGIDList().begin()) );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnSVarVsrcGIDVec
// Purpose       : Generate ordered vector of variables connected to vsrcs.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
void Topology::returnSVarVsrcGIDVec( std::vector<int> & SVarVsrcGIDVec )
{
  setOrderedNodeList();

  SVarVsrcGIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    int owned = (*it_cnL)->get_IsOwned();
    if (owned)
    {
      int type = (*it_cnL)->type();
      if (type == _DNODE)
      {
        const std::string & id = (*it_cnL)->get_id();
        std::string::size_type col = id.find_first_of(':');

        if ( id[col+1] == 'V' || id[col+1] == 'v' )
        {
          std::list<int>::const_iterator iterIL = (*it_cnL)->get_SolnVarGIDList().begin();
          std::list<int>::const_iterator endIL = (*it_cnL)->get_SolnVarGIDList().end();
          for( ; iterIL != endIL; ++iterIL )
          {
            SVarVsrcGIDVec.push_back( *iterIL );
          }
          std::list<int>::const_iterator iterEIL = (*it_cnL)->get_ExtSolnVarGIDList().begin();
          std::list<int>::const_iterator endEIL = (*it_cnL)->get_ExtSolnVarGIDList().end();
          for( ; iterEIL != endEIL; ++iterEIL )
          {
            SVarVsrcGIDVec.push_back( *iterEIL );
          }
        }
      }
    }
  }

#ifdef Xyce_DEBUG_TOPOLOGY
    int isize= SVarVsrcGIDVec.size();
    for (int ieric=0;ieric<isize;++ieric)
    {
      Xyce::dout() << "SVarVsrcGIDVec["<<ieric<<"] = " << SVarVsrcGIDVec[ieric] << std::endl;
    }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Topology::returnSVarNoDCPathIDVec
// Purpose       : Generate ordered list of solution IDs with no dc path to
//                 ground.
// Special Notes :
// Scope         : public
// Creator       : Keith Santarelli, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/27/07
//-----------------------------------------------------------------------------
void Topology::returnSVarNoDCPathIDVec( std::vector<std::string> &
					       SVarNoDCPathIDVec )
{
  setOrderedNodeList();

  SVarNoDCPathIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 )
         && ( (*it_cnL)->type() == _VNODE ) )
    {
      if ((*it_cnL)->getNoDCPathVar())
      {
        SVarNoDCPathIDVec.push_back((*it_cnL)->get_id());
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::returnSVarConnToOneTermIDVec
// Purpose       : Generate ordered list of solution vars that are only
//                 connected to one terminal.
// Scope         : public
// Creator       : Keith Santarelli, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/27/07
//-----------------------------------------------------------------------------
void Topology::returnSVarConnToOneTermIDVec( std::vector<std::string> &
					       SVarConnToOneTermIDVec )
{
  setOrderedNodeList();

  SVarConnToOneTermIDVec.clear();

  std::list<CktNode*>::iterator it_cnL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator it_cnL_end = orderedNodeListPtr_->end();

  for( ; it_cnL != it_cnL_end; ++it_cnL )
  {
    if( (*it_cnL)->get_IsOwned() && ( (*it_cnL)->get_gID() != -1 )
         && ( (*it_cnL)->type() == _VNODE ) )
    {
      if ((*it_cnL)->getConnToOneTermVar())
      {
        SVarConnToOneTermIDVec.push_back((*it_cnL)->get_id());
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::getNodeSVarGIDs
// Purpose       : Return list of solution var indices for named node.
// Special Notes : returns false if node not owned or not local.
// Scope         : public
// Creator       : Rob  Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/11/00
//-----------------------------------------------------------------------------
bool Topology::getNodeSVarGIDs( const NodeID& id,
		                      std::list<int> & sVarGIDList,
		                      std::list<int> & extSVarGIDList,
                                      char & type )
{
  CktNode * cnPtr = mainGraphPtr_->FindCktNode( id );

  if( cnPtr != NULL )
  {
    if( cnPtr->type() == _DNODE ) type = 'D';
    else                          type = 'V';

    sVarGIDList.assign( cnPtr->get_SolnVarGIDList().begin(),
                        cnPtr->get_SolnVarGIDList().end() );
    extSVarGIDList.assign( cnPtr->get_ExtSolnVarGIDList().begin(),
                        cnPtr->get_ExtSolnVarGIDList().end() );
    if( cnPtr->get_IsOwned() )
    {
      return true;
    }
    else
    {
      sVarGIDList.clear();
      return false;
    }
  }
  else
    return false;
}

//-----------------------------------------------------------------------------
// Function      : Topology::extractMigrateNodes
// Purpose       : Extracts nodes and devices to be migrated.
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
void Topology::extractMigrateNodes( const int & num,
			const int * nodeGIDs, const int * procs )
{
  CktNode * nodePtr;
  NodeBlock * nodeBlockPtr;
  Device::InstanceBlock * devBlockPtr;
  NodeDevBlock * migrateNodePtr;

  bool devFlag = false;

  for( int i = 0; i < num; ++i )
  {
    devFlag = false;

    nodePtr = mainGraphPtr_->FindCktNode( nodeGIDs[i] );

    nodePtr->set_ProcNum( procs[i] );
    nodeBlockPtr = nodePtr->extractNodeBlock();
    nodePtr->set_IsOwned( false );

    if( devInstBlockMap_.find( nodePtr->get_id() ) != devInstBlockMap_.end() )
      devFlag = true;

    std::list<int> gidList, svGIDList, procList;
    std::list<NodeID> idList;
    mainGraphPtr_->returnAdjNodesWithGround( nodeGIDs[i], gidList,
		svGIDList, procList, idList );

    if( devFlag )
    {
      devBlockPtr = devInstBlockMap_[ nodePtr->get_id() ];

      std::list<tagged_param> nList, npList;
      std::list<NodeID>::iterator it_idL = idList.begin();
      std::list<int>::iterator it_gidL = gidList.begin();
      std::list<int>::iterator it_pL = procList.begin();
      for( ; it_idL != idList.end(); ++it_idL, ++it_gidL, ++it_pL )
      {
        nList.push_back( tagged_param( (*it_idL).first, *it_gidL ) );
        npList.push_back( tagged_param( (*it_idL).first, *it_pL ) );
      }

      nodeBlockPtr->set_NodeList( nList );
      nodeBlockPtr->set_NodeProcList( npList );

      migrateNodePtr = new NodeDevBlock( *nodeBlockPtr, *devBlockPtr );
    }
    else
    {
      migrateNodePtr = new NodeDevBlock( *nodeBlockPtr );
      migrateNodePtr->getDevBlock().setName("");
    }

    if( migrateNodeMap_[ nodeGIDs[i] ].find( nodePtr->get_gID() )
	       == migrateNodeMap_[ nodeGIDs[i] ].end() )
    {
      migrateNodeMap_[ nodeGIDs[i] ][ nodePtr->get_gID() ] = migrateNodePtr;
    }
    else
    {
      delete migrateNodePtr;
    }

    delete nodeBlockPtr;

    //Add neighbors to migrateNodeMap_
    std::list<int>::iterator it_iL = gidList.begin();
    for( ; it_iL != gidList.end(); ++it_iL )
    {
      nodePtr = mainGraphPtr_->FindCktNode( *it_iL );

      if( *it_iL == -1 )
        nodeBlockPtr = new NodeBlock( "0", std::list<tagged_param>(),
                                    std::list<tagged_param>(), false, -1 );
      else
        nodeBlockPtr = nodePtr->extractNodeBlock();

      nodeBlockPtr->set_IsOwned( false );

      std::list<int> gidList2;

      if( !devFlag )
      {
        devBlockPtr = devInstBlockMap_[ nodePtr->get_id() ];

        std::list<int> svGIDList2, procList2;
        std::list<NodeID> idList2;
        mainGraphPtr_->returnAdjNodesWithGround( *it_iL, gidList2,
            svGIDList2, procList2, idList2 );

        std::list<tagged_param> nList2, npList2;
        std::list<NodeID>::iterator it_idL2 = idList2.begin();
        std::list<int>::iterator it_gidL2 = gidList2.begin();
        std::list<int>::iterator it_pL2 = procList2.begin();
        for( ; it_idL2 != idList2.end(); ++it_idL2, ++it_gidL2, ++it_pL2 )
        {
          nList2.push_back( tagged_param( (*it_idL2).first, *it_gidL2 ) );
          npList2.push_back( tagged_param( (*it_idL2).first, *it_pL2 ) );
        }

        nodeBlockPtr->set_NodeList( nList2 );
        nodeBlockPtr->set_NodeProcList( npList2 );

        migrateNodePtr = new NodeDevBlock( *nodeBlockPtr, *devBlockPtr );
      }
      else
      {
        migrateNodePtr = new NodeDevBlock( *nodeBlockPtr );
        migrateNodePtr->getDevBlock().setName("");
      }

      if( migrateNodeMap_[ nodeGIDs[i] ].find( nodePtr->get_gID() )
          == migrateNodeMap_[ nodeGIDs[i] ].end() )
      {
        migrateNodeMap_[ nodeGIDs[i] ][ nodePtr->get_gID() ] = migrateNodePtr;
      }
      else
      {
        delete migrateNodePtr;
      }

      delete nodeBlockPtr;

      //if initial node is voltage node, add neighbors neighbors
      if( !devFlag )
      {
        for( std::list<int>::iterator it_iL2 = gidList2.begin();
              it_iL2 != gidList2.end(); ++it_iL2 )
        {
          nodePtr = mainGraphPtr_->FindCktNode( *it_iL2 );

          nodeBlockPtr = nodePtr->extractNodeBlock();
          nodeBlockPtr->set_IsOwned( false );

          migrateNodePtr = new NodeDevBlock( *nodeBlockPtr );
          migrateNodePtr->getDevBlock().setName("");

          if( migrateNodeMap_[ nodeGIDs[i] ].find( nodePtr->get_gID() )
	      == migrateNodeMap_[ nodeGIDs[i] ].end() )
            migrateNodeMap_[ nodeGIDs[i] ][ nodePtr->get_gID() ] =
						migrateNodePtr;
          else
            delete migrateNodePtr;

          delete nodeBlockPtr;
        }
      }
    }

    std::map<int,NodeDevBlock*>::iterator it_ndbM =
	    migrateNodeMap_[ nodeGIDs[i] ].begin();
    std::map<int,NodeDevBlock*>::iterator it_ndbM_end =
      migrateNodeMap_[nodeGIDs[i] ].end();

    Xyce::dout() << "Extracting Migrate Node: " << nodeGIDs[i] << std::endl;

    for( ; it_ndbM != it_ndbM_end; ++it_ndbM )
    {
      Xyce::dout() << it_ndbM->second->getNodeBlock() << std::endl;
    }
    Xyce::dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::insertMigrateNodes
// Purpose       : Inserts nodes and devices from migration Map.
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
void Topology::insertMigrateNodes( const int & num,
			const int * nodeGIDs, const int * procs )
{
  std::map<std::string,int> foundMap;
  std::map<std::string,int> ownMap;

  for( int i = 0; i < num; ++i )
  {
    std::map<int,NodeDevBlock *>::iterator it_ndbM =
          migrateNodeMap_[ nodeGIDs[i] ].begin();
    std::map<int,NodeDevBlock *>::iterator it_ndbM_end =
          migrateNodeMap_[ nodeGIDs[i] ].end();

    for( ; it_ndbM != it_ndbM_end ; ++it_ndbM )
    {
      if( it_ndbM->second->isDevice() )
        if( it_ndbM->second->getNodeBlock().get_IsOwned() )
          if( ownMap.find( it_ndbM->second->getID() ) == ownMap.end() )
            ownMap[it_ndbM->second->getID()] = 1;
    }
  }

  for( int i = 0; i < num; ++i )
  {
    Xyce::dout() << "Inserting Migrate Node: " << nodeGIDs[i] << std::endl;

    std::map<int,NodeDevBlock *>::iterator it_ndbM =
      migrateNodeMap_[ nodeGIDs[i] ].begin();

    std::map<int,NodeDevBlock *>::iterator it_ndbM_end =
      migrateNodeMap_[ nodeGIDs[i] ].end();

    for( ; it_ndbM != it_ndbM_end; ++it_ndbM )
    {
      Xyce::dout() << it_ndbM->second->getNodeBlock() << std::endl;
      Xyce::dout() << it_ndbM->second->isDevice() << std::endl;
      if( it_ndbM->second->isDevice() )
      {
        if( ownMap.find( it_ndbM->second->getID() ) != ownMap.end() &&
            it_ndbM->second->getNodeBlock().get_IsOwned() )
        {
          Teuchos::RCP<Device::InstanceBlock> ibRcp = rcp(  new Device::InstanceBlock(it_ndbM->second->getDevBlock()) );
          addDevice( it_ndbM->second->getNodeBlock(),ibRcp );
        }
        else if( foundMap.find( it_ndbM->second->getID() ) == foundMap.end() )
        {
          foundMap[ it_ndbM->second->getID() ] = 1;
          Teuchos::RCP<Device::InstanceBlock> ibRcp = rcp( new Device::InstanceBlock(it_ndbM->second->getDevBlock()) );
          addDevice( it_ndbM->second->getNodeBlock(), ibRcp );
        }
      }
      else
      {
        Xyce::dout() << "Adding Voltage Node: " << it_ndbM->second->getID() << std::endl;
        addVoltageNode( it_ndbM->second->getNodeBlock() );
      }
    }
    //Xyce::dout() << std::endl << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Topology::clearMigrateNodeMap
// Purpose       : Clears migration map.
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
void Topology::clearMigrateNodeMap()
{
  for( std::map<int,std::map<int,NodeDevBlock *> >::iterator it_ndbM =
	migrateNodeMap_.begin(); it_ndbM != migrateNodeMap_.end(); ++it_ndbM )
    for( std::map<int,NodeDevBlock *>::iterator it_ndbM2 =
	it_ndbM->second.begin(); it_ndbM2 != it_ndbM->second.end(); ++it_ndbM2 )
      if( it_ndbM2->second != NULL ) delete it_ndbM2->second;
}

//-----------------------------------------------------------------------------
// Function      : Topology::addMigrateNode
// Purpose       : Add node to migration node map.
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
void Topology::addMigrateNode( const int & id,
		std::map<int,NodeDevBlock *> & ndbL )
{
  migrateNodeMap_[id] = ndbL;
}

//-----------------------------------------------------------------------------
// Function      : Topology::getMigrateNode
// Purpose       : Get node from migration node map.
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
std::map<int,NodeDevBlock *> & Topology::getMigrateNode( const int & id )
{
  return migrateNodeMap_[id];
}

//-----------------------------------------------------------------------------
// Function      : Topology::pruneCkt
// Purpose       : Prune out unnecessary nodes.
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
void Topology::pruneCkt( const int & num, const int * nodeGIDs )
{
  CktNode * nodePtr;

  for( int i = 0; i < num; ++i )
  {
    nodePtr = mainGraphPtr_->FindCktNode( nodeGIDs[i] );

    if( nodePtr != NULL )
    {
      if( devInstBlockMap_.find( nodePtr->get_id() ) != devInstBlockMap_.end() )
      {
        pruneDevNode( nodeGIDs[i] );
      }
      else
      {
        std::list<int> gidList, svGIDList, procList;
        std::list<NodeID> idList;
        mainGraphPtr_->returnAdjNodes( nodeGIDs[i], gidList, svGIDList,
                                       procList, idList );

        std::list<int>::iterator it_iL = gidList.begin();
        std::list<int>::iterator it_iL_end = gidList.end();
        for( ; it_iL != it_iL_end; ++it_iL )
        {
          if( !(mainGraphPtr_->FindCktNode( *it_iL )->get_IsOwned()) )
          {
            pruneDevNode( *it_iL );
          }
        }
      }
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Topology::pruneDevNode
// Purpose       : prune out unnecessary nodes
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
void Topology::pruneDevNode( const int & node )
{
  CktNode * nodePtr;

  std::list<int> gidList, svGIDList, procList;
  std::list<NodeID> idList;
  mainGraphPtr_->returnAdjNodes( node, gidList, svGIDList, procList,
		idList );

  bool needed = false;

  nodePtr = mainGraphPtr_->FindCktNode( node );

  if( nodePtr != NULL )
  {
    std::list<int>::iterator it_iL = gidList.begin();
    std::list<int>::iterator it_iL_end = gidList.end();
    for( ; it_iL != it_iL_end; ++it_iL )
    {
      if( mainGraphPtr_->FindCktNode( *it_iL )->get_IsOwned() ) needed = true;
    }

    if( !needed )
    {
      nodePtr = mainGraphPtr_->ExtractNode( node );
      devIntPtr_->deleteDeviceInstance( nodePtr->get_id() );
      devInstMap_.erase( nodePtr->get_id() );
      delete devInstBlockMap_[nodePtr->get_id()];
      devInstBlockMap_.erase(nodePtr->get_id());
      delete nodePtr;

      //Check Neighboring V-nodes
      std::list<int>::iterator it_iL = gidList.begin();
      std::list<int>::iterator it_iL_end = gidList.end();
      for( ; it_iL != it_iL_end; ++it_iL )
      {
        if( mainGraphPtr_->returnNumEdges( *it_iL ) == 0 )
        {
          nodePtr = mainGraphPtr_->ExtractNode( *it_iL );
          delete nodePtr;
        }
      }
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Topology::regenerateGIDNodeMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/6/01
//-----------------------------------------------------------------------------
void Topology::regenerateGIDNodeMap()
{
  mainGraphPtr_->regenerateGIDNodeMap();
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : implementation for debugging
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const Topology &topo)
{
  std::list<CktGraph*>::const_iterator it_gL, end_gL;

  for ( it_gL = topo.graphList_.begin(), end_gL = topo.graphList_.end();
				 it_gL != end_gL; ++it_gL )
    os << *(*it_gL) << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Topology::generateICLoader
// Purpose       : SuperNode all nodes associated with voltage sources
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/20/01
//-----------------------------------------------------------------------------
bool Topology::generateICLoader()
{
  std::vector< std::pair<int,double> > * vecP = NULL;
  if( icSettings_ != NULL )
  {
    Xyce::dout() << *icSettings_;

    vecP = new std::vector< std::pair<int,double> >();

    CktNode * cnp;
    for( std::list<N_UTL_Param>::iterator iterPL = icSettings_->getParams().begin();
         iterPL != icSettings_->getParams().end(); ++iterPL )
    {
      cnp = mainGraphPtr_->FindCktNode( NodeID( iterPL->tag(), -1 ) );
      if( cnp != NULL )
        if( cnp->get_IsOwned() )
          vecP->push_back( std::pair<int,double>(
		*(cnp->get_SolnVarGIDList().begin()), iterPL->getImmutableValue<double>() ) );
    }

    return devIntPtr_->registerICLoads( vecP );
  }
  else
    return true;
}

//-----------------------------------------------------------------------------
// Function      : Topology::getRestartNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/29/01
//-----------------------------------------------------------------------------
bool Topology::getRestartNodes( std::vector<IO::RestartNode*> & nodeVec )
{
  if( orderedNodeListPtr_ == NULL ) return false;

  int count = 0;
  std::list<CktNode*>::iterator iterCN = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator endCN = orderedNodeListPtr_->end();
  for( ; iterCN != endCN; ++iterCN )
    if( (*iterCN)->get_IsOwned() && (*iterCN)->get_gID() != -1 ) ++count;

  nodeVec.resize(count);

  CktNode * cnP;
  iterCN = orderedNodeListPtr_->begin();
  count = 0;
  for( iterCN = orderedNodeListPtr_->begin(); iterCN != endCN; ++iterCN )
    if( (*iterCN)->get_IsOwned() && (*iterCN)->get_gID() != -1 )
    {
      cnP = *iterCN;
      nodeVec[count] = new IO::RestartNode( cnP->get_id(), cnP->type() );

      int i = 0;
      const std::list<int> & gidList = cnP->get_SolnVarGIDList();
      nodeVec[count]->solnVarData.resize( gidList.size() );
      for( std::list<int>::const_iterator iterIC = gidList.begin();
           iterIC != gidList.end(); ++iterIC, ++i )
        anaIntPtr_->getSolnVarData( *iterIC, nodeVec[count]->solnVarData[i] );

      i = 0;
      const std::list<int> & gidList2 = cnP->get_StateVarGIDList();
      nodeVec[count]->stateVarData.resize( gidList2.size() );
      for( std::list<int>::const_iterator iterIC = gidList2.begin();
           iterIC != gidList2.end(); ++iterIC, ++i )
        anaIntPtr_->getStateVarData( *iterIC, nodeVec[count]->stateVarData[i] );

      i = 0;
      const std::list<int> & gidList3 = cnP->get_StoreVarGIDList();
      nodeVec[count]->storeVarData.resize( gidList3.size() );
      for( std::list<int>::const_iterator iterIC = gidList3.begin();
           iterIC != gidList3.end(); ++iterIC, ++i )
        anaIntPtr_->getStoreVarData( *iterIC, nodeVec[count]->storeVarData[i] );

      if( cnP->type() == _DNODE )
        nodeVec[count]->devState =
         (dynamic_cast<CktNode_Dev*>(cnP))->getDevState();

      ++count;
    }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Topology::restoreRestartNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/29/01
//-----------------------------------------------------------------------------
bool Topology::restoreRestartNodes( std::vector<IO::RestartNode*> & nodeVec )
{
  CktNode * cnP;

  for( unsigned int i = 0; i < nodeVec.size(); ++i )
  {
    cnP = mainGraphPtr_->FindCktNode( NodeID(nodeVec[i]->ID,nodeVec[i]->type) );

    if( cnP != NULL )
    {
#ifdef Xyce_DEBUG_RESTART
      Xyce::dout() << "Restoring Node: " << nodeVec[i]->ID << std::endl;
#endif

      const std::list<int> & gidList = cnP->get_SolnVarGIDList();
      int pos = 0;
      for( std::list<int>::const_iterator iterIC = gidList.begin();
           iterIC != gidList.end(); ++iterIC, ++pos )
        anaIntPtr_->setSolnVarData( *iterIC, nodeVec[i]->solnVarData[pos] );

      const std::list<int> & gidList2 = cnP->get_StateVarGIDList();
      pos = 0;
      for( std::list<int>::const_iterator iterIC = gidList2.begin();
           iterIC != gidList2.end(); ++iterIC, ++pos )
        anaIntPtr_->setStateVarData( *iterIC, nodeVec[i]->stateVarData[pos] );

      const std::list<int> & gidList3 = cnP->get_StoreVarGIDList();
      pos = 0;
      for( std::list<int>::const_iterator iterIC = gidList3.begin();
           iterIC != gidList3.end(); ++iterIC, ++pos )
        anaIntPtr_->setStoreVarData( *iterIC, nodeVec[i]->storeVarData[pos] );

      if( nodeVec[i]->devState != NULL && nodeVec[i]->devState != 0 )
        (dynamic_cast<CktNode_Dev*>(cnP))->setDevState( *nodeVec[i]->devState );
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Topology::generateDirectory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/29/01
//-----------------------------------------------------------------------------
bool Topology::generateDirectory()
{
  dirPtr_ = new Directory( this, pdsMgrPtr_ );
  dirPtr_->generateDirectory();

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Topology::outputNameFile
// Purpose       : This is a kludgy function designed to output all the
//                 solution variable indices and their respective names
//                 to a file.
//
//                 The original point was to create a file to compare
//                 with SPICE, so the names needed to be as similar
//                 as possible to SPICE's naming convention.
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/07/01
//-----------------------------------------------------------------------------
bool Topology::outputNameFile (bool overRideOutput)
{
  if( lsUtilPtr_->namesFileFlag() || overRideOutput )
  {
    std::string netListFile("");
    if (commandLine_.getArgumentValue("netlist") != "")
    {
      netListFile = commandLine_.getArgumentValue("netlist");
    }

    std::string fileName("namesMap_"+netListFile+".txt");

    int masterRank = 0;
    int numProcs = 1;
    int thisProc = 0;

#ifdef Xyce_PARALLEL_MPI
    N_PDS_Comm * commPtr = pdsMgrPtr_->getPDSComm();
    numProcs = commPtr->numProc();
    thisProc = commPtr->procID();
#endif

    for( int p = 0; p < numProcs; ++p )
    {
#ifdef Xyce_PARALLEL_MPI
      commPtr->barrier();
#endif

      if (p==thisProc)
      {
        std::ostream * outStreamPtr_;
        if( masterRank == thisProc )
        {
          outStreamPtr_ = new std::ofstream(fileName.c_str());
        }
        else
        {
          outStreamPtr_ = new std::ofstream(fileName.c_str(),std::ios_base::app);
        }

        if(outStreamPtr_->fail())
        {
          if (p==0)
          {
            Report::UserWarning() << "Unable to open names file" <<std::endl;
          }
        }
        else
        {
          if (p==0)
          {
            // to line up with the outputted files of N_LAS_Vector, which start with N.
            (*outStreamPtr_) << "HEADER" << std::endl;
          }

          // Loop over the graph.  If this is a voltage node, just output
          // the name and solution variable index.  If a device node,
          // call the device package to print out names of internal nodes.

          std::list<CktNode*>::iterator iterONL = orderedNodeListPtr_->begin();
          std::list<CktNode*>::iterator endONL  = orderedNodeListPtr_->end  ();

          for( ; iterONL != endONL; ++iterONL )
          {
            int owned = (*iterONL)->get_IsOwned();
            if (owned)
            {
              // get node type:
              int type = (*iterONL)->type();

              // if voltage node, just output first GID and the name.
              if (type == _VNODE )
              {

                // get GID list:
                std::list<int> gidList = (*iterONL)->get_SolnVarGIDList ();
                std::list<int>::iterator iterGID = gidList.begin ();
                std::list<int>::iterator endGID  = gidList.end   ();

                if ( (*iterGID) >= 0 )
                {
                  (*outStreamPtr_) << "\t";
                  // index
                  (*outStreamPtr_) << (*iterGID) << "\t";
                  // name
                  outStreamPtr_->width(12);
                  (*outStreamPtr_) << ExtendedString((*iterONL)->get_id()).toLower();
                  (*outStreamPtr_) <<std::endl;
                }
              }
              // if device node, call the device and have it output indices
              // of the internal variables, plus names.  (only the device
              // itself knows these names, I think).
              else if (type == _DNODE)
              {
#ifdef Xyce_PARALLEL_MPI
                Indexor indexor( pdsMgrPtr_ );
#endif
                std::string mapName("SOLUTION_OVERLAP_GND");
                std::vector<int> convIntVec(1,0);

                std::map<int,std::string> & inMap = (*iterONL)->getIntNameMap ();

                std::map<int,std::string>::iterator iterMap = inMap.begin();
                std::map<int,std::string>::iterator  endMap = inMap.end  ();

                for ( ;iterMap != endMap ; ++iterMap )
                {
                  convIntVec[0] = iterMap->first;
#ifdef Xyce_PARALLEL_MPI
                  bool success = indexor.localToGlobal( mapName, convIntVec );
#endif
                  if ( convIntVec[0] != -1 )
                  //if ( (iterMap->first) != -1 )
                  {
                    (*outStreamPtr_) << "\t";
                    (*outStreamPtr_) << convIntVec[0] << "\t"; // index
                    outStreamPtr_->width(12);
                    (*outStreamPtr_) <<ExtendedString(iterMap->second).toLower(); // name
                    (*outStreamPtr_) <<std::endl;
                  }
                }
              }
              else
                N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
                  "Trying to output something other than a vnode or dnode!\n" );
            }
          }

        }
        delete outStreamPtr_;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Topology::getNodeNames
// Purpose       : get node names for nodes owned by this processor
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date :
//-----------------------------------------------------------------------------
void Topology::getNodeNames (std::map<std::string, std::pair<int, double>, Xyce::LessNoCase > & nodeNames)
{
  std::list<CktNode*>::iterator iterONL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator endONL  = orderedNodeListPtr_->end  ();
  Indexor indexor( pdsMgrPtr_ );
  std::string mapName("SOLUTION_OVERLAP_GND");
  std::vector<int> intVec(1,0);

  // two other functions will do identical searches over the node list looking
  // for devices.  At this time we will fill up another list with just devices
  // to speed that later work.

  bool fillDeviceNodeList=false;
  if( deviceNodeListPtr_.empty() )
  {
    fillDeviceNodeList=true;
  }

  for( ; iterONL != endONL; ++iterONL )
  {
    int owned = (*iterONL)->get_IsOwned();
    if (owned)
    {
      int type = (*iterONL)->type();
      if (type == _VNODE )
      {
        std::list<int> gidList = (*iterONL)->get_SolnVarGIDList ();
        std::list<int>::iterator iterGID = gidList.begin ();
        if (*iterGID >= 0)
        {
          intVec[0] = *iterGID;
          indexor.globalToLocal( mapName, intVec );

          nodeNames[(*iterONL)->get_id()].first = intVec[0];
        }
      }
      else if (type == _DNODE)
      {
        // if the device node list was empty at the start,
        // store the pointers to the devices.
        if( fillDeviceNodeList )
        {
          deviceNodeListPtr_.push_back( *iterONL);
        }

        std::map<int,std::string> & inMap = (*iterONL)->getIntNameMap ();

        std::map<int,std::string>::iterator iterMap = inMap.begin();
        std::map<int,std::string>::iterator  endMap = inMap.end  ();

        for ( ;iterMap != endMap ; ++iterMap )
        {
          if( iterMap->first != -1 )
          {
            nodeNames[iterMap->second].first = iterMap->first;
          }
        }
      }
    }
  }

  return;
}


//-----------------------------------------------------------------------------
// Function      : Topology::getStateNodeNames
// Purpose       : get node names for state vector elements owned by this processor
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, SNL, Electrical Systems Modeling
// Creation Date : 08/01/2012
//-----------------------------------------------------------------------------
void Topology::getStateNodeNames (std::map<std::string, std::pair<int, double>, Xyce::LessNoCase > & nodeNames)
{
  if( deviceNodeListPtr_.empty() )
  {
    std::list<CktNode*>::iterator iterONL = orderedNodeListPtr_->begin();
    std::list<CktNode*>::iterator endONL  = orderedNodeListPtr_->end  ();

    for( ; iterONL != endONL; ++iterONL )
    {
      int owned = (*iterONL)->get_IsOwned();
      if (owned)
      {
        int type = (*iterONL)->type();
        // Only device nodes have state/store components.  Thus only querry those
        //if (type == _VNODE )
        //{
        //}
        if (type == _DNODE)
        {
          std::map<int,std::string> & inMap = (*iterONL)->getStateNameMap ();

          std::map<int,std::string>::iterator iterMap = inMap.begin();
          std::map<int,std::string>::iterator  endMap = inMap.end  ();

          for ( ;iterMap != endMap ; ++iterMap )
          {
            if( iterMap->first != -1 )
            {
              nodeNames[iterMap->second].first = iterMap->first;
            }
          }
        }
      }
    }
  }
  else
  {
    // we can do a slightly simpler loop here as we already have a list of
    // owned, device nodes.
    std::list<CktNode*>::iterator iterONL = deviceNodeListPtr_.begin();
    std::list<CktNode*>::iterator endONL  = deviceNodeListPtr_.end  ();

    for( ; iterONL != endONL; ++iterONL )
    {
      std::map<int,std::string> & inMap = (*iterONL)->getStateNameMap ();

      std::map<int,std::string>::iterator iterMap = inMap.begin();
      std::map<int,std::string>::iterator  endMap = inMap.end  ();

      for ( ;iterMap != endMap ; ++iterMap )
      {
        if( iterMap->first != -1 )
        {
          nodeNames[iterMap->second].first = iterMap->first;
        }
      }
    }
  }

  return;
}


//-----------------------------------------------------------------------------
// Function      : Topology::getStoreNodeNames
// Purpose       : get node names for store vector elements owned by this processor
// Special Notes :
// Scope         :
// Creator       : Richard Schiek, SNL, Electrical Systems Modeling
// Creation Date : 08/01/2012
//-----------------------------------------------------------------------------
void Topology::getStoreNodeNames (std::map<std::string, std::pair<int, double>, Xyce::LessNoCase > & nodeNames)
{
  if( deviceNodeListPtr_.size() == 0 )
  {
    std::list<CktNode*>::iterator iterONL = orderedNodeListPtr_->begin();
    std::list<CktNode*>::iterator endONL  = orderedNodeListPtr_->end  ();

    for( ; iterONL != endONL; ++iterONL )
    {
      int owned = (*iterONL)->get_IsOwned();
      if (owned)
      {
        int type = (*iterONL)->type();
        // Only device nodes have state/store components.  Thus only querry those
        //if (type == _VNODE )
        //{
        //}
        if (type == _DNODE)
        {
          std::map<int,std::string> & inMap = (*iterONL)->getStoreNameMap ();

          std::map<int,std::string>::iterator iterMap = inMap.begin();
          std::map<int,std::string>::iterator  endMap = inMap.end  ();

          for ( ;iterMap != endMap ; ++iterMap )
          {
            if( iterMap->first != -1 )
            {
              nodeNames[iterMap->second].first = iterMap->first;
            }
          }
        }
      }
    }
  }
  else
  {
    // we can do a slightly simpler loop here as we already have a list of
    // owned, device nodes.
    std::list<CktNode*>::iterator iterONL = deviceNodeListPtr_.begin();
    std::list<CktNode*>::iterator endONL  = deviceNodeListPtr_.end  ();

    for( ; iterONL != endONL; ++iterONL )
    {
      std::map<int,std::string> & inMap = (*iterONL)->getStoreNameMap ();

      std::map<int,std::string>::iterator iterMap = inMap.begin();
      std::map<int,std::string>::iterator  endMap = inMap.end  ();

      for ( ;iterMap != endMap ; ++iterMap )
      {
        if( iterMap->first != -1 )
        {
          nodeNames[iterMap->second].first = iterMap->first;
        }
      }
    }
  }
  return;

}

//-----------------------------------------------------------------------------
// Function      : Topology::getExternNodeNames
// Purpose       : get node names for nodes owned by this processor
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter
// Creation Date : 4/25/2012
//-----------------------------------------------------------------------------
void Topology::getExternNodeNames (std::map<std::string, std::pair<int, double>, Xyce::LessNoCase > & nodeNames)
{
  std::list<CktNode*>::iterator iterONL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator endONL  = orderedNodeListPtr_->end  ();
  Indexor indexor( pdsMgrPtr_ );
  std::string mapName("SOLUTION_OVERLAP_GND");
  std::vector<int> intVec(1,0);

  for( ; iterONL != endONL; ++iterONL )
  {
    int owned = (*iterONL)->get_IsOwned();
    if (owned)
    {
      int type = (*iterONL)->type();
      if (type == _VNODE )
      {
        std::list<int> gidList = (*iterONL)->get_SolnVarGIDList ();
        std::list<int>::iterator iterGID = gidList.begin ();
        if (*iterGID >= 0)
        {
          intVec[0] = *iterGID;
          indexor.globalToLocal( mapName, intVec );

          nodeNames[(*iterONL)->get_id()].first = intVec[0];
        }
      }
      else if (type == _DNODE)
      {
        // do nothing as we only want externs.
      }
    }
  }

  return;
}


//-----------------------------------------------------------------------------
// Function      : Topology::getVsrcNodes
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void Topology::getVsrcNodes (std::set<std::string> & vsrcSet)
{

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "In Topology::getVsrcNodes" << std::endl;
#endif

  if( deviceNodeListPtr_.size() == 0 )
  {
    std::list<CktNode*>::iterator iterONL = orderedNodeListPtr_->begin();
    std::list<CktNode*>::iterator endONL  = orderedNodeListPtr_->end  ();

    for( ; iterONL != endONL; ++iterONL )
    {
      int owned = (*iterONL)->get_IsOwned();
      if (owned)
      {
        int type = (*iterONL)->type();
        if (type == _DNODE)
        {
          std::vector<NodeID> adj_ids;
          const std::string & id = (*iterONL)->get_id();
          std::string::size_type col = id.find_first_of(':');

          if ( id[col+1] == 'V' || id[col+1] == 'v' )
          {
  #ifdef Xyce_DEBUG_TOPOLOGY
            Xyce::dout() << "Getting adjacent nodes for: " << id << std::endl;
  #endif
            mainGraphPtr_->returnAdjIDs( NodeID(id,type), adj_ids );
            int adjSize = adj_ids.size();
            for( int i = 0; i < adjSize; ++i )
            {
  #ifdef Xyce_DEBUG_TOPOLOGY
              Xyce::dout() << "adj_ids["<<i<<"] = " << adj_ids[i] << std::endl;
  #endif
              if( adj_ids[i].first != "0" )
              {
                vsrcSet.insert(adj_ids[i].first);
              }
            }
          }
        }
      }
    }
  }
  else
  {
    // slightly simpler loop as we already have a list of owned, device nodes
    std::list<CktNode*>::iterator iterONL = deviceNodeListPtr_.begin();
    std::list<CktNode*>::iterator endONL  = deviceNodeListPtr_.end  ();
    for( ; iterONL != endONL; ++iterONL )
    {
      std::vector<NodeID> adj_ids;
      const std::string & id = (*iterONL)->get_id();
      int type = (*iterONL)->type();
      std::string::size_type col = id.find_first_of(':');

      if ( id[col+1] == 'V' || id[col+1] == 'v' )
      {
#ifdef Xyce_DEBUG_TOPOLOGY
        Xyce::dout() << "Getting adjacent nodes for: " << id << std::endl;
#endif
        mainGraphPtr_->returnAdjIDs( NodeID(id,type), adj_ids );
        int adjSize = adj_ids.size();
        for( int i = 0; i < adjSize; ++i )
        {
#ifdef Xyce_DEBUG_TOPOLOGY
          Xyce::dout() << "adj_ids["<<i<<"] = " << adj_ids[i] << std::endl;
#endif
          if( adj_ids[i].first != "0" )
          {
            vsrcSet.insert(adj_ids[i].first);
          }
        }
      }
    }
  }
  return;
}


//-----------------------------------------------------------------------------
// Function      : Topology::getRawData
// Purpose       : get node ids, names, and types
// Special Notes :
// Scope         :
// Creator       : Eric Rankin
// Creation Date :
//-----------------------------------------------------------------------------
bool Topology::getRawData( std::map< int, std::string > & nRef, std::vector< char > & tRef )
{
  // collect internal soln var types
  std::vector< char > varTypeVec;
  returnVarTypeVec( tRef );

  std::list<CktNode*>::iterator iterONL = orderedNodeListPtr_->begin();
  std::list<CktNode*>::iterator endONL  = orderedNodeListPtr_->end  ();

  for( ; iterONL != endONL; ++iterONL )
  {
    // get node type:
    int type = (*iterONL)->type();

    // if voltage node, just retrieve the first GID and the name.
    if( type == _VNODE )
    {
      std::list<int> gidList = (*iterONL)->get_SolnVarGIDList ();

      if ( *gidList.begin() != -1 )
      {
        nRef[*gidList.begin()] = ExtendedString( (*iterONL)->get_id() ).toLower();
      }
    }

    // if device node, call the device and have it output indices
    // of the internal variables, plus names.
    else if( type == _DNODE )
    {
      std::map< int,std::string > & inMap = (*iterONL)->getIntNameMap();

      std::map<int,std::string>::iterator iterMap = inMap.begin();
      std::map<int,std::string>::iterator  endMap = inMap.end  ();

      for( ; iterMap != endMap ; ++iterMap )
      {
        if( iterMap->first != -1 )
        {
          nRef[iterMap->first] = ExtendedString( iterMap->second ).toLower();
        }
      }
    }

    else
    {
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
       "Trying to output something other than a vnode or dnode!\n" );
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Topology::addResistors
// Purpose       : Adds resistors (between ground and nodes which are connected
//                 to only one device terminal) to a copy of the netlist file.
// Special Notes :
// Scope         : public
// Creator       : Keith Santarelli, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/11/07
//-----------------------------------------------------------------------------
void Topology::addResistors
 (const std::vector<std::string> & inputVec,
  const std::string & resValue,
  const std::string & netlistFile, bool oneTermNotNoDCPath)
{
  std::string netlistCopy(netlistFile);
  netlistCopy += "_xyce.cir";
  std::ofstream copyFile;
  copyFile.open(netlistCopy.c_str(),std::ios::app);  //put in append mode

  std::string msg("");

#ifdef Xyce_DEBUG_IO
  if (oneTermNotNoDCPath)
  {
    msg = "Adding resistors of value ";
    //msg += resvalue;
    msg += " between ground and nodes connected to only one device";
    msg += "terminal in file ";
    msg += netlistCopy;
  }
  else
  {
    msg = "Adding resistors of value ";
    //msg += resvalue;
    msg += " between ground and nodes with no DC path to ground in file ";
    msg += netlistCopy;
  }

  Xyce::dout() << msg << std::endl;
#endif

  //Some error checking in case we can't open the file.
  if(copyFile.fail())
  {
    if (oneTermNotNoDCPath)
    {
      msg = "Error in adding resistors between ground and nodes connected";
      msg += " to only one device terminal:  cannot open file ";
      msg += netlistCopy;
    }
    else
    {
      msg = "Error in adding resistors between ground and nodes with no";
      msg += " DC path to ground:  cannot open file ";
      msg += netlistCopy;
    }

    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg);
  }


  //Create a header banner to say what is is we're adding.  Only add the banner
  // on proc 0.

  //#ifdef Xyce_PARALLEL_MPI
  //N_PDS_Comm & comm = *(pdsMgrPtr_->getPDSComm());
  //int procCnt = comm.numProc();
  //int procID = comm.procID();

  //if (procID == 0)
  //  {
  //#endif

  std::string banner("");

  if (oneTermNotNoDCPath)
  {
    banner = "*XYCE-GENERATED OUTPUT:  Adding resistors between ground and ";
    banner += "nodes connected to only 1 device terminal:";
  }
  else
  {
    banner = "*XYCE-GENERATED OUTPUT:  Adding resistors between ground and ";
    banner += "nodes with no DC path to ground:";
  }

  copyFile << std::endl << std::endl << banner << std::endl << std::endl;

  //#ifdef Xyce_PARALLEL_MPI
  //}
  //#endif

  //Now, loop through the ids in inputVec and add the resistors.

  std::vector<std::string>::const_iterator inputit = inputVec.begin();
  std::vector<std::string>::const_iterator inputend = inputVec.end();
  std::string node("");
  int count = 0;

  while (inputit != inputend)
  {
    node = *inputit;

    //If node is blank (""), then something has gone wrong---we're trying to
    //access a node that's not of type _VNODE (this *shouldn't* happen,
    //assuming that everything here is done correctly, but, in case it does,
    //we don't want to add junk lines to the copy of the netlist file.

    if (node == "")
    {
      msg = "Error in netlist copy:  attempt to access circuit node not";
      msg += " of type _VNODE.  No line will be printed in netlist copy";
      msg += " file.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg);
    }
    else
    {
      std::string resname("R");

      if (oneTermNotNoDCPath)
      {
        resname += "ONETERM";
      }
      else
      {
        resname += "NODCPATH";
      }

      copyFile << resname;
      copyFile << count+1 << " " << node << " 0 " << resValue;
      copyFile << std::endl;
    }
    inputit++;
    count++;
  }
  copyFile.close();
}

//-----------------------------------------------------------------------------
// Function      : Topology::appendEndStatement
// Purpose       : Adds ".END" to a copy of the netlist file.
// Special Notes :
// Scope         : public
// Creator       : Keith Santarelli, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/11/07
//-----------------------------------------------------------------------------
void Topology::appendEndStatement(const std::string & netlistFile)
{
  std::string netlistCopy(netlistFile);
  netlistCopy += "_xyce.cir";
  std::ofstream copyFile;
  copyFile.open(netlistCopy.c_str(),std::ios::app);  //put in append mode

  std::string msg("");

  //Some error checking in case we can't open the file.
  if(copyFile.fail())
  {
    msg = "Error in attempt to append .END statement as part of netlist";
    msg += " copy procedure:  Cannot open file ";
    msg += netlistCopy;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg);
  }

  copyFile << std::endl << ".END" << std::endl;
  copyFile.close();
}

} // namespace Topo
} // namespace Xyce
