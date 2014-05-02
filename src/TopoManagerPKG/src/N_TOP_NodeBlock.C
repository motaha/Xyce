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
// Filename       : $RCSfile: N_TOP_NodeBlock.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/23/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TOP_NodeBlock.h>
#include <N_PDS_Comm.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : NodeBlock::operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
std::ostream & operator<< ( std::ostream & os, const NodeBlock & nb )
{
  os << "NodeBlock: " << nb.id_ << std::endl;
  os << " Connected Nodes: ";
  std::list<tagged_param>::const_iterator it_tpL = nb.nodeList_.begin(); 
  std::list<tagged_param>::const_iterator it_tpL_end = nb.nodeList_.end();
  for( ; it_tpL != it_tpL_end; ++it_tpL)
    os << "   " << it_tpL->tag << " " << it_tpL->param;
  os << std::endl;

  os << " Node-Proc List: ";
  it_tpL = nb.nodeProcList_.begin(); 
  it_tpL_end = nb.nodeProcList_.end();
  for( ; it_tpL != it_tpL_end; ++it_tpL)
    os << "   " << it_tpL->tag << " " << it_tpL->param;
  os << std::endl;

  os << "  Proc Num: " << nb.procNum_ << std::endl;
  os << "  Global ID: " << nb.gID_ << std::endl;
  os << "  SV Global IDs: ";
  std::list<int>::const_iterator it_iL = nb.solnVarGIDList_.begin(); 
  std::list<int>::const_iterator it_iL_end = nb.solnVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL)
    os << "   " << *it_iL;
  os << std::endl;

  os << "  StateVar Global IDs: ";
  it_iL = nb.stateVarGIDList_.begin(); 
  it_iL_end = nb.stateVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL)
    os << "   " << *it_iL;
  os << std::endl;

  os << "  StoreVar Global IDs: ";
  it_iL = nb.storeVarGIDList_.begin(); 
  it_iL_end = nb.storeVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL)
    os << "   " << *it_iL;
  os << std::endl;

  os << " isOwned: " << nb.isOwned_ << std::endl;

  return os << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
Packable * NodeBlock::instance() const
{
  return new NodeBlock();
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::packedByteCount
// Purpose       : Counts number of bytes needed to pack object
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
int NodeBlock::packedByteCount() const
{
  int byteCount = 0;

  int length, size, i;

  std::list<tagged_param>::const_iterator it_tpL;

  //----- count id_
  length = id_.length();
  byteCount += sizeof( int );
  byteCount += length * sizeof( char );

  //----- count nodeList_
  size = nodeList_.size();
  byteCount += sizeof( int );
  it_tpL = nodeList_.begin();
  for( i = 0; i < size; ++i, ++it_tpL )
  {
    length = it_tpL->tag.length();
    byteCount += sizeof( int );
    byteCount += length * sizeof( char );
    byteCount += sizeof( double );
  }

  //----- count nodeProcList_
  size = nodeProcList_.size();
  byteCount += sizeof( int );
  it_tpL = nodeProcList_.begin();
  for( i = 0; i < size; ++i, ++it_tpL )
  {
    length = it_tpL->tag.length();
    byteCount += sizeof( int );
    byteCount += length;
    byteCount += sizeof( double );
  }

  //----- count isOwned_
  byteCount += sizeof( int );

  //----- count gID_
  byteCount += sizeof( int );

  //----- count procNum__
  byteCount += sizeof( int );

  //----- count solnVarGIDList_
  size = solnVarGIDList_.size();
  byteCount += sizeof( int );
  byteCount += ( size * sizeof( int ) );

  //----- pack stateVarGIDList_
  size = stateVarGIDList_.size();
  byteCount += sizeof( int );
  byteCount += ( size * sizeof( int ) );

  //----- pack storeVarGIDList_
  size = storeVarGIDList_.size();
  byteCount += sizeof( int );
  byteCount += ( size * sizeof( int ) );

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::pack
// Purpose       : Packs NodeBlock into char buffer using MPI_PACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/23/00
//-----------------------------------------------------------------------------
void NodeBlock::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
{

  int size, length;
  int i;
  std::list<tagged_param>::const_iterator it_tpL;
  std::list<int>::const_iterator it_iL;

  //----- pack id_
  length = id_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( id_.c_str(), length, buf, bsize, pos );

  //----- pack nodeList_
  size = nodeList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  it_tpL = nodeList_.begin();
  for( i = 0; i < size; ++i, ++it_tpL )
  {
    length = it_tpL->tag.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( it_tpL->tag.c_str(), length, buf, bsize, pos );
    comm->pack( &(it_tpL->param), 1, buf, bsize, pos );
  }

  //----- pack nodeProcList_
  size = nodeProcList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  it_tpL = nodeProcList_.begin();
  for( i = 0; i < size; ++i, ++it_tpL )
  {
    length = it_tpL->tag.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( it_tpL->tag.c_str(), length, buf, bsize, pos );
    comm->pack( &(it_tpL->param), 1, buf, bsize, pos );
  }

  //----- pack isOwned_
  i = isOwned_;
  comm->pack( &i, 1, buf, bsize, pos );

  //----- pack gID_
  comm->pack( &gID_, 1, buf, bsize, pos );

  //----- pack procNum__
  comm->pack( &procNum_, 1, buf, bsize, pos );

  int tmpInt;

  //----- pack solnVarGIDList_
  size = solnVarGIDList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  it_iL = solnVarGIDList_.begin();
  for( i = 0; i < size; ++i, ++it_iL )
  {
    tmpInt = (*it_iL);
    comm->pack( &tmpInt, 1, buf, bsize, pos );
  }

  //----- pack stateVarGIDList_
  size = stateVarGIDList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  it_iL = stateVarGIDList_.begin();
  for( i = 0; i < size; ++i, ++it_iL )
  {
    tmpInt = (*it_iL);
    comm->pack( &tmpInt, 1, buf, bsize, pos );
  }

  //----- pack storeVarGIDList_
  size = storeVarGIDList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  it_iL = storeVarGIDList_.begin();
  for( i = 0; i < size; ++i, ++it_iL )
  {
    tmpInt = (*it_iL);
    comm->pack( &tmpInt, 1, buf, bsize, pos );
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Packed " << pos << " bytes for NodeBlock: " << id_ << std::endl;
#endif

}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::unpack
// Purpose       : Unpacks NodeBlock from char buffer using MPI_UNPACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/23/00
//-----------------------------------------------------------------------------
void NodeBlock::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{

  int size, length;
  double val;
  int i, j;
  std::list<tagged_param>::iterator it_tpL;
  std::list<int>::iterator it_iL;
  std::string tmpStr;

  //----- unpack id_
  comm->unpack( pB, bsize, pos, &length, 1 );
  id_ = std::string( (pB+pos), length );
  pos += length;

  //----- unpack nodeList_
  comm->unpack( pB, bsize, pos, &size, 1 );
  nodeList_.clear();
  for( i = 0 ; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    tmpStr = std::string( (pB+pos), length );
    pos += length;
    comm->unpack( pB, bsize, pos, &val, 1 );
    nodeList_.push_back( tagged_param( tmpStr, val ) );
  }

  //----- unpack nodeProcList_
  comm->unpack( pB, bsize, pos, &size, 1 );
  nodeProcList_.clear();
  for( i = 0 ; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    tmpStr = std::string( (pB+pos), length );
    pos += length;
    comm->unpack( pB, bsize, pos, &val, 1 );
    nodeProcList_.push_back( tagged_param( tmpStr, val ) );
  }

  //----- unpack isOwned_
  comm->unpack( pB, bsize, pos, &i, 1 );
  isOwned_ = ( i != 0 );

  //----- unpack gID_
  comm->unpack( pB, bsize, pos, &gID_, 1 );

  //----- unpack procNum__
  comm->unpack( pB, bsize, pos, &procNum_, 1 );

  //----- unpack solnVarGIDList_
  comm->unpack( pB, bsize, pos, &size, 1 );
  solnVarGIDList_.clear();
  for( i = 0 ; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &j, 1 );
    solnVarGIDList_.push_back( j );
  }

  //----- unpack stateVarGIDList_
  comm->unpack( pB, bsize, pos, &size, 1 );
  stateVarGIDList_.clear();
  for( i = 0 ; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &j, 1 );
    stateVarGIDList_.push_back( j );
  }

  //----- unpack storeVarGIDList_
  comm->unpack( pB, bsize, pos, &size, 1 );
  storeVarGIDList_.clear();
  for( i = 0 ; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &j, 1 );
    storeVarGIDList_.push_back( j );
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Unpacked " << pos << " bytes for NodeBlock: " << id_ << std::endl;
#endif
}

} // namespace Topo
} // namespace Xyce
