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
// Filename       : $RCSfile: N_TOP_NodeDevBlock.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 3/3/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.19 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TOP_NodeDevBlock.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : NodeDevBlock::instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
Packable * NodeDevBlock::instance() const
{
  return new NodeDevBlock();
}


//-----------------------------------------------------------------------------
// Function      : NodeDevBlock::clear
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 2/8/2010
//-----------------------------------------------------------------------------
void NodeDevBlock::clear()
{
  id_.erase();
  gID_=0;
  devBlock_.clear();
  nodeBlock_.clear();
}


//-----------------------------------------------------------------------------
// Function      : NodeBlock::packedByteCount
// Purpose       : Counts number of bytes needed to pack object
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
int NodeDevBlock::packedByteCount() const
{
  int byteCount = 0;

  byteCount += nodeBlock_.packedByteCount();

  //flag for device
  byteCount += sizeof(int);

  if( isDevice() ) byteCount += devBlock_.packedByteCount();

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : NodeDevBlock::pack
// Purpose       : Packs NodeDevBlock into char buffer using MPI_PACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
void NodeDevBlock::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
{

  int flag;

  //----- pack nodeBlock_
  nodeBlock_.pack( buf, bsize, pos, comm );

  //----- pack nodeBlock_ flag
  flag = isDevice() ? 1 : 0;
  comm->pack( &flag, 1, buf, bsize, pos );

  //----- pack nodeBlock_
  if( isDevice() ) devBlock_.pack( buf, bsize, pos, comm );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Packed " << pos << " bytes for NodeDevBlock: " <<
	id_ << std::endl;
#endif

}

//-----------------------------------------------------------------------------
// Function      : NodeDevBlock::unpack
// Purpose       : Unpacks NodeDevBlock from char buffer using MPI_UNPACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
void NodeDevBlock::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{

  int flag;

  //----- unpack nodeBlock_
  nodeBlock_.unpack( pB, bsize, pos, comm );
  id_ = nodeBlock_.get_id();
  gID_ = nodeBlock_.get_gID();

  comm->unpack( pB, bsize, pos, &flag, 1 );
  if( flag == 1 ) devBlock_.unpack( pB, bsize, pos, comm );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Unpacked " << pos << " bytes for NodeDevBlock: " <<
	id_ << std::endl;
#endif
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
std::ostream & operator<< ( std::ostream & os, const NodeDevBlock & ndb )
{
  os << "NodeDevBlock: " << ndb.id_ << std::endl;
  os << ndb.nodeBlock_ << std::endl;
  if( ndb.isDevice() ) os << ndb.devBlock_ << std::endl;
  return os << std::endl;
}

} // namespace Topo
} // namespace Xyce
