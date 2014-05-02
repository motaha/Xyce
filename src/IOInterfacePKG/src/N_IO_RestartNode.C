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
// Filename       : $RCSfile: N_IO_RestartNode.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 7/19/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.30 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <iterator>

#include <N_IO_RestartNode.h>

#include <N_DEV_DeviceState.h>

#include <N_PDS_Comm.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : RestartNode::RestartNode
// Purpose       : Copy
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
RestartNode::RestartNode(const RestartNode & right)
: ID(right.ID),
  type(right.type),
  solnVarData(right.solnVarData),
  stateVarData(right.stateVarData),
  storeVarData(right.storeVarData),
  devState(0)
{
  if( right.devState ) devState = new Device::DeviceState( *(right.devState) );
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::operator=
// Purpose       : Assign
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
RestartNode & RestartNode::operator=(const RestartNode & right)
{
  ID = right.ID;
  type = right.type;

  solnVarData = right.solnVarData;
  stateVarData = right.stateVarData;
  storeVarData = right.storeVarData;
  if( right.devState ) devState = new Device::DeviceState( *(right.devState) );

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::~RestartNode
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
RestartNode::~RestartNode()
{
  if( devState != 0 ) delete devState;
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::packedByteCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
int RestartNode::packedByteCount() const
{
  int byteCount = ID.length() + 2*sizeof(int); // ID + length + type

  int svdSize = solnVarData.size();

  byteCount += sizeof(int) * (1 + svdSize);
  for (int i = 0; i < svdSize; ++i)
    byteCount += solnVarData[i].size() * sizeof(double);

  svdSize = stateVarData.size();
  byteCount += sizeof(int) * (1 + svdSize);
  for (int i = 0; i < svdSize; ++i)
    byteCount += stateVarData[i].size() * sizeof(double);

  svdSize = storeVarData.size();
  byteCount += sizeof(int) * (1 + svdSize);
  for (int i = 0; i < svdSize; ++i)
    byteCount += storeVarData[i].size() * sizeof(double);

  byteCount += sizeof(int); //devState flag
  if (devState) byteCount += devState->packedByteCount();

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
void RestartNode::pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const
{
  int size, size2;
  int predictedPos = pos+packedByteCount();

  //pack tag
  size = ID.length();
  comm->pack( &size, 1, buf, bsize, pos );
  comm->pack( ID.c_str(), size, buf, bsize, pos );
  comm->pack( &type, 1, buf, bsize, pos );

  size = solnVarData.size();
  comm->pack( &size, 1, buf, bsize, pos );

  for( int i = 0; i < size; ++i )
  {
    size2 = solnVarData[i].size();
    comm->pack( &size2, 1, buf, bsize, pos );
    for( int ii = 0; ii < size2; ++ii )
      comm->pack( &(solnVarData[i][ii]), 1, buf, bsize, pos );
  }

  size = stateVarData.size();
  comm->pack( &size, 1, buf, bsize, pos );

  for( int i = 0; i < size; ++i )
  {
    size2 = stateVarData[i].size();
    comm->pack( &size2, 1, buf, bsize, pos );
    for( int ii = 0; ii < size2; ++ii )
      comm->pack( &(stateVarData[i][ii]), 1, buf, bsize, pos );
  }

  size = storeVarData.size();
  comm->pack( &size, 1, buf, bsize, pos );

  for( int i = 0; i < size; ++i )
  {
    size2 = storeVarData[i].size();
    comm->pack( &size2, 1, buf, bsize, pos );
    for( int ii = 0; ii < size2; ++ii )
      comm->pack( &(storeVarData[i][ii]), 1, buf, bsize, pos );
  }

  int flag = devState?1:0;
  comm->pack( &flag, 1, buf, bsize, pos );
  if( devState ) devState->pack( buf, bsize, pos, comm );
  if (pos != predictedPos)
  {
    std::string msg("Predicted pos does not match actual pos in RestartNode::pack");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
  }
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
void RestartNode::unpack(char * buf, int bsize, int & pos,
  N_PDS_Comm * comm)
{
  int size, size2;
  comm->unpack( buf, bsize, pos, &size, 1 );
  ID = std::string( (buf+pos), size);
  pos += size;
  comm->unpack( buf, bsize, pos, &type, 1 );

  comm->unpack( buf, bsize, pos, &size, 1 );
  solnVarData.resize(size);
  for( int i = 0; i < size; ++i )
  {
    comm->unpack( buf, bsize, pos, &size2, 1 );
    solnVarData[i].resize(size2);
    for( int ii = 0; ii < size2; ++ii )
      comm->unpack( buf, bsize, pos, &(solnVarData[i][ii]), 1 );
  }

  comm->unpack( buf,  bsize, pos, &size, 1 );
  stateVarData.resize(size);
  for( int i = 0; i < size; ++i )
  {
    comm->unpack( buf, bsize, pos, &size2, 1 );
    stateVarData[i].resize(size2);
    for( int ii = 0; ii < size2; ++ii )
      comm->unpack( buf, bsize, pos, &(stateVarData[i][ii]), 1 );
  }

  comm->unpack( buf,  bsize, pos, &size, 1 );
  storeVarData.resize(size);
  for( int i = 0; i < size; ++i )
  {
    comm->unpack( buf, bsize, pos, &size2, 1 );
    storeVarData[i].resize(size2);
    for( int ii = 0; ii < size2; ++ii )
      comm->unpack( buf, bsize, pos, &(storeVarData[i][ii]), 1 );
  }

  int flag;
  comm->unpack( buf, bsize, pos, &flag, 1 );
  if( flag )
  {
    devState = new Device::DeviceState();
    devState->unpack( buf, bsize, pos, comm );
  }
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const RestartNode & rn)
{
  os << Xyce::subsection_divider << std::endl
     << "Restart Node: " << rn.ID << " ( type = " << rn.type << " )" << std::endl;
  std::ostream_iterator<double> out( os, " " );
  if( !rn.solnVarData.empty() )
  {
    os << " SolnVarData: " << std::endl;
    for( unsigned int i = 0; i < rn.solnVarData.size(); ++i )
    {
      os << " " << i << " ";
      copy( rn.solnVarData[i].begin(), rn.solnVarData[i].end(), out );
      os << std::endl;
    }
    os << std::endl;
  }
  if( !rn.stateVarData.empty() )
  {
    os << " StateVarData: " << std::endl;
    for( unsigned int i = 0; i < rn.stateVarData.size(); ++i )
    {
      os << " " << i << " ";
      copy( rn.stateVarData[i].begin(), rn.stateVarData[i].end(), out );
      os << std::endl;
    }
    os << std::endl;
  }

  if( !rn.storeVarData.empty() )
  {
    os << " StoreVarData: " << std::endl;
    for( unsigned int i = 0; i < rn.storeVarData.size(); ++i )
    {
      os << " " << i << " ";
      copy( rn.storeVarData[i].begin(), rn.storeVarData[i].end(), out );
      os << std::endl;
    }
    os << std::endl;
  }

  if( rn.devState ) os << *(rn.devState) << std::endl;

  os << Xyce::subsection_divider << std::endl;
  return os;
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::dump
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 7/23/03
//-----------------------------------------------------------------------------
void RestartNode::dump( std::ostream & os )
{
  os << ID << " ";
  os << type << " ";

  int dsize = solnVarData.size();
  os << dsize << " ";
  for( int i = 0; i < dsize; ++i )
  {
    int nsize = solnVarData[i].size();
    os << nsize << " ";
    for( int j = 0; j < nsize; ++j )
      os << solnVarData[i][j] << " ";
  }

  dsize = stateVarData.size();
  os << dsize << " ";
  for( int i = 0; i < dsize; ++i )
  {
    int nsize = stateVarData[i].size();
    os << nsize << " ";
    for( int j = 0; j < nsize; ++j )
      os << stateVarData[i][j] << " ";
  }

  dsize = storeVarData.size();
  os << dsize << " ";
  for( int i = 0; i < dsize; ++i )
  {
    int nsize = storeVarData[i].size();
    os << nsize << " ";
    for( int j = 0; j < nsize; ++j )
      os << storeVarData[i][j] << " ";
  }

  int flag = 1;
  if( devState )
  {
    os << flag << " ";
    devState->dump( os );
  }
  else
  {
    flag = 0;
    os << flag << " ";
  }
}

//-----------------------------------------------------------------------------
// Function      : RestartNode::restore
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 7/23/03
//-----------------------------------------------------------------------------
void RestartNode::restore( std::istream & is )
{
  is >> ID;
  is >> type;

  int dsize;
  is >> dsize;
  solnVarData.resize(dsize);
  for( int i = 0; i < dsize; ++i )
  {
    int nsize;
    is >> nsize;
    solnVarData[i].resize(nsize);
    for( int j = 0; j < nsize; ++j )
      is >> solnVarData[i][j];
  }

  is >> dsize;
  stateVarData.resize(dsize);
  for( int i = 0; i < dsize; ++i )
  {
    int nsize;
    is >> nsize;
    stateVarData[i].resize(nsize);
    for( int j = 0; j < nsize; ++j )
      is >> stateVarData[i][j];
  }

  is >> dsize;
  storeVarData.resize(dsize);
  for( int i = 0; i < dsize; ++i )
  {
    int nsize;
    is >> nsize;
    storeVarData[i].resize(nsize);
    for( int j = 0; j < nsize; ++j )
      is >> storeVarData[i][j];
  }

  int flag;
  is >> flag;
  if( flag == 1 )
  {
    devState = new Device::DeviceState();
    devState->restore(is);
  }
}

} // namespace IO
} // namespace Xyce
