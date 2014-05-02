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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DeviceState.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 09/02/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.15 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <iostream>

#include <N_DEV_DeviceState.h>

#include <N_PDS_Comm.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceState::packedByteCount
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoektra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
int DeviceState::packedByteCount() const
{
  int bCnt = sizeof(int);  //ID length
  bCnt += ID.length();

  bCnt += sizeof(int);  //data double length
  bCnt += data.size() * sizeof(double);

  bCnt += sizeof(int);  //data int length
  bCnt += dataInt.size() * sizeof(int);

  bCnt += sizeof(int);  //data size_t length
  bCnt += dataSizeT.size() * sizeof(size_t);

  return bCnt;
}

//-----------------------------------------------------------------------------
// Function      : DeviceState::pack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
void DeviceState::pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const
{
  int length;

  //----- pack ID
  length = ID.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( ID.c_str(), length, buf, bsize, pos );

  //----- pack double data
  length = data.size();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( &(data[0]), length, buf, bsize, pos );

  //----- pack int data
  length = dataInt.size();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( &(dataInt[0]), length, buf, bsize, pos );

  //----- pack size_t data
  length = dataSizeT.size();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( &(dataSizeT[0]), length, buf, bsize, pos );
}

//-----------------------------------------------------------------------------
// Function      : DeviceState::unpack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
void DeviceState::unpack(char * buf, int bsize, int & pos, N_PDS_Comm * comm)
{
  int length;

  //----- unpack ID
  comm->unpack( buf, bsize, pos, &length, 1 );
  ID = std::string( (buf+pos), length);
  pos += length;

  //----- unpack data
  comm->unpack( buf, bsize, pos, &length, 1 );
  data.resize(length);
  comm->unpack( buf, bsize, pos, &(data[0]), length );

  //----- unpack int data
  comm->unpack( buf, bsize, pos, &length, 1 );
  dataInt.resize(length);
  comm->unpack( buf, bsize, pos, &(dataInt[0]), length );

  //----- unpack size_t data
  comm->unpack( buf, bsize, pos, &length, 1 );
  dataSizeT.resize(length);
  comm->unpack( buf, bsize, pos, &(dataSizeT[0]), length );
}

//-----------------------------------------------------------------------------
// Function      : DeviceState::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
std::ostream & operator<<( std::ostream & os, const DeviceState & ds )
{
  os << "Device State: " << ds.ID << std::endl;
  os << " -------------" << std::endl;
  for( int i = 0; i < ds.data.size(); ++i )
    os << " " << i << ": " << ds.data[i] << std::endl;
  os << " -------------" << std::endl;
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : DeviceState::dump
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/03/04
//-----------------------------------------------------------------------------
void DeviceState::dump( std::ostream & os )
{
  os << ID << " ";

  int size = data.size();
  os << size << " ";
  for( int i = 0; i < size; ++i )
    os << data[i] << " ";


  size = dataInt.size();
  os << size << " ";
  for( int i = 0; i < size; ++i )
    os << dataInt[i] << " ";


  size = dataSizeT.size();
  os << size << " ";
  for( int i = 0; i < size; ++i )
    os << dataSizeT[i] << " ";
}

//-----------------------------------------------------------------------------
// Function      : DeviceState::restore
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/03/04
//-----------------------------------------------------------------------------
void DeviceState::restore( std::istream & is )
{
  is >> ID;

  int size;
  is >> size;
  data.resize(size);
  for( int i = 0; i < size; ++i )
    is >> data[i];

  is >> size;
  dataInt.resize(size);
  for( int i = 0; i < size; ++i )
    is >> dataInt[i];

  is >> size;
  dataSizeT.resize(size);
  for( int i = 0; i < size; ++i )
    is >> dataSizeT[i];
}

} // namespace Device
} // namespace Xyce
