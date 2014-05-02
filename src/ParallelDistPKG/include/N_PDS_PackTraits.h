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
// Filename       : $RCSfile: N_PDS_PackTraits.h,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/11/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_PackTraits_h
#define Xyce_N_PDS_PackTraits_h

// ---------- Standard Includes ----------

#include <string>
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_PDS_Comm.h>
#include <N_TOP_Misc.h>

// ----------  Other Includes   ----------

// ----------  Fwd Declarations ----------

namespace Xyce {
namespace Parallel {

typedef N_PDS_Comm Comm;

//-----------------------------------------------------------------------------
// Class         : Xyce::Parallel::PackTraits
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 08/11/03
//-----------------------------------------------------------------------------
template <typename T>
struct PackTraits
{
  static int size( T const & object )
  { return const_cast<T&>(object).packedByteCount(); }

  static void pack( T const & object, char * buf, int size, int & pos, Comm & comm )
  { const_cast<T&>(object).pack( buf, size, pos, comm ); }

  static void unpack( T & object, char * buf, int size, int & pos, Comm & comm )
  { object.unpack( buf, size, pos, comm ); }
};

template <>
struct PackTraits<std::string>
{
  static int size( std::string const & object )
  { return object.length() + sizeof(int); }

  static void pack( std::string const & object, char * buf, int size, int & pos, Comm & comm )
  {
    int len = object.length();
    comm.pack( &len, 1, buf, size, pos );
    comm.pack( object.c_str(), len, buf, size, pos );
  }

  static void unpack( std::string & object, char * buf, int size, int & pos, Comm & comm )
  {
    int len;
    comm.unpack( buf, size, pos, &len, 1 );
    object = std::string( buf+pos, len );
    pos += len;
  }
};

template <>
struct PackTraits< std::vector<int> >
{
  static int size( std::vector<int> const & object )
  { return (object.size()+1) * sizeof(int); }

  static void pack( std::vector<int> const & object, char * buf, int size, int & pos, Comm & comm )
  {
    int len = object.size();
    comm.pack( &len, 1, buf, size, pos );
    for( int i = 0; i < len; ++i )
      comm.pack( &object[i], 1, buf, size, pos );
  }

  static void unpack( std::vector<int> & object, char * buf, int size, int & pos, Comm & comm )
  {
    int len;
    comm.unpack( buf, size, pos, &len, 1 );
    object.resize(len);
    for( int i = 0; i < len; ++i )
      comm.unpack( buf, size, pos, &object[i], 1 );
  }
};

template <>
struct PackTraits< NodeID >
{
  static int size( NodeID const & object )
  { return (object.first.length() + 2*sizeof(int)); }

  static void pack( NodeID const & object, char * buf, int size, int & pos, Comm & comm )
  {
    int len = object.first.length();
    comm.pack( &len, 1, buf, size, pos );
    comm.pack( object.first.c_str(), len, buf, size, pos );
    comm.pack( &object.second, 1, buf, size, pos );
  }

  static void unpack( NodeID & object, char * buf, int size, int & pos, Comm & comm )
  {
    int len;
    comm.unpack( buf, size, pos, &len, 1 );
    object.first = std::string( buf+pos, len );
    pos += len;
    comm.unpack( buf, size, pos, &object.second, 1 );
  }
};

} //namespace Parallel
} //namespace Xyce

#endif
