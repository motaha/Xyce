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
// Filename       : $RCSfile: N_TOP_NodeDevBlock.h,v $
//
// Purpose        : Block of Node and Dev data for passing between pkgs
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
// Revision Number: $Revision: 1.18 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_NodeDevBlock_h
#define N_TOP_NodeDevBlock_h 1

#include <iosfwd>

#include <N_UTL_Packable.h>
#include <N_TOP_NodeBlock.h>
#include <N_DEV_DeviceBlock.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : NodeDevBlock
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
class NodeDevBlock : public Packable
{

public:

  // Constructor
  NodeDevBlock( NodeBlock nb = NodeBlock(),
  	              Device::InstanceBlock ib = Device::InstanceBlock() )
  : id_(nb.get_id()),
    gID_(nb.get_gID()),
    nodeBlock_(nb),
    devBlock_(ib)
  {}

  // Copy constructor.
  NodeDevBlock(const NodeDevBlock & right)
  : id_(right.id_),
    gID_(right.gID_),
    nodeBlock_(right.nodeBlock_),
    devBlock_(right.devBlock_)
  {}

  // Destructor
  ~NodeDevBlock()
  {}

  // Assignment operator (autogen)
  //NodeDevBlock & operator = (const NodeDevBlock & right);

  // claer data for reuse of this object
  void clear();

  // Equality operator
  bool operator==( const NodeDevBlock & right ) const
  { return id_ == right.id_; }

  // Non-equality operator
  bool operator!=( const NodeDevBlock & right ) const
  { return id_ != right.id_; }

  // Accessors for private attributes.

  // Get the node ID.
  const std::string & getID() const { return id_; }
  // Set the node ID.
  void setID( const std::string & id ) { id_ = id; }

  // Get the node global ID.
  const int & getGID() const { return gID_; }
  // Set the node global ID.
  void setGID( const int & gid ) { gID_ = gid; }

  NodeBlock & getNodeBlock() { return nodeBlock_; }
  const NodeBlock & getNodeBlock() const { return nodeBlock_; }
  Device::InstanceBlock & getDevBlock() { return devBlock_; }
  const Device::InstanceBlock & getDevBlock() const { return devBlock_; }

  bool isDevice() const { return devBlock_.getName() != ""; }

  // For packing and unpacking using MPI.
  Packable * instance() const;

  // Counts number of bytes needed to pack object.
  int packedByteCount() const;

  // Packs NodeBlock into char buffer using MPI_PACK.
  void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;
  // Unpacks NodeBlock from char buffer using MPI_UNPACK.
  void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);

protected:

  std::string id_;

  // Global ID.
  int gID_;

  NodeBlock nodeBlock_;
  Device::InstanceBlock devBlock_;

    friend std::ostream & operator << (std::ostream & os, const NodeDevBlock & ndb);

};

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::NodeDevBlock N_TOP_NodeDevBlock;

#endif
