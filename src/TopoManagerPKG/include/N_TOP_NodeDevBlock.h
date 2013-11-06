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
// Revision Number: $Revision: 1.12.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_NodeDevBlock_h
#define N_TOP_NodeDevBlock_h 1

//----------- Std includes ---------
#include <N_UTL_Misc.h>
#include <iosfwd>

//----------- Xyce includes ---------


#include <N_UTL_Packable.h>

#include <N_TOP_NodeBlock.h>
#include <N_DEV_DeviceBlock.h>

//----------- Other includes ---------

//----------- Fwd Declares ----------

//-----------------------------------------------------------------------------
// Class         : N_TOP_NodeDevBlock
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/3/01
//-----------------------------------------------------------------------------
class N_TOP_NodeDevBlock : public Packable
{

public:

  // Constructor
  N_TOP_NodeDevBlock( N_TOP_NodeBlock nb = N_TOP_NodeBlock(),
  	              N_DEV_InstanceBlock ib = N_DEV_InstanceBlock() )
  : id_(nb.get_id()),
    gID_(nb.get_gID()),
    nodeBlock_(nb),
    devBlock_(ib)
  {}

  // Copy constructor.
  N_TOP_NodeDevBlock(const N_TOP_NodeDevBlock & right)
  : id_(right.id_),
    gID_(right.gID_),
    nodeBlock_(right.nodeBlock_),
    devBlock_(right.devBlock_)
  {}

  // Destructor
  ~N_TOP_NodeDevBlock()
  {}

  // Assignment operator (autogen)
  //N_TOP_NodeDevBlock & operator = (const N_TOP_NodeDevBlock & right);

  // claer data for reuse of this object
  void clear();

  // Equality operator
  bool operator==( const N_TOP_NodeDevBlock & right ) const
  { return id_ == right.id_; }

  // Non-equality operator
  bool operator!=( const N_TOP_NodeDevBlock & right ) const
  { return id_ != right.id_; }

  // Accessors for private attributes.

  // Get the node ID.
  const string & getID() const { return id_; }
  // Set the node ID.
  void setID( const string & id ) { id_ = id; }

  // Get the node global ID.
  const int & getGID() const { return gID_; }
  // Set the node global ID.
  void setGID( const int & gid ) { gID_ = gid; }

  N_TOP_NodeBlock & getNodeBlock() { return nodeBlock_; }
  N_DEV_InstanceBlock & getDevBlock() { return devBlock_; }

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

  string id_;

  // Global ID.
  int gID_;

  N_TOP_NodeBlock nodeBlock_;
  N_DEV_InstanceBlock devBlock_;

  friend ostream & operator << (ostream & os, const N_TOP_NodeDevBlock & ndb);

};

#endif
