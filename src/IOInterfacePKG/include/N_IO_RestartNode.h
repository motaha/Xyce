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
// Filename       : $RCSfile: N_IO_RestartNode.h,v $
//
// Purpose        : Node storing restart info associated with an ID
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 8/22/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.20 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_RestartNode_h
#define Xyce_N_IO_RestartNode_h

#include <vector>
#include <iosfwd>

#include <N_DEV_fwd.h>
#include <N_UTL_Packable.h>
#include <N_TOP_Misc.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : RestartNode
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
class RestartNode : public Packable
{
public:

  // Constructor
  RestartNode(const std::string & id = "", const int inType = _VNODE ) : ID(id), type(inType), devState(0) { }

  // Destructor
  ~RestartNode();

  // Copy constructor
  RestartNode(const RestartNode & right);

  // Assignment operator
  RestartNode & operator = (const RestartNode & right);

  bool operator == (const RestartNode & right) { return ID == right.ID; }

  bool operator != (const RestartNode & right) { return ID != right.ID; }

  bool operator < (const RestartNode & right) { return ID < right.ID; }

  Packable * instance() const { return new RestartNode(); }

  int packedByteCount() const;

  void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;
  void unpack(char * buf, int bsize, int & pos, N_PDS_Comm * comm);

  void dump( std::ostream & os );
  void restore( std::istream & is );

  std::string ID;
  int type;

  std::vector< std::vector< double > > solnVarData;
  std::vector< std::vector< double > > stateVarData;
  std::vector< std::vector< double > > storeVarData;

  Device::DeviceState * devState;

  friend std::ostream & operator << (std::ostream & os, const RestartNode & rn);
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::RestartNode N_IO_RestartNode;

#endif
