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
// Revision Number: $Revision: 1.14.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_RestartNode_h
#define Xyce_N_IO_RestartNode_h

// ---------- Standard Includes ----------
#include <vector>
#include <iosfwd>

// ----------   Xyce Includes   ----------

#include <N_DEV_fwd.h>
#include <N_UTL_Packable.h>
#include <N_TOP_Misc.h>

// ---------- Forward Declarations ----------

//-----------------------------------------------------------------------------
// Class         : N_IO_RestartNode
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
class N_IO_RestartNode : public Packable
{
public:

  // Constructor
  N_IO_RestartNode(const string & id = "", const int inType = _VNODE ) : ID(id), type(inType), devState(0) { }

  // Destructor
  ~N_IO_RestartNode();

  // Copy constructor
  N_IO_RestartNode(const N_IO_RestartNode & right);

  // Assignment operator
  N_IO_RestartNode & operator = (const N_IO_RestartNode & right);

  bool operator == (const N_IO_RestartNode & right) { return ID == right.ID; }

  bool operator != (const N_IO_RestartNode & right) { return ID != right.ID; }

  bool operator < (const N_IO_RestartNode & right) { return ID < right.ID; }

  Packable * instance() const { return new N_IO_RestartNode(); }

  int packedByteCount() const;

  void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;
  void unpack(char * buf, int bsize, int & pos, N_PDS_Comm * comm);

  void dump( ostream & os );
  void restore( istream & is );

  string ID;
  int type;

  vector < vector < double > > solnVarData;
  vector < vector < double > > stateVarData;
  vector < vector < double > > storeVarData;

  N_DEV_DeviceState * devState;

  friend ostream & operator << (ostream & os, const N_IO_RestartNode & rn);
};

#endif
