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
// Filename       : $RCSfile: N_UTL_Packable.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  _PACKABLE_H
#define  _PACKABLE_H

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_PDS_fwd.h>

//-----------------------------------------------------------------------------
// Class         : Packable
// Purpose       : Abstract class for packable objects.
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/23/01
//-----------------------------------------------------------------------------
class Packable
{
public:
  virtual Packable * instance() const = 0;

  // Counts bytes needed to pack block (abstract).
  virtual int packedByteCount() const = 0;

  // Packs OptionBlock into char buffer using MPI_PACK (abstract).
  virtual void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const = 0;

  // Unpacks OptionBlock from char buffer using MPI_UNPACK (abstract).
  virtual void unpack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) = 0;
};

#endif
