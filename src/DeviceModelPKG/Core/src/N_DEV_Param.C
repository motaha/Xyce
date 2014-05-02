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
// Filename      : $RCSfile: N_DEV_Param.C,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Robert Hoekstra, SNL
//
// Creation Date : 5/15/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.30 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_Param.h>
#include <N_PDS_Comm.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : packedByteCount
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
int Param::packedByteCount() const
{
  // Util::Param info
  int byteCount = Util::Param::packedByteCount();

  // given & default
  byteCount += sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void Param::pack( char * buf, const int bsize, int & pos,
		N_PDS_Comm * comm ) const
{
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  Util::Param::pack( buf, bsize, pos, comm );

  //pack given_
  int dg = (isGiven_?1:0) + 2*(isDefault_?1:0);
  comm->pack( &dg, 1, buf, bsize, pos );

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    DevelFatal(*this, "Param::pack") << "Predicted pos does not match actual pos";
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void Param::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{
  Util::Param::unpack( pB, bsize, pos, comm );

  //unpack given_
  int dg;
  comm->unpack( pB, bsize, pos, &dg, 1 );
  isGiven_ = ( dg%2 != 0 );
  isDefault_ = ( dg >= 2 );
}

} // namespace Device
} // namespace Xyce
