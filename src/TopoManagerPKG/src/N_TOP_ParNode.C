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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_TOP_ParNode.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/17/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>



// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_TOP_ParNode.h>
#include <N_PDS_Comm.h>

//-----------------------------------------------------------------------------
// Function      : N_TOP_ParNode::pack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/17/01
//-----------------------------------------------------------------------------
void N_TOP_ParNode::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
{
  N_TOP_Node::pack( buf, bsize, pos, comm );

  comm->pack( &proc_, 1, buf, bsize, pos );

  int var = ghosts_.size();
  comm->pack( &var, 1, buf, bsize, pos );
  set<int>::const_iterator iterGS = ghosts_.begin();
  set<int>::const_iterator iterGS_end = ghosts_.end();
  for( ; iterGS != iterGS_end; ++iterGS )
  {
    var = *iterGS;
    comm->pack( &var, 1, buf, bsize, pos );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_ParNode::unpack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/17/01
//-----------------------------------------------------------------------------
void N_TOP_ParNode::unpack( char * buf, int bsize, int & pos, N_PDS_Comm * comm )
{
  N_TOP_Node::unpack( buf, bsize, pos, comm );

  comm->unpack( buf, bsize, pos, &proc_, 1 );

  ghosts_.clear();
  int size, var;
  comm->unpack( buf, bsize, pos, &size, 1 );
  for( int i = 0; i < size; ++i )
  {
    comm->unpack( buf, bsize, pos, &var, 1 );
    ghosts_.insert( var );
  }
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/17/01
//-----------------------------------------------------------------------------
ostream & operator<<( ostream & os, const N_TOP_ParNode & pn )
{
  pn.put(os);
  return os;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_ParNode::put
// Purpose       : Constructor
// Scope         : public
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/01
//-----------------------------------------------------------------------------
ostream & N_TOP_ParNode::put(ostream & os) const
{
  os << "-------------------------" << endl;
  os << "PARALLEL Node" << endl;
  N_TOP_Node::put(os);
  os << "Proc Owner:\t" << proc_ << endl;

  if (!ghosts_.empty())
  {
    os << "Ghosting Procs:";
    set < int >::const_iterator iterGS = ghosts_.begin();
    set < int >::const_iterator iterGS_end = ghosts_.end();
    for ( ; iterGS != iterGS_end; ++iterGS)
    {
      os << " " << * iterGS;
    }
    os << endl;
  }
  os << "-------------------------" << endl;
  os << endl;
  return os;
}
