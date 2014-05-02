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
// Filename       : $RCSfile: N_TOP_CktNode_Ckt.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/20/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_TOP_CktNode_Ckt.h>
#include <N_TOP_CktGraph.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : CktNode_Ckt::put
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
std::ostream& CktNode_Ckt::put(std::ostream& os) const
{
  os << "CN_Ckt: " << id_ << std::endl;
  os << "   GID: " << gID_ << "  Proc: " << procNum_ << std::endl;
  os << "   Soln Var GID List: ";
  std::list<int>::const_iterator it_iL = solnVarGIDList_.begin();
  std::list<int>::const_iterator it_iL_end = solnVarGIDList_.end();
  for( ; it_iL != it_iL_end; ++it_iL )
  {
    os << *it_iL << "  ";
  }
  os << std::endl;

  os << *cktGphPtr_ << std::endl;

  return os;
}

} // namespace Topo
} // namespace Xyce
