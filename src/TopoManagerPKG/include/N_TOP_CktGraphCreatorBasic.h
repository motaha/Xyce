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
// Filename       : $RCSfile: N_TOP_CktGraphCreatorBasic.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/10/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktGraphCreatorBasic_h
#define N_TOP_CktGraphCreatorBasic_h 1

#include <N_TOP_CktGraphCreator.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktGraphCreatorBasic
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
class CktGraphCreatorBasic : public CktGraphCreator
{

public:

  // Default constructor.
  CktGraphCreatorBasic(const int maxTries) { maxTries_ = maxTries; }

  // Copy constructor.
  CktGraphCreatorBasic(const CktGraphCreatorBasic & right) { maxTries_ = right.maxTries_; }

  // Destructor
  ~CktGraphCreatorBasic() { }

  // Method to create a new 'Basic' based circuit graph.
  CktGraph * create(const std::string & cgID);

  // Method to create a new 'Basic' based circuit graph.
  CktGraph * create(const std::string & cgID,
                          const std::list<NodeID> & nodeList);

private:
  // Max number of attempts to compute the graph center.
  int maxTries_;
};

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::CktGraphCreatorBasic N_TOP_CktGraphCreatorBasic;

#endif
