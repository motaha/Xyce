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
// Revision Number: $Revision: 1.4.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:50 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktGraphCreatorBasic_h
#define N_TOP_CktGraphCreatorBasic_h 1

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_TOP_CktGraphCreator.h>

// ---------- Forward Declarations ----------

class N_TOP_CktGraph;

//-----------------------------------------------------------------------------
// Class         : N_DEV_CktGraphCreatorBasic
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
class N_TOP_CktGraphCreatorBasic : public N_TOP_CktGraphCreator
{

public:

  // Default constructor.
  N_TOP_CktGraphCreatorBasic(const int maxTries) { maxTries_ = maxTries; }

  // Copy constructor.
  N_TOP_CktGraphCreatorBasic(const N_TOP_CktGraphCreatorBasic & right) { maxTries_ = right.maxTries_; }

  // Destructor
  ~N_TOP_CktGraphCreatorBasic() { }

  // Method to create a new 'Basic' based circuit graph.
  N_TOP_CktGraph * create(const string & cgID);

  // Method to create a new 'Basic' based circuit graph.
  N_TOP_CktGraph * create(const string & cgID,
                          const list <NodeID> & nodeList);

private:
  // Max number of attempts to compute the graph center.
  int maxTries_;
};

#endif
