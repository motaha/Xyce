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
// Filename       : $RCSfile: N_TOP_CktNode_Ckt.h,v $
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
// Revision Number: $Revision: 1.19 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktNode_Ckt_h
#define N_TOP_CktNode_Ckt_h 1

#include <N_TOP_CktNode.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktNode_Ckt
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktNode_Ckt : public CktNode
{

private:

  // Default constructor (private).
  CktNode_Ckt() { }

  // Assignment operator (private).
  CktNode_Ckt & operator = (const CktNode_Ckt & right);

public:

  // Constructor
  CktNode_Ckt(CktGraph * grphPtr, const std::string & ID = std::string(""),
                    const int & gID = 0,
                    const std::list< int > & varGIDList = std::list< int > (),
                    const std::list< int > & statevarGIDList = std::list< int > (),
                    const std::list< int > & storevarGIDList = std::list< int > (),
                    const int & pNum = 0, const bool & owned = true)
    :
    CktNode(ID, gID, varGIDList, statevarGIDList, storevarGIDList, pNum, owned),
    cktGphPtr_(grphPtr)
    { }

  // Constructor
  CktNode_Ckt(CktGraph * grphPtr, const NodeBlock & nb)
    :
    CktNode(nb), cktGphPtr_(grphPtr)
    { }

  // Destructor
  ~CktNode_Ckt() { }

  int type() const { return _CNODE; }

  bool getNoDCPathVar() {return false;} 

  bool getConnToOneTermVar() {return false;}

  void setTrueNoDCPathVar() {}
  
  void setTrueConnToOneTermVar() {}

private:

  // Pointer to the circuit graph.
  CktGraph * cktGphPtr_;

public:

  std::ostream & put(std::ostream & os) const;

};

//-----------------------------------------------------------------------------
// Function      : CktNode_Ckt::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
inline CktNode_Ckt & CktNode_Ckt::operator = (
  const CktNode_Ckt & right)
{
  return * this;
}

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::CktNode_Ckt N_TOP_CktNode_Ckt;

#endif
