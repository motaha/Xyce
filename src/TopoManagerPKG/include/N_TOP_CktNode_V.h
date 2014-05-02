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
// Filename       : $RCSfile: N_TOP_CktNode_V.h,v $
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
// Revision Number: $Revision: 1.21 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktNode_V_h
#define N_TOP_CktNode_V_h 1

#include <N_TOP_CktNode.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktNode_V
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktNode_V : public CktNode
{

public:

  // Constructor
  CktNode_V(const std::string & nodeID = std::string(""), const int & globalID = 0,
                  const std::list< int > & varGIDList = std::list< int > (),
                  const std::list< int > & statevarGIDList = std::list< int > (),
                  const std::list< int > & storevarGIDList = std::list< int > (),
                  const int & pNum = 0, const bool & owned = true)
    : CktNode(nodeID, globalID, varGIDList, statevarGIDList, storevarGIDList, pNum, owned),
    noDCPath_(false),
    connToOneTerm_(false)
    {}

  // Constructor
  CktNode_V(const NodeBlock & nb)
    : CktNode(nb),
    noDCPath_(false),
    connToOneTerm_(false)
    {}

  // Copy constructor
  CktNode_V(const CktNode_V & right)
  : CktNode(right.id_, right.gID_, right.solnVarGIDList_,
                  right.stateVarGIDList_, 
                  right.storeVarGIDList_, 
                  right.procNum_, right.isOwned_),
    noDCPath_(right.noDCPath_),
    connToOneTerm_(right.connToOneTerm_)
  { }

  // Destructor
  ~CktNode_V() { }

  int type() const { return _VNODE; }

  int solnVarCount() { return 1; }

  bool getNoDCPathVar() {return noDCPath_;} //See explanation in "private" 
                                            //section
  bool getConnToOneTermVar() {return connToOneTerm_;}

  void setTrueNoDCPathVar() { noDCPath_ = true;}
  
  void setTrueConnToOneTermVar() { connToOneTerm_ = true;}

  std::ostream & put(std::ostream & os) const;

private:
  //These new boolean variables are being introduced to detect nodes for 
  //which there is no DC path to ground or nodes which are only connected to
  //one device terminal.  Detection of these situations is being used to 
  //(optionally) add large resistors between such nodes and ground  

  bool noDCPath_;
  bool connToOneTerm_;

};

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::CktNode_V N_TOP_CktNode_V;

#endif
