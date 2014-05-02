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
// Filename       : $RCSfile: N_TOP_NodeBlock.h,v $
//
// Purpose        : Block of Node data for passing between pkgs
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/18/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.23 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_NodeBlock_h
#define N_TOP_NodeBlock_h 1

#include <string>
#include <iosfwd>
#include <list>

#include <N_UTL_Misc.h>
#include <N_UTL_Packable.h>

namespace Xyce {
namespace Topo {

class NodeBlock : public Packable
{

public:

  // Constructor
  NodeBlock(const std::string & nodeID = std::string(),
  	const std::list< tagged_param > & nList = std::list< tagged_param > (),
  	const std::list< tagged_param > & npList = std::list< tagged_param > (),
                  const bool & own = false, const int & globalID = 0,
                  const int & pNum = 0,
                  const std::list< int > & varGIDList = std::list< int > (),
                  const std::list< int > & statevarGIDList = std::list< int > (),
                  const std::list< int > & storevarGIDList = std::list< int > ()
                  );

  // Destructor
  ~NodeBlock();

  void clear();

  // Equality operator
  bool operator == (const NodeBlock & right) const;

  // Non-equality operator
  bool operator != (const NodeBlock & right) const;

  //------- accessors for private attributes

  // Get id
  const std::string & get_id() const;
  // Set id
  void set_id(const std::string & nodeID);

  void addNode(const tagged_param &);
  const std::list< tagged_param > & get_NodeList() const;
  void set_NodeList(const std::list< tagged_param > & nList);

  const std::list< tagged_param > & get_NodeProcList() const;
  void set_NodeProcList(const std::list< tagged_param > & npList);

  const bool & get_IsOwned() const;
  void set_IsOwned(const bool & own);

  // Get the global ID.
  const int & get_gID() const;
  // Set the global ID.
  void set_gID(const int & globalID);

  // Get the processor number.
  const int & get_ProcNum() const;
  // Set the processor number.
  void set_ProcNum(const int & pNum);

  const std::list< int > & get_SolnVarGIDList() const;
  void set_SolnVarGIDList(const std::list< int > & svGList);

  const std::list< int > & get_StateVarGIDList() const;
  void set_StateVarGIDList(const std::list< int > & svGList);

  const std::list< int > & get_StoreVarGIDList() const;
  void set_StoreVarGIDList(const std::list< int > & svGList);

  // For packing and unpacking using MPI.
  Packable * instance() const;

  // Counts number of bytes needed to pack object.
  int packedByteCount() const;

  // Packs NodeBlock into char buffer using MPI_PACK.
  void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;
  // Unpacks NodeBlock from char buffer using MPI_UNPACK.
  void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);

protected:

  std::string id_;

  // tag = adj node name; param = adj node global index.
  std::list< tagged_param > nodeList_;

  // tag = adj node name; param = adj node owning processor.
  std::list< tagged_param > nodeProcList_;

  // Owned by this processor, otherwise a ghosted node.
  bool isOwned_;

  // Global index.
  int gID_;

  // Owning processor.
  int procNum_;

  // List of associated solution variables.
  std::list< int > solnVarGIDList_;

  // List of associated state variables.
  std::list< int > stateVarGIDList_;

  // List of associated store variables.
  std::list< int > storeVarGIDList_;

    friend std::ostream & operator << (std::ostream & os, const NodeBlock & nb);

};

//-----------------------------------------------------------------------------
// Function      : NodeBlock::NodeBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline NodeBlock::NodeBlock(const std::string & nodeID,
                                        const std::list< tagged_param > & nList,
	const std::list< tagged_param > & npList, const bool & own,
                                        const int & globalID, const int & pNum,
                                        const std::list< int > & varGIDList,
                                        const std::list< int > & statevarGIDList,
                                        const std::list< int > & storevarGIDList
                                        )
  :
  id_(nodeID), isOwned_(own), gID_(globalID), procNum_(pNum)
{
  solnVarGIDList_.assign(varGIDList.begin(), varGIDList.end());
  stateVarGIDList_.assign(statevarGIDList.begin(), statevarGIDList.end());
  storeVarGIDList_.assign(storevarGIDList.begin(), storevarGIDList.end());
  nodeList_.assign(nList.begin(), nList.end());
  nodeProcList_.assign(npList.begin(), npList.end());
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::~NodeBlock
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline NodeBlock::~NodeBlock()
{
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::clear
// Purpose       : clears data from block
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
inline void NodeBlock::clear()
{
  id_ = "";
  gID_ = 0;
  procNum_ = 0;
  isOwned_ = false;

  solnVarGIDList_.erase(solnVarGIDList_.begin(), solnVarGIDList_.end());
  stateVarGIDList_.erase(stateVarGIDList_.begin(), stateVarGIDList_.end());
  storeVarGIDList_.erase(storeVarGIDList_.begin(), storeVarGIDList_.end());
  nodeList_.erase(nodeList_.begin(), nodeList_.end());
  nodeProcList_.erase(nodeProcList_.begin(), nodeProcList_.end());
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::operator==
// Purpose       : equality operator
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline bool NodeBlock::operator == (const NodeBlock & right) const
{
  return id_ == right.id_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::operator!=
// Purpose       : inequality operator
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline bool NodeBlock::operator != (const NodeBlock & right) const
{
  return id_ != right.id_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_id
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline const std::string & NodeBlock::get_id() const
{
  return id_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_id
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline void NodeBlock::set_id(const std::string & nodeID)
{
  id_ = nodeID;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::addNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 03/16/06
//-----------------------------------------------------------------------------
inline void NodeBlock::addNode(const tagged_param & p)
{
  nodeList_.push_back(p);
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_NodeList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline const std::list< tagged_param > & NodeBlock::get_NodeList() const
{
  return nodeList_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_NodeList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline void NodeBlock::set_NodeList(const std::list< tagged_param > & nList)
{
  nodeList_.assign(nList.begin(), nList.end());
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_NodeProcList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline const std::list< tagged_param > & NodeBlock::get_NodeProcList() const
{
  return nodeProcList_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_NodeProcList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline void NodeBlock::set_NodeProcList(
  const std::list< tagged_param > & npList)
{
  nodeProcList_.assign(npList.begin(), npList.end());
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_IsOwned
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline const bool & NodeBlock::get_IsOwned() const
{
  return isOwned_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_IsOwned
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline void NodeBlock::set_IsOwned(const bool & own)
{
  isOwned_ = own;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_gID
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline const int & NodeBlock::get_gID() const
{
  return gID_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_gID
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline void NodeBlock::set_gID(const int & globalID)
{
  gID_ = globalID;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_ProcNum
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline const int & NodeBlock::get_ProcNum() const
{
  return procNum_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_ProcNum
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline void NodeBlock::set_ProcNum(const int & pNum)
{
  procNum_ = pNum;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_SolnVarGIDList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline const std::list< int > & NodeBlock::get_SolnVarGIDList() const
{
  return solnVarGIDList_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_SolnVarGIDList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
inline void NodeBlock::set_SolnVarGIDList(const std::list< int > & svGList)
{
  solnVarGIDList_.assign(svGList.begin(), svGList.end());
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_StateVarGIDList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
inline const std::list< int > & NodeBlock::get_StateVarGIDList() const
{
  return stateVarGIDList_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_StateVarGIDList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
inline void NodeBlock::set_StateVarGIDList(const std::list< int > & svGList)
{
  stateVarGIDList_.assign(svGList.begin(), svGList.end());
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::get_StoreVarGIDList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
inline const std::list< int > & NodeBlock::get_StoreVarGIDList() const
{
  return storeVarGIDList_;
}

//-----------------------------------------------------------------------------
// Function      : NodeBlock::set_StoreVarGIDList
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
inline void NodeBlock::set_StoreVarGIDList(const std::list< int > & svGList)
{
  storeVarGIDList_.assign(svGList.begin(), svGList.end());
}

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::NodeBlock N_TOP_NodeBlock;

#endif

