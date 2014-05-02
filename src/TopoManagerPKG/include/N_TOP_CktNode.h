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
// Filename       : $RCSfile: N_TOP_CktNode.h,v $
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
// Revision Number: $Revision: 1.49 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktNode_h
#define N_TOP_CktNode_h 1

#include <iosfwd>
#include <string>
#include <list>
#include <vector>
#include <map>

#include <N_TOP_fwd.h>

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_TOP_Misc.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktNode
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktNode
{

public:

  // Constructor
  CktNode(const std::string & nodeID = std::string(""),
                const int & globalID = 0,
                const std::list< int > & varGIDList = std::list< int > (),
                const std::list< int > & statevarGIDList = std::list< int > (),
                const std::list< int > & storevarGIDList = std::list< int > (),
                const int & pNum = 0, const bool & owned = true,
                GraphNode * gN = 0 )
  : id_(nodeID),
    gID_(globalID),
    procNum_(pNum),
    isOwned_(owned),
    Offset_(-1),
    graphNodePtr_(gN),
    solnVarGIDList_(varGIDList),
    stateVarGIDList_(statevarGIDList),
    storeVarGIDList_(storevarGIDList)
  {}

  // Constructor
  CktNode( const NodeBlock & nb,
                 GraphNode * gN = 0 );

  // Auto generated copy construct.

  // Destructor
  virtual ~CktNode() { }

  bool operator == (const CktNode & right) const { return id_ == right.id_; }
  bool operator != (const CktNode & right) const { return id_ != right.id_; }

  //-------- get and set methods for attributes

  // Get circuit node id.
  const std::string & get_id() const { return id_; }
  // Set circuit node id.
  void set_id(const std::string & value) { id_ = value; }

  // Get circuit node global id.
  const int & get_gID() const { return gID_; }
  // Set circuit node global id.
  void set_gID(const int & globalID) { gID_ = globalID; }

  virtual int type() const = 0;

  //The following get/set functions are only used in CktNode_V class
  //to access and change the noDCPath_ and connToOneTermVar_ internal variables
  //These boolean variables are not defined in the other derived classes, so
  //we create dummy functions there that do nothing.

  virtual bool getNoDCPathVar() = 0;

  virtual bool getConnToOneTermVar() = 0;

  virtual void setTrueNoDCPathVar() = 0;

  virtual void setTrueConnToOneTermVar() = 0;



  const std::list<int> & get_SolnVarGIDList() const
  { return solnVarGIDList_; }
  void set_SolnVarGIDList(const std::list<int> & svGIDList)
  { solnVarGIDList_ = svGIDList; }
  void set_SolnVarGIDVec(const std::vector<int> & svGIDVec)
  { solnVarGIDList_.assign( svGIDVec.begin(), svGIDVec.end() ); }

  const std::list<int> & get_ExtSolnVarGIDList() const
  { return extSolnVarGIDList_; }
  void set_ExtSolnVarGIDList(const std::list<int> & svGIDList)
  { extSolnVarGIDList_ = svGIDList; }
  void set_ExtSolnVarGIDVec(const std::vector<int> & svGIDVec)
  { extSolnVarGIDList_.assign( svGIDVec.begin(), svGIDVec.end() ); }

  const std::list< int > & get_StateVarGIDList() const
  { return stateVarGIDList_; }
  void set_StateVarGIDList(const std::list< int > & svGIDList)
  { stateVarGIDList_ = svGIDList; }
  void set_StateVarGIDVec(const std::vector< int > & svGIDVec)
  { stateVarGIDList_.assign( svGIDVec.begin(), svGIDVec.end() ); }

  const std::list< int > & get_StoreVarGIDList() const
  { return storeVarGIDList_; }
  void set_StoreVarGIDList(const std::list< int > & svGIDList)
  { storeVarGIDList_ = svGIDList; }
  void set_StoreVarGIDVec(const std::vector< int > & svGIDVec)
  { storeVarGIDList_.assign( svGIDVec.begin(), svGIDVec.end() ); }

  const std::vector< std::vector<int> > & get_DepSolnGIDVec() const
  { return depSolnGIDVec_; }
  void set_DepSolnGIDVec(const std::vector< std::vector<int> > & dsGIDs)
  { depSolnGIDVec_ = dsGIDs; }
  const std::vector< std::vector<int> > & get_DepStateGIDVec() const
  { return depStateGIDVec_; }
  void set_DepStateGIDVec(const std::vector< std::vector<int> > & dsGIDs)
  { depStateGIDVec_ = dsGIDs; }
  const std::vector< std::vector<int> > & get_DepStoreGIDVec() const
  { return depStoreGIDVec_; }
  void set_DepStoreGIDVec(const std::vector< std::vector<int> > & dsGIDs)
  { depStoreGIDVec_ = dsGIDs; }

  // Get the processor number.
  const int & get_ProcNum() const { return procNum_; }
  // Set the processor number.
  void set_ProcNum(const int & pNum) { procNum_ = pNum; }

  // Get the processor ownership flag.
  const bool & get_IsOwned() const { return isOwned_; }
  // Set the processor ownership flag.
  void set_IsOwned(const bool & owned) { isOwned_ = owned; }

  const int & get_Offset() const { return Offset_; }
  void set_Offset(const int & Offset) { Offset_ = Offset; }

  GraphNode * get_GraphNode() { return graphNodePtr_; }
  void set_GraphNode(GraphNode * gN) { graphNodePtr_ = gN; }

  NodeBlock * extractNodeBlock();

  virtual const std::vector< std::vector<int> > & jacobianStamp() const
  { static std::vector< std::vector<int> > dummy; return dummy; }
  virtual void registerJacLIDswithDev( const std::vector< std::vector<int> > & jacLIDVec ) {}

  virtual void registerGIDDataWithDev(
               const std::vector<int> & counts,
               const std::vector<int> & GIDs,
               const std::vector< std::vector<int> > & jacGIDs ) {}

  //------- Added for use with the outputFileName function
  virtual std::map<int,std::string> & getIntNameMap() { return intNameMap_; }
  virtual std::map<int,std::string> & getStateNameMap() { return stateNameMap_; }
  virtual std::map<int,std::string> & getStoreNameMap() { return storeNameMap_; }

  virtual void varTypeList( std::vector<char> & varTypeVec ) {}

protected:

  // Node id.
  std::string id_;
  // Global node id.
  int gID_;

  // Processor number.
  int procNum_;
  // Processor ownership flag.
  bool isOwned_;

  int Offset_;

  GraphNode * graphNodePtr_;

  // Solution variable global id list.
  std::list<int> solnVarGIDList_;
  std::list<int> extSolnVarGIDList_;
  // State variable global id list.
  std::list<int> stateVarGIDList_;
  std::list<int> storeVarGIDList_;

  std::vector< std::vector<int> > depSolnGIDVec_;
  std::vector< std::vector<int> > depStateGIDVec_;
  std::vector< std::vector<int> > depStoreGIDVec_;

  std::vector<int> depSolnGIDJacVec_;

  std::map<int,std::string> intNameMap_;
  std::map<int,std::string> stateNameMap_;
  std::map<int,std::string> storeNameMap_;

public:

  //------- registration of global ids with device instances

  virtual int solnVarCount() { return solnVarGIDList_.size(); }
  virtual int extSolnVarCount() { return extSolnVarGIDList_.size(); }
  virtual int stateVarCount() { return stateVarGIDList_.size(); }
  virtual int depSolnVarCount() { return depSolnGIDVec_.size(); }
  virtual int depStateVarCount() { return depStateGIDVec_.size(); }

  virtual int storeVarCount() { return storeVarGIDList_.size(); }
  virtual int depStoreVarCount() { return depStoreGIDVec_.size(); }

  virtual void leadConnect(std::vector<int> &) {return;}

  virtual void registerGIDswithDev( const std::list<index_pair> & intGIDList,
  	                            const std::list<index_pair> & extGIDList) {}

  virtual void registerStateGIDswithDev( const std::list<index_pair> & stateGIDList) {}
  virtual void registerStoreGIDswithDev( const std::list<index_pair> & storeGIDList) {}

  virtual void registerLIDswithDev( const std::vector<int> & intLIDVec,
                                    const std::vector<int> & extLIDVec ) {}
  virtual void registerStateLIDswithDev( const std::vector<int> & stateLIDVec ) {}
  virtual void registerStoreLIDswithDev( const std::vector<int> & storeLIDVec ) {}

  virtual void registerDepLIDswithDev( const std::vector< std::vector<int> > & depLIDVec ) {}
  virtual void registerDepStateLIDswithDev( const std::vector< std::vector<int> > & depStateLIDVec ) {}
  virtual void registerDepStoreLIDswithDev( const std::vector< std::vector<int> > & depStoreLIDVec ) {}

  virtual void getDepSolnVars( std::vector< NodeID >& dsVars );
  virtual void registerDepSolnGIDs( const std::vector< std::vector<int> > & dsGIDs) {}
  virtual void getDepStateVars( std::vector< NodeID >& dsVars );
  virtual void registerDepStateGIDs(const std::vector< std::vector< int > > & dsGIDs) {}
  virtual void getDepStoreVars( std::vector< NodeID >& dsVars );
  virtual void registerDepStoreGIDs(const std::vector< std::vector< int > > & dsGIDs) {}

  virtual void getRowColPairs(std::list<index_pair> & rcList) {}

  const std::vector<int> & get_DepSolnGIDJacVec() { return depSolnGIDJacVec_; }

  virtual std::ostream & put(std::ostream & os) const = 0;

  friend std::ostream & operator << (std::ostream & os, const CktNode & cn);

};

//-----------------------------------------------------------------------------
// Function      : CktNode::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/01/01
//-----------------------------------------------------------------------------
inline void CktNode::getDepSolnVars( std::vector< NodeID >& dsVars )
{
  dsVars.clear();
}

//-----------------------------------------------------------------------------
// Function      : CktNode::getDepStateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/01/01
//-----------------------------------------------------------------------------
inline void CktNode::getDepStateVars( std::vector< NodeID >& dsVars )
{
  dsVars.clear();
}

//-----------------------------------------------------------------------------
// Function      : CktNode::getDepStoreVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
inline void CktNode::getDepStoreVars( std::vector< NodeID >& dsVars )
{
  dsVars.clear();
}

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::CktNode N_TOP_CktNode;

#endif
