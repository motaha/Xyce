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
// Revision Number: $Revision: 1.44.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktNode_h
#define N_TOP_CktNode_h 1

// ---------- Standard Includes ----------
#include <iosfwd>
#include <string>
#include <list>
#include <vector>
#include <map>

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_TOP_Misc.h>

// ---------- Forward Declarations ----------

class N_TOP_Topology;
class N_TOP_NodeBlock;

class N_TOP_GraphNode;

//-----------------------------------------------------------------------------
// Class         : N_TOP_CktNode
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class N_TOP_CktNode
{

public:

  // Constructor
  N_TOP_CktNode(const string & nodeID = string(""),
                const int & globalID = 0,
                const list < int > & varGIDList = list < int > (),
                const list < int > & statevarGIDList = list < int > (),
                const list < int > & storevarGIDList = list < int > (),
                const int & pNum = 0, const bool & owned = true,
                N_TOP_GraphNode * gN = 0 )
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
  N_TOP_CktNode( const N_TOP_NodeBlock & nb,
                 N_TOP_GraphNode * gN = 0 );

  // Auto generated copy construct.

  // Destructor
  virtual ~N_TOP_CktNode() { }

  bool operator == (const N_TOP_CktNode & right) const { return id_ == right.id_; }
  bool operator != (const N_TOP_CktNode & right) const { return id_ != right.id_; }

  //-------- get and set methods for attributes

  // Get circuit node id.
  const string & get_id() const { return id_; }
  // Set circuit node id.
  void set_id(const string & value) { id_ = value; }

  // Get circuit node global id.
  const int & get_gID() const { return gID_; }
  // Set circuit node global id.
  void set_gID(const int & globalID) { gID_ = globalID; }

  virtual int type() const = 0;

  //The following get/set functions are only used in N_TOP_CktNode_V class
  //to access and change the noDCPath_ and connToOneTermVar_ internal variables
  //These boolean variables are not defined in the other derived classes, so
  //we create dummy functions there that do nothing.

  virtual bool getNoDCPathVar() = 0;

  virtual bool getConnToOneTermVar() = 0;

  virtual void setTrueNoDCPathVar() = 0;

  virtual void setTrueConnToOneTermVar() = 0;



  const list<int> & get_SolnVarGIDList() const
  { return solnVarGIDList_; }
  void set_SolnVarGIDList(const list<int> & svGIDList)
  { solnVarGIDList_ = svGIDList; }
  void set_SolnVarGIDVec(const vector<int> & svGIDVec)
  { solnVarGIDList_.assign( svGIDVec.begin(), svGIDVec.end() ); }

  const list<int> & get_ExtSolnVarGIDList() const
  { return extSolnVarGIDList_; }
  void set_ExtSolnVarGIDList(const list<int> & svGIDList)
  { extSolnVarGIDList_ = svGIDList; }
  void set_ExtSolnVarGIDVec(const vector<int> & svGIDVec)
  { extSolnVarGIDList_.assign( svGIDVec.begin(), svGIDVec.end() ); }

  const list < int > & get_StateVarGIDList() const
  { return stateVarGIDList_; }
  void set_StateVarGIDList(const list < int > & svGIDList)
  { stateVarGIDList_ = svGIDList; }
  void set_StateVarGIDVec(const vector < int > & svGIDVec)
  { stateVarGIDList_.assign( svGIDVec.begin(), svGIDVec.end() ); }

  const list < int > & get_StoreVarGIDList() const
  { return storeVarGIDList_; }
  void set_StoreVarGIDList(const list < int > & svGIDList)
  { storeVarGIDList_ = svGIDList; }
  void set_StoreVarGIDVec(const vector < int > & svGIDVec)
  { storeVarGIDList_.assign( svGIDVec.begin(), svGIDVec.end() ); }

  const vector< vector<int> > & get_DepSolnGIDVec() const
  { return depSolnGIDVec_; }
  void set_DepSolnGIDVec(const vector< vector<int> > & dsGIDs)
  { depSolnGIDVec_ = dsGIDs; }
  const vector< vector<int> > & get_DepStateGIDVec() const
  { return depStateGIDVec_; }
  void set_DepStateGIDVec(const vector< vector<int> > & dsGIDs)
  { depStateGIDVec_ = dsGIDs; }
  const vector< vector<int> > & get_DepStoreGIDVec() const
  { return depStoreGIDVec_; }
  void set_DepStoreGIDVec(const vector< vector<int> > & dsGIDs)
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

  N_TOP_GraphNode * get_GraphNode() { return graphNodePtr_; }
  void set_GraphNode(N_TOP_GraphNode * gN) { graphNodePtr_ = gN; }

  N_TOP_NodeBlock * extractNodeBlock();

  virtual const vector< vector<int> > & jacobianStamp() const
  { static vector< vector<int> > dummy; return dummy; }
  virtual void registerJacLIDswithDev( const vector< vector<int> > & jacLIDVec ) {}

  virtual void registerGIDDataWithDev(
               const vector<int> & counts,
               const vector<int> & GIDs,
               const vector< vector<int> > & jacGIDs ) {}

  //------- Added for use with the outputFileName function
  virtual map<int,string> & getIntNameMap() { return intNameMap_; }
  virtual map<int,string> & getStateNameMap() { return stateNameMap_; }
  virtual map<int,string> & getStoreNameMap() { return storeNameMap_; }

  virtual void varTypeList( vector<char> & varTypeVec ) {}

protected:

  // Node id.
  string id_;
  // Global node id.
  int gID_;

  // Processor number.
  int procNum_;
  // Processor ownership flag.
  bool isOwned_;

  int Offset_;

  N_TOP_GraphNode * graphNodePtr_;

  // Solution variable global id list.
  list<int> solnVarGIDList_;
  list<int> extSolnVarGIDList_;
  // State variable global id list.
  list<int> stateVarGIDList_;
  list<int> storeVarGIDList_;

  vector< vector<int> > depSolnGIDVec_;
  vector< vector<int> > depStateGIDVec_;
  vector< vector<int> > depStoreGIDVec_;

  vector<int> depSolnGIDJacVec_;

  map<int,string> intNameMap_;
  map<int,string> stateNameMap_;
  map<int,string> storeNameMap_;

public:

  //------- registration of global ids with device instances

  virtual int solnVarCount() { return solnVarGIDList_.size(); }
  virtual int extSolnVarCount() { return extSolnVarGIDList_.size(); }
  virtual int stateVarCount() { return stateVarGIDList_.size(); }
  virtual int depSolnVarCount() { return depSolnGIDVec_.size(); }
  virtual int depStateVarCount() { return depStateGIDVec_.size(); }

  virtual int storeVarCount() { return storeVarGIDList_.size(); }
  virtual int depStoreVarCount() { return depStoreGIDVec_.size(); }

  virtual void leadConnect(vector<int> &) {return;}

  virtual void registerGIDswithDev( const list<index_pair> & intGIDList,
  	                            const list<index_pair> & extGIDList) {}

  virtual void registerStateGIDswithDev( const list<index_pair> & stateGIDList) {}
  virtual void registerStoreGIDswithDev( const list<index_pair> & storeGIDList) {}

  virtual void registerLIDswithDev( const vector<int> & intLIDVec,
                                    const vector<int> & extLIDVec ) {}
  virtual void registerStateLIDswithDev( const vector<int> & stateLIDVec ) {}
  virtual void registerStoreLIDswithDev( const vector<int> & storeLIDVec ) {}

  virtual void registerDepLIDswithDev( const vector< vector<int> > & depLIDVec ) {}
  virtual void registerDepStateLIDswithDev( const vector< vector<int> > & depStateLIDVec ) {}
  virtual void registerDepStoreLIDswithDev( const vector< vector<int> > & depStoreLIDVec ) {}

  virtual void getDepSolnVars( vector< NodeID >& dsVars );
  virtual void registerDepSolnGIDs( const vector< vector<int> > & dsGIDs) {}
  virtual void getDepStateVars( vector< NodeID >& dsVars );
  virtual void registerDepStateGIDs(const vector < vector < int > > & dsGIDs) {}
  virtual void getDepStoreVars( vector< NodeID >& dsVars );
  virtual void registerDepStoreGIDs(const vector < vector < int > > & dsGIDs) {}

  virtual void getRowColPairs(list<index_pair> & rcList) {}

  const vector<int> & get_DepSolnGIDJacVec() { return depSolnGIDJacVec_; }

  virtual ostream & put(ostream & os) const = 0;

  friend ostream & operator << (ostream & os, const N_TOP_CktNode & cn);

};

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/01/01
//-----------------------------------------------------------------------------
inline void N_TOP_CktNode::getDepSolnVars( vector< NodeID >& dsVars )
{
  dsVars.clear();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode::getDepStateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/01/01
//-----------------------------------------------------------------------------
inline void N_TOP_CktNode::getDepStateVars( vector< NodeID >& dsVars )
{
  dsVars.clear();
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode::getDepStoreVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
inline void N_TOP_CktNode::getDepStoreVars( vector< NodeID >& dsVars )
{
  dsVars.clear();
}

#endif
