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
// Filename       : $RCSfile: N_TOP_CktNode_Dev.h,v $
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
// Revision Number: $Revision: 1.37.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktNode_Dev_h
#define N_TOP_CktNode_Dev_h 1

// ---------- Standard Includes ----------
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_TOP_CktNode.h>
#include <N_DEV_fwd.h>

// ---------- Forward Declarations ----------

//-----------------------------------------------------------------------------
// Class         : N_TOP_CktNode_Dev
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class N_TOP_CktNode_Dev : public N_TOP_CktNode
{

private:


  // Default constructor (private).
  N_TOP_CktNode_Dev() {}

  // Assignment operator (private).
  N_TOP_CktNode_Dev & operator=(const N_TOP_CktNode_Dev & right) { return *this; }

public:

  // Constructor
  N_TOP_CktNode_Dev(N_DEV_DeviceInstance * devinstPtr,
                    const string & ID = string(""), const int & gID = 0,
                    const list < int > & varGIDList = list < int > (),
                    const list < int > & statevarGIDList = list < int > (),
                    const list < int > & storevarGIDList = list < int > (),
                    const int & pNum = 0, const bool & owned = true)
  : N_TOP_CktNode(ID, gID, varGIDList, statevarGIDList, storevarGIDList, pNum, owned),
    devPtr_(devinstPtr)
  {}

  // Constructor
  N_TOP_CktNode_Dev(N_DEV_DeviceInstance * devinstPtr,
                    const N_TOP_NodeBlock & nb)
  : N_TOP_CktNode(nb),
    devPtr_(devinstPtr)
  {}

  N_TOP_CktNode_Dev( const N_TOP_NodeBlock & nb,
		     const Teuchos::RefCountPtr<N_DEV_InstanceBlock> ibPtr,
		     N_DEV_DeviceInterface & devIF )
  : N_TOP_CktNode(nb),
    devPtr_(0),
    iface_(&devIF),
    instance_( ibPtr )
  {}

  // Destructor
  ~N_TOP_CktNode_Dev();

  int type() const { return _DNODE; }

  bool getNoDCPathVar() {return false;}

  bool getConnToOneTermVar() {return false;}

  void setTrueNoDCPathVar() {}

  void setTrueConnToOneTermVar() {}

  bool instantiated() const { return (devPtr_!=0); }
  bool instantiate();

  Teuchos::RCP<N_DEV_InstanceBlock> devBlock() { return instance_; }

  // Get's the device state object.
  N_DEV_DeviceState * getDevState();

  // Set's the device state object.
  bool setDevState(const N_DEV_DeviceState & state);

  // Registers int. and ext. global ids with dev instance.
  void registerGIDswithDev(const list < index_pair > & intGIDList,
  	const list < index_pair > & extGIDList);

  // Registers state global ids with dev instance.
  void registerStateGIDswithDev(const list < index_pair > & stateGIDList);
  void registerStoreGIDswithDev(const list < index_pair > & storeGIDList);

  void registerLIDswithDev( const vector<int> & intLIDVec,
                            const vector<int> & extLIDVec );
  void registerStateLIDswithDev( const vector<int> & stateLIDVec );
  void registerStoreLIDswithDev( const vector<int> & storeLIDVec );

  void registerDepLIDswithDev( const vector< vector<int> > & depLIDVec );
  void registerDepStateLIDswithDev( const vector< vector<int> > & depStateLIDVec );
  void registerDepStoreLIDswithDev( const vector< vector<int> > & depStoreLIDVec );

  // Setup secondary dependencies.
  void getDepSolnVars( vector< NodeID >& dsVars );
  void registerDepSolnGIDs(const vector < vector < int > > & dsGIDs);
  void getDepStateVars( vector< NodeID >& dsVars );
  void registerDepStateGIDs(const vector < vector < int > > & dsGIDs);
  void getDepStoreVars( vector< NodeID >& dsVars );
  void registerDepStoreGIDs(const vector < vector < int > > & dsGIDs);

  const vector<string> & getDepStoreVars();


  // Get RowCol pairs from devices.
  void getRowColPairs(list < index_pair > & rcList);

  int solnVarCount();
  int stateVarCount();
  int storeVarCount();
  void leadConnect(vector<int> &);

  // Added for use with the outputFileName function.
  map < int, string > & getIntNameMap();
  map < int, string > & getStateNameMap();
  map < int, string > & getStoreNameMap();

  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDswithDev( const vector< vector<int> > & jacLIDVec );

  void registerGIDDataWithDev(
        const vector<int> & counts,
        const vector<int> & GIDs,
  	const vector< vector<int> > & jacGIDs );

  void varTypeList( vector<char> & varTypeVec );

private:

  // Pointer to a device instance.
  N_DEV_DeviceInstance * devPtr_;

  N_DEV_DeviceInterface * iface_;
  Teuchos::RefCountPtr<N_DEV_InstanceBlock> instance_;

public:

  ostream & put(ostream & os) const;

};

#endif
