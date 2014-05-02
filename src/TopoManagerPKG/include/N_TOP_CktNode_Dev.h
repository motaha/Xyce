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
// Revision Number: $Revision: 1.42 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktNode_Dev_h
#define N_TOP_CktNode_Dev_h 1

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_TOP_CktNode.h>
#include <N_DEV_fwd.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktNode_Dev
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktNode_Dev : public CktNode
{

private:


  // Default constructor (private).
  CktNode_Dev() {}

  // Assignment operator (private).
  CktNode_Dev & operator=(const CktNode_Dev & right) { return *this; }

public:

  // Constructor
  CktNode_Dev(Device::DeviceInstance * devinstPtr,
                    const std::string & ID = std::string(""), const int & gID = 0,
                    const std::list< int > & varGIDList = std::list< int > (),
                    const std::list< int > & statevarGIDList = std::list< int > (),
                    const std::list< int > & storevarGIDList = std::list< int > (),
                    const int & pNum = 0, const bool & owned = true)
  : CktNode(ID, gID, varGIDList, statevarGIDList, storevarGIDList, pNum, owned),
    devPtr_(devinstPtr)
  {}

  // Constructor
  CktNode_Dev(Device::DeviceInstance * devinstPtr,
                    const NodeBlock & nb)
  : CktNode(nb),
    devPtr_(devinstPtr)
  {}

  CktNode_Dev( const NodeBlock & nb,
		     const Teuchos::RefCountPtr<Device::InstanceBlock> ibPtr,
		     Device::DeviceInterface & devIF )
  : CktNode(nb),
    devPtr_(0),
    iface_(&devIF),
    instance_( ibPtr )
  {}

  // Destructor
  ~CktNode_Dev();

  int type() const { return _DNODE; }

  bool getNoDCPathVar() {return false;}

  bool getConnToOneTermVar() {return false;}

  void setTrueNoDCPathVar() {}

  void setTrueConnToOneTermVar() {}

  bool instantiated() const { return (devPtr_!=0); }
  bool instantiate();

  Teuchos::RCP<Device::InstanceBlock> devBlock() { return instance_; }

  // Get's the device state object.
  Device::DeviceState * getDevState();

  // Set's the device state object.
  bool setDevState(const Device::DeviceState & state);

  // Registers int. and ext. global ids with dev instance.
  void registerGIDswithDev(const std::list< index_pair > & intGIDList,
  	const std::list< index_pair > & extGIDList);

  // Registers state global ids with dev instance.
  void registerStateGIDswithDev(const std::list< index_pair > & stateGIDList);
  void registerStoreGIDswithDev(const std::list< index_pair > & storeGIDList);

  void registerLIDswithDev( const std::vector<int> & intLIDVec,
                            const std::vector<int> & extLIDVec );
  void registerStateLIDswithDev( const std::vector<int> & stateLIDVec );
  void registerStoreLIDswithDev( const std::vector<int> & storeLIDVec );

  void registerDepLIDswithDev( const std::vector< std::vector<int> > & depLIDVec );
  void registerDepStateLIDswithDev( const std::vector< std::vector<int> > & depStateLIDVec );
  void registerDepStoreLIDswithDev( const std::vector< std::vector<int> > & depStoreLIDVec );

  // Setup secondary dependencies.
  void getDepSolnVars( std::vector< NodeID >& dsVars );
  void registerDepSolnGIDs(const std::vector< std::vector< int > > & dsGIDs);
  void getDepStateVars( std::vector< NodeID >& dsVars );
  void registerDepStateGIDs(const std::vector< std::vector< int > > & dsGIDs);
  void getDepStoreVars( std::vector< NodeID >& dsVars );
  void registerDepStoreGIDs(const std::vector< std::vector< int > > & dsGIDs);

  const std::vector<std::string> & getDepStoreVars();


  // Get RowCol pairs from devices.
  void getRowColPairs(std::list< index_pair > & rcList);

  int solnVarCount();
  int stateVarCount();
  int storeVarCount();
  void leadConnect(std::vector<int> &);

  // Added for use with the outputFileName function.
  std::map< int, std::string > & getIntNameMap();
  std::map< int, std::string > & getStateNameMap();
  std::map< int, std::string > & getStoreNameMap();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDswithDev( const std::vector< std::vector<int> > & jacLIDVec );

  void registerGIDDataWithDev(
        const std::vector<int> & counts,
        const std::vector<int> & GIDs,
  	const std::vector< std::vector<int> > & jacGIDs );

  void varTypeList( std::vector<char> & varTypeVec );

private:

  // Pointer to a device instance.
  Device::DeviceInstance * devPtr_;

  Device::DeviceInterface * iface_;
  Teuchos::RefCountPtr<Device::InstanceBlock> instance_;

public:

    std::ostream & put(std::ostream & os) const;

};

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::CktNode_Dev N_TOP_CktNode_Dev;

#endif
