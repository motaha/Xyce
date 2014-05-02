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
// Filename      : NetlistImportTool.C
//
// Purpose       : Implement the interface to to read and parse a netlist for
//                 an electrical circuit.
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.92 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <N_UTL_fwd.h>

#include <N_IO_NetlistImportTool.h>

#include <N_IO_DistributionTool.h>

#include <N_PDS_ParallelMachine.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

#include <N_IO_DeviceBlock.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_CircuitBlock.h>

#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInterface.h>

#include <N_IO_CmdParse.h>

#include <N_ERH_ErrorMgr.h>

#include <N_TOP_NodeBlock.h>
#include <N_TOP_NodeDevBlock.h>
#include <N_TOP_InsertionTool.h>
#include <N_TOP_TopologyMgr.h>
#include <N_DEV_SourceData.h>
#include <N_UTL_Xyce.h>

#include <N_UTL_Expression.h>

namespace Xyce {
namespace IO {

//-------------------------------------------------------------------------
// Function      : NetlistImportTool::factory
// Purpose       : This function returns a pointer to the only
//                 instance of this class. (for the current library
//                 instance, hopefully)
//
// Special Notes : ERK.  10/16/2005.  This used to be a singleton (ie a
//                 static pointer was returned) but had to be changed
//                 so that the library version of Xyce would work
//                 correctly.
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool* NetlistImportTool::factory(IO::CmdParse & cp,
    Xyce::Topo::Manager & tm)
{
  NetlistImportTool * NIT_ptr = new NetlistImportTool(cp,tm);
  return NIT_ptr;
}



//-------------------------------------------------------------------------
// Function      : NetlistImportTool::~NetlistImportTool
// Purpose       : Destructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool::~NetlistImportTool()
{
    delete circuitBlock_;
    delete distToolPtr_;
}

//-----------------------------------------------------------------------------
// Function      : NetlistImportTool::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerPkgOptionsMgr( PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  return true;
}

typedef std::list<Util::Param> ParameterList;

namespace { // <unnamed>

void setLeadCurrentDevices(const ParameterList &variableList_, Device::DeviceInterface * devPtr_)
{
  std::set<std::string> devicesNeedingLeadCurrents;

  for (ParameterList::const_iterator iterParam = variableList_.begin() ; iterParam != variableList_.end(); ++iterParam)
  {
    std::string varType(iterParam->tag());

    if (Util::hasExpressionTag(*iterParam))
    {
      std::vector<std::string> leads;

      Util::Expression exp;
      exp.set(iterParam->tag());
      exp.get_names(XEXP_LEAD, leads);

      // any lead currents found in this expression need to be communicated to the device manager.
      // Multi terminal devices have an extra designator on the name as in name{lead_name}
      // need to remove any {} in the name.
      for (std::vector<std::string>::const_iterator currLeadItr = leads.begin(); currLeadItr != leads.end(); ++currLeadItr)
      {
        size_t leadDesignator = currLeadItr->find_first_of("{");
        devicesNeedingLeadCurrents.insert( currLeadItr->substr(0, leadDesignator));
      }
    }
    else
    {
      int numIndices = iterParam->getImmutableValue<int>();
      if (numIndices > 0 &&(varType == "I" ||(varType.size() == 2 && varType[0] == 'I')))
      {
        // any devices found in this I(xxx) structure need to be communicated to the device manager
        // so that the lead currents can be calculated
        if (numIndices != 1)
        {
          Report::DevelFatal0() << "Only one device argument allowed in I() in .print";
        }
        ++iterParam;
        if (varType.size() == 2)
        {
          devicesNeedingLeadCurrents.insert(iterParam->tag());
        }
        else
        {
          devicesNeedingLeadCurrents.insert(iterParam->tag());
        }
      }
    }
  }

  // for debugging output the list of devices that need lead currents.
  if ( Xyce::DEBUG_IO && !devicesNeedingLeadCurrents.empty())
  {
    std::set<std::string>::iterator currentDeviceNameItr = devicesNeedingLeadCurrents.begin();
    std::set<std::string>::iterator endDeviceNameItr = devicesNeedingLeadCurrents.end();
    Xyce::dout() << "OutputMgr::printLineDiagnostics Devices for which lead currents were requested: ";
    while ( currentDeviceNameItr != endDeviceNameItr)
    {
      Xyce::dout() << *currentDeviceNameItr << "  ";
      currentDeviceNameItr++;
    }
  }

// This is the last call on the output manager before devices are constructed
  // So it's the last time I can isolate lead currents. However it won't be sufficient
  // as we haven't parsed expressions in devices yet.  I'll need to rethink
  // when this call is made to the device manager.  Sometime just after device instance
  // construction but before the store vector is allocated.
  devPtr_->setLeadCurrentRequests(devicesNeedingLeadCurrents);
}

} // namespace <unnamed>

//-------------------------------------------------------------------------
// Function      : NetlistImportTool::constructCircuitFromNetlist
// Purpose       : Construct a circuit from a netlist.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
int NetlistImportTool::constructCircuitFromNetlist(const std::string &netlistFile,
    const std::vector< std::pair< std::string, std::string> > & externalNetlistParams )
{
    std::map<std::string,RefCountPtr<Device::InstanceBlock> > dNames;
    std::set<std::string> nNames;

    // build metadata
    metadata_.buildMetadata();

    // Create a circuitBlock instance to hold netlist circuit data.
    circuitBlock_ = new IO::CircuitBlock(
        netlistFile,
        commandLine_,
        metadata_,
        modelNames_ ,
        ssfMap_,
        circuitContext_,
        outputMgrPtr_,
        useCount_,
        globalParamsInserted_,
        dNames,
        nNames,
        externalNetlistParams );

    // Parse the netlist file.
    if (Xyce::DEBUG_IO)
      Xyce::dout() << "Starting netlist parsing." << std::endl;

    bool result;

    // Get the device insertion tool. Note: the argument determines
    // the insertion tool returned, since only the device insertion
    // is available at this point (08/27/2003) only it is available
    // and the argument is does nothing. This may need to be changed
    // later.

    Xyce::Topo::InsertionTool* insertToolPtr = topMgr_.getInsertionTool("");
    circuitBlock_->registerInsertionTool(insertToolPtr);

    // Regsiter the device interface with the circuit block.
    circuitBlock_->registerDeviceInterface(devIntPtr_);

    distToolPtr_ = new IO::DistributionTool(circuitBlock_, commandLine_, pdsCommPtr_);

    distToolPtr_->registerPkgOptionsMgr(pkgOptMgrPtr_);

    N_ERH_ErrorMgr::safeBarrier(comm_);  // All procs call (2)

    // Issue the distribution tool start. Proc 0 will return right away, the
    // other procs will return when proc 0 issues the done signal.
    result = distToolPtr_->start();

    bool ok = true;

    if (Parallel::rank(comm_) == 0)
    {
      circuitBlock_->registerDistributionTool(distToolPtr_);
      ok = circuitBlock_->parseNetlistFilePass1();

      if (ok)
        ok = circuitBlock_->parseNetlistFilePass2();

      if (!ok)
        N_ERH_ErrorMgr::safeBarrier(comm_);
    }

    if (ok) {
      distToolPtr_->done();

      if (Xyce::DEBUG_IO)
        Xyce::dout() << "Completed netlist parsing. ";

      circuitBlock_->fixupYDeviceNames();
      if (outputMgrPtr_->getVariableList())
        setLeadCurrentDevices(*outputMgrPtr_->getVariableList(), devIntPtr_);
      outputMgrPtr_->printLineDiagnostics();
    }

  return 1;
}



//-------------------------------------------------------------------------
// Function      : NetlistImportTool::registerDevMgr
// Purpose       :
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
bool NetlistImportTool::registerDevMgr(Device::DeviceInterface * devPtr)
{
   return ( (devIntPtr_ = devPtr) != NULL );
}

//-------------------------------------------------------------------------
// Function      : NetlistImportTool::registerOutputMgr
// Purpose       :
//
// Special Notes :
//
// Creator       : Dave Shirley, PSSI
//
// Creation Date : 08/15/06
//-------------------------------------------------------------------------
bool NetlistImportTool::registerOutputMgr(IO::OutputMgr * outputMgrPtr)
{
   return ( (outputMgrPtr_ = outputMgrPtr) != NULL );
}

//-----------------------------------------------------------------------------
// Function      : NetlistImportTool::registerParallelServices
// Purpose       : Registers N_PDS_Comm object for parallel communication.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/17/2003
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerParallelServices(N_PDS_Comm * tmp_pds_ptr)
{
  if( !tmp_pds_ptr )
    return false;

  pdsCommPtr_ = tmp_pds_ptr;
  comm_ = pdsCommPtr_->comm();

  if (Xyce::DEBUG_DISTRIBUTION)
    Xyce::dout() << "Number of Processors:   " << Parallel::size(comm_) << std::endl
                 << "Processor ID:           " << Parallel::rank(comm_) << std::endl;

  return true;
}

//-------------------------------------------------------------------------
// Name          : NetlistImportTool::NetlistImportTool
// Purpose       : Default constructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool::NetlistImportTool(
    IO::CmdParse & cp,
    Xyce::Topo::Manager & tm)
  : circuitBlock_(NULL),
    commandLine_(cp),
    topMgr_(tm),
    distToolPtr_(NULL),
    useCount_(0),
    globalParamsInserted_(false),
    currentContextPtr_(NULL),
    circuitContext_( metadata_, contextList_, currentContextPtr_ ),
    comm_(0),
    pdsCommPtr_(0)
{}

} // namespace IO
} // namespace Xyce
