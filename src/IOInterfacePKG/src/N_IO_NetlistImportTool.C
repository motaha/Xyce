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
// Revision Number: $Revision: 1.66.2.4 $
//
// Revision Date  : $Date: 2014/01/29 18:42:10 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif

#include <iostream>
#include <list>
#include <string>
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_IO_NetlistImportTool.h>

#include <N_IO_DistributionTool.h>

#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_Comm.h>
#endif

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

// ----------   Macro Definitions ---------

#define RETURN_ON_FAILURE(result) \
   if ( false == (result) )       \
   {                              \
      return false;               \
   }

// ---------- EXTERN DECLARATIONS ----------


// this is defined in N_DEV_SourceData.h
// device indices:
//enum Src_index{
//    _SIN_DATA,      // 0
//    _EXP_DATA,      // 1
//    _PULSE_DATA,    // 2
//    _PWL_DATA,      // 3
//    _SFFM_DATA,     // 4
//    _DC_DATA,    // 5
//    _SMOOTH_PULSE_DATA, // 6
//    _AC_DATA,     //tmei: 05/02
//    _NUM_DATA       // total number of data types
//};

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
NetlistImportTool* NetlistImportTool::factory(N_IO_CmdParse & cp,
    Xyce::Topology::Manager & tm)
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
bool NetlistImportTool::registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  return true;
}

typedef std::list<N_UTL_Param> ParameterList;

namespace { // <unnamed>

void setLeadCurrentDevices(const ParameterList &variableList_, N_DEV_DeviceInterface * devPtr_)
{
  std::set<std::string> devicesNeedingLeadCurrents;

  for (ParameterList::const_iterator iterParam = variableList_.begin() ; iterParam != variableList_.end(); ++iterParam)
  {
    std::string varType(iterParam->tag());

    if (iterParam->hasExpressionTag())
    {
      std::vector<std::string> leads;

      N_UTL_Expression exp;
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
      int numIndices = iterParam->iVal();
      if (numIndices > 0 &&(varType == "I" ||(varType.size() == 2 && varType[0] == 'I')))
      {
        // any devices found in this I(xxx) structure need to be communicated to the device manager
        // so that the lead currents can be calculated
        if (numIndices != 1)
        {
          std::string msg("Only one device argument allowed in I() in .print");
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
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

#ifdef Xyce_DEBUG_IO
  // for debugging output the list of devices that need lead currents.
  if (! devicesNeedingLeadCurrents.empty())
  {
    std::set<std::string>::iterator currentDeviceNameItr = devicesNeedingLeadCurrents.begin();
    std::set<std::string>::iterator endDeviceNameItr = devicesNeedingLeadCurrents.end();
    std::cout << "OutputMgr::printLineDiagnostics Devices for which lead currents were requested: ";
    while ( currentDeviceNameItr != endDeviceNameItr)
    {
      std::cout << *currentDeviceNameItr << "  ";
      currentDeviceNameItr++;
    }
  }
#endif

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
int NetlistImportTool::constructCircuitFromNetlist( string const & netlistFile,
    const vector< pair< string, string> > & externalNetlistParams )
{
    map<string,RefCountPtr<N_DEV_InstanceBlock> > dNames;
    set<string> nNames;

    // build metadata
    metadata_.buildMetadata();

    // Create a circuitBlock instance to hold netlist circuit data.
    circuitBlock_ = new N_IO_CircuitBlock(
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
#ifdef Xyce_DEBUG_IO
    cout << "Starting netlist parsing." << endl;
#endif

    bool result;

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();  // All procs call (2)
#endif

    // Get the device insertion tool. Note: the argument determines
    // the insertion tool returned, since only the device insertion
    // is available at this point (08/27/2003) only it is available
    // and the argument is does nothing. This may need to be changed
    // later.

    Xyce::Topology::InsertionTool* insertToolPtr = topMgr_.getInsertionTool("");
    circuitBlock_->registerInsertionTool(insertToolPtr);

    // Regsiter the device interface with the circuit block.
    circuitBlock_->registerDeviceInterface(devIntPtr_);

#ifdef Xyce_PARALLEL_MPI

    // Create the distribution tool (parallel)
    if (pdsCommPtr_ != NULL)
    {
      distToolPtr_ = new N_IO_DistributionTool(circuitBlock_, commandLine_, pdsCommPtr_);
    }
    else
    {
      string msg("Parallel Services not registered with NetlistImportTool\n");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
    }

#else

    // Create the distribution tool (serial)
    distToolPtr_ = new N_IO_DistributionTool( circuitBlock_, commandLine_, NULL );

#endif
    distToolPtr_->registerPkgOptionsMgr(pkgOptMgrPtr_);

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);  // All procs call (2)
#endif

    // Issue the distribution tool start. Proc 0 will return right away, the
    // other procs will return when proc 0 issues the done signal.
    result = distToolPtr_->start();


    bool ok = true;
#ifdef Xyce_PARALLEL_MPI

    if( isSerial_ || pdsCommPtr_->procID() == 0 )
    {

#endif

      circuitBlock_->registerDistributionTool(distToolPtr_);
      ok = circuitBlock_->parseNetlistFilePass1();

      if (ok)
        ok = circuitBlock_->parseNetlistFilePass2();

#ifdef Xyce_PARALLEL_MPI

    }
    N_ERH_ErrorMgr::safeBarrier(0);   // All procs call (END)
    int stat;
    if (!ok)
    {
      stat = -1;
      pdsCommPtr_->bcast (&stat, 1, 0);
    }

    stat = ok?1:0;
    pdsCommPtr_->bcast (&stat, 1, 0);
    ok = (stat != 0);
#endif

    if (!ok)
    {
      string msg("Unrecoverable error(s) occurred in parsing");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
    }

    distToolPtr_->done();

#ifdef Xyce_DEBUG_IO
    cout << "Completed netlist parsing. ";
#endif

  circuitBlock_->fixupYDeviceNames();

  if (outputMgrPtr_->getVariableList())
    setLeadCurrentDevices(*outputMgrPtr_->getVariableList(), devIntPtr_);

  outputMgrPtr_->printLineDiagnostics();

  if ( commandLine_.argExists("-syntax") )
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                         " ***** Netlist syntax OK\n");

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                         " ***** Device Type Counts ...\n");
    outputMgrPtr_->printDeviceCounts();

    Xyce_exit( 0 );
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
bool NetlistImportTool::registerDevMgr(N_DEV_DeviceInterface * devPtr)
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
bool NetlistImportTool::registerOutputMgr(N_IO_OutputMgr * outputMgrPtr)
{
   return ( (outputMgrPtr_ = outputMgrPtr) != NULL );
}

#ifdef Xyce_PARALLEL_MPI
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
  pdsCommPtr_ = tmp_pds_ptr;

  if( !pdsCommPtr_ ) return false;

  isSerial_ = ( pdsCommPtr_->numProc() == 1 );

#ifdef Xyce_DEBUG_DISTRIBUTION
  const string numProcsMsg("Number of Processors:\t");
  const string procIDMsg("Processor ID:\t\t");

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG, numProcsMsg,
			pdsCommPtr_->numProc(), "\n" );
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG, procIDMsg,
			pdsCommPtr_->procID(), "\n" );
#endif

  return true;
}
#endif

// ********** END PUBLIC FUNCTIONS            **********

// ********** BEGIN PROTECTED FUNCTIONS       **********

//-------------------------------------------------------------------------
// Name          : NetlistImportTool::NetlistImportTool
// Purpose       : Default constructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool::NetlistImportTool(
    N_IO_CmdParse & cp,
    Xyce::Topology::Manager & tm)
  : circuitBlock_(NULL),
    commandLine_(cp),
    topMgr_(tm),
    isSerial_(true),
    distToolPtr_(NULL),
    useCount_(0),
    globalParamsInserted_(false),
    currentContextPtr_(NULL),
    circuitContext_( metadata_, contextList_, currentContextPtr_ )
#ifdef Xyce_PARALLEL_MPI
    ,
    pdsCommPtr_(NULL)
#endif
{
}

// ********** END PROTECTED FUNCTIONS         **********
