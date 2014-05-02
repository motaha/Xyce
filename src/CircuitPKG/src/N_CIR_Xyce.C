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
// Filename       : $RCSfile: N_CIR_Xyce.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/27/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.276.2.3 $
//
// Revision Date  : $Date: 2014/03/10 16:15:20 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <stdexcept>

// ----------   Xyce Includes   ----------

#include <N_CIR_Xyce.h>

#include <N_DEV_fwd.h>
#include <N_DEV_RegisterDevices.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DAC.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternalSimulationData.h>
#include <N_DEV_Print.h>

#include <N_ERH_ErrorMgr.h>

#include <N_IO_NetlistImportTool.h>

#include <N_IO_OutputMgr.h>
#include <N_IO_RestartMgr.h>
#include <N_IO_PkgOptionsMgr.h>

#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_LOA_LoaderMgr.h>
#include <N_LOA_Loader.h>
#include <N_LOA_CktLoader.h>
#include <N_LOA_NonlinearEquationLoader.h>

#include <N_NLS_Manager.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_ANP_AnalysisInterface.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TwoLevelError.h>

#include <N_TOP_Topology.h>
#include <N_TOP_TopologyMgr.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Timer.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_Expression.h>
#include <N_UTL_LogStream.h>

#include <N_UTL_Version.h>
#include <N_UTL_BreakPoint.h>


//--------  Forward Declarations ------------

class N_LAS_QueryUtil;

class N_LAS_MultiVector;

namespace Xyce {
namespace Circuit {

//--------  Global Declarations ------------
void report_handler(const char *message, unsigned report_mask);

//-----------------------------------------------------------------------------
// Function      : Simulator::Simulator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
Simulator::Simulator(Xyce::Parallel::Machine comm)
:
  comm_(comm),
  devIntPtr_(0),
  topPtr_(0),
  topMgrPtr_(0),
  netlistImportToolPtr_(0),
  outMgrPtr_(0),
  lasSysPtr_(0),
  lasBuilderPtr_(0),
  anaIntPtr_(0),
  nlsMgrPtr_(0),
  loaderMgrPtr_(0),
  cktLoaderPtr_(0),
  nonlinearEquationLoaderPtr_(0),
  parMgrPtr_(0),
  resMgrPtr_(0),
  XyceTimerPtr_(0),
  ElapsedTimerPtr_(0),
  multiThreading_(false),
  numThreads_(0),
  initializeAllFlag_(false)
{
  previousReportHandler_ = Xyce::set_report_handler(report_handler);

  Xyce::Device::registerDevices();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::~Xyce
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
Simulator::~Simulator()
{
  Xyce::set_report_handler(previousReportHandler_);
}

//---------------------------------------------------------------------------
// Function      : Simulator::setNetlistParameters
// Purpose       : This passes a vector of pairs "key" "value" that will
//                 be substituted during the processing of the netlist.  This
//                 more easily allows Dakota to change any netlist parameter
//                 during netlist setup.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 10/9/2008
//---------------------------------------------------------------------------
void Simulator::setNetlistParameters( const std::vector< std::pair< std::string, std::string > > & externalParams )
{
  externalNetlistParams_ = externalParams;
}

//---------------------------------------------------------------------------
// Function      : Simulator::setNetlistParameters
// Purpose       : Call through to the output manager to set the suffix to
//                 be used on the output file, as in circuit + suffix + prn
//                 This is useful in Dakota controlled runs to keep each
//                 simulation from overwritting the last one.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 10/9/2008
//---------------------------------------------------------------------------
void Simulator::setOutputFileSuffix( const std::string newSuffix )
{
  if( outMgrPtr_ )
  {
    outMgrPtr_->setOutputFilenameSuffix( newSuffix );
  }
}
//-----------------------------------------------------------------------------
// Function      : Simulator::setupParMgr_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool Simulator::setupParMgr_( int iargs, char **cargs)
{
  bool bsuccess = true;

  // Setup the Parallel Mgr. with a default load-balance based on the numProc
  // value.
  parMgrPtr_ = new N_PDS_Manager( isSerialFlag_, procZeroFlag_, iargs, cargs, comm_);
  if (comm_ == 0)
    comm_ = parMgrPtr_->getPDSComm()->comm();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doAllocations_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::doAllocations_()

{
  bool bsuccess;

  std::string topotype = "Basic";

  // Allocate device manager:
  devIntPtr_  = Xyce::Device::DeviceInterface::factory(commandLine);

  // Allocate Topology:
  topMgrPtr_ = Xyce::Topo::Manager::instance();
  Xyce::Topo::Manager & topMgr = (*topMgrPtr_);

  topPtr_ = topMgrPtr_->createTopology(commandLine);

  // Allocate distribution mgr:
  netlistImportToolPtr_ = Xyce::IO::NetlistImportTool::factory(commandLine, topMgr);

  // Allocate output mgr:
  outMgrPtr_ = N_IO_OutputMgr::factory(commandLine);

  // Allocate restart mgr:
  resMgrPtr_ = N_IO_RestartMgr::factory(commandLine);

  // Linear Algebra allocations:
  lasSysPtr_       = new N_LAS_System();
  lasBuilderPtr_   = new N_LAS_Builder();

  // Allocate Time Integration:
  anaIntPtr_    = new N_ANP_AnalysisInterface(commandLine);

  // Allocate nonlinear solver:
  nlsMgrPtr_ = new N_NLS_Manager(commandLine);

  // Allocate loader:
  loaderMgrPtr_ = new N_LOA_LoaderMgr();

  cktLoaderPtr_ = loaderMgrPtr_->createCktLoader();
  nonlinearEquationLoaderPtr_ = loaderMgrPtr_->createNonlinearEquationLoader();

  pkgOptionsMgrPtr_ = new N_IO_PkgOptionsMgr();

  bsuccess = true;
  bsuccess = bsuccess && (devIntPtr_ != 0);

  bsuccess = bsuccess && (topPtr_ != 0);

  bsuccess = bsuccess && (netlistImportToolPtr_ != 0 );

  bsuccess = bsuccess && (outMgrPtr_ != 0);

  bsuccess = bsuccess && (lasSysPtr_ != 0);
  bsuccess = bsuccess && (lasBuilderPtr_ != 0);

  bsuccess = bsuccess && (anaIntPtr_ != 0);

  bsuccess = bsuccess && (nlsMgrPtr_ != 0);

  bsuccess = bsuccess && (loaderMgrPtr_ != 0);
  bsuccess = bsuccess && (cktLoaderPtr_ != 0);
  bsuccess = bsuccess && (nonlinearEquationLoaderPtr_ != 0);

  bsuccess = bsuccess && (parMgrPtr_ != 0);

  bsuccess = bsuccess && (pkgOptionsMgrPtr_ != 0);

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : Simulator::doDeAllocations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::doDeAllocations_()

{
  // de-allocate the device manager:

  delete devIntPtr_;
  delete nlsMgrPtr_;
  delete anaIntPtr_;
  delete loaderMgrPtr_;
  delete cktLoaderPtr_;
  delete nonlinearEquationLoaderPtr_;
  delete outMgrPtr_;
  delete lasSysPtr_;
  delete lasBuilderPtr_;
  delete XyceTimerPtr_;
  delete ElapsedTimerPtr_;
  delete parMgrPtr_;
  delete topMgrPtr_;
  delete resMgrPtr_;
  delete pkgOptionsMgrPtr_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doRegistrations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::doRegistrations_()

{
  bool bsuccess = true;
  bool bs1 = true;


  // ParallelDist Manager registrations
  bs1 = parMgrPtr_->registerTopology(topPtr_);
  bsuccess = bsuccess && bs1;

  // Device Manager registrations
  bs1 = devIntPtr_->registerLinearSystem(lasSysPtr_);    bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerNonlinearSolver(nlsMgrPtr_); bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerAnalysisInterface(anaIntPtr_);  bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerParallelMgr(parMgrPtr_);     bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerOutputMgr(outMgrPtr_);       bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

  // Topology registrations:
  bs1 = topPtr_->registerDeviceInterface(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = topPtr_->registerParallelMgr(parMgrPtr_);     bsuccess = bsuccess && bs1;
  bs1 = topPtr_->registeranaInt(anaIntPtr_);          bsuccess = bsuccess && bs1;
  bs1 = topPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

  // Distribution manager registrations:
  bs1 = netlistImportToolPtr_->registerParallelServices(parMgrPtr_->getPDSComm());
  bsuccess = bsuccess && bs1;

  bs1 = netlistImportToolPtr_->registerDevMgr(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = netlistImportToolPtr_->registerOutputMgr(outMgrPtr_); bsuccess = bsuccess && bs1;
  bs1 = netlistImportToolPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

  // Restart manager registrations:
  bs1 = resMgrPtr_->registerTopology(topPtr_);           bsuccess = bsuccess && bs1;
  bs1 = resMgrPtr_->registerDeviceInterface(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = resMgrPtr_->registeranaInt(anaIntPtr_);          bsuccess = bsuccess && bs1;
  bs1 = resMgrPtr_->registerParallelServices(parMgrPtr_);bsuccess = bsuccess && bs1;
  bs1 = resMgrPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

  // Output manager registrations:
  bs1 = outMgrPtr_->registerTopology(topPtr_);           bsuccess = bsuccess && bs1;
  bs1 = outMgrPtr_->registerDeviceInterface(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = outMgrPtr_->registerAnalysisInterface(anaIntPtr_);          bsuccess = bsuccess && bs1;
  bs1 = outMgrPtr_->registerPkgOptionsMgr( *pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;

//#ifdef Xyce_PARALLEL_MPI
  bs1 = outMgrPtr_->registerParallelServices( parMgrPtr_->getPDSComm() );
  bsuccess = bsuccess && bs1;
//#endif

  // Analysis manager registrations:
  bs1 = anaIntPtr_->registerLinearSystem(lasSysPtr_); bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerNLSManager(nlsMgrPtr_);   bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerLoader(cktLoaderPtr_);    bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerOutputMgr(outMgrPtr_);    bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerElapsedTimer(ElapsedTimerPtr_); bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerParallelServices(parMgrPtr_); bsuccess = bsuccess && bs1;

  bs1 = anaIntPtr_->registerRestartMgr(resMgrPtr_); bsuccess = bsuccess && bs1;

  // Nonlinear Solver registrations:
  bs1 = nlsMgrPtr_->registerLoader(nonlinearEquationLoaderPtr_);      bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerLinearSystem(lasSysPtr_);   bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerAnalysisInterface(anaIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerOutputMgr(outMgrPtr_);      bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerPkgOptionsMgr( pkgOptionsMgrPtr_ ); bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerTopology(topPtr_); bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->registerParallelMgr(parMgrPtr_); bsuccess = bsuccess && bs1;

  // Loader(s) registrations:
  bs1 = cktLoaderPtr_->registerDeviceInterface(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = nonlinearEquationLoaderPtr_->registerDeviceInterface(devIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = nonlinearEquationLoaderPtr_->registerTIA (anaIntPtr_);            bsuccess = bsuccess && bs1;

  // Linear Solver registrations:
  bs1 = lasSysPtr_->registerPDSManager(parMgrPtr_);    bsuccess = bsuccess && bs1;
  bs1 = lasSysPtr_->registerBuilder( lasBuilderPtr_ ); bsuccess = bsuccess && bs1;
  bs1 = lasSysPtr_->registerQueryUtil(
    (N_LAS_QueryUtil *) topPtr_->get_LinSolvUtil() );  bsuccess = bsuccess && bs1;
  bs1 = lasSysPtr_->registerANPInterface(anaIntPtr_);    bsuccess = bsuccess && bs1;

  bs1 = lasBuilderPtr_->registerPDSManager(parMgrPtr_); bsuccess = bsuccess && bs1;
  bs1 = lasBuilderPtr_->registerSystem(lasSysPtr_);     bsuccess = bsuccess && bs1;
  bs1 = lasBuilderPtr_->registerQueryUtil(
    (N_LAS_QueryUtil *) topPtr_->get_LinSolvUtil() );
  bsuccess = bsuccess && bs1;

  // Additional ANP registrations for MPDE
  bs1 = anaIntPtr_->registerDeviceInterface(devIntPtr_);  bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerTopology(topPtr_);            bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerNonlinearEquationLoader(nonlinearEquationLoaderPtr_);     bsuccess = bsuccess && bs1;
  bs1 = anaIntPtr_->registerApplicationBuilder(lasBuilderPtr_); bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::setUpTopology_
// Purpose       : This function handles a lot of the initial setup.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/00
//-----------------------------------------------------------------------------
bool Simulator::setUpTopology_()
{
  std::string netListFile;

  if( procZeroFlag_ )
  {
    if (commandLine.getArgumentValue("netlist") != "")
      netListFile = commandLine.getArgumentValue("netlist");
    else
      netListFile = "xyce.in";

    FILE * testFile;
    if ( (testFile = fopen(netListFile.c_str(), "r")) == 0)
    {
      Xyce::Report::UserError() << netListFile << " file not found";
    }
    else
    {
      fclose( testFile );
    }
  }

  Xyce::Report::safeBarrier(comm_);

  Xyce::lout() << "***** Reading and parsing netlist..." << std::endl;

  netlistImportToolPtr_->constructCircuitFromNetlist(netListFile, externalNetlistParams_);

  Xyce::Report::safeBarrier(comm_);

  if ( commandLine.argExists("-syntax") )
  {
    Xyce::lout() << "***** Netlist syntax OK\n" << std::endl
                 << "***** Device Type Counts ...\n" << std::endl;

    Xyce::IO::printGlobalDeviceCounts(Xyce::lout(), comm_, devIntPtr_->getDeviceCountMap());
    if (outMgrPtr_->getDetailedDeviceFlag())
      Xyce::IO::printLocalDeviceCount(Xyce::lout(), comm_, devIntPtr_->getDeviceCountMap());

    Xyce::lout() << std::endl;

    reportTotalElapsedTime ();

    Xyce_exit(0);
  }

  delete netlistImportToolPtr_;

  Xyce::lout() << "***** Setting up topology...\n" << std::endl;

  // topology query's device manager to see if any devices are bad (i.e. a resistor with zero resistance)
  // if so, a list of nodes to be supernoded is created
  topPtr_->verifyNodesAndDevices();

#ifdef Xyce_PARALLEL_MPI
  // create a union of the supernode list on all processors
  topPtr_->mergeOffProcTaggedNodesAndDevices();
#endif

// combine nodes into supernodes and remove now redundant devices (i.e. those only connected to 1 processor )
  topPtr_->removeTaggedNodesAndDevices();

  // if "-remeasure" was on the command line, then we don't need to
  // instantiate the devices.
  if (commandLine.argExists("-remeasure"))
  {
    outMgrPtr_->remeasure();
    Xyce::lout() << "***** Remeasure analysis complete\n" << std::endl;
    return false;
  }
  topPtr_->instantiateDevices();

  outMgrPtr_->delayedPrintLineDiagnostics();

  Xyce::Report::safeBarrier(comm_);

#ifdef Xyce_PARALLEL_MPI
  devIntPtr_->setGlobalFlags();
#endif

  // Setup of indices including global reordering.
  topPtr_->setupGlobalIndices();

  Xyce::Report::safeBarrier(comm_);

#ifdef Xyce_DEBUG_DEVICE
  for (Xyce::Device::EntityTypeIdDeviceMap::const_iterator it = devIntPtr_->getDeviceMap().begin(); it != devIntPtr_->getDeviceMap().end(); ++it)
    print(Xyce::lout(), *(*it).second);
#endif

  return true;

}

//-----------------------------------------------------------------------------
// Function      : Simulator::setUpMatrixStructure_
// Purpose       : This function needs to set up the various linear algebra
//                 entities which are owned by the LAS system class.  This
//                 includes, most importantly, the Jacobian matrix and the
//                 right hand side vector.  It should also set the solution
//                 vector size and the state vector size.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/00
//-----------------------------------------------------------------------------
bool Simulator::setUpMatrixStructure_()
{
  lasBuilderPtr_->generateParMaps();
  lasBuilderPtr_->generateGraphs();

  lasSysPtr_->initializeSystem();

  topPtr_->registerLIDswithDevs();
  topPtr_->registerJacLIDswithDevs();

#ifdef Xyce_EXTDEV
  devIntPtr_->setupExternalDevices();
#endif

  int lasSize = lasSysPtr_->getGlobalSolutionSize();
  Xyce::lout() << "***** Number of Unknowns = " << lasSize << std::endl;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::doInitializations_
// Purpose       : This function calls "initializeAll" functions in all the
//                 packages which have them.  Packages that have the LAS system
//                 class registered with them will have one of these functions.
//
// Special Notes : This is called once the size of the linear system is known,
//                 as various packages will need to allocate vectors, matrices,
//                 etc.
//
//                 These probably can be done in any order, but to be safe,
//                 make sure for now that the time integrator's initializeAll
//                 function is called first.  The other two are essentially
//                 secondary registrations, while the TIA one includes a lot of
//                 allocations.
//
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool Simulator::doInitializations_()
{
  bool bsuccess = true;
  bool bs1 = true;

  bs1 = anaIntPtr_->initializeAll();  bsuccess = bsuccess && bs1;
  bs1 = devIntPtr_->initializeAll();  bsuccess = bsuccess && bs1;
  bs1 = nlsMgrPtr_->initializeAll();  bsuccess = bsuccess && bs1;

  if( resMgrPtr_->isRestart() ) resMgrPtr_->restoreRestartData();

  topPtr_->generateICLoader();

  initializeAllFlag_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::runSolvers_
// Purpose       : This function runs the solvers.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool Simulator::runSolvers_()
{
  return anaIntPtr_->run();
}


//-----------------------------------------------------------------------------
// Function      : Simulator::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
bool Simulator::run(int iargs_tmp, char *cargs_tmp[])
{
  bool bsuccess = true;

  try 
  {
    bool bs1 = initialize(iargs_tmp, cargs_tmp);
    bsuccess = bsuccess && bs1;
  }
  catch (std::exception &x) 
  {
    Xyce::lout() << "Exception " << x.what() << std::endl;
    bsuccess = false;
  }

  if (!bsuccess)
  {
    reportTotalElapsedTime ();
    Xyce::lout() << "Xyce Initialization Phase failed." << std::endl;
    Xyce_exit(-1);
  }

  try {
    bool bs1 = runSimulation();
    bsuccess = bsuccess && bs1;

    bs1 = finalize();
    bsuccess = bsuccess && bs1;
  }
  catch (std::exception &x) {
    Xyce::lout() << "Exception " << x.what() << std::endl;
    bsuccess = false;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Simulator::runSimulation
// Purpose       : Main simulation driver.
// Special Notes : Not private as this is also called from N_DAK_DakotaInterface
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool Simulator::runSimulation()
{
  return runSolvers_();
}

// ---------------------------------------------------------------------------
// API METHODS NEEDED FOR MIXED-SIGNAL and other external applications
//
//-----------------------------------------------------------------------------
// Function      : Simulator::initialize
// Purpose       : capture all "initialization-type" activities in one
//                 method
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 02/17/2004
//-----------------------------------------------------------------------------
bool Simulator::initialize( int iargs_tmp, char **cargs_tmp)
{
  iargs = iargs_tmp;
  cargs = cargs_tmp;

  // Setup the Parallel Manager first to give us a N_PDS_Comm object so the
  // reporting package will work.
  bool bsuccess = setupParMgr_( iargs, cargs );
  bool b2;

  // Start the solvers timer.
  ElapsedTimerPtr_ = new N_UTL_Timer( *(parMgrPtr_->getPDSComm()) );

#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_DEBUG_ALL_PROCS_SAME_WD
    N_PDS_Comm * commPtr  = parMgrPtr_->getPDSComm();
    size_t dirsize = 512;
    char directory[dirsize];  for (int idir=0;idir<dirsize;++idir)  directory[idir] = 0;
#endif
#endif

  // register parallel mgr to allow parallel support
    commandLine.registerParallelMgr( parMgrPtr_ );

    Xyce::initializeLogStream(parMgrPtr_->getPDSComm()->procID(), parMgrPtr_->getPDSComm()->numProc());

  // read in command line arguments
    int status = commandLine.parseCommandLine( iargs, cargs );
    if (status != 0)
      Xyce_exit(status == -1 ? 0 : status);

  // Set the output stream of the "-l" flag exists
  if (commandLine.argExists("-l"))
  {
    Xyce::openLogFile(commandLine.getArgumentValue("-l"), commandLine.argExists("-per-processor"));
  }

  if( procZeroFlag_ )
  {
    Xyce::lout() << "\n"
                 << "*****\n"
                 << "***** Welcome to the Xyce(TM) Parallel Electronic Simulator\n"
                 << "*****\n"
                 << "***** This is version " << N_UTL_Version::getFullVersionString() << "\n\n\n"
                 << "***** Executing netlist " << (commandLine.getArgumentValue("netlist") == "" ? "MISSING!" : commandLine.getArgumentValue("netlist")) << "\n\n";

    // Don't bother doing anything else if we weren't given a netlist!
    if (commandLine.getArgumentValue("netlist") == "") 
    {
      Xyce::lout() << "  Usage: " << cargs[0] << " netlist\n" << std::endl;
      return false;
    }

    // check if a parameter file was specified for param substitution during parsing
    if (commandLine.argExists(std::string("-prf")))
    {
      std::string parameterFile = commandLine.getArgumentValue("-prf");
      readExternalParamsFromFile( parameterFile, iterationSuffix_, externalNetlistParams_ );
#ifdef Xyce_DEBUG
      // debug print out externNetlistParams_
      std::vector<std::pair<std::string,std::string> >::iterator currentPair = externalNetlistParams_.begin();
      std::vector<std::pair<std::string,std::string> >::iterator endPair = externalNetlistParams_.end();
      while( currentPair != endPair )
      {
        Xyce::dout() << "\"" << currentPair->first << "\" = \"" << currentPair->second << "\"" << std::endl;
        currentPair++;
      }
#endif
    }

  }

#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_DEBUG_ALL_PROCS_SAME_WD

  if( procZeroFlag_ )
  {
    (void) getcwd(directory, dirsize);
  }
  commPtr->bcast(directory, dirsize, 0);
  if( !procZeroFlag_ )
  {
    chdir(directory);
  }

#endif
#endif

  // Allocate all the various packages:

  b2 = doAllocations_();
  bsuccess = bsuccess && b2;

  if (b2) {
    if (Xyce::DEBUG_CIRCUIT)
      Xyce::dout() << "Allocation was successful.";
  }
  else
    Xyce::Report::DevelFatal() << "Allocation was NOT successful.";

  // Register the external package pointers:
  b2 = doRegistrations_();
  bsuccess = bsuccess && b2;

  if (b2) {
    if (Xyce::DEBUG_CIRCUIT)
      Xyce::dout() << "Registration was successful.";
  }
  else
    Xyce::Report::DevelFatal() << "Registration was NOT successful.";

  b2 = setUpTopology_();
  if( !b2 )
  {
    return false;
  }
  bsuccess = bsuccess && b2;

  Xyce::lout() << "***** Device Count Summary ..." << std::endl;

  Xyce::IO::printGlobalDeviceCounts(Xyce::lout(), comm_, devIntPtr_->getDeviceCountMap());
  if (outMgrPtr_->getDetailedDeviceFlag())
    Xyce::IO::printLocalDeviceCount(Xyce::lout(), comm_, devIntPtr_->getDeviceCountMap());

  if( commandLine.argExists( "-norun" ) || 
      commandLine.argExists( "-namesfile" )  )
  {
    Xyce::lout() << "\n***** Syntax and topology analysis complete" << std::endl;

    if( commandLine.argExists( "-namesfile" ) )
    {
      setUpMatrixStructure_();
      topPtr_->outputNameFile(true);
    }

    return false;
  }
  else
  {
    Xyce::lout() << "\n***** Setting up matrix structure..." << std::endl;
    b2 = setUpMatrixStructure_();
    bsuccess = bsuccess && b2;

    Xyce::lout() << "***** Initializing...\n" << std::endl;
    b2 = doInitializations_();
    bsuccess = bsuccess && b2;

    // optional diagnostic output file:
    topPtr_->outputNameFile();

    // if we loaded a parameter file from the command line, then the output manager
    // should scan the results as well as such files can specify response functions
    // than need to be reported by the output manager.
    if (commandLine.argExists(std::string("-prf")))
    {
      outMgrPtr_->setExternalNetlistParams( externalNetlistParams_ );
      if(iterationSuffix_.length() > 0)
      {
        outMgrPtr_->setOutputFilenameSuffix( iterationSuffix_ );
      }
    }
    if (commandLine.argExists(std::string("-rsf")))
    {
      std::string responseFile = commandLine.getArgumentValue("-rsf");
      outMgrPtr_->setResponseFilename( responseFile );
    }
    // Start the solvers timer.
    XyceTimerPtr_ = new N_UTL_Timer( *(parMgrPtr_->getPDSComm()) );
  }

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : Simulator::getDACDeviceNames
// Purpose        : Gets the (stripped) names of the DAC devices
//                  in the circuit.
// Special Notes  :
// Scope          :
// Creator        : Lisa Maynes
// Creation Date  : 06/13/2003
//----------------------------------------------------------------------------
bool Simulator::getDACDeviceNames(std::vector< std::string >& dacNames)
{
  bool bsuccess = true;
  bsuccess = devIntPtr_ -> getDACDeviceNames( dacNames );
  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : getADCMap
// Purpose        : Gets the (stripped) names of the ADC devices
//                 in the circuit(as key of map) and map of parameters
//                 (keyed by parameter name) for each device
// Special Notes  :
// Scope          :
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Simulator::getADCMap(std::map<std::string, std::map<std::string,double> >&ADCMap)
{
  bool bsuccess = true;
  bsuccess = devIntPtr_ -> getADCMap(ADCMap);
  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : updateTimeVoltagePairs
// Purpose        : Update the DAC devices in a circuit by adding the set
//                  of time and voltage pairs built up on the "digital side"
//                  since the last update and by removing the time-voltage
//                  pairs for times that pre-date the given simulation time.
// Special Notes  : The current method for locating DAC instances
//                  works for the serial case only. The parallel
//                  case will be added in later modifications.
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 06/09/2003
//----------------------------------------------------------------------------
bool Simulator::updateTimeVoltagePairs(
   std::map< std::string, std::vector< std::pair<double,double> >* > const & timeVoltageUpdateMap)
{
  bool bsuccess = true;

  bsuccess=devIntPtr_->updateTimeVoltagePairs(timeVoltageUpdateMap);

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : getTimeVoltagePairs
// Purpose        : get a map of all time-voltage pairs from all ADC instances
//
// Special Notes  : The current method for locating ADC instances
//                  works for the serial case only. The parallel
//                  case will be added in later modifications.
// Scope          : public
// Creator        : Tom Russo
// Creation Date  : 05/10/2004
//----------------------------------------------------------------------------
bool Simulator::getTimeVoltagePairs(
   std::map< std::string, std::vector< std::pair<double,double> > > & timeVoltageUpdateMap)
{
  bool bsuccess = true;

  bsuccess=devIntPtr_->getTimeVoltagePairs(timeVoltageUpdateMap);

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : setADCWidths
// Purpose        : Update the ADC devices in a circuit by informing them
//                  of the width of their bitvector output on the
//                  "digital side"
// Special Notes  :
// Scope          :
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Simulator::setADCWidths(
   std::map< std::string, int > const & ADCWidthMap)
{
  bool bsuccess = true;

  bsuccess=devIntPtr_->setADCWidths(ADCWidthMap);

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::simulateUntil
// Purpose       : To continue the existing analog circuit simulation
//                 until either the given <requestedUntilTime> is reached
//                 or the simulation termination criterion is met.
//                 Return a Boolean indicating whether the simulation
//                 run was successful. (Note that the run is successful
//                 even when the given <requestedUntilTime> is not reached,
//                 so long as the run completed normally.)
// Special Notes : The time variables are in units of seconds.
// Scope         : public
// Creator       : Lon Waters, SNL
//               : Tom Russo, SNL, Component Information and Models
// Creation Date : 06/03/2003
//---------------------------------------------------------------------------
bool Simulator::simulateUntil(double requestedUntilTime,
                               double& completedUntilTime)
{
  bool bsuccess = false;
  double currentTimeBeforeSim = anaIntPtr_->getTime();
  double finalTime=anaIntPtr_->getFinalTime();
  double initialTime = anaIntPtr_->getInitialTime();

  // Silence the "Percent Complete" noise, we don't want it when using this
  // interface

  anaIntPtr_->silenceProgress();

#ifdef Xyce_DEBUG_CIRCUIT
  Xyce::dout() << "Xyce::simulateUntil: ";
  Xyce::dout() << "finalTime = " << finalTime
            << ", currentTimeBeforeSim = " << currentTimeBeforeSim << std::endl;
#endif

  if (currentTimeBeforeSim >= finalTime)
  {
    // We have already simulated as far as the netlist requested.
    bsuccess = true;
    completedUntilTime = currentTimeBeforeSim;
#ifdef Xyce_DEBUG_CIRCUIT
    Xyce::dout() << "Case1: completedUntilTime = " << completedUntilTime;
#endif
  }
  else
  {
    // We are not already past the end of the netlist time
    anaIntPtr_->setPauseTime(Xycemin(requestedUntilTime,finalTime));
#ifdef Xyce_DEBUG_CIRCUIT
    Xyce::dout() << "Xyce::simulateUntil currentTimeBeforeSim = " << currentTimeBeforeSim << "  initialTime = " << initialTime << std::endl;
#endif
    if (currentTimeBeforeSim > initialTime)
    {
      anaIntPtr_->resumeSimulation();
    }

#ifdef Xyce_DEBUG_CIRCUIT
    Xyce::dout() << "Xyce::simulateUntil: Case2: requestedUntilTime = " << requestedUntilTime
              << ", pauseTime = " << anaIntPtr_->getPauseTime() << std::endl;
#endif

    bsuccess = runSolvers_();
    completedUntilTime = anaIntPtr_->getTime();
#ifdef Xyce_DEBUG_CIRCUIT
    Xyce::dout() << "Xyce::simulateUntil: Case2: completedUntilTime = " << completedUntilTime << std::endl;
#endif
  }

#ifdef Xyce_DEBUG_CIRCUIT
  Xyce::dout() << std::endl;
#endif
  //  return true;
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::finalize
// Purpose       : To clean up after driving Xyce with the SIMBUS
//                 simulation backplane. This includes the following:
//                    Free any dynamically allocated memory...
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
//               : Lisa Maynes, CoMeT Solutions, Inc.
// Creation Date : 06/03/2003
//---------------------------------------------------------------------------
bool Simulator::finalize()
{
  bool bsuccess = true;

  Xyce::lout() << "\n***** Solution Summary *****"  << std::endl;

  anaIntPtr_->outputSummary();
  outMgrPtr_->outputMacroResults();

  {
    Xyce::lout() << std::endl
                 << "***** Total Simulation Solvers Run Time: " << XyceTimerPtr_->elapsedTime() << " seconds" << std::endl
                 << "***** Total Elapsed Run Time:            " << ElapsedTimerPtr_->elapsedTime() << " seconds" << std::endl
                 << "*****" << std::endl
                 << "***** End of Xyce(TM) Simulation" << std::endl
                 << "*****" << std::endl;
  }

  // Close the output stream:
  Xyce::closeLogFile();

  bsuccess = doDeAllocations_();

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::reportTotalElapsedTime ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/01/2009
//---------------------------------------------------------------------------
void Simulator::reportTotalElapsedTime ()
{
  Xyce::lout() <<  "\n***** Total Elapsed Run Time: " << ElapsedTimerPtr_->elapsedTime() << " seconds" << std::endl;
}

//---------------------------------------------------------------------------
// Function      : Simulator::simulationComplete
// Purpose       : Simply report whether we've reached the end of the
//                 simulation
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//---------------------------------------------------------------------------
bool Simulator::simulationComplete()
{
  return anaIntPtr_->simulationComplete();
}

//
// new mixed-signal functions:
// These are provisional!

//---------------------------------------------------------------------------
// Function      : Simulator::provisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/2009
//---------------------------------------------------------------------------
bool Simulator::provisionalStep
  (double maxTimeStep,
   double &timeStep,
   std::map< std::string, std::vector< std::pair<double,double> > > & timeVoltageUpdateMap)
{
  bool bsuccess=true;

  bool b1 = anaIntPtr_->provisionalStep(maxTimeStep, timeStep);
  bsuccess = bsuccess && b1;

  b1=getTimeVoltagePairs(timeVoltageUpdateMap);

  bsuccess = bsuccess && b1;

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::getFinalTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/25/2009
//---------------------------------------------------------------------------
double Simulator::getFinalTime()
{
  double ft=0.0;
  if (anaIntPtr_!=0)
  {
    ft = anaIntPtr_->getFinalTime();
  }
  return ft;
}

//---------------------------------------------------------------------------
// Function      : Simulator::getTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//---------------------------------------------------------------------------
double Simulator::getTime()
{
  double t=0.0;
  if (anaIntPtr_!=0)
  {
    t = anaIntPtr_->getTime();
  }
  return t;
}

//---------------------------------------------------------------------------
// Function      : Simulator::acceptProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//---------------------------------------------------------------------------
void Simulator::acceptProvisionalStep()
{
  anaIntPtr_->acceptProvisionalStep ();
}

//---------------------------------------------------------------------------
// Function      : Simulator::rejectProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/2009
//---------------------------------------------------------------------------
void Simulator::rejectProvisionalStep()
{
  anaIntPtr_->rejectProvisionalStep();
}

// ---------------------------------------------------------------------------
// API METHODS NEEDED FOR Two-level Functions:

//---------------------------------------------------------------------------
// Function      : Simulator::simulateStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/2006
//---------------------------------------------------------------------------
bool Simulator::simulateStep
      ( const N_DEV_SolverState & solState,
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError)
{
  bool bsuccess = false;

#ifdef Xyce_DEBUG_CIRCUIT
  Xyce::dout() << "\nXyce::simulateStep: " << std::endl;
#endif

  // Apply the input voltages to their appropriate sources:
  std::map<std::string,double>::const_iterator iterM = inputMap.begin();
  std::map<std::string,double>::const_iterator  endM = inputMap.end  ();
  int i=0;
  for (; iterM != endM;++iterM, ++i)
  {
    bool found = true;
    std::string name = iterM->first;
    double val  = iterM->second;
    devIntPtr_->setParam (name,val);
  }

  cktLoaderPtr_->setExternalSolverState (solState);

  bsuccess = anaIntPtr_->runStep (solState.tiInfo, tlError);

  // calculate the conductance:
  nlsMgrPtr_->obtainConductances(
        inputMap,
        outputVector,
        jacobian
    );

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::simulateStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 07/16/2009
//---------------------------------------------------------------------------
bool Simulator::simulateStep
      (const N_DEV_ExternalSimulationData & ext_data,
       const std::map<std::string,double> & inputMap,
       std::vector<double> & outputVector,
       std::vector< std::vector<double> > & jacobian,
       N_TIA_TwoLevelError & tlError)
{
  bool bsuccess = false;

#ifdef Xyce_DEBUG_CIRCUIT
  Xyce::dout() << "\nXyce::simulateStep: " << std::endl;
#endif

  // Create N_DEX_SolverState object
  N_DEV_SolverState state;

  if (ext_data.is_transient)
  {
    Xyce::lout() << "Xychron (mixed-model): Transient Xyce Step" << std::endl;

    state.currTimeStep = ext_data.current_time_step_size;
    state.lastTimeStep = ext_data.previous_time_step_size;
    state.currTime = ext_data.current_time;
    state.finalTime = ext_data.final_time;
    state.timeStepNumber = ext_data.time_step_number;
    state.startingTimeStep = ext_data.current_time_step_size;
    state.tranopFlag = false;
    state.transientFlag = true;
    state.dcopFlag = false;
    state.dcsweepFlag = false;
    // tiInfo
    state.tiInfo.nextTimeStep = ext_data.current_time_step_size;
    state.tiInfo.currTimeStep = ext_data.current_time_step_size;
    state.tiInfo.nextTime = ext_data.current_time;
    state.tiInfo.finalTime = ext_data.final_time;
    state.tiInfo.dcopFlag = false;
    state.tiInfo.tranopFlag = false;
    state.tiInfo.transientFlag = true;
    state.tiInfo.dcsweepFlag = false;
    state.tiInfo.timeStepNumber = ext_data.time_step_number;
    state.tiInfo.timeIntMode = 6;

    // New changes
    // pdt = 1/dt
    state.pdt = 1.0/ext_data.current_time_step_size;
    // bptol - ignore for now
    state.doubleDCOPEnabled = false;
    state.doubleDCOPStep = 1;
    // newtonIteration - doens't matter for particular devices.  For volt lim

    if (state.timeStepNumber == 0)
    {
      state.initTranFlag = true;
      state.tiInfo.initTranFlag = true;
      state.tiInfo.beginIntegrationFlag = true;
    }
    else
    {
      state.initTranFlag = false;
      state.tiInfo.initTranFlag = false;
      state.tiInfo.beginIntegrationFlag = false;
    }

    // tia pdt - set this
    state.tiInfo.pdt = state.pdt;
    state.tiInfo.currentOrder = 1;
  }
  else
  {
    state.tranopFlag = true;  // for a DCOP only, this should be false
    state.transientFlag = false;
    state.dcopFlag = false;
    state.dcsweepFlag = false;
    state.tiInfo.dcopFlag = true;
    state.tiInfo.tranopFlag = true;
    state.tiInfo.transientFlag = false;
    state.tiInfo.dcsweepFlag = false;
    state.tiInfo.timeIntMode = 0;
  }

  bsuccess = this->simulateStep(state,inputMap,outputVector,jacobian,tlError);

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::startupSolvers
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/2006
//---------------------------------------------------------------------------
bool Simulator::startupSolvers ()
{
  bool bsuccess = true;
  bsuccess = anaIntPtr_->startupSolvers ();
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::finishSolvers
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/2006
//---------------------------------------------------------------------------
bool Simulator::finishSolvers ()
{
  bool bsuccess = true;
  bsuccess = anaIntPtr_->finishSolvers ();
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::homotopyStepSuccess
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//---------------------------------------------------------------------------
void Simulator::homotopyStepSuccess
    (const std::vector<std::string> & paramNames,
     const std::vector<double> & paramVals)
{
  anaIntPtr_->homotopyStepSuccess ( paramNames, paramVals);
  return;
}

//---------------------------------------------------------------------------
// Function      : Simulator::homotopyStepFailure
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/2006
//---------------------------------------------------------------------------
void Simulator::homotopyStepFailure ()
{
  anaIntPtr_->homotopyStepFailure ();
  return;
}

//---------------------------------------------------------------------------
// Function      : Simulator::stepSuccess
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
void Simulator::stepSuccess (int analysis)
{
  anaIntPtr_->stepSuccess (analysis);
  return;
}

//---------------------------------------------------------------------------
// Function      : Simulator::stepFailure
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
void Simulator::stepFailure (int analysis)
{
  anaIntPtr_->stepFailure (analysis);
  return;
}


//---------------------------------------------------------------------------
// Function      : Simulator::getInitialQnorm
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
bool Simulator::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  bool bsuccess = true;

  if (initializeAllFlag_)
  {
    bsuccess = anaIntPtr_->getInitialQnorm (tle);
  }

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::getBreakPoints
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
bool Simulator::getBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes)
{
  bool bsuccess = true;

  if (initializeAllFlag_)
  {
    bsuccess = anaIntPtr_->getBreakPoints (breakPointTimes);
  }

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::updateStateArrays
// Purpose       :
// Special Notes : Used for 2-level Newton solves, with LOCA.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/18/2006
//---------------------------------------------------------------------------
bool Simulator::updateStateArrays ()
{
  bool bsuccess = true;
#if 0
  if (initializeAllFlag_)
  {
    bsuccess = anaIntPtr_->updateStateArrays ();
  }
#endif
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::setInternalParam
// Purpose       :
// Special Notes : Used for 2-level Newton solves, with LOCA.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/18/2006
//---------------------------------------------------------------------------
bool Simulator::setInternalParam (std::string & name, double val)
{
  cktLoaderPtr_->setParam (name, val);
  return true;
}

//---------------------------------------------------------------------------
// Function      : Simulator::startTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//---------------------------------------------------------------------------
bool Simulator::startTimeStep (const N_TIA_TimeIntInfo & tiInfo)
{
  bool bsuccess = anaIntPtr_->startTimeStep(tiInfo);
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::startTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 08/28/2009
//---------------------------------------------------------------------------
bool Simulator::startTimeStep (const N_DEV_ExternalSimulationData & ext_data)
{
  N_TIA_TimeIntInfo tiInfo;

  if (ext_data.is_transient)
  {
    tiInfo.nextTimeStep = ext_data.current_time_step_size;
    tiInfo.currTimeStep = ext_data.current_time_step_size;
    tiInfo.nextTime = ext_data.current_time;
    tiInfo.finalTime = ext_data.final_time;
    tiInfo.dcopFlag = false;
    tiInfo.tranopFlag = false;
    tiInfo.transientFlag = true;
    tiInfo.dcsweepFlag = false;
    tiInfo.timeStepNumber = ext_data.time_step_number;
    tiInfo.timeIntMode = 6;

    if (ext_data.time_step_number == 0)
    {
      tiInfo.initTranFlag = true;
      tiInfo.beginIntegrationFlag = true;
    }
    else
    {
      tiInfo.initTranFlag = false;
      tiInfo.beginIntegrationFlag = false;
    }

    // tia pdt - set this
    tiInfo.pdt = 1.0/ext_data.current_time_step_size;
    tiInfo.currentOrder = 1;
  }
  else
  {
    tiInfo.dcopFlag = true;
    tiInfo.tranopFlag = true;
    tiInfo.transientFlag = false;
    tiInfo.dcsweepFlag = false;
    tiInfo.timeIntMode = 0;
  }

  bool bsuccess = anaIntPtr_->startTimeStep(tiInfo);
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::endTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Russell Hooper, SNL
// Creation Date : 10/16/2012
//---------------------------------------------------------------------------
bool Simulator::endTimeStep (N_DEV_ExternalSimulationData & ext_data)
{
  // We could opt to obtain selected time integration data here or
  // allow it to be obtained higher up, eg in charon::sc::CircuitDriver
  bool bsuccess = true;
  N_TIA_TimeIntInfo tiInfo;
  anaIntPtr_->getTimeIntInfo(tiInfo);

  // Copy N_TIA_TimeIntInfo into our ext_data container
  ext_data.currentOrder         = tiInfo.currentOrder        ;
  ext_data.numberOfSteps        = tiInfo.numberOfSteps       ;
  ext_data.usedOrder            = tiInfo.usedOrder           ;
  ext_data.nscsco               = tiInfo.nscsco              ;
  ext_data.pdt                  = tiInfo.pdt                 ;
  ext_data.nextTimeStep         = tiInfo.nextTimeStep        ;
  ext_data.currTimeStep         = tiInfo.currTimeStep        ;
  ext_data.currentTime          = tiInfo.currentTime         ;
  ext_data.nextTime             = tiInfo.nextTime            ;
  ext_data.finalTime            = tiInfo.finalTime           ;
  ext_data.startingTimeStep     = tiInfo.startingTimeStep    ;
  ext_data.bpTol                = tiInfo.bpTol               ;
  ext_data.dcopFlag             = tiInfo.dcopFlag            ;
  ext_data.acopFlag             = tiInfo.acopFlag            ;
  ext_data.inputOPFlag          = tiInfo.inputOPFlag         ;
  ext_data.tranopFlag           = tiInfo.tranopFlag          ;
  ext_data.transientFlag        = tiInfo.transientFlag       ;
  ext_data.dcsweepFlag          = tiInfo.dcsweepFlag         ;
  ext_data.timeStepNumber       = tiInfo.timeStepNumber      ;
  ext_data.initTranFlag         = tiInfo.initTranFlag        ;
  ext_data.beginIntegrationFlag = tiInfo.beginIntegrationFlag;
  ext_data.doubleDCOPStep       = tiInfo.doubleDCOPStep      ;
  ext_data.doubleDCOPEnabled    = tiInfo.doubleDCOPEnabled   ;
  ext_data.stepLoopIter         = tiInfo.stepLoopIter        ;
  ext_data.timeIntMode          = tiInfo.timeIntMode         ;
  ext_data.sweepSourceResetFlag = tiInfo.sweepSourceResetFlag;

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : Simulator::endTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Russell Hooper, SNL
// Creation Date : 10/16/2012
//---------------------------------------------------------------------------
void Simulator::enable_lte_analysis()
{
  N_TIA_TIAParams & tiaParams = anaIntPtr_->getTIAParams();
  tiaParams.errorAnalysisOption = 0; // use local truncation error estimates
}

//---------------------------------------------------------------------------
// Function      : Simulator::readExternalParamsFromFile
// Purpose       :
// Special Notes : Used to read parameters "tag" = "value" from an external
//                 file.  Any "tag"s found while reading the netlist during
//                 parsing will be replaced by "value".
// Scope         : public
// Creator       : Richard Schiek, 1437 Electrical Systems Modeling
// Creation Date : 07/24/2012
//---------------------------------------------------------------------------
void Simulator::readExternalParamsFromFile( std::string filename, std::string & iterationSuffix,
  std::vector< std::pair< std::string, std::string > > & paramList )
{
  // attempt to get parameter file suffix so that if it is  part of a unique identifier
  // like an iteration or realization number it can be attached to the response output
  // file as well. (i.e. params.in.xxxx)  For now I'll assume that the last ".xxx" is
  // significant.  if we just get ".in" or ".dat" then that will be reflected on the output
  // file.

  std::string::size_type suffixStart=filename.find_last_of( "." );
  if( suffixStart != std::string::npos )
  {
    iterationSuffix.assign( filename, suffixStart, filename.length()-suffixStart );
    // check for trivial results of ".dat" ".txt" and ".in"
    if( iterationSuffix==".dat" || iterationSuffix==".txt" || iterationSuffix==".in" )
      iterationSuffix.erase();
  }


  // at this stage just support the Dakota params.in format of "value" = "tag".
  // we could support other formats as well.

  const std::string allowedChars("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_.$");
  const std::string whiteSpace(" \t=\n\r");    // note we will treat "=" as whitespace
  const std::string commentChars("*#;");  // if we find any of these then the rest of the line is a comment

  // attempt to open the params file
  std::ifstream paramFile(filename.c_str(), std::ios::in);
  if( paramFile )
  {
    // loop over the file trying to gather names / values and response functions required.
    // general format:
    // white space tag/value white space or = tag/value
    //
    std::string aLine;
    getline(paramFile, aLine);
    while ( paramFile.good() )
    {
      // getline worked, so try to parse the line
      // procedure (1) discard any comments
      //           (2) find tag/value boundaries by looking for allowedChars followed by whitespace
      //           (3) sort out order of tag/value pair
      //           (4) store pair in paramList
      //           (5) try to get another line

      // don't bother with lines that are blank or essentially blank (i.e. just \r\n)
      // shortest line we could have is x=1 or 3 chars.
      if( aLine.length() > 2 )
      {
        std::string::size_type commentLoc = aLine.find_first_of( commentChars, 0 );
        if( commentLoc != std::string::npos )
        {
          // chop off comment to end of line.
          aLine.erase( commentLoc, aLine.length()-commentLoc );
        }
        //check overall line length again.  This could have been just a comment line
        if( aLine.length() > 2 )
        {
          std::string::size_type word1Start = aLine.find_first_of( allowedChars, 0 );
          // check that we found a valid word1Start otherwise stop trying to parse this line.
          if( word1Start != std::string::npos )
          {
            std::string::size_type word1End   = aLine.find_first_of( whiteSpace, word1Start );
            // check that we found a valid word1End otherwise stop trying to parse this line.
            if( word1End != std::string::npos )
            {
              std::string::size_type word2Start = aLine.find_first_of( allowedChars, word1End );
              // check that we found a valid word2Start otherwise stop trying to parse this line.
              if( word2Start != std::string::npos )
              {
                std::string::size_type word2End   = aLine.find_first_of( whiteSpace, word2Start );
                // check that we found a valid word2End
                if( word2End == std::string::npos )
                  word2End = aLine.length();
                // if we get here then we have valid start,end indicies for word1 and word2

                std::string word1=aLine.substr( word1Start, (word1End - word1Start) );
                std::string word2=aLine.substr( word2Start, (word2End - word2Start) );

                // sort out tag/value ordering.
                // if word1=number assume format is value = tag
                // otherwise assume format is tag = value
                std::stringstream converter;
                converter << word1;
                double testvalue;
                converter >> testvalue;
                if( converter.fail() )
                {
                  // couldn't convert tag1 to a double value so assume format is word1=tag word2=value
                  paramList.push_back( std::pair<std::string,std::string>(word1,word2) );
                }
                else
                {
                  // tag1 was successfully converted to a value so assume format is word1=value word2=tag
                  paramList.push_back( std::pair<std::string,std::string>(word2,word1) );
                }

              }  // if ( word2Start != std::string::npos )
            }    // if( word1End != std::string::npos )
          }      // if( word1Start != std::string::npos )
        }        // if( aLine.length() > 2 ) -- check after removing comments
      }          // if( aLine.lenght() > 2 ) -- outer check
      // try to get another line
      getline(paramFile, aLine);
    }
  }
  else
  {
    Xyce::Report::UserWarning() << "Could not open parameter file: " + filename + ". Attempting to continue.";
  }

  // for debug purposes.  output the params as read
  Xyce::dout() << "Parameters read from \"" << filename << "\"" << std::endl;
  std::vector< std::pair< std::string, std::string > >::iterator listitr = paramList.begin();
  std::vector< std::pair< std::string, std::string > >::iterator enditr = paramList.end();
  while( listitr != enditr )
  {
    Xyce::dout() << "  " << listitr->first << " , " << listitr->second << std::endl;
    listitr++;
  }
}

//---------------------------------------------------------------------------
// Function      : Simulator::registerResponseVars
// Purpose       :
// Special Notes : Used when Dakota or an external program calls Xyce to tell
//                 Xyce what variables and expressions should be held as
//                 response funcctions to pass back to Dakota.  This
//                 is just pass through to the output Manager
// Scope         : public
// Creator       : Richard Schiek, 1437 Electrical and MEMS modeling
// Creation Date : 10/22/2008
//---------------------------------------------------------------------------
bool Simulator::registerResponseVars (std::string objString, RCP<std::vector< double > > varVectorPtr )
{
  if( outMgrPtr_ == 0)
  {
    Xyce::Report::DevelFatal0() << "Xyce::registerResponseVars outMgrPtr_ is null";
  }
  bool result = outMgrPtr_->registerResponseVars( objString, varVectorPtr );

  return result;

}

//---------------------------------------------------------------------------
// Function      : Simulator::registerResponseVars
// Purpose       :
// Special Notes : Used when Dakota or an external program calls Xyce to tell
//                 Xyce what variables and expressions should be held as
//                 response funcctions to pass back at end of simulation.
//                 This is just pass through to the output Manager
// Scope         : public
// Creator       : Richard Schiek, 1437 Electrical and MEMS modeling
// Creation Date : 10/22/2008
//---------------------------------------------------------------------------
void Simulator::finalizeResponseVars()
{
  if( outMgrPtr_ == 0)
  {
    Xyce::Report::DevelFatal0() << "Xyce::finalizeResponseVars outMgrPtr_ is null";
  }
  outMgrPtr_->finalizeResponseVars();
}

//---------------------------------------------------------------------------
// Function      : Simulator::initializeTransientModel
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
void Simulator::initializeTransientModel()
{
  anaIntPtr_->initializeTransientModel();
}

//---------------------------------------------------------------------------
// Function      : Simulator::evalTransientModel
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
bool Simulator::evalTransientModel
    (
     double t,
     N_LAS_Vector * SolVectorPtr,
     N_LAS_Vector * CurrSolVectorPtr,
     N_LAS_Vector * LasSolVectorPtr,
     N_LAS_Vector * StaVectorPtr,
     N_LAS_Vector * CurrStaVectorPtr,
     N_LAS_Vector * LasStaVectorPtr,
     N_LAS_Vector * StaDerivVectorPtr,
     N_LAS_Vector * StoVectorPtr,
     N_LAS_Vector * CurrStoVectorPtr,
     N_LAS_Vector * LasStoVectorPtr,
     N_LAS_Vector * stoLeadCurrQCompVectorPtr,
     N_LAS_Vector * QVectorPtr,
     N_LAS_Vector * FVectorPtr,
     N_LAS_Vector * dFdxdVpVectorPtr,
     N_LAS_Vector * dQdxdVpVectorPtr,
     N_LAS_Matrix * dQdxMatrixPtr,
     N_LAS_Matrix * dFdxMatrixPtr
    )
{
  return (
      anaIntPtr_->evalTransientModel(
      t,
      SolVectorPtr,
      CurrSolVectorPtr,
      LasSolVectorPtr,

      StaVectorPtr,
      CurrStaVectorPtr,
      LasStaVectorPtr,
      StaDerivVectorPtr,
      StoVectorPtr,
      CurrStoVectorPtr,
      LasStoVectorPtr,
      stoLeadCurrQCompVectorPtr,
      QVectorPtr,
      FVectorPtr,
      dFdxdVpVectorPtr,
      dQdxdVpVectorPtr,
      dQdxMatrixPtr,
      dFdxMatrixPtr
      )
    );
}

//---------------------------------------------------------------------------
// Function      : Simulator::evalTransientModelState
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//---------------------------------------------------------------------------
bool Simulator::evalTransientModelState
    (
     double t,
     N_LAS_Vector * SolVectorPtr,
     N_LAS_Vector * StaVectorPtr,
     N_LAS_Vector * StoVectorPtr
    )
{
  return (
      anaIntPtr_->evalTransientModelState(
      t,
      SolVectorPtr,
      StaVectorPtr,
      StoVectorPtr
      )
    );
}

//---------------------------------------------------------------------------
// Function      : Simulator::getMapsAndGraphs
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
void Simulator::getMapsAndGraphs
  (
   RCP<N_PDS_ParMap> & x_map,
   RCP<N_PDS_ParMap> & x_map_ognd,
   RCP<N_PDS_ParMap> & s_map,
   RCP<N_PDS_ParMap> & store_map,
   RCP<Epetra_CrsGraph> & dQdx_graph,
   RCP<Epetra_CrsGraph> & dQdx_graph_ognd,
   RCP<Epetra_CrsGraph> & dFdx_graph,
   RCP<Epetra_CrsGraph> & dFdx_graph_ognd
  )
{
  lasBuilderPtr_->getSolutionMaps(x_map,x_map_ognd);
  lasBuilderPtr_->getStateMap(s_map);
  lasBuilderPtr_->getStoreMap(store_map);
  lasBuilderPtr_->getdQdxGraphs(dQdx_graph,dQdx_graph_ognd);
  lasBuilderPtr_->getdFdxGraphs(dFdx_graph,dFdx_graph_ognd);
}


//---------------------------------------------------------------------------
// Function      : Simulator::getVariableNames
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey
// Creation Date : 07/28/09
//---------------------------------------------------------------------------
std::vector<std::string> Simulator::getVariableNames()
{
  return outMgrPtr_->getVariableNames();
}

void
report_handler(
  const char *  message,
  unsigned      report_mask)
{
  // if ( !comm_->isSerial() && !(report_mask & MSG_SYMMETRIC))
  //   os << "P" << comm_->procID() << " - ";

  std::ostringstream oss;
  Xyce::Util::word_wrap(oss, message, 78, " ", "");

  // If symmetric then all processors are getting the same message, only write to p0 and not to ~p0 backlog.  If
  // asymetric then one processor is getting the message, write to per processor stream which writes to per processor
  // log file and to backlog.
  if (report_mask & Xyce::Report::MSG_SYMMETRIC)
    Xyce::lout() << oss.str();
  else
    Xyce::pout() << oss.str();

  // If fatal error also send the message to the standard error file:
  // Also save it for output on proc 0 if running in parallel
  if (report_mask & Xyce::Report::MSG_TERMINATE)
  {
    std::cerr << oss.str() << std::endl;
    Xyce::Report::abort();
  }
}

} //namespace Circuit
} // namespace Xyce
