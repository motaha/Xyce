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
// Revision Number: $Revision: 1.241.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:32 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

// ----------   Xyce Includes   ----------

#include <N_CIR_Xyce.h>

#include <N_DEV_fwd.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DAC.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternalSimulationData.h>

#include <N_ERH_ErrorMgr.h>

#include <IO_NetlistImportTool.h>

#include <N_IO_OutputMgr.h>
#include <N_IO_CmdParse.h>
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

#include <N_UTL_Version.h>
#include <N_UTL_BreakPoint.h>


//--------  Forward Declarations ------------

class N_LAS_QueryUtil;
class N_TOP_TopoLSUtil;

class N_LAS_MultiVector;

//--------  Global Declarations ------------

// Create a global instance of the command line parser that can be
// accessed from anywhere in Xyce.
//N_IO_CmdParse commandLine;

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::N_CIR_Xyce
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
N_CIR_Xyce::N_CIR_Xyce(
#ifdef Xyce_PARALLEL_MPI
                       MPI_Comm * comm
#endif
                      )
: devIntPtr_(0),
#ifdef Xyce_PARALLEL_MPI
  commPtr_(comm),
#endif
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
  Xyce::Device::registerDevices();
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::~N_CIR_Xyce
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
N_CIR_Xyce::~N_CIR_Xyce()
{
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::setNetlistParameters
// Purpose       : This passes a vector of pairs "key" "value" that will
//                 be substituted during the processing of the netlist.  This
//                 more easily allows Dakota to change any netlist parameter
//                 during netlist setup.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 10/9/2008
//---------------------------------------------------------------------------
void N_CIR_Xyce::setNetlistParameters( const vector< pair< string, string > > & externalParams )
{
  externalNetlistParams_ = externalParams;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::setNetlistParameters
// Purpose       : Call through to the output manager to set the suffix to
//                 be used on the output file, as in circuit + suffix + prn
//                 This is useful in Dakota controlled runs to keep each
//                 simulation from overwritting the last one.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and MEMS Modeling
// Creation Date : 10/9/2008
//---------------------------------------------------------------------------
void N_CIR_Xyce::setOutputFileSuffix( const string newSuffix )
{
  if( outMgrPtr_ )
  {
    outMgrPtr_->setOutputFilenameSuffix( newSuffix );
  }
}
//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::setupParMgr_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::setupParMgr_( int iargs, char *cargs[] )
{
  bool bsuccess = true;

  // Setup the Parallel Mgr. with a default load-balance based on the numProc
  // value.
  parMgrPtr_ = new N_PDS_Manager( isSerialFlag_, procZeroFlag_, iargs, cargs
#ifdef Xyce_PARALLEL_MPI
                                  , commPtr_
#endif
                                );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::doAllocations_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::doAllocations_()

{
  bool bsuccess;

  string topotype = "Basic";

  // Allocate device manager:
  devIntPtr_  = N_DEV_DeviceInterface::factory(commandLine);

  // Allocate Topology:
  topMgrPtr_ = Xyce::Topology::Manager::instance();
  Xyce::Topology::Manager & topMgr = (*topMgrPtr_);

  topPtr_ = topMgrPtr_->createTopology(commandLine);

  // Allocate distribution mgr:
  netlistImportToolPtr_ = NetlistImportTool::factory(commandLine, topMgr);

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

  pkgOptionsMgrPtr_ = rcp( new N_IO_PkgOptionsMgr() );

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

  bsuccess = bsuccess && ( !(Teuchos::is_null(pkgOptionsMgrPtr_)) );

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::doDeAllocations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::doDeAllocations_()

{
  // de-allocate the device manager:

  if (devIntPtr_       != 0) delete devIntPtr_;


  if (nlsMgrPtr_       != 0) delete nlsMgrPtr_;
  if (anaIntPtr_       != 0) delete anaIntPtr_;

  if (loaderMgrPtr_    != 0) delete loaderMgrPtr_;
  if (cktLoaderPtr_    != 0) delete cktLoaderPtr_;
  if (nonlinearEquationLoaderPtr_    != 0) delete nonlinearEquationLoaderPtr_;
  if (outMgrPtr_       != 0) delete outMgrPtr_;
  if (lasSysPtr_       != 0) delete lasSysPtr_;
  if (lasBuilderPtr_   != 0) delete lasBuilderPtr_;
  if (XyceTimerPtr_    != 0) delete XyceTimerPtr_;
  if (ElapsedTimerPtr_ != 0) delete ElapsedTimerPtr_;
  if (parMgrPtr_       != 0) delete parMgrPtr_;

  if (N_ERH_ErrorMgr::output != 0) delete N_ERH_ErrorMgr::output;

  if (topMgrPtr_       != 0) delete topMgrPtr_;

  if (resMgrPtr_       != 0) delete resMgrPtr_;

  pkgOptionsMgrPtr_ = Teuchos::null;

  return true;

}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::doRegistrations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::doRegistrations_()

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


#ifdef Xyce_PARALLEL_MPI
  // Distribution manager registrations:
  bs1 = netlistImportToolPtr_->registerParallelServices(parMgrPtr_->getPDSComm());
  bsuccess = bsuccess && bs1;
#endif

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
// Function      : N_CIR_Xyce::setUpTopology_
// Purpose       : This function handles a lot of the initial setup.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/00
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::setUpTopology_()
{
  string netListFile;

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();
#endif
  if( procZeroFlag_ )
  {
    if (commandLine.getArgumentValue("netlist") != "")
      netListFile = commandLine.getArgumentValue("netlist");
    else
      netListFile = "xyce.in";

    FILE * testFile;
    if ( (testFile=fopen(netListFile.c_str(), "r")) == 0)
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, netListFile +
                           " file not found.");
    int ierr = fclose( testFile );
  }

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);
#endif

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                         " ***** Reading and parsing netlist...\n");
  netlistImportToolPtr_->constructCircuitFromNetlist(netListFile, externalNetlistParams_);

  delete netlistImportToolPtr_;

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                         " ***** Setting up topology...\n");
#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();
#endif

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
  if (commandLine.argExists(string("-remeasure")) )
  {
    outMgrPtr_->remeasure();
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                         " ***** Remeasure analysis complete\n");
    return false;
  }
  topPtr_->instantiateDevices();

  outMgrPtr_->delayedPrintLineDiagnostics();

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);
#endif

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();
#endif

#ifdef Xyce_PARALLEL_MPI
  devIntPtr_->setGlobalFlags();
#endif

  // Setup of indices including global reordering.
  topPtr_->setupGlobalIndices();

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);
#endif

#ifdef Xyce_DEBUG_DEVICE
  devIntPtr_->printOutLists();
#endif

  return true;

}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::setUpMatrixStructure_
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
bool N_CIR_Xyce::setUpMatrixStructure_()
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
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
		       " ***** Number of Unknowns = ",lasSize,"\n");

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::doInitializations_
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
bool N_CIR_Xyce::doInitializations_()
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
// Function      : N_CIR_Xyce::runSolvers_
// Purpose       : This function runs the solvers.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::runSolvers_()
{
  return anaIntPtr_->run();
}


//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::run(int iargs_tmp, char *cargs_tmp[])
{
  bool bsuccess = true;
  bool bs1 = true;

  bs1 = initialize(iargs_tmp, cargs_tmp);
  bsuccess = bsuccess && bs1;

  if (bsuccess)
  {
    bs1 = runSimulation();
    bsuccess = bsuccess && bs1;

    bs1 = finalize();
    bsuccess = bsuccess && bs1;
  }
  else
  {
    reportTotalElapsedTime ();
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "Xyce Initialization Phase failed.");
    Xyce_exit(0);
  }

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::runSimulation
// Purpose       : Main simulation driver.
// Special Notes : Not private as this is also called from N_DAK_DakotaInterface
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::runSimulation()
{
  return runSolvers_();
}

// ---------------------------------------------------------------------------
// API METHODS NEEDED FOR MIXED-SIGNAL and other external applications
//
//-----------------------------------------------------------------------------
// Function      : N_CIR_Xyce::initialize
// Purpose       : capture all "initialization-type" activities in one
//                 method
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 02/17/2004
//-----------------------------------------------------------------------------
bool N_CIR_Xyce::initialize( int iargs_tmp, char *cargs_tmp[] )
{
  iargs = iargs_tmp;
  cargs = cargs_tmp;

  // Setup the Parallel Manager first to give us a N_PDS_Comm object so the
  // reporting package will work.
  bool bsuccess = setupParMgr_( iargs, cargs );
  bool b2;

  // Start the solvers timer.
  ElapsedTimerPtr_ = new N_UTL_Timer( *(parMgrPtr_->getPDSComm()) );

  // Banner
  static const string bannerHdg("*****");
  static const string bannerMsg("***** Welcome to the Xyce(TM) "
                                    "Parallel Electronic Simulator\n");

  const string versionMsg = "***** This is version " +
   N_UTL_Version::getFullVersionString() + "\n";

  static const string executeMsg("***** Executing ");
  string netlistMsg = "netlist " ;
  static const string crMsg = ("\n");

#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_DEBUG_ALL_PROCS_SAME_WD
    N_PDS_Comm * commPtr  = parMgrPtr_->getPDSComm();
    size_t dirsize = 512;
    char directory[dirsize];  for (int idir=0;idir<dirsize;++idir)  directory[idir] = 0;
#endif
#endif

#ifdef Xyce_DEBUG_CIRCUIT
  string msg;
#endif

  // register parallel mgr to allow parallel support
  commandLine.registerParallelMgr( parMgrPtr_ );

  // read in command line arguments
  commandLine.parseCommandLine( iargs, cargs );


  if( procZeroFlag_ )
  {
    // This utilizes the "-l <file>" command line option to output to a file
    // instead of stdout.  Currently only works for serial (proc. 0) output.
    // Parallel output will be directed to stdout - SAH, 27 Aug. 2002.

    // Set the output stream of the "-l" flag exists
    if (commandLine.argExists(string("-l")))
    {
      string outputFile;
      outputFile = commandLine.getArgumentValue("-l");

      // Allocate the output stream
      N_ERH_ErrorMgr::output = new ofstream();

      N_ERH_ErrorMgr::output->open(outputFile.c_str());

      if (N_ERH_ErrorMgr::output->fail())
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, "Unable to open output file: ");

      if (!isSerialFlag_)
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, crMsg +
                               "-l <file> syntax not supported in parallel.  "
                               "Parallel output will be directed to stdout "
                               "but serial output directed to processor 0 "
                               "will be sent to the requested file");
    }

    if (commandLine.getArgumentValue("netlist") != "")
      netlistMsg += commandLine.getArgumentValue("netlist") + "...\n";
    else
      netlistMsg = netlistMsg + "MISSING!  Usage: " + string(cargs[0]) +
                   " netlist \n";

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, crMsg + bannerHdg +
                         crMsg + bannerMsg + bannerHdg + crMsg + versionMsg +
                         crMsg + crMsg + executeMsg + netlistMsg + crMsg +
                         crMsg);

    // Don't bother doing anything else if we weren't given a netlist!
    if (commandLine.getArgumentValue("netlist") == "") return (STATUS_FAILURE);

    // check if a parameter file was specified for param substitution during parsing
    if (commandLine.argExists(string("-prf")))
    {
      string parameterFile = commandLine.getArgumentValue("-prf");
      readExternalParamsFromFile( parameterFile, iterationSuffix_, externalNetlistParams_ );
#ifdef Xyce_DEBUG
      // debug print out externNetlistParams_
      vector<pair<string,string> >::iterator currentPair = externalNetlistParams_.begin();
      vector<pair<string,string> >::iterator endPair = externalNetlistParams_.end();
      while( currentPair != endPair )
      {
        std::cout << "\"" << currentPair->first << "\" = \"" << currentPair->second << "\"" << endl;
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

#ifdef Xyce_DEBUG_CIRCUIT
  if (b2 == true)
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg +
                           "Allocation was successful.");
  else
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg +
                           "Allocation was NOT successful.");
#endif

  // Register the external package pointers:
  b2 = doRegistrations_();
  bsuccess = bsuccess && b2;

#ifdef Xyce_DEBUG_CIRCUIT
  if (b2 == true)
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg +
                           "Registration was successful.");
  else
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg +
                           "Registration was NOT successful.");
#endif

  b2 = setUpTopology_();
  if( !b2 )
  {
    return false;
  }
  bsuccess = bsuccess && b2;

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                         " ***** Device Count Summary ...\n");
  outMgrPtr_->printDeviceCounts();

  if( commandLine.argExists( "-norun" ) )
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                         " ***** Syntax and topology analysis complete\n");
    return false;
  }
  else
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           " ***** Setting up matrix structure...\n");
    b2 = setUpMatrixStructure_();
    bsuccess = bsuccess && b2;

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           " ***** Initializing...\n");
    b2 = doInitializations_();
    bsuccess = bsuccess && b2;

#ifdef Xyce_TEST_SOLN_VAR_MAP
    topPtr_->outputNameFile();
#endif

    // if we loaded a parameter file from the command line, then the output manager
    // should scan the results as well as such files can specify response functions
    // than need to be reported by the output manager.
    if (commandLine.argExists(string("-prf")))
    {
      outMgrPtr_->setExternalNetlistParams( externalNetlistParams_ );
      if(iterationSuffix_.length() > 0)
      {
        outMgrPtr_->setOutputFilenameSuffix( iterationSuffix_ );
      }
    }
    if (commandLine.argExists(string("-rsf")))
    {
      string responseFile = commandLine.getArgumentValue("-rsf");
      outMgrPtr_->setResponseFilename( responseFile );
    }
    // Start the solvers timer.
    XyceTimerPtr_ = new N_UTL_Timer( *(parMgrPtr_->getPDSComm()) );
  }

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : getDACDeviceNames
// Purpose        : Gets the (stripped) names of the DAC devices
//                  in the circuit.
// Special Notes  :
// Scope          :
// Creator        : Lisa Maynes
// Creation Date  : 06/13/2003
//----------------------------------------------------------------------------
bool N_CIR_Xyce::getDACDeviceNames(vector< string >& dacNames)
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
bool N_CIR_Xyce::getADCMap(map<string, map<string,double> >&ADCMap)
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
bool N_CIR_Xyce::updateTimeVoltagePairs(
   map< string, vector< pair<double,double> >* > const & timeVoltageUpdateMap)
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
bool N_CIR_Xyce::getTimeVoltagePairs(
   map< string, vector< pair<double,double> > > & timeVoltageUpdateMap)
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
bool N_CIR_Xyce::setADCWidths(
   map< string, int > const & ADCWidthMap)
{
  bool bsuccess = true;

  bsuccess=devIntPtr_->setADCWidths(ADCWidthMap);

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::simulateUntil
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
bool N_CIR_Xyce::simulateUntil(double requestedUntilTime,
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
  std::cout << "N_CIR_Xyce::simulateUntil: ";
  std::cout << "finalTime = " << finalTime
            << ", currentTimeBeforeSim = " << currentTimeBeforeSim << std::endl;
#endif

  if (currentTimeBeforeSim >= finalTime)
  {
    // We have already simulated as far as the netlist requested.
    bsuccess = true;
    completedUntilTime = currentTimeBeforeSim;
#ifdef Xyce_DEBUG_CIRCUIT
    std::cout << "Case1: completedUntilTime = " << completedUntilTime;
#endif
  }
  else
  {
    // We are not already past the end of the netlist time
    anaIntPtr_->setPauseTime(Xycemin(requestedUntilTime,finalTime));
#ifdef Xyce_DEBUG_CIRCUIT
    std::cout << "N_CIR_Xyce::simulateUntil currentTimeBeforeSim = " << currentTimeBeforeSim << "  initialTime = " << initialTime << std::endl;
#endif
    if (currentTimeBeforeSim > initialTime)
    {
      anaIntPtr_->resumeSimulation();
    }

#ifdef Xyce_DEBUG_CIRCUIT
    std::cout << "N_CIR_Xyce::simulateUntil: Case2: requestedUntilTime = " << requestedUntilTime
              << ", pauseTime = " << anaIntPtr_->getPauseTime() << std::endl;
#endif

    bsuccess = runSolvers_();
    completedUntilTime = anaIntPtr_->getTime();
#ifdef Xyce_DEBUG_CIRCUIT
    std::cout << "N_CIR_Xyce::simulateUntil: Case2: completedUntilTime = " << completedUntilTime << std::endl;
#endif
  }

#ifdef Xyce_DEBUG_CIRCUIT
  std::cout << std::endl;
#endif
  //  return true;
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::finalize
// Purpose       : To clean up after driving Xyce with the SIMBUS
//                 simulation backplane. This includes the following:
//                    Free any dynamically allocated memory...
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
//               : Lisa Maynes, CoMeT Solutions, Inc.
// Creation Date : 06/03/2003
//---------------------------------------------------------------------------
bool N_CIR_Xyce::finalize()
{
  bool bsuccess = true;

  static const string crMsg = ("\n");
  static const string bannerHdg("*****");
  static const string timeMsgPre = ("\n\tTotal Simulation Solvers Run "
                                    "Time:\t");
  static const string totTimeMsgPre = ("\n***** Total Elapsed Run Time: ");
  static const string timeMsgPost = (" seconds");
  string msg;
  msg = ("\n***** Solution Summary *****");
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);


  anaIntPtr_->outputSummary();
  outMgrPtr_->outputMacroResults();

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, timeMsgPre,
                         XyceTimerPtr_->elapsedTime(), timeMsgPost);

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, totTimeMsgPre,
                         ElapsedTimerPtr_->elapsedTime(), timeMsgPost);

  // Closing banner
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, crMsg + bannerHdg +
                         crMsg + "***** End of Xyce(TM) Simulation " +
                         crMsg + bannerHdg);

  // Close the output stream:
  if (N_ERH_ErrorMgr::output != 0 && N_ERH_ErrorMgr::output->is_open())
    N_ERH_ErrorMgr::output->close();

  bsuccess = doDeAllocations_();

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::reportTotalElapsedTime ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/01/2009
//---------------------------------------------------------------------------
void N_CIR_Xyce::reportTotalElapsedTime ()
{
  static const string totTimeMsgPre = ("\n***** Total Elapsed Run Time: ");
  static const string timeMsgPost = (" seconds");
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, totTimeMsgPre,
                         ElapsedTimerPtr_->elapsedTime(), timeMsgPost);
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::simulationComplete
// Purpose       : Simply report whether we've reached the end of the
//                 simulation
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//---------------------------------------------------------------------------
bool N_CIR_Xyce::simulationComplete()
{
  return anaIntPtr_->simulationComplete();
}

//
// new mixed-signal functions:
// These are provisional!

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::provisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/2009
//---------------------------------------------------------------------------
bool N_CIR_Xyce::provisionalStep
  (double maxTimeStep,
   double &timeStep,
   map< string, vector< pair<double,double> > > & timeVoltageUpdateMap)
{
  bool bsuccess=true;

  bool b1 = anaIntPtr_->provisionalStep(maxTimeStep, timeStep);
  bsuccess = bsuccess && b1;

  b1=getTimeVoltagePairs(timeVoltageUpdateMap);

  bsuccess = bsuccess && b1;

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::getFinalTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/25/2009
//---------------------------------------------------------------------------
double N_CIR_Xyce::getFinalTime()
{
  double ft=0.0;
  if (anaIntPtr_!=0)
  {
    ft = anaIntPtr_->getFinalTime();
  }
  return ft;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::getTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//---------------------------------------------------------------------------
double N_CIR_Xyce::getTime()
{
  double t=0.0;
  if (anaIntPtr_!=0)
  {
    t = anaIntPtr_->getTime();
  }
  return t;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::acceptProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//---------------------------------------------------------------------------
void N_CIR_Xyce::acceptProvisionalStep()
{
  anaIntPtr_->acceptProvisionalStep ();
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::rejectProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/2009
//---------------------------------------------------------------------------
void N_CIR_Xyce::rejectProvisionalStep()
{
  anaIntPtr_->rejectProvisionalStep();
}

// ---------------------------------------------------------------------------
// API METHODS NEEDED FOR Two-level Functions:

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::simulateStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/2006
//---------------------------------------------------------------------------
bool N_CIR_Xyce::simulateStep
      ( const N_DEV_SolverState & solState,
        const map<string,double> & inputMap,
        vector<double> & outputVector,
        vector< vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError)
{
  bool bsuccess = false;

#ifdef Xyce_DEBUG_CIRCUIT
  cout << "\nN_CIR_Xyce::simulateStep: " << endl;
#endif

  // Apply the input voltages to their appropriate sources:
  map<string,double>::const_iterator iterM = inputMap.begin();
  map<string,double>::const_iterator  endM = inputMap.end  ();
  int i=0;
  for (; iterM != endM;++iterM, ++i)
  {
    bool found = true;
    string name = iterM->first;
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
// Function      : N_CIR_Xyce::simulateStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 07/16/2009
//---------------------------------------------------------------------------
bool N_CIR_Xyce::simulateStep
      (const N_DEV_ExternalSimulationData & ext_data,
       const map<string,double> & inputMap,
       vector<double> & outputVector,
       vector< vector<double> > & jacobian,
       N_TIA_TwoLevelError & tlError)
{
  bool bsuccess = false;

#ifdef Xyce_DEBUG_CIRCUIT
  cout << "\nN_CIR_Xyce::simulateStep: " << endl;
#endif

  // Create N_DEX_SolverState object
  N_DEV_SolverState state;

  if (ext_data.is_transient)
  {
	  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO, "Xychron (mixed-model): Transient Xyce Step" );

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
// Function      : N_CIR_Xyce::startupSolvers
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/2006
//---------------------------------------------------------------------------
bool N_CIR_Xyce::startupSolvers ()
{
  bool bsuccess = true;
  bsuccess = anaIntPtr_->startupSolvers ();
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::finishSolvers
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/2006
//---------------------------------------------------------------------------
bool N_CIR_Xyce::finishSolvers ()
{
  bool bsuccess = true;
  bsuccess = anaIntPtr_->finishSolvers ();
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::homotopyStepSuccess
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//---------------------------------------------------------------------------
void N_CIR_Xyce::homotopyStepSuccess
    (const vector<string> & paramNames,
     const vector<double> & paramVals)
{
  anaIntPtr_->homotopyStepSuccess ( paramNames, paramVals);
  return;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::homotopyStepFailure
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/2006
//---------------------------------------------------------------------------
void N_CIR_Xyce::homotopyStepFailure ()
{
  anaIntPtr_->homotopyStepFailure ();
  return;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::stepSuccess
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
void N_CIR_Xyce::stepSuccess (int analysis)
{
  anaIntPtr_->stepSuccess (analysis);
  return;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::stepFailure
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
void N_CIR_Xyce::stepFailure (int analysis)
{
  anaIntPtr_->stepFailure (analysis);
  return;
}


//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::getInitialQnorm
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
bool N_CIR_Xyce::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  bool bsuccess = true;

  if (initializeAllFlag_)
  {
    bsuccess = anaIntPtr_->getInitialQnorm (tle);
  }

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::getBreakPoints
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
bool N_CIR_Xyce::getBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes)
{
  bool bsuccess = true;

  if (initializeAllFlag_)
  {
    bsuccess = anaIntPtr_->getBreakPoints (breakPointTimes);
  }

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::updateStateArrays
// Purpose       :
// Special Notes : Used for 2-level Newton solves, with LOCA.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/18/2006
//---------------------------------------------------------------------------
bool N_CIR_Xyce::updateStateArrays ()
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
// Function      : N_CIR_Xyce::setInternalParam
// Purpose       :
// Special Notes : Used for 2-level Newton solves, with LOCA.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/18/2006
//---------------------------------------------------------------------------
bool N_CIR_Xyce::setInternalParam (string & name, double val)
{
  cktLoaderPtr_->setParam (name, val);
  return true;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::startTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//---------------------------------------------------------------------------
bool N_CIR_Xyce::startTimeStep (const N_TIA_TimeIntInfo & tiInfo)
{
  bool bsuccess = anaIntPtr_->startTimeStep(tiInfo);
  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::startTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 08/28/2009
//---------------------------------------------------------------------------
bool N_CIR_Xyce::startTimeStep (const N_DEV_ExternalSimulationData & ext_data)
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
// Function      : N_CIR_Xyce::endTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Russell Hooper, SNL
// Creation Date : 10/16/2012
//---------------------------------------------------------------------------
bool N_CIR_Xyce::endTimeStep (N_DEV_ExternalSimulationData & ext_data)
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
// Function      : N_CIR_Xyce::endTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Russell Hooper, SNL
// Creation Date : 10/16/2012
//---------------------------------------------------------------------------
void N_CIR_Xyce::enable_lte_analysis()
{
  N_TIA_TIAParams & tiaParams = anaIntPtr_->getTIAParams();
  tiaParams.errorAnalysisOption = 0; // use local truncation error estimates
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::readExternalParamsFromFile
// Purpose       :
// Special Notes : Used to read parameters "tag" = "value" from an external
//                 file.  Any "tag"s found while reading the netlist during
//                 parsing will be replaced by "value".
// Scope         : public
// Creator       : Richard Schiek, 1437 Electrical Systems Modeling
// Creation Date : 07/24/2012
//---------------------------------------------------------------------------
void N_CIR_Xyce::readExternalParamsFromFile( string filename, string & iterationSuffix,
  vector< pair< string, string > > & paramList )
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

  const string allowedChars("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_.$");
  const string whiteSpace(" \t=\n\r");    // note we will treat "=" as whitespace
  const string commentChars("*#;");  // if we find any of these then the rest of the line is a comment

  // attempt to open the params file
  ifstream paramFile(filename.c_str(), ifstream::in);
  if( paramFile )
  {
    // loop over the file trying to gather names / values and response functions required.
    // general format:
    // white space tag/value white space or = tag/value
    //
    string aLine;
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
        string::size_type commentLoc = aLine.find_first_of( commentChars, 0 );
        if( commentLoc != string::npos )
        {
          // chop off comment to end of line.
          aLine.erase( commentLoc, aLine.length()-commentLoc );
        }
        //check overall line length again.  This could have been just a comment line
        if( aLine.length() > 2 )
        {
          string::size_type word1Start = aLine.find_first_of( allowedChars, 0 );
          // check that we found a valid word1Start otherwise stop trying to parse this line.
          if( word1Start != string::npos )
          {
            string::size_type word1End   = aLine.find_first_of( whiteSpace, word1Start );
            // check that we found a valid word1End otherwise stop trying to parse this line.
            if( word1End != string::npos )
            {
              string::size_type word2Start = aLine.find_first_of( allowedChars, word1End );
              // check that we found a valid word2Start otherwise stop trying to parse this line.
              if( word2Start != string::npos )
              {
                string::size_type word2End   = aLine.find_first_of( whiteSpace, word2Start );
                // check that we found a valid word2End
                if( word2End == string::npos )
                  word2End = aLine.length();
                // if we get here then we have valid start,end indicies for word1 and word2

                string word1=aLine.substr( word1Start, (word1End - word1Start) );
                string word2=aLine.substr( word2Start, (word2End - word2Start) );

                // sort out tag/value ordering.
                // if word1=number assume format is value = tag
                // otherwise assume format is tag = value
                stringstream converter;
                converter << word1;
                double testvalue;
                converter >> testvalue;
                if( converter.fail() )
                {
                  // couldn't convert tag1 to a double value so assume format is word1=tag word2=value
                  paramList.push_back( pair<string,string>(word1,word2) );
                }
                else
                {
                  // tag1 was successfully converted to a value so assume format is word1=value word2=tag
                  paramList.push_back( pair<string,string>(word2,word1) );
                }

              }  // if ( word2Start != string::npos )
            }    // if( word1End != string::npos )
          }      // if( word1Start != string::npos )
        }        // if( aLine.length() > 2 ) -- check after removing comments
      }          // if( aLine.lenght() > 2 ) -- outer check
      // try to get another line
      getline(paramFile, aLine);
    }

  }
  else
  {
    // emit warning that param file could not be found
    string message = "Could not open parameter file: " + filename + ". Attempting to continue.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING, message );
  }

  // for debug purposes.  output the params as read
  std::cout << "Parameters read from \"" << filename << "\"" << std::endl;
  vector< pair< string, string > >::iterator listitr = paramList.begin();
  vector< pair< string, string > >::iterator enditr = paramList.end();
  while( listitr != enditr )
  {
    std::cout << "  " << listitr->first << " , " << listitr->second << std::endl;
    listitr++;
  }



}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::registerResponseVars
// Purpose       :
// Special Notes : Used when Dakota or an external program calls Xyce to tell
//                 Xyce what variables and expressions should be held as
//                 response funcctions to pass back to Dakota.  This
//                 is just pass through to the output Manager
// Scope         : public
// Creator       : Richard Schiek, 1437 Electrical and MEMS modeling
// Creation Date : 10/22/2008
//---------------------------------------------------------------------------
bool N_CIR_Xyce::registerResponseVars (string objString, RCP<vector< double > > varVectorPtr )
{
  if( outMgrPtr_ == 0)
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0, "N_CIR_Xyce::registerResponseVars outMgrPtr_ is null");
  }
  bool result = outMgrPtr_->registerResponseVars( objString, varVectorPtr );
  return result;

}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::registerResponseVars
// Purpose       :
// Special Notes : Used when Dakota or an external program calls Xyce to tell
//                 Xyce what variables and expressions should be held as
//                 response funcctions to pass back at end of simulation.
//                 This is just pass through to the output Manager
// Scope         : public
// Creator       : Richard Schiek, 1437 Electrical and MEMS modeling
// Creation Date : 10/22/2008
//---------------------------------------------------------------------------
void N_CIR_Xyce::finalizeResponseVars()
{
  if( outMgrPtr_ == 0)
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0, "N_CIR_Xyce::finalizeResponseVars outMgrPtr_ is null");
  }
  outMgrPtr_->finalizeResponseVars();
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::initializeTransientModel
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
void N_CIR_Xyce::initializeTransientModel()
{
  anaIntPtr_->initializeTransientModel();
}

//---------------------------------------------------------------------------
// Function      : N_CIR_Xyce::evalTransientModel
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
bool N_CIR_Xyce::evalTransientModel
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
// Function      : N_CIR_Xyce::evalTransientModelState
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//---------------------------------------------------------------------------
bool N_CIR_Xyce::evalTransientModelState
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
// Function      : N_CIR_Xyce::getMapsAndGraphs
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//---------------------------------------------------------------------------
void N_CIR_Xyce::getMapsAndGraphs
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
// Function      : N_CIR_Xyce::getVariableNames
// Purpose       :
// Special Notes : Used for ModelEvaluator interface
// Scope         : public
// Creator       : Coffey
// Creation Date : 07/28/09
//---------------------------------------------------------------------------
std::vector<std::string> N_CIR_Xyce::getVariableNames()
{
  return outMgrPtr_->getVariableNames();
}

