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
// Filename      : $RCSfile: N_MPDE_Manager.h,v $
//
// Purpose       : This file defines the manager class for the MPDE package
//
// Special Notes :
//
// Creator       : Robert Hoekstra, 9233, Computational Sciences
//
// Creation Date : 3/11/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.74.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:46 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_MANAGER_H
#define Xyce_MPDE_MANAGER_H

// ---------- Standard Includes ----------

#include <string>
#include <map>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_IO_fwd.h>

#include <N_LAS_BlockVector.h>

#include <N_IO_PkgOptionsMgr.h>

#include <N_MPDE_State.h>
#include <N_TIA_TIAParams.h>
#include <N_MPDE_WarpedPhaseCondition.h>
#include <N_MPDE_DeviceInterface.h>

// ----------   Trilinos Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Forward Declarations ----------

class N_MPDE_Loader;
class N_MPDE_Builder;
class N_MPDE_Discretization;

class N_ANP_AnalysisInterface;
class N_TIA_MPDEInterface;
class N_NLS_Manager;
class N_PDS_Manager;
class N_TOP_Topology;
class N_IO_RestartMgr;
class N_IO_CmdParse;

class N_LOA_Loader;
class N_LOA_NonlinearEquationLoader;
class N_LAS_Builder;
class N_LAS_System;

class N_LAS_PrecondFactory;

class N_UTL_Timer;

// ---------- Enum Definitions ----------

//-----------------------------------------------------------------------------
// Class         : N_MPDE_Manager
// Purpose       : MPDE Manager Class
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
class N_MPDE_Manager
{
 public:

  enum MPDE_IC
  {
    MPDE_IC_DCOP,         // 0  -- use DC sweep to approx IC
    MPDE_IC_SAWTOOTH,     // 1  -- use many small transients to approx IC.
    MPDE_IC_TRANSIENT,    // 2  -- use a transient sim as an approx IC.
    MPDE_IC_TRAN_MAP,     // 3  -- use many transient sims map out an IC.
    MPDE_IC_TWO_PERIOD    // 4  -- use two periods to interpolate IC.
  };


  // Default constructor
  N_MPDE_Manager(N_IO_CmdParse & cp);

  // Destructor
  ~N_MPDE_Manager() {}

  //Runs MPDE analysis
  bool run();

  // MPDE analysis flag
  bool blockAnalysisFlag() const { return blockAnalysisFlag_; }

  //Registrations
  void registerAnalysisInterface( Teuchos::RCP<N_ANP_AnalysisInterface> anaIntPtr );

  void registerTIAMPDEInterface( Teuchos::RCP<N_TIA_MPDEInterface> tiaMPDEIfacePtr );

  void registerNonlinearSolver( Teuchos::RCP<N_NLS_Manager> nlsMgrPtr );

  void registerDeviceInterface ( Teuchos::RCP<N_DEV_DeviceInterface> devInterfacePtr );

  void registerParallelManager( Teuchos::RCP<N_PDS_Manager> pdsMgrPtr );

  void registerTopology( Teuchos::RCP<N_TOP_Topology> topoMgrPtr );

  void registerRestartManager( Teuchos::RCP<N_IO_RestartMgr> resMgrPtr );

  void registerOutputManager( Teuchos::RCP<N_IO_OutputMgr> outMgrPtr );

  void registerElapsedTimer ( Teuchos::RCP<N_UTL_Timer> );

  //registrations of application system builder and loader
  void registerApplicationLoader( Teuchos::RCP<N_LOA_Loader> appLoaderPtr );

  void registerNonlinearEquationLoader( Teuchos::RCP<N_LOA_NonlinearEquationLoader> appLoaderPtr );

  void registerApplicationBuilder( Teuchos::RCP<N_LAS_Builder> appBuilderPtr );

  void registerLinearSystem( Teuchos::RCP<N_LAS_System> lasSysPtr );

  // Method to register the utility options.
  bool setMPDEAnalysisParams(const N_UTL_OptionBlock & OB);

  // Method to register the utility options.
  bool setMPDEOptions(const N_UTL_OptionBlock & OB);

  // Method to register the utility options.
  bool registerTranMPDEOptions(const N_UTL_OptionBlock & OB);

  // "get" function for the fast time scale.
  //double getFastTime ();

  // "get" function for MPDE flag. (true if not IC)
  bool getMPDEFlag ();

  bool getMPDEStartupFlag();
  // get function for MPDE initial condition flag (true if MPDE & IC)
  bool getMPDEIcFlag();

  // "get" function for WaMPDE flag. (true if not IC)
  bool getWaMPDEFlag ();

    const vector<double> & getFastTimePoints () const;

    const vector<double> & getFreqPoints () const;

  // "get" function for GID of phi variable in Warped MPDE case
  int getPhiGID ();

 private :
  void printParams_ ();

 public:
  int debugLevel;

 private:

  //Run Stages
  bool initializeAll_();
  bool initializeOscOut_();
  bool runInitialCondition_();
  bool runDCOP_();
  bool runStartupPeriods_();
  bool runTransientIC_();
  bool filterFastTimePoints_();
  bool setupMPDEProblem_();
  bool runMPDEProblem_();
  bool runTests_();

  void setMPDEOnFlag( bool flagVal );
  void setRunFlag( bool flagVal );
  void setBlockAnalysisFlag_( bool flagVal );

  // diagnostic function
  double checkPeriodicity_();

  N_IO_CmdParse & commandLine_;

  Teuchos::RCP<N_ANP_AnalysisInterface> anaIntPtr_;
  Teuchos::RCP<N_TIA_MPDEInterface> tiaMPDEIfacePtr_;

  N_TIA_TIAParams tiaMPDEParams_;

  N_UTL_OptionBlock saved_tranMPDEOB_;

  Teuchos::RCP<N_NLS_Manager  > nlsMgrPtr_;

  Teuchos::RCP< N_DEV_DeviceInterface > devInterfacePtr_;
  Teuchos::RCP< N_MPDE_DeviceInterface > mpdeDeviceInterfacePtr_;

  Teuchos::RCP<N_PDS_Manager  > pdsMgrPtr_;

  Teuchos::RCP<N_TOP_Topology > topoMgrPtr_;

  Teuchos::RCP<N_IO_RestartMgr> resMgrPtr_;
  Teuchos::RCP<N_IO_OutputMgr > outMgrPtr_;

  Teuchos::RCP<N_LOA_Loader > appLoaderPtr_;
  Teuchos::RCP<N_LOA_NonlinearEquationLoader > nonlinearEquationLoaderPtr_;
  Teuchos::RCP<N_LAS_Builder> appBuilderPtr_;
  Teuchos::RCP<N_LAS_System> lasSysPtr_;
  Teuchos::RCP<N_LAS_System> lasMPDESysPtr_;

  N_MPDE_State mpdeState_;
  Teuchos::RCP<N_MPDE_Loader> mpdeLoaderPtr_;
  Teuchos::RCP<N_MPDE_Builder> mpdeBuilderPtr_;
  Teuchos::RCP<N_MPDE_Discretization> mpdeDiscPtr_;

  Teuchos::RCP<N_UTL_Timer> ElapsedTimerPtr_;

  Teuchos::RCP<N_LAS_Vector> dcOpSolVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpStateVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpQVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpStoreVecPtr_;

  Teuchos::RCP<N_LAS_Vector> endIcSolVecPtr_;
  Teuchos::RCP<N_LAS_Vector> endIcStateVecPtr_;
  Teuchos::RCP<N_LAS_Vector> endIcQVecPtr_;
  Teuchos::RCP<N_LAS_Vector> endIcStoreVecPtr_;

  // Do we have a .MPDE analysis line
  bool blockAnalysisFlag_;

  //Testing Flag
  bool test_;

  //MPDE Problem Size Factor
  int size_;
  bool tranRunForSize_;      // use an initial transient run to calculate size_
  int maxCalcSize_;          // max size to use from transient run.
  bool maxCalcSizeGiven_;

  //MPDE Fast Src driving the fast time scale oscillations
  string fastSrc_;
  bool fastSrcGiven_;
  vector<string> srcVec_;

  // Independent variable for warped MPDE.
  string oscOut_;
  bool oscOutGiven_;

  // Number of steps to take during the start of the full MPDE
  // calculation before turning on truncation errror control
  int nonLteSteps_;
  bool nonLteStepsGiven_;

  //MPDE Fast Time Scale Period
  double period_;
  bool periodGiven_;

  //MPDE number of fast time periods to integrate over and IGNORE before
  // getting initial conditions for MPDE.  Default is zero.
  int startUpPeriods_;
  bool startUpPeriodsGiven_;
  bool startUpPeriodsFinished_;
  bool saveIcData_;

  // MPDE initial condition
  Teuchos::RCP<N_LAS_BlockVector> mpdeICVectorPtr_;

  // MPDE initial state condition
  Teuchos::RCP<N_LAS_BlockVector> mpdeICStateVectorPtr_;

  // MPDE initial Q condition
  Teuchos::RCP<N_LAS_BlockVector> mpdeICQVectorPtr_;

  // MPDE initial store condition
  Teuchos::RCP<N_LAS_BlockVector> mpdeICStoreVectorPtr_;

  //MPDE fast time points
  // 12/5/06 tscoffe:  Note, the period T2 is the last element in fastTimes.
  // This means that the number of fast time points is fastTimes_.size()-1
  vector<double> fastTimes_;

  vector<double> freqPoints_;

  // Fast time discretization
  int fastTimeDisc_;
  int fastTimeDiscOrder_;

  //if we pull data directly from an initial transient run, keep
  // a list of the indices we used so that we can pull out the solution
  // and state data too.
  vector<int> indicesUsed_;
  vector<bool> nonPeriodicFlags;
  // warped MPDE setting in netlist
  bool warpMPDE_;

  // WaMPDE phase condition class
  Teuchos::RCP<N_MPDE_WarpedPhaseCondition> warpMPDEPhasePtr_;

  // WaMPDE OSCOUT GID:
  int warpMPDEOSCOUT_;

  // WaMPDE Phase equation:
  int warpPhase_;
  bool warpPhaseGiven_;

  // WaMPDE Phase equation constant:
  double warpPhaseCoeff_;
  bool warpPhaseCoeffGiven_;

  // frequency domain flag
  bool fftFlag_;

  // Number of fast periods used for initial condition.
  int icPer_;

  // Initial condition strategy.
  int initialCondition_;

  // MPDE mode flag.  if false, run initial condition, if true run MPDE.
  bool mpdeOnFlag_;

  // MPDE IC flag; if true then we're calculating ic's for MPDE
  bool mpdeICFlag_;

  // debug flags:
  bool dcopExitFlag_;
  bool icExitFlag_;
  int  exitSawtoothStep_;

  // An analysis-dependent preconditioner factory.
  Teuchos::RCP<N_LAS_PrecondFactory> precFactory_;

};

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::getMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline bool N_MPDE_Manager::getMPDEFlag()
{
  return mpdeOnFlag_;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::getMPDEStartupFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline bool N_MPDE_Manager::getMPDEStartupFlag()
{
    return startUpPeriodsFinished_;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::getMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437, Elec. Modeling
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline bool N_MPDE_Manager::getMPDEIcFlag()
{
  return mpdeICFlag_;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::getWaMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline bool N_MPDE_Manager::getWaMPDEFlag()
{
  return warpMPDE_;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerAnalysisInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerAnalysisInterface( Teuchos::RCP<N_ANP_AnalysisInterface> anaIntPtr )
{
  anaIntPtr_ = anaIntPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerTIAMPDEInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerTIAMPDEInterface( Teuchos::RCP<N_TIA_MPDEInterface> tiaMPDEIfacePtr )
{
  tiaMPDEIfacePtr_ = tiaMPDEIfacePtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerNonlinearSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerNonlinearSolver( Teuchos::RCP<N_NLS_Manager> nlsMgrPtr )
{
  nlsMgrPtr_ = nlsMgrPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerDeviceInterface( Teuchos::RCP<N_DEV_DeviceInterface> devInterfacePtr )
{
  devInterfacePtr_ = devInterfacePtr;
  mpdeDeviceInterfacePtr_ = Teuchos::rcp( new N_MPDE_DeviceInterface() );
  mpdeDeviceInterfacePtr_->registerDeviceInterface( devInterfacePtr );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerParallelManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerParallelManager( Teuchos::RCP<N_PDS_Manager> pdsMgrPtr )
{
  pdsMgrPtr_ = pdsMgrPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerTopology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerTopology( Teuchos::RCP<N_TOP_Topology> topoMgrPtr )
{
  topoMgrPtr_ = topoMgrPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerRestartManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerRestartManager( Teuchos::RCP<N_IO_RestartMgr> resMgrPtr )
{
  resMgrPtr_ = resMgrPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerOutputManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerOutputManager( Teuchos::RCP<N_IO_OutputMgr> outMgrPtr )
{
  outMgrPtr_ = outMgrPtr;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerApplicationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerApplicationLoader( Teuchos::RCP<N_LOA_Loader> appLoaderPtr )
{
  appLoaderPtr_ = appLoaderPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/06/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerNonlinearEquationLoader( Teuchos::RCP<N_LOA_NonlinearEquationLoader> nonlinearEquationLoaderPtr )
{
  nonlinearEquationLoaderPtr_ = nonlinearEquationLoaderPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerApplicationBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerApplicationBuilder( Teuchos::RCP<N_LAS_Builder> appBuilderPtr )
{
  appBuilderPtr_ = appBuilderPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerLinearSystem( Teuchos::RCP<N_LAS_System> lasSysPtr )
{
  lasSysPtr_ = lasSysPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerElapsedTimer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 08/25/06
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerElapsedTimer (RCP<N_UTL_Timer> et)
{
  ElapsedTimerPtr_ = et;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::getFreqTimePoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 09/02/08
//-----------------------------------------------------------------------------
inline const vector<double> & N_MPDE_Manager::getFreqPoints() const
{
  return freqPoints_;
}



//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::getFastTimePoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 04/06/04
//-----------------------------------------------------------------------------
inline const vector<double> & N_MPDE_Manager::getFastTimePoints() const
{
  return fastTimes_;
}



//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::getPhiGID
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey 1414
// Creation Date : 01/10/07
//-----------------------------------------------------------------------------
inline int N_MPDE_Manager::getPhiGID()
{
  return warpMPDEPhasePtr_->getPhiGID();
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setMPDEOnFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey 1414
// Creation Date : 01/10/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::setMPDEOnFlag( bool flagVal )
{
  mpdeOnFlag_ = flagVal;
  mpdeDeviceInterfacePtr_->setMPDEFlag( flagVal );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey 1414
// Creation Date : 01/10/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::setBlockAnalysisFlag_( bool flagVal )
{
  blockAnalysisFlag_ = flagVal;
  mpdeDeviceInterfacePtr_->setBlockAnalysisFlag( flagVal );
}


#endif //Xyce_MPDE_MANAGER_H
