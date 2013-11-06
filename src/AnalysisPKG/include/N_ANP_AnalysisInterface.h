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
// Filename      : $RCSfile: N_ANP_AnalysisInterface.h,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 01/21/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.36.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_ANALYSIS_INTERFACE_H
#define Xyce_ANALYSIS_INTERFACE_H

// ---------- Standard Includes ----------

#include <list>
#include <vector>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_DEV_fwd.h>
#include <N_UTL_fwd.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Timer.h>
#include <N_NLS_Manager.h>

#include <N_IO_fwd.h>
#include <N_IO_PkgOptionsMgr.h>

// ---------- Forward Declarations ----------
class N_MPDE_Manager;

class N_TIA_TIAParams;
class N_TIA_TimeIntInfo;
class N_TIA_TwoLevelError;

class N_NLS_Manager;

class N_LAS_System;

class N_IO_CmdParse;

class N_LOA_Loader;

class N_IO_RestartMgr;
class N_IO_PkgOptionsMgr;

class N_PDS_Comm;
class N_PDS_Manager;
class N_ANP_AnalysisManager;
class N_TIA_MPDEInterface;
class N_UTL_BreakPoint;

class N_LOA_NonlinearEquationLoader;
class N_TOP_Topology;
class N_LAS_Builder;

// ---------- Enum Definitions ----------
enum Vector_Tag
{
  SOLUTION_VECTOR, // 0
  STATE_VECTOR, // 1
  NUM_VECTOR_TYPES // 2
};

enum ANP_Analysis_Mode
{
  ANP_MODE_INVALID,
  ANP_MODE_DC_OP,
  ANP_MODE_DC_SWEEP,
  ANP_MODE_TRANSIENT,
  ANP_MODE_MPDE,
  ANP_MODE_HB,
  ANP_MODE_AC,
  ANP_MODE_MOR,
  ANP_MODE_NUM_MODES
};

//-----------------------------------------------------------------------------
// Class         :
// Purpose       : This function converts between N_NLS_Manager.h AnalysisMode
//               : and N_ANP_AnalysisInterface.h ANP_Analysis_Mode
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
AnalysisMode anpAnalysisModeToNLS(ANP_Analysis_Mode mode);




//-----------------------------------------------------------------------------
// Class         : N_ANP_AnalysisInterface
// Purpose       :
// Special Notes :
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 11/09/00
//-----------------------------------------------------------------------------
class N_ANP_AnalysisInterface
{
public:

  // Default donstructor
  N_ANP_AnalysisInterface(N_IO_CmdParse & cp);

  // Destructor
  ~N_ANP_AnalysisInterface();


  // Execution functions:
  void resetAll();

  // Execute the control loop for the set analysis type.
  bool run();

  // Method to register the TIA parameters object.
  bool registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp);

  // Method to register the linear system object pointer.
  bool registerLinearSystem(N_LAS_System * lasSysPtr);

  // Method to register the nonlinear solver manager pointer.
  bool registerNLSManager(N_NLS_Manager * nlsMgrPtr);

  // Method to register the nonlinear loader object pointer.
  bool registerLoader(N_LOA_Loader * loaderPtr);

  // Method to register the ?DAE? Loader pointer
  bool registerNonlinearEquationLoader( N_LOA_NonlinearEquationLoader * nonlinearEquationLoaderPtr );

  // Method to register the output manager pointer.
  bool registerOutputMgr(N_IO_OutputMgr * outputPtr);

  // Method to register the restart manager pointer.
  bool registerRestartMgr(N_IO_RestartMgr * restartPtr);

  // Method to register the Device inerfaace pointer
  bool registerDeviceInterface ( N_DEV_DeviceInterface * devInterfacePtr );

  // Method to register the Parallel Services Manager pointer
  bool registerParallelManager( N_PDS_Manager * pdsMgrPtr );

  // Method to register the Topology pointer
  bool registerTopology( N_TOP_Topology * topoMgrPtr );

  // Method to register the Application Builder pointer
  bool registerApplicationBuilder( N_LAS_Builder * appBuilderPtr );

  // Method to register the package options manager
  bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

  // Method to set the transient solution parameters.
  bool setTranAnalysisParams(const N_UTL_OptionBlock & tiaParamBlock);

  // Method to set the DC sweep solution parameters.
  bool setDCAnalysisParams(const N_UTL_OptionBlock & tiaParamBlock);

  // Method to handle OP statements.
  bool setOPAnalysisParams(const N_UTL_OptionBlock & tiaParamBlock);

  // Method to set the STEP solution parameters.
  bool setSTEPAnalysisParams(const N_UTL_OptionBlock & tiaParamBlock);

  // Method to handle SAVE statements
  bool setSaveOptions(const N_UTL_OptionBlock & tiaParamBlock);

  // Method to handle DCOP restart statements
  bool setDCOPRestartParams(const N_UTL_OptionBlock & tiaParamBlock);

  // Method to set a pause time (goes down to AnalysisManager)
  void setPauseTime(const double pauseTime);

  // Method to get a pause time (goes down to AnalysisManager)
  double getPauseTime();

  // Method to check if a simulation is currently paused
  bool isPaused();

  // Method to signal that simulation is resuming
  void resumeSimulation();

  // Method to register the AC Analysis options.
  bool setACAnalysisParams(const N_UTL_OptionBlock & OB);

  // Method to register the MOR Analysis options.
  bool setMORAnalysisParams(const N_UTL_OptionBlock & OB);

  // Method to register the utility options.
  bool setMOROptions(const N_UTL_OptionBlock & OB);

  // Method to register the utility options.
  bool setTranOptions(const N_UTL_OptionBlock & OB);

  // Method to register the MPDE Analysis options.
  bool setMPDEAnalysisParams(const N_UTL_OptionBlock & OB);

  // Method to register the MPDE utility options.
  bool setMPDEOptions(const N_UTL_OptionBlock & OB);

  // Method to register the MPDE utility options.
  bool setTRANMPDEOptions(const N_UTL_OptionBlock & OB);

  // Method to register the HB Analysis options.
  bool setHBAnalysisParams(const N_UTL_OptionBlock & OB);

  // Method to register the HBINT options.
  bool setHBOptions(const N_UTL_OptionBlock & OB);

  // Method to register the LINSOL options.
  bool setLinSol(const N_UTL_OptionBlock & OB);

  // Method to register the HB-LINSOL options.
  bool setHBLinSol(const N_UTL_OptionBlock & OB);

  // Method to register the sensitivity options.
  bool setSensOptions(const N_UTL_OptionBlock & OB);

  // Method to register the sensitivity parameter list options.
  bool setSensParamsOptions(const N_UTL_OptionBlock & OB);

  // Method to register the elapsed time timer
  bool registerElapsedTimer (N_UTL_Timer *);

  // Gets flag for MPDE or HB block analysis
  bool getBlockAnalysisFlag() const;

  // Gets flag for transient analysis
  bool getTransientFlag() const;

  // This function performs the final initializations.  Mainly, it initializes
  // the two data container classes, DataStore and StepErrorControl.  It also
  // registers the neccessary vectors with the LAS system class.
  bool initializeAll();

  // Gets the next time value.
  double getTime() const;

  // Gets the next time value.
  double getCurrentTime() const;

  // Gets the final time value.
  double getFinalTime() const;

  // Gets the initial time value.
  double getInitialTime() const;

  // Updates the divided difference values.
  bool updateDivDiffs();

  // Calls the time int. method to update the corrector derivatives.
  bool updateDerivs();

  // updates state vectors after initial step of starting with a previous OP
  bool completeOPStartStep();

  // updates the State vectors. This function is called from the LOCA interface.
  bool completeHomotopyStep
    ( const vector<string> & paramNames,
      const vector<double> & paramVals,
      N_LAS_Vector * solnVecPtr );

  // Called from the LOCA interface:
  bool failHomotopyStep ();

  bool startTimeStep (const N_TIA_TimeIntInfo & tiInfo);

  // This function equates the 6 temporary vectors with their "next" vector
  // equivalents.  This function is neccessary for the nonlinear solver damping
  // loop.
  bool equateTmpVectors();

  bool updateDerivsBlock(const list < index_pair > & solGIDList,
                         const list < index_pair > & staGIDlist);

  // Gets the size of the restart data (bytes?).
  int restartDataSize( bool pack = true );

  // Dumps the restart data to a file.
  bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack = true );

  // Restores the restart data from a file.
  bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack = true );

  // Gets the solution variable data.
  bool getSolnVarData(const int & gid, vector < double > & varData);

  // Gets the state variable data.
  bool getStateVarData(const int & gid, vector < double > & varData);

  // Gets the store variable data.
  bool getStoreVarData(const int & gid, vector < double > & varData);

  // Sets the solution variable data.
  bool setSolnVarData(const int & gid, const vector < double > & varData);

  // Sets the state variable data.
  bool setStateVarData(const int & gid, const vector < double > & varData);

  // Sets the store variable data.
  bool setStoreVarData(const int & gid, const vector < double > & varData);

  // set begining integration to true
  void setBeginningIntegrationFlag(bool bif);

  // Outputs summary information for time-integration function.
  bool outputSummary();

  // Registers the pointer to the Parallel Distribution Manager.
  bool registerParallelServices(N_PDS_Manager * pds);

  // Accessor so that an outsider can use a different solution vector.
  bool setNextSolVectorPtr (N_LAS_Vector * solVecPtr);

  // Accessor to the N_TIA_TIAParams object.
  N_TIA_TIAParams& getTIAParams();

  // .TRAN
  struct N_ANP_AnalysisInterface_TransAnalysisReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_TransAnalysisReg( N_ANP_AnalysisInterface * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setTranAnalysisParams( options ); }

    N_ANP_AnalysisInterface * Mgr;
  };

  // .options TIMEINT
  struct N_ANP_AnalysisInterface_TranOptionsReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_TranOptionsReg( N_ANP_AnalysisInterface * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setTranOptions( options ); }

    N_ANP_AnalysisInterface * Mgr;
  };

  // .DC
  struct N_ANP_AnalysisInterface_DCAnalysisReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_DCAnalysisReg( N_ANP_AnalysisInterface * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setDCAnalysisParams( options ); }

    N_ANP_AnalysisInterface * Mgr;
  };

  // .options OP_IO
  struct N_ANP_AnalysisInterface_DCOPOptionsReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_DCOPOptionsReg( N_ANP_AnalysisInterface * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setDCOPRestartParams( options ); }

    N_ANP_AnalysisInterface * Mgr;
  };


  // .OP
  struct N_ANP_AnalysisInterface_OPAnalysisReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_OPAnalysisReg( N_ANP_AnalysisInterface * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setOPAnalysisParams( options ); }

    N_ANP_AnalysisInterface * Mgr;
  };

  // .STEP
  struct N_ANP_AnalysisInterface_STEPAnalysisReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_STEPAnalysisReg( N_ANP_AnalysisInterface * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setSTEPAnalysisParams( options ); }

    N_ANP_AnalysisInterface * Mgr;
  };

  // .options SAVE
  struct N_ANP_AnalysisInterface_SaveOptionsReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_SaveOptionsReg( N_ANP_AnalysisInterface * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setSaveOptions( options ); }

    N_ANP_AnalysisInterface * Mgr;
  };

  // .MPDE
  struct N_ANP_AnalysisInterface_MPDE_AnalysisReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_MPDE_AnalysisReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setMPDEAnalysisParams( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .HB
  struct N_ANP_AnalysisInterface_HB_AnalysisReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_HB_AnalysisReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setHBAnalysisParams( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .options HBINT
  struct N_ANP_AnalysisInterface_HB_OptionsReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_HB_OptionsReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setHBOptions( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .options LINSOL
  struct N_ANP_AnalysisInterface_LinSolReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_LinSolReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setLinSol( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .options LINSOL-HB
  struct N_ANP_AnalysisInterface_HB_LinSolReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_HB_LinSolReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setHBLinSol( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .options MPDEINT
  struct N_ANP_AnalysisInterface_MPDE_OptionsReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_MPDE_OptionsReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setMPDEOptions( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .options TIMEINT-MPDE
  struct N_ANP_AnalysisInterface_MPDE_TranMPDEOptionsReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_MPDE_TranMPDEOptionsReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setTRANMPDEOptions( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .AC
  struct N_ANP_AnalysisInterface_AC_AnalysisReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_AC_AnalysisReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setACAnalysisParams( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .MOR
  struct N_ANP_AnalysisInterface_MOR_AnalysisReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_MOR_AnalysisReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setMORAnalysisParams( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // .options MOR_OPTS
  struct N_ANP_AnalysisInterface_MOR_OptionsReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_MOR_OptionsReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setMOROptions( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  
  // .options SENS 
  struct N_ANP_AnalysisInterface_SensOptionsReg : public N_IO_PkgOptionsReg
  {
    N_ANP_AnalysisInterface_SensOptionsReg( N_ANP_AnalysisInterface * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setSensOptions( options ); }

    N_ANP_AnalysisInterface * anpInt_;
  };

  // Function calls relating to residual and jacobian assembly
  // These all come from the NonlinearEquationLoader class.
  bool loadRHS ();
  bool loadJacobian ();
  bool applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result); // for HB Only (matrix free)
  double getResidualTime();
  double getJacobianTime();
  // End of NonlinearEquationLoader function calls.

  bool simulationComplete();

  // Execute the control loop for the
  // set analysis type, for a set number of steps.
  bool runStep
    (const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError);

  void condutanceTest ();
  bool startupSolvers ();
  bool finishSolvers ();

  bool provisionalStep (double maxTimeStep, double &currTimeStep);

  void acceptProvisionalStep ();
  void rejectProvisionalStep ();

  void homotopyStepSuccess
    ( const vector<string> & paramNames,
      const vector<double> & paramVals);

  void homotopyStepFailure ();

  void stepSuccess (int analysis);
  void stepFailure (int analysis);
  bool getInitialQnorm (N_TIA_TwoLevelError & tle);
  bool getBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes);

#ifdef Xyce_Dakota
  // routines to get/set Dakota run flags and actually run a Dakota iteration
  // calles underlying ControlAlgorithm functions
  bool getDakotaRunFlag();
  void setDakotaRunFlag( bool flag );
  int getDakotaIteration();
  void setDakotaIteration( int iterNumber );
#endif

  // Get all the various solver info:
  void getTimeIntInfo(N_TIA_TimeIntInfo & tiInfo);

  // manage verbosity of "Percent Progress" printing
  void silenceProgress();
  void enableProgress();

  N_ANP_AnalysisManager * getAnalysisMgr () const { return &*anaManagerPtr_; }

  void initializeTransientModel();
  bool evalTransientModel(
      double t,
      N_LAS_Vector * SolVectorPtr,
      N_LAS_Vector * CurrSolVectorPtr,
      N_LAS_Vector * LastSolVectorPtr,
      N_LAS_Vector * StateVectorPtr,
      N_LAS_Vector * CurrStateVectorPtr,
      N_LAS_Vector * LastStateVectorPtr,
      N_LAS_Vector * StateDerivVectorPtr,
      N_LAS_Vector * StoreVectorPtr,
      N_LAS_Vector * CurrStoreVectorPtr,
      N_LAS_Vector * LastStoreVectorPtr,
      N_LAS_Vector * stoLeadCurrQCompVectorPtr,
      N_LAS_Vector * QVectorPtr,
      N_LAS_Vector * FVectorPtr,
      N_LAS_Vector * dFdxdVpVectorPtr,
      N_LAS_Vector * dQdxdVpVectorPtr,
      N_LAS_Matrix * dQdxMatrixPtr,
      N_LAS_Matrix * dFdxMatrixPtr
      );
  bool evalTransientModelState(
      double t,
      N_LAS_Vector * SolVectorPtr,
      N_LAS_Vector * StaVectorPtr ,
      N_LAS_Vector * StoVectorPtr
      );

private :

  // Pointer to the TIA AnalysisManager.
  RefCountPtr<N_ANP_AnalysisManager> anaManagerPtr_;

  // commandline object
  N_IO_CmdParse & commandLine_;

  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;

};

#endif //_TIME_MANAGER_H

