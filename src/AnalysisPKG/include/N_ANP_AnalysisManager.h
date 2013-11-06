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
// Filename      : $RCSfile: N_ANP_AnalysisManager.h,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.49.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIME_INTEG_ALG_H
#define Xyce_N_TIME_INTEG_ALG_H

// ---------- Standard Includes ----------

#include <list>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>
#include <N_IO_fwd.h>
#include <N_UTL_OptionBlock.h>

#include <N_TIA_TIAParams.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TimeIntegrationEnum.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TimeIntInfo.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>

#include <N_LOA_NonlinearEquationLoader.h>
#include <N_DEV_DeviceInterface.h>
#include <N_TOP_Topology.h>
#include <N_LAS_Builder.h>
#include <N_NLS_Manager.h>

// ---------- Forward Declarations ----------

class N_TIA_Assembler;
class N_TIA_TimeIntegrationMethod;
class N_TIA_TwoLevelError;
class N_TIA_MPDEInterface;

class N_NLS_Manager;

class N_LAS_System;

class N_LOA_Loader;

class N_UTL_Timer;

class N_IO_RestartMgr;
class N_IO_CmdParse;

class N_PDS_Comm;
class N_PDS_Manager;

class N_MPDE_Manager;

//-----------------------------------------------------------------------------
// Class         : N_ANP_AnalysisManager
//
// Purpose       : This class manages, allocates, and sets up the
//                 various analysis types, such as DC, Tran, HB, etc.
//
// Special Notes : Some of this class was once in the N_TIA_ControlAlgorithm
//                 class, which was set up back when Xyce only did transient
//                 simulations.  As we added analysis types it became necessary
//                 to refactor the code so that each analysis type had its own
//                 set of classes.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 6/01/00 (N_TIA_ControlAlgorithm, now deprecated)
// Creation Date : 1/24/08 (for this version of the class. Date is approximate)
//-----------------------------------------------------------------------------
class N_ANP_AnalysisManager
{
  public:

    // Default constructor.
    N_ANP_AnalysisManager(N_IO_CmdParse & cp, N_ANP_AnalysisInterface * anaIntPtr);

    // Destructor
    ~N_ANP_AnalysisManager();

    // Execution functions:
    void resetAll();

    // Execute the control loop for the set analysis type.
    bool run();

    // This function performs the final initializations.  Mainly, it initializes
    // the two data container classes, DataStore and StepErrorControl.  It also
    // registers the neccessary vectors with the LAS system class.
    bool initializeAll();

    // Returns the current partial time derivative
    double partialTimeDerivative();

    // Gets the next time-step value.
    double getTime() const;

    // Gets the current time-step value.
    double getCurrentTime() const;

    // Gets the final time-step value.
    double getFinalTime() const;

    // Gets the initial time-step value.
    double getInitialTime() const;

    // Gets the starting time-step value.
    double getStartingTimeStep();

    // Updates the divided difference values.
    bool updateDivDiffs();

    // Calls the time int. method to update the corrector derivatives.
    bool updateDerivs();

    // updates state vectors following initial solve with previous operating point constraints
    bool completeOPStartStep();

    // updates the State vectors. This function is called from the LOCA interface.
    bool completeHomotopyStep
      ( const vector<string> & paramNames,
        const vector<double> & paramVals,
        N_LAS_Vector * solnVecPtr );

    bool failHomotopyStep ();

    // This function equates the 6 temporary vectors with their "next" vector
    // equivalents.  This function is neccessary for the nonlinear solver damping
    // loop.
    bool equateTmpVectors();

    // Compute an estimate of the error in the integration step.
    bool updateDerivsBlock(const list < index_pair > & solGIDList,
                           const list < index_pair > & staGIDList);

    // Prints out time loop information.
    bool printLoopInfo(int start, int finish);

    // Gets the time-integration method.
    int getTimeIntMode();

    // Get the steady-state flag (true if time int mode is none)
    bool getSteadyStateFlag ();

    // Get the dcop flag for transient.(true if doing DCOP for transient initialization)
    bool getTranOPFlag ();

    // Get the dcop flag for AC.(true if doing DCOP for AC initialization)
    bool getACOPFlag ();

    bool getDCSweepFlag ();
    
    bool getDotOpFlag() {return dotOpSpecified_;};
    bool getStepFlag() {return stepLoopFlag_;};

    bool getSweepSourceResetFlag () {return sweepSourceResetFlag_;};
    void setSweepSourceResetFlag (bool ssrf) { sweepSourceResetFlag_=ssrf;} ;

    bool getTransientFlag () const;

    // Is the doubleDCOP algorithm enabled?
    bool getDoubleDCOPEnabled ();

    // Get block analysis information for HB
    void setBlockAnalysisFlag( bool flagVal ) { blockAnalysisFlag_ = flagVal; }
    bool getBlockAnalysisFlag() const;  
    
    void setHBFlag( bool flagVal ) { hbFlag_ = flagVal; }
    bool getHBFlag () {return hbFlag_; }

    // gets the index of the DCOP step.  
    // 0 = nonlinear poisson, 1=full TCAD 
    int getDoubleDCOPStep();

    // Gets/sets the step number.
    int getStepNumber();
    int getTranStepNumber();

    void setStepNumber(int step);
    void setTranStepNumber(int step);

    // This is true only at the beginning of integration, not at a breakpoint.
    bool getInitTranFlag();

    // Gets the step number.
    //int getStepLoopIter () { return stepLoopIter_; }

    // Returns the time integration order
    int getOrder ();

    int getNumberOfSteps ();
    int getUsedOrder ();
    int getNscsco ();

    // Returns the "current" time step size.
    double getCurrentStepSize();

    // Returns the "last" time step size.
    double getLastStepSize();

    // Returns the breakpoint tolerance.
    double getBreakpointTol();
    // Sets the breakpoint tolerance.
    void setBreakpointTol(double bptol);

    // returns whether transient analysis is completed
    bool simulationComplete();

    // Gets the size of the restart data (bytes?).
    int restartDataSize( bool pack );

    // Sets the transient calculations parameters
    bool setTranAnalysisParams(const N_UTL_OptionBlock & paramsBlock);

    // Sets the DC sweep calculation parameters.
    bool setDCAnalysisParams(const N_UTL_OptionBlock & paramsBlock);

    // Method to handle OP statements.
    bool setOPAnalysisParams(const N_UTL_OptionBlock & paramsBlock);

    // Sets the STEP calculation parameters.
    bool setSTEPAnalysisParams(const N_UTL_OptionBlock & paramsBlock);

    // Sets the SAVE parameters.
    bool setSaveOptions(const N_UTL_OptionBlock & OB);

    // Sets the DCOP restart parameters.
    bool setDCOPRestartParams(const N_UTL_OptionBlock & OB);

    // Method to register the AC Analysis options.
    bool setACAnalysisParams(const N_UTL_OptionBlock & OB);

    // Method to register the MOR Analysis options.
    bool setMORAnalysisParams(const N_UTL_OptionBlock & OB);

    // Method to register the MOR utility options.
    bool setMOROptions(const N_UTL_OptionBlock & OB);

    // sets a time at which to pause the simulation
    void setPauseTime(double pauseTime);
    // returns time at which to pause the simulation
    double getPauseTime();

    // returns true if the simulation is currently paused
    bool isPaused();

    // signal that simulation is being resumed from previously paused state
    void resumeSimulation();
    // reset the resume flag to false.
    void unset_resumeSimulation();

    // Registers the options block reference.
    bool setTranOptions(const N_UTL_OptionBlock & OB);

    // Method to register the MPDE Analysis options.
    bool setMPDEAnalysisParams(const N_UTL_OptionBlock & OB);

    // Method to register the MPDE utility options.
    bool setMPDEOptions(const N_UTL_OptionBlock & OB);

    // Method to register the HB Analysis options.
    bool setHBAnalysisParams(const N_UTL_OptionBlock & OB);

    // Method to register the HB utility options.
    bool setHBOptions(const N_UTL_OptionBlock & OB);

    // Method to register the linear solver / preconditioning options.
    bool setLinSol(const N_UTL_OptionBlock & OB);

    // Method to register the HB linear solver / preconditioning options.
    bool setHBLinSol(const N_UTL_OptionBlock & OB);

    // Method to register the MPDE utility options.
    bool setTRANMPDEOptions(const N_UTL_OptionBlock & OB);

    // Method to register the sensitivity options.
    bool setSensOptions(const N_UTL_OptionBlock & OB);

    // Registers the TIA parameters block reference.
    bool registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp);

    // Registers the TIA MPDE Interface reference.
    bool registerMPDEInterface(N_TIA_MPDEInterface * tiaMPDEIface_tmp);

    // Registers the linear system pointer.
    bool registerLinearSystem(N_LAS_System * lasSysPtr_tmp);

    // Registers the nonlinear system manager pointer.
    bool registerNLSManager(N_NLS_Manager * nlsMgrPtr_tmp);

    // Registers the nonlinear loader pointer.
    bool registerLoader(N_LOA_Loader * loaderPtr_tmp);

    // Registers the output manager pointer.
    bool registerOutputMgr(N_IO_OutputMgr * outputPtr_tmp);

    // Registers the restart manager pointer.
    bool registerRestartMgr(N_IO_RestartMgr * restartPtr_tmp);

    // Method to register the ?DAE? Loader pointer
    bool registerNonlinearEquationLoader( N_LOA_NonlinearEquationLoader * nonlinearEquationLoaderPtr );

    // Method to register the Device inerfaace pointer
    bool registerDeviceInterface ( N_DEV_DeviceInterface * devInterfacePtr );

    // Method to register the Topology pointer
    bool registerTopology( N_TOP_Topology * topoMgrPtr );

    // Method to register the Restart Manager pointer
    bool registerRestartManager( N_IO_RestartMgr * resMgrPtr );

    // Method to register the Output Manager pointer
    bool registerOutputManager( N_IO_OutputMgr * outMgrPtr );

    // Method to register the Application Builder pointer
    bool registerApplicationBuilder( N_LAS_Builder * appBuilderPtr );

    // Registers the parallel services manager pointer.
    bool registerParallelServices(N_PDS_Manager * pds_tmp);

    // Registers the restart intervals.
    bool registerRestartIntervals();

    // Registers the restart output intervals.
    bool registerOutputIntervals();

    // Registers the elapsed time timer
    bool registerElapsedTimer(N_UTL_Timer *);

    // Writes-out the restart data.
    bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack);

    // Restores the restart data.
    bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

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

    // set/get beginning integration flag
    void setBeginningIntegrationFlag(bool bif);
    bool getBeginningIntegrationFlag();

    // set/get integration method
    void setIntegrationMethod(int im);
    unsigned int getIntegrationMethod();

    // Gets the total time spent in the linear solvers.
    double getTotalLinearSolutionTime() const;

    // Gets the total time spent in the residual load calculations.
    double getTotalResidualLoadTime() const;

    // Gets the total time spent in the Jacobian load calculations.
    double getTotalJacobianLoadTime() const;

    // set the next solution vector pointer.  Needed for NOX...
    bool setNextSolVectorPtr (N_LAS_Vector * solVecPtr);

    // Habanero API mixed signal functions:
    bool provisionalStep (double maxTimeStep, double &currTimeStep);
    void acceptProvisionalStep ();
    void rejectProvisionalStep ();

    // Two-level Newton API functions:
    // Execute the control loop for the set analysis type,
    // for a set number of steps.
    bool runStep
      (const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError);

    void conductanceTest ();
    bool startupSolvers ();
    bool finishSolvers ();

    void homotopyStepSuccess
      ( const vector<string> & paramNames,
        const vector<double> & paramVals);

    void homotopyStepFailure ();

    void stepSuccess (int analysisUpper);
    void stepFailure (int analysisUpper);
    bool getInitialQnorm (N_TIA_TwoLevelError & tle);
    bool getBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes);
    bool startTimeStep (const N_TIA_TimeIntInfo & tiInfo);


    // routines to get/set Dakota run flags and actually run a Dakota iteration
    bool getDakotaRunFlag();
    void setDakotaRunFlag( bool flag );
    int getDakotaIteration();
    void setDakotaIteration( int iterNumber );
    void getTimeIntInfo (N_TIA_TimeIntInfo & tiInfo);

    void initializeTransientModel();
    bool evalTransientModel(
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
        );
    bool evalTransientModelState(
        double t,
        N_LAS_Vector * SolVectorPtr,
        N_LAS_Vector * StaVectorPtr,
        N_LAS_Vector * StoVectorPtr
        );

  protected:

  private :

    bool getInputOPFlag ();

    // Allocate and register the MPDE Manager
    void setupMPDEMgr_();

    // allocate analysis objects:
    void allocateAnalysisObject_();

    void initializeIntegrationProcess_();

    void computeDividedDifferences_();

    // Sets the nonlinear solver solution parameters.
    void setNLSParams_();

    // Two level Newton API:
    bool calledBeforeTwoLevelTran_;
    Teuchos::RefCountPtr<N_ANP_AnalysisBase> twoLevelAnalysisObject_;

    // Habanero mixed-signal API:
    Teuchos::RefCountPtr<N_ANP_AnalysisBase> mixedSignalAnalysisObject_;

    // Queries about output and restart times
    bool outputIntervalSpecified_();
    bool testRestartSaveTime_();

    // Queries to the MPDE Manager.
    // keep there like this so we don't have to
    // expose the MPDE Manager pointer
    void setMPDEFlag( bool flagVal ) { mpdeFlag_ = flagVal; }
    bool getMPDEFlag ();    // "get" function for MPDE flag. (true if not IC)
    bool getMPDEIcFlag();   // get function for MPDE initial condition flag (true if MPDE & IC)
    bool getMPDEStartupFlag();   // True if you have done an initial transient simulation before starting MPDE IC calculation
    bool getWaMPDEFlag ();  // "get" function for WaMPDE flag. (true if not IC)

    // For DCOP restart and/or .SAVE files.
    bool testDCOPOutputTime_();
    bool testSaveOutputTime_();

  public:

    unsigned int breakPointRestartStep;

    // The current time-integration method.
    RefCountPtr<N_TIA_WorkingIntegrationMethod> wimPtr;

    // TIA params object.
    N_TIA_TIAParams tiaParams;

    // Pointer to the ANP Analysis Interface
    RefCountPtr<N_ANP_AnalysisInterface> anaIntPtr;

    // Pointer to the linear system information and containers.
    RefCountPtr<N_LAS_System> lasSysPtr;

    // Pointer to the nonlinear solver manager.
    RefCountPtr<N_NLS_Manager> nlsMgrPtr;

    // Pointer to the nonlinear loader object.
    RefCountPtr<N_LOA_Loader> loaderPtr;

    // Pointer to residual assembler object
    RefCountPtr<N_TIA_Assembler> assemblerPtr;

    // Pointer to the restart manager.
    RefCountPtr<N_IO_RestartMgr> restartPtr;

    // Pointer to the application loader
    RefCountPtr<N_LOA_NonlinearEquationLoader> nonlinearEquationLoaderPtr;

    // Pointer to the device interface
    RefCountPtr<N_DEV_DeviceInterface> devInterfacePtr;

    // Pointer to the topology manager
    RefCountPtr<N_TOP_Topology> topoMgrPtr;

    // Pointer to the output manager
    RefCountPtr<N_IO_OutputMgr> outMgrPtr;

    // Pointer to the applicaiton builder
    RefCountPtr<N_LAS_Builder> appBuilderPtr;

    // Pointer to the parallel services manager.
    RefCountPtr<N_PDS_Manager> pdsMgrPtr;

    ANP_Analysis_Mode analysis;

    bool analysisParamsRegistered;

    // 0 = tranop
    // 1 = transient (time stepping)
    // 2 = dc sweep.
    int currentMode_;

    bool firstTime;
    double oldPercentComplete;
    double startSimTime;

    N_MPDE_Manager * getMPDEManager() const { return &*mpdeMgrPtr_; }
    N_TIA_DataStore * getTIADataStore() const { return &*dsPtr_; }
    Teuchos::RefCountPtr<const N_ANP_AnalysisBase> getAnalysisObject() const { return primaryAnalysisObject_; }

    void silenceProgress();
    void enableProgress();

  protected:
  private:

#if 0
    // Current time-integration method flag.
    unsigned int integrationMethod_;
#endif

    // Switch the integration flag.
    bool switchIntegrator_;

    // DC Operating Point flag.
    bool initializeAllFlag_;      // true if the initializeAll function has been
                                  // called once.

    //
    // 05/26/09 Coffey,Schiek,Mei:  We considered the data members down to this point for moving to N_ANP_Transient
    //
    double startTRANtime;

    bool stepLoopFlag_;           // true if there is an external
                                  // parameter stepper loop around the dcop or
                                  // transient simulation loops.

    bool stepLoopInitialized_;    // true if the step loop has been set up.
    bool dcLoopInitialized_;      // true if the dc sweep loop has been set up.
    bool gui_;                    // command line arg -gui is present


    bool daeStateDerivFlag_;   // true if running new-DAE and need
                              // state derivative.  true by default.
                              // If set to true, it breaks MPDE.

    bool initializeSolvers_mixedSignal_;

    int dcLoopSize_;

    bool sweepSourceResetFlag_;

    // Flag to decide whether to print progress
    bool progressFlag_;

    // Xyce timing utility for timing the transient simulation CPU time.
    RefCountPtr<N_UTL_Timer> xyceTranTimerPtr_;

    // Xyce timing utility for timing elapsed run time
    RefCountPtr<N_UTL_Timer> elapsedTimerPtr_;

    double solverStartTime_;

    // pointer to the output manager adapter
    RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr_;

    // Pointer to the TIA data store object.
    RefCountPtr<N_TIA_DataStore> dsPtr_;

    // Pointer to the TIA step-error control object.
    RefCountPtr<N_TIA_StepErrorControl> secPtr_;

    bool dakotaRunFlag_;
    int dakotaIterationNumber_;

    // for .SAVE and/or DCOP restart.
    double saveTime_;
    bool saveTimeGiven_;
    bool saveFlag_;
    bool savedAlready_;
    bool dcopRestartFlag_;

    // .OP flag(s)
    bool dotOpSpecified_;

    // output and restart interval info
    double initialOutputInterval_;
    vector < pair < double, double > > outputIntervals_;
    double nextOutputTime_;

    double initialRestartInterval_;
    vector < pair < double, double > > restartIntervals_;
    double nextRestartSaveTime_;

    // for HB, MPDE, or any other block analysis type
    bool blockAnalysisFlag_;
    bool hbFlag_;
    bool mpdeFlag_;

    // sensitivity flag(s)
    bool sensFlag_;

    // for handling expressions given for a time dependent max time step
    //bool maxTimeStepExpressionGiven_;
    //string maxTimeStepExpressionAsString_;

    // ref counted pointers for various analyses. Not all are used in every simulation
    Teuchos::RefCountPtr<N_ANP_AnalysisBase> analysisObject_;
    Teuchos::RefCountPtr<N_ANP_AnalysisBase> stepAnalysisTarget_;
    Teuchos::RefCountPtr<N_ANP_AnalysisBase> dakotaAnalysisTarget_;
    Teuchos::RefCountPtr<N_ANP_AnalysisBase> primaryAnalysisObject_;

    // Friended class for lookup and set of restart data.
    friend class N_TIA_StepErrorControl;
    friend class N_ANP_AnalysisInterface;
    friend class N_TIA_DAE_Assembler;
    friend class N_TIA_Assembler;
    friend class N_ANP_AnalysisBase;
    friend class N_ANP_Transient;
    friend class N_ANP_MPDE;
    friend class N_ANP_AC;
    friend class N_ANP_MOR;
    friend class N_ANP_HB;
    friend class N_ANP_DCSweep;
    friend class N_ANP_Step;
    friend class N_ANP_Dakota;

    // command line object
    N_IO_CmdParse & commandLine_;

    // MPDE stuff:
    RefCountPtr<N_MPDE_Manager> mpdeMgrPtr_;
    RefCountPtr<N_TIA_MPDEInterface> tiaMPDEIfacePtr_;

    // parameter objects for each analysis type.  These are saved until the
    // analysis object of choice is allocated.
    N_UTL_OptionBlock tranParamsBlock;
    N_UTL_OptionBlock opParamsBlock;
    N_UTL_OptionBlock acParamsBlock;
    N_UTL_OptionBlock morParamsBlock;

    // Different block passed in for each sweep variable, so need to have
    // these containers be vectors.
    vector<N_UTL_OptionBlock> dcParamsBlockVec;
    vector<N_UTL_OptionBlock> stepParamsBlockVec;

    N_UTL_OptionBlock mpdeParamsBlock;
    N_UTL_OptionBlock hbParamsBlock;
    N_UTL_OptionBlock hbOptionsBlock;
    N_UTL_OptionBlock hbLinSolBlock;
    N_UTL_OptionBlock linSolBlock;
    N_UTL_OptionBlock dakotaParamsBlock;
};

#endif
