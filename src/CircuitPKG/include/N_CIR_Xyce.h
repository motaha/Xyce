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
// Filename       : $RCSfile: N_CIR_Xyce.h,v $
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
// Revision Number: $Revision: 1.90.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:31 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_CIRCUIT_h
#define Xyce_N_CIR_CIRCUIT_h

// ---------- Standard Includes ----------

#include <string>
#include <map>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

// ---------- Forward Declarations ----------

class N_TOP_Topology;

namespace Xyce
{
  namespace Topology
  {
    class Manager;
  }
}

class N_LAS_System;
class N_LAS_Builder;

class N_ANP_AnalysisInterface;
class N_TIA_TimeIntInfo;
class N_TIA_TwoLevelError;

class N_PDS_Manager;
class N_PDS_ParMap;

class N_ERH_ErrorMgr;

class NetlistImportTool;

class N_IO_RestartMgr;
class N_IO_PkgOptionsMgr;

class N_NLS_Manager;

class N_LOA_Loader;
class N_LOA_CktLoader;
class N_LOA_NonlinearEquationLoader;
class N_LOA_LoaderMgr;

class N_UTL_Timer;

class N_UTL_BreakPoint;

class N_LAS_Vector;
class N_LAS_Matrix;
class Epetra_CrsGraph;

#include <N_IO_CmdParse.h>

//-----------------------------------------------------------------------------
// Class         : N_CIR_Xyce
// Purpose       : This is the main "top level" class for Xyce.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
class N_CIR_Xyce
{
  // Methods
 public:

  // Default constructor
  N_CIR_Xyce(
#ifdef Xyce_PARALLEL_MPI
             MPI_Comm * comm = 0
#endif
            );

  // Default destructor (virtual so we derive other classes from this properly)
  virtual ~N_CIR_Xyce();

  // These are all the API calls that we are suppose to be making available
  // for external programs and/or other objects

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::setNetlistParameters
  // Purpose       : This passes a vector of pairs "key", "value" that will
  //                 be substituted during the processing of the netlist.  This
  //                 more easily allows Dakota to change any netlist parameter
  //                 during netlist setup.
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical and MEMS Modeling
  // Creation Date : 10/9/2008
  //---------------------------------------------------------------------------
  void setNetlistParameters( const vector< pair< string, string > > & externalParams );


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
  void setOutputFileSuffix( const string newSuffix );

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::run
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 02/19/01
  //---------------------------------------------------------------------------
  bool run(int iargs_tmp, char *cargs_tmp[]);

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::initialize
  // Purpose       : To initialize Xyce to be driven by the SAX
  //                 simulation backplane. This includes the following:
  //                    Set up and register the Parallel Manager.
  //                    Parse the command-line arguments.
  //                    Redirect the output stream of processor 1,
  //                    if requested.
  //                    Read in the Netlist.
  //                    Allocate and register the external packages.
  //                    Set up the representation of the circuit topology.
  //                    Set up the matrix structures.
  //                    Initialize the solvers.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/28/03
  //---------------------------------------------------------------------------
  bool initialize(int iargs_tmp, char *cargs_tmp[]);

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::runSimulation
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 5/27/00
  //---------------------------------------------------------------------------
  bool runSimulation();

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::getDACDeviceNames
  // Purpose       : Gets the (stripped) names of the DAC devices
  //                 in the circuit.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 06/13/03
  //---------------------------------------------------------------------------
  bool getDACDeviceNames(std::vector< string >& dacNames);

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::getADCMap
  // Purpose       : Gets the (stripped) names of the ADC devices
  //                 in the circuit(as key of map) and map of parameters
  //                 (keyed by parameter name) for each device
  // Special Notes :
  // Scope         : public
  // Creator       : Tom Russo, SNL, Component Information and Models
  // Creation Date : 05/06/2004
  //---------------------------------------------------------------------------
  bool getADCMap(map<string,map<string,double> >& ADCMap);

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::updateTimeVoltagePairs
  // Purpose       : Update the DAC devices in a circuit by adding the set
  //                 of time and voltage pairs built up on the "digital side"
  //                 since the last update and by removing the time-voltage
  //                 pairs for times that pre-date the given simulation time.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 06/10/03
  //---------------------------------------------------------------------------
  bool updateTimeVoltagePairs(
        map< string, vector< pair<double,double> >* > const&
        timeVoltageUpdateMap);

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::getTimeVoltagePairs
  // Purpose       : query the DAC devices in a circuit for the set
  //                 of time and voltage pairs
  // Special Notes :
  // Scope         : public
  // Creator       : Tom Russo, SNL ComponentInformation and Models
  // Creation Date : 05/10/2004
  //---------------------------------------------------------------------------
  bool getTimeVoltagePairs(
        map< string, vector< pair<double,double> > > &
        timeVoltageUpdateMap);

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
  bool setADCWidths(map< string, int > const& ADCWidthMap);

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
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/28/03
  //---------------------------------------------------------------------------
  bool simulateUntil(double requestedUntilTime, double& completedUntilTime);

  //---------------------------------------------------------------------------
  // Function      : N_CIR_Xyce::finalize
  // Purpose       : To clean up after driving Xyce with the SIMBUS
  //                 simulation backplane. This includes the following:
  //                    Free any dynamically allocated memory...
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/29/03
  //---------------------------------------------------------------------------
  bool finalize();

  void reportTotalElapsedTime ();

  void readExternalParamsFromFile( string filename, string & iterationSuffix, vector< pair< string, string > > & paramList );

  // accessor routines for external program controlled simulations that
  // need a response function
  bool registerResponseVars (string objString, RCP<vector< double > > varVectorPtr );
  void finalizeResponseVars();

  // report on whether simulation is finished or not
  bool simulationComplete();

  // new mixed-signal functions:
  bool provisionalStep
    (double maxtimeStep,
     double &timeStep,
     map< string, vector< pair<double,double> > > & timeVoltageUpdateMap);

  void acceptProvisionalStep();
  void rejectProvisionalStep();

  double getFinalTime();
  double getTime();

  // 2-level, power node functions:
  bool startupSolvers ();
  bool finishSolvers ();

  void homotopyStepSuccess
    (const vector<string> & paramNames,
     const vector<double> & paramVals);

  void homotopyStepFailure ();

  void stepSuccess (int analysis);
  void stepFailure (int analysis);
  bool getBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes);
  bool updateStateArrays ();
  bool startTimeStep (const N_TIA_TimeIntInfo & tiInfo);
  bool startTimeStep (const N_DEV_ExternalSimulationData & ext_data);
  bool endTimeStep (N_DEV_ExternalSimulationData & ext_data);
  void enable_lte_analysis();
  bool setInternalParam (string & name, double val);

  bool getInitialQnorm (N_TIA_TwoLevelError & tle);

  bool simulateStep
      ( const N_DEV_SolverState & solState,
        const map <string,double> & inputMap,
        vector<double> & outputVector,
        vector< vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError);

  // 05/27/09 ModelEvaluator Interface functions
  void initializeTransientModel();

  bool evalTransientModel (
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

  bool evalTransientModelState (
     double t,
     N_LAS_Vector * SolVectorPtr,
     N_LAS_Vector * StaVectorPtr,
     N_LAS_Vector * StoVectorPtr
    );

  void getMapsAndGraphs (
    RCP<N_PDS_ParMap> & x_map,
    RCP<N_PDS_ParMap> & x_map_ognd,
    RCP<N_PDS_ParMap> & s_map,
    RCP<N_PDS_ParMap> & store_map,
    RCP<Epetra_CrsGraph> & dFdx_graph,
    RCP<Epetra_CrsGraph> & dFdx_graph_ognd,
    RCP<Epetra_CrsGraph> & dQdx_graph,
    RCP<Epetra_CrsGraph> & dQdx_graph_ognd
    );

  std::vector<std::string> getVariableNames();

  bool simulateStep
      ( const N_DEV_ExternalSimulationData & ext_data,
        const map <string,double> & inputMap,
        vector<double> & outputVector,
        vector< vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError);

 private:

  bool setupParMgr_(int iargs, char *cargs[]);
  bool doAllocations_();
  bool doInitializations_();
  bool doRegistrations_();
  bool doDeAllocations_();
  bool setUpTopology_();
  bool setUpMatrixStructure_();
  bool runSolvers_();

 public:

protected:

private:

  // Attributes
public:

protected:
  // We put the device interface pointer as a protected method so we can
  // access it from derived classes without kludges.

  // Device manager
  N_DEV_DeviceInterface          * devIntPtr_;

private:

  // Topology manager
  Xyce::Topology::Manager        * topMgrPtr_;
  N_TOP_Topology                 * topPtr_;

  // Linear algebra system
  N_LAS_System                   * lasSysPtr_;

  // Linear algebra system builder
  N_LAS_Builder                  * lasBuilderPtr_;

  // Time integration manager
  N_ANP_AnalysisInterface        * anaIntPtr_;

  // Nonlinear solver manager
  N_NLS_Manager                  * nlsMgrPtr_;

  // Loader (residual and Jacobian) manager
  N_LOA_LoaderMgr                * loaderMgrPtr_;

  // Loader (residual and Jacobian) pointer
  //N_LOA_Loader                   * loaderPtr_;

  // Circuit loader (residual and Jacobian) pointer
  N_LOA_CktLoader                * cktLoaderPtr_;

  // Dae loader (residual and Jacobian) pointer
  N_LOA_NonlinearEquationLoader                * nonlinearEquationLoaderPtr_;

  // Parallel distribution manager
  N_PDS_Manager                  * parMgrPtr_;

  // Netlist import tool.
  NetlistImportTool              * netlistImportToolPtr_;

  // Output manager
  N_IO_OutputMgr                 * outMgrPtr_;

  // Restart manager
  N_IO_RestartMgr                * resMgrPtr_;

  // Xyce solver timing utility
  N_UTL_Timer                    * XyceTimerPtr_;

  // Elapsed time from beginning of run
  N_UTL_Timer                    * ElapsedTimerPtr_;

  // package options manager
  RCP<N_IO_PkgOptionsMgr>          pkgOptionsMgrPtr_;

  int iargs;
  char **cargs;

  // if the user is providing an externa file with parameters in it,
  // (as in x=1, y=4) the filename may also have an iteration / realization suffix
  // such as params.in.1234.  iterationSuffix_ will attempt to capture that unique
  // identifier while externalNetlistParams_ will map params to values.  The values
  // are treated as strings at this point but later converted to numbers in parsing.
  string iterationSuffix_;
  vector< pair< string, string> > externalNetlistParams_;

  bool isSerialFlag_;
  bool procZeroFlag_;
  bool multiThreading_;
  int numThreads_;

  bool initializeAllFlag_;

  N_IO_CmdParse commandLine;

#ifdef Xyce_PARALLEL_MPI
  MPI_Comm * commPtr_;
#endif

};

#endif

