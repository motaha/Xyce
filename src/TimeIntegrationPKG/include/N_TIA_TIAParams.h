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
// Filename      : $RCSfile: N_TIA_TIAParams.h,v $
//
// Purpose       : This file identifies the class associated with all user
//                 specified parameters which relate to the time integration
//                 algorithms and problem definition.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.72.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_PARAMS_H
#define Xyce_N_TIA_PARAMS_H

// ----------   Standard Includes   ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>

// ---------- Forward Declarations ----------

class N_TIA_DataStore;
class N_IO_CmdParse;

//-----------------------------------------------------------------------------
// Class         : N_TIA_TIAParams::N_TIA_TIAParams
// Purpose       : This is a class that sets up time integration information,
//                 associated with user related input only, that defines the
//                 problem setup.
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class N_TIA_TIAParams
{
  public:

    // tStart is the time at which printing starts. 
    double tStart;
    bool tStartGiven;

    // Beginning time for the time-integrator.
    double initialTime;

    // Last time step value.
    double finalTime;

    // Time step value at which to "pause" the simulation.
    double pauseTime;
    
    // flag used to indicate that a puase was specifically set at time zero 
    // and thus should not be ignored.
    bool pauseSetAtZero;

    // Initial time-step value.
    double userSpecified_startingTimeStep;

    // User-specified maximum time step.
    double maxTimeStep;
    bool maxTimeStepGiven;

    // User-specified mininum number of steps between breakpoints
    int minTimeStepsBP;
    bool minTimeStepsBPGiven;
    
    double userSpecMinTimeStep;
    bool userSpecMinTimeStepGiven;

    // The time-integrator will exit when it exceeds this time.
    double exitTime;

    // Constant time step-size integration flag.
    bool constantStepSize;

    // Flag to set if we are letting the devices impose a time step maximum
    bool useDeviceTimeStepMax;

    // error analysis option
    int errorAnalysisOption;
    
    // iteration count down to reset errorAnalysisOption to zero
    int errorAnalysisOptionResetIter;

    // Restarting integration from a steady-state solution flag.
    bool restartingIntegrationFromSSS;

    // No operating point flag, if true, the operating point calculation will
    // be skipped
    bool NOOP;

    // Continue a paused calculation.
    bool resume;

    // Enable breakpoints flag.
    bool bpEnable;

    // This is the extent to which the time step is 
    // scaled coming out of a breakpoint.
    double restartTimeStepScale;

  
    //iteration count algorithm parameters
    int NLmin;
    int NLmax;
    bool TimeStepLimitedbyBP;
    double delmax;
    bool delmaxGiven;
    bool timestepsReversal;
    bool testFirstStep;
    
    // Time-integration method flag.
    unsigned int integrationMethod;
    unsigned int sweepSteps;
    unsigned int solutionSize;
    unsigned int stateSize;
    
    bool newBPStepping; 
    bool newLte; 
    
    int doubleDCOPStep;
    int firstDCOPStep;
    int lastDCOPStep;

    // The time-integrator will exit after taking this many steps.
    int exitStep;

    bool doubleDCOPAll;

    // Error Tolerances:

    // Relative error tolerance.  This value should be selected to be 10^(-(m+1))
    // where m is the desired number of significant digits in the solution.
    double relErrorTol;
    
    bool relErrorTolGiven;
    
    // Absolute error tolerance.  In general, this should be selected to be small
    // enough such that any solution value smaller than this tolerance will be
    // considered "insignificant".
    double absErrorTol;

    // Error acceptance tolerance.
    double errTolAcceptance;
    
    bool scalarTolerances;

    // Vector containing user-specified time-integration break points.
    //vector < double > userSpecBreakPoints;

    // Debug Level
    int debugLevel;

    // NL solver "near convergence" flag.
    bool nlNearConvFlag;

    // NL solver "small update" flag.
    bool nlSmallUpdateFlag;

    // Jacobian limit - to prevent capacitive spiral of death!
    bool jacLimitFlag;
    double jacLimit;

    // Maximum/Minimum order desired (BackwardDifferentiation15 specific)
    int maxOrder;
    int minOrder;

    // flag for interpolating the MPDE output.
    // Sometimes, the interpolation is the hard part.
    bool outputInterpMPDE;

    // Period of oscillation used by HB
    double freq;
    bool freqGiven; 
     
    // AC  
    string type; 
    double np;  
    double fStart; 
    double fStop; 
   
    // MOR
    int ROMsize;
    std::string morMethod;
    bool morSaveRedSys;
    bool morCompOrigTF;
    bool morCompRedTF;
    std::string morCompType;
    int morCompNP;
    double morCompFStart;
    double morCompFStop;
    double morExpPoint;
    double morScaleFactor;
    int morScaleType;
    double morScaleFactor1;
    int morSparsificationType;
 
    // flag for interpolated output
    bool interpOutputFlag;

    // flag for conductance test
    bool condTestFlag;

    // flag to save timestpes in data store for later use
    bool saveTimeStepsFlag;
    
    // option to pass some non-linear solver failures
    bool passNLStall;
    
    // if we have rejected several time steps and are about to fail due to a 
    // time step too smal error, we'll go back and accept the step that had
    // the minimum estimated error over tol if this flag is true.
    // (set by the user in the netlist via .options timeint mintimesteprecovery=<int>)
    int minTimeStepRecoveryCounter;
    
    // if fastTests == true, then we'll consider voltages and currents
    // below volAbsTol and currAbsTol to be converged in calculating the 
    // weights.
    bool fastTests;
    double voltZeroTol;
    double currZeroTol;

    // Xyce tracks how the last few timesteps have gone for error reporting
    // if Xyce has to exit.  This var specifies how many time steps (both
    // passing and failing) to keep in its history.  If it's zero, then
    // history tracking is off.
    int historyTrackingDepth;

    // list of names for conductance test
    std::list< std::string > condTestDeviceNames;

    // commandline parser object
    N_IO_CmdParse & commandLine;

  protected:
  private :

  public :
    // Default constructor
    N_TIA_TIAParams(N_IO_CmdParse & cp);

    // Destructor
    ~N_TIA_TIAParams();

    N_TIA_TIAParams & operator=(const N_TIA_TIAParams & right);

#ifdef Xyce_VERBOSE_TIME
    // Print out time-integration parameters.
    void printParams(int analysis);
#endif

};

#endif // Xyce_N_TIA_TIAParams_H
