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
// Filename      : $RCSfile: N_TIA_TIAParams.C,v $
//
// Purpose       : This file implements the class associated with all user
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
// Revision Number: $Revision: 1.121 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif

#include <string>

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisInterface.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_TIAParams.h>
#include <N_ERH_ErrorMgr.h>

// -------- Forward declaration -----------
#include <N_IO_CmdParse.h>
 
//-----------------------------------------------------------------------------
// Function      : N_TIA_TIAParams::N_TIA_TIAParams
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------

N_TIA_TIAParams::N_TIA_TIAParams(N_IO_CmdParse & cp)
  :
  commandLine(&cp),
  tStart(0.0),
  tStartGiven(false),
  initialTime(0.0),
  finalTime(0.0),
  pauseTime(0.0),
  pauseSetAtZero(false),
  integrationMethod(7),
  userSpecified_startingTimeStep(1.0e-10),
  maxTimeStep(1.0e+99),
  maxTimeStepGiven(false),
  userSpecMinTimeStep(0.0),
  userSpecMinTimeStepGiven(false),
  constantStepSize(false),
  useDeviceTimeStepMax(true),
  restartingIntegrationFromSSS(false),
  newBPStepping(true),
//  newBPStepping(false),
  minTimeStepsBP(10),
//  minTimeStepsBPGiven(true),
  minTimeStepsBPGiven(false),
//  newLte(false),
  newLte(true),
  errorAnalysisOption(0),
  NLmin(3),
  NLmax(8),
  TimeStepLimitedbyBP(false),
  delmax(1.0e+99),
  delmaxGiven(false),
  errorAnalysisOptionResetIter(0),
  timestepsReversal(false), 
  testFirstStep(false),
  solutionSize(0),
  stateSize(0),
  relErrorTol(1.0e-3),
//  relErrorTol(1.0e-2),
  absErrorTol(1.0e-6),
  relErrorTolGiven(false),
  errTolAcceptance(1.0),
//  errTolAcceptance(0.6859),
  NOOP(false),
  resume(false),
//#ifdef Xyce_MPDE
//  sweepSteps(1),
//#else
  sweepSteps(0),
//#endif // Xyce_MPDE
  doubleDCOPStep(0),
  firstDCOPStep(0),
  lastDCOPStep(1),
  doubleDCOPAll(false),
  exitTime(0.0),
  exitStep(-1),
#ifdef Xyce_DEBUG_TIME
  debugLevel(1),
#else
  debugLevel(0),
#endif
  bpEnable(true),
  restartTimeStepScale(0.005),
//  nlNearConvFlag(true),
  nlNearConvFlag(false),
  nlSmallUpdateFlag(true),
  jacLimitFlag(false),
  jacLimit(1.0e+17),
  maxOrder(5),
  minOrder(1),
//  freq(0.0),
//  freq2(0.0),
//  freq2Given(false),
  freqGiven(false),
  type("DEC"),
  np(10.0),
  fStart(1.0),
  fStop(1.0), 
  ROMsize(-1),
  morMethod("PRIMA"),
  morSaveRedSys(false),
  morCompOrigTF(false), 
  morCompRedTF(false), 
  morCompType("DEC"), 
  morCompNP(10),
  morCompFStart(1.0), 
  morCompFStop(1.0), 
  morExpPoint(0.0), 
  morScaleFactor(1.0),
  morScaleType(0),
  morScaleFactor1(0.01), 
  morSparsificationType(0),
  outputInterpMPDE(true),
  interpOutputFlag(true),
  condTestFlag(false),
  saveTimeStepsFlag(false),
  passNLStall(false),
  minTimeStepRecoveryCounter(0),
  fastTests(false),
  voltZeroTol(1.0e-6),
  currZeroTol(1.0e-6),
  historyTrackingDepth(25)
{
  scalarTolerances = true;

  freqs.clear();
#ifdef Xyce_DEBUG_TIME
  // override default debug level based on command line options
  if ( commandLine->argExists( "-tdl" ) )
  {
    debugLevel = atoi( commandLine->getArgumentValue( "-tdl" ).c_str() );
  }
#endif 
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TIAParams::N_TIA_TIAParams
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------

N_TIA_TIAParams::~N_TIA_TIAParams()
{

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TIAParams::printParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233.
// Creation Date : 7/12/01
//-----------------------------------------------------------------------------
void N_TIA_TIAParams::printParams(std::ostream &os, int analysis)
{
  os << "\n" << std::endl;
  os << Xyce::section_divider << std::endl;
  os << "\n***** Time Integration solver options:\n" << std::endl;

  if (analysis == TRANSIENT)
  {
    os << "\tAnalysis:\t\t\tTRANSIENT" << std::endl
                  << "\tInitial Time (sec):\t\t" <<  initialTime << std::endl
                 << "\tFinal Time (sec):\t\t" <<  finalTime << std::endl
                 << "\tStarting Time Step(sec):\t" <<  userSpecified_startingTimeStep << std::endl
                 << "\tRestart Time Step Scale:\t" <<  restartTimeStepScale << std::endl
                 << "\tError Analysis option:\t" <<  errorAnalysisOption << std::endl
                 << "\ttStart (initial outputTime):\t" <<  tStart << std::endl;

    switch (integrationMethod)
    {
    case TIAMethod_BACKWARD_EULER:
      os << "\tTime Integration method:\tBACKWARD EULER" << std::endl;
      break;

    case TIAMethod_BACKWARD_DIFFERENTIATION_2:
      os << "\tTime Integration method:\tBACKWARD DIFFERENTIATION ORDER 2" << std::endl;
      break;

    case TIAMethod_BACKWARD_DIFFERENTIATION_15:
      os << "\tTime Integration method:\tBACKWARD DIFFERENTIATION ORDER 15" << std::endl;
      break;

    case TIAMethod_GEAR_12:
      os << "\tTime Integration method:\tGEAR 12" << std::endl;
      break;

    case TIAMethod_ONESTEP:
      os << "\tTime Integration method:\tONESTEP" << std::endl;
      break;

    case TIAMethod_TRAPEZOIDAL:
      os << "\tTime Integration method:\tTRAPEZOIDAL" << std::endl;
      break;

    case TIAMethod_VARIABLE_THETA:
      os << "\tTime Integration method:\tVARIABLE THETA" << std::endl;
      break;

    case TIAMethod_A_CONTRACTIVE_2:
      os << "\tTime Integration method:\tA-CONTRACTIVE order 2" << std::endl;
      break;
    }

    os << (constantStepSize ? "\tUsing Constant Step Size" : "\tUsing Variable Step Size") << std::endl
                 << (useDeviceTimeStepMax ? "\tUsing Device specified maximum stepsize" : "\tNOT using Device specified maximum stepsize") << std::endl
                 << (fastTests ? "\tTime integration FastTests is ON" : "\tTime integration FastTests is OFF") << std::endl
                 << (nlNearConvFlag ? "\tNL Near Convergence Flag is ON" : "\tNL Near Convergence Flag is OFF") << std::endl
                 << (passNLStall ? "\tNL Pass Non-linear Stalls is ON" : "\tNL Pass Non-linear Stalls is OFF") << std::endl;
  }
  else
  {
    os << "\tAnalysis:\t\t\tDC SWEEP" << std::endl
                 <<"\tTotal DC Sweep Steps:\t\t" << sweepSteps << std::endl;

  }

  os << "\tabsErrorTol:\t\t\t" <<  absErrorTol  << std::endl
               << "\trelErrorTol:\t\t\t" <<  relErrorTol  << std::endl
               << "\texitTime:\t\t\t" <<  exitTime  << std::endl
               << "\texitStep:\t\t\t" <<  exitStep  << std::endl
               << "\tdebugLevel:\t\t\t" <<  debugLevel << std::endl
               << "\tMaximum Order:\t\t\t" <<  maxOrder << std::endl
               << "\tMinimum Order:\t\t\t" <<  minOrder << std::endl
               << "\tInterpolated Output Flag:\t\t " << (interpOutputFlag ? "true": "false") << std::endl
               << "\tConductance Test Flag:\t\t" << (condTestFlag  ? "true": "false") << std::endl;
  
  if (condTestDeviceNames.empty())
  {
    os <<
		       "\tConductance Test Device Name List is:\t\tEMPTY" << std::endl;
  }
  else
  {
    os << "\tConductance Test Device Name List contains:  " << std::endl;
    for (std::list< std::string >::iterator it = condTestDeviceNames.begin(); it != condTestDeviceNames.end(); ++it)
    {
      os << " \"" << *it << "\"";
    }
  }
  os << Xyce::section_divider << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TIAParams::operator=
// Purpose       : "=" operator.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/17/05
//-----------------------------------------------------------------------------
N_TIA_TIAParams & N_TIA_TIAParams::operator=(const N_TIA_TIAParams & right)
{
  commandLine   = right.commandLine;

  tStart = right.tStart;
  initialTime = right.initialTime;
  finalTime = right.finalTime;
  pauseTime = right.pauseTime;
  pauseSetAtZero = right.pauseSetAtZero;
  userSpecified_startingTimeStep = right.userSpecified_startingTimeStep;
  maxTimeStep = right.maxTimeStep;
  maxTimeStepGiven = right.maxTimeStepGiven;
  userSpecMinTimeStep = right.userSpecMinTimeStep;
  userSpecMinTimeStepGiven = right.userSpecMinTimeStepGiven;
  exitTime = right.exitTime;
  constantStepSize = right.constantStepSize;
  useDeviceTimeStepMax = right.useDeviceTimeStepMax;
  errorAnalysisOption = right.errorAnalysisOption;
  restartingIntegrationFromSSS = right.restartingIntegrationFromSSS;
  NOOP = right.NOOP;
  resume = right.resume;
  bpEnable = right.bpEnable;
  restartTimeStepScale = right.restartTimeStepScale;
  integrationMethod = right.integrationMethod;
  sweepSteps = right.sweepSteps;
  solutionSize = right.solutionSize;
  stateSize = right.stateSize;
  doubleDCOPStep = right.doubleDCOPStep;
  firstDCOPStep = right.firstDCOPStep;
  lastDCOPStep = right.lastDCOPStep;
  exitStep = right.exitStep;
  doubleDCOPAll = right.doubleDCOPAll;
  relErrorTol = right.relErrorTol;
  absErrorTol = right.absErrorTol;
  errTolAcceptance = right.errTolAcceptance;
  scalarTolerances = right.scalarTolerances;
  debugLevel = right.debugLevel;
  nlNearConvFlag = right.nlNearConvFlag;
  nlSmallUpdateFlag = right.nlSmallUpdateFlag;
  jacLimitFlag = right.jacLimitFlag;
  jacLimit = right.jacLimit;
  maxOrder = right.maxOrder;
  minOrder = right.minOrder;
  freqs = right.freqs;
  freqGiven = right.freqGiven;
  interpOutputFlag = right.interpOutputFlag;
  condTestFlag = right.condTestFlag;
  condTestDeviceNames = right.condTestDeviceNames;
  saveTimeStepsFlag = right.saveTimeStepsFlag;
  passNLStall = right.passNLStall;
  fastTests = right.fastTests;
  voltZeroTol = right.voltZeroTol;
  currZeroTol = right.currZeroTol;

  return *this;
}


