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
// Filename       : $RCSfile: N_NLS_NLParams.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/13/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.74 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include <N_NLS_NLParams.h>
#include <N_ERH_ErrorMgr.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_DampedNewton.h>
#include <N_UTL_Param.h>
#include <N_UTL_OptionBlock.h>
#include <N_NLS_TwoLevelNewton.h>
#include <N_IO_CmdParse.h>

// ---------- Forward Declarations ----------

// ---------- Static Initializations ----------

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::N_NLS_NLParams
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/09/00
//-----------------------------------------------------------------------------
N_NLS_NLParams::N_NLS_NLParams(AnalysisMode mode, N_IO_CmdParse & cp)
  : commandLine_(&cp),
    modeToggled_(true),
#ifdef Xyce_VERBOSE_NONLINEAR
    printParamsFlag_(true),
#endif
#ifdef Xyce_DEBUG_NONLINEAR
    debugLevel_(1),
    debugMinTimeStep_(0),
    debugMaxTimeStep_(N_UTL_MachineDependentParams::IntMax()),
    debugMinTime_(0.0),
    debugMaxTime_(N_UTL_MachineDependentParams::DoubleMax()),
    screenOutputFlag_(false),
#endif
    analysisMode_(mode)
{

//  lasSolverPtr = 0;
  // Set defaults

  // Set default update norm tolerance
  resetDeltaXTol();

  // Set default small update norm tolerance
  resetSmallUpdateTol ();

  // Set default residual norm tolerance
  resetRHSTol();

  // Set default absolute tolerance value for use in weighted norm
  resetAbsTol();

  // Set default relative tolerance value for use in weighted norm
  resetRelTol();

  resetEnforceDeviceConvFlag ();

  resetSearchMethod();
  resetDirection();
  resetNLStrategy();
  resetMaxNewtonStep();
  resetMaxSearchStep();
  resetForcingFlag();
  resetForcingTerm();
  resetNormLevel();
  resetLinearOpt();
  resetConstraintBT();
  resetGlobalBTMax();
  resetGlobalBTMin();
  resetGlobalBTChange();

  // Set the default parameters for transient, if the specified mode
  // is TRANSIENT.
  // The defaults set in the NLParams constructor are for DC_OP,
  // so the transient ones should be reset here.
  if  (mode==TRANSIENT)
  {
    setNLStrategy(NEWTON);
    setSearchMethod(FULL);
    setMaxSearchStep(2);
    setMaxNewtonStep(20);
    setDeltaXTol(0.33);
    setForcingFlag(false);
    setAbsTol(1.0e-06);
    setRelTol(1.0e-02);
    setRHSTol(1.0e-02);
  }

  if  (mode == HB_MODE)
  {
    setNLStrategy(NEWTON);
    setSearchMethod(FULL);
//    setMaxSearchStep(2);
//    setMaxNewtonStep(10);
//    setDeltaXTol(0.33);
//    setForcingFlag(false);
    setAbsTol(1.0e-9);
//    setRelTol(1.0e-02);
    setRHSTol(1.0e-4);
  }


  setCmdLineOptions ();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::N_NLS_NLParams
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/09/00
//-----------------------------------------------------------------------------
N_NLS_NLParams::N_NLS_NLParams(const N_NLS_NLParams & right)
  : commandLine_(right.commandLine_),
    absTol_(right.absTol_),
    relTol_(right.relTol_),
    INForcingFlag_(right.INForcingFlag_),
    eta_(right.eta_),
    normLevel_(right.normLevel_),
    linearOptimization_(right.linearOptimization_),
    modeToggled_(right.modeToggled_),
#ifdef Xyce_VERBOSE_NONLINEAR
    printParamsFlag_(right.printParamsFlag_),
#endif
    analysisMode_(right.analysisMode_),
    maxNewtonStep_(right.maxNewtonStep_),
    maxSearchStep_(right.maxSearchStep_),
    nlStrategy_(right.nlStrategy_),
    searchMethod_(right.searchMethod_),
    direction_(right.direction_),
    deltaXTol_(right.deltaXTol_),
    RHSTol_(right.RHSTol_),
    constraintBT_(right.constraintBT_),
    globalBTMax_(right.globalBTMax_),
    globalBTMin_(right.globalBTMin_),
    globalBTChange_(right.globalBTChange_)
#ifdef Xyce_DEBUG_NONLINEAR
    ,
    debugLevel_(right.debugLevel_),
    debugMinTimeStep_(right.debugMinTimeStep_),
    debugMaxTimeStep_(right.debugMaxTimeStep_),
    debugMinTime_(right.debugMinTime_),
    debugMaxTime_(right.debugMaxTime_)
#endif
{
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::~N_NLS_NLParams
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/09/00
//-----------------------------------------------------------------------------
N_NLS_NLParams::~N_NLS_NLParams()
{
//  if (lasSolverPtr != 0) delete lasSolverPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setOptions
// Purpose       : This function takes an .options statement option block,
//                 and uses it to set the various parameters of
//                 the NLParams class.
//
// Special Notes : This was originally in N_NLS_DampedNewton, but it makes
//                 more sense to have it here.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------

bool N_NLS_NLParams::setOptions (const N_UTL_OptionBlock & OB)
{
  for (std::list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
       it_tpL != OB.getParams().end(); ++ it_tpL)
  {
    if (it_tpL->uTag() == "ABSTOL")
    {
      setAbsTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "RELTOL")
    {
      setRelTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "DELTAXTOL")
    {
      setDeltaXTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "SMALLUPDATETOL")
    {
      setSmallUpdateTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "ENFORCEDEVICECONV")
    {
      setEnforceDeviceConvFlag(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "RHSTOL")
    {
      setRHSTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "MAXSTEP")
    {
      setMaxNewtonStep(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "LINOPT")
    {
      setLinearOpt(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "CONSTRAINTBT")
    {
      setConstraintBT(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "CONSTRAINTMAX")
    {
      setGlobalBTMax(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "CONSTRAINTMIN")
    {
      setGlobalBTMin(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "CONSTRAINTCHANGE")
    {
      setGlobalBTChange(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "NLSTRATEGY")
    {
      setNLStrategy(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "SEARCHMETHOD")
    {
      setSearchMethod(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "MAXSEARCHSTEP")
    {
      setMaxSearchStep(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "IN_FORCING")
    {
      setForcingFlag(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "NORMLVL")
    {
      setNormLevel(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "NOX")
    {
      // do nothing.
    }
    else if (it_tpL->uTag() == "MATRIXMARKET")
    {
#ifdef Xyce_DEBUG_NONLINEAR
      setMMFormat (static_cast<bool>(it_tpL->getImmutableValue<double>()));
#endif
    }
    else if (it_tpL->uTag() == "DEBUGLEVEL")
    {
#ifdef Xyce_DEBUG_NONLINEAR
      setDebugLevel(it_tpL->getImmutableValue<int>());
#endif
    }
    else if (it_tpL->uTag() == "DEBUGMINTIMESTEP")
    {
#ifdef Xyce_DEBUG_NONLINEAR
      setDebugMinTimeStep(it_tpL->getImmutableValue<int>());
#endif
    }
    else if (it_tpL->uTag() == "DEBUGMAXTIMESTEP")
    {
#ifdef Xyce_DEBUG_NONLINEAR
      setDebugMaxTimeStep(it_tpL->getImmutableValue<int>());
#endif
    }
    else if (it_tpL->uTag() == "DEBUGMINTIME")
    {
#ifdef Xyce_DEBUG_NONLINEAR
      setDebugMinTime(it_tpL->getImmutableValue<double>());
#endif
    }
    else if (it_tpL->uTag() == "DEBUGMAXTIME")
    {
#ifdef Xyce_DEBUG_NONLINEAR
      setDebugMaxTime(it_tpL->getImmutableValue<double>());
#endif
    }
    else if (it_tpL->uTag() == "SCREENOUTPUT")
    {
#ifdef Xyce_DEBUG_NONLINEAR
      setScreenOutputFlag (static_cast<bool>(it_tpL->getImmutableValue<double>()));
#endif
    }
    else
    {
      std::string tmp = it_tpL->uTag() +
        " is not a recognized nonlinear solver option.\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, tmp);
    }
  }

  setCmdLineOptions ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::printParams
// Purpose       : Print out the nonlinear solver parameter values.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
void N_NLS_NLParams::printParams(std::ostream &os)
{
  os << "\n" << std::endl
               << Xyce::section_divider << std::endl;
  os << "\n***** Nonlinear solver options:\n" << std::endl
               << "\tabsTol:\t\t\t" << getAbsTol()
               << "\trelTol:\t\t\t" << getRelTol()
               << "\tdeltaXTol (weighted):\t" << getDeltaXTol()
               << "\tRHSTol:\t\t\t" << getRHSTol()
               << "\tSmall Update Tol:\t" << getSmallUpdateTol()
               << "\tmax NL Steps:\t\t" << getMaxNewtonStep();

  if (analysisMode_ == DC_OP)
    os << "\tAnalysis Mode:\t\t" << analysisMode_ << "\t(DC Op)" << std::endl;
  else if (analysisMode_ == DC_SWEEP)
    os << "\tAnalysis Mode:\t\t" << analysisMode_ <<  "\t(DC Sweep)" << std::endl;
  else if (analysisMode_ == TRANSIENT)
    os << "\tAnalysis Mode:\t\t" << analysisMode_ << "\t(Transient)" << std::endl;

  NLStrategy strategy = getNLStrategy();
  if (strategy == NEWTON)
    os << "\tNL Strategy:\t\t" << strategy << "\t(None => Newton)" << std::endl;
  else if (strategy == GRADIENT)
    os << "\tNL Strategy:\t\t" << strategy << "\t(Gradient)" << std::endl;
  else if (strategy == NEWTON_GRADIENT)
    os << "\tNL Strategy:\t\t" << strategy << "\t(Newton/Gradient)" << std::endl;
  else if (strategy == MOD_NEWTON)
    os << "\tNL Strategy:\t\t" << strategy << "\t(Modified-Newton)" << std::endl;
  else if (strategy == MOD_NEWTON_GRADIENT)
    os << "\tNL Strategy:\t\t" << strategy << "\t(Modified-Newton/Gradient)" << std::endl;

  LineSearchMethod searchMethod = getSearchMethod();
  if (searchMethod == FULL)
    os << "\tsearch method:\t\t" << searchMethod << "\t(None => Full Newton Steps)" << std::endl;

  else if (searchMethod == DIVIDE)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Divide)" << std::endl;

  else if (searchMethod == BACKTRACK)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Backtrack)" << std::endl;

  else if (searchMethod == SIMPLE_BACKTRACK)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Simple Backtrack)" << std::endl;

  else if (searchMethod == BANK_ROSE)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Bank and Rose Algorithm)" << std::endl;

  else if (searchMethod == DESCENT)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Line Search)" << std::endl;

  os << "\tmax search steps:\t" << getMaxSearchStep()
               << "\tinexact-Newton forcing:\t" << getForcingFlag()
               << "\tnorm level:\t\t" << getNormLevel()
               << "\tlinear optimization:\t" << getLinearOpt()
               << "\tconstraint backtrack:\t" << getConstraintBT()
#ifdef Xyce_DEBUG_NONLINEAR
               << "\tdebugLevel:\t\t" << getDebugLevel ()
               << "\tdebugMinTimeStep:\t" << getDebugMinTimeStep ()
               << "\tdebugMaxTimeStep:\t" << getDebugMaxTimeStep ()
#endif
               << Xyce::section_divider << "\n" << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::operator=
// Purpose       : "=" operator.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 09/05/01
//-----------------------------------------------------------------------------
N_NLS_NLParams & N_NLS_NLParams::operator=(const N_NLS_NLParams & right)
{
  commandLine_   = right.commandLine_;
  nlStrategy_    = right.nlStrategy_;
  searchMethod_  = right.searchMethod_;
  direction_     = right.direction_;
  deltaXTol_     = right.deltaXTol_;
  RHSTol_        = right.RHSTol_;
  absTol_        = right.absTol_;
  relTol_        = right.relTol_;
  maxNewtonStep_ = right.maxNewtonStep_;
  maxSearchStep_ = right.maxSearchStep_;
  INForcingFlag_ = right.INForcingFlag_;
  eta_           = right.eta_;
  normLevel_     = right.normLevel_;

  analysisMode_  = right.analysisMode_;

  linearOptimization_ = right.linearOptimization_;

  constraintBT_   = right.constraintBT_;
  globalBTMax_    = right.globalBTMax_;
  globalBTMin_    = right.globalBTMin_;
  globalBTChange_ = right.globalBTChange_;

#ifdef Xyce_DEBUG_NONLINEAR
  // Debug output options:
  debugLevel_       = right.debugLevel_;
  debugMinTimeStep_ = right.debugMinTimeStep_;
  debugMaxTimeStep_ = right.debugMaxTimeStep_;
  debugMinTime_     = right.debugMinTime_;
  debugMaxTime_     = right.debugMaxTime_;
#endif

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setCmdLineOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 05/24/05
//-----------------------------------------------------------------------------

bool N_NLS_NLParams::setCmdLineOptions ()
{
#ifdef Xyce_DEBUG_NONLINEAR
  // set (or override) debug levels based on command line options
  if ( commandLine_->argExists( "-ndl" ) )
  {
    setDebugLevel( atoi( commandLine_->getArgumentValue( "-ndl" ).c_str() ) );
  }
#endif
  return true;
}


