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
// Filename       : $RCSfile: N_NLS_NLParams.h,v $
//
// Purpose        : Basic user-specifid parameters data structure.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
// Revision Information:
// ----------------------
//
// Revision Number: $Revision: 1.82 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NLParams_h
#define Xyce_N_NLS_NLParams_h

// ----------   Standard Includes   ----------

#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <assert.h>

#ifdef HAVE_CLIMITS
#include <climits>
#else
#include <limits.h>
#endif

#include <N_UTL_Xyce.h>
#include <N_IO_fwd.h>
#include <N_NLS_Manager.h>
#include <N_ERH_ErrorMgr.h>

// ---------- Enumerated Types ----------

// Support line-search methods.
enum LineSearchMethod {
  FULL,            // Newton's Method, undamped.
  DIVIDE,          // Reduce step-size by successively halving the previous
                   // step.
  BACKTRACK,       // Backtrack using interpolation (ala Dennis & Schnabel).
  BANK_ROSE,       // Bank and Rose algorithm - see method notes.
  DESCENT,         // More'-Thuente line search
  SIMPLE_BACKTRACK // Simple backtracking method
};

// Search directions.
enum Direction {
  NEWTON_DIR,            // Newton's direction
  GRADIENT_DIR,          // Steepest descent direction
  MOD_NEWTON_DIR         // Modified Newton direction
};

// Nonlinear solution "strategies".
enum NLStrategy {
  NEWTON,             // Pure Newton's method
  GRADIENT,           // Pure gradient method
  NEWTON_GRADIENT,    // Combined Newton/gradient method
  MOD_NEWTON,         // Pure modified Newton's method
  MOD_NEWTON_GRADIENT // Combined modified-Newton/gradient method
};

//-----------------------------------------------------------------------------
// Class         : N_NLS_NLParams
// Purpose       : Stores nonlinear solver user settings.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
class N_NLS_NLParams
{
public:
  //N_NLS_NLParams(AnalysisMode mode=DC_OP);
  N_NLS_NLParams(AnalysisMode mode, N_IO_CmdParse & cp);
  N_NLS_NLParams(const N_NLS_NLParams & right);

  ~N_NLS_NLParams();

  N_NLS_NLParams & operator=(const N_NLS_NLParams & right);

  // ***** Accessor functions *****

  void setPrintParamsFlag();
  void clearPrintParamsFlag();
  bool getPrintParamsFlag() const;

  bool setOptions(const N_UTL_OptionBlock & OB);
  bool setCmdLineOptions ();

  inline void setNLStrategy(NLStrategy strategy);
  inline void setNLStrategy(int strategy);
  inline void resetNLStrategy();
  inline NLStrategy getNLStrategy() const;

  inline void setSearchMethod(LineSearchMethod method);
  inline void setSearchMethod(int method);
  inline void resetSearchMethod();
  inline LineSearchMethod getSearchMethod() const;

  inline void setDirection(Direction value);
  inline void resetDirection();
  inline Direction getDirection() const;

  inline void   setDeltaXTol(double Tolerance);
  inline void   resetDeltaXTol();
  inline double getDeltaXTol() const;

  inline void   setSmallUpdateTol (double Tolerance);
  inline void   resetSmallUpdateTol ();
  inline double getSmallUpdateTol () const;

  inline void   setEnforceDeviceConvFlag (bool flag);
  inline void   resetEnforceDeviceConvFlag ();
  inline bool   getEnforceDeviceConvFlag () const;

  inline void   setRHSTol(double Tolerance);
  inline void   resetRHSTol();
  inline double getRHSTol() const;

  inline void   setAbsTol(double Tolerance);
  inline void   resetAbsTol();
  inline double getAbsTol() const;

  inline void   setRelTol(double Tolerance);
  inline void   resetRelTol();
  inline double getRelTol() const;

  inline void     setMaxNewtonStep(unsigned maxNewtonStep);
  inline void     resetMaxNewtonStep();
  inline unsigned getMaxNewtonStep() const;

  inline void     setMaxSearchStep(unsigned maxSearchStep);
  inline void     resetMaxSearchStep();
  inline unsigned getMaxSearchStep() const;

  inline void     setForcingFlag(bool flag);
  inline void     resetForcingFlag();
  inline bool     getForcingFlag() const;

  inline void     setForcingTerm(double value);
  inline void     resetForcingTerm();
  inline double   getForcingTerm() const;

  inline void     setNormLevel(int level);
  inline void     resetNormLevel();
  inline int      getNormLevel() const;

  inline void     setLinearOpt(bool flag);
  inline void     resetLinearOpt();
  inline bool     getLinearOpt() const;

  inline void     setConstraintBT(bool flag);
  inline void     resetConstraintBT();
  inline bool     getConstraintBT() const;

  inline void     setGlobalBTMax(double value);
  inline void     resetGlobalBTMax();
  inline double   getGlobalBTMax() const;

  inline void     setGlobalBTMin(double value);
  inline void     resetGlobalBTMin();
  inline double   getGlobalBTMin() const;

  inline void     setGlobalBTChange(double value);
  inline void     resetGlobalBTChange();
  inline double   getGlobalBTChange() const;

    void printParams(std::ostream &os);

  inline void setDebugLevel(int value);
  inline void resetDebugLevel();
  inline int  getDebugLevel() const;

  inline void setDebugMinTimeStep(int value);
  inline void resetDebugMinTimeStep();
  inline int getDebugMinTimeStep() const;

  inline void setDebugMaxTimeStep(int value);
  inline void resetDebugMaxTimeStep();
  inline int getDebugMaxTimeStep() const;

  inline void setDebugMinTime(double value);
  inline void resetDebugMinTime();
  inline double getDebugMinTime() const;

  inline void setDebugMaxTime(double value);
  inline void resetDebugMaxTime();
  inline double getDebugMaxTime() const;

  inline void setScreenOutputFlag (bool bval);
  inline void resetScreenOutputFlag ();
  inline bool getScreenOutputFlag () const;

  inline void setMMFormat(bool value);
  inline void resetMMFormat();
  inline bool getMMFormat() const;

protected:
private:

public:

protected:
  // Print control flag
  bool printParamsFlag_;

  // command line parser
  N_IO_CmdParse * commandLine_;

  // Calling solution method (e.g., DC Operation Point, Transient, etc.)
  AnalysisMode analysisMode_;
  bool         modeToggled_;

  // Nonlinear solution strategy
  NLStrategy nlStrategy_;

  // Damping method (e.g., Bank and Rose)
  LineSearchMethod searchMethod_;

  // Direction flag - dictates which direction to take in advancing the
  // nonlinear solver.
  Direction direction_;

  // Absolute convergence tolerance for the norm of the residual.
  double absTol_;

  // Relative convergence tolerance.
  double relTol_;

  // Weighted deltaX (update) norm tolerance (Petzold)
  double deltaXTol_;

  // Special deltaX (update) norm tolerance, used for the 
  // "small update" test.   
  double smallUpdateTol_;

  // Weighted RHS (residual) norm tolerance
  double RHSTol_;

  // Maximum number of solution steps to attempt.
  unsigned maxNewtonStep_;

  // maximum number of damping steps:
  unsigned maxSearchStep_;

  // inexact-Newton forcing flag:
  bool INForcingFlag_;

  // check device convergence flag:
  bool enforceDeviceConvFlag_;

  // inexact-Newton forcing term:
  double eta_;

  // norm level:
  int normLevel_;

  // linear optimization flag
  bool linearOptimization_;

  // Constraint backtracking flag
  bool constraintBT_;

  // Backtracking constraint scalar values.  Note that these are used to set
  // vectors of values to these defaults but, in general, the constraints can
  // be individually applied to all solution variables (see
  // N_NLS_ConstraintBT.h).
  double globalBTMax_;
  double globalBTMin_;
  double globalBTChange_;

  // Debug output options:
  int debugLevel_;
  int debugMinTimeStep_;
  int debugMaxTimeStep_;
  double debugMinTime_;
  double debugMaxTime_;
  bool screenOutputFlag_;
  bool matrixMarketFormat_;
};

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setPrintParamsFlag
// Purpose       : Sets the flag which controls the printing of the nonlinear
//                 solver parameters list.  This flag should be set anytime one
//                 of the parameters is modified as the code runs.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setPrintParamsFlag()
{
  printParamsFlag_ = true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::clearPrintParamsFlag
// Purpose       : Clears the flag which controls the printing of the nonlinear
//                 solver parameters list.  This flag is cleared once the
//                 output is printed.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::clearPrintParamsFlag()
{
  printParamsFlag_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getPrintParamsFlag
// Purpose       : Returns the flag value
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
inline bool N_NLS_NLParams::getPrintParamsFlag() const
{
  return printParamsFlag_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setNLStrategy
// Purpose       : Accessor method to set the nonlinear solver strategy.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setNLStrategy(NLStrategy strategy)
{
  nlStrategy_ = strategy;

  // Make sure the method selected is supported, otherwise revert to the
  // default.
  if (nlStrategy_ != NEWTON &&
      nlStrategy_ != GRADIENT &&
      nlStrategy_ != NEWTON_GRADIENT &&
      nlStrategy_ != MOD_NEWTON &&
      nlStrategy_ != MOD_NEWTON_GRADIENT)
    resetNLStrategy();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setNLStrategy
// Purpose       : Accessor method to set the nonlinear solver strategy.
// Special Notes : Takes an integer argument
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setNLStrategy(int strategy)
{
  nlStrategy_ = static_cast<NLStrategy> (strategy);

  // Make sure the method selected is supported, otherwise revert to the
  // default.
  if (nlStrategy_ != NEWTON &&
      nlStrategy_ != GRADIENT &&
      nlStrategy_ != NEWTON_GRADIENT &&
      nlStrategy_ != MOD_NEWTON &&
      nlStrategy_ != MOD_NEWTON_GRADIENT)
    resetNLStrategy();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetNLStrategy
// Purpose       : Accessor method to reset the nonlinear solver strategy to
//                 its default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetNLStrategy()
{
  nlStrategy_ = NEWTON;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getNLStrategy
// Purpose       : Accessor method to return the nonlinear solver strategy.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
inline NLStrategy N_NLS_NLParams::getNLStrategy() const
{
  return nlStrategy_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setSearchMethod
// Purpose       : Accessor method to set the line-search method
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/13/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setSearchMethod(LineSearchMethod method)
{
  searchMethod_ = method;

  // Revert to the default if the method selected is not supported
  if (searchMethod_ != FULL &&
      searchMethod_ != DIVIDE &&
      searchMethod_ != BACKTRACK &&
      searchMethod_ != SIMPLE_BACKTRACK &&
      searchMethod_ != BANK_ROSE &&
      searchMethod_ != DESCENT)
    searchMethod_ = FULL;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setSearchMethod
// Purpose       : Accessor method to set the line-search method
// Special Notes : Takes and integer argument
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/13/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setSearchMethod(int method)
{
  searchMethod_ = static_cast<LineSearchMethod> (method);

  // Revert to the default if the method selected is not supported
  if (searchMethod_ != FULL &&
      searchMethod_ != DIVIDE &&
      searchMethod_ != BACKTRACK &&
      searchMethod_ != SIMPLE_BACKTRACK &&
      searchMethod_ != BANK_ROSE &&
      searchMethod_ != DESCENT)
    searchMethod_ = FULL;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetSearchMethod
// Purpose       : Accessor method to reset the default line-search method
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetSearchMethod()
{
  searchMethod_ = FULL;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getSearchMethod
// Purpose       : Accessor method to return the line-search method
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/13/01
//-----------------------------------------------------------------------------
inline LineSearchMethod N_NLS_NLParams::getSearchMethod() const
{
  return searchMethod_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setDirection
// Purpose       : Accessor method to set the nonlinear direction flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setDirection(Direction value)
{

  // Make sure we have a valid value here.
  assert(value >= 0);
  direction_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetDirection
// Purpose       : Accessor method to reset the default nonlinear direction
//                 flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetDirection()
{
  direction_ = NEWTON_DIR;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getDirection
// Purpose       : Accessor method to return the nonlinear diection flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
inline Direction N_NLS_NLParams::getDirection() const
{
  return direction_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setDeltaXTol
// Purpose       : Accessor method to set the deltaX (update) tolerance.
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/20/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setDeltaXTol(double value)
{
  deltaXTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetDeltaXTol
// Purpose       : Accessor method to reset the deltaX (update) tolerance to
//                 the default value.
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/20/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetDeltaXTol()
{
  deltaXTol_ = 1.0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getDeltaXTol
// Purpose       : Accessor method to return the deltaX (update) tolerance
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/20/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getDeltaXTol() const
{
  return deltaXTol_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setSmallUpdateTol
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/02/03
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setSmallUpdateTol (double value)
{
  smallUpdateTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetSmallUpdateTol
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/02/03
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetSmallUpdateTol ()
{
  smallUpdateTol_ = 1.0e-6;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getSmallUpdateTol
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/02/03
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getSmallUpdateTol () const
{
  return smallUpdateTol_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setRHSTol
// Purpose       : Accessor method to set the RHS (residual) tolerance.
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setRHSTol(double value)
{
  RHSTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetRHSTol
// Purpose       : Accessor method to reset the RHS (residual) tolerance to
//                 the default value.
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetRHSTol()
{
  RHSTol_ = 1.0e-06;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getRHSTol
// Purpose       : Accessor method to return the RHS (residual) tolerance
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getRHSTol() const
{
  return RHSTol_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setAbsTol
// Purpose       : Accessor method to set the absTol tolerance.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/22/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setAbsTol(double value)
{
  absTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetAbsTol
// Purpose       : Accessor method to reset the absTol tolerance to
//                 the default value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/21/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetAbsTol()
{
  absTol_ = 1.0e-12;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getAbsTol
// Purpose       : Accessor method to return the absTol tolerance
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/21/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getAbsTol() const
{
  return absTol_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setRelTol
// Purpose       : Accessor method to set the relTol tolerance.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/22/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setRelTol(double value)
{
  relTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetRelTol
// Purpose       : Accessor method to reset the relTol tolerance to
//                 the default value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/21/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetRelTol()
{
  relTol_ = 1.0e-03;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getRelTol
// Purpose       : Accessor method to return the relTol tolerance
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/21/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getRelTol() const
{
  return relTol_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setMaxNewtonStep
// Purpose       : Accessor method to set the maximum number of Newton steps.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setMaxNewtonStep(unsigned value)
{
  maxNewtonStep_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetMaxNewtonStep
// Purpose       : Accessor method to reset the maximum number of Newton steps
//                 to the default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetMaxNewtonStep()
{
  maxNewtonStep_ = 200;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getMaxNewtonStep
// Purpose       : Accessor method to get the maximum number of Newton steps.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline unsigned N_NLS_NLParams::getMaxNewtonStep() const
{
  return maxNewtonStep_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setMaxSearchStep
// Purpose       : Accessor method to set the maximum number of line-search
//                 steps.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setMaxSearchStep(unsigned value)
{
  maxSearchStep_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetMaxSearchStep
// Purpose       : Accessor method to reset the maximum number of line-search
//                 steps to the default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetMaxSearchStep()
{
  maxSearchStep_ = 9;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getMaxSearchStep
// Purpose       : Accessor method to get the maximum number of line-search
//                 steps.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline unsigned N_NLS_NLParams::getMaxSearchStep() const
{
  return maxSearchStep_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setForcingFlag
// Purpose       : Accessor method to set the inexact-Newton forcing flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/04/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setForcingFlag(bool value)
{

  // Reset the linear solver tolerance if we're switching forcing off
  // NOTE:  I couldn't get this to work because we're instantiating several
  // NLParams objects and the one for which the lasSolverPtr is getting set
  // isn't the one needed.  I've found a work-around for now but we'll need to
  // fix this - SAH, 6/25/01
//   if (INForcingFlag_ && !value)
//     lasSolverPtr_->resetTolerance();

  INForcingFlag_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetForcingFlag
// Purpose       : Accessor method to reset the default inexact-Newton forcing
//                 flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/04/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetForcingFlag()
{
    INForcingFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getForcingFlag
// Purpose       : Accessor method to get the inexact-Newton forcing flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/04/01
//-----------------------------------------------------------------------------
inline bool N_NLS_NLParams::getForcingFlag() const
{
  return INForcingFlag_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setEnforceDeviceConvFlag
// Purpose       : Accessor method to set the inexact-Newton forcing term.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 06/30/05
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setEnforceDeviceConvFlag (bool flag)
{
  enforceDeviceConvFlag_ = flag;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetEnforceDeviceConvFlag
// Purpose       : Accessor method to reset the default inexact-Newton forcing
//                 value.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 06/30/05
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetEnforceDeviceConvFlag ()
{
  enforceDeviceConvFlag_ = true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getEnforceDeviceConvFlag
// Purpose       : Accessor method to get the inexact-Newton forcing term.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 06/30/05
//-----------------------------------------------------------------------------
inline bool N_NLS_NLParams::getEnforceDeviceConvFlag () const
{
  return enforceDeviceConvFlag_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setForcingTerm
// Purpose       : Accessor method to set the inexact-Newton forcing term.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/29/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setForcingTerm(double value)
{
  eta_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetForcingTerm
// Purpose       : Accessor method to reset the default inexact-Newton forcing
//                 value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/29/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetForcingTerm()
{
  eta_ = 1.0e-01;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getForcingTerm
// Purpose       : Accessor method to get the inexact-Newton forcing term.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 06/29/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getForcingTerm() const
{
  return eta_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setNormLevel
// Purpose       : Accessor method to set the lp norm level (p).
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setNormLevel(int value)
{
  normLevel_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetNormLevel
// Purpose       : Accessor method to reset the lp norm level (p) to the
//                 default value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetNormLevel()
{
  normLevel_ = 2;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getNormLevel
// Purpose       : Accessor method to return the lp norm level (p).
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline int N_NLS_NLParams::getNormLevel() const
{
  return normLevel_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setLinearOpt
// Purpose       : Accessor method to set the linear optimization flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setLinearOpt(bool flag)
{
  linearOptimization_ = flag;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetLinearOpt
// Purpose       : Accessor method to reset the default linear optimization
//                 flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetLinearOpt()
{
  linearOptimization_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getLinearOpt
// Purpose       : Accessor method to return the linear optimization flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline bool N_NLS_NLParams::getLinearOpt() const
{
  return linearOptimization_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setConstraintBT
// Purpose       : Accessor method to set the constraint backtracking flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/06/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setConstraintBT(bool flag)
{
  constraintBT_ = flag;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetConstraintBT
// Purpose       : Accessor method to reset the default constraint backtracking
//                 flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/06/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetConstraintBT()
{
  constraintBT_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getConstraintBT
// Purpose       : Accessor method to return the constraint backtracking flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/06/01
//-----------------------------------------------------------------------------
inline bool N_NLS_NLParams::getConstraintBT() const
{
  return constraintBT_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setGlobalBTMax
// Purpose       : Accessor method to set the constraint maximum value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setGlobalBTMax(double value)
{
  globalBTMax_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetGlobalBTMax
// Purpose       : Accessor method to set the constraint maximum default value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetGlobalBTMax()
{
  globalBTMax_ = N_UTL_MachineDependentParams::DoubleMax();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getGlobalBTMax
// Purpose       : Accessor method to return the constraint maximum value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getGlobalBTMax() const
{
  return globalBTMax_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setGlobalBTMin
// Purpose       : Accessor method to set the constraint minimum value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setGlobalBTMin(double value)
{
  globalBTMin_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetGlobalBTMin
// Purpose       : Accessor method to set the constraint minimum default value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetGlobalBTMin()
{
  globalBTMin_ = -N_UTL_MachineDependentParams::DoubleMax();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getGlobalBTMin
// Purpose       : Accessor method to return the constraint minimum value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getGlobalBTMin() const
{
  return globalBTMin_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setGlobalBTChange
// Purpose       : Accessor method to set the constraint percentage change
//                 value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setGlobalBTChange(double value)
{
  globalBTChange_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetGlobalBTChange
// Purpose       : Accessor method to set the constraint percentage change
//                 default value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetGlobalBTChange()
{
  globalBTChange_ = sqrt(N_UTL_MachineDependentParams::DoubleMax());
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getGlobalBTChange
// Purpose       : Accessor method to return the constraint percentage change
//                 value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getGlobalBTChange() const
{
  return globalBTChange_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setDebugLevel(int value)
{
  debugLevel_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetDebugLevel()
{
  debugLevel_ = 1;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline int  N_NLS_NLParams::getDebugLevel() const
{
  return debugLevel_;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setDebugMinTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setDebugMinTimeStep(int value)
{
  debugMinTimeStep_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetDebugMinTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetDebugMinTimeStep()
{
  debugMinTimeStep_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getDebugMinTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline int N_NLS_NLParams::getDebugMinTimeStep() const
{
  return debugMinTimeStep_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setDebugMaxTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setDebugMaxTimeStep(int value)
{
  debugMaxTimeStep_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetDebugMaxTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetDebugMaxTimeStep()
{
  debugMaxTimeStep_ = N_UTL_MachineDependentParams::IntMax(); 
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getDebugMaxTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline int N_NLS_NLParams::getDebugMaxTimeStep() const
{
  return debugMaxTimeStep_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setDebugMinTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setDebugMinTime(double value)
{
  debugMinTime_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetDebugMinTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetDebugMinTime()
{
  debugMinTime_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getDebugMinTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getDebugMinTime() const
{
  return debugMinTime_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setDebugMaxTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setDebugMaxTime(double value)
{
  debugMaxTime_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetDebugMaxTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetDebugMaxTime()
{
  debugMaxTime_ = N_UTL_MachineDependentParams::DoubleMax();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getDebugMaxTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline double N_NLS_NLParams::getDebugMaxTime() const
{
  return debugMaxTime_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 04/25/04
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setScreenOutputFlag (bool bval)
{
  screenOutputFlag_ = bval;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 04/25/04
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetScreenOutputFlag ()
{
  screenOutputFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 04/25/04
//-----------------------------------------------------------------------------
inline bool N_NLS_NLParams::getScreenOutputFlag () const
{
  return screenOutputFlag_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::setMMFormat
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::setMMFormat(bool value)
{
  matrixMarketFormat_=value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::resetMMFormat
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_NLS_NLParams::resetMMFormat()
{
  matrixMarketFormat_=false;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NLParams::getMMFormat
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
inline bool N_NLS_NLParams::getMMFormat() const
{
  return matrixMarketFormat_;
}

#endif // Xyce_N_NLS_NLParams_h
