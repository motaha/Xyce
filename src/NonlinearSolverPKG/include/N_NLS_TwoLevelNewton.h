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
// Filename       : $RCSfile: N_NLS_TwoLevelNewton.h,v $
//
// Purpose        : Defines N_NLS_TwoLevelNewton class.
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/20/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.72.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_TwoLevelNewton_h
#define Xyce_N_NLS_TwoLevelNewton_h

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>
#include <N_NLS_NLParams.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_TwoLevelEnum.h>

// ---------- Forward Declarations ----------
// Time Integrator
class N_ANP_AnalysisInterface;

// Loader for RHS and Jacobian
class N_LOA_Loader;

// Linear Algebra Support
class N_LAS_Matrix;
class N_LAS_Vector;
class N_LAS_System;

// Options
class N_UTL_Param;

//-----------------------------------------------------------------------------
// Class         : N_NLS_TwoLevelNewton
// Purpose       : This class is the manager for the two level newton
//                 algorithm.  Mostly, it will contain a control loop
//                 which repeatedly calls the "true" nonlinear solver
//                 (either N_NLS_DampedNewton's solve, or NOX, depending
//                 on the options specified).
//
//                 The idea of two level newton is divide the problem
//                 up into sub-problems and solve them separately
//                 for part or all of the solve.
//
//                 The control loop contained in this class repeatedly
//                 loops over all the sub-problems, and also determines
//                 if the total problem has converged.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------

class N_NLS_TwoLevelNewton : public N_NLS_NonLinearSolver
{
public:
  N_NLS_TwoLevelNewton(bool noxFlag, bool noxFlagInner, N_IO_CmdParse & cp);
  ~N_NLS_TwoLevelNewton();

  int getNumIterations() const;
#ifdef Xyce_DEBUG_NONLINEAR
  int getDebugLevel () const;
  bool getScreenOutputFlag () const;
  double getDebugMinTime() const;
  double getDebugMaxTime() const;
  int getDebugMinTimeStep() const;
  int getDebugMaxTimeStep() const;
  bool getMMFormat () const;
#endif
  int getContinuationStep() const;
  int getParameterNumber() const;
  bool isFirstContinuationParam() const;
  bool isFirstSolveComplete() const;

  int solve (N_NLS_NonLinearSolver * nlsTmpPtr=NULL);
  bool setOptions     (const N_UTL_OptionBlock & OB);
  bool setTranOptions (const N_UTL_OptionBlock & OB);
  bool setHBOptions (const N_UTL_OptionBlock & OB);
  bool setTwoLevelOptions     (const N_UTL_OptionBlock & OB);
  bool setTwoLevelTranOptions (const N_UTL_OptionBlock & OB);
  bool setLocaOptions (const N_UTL_OptionBlock& OB);
  bool setTwoLevelLocaOptions (const N_UTL_OptionBlock& OB);
  bool registerLinearSystem(N_LAS_System* ptr);
  bool registerAnalysisInterface (N_ANP_AnalysisInterface * tiaPtr_tmp);
  bool registerLoader(N_LOA_Loader* ptr);
  bool registerOutputMgr (N_IO_OutputMgr * ptr);
  bool initializeAll ();
  TwoLevelNewtonMode  getCouplingMode ();
  void setAnalysisMode (AnalysisMode mode);
  bool setPetraOptions (const N_UTL_OptionBlock & OB);
  int  getContStepNumber ();
  bool enableSensitivity ();

  // solver statistics:
  virtual int getNumResidualLoads();
  virtual int getNumJacobianLoads();
  virtual int getNumLinearSolves();
  virtual int getNumFailedLinearSolves();
  virtual int getNumJacobianFactorizations();
  virtual unsigned int getTotalNumLinearIters();
  virtual double getTotalLinearSolveTime();
  virtual double getTotalResidualLoadTime();
  virtual double getTotalJacobianLoadTime();

  double getMaxNormF() const;
  int getMaxNormFindex() const;
private:
  N_NLS_TwoLevelNewton();

#ifdef Xyce_VERBOSE_NONLINEAR

  void printHeader_();
  void printStepInfo_(int step, int success, TwoLevelNewtonMode solveType);

#endif
  void zeroInnerLoopStatistics_ ();
  void calcInnerLoopStatistics_ ();

  void calcOuterLoopStatistics_ ();

  bool calcCouplingTerms_ ();

  int continuationLoop_ ();
  int locaLoop_ ();

  int algorithm0_();
  int algorithm1_();
  int algorithm2_();
  int algorithm3_();
  int algorithm4_();
  int algorithm5_();

public:

private:
  // Pointer to the outer loop nonlinear solver that is being used.
  N_NLS_NonLinearSolver *nlsOuterPtr_;

  // Pointer to the inner loop nonlinear solver that is being used.
  N_NLS_NonLinearSolver *nlsInnerPtr_;

  // Pointer to the time integrator that is being used.
  N_ANP_AnalysisInterface *tiaPtr_;

  // number of iterations for the "outer" loop, if it exists. (algorithm>0)
  int maxOuterSteps_;

  // number of iterations for the inner "voltage stepper" loop.
  // This parameter is kind of irrelevant, if running with variable step
  // size.
  int maxContSteps_;
  int maxContStepsTran_;
  int contStep_;

  double increaseContScalar_;
  double decreaseContScalar_;

  // this variable determines which variant of two level newton
  // is being used.
  //
  // Here's the key:
  //
  // 0 = full newton , one level (as though two-level weren't being used)
  //      This is probably the best algorithm for transient.
  //
  // 1 = full newton outter loop, device pde newton inner loop.
  //
  // 2 = full newton outter loop, with pde device inner loop and
  //     continuation.  Sort of a 3-level.
  //
  // 3 = ckt only outter loop, with pde device inner loop, which uses
  //     continuation.  This is sort of a 3-level Newton.  This is
  //     the best algorithm for DCOP.
  //
  int twoLevelAlgorithm_;
  int twoLevelAlgorithmTran_;  // same, but for transient.

  // this flag indicates if the current solve iteration is of the
  // outer loop or the inner loop.
  bool outerLoopActiveFlag_;

  AnalysisMode externalAnalysisMode;

  bool setupOuterLoopParamsFlag_;
  bool setupTranParamsFlag_;
  bool noxFlag_;
  bool noxFlagInner_;

  // inner loop statistics vars:
  int numResidualLoads_;
  int numJacobianLoads_;
  int numLinearSolves_;
  int numFailedLinearSolves_;
  int numJacobianFactorizations_;
  unsigned int totalNumLinearIters_;
  double totalLinearSolveTime_;
  double totalResidualLoadTime_;
  double totalJacobianLoadTime_;

  bool numInterfaceNodesSetup_;

  vector<int> numInterfaceNodes_;

  TwoLevelNewtonMode twoLevelCouplingMode_;

  N_LAS_Vector *savedRHSPtr_;
  N_LAS_Vector *savedNextSolPtr_;
  N_LAS_Vector *jdxpVectorPtr_;

  int numSubProblems_;
  int continuationType_;
  bool innerLoopFailFatal_;
  bool totalSolveFailFatal_;
  bool doFullNewtonFinalEnforcement_;

  N_NLS_NonLinearSolver * nlsPassingPtr_;

  bool firstDCOPFlag_;

  bool continuationCalledBefore_;

  // continuation parameters used for algorithm 4.  (not the
  // continuationLoop_ function, that's different)
  vector <string> paramNameList;
  vector <double> paramFinalVal;
  vector <double> paramCurrentVal;

  // These options are saved in case they have to be modified
  // over the course of the solve.
  N_UTL_OptionBlock innerSolverOptions_;
  N_UTL_OptionBlock innerLocaOptions_;
  N_UTL_OptionBlock outerLocaOptions_;

  // voltage limiter tolerance:
  double voltLimTol_;
};

//-----------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::getContStepNumber
// Purpose       : returns the current continuation step number.
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Compuational Sciences
// Creation Date : 10/24/02
//-----------------------------------------------------------------------------
inline int N_NLS_TwoLevelNewton::getContStepNumber ()
{
  return contStep_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::isFirstSolveComplete
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/22/06
//-----------------------------------------------------------------------------
inline bool N_NLS_TwoLevelNewton::isFirstSolveComplete () const
{
  return true;
}

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::getDebugLevel
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
inline int N_NLS_TwoLevelNewton::getDebugLevel () const
{
  return -100;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::getScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
inline bool N_NLS_TwoLevelNewton::getScreenOutputFlag () const
{
  return false;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::getDebugMinTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double N_NLS_TwoLevelNewton::getDebugMinTime() const
{
  return 0.0;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::getDebugMaxTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double N_NLS_TwoLevelNewton::getDebugMaxTime() const
{
  return N_UTL_MachineDependentParams::DoubleMax();
}

//---------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::getDebugMinTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int N_NLS_TwoLevelNewton::getDebugMinTimeStep() const
{
  return 0;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::getDebugMaxTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int N_NLS_TwoLevelNewton::getDebugMaxTimeStep() const
{
  return N_UTL_MachineDependentParams::IntMax();
}

//---------------------------------------------------------------------------
// Function      : N_NLS_TwoLevelNewton::getMMFormat
//
// Return Type   : bool
//---------------------------------------------------------------------------
inline bool N_NLS_TwoLevelNewton::getMMFormat () const
{
  return false;
}

#endif // debug nonlin
#endif
