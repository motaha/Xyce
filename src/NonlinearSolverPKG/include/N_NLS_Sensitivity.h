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
// Filename       : $RCSfile: N_NLS_Sensitivity.h,v $
//
// Purpose        : This file contains the sensitivity class.   It mostly
//                  manages the calculations of direct (and possibly later,
//                  adjoint) sensitivities.
//
// Special Notes  : The main reason that this class is derived from
//                  N_NLS_NonLinearSolver is that this class needs to
//                  do a series of linear solves, using the jacobian
//                  matrix.  This seemed similar enough to the requirements
//                  of a nonlinear solver to have one derived off the other.
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/30/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.50 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_Sensitivity_h
#define Xyce_N_NLS_Sensitivity_h

#include<vector>

#include <N_UTL_fwd.h>
#include <N_NLS_NLParams.h>
#include <N_NLS_NonLinearSolver.h>

enum sensDiffMode
{
  SENS_FWD,
  SENS_REV,
  SENS_CNT,
  NUM_DIFF_MODES
};

//-----------------------------------------------------------------------------
// Class         : N_NLS_Sensitivity
// Purpose       :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------

class N_NLS_Sensitivity : public N_NLS_NonLinearSolver
{
public:
  N_NLS_Sensitivity ( N_NLS_NonLinearSolver & nls_, 
                      N_TOP_Topology & top_,
                      N_IO_CmdParse & cp);

  ~N_NLS_Sensitivity ();

   int solve (N_NLS_NonLinearSolver * nlsTmpPtr=NULL) {return -1;};
   int solve (
       std::vector<double> & objectiveVec,
       std::vector<double> & dOdpVec, 
       std::vector<double> & dOdpAdjVec,
       std::vector<double> & scaled_dOdpVec, 
       std::vector<double> & scaled_dOdpAdjVec);

   int solveDirect  ();
   int solveAdjoint ();

   void stdOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities
       );

   void fileOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities
       );

   void dakOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities
       );

   bool calcSensitivities ();

   bool calcObjFuncDerivs ();

   bool setOptions(const N_UTL_OptionBlock& OB);
   bool setSensitivityOptions(const N_UTL_OptionBlock& OB);
   bool setTranOptions(const N_UTL_OptionBlock& OB);
   bool setHBOptions(const N_UTL_OptionBlock& OB);
   
   // Note, many of the following are here b/c they are purely
   // virtual functions of the nonlinear solver class.
   // They don't have much meaning here.  It may turn out that
   // having this class derive off of the N_NLS_NonLinearSolver
   // class doesn't make much sense.  If so, I'll change it later. ERK

   int getNumIterations() const;
#ifdef Xyce_DEBUG_NONLINEAR
   int getDebugLevel() const;
   bool getScreenOutputFlag() const;
   double getDebugMinTime() const;
   double getDebugMaxTime() const;
   int getDebugMinTimeStep() const;
   int getDebugMaxTimeStep() const;
   bool getMMFormat () const;
#endif
   double getMaxNormF() const;
   int getMaxNormFindex() const;

   int getContinuationStep() const;
   int getParameterNumber() const;
   bool isFirstContinuationParam() const;
   bool isFirstSolveComplete() const;
   void setAnalysisMode(AnalysisMode mode);

protected:
private:

public:
protected:
private:

  bool allocateddXVec_;
  int debugLevel_;
  int solutionSize_;
  bool solveDirectFlag_;
  bool solveAdjointFlag_;
  bool outputScaledFlag_; // include scaled sensitivities in IO 
  bool outputUnscaledFlag_; // include unscaled sensitivities in IO
  int maxParamStringSize_;

  bool stdOutputFlag_;
  bool fileOutputFlag_;
  bool dakotaFileOutputFlag_;
  int numSolves_;

  // expression related stuff:
  int difference;
  bool objFuncGiven_;
  bool objFuncGIDsetup_;
  int            expNumVars_;
  std::vector<std::string> expVarNames_;
  std::vector<int>    expVarGIDs_;
  std::vector<int>    expVarLocal_;
  std::vector<double> expVarVals_;
  std::vector<double> expVarDerivs_;
  double         expVal_;
  std::string objFuncString_;

  double curValue_;   // current value of the variable.
  double objFuncEval_;// value of the evaluated objective function.
  double dOdp_;
  double sqrtEta_;
  bool sqrtEtaGiven_;

  std::vector<N_LAS_Vector*> dfdpPtrVector_;
  std::vector<N_LAS_Vector*> dXdpPtrVector_;

  N_LAS_Vector* dOdXVectorPtr_; // size of solution std::vector.
  std::vector<double> dOdpVec_; // size = number of sensitivity params.
  std::vector<double> dOdpAdjVec_; // size = number of sensitivity params.

  std::vector<double> scaled_dOdpVec_; // size = number of sensitivity params.
  std::vector<double> scaled_dOdpAdjVec_; // size = number of sensitivity params.

  std::vector<double> paramOrigVals_; // size = number of sensitivity params.

  N_LAS_Vector * lambdaVectorPtr_;
  N_LAS_Vector * savedRHSVectorPtr_;
  N_LAS_Vector * savedNewtonVectorPtr_;
  N_LAS_Vector * origFVectorPtr_;
  N_LAS_Vector * pertFVectorPtr_;
  N_LAS_Vector * testFVectorPtr_;

  N_NLS_NonLinearSolver & nls_;

  N_TOP_Topology & top_;

  N_UTL_Expression * expPtr_;

  int numSensParams_;
  std::vector<std::string> paramNameVec_;
};

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getNumIterations
// Purpose       : doesn't do anything, is just a placeholder.
// Special Notes : This one may be needed later, I'm not sure.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
inline int N_NLS_Sensitivity::getNumIterations() const
{
  return 0;
}


//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline int N_NLS_Sensitivity::getContinuationStep() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline int N_NLS_Sensitivity::getParameterNumber() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline bool N_NLS_Sensitivity::isFirstContinuationParam() const
{
  return true;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline bool N_NLS_Sensitivity::isFirstSolveComplete() const
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::setAnalysisMode
// Purpose       : doesn't do anything, is just a placeholder.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
inline void N_NLS_Sensitivity::setAnalysisMode(AnalysisMode mode)
{

}

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getDebugLevel
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
inline int N_NLS_Sensitivity::getDebugLevel() const
{
  return -100;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getScreenOutputFlag
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
inline bool N_NLS_Sensitivity::getScreenOutputFlag () const
{
  return false;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getDebugMinTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double N_NLS_Sensitivity::getDebugMinTime() const
{
  return 0.0;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getDebugMaxTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double N_NLS_Sensitivity::getDebugMaxTime() const
{
  return N_UTL_MachineDependentParams::DoubleMax();
}

//---------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getDebugMinTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int N_NLS_Sensitivity::getDebugMinTimeStep() const
{
  return 0;
}

//---------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getDebugMaxTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int N_NLS_Sensitivity::getDebugMaxTimeStep() const
{
  return N_UTL_MachineDependentParams::IntMax();
}

//---------------------------------------------------------------------------
// Function      : N_NLS_Sensitivity::getMMFormat
//
// Return Type   : bool
//---------------------------------------------------------------------------
inline bool N_NLS_Sensitivity::getMMFormat () const
{
  return false;
}


#endif // debug nonlin

#endif

