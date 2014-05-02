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
// Filename       : $RCSfile: N_LAS_ShyLUSolver.h,v $
//
// Purpose        : ShyLU Interface
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 05/18/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_ShyLUSolver_h
#define Xyce_N_LAS_ShyLUSolver_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <Teuchos_RCP.hpp>

#ifdef Xyce_SHYLU
#include <Ifpack_ShyLU.h>
#endif

#include <N_LAS_Solver.h>

// ---------  Other Includes  -----------

// ---------  Fwd Declares --------------

class Epetra_LinearProblem;
class Epetra_MultiVector;
class Epetra_Operator;

class N_LAS_Problem;
//class N_UTL_Timer;

namespace Teuchos {
  class ParameterList;
}

#include <N_UTL_Param.h>
#include <N_LAS_TransformTool.h>

//-----------------------------------------------------------------------------
// Class         : N_LAS_ShyLUSolver
// Purpose       : ShyLU Interface
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
class N_LAS_ShyLUSolver : public N_LAS_Solver
{

public:

  //Constructors
  N_LAS_ShyLUSolver( N_LAS_Problem & problem, N_UTL_OptionBlock & options );

  // Destructor
  ~N_LAS_ShyLUSolver();

  //Control of solver options
  bool setOptions( const N_UTL_OptionBlock & OB );
  bool setDefaultOptions();
  bool setDefaultOption( const string & option );
  bool setParam( const N_UTL_Param & param );

  //Retrieve of solver information
  bool getInfo( N_UTL_Param & info );

  //Actual Solve Call
  int solve( bool ReuseFactors );

  //Residual Value
  double residual() { return linearResidual_; }

  // Output Flag
  const int & getSymmetry() const
  { return symmetry_; }
  void setSymmetry(const int & value)
  { symmetry_ = value; shyluParams_->set("Symmetry", symmetry_); }
  void resetSymmetry()
  { setSymmetry( symmetryDefault_ ); }

  // Maximum Iterations
  const int & getInnerMaxIter() const
  { return innerMaxIter_; }
  void setInnerMaxIter(const int & value)
  { innerMaxIter_ = value; shyluParams_->set("Inner Solver MaxIters", innerMaxIter_); }
  void resetInnerMaxIter()
  { setInnerMaxIter( innerMaxIterDefault_ ); }

  // Linear Convergence Tolerance
  const double & getInnerTol() const
  { return innerTol_; }
  void setInnerTol(const double & value)
  { innerTol_ = value; shyluParams_->set("Inner Solver Tolerance", innerTol_); }
  void resetInnerTol()
  { setInnerTol( innerTolDefault_ ); }

  // Relative Threshold 
  const double & getRelThresh() const
  { return relThresh_; }
  void setRelThresh(const double & value)
  { relThresh_ = value; shyluParams_->set("Relative Threshold", relThresh_); }
  void resetRelThresh()
  { setRelThresh( relThreshDefault_ ); }

  // Diagonal Factor
  const double & getDiagFactor() const
  { return diagFactor_; }
  void setDiagFactor(const double & value)
  { diagFactor_ = value; shyluParams_->set("Diagonal Factor", diagFactor_); }
  void resetDiagFactor()
  { setDiagFactor( diagFactorDefault_ ); }

  // Outer Solver
  const std::string & getOuterSolver() const
  { return outerSolver_; }
  void setOuterSolver(const std::string& value)
  { outerSolver_ = value; shyluParams_->set("Outer Solver Library", outerSolver_); }
  void resetOuterSolver()
  { setOuterSolver( outerSolverDefault_ ); }

  // Separator Type
  const std::string & getSepType() const
  { return separatorType_; }
  void setSepType(const std::string& value)
  { separatorType_ = value; shyluParams_->set("Separator Type", separatorType_); }
  void resetSepType()
  { setSepType( separatorTypeDefault_ ); }

  // Schur Approximation Method
  const std::string & getSchurApproxType() const
  { return schurApproxType_; }
  void setSchurApproxType(const std::string& value)
  { schurApproxType_ = value; shyluParams_->set("Schur Approximation Method", schurApproxType_); }
  void resetSchurApproxType()
  { setSchurApproxType( schurApproxTypeDefault_ ); }

  // Schur Complement Solver
  const std::string & getSchurSolver() const
  { return schurSolver_; }
  void setSchurSolver(const std::string& value)
  { schurSolver_ = value; shyluParams_->set("Schur Complement Solver", schurSolver_); }
  void resetSchurSolver()
  { setSchurSolver( schurSolverDefault_ ); }

private:

  // Number of iterations
  int numLinearIters_;
  // Symmetry
  int symmetry_;
  // Max iterations
  int innerMaxIter_;
  // Tolerance for inner solve
  double innerTol_;
  // Relative threshold
  double relThresh_;
  // Diagonal factor
  double diagFactor_;
  // ShyLU outer solver
  std::string outerSolver_;
  // ShyLU separator type
  std::string separatorType_;
  // ShyLU Schur complement solver
  std::string schurSolver_ ;
  // ShyLU Schur approximation type
  std::string schurApproxType_ ; 

  // Default parameters for solver.
  const int symmetryDefault_;
  const int innerMaxIterDefault_;
  const double innerTolDefault_;
  const double relThreshDefault_;
  const double diagFactorDefault_;
  const std::string outerSolverDefault_;
  const std::string separatorTypeDefault_;
  const std::string schurSolverDefault_;
  const std::string schurApproxTypeDefault_;

  // Set ShyLU options
  bool setShyLUOption_(const char * paramName, const int val);
  // Set ShyLU params
  bool setShyLUParam_(const char * paramName, const double val);

  // Output linear system every outputLS_ or outputBaseLS_ calls
  int outputLS_;
  int outputBaseLS_;

  // Output solver parameters
  void printParams_() const;

  // Linear Problem
  N_LAS_Problem & lasProblem_;
  Epetra_LinearProblem * problem_;

  // Option block pointer
  Teuchos::RCP<N_UTL_OptionBlock> options_;

  // Flag for change of params
  bool updatedParams_;

  // Linear residual returned from the linear solver
  double linearResidual_;

  // Pointer to the ShyLU linear solver
  Teuchos::RCP<Ifpack_ShyLU> solver_;

  // Timing tool
  Teuchos::RCP<N_UTL_Timer> timer_;

  // Transform support
  Teuchos::RCP<N_LAS_Transform> transform_;
  Epetra_LinearProblem * tProblem_;

  // ShyLU parameter list.
  Teuchos::RCP<Teuchos::ParameterList> shyluParams_;

  bool setShyLUCntl_( const N_UTL_Param & param ) { return true; }
};

#endif



