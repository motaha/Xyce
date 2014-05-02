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
// Filename       : $RCSfile: N_LAS_BelosSolver.h,v $
//
// Purpose        : Belos Interface
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
// Revision Number: $Revision: 1.18 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_BelosSolver_h
#define Xyce_N_LAS_BelosSolver_h

#include <Teuchos_RCP.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosSolverManager.hpp>

#include <N_UTL_fwd.h>
#include <N_LAS_Solver.h>

class Epetra_LinearProblem;
class Epetra_MultiVector;
class Epetra_Operator;

class Ifpack_CrsRiluk;
class Ifpack_IlukGraph;

class N_LAS_Problem;

namespace Teuchos {
  class ParameterList;
}

class N_LAS_Preconditioner;
class N_LAS_Problem;

#include <N_UTL_Param.h>
#include <N_LAS_TransformTool.h>

//-----------------------------------------------------------------------------
// Class         : N_LAS_BelosSolver
// Purpose       : Belos Interface
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
class N_LAS_BelosSolver : public N_LAS_Solver
{

public:

  //Constructors
  N_LAS_BelosSolver( N_LAS_Problem & problem, N_UTL_OptionBlock & options );

  // Destructor
  ~N_LAS_BelosSolver();

  //Control of solver options
  bool setOptions( const N_UTL_OptionBlock & OB );
  bool setDefaultOptions();
  bool setDefaultOption( const std::string & option );
  bool setParam( const N_UTL_Param & param );
  bool setPreconditioner( const Teuchos::RCP<N_LAS_Preconditioner>& precond )
  { precond_ = precond; isPrecSet_ = false; return true; }

  //Retrieve of solver information
  bool getInfo( N_UTL_Param & info );

  //Actual Solve Call
  int solve( bool ReuseFactors );

  //Residual Value
  double residual() { return linearResidual_; }

  // Output Flag
  const int & getOutput() const
  { return output_; }
  void setOutput(const int & value)
  { output_ = value; belosParams_->set("Verbosity", output_); }
  void resetOutput()
  { setOutput( outputDefault_ ); }

  // Maximum Iterations
  const int & getMaxIter() const
  { return maxIter_; }
  void setMaxIter(const int & value)
  { maxIter_ = value; belosParams_->set("Maximum Iterations", maxIter_); }
  void resetMaxIter()
  { setMaxIter( maxIterDefault_ ); }

  // Krylov-subspace size for GMRES
  const int & getKSpace() const
  { return KSpace_; }
  void setKSpace(const int & value)
  { KSpace_ = value; belosParams_->set("Num Blocks", KSpace_ ); }
  void resetKSpace()
  { setKSpace( KSpaceDefault_ ); }

  // Linear Convergence Tolerance
  const double & getTolerance() const
  { return tolerance_; }
  void setTolerance(const double & value)
  { tolerance_ = value; belosParams_->set("Convergence Tolerance", tolerance_ ); }
  void resetTolerance()
  { setTolerance( toleranceDefault_ ); }

  // Recycle space size
  const int & getRecycleSize() const
  { return recycle_; }
  void setRecycleSize(const int & value)
  { recycle_ = value; belosParams_->set("Recycle Size", recycle_); }
  void resetRecycleSize()
  { setRecycleSize( recycleDefault_ ); }

private:

  // Output type value.
  int output_;
  // Max iterations
  int maxIter_;
  // Krylov subspace size for restarted GMRES
  int KSpace_;
  // Tolerance for convergence tests
  double tolerance_;
  // Number of iterations
  int numLinearIters_;
  // Number of recycled vectors
  int recycle_;
  // Belos solver
  std::string belosSolver_;

  // Default parameters for solver.
  const int outputDefault_;
  const int maxIterDefault_;
  const int KSpaceDefault_;
  const double toleranceDefault_;
  const int recycleDefault_;
  const std::string belosSolverDefault_;

  // Set Belos options
  bool setBelosOption_(const char * paramName, const int val);
  // Set Belos params
  bool setBelosParam_(const char * paramName, const double val);

  // Output linear system every outputLS_ or outputBaseLS_ calls
  int outputLS_;
  int outputBaseLS_;

  // Output solver parameters
  void printParams_() const;

  // Linear Problem
  N_LAS_Problem & lasProblem_;
  Epetra_LinearProblem * problem_;
  Teuchos::RCP<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > belosProblem_;

  // Option block pointer
  Teuchos::RCP<N_UTL_OptionBlock> options_;

  // Flag for change of params
  bool updatedParams_;

  // Linear residual returned from the linear solver
  double linearResidual_;

  // Pointer to the Belos linear solver
  Teuchos::RCP<Belos::SolverManager<double,Epetra_MultiVector,Epetra_Operator> > solver_;

  // Timing tool
  Teuchos::RCP<N_UTL_Timer> timer_;

  // Transform support
  Teuchos::RCP<N_LAS_Transform> transform_;
  Epetra_LinearProblem * tProblem_;

  bool isPrecSet_;
  Teuchos::RCP<N_LAS_Preconditioner> precond_;

  // Belos preconditioning object.
  Teuchos::RCP<Epetra_Operator> belosPrecond_;

  // Belos parameter list.
  Teuchos::RCP<Teuchos::ParameterList> belosParams_;

  bool setBelosCntl_( const N_UTL_Param & param ) { return true; }
};

#endif



