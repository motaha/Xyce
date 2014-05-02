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
// Filename       : $RCSfile: N_LAS_AmesosSolver.h,v $
//
// Purpose        : Amesos Direct Linear Solver Interface
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/24/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.20 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_AmesosSolver_h
#define Xyce_N_LAS_AmesosSolver_h

#include <string>

#include <N_UTL_fwd.h>

#include <N_LAS_Solver.h>
#include <N_LAS_TransformTool.h>
#include <Teuchos_RCP.hpp>

class N_LAS_Problem;

class Amesos_BaseSolver;
class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_Export;

//-----------------------------------------------------------------------------
// Class         : N_LAS_AmesosSolver
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class N_LAS_AmesosSolver : public N_LAS_Solver
{

public:
  // Constructor
  N_LAS_AmesosSolver( const std::string & type,
                      N_LAS_Problem & prob,
                      N_UTL_OptionBlock & options );

  // Destructor
  ~N_LAS_AmesosSolver();

  // Set the solver options
  bool setOptions(const N_UTL_OptionBlock & OB);
  bool setDefaultOptions();
  bool setDefaultOption( const std::string & option );

  // Set individual options
  bool setParam( const N_UTL_Param & param );

  // Get info such as Num Iterations, Residual, etc.
  bool getInfo( N_UTL_Param & info );

  // Solve function: x = A^(-1) b.
  // input parameter 'ReuseFactors': If 'true', do not factor A, rather reuse
  // factors from previous solve.  Useful for inexact nonlinear techniques and
  // multiple RHS solves.
  int solve( bool ReuseFactors = false );

private:

  //Solver Type
  const std::string type_;

  //Primary problem access
  N_LAS_Problem & lasProblem_;
  Epetra_LinearProblem & problem_;

  //Wrapped solver object
  Amesos_BaseSolver * solver_;

  //Repivot every time or use static pivoting
  bool repivot_;
  
  //Have Amesos reindex the linear problem
  bool reindex_;

  //Output linear system every outputLS_ calls
  int outputLS_;
  int outputBaseLS_;
  int outputFailedLS_;

  // Transform Support
  Teuchos::RCP<N_LAS_Transform> transform_;
  Epetra_LinearProblem * tProblem_;

  // Optimized Matrix
  Epetra_LinearProblem * optProb_;
  Epetra_CrsMatrix * optMat_;
  Epetra_CrsMatrix * origMat_;
  Epetra_Export * optExporter_;

  //Options
  N_UTL_OptionBlock * options_;

  //Timer
  N_UTL_Timer * timer_;

};

#endif

