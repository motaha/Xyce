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
// Filename       : $RCSfile: N_LAS_KSparseSolver.h,v $
//
// Purpose        : KSparse Direct Linear Solver Interface
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
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_KSparseSolver_h
#define Xyce_N_LAS_KSparseSolver_h

// ---------- Standard Includes ----------

#include <string>
#include <N_UTL_fwd.h>

class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_Export;
class Epetra_CrsKundertSparse;
class Epetra_Import;
class Epetra_Map;
class Epetra_MultiVector;

// ----------   Xyce Includes   ----------

#include <N_LAS_Solver.h>
#include <Teuchos_RCP.hpp>

class N_LAS_Problem;
class N_LAS_Transform;

//-----------------------------------------------------------------------------
// Class         : N_LAS_KSparseSolver
// Purpose       : 
// Special Notes : 
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class N_LAS_KSparseSolver : public N_LAS_Solver
{

public:
  // Constructor
  N_LAS_KSparseSolver( N_LAS_Problem & prob,
                       N_UTL_OptionBlock & options );

  // Destructor
  ~N_LAS_KSparseSolver();

  // Set the solver options
  bool setOptions(const N_UTL_OptionBlock & OB);
  bool setDefaultOptions();
  bool setDefaultOption( const string & option );
   
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

  //Primary problem access
  N_LAS_Problem & lasProblem_;
  Epetra_LinearProblem & problem_;

  //Repivot every time or use static pivoting
  bool repivot_;

  //Output linear system every outputLS_ calls
  int outputLS_;
  int outputBaseLS_;
  int outputFailedLS_;

  // Transform Support
  Teuchos::RCP<N_LAS_Transform> transform_;
  Epetra_LinearProblem * tProblem_;

  // Serialized Matrix (if using parallel load serial solve scenario) 
  Teuchos::RCP<Epetra_Map> serialMap_;
  Teuchos::RCP<Epetra_LinearProblem> serialProblem_;
  Teuchos::RCP<Epetra_CrsMatrix> serialMat_;
  Teuchos::RCP<Epetra_MultiVector> serialLHS_, serialRHS_;
  Teuchos::RCP<Epetra_Import> serialImporter_;

  // Import and export matrices in parallel load serial solve scenario
  Epetra_LinearProblem * importToSerial();
  int exportToGlobal();

  // KSparse Solver
  Teuchos::RCP<Epetra_CrsKundertSparse> solver_;

  //Options
  N_UTL_OptionBlock * options_;

  //Timer
  N_UTL_Timer * timer_;

};

#endif

