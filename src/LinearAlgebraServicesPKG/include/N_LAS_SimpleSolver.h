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
// Filename       : $RCSfile: N_LAS_SimpleSolver.h,v $
//
// Purpose        : Simple Direct Linear Solver Interface
//
// Special Notes  : This direct solver is used for trivial 1x1 linear systems
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 03/07/2013
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

#ifndef Xyce_N_LAS_SimpleSolver_h
#define Xyce_N_LAS_SimpleSolver_h

#include <string>

#include <N_LAS_Solver.h>
#include <N_UTL_fwd.h>
#include <Teuchos_RCP.hpp>

class N_LAS_Problem;
class N_LAS_Transform;

class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_Export;

//-----------------------------------------------------------------------------
// Class         : N_LAS_SimpleSolver
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class N_LAS_SimpleSolver : public N_LAS_Solver
{

public:
  // Constructor
  N_LAS_SimpleSolver( N_LAS_Problem & prob,
                      N_UTL_OptionBlock & options );

  // Destructor
  ~N_LAS_SimpleSolver();

  // Set the solver options
  bool setOptions(const N_UTL_OptionBlock & OB);
  bool setDefaultOptions();
  bool setDefaultOption( const std::string & option );

  // Set individual options
  bool setParam( const N_UTL_Param & param );

  // Get info such as Num Iterations, Residual, etc.
  bool getInfo( N_UTL_Param & info );

  // Solve function: x = A^(-1) b.
  // This class is only used when A is a 1x1 matrix so x = b / A(1,1)
  int solve( bool ReuseFactors = false );

private:

  //Primary problem access
  N_LAS_Problem & lasProblem_;
  Epetra_LinearProblem & problem_;

  //Output linear system every outputLS_ calls
  int outputLS_;
  int outputBaseLS_;
  int outputFailedLS_;

  //Options
  N_UTL_OptionBlock * options_;

  //Timer
  N_UTL_Timer * timer_;

};

#endif

