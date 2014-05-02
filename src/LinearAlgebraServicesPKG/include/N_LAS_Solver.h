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
// Filename       : $RCSfile: N_LAS_Solver.h,v $
//
// Purpose        : Abstract interface to linear solver type.
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/17/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Solver_h
#define Xyce_N_LAS_Solver_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_UTL_fwd.h>

#include <Teuchos_RCP.hpp>

// ----------  Fwd Declares     ----------

class N_LAS_Problem;
class N_LAS_Preconditioner;

//-----------------------------------------------------------------------------
// Class         : N_LAS_Solver
// Purpose       : Abstract interface to linear solver type.
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/17/04
//-----------------------------------------------------------------------------
class N_LAS_Solver
{

public:
  //Constructors
  N_LAS_Solver( bool isIterative )
  : solutionTime_(0.0),
    isIterative_(isIterative)
  {}

  //Destructor
  virtual ~N_LAS_Solver() {}

  // Set the solver options
  virtual bool setOptions(const N_UTL_OptionBlock & OB) = 0;
  virtual bool setDefaultOptions() = 0;
  virtual bool setDefaultOption( const std::string & option ) = 0;
  virtual bool setPreconditioner( const Teuchos::RCP<N_LAS_Preconditioner>& precond ) { return true; }

  // Set individual options
  virtual bool setParam( const N_UTL_Param & param ) = 0;

  // Get info such as Num Iterations, Residual, etc.
  virtual bool getInfo( N_UTL_Param & info ) = 0;

  //Residual
  virtual double residual() { return 0.0; }
  virtual void setTolerance( const double & tol ) {}

  // Solve function: x = A^(-1) b.
  // input parameter 'ReuseFactors': If 'true', do not factor A, rather reuse
  // factors from previous solve.  Useful for inexact nonlinear techniques and
  // multiple RHS solves.
  virtual int solve( bool ReuseFactors = false ) = 0;

  double solutionTime() { return solutionTime_; }

  bool isIterative() { return isIterative_; }

protected:
  double solutionTime_;

private:
  const bool isIterative_;

  //No copying
  N_LAS_Solver(const N_LAS_Solver & right);
  N_LAS_Solver & operator=(const N_LAS_Solver & right);

  //No comparison
  bool operator==(const N_LAS_Solver & right) const;
  bool operator!=(const N_LAS_Solver & right) const;

};

#endif
