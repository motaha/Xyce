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
// Filename       : $RCSfile: N_LAS_Problem.h,v $
//
// Purpose        : interface to linear problem
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
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Problem_h
#define Xyce_N_LAS_Problem_h


#include <Teuchos_RCP.hpp>

#include <N_UTL_Xyce.h>
#include <N_PDS_fwd.h>
#include <N_PDS_ParMap.h>
#include <N_LAS_MultiVector.h>

class N_LAS_Matrix;
class N_LAS_Vector;

class Epetra_LinearProblem;
class Epetra_Operator;

using Teuchos::RCP;
using Teuchos::rcp;

//-----------------------------------------------------------------------------
// Class         : N_LAS_Problem
// Purpose       : interface to linear problem
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/17/04
//-----------------------------------------------------------------------------
class N_LAS_Problem
{

public:
  //Constructors
  N_LAS_Problem( const RCP<N_LAS_Matrix> & A, const RCP<N_LAS_MultiVector> & x, const RCP<N_LAS_MultiVector> & b );
  N_LAS_Problem( const RCP<Epetra_Operator> & Op, const RCP<N_LAS_MultiVector> & x, const RCP<N_LAS_MultiVector> & b );
  N_LAS_Problem( const RCP<Epetra_LinearProblem> & epetraProblem );                                                                                        
  Epetra_LinearProblem & epetraObj() { return *epetraProblem_; }
                                                                                          
  //Destructor
  ~N_LAS_Problem();
                                                                                          
  bool matrixFree() const { return(matrixFreeFlag_); }

private:
                                                                                          
  RCP<N_LAS_Matrix> A_;
  RCP<Epetra_Operator> Op_;
  RCP<N_LAS_MultiVector> x_;
  RCP<N_LAS_MultiVector> b_;
                                                                                          
  RCP<Epetra_LinearProblem> epetraProblem_;
                                                                                          
  bool matrixFreeFlag_;
};
                                                                                          
#endif
