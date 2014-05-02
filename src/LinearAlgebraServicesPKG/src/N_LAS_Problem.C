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
// Filename       : $RCSfile: N_LAS_Problem.C,v $
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>

// ----------   Xyce Includes   ----------

#include <N_LAS_Problem.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>

//-----------------------------------------------------------------------------
// Function      : N_LAS_Problem::N_LAS_Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
N_LAS_Problem::N_LAS_Problem( const RCP<N_LAS_Matrix> & A, const RCP<N_LAS_MultiVector> & x, const RCP<N_LAS_MultiVector> & b )
 : A_(A),
   x_(x),
   b_(b),
   epetraProblem_( rcp( new Epetra_LinearProblem( dynamic_cast<Epetra_RowMatrix*>(&(A_->epetraObj())),
                                             &(x_->epetraObj()),
                                             &(b_->epetraObj()) ) ) )
{
  matrixFreeFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Problem::N_LAS_Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/04/08
//-----------------------------------------------------------------------------
N_LAS_Problem::N_LAS_Problem( const RCP<Epetra_Operator> & Op, const RCP<N_LAS_MultiVector> & x, const RCP<N_LAS_MultiVector> & b )
 : Op_(Op),
   x_(x),
   b_(b),
   epetraProblem_( rcp( new Epetra_LinearProblem( &*Op,
                                             &(x_->epetraObj()),
                                             &(b_->epetraObj()) ) ) )
{
  matrixFreeFlag_ = true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Problem::N_LAS_Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 10/08/08
//-----------------------------------------------------------------------------
N_LAS_Problem::N_LAS_Problem( const RCP<Epetra_LinearProblem> & epetraProblem )
 : x_(Teuchos::rcp(new N_LAS_MultiVector(epetraProblem->GetLHS(), false))),
   b_(Teuchos::rcp(new N_LAS_MultiVector(epetraProblem->GetRHS(), false))),
   epetraProblem_(epetraProblem)
{
  if (epetraProblem_->GetMatrix()) {
    matrixFreeFlag_ = false;
    A_ = Teuchos::rcp(new N_LAS_Matrix(dynamic_cast<Epetra_CrsMatrix *>(epetraProblem_->GetMatrix()), false));
  }
  else {
    matrixFreeFlag_ = true;
    Op_ = Teuchos::rcp(epetraProblem_->GetOperator(), false);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Problem::~N_LAS_Problem
// Purpose       : 
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
N_LAS_Problem::~N_LAS_Problem()
{
}

