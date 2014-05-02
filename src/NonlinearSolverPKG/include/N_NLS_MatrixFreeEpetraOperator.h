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
// Filename       : $RCSfile: N_NLS_MatrixFreeEpetraOperator.h,v $
//
// Purpose        : This is an HB specific class that derives off of
// Epetra_Operator and supports the matrix-free load method that we need for
// HB.  It takes a pointer to the NonLinearSolver so that it can correctly do
// the apply function, which calls applyJacobian.
//
// Creator        : Todd Coffey, 1414
//
// Creation Date  : 09/04/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_MatrixFreeEpetraOperator_h
#define Xyce_N_NLS_MatrixFreeEpetraOperator_h

// ---------- Standard Includes ----------
#include <Teuchos_RefCountPtr.hpp>

// ---------- Trilinos Includes ----------

#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

// ----------   Xyce Includes   ----------

#include <N_LAS_MultiVector.h>

// ---------- Forward Declarations ----------
class N_NLS_NonLinearSolver;
class N_LAS_Vector;

// ---------- Using Declarations ------------
using Teuchos::RefCountPtr;
using Teuchos::rcp;

//-----------------------------------------------------------------------------
// Class         : N_NLS_MatrixFreeEpetraOperator
// Purpose       : Matrix Free Epetra Operator concrete class
// Special Notes : 
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------

class N_NLS_MatrixFreeEpetraOperator : virtual public Epetra_Operator
{

public:

  N_NLS_MatrixFreeEpetraOperator();

  virtual ~N_NLS_MatrixFreeEpetraOperator();

  void initialize(
      RefCountPtr<N_NLS_NonLinearSolver> nonlinearSolver,
      RefCountPtr<N_LAS_Vector> solVector, 
      RefCountPtr<N_LAS_Vector> rhsVector,
      RefCountPtr<const Epetra_Map> solutionMap
      );

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param[in]
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
    */
  int SetUseTranspose(bool UseTranspose);

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*
    \param[in]
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param[out]
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int Apply(const N_LAS_MultiVector& X, N_LAS_MultiVector& Y) const;

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*
    \param[in]
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param[out]
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must 
              support the case where X and Y are the same object.
    */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int ApplyInverse(const N_LAS_MultiVector& X, N_LAS_MultiVector& Y) const;

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

       \warning This method must not be called unless HasNormInf() returns true.
    */ 
  double NormInf() const;

  //! Returns a character string describing the operator
  const char * Label() const;

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const;

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const;

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const;

private:
  bool isInitialized_;
  Teuchos::RefCountPtr<N_LAS_Vector> solVectorRCPtr_;
  Teuchos::RefCountPtr<N_LAS_Vector> rhsVectorRCPtr_;
  Teuchos::RefCountPtr<N_NLS_NonLinearSolver> nonlinearSolverRCPtr_;
  Teuchos::RefCountPtr<const Epetra_Map> solutionMap_;
};

// Non-member constructor
RefCountPtr<N_NLS_MatrixFreeEpetraOperator> matrixFreeEpetraOperator(
    RefCountPtr<N_NLS_NonLinearSolver> nonlinearSolver,
    RefCountPtr<N_LAS_Vector> solVector,
    RefCountPtr<N_LAS_Vector> rhsVector,
    RefCountPtr<const Epetra_Map> solutionMap
    );

#endif // Xyce_N_NLS_MatrixFreeEpetraOperator_h

