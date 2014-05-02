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
// Filename       : $RCSfile: N_NLS_MatrixFreeEpetraOperator.C,v $
//
// Purpose        :
//
// Creator        : Todd Coffey, 1414
//
// Creation Date  : 09/04/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_MatrixFreeEpetraOperator.h>
#include <N_LAS_Vector.h>
#include <N_PDS_ParMap.h>

//-----------------------------------------------------------------------------
// Function      : matrixFreeEpetraOperator
// Purpose       : non-member constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
RefCountPtr<N_NLS_MatrixFreeEpetraOperator> matrixFreeEpetraOperator(
    RefCountPtr<N_NLS_NonLinearSolver> nonlinearSolver,
    RefCountPtr<N_LAS_Vector> solVector,
    RefCountPtr<N_LAS_Vector> rhsVector,
    RefCountPtr<const Epetra_Map> solutionMap
    )
{
  RefCountPtr<N_NLS_MatrixFreeEpetraOperator> epetraOperator =
    rcp(new N_NLS_MatrixFreeEpetraOperator);
  epetraOperator->initialize(nonlinearSolver,
      solVector,
      rhsVector,
      solutionMap
      );
  return epetraOperator;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::N_NLS_MatrixFreeEpetraOperator
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
N_NLS_MatrixFreeEpetraOperator::N_NLS_MatrixFreeEpetraOperator()
{
  isInitialized_ = false;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::N_NLS_MatrixFreeEpetraOperator
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
N_NLS_MatrixFreeEpetraOperator::~N_NLS_MatrixFreeEpetraOperator()
{
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::initialize
// Purpose       : Initialization
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
void N_NLS_MatrixFreeEpetraOperator::initialize(
      RefCountPtr<N_NLS_NonLinearSolver> nonlinearSolver,
      RefCountPtr<N_LAS_Vector> solVector,
      RefCountPtr<N_LAS_Vector> rhsVector,
      RefCountPtr<const Epetra_Map> solutionMap
    )
{
  nonlinearSolverRCPtr_ = nonlinearSolver;
  solVectorRCPtr_ = solVector;
  rhsVectorRCPtr_ = rhsVector;
  solutionMap_ = solutionMap;
  isInitialized_ = true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::SetUseTranspose
// Purpose       : Define if transpose Apply and ApplyInverse is to be used.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int N_NLS_MatrixFreeEpetraOperator::SetUseTranspose(bool UseTranspose)
{
  // This is not supported for the HB load layers.
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::Apply
// Purpose       : Apply matrix free operator with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int N_NLS_MatrixFreeEpetraOperator::Apply(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  // Convert these to N_LAS_MultiVectors and call the other Apply

  // Cast away the const until the Apply which will enforce it.
  // This is necessary because there is no const view function in N_LAS_MultiVector
//  Epetra_MultiVector* Xptr = const_cast<Epetra_MultiVector*>(&X);
//  N_LAS_MultiVector las_X(Xptr);  // This is the wrong thing to do, when it goes out of scope, it deletes the Epetra_MultiVector Ptr.
//  N_LAS_MultiVector las_Y(&Y);

  // COPY the multi-vector data into new objects on the stack.
  Epetra_MultiVector* Xcopy = new Epetra_MultiVector(X); // This gets deleted by the N_LAS_MultiVector below
  Epetra_MultiVector* Ycopy = new Epetra_MultiVector(Y); // This gets deleted by the N_LAS_MultiVector below
  N_LAS_MultiVector las_X(Xcopy); // this co-ops the Epetra_MultiVector and uses (and owns) its memory
  N_LAS_MultiVector las_Y(Ycopy); // this co-ops the Epetra_MultiVector and uses (and owns) its memory
  int status = Apply(las_X,las_Y);
  // COPY the Ycopy data back into Y
  Y = las_Y.epetraObj();
  return(status);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::Apply
// Purpose       : Apply matrix free operator with N_LAS_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int N_NLS_MatrixFreeEpetraOperator::Apply(
  const N_LAS_MultiVector& X,
  N_LAS_MultiVector& Y
  ) const
{
  if (!isInitialized_)
  {
    std::string msg = "N_NLS_MatrixFreeEpetraOperator::Apply:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  bool status = true;
  for (int i=0 ; i<X.numVectors() ; ++i)
  {
    const N_LAS_Vector x(X.epetraVector(i));
    N_LAS_Vector y(Y.epetraVector(i));
    bool localStatus = nonlinearSolverRCPtr_->applyJacobian(x,y);
    status = status && localStatus;
  }
  if (status)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}
//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::ApplyInverse
// Purpose       : Apply inverse of matrix free operator with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int N_NLS_MatrixFreeEpetraOperator::ApplyInverse(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  std::string msg = "N_NLS_MatrixFreeEpetraOperator::ApplyInverse is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::ApplyInverse
// Purpose       : Apply inverse of matrix free operator with N_LAS_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int N_NLS_MatrixFreeEpetraOperator::ApplyInverse(
  const N_LAS_MultiVector& X,
  N_LAS_MultiVector& Y
  ) const
{
  std::string msg = "N_NLS_MatrixFreeEpetraOperator::ApplyInverse is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::NormInf
// Purpose       : Norm Inf of matrix
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
double N_NLS_MatrixFreeEpetraOperator::NormInf() const
{
  std::string msg = "N_NLS_MatrixFreeEpetraOperator::NormInf is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1.0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::Label
// Purpose       : Label for operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const char * N_NLS_MatrixFreeEpetraOperator::Label() const
{
  return "Matrix Free Harmonic Balance Epetra Operator";
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::UseTranspose
// Purpose       : Query for useTranspose setting
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
bool N_NLS_MatrixFreeEpetraOperator::UseTranspose() const
{
  // Use Transpose is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::HasNormInf
// Purpose       : Query for normInf support
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
bool N_NLS_MatrixFreeEpetraOperator::HasNormInf() const
{
  // Norm Inf is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::Comm
// Purpose       : Return Epetra_Comm object
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Comm & N_NLS_MatrixFreeEpetraOperator::Comm() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_NLS_MatrixFreeEpetraOperator::Comm:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return(rhsVectorRCPtr_->epetraObj().Comm());
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::OperatorDomainMap
// Purpose       : Return Epetra_Map corresponding to domain of operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Map & N_NLS_MatrixFreeEpetraOperator::OperatorDomainMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_NLS_MatrixFreeEpetraOperator::OperatorDomainMap:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return(*solutionMap_);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_MatrixFreeEpetraOperator::OperatorRangeMap
// Purpose       : Return Epetra_Map corresponding to range of operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Map & N_NLS_MatrixFreeEpetraOperator::OperatorRangeMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_NLS_MatrixFreeEpetraOperator::OperatorRangeMap:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  const Epetra_Map* emap = dynamic_cast<const Epetra_Map*>(&rhsVectorRCPtr_->epetraObj().Map());
  return(*solutionMap_);
}

