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
// Filename       : $RCSfile: N_LAS_HBBlockJacobiEpetraOperator.C,v $
//
// Purpose        :
//
// Creator        : Heidi Thornquist, 1437
//
// Creation Date  : 09/04/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_HBBlockJacobiEpetraOperator.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_PDS_ParMap.h>

// ----------   Trilinos Includes   ----------

#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Amesos_BaseSolver.h>

using Teuchos::RCP;

//-----------------------------------------------------------------------------
// Function      : matrixFreeEpetraOperator
// Purpose       : non-member constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
RCP<N_LAS_HBBlockJacobiEpetraOperator> blockJacobiOperator(
    const std::vector<Teuchos::RCP<Epetra_LinearProblem> >& epetraProblems,
    const std::vector<Teuchos::RCP<Amesos_BaseSolver> >& amesosSolvers,
    const Teuchos::RCP<N_LAS_HBBuilder>& hbBuilder
    )
{
  RCP<N_LAS_HBBlockJacobiEpetraOperator> epetraOperator =
    rcp(new N_LAS_HBBlockJacobiEpetraOperator);
  epetraOperator->initialize(epetraProblems,
      amesosSolvers,
      hbBuilder
      );
  return epetraOperator;
}


//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::N_LAS_HBBlockJacobiEpetraOperator
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
N_LAS_HBBlockJacobiEpetraOperator::N_LAS_HBBlockJacobiEpetraOperator()
{
  isInitialized_ = false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::N_LAS_HBBlockJacobiEpetraOperator
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
N_LAS_HBBlockJacobiEpetraOperator::~N_LAS_HBBlockJacobiEpetraOperator()
{
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::initialize
// Purpose       : Initialization
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
void N_LAS_HBBlockJacobiEpetraOperator::initialize(
      const std::vector<Teuchos::RCP<Epetra_LinearProblem> >& epetraProblems,
      const std::vector<Teuchos::RCP<Amesos_BaseSolver> >& amesosSolvers,
      const Teuchos::RCP<N_LAS_HBBuilder>& hbBuilder
    )
{
  epetraProblems_ = epetraProblems;
  amesosSolvers_ = amesosSolvers;
  hbBuilder_ = hbBuilder;
  N_ = epetraProblems_.size();
  M_ = (N_-1)/2;
  isInitialized_ = true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::SetUseTranspose
// Purpose       : Define if transpose Apply and ApplyInverse is to be used.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBBlockJacobiEpetraOperator::SetUseTranspose(bool UseTranspose)
{
  // This is not supported for the HB load layers.
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::ApplyInverse
// Purpose       : Apply matrix free preconditioner with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBBlockJacobiEpetraOperator::ApplyInverse(
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
  int status = ApplyInverse(las_X,las_Y);
  // COPY the Ycopy data back into Y
  Y = las_Y.epetraObj();
  return(status);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::ApplyInverse
// Purpose       : Apply matrix free preconditioner with N_LAS_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBBlockJacobiEpetraOperator::ApplyInverse(
  const N_LAS_MultiVector& X,
  N_LAS_MultiVector& Y
  ) const
{
  if (!isInitialized_)
  {
    std::string msg = "N_LAS_HBBlockJacobiEpetraOperator::ApplyInverse:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  bool status = true;

  // Determine number of time-domain variables.
  int n = X.localLength() / (2*N_);

  // Apply the j-th block to all the vectors in the multivector.
  for (int nB=0; nB<N_; ++nB) {

    Epetra_MultiVector * nB_RHS = epetraProblems_[nB]->GetRHS();
    double * nB_Soln = epetraProblems_[nB]->GetLHS()->Values();

    // Apply preconditioner to each vector
    for (int i=0 ; i<X.numVectors() ; ++i) {

      const N_LAS_Vector x(X.epetraVector(i));
      N_LAS_Vector y(Y.epetraVector(i));

      // Load the real and imaginary part of the RHS.
      // NOTE:  We have to convert from complex pair storage to real-equivalent form.
      for (int j=0; j<n; ++j) {
	nB_RHS->ReplaceMyValue(j, 0, x[j*(2*N_)+2*nB]);        // real
	nB_RHS->ReplaceMyValue(n+j, 0, x[j*(2*N_)+2*nB+1]);    // imaginary
      }

      // Solve the linear system
      amesosSolvers_[nB]->Solve();

      // Copy the solution back into y
      for (int j=0; j<n; ++j) {
	y[j*(2*N_)+2*nB]   = nB_Soln[j];
	y[j*(2*N_)+2*nB+1] = nB_Soln[n+j];
      }
 
    }
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
// Function      : N_LAS_HBBlockJacobiEpetraOperator::Apply
// Purpose       : Apply inverse of matrix free preconditioner with Epetra_MultiVectors
// Special Notes : Not supported!
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBBlockJacobiEpetraOperator::Apply(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  std::string msg = "N_LAS_HBBlockJacobiEpetraOperator::Apply is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::Apply
// Purpose       : Apply inverse of matrix free preconditioner with N_LAS_MultiVectors
// Special Notes : Not supported!
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBBlockJacobiEpetraOperator::Apply(
  const N_LAS_MultiVector& X,
  N_LAS_MultiVector& Y
  ) const
{
  std::string msg = "N_LAS_HBBlockJacobiEpetraOperator::Apply is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::NormInf
// Purpose       : Norm Inf of matrix
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
double N_LAS_HBBlockJacobiEpetraOperator::NormInf() const
{
  std::string msg = "N_LAS_HBBlockJacobiEpetraOperator::NormInf is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1.0;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::Label
// Purpose       : Label for operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const char * N_LAS_HBBlockJacobiEpetraOperator::Label() const
{
  return "Matrix Free Harmonic Balance Block Jacobi Preconditioner";
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::UseTranspose
// Purpose       : Query for useTranspose setting
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiEpetraOperator::UseTranspose() const
{
  // Use Transpose is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::HasNormInf
// Purpose       : Query for normInf support
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
bool N_LAS_HBBlockJacobiEpetraOperator::HasNormInf() const
{
  // Norm Inf is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::Comm
// Purpose       : Return Epetra_Comm object
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Comm & N_LAS_HBBlockJacobiEpetraOperator::Comm() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_LAS_HBBlockJacobiEpetraOperator::Comm:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return((epetraProblems_[0])->GetMatrix()->Comm());
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::OperatorDomainMap
// Purpose       : Return Epetra_Map corresponding to domain of operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Map & N_LAS_HBBlockJacobiEpetraOperator::OperatorDomainMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_LAS_HBBlockJacobiEpetraOperator::OperatorDomainMap:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return(*(hbBuilder_->getSolutionMap()));
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBlockJacobiEpetraOperator::OperatorRangeMap
// Purpose       : Return Epetra_Map corresponding to range of operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Map & N_LAS_HBBlockJacobiEpetraOperator::OperatorRangeMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_LAS_HBBlockJacobiEpetraOperator::OperatorRangeMap:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return(*(hbBuilder_->getSolutionMap()));
}

