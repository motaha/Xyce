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
// Filename       : $RCSfile: N_LAS_HBFDJacobianEpetraOperator.C,v $
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
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_HBFDJacobianEpetraOperator.h>
#include <N_LAS_HBBuilder.h>
#include <N_LOA_HBLoader.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_PDS_ParMap.h>

// ----------   Trilinos Includes   ----------

#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
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
RCP<N_LAS_HBFDJacobianEpetraOperator> fdJacobianOperator(
    const std::vector<Teuchos::RCP<Epetra_LinearProblem> >& epetraProblems,
    const std::vector<Teuchos::RCP<Amesos_BaseSolver> >& amesosSolvers,
    const Teuchos::RCP<N_LAS_HBBuilder>& hbBuilder,
    const Teuchos::RCP<N_LOA_HBLoader>& hbLoader,
    const std::vector<double>& timeSteps
    )
{
  RCP<N_LAS_HBFDJacobianEpetraOperator> epetraOperator =
    rcp(new N_LAS_HBFDJacobianEpetraOperator);
  epetraOperator->initialize(epetraProblems,
      amesosSolvers,
      hbBuilder,
      hbLoader,
      timeSteps
      );
  return epetraOperator;
}


//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::N_LAS_HBFDJacobianEpetraOperator
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
N_LAS_HBFDJacobianEpetraOperator::N_LAS_HBFDJacobianEpetraOperator()
{
  isInitialized_ = false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::N_LAS_HBFDJacobianEpetraOperator
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
N_LAS_HBFDJacobianEpetraOperator::~N_LAS_HBFDJacobianEpetraOperator()
{
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::initialize
// Purpose       : Initialization
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
void N_LAS_HBFDJacobianEpetraOperator::initialize(
      const std::vector<Teuchos::RCP<Epetra_LinearProblem> >& epetraProblems,
      const std::vector<Teuchos::RCP<Amesos_BaseSolver> >& amesosSolvers,
      const Teuchos::RCP<N_LAS_HBBuilder>& hbBuilder,
      const Teuchos::RCP<N_LOA_HBLoader>& hbLoader,
      const std::vector<double>& timeSteps
    )
{
  epetraProblems_ = epetraProblems;
  amesosSolvers_ = amesosSolvers;
  hbBuilder_ = hbBuilder;
  hbLoader_ = hbLoader;
  timeSteps_ = timeSteps;
  N_ = epetraProblems_.size();
  isInitialized_ = true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::SetUseTranspose
// Purpose       : Define if transpose Apply and ApplyInverse is to be used.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBFDJacobianEpetraOperator::SetUseTranspose(bool UseTranspose)
{
  // This is not supported for the HB load layers.
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::ApplyInverse
// Purpose       : Apply matrix free preconditioner with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBFDJacobianEpetraOperator::ApplyInverse(
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
// Function      : N_LAS_HBFDJacobianEpetraOperator::ApplyInverse
// Purpose       : Apply matrix free preconditioner with N_LAS_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBFDJacobianEpetraOperator::ApplyInverse(
  const N_LAS_MultiVector& X,
  N_LAS_MultiVector& Y
  ) const
{
  if (!isInitialized_)
  {
    std::string msg = "N_LAS_HBFDJacobianEpetraOperator::ApplyInverse:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  bool status = true;

  // Assume one column in the multivector.
  N_LAS_BlockVector bXf( *(X.getVectorView( 0 )), 2*N_ );

  // Create frequency and time domain work vectors.
  Teuchos::RCP<N_LAS_BlockVector> bXtPtr = hbBuilder_->createTimeDomainBlockVector();
  N_LAS_BlockVector bYf( bXf ), bYtPtr( *bXtPtr );
  N_LAS_Vector tmpXt( bXtPtr->block(0) );
  bYf.putScalar( 0.0 );
  bYtPtr.putScalar( 0.0 );
  tmpXt.putScalar( 0.0 );

  // Permute the input vector from the frequency to time domain, since this
  // is a time domain preconditioner.
  hbLoader_->permutedIFT(bXf, &*bXtPtr);  

  // Now apply preconditioner, which consists of a forward solve since we are using
  // Backward Euler.  The matrix representation of BE is:
  //
  // [ C_1/h+G_1    0        0     ... C_M/h ]
  // [  -C_1/h  C_2/h+G_2    0     ...   0   ]
  // [     0     -C_2/h  C_3/h+G_3 ...   0   ]
  // [ and so on ... ]
  //
  // If the entry in the upper right corner of the matrix is ignored, this is a
  // lower triangular, banded matrix which can be solved with forward solve.

  // Get the stored C matrices from the HB loader.
  std::vector<Teuchos::RCP<N_LAS_Matrix> > vecAppdQdx = hbLoader_->getStoredQdx();
  
  for (int nB=0; nB<N_; ++nB)
  {
    // Get time step.
    double h = timeSteps_[0];
    if (timeSteps_.size() > 1)
      h = timeSteps_[nB];

    // Get the pointer to the RHS and solution vector, set them for row nB.
    Epetra_MultiVector * nB_RHS = epetraProblems_[nB]->GetRHS();
    Epetra_MultiVector * nB_Soln = epetraProblems_[nB]->GetLHS(); 

    // Subtract off the previous solutions
    // Only need to remove one solution since the matrix is banded.
    if (nB-1 >= 0)
    {
      // Compute C_i-1 * x_i-1.
      (vecAppdQdx[nB-1]->epetraObj()).Apply( (bXtPtr->block(nB-1)).epetraObj(),
                                             tmpXt.epetraObj() );
 
      // Compute updated RHS:  y_i - C_i-1/h * x_i-1 
      (bXtPtr->block(nB)).addVec( -1.0/h, tmpXt );
    }

    // Set the RHS of the linear system
    (*nB_RHS) = (bXtPtr->block(nB)).epetraObj();

    // Solve the linear system
    amesosSolvers_[nB]->Solve();

    // Copy the solution back to the appropriate bXtPtr block.
    (bXtPtr->block(nB)).epetraObj() = (*nB_Soln);
  }

  // Permute the output vector from the time domain to the frequency domain
  hbLoader_->permutedFFT(*bXtPtr, &bYf);

  // Set Y = bYf. 
  Y = bYf;

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
// Function      : N_LAS_HBFDJacobianEpetraOperator::Apply
// Purpose       : Apply inverse of matrix free preconditioner with Epetra_MultiVectors
// Special Notes : Not supported!
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBFDJacobianEpetraOperator::Apply(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  std::string msg = "N_LAS_HBFDJacobianEpetraOperator::Apply is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::Apply
// Purpose       : Apply inverse of matrix free preconditioner with N_LAS_MultiVectors
// Special Notes : Not supported!
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int N_LAS_HBFDJacobianEpetraOperator::Apply(
  const N_LAS_MultiVector& X,
  N_LAS_MultiVector& Y
  ) const
{
  std::string msg = "N_LAS_HBFDJacobianEpetraOperator::Apply is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::NormInf
// Purpose       : Norm Inf of matrix
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
double N_LAS_HBFDJacobianEpetraOperator::NormInf() const
{
  std::string msg = "N_LAS_HBFDJacobianEpetraOperator::NormInf is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1.0;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::Label
// Purpose       : Label for operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const char * N_LAS_HBFDJacobianEpetraOperator::Label() const
{
  return "Matrix Free Harmonic Balance Finite Difference Preconditioner";
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::UseTranspose
// Purpose       : Query for useTranspose setting
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianEpetraOperator::UseTranspose() const
{
  // Use Transpose is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::HasNormInf
// Purpose       : Query for normInf support
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianEpetraOperator::HasNormInf() const
{
  // Norm Inf is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::Comm
// Purpose       : Return Epetra_Comm object
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Comm & N_LAS_HBFDJacobianEpetraOperator::Comm() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_LAS_HBFDJacobianEpetraOperator::Comm:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return((epetraProblems_[0])->GetMatrix()->Comm());
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::OperatorDomainMap
// Purpose       : Return Epetra_Map corresponding to domain of operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Map & N_LAS_HBFDJacobianEpetraOperator::OperatorDomainMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_LAS_HBFDJacobianEpetraOperator::OperatorDomainMap:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return(*(hbBuilder_->getSolutionMap()));
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianEpetraOperator::OperatorRangeMap
// Purpose       : Return Epetra_Map corresponding to range of operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Map & N_LAS_HBFDJacobianEpetraOperator::OperatorRangeMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "N_LAS_HBFDJacobianEpetraOperator::OperatorRangeMap:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return(*(hbBuilder_->getSolutionMap()));
}

