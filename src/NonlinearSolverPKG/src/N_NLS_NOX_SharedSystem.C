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
// Filename       : $RCSfile: N_NLS_NOX_SharedSystem.C,v $
//
// Purpose        : Interface to Xyce vectors for NOX.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.35 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>

#include "N_NLS_NOX_SharedSystem.h"
#include "N_NLS_NOX_Interface.h"
#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "N_LAS_System.h"
#include "N_LAS_Builder.h"
#include "N_ERH_ErrorMgr.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

// ----------   NOX Includes   ----------

// ---------- Namespaces ---------------
using namespace N_NLS_NOX;

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::SharedSystem
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
SharedSystem::SharedSystem(N_LAS_Vector& soln,
			   N_LAS_Vector& f,
			   N_LAS_Matrix& jacobian,
			   N_LAS_Vector& newton,
			   N_LAS_Vector& gradient,
			   N_LAS_System& lasSys,
			   N_NLS_NOX::Interface& interface) :

  xyceSolnPtr_(0),
  xyceFPtr_(0),
  xyceJacobianPtr_(0),
  xyceNewtonPtr_(0),
  xyceGradientPtr_(0),
  xyceLasSysPtr_(0),
  xyceInterfacePtr_(0),
  matrixFreeFlag_(false),
  ownerOfJacobian_(0),
  ownerOfStateVectors_(0),
  ifpackGraphPtr_(0),
  ifpackPreconditionerPtr_(0)
{
  reset(soln, f, jacobian, newton, gradient, lasSys, interface);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::~SharedSystem
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
SharedSystem::~SharedSystem()
{
  deletePreconditioner();
  delete xyceSolnPtr_;
  delete xyceFPtr_;
  delete xyceNewtonPtr_;
  delete xyceGradientPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::reset
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::reset(N_LAS_Vector& x,
			 N_LAS_Vector& f,
			 N_LAS_Matrix& jacobian,
			 N_LAS_Vector& newton,
			 N_LAS_Vector& gradient,
			 N_LAS_System& lasSys,
			 N_NLS_NOX::Interface& interface)
{
  // Clear out old views
  delete xyceSolnPtr_;
  delete xyceFPtr_;
  delete xyceNewtonPtr_;
  delete xyceGradientPtr_;

  xyceJacobianPtr_ = &jacobian;
  xyceLasSysPtr_ = &lasSys;
  xyceInterfacePtr_ = &interface;
  matrixFreeFlag_ = xyceInterfacePtr_->getMatrixFreeFlag();

  // Create views of the data used for fills in xyce
  xyceSolnPtr_ = new N_NLS_NOX::Vector(x, lasSys);
  xyceFPtr_ = new N_NLS_NOX::Vector(f, lasSys);
  xyceNewtonPtr_ = new N_NLS_NOX::Vector(newton, lasSys);
  xyceGradientPtr_ = new N_NLS_NOX::Vector(gradient, lasSys);

  // Wipe the preconditioner clean
  deletePreconditioner();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::computeF
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeF(const Vector& solution, Vector& F,
			    const Group* grp)
{
  ownerOfStateVectors_ = grp;

#ifdef Xyce_NOX_USE_VECTOR_COPY
  *xyceSolnPtr_ = solution;
  bool status = xyceInterfacePtr_->computeF();
#else
  bool status = xyceInterfacePtr_->computeF(F, solution);
#endif

  if (status == false) {
    const string message = "Error: N_NLS_NOX::SharedSystem::computeF() - compute F failed!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

#ifdef Xyce_NOX_USE_VECTOR_COPY
  F = *xyceFPtr_;
#endif
  return status;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::computeJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeJacobian(Group* grp)
{
  ownerOfJacobian_ = grp;

#ifdef Xyce_NOX_USE_VECTOR_COPY
  *xyceSolnPtr_ = grp->getX();
#endif

  if (!areStateVectors(grp)) {
#ifdef Xyce_VERBOSE_NOX
    if (1) { //RPP: Need to add priting utilities to group ctor
      cout << "Warning: N_NLS_NOX::SharedSystem::computeJacobian() - State "
	   << "Vectors are not valid wrt solution!" << endl;
      cout << "Calling computeResidual to fix this!" << endl;
    }
#endif
    // RPP: This line is not needed since we now call the group
    //ownerOfStateVectors_ = grp;
    
    NOX::Abstract::Group::ReturnType status = grp->computeF();

    if (status != NOX::Abstract::Group::Ok) {
      const string message = "N_NLS_NOX::SharedSystem::computeJacobian() - compute F failed!";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
    }

  }

#ifdef Xyce_NOX_USE_VECTOR_COPY
  bool status = xyceInterfacePtr_->computeJacobian();
#else
  bool status = xyceInterfacePtr_->computeJacobian
    (dynamic_cast<const N_NLS_NOX::Vector &> (grp->getX()));
#endif

  if (status == false) {
    const string message = "N_NLS_NOX::SharedSystem::computeJacobian() - compute Jac failed!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::computeNewton
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeNewton(const Vector& F, Vector& Newton,
				 Teuchos::ParameterList& params)
{
  *xyceFPtr_ = F;
  // Zero out the Newton vector
  xyceNewtonPtr_->scale(0.0);
  bool status = xyceInterfacePtr_->computeNewton(params);
  Newton = *xyceNewtonPtr_;

  return status;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::computeGradient
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computeGradient(const Vector& F, Vector& Gradient)
{
  *xyceFPtr_ = F;
  xyceGradientPtr_->scale(0.0);
  bool status = xyceInterfacePtr_->computeGradient();
  Gradient = *xyceGradientPtr_;
  return status;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::applyJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::applyJacobian(const Vector& input, Vector& result) const
{
  if (!matrixFreeFlag_) {
    bool NoTranspose = false;
    xyceJacobianPtr_->matvec(NoTranspose, input.getNativeVectorRef(), result.getNativeVectorRef());
  } else {
    // tscoffe/tmei HB 07/29/08
#ifndef Xyce_NOX_USE_VECTOR_COPY
    const string message = "N_NLS_NOX::SharedSystem::applyJacobian() - ERROR, Xyce_NOX_USE_VECTOR_COPY required";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
#endif
    bool status = xyceInterfacePtr_->applyJacobian(input.getNativeVectorRef(), result.getNativeVectorRef());
    if (status == false) {
      const string message = "N_NLS_NOX::SharedSystem::applyJacobian() - apply Jac failed!";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::applyJacobianTranspose
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  if (!matrixFreeFlag_) {
    bool Transpose = true;
    xyceJacobianPtr_->matvec(Transpose, input.getNativeVectorRef(), result.getNativeVectorRef());
  } else {
      const string message = "N_NLS_NOX::SharedSystem::applyJacobianTranspose() - Not Supported for Matrix Free Loads!";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::computePreconditioner
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::computePreconditioner()
{
  Epetra_CrsMatrix* crs = dynamic_cast<Epetra_CrsMatrix*>
    (&(xyceJacobianPtr_->epetraObj()));

  if (crs == 0) {
    cout << "N_NLS_NOX::SharedSystem::computePreconditioner() - " 
	 << "Dynamic cast to CRS Matrix failed!" << endl;
  }
  
  deletePreconditioner();
  ifpackGraphPtr_ = new Ifpack_IlukGraph(crs->Graph(),
					 1,
					 0);
  ifpackGraphPtr_->ConstructFilledGraph();
  ifpackPreconditionerPtr_ = new Ifpack_CrsRiluk(*ifpackGraphPtr_);
  ifpackPreconditionerPtr_->InitValues(*crs);
  ifpackPreconditionerPtr_->Factor();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::deletePreconditioner
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::deletePreconditioner()
{
  delete ifpackPreconditionerPtr_;
  delete ifpackGraphPtr_;
  ifpackPreconditionerPtr_ = 0;
  ifpackGraphPtr_ = 0;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::applyRightPreconditioning
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool SharedSystem::applyRightPreconditioning(bool useTranspose, 
					     Teuchos::ParameterList& params,
					     const Vector& input, 
					     Vector& result)
{
  if (ifpackPreconditionerPtr_ == 0) {
    const string message = "N_NLS_NOX::SharedSystem::applyRightPreconditioning - Preconditioner is 0!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  if (useTranspose)
    ifpackPreconditionerPtr_->SetUseTranspose(useTranspose);
  
  N_LAS_Vector& nonConstInput = 
    const_cast<N_LAS_Vector&>(input.getNativeVectorRef_());
  Epetra_MultiVector& epVecInput = 
    const_cast<Epetra_MultiVector&>(nonConstInput.epetraObj());
  
  N_LAS_Vector& nonConstResult = 
    const_cast<N_LAS_Vector&>(result.getNativeVectorRef_());
  Epetra_MultiVector& epVecResult = 
    const_cast<Epetra_MultiVector&>(nonConstResult.epetraObj());
  
  int errorCode = ifpackPreconditionerPtr_->
    ApplyInverse(epVecInput, epVecResult);
  
  // Unset the transpose call
  if (useTranspose)
    ifpackPreconditionerPtr_->SetUseTranspose(false);    

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::getSolutionVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
N_NLS_NOX::Vector& SharedSystem::getSolutionVector()
{
  return *xyceSolnPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::getJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const N_LAS_Matrix& SharedSystem::getJacobian() const
{
  return *xyceJacobianPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::getJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
N_LAS_Matrix& SharedSystem::getJacobian(const Group* grp)
{
  ownerOfJacobian_ = grp;
  return *xyceJacobianPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::getStateVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::getStateVectors(const Group* grp)
{
  ownerOfStateVectors_ = grp;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::getLasSystem
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
N_LAS_System* SharedSystem::getLasSystem()
{
  return xyceLasSysPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::cloneSolutionVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
N_NLS_NOX::Vector* SharedSystem::cloneSolutionVector() const
{
  N_NLS_NOX::Vector* tmpVectorPtr = 0;
  tmpVectorPtr = 
    dynamic_cast<N_NLS_NOX::Vector*>(xyceSolnPtr_->clone(NOX::DeepCopy).release().get());

  if (tmpVectorPtr == 0) {
    const string message = 
      "N_NLS_NOX::SharedSystem::cloneSolutionVector() - dynamic cast/ memory allocation failure!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return (tmpVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem:: getNewtonVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const N_NLS_NOX::Vector & SharedSystem::getNewtonVector() const
{
  return *xyceNewtonPtr_;                                                       
}                                                                               

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::debugOutput1
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::debugOutput1
   (N_LAS_Matrix & jacobian, N_LAS_Vector & rhs)
{
  xyceInterfacePtr_->debugOutput1(jacobian, rhs);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::SharedSystem::debugOutput3
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void SharedSystem::debugOutput3 
   (N_LAS_Vector & dxVector, N_LAS_Vector & xVector)
{
  xyceInterfacePtr_->debugOutput3(dxVector, xVector);
}
#endif

