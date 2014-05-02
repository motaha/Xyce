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
// Filename       : $RCSfile: N_NLS_NOX_Group.C,v $
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
// Revision Number: $Revision: 1.32 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX_Group.h"
#include "N_NLS_NOX.h"
#include "N_NLS_NOX_Vector.h"
#include "N_NLS_NOX_SharedSystem.h"
#include "N_LAS_Vector.h"
#include "N_ERH_ErrorMgr.h"

// ----------   NOX Includes   ----------
#include "NOX_Abstract_Vector.H"

using namespace N_NLS_NOX;

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::Group
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Group::Group(SharedSystem& s) :
  sharedSystemPtr_(&s),
  xVecPtr_(Teuchos::rcp(dynamic_cast<N_NLS_NOX::Vector*>(s.cloneSolutionVector()))),
  xVec_(*xVecPtr_),
  fVecPtr_(Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>(xVec_.clone(NOX::ShapeCopy))),
  fVec_(*fVecPtr_),
  normF_(0.0),
  haveSolverFactors_(false)
{
  resetIsValid_();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::Group
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Group::Group(const Group& source, NOX::CopyType type) :
  sharedSystemPtr_(source.sharedSystemPtr_),
  xVecPtr_(Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>(source.getX().clone(type))),
  xVec_(*xVecPtr_),
  fVecPtr_(Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>(xVec_.clone(NOX::ShapeCopy))),
  fVec_(*fVecPtr_),
  normF_(0.0),
  haveSolverFactors_(false)
{
  // Default the is valid flags to "false"
  resetIsValid_();

  switch(type)
  {
  case NOX::DeepCopy:

    // Copy F
    if (source.isF()) {
      isValidF_ = true;
      fVec_ = source.fVec_;
      normF_ = source.normF_;
      if (sharedSystemPtr_->areStateVectors(&source))
        sharedSystemPtr_->getStateVectors(this);
    }

    // Take ownership of the Jacobian.
    if (source.isJacobian()) {
      //isValidJacobian_ = false;
      isValidJacobian_ = true;
      sharedSystemPtr_->getJacobian(this);
      haveSolverFactors_ = source.haveSolverFactors_;
    }

    // Copy Gradient Vector
    if (source.isGradient()) {

      if (Teuchos::is_null(gradVecPtr_)) {
        gradVecPtr_ = Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>
        (source.gradVecPtr_->clone(NOX::DeepCopy));
      }
      else
        *gradVecPtr_ = *(source.gradVecPtr_);
      
      isValidGradient_ = true;
    }

    // Copy Newton Vector
    if (source.isNewton()) {

      if (Teuchos::is_null(newtonVecPtr_)) {
        newtonVecPtr_ = Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>
        (source.newtonVecPtr_->clone(NOX::DeepCopy));
      }
      else
        *newtonVecPtr_ = *(source.newtonVecPtr_);
      
      isValidNewton_ = true;
    }

    break;

  case NOX::ShapeCopy:
    resetIsValid_();
    break;

  default:
    const std::string message =  "N_NLS::NOX::Group::Group() - Invalid ConstructorType for group copy constructor";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::~Group
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Group::~Group()
{
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::operator=
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group& Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&>(source));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::operator=
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group& Group::operator=(const Group& source)
{
  resetIsValid_();

  // Copy solution vector
  xVec_ = source.xVec_;
  
  // Copy residual vector and take ownership of state vectors
  if (source.isF()) {
    isValidF_ = true;
    fVec_ = source.fVec_;
    normF_ = source.normF_;
    if (source.sharedSystemPtr_->areStateVectors(&source))
      sharedSystemPtr_->getStateVectors(this);
  }
  
  // Copy the Jacobian by taking ownership of the shared system
  if (source.isJacobian()) {
    isValidJacobian_ = true;
    sharedSystemPtr_->getJacobian(this);
    haveSolverFactors_ = source.haveSolverFactors_;
  }

  // Copy Gradient Vector
  if ( source.isGradient() ) {

    if (Teuchos::is_null(gradVecPtr_)) {
      gradVecPtr_ = Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>
	(source.gradVecPtr_->clone(NOX::DeepCopy));
    }
    else
      *gradVecPtr_ = *(source.gradVecPtr_);
    
    isValidGradient_ = true;
  }
  
  // Copy Newton Vector
  if (source.isNewton()) {
    if (Teuchos::is_null(newtonVecPtr_)) {
      newtonVecPtr_ = Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>
        (source.newtonVecPtr_->clone(NOX::DeepCopy));
    }
    else
      *newtonVecPtr_ = *(source.newtonVecPtr_);
    
    isValidNewton_ = true;
  }
  
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setX
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void Group::setX(const NOX::Abstract::Vector& input)
{
  setX(dynamic_cast<const Vector&>(input));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setX
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void Group::setX(const Vector& input)
{
  resetIsValid_();
  xVec_ = input;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::computeX
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void Group::computeX(const NOX::Abstract::Group& grp, 
		     const NOX::Abstract::Vector& d, 
		     double step)
{
  computeX(dynamic_cast<const Group&>(grp), dynamic_cast<const Vector&>(d), 
	   step);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::computeX
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void Group::computeX(const Group& grp, const Vector& d, double step)
{
  resetIsValid_();
  xVec_.update(step, d, 1.0, grp.getX(), 0.0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::computeF
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::computeF()
{
  if( isF() ) return Ok;

  isValidF_ = sharedSystemPtr_->computeF(xVec_, fVec_, this);

  // Xyce computes "-f" for efficiency of Newton solve:
  // "Js = f" instead of "Js = -f" We need the real F!
  fVec_.scale(-1.0);

  normF_ = fVec_.norm();

  return (isF()?Ok:Failed);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::computeJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::computeJacobian()
{
  if (isJacobian()) return Ok;

  isValidJacobian_ = sharedSystemPtr_->computeJacobian(this);

  haveSolverFactors_ = false;

  return (isJacobian()?Ok:Failed);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::computeGradient
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::computeGradient()
{
  if (isGradient()) return Ok;

  if (!isF()) 
    throwError("computeGradient", "F is not Valid!");

  if (!isJacobian())
    throwError("computeGradient", "Jacobian is not Valid!");

  if (Teuchos::is_null(gradVecPtr_))
    gradVecPtr_ = Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>
      (fVec_.clone(NOX::ShapeCopy));

  NOX::Abstract::Group::ReturnType status = 
    applyJacobianTranspose(fVec_, *gradVecPtr_);

  //isValidGradient_ = sharedSystemPtr_->computeGradient(fVec_, *gradVecPtr_);
  if (status == Ok)
    isValidGradient_ = true;
  else 
    isValidGradient_ = false;

  return (isGradient()?Ok:Failed);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::computeNewton
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::computeNewton(Teuchos::ParameterList& params)
{
  if (isNewton()) return Ok;

  if (!isF())
    throwError("computeNewton", "F is not Valid!");

  if (!isJacobian())
    throwError("computeNewton", "Jacobian is not Valid!");

  if (Teuchos::is_null(newtonVecPtr_)) 
    newtonVecPtr_ = Teuchos::rcp_dynamic_cast<N_NLS_NOX::Vector>
      (fVec_.clone(NOX::ShapeCopy));

  isValidNewton_ = sharedSystemPtr_->computeNewton(fVec_, *newtonVecPtr_, params);

  haveSolverFactors_ = true;

  newtonVecPtr_->scale(-1.0);

  return (isNewton()?Ok:Failed);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::applyJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
{
  return applyJacobian(dynamic_cast<const Vector&>(input), dynamic_cast<Vector&>(result));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::applyJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::applyJacobian(const Vector& input, Vector& result) const
{
  if (!isJacobian()) {
    throwError("applyJacobian", "Jacobian is not Valid!");

    //RPP Hack to get Homotopy working!!
    //cout << "RPP: Inefficient Jacobian computation!!" << endl; 
    //(const_cast<N_NLS_NOX::Group*>(this))->computeF();
    //(const_cast<N_NLS_NOX::Group*>(this))->computeJacobian();
  }

  return (sharedSystemPtr_->applyJacobian(input, result)?Ok:Failed);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::applyJacobianTranspose
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::applyJacobianTranspose(const NOX::Abstract::Vector& input,
				   NOX::Abstract::Vector& result) const
{
  return applyJacobianTranspose(dynamic_cast<const Vector&>(input),
				dynamic_cast<Vector&>(result));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::applyJacobianTranspose
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  if (!isJacobian())
    throwError("applyJacobianTranspose", "Jacobian is not Valid!");

  return (sharedSystemPtr_->applyJacobianTranspose(input, result)?Ok:Failed);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::applyJacobianInverse
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType 
Group::applyJacobianInverse(Teuchos::ParameterList& params,
			    const NOX::Abstract::Vector& input,
			    NOX::Abstract::Vector& result) const
{
  return applyJacobianInverse(params, dynamic_cast<const Vector&>(input),
			      dynamic_cast<Vector&>(result));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::applyJacobianInverse
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType 
Group::applyJacobianInverse(Teuchos::ParameterList& params,
			    const Vector& input, Vector& result) const
{
  if (!isJacobian())
    throwError("applyJacobianInverse", "Jacobian is not Valid!");

  bool status = sharedSystemPtr_->computeNewton(input, result, params);

  haveSolverFactors_ = true;

  return (status ? Ok : Failed);
}
  
//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::applyRightPreconditioning
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType 
Group::applyRightPreconditioning(bool useTranspose,
				 Teuchos::ParameterList& params,
				 const NOX::Abstract::Vector& input, 
				 NOX::Abstract::Vector& result) const
{
  return applyRightPreconditioning(useTranspose, params,
				   dynamic_cast<const Vector&>(input), 
				   dynamic_cast<Vector&> (result));
}
 
//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::applyRightPreconditioning
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType 
Group::applyRightPreconditioning(bool useTranspose,
				 Teuchos::ParameterList& params,
				 const Vector& input, 
				 Vector& result) const
{
  if (!isJacobian()) {
    // throw error - Jacobian is not owned by this group 
  }

  if (!isValidPreconditioner_) {
    sharedSystemPtr_->computePreconditioner();
    isValidPreconditioner_ = true;
  }
  sharedSystemPtr_->applyRightPreconditioning(useTranspose, params,
					      input, result);
  
  return NOX::Abstract::Group::Ok;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::isF
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool Group::isF() const
{
  return isValidF_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::isJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool Group::isJacobian() const
{
  return (isValidJacobian_ && sharedSystemPtr_->isJacobianOwner(this));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::isGradient
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool Group::isGradient() const
{
  return isValidGradient_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::isNewton
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool Group::isNewton() const
{
  return isValidNewton_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getX
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const NOX::Abstract::Vector& Group::getX() const
{
  return xVec_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getF
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const NOX::Abstract::Vector& Group::getF() const
{
  if (!isF()) {
    throwError("getF", 
	       "F is not current with respect to the solution vector!");
  }

  return fVec_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getNormF
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
double Group::getNormF() const
{
  if (!isF()) {
    //cout << "N_NLS_NOX::Group::getNormF() - F is not current "
    // << "with respect to the solution vector!" << endl;
    //throw "NOX Error";
    (const_cast<N_NLS_NOX::Group*>(this))->computeF();
  }

  return normF_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getGradient
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const NOX::Abstract::Vector& Group::getGradient() const
{
  if (Teuchos::is_null(gradVecPtr_))
    throwError("getGradient", "gradVecPtr_ is 0!");
  
  return *gradVecPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getNewton
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
const NOX::Abstract::Vector& Group::getNewton() const
{
  if (Teuchos::is_null(newtonVecPtr_))
  {
   throwError("getNewton", "newtonVecPtr_ is 0!");
  }
  return *newtonVecPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::clone
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<NOX::Abstract::Group> 
Group::clone(NOX::CopyType type) const
{
  Teuchos::RefCountPtr<Group> ptr = Teuchos::rcp(new Group(*this, type));
  return ptr;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::resetIsValid
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
// resets the isValid flags to false
void Group::resetIsValid_()
{
  isValidF_ = false;
  isValidGradient_ = false;
  isValidJacobian_ = false;
  isValidNewton_ = false;
  isValidPreconditioner_ = false;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::throwError
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void Group::throwError(std::string method, std::string message) const
{
  const std::string leader = "N_NLS::NOX::Group::";
  const std::string fcn = "() - ";

  std::string error = leader + method + fcn + message;

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, error);
}


