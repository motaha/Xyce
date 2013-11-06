//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Filename       : $$
//
// Purpose        : Constraint Backtracking Class.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 01/26/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $$
//
// Revision Date  : $$
//
// Current Owner  : $$
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------

#include <N_NLS_ConstraintBT.h>

#include <N_LAS_Vector.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>
#include <N_NLS_NLParams.h>

#ifdef Xyce_DEBUG_NONLINEAR
#include <N_ERH_ErrorMgr.h>
#endif

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::N_NLS_ConstraintBT
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
N_NLS_ConstraintBT::N_NLS_ConstraintBT()
 : constraintMinVector_(0),
   constraintMaxVector_(0),
   constraintChangeVector_(0),
   constraintTempVector_(0)
{

  // Initialize protected data
  resetThetaBoundNeg();
  resetThetaBoundPos();
  resetThetaChange();

}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::N_NLS_ConstraintBT
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
N_NLS_ConstraintBT::N_NLS_ConstraintBT(const N_NLS_ConstraintBT & right)
{

}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::~N_NLS_ConstraintBT
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
N_NLS_ConstraintBT::~N_NLS_ConstraintBT()
{
  if (constraintMinVector_ != 0)    delete constraintMinVector_;
  if (constraintMaxVector_ != 0)    delete constraintMaxVector_;
  if (constraintChangeVector_ != 0) delete constraintChangeVector_;
  if (constraintTempVector_ != 0)   delete constraintTempVector_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::operator=
// Purpose       : Assignment operator
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
N_NLS_ConstraintBT & N_NLS_ConstraintBT::operator=(
  const N_NLS_ConstraintBT & right)
{
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::operator==
// Purpose       : Equal operator
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
int N_NLS_ConstraintBT::operator==(const N_NLS_ConstraintBT & right) const

{
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::operator!=
// Purpose       : Not-Equal operator
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
int N_NLS_ConstraintBT::operator!=(const N_NLS_ConstraintBT & right) const

{
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::initializeAll
// Purpose       : Not-Equal operator
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/29/02
//-----------------------------------------------------------------------------
bool N_NLS_ConstraintBT::initializeAll 
    (N_LAS_System * lasSysPtr, const N_NLS_NLParams & nlParams)
{
  // create and initialize constraint backtracking vectors:
  constraintMinVector_ = lasSysPtr->builder().createVector();
  constraintMaxVector_ = lasSysPtr->builder().createVector();

  constraintMinVector_->putScalar(nlParams.getGlobalBTMin());
  constraintMaxVector_->putScalar(nlParams.getGlobalBTMax());

  constraintChangeVector_ = lasSysPtr->builder().createVector();

  constraintChangeVector_->putScalar (nlParams.getGlobalBTChange());

  constraintTempVector_ = lasSysPtr->builder().createVector();

  constraintTempVector_->putScalar(0.0);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::updateThetaBoundNeg
// Purpose       : Updates the minimum bound value for the backtracking
//                 algorithm.
// Special Notes : This is implemented according to a internal
//                 communication with John Shadid (SNL) on their MPSalsa
//                 constraint backtracking implementation.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
void N_NLS_ConstraintBT::updateThetaBoundNeg(const N_LAS_Vector * oldSoln,
                                             const N_LAS_Vector * solnUpdate)

{

#ifdef Xyce_DEBUG_NONLINEAR
  string msg = "N_NLS_ConstraintBT::updateThetaBoundNeg: ";
#endif

  N_LAS_Vector * solnPtr;
  N_LAS_Vector * solnUpdatePtr;

  // Initialize
  solnPtr       = const_cast<N_LAS_Vector *> (oldSoln);
  solnUpdatePtr = const_cast<N_LAS_Vector *> (solnUpdate);

  // First, form a vector of constraints...
  for (int i = 0; i < solnPtr->localLength(); ++i)
  {
    if ((*(solnUpdatePtr))[i] < 0.0)
      (*(constraintTempVector_))[i] = ((*(constraintMinVector_))[i] -
      (*(solnPtr))[i]) / (*(solnUpdatePtr))[i];
    else
      (*(constraintTempVector_))[i] = N_UTL_MachineDependentParams::DoubleMax();

#ifdef Xyce_DEBUG_NONLINEAR
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "min: ",
                           (*(constraintMinVector_))[i], "\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "soln: ",
                           (*(solnPtr))[i], "\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "solnUpdate: ",
                           (*(solnUpdatePtr))[i], "\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "constraint: ",
                           (*(constraintTempVector_))[i], "\n");
#endif
  }

  // Find minimum
  constraintTempVector_->minValue(&thetaBoundNeg_);

#ifdef Xyce_DEBUG_NONLINEAR
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg +
                         "thetaBoundNeg_: ", thetaBoundNeg_, "\n");
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::updateThetaBoundPos
// Purpose       : Updates the maximum bound value for the backtracking
//                 algorithm.
// Special Notes : This is implemented according to a internal
//                 communication with John Shadid (SNL) on their MPSalsa
//                 constraint backtracking implementation.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
void N_NLS_ConstraintBT::updateThetaBoundPos(const N_LAS_Vector * oldSoln,
                                             const N_LAS_Vector * solnUpdate)

{

#ifdef Xyce_DEBUG_NONLINEAR
  string msg = "N_NLS_ConstraintBT::updateThetaBoundPos: ";
#endif

  N_LAS_Vector * solnPtr;
  N_LAS_Vector * solnUpdatePtr;

  // Initialize
  solnPtr       = const_cast<N_LAS_Vector *> (oldSoln);
  solnUpdatePtr = const_cast<N_LAS_Vector *> (solnUpdate);

  // First, form a vector of constraints...
  for (int i = 0; i < solnPtr->localLength(); ++i)
  {
    if ((*(solnUpdatePtr))[i] > 0.0)
      (*(constraintTempVector_))[i] = ((*(constraintMaxVector_))[i] -
      (*(solnPtr))[i]) / (*(solnUpdatePtr))[i];
    else
      (*(constraintTempVector_))[i] = N_UTL_MachineDependentParams::DoubleMax();

#ifdef Xyce_DEBUG_NONLINEAR
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "max: ",
                           (*(constraintMaxVector_))[i], "\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "soln: ",
                           (*(solnPtr))[i], "\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "solnUpdate: ",
                           (*(solnUpdatePtr))[i], "\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "constraint: ",
                           (*(constraintTempVector_))[i], "\n");
#endif
  }

  // Find minimum
  constraintTempVector_->minValue(&thetaBoundPos_);

#ifdef Xyce_DEBUG_NONLINEAR
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg +
                         "thetaBoundPos_: ", thetaBoundPos_, "\n");
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_NLS_ConstraintBT::updateThetaChange
// Purpose       : Updates the percentage-change bound value for the
//                 backtracking algorithm.
// Special Notes : This is implemented according to a internal
//                 communication with John Shadid (SNL) on their MPSalsa
//                 constraint backtracking implementation.
//                 This function returns:
//
//                 theta_u = min {gamma_i | oldSoln_i | / | solnUpdate_i | }
//                            i
//
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
void N_NLS_ConstraintBT::updateThetaChange(const N_LAS_Vector * oldSoln,
                                           const N_LAS_Vector * solnUpdate)

{

#ifdef Xyce_DEBUG_NONLINEAR
  string msg = "N_NLS_ConstraintBT::updateThetaChange: ";
#endif

  N_LAS_Vector * solnPtr;
  N_LAS_Vector * solnUpdatePtr;

  // Initialize
  solnPtr       = const_cast<N_LAS_Vector *> (oldSoln);
  solnUpdatePtr = const_cast<N_LAS_Vector *> (solnUpdate);

  // First, form a vector of constraints...
  for (int i = 0; i < solnPtr->localLength(); ++i)
  {
    if (fabs((*(solnUpdatePtr))[i]) > N_UTL_MachineDependentParams::DoubleMin() &&
        fabs((*(solnPtr))[i]) > 0.0)
      (*(constraintTempVector_))[i] = (*(constraintChangeVector_))[i] *
        fabs((*(solnPtr))[i]) / fabs((*(solnUpdatePtr))[i]);

    else
      (*(constraintTempVector_))[i] = N_UTL_MachineDependentParams::DoubleMax();

#ifdef Xyce_DEBUG_NONLINEAR
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "change: ",
                           (*(constraintChangeVector_))[i], "\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "constraint: ",
                           (*(constraintTempVector_))[i], "\n");
#endif
  }

  // Find minimum
  constraintTempVector_->minValue(&thetaChange_);

#ifdef Xyce_DEBUG_NONLINEAR
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_DEBUG_0, msg + "thetaChange_: ",
                         thetaChange_, "\n");
#endif

}

