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
// Filename      : $RCSfile: N_TIA_NoTimeIntegration.C,v $
//
// Purpose       : This file contains the functions which define the
//		   time integration methods classes.
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 7/21/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.43 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_NoTimeIntegration.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TwoLevelError.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>

//Need Epetra Stuff for now, since directly manipulating Epetra Objs for
//obtainJacobian.  This will be abstracted later
#include <Epetra_Import.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::N_TIA_NoTimeIntegration
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
N_TIA_NoTimeIntegration::N_TIA_NoTimeIntegration
    (N_TIA_TIAParams & tiaP,
     N_TIA_StepErrorControl & secTmp,
     N_TIA_DataStore & dsTmp)
     
  : N_TIA_TimeIntegrationMethod(tiaP,secTmp,dsTmp)
{
  leadingCoeff = 1.0;
  alphas = -1.0;
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::~N_TIA_NoTimeIntegration
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
N_TIA_NoTimeIntegration::~N_TIA_NoTimeIntegration() {}


//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::obtainCorrectorDeriv
// Purpose       : Evaluate the Corrector Derivative Formula
// Special Notes : For "no integration" the derivatives should always be zero.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
void N_TIA_NoTimeIntegration::obtainCorrectorDeriv()
{
  ds.nextSolutionDerivPtr->putScalar(0.0);
  ds.nextStateDerivPtr->putScalar(0.0);
  ds.nextStoreLeadCurrQCompDerivPtr->putScalar(0.0);
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::obtainResidual
//
// Purpose       : This function returns the residual for the steady state
//                 case.
//
// Special Notes : For "no integration" the derivatives should always be zero,
//                 so don't add in the dqdt term.
//
//                 This function is only called in the new-DAE case.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/09/04
//-----------------------------------------------------------------------------
void N_TIA_NoTimeIntegration::obtainResidual ()
{
  ds.RHSVectorPtr->linearCombo
    (0.0,*ds.RHSVectorPtr,+1.0,*ds.daeFVectorPtr);

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  ds.RHSVectorPtr->scale(-1.0);

  // if voltage limiting is on, add it in:
  if (ds.limiterFlag)
  {
    (ds.RHSVectorPtr)->daxpy(
      *(ds.RHSVectorPtr), +1.0, *(ds.dFdxdVpVectorPtr));
  }
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::obtainJacobian
//
// Purpose       : Returns the full Jacobian matrix for the steady state
//                 case.  
//
// Special Notes : For "no integration" the derivatives should always be zero.
//                 However, to prevent singular matrices, the dQdt term is
//                 added in anyway, with an assumed very large time step.
//
//                 There may be a better way to handle this - the assumed
//                 large time step leads to a very small alpha/dt, which
//                 could (maybe?) have an adverse effect on the matrix
//                 conditioning.
//
//                 The singular matrix problem will happen for the case of
//                 capacitors in serial.  Any voltage node which is *only*
//                 connected to capacitors will not really be part of the
//                 system of equations.  Ideally nodes like this would just
//                 be removed from the system altogether.  No current is
//                 passing through them, and they are really just floating
//                 in space.
//
//                 For the time being, however, what we do is put some
//                 bogus C*alpha/dt terms into the Jacobian.  If dt is
//                 large, then the C*alpha/dt is very small, and doesn't
//                 significantly change any Jacobian entries, except for
//                 zero entries.
//
//                 This function is only called in the new-DAE case.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/09/04
//-----------------------------------------------------------------------------
void N_TIA_NoTimeIntegration::obtainJacobian ()
{
  N_LAS_Matrix & dQdx = *(ds.dQdxMatrixPtr);
  N_LAS_Matrix & dFdx = *(ds.dFdxMatrixPtr);
  N_LAS_Matrix & Jac = *(ds.JMatrixPtr);

  Jac.linearCombo( 1.0e-20, dQdx, 1.0, dFdx );
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::applyJacobian
//
// Purpose       : Applies the Jacobian operator for the steady state
//                 case.  
//
// Special Notes : 
//                 This function is only called in the new-DAE HB (matrix-free) case.
//
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
void N_TIA_NoTimeIntegration::applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result)
{
  N_LAS_Vector & dQdxV = *(ds.dQdxVecVectorPtr);
  N_LAS_Vector & dFdxV = *(ds.dFdxVecVectorPtr);
  result.linearCombo( 1.0e-20, dQdxV, 1.0, dFdxV );
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::getInitialQnorm
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
void N_TIA_NoTimeIntegration::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  tle.q1HistorySum = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::setupTwoLevelError
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
void N_TIA_NoTimeIntegration::setupTwoLevelError(N_TIA_TwoLevelError & tle)
{
  tle.xErrorSum    = 0.0;
  tle.qErrorSum    = 0.0;
  tle.xErrorSum_m1 = 0.0;
  tle.xErrorSum_m2 = 0.0;
  tle.innerSize    = ds.globalLength ();
}

