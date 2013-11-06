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

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_TIA_TimeIntegrationMethods.C,v $
//
// Purpose       : This file contains the functions which define the
//		             time integration methods classes.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.47.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:50 $
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

#include <stdio.h>

// ----------   Xyce Includes   ----------
#include <N_TIA_TimeIntegrationMethods.h>

#include <N_TIA_BackwardDifferentiation15.h>
#include <N_TIA_Gear12.h>
#include <N_TIA_OneStep.h>
#include <N_TIA_NoTimeIntegration.h>

#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TIAParams.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_System.h>

// ---------- Static Initializations ----------
//-----------------------------------------------------------------------------
// Function      : N_TIA_WorkingIntegrationMethod::N_TIA_WorkingIntegrationMethod
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
N_TIA_WorkingIntegrationMethod::
  N_TIA_WorkingIntegrationMethod
   (N_TIA_TIAParams & tiaP,
    N_TIA_StepErrorControl & secTmp,
    N_TIA_DataStore & dsTmp)
  : integMethodPtr(0),
    tiaParams(tiaP),
    sec(secTmp),
    ds(dsTmp)
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_WorkingIntegrationMethod::N_TIA_WorkingIntegrationMethod
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
N_TIA_WorkingIntegrationMethod::N_TIA_WorkingIntegrationMethod(
  const unsigned int integration_method,
  N_TIA_TIAParams & tiaP,
  N_TIA_StepErrorControl & secTmp,
  N_TIA_DataStore & dsTmp)
  :
  workingIntegMethod(integration_method),
  integMethodPtr(0),
  tiaParams(tiaP),
  sec(secTmp),
  ds(dsTmp)
{
  createTimeIntegMethod(integration_method);
  return;
}

//-----------------------------------------------------------------------------
// Function      :
//            N_TIA_WorkingIntegrationMethod::~N_TIA_WorkingIntegrationMethod()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
N_TIA_WorkingIntegrationMethod::~N_TIA_WorkingIntegrationMethod()
{
  if (integMethodPtr != 0) 
  {
    delete integMethodPtr;
    integMethodPtr = 0; 
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_WorkingIntegrationMethod::createTimeIntegMethod
// Purpose       : Creates the time integration method class --- assigning a
//                 pointer and the Leading Coefficient value of the method.
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
void N_TIA_WorkingIntegrationMethod::createTimeIntegMethod(
  const unsigned int integration_method)
{
#ifdef Xyce_DEBUG_TIME
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "\n  ********** createTimeIntegMethod Function Called");
#endif

  N_TIA_TimeIntegrationMethod* pTIM;
  char ch_msg[256];

  workingIntegMethod = integration_method;

  switch(integration_method)
  {
    case TIAMethod_NONE:
      pTIM = N_TIA_NoTimeIntegration::factory(tiaParams,sec,ds);
      break;
    case TIAMethod_BACKWARD_DIFFERENTIATION_15:
      pTIM = N_TIA_BackwardDifferentiation15::factory(tiaParams,sec,ds);
      break;
    case TIAMethod_GEAR_12:
      pTIM = N_TIA_Gear12::factory(tiaParams,sec,ds);
      break;
    case TIAMethod_ONESTEP:
      pTIM = N_TIA_OneStep::factory(tiaParams,sec,ds);
      break;

    // deprecated old-DAE methods:
    case TIAMethod_BACKWARD_EULER:
    case TIAMethod_BACKWARD_DIFFERENTIATION_2:
    case TIAMethod_TRAPEZOIDAL:
      sprintf(ch_msg,
        "N_TIA_WorkingIntegrationMethod::createTimeIntegMethod.  Non-valid method specified.  integration_method = %d",
        integration_method);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, string(ch_msg));
      break;

    default:
      sprintf(ch_msg,
        "N_TIA_WorkingIntegrationMethod::createTimeIntegMethod.  Non-valid method specified.  integration_method = %d",
        integration_method);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, string(ch_msg));
      break;
  }

  if( integMethodPtr ) { delete integMethodPtr; integMethodPtr = 0; }

  integMethodPtr = pTIM;

#ifdef Xyce_VERBOSE_TIME
  printWorkingIntegMethod();
#endif

  return;
}

#ifdef Xyce_VERBOSE_TIME
//-----------------------------------------------------------------------------
// Function      : N_TIA_WorkingIntegrationMethod::printWorkingIntegMethod
// Purpose       : This function is a debug output function.  It prints
//                 to the screen the current integration method.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void N_TIA_WorkingIntegrationMethod::printWorkingIntegMethod()
{

  string msg = "  Integration method = ";
  
  switch (workingIntegMethod)
  {
    case TIAMethod_NONE:
      msg = "  Integration method = None\n";
      break;
    case TIAMethod_BACKWARD_DIFFERENTIATION_15:
      msg = "  Integration method = Backward Differentiation 15\n";
      break;
    case TIAMethod_ONESTEP:
      msg = "Onestep: Trapezoidal\n";
      break;
    case TIAMethod_GEAR_12:
      msg = "  Integration method = Gear 12\n";
      break;
    default:
      msg = "N_TIA_WorkingIntegrationMethod::printWorkingIntegMethod  ";
      msg += "Time Integration method not specified correctly.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      break;
  }

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);
}
#endif

//*****************************************************************************
//************* Functions for Time Integration Method Base class  *************
//*****************************************************************************

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::N_TIA_TimeIntegrationMethod
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, 9233: Computational Sciences
// Creation Date : 9/11/00
//-----------------------------------------------------------------------------

N_TIA_TimeIntegrationMethod::N_TIA_TimeIntegrationMethod
  ( N_TIA_TIAParams & tiaP, 
    N_TIA_StepErrorControl & secTmp,
    N_TIA_DataStore & dsTmp) 
  : leadingCoeff(1.0),
    tiaParams(tiaP),
    sec(secTmp),
    ds(dsTmp)
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::~N_TIA_TimeIntegrationMethod
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
N_TIA_TimeIntegrationMethod::~N_TIA_TimeIntegrationMethod() 
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::partialTimeDeriv
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 1/20/07
//-----------------------------------------------------------------------------
double N_TIA_TimeIntegrationMethod::partialTimeDeriv()
{
  if (sec.currentTimeStep < 1e-30) 
  {
    string msg = 
      "Excessively small current time step in N_TIA_TimeIntegrationMethods.h, incorrectly returning with large value";

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);
    return (leadingCoeff * 1.e+30);
  }
  return (leadingCoeff / sec.currentTimeStep);
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::obtainResidual
// Purpose       : Evaluate residual for nonlinear solver
// Special Notes : 
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/08/04
//-----------------------------------------------------------------------------
void N_TIA_TimeIntegrationMethod::obtainResidual()
{
  string msg = "N_TIA_ControlMethod::obtainResidual";
  msg += " The current algorithm does not have an implemented";
  msg += " obtainResidual function.\n";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::rejectStep
// Purpose       : restore history & choose new step-size & order
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
void N_TIA_TimeIntegrationMethod::rejectStep()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::completeStep
// Purpose       : update history & choose new step-size & order
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
void N_TIA_TimeIntegrationMethod::completeStep()
{
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::printOutputSolution()
// Purpose       : Print output that is dumbed down in terms of order.
//
// Special Notes : This routine picks smaller time steps to approximate first
//                 order integration from the perspective of the output.
//                 
//                 For the old method classes, this function does not do
//                 any interpolation (none is possible), and just calls 
//                 the output manager adapter directly.
//
//                 ERK:  Note, the old methods (old-DAE) have all been 
//                 removed from Xyce, so possibly this function isn't needed
//                 anymore.
//
// Scope         : public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 11/22/05
//-----------------------------------------------------------------------------
bool N_TIA_TimeIntegrationMethod::printOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr, 
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const bool doNotInterpolate,
                  const vector <double> & outputInterpolationTimes,
                  bool skipPrintLineOutput )
{
#ifdef Xyce_DEBUG_TIME
  cout << "Calling conventional outputs!" << endl;
#endif

  outputMgrAdapterRCPtr->tranOutput( time, *solnVecPtr, *ds.currStatePtr, *ds.currStorePtr, skipPrintLineOutput  );

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::saveOutputSolution
// Purpose       : 
// Special Notes : For the old method functions (old-DAE/ODE) no interpolation
//                 is possible, so this function calls directly through to the
//                 output manager.
// Scope         : public
// Creator       : Eric Keiter, SNL, 1437
// Creation Date : 10/25/07
//-----------------------------------------------------------------------------
bool N_TIA_TimeIntegrationMethod::saveOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  N_LAS_Vector * solnVecPtr,
                  const double saveTime,
                  const bool doNotInterpolate) 
{
  outputMgrAdapterRCPtr->outputDCOP( *(solnVecPtr) );
  return true;
}

