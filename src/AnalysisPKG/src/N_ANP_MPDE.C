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
// Filename      : $RCSfile: N_ANP_MPDE.C,v $
// Purpose       : MPDE analysis functions.
// Special Notes :
// Creator       : Todd Coffey, 1414, Ting Mei, 1437
// Creation Date : 07/23/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.9.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:31 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------


// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_TIA_DataStore.h>
#include <N_LOA_Loader.h>
#include <N_MPDE_Manager.h>
#include <N_IO_OutputMgr.h>

#include <N_IO_CmdParse.h>

#include<N_ANP_MPDE.h>

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::N_ANP_MPDE( N_ANP_AnalysisManager * )
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
N_ANP_MPDE::N_ANP_MPDE( N_ANP_AnalysisManager * anaManagerPtr ) :
  N_ANP_AnalysisBase(anaManagerPtr),
  isPaused(false),
  startDCOPtime(0.0),  
  endTRANtime(0.0)
{  
  devInterfacePtr_ = anaManagerRCPtr_->devInterfacePtr;
  topoMgrPtr_ = anaManagerRCPtr_->topoMgrPtr;
  nonlinearEquationLoaderPtr_ = anaManagerRCPtr_->nonlinearEquationLoaderPtr;
  appBuilderPtr_ = anaManagerRCPtr_->appBuilderPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::run()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::run()
{
  anaManagerRCPtr_->mpdeMgrPtr_->run();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::init()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::init()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::loopProcess()
// Purpose       : Conduct the time stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::loopProcess()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::processSuccessfulDCOP()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::processSuccessfulDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::processSuccessfulStep()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::processSuccessfulStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::processFailedStep
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::processFailedStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::processFailedDCOP
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::processFailedDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::finish
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::finish()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::resetForStepAnalysis() 
// Purpose       : When doing a .STEP sweep, some data must be reset to its 
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::resetForStepAnalysis()
{
  totalNumberSuccessStepsThisParameter_ = 0;
  return false;
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::finalVerboseOutput
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool N_ANP_MPDE::finalVerboseOutput()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::takeAnIntegrationStep_
// Purpose       : Take a transient integration step.
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void N_ANP_MPDE::takeAnIntegrationStep_()
{
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::printStepHeader()
// Purpose       : Prints out time step information.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
void N_ANP_MPDE::printStepHeader() 
{
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_MPDE::printProgress()
// Purpose       : Outputs run completion percentage and estimated
//                 time-to-completion.
//
// Special Notes : This will need some fixing to work with .STEP.
//
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 06/07/2002
//-----------------------------------------------------------------------------
void N_ANP_MPDE::printProgress() 
{
}

