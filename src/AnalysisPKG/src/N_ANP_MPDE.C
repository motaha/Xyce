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
// Filename      : $RCSfile: N_ANP_MPDE.C,v $
// Purpose       : MPDE analysis functions.
// Special Notes :
// Creator       : Todd Coffey, 1414, Ting Mei, 1437
// Creation Date : 07/23/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.14 $
// Revision Date  : $Date: 2014/02/24 23:49:12 $
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

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : MPDE::MPDE( AnalysisManager * )
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
MPDE::MPDE( AnalysisManager * anaManagerPtr ) :
  AnalysisBase(anaManagerPtr),
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
// Function      : MPDE::run()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::run()
{
  anaManagerRCPtr_->getMPDEManager()->run();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::init()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::init()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::loopProcess()
// Purpose       : Conduct the time stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::loopProcess()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::processSuccessfulDCOP()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::processSuccessfulDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::processSuccessfulStep()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::processSuccessfulStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::processFailedStep
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::processFailedStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::processFailedDCOP
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::processFailedDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::finish
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::finish()
{
  return false;
}

bool MPDE::handlePredictor()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::resetForStepAnalysis() 
// Purpose       : When doing a .STEP sweep, some data must be reset to its 
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool MPDE::resetForStepAnalysis()
{
  totalNumberSuccessStepsThisParameter_ = 0;
  return false;
}


//-----------------------------------------------------------------------------
// Function      : MPDE::finalVerboseOutput
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date : 3/11/06
//-----------------------------------------------------------------------------
bool MPDE::finalVerboseOutput()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MPDE::takeAnIntegrationStep_
// Purpose       : Take a transient integration step.
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void MPDE::takeAnIntegrationStep_()
{
}

} // namespace Analysis
} // namespace Xyce
