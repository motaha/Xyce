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
// Filename      : $RCSfile: N_MPDE_DeviceInterface.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.6.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:47 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_DEV_DeviceInterface.h>

#include <N_MPDE_DeviceInterface.h>

//-----------------------------------------------------------------------------
// Class         : N_MPDE_DeviceInterface
// Purpose       :
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
N_MPDE_DeviceInterface::N_MPDE_DeviceInterface()
{
}


//-----------------------------------------------------------------------------
// Class         : N_MPDE_DeviceInterface
// Purpose       :
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::registerDeviceInterface( RefCountPtr<N_DEV_DeviceInterface> devIntPtr )
{
  devInterfacePtr_ = devIntPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::registerLinearSystem (RefCountPtr<N_LAS_System> tmp_system_ptr)
{
  devInterfacePtr_->registerLinearSystem(&*tmp_system_ptr);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::registerAnalysisInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::registerAnalysisInterface (RefCountPtr<N_ANP_AnalysisInterface> tmp_anaIntPtr)
{
  devInterfacePtr_->registerAnalysisInterface(&*tmp_anaIntPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::registerNonlinearSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::registerNonlinearSolver (RefCountPtr<N_NLS_Manager> tmp_nlsMgrPtr)
{
  devInterfacePtr_->registerNonlinearSolver(&*tmp_nlsMgrPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::getFastSourcePeriod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
vector<double> N_MPDE_DeviceInterface::getFastSourcePeriod (vector<string>& sourceNames)
{
  return devInterfacePtr_->getFastSourcePeriod (sourceNames);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::registerFastSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
vector<double> N_MPDE_DeviceInterface::registerFastSources(vector<string>& sourceNames)
{
  return devInterfacePtr_->registerFastSources (sourceNames);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::deactivateSlowSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::deactivateSlowSources() 
{
  devInterfacePtr_->deactivateSlowSources();
}
  
 
//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::activateSlowSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::activateSlowSources() 
{
  devInterfacePtr_->activateSlowSources();
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::setMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::setMPDEFlag( bool flagVal ) 
{
  devInterfacePtr_->setMPDEFlag( flagVal );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::setBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::setBlockAnalysisFlag( bool flagVal ) 
{
//  devInterfacePtr_->setMPDERunFlag( flagVal );
  devInterfacePtr_->setBlockAnalysisFlag( flagVal );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_DeviceInterface::setFastTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
void N_MPDE_DeviceInterface::setFastTime( double timeVal ) 
{
  devInterfacePtr_->setFastTime( timeVal );
}

