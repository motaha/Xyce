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
// Filename       : $RCSfile: N_MPDE_DeviceInterface.h,v $
//
// Purpose        : MPDE Specific Loader
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 07/21/2008
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:46 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_MPDE_DeviceInterface_H
#define Xyce_N_MPDE_DeviceInterface_H

// ---------- Standard Includes ----------
#include <string>
#include <vector>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>

// ---------- Forward declarations --------
class N_ANP_AnalysisInterface;
class N_NLS_Manager;
class N_LAS_System;


//-----------------------------------------------------------------------------
// Class         : N_MPDE_DeviceInterface
// Purpose       :
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
class N_MPDE_DeviceInterface
{

public:
  N_MPDE_DeviceInterface();

  void registerDeviceInterface( RefCountPtr<N_DEV_DeviceInterface> devIntPtr );
  void registerLinearSystem(RefCountPtr<N_LAS_System> tmp_system_ptr);
  void registerAnalysisInterface(RefCountPtr<N_ANP_AnalysisInterface> tmp_anaIntPtr);
  void registerNonlinearSolver (RefCountPtr<N_NLS_Manager> tmp_nlsMgrPtr);

  vector<double> getFastSourcePeriod (vector<string>& sourceNames);
  vector<double> registerFastSources (vector<string> & sourceNames);

  void deactivateSlowSources();
  void activateSlowSources();

  void setMPDEFlag( bool flagVal );
  void setBlockAnalysisFlag( bool flagVal );
  void setFastTime( double timeVal );

private:
  RefCountPtr< N_DEV_DeviceInterface > devInterfacePtr_;

};

#endif

