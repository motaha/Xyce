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
// Filename       : $RCSfile: N_LAS_HBPrecondFactory.h,v $
//
// Purpose        : Preconditioner Factory for HB
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 10/01/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:44 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_HBPrecondFactory_h
#define Xyce_N_LAS_HBPrecondFactory_h

// ---------- Standard Includes ----------

#include <string>

// ----------   Xyce Includes   ----------

#include <N_LAS_PrecondFactory.h>
#include <N_ERH_ErrorMgr.h>
#include <N_DEV_DeviceInterface.h>
#include <N_UTL_OptionBlock.h>

// ----------  Fwd Declares  -------------

class N_LAS_Preconditioner;
class N_LAS_Problem;
class N_LOA_Loader;
class N_MPDE_State;
class N_LAS_HBBuilder;
class N_LOA_HBLoader;
class N_LAS_Builder;
class N_LAS_System;

//-----------------------------------------------------------------------------
// Class         : N_LAS_HBPrecondFactory
// Purpose       : 
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
class N_LAS_HBPrecondFactory : public N_LAS_PrecondFactory
{
public:

  // Default Constructor, sets null preconditioner.
  N_LAS_HBPrecondFactory();

  // Basic Constructor, sets preconditioner factory options.
  N_LAS_HBPrecondFactory( const N_UTL_OptionBlock & OB );

  // Destructor
  virtual ~N_LAS_HBPrecondFactory() {} 

  // Creates a new preconditioner (matrix based).
  // NOTE:  This type of creation is not supported by this preconditioner factory.
  Teuchos::RCP<N_LAS_Preconditioner> create( const Teuchos::RCP<N_LAS_Problem> & problem )
  {
    std::string msg = "N_LAS_HBPrecondFactory::create() using N_LAS_Problem is not supported!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR_0, msg);
    return Teuchos::null;
  }

  // Creates a new preconditioner (matrix free).
  Teuchos::RCP<N_LAS_Preconditioner> create( const Teuchos::RCP<N_LAS_System> & lasSystem );

  // Set the fast times being used in the HB analysis.
  void setFastTimes( const std::vector<double> & times )
    { times_ = times; }

  // Register MPDE state
  void registerMPDEState( const Teuchos::RCP<N_MPDE_State>& mpdeState )
    { statePtr_ = mpdeState; }

  // Register the application system loader
  void registerAppLoader( const Teuchos::RCP<N_LOA_Loader>& appLoaderPtr ) 
    { appLoaderPtr_ = appLoaderPtr; }

  // Register the application system builder
  void registerAppBuilder( const Teuchos::RCP<N_LAS_Builder>& appBuilderPtr ) 
    { appBuilderPtr_ = appBuilderPtr; }

  // Register the application system loader
  void registerHBLoader( const Teuchos::RCP<N_LOA_HBLoader>& hbLoaderPtr ) 
    { hbLoaderPtr_ = hbLoaderPtr; }

  // Register the HB builder 
  void registerHBBuilder( const Teuchos::RCP<N_LAS_HBBuilder>& hbBuilder ) 
    { hbBuilderPtr_ = hbBuilder; }

  // Register the interface to the devices
  void registerDeviceInterface( const Teuchos::RCP< N_DEV_DeviceInterface >& devInterfacePtr )
    { devInterfacePtr_ = devInterfacePtr; }

private:

  std::vector<double> times_;
  std::string precType_;
  Teuchos::RCP<N_MPDE_State> statePtr_;
  Teuchos::RCP<N_LOA_Loader> appLoaderPtr_;
  Teuchos::RCP<N_LAS_Builder> appBuilderPtr_;
  Teuchos::RCP<N_LAS_System> lasSysPtr_;
  Teuchos::RCP<N_LOA_HBLoader> hbLoaderPtr_;
  Teuchos::RCP<N_LAS_HBBuilder> hbBuilderPtr_;
  Teuchos::RCP<N_DEV_DeviceInterface> devInterfacePtr_;
  Teuchos::RCP<const N_UTL_OptionBlock> OB_;

  // Copy constructor.
  N_LAS_HBPrecondFactory( const N_LAS_HBPrecondFactory& pf );
};


#endif

