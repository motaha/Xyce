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
// Filename       : $RCSfile: N_LOA_LoaderMgr.C,v $
//
// Purpose        : This file contains functions for the loader
//                  manager class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Xyce Includes   ----------
#include <N_LOA_LoaderMgr.h>
#include <N_LOA_Loader.h>
#include <N_LOA_CktLoader.h>
#include <N_LOA_NonlinearEquationLoader.h>

//-----------------------------------------------------------------------------
// Function      : N_LOA_LoaderMgr::N_LOA_LoaderMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
N_LOA_LoaderMgr::N_LOA_LoaderMgr ()
{

}

//-----------------------------------------------------------------------------
// Function      : N_LOA_LoaderMgr::~N_LOA_LoaderMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
N_LOA_LoaderMgr::~N_LOA_LoaderMgr ()
{

}

//-----------------------------------------------------------------------------
// Function      : N_LOA_LoaderMgr::createLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
N_LOA_Loader * N_LOA_LoaderMgr::createLoader (int loaderIndex)
{
  // for now, only option is the circuit loader:
  loaderPtr_ = new N_LOA_CktLoader ();

  return loaderPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_LoaderMgr::createCktLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
N_LOA_CktLoader * N_LOA_LoaderMgr::createCktLoader ()
{
  // for now, only option is the circuit loader:
  cktLoaderPtr_ = new N_LOA_CktLoader ();
  loaderPtr_ = cktLoaderPtr_;

  return cktLoaderPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_LoaderMgr::createNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/18/05
//-----------------------------------------------------------------------------
N_LOA_NonlinearEquationLoader * N_LOA_LoaderMgr::createNonlinearEquationLoader ()
{
  // for now, only option is the circuit loader:
  nonlinearEquationLoaderPtr_ = new N_LOA_NonlinearEquationLoader ();
  loaderPtr_ = nonlinearEquationLoaderPtr_;

  return nonlinearEquationLoaderPtr_;
}

