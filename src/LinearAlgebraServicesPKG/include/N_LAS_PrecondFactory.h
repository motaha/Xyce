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
// Filename       : $RCSfile: N_LAS_PrecondFactory.h,v $
//
// Purpose        : Preconditioner Factory
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
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_PrecondFactory_h
#define Xyce_N_LAS_PrecondFactory_h

// ---------- Standard Includes ----------

#include <string>

#include <N_UTL_fwd.h>

class N_LAS_Preconditioner;
class N_LAS_Problem;
class N_LAS_System;

#include <Teuchos_RCP.hpp>

//-----------------------------------------------------------------------------
// Class         : N_LAS_PrecondFactory
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
class N_LAS_PrecondFactory
{
public:
  // Default Constructor
  N_LAS_PrecondFactory() {}

  // Basic Constructor, sets preconditioner factory options.
  N_LAS_PrecondFactory( const N_UTL_OptionBlock & OB ) {}

  // Destructor
  virtual ~N_LAS_PrecondFactory() {}

  // Creates a new preconditioner (matrix based).
  virtual Teuchos::RCP<N_LAS_Preconditioner> create( const Teuchos::RCP<N_LAS_Problem> & problem ) = 0;

  // Creates a new preconditioner (matrix free).
  virtual Teuchos::RCP<N_LAS_Preconditioner> create( const Teuchos::RCP<N_LAS_System> & lasSystem ) = 0;
};

#endif

