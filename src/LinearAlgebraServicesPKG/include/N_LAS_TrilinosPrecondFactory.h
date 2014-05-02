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
// Filename       : $RCSfile: N_LAS_TrilinosPrecondFactory.h,v $
//
// Purpose        : Preconditioner Factory for Trilinos
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

#ifndef Xyce_N_LAS_TrilinosPrecondFactory_h
#define Xyce_N_LAS_TrilinosPrecondFactory_h

// ---------- Standard Includes ----------

#include <string>
#include <N_UTL_fwd.h>

// ----------   Xyce Includes   ----------

#include <N_LAS_PrecondFactory.h>
#include <N_ERH_ErrorMgr.h>

// ----------  Fwd Declares  -------------

class N_LAS_Preconditioner;
class N_LAS_Problem;
class N_LAS_System;

//-----------------------------------------------------------------------------
// Class         : N_LAS_PrecondFactory
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
class N_LAS_TrilinosPrecondFactory : public N_LAS_PrecondFactory
{
public:
  // Basic Constructor, sets preconditioner factory options.
  N_LAS_TrilinosPrecondFactory( const N_UTL_OptionBlock & OB );

  // Destructor
  virtual ~N_LAS_TrilinosPrecondFactory() {}

  // Creates a new preconditioner (matrix based).
  Teuchos::RCP<N_LAS_Preconditioner> create( const Teuchos::RCP<N_LAS_Problem> & problem );

  // Creates a new preconditioner (matrix free)
  // NOTE:  This creation type is not supported by the Trilinos factory at this time.
  Teuchos::RCP<N_LAS_Preconditioner> create( const Teuchos::RCP<N_LAS_System> & lasSystem )
  {
    Xyce::Report::DevelFatal0().in("N_LAS_TrilinosPrecondFactory::create()") << "using N_LAS_System is not supported!";
    return Teuchos::null;
  }

private:

  std::string precType_;
  Teuchos::RCP<const N_UTL_OptionBlock> OB_;

  // Default constructor and copy constructor.
  N_LAS_TrilinosPrecondFactory();
  N_LAS_TrilinosPrecondFactory( const N_LAS_TrilinosPrecondFactory& pf );
};


#endif

