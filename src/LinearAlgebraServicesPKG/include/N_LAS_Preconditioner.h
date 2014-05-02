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
// Filename       : $RCSfile: N_LAS_Preconditioner.h,v $
//
// Purpose        : interface to preconditioner
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 09/27/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Preconditioner_h
#define Xyce_N_LAS_Preconditioner_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_UTL_fwd.h>

// ---------- Trilinos Includes ----------

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>

// ----------  Fwd Declares     ----------

class N_LAS_Problem;
class N_LAS_MultiVector;

class Epetra_Operator;

//-----------------------------------------------------------------------------
// Class         : N_LAS_Preconditioner
// Purpose       : interface to preconditioner
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
class N_LAS_Preconditioner : public Teuchos::Describable
{

public:
  // Constructors
  N_LAS_Preconditioner() {}

  // Destructor
  virtual ~N_LAS_Preconditioner() {}

  // Set the preconditioner options
  virtual bool setOptions(const N_UTL_OptionBlock & OB) = 0;
  virtual bool setDefaultOptions() = 0;
  virtual bool setDefaultOption( const std::string & option ) = 0;

  // Set individual preconditioner options
  virtual bool setParam( const N_UTL_Param & param ) = 0;

  // Set the matrix pattern for the preconditioner
  virtual bool initGraph( const Teuchos::RCP<N_LAS_Problem> & problem ) = 0;

  // Set the matrix values for the preconditioner
  virtual bool initValues( const Teuchos::RCP<N_LAS_Problem> & problem ) = 0;

  // Compute the preconditioner using the current matrix values.
  virtual bool compute() = 0;

  // Apply the preconditioner; y = M*x.
  virtual int apply( N_LAS_MultiVector & x, N_LAS_MultiVector & y ) = 0;

  // Return the preconditioner as an Epetra_Operator object.
  virtual Teuchos::RCP<Epetra_Operator> epetraObj() = 0;

private:

  // No copying
  N_LAS_Preconditioner(const N_LAS_Preconditioner & right);
  N_LAS_Preconditioner & operator=(const N_LAS_Preconditioner & right);

  // No comparison
  bool operator==(const N_LAS_Preconditioner & right) const;
  bool operator!=(const N_LAS_Preconditioner & right) const;

};

#endif
