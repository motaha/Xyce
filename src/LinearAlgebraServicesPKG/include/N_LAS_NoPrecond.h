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
// Filename       : $RCSfile: N_LAS_NoPrecond.h,v $
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
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_NoPrecond_h
#define Xyce_N_LAS_NoPrecond_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_UTL_fwd.h>
#include <N_LAS_Preconditioner.h>

// ---------- Trilinos Includes ----------

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>

// ----------  Fwd Declares     ----------

class N_LAS_Problem;
class N_LAS_MultiVector;

class Epetra_Operator;

//-----------------------------------------------------------------------------
// Class         : N_LAS_NoPrecond
// Purpose       : interface to preconditioner
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
class N_LAS_NoPrecond : public N_LAS_Preconditioner
{

public:
  // Constructors
  N_LAS_NoPrecond() {}

  // Destructor
  virtual ~N_LAS_NoPrecond() {}

  // Set the preconditioner options
  virtual bool setOptions(const N_UTL_OptionBlock & OB) { return true; }
  virtual bool setDefaultOptions() { return true; }
  virtual bool setDefaultOption( const std::string & option ) { return true; }

  // Set individual preconditioner options
  virtual bool setParam( const N_UTL_Param & param ) { return true; }

  // Set the matrix pattern for the preconditioner
  virtual bool initGraph( const Teuchos::RCP<N_LAS_Problem> & problem ) { return true; }

  // Set the matrix values for the preconditioner
  virtual bool initValues( const Teuchos::RCP<N_LAS_Problem> & problem ) { return true; }

  // Compute the preconditioner using the current matrix values.
  virtual bool compute() { return true; }

  // Apply the preconditioner; y = M*x.
  virtual int apply( N_LAS_MultiVector & x, N_LAS_MultiVector & y ) { return 0; }

  // Return the preconditioner as an Epetra_Operator object.
  virtual Teuchos::RCP<Epetra_Operator> epetraObj() { return Teuchos::null; }

private:

  // No copying
  N_LAS_NoPrecond(const N_LAS_NoPrecond & right);
  N_LAS_NoPrecond & operator=(const N_LAS_NoPrecond & right);

  // No comparison
  bool operator==(const N_LAS_NoPrecond & right) const;
  bool operator!=(const N_LAS_NoPrecond & right) const;

};

#endif
