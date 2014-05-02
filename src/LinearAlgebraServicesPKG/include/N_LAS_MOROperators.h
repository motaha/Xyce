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
// Filename       : $RCSfile: N_LAS_MOROperators.h,v $
//
// Purpose        : Specification file for operators that are necessary for 
//                  performing model-order reduction.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL
//
// Creation Date  : 06/04/12
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_MOROperators_h
#define Xyce_N_LAS_MOROperators_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

// ----------  Other Includes   ----------

// Include header for Epetra compressed-row storage matrix and linear problem
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h" 

// Include header for Amesos solver and solver interface for Epetra_Operator
#include "Amesos.h"

// Include selected communicator class and map required by Epetra objects
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//-----------------------------------------------------------------------------
// Class         : N_LAS_AmesosGenOp
// Purpose       : Implementation of the Epetra_Operator interface to define
//                 inv(A)*B, where inv(A) is computed by Amesos.
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 06/04/12
//-----------------------------------------------------------------------------
class N_LAS_AmesosGenOp : public virtual Epetra_Operator
{
public:
  // Basic constructor
  N_LAS_AmesosGenOp( const Teuchos::RCP<Amesos_BaseSolver>& solver,
                     const Teuchos::RCP<Epetra_Operator>& B,
                     bool useTranspose = false );
  // Destructor
  ~N_LAS_AmesosGenOp() {};

  // Methods for supporting Epetra_Operator interface
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;
  const char* Label() const { return "Amesos direct solver for applying A^{-1}B"; }
  bool UseTranspose() const { return useTranspose_; }
  int SetUseTranspose( bool useTranspose ) { useTranspose_ = useTranspose; return 0; }
  const Epetra_Comm& Comm() const { return solver_->Comm(); };
  const Epetra_Map& OperatorDomainMap() const { return B_->OperatorDomainMap(); }
  const Epetra_Map& OperatorRangeMap() const { return B_->OperatorRangeMap(); }

  // Epetra_Operator interface methods that are not supported.
  // Note:  ApplyInverse not defined because M not guaranteed to have an inverse.
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const { return -1; }
  bool HasNormInf() const { return false; };
  double NormInf() const { return -1.0; };
   
private:  
  // Default constructor
  N_LAS_AmesosGenOp () {};
  
  // Copy constructor 
  N_LAS_AmesosGenOp ( const N_LAS_AmesosGenOp& genOp ) {};
  
  // Epetra_LinearProblem contained in the Amesos_BaseSolver
  bool useTranspose_;
  Teuchos::RCP<Amesos_BaseSolver> solver_;
  Teuchos::RCP<Epetra_Operator> B_;
  Teuchos::RCP<Epetra_LinearProblem> problem_;
  
};

#endif


