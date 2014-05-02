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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys_IC.h,v $
//
// Purpose        : Concrete class for augmenting the Jacobian for
//                  .IC simulations (initial condition)
//
// Special Notes  :
//
// Creator        : Eric R. Keiter
//
// Creation Date  : 09/15/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_IC_h
#define Xyce_N_NLS_NOX_AugmentLinSys_IC_h

#include "Teuchos_RefCountPtr.hpp"
#include "N_NLS_NOX_AugmentLinSys.h"
#include "N_PDS_ParMap.h"

#include <N_UTL_fwd.h>

#include <set>

class Epetra_MapColoring;

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::AugmentLinSysIC
// Purpose       : Handles matrix augmentation to support .IC statements.
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems modeling.
// Creation Date : 09/15/2007
//-----------------------------------------------------------------------------
namespace N_NLS_NOX {

class AugmentLinSysIC : public N_NLS_NOX::AugmentLinSys {

public:
  //! Ctor.
  AugmentLinSysIC(Xyce::NodeNamePairMap & op_in,
      const Teuchos::RefCountPtr <Epetra_MapColoring>& color_map,
      N_LAS_Vector* cloneVector);

  //! Dtor.
  ~AugmentLinSysIC();

  void setProgressVariable(double dummy) {return;}

  void augmentResidual(const N_LAS_Vector * solution,
		       N_LAS_Vector * residual_vector);

  void augmentJacobian(N_LAS_Matrix * jacobian);

 private:

  //! map of specified variables
  Xyce::NodeNamePairMap & op_;

  //! Color 0 are the voltage unknowns.
  Teuchos::RefCountPtr<Epetra_MapColoring> color_map_;

  //! Temporary vector used to store diagonal.
  N_LAS_Vector* tmp_vector_ptr_;

};

}

#endif

