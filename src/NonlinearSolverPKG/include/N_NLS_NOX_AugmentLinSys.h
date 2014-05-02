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
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys.h,v $
//
// Purpose        : Strategy pattern for allowing algorithms to augment
//                  the Jacobian and residual.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, NLS, 9233
//
// Creation Date  : 3/6/2006
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_h
#define Xyce_N_NLS_NOX_AugmentLinSys_h

//-----------------------------------------------------------------------------
// Class         : N_NLS_NOX::AugmentLinSys
// Purpose       :
// Creator       : Roger Pawlowski, SNL, 9233
// Creation Date : 3/6/2006
//-----------------------------------------------------------------------------

/*! \brief Pure virtual class to augment a linear system.

    Strategy pattern for allowing algorithms to augment the Jacobian
    and residual.  This class is used to provide an algorithm for
    augmenting the linear system for various solution techniques.
    Homotopy, pseudo-transient, and "gmin" stepping are examples.  In
    each case, the Jacobian's diagonal is changed/added to.  Some
    algorithms also augment the residual.
  
*/

// Forward Decalarations
class N_LAS_Vector;
class N_LAS_Matrix;

namespace N_NLS_NOX {

class AugmentLinSys {

public:

  //! Ctor.
  AugmentLinSys(){};

  //! Dtor.
  ~AugmentLinSys() {};

  //! Set the progress variable (time step size for pseudo transient).
  virtual void setProgressVariable(double value) = 0;

  //! Augments the Residual.
  virtual void augmentResidual(const N_LAS_Vector * solution,
			       N_LAS_Vector * residual_vector) = 0;
  
  //! Augments the Jacobian.
  virtual void augmentJacobian(N_LAS_Matrix * jacobian) = 0;
  
};
 
} 

#endif 

