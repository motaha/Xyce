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
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys_GStepping.h,v $
//
// Purpose        : Concrete class for augmenting the Jacobian for
//                  pseudo-transient solves.
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
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_GStepping_h
#define Xyce_N_NLS_NOX_AugmentLinSys_GStepping_h


#include <vector>
#include "Teuchos_RefCountPtr.hpp"
#include "N_NLS_NOX_AugmentLinSys.h"          

class Epetra_MapColoring;

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::GStepping
// Purpose       : Handles matrix augmentation to support GMIN stepping.
// Creator       : Roger Pawlowski, SNL, 9233
// Creation Date : 3/6/2006
//-----------------------------------------------------------------------------
namespace N_NLS_NOX {
  
  class GStepping : public N_NLS_NOX::AugmentLinSys {
    
  public:

    enum NodeListType {
      NLT_VoltageNodes,
      NLT_AllVoltageUnknowns
    };

  public:
    
    //! Ctor for the voltage nodes as a GID list.
    GStepping(const std::vector<int>& vnodeGIDVec,
	      N_LAS_Vector* cloneVector,
	      double endValue,
        double residCond=0);
    
    //! Ctor with an EpetraMapColoring for all voltage unknowns.
    GStepping(const Teuchos::RefCountPtr
	      <Epetra_MapColoring>& color_map,
	      N_LAS_Vector* cloneVector,
	      double endValue,
        double residCond=0);
    
    //! Dtor.
    ~GStepping();
    
    void setProgressVariable(double time_step_size);
    
    void augmentResidual(const N_LAS_Vector * solution,
			 N_LAS_Vector * residual_vector);
    
    void augmentJacobian(N_LAS_Matrix * jacobian);

    inline void setResidualConductance(double c) {residualConductance_=c;};

  private:

    //! Type of list we are using.
    NodeListType node_list_type_;

    //! Conductance.
    double conductance_;
    
    //! List of voltage node GIDs.
    const std::vector<int> vnodeGIDVec_;
    
    //! Color 0 are the voltage unknowns.
    Teuchos::RefCountPtr<Epetra_MapColoring> color_map_;
  
    //! Temporary vector used to store diagonal.
    N_LAS_Vector* tmp_vector_ptr_;
    
    //! low end of the exponential term.
    double scaled_end_value_;

    //! residual value of the conductance.  Should almost always be zero
    double residualConductance_;

  };
  
} 

#endif 

