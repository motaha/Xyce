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
//    MERCHANTABILITY or FITNESS FOR A PARTIC_GminULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys_IC_Gmin.h,v $
//
// Purpose        : Concrete class for augmenting the Jacobian for
//                  .IC simulations (initial condition) with gmin stepping.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter
//
// Creation Date  : 04/29/12
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

#ifndef Xyce_N_NLS_NOX_AugmentLinSys_IC_Gmin_h
#define Xyce_N_NLS_NOX_AugmentLinSys_IC_Gmin_h

#include <vector>
#include <set>

#include "Teuchos_RefCountPtr.hpp"
#include "N_NLS_NOX_AugmentLinSys.h"

#include <N_UTL_fwd.h>

class Epetra_MapColoring;

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::AugmentLinSysIC_Gmin
// Purpose       : Handles matrix augmentation to support .IC statements, if
//                 gmin stepping is also being applied.
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems modeling.
// Creation Date : 04/29/2012
//-----------------------------------------------------------------------------
namespace N_NLS_NOX {

class AugmentLinSysIC_Gmin : public N_NLS_NOX::AugmentLinSys {

  public:
    enum NodeListType {
      NLT_VoltageNodes,
      NLT_AllVoltageUnknowns
    };

  public:
    //! Ctor.
    AugmentLinSysIC_Gmin
    ( Xyce::NodeNamePairMap & op_in,
      const Teuchos::RefCountPtr <Epetra_MapColoring>& ICcolor_map,
      const std::vector<int>& vnodeGIDVec,
      N_LAS_Vector* cloneVector,
			double scaledEndValue,
      double resCond);

    //! Ctor.
    AugmentLinSysIC_Gmin
    ( Xyce::NodeNamePairMap & op_in,
      const Teuchos::RefCountPtr <Epetra_MapColoring>& ICcolor_map,
      const Teuchos::RefCountPtr <Epetra_MapColoring>& GMINcolor_map,
      N_LAS_Vector* cloneVector,
			double scaledEndValue,
      double resCond);

  //! Dtor.
    ~AugmentLinSysIC_Gmin();

    void setProgressVariable(double dummy);

    void augmentResidual(const N_LAS_Vector * solution,
             N_LAS_Vector * residual_vector);

    void augmentJacobian(N_LAS_Matrix * jacobian);

  private:

    //! Type of list we are using.
    NodeListType node_list_type_;

    //! Conductance.
    double conductance_;

    //! low end of the exponential term.
    double scaled_end_value_;

    //! residual value of the conductance.  Should almost always be zero
    double residualConductance_;

    //! List of voltage node GIDs.
    const std::vector<int> vnodeGIDVec_;

    //! map of specified variables
    Xyce::NodeNamePairMap & op_;

    //! Color 0 are the voltage unknowns.
    //! For the IC color map, the voltage nodes attached to
    //! independent voltage sources are not included.
    Teuchos::RefCountPtr<Epetra_MapColoring> ICcolor_map_;
    Teuchos::RefCountPtr<Epetra_MapColoring> GMINcolor_map_;

    //! Temporary vectors used to store diagonal.
    N_LAS_Vector* vecptr1_;
    N_LAS_Vector* vecptr2_;

};

}

#endif

