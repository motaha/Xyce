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
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys_GStepping.C,v $
//
// Purpose        : Algorithm for augmenting the Jacobian for pseudo
//                  transient solves using vnode conductance.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 03/07/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "Epetra_MapColoring.h"
#include "N_NLS_NOX_AugmentLinSys_GStepping.h"


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::GStepping::GStepping
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
N_NLS_NOX::GStepping::GStepping(
        const std::vector<int>& vnodeGIDVec,
				N_LAS_Vector* cloneVector,
				double scaledEndValue,
        double resCond) :
  node_list_type_(NLT_VoltageNodes),
  vnodeGIDVec_(vnodeGIDVec),
  tmp_vector_ptr_(0),
  scaled_end_value_(scaledEndValue),
  residualConductance_(resCond)
{
  tmp_vector_ptr_ = new N_LAS_Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::GStepping::GStepping
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
N_NLS_NOX::GStepping::
GStepping(const Teuchos::RefCountPtr<Epetra_MapColoring>& color_map,
	  N_LAS_Vector* cloneVector,
	  double scaledEndValue,
    double resCond) :
  node_list_type_(NLT_AllVoltageUnknowns),
  scaled_end_value_(scaledEndValue),
  residualConductance_(resCond)
{
  color_map_ = color_map;
  tmp_vector_ptr_ = new N_LAS_Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::GStepping::~GStepping
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
N_NLS_NOX::GStepping::~GStepping()
{
  delete tmp_vector_ptr_;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::GStepping::setProgressVariable
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
void N_NLS_NOX::GStepping::setProgressVariable(double conductance)
{
  // Direct continuation of conductance (con param goes from 1.0e4 -> 0.0
  //conductance_ = conductance;
  
  // Exponential Continuation (con param goes from +4 -> -log10(endValue))
  conductance_ = pow(10.0, conductance) - pow(10.0, scaled_end_value_) + residualConductance_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::GStepping::augmentResidual
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
void N_NLS_NOX::GStepping::augmentResidual(const N_LAS_Vector * solution,
					   N_LAS_Vector * residualVector)
{
  if (node_list_type_ == NLT_VoltageNodes) 
  {
    std::vector<int>::const_iterator i = vnodeGIDVec_.begin();
    std::vector<int>::const_iterator stop = vnodeGIDVec_.end();
    for ( ; i < stop; ++i) 
    {
      double value = conductance_ * 
        (const_cast<N_LAS_Vector*>(solution))->getElementByGlobalIndex(*i);

      residualVector->sumElementByGlobalIndex(*i, value);
    }
  }
  else 
  {
    for (std::size_t i = 0; i <  tmp_vector_ptr_->localLength(); ++i) 
    {
      if ( (*color_map_)[i] == 0)
      {
        (*residualVector)[i] += conductance_ * (const_cast<N_LAS_Vector&>(*solution))[i]; 
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::GStepping::augmentJacobian
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date : 
//-----------------------------------------------------------------------------
void N_NLS_NOX::GStepping::augmentJacobian(N_LAS_Matrix * jacobian)
{
  jacobian->getDiagonal(*tmp_vector_ptr_);
  
  if (node_list_type_ == NLT_VoltageNodes) 
  {
    std::vector<int>::const_iterator i = vnodeGIDVec_.begin();
    std::vector<int>::const_iterator stop = vnodeGIDVec_.end();
    for ( ; i < stop; ++i) 
    {
      tmp_vector_ptr_->sumElementByGlobalIndex(*i, conductance_);
    }
  }
  else 
  {
    for (std::size_t i = 0; i <  tmp_vector_ptr_->localLength(); ++i) 
    {
      if ( (*color_map_)[i] == 0)
      {
        (*tmp_vector_ptr_)[i] += conductance_; 
      }
    }
  }

  jacobian->replaceDiagonal(*tmp_vector_ptr_);
}

