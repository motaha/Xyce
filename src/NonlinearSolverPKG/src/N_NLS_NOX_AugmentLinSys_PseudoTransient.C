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
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys_PseudoTransient.C,v $
//
// Purpose        : Algorithm for augmenting the Jacobian for pseudo
//                  transient solves.
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
// Revision Number: $Revision: 1.8 $
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
#include "N_NLS_NOX_AugmentLinSys_PseudoTransient.h"


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysPseudoTransient::AugmentLinSysPseudoTransient
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
N_NLS_NOX::AugmentLinSysPseudoTransient::AugmentLinSysPseudoTransient
    (const Teuchos::RefCountPtr<Epetra_MapColoring>&
			     color_map,
			     N_LAS_Vector* cloneVector,
			     bool useVoltageScaleFactor,
			     double voltageScaleFactor)
{
  use_voltage_scale_factor_ = useVoltageScaleFactor;
  voltage_scale_factor_ = voltageScaleFactor;
  color_map_ = color_map;
  tmp_vector_ptr_ = new N_LAS_Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysPseudoTransient::~AugmentLinSysPseudoTransient
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
N_NLS_NOX::AugmentLinSysPseudoTransient::~AugmentLinSysPseudoTransient()
{
  delete tmp_vector_ptr_;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysPseudoTransient::setProgressVariable
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_NOX::AugmentLinSysPseudoTransient::setProgressVariable
  (double time_step_size)
{
  time_step_size_ = time_step_size;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysPseudoTransient::augmentResidual
// Purpose       :
// Special Notes : no-op for pseudo-transient.
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_NOX::AugmentLinSysPseudoTransient::augmentResidual
  (const N_LAS_Vector * solution, N_LAS_Vector * residual_vector)
{
  // Nothing to do for pseudo transient!!
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysPseudoTransient::augmentJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_NOX::AugmentLinSysPseudoTransient::augmentJacobian
  (N_LAS_Matrix * jacobian)
{
  //cout << "Augmenting Jacobian for Pseudo Transient" << endl;
  //cout << "Pseudo Trans Step Size = " << pseudoTransientTimeStep_ << endl;
  
  //color_map_->Print(cout);
    
  //jacobian->printPetraObject();

  //jacobian->scale(conParamValue);

  jacobian->getDiagonal(*tmp_vector_ptr_);
  
  //tmp_vector_ptr_->printPetraObject();
    
  double value = 1.0 / time_step_size_;
    
  //cout << "Pseudo Transient Time Step Size = " << time_step_size_ << endl;

  if (!use_voltage_scale_factor_) 
  {
    tmp_vector_ptr_->addScalar(value);
  }
  else 
  {
    for (std::size_t i = 0; i <  tmp_vector_ptr_->localLength(); ++i) 
    {
      if ( (*color_map_)[i] == 0)
      {
        (*tmp_vector_ptr_)[i] += value * voltage_scale_factor_; 
      }
      else
      {
        (*tmp_vector_ptr_)[i] += value;
      }
    }
    //RPP Might need to export local values for tmp_vector_ptr_ here
    //for parallel.
  }

  jacobian->replaceDiagonal(*tmp_vector_ptr_);

  //jacobian->printPetraObject();
   
  //cout << "Time step size = " << time_step_size_ << endl;

  //if (use_voltage_scale_factor_)
  //cout << "Voltage scale factor = " << voltage_scale_factor_ << endl;
}

