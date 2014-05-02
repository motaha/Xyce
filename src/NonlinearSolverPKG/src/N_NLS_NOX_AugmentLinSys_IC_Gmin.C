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
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys_IC_Gmin.C,v $
//
// Purpose        : Algorithm for augmenting the Jacobian for .IC_Gmin
//                  operating point solves.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 4/29/2012
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------

#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "Epetra_MapColoring.h"
#include "N_ERH_ErrorMgr.h"
#include "N_NLS_NOX_AugmentLinSys_IC_Gmin.h"

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC_Gmin::AugmentLinSysIC_Gmin
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
N_NLS_NOX::AugmentLinSysIC_Gmin::AugmentLinSysIC_Gmin(
  Xyce::NodeNamePairMap & op_in,
        const Teuchos::RefCountPtr <Epetra_MapColoring>& ICcolor_map,
        const std::vector<int>& vnodeGIDVec,
        N_LAS_Vector* cloneVector,
				double scaledEndValue,
        double resCond) :
    op_       (op_in),
    vecptr1_(0),
    vecptr2_(0),
    node_list_type_(NLT_VoltageNodes),
    vnodeGIDVec_(vnodeGIDVec),
    scaled_end_value_(scaledEndValue),
    residualConductance_(resCond)
{
  ICcolor_map_ = ICcolor_map;
  vecptr1_ = new N_LAS_Vector(*cloneVector);
  vecptr2_ = new N_LAS_Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC_Gmin::AugmentLinSysIC_Gmin
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
N_NLS_NOX::AugmentLinSysIC_Gmin::AugmentLinSysIC_Gmin(
  Xyce::NodeNamePairMap & op_in,
        const Teuchos::RefCountPtr <Epetra_MapColoring>& ICcolor_map,
        const Teuchos::RefCountPtr <Epetra_MapColoring>& GMINcolor_map,
        N_LAS_Vector* cloneVector,
				double scaledEndValue,
        double resCond) :
    op_       (op_in),
    vecptr1_(0),
    vecptr2_(0),
    node_list_type_(NLT_VoltageNodes),
    scaled_end_value_(scaledEndValue),
    residualConductance_(resCond)
{
  ICcolor_map_ = ICcolor_map;
  GMINcolor_map_ = GMINcolor_map;
  vecptr1_ = new N_LAS_Vector(*cloneVector);
  vecptr2_ = new N_LAS_Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC_Gmin::~AugmentLinSysIC_Gmin
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
N_NLS_NOX::AugmentLinSysIC_Gmin::~AugmentLinSysIC_Gmin()
{
  delete vecptr1_;
  delete vecptr2_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::IC_Gmin::setProgressVariable
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_NOX::AugmentLinSysIC_Gmin::setProgressVariable(double conductance)
{
  // Direct continuation of conductance (con param goes from 1.0e4 -> 0.0
  //conductance_ = conductance;

  // Exponential Continuation (con param goes from +4 -> -log10(endValue))
  conductance_ = pow(10.0, conductance) - pow(10.0, scaled_end_value_) + residualConductance_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC_Gmin::augmentResidual
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
void N_NLS_NOX::AugmentLinSysIC_Gmin::augmentResidual
   (const N_LAS_Vector * solution, N_LAS_Vector * residual_vector)
{
#ifdef Xyce_DEBUG_IC_Gmin
  cout << "Inside AugmentLinSysIC_Gmin::augmentResidual:"  << endl;
#endif

  // GMIN portion
   if (node_list_type_ == NLT_VoltageNodes)
  {
    std::vector<int>::const_iterator i = vnodeGIDVec_.begin();
    std::vector<int>::const_iterator stop = vnodeGIDVec_.end();
    for ( ; i < stop; ++i)
    {
      double value = conductance_ *
        (const_cast<N_LAS_Vector*>(solution))->getElementByGlobalIndex(*i);

      residual_vector->sumElementByGlobalIndex(*i, value);
    }
  }
  else
  {
    for (std::size_t i = 0; i <  vecptr1_->localLength(); ++i)
    {
      if ( (*GMINcolor_map_)[i] == 0)
      {
        (*residual_vector)[i] += conductance_ * (const_cast<N_LAS_Vector&>(*solution))[i];
      }
    }
  }

  // IC portion
  Xyce::NodeNamePairMap::iterator op_i = op_.begin();
  Xyce::NodeNamePairMap::iterator op_end = op_.end();
  for ( ; op_i != op_end ; ++op_i)
  {
    int row = (*op_i).second.first;
    int global_row(row);

    if ( (*ICcolor_map_)[row] == 0)
    {
      (*residual_vector)[row]  = 0.0;
    }
  }



  return;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC_Gmin::augmentJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/29/2012
//-----------------------------------------------------------------------------
void N_NLS_NOX::AugmentLinSysIC_Gmin::augmentJacobian(N_LAS_Matrix * jacobian)
{
#ifdef Xyce_DEBUG_IC_Gmin
  cout << "Inside AugmentLinSysIC_Gmin::augmentJacobian:"  << endl;
#endif

  // GMIN portion
  jacobian->getDiagonal(*vecptr1_);
  if (node_list_type_ == NLT_VoltageNodes)
  {
    std::vector<int>::const_iterator i = vnodeGIDVec_.begin();
    std::vector<int>::const_iterator stop = vnodeGIDVec_.end();
    for ( ; i < stop; ++i)
    {
      vecptr1_->sumElementByGlobalIndex(*i, conductance_);
    }
  }
  else
  {
    for (std::size_t i = 0; i <  vecptr1_->localLength(); ++i)
    {
      if ( (*GMINcolor_map_)[i] == 0)
      {
        (*vecptr1_)[i] += conductance_;
      }
    }
  }
  jacobian->replaceDiagonal(*vecptr1_);

  // IC portion
  std::vector<int> col;
  std::vector<double> val;
  Xyce::NodeNamePairMap::iterator op_i = op_.begin();
  Xyce::NodeNamePairMap::iterator op_end = op_.end();

  jacobian->getDiagonal(*vecptr2_);

  for ( ; op_i != op_end ; ++op_i)
  {
    int row = (*op_i).second.first;
    int rowLen(0);
    int numEntries(0);

    if ( (*ICcolor_map_)[row] == 0)
    {
      rowLen = jacobian->getLocalRowLength(row);

      col.resize(rowLen,0);
      val.resize(rowLen,0.0);
      jacobian->getLocalRowCopy(row, rowLen, numEntries, &val[0], &col[0]);

      // zero out the entire row.
      for (int i=0;i<val.size();++i) val[i] = 0.0;

      jacobian->putLocalRow(row, rowLen, &val[0], &col[0]);

      // set the diagonal to 1.0.
      (*vecptr2_)[row] = 1.0;
    }
  }

  jacobian->replaceDiagonal(*vecptr2_);

  return;
}

