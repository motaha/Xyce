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
// Filename       : $RCSfile: N_NLS_NOX_AugmentLinSys_IC.C,v $
//
// Purpose        : Algorithm for augmenting the Jacobian for .IC
//                  operating point solves.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 09/15/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.16 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_fwd.h>
#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------

#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "Epetra_MapColoring.h"
#include "N_ERH_ErrorMgr.h"
#include "N_NLS_NOX_AugmentLinSys_IC.h"

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC::AugmentLinSysIC
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/15/07
//-----------------------------------------------------------------------------
N_NLS_NOX::AugmentLinSysIC::AugmentLinSysIC
  (Xyce::NodeNamePairMap & op_in,
  const Teuchos::RefCountPtr <Epetra_MapColoring>& color_map,
  N_LAS_Vector* cloneVector
  )
  : op_       (op_in),
    tmp_vector_ptr_(0)
{
  color_map_ = color_map;
  tmp_vector_ptr_ = new N_LAS_Vector(*cloneVector);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC::~AugmentLinSysIC
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/15/07
//-----------------------------------------------------------------------------
N_NLS_NOX::AugmentLinSysIC::~AugmentLinSysIC()
{
  delete tmp_vector_ptr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC::augmentResidual
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/15/07
//-----------------------------------------------------------------------------
void N_NLS_NOX::AugmentLinSysIC::augmentResidual
   (const N_LAS_Vector * solution, N_LAS_Vector * residual_vector)
{
#ifdef Xyce_DEBUG_IC
  cout << "Inside AugmentLinSysIC::augmentResidual:"  << endl;
#endif

  Xyce::NodeNamePairMap::iterator op_i = op_.begin();
  Xyce::NodeNamePairMap::iterator op_end = op_.end();
  for ( ; op_i != op_end ; ++op_i)
  {
    int row = (*op_i).second.first;
    int global_row(row);

    if ( (*color_map_)[row] == 0)
    {
      (*residual_vector)[row]  = 0.0;
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::AugmentLinSysIC::augmentJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/15/07
//-----------------------------------------------------------------------------
void N_NLS_NOX::AugmentLinSysIC::augmentJacobian(N_LAS_Matrix * jacobian)
{
#ifdef Xyce_DEBUG_IC
  cout << "Inside AugmentLinSysIC::augmentJacobian:"  << endl;
#endif

  std::vector<int> col;
  std::vector<double> val;
  Xyce::NodeNamePairMap::iterator op_i = op_.begin();
  Xyce::NodeNamePairMap::iterator op_end = op_.end();

  jacobian->getDiagonal(*tmp_vector_ptr_);

  for ( ; op_i != op_end ; ++op_i)
  {
    int row = (*op_i).second.first;
    int rowLen(0);
    int numEntries(0);

    if ( (*color_map_)[row] == 0)
    {
      rowLen = jacobian->getLocalRowLength(row);

      col.resize(rowLen,0);
      val.resize(rowLen,0.0);
      jacobian->getLocalRowCopy(row, rowLen, numEntries, &val[0], &col[0]);

      // zero out the entire row.
      for (int i=0;i<val.size();++i) val[i] = 0.0;

      jacobian->putLocalRow(row, rowLen, &val[0], &col[0]);

      // set the diagonal to 1.0.
      (*tmp_vector_ptr_)[row] = 1.0;
    }
  }

  jacobian->replaceDiagonal(*tmp_vector_ptr_);

  return;
}

