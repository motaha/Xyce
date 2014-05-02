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
// Filename       : $RCSfile: N_LAS_Util.h,v $
//
// Purpose        : Linear algebra utilities.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Computational Sciences
//
// Creation Date  : 8/9/12
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Util_h
#define Xyce_N_LAS_Util_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

// ----------  Other Includes   ----------

// --------  Forward Declarations --------

//-----------------------------------------------------------------------------
// Class         : N_LAS_Util
// Purpose       : Provides some basic utilities for linear algebra computations.
// Special Notes :
// Creator       : Heidi K. Thornquist, SNL, Computational Sciences
// Creation Date : 8/9/12
//-----------------------------------------------------------------------------
namespace N_LAS_Util
{

  // y = alpha*A*x + beta*y
  void crsAxpy( int Nrows, double alpha, double* Aval, int* ArowPtr, int* AcolInd, double* xval, double beta, double* yval )
  {
    int nnz = ArowPtr[Nrows];

    // y *= beta
    for (int i=0; i<Nrows; i++)
      yval[i] *= beta;

    for (int i=0; i<Nrows; i++)
    {
      double sum=0.0;
      for (int j=ArowPtr[i]; j<ArowPtr[i+1]; j++)
      {
        sum += Aval[j]*xval[AcolInd[j]];
      }
      yval[i] += alpha*sum;
    }
  }

}

#endif

