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
//
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DiodePDE_DAE.C,v $
//
// Purpose        : One dimensional PDE device, new-DAE functions.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/13/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.29 $
//
// Revision Date  : $Date: 2014/02/24 23:49:18 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#include <N_UTL_Misc.h>
#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_DiodePDE.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

// default number of mesh points:
#define NUM_MESH_POINTS 11

// default maximum number of nonzero entries in a matrix row
#define MAX_COLS_PER_ROW 10

namespace Xyce {
namespace Device {
namespace DiodePDE {

//-----------------------------------------------------------------------------
// Function      : N_DEV_Instance::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;
  bool bs1;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "loadDAEFVector:  doubleDCOPStep="<<getSolverState().doubleDCOPStep;
    if (getSolverState().dcopFlag) Xyce::dout() << "  DCOP load" << std::endl;
    else                   Xyce::dout() << "  Transient load" << std::endl;
  }
#endif

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep == 0)
  {
    equationSet = 0;
    bs1 = loadDAEFNonlinPoisson ();
  }
  else
  {
    equationSet = 1;

    if (getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM ||
        getSolverState().twoLevelNewtonCouplingMode==FULL_PROBLEM)
    {
      bs1 = loadDAEFDDFormulation ();
    }
    else if (getSolverState().twoLevelNewtonCouplingMode==OUTER_PROBLEM)
    {
      bs1 = loadDAEFExtractedConductance ();
    }
    else
    {
      std::string msg = "Instance::loadDAEFVector."
                   "Invalid coupling Mode";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
  }

  return (bsuccess && bs1);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFNonlinPoisson
// Purpose       : Loads the F-vector the nonlinear poisson calculation.
//
// Special Notes : This should be identical to the original loadRHS function,
//                 for the nonlinear poisson.  All of that goes into the
//                 FVector, but with the opposite sign.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFNonlinPoisson ()
{
  double * fvec = extData.daeFVectorRawPtr;

  bool bsuccess = loadVecNLPoisson ( fvec );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFDDFormulation
// Purpose       : This function should be called from the loadDAEFVector
//                 function when solving the drift-diffusion equations.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFDDFormulation ()
{
  double * fvec = extData.daeFVectorRawPtr;

  bool bsuccess = loadVecDDForm ( fvec );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFExtractedConductance ()
{
  double coef;

  N_LAS_Vector & fVec = *(extData.daeFVectorPtr) ;

  // KCL equations for the various connecting terminals:
  int iRow = 0;
  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC, ++iRow)
  {
    coef = bcVec[iBC].currentSum;

    if (iBC==1) coef *= -1.0;

    double voltLimFac = 0.0;
    if (getDeviceOptions().voltageLimiterFlag && voltLimFlag)
    {
      int iCol;
      for (iCol=0; iCol < numElectrodes ; ++iCol)
      {
        if (bcVec[iCol].gid == -1) continue;

        double vdiff = bcVec[iCol].Vckt_final - bcVec[iCol].Vckt_orig;

        vdiff *= scalingVars.V0;

        voltLimFac += vdiff * condVec[iRow][iCol];
      }
    }

    fVec[bcVec[iBC].lid] += -coef + voltLimFac;
  }

  // Load in the zeros...
  for (int i=0;i<NX;++i)
  {
    fVec[ li_Vrowarray[i] ] = 0.0;
    fVec[ li_Nrowarray[i] ] = 0.0;
    fVec[ li_Prowarray[i] ] = 0.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;
  bool bs1;

  N_LAS_Vector & qvec = *(extData.daeQVectorPtr);

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep == 0)
  {
    equationSet = 0;

    bs1 = true; // no-op
  }
  else
  {
    equationSet = 1;

    if (getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM ||
        getSolverState().twoLevelNewtonCouplingMode==FULL_PROBLEM)
    {
      bs1 = loadDAEQDDFormulation ();
    }
    else if (getSolverState().twoLevelNewtonCouplingMode==OUTER_PROBLEM)
    {
      bs1 = loadDAEQExtractedConductance ();
    }
    else
    {
      std::string msg = "Instance::loadDAEQVector."
                   "Invalid coupling Mode";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
  }

  return (bsuccess && bs1);
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQDDFormulation
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQDDFormulation ()
{
  N_LAS_Vector & qvec = *(extData.daeQVectorPtr);
  int i;

  // mesh points for the PDE problem:
  // only do the interior mesh points - boundaries have no dndt terms.
  for (i=1;i<NX-1;++i)
  {
    qvec[ li_Nrowarray[i] ] = -nnVec[i]*scalingVars.t0;
    qvec[ li_Prowarray[i] ] = -npVec[i]*scalingVars.t0;
  } // ip_iter, row loop...

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQExtractedConductance ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx
// Purpose       : This function performs an analytic Jacobian matrix load for
//                 the diode-pde class, for the case of solving a nonlinear
//                 poisson equation.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess;

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep == 0)
  {
    bsuccess = loadDAEdFdxNonlinPoisson ();
  }
  else
  {
    if (getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM ||
        getSolverState().twoLevelNewtonCouplingMode==FULL_PROBLEM)
    {
      bsuccess = loadDAEdFdxDDFormulation ();
    }
    else if (getSolverState().twoLevelNewtonCouplingMode==OUTER_PROBLEM)
    {
      bsuccess = loadDAEdFdxExtractedConductance ();
    }
    else
    {
      Report::DevelFatal().in("Instance::loadDAEdFdx") << "Invalid coupling Mode" << numElectrodes;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxNonlinPoisson
// Purpose       : This function performs an analytic Jacobian matrix load for
//                 the diode-pde class, for the case of solving a nonlinear
//                 poisson equation.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxNonlinPoisson ()
{
  return loadMatNLPoisson( *(extData.dFdxMatrixPtr) );
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxDDFormulation
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxDDFormulation ()
{
  return loadMatDDForm ( *(extData.dFdxMatrixPtr));
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxExtractedConductance ()
{
  bool bsuccess = true;
  bool bs1 = true;

  N_LAS_Matrix & mat = *(extData.dFdxMatrixPtr);

  int iRow, iCol;
  int colsLID[10];
  double valsLID[10];
  int i;
  for (i=0;i<10;++i)
  {
    colsLID[i] = -1;
    valsLID[i] = 0.0;
  }

  // first put 1's on the diagonals of all the mesh-rows:
  for (i=0;i<numMeshPoints;++i)
  {
    valsLID[0] = 1.0;

    iRow = li_Vrowarray[i];
    if (iRow != -1)
    {
      colsLID[0] = iRow;
      bs1 = mat.putLocalRow(iRow,1, valsLID, colsLID);
      bsuccess = bsuccess && bs1;
    }

    iRow = li_Nrowarray[i];
    if (iRow != -1)
    {
      colsLID[0] = iRow;
      bs1 = mat.putLocalRow(iRow,1, valsLID, colsLID);
      bsuccess = bsuccess && bs1;
    }

    iRow = li_Prowarray[i];
    if (iRow != -1)
    {
      colsLID[0] = iRow;
      bs1 = mat.putLocalRow(iRow,1, valsLID, colsLID);
      bsuccess = bsuccess && bs1;
    }
  }

  // now load the equivalent conductances.
  for (iRow=0; iRow < numElectrodes ; ++iRow)
  {

    int iRowLID = bcVec[iRow].gid;

    if (iRowLID == -1) continue;

    int count = 0;

    for (iCol=0; iCol < numElectrodes ; ++iCol)
    {
      if (bcVec[iCol].gid == -1) continue;

      colsLID[count] = bcVec[iCol].gid;
      valsLID[count] = condVec[iRow][iCol];
      ++count;
    }

    bs1 = mat.sumIntoLocalRow (iRowLID, count, valsLID, colsLID);
    bsuccess = bsuccess && bs1;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -3)
  {
    int iE1, iE2;
    char tmpString[128]; for (int i=0;i<128;++i) tmpString[i] = 0;

    Xyce::dout() << std::endl;
    sprintf(tmpString,"ConArray:\t           "); Xyce::dout() << std::string(tmpString);
    for (iE2 = 0; iE2 < numElectrodes; ++ iE2)
    {
      sprintf(tmpString,"\t%14s",bcVec[iE2].eName.c_str()); Xyce::dout() << std::string(tmpString);
    }
    Xyce::dout() << std::endl;

    for (iE1 = 0; iE1 < numElectrodes; ++iE1)
    {
      sprintf(tmpString,"ConArray:\t%14s",bcVec[iE1].eName.c_str()); Xyce::dout() << std::string(tmpString);
      for (iE2 = 0; iE2 < numElectrodes; ++ iE2)
      {
        sprintf(tmpString,"\t%14.4e",condVec[iE1][iE2]); Xyce::dout() << std::string(tmpString);
      }
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess;
  N_LAS_Matrix & dQdxMat = *(extData.dQdxMatrixPtr);

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep == 0)
  {
    bsuccess = true; // no-op
  }
  else
  {
    if (getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM ||
        getSolverState().twoLevelNewtonCouplingMode==FULL_PROBLEM)
    {
      bsuccess = loadDAEdQdxDDFormulation ();
    }
    else if (getSolverState().twoLevelNewtonCouplingMode==OUTER_PROBLEM)
    {
      bsuccess = loadDAEdQdxExtractedConductance ();
    }
    else
    {
      Report::DevelFatal().in("Instance::loadDAEdQdx") << "Invalid coupling Mode " << numElectrodes;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdxDDFormulation
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdxDDFormulation ()
{
  bool bsuccess = true;

  N_LAS_Matrix & dQdxMat = *(extData.dQdxMatrixPtr);

  // load the rows associated with the PDE mesh:
  // use direct matrix access:
  for (int i=1;i<NX-1;++i)
  {
    int li_Nrow = li_Nrowarray[i];
    int li_Prow = li_Prowarray[i];

    // Note::these terms need to be the SAME sign as the terms
    // in loadDAEQVector!

    // electron continuity row:
    // derivative w.r.t. nnVec[i  ]:
    dQdxMat[li_Nrow][li_Ncolarray[i][1]] = -scalingVars.t0;

    // hole continuity row:
    // derivative w.r.t. npVec[i  ]:
    dQdxMat[li_Prow][li_Pcolarray[i][1]] = -scalingVars.t0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdxExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdxExtractedConductance ()
{
  bool bsuccess = true;
  return bsuccess;
}

} // namespace DiodePDE
} // namespace Device
} // namespace Xyce
