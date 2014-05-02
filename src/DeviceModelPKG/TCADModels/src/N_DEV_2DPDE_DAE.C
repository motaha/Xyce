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
// Filename       : $RCSfile: N_DEV_2DPDE_DAE.C,v $
//
// Purpose        : This file contains a lot of the
//                  implementation of the instance class for the two
//                  dimensional PDE based semiconductor device.
//
//                  Functions pertaining to the initial setup are in other
//                  files, as are functions relating to mesh handling and
//                  parameter handling.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/05/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26 $
//
// Revision Date  : $Date: 2014/02/24 23:49:18 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_2DPDE.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_PDE_2DMesh.h>
#include <N_DEV_SourceData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {
namespace TwoDPDE {

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;
  bool bs1;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "loadDAEFVector:  doubleDCOPStep="<<getSolverState().doubleDCOPStep;
    if (getSolverState().dcopFlag) Xyce::dout() << "  DCOP load" <<std::endl;
    else                   Xyce::dout() << "  Transient load" <<std::endl;
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
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFNonlinPoisson ()
{
  return loadVecNLPoisson ( -1.0, extData.daeFVectorPtr );
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFDDFormulation
// Purpose       : This function should be called from the loadDAEFVector
//                 function when solving the drift-diffusion equations.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFDDFormulation ()
{
  // calcTerminalCurrents needs to be here b/c it is otherwise
  // called from updateSecondaryState
  calcTerminalCurrents ();
  return loadVecDDForm (-1.0,0.0, extData.daeFVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFExtractedConductance ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;
  bool bs1;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "loadDAEQVector:  doubleDCOPStep="<<getSolverState().doubleDCOPStep;
    if (getSolverState().dcopFlag) Xyce::dout() << "  DCOP load" <<std::endl;
    else                   Xyce::dout() << "  Transient load" <<std::endl;
  }
#endif

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
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQDDFormulation ()
{
  bool bsuccess = true;
  bool bs1 = true;

  N_LAS_Vector * vecPtr = extData.daeQVectorPtr;

  int i;
  int Nrow, Prow;

  // mesh points for the PDE problem:
  for (i=0;i<numMeshPoints;++i)
  {

    // skip boundary points, for both NEW_BC and old.
    if (boundarySten[i]) continue;

    if( useVectorGIDFlag )
    {
      Nrow = Nrowarray[i];
      Prow = Prowarray[i];

      bs1 = vecPtr->sumElementByGlobalIndex(Nrow, (-nnVec[i]*scalingVars.t0), 0);
      bsuccess = bsuccess && bs1;

      bs1 = vecPtr->sumElementByGlobalIndex(Prow, (-npVec[i]*scalingVars.t0), 0);
      bsuccess = bsuccess && bs1;
    }
    else
    {
      Nrow = li_Nrowarray[i];
      Prow = li_Prowarray[i];
      (*vecPtr)[Nrow] = -nnVec[i]*scalingVars.t0;
      (*vecPtr)[Prow] = -npVec[i]*scalingVars.t0;
    }

  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQExtractedConductance ()
{
  bool bsuccess = true;

  return bsuccess;
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
// Creation Date : 6/21/05
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
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxNonlinPoisson ()
{
  return loadMatNLPoisson( extData.dFdxMatrixPtr );
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxDDFormulation
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxDDFormulation ()
{
  bool bsuccess = true;
  bool bs1 = true;
  // set up some of the partial derivative arrays:
  bs1 = pdRecombination ();   bsuccess = bsuccess && bs1;
  bs1 = pdElectronCurrent (); bsuccess = bsuccess && bs1;
  bs1 = pdHoleCurrent ();     bsuccess = bsuccess && bs1;
  bs1 = pdTerminalCurrents ();bsuccess = bsuccess && bs1;

  if ( !(getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM))
  {
    bs1 = loadMatKCLDDForm ( extData.dFdxMatrixPtr );
    bsuccess = bsuccess && bs1;
  }
  else
  {
    bs1 = loadMatCktTrivial ( extData.dFdxMatrixPtr );
    bsuccess = bsuccess && bs1;
  }

  bs1 = loadMatDDForm ( 0.0, extData.dFdxMatrixPtr );
  bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdxExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdxExtractedConductance ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::loadDAEdFdxExtractedConductance" << std::endl;
  }

#endif

#if 0            // Not yet implemented -- commented out to remove warnings
  int colsGID[10];
  double valsGID[10];
  int i;
  for (i=0;i<10;++i)
  {
    colsGID[i] = -1;
    valsGID[i] = 0.0;
  }

  // first put 1's on the diagonals of all the mesh-rows:
#endif

  // now load the equivalent conductances.

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess;

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
      Report::DevelFatal().in("Instance::loadDAEdQdx") << "Invalid coupling Mode" << numElectrodes;
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
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdxDDFormulation ()
{
  bool bsuccess = true;
  bool bs1 = true;
  int i,j,count;
  int Nrow, Prow;

  int numCol = cols.size();
  if (vals.size () < cols.size()) numCol = vals.size();

  N_LAS_Matrix & QMatrix = (*extData.dQdxMatrixPtr);

  // load the rows associated with the PDE mesh:
  if( useMatrixGIDFlag )
  {
    for (i=0;i<numMeshPoints;++i)
    {
      if (boundarySten[i]) continue;

      Nrow = Nrowarray[i];
      Prow = Prowarray[i];

      // electron continuity row:
      for (j=0;j<numCol;++j) { cols[j] = -1; vals[j] = 0.0; }
      count = 0;

      // center point
      vals[count] = -scalingVars.t0;
      cols[count] = Ncolarray[i][0];
      ++count;

      if (Nrow != -1 && count > 0)
      {
        bs1 = QMatrix.sumIntoRow (Nrow,count,&vals[0],&cols[0]);
        bsuccess = bsuccess && bs1;
      }
      else
      {
        Xyce::dout() << "OOOPS 1!" <<std::endl;
        exit(0);
      }

      // hole continuity row:
      for (j=0;j<numCol;++j) { cols[j] = -1; vals[j] = 0.0; }
      count = 0;

      // center point
      vals[count] = -scalingVars.t0;
      cols[count] = Pcolarray[i][0];
      ++count;

      if (Prow != -1 && count > 0)
      {
        bs1 = QMatrix.sumIntoRow (Prow,count,&vals[0],&cols[0]);
        bsuccess = bsuccess && bs1;
      }
      else
      {
        Xyce::dout() << "OOOPS 2!" <<std::endl;
        exit(0);
      }

    } // mesh loop.
  }
  else // direct matrix access:
  {
    for (i=0;i<numMeshPoints;++i)
    {
      if (boundarySten[i]) continue;

      Nrow = li_Nrowarray[i];
      Prow = li_Prowarray[i];

      std::vector<int> & Noff = li_NoffsetArray[i];
      std::vector<int> & Poff = li_PoffsetArray[i];

      QMatrix[Nrow][Noff[0]] += -scalingVars.t0;
      QMatrix[Prow][Poff[0]] += -scalingVars.t0;

    } // mesh loop.
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdxExtractedConductance
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdxExtractedConductance ()
{
  bool bsuccess = true;
  return bsuccess;
}

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce

