//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_ExternData.h,v $
//
// Purpose        : This class is a container class for holding pointers,
//                  references, etc., to external linear solver data
//                  structures.  Occasionally, it will hold pointers
//                  to other things, but that is not the primary intention.
//
//                  In general, stuff that goes into this class should
//                  be stuff needed by more than one device instance type.
//
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.42.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_ExternData_h
#define Xyce_N_DEV_ExternData_h

#include <map>
#include <vector>
#include <N_UTL_Xyce.h>
#include <N_DEV_fwd.h>

// ---------- Forward Declarations ----------
class N_LAS_Matrix;
class N_LAS_MultiVector;
class N_LAS_Vector;
class N_LAS_System;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : ExternData
// Purpose       : Container for references, pointers, etc. to data structures
//                 outside the device package.  Mostly these are linear
//                 algebra objects like the Jacobian.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/25/02
//-----------------------------------------------------------------------------
class ExternData
{

  public:
    ExternData () :
      lasSysPtr(0),
#ifdef Xyce_DEBUG_DEVICE
      fVectorPtr(0),
      JdxpVectorPtr(0),
#endif
#ifdef Xyce_DEBUG_VOLTLIM
      dxVoltlimVectorPtr(0),

      JTestMatrixPtr(0),
      Jdx2VectorPtr(0),

      FTestMatrixPtr(0),
      Fdx2VectorPtr(0),
      QTestMatrixPtr(0),
      Qdx2VectorPtr(0),
#endif
      tmpdIdXPtr(0),
      tmpdQdXPtr(0),
      daeQVectorPtr(0),
      daeFVectorPtr(0),

      bVecRealPtr(0),
      bVecImagPtr(0),

      dFdxdVpVectorPtr(0),
      dQdxdVpVectorPtr(0),
      dQdxMatrixPtr(0),
      dFdxMatrixPtr(0),

      currSolVectorPtr(0),
      nextSolVectorPtr(0),
      lastSolVectorPtr(0),

      currStaVectorPtr(0),
      nextStaVectorPtr(0),
      lastStaVectorPtr(0),

      currStoVectorPtr(0),
      nextStoVectorPtr(0),
      lastStoVectorPtr(0),
      storeLeadCurrQCompPtr(0),
      flagSolVectorPtr(0),

      perturbVectorPtr(0),
      numJacRHSVectorPtr(0),
      numJacFVectorPtr(0),
      numJacQVectorPtr(0),
      numJacLoadFlagPtr(0),

#if 0
      nextSolDerivVectorPtr(0),
#endif
      nextStaDerivVectorPtr(0),

      deviceMaskVectorPtr(0),

      daeQVectorRawPtr(0),
      daeFVectorRawPtr(0),
      dFdxdVpVectorRawPtr(0),
      dQdxdVpVectorRawPtr(0),
      nextSolVectorRawPtr(0),
      currSolVectorRawPtr(0),
      lastSolVectorRawPtr(0),
      nextStaVectorRawPtr(0),
      currStaVectorRawPtr(0),
      lastStaVectorRawPtr(0),
      nextStaDerivVectorRawPtr(0),
      nextStoVectorRawPtr(0),
      currStoVectorRawPtr(0),
      lastStoVectorRawPtr(0),
      storeLeadCurrQCompRawPtr(0),
      bVecRealRawPtr(0),
      bVecImagRawPtr(0),

      devMgrPtr(0),
      solDevInstMap(0),
      initializeAllFlag(false)
    { };

  protected:

  private:

  public:
    N_LAS_System * lasSysPtr;

#ifdef Xyce_DEBUG_DEVICE
    N_LAS_Vector * fVectorPtr;
    N_LAS_Vector * JdxpVectorPtr;
#endif
#ifdef Xyce_DEBUG_VOLTLIM
    // voltlim DX vector:
    N_LAS_Vector * dxVoltlimVectorPtr;

    // old DAE matrix
    N_LAS_Vector * Jdx2VectorPtr;
    N_LAS_Matrix * JTestMatrixPtr;

    // new DAE matrices
    N_LAS_Vector * Fdx2VectorPtr;
    N_LAS_Matrix * FTestMatrixPtr;
    N_LAS_Vector * Qdx2VectorPtr;
    N_LAS_Matrix * QTestMatrixPtr;
#endif
    N_LAS_Vector * tmpdIdXPtr;
    N_LAS_Vector * tmpdQdXPtr;

    // DAE formulation vectors


    N_LAS_Vector * bVecRealPtr;
    N_LAS_Vector * bVecImagPtr;

    N_LAS_Vector * daeQVectorPtr;
    N_LAS_Vector * daeFVectorPtr;

    N_LAS_Vector *  dFdxdVpVectorPtr;
    N_LAS_Vector *  dQdxdVpVectorPtr;

    // DAE formulation matrices
    N_LAS_Matrix * dQdxMatrixPtr;
    N_LAS_Matrix * dFdxMatrixPtr;

    N_LAS_Vector * currSolVectorPtr;
    N_LAS_Vector * nextSolVectorPtr;
    N_LAS_Vector * lastSolVectorPtr;

    N_LAS_Vector * currStaVectorPtr;
    N_LAS_Vector * nextStaVectorPtr;
    N_LAS_Vector * lastStaVectorPtr;

    N_LAS_Vector * currStoVectorPtr;
    N_LAS_Vector * nextStoVectorPtr;
    N_LAS_Vector * lastStoVectorPtr;
    N_LAS_Vector * storeLeadCurrQCompPtr;

    N_LAS_Vector * flagSolVectorPtr;

    N_LAS_Vector *  perturbVectorPtr;
    N_LAS_Vector *  numJacRHSVectorPtr;
    N_LAS_Vector *  numJacFVectorPtr;
    N_LAS_Vector *  numJacQVectorPtr;
    N_LAS_Vector *  numJacLoadFlagPtr;

#if 0
    N_LAS_Vector  * nextSolDerivVectorPtr;
#endif
    N_LAS_Vector  * nextStaDerivVectorPtr;

    N_LAS_Vector  * deviceMaskVectorPtr;

    // raw pointers (to internal vector data):
    double * daeQVectorRawPtr;
    double * daeFVectorRawPtr;
    double * dFdxdVpVectorRawPtr;
    double * dQdxdVpVectorRawPtr;

    double * nextSolVectorRawPtr;
    double * currSolVectorRawPtr;
    double * lastSolVectorRawPtr;

    double * nextStaVectorRawPtr;
    double * currStaVectorRawPtr;
    double * lastStaVectorRawPtr;

    double * nextStoVectorRawPtr;
    double * currStoVectorRawPtr;
    double * lastStoVectorRawPtr;
    double * storeLeadCurrQCompRawPtr;

    double * nextStaDerivVectorRawPtr;

    double * bVecRealRawPtr;
    double * bVecImagRawPtr;

    // This device class manager pointer should go somewhere
    // else eventually.
    N_DEV_DeviceMgr * devMgrPtr;

    // This map should go elsewhere eventually also:
    multimap <int, DeviceInstance*> * solDevInstMap;

    bool initializeAllFlag;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ExternData N_DEV_ExternData;

#endif

