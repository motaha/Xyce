//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: HB_Loader_Helpers.h,v $
// Purpose       : This file contains some helper functions for create N_HB_Loaders.
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/10/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.2 $
// Revision Date  : $Date: 2008/09/18 17:04:27 $
// Current Owner  : $Author: tscoffe $
//-----------------------------------------------------------------------------

#ifndef HB_LOADER_HELPERS_H
#define HB_LOADER_HELPERS_H

#include <Teuchos_RefCountPtr.hpp>
#include <N_LAS_Vector.h>
#include <N_LOA_Loader.h>
#include <N_LAS_Matrix.h>

using Teuchos::RefCountPtr;

class N_HB_Loader;
class N_HB_Builder;
class N_LOA_Loader;
class N_LAS_Builder;

class SinCosLoader : public virtual N_LOA_Loader
{
  public:
    SinCosLoader() { period_ = 1.0; pi_ = 4.0*atan(1.0); }
    virtual ~SinCosLoader() {}
    bool loadJacobian () { return false; }
    bool loadRHS () { return false; }
    bool loadDAEMatrices(
        N_LAS_Vector * tmpSolVectorPtr,
        N_LAS_Vector * tmpStaVectorPtr,
        N_LAS_Vector * tmpStaDerivVectorPtr,
        N_LAS_Matrix * tmpdQdxMatrixPtr,
        N_LAS_Matrix * tmpdFdxMatrixPtr
        )
    {
      N_LAS_Vector& xVector = *tmpSolVectorPtr;
      N_LAS_Matrix& qMatrix = *tmpdQdxMatrixPtr;
      N_LAS_Matrix& fMatrix = *tmpdFdxMatrixPtr;
      qMatrix.put(0.0);
      fMatrix.put(0.0);

      double coeffs[2];
      int colIndices[2];
      colIndices[0] = 0;
      colIndices[1] = 1;
      // Q matrix
      coeffs[0] = 1.0;
      coeffs[1] = 0.0;
      qMatrix.insertRow(0,2,coeffs,colIndices); 
      coeffs[0] = 0.0;
      coeffs[1] = 1.0;
      qMatrix.insertRow(1,2,coeffs,colIndices);
      // F matrix
      coeffs[0] = 0.0;
      coeffs[1] = 2.0*pi_/period_;
      fMatrix.insertRow(0,2,coeffs,colIndices);
      coeffs[0] = -2.0*pi_/period_;
      coeffs[1] = 0.0;
      fMatrix.insertRow(1,2,coeffs,colIndices);

      qMatrix.fillComplete();
      fMatrix.fillComplete();
      return true;
    }
    bool loadDAEVectors(
        N_LAS_Vector * tmpSolVectorPtr,
        N_LAS_Vector * tmpStaVectorPtr,
        N_LAS_Vector * tmpcurrStaVectorPtr,
        N_LAS_Vector * tmplastStaVectorPtr,
        N_LAS_Vector * tmpStaDerivVectorPtr,
        N_LAS_Vector * tmpQVectorPtr,
        N_LAS_Vector * tmpFVectorPtr,
        N_LAS_Vector * tmpdFdxdVpVectorPtr,
        N_LAS_Vector * tmpdQdxdVpVectorPtr
        )
    { 
      N_LAS_Vector& fVector = *tmpFVectorPtr;
      N_LAS_Vector& qVector = *tmpQVectorPtr;
      N_LAS_Vector& xVector = *tmpSolVectorPtr;
      fVector[0] = (2.0*pi_/period_)*xVector[1];
      fVector[1] = -(2.0*pi_/period_)*xVector[0];
      qVector[0] = xVector[0];
      qVector[1] = xVector[1];
      return true;
    }
    bool applyDAEMatrices(
        N_LAS_Vector * tmpSolVectorPtr,
        N_LAS_Vector * tmpStaVectorPtr,
        N_LAS_Vector * tmpStaDerivVectorPtr,
        const N_LAS_Vector & tmpVecVectorPtr,
        N_LAS_Vector * tmpdQdxVecVectorPtr,
        N_LAS_Vector * tmpdFdxVecVectorPtr
        )
    {
      const N_LAS_Vector& vVector = tmpVecVectorPtr;
      N_LAS_Vector& dQdxV = *tmpdQdxVecVectorPtr;
      N_LAS_Vector& dFdxV = *tmpdFdxVecVectorPtr;
      dQdxV[0] = vVector[0];
      dQdxV[1] = vVector[1];
      dFdxV[0] = (2.0*pi_/period_)*vVector[1];
      dFdxV[1] = -(2.0*pi_/period_)*vVector[0];
      return true;
    }
    void setPeriod(double period) { period_ = period; }
    double getPeriod() { return period_; }
  private:
    double period_;
    double pi_;
};

RefCountPtr<N_HB_Loader> createHBLoader(RefCountPtr<N_LOA_Loader> appLoader, bool matrixFreeFlag);
RefCountPtr<N_LOA_Loader> createAppLoader();
RefCountPtr<N_LAS_Builder> createAppBuilder(int numSolutionVars, int numStateVars);
void registerVectorsOnHBLoader(N_HB_Loader* hbLoader, const N_HB_Builder& hbBuilder, const N_LAS_Builder& appBuilder);

RefCountPtr<N_LAS_MultiVector> loadJacobian(N_LOA_Loader& loader, const N_LAS_Builder& builder);
RefCountPtr<N_LAS_MultiVector> diffJacobian(N_LOA_Loader& loader, const N_LAS_Builder& builder);
RefCountPtr<N_LAS_MultiVector> applyJacobian(N_LOA_Loader& loader, const N_LAS_Builder& builder);

#endif // HB_LOADER_HELPERS_H

