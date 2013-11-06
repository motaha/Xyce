//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_FFTW_Interface.hpp,v $
//
// Purpose        : This class acts as an abstract interface to an FFT library
//                  for FFT and IFT calculations.  This class should isolate
//                  Xyce from the specifics of a given FFT library so 
//                  that multiple libraries can be used.  
//
// Special Notes  : 
//
// Creator        : Heidi Thornquist
//
// Creation Date  : 5/27/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2013/06/06 22:42:52 $
//
// Current Owner  : $Author: hkthorn $
//-------------------------------------------------------------------------
#ifndef N_UTL_FFTW_INTERFACE_HPP
#define N_UTL_FFTW_INTERFACE_HPP


// ---------- Standard Includes ----------

#include <N_UTL_FFTInterfaceDecl.hpp>

// ----------   Other Includes   ----------

#include <fftw3.h>
#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterface for FFTW
// Purpose       : This class acts as an abstract interface to an FFT library
//                 for FFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given FFT library so 
//                 that multiple libraries can be used.  
// Special Notes :
// Creator       : Heidi Thornquist
// Creation Date : 11/11/10
// Last Modified : 11/11/10
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_FFTW_Interface: public N_UTL_FFTInterfaceDecl<VectorType>
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[FFT/IFT] methods.
    N_UTL_FFTW_Interface( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
      : N_UTL_FFTInterfaceDecl<VectorType>(length, numSignals, reqStride, overwrite),
        firstForwardFFT_(true),
        firstInverseFFT_(true)
    {}

    // Basic destructor 
    virtual ~N_UTL_FFTW_Interface() 
    {
      if (!firstForwardFFT_)
        fftw_destroy_plan(forwardPlan_);
      if (!firstInverseFFT_)
        fftw_destroy_plan(inversePlan_);
      fftw_cleanup();
    }

    // Register new vectors for the FFT/IFT interface to use.
    void registerVectors( const Teuchos::RCP<const VectorType>& fftInData, const Teuchos::RCP<VectorType>& fftOutData,
                          const Teuchos::RCP<const VectorType>& iftInData, const Teuchos::RCP<VectorType>& iftOutData )
    { 
      this->fftInData_ = fftInData; 
      this->fftOutData_ = fftOutData; 
      this->iftInData_ = iftInData; 
      this->iftOutData_ = iftOutData; 
      bool cleanup = false;
      if (!firstForwardFFT_)
      {
        fftw_destroy_plan(forwardPlan_);
        firstForwardFFT_ = true;
        cleanup = true;
      }
      if (!firstInverseFFT_)
      {
        fftw_destroy_plan(inversePlan_);
        firstInverseFFT_ = true;
        cleanup = true;
      }
      if (cleanup)
        fftw_cleanup();
    }

    // Calculate FFT with new vectors, not the ones that have been registered or used in the constructor.
    void calculateFFT( const Teuchos::RCP<const VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      // Check if the vectors are new, if so destroy the current plan an rebuild it in the next call to calculateFFT()
      if ((inData != this->fftInData_) || (outData != this->fftOutData_))
      {
        this->fftInData_ = inData;
        this->fftOutData_ = outData;
        if (!firstForwardFFT_)
        {
          fftw_destroy_plan(forwardPlan_);
          firstForwardFFT_ = true;
        }
      }
      calculateFFT();
    }

    // Calculate IFT with new vectors, not the ones that have been registered or used in the constructor.
    void calculateIFT( const Teuchos::RCP<const VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      // Check if the vectors are new, if so destroy the current plan an rebuild it in the next call to calculateIFT()
      if ((inData != this->iftInData_) || (outData != this->iftOutData_))
      {
        this->iftInData_ = inData;
        this->iftOutData_ = outData;
        if (!firstInverseFFT_)
        {
          fftw_destroy_plan(inversePlan_);
          firstInverseFFT_ = true;
        }
      }
      calculateIFT();
    }

    // Calculate FFT with the vectors that have been registered.
    // NOTE:  This method must be specialized for each type of vector used by this class,
    //        or the lack of method definition will result in a build failure.
    void calculateFFT();
    // Calculate IFT with the vectors that have been registered.
    // NOTE:  This method must be specialized for each type of vector used by this class,
    //        or the lack of method definition will result in a build failure.
    void calculateIFT();

  private:
    bool firstForwardFFT_, firstInverseFFT_;
    fftw_plan forwardPlan_;
    fftw_plan inversePlan_;
};

#endif
