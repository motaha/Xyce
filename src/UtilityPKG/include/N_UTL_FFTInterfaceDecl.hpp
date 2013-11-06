//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_FFTInterfaceDecl.hpp,v $
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
// Creation Date  : 11/11/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2010/11/18 00:08:37 $
//
// Current Owner  : $Author: hkthorn $
//-------------------------------------------------------------------------
#ifndef N_UTL_FFTINTERFACE_DECL_HPP
#define N_UTL_FFTINTERFACE_DECL_HPP


// ---------- Standard Includes ----------

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterfaceDecl
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
class N_UTL_FFTInterfaceDecl
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[FFT/IFT] methods.
    N_UTL_FFTInterfaceDecl( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
      { signalLength_ = length; numberSignals_ = numSignals; stride_ = reqStride; overwrite_ = overwrite; }
    
    // Basic destructor 
    virtual ~N_UTL_FFTInterfaceDecl() {};

    // Register new vectors for the FFT/IFT interface to use.
    virtual void registerVectors( const Teuchos::RCP<const VectorType>& fftInData, const Teuchos::RCP<VectorType>& fftOutData,
                                  const Teuchos::RCP<const VectorType>& iftInData, const Teuchos::RCP<VectorType>& iftOutData ) = 0;
    
    // Calculate FFT with new vectors, not the ones that have been registered or used in the constructor.
    virtual void calculateFFT( const Teuchos::RCP<const VectorType>& inData, const Teuchos::RCP<VectorType>& outData ) = 0;
    // Calculate IFT with new vectors, not the ones that have been registered or used in the constructor.
    virtual void calculateIFT( const Teuchos::RCP<const VectorType>& inData, const Teuchos::RCP<VectorType>& outData ) = 0;

    // Calculate FFT with the vectors that have been registered.
    virtual void calculateFFT() = 0;
    // Calculate IFT with the vectors that have been registered.
    virtual void calculateIFT() = 0;

  protected:
    // this is the length of the real array.  The complex result will be signalLength+2 long
    // if signalLength is even and signalLength+1 if signalLength is odd
    int signalLength_;

    // this is the number of signals on which we will take an fft/ift
    int numberSignals_;

    // If the signals are grouped by blocks at the same time (say x0, x1 ... xn at t0) and
    // then (x0, x1 ... xn at t1). Then stride is the spacing from one x0 at t0 to the next
    // x0 at t1.  This lets one take ffts/ifts of data that is blocked by time
    int stride_;

    // Whether the input and output vectors should be expected to be the same, save space if possible.
    bool overwrite_;

    // FFT/IFT vectors
    Teuchos::RCP<const VectorType> fftInData_, iftInData_;
    Teuchos::RCP<VectorType> fftOutData_, iftOutData_;

};

#endif

