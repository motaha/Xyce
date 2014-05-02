//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_FFTInterface.hpp,v $
//
// Purpose        : This class acts as an interface to an FFT library
//                  for FFT and IFT calculations.  This class should isolate
//                  Xyce from the specifics of a given FFT library so 
//                  that multiple libraries can be used.  It is originally
//                  implemented for Intel's Math Library but may be extended
//                  to FFTW at some time in the future.
//
// Special Notes  : 
//
// Creator        : Richard Schiek 
//
// Creation Date  : 5/27/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6 $
//
// Revision Date  : $Date: 2014/01/28 19:03:42 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#ifndef N_UTL_FFTINTERFACE_HPP
#define N_UTL_FFTINTERFACE_HPP 1


// ---------- Standard Includes ----------
#include <Xyce_config.h>

#ifdef Xyce_USE_INTEL_FFT
#include <N_UTL_IntelFFT_Interface.hpp>
#endif

#ifdef Xyce_USE_FFTW
#include <N_UTL_FFTW_Interface.hpp>
#endif

#include <N_ERH_ErrorMgr.h>

// ----------   Other Includes   ----------

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterface
// Purpose       : This class acts as a templated interface to any FFT library
//                 for FFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given FFT library so 
//                 that multiple libraries can be used.  It is originally
//                 implemented for Intel's Math Library but may be extended
//                 to FFTW at some time in the future.
// Special Notes :
// Creator       : Richard Schiek (templated by Heidi Thornquist)
// Creation Date : 5/27/08
// Last Modified : 11/17/10
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_FFTInterface
{
  public:
    N_UTL_FFTInterface( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
#ifdef Xyce_USE_INTEL_FFT
      : intelfftInterface( length, numSignals, reqStride, overwrite )
#elif defined(Xyce_USE_FFTW)
      : fftwInterface( length, numSignals, reqStride, overwrite )
#endif
    {
#ifndef Xyce_USE_FFT
      std::string msg = "Xyce has not been configured with an FFT library! Please reconfigure with FFT enabled to perform any frequency analysis!\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
#endif
    }
    
    virtual ~N_UTL_FFTInterface() {}
   
    // Register new vectors for the FFT/IFT interface to use.
    void registerVectors( const VectorType& fftInData, VectorType* fftOutData,
                          const VectorType& iftInData, VectorType* iftOutData )
    {
#ifdef Xyce_USE_INTEL_FFT
      intelfftInterface.registerVectors( Teuchos::rcp( &fftInData, false ), Teuchos::rcp( fftOutData, false ),
                                         Teuchos::rcp( &iftInData, false ), Teuchos::rcp( iftOutData, false ) );
#elif defined(Xyce_USE_FFTW)
      fftwInterface.registerVectors( Teuchos::rcp( &fftInData, false ), Teuchos::rcp( fftOutData, false ),
                                     Teuchos::rcp( &iftInData, false ), Teuchos::rcp( iftOutData, false ) );
#endif
    }
 
    void calculateFFT( const VectorType& inData, VectorType* outResult)
    {
#ifdef Xyce_USE_INTEL_FFT
      intelfftInterface.calculateFFT( Teuchos::rcp( &inData, false ), Teuchos::rcp( outResult, false ) );
#elif defined(Xyce_USE_FFTW)
      fftwInterface.calculateFFT( Teuchos::rcp( &inData, false ), Teuchos::rcp( outResult, false ) );
#endif
    }
    
    void calculateIFT( const VectorType& inData, VectorType* outResult)
    {
#ifdef Xyce_USE_INTEL_FFT
      intelfftInterface.calculateIFT( Teuchos::rcp( &inData, false ), Teuchos::rcp( outResult, false ) );
#elif defined(Xyce_USE_FFTW)
      fftwInterface.calculateIFT( Teuchos::rcp( &inData, false ), Teuchos::rcp( outResult, false ) );
#endif
    }
   
    void calculateFFT()
    {
#ifdef Xyce_USE_INTEL_FFT
      intelfftInterface.calculateFFT();
#elif defined(Xyce_USE_FFTW)
      fftwInterface.calculateFFT();
#endif
    }
    
    void calculateIFT()
    {
#ifdef Xyce_USE_INTEL_FFT
      intelfftInterface.calculateIFT();
#elif defined(Xyce_USE_FFTW)
      fftwInterface.calculateIFT();
#endif
    }

  private:

#ifdef Xyce_USE_INTEL_FFT
    N_UTL_IntelFFT_Interface<VectorType> intelfftInterface;
#elif defined(Xyce_USE_FFTW)
    N_UTL_FFTW_Interface<VectorType> fftwInterface;
#endif    
};

#endif

