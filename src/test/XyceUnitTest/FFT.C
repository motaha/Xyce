//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: FFT.C,v $
// Purpose       : This file contains unit tests for the FFT interfaces
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/16/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.7 $
// Revision Date  : $Date: 2012/11/02 21:35:07 $
// Current Owner  : $Author: eczeek $
//-----------------------------------------------------------------------------

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <N_UTL_FFTInterface.hpp>
#include <N_LAS_Vector.h>
#include <FFT_Helpers.h>
using Teuchos::RefCountPtr;
using Teuchos::rcp;
using Teuchos::is_null;

TEUCHOS_UNIT_TEST( FFTInterface, create ) {
  int numPts = 11;
  RefCountPtr<N_UTL_FFTInterface<std::vector<double> > > myTransform = rcp(new N_UTL_FFTInterface<std::vector<double> >(numPts) );
  TEST_EQUALITY_CONST( is_null(myTransform), false );
}

TEUCHOS_UNIT_TEST( FFTInterface, std_vector_simple ) {
  StdSignalTransform signal = createStdSimpleSineWave();
  int numPts = signal.numPts;
  int lengthTransformedSignal = signal.numFreqPts;
  N_UTL_FFTInterface<std::vector<double> > myTransform(numPts);
  std::vector<double>& input = *signal.inputSignal;
  std::vector<double> output(lengthTransformedSignal,0.0);
  std::vector<double>& mask = *signal.outputSignalMask;
  std::vector<double>& sol = *signal.outputSignal;
  myTransform.calculateFFT( input, &output );
  double tol=1.0e-10;
  for (int i=0 ; i<lengthTransformedSignal; ++i) {
    TEST_FLOATING_EQUALITY( output[i]+mask[i], sol[i]+mask[i], tol );
  }
}

TEUCHOS_UNIT_TEST( FFTInterface, LAS_Vector_simple ) {
  LAS_SignalTransform signal = createLASSimpleSineWave();
  int numPts = signal.numPts;
  int lengthTransformedSignal = signal.numFreqPts;
  N_UTL_FFTInterface<N_LAS_Vector> myTransform(numPts);
  N_LAS_Vector& input = *signal.inputSignal;
  N_LAS_Vector output(*signal.outputSignalMask);
  output.putScalar(0.0);
  N_LAS_Vector& mask = *signal.outputSignalMask;
  N_LAS_Vector& sol = *signal.outputSignal;
  myTransform.calculateFFT( input, &output );
  double tol=1.0e-10;
  for (int i=0 ; i<lengthTransformedSignal; ++i) {
    TEST_FLOATING_EQUALITY( output[i]+mask[i], sol[i]+mask[i], tol );
  }
}

TEUCHOS_UNIT_TEST( FFTInterface, MultiVec_Vector_simple ) {
  MultiVec_SignalTransform signal = createMultiVecSimpleSineWave();
  int numPts = signal.numPts;
  int lengthTransformedSignal = signal.numFreqPts;
  int numVecs = signal.numVecs;
  N_UTL_FFTInterface<N_LAS_MultiVector> myTransform(numPts);
  N_LAS_MultiVector& input = *signal.inputSignal;
  N_LAS_MultiVector output(*signal.outputSignalMask);
  output.putScalar(0.0);
  N_LAS_MultiVector& mask = *signal.outputSignalMask;
  N_LAS_MultiVector& sol = *signal.outputSignal;
  myTransform.calculateFFT( input, &output );
  double tol=1.0e-10;
  for (int j=0 ; j<numVecs ; ++j )
  {
    for (int i=0 ; i<lengthTransformedSignal; ++i) {
      TEST_FLOATING_EQUALITY( output[j][i]+mask[j][i], sol[j][i]+mask[j][i], tol );
    }
  }
}

TEUCHOS_UNIT_TEST( FFTInterface, LAS_Vector_nonSimple ) {
  LAS_SignalTransform signal = createLASModifiedSineWave();
  int numPts = signal.numPts;
  int lengthTransformedSignal = signal.numFreqPts;
  N_UTL_FFTInterface<N_LAS_Vector> myTransform(numPts);
  N_LAS_Vector& input = *signal.inputSignal;
  N_LAS_Vector output(*signal.outputSignalMask);
  output.putScalar(0.0);
  N_LAS_Vector& mask = *signal.outputSignalMask;
  N_LAS_Vector& sol = *signal.outputSignal;
  myTransform.calculateFFT( input, &output );
  double tol=1.0e-10;
  for (int i=0 ; i<lengthTransformedSignal; ++i) {
    TEST_FLOATING_EQUALITY( output[i]+mask[i], sol[i]+mask[i], tol );
  }
}

TEUCHOS_UNIT_TEST( FFTInterface, std_vector_nonSimple ) {
  StdSignalTransform signal = createStdModifiedSineWave();
  int numPts = signal.numPts;
  int lengthTransformedSignal = signal.numFreqPts;
  N_UTL_FFTInterface<std::vector<double> > myTransform(numPts);
  std::vector<double>& input = *signal.inputSignal;
  std::vector<double> output(lengthTransformedSignal);
  std::vector<double>& mask = *signal.outputSignalMask;
  std::vector<double>& sol = *signal.outputSignal;
  myTransform.calculateFFT( input, &output );
  double tol=1.0e-10;
  for (int i=0 ; i<lengthTransformedSignal; ++i) {
    TEST_FLOATING_EQUALITY( output[i]+mask[i], sol[i]+mask[i], tol );
  }
}

