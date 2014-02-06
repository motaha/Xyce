//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: FFT_Helpers.C,v $
// Purpose       : This file contains some helper functions for the FFT interface
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/16/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.4 $
// Revision Date  : $Date: 2012/11/02 21:35:07 $
// Current Owner  : $Author: eczeek $
//-----------------------------------------------------------------------------

#include <FFT_Helpers.h>
#include <N_LAS_Vector.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <N_PDS_SerialComm.h>

StdSignalTransform createStdSimpleSineWave()
{
  StdSignalTransform signal;
  signal.numPts = 11;
  signal.numFreqPts = 12;
  signal.inputSignal = rcp( new std::vector<double>(signal.numPts, 0.0) );
  signal.outputSignal = rcp( new std::vector<double>(signal.numFreqPts, 0.0) );
  signal.outputSignalMask = rcp( new std::vector<double>(signal.numFreqPts, 1.0) );
  double pi = 4.0*atan(1.0);
  for (int i=0 ; i<signal.numPts ; ++i) {
    (*signal.inputSignal)[i] = std::sin(2*pi*i/signal.numPts); 
  }
  (*signal.outputSignal)[3] = -5.5;
  (*signal.outputSignalMask)[3] = 0.0;
  return signal;
}

LAS_SignalTransform createLASSimpleSineWave()
{
  StdSignalTransform stdSignal = createStdSimpleSineWave();
  LAS_SignalTransform lasSignal;
  lasSignal.numPts = stdSignal.numPts;
  lasSignal.numFreqPts = stdSignal.numFreqPts;
  RefCountPtr<N_PDS_Comm> pdsComm = rcp(new N_PDS_SerialComm);
  RefCountPtr<Epetra_Comm> comm = rcp(pdsComm->petraComm(),false);
  lasSignal.emap_input = rcp(new Epetra_Map(lasSignal.numPts,0,*comm));
  lasSignal.emap_output = rcp(new Epetra_Map(lasSignal.numFreqPts,0,*comm));
  RefCountPtr<Epetra_Vector> evec_input = rcp(new Epetra_Vector(*lasSignal.emap_input),false);
  RefCountPtr<Epetra_Vector> evec_output = rcp(new Epetra_Vector(*lasSignal.emap_output),false);
  RefCountPtr<Epetra_Vector> evec_outputMask = rcp(new Epetra_Vector(*lasSignal.emap_output),false);
  lasSignal.inputSignal = rcp(new N_LAS_Vector(&*evec_input,true));
  lasSignal.outputSignal = rcp(new N_LAS_Vector(&*evec_output,true));
  lasSignal.outputSignalMask = rcp(new N_LAS_Vector(&*evec_outputMask,true));
  for (int i=0 ; i<lasSignal.numPts ; ++i) {
    (*lasSignal.inputSignal)[i] = (*stdSignal.inputSignal)[i];
  }
  for (int i=0 ; i<lasSignal.numFreqPts ; ++i) {
    (*lasSignal.outputSignal)[i] = (*stdSignal.outputSignal)[i];
    (*lasSignal.outputSignalMask)[i] = (*stdSignal.outputSignalMask)[i];
  }
  return lasSignal;
}

MultiVec_SignalTransform createMultiVecSimpleSineWave()
{
  StdSignalTransform stdSignal = createStdSimpleSineWave();
  MultiVec_SignalTransform multiVecSignal;
  multiVecSignal.numPts = stdSignal.numPts;
  multiVecSignal.numFreqPts = stdSignal.numFreqPts;
  multiVecSignal.numVecs = 10;   //the number of vecs in this multivector
  
  // make the epetra objects that can be used to make the multivectors
  RefCountPtr<N_PDS_Comm> pdsComm = rcp(new N_PDS_SerialComm);
  RefCountPtr<Epetra_Comm> comm = rcp(pdsComm->petraComm(),false);
  multiVecSignal.emap_input = rcp(new Epetra_Map(multiVecSignal.numPts,0,*comm));
  multiVecSignal.emap_output = rcp(new Epetra_Map(multiVecSignal.numFreqPts,0,*comm));
  RefCountPtr<Epetra_MultiVector> evec_input = rcp(new Epetra_MultiVector(*multiVecSignal.emap_input, multiVecSignal.numVecs),false);
  RefCountPtr<Epetra_MultiVector> evec_output = rcp(new Epetra_MultiVector(*multiVecSignal.emap_output, multiVecSignal.numVecs),false);
  RefCountPtr<Epetra_MultiVector> evec_outputMask = rcp(new Epetra_MultiVector(*multiVecSignal.emap_output, multiVecSignal.numVecs),false);
  
  // now make the multivectors
  multiVecSignal.inputSignal = rcp(new N_LAS_MultiVector(&*evec_input,true));
  multiVecSignal.outputSignal = rcp(new N_LAS_MultiVector(&*evec_output,true));
  multiVecSignal.outputSignalMask = rcp(new N_LAS_MultiVector(&*evec_outputMask,true));
  for (int j=0; j<multiVecSignal.numVecs; ++j)
  {
    for (int i=0 ; i<multiVecSignal.numPts ; ++i) {
      (*multiVecSignal.inputSignal)[j][i] = (*stdSignal.inputSignal)[i];
    }
    for (int i=0 ; i<multiVecSignal.numFreqPts ; ++i) {
      (*multiVecSignal.outputSignal)[j][i] = (*stdSignal.outputSignal)[i];
      (*multiVecSignal.outputSignalMask)[j][i] = (*stdSignal.outputSignalMask)[i];
    }
  }
  return multiVecSignal;
}


StdSignalTransform createStdModifiedSineWave()
{
  StdSignalTransform signal;
  signal.numPts = 11;
  signal.numFreqPts = 12;
  signal.inputSignal = rcp( new std::vector<double>(signal.numPts, 0.0) );
  signal.outputSignal = rcp( new std::vector<double>(signal.numFreqPts, 0.0) );
  signal.outputSignalMask = rcp( new std::vector<double>(signal.numFreqPts, 0.0) );
  double pi = 4.0*atan(1.0);
  for (int i=0 ; i<signal.numPts ; ++i) {
    (*signal.inputSignal)[i] = std::sin(2*pi*i/signal.numPts)+pow(1.0*i/signal.numPts,2.0);
  }
  (*signal.outputSignal)[0] = 3.18181818181818;
  (*signal.outputSignal)[1] = 0.0;
  (*signal.outputSignal)[2] = 0.07266843496059;
  (*signal.outputSignal)[3] = -3.79715638055538;
  (*signal.outputSignal)[4] = -0.34448952174445;
  (*signal.outputSignal)[5] = 0.77801518648772;
  (*signal.outputSignal)[6] = -0.42041678190357;
  (*signal.outputSignal)[7] = 0.43325246627153;
  (*signal.outputSignal)[8] = -0.44506541303185;
  (*signal.outputSignal)[9] = 0.22834234895178;
  (*signal.outputSignal)[10] = -0.45360580918981;
  (*signal.outputSignal)[11] = 0.07188914699749;
  (*signal.outputSignalMask)[1] = 1.0;
  return signal;
}

LAS_SignalTransform createLASModifiedSineWave()
{
  StdSignalTransform stdSignal = createStdModifiedSineWave();
  LAS_SignalTransform lasSignal = createLASSimpleSineWave();
  int numPts = lasSignal.numPts;
  int numFreqPts = lasSignal.numFreqPts;
  for (int i=0 ; i<numPts ; ++i) {
    (*lasSignal.inputSignal)[i] = (*stdSignal.inputSignal)[i];
  }
  for (int i=0 ; i<numFreqPts ; ++i) {
    (*lasSignal.outputSignal)[i] = (*stdSignal.outputSignal)[i];
    (*lasSignal.outputSignalMask)[i] = (*stdSignal.outputSignalMask)[i];
  }
  return lasSignal;
}

MultiVec_SignalTransform createMultiVecModifiedSineWave()
{
  StdSignalTransform stdSignal = createStdModifiedSineWave();
  MultiVec_SignalTransform multiVecSignal = createMultiVecSimpleSineWave();
  int numPts = multiVecSignal.numPts;
  int numFreqPts = multiVecSignal.numFreqPts;
  int numVecs = multiVecSignal.numVecs;
  for (int j=0 ; j<numVecs; ++j) {
    for (int i=0 ; i<numPts ; ++i) {
      (*multiVecSignal.inputSignal)[j][i] = (*stdSignal.inputSignal)[i];
    }
    for (int i=0 ; i<numFreqPts ; ++i) {
      (*multiVecSignal.outputSignal)[j][i] = (*stdSignal.outputSignal)[i];
      (*multiVecSignal.outputSignalMask)[j][i] = (*stdSignal.outputSignalMask)[i];
    }
  }
  return multiVecSignal;
}
