//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: FFT_Helpers.h,v $
// Purpose       : This file contains some helper functions for the FFT Interface
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/16/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.3 $
// Revision Date  : $Date: 2008/09/25 21:32:34 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#ifndef FFT_HELPERS_H
#define FFT_HELPERS_H

#include <Teuchos_RefCountPtr.hpp>
#include <vector>

using Teuchos::RefCountPtr;
using Teuchos::rcp;

class N_LAS_Vector;
class N_LAS_MultiVector;
class Epetra_Map;

struct StdSignalTransform
{
  int numPts;
  int numFreqPts;
  RefCountPtr<std::vector<double> > inputSignal;
  RefCountPtr<std::vector<double> > outputSignal;
  // the Mask is necessary for getting around the relative error test in TEST_FLOATING_EQUALITY
  RefCountPtr<std::vector<double> > outputSignalMask;
};

StdSignalTransform createStdSimpleSineWave();
StdSignalTransform createStdModifiedSineWave();

struct LAS_SignalTransform
{
  int numPts;
  int numFreqPts;
  RefCountPtr<N_LAS_Vector> inputSignal;
  RefCountPtr<N_LAS_Vector> outputSignal;
  // the Mask is necessary for getting around the relative error test in TEST_FLOATING_EQUALITY
  RefCountPtr<N_LAS_Vector> outputSignalMask;
  RefCountPtr<const Epetra_Map> emap_input;
  RefCountPtr<const Epetra_Map> emap_output;
};

LAS_SignalTransform createLASSimpleSineWave();
LAS_SignalTransform createLASModifiedSineWave();

struct MultiVec_SignalTransform
{
  int numPts;
  int numFreqPts;
  int numVecs;
  RefCountPtr<N_LAS_MultiVector> inputSignal;
  RefCountPtr<N_LAS_MultiVector> outputSignal;
  // the Mask is necessary for getting around the relative error test in TEST_FLOATING_EQUALITY
  RefCountPtr<N_LAS_MultiVector> outputSignalMask;
  RefCountPtr<const Epetra_Map> emap_input;
  RefCountPtr<const Epetra_Map> emap_output;
};

MultiVec_SignalTransform createMultiVecSimpleSineWave();
MultiVec_SignalTransform createMultiVecModifiedSineWave();

#endif // FFT_HELPERS_H

