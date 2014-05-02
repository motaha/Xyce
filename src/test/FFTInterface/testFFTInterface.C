//
// test the N_UTL_FFTInterface class
//

#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_CommFactory.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Vector.h>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>

#include <iostream>
#include <vector>
#include <cmath>

int main(int argc, char* argv[])
{
  //
  // This first part of the code tests out making a signal x(t) and then
  // taking its fourier transform fft( x(t) ) and then its back 
  // transform ift( fft( x(t) ) ).  This involves running the FFT interface
  // with a single input vector (either a vector<double> or an N_LAS_Vector
  // object.  In the second half of this code we try doing multiple FFTs at
  // the same time using an N_LAS_BlockVector object
  //
  
  // create arrays to hold some data for testing
  //int numPts = 128;
  int numPts = 11;
  int lengthTransformedSignal = numPts + 2;
  if( numPts % 2 )
  {
    // numPts was odd so the actual number of elements in the transformed signal is numPts+1
    lengthTransformedSignal = numPts + 1;
  }
  std::cout << "numPts = " << numPts << std::endl;
  std::cout << "lengthTransformedSignal = " << lengthTransformedSignal << std::endl;
  
  double timeStart = 0.0;
  double timeStop = 1;
  double deltaTime = (timeStop - timeStart)/numPts;
  double freqDelta = 1.0 / (timeStop - timeStart);
  
  // these are to test the vector of double accecess to FFT interface 
  std::vector<double> time(numPts, 0.0);
  std::vector<double> inputSignal(numPts, 0.0);
  std::vector<double> outputSignal(lengthTransformedSignal, 0.0);
  std::vector<double> backSignal(numPts, 0.0);
  
  // these is to make N_LAS_MultiVector objects to test that part of
  // the interface since N_LAS_Vector derives from MultiVector, one
  // can use this same interface for N_LAS_Vector objects
  Teuchos::RCP<N_PDS_Comm> pdsComm = Teuchos::rcp( N_PDS_CommFactory::create( argc, argv ) );
  std::vector<int> realLbMap(numPts, 0.0);
  std::vector<int> complexLbMap(lengthTransformedSignal, 0.0);
  N_PDS_ParMap parMapForReal( numPts, numPts, realLbMap, 0, pdsComm.get() );
  N_PDS_ParMap parMapForComplex( lengthTransformedSignal, lengthTransformedSignal, complexLbMap, 0, pdsComm.get() );
  
  N_LAS_Vector timeMV( parMapForReal );
  N_LAS_Vector inputSignalMV( parMapForReal );
  N_LAS_Vector outputSignalMV( parMapForComplex );
  N_LAS_Vector backSignalMV( parMapForReal );
  
  // Fill in the vectors
  for(int i=0; i<numPts; i++)
  {
    time[i] = timeStart + i * deltaTime;
    
    inputSignal[i] = std::sin( 2.0 * M_PI * 1 * time[i] ) + time[i]*time[i];      }

  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > myTransform( numPts );
  
  // try taking an FFT
  myTransform.calculateFFT( inputSignal, &outputSignal );
  myTransform.calculateIFT( outputSignal, &backSignal );
  
  std::cout.precision(6);
  Xyce::dout() << "time\tInputSignal\tBackSignal" << std::endl;
  for(int i=0; i<numPts; i++)
  {
    Xyce::dout() << time[i] << "\t" << inputSignal[i] << "\t" << backSignal[i] << std::endl;
     
  }
  Xyce::dout() << "Frequency\tReal + Imaginary" << std::endl;
  for (int i=0 ; i<lengthTransformedSignal/2  ; ++i) {
      Xyce::dout() << i*freqDelta << "\t" << outputSignal[2*i] << " + "
        << outputSignal[2*i+1] << "i" << std::endl;
  }
  
  return 0;
}
