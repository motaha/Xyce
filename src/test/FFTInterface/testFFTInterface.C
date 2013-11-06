//
// test the N_UTL_FFTInterface class
//

#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Misc.h>
#include <N_PDS_ParMap.h>
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
  cout << "numPts = " << numPts << std::endl;
  cout << "lengthTransformedSignal = " << lengthTransformedSignal << std::endl;
  
  double timeStart = 0.0;
  double timeStop = 1;
  double deltaTime = (timeStop - timeStart)/numPts;
  double freqDelta = 1.0 / (timeStop - timeStart);
  
  // these are to test the vector of double accecess to FFT interface 
  std::vector<double> time(numPts, 0.0);
  std::vector<double> inputSignal(numPts, 0.0);
  std::vector<double> outputSignal(lengthTransformedSignal, 0.0);
  std::vector<double> backSignal(numPts, 0.0);
  
  // these is to make N_LAS_MultiVector objects to test that part of the interface
  // since N_LAS_Vector derrives from MultiVector, one can use this same interface for
  // N_LAS_Vector objects
  std::vector<int> realLbMap(numPts, 0.0);
  std::vector<int> complexLbMap(lengthTransformedSignal, 0.0);
  N_PDS_ParMap parMapForReal( numPts, numPts,realLbMap);
  N_PDS_ParMap parMapForComplex( lengthTransformedSignal, lengthTransformedSignal, complexLbMap);
  
  //N_LAS_MultiVector timeMV( parMapForReal );
  //N_LAS_MultiVector inputSignalMV( parMapForReal );
  //N_LAS_MultiVector outputSignalMV( parMapForComplex );
  //N_LAS_MultiVector backSignalMV( parMapForReal );
  
  N_LAS_Vector timeMV( parMapForReal );
  N_LAS_Vector inputSignalMV( parMapForReal );
  N_LAS_Vector outputSignalMV( parMapForComplex );
  N_LAS_Vector backSignalMV( parMapForReal );
  
  // Fill in the vectors
  for(int i=0; i<numPts; i++)
  {
    time[i] = timeStart + i * deltaTime;
    
    inputSignal[i] = std::sin( 2.0 * M_PI * 1 * time[i] ) + time[i]*time[i];       // period = 1
    //inputSignal[i] = std::sin( 2.0 * M_PI * 10.0 * time[i] )       // frequency = 10.0
    //    + std::cos( 2.0 * M_PI * 75.0 * time[i] )                  // frequency = 75.0
	  //  + std::sin( 2.0 * M_PI * 103.0 * time[i] );                // frequency = 103.0
     
      
//    timeMV[i] = timeStart + i * deltaTime;
//    inputSignalMV[i] = inputSignal[i];
    
    //inputSignal[i] = std::sin( 2.0 * M_PI * time[i] * time[i] );     // frequency = 1.0
  }

  // create the interface  
  N_UTL_FFTInterface<std::vector<double> > myTransform( numPts );
  
  // try taking an FFT
  myTransform.calculateFFT( inputSignal, &outputSignal );
  myTransform.calculateIFT( outputSignal, &backSignal );
  
//  myTransform.calculateFFT( inputSignalMV, &outputSignalMV );
//  myTransform.calculateIFT( outputSignalMV, &backSignalMV );
  
  cout.precision(6);
  std::cout << "time\tInputSignal\tBackSignal" << std::endl;
  for(int i=0; i<numPts; i++)
  {
    std::cout << time[i] << "\t" << inputSignal[i] << "\t" << backSignal[i] << std::endl;
//      << "\t" << timeMV[i] << "\t" << inputSignalMV[i] << "\t" << backSignalMV[i];
     
  }
  std::cout << "Frequency\tReal + Imaginary" << std::endl;
  for (int i=0 ; i<lengthTransformedSignal/2  ; ++i) {
      std::cout << i*freqDelta << "\t" << outputSignal[2*i] << " + "
        << outputSignal[2*i+1] << "i" << std::endl;
  }
/*    if( i <= (numPts/2) )
    {
      std::cout << "\t" << i*freqDelta << "\t" << outputSignal[2*i] << "\t" << outputSignal[2*i+1] ;
//        << "\t" << outputSignalMV[2*i] << "\t" << outputSignalMV[2*i+1];
    }
    else 
    {
      std::cout << "\t" << i*freqDelta << "\t" << outputSignal[2*lengthTransformedSignal - 2*i] 
        << "\t" << outputSignal[2*lengthTransformedSignal - 2*i + 1];
//        << "\t" << outputSignalMV[2*lengthTransformedSignal - 2*i] 
//        << "\t" << outputSignalMV[2*lengthTransformedSignal - 2*i + 1];
    }
    std::cout << std::endl;
*/
  
  // }
 
  
  //
  // Now we'll try doing multiple FFT's at the same time using the N_LAS_BlockVectors
  // Each block will contain all of the problem's unknowns (10 in this case) and the
  // total number of blocks will be the total number of time points.
  /*
  
  int numUnknowns = 4;
  Epetra_SerialComm aSerialComm;
  Epetra_Map baseMap( numUnknowns, 0, aSerialComm );
  Epetra_Map wholeMap( numPts*numUnknowns, numPts*numUnknowns, 0, aSerialComm );
  Epetra_Map complexWholeMap( lengthTransformedSignal*numUnknowns, lengthTransformedSignal*numUnknowns, 0, aSerialComm );
  
  N_LAS_BlockVector realSolution( numPts, baseMap, wholeMap );
  N_LAS_BlockVector complexSolution( lengthTransformedSignal, baseMap, complexWholeMap );
  N_LAS_BlockVector backSolution( numPts, baseMap, wholeMap );
  
  
  // Next, fill up the vectors with some data
  for( int i=0; i<numPts; i++)
  {
    // get the block that we're going to work on
    N_LAS_Vector vec = realSolution.block(i);
    // for this block, time is time[i]
    for( int j=0; j<numUnknowns; j++ )
    {
      double ans = 0.0;
      switch (j) 
      {  
        case 0:
          ans = std::sin( 2.0 * M_PI * 10.0 * time[i] );        // frequency = 10.0
          break;
        case 1: 
          ans = std::cos( 2.0 * M_PI * 75.0 * time[i] );        // frequency = 75.0
          break;
        case 2:
          ans = std::sin( 2.0 * M_PI * 103.0 * time[i] );       // frequency = 103.0
          break;
        case 3:
          ans = std::sin( 2.0 * M_PI * 10.0 * time[i] )         // frequency = 10.0
              + std::cos( 2.0 * M_PI * 75.0 * time[i] )         // frequency = 75.0
              + std::sin( 2.0 * M_PI * 103.0 * time[i] );       // frequency = 103.0
          break;
        case 4: 
          ans = std::sin( 2.0 * M_PI * time[i] * time[i] );
          break;
        case 5:
          ans = std::cos( 2.0 * M_PI * time[i] * time[i] );
          break;
        case 6:
          ans = std::exp( -10.0 * time[i] * time[i] );
          break;
        case 7:
          ans = 1.0 / (std::sqrt( std::fabs( time[i] + 1.0e-10 ) ) );   // don't divide by zero
          break;
        case 8:
          ans = 1.0;
          break;
        case 9:
          ans = time[i];
          break;
        default:
          break;
      }
      vec[j] = ans;
    }
  }
  
  // now set up a new fft interface with strides set up 
  N_UTL_FFTInterface blockTransform( numPts, numUnknowns, numUnknowns );
   
  blockTransform.calculateFFT( realSolution, &complexSolution );
  //blockTransform.calculateIFT( complexSolution, &backSolution );
  
  //char aChar;
  //std::cout << "Hit a key to quit";
  //std::cin >> aChar;
  
  */
  
  return 0;
}
