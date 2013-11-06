//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: bug1574.C,v $
// Purpose       : This file contains functions to call Xyce as a library
//                 and test the functionality of the Xyce library interface.
// Special Notes :
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/15/200901/15/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.4 $
// Revision Date  : $Date: 2009/05/05 19:18:20 $
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

// ---------- Standard Includes ----------
#include <iostream>
#include <ctime>
#include <cmath>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_CIR_Xyce.h>

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/02/09
//-----------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
  
  // make a Xyce object
  RefCountPtr< N_CIR_Xyce > xycePtr_ = rcp( new N_CIR_Xyce() );
  
  // Initialize Xyce.
  if ( ! xycePtr_->initialize(argc, argv) )
  {
    std::cerr << "Failed in N_CIR_Xyce::initialize" << std::endl;
    //    return;
    exit(-1);
  }
 
#if 0 
  // Set up the ADC's
  // map of <ADCNAME, <parametername,parametervalue>>
  std::map<string, std::map<string,double> > ADCParamsMap_;
  xycePtr_->getADCMap(ADCParamsMap_);
  
  // Setting up ADC width
  // must construct a map of ADC name/bit vector widths
  map<string,int> ADCWidthMap;
  std::map<string, std::map<string,double> > ::iterator currentADC = ADCParamsMap_.begin();
  std::map<string, std::map<string,double> > ::iterator endADC = ADCParamsMap_.end();
  const int defaultWidth = 8;
  while( currentADC != endADC )
  {
    ADCWidthMap[currentADC->first] = defaultWidth;
    currentADC++;
  }
  xycePtr_->setADCWidths(ADCWidthMap);
#endif
  
  // Set up time voltage pairs
  map< string, vector< pair<double,double> > > xyceTimeVoltageUpdateMap;
  xycePtr_->getTimeVoltagePairs(xyceTimeVoltageUpdateMap);
  map< string, vector< pair<double,double> > >::iterator tvCurrentPair = xyceTimeVoltageUpdateMap.begin();
  map< string, vector< pair<double,double> > >::iterator tvEndPair = xyceTimeVoltageUpdateMap.end();
  
  map< string, vector< pair<double,double> >* > timeVoltageUpdateMap;

  // Get the names of all DACs in the analog circuit.
  vector< string > dacNames;
  if ( ! xycePtr_->getDACDeviceNames(dacNames) )
  {
    std::cerr << "Failed to get the DAC device names" << std::endl;
    //    return;
    exit(-1);
  }

  cout.width(14);cout.precision(8);cout.setf(ios::scientific);
  
  // simulate the circuit
  bool stepStatus=false;
  double oldTime = 0.0;
  double actualTime = 0.0;
  double finalTime = xycePtr_->getFinalTime ();
    
    //0.25e-3;  // there is still a bug in Xyce if this is called with 0.0 and later resumed. RLS

  double initialTime = 0.0;
  double dt = (finalTime-initialTime)/100.0;
  double curr_computed_dt = dt;

  // For each DAC being simulated...
  vector<string>::iterator nameIter = dacNames.begin();
  vector<string>::iterator endIter = dacNames.end();

  // this is just the amount of time ahead of the current time where 
  // we should place the next breakpoint.  Making this larger makes 
  // this test faster because the simulation is paused fewer times.
  double dacValue = 1.0;
  double time=0.0;
  double xycetime=0.0;
  for( ; nameIter != endIter; nameIter++ )
  {
    vector< pair<double,double> >* timeVoltageUpdatesPtr(0);
    timeVoltageUpdatesPtr = new vector< pair<double,double> >;

    // update this dac
    for (int i=0;i<10;++i)
    {
      double scale = time/finalTime;
      dacValue = sin(M_PI*scale);

      timeVoltageUpdatesPtr->push_back( make_pair( time, dacValue ) );
      cout.width(14);cout.precision(8);cout.setf(ios::scientific);
      cout << "bug1574:134:" << *nameIter<< ": time,value="<<time<<","<<dacValue<<endl;

      time += 10*dt;
    }
    cout << "size of the vector: " << timeVoltageUpdatesPtr->size() << endl;
    timeVoltageUpdateMap[*nameIter] = timeVoltageUpdatesPtr;
  }
  cout << "double-checking the time,value pairs:" << endl;
  nameIter = dacNames.begin();
  for( ; nameIter != endIter; nameIter++ )
  {
    vector< pair<double,double> >* tVUPtr(0);

    tVUPtr = timeVoltageUpdateMap[*nameIter];
    for (int i=0;i<tVUPtr->size();++i)
    {
      double time = (*tVUPtr)[i].first;
      double value = (*tVUPtr)[i].second;
      cout << "bug1574:151:" << *nameIter << ": time,value="<<time<<","<<value<<endl;
    }
  }

  bool bs1 = true;
  cout.width(14);cout.precision(8);cout.setf(ios::scientific);
  cout << "Calling dcop provisionalStep.  actualTime = " << actualTime <<endl;
  cout << "curr_computed_dt = " << curr_computed_dt << endl;
  stepStatus = xycePtr_->provisionalStep(0.0, curr_computed_dt, timeVoltageUpdateMap);

  //(double maxTimeStep, 
   //double &timeStep, 
   //map< string, vector< pair<double,double> > > & timeVoltageUpdateMap)

  xycePtr_->acceptProvisionalStep();
  actualTime = xycePtr_->getTime();

  double maxdt = 1.0e-7;

  do 
  {
    cout.width(14);cout.precision(8);cout.setf(ios::scientific);
    cout << "Calling dcop provisionalStep.  actualTime = " << actualTime <<endl;

    stepStatus = xycePtr_->provisionalStep(maxdt, curr_computed_dt, timeVoltageUpdateMap);

    cout.width(14);cout.precision(8);cout.setf(ios::scientific);
    cout << "actual dt = "<<curr_computed_dt<<endl;

    if (curr_computed_dt > maxdt)
    {
      xycePtr_->rejectProvisionalStep();
    }
    else
    {
      xycePtr_->acceptProvisionalStep();
    }

    actualTime = xycePtr_->getTime ();

  } while ( actualTime < finalTime );
    
  if ( ! xycePtr_->finalize() )
  {
    cerr << "Failed to N_CIR_Xyce::finalize\n";
    exit(1);
  }
  
  exit(0);
}
