//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, 2013, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: XyceLibTest.C,v $
// Purpose       : This file contains functions to call Xyce as a library
//                 and test the functionality of the Xyce library interface.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 02/14/2008
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.6 $
// Revision Date  : $Date: 2013/09/26 20:14:46 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

// 
// To make this a better test for specific bugs in the bug repository,
// I'll have the code check for specific functionality and emit an error
// only if that funciton fails.  The bugs checked are:
//
// BUG 1466  -- simulateUntil() fails to return upon ADC output change
// BUG 1499  -- Use raw name for ADC and DAC
// BUG 1500  -- Breakpoint needed at end of ramp up/down DACC
// BUG 1628  -- Initial ADC values needed from getTimeVoltagePairs()
//
// Not checked: 
// BUG 1467  -- ADC / DAC are not compatible with newdae
//
// From here I can't check if Xyce is running in new dae mode
// however, the netlist used for this test ADC_DACtest.cir
// does not force the time integrator to use the old dae load.



// ---------- Standard Includes ----------
#include <iostream>
#include <ctime>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_CIR_Xyce.h>

int main( int argc, char * argv[] )
{
  // flags for the bug's we're checking
  bool pass_bug_1466 = false;
  bool pass_bug_1499 = false;
  bool pass_bug_1500 = false;
  bool pass_bug_1628 = false;
  
  // this is just the amount of time ahead of the current time where 
  // we should place the next breakpoint.  Making this larger makes 
  // this test faster because the simulation is paused fewer times.
  const double deltaTimeForPause = 1.0e-5;
  
  // make a Xyce object
  RefCountPtr< N_CIR_Xyce > xycePtr_ = rcp( new N_CIR_Xyce() );
  
  // Initialize Xyce.
  if ( ! xycePtr_->initialize(argc, argv) )
  {
    std::cerr << "Failed in N_CIR_Xyce::initialize" << std::endl;
    //    return;
    exit(-1);
  }
  
  // Get the names of all DACs in the analog circuit.
  vector< string > dacNames;
  if ( ! xycePtr_->getDACDeviceNames(dacNames) )
  {
    std::cerr << "Failed to get the DAC device names" << std::endl;
    //    return;
    exit(-1);
  }
  
  // now check that the names start with "Y%" this idicates that we have 
  // the raw device names.
  pass_bug_1499 = true;
  vector<string>::iterator nameIter = dacNames.begin();
  vector<string>::iterator endIter = dacNames.end();
  for( ; nameIter != endIter; nameIter++ )
  {
   if( nameIter->find("Y%") == string::npos )
   {
     pass_bug_1499 = false;
   }
  }
  
  // Setting up ADC's 
  // set up the ADC's
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
  
  // Set up time voltage pairs
  map< string, vector< pair<double,double> > > xyceTimeVoltageUpdateMap;
  xycePtr_->getTimeVoltagePairs(xyceTimeVoltageUpdateMap);
  map< string, vector< pair<double,double> > >::iterator tvCurrentPair = xyceTimeVoltageUpdateMap.begin();
  map< string, vector< pair<double,double> > >::iterator tvEndPair = xyceTimeVoltageUpdateMap.end();
  
  //
  // for debugging print out the set of time voltage pairs
  //
  //   Xyce::dout() << "xyceTimeVoltageUpdateMap = ";
  //   for ( ; tvCurrentPair != tvEndPair; ++tvCurrentPair )
  //   {
  //     Xyce::dout() << "\"" << tvCurrentPair->first << "\": ";
  //     vector< pair<double,double> >::iterator currentPair = (tvCurrentPair->second).begin();
  //     vector< pair<double,double> >::iterator endPair = (tvCurrentPair->second).end();
  //     for( ;currentPair != endPair; ++currentPair)
  //     {
  //       Xyce::dout() << "( " << currentPair->first << ", " << currentPair->second << " ) ";
  //     }
  //   }
  //   Xyce::dout() << std::endl;
  
  map< string, vector< pair<double,double> >* > timeVoltageUpdateMap;

  // For each DAC being simulated...
  nameIter = dacNames.begin();
  double dacValue = 1.0;
  for( ; nameIter != endIter; nameIter++ )
  {
    // update this dac
    vector< pair<double,double> >* timeVoltageUpdatesPtr = new vector< pair<double,double> >;
    timeVoltageUpdatesPtr->push_back( make_pair( deltaTimeForPause, dacValue ) );
    timeVoltageUpdateMap[*nameIter] = timeVoltageUpdatesPtr;
  }
  
  //
  // for debugging print out the new set of time voltage pairs
  //
  //   Xyce::dout() << "New time voltage pairs " << std::endl;
  //   map< string, vector< pair<double,double> >* >::iterator ntvCurrentPair = timeVoltageUpdateMap.begin();
  //   map< string, vector< pair<double,double> >* >::iterator ntvEndPair = timeVoltageUpdateMap.end();
  //   
  //   Xyce::dout() << "timeVoltageUpdateMap = ";
  //   for ( ; ntvCurrentPair != ntvEndPair; ++ntvCurrentPair )
  //   {
  //     Xyce::dout() << "\"" << ntvCurrentPair->first << "\": ";
  //     vector< pair<double,double> >::iterator currentPair = ntvCurrentPair->second->begin();
  //     vector< pair<double,double> >::iterator endPair = ntvCurrentPair->second->end();
  //     for( ;currentPair != endPair; ++currentPair)
  //     {
  //       Xyce::dout() << "( " << currentPair->first << ", " << currentPair->second << " ) ";
  //     }
  //   }
  //   Xyce::dout() << std::endl;

  // Update the time-voltage pairs in Xyce
  if (timeVoltageUpdateMap.size() != 0)
  {
    //Xyce::dout() << "Calling updateTimeVoltagePairs " << std::endl;
    xycePtr_->updateTimeVoltagePairs(timeVoltageUpdateMap);
  }
  
  // simulate the circuit
  bool simStatus=false;
  double actualTime = 0.0;

  double requstedTime = 1.0e-4;  // there is still a bug in Xyce if this is called with 0.0 and later resumed. RLS
  double nextDacPauseTime = deltaTimeForPause;
  do
  {

    //Xyce::dout() << "Calling simulateUntil( " << requstedTime << ", " << actualTime << ")" << std::endl;
    simStatus = xycePtr_->simulateUntil(requstedTime, actualTime);
    //Xyce::dout() << "Simulation Paused at " << actualTime << std::endl;
    if( !simStatus )
    {
      Xyce::dout() << "Simulation exited early. " << std::endl;
      break;
    } 
    
    if( !pass_bug_1628 )
    {
      // Need to get the time voltage pairs and see if we get the dc op value too
      xyceTimeVoltageUpdateMap.clear();
      xycePtr_->getTimeVoltagePairs(xyceTimeVoltageUpdateMap);
      tvCurrentPair = xyceTimeVoltageUpdateMap.begin();
      tvEndPair = xyceTimeVoltageUpdateMap.end();
      
      //Xyce::dout() << "BUG 1628 Check: xyceTimeVoltageUpdateMap = ";
      for ( ; tvCurrentPair != tvEndPair; ++tvCurrentPair )
      {
        //Xyce::dout() << "\"" << tvCurrentPair->first << "\": ";
        vector< pair<double,double> >::iterator currentPair = (tvCurrentPair->second).begin();
        vector< pair<double,double> >::iterator endPair = (tvCurrentPair->second).end();
        for( ;currentPair != endPair; ++currentPair)
        {
          //Xyce::dout() << "( " << currentPair->first << ", " << currentPair->second << " ) ";
          if( currentPair->first == 0.0 )
          {
            // we have a dc op value so but 1628 is considered a pass
            pass_bug_1628 = true;
          }
        }
      }
    }
    
    if( actualTime > nextDacPauseTime )
    {
      //Xyce::dout() << "actualTime > nextDacPauseTime " << actualTime << " > " << nextDacPauseTime << std::endl;
      if( !pass_bug_1500 )
      {
        // On bug 1500 the ADC/DACs were not setting breakpoints that Xyce could use
        // to pause the simulation.  If we get this this part of the test run, then
        // the ADC/DACs are setting breakpoints, so this bug is considered passed
        pass_bug_1500 = true;
      }
      
      // try and add a new state change on the DAC 
      if( dacValue > 0.0 )
      {
        dacValue = 0.0;
      }
      else
      {
        dacValue = 1.0;
        if( !pass_bug_1466 )
        {
          // On bug 1466 a change in voltage of a DAC was not setting a breakpiont
          // if we make it here then the default of DAC value of 1.0 set before this
          // loop has been set to zero which caused a break point to be set and brought
          // back here to puch the DAC back to 1.0.  Thus, the dac's are functioning
          // properly and bug 1466 is passing
          pass_bug_1466 = true;
        }
        
      }
      nameIter = dacNames.begin();
      for( ; nameIter != endIter; nameIter++ )
      {
        // update this dac
        vector< pair<double,double> >* timeVoltageUpdatesPtr = new vector< pair<double,double> >;
        timeVoltageUpdatesPtr->push_back( make_pair( actualTime+deltaTimeForPause, dacValue ) );
        timeVoltageUpdateMap[*nameIter] = timeVoltageUpdatesPtr;
      }
    
      // Update the time-voltage pairs.
      if (timeVoltageUpdateMap.size() != 0)
      {

        //Xyce::dout() << "Calling updateTimeVoltagePairs " << std::endl;
        xycePtr_->updateTimeVoltagePairs(timeVoltageUpdateMap);
      }
      nextDacPauseTime += deltaTimeForPause;
    }
    
    requstedTime = 1.0e-4;
  } while ( actualTime < requstedTime );
    
  if ( ! xycePtr_->finalize() )
  {
    cerr << "Failed to N_CIR_Xyce::finalize\n";
    //    return;
    exit(1);
  }
  
  //
  // print out status:
  //
  Xyce::dout() << "BUG  1466 = ";
  if( pass_bug_1466 ) 
    Xyce::dout() << "pass" << std::endl;
  else
    Xyce::dout() << "FAIL" << std::endl;
    
    Xyce::dout() << "BUG  1499 = ";
  if( pass_bug_1499 ) 
    Xyce::dout() << "pass" << std::endl;
  else
    Xyce::dout() << "FAIL" << std::endl;
  
  Xyce::dout() << "BUG  1500 = ";
  if( pass_bug_1500 ) 
    Xyce::dout() << "pass" << std::endl;
  else
    Xyce::dout() << "FAIL" << std::endl;
    
  Xyce::dout() << "BUG  1628 = ";
  if( pass_bug_1628 ) 
    Xyce::dout() << "pass" << std::endl;
  else
    Xyce::dout() << "FAIL" << std::endl;
    
  int returnStatus = 0;
  if( !(pass_bug_1466 && pass_bug_1499 && pass_bug_1500 && pass_bug_1628) )
    returnStatus = -1;


  return returnStatus;
}
