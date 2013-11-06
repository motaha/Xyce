//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_IO_MeasureManager.C,v $
// Purpose       : This file contains the functions to manage measure objects
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.13.2.3 $
// Revision Date  : $Date: 2013/10/03 17:23:42 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_IO_MeasureManager.h>
#include <N_IO_MeasureRiseFallDelay.h>
#include <N_IO_MeasureAverage.h>
#include <N_IO_MeasureMax.h>
#include <N_IO_MeasureMin.h>
#include <N_IO_MeasurePeakToPeak.h>
#include <N_IO_MeasureRMS.h>
#include <N_IO_MeasureFrequency.h>
#include <N_IO_MeasureDuty.h>
#include <N_IO_MeasureOnTime.h>
#include <N_IO_MeasureOffTime.h>
#include <N_IO_MeasureFindWhen.h>
#include <N_IO_MeasureEquationEvaluation.h>
#include <N_IO_MeasureDerivativeEvaluation.h>
#include <N_IO_MeasureIntegralEvaluation.h>
#include <N_IO_MeasureRelativeError.h>
#include <N_IO_MeasureFourier.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>

//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureManager::N_IO_MeasureManager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
N_IO_MeasureManager::N_IO_MeasureManager(N_IO_OutputMgr &outputManager)
  : outputManager_(outputManager)
{}

//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureManager::N_IO_MeasureManager
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
N_IO_MeasureManager::~N_IO_MeasureManager()
{

}

//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureManager::addMeasure
// Purpose       : Entery point when .measure lines are pased in the netlist
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool N_IO_MeasureManager::addMeasure(const N_UTL_OptionBlock & measureBlock )
{
  // based on what's in the option block passed in, we
  // create the needed measure instance
#ifdef Xyce_DEBUG_IO
  std::cout << "In N_IO_MeasureManager::addMeasure" << std::endl;
  std::cout << "measureLine passed it was: " << std::endl << measureBlock << std::endl;
#endif

  // here's the base data we should pull from the option block while
  string type;

  if( measureBlock.tagExists( "TYPE" ) )
  {
    type = measureBlock.getTagValueAsString("TYPE");
  }
  else
  {
    // this shouldn't happen, but catch it if does
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Missing TYPE in N_IO_MeasureManager");
  }

  RCP< N_IO_MeasureBase > theMeasureObject;
  // have enough info to make the correct measure class
  if( type=="TRIG" || type=="TARG" )
  {
    theMeasureObject = rcp( new N_IO_MeasureRiseFallDelay( measureBlock, outputManager_ ));
  }
  else if( type=="AVG" )
  {
    theMeasureObject = rcp( new N_IO_MeasureAverage( measureBlock, outputManager_ ));
  }
  else if( type=="MAX" )
  {
    theMeasureObject = rcp( new N_IO_MeasureMax( measureBlock, outputManager_ ));
  }
  else if( type=="MIN" )
  {
    theMeasureObject = rcp( new N_IO_MeasureMin( measureBlock, outputManager_ ));
  }
  else if( type=="PP" )
  {
    theMeasureObject = rcp( new N_IO_MeasurePeakToPeak( measureBlock, outputManager_ ));
  }
  else if( type=="RMS" )
  {
    theMeasureObject = rcp( new N_IO_MeasureRMS( measureBlock, outputManager_ ));
  }
  else if( type=="FREQ" )
  {
    theMeasureObject = rcp( new N_IO_MeasureFrequency( measureBlock, outputManager_ ));
  }
  else if( type=="DUTY" )
  {
    theMeasureObject = rcp( new N_IO_MeasureDuty( measureBlock, outputManager_ ));
  }
  else if( type=="ON_TIME" )
  {
    theMeasureObject = rcp( new N_IO_MeasureOnTime( measureBlock, outputManager_ ));
  }
  else if( type=="OFF_TIME" )
  {
    theMeasureObject = rcp( new N_IO_MeasureOffTime( measureBlock, outputManager_ ));
  }
  else if( type=="FIND" || type=="WHEN" )
  {
    theMeasureObject = rcp( new N_IO_MeasureFindWhen( measureBlock, outputManager_ ));
  }
  else if( type=="PARAM" )
  {
    theMeasureObject = rcp( new N_IO_MeasureEquationEvaluation( measureBlock, outputManager_ ));
  }
  else if( type=="DERIVATIVE" || type=="DERIV" )
  {
    theMeasureObject = rcp( new N_IO_MeasureDerivativeEvaluation( measureBlock, outputManager_ ));
  }
  else if( type=="INTEGRAL" || type=="INTEG" )
  {
    theMeasureObject = rcp( new N_IO_MeasureIntegralEvaluation( measureBlock, outputManager_ ));
  }
  else if( type=="ERROR" )
  {
    theMeasureObject = rcp( new N_IO_MeasureRelativeError( measureBlock, outputManager_ ));
  }
  else if( type=="FOUR" )
  {
    theMeasureObject = rcp( new N_IO_MeasureFourier( measureBlock, outputManager_ ));
  }
  else
  {
    // unknown type issue warning.
    string msg = "Unknown MEASURE type requested, \"" + type + "\" This measure will be ignored.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING, msg);
  }

  // if the measure object is supported, then add it to the active and all lists
  // if it's not supported, then the ref count pointer will let it go on exiting
  // this routine.
  if( theMeasureObject->typeSupported_ )
  {
    allMeasuresList_.push_back( theMeasureObject );
    activeMeasuresList_.push_back( theMeasureObject );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureManager::updateDcMeasures
// Purpose       : Called during the simulation to update the measure objects
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureManager::updateTranMeasures( const double circuitTime, RCP< N_LAS_Vector > solnVecRCP)
{
  // loop over active masure objects and get them to update themselves.
  list<RCP<N_IO_MeasureBase> >::iterator currentMeasure = activeMeasuresList_.begin();
  list<RCP<N_IO_MeasureBase> >::iterator endMeasure = activeMeasuresList_.end();
  while( currentMeasure != endMeasure )
  {
    (*currentMeasure)->updateTran(circuitTime, solnVecRCP);
    if( (*currentMeasure)->finishedCalculation() )
    {
      // remove this one from the active list
      list<RCP<N_IO_MeasureBase> >::iterator compleatedMeasure = currentMeasure;
      currentMeasure++;
      activeMeasuresList_.erase( compleatedMeasure );
    }
    else
    {
      currentMeasure++;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureManager::updateDcMeasures
// Purpose       : Called during the simulation to update the measure objects
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureManager::updateDcMeasures( const vector<N_ANP_SweepParam> & dcParamsVec, RCP< N_LAS_Vector > solnVecRCP)
{
  // loop over active measrue objects and get them to update themselves.
  list<RCP<N_IO_MeasureBase> >::iterator currentMeasure = activeMeasuresList_.begin();
  list<RCP<N_IO_MeasureBase> >::iterator endMeasure = activeMeasuresList_.end();
  while( currentMeasure != endMeasure )
  {
    (*currentMeasure)->updateDC(dcParamsVec, solnVecRCP);
    if( (*currentMeasure)->finishedCalculation() )
    {
      // remove this one from the active list
      list<RCP<N_IO_MeasureBase> >::iterator compleatedMeasure = currentMeasure;
      currentMeasure++;
      activeMeasuresList_.erase( compleatedMeasure );
    }
    else
    {
      currentMeasure++;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureManager::outputResults
// Purpose       : Output measure results at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void N_IO_MeasureManager::outputResults( std::ostream& outputStream, bool printHeader )
{

  if (!allMeasuresList_.empty())
  {
    if (printHeader)
    {
      outputStream << std::endl
        << " ***** Measure Functions ***** " << std::endl
        << std::endl;
    }

    // loop over active measure objects and get the results
    list<RCP<N_IO_MeasureBase> >::iterator currentMeasure = allMeasuresList_.begin();
    list<RCP<N_IO_MeasureBase> >::iterator endMeasure = allMeasuresList_.end();
    while( currentMeasure != endMeasure )
    {
      (*currentMeasure)->printMeasureResult( outputStream );
      currentMeasure++;
    }
  }

}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureManager::getMeasureValue
// Purpose       : Get the value of a .measure test
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/29/2010
//-----------------------------------------------------------------------------
void N_IO_MeasureManager::getMeasureValue (const string &name, double &value, bool &found)
{
  found = false;
  if (allMeasuresList_.empty())
    return;
  list<RCP<N_IO_MeasureBase> >::iterator currentMeasure = allMeasuresList_.begin();
  list<RCP<N_IO_MeasureBase> >::iterator endMeasure = allMeasuresList_.end();
  while( currentMeasure != endMeasure )
  {
    if ((*currentMeasure)->name_ == name) {
      found = true;
      value = (*currentMeasure)->getMeasureResult();
      return;
    }
    ++currentMeasure;
  }
  return;
}
