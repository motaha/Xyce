//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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
// Revision Number: $Revision: 1.28 $
// Revision Date  : $Date: 2014/02/24 23:49:20 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <utility>

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

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Manager::Manager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
Manager::Manager(N_IO_OutputMgr &outputManager)
  : outputManager_(outputManager)
{}

//-----------------------------------------------------------------------------
// Function      : Manager::Manager
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
Manager::~Manager()
{
  for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
    delete (*it);
}

void
Manager::fixupMeasureParameters() 
{
  for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
    (*it)->fixupMeasureParameters();
}


//-----------------------------------------------------------------------------
// Function      : Manager::addMeasure
// Purpose       : Entery point when .measure lines are pased in the netlist
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool Manager::addMeasure(const Util::OptionBlock & measureBlock )
{
  // based on what's in the option block passed in, we
  // create the needed measure instance
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "In Manager::addMeasure" << std::endl;
  Xyce::dout() << "measureLine passed it was: " << std::endl << measureBlock << std::endl;
#endif

  // here's the base data we should pull from the option block while
  std::string type;

  if( measureBlock.tagExists( "TYPE" ) )
  {
    type = measureBlock.getTagValueAsString("TYPE");
  }
  else
  {
    // this shouldn't happen, but catch it if does
    Report::UserError0() << "Missing TYPE in Manager";
  }

  N_IO_MeasureBase  * theMeasureObject = 0;
  // have enough info to make the correct measure class
  if( type=="TRIG" || type=="TARG" )
  {
    theMeasureObject = new N_IO_MeasureRiseFallDelay( measureBlock, outputManager_ );
  }
  else if( type=="AVG" )
  {
    theMeasureObject = new N_IO_MeasureAverage( measureBlock, outputManager_ );
  }
  else if( type=="MAX" )
  {
    theMeasureObject = new N_IO_MeasureMax( measureBlock, outputManager_ );
  }
  else if( type=="MIN" )
  {
    theMeasureObject = new N_IO_MeasureMin( measureBlock, outputManager_ );
  }
  else if( type=="PP" )
  {
    theMeasureObject = new N_IO_MeasurePeakToPeak( measureBlock, outputManager_ );
  }
  else if( type=="RMS" )
  {
    theMeasureObject = new N_IO_MeasureRMS( measureBlock, outputManager_ );
  }
  else if( type=="FREQ" )
  {
    theMeasureObject = new N_IO_MeasureFrequency( measureBlock, outputManager_ );
  }
  else if( type=="DUTY" )
  {
    theMeasureObject = new N_IO_MeasureDuty( measureBlock, outputManager_ );
  }
  else if( type=="ON_TIME" )
  {
    theMeasureObject = new N_IO_MeasureOnTime( measureBlock, outputManager_ );
  }
  else if( type=="OFF_TIME" )
  {
    theMeasureObject = new N_IO_MeasureOffTime( measureBlock, outputManager_ );
  }
  else if( type=="FIND" || type=="WHEN" )
  {
    theMeasureObject = new N_IO_MeasureFindWhen( measureBlock, outputManager_ );
  }
  else if( type=="PARAM" || type=="EQN"  )
  {
    theMeasureObject = new N_IO_MeasureEquationEvaluation( measureBlock, outputManager_ );
  }
  else if( type=="DERIVATIVE" || type=="DERIV" )
  {
    theMeasureObject = new N_IO_MeasureDerivativeEvaluation( measureBlock, outputManager_ );
  }
  else if( type=="INTEGRAL" || type=="INTEG" )
  {
    theMeasureObject = new N_IO_MeasureIntegralEvaluation( measureBlock, outputManager_ );
  }
  else if( type=="ERROR" )
  {
    theMeasureObject = new N_IO_MeasureRelativeError( measureBlock, outputManager_ );
  }
  else if( type=="FOUR" )
  {
    theMeasureObject = new N_IO_MeasureFourier( measureBlock, outputManager_ );
  }
  else
  {
    // unknown type issue warning.
    Xyce::Report::UserWarning() << "Unknown MEASURE type requested, \"" << type << "\".  This measure will be ignored";
  }

  // if the measure object is supported, then add it to the active and all lists
  if (theMeasureObject && theMeasureObject->typeSupported_ )
  {
    allMeasuresList_.push_back( theMeasureObject );
    activeMeasuresList_.push_back( theMeasureObject );
  }
  else
    delete theMeasureObject;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::updateTranMeasures
// Purpose       : Called during the simulation to update the measure objects
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void Manager::updateTranMeasures( const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  // loop over active masure objects and get them to update themselves.
  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateTran(circuitTime, solnVec, stateVec, storeVec);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), std::mem_fun(&Measure::Base::finishedCalculation)),
                            activeMeasuresList_.end());
}


//-----------------------------------------------------------------------------
// Function      : Manager::updateDcMeasures
// Purpose       : Called during the simulation to update the measure objects
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void Manager::updateDcMeasures( const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateDC(dcParamsVec, solnVec, stateVec, storeVec);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), std::mem_fun(&Measure::Base::finishedCalculation)),
                            activeMeasuresList_.end());
}

//-----------------------------------------------------------------------------
// Function      : Manager::updateAcMeasures
// Purpose       : Called during the simulation to update the measure objects
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/21/2014
//-----------------------------------------------------------------------------
void Manager::updateAcMeasures( const double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateAC(frequency, real_solution_vector, imaginary_solution_vector);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), std::mem_fun(&Measure::Base::finishedCalculation)),
                            activeMeasuresList_.end()); }

//-----------------------------------------------------------------------------
// Function      : Manager::outputResults
// Purpose       : Output measure results at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
std::ostream &Manager::outputResults( std::ostream& outputStream, bool printHeader )
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
    for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
    {
      (*it)->printMeasureResult( outputStream );
    }
  }
  
  return outputStream;
}


//-----------------------------------------------------------------------------
// Function      : Manager::getMeasureValue
// Purpose       : Get the value of a .measure test
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/29/2010
//-----------------------------------------------------------------------------
void Manager::getMeasureValue(const std::string &name, double &value, bool &found) const
{
  found = false;

  for (MeasurementVector::const_iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
  {
    if ((*it)->name_ == name) {
      value = (*it)->getMeasureResult();
      found = true;
      return;
    }
  }
}

const Base *Manager::find(const std::string &name) const
{
  for (MeasurementVector::const_iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
    if ((*it)->name_ == name)
      return *it;

  return 0;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce
