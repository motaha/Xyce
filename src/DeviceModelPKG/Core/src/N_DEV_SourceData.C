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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_SourceData.C,v $
//
// Purpose        : This file  contains the member functions of the
//                  N_DEV_SourceData class, which is used by the Vsrc
//                  and ISRC devices.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.124 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_SourceData.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : SourceData::SourceData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
SourceData::SourceData(const SolverState & ss1,
                       const DeviceOptions & do1)
  : SourceValue(0.0),
    initializeFlag_(false),
    resetFlag_(true),
    solState_(ss1),
    time(0.0),
    devOptions_(do1),
    typeName_(""),
    sourceName_(""),
    defaultParamName_("") ,
    fastTimeScaleFlag_(false),
    realFlag_(true)
{

}

//-----------------------------------------------------------------------------
// Function      : SourceData::SourceData
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
SourceData::SourceData(const SourceData & right)
  : SourceValue (right.SourceValue),
    initializeFlag_(right.initializeFlag_),
    resetFlag_(right.resetFlag_),
    solState_(right.solState_),
    time(right.time),
    devOptions_(right.devOptions_),
    typeName_(right.typeName_),
    sourceName_(right.sourceName_),
    defaultParamName_(right.defaultParamName_)
{

}

//-----------------------------------------------------------------------------
// Function      : SourceData::~SourceData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
SourceData::~SourceData()
{

}

//-----------------------------------------------------------------------------
// Function      : SourceData::initializeSource
// Purpose       : Base class initialization function.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
bool SourceData::initializeSource ()
{
  initializeFlag_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::updateSource
// Purpose       : Base class update function - flags an error.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
bool SourceData::updateSource ()
{
  std::string msg;

  msg = "SourceData::updateSource:\n";
  msg += "\tAttempted to update the source for a source type that is not";
  msg += " yet implemented";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);

  return  false;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::returnSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
double SourceData::returnSource()
{
  return SourceValue;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::returnSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
std::string SourceData::getSourceTypeName()
{
  return typeName_;
}



#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : SourceData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void SourceData::printOutParams()
{

}
#endif

//-----------------------------------------------------------------------------
// Function      : SourceData::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/01/02
//-----------------------------------------------------------------------------
double SourceData::getMaxTimeStepSize ()
{
  return devOptions_.defaultMaxTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::getTime_
// Purpose       :
// Special Notes :
// Scope         : protected
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/25/04
//-----------------------------------------------------------------------------
double   SourceData::getTime_()
{
  double tmpTime = 0.0;

//#ifdef Xyce_MPDE
  if (fastTimeScaleFlag_)
    tmpTime = solState_.currFastTime;
  else
//#endif // Xyce_MPDE
    tmpTime = solState_.currTime;

#ifdef Xyce_DEBUG_DEVICE
//#ifdef Xyce_MPDE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "SourceData::getTime.  time = ";
    Xyce::dout() << tmpTime << std::endl;
    Xyce::dout() << "SourceData::getTime.  currFastTime = ";
    Xyce::dout() << solState_.currFastTime<< std::endl;
    Xyce::dout() << "SourceData::getTime.  currTime = ";
    Xyce::dout() << solState_.currTime << std::endl;
  }
//#endif // Xyce_MPDE
#endif
  return tmpTime;
}

// Class SinData

//-----------------------------------------------------------------------------
// Function      : SinData::SinData
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
SinData::SinData(const SinData & right)
  : SourceData(right),
    V0(right.V0),
    VA(right.VA),
    FREQ(right.FREQ),
    TD(right.TD),
    THETA(right.THETA),
    PHASE(right.PHASE),
    V0given(right.V0given),
    VAgiven(right.VAgiven),
    FREQgiven(right.FREQgiven),
    TDgiven(right.TDgiven),
    THETAgiven(right.THETAgiven),
    PHASEgiven(right.PHASEgiven)
{

}

//-----------------------------------------------------------------------------
// Function      : SinData::SinData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//
// From 3f5:  remember to incorporate this later:
//
//-----------------------------------------------------------------------------
#if 0
  // this is from 3f5:
 #define FREQ  (((here->VSRCfunctionOrder >=3) && (*(here->VSRCcoeffs+2)))? \
    (*(here->VSRCcoeffs+2)):(1/ckt->CKTfinalTime))
#endif

SinData::SinData( const std::vector<Param> & paramRef,
                  const SolverState   & ss1,
                  const DeviceOptions & do1)
  : SourceData(ss1,do1),
    V0(0.0),
    VA(0.0),
    FREQ(0.0),
    TD(0.0),
    THETA(0.0),
    PHASE(0.0),
    V0given(false),
    VAgiven(false),
    FREQgiven(false),
    TDgiven(false),
    THETAgiven(false),
    PHASEgiven(false)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V0")    { V0    = iter->getImmutableValue<double>(); V0given = iter->given();}
    if (tmpname == "VA")    { VA    = iter->getImmutableValue<double>(); VAgiven = iter->given();}
    if (tmpname == "FREQ")  { FREQ  = iter->getImmutableValue<double>(); FREQgiven = iter->given();}
    if (tmpname == "TD")    { TD    = iter->getImmutableValue<double>(); TDgiven = iter->given();}
    if (tmpname == "THETA") { THETA = iter->getImmutableValue<double>(); THETAgiven =
                                                      iter->given();}
    if (tmpname == "PHASE") { PHASE = iter->getImmutableValue<double>();
    PHASEgiven = iter->given();}
  }

  typeName_ = "SIN";
  defaultParamName_ = "V0";
}

//-----------------------------------------------------------------------------
// Function      : SinData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
bool SinData::initializeSource ()
{

  // If neccessary, set defaults:
  double tstop = solState_.finalTime;

  if (!FREQgiven)  FREQ  = 1.0/tstop;
  if (!TDgiven)    TD    = 0.0;
  if (!THETAgiven) THETA = 0.0;

  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SinData::~SinData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
SinData::~SinData()
{

}

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : SinData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------

void SinData::printOutParams()

{
  Xyce::dout() << "SinData:\n";
  Xyce::dout() << "V0 = "    << V0 << std::endl;
  Xyce::dout() << "VA = "    << VA << std::endl;
  Xyce::dout() << "FREQ = "  << FREQ << std::endl;
  Xyce::dout() << "TD = "    << TD << std::endl;
  Xyce::dout() << "THETA = " << THETA << std::endl;
  Xyce::dout() << "PHASE = " << PHASE << std::endl;
}
#endif

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : SinData::updateSource
// Purpose       : Update the sinwave source.
// Special Notes :
//
//   V0    - offset  (V or A)
//   VA    - Amplitude  (V or A)
//   FREQ  - frequency in Hz
//   TD    - delay in seconds
//   THETA - damping factor (Hz).
//   PHASE - phase (degrees)
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
bool SinData::updateSource()
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

//#ifdef Xyce_MPDE
  time = getTime_();
//#else
//  time = solState_.currTime;
//#endif // Xyce_MPDE

  time -= TD;
  double mpi = M_PI;
  if (time <= 0)
  {
    //SourceValue = V0;
    SourceValue = V0 + VA * sin (2.0*mpi*(PHASE/360)) ;
  }
  else
  {
    // 2PI to convert from hz to radians/sec
    SourceValue = V0 + VA * sin (2.0*mpi*(FREQ*time + PHASE/360)) *
      exp( -(time*THETA));
  }

  resetFlag_ = false;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : SinData::getParams
// Purpose       : Pass back the sine source params.
// Special Notes : TD and FREQ are interchanged from their normal order
//                 to accomodate the storage layout in the device classes
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void SinData::getParams(double *params)
{
  params[0] = V0;
  params[1] = VA;
  params[2] = TD;
  params[3] = FREQ;
  params[4] = THETA;
  params[5] = PHASE;
  return;
}

//-----------------------------------------------------------------------------
// Function      : SinData::setParams
// Purpose       : Update the sine source params.
// Special Notes : TD and FREQ are interchanged from their normal order
//                 to accomodate the storage layout in the device classes
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void SinData::setParams(double *params)
{
  bool reset=false;
  if (V0 != params[0])
  {
    V0 = params[0];
    reset = true;
  }
  if (VA != params[1])
  {
    VA = params[1];
    reset = true;
  }
  if (TD != params[2])
  {
    TD = params[2];
    reset = true;
  }
  if (FREQ != params[3])
  {
    FREQ = params[3];
    reset = true;
  }
  if (THETA != params[4])
  {
    THETA = params[4];
    reset = true;
  }
  if (PHASE != params[5])
  {
    PHASE = params[5];
    reset = true;
  }
  if (reset)
  {
    updateSource();
  }

  return;
}

// Class ExpData

//-----------------------------------------------------------------------------
// Function      : ExpData::ExpData
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ExpData::ExpData(const ExpData &right)
  : SourceData(right),
    V1   (right.V1),
    V2   (right.V2),
    TD1  (right.TD1),
    TAU1 (right.TAU1),
    TD2  (right.TD2),
    TAU2 (right.TAU2),
    V1given (right.V1given),
    V2given (right.V2given),
    TD1given (right.TD1given),
    TAU1given (right.TAU1given),
    TD2given (right.TD2given),
    TAU2given (right.TAU2given)
{

}

//-----------------------------------------------------------------------------
// Function      : ExpData::ExpData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ExpData::ExpData(const std::vector<Param> & paramRef,
                 const SolverState   & ss1,
                 const DeviceOptions & do1)
  : SourceData(ss1,do1),
    V1   (0.0),
    V2   (0.0),
    TD1  (0.0),
    TAU1 (0.0),
    TD2  (0.0),
    TAU2 (0.0),
    V1given (false),
    V2given (false),
    TD1given (false),
    TAU1given (false),
    TD2given (false),
    TAU2given (false)
{
  // Set the user-defined params:
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V1")    { V1    = iter->getImmutableValue<double>(); V1given   = iter->given();}
    if (tmpname == "V2")    { V2    = iter->getImmutableValue<double>(); V2given   = iter->given();}
    if (tmpname == "TD1")   { TD1   = iter->getImmutableValue<double>(); TD1given  = iter->given();}
    if (tmpname == "TAU1")  { TAU1  = iter->getImmutableValue<double>(); TAU1given = iter->given();}
    if (tmpname == "TD2")   { TD2   = iter->getImmutableValue<double>(); TD2given  = iter->given();}
    if (tmpname == "TAU2")  { TAU2  = iter->getImmutableValue<double>(); TAU2given = iter->given();}
  }

  typeName_ = "EXP";
  defaultParamName_ = "V1";
}

//-----------------------------------------------------------------------------
// Function      : ExpData::~ExpData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ExpData::~ExpData()
{

}

// Additional Declarations

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : ExpData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void ExpData::printOutParams()
{
  Xyce::dout() << "ExpData:\n";

  Xyce::dout() << "V1 = " << V1 << std::endl;
  Xyce::dout() << "V2 = " << V2 << std::endl;

  Xyce::dout() << "TD1 = " << TD1 << std::endl;
  Xyce::dout() << "TAU1 = " << TAU1 << std::endl;

  Xyce::dout() << "TD2 = " << TD2 << std::endl;
  Xyce::dout() << "TAU2 = " << TAU2 << std::endl;

}
#endif

//-----------------------------------------------------------------------------
// Function      : ExpData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
bool ExpData::initializeSource ()
{
  // If neccessary, set defaults:
  double tstep = solState_.startingTimeStep;

  if (!TD1given)  TD1 = 0.0;
  if (!TAU1given) TAU1 = tstep;
  if (!TD2given)  TD2 = TD1 + tstep;
  if (!TAU2given) TAU2 = tstep;

  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ExpData::updateSource
// Purpose       : Updates an exponential source:
// Special Notes :
//
//    V1    -  Initial value (V or A)
//    V2    -  Pulsed value (V or A).
//    TD1   -  Rise delay time (seconds).
//    TAU1  -  Rise time constant (seconds)
//    TD2   -  Fall delay time (seconds).
//    TAU2  -  Fall time constant (seconds)
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------

bool ExpData::updateSource()

{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

//#ifdef Xyce_MPDE
  time = getTime_();
//#else
//  time = solState_.currTime;
//#endif // Xyce_MPDE

  if (time <= TD1)
  {
    SourceValue = V1;
  }
  else if (time <= TD2)
  {
    SourceValue = V1 + (V2-V1)*(1-exp(-(time-TD1)/TAU1));
  }
  else
  {
    SourceValue = V1 + (V2-V1)*(1-exp(-(time-TD1)/TAU1)) +
      (V1-V2)*(1-exp(-(time-TD2)/TAU2)) ;
  }

  resetFlag_ = false;

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : ExpData::getParams
// Purpose       : Pass back the exponential source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void ExpData::getParams(double *params)
{
  params[0] = V1;
  params[1] = V2;
  params[2] = TD1;
  params[3] = TAU1;
  params[4] = TD2;
  params[5] = TAU2;
  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpData::setParams
// Purpose       : Update the exponential source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void ExpData::setParams(double *params)
{
  bool reset=false;
  if (V1 != params[0])
  {
    V1 = params[0];
    reset = true;
  }
  if (V2 != params[1])
  {
    V2 = params[1];
    reset = true;
  }
  if (TD1 != params[2])
  {
    TD1 = params[2];
    reset = true;
  }
  if (TAU1 != params[3])
  {
    TAU1 = params[3];
    reset = true;
  }
  if (TD2 != params[4])
  {
    TD2 = params[4];
    reset = true;
  }
  if (TAU2 != params[5])
  {
    TAU2 = params[5];
    reset = true;
  }
  if (reset)
  {
    updateSource();
  }

  return;
}

// Class PulseData
//-----------------------------------------------------------------------------
// Function      : PulseData::PulseData
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PulseData::PulseData(const PulseData &right)
  : SourceData(right),
    V1  (right.V1),
    V2  (right.V2),
    TD  (right.TD),
    TR  (right.TR),
    TF  (right.TF),
    PW  (right.PW),
    PER (right.PER),
    V1given (right.V1given),
    V2given (right.V2given),
    TDgiven (right.TDgiven),
    TRgiven (right.TRgiven),
    TFgiven (right.TFgiven),
    PWgiven (right.PWgiven),
    PERgiven (right.PERgiven)
{

}

//-----------------------------------------------------------------------------
// Function      : PulseData::PulseData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PulseData::PulseData( const std::vector<Param> & paramRef,
                      const SolverState   & ss1,
                      const DeviceOptions & do1)
  : SourceData (ss1,do1),
    V1  (0.0),
    V2  (0.0),
    TD  (0.0),
    TR  (0.0),
    TF  (0.0),
    PW  (0.0),
    PER (0.0),
    V1given (false),
    V2given (false),
    TDgiven (false),
    TRgiven (false),
    TFgiven (false),
    PWgiven (false),
    PERgiven (false)
{

  // Get the user-defined values:
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V1")  { V1    = iter->getImmutableValue<double>(); V1given = iter->given();}
    if (tmpname == "V2")  { V2    = iter->getImmutableValue<double>(); V2given = iter->given();}
    if (tmpname == "TD")  { TD    = iter->getImmutableValue<double>(); TDgiven = iter->given();}
    if (tmpname == "TR")  { TR    = iter->getImmutableValue<double>(); TRgiven = iter->given();}
    if (tmpname == "TF")  { TF    = iter->getImmutableValue<double>(); TFgiven = iter->given();}
    if (tmpname == "PW")  { PW    = iter->getImmutableValue<double>(); PWgiven = iter->given();}
    if (tmpname == "PER") { PER   = iter->getImmutableValue<double>(); PERgiven = iter->given();}
  }

  typeName_ = "PULSE";
  defaultParamName_ = "V2";
}

//-----------------------------------------------------------------------------
// Function      : PulseData::~PulseData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PulseData::~PulseData()
{

}

// Additional Declarations

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : PulseData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void PulseData::printOutParams()
{
  Xyce::dout() << std::endl;
  Xyce::dout() << "  PulseData::printOutParams\n";
  Xyce::dout() << "  V1  = "    << V1 << std::endl;
  Xyce::dout() << "  V2  = "    << V2 << std::endl;

  Xyce::dout() << "  TD  = "    << TD << std::endl;
  Xyce::dout() << "  TR  = "    << TR << std::endl;
  Xyce::dout() << "  TF  = "    << TF << std::endl;
  Xyce::dout() << "  PW  = "    << PW << std::endl;
  Xyce::dout() << "  PER = "    << PER << std::endl;
  Xyce::dout() << std::endl;

}
#endif

//-----------------------------------------------------------------------------
// Function      : PulseData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
bool PulseData::initializeSource ()
{
  // If neccessary, set the defaults:

  double tstep = solState_.startingTimeStep;
  double tstop = solState_.finalTime;

  if (!TDgiven)  TD  = 0.0;
  if (!TRgiven)  TR  = tstep;
  if (!TFgiven)  TF  = tstep;
  if (!PWgiven)  PW  = tstop;
  if (!PERgiven) PER = tstop;

  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::updateSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
bool PulseData::updateSource()

{
//notused:  double tstep = solState_.startingTimeStep;
//notused:  double tstop = solState_.finalTime;
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  double basetime = 0;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  PulseData::updateSources\n";
    printOutParams();
  }
#endif

//#ifdef Xyce_MPDE
  time = getTime_();
//#else
//  time = solState_.currTime;
//#endif // Xyce_MPDE

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  Time = " << time << std::endl;
  }
#endif

  time -= TD;

  if (time > PER && PER != 0.0)
  {
    // repeating signal - figure out where we are in period
    basetime = PER * floor(time/PER);
    time -= basetime;
  }

  // This section got ugly because of a nasty roundoff bug.
  // Instead of doing "time > X" you need also check that time
  // is not within bptol of X.
  // So the following translation is used:
  // Instead of:                           we do:
  //  time > X                            time>X && fabs(time-x)>bptol
  //  time <= X                           time<X || fabs(time-x)<bptol

  if (time <= 0 || (time > (TR + PW + TF) &&
		    (fabs (time - (TR+PW+TF)) > solState_.bpTol) ) )
  {
    SourceValue = V1;
  }
  else if ((time > TR && fabs(time-TR) > solState_.bpTol)
	   && (time < (TR + PW) || fabs (time-(TR+PW))<solState_.bpTol) )
  {
    SourceValue = V2;
  }
  else if (time > 0 &&
	   (time < TR || fabs(time-TR) < solState_.bpTol))
  {
    if (TR != 0.0)
    {
      SourceValue = V1 + (V2 - V1) * (time) / TR;
    }
    else
    {
      SourceValue = V1;
    }
  }
  else
  { // time > (TR + PW) && <= (TR + PW + TF)
    if (TF != 0.0)
    {
      SourceValue = V2 + (V1 - V2) * (time - (TR + PW)) / TF;
    }
    else
    {
      SourceValue = V2; //      SourceValue = 0.5 * (V1 + V2);
    }
  }

#ifdef  Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  SourceValue = " << SourceValue << std::endl;
  }
#endif

  resetFlag_ = false;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::getParams
// Purpose       : Pass back the pulse source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void PulseData::getParams(double *params)
{
  params[0] = V1;
  params[1] = V2;
  params[2] = TD;
  params[3] = TR;
  params[4] = TF;
  params[5] = PW;
  params[6] = PER;
  return;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::setParams
// Purpose       : Update the pulse source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void PulseData::setParams(double *params)
{
  bool reset=false;
  if (V1 != params[0])
  {
    V1 = params[0];
    reset = true;
  }
  if (V2 != params[1])
  {
    V2 = params[1];
    reset = true;
  }
  if (TD != params[2])
  {
    TD = params[2];
    reset = true;
  }
  if (TR != params[3])
  {
    TR = params[3];
    reset = true;
  }
  if (TF != params[4])
  {
    TF = params[4];
    reset = true;
  }
  if (PW != params[5])
  {
    PW = params[5];
    reset = true;
  }
  if (PER != params[6])
  {
    PER = params[6];
    reset = true;
  }
  if (reset)
  {
    updateSource();
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::getBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
//                 It does not bother to check them in any way, or put them
//                 in order.  It only adds them in.
//
// Special Notes : Like much of this file, this is adapted from spice 3f5.
//                 Some of the stuff in it is a little hokey, and may get
//                 removed or modified later.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/08/01
//-----------------------------------------------------------------------------
bool PulseData::getBreakPoints
  (std::vector<N_UTL_BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  // currPeriodIndex is the integer index representing the period of
  // the current circuit time.
  int currPeriodIndex = 0;

  // subtract out the delay.
  double basetime = 0.0;

//#ifdef Xyce_MPDE
  time = getTime_() - TD;
//#else
//  time = solState_.currTime - TD;
//#endif // Xyce_MPDE

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "  In PulseData::getBreakPoints\n";
    Xyce::dout() << "  time = " << time << std::endl;
    Xyce::dout() << "  TD   = " << TD  << std::endl;
    Xyce::dout() << "  PER  = " << PER << std::endl;
  }
#endif

  // repeating signal - figure out where we are in period
  if(time >= PER)
  {
    if (PER != 0.0)
    {
      currPeriodIndex = (static_cast<int> (floor(time/PER)));
      basetime = PER * (static_cast<double> (currPeriodIndex));
      time -= basetime;
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  time            = " << time << std::endl;
    Xyce::dout() << " basetime         = " << basetime << std::endl;
    Xyce::dout() << "  currPeriodIndex = " << currPeriodIndex << std::endl;
    Xyce::dout() << std::endl;
  }
#endif

   // now that we know which period this is, push_back all breakpoints
   // in this period and the next.  If we are still in the delay, then
   // just use first two periods.

   // current period:
   breakPointTimes.push_back(basetime+TD);
   breakPointTimes.push_back(basetime+TD+TR);
   breakPointTimes.push_back(basetime+TD+TR+PW);
   breakPointTimes.push_back(basetime+TD+TR+PW+TF);

   if (PER != 0.0)
   {
     breakPointTimes.push_back(basetime+TD+PER);

     // next period:
     breakPointTimes.push_back(basetime+TD+PER+TR);
     breakPointTimes.push_back(basetime+TD+PER+TR+PW);
     breakPointTimes.push_back(basetime+TD+PER+TR+PW+TF);
     breakPointTimes.push_back(basetime+TD+PER+PER);
   }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/25/01
//-----------------------------------------------------------------------------
double PulseData::getMaxTimeStepSize ()
{
  double maxTimeStep = devOptions_.defaultMaxTimeStep;

  // check if we are still in the delay or not.
//#ifdef Xyce_MPDE
  time = getTime_();
//#else
//  time = solState_.currTime;
//#endif // Xyce_MPDE

  if (time < TD) maxTimeStep = (0.1*TD );
  else           maxTimeStep = (0.1*PER);

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "\nIn PulseData::getMaxTimeStepSize.  ";
    Xyce::dout() << " maxTimeStep = "<< maxTimeStep;
    Xyce::dout() << "  TD = " << TD << "  PER = " <<PER;
    Xyce::dout() << "  time = "<< time << std::endl;
  }
#endif

  return maxTimeStep;
}

// Class PWLinData

//-----------------------------------------------------------------------------
// Function      : PWLinData::PWLinData
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PWLinData::PWLinData(const PWLinData &right)
  : SourceData(right),
    NUM(right.NUM),
    REPEAT(right.REPEAT),
    REPEATTIME(right.REPEATTIME),
    TD(right.TD),
    TVVEC(right.TVVEC),
    loc_(right.loc_),
    starttime_(right.starttime_)
{
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::PWLinData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PWLinData::PWLinData( const std::vector<Param> & paramRef,
                      const SolverState & ss1,
                      const DeviceOptions & do1)
  : SourceData(ss1,do1),
    NUM(0),
    REPEAT(false),
    REPEATTIME(0.0),
    TD(0.0),
    loc_(0),
    starttime_(0.0)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "NUM")        NUM        = iter->getImmutableValue<int>();
    if (tmpname == "REPEAT")     REPEAT     = iter->getImmutableValue<bool>();
    if (tmpname == "REPEATTIME") REPEATTIME = iter->getImmutableValue<double>();
    if (tmpname == "TD")         TD         = iter->getImmutableValue<double>();

    if ( tmpname == "T" && iter->given() )
    {
      time = iter->getImmutableValue<double>();
      ++iter;

      TVVEC.push_back(std::pair<double,double>(time, iter->getImmutableValue<double>()));
    }
  }

  typeName_ = "PWL";
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::PWLinData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PWLinData::~PWLinData()
{
}

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : PWLinData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void PWLinData::printOutParams()
{
  Xyce::dout() << std::endl;
  Xyce::dout() << "  NUM  = "    << NUM << std::endl;
  Xyce::dout() << "  REPEAT  = "    << REPEAT << std::endl;
  Xyce::dout() << "  REPEATTIME  = "    << REPEATTIME << std::endl;
  Xyce::dout() << "  TD  = "    << TD << std::endl;
  Xyce::dout() << "  loc_  = "    << loc_ << std::endl;
  Xyce::dout() << "  starttime_  = "    << starttime_ << std::endl;

  Xyce::dout() << "  Time    Voltage" << std::endl;
  for( int i = 0; i < NUM; ++i )
    Xyce::dout() << " " << TVVEC[i].first << "  " << TVVEC[i].second << std::endl;

  Xyce::dout() << std::endl;
}
#endif

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : PWLinData::updateSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
bool PWLinData::updateSource()
{

  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "  PWLinData::updateSource\n";
    printOutParams();
  }
#endif

//#ifdef Xyce_MPDE
  time = getTime_();
//#else
//  time = solState_.currTime;
//#endif // Xyce_MPDE

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "  Time = " << time << std::endl;
  }
#endif

  double time1, time2, voltage1, voltage2;

  if( time >= TD )
  {
    time -= TD;

    if( time <= TVVEC[NUM-1].first )
    {
      for( int i = 0; i < NUM; ++i )
        if( time < TVVEC[i].first ) {loc_ = i;break;}

      if( loc_ == 0 )
      {
        time1 = 0.0;
        voltage1 = 0.0;
      }
      else
      {
        time1 = TVVEC[loc_-1].first;
        voltage1 = TVVEC[loc_-1].second;
      }
      time2 = TVVEC[loc_].first;
      voltage2 = TVVEC[loc_].second;

    }
    else if( !REPEAT )
    {
      time1 = 0.0;
      time2 = 1.0;
      voltage1 = voltage2 = TVVEC[NUM-1].second;
    }
    else
    {
      double looptime = TVVEC[NUM-1].first - REPEATTIME;

      time -= TVVEC[NUM-1].first;
      time -= looptime * floor(time / looptime);
      time += REPEATTIME;


      for( int i = 0; i < NUM; ++i )
      {
        if( time < TVVEC[i].first ) {loc_ = i;break;}
      }

      if (time == REPEATTIME)
      {
        time1 = 0.0;
        time2 = 1.0;
        voltage1 = voltage2 = TVVEC[NUM-1].second;
      }
      else
      {
        if( loc_ == 0 )
        {
          time1 = REPEATTIME;
          voltage1 = TVVEC[NUM-1].second;
        }
        else
        {
          time1 = TVVEC[loc_-1].first;
          voltage1 = TVVEC[loc_-1].second;
        }
        time2 = TVVEC[loc_].first;
        voltage2 = TVVEC[loc_].second;
      }

    }

    if( time1 == time2 )
      SourceValue = voltage2;
    else
    {
      double length = time2 - time1;
      SourceValue = ( time2 - time ) * voltage1 / length;
      SourceValue += ( -time1 + time ) * voltage2 / length;
    }

#ifdef Xyce_DEBUG_DEVICE
   if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
   {
      Xyce::dout() << "time: " << time << std::endl;
      Xyce::dout() << "time1: " << time1 << std::endl;
      Xyce::dout() << "time2: " << time2 << std::endl;
      Xyce::dout() << "voltage1: " << voltage1 << std::endl;
      Xyce::dout() << "voltage2: " << voltage2 << std::endl;
      Xyce::dout() << "Src: " << SourceValue << std::endl;
      Xyce::dout() << Xyce::section_divider << std::endl;
    }
#endif

  }
  else
  {
    SourceValue = 0.0;
  }


  resetFlag_ = false;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/08/01
//-----------------------------------------------------------------------------
bool PWLinData::getBreakPoints
  (std::vector<N_UTL_BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  In PWLinData::getBreakPoints\n";
  }
#endif

  if (!initializeFlag_) bsuccess = initializeSource ();

  time = solState_.currTime;

  time -= TD;

  // if it's a repeating signal, figure out what period we are in
  if (REPEAT && time >= TVVEC[NUM - 1].first)
  {
    double loopBaseTime = 0.0;
    double bp_time;

    double looptime = TVVEC[NUM-1].first - REPEATTIME;
    loopBaseTime = looptime * (1.0 + floor((time - TVVEC[NUM - 1].first) / looptime));
#ifdef Xyce_DEBUG_DEVICE
    if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
    {
      Xyce::dout() << "loopBaseTime: " << loopBaseTime << std::endl;
      Xyce::dout() << "floor function: " << floor((time - TVVEC[NUM - 1].first) / (TVVEC[NUM - 1].first - REPEATTIME)) << std::endl;
    }
#endif
    
    for (int i = 0; i < NUM; ++i)
    {
      bp_time = TVVEC[i].first;
      if (bp_time >= REPEATTIME)
      {
        breakPointTimes.push_back(bp_time + loopBaseTime);
#ifdef Xyce_DEBUG_DEVICE
        if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
        {
          Xyce::dout() << "bp_time + loopBaseTime: " << bp_time + loopBaseTime << std::endl;
        }
#endif
      }
    }
  }
  else
  {
    for (int i = 0; i < NUM; ++i)
    {
      double bp_time = TVVEC[i].first;
      breakPointTimes.push_back(bp_time);
#ifdef Xyce_DEBUG_DEVICE
      if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
      {
        Xyce::dout() << "bp_time: " << bp_time << std::endl;
      }
#endif
    }
  }
  // now that we know which period this is, push_back all breakpoints
  // in this period and the next.

  return bsuccess;
}

// Class SFFMData

//-----------------------------------------------------------------------------
// Function      : SFFMData::SFFMData
// Purpose       : copy constructor
// Special Notes : SFFM = spice frequency modulation
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

SFFMData::SFFMData(const SFFMData &right)
  : SourceData(right),
    V0(right.V0),
    VA(right.VA),
    FC(right.FC),
    MDI(right.MDI),
    FS(right.FS),
    V0given  (right.V0given),
    VAgiven  (right.VAgiven),
    FCgiven  (right.FCgiven),
    MDIgiven (right.MDIgiven),
    FSgiven  (right.FSgiven)
{

}

//-----------------------------------------------------------------------------
// Function      : SFFMData::SFFMData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

SFFMData::SFFMData( const std::vector<Param> & paramRef,
                    const SolverState   & ss1,
                    const DeviceOptions & do1)
  : SourceData(ss1,do1),
    V0  (0.0),
    VA  (0.0),
    FC  (0.0),
    MDI (0.0),
    FS  (0.0),
    V0given  (false),
    VAgiven  (false),
    FCgiven  (false),
    MDIgiven (false),
    FSgiven  (false)

{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V0")  { V0    = iter->getImmutableValue<double>(); V0given = iter->given(); }
    if (tmpname == "VA")  { VA    = iter->getImmutableValue<double>(); VAgiven = iter->given(); }
    if (tmpname == "FC")  { FC    = iter->getImmutableValue<double>(); FCgiven = iter->given(); }
    if (tmpname == "MDI") { MDI   = iter->getImmutableValue<double>(); MDIgiven = iter->given(); }
    if (tmpname == "FS")  { FS    = iter->getImmutableValue<double>(); FSgiven = iter->given(); }
  }

  typeName_ = "SFFM";
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::~SFFMData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

SFFMData::~SFFMData()
{

}

// Additional Declarations

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : SFFMData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------

void SFFMData::printOutParams()
{
  Xyce::dout() << "SFFMData:\n";
  Xyce::dout() << "V0 = " << V0 << std::endl;
  Xyce::dout() << "VA = " << VA << std::endl;
  Xyce::dout() << "FC = " << FC << std::endl;
  Xyce::dout() << "MDI = " << MDI << std::endl;
  Xyce::dout() << "FS = " << FS << std::endl;
}
#endif

//-----------------------------------------------------------------------------
// Function      : SFFMData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------

bool SFFMData::initializeSource ()
{
  // If neccessary, set the defaults:
  double tstop = solState_.finalTime;

  if (!FCgiven) FC = 1.0/tstop;
  if (!FSgiven) FS = 1.0/tstop;

  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::updateSource
// Purpose       :
// Special Notes :
//
//   V0    - offset  (V or A)
//   VA    - Amplitude  (V or A)
//   FC    - carrier frequency
//   MDI   - modulation index
//   FS    - signal frequency
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------

bool SFFMData::updateSource()
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

//#ifdef Xyce_MPDE
  time = getTime_();
//#else
//  time = solState_.currTime;
//#endif // Xyce_MPDE

  double mpi = M_PI;
  SourceValue = V0 + VA * sin((2 * mpi * FC * time) +
                              MDI * sin (2 * mpi * FS * time));

  resetFlag_ = false;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::getParams
// Purpose       : Pass back the SFFM source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void SFFMData::getParams(double *params)
{
  params[0] = V0;
  params[1] = VA;
  params[2] = FC;
  params[3] = MDI;
  params[4] = FS;
  return;
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::setParams
// Purpose       : Update the SFFM source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void SFFMData::setParams(double *params)
{
  bool reset=false;
  if (V0 != params[0])
  {
    V0 = params[0];
    reset = true;
  }
  if (VA != params[1])
  {
    VA = params[1];
    reset = true;
  }
  if (FC != params[2])
  {
    FC = params[2];
    reset = true;
  }
  if (MDI != params[3])
  {
    MDI = params[3];
    reset = true;
  }
  if (FS != params[4])
  {
    FS = params[4];
    reset = true;
  }
  if (reset)
    updateSource();

  return;
}

// Class ACData

//-----------------------------------------------------------------------------
// Function      : ACData::ACData
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------

ACData::ACData(const ACData & right)
  : SourceData(right),
    ACMAG(right.ACMAG),
    ACPHASE(right.ACPHASE),
    ACMAGgiven(right.ACMAGgiven),
    ACPHASEgiven(right.ACPHASEgiven)
{

}


//-----------------------------------------------------------------------------
// Function      : ACData::ACData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------

ACData::ACData( const std::vector<Param> & paramRef,
                const SolverState   & ss1,
                const DeviceOptions & do1)
  : SourceData(ss1,do1),
    ACMAG(1.0),
    ACPHASE(0.0),
    ACMAGgiven(false),
    ACPHASEgiven(false)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "ACMAG")    { ACMAG    = iter->getImmutableValue<double>(); ACMAGgiven = iter->given();}
    if (tmpname == "ACPHASE") { ACPHASE = iter->getImmutableValue<double>(); ACPHASEgiven = iter->given();}
  }

  typeName_ = "AC";
  defaultParamName_ = "ACMAG";
}

//-----------------------------------------------------------------------------
// Function      : ACData::~ACData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

ACData::~ACData()
{

}

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : ACData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------

void ACData::printOutParams()

{
  Xyce::dout() << "ACData:\n";
  Xyce::dout() << "ACMAG = " << ACMAG << std::endl;
  Xyce::dout() << "ACPHASE = " << ACPHASE << std::endl;
}
#endif

//-----------------------------------------------------------------------------
// Function      : ACData::updateSource
// Purpose       : Update the sinwave source.
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------

bool ACData::updateSource()
{
  bool bsuccess = true;

  double mpi = M_PI;

  if (!initializeFlag_) bsuccess = initializeSource ();

    if (realFlag_)
    {
      SourceValue = ACMAG * cos(2.0*mpi*ACPHASE/360);
    }
    else
    {
      SourceValue = ACMAG * sin(2.0*mpi*ACPHASE/360);
    }

#ifdef  Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "  SourceValue = " << SourceValue << std::endl;
  }
#endif

  resetFlag_ = false;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ACData::getParams
// Purpose       : Pass back the AC source params.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date :
//-----------------------------------------------------------------------------
void ACData::getParams(double *params)
{
  params[0] = ACMAG;
  params[1] = ACPHASE;
  return;
}

//-----------------------------------------------------------------------------
// Function      : ACData::setParams
// Purpose       : Update the AC source params.
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei
// Creation Date :
//-----------------------------------------------------------------------------
void ACData::setParams(double *params)
{
  bool reset=false;
  if (ACMAG!=  params[0])
  {
    ACMAG = params[0];
    reset = true;
  }
  if ( ACPHASE != params[1])
  {
    ACPHASE = params[1];
    reset = true;
  }
  if (reset)
    updateSource();

  return;
}


// Class ConstData

//-----------------------------------------------------------------------------
// Function      : ConstData::ConstData
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------
ConstData::ConstData(const ConstData & right)
  : SourceData(right),
    V0(right.V0)
{

}

//-----------------------------------------------------------------------------
// Function      : ConstData::ConstData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------

ConstData::ConstData( const std::vector<Param> & paramRef,
                      const SolverState   & ss1,
                      const DeviceOptions & do1)
  : SourceData(ss1,do1),
    V0(0.0)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();
    if (tmpname == "DCV0")
    {
      V0 = iter->getImmutableValue<double>();
    }
  }

  typeName_ = "CONST";
  defaultParamName_ = "DCV0";
//  SourceValue = V0; // updateSource function is a no-op, essentially.
}

//-----------------------------------------------------------------------------
// Function      : ConstData::~ConstData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------
ConstData::~ConstData()
{

}

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : ConstData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------

void ConstData::printOutParams()

{
  Xyce::dout() << "ConstData:\n";
  Xyce::dout() << "V0: " << V0 << std::endl;
}
#endif

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : ConstData::updateSource
// Purpose       : Update the const source.
// Special Notes : ERK: this is now a no-op, as the source value is set in
//                 the constructor.  The value for this source never changes,
//                 so there isn't any point is re-setting the same value
//                 over and over.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------
bool ConstData::updateSource()
{

  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  SourceValue = V0;

  resetFlag_ = false;

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : ConstData::getParams
// Purpose       : Pass back the const source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/13/05
//-----------------------------------------------------------------------------
void ConstData::getParams(double *params)
{
  params[0] = V0;
  return;
}

//-----------------------------------------------------------------------------
// Function      : ConstData::setParams
// Purpose       : Update the const source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/13/05
//-----------------------------------------------------------------------------
void ConstData::setParams(double *params)
{
  if (V0 != params[0])
  {
    V0 = params[0];
//    SourceValue = V0;
    updateSource();
  }
  return;
}


// Class SmoothPulseData
//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::SmoothPulseData
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------

SmoothPulseData::SmoothPulseData(const SmoothPulseData &right)
  : SourceData(right),
    V1  (right.V1),
    V2  (right.V2),
    TD  (right.TD),
    TR  (right.TR),
    TF  (right.TF),
    PW  (right.PW),
    PER (right.PER),
    riseScaleFactor_(right.riseScaleFactor_),
    fallScaleFactor_(right.fallScaleFactor_),
    functionScaleFactor_(right.functionScaleFactor_),
    V1given (right.V1given),
    V2given (right.V2given),
    TDgiven (right.TDgiven),
    TRgiven (right.TRgiven),
    TFgiven (right.TFgiven),
    PWgiven (right.PWgiven),
    PERgiven (right.PERgiven),
    functionScaleFactorGiven_(right.functionScaleFactorGiven_)
{

}

//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::SmoothPulseData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------

SmoothPulseData::SmoothPulseData( const std::vector<Param> & paramRef,
                                  const SolverState   & ss1,
                                  const DeviceOptions & do1)
  : SourceData (ss1,do1),
    V1  (0.0),
    V2  (0.0),
    TD  (0.0),
    TR  (0.0),
    TF  (0.0),
    PW  (0.0),
    PER (0.0),
    riseScaleFactor_(0.0),
    fallScaleFactor_(0.0),
    functionScaleFactor_(20.0),
    V1given (false),
    V2given (false),
    TDgiven (false),
    TRgiven (false),
    TFgiven (false),
    PWgiven (false),
    PERgiven (false),
    functionScaleFactorGiven_(false)
{

  // Get the user-defined values:
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V1")  { V1    = iter->getImmutableValue<double>(); V1given = iter->given();}
    if (tmpname == "V2")  { V2    = iter->getImmutableValue<double>(); V2given = iter->given();}
    if (tmpname == "TD")  { TD    = iter->getImmutableValue<double>(); TDgiven = iter->given();}
    if (tmpname == "TR")  { TR    = iter->getImmutableValue<double>(); TRgiven = iter->given();}
    if (tmpname == "TF")  { TF    = iter->getImmutableValue<double>(); TFgiven = iter->given();}
    if (tmpname == "PW")  { PW    = iter->getImmutableValue<double>(); PWgiven = iter->given();}
    if (tmpname == "PER") { PER   = iter->getImmutableValue<double>(); PERgiven = iter->given();}
    if (tmpname == "SF")  { functionScaleFactor_ = iter->getImmutableValue<double>(); functionScaleFactorGiven_ = iter->given();}
  }

  typeName_ = "SMOOTHPULSE";
  defaultParamName_ = "V2";
}

//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::~SmoothPulseData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------

SmoothPulseData::~SmoothPulseData()
{

}

// Additional Declarations

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------

void SmoothPulseData::printOutParams()
{

  Xyce::dout() << std::endl;
  Xyce::dout() << "  SmoothPulseData::printOutParams\n";
  Xyce::dout() << "  V1  = "    << V1 << std::endl;
  Xyce::dout() << "  V2  = "    << V2 << std::endl;

  Xyce::dout() << "  TD  = "    << TD << std::endl;
  Xyce::dout() << "  TR  = "    << TR << std::endl;
  Xyce::dout() << "  TF  = "    << TF << std::endl;
  Xyce::dout() << "  PW  = "    << PW << std::endl;
  Xyce::dout() << "  PER = "    << PER << std::endl;
  Xyce::dout() << "  SF  = "    << functionScaleFactor_ << std::endl;
  Xyce::dout() << std::endl;

}
#endif

//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------

bool SmoothPulseData::initializeSource ()
{

  // If neccessary, set the defaults:

  double tstep = solState_.startingTimeStep;
  double tstop = solState_.finalTime;

  if (!TDgiven)  TD  = 0.0;
  if (!TRgiven)  TR  = tstep;
  if (!TFgiven)  TF  = tstep;
  if (!PWgiven)  PW  = tstop;
  if (!PERgiven) PER = tstop;
  // scale the amplitude of the response so it better joins the low and hi value
  // riseFallScaleFactor_ = M_PI;
  riseScaleFactor_ = 2.0*fabs( atan(M_PI * functionScaleFactor_ * (0.5*TR) / TR) );
  fallScaleFactor_ = 2.0*fabs( atan(M_PI * functionScaleFactor_ * (0.5*TF) / TF) );
  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::updateSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------

bool SmoothPulseData::updateSource()

{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  double basetime = 0;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  SmoothPulseData::updateSources\n";
    printOutParams();
  }
#endif

  time = getTime_();

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  Time = " << time << std::endl;
  }
#endif

  time -= TD;

  if (time > PER && PER != 0.0)
  {
    // repeating signal - figure out where we are in period
    basetime = PER * floor(time/PER);
    time -= basetime;
  }

  // This section got ugly because of a nasty roundoff bug.
  // Instead of doing "time > X" you need also check that time
  // is not within bptol of X.
  // So the following translation is used:
  // Instead of:                           we do:
  //  time > X                            time>X && fabs(time-x)>bptol
  //  time <= X                           time<X || fabs(time-x)<bptol

  if (time <= 0 || (time > (TR + PW + TF) &&
		    (fabs (time - (TR+PW+TF)) > solState_.bpTol) ) )
  {
    SourceValue = V1;
  }
  else if ((time > TR && fabs(time-TR) > solState_.bpTol)
	   && (time < (TR + PW) || fabs (time-(TR+PW))<solState_.bpTol) )
  {
    SourceValue = V2;
  }
  else if (time > 0 &&
	   (time < TR || fabs(time-TR) < solState_.bpTol))
  {
    if (TR != 0.0)
      SourceValue = V1 + (V2 - V1) *
        (((atan(M_PI * functionScaleFactor_ * (time-0.5*TR) / TR))/riseScaleFactor_) + 0.5);
    else
      SourceValue = V1;
  }
  else
  { // time > (TR + PW) && <= (TR + PW + TF)
    if (TF != 0.0)
      SourceValue = V2 + (V1 - V2) *
        (((atan(M_PI * functionScaleFactor_ * (time - (TR + PW + 0.5*TF)) / TF))/fallScaleFactor_) + 0.5);
    else
      SourceValue = V2;
  }

#ifdef  Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  SourceValue = " << SourceValue << std::endl;
  }
#endif

  resetFlag_ = false;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::getParams
// Purpose       : Pass back the pulse source params.
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------
void SmoothPulseData::getParams(double *params)
{
  params[0] = V1;
  params[1] = V2;
  params[2] = TD;
  params[3] = TR;
  params[4] = TF;
  params[5] = PW;
  params[6] = PER;
  return;
}

//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::setParams
// Purpose       : Update the pulse source params.
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------
void SmoothPulseData::setParams(double *params)
{
  bool reset=false;
  if (V1 != params[0])
  {
    V1 = params[0];
    reset = true;
  }
  if (V2 != params[1])
  {
    V2 = params[1];
    reset = true;
  }
  if (TD != params[2])
  {
    TD = params[2];
    reset = true;
  }
  if (TR != params[3])
  {
    TR = params[3];
    reset = true;
  }
  if (TF != params[4])
  {
    TF = params[4];
    reset = true;
  }
  if (PW != params[5])
  {
    PW = params[5];
    reset = true;
  }
  if (PER != params[6])
  {
    PER = params[6];
    reset = true;
  }
  if (reset)
    updateSource();

  return;
}

//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::getBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
//                 It does not bother to check them in any way, or put them
//                 in order.  It only adds them in.
//
// Special Notes : Like much of this file, this is adapted from spice 3f5.
//                 Some of the stuff in it is a little hokey, and may get
//                 removed or modified later.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------
bool SmoothPulseData::getBreakPoints
  (std::vector<N_UTL_BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : SmoothPulseData::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/13/11
//-----------------------------------------------------------------------------
double SmoothPulseData::getMaxTimeStepSize ()
{
  double maxTimeStep = devOptions_.defaultMaxTimeStep;

  // check if we are still in the delay or not.
  time = getTime_();

  if (time < TD) maxTimeStep = (0.1*TD );
  else           maxTimeStep = (0.1*PER);

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    Xyce::dout() << "\nIn SmoothPulseData::getMaxTimeStepSize.  ";
    Xyce::dout() << " maxTimeStep = "<< maxTimeStep;
    Xyce::dout() << "  TD = " << TD << "  PER = " <<PER;
    Xyce::dout() << "  time = "<< time << std::endl;
  }
#endif

  return maxTimeStep;
}

} // namespace Device
} // namespace Xyce
