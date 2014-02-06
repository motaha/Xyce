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
// Filename       : $RCSfile: N_DEV_SourceData.h,v $
//
// Purpose        : Source data containers.  Used by the vsrc and isrc
//                  devices.
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
// Revision Number: $Revision: 1.44.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_SourceData_h
#define Xyce_N_DEV_SourceData_h

// ---------- Standard Includes ----------
#include <vector>
#include <list>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceBlock.h>
#include <N_UTL_BreakPoint.h>

// enum statements:
// device indices:
enum Src_index{
    _SIN_DATA,      // 0
    _EXP_DATA,      // 1
    _PULSE_DATA,    // 2
    _PWL_DATA,      // 3
    _SFFM_DATA,     // 4
    _DC_DATA,    // 5
    _SMOOTH_PULSE_DATA, // 6
    _AC_DATA,     //tmei: 05/02
    //_DISTOF1_DATA,
    //_DISTOF2_DATA,
    _NUM_SRC_DATA       // total number of data types
};


namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : SourceData
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------

class SourceData
{

public:
  SourceData(SolverState & ss1,
             DeviceOptions & do1);
  SourceData(const SourceData &right);

  virtual ~SourceData();

  virtual bool initializeSource ();

  virtual bool updateSource   ();
  virtual bool getBreakPoints (vector<N_UTL_BreakPoint> & breakPointTimes)
  { return true; }

  virtual double getMaxTimeStepSize ();

  virtual void setRealFlag(bool flag) { realFlag_ = true;}

  virtual double period() { return 0.0; }

  double returnSource ();

  bool getResetFlag ()
  { return resetFlag_; }

  string getSourceTypeName ();

  virtual void getParams (double *) {}
  virtual void setParams (double *) {}

#ifdef Xyce_DEBUG_DEVICE
  virtual void printOutParams ();
#endif

  bool getFastTimeScaleFlag() const {
    return fastTimeScaleFlag_;
  }
  void setFastTimeScaleFlag(const bool &fastTimeScaleFlag) {
    fastTimeScaleFlag_ = fastTimeScaleFlag;
  }

protected:
  double getTime_();

private:
  SourceData ();

public:
protected:

  string sourceName_;
  string typeName_;
  string defaultParamName_;

  double time;
  double SourceValue;

  bool initializeFlag_;

  bool resetFlag_;

  SolverState & solState_;
  DeviceOptions & devOptions_;

  bool fastTimeScaleFlag_;

  bool realFlag_;

private:
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
  friend class SourceInstance;
};

//-----------------------------------------------------------------------------
// Class         : SinData
// Purpose       : This class contains data and functions associated with
//                 sinusoidal independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

class SinData : public SourceData
{

public:
  SinData(const SinData &right);
  SinData(const vector<Param> & paramRef,
          SolverState   & ss1,
          DeviceOptions & do1);

  ~SinData();

  bool initializeSource ();
  bool updateSource ();
  void getParams (double *);
  void setParams (double *);

#ifdef Xyce_DEBUG_DEVICE
  void printOutParams();
#endif

  double getMaxTimeStepSize () { return (0.1/FREQ); }

  double period() { return (1.0/FREQ); }

protected:
private:

public:
protected:
private:
  // Data Members for Class Attributes

  double V0;    // Offset (V or A)
  double VA;    // Amplitude (V or A)
  double FREQ;  // Frequency (Hz)
  double TD;    // Delay (seconds)
  double THETA; // Damping factor (1/seconds)
  double PHASE; // Phase (degrees)

  bool   V0given;
  bool   VAgiven;
  bool   FREQgiven;
  bool   TDgiven;
  bool   THETAgiven;
  bool   PHASEgiven;

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};

//-----------------------------------------------------------------------------
// Class         : ExpData
// Purpose       : This class contains data and functions associated with
//                 exponential independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class ExpData : public SourceData
{
public:
  ExpData(const ExpData & right);
  ExpData(const vector<Param> & paramRef,
          SolverState & ss1,
          DeviceOptions & do1);

  ~ExpData();

  bool initializeSource ();
  bool updateSource ();
  void getParams (double *);
  void setParams (double *);

#ifdef Xyce_DEBUG_DEVICE
  void printOutParams ();
#endif

protected:
private:

public:
protected:
private:
  double V1;   // Initial value (V or A)
  double V2;   // Pulsed value (V or A).
  double TD1;  // Rise delay time (seconds).
  double TAU1; // Rise time constant (seconds)
  double TD2;  // Fall delay time (seconds).
  double TAU2; // Fall time constant (seconds)

  bool  V1given;
  bool  V2given;
  bool  TD1given;
  bool  TAU1given;
  bool  TD2given;
  bool  TAU2given;

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};


//-----------------------------------------------------------------------------
// Class         : ACData
// Purpose       : This class contains data and functions associated with
//                 AC independent sources.
// Special Notes :
// Creator       : Ting Mei
// Creation Date :
//-----------------------------------------------------------------------------

class ACData : public SourceData
{

public:
  ACData(const ACData &right);
  ACData(const vector<Param> & paramRef,
         SolverState   & ss1,
         DeviceOptions & do1);

  ~ACData();

  bool updateSource ();
  void getParams (double *);
  void setParams (double *);

#ifdef Xyce_DEBUG_DEVICE
  void printOutParams();
#endif

  void setRealFlag(bool flag) { realFlag_ = flag; }

protected:
private:

public:
protected:
private:

  double ACMAG;    // Amplitude (V or A)
  double ACPHASE; // Phase (degrees)

  bool ACMAGgiven;
  bool ACPHASEgiven;

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};


//-----------------------------------------------------------------------------
// Class         : PulseData
// Purpose       : This class contains data and functions associated with
//                 pulsed independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class PulseData : public SourceData
{
public:
  PulseData(const PulseData  & right);
  PulseData(const vector<Param> & paramRef,
            SolverState   & ss1,
            DeviceOptions & do1);

  ~PulseData();

  bool initializeSource();
  bool updateSource();
  void getParams (double *);
  void setParams (double *);
  bool getBreakPoints(vector<N_UTL_BreakPoint> & breakPointTimes);

#ifdef Xyce_DEBUG_DEVICE
  void printOutParams();
#endif

  double getMaxTimeStepSize ();

  double period() { return PER; }

protected:
private:

public:
  double V1;  // Initial value  (for a voltage source, units are Volts,
  // For a current source, units are Amps)
  double V2;  // Pulsed value. (Volts or Amps)
  double TD;  // Delay time  (seconds)
  double TR;  // Rise time (seconds)
  double TF;  // Fall Time  (seconds)
  double PW;  // Pulse Width (seconds)
  double PER; // Period (seconds)

  bool V1given;
  bool V2given;
  bool TDgiven;
  bool TRgiven;
  bool TFgiven;
  bool PWgiven;
  bool PERgiven;

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};


//-----------------------------------------------------------------------------
// Class         : PWLinData
// Purpose       : This class contains the data and functions associated
//                 with piece-wise linear independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class PWLinData : public SourceData
{
public:
  PWLinData(const PWLinData &right);
  PWLinData(const vector<Param> & paramRef,
            SolverState   & ss1,
            DeviceOptions & do1);

  ~PWLinData();

  bool updateSource();
  bool getBreakPoints( vector<N_UTL_BreakPoint> & breakPointTimes);

#ifdef Xyce_DEBUG_DEVICE
  void printOutParams ();
#endif

protected:
private:

public:
protected:
private:
  // Data Members for Class Attributes
  int NUM; //number of time,voltage pairs
  bool REPEAT; //repeat cycle?
  double REPEATTIME; //start time in cycle for repeat
  double TD; //time delay
  vector< pair<double,double> > TVVEC;  // Array (time,voltage)

  int loc_; //current location in time vector
  double starttime_; //absolute start time of current cycle

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};


//-----------------------------------------------------------------------------
// Class         : SFFMData
// Purpose       : This class contains data and functions associated with
//                 Single-frequency FM independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class SFFMData : public SourceData
{
public:
  SFFMData(const SFFMData   & right);
  SFFMData(const vector<Param> & paramRef,
           SolverState   & ss1,
           DeviceOptions & do1);

  ~SFFMData();

  bool initializeSource ();
  bool updateSource ();
  void getParams (double *);
  void setParams (double *);

#ifdef Xyce_DEBUG_DEVICE
  void printOutParams ();
#endif

protected:
private:

public:
protected:
private:
  double V0;  // Offset. (V or A)
  double VA;  // Amplitude (V or A)
  double FC;  // Carrier frequency (Hz)
  double MDI; // Modulation index
  double FS;  // Signal frequency (Hz)

  bool V0given;
  bool VAgiven;
  bool FCgiven;
  bool MDIgiven;
  bool FSgiven;

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};

//-----------------------------------------------------------------------------
// Class         : ConstData
// Purpose       : This class contains data and functions associated with
//                 const DC sources.  It is not yet implemented.
// Special Notes :
// Creator       : Robert Hoekstra
// Creation Date : 10/4/00
//-----------------------------------------------------------------------------
class ConstData : public SourceData
{
public:
  ConstData(const ConstData   & right);
  ConstData(const vector<Param> & paramRef,
            SolverState   & ss1,
            DeviceOptions & do1);

  ~ConstData();

  bool updateSource();
  void getParams (double *);
  void setParams (double *);

#ifdef Xyce_DEBUG_DEVICE
  void printOutParams();
#endif
protected:
private:

public:
protected:
private:
  double V0;

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};


//-----------------------------------------------------------------------------
// Class         : PulseData
// Purpose       : This class contains data and functions associated with
//                 a smooth pulsed independent source.
// Special Notes :
// Creator       : Richard Schiek
// Creation Date : 4/13/11
//-----------------------------------------------------------------------------
class SmoothPulseData : public SourceData
{
public:
  SmoothPulseData(const SmoothPulseData  & right);
  SmoothPulseData(const vector<Param> & paramRef,
                  SolverState   & ss1,
                  DeviceOptions & do1);

  ~SmoothPulseData();

  bool initializeSource();
  bool updateSource();
  void getParams (double *);
  void setParams (double *);
  bool getBreakPoints(vector<N_UTL_BreakPoint> & breakPointTimes);

#ifdef Xyce_DEBUG_DEVICE
  void printOutParams();
#endif

  double getMaxTimeStepSize ();

  double period() { return PER; }

protected:
private:

public:
  double V1;  // Initial value  (for a voltage source, units are Volts,
  // For a current source, units are Amps)
  double V2;  // Pulsed value. (Volts or Amps)
  double TD;  // Delay time  (seconds)
  double TR;  // Rise time (seconds)
  double TF;  // Fall Time  (seconds)
  double PW;  // Pulse Width (seconds)
  double PER; // Period (seconds)
  double riseScaleFactor_; // scaling of amplitude of rise smoothing function so that
  // it joins the low and high signal appropriately.
  double fallScaleFactor_; // scaling of amplitude of rise smoothing function so that
  // it joins the low and high signal appropriately.
  double functionScaleFactor_; // a scale factor for the smooth function used to fill in the
  // rise and fall segment

  bool V1given;
  bool V2given;
  bool TDgiven;
  bool TRgiven;
  bool TFgiven;
  bool PWgiven;
  bool PERgiven;
  bool functionScaleFactorGiven_;

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::SourceData N_DEV_SourceData;
typedef Xyce::Device::SmoothData N_DEV_SmoothData;
typedef Xyce::Device::SinData N_DEV_SinData;
typedef Xyce::Device::ExpData N_DEV_ExpData;
typedef Xyce::Device::ACData N_DEV_ACData;
typedef Xyce::Device::PulseData N_DEV_PulseData;
typedef Xyce::Device::PWLinData N_DEV_PWLinData;
typedef Xyce::Device::SFFMData N_DEV_SFFMData;\
typedef Xyce::Device::ConstData N_DEV_ConstData;
typedef Xyce::Device::SmoothPulseData N_DEV_SmoothPulseData;

#endif
