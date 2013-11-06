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
// Filename       : $RCSfile: N_DEV_DeviceOptions.h,v $
//
// Purpose        :
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
// Revision Number: $Revision: 1.58.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DeviceOptions_h
#define Xyce_N_DEV_DeviceOptions_h

#include <string>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_UTL_Param.h>

// ----------   Forward Declarations   ----------
class N_IO_CmdParse;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceOptions
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DeviceOptions
{
  public:
    DeviceOptions(N_IO_CmdParse & cp);
    ~DeviceOptions();

    bool setupDefaultOptions ();
    bool applyCmdLineOptions ();
    bool registerOptions (const N_UTL_OptionBlock & OB);

    void setBlockAnalysisFlag(bool flagVal);

    // Assignment operator:
    DeviceOptions & operator=(DeviceOptions const & rhs);

#ifdef Xyce_DEBUG_DEVICE
    friend ostream & operator<<(ostream & os, const DeviceOptions & devOp);
#endif

  protected:

  private:

  public:

    // some general MOS parameters:
    double defad;  // MOS drain diffusion area.
    double defas;  // MOS source diffusion area.
    double defl;   // MOS channel length.
    double defw;   // MOS channel width.

    double abstol; // absolute current error tolerance.
    double reltol; // relative current error tolerance.
    double chgtol; // absolute charge error tolerance.

    double gmin;      // minimum allowed conductance.
    double gmin_orig; // this is needed for gmin-homotopy.
    double gmin_init; // this is needed for gmin-homotopy.
    double gmin_scalar; // this is needed for gmin-homotopy.

    double gmax;   // maximum allowed conductance.

    double testJac_relTol; // reltol for num. jacobian diagnostic
    double testJac_absTol; // abstol for num. jacobian diagnostic.
    double testJac_SqrtEta;// dx = numJacSqrtEta * (1.0 + fabs(soln[i]));
    double deviceSens_dp;  // similar to eta, but for numerical device sensitivities

    double tnom;   // nominal temperature for device params.
    N_UTL_Param temp;   // operating temperature of ckt.

    double scale_src; // scaling for source loads

    bool numericalJacobianFlag;
    bool testJacobianFlag;
    int testJacStartStep;
    int testJacStopStep;
    bool testJacWarn;
    bool testJacDeviceNameGiven;
    string testJacDeviceName;
    bool voltageLimiterFlag;
    int lambertWFlag;

    bool newMeyerFlag;

    double icMultiplier;

    double defaultMaxTimeStep;

    // mosfet homotopy:
    double vdsScaleMin;
    double vgstConst;
    int numGainScaleBlocks;
    bool staggerGainScale;
    bool randomizeVgstConst;
    double length0;
    double width0;
    double tox0;
    double minRes;
    double minCap;
    double exp_order;

    // tolerance on resistance below which it will be treated as zero
    double zeroResistanceTol;
    bool checkForZeroResistance;
    bool detailedDeviceCounts;

#ifdef Xyce_DEBUG_DEVICE
    int    sensDebugLevel;
    int    debugLevel;
    int    debugMinTimestep;
    int    debugMaxTimestep;
    double debugMinTime;
    double debugMaxTime;
#endif
    int    verboseLevel;

    bool blockAnalysisFlag;       // This indicates an MPDE/HB run.  This is true during both IC and MPDE/HB phase.
    // It is toggled by the presence of the .mpde analysis statement in the netlist.

    bool newExcessPhase;
    bool defaultNewExcessPhase;   // default is true for MPDE, false for non-MPDE.

    double excessPhaseScalar1;
    double excessPhaseScalar2;

    unsigned int randomSeed;      // seed for random number generator used by some devices.
    // note: each device gets its own random number generator so
    // it must initialize thing correctly. (See N_DEV_Synapse3 for an
    // example)

    bool tryToCompact;            // Try to compact past history for LTRA device(s).

    bool calculateAllLeadCurrents;   // class level flag to have the device manager
                                     // configure all devices to load lead current.
                                     // data in store and storeQvec.

    N_IO_CmdParse & commandLine;
};

inline void DeviceOptions::setBlockAnalysisFlag(bool flagVal)
{
  blockAnalysisFlag = flagVal;

  // tscoffe/tmei 07/30/08:  Note, if we call this with "false", then the
  // newExcessPhase and defaultNewExcessPhase flags will not go back to their
  // "user-set" values, they will be set to false.
  newExcessPhase = flagVal;
  defaultNewExcessPhase = flagVal;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceOptions N_DEV_DeviceOptions;

#endif // Xyce_N_DEV_DeviceOptions_h
