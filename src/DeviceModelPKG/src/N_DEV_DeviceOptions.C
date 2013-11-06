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
// Filename       : $RCSfile: N_DEV_DeviceOptions.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/06/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.76.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

#ifdef HAVE_CLIMITS
#include <climits>
#else
#include <limits.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_UTL_OptionBlock.h>
#include <N_ERH_ErrorMgr.h>

#include <N_IO_CmdParse.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::DeviceOptions
// Purpose       : constructor
//
// Special Notes : The are initialized to default values, but are reset
//                 in the setupDefaultOptions function.  Confusing I know,
//                 but I had a reason for doing this, I think.
//
//                 Consider  setupDefaultOptions function to be the
//                 ultimate authority.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/01/02
//-----------------------------------------------------------------------------
DeviceOptions::DeviceOptions(N_IO_CmdParse & cp) :

  defad (0.0e+0),  // MOS drain diffusion area.
  defas (0.0e+0),  // MOS source diffusion area.
  defl  (1.0e-4),  // MOS channel length.
  defw  (1.0e-4),  // MOS channel width.
  abstol(1.0e-12), // absolute current error tol.
  reltol(1.0e-4),  // relative current error tol.
  chgtol(1.0e-12), // absolute charge error tol.
  gmin  (1.0e-12), // minimum allowed conductance.
  gmin_orig (1.0e-12), // minimum allowed conductance, final
  gmin_init (1.0e-02), // minimum allowed conductance, initial
  gmin_scalar(1.0e10),
  gmax  (1.0e20),   // maximum allowed conductance.
  testJac_relTol(0.01),   // reltol for num. jacobian diagnostic
  testJac_absTol(1.0e-8), // abstol for num. jacobian diagnostic.
  testJac_SqrtEta(1.0e-8), // dx = numJacSqrtEta * (1.0 + fabs(soln[i]));
  deviceSens_dp(1.0e-8),  // similar to eta, but for numerical device sensitivities
  tnom  (CONSTREFTEMP),
  scale_src (0.0),
  numericalJacobianFlag (false),
  testJacobianFlag (false),
  testJacStartStep(0),
  testJacStopStep(N_UTL_MachineDependentParams::IntMax()),
  testJacWarn (false),
  testJacDeviceName(""),
  testJacDeviceNameGiven( false ),
  voltageLimiterFlag (true),
  lambertWFlag (0),
  icMultiplier (10000.0),
  defaultMaxTimeStep (1.0e99),
  vdsScaleMin(0.3),
  vgstConst(4.5),
  numGainScaleBlocks(1),
  staggerGainScale(false),
  randomizeVgstConst(false),
  length0(5.0e-6), // used in mosfet "size" homotopy
  width0(200.0e-6), // used in mosfet "size" homotopy
  tox0(6.0e-8), // used in mosfet "size" homotopy
  minRes(0.0),
  minCap(0.0),
  exp_order(100.0),
  zeroResistanceTol(1.0e-100),
  checkForZeroResistance(true),
  detailedDeviceCounts(false),
#ifdef Xyce_DEBUG_DEVICE
  sensDebugLevel(0),
  debugLevel (1),
  debugMinTimestep (0),
  debugMaxTimestep (N_UTL_MachineDependentParams::IntMax()), // later, this should be MAX_INT
  debugMinTime (0),
  debugMaxTime (N_UTL_MachineDependentParams::DoubleMax()),
#endif
  verboseLevel (0),
  blockAnalysisFlag (false),
#ifndef Xyce_NEW_EXCESS_PHASE
  newExcessPhase    (false),
  defaultNewExcessPhase    (false),
#else
  newExcessPhase    (true),
  defaultNewExcessPhase    (true),
#endif
  excessPhaseScalar1 (1.0),
  excessPhaseScalar2 (1.0),
  randomSeed (0),
  tryToCompact (false),
  calculateAllLeadCurrents (false),
  commandLine (cp),
  newMeyerFlag (false)
{
  temp = N_UTL_Param("TEMP",CONSTREFTEMP);
  setupDefaultOptions ();
  applyCmdLineOptions ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::~DeviceOptions
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/01/02
//-----------------------------------------------------------------------------
DeviceOptions::~DeviceOptions ()
{
}

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::setupDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/13/00
//-----------------------------------------------------------------------------
bool DeviceOptions::setupDefaultOptions ()
{

  defad = 0.0e+0;  // MOS drain diffusion area.
  defas = 0.0e+0;  // MOS source diffusion area.
  defl  = 1.0e-4;  // MOS channel length.
  defw  = 1.0e-4;  // MOS channel width.
  abstol= 1.0e-12; // absolute current error tol.
  reltol = 1.0e-4; // relative current error tol.
  chgtol= 1.0e-12; // absolute charge error tol.
  gmin  = 1.0e-12; // minimum allowed conductance.
  gmin_orig  = 1.0e-12; // minimum allowed conductance, reference
  gmin_init = 1.0e-02;  // minimum allowed conductance, initial
  gmin_scalar = 1.0e10;
  gmax  = 1.0e20;   // maximum allowed conductance.

  testJac_relTol = 0.01;   // reltol for num. jacobian diagnostic
  testJac_absTol = 1.0e-8; // abstol for num. jacobian diagnostic.
  testJac_SqrtEta = 1.0e-8; // dx = numJacSqrtEta * (1.0 + fabs(soln[i]));
  deviceSens_dp = 1.0e-8;  // similar to eta, but for numerical device sensitivities

  tnom  = CONSTREFTEMP; // nominal temp. for device params.
  temp  = N_UTL_Param("TEMP",CONSTREFTEMP); // operating temp. of ckt.

  scale_src = 0;

  numericalJacobianFlag = false;
  testJacobianFlag      = false;
  testJacStartStep      = 0;
  testJacStopStep       = N_UTL_MachineDependentParams::IntMax();
  testJacWarn           = false;
  voltageLimiterFlag    = true;
  lambertWFlag = 0;

  newMeyerFlag = false;

  icMultiplier = 10000.0;

  defaultMaxTimeStep = 1.0e99;

#ifdef Xyce_DEBUG_DEVICE
  sensDebugLevel = 0;
  debugLevel = 1;
  debugMinTimestep=0;
  debugMaxTimestep=N_UTL_MachineDependentParams::IntMax();
  debugMinTime=0;
  debugMaxTime=N_UTL_MachineDependentParams::DoubleMax();
#endif
  verboseLevel = 0;

  newExcessPhase     = defaultNewExcessPhase; // true if MPDE, false otherwise.
  excessPhaseScalar1  = 1.0;
  excessPhaseScalar2  = 1.0;

  if( commandLine.getArgumentValue( "-dva" ) == "off" )
  {
    string tmp(
      "  Warning: -dva is no longer a recognized option.\n");
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0,  tmp);
  }

  if( commandLine.getArgumentValue( "-dma" ) == "off" )
  {
    string tmp(
      "  Warning: -dma is no longer a recognized option.\n");
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0,  tmp);
  }

  return true;

}

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::registerOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
bool DeviceOptions::registerOptions(const N_UTL_OptionBlock & OB)
{
  list<N_UTL_Param>::const_iterator iter = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator end   = OB.getParams().end();

  for ( ; iter !=  end;  ++iter)
  {
    string tag(iter->lTag());

    if (tag == "defad")         defad     = iter->dVal();
    else if (tag == "defas")    defas     = iter->dVal();
    else if (tag == "defl")     defl      = iter->dVal();
    else if (tag == "defw")     defw      = iter->dVal();
    else if (tag == "abstol")   abstol    = iter->dVal();
    else if (tag == "reltol")   reltol    = iter->dVal();
    else if (tag == "chgtol")   chgtol    = iter->dVal();
    else if (tag == "gmin")     gmin      = iter->dVal();
    else if (tag == "gminscalar") gmin_scalar = iter->dVal();
    else if (tag == "gmax")     gmax      = iter->dVal();
    else if (tag == "tjreltol") testJac_relTol = iter->dVal();
    else if (tag == "tjabstol") testJac_absTol = iter->dVal();
    else if (tag == "tjsqrteta") testJac_SqrtEta = iter->dVal();
    else if (tag == "sensdp") deviceSens_dp = iter->dVal();
    else if (tag == "tnom")     tnom      = iter->dVal()+CONSTCtoK;
    else if (tag == "temp")
    {
      temp = *iter;
      if ( !iter->isTimeDependent() ) temp.setVal(iter->dVal()+CONSTCtoK);
      else { temp.setVal( iter->sVal() ); temp.setTimeDependent( true ); }
    }
    else if (tag == "scalesrc") scale_src = iter->dVal();
    else if (tag == "numjac")
      numericalJacobianFlag = static_cast<bool>(iter->iVal());
    else if (tag == "testjac")
      testJacobianFlag = static_cast<bool>(iter->iVal());
    else if (tag == "testjacstartstep")
      testJacStartStep = iter->iVal();
    else if (tag == "testjacstopstep")
      testJacStopStep = iter->iVal();
    else if (tag == "testjacwarn")
      testJacWarn = static_cast<bool> (iter->iVal());
    else if (tag == "testjacdevicename")
    {
      testJacDeviceName = iter->sVal();
      testJacDeviceNameGiven = true;
    }
    else if (tag == "voltlim")
      voltageLimiterFlag    = static_cast<bool>(iter->iVal());
    else if (tag == "lambertw")  lambertWFlag = static_cast<int>(iter->iVal());
    else if (tag == "icfac" )   icMultiplier = iter->dVal();
    else if (tag == "maxtimestep" ) defaultMaxTimeStep = iter->dVal();

    else if (tag == "vdsscalemin" ) vdsScaleMin = iter->dVal();
    else if (tag == "vgstconst" ) vgstConst  = iter->dVal();
    else if (tag == "numgainscaleblocks" )
      numGainScaleBlocks = static_cast<int>(iter->iVal());
    else if (tag == "staggergainscale")
      staggerGainScale = static_cast<bool>(iter->iVal());
    else if (tag == "randomizevgstconst")
      randomizeVgstConst = static_cast<bool>(iter->iVal());
    else if (tag == "length0" ) length0 = iter->dVal();
    else if (tag == "width0" )  width0  = iter->dVal();
    else if (tag == "tox0" )    tox0    = iter->dVal();
    else if (tag == "minres" )  minRes  = iter->dVal();
    else if (tag == "mincap" )  minCap  = iter->dVal();
    else if (tag == "newmeyer" )
      newMeyerFlag = iter->bVal();
    else if (tag == "sensdebuglevel")
    {
#ifdef Xyce_DEBUG_DEVICE
      sensDebugLevel   = (iter->iVal());
#endif
    }
    else if (tag == "debuglevel")
    {
#ifdef Xyce_DEBUG_DEVICE
      debugLevel       = (iter->iVal());
#endif
    }
    else if (tag == "verboselevel")
    {
      verboseLevel     = (iter->iVal());
    }
    else if (tag == "debugmintimestep")
    {
#ifdef Xyce_DEBUG_DEVICE
      debugMinTimestep = (iter->iVal());
#endif
    }
    else if (tag == "debugmaxtimestep")
    {
#ifdef Xyce_DEBUG_DEVICE
      debugMaxTimestep = (iter->iVal());
#endif
    }
    else if (tag == "debugmintime")
    {
#ifdef Xyce_DEBUG_DEVICE
      debugMinTime     = (iter->dVal());
#endif
    }
    else if (tag == "debugmaxtime")
    {
#ifdef Xyce_DEBUG_DEVICE
      debugMaxTime     = (iter->dVal());
#endif
    }
    else if (tag == "newexcessphase")
    {
      newExcessPhase = static_cast<bool> (iter->iVal());
    }
    else if (tag == "excessphasescalar1")
    {
      excessPhaseScalar1 = (iter->dVal());
    }
    else if (tag == "excessphasescalar2")
    {
      excessPhaseScalar2 = (iter->dVal());
    }
    else if (tag == "zeroresistancetol")
    {
      zeroResistanceTol = (iter->dVal());
    }
    else if (tag == "checkforzeroresistance")
    {
      checkForZeroResistance = (iter->bVal());
    }
    else if (tag == "detailed_device_counts" )
    {
      detailedDeviceCounts = (iter->bVal());
    }
    else if (tag == "randomseed" )
    {
      randomSeed = (iter->lVal());
    }
    else if (tag == "trytocompact")
    {
      tryToCompact = static_cast<bool> (iter->iVal());
    }
    else if (tag == "calculateallleadcurrents")
    {
      calculateAllLeadCurrents = (iter->bVal());
    }
    else
    {
      string tmp(
        "  ***** ERROR: " + tag +
        " is not a recognized device package option.\n");
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0,  tmp);
    }
  }

  applyCmdLineOptions ();

  gmin_orig = gmin;
  gmin_init = gmin*gmin_scalar;  // by default, 10 orders of magnitude larger.

#ifdef Xyce_DEBUG_DEVICE
  ostringstream ost;

  ost << *this;
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, ost.str());
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::applyCmdLineOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/24/05
//-----------------------------------------------------------------------------
bool DeviceOptions::applyCmdLineOptions ()
{
#ifdef Xyce_DEBUG_DEVICE
  // set (or override) debug levels based on command line options
  if ( commandLine.argExists( "-sdl" ) )
  {
    sensDebugLevel = atoi( commandLine.getArgumentValue( "-sdl" ).c_str() );
  }

  if ( commandLine.argExists( "-ddl" ) )
  {
    debugLevel = atoi( commandLine.getArgumentValue( "-ddl" ).c_str() );
  }
#endif

  if (commandLine.argExists("-jacobian_test"))
  {
    testJacobianFlag = true;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DeviceOptions::operator=
// Purpose       : assignment operator
// Special Notes :
// Scope         :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/16/05
//-----------------------------------------------------------------------------
DeviceOptions & DeviceOptions::operator=(DeviceOptions const & rhs)
{
  commandLine = rhs.commandLine;

  defad = rhs.defad;
  defas = rhs.defas;
  defl = rhs.defl;
  defw = rhs.defw;
  abstol = rhs.abstol;
  reltol = rhs.reltol;
  chgtol = rhs.chgtol;
  gmin = rhs.gmin;
  gmin_orig = rhs.gmin_orig;
  gmin_init = rhs.gmin_init;
  gmin_scalar = rhs.gmin_scalar;
  gmax = rhs.gmax;
  testJac_relTol = rhs.testJac_relTol;
  testJac_absTol = rhs.testJac_absTol;
  testJac_SqrtEta = rhs.testJac_SqrtEta;
  deviceSens_dp = rhs.deviceSens_dp;
  tnom = rhs.tnom;
  temp = rhs.temp;
  scale_src = rhs.scale_src;
  numericalJacobianFlag = rhs.numericalJacobianFlag;
  testJacobianFlag = rhs.testJacobianFlag;
  testJacStartStep = rhs.testJacStartStep;
  testJacStopStep = rhs.testJacStopStep;
  testJacWarn = rhs.testJacWarn;
  voltageLimiterFlag = rhs.voltageLimiterFlag;
  lambertWFlag = rhs.lambertWFlag;
  icMultiplier = rhs.icMultiplier;
  defaultMaxTimeStep = rhs.defaultMaxTimeStep;
  vdsScaleMin = rhs.vdsScaleMin;
  vgstConst = rhs.vgstConst;
  numGainScaleBlocks = rhs.numGainScaleBlocks;
  staggerGainScale = rhs.staggerGainScale;
  randomizeVgstConst = rhs.randomizeVgstConst;
  length0 = rhs.length0;
  width0 = rhs.width0;
  tox0 = rhs.tox0;
  minRes = rhs.minRes;
  minCap = rhs.minCap;
#ifdef Xyce_DEBUG_DEVICE
  sensDebugLevel = rhs.sensDebugLevel;
  debugLevel = rhs.debugLevel;
  debugMinTimestep = rhs.debugMinTimestep;
  debugMaxTimestep = rhs.debugMaxTimestep;
  debugMinTime = rhs.debugMinTime;
  debugMaxTime = rhs.debugMaxTime;
#endif
  verboseLevel = rhs.verboseLevel;

  blockAnalysisFlag = rhs.blockAnalysisFlag;
  newExcessPhase = rhs.newExcessPhase;
  defaultNewExcessPhase = rhs.defaultNewExcessPhase;
  excessPhaseScalar1 = rhs.excessPhaseScalar1;
  excessPhaseScalar2 = rhs.excessPhaseScalar2 ;
  newMeyerFlag=rhs.newMeyerFlag;
  return *this;
}


#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : DeviceOptions::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/23/05
//-----------------------------------------------------------------------------
ostream & operator<<(ostream & os, const DeviceOptions & devOp)
{
  os << "\n\n-----------------------------------------\n";
  os << "\tDevice Options:\n";
  os << "\t\tdefad                 = " << devOp.defad <<"\n";
  os << "\t\tdefas                 = " << devOp.defas <<"\n";
  os << "\t\tdefl                  = " << devOp.defl <<"\n";
  os << "\t\tdefw                  = " << devOp.defw <<"\n";
  os << "\t\tabstol                = " << devOp.abstol <<"\n";
  os << "\t\treltol                = " << devOp.reltol <<"\n";
  os << "\t\tchgtol                = " << devOp.chgtol <<"\n";
  os << "\t\tgmin                  = " << devOp.gmin <<"\n";
  os << "\t\tgmin_orig             = " << devOp.gmin_orig <<"\n";
  os << "\t\tgmin_init             = " << devOp.gmin_init <<"\n";
  os << "\t\tgmin_scalar           = " << devOp.gmin_scalar <<"\n";
  os << "\t\tgmax                  = " << devOp.gmax <<"\n";
  os << "\t\ttnom                  = " << devOp.tnom <<"\n";
  //os << "\t\tN_UTL_Param temp    = " << devOp.N_UTL_Param temp <<"\n";
  os << "\t\tscale_src             = " << devOp.scale_src <<"\n";
  os << "\t\tnumericalJacobianFlag = " << devOp.numericalJacobianFlag <<"\n";
  os << "\t\ttestJacobianFlag      = " << devOp.testJacobianFlag << "\n";
  os << "\t\ttestJacStartStep      = " << devOp.testJacStartStep << "\n";
  os << "\t\ttestJacStopStep       = " << devOp.testJacStopStep << "\n";
  os << "\t\ttestJacWarn           = " << devOp.testJacWarn     << "\n";
  os << "\t\ttestJacDeviceName     = " << devOp.testJacDeviceName << "\n";

  os << "\t\ttestJac_relTol        = " << devOp.testJac_relTol << "\n";
  os << "\t\ttestJac_absTol        = " << devOp.testJac_absTol << "\n";
  os << "\t\ttestJac_SqrtEta       = " << devOp.testJac_SqrtEta << "\n";
  os << "\t\tdeviceSens_dp         = " << devOp.deviceSens_dp << "\n";

  os << "\t\tvoltageLimiterFlag    = " << devOp.voltageLimiterFlag <<"\n";
  os << "\t\tlambertWFlag          = " << devOp.lambertWFlag <<"\n";
  os << "\t\ticMultiplier          = " << devOp.icMultiplier <<"\n";
  os << "\t\tdefaultMaxTimeStep    = " << devOp.defaultMaxTimeStep <<"\n";
  os << "\t\tvdsScaleMin           = " << devOp.vdsScaleMin <<"\n";
  os << "\t\tvgstConst             = " << devOp.vgstConst <<"\n";
  os << "\t\tnumGainScaleBlocks    = " << devOp.numGainScaleBlocks <<"\n";
  os << "\t\tstaggerGainScale      = " << devOp.staggerGainScale <<"\n";
  os << "\t\trandomizeVgstConst    = " << devOp.randomizeVgstConst <<"\n";
  os << "\t\tlength0               = " << devOp.length0 <<"\n";
  os << "\t\twidth0                = " << devOp.width0 <<"\n";
  os << "\t\ttox0                  = " << devOp.tox0 <<"\n";
  os << "\t\tminres                = " << devOp.minRes <<"\n";
  os << "\t\tmincap                = " << devOp.minCap <<"\n";
#ifdef Xyce_DEBUG_DEVICE
  os << "\t\tsensDebugLevel        = " << devOp.sensDebugLevel <<"\n";
  os << "\t\tdebugLevel            = " << devOp.debugLevel <<"\n";
  os << "\t\tdebugMinTimestep      = " << devOp.debugMinTimestep <<"\n";
  os << "\t\tdebugMaxTimestep      = " << devOp.debugMaxTimestep <<"\n";
  os << "\t\tdebugMinTime          = " << devOp.debugMinTime <<"\n";
  os << "\t\tdebugMaxTime          = " << devOp.debugMaxTime <<"\n";
#endif
  os << "\t\tverboseLevel          = " << devOp.verboseLevel <<"\n";
  os << "\t\tblockAnalysisFlag     = " << devOp.blockAnalysisFlag << "\n";
  os << "\t\tnewExcessPhase        = " << devOp.newExcessPhase << "\n";
  os << "\t\tdefaultNewExcessPhase = " << devOp.defaultNewExcessPhase << "\n";
  os << "\t\texcessPhaseScalar1    = " << devOp.excessPhaseScalar1 << "\n";
  os << "\t\texcessPhaseScalar2    = " << devOp.excessPhaseScalar2 << "\n";
  os << "\t\tnewMeyerFlag    = " << devOp.newMeyerFlag << "\n";
  os << "-----------------------------------------\n";
  os << endl;

  return os;
}
#endif

} // namespace Device
} // namespace Xyce
