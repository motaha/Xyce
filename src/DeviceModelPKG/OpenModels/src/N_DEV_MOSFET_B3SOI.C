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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_MOSFET_B3SOI.C,v $
//
// Purpose        : This file implements the B3SOI MOSFET model.  It
//                  is intended to be compatible with the Berkeley SPICE
//                  (3f5) version, B3SOI version 3.2.
//
// Special Notes  :
//
//
// Creator        : Dave Shirley
//
// Creation Date  : 05/20/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.287.2.4 $
//
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <sstream>

// ----------   Xyce Includes   ----------

#include <N_DEV_Const.h>
#include <N_DEV_MOSFET_B3SOI.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<MOSFET_B3SOI::Instance>::ParametricData()
{
    // Set up configuration constants:
    setNumNodes(4);
    setNumOptionalNodes(3);
    setNumFillNodes(0);
    setModelRequired(1);
    addModelType("NMOS");
    addModelType("PMOS");

    // Set up double precision variables:
    addPar ("L", 5.0e-6, false, NO_DEP,
      &MOSFET_B3SOI::Instance::l,
      NULL, U_METER,  CAT_GEOMETRY, "Channel length");

    addPar ("W", 5.0e-6, false, NO_DEP,
      &MOSFET_B3SOI::Instance::w,
      NULL, U_METER,  CAT_GEOMETRY, "Channel width");

    addPar ("AD", 0.0, false,   NO_DEP,
      &MOSFET_B3SOI::Instance::drainArea,
      NULL, U_METER2,  CAT_GEOMETRY, "Drain diffusion area");

    addPar ("AS", 0.0, false,   NO_DEP,
      &MOSFET_B3SOI::Instance::sourceArea,
      NULL, U_METER2,  CAT_GEOMETRY, "Source diffusion area");

    addPar ("NRD", 1.0, false,  NO_DEP,
      &MOSFET_B3SOI::Instance::drainSquares,
      NULL, U_SQUARES, CAT_GEOMETRY, "Multiplier for RSH to yield parasitic resistance of drain");

    addPar ("NRS", 1.0, false,  NO_DEP,
      &MOSFET_B3SOI::Instance::sourceSquares,
      NULL, U_SQUARES, CAT_GEOMETRY, "Multiplier for RSH to yield parasitic resistance of source");

    addPar ("PD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Instance::drainPerimeter,
      NULL, U_METER, CAT_GEOMETRY, "Drain diffusion perimeter");

    addPar ("PS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Instance::sourcePerimeter,
      NULL, U_METER, CAT_GEOMETRY, "Source diffusion perimeter");

    addPar ("IC1", 0.0, false,NO_DEP,
      &MOSFET_B3SOI::Instance::icVDS,
      &MOSFET_B3SOI::Instance::icVDSGiven,
       U_VOLT, CAT_VOLT, "Initial condition on Vds");

    addPar ("IC2", 0.0, false,NO_DEP,
      &MOSFET_B3SOI::Instance::icVGS,
      &MOSFET_B3SOI::Instance::icVGSGiven,
      U_VOLT, CAT_VOLT, "Initial condition on Vgs");

    addPar ("IC3", 0.0, false,NO_DEP,
      &MOSFET_B3SOI::Instance::icVBS,
      &MOSFET_B3SOI::Instance::icVBSGiven,
       U_VOLT, CAT_VOLT, "Initial condition on Vbs");

    addPar ("IC4", 0.0, false,NO_DEP,
      &MOSFET_B3SOI::Instance::icVES,
      &MOSFET_B3SOI::Instance::icVESGiven,
       U_VOLT, CAT_VOLT, "Initial condition on Ves");

    addPar ("IC5", 0.0, false,NO_DEP,
      &MOSFET_B3SOI::Instance::icVPS,
      &MOSFET_B3SOI::Instance::icVPSGiven,
       U_VOLT, CAT_VOLT, "Initial condition on Vps");

    addPar ("TEMP", 27.0, false,TIME_DEP,
      &MOSFET_B3SOI::Instance::temp,
      NULL,   STANDARD, CAT_NONE, "");

    addPar ("RTH0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Instance::rth0,
      NULL, U_OHM, CAT_TEMP, "normalized thermal resistance");

    addPar ("CTH0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Instance::cth0,
      NULL, U_FARAD, CAT_TEMP, "Thermal capacitance");

    addPar ("NRB", 1.0, false,  NO_DEP,
      &MOSFET_B3SOI::Instance::bodySquares,
      NULL, U_NONE, CAT_GEOMETRY, "Number of squares in body");

    addPar ("FRBODY", 1.0, false,NO_DEP,
      &MOSFET_B3SOI::Instance::frbody,
      NULL, U_NONE, CAT_GEOMETRY, "Layout dependent body-resistance coefficient");

    addPar ("NBC", 0.0, false,   NO_DEP,
      &MOSFET_B3SOI::Instance::nbc,
      NULL, U_NONE, CAT_GEOMETRY, "Number of body contact isolation edge");

    addPar ("NSEG", 1.0, false,  NO_DEP,
      &MOSFET_B3SOI::Instance::nseg,
      NULL, U_NONE, CAT_GEOMETRY, "Number segments for width partitioning");

    addPar ("PDBCP", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Instance::pdbcp,
      NULL, U_METER, CAT_GEOMETRY, "Perimeter length for bc parasitics at drain side");

    addPar ("PSBCP", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Instance::psbcp,
      NULL, U_METER, CAT_GEOMETRY, "Perimeter length for bc parasitics at source side");

    addPar ("AGBCP", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Instance::agbcp,
      NULL, U_METER2, CAT_GEOMETRY, "Gate to body overlap area for bc parasitics");

    addPar ("AEBCP", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Instance::aebcp,
      NULL, U_METER2, CAT_GEOMETRY, "Substrate to body overlap area for bc prasitics");

    addPar ("VBSUSR", 0.0, false,NO_DEP,
      &MOSFET_B3SOI::Instance::vbsusr,
      NULL, U_VOLT, CAT_DC, "Vbs specified by user");

    addPar ("M", 1.0, false,   NO_DEP,
      &MOSFET_B3SOI::Instance::numberParallel,
      NULL, U_NONE, CAT_CONTROL,  "Multiplier for M devices connected in parallel");

    // Set up non-double precision variables:
    addPar ("OFF",false,false, NO_DEP,
            &MOSFET_B3SOI::Instance::OFF,
            NULL, U_LOGIC, CAT_VOLT,
            "Initial condition of no voltage drops accross device");

    addPar ("BJTOFF", 0, false, NO_DEP,
            &MOSFET_B3SOI::Instance::bjtoff, NULL,
           U_LOGIC, CAT_NONE, "BJT on/off flag");
    addPar ("DEBUG", 0, false, NO_DEP,
            &MOSFET_B3SOI::Instance::debugMod, NULL,
           U_LOGIC, CAT_NONE, "BJT on/off flag");
    addPar ("SOIMOD", 0, false, NO_DEP,
            &MOSFET_B3SOI::Instance::soiMod, NULL,
           U_NONE, CAT_CONTROL, "SIO model selector, SOIMOD=0: BSIMPD, SOIMOD=1: undefined model for PD and FE, SOIMOD=2: ideal FD");
    addPar ("TNODEOUT", 0, false, NO_DEP,
            &MOSFET_B3SOI::Instance::tnodeout, NULL,
           U_LOGIC, CAT_NONE, "Flag indicating external temp node");
    addPar ("RGATEMOD", 0, false, NO_DEP,
            &MOSFET_B3SOI::Instance::rgateMod, NULL,
           U_NONE, CAT_RF, "Gate resistance model selector");
    addPar ("VLDEBUG", false, false, NO_DEP,
            &MOSFET_B3SOI::Instance::vlDebug, NULL,
           U_LOGIC, CAT_NONE, "");
    // This tells the parser that IC1 - IC5 are to be input as a vector of "IC"
    makeVector ("IC", 5);
}

template<>
ParametricData<MOSFET_B3SOI::Model>::ParametricData()
{
    // Set up double precision variables:
    addPar ("TOX", 100.0e-10, false, NO_DEP,
      &MOSFET_B3SOI::Model::tox,
      NULL, U_METER, CAT_GEOMETRY, "Gate oxide thickness");

    addPar ("TOXM", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::toxm,
      NULL,  U_METER, CAT_PROCESS, "Gate oxide thickness used in extraction");

    addPar ("DTOXCV", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dtoxcv,
      NULL, U_METER, CAT_NONE, "Delta oxide thickness in meters in CapMod3");

    addPar ("CDSC", 2.4e-4, false, NO_DEP,
      &MOSFET_B3SOI::Model::cdsc,
      NULL,  U_FARADMM2, CAT_DC, "Drain/source to channel coupling capacitance");

    addPar ("CDSCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cdscb,
      NULL,  U_FVM1MM2, CAT_DC, "Body-bias sensitivity of CDSC");

    addPar ("CDSCD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cdscd,
      NULL,  U_FVM1MM2, CAT_DC, "Drain-bias sensitivity of CDSC");

    addPar ("CIT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cit,
      NULL,  U_FARADMM2, CAT_DC, "Interface trap capacitance");

    addPar ("NFACTOR", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::nfactor,
      NULL,  U_NONE, CAT_DC, "Subthreshold swing factor");

    addPar ("VSAT", 8.0e4, false, NO_DEP,
      &MOSFET_B3SOI::Model::vsat,
      NULL,  U_MSM1, CAT_DC, "Saturation velocity at temp = TNOM");

    addPar ("AT", 3.3e4, false, NO_DEP,
      &MOSFET_B3SOI::Model::at,
      NULL,  U_MSM1, CAT_TEMP, "Temperature coefficient for saturation velocity");

    addPar ("A0", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::a0,
      NULL,  U_NONE, CAT_DC, "Bulk charge effect coefficient for channel length");

    addPar ("AGS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ags,
      NULL,  U_VOLTM1, CAT_DC, "Gate-bias coefficient of abulk");

    addPar ("A1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::a1,
      NULL,  U_VOLTM1, CAT_DC, "First non-saturation effect parameter");

    addPar ("A2", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::a2,
      NULL,  U_NONE, CAT_DC, "Second non-saturation factor");

    addPar ("KETA", -0.6, false, NO_DEP,
      &MOSFET_B3SOI::Model::keta,
      NULL,  U_VOLTM1, CAT_DC, "Body-bias coefficient of bulk charge effect");

    addPar ("NSUB", 6.0e16, false, NO_DEP,
      &MOSFET_B3SOI::Model::nsub,
      NULL, U_CMM3, CAT_DOPING, "Substrate doping density");

    addPar ("NCH", 1.7e17, false, NO_DEP,
      &MOSFET_B3SOI::Model::npeak,
      &MOSFET_B3SOI::Model::npeakGiven,    
       U_CMM3, CAT_PROCESS, "Channel doping concentration");

    addPar ("NGATE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ngate,
      NULL,  U_CMM3, CAT_DC, "Poly gate doping concentration");

    addPar ("GAMMA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::gamma1,
      &MOSFET_B3SOI::Model::gamma1Given,
       U_VOLTH, CAT_PROCESS, "Body effect coefficient near the surface");

    addPar ("GAMMA2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::gamma2,
      &MOSFET_B3SOI::Model::gamma1Given,
       U_VOLTH, CAT_PROCESS, "Body effect coefficient in the bulk");

    addPar ("VBX", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vbx,
      &MOSFET_B3SOI::Model::vbxGiven,
       U_VOLT, CAT_PROCESS, "Vbs at which the depetion region = XT");

    addPar ("VBM", -3.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vbm,
      &MOSFET_B3SOI::Model::vbmGiven,
       U_VOLT, CAT_DC, "Maximum applied body-bias in threshold voltage calculation");

    addPar ("XT", 1.55e-7, false, NO_DEP,
      &MOSFET_B3SOI::Model::xt,
      &MOSFET_B3SOI::Model::xtGiven,
       U_METER, CAT_PROCESS, "Doping depth");

    addPar ("K1", 0.53, false, NO_DEP,
      &MOSFET_B3SOI::Model::k1,
      &MOSFET_B3SOI::Model::k1Given,
       U_VOLTH, CAT_DC, "First-order body effect coefficient");

    addPar ("KT1", -0.11, false, NO_DEP,
      &MOSFET_B3SOI::Model::kt1,
      NULL,  U_VOLT, CAT_TEMP, "Themperature coefficient for threshold voltage");

    addPar ("KT1L", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::kt1l,
      NULL,  U_VM, CAT_TEMP, "Channel length dependence of the temerature coefficient for the threshold voltage");

    addPar ("KT2", 0.022, false, NO_DEP,
      &MOSFET_B3SOI::Model::kt2,
      NULL,  U_NONE, CAT_TEMP, "Body-bias coefficient fo the threshold voltage temperature effect");

    addPar ("K2", -0.0186, false, NO_DEP,
      &MOSFET_B3SOI::Model::k2,
      &MOSFET_B3SOI::Model::k2Given,
       U_NONE, CAT_DC, "second-order body effect coefficient");

    addPar ("K3", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::k3,
      NULL,  U_NONE, CAT_DC, "Narrow width coefficient");

    addPar ("K3B", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::k3b,
      NULL,  U_VOLTM1, CAT_DC, "Body effect coefficient of K3");

    addPar ("W0", 2.5e-6, false, NO_DEP,
      &MOSFET_B3SOI::Model::w0,
      NULL,  U_METER, CAT_DC, "Narrow-width paameter");

    addPar ("NLX", 1.74e-7, false, NO_DEP,
      &MOSFET_B3SOI::Model::nlx,
      NULL,  U_METER, CAT_DC, "Lateral non-uniform doping parameter");

    addPar ("DVT0", 2.2, false, NO_DEP,
      &MOSFET_B3SOI::Model::dvt0,
      NULL,  U_NONE, CAT_DC, "First coefficient of short-channel effect effect on threshold voltage");

    addPar ("DVT1", 0.53, false, NO_DEP,
      &MOSFET_B3SOI::Model::dvt1,
      NULL,  	 U_NONE, CAT_DC, "Second coefficient of short-channel effect effect on threshold voltage");

    addPar ("DVT2", -0.032, false, NO_DEP,
      &MOSFET_B3SOI::Model::dvt2,
      NULL,  	 U_VOLTM1, CAT_DC, "Body-bias coefficient of short-channel effect effect on threshold voltage");

    addPar ("DVT0W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dvt0w,
      NULL,  	 U_METERM1, CAT_DC, "First coefficient of narrow-width effect effect on threshold voltage for small channel length");

    addPar ("DVT1W", 5.3e6, false, NO_DEP,
      &MOSFET_B3SOI::Model::dvt1w,
      NULL,  	 U_METERM1, CAT_DC, "Second coefficient of narrow-width effect effect on threshold voltage for small channel length");

    addPar ("DVT2W", -0.032, false, NO_DEP,
      &MOSFET_B3SOI::Model::dvt2w,
      NULL,  	 U_VOLTM1, CAT_DC, "Body-bias coefficient of narrow-width effect effect on threshold voltage for small channel length");

    addPar ("DROUT", 0.56, false, NO_DEP,
      &MOSFET_B3SOI::Model::drout,
      NULL,  U_NONE, CAT_DC, "L-depedance Coefficient of the DIBL correction parameter in Rout");

    addPar ("DSUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dsub,
      NULL,  U_NONE, CAT_DC, "DIBL coefficient exponent in subthreshhold region");

    addPar ("UA", 2.25e-9, false, NO_DEP,
      &MOSFET_B3SOI::Model::ua,
      NULL,  U_MVM1, CAT_DC, "First-order mobility degradation coefficient");

    addPar ("UA1", 4.31e-9, false, NO_DEP,
      &MOSFET_B3SOI::Model::ua1,
      NULL,  U_MVM1, CAT_TEMP, "Temperature coefficient for UA");

    addPar ("UB", 5.87e-19, false, NO_DEP,
      &MOSFET_B3SOI::Model::ub,
      NULL,  U_M2VM2, CAT_DC, "First-order mobility degradation coefficient");

    addPar ("UB1", -7.61e-18, false, NO_DEP,
      &MOSFET_B3SOI::Model::ub1,
      NULL,  U_M2VM2, CAT_TEMP, "Temperature coefficient for UB");

    addPar ("UC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::uc,
      NULL,  U_MVM2, CAT_DC, "Body effect of mobility degridation coefficient");

    addPar ("UC1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::uc1,
      NULL,  U_MVM2DEGCM1, CAT_TEMP, "Temperature coefficient for UC");

    addPar ("U0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::u0,
      NULL,  
       U_CMM2VM1SM1, CAT_PROCESS, "Surface mobility");

    addPar ("UTE", -1.5, false, NO_DEP,
      &MOSFET_B3SOI::Model::ute,
      NULL,  U_NONE, CAT_TEMP, "Mobility temerature exponent");

    addPar ("VOFF", -0.08, false, NO_DEP,
      &MOSFET_B3SOI::Model::voff,
      NULL,  U_VOLT, CAT_DC, "Offset voltage in the subthreshold region at large W and L");

    addPar ("TNOM", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::tnom,
      NULL, STANDARD, CAT_NONE, "");

    addPar ("CGSO", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cgso,
      NULL,  U_FARADMM1, CAT_CAP, "Non-LLD region source-gate overlap capacitance per unit channel length");

    addPar ("CGDO", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cgdo,
      NULL,  U_FARADMM1, CAT_CAP, "Non-LLD region drain-gate overlap capacitance per unit channel length");

    addPar ("XPART", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xpart,
      NULL,  U_NONE, CAT_CAP, "Charge partitioning rate flag");

    addPar ("DELTA", 0.01, false, NO_DEP,
      &MOSFET_B3SOI::Model::delta,
      NULL,  U_VOLT, CAT_DC, "Effective Vds parameter");

    addPar ("RSH", 0.0, false, MIN_RES,
      &MOSFET_B3SOI::Model::sheetResistance,
      NULL, U_OHM, CAT_RES, "Drain, source diffusion sheet resistance");

    addPar ("RDSW", 100.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::rdsw,
      NULL,  U_OHMMICRON, CAT_DC, "Parasitic resistance per unit width");

    addPar ("PRWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::prwg,
      NULL,  U_VOLTM1, CAT_DC, "Gate-bias effect coefficient of RDSW");

    addPar ("PRWB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::prwb,
      NULL,  U_VOLTMH, CAT_DC, "Body effect coefficient of RDSW");

    addPar ("PRT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::prt,
      NULL,  U_OHMMICRON, CAT_TEMP, "Temerature coefficient for RDSW");

    addPar ("ETA0", 0.08, false, NO_DEP,
      &MOSFET_B3SOI::Model::eta0,
      NULL,  U_NONE, CAT_DC, "DIBL coefficient in subthreshold region");

    addPar ("ETAB", -0.07, false, NO_DEP,
      &MOSFET_B3SOI::Model::etab,
      NULL,  U_VOLTM1, CAT_DC, "Body-bias coefficient for the subthreshold DIBL effect");

    addPar ("PCLM", 1.3, false, NO_DEP,
      &MOSFET_B3SOI::Model::pclm,
      NULL,  U_NONE, CAT_DC, "Channel length modulation parameter");

    addPar ("PDIBLC1", 0.39, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdibl1,
      NULL,  U_NONE, CAT_DC, "First output resistance DIBL effect correction parameter");

    addPar ("PDIBLC2", 0.0086, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdibl2,
      NULL,  U_NONE, CAT_DC, "Second output resistance DIBL effect correction parameter");

    addPar ("PDIBLCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdiblb,
      NULL,  U_VOLTM1, CAT_DC, "Body effect coefficient of DIBL correction parameter");

    addPar ("PVAG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvag,
      NULL,  U_NONE, CAT_DC, "Gate dependence of early voltage");

    addPar ("SHMOD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::shMod,
      NULL, U_NONE, CAT_CONTROL, "Flag for self-heating, 0-no self-heating, 1-self-heating");

    addPar ("TBOX", 3.0e-7, false, NO_DEP,
      &MOSFET_B3SOI::Model::tbox,
      NULL, U_METER, CAT_PROCESS, "Buried oxide thickness");

    addPar ("TSI", 1.0e-7, false, NO_DEP,
      &MOSFET_B3SOI::Model::tsi,
      NULL, U_METER, CAT_PROCESS, "Silicon film thickness");

    addPar ("XJ", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xj,
      NULL,  U_METER, CAT_GEOMETRY, "Junction depth");

    addPar ("RTH0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::rth0,
      NULL, U_OHMMM1, CAT_TEMP, "Thermal resistance per unit width");

    addPar ("CTH0", 1.0e-5, false, NO_DEP,
      &MOSFET_B3SOI::Model::cth0,
      NULL, U_FARADMM1, CAT_TEMP, "Thermal capacitance per unit width");

    addPar ("NGIDL", 1.2, false, NO_DEP,
      &MOSFET_B3SOI::Model::ngidl,
      NULL, U_VOLT, CAT_DC, "GIDL Vds enhancement coefficient");

    addPar ("AGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::agidl,
      NULL, U_OHMM1, CAT_DC, "GIDL constant");

    addPar ("BGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::bgidl,
      NULL, U_VMM1, CAT_DC, "GIDL exponential coefficient");

    addPar ("NDIODE", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ndiode,
      NULL, U_NONE, CAT_DC, "Diode non-ideality factor");

    addPar ("XBJT", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xbjt,
      NULL, U_NONE, CAT_TEMP, "Power dependence of JBJT on temperature");

    addPar ("XDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xdif,
      NULL, U_NONE, CAT_TEMP, "Power dependence of JDIF on temperature");

    addPar ("XREC", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xrec,
      NULL, U_NONE, CAT_TEMP, "Power dependence of JREC on temperature");

    addPar ("XTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xtun,
      NULL, U_NONE, CAT_TEMP, "Power dependence of JTUN on temperature");

    addPar ("PBSWG", 0.7, false, NO_DEP,
      &MOSFET_B3SOI::Model::GatesidewallJctPotential,
      NULL,  U_VOLT, CAT_CAP, "Source/drain gate sidewall junction built-in potential");

    addPar ("MJSWG", 0.5, false, NO_DEP,
      &MOSFET_B3SOI::Model::bodyJctGateSideGradingCoeff,
      NULL,  U_NONE, CAT_CAP, "Source/grain gate sidewall junction capacitance grading coeficient");

    addPar ("CJSWG", 1.0e-10, false, NO_DEP,
      &MOSFET_B3SOI::Model::unitLengthGateSidewallJctCap,
      NULL,  U_FARADMM1, CAT_CAP, "Source/grain gate sidewall junction capacitance per unit width");

    addPar ("LINT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lint,
      NULL,  U_METER, CAT_DC, "Length of offset fiting parameter from I-V without bias");

    addPar ("LL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Ll,
      NULL,  U_MEXPLL, CAT_GEOMETRY, "Coefficient of length dependence for length offset");

    addPar ("LLC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Llc,
      NULL,  U_MEXPLL, CAT_GEOMETRY, "Coefficient of length dependence for CV channel length offset");

    addPar ("LLN", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lln,
      NULL,  U_NONE, CAT_GEOMETRY, "Power of length dependence for length offset");

    addPar ("LW", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lw,
      NULL,  U_MEXPLW, CAT_GEOMETRY, "Coefficient of width dependence for length offset");

    addPar ("LWC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lwc,
      NULL,  U_MEXPLW, CAT_GEOMETRY, "Coefficient of width dependence for channel length offset");

    addPar ("LWN", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lwn,
      NULL,  U_NONE, CAT_GEOMETRY, "Power of width dependence for length offset");

    addPar ("LWL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lwl,
      NULL,  U_MEXPLLLW, CAT_GEOMETRY, "Coefficient of length and width cross term for length offset");

    addPar ("LWLC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lwlc,
      NULL,  U_MEXPLLLW, CAT_GEOMETRY, "Coefficient of length and width dependence for CV channel length offset");

    addPar ("WR", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wr,
      NULL,  U_NONE, CAT_DC, "Width offset from Weff for Rds Calculation");

    addPar ("WINT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wint,
      NULL,  U_METER, CAT_DC, "Width-offset fitting parameter from I-V without bias");

    addPar ("DWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dwg,
      NULL,  U_MVMH, CAT_DC, "Coefficient of gate depedence of Weff");

    addPar ("DWB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dwb,
      NULL,  U_MVMH, CAT_DC, "Coefficient of substrate body bias dependence of Weff");

    addPar ("WL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wl,
      NULL,  U_MEXPWL, CAT_GEOMETRY, "Coefficient of length dependence for width offset");

    addPar ("WLC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wlc,
      NULL,  U_MEXPWL, CAT_GEOMETRY, "Coefficient of length dependence for CV channel width offset");

    addPar ("WLN", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wln,
      NULL,  U_NONE, CAT_GEOMETRY, "Power of length dependece of width offset");

    addPar ("WW", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Ww,
      NULL,  U_MEXPWW, CAT_GEOMETRY, "Coefficient of width dependence for width offset");

    addPar ("WWC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wwc,
      NULL,  U_MEXPWW, CAT_GEOMETRY, "Coefficient of width dependence for CV channel width offset");

    addPar ("WWN", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wwn,
      NULL,  U_NONE, CAT_GEOMETRY, "Power of width dependence of width offset");

    addPar ("WWL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wwl,
      NULL,  U_MEXPWLWW, CAT_GEOMETRY, "Coefficient of length and width cross term for width offset");

    addPar ("WWLC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wwlc,
      NULL,  U_MEXPWLWW, CAT_GEOMETRY, "Coefficient of length and width dependence for CV channel width offset");

    addPar ("B0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::b0,
      NULL,  U_METER, CAT_DC, "Bulk charge effect coefficient for channel width");

    addPar ("B1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::b1,
      NULL,  U_METER, CAT_DC, "Bulk charge effect offset");

    addPar ("CGSL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cgsl,
      NULL,  U_FARADMM1, CAT_CAP, "Light-doped source-gate region overlap capacitance");

    addPar ("CGDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cgdl,
      NULL,  U_FARADMM1, CAT_CAP, "Light-doped drain-gate region overlap capacitance");

    addPar ("CKAPPA", 0.6, false, NO_DEP,
      &MOSFET_B3SOI::Model::ckappa,
      NULL,  U_FARADMM1, CAT_CAP, "Coefficient for lightly doped region overlap capacitance fireing field capacitance");

    addPar ("CF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cf,
      NULL,  U_FARADMM1, CAT_CAP, "Firing field capacitance");

    addPar ("CLC", 0.1e-7, false,NO_DEP,
      &MOSFET_B3SOI::Model::clc,
      NULL,  U_METER, CAT_CAP, "Constant term for short-channel model");

    addPar ("CLE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cle,
      NULL,  U_NONE, CAT_CAP, "Exponetial term for the short-channel model");

    addPar ("DWC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dwc,
      NULL,  U_METER, CAT_CAP, "Width offset fitting parameter from C-V");

    addPar ("DLC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dlc,
      NULL,  U_METER, CAT_CAP, "Length offset fitting parameter from C-V");

    addPar ("ALPHA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::alpha0,
      NULL,  U_MVM1, CAT_DC, "First parameter of impact-ionization current");

    addPar ("NOIA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::oxideTrapDensityA,
      NULL,  U_NONE, CAT_FLICKER, "Noise parameter a");

    addPar ("NOIB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::oxideTrapDensityB,
      NULL,  U_NONE, CAT_FLICKER, "Noise parameter b");

    addPar ("NOIC", 8.75e9, false, NO_DEP,
      &MOSFET_B3SOI::Model::oxideTrapDensityC,
      NULL,  U_NONE, CAT_FLICKER, "Noise parameter c");

    addPar ("FNOIMOD", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::fnoiMod,
      NULL, U_NONE, CAT_NONE, "Flicker noise model selector");

    addPar ("TNOIMOD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::tnoiMod,
      NULL, U_NONE, CAT_NONE, "Thermal noise model selector");

    addPar ("TNOIA", 1.5, false, NO_DEP,
      &MOSFET_B3SOI::Model::tnoia,
      NULL, U_NONE, CAT_NONE, "Thermal noise parameter");

    addPar ("TNOIB", 3.5, false, NO_DEP,
      &MOSFET_B3SOI::Model::tnoib,
      NULL, U_NONE, CAT_NONE, "Thermal noise parameter");

    addPar ("RNOIA", 0.577, false, NO_DEP,
      &MOSFET_B3SOI::Model::rnoia,
      NULL, U_NONE, CAT_NONE, "Thermal noise coefficient");

    addPar ("RNOIB", 0.37, false, NO_DEP,
      &MOSFET_B3SOI::Model::rnoib,
      NULL, U_NONE, CAT_NONE, "Thermal noise coefficient");

    addPar ("NTNOI", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ntnoi,
      NULL, U_NONE, CAT_NONE, "Thermal noise parameter");

    addPar ("EM", 4.1e7, false, NO_DEP,
      &MOSFET_B3SOI::Model::em,
      NULL,  U_VMM1, CAT_FLICKER, "Saturation field");

    addPar ("EF", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ef,
      NULL,  U_NONE, CAT_FLICKER, "Flicker exponent");

    addPar ("AF", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::af,
      NULL, U_NONE, CAT_FLICKER, "Flicker noise exponent");

    addPar ("KF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::kf,
      NULL, U_NONE, CAT_FLICKER, "Flicker noise coefficient");

    addPar ("NOIF", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::noif,
      NULL, U_NONE, CAT_NONE, "Floating body excess noise ideality factor");

    addPar ("K1W1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::k1w1,
      NULL, U_METER, CAT_DC, "First body effect width depenent parameter");

    addPar ("K1W2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::k1w2,
      NULL, U_METER, CAT_DC, "Second body effect width depenent parameter");

    addPar ("KETAS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ketas,
      NULL, U_VOLT, CAT_DC, "Surface potential adjustment for bulk charge effect");

    addPar ("DWBC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dwbc,
      NULL, U_METER, CAT_DC, "Width offset for body contact isolation edge");

    addPar ("BETA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::beta0,
      NULL,  U_VOLT, CAT_DC, "Second parameter of impact-ionization current");

    addPar ("BETA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::beta1,
      NULL, U_NONE, CAT_DC, "Second Vds dependent parameter of impact ionizatin current");

    addPar ("BETA2", 0.1, false, NO_DEP,
      &MOSFET_B3SOI::Model::beta2,
      NULL, U_VOLT, CAT_DC, "Third Vds dependent parameter of impact ionizatin current");

    addPar ("VDSATII0", 0.9, false, NO_DEP,
      &MOSFET_B3SOI::Model::vdsatii0,
      NULL, U_VOLT, CAT_DC, "Normal drain saturatio voltage at threshold for impact ionization current");

    addPar ("TII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::tii,
      NULL, U_NONE, CAT_DC, "Temperature dependent parameter for impact ionization current");

    addPar ("LII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lii,
      NULL, U_NONE, CAT_DC, "Channel length dependent parameter at threshold for impact ionization current");

    addPar ("SII0", 0.5, false, NO_DEP,
      &MOSFET_B3SOI::Model::sii0,
      NULL, U_VOLTM1, CAT_DC, "First Vgs dependent parameter of impact ionizatin current");

    addPar ("SII1", 0.1, false, NO_DEP,
      &MOSFET_B3SOI::Model::sii1,
      NULL, U_VOLTM1, CAT_DC, "Second Vgs dependent parameter of impact ionizatin current");

    addPar ("SII2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::sii2,
      NULL, U_NONE, CAT_DC, "Third Vgs dependent parameter of impact ionizatin current");

    addPar ("SIID", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::siid,
      NULL, U_VOLTM1, CAT_DC, "Vds dependent parameter of drain saturation voltage for impact ionizatin current");

    addPar ("FBJTII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::fbjtii,
      NULL, U_NONE, CAT_DC, "Fraction of bipolar current affecting the impact ionization");

    addPar ("ESATII", 1.0e7, false, NO_DEP,
      &MOSFET_B3SOI::Model::esatii,
      NULL, U_VMM1, CAT_DC, "Saturation channel electric field for impact ionization current");

    addPar ("NTUN", 10.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ntun,
      NULL, U_NONE, CAT_DC, "Reverse tunneling non-ideality factor");

    addPar ("NRECF0", 2.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::nrecf0,
      NULL, U_NONE, CAT_DC, "Recombination non-ideality factor at foward bias");

    addPar ("NRECR0", 10.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::nrecr0,
      NULL, U_NONE, CAT_DC, "Recombination non-ideality factor at reverse bias");

    addPar ("ISBJT", 1.0e-6, false, NO_DEP,
      &MOSFET_B3SOI::Model::isbjt,
      NULL, U_AMPMM2, CAT_DC, "BJT injection saturation current");

    addPar ("ISDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::isdif,
      NULL, U_AMPMM2, CAT_DC, "BOdy to source/drain injection saturation current");

    addPar ("ISREC", 1.0e-5, false, NO_DEP,
      &MOSFET_B3SOI::Model::isrec,
      NULL, U_AMPMM2, CAT_DC, "Recombinatin in depletion saturation current");

    addPar ("ISTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::istun,
      NULL, U_AMPMM2, CAT_DC, "Reverse tunneling saturation current");

    addPar ("LN", 2.0e-6, false, NO_DEP,
      &MOSFET_B3SOI::Model::ln,
      NULL, U_METER, CAT_DC, "Electron/hole diffusion length");

    addPar ("VREC0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vrec0,
      NULL, U_VOLT, CAT_DC, "Voltage dependent parameter for recombination current");

    addPar ("VTUN0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vtun0,
      NULL, U_VOLT, CAT_DC, "Voltage dependent parameter for tunneling current");

    addPar ("NBJT", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::nbjt,
      NULL, U_NONE, CAT_DC, "Power coefficient of channel length");

    addPar ("LBJT0", 2.0e-7, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbjt0,
      NULL, U_METER, CAT_DC, "Reference channel length for bipolar current");

    addPar ("LDIF0", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldif0,
      NULL, U_NONE, CAT_CAP, "Channel length dependency coefficient of diffusion capacitance");

    addPar ("VABJT", 10.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vabjt,
      NULL, U_VOLT, CAT_DC, "Early voltage for bipolar current");

    addPar ("AELY", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::aely,
      NULL, U_VMM1, CAT_DC, "Channel length dependency of early voltage for bipolar current");

    addPar ("AHLI", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ahli,
      NULL, U_NONE, CAT_DC, "High level injection parameter for bipolar current");

    addPar ("RBODY", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::rbody,
      NULL, U_OSQM1, CAT_DC, "Intrinsic body contact sheet resistance");

    addPar ("RBSH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::rbsh,
      NULL, U_OSQM1, CAT_DC, "Intrinsic body contact sheet resistance");

    addPar ("CGEO", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cgeo,
      NULL, U_FARADMM1, CAT_CAP, "Gate substrate overlap capacitance per unit channel length");

    addPar ("TT", 1.0e-12, false, NO_DEP,
      &MOSFET_B3SOI::Model::tt,
      NULL, U_SECOND, CAT_CAP, "Diffusion capacitance transit time coefficient");

    addPar ("NDIF", -1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ndif,
      NULL, U_NONE, CAT_CAP, "Power coefficient of channel length dependency for diffusion capacitance");

    addPar ("VSDFB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vsdfb,
      &MOSFET_B3SOI::Model::vsdfbGiven,
    U_VOLT, CAT_CAP, "Sorce/Drain bottom diffusion capacitance flatband voltage");

    addPar ("VSDTH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vsdth,
      &MOSFET_B3SOI::Model::vsdthGiven,
    U_VOLT, CAT_CAP, "Sorce/Drain bottom diffusion capacitance threshold voltage");

    addPar ("CSDMIN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::csdmin,
      &MOSFET_B3SOI::Model::csdminGiven,
    U_VOLT, CAT_CAP, "Sorce/Drain bottom diffusion minimum capacitance");

    addPar ("ASD", 0.3, false, NO_DEP,
      &MOSFET_B3SOI::Model::asd,
      NULL, U_NONE, CAT_CAP, "Sorce/Drain bottom diffusion smoothing parameter");

    addPar ("CSDESW", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::csdesw,
      NULL, U_FARADMM1, CAT_CAP, "Sorce/Drain sidewall fringing capacitance per unit length");

    addPar ("NTRECF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ntrecf,
      NULL, U_NONE, CAT_TEMP, "Temperature coefficient for NRECF");

    addPar ("NTRECR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ntrecr,
      NULL, U_NONE, CAT_TEMP, "Temperature coefficient for NRECR");

    addPar ("DLCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dlcb,
      NULL, U_METER, CAT_CAP, "Length offset fitting parameter for body charge");

    addPar ("FBODY", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::fbody,
      NULL, U_NONE, CAT_CAP, "Scaling factor for body charge");

    addPar ("TCJSWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::tcjswg,
      NULL,  U_KM1, CAT_TEMP, "Temperature coefficient of Cjswg");

    addPar ("TPBSWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::tpbswg,
      NULL,  U_VKM1, CAT_TEMP, "Temperature coefficient of Pbswg");

    addPar ("ACDE", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::acde,
      NULL,  U_MVM1, CAT_CAP, "Exponetial coefficient for charge thickness in capmod = 3 for accumulation and depletion regions");

    addPar ("MOIN", 15.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::moin,
      NULL,  U_NONE, CAT_CAP, "Coefficient for the gate-bias dependent surface potential");

    addPar ("NOFF", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::noff,
      NULL,  U_NONE, CAT_CAP, "CV parameter in Vgsteff, CV for weak to strong inversion");

    addPar ("DELVT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::delvt,
      NULL, U_VOLT, CAT_CAP, "Threshold voltage adjust for C-V");

    addPar ("KB1", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::kb1,
      NULL, U_NONE, CAT_NONE, "Scaling factor for backgate charge");

    addPar ("DLBG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dlbg,
      NULL, U_METER, CAT_CAP, "Length offset fitting parameter for backgate charge");

    addPar ("IGCMOD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::igcMod,
      NULL, U_NONE, CAT_NONE, "Gate-channel tunneling current model selector");

    addPar ("TOXQM", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::toxqm,
      NULL, U_METER, CAT_TUNNEL, "Oxide thickness for Igb calculation");

    addPar ("WTH0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wth0,
      NULL, U_METER, CAT_TEMP, "Minimum width for thermal resistance calculation");

    addPar ("RHALO", 1.0e15, false, NO_DEP,
      &MOSFET_B3SOI::Model::rhalo,
      NULL, U_OMM1, CAT_DC, "Body halo sheet resistance");

    addPar ("NTOX", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ntox,
      NULL, U_NONE, CAT_TUNNEL, "Power term of gate current");

    addPar ("TOXREF", 2.5e-9, false, NO_DEP,
      &MOSFET_B3SOI::Model::toxref,
      NULL, U_METER, CAT_TUNNEL, "Target oxide thickness");

    addPar ("EBG", 1.2, false, NO_DEP,
      &MOSFET_B3SOI::Model::ebg,
      NULL, U_VOLT, CAT_TUNNEL, "Effective bandgap in gate current calculation");

    addPar ("VEVB", 0.075, false, NO_DEP,
      &MOSFET_B3SOI::Model::vevb,
      NULL, U_NONE, CAT_TUNNEL, "Vaux parameter for valence band electron tunneling");

    addPar ("ALPHAGB1", 0.35, false, NO_DEP,
      &MOSFET_B3SOI::Model::alphaGB1,
      NULL, U_VOLTM1, CAT_TUNNEL, "First Vox dependent parameter for gate current in inversion");

    addPar ("BETAGB1", 0.03, false, NO_DEP,
      &MOSFET_B3SOI::Model::betaGB1,
      NULL, U_VOLTM2, CAT_TUNNEL, "Second Vox dependent parameter for gate current in inversion");

    addPar ("VGB1", 300.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vgb1,
      NULL, U_VOLT, CAT_TUNNEL, "Third Vox dependent parameter for gate current in inversion");

    addPar ("VECB", 0.026, false, NO_DEP,
      &MOSFET_B3SOI::Model::vecb,
      NULL, U_NONE, CAT_TUNNEL, "Vaux parameter for conduction band electron tunneling");

    addPar ("ALPHAGB2", 0.43, false, NO_DEP,
      &MOSFET_B3SOI::Model::alphaGB2,
      NULL, U_VOLTM1, CAT_TUNNEL, "First Vox dependent parameter for gate current in accumulation");

    addPar ("BETAGB2", 0.05, false, NO_DEP,
      &MOSFET_B3SOI::Model::betaGB2,
      NULL, U_VOLTM2, CAT_TUNNEL, "First Vox dependent parameter for gate current in accumulation");

    addPar ("VGB2", 17.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vgb2,
      NULL, U_VOLT, CAT_TUNNEL, "Third Vox dependent parameter for gate current in accumulation");

    addPar ("VOXH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::voxh,
      NULL, U_NONE, CAT_NONE, "The limit of Vox in gate current calculation");

    addPar ("DELTAVOX", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::deltavox,
      NULL, U_NONE, CAT_NONE, "The smoothing parameter in the Vox smoothing function");

    addPar ("AIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::aigc,
      NULL, U_FHGMHSMVM1, CAT_CURRENT, "Parameter for Igc");

    addPar ("BIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::bigc,
      NULL, U_FHGMHSMVM1, CAT_CURRENT, "Parameter for Igc");

    addPar ("CIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cigc,
      NULL, U_VOLTM1, CAT_CURRENT, "Parameter for Igc");

    addPar ("AIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::aigsd,
      NULL, U_FHGMHSMVM1, CAT_CURRENT, "Parameter for Igs,d");

    addPar ("BIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::bigsd,
      NULL, U_FHGMHSMVM1, CAT_CURRENT, "Parameter for Igs,d");

    addPar ("CIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::cigsd,
      NULL, U_VOLTM1, CAT_CURRENT, "Parameter for Igs,d");

    addPar ("NIGC", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::nigc,
      NULL, U_NONE, CAT_CURRENT, "Parameter for Igc slope");

    addPar ("PIGCD", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pigcd,
      NULL, U_NONE, CAT_CURRENT, "Parameter for Igc partition");

    addPar ("POXEDGE", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::poxedge,
      NULL, U_NONE, CAT_NONE, "Factor for the gate edge Tox");

    addPar ("DLCIG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dlcig,
      NULL, U_VOLTM1, CAT_CURRENT, "Delta L for Ig model");

    addPar ("VBS0PD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vbs0pd,
      NULL, U_NONE, CAT_NONE, "Upper bound of built-in potential lowering for FD operation");

    addPar ("VBS0FD", 0.5, false, NO_DEP,
      &MOSFET_B3SOI::Model::vbs0fd,
      NULL, U_VOLT, CAT_NONE, "Lower bound of built-in potential lowering for FD operation");

    addPar ("VBSA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vbsa,
      NULL, U_VOLT, CAT_VBI, "Offset voltage due to non-idealities");

    addPar ("NOFFFD", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::nofffd,
      NULL, U_NONE, CAT_VBI, "Smoothing parameter in FD module");

    addPar ("VOFFFD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vofffd,
      NULL, U_VOLT, CAT_VBI, "Smoothing parameter in FD module");

    addPar ("K1B", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::k1b,
      NULL, U_NONE, CAT_VBI, "First backgate body effect parameter");

    addPar ("K2B", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::k2b,
      NULL, U_NONE, CAT_VBI, "Second backgate body effect parameter for short channel effect");

    addPar ("DK2B", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dk2b,
      NULL, U_NONE, CAT_VBI, "Third backgate body effect parameter for short channel effect");

    addPar ("DVBD0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dvbd0,
      NULL, U_NONE, CAT_VBI, "First short channel effect parameter in FD module");

    addPar ("DVBD1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::dvbd1,
      NULL, U_NONE, CAT_VBI, "Second short channel effect parameter in FD module");

    addPar ("MOINFD", 1e3, false, NO_DEP,
      &MOSFET_B3SOI::Model::moinFD,
      NULL, U_NONE, CAT_VBI, "Gate bias dependance coefficient of surface potential in FD module");

    addPar ("XRCRG1", 12.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xrcrg1,
      NULL, U_NONE, CAT_RF, "Parameter for distributed channel resistance effect for intrinsic input resistance");

    addPar ("XRCRG2", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xrcrg2,
      NULL, U_NONE, CAT_RF, "Parameter to account for the excess channel diffusion resistance for intrinsic input resistance");

    addPar ("RSHG", 0.1, false, NO_DEP,
      &MOSFET_B3SOI::Model::rshg,
      NULL, U_NONE, CAT_NONE, "Gate sheet resistance");

    addPar ("NGCON", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ngcon,
      NULL, U_NONE, CAT_RF, "Number of gate contacts");

    addPar ("XGW", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xgw,
      NULL, U_METER, CAT_RF, "Distance from the gate contact to the channel edge");

    addPar ("XGL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::xgl,
      NULL, U_METER, CAT_RF, "Offset of the gate length due to variations in patterning");

    addPar ("LXJ", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lxj,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of xj");

    addPar ("LALPHAGB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lalphaGB1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ALPHAGB1");

    addPar ("LBETAGB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbetaGB1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of BETAGB1");

    addPar ("LALPHAGB2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lalphaGB2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ALPHAGB2");

    addPar ("LBETAGB2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbetaGB2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of BETAGB2");

    addPar ("LCGSL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lcgsl,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of CGSL");

    addPar ("LCGDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lcgdl,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of CGDL");

    addPar ("LCKAPPA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lckappa,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of CKAPPA");

    addPar ("LNDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lndif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NDIF");

    addPar ("LUTE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lute,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of UTE");

    addPar ("LKT1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lkt1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of KT1");

    addPar ("LKT1L", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lkt1l,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of KT1L");

    addPar ("LKT2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lkt2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of KT2");

    addPar ("LUA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lua1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of UA1");

    addPar ("LUB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lub1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of UB1");

    addPar ("LUC1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::luc1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of UC1");

    addPar ("LAT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lat,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of at");

    addPar ("LPRT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lprt,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of PRT");

    addPar ("LNTRECF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lntrecf,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NTRECF");

    addPar ("LNTRECR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lntrecr,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NTRECR");

    addPar ("LXBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lxbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of XBJT");

    addPar ("LXDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lxdif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CDIF");

    addPar ("LXREC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lxrec,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of XREC");

    addPar ("LXTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lxtun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of XTUN");

    addPar ("LAIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::laigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of AIGC");

    addPar ("LBIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of BIGC");

    addPar ("LCIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lcigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CIGC");

    addPar ("LAIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::laigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of AIGSD");

    addPar ("LBIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of BIGSD");

    addPar ("LCIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lcigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CIGSD");

    addPar ("LNIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnigc,
      NULL, U_METER, CAT_DEPENDENCY, "Length dependence of NIGC");

    addPar ("LPIGCD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lpigcd,
      NULL, U_NONE, CAT_DEPENDENCY, "Length dependence of PIGCD");

    addPar ("LPOXEDGE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lpoxedge,
      NULL, U_METER, CAT_DEPENDENCY, "Length dependence of POXEDGE");

    addPar ("LNCH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnpeak,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of nch");

    addPar ("LNSUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnsub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of nsub");

    addPar ("LNGATE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lngate,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of ngate");

    addPar ("LVTH0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvth0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of VT0");

    addPar ("LK1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lk1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of K1");

    addPar ("LK1W1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lk1w1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of K1W1");

    addPar ("LK1W2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lk1w2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of K1W2");

    addPar ("LK2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lk2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of K2");

    addPar ("LK3", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lk3,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of K3");

    addPar ("LK3B", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lk3b,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of K3B");

    addPar ("LKB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lkb1,
      NULL, U_METER, CAT_DEPENDENCY, "Length dependence of KB1");

    addPar ("LW0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lw0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of W0");

    addPar ("LNLX", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnlx,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of NLX");

    addPar ("LDVT0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldvt0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT0");

    addPar ("LDVT1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldvt1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT1");

    addPar ("LDVT2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldvt2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT2");

    addPar ("LDVT0W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldvt0w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT0W");

    addPar ("LDVT1W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldvt1w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT1W");

    addPar ("LDVT2W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldvt2w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT2W");

    addPar ("LU0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lu0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of U0");

    addPar ("LUA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lua,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of UA");

    addPar ("LUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of UB");

    addPar ("LUC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::luc,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of UC");

    addPar ("LVSAT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvsat,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of vsat");

    addPar ("LA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::la0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of a0");

    addPar ("LAGS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lags,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of ags");

    addPar ("LB0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lb0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of B0");

    addPar ("LB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lb1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of B1");

    addPar ("LKETA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lketa,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of keta");

    addPar ("LKETAS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lketas,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of KETAS");

    addPar ("LA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::la1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of a1");

    addPar ("LA2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::la2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of a2");

    addPar ("LRDSW", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lrdsw,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of RDSW");

    addPar ("LPRWB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lprwb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of PRWB");

    addPar ("LPRWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lprwg,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of PRWG");

    addPar ("LWR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lwr,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of WR");

    addPar ("LNFACTOR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnfactor,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of nfactor");

    addPar ("LDWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldwg,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DWG");

    addPar ("LDWB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldwb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DWB");

    addPar ("LVOFF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvoff,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of VOFF");

    addPar ("LETA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::leta0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of ETA0");

    addPar ("LETAB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::letab,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of ETAB");

    addPar ("LDSUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldsub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of LDSUB");

    addPar ("LCIT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lcit,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of cit");

    addPar ("LCDSC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lcdsc,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of cdsc");

    addPar ("LCDSCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lcdscb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of cdscb");

    addPar ("LCDSCD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lcdscd,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of cdscd");

    addPar ("LPCLM", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lpclm,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of PCLM");

    addPar ("LPDIBLC1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lpdibl1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of PDIBLC1");

    addPar ("LPDIBLC2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lpdibl2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of PDIBLC2");

    addPar ("LPDIBLCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lpdiblb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of PDIBLCB");

    addPar ("LDROUT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldrout,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DROUT");

    addPar ("LPVAG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lpvag,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of PVAG");

    addPar ("LDELTA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldelta,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of DELTA");

    addPar ("LALPHA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lalpha0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of ALPHA0");

    addPar ("LFBJTII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lfbjtii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of FBJTII");

    addPar ("LBETA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbeta0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of BETA0");

    addPar ("LBETA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbeta1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of BETA1");

    addPar ("LBETA2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbeta2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of BETA2");

    addPar ("LVDSATII0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvdsatii0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VDSATII0");

    addPar ("LLII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::llii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of LII");

    addPar ("LESATII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lesatii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ESATII");

    addPar ("LSII0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lsii0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of SII0");

    addPar ("LSII1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lsii1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of SII1");

    addPar ("LSII2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lsii2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of SII2");

    addPar ("LSIID", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lsiid,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of SIID");

    addPar ("LAGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lagidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of AGIDL");

    addPar ("LBGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lbgidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of BGIDL");

    addPar ("LNGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lngidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NGIDL");

    addPar ("LNTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lntun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NTUN");

    addPar ("LNDIODE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lndiode,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NDIODE");

    addPar ("LNRECF0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnrecf0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NRECF0");

    addPar ("LNRECR0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnrecr0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NRECR0");

    addPar ("LISBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lisbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ISBJT");

    addPar ("LISDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lisdif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ISDIF");

    addPar ("LISREC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lisrec,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ISREC");

    addPar ("LISTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::listun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ISTUN");

    addPar ("LVREC0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvrec0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VREC0");

    addPar ("LVTUN0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvtun0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VTUN0");

    addPar ("LNBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NBJT");

    addPar ("LLBJT0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::llbjt0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of LBJT0");

    addPar ("LVABJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvabjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VABJT");

    addPar ("LAELY", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::laely,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of AELY");

    addPar ("LAHLI", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lahli,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of AHLI");

    addPar ("LVSDFB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvsdfb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VSDFB");

    addPar ("LVSDTH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lvsdth,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VSDTH");

    addPar ("LDELVT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ldelvt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DELVT");

    addPar ("LACDE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lacde,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of ACDE");

    addPar ("LMOIN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lmoin,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of MOIN");

    addPar ("LNOFF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lnoff,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Length dependence of NOFF");

    addPar ("LXRCRG1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lxrcrg1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of XRCRG1");

    addPar ("LXRCRG2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::lxrcrg2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of XRCRG2");

    addPar ("WXJ", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wxj,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of XJ");

    addPar ("WALPHAGB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::walphaGB1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ALPHAGB1");

    addPar ("WBETAGB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wbetaGB1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of BETAGB1");

    addPar ("WALPHAGB2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::walphaGB2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ALPHAGB2");

    addPar ("WBETAGB2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wbetaGB2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of BETAGB2");

    addPar ("WCGSL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wcgsl,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of CGSL");

    addPar ("WCGDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wcgdl,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of CGDL");

    addPar ("WCKAPPA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wckappa,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of CKAPPA");

    addPar ("WNDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wndif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NDIF");

    addPar ("WUTE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wute,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of UTE");

    addPar ("WKT1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wkt1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of KT1");

    addPar ("WKT1L", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wkt1l,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of KT1L");

    addPar ("WKT2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wkt2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of KT2");

    addPar ("WUA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wua1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of UA1");

    addPar ("WUB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wub1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of UB1");

    addPar ("WUC1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wuc1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of UC1");

    addPar ("WAT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wat,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of AT");

    addPar ("WPRT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wprt,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of PRT");

    addPar ("WNTRECF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wntrecf,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NTRECF");

    addPar ("WNTRECR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wntrecr,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NTRECR");

    addPar ("WXBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wxbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of XBIT");

    addPar ("WXDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wxdif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of XDIF");

    addPar ("WXREC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wxrec,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of XREC");

    addPar ("WXTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wxtun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of XTUN");

    addPar ("WAIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::waigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of AIGC");

    addPar ("WBIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wbigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of BIGC");

    addPar ("WCIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wcigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CIGC");

    addPar ("WAIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::waigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of AIGSD");

    addPar ("WBIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wbigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of BIGSD");

    addPar ("WCIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wcigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CIGSD");

    addPar ("WNIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnigc,
      NULL, U_METER, CAT_DEPENDENCY, "Width dependence of NIGC");

    addPar ("WPIGCD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wpigcd,
      NULL, U_NONE, CAT_DEPENDENCY, "Width dependence of PIGCD");

    addPar ("WPOXEDGE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wpoxedge,
      NULL, U_METER, CAT_DEPENDENCY, "Width dependence of POXEDGE");

    addPar ("WNCH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnpeak,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of NCH");

    addPar ("WNSUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnsub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of NSUB");

    addPar ("WNGATE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wngate,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of NGATE");

    addPar ("WVTH0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvth0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of VTO");

    addPar ("WK1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wk1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of K1");

    addPar ("WK1W1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wk1w1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of K1W1");

    addPar ("WK1W2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wk1w2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of K1W2");

    addPar ("WK2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wk2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of K2");

    addPar ("WK3", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wk3,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of K3");

    addPar ("WK3B", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wk3b,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of K3B");

    addPar ("WKB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wkb1,
      NULL, U_METER, CAT_DEPENDENCY, "Width dependence of KB1");

    addPar ("WW0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ww0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of W0");

    addPar ("WNLX", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnlx,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of NLX");

    addPar ("WDVT0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdvt0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT0");

    addPar ("WDVT1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdvt1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT1");

    addPar ("WDVT2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdvt2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT2");

    addPar ("WDVT0W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdvt0w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT0W");

    addPar ("WDVT1W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdvt1w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT1W");

    addPar ("WDVT2W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdvt2w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT2W");

    addPar ("WU0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wu0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of U0");

    addPar ("WUA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wua,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of UA");

    addPar ("WUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of UB");

    addPar ("WUC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wuc,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of UC");

    addPar ("WVSAT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvsat,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of VSAT");

    addPar ("WA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wa0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of A0");

    addPar ("WAGS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wags,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of AGS");

    addPar ("WB0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wb0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of B0");

    addPar ("WB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wb1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of B1");

    addPar ("WKETA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wketa,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of KETA");

    addPar ("WKETAS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wketas,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of KETAS");

    addPar ("WA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wa1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of A1");

    addPar ("WA2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wa2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of A2");

    addPar ("WRDSW", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wrdsw,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of RDSW");

    addPar ("WPRWB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wprwb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of PRWB");

    addPar ("WPRWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wprwg,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of PRWG");

    addPar ("WWR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wwr,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of WR");

    addPar ("WNFACTOR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnfactor,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of NFACTOR");

    addPar ("WDWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdwg,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of WG");

    addPar ("WDWB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdwb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DWB");

    addPar ("WVOFF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvoff,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of VOFF");

    addPar ("WETA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::weta0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of ETA0");

    addPar ("WETAB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wetab,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of ETAB");

    addPar ("WDSUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdsub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DSUB");

    addPar ("WCIT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wcit,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of CIT");

    addPar ("WCDSC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wcdsc,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of CDSC");

    addPar ("WCDSCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wcdscb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of CDSCB");

    addPar ("WCDSCD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wcdscd,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of CDSCD");

    addPar ("WPCLM", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wpclm,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of PCLM");

    addPar ("WPDIBLC1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wpdibl1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of PDIBLC1");

    addPar ("WPDIBLC2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wpdibl2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of PDIBLC2");

    addPar ("WPDIBLCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wpdiblb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of PDIBLCB");

    addPar ("WDROUT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdrout,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DROUT");

    addPar ("WPVAG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wpvag,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of PVAG");

    addPar ("WDELTA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdelta,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of DELTA");

    addPar ("WALPHA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::walpha0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of ALPHA0");

    addPar ("WFBJTII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wfbjtii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of FBJTII");

    addPar ("WBETA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wbeta0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of BETA0");

    addPar ("WBETA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wbeta1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of BETA1");

    addPar ("WBETA2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wbeta2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of BETA2");

    addPar ("WVDSATII0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvdsatii0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VDSATII0");

    addPar ("WLII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wlii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of LII");

    addPar ("WESATII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wesatii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ESATII");

    addPar ("WSII0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wsii0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of SII0");

    addPar ("WSII1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wsii1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of SII1");

    addPar ("WSII2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wsii2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of SII2");

    addPar ("WSIID", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wsiid,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of SIID");

    addPar ("WAGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wagidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of AGIDL");

    addPar ("WBGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wbgidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of BGIDL");

    addPar ("WNGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wngidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NGIDL");

    addPar ("WNTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wntun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NTUN");

    addPar ("WNDIODE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wndiode,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NDIODE");

    addPar ("WNRECF0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnrecf0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NRECF0");

    addPar ("WNRECR0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnrecr0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NRECR0");

    addPar ("WISBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wisbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ISBJT");

    addPar ("WISDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wisdif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ISDIF");

    addPar ("WISREC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wisrec,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ISREC");

    addPar ("WISTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wistun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ISTUN");

    addPar ("WVREC0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvrec0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VREC0");

    addPar ("WVTUN0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvtun0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VTUN0");

    addPar ("WNBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NBJT");

    addPar ("WLBJT0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wlbjt0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of LBJT0");

    addPar ("WVABJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvabjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VABJT");

    addPar ("WAELY", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::waely,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of AELY");

    addPar ("WAHLI", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wahli,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of AHLI");

    addPar ("WVSDFB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvsdfb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VSDFB");

    addPar ("WVSDTH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wvsdth,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VSDTH");

    addPar ("WDELVT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wdelvt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DELVT");

    addPar ("WACDE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wacde,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of ACDE");

    addPar ("WMOIN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wmoin,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of MOIN");

    addPar ("WNOFF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wnoff,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Width dependence of NOFF");

    addPar ("WXRCRG1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wxrcrg1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of XRCRG1");

    addPar ("WXRCRG2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::wxrcrg2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of XRCRG2");

    addPar ("PXJ", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pxj,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XJ");

    addPar ("PALPHAGB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::palphaGB1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ALPHAGB1");

    addPar ("PBETAGB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pbetaGB1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BETAGB1");

    addPar ("PALPHAGB2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::palphaGB2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ALPHAGB2");

    addPar ("PBETAGB2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pbetaGB2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BETAGB2");

    addPar ("PCGSL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pcgsl,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CGSL");

    addPar ("PCGDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pcgdl,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CGDL");

    addPar ("PCKAPPA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pckappa,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CKAPPA");

    addPar ("PNDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pndif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NDIF");

    addPar ("PUTE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pute,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UTE");

    addPar ("PKT1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pkt1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KT1");

    addPar ("PKT1L", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pkt1l,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KT1L");

    addPar ("PKT2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pkt2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KT2");

    addPar ("PUA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pua1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UA1");

    addPar ("PUB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pub1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UB1");

    addPar ("PUC1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::puc1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UC1");

    addPar ("PAT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pat,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AT");

    addPar ("PPRT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pprt,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PRT");

    addPar ("PNTRECF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pntrecf,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NTRECF");

    addPar ("PNTRECR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pntrecr,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NTRECR");

    addPar ("PXBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pxbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XBJT");

    addPar ("PXDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pxdif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XDIF");

    addPar ("PXREC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pxrec,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XREC");

    addPar ("PXTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pxtun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XTUN");

    addPar ("PAIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::paigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AIGC");

    addPar ("PBIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pbigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BIGC");

    addPar ("PCIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pcigc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CIGC");

    addPar ("PAIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::paigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AIGSD");

    addPar ("PBIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pbigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BIGSD");

    addPar ("PCIGSD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pcigsd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CIGSD");

    addPar ("PNIGC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnigc,
      NULL, U_METER2, CAT_DEPENDENCY, "Cross-term dependence of NIGC");

    addPar ("PPIGCD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ppigcd,
      NULL, U_NONE, CAT_DEPENDENCY, "Cross-term dependence of PIGCD");

    addPar ("PPOXEDGE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ppoxedge,
      NULL, U_METER2, CAT_DEPENDENCY, "Cross-term dependence of OXEDGE");

    addPar ("PNCH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnpeak,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NCH");

    addPar ("PNSUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnsub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NSUB");

    addPar ("PNGATE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pngate,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NGATE");

    addPar ("PVTH0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvth0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VT0");

    addPar ("PK1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pk1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K1");

    addPar ("PK1W1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pk1w1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K1W1");

    addPar ("PK1W2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pk1w2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K1W2");

    addPar ("PK2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pk2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K2");

    addPar ("PK3", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pk3,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K3");

    addPar ("PK3B", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pk3b,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K3B");

    addPar ("PKB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pkb1,
      NULL, U_METER2, CAT_DEPENDENCY, "Cross-term dependence of KB1");

    addPar ("PW0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pw0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of W0");

    addPar ("PNLX", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnlx,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NLX");

    addPar ("PDVT0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdvt0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT0");

    addPar ("PDVT1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdvt1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT1");

    addPar ("PDVT2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdvt2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT2");

    addPar ("PDVT0W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdvt0w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT0W");

    addPar ("PDVT1W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdvt1w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT1W");

    addPar ("PDVT2W", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdvt2w,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT2W");

    addPar ("PU0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pu0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of U0");

    addPar ("PUA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pua,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UA");

    addPar ("PUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UB");

    addPar ("PUC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::puc,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UC");

    addPar ("PVSAT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvsat,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VSAT");

    addPar ("PA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pa0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of A0");

    addPar ("PAGS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pags,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AGS");

    addPar ("PB0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pb0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of B0");

    addPar ("PB1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pb1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of B1");

    addPar ("PKETA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pketa,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KETA");

    addPar ("PKETAS", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pketas,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KETAS");

    addPar ("PA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pa1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of A1");

    addPar ("PA2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pa2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of A2");

    addPar ("PRDSW", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::prdsw,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of RDSW");

    addPar ("PPRWB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pprwb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PRWB");

    addPar ("PPRWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pprwg,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PRWG");

    addPar ("PWR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pwr,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of WR");

    addPar ("PNFACTOR", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnfactor,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NFACTOR");

    addPar ("PDWG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdwg,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DWG");

    addPar ("PDWB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdwb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DWB");

    addPar ("PVOFF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvoff,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VOFF");

    addPar ("PETA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::peta0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ETA0");

    addPar ("PETAB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::petab,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ETAB");

    addPar ("PDSUB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdsub,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DSUB");

    addPar ("PCIT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pcit,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CIT");

    addPar ("PCDSC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pcdsc,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CDSC");

    addPar ("PCDSCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pcdscb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CDSCB");

    addPar ("PCDSCD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pcdscd,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CDSCD");

    addPar ("PPCLM", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ppclm,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PCLM");

    addPar ("PPDIBLC1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ppdibl1,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PDIBLC1");

    addPar ("PPDIBLC2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ppdibl2,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PDIBLC2");

    addPar ("PPDIBLCB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ppdiblb,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PDIBLCB");

    addPar ("PDROUT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdrout,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DROUT");

    addPar ("PPVAG", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::ppvag,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PVAG");

    addPar ("PDELTA", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdelta,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DELTA");

    addPar ("PALPHA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::palpha0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ALPHA0");

    addPar ("PFBJTII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pfbjtii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of FBJTII");

    addPar ("PBETA0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pbeta0,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BETA0");

    addPar ("PBETA1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pbeta1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BETA1");

    addPar ("PBETA2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pbeta2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BETA2");

    addPar ("PVDSATII0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvdsatii0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VDSATII0");

    addPar ("PLII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::plii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of LII");

    addPar ("PESATII", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pesatii,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ESATII");

    addPar ("PSII0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::psii0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of SII0");

    addPar ("PSII1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::psii1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of SII1");

    addPar ("PSII2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::psii2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of SII2");

    addPar ("PSIID", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::psiid,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of SIID");

    addPar ("PAGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pagidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AGIDL");

    addPar ("PBGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pbgidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BGIDL");

    addPar ("PNGIDL", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pngidl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NGIDL");

    addPar ("PNTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pntun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NTON");

    addPar ("PNDIODE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pndiode,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NDIODE");

    addPar ("PNRECF0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnrecf0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NRECF0");

    addPar ("PNRECR0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnrecr0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NRECR0");

    addPar ("PISBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pisbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ISBJT");

    addPar ("PISDIF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pisdif,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ISDIF");

    addPar ("PISREC", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pisrec,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ISREC");

    addPar ("PISTUN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pistun,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ISTUN");

    addPar ("PVREC0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvrec0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VREC0");

    addPar ("PVTUN0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvtun0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VTUN0");

    addPar ("PNBJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnbjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NBJT");

    addPar ("PLBJT0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::plbjt0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of LBJT0");

    addPar ("PVABJT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvabjt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VABJT");

    addPar ("PAELY", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::paely,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AELY");

    addPar ("PAHLI", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pahli,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AHLI");

    addPar ("PVSDFB", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvsdfb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VSDFB");

    addPar ("PVSDTH", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pvsdth,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VSDTH");

    addPar ("PDELVT", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pdelvt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DELVT");

    addPar ("PACDE", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pacde,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ACDE");

    addPar ("PMOIN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pmoin,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of MOIN");

    addPar ("PNOFF", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pnoff,
      NULL,  U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NOFF");

    addPar ("PXRCRG1", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pxrcrg1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XRCRG1");

    addPar ("PXRCRG2", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::pxrcrg2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XRCRG2");

    addPar ("L", 5.0e-6, false, NO_DEP,
      &MOSFET_B3SOI::Model::model_l,
      NULL, U_METER,  CAT_GEOMETRY, "Channel length");

    addPar ("W", 5.0e-6, false, NO_DEP,
      &MOSFET_B3SOI::Model::model_w,
      NULL, U_METER,  CAT_GEOMETRY, "Channel width");

    addPar ("LMAX", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lmax,
      NULL,  U_METER, CAT_BIN, "Maximum channel length");

    addPar ("LMIN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Lmin,
      NULL,  U_METER, CAT_BIN, "Minimum channel length");

    addPar ("WMAX", 1.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wmax,
      NULL,  U_METER, CAT_BIN, "Maximum channel width");

    addPar ("WMIN", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::Wmin,
      NULL,  U_METER, CAT_BIN, "Minimum channel width");

    addPar ("IGMOD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::igbMod,
      &MOSFET_B3SOI::Model::igbModGiven,
    U_NONE, CAT_TUNNEL, "Gate current model selector");

    addPar ("IGBMOD", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::igbMod,
      &MOSFET_B3SOI::Model::igbModGiven,
    U_NONE, CAT_NONE, "Flicker noise model selector");

    addPar ("VTHO", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vth0,
      &MOSFET_B3SOI::Model::vth0Given,
    U_NONE, CAT_NONE, "Threshold voltage");

    addPar ("VTH0", 0.0, false, NO_DEP,
      &MOSFET_B3SOI::Model::vth0,
      &MOSFET_B3SOI::Model::vth0Given,
    
       U_VOLT, CAT_DC, "Threshold voltage at Vbs = 0 for large L");

    // Set up non-double precision variables:
    addPar ("VERSION", "3.2", false, NO_DEP,
           &MOSFET_B3SOI::Model::version, NULL,
           
            U_NONE, CAT_CONTROL, "Version number");
    addPar ("CAPMOD", 2, false, NO_DEP,
           &MOSFET_B3SOI::Model::capMod, NULL,
           
            U_NONE, CAT_CONTROL, "Flag for capacitance models");
    addPar ("MOBMOD", 1, false, NO_DEP,
           &MOSFET_B3SOI::Model::mobMod, NULL,
           
            U_NONE, CAT_CONTROL, "Mobility model selector");
    addPar ("PARAMCHK", 0, false, NO_DEP,
           &MOSFET_B3SOI::Model::paramChk, NULL,
           
            U_NONE, CAT_CONTROL, "Parameter value check");
    addPar ("BINUNIT", 1, false, NO_DEP,
           &MOSFET_B3SOI::Model::binUnit, NULL,
           
            U_NONE, CAT_CONTROL, "Binning unit selector");
    addPar ("SOIMOD", 0, false, NO_DEP,
           &MOSFET_B3SOI::Model::soiMod, NULL,
           U_NONE, CAT_VBI, "SIO model selector, SOIMOD=0: BSIMPD, SOIMOD=1: undefined model for PD and FE, SOIMOD=2: ideal FD");
    addPar ("RGATEMOD", 0, false, NO_DEP,
           &MOSFET_B3SOI::Model::rgateMod, NULL,
           U_NONE, CAT_RF, "Gate resistance model selector");
    addPar ("BUG1830FIX", 0, false, NO_DEP,
           &MOSFET_B3SOI::Model::bug1830fix, NULL,
           U_NONE, CAT_RF, "Voltage limter fix for bug 1830");

    // Thermal model setup:
    DeviceModel::initThermalModel(*this);
}

namespace MOSFET_B3SOI {
// static jacobian maps:

vector< vector< vector<int> > > Instance::jacStamp_v;
vector< vector<int> > Instance::jacMap_v;
vector< vector< vector<int> > > Instance::jacMap2_v;

ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

// pre-processor macros.  Note:  usually, numerical constants such as
// EPSOX come from  Const.h.  The constants listed here are intended
// to be used for debugging and comparing to the original spice3f5 code.
//
// Note2: While the initial intent of these internally defined constants was
// SOI debugging, the SOI is relatively unstable device, numerically.
// For now, leaving the constants the same as in Spice3, although this
// will cause some compiler complaints.
//
// Once someone takes the time to check the SOI results and stability in
// detail, the spice3 constants will be commented out.

#define Xyce_USE_BSIMSOI_CONST 1

#ifdef Xyce_USE_BSIMSOI_CONST
#define CONSTEPSOX  (3.453133e-11)
#define CONSTEPSSI  (1.03594e-10)

#ifdef CONSTQ   // if this is defined, we need to undef it to avoid compiler warnings
#undef CONSTQ
#endif
#define CONSTQ      (1.60219e-19)

#ifdef  CONSTKoverQ  // if this is defined, we need to undef it to avoid compiler warnings
#undef CONSTKoverQ
#endif
#define CONSTKoverQ (8.617087e-5)  //  Kb / q

#ifdef CONSTEg300   // if this is defined, we need to undef it to avoid compiler warnings
#undef CONSTEg300
#endif
#define CONSTEg300  (1.115)   //  energy gap at 300K

#define CONSTboltz  (1.3806226e-23)

#ifdef M_PI   // if this is defined, we need to undef it to avoid compiler warnings
#undef M_PI
#endif
#define M_PI        (3.141592654)

#endif     // Xyce_USE_BSIMSOI_CONST

#define DELTA_1 0.02
#define DELTA_2 0.02
#define DELTA_3 0.02
// Original is 0.02, for matching IBM model, change to 0.08
#define DELTA_3_SOI 0.08
#define DELTA_4 0.02
#define DELT_Vbseff  0.005
#define DELTA_VFB  0.02
#define OFF_Vbsitf 0.02   // v3.1

#define MAX_EXPL 2.688117142e+43
#define MIN_EXPL 3.720075976e-44
#define EXPL_THRESHOLD 100.0
#define DEXP(A,B) {                                          \
        if (A > EXPL_THRESHOLD) {                            \
            B = MAX_EXPL*(1.0+(A)-EXPL_THRESHOLD);           \
        } else if (A < -EXPL_THRESHOLD)  {                   \
            B = MIN_EXPL;                                    \
        } else   {                                           \
            B = exp(A);                                      \
        }                                                    \
    }
#define CEXP(A,B,C) {                                        \
        if (A > EXPL_THRESHOLD) {                            \
            B = MAX_EXPL*(1.0+(A)-EXPL_THRESHOLD);           \
            C = MAX_EXPL;                                    \
        } else if (A < -EXPL_THRESHOLD)  {                   \
            B = MIN_EXPL;                                    \
            C = 0;                                           \
        } else   {                                           \
            B = exp(A);                                      \
            C = B;                                           \
        }                                                    \
    }

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{

  // now set the temperature related stuff.
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Instance::Instance
( InstanceBlock & IB,
  Model &model,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)
  : DeviceInstance              (IB,mlData1,ss1,ed1,do1),
    model_                              (model),
    paramPtr                              (NULL),
    dNode                                 (0),
    gNode                                 (0),
    sNode                                 (0),
    eNode                                 (0),
    pNode                                 (0),
    pNodeMappedToB                        (false),
    bNode                                 (0),
    tNode                                 (0),
    dNodePrime                            (0),
    sNodePrime                            (0),
    gNodePrime                            (0),
    gNodeMid                              (0),
    jacID                                 (0),
    ueff                                  (0.0),
    thetavth                              (0.0),
    von                                   (0.0),
    vdsat                                 (0.0),
    cgdo                                  (0.0),
    cgso                                  (0.0),
    l                                     (0.0),
    w                                     (0.0),
    numberParallel                        (0.0),
    drainArea                             (0.0),
    sourceArea                            (0.0),
    drainSquares                          (1.0),
    sourceSquares                         (1.0),
    drainPerimeter                        (0.0),
    sourcePerimeter                       (0.0),
    drainConductance                      (0.0),
    FwdSum                                (0.0),
    gbbb                                  (0.0),
    gbbdp                                 (0.0),
    gbbe                                  (0.0),
    gbbg                                  (0.0),
    gbbsp                                 (0.0),
    gbbT                                  (0.0),
    gcbdb                                 (0.0),
    gcbeb                                 (0.0),
    gcbgb                                 (0.0),
    gcbsb                                 (0.0),
    gcbT                                  (0.0),
    gcddb                                 (0.0),
    gcdeb                                 (0.0),
    gcdgb                                 (0.0),
    gcdgmb                                (0.0),
    gcdsb                                 (0.0),
    gcdT                                  (0.0),
    gcedb                                 (0.0),
    gceeb                                 (0.0),
    gcegb                                 (0.0),
    gcegmb                                (0.0),
    gcesb                                 (0.0),
    gceT                                  (0.0),
    gcgbb                                 (0.0),
    gcgdb                                 (0.0),
    gcgeb                                 (0.0),
    gcggb                                 (0.0),
    gcgmdb                                (0.0),
    gcgmeb                                (0.0),
    gcgmgmb                               (0.0),
    gcgmsb                                (0.0),
    gcgsb                                 (0.0),
    gcgT                                  (0.0),
    gcrg                                  (0.0),
    gcrgb                                 (0.0),
    gcrgd                                 (0.0),
    gcrgg                                 (0.0),
    gcrgs                                 (0.0),
    gcrg_jac                              (0.0),
    gcrgb_jac                             (0.0),
    gcrgd_jac                             (0.0),
    gcrgg_jac                             (0.0),
    gcrgs_jac                             (0.0),
    gcsdb                                 (0.0),
    gcseb                                 (0.0),
    gcsgb                                 (0.0),
    gcsgmb                                (0.0),
    gcssb                                 (0.0),
    gcsT                                  (0.0),
    gcTt                                  (0.0),
    gddpb                                 (0.0),
    gddpdp                                (0.0),
    gddpe                                 (0.0),
    gddpg                                 (0.0),
    gddpsp                                (0.0),
    gddpT                                 (0.0),
    gds                                   (0.0),
    geltd                                 (0.0),
    gigb                                  (0.0),
    gigd                                  (0.0),
    gige                                  (0.0),
    gigg                                  (0.0),
    gigs                                  (0.0),
    gigT                                  (0.0),
    gigb_jac                              (0.0),
    gigd_jac                              (0.0),
    gige_jac                              (0.0),
    gigg_jac                              (0.0),
    gigs_jac                              (0.0),
    gigT_jac                              (0.0),
    gIdtotb                               (0.0),
    gIgtotb                               (0.0),
    gIgtotd                               (0.0),
    gIgtotg                               (0.0),
    gIgtots                               (0.0),
    gIstotb                               (0.0),
    gIstotd                               (0.0),
    gIstotg                               (0.0),
    gIstots                               (0.0),
    gppb                                  (0.0),
    gppp                                  (0.0),
    gsspb                                 (0.0),
    gsspdp                                (0.0),
    gsspe                                 (0.0),
    gsspg                                 (0.0),
    gsspsp                                (0.0),
    gsspT                                 (0.0),
    gTtb                                  (0.0),
    gTtdp                                 (0.0),
    gTte                                  (0.0),
    gTtg                                  (0.0),
    gTtsp                                 (0.0),
    gTtt                                  (0.0),
    Gm                                    (0.0),
    Gmbs                                  (0.0),
    Gme                                   (0.0),
    Gmin                                  (0.0),
    GmT                                   (0.0),
    RevSum                                (0.0),
    sourceConductance                     (0.0),
    gIdtotg                               (0.0),
    gIdtotd                               (0.0),
    gIdtots                               (0.0),
    rbodyext                              (0.0),
    csesw                                 (0.0),
    dt4                                   (0.0),
    st4                                   (0.0),
    cdmin                                 (0.0),
    cdbox                                 (0.0),
    csbox                                 (0.0),
    csmin                                 (0.0),
    grgeltd                               (0.0),
    phi                                   (0.0),
    cdesw                                 (0.0),
    mode                                  (0),
    bjtoff                                (0),
    debugMod                              (0),
    OFF                                   (false),
    rth0                                  (0.0),
    cth0                                  (0.0),
    bodySquares                           (0.0),
    frbody                                (0.0),
    soiMod                                (0),
    nbc                                   (0.0),
    nseg                                  (0.0),
    pdbcp                                 (0.0),
    psbcp                                 (0.0),
    agbcp                                 (0.0),
    aebcp                                 (0.0),
    vbsusr                                (0.0),
    tnodeout                              (0),
    rgateMod                              (0),
    cdrain                                (0.0),
    gIgsg                                 (0.0),
    gIgss                                 (0.0),
    gIgcdg                                (0.0),
    gIgcds                                (0.0),
    gIgcdd                                (0.0),
    gIgcdb                                (0.0),
    Igs                                   (0.0),
    Igcd                                  (0.0),
    gIgdg                                 (0.0),
    gIgcsg                                (0.0),
    gIgdd                                 (0.0),
    gIgcss                                (0.0),
    gIgcsd                                (0.0),
    gIgcsb                                (0.0),
    Igd                                   (0.0),
    Igcs                                  (0.0),
    qinv                                  (0.0),
    cb                                    (0.0),
    cd                                    (0.0),
    cbd                                   (0.0),
    gm                                    (0.0),
    gmbs                                  (0.0),
    gbbs                                  (0.0),
    gbgs                                  (0.0),
    gbds                                  (0.0),
    cggb                                  (0.0),
    cgdb                                  (0.0),
    cgsb                                  (0.0),
    cbgb                                  (0.0),
    cbdb                                  (0.0),
    cbsb                                  (0.0),
    cdgb                                  (0.0),
    cddb                                  (0.0),
    cdsb                                  (0.0),
    qbulk                                 (0.0),
    icVDS                                 (0.0),
    icVGS                                 (0.0),
    icVBS                                 (0.0),
    icVES                                 (0.0),
    icVPS                                 (0.0),

    ChargeComputationNeeded               (true),
    selfheat                              (false),
    bodyMod                               (0),
    floating                              (0),
    dxpart                                (0.0),
    sxpart                                (0.0),
    cdreq                                 (0.0),
    ceqbd                                 (0.0),
    ceqbs                                 (0.0),
    qgdo                                  (0.0),
    qgso                                  (0.0),
    qgd                                   (0.0),
    qgs                                   (0.0),
    qge                                   (0.0),
    qgme                                  (0.0),
    qgate                                 (0.0),
    qbody                                 (0.0),
    qdrn                                  (0.0),
    qsub                                  (0.0),
    qsrc                                  (0.0),
    ceqbody                               (0.0),
    ceqgate                               (0.0),
    ceqgcrg                               (0.0),
    ceqqe                                 (0.0),
    ceqqgmid                              (0.0),
    ceqbodcon                             (0.0),
    ceqth                                 (0.0),
    ceqqth                                (0.0),
    Igtoteq                               (0.0),
    Idtoteq                               (0.0),
    Istoteq                               (0.0),
    dVgst_dVg                             (0.0),
    dVgst_dVb                             (0.0),
    dVgs_eff_dVg                          (0.0),
    dDeltaPhi_dVg                         (0.0),
    dDeltaPhi_dVd                         (0.0),
    dDeltaPhi_dVb                         (0.0),
    cqdrn                                 (0.0),
    cqgate                                (0.0),
    cqsub                                 (0.0),
    cqbody                                (0.0),
    cqtemp                                (0.0),
    vtm                                   (0.0),
    temp                                  (0.0),
    Vd                                    (0.0),
    Vg                                    (0.0),
    Vs                                    (0.0),
    Ve                                    (0.0),
    Vp                                    (0.0),
    Vb                                    (0.0),
    Vsp                                   (0.0),
    Vdp                                   (0.0),
    Vgp                                   (0.0),
    Vgm                                   (0.0),
    Idrain                                (0.0),
    Isource                               (0.0),
    Igate                                 (0.0),
    IgateMid                              (0.0),
    vgb                                   (0.0),
    vgd                                   (0.0),
    ceqqd                                 (0.0),
    ceqqb                                 (0.0),
    ceqqg                                 (0.0),
    ceqqg_Jdxp                            (0.0),
    ceqqb_Jdxp                            (0.0),
    ceqqd_Jdxp                            (0.0),
    ceqqe_Jdxp                            (0.0),
    ceqqth_Jdxp                           (0.0),
    ceqqgmid_Jdxp                         (0.0),
    cdreq_Jdxp                            (0.0),
    Idtoteq_Jdxp                          (0.0),
    Istoteq_Jdxp                          (0.0),
    Igtoteq_Jdxp                          (0.0),
    cbodcon_Jdxp                          (0.0),
    ceqbody_Jdxp                          (0.0),
    ceqgate_Jdxp                          (0.0),
    ceqbs_Jdxp                            (0.0),
    ceqbd_Jdxp                            (0.0),
    ceqbodcon_Jdxp                        (0.0),
    ceqgcrg_Jdxp                          (0.0),
    Idrain_Jdxp                           (0.0),
    Isource_Jdxp                          (0.0),
    Igate_Jdxp                            (0.0),
    IgateMid_Jdxp                         (0.0),
    ceqth_Jdxp                            (0.0),
    vbd                                   (0.0),
    vbs                                   (0.0),
    vps                                   (0.0),
    vpd                                   (0.0),
    ved                                   (0.0),
    veb                                   (0.0),
    ves                                   (0.0),
    vgs                                   (0.0),
    vge                                   (0.0),
    vds                                   (0.0),
    vged                                  (0.0),
    vgmd                                  (0.0),
    vgme                                  (0.0),
    vgmb                                  (0.0),
    vg                                    (0.0),
    vd                                    (0.0),
    vs                                    (0.0),
    vp                                    (0.0),
    ve                                    (0.0),
    deltemp                               (0.0),
    delTemp                               (0.0),
    TempRatioMinus1                       (0.0),
    vges                                  (0.0),
    vgms                                  (0.0),
    vbd_orig                              (0.0),
    vbs_orig                              (0.0),
    vps_orig                              (0.0),
    vpd_orig                              (0.0),
    ves_orig                              (0.0),
    ved_orig                              (0.0),
    vgs_orig                              (0.0),
    vds_orig                              (0.0),
    delTemp_orig                          (0.0),
    vges_orig                             (0.0),
    vgms_orig                             (0.0),
    vgd_orig                              (0.0),
    Vds                                   (0.0),
    Vgs                                   (0.0),
    Vbs                                   (0.0),
    Vbd                                   (0.0),
    Ves                                   (0.0),
    Vps                                   (0.0),
    Vds_orig                              (0.0),
    Vgs_orig                              (0.0),
    Vbs_orig                              (0.0),
    Vbd_orig                              (0.0),
    Ves_orig                              (0.0),
    Vps_orig                              (0.0),
    Vd_orig                               (0.0),
    Vg_orig                               (0.0),
    Vs_orig                               (0.0),
    Ve_orig                               (0.0),
    Vb_orig                               (0.0),
    Vp_orig                               (0.0),
    Vsp_orig                              (0.0),
    Vdp_orig                              (0.0),
    Vgp_orig                              (0.0),
    Vgm_orig                              (0.0),
    newtonIterOld                         (0),
    Vgsteff                               (0.0),
    Vdseff                                (0.0),
    ni                                    (0.0),
    Abulk                                 (0.0),
    vbseff                                (0.0),
    nstar                                 (0.0),
    rds                                   (0.0),
    AbovVgst2Vtm                          (0.0),
    ids                                   (0.0),
    igidl                                 (0.0),
    ic                                    (0.0),
    ig                                    (0.0),
    itun                                  (0.0),
    ibs                                   (0.0),
    ibd                                   (0.0),
    iii                                   (0.0),
    ibp                                   (0.0),
    gbpbs                                 (0.0),
    gbpps                                 (0.0),
    gbpT                                  (0.0),
    cbodcon                               (0.0),
    gme                                   (0.0),
    gmT                                   (0.0),
    gtempg                                (0.0),
    gtempb                                (0.0),
    gtempe                                (0.0),
    gtempT                                (0.0),
    gtempd                                (0.0),
    cth                                   (0.0),
    cjs                                   (0.0),
    cjd                                   (0.0),
    cbody                                 (0.0),
    cgate                                 (0.0),
    gjdb                                  (0.0),
    gjdd                                  (0.0),
    gjdg                                  (0.0),
    gjde                                  (0.0),
    gjdT                                  (0.0),
    gjsb                                  (0.0),
    gjsd                                  (0.0),
    gjsg                                  (0.0),
    gjsT                                  (0.0),
    gbes                                  (0.0),
    gbps                                  (0.0),
    gbT                                   (0.0),
    cgT                                   (0.0),
    cbT                                   (0.0),
    ceT                                   (0.0),
    cdT                                   (0.0),
    cbeb                                  (0.0),
    ceeb                                  (0.0),
    cdeb                                  (0.0),
    qse                                   (0.0),
    qde                                   (0.0),
    qbf                                   (0.0),
    qjs                                   (0.0),
    qjd                                   (0.0),
    cbb                                   (0.0),
    cbg                                   (0.0),
    gcse                                  (0.0),
    gcde                                  (0.0),
    qb                                    (0.0),
    qg                                    (0.0),
    qd                                    (0.0),
    qe                                    (0.0),
    qgmid                                 (0.0),
    qth                                   (0.0),
    wdiosCV_NoSwap                        (0.0),
    wdiodCV_NoSwap                        (0.0),
    CAPcgmgmb     (0.0),
    CAPcgmdb      (0.0),
    CAPcgmsb      (0.0),
    CAPcgmeb      (0.0),
    CAPcdgmb      (0.0),
    CAPcsgmb      (0.0),
    CAPcegmb      (0.0),
    CAPcggb       (0.0),
    CAPcgdb       (0.0),
    CAPcgsb       (0.0),
    CAPcgeb       (0.0),
    CAPcgbb       (0.0),
    CAPcdgb       (0.0),
    CAPcegb       (0.0),
    CAPcsgb       (0.0),
    CAPcbgb       (0.0),
    CAPcddb       (0.0),
    CAPcdsb       (0.0),
    CAPcdeb       (0.0),
    CAPcdT        (0.0),
    CAPcsdb       (0.0),
    CAPcssb       (0.0),
    CAPcseb       (0.0),
    CAPcsT        (0.0),
    CAPcgT        (0.0),
    CAPcbdb       (0.0),
    CAPcbsb       (0.0),
    CAPcbeb       (0.0),
    CAPcbT        (0.0),
    CAPcedb       (0.0),
    CAPcesb       (0.0),
    CAPceeb       (0.0),
    CAPceT        (0.0),
    CAPcTt        (0.0),
    Qeqqg         (0.0),
    Qeqqb         (0.0),
    Qeqqd         (0.0),
    Qeqqe         (0.0),
    Qeqqth        (0.0),
    Qeqqgmid      (0.0),
    Qeqqg_Jdxp    (0.0),
    Qeqqb_Jdxp    (0.0),
    Qeqqd_Jdxp    (0.0),
    Qeqqe_Jdxp    (0.0),
    Qeqqth_Jdxp    (0.0),
    Qeqqgmid_Jdxp    (0.0),
    li_store_vbd                          (0),
    li_store_vbs                          (0),
    li_store_vgs                          (0),
    li_store_vds                          (0),
    li_store_ves                          (0),
    li_store_vps                          (0),
    li_store_vg                           (0),
    li_store_vd                           (0),
    li_store_vs                           (0),
    li_store_vp                           (0),
    li_store_ve                           (0),
    li_store_deltemp                      (0),
    li_store_vges                         (0),
    li_store_vgms                         (0),
    li_store_vgp                          (0),
    li_store_vgm                          (0),
    li_store_dev_id                       (0),
    li_store_dev_ig                       (0),
    li_store_dev_is                       (0),
    li_store_dev_ie                       (0),
    li_store_dev_ib                       (0),

    li_state_qb                           (0),
    li_state_qg                           (0),
    li_state_qd                           (0),
    li_state_qe                           (0),
    li_state_qgmid                        (0),
    li_state_qth                          (0),
    li_Drain                              (0),
    li_Gate                               (0),
    li_Source                             (0),
    li_Substrate                          (0),
    li_ExtBody                            (0),
    li_Body                               (0),
    li_Temperature                        (0),
    li_DrainPrime                         (0),
    li_SourcePrime                        (0),
    li_GatePrime                          (0),
    li_GateMid                            (0),
    li_Ids                                (0),
    li_Igs                                (0),
    li_Ibs                                (0),
    li_Ies                                (0),
    li_Ips                                (0),
    ADrainEquDrainNodeOffset              (0),
    ADrainEquDrainPrimeNodeOffset         (0),
    ADrainEquIdsOffset                    (0),
    AGateEquGateNodeOffset                (0),
    AGateEquBodyNodeOffset                (0),
    AGateEquDrainPrimeNodeOffset          (0),
    AGateEquSourcePrimeNodeOffset         (0),
    AGateEquGatePrimeNodeOffset           (0),
    AGateEquGateMidNodeOffset             (0),
    AGateEquIgsOffset                     (0),
    ASourceEquSourceNodeOffset            (0),
    ASourceEquSourcePrimeNodeOffset       (0),
    ASourceEquIdsOffset                   (0),
    ASourceEquIgsOffset                   (0),
    ASourceEquIbsOffset                   (0),
    ASourceEquIesOffset                   (0),
    ASourceEquIpsOffset                   (0),
    ASubstrateEquSubstrateNodeOffset      (0),
    ASubstrateEquBodyNodeOffset           (0),
    ASubstrateEquTemperatureNodeOffset    (0),
    ASubstrateEquDrainPrimeNodeOffset     (0),
    ASubstrateEquSourcePrimeNodeOffset    (0),
    ASubstrateEquGatePrimeNodeOffset      (0),
    ASubstrateEquGateMidNodeOffset        (0),
    ASubstrateEquIesOffset                (0),
    AExtBodyEquExtBodyNodeOffset          (0),
    AExtBodyEquBodyNodeOffset             (0),
    AExtBodyEquIpsOffset                  (0),
    ABodyEquSubstrateNodeOffset           (0),
    ABodyEquExtBodyNodeOffset             (0),
    ABodyEquBodyNodeOffset                (0),
    ABodyEquTemperatureNodeOffset         (0),
    ABodyEquDrainPrimeNodeOffset          (0),
    ABodyEquSourcePrimeNodeOffset         (0),
    ABodyEquGatePrimeNodeOffset           (0),
    ABodyEquIbsOffset                     (0),
    ATemperatureEquSubstrateNodeOffset    (0),
    ATemperatureEquBodyNodeOffset         (0),
    ATemperatureEquTemperatureNodeOffset  (0),
    ATemperatureEquDrainPrimeNodeOffset   (0),
    ATemperatureEquSourcePrimeNodeOffset  (0),
    ATemperatureEquGatePrimeNodeOffset    (0),
    ADrainPrimeEquDrainNodeOffset         (0),
    ADrainPrimeEquSubstrateNodeOffset     (0),
    ADrainPrimeEquBodyNodeOffset          (0),
    ADrainPrimeEquTemperatureNodeOffset   (0),
    ADrainPrimeEquDrainPrimeNodeOffset    (0),
    ADrainPrimeEquSourcePrimeNodeOffset   (0),
    ADrainPrimeEquGatePrimeNodeOffset     (0),
    ADrainPrimeEquGateMidNodeOffset       (0),
    ASourcePrimeEquSourceNodeOffset       (0),
    ASourcePrimeEquSubstrateNodeOffset    (0),
    ASourcePrimeEquBodyNodeOffset         (0),
    ASourcePrimeEquTemperatureNodeOffset  (0),
    ASourcePrimeEquDrainPrimeNodeOffset   (0),
    ASourcePrimeEquSourcePrimeNodeOffset  (0),
    ASourcePrimeEquGatePrimeNodeOffset    (0),
    ASourcePrimeEquGateMidNodeOffset      (0),
    AGatePrimeEquGateNodeOffset           (0),
    AGatePrimeEquSubstrateNodeOffset      (0),
    AGatePrimeEquBodyNodeOffset           (0),
    AGatePrimeEquTemperatureNodeOffset    (0),
    AGatePrimeEquDrainPrimeNodeOffset     (0),
    AGatePrimeEquSourcePrimeNodeOffset    (0),
    AGatePrimeEquGatePrimeNodeOffset      (0),
    AGatePrimeEquGateMidNodeOffset        (0),
    AGateMidEquGateNodeOffset             (0),
    AGateMidEquSubstrateNodeOffset        (0),
    AGateMidEquBodyNodeOffset             (0),
    AGateMidEquDrainPrimeNodeOffset       (0),
    AGateMidEquSourcePrimeNodeOffset      (0),
    AGateMidEquGatePrimeNodeOffset        (0),
    AGateMidEquGateMidNodeOffset          (0),
    icVDSEquVdOffset                      (0),
    icVDSEquVsOffset                      (0),
    icVDSEquIdsOffset                     (0),
    icVGSEquVgOffset                      (0),
    icVGSEquVsOffset                      (0),
    icVGSEquIgsOffset                     (0),
    icVBSEquVsOffset                      (0),
    icVBSEquVbOffset                      (0),
    icVBSEquIbsOffset                     (0),
    icVESEquVeOffset                      (0),
    icVESEquVsOffset                      (0),
    icVESEquIesOffset                     (0),
    icVPSEquVpOffset                      (0),
    icVPSEquVsOffset                      (0),
    icVPSEquIpsOffset                     (0),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // Jacobian Matrix f-Ptrs:
//  drain row:
    f_DrainEquDrainNodePtr(0),
    f_DrainEquDrainPrimeNodePtr(0),
    f_DrainEquIdsPtr(0),

//  gate row:
    f_GateEquGateNodePtr(0),
    f_GateEquBodyNodePtr(0),
    f_GateEquDrainPrimeNodePtr(0),
    f_GateEquSourcePrimeNodePtr(0),
    f_GateEquGatePrimeNodePtr(0),
    f_GateEquGateMidNodePtr(0),
    f_GateEquIgsPtr(0),

//  source row:
    f_SourceEquSourceNodePtr(0),
    f_SourceEquSourcePrimeNodePtr(0),
    f_SourceEquIdsPtr(0),
    f_SourceEquIgsPtr(0),
    f_SourceEquIbsPtr(0),
    f_SourceEquIesPtr(0),
    f_SourceEquIpsPtr(0),

//  substrate row:
    f_SubstrateEquSubstrateNodePtr(0),
    f_SubstrateEquBodyNodePtr(0),
    f_SubstrateEquTemperatureNodePtr(0),
    f_SubstrateEquDrainPrimeNodePtr(0),
    f_SubstrateEquSourcePrimeNodePtr(0),
    f_SubstrateEquGatePrimeNodePtr(0),
    f_SubstrateEquGateMidNodePtr(0),
    f_SubstrateEquIesPtr(0),

// external body row:
    f_ExtBodyEquExtBodyNodePtr(0),
    f_ExtBodyEquBodyNodePtr(0),
    f_ExtBodyEquIpsPtr(0),

// body row:
    f_BodyEquSubstrateNodePtr(0),
    f_BodyEquExtBodyNodePtr(0),
    f_BodyEquBodyNodePtr(0),
    f_BodyEquTemperatureNodePtr(0),
    f_BodyEquDrainPrimeNodePtr(0),
    f_BodyEquSourcePrimeNodePtr(0),
    f_BodyEquGatePrimeNodePtr(0),
    f_BodyEquIbsPtr(0),

// temperature row:
    f_TemperatureEquSubstrateNodePtr(0),
    f_TemperatureEquBodyNodePtr(0),
    f_TemperatureEquTemperatureNodePtr(0),
    f_TemperatureEquDrainPrimeNodePtr(0),
    f_TemperatureEquSourcePrimeNodePtr(0),
    f_TemperatureEquGatePrimeNodePtr(0),

// drain' row:
    f_DrainPrimeEquDrainNodePtr(0),
    f_DrainPrimeEquSubstrateNodePtr(0),
    f_DrainPrimeEquBodyNodePtr(0),
    f_DrainPrimeEquTemperatureNodePtr(0),
    f_DrainPrimeEquDrainPrimeNodePtr(0),
    f_DrainPrimeEquSourcePrimeNodePtr(0),
    f_DrainPrimeEquGatePrimeNodePtr(0),
    f_DrainPrimeEquGateMidNodePtr(0),

// source' row:
    f_SourcePrimeEquSourceNodePtr(0),
    f_SourcePrimeEquSubstrateNodePtr(0),
    f_SourcePrimeEquBodyNodePtr(0),
    f_SourcePrimeEquTemperatureNodePtr(0),
    f_SourcePrimeEquDrainPrimeNodePtr(0),
    f_SourcePrimeEquSourcePrimeNodePtr(0),
    f_SourcePrimeEquGatePrimeNodePtr(0),
    f_SourcePrimeEquGateMidNodePtr(0),

// gate' row:
    f_GatePrimeEquGateNodePtr(0),
    f_GatePrimeEquSubstrateNodePtr(0),
    f_GatePrimeEquBodyNodePtr(0),
    f_GatePrimeEquTemperatureNodePtr(0),
    f_GatePrimeEquDrainPrimeNodePtr(0),
    f_GatePrimeEquSourcePrimeNodePtr(0),
    f_GatePrimeEquGatePrimeNodePtr(0),
    f_GatePrimeEquGateMidNodePtr(0),

// gate mid row:
    f_GateMidEquGateNodePtr(0),
    f_GateMidEquSubstrateNodePtr(0),
    f_GateMidEquBodyNodePtr(0),
    f_GateMidEquDrainPrimeNodePtr(0),
    f_GateMidEquSourcePrimeNodePtr(0),
    f_GateMidEquGatePrimeNodePtr(0),
    f_GateMidEquGateMidNodePtr(0),

    // These offset are for the voltage sources that represent initial
    // conditions on Vds, Vgs, Vbs, Ves and Vps

    // f_icVDS
    f_icVDSEquVsPtr(0),
    f_icVDSEquVdPtr(0),
    f_icVDSEquIdsPtr(0),

    // f_icVGS
    f_icVGSEquVsPtr(0),
    f_icVGSEquVgPtr(0),
    f_icVGSEquIgsPtr(0),

    // f_icVBS
    f_icVBSEquVsPtr(0),
    f_icVBSEquVbPtr(0),
    f_icVBSEquIbsPtr(0),

    // f_icVES
    f_icVESEquVsPtr(0),
    f_icVESEquVePtr(0),
    f_icVESEquIesPtr(0),

    // f_icVPS
    f_icVPSEquVsPtr(0),
    f_icVPSEquVpPtr(0),
    f_icVPSEquIpsPtr(0),

    // Jacobian Matrix q-Ptrs:

//  drain row:
    q_DrainEquDrainNodePtr(0),
    q_DrainEquDrainPrimeNodePtr(0),
    q_DrainEquIdsPtr(0),

//  gate row:
    q_GateEquGateNodePtr(0),
    q_GateEquBodyNodePtr(0),
    q_GateEquDrainPrimeNodePtr(0),
    q_GateEquSourcePrimeNodePtr(0),
    q_GateEquGatePrimeNodePtr(0),
    q_GateEquGateMidNodePtr(0),
    q_GateEquIgsPtr(0),

//  source row:
    q_SourceEquSourceNodePtr(0),
    q_SourceEquSourcePrimeNodePtr(0),
    q_SourceEquIdsPtr(0),
    q_SourceEquIgsPtr(0),
    q_SourceEquIbsPtr(0),
    q_SourceEquIesPtr(0),
    q_SourceEquIpsPtr(0),

//  substrate row:
    q_SubstrateEquSubstrateNodePtr(0),
    q_SubstrateEquBodyNodePtr(0),
    q_SubstrateEquTemperatureNodePtr(0),
    q_SubstrateEquDrainPrimeNodePtr(0),
    q_SubstrateEquSourcePrimeNodePtr(0),
    q_SubstrateEquGatePrimeNodePtr(0),
    q_SubstrateEquGateMidNodePtr(0),
    q_SubstrateEquIesPtr(0),

// external body row:
    q_ExtBodyEquExtBodyNodePtr(0),
    q_ExtBodyEquBodyNodePtr(0),
    q_ExtBodyEquIpsPtr(0),

// body row:
    q_BodyEquSubstrateNodePtr(0),
    q_BodyEquExtBodyNodePtr(0),
    q_BodyEquBodyNodePtr(0),
    q_BodyEquTemperatureNodePtr(0),
    q_BodyEquDrainPrimeNodePtr(0),
    q_BodyEquSourcePrimeNodePtr(0),
    q_BodyEquGatePrimeNodePtr(0),
    q_BodyEquIbsPtr(0),

// temperature row:
    q_TemperatureEquSubstrateNodePtr(0),
    q_TemperatureEquBodyNodePtr(0),
    q_TemperatureEquTemperatureNodePtr(0),
    q_TemperatureEquDrainPrimeNodePtr(0),
    q_TemperatureEquSourcePrimeNodePtr(0),
    q_TemperatureEquGatePrimeNodePtr(0),

// drain' row:
    q_DrainPrimeEquDrainNodePtr(0),
    q_DrainPrimeEquSubstrateNodePtr(0),
    q_DrainPrimeEquBodyNodePtr(0),
    q_DrainPrimeEquTemperatureNodePtr(0),
    q_DrainPrimeEquDrainPrimeNodePtr(0),
    q_DrainPrimeEquSourcePrimeNodePtr(0),
    q_DrainPrimeEquGatePrimeNodePtr(0),
    q_DrainPrimeEquGateMidNodePtr(0),

// source' row:
    q_SourcePrimeEquSourceNodePtr(0),
    q_SourcePrimeEquSubstrateNodePtr(0),
    q_SourcePrimeEquBodyNodePtr(0),
    q_SourcePrimeEquTemperatureNodePtr(0),
    q_SourcePrimeEquDrainPrimeNodePtr(0),
    q_SourcePrimeEquSourcePrimeNodePtr(0),
    q_SourcePrimeEquGatePrimeNodePtr(0),
    q_SourcePrimeEquGateMidNodePtr(0),

// gate' row:
    q_GatePrimeEquGateNodePtr(0),
    q_GatePrimeEquSubstrateNodePtr(0),
    q_GatePrimeEquBodyNodePtr(0),
    q_GatePrimeEquTemperatureNodePtr(0),
    q_GatePrimeEquDrainPrimeNodePtr(0),
    q_GatePrimeEquSourcePrimeNodePtr(0),
    q_GatePrimeEquGatePrimeNodePtr(0),
    q_GatePrimeEquGateMidNodePtr(0),

// gate mid row:
    q_GateMidEquGateNodePtr(0),
    q_GateMidEquSubstrateNodePtr(0),
    q_GateMidEquBodyNodePtr(0),
    q_GateMidEquDrainPrimeNodePtr(0),
    q_GateMidEquSourcePrimeNodePtr(0),
    q_GateMidEquGatePrimeNodePtr(0),
    q_GateMidEquGateMidNodePtr(0),

    // These offset are for the voltage sources that represent initial
    // conditions on Vds, Vgs, Vbs, Ves and Vps

    // q_icVDS
    q_icVDSEquVsPtr(0),
    q_icVDSEquVdPtr(0),
    q_icVDSEquIdsPtr(0),

    // q_icVGS
    q_icVGSEquVsPtr(0),
    q_icVGSEquVgPtr(0),
    q_icVGSEquIgsPtr(0),

    // q_icVBS
    q_icVBSEquVsPtr(0),
    q_icVBSEquVbPtr(0),
    q_icVBSEquIbsPtr(0),

    // q_icVES
    q_icVESEquVsPtr(0),
    q_icVESEquVePtr(0),
    q_icVESEquIesPtr(0),

    // q_icVPS
    q_icVPSEquVsPtr(0),
    q_icVPSEquVpPtr(0),
    q_icVPSEquIpsPtr(0),
#endif

    vlDebug                               (false),
    blockHomotopyID                       (0),
    randomPerturb                         (0.0)
{
  int i,j,k;

  numExtVars   = IB.numExtVars;

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = devOptions.temp.dVal();
  if (!given("L"))
    l =model_.model_l;
  if (!given("W"))
    w =model_.model_w;
  if (!given("RGATEMOD"))
    rgateMod = model_.rgateMod;
  if (!given("SOIMOD"))
    soiMod = model_.soiMod;

  if (!given("RTH0"))
    rth0 = model_.rth0;
  if (!given("CTH0"))
    cth0 = model_.cth0;

  double Qsi, Vbs0t;
  Qsi = CONSTQ * model_.npeak
      * (1.0 + model_.nlx / l) * 1e6 * model_.tsi;
  Vbs0t = 0.8 - 0.5 * Qsi / model_.csi + model_.vbsa;

  selfheat = (model_.shMod == 1) && (rth0 != 0);

  if (soiMod == 3) /* auto selection */
     if (Vbs0t > model_.vbs0fd)
        soiMod = 2; /* ideal FD mode */
     else
        if (Vbs0t < model_.vbs0pd)
           soiMod = 0; /* BSIMPD */
        else
           soiMod = 1;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  //
  if (devOptions.verboseLevel > 0 && (l > model_.Lmax || l < model_.Lmin))
  {
    ostringstream msg;
    msg << ":: instance processParams \n"
        << "channel length " << l << " out of specified range ("
        << model_.Lmin << " , " << model_.Lmax << ")" << endl;
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  if (devOptions.verboseLevel > 0 && (w > model_.Wmax || w < model_.Wmin))
  {
    ostringstream msg;
    msg << ":: instance processParams \n"
        << "channel width " << w << " out of specified range ("
        << model_.Wmin << " , " << model_.Wmax << ")" << endl;
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  // The next block of code determines which of the nodes (solution variables)
  // need to be present.  The "Node" variables, like dNode and gNode, use
  // the following convention:
  //
  // 0 - does not exist
  // 1 - exists as an internal node
  // 2 - exists as an external node
  //
  if (numExtVars < 4 || (numExtVars < 5 && given("TNODEOUT")))
  {
    string msg;
    if (given("TNODEOUT"))
      msg += "Less than 5 external nodes with tnodeout set";
    else
      msg += "Less than 4 external nodes without tnodeout set";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }
  if (numExtVars > 7 || (numExtVars > 6 && !given("TNODEOUT")))
  {
    string msg;
    if (given("TNODEOUT"))
      msg += "Over 7 nodes with tnodeout set";
    else
      msg += "Over 6 nodes without tnodeout set";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }
  numIntVars   = 0;
  floating = 0;
  dNode = 2;
  gNode = 2;
  sNode = 2;
  eNode = 2;
  bNode = 0;
  pNode = 0;
  pNodeMappedToB = false;
  tNode = 0;
  dNodePrime = 0;
  sNodePrime = 0;
  gNodePrime = 0;
  gNodeMid = 0;

  processParams ();

  if (model_.sheetResistance > 0.0)
  {
    if (sourceSquares > 0.0)
    {
      sNodePrime = 1;
      ++numIntVars;
    }
    if (drainSquares > 0.0)
    {
      dNodePrime = 1;
      ++numIntVars;
    }
  }

  if (rgateMod >= 1)
  {
    gNodePrime = 1;
    ++numIntVars;
  }

  if (rgateMod == 3)
  {
    gNodeMid = 1;
    ++numIntVars;
  }

  if (soiMod == 2) // v3.2
  {
    bodyMod = 0;
    bNode = 0;   // body node, does not exist
    pNode = 0;   // p-node(ExtBody) does not exist.
  } // For ideal FD, body contact is disabled and no body node
  else
  {
    if (numExtVars == 4) //  floating body case -- 4-node
    {
      bNode = 1; // internal
      pNode = 0; // does not exist
      ++numIntVars;
      floating = 1;
      bodyMod = 0;
    }
    else // the 5th Node has been assigned
    {
      if (!given("TNODEOUT"))
      {
        if (numExtVars == 5) // 5-node body tie, bNode has not been assigned
        {
          if ((model_.rbody == 0.0) && (model_.rbsh == 0.0))
          {                 // ideal body tie, pNode is not used
            bodyMod = 2;
            bNode = 2;   // body node, exists as external variable
            pNode = 0;   // p-node(ExtBody) does not exist.
            pNodeMappedToB = true; // but is mapped onto B
          }
          else
          {                 // nonideal body tie
            bodyMod = 1;
            bNode = 1;   // body node, exists as internal variable
            pNode = 2;   // p-node(ExtBody) exists as external variable
            ++numIntVars;
          }
        }
        else          // 6-node body tie, bNode has been assigned
        {
          bodyMod = 1;
          bNode = 2;  // body node exists and is internal
          pNode = 2;  // pnode (ExtBody) exists and is external

          // error test - make sure body resistor has a nonzero value
          if ((model_.rbody == 0.0) && (model_.rbsh == 0.0))
          {
             string msg =
             "Instance::Instance:";
             msg += "model parameter rbody=0!\n";
             N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
             model_.rbody = 1.;
          }
        }
      }
      else
      {                        // t-node assigned
        tNode = 2;
        if (numExtVars == 5)
        {                     // 4 nodes & t-node, floating body
           bNode = 1;  // internal
           pNode = 0;  // does not exist
           ++numIntVars;
           floating = 1;
           bodyMod = 0;
        }
        else
        {                     // 5 or 6 nodes & t-node, body-contact device
          if (numExtVars == 6)
          {                   // 5 nodes & tnode
            if ((model_.rbody == 0.0) && (model_.rbsh == 0.0))
            {             // ideal body tie, pNode is not used
               bNode = 2; // external
               pNode = 0; // doesn't exist
               pNodeMappedToB = true; // P and B nodes are the same
               bodyMod = 2;
            }
            else
            {             // nonideal body tie
              bNode = 1;  // internal
              pNode = 2;  // external
              ++numIntVars;
              bodyMod = 1;
            }
          }
          else
          {  // 6 nodes & t-node
            bodyMod = 1;
            bNode = 2; // external
            pNode = 2; // external

            // error test - make sure body resistor has a nonzero value
            if ((model_.rbody == 0.0) && (model_.rbsh == 0.0))
            {
               string msg =
                 "Instance::Instance:";
               msg += "model parameter rbody=0!\n";
               N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
               model_.rbody = 1.;
            }
          }
        }
      }
    }
  }

  if (selfheat && tNode == 0)
  {
    tNode = 1;
    ++numIntVars;
  }

  if (tNode == 2 && bNode == 1) {
    if (pNode == 0)
    {
      P_index = 6;
      B_index = 5;
      T_index = 4;
    }
    else
    {
      P_index = 4;
      B_index = 6;
      T_index = 5;
    }
  }
  else
  {
    P_index = 4;
    B_index = 5;
    T_index = 6;
  }

// The T node is either present or not, and is really not a normal node, but
// rather is a stand alone node that is connected to ground by a thermal
// resistance and capacitance.  It is treated specially and is the highest
// order term in the indexing of the various jacobians.

  if (tNode == 0)
  {
    i = 24;
  }
  else
  {
    if (B_index == 6)
      i = 72;
    else if (P_index == 6)
      i = 84;
    else
      i = 0;
  }

// Similarly, the P or B nodes can either be present or not.  In the fully
// depleted model they are not present.  Their presence or lack there of
// is used to indicate which chunk of jacobians/maps are to be used.

  if (pNode == 0 && bNode == 0)
  {
    j = 4;
    i = 48+i/2;
  }
  else
  {
    if (T_index == 6)
      j = 8;
    else
      j = 4;
  }

// Note that gNodeMid can only be present when gNodePrime is present so there
// are three possibilities for the gate node configuration.  The others which
// follow are the standard presence or absense of the prime node depending
// on whether a lead resistance exists.

  if (gNodeMid == 0)
  {
    i += j;
    if (gNodePrime == 0)
      i += j;
  }
  j /= 2;
  if (sNodePrime == 0)
    i += j;
  j /= 2;
  if (dNodePrime == 0)
    i += j;
  j /= 2;
  if (pNode == 0 || bNode == 0)
    i += j;

  jacID = i;

  // update numIntVars for any specified initial conditions
  if (icVDSGiven) ++numIntVars;
  if (icVGSGiven) ++numIntVars;
  if (icVBSGiven) ++numIntVars;
  if (icVESGiven) ++numIntVars;
  // Spice model documentation says to ignore icVPS if
  // there isn't a p terminal.
  if (pNode != 0)
  {
    if (icVPSGiven) ++numIntVars;
  }

  if (!given("AD")) {drainArea = devOptions.defad;}
  if (!given("AS")) {sourceArea = devOptions.defas;}

  // process source/drain series resistance
  drainConductance = model_.sheetResistance * drainSquares;
  if (drainConductance > 0.0)
    drainConductance = 1.0 / drainConductance;
  else
    drainConductance = 0.0;

  sourceConductance = model_.sheetResistance * sourceSquares;
  if (sourceConductance > 0.0)
    sourceConductance = 1.0 / sourceConductance;
  else
    sourceConductance = 0.0;

  devConMap.resize(numExtVars);
  devConMap[0] = 3;
  devConMap[1] = 2;
  devConMap[2] = 3;
  devConMap[3] = 1;
  k = 4;
  if (bNode == 2)
  {
    devConMap[k++] = 3;
    if (pNode ==2)
      devConMap[k++] = 3;
    if (tNode == 2)
      devConMap[k++] = 4;
  }
  else
  {
    if (tNode == 2)
      devConMap[k++] = 4;
    if (pNode == 2)
      devConMap[k++] = 3;
  }
  if (k != numExtVars)
  {
    cerr << endl;
    cerr << "numExtVars  = " << numExtVars << endl;
    cerr << "numIntVars  = " << numIntVars << endl;
    cerr << "dNode = " << dNode << endl;
    cerr << "gNode = " << gNode << endl;
    cerr << "sNode = " << sNode << endl;
    cerr << "eNode = " << eNode << endl;
    cerr << "bNode = " << bNode << endl;
    cerr << "pNode = " << pNode << endl;
    cerr << "tNode = " << tNode << endl;
    cerr << "dNodePrime = " << dNodePrime << endl;
    cerr << "sNodePrime = " << sNodePrime << endl;
    cerr << "gNodePrime = " << gNodePrime << endl;
    cerr << "gNodeMid = " << gNodeMid << endl;
    string msg = "Instance::Instance: Internal error in lead connectivity";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  setNumStoreVars(16);
  numStateVars = 6;
  numLeadCurrentStoreVars = 5; // drain, gate, source, ext & base lead currents
  setName(IB.getName());
  setModelName(model_.getName());

  blockHomotopyID =
    devSupport.getGainScaleBlockID(devOptions.numGainScaleBlocks);
  randomPerturb =
    devSupport.getRandomPerturbation();

  setupJacStamp();

#ifdef Xyce_DEBUG_DEVICE
  debugOutputModelParams();
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupJacStamp
// Purpose       : Sets up the static jac stamps.
//
// Special Notes : This stuff was originally in the constructor, but it is
//                 complex enough to require its own function.  Also, it
//                 may be a better choice to call this at a later time than
//                 the constructor, especially as the constructor will
//                 get called multiple times (once by the IO package, and once
//                 by topology).
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/13/08
//-----------------------------------------------------------------------------
void Instance::setupJacStamp ()
{
  int c, d, off, i, j, k, m, i_b, j_b, k_b, sz, t;
  int row1, row2, dest, n_swap, top;

// The jacobian stamps come in 96 varieties.  These are derived from the maximum
// sized jacobian and stored in the stamp and map vectors in the following
// chunks, based on certain nodes being present or missing:

// 0-23  Full jacobian
// 24-47 Temperature node missing (no self heating)
// 48-59 P and B nodes missing (fully depleted, soiMod = 2)
// 60-71 Temperature, P, and B nodes missing (no self heating, fully depleted)
// 72-83 Temperature external, Body is internal, Ext Body is external.
//       Special case where the positions of T and B must be switched so that
//       all the external nodes come before the internal nodes.
// 84-95 Temperature external, Body is internal, Ext Body is missing.
//       Special case where the positions of T and P must be switched so that
//       all the external nodes come before the internal nodes.

// Each chunk is then broken into three sub chunks based on which gate model
// is used.  The gate model defines which gate nodes are present.  The three
// possibilities are shown with the range of jacobian indicies given for the
// full jacobian base as an example:

// 0-7   Gmid and Gprime both present
// 8-15  Gprime present, Gmid merged with G node
// 16-23 Gmid and Gprime both merged with G node

// These in turn are similarly subdivided into two chunks, the higher of which
// has the Sprime node merged with S

// These in turn are similarly subdivided into two chunks, the higher of which
// has the Dprime node merged with D

// Finally, for the chunks with 24 members there is a final subdivision
// depending on whether P and B are to be merged.  For the chunks starting at
// 48, 60, and 72 P and B are not merged.  For the chunk starting at 84, P and
// B are always merged.

  if (jacStamp_v.empty() )
  {

// ----------------------------------------------------------------------------
// |NZ|      |                                                                |
// |  |      |V_d   V_g   V_s   V_e   V_p   V_b   Temp  V_d'  V_s'  V_g' V_gm |
// ----------------------------------------------------------------------------
// | 2|KCL_d | a                                         b                    |
// | 6|KCL_g |       c                       d           e     f     g    h   |
// | 2|KCL_s |             i                                   j              |
// | 7|KCL_e |                   k           l     m     n     o     p    q   |
// | 2|KCL_p |                         r     s                                |
// | 7|KCL_b |                   t     u     v     w     x     y     z        |
// | 6|Temp  |                   A           B     C     D     E     F        |
// | 8|KCL_d'| G                 H           I     J     K     L     M    N   |
// | 8|KCL_s'|             O     P           Q     R     S     T     U    V   |
// | 8|KCL_g'|       W           X           Y     Z     aa    bb    cc   dd  |
// | 7|KCL_gm|       ee          ff          gg          hh    ii    jj   kk  |
// ----------------------------------------------------------------------------

    jacMap_v.resize(96);
    jacMap2_v.resize(96);

    jacStamp_v.resize(96);
    jacStamp_v[0].resize(11);
    jacStamp_v[0][0].resize(2);
    jacStamp_v[0][0][0] = 0;
    jacStamp_v[0][0][1] = 7;
    jacStamp_v[0][1].resize(6);
    jacStamp_v[0][1][0] = 1;
    jacStamp_v[0][1][1] = 5;
    jacStamp_v[0][1][2] = 7;
    jacStamp_v[0][1][3] = 8;
    jacStamp_v[0][1][4] = 9;
    jacStamp_v[0][1][5] = 10;
    jacStamp_v[0][2].resize(2);
    jacStamp_v[0][2][0] = 2;
    jacStamp_v[0][2][1] = 8;
    jacStamp_v[0][3].resize(7);
    jacStamp_v[0][3][0] = 3;
    jacStamp_v[0][3][1] = 5;
    jacStamp_v[0][3][2] = 6;
    jacStamp_v[0][3][3] = 7;
    jacStamp_v[0][3][4] = 8;
    jacStamp_v[0][3][5] = 9;
    jacStamp_v[0][3][6] = 10;
    jacStamp_v[0][4].resize(2);
    jacStamp_v[0][4][0] = 4;
    jacStamp_v[0][4][1] = 5;
    jacStamp_v[0][5].resize(7);
    jacStamp_v[0][5][0] = 3;
    jacStamp_v[0][5][1] = 4;
    jacStamp_v[0][5][2] = 5;
    jacStamp_v[0][5][3] = 6;
    jacStamp_v[0][5][4] = 7;
    jacStamp_v[0][5][5] = 8;
    jacStamp_v[0][5][6] = 9;
    jacStamp_v[0][6].resize(6);
    jacStamp_v[0][6][0] = 3;
    jacStamp_v[0][6][1] = 5;
    jacStamp_v[0][6][2] = 6;
    jacStamp_v[0][6][3] = 7;
    jacStamp_v[0][6][4] = 8;
    jacStamp_v[0][6][5] = 9;
    jacStamp_v[0][7].resize(8);
    jacStamp_v[0][7][0] = 0;
    jacStamp_v[0][7][1] = 3;
    jacStamp_v[0][7][2] = 5;
    jacStamp_v[0][7][3] = 6;
    jacStamp_v[0][7][4] = 7;
    jacStamp_v[0][7][5] = 8;
    jacStamp_v[0][7][6] = 9;
    jacStamp_v[0][7][7] = 10;
    jacStamp_v[0][8].resize(8);
    jacStamp_v[0][8][0] = 2;
    jacStamp_v[0][8][1] = 3;
    jacStamp_v[0][8][2] = 5;
    jacStamp_v[0][8][3] = 6;
    jacStamp_v[0][8][4] = 7;
    jacStamp_v[0][8][5] = 8;
    jacStamp_v[0][8][6] = 9;
    jacStamp_v[0][8][7] = 10;
    jacStamp_v[0][9].resize(8);
    jacStamp_v[0][9][0] = 1;
    jacStamp_v[0][9][1] = 3;
    jacStamp_v[0][9][2] = 5;
    jacStamp_v[0][9][3] = 6;
    jacStamp_v[0][9][4] = 7;
    jacStamp_v[0][9][5] = 8;
    jacStamp_v[0][9][6] = 9;
    jacStamp_v[0][9][7] = 10;
    jacStamp_v[0][10].resize(7);
    jacStamp_v[0][10][0] = 1;
    jacStamp_v[0][10][1] = 3;
    jacStamp_v[0][10][2] = 5;
    jacStamp_v[0][10][3] = 7;
    jacStamp_v[0][10][4] = 8;
    jacStamp_v[0][10][5] = 9;
    jacStamp_v[0][10][6] = 10;

// Map zero now contains the maximum sized jacobian stamp.  The others will be
// generated by removal of row/cols and by using jacStampMap to combine
// row/cols corresponding to combining nodes in the usual way.

// Stamp 24 omits rol/col 6 (T) from the original

    jacStamp_v[24].resize(10);
    j = 0;
    for (i=0 ; i<11 ; ++i)
    {
      if (i != 6)
      {
        jacStamp_v[24][j].clear();
        for (k=0 ; k<jacStamp_v[0][i].size() ; ++k)
        {
          if (jacStamp_v[0][i][k] < 6)
            jacStamp_v[24][j].push_back(jacStamp_v[0][i][k]);
          else if (jacStamp_v[0][i][k] > 6)
            jacStamp_v[24][j].push_back(jacStamp_v[0][i][k]-1);
        }
        ++j;
      }
    }

// Stamp 72 swaps row/col 5 (B) with 6 (T) to make external vars come first
// Stamp 84 swaps row/col 4 (P) with 6 (T) to make external vars come first
// Stamp 73 is used as a temporary for this case

    for (c=0 ; c<2 ; ++c)
    {
      row1 = 5-c;
      row2 = 6;
      dest = 72+c;
      jacStamp_v[dest].resize(11);
      for (i=0 ; i<11 ; ++i)
      {
        if (i == row1)
          j = row2;
        else if (i == row2)
          j = row1;
        else
          j = i;
        jacStamp_v[dest][i].clear();
        for (k=0 ; k<jacStamp_v[0][j].size() ; ++k)
        {
          if (jacStamp_v[0][j][k] == row1)
            jacStamp_v[dest][i].push_back(row2);
          else if (jacStamp_v[0][j][k] == row2)
            jacStamp_v[dest][i].push_back(row1);
          else
            jacStamp_v[dest][i].push_back(jacStamp_v[0][j][k]);
        }
        n_swap = 1;
        while (n_swap > 0)
        {
          n_swap = 0;
          for (k=0 ; k<jacStamp_v[0][j].size()-1 ; ++k)
          {
            if (jacStamp_v[dest][i][k] > jacStamp_v[dest][i][k+1])
            {
              t = jacStamp_v[dest][i][k];
              jacStamp_v[dest][i][k] = jacStamp_v[dest][i][k+1];
              jacStamp_v[dest][i][k+1] = t;
              ++n_swap;
            }
          }
        }
      }
    }

// Now, 84 is generated from 73.  This is done by combining 5 (B) and 6 (P)

    jacMap_v[73].clear();
    jacStampMap (jacStamp_v[73], jacMap_v[73], jacMap2_v[73],
           jacStamp_v[84], jacMap_v[84], jacMap2_v[84], 6, 5, 11);

// Finally, clean up temporary back to virgin state in case this matters
    for (i=0 ; i<11 ; ++i)
      jacStamp_v[73][i].clear();
    jacStamp_v[73].clear();
    jacMap_v[73].clear();
    jacMap2_v[73].clear();

// Stamp 48 omits 4 (P) and 5 (B) from the original
// Stamp 60 omits 4 (P) and 5 (B) from 24, which omitted 6 (T)

    for (c=0 ; c<2 ; ++c)
    {
      jacStamp_v[48+12*c].resize(9-c);
      j = 0;
      for (i=0 ; i<11 ; ++i)
      {
        if (i != 4 && i != 5 && (c != 1 || i != 10))
        {
          jacStamp_v[48+12*c][j].clear();
          for (k=0 ; k<jacStamp_v[24*c][i].size() ; ++k)
          {
            if (jacStamp_v[24*c][i][k] < 4)
              jacStamp_v[48+12*c][j].push_back(jacStamp_v[24*c][i][k]);
            else if (jacStamp_v[24*c][i][k] > 5)
              jacStamp_v[48+12*c][j].push_back(jacStamp_v[24*c][i][k]-2);
          }
          ++j;
        }
      }
    }

    for (c=0 ; c<6 ; ++c)
    {
      if (c >= 2)
      {
        off = 48+(c-2)*12;
        d = 2;
      }
      else
      {
        off = c*24;
        d = 1;
      }
      if (c == 4)
      {
        sz = 11;
        top = sz;
      }
      else if (c == 5)
      {
        sz = 11;
        top = 10;
      }
      else
      {
        sz = 11-c;
        top = sz;
      }
      if (c != 5)
        jacMap_v[off].clear();
      for (i=0 ; i<3 ; ++i)
      {
        i_b = 8*i/d;
        if (i >= 1)
          jacStampMap (jacStamp_v[off+i_b-8/d],
                       jacMap_v[off+i_b-8/d], jacMap2_v[off+i_b-8/d],
                       jacStamp_v[off+i_b],
                       jacMap_v[off+i_b], jacMap2_v[off+i_b], top-i, 1, sz);
        for (j=0 ; j<2 ; ++j)
        {
          if (j == 1)
            jacStampMap (jacStamp_v[off+i_b],
                         jacMap_v[off+i_b], jacMap2_v[off+i_b],
                         jacStamp_v[off+i_b+4/d],
                         jacMap_v[off+i_b+4/d], jacMap2_v[off+i_b+4/d],
                         top-3, 2, sz);
          j_b = i_b+4*j/d;
          for (k=0 ; k<2 ; ++k)
          {
            if (k == 1)
            {
              jacStampMap (jacStamp_v[off+j_b],
                           jacMap_v[off+j_b], jacMap2_v[off+j_b],
                           jacStamp_v[off+j_b+2/d],
                           jacMap_v[off+j_b+2/d], jacMap2_v[off+j_b+2/d],
                           top-4, 0, sz);
            }
            if (d == 1) {
              k_b = j_b+2*k;
              for (m=0 ; m<2 ; ++m)
              {
                if (m == 1)
                  jacStampMap (jacStamp_v[off+k_b],
                               jacMap_v[off+k_b], jacMap2_v[off+k_b],
                               jacStamp_v[off+k_b+1],
                               jacMap_v[off+k_b+1], jacMap2_v[off+k_b+1],
                               5, 4, sz);
              }
            }
          }
        }
      }
    }
  }

  // If the user has specified an initial condition for this instance,
  // then we'll make an instance specific jacobian stamp (otherwise
  // we have to generate 5 times the number of jacobian stamps already
  // done above.

  if( icVDSGiven || icVGSGiven || icVBSGiven || icVESGiven || icVPSGiven )
  {
    // allocate all the rows we'll need
    int numRows = jacStamp_v[ jacID ].size();
    int numRowsMap = jacMap_v[ jacID ].size();

    // since not all the initial conditions may be present, we need
    // to calculate what the row/column numbers will be for the new
    // variables.
    int icVDSRow = 0;
    int icVGSRow = 0;
    int icVBSRow = 0;
    int icVESRow = 0;
    int icVPSRow = 0;
    int lastRow = numRows;
    int numExtraRows = 0;

    // we can update numRows for the new total number of rows
    numExtraRows += ( icVDSGiven ? 1 : 0 ) +
      ( icVGSGiven ? 1 : 0 ) +
      ( icVBSGiven ? 1 : 0 ) +
      ( icVESGiven ? 1 : 0 ) +
      ( icVPSGiven ? 1 : 0 );
    numRows += numExtraRows;
    numRowsMap += numExtraRows;

    // While we at it we'll add up how many extra terms we'll need in
    // each row
    vector< int > additionalValues( numRows, 0 );

    if( icVDSGiven )
    {
      icVDSRow = lastRow;
      ++lastRow;
      additionalValues[ 0 ]        += 1;  // KCL_d connection
      additionalValues[ 2 ]        += 1;  // KCL_s connection
      additionalValues[ icVDSRow ] += 3;  // branch equation terms
    }

    if( icVGSGiven )
    {
      icVGSRow = lastRow;
      ++lastRow;
      additionalValues[ 1 ]        += 1;  // KCL_g connection
      additionalValues[ 2 ]        += 1;  // KCL_s connection
      additionalValues[ icVGSRow ] += 3;  // branch equation terms
    }

    if( icVBSGiven )
    {
      icVBSRow = lastRow;
      ++lastRow;
      additionalValues[ 5 ]        += 1;  // KCL_b connection
      additionalValues[ 2 ]        += 1;  // KCL_s connection
      additionalValues[ icVBSRow ] += 3;  // branch equation terms
    }

    if( icVESGiven )
    {
      icVESRow = lastRow;
      ++lastRow;
      additionalValues[ 3 ]        += 1;  // KCL_e connection
      additionalValues[ 2 ]        += 1;  // KCL_s connection
      additionalValues[ icVESRow ] += 3;  // branch equation terms
    }

    if( icVPSGiven )
    {
      icVPSRow = lastRow;
      ++lastRow;
      additionalValues[ 4 ]        += 1;  // KCL_p connection
      additionalValues[ 2 ]        += 1;  // KCL_s connection
      additionalValues[ icVPSRow ] += 3;  // branch equation terms
    }

    // now resize the jacStampIC for the right number of rows.
    jacStampIC.resize( numRows );

    // resize each row and copy the normal jacstamp over
    for ( int rw=0; rw < numRows; ++rw )
    {
      int numNonZeros = 0;
      if( rw < jacStamp_v[ jacID ].size() )
      {
        numNonZeros = jacStamp_v[jacID][rw].size() +
        additionalValues[rw];

        // copy this row of the jacStamp
        jacStampIC[rw].resize( numNonZeros );
        copy(jacStamp_v[jacID][rw].begin(), jacStamp_v[jacID][rw].end(),
        jacStampIC[rw].begin() );

        // sanity check to ensure that the copy stl function didn't act
        // as an insert operation
        assert( numNonZeros == jacStampIC[rw].size() );
      }
      else
      {
        numNonZeros = additionalValues[ rw ];
        jacStampIC[rw].resize( numNonZeros );
      }
    }

    // resize each map row and copy current map over
    // can't do this during the loop above because the map's may have extra
    // rows in them if some rows and columns have been combined.
    jacMapIC.resize( numRowsMap );
    jacMapIC2.resize( numRowsMap );

    // this offset lets us easily add to the jacStamp and the maps
    // when the maps have more rows than the jacStamp.
    int mapRowOffset = jacMapIC2.size() - jacStampIC.size() ;

    for( int rw=0; rw < numRowsMap; ++rw )
    {
      int numNonZeros = 0;
      if( rw < jacMap_v[ jacID ].size() )
      {
        numNonZeros = jacMap2_v[ jacID ][rw].size();
        if( rw < jacStamp_v[ jacID ].size() )
        numNonZeros += additionalValues[ rw ];

        // copy this row of the jacStampMap2
        jacMapIC2[rw].resize( numNonZeros );
        copy(jacMap2_v[jacID][rw].begin(), jacMap2_v[jacID][rw].end(),
        jacMapIC2[rw].begin() );

        // copy jacMap_v[ jacID ][rw] as well
        jacMapIC[rw] = jacMap_v[jacID][rw];
      }
      else
      {
        numNonZeros = 3;
        jacMapIC2[rw].resize( numNonZeros );
        jacMapIC[rw] = rw  - mapRowOffset;
      }
    }

    // fix up rows with new data in them, and the map2 as well
    if( icVDSGiven )
    {
      int offset = jacStampIC[0].size() - 1;   // only one term could have
      int mapOffset = jacMapIC2[0].size() - 1; // been added here.

      // coupling to KCL_d row
      jacStampIC[0][ offset ] = icVDSRow;
      jacMapIC2 [0][ mapOffset ] = offset;

      // coupling to KCL_s row
      offset = jacStampIC[2].size() - 1
      - ( icVGSGiven ? 1 : 0 ) - ( icVBSGiven ? 1 : 0 )
      - ( icVESGiven ? 1 : 0 ) - ( icVPSGiven ? 1 : 0 );
      mapOffset = offset - jacStampIC[2].size() + jacMapIC2[2].size();
      jacStampIC[2][ offset ] = icVDSRow;
      jacMapIC2 [2][ mapOffset ] = offset;

      // extra row for icVDS
      jacStampIC[icVDSRow][0] = 0;
      jacStampIC[icVDSRow][1] = 2;
      jacStampIC[icVDSRow][2] = icVDSRow;
      jacMapIC2 [icVDSRow + mapRowOffset ][0] = 0;
      jacMapIC2 [icVDSRow + mapRowOffset ][1] = 1;
      jacMapIC2 [icVDSRow + mapRowOffset ][2] = 2;
    }

    if( icVGSGiven )
    {
      int offset = jacStampIC[1].size() - 1;   // only one term could have
      int mapOffset = jacMapIC2[1].size() - 1; // been added here.

      // coupling to KCL_g row
      jacStampIC[1][ offset ] = icVGSRow;
      jacMapIC2 [1][ mapOffset ] = offset;

      // coupling to KCL_s row
      offset = jacStampIC[2].size() - 1 - ( icVBSGiven ? 1 : 0 )
      - ( icVESGiven ? 1 : 0 ) - ( icVPSGiven ? 1 : 0 );
      mapOffset = offset - jacStampIC[2].size() + jacMapIC2[2].size();
      jacStampIC[2][ offset ] = icVGSRow;
      jacMapIC2 [2][ mapOffset ] = offset;

      // extra row for icVGS
      jacStampIC[icVGSRow][0] = 1;
      jacStampIC[icVGSRow][1] = 2;
      jacStampIC[icVGSRow][2] = icVGSRow;
      jacMapIC2 [icVGSRow + mapRowOffset ][0] = 0;
      jacMapIC2 [icVGSRow + mapRowOffset ][1] = 1;
      jacMapIC2 [icVGSRow + mapRowOffset ][2] = 2;
    }

    if( icVBSGiven )
    {
      int offset = jacStampIC[5].size() - 1;   // only one term could have
      int mapOffset = jacMapIC2[5].size() - 1;  // been added here.

      // coupling to KCL_b row
      jacStampIC[5][ offset ] = icVBSRow;
      jacMapIC2 [5][ mapOffset ] = offset;

      // coupling to KCL_s row
      offset = jacStampIC[2].size() - 1 -
      ( icVESGiven ? 1 : 0) + ( icVPSGiven ? 1 : 0 );
      mapOffset = offset - jacStampIC[2].size() + jacMapIC2[2].size();
      jacStampIC[2][ offset ] = icVBSRow;
      jacMapIC2 [2][ mapOffset ] = offset;

      // extra row for icVBS
      jacStampIC[icVBSRow][0] = 5;
      jacStampIC[icVBSRow][1] = 2;
      jacStampIC[icVBSRow][2] = icVBSRow;
      jacMapIC2 [icVBSRow + mapRowOffset ][0] = 0;
      jacMapIC2 [icVBSRow + mapRowOffset ][1] = 1;
      jacMapIC2 [icVBSRow + mapRowOffset ][2] = 2;
    }

    if( icVESGiven )
    {
      int offset = jacStampIC[3].size() - 1;    // only one term could have
      int mapOffset = jacMapIC2[3].size() - 1;  // been added here.

      // coupling to KCL_e row
      jacStampIC[3][ offset ] = icVESRow;
      jacMapIC2 [3][ mapOffset ] = offset;

      // coupling to KCL_s row
      offset = jacStampIC[2].size() - 1 - ( icVPSGiven ? 1 : 0 );
      mapOffset = offset - jacStampIC[2].size() + jacMapIC2[2].size();
      jacStampIC[2][ offset ] = icVESRow;
      jacMapIC2 [2][ mapOffset ] = offset;

      // extra row for icVES
      jacStampIC[icVESRow][0] = 3;
      jacStampIC[icVESRow][1] = 2;
      jacStampIC[icVESRow][2] = icVESRow;
      jacMapIC2 [icVESRow + mapRowOffset ][0] = 0;
      jacMapIC2 [icVESRow + mapRowOffset ][1] = 1;
      jacMapIC2 [icVESRow + mapRowOffset ][2] = 2;
    }

    if( icVPSGiven )
    {
      int offset = jacStampIC[4].size() - 1;    // only one term could have
      int mapOffset = jacMapIC2[4].size() - 1;  // been added here.

      // coupling to KCL_p row
      jacStampIC[4][ offset ] = icVPSRow;
      jacMapIC2 [4][ mapOffset ] = offset;

      // coupling to KCL_s row
      offset = jacStampIC[2].size() - 1;
      mapOffset = offset - jacStampIC[2].size() + jacMapIC2[2].size();
      jacStampIC[2][ offset ] = icVPSRow;
      jacMapIC2 [2][ mapOffset ] = offset;

      // extra row for icVPS
      jacStampIC[icVPSRow][0] = 4;
      jacStampIC[icVPSRow][1] = 2;
      jacStampIC[icVPSRow][2] = icVPSRow;
      jacMapIC2 [icVGSRow + mapRowOffset ][0] = 0;
      jacMapIC2 [icVGSRow + mapRowOffset ][1] = 1;
      jacMapIC2 [icVGSRow + mapRowOffset ][2] = 2;
    }

#ifdef Xyce_DEBUG_DEVICE
    if (devOptions.debugLevel > 0)
    {
      cout << "Original jacobian stamp before IC's added" << endl;
      for( int rw=0; rw < jacStamp_v[ jacID ].size() ; ++rw  )
      {
        cout << "jacStamp_v[ " << jacID << " ][ " << rw << "] = { " ;
        for( int cl=0; cl < jacStamp_v[ jacID ][rw].size(); ++cl )
        {
          cout << jacStamp_v[ jacID ][rw][cl];
          if( cl != (jacStamp_v[ jacID ][rw].size()-1) )
          {
            cout << ", ";
          }
        }
        cout << "}" << endl;
      }
      cout << endl;

      cout << "And as viewed through the maps"  << endl;
      for( int rw=0; rw < jacMap_v[ jacID ].size() ; ++rw  )
      {
        cout << "jacStampIC[ " << jacID << " ][ " << jacMap_v[jacID][rw]
        << "] = { " ;
        for( int cl=0; cl < jacMap2_v[ jacID ][rw].size(); ++cl )
        {
          cout << jacStamp_v[jacID][jacMap_v[jacID][rw]][jacMap2_v[jacID][rw][cl]];
          if( cl != (jacMap2_v[ jacID ][rw].size()-1) )
          {
            cout << ", ";
          }
        }
        cout << "}" << endl;
      }
      cout << endl;


      cout << "jacobian stamp including initial conditions: " << endl
      << "icVDSRow = " << icVDSRow << endl
      << "icVGSRow = " << icVGSRow << endl
      << "icVBSRow = " << icVBSRow << endl
      << "icVESRow = " << icVESRow << endl
      << "icVPSRow = " << icVPSRow << endl;
      for( int rw=0; rw < jacStampIC.size() ; ++rw  )
      {
        cout << "jacStampIC[ " << rw << "] = { " ;
        for( int cl=0; cl < jacStampIC[rw].size(); ++cl )
        {
          cout << jacStampIC[rw][cl];
          if( cl != (jacStampIC[rw].size()-1) )
          {
            cout << ", ";
          }
        }
        cout << "}" << endl;
      }
      cout << endl;

      cout << "And as viewed through the maps"  << endl;
      for( int rw=0; rw < jacMapIC.size() ; ++rw  )
      {
        cout << "jacStampIC[ " << jacMapIC[rw] << "] = { " ;
        for( int cl=0; cl < jacMapIC2[rw].size(); ++cl )
        {
          cout << jacStampIC[ jacMapIC[rw] ][ jacMapIC2[rw][cl] ];
          if( cl != (jacMapIC2[rw].size()-1) )
          {
            cout << ", ";
          }
        }
        cout << "}" << endl;
      }
      cout << endl;
    }

#endif // Xyce_DEBUG_DEVICE

  }
}

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : Instance::debugOutputModelParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/13/08
//-----------------------------------------------------------------------------
void Instance::debugOutputModelParams()
{
  if (devOptions.debugLevel > 0)
  {
    cout << "------------------------------------------------------------"
         <<endl;
    cout << "Instance: constructor: "<< getName() <<endl;
    cout << "Model Params:  " << endl;
    cout << "capMod =                       " << model_.capMod << endl;
    cout << "mobMod =                       " << model_.mobMod << endl;
    cout << "paramChk =                     " << model_.paramChk << endl;
    cout << "binUnit =                      " << model_.binUnit << endl;
    cout << "version =                      " << model_.version << endl;
    cout << "tox =                          " << model_.tox << endl;
    cout << "toxm =                         " << model_.toxm << endl;
    cout << "dtoxcv =                       " << model_.dtoxcv << endl;
    cout << "cdsc =                         " << model_.cdsc << endl;
    cout << "cdscb =                        " << model_.cdscb << endl;
    cout << "cdscd =                        " << model_.cdscd << endl;
    cout << "cit =                          " << model_.cit << endl;
    cout << "nfactor =                      " << model_.nfactor << endl;
    cout << "vsat =                         " << model_.vsat << endl;
    cout << "at =                           " << model_.at << endl;
    cout << "a0 =                           " << model_.a0 << endl;
    cout << "ags =                          " << model_.ags << endl;
    cout << "a1 =                           " << model_.a1 << endl;
    cout << "a2 =                           " << model_.a2 << endl;
    cout << "keta =                         " << model_.keta << endl;
    cout << "nsub =                         " << model_.nsub << endl;
    cout << "npeak =                        " << model_.npeak << endl;
    cout << "ngate =                        " << model_.ngate << endl;
    cout << "gamma1 =                       " << model_.gamma1 << endl;
    cout << "gamma2 =                       " << model_.gamma2 << endl;
    cout << "vbx =                          " << model_.vbx << endl;
    cout << "vbm =                          " << model_.vbm << endl;
    cout << "xt =                           " << model_.xt << endl;
    cout << "k1 =                           " << model_.k1 << endl;
    cout << "kt1 =                          " << model_.kt1 << endl;
    cout << "kt1l =                         " << model_.kt1l << endl;
    cout << "kt2 =                          " << model_.kt2 << endl;
    cout << "k2 =                           " << model_.k2 << endl;
    cout << "k3 =                           " << model_.k3 << endl;
    cout << "k3b =                          " << model_.k3b << endl;
    cout << "w0 =                           " << model_.w0 << endl;
    cout << "nlx =                          " << model_.nlx << endl;
    cout << "dvt0 =                         " << model_.dvt0 << endl;
    cout << "dvt1 =                         " << model_.dvt1 << endl;
    cout << "dvt2 =                         " << model_.dvt2 << endl;
    cout << "dvt0w =                        " << model_.dvt0w << endl;
    cout << "dvt1w =                        " << model_.dvt1w << endl;
    cout << "dvt2w =                        " << model_.dvt2w << endl;
    cout << "drout =                        " << model_.drout << endl;
    cout << "dsub =                         " << model_.dsub << endl;
    cout << "vth0=                          " << model_.vth0 << endl;
    cout << "ua =                           " << model_.ua << endl;
    cout << "ua1 =                          " << model_.ua1 << endl;
    cout << "ub =                           " << model_.ub << endl;
    cout << "ub1 =                          " << model_.ub1 << endl;
    cout << "uc =                           " << model_.uc << endl;
    cout << "uc1 =                          " << model_.uc1 << endl;
    cout << "u0 =                           " << model_.u0 << endl;
    cout << "ute =                          " << model_.ute << endl;
    cout << "voff =                         " << model_.voff << endl;
    cout << "tnom =                         " << model_.tnom << endl;
    cout << "cgso =                         " << model_.cgso << endl;
    cout << "cgdo =                         " << model_.cgdo << endl;
    cout << "xpart =                        " << model_.xpart << endl;
    cout << "delta =                        " << model_.delta << endl;
    cout << "sheetResistance =              " << model_.sheetResistance << endl;
    cout << "rdsw =                         " << model_.rdsw << endl;
    cout << "prwg =                         " << model_.prwg << endl;
    cout << "prwb =                         " << model_.prwb << endl;
    cout << "prt =                          " << model_.prt << endl;
    cout << "eta0 =                         " << model_.eta0 << endl;
    cout << "etab =                         " << model_.etab << endl;
    cout << "pclm =                         " << model_.pclm << endl;
    cout << "pdibl1 =                       " << model_.pdibl1 << endl;
    cout << "pdibl2 =                       " << model_.pdibl2 << endl;
    cout << "pdiblb =                       " << model_.pdiblb << endl;
    cout << "pvag =                         " << model_.pvag << endl;
    cout << "shMod =                        " << model_.shMod << endl;
    cout << "tbox =                         " << model_.tbox << endl;
    cout << "tsi =                          " << model_.tsi << endl;
    cout << "xj =                           " << model_.xj << endl;
    cout << "rth0 =                         " << model_.rth0 << endl;
    cout << "cth0 =                         " << model_.cth0 << endl;
    cout << "ngidl =                        " << model_.ngidl << endl;
    cout << "agidl =                        " << model_.agidl << endl;
    cout << "bgidl =                        " << model_.bgidl << endl;
    cout << "ndiode =                       " << model_.ndiode << endl;
    cout << "xbjt =                         " << model_.xbjt << endl;
    cout << "xdif =                         " << model_.xdif << endl;
    cout << "xrec =                         " << model_.xrec << endl;
    cout << "xtun =                         " << model_.xtun << endl;
    cout << "GatesidewallJctPotential =     " << model_.GatesidewallJctPotential << endl;
    cout << "bodyJctGateSideGradingCoeff =  " << model_.bodyJctGateSideGradingCoeff << endl;
    cout << "unitLengthGateSidewallJctCap = " << model_.unitLengthGateSidewallJctCap << endl;
    cout << "Lint =                         " << model_.Lint << endl;
    cout << "Ll =                           " << model_.Ll << endl;
    cout << "Llc =                          " << model_.Llc << endl;
    cout << "Lln =                          " << model_.Lln << endl;
    cout << "Lw =                           " << model_.Lw << endl;
    cout << "Lwc =                          " << model_.Lwc << endl;
    cout << "Lwn =                          " << model_.Lwn << endl;
    cout << "Lwl =                          " << model_.Lwl << endl;
    cout << "Lwlc =                         " << model_.Lwlc << endl;
    cout << "wr =                           " << model_.wr << endl;
    cout << "Wint =                         " << model_.Wint << endl;
    cout << "dwg =                          " << model_.dwg << endl;
    cout << "dwb =                          " << model_.dwb << endl;
    cout << "Wl =                           " << model_.Wl << endl;
    cout << "Wlc =                          " << model_.Wlc << endl;
    cout << "Wln =                          " << model_.Wln << endl;
    cout << "Ww =                           " << model_.Ww << endl;
    cout << "Wwc =                          " << model_.Wwc << endl;
    cout << "Wwn =                          " << model_.Wwn << endl;
    cout << "Wwl =                          " << model_.Wwl << endl;
    cout << "Wwlc =                         " << model_.Wwlc << endl;
    cout << "b0 =                           " << model_.b0 << endl;
    cout << "b1 =                           " << model_.b1 << endl;
    cout << "cgsl =                         " << model_.cgsl << endl;
    cout << "cgdl =                         " << model_.cgdl << endl;
    cout << "ckappa =                       " << model_.ckappa << endl;
    cout << "cf =                           " << model_.cf << endl;
    cout << "clc =                          " << model_.clc << endl;
    cout << "cle =                          " << model_.cle << endl;
    cout << "dwc =                          " << model_.dwc << endl;
    cout << "dlc =                          " << model_.dlc << endl;
    cout << "alpha0 =                       " << model_.alpha0 << endl;
    cout << "oxideTrapDensityA =            " << model_.oxideTrapDensityA << endl;
    cout << "oxideTrapDensityB =            " << model_.oxideTrapDensityB << endl;
    cout << "oxideTrapDensityC =            " << model_.oxideTrapDensityC << endl;
    cout << "fnoiMod =                      " << model_.fnoiMod << endl;
    cout << "tnoiMod =                      " << model_.tnoiMod << endl;
    cout << "tnoia =                        " << model_.tnoia << endl;
    cout << "tnoib =                        " << model_.tnoib << endl;
    cout << "rnoia =                        " << model_.rnoia << endl;
    cout << "rnoib =                        " << model_.rnoib << endl;
    cout << "ntnoi =                        " << model_.ntnoi << endl;
    cout << "em =                           " << model_.em << endl;
    cout << "ef =                           " << model_.ef << endl;
    cout << "af =                           " << model_.af << endl;
    cout << "kf =                           " << model_.kf << endl;
    cout << "noif =                         " << model_.noif << endl;
    cout << "k1w1 =                         " << model_.k1w1 << endl;
    cout << "k1w2 =                         " << model_.k1w2 << endl;
    cout << "ketas =                        " << model_.ketas << endl;
    cout << "dwbc =                         " << model_.dwbc << endl;
    cout << "beta0 =                        " << model_.beta0 << endl;
    cout << "beta1 =                        " << model_.beta1 << endl;
    cout << "beta2 =                        " << model_.beta2 << endl;
    cout << "vdsatii0 =                     " << model_.vdsatii0 << endl;
    cout << "tii =                          " << model_.tii << endl;
    cout << "lii =                          " << model_.lii << endl;
    cout << "sii0 =                         " << model_.sii0 << endl;
    cout << "sii1 =                         " << model_.sii1 << endl;
    cout << "sii2 =                         " << model_.sii2 << endl;
    cout << "siid =                         " << model_.siid << endl;
    cout << "fbjtii =                       " << model_.fbjtii << endl;
    cout << "esatii =                       " << model_.esatii << endl;
    cout << "ntun =                         " << model_.ntun << endl;
    cout << "nrecf0 =                       " << model_.nrecf0 << endl;
    cout << "nrecr0 =                       " << model_.nrecr0 << endl;
    cout << "isbjt =                        " << model_.isbjt << endl;
    cout << "isdif =                        " << model_.isdif << endl;
    cout << "isrec =                        " << model_.isrec << endl;
    cout << "istun =                        " << model_.istun << endl;
    cout << "ln =                           " << model_.ln << endl;
    cout << "vrec0 =                        " << model_.vrec0 << endl;
    cout << "vtun0 =                        " << model_.vtun0 << endl;
    cout << "nbjt =                         " << model_.nbjt << endl;
    cout << "lbjt0 =                        " << model_.lbjt0 << endl;
    cout << "ldif0 =                        " << model_.ldif0 << endl;
    cout << "vabjt =                        " << model_.vabjt << endl;
    cout << "aely =                         " << model_.aely << endl;
    cout << "ahli =                         " << model_.ahli << endl;
    cout << "rbody =                        " << model_.rbody << endl;
    cout << "rbsh =                         " << model_.rbsh << endl;
    cout << "cgeo =                         " << model_.cgeo << endl;
    cout << "tt =                           " << model_.tt << endl;
    cout << "ndif =                         " << model_.ndif << endl;
    cout << "vsdfb =                        " << model_.vsdfb << endl;
    cout << "vsdth =                        " << model_.vsdth << endl;
    cout << "csdmin =                       " << model_.csdmin << endl;
    cout << "asd =                          " << model_.asd << endl;
    cout << "csdesw =                       " << model_.csdesw << endl;
    cout << "ntrecf =                       " << model_.ntrecf << endl;
    cout << "ntrecr =                       " << model_.ntrecr << endl;
    cout << "dlcb =                         " << model_.dlcb << endl;
    cout << "fbody =                        " << model_.fbody << endl;
    cout << "tcjswg =                       " << model_.tcjswg << endl;
    cout << "tpbswg =                       " << model_.tpbswg << endl;
    cout << "acde =                         " << model_.acde << endl;
    cout << "moin =                         " << model_.moin << endl;
    cout << "noff =                         " << model_.noff << endl;
    cout << "delvt =                        " << model_.delvt << endl;
    cout << "kb1 =                          " << model_.kb1 << endl;
    cout << "dlbg =                         " << model_.dlbg << endl;
    cout << "igbMod =                       " << model_.igbMod << endl;
    cout << "igcMod=                        " << model_.igcMod << endl;
    cout << "toxqm =                        " << model_.toxqm << endl;
    cout << "wth0 =                         " << model_.wth0 << endl;
    cout << "rhalo =                        " << model_.rhalo << endl;
    cout << "ntox =                         " << model_.ntox << endl;
    cout << "toxref =                       " << model_.toxref << endl;
    cout << "ebg =                          " << model_.ebg << endl;
    cout << "vevb =                         " << model_.vevb << endl;
    cout << "alphaGB1 =                     " << model_.alphaGB1 << endl;
    cout << "betaGB1 =                      " << model_.betaGB1 << endl;
    cout << "vgb1 =                         " << model_.vgb1 << endl;
    cout << "vecb =                         " << model_.vecb << endl;
    cout << "alphaGB2 =                     " << model_.alphaGB2 << endl;
    cout << "betaGB2 =                      " << model_.betaGB2 << endl;
    cout << "vgb2 =                         " << model_.vgb2 << endl;
    cout << "voxh =                         " << model_.voxh << endl;
    cout << "deltavox =                     " << model_.deltavox << endl;
    cout << "aigc =                         " << model_.aigc << endl;
    cout << "bigc =                         " << model_.bigc << endl;
    cout << "cigc =                         " << model_.cigc << endl;
    cout << "aigsd =                        " << model_.aigsd << endl;
    cout << "bigsd =                        " << model_.bigsd << endl;
    cout << "cigsd =                        " << model_.cigsd << endl;
    cout << "nigc =                         " << model_.nigc << endl;
    cout << "pigcd =                        " << model_.pigcd << endl;
    cout << "poxedge =                      " << model_.poxedge << endl;
    cout << "dlcig =                        " << model_.dlcig << endl;
    cout << "soiMod =                       " << model_.soiMod << endl;
    cout << "vbs0pd =                       " << model_.vbs0pd << endl;
    cout << "vbs0fd =                       " << model_.vbs0fd << endl;
    cout << "vbsa =                         " << model_.vbsa << endl;
    cout << "nofffd =                       " << model_.nofffd << endl;
    cout << "vofffd =                       " << model_.vofffd << endl;
    cout << "k1b =                          " << model_.k1b << endl;
    cout << "k2b =                          " << model_.k2b << endl;
    cout << "dk2b =                         " << model_.dk2b << endl;
    cout << "dvbd0 =                        " << model_.dvbd0 << endl;
    cout << "dvbd1 =                        " << model_.dvbd1 << endl;
    cout << "moinFD =                       " << model_.moinFD << endl;
    cout << "rgateMod =                     " << model_.rgateMod << endl;
    cout << "xrcrg1 =                       " << model_.xrcrg1 << endl;
    cout << "xrcrg2 =                       " << model_.xrcrg2 << endl;
    cout << "rshg =                         " << model_.rshg << endl;
    cout << "ngcon =                        " << model_.ngcon << endl;
    cout << "xgw =                          " << model_.xgw << endl;
    cout << "xgl =                          " << model_.xgl << endl;
    cout << "lxj =                          " << model_.lxj << endl;
    cout << "lalphaGB1 =                    " << model_.lalphaGB1 << endl;
    cout << "lbetaGB1 =                     " << model_.lbetaGB1 << endl;
    cout << "lalphaGB2 =                    " << model_.lalphaGB2 << endl;
    cout << "lbetaGB2 =                     " << model_.lbetaGB2 << endl;
    cout << "lcgsl =                        " << model_.lcgsl << endl;
    cout << "lcgdl =                        " << model_.lcgdl << endl;
    cout << "lckappa =                      " << model_.lckappa << endl;
    cout << "lndif =                        " << model_.lndif << endl;
    cout << "lute =                         " << model_.lute << endl;
    cout << "lkt1 =                         " << model_.lkt1 << endl;
    cout << "lkt1l =                        " << model_.lkt1l << endl;
    cout << "lkt2 =                         " << model_.lkt2 << endl;
    cout << "lua1 =                         " << model_.lua1 << endl;
    cout << "lub1 =                         " << model_.lub1 << endl;
    cout << "luc1 =                         " << model_.luc1 << endl;
    cout << "lat =                          " << model_.lat << endl;
    cout << "lprt =                         " << model_.lprt << endl;
    cout << "lntrecf =                      " << model_.lntrecf << endl;
    cout << "lntrecr =                      " << model_.lntrecr << endl;
    cout << "lxbjt =                        " << model_.lxbjt << endl;
    cout << "lxdif =                        " << model_.lxdif << endl;
    cout << "lxrec =                        " << model_.lxrec << endl;
    cout << "lxtun =                        " << model_.lxtun << endl;
    cout << "laigc =                        " << model_.laigc << endl;
    cout << "lbigc =                        " << model_.lbigc << endl;
    cout << "lcigc =                        " << model_.lcigc << endl;
    cout << "laigsd =                       " << model_.laigsd << endl;
    cout << "lbigsd =                       " << model_.lbigsd << endl;
    cout << "lcigsd =                       " << model_.lcigsd << endl;
    cout << "lnigc =                        " << model_.lnigc << endl;
    cout << "lpigcd =                       " << model_.lpigcd << endl;
    cout << "lpoxedge =                     " << model_.lpoxedge << endl;
    cout << "lnpeak =                       " << model_.lnpeak << endl;
    cout << "lnsub =                        " << model_.lnsub << endl;
    cout << "lngate =                       " << model_.lngate << endl;
    cout << "lvth0 =                        " << model_.lvth0 << endl;
    cout << "lk1 =                          " << model_.lk1 << endl;
    cout << "lk1w1 =                        " << model_.lk1w1 << endl;
    cout << "lk1w2 =                        " << model_.lk1w2 << endl;
    cout << "lk2 =                          " << model_.lk2 << endl;
    cout << "lk3 =                          " << model_.lk3 << endl;
    cout << "lk3b =                         " << model_.lk3b << endl;
    cout << "lkb1 =                         " << model_.lkb1 << endl;
    cout << "lw0 =                          " << model_.lw0 << endl;
    cout << "lnlx =                         " << model_.lnlx << endl;
    cout << "ldvt0 =                        " << model_.ldvt0 << endl;
    cout << "ldvt1 =                        " << model_.ldvt1 << endl;
    cout << "ldvt2 =                        " << model_.ldvt2 << endl;
    cout << "ldvt0w =                       " << model_.ldvt0w << endl;
    cout << "ldvt1w =                       " << model_.ldvt1w << endl;
    cout << "ldvt2w =                       " << model_.ldvt2w << endl;
    cout << "lu0 =                          " << model_.lu0 << endl;
    cout << "lua =                          " << model_.lua << endl;
    cout << "lub =                          " << model_.lub << endl;
    cout << "luc =                          " << model_.luc << endl;
    cout << "lvsat =                        " << model_.lvsat << endl;
    cout << "la0 =                          " << model_.la0 << endl;
    cout << "lags =                         " << model_.lags << endl;
    cout << "lb0 =                          " << model_.lb0 << endl;
    cout << "lb1 =                          " << model_.lb1 << endl;
    cout << "lketa =                        " << model_.lketa << endl;
    cout << "lketas =                       " << model_.lketas << endl;
    cout << "la1 =                          " << model_.la1 << endl;
    cout << "la2 =                          " << model_.la2 << endl;
    cout << "lrdsw =                        " << model_.lrdsw << endl;
    cout << "lprwb =                        " << model_.lprwb << endl;
    cout << "lprwg =                        " << model_.lprwg << endl;
    cout << "lwr =                          " << model_.lwr << endl;
    cout << "lnfactor =                     " << model_.lnfactor << endl;
    cout << "ldwg =                         " << model_.ldwg << endl;
    cout << "ldwb =                         " << model_.ldwb << endl;
    cout << "lvoff =                        " << model_.lvoff << endl;
    cout << "leta0 =                        " << model_.leta0 << endl;
    cout << "letab =                        " << model_.letab << endl;
    cout << "ldsub =                        " << model_.ldsub << endl;
    cout << "lcit =                         " << model_.lcit << endl;
    cout << "lcdsc =                        " << model_.lcdsc << endl;
    cout << "lcdscb =                       " << model_.lcdscb << endl;
    cout << "lcdscd =                       " << model_.lcdscd << endl;
    cout << "lpclm =                        " << model_.lpclm << endl;
    cout << "lpdibl1 =                      " << model_.lpdibl1 << endl;
    cout << "lpdibl2 =                      " << model_.lpdibl2 << endl;
    cout << "lpdiblb =                      " << model_.lpdiblb << endl;
    cout << "ldrout =                       " << model_.ldrout << endl;
    cout << "lpvag =                        " << model_.lpvag << endl;
    cout << "ldelta =                       " << model_.ldelta << endl;
    cout << "lalpha0 =                      " << model_.lalpha0 << endl;
    cout << "lfbjtii =                      " << model_.lfbjtii << endl;
    cout << "lbeta0 =                       " << model_.lbeta0 << endl;
    cout << "lbeta1 =                       " << model_.lbeta1 << endl;
    cout << "lbeta2 =                       " << model_.lbeta2 << endl;
    cout << "lvdsatii0 =                    " << model_.lvdsatii0 << endl;
    cout << "llii =                         " << model_.llii << endl;
    cout << "lesatii =                      " << model_.lesatii << endl;
    cout << "lsii0 =                        " << model_.lsii0 << endl;
    cout << "lsii1 =                        " << model_.lsii1 << endl;
    cout << "lsii2 =                        " << model_.lsii2 << endl;
    cout << "lsiid =                        " << model_.lsiid << endl;
    cout << "lagidl =                       " << model_.lagidl << endl;
    cout << "lbgidl =                       " << model_.lbgidl << endl;
    cout << "lngidl =                       " << model_.lngidl << endl;
    cout << "lntun =                        " << model_.lntun << endl;
    cout << "lndiode =                      " << model_.lndiode << endl;
    cout << "lnrecf0 =                      " << model_.lnrecf0 << endl;
    cout << "lnrecr0 =                      " << model_.lnrecr0 << endl;
    cout << "lisbjt =                       " << model_.lisbjt << endl;
    cout << "lisdif =                       " << model_.lisdif << endl;
    cout << "lisrec =                       " << model_.lisrec << endl;
    cout << "listun =                       " << model_.listun << endl;
    cout << "lvrec0 =                       " << model_.lvrec0 << endl;
    cout << "lvtun0 =                       " << model_.lvtun0 << endl;
    cout << "lnbjt =                        " << model_.lnbjt << endl;
    cout << "llbjt0 =                       " << model_.llbjt0 << endl;
    cout << "lvabjt =                       " << model_.lvabjt << endl;
    cout << "laely =                        " << model_.laely << endl;
    cout << "lahli =                        " << model_.lahli << endl;
    cout << "lvsdfb =                       " << model_.lvsdfb << endl;
    cout << "lvsdth =                       " << model_.lvsdth << endl;
    cout << "ldelvt =                       " << model_.ldelvt << endl;
    cout << "lacde =                        " << model_.lacde << endl;
    cout << "lmoin =                        " << model_.lmoin << endl;
    cout << "lnoff =                        " << model_.lnoff << endl;
    cout << "lxrcrg1 =                      " << model_.lxrcrg1 << endl;
    cout << "lxrcrg2 =                      " << model_.lxrcrg2 << endl;
    cout << "wxj =                          " << model_.wxj << endl;
    cout << "walphaGB1 =                    " << model_.walphaGB1 << endl;
    cout << "wbetaGB1 =                     " << model_.wbetaGB1 << endl;
    cout << "walphaGB2 =                    " << model_.walphaGB2 << endl;
    cout << "wbetaGB2 =                     " << model_.wbetaGB2 << endl;
    cout << "wcgsl =                        " << model_.wcgsl << endl;
    cout << "wcgdl =                        " << model_.wcgdl << endl;
    cout << "wckappa =                      " << model_.wckappa << endl;
    cout << "wndif =                        " << model_.wndif << endl;
    cout << "wute =                         " << model_.wute << endl;
    cout << "wkt1 =                         " << model_.wkt1 << endl;
    cout << "wkt1l =                        " << model_.wkt1l << endl;
    cout << "wkt2 =                         " << model_.wkt2 << endl;
    cout << "wua1 =                         " << model_.wua1 << endl;
    cout << "wub1 =                         " << model_.wub1 << endl;
    cout << "wuc1 =                         " << model_.wuc1 << endl;
    cout << "wat =                          " << model_.wat << endl;
    cout << "wprt =                         " << model_.wprt << endl;
    cout << "wntrecf =                      " << model_.wntrecf << endl;
    cout << "wntrecr =                      " << model_.wntrecr << endl;
    cout << "wxbjt =                        " << model_.wxbjt << endl;
    cout << "wxdif =                        " << model_.wxdif << endl;
    cout << "wxrec =                        " << model_.wxrec << endl;
    cout << "wxtun =                        " << model_.wxtun << endl;
    cout << "waigc =                        " << model_.waigc << endl;
    cout << "wbigc =                        " << model_.wbigc << endl;
    cout << "wcigc =                        " << model_.wcigc << endl;
    cout << "waigsd =                       " << model_.waigsd << endl;
    cout << "wbigsd =                       " << model_.wbigsd << endl;
    cout << "wcigsd =                       " << model_.wcigsd << endl;
    cout << "wnigc =                        " << model_.wnigc << endl;
    cout << "wpigcd =                       " << model_.wpigcd << endl;
    cout << "wpoxedge =                     " << model_.wpoxedge << endl;
    cout << "wnpeak =                       " << model_.wnpeak << endl;
    cout << "wnsub =                        " << model_.wnsub << endl;
    cout << "wngate =                       " << model_.wngate << endl;
    cout << "wvth0 =                        " << model_.wvth0 << endl;
    cout << "wk1 =                          " << model_.wk1 << endl;
    cout << "wk1w1 =                        " << model_.wk1w1 << endl;
    cout << "wk1w2 =                        " << model_.wk1w2 << endl;
    cout << "wk2 =                          " << model_.wk2 << endl;
    cout << "wk3 =                          " << model_.wk3 << endl;
    cout << "wk3b =                         " << model_.wk3b << endl;
    cout << "wkb1 =                         " << model_.wkb1 << endl;
    cout << "ww0 =                          " << model_.ww0 << endl;
    cout << "wnlx =                         " << model_.wnlx << endl;
    cout << "wdvt0 =                        " << model_.wdvt0 << endl;
    cout << "wdvt1 =                        " << model_.wdvt1 << endl;
    cout << "wdvt2 =                        " << model_.wdvt2 << endl;
    cout << "wdvt0w =                       " << model_.wdvt0w << endl;
    cout << "wdvt1w =                       " << model_.wdvt1w << endl;
    cout << "wdvt2w =                       " << model_.wdvt2w << endl;
    cout << "wu0 =                          " << model_.wu0 << endl;
    cout << "wua =                          " << model_.wua << endl;
    cout << "wub =                          " << model_.wub << endl;
    cout << "wuc =                          " << model_.wuc << endl;
    cout << "wvsat =                        " << model_.wvsat << endl;
    cout << "wa0 =                          " << model_.wa0 << endl;
    cout << "wags =                         " << model_.wags << endl;
    cout << "wb0 =                          " << model_.wb0 << endl;
    cout << "wb1 =                          " << model_.wb1 << endl;
    cout << "wketa =                        " << model_.wketa << endl;
    cout << "wketas =                       " << model_.wketas << endl;
    cout << "wa1 =                          " << model_.wa1 << endl;
    cout << "wa2 =                          " << model_.wa2 << endl;
    cout << "wrdsw=                         " << model_.wrdsw << endl;
    cout << "wprwb =                        " << model_.wprwb << endl;
    cout << "wprwg =                        " << model_.wprwg << endl;
    cout << "wwr =                          " << model_.wwr << endl;
    cout << "wnfactor =                     " << model_.wnfactor << endl;
    cout << "wdwg =                         " << model_.wdwg << endl;
    cout << "wdwb =                         " << model_.wdwb << endl;
    cout << "wvoff =                        " << model_.wvoff << endl;
    cout << "weta0 =                        " << model_.weta0 << endl;
    cout << "wetab =                        " << model_.wetab << endl;
    cout << "wdsub =                        " << model_.wdsub << endl;
    cout << "wcit =                         " << model_.wcit << endl;
    cout << "wcdsc =                        " << model_.wcdsc << endl;
    cout << "wcdscb =                       " << model_.wcdscb << endl;
    cout << "wcdscd =                       " << model_.wcdscd << endl;
    cout << "wpclm =                        " << model_.wpclm << endl;
    cout << "wpdibl1 =                      " << model_.wpdibl1 << endl;
    cout << "wpdibl2 =                      " << model_.wpdibl2 << endl;
    cout << "wpdiblb =                      " << model_.wpdiblb << endl;
    cout << "wdrout =                       " << model_.wdrout << endl;
    cout << "wpvag =                        " << model_.wpvag << endl;
    cout << "wdelta =                       " << model_.wdelta << endl;
    cout << "walpha0 =                      " << model_.walpha0 << endl;
    cout << "wfbjtii =                      " << model_.wfbjtii << endl;
    cout << "wbeta0 =                       " << model_.wbeta0 << endl;
    cout << "wbeta1 =                       " << model_.wbeta1 << endl;
    cout << "wbeta2 =                       " << model_.wbeta2 << endl;
    cout << "wvdsatii0 =                    " << model_.wvdsatii0 << endl;
    cout << "wlii =                         " << model_.wlii << endl;
    cout << "wesatii =                      " << model_.wesatii << endl;
    cout << "wsii0 =                        " << model_.wsii0 << endl;
    cout << "wsii1 =                        " << model_.wsii1 << endl;
    cout << "wsii2 =                        " << model_.wsii2 << endl;
    cout << "wsiid =                        " << model_.wsiid << endl;
    cout << "wagidl =                       " << model_.wagidl << endl;
    cout << "wbgidl =                       " << model_.wbgidl << endl;
    cout << "wngidl =                       " << model_.wngidl << endl;
    cout << "wntun =                        " << model_.wntun << endl;
    cout << "wndiode =                      " << model_.wndiode << endl;
    cout << "wnrecf0 =                      " << model_.wnrecf0 << endl;
    cout << "wnrecr0 =                      " << model_.wnrecr0 << endl;
    cout << "wisbjt =                       " << model_.wisbjt << endl;
    cout << "wisdif =                       " << model_.wisdif << endl;
    cout << "wisrec =                       " << model_.wisrec << endl;
    cout << "wistun =                       " << model_.wistun << endl;
    cout << "wvrec0 =                       " << model_.wvrec0 << endl;
    cout << "wvtun0 =                       " << model_.wvtun0 << endl;
    cout << "wnbjt =                        " << model_.wnbjt << endl;
    cout << "wlbjt0 =                       " << model_.wlbjt0 << endl;
    cout << "wvabjt =                       " << model_.wvabjt << endl;
    cout << "waely =                        " << model_.waely << endl;
    cout << "wahli =                        " << model_.wahli << endl;
    cout << "wvsdfb =                       " << model_.wvsdfb << endl;
    cout << "wvsdth =                       " << model_.wvsdth << endl;
    cout << "wdelvt =                       " << model_.wdelvt << endl;
    cout << "wacde =                        " << model_.wacde << endl;
    cout << "wmoin =                        " << model_.wmoin << endl;
    cout << "wnoff =                        " << model_.wnoff << endl;
    cout << "wxrcrg1 =                      " << model_.wxrcrg1 << endl;
    cout << "wxrcrg2 =                      " << model_.wxrcrg2 << endl;
    cout << "pxj =                          " << model_.pxj << endl;
    cout << "palphaGB1 =                    " << model_.palphaGB1 << endl;
    cout << "pbetaGB1 =                     " << model_.pbetaGB1 << endl;
    cout << "palphaGB2 =                    " << model_.palphaGB2 << endl;
    cout << "pbetaGB2 =                     " << model_.pbetaGB2 << endl;
    cout << "pcgsl =                        " << model_.pcgsl << endl;
    cout << "pcgdl =                        " << model_.pcgdl << endl;
    cout << "pckappa =                      " << model_.pckappa << endl;
    cout << "pndif =                        " << model_.pndif << endl;
    cout << "pute =                         " << model_.pute << endl;
    cout << "pkt1 =                         " << model_.pkt1 << endl;
    cout << "pkt1l =                        " << model_.pkt1l << endl;
    cout << "pkt2 =                         " << model_.pkt2 << endl;
    cout << "pua1 =                         " << model_.pua1 << endl;
    cout << "pub1 =                         " << model_.pub1 << endl;
    cout << "puc1 =                         " << model_.puc1 << endl;
    cout << "pat =                          " << model_.pat << endl;
    cout << "pprt =                         " << model_.pprt << endl;
    cout << "pntrecf =                      " << model_.pntrecf << endl;
    cout << "pntrecr =                      " << model_.pntrecr << endl;
    cout << "pxbjt =                        " << model_.pxbjt << endl;
    cout << "pxdif =                        " << model_.pxdif << endl;
    cout << "pxrec =                        " << model_.pxrec << endl;
    cout << "pxtun =                        " << model_.pxtun << endl;
    cout << "paigc =                        " << model_.paigc << endl;
    cout << "pbigc =                        " << model_.pbigc << endl;
    cout << "pcigc =                        " << model_.pcigc << endl;
    cout << "paigsd =                       " << model_.paigsd << endl;
    cout << "pbigsd =                       " << model_.pbigsd << endl;
    cout << "pcigsd =                       " << model_.pcigsd << endl;
    cout << "pnigc =                        " << model_.pnigc << endl;
    cout << "ppigcd =                       " << model_.ppigcd << endl;
    cout << "ppoxedge =                     " << model_.ppoxedge << endl;
    cout << "pnpeak =                       " << model_.pnpeak << endl;
    cout << "pnsub =                        " << model_.pnsub << endl;
    cout << "pngate =                       " << model_.pngate << endl;
    cout << "pvth0 =                        " << model_.pvth0 << endl;
    cout << "pk1 =                          " << model_.pk1 << endl;
    cout << "pk1w1 =                        " << model_.pk1w1 << endl;
    cout << "pk1w2 =                        " << model_.pk1w2 << endl;
    cout << "pk2 =                          " << model_.pk2 << endl;
    cout << "pk3 =                          " << model_.pk3 << endl;
    cout << "pk3b =                         " << model_.pk3b << endl;
    cout << "pkb1 =                         " << model_.pkb1 << endl;
    cout << "pw0 =                          " << model_.pw0 << endl;
    cout << "pnlx =                         " << model_.pnlx << endl;
    cout << "pdvt0 =                        " << model_.pdvt0 << endl;
    cout << "pdvt1 =                        " << model_.pdvt1 << endl;
    cout << "pdvt2 =                        " << model_.pdvt2 << endl;
    cout << "pdvt0w =                       " << model_.pdvt0w << endl;
    cout << "pdvt1w =                       " << model_.pdvt1w << endl;
    cout << "pdvt2w =                       " << model_.pdvt2w << endl;
    cout << "pu0 =                          " << model_.pu0 << endl;
    cout << "pua =                          " << model_.pua << endl;
    cout << "pub =                          " << model_.pub << endl;
    cout << "puc =                          " << model_.puc << endl;
    cout << "pvsat =                        " << model_.pvsat << endl;
    cout << "pa0 =                          " << model_.pa0 << endl;
    cout << "pags =                         " << model_.pags << endl;
    cout << "pb0 =                          " << model_.pb0 << endl;
    cout << "pb1 =                          " << model_.pb1 << endl;
    cout << "pketa =                        " << model_.pketa << endl;
    cout << "pketas =                       " << model_.pketas << endl;
    cout << "pa1 =                          " << model_.pa1 << endl;
    cout << "pa2 =                          " << model_.pa2 << endl;
    cout << "prdsw =                        " << model_.prdsw << endl;
    cout << "pprwb =                        " << model_.pprwb << endl;
    cout << "pprwg =                        " << model_.pprwg << endl;
    cout << "pwr =                          " << model_.pwr << endl;
    cout << "pnfactor =                     " << model_.pnfactor << endl;
    cout << "pdwg =                         " << model_.pdwg << endl;
    cout << "pdwb =                         " << model_.pdwb << endl;
    cout << "pvoff=                         " << model_.pvoff << endl;
    cout << "peta0 =                        " << model_.peta0 << endl;
    cout << "petab =                        " << model_.petab << endl;
    cout << "pdsub =                        " << model_.pdsub << endl;
    cout << "pcit =                         " << model_.pcit << endl;
    cout << "pcdsc =                        " << model_.pcdsc << endl;
    cout << "pcdscb =                       " << model_.pcdscb << endl;
    cout << "pcdscd =                       " << model_.pcdscd << endl;
    cout << "ppclm =                        " << model_.ppclm << endl;
    cout << "ppdibl1 =                      " << model_.ppdibl1 << endl;
    cout << "ppdibl2 =                      " << model_.ppdibl2 << endl;
    cout << "ppdiblb =                      " << model_.ppdiblb << endl;
    cout << "pdrout =                       " << model_.pdrout << endl;
    cout << "ppvag =                        " << model_.ppvag << endl;
    cout << "pdelta =                       " << model_.pdelta << endl;
    cout << "palpha0 =                      " << model_.palpha0 << endl;
    cout << "pfbjtii =                      " << model_.pfbjtii << endl;
    cout << "pbeta0 =                       " << model_.pbeta0 << endl;
    cout << "pbeta1 =                       " << model_.pbeta1 << endl;
    cout << "pbeta2 =                       " << model_.pbeta2 << endl;
    cout << "pvdsatii0 =                    " << model_.pvdsatii0 << endl;
    cout << "plii =                         " << model_.plii << endl;
    cout << "pesatii =                      " << model_.pesatii << endl;
    cout << "psii0 =                        " << model_.psii0 << endl;
    cout << "psii1 =                        " << model_.psii1 << endl;
    cout << "psii2 =                        " << model_.psii2 << endl;
    cout << "psiid =                        " << model_.psiid << endl;
    cout << "pagidl =                       " << model_.pagidl << endl;
    cout << "pbgidl =                       " << model_.pbgidl << endl;
    cout << "pngidl=                        " << model_.pngidl << endl;
    cout << "pntun =                        " << model_.pntun << endl;
    cout << "pndiode =                      " << model_.pndiode << endl;
    cout << "pnrecf0 =                      " << model_.pnrecf0 << endl;
    cout << "pnrecr0 =                      " << model_.pnrecr0 << endl;
    cout << "pisbjt =                       " << model_.pisbjt << endl;
    cout << "pisdif =                       " << model_.pisdif << endl;
    cout << "pisrec =                       " << model_.pisrec << endl;
    cout << "pistun =                       " << model_.pistun << endl;
    cout << "pvrec0 =                       " << model_.pvrec0 << endl;
    cout << "pvtun0 =                       " << model_.pvtun0 << endl;
    cout << "pnbjt =                        " << model_.pnbjt << endl;
    cout << "plbjt0 =                       " << model_.plbjt0 << endl;
    cout << "pvabjt =                       " << model_.pvabjt << endl;
    cout << "paely =                        " << model_.paely << endl;
    cout << "pahli =                        " << model_.pahli << endl;
    cout << "pvsdfb =                       " << model_.pvsdfb << endl;
    cout << "pvsdth =                       " << model_.pvsdth << endl;
    cout << "pdelvt =                       " << model_.pdelvt << endl;
    cout << "pacde =                        " << model_.pacde << endl;
    cout << "pmoin =                        " << model_.pmoin << endl;
    cout << "pnoff =                        " << model_.pnoff << endl;
    cout << "pxrcrg1 =                      " << model_.pxrcrg1 << endl;
    cout << "pxrcrg2 =                      " << model_.pxrcrg2 << endl;

    cout << endl;
    cout << "Instance Params: " << endl;
    cout << "l:               " << l << endl;
    cout << "w:               " << w << endl;
    cout << "drainArea:       " << drainArea << endl;
    cout << "sourceArea:      " << sourceArea << endl;
    cout << "drainSquares:    " << drainSquares << endl;
    cout << "sourceSquares:   " << sourceSquares << endl;
    cout << "drainPerimeter:  " << drainPerimeter << endl;
    cout << "sourcePerimeter: " << sourcePerimeter << endl;
    cout << "icVBS:           " << icVBS << endl;
    cout << "icVDS:           " << icVDS << endl;
    cout << "icVGS:           " << icVGS << endl;
    cout << "bjtoff:          " << bjtoff << endl;
    cout << "debugMod:        " << debugMod << endl;
    cout << "rth0:            " << rth0 << endl;
    cout << "cth0:            " << cth0 << endl;
    cout << "bodySquares:     " << bodySquares << endl;
    cout << "frbody:          " << frbody << endl;
    cout << "soiMod:          " << soiMod << endl;
    cout << "nbc:             " << nbc << endl;
    cout << "nseg:            " << nseg << endl;
    cout << "pdbcp:           " << pdbcp << endl;
    cout << "psbcp:           " << psbcp << endl;
    cout << "agbcp:           " << agbcp << endl;
    cout << "aebcp:           " << aebcp << endl;
    cout << "vbsusr:          " << vbsusr << endl;
    cout << "tnodeout:        " << tnodeout << endl;
    cout << "rgateMod:        " << rgateMod << endl;
    cout << "numberParallel:  " << numberParallel << endl;
  }
}
#endif


//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
#ifndef REUSE_PARAMPTR
  if (paramPtr != static_cast<SizeDependParam *> (NULL))
    delete paramPtr;
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                                            const vector<int> & extLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (devOptions.debugLevel > 0)
  {
    cout << dashedline << endl;
    cout << "  In Instance::register LIDs\n\n";
    cout << "  name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "  number of internal variables: " << numInt << endl;
    cout << "  number of external variables: " << numExt << endl;

    int i1;
    for (i1=0;i1<numInt;++i1)
    {
      cout << "int["<<i1<<"] = " << intLIDVecRef[i1] << endl;
    }
    cout << endl;
    for (i1=0;i1<numExt;++i1)
    {
      cout << "ext["<<i1<<"] = " << extLIDVecRef[i1] << endl;
    }
  }
#endif

  numIntVars = 0;
  if (tNode == 1)
    ++numIntVars;
  if (bNode == 1)
    ++numIntVars;
  if (dNodePrime == 1)
    ++numIntVars;
  if (sNodePrime == 1)
    ++numIntVars;
  if (gNodePrime == 1)
    ++numIntVars;
  if (gNodeMid == 1)
    ++numIntVars;

  if (icVDSGiven) ++numIntVars;
  if (icVGSGiven) ++numIntVars;
  if (icVBSGiven) ++numIntVars;
  if (icVESGiven) ++numIntVars;
  // Spice model documentation says to ignore icVPS if
  // there isn't a p terminal.
  if (pNode != 0)
  {
    if (icVPSGiven) ++numIntVars;
  }

  if ( numIntVars !=  numInt )
  {
    cout << numIntVars << " != " << numInt << endl;
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  int intLoc = 0;
  int extLoc = 0;

  li_Drain = extLIDVec[extLoc++];
  li_Gate = extLIDVec[extLoc++];
  li_Source = extLIDVec[extLoc++];
  li_Substrate = extLIDVec[extLoc++];

  // Note: In SPICE, the "extBody" is called the p-node.
  if (pNode == 2)
    li_ExtBody = extLIDVec[extLoc++];
  else
    li_ExtBody = -1;

  if (bNode == 1)
    li_Body = intLIDVec[intLoc++];
  else if (bNode == 2)
    li_Body = extLIDVec[extLoc++];
  else
    li_Body = -1;

  // correct mistake in original Xyce implementation.
  // Use of orthogonal flag necessary to avoid impacting jacstamp set-up.
  if (pNodeMappedToB)
    li_ExtBody = li_Body;

  if (tNode == 1)
    li_Temperature = intLIDVec[intLoc++];
  else if (tNode == 2)
    li_Temperature = extLIDVec[extLoc++];
  else
    li_Temperature = -1;

  if (dNodePrime == 1)
    li_DrainPrime = intLIDVec[intLoc++];
  else
    li_DrainPrime = li_Drain;

  if (sNodePrime == 1)
    li_SourcePrime = intLIDVec[intLoc++];
  else
    li_SourcePrime = li_Source;

  if (gNodePrime == 1)
    li_GatePrime = intLIDVec[intLoc++];
  else
    li_GatePrime = li_Gate;

  if (gNodeMid == 1)
    li_GateMid = intLIDVec[intLoc++];
  else
    li_GateMid = li_Gate;

  if( icVDSGiven )
  {
    if( li_Drain == li_Source )
    {
      msg = "Instance::registerLIDs:";
      msg += "Tried to specify an initial condition on V_Drain_Source ";
      msg += "when Drain and Source nodes are the same node.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    li_Ids = intLIDVec[intLoc++];
  }

  if( icVGSGiven )
  {
    if( li_Gate == li_Source )
    {
      msg = "Instance::registerLIDs:";
      msg += "Tried to specify an initial condition on V_Gate_Source ";
      msg += "when Gate and Source nodes are the same node.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    li_Igs = intLIDVec[intLoc++];
  }

  if( icVBSGiven )
  {
    if( (li_Body == li_Source) || (li_Body == -1) )
    {
      msg = "Instance::registerLIDs:";
      msg += "Tried to specify an initial condition on V_Body_Source ";
      msg += "when Body and Source nodes are the same node, or";
      msg += "when Body node does not exist.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    li_Ibs = intLIDVec[intLoc++];
  }

  if( icVESGiven )
  {
    if( li_Substrate == li_Source )
    {
      msg = "Instance::registerLIDs:";
      msg += "Tried to specify an initial condition on V_Substrate_Source ";
      msg += "when Substrate and Source nodes are the same node.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    li_Ies = intLIDVec[intLoc++];
  }
  if( icVPSGiven )
  {
    if( (li_ExtBody == li_Source) || (li_ExtBody == -1) )
    {
      msg = "Instance::registerLIDs:";
      msg += "Tried to specify an initial condition on V_ExtBody_Source ";
      msg += "when External Body and Source nodes are the same node, ";
      msg += "or when External Body node does not exist.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    li_Ips = intLIDVec[intLoc++];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "\n  local variable indices:\n";
    cout << "  li_Drain         = " << li_Drain << endl;
    cout << "  li_Gate          = " << li_Gate << endl;
    cout << "  li_Source        = " << li_Source << endl;
    cout << "  li_Substrate     = " << li_Substrate << endl;
    cout << "  li_ExtBody       = " << li_ExtBody << endl;
    cout << "  li_Body          = " << li_Body << endl;
    cout << "  li_Temperature   = " << li_Temperature << endl;
    cout << "  li_DrainPrime    = " << li_DrainPrime << endl;
    cout << "  li_SourcePrime   = " << li_SourcePrime << endl;
    cout << "  li_GatePrime     = " << li_GatePrime << endl;
    cout << "  li_GateMid       = " << li_GateMid << endl;
    if (icVDSGiven)
      cout << "  li_Ids         = " << li_Ids << endl;

    if (icVGSGiven)
      cout << "  li_Igs         = " << li_Igs << endl;

    if (icVBSGiven)
      cout << "  li_Ibs         = " << li_Ibs << endl;

    if (icVESGiven)
      cout << "  li_Ies         = " << li_Ies << endl;

    if (icVPSGiven)
      cout << "  li_Ips         = " << li_Ips << endl;

    cout << dashedline << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    // set up the internal names map
    string tmpstr;
    if (bNode == 1)
    {
      tmpstr = getName()+"_Body";
      spiceInternalName (tmpstr);
      intNameMap[li_Body] = tmpstr;
    }
    if (tNode == 1)
    {
      tmpstr = getName()+"_Temp";
      spiceInternalName (tmpstr);
      intNameMap[li_Temperature] = tmpstr;
    }
    if (dNodePrime == 1)
    {
      tmpstr = getName()+"_drain";
      spiceInternalName (tmpstr);
      intNameMap[li_DrainPrime] = tmpstr;
    }
    if (sNodePrime == 1)
    {
      tmpstr = getName()+"_source";
      spiceInternalName (tmpstr);
      intNameMap[li_SourcePrime] = tmpstr;
    }
    if (gNodePrime == 1)
    {
      tmpstr = getName()+"_gate";
      spiceInternalName (tmpstr);
      intNameMap[li_GatePrime] = tmpstr;
    }
    if (gNodeMid == 1)
    {
      tmpstr = getName()+"_midgate";
      spiceInternalName (tmpstr);
      intNameMap[li_GateMid] = tmpstr;
    }
    if (icVDSGiven)
    {
      tmpstr = getName()+"_branch_DS";
      spiceInternalName (tmpstr);
      intNameMap[li_Ids] = tmpstr;
    }
    if (icVGSGiven)
    {
      tmpstr = getName()+"_branch_GS";
      spiceInternalName (tmpstr);
      intNameMap[li_Igs] = tmpstr;
    }
    if (icVBSGiven)
    {
      tmpstr = getName()+"_branch_BS";
      spiceInternalName (tmpstr);
      intNameMap[li_Ibs] = tmpstr;
    }
    if (icVESGiven)
    {
      tmpstr = getName()+"_branch_ES";
      spiceInternalName (tmpstr);
      intNameMap[li_Ies] = tmpstr;
    }
  // Spice model documentation says to ignore icVPS if
  // there isn't a p terminal.
    if (pNode != 0)
    {
      if (icVPSGiven)
      {
        tmpstr = getName()+"_branch_PS";
        spiceInternalName (tmpstr);
        intNameMap[li_Ips] = tmpstr;
      }
    }

  }

  return intNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 4/4/2013
//-----------------------------------------------------------------------------
map<int,string> & Instance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if( loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    string modName(getName());
    spiceInternalName(modName);
    string tmpstr;
    tmpstr = modName+":DEV_ID";
    storeNameMap[ li_store_dev_id ] = tmpstr;
    tmpstr = modName+":DEV_IG";
    storeNameMap[ li_store_dev_ig ] = tmpstr;
    tmpstr = modName+":DEV_IS";
    storeNameMap[ li_store_dev_is ] = tmpstr;
    tmpstr = modName+":DEV_IE";
    storeNameMap[ li_store_dev_ie ] = tmpstr;
    tmpstr = modName+":DEV_IB";
    storeNameMap[ li_store_dev_ib ] = tmpstr;
  }

  return storeNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(
   const vector<int> & staLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (devOptions.debugLevel > 0)
  {
    cout << endl;
    cout << dashedline << endl;
    cout << "  In Instance::registerStateLIDs\n\n";
    cout << "  name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "  Number of State LIDs: " << numSta << endl;
  }
#endif

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int n = 0;
  // intrinsic capacitors
  li_state_qb       = staLIDVec[n++];
  li_state_qg       = staLIDVec[n++];
  li_state_qd       = staLIDVec[n++];
  li_state_qe       = staLIDVec[n++];
  li_state_qgmid    = staLIDVec[n++];
  li_state_qth      = staLIDVec[n++];

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "  Local State indices:" << endl;
    cout << endl;

    cout << "  li_state_qb           = " << li_state_qb << endl;
    cout << "  li_state_qg           = " << li_state_qg << endl;
    cout << "  li_state_qd           = " << li_state_qd << endl;
    cout << "  li_state_qe           = " << li_state_qe << endl;
    cout << "  li_state_qgmid        = " << li_state_qgmid << endl;
    cout << "  li_state_qth          = " << li_state_qth << endl;
    cout << endl;
    cout << dashedline << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(
   const vector<int> & stoLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (devOptions.debugLevel > 0)
  {
    cout << endl;
    cout << dashedline << endl;
    cout << "  In Instance::registerStoreLIDs\n\n";
    cout << "  name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSto = stoLIDVecRef.size();

  if (numSto != getNumStoreVars())
  {
    msg = "Instance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "  Number of Store LIDs: " << numSto << endl;
  }
#endif

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int n = 0;
    // voltage drops
  li_store_vbd      = stoLIDVec[n++];
  li_store_vbs      = stoLIDVec[n++];
  li_store_vgs      = stoLIDVec[n++];
  li_store_vds      = stoLIDVec[n++];
  li_store_ves      = stoLIDVec[n++];
  li_store_vps      = stoLIDVec[n++];

  li_store_vg       = stoLIDVec[n++];
  li_store_vd       = stoLIDVec[n++];
  li_store_vs       = stoLIDVec[n++];
  li_store_vp       = stoLIDVec[n++];
  li_store_ve       = stoLIDVec[n++];
  li_store_vgp     = stoLIDVec[n++];
  li_store_vgm     = stoLIDVec[n++];
  li_store_deltemp  = stoLIDVec[n++];

  li_store_vges     = stoLIDVec[n++];
  li_store_vgms     = stoLIDVec[n++];

  if( loadLeadCurrent )
  {
    li_store_dev_id = stoLIDVec[n++];
    li_store_dev_ig = stoLIDVec[n++];
    li_store_dev_is = stoLIDVec[n++];
    li_store_dev_ie = stoLIDVec[n++];
    li_store_dev_ib = stoLIDVec[n++];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "  Local Store indices:" << endl;
    cout << endl;
    cout << "  li_store_vbd          = " << li_store_vbd << endl;
    cout << "  li_store_vbs          = " << li_store_vbs << endl;
    cout << "  li_store_vgs          = " << li_store_vgs << endl;
    cout << "  li_store_vds          = " << li_store_vds << endl;
    cout << "  li_store_ves          = " << li_store_ves << endl;
    cout << "  li_store_vps          = " << li_store_vps << endl;
    cout << "  li_store_vg           = " << li_store_vg  << endl;
    cout << "  li_store_vd           = " << li_store_vd << endl;
    cout << "  li_store_vs           = " << li_store_vs << endl;
    cout << "  li_store_vp           = " << li_store_vp << endl;
    cout << "  li_store_ve           = " << li_store_ve << endl;
    cout << "  li_store_vgp          = " << li_store_vgp << endl;
    cout << "  li_store_vgm          = " << li_store_vgm << endl;
    cout << "  li_store_deltemp      = " << li_store_deltemp << endl;
    cout << "  li_store_vges         = " << li_store_vges << endl;
    cout << "  li_store_vgms         = " << li_store_vgms << endl;
    cout << endl;
    cout << dashedline << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  if ( icVDSGiven || icVGSGiven || icVBSGiven || icVESGiven || icVPSGiven )
  {
    return jacStampIC;
  }

  return jacStamp_v[jacID];
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs(
   const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  vector<int> map;
  vector< vector<int> > map2;

  // these are used to helpe us pluck the LID's from the right spot
  // without adding a lot of new code to each jacID condition below
  int lastDrainIndex     = 0;
  int lastGateIndex      = 0;
  int lastBodyIndex      = 0;
  int lastSourceIndex    = 0;
  int lastSubstrateIndex = 0;
  int lastExtBodyIndex   = 0;
  int lastRowUsed        = 0;

  if ( icVDSGiven || icVGSGiven || icVBSGiven || icVESGiven || icVPSGiven )
  {
    map = jacMapIC;
    map2 = jacMapIC2;
  }
  else
  {
    map = jacMap_v[jacID];
    map2 = jacMap2_v[jacID];
  }

  if (jacID < 24)
  {
    ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
    ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

    AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
    AGateEquBodyNodeOffset               = jacLIDVec[map[1]][map2[1][1]];
    AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][2]];
    AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][3]];
    AGateEquGatePrimeNodeOffset          = jacLIDVec[map[1]][map2[1][4]];
    AGateEquGateMidNodeOffset            = jacLIDVec[map[1]][map2[1][5]];

    ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
    ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

    ASubstrateEquSubstrateNodeOffset     = jacLIDVec[map[3]][map2[3][0]];
    ASubstrateEquBodyNodeOffset          = jacLIDVec[map[3]][map2[3][1]];
    ASubstrateEquTemperatureNodeOffset   = jacLIDVec[map[3]][map2[3][2]];
    ASubstrateEquDrainPrimeNodeOffset    = jacLIDVec[map[3]][map2[3][3]];
    ASubstrateEquSourcePrimeNodeOffset   = jacLIDVec[map[3]][map2[3][4]];
    ASubstrateEquGatePrimeNodeOffset     = jacLIDVec[map[3]][map2[3][5]];
    ASubstrateEquGateMidNodeOffset       = jacLIDVec[map[3]][map2[3][6]];

    AExtBodyEquExtBodyNodeOffset         = jacLIDVec[map[4]][map2[4][0]];
    AExtBodyEquBodyNodeOffset            = jacLIDVec[map[4]][map2[4][1]];

    ABodyEquSubstrateNodeOffset          = jacLIDVec[map[5]][map2[5][0]];
    ABodyEquExtBodyNodeOffset            = jacLIDVec[map[5]][map2[5][1]];
    ABodyEquBodyNodeOffset               = jacLIDVec[map[5]][map2[5][2]];
    ABodyEquTemperatureNodeOffset        = jacLIDVec[map[5]][map2[5][3]];
    ABodyEquDrainPrimeNodeOffset         = jacLIDVec[map[5]][map2[5][4]];
    ABodyEquSourcePrimeNodeOffset        = jacLIDVec[map[5]][map2[5][5]];
    ABodyEquGatePrimeNodeOffset          = jacLIDVec[map[5]][map2[5][6]];

    ATemperatureEquSubstrateNodeOffset   = jacLIDVec[map[6]][map2[6][0]];
    ATemperatureEquBodyNodeOffset        = jacLIDVec[map[6]][map2[6][1]];
    ATemperatureEquTemperatureNodeOffset = jacLIDVec[map[6]][map2[6][2]];
    ATemperatureEquDrainPrimeNodeOffset  = jacLIDVec[map[6]][map2[6][3]];
    ATemperatureEquSourcePrimeNodeOffset = jacLIDVec[map[6]][map2[6][4]];
    ATemperatureEquGatePrimeNodeOffset   = jacLIDVec[map[6]][map2[6][5]];

    ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[7]][map2[7][0]];
    ADrainPrimeEquSubstrateNodeOffset    = jacLIDVec[map[7]][map2[7][1]];
    ADrainPrimeEquBodyNodeOffset         = jacLIDVec[map[7]][map2[7][2]];
    ADrainPrimeEquTemperatureNodeOffset  = jacLIDVec[map[7]][map2[7][3]];
    ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[7]][map2[7][4]];
    ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[7]][map2[7][5]];
    ADrainPrimeEquGatePrimeNodeOffset    = jacLIDVec[map[7]][map2[7][6]];
    ADrainPrimeEquGateMidNodeOffset      = jacLIDVec[map[7]][map2[7][7]];

    ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[8]][map2[8][0]];
    ASourcePrimeEquSubstrateNodeOffset   = jacLIDVec[map[8]][map2[8][1]];
    ASourcePrimeEquBodyNodeOffset        = jacLIDVec[map[8]][map2[8][2]];
    ASourcePrimeEquTemperatureNodeOffset = jacLIDVec[map[8]][map2[8][3]];
    ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[8]][map2[8][4]];
    ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[8]][map2[8][5]];
    ASourcePrimeEquGatePrimeNodeOffset   = jacLIDVec[map[8]][map2[8][6]];
    ASourcePrimeEquGateMidNodeOffset     = jacLIDVec[map[8]][map2[8][7]];

    AGatePrimeEquGateNodeOffset          = jacLIDVec[map[9]][map2[9][0]];
    AGatePrimeEquSubstrateNodeOffset     = jacLIDVec[map[9]][map2[9][1]];
    AGatePrimeEquBodyNodeOffset          = jacLIDVec[map[9]][map2[9][2]];
    AGatePrimeEquTemperatureNodeOffset   = jacLIDVec[map[9]][map2[9][3]];
    AGatePrimeEquDrainPrimeNodeOffset    = jacLIDVec[map[9]][map2[9][4]];
    AGatePrimeEquSourcePrimeNodeOffset   = jacLIDVec[map[9]][map2[9][5]];
    AGatePrimeEquGatePrimeNodeOffset     = jacLIDVec[map[9]][map2[9][6]];
    AGatePrimeEquGateMidNodeOffset       = jacLIDVec[map[9]][map2[9][7]];

    AGateMidEquGateNodeOffset            = jacLIDVec[map[10]][map2[10][0]];
    AGateMidEquSubstrateNodeOffset       = jacLIDVec[map[10]][map2[10][1]];
    AGateMidEquBodyNodeOffset            = jacLIDVec[map[10]][map2[10][2]];
    AGateMidEquDrainPrimeNodeOffset      = jacLIDVec[map[10]][map2[10][3]];
    AGateMidEquSourcePrimeNodeOffset     = jacLIDVec[map[10]][map2[10][4]];
    AGateMidEquGatePrimeNodeOffset       = jacLIDVec[map[10]][map2[10][5]];
    AGateMidEquGateMidNodeOffset         = jacLIDVec[map[10]][map2[10][6]];

    lastDrainIndex     =  1;
    lastGateIndex      =  5;
    lastBodyIndex      =  6;
    lastSourceIndex    =  1;
    lastSubstrateIndex =  6;
    lastExtBodyIndex   =  1;
    lastRowUsed        = 10;
  }
  else if (jacID >= 24 && jacID < 48)
  {
    ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
    ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

    AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
    AGateEquBodyNodeOffset               = jacLIDVec[map[1]][map2[1][1]];
    AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][2]];
    AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][3]];
    AGateEquGatePrimeNodeOffset          = jacLIDVec[map[1]][map2[1][4]];
    AGateEquGateMidNodeOffset            = jacLIDVec[map[1]][map2[1][5]];

    ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
    ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

    ASubstrateEquSubstrateNodeOffset     = jacLIDVec[map[3]][map2[3][0]];
    ASubstrateEquBodyNodeOffset          = jacLIDVec[map[3]][map2[3][1]];
    ASubstrateEquDrainPrimeNodeOffset    = jacLIDVec[map[3]][map2[3][2]];
    ASubstrateEquSourcePrimeNodeOffset   = jacLIDVec[map[3]][map2[3][3]];
    ASubstrateEquGatePrimeNodeOffset     = jacLIDVec[map[3]][map2[3][4]];
    ASubstrateEquGateMidNodeOffset       = jacLIDVec[map[3]][map2[3][5]];

    AExtBodyEquExtBodyNodeOffset         = jacLIDVec[map[4]][map2[4][0]];
    AExtBodyEquBodyNodeOffset            = jacLIDVec[map[4]][map2[4][1]];

    ABodyEquSubstrateNodeOffset          = jacLIDVec[map[5]][map2[5][0]];
    ABodyEquExtBodyNodeOffset            = jacLIDVec[map[5]][map2[5][1]];
    ABodyEquBodyNodeOffset               = jacLIDVec[map[5]][map2[5][2]];
    ABodyEquDrainPrimeNodeOffset         = jacLIDVec[map[5]][map2[5][3]];
    ABodyEquSourcePrimeNodeOffset        = jacLIDVec[map[5]][map2[5][4]];
    ABodyEquGatePrimeNodeOffset          = jacLIDVec[map[5]][map2[5][5]];

    ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[6]][map2[6][0]];
    ADrainPrimeEquSubstrateNodeOffset    = jacLIDVec[map[6]][map2[6][1]];
    ADrainPrimeEquBodyNodeOffset         = jacLIDVec[map[6]][map2[6][2]];
    ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[6]][map2[6][3]];
    ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[6]][map2[6][4]];
    ADrainPrimeEquGatePrimeNodeOffset    = jacLIDVec[map[6]][map2[6][5]];
    ADrainPrimeEquGateMidNodeOffset      = jacLIDVec[map[6]][map2[6][6]];

    ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[7]][map2[7][0]];
    ASourcePrimeEquSubstrateNodeOffset   = jacLIDVec[map[7]][map2[7][1]];
    ASourcePrimeEquBodyNodeOffset        = jacLIDVec[map[7]][map2[7][2]];
    ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[7]][map2[7][3]];
    ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[7]][map2[7][4]];
    ASourcePrimeEquGatePrimeNodeOffset   = jacLIDVec[map[7]][map2[7][5]];
    ASourcePrimeEquGateMidNodeOffset     = jacLIDVec[map[7]][map2[7][6]];

    AGatePrimeEquGateNodeOffset          = jacLIDVec[map[8]][map2[8][0]];
    AGatePrimeEquSubstrateNodeOffset     = jacLIDVec[map[8]][map2[8][1]];
    AGatePrimeEquBodyNodeOffset          = jacLIDVec[map[8]][map2[8][2]];
    AGatePrimeEquDrainPrimeNodeOffset    = jacLIDVec[map[8]][map2[8][3]];
    AGatePrimeEquSourcePrimeNodeOffset   = jacLIDVec[map[8]][map2[8][4]];
    AGatePrimeEquGatePrimeNodeOffset     = jacLIDVec[map[8]][map2[8][5]];
    AGatePrimeEquGateMidNodeOffset       = jacLIDVec[map[8]][map2[8][6]];

    AGateMidEquGateNodeOffset            = jacLIDVec[map[9]][map2[9][0]];
    AGateMidEquSubstrateNodeOffset       = jacLIDVec[map[9]][map2[9][1]];
    AGateMidEquBodyNodeOffset            = jacLIDVec[map[9]][map2[9][2]];
    AGateMidEquDrainPrimeNodeOffset      = jacLIDVec[map[9]][map2[9][3]];
    AGateMidEquSourcePrimeNodeOffset     = jacLIDVec[map[9]][map2[9][4]];
    AGateMidEquGatePrimeNodeOffset       = jacLIDVec[map[9]][map2[9][5]];
    AGateMidEquGateMidNodeOffset         = jacLIDVec[map[9]][map2[9][6]];

    lastDrainIndex     =  1;
    lastGateIndex      =  5;
    lastBodyIndex      =  5;
    lastSourceIndex    =  1;
    lastSubstrateIndex =  5;
    lastExtBodyIndex   =  1;
    lastRowUsed        =  9;
  }
  else if (jacID >= 48 && jacID < 60)
  {
    ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
    ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

    AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
    AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][1]];
    AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][2]];
    AGateEquGatePrimeNodeOffset          = jacLIDVec[map[1]][map2[1][3]];
    AGateEquGateMidNodeOffset            = jacLIDVec[map[1]][map2[1][4]];

    ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
    ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

    ASubstrateEquSubstrateNodeOffset     = jacLIDVec[map[3]][map2[3][0]];
    ASubstrateEquTemperatureNodeOffset   = jacLIDVec[map[3]][map2[3][1]];
    ASubstrateEquDrainPrimeNodeOffset    = jacLIDVec[map[3]][map2[3][2]];
    ASubstrateEquSourcePrimeNodeOffset   = jacLIDVec[map[3]][map2[3][3]];
    ASubstrateEquGatePrimeNodeOffset     = jacLIDVec[map[3]][map2[3][4]];
    ASubstrateEquGateMidNodeOffset       = jacLIDVec[map[3]][map2[3][5]];

    ATemperatureEquSubstrateNodeOffset   = jacLIDVec[map[4]][map2[4][0]];
    ATemperatureEquTemperatureNodeOffset = jacLIDVec[map[4]][map2[4][1]];
    ATemperatureEquDrainPrimeNodeOffset  = jacLIDVec[map[4]][map2[4][2]];
    ATemperatureEquSourcePrimeNodeOffset = jacLIDVec[map[4]][map2[4][3]];
    ATemperatureEquGatePrimeNodeOffset   = jacLIDVec[map[4]][map2[4][4]];

    ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[5]][map2[5][0]];
    ADrainPrimeEquSubstrateNodeOffset    = jacLIDVec[map[5]][map2[5][1]];
    ADrainPrimeEquTemperatureNodeOffset  = jacLIDVec[map[5]][map2[5][2]];
    ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[5]][map2[5][3]];
    ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[5]][map2[5][4]];
    ADrainPrimeEquGatePrimeNodeOffset    = jacLIDVec[map[5]][map2[5][5]];
    ADrainPrimeEquGateMidNodeOffset      = jacLIDVec[map[5]][map2[5][6]];

    ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[6]][map2[6][0]];
    ASourcePrimeEquSubstrateNodeOffset   = jacLIDVec[map[6]][map2[6][1]];
    ASourcePrimeEquTemperatureNodeOffset = jacLIDVec[map[6]][map2[6][2]];
    ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[6]][map2[6][3]];
    ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[6]][map2[6][4]];
    ASourcePrimeEquGatePrimeNodeOffset   = jacLIDVec[map[6]][map2[6][5]];
    ASourcePrimeEquGateMidNodeOffset     = jacLIDVec[map[6]][map2[6][6]];

    AGatePrimeEquGateNodeOffset          = jacLIDVec[map[7]][map2[7][0]];
    AGatePrimeEquSubstrateNodeOffset     = jacLIDVec[map[7]][map2[7][1]];
    AGatePrimeEquTemperatureNodeOffset   = jacLIDVec[map[7]][map2[7][2]];
    AGatePrimeEquDrainPrimeNodeOffset    = jacLIDVec[map[7]][map2[7][3]];
    AGatePrimeEquSourcePrimeNodeOffset   = jacLIDVec[map[7]][map2[7][4]];
    AGatePrimeEquGatePrimeNodeOffset     = jacLIDVec[map[7]][map2[7][5]];
    AGatePrimeEquGateMidNodeOffset       = jacLIDVec[map[7]][map2[7][6]];

    AGateMidEquGateNodeOffset            = jacLIDVec[map[8]][map2[8][0]];
    AGateMidEquSubstrateNodeOffset       = jacLIDVec[map[8]][map2[8][1]];
    AGateMidEquDrainPrimeNodeOffset      = jacLIDVec[map[8]][map2[8][2]];
    AGateMidEquSourcePrimeNodeOffset     = jacLIDVec[map[8]][map2[8][3]];
    AGateMidEquGatePrimeNodeOffset       = jacLIDVec[map[8]][map2[8][4]];
    AGateMidEquGateMidNodeOffset         = jacLIDVec[map[8]][map2[8][5]];

    lastDrainIndex     =  1;
    lastGateIndex      =  4;
    lastBodyIndex      =  0;
    lastSourceIndex    =  1;
    lastSubstrateIndex =  5;
    lastExtBodyIndex   =  0;
    lastRowUsed        =  8;
  }
  else if (jacID >= 60 && jacID < 72)
  {
    ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
    ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

    AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
    AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][1]];
    AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][2]];
    AGateEquGatePrimeNodeOffset          = jacLIDVec[map[1]][map2[1][3]];
    AGateEquGateMidNodeOffset            = jacLIDVec[map[1]][map2[1][4]];

    ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
    ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

    ASubstrateEquSubstrateNodeOffset     = jacLIDVec[map[3]][map2[3][0]];
    ASubstrateEquDrainPrimeNodeOffset    = jacLIDVec[map[3]][map2[3][1]];
    ASubstrateEquSourcePrimeNodeOffset   = jacLIDVec[map[3]][map2[3][2]];
    ASubstrateEquGatePrimeNodeOffset     = jacLIDVec[map[3]][map2[3][3]];
    ASubstrateEquGateMidNodeOffset       = jacLIDVec[map[3]][map2[3][4]];

    ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[4]][map2[4][0]];
    ADrainPrimeEquSubstrateNodeOffset    = jacLIDVec[map[4]][map2[4][1]];
    ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[4]][map2[4][2]];
    ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[4]][map2[4][3]];
    ADrainPrimeEquGatePrimeNodeOffset    = jacLIDVec[map[4]][map2[4][4]];
    ADrainPrimeEquGateMidNodeOffset      = jacLIDVec[map[4]][map2[4][5]];

    ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[5]][map2[5][0]];
    ASourcePrimeEquSubstrateNodeOffset   = jacLIDVec[map[5]][map2[5][1]];
    ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[5]][map2[5][2]];
    ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[5]][map2[5][3]];
    ASourcePrimeEquGatePrimeNodeOffset   = jacLIDVec[map[5]][map2[5][4]];
    ASourcePrimeEquGateMidNodeOffset     = jacLIDVec[map[5]][map2[5][5]];

    AGatePrimeEquGateNodeOffset          = jacLIDVec[map[6]][map2[6][0]];
    AGatePrimeEquSubstrateNodeOffset     = jacLIDVec[map[6]][map2[6][1]];
    AGatePrimeEquDrainPrimeNodeOffset    = jacLIDVec[map[6]][map2[6][2]];
    AGatePrimeEquSourcePrimeNodeOffset   = jacLIDVec[map[6]][map2[6][3]];
    AGatePrimeEquGatePrimeNodeOffset     = jacLIDVec[map[6]][map2[6][4]];
    AGatePrimeEquGateMidNodeOffset       = jacLIDVec[map[6]][map2[6][5]];

    AGateMidEquGateNodeOffset            = jacLIDVec[map[7]][map2[7][0]];
    AGateMidEquSubstrateNodeOffset       = jacLIDVec[map[7]][map2[7][1]];
    AGateMidEquDrainPrimeNodeOffset      = jacLIDVec[map[7]][map2[7][2]];
    AGateMidEquSourcePrimeNodeOffset     = jacLIDVec[map[7]][map2[7][3]];
    AGateMidEquGatePrimeNodeOffset       = jacLIDVec[map[7]][map2[7][4]];
    AGateMidEquGateMidNodeOffset         = jacLIDVec[map[7]][map2[7][5]];

    lastDrainIndex     =  1;
    lastGateIndex      =  4;
    lastBodyIndex      =  0;
    lastSourceIndex    =  1;
    lastSubstrateIndex =  4;
    lastExtBodyIndex   =  0;
    lastRowUsed        =  7;
  }
  else if (jacID >= 72 && jacID < 84)
  {
    ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
    ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

    AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
    AGateEquBodyNodeOffset               = jacLIDVec[map[1]][map2[1][1]];
    AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][2]];
    AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][3]];
    AGateEquGatePrimeNodeOffset          = jacLIDVec[map[1]][map2[1][4]];
    AGateEquGateMidNodeOffset            = jacLIDVec[map[1]][map2[1][5]];

    ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
    ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

    ASubstrateEquSubstrateNodeOffset     = jacLIDVec[map[3]][map2[3][0]];
    ASubstrateEquTemperatureNodeOffset   = jacLIDVec[map[3]][map2[3][1]];
    ASubstrateEquBodyNodeOffset          = jacLIDVec[map[3]][map2[3][2]];
    ASubstrateEquDrainPrimeNodeOffset    = jacLIDVec[map[3]][map2[3][3]];
    ASubstrateEquSourcePrimeNodeOffset   = jacLIDVec[map[3]][map2[3][4]];
    ASubstrateEquGatePrimeNodeOffset     = jacLIDVec[map[3]][map2[3][5]];
    ASubstrateEquGateMidNodeOffset       = jacLIDVec[map[3]][map2[3][6]];

    AExtBodyEquExtBodyNodeOffset         = jacLIDVec[map[4]][map2[4][0]];
    AExtBodyEquBodyNodeOffset            = jacLIDVec[map[4]][map2[4][1]];

    ATemperatureEquSubstrateNodeOffset   = jacLIDVec[map[5]][map2[5][0]];
    ATemperatureEquTemperatureNodeOffset = jacLIDVec[map[5]][map2[5][1]];
    ATemperatureEquBodyNodeOffset        = jacLIDVec[map[5]][map2[5][2]];
    ATemperatureEquDrainPrimeNodeOffset  = jacLIDVec[map[5]][map2[5][3]];
    ATemperatureEquSourcePrimeNodeOffset = jacLIDVec[map[5]][map2[5][4]];
    ATemperatureEquGatePrimeNodeOffset   = jacLIDVec[map[5]][map2[5][5]];

    ABodyEquSubstrateNodeOffset          = jacLIDVec[map[6]][map2[6][0]];
    ABodyEquExtBodyNodeOffset            = jacLIDVec[map[6]][map2[6][1]];
    ABodyEquTemperatureNodeOffset        = jacLIDVec[map[6]][map2[6][2]];
    ABodyEquBodyNodeOffset               = jacLIDVec[map[6]][map2[6][3]];
    ABodyEquDrainPrimeNodeOffset         = jacLIDVec[map[6]][map2[6][4]];
    ABodyEquSourcePrimeNodeOffset        = jacLIDVec[map[6]][map2[6][5]];
    ABodyEquGatePrimeNodeOffset          = jacLIDVec[map[6]][map2[6][6]];

    ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[7]][map2[7][0]];
    ADrainPrimeEquSubstrateNodeOffset    = jacLIDVec[map[7]][map2[7][1]];
    ADrainPrimeEquTemperatureNodeOffset  = jacLIDVec[map[7]][map2[7][2]];
    ADrainPrimeEquBodyNodeOffset         = jacLIDVec[map[7]][map2[7][3]];
    ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[7]][map2[7][4]];
    ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[7]][map2[7][5]];
    ADrainPrimeEquGatePrimeNodeOffset    = jacLIDVec[map[7]][map2[7][6]];
    ADrainPrimeEquGateMidNodeOffset      = jacLIDVec[map[7]][map2[7][7]];

    ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[8]][map2[8][0]];
    ASourcePrimeEquSubstrateNodeOffset   = jacLIDVec[map[8]][map2[8][1]];
    ASourcePrimeEquTemperatureNodeOffset = jacLIDVec[map[8]][map2[8][2]];
    ASourcePrimeEquBodyNodeOffset        = jacLIDVec[map[8]][map2[8][3]];
    ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[8]][map2[8][4]];
    ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[8]][map2[8][5]];
    ASourcePrimeEquGatePrimeNodeOffset   = jacLIDVec[map[8]][map2[8][6]];
    ASourcePrimeEquGateMidNodeOffset     = jacLIDVec[map[8]][map2[8][7]];

    AGatePrimeEquGateNodeOffset          = jacLIDVec[map[9]][map2[9][0]];
    AGatePrimeEquSubstrateNodeOffset     = jacLIDVec[map[9]][map2[9][1]];
    AGatePrimeEquTemperatureNodeOffset   = jacLIDVec[map[9]][map2[9][2]];
    AGatePrimeEquBodyNodeOffset          = jacLIDVec[map[9]][map2[9][3]];
    AGatePrimeEquDrainPrimeNodeOffset    = jacLIDVec[map[9]][map2[9][4]];
    AGatePrimeEquSourcePrimeNodeOffset   = jacLIDVec[map[9]][map2[9][5]];
    AGatePrimeEquGatePrimeNodeOffset     = jacLIDVec[map[9]][map2[9][6]];
    AGatePrimeEquGateMidNodeOffset       = jacLIDVec[map[9]][map2[9][7]];

    AGateMidEquGateNodeOffset            = jacLIDVec[map[10]][map2[10][0]];
    AGateMidEquSubstrateNodeOffset       = jacLIDVec[map[10]][map2[10][1]];
    AGateMidEquBodyNodeOffset            = jacLIDVec[map[10]][map2[10][2]];
    AGateMidEquDrainPrimeNodeOffset      = jacLIDVec[map[10]][map2[10][3]];
    AGateMidEquSourcePrimeNodeOffset     = jacLIDVec[map[10]][map2[10][4]];
    AGateMidEquGatePrimeNodeOffset       = jacLIDVec[map[10]][map2[10][5]];
    AGateMidEquGateMidNodeOffset         = jacLIDVec[map[10]][map2[10][6]];

    lastDrainIndex     =  1;
    lastGateIndex      =  5;
    lastBodyIndex      =  6;
    lastSourceIndex    =  1;
    lastSubstrateIndex =  6;
    lastExtBodyIndex   =  1;
    lastRowUsed        = 10;
  }
  else if (jacID >= 84 && jacID < 96)
  {
    ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
    ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

    AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
    AGateEquBodyNodeOffset               = jacLIDVec[map[1]][map2[1][1]];
    AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][2]];
    AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][3]];
    AGateEquGatePrimeNodeOffset          = jacLIDVec[map[1]][map2[1][4]];
    AGateEquGateMidNodeOffset            = jacLIDVec[map[1]][map2[1][5]];

    ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
    ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

    ASubstrateEquSubstrateNodeOffset     = jacLIDVec[map[3]][map2[3][0]];
    ASubstrateEquTemperatureNodeOffset   = jacLIDVec[map[3]][map2[3][1]];
    ASubstrateEquBodyNodeOffset          = jacLIDVec[map[3]][map2[3][2]];
    ASubstrateEquDrainPrimeNodeOffset    = jacLIDVec[map[3]][map2[3][3]];
    ASubstrateEquSourcePrimeNodeOffset   = jacLIDVec[map[3]][map2[3][4]];
    ASubstrateEquGatePrimeNodeOffset     = jacLIDVec[map[3]][map2[3][5]];
    ASubstrateEquGateMidNodeOffset       = jacLIDVec[map[3]][map2[3][6]];

    ATemperatureEquSubstrateNodeOffset   = jacLIDVec[map[4]][map2[4][0]];
    ATemperatureEquTemperatureNodeOffset = jacLIDVec[map[4]][map2[4][1]];
    ATemperatureEquBodyNodeOffset        = jacLIDVec[map[4]][map2[4][2]];
    ATemperatureEquDrainPrimeNodeOffset  = jacLIDVec[map[4]][map2[4][3]];
    ATemperatureEquSourcePrimeNodeOffset = jacLIDVec[map[4]][map2[4][4]];
    ATemperatureEquGatePrimeNodeOffset   = jacLIDVec[map[4]][map2[4][5]];

    ABodyEquSubstrateNodeOffset          = jacLIDVec[map[5]][map2[5][0]];
    ABodyEquTemperatureNodeOffset        = jacLIDVec[map[5]][map2[5][1]];
    ABodyEquBodyNodeOffset               = jacLIDVec[map[5]][map2[5][2]];

    ABodyEquExtBodyNodeOffset            = jacLIDVec[map[5]][map2[5][3]];
    ABodyEquDrainPrimeNodeOffset         = jacLIDVec[map[5]][map2[5][4]];
    ABodyEquSourcePrimeNodeOffset        = jacLIDVec[map[5]][map2[5][5]];
    ABodyEquGatePrimeNodeOffset          = jacLIDVec[map[5]][map2[5][6]];

    AExtBodyEquBodyNodeOffset            = jacLIDVec[map[6]][map2[6][0]];
    AExtBodyEquExtBodyNodeOffset         = jacLIDVec[map[6]][map2[6][1]];

    ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[7]][map2[7][0]];
    ADrainPrimeEquSubstrateNodeOffset    = jacLIDVec[map[7]][map2[7][1]];
    ADrainPrimeEquTemperatureNodeOffset  = jacLIDVec[map[7]][map2[7][2]];
    ADrainPrimeEquBodyNodeOffset         = jacLIDVec[map[7]][map2[7][3]];
    ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[7]][map2[7][4]];
    ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[7]][map2[7][5]];
    ADrainPrimeEquGatePrimeNodeOffset    = jacLIDVec[map[7]][map2[7][6]];
    ADrainPrimeEquGateMidNodeOffset      = jacLIDVec[map[7]][map2[7][7]];

    ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[8]][map2[8][0]];
    ASourcePrimeEquSubstrateNodeOffset   = jacLIDVec[map[8]][map2[8][1]];
    ASourcePrimeEquTemperatureNodeOffset = jacLIDVec[map[8]][map2[8][2]];
    ASourcePrimeEquBodyNodeOffset        = jacLIDVec[map[8]][map2[8][3]];
    ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[8]][map2[8][4]];
    ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[8]][map2[8][5]];
    ASourcePrimeEquGatePrimeNodeOffset   = jacLIDVec[map[8]][map2[8][6]];
    ASourcePrimeEquGateMidNodeOffset     = jacLIDVec[map[8]][map2[8][7]];

    AGatePrimeEquGateNodeOffset          = jacLIDVec[map[9]][map2[9][0]];
    AGatePrimeEquSubstrateNodeOffset     = jacLIDVec[map[9]][map2[9][1]];
    AGatePrimeEquTemperatureNodeOffset   = jacLIDVec[map[9]][map2[9][2]];
    AGatePrimeEquBodyNodeOffset          = jacLIDVec[map[9]][map2[9][3]];
    AGatePrimeEquDrainPrimeNodeOffset    = jacLIDVec[map[9]][map2[9][4]];
    AGatePrimeEquSourcePrimeNodeOffset   = jacLIDVec[map[9]][map2[9][5]];
    AGatePrimeEquGatePrimeNodeOffset     = jacLIDVec[map[9]][map2[9][6]];
    AGatePrimeEquGateMidNodeOffset       = jacLIDVec[map[9]][map2[9][7]];

    AGateMidEquGateNodeOffset            = jacLIDVec[map[10]][map2[10][0]];
    AGateMidEquSubstrateNodeOffset       = jacLIDVec[map[10]][map2[10][1]];
    AGateMidEquBodyNodeOffset            = jacLIDVec[map[10]][map2[10][2]];
    AGateMidEquDrainPrimeNodeOffset      = jacLIDVec[map[10]][map2[10][3]];
    AGateMidEquSourcePrimeNodeOffset     = jacLIDVec[map[10]][map2[10][4]];
    AGateMidEquGatePrimeNodeOffset       = jacLIDVec[map[10]][map2[10][5]];
    AGateMidEquGateMidNodeOffset         = jacLIDVec[map[10]][map2[10][6]];

    lastDrainIndex     =  1;
    lastGateIndex      =  5;
    lastBodyIndex      =  6;
    lastSourceIndex    =  1;
    lastSubstrateIndex =  6;
    lastExtBodyIndex   =  1;
    lastRowUsed        = 10;
  }
  else
  {
    string msg = "Instance::registerJacLIDs:";
    msg += "jacID out of supported range";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  int extraVars = 0;
  if( icVDSGiven )
  {
    ++extraVars;
    ADrainEquIdsOffset     = jacLIDVec[map[0]][map2[0][lastDrainIndex+1]];
    ASourceEquIdsOffset    = jacLIDVec[map[2]][map2[2][lastSourceIndex+extraVars]];;
    icVDSEquVdOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][0]];
    icVDSEquVsOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][1]];
    icVDSEquIdsOffset      = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][2]];
  }
  if( icVGSGiven )
  {
    ++extraVars;
    AGateEquIgsOffset      = jacLIDVec[map[1]][map2[1][lastGateIndex+1]];
    ASourceEquIgsOffset    = jacLIDVec[map[2]][map2[2][lastSourceIndex+extraVars]];
    icVGSEquVgOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][0]];
    icVGSEquVsOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][1]];
    icVGSEquIgsOffset      = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][2]];
  }
  if( icVBSGiven )
  {
    ++extraVars;
    ABodyEquIbsOffset      = jacLIDVec[map[5]][map2[5][lastBodyIndex+1]];
    ASourceEquIbsOffset    = jacLIDVec[map[2]][map2[2][lastSourceIndex+extraVars]];
    icVBSEquVbOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][0]];
    icVBSEquVsOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][1]];
    icVBSEquIbsOffset      = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][2]];
  }
  if( icVESGiven )
  {
    ++extraVars;
    ASubstrateEquIesOffset = jacLIDVec[map[3]][map2[3][lastSubstrateIndex+1]];
    ASourceEquIesOffset    = jacLIDVec[map[2]][map2[2][lastSourceIndex+extraVars]];
    icVESEquVeOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][0]];
    icVESEquVsOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][1]];
    icVESEquIesOffset      = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][2]];
  }
  if( icVPSGiven )
  {
    ++extraVars;
    AExtBodyEquIpsOffset   = jacLIDVec[map[4]][map2[4][lastExtBodyIndex+1]];
    ASourceEquIpsOffset    = jacLIDVec[map[2]][map2[2][lastSourceIndex+extraVars]];
    icVPSEquVpOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][0]];
    icVPSEquVsOffset       = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][1]];
    icVPSEquIpsOffset      = jacLIDVec[map[lastRowUsed+extraVars]][map2[lastRowUsed+extraVars][2]];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/07/09
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  // F-matrix:
  if (rgateMod == 1)
  {
    f_GateEquGateNodePtr = &(dFdx[li_Gate][AGateEquGateNodeOffset]);
    f_GatePrimeEquGateNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquGateNodeOffset]);
    f_GateEquGatePrimeNodePtr = &(dFdx[li_Gate][AGateEquGatePrimeNodeOffset]);
    // It seems that this should be here, but it shouldn't.
    // The mi.geltd term is added into the G'-G' element later down in this
    // routine
    //
    //  f_GatePrime][AGatePrimeEquGatePrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]);
  }
  else if (rgateMod == 2)
  {
    f_GateEquGateNodePtr = &(dFdx[li_Gate][AGateEquGateNodeOffset]);
    f_GateEquGatePrimeNodePtr = &(dFdx[li_Gate][AGateEquGatePrimeNodeOffset]);
    f_GateEquDrainPrimeNodePtr = &(dFdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
    f_GateEquSourcePrimeNodePtr = &(dFdx[li_Gate][AGateEquSourcePrimeNodeOffset]);
    f_GatePrimeEquGateNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquGateNodeOffset]);
    if (soiMod != 2)
      f_GateEquBodyNodePtr = &(dFdx[li_Gate][AGateEquBodyNodeOffset]);
  }
  else if (rgateMod == 3)
  {
    f_GateEquGateNodePtr = &(dFdx[li_Gate][AGateEquGateNodeOffset]);
    f_GateEquGateMidNodePtr = &(dFdx[li_Gate][AGateEquGateMidNodeOffset]);
    f_GateMidEquGateNodePtr = &(dFdx[li_GateMid][AGateMidEquGateNodeOffset]);
    f_GateMidEquGateMidNodePtr = &(dFdx[li_GateMid][AGateMidEquGateMidNodeOffset]);
    f_GateMidEquDrainPrimeNodePtr = &(dFdx[li_GateMid][AGateMidEquDrainPrimeNodeOffset]);
    f_GateMidEquGatePrimeNodePtr = &(dFdx[li_GateMid][AGateMidEquGatePrimeNodeOffset]);
    f_GateMidEquSourcePrimeNodePtr = &(dFdx[li_GateMid][AGateMidEquSourcePrimeNodeOffset]);

    if (soiMod != 2)
    {
      f_GateMidEquBodyNodePtr = &(dFdx[li_GateMid][AGateMidEquBodyNodeOffset]);
    }

    f_GatePrimeEquGateMidNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquGateMidNodeOffset]);
  }
  if (soiMod != 0)
  {
    f_DrainPrimeEquSubstrateNodePtr = &(dFdx[li_DrainPrime][ADrainPrimeEquSubstrateNodeOffset]);
    f_SourcePrimeEquSubstrateNodePtr = &(dFdx[li_SourcePrime][ASourcePrimeEquSubstrateNodeOffset]);
    if (soiMod != 2)
    {
      f_GatePrimeEquSubstrateNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquSubstrateNodeOffset]);
      f_BodyEquSubstrateNodePtr = &(dFdx[li_Body][ABodyEquSubstrateNodeOffset]);
    }
  }

  if (soiMod != 2)
  {
    if (rgateMod == 0 || rgateMod == 1)
    {
      f_GatePrimeEquBodyNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquBodyNodeOffset]);
    }
    else
    {
      f_GatePrimeEquBodyNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquBodyNodeOffset]);
    }
    f_DrainPrimeEquBodyNodePtr = &(dFdx[li_DrainPrime][ADrainPrimeEquBodyNodeOffset]);
    f_SourcePrimeEquBodyNodePtr = &(dFdx[li_SourcePrime][ASourcePrimeEquBodyNodeOffset]);
    f_BodyEquSubstrateNodePtr = &(dFdx[li_Body][ABodyEquSubstrateNodeOffset]);
    f_BodyEquGatePrimeNodePtr = &(dFdx[li_Body][ABodyEquGatePrimeNodeOffset]);
    f_BodyEquDrainPrimeNodePtr = &(dFdx[li_Body][ABodyEquDrainPrimeNodeOffset]);
    f_BodyEquSourcePrimeNodePtr = &(dFdx[li_Body][ABodyEquSourcePrimeNodeOffset]);
    f_BodyEquBodyNodePtr = &(dFdx[li_Body][ABodyEquBodyNodeOffset]);
  }
  if (rgateMod == 0)
  {
    f_GatePrimeEquGatePrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]);
    f_GatePrimeEquDrainPrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]);
    f_GatePrimeEquSourcePrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]);
  }
  else if (rgateMod == 1)
  {
    f_GatePrimeEquGatePrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]);
    f_GatePrimeEquDrainPrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]);
    f_GatePrimeEquSourcePrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]);
  }
  else
  {
    f_GatePrimeEquGatePrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]);
    f_GatePrimeEquDrainPrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]);
    f_GatePrimeEquSourcePrimeNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]);
  }
  f_DrainPrimeEquGatePrimeNodePtr = &(dFdx[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset]);
  f_DrainPrimeEquDrainPrimeNodePtr = &(dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  f_DrainPrimeEquSourcePrimeNodePtr = &(dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);
  f_DrainPrimeEquDrainNodePtr = &(dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  f_SourcePrimeEquGatePrimeNodePtr = &(dFdx[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset]);
  f_SourcePrimeEquDrainPrimeNodePtr = &(dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  f_SourcePrimeEquSourcePrimeNodePtr = &(dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);
  f_SourcePrimeEquSourceNodePtr = &(dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  f_DrainEquDrainNodePtr = &(dFdx[li_Drain][ADrainEquDrainNodeOffset]);
  f_DrainEquDrainPrimeNodePtr = &(dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);
  f_SourceEquSourceNodePtr = &(dFdx[li_Source][ASourceEquSourceNodeOffset]);
  f_SourceEquSourcePrimeNodePtr = &(dFdx[li_Source][ASourceEquSourcePrimeNodeOffset]);
  if (bodyMod == 1)
  {
    f_BodyEquExtBodyNodePtr = &(dFdx[li_Body][ABodyEquExtBodyNodeOffset]);
    f_ExtBodyEquBodyNodePtr = &(dFdx[li_ExtBody][AExtBodyEquBodyNodeOffset]);
    f_ExtBodyEquExtBodyNodePtr = &(dFdx[li_ExtBody][AExtBodyEquExtBodyNodeOffset]);
  }
  if (selfheat)
  {
    f_DrainPrimeEquTemperatureNodePtr = &(dFdx[li_DrainPrime][ADrainPrimeEquTemperatureNodeOffset]);
    f_SourcePrimeEquTemperatureNodePtr = &(dFdx[li_SourcePrime][ASourcePrimeEquTemperatureNodeOffset]);

    f_GatePrimeEquTemperatureNodePtr = &(dFdx[li_GatePrime][AGatePrimeEquTemperatureNodeOffset]);
    f_TemperatureEquTemperatureNodePtr = &(dFdx[li_Temperature][ATemperatureEquTemperatureNodeOffset]);
    f_TemperatureEquGatePrimeNodePtr = &(dFdx[li_Temperature][ATemperatureEquGatePrimeNodeOffset]);
    f_TemperatureEquDrainPrimeNodePtr = &(dFdx[li_Temperature][ATemperatureEquDrainPrimeNodeOffset]);
    f_TemperatureEquSourcePrimeNodePtr = &(dFdx[li_Temperature][ATemperatureEquSourcePrimeNodeOffset]);
    if (soiMod != 0)
    {
      f_TemperatureEquSubstrateNodePtr = &(dFdx[li_Temperature][ATemperatureEquSubstrateNodeOffset]);
    }
    if (bNode > 0)
    {
      f_BodyEquTemperatureNodePtr = &(dFdx[li_Body][ABodyEquTemperatureNodeOffset]);
      f_TemperatureEquBodyNodePtr = &(dFdx[li_Temperature][ATemperatureEquBodyNodeOffset]);
    }
  }
  if( icVDSGiven )
  {


    f_DrainEquIdsPtr = &(dFdx[li_Drain][ADrainEquIdsOffset]);
    f_SourceEquIdsPtr = &(dFdx[li_Source][ASourceEquIdsOffset]);
    f_icVDSEquVdPtr = &(dFdx[li_Ids][icVDSEquVdOffset]);
    f_icVDSEquVsPtr = &(dFdx[li_Ids][icVDSEquVsOffset]);



    f_icVDSEquIdsPtr = &(dFdx[li_Ids][icVDSEquIdsOffset]);

  }

  if( icVGSGiven )
  {


    f_GateEquIgsPtr = &(dFdx[li_Gate][AGateEquIgsOffset]);
    f_SourceEquIgsPtr = &(dFdx[li_Source][ASourceEquIgsOffset]);
    f_icVGSEquVgPtr = &(dFdx[li_Igs][icVGSEquVgOffset]);
    f_icVGSEquVsPtr = &(dFdx[li_Igs][icVGSEquVsOffset]);



    f_icVGSEquIgsPtr = &(dFdx[li_Igs][icVGSEquIgsOffset]);

  }

  if( icVBSGiven )
  {


    f_BodyEquIbsPtr = &(dFdx[li_Body][ABodyEquIbsOffset]);
    f_SourceEquIbsPtr = &(dFdx[li_Source][ASourceEquIbsOffset]);
    f_icVBSEquVbPtr = &(dFdx[li_Ibs][icVBSEquVbOffset]);
    f_icVBSEquVsPtr = &(dFdx[li_Ibs][icVBSEquVsOffset]);



    f_icVBSEquIbsPtr = &(dFdx[li_Ibs][icVBSEquIbsOffset]);

  }

  if( icVESGiven )
  {


    f_SubstrateEquIesPtr = &(dFdx[li_Substrate][ASubstrateEquIesOffset]);
    f_SourceEquIesPtr = &(dFdx[li_Source][ASourceEquIesOffset]);
    f_icVESEquVePtr = &(dFdx[li_Ies][icVESEquVeOffset]);
    f_icVESEquVsPtr = &(dFdx[li_Ies][icVESEquVsOffset]);



    f_icVESEquIesPtr = &(dFdx[li_Ies][icVESEquIesOffset]);

  }

  if( icVPSGiven )
  {


    f_ExtBodyEquIpsPtr = &(dFdx[li_ExtBody][AExtBodyEquIpsOffset]);
    f_SourceEquIpsPtr = &(dFdx[li_Source][ASourceEquIpsOffset]);
    f_icVPSEquVpPtr = &(dFdx[li_Ips][icVPSEquVpOffset]);
    f_icVPSEquVsPtr = &(dFdx[li_Ips][icVPSEquVsOffset]);



    f_icVPSEquIpsPtr = &(dFdx[li_Ips][icVPSEquIpsOffset]);

  }


  // Q-matrix:
  if (rgateMod == 3)
  {
    q_GateMidEquGateMidNodePtr = &(dQdx[li_GateMid][AGateMidEquGateMidNodeOffset]);
    q_GateMidEquDrainPrimeNodePtr = &(dQdx[li_GateMid][AGateMidEquDrainPrimeNodeOffset]);

    q_GateMidEquSourcePrimeNodePtr = &(dQdx[li_GateMid][AGateMidEquSourcePrimeNodeOffset]);
    q_GateMidEquSubstrateNodePtr = &(dQdx[li_GateMid][AGateMidEquSubstrateNodeOffset]);

    q_DrainPrimeEquGateMidNodePtr = &(dQdx[li_DrainPrime][ADrainPrimeEquGateMidNodeOffset]);

    q_SourcePrimeEquGateMidNodePtr = &(dQdx[li_SourcePrime][ASourcePrimeEquGateMidNodeOffset]);
    q_SubstrateEquGateMidNodePtr = &(dQdx[li_Substrate][ASubstrateEquGateMidNodeOffset]);
  }

  q_SubstrateEquDrainPrimeNodePtr = &(dQdx[li_Substrate][ASubstrateEquDrainPrimeNodeOffset]);

  q_DrainPrimeEquSubstrateNodePtr = &(dQdx[li_DrainPrime][ADrainPrimeEquSubstrateNodeOffset]);
  q_SourcePrimeEquSubstrateNodePtr = &(dQdx[li_SourcePrime][ASourcePrimeEquSubstrateNodeOffset]);
  q_SubstrateEquGatePrimeNodePtr = &(dQdx[li_Substrate][ASubstrateEquGatePrimeNodeOffset]);
  q_GatePrimeEquSubstrateNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquSubstrateNodeOffset]);
  if (soiMod != 2)
  {
    q_SubstrateEquBodyNodePtr = &(dQdx[li_Substrate][ASubstrateEquBodyNodeOffset]);
    if (rgateMod == 0 || rgateMod == 1)
    {
      q_GatePrimeEquBodyNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquBodyNodeOffset]);
    }
    else
    {
      q_GatePrimeEquBodyNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquBodyNodeOffset]);
    }
    q_DrainPrimeEquBodyNodePtr = &(dQdx[li_DrainPrime][ADrainPrimeEquBodyNodeOffset]);
    q_SourcePrimeEquBodyNodePtr = &(dQdx[li_SourcePrime][ASourcePrimeEquBodyNodeOffset]);
    q_BodyEquSubstrateNodePtr = &(dQdx[li_Body][ABodyEquSubstrateNodeOffset]);
    q_BodyEquGatePrimeNodePtr = &(dQdx[li_Body][ABodyEquGatePrimeNodeOffset]);
    q_BodyEquDrainPrimeNodePtr = &(dQdx[li_Body][ABodyEquDrainPrimeNodeOffset]);
    q_BodyEquSourcePrimeNodePtr = &(dQdx[li_Body][ABodyEquSourcePrimeNodeOffset]);
    q_BodyEquBodyNodePtr = &(dQdx[li_Body][ABodyEquBodyNodeOffset]);
  }
  q_SubstrateEquSubstrateNodePtr = &(dQdx[li_Substrate][ASubstrateEquSubstrateNodeOffset]);
  if (rgateMod == 0)
  {
    q_GatePrimeEquGatePrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]);
    q_GatePrimeEquDrainPrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]);
    q_GatePrimeEquSourcePrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]);
  }
  else if (rgateMod == 1)
  {
    q_GatePrimeEquGatePrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]);
    q_GatePrimeEquDrainPrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]);
    q_GatePrimeEquSourcePrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]);
  }
  else
  {
    q_GatePrimeEquGatePrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]);
    q_GatePrimeEquDrainPrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]);
    q_GatePrimeEquSourcePrimeNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]);
  }
  q_DrainPrimeEquGatePrimeNodePtr = &(dQdx[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset]);
  q_DrainPrimeEquDrainPrimeNodePtr = &(dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  q_DrainPrimeEquSourcePrimeNodePtr = &(dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  q_SourcePrimeEquGatePrimeNodePtr = &(dQdx[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset]);
  q_SourcePrimeEquDrainPrimeNodePtr = &(dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  q_SourcePrimeEquSourcePrimeNodePtr = &(dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);

  if (selfheat)
  {
    q_DrainPrimeEquTemperatureNodePtr = &(dQdx[li_DrainPrime][ADrainPrimeEquTemperatureNodeOffset]);
    q_SourcePrimeEquTemperatureNodePtr = &(dQdx[li_SourcePrime][ASourcePrimeEquTemperatureNodeOffset]);
    q_SubstrateEquTemperatureNodePtr = &(dQdx[li_Substrate][ASubstrateEquTemperatureNodeOffset]);
    q_GatePrimeEquTemperatureNodePtr = &(dQdx[li_GatePrime][AGatePrimeEquTemperatureNodeOffset]);
    q_TemperatureEquTemperatureNodePtr = &(dQdx[li_Temperature][ATemperatureEquTemperatureNodeOffset]);

    if (bNode > 0)
    {
      q_BodyEquTemperatureNodePtr = &(dQdx[li_Body][ABodyEquTemperatureNodeOffset]);
    }
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       : This updates all the instance-owned paramters which
//                 are temperature dependent.
//
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool Instance::updateTemperature (const double & temp_tmp)
{

double tmp, tmp1, tmp2, T0, T1, T2, T3, T4, T5, Ldrn, Wdrn;
double TempRatio, Inv_L, Inv_W, Inv_LW, Tnom;
double SDphi, SDgamma;
double tmp3, T7, Eg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout << endl << dashedline2 << endl;
    cout << "Instance::updateTemperature\n";
    cout << "name = " << getName() << endl;
  }
#endif

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) temp = temp_tmp;

  if (model_.interpolateTNOM(temp))
  {
    // make sure interpolation doesn't take any resistance negative
    if (model_.sheetResistance < 0.0) model_.sheetResistance = 0.0;
    if (model_.rshg < 0.0) model_.rshg = 0.0;

    // some params may have changed during interpolation
    model_.processParams();
  }
// PMC
//  model_.outputParams(0);

  Tnom = model_.tnom;
  TempRatio = temp/Tnom;

  vtm = CONSTKoverQ * temp;
  Eg = 1.16 - 7.02e-4 * temp * temp / (temp + 1108.0);
  ni   = CONSTNi0 * (temp / CONSTREFTEMP) * sqrt(temp / CONSTREFTEMP)
    * exp(21.5565981 - Eg / (2.0 * vtm));

  rbodyext = bodySquares * model_.rbsh;

#ifdef REUSE_PARAMPTR
  // This next block determines whether  or not to use the model
  // size dependent parameters, or those of the instance.  If we are
  // using the instance, we may have to allocate the structure.

  list<SizeDependParam*>::iterator it_dpL =
    model_.sizeDependParamList.begin();
  list<SizeDependParam*>::iterator end_dpL =
    model_.sizeDependParamList.end();

  paramPtr = NULL;

  for( ; it_dpL != end_dpL; ++it_dpL ) {
    if( (*it_dpL)->Length == l && (*it_dpL)->Width == w &&
        (*it_dpL)->Rth0 == rth0 && (*it_dpL)->Cth0 == cth0) {
      paramPtr = (*it_dpL);
      break;
    }
  }

  if ( paramPtr != NULL )
  {
  }
  else
  {
    paramPtr = new SizeDependParam ();

    model_.sizeDependParamList.push_back( paramPtr );
#else
    if (paramPtr == static_cast<SizeDependParam *> (NULL))
      paramPtr = new SizeDependParam ();
#endif

    paramPtr->referenceTemperature = temp_tmp;

      Ldrn = l;
      Wdrn = w;
      paramPtr->Length = Ldrn;
      paramPtr->Width = Wdrn;
      paramPtr->Rth0 = rth0;
      paramPtr->Cth0 = cth0;

      T0 = pow(Ldrn, model_.Lln);
      T1 = pow(Wdrn, model_.Lwn);
      tmp1 = model_.Ll / T0 + model_.Lw / T1
           + model_.Lwl / (T0 * T1);
      paramPtr->dl = model_.Lint + tmp1;

// v2.2.3
      tmp1 = model_.Llc / T0 + model_.Lwc / T1
           + model_.Lwlc / (T0 * T1);
      paramPtr->dlc = model_.dlc + tmp1;

// v3.0
      paramPtr->dlcig = model_.dlcig + tmp1;


      T2 = pow(Ldrn, model_.Wln);
      T3 = pow(Wdrn, model_.Wwn);
      tmp2 = model_.Wl / T2 + model_.Ww / T3
           + model_.Wwl / (T2 * T3);
      paramPtr->dw = model_.Wint + tmp2;

// v2.2.3
      tmp2 = model_.Wlc / T2 + model_.Wwc / T3
           + model_.Wwlc / (T2 * T3);
      paramPtr->dwc = model_.dwc + tmp2;


      paramPtr->leff = l - 2.0 * paramPtr->dl;
      if (paramPtr->leff <= 0.0)
      {
        string msg = " In device: ";
        msg += getName();
        msg += ", Effective channel length <= 0";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }

      paramPtr->weff = w - nbc * model_.dwbc
         - (2.0 - nbc) * paramPtr->dw;
      if (paramPtr->weff <= 0.0)
      {
        string msg = " In device: ";
        msg += getName();
        msg += ", Effective channel length <= 0";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }

      paramPtr->wdiod = paramPtr->weff / nseg + pdbcp;
      paramPtr->wdios = paramPtr->weff / nseg + psbcp;

      paramPtr->leffCV = l - 2.0 * paramPtr->dlc;
      if (paramPtr->leffCV <= 0.0)
      {
        string msg = " In device: ";
        msg += getName();
        msg += ", Effective channel length for C-V <= 0";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }

      paramPtr->weffCV = w - nbc * model_.dwbc
         - (2.0 - nbc) * paramPtr->dwc;
      if (paramPtr->weffCV <= 0.0)
      {
        string msg = " In device: ";
        msg += getName();
        msg += ", Effective channel length for C-V <= 0";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }

      paramPtr->wdiodCV = paramPtr->weffCV / nseg + pdbcp;
      paramPtr->wdiosCV = paramPtr->weffCV / nseg + psbcp;

      paramPtr->leffCVb = l - 2.0 * paramPtr->dlc - model_.dlcb;
      if (paramPtr->leffCVb <= 0.0)
      {
        string msg = " In device: ";
        msg += getName();
        msg += ", Effective channel length for C-V (body) <= 0";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }

      paramPtr->leffCVbg = paramPtr->leffCVb + 2 * model_.dlbg;
      if (paramPtr->leffCVbg <= 0.0)
      {
        string msg = " In device: ";
        msg += getName();
        msg += ", Effective channel length for C-V (backgate) <= 0";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }


      /* Not binned - START */
      paramPtr->gamma1 = model_.gamma1;
      paramPtr->gamma2 = model_.gamma2;
      paramPtr->vbx = model_.vbx;
      paramPtr->vbm = model_.vbm;
      paramPtr->xt = model_.xt;
      /* Not binned - END */

      /* CV model */
      paramPtr->cf = model_.cf;
      paramPtr->clc = model_.clc;
      paramPtr->cle = model_.cle;

      paramPtr->abulkCVfactor = 1.0 + pow((paramPtr->clc / paramPtr->leff),
                                 paramPtr->cle);

      /* Added for binning - START */
      if (model_.binUnit == 1)
      {   Inv_L = 1.0e-6 / paramPtr->leff;
          Inv_W = 1.0e-6 / paramPtr->weff;
          Inv_LW = 1.0e-12 / (paramPtr->leff
                 * paramPtr->weff);
      }
      else
      {   Inv_L = 1.0 / paramPtr->leff;
          Inv_W = 1.0 / paramPtr->weff;
          Inv_LW = 1.0 / (paramPtr->leff
                 * paramPtr->weff);
      }
      paramPtr->npeak = model_.npeak
                         + model_.lnpeak * Inv_L
                         + model_.wnpeak * Inv_W
                         + model_.pnpeak * Inv_LW;
      paramPtr->nsub = model_.nsub
                        + model_.lnsub * Inv_L
                        + model_.wnsub * Inv_W
                        + model_.pnsub * Inv_LW;
      paramPtr->ngate = model_.ngate
                         + model_.lngate * Inv_L
                         + model_.wngate * Inv_W
                         + model_.pngate * Inv_LW;
      paramPtr->vth0 = model_.vth0
                        + model_.lvth0 * Inv_L
                        + model_.wvth0 * Inv_W
                        + model_.pvth0 * Inv_LW;
      paramPtr->k1 = model_.k1
                      + model_.lk1 * Inv_L
                      + model_.wk1 * Inv_W
                      + model_.pk1 * Inv_LW;
      paramPtr->k2 = model_.k2
                      + model_.lk2 * Inv_L
                      + model_.wk2 * Inv_W
                      + model_.pk2 * Inv_LW;
      paramPtr->k1w1 = model_.k1w1
                      + model_.lk1w1 * Inv_L
                      + model_.wk1w1 * Inv_W
                      + model_.pk1w1 * Inv_LW;
      paramPtr->k1w2 = model_.k1w2
                      + model_.lk1w2 * Inv_L
                      + model_.wk1w2 * Inv_W
                      + model_.pk1w2 * Inv_LW;
      paramPtr->k3 = model_.k3
                      + model_.lk3 * Inv_L
                      + model_.wk3 * Inv_W
                      + model_.pk3 * Inv_LW;
      paramPtr->k3b = model_.k3b
                       + model_.lk3b * Inv_L
                       + model_.wk3b * Inv_W
                       + model_.pk3b * Inv_LW;
      paramPtr->kb1 = model_.kb1
                       + model_.lkb1 * Inv_L
                       + model_.wkb1 * Inv_W
                       + model_.pkb1 * Inv_LW;
      paramPtr->w0 = model_.w0
                      + model_.lw0 * Inv_L
                      + model_.ww0 * Inv_W
                      + model_.pw0 * Inv_LW;
      paramPtr->nlx = model_.nlx
                       + model_.lnlx * Inv_L
                       + model_.wnlx * Inv_W
                       + model_.pnlx * Inv_LW;
      paramPtr->dvt0 = model_.dvt0
                        + model_.ldvt0 * Inv_L
                        + model_.wdvt0 * Inv_W
                        + model_.pdvt0 * Inv_LW;
      paramPtr->dvt1 = model_.dvt1
                        + model_.ldvt1 * Inv_L
                        + model_.wdvt1 * Inv_W
                        + model_.pdvt1 * Inv_LW;
      paramPtr->dvt2 = model_.dvt2
                        + model_.ldvt2 * Inv_L
                        + model_.wdvt2 * Inv_W
                        + model_.pdvt2 * Inv_LW;
      paramPtr->dvt0w = model_.dvt0w
                        + model_.ldvt0w * Inv_L
                        + model_.wdvt0w * Inv_W
                        + model_.pdvt0w * Inv_LW;
      paramPtr->dvt1w = model_.dvt1w
                        + model_.ldvt1w * Inv_L
                        + model_.wdvt1w * Inv_W
                        + model_.pdvt1w * Inv_LW;
      paramPtr->dvt2w = model_.dvt2w
                        + model_.ldvt2w * Inv_L
                        + model_.wdvt2w * Inv_W
                        + model_.pdvt2w * Inv_LW;
      paramPtr->u0 = model_.u0
                      + model_.lu0 * Inv_L
                      + model_.wu0 * Inv_W
                      + model_.pu0 * Inv_LW;
      paramPtr->ua = model_.ua
                      + model_.lua * Inv_L
                      + model_.wua * Inv_W
                      + model_.pua * Inv_LW;
      paramPtr->ub = model_.ub
                      + model_.lub * Inv_L
                      + model_.wub * Inv_W
                      + model_.pub * Inv_LW;
      paramPtr->uc = model_.uc
                      + model_.luc * Inv_L
                      + model_.wuc * Inv_W
                      + model_.puc * Inv_LW;
      paramPtr->vsat = model_.vsat
                        + model_.lvsat * Inv_L
                        + model_.wvsat * Inv_W
                        + model_.pvsat * Inv_LW;
      paramPtr->a0 = model_.a0
                      + model_.la0 * Inv_L
                      + model_.wa0 * Inv_W
                      + model_.pa0 * Inv_LW;
      paramPtr->ags = model_.ags
                      + model_.lags * Inv_L
                      + model_.wags * Inv_W
                      + model_.pags * Inv_LW;
      paramPtr->b0 = model_.b0
                      + model_.lb0 * Inv_L
                      + model_.wb0 * Inv_W
                      + model_.pb0 * Inv_LW;
      paramPtr->b1 = model_.b1
                      + model_.lb1 * Inv_L
                      + model_.wb1 * Inv_W
                      + model_.pb1 * Inv_LW;
      paramPtr->keta = model_.keta
                        + model_.lketa * Inv_L
                        + model_.wketa * Inv_W
                        + model_.pketa * Inv_LW;
      paramPtr->ketas = model_.ketas
                        + model_.lketas * Inv_L
                        + model_.wketas * Inv_W
                        + model_.pketas * Inv_LW;
      paramPtr->a1 = model_.a1
                      + model_.la1 * Inv_L
                      + model_.wa1 * Inv_W
                      + model_.pa1 * Inv_LW;
      paramPtr->a2 = model_.a2
                      + model_.la2 * Inv_L
                      + model_.wa2 * Inv_W
                      + model_.pa2 * Inv_LW;
      paramPtr->rdsw = model_.rdsw
                        + model_.lrdsw * Inv_L
                        + model_.wrdsw * Inv_W
                        + model_.prdsw * Inv_LW;
      paramPtr->prwb = model_.prwb
                        + model_.lprwb * Inv_L
                        + model_.wprwb * Inv_W
                        + model_.pprwb * Inv_LW;
      paramPtr->prwg = model_.prwg
                        + model_.lprwg * Inv_L
                        + model_.wprwg * Inv_W
                        + model_.pprwg * Inv_LW;
      paramPtr->wr = model_.wr
                      + model_.lwr * Inv_L
                      + model_.wwr * Inv_W
                      + model_.pwr * Inv_LW;
      paramPtr->nfactor = model_.nfactor
                           + model_.lnfactor * Inv_L
                           + model_.wnfactor * Inv_W
                           + model_.pnfactor * Inv_LW;
      paramPtr->dwg = model_.dwg
                       + model_.ldwg * Inv_L
                       + model_.wdwg * Inv_W
                       + model_.pdwg * Inv_LW;
      paramPtr->dwb = model_.dwb
                       + model_.ldwb * Inv_L
                       + model_.wdwb * Inv_W
                       + model_.pdwb * Inv_LW;
      paramPtr->voff = model_.voff
                        + model_.lvoff * Inv_L
                        + model_.wvoff * Inv_W
                        + model_.pvoff * Inv_LW;
      paramPtr->eta0 = model_.eta0
                        + model_.leta0 * Inv_L
                        + model_.weta0 * Inv_W
                        + model_.peta0 * Inv_LW;
      paramPtr->etab = model_.etab
                        + model_.letab * Inv_L
                        + model_.wetab * Inv_W
                        + model_.petab * Inv_LW;
      paramPtr->dsub = model_.dsub
                        + model_.ldsub * Inv_L
                        + model_.wdsub * Inv_W
                        + model_.pdsub * Inv_LW;
      paramPtr->cit = model_.cit
                       + model_.lcit * Inv_L
                       + model_.wcit * Inv_W
                       + model_.pcit * Inv_LW;
      paramPtr->cdsc = model_.cdsc
                        + model_.lcdsc * Inv_L
                        + model_.wcdsc * Inv_W
                        + model_.pcdsc * Inv_LW;
      paramPtr->cdscb = model_.cdscb
                         + model_.lcdscb * Inv_L
                         + model_.wcdscb * Inv_W
                         + model_.pcdscb * Inv_LW;
      paramPtr->cdscd = model_.cdscd
                         + model_.lcdscd * Inv_L
                         + model_.wcdscd * Inv_W
                         + model_.pcdscd * Inv_LW;
      paramPtr->pclm = model_.pclm
                        + model_.lpclm * Inv_L
                        + model_.wpclm * Inv_W
                        + model_.ppclm * Inv_LW;
      paramPtr->pdibl1 = model_.pdibl1
                          + model_.lpdibl1 * Inv_L
                          + model_.wpdibl1 * Inv_W
                          + model_.ppdibl1 * Inv_LW;
      paramPtr->pdibl2 = model_.pdibl2
                          + model_.lpdibl2 * Inv_L
                          + model_.wpdibl2 * Inv_W
                          + model_.ppdibl2 * Inv_LW;
      paramPtr->pdiblb = model_.pdiblb
                          + model_.lpdiblb * Inv_L
                          + model_.wpdiblb * Inv_W
                          + model_.ppdiblb * Inv_LW;
      paramPtr->drout = model_.drout
                         + model_.ldrout * Inv_L
                         + model_.wdrout * Inv_W
                         + model_.pdrout * Inv_LW;
      paramPtr->pvag = model_.pvag
                        + model_.lpvag * Inv_L
                        + model_.wpvag * Inv_W
                        + model_.ppvag * Inv_LW;
      paramPtr->delta = model_.delta
                         + model_.ldelta * Inv_L
                         + model_.wdelta * Inv_W
                         + model_.pdelta * Inv_LW;
      paramPtr->alpha0 = model_.alpha0
                          + model_.lalpha0 * Inv_L
                          + model_.walpha0 * Inv_W
                          + model_.palpha0 * Inv_LW;
      paramPtr->fbjtii = model_.fbjtii
                          + model_.lfbjtii * Inv_L
                          + model_.wfbjtii * Inv_W
                          + model_.pfbjtii * Inv_LW;
      paramPtr->beta0 = model_.beta0
                         + model_.lbeta0 * Inv_L
                         + model_.wbeta0 * Inv_W
                         + model_.pbeta0 * Inv_LW;
      paramPtr->beta1 = model_.beta1
                         + model_.lbeta1 * Inv_L
                         + model_.wbeta1 * Inv_W
                         + model_.pbeta1 * Inv_LW;
      paramPtr->beta2 = model_.beta2
                         + model_.lbeta2 * Inv_L
                         + model_.wbeta2 * Inv_W
                         + model_.pbeta2 * Inv_LW;
      paramPtr->vdsatii0 = model_.vdsatii0
                          + model_.lvdsatii0 * Inv_L
                          + model_.wvdsatii0 * Inv_W
                          + model_.pvdsatii0 * Inv_LW;
      paramPtr->lii = model_.lii
                          + model_.llii * Inv_L
                          + model_.wlii * Inv_W
                          + model_.plii * Inv_LW;
      paramPtr->esatii = model_.esatii
                          + model_.lesatii * Inv_L
                          + model_.wesatii * Inv_W
                          + model_.pesatii * Inv_LW;
      paramPtr->sii0 = model_.sii0
                          + model_.lsii0 * Inv_L
                          + model_.wsii0 * Inv_W
                          + model_.psii0 * Inv_LW;
      paramPtr->sii1 = model_.sii1
                          + model_.lsii1 * Inv_L
                          + model_.wsii1 * Inv_W
                          + model_.psii1 * Inv_LW;
      paramPtr->sii2 = model_.sii2
                          + model_.lsii2 * Inv_L
                          + model_.wsii2 * Inv_W
                          + model_.psii2 * Inv_LW;
      paramPtr->siid = model_.siid
                          + model_.lsiid * Inv_L
                          + model_.wsiid * Inv_W
                          + model_.psiid * Inv_LW;
      paramPtr->agidl = model_.agidl
                          + model_.lagidl * Inv_L
                          + model_.wagidl * Inv_W
                          + model_.pagidl * Inv_LW;
      paramPtr->bgidl = model_.bgidl
                          + model_.lbgidl * Inv_L
                          + model_.wbgidl * Inv_W
                          + model_.pbgidl * Inv_LW;
      paramPtr->ngidl = model_.ngidl
                          + model_.lngidl * Inv_L
                          + model_.wngidl * Inv_W
                          + model_.pngidl * Inv_LW;
      paramPtr->ntun = model_.ntun
                          + model_.lntun * Inv_L
                          + model_.wntun * Inv_W
                          + model_.pntun * Inv_LW;
      paramPtr->ndiode = model_.ndiode
                          + model_.lndiode * Inv_L
                          + model_.wndiode * Inv_W
                          + model_.pndiode * Inv_LW;
      paramPtr->nrecf0 = model_.nrecf0
                      + model_.lnrecf0 * Inv_L
                      + model_.wnrecf0 * Inv_W
                      + model_.pnrecf0 * Inv_LW;
      paramPtr->nrecr0 = model_.nrecr0
                      + model_.lnrecr0 * Inv_L
                      + model_.wnrecr0 * Inv_W
                      + model_.pnrecr0 * Inv_LW;
      paramPtr->isbjt = model_.isbjt
                      + model_.lisbjt * Inv_L
                      + model_.wisbjt * Inv_W
                      + model_.pisbjt * Inv_LW;
      paramPtr->isdif = model_.isdif
                      + model_.lisdif * Inv_L
                      + model_.wisdif * Inv_W
                      + model_.pisdif * Inv_LW;
      paramPtr->isrec = model_.isrec
                      + model_.lisrec * Inv_L
                      + model_.wisrec * Inv_W
                      + model_.pisrec * Inv_LW;
      paramPtr->istun = model_.istun
                      + model_.listun * Inv_L
                      + model_.wistun * Inv_W
                      + model_.pistun * Inv_LW;
      paramPtr->vrec0 = model_.vrec0
                      + model_.lvrec0 * Inv_L
                      + model_.wvrec0 * Inv_W
                      + model_.pvrec0 * Inv_LW;
      paramPtr->vtun0 = model_.vtun0
                      + model_.lvtun0 * Inv_L
                      + model_.wvtun0 * Inv_W
                      + model_.pvtun0 * Inv_LW;
      paramPtr->nbjt = model_.nbjt
                      + model_.lnbjt * Inv_L
                      + model_.wnbjt * Inv_W
                      + model_.pnbjt * Inv_LW;
      paramPtr->lbjt0 = model_.lbjt0
                      + model_.llbjt0 * Inv_L
                      + model_.wlbjt0 * Inv_W
                      + model_.plbjt0 * Inv_LW;
      paramPtr->vabjt = model_.vabjt
                      + model_.lvabjt * Inv_L
                      + model_.wvabjt * Inv_W
                      + model_.pvabjt * Inv_LW;
      paramPtr->aely = model_.aely
                      + model_.laely * Inv_L
                      + model_.waely * Inv_W
                      + model_.paely * Inv_LW;
      paramPtr->ahli = model_.ahli
                      + model_.lahli * Inv_L
                      + model_.wahli * Inv_W
                      + model_.pahli * Inv_LW;


      paramPtr->xj = model_.xj
                         + model_.lxj * Inv_L
                         + model_.wxj * Inv_W
                         + model_.pxj * Inv_LW;
      paramPtr->alphaGB1 = model_.alphaGB1
                         + model_.lalphaGB1 * Inv_L
                         + model_.walphaGB1 * Inv_W
                         + model_.palphaGB1 * Inv_LW;
      paramPtr->alphaGB2 = model_.alphaGB2
                         + model_.lalphaGB2 * Inv_L
                         + model_.walphaGB2 * Inv_W
                         + model_.palphaGB2 * Inv_LW;
      paramPtr->betaGB1 = model_.betaGB1
                         + model_.lbetaGB1* Inv_L
                         + model_.wbetaGB1 * Inv_W
                         + model_.pbetaGB1 * Inv_LW;
      paramPtr->betaGB2 = model_.betaGB2
                         + model_.lbetaGB2 * Inv_L
                         + model_.wbetaGB2 * Inv_W
                         + model_.pbetaGB2 * Inv_LW;
      paramPtr->ndif = model_.ndif
                         + model_.lndif * Inv_L
                         + model_.wndif * Inv_W
                         + model_.pndif * Inv_LW;
      paramPtr->ntrecf = model_.ntrecf
                         + model_.lntrecf* Inv_L
                         + model_.wntrecf * Inv_W
                         + model_.pntrecf * Inv_LW;
      paramPtr->ntrecr = model_.ntrecr
                         + model_.lntrecr * Inv_L
                         + model_.wntrecr * Inv_W
                         + model_.pntrecr * Inv_LW;
      paramPtr->xbjt = model_.xbjt
                         + model_.lxbjt * Inv_L
                         + model_.wxbjt * Inv_W
                         + model_.pxbjt * Inv_LW;
      paramPtr->xdif = model_.xdif
                         + model_.lxdif* Inv_L
                         + model_.wxdif * Inv_W
                         + model_.pxdif * Inv_LW;
      paramPtr->xrec = model_.xrec
                         + model_.lxrec * Inv_L
                         + model_.wxrec * Inv_W
                         + model_.pxrec * Inv_LW;
      paramPtr->xtun = model_.xtun
                         + model_.lxtun * Inv_L
                         + model_.wxtun * Inv_W
                         + model_.pxtun * Inv_LW;
      paramPtr->cgdl = model_.cgdl
                         + model_.lcgdl * Inv_L
                         + model_.wcgdl * Inv_W
                         + model_.pcgdl * Inv_LW;
      paramPtr->cgsl = model_.cgsl
                         + model_.lcgsl * Inv_L
                         + model_.wcgsl * Inv_W
                         + model_.pcgsl * Inv_LW;
      paramPtr->ckappa = model_.ckappa
                         + model_.lckappa * Inv_L
                         + model_.wckappa * Inv_W
                         + model_.pckappa * Inv_LW;
      paramPtr->ute = model_.ute
                         + model_.lute * Inv_L
                         + model_.wute * Inv_W
                         + model_.pute * Inv_LW;
      paramPtr->kt1 = model_.kt1
                         + model_.lkt1 * Inv_L
                         + model_.wkt1 * Inv_W
                         + model_.pkt1 * Inv_LW;
      paramPtr->kt2 = model_.kt2
                         + model_.lkt2 * Inv_L
                         + model_.wkt2 * Inv_W
                         + model_.pkt2 * Inv_LW;
      paramPtr->kt1l = model_.kt1l
                         + model_.lkt1l * Inv_L
                         + model_.wkt1l * Inv_W
                         + model_.pkt1l * Inv_LW;
      paramPtr->ua1 = model_.ua1
                         + model_.lua1 * Inv_L
                         + model_.wua1 * Inv_W
                         + model_.pua1 * Inv_LW;
      paramPtr->ub1 = model_.ub1
                         + model_.lub1* Inv_L
                         + model_.wub1 * Inv_W
                         + model_.pub1 * Inv_LW;
      paramPtr->uc1 = model_.uc1
                         + model_.luc1 * Inv_L
                         + model_.wuc1 * Inv_W
                         + model_.puc1 * Inv_LW;
      paramPtr->at = model_.at
                         + model_.lat * Inv_L
                         + model_.wat * Inv_W
                         + model_.pat * Inv_LW;
      paramPtr->prt = model_.prt
                         + model_.lprt * Inv_L
                         + model_.wprt * Inv_W
                         + model_.pprt * Inv_LW;


      paramPtr->nigc = model_.nigc
                         + model_.lnigc * Inv_L
                         + model_.wnigc * Inv_W
                         + model_.pnigc * Inv_LW;
      paramPtr->aigc = model_.aigc
                         + model_.laigc * Inv_L
                         + model_.waigc * Inv_W
                         + model_.paigc * Inv_LW;
      paramPtr->bigc = model_.bigc
                         + model_.lbigc * Inv_L
                         + model_.wbigc * Inv_W
                         + model_.pbigc * Inv_LW;
      paramPtr->cigc = model_.cigc
                         + model_.lcigc * Inv_L
                         + model_.wcigc * Inv_W
                         + model_.pcigc * Inv_LW;
      paramPtr->aigsd = model_.aigsd
                         + model_.laigsd * Inv_L
                         + model_.waigsd * Inv_W
                         + model_.paigsd * Inv_LW;
      paramPtr->bigsd = model_.bigsd
                         + model_.lbigsd * Inv_L
                         + model_.wbigsd * Inv_W
                         + model_.pbigsd * Inv_LW;
      paramPtr->cigsd = model_.cigsd
                         + model_.lcigsd * Inv_L
                         + model_.wcigsd * Inv_W
                         + model_.pcigsd * Inv_LW;
      paramPtr->pigcd = model_.pigcd
                         + model_.lpigcd * Inv_L
                         + model_.wpigcd * Inv_W
                         + model_.ppigcd * Inv_LW;
      paramPtr->poxedge = model_.poxedge
                           + model_.lpoxedge * Inv_L
                           + model_.wpoxedge * Inv_W
                           + model_.ppoxedge * Inv_LW;
// v3.0

// v3.1 wanh added for RF
      paramPtr->xrcrg1 = model_.xrcrg1
                          + model_.lxrcrg1 * Inv_L
                          + model_.wxrcrg1 * Inv_W
                          + model_.pxrcrg1 * Inv_LW;
      paramPtr->xrcrg2 = model_.xrcrg2
                          + model_.lxrcrg2 * Inv_L
                          + model_.wxrcrg2 * Inv_W
                          + model_.pxrcrg2 * Inv_LW;
// v3.1 wanh added for RF end


      // CV model
      paramPtr->vsdfb = model_.vsdfb
                      + model_.lvsdfb * Inv_L
                      + model_.wvsdfb * Inv_W
                      + model_.pvsdfb * Inv_LW;
      paramPtr->vsdth = model_.vsdth
                      + model_.lvsdth * Inv_L
                      + model_.wvsdth * Inv_W
                      + model_.pvsdth * Inv_LW;
      paramPtr->delvt = model_.delvt
                      + model_.ldelvt * Inv_L
                      + model_.wdelvt * Inv_W
                      + model_.pdelvt * Inv_LW;
      paramPtr->acde = model_.acde
                      + model_.lacde * Inv_L
                      + model_.wacde * Inv_W
                      + model_.pacde * Inv_LW;
      paramPtr->acde = paramPtr->acde *
                        pow((paramPtr->npeak / 2.0e16), -0.25);
      // v3.2 bug fix

      paramPtr->moin = model_.moin
                      + model_.lmoin * Inv_L
                      + model_.wmoin * Inv_W
                      + model_.pmoin * Inv_LW;
      paramPtr->noff = model_.noff
                        + model_.lnoff * Inv_L
                        + model_.wnoff * Inv_W
                        + model_.pnoff * Inv_LW; // v3.2
      /* Added for binning - END */

      T0 = (TempRatio - 1.0);

      paramPtr->uatemp = paramPtr->ua;  /*  save ua, ub, and uc for b3soild.c */
      paramPtr->ubtemp = paramPtr->ub;
      paramPtr->uctemp = paramPtr->uc;
      paramPtr->rds0denom = pow(paramPtr->weff * 1E6, paramPtr->wr);


/* v2.2 release */
      paramPtr->rth = rth0 / (paramPtr->weff + model_.wth0)
                       * nseg;
      paramPtr->cth = cth0 * (paramPtr->weff + model_.wth0)
                       / nseg;

/* v2.2.2 adding layout-dependent Frbody multiplier */
      paramPtr->rbody = frbody *model_.rbody * model_.rhalo
        / (2 * model_.rbody + model_.rhalo * paramPtr->leff)
        * paramPtr->weff / nseg;

      paramPtr->oxideRatio = pow(model_.toxref/model_.toxqm,
                                 model_.ntox) /
        model_.toxqm/model_.toxqm;
/* v2.2 release */


      paramPtr->ua = paramPtr->ua + paramPtr->ua1 * T0;
      paramPtr->ub = paramPtr->ub + paramPtr->ub1 * T0;
      paramPtr->uc = paramPtr->uc + paramPtr->uc1 * T0;
      if (paramPtr->u0 > 1.0)
          paramPtr->u0 = paramPtr->u0 / 1.0e4;

      paramPtr->u0temp = paramPtr->u0
                          * pow(TempRatio, paramPtr->ute);
      paramPtr->vsattemp = paramPtr->vsat - paramPtr->at
                            * T0;
      paramPtr->rds0 = (paramPtr->rdsw + paramPtr->prt * T0)
                        / pow(paramPtr->weff * 1E6, paramPtr->wr);

      if (!checkModel ())
      {
        string msg = " In device: ";
        msg += getName();
        msg += ", detected during B3SOIV3 parameter check";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }

      paramPtr->cgdo = (model_.cgdo + paramPtr->cf)
                        * paramPtr->wdiodCV;
      paramPtr->cgso = (model_.cgso + paramPtr->cf)
                        * paramPtr->wdiosCV;

      paramPtr->cgeo = model_.cgeo
                        * paramPtr->leffCV;


      if (!model_.npeakGiven && model_.gamma1Given)
      {   T0 = paramPtr->gamma1 * model_.cox;
          paramPtr->npeak = 3.021E22 * T0 * T0;
      }


      T4 = CONSTEg300 / vtm * (TempRatio - 1.0);
      T7 = paramPtr->xbjt * T4 / paramPtr->ndiode;
      DEXP(T7, T0);
      T7 = paramPtr->xdif * T4 / paramPtr->ndiode;
      DEXP(T7, T1);
      T7 = paramPtr->xrec * T4 / paramPtr->nrecf0;
      DEXP(T7, T2);

      /* v2.2.2 bug fix */
      paramPtr->ahli0 = paramPtr->ahli * T0;

      paramPtr->jbjt = paramPtr->isbjt * T0;
      paramPtr->jdif = paramPtr->isdif * T1;
      paramPtr->jrec = paramPtr->isrec * T2;

      T7 = paramPtr->xtun * (TempRatio - 1);
      DEXP(T7, T0);
      paramPtr->jtun = paramPtr->istun * T0;


      if (paramPtr->nsub > 0)
         paramPtr->vfbb = -model_.dtype * vtm *
                    log(paramPtr->npeak/ paramPtr->nsub);
      else
         paramPtr->vfbb = -model_.dtype * vtm *
                    log(-paramPtr->npeak* paramPtr->nsub/
                    (ni*ni));

      if (!model_.vsdfbGiven)
      {
         if (paramPtr->nsub > 0)
            paramPtr->vsdfb = -model_.dtype * (vtm*log(1e20 *
                            paramPtr->nsub /
                           (ni*ni)) - 0.3);
         else if (paramPtr->nsub < 0)
            paramPtr->vsdfb = -model_.dtype * (vtm*log(-1e20 /
                                paramPtr->nsub) + 0.3);
      }

      /* Phi  & Gamma */
      SDphi = 2.0*vtm*log(fabs(paramPtr->nsub) / ni);
      SDgamma = 5.753e-12 * sqrt(fabs(paramPtr->nsub)) / model_.cbox;

      if (!model_.vsdthGiven)
      {
         if ( ((paramPtr->nsub > 0) && (model_.dtype > 0)) ||
              ((paramPtr->nsub < 0) && (model_.dtype < 0)) )
            paramPtr->vsdth = paramPtr->vsdfb + SDphi +
                                SDgamma * sqrt(SDphi);
         else
            paramPtr->vsdth = paramPtr->vsdfb - SDphi -
                                SDgamma * sqrt(SDphi);
      }

      if (!model_.csdminGiven)
      {
         /* Cdmin */
         tmp = sqrt(2.0 * CONSTEPSSI * SDphi / (CONSTQ *
                    fabs(paramPtr->nsub) * 1.0e6));
         tmp1 = CONSTEPSSI / tmp;
         model_.csdmin = tmp1 * model_.cbox /
                              (tmp1 + model_.cbox);
      }


      paramPtr->phi = 2.0 * vtm * log(paramPtr->npeak / ni);

      paramPtr->sqrtPhi = sqrt(paramPtr->phi);
      paramPtr->phis3 = paramPtr->sqrtPhi * paramPtr->phi;

      paramPtr->Xdep0 = sqrt(2.0 * CONSTEPSSI / (CONSTQ
                         * paramPtr->npeak * 1.0e6))
                         * paramPtr->sqrtPhi;
      paramPtr->sqrtXdep0 = sqrt(paramPtr->Xdep0);
      paramPtr->litl = sqrt(3.0 * paramPtr->xj
                        * model_.tox);
      paramPtr->vbi = vtm * log(1.0e20
                       * paramPtr->npeak / (ni*ni));
      paramPtr->cdep0 = sqrt(CONSTQ * CONSTEPSSI
                         * paramPtr->npeak * 1.0e6 / 2.0
                         / paramPtr->phi);

// v3.0
      if (paramPtr->ngate > 0.0)
      {   paramPtr->vfbsd = model_.Vtm0 * log(paramPtr->ngate
                             / 1.0e20);
      }
      else
          paramPtr->vfbsd = 0.0;

      paramPtr->ToxRatio = exp(model_.ntox
                            * log(model_.toxref /model_.toxqm))
                            /model_.toxqm /model_.toxqm;
      paramPtr->ToxRatioEdge = exp(model_.ntox
                                * log(model_.toxref
                                / (model_.toxqm * paramPtr->poxedge)))
                                / model_.toxqm / model_.toxqm
                                / paramPtr->poxedge / paramPtr->poxedge;
      paramPtr->Aechvb = (model_.dtype==CONSTNMOS)?4.97232e-7:3.42537e-7;
      paramPtr->Bechvb = (model_.dtype==CONSTNMOS)?7.45669e11:1.16645e12;
      paramPtr->AechvbEdge = paramPtr->Aechvb * paramPtr->weff/nseg
        * paramPtr->dlcig * paramPtr->ToxRatioEdge; // v3.1 bug fix
      paramPtr->BechvbEdge = -paramPtr->Bechvb
                              * model_.toxqm * paramPtr->poxedge;
      paramPtr->Aechvb *= paramPtr->weff/nseg * paramPtr->leff
                           * paramPtr->ToxRatio; // v3.1 bug fix
      paramPtr->Bechvb *= -model_.toxqm;
// v3.0

      if (model_.k1Given || model_.k2Given)
      {
          if (!model_.k1Given)
          {   string msg = "Warning: k1 should be specified with k2.";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
              paramPtr->k1 = 0.53;
          }
          if (!model_.k2Given)
          {   string msg = "Warning: k2 should be specified with k1.";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
              paramPtr->k2 = -0.0186;
          }
          if (model_.xtGiven)
          {
              string msg = "Warning: xt is ignored because k1 or k2 is given.";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
          }
          if (model_.vbxGiven)
          {
              string msg ="Warning: vbx is ignored because k1 or k2 is given.";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
          }
          if (model_.vbmGiven)
          {
              string msg ="Warning: vbm is ignored because k1 or k2 is given.";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
          }
          if (model_.gamma1Given)
          {
              string msg="Warning: gamma1 ignored because k1 or k2 is given.";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
          }
          if (model_.gamma2Given)
          {
              string msg="Warning: gamma2 ignored because k1 or k2 is given.";
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
          }
      }
      else
      {
          if (!model_.vbxGiven)
              paramPtr->vbx = paramPtr->phi - 7.7348e-4
                               * paramPtr->npeak
                               * paramPtr->xt * paramPtr->xt;
          if (paramPtr->vbx > 0.0)
              paramPtr->vbx = -paramPtr->vbx;
          if (paramPtr->vbm > 0.0)
              paramPtr->vbm = -paramPtr->vbm;

          if (!model_.gamma1Given)
              paramPtr->gamma1 = 5.753e-12
                                  * sqrt(paramPtr->npeak)
                                  / model_.cox;
          if (!model_.gamma2Given)
              paramPtr->gamma2 = 5.753e-12
                                  * sqrt(paramPtr->nsub)
                                  / model_.cox;

          T0 = paramPtr->gamma1 - paramPtr->gamma2;
          T1 = sqrt(paramPtr->phi - paramPtr->vbx)
             - paramPtr->sqrtPhi;
          T2 = sqrt(paramPtr->phi * (paramPtr->phi
             - paramPtr->vbm)) - paramPtr->phi;
          paramPtr->k2 = T0 * T1 / (2.0 * T2 + paramPtr->vbm);
          paramPtr->k1 = paramPtr->gamma2 - 2.0
                          * paramPtr->k2 * sqrt(paramPtr->phi
                          - paramPtr->vbm);
      }

      if (paramPtr->k2 < 0.0)
      {   T0 = 0.5 * paramPtr->k1 / paramPtr->k2;
          paramPtr->vbsc = 0.9 * (paramPtr->phi - T0 * T0);
          if (paramPtr->vbsc > -3.0)
              paramPtr->vbsc = -3.0;
          else if (paramPtr->vbsc < -30.0)
              paramPtr->vbsc = -30.0;
      }
      else
      {   paramPtr->vbsc = -30.0;
      }
      if (paramPtr->vbsc > paramPtr->vbm)
          paramPtr->vbsc = paramPtr->vbm;

      if ((T0 = paramPtr->weff + paramPtr->k1w2) < 1e-8)
         T0 = 1e-8;
      paramPtr->k1eff = paramPtr->k1 * (1 + paramPtr->k1w1/T0);

      if (model_.vth0Given)
      {   paramPtr->vfb = model_.dtype * paramPtr->vth0
                           - paramPtr->phi - paramPtr->k1eff
                           * paramPtr->sqrtPhi;
      }
      else
      {   paramPtr->vfb = -1.0;
          paramPtr->vth0 = model_.dtype * (paramPtr->vfb
                            + paramPtr->phi + paramPtr->k1eff
                            * paramPtr->sqrtPhi);
      }

// v3.2
      paramPtr->k1eff *= model_.tox / model_.toxm;
      paramPtr->k2 *= model_.tox / model_.toxm;


      T1 = sqrt(CONSTEPSSI / CONSTEPSOX * model_.tox
         * paramPtr->Xdep0);
      T0 = exp(-0.5 * paramPtr->dsub * paramPtr->leff / T1);
      paramPtr->theta0vb0 = (T0 + 2.0 * T0 * T0);

      T0 = exp(-0.5 * paramPtr->drout * paramPtr->leff / T1);
      T2 = (T0 + 2.0 * T0 * T0);
      paramPtr->thetaRout = paramPtr->pdibl1 * T2
                             + paramPtr->pdibl2;
#ifdef REUSE_PARAMPTR
  }
#endif

  csbox = model_.cbox*sourceArea;
  csmin = model_.csdmin*sourceArea;
  cdbox = model_.cbox*drainArea;
  cdmin = model_.csdmin*drainArea;

  if ( ((paramPtr->nsub > 0) && (model_.dtype > 0)) ||
       ((paramPtr->nsub < 0) && (model_.dtype < 0)) )
  {
     T0 = paramPtr->vsdth - paramPtr->vsdfb;
     paramPtr->sdt1 = paramPtr->vsdfb + model_.asd * T0;
     T1 = csbox - csmin;
     T2 = T1 / T0 / T0;
     paramPtr->st2 = T2 / model_.asd;
     paramPtr->st3 = T2 /( 1 - model_.asd);
     st4 =  T0 * T1 * (1 + model_.asd) / 3
                      - csmin * paramPtr->vsdfb;

     T1 = cdbox - cdmin;
     T2 = T1 / T0 / T0;
     paramPtr->dt2 = T2 / model_.asd;
     paramPtr->dt3 = T2 /( 1 - model_.asd);
     dt4 =  T0 * T1 * (1 + model_.asd) / 3
                      - cdmin * paramPtr->vsdfb;
  } else
  {
     T0 = paramPtr->vsdfb - paramPtr->vsdth;
     paramPtr->sdt1 = paramPtr->vsdth + model_.asd * T0;
     T1 = csmin - csbox;
     T2 = T1 / T0 / T0;
     paramPtr->st2 = T2 / model_.asd;
     paramPtr->st3 = T2 /( 1 - model_.asd);
     st4 =  T0 * T1 * (1 + model_.asd) / 3
                      - csbox * paramPtr->vsdth;

     T1 = cdmin - cdbox;
     T2 = T1 / T0 / T0;
     paramPtr->dt2 = T2 / model_.asd;
     paramPtr->dt3 = T2 /( 1 - model_.asd);
     dt4 =  T0 * T1 * (1 + model_.asd) / 3
                          - cdbox * paramPtr->vsdth;
  }

  /* v2.2.2 bug fix */
  T0 = model_.csdesw * log(1 + model_.tsi /
     model_.tbox);
  T1 = sourcePerimeter - w;
  if (T1 > 0.0)
     csesw = T0 * T1;
  else
     csesw = 0.0;
  T1 = drainPerimeter - w;
  if (T1 > 0.0)
     cdesw = T0 * T1;
  else
     cdesw = 0.0;

  phi = paramPtr->phi;
  cgso = paramPtr->cgso;
  cgdo = paramPtr->cgdo;

/* v2.0 release */
  if (model_.ln < 1e-15) model_.ln = 1e-15;
  T0 = -0.5 * paramPtr->leff * paramPtr->leff /
    model_.ln / model_.ln;
  DEXP(T0,T1);
  paramPtr->arfabjt = T1;

  T0 = paramPtr->lbjt0 * (1.0 / paramPtr->leff + 1.0 / model_.ln);
  paramPtr->lratio = pow(T0,paramPtr->nbjt);
  paramPtr->lratiodif = 1.0 + model_.ldif0 * pow(T0,paramPtr->ndif);

  if ((paramPtr->vearly = paramPtr->vabjt + paramPtr->aely*paramPtr->leff) < 1)
     paramPtr->vearly = 1;

  /* vfbzb calculation for capMod 3 */
  tmp = sqrt(paramPtr->Xdep0);
  tmp1 = paramPtr->vbi - paramPtr->phi;
  tmp2 = model_.factor1 * tmp;

  T0 = -0.5 * paramPtr->dvt1w * paramPtr->weff
     * paramPtr->leff / tmp2;
  if (T0 > -EXPL_THRESHOLD)
  {   T1 = exp(T0);
      T2 = T1 * (1.0 + 2.0 * T1);
  }
  else
  {   T1 = MIN_EXPL;
      T2 = T1 * (1.0 + 2.0 * T1);
  }
  T0 = paramPtr->dvt0w * T2;
  T2 = T0 * tmp1;

  T0 = -0.5 * paramPtr->dvt1 * paramPtr->leff / tmp2;
  if (T0 > -EXPL_THRESHOLD)
  {   T1 = exp(T0);
      T3 = T1 * (1.0 + 2.0 * T1);
  }
  else
  {   T1 = MIN_EXPL;
      T3 = T1 * (1.0 + 2.0 * T1);
  }
  T3 = paramPtr->dvt0 * T3 * tmp1;

/* v2.2.3 */
  T4 = (model_.tox - model_.dtoxcv) * paramPtr->phi
     / (paramPtr->weff + paramPtr->w0);

  T0 = sqrt(1.0 + paramPtr->nlx / paramPtr->leff);
  T5 = paramPtr->k1eff * (T0 - 1.0) * paramPtr->sqrtPhi
     + (paramPtr->kt1 + paramPtr->kt1l / paramPtr->leff)
     * (TempRatio - 1.0);

  tmp3 = model_.dtype * paramPtr->vth0
       - T2 - T3 + paramPtr->k3 * T4 + T5;
  paramPtr->vfbzb = tmp3 - paramPtr->phi - paramPtr->k1eff
                     * paramPtr->sqrtPhi;
  /* End of vfbzb */


// v3.2
  paramPtr->qsi = CONSTQ * model_.npeak
                   * (1.0 + paramPtr->nlx / paramPtr->leff)
                   * 1e6 * model_.tsi;


// v3.1 wanh added for RF
  grgeltd = model_.rshg * (model_.xgw
            +paramPtr->weff/nseg / 3.0 / model_.ngcon) /
            (model_.ngcon * (l - model_.xgl));
  if (grgeltd > 0.0)
    grgeltd = 1.0 / grgeltd;
  else
  { grgeltd = 1.0e3; // mho
    if (rgateMod !=0)
    {
      string msg = "Warning: The gate conductance reset to 1.0e3 mho.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
    }
  }
// v3.1 wanh added for RF end

  paramPtr->ldeb = sqrt(CONSTEPSSI * model_.Vtm0 /
                        (CONSTQ * paramPtr->npeak * 1.0e6)) / 3.0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

  double VgstNVt, ExpVgst;

  double arg;

  double Vfbeff, dVfbeff_dVd, dVfbeff_dVg, dVfbeff_dVrg, dVfbeff_dVb;
  double dVfbeff_dT, V3, V4;
  double Vgs_eff, Vfb, dVfb_dVb, dVfb_dVd, dVfb_dT;

  double MJSWG;

  double qinoi;
  double Phis, dPhis_dVb, sqrtPhis, dsqrtPhis_dVb, Vth, dVth_dVb, dVth_dVd;
  double sqrtPhisExt, dsqrtPhisExt_dVb;
  double Vgst;

  double Vtm;
  double n, dn_dVb, dn_dVd, noff, dnoff_dVd, dnoff_dVb;
  double ExpArg, V0, CoxWLcen, QovCox, LINK;
  double CoxWLb, CoxWLcenb;
  double DeltaPhi;

  double Cox, Tox, Tcen, dTcen_dVg, dTcen_dVd, dTcen_dVb;
  double Ccen, Coxeff, dCoxeff_dVg, dCoxeff_dVd, dCoxeff_dVb;
  double dTcen_dT, dCoxeff_dT, dCoxWLcenb_dT;
  double Denomi ,dDenomi_dVg ,dDenomi_dVd ,dDenomi_dVb ,dDenomi_dT;
  double Ce1b ,Ce1e, Ce1T;
  double CbT, CsT, CgT;

  double Giie, dRatio_dVe;
  double Vdsatii0, dVdsatii0_dT;
  double VgsStep, dVgsStep_dT;
  double Vdsatii;
  double Vdiff, dVdiff_dVg, dVdiff_dVb, dVdiff_dVd, dVdiff_dVe, dVdiff_dT;
  double Ratio, dRatio_dVg, dRatio_dVb, dRatio_dVd, dRatio_dT;
  double Gbpbs, Gbpps;
  double qinv_local;
  double qjs_local, gcjsbs, gcjsT;
  double qjd_local, gcjdbs, gcjdds, gcjdT;
  double cjsbs, cjdbs, dcjdbs_dT, dcjsbs_dT;
  double PhiBSWG, dPhiBSWG_dT;
  double darg_dT, ddT3_dVb_dT;

  double gbbp;

  double OxideRatio, Vaux, dVaux_dVg, dVaux_dVd, dVaux_dVb;
  double Igb, dIgb_dVg, dIgb_dVd, dIgb_dVb, dIgb_dVe, dIgb_dT;
  double Ibs, Ibd;

  double Ibs1 ,dIbs1_dVb ,dIbs1_dT;
  double Ibs2 ,dIbs2_dVb ,dIbs2_dT;
  double Ibs3 ,dIbs3_dVb ,dIbs3_dVd, dIbs3_dT;
  double Ibs4 ,dIbs4_dVb ,dIbs4_dT;
  double Ibd1 ,dIbd1_dVb ,dIbd1_dVd ,dIbd1_dT;
  double Ibd2 ,dIbd2_dVb ,dIbd2_dVd ,dIbd2_dT;
  double Ibd3 ,dIbd3_dVb ,dIbd3_dVd ,dIbd3_dT;
  double Ibd4 ,dIbd4_dVb ,dIbd4_dVd ,dIbd4_dT;
  double Igc, dIgc_dVg, dIgc_dVd, dIgc_dVb;
  double Igs_local, dIgs_dVg, dIgs_dVs, Igd_local, dIgd_dVg, dIgd_dVd;
  double Igcs_local, dIgcs_dVg, dIgcs_dVd, dIgcs_dVb;
  double Igcd_local, dIgcd_dVg, dIgcd_dVd, dIgcd_dVb;
  double Igb1, dIgb1_dVg, dIgb1_dVd, dIgb1_dVb, dIgb1_dT, dIgb1_dVe;
  double Igb2, dIgb2_dVg, dIgb2_dVd, dIgb2_dVb, dIgb2_dT;

  double vgs_eff, dvgs_eff_dvg, vgd_eff, dvgd_eff_dvg;
  double dT0_dVox, Voxeff, dVoxeff_dVox, dVox_dT, dVaux_dT;

  double Vgb, dVgb_dVg, dVgb_dVb, Vox, dVox_dVg, dVox_dVd, dVox_dVb;
  double dT1_dVe, dT5_dVe, dVox_dVe, dVoxdepinv_dVe, dVaux_dVe;
  double Voxacc, dVoxacc_dVg, dVoxacc_dVd, dVoxacc_dVb;
  double Voxdepinv, dVoxdepinv_dVg, dVoxdepinv_dVb, dVoxdepinv_dVd;
  double dVoxdepinv_dT, VxNVt, ExpVxNVt;

  double Gjsd=0.0, Gjsb=0.0, GjsT=0.0, Gjdd=0.0, Gjdb=0.0, GjdT=0.0;
  double Ien, dIen_dT, Iendif, dIendif_dT;
  double Ibsdif, dIbsdif_dVb, dIbsdif_dT;
  double Ibddif, dIbddif_dVb, dIbddif_dVd, dIbddif_dT;
  double Ehlis, dEhlis_dVb, dEhlis_dT;
  double EhlisFactor, dEhlisFactor_dVb, dEhlisFactor_dT;
  double Ehlid, dEhlid_dVb, dEhlid_dVd, dEhlid_dT;
  double EhlidFactor, dEhlidFactor_dVb, dEhlidFactor_dVd, dEhlidFactor_dT;
  double E2ndFactor, dE2ndFactor_dVb, dE2ndFactor_dVd, dE2ndFactor_dT;

  double ExpVbsNVtm, dExpVbsNVtm_dVb, dExpVbsNVtm_dT;
  double ExpVbdNVtm, dExpVbdNVtm_dVb, dExpVbdNVtm_dVd, dExpVbdNVtm_dT;

  double ueff_local, dueff_dVg, dueff_dVd, dueff_dVb, dueff_dT;
  double Esat;

  double Vdsat;

  double EsatL, dEsatL_dVg, dEsatL_dVd, dEsatL_dVb, dEsatL_dT;

  double dVdsat_dVg, dVdsat_dVb, dVdsat_dVd, dVdsat_dT, Vasat;
  double dVasat_dVg, dVasat_dVb, dVasat_dVd, dVasat_dT;
  double Va, dVa_dVd, dVa_dVg, dVa_dVb, dVa_dT;
  double WTsi, NVtm1, NVtm2;
  double Ic;
  double dNVtm1_dT;
  double WsTsi, WdTsi;
  double NVtmf, NVtmr, dNVtmf_dT, dNVtmr_dT;

  double Ibp, Iii, Gcd, Gcb, GcT;
  double Giid=0.0, Giig=0.0, Giib=0.0, GiiT=0.0;

  double Vbseff, dVbseff_dVb;

  double CoxWL;
  double Idgidl, Gdgidld, Gdgidlg, Isgidl, Gsgidlg;

  double K1, Cbox, CboxWL, Vesfb;
  double T0, dT0_dVg, dT0_dVd, dT0_dVb, dT0_dVe;
  double T1, dT1_dVg, dT1_dVd, dT1_dVb;
  double T2, dT2_dVg, dT2_dVd, dT2_dVb, dT2_dVe;
  double T3, dT3_dVg, dT3_dVd, dT3_dVb, dT3_dT;
  double T4, dT4_dVd, dT4_dVb, dT4_dT;
  double T5, dT5_dVg, dT5_dVd, dT5_dVb, dT5_dT;
  double T6, dT6_dVg, dT6_dVd, dT6_dVb, dT6_dT;
  double T7;
  double T8, dT8_dVd;
  double T9, dT9_dVd;
  double T10, dT10_dVb, dT10_dVd;
  double T11, T12;
  double dT10_dT, dT11_dT, DioMax;

  double T13, T14;
  double dT11_dVb, dT13_dVb, dT14_dVb;

  double tmp, Abulk_local, dAbulk_dVb, Abulk0, dAbulk0_dVb;

  double VACLM, dVACLM_dVg, dVACLM_dVd, dVACLM_dVb, dVACLM_dT;
  double VADIBL, dVADIBL_dVg, dVADIBL_dVd, dVADIBL_dVb, dVADIBL_dT;

  double Xdep, dXdep_dVb, lt1, dlt1_dVb, ltw, dltw_dVb;
  double Delt_vth, dDelt_vth_dVb;

  double Theta0, dTheta0_dVb;
  double T3zb, lt1zb, ltwzb, Theta0zb;
  double Delt_vthzb, dDelt_vthzb_dT;
  double DeltVthwzb, dDeltVthwzb_dT;
  double DeltVthtempzb, dDeltVthtempzb_dT;
  double Vthzb, dVthzb_dT, Vfbzb, dVfbzb_dT;


  double Temp, TempRatio, dTempRatio_dT, tmp1, tmp2, tmp3, tmp4;

  // Temp-dependent values for self heating
  double phi_local;
  double ni_local, Eg, vbi, dvbi_dT, vfbb, dvfbb_dT, sqrtPhi, Xdep0;
  double Ahli, dAhli_dT, jbjt, jdif, jrec, djbjt_dT, djdif_dT, djrec_dT;
  double jtun, djtun_dT, u0temp, du0temp_dT, vsattemp, dvsattemp_dT;
  double rds0, drds0_dT, ua, ub, uc, dua_dT, dub_dT, duc_dT;
  double dni_dT, dT7_dT, dT0_dT, dT0_dT7, dT1_dT, dT1_dT7;
  double dT2_dT, dT2_dT7;

  double Vbs0, dVbs0_dVg, dVbs0_dVd, dVbs0_dVe, dVbs0_dT;
  double Vbs0mos, dVbs0mos_dVe, dVbs0mos_dT;
  double Vbsmos, dVbsmos_dVg, dVbsmos_dVd, dVbsmos_dVb, dVbsmos_dVe, dVbsmos_dT;

  double wdios, wdiod, wdiosCV, wdiodCV;

  double Vbp, dVbp_dVb, dVtm_dT;
  double Vbsh, dVbsh_dVb;
  double dDelt_vth_dT;
  double DeltVthw, dDeltVthw_dVb, dDeltVthw_dT, DeltVthtemp, dDeltVthtemp_dT;
  double VthFD, dVthFD_dVd, dVthFD_dVb, dVthFD_dVe, dVthFD_dT;
  double VtgsFD, ExpVtgsFD, VgstFD, ExpVgstFD;
  double VtgseffFD, dVtgseffFD_dVd, dVtgseffFD_dVg, dVtgseffFD_dVe;
  double dVtgseffFD_dT;
  double Vbsitf, dVbsitf_dVg, dVbsitf_dVd, dVbsitf_dVb, dVbsitf_dVe, dVbsitf_dT;
  double PhiFD, dPhiFD_dVg, dPhiFD_dVd, dPhiFD_dVe, dPhiFD_dT;
  double PhiON, dPhiON_dVg, dPhiON_dVd, dPhiON_dVe, dPhiON_dT;
  double dVbsh_dVg, dVbsh_dVd, dVbsh_dVe, dVbsh_dT;
  double VgsteffFD, dVgsteffFD_dVd, dVgsteffFD_dVg, dVgsteffFD_dVe;
  double dVgsteffFD_dT;
  double dVgsteff_dVe, dVbseff_dVg, dVbseff_dVd, dVbseff_dVe, dVbseff_dT;
  double Vbs0t, dVbs0t_dVg, dVbs0t_dVd, dVbs0t_dVe, dVbs0t_dT;
  double Vgsteff_local, dVgsteff_dVg, dVgsteff_dVd, dVgsteff_dVb, dVgsteff_dT;
  double Vdseff_local, dVdseff_dVg, dVdseff_dVd, dVdseff_dVb, dVdseff_dT;
  double dVgst_dVd, dVth_dT, dVgst2Vtm_dT;

  double DIBL_Sft, dDIBL_Sft_dVd, dDIBL_Sft_dVb;

  double Lambda, dLambda_dVg;

  double a1;

  double VdseffCV, dVdseffCV_dVg, dVdseffCV_dVd, dVdseffCV_dVb;
  double diffVds;

  double dAbulk_dVg;
  double beta, dbeta_dVg, dbeta_dVd, dbeta_dVb, dbeta_dT;
  double gche, dgche_dVg, dgche_dVd, dgche_dVb, dgche_dT;
  double fgche1, dfgche1_dVg, dfgche1_dVd, dfgche1_dVb, dfgche1_dT;
  double fgche2, dfgche2_dVg, dfgche2_dVd, dfgche2_dVb, dfgche2_dT;
  double Idl, dIdl_dVg, dIdl_dVd, dIdl_dVb, dIdl_dT;
  double Gm0, Gds0, Gmb0, GmT0;
  double Ids;

  double Gds, Gmb;

  double CoxWovL;
  double Rds, dRds_dVg, dRds_dVb, dRds_dT, WVCox, WVCoxRds;
  double Vgst2Vtm, VdsatCV;

//not used  double dVdsatCV_dVg, dVdsatCV_dVb;

  double Leff, Weff, dWeff_dVg, dWeff_dVb;
  double AbulkCV, dAbulkCV_dVb;

  double Qe1 , dQe1_dVb, dQe1_dVe, dQe1_dT;

  double Csg(0.0), Csd(0.0), Csb(0.0), Cbg(0.0), Cbd(0.0), Cbb(0.0);
  double Cgb(0.0), Cgg(0.0), Cgd(0.0);
  double Cgg1(0.0), Cgb1(0.0), Cgd1(0.0), Cbg1(0.0), Cbb1(0.0), Cbd1(0.0);
  double Csg1(0.0), Csd1(0.0), Csb1(0.0);

  double Qac0, Qsub0;
  double dQac0_dVg, dQac0_dVb, dQac0_dVd, dQac0_dVrg, dQac0_dT;
  double dQsub0_dVg, dQsub0_dVd, dQsub0_dVb, dQsub0_dVrg, dQsub0_dT;

  // These are needed for voltage limiting, no need for them to be instance
  // variables:
  double vg_old, vd_old, vp_old, ve_old, vgp_old, vgm_old, vs_old, delTemp_old;
  double vbd_old, vbs_old, vds_old;


  // Don't do charge computations in DC sweeps.
  if (solState.tranopFlag || solState.acopFlag || solState.transientFlag)
  {
    ChargeComputationNeeded = true;
  }
  else
  {
    ChargeComputationNeeded = false;
  }

  // this block of variables were originally set up as local function variables
  // but they've been moved to the instance class.
  qgdo = 0.0;
  qgso = 0.0;
  qgd = 0.0;
  qgs = 0.0;
  qge = 0.0;
  qgme = 0.0;
  qgate = 0.0;
  qbody = 0.0;
  qdrn = 0.0;
  qsub = 0.0;
  qsrc = 0.0;

  // end of b3ld.c parameters.

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  In updateIntermediateVars\n";
    cout << "  name = " << getName();
    cout << "  model name = " << model_.getName();
    cout <<"   dtype is " << model_.dtype << endl;
    cout.width(21); cout.precision(13); cout.setf(ios::scientific);
    cout << "  " << endl;
  }
#endif

  int Check = 0; // The limiter function sets this to 1 if it changes things

  // The first block of code in b3ld.c basically sets up, locally,
  // what the load function should use as values for the various solution
  // variables.  There is a series of IF statements which are dependent
  // upon the mode.  (transient, initializing transient, operating point,
  // small signal, etc.).  Xyce treats the operating point and transient
  // calculation in the same way, from the device's point of view, and
  // we don't support any of the other modes.  Therefore most of these
  // mode options are not here - only the transient mode stuff.

  // First get some of the needed solution variables:
  Vd     = 0.0;
  Vg     = 0.0;
  Vs     = 0.0;
  Ve     = 0.0;
  Vp     = 0.0;
  Vb     = 0.0;
  delTemp = 0.0;
  Vsp    = 0.0;
  Vdp    = 0.0;
  Vgp    = 0.0;
  Vgm    = 0.0;

  Vd = (extData.nextSolVectorRawPtr)[li_Drain];
  // possible source of confusion: Original SPICE b3soi calls this "GE" for
  // "Gate External", we call it Gate
  Vg = (extData.nextSolVectorRawPtr)[li_Gate];
  Vs = (extData.nextSolVectorRawPtr)[li_Source];
  Ve = (extData.nextSolVectorRawPtr)[li_Substrate];

  if(li_Body != -1)  {
     Vb = (extData.nextSolVectorRawPtr)[li_Body];
  }
  if(li_ExtBody != -1)
     Vp = (extData.nextSolVectorRawPtr)[li_ExtBody];
  if(li_Temperature != -1)
     delTemp = (extData.nextSolVectorRawPtr)[li_Temperature];
  Vsp = (extData.nextSolVectorRawPtr)[li_SourcePrime];
  Vdp = (extData.nextSolVectorRawPtr)[li_DrainPrime];
  // Possible source of confusion: Original Spice calls this Gate, we call
  // it GatePrime to be consistent with our usage in other devices
  // (Prime means internal, unadorned means external)
  Vgp = (extData.nextSolVectorRawPtr)[li_GatePrime];
  Vgm = (extData.nextSolVectorRawPtr)[li_GateMid];

  // orig variables, nodal version:
  Vd_orig = Vd;
  Vg_orig = Vg;
  Vs_orig = Vs;
  Ve_orig = Ve;
  Vb_orig = Vb;
  Vp_orig = Vp;
  Vsp_orig = Vsp;
  Vdp_orig = Vdp;
  Vgp_orig = Vgp;
  Vgm_orig = Vgm;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout.precision(18);
    cout << endl;
    cout << getName() << "  Blim:     Vd = " << Vd << endl;
    cout << getName() << "  Blim:     Vg = " << Vg << endl;
    cout << getName() << "  Blim:     Vs = " << Vs << endl;
    cout << getName() << "  Blim:     Ve = " << Ve << endl;
    cout << getName() << "  Blim:     Vb = " << Vb << endl;
    cout << getName() << "  Blim:     Vp = " << Vp << endl;
    cout << getName() << "  Blim:delTemp = " << delTemp << endl;
    cout << getName() << "  Blim:    Vsp = " << Vsp << endl;
    cout << getName() << "  Blim:    Vdp = " << Vdp << endl;
    cout << getName() << "  Blim:    Vgp = " << Vgp << endl;
    cout << getName() << "  Blim:    Vgm = " << Vgm << endl;
    cout <<  endl;
  }
#endif

  // modified from b3ld:  (see lines 221-230)
  vbs  = model_.dtype * (Vb - Vsp);
  vps  = model_.dtype * (Vp - Vsp);
  vgs  = model_.dtype * (Vgp- Vsp);
  ves  = model_.dtype * (Ve - Vsp);
  vds  = model_.dtype * (Vdp- Vsp);
  vges = model_.dtype * (Vg - Vsp);
  vgms = model_.dtype * (Vgm - Vsp);

  vbd  = vbs - vds;
  vgd  = vgs - vds;
  ved  = ves - vds;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout.precision(18);
    cout << endl;
    cout << getName() << "  Blim:    vbs = " << vbs << endl;
    cout << getName() << "  Blim:    vbd = " << vbd << endl;
    cout << getName() << "  Blim:    vps = " << vps << endl;
    cout << getName() << "  Blim:    vpd = " << vpd << endl;
    cout << getName() << "  Blim:    vgs = " << vgs << endl;
    cout << getName() << "  Blim:    vds = " << vds << endl;
    cout << getName() << "  Blim:    ves = " << ves << endl;
    cout << getName() << "  Blim:    ved = " << ved << endl;
    cout << getName() << "  Blim:    vgd = " << vgd << endl;
    cout << getName() << "  Blim:   vges = " << vges<< endl;
    cout << getName() << "  Blim:   vgms = " << vgms<< endl;
  }
#endif

  origFlag = 1;
  vbs_orig = vbs;
  vbd_orig = vbd;
  vps_orig = vps;
  vpd_orig = vpd;
  vgs_orig = vgs;
  vds_orig = vds;
  ves_orig = ves;
  ved_orig = ved;
  vgd_orig = vgd;
  vges_orig = vges;
  vgms_orig = vgms;
  delTemp_orig = delTemp;

  // note initJctFlag will only be true for dcop.
  if (solState.initJctFlag && !OFF && devOptions.voltageLimiterFlag)
  {
    if (solState.inputOPFlag)
    {
      N_LAS_Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
      if ((*flagSolVectorPtr)[li_Drain] == 0 || (*flagSolVectorPtr)[li_Gate] == 0 ||
          (*flagSolVectorPtr)[li_Source] == 0 || (*flagSolVectorPtr)[li_Substrate] == 0 ||
          (*flagSolVectorPtr)[li_DrainPrime] == 0 || (*flagSolVectorPtr)[li_GatePrime] == 0 ||
          (*flagSolVectorPtr)[li_SourcePrime] == 0 || (*flagSolVectorPtr)[li_GateMid] == 0 ||
          (li_Body != -1 && (*flagSolVectorPtr)[li_Body] == 0) ||
          (li_ExtBody != -1 && (*flagSolVectorPtr)[li_ExtBody] == 0) ||
          (li_Temperature != -1 && (*flagSolVectorPtr)[li_Temperature] == 0) )
      {
        vbs = 0.0;
        vgs = model_.dtype*0.1 + paramPtr->vth0;
        vds = 0.0;
        ves = 0.0;
        vps = 0.0;
        vpd = 0.0;
        vges = vgms = vgs;
        origFlag = 0;
      }
    }
    else
    {
      vgs = model_.dtype*0.1 + paramPtr->vth0;
      vges = vgms = vgs;
      origFlag = 0;
    }
  }
  else if ((solState.initFixFlag || solState.initJctFlag) && OFF)
  {
    delTemp = vps = vbs = vgs = vds = ves = 0.0;
    Vg = Vd = Vs = Vp = Ve = 0.0;

    vges = vgms = 0.0;
  }

  if (solState.newtonIter == 0)
  {
    newtonIterOld = 0;

    if (!solState.dcopFlag || (solState.locaEnabledFlag && solState.dcopFlag))
    // ie, first newton step of a transient time step or DCOP continuation step.
    {
      // if not dcop, then state vector has final drops of last nonlinear step
      // of most recent successful time step.  Use those as our "old" values
      // for limiting
      vg_old = (extData.currStoVectorRawPtr)[li_store_vg];
      vd_old = (extData.currStoVectorRawPtr)[li_store_vd];
      vs_old = (extData.currStoVectorRawPtr)[li_store_vs];
      vp_old = (extData.currStoVectorRawPtr)[li_store_vp];
      ve_old = (extData.currStoVectorRawPtr)[li_store_ve];
      vgp_old = (extData.currStoVectorRawPtr)[li_store_vgp];
      vgm_old = (extData.currStoVectorRawPtr)[li_store_vgm];
      delTemp_old = (extData.currStoVectorRawPtr)[li_store_deltemp];

      // old voltage drops  needed for second round of limiting
      vds_old = (extData.currStoVectorRawPtr)[li_store_vds];
      vbs_old = (extData.currStoVectorRawPtr)[li_store_vbs];
      vbd_old = (extData.currStoVectorRawPtr)[li_store_vbd];
    }
    else // first newton step of DCOP.
    {
      vg_old = Vg;
      vd_old = Vdp;
      vs_old = Vsp;
      vp_old = Vp;
      ve_old = Ve;
      vgp_old = Vgp;
      vgm_old = Vgm;
      delTemp_old = delTemp;

      // old voltage drops
      vds_old = vds;
      vbs_old = vbs;
      vbd_old = vbd;
    }
  }
  else
  {
    vg_old = (extData.nextStoVectorRawPtr)[li_store_vg];
    vd_old = (extData.nextStoVectorRawPtr)[li_store_vd];
    vs_old = (extData.nextStoVectorRawPtr)[li_store_vs];
    vp_old = (extData.nextStoVectorRawPtr)[li_store_vp];
    ve_old = (extData.nextStoVectorRawPtr)[li_store_ve];
    vgp_old = (extData.nextStoVectorRawPtr)[li_store_vgp];
    vgm_old = (extData.nextStoVectorRawPtr)[li_store_vgm];
    delTemp_old = (extData.nextStoVectorRawPtr)[li_store_deltemp];

    // old voltage drops  needed for second round of limiting
    vds_old = (extData.nextStoVectorRawPtr)[li_store_vds];
    vbs_old = (extData.nextStoVectorRawPtr)[li_store_vbs];
    vbd_old = (extData.nextStoVectorRawPtr)[li_store_vbd];
  }

  // B3SOI does limiting in a very different manner than any of the other
  // SPICE/Xyce devices did.  Instead of limiting junction voltages and
  // preventing them from swinging too much, B3SOI limits actual voltage
  // nodes directly.
  // The impact of this on Xyce voltage limiting is significant:
  // The most straightforward way to apply Xyce style voltage limiting is
  // to apply the jacobian directly to the delta-x resulting from the
  // voltage node changes.  But since most of the B3SOI code that could work
  // with limiting was written to use voltage drop changes instead, we're
  // gonna try to cheat.
  // We'll limit the individual nodes as B3SOI does, then recalculate the
  // drops.  We have already saved the drops from the input solution vector
  // in the _orig variables.  So now everywhere that spice does something
  // like:
  // RightHandSideTerm = current-conductance*voltagedrop
  // we can do
  // RightHandSideTerm = current;
  // RightHandSideTerm_Jdxp = conductance*(voltagedrop-voltagedrop_orig);

  if (devOptions.voltageLimiterFlag && !(solState.initJctFlag)
      && !(solState.initFixFlag && OFF))
  {
    Vg  = B3SOIlimit(Vg,  vg_old,  3.0, &Check);
    Vdp = B3SOIlimit(Vdp, vd_old,  3.0, &Check);
    Vsp = B3SOIlimit(Vsp, vs_old,  3.0, &Check);
    Vp  = B3SOIlimit(Vp,  vp_old,  3.0, &Check);
    Ve  = B3SOIlimit(Ve,  ve_old,  3.0, &Check);
    Vgp = B3SOIlimit(Vgp, vgp_old, 3.0, &Check);
    Vgm = B3SOIlimit(Vgm, vgm_old, 3.0, &Check);
    if (Check == 1)  // something changed, recalculate the drops
    {

      // Note that vbs is NOT recalculated!  In SPICE b3soi vbs is calculated
      // from the unlimited b and sp nodes
      vps  = model_.dtype * (Vp - Vsp);
      vgs  = model_.dtype * (Vgp- Vsp);
      ves  = model_.dtype * (Ve - Vsp);
      vds  = model_.dtype * (Vdp- Vsp);
      vges = model_.dtype * (Vg - Vsp);
      vgms = model_.dtype * (Vgm - Vsp);

      vbd  = vbs - vds;
      vgd  = vgs - vds;
      ved  = ves - vds;
      origFlag = 0;
    }

  } // devOptions.voltageLimiterFlag

  // ALL Bypass code removed by TVR --- we never, ever use it, and it just
  // clutters the code

  // There is now a second round of limiting-related stuff in the SPICE B3SOI:
  // In spice3F5 version there are additional lines here but they calculate
  // junk that is only used in bypass, which we never do

  if (devOptions.voltageLimiterFlag && !(solState.initJctFlag))
  {
    if (vds_old >= 0)  // normal mode
    {
      T0 = vbs_old;

      if ( model_.bug1830fix != 0)
      {
        T1 = vbd_old;
      }

    }
    else  // reverse
    {
      T0 = vbd_old;

      if ( model_.bug1830fix != 0)
      {
        T1 = vbs_old;
      }

    }

    if (vds >= 0)
    {
      vbs = B3SOIlimit(vbs, T0, 0.2, &Check);
      vbd = vbs - vds;

      if ( model_.bug1830fix != 0)
      {
        //vbd = B3SOIlimit(vbd, T1, 0.2, &Check);
        vbd = B3SOIlimit(vbd, vbd_old, 0.2, &Check);
      }

      Vb = model_.dtype*vbs+Vsp;
    }
    else
    {
      vbd = B3SOIlimit(vbd,T0,0.2,&Check);
      vbs = vbd+vds;

      if ( model_.bug1830fix != 0)
      {
        //vbs = B3SOIlimit(vbs, T1, 0.2, &Check);
        vbs = B3SOIlimit(vbs, vbs_old, 0.2, &Check);
      }

      Vb = model_.dtype*vbs+Vdp;
    }

    delTemp = B3SOIlimit(delTemp, delTemp_old, 5.0, &Check);

    if (Check == 1) origFlag = 0;

#ifdef Xyce_DEBUG_DEVICE
    if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
    {
      cout << endl;
      cout.precision(18);
      cout << endl;
      cout << getName() << "  Alim:     Vd = " << Vd << endl;
      cout << getName() << "  Alim:     Vg = " << Vg << endl;
      cout << getName() << "  Alim:     Vs = " << Vs << endl;
      cout << getName() << "  Alim:     Ve = " << Ve << endl;
      cout << getName() << "  Alim:     Vb = " << Vb << endl;
      cout << getName() << "  Alim:     Vp = " << Vp << endl;
      cout << getName() << "  Alim:delTemp = " << delTemp << endl;
      cout << getName() << "  Alim:    Vsp = " << Vsp << endl;
      cout << getName() << "  Alim:    Vdp = " << Vdp << endl;
      cout << getName() << "  Alim:    Vgp = " << Vgp << endl;
      cout << getName() << "  Alim:    Vgm = " << Vgm << endl;

      cout << getName() << "  Alim:    vbs = " << vbs << endl;
      cout << getName() << "  Alim:    vps = " << vps << endl;
      cout << getName() << "  Alim:    vgs = " << vgs << endl;
      cout << getName() << "  Alim:    ves = " << ves << endl;
      cout << getName() << "  Alim:    vds = " << vds << endl;
      cout << getName() << "  Alim:   vges = " <<vges << endl;
      cout << getName() << "  Alim:   vgms = " <<vgms << endl;

      cout << getName() << "  Alim:    vbd = " << vbd << endl;

      cout << getName() << "  Alim:     T0 = " << vbd << endl;
      cout <<  endl;
      cout <<  endl;

      cout <<getName()<<"  Vg_orig      = " << Vg_orig  << endl;
      cout <<getName()<<"  Vs_orig      = " << Vs_orig  << endl;
      cout <<getName()<<"  Ve_orig      = " << Ve_orig  << endl;
      cout <<getName()<<"  Vb_orig      = " << Vb_orig  << endl;
      cout <<getName()<<"  Vp_orig      = " << Vp_orig  << endl;
      cout <<getName()<<"  Vsp_orig     = " << Vsp_orig  << endl;
      cout <<getName()<<"  Vdp_orig     = " << Vdp_orig  << endl;
      cout <<getName()<<"  Vgp_orig     = " << Vgp_orig  << endl;
      cout <<getName()<<"  Vgm_orig     = " << Vgm_orig  << endl;
      cout << endl;
      cout <<getName()<<"  vbs_orig     = " << vbs_orig  << endl;
      cout <<getName()<<"  vbd_orig     = " << vbd_orig  << endl;
      cout <<getName()<<"  vps_orig     = " << vps_orig  << endl;
      cout <<getName()<<"  vpd_orig     = " << vpd_orig  << endl;
      cout <<getName()<<"  vgs_orig     = " << vgs_orig  << endl;
      cout <<getName()<<"  vds_orig     = " << vds_orig  << endl;
      cout <<getName()<<"  ves_orig     = " << ves_orig  << endl;
      cout <<getName()<<"  ved_orig     = " << ved_orig  << endl;
      cout <<getName()<<"  vgd_orig     = " << vgd_orig  << endl;
      cout <<getName()<<"  vges_orig    = " << vges_orig  << endl;
      cout <<getName()<<"  vgms_orig    = " << vgms_orig  << endl;
      cout <<getName()<<"  delTemp_orig = " << delTemp_orig  << endl;
    }
#endif
  }

  // update the "old" newton iteration number.
  if (solState.newtonIter != 0 && solState.newtonIter != newtonIterOld)
  {
    newtonIterOld = solState.newtonIter;
  }

  // Finished with what would have been the series of CKTmode
  // IF statements...

  // Calculate temperature dependent values for self-heating effect
  Temp = delTemp + temp;
  dTempRatio_dT = 1 / model_.tnom;
  TempRatio = Temp * dTempRatio_dT;

  if (selfheat)
  {
    Vtm = CONSTKoverQ * Temp;

    T0 = 1108.0 + Temp;
    T5 = Temp * Temp;
    Eg = 1.16 - 7.02e-4 * T5 / T0;
    T1 = ((7.02e-4 * T5) - T0 * (14.04e-4 * Temp)) / T0 / T0;
    //  T1 = dEg / dT

    T2 = 1.9230584e-4;
    //  T2 = 1 / 300.15^(3/2)
    T5 = sqrt(Temp);
    T3 = 1.45e10 * Temp * T5 * T2;
    T4 = exp(21.5565981 - Eg / (2.0 * Vtm));
    ni_local = T3 * T4;
    dni_dT = 2.175e10 * T2 * T5 * T4 + T3 * T4 *
                        (-Vtm * T1 + Eg * CONSTKoverQ) / (2.0 * Vtm * Vtm);

    T0 = log(1.0e20 * paramPtr->npeak / (ni_local * ni_local));
    vbi = Vtm * T0;
    dvbi_dT = CONSTKoverQ * T0 + Vtm * (-2.0 * dni_dT / ni_local);

    if (paramPtr->nsub > 0)
    {
      T0 = log(paramPtr->npeak / paramPtr->nsub);
      vfbb = -model_.dtype * Vtm * T0;
      dvfbb_dT = -model_.dtype * CONSTKoverQ * T0;
    }
    else
    {
      T0 = log(-paramPtr->npeak * paramPtr->nsub / ni_local / ni_local);
      vfbb = -model_.dtype * Vtm * T0;
      dvfbb_dT = -model_.dtype *
        (CONSTKoverQ * T0 - Vtm * 2.0 * dni_dT / ni_local);
    }

    //    phi_local = 2.0 * Vtm * log(paramPtr->npeak / ni_local);
    phi_local = phi;
    sqrtPhi = sqrt(phi_local);
    Xdep0 = sqrt(2.0 * CONSTEPSSI/(CONSTQ * paramPtr->npeak * 1.0e6))
          * sqrtPhi;

    //  Save the values below for phi_local calculation in B3SOIaccept()
    vtm = Vtm;
    ni = ni_local;

    T3 = TempRatio - 1.0;
    T8 = 1/ model_.tnom;
    T4 = CONSTEg300 / Vtm * T3;
    dT4_dT = CONSTEg300 / Vtm / Vtm * (Vtm * T8 - T3 * CONSTKoverQ);

    T7 = paramPtr->xbjt * T4 / paramPtr->ndiode;
    dT7_dT = paramPtr->xbjt * dT4_dT / paramPtr->ndiode;
    CEXP(T7, T0, dT0_dT7);
    dT0_dT = dT0_dT7 * dT7_dT;

    if (paramPtr->xbjt == paramPtr->xdif)
    {
      T1 = T0;
      dT1_dT = dT0_dT;
    }
    else
    {
      T7 = paramPtr->xdif * T4 / paramPtr->ndiode;
      dT7_dT = paramPtr->xdif * dT4_dT / paramPtr->ndiode;
      CEXP(T7, T1, dT1_dT7);
      dT1_dT = dT1_dT7 * dT7_dT;
    }

    T7 = paramPtr->xrec * T4 / paramPtr->nrecf0;
    dT7_dT = paramPtr->xrec * dT4_dT / paramPtr->nrecf0;
    CEXP(T7, T2, dT2_dT7);
    dT2_dT = dT2_dT7 * dT7_dT;

    /* high level injection */
    Ahli = paramPtr->ahli * T0;
    dAhli_dT = paramPtr->ahli * dT0_dT;

    jbjt = paramPtr->isbjt * T0;
    jdif = paramPtr->isdif * T1;
    jrec = paramPtr->isrec * T2;
    djbjt_dT = paramPtr->isbjt * dT0_dT;
    djdif_dT = paramPtr->isdif * dT1_dT;
    djrec_dT = paramPtr->isrec * dT2_dT;

    T7 = paramPtr->xtun * T3;
    dT7_dT = paramPtr->xtun * T8;
    CEXP(T7, T0, dT0_dT7);
    dT0_dT = dT0_dT7 * dT7_dT;
    jtun = paramPtr->istun * T0;
    djtun_dT = paramPtr->istun * dT0_dT;

    u0temp = paramPtr->u0 * pow(TempRatio, paramPtr->ute);
    du0temp_dT = paramPtr->u0 * paramPtr->ute *
                            pow(TempRatio, paramPtr->ute - 1.0) * T8;

    vsattemp = paramPtr->vsat - paramPtr->at * T3;
    dvsattemp_dT = -paramPtr->at * T8;

    rds0 = (paramPtr->rdsw + paramPtr->prt * T3) / paramPtr->rds0denom;
    drds0_dT = paramPtr->prt / paramPtr->rds0denom * T8;

    ua = paramPtr->uatemp + paramPtr->ua1 * T3;
    ub = paramPtr->ubtemp + paramPtr->ub1 * T3;
    uc = paramPtr->uctemp + paramPtr->uc1 * T3;
    dua_dT = paramPtr->ua1 * T8;
    dub_dT = paramPtr->ub1 * T8;
    duc_dT = paramPtr->uc1 * T8;
  }
  else
  {
    Vtm = CONSTKoverQ * temp;
    vbi = paramPtr->vbi;
    vfbb = paramPtr->vfbb;
    phi_local = paramPtr->phi;
    sqrtPhi = paramPtr->sqrtPhi;
    Xdep0 = paramPtr->Xdep0;
    jbjt = paramPtr->jbjt;
    jdif = paramPtr->jdif;
    jrec = paramPtr->jrec;
    jtun = paramPtr->jtun;

    // v2.2.2 bug fix
    Ahli = paramPtr->ahli0;

    u0temp = paramPtr->u0temp;
    vsattemp = paramPtr->vsattemp;
    rds0 = paramPtr->rds0;
    ua = paramPtr->ua;
    ub = paramPtr->ub;
    uc = paramPtr->uc;
    dni_dT = dvbi_dT = dvfbb_dT = djbjt_dT = djdif_dT = 0.0;
    djrec_dT = djtun_dT = du0temp_dT = dvsattemp_dT = 0.0;
    drds0_dT = dua_dT = dub_dT = duc_dT = 0.0;
    dAhli_dT = 0;
  }

  // TempRatio used for Vth and mobility
  if (selfheat)
  {
    TempRatioMinus1 = Temp / model_.tnom - 1.0;
  }
  else
  {
    TempRatioMinus1 =  temp / model_.tnom - 1.0;
  }

  // file: b3ld.c   line: 937
  // determine DC current and derivatives
  vbd = vbs - vds;
  vgd = vgs - vds;
  vgb = vgs - vbs;
  ved = ves - vds;
  veb = ves - vbs;
  vge = vgs - ves;
  vpd = vps - vds;

  vged = vges - vds;
  vgmd = vgms - vds;
  vgme = vgms - ves;
  vgmb = vgms - vbs;    // v3.2 bug fix

// v3.1 bug fix
  wdiosCV_NoSwap = paramPtr->wdiosCV;
  wdiodCV_NoSwap = paramPtr->wdiodCV;

  // This stuff impacts voltage limiting!  Be careful!
  // It changes the meanings of things under the covers, so we must be
  // consistent
  if (vds >= 0.0)
  {   // normal mode
      mode = 1;
      Vds = vds;
      Vds_orig = vds_orig;
      Vgs = vgs;
      Vgs_orig = vgs_orig;
      Vbs = vbs;
      Vbs_orig = vbs_orig;
      Vbd = vbd;
      Vbd_orig = vbd_orig;
      Ves = ves;
      Ves_orig = ves_orig;
      Vps = vps;
      Vps_orig = Vps_orig;

      wdios = paramPtr->wdios;
      wdiod = paramPtr->wdiod;
      wdiosCV = paramPtr->wdiosCV;
      wdiodCV = paramPtr->wdiodCV;
  }
  else
  {   // inverse mode
      mode = -1;
      Vds = -vds;
      Vds_orig = -vds_orig;
      Vgs = vgd;
      Vgs_orig = vgd_orig;
      Vbs = vbd;
      Vbs_orig = vbd_orig;
      Vbd = vbs;
      Vbd_orig = vbs_orig;
      Ves = ved;
      Ves_orig = ved_orig;
      Vps = vpd;
      Vps_orig = vpd_orig;

      wdios = paramPtr->wdiod;
      wdiod = paramPtr->wdios;
      wdiosCV = paramPtr->wdiodCV;
      wdiodCV = paramPtr->wdiosCV;
  }

 // mosfet continuation.
  // This idea is based, loosely, on a paper by Jaijeet
  // Rosychowdhury.  If the artificial parameter flag has been enabled,
  // modify Vds and Vgs.
#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout << "HOMOTOPY INFO: gainscale   = "
         << solState.gainScale[blockHomotopyID] << endl;
    cout << "HOMOTOPY INFO: before vds  = " << Vds << endl;
    cout << "HOMOTOPY INFO: before vgst = " << Vgs << endl;
    cout << "vgstConst= " << devOptions.vgstConst << endl;
  }
#endif
  if (solState.artParameterFlag)
  {
    double alpha = solState.gainScale[blockHomotopyID];
    if (devOptions.staggerGainScale)
    {
      alpha *= (0.3 * randomPerturb + 1.0);
      if (alpha > 1.0)
      {
        alpha = 1.0;
      }
    }
    double vgstConst = devOptions.vgstConst;
    if (devOptions.randomizeVgstConst)
    {
      vgstConst *= randomPerturb;
    }

    Vds = devSupport.contVds (Vds, solState.nltermScale,
                              devOptions.vdsScaleMin);
    Vgs = devSupport.contVgst(Vgs, alpha, vgstConst);
  }
#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout << "HOMOTOPY INFO: after vds   = " << Vds << endl;
    cout << "HOMOTOPY INFO: after vgst  = " << Vgs << endl;
  }
#endif
  // end of mosfet continuation block.

  Vesfb = Ves - vfbb;
  Cbox  = model_.cbox;
  K1    = paramPtr->k1eff;

// Poly Gate Si Depletion Effect
  T0 = paramPtr->vfb + phi_local;
  if ((paramPtr->ngate > 1.e18) && (paramPtr->ngate < 1.e25) && (Vgs > T0))
      // added to avoid the problem caused by ngate
  {   T1 = 1.0e6 * CONSTQ * CONSTEPSSI * paramPtr->ngate
         / (model_.cox * model_.cox);
      T4 = sqrt(1.0 + 2.0 * (Vgs - T0) / T1);
      T2 = T1 * (T4 - 1.0);
      T3 = 0.5 * T2 * T2 / T1;   // T3 = Vpoly
      T7 = 1.12 - T3 - 0.05;
      T6 = sqrt(T7 * T7 + 0.224);
      T5 = 1.12 - 0.5 * (T7 + T6);
      Vgs_eff = Vgs - T5;
      dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
  }
  else
  {   Vgs_eff = Vgs;
      dVgs_eff_dVg = 1.0;
  }

  Leff = paramPtr->leff;

  if (selfheat)
  {
      Vtm = CONSTKoverQ * Temp;
      dVtm_dT = CONSTKoverQ;
  }
  else {
      dVtm_dT = 0.0;
  }

  V0 = vbi - phi_local;


// begin of v3.0 block addition
// B/S built-in potential lowering calculation
  if (soiMod == 0) // BSIMPD -  v3.2
  {
     Vbsmos = Vbs;
     dVbsmos_dVg = 0.0;
     dVbsmos_dVd = 0.0;
     dVbsmos_dVb = 1.0;
     dVbsmos_dVe = 0.0;
     if (selfheat)  dVbsmos_dT = 0.0;
     else  dVbsmos_dT = 0.0;

     Vbp = Vbs - Vps;
     dVbp_dVb = 1;
  }
  else   // soiMod = 1 or 2: adding FD module on top of BSIMPD
  {
     // prepare Vbs0 & Vbs0mos for VthFD calculation
     T0 = -model_.dvbd1 * paramPtr->leff / paramPtr->litl;
     T1 = model_.dvbd0 * (exp(0.5*T0) + 2*exp(T0));
     T2 = T1 * (vbi - phi_local);
     T3 = 0.5 * paramPtr->qsi / model_.csi; // v3.2
     Vbs0t = phi_local - T3 + model_.vbsa + T2;
     if (selfheat)
        dVbs0t_dT = T1 * dvbi_dT;
     else
        dVbs0t_dT = 0.0;

     T0 = 1 + model_.csi / Cbox;
     T3 = -model_.dk2b * paramPtr->leff / paramPtr->litl;
     T5 = model_.k2b * (exp(0.5*T3) + 2*exp(T3));
     T1 = (model_.k1b - T5) / T0;
     T2 = T1 * Vesfb;
     T4 = 1.0/(1 + Cbox / model_.csi);
     Vbs0 = T4 * Vbs0t + T2;
     dVbs0_dVe = T1;
     if (selfheat)
        dVbs0_dT = T4 * dVbs0t_dT - T1 * dvfbb_dT;
     else
        dVbs0_dT = 0.0;


     // zero field body potential calc
     T1 = Vbs0t - Vbs0 - 0.005;
     T2 = sqrt(T1 * T1 + (2.5e-5));
     T3 = 0.5 * (T1 + T2);
     T4 = T3 * model_.csi / paramPtr->qsi; // v3.2
     Vbs0mos = Vbs0 - 0.5 * T3 * T4;
     T5 = 0.5 * T4 * (1 + T1 / T2);
     dVbs0mos_dVe = dVbs0_dVe * (1 + T5);
     if (selfheat)
        dVbs0mos_dT = dVbs0_dT * (1 + T5) - T5 * dVbs0t_dT;
     else
        dVbs0mos_dT = 0.0;


     // set the upperbound of Vbs0mos to be phi_local for sqrt calc.
     T1 = phi_local - 0.02;
     T2 = T1 - Vbs0mos - 0.005;
     T3 = sqrt(T2 * T2 + 4.0 * 0.005);
     Vbs0mos = T1 - 0.5 * (T2 + T3);
     T4 = 0.5 * (1 + T2 / T3);
     dVbs0mos_dVe = T4 * dVbs0mos_dVe;
     if (selfheat)
        dVbs0mos_dT = T4 * dVbs0mos_dT;
     else  dVbs0mos_dT = 0.0;


     // VthFD calculation
     Phis = phi_local - Vbs0mos;
     dPhis_dVb = -1;   // w.r.t Vbs0mos
     sqrtPhis = sqrt(Phis);
     dsqrtPhis_dVb = -0.5 / sqrtPhis;
     Xdep = Xdep0 * sqrtPhis / sqrtPhi;
     dXdep_dVb = (Xdep0 / sqrtPhi) * dsqrtPhis_dVb;
     T3 = sqrt(Xdep);

     T0 = paramPtr->dvt2 * Vbs0mos;
     if (T0 >= - 0.5)
     {   T1 = 1.0 + T0;
         T2 = paramPtr->dvt2 ;
     }
     else   // Added to avoid any discontinuity problems caused by dvt2
     {   T4 = 1.0 / (3.0 + 8.0 * T0);
         T1 = (1.0 + 3.0 * T0) * T4;
         T2 = paramPtr->dvt2 * T4 * T4 ;
     }
     lt1 = model_.factor1 * T3 * T1;
     dlt1_dVb =model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

     T0 = paramPtr->dvt2w * Vbs0mos;
     if (T0 >= - 0.5)
     {   T1 = 1.0 + T0;
         T2 = paramPtr->dvt2w ;
     }
     else    // Added to avoid any discontinuity problems caused by dvt2w
     {   T4 = 1.0 / (3.0 + 8.0 * T0);
         T1 = (1.0 + 3.0 * T0) * T4;
         T2 = paramPtr->dvt2w * T4 * T4 ;
     }
     ltw= model_.factor1 * T3 * T1;
     dltw_dVb=model_.factor1*(0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

     T0 = -0.5 * paramPtr->dvt1 * Leff / lt1;
     if (T0 > -EXPL_THRESHOLD)
     {   T1 = exp(T0);
         Theta0 = T1 * (1.0 + 2.0 * T1);
         dT1_dVb = -T0 / lt1 * T1 * dlt1_dVb;
         dTheta0_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
     }
     else
     {   T1 = MIN_EXPL;
         Theta0 = T1 * (1.0 + 2.0 * T1);
         dTheta0_dVb = 0.0;
     }

     thetavth = paramPtr->dvt0 * Theta0;
     Delt_vth = thetavth * V0;
     dDelt_vth_dVb = paramPtr->dvt0 * dTheta0_dVb * V0;
     if (selfheat)  dDelt_vth_dT = thetavth * dvbi_dT;
     else  dDelt_vth_dT = 0.0;

     T0 = -0.5 * paramPtr->dvt1w * paramPtr->weff * Leff / ltw;
     if (T0 > -EXPL_THRESHOLD)
     {   T1 = exp(T0);
         T2 = T1 * (1.0 + 2.0 * T1);
         dT1_dVb = -T0 / ltw * T1 * dltw_dVb;
         dT2_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
     }
     else
     {   T1 = MIN_EXPL;
         T2 = T1 * (1.0 + 2.0 * T1);
         dT2_dVb = 0.0;
     }

     T0 = paramPtr->dvt0w * T2;
     DeltVthw = T0 * V0;
     dDeltVthw_dVb = paramPtr->dvt0w * dT2_dVb * V0;
     if (selfheat)   dDeltVthw_dT = T0 * dvbi_dT;
     else   dDeltVthw_dT = 0.0;

     T0 = sqrt(1.0 + paramPtr->nlx / Leff);
     T1 = (paramPtr->kt1 + paramPtr->kt1l / Leff
           + paramPtr->kt2 * Vbs0mos);
     DeltVthtemp = paramPtr->k1eff * (T0 - 1.0) * sqrtPhi + T1
                 * TempRatioMinus1;
     if (selfheat)
        dDeltVthtemp_dT = T1 / model_.tnom;
     else
        dDeltVthtemp_dT = 0.0;

     tmp2 = model_.tox * phi_local
          / (paramPtr->weff + paramPtr->w0);

     T3 = paramPtr->eta0 + paramPtr->etab * Vbs0mos;
     if (T3 < 1.0e-4)   // avoid  discontinuity problems caused by etab
     {   T9 = 1.0 / (3.0 - 2.0e4 * T3);
         T3 = (2.0e-4 - T3) * T9;
         T4 = T9 * T9 * paramPtr->etab;
         dT3_dVb = T4 ;
     }
     else
     {
         dT3_dVb = paramPtr->etab ;
     }
     DIBL_Sft = T3 * paramPtr->theta0vb0 * Vds;
     dDIBL_Sft_dVd = paramPtr->theta0vb0 * T3;
     dDIBL_Sft_dVb = paramPtr->theta0vb0 * Vds * dT3_dVb;

     VthFD = model_.dtype * paramPtr->vth0 + paramPtr->k1eff
         * (sqrtPhis - sqrtPhi) - paramPtr->k2
         * Vbs0mos- Delt_vth - DeltVthw +(paramPtr->k3 + paramPtr->k3b
         * Vbs0mos) * tmp2 + DeltVthtemp - DIBL_Sft;

     T6 = paramPtr->k3b * tmp2 - paramPtr->k2
          + paramPtr->kt2 * TempRatioMinus1;
     dVthFD_dVb = paramPtr->k1eff * dsqrtPhis_dVb
              - dDelt_vth_dVb - dDeltVthw_dVb
              + T6 - dDIBL_Sft_dVb;   //  this is actually dVth_dVbs0mos
     dVthFD_dVe = dVthFD_dVb * dVbs0mos_dVe;
     dVthFD_dVd = -dDIBL_Sft_dVd;
     if (selfheat)
        dVthFD_dT = dDeltVthtemp_dT - dDelt_vth_dT - dDeltVthw_dT
                  + dVthFD_dVb * dVbs0mos_dT;
     else  dVthFD_dT = 0.0;


     // VtgseffFD calculation for PhiFD
     VtgsFD = VthFD - Vgs_eff;
     T10 = model_.nofffd * Vtm;
     CEXP((VtgsFD - model_.vofffd)/ T10, ExpVtgsFD, T0);
     VtgseffFD = T10 * log(1.0 + ExpVtgsFD);
     T0 /= (1.0 + ExpVtgsFD);
     dVtgseffFD_dVd = T0 * dVthFD_dVd;
     dVtgseffFD_dVg = -T0 * dVgs_eff_dVg;
     dVtgseffFD_dVe = T0 * dVthFD_dVe;
     if (selfheat)
        dVtgseffFD_dT = T0 * (dVthFD_dT - (VtgsFD - model_.vofffd)
                        / Temp) + VtgseffFD/Temp;
     else dVtgseffFD_dT = 0.0;


     // surface potential modeling at strong inversion: PhiON
     VgstFD = Vgs_eff - VthFD;
     CEXP((VgstFD - model_.vofffd)/ T10, ExpVgstFD, T0);
     VgsteffFD = T10 * log(1.0 + ExpVgstFD);
     T0 /= (1.0 + ExpVgstFD);
     dVgsteffFD_dVd = -T0 * dVthFD_dVd;
     dVgsteffFD_dVg = T0 * dVgs_eff_dVg;
     dVgsteffFD_dVe = -T0 * dVthFD_dVe;
     if (selfheat)
        dVgsteffFD_dT = T0 * (-dVthFD_dT - (VgstFD - model_.vofffd)
                        / Temp) + VgsteffFD/Temp;
     else dVgsteffFD_dT = 0.0;


     T1 = model_.moinFD*paramPtr->k1eff*Vtm*Vtm;
     if (selfheat) dT1_dT = 2*T1/Temp;
     else dT1_dT=0.0;

     T2 = VgsteffFD+ 2*paramPtr->k1eff*sqrt(phi_local);
     dT2_dVg = dVgsteffFD_dVg;
     dT2_dVd = dVgsteffFD_dVd;
     dT2_dVe = dVgsteffFD_dVe;
     if (selfheat) dT2_dT = dVgsteffFD_dT;
     else dT2_dT = 0.0;

     T0 = 1+ VgsteffFD * T2 / T1;
     dT0_dVg = (VgsteffFD * dT2_dVg + T2 * dVgsteffFD_dVg) / T1;
     dT0_dVd = (VgsteffFD * dT2_dVd + T2 * dVgsteffFD_dVd) / T1;
     dT0_dVe = (VgsteffFD * dT2_dVe + T2 * dVgsteffFD_dVe) / T1;
     if (selfheat)
        dT0_dT = (VgsteffFD * (dT2_dT - T2/T1 * dT1_dT) + T2
                * dVgsteffFD_dT) / T1;
     else dT0_dT = 0.0;


     PhiON = phi_local + Vtm* log(T0) ;
     dPhiON_dVg = Vtm* dT0_dVg/T0 ;
     dPhiON_dVd = Vtm* dT0_dVd/T0 ;
     dPhiON_dVe = Vtm* dT0_dVe/T0 ;
     if (selfheat)
        dPhiON_dT = Vtm* dT0_dT/T0 + (PhiON-phi_local)/Temp ;
     else dPhiON_dT = 0.0;


     // surface potential from subthreshold to inversion: PhiFD
     T0 = model_.cox / (model_.cox +
                             1.0/(1.0/model_.csi + 1.0/Cbox));
     PhiFD = PhiON - T0 * VtgseffFD;
     dPhiFD_dVg = dPhiON_dVg - T0 * dVtgseffFD_dVg;
     dPhiFD_dVd = dPhiON_dVd - T0 * dVtgseffFD_dVd;
     dPhiFD_dVe = dPhiON_dVe - T0 * dVtgseffFD_dVe;
     if (selfheat)
        dPhiFD_dT = dPhiON_dT - T0 * dVtgseffFD_dT;
     else dPhiFD_dT = 0;


     // built-in potential lowering: Vbs0
     T0 = -model_.dvbd1 * paramPtr->leff / paramPtr->litl;
     T1 = model_.dvbd0 * (exp(0.5*T0) + 2*exp(T0));
     T2 = T1 * (vbi - phi_local);
     T3 = 0.5 * paramPtr->qsi / model_.csi;   // v3.2
     Vbs0t = PhiFD - T3 + model_.vbsa + T2;
     dVbs0t_dVg = dPhiFD_dVg;
     dVbs0t_dVd = dPhiFD_dVd;
     dVbs0t_dVe = dPhiFD_dVe;
     if (selfheat)
        dVbs0t_dT = dPhiFD_dT + T1 * dvbi_dT;
     else dVbs0t_dT = 0;


     T0 = 1 + model_.csi / Cbox;
     T3 = -model_.dk2b * paramPtr->leff / paramPtr->litl;
     T5 = model_.k2b * (exp(0.5*T3) + 2*exp(T3));
     T1 = (model_.k1b - T5) / T0;
     T2 = T1 * Vesfb;
     T0 = 1.0/(1 + Cbox / model_.csi);
     Vbs0 = T0 * Vbs0t + T2;
     dVbs0_dVg = T0 * dVbs0t_dVg;
     dVbs0_dVd = T0 * dVbs0t_dVd;
     dVbs0_dVe = T0 * dVbs0t_dVe + T1;
     if (selfheat)
        dVbs0_dT =  T0 * dVbs0t_dT - T1 * dvfbb_dT;
     else
        dVbs0_dT = 0.0;


     // set lowerbound of Vbs (from SPICE) to Vbs0: Vbsitf
     //     (Vbs at back interface)
     if (soiMod == 2) // v3.2  -  v3.1 ideal FD: Vbsitf is pinned at Vbs0
     {
        Vbs = Vbsitf = Vbs0 + OFF_Vbsitf;
        dVbsitf_dVg = dVbs0_dVg;
        dVbsitf_dVd = dVbs0_dVd;
        dVbsitf_dVe = dVbs0_dVe;
        dVbsitf_dVb = 0.0;
        if (selfheat) dVbsitf_dT = dVbs0_dT;
        else dVbsitf_dT = 0;
     }
     else // soiMod = 1
     {
        T1 = Vbs - (Vbs0 + OFF_Vbsitf) - 0.01;
        T2 = sqrt(T1*T1 + 0.0001);
        T3 = 0.5 * (1 + T1/T2);
        Vbsitf = (Vbs0 + OFF_Vbsitf) + 0.5 * (T1 + T2);
        dVbsitf_dVg = (1 - T3) * dVbs0_dVg;
        dVbsitf_dVd = (1 - T3) * dVbs0_dVd;
        dVbsitf_dVe = (1 - T3) * dVbs0_dVe;
        dVbsitf_dVb = T3 ;
        if (selfheat)  dVbsitf_dT = (1 - T3) * dVbs0_dT;
        else  dVbsitf_dT = 0.0;
     }

  // Based on Vbsitf, calculate zero-field body potential for MOS: Vbsmos
     T1 = Vbs0t - Vbsitf - 0.005;
     T2 = sqrt(T1 * T1 + (2.5e-5));
     T3 = 0.5 * (T1 + T2);
     T4 = T3 * model_.csi / paramPtr->qsi; // v3.2
     Vbsmos = Vbsitf - 0.5 * T3 * T4;
     T5 = 0.5 * T4 * (1 + T1 / T2);
     dVbsmos_dVg = dVbsitf_dVg * (1 + T5) - T5 * dVbs0t_dVg;
     dVbsmos_dVd = dVbsitf_dVd * (1 + T5) - T5 * dVbs0t_dVd;
     dVbsmos_dVb = dVbsitf_dVb * (1 + T5);
     dVbsmos_dVe = dVbsitf_dVe * (1 + T5) - T5 * dVbs0t_dVe;
     if (selfheat)
        dVbsmos_dT = dVbsitf_dT * (1 + T5) - T5 * dVbs0t_dT;
     else
        dVbsmos_dT = 0.0;
     // Vbsmos should be used in MOS after some limiting (Vbseff)


     Vbp = Vbs - Vps;
     dVbp_dVb = 1;
  }
// end of v3.0 block edition

// v3.0 modification
  // T2 is Vbsmos limited above Vbsc=-5
  T0 = Vbsmos + 5 - 0.001;
  T1 = sqrt(T0 * T0 - 0.004 * (-5));
  T2 = (-5) + 0.5 * (T0 + T1);
  dT2_dVb = (0.5 * (1.0 + T0 / T1)) * dVbsmos_dVb;
  dT2_dVg = (0.5 * (1.0 + T0 / T1)) * dVbsmos_dVg;
  dT2_dVd = (0.5 * (1.0 + T0 / T1)) * dVbsmos_dVd;
  dT2_dVe = (0.5 * (1.0 + T0 / T1)) * dVbsmos_dVe;
  if (selfheat) dT2_dT = (0.5 * (1.0 + T0 / T1)) * dVbsmos_dT;
  else  dT2_dT = 0.0;

  // Vbsh is T2 limited below 1.5
  T0 = 1.5;
  T1 = T0 - T2 - 0.002;
  T3 = sqrt(T1 * T1 + 0.008 * T0);
  Vbsh = T0 - 0.5 * (T1 + T3);
  dVbsh_dVb = 0.5 * (1.0 + T1 / T3) * dT2_dVb;
  dVbsh_dVg = 0.5 * (1.0 + T1 / T3) * dT2_dVg;
  dVbsh_dVd = 0.5 * (1.0 + T1 / T3) * dT2_dVd;
  dVbsh_dVe = 0.5 * (1.0 + T1 / T3) * dT2_dVe;
  if (selfheat) dVbsh_dT = 0.5 * (1.0 + T1 / T3) * dT2_dT;
  else  dVbsh_dT = 0.0;

  // Vbseff is Vbsh limited to 0.95*phi_local
  T0 = 0.95 * phi_local;
  T1 = T0 - Vbsh - 0.002;
  T2 = sqrt(T1 * T1 + 0.008 * T0);
  Vbseff = T0 - 0.5 * (T1 + T2);
  dVbseff_dVb = 0.5 * (1.0 + T1 / T2) * dVbsh_dVb;
  dVbseff_dVg = 0.5 * (1.0 + T1 / T2) * dVbsh_dVg;
  dVbseff_dVd = 0.5 * (1.0 + T1 / T2) * dVbsh_dVd;
  dVbseff_dVe = 0.5 * (1.0 + T1 / T2) * dVbsh_dVe;
  if (selfheat)  dVbseff_dT = 0.5 * (1.0 + T1 / T2) * dVbsh_dT;
  else  dVbseff_dT = 0.0;
  vbseff = Vbseff;   // SPICE sol.

// end of v3.0 modification

  // Below all the variables refer to Vbseff
  if (dVbseff_dVb < 1e-20) {
     dVbseff_dVb = 1e-20;
     dVbsh_dVb *= 1e20;
  }
  else
     dVbsh_dVb /= dVbseff_dVb;

  Phis = phi_local - Vbseff;
  dPhis_dVb = -1;
  sqrtPhis = sqrt(Phis);
  dsqrtPhis_dVb = -0.5 / sqrtPhis;

  Xdep = Xdep0 * sqrtPhis / sqrtPhi;
  dXdep_dVb = (Xdep0 / sqrtPhi) * dsqrtPhis_dVb;

  // Vth Calculation
  T3 = sqrt(Xdep);

  T0 = paramPtr->dvt2 * Vbseff;
  if (T0 >= - 0.5)
  {   T1 = 1.0 + T0;
      T2 = paramPtr->dvt2 ;
  }
  else
  // Added to avoid any discontinuity problems caused by dvt2
  {   T4 = 1.0 / (3.0 + 8.0 * T0);
      T1 = (1.0 + 3.0 * T0) * T4;
      T2 = paramPtr->dvt2 * T4 * T4 ;
  }
  lt1 = model_.factor1 * T3 * T1;
  dlt1_dVb =model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = paramPtr->dvt2w * Vbseff;
  if (T0 >= - 0.5)
  {   T1 = 1.0 + T0;
      T2 = paramPtr->dvt2w ;
  }
  else
  // Added to avoid any discontinuity problems caused by dvt2w
  {   T4 = 1.0 / (3.0 + 8.0 * T0);
      T1 = (1.0 + 3.0 * T0) * T4;
      T2 = paramPtr->dvt2w * T4 * T4 ;
  }
  ltw= model_.factor1 * T3 * T1;
  dltw_dVb=model_.factor1*(0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = -0.5 * paramPtr->dvt1 * Leff / lt1;
  if (T0 > -EXPL_THRESHOLD)
  {   T1 = exp(T0);
      Theta0 = T1 * (1.0 + 2.0 * T1);
      dT1_dVb = -T0 / lt1 * T1 * dlt1_dVb;
      dTheta0_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
  }
  else
  {   T1 = MIN_EXPL;
      Theta0 = T1 * (1.0 + 2.0 * T1);
      dTheta0_dVb = 0.0;
  }

  thetavth = paramPtr->dvt0 * Theta0;
  Delt_vth = thetavth * V0;
  dDelt_vth_dVb = paramPtr->dvt0 * dTheta0_dVb * V0;
  if (selfheat)  dDelt_vth_dT = thetavth * dvbi_dT;
  else  dDelt_vth_dT = 0.0;

  T0 = -0.5 * paramPtr->dvt1w * paramPtr->weff * Leff / ltw;
  if (T0 > -EXPL_THRESHOLD)
  {   T1 = exp(T0);
      T2 = T1 * (1.0 + 2.0 * T1);
      dT1_dVb = -T0 / ltw * T1 * dltw_dVb;
      dT2_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
  }
  else
  {   T1 = MIN_EXPL;
      T2 = T1 * (1.0 + 2.0 * T1);
      dT2_dVb = 0.0;
  }

  T0 = paramPtr->dvt0w * T2;
  DeltVthw = T0 * V0;
  dDeltVthw_dVb = paramPtr->dvt0w * dT2_dVb * V0;
  if (selfheat)   dDeltVthw_dT = T0 * dvbi_dT;
  else   dDeltVthw_dT = 0.0;

  T0 = sqrt(1.0 + paramPtr->nlx / Leff);
  T1 = (paramPtr->kt1 + paramPtr->kt1l / Leff
        + paramPtr->kt2 * Vbseff);
  DeltVthtemp = paramPtr->k1eff * (T0 - 1.0) * sqrtPhi
                                          + T1 * TempRatioMinus1;
  if (selfheat)
     dDeltVthtemp_dT = T1 / model_.tnom;
  else
     dDeltVthtemp_dT = 0.0;

  tmp2 = model_.tox * phi_local
       / (paramPtr->weff + paramPtr->w0);

  T3 = paramPtr->eta0 + paramPtr->etab * Vbseff;
  if (T3 < 1.0e-4)   // avoid  discontinuity problems caused by etab
  {   T9 = 1.0 / (3.0 - 2.0e4 * T3);
      T3 = (2.0e-4 - T3) * T9;
      T4 = T9 * T9 * paramPtr->etab;
      dT3_dVb = T4 ;
  }
  else
  {
      dT3_dVb = paramPtr->etab ;
  }
  DIBL_Sft = T3 * paramPtr->theta0vb0 * Vds;
  dDIBL_Sft_dVd = paramPtr->theta0vb0 * T3;
  dDIBL_Sft_dVb = paramPtr->theta0vb0 * Vds * dT3_dVb;

  T9 =  2.2361 / sqrtPhi;
  sqrtPhisExt = sqrtPhis - T9 * (Vbsh - Vbseff);
  dsqrtPhisExt_dVb = dsqrtPhis_dVb - T9 * (dVbsh_dVb - 1);

  Vth = model_.dtype * paramPtr->vth0 + paramPtr->k1eff
      * (sqrtPhisExt - sqrtPhi) - paramPtr->k2
      * Vbseff- Delt_vth - DeltVthw +(paramPtr->k3 + paramPtr->k3b
      * Vbseff) * tmp2 + DeltVthtemp - DIBL_Sft;
  von = Vth;

  T6 = paramPtr->k3b * tmp2 - paramPtr->k2
       + paramPtr->kt2 * TempRatioMinus1;
  dVth_dVb = paramPtr->k1eff * dsqrtPhisExt_dVb
           - dDelt_vth_dVb - dDeltVthw_dVb
           + T6 - dDIBL_Sft_dVb;
  //  this is actually dVth_dVbseff

  dVth_dVd = -dDIBL_Sft_dVd;
  if (selfheat)
     dVth_dT = dDeltVthtemp_dT - dDelt_vth_dT - dDeltVthw_dT;
  else
     dVth_dT = 0.0;

  // dVthzb_dT calculation
  if ((model_.capMod == 3) && (selfheat == 1)) {
     T3zb = sqrt(Xdep0);
     ltwzb = lt1zb = model_.factor1 * T3zb;

     T0 = -0.5 * paramPtr->dvt1 * Leff / lt1zb;
     if (T0 > -EXPL_THRESHOLD)
     {   T1 = exp(T0);
         Theta0zb = T1 * (1.0 + 2.0 * T1);
     }
     else
     {   T1 = MIN_EXPL;
         Theta0zb = T1 * (1.0 + 2.0 * T1);
     }
     Delt_vthzb = paramPtr->dvt0 * Theta0zb * V0;
     dDelt_vthzb_dT = paramPtr->dvt0 * Theta0zb * dvbi_dT;

     T0 = -0.5 * paramPtr->dvt1w * paramPtr->weff * Leff / ltwzb;
     if (T0 > -EXPL_THRESHOLD)
     {   T1 = exp(T0);
         T2 = T1 * (1.0 + 2.0 * T1);
     }
     else
     {   T1 = MIN_EXPL;
         T2 = T1 * (1.0 + 2.0 * T1);
     }
     T0 = paramPtr->dvt0w * T2;
     DeltVthwzb = T0 * V0;
     dDeltVthwzb_dT = T0 * dvbi_dT;

     T0 = sqrt(1.0 + paramPtr->nlx / Leff);
     T1 = (paramPtr->kt1 + paramPtr->kt1l / Leff);
     DeltVthtempzb = paramPtr->k1eff * (T0 - 1.0) * sqrtPhi
                   + T1 * TempRatioMinus1;
     dDeltVthtempzb_dT = T1 / model_.tnom;

     Vthzb = model_.dtype * paramPtr->vth0
           - Delt_vthzb - DeltVthwzb + paramPtr->k3 * tmp2
           + DeltVthtempzb;
     dVthzb_dT = dDeltVthtempzb_dT - dDelt_vthzb_dT - dDeltVthwzb_dT;
  }

// Calculate nstar v3.2
  nstar = vtm / CONSTQ * (model_.cox + CONSTEPSSI /
                          Xdep + paramPtr->cit);
// Calculate n
  T2 = paramPtr->nfactor * CONSTEPSSI / Xdep;
  dT2_dVb = - T2 / Xdep * dXdep_dVb;

  T3 = paramPtr->cdsc + paramPtr->cdscb * Vbseff
       + paramPtr->cdscd * Vds;
  dT3_dVb = paramPtr->cdscb;
  dT3_dVd = paramPtr->cdscd;

  T4 = (T2 + T3 * Theta0 + paramPtr->cit) / model_.cox;
  dT4_dVb = (dT2_dVb + Theta0 * dT3_dVb + dTheta0_dVb * T3)
            / model_.cox;
  dT4_dVd = Theta0 * dT3_dVd / model_.cox;

  if (T4 >= -0.5)
  {   n = 1.0 + T4;
      dn_dVb = dT4_dVb;
      dn_dVd = dT4_dVd;
  }
  else
   // avoid  discontinuity problems caused by T4 *
  {   T0 = 1.0 / (3.0 + 8.0 * T4);
      n = (1.0 + 3.0 * T4) * T0;
      T0 *= T0;
      dn_dVb = T0 * dT4_dVb;
      dn_dVd = T0 * dT4_dVd;
  }

// Effective Vgst (Vgsteff_local) Calculation

  Vgst = Vgs_eff - Vth;
  dVgst_dVg = dVgs_eff_dVg;
  dVgst_dVd = -dVth_dVd;
  dVgst_dVb = -dVth_dVb;

  T10 = 2.0 * n * Vtm;
  VgstNVt = Vgst / T10;
  ExpArg = (2.0 * paramPtr->voff - Vgst) / T10;


  // MCJ: Very small Vgst
  if (VgstNVt > EXPL_THRESHOLD)
  {   Vgsteff_local = Vgst;
      // T0 is dVgsteff_dVbseff
      T0 = -dVth_dVb;
      dVgsteff_dVg = dVgs_eff_dVg + T0 * dVbseff_dVg; // v3.0
      dVgsteff_dVd = -dVth_dVd + T0 * dVbseff_dVd;    // v3.0
      dVgsteff_dVb = T0 * dVbseff_dVb;
      dVgsteff_dVe = T0 * dVbseff_dVe; // v3.0
      if (selfheat)
         dVgsteff_dT  = -dVth_dT + T0 * dVbseff_dT;   // v3.0
      else
         dVgsteff_dT = 0.0;
  }
  else if (ExpArg > EXPL_THRESHOLD)
  {   T0 = (Vgst - paramPtr->voff) / (n * Vtm);
      ExpVgst = exp(T0);
      Vgsteff_local = Vtm * paramPtr->cdep0 / model_.cox * ExpVgst;
      T3 = Vgsteff_local / (n * Vtm) ;
      // T1 is dVgsteff_dVbseff
      T1  = -T3 * (dVth_dVb + T0 * Vtm * dn_dVb);
      dVgsteff_dVg = T3 * dVgs_eff_dVg+ T1 * dVbseff_dVg; // v3.0
      dVgsteff_dVd = -T3 * (dVth_dVd + T0 * Vtm * dn_dVd)+
                      T1 * dVbseff_dVd; // v3.0
      dVgsteff_dVe = T1 * dVbseff_dVe;  // v3.0
      dVgsteff_dVb = T1 * dVbseff_dVb;
      if (selfheat)
        dVgsteff_dT = -T3 * (dVth_dT + T0 * dVtm_dT * n)
                      + Vgsteff_local / Temp+ T1 * dVbseff_dT; // v3.0
      else
        dVgsteff_dT = 0.0;
  }
  else
  {   ExpVgst = exp(VgstNVt);
      T1 = T10 * log(1.0 + ExpVgst);

      dT1_dVg = ExpVgst / (1.0 + ExpVgst);
      dT1_dVb = -dT1_dVg * (dVth_dVb + Vgst / n * dn_dVb)
              + T1 / n * dn_dVb;
      dT1_dVd = -dT1_dVg * (dVth_dVd + Vgst / n * dn_dVd)
              + T1 / n * dn_dVd;
      T3 = (1.0 / Temp);
      if (selfheat)
         dT1_dT = -dT1_dVg * (dVth_dT + Vgst * T3) + T1 * T3;
      else
         dT1_dT = 0.0;

      dT2_dVg = -model_.cox / (Vtm * paramPtr->cdep0) * exp(ExpArg);
      T2 = 1.0 - T10 * dT2_dVg;
      dT2_dVd = -dT2_dVg * (dVth_dVd - 2.0 * Vtm * ExpArg * dn_dVd)
              + (T2 - 1.0) / n * dn_dVd;
      dT2_dVb = -dT2_dVg * (dVth_dVb - 2.0 * Vtm * ExpArg * dn_dVb)
              + (T2 - 1.0) / n * dn_dVb;
      if (selfheat)
         dT2_dT = -dT2_dVg * (dVth_dT - ExpArg * T10 * T3);
      else
         dT2_dT = 0.0;

      Vgsteff_local = T1 / T2;
      T3 = T2 * T2;
      //  T4 is dVgsteff_dVbseff
      T4 = (T2 * dT1_dVb - T1 * dT2_dVb) / T3;
      dVgsteff_dVb = T4 * dVbseff_dVb;
      dVgsteff_dVe = T4 * dVbseff_dVe; // v3.0
      dVgsteff_dVg = (T2 * dT1_dVg - T1 * dT2_dVg) / T3 * dVgs_eff_dVg
                                     + T4 * dVbseff_dVg; // v3.0
      dVgsteff_dVd = (T2 * dT1_dVd - T1 * dT2_dVd) / T3 +
                                       T4 * dVbseff_dVd; // v3.0
      if (selfheat)
         dVgsteff_dT = (T2 * dT1_dT - T1 * dT2_dT) / T3 +
                                       T4 * dVbseff_dT;  // v3.0
      else
         dVgsteff_dT = 0.0;
  }
  Vgst2Vtm = Vgsteff_local + 2.0 * Vtm;
  if (selfheat)
     dVgst2Vtm_dT = dVgsteff_dT + 2.0 * dVtm_dT; // v3.1.1 bug fix
  else  dVgst2Vtm_dT = 0.0;
  Vgsteff = Vgsteff_local; // v2.2.3 bug fix

// Calculate Effective Channel Geometry
  T9 = sqrtPhis - sqrtPhi;
  Weff = paramPtr->weff - (2.0 - nbc) *
                   (paramPtr->dwg * Vgsteff_local + paramPtr->dwb * T9);
  dWeff_dVg = -(2.0 - nbc) * paramPtr->dwg;
  dWeff_dVb = -(2.0 - nbc) * paramPtr->dwb * dsqrtPhis_dVb;

  if (Weff < 2.0e-8) // to avoid the discontinuity problem due to Weff
  {   T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
      Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
      T0 *= T0 * 4.0e-16;
      dWeff_dVg *= T0;
      dWeff_dVb *= T0;
  }

  T0 = paramPtr->prwg * Vgsteff_local + paramPtr->prwb * T9;
  if (T0 >= -0.9)
  {   Rds = rds0 * (1.0 + T0);
      dRds_dVg = rds0 * paramPtr->prwg;
      dRds_dVb = rds0 * paramPtr->prwb * dsqrtPhis_dVb;

      if (selfheat && (Rds!=0.0))  dRds_dT = (1.0 + T0) * drds0_dT;
      else  dRds_dT = 0.0;
  }
  else
   // to avoid the discontinuity problem due to prwg and prwb
  {   T1 = 1.0 / (17.0 + 20.0 * T0);
      Rds = rds0 * (0.8 + T0) * T1;
      T1 *= T1;
      dRds_dVg = rds0 * paramPtr->prwg * T1;
      dRds_dVb = rds0 * paramPtr->prwb * dsqrtPhis_dVb * T1;

      if (selfheat && (Rds!=0.0))  dRds_dT = (0.8 + T0) * T1 * drds0_dT;
      else  dRds_dT = 0.0;
  }
  rds = Rds; // v2.2.3 bug fix

// Calculate Abulk_local
  if (paramPtr->a0 == 0.0) {
     Abulk0 = Abulk_local = 1.0;
     dAbulk0_dVb = dAbulk_dVg = dAbulk_dVb = 0.0;
  }
  else {
     T10 = paramPtr->keta * Vbsh;
     if (T10 >= -0.9) {
        T11 = 1.0 / (1.0 + T10);
        dT11_dVb = -paramPtr->keta * T11 * T11 * dVbsh_dVb;
     }
     else { // added to avoid the problems caused by Keta
        T12 = 1.0 / (0.8 + T10);
        T11 = (17.0 + 20.0 * T10) * T12;
        dT11_dVb = -paramPtr->keta * T12 * T12 * dVbsh_dVb;
     }

// v3.0 bug fix
     T10 = phi_local + paramPtr->ketas;

     T13 = (Vbsh * T11) / T10;
     dT13_dVb = (Vbsh * dT11_dVb + T11 * dVbsh_dVb) / T10;

     // limit 1/sqrt(1-T13) to 6, starting at T13=0.96
     if (T13 < 0.96) {
        T14 = 1 / sqrt(1-T13);
        T10 = 0.5 * T14 / (1-T13);
        dT14_dVb = T10 * dT13_dVb;
     }
     else {
        T11 = 1.0 / (1.0 - 1.043406*T13);
        T14 = (6.00167 - 6.26044 * T13) * T11;
        T10 = 0.001742 * T11 * T11;
        dT14_dVb = T10 * dT13_dVb;
     }

// v3.0 bug fix
     T10 = 0.5 * paramPtr->k1eff / sqrt(phi_local + paramPtr->ketas);

     T1 = T10 * T14;
     dT1_dVb = T10 * dT14_dVb;

     T9 = sqrt(paramPtr->xj * Xdep);
     tmp1 = Leff + 2.0 * T9;
     T5 = Leff / tmp1;
     tmp2 = paramPtr->a0 * T5;
     tmp3 = paramPtr->weff + paramPtr->b1;
     tmp4 = paramPtr->b0 / tmp3;
     T2 = tmp2 + tmp4;
     dT2_dVb = -T9 * tmp2 / tmp1 / Xdep * dXdep_dVb;
     T6 = T5 * T5;
     T7 = T5 * T6;

     Abulk0 = 1 + T1 * T2;
     dAbulk0_dVb = T1 * dT2_dVb + T2 * dT1_dVb;

     T8 = paramPtr->ags * paramPtr->a0 * T7;
     dAbulk_dVg = -T1 * T8;
     Abulk_local = Abulk0 + dAbulk_dVg * Vgsteff_local;

     dAbulk_dVb = dAbulk0_dVb
           - T8 * Vgsteff_local * (dT1_dVb + 3.0 * T1 * dT2_dVb / tmp2);
  }

  if (Abulk0 < 0.01)
  {
     T9 = 1.0 / (3.0 - 200.0 * Abulk0);
     Abulk0 = (0.02 - Abulk0) * T9;
     dAbulk0_dVb *= T9 * T9;
  }

  if (Abulk_local < 0.01)
  {
     T9 = 1.0 / (3.0 - 200.0 * Abulk_local);
     Abulk_local = (0.02 - Abulk_local) * T9;
     dAbulk_dVb *= T9 * T9;
     T10 = T9 * T9;        // 3.2 bug fix
     dAbulk_dVg *= T10;    // 3.2 bug fix
  }

  Abulk = Abulk_local;    // v3.2 for noise

// Mobility calculation
  if (model_.mobMod == 1)
  {   T0 = Vgsteff_local + Vth + Vth;
      T2 = ua + uc * Vbseff;
      T3 = T0 / model_.tox;
      T5 = T3 * (T2 + ub * T3);
      dDenomi_dVg = (T2 + 2.0 * ub * T3) / model_.tox;
      dDenomi_dVd = dDenomi_dVg * 2 * dVth_dVd;
      dDenomi_dVb = dDenomi_dVg * 2 * dVth_dVb + uc * T3 ;
      if (selfheat)
        dDenomi_dT = dDenomi_dVg * 2 * dVth_dT
                     + (dua_dT + Vbseff * duc_dT + dub_dT * T3 ) * T3;
       else
        dDenomi_dT = 0.0;
  }
  else if (model_.mobMod == 2)
  {   T5 = Vgsteff_local / model_.tox * (ua + uc * Vbseff
                             + ub * Vgsteff_local / model_.tox);
      dDenomi_dVg = (ua + uc * Vbseff + 2.0 * ub * Vgsteff_local
                                / model_.tox) / model_.tox;
      dDenomi_dVd = 0.0;
      dDenomi_dVb = Vgsteff_local * uc / model_.tox ;
      if (selfheat)
        dDenomi_dT = Vgsteff_local / model_.tox * (dua_dT + Vbseff
                     * duc_dT + dub_dT * Vgsteff_local / model_.tox);
      else
        dDenomi_dT = 0.0;
  }
  else  //  mobMod == 3
  {   T0 = Vgsteff_local + Vth + Vth;
      T2 = 1.0 + uc * Vbseff;
      T3 = T0 / model_.tox;
      T4 = T3 * (ua + ub * T3);
      T5 = T4 * T2;
      dDenomi_dVg = (ua + 2.0 * ub * T3) * T2 / model_.tox;
      dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
      dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + uc * T4 ;
      if (selfheat)
        dDenomi_dT = dDenomi_dVg * 2.0 * dVth_dT + (dua_dT + dub_dT * T3)
                                 * T3 * T2 + T4 * Vbseff * duc_dT;
      else
        dDenomi_dT = 0.0;
  }

  if (T5 >= -0.8)
  {   Denomi = 1.0 + T5;
  }
  else // Added to avoid the discontinuity problem caused by ua and ub
  {   T9 = 1.0 / (7.0 + 10.0 * T5);
      Denomi = (0.6 + T5) * T9;
      T9 *= T9;
      dDenomi_dVg *= T9;
      dDenomi_dVd *= T9;
      dDenomi_dVb *= T9;
      if (selfheat)  dDenomi_dT *= T9;
      else   dDenomi_dT = 0.0;
  }

  ueff = ueff_local = u0temp / Denomi;
  T9 = -ueff_local / Denomi;
  dueff_dVg = T9 * dDenomi_dVg;
  dueff_dVd = T9 * dDenomi_dVd;
  dueff_dVb = T9 * dDenomi_dVb;
  if (selfheat)  dueff_dT = T9 * dDenomi_dT + du0temp_dT / Denomi;
  else  dueff_dT = 0.0;

// Saturation Drain Voltage  Vdsat
  WVCox = Weff * vsattemp * model_.cox;
  WVCoxRds = WVCox * Rds;

//dWVCoxRds_dT = WVCox * dRds_dT
//           + Weff * model_.cox * Rds * dvsattemp_dT;

  Esat = 2.0 * vsattemp / ueff_local;
  EsatL = Esat * Leff;
  T0 = -EsatL /ueff_local;
  dEsatL_dVg = T0 * dueff_dVg;
  dEsatL_dVd = T0 * dueff_dVd;
  dEsatL_dVb = T0 * dueff_dVb;
  if (selfheat)
    dEsatL_dT = T0 * dueff_dT + EsatL / vsattemp * dvsattemp_dT;
  else
    dEsatL_dT = 0.0;

  // Sqrt()
  a1 = paramPtr->a1;
  if (a1 == 0.0)
  {   Lambda = paramPtr->a2;
      dLambda_dVg = 0.0;
  }
  else if (a1 > 0.0)
// Added to avoid the discontinuity problem caused by a1 and a2 (Lambda)
  {   T0 = 1.0 - paramPtr->a2;
      T1 = T0 - paramPtr->a1 * Vgsteff_local - 0.0001;
      T2 = sqrt(T1 * T1 + 0.0004 * T0);
      Lambda = paramPtr->a2 + T0 - 0.5 * (T1 + T2);
      dLambda_dVg = 0.5 * paramPtr->a1 * (1.0 + T1 / T2);
  }
  else
  {   T1 = paramPtr->a2 + paramPtr->a1 * Vgsteff_local - 0.0001;
      T2 = sqrt(T1 * T1 + 0.0004 * paramPtr->a2);
      Lambda = 0.5 * (T1 + T2);
      dLambda_dVg = 0.5 * paramPtr->a1 * (1.0 + T1 / T2);
  }

  AbovVgst2Vtm = Abulk_local /Vgst2Vtm; // v2.2.3 bug fix

  if (Rds > 0)
  {   tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
      tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
  }
  else
  {   tmp2 = dWeff_dVg / Weff;
      tmp3 = dWeff_dVb / Weff;
  }
  if ((Rds == 0.0) && (Lambda == 1.0))
  {   T0 = 1.0 / (Abulk_local * EsatL + Vgst2Vtm);
      tmp1 = 0.0;
      T1 = T0 * T0;
      T2 = Vgst2Vtm * T0;
      T3 = EsatL * Vgst2Vtm;
      Vdsat = T3 * T0;

      dT0_dVg = -(Abulk_local * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
      dT0_dVd = -(Abulk_local * dEsatL_dVd) * T1;
      dT0_dVb = -(Abulk_local * dEsatL_dVb + EsatL * dAbulk_dVb) * T1;
                      if (selfheat)
         dT0_dT  = -(Abulk_local * dEsatL_dT + dVgst2Vtm_dT) * T1;
                      else dT0_dT  = 0.0;

      dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
      dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
      dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
      if (selfheat)
        dVdsat_dT  = T3 * dT0_dT  + T2 * dEsatL_dT
                                       + EsatL * T0 * dVgst2Vtm_dT;
      else dVdsat_dT  = 0.0;
  }
  else
  {   tmp1 = dLambda_dVg / (Lambda * Lambda);
      T9 = Abulk_local * WVCoxRds;
      T8 = Abulk_local * T9;
      T7 = Vgst2Vtm * T9;
      T6 = Vgst2Vtm * WVCoxRds;
      T0 = 2.0 * Abulk_local * (T9 - 1.0 + 1.0 / Lambda);
      dT0_dVg = 2.0 * (T8 * tmp2 - Abulk_local * tmp1
              + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);
//  this is equivalent to one below, but simpler
//  dT0_dVb = 2.0 * (T8 * tmp3 + (2.0 * T9 + 1.0 /
//                                   Lambda - 1.0) * dAbulk_dVg);
      dT0_dVb = 2.0 * (T8 * (2.0 / Abulk_local * dAbulk_dVb + tmp3)
              + (1.0 / Lambda - 1.0) * dAbulk_dVb);
      dT0_dVd = 0.0;

      if (selfheat)
      {
         if (Rds!=0.0)
            tmp4 = dRds_dT / Rds + dvsattemp_dT / vsattemp;
         else
            tmp4 = dvsattemp_dT / vsattemp;

         dT0_dT = 2.0 * T8 * tmp4;

       } else tmp4 = dT0_dT = 0.0;

      T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk_local * EsatL + 3.0 * T7;

      dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
              + Abulk_local * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
              + T7 * tmp2 + T6 * dAbulk_dVg);
      dT1_dVb = Abulk_local * dEsatL_dVb + EsatL * dAbulk_dVb
              + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);
      dT1_dVd = Abulk_local * dEsatL_dVd;

      if (selfheat)
      {
         tmp4 += dVgst2Vtm_dT / Vgst2Vtm;
         dT1_dT  = (2.0 / Lambda - 1.0) * dVgst2Vtm_dT + Abulk_local
                                            * dEsatL_dT + 3.0 * T7 * tmp4;
      } else
         dT1_dT = 0.0;

      T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
      dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
              + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);
      dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
      dT2_dVd = Vgst2Vtm * dEsatL_dVd;
      if (selfheat)
         dT2_dT  = Vgst2Vtm * dEsatL_dT + EsatL * dVgst2Vtm_dT
                 + 2.0 * T6 * (dVgst2Vtm_dT + Vgst2Vtm * tmp4);
      else
         dT2_dT  = 0.0;

      T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
      Vdsat = (T1 - T3) / T0;

      dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
                 - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;
      dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
                 - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;
      dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
      if (selfheat)
         dVdsat_dT  = (dT1_dT - (T1 * dT1_dT - dT0_dT * T2
                    - T0 * dT2_dT) / T3 - Vdsat * dT0_dT) / T0;
      else dVdsat_dT  = 0.0;
  }
  vdsat = Vdsat;

// Effective Vds (Vdseff_local) Calculation
  T1 = Vdsat - Vds - paramPtr->delta;
  dT1_dVg = dVdsat_dVg;
  dT1_dVd = dVdsat_dVd - 1.0;
  dT1_dVb = dVdsat_dVb;
  dT1_dT  = dVdsat_dT;

  T2 = sqrt(T1 * T1 + 4.0 * paramPtr->delta * Vdsat);
  T0 = T1 / T2;
  T3 = 2.0 * paramPtr->delta / T2;
  dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
  dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
  dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;
  if (selfheat)
     dT2_dT  = T0 * dT1_dT  + T3 * dVdsat_dT;
  else
     dT2_dT  = 0.0;

  Vdseff_local = Vdsat - 0.5 * (T1 + T2);
  dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
  dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
  dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);

  if (selfheat)
     dVdseff_dT  = dVdsat_dT  - 0.5 * (dT1_dT  + dT2_dT);
  else
     dVdseff_dT  = 0.0;

  if (Vdseff_local > Vds)
     Vdseff_local = Vds;  // This code is added to fixed the problem
                          // caused by computer precision when
                          // Vds is very close to Vdseff_local.
  diffVds = Vds - Vdseff_local;
  Vdseff = Vdseff_local; // v2.2.3 bug fix

// Calculate VAsat
  tmp4 = 1.0 - 0.5 * Abulk_local * Vdsat / Vgst2Vtm;
  T9 = WVCoxRds * Vgsteff_local;
  T8 = T9 / Vgst2Vtm;
  T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

  T7 = 2.0 * WVCoxRds * tmp4;
  dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff_local)
          - T8 * (Abulk_local * dVdsat_dVg - Abulk_local * Vdsat / Vgst2Vtm
          + Vdsat * dAbulk_dVg);
  dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff_local
          - T8 * (dAbulk_dVb * Vdsat + Abulk_local * dVdsat_dVb);
  dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk_local * dVdsat_dVd;

  if (selfheat)
  {
     if (Rds!=0.0)
        tmp4 = dRds_dT / Rds + dvsattemp_dT / vsattemp;
     else
        tmp4 = dvsattemp_dT / vsattemp;
     dT0_dT  = dEsatL_dT + dVdsat_dT + T7 * tmp4 * Vgsteff_local - T8 *
           (Abulk_local * dVdsat_dT - Abulk_local * Vdsat * dVgst2Vtm_dT
                                                        / Vgst2Vtm);
  } else
     dT0_dT = 0.0;

  T9 = WVCoxRds * Abulk_local;
  T1 = 2.0 / Lambda - 1.0 + T9;
  dT1_dVg = -2.0 * tmp1 +  WVCoxRds * (Abulk_local * tmp2 + dAbulk_dVg);
  dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;
  if (selfheat)
     dT1_dT  = T9 * tmp4;
  else
     dT1_dT  = 0.0;

  Vasat = T0 / T1;
  dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
  dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
  dVasat_dVd = dT0_dVd / T1;
  if (selfheat)
    dVasat_dT  = (dT0_dT  - Vasat * dT1_dT)  / T1;
  else
    dVasat_dT  = 0.0;

// Calculate VACLM
  if ((paramPtr->pclm > 0.0) && (diffVds > 1.0e-10))
  {   T0 = 1.0 / (paramPtr->pclm * Abulk_local * paramPtr->litl);
      dT0_dVb = -T0 / Abulk_local * dAbulk_dVb;
      dT0_dVg = -T0 / Abulk_local * dAbulk_dVg;

      T2 = Vgsteff_local / EsatL;
      T1 = Leff * (Abulk_local + T2);
      dT1_dVg = Leff * ((1.0 - T2 * dEsatL_dVg) / EsatL + dAbulk_dVg);
      dT1_dVb = Leff * (dAbulk_dVb - T2 * dEsatL_dVb / EsatL);
      dT1_dVd = -T2 * dEsatL_dVd / Esat;
      if (selfheat)
         dT1_dT  = -T2 * dEsatL_dT / Esat;
      else
         dT1_dT  = 0.0;

      T9 = T0 * T1;
      VACLM = T9 * diffVds;
      dVACLM_dVg = T0 * dT1_dVg * diffVds - T9 * dVdseff_dVg
                 + T1 * diffVds * dT0_dVg;
      dVACLM_dVb = (dT0_dVb * T1 + T0 * dT1_dVb) * diffVds
                 - T9 * dVdseff_dVb;
      dVACLM_dVd = T0 * dT1_dVd * diffVds + T9 * (1.0 - dVdseff_dVd);
      if (selfheat)
         dVACLM_dT  = T0 * dT1_dT * diffVds - T9 * dVdseff_dT;
      else
         dVACLM_dT  = 0.0;
  }
  else
  {   VACLM = MAX_EXPL;
      dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = dVACLM_dT = 0.0;
  }

// Calculate VADIBL
  if (paramPtr->thetaRout > 0.0)
  {   T8 = Abulk_local * Vdsat;
      T0 = Vgst2Vtm * T8;
      T1 = Vgst2Vtm + T8;
      dT0_dVg = Vgst2Vtm * Abulk_local * dVdsat_dVg + T8
              + Vgst2Vtm * Vdsat * dAbulk_dVg;
      dT1_dVg = 1.0 + Abulk_local * dVdsat_dVg + Vdsat * dAbulk_dVg;
      dT1_dVb = dAbulk_dVb * Vdsat + Abulk_local * dVdsat_dVb;
      dT0_dVb = Vgst2Vtm * dT1_dVb;
      dT1_dVd = Abulk_local * dVdsat_dVd;
      dT0_dVd = Vgst2Vtm * dT1_dVd;
      if (selfheat)
      {
         dT0_dT  = dVgst2Vtm_dT * T8 + Abulk_local * Vgst2Vtm * dVdsat_dT;
         dT1_dT  = dVgst2Vtm_dT + Abulk_local * dVdsat_dT;
      } else
         dT0_dT = dT1_dT = 0.0;

      T9 = T1 * T1;
      T2 = paramPtr->thetaRout;
      VADIBL = (Vgst2Vtm - T0 / T1) / T2;
      dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
      dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
      dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;
      if (selfheat)
         dVADIBL_dT = (dVgst2Vtm_dT - dT0_dT/T1 + T0*dT1_dT/T9) / T2;
      else
         dVADIBL_dT = 0.0;

      T7 = paramPtr->pdiblb * Vbseff;
      if (T7 >= -0.9)
      {   T3 = 1.0 / (1.0 + T7);
          VADIBL *= T3;
          dVADIBL_dVg *= T3;
          dVADIBL_dVb = (dVADIBL_dVb - VADIBL * paramPtr->pdiblb) * T3;
          dVADIBL_dVd *= T3;
          if (selfheat)  dVADIBL_dT  *= T3;
          else  dVADIBL_dT  = 0.0;
      }
      else
// Added to avoid the discontinuity problem caused by pdiblcb
      {   T4 = 1.0 / (0.8 + T7);
          T3 = (17.0 + 20.0 * T7) * T4;
          dVADIBL_dVg *= T3;
          dVADIBL_dVb = dVADIBL_dVb * T3
                      - VADIBL * paramPtr->pdiblb * T4 * T4;
          dVADIBL_dVd *= T3;
          if (selfheat)  dVADIBL_dT  *= T3;
          else  dVADIBL_dT  = 0.0;
          VADIBL *= T3;
      }
  }
  else
  {   VADIBL = MAX_EXPL;
      dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = dVADIBL_dT = 0.0;
  }

// Calculate VA

  T8 = paramPtr->pvag / EsatL;
  T9 = T8 * Vgsteff_local;
  if (T9 > -0.9)
  {   T0 = 1.0 + T9;
      dT0_dVg = T8 * (1.0 - Vgsteff_local * dEsatL_dVg / EsatL);
      dT0_dVb = -T9 * dEsatL_dVb / EsatL;
      dT0_dVd = -T9 * dEsatL_dVd / EsatL;
      if (selfheat)
         dT0_dT  = -T9 * dEsatL_dT / EsatL;
      else
         dT0_dT  = 0.0;
  }
  else // Added to avoid the discontinuity problems caused by pvag
  {   T1 = 1.0 / (17.0 + 20.0 * T9);
      T0 = (0.8 + T9) * T1;
      T1 *= T1;
      dT0_dVg = T8 * (1.0 - Vgsteff_local * dEsatL_dVg / EsatL) * T1;

      T9 *= T1 / EsatL;
      dT0_dVb = -T9 * dEsatL_dVb;
      dT0_dVd = -T9 * dEsatL_dVd;
      if (selfheat)
         dT0_dT  = -T9 * dEsatL_dT;
      else
         dT0_dT  = 0.0;
  }

  tmp1 = VACLM * VACLM;
  tmp2 = VADIBL * VADIBL;
  tmp3 = VACLM + VADIBL;

  T1 = VACLM * VADIBL / tmp3;
  tmp3 *= tmp3;
  dT1_dVg = (tmp1 * dVADIBL_dVg + tmp2 * dVACLM_dVg) / tmp3;
  dT1_dVd = (tmp1 * dVADIBL_dVd + tmp2 * dVACLM_dVd) / tmp3;
  dT1_dVb = (tmp1 * dVADIBL_dVb + tmp2 * dVACLM_dVb) / tmp3;
  if (selfheat)
     dT1_dT  = (tmp1 * dVADIBL_dT  + tmp2 * dVACLM_dT ) / tmp3;
  else
     dT1_dT  = 0.0;

  Va = Vasat + T0 * T1;
  dVa_dVg = dVasat_dVg + T1 * dT0_dVg + T0 * dT1_dVg;
  dVa_dVd = dVasat_dVd + T1 * dT0_dVd + T0 * dT1_dVd;
  dVa_dVb = dVasat_dVb + T1 * dT0_dVb + T0 * dT1_dVb;
  if (selfheat)
     dVa_dT  = dVasat_dT  + T1 * dT0_dT  + T0 * dT1_dT;
  else
     dVa_dT  = 0.0;

// Calculate Ids
  CoxWovL = model_.cox * Weff / Leff;
  beta = ueff_local * CoxWovL;
  dbeta_dVg = CoxWovL * dueff_dVg + beta * dWeff_dVg / Weff ;
  dbeta_dVd = CoxWovL * dueff_dVd;
  dbeta_dVb = CoxWovL * dueff_dVb + beta * dWeff_dVb / Weff ;
  if (selfheat)  dbeta_dT  = CoxWovL * dueff_dT;
  else  dbeta_dT  = 0.0;

  T0 = 1.0 - 0.5 * Abulk_local * Vdseff_local / Vgst2Vtm;
  dT0_dVg = -0.5 * (Abulk_local * dVdseff_dVg
          - Abulk_local * Vdseff_local / Vgst2Vtm + Vdseff_local
                      * dAbulk_dVg) / Vgst2Vtm;
  dT0_dVd = -0.5 * Abulk_local * dVdseff_dVd / Vgst2Vtm;
  dT0_dVb = -0.5 * (Abulk_local * dVdseff_dVb + dAbulk_dVb *
                                        Vdseff_local) / Vgst2Vtm;
  if (selfheat)
     dT0_dT  = -0.5 * (Abulk_local * dVdseff_dT
             - Abulk_local * Vdseff_local / Vgst2Vtm
             * dVgst2Vtm_dT) / Vgst2Vtm;
  else
     dT0_dT = 0.0;

  fgche1 = Vgsteff_local * T0;
  dfgche1_dVg = Vgsteff_local * dT0_dVg + T0;
  dfgche1_dVd = Vgsteff_local * dT0_dVd;
  dfgche1_dVb = Vgsteff_local * dT0_dVb;
  if (selfheat)  dfgche1_dT  = Vgsteff_local * dT0_dT;
  else  dfgche1_dT  = 0.0;

  T9 = Vdseff_local / EsatL;
  fgche2 = 1.0 + T9;
  dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
  dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
  dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;
  if (selfheat)
     dfgche2_dT  = (dVdseff_dT  - T9 * dEsatL_dT)  / EsatL;
  else
     dfgche2_dT  = 0.0;

  gche = beta * fgche1 / fgche2;

  dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
            - gche * dfgche2_dVg) / fgche2;
  dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
            - gche * dfgche2_dVd) / fgche2;
  dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
            - gche * dfgche2_dVb) / fgche2;
  if (selfheat)
     dgche_dT  = (beta * dfgche1_dT  + fgche1 * dbeta_dT
               - gche * dfgche2_dT)  / fgche2;
  else
     dgche_dT  = 0.0;

  T0 = 1.0 + gche * Rds;
  T9 = Vdseff_local / T0;
  Idl = gche * T9;

//  Whoa, these formulas for the derivatives of Idl are convoluted,
//  but I verified them to be correct

  dIdl_dVg = (gche * dVdseff_dVg + T9 * dgche_dVg) / T0
           - Idl * gche / T0 * dRds_dVg ;
  dIdl_dVd = (gche * dVdseff_dVd + T9 * dgche_dVd) / T0;
  dIdl_dVb = (gche * dVdseff_dVb + T9 * dgche_dVb
           - Idl * dRds_dVb * gche) / T0;
  if (selfheat)
     dIdl_dT  = (gche * dVdseff_dT + T9 * dgche_dT
              - Idl * dRds_dT * gche) / T0;
  else
     dIdl_dT  = 0.0;

  T9 =  diffVds / Va;
  T0 =  1.0 + T9;
  ids = Ids = Idl * T0 / nseg;

  Gm0  = T0 * dIdl_dVg - Idl * (      dVdseff_dVg + T9 * dVa_dVg) / Va;
  Gmb0 = T0 * dIdl_dVb - Idl * (      dVdseff_dVb + T9 * dVa_dVb) / Va;
  Gds0 = T0 * dIdl_dVd + Idl * (1.0 - dVdseff_dVd - T9 * dVa_dVd) / Va;

  // v3.1 wanh added for RF
  tmp1 = Gds0 + Gm0 * dVgsteff_dVd;
  tmp2 = Gmb0 + Gm0 * dVgsteff_dVb;
  tmp3 = Gm0;
  // v3.1 wanh added for RF end

  if (selfheat)
     GmT0 = T0 * dIdl_dT - Idl * (dVdseff_dT + T9 * dVa_dT) / Va;
  else
     GmT0 = 0.0;

// This includes all dependencies from Vgsteff_local, Vbseff

  Gm = (Gm0 * dVgsteff_dVg+ Gmb0 * dVbseff_dVg) / nseg; // v3.0
  Gmb = (Gm0 * dVgsteff_dVb + Gmb0 * dVbseff_dVb) / nseg;
  Gds = (Gm0 * dVgsteff_dVd+ Gmb0 * dVbseff_dVd + Gds0) / nseg;
  Gme = (Gm0 * dVgsteff_dVe + Gmb0 * dVbseff_dVe) / nseg; // v3.0
  if (selfheat)
     GmT = (Gm0 * dVgsteff_dT+ Gmb0 * dVbseff_dT + GmT0) / nseg;
  else
     GmT = 0.0;

// v3.1
  if (soiMod != 2) // v3.2
  {
     //  calculate GIDL current
     T0 = 3 * model_.tox;
     // For drain side
     T1 = (Vds - Vgs_eff - paramPtr->ngidl) / T0;
     if ((paramPtr->agidl <= 0.0) || (paramPtr->bgidl <= 0.0) ||
         (T1 <= 0.0))
     {   Idgidl = Gdgidld = Gdgidlg = 0.0;
     }
     else {
        dT1_dVd = 1 / T0;
        dT1_dVg = - dT1_dVd * dVgs_eff_dVg;
        T2 = paramPtr->bgidl / T1;
        if (T2 < EXPL_THRESHOLD)
        {
           Idgidl = wdiod * paramPtr->agidl * T1 * exp(-T2);
           T3 = Idgidl / T1 * (T2 + 1);
           Gdgidld = T3 * dT1_dVd;
           Gdgidlg = T3 * dT1_dVg;
        } else
        {
           T3 = wdiod * paramPtr->agidl * MIN_EXPL;
           Idgidl = T3 * T1 ;
           Gdgidld  = T3 * dT1_dVd;
           Gdgidlg  = T3 * dT1_dVg;
        }
     }
     igidl = Idgidl;

     // For source side
     T1 = (- Vgs_eff - paramPtr->ngidl) / T0;
     if ((paramPtr->agidl <= 0.0) || (paramPtr->bgidl <= 0.0)
                           || (T1 <= 0.0))
     {   Isgidl = Gsgidlg = 0;
     }
     else
     {
        dT1_dVg = - dVgs_eff_dVg / T0;
        T2 = paramPtr->bgidl / T1;
        if (T2 < EXPL_THRESHOLD)
        {
           Isgidl = wdios * paramPtr->agidl * T1 * exp(-T2);
           T3 = Isgidl / T1 * (T2 + 1);
           Gsgidlg = T3 * dT1_dVg;
        } else
        {
           T3 = wdios * paramPtr->agidl * MIN_EXPL;
           Isgidl = T3 * T1 ;
           Gsgidlg = T3 * dT1_dVg;
        }
     }

     // calculate diode and BJT current
     WsTsi = wdios * model_.tsi;
     WdTsi = wdiod * model_.tsi;

     NVtm1 = Vtm * paramPtr->ndiode;
     if (selfheat)
        dNVtm1_dT = paramPtr->ndiode * dVtm_dT;
     else
        dNVtm1_dT = 0;

     T0 = Vbs / NVtm1;
     dT0_dVb = 1.0 / NVtm1;
     if (selfheat)
        dT0_dT = -Vbs / NVtm1 / NVtm1 * dNVtm1_dT;
     else
        dT0_dT = 0;
     CEXP(T0, ExpVbsNVtm, T1);
     dExpVbsNVtm_dVb = T1 * dT0_dVb;
     if (selfheat)
        dExpVbsNVtm_dT = T1 * dT0_dT;
     else
        dExpVbsNVtm_dT = 0;

     T0 = Vbd / NVtm1;
     dT0_dVb = 1.0 / NVtm1;
     dT0_dVd = -dT0_dVb;
     if (selfheat)
        dT0_dT = -Vbd / NVtm1 / NVtm1 * dNVtm1_dT;
     else
        dT0_dT = 0;
     CEXP(T0, ExpVbdNVtm, T1);
     dExpVbdNVtm_dVb = T1 * dT0_dVb;
     dExpVbdNVtm_dVd = -dExpVbdNVtm_dVb;
     if (selfheat)
        dExpVbdNVtm_dT = T1 * dT0_dT;
     else
        dExpVbdNVtm_dT = 0;

     // Ibs1 / Ibd1 : diffusion current
     if (jdif == 0) {
        Ibs1 = dIbs1_dVb = dIbs1_dT = Ibd1 = dIbd1_dVb
                                        = dIbd1_dVd = dIbd1_dT = 0;
     }
     else {
        T0 = WsTsi * jdif;
        if (selfheat)
           dT0_dT = WsTsi * djdif_dT;
        else
           dT0_dT = 0;
        Ibs1 = T0 * (ExpVbsNVtm - 1);
        dIbs1_dVb = T0 * dExpVbsNVtm_dVb;
        if (selfheat)
           dIbs1_dT = T0 * dExpVbsNVtm_dT + (ExpVbsNVtm - 1) * dT0_dT;
        else
           dIbs1_dT = 0;

        T0 = WdTsi * jdif;
        if (selfheat)
           dT0_dT = WdTsi * djdif_dT;
        else
           dT0_dT = 0;
        Ibd1 = T0 * (ExpVbdNVtm - 1);
        dIbd1_dVb = T0 * dExpVbdNVtm_dVb;
        dIbd1_dVd = -dIbd1_dVb;
        if (selfheat)
           dIbd1_dT = T0 * dExpVbdNVtm_dT + (ExpVbdNVtm -1) * dT0_dT;
        else
           dIbd1_dT = 0;
     }

     // Ibs2:recombination/trap-assisted tunneling current
     NVtmf = 0.026 * paramPtr->nrecf0
           * (1 + paramPtr->ntrecf * (TempRatio - 1));
     NVtmr = 0.026 * paramPtr->nrecr0 // v2.2.2 bug fix
           * (1 + paramPtr->ntrecr * (TempRatio - 1));
     if (selfheat) {
        dNVtmf_dT = paramPtr->nrecf0 * 0.026
                  * paramPtr->ntrecf * dTempRatio_dT;
        dNVtmr_dT = paramPtr->nrecr0 * 0.026 // v2.2.2 bug fix
                  * paramPtr->ntrecr * dTempRatio_dT;
     }
     else
        dNVtmf_dT = dNVtmr_dT = 0;

     if (jrec == 0) {
        Ibs2 = dIbs2_dVb = dIbs2_dT = 0;
        Ibd2 = dIbd2_dVb = dIbd2_dVd = dIbd2_dT = 0;
     }
     else {
        // forward bias
        T0 = Vbs / NVtmf;
        CEXP(T0,T10,T2);
        T4 = 1 / NVtmf;
        dT10_dVb = T4 * T2;
        if (selfheat)
           dT10_dT  = - T4 * T2 * Vbs / NVtmf * dNVtmf_dT ;
        else   dT10_dT  = 0.0;

        // reverse bias
        if ((paramPtr->vrec0 - Vbs) < 1e-3) {

        // v2.2.3 bug fix
           T1 = 1e3;
           T0 = -Vbs / NVtmr * paramPtr->vrec0 * T1;
           T11 = -exp(T0);

           dT11_dVb = dT11_dT = 0;
        }
        else {
           T1 = 1 / (paramPtr->vrec0 - Vbs);
           T0 = -Vbs / NVtmr * paramPtr->vrec0 * T1;
           dT0_dVb = -paramPtr->vrec0 / NVtmr * (T1 + Vbs * T1 * T1) ;
           if (selfheat)
              dT0_dT = -T0 / NVtmr * dNVtmr_dT;
           else   dT0_dT = 0;

           CEXP(T0, T11, T2);
           T11 = -T11;
           dT11_dVb = -T2 * dT0_dVb;
           if (selfheat)
              dT11_dT = -T2 * dT0_dT;
           else   dT11_dT = 0;
        }
        T3 = WsTsi * jrec;
        Ibs2 = T3 * (T10 + T11);
        dIbs2_dVb = T3 * (dT10_dVb + dT11_dVb);
        if (selfheat)
           dIbs2_dT = T3 * (dT10_dT + dT11_dT) + WsTsi * (T10 + T11)
                                                       * djrec_dT;
        else   dIbs2_dT = 0;

        // Ibd2
        T0 = Vbd / NVtmf;
        CEXP(T0,T10,T2);
        T4 = 1 / NVtmf;
        dT10_dVb = T4 * T2;
        if (selfheat)
        {
           dT10_dT  = - T4 * T2 * Vbd / NVtmf * dNVtmf_dT ;
        }
        else
        {
          dT10_dT  = 0.0;
        }

        if ((paramPtr->vrec0 - Vbd) < 1e-3) {

        // v2.2.3 bug fix
           T1 = 1e3;
           T0 = -Vbd / NVtmr * paramPtr->vrec0 * T1;
           T11 = -exp(T0);

           dT11_dVb = dT11_dT = 0;
        }
        else {
           T1 = 1 / (paramPtr->vrec0 - Vbd);
           T0 = -Vbd / NVtmr * paramPtr->vrec0 * T1;
           dT0_dVb = -paramPtr->vrec0 / NVtmr * (T1 + Vbd * T1 * T1) ;
           if (selfheat)
              dT0_dT = -T0 / NVtmr * dNVtmr_dT;
           else
              dT0_dT = 0;
           CEXP(T0, T11, T2);
           T11 = - T11;
           dT11_dVb = -T2 * dT0_dVb;
           if (selfheat)
              dT11_dT = -T2 * dT0_dT;
           else
              dT11_dT = 0;
        }
        T3 = WdTsi * jrec;
        Ibd2 = T3 * (T10 + T11);
        dIbd2_dVb = T3 * (dT10_dVb + dT11_dVb);

        dIbd2_dVd = -dIbd2_dVb;
        if (selfheat)
           dIbd2_dT = T3 * (dT10_dT + dT11_dT) + WdTsi * (T10 + T11)
                                                       * djrec_dT;
        else
           dIbd2_dT = 0;
     }

     // Ibs3/Ibd3:  recombination current in neutral body
     WTsi = paramPtr->weff / nseg * model_.tsi;
     if (jbjt == 0.0)
     {
        Ibs3 = dIbs3_dVb = dIbs3_dVd = dIbs3_dT = 0.0;
        Ibd3 = dIbd3_dVb = dIbd3_dVd = dIbd3_dT = 0.0;
        Ibsdif = dIbsdif_dVb = dIbsdif_dT = 0;
        Ibddif = dIbddif_dVb = dIbddif_dVd = dIbddif_dT = 0;
        ic = Ic = Gcd = Gcb = GcT = 0.0;
     }
     else {
        Ien = WTsi * jbjt * paramPtr->lratio;
        if (selfheat)
           dIen_dT = WTsi * djbjt_dT * paramPtr->lratio;
        else
           dIen_dT = 0;

        // high level injection of source side
        if ((Ehlis = Ahli * (ExpVbsNVtm - 1)) < 1e-5) {
           Ehlis = dEhlis_dVb = dEhlis_dT = 0;
           EhlisFactor = 1;
           dEhlisFactor_dVb = dEhlisFactor_dT = 0;
        }
        else {
           dEhlis_dVb = Ahli * dExpVbsNVtm_dVb;
           if (selfheat)
              dEhlis_dT = Ahli * dExpVbsNVtm_dT + (ExpVbsNVtm - 1) * dAhli_dT;
           else
              dEhlis_dT = 0;
           EhlisFactor = 1.0 / sqrt(1 + Ehlis);
           T0 = -0.5 * EhlisFactor / (1 + Ehlis);
           dEhlisFactor_dVb = T0 * dEhlis_dVb;
           if (selfheat)
              dEhlisFactor_dT = T0 * dEhlis_dT;
           else
              dEhlisFactor_dT = 0;
        }

        // high level injection of drain side
        if ((Ehlid = Ahli * (ExpVbdNVtm - 1)) < 1e-5) {
           Ehlid = dEhlid_dVb = dEhlid_dVd = dEhlid_dT = 0;
           EhlidFactor = 1;
           dEhlidFactor_dVb = dEhlidFactor_dVd = dEhlidFactor_dT = 0;
        }
        else {
           dEhlid_dVb = Ahli * dExpVbdNVtm_dVb;
           dEhlid_dVd = -dEhlid_dVb;
           if (selfheat)
              dEhlid_dT = Ahli * dExpVbdNVtm_dT + (ExpVbdNVtm - 1) * dAhli_dT;
           else
              dEhlid_dT = 0;
           EhlidFactor = 1.0 / sqrt(1 + Ehlid);
           T0 = -0.5 * EhlidFactor / (1 + Ehlid);
           dEhlidFactor_dVb = T0 * dEhlid_dVb;
           dEhlidFactor_dVd = -dEhlidFactor_dVb;
           if (selfheat)
              dEhlidFactor_dT = T0 * dEhlid_dT;
           else
              dEhlidFactor_dT = 0;
        }


        // v3.1.1 bug fix for Ibjt(L) discontinuity
        T0 = 1 - paramPtr->arfabjt;
        T1 = T0 * Ien;
        if (selfheat)
           dT1_dT = T0 * dIen_dT;
        else
           dT1_dT = 0;

        Ibs3 = T1 * (ExpVbsNVtm - 1) * EhlisFactor;
        dIbs3_dVb = T1 * (dExpVbsNVtm_dVb * EhlisFactor
                  + (ExpVbsNVtm - 1) * dEhlisFactor_dVb);
        dIbs3_dVd = 0;
        if (selfheat)
           dIbs3_dT = dT1_dT * (ExpVbsNVtm - 1) * EhlisFactor
                    + T1 * (dExpVbsNVtm_dT * EhlisFactor
                    + (ExpVbsNVtm - 1) * dEhlisFactor_dT);
        else
           dIbs3_dT = 0.0;

        Ibd3 = T1 * (ExpVbdNVtm - 1) * EhlidFactor;
        dIbd3_dVb = T1 * (dExpVbdNVtm_dVb * EhlidFactor
                  + (ExpVbdNVtm - 1) * dEhlidFactor_dVb);
        dIbd3_dVd = -dIbd3_dVb;
        if (selfheat)
           dIbd3_dT = dT1_dT * (ExpVbdNVtm - 1) * EhlidFactor
                    + T1 * (dExpVbdNVtm_dT * EhlidFactor
                    + (ExpVbdNVtm - 1) * dEhlidFactor_dT);
        else
           dIbd3_dT = 0.0;


        // effective diffusion current for capacitance calcu.
        Iendif = WTsi * jbjt * paramPtr->lratiodif;
        if (selfheat)
           dIendif_dT = WTsi * djbjt_dT * paramPtr->lratiodif;
        else
           dIendif_dT = 0;

        Ibsdif = Iendif * (ExpVbsNVtm - 1) * EhlisFactor;
        dIbsdif_dVb = Iendif * (dExpVbsNVtm_dVb * EhlisFactor
                    + (ExpVbsNVtm - 1) * dEhlisFactor_dVb);
        if (selfheat)
           dIbsdif_dT = dIendif_dT * (ExpVbsNVtm - 1) * EhlisFactor
                      + Iendif * (dExpVbsNVtm_dT * EhlisFactor
                      + (ExpVbsNVtm - 1) * dEhlisFactor_dT);
        else
           dIbsdif_dT = 0;

        Ibddif = Iendif * (ExpVbdNVtm - 1) * EhlidFactor;
        dIbddif_dVb = Iendif * (dExpVbdNVtm_dVb * EhlidFactor
                    + (ExpVbdNVtm - 1) * dEhlidFactor_dVb);
        dIbddif_dVd = -dIbddif_dVb;
        if (selfheat)
           dIbddif_dT = dIendif_dT * (ExpVbdNVtm - 1) * EhlidFactor
                      + Iendif * (dExpVbdNVtm_dT * EhlidFactor
                      + (ExpVbdNVtm - 1) * dEhlidFactor_dT);
        else
           dIbddif_dT = 0;

        // Ic: Bjt collector current
        if ((bjtoff == 1) || (Vds == 0.0)) {
           ic = Ic = Gcd = Gcb = GcT = 0.0;
        }
        else {
           // second order effects
           T0 = 1 + (Vbs + Vbd) / paramPtr->vearly;
           dT0_dVb = 2.0 / paramPtr->vearly;
           dT0_dVd = -1.0 / paramPtr->vearly;

           T1 = Ehlis + Ehlid;
           dT1_dVb = dEhlis_dVb + dEhlid_dVb;
           dT1_dVd = dEhlid_dVd;
           if (selfheat)
              dT1_dT = dEhlis_dT + dEhlid_dT;
           else
              dT1_dT = 0;

           T3 = sqrt(T0 * T0 + 4 * T1);
           dT3_dVb = 0.5 / T3 * (2 * T0 * dT0_dVb + 4 * dT1_dVb);
           dT3_dVd = 0.5 / T3 * (2 * T0 * dT0_dVd + 4 * dT1_dVd);
           if (selfheat)
              dT3_dT = 2 * dT1_dT / T3;
           else
              dT3_dT = 0;

           T2 = (T0 + T3) / 2.0;
           dT2_dVb = (dT0_dVb + dT3_dVb) / 2.0;
           dT2_dVd = (dT0_dVd + dT3_dVd) / 2.0;
           if (selfheat)
              dT2_dT = dT3_dT /2.0;
           else
              dT2_dT = 0;

           if (T2 < .1)
           {
              E2ndFactor = 10.0;
              dE2ndFactor_dVb = dE2ndFactor_dVd = dE2ndFactor_dT = 0;
           }

           else {
              E2ndFactor = 1.0 / T2;
              dE2ndFactor_dVb = -E2ndFactor / T2 * dT2_dVb;
              dE2ndFactor_dVd = -E2ndFactor / T2 * dT2_dVd;
              if (selfheat)
                 dE2ndFactor_dT = -E2ndFactor / T2 * dT2_dT;
              else
                 dE2ndFactor_dT = 0;
           }

           T0 = paramPtr->arfabjt * Ien;
           if (selfheat)
              dT0_dT = paramPtr->arfabjt * dIen_dT;
           else
              dT0_dT = 0;
           ic = Ic
                         = T0 * (ExpVbsNVtm - ExpVbdNVtm) * E2ndFactor;
           Gcb = T0 * ((dExpVbsNVtm_dVb - dExpVbdNVtm_dVb) * E2ndFactor
               + (ExpVbsNVtm - ExpVbdNVtm) * dE2ndFactor_dVb);
           Gcd = T0 * (-dExpVbdNVtm_dVd * E2ndFactor
               + (ExpVbsNVtm - ExpVbdNVtm) * dE2ndFactor_dVd);
           if (selfheat)
              GcT = T0 * (dExpVbsNVtm_dT - dExpVbdNVtm_dT) * E2ndFactor
                  + dT0_dT * (ExpVbsNVtm - ExpVbdNVtm) * E2ndFactor
                  + T0 * (ExpVbsNVtm - ExpVbdNVtm) * dE2ndFactor_dT;
           else
              GcT = 0;
        }
     }

     // Ibs4/Ibd4 : tunneling
     NVtm2 = 0.026 * paramPtr->ntun;
     if (jtun == 0)
     {  Ibs4 = Ibd4 = dIbs4_dVb = dIbs4_dT = dIbd4_dVb
                                          = dIbd4_dVd = dIbd4_dT = 0;
     } else
     {
        if ((paramPtr->vtun0 - Vbs) < 1e-3)
        {
           // v2.2.3 bug fix
           T1=1e3;
           T0 = -Vbs / NVtm2 * paramPtr->vtun0 * T1;
           T1 = exp(T0);
           T3 = WsTsi * jtun;
           Ibs4 = T3 * (1- T1);

           dIbs4_dVb = dIbs4_dT = 0;
        }
        else {
           T1 = 1 / (paramPtr->vtun0 - Vbs);
           T0 = -Vbs / NVtm2 * paramPtr->vtun0 * T1;
           dT0_dVb = -paramPtr->vtun0 / NVtm2 * (T1 + Vbs * T1 * T1) ;

           CEXP(T0, T1, T2);
           T3 = WsTsi * jtun;
           Ibs4 =  T3 * (1- T1);
           dIbs4_dVb = -T3 * T2 * dT0_dVb;
           if (selfheat)
              dIbs4_dT = (1 - T1) * WsTsi * djtun_dT;
           else
              dIbs4_dT = 0;
        }

        if ((paramPtr->vtun0 - Vbd) < 1e-3)
        {
           // v2.2.3 bug fix
           T1=1e3;
           T0 = -Vbd / NVtm2 * paramPtr->vtun0 * T1;
           T1 = exp(T0);
           T3 = WdTsi * jtun;
           Ibd4 = T3 * (1- T1);

           dIbd4_dVb = dIbd4_dT = 0;
           dIbd4_dVd = 0;
        }
        else
        {
           T1 = 1 / (paramPtr->vtun0 - Vbd);
           T0 = -Vbd / NVtm2 * paramPtr->vtun0 * T1;
           dT0_dVb = -paramPtr->vtun0 / NVtm2 * (T1 + Vbd * T1 * T1) ;

           CEXP(T0, T1, T2);
           T3 = WdTsi * jtun;
           Ibd4 =  T3 * (1- T1);
           dIbd4_dVb = -T3 * T2 * dT0_dVb;

           dIbd4_dVd = -dIbd4_dVb;

           if (selfheat)
              dIbd4_dT = (1 - T1) * WdTsi * djtun_dT;
           else   dIbd4_dT = 0;
        }
     }

     itun = - Ibd3 - Ibd4;
     ibs = Ibs = Ibs1 + Ibs2 + Ibs3 + Ibs4;
     ibd = Ibd = Ibd1 + Ibd2 + Ibd3 + Ibd4;

     Gjsb = dIbs1_dVb + dIbs2_dVb + dIbs3_dVb + dIbs4_dVb;
     Gjsd = dIbs3_dVd;
     if (selfheat)
        GjsT = dIbs1_dT + dIbs2_dT + dIbs3_dT + dIbs4_dT;
     else
        GjsT = 0.0;

     Gjdb = dIbd1_dVb + dIbd2_dVb + dIbd3_dVb + dIbd4_dVb;

     Gjdd = dIbd1_dVd + dIbd2_dVd + dIbd3_dVd + dIbd4_dVd;
     if (selfheat)
        GjdT = dIbd1_dT  + dIbd2_dT + dIbd3_dT + dIbd4_dT;
     else
        GjdT = 0.0;

  }
  else // v3.1 soiMod=2: ideal FD
  {
     igidl= Idgidl = Gdgidld = Gdgidlg = 0.0;
     Isgidl = Gsgidlg = 0;

     itun = 0;
     ibs = Ibs = 0;
     ibd = Ibd = 0;
     ic = Ic = Gcd = Gcb = GcT = 0.0;

     Gjsb = 0;
     Gjsd = 0;
     GjsT = 0;

     Gjdb = 0;
     Gjdd = 0;
     GjdT = 0;
  }

  // v3.0: gate-tunneling
  if ((model_.igbMod != 0) || (model_.igcMod != 0)) {
     Vgb = Vgs_eff - Vbs;
     dVgb_dVg = dVgs_eff_dVg;
     dVgb_dVb = -1;

     // Calculate Vox first
     Vfb = model_.dtype * paramPtr->vth0 - phi_local
                                      - paramPtr->k1eff * sqrtPhi;

     T3 = Vfb - Vgs_eff + Vbs - DELTA_3;
     dT3_dVg = -dVgs_eff_dVg;
     dT3_dVd = 0;
     dT3_dVb = 1;

     if (Vfb <= 0.0) {
        T0 = sqrt(T3 * T3 - 4.0 * DELTA_3 * Vfb);
        dT0_dVg = 1.0/(2.0 * T0) * 2.0*T3 * dT3_dVg;
        dT0_dVb = 0.5*(1.0/T0) * 2.0*T3 * dT3_dVb;
     }
     else {
        T0 = sqrt(T3 * T3 + 4.0 * DELTA_3 * Vfb);
        dT0_dVg = 1.0/(2.0 * T0) * 2.0*T3 * dT3_dVg;
        dT0_dVb = 0.5*(1.0/T0) * 2.0*T3 * dT3_dVb;
     }

     Vfbeff = Vfb - 0.5 * (T3 + T0);
     dVfbeff_dVg = -0.5 * (dT3_dVg + dT0_dVg);
     dVfbeff_dVb = -0.5 * (dT3_dVb + dT0_dVb);

     Voxacc = Vfb - Vfbeff;
     dVoxacc_dVg = -dVfbeff_dVg;
     dVoxacc_dVd = 0.0;
     dVoxacc_dVb = -dVfbeff_dVb;
     if (Voxacc < 0.0)
        Voxacc = dVoxacc_dVg = dVoxacc_dVb = 0.0;

     T0 = Vgs_eff - Vgsteff_local - Vfbeff - Vbseff;
     dT0_dVg = dVgs_eff_dVg - dVgsteff_dVg - dVfbeff_dVg
                             - dVbseff_dVg; // v3.0
     dT0_dVd = -dVgsteff_dVd - dVbseff_dVd; // v3.0
     dT0_dVb = -dVgsteff_dVb - dVfbeff_dVb - dVbseff_dVb;
     dT0_dVe = -dVgsteff_dVe - dVbseff_dVe;

     if (selfheat)
        dT0_dT = -dVgsteff_dT - dVbseff_dT; // v3.0

     if (paramPtr->k1eff == 0.0) {
        Voxdepinv = dVoxdepinv_dVg
                  = dVoxdepinv_dVd
                  = dVoxdepinv_dVb
                  = dVoxdepinv_dT
                  = 0.0;
     } else {
        if (T0 < 0.0) {
           T1 = T0/paramPtr->k1eff;
           dT1_dVg = dT0_dVg/paramPtr->k1eff;
           dT1_dVd = dT0_dVd/paramPtr->k1eff;
           dT1_dVb = dT0_dVb/paramPtr->k1eff;
           dT1_dVe = dT0_dVe/paramPtr->k1eff; // v3.0
           if (selfheat) dT1_dT = dT0_dT/paramPtr->k1eff;
        }
        else {
           T1 = paramPtr->k1eff/2*(-1 + sqrt(1 +
                4*T0/paramPtr->k1eff/paramPtr->k1eff));
           T2 = paramPtr->k1eff/2 *
                0.5/sqrt(1 + 4*T0/paramPtr->k1eff
                           /paramPtr->k1eff) *
                4/paramPtr->k1eff/paramPtr->k1eff;
           dT1_dVg = T2 * dT0_dVg;
           dT1_dVd = T2 * dT0_dVd;
           dT1_dVb = T2 * dT0_dVb;
           dT1_dVe = T2 * dT0_dVe; // v3.0
           if (selfheat)
              dT1_dT = T2 * dT0_dT;
         }

        Voxdepinv = Vgs_eff - (T1*T1 + Vbs) - Vfb;
        dVoxdepinv_dVg = dVgs_eff_dVg - (2.0*T1*dT1_dVg);
        dVoxdepinv_dVd = -(2.0*T1*dT1_dVd);
        dVoxdepinv_dVb = -(2.0*T1*dT1_dVb + 1);
        dVoxdepinv_dVe = -(2.0*T1*dT1_dVe); // v3.0
        if (selfheat)
           dVoxdepinv_dT = -(2.0*T1*dT1_dT);
     }
  }


  // gate-channel tunneling component
  if (model_.igcMod)
  {   T0 = Vtm * paramPtr->nigc;
      VxNVt = (Vgs_eff - model_.dtype * paramPtr->vth0) / T0;
                            // Vth instead of Vth0 may be used
      if (VxNVt > EXPL_THRESHOLD)
      {   Vaux = Vgs_eff - model_.dtype * paramPtr->vth0;
          dVaux_dVg = dVgs_eff_dVg;
          dVaux_dVd = 0.0;
          dVaux_dVb = 0.0;
      }
      else if (VxNVt < -EXPL_THRESHOLD)
      {   Vaux = T0 * log(1.0 + MIN_EXPL);
          dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
      }
      else
      {   ExpVxNVt = exp(VxNVt);
          Vaux = T0 * log(1.0 + ExpVxNVt);
          dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
          dVaux_dVd = -dVaux_dVg * 0.0;
          dVaux_dVb = -dVaux_dVg * 0.0;
          dVaux_dVg *= dVgs_eff_dVg;
      }

      T2 = Vgs_eff * Vaux;
      dT2_dVg = dVgs_eff_dVg * Vaux + Vgs_eff * dVaux_dVg;
      dT2_dVd = Vgs_eff * dVaux_dVd;
      dT2_dVb = Vgs_eff * dVaux_dVb;

      T11 = paramPtr->Aechvb;
      T12 = paramPtr->Bechvb;
      T3 = paramPtr->aigc * paramPtr->cigc - paramPtr->bigc;
      T4 = paramPtr->bigc * paramPtr->cigc;
      T5 = T12 * (paramPtr->aigc + T3 * Voxdepinv
         - T4 * Voxdepinv * Voxdepinv);

      if (T5 > EXPL_THRESHOLD)
      {   T6 = MAX_EXPL;
          dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
      }
      else if (T5 < -EXPL_THRESHOLD)
      {   T6 = MIN_EXPL;
          dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
      }
      else
      {   T6 = exp(T5);
          dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
          dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
          dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
          dT6_dVg *= dVoxdepinv_dVg;
      }

      Igc = T11 * T2 * T6;
      dIgc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
      dIgc_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
      dIgc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

      T7 = -paramPtr->pigcd * Vds;
      T8 = T7 * T7 + 2.0e-4;
      dT8_dVd = -2.0 * paramPtr->pigcd * T7;
      if (T7 > EXPL_THRESHOLD)
      {   T9 = MAX_EXPL;
          dT9_dVd = 0.0;
      }
      else if (T7 < -EXPL_THRESHOLD)
      {   T9 = MIN_EXPL;
          dT9_dVd = 0.0;
      }
      else
      {   T9 = exp(T7);
          dT9_dVd = -T9 * paramPtr->pigcd;
      }

      T0 = T8 * T8;
      T1 = T9 - 1.0 + 1.0e-4;
      T10 = (T1 - T7) / T8;
      dT10_dVd = ((paramPtr->pigcd + dT9_dVd) * T8
               - (T1 - T7) * dT8_dVd) / T0;
      Igcs_local = Igc * T10;
      dIgcs_dVg = dIgc_dVg * T10;
      dIgcs_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
      dIgcs_dVb = dIgc_dVb * T10;

      T1 = T9 - 1.0 - 1.0e-4;
      T10 = (T7 * T9 - T1) / T8;
      dT10_dVd = (-paramPtr->pigcd * T9 + (T7 - 1.0)
               * dT9_dVd - T10 * dT8_dVd) / T8;
      Igcd_local = Igc * T10;
      dIgcd_dVg = dIgc_dVg * T10;
      dIgcd_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
      dIgcd_dVb = dIgc_dVb * T10;

      Igcs = Igcs_local;
      gIgcsg = dIgcs_dVg;
      gIgcsd = dIgcs_dVd;
      gIgcsb =  dIgcs_dVb * dVbseff_dVb;
      Igcd = Igcd_local;
      gIgcdg = dIgcd_dVg;
      gIgcdd = dIgcd_dVd;
      gIgcdb = dIgcd_dVb * dVbseff_dVb;


      T0 = vgs - paramPtr->vfbsd;
      vgs_eff = sqrt(T0 * T0 + 1.0e-4);
      dvgs_eff_dvg = T0 / vgs_eff;

      T2 = vgs * vgs_eff;
      dT2_dVg = vgs * dvgs_eff_dvg + vgs_eff;
      T11 = paramPtr->AechvbEdge;
      T12 = paramPtr->BechvbEdge;
      T3 = paramPtr->aigsd * paramPtr->cigsd
         - paramPtr->bigsd;
      T4 = paramPtr->bigsd * paramPtr->cigsd;
      T5 = T12 * (paramPtr->aigsd + T3 * vgs_eff - T4 *
                                         vgs_eff * vgs_eff);
      if (T5 > EXPL_THRESHOLD)
      {   T6 = MAX_EXPL;
          dT6_dVg = 0.0;
      }
      else if (T5 < -EXPL_THRESHOLD)
      {   T6 = MIN_EXPL;
          dT6_dVg = 0.0;
      }
      else
      {   T6 = exp(T5);
          dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgs_eff)
                  * dvgs_eff_dvg;
      }
      Igs_local = T11 * T2 * T6;
      dIgs_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
      dIgs_dVs = -dIgs_dVg;


      T0 = vgd - paramPtr->vfbsd;
      vgd_eff = sqrt(T0 * T0 + 1.0e-4);
      dvgd_eff_dvg = T0 / vgd_eff;

      T2 = vgd * vgd_eff;
      dT2_dVg = vgd * dvgd_eff_dvg + vgd_eff;
      T5 = T12 * (paramPtr->aigsd + T3 * vgd_eff
                                - T4 * vgd_eff * vgd_eff);
      if (T5 > EXPL_THRESHOLD)
      {   T6 = MAX_EXPL;
          dT6_dVg = 0.0;
      }
      else if (T5 < -EXPL_THRESHOLD)
      {   T6 = MIN_EXPL;
          dT6_dVg = 0.0;
      }
      else
      {   T6 = exp(T5);
          dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgd_eff)
                  * dvgd_eff_dvg;
      }
      Igd_local = T11 * T2 * T6;
      dIgd_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
      dIgd_dVd = -dIgd_dVg;

      Igs = Igs_local;
      gIgsg = dIgs_dVg;
      gIgss = dIgs_dVs;
      Igd = Igd_local;
      gIgdg = dIgd_dVg;
      gIgdd = dIgd_dVd;
  }
  else
  {   Igcs = gIgcsg = gIgcsd = gIgcsb = 0.0;
      Igcd = gIgcdg = gIgcdd = gIgcdb = 0.0;
      Igs = gIgsg = gIgss = 0.0;
      Igd = gIgdg = gIgdd = 0.0;
  }

  gIgcss = -(gIgcsg + gIgcsd + gIgcsb);
  gIgcds = -(gIgcdg + gIgcdd + gIgcdb);


  // gate-body tunneling component
  if ((model_.igbMod!= 0) && (soiMod != 2))  // v3.2
  // v3.1: the Igb calculation is skipped for the ideal FD mode
  {
     OxideRatio = paramPtr->oxideRatio;

     Vox = Voxdepinv;
     // Voxeff is Vox limited below Voxh
     T0 = model_.voxh;
     T1 = T0 - Vox - model_.deltavox;
     T3 = sqrt(T1 * T1 + 4*model_.deltavox * T0);
     Voxeff = T0 - 0.5 * (T1 + T3);
     dVoxeff_dVox = 0.5 * (1.0 + T1 / T3);

     Vox = Voxeff;
     dVox_dVg = dVoxdepinv_dVg * dVoxeff_dVox;
     dVox_dVd = dVoxdepinv_dVd * dVoxeff_dVox;
     dVox_dVb = dVoxdepinv_dVb * dVoxeff_dVox;
     dVox_dVe = dVoxdepinv_dVe * dVoxeff_dVox; // v3.0
     dVox_dT = dVoxdepinv_dT * dVoxeff_dVox;


     T0 = (Vox - model_.ebg)/model_.vevb;
     if (selfheat)
        dT0_dT = dVox_dT /model_.vevb;

     CEXP(T0, T1, T2); // T1=exp(T0), T2=dT1_dT0
     if (selfheat)
        dT1_dT = T2 * dT0_dT;

     Vaux = model_.vevb * log(1 + T1);
     dVaux_dVg = T2 / (1 + T1) * dVox_dVg;
     dVaux_dVd = T2 / (1 + T1) * dVox_dVd;
     dVaux_dVb = T2 / (1 + T1) * dVox_dVb;
     dVaux_dVe = T2 / (1 + T1) * dVox_dVe; // v3.0
     if (selfheat)
        dVaux_dT = T2 / (1 + T1) * dVox_dT;

     if (model_.vgb1 != 0) {
        T0 = 1 - Vox / model_.vgb1;
        dT0_dVox = -1.0/model_.vgb1;
        if (selfheat)
           dT0_dT = -dVox_dT / model_.vgb1;
     } else {
          T0 = 1;
          dT0_dVox = dT0_dT = 0.0;
       }

     if (T0 < 0.01) {
        T0 = 0.01;
        dT0_dVox = dT0_dT = 0.0;
     }

     // v2.2.3 bug fix
     T1 = Leff * Weff * 3.7622e-7 * OxideRatio / nseg;

     T2 = -3.1051e10 * model_.toxqm;
     T3 = paramPtr->alphaGB1;
     T4 = paramPtr->betaGB1;

     T6 = T2*(T3 - T4 * Vox) / T0;
     if (selfheat) dT6_dT = -T2 * T4 * dVox_dT / T0 - T6/T0 * dT0_dT;

     CEXP(T6, T5, T7); /* T5=exp(T6), T7=dT5_dT6 */
     dT5_dVg = -T7 * dVox_dVg * T2 / T0 * (T4 + (T3 - T4 * Vox)
                                   / T0 * dT0_dVox);
     dT5_dVd = -T7 * dVox_dVd * T2 / T0 * (T4 + (T3 - T4 * Vox)
                                   / T0 * dT0_dVox);
     dT5_dVb = -T7 * dVox_dVb * T2 / T0 * (T4 + (T3 - T4 * Vox)
                                   / T0 * dT0_dVox);
     dT5_dVe = -T7 * dVox_dVe * T2 / T0 * (T4 + (T3 - T4 * Vox)
                                   / T0 * dT0_dVox); // v3.0
     if (selfheat)
        dT5_dT = T7 * dT6_dT;

     Igb1 = T1 * Vgb * Vaux * T5;
     dIgb1_dVg = T1 * (Vgb*Vaux*dT5_dVg + dVgb_dVg*Vaux*T5 +
                 Vgb*T5*dVaux_dVg);
     dIgb1_dVd = T1 * (Vgb*Vaux*dT5_dVd + Vgb*T5*dVaux_dVd);
     dIgb1_dVb = T1 * (Vgb*Vaux*dT5_dVb + dVgb_dVb*Vaux*T5 +
                 Vgb*T5*dVaux_dVb);
     dIgb1_dVe = T1 * (Vgb*Vaux*dT5_dVe + Vgb*T5*dVaux_dVe); // v3.0
     if (selfheat)
        dIgb1_dT = T1 * Vgb * (Vaux*dT5_dT + T5*dVaux_dT);
     else dIgb1_dT = 0.0;


     Vox = Voxacc;
     // Voxeff is Vox limited below Voxh
     T0 = model_.voxh;
     T1 = T0 - Vox - model_.deltavox;
     T3 = sqrt(T1 * T1 + 4*model_.deltavox * T0);
     Voxeff = T0 - 0.5 * (T1 + T3);
     dVoxeff_dVox = 0.5 * (1.0 + T1 / T3);

     Vox = Voxeff;
     dVox_dVg = dVoxacc_dVg * dVoxeff_dVox;
     dVox_dVd = dVoxacc_dVd * dVoxeff_dVox;
     dVox_dVb = dVoxacc_dVb * dVoxeff_dVox;
     dVox_dT = 0;

     T0 = (-Vgb+(Vfb))/model_.vecb;
     if (selfheat)
        dT0_dT = 0;

     CEXP(T0, T1, T2); // T1=exp(T0), T2=dT1_dT0
     if (selfheat)
        dT1_dT = 0;

     Vaux = model_.vecb* log(1 + T1);
     dVaux_dVg = -T2 / (1 + T1);
     dVaux_dVd = 0;
     dVaux_dVb = -dVaux_dVg;
     if (selfheat)
        dVaux_dT = 0;

     if (model_.vgb2 != 0) {
        T0 = 1 - Vox / model_.vgb2;
        dT0_dVox = -1.0/model_.vgb2;
        if (selfheat) dT0_dT = -dVox_dT / model_.vgb2;
     } else {
          T0 = 1;
          dT0_dVox = dT0_dT =0.0;
       }

     if (T0 < 0.01) {
        T0 = 0.01;
        dT0_dVox = dT0_dT =0.0;
     }

     // v2.2.3 bug fix
     T1 = Leff * Weff * 4.9758e-7  * OxideRatio / nseg;

     T2 = -2.357e10 * model_.toxqm;
     T3 = paramPtr->alphaGB2;
     T4 = paramPtr->betaGB2;

     T6 = T2*(T3 - T4 * Vox) / T0;
     if (selfheat) dT6_dT = -T2 * T4 * dVox_dT / T0 - T6/T0 * dT0_dT;

     CEXP(T6, T5, T7); /* T5=exp(T6), T7=dT5_dT6 */
     dT5_dVg = -T7 * dVox_dVg * T2 / T0 * (T4 + (T3 - T4 * Vox)
                                   / T0 * dT0_dVox);
     dT5_dVd = -T7 * dVox_dVd * T2 / T0 * (T4 + (T3 - T4 * Vox)
                                   / T0 * dT0_dVox);
     dT5_dVb = -T7 * dVox_dVb * T2 / T0 * (T4 + (T3 - T4 * Vox)
                                   / T0 * dT0_dVox);
     if (selfheat)
        dT5_dT = T7 * dT6_dT;

     Igb2 = T1 * Vgb * Vaux * T5;
     dIgb2_dVg = T1 * (Vgb*Vaux*dT5_dVg + dVgb_dVg*Vaux*T5 +
                 Vgb*T5*dVaux_dVg);
     dIgb2_dVd = T1 * (Vgb*Vaux*dT5_dVd + Vgb*T5*dVaux_dVd);
     dIgb2_dVb = T1 * (Vgb*Vaux*dT5_dVb + dVgb_dVb*Vaux*T5 +
                 Vgb*T5*dVaux_dVb);
     if (selfheat)
        dIgb2_dT = T1 * Vgb * (Vaux*dT5_dT + T5*dVaux_dT);
     else dIgb2_dT = 0.0;


     // Igb1 dominates in inversion region,
     // while Igb2 doninates in accumulation
     // v2.2.3 bug fix for residue at low Vgb
     if (Vgb >= 0)
     {
        Igb = Igb1;
        dIgb_dVg = dIgb1_dVg;
        dIgb_dVd = dIgb1_dVd;
        dIgb_dVb = dIgb1_dVb;
        dIgb_dVe = dIgb1_dVe; // v3.0
        dIgb_dT = dIgb1_dT;
     }
     else
     {
        Igb = Igb2;
        dIgb_dVg = dIgb2_dVg;
        dIgb_dVd = dIgb2_dVd;
        dIgb_dVb = dIgb2_dVb;
        dIgb_dVe = 0; // v3.0
        dIgb_dT = dIgb2_dT;
     }

  }
  else {
     Igb = 0.0;
     dIgb_dVg = 0.0;
     dIgb_dVd = 0.0;
     dIgb_dVb = 0.0;
     dIgb_dVe = 0.0; // v3.0
     dIgb_dT = 0.0;
  }

  ig = Igb;
  gigg = dIgb_dVg;
  gigd = dIgb_dVd;
  gigb = dIgb_dVb;
  gige = dIgb_dVe; // v3.0
  gigs = -(dIgb_dVg + dIgb_dVd + dIgb_dVb + dIgb_dVe);
  gigT = dIgb_dT;
  // end of gate-body tunneling
  // end of v3.0 gate-tunneling

  // v3.1

  if (soiMod != 2)  // v3.2
  {
    // calculate substrate current Iii
    if (paramPtr->alpha0 <= 0.0) {
       Giig = Giib = Giid = GiiT = 0.0;
       Giie = 0; // v3.0
       iii = Iii = 0.0;
    }
    else {
      Vdsatii0 = paramPtr->vdsatii0 * (1 + model_.tii * (TempRatio-1.0))
         - paramPtr->lii / Leff;
      if (selfheat)
         dVdsatii0_dT = paramPtr->vdsatii0 * model_.tii * dTempRatio_dT;
      else
         dVdsatii0_dT = 0;

      // Calculate VgsStep
      T0 = paramPtr->esatii * Leff;
               // v3.0 bug fix: T0 is dimentionless (i.e., scaled by 1V)
      T1 = paramPtr->sii0 * T0 / (1.0 + T0);

      T0 = 1 / (1 + paramPtr->sii1 * Vgsteff_local);
      if (selfheat)
         dT0_dT = - paramPtr->sii1 * T0 * T0 *dVgsteff_dT;
      else
         dT0_dT = 0;
      T3 = T0 + paramPtr->sii2;
      T4 = Vgst * paramPtr->sii1 * T0 * T0;
      T2 = Vgst * T3;
      dT2_dVg = T3 * (dVgst_dVg - dVth_dVb * dVbseff_dVg) - T4
                                         * dVgsteff_dVg; // v3.0
      dT2_dVb = T3 * dVgst_dVb * dVbseff_dVb - T4 * dVgsteff_dVb;
      dT2_dVe = T3 * dVgst_dVb * dVbseff_dVe - T4 * dVgsteff_dVe; // v3.0
      dT2_dVd = T3 * (dVgst_dVd - dVth_dVb * dVbseff_dVd) - T4
                                                  * dVgsteff_dVd; // v3.0
      if (selfheat)
         dT2_dT = -(dVth_dT + dVth_dVb * dVbseff_dT) * T3
                                             + Vgst * dT0_dT;     // v3.0
      else dT2_dT = 0;


      T3 = 1 / (1 + paramPtr->siid * Vds);
      dT3_dVd = - paramPtr->siid * T3 * T3;

      VgsStep = T1 * T2 * T3;
      if (selfheat)
         dVgsStep_dT = T1 * T3 * dT2_dT;
      else dVgsStep_dT = 0;
      Vdsatii = Vdsatii0 + VgsStep;
      Vdiff = Vds - Vdsatii;
      dVdiff_dVg = - T1 * T3 * dT2_dVg;
      dVdiff_dVb = - T1 * T3 * dT2_dVb;
      dVdiff_dVe = - T1 * T3 * dT2_dVe; // v3.0
      dVdiff_dVd = 1.0 - T1 * (T3 * dT2_dVd + T2 * dT3_dVd);
      if (selfheat)
         dVdiff_dT  = -(dVdsatii0_dT + dVgsStep_dT);
      else dVdiff_dT = 0;

      T0 = paramPtr->beta2 + paramPtr->beta1 * Vdiff
         + paramPtr->beta0 * Vdiff * Vdiff;
      if (T0 < 1e-5)
      {
         T0 = 1e-5;
         dT0_dVg = dT0_dVd = dT0_dVb = dT0_dT = 0.0;
         dT0_dVe = 0; /* v3.0 */
      }
      else
      {
         T1 = paramPtr->beta1 + 2 * paramPtr->beta0 * Vdiff;
         dT0_dVg = T1 * dVdiff_dVg;
         dT0_dVb = T1 * dVdiff_dVb;
         dT0_dVd = T1 * dVdiff_dVd;
         dT0_dVe = T1 * dVdiff_dVe; // v3.0
         if (selfheat)
            dT0_dT = T1 * dVdiff_dT;
         else
            dT0_dT = 0;
      }

      if ((T0 < Vdiff / EXPL_THRESHOLD) && (Vdiff > 0.0)) {
         Ratio = paramPtr->alpha0 * MAX_EXPL;
         dRatio_dVg = dRatio_dVb = dRatio_dVd = dRatio_dT = 0.0;
         dRatio_dVe = 0; /* v3.0 */
      }
      else if ((T0 < -Vdiff / EXPL_THRESHOLD) && (Vdiff < 0.0)) {
         Ratio = paramPtr->alpha0 * MIN_EXPL;
         dRatio_dVg = dRatio_dVb = dRatio_dVd = dRatio_dT = 0.0;
         dRatio_dVe = 0; /* v3.0 */
      }
      else {
         Ratio = paramPtr->alpha0 * exp(Vdiff / T0);
         T1 = Ratio / T0 / T0;
         dRatio_dVg = T1 * (T0 * dVdiff_dVg - Vdiff * dT0_dVg);
         dRatio_dVb = T1 * (T0 * dVdiff_dVb - Vdiff * dT0_dVb);
         dRatio_dVd = T1 * (T0 * dVdiff_dVd - Vdiff * dT0_dVd);
         /* v3.0 */
         dRatio_dVe = T1 * (T0 * dVdiff_dVe - Vdiff * dT0_dVe);

         if (selfheat)
            dRatio_dT = T1 * (T0 * dVdiff_dT - Vdiff * dT0_dT);
         else
            dRatio_dT = 0;
      }

      // Avoid too high ratio
      if (Ratio > 10.0) {
         Ratio = 10.0;
         dRatio_dVg = dRatio_dVb = dRatio_dVd = dRatio_dT = 0.0;
         dRatio_dVe = 0; /* v3.0 */
      }

      T0 = Ids + paramPtr->fbjtii * Ic;
      iii = Iii = Ratio * T0;
      Giig = Ratio * Gm + T0 * dRatio_dVg;
      Giib = Ratio * (Gmb + paramPtr->fbjtii * Gcb)
           + T0 * dRatio_dVb;
      Giid = Ratio * (Gds + paramPtr->fbjtii * Gcd)
           + T0 * dRatio_dVd;
      // v3.0
      Giie = Ratio * Gme + T0 * dRatio_dVe;

      if (selfheat)
         GiiT = Ratio * (GmT + paramPtr->fbjtii * GcT)
              + T0 * dRatio_dT;
      else
         GiiT = 0.0;
    }   // alpha0 <= 0


    // Current through body resistor
    // Current going out is +ve
    if ((bodyMod == 0) || (bodyMod == 2))
    {
       Ibp = Gbpbs = Gbpps = 0.0;
    }
    else
    {   // bodyMod == 1
       if (paramPtr->rbody < 1e-3)         // 3.2 bug fix
       {
          if (rbodyext <= 1e-3) // 3.2 bug fix
             T0 = 1.0 / 1e-3;              // 3.2 bug fix
          else
             T0 = 1.0 / rbodyext;
          Ibp = Vbp * T0;
          Gbpbs = T0 * dVbp_dVb;
          Gbpps = -T0 * dVbp_dVb;
       }
       else
       {
           Gbpbs = 1.0 / (paramPtr->rbody + rbodyext);
           Ibp = Vbp * Gbpbs;
           Gbpps = - Gbpbs;
       }

    }  //  soiMod != 2

    ibp = Ibp;
    gbpbs = Gbpbs;
    gbpps = Gbpps;
    gbpT = 0.0;

//  spice has:
//  cbodcon = Ibp - (Gbpbs * Vbs + Gbpps * Vps);

    cbodcon = Ibp;
    cbodcon_Jdxp = 0.0;
    if (devOptions.voltageLimiterFlag && !origFlag)
    {
      tmp = - (Gbpbs * (Vbs-Vbs_orig) + Gbpps * (Vps-Vps_orig));
      cbodcon_Jdxp = tmp;
    }
  }
  else // v3.1 soiMod=2: ideal FD
  {
    Giig = Giib = Giid = Giie = GiiT = 0.0;
    iii = Iii = 0.0;

    ibp = Ibp = 0.0;
    gbpbs = 0.0;
    gbpps = gbpT = cbodcon = 0.0;
    Gbpbs = Gbpps = 0.0;
  }
  // v3.1

  //  Current going out of drainprime node into the drain of device
  //  "node" means the SPICE circuit node


  cdrain = Ids + Ic;
  cd = Ids + Ic - Ibd + Iii + Idgidl;
  cb = Ibs + Ibd + Ibp - Iii - Idgidl - Isgidl - Igb;
  gds = Gds + Gcd;
  gm = Gm;
  gmbs = Gmb + Gcb;
  // v3.0
  gme = Gme;

  // v3.1 wanh added Rg for RF
  // Calculate Rg
  if (rgateMod >1)
  {  T9 = paramPtr->xrcrg2 * vtm;
     T0 = T9 *beta;
     dT0_dVd = (dbeta_dVd + dbeta_dVg * dVgsteff_dVd) * T9;
     dT0_dVb = (dbeta_dVb + dbeta_dVg * dVgsteff_dVb) * T9;
     dT0_dVg = dbeta_dVg * T9;

     gcrg = paramPtr->xrcrg1 * (T0 + Ids);
     gcrgd = paramPtr->xrcrg1 * (dT0_dVd +tmp1);
     gcrgb = paramPtr->xrcrg1 * (dT0_dVb +tmp2)
                      * dVbseff_dVb;
     gcrgg = paramPtr->xrcrg1 * (dT0_dVg +tmp3)
                      * dVgsteff_dVg;

     if (rgateMod == 2)
     {   T10 = grgeltd * grgeltd;
         T11 = grgeltd + gcrg;
         gcrg = grgeltd * gcrg / T11;
         T12 = T10 / T11 /T11;
         gcrgg *= T12;
         gcrgd *= T12;
         gcrgb *= T12;
     }

     gcrgs = -(gcrgg + gcrgd + gcrgb);
  }
  // v3.1 wanh added Rg for RF end

  if (selfheat)
     gmT = GmT + GcT;
  else
     gmT = 0.0;

  //  note that sign is switched because power flows out
  //  of device into the temperature node.
  //  Currently ommit self-heating due to bipolar current
  //  because it can cause convergence problem

  gtempg = -Gm  * Vds;
  gtempb = -Gmb * Vds;
  // v3.0
  gtempe = -Gme * Vds;

  gtempT = -GmT * Vds;
  gtempd = -Gds * Vds - Ids;

//  spice has:
//  cth = - Ids * Vds - model_.dtype * (gtempg * Vgs + gtempb * Vbs
//                  + gtempe * Ves + gtempd * Vds) - gtempT * delTemp; // v3.0

  cth = - Ids * Vds;
  cth_Jdxp = 0.0;
  if (devOptions.voltageLimiterFlag && !origFlag)
  {
    tmp = -model_.dtype *(gtempg * (Vgs-Vgs_orig)
                                +gtempb * (Vbs-Vbs_orig)
                                +gtempe * (Ves-Ves_orig)
                                +gtempd * (Vds-Vds_orig))
      -gtempT * (delTemp-delTemp_orig);
    cth_Jdxp = tmp;
  }


  //  Body current which flows into drainprime node from the drain of device

  gjdb = Gjdb - Giib;
  gjdd = Gjdd - (Giid + Gdgidld);
  gjdg = - (Giig + Gdgidlg);
  // v3.0
  gjde = - Giie;

  if (selfheat) gjdT = GjdT - GiiT;
  else gjdT = 0.0;

  cjd = Ibd - Iii - Idgidl;
  cjd_Jdxp = 0.0;
  if (devOptions.voltageLimiterFlag && !origFlag)
  {
    tmp = -(gjdb*(Vbs-Vbs_orig)
            + gjdd*(Vds-Vds_orig)
            + gjdg*(Vgs-Vgs_orig)
            + gjde*(Ves-Ves_orig)
            + gjdT*(delTemp-delTemp_orig));
    cjd_Jdxp = tmp;
  }

  //  Body current which flows into sourceprime node from the source of device
  gjsb = Gjsb;
  gjsd = Gjsd;
  gjsg = - Gsgidlg;
  if (selfheat) gjsT = GjsT;
  else gjsT = 0.0;

  cjs = Ibs - Isgidl;
  cjs_Jdxp = 0.0;
  if (devOptions.voltageLimiterFlag && !origFlag)
  {
    tmp = -(gjsb*(Vbs-Vbs_orig)
           + gjsd*(Vds-Vds_orig)
           + gjsg*(Vgs-Vgs_orig)
           + gjsT*(delTemp-delTemp_orig));
    cjs_Jdxp = tmp;
  }

  //  Current flowing into body node

  gbbs = Giib - Gjsb - Gjdb - Gbpbs;
  gbgs = Giig + Gdgidlg + Gsgidlg;
  gbds = Giid + Gdgidld - Gjsd - Gjdd;
  // v3.0
  gbes = Giie;

  gbps = - Gbpps;
  if (selfheat) gbT = GiiT - GjsT - GjdT;
  else gbT = 0.0;

//  spice has:
//  cbody = Iii + Idgidl + Isgidl - Ibs - Ibd - Ibp + Igb
//                   - ( (gbbs + dIgb_dVb) * Vbs
//                   + (gbgs + dIgb_dVg) * Vgs
//                   + (gbds + dIgb_dVd) * Vds
//                   + gbps * Vps
//                   + (gbes + dIgb_dVe) * Ves
//                   + (gbT + dIgb_dT) * delTemp); // v3.0
//
  cbody_Jdxp = 0.0;
  cbody = Iii + Idgidl + Isgidl - Ibs - Ibd - Ibp + Igb;
  if (devOptions.voltageLimiterFlag && !origFlag)
  {
    tmp = -( (gbbs + dIgb_dVb) * (Vbs - Vbs_orig)
            +(gbgs + dIgb_dVg) * (Vgs - Vgs_orig)
            +(gbds + dIgb_dVd) * (Vds - Vds_orig)
            +gbps              * (Vps - Vps_orig)
            +(gbes + dIgb_dVe) * (Ves - Ves_orig)
            +(gbT + dIgb_dT)   * (delTemp - delTemp_orig));
    cbody_Jdxp  = tmp;
  }

//  spice has:
//  cgate = Igb - (dIgb_dVb * Vbs + dIgb_dVe * Ves + dIgb_dVg
//                   * Vgs + dIgb_dVd * Vds + dIgb_dT * delTemp); // v3.0

  cgate_Jdxp = 0.0;
  cgate = Igb;
  if (devOptions.voltageLimiterFlag && !origFlag)
  {
    tmp = -(dIgb_dVb * (Vbs - Vbs_orig) +
            dIgb_dVe * (Ves - Ves_orig) +
            dIgb_dVg * (Vgs - Vgs_orig) +
            dIgb_dVd * (Vds - Vds_orig) +
            dIgb_dT  * (delTemp - delTemp_orig));
    cgate_Jdxp  = tmp;
  }

  // Calculate Qinv for Noise analysis

  T1 = Vgsteff_local * (1.0 - 0.5 * Abulk_local * Vdseff_local / Vgst2Vtm);
  qinv = -model_.cox * paramPtr->weff * Leff * T1;

  //  Begin CV (charge) model

  if ((model_.xpart < 0) || (!ChargeComputationNeeded))
  {    qgate  = qdrn = qsrc = qbody = qsub = 0.0; // v2.2.3 bug fix
       cggb = cgsb = cgdb = 0.0;
       cdgb = cdsb = cddb = 0.0;
       cbgb = cbsb = cbdb = 0.0;
       Qsub0 = 0;
       Qac0 = 0;
       qjs_local = 0;
       qjd_local = 0;
       goto finished;
  }
  else
  {
       CoxWL  = model_.cox * (paramPtr->weffCV / nseg
                       * paramPtr->leffCV + agbcp);
       CoxWLb = model_.fbody * model_.cox
                       * (paramPtr->weffCV / nseg
                       * paramPtr->leffCVb + agbcp);


       // v3.2 Seperate VgsteffCV with noff
       noff = n * paramPtr->noff;
       dnoff_dVd = paramPtr->noff * dn_dVd;
       dnoff_dVb = paramPtr->noff * dn_dVb;

       if ((VgstNVt > -EXPL_THRESHOLD) && (VgstNVt < EXPL_THRESHOLD))
       {
           ExpVgst *= ExpVgst;
           ExpVgst *= exp( -(paramPtr->delvt / (noff * Vtm)));
           Vgsteff_local = noff * Vtm * log(1.0 + ExpVgst);
           T0 = ExpVgst / (1.0 + ExpVgst);
           T1 = -T0 * (dVth_dVb + (Vgst-paramPtr->delvt) / noff * dnoff_dVb)
                                + Vgsteff_local / noff * dnoff_dVb;
           dVgsteff_dVd = -T0 * (dVth_dVd + dVth_dVb*dVbseff_dVd + Vgst / noff
                          * dnoff_dVd) + Vgsteff_local / noff * dnoff_dVd;
           dVgsteff_dVg = T0 * (dVgs_eff_dVg - dVth_dVb*dVbseff_dVg);
                           dVgsteff_dVb = T1 * dVbseff_dVb;
                           dVgsteff_dVe = T1 * dVbseff_dVe;
           if (selfheat)
               dVgsteff_dT = -T0 * (dVth_dT+dVth_dVb*dVbseff_dT
                       + (Vgst - paramPtr->delvt) / Temp)
                                               + Vgsteff_local / Temp;
            else
               dVgsteff_dT  = 0.0;
       }
       // v3.2


       if (model_.capMod == 2)
       {

          // v3.1
         if (soiMod == 2) // v3.2 - ideal FD
         {
             Qac0 = dQac0_dVrg = dQac0_dVd = dQac0_dVb = dQac0_dT = 0.0;
             Qsub0 = dQsub0_dVrg = dQsub0_dVg = dQsub0_dVd
                                        = dQsub0_dVb = dQsub0_dT = 0.0;
         }
         else // soiMod = 0 or 1
         {
             Vfb = Vth - phi_local - paramPtr->k1eff * sqrtPhis
                                          + paramPtr->delvt;
             dVfb_dVb = dVth_dVb - paramPtr->k1eff * dsqrtPhis_dVb;
             dVfb_dVd = dVth_dVd;
             dVfb_dT  = dVth_dT;

             V3 = Vfb - Vgs_eff + Vbseff - DELTA_3_SOI;
             if (Vfb <= 0.0)
             {   T0 = sqrt(V3 * V3 - 4.0 * DELTA_3_SOI * Vfb);
                 T2 = -DELTA_3_SOI / T0;
             }
             else
             {   T0 = sqrt(V3 * V3 + 4.0 * DELTA_3_SOI * Vfb);
                 T2 = DELTA_3_SOI / T0;
             }

             T1 = 0.5 * (1.0 + V3 / T0);
             Vfbeff = Vfb - 0.5 * (V3 + T0);
             dVfbeff_dVd = (1.0 - T1 - T2) * dVfb_dVd;
             dVfbeff_dVb = (1.0 - T1 - T2) * dVfb_dVb - T1;
             dVfbeff_dVrg = T1 * dVgs_eff_dVg;
             if (selfheat)
                dVfbeff_dT = (1.0 - T1 - T2) * dVfb_dT;
             else
                dVfbeff_dT = 0.0;

             Qac0 = CoxWLb * (Vfbeff - Vfb);
             dQac0_dVrg = CoxWLb * dVfbeff_dVrg;
             dQac0_dVd = CoxWLb * (dVfbeff_dVd - dVfb_dVd);
             dQac0_dVb = CoxWLb * (dVfbeff_dVb - dVfb_dVb);
             if (selfheat)
                dQac0_dT = CoxWLb * (dVfbeff_dT - dVfb_dT);
             else
                dQac0_dT = 0.0;

             T0 = 0.5 * K1;
             T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff_local;
             if (paramPtr->k1eff == 0.0)
             {   T1 = 0.0;
                 T2 = 0.0;
             }
             else if (T3 < 0.0)
             {   T1 = T0 + T3 / paramPtr->k1eff;
                 T2 = CoxWLb;
             }
             else
             {   T1 = sqrt(T0 * T0 + T3);
                 T2 = CoxWLb * T0 / T1;
             }

             Qsub0 = CoxWLb * K1 * (T1 - T0);
             dQsub0_dVrg = T2 * (dVgs_eff_dVg - dVfbeff_dVrg);
             dQsub0_dVg = -T2;
             dQsub0_dVd = -T2 * dVfbeff_dVd;
             dQsub0_dVb = -T2 * (dVfbeff_dVb + 1);
             if (selfheat)
                dQsub0_dT  = -T2 * dVfbeff_dT;
             else
                dQsub0_dT = 0.0;
         }
         // v3.1


         AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
         dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;

         VdsatCV = Vgsteff_local / AbulkCV;
//not used         dVdsatCV_dVg = 1.0 / AbulkCV;
//not used         dVdsatCV_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;

         V4 = VdsatCV - Vds - DELTA_4;
         T0 = sqrt(V4 * V4 + 4.0 * DELTA_4 * VdsatCV);
         VdseffCV = VdsatCV - 0.5 * (V4 + T0);
         T1 = 0.5 * (1.0 + V4 / T0);
         T2 = DELTA_4 / T0;
         T3 = (1.0 - T1 - T2) / AbulkCV;
         dVdseffCV_dVg = T3;
         dVdseffCV_dVd = T1;
         dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;


         // v3.1
         if (soiMod == 2) // v3.2 - ideal FD
         {
            qbulk = Cbg1 = Cbd1 = Cbb1 = 0;
         }
         else
         {
             T0 = AbulkCV * VdseffCV;
             T1 = 12.0 * (Vgsteff_local - 0.5 * T0 + 1e-20);
             T2 = VdseffCV / T1;
             T3 = T0 * T2;
             T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
             T5 = (6.0 * T0 * (4.0 * Vgsteff_local- T0) / (T1 * T1) - 0.5);
             T6 = 12.0 * T2 * T2 * Vgsteff_local;

             T7 = 1.0 - AbulkCV;
             qbulk = CoxWLb * T7 * (0.5 * VdseffCV - T3);
             T4 = -T7 * (T4 - 1.0);
             T5 = -T7 * T5;
             T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));

             Cbg1 = CoxWLb * (T4 + T5 * dVdseffCV_dVg);
             Cbd1 = CoxWLb * T5 * dVdseffCV_dVd ;
             Cbb1 = CoxWLb * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb);
         }
         // v3.1


         // Total inversion charge
         T0 = AbulkCV * VdseffCV;
         T1 = 12.0 * (Vgsteff_local - 0.5 * T0 + 1e-20);
//       T2 = VdseffCV / T1;
         T2 = T0 / T1;
         T3 = T0 * T2;

//       T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
//       T5 = (6.0 * T0 * (4.0 * Vgsteff_local - T0) / (T1 * T1) - 0.5);
//       T6 = 12.0 * T2 * T2 * Vgsteff_local;
         T4 = (1.0 - 12.0 * T2 * T2);      //bug fix by wanh
         T7 = T2 * (2.0 + 6.0 * T2) - 0.5; //bug fix by wanh

         T5 = T7 * AbulkCV;
         T6 = T7 * VdseffCV;

//       qinv_local = CoxWL * (Vgsteff_local - 0.5 * VdseffCV + T3);
         qinv_local = CoxWL * (Vgsteff_local - 0.5 * T0 + T3);

         qinv = -qinv_local; // for noise v3.2

         Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
         Cgd1 = CoxWL * T5 * dVdseffCV_dVd;
         Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb);

         // Inversion charge partitioning into S / D
         if (model_.xpart > 0.5)
         {   // 0/100 Charge partition model
             T1 = T1 + T1;
             qsrc = -CoxWL * (0.5 * Vgsteff_local + 0.25 * T0
                  - T0 * T0 / T1);
             T7 = (4.0 * Vgsteff_local - T0) / (T1 * T1);
             T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
             T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
             T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
             Csg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
             Csd1 = CoxWL * T5 * dVdseffCV_dVd;
             Csb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb);

         }
         else if (model_.xpart < 0.5)
         {   // 40/60 Charge partition model
             T1 = T1 / 12.0;
             T2 = 0.5 * CoxWL / (T1 * T1);
             T3 = Vgsteff_local * (2.0 * T0 * T0 / 3.0 + Vgsteff_local
                  * (Vgsteff_local - 4.0 * T0 / 3.0))
                         - 2.0 * T0 * T0 * T0 / 15.0;
             qsrc = -T2 * T3;
             T7 = 4.0 / 3.0 * Vgsteff_local * (Vgsteff_local - T0)
                           + 0.4 * T0 * T0;
             T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff_local * (3.0
                  * Vgsteff_local - 8.0 * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
             T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
             T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
             Csg1 = T4 + T5 * dVdseffCV_dVg;
             Csd1 = T5 * dVdseffCV_dVd;
             Csb1 = T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb;
         }
         else
         {   // 50/50 Charge partition model
             qsrc = - 0.5 * (qinv_local + qbulk);
             Csg1 = - 0.5 * (Cgg1 + Cbg1);
             Csb1 = - 0.5 * (Cgb1 + Cbb1);
             Csd1 = - 0.5 * (Cgd1 + Cbd1);
         }



         // Backgate charge
         // v3.1
         if (soiMod == 2) // v3.2 - ideal FD
         {
            Qe1 = dQe1_dVb = dQe1_dVe = dQe1_dT = 0;
         }
         else // soiMod = 0 or 1
         {
             CboxWL = paramPtr->kb1 * model_.fbody * Cbox
                                    * (paramPtr->weffCV / nseg
                                    * paramPtr->leffCVbg + aebcp);
             Qe1 = CboxWL * (Vesfb - Vbs);
             dQe1_dVb = -CboxWL;
             dQe1_dVe = CboxWL;
             if (selfheat)
                dQe1_dT = -CboxWL * dvfbb_dT;
             else
                dQe1_dT = 0;
         }
         // v3.1


         qgate = qinv_local + Qac0 + Qsub0;
         qbody = (qbulk - Qac0 - Qsub0 - Qe1);
         qsub = Qe1;
         qdrn = -(qgate + qsrc + qbody + qsub);

         // This transform all the dependency on Vgsteff_local, Vbseff
         //    into real ones
         Ce1b = dQe1_dVb;
         Ce1e = dQe1_dVe;

         Csg = Csg1 * dVgsteff_dVg;
         Csd = Csd1 + Csg1 * dVgsteff_dVd;
         Csb = Csg1 * dVgsteff_dVb + Csb1 * dVbseff_dVb;
         if (selfheat)
            CsT = Csg1 * dVgsteff_dT;
         else
            CsT = 0.0;

         Cgg = (Cgg1 + dQsub0_dVg) * dVgsteff_dVg
                                + dQac0_dVrg + dQsub0_dVrg;
         Cgd = (Cgg1 + dQsub0_dVg) * dVgsteff_dVd + Cgd1
                                + dQac0_dVd + dQsub0_dVd;
         Cgb = (Cgg1 + dQsub0_dVg) * dVgsteff_dVb
                          + (Cgb1 + dQsub0_dVb + dQac0_dVb) * dVbseff_dVb;
         if (selfheat)
            CgT = (Cgg1 + dQsub0_dVg) * dVgsteff_dT + dQac0_dT + dQsub0_dT;
         else
            CgT = 0.0;

         Cbg = (Cbg1 - dQsub0_dVg) * dVgsteff_dVg - dQac0_dVrg - dQsub0_dVrg;
         Cbd = (Cbg1 - dQsub0_dVg) * dVgsteff_dVd + Cbd1 - dQac0_dVd
                                                    - dQsub0_dVd;
         Cbb = (Cbg1 - dQsub0_dVg) * dVgsteff_dVb - dQe1_dVb
                          + (Cbb1 - dQsub0_dVb - dQac0_dVb) * dVbseff_dVb;
         if (selfheat)
            CbT = (Cbg1 - dQsub0_dVg) * dVgsteff_dT - dQac0_dT - dQsub0_dT
                                                   - dQe1_dT;
         else
            CbT = 0.0;


         cggb = Cgg ;
         cgsb = - (Cgg  + Cgd  + Cgb);
         cgdb = Cgd;
         cgT = CgT;

         cbgb = Cbg;
         cbsb = -(Cbg  + Cbd  + Cbb) + Ce1e;
         cbdb = Cbd;
         cbeb = - Ce1e ;
         cbT = CbT;

         ceeb = Ce1e ;
         ceT = dQe1_dT;

         cdgb = -(Cgg + Cbg + Csg);
         cddb = -(Cgd + Cbd + Csd);
         cdeb = 0;
         cdT = -(CgT + CbT + CsT) - dQe1_dT;
         cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                          + Csg + Csd + Csb) + Ce1b;
     } // End of if capMod == 2

     else if (model_.capMod == 3)
     {

        dVgsteff_dVb /= dVbseff_dVb;
        Cox = 3.453133e-11 / (model_.tox - model_.dtoxcv); //v2.2.3
        CoxWL *= model_.tox/ (model_.tox - model_.dtoxcv);
        CoxWLb *= model_.tox/ (model_.tox - model_.dtoxcv);
        Tox = 1.0e8 * (model_.tox - model_.dtoxcv);


        // v3.1
        if (soiMod == 2) // v3.2 - ideal FD
        {
           Qac0 = dQac0_dVg = dQac0_dVb = dQac0_dT = 0.0;
           Qsub0 = dQsub0_dVg = dQsub0_dVd = dQsub0_dVb = dQsub0_dT = 0.0;
        }
        else // soiMod = 0 or 1
        {
           if (selfheat) {
              Vfbzb = Vthzb - phi_local - paramPtr->k1eff * sqrtPhi
                    + paramPtr->delvt;
              dVfbzb_dT = dVthzb_dT;
           }
           else {
              Vfbzb = paramPtr->vfbzb + paramPtr->delvt;
              dVfbzb_dT = 0;
           }

           V3 = Vfbzb - Vgs_eff + Vbseff - DELTA_3;
           if (Vfbzb <= 0.0)
           {   T0 = sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfbzb);
               T2 = -DELTA_3 / T0;
           }
           else
           {   T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfbzb);
               T2 = DELTA_3 / T0;
           }

           T1 = 0.5 * (1.0 + V3 / T0);
           Vfbeff = Vfbzb - 0.5 * (V3 + T0);
           dVfbeff_dVg = T1 * dVgs_eff_dVg;
           dVfbeff_dVb = -T1;
           if (selfheat)
              dVfbeff_dT = (1.0 - T1 - T2) * dVfbzb_dT;
           else
              dVfbeff_dT = 0.0;


           T0 = (Vgs_eff - Vbseff - Vfbzb) / Tox;
           dT0_dVg = dVgs_eff_dVg / Tox;
           dT0_dVb = -1.0 / Tox;

           tmp = T0 * paramPtr->acde;
           if ((-EXPL_THRESHOLD < tmp) && (tmp < EXPL_THRESHOLD))
           {   Tcen = paramPtr->ldeb * exp(tmp);
               dTcen_dVg = paramPtr->acde * Tcen;
               dTcen_dVb = dTcen_dVg * dT0_dVb;
               dTcen_dVg *= dT0_dVg;
               if (selfheat)
                  dTcen_dT = -Tcen * paramPtr->acde * dVfbzb_dT / Tox;
               else
                  dTcen_dT = 0;
           }
           else if (tmp <= -EXPL_THRESHOLD)
           {   Tcen = paramPtr->ldeb * MIN_EXPL;
               dTcen_dVg = dTcen_dVb = dTcen_dT = 0.0;
           }
           else
           {   Tcen = paramPtr->ldeb * MAX_EXPL;
               dTcen_dVg = dTcen_dVb = dTcen_dT = 0.0;
           }

           LINK = 1.0e-3 * (model_.tox - model_.dtoxcv); // v2.2.3
           V3 = paramPtr->ldeb - Tcen - LINK;
           V4 = sqrt(V3 * V3 + 4.0 * LINK * paramPtr->ldeb);
           Tcen = paramPtr->ldeb - 0.5 * (V3 + V4);
           T1 = 0.5 * (1.0 + V3 / V4);
           dTcen_dVg *= T1;
           dTcen_dVb *= T1;
           if (selfheat)
              dTcen_dT *= T1;
           else
              dTcen_dT = 0;

           Ccen = CONSTEPSSI / Tcen;
           T2 = Cox / (Cox + Ccen);
           Coxeff = T2 * Ccen;
           T3 = -Ccen / Tcen;
           dCoxeff_dVg = T2 * T2 * T3;
           dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
           dCoxeff_dVg *= dTcen_dVg;
           if (selfheat)
              dCoxeff_dT = T3 * dTcen_dT * (T2 - Coxeff / (Cox + Ccen));
           else
              dCoxeff_dT = 0;
           CoxWLcenb = CoxWLb * Coxeff / Cox;
           if (selfheat)
              dCoxWLcenb_dT = CoxWLb * dCoxeff_dT / Cox;
           else
              dCoxWLcenb_dT = 0;

           Qac0 = CoxWLcenb * (Vfbeff - Vfbzb);
           QovCox = Qac0 / Coxeff;
           dQac0_dVg = CoxWLcenb * dVfbeff_dVg
                     + QovCox * dCoxeff_dVg;
           dQac0_dVb = CoxWLcenb * dVfbeff_dVb
                     + QovCox * dCoxeff_dVb;
           if (selfheat)
              dQac0_dT = CoxWLcenb * (dVfbeff_dT - dVfbzb_dT)
                                  + dCoxWLcenb_dT * (Vfbeff - Vfbzb);
           else
              dQac0_dT = 0.0;

           T0 = 0.5 * paramPtr->k1eff;
           T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff_local;
           if (paramPtr->k1eff == 0.0)
           {   T1 = 0.0;
               T2 = 0.0;
           }
           else if (T3 < 0.0)
           {   T1 = T0 + T3 / paramPtr->k1eff;
               T2 = CoxWLcenb;
           }
           else
           {   T1 = sqrt(T0 * T0 + T3);
               T2 = CoxWLcenb * T0 / T1;
           }

           Qsub0 = CoxWLcenb * paramPtr->k1eff * (T1 - T0);
           QovCox = Qsub0 / Coxeff;
           dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                      + QovCox * dCoxeff_dVg;
           dQsub0_dVd = -T2 * dVgsteff_dVd;
           dQsub0_dVb = -T2 * (dVfbeff_dVb + 1 + dVgsteff_dVb)
                      + QovCox * dCoxeff_dVb;
           if (selfheat)
              dQsub0_dT = -T2 * (dVfbeff_dT + dVgsteff_dT)
                        + dCoxWLcenb_dT * paramPtr->k1eff * (T1 - T0);
           else
              dQsub0_dT = 0.0;
        }
        // v3.1


        // Gate-bias dependent delta Phis begins
        if (paramPtr->k1eff <= 0.0)
        {   Denomi = 0.25 * paramPtr->moin * Vtm;
            T0 = 0.5 * paramPtr->sqrtPhi;
        }
        else
        {   Denomi = paramPtr->moin * Vtm
                   * paramPtr->k1eff * paramPtr->k1eff;
            T0 = paramPtr->k1eff * paramPtr->sqrtPhi;
        }
        T1 = 2.0 * T0 + Vgsteff_local;

        DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff_local / Denomi);
        dDeltaPhi_dVg = 2.0 * Vtm * (T1 -T0) / (Denomi + T1 * Vgsteff_local);
        dDeltaPhi_dVd = dDeltaPhi_dVg * dVgsteff_dVd;
        dDeltaPhi_dVb = dDeltaPhi_dVg * dVgsteff_dVb;
        // End of delta Phis


        // v3.1.1 bug fix for discontinuity
        T3 = 4.0 * (Vth - Vfbzb - phi_local);
        T2 = sqrt(T3*T3 + 0.0001);
        T5 = 0.5 * (1 + T3/T2);
        T4 = 0.5 * (T3 + T2);

        Tox += Tox;
        T0 = (Vgsteff_local + T4) / Tox;
        tmp = exp(0.7 * log(T0));
        T1 = 1.0 + tmp;
        T2 = 0.7 * tmp / (T0 * Tox);
        Tcen = 1.9e-9 / T1;
        dTcen_dVg = -Tcen * T2 / T1;
        dTcen_dVd = dTcen_dVg * (T5 * 4.0 * dVth_dVd + dVgsteff_dVd);
        dTcen_dVb = dTcen_dVg * (T5 * 4.0 * dVth_dVb + dVgsteff_dVb);
        dTcen_dVg *= dVgsteff_dVg;
        if (selfheat)
           dTcen_dT = -Tcen * T2 / T1
                    * (T5 * 4.0 * (dVth_dT - dVfbzb_dT) + dVgsteff_dT);
        else
           dTcen_dT = 0;


        Ccen = CONSTEPSSI / Tcen;
        T0 = Cox / (Cox + Ccen);
        Coxeff = T0 * Ccen;
        T1 = -Ccen / Tcen;
        dCoxeff_dVg = T0 * T0 * T1;
        dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
        dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
        dCoxeff_dVg *= dTcen_dVg;
        if (selfheat)
           dCoxeff_dT = T1 * dTcen_dT * (T0 - Coxeff / (Cox + Ccen));
        else dCoxeff_dT = 0;
        CoxWLcen = CoxWL * Coxeff / Cox;
        CoxWLcenb = CoxWLb * Coxeff / Cox;

        AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
        dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
        VdsatCV = (Vgsteff_local - DeltaPhi) / AbulkCV;
        V4 = VdsatCV - Vds - DELTA_4;
        T0 = sqrt(V4 * V4 + 4.0 * DELTA_4 * VdsatCV);
        VdseffCV = VdsatCV - 0.5 * (V4 + T0);
        T1 = 0.5 * (1.0 + V4 / T0);
        T2 = DELTA_4 / T0;
        T3 = (1.0 - T1 - T2) / AbulkCV;
        T4 = T3 * ( 1.0 - dDeltaPhi_dVg);
        dVdseffCV_dVg = T4;
        dVdseffCV_dVd = T1;
        dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;

        T0 = AbulkCV * VdseffCV;
        T1 = Vgsteff_local - DeltaPhi;
        T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
        T3 = T0 / T2;
        T4 = 1.0 - 12.0 * T3 * T3;
        T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
        T6 = T5 * VdseffCV / AbulkCV;

        qinv_local = qgate = qinoi = CoxWLcen * (T1 - T0 * (0.5 - T3));
        QovCox = qgate / Coxeff;
        Cgg1 = CoxWLcen * (T4 * (1.0 - dDeltaPhi_dVg)
             + T5 * dVdseffCV_dVg);
        Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
             * dVgsteff_dVd + QovCox * dCoxeff_dVd;
        Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
             + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
        Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


        // v3.1
        if (soiMod == 2) // v3.2 - ideal FD
        {
           qbulk = Cbg1 = Cbd1 = Cbb1 = Cbg1 = 0;
        }
        else // soiMod = 0 or 1
        {
           T7 = 1.0 - AbulkCV;
           T8 = T2 * T2;
           T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
           T10 = T9 * (1.0 - dDeltaPhi_dVg);
           T11 = -T7 * T5 / AbulkCV;
           T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

           qbulk = CoxWLcenb * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
           QovCox = qbulk / Coxeff;
           Cbg1 = CoxWLcenb * (T10 + T11 * dVdseffCV_dVg);
           Cbd1 = CoxWLcenb * T11 * dVdseffCV_dVd + Cbg1
                * dVgsteff_dVd + QovCox * dCoxeff_dVd;
           Cbb1 = CoxWLcenb * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
                + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
           Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;
        }
        // v3.1


        if (model_.xpart > 0.5)
        {   // 0/100 partition
            qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0
                 - 0.5 * T0 * T0 / T2);
            QovCox = qsrc / Coxeff;
            T2 += T2;
            T3 = T2 * T2;
            T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
            T4 = -(0.5 + 24.0 * T0 * T0 / T3) * (1.0 - dDeltaPhi_dVg);
            T5 = T7 * AbulkCV;
            T6 = T7 * VdseffCV;

            Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
            Csd = CoxWLcen * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd
                + QovCox * dCoxeff_dVd;
            Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                + Csg * dVgsteff_dVb + QovCox * dCoxeff_dVb;
            Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
        }
        else if (model_.xpart < 0.5)
        {   // 40/60 partition
            T2 = T2 / 12.0;
            T3 = 0.5 * CoxWLcen / (T2 * T2);
            T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0
               * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
            qsrc = -T3 * T4;
            QovCox = qsrc / Coxeff;
            T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
            T5 = -2.0 * qsrc / T2 - T3 * (T1 * (3.0 * T1 - 8.0
               * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
            T6 = AbulkCV * (qsrc / T2 + T3 * T8);
            T7 = T6 * VdseffCV / AbulkCV;

            Csg = T5 * (1.0 - dDeltaPhi_dVg) + T6 * dVdseffCV_dVg;
            Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd
                + QovCox * dCoxeff_dVd;
            Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb
                + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb;
            Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
        }
        else
        {   // 50/50 partition
            qsrc = -0.5 * qgate;
            Csg = -0.5 * Cgg1;
            Csd = -0.5 * Cgd1;
            Csb = -0.5 * Cgb1;
        }


        // Backgate charge
        // v3.1
        if (soiMod == 2) // v3.2 - ideal FD
        {
           Qe1 = Ce1b = Ce1e = Ce1T = dQe1_dT = 0;
        }
        else // soiMod = 0 or 1
        {
           CboxWL = paramPtr->kb1 * model_.fbody * Cbox
                  * (paramPtr->weffCV / nseg
                  * paramPtr->leffCVbg + aebcp);
           Qe1 = CboxWL * (Vesfb - Vbs);
           Ce1b = dQe1_dVb = -CboxWL;
           Ce1e = dQe1_dVe = CboxWL;
           if (selfheat)
              Ce1T = dQe1_dT = -CboxWL * dvfbb_dT;
           else
              Ce1T = dQe1_dT = 0.0;
        }
        // v3.1


        qgate += Qac0 + Qsub0 - qbulk;
        qbody = qbulk - Qac0 - Qsub0 - Qe1;
        qsub = Qe1;
        qdrn = -(qgate + qbody + qsub + qsrc);

        Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
        Cbd = Cbd1 - dQsub0_dVd;
        Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb - Ce1b / dVbseff_dVb;
        if (selfheat)
           CbT = Cbg1 * dVgsteff_dT - dQac0_dT
               - dQsub0_dT - dQe1_dT;
        else
           CbT = 0.0;

        Cgg = Cgg1 - Cbg;
        Cgd = Cgd1 - Cbd;
        Cgb = Cgb1 - Cbb - Ce1b / dVbseff_dVb;
        if (selfheat)
           CgT = Cgg1 * dVgsteff_dT + dQac0_dT + dQsub0_dT;
        else
           CgT = 0.0;

        Cgb *= dVbseff_dVb;
        Cbb *= dVbseff_dVb;
        Csb *= dVbseff_dVb;
        if (selfheat)
           CsT = Csg * dVgsteff_dT;
        else
           CsT = 0.0;

        cggb = Cgg;
        cgsb = -(Cgg + Cgd + Cgb);
        cgdb = Cgd;
        cgT  = CgT;

        cbgb = Cbg;
        cbsb = -(Cbg + Cbd + Cbb)
                        + Ce1e;
        cbdb = Cbd;
        cbeb = -Ce1e;
        cbT  = CbT;

        ceT = Ce1T;
        ceeb = Ce1e ;

        cdgb = -(Cgg + Cbg + Csg);
        cddb = -(Cgd + Cbd + Csd);
        cdeb = 0;
        cdT   = -(CgT+CbT+CsT) - Ce1T;
        cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                        + Csg + Csd + Csb) + Ce1b;
        qinv = -qinoi;

     } // End of if capMod ==3
  }

  finished: // returning Values to Calling Routine

  //
  //  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
  //


  if (ChargeComputationNeeded)
  {
      // Intrinsic S/D junction charge
      // v3.1
      if (soiMod == 2) // v3.2 - ideal FD
      {
         qjs_local = qjd_local = 0.0;
         gcjdds = gcjdbs = gcjdT = 0.0;
         gcjsbs = gcjsT = 0.0;
      }
      else // soiMod = 0 or 1
      {
         PhiBSWG = model_.GatesidewallJctPotential;
         dPhiBSWG_dT = -model_.tpbswg;
         PhiBSWG += dPhiBSWG_dT * (Temp - model_.tnom);
         MJSWG = model_.bodyJctGateSideGradingCoeff;

         cjsbs = model_.unitLengthGateSidewallJctCap
                           * wdiosCV * model_.tsi / 1e-7;
         dcjsbs_dT = cjsbs * model_.tcjswg;
         cjsbs += dcjsbs_dT * (Temp - model_.tnom);

         cjdbs = model_.unitLengthGateSidewallJctCap
                           * wdiodCV * model_.tsi / 1e-7;
         dcjdbs_dT = cjdbs * model_.tcjswg;
         cjdbs += dcjdbs_dT * (Temp - model_.tnom);

         DioMax = 0.9 * (PhiBSWG);

         arg = 1.0 - (Vbs > DioMax ? DioMax : Vbs) / PhiBSWG;

         if (selfheat)
            darg_dT = (1 - arg) / PhiBSWG * dPhiBSWG_dT;

         if (MJSWG == 0.5) {
            dT3_dVb = 1.0 / sqrt(arg);

            if (selfheat) ddT3_dVb_dT = -0.5 * dT3_dVb / arg * darg_dT;
         }
         else {
            dT3_dVb = exp(-MJSWG * log(arg));

            if (selfheat) ddT3_dVb_dT = -MJSWG * dT3_dVb / arg * darg_dT;
         }
         T3 = (1.0 - arg * dT3_dVb) * PhiBSWG / (1.0 - MJSWG);

         if (selfheat)
            dT3_dT = (1.0 - arg * dT3_dVb) * dPhiBSWG_dT / (1.0 - MJSWG)
                    - (arg * ddT3_dVb_dT + darg_dT * dT3_dVb) * PhiBSWG
                                                / (1.0 - MJSWG);

         if (Vbs > DioMax)
            T3 += dT3_dVb * (Vbs - DioMax);

         qjs_local = cjsbs * T3 + model_.tt * Ibsdif;
         gcjsbs = cjsbs * dT3_dVb + model_.tt * dIbsdif_dVb;

         if (selfheat)
            gcjsT = model_.tt * dIbsdif_dT + dcjsbs_dT * T3
                                       + dT3_dT * cjsbs;
         else
            gcjsT = 0.0;


         arg = 1.0 - (Vbd > DioMax ? DioMax : Vbd) / PhiBSWG;

         if (selfheat)
            darg_dT = (1 - arg) / PhiBSWG * dPhiBSWG_dT;

         if (MJSWG == 0.5) {
            dT3_dVb = 1.0 / sqrt(arg);

            if (selfheat) ddT3_dVb_dT = -0.5 * dT3_dVb / arg * darg_dT;
         }
         else {
            dT3_dVb = exp(-MJSWG * log(arg));

            if (selfheat) ddT3_dVb_dT = -MJSWG * dT3_dVb / arg * darg_dT;
         }
         T3 = (1.0 - arg * dT3_dVb) * PhiBSWG / (1.0 - MJSWG);

         if (selfheat)
            dT3_dT = (1.0 - arg * dT3_dVb) * dPhiBSWG_dT / (1.0 - MJSWG)
            - (arg * ddT3_dVb_dT + darg_dT * dT3_dVb) * PhiBSWG
                                                / (1.0 - MJSWG);

         if (Vbd > DioMax)
            T3 += dT3_dVb * (Vbd - DioMax);

         dT3_dVd = -dT3_dVb;

         qjd_local = cjdbs * T3 + model_.tt * Ibddif;
         gcjdbs = cjdbs * dT3_dVb + model_.tt * dIbddif_dVb;
         gcjdds = cjdbs * dT3_dVd + model_.tt * dIbddif_dVd;

         if (selfheat)
            gcjdT = model_.tt * dIbddif_dT + dcjdbs_dT * T3
                                              + dT3_dT * cjdbs;
         else
            gcjdT = 0.0;
      }
      // v3.1

      qdrn -= qjd_local;
      qbody += (qjs_local + qjd_local);
      qsrc = -(qgate + qbody + qdrn + qsub);

      // Update the conductance
      cddb -= gcjdds;
      cdT -= gcjdT;
      cdsb += gcjdds + gcjdbs;

      cbdb += (gcjdds);
      cbT += (gcjdT + gcjsT);
      cbsb -= (gcjdds + gcjdbs + gcjsbs);

      // Extrinsic Bottom S/D to substrate charge
      T10 = -model_.dtype * ves;
      // T10 is vse without dtype conversion
      T11 = model_.dtype * (vds - ves);
      // T11 is vde without dtype conversion

      if (model_.csdmin != 0.0)
      {
         if ( ((paramPtr->nsub > 0) && (model_.dtype > 0)) ||
              ((paramPtr->nsub < 0) && (model_.dtype < 0)) )
         {
            if (T10 < paramPtr->vsdfb)
            {  qse = csbox * (T10 - paramPtr->vsdfb);
               gcse = csbox;
            }
            else if (T10 < paramPtr->sdt1)
            {  T0 = T10 - paramPtr->vsdfb;
               T1 = T0 * T0;
               qse = T0 * (csbox - paramPtr->st2 / 3 * T1) ;
               gcse = csbox - paramPtr->st2 * T1;
            }
            else if (T10 < paramPtr->vsdth)
            {  T0 = T10 - paramPtr->vsdth;
               T1 = T0 * T0;
               qse = csmin * T10 + st4 + paramPtr->st3 / 3 * T0 * T1;
               gcse = csmin + paramPtr->st3 * T1;
            }
            else
            {  qse = csmin * T10 + st4;
               gcse = csmin;
            }
         } else
         {
            if (T10 < paramPtr->vsdth)
            {  qse = csmin * (T10 - paramPtr->vsdth);
               gcse = csmin;
            }
            else if (T10 < paramPtr->sdt1)
            {  T0 = T10 - paramPtr->vsdth;
               T1 = T0 * T0;
               qse = T0 * (csmin - paramPtr->st2 / 3 * T1) ;
               gcse = csmin - paramPtr->st2 * T1;
            }
            else if (T10 < paramPtr->vsdfb)
            {  T0 = T10 - paramPtr->vsdfb;
               T1 = T0 * T0;
               qse = csbox * T10 + st4 + paramPtr->st3 / 3 * T0 * T1;
               gcse = csbox + paramPtr->st3 * T1;
            }
            else
            {  qse = csbox * T10 + st4;
               gcse = csbox;
            }
         }

         if ( ((paramPtr->nsub > 0) && (model_.dtype > 0)) ||
              ((paramPtr->nsub < 0) && (model_.dtype < 0)) )
         {
            if (T11 < paramPtr->vsdfb)
            {  qde = cdbox * (T11 - paramPtr->vsdfb);
               gcde = cdbox;
            }
            else if (T11 < paramPtr->sdt1)
            {  T0 = T11 - paramPtr->vsdfb;
               T1 = T0 * T0;
               qde = T0 * (cdbox - paramPtr->dt2 / 3 * T1) ;
               gcde = cdbox - paramPtr->dt2 * T1;
            }
            else if (T11 < paramPtr->vsdth)
            {  T0 = T11 - paramPtr->vsdth;
               T1 = T0 * T0;
               qde = cdmin * T11 + dt4 + paramPtr->dt3 / 3 * T0 * T1;
               gcde = cdmin + paramPtr->dt3 * T1;
            }
            else
            {  qde = cdmin * T11 + dt4;
               gcde = cdmin;
            }
         } else
         {
            if (T11 < paramPtr->vsdth)
            {  qde = cdmin * (T11 - paramPtr->vsdth);
               gcde = cdmin;
            }
            else if (T11 < paramPtr->sdt1)
            {  T0 = T11 - paramPtr->vsdth;
               T1 = T0 * T0;
               qde = T0 * (cdmin - paramPtr->dt2 / 3 * T1) ;
               gcde = cdmin - paramPtr->dt2 * T1;
            }
            else if (T11 < paramPtr->vsdfb)
            {  T0 = T11 - paramPtr->vsdfb;
               T1 = T0 * T0;
               qde = cdbox * T11 + dt4 + paramPtr->dt3 / 3 * T0 * T1;
               gcde = cdbox + paramPtr->dt3 * T1;
            }
            else
            {  qde = cdbox * T11 + dt4;
               gcde = cdbox;
            }
         }
      }
      else {
         qse = csbox * T10;
         gcse = csbox;
         qde = cdbox * T11;
         gcde = cdbox;
      }

      // Extrinsic : Sidewall fringing S/D charge
      qse += csesw * T10;
      gcse += csesw;
      qde += cdesw * T11;
      gcde += cdesw;

      // All charge are mutliplied with dtype at the end, but qse and qde
      //   have true polarity => so pre-mutliplied with dtype
      qse *= model_.dtype;
      qde *= model_.dtype;
  }

  cbb = Cbb;
  cbd = Cbd;
  cbg = Cbg;
  qbf = -Qsub0 - Qac0;
  qjs = qjs_local;
  qjd = qjd_local;

  // Setting up a few currents for the RHS load:
  // Note:  Vdp and Vsp get limited, but Vd and Vs do not.
  if (dNodePrime == 1)  // if drain resistor exists
  {
    Idrain   = drainConductance  * (Vd - Vdp);
    Idrain_Jdxp = 0.0;

    if (devOptions.voltageLimiterFlag && !origFlag)
    {                       // no limit on Vd, just Vdp
      tmp = - drainConductance  * ( - (Vdp-Vdp_orig));
      Idrain_Jdxp  = tmp;
    }
  }

  if (sNodePrime == 1)  // if source resistor exists
  {
    Isource  = sourceConductance * (Vs - Vsp);
    Isource_Jdxp = 0.0;

    if (devOptions.voltageLimiterFlag && !origFlag)
    {                       // no limit on Vs, just Vsp
      tmp = - sourceConductance  * ( - (Vsp-Vsp_orig));
      Isource_Jdxp  = tmp;
    }
  }

  if (selfheat)
    Ith = delTemp/paramPtr->rth;
  else
    Ith = 0;

  Igate = IgateMid = 0.0;
  if(rgateMod == 1)
  {
    Igate = grgeltd * (Vg - Vgp);
    Igate_Jdxp = 0.0;

    if (devOptions.voltageLimiterFlag && !origFlag)
    {
      tmp = - grgeltd * ((Vg-Vg_orig) - (Vgp-Vgp_orig));
      Igate_Jdxp  = tmp;
     }
  }
  else if(rgateMod == 2)
  {
    // gcrg is already the combined conductance of Rgeltd and Rii
    // No need to calculate jdxp here, because in fact ceqgcrg's jdxp
    // already handles it.  This is not the case for rgateMod==1
    Igate = (gcrg) * (Vg - Vgp);
    Igate_Jdxp = 0.0;
  }
  else if(rgateMod == 3)
  {
    // Need to check this --- might need some Jdxp terms for Igate, but
    // not for IgateMid, whose Jdxp terms are already in ceqgcrg
    Igate = grgeltd * (Vg - Vgm);
    IgateMid = gcrg * (Vgm - Vgp);
    Igate_Jdxp = 0.0;
  }

  // There is a spice3f5 convergence check that would happen here.
  // (line 2404) skipping...

  // In 3f5, loading a bunch of things into the state vector at this point.
  // (line 2433) skipping...

  // The following lines of charge computation are necessary in this routine,
  // and can't be deferred to auxChargeCalculations:  TVR ADDITION

  setupCapacitors_newDAE ();

  if (mode >= 0)
  {   Gm = gm;
      Gmbs = gmbs;
      Gme = gme;

      GmT = model_.dtype * gmT;

      FwdSum = Gm + Gmbs + Gme;   // v3.0
      RevSum = 0.0;

// spice has:
//      cdreq = model_.dtype * (cdrain - gds * vds
//            - Gm * vgs - Gmbs * vbs - Gme * ves) - GmT * delTemp;

      cdreq_Jdxp = 0.0;
      cdreq = model_.dtype * cdrain;
      if (devOptions.voltageLimiterFlag && !origFlag)
      {
        tmp = model_.dtype *
               (-gds *  (vds-vds_orig)
                -Gm  *  (vgs-vgs_orig)
                -Gmbs * (vbs-vbs_orig)
                -Gme  * (ves-ves_orig))
                -GmT * (delTemp-delTemp_orig);
        cdreq_Jdxp  = tmp;
      }

      // ceqbs now is compatible with cdreq, ie. going in is +ve
      // Equivalent current source from the diode
      ceqbs = cjs;
      ceqbd = cjd;

      ceqbs_Jdxp = cjs_Jdxp;
      ceqbd_Jdxp = cjd_Jdxp;

      // Current going in is +ve
      ceqbody = -cbody;
      ceqbody_Jdxp = -cbody_Jdxp;

      ceqgate = cgate;
      ceqgate_Jdxp = cgate_Jdxp;

      gigg_jac = gigg;
      gigb_jac = gigb;
      gige_jac = gige; // v3.0
      gigs_jac = gigs;
      gigd_jac = gigd;
      gigT_jac = model_.dtype * gigT;

      ceqth = cth;
      ceqth_Jdxp = cth_Jdxp;
      ceqbodcon = cbodcon;
      ceqbodcon_Jdxp = cbodcon_Jdxp;

      gbbg  = -gbgs;
      gbbdp = -gbds;
      gbbb  = -gbbs;
      gbbp  = -gbps;
      gbbT  = -model_.dtype * gbT; // v3.0
      gbbe  = -gbes;
      gbbsp = - ( gbbg + gbbdp + gbbb + gbbp + gbbe);

      gddpg  = -gjdg;
      gddpdp = -gjdd;
      gddpb  = -gjdb;
      gddpT  = -model_.dtype * gjdT; // v3.0
      gddpe  = -gjde;
      gddpsp = - ( gddpg + gddpdp + gddpb + gddpe);

      gsspg  = -gjsg;
      gsspdp = -gjsd;
      gsspb  = -gjsb;
      gsspT  = -model_.dtype * gjsT; // v3.0
      gsspe  = 0.0;
      gsspsp = - (gsspg + gsspdp + gsspb + gsspe);

      gppb = -gbpbs;
      gppp = -gbpps;

      gTtg  = gtempg;
      gTtb  = gtempb;
      gTtdp = gtempd;
      gTtt  = gtempT;

      // v3.0
      gTte  = gtempe;
      gTtsp = - (gTtg + gTtb + gTtdp + gTte);


      // v3.0
      if (model_.igcMod)
      {   gIstotg = gIgsg + gIgcsg;
          gIstotd = gIgcsd;
          gIstots = gIgss + gIgcss;
          gIstotb = gIgcsb;
//          Istoteq = model_.dtype * (Igs + Igcs
//                  - gIstotg * vgs - gIgcsd * vds - gIgcsb * vbs);

          Istoteq_Jdxp = 0.0;
          Istoteq = model_.dtype * (Igs + Igcs);
          if (devOptions.voltageLimiterFlag && !origFlag)
          {
            // NOTE: mimicing spice, we don't use the "Vgs/Vds/Vbs" versions
            tmp = model_.dtype *
                 (- gIstotg * (vgs-vgs_orig)
                  - gIgcsd  * (vds-vds_orig)
                  - gIgcsb  * (vbs-vbs_orig));
            Istoteq_Jdxp  = tmp;
          }

          gIdtotg = gIgdg + gIgcdg;
          gIdtotd = gIgdd + gIgcdd;
          gIdtots = gIgcds;
          gIdtotb = gIgcdb;
//          Idtoteq = model_.dtype * (Igd + Igcd
//                 - gIgdg * vgd - gIgcdg * vgs - gIgcdd * vds - gIgcdb * vbs);

          Idtoteq_Jdxp = 0.0;
          Idtoteq = model_.dtype * (Igd + Igcd);
          if (devOptions.voltageLimiterFlag && !origFlag)
          {
            // NOTE: mimicing SPICE, we use the original vgd/vgs/vds not
            // the "mode aware" "Vgs/Vds" etc.
            tmp = model_.dtype *
                  (- gIgdg  * (vgd-vgd_orig)
                   - gIgcdg * (vgs-vgs_orig)
                   - gIgcdd * (vds-vds_orig)
                   - gIgcdb * (vbs-vbs_orig));
            Idtoteq_Jdxp  = tmp;
          }

          gIgtotg = gIstotg + gIdtotg;
          gIgtotd = gIstotd + gIdtotd;
          gIgtots = gIstots + gIdtots;
          gIgtotb = gIstotb + gIdtotb;
          Igtoteq = Istoteq + Idtoteq;
          Igtoteq_Jdxp = 0;
          if (devOptions.voltageLimiterFlag && !origFlag)
          {
            Igtoteq_Jdxp = Istoteq_Jdxp + Idtoteq_Jdxp;
          }
      }
      else
      {   gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
          gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;

          gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;
      }

      // v3.1 wanh added for RF
      if (rgateMod == 2)
         T0 = vges - vgs;
      else if (rgateMod == 3)
         T0 = vgms - vgs;
      if (rgateMod > 1)
      {
         gcrgd_jac = gcrgd * T0;
         gcrgg_jac = gcrgg * T0;
         gcrgs_jac = gcrgs * T0;
         gcrgb_jac = gcrgb * T0;

//         ceqgcrg = -(gcrgd_jac * vds + gcrgg_jac * vgs + gcrgb_jac * vbs);

         ceqgcrg_Jdxp = 0.0;
         ceqgcrg = 0.0;
         if (devOptions.voltageLimiterFlag && !origFlag)
         {
           tmp = gcrgd_jac * (vds - vds_orig) +
                 gcrgg_jac * (vgs - vgs_orig) +
                 gcrgb_jac * (vbs - vbs_orig);
           ceqgcrg_Jdxp -= tmp;
         }

         gcrgg_jac -= gcrg;
         gcrg_jac = gcrg;
      }
      else
         ceqgcrg = gcrg_jac = gcrgd_jac = gcrgg_jac
                 = gcrgs_jac = gcrgb_jac = 0.0;
      // v3.1 wanh added for RF end

  } // end of mode>=0

  else
  {   Gm = -gm;
      Gmbs = -gmbs;
      Gme = -gme; // v3.0

      GmT = -model_.dtype * gmT;
      FwdSum = 0.0;
      RevSum = -(Gm + Gmbs + Gme); // v3.0


      // v3.1 bug fix
//      cdreq = -model_.dtype * (cdrain + gds*vds + Gm * vgd
//                    + Gmbs * vbd + Gme * (ves - vds)) - GmT * delTemp;

      cdreq_Jdxp = 0.0;
      cdreq = -model_.dtype * cdrain;
      if (devOptions.voltageLimiterFlag && !origFlag)
      {
        tmp = -model_.dtype * (gds   * (vds-vds_orig)
                                     +Gm   * (vgd-vgd_orig)
                                     +Gmbs * (vbd-vbd_orig)
                                     +Gme  * (ves-ves_orig-(vds-vds_orig))
                                     -GmT  * (delTemp-delTemp_orig));
        cdreq_Jdxp  = tmp;
      }

      ceqbs = cjd;
      ceqbd = cjs;

      ceqbs_Jdxp = cjd_Jdxp;
      ceqbd_Jdxp = cjs_Jdxp;

      // Current going in is +ve
      ceqbody = -cbody;
      ceqbody_Jdxp = -cbody_Jdxp;

      ceqgate = cgate;
      ceqgate_Jdxp = cgate_Jdxp;

      gigg_jac = gigg;
      gigb_jac = gigb;
      gige_jac = gige; // v3.0
      gigs_jac = gigd;
      gigd_jac = gigs;
      gigT_jac = model_.dtype * gigT;

      ceqth = cth;
      ceqth_Jdxp = cth_Jdxp;
      ceqbodcon = cbodcon;
      ceqbodcon_Jdxp = cbodcon_Jdxp;

      gbbg  = -gbgs;
      gbbb  = -gbbs;
      gbbp  = -gbps;
      gbbsp = -gbds;
      gbbT  = -model_.dtype * gbT;
      gbbe  = -gbes; // v3.0
      gbbdp = - ( gbbg + gbbsp + gbbb + gbbp + gbbe);

      gddpg  = -gjsg;
      gddpsp = -gjsd;
      gddpb  = -gjsb;
      gddpT  = -model_.dtype * gjsT;
      gddpe  = 0.0; // v3.0
      gddpdp = - (gddpg + gddpsp + gddpb + gddpe);

      gsspg  = -gjdg;
      gsspsp = -gjdd;
      gsspb  = -gjdb;
      gsspT  = -model_.dtype * gjdT; // v3.0
      gsspe  = -gjde;
      gsspdp = - ( gsspg + gsspsp + gsspb + gsspe);

      gppb = -gbpbs;
      gppp = -gbpps;

      gTtg  = gtempg;
      gTtb  = gtempb;
      gTtsp = gtempd;
      gTtt  = gtempT;

      // v3.0
      gTte  = gtempe;
      gTtdp = - (gTtg + gTtb + gTtsp + gTte);

      // v3.0
      if (model_.igcMod)
      {   gIstotg = gIgsg + gIgcdg;
          gIstotd = gIgcds;
          gIstots = gIgss + gIgcdd;
          gIstotb = gIgcdb;
//          Istoteq = model_.dtype * (Igs + Igcd
//                 - gIgsg * vgs - gIgcdg * vgd + gIgcdd * vds - gIgcdb * vbd);

          Istoteq_Jdxp = 0.0;
          Istoteq = model_.dtype * (Igs + Igcd);
          if (devOptions.voltageLimiterFlag && !origFlag)
          {
            tmp = model_.dtype *
                 (- gIgsg  * (vgs-vgs_orig)
                  - gIgcdg * (vgd-vgd_orig)
                  + gIgcdd * (vds-vds_orig)
                  - gIgcdb * (vbd-vbd_orig));
            Istoteq_Jdxp  = tmp;
          }


          gIdtotg = gIgdg + gIgcsg;
          gIdtotd = gIgdd + gIgcss;
          gIdtots = gIgcsd;
          gIdtotb = gIgcsb;
//          Idtoteq = model_.dtype * (Igd + Igcs
//                  - (gIgdg + gIgcsg) * vgd + gIgcsd * vds - gIgcsb * vbd);

          Idtoteq_Jdxp = 0.0;
          Idtoteq = model_.dtype * (Igd + Igcs);
          if (devOptions.voltageLimiterFlag && !origFlag)
          {
            tmp = model_.dtype *
                  (- (gIgdg + gIgcsg) * (vgd-vgd_orig)
                   + gIgcsd * (vds-vds_orig)
                   - gIgcsb * (vbd-vbd_orig));
            Idtoteq_Jdxp  = tmp;
          }

          gIgtotg = gIstotg + gIdtotg;
          gIgtotd = gIstotd + gIdtotd;
          gIgtots = gIstots + gIdtots;
          gIgtotb = gIstotb + gIdtotb;
          Igtoteq = Istoteq + Idtoteq;

          Igtoteq_Jdxp = 0;
          if (devOptions.voltageLimiterFlag && !origFlag)
          {
            Igtoteq_Jdxp = Istoteq_Jdxp + Idtoteq_Jdxp;
          }
      }
      else
      {   gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
          gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;

          gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;
      }

      // v3.1 wanh added for RF
      if (rgateMod == 2)
         T0 = vges - vgs;
      else if (rgateMod == 3)
         T0 = vgms - vgs;
      if (rgateMod > 1)
      {
         gcrgd_jac = gcrgs * T0;
         gcrgg_jac = gcrgg * T0;
         gcrgs_jac = gcrgd * T0;
         gcrgb_jac = gcrgb * T0;

         // ceqgcrg = -(gcrgg_jac * vgd - gcrgs_jac * vds + gcrgb_jac * vbd);
         ceqgcrg_Jdxp = 0.0;
         ceqgcrg = 0.0;
         if (devOptions.voltageLimiterFlag && !origFlag)
         {
           tmp = gcrgg_jac * (vgd - vgd_orig) -
                 gcrgs_jac * (vds - vds_orig) +
                 gcrgb_jac * (vbd - vbd_orig);

           ceqgcrg_Jdxp -= tmp;
       }

         gcrgg_jac -= gcrg;
         gcrg_jac = gcrg;
      }
      else
      {
         ceqgcrg = gcrg_jac = gcrgd_jac = gcrgg_jac
                 = gcrgs_jac = gcrgb_jac = 0.0;
      }
      // v3.1 wanh added for RF end

  } // end of soimod<0

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::B3SOIlimit
// Purpose       : Limits the per-iteration change of absolute voltage values
//
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
double Instance::B3SOIlimit(double vnew, double vold,
                                                     double limit, int *check)
{
  double T0, T1;
  T0 = vnew - vold;
  T1 = fabs(T0);
  if (T1 > limit)
  {
    if (T0 > 0.0)
    {
      vnew = vold + limit;
    }
    else
    {
      vnew = vold - limit;
    }
    *check = 1;
  }
  return vnew;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupCapacitors_oldDAE
// Purpose       :
//
// Special Notes : This was taken from updateIntermediateVars.
//                 This function calculates capacitance *conductances*.
//                 (C/dt terms), and also the charge variables:
//                 qb, qg, qe, qd, qth, qgmid
//
//                 Another way of looking at this function - this is analogous
//                 to all the capacitor-related code that gets executed
//                 in the SPICE3 SOI device, prior to the NIintegrate
//                 functions being called.
//
//                 In the old-dae form, this function needs to be called
//                 before charges can be put into the state vector
//                 (in updatePrimaryState).
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/25/05
//-----------------------------------------------------------------------------
bool Instance::setupCapacitors_oldDAE ()
{
  bool bsuccess = true;

  // these are temporary local variables.
  double T0, T1, T2, T3, T4;
  double cgdo_local;
  double cgso_local;
  double ag0 = solState.pdt;

  // ERK. 12/17/2006.
  // It is necessary to set ag0=0.0, because for the first time step out of
  // the DCOP, all the time derivatives are forced to be zero.  Thus, all
  // their derivatives should also be zero.  If it wasn't for that, then ag0
  // could always be pdt.  (it used to be, before the -jacobian_test capability).
  if (!(solState.dcopFlag) && solState.initTranFlag && solState.newtonIter==0)
  {
    ag0 = 0.0;
  }

  if (ChargeComputationNeeded)
  {

    T0 = vgd + DELTA_1;
    if (rgateMod == 3) T0 = vgmd + DELTA_1; /* v3.2 bug fix */
    T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
    T2 = 0.5 * (T0 - T1);

    /* v2.2.3 bug fix */
    T3 = wdiodCV_NoSwap * paramPtr->cgdl; /* v3.1 bug fix */

    T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappa);
    cgdo_local = paramPtr->cgdo + T3 - T3 * (1.0 - 1.0 / T4)
      * (0.5 - 0.5 * T0 / T1);
    qgdo = (paramPtr->cgdo + T3) * vgd - T3 * (T2
            + 0.5 * paramPtr->ckappa * (T4 - 1.0));

    if (rgateMod == 3) {
      qgdo = (paramPtr->cgdo + T3) * vgmd - T3 * (T2
               + 0.5 * paramPtr->ckappa * (T4 - 1.0));
    } /* v3.2 bug fix */

    T0 = vgs + DELTA_1;
    if (rgateMod == 3) T0 = vgms + DELTA_1; /* v3.2 bug fix */
    T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
    T2 = 0.5 * (T0 - T1);

    /* v2.2.3 bug fix */
    T3 = wdiosCV_NoSwap * paramPtr->cgsl; /* v3.1 bug fix */

    T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappa);
    cgso_local = paramPtr->cgso + T3 - T3 * (1.0 - 1.0 / T4)
      * (0.5 - 0.5 * T0 / T1);
    qgso = (paramPtr->cgso + T3) * vgs - T3 * (T2
      + 0.5 * paramPtr->ckappa * (T4 - 1.0));

    if (rgateMod == 3) {
      qgso = (paramPtr->cgso + T3) * vgms - T3 * (T2
      + 0.5 * paramPtr->ckappa * (T4 - 1.0));
    } /* v3.2 bug fix */

    if (mode > 0)
    {

      /* v3.1 wanh added for RF */
      if (rgateMod == 3)
      {
        gcgmgmb = (cgdo_local + cgso_local + paramPtr->cgeo) * ag0;

        gcgmdb = -cgdo_local * ag0;
        gcgmsb = -cgso_local * ag0;
        gcgmeb = -paramPtr->cgeo * ag0;
        gcdgmb = gcgmdb;
        gcsgmb = gcgmsb;
        gcegmb = gcgmeb;

        gcggb = cggb * ag0;
        gcgdb = cgdb * ag0;
        gcgsb = cgsb * ag0;
        gcgeb = 0 ;/*v3.1 wanh changed*/
        gcgbb = -(gcggb + gcgdb + gcgsb + gcgeb);

        gcdgb = cdgb * ag0;
        gcegb = gcgeb; /*v3.1 wanh added*/
        gcsgb = -(cggb + cbgb + cdgb) * ag0 - gcegb;
        gcbgb = cbgb * ag0;

        qgd = qgdo;
        qgs = qgso;
        qge = 0; /* v3.1 wanh changed */

        qgme = paramPtr->cgeo * vgme;
        qgmid = qgdo + qgso + qgme;
        qgate += qge;
        qbody -= 0;
        qdrn += qde - qgd;
        qsub -= qgme + qse + qde;
        qsrc = -(qgate + qgmid + qbody + qdrn + qsub);
      }
      else
      {
        gcggb = (cggb + cgdo_local + cgso_local + paramPtr->cgeo) * ag0;
        gcgdb = (cgdb - cgdo_local) * ag0;
        gcgsb = (cgsb - cgso_local) * ag0;
        gcgeb = (-paramPtr->cgeo) *ag0;
        gcgbb = -(gcggb + gcgdb + gcgsb + gcgeb);

        gcegb = (- paramPtr->cgeo) * ag0;
        gcdgb = (cdgb - cgdo_local) * ag0;
        gcsgb = -(cggb + cbgb
                  + cdgb + cgso_local) * ag0;
        gcbgb = cbgb * ag0;

        gcdgmb = gcsgmb = gcegmb = 0.0;
        gcgmdb = gcgmsb = gcgmeb = 0.0;

        /* Lump the overlap capacitance and S/D parasitics */
        qgd = qgdo;
        qgs = qgso;
        qge = paramPtr->cgeo * vge;
        qgate += qgd + qgs + qge;
        qdrn += qde - qgd;
        qsub -= qge + qse + qde;
        qsrc = -(qgate + qbody + qdrn + qsub);
      }

      gcddb = (cddb + cgdo_local + gcde) * ag0;
      gcdsb = cdsb * ag0;
      gcdeb = (cdeb - gcde) * ag0;
      gcdT = model_.dtype * cdT * ag0;

      gcsdb = -(cgdb + cbdb + cddb) * ag0;
      gcssb = (cgso_local + gcse - (cgsb + cbsb + cdsb)) * ag0;
      gcseb = -(gcse + cbeb + cdeb + ceeb) * ag0;
      gcsT = - model_.dtype * (cgT + cbT + cdT + ceT) * ag0;

      gcgT = model_.dtype * cgT * ag0;

      gcbdb = cbdb * ag0;
      gcbsb = cbsb * ag0;
      gcbeb = cbeb * ag0;
      gcbT = model_.dtype * cbT * ag0;

      gcedb = (- gcde) * ag0;
      gcesb = (- gcse) * ag0;
      gceeb = (gcse + gcde +
               ceeb + paramPtr->cgeo) * ag0;

      gceT = model_.dtype * ceT * ag0;

      gcTt = paramPtr->cth * ag0;

      sxpart = 0.6;
      dxpart = 0.4;


      /* v3.1 wanh moved the following original code ahead */
      /* Lump the overlap capacitance and S/D parasitics */
      /*                  qgd = qgdo;
                          qgs = qgso;
                          qge = paramPtr->cgeo * vge;
                          qgate += qgd + qgs + qge;
                          qdrn += qde - qgd;
                          qsub -= qge + qse + qde;
                          qsrc = -(qgate + qbody + qdrn + qsub);
      */
      /* v3.1 wanh end */

    }

    else
    {
      if (rgateMod == 3)
      {
        gcgmgmb = (cgdo_local + cgso_local + paramPtr->cgeo) * ag0;
        gcgmdb = -cgdo_local * ag0;
        gcgmsb = -cgso_local * ag0;
        gcgmeb = -paramPtr->cgeo * ag0;
        gcdgmb = gcgmdb;
        gcsgmb = gcgmsb;
        gcegmb = gcgmeb;

        gcggb = cggb * ag0;
        gcgsb = cgdb * ag0;
        gcgdb = cgsb * ag0;
        gcgeb = 0; /* v3.1 wanh changed */
        gcgbb = -(gcggb + gcgdb + gcgsb + gcgeb); /* v3.1 wanh added gcgeb */

        gcsgb = cdgb * ag0;
        gcegb = gcgeb; /* v3.1 wanh added */
        gcdgb = -(cggb + cbgb
                  + cdgb) * ag0 - gcegb; /*v3.1 wanh added gcegb*/
        gcbgb = cbgb * ag0;

        qgd = qgdo;
        qgs = qgso;
        qge = 0; /* v3.1 wanh changed */
        qgme = paramPtr->cgeo * vgme;
        qgmid = qgdo + qgso + qgme;
        qgate += qge;
        qbody -= 0;
        qsrc = qdrn - qgs + qse;
        qsub -= qgme + qse + qde;
        qdrn = -(qgate + qgmid + qbody + qsrc + qsub);
      }
      else
      {
        gcggb = (cggb + cgdo_local + cgso_local + paramPtr->cgeo) * ag0;
        gcgdb = (cgsb - cgdo_local) * ag0;
        gcgsb = (cgdb - cgso_local) * ag0;
        gcgeb = (- paramPtr->cgeo) * ag0;
        gcgbb = -(gcggb + gcgdb + gcgsb + gcgeb); /*wanh added gcgbb*/

        gcegb = gcgeb; /* v3.1 wanh added */
        gcsgb = (cdgb - cgso_local) * ag0;
        gcdgb = -(cggb + cbgb + cdgb + cgdo_local) * ag0;
        gcdgmb = gcsgmb = gcegmb = 0.0;
        gcgmdb = gcgmsb = gcgmeb = 0.0;

        /* Lump the overlap capacitance and S/D parasitics */
        qgd = qgdo;
        qgs = qgso;
        qge = paramPtr->cgeo * vge;
        qgate += qgd + qgs + qge;
        qsrc = qdrn - qgs + qse;
        qsub -= qge + qse + qde;
        qdrn = -(qgate + qbody + qsrc + qsub);
      }

      gcssb = (cddb + cgso_local + gcse) * ag0;
      gcsdb = cdsb * ag0;
      gcseb = (cdeb - gcse) * ag0;
      gcsT = model_.dtype * cdT * ag0;
      gcdsb = -(cgdb + cbdb + cddb) * ag0;
      gcddb = (cgdo_local + gcde - (cgsb + cbsb
                                    + cdsb)) * ag0;
      gcdeb = -(gcde + cbeb + cdeb
                + ceeb) * ag0;
      gcdT = - model_.dtype * (cgT + cbT
                                     + cdT + ceT) * ag0;

      gcgT = model_.dtype * cgT * ag0;

      gcbgb = cbgb * ag0;
      gcbsb = cbdb * ag0;
      gcbdb = cbsb * ag0;
      gcbeb = cbeb * ag0;
      gcbT = model_.dtype * cbT * ag0;

      /*                  gcegb = (-paramPtr->cgeo) * ag0; V3.2 bug fix */
      gcesb = (- gcse) * ag0;
      gcedb = (- gcde) * ag0;
      gceeb = (ceeb + paramPtr->cgeo +
               gcse + gcde) * ag0;
      gceT = model_.dtype * ceT * ag0;

      gcTt = paramPtr->cth * ag0;

      dxpart = 0.6;
      sxpart = 0.4;


      /* v3.1 wanh moved the following code ahead */
      /* Lump the overlap capacitance */
      /*
        qgd = qgdo;
        gs = qgso;
        qge = paramPtr->cgeo * vge;
        qgate += qgd + qgs + qge;
        qsrc = qdrn - qgs + qse;
        qsub -= qge + qse + qde;
        qdrn = -(qgate + qbody + qsrc + qsub);
      */
      /* v3.1 wanh end*/


    }

    cgdo = cgdo_local;
    cgso = cgso_local;

    qe = qsub;
    qg = qgate;
    qd = qdrn;
    qb = qbody;
    if ((model_.shMod == 1) && (rth0!=0.0))
      qth = paramPtr->cth * delTemp;
    if (rgateMod == 3)
      qgmid = qgmid;
    /* 3.1 bug fix */
  }
  else // !chargeComputationNeeded
  {
    //     initialize some stuff to 0
    // Original spice comment:  v3.1 wanh added for RF
    gcgmgmb = gcgmdb = gcgmsb = gcgmeb = 0.0;
    gcdgmb = gcsgmb = gcegmb = 0.0;
    gcgbb = 0.0;
  }
  // END TVR ADDITION

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
//
// Purpose       : This function sets up the primaray state variables into
//                 the primary state vector.
//
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  double * staVec = extData.nextStaVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;

  bsuccess = updateIntermediateVars ();

  // voltage drops:
  stoVec[li_store_vbd] = vbd;
  stoVec[li_store_vbs] = vbs;
  stoVec[li_store_vgs] = vgs;
  stoVec[li_store_vds] = vds;
  stoVec[li_store_ves] = ves;
  stoVec[li_store_vps] = vps;

  stoVec[li_store_vgp] = Vgp;
  stoVec[li_store_vd] = Vdp;
  stoVec[li_store_vs] = Vsp;
  stoVec[li_store_vp] = Vp;
  stoVec[li_store_ve] = Ve;
  stoVec[li_store_vg] = Vg;
  stoVec[li_store_vgm] = Vgm;
  stoVec[li_store_deltemp] = delTemp;

  stoVec[li_store_vges] = vges;
  stoVec[li_store_vgms] = vgms;

  // intrinsic capacitors:
  staVec[li_state_qb] = qb;
  staVec[li_state_qg] = qg;
  staVec[li_state_qd] = qd;
  staVec[li_state_qe] = qe;
  staVec[li_state_qgmid] = qgmid;
  staVec[li_state_qth] = qth;

  // if this is the first newton step of the first time step
  // of the transient simulation, we need to enforce that the
  // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
  // compatibility.  ERK.

  // Note:  I think this kind of thing is enforced (or should be enforced,
  // anyway) at the time integration level.  So I'm not sure this step is
  // really needed, at least for new-DAE.  Derivatives out of the DCOP
  // are supposed to be zero at the first newton step.

  if (!(solState.dcopFlag) && solState.initTranFlag && solState.newtonIter==0)
  {
    // re-set the state vector pointer that we are using to the "current"
    // pointer, rather than the "next" pointer.
    double * currStaVec = extData.currStaVectorRawPtr;

    // intrinsic capacitors:
    currStaVec[li_state_qb] = qb;
    currStaVec[li_state_qg] = qg;
    currStaVec[li_state_qd] = qd;
    currStaVec[li_state_qe] = qe;
    currStaVec[li_state_qgmid] = qgmid;
    currStaVec[li_state_qth] = qth;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 voltage source instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) + B(t) = 0
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  auxChargeCalculations ();

  double * qVec = extData.daeQVectorRawPtr;
  double * dQdxdVpPtr = extData.dQdxdVpVectorRawPtr;

  double Coef_body=0.0;
  double Coef_extBody=0.0;
  double Coef_gate=0.0;
  double Coef_gatePrime=0.0;
  double Coef_gateMid=0.0;
  double Coef_drainPrime=0.0;
  double Coef_sourcePrime=0.0;
  double Coef_substrate=0.0;
  double Coef_temp=0.0;

  if (soiMod != 2)
  {
    Coef_body  -= (model_.dtype*(Qeqqb))* numberParallel;
  }

  Coef_gatePrime  -= (model_.dtype*(Qeqqg))* numberParallel;

  Coef_drainPrime  += (model_.dtype*(-Qeqqd)) * numberParallel;

  Coef_sourcePrime  += ( + model_.dtype*(+ Qeqqg + Qeqqb
                    + Qeqqd + Qeqqe + Qeqqgmid)) * numberParallel;

  Coef_substrate  -= model_.dtype*Qeqqe * numberParallel;

  if (rgateMod == 3)
  {
    Coef_gateMid -= model_.dtype*(Qeqqgmid) * numberParallel;
  }

  if (selfheat)
  {
    Coef_temp -= (Qeqqth) * numberParallel;
  }

  ////////////////////////////////////////////////////////////////////////
  // Now place all the calculated values in the Q vector.

  if(li_Body != -1)
  {
    qVec[li_Body       ] -= Coef_body;
  }
  if(li_ExtBody != -1)
  {
    qVec[li_ExtBody    ] -= Coef_extBody;
  }

  qVec[li_Gate       ] -= Coef_gate;
  qVec[li_GatePrime  ] -= Coef_gatePrime;
  qVec[li_GateMid    ] -= Coef_gateMid;
  qVec[li_DrainPrime ] -= Coef_drainPrime;
  qVec[li_SourcePrime] -= Coef_sourcePrime;
  qVec[li_Substrate  ] -= Coef_substrate;

  if(li_Temperature != -1)
  {
    qVec[li_Temperature] -= Coef_temp;
  }

  if( loadLeadCurrent )
  {
    double * storeLeadQ = extData.storeLeadCurrQCompRawPtr;
    storeLeadQ[li_store_dev_id] = 0.0;
    storeLeadQ[li_store_dev_ig] = -Coef_gate;
    storeLeadQ[li_store_dev_is] = 0.0;
    storeLeadQ[li_store_dev_ie] = -Coef_substrate;
    if(li_Body != -1)
    {
      storeLeadQ[li_store_dev_ib] = -Coef_body;
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // Voltage limiting section:
  double Coef_body_Jdxp        = 0.0;
  double Coef_extBody_Jdxp     = 0.0;
  double Coef_gate_Jdxp        = 0.0;
  double Coef_gatePrime_Jdxp   = 0.0;
  double Coef_gateMid_Jdxp     = 0.0;
  double Coef_drainPrime_Jdxp  = 0.0;
  double Coef_sourcePrime_Jdxp = 0.0;
  double Coef_substrate_Jdxp   = 0.0;
  double Coef_temp_Jdxp        = 0.0;

  if (devOptions.voltageLimiterFlag)
  {
    if (soiMod != 2)
    {
      Coef_body_Jdxp  -=
      (model_.dtype*(Qeqqb_Jdxp))* numberParallel;
    }

    Coef_gatePrime_Jdxp  -=
      (model_.dtype*(Qeqqg_Jdxp))* numberParallel;

    Coef_drainPrime_Jdxp  +=
      (model_.dtype*(- Qeqqd_Jdxp)) * numberParallel;

    Coef_sourcePrime_Jdxp  += (+model_.dtype*(Qeqqg_Jdxp + Qeqqb_Jdxp
                        + Qeqqd_Jdxp + Qeqqe_Jdxp + Qeqqgmid_Jdxp)
                         ) * numberParallel;

    Coef_substrate_Jdxp  -= model_.dtype*Qeqqe_Jdxp * numberParallel;

    if (rgateMod == 3)
    {
      Coef_gateMid_Jdxp -=
      model_.dtype*(Qeqqgmid_Jdxp) * numberParallel;
    }

    if (selfheat)
    {
      Coef_temp_Jdxp -= (Qeqqth_Jdxp ) * numberParallel;
    }

    if(li_Body != -1)
    {
      dQdxdVpPtr[li_Body       ] += Coef_body_Jdxp;
    }
    if(li_ExtBody != -1)
    {
      dQdxdVpPtr[li_ExtBody    ] += Coef_extBody_Jdxp;
    }
    dQdxdVpPtr[li_Gate       ] += Coef_gate_Jdxp;
    dQdxdVpPtr[li_GatePrime  ] += Coef_gatePrime_Jdxp;
    dQdxdVpPtr[li_GateMid    ] += Coef_gateMid_Jdxp;
    dQdxdVpPtr[li_DrainPrime ] += Coef_drainPrime_Jdxp;
    dQdxdVpPtr[li_SourcePrime] += Coef_sourcePrime_Jdxp;
    dQdxdVpPtr[li_Substrate  ] += Coef_substrate_Jdxp;

    if(li_Temperature != -1)
    {
      dQdxdVpPtr[li_Temperature] += Coef_temp_Jdxp;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 B3SOI instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * dFdxdVpPtr = extData.dFdxdVpVectorRawPtr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Instance::loadDAEFVector" << endl;
    cout << "  name = " << getName() <<endl;
    cout.width(28); cout.precision(20); cout.setf(ios::scientific);
    cout << " " << endl;
  }
#endif

  double Coef_body=0.0;
  double Coef_extBody=0.0;
  double Coef_gate=0.0;
  double Coef_gatePrime=0.0;
  double Coef_gateMid=0.0;
  double Coef_drain=0.0;
  double Coef_drainPrime=0.0;
  double Coef_source=0.0;
  double Coef_sourcePrime=0.0;
  double Coef_substrate=0.0;
  double Coef_temp=0.0;

  Gmin = devOptions.gmin * 1e-6;
  geltd = grgeltd;

  // I have elected to add in the Gmin-based Jdxp terms right here.
  // In adding them, I've made use of the fact that whether we need voltage
  // limiting or not, we have to add Gmin*somedrop.  If we need voltage
  // limiting, we have to add in Gmin*(somedrop-somedrop_orig) with the
  // opposite sign.  I've simplified that down to the single term that
  // remains, which is correct even if voltage limiting is not needed
  // --- TVR

  double iGmin_bs = model_.dtype*(Gmin*vbs_orig);
  double iGmin_gd = model_.dtype*(Gmin*vgd_orig);

  if (soiMod != 2)
  {
    Coef_body  -= (model_.dtype*(ceqbody)
                   + iGmin_bs)* numberParallel;
  }

  Coef_gatePrime  -= (model_.dtype*(ceqgate - ceqgcrg)
                      + Igtoteq + iGmin_gd)* numberParallel;

  Coef_drainPrime  += (model_.dtype*(ceqbd)  - cdreq + Idtoteq
                       + Idrain + iGmin_gd) * numberParallel;

  Coef_sourcePrime  += (cdreq
                        + model_.dtype*(ceqbs)
                        + Istoteq + Isource + iGmin_bs) * numberParallel;

  if (rgateMod == 2)
  {
    Coef_gate -= model_.dtype*ceqgcrg * numberParallel;
  }
  else if (rgateMod == 3)
  {
    Coef_gateMid -= model_.dtype*(ceqgcrg) * numberParallel;
  }

  if (bodyMod == 1)
  {
    Coef_extBody += model_.dtype*ceqbodcon * numberParallel;
  }

  if (selfheat)
  {
    Coef_temp -= (ceqth + Ith) * numberParallel;
  }

  ////////////////////////////////////////////////////////////////////////
  // This next section deals with linear resistor currents (mostly) that
  // would not be part of the spice3f5 RHS load.
  //
  // These include the source and drain load resistors, the body tie
  // resistor and the various options for gate resistors.
  //
  // Isource and Idrain correspond to the linear load resistors.

  if (sNodePrime == 1) Coef_source  -= Isource * numberParallel;
  if (dNodePrime == 1) Coef_drain -= Idrain * numberParallel;

  // Note: extra body nodes were already handled, above.

  // Handle extra gate nodes, if they exist.
  // rgateMod==0   no gate resistor.
  // rgateMod==1   linear gate resistor
  // rgateMod==2   nonlinear gate resistor
  // rgateMod==3   2 gate resistors, in series.
  //
  if (rgateMod > 0)
  {
    Coef_gate -= Igate * numberParallel;
    if (rgateMod == 3)
    {
      Coef_gateMid -= (IgateMid - Igate) * numberParallel;
      Coef_gatePrime += IgateMid * numberParallel;
    }
    else
    {
      Coef_gatePrime += Igate * numberParallel;
    }
  }

  // Now place all the calculated values in the RHS vector.
  if(li_Body != -1)
  {
    fVec[li_Body       ] -= Coef_body;
  }
  if(li_ExtBody != -1)
  {
    fVec[li_ExtBody    ] -= Coef_extBody;
  }
  fVec[li_Gate       ] -= Coef_gate;
  fVec[li_GatePrime  ] -= Coef_gatePrime;
  fVec[li_GateMid    ] -= Coef_gateMid;
  fVec[li_Drain      ] -= Coef_drain;
  fVec[li_DrainPrime ] -= Coef_drainPrime;
  fVec[li_Source     ] -= Coef_source;
  fVec[li_SourcePrime] -= Coef_sourcePrime;
  fVec[li_Substrate  ] -= Coef_substrate;

  if(li_Temperature != -1)
  {
    fVec[li_Temperature] -= Coef_temp;
  }

  if( loadLeadCurrent )
  {
    double * storeLeadF = extData.nextStoVectorRawPtr;
    storeLeadF[li_store_dev_id] = -Coef_drain;
    storeLeadF[li_store_dev_ig] = -Coef_gate;
    storeLeadF[li_store_dev_is] = -Coef_source;
    storeLeadF[li_store_dev_ie] = -Coef_substrate;
    if(li_Body != -1)
    {
      storeLeadF[li_store_dev_ib] = -Coef_body;
    }
    else
    {
      storeLeadF[li_store_dev_ib] = 0.0;
    }
  }

  if( solState.dcopFlag && icVDSGiven )
  {
    if ( icVDSGiven )
    {
      if ( icVDSGiven )
      {
        double coef = (extData.nextSolVectorRawPtr)[li_Ids];
        fVec[li_Drain] += coef;
        fVec[li_Source] += -coef;
        if( loadLeadCurrent )
        {
          double * storeLeadF = extData.nextStoVectorRawPtr;
          storeLeadF[li_store_dev_id] = coef;
          storeLeadF[li_store_dev_is] = coef;
        }
        double cVs = (extData.nextSolVectorRawPtr)[li_Source];
        double cVd = (extData.nextSolVectorRawPtr)[li_Drain];
        fVec[li_Ids] += (cVd - cVs - icVDS);
      }
      if ( icVGSGiven )
      {
        double coef = (extData.nextSolVectorRawPtr)[li_Igs];
        fVec[li_Gate] += coef;
        fVec[li_Source] += -coef;
        if( loadLeadCurrent )
        {
          double * storeLeadF = extData.nextStoVectorRawPtr;
          storeLeadF[li_store_dev_ig] = coef;
          storeLeadF[li_store_dev_is] = coef;
        }
        double cVs = (extData.nextSolVectorRawPtr)[li_Source];
        double cVg = (extData.nextSolVectorRawPtr)[li_Gate];
        fVec[li_Igs] += (cVg - cVs - icVGS);
      }
      if ( icVBSGiven )
      {
        double coef = (extData.nextSolVectorRawPtr)[li_Ibs];
        fVec[li_Body] += coef;
        fVec[li_Source] += -coef;
        if( loadLeadCurrent )
        {
          double * storeLeadF = extData.nextStoVectorRawPtr;
          storeLeadF[li_store_dev_ib] = coef;
          storeLeadF[li_store_dev_is] = coef;
        }
        double cVs = (extData.nextSolVectorRawPtr)[li_Source];
        double cVb = (extData.nextSolVectorRawPtr)[li_Body];
        fVec[li_Ibs] += (cVb - cVs - icVBS);
      }
      if ( icVESGiven )
      {
        double coef = (extData.nextSolVectorRawPtr)[li_Ies];
        fVec[li_Substrate] += coef;
        fVec[li_Source] += -coef;
        if( loadLeadCurrent )
        {
          double * storeLeadF = extData.nextStoVectorRawPtr;
          storeLeadF[li_store_dev_ie] = coef;
          storeLeadF[li_store_dev_is] = coef;
        }
        double cVs = (extData.nextSolVectorRawPtr)[li_Source];
        double cVe = (extData.nextSolVectorRawPtr)[li_Substrate];
        fVec[li_Ies] += (cVe - cVs - icVES);
      }
      if ( icVPSGiven )
      {
        double coef = (extData.nextSolVectorRawPtr)[li_Ips];
        fVec[li_ExtBody] += coef;
        fVec[li_Source] += -coef;
        if( loadLeadCurrent )
        {
          double * storeLeadF = extData.nextStoVectorRawPtr;
          // this may need a separate store var if the internal and external body
          // nodes are not just connected by a resistor
          storeLeadF[li_store_dev_ib] = coef;
          storeLeadF[li_store_dev_is] = coef;
        }
        double cVs = (extData.nextSolVectorRawPtr)[li_Source];
        double cVp = (extData.nextSolVectorRawPtr)[li_ExtBody];
        fVec[li_Ies] += (cVp - cVs - icVPS);
      }
    }
  }

  // Set up the Jdxp vector:
  double Coef_body_Jdxp        = 0.0;
  double Coef_extBody_Jdxp     = 0.0;
  double Coef_gate_Jdxp        = 0.0;
  double Coef_gatePrime_Jdxp   = 0.0;
  double Coef_gateMid_Jdxp     = 0.0;
  double Coef_drain_Jdxp       = 0.0;
  double Coef_drainPrime_Jdxp  = 0.0;
  double Coef_source_Jdxp      = 0.0;
  double Coef_sourcePrime_Jdxp = 0.0;
  double Coef_substrate_Jdxp   = 0.0;
  double Coef_temp_Jdxp        = 0.0;

  double iGmin_bs_Jdxp = model_.dtype*Gmin*(vbs-vbs_orig);
  double iGmin_gd_Jdxp = model_.dtype*Gmin*(vgd-vgd_orig);

  if (devOptions.voltageLimiterFlag)
  {
    if (soiMod != 2)
    {
      Coef_body_Jdxp  -=
      (model_.dtype*(ceqbody_Jdxp )
         - iGmin_bs_Jdxp)* numberParallel;
    }

    Coef_gatePrime_Jdxp  -=
        (model_.dtype*(ceqgate_Jdxp - ceqgcrg_Jdxp)
         + Igtoteq_Jdxp - iGmin_gd_Jdxp)* numberParallel;

    Coef_drainPrime_Jdxp  +=
        (model_.dtype*(ceqbd_Jdxp)
         - cdreq_Jdxp + Idtoteq_Jdxp
         + Idrain_Jdxp - iGmin_gd_Jdxp) * numberParallel;

    Coef_sourcePrime_Jdxp  += (cdreq_Jdxp
      + model_.dtype*(ceqbs_Jdxp)
                               + Istoteq_Jdxp
                               + Isource_Jdxp - iGmin_bs_Jdxp
                          ) * numberParallel;

    Coef_substrate_Jdxp  -= model_.dtype*ceqqe_Jdxp * numberParallel;

    if (rgateMod == 2)
    {
      Coef_gate_Jdxp -=
        model_.dtype*ceqgcrg_Jdxp * numberParallel;
    }
    else if (rgateMod == 3)
    {
      Coef_gateMid_Jdxp -=
      model_.dtype*(ceqgcrg_Jdxp) * numberParallel;
    }

    if (bodyMod == 1)
    {
      Coef_extBody_Jdxp += model_.dtype*ceqbodcon_Jdxp * numberParallel;
    }

    if (selfheat)
    {
      Coef_temp_Jdxp -= ( ceqth_Jdxp ) * numberParallel;
    }

    if (sNodePrime == 1) Coef_source_Jdxp  -= Isource_Jdxp * numberParallel;
    if (dNodePrime == 1) Coef_drain_Jdxp -= Idrain_Jdxp * numberParallel;

    if (rgateMod > 0)
    {
      Coef_gate_Jdxp -= Igate_Jdxp * numberParallel;
      if (rgateMod == 3)
      {
        Coef_gateMid_Jdxp -= (IgateMid_Jdxp - Igate_Jdxp) * numberParallel;
        Coef_gatePrime_Jdxp += IgateMid_Jdxp * numberParallel;
      }
      else
      {
        Coef_gatePrime_Jdxp += Igate_Jdxp * numberParallel;
      }
    }

    if(li_Body != -1)
    {
      dFdxdVpPtr[li_Body       ] += Coef_body_Jdxp;
    }
    if(li_ExtBody != -1)
    {
      dFdxdVpPtr[li_ExtBody    ] += Coef_extBody_Jdxp;
    }
    dFdxdVpPtr[li_Gate       ] += Coef_gate_Jdxp;
    dFdxdVpPtr[li_GatePrime  ] += Coef_gatePrime_Jdxp;
    dFdxdVpPtr[li_GateMid    ] += Coef_gateMid_Jdxp;
    dFdxdVpPtr[li_Drain      ] += Coef_drain_Jdxp;
    dFdxdVpPtr[li_DrainPrime ] += Coef_drainPrime_Jdxp;
    dFdxdVpPtr[li_Source     ] += Coef_source_Jdxp;
    dFdxdVpPtr[li_SourcePrime] += Coef_sourcePrime_Jdxp;
    dFdxdVpPtr[li_Substrate  ] += Coef_substrate_Jdxp;

    if(li_Temperature != -1)
    {
      dFdxdVpPtr[li_Temperature] += Coef_temp_Jdxp;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupCapacitors_newDAE
// Purpose       :
//
// Special Notes : This function is identical to setupCapacitors_oldDAE,
//                 except that in that function, capacitance
//                 *conductances*. (C/dt terms) are set up, whereas in
//                 this function, just capacitances (C).
//
//                 As a result, a lot of gc* variables have been replaced
//                 with a lot of CAPc* variables.  Also, ag0 isn't used
//                 in this function, whereas it is used in the old-DAE
//                 version.
//
//                 This function also sets up the charge variables:
//                 qb, qg, qe, qd, qth, qgmid
//
//                 Another way of looking at this function - this is analogous
//                 to all the capacitor-related code that gets executed
//                 in the SPICE3 SOI device, prior to the NIintegrate
//                 functions being called.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/23/05
//-----------------------------------------------------------------------------
bool Instance::setupCapacitors_newDAE ()
{
  // these are temporary local variables.
  double T0, T1, T2, T3, T4;
  double cgdo_local;
  double cgso_local;

  if (ChargeComputationNeeded)
  {
    T0 = vgd + DELTA_1;
    if (rgateMod == 3) T0 = vgmd + DELTA_1; /* v3.2 bug fix */
    T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
    T2 = 0.5 * (T0 - T1);

    /* v2.2.3 bug fix */
    T3 = wdiodCV_NoSwap * paramPtr->cgdl; /* v3.1 bug fix */

    T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappa);
    cgdo_local = paramPtr->cgdo + T3 - T3 * (1.0 - 1.0 / T4)
      * (0.5 - 0.5 * T0 / T1);
    qgdo = (paramPtr->cgdo + T3) * vgd - T3 * (T2
            + 0.5 * paramPtr->ckappa * (T4 - 1.0));

    if (rgateMod == 3) {
      qgdo = (paramPtr->cgdo + T3) * vgmd - T3 * (T2
               + 0.5 * paramPtr->ckappa * (T4 - 1.0));
    }   /* v3.2 bug fix */

    T0 = vgs + DELTA_1;
    if (rgateMod == 3) T0 = vgms + DELTA_1; /* v3.2 bug fix */
    T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
    T2 = 0.5 * (T0 - T1);

    /* v2.2.3 bug fix */
    T3 = wdiosCV_NoSwap * paramPtr->cgsl; /* v3.1 bug fix */

    T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappa);
    cgso_local = paramPtr->cgso + T3 - T3 * (1.0 - 1.0 / T4)
      * (0.5 - 0.5 * T0 / T1);
    qgso = (paramPtr->cgso + T3) * vgs - T3 * (T2
      + 0.5 * paramPtr->ckappa * (T4 - 1.0));

    if (rgateMod == 3) {
      qgso = (paramPtr->cgso + T3) * vgms - T3 * (T2
      + 0.5 * paramPtr->ckappa * (T4 - 1.0));
    }   /* v3.2 bug fix */

    if (mode > 0)
    {

      /* v3.1 wanh added for RF */
      if (rgateMod == 3)
      {
        CAPcgmgmb = (cgdo_local + cgso_local + paramPtr->cgeo);

        CAPcgmdb = -cgdo_local;
        CAPcgmsb = -cgso_local;
        CAPcgmeb = -paramPtr->cgeo;
        CAPcdgmb = CAPcgmdb;
        CAPcsgmb = CAPcgmsb;
        CAPcegmb = CAPcgmeb;

        CAPcggb = cggb;
        CAPcgdb = cgdb;
        CAPcgsb = cgsb;
        CAPcgeb = 0 ;/*v3.1 wanh changed*/
        CAPcgbb = -(CAPcggb + CAPcgdb + CAPcgsb + CAPcgeb);

        CAPcdgb = cdgb;
        CAPcegb = CAPcgeb; /*v3.1 wanh added*/
        CAPcsgb = -(cggb + cbgb + cdgb) - CAPcegb;
        CAPcbgb = cbgb;

        qgd = qgdo;
        qgs = qgso;
        qge = 0; /* v3.1 wanh changed */

        qgme = paramPtr->cgeo * vgme;
        qgmid = qgdo + qgso + qgme;
        qgate += qge;
        qbody -= 0;
        qdrn += qde - qgd;
        qsub -= qgme + qse + qde;
        qsrc = -(qgate + qgmid + qbody + qdrn + qsub);
      }
      else
      {
        CAPcggb = (cggb + cgdo_local + cgso_local + paramPtr->cgeo);
        CAPcgdb = (cgdb - cgdo_local);
        CAPcgsb = (cgsb - cgso_local);
        CAPcgeb = (-paramPtr->cgeo);
        CAPcgbb = -(CAPcggb + CAPcgdb + CAPcgsb + CAPcgeb);

        CAPcegb = (- paramPtr->cgeo);
        CAPcdgb = (cdgb - cgdo_local);
        CAPcsgb = -(cggb + cbgb
                  + cdgb + cgso_local);
        CAPcbgb = cbgb;

        CAPcdgmb = CAPcsgmb = CAPcegmb = 0.0;
        CAPcgmdb = CAPcgmsb = CAPcgmeb = 0.0;

        /* Lump the overlap capacitance and S/D parasitics */
        qgd = qgdo;
        qgs = qgso;
        qge = paramPtr->cgeo * vge;
        qgate += qgd + qgs + qge;
        qdrn += qde - qgd;
        qsub -= qge + qse + qde;
        qsrc = -(qgate + qbody + qdrn + qsub);
      }

      CAPcddb = (cddb + cgdo_local + gcde);
      CAPcdsb = cdsb;
      CAPcdeb = (cdeb - gcde);
      CAPcdT = model_.dtype * cdT;

      CAPcsdb = -(cgdb + cbdb + cddb);
      CAPcssb = (cgso_local + gcse - (cgsb + cbsb + cdsb));
      CAPcseb = -(gcse + cbeb + cdeb + ceeb);
      CAPcsT = - model_.dtype * (cgT + cbT + cdT + ceT);

      CAPcgT = model_.dtype * cgT;

      CAPcbdb = cbdb;
      CAPcbsb = cbsb;
      CAPcbeb = cbeb;
      CAPcbT = model_.dtype * cbT;

      CAPcedb = (- gcde);
      CAPcesb = (- gcse);
      CAPceeb = (gcse + gcde +
               ceeb + paramPtr->cgeo);

      CAPceT = model_.dtype * ceT;

      CAPcTt = paramPtr->cth;

      sxpart = 0.6;
      dxpart = 0.4;


      /* v3.1 wanh moved the following original code ahead */
      /* Lump the overlap capacitance and S/D parasitics */
      /*                  qgd = qgdo;
                          qgs = qgso;
                          qge = paramPtr->cgeo * vge;
                          qgate += qgd + qgs + qge;
                          qdrn += qde - qgd;
                          qsub -= qge + qse + qde;
                          qsrc = -(qgate + qbody + qdrn + qsub);
      */
      /* v3.1 wanh end */

    }

    else
    {
      if (rgateMod == 3)
      {
        CAPcgmgmb = (cgdo_local + cgso_local + paramPtr->cgeo);
        CAPcgmdb = -cgdo_local;
        CAPcgmsb = -cgso_local;
        CAPcgmeb = -paramPtr->cgeo;
        CAPcdgmb = CAPcgmdb;
        CAPcsgmb = CAPcgmsb;
        CAPcegmb = CAPcgmeb;

        CAPcggb = cggb;
        CAPcgsb = cgdb;
        CAPcgdb = cgsb;
        CAPcgeb = 0; /* v3.1 wanh changed */
        CAPcgbb = -(CAPcggb + CAPcgdb + CAPcgsb + CAPcgeb); /* v3.1 wanh added gcgeb */

        CAPcsgb = cdgb;
        CAPcegb = CAPcgeb; /* v3.1 wanh added */
        CAPcdgb = -(cggb + cbgb
                  + cdgb) - CAPcegb; /*v3.1 wanh added gcegb*/
        CAPcbgb = cbgb;

        qgd = qgdo;
        qgs = qgso;
        qge = 0; /* v3.1 wanh changed */
        qgme = paramPtr->cgeo * vgme;
        qgmid = qgdo + qgso + qgme;
        qgate += qge;
        qbody -= 0;
        qsrc = qdrn - qgs + qse;
        qsub -= qgme + qse + qde;
        qdrn = -(qgate + qgmid + qbody + qsrc + qsub);
      }
      else
      {
        CAPcggb = (cggb + cgdo_local + cgso_local + paramPtr->cgeo);
        CAPcgdb = (cgsb - cgdo_local);
        CAPcgsb = (cgdb - cgso_local);
        CAPcgeb = (- paramPtr->cgeo);
        CAPcgbb = -(CAPcggb + CAPcgdb + CAPcgsb + CAPcgeb); /*wanh added gcgbb*/

        CAPcegb = CAPcgeb; /* v3.1 wanh added */
        CAPcsgb = (cdgb - cgso_local);
        CAPcdgb = -(cggb + cbgb + cdgb + cgdo_local);
        CAPcdgmb = CAPcsgmb = CAPcegmb = 0.0;
        CAPcgmdb = CAPcgmsb = CAPcgmeb = 0.0;

        /* Lump the overlap capacitance and S/D parasitics */
        qgd = qgdo;
        qgs = qgso;
        qge = paramPtr->cgeo * vge;
        qgate += qgd + qgs + qge;
        qsrc = qdrn - qgs + qse;
        qsub -= qge + qse + qde;
        qdrn = -(qgate + qbody + qsrc + qsub);
      }

      CAPcssb = (cddb + cgso_local + gcse);
      CAPcsdb = cdsb;
      CAPcseb = (cdeb - gcse);
      CAPcsT = model_.dtype * cdT;
      CAPcdsb = -(cgdb + cbdb + cddb);
      CAPcddb = (cgdo_local + gcde - (cgsb + cbsb
                                    + cdsb));
      CAPcdeb = -(gcde + cbeb + cdeb + ceeb);
      CAPcdT = - model_.dtype * (cgT + cbT
                                     + cdT + ceT);

      CAPcgT = model_.dtype * cgT;

      CAPcbgb = cbgb;
      CAPcbsb = cbdb;
      CAPcbdb = cbsb;
      CAPcbeb = cbeb;
      CAPcbT = model_.dtype * cbT;

      /*                  CAPcegb = (-paramPtr->cgeo); V3.2 bug fix */
      CAPcesb = (- gcse);
      CAPcedb = (- gcde);
      CAPceeb = (ceeb + paramPtr->cgeo + gcse + gcde);
      CAPceT = model_.dtype * ceT;

      CAPcTt = paramPtr->cth;

      dxpart = 0.6;
      sxpart = 0.4;


      /* v3.1 wanh moved the following code ahead */
      /* Lump the overlap capacitance */
      /*
        qgd = qgdo;
        gs = qgso;
        qge = paramPtr->cgeo * vge;
        qgate += qgd + qgs + qge;
        qsrc = qdrn - qgs + qse;
        qsub -= qge + qse + qde;
        qdrn = -(qgate + qbody + qsrc + qsub);
      */
      /* v3.1 wanh end*/


    }

    cgdo = cgdo_local;
    cgso = cgso_local;

    qe = qsub;
    qg = qgate;
    qd = qdrn;
    qb = qbody;
    if ((model_.shMod == 1) && (rth0!=0.0))
      qth = paramPtr->cth * delTemp;
    if (rgateMod == 3)
      qgmid = qgmid;
    /* 3.1 bug fix */
  }
  else
  {
    //     initialize some stuff to 0
    // Original spice comment:  v3.1 wanh added for RF
    CAPcgmgmb = CAPcgmdb = CAPcgmsb = CAPcgmeb = 0.0;
    CAPcdgmb = CAPcsgmb = CAPcegmb = 0.0;
    CAPcgbb = 0.0;
  }
  // END TVR ADDITION

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::auxChargeCalculations
// Purpose       :
//
// Special Notes :
//                 Another way of looking at this function - this is analogous
//                 to the capacitor-related code that gets executed in the SPICE3
//                 SOI device, after the NIintegrate functions have been called.
//
//                 About all this function really does is set up some
//                 voltlim terms.
//                 It really is only here to provide capacitor-related
//                 voltage limiter terms for the QVector load.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/23/05
//-----------------------------------------------------------------------------
bool Instance::auxChargeCalculations ()
{
  bool bsuccess = true;
  double vgb_orig, vgmb_orig, veb_orig;

  if (!ChargeComputationNeeded)
  {
    sxpart = (1.0 - (dxpart = (mode > 0) ? 0.4 : 0.6));
  }
  else  // ChargeComputation is needed
  {
    Qeqqg = qg;
    Qeqqb = qb;
    Qeqqd = qd;
    Qeqqe = qe;
    Qeqqth = qth;
    Qeqqgmid = qgmid;

    // Note - in the SPICE version of the SOI, the variables B3SOIcd and
    // B3SOIcb contain capacitor current contributions.  So, one might
    // think that equivalent charge variables should be set up here
    // (Qcd, Qcb, etc).  This isn't necessary, it turns out.  cd and
    // cb are only used in SPICE's convergence and bypass analysis, and
    // are not ever loaded (directly or indirectly) into the SPICE
    // rhs vector.

    // evaluate equivalent charge current
    vgb_orig = vgs_orig - vbs_orig;
    vgmb_orig = vgms_orig - vbs_orig;
    veb_orig = ves_orig - vbs_orig;

//    ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs
//            - gcgeb * veb - gcgT * delTemp;

    Qeqqg_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqg_Jdxp  = - CAPcggb * (vgb-vgb_orig)
                    + CAPcgdb * (vbd-vbd_orig)
                    + CAPcgsb * (vbs-vbs_orig)
                    + CAPcgeb * (veb-veb_orig)
                    + CAPcgT * (delTemp-delTemp_orig);
    }

//    ceqqb = cqbody - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs
//            - gcbeb * veb - gcbT * delTemp;

    Qeqqb_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqb_Jdxp  = - CAPcbgb * (vgb-vgb_orig)
                    + CAPcbdb * (vbd-vbd_orig)
                    + CAPcbsb * (vbs-vbs_orig)
                    - CAPcbeb * (veb-veb_orig)
                    -  CAPcbT * (delTemp-delTemp_orig);
    }

//    ceqqd = cqdrn - gcdgb * vgb + gcddb * vbd + gcdsb * vbs
//            - gcdeb * veb - gcdT * delTemp -gcdgmb * vgmb;

    Qeqqd_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqd_Jdxp  = - CAPcdgb * (vgb-vgb_orig)
                    + CAPcddb * (vbd-vbd_orig)
                    + CAPcdsb * (vbs-vbs_orig)
                    - CAPcdeb * (veb-veb_orig)
                    -  CAPcdT * (delTemp-delTemp_orig)
                    -CAPcdgmb * (vgmb-vgmb_orig);
    }

//    ceqqe = cqsub - gcegb * vgb + gcedb * vbd + gcesb * vbs
//            - gceeb * veb - gceT * delTemp - gcegmb * vgmb;

    Qeqqe_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqe_Jdxp  = - CAPcegb * (vgb-vgb_orig)
                    + CAPcedb * (vbd-vbd_orig)
                    + CAPcesb * (vbs-vbs_orig)
                    - CAPceeb * (veb-veb_orig)
                    -  CAPceT * (delTemp-delTemp_orig)
                    -CAPcegmb * (vgmb-vgmb_orig);
    }

//    ceqqth = cqtemp - gcTt * delTemp;

    Qeqqth_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqth_Jdxp  = -CAPcTt * (delTemp-delTemp_orig);
    }

    if (rgateMod == 3)
    {
//    ceqqgmid = *(ckt->CKTstate0 + cqgmid)
//        + gcgmdb * vbd + gcgmsb * vbs - gcgmgmb * vgmb;

      Qeqqgmid_Jdxp = 0.0;
      if (!origFlag)
      {
        Qeqqgmid_Jdxp  =   CAPcgmdb  * (vbd-vbd_orig)
                         + CAPcgmsb  * (vbs-vbs_orig)
                         - CAPcgmgmb * (vgmb-vgmb_orig);
      }
    }
    else
    {
      Qeqqgmid_Jdxp = 0.0;
    }
  } // !ChargeComputationNeeded

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 B3SOI instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  if (rgateMod == 3)
  {
    dQdx[li_GateMid][AGateMidEquGateMidNodeOffset]
      += (CAPcgmgmb)*numberParallel;
    dQdx[li_GateMid][AGateMidEquDrainPrimeNodeOffset]
      += (CAPcgmdb)*numberParallel;

    dQdx[li_GateMid][AGateMidEquSourcePrimeNodeOffset]
      += (CAPcgmsb)*numberParallel;
    dQdx[li_GateMid][AGateMidEquSubstrateNodeOffset]
      += CAPcgmeb*numberParallel;

    dQdx[li_DrainPrime][ADrainPrimeEquGateMidNodeOffset]
      += CAPcdgmb*numberParallel;

    dQdx[li_SourcePrime][ASourcePrimeEquGateMidNodeOffset]
      += CAPcsgmb*numberParallel;
    dQdx[li_Substrate][ASubstrateEquGateMidNodeOffset]
      += CAPcegmb*numberParallel;
  }

  dQdx[li_Substrate][ASubstrateEquDrainPrimeNodeOffset]
      += CAPcedb*numberParallel;

  dQdx[li_DrainPrime][ADrainPrimeEquSubstrateNodeOffset]
      += CAPcdeb*numberParallel;
  dQdx[li_SourcePrime][ASourcePrimeEquSubstrateNodeOffset]
      += CAPcseb*numberParallel;
  dQdx[li_Substrate][ASubstrateEquGatePrimeNodeOffset]
      += CAPcegb*numberParallel;
  dQdx[li_GatePrime][AGatePrimeEquSubstrateNodeOffset]
      += CAPcgeb*numberParallel;
  if (soiMod != 2)
  {
    dQdx[li_Substrate][ASubstrateEquBodyNodeOffset]
      -= (CAPcegb + CAPcedb + CAPcesb + CAPceeb + CAPcegmb)*numberParallel;
    if (rgateMod == 0 || rgateMod == 1)
    {
      dQdx[li_GatePrime][AGatePrimeEquBodyNodeOffset]
      -= (+CAPcggb + CAPcgdb + CAPcgsb + CAPcgeb)*numberParallel;
    }
    else
    {
      dQdx[li_GatePrime][AGatePrimeEquBodyNodeOffset]
      += (+ CAPcgbb)*numberParallel;
    }
    dQdx[li_DrainPrime][ADrainPrimeEquBodyNodeOffset]
      -= ((+ CAPcdgb + CAPcddb + CAPcdeb + CAPcdsb) + CAPcdgmb)*numberParallel;
    dQdx[li_SourcePrime][ASourcePrimeEquBodyNodeOffset]
      -= ((+ CAPcsgb + CAPcsdb + CAPcseb + CAPcssb) + CAPcsgmb)*numberParallel;
    dQdx[li_Body][ABodyEquSubstrateNodeOffset]
      += (+ CAPcbeb)*numberParallel;
    dQdx[li_Body][ABodyEquGatePrimeNodeOffset]
      += (+ CAPcbgb)*numberParallel;
    dQdx[li_Body][ABodyEquDrainPrimeNodeOffset]
      += (+ CAPcbdb)*numberParallel;
    dQdx[li_Body][ABodyEquSourcePrimeNodeOffset]
      += (CAPcbsb)*numberParallel;
    dQdx[li_Body][ABodyEquBodyNodeOffset]
      += (-CAPcbgb - CAPcbdb - CAPcbsb - CAPcbeb)*numberParallel;
  }
  dQdx[li_Substrate][ASubstrateEquSubstrateNodeOffset]
    += CAPceeb*numberParallel;
  if (rgateMod == 0)
  {
    dQdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (CAPcggb)*numberParallel;
    dQdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (CAPcgdb)*numberParallel;
    dQdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (CAPcgsb)*numberParallel;
  }
  else if (rgateMod == 1)
  {
    dQdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (CAPcggb)*numberParallel;
    dQdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (CAPcgdb)*numberParallel;
    dQdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (CAPcgsb)*numberParallel;
  }
  else
  {
    dQdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (CAPcggb)*numberParallel;
    dQdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (CAPcgdb)*numberParallel;
    dQdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (CAPcgsb)*numberParallel;
  }
  dQdx[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset]
    += ((CAPcdgb))*numberParallel;
  dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]
    += ((CAPcddb) )*numberParallel;
  dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]
    -= ((-CAPcdsb) )*numberParallel;

  dQdx[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset]
    += (CAPcsgb)*numberParallel;
  dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]
    -= (-CAPcsdb )*numberParallel;
  dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
    += ((CAPcssb) )*numberParallel;

  if (selfheat)
  {
    dQdx[li_DrainPrime][ADrainPrimeEquTemperatureNodeOffset]
      += (CAPcdT)*numberParallel;
    dQdx[li_SourcePrime][ASourcePrimeEquTemperatureNodeOffset]
      += (CAPcsT)*numberParallel;
    dQdx[li_Substrate][ASubstrateEquTemperatureNodeOffset]
      += CAPceT*numberParallel;
    dQdx[li_GatePrime][AGatePrimeEquTemperatureNodeOffset]
      += (CAPcgT)*numberParallel;
    dQdx[li_Temperature][ATemperatureEquTemperatureNodeOffset]
      += (CAPcTt)*numberParallel;

    if (bNode > 0)
    {
      dQdx[li_Body][ABodyEquTemperatureNodeOffset]
      += (CAPcbT)*numberParallel;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes : The F-vector is an algebraic constaint.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if (rgateMod == 1)
  {
    dFdx[li_Gate][AGateEquGateNodeOffset]
      += geltd*numberParallel;
    dFdx[li_GatePrime][AGatePrimeEquGateNodeOffset]
      -= geltd*numberParallel;
    dFdx[li_Gate][AGateEquGatePrimeNodeOffset]
      -= geltd*numberParallel;
    // It seems that this should be here, but it shouldn't.
    // The geltd term is added into the G'-G' element later down in this
    // routine
    //
    //    dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
    //      += geltd*numberParallel;
  }
  else if (rgateMod == 2)
  {
    dFdx[li_Gate][AGateEquGateNodeOffset]
      += gcrg_jac*numberParallel;
    dFdx[li_Gate][AGateEquGatePrimeNodeOffset]
      += gcrgg_jac*numberParallel;
    dFdx[li_Gate][AGateEquDrainPrimeNodeOffset]
      += gcrgd_jac*numberParallel;
    dFdx[li_Gate][AGateEquSourcePrimeNodeOffset]
      += gcrgs_jac*numberParallel;
    dFdx[li_GatePrime][AGatePrimeEquGateNodeOffset]
      -= gcrg_jac*numberParallel;
    if (soiMod != 2)
      dFdx[li_Gate][AGateEquBodyNodeOffset]
      += gcrgb_jac*numberParallel;
  }
  else if (rgateMod == 3)
  {
    dFdx[li_Gate][AGateEquGateNodeOffset]
      += geltd*numberParallel;
    dFdx[li_Gate][AGateEquGateMidNodeOffset]
      += -geltd*numberParallel;
    dFdx[li_GateMid][AGateMidEquGateNodeOffset]
      += -geltd*numberParallel;
    dFdx[li_GateMid][AGateMidEquGateMidNodeOffset]
      += (geltd + gcrg_jac)*numberParallel;
    dFdx[li_GateMid][AGateMidEquDrainPrimeNodeOffset]
      += (gcrgd_jac)*numberParallel;
    dFdx[li_GateMid][AGateMidEquGatePrimeNodeOffset]
      += gcrgg_jac*numberParallel;
    dFdx[li_GateMid][AGateMidEquSourcePrimeNodeOffset]
      += (gcrgs_jac)*numberParallel;

    if (soiMod != 2)
    {
      dFdx[li_GateMid][AGateMidEquBodyNodeOffset]
      += gcrgb_jac*numberParallel;
    }

    dFdx[li_GatePrime][AGatePrimeEquGateMidNodeOffset]
      -= gcrg_jac*numberParallel;
  }
  if (soiMod != 0)
  {
    dFdx[li_DrainPrime][ADrainPrimeEquSubstrateNodeOffset]
      += (Gme + gddpe)*numberParallel;
    dFdx[li_SourcePrime][ASourcePrimeEquSubstrateNodeOffset]
      += (gsspe - Gme)*numberParallel;
    if (soiMod != 2)
    {
      dFdx[li_GatePrime][AGatePrimeEquSubstrateNodeOffset]
      += gige_jac*numberParallel;
      dFdx[li_Body][ABodyEquSubstrateNodeOffset]
      -= gige_jac*numberParallel;
    }
  }

  if (soiMod != 2)
  {
    if (rgateMod == 0 || rgateMod == 1)
    {
      dFdx[li_GatePrime][AGatePrimeEquBodyNodeOffset]
      -= (-gigb_jac - gIgtotb)*numberParallel;
    }
    else
    {
      dFdx[li_GatePrime][AGatePrimeEquBodyNodeOffset]
      += (gigb_jac +gIgtotb - gcrgb_jac)*numberParallel;
    }
    dFdx[li_DrainPrime][ADrainPrimeEquBodyNodeOffset]
      -= ((-gddpb - Gmbs) + gIdtotb)*numberParallel;
    dFdx[li_SourcePrime][ASourcePrimeEquBodyNodeOffset]
      -= ((-gsspb + Gmbs) + Gmin + gIstotb)*numberParallel;
    dFdx[li_Body][ABodyEquSubstrateNodeOffset]
      += (gbbe)*numberParallel;
    dFdx[li_Body][ABodyEquGatePrimeNodeOffset]
      += (-gigg + gbbg)*numberParallel;
    dFdx[li_Body][ABodyEquDrainPrimeNodeOffset]
      += (-gigd_jac + gbbdp)*numberParallel;
    dFdx[li_Body][ABodyEquSourcePrimeNodeOffset]
      += (gbbsp - Gmin - gigs)*numberParallel;
    dFdx[li_Body][ABodyEquBodyNodeOffset]
      += (-gigb_jac + gbbb + Gmin)*numberParallel;
  }
  if (rgateMod == 0)
  {
    dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (gigg + Gmin + gIgtotg)*numberParallel;
    dFdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (gigd_jac - Gmin + gIgtotd)*numberParallel;
    dFdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (gigs + gIgtots)*numberParallel;
  }
  else if (rgateMod == 1)
  {
    dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (gigg_jac + Gmin + gIgtotg + geltd)*numberParallel;
    dFdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (gigd - Gmin + gIgtotd)*numberParallel;
    dFdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (gigs_jac + gIgtots)*numberParallel;
  }
  else
  {
    dFdx[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (gigg + Gmin + gIgtotg - gcrgg_jac)*numberParallel;
    dFdx[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (gigd_jac - Gmin + gIgtotd - gcrgd_jac)*numberParallel;
    dFdx[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (gigs + gIgtots - gcrgs_jac)*numberParallel;
  }
  dFdx[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset]
    += ((Gm) + gddpg - Gmin - gIdtotg)*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]
    += ((drainConductance + gds + gddpdp + RevSum) + Gmin - gIdtotd)*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]
    -= ((-gddpsp + gds + FwdSum) + gIdtots)*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]
    += -drainConductance*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset]
    += (-Gm + gsspg - gIstotg)*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]
    -= (gds - gsspdp + RevSum + gIstotd)*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
    += ((sourceConductance + gds + gsspsp + FwdSum) + Gmin - gIstots)*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]
    -= sourceConductance*numberParallel;
  dFdx[li_Drain][ADrainEquDrainNodeOffset]
    += drainConductance*numberParallel;
  dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset]
    -= drainConductance*numberParallel;
  dFdx[li_Source][ASourceEquSourceNodeOffset]
    += sourceConductance*numberParallel;
  dFdx[li_Source][ASourceEquSourcePrimeNodeOffset]
    -= sourceConductance*numberParallel;
  if (bodyMod == 1)
  {
    dFdx[li_Body][ABodyEquExtBodyNodeOffset]
      -= gppp*numberParallel;
    dFdx[li_ExtBody][AExtBodyEquBodyNodeOffset]
      += gppb*numberParallel;
    dFdx[li_ExtBody][AExtBodyEquExtBodyNodeOffset]
      += gppp*numberParallel;
  }
  if (selfheat)
  {
    dFdx[li_DrainPrime][ADrainPrimeEquTemperatureNodeOffset]
      += (GmT + gddpT)*numberParallel;
    dFdx[li_SourcePrime][ASourcePrimeEquTemperatureNodeOffset]
      += (-GmT + gsspT)*numberParallel;

    dFdx[li_GatePrime][AGatePrimeEquTemperatureNodeOffset]
      += (gigT_jac)*numberParallel;
    dFdx[li_Temperature][ATemperatureEquTemperatureNodeOffset]
      += (gTtt  + 1/paramPtr->rth)*numberParallel;
    dFdx[li_Temperature][ATemperatureEquGatePrimeNodeOffset]
      += gTtg*numberParallel;
    dFdx[li_Temperature][ATemperatureEquDrainPrimeNodeOffset]
      += gTtdp*numberParallel;
    dFdx[li_Temperature][ATemperatureEquSourcePrimeNodeOffset]
      += gTtsp*numberParallel;
    if (soiMod != 0)
    {
      dFdx[li_Temperature][ATemperatureEquSubstrateNodeOffset]
        += gTte*numberParallel;
    }
    if (bNode > 0)
    {
      dFdx[li_Body][ABodyEquTemperatureNodeOffset]
        += (gbbT - gigT_jac)*numberParallel;
      dFdx[li_Temperature][ATemperatureEquBodyNodeOffset]
        += gTtb*numberParallel;
    }
  }
  if( icVDSGiven )
  {
    if( solState.dcopFlag  )
    {
      dFdx[li_Drain][ADrainEquIdsOffset] += 1.0;
      dFdx[li_Source][ASourceEquIdsOffset] -= 1.0;
      dFdx[li_Ids][icVDSEquVdOffset] += 1.0;
      dFdx[li_Ids][icVDSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Ids][icVDSEquIdsOffset] += 1.0;
    }
  }

  if( icVGSGiven )
  {
    if( solState.dcopFlag  )
    {
      dFdx[li_Gate][AGateEquIgsOffset] += 1.0;
      dFdx[li_Source][ASourceEquIgsOffset] -= 1.0;
      dFdx[li_Igs][icVGSEquVgOffset] += 1.0;
      dFdx[li_Igs][icVGSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Igs][icVGSEquIgsOffset] += 1.0;
    }
  }

  if( icVBSGiven )
  {
    if( solState.dcopFlag  )
    {
      dFdx[li_Body][ABodyEquIbsOffset] += 1.0;
      dFdx[li_Source][ASourceEquIbsOffset] -= 1.0;
      dFdx[li_Ibs][icVBSEquVbOffset] += 1.0;
      dFdx[li_Ibs][icVBSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Ibs][icVBSEquIbsOffset] += 1.0;
    }
  }

  if( icVESGiven )
  {
    if( solState.dcopFlag  )
    {
      dFdx[li_Substrate][ASubstrateEquIesOffset] += 1.0;
      dFdx[li_Source][ASourceEquIesOffset] -= 1.0;
      dFdx[li_Ies][icVESEquVeOffset] += 1.0;
      dFdx[li_Ies][icVESEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Ies][icVESEquIesOffset] += 1.0;
    }
  }

  if( icVPSGiven )
  {
    if( solState.dcopFlag  )
    {
      dFdx[li_ExtBody][AExtBodyEquIpsOffset] += 1.0;
      dFdx[li_Source][ASourceEquIpsOffset] -= 1.0;
      dFdx[li_Ips][icVPSEquVpOffset] += 1.0;
      dFdx[li_Ips][icVPSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Ips][icVPSEquIpsOffset] += 1.0;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance:loadMatrix
//
// Purpose       : The contents of this function were in the
//                 loadAnalyticJacobian function previously.  I moved them to
//                 this new function as a debug exercise, b/c I needed to
//                 load more than one matrix, as a debug tool.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, 9233
// Creation Date : 01/06/05
//-----------------------------------------------------------------------------
bool Instance::loadMatrix (N_LAS_Matrix & JMat)
{
  bool bsuccess = true;

  if (rgateMod == 1)
  {
    (JMat)[li_Gate][AGateEquGateNodeOffset]
      += geltd*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquGateNodeOffset]
      -= geltd*numberParallel;
    (JMat)[li_Gate][AGateEquGatePrimeNodeOffset]
      -= geltd*numberParallel;
    // It seems that this should be here, but it shouldn't.
    // The geltd term is added into the G'-G' element later down in this
    // routine
    //
    //    (JMat)[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
    //      += geltd*numberParallel;
  }
  else if (rgateMod == 2)
  {
    (JMat)[li_Gate][AGateEquGateNodeOffset]
      += gcrg_jac*numberParallel;
    (JMat)[li_Gate][AGateEquGatePrimeNodeOffset]
      += gcrgg_jac*numberParallel;
    (JMat)[li_Gate][AGateEquDrainPrimeNodeOffset]
      += gcrgd_jac*numberParallel;
    (JMat)[li_Gate][AGateEquSourcePrimeNodeOffset]
      += gcrgs_jac*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquGateNodeOffset]
      -= gcrg_jac*numberParallel;
    if (soiMod != 2)
      (JMat)[li_Gate][AGateEquBodyNodeOffset]
        += gcrgb_jac*numberParallel;
  }
  else if (rgateMod == 3)
  {
    (JMat)[li_Gate][AGateEquGateNodeOffset]
      += geltd*numberParallel;
    (JMat)[li_Gate][AGateEquGateMidNodeOffset]
      += -geltd*numberParallel;
    (JMat)[li_GateMid][AGateMidEquGateNodeOffset]
      += -geltd*numberParallel;
    (JMat)[li_GateMid][AGateMidEquGateMidNodeOffset]
      += (geltd + gcrg_jac + gcgmgmb)*numberParallel;
    (JMat)[li_GateMid][AGateMidEquDrainPrimeNodeOffset]
      += (gcrgd_jac + gcgmdb)*numberParallel;
    (JMat)[li_GateMid][AGateMidEquGatePrimeNodeOffset]
      += gcrgg_jac*numberParallel;
    (JMat)[li_GateMid][AGateMidEquSourcePrimeNodeOffset]
      += (gcrgs_jac + gcgmsb)*numberParallel;
    (JMat)[li_GateMid][AGateMidEquSubstrateNodeOffset]
      += gcgmeb*numberParallel;
    if (soiMod != 2)
    {
      (JMat)[li_GateMid][AGateMidEquBodyNodeOffset]
        += gcrgb_jac*numberParallel;
    }
    (JMat)[li_DrainPrime][ADrainPrimeEquGateMidNodeOffset]
      += gcdgmb*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquGateMidNodeOffset]
      -= gcrg_jac*numberParallel;
    (JMat)[li_SourcePrime][ASourcePrimeEquGateMidNodeOffset]
      += gcsgmb*numberParallel;
    (JMat)[li_Substrate][ASubstrateEquGateMidNodeOffset]
      += gcegmb*numberParallel;
  }
  if (soiMod != 0)
  {
    (JMat)[li_DrainPrime][ADrainPrimeEquSubstrateNodeOffset]
      += (Gme + gddpe)*numberParallel;
    (JMat)[li_SourcePrime][ASourcePrimeEquSubstrateNodeOffset]
      += (gsspe - Gme)*numberParallel;
    if (soiMod != 2)
    {
      (JMat)[li_GatePrime][AGatePrimeEquSubstrateNodeOffset]
      += gige_jac*numberParallel;
      (JMat)[li_Body][ABodyEquSubstrateNodeOffset]
      -= gige_jac*numberParallel;
    }
  }
  (JMat)[li_Substrate][ASubstrateEquDrainPrimeNodeOffset]
      += gcedb*numberParallel;
  (JMat)[li_Substrate][ASubstrateEquSourcePrimeNodeOffset]
      += gcesb*numberParallel;
  (JMat)[li_DrainPrime][ADrainPrimeEquSubstrateNodeOffset]
      += gcdeb*numberParallel;
  (JMat)[li_SourcePrime][ASourcePrimeEquSubstrateNodeOffset]
      += gcseb*numberParallel;
  (JMat)[li_Substrate][ASubstrateEquGatePrimeNodeOffset]
      += gcegb*numberParallel;
  (JMat)[li_GatePrime][AGatePrimeEquSubstrateNodeOffset]
      += gcgeb*numberParallel;
  if (soiMod != 2)
  {
    (JMat)[li_Substrate][ASubstrateEquBodyNodeOffset]
      -= (gcegb + gcedb + gcesb + gceeb + gcegmb)*numberParallel;
    if (rgateMod == 0 || rgateMod == 1)
    {
      (JMat)[li_GatePrime][AGatePrimeEquBodyNodeOffset]
        -= (-gigb_jac + gcggb + gcgdb + gcgsb + gcgeb - gIgtotb)*numberParallel;
    }
    else
    {
      (JMat)[li_GatePrime][AGatePrimeEquBodyNodeOffset]
        += (gigb_jac + gcgbb +gIgtotb - gcrgb_jac)*numberParallel;
    }
    (JMat)[li_DrainPrime][ADrainPrimeEquBodyNodeOffset]
      -= ((-gddpb - Gmbs + gcdgb + gcddb + gcdeb + gcdsb) + gcdgmb + gIdtotb)*numberParallel;
    (JMat)[li_SourcePrime][ASourcePrimeEquBodyNodeOffset]
      -= ((-gsspb + Gmbs + gcsgb + gcsdb + gcseb + gcssb) + gcsgmb + Gmin + gIstotb)*numberParallel;
    (JMat)[li_Body][ABodyEquSubstrateNodeOffset]
      += (gbbe + gcbeb)*numberParallel;
    (JMat)[li_Body][ABodyEquGatePrimeNodeOffset]
      += (-gigg_jac + gcbgb + gbbg)*numberParallel;
    (JMat)[li_Body][ABodyEquDrainPrimeNodeOffset]
      += (-gigd_jac + gcbdb + gbbdp)*numberParallel;
    (JMat)[li_Body][ABodyEquSourcePrimeNodeOffset]
      += (gcbsb + gbbsp - Gmin - gigs_jac)*numberParallel;
    (JMat)[li_Body][ABodyEquBodyNodeOffset]
      += (-gigb_jac + gbbb - gcbgb - gcbdb - gcbsb - gcbeb + Gmin)*numberParallel;
  }
  (JMat)[li_Substrate][ASubstrateEquSubstrateNodeOffset]
    += gceeb*numberParallel;
  if (rgateMod == 0)
  {
    (JMat)[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (gigg_jac + gcggb + Gmin + gIgtotg)*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (gigd_jac + gcgdb - Gmin + gIgtotd)*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (gcgsb + gigs_jac + gIgtots)*numberParallel;
  }
  else if (rgateMod == 1)
  {
    (JMat)[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (gigg_jac + gcggb + Gmin + gIgtotg + geltd)*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (gigd_jac + gcgdb - Gmin + gIgtotd)*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (gcgsb + gigs_jac + gIgtots)*numberParallel;
  }
  else
  {
    (JMat)[li_GatePrime][AGatePrimeEquGatePrimeNodeOffset]
      += (gigg_jac + gcggb + Gmin + gIgtotg - gcrgg_jac)*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquDrainPrimeNodeOffset]
      += (gigd_jac + gcgdb - Gmin + gIgtotd - gcrgd_jac)*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquSourcePrimeNodeOffset]
      += (gcgsb + gigs_jac + gIgtots - gcrgs_jac)*numberParallel;
  }
  (JMat)[li_DrainPrime][ADrainPrimeEquGatePrimeNodeOffset]
    += ((Gm + gcdgb) + gddpg - Gmin - gIdtotg)*numberParallel;
  (JMat)[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]
    += ((drainConductance + gds + gddpdp + RevSum + gcddb) + Gmin - gIdtotd)*numberParallel;
  (JMat)[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]
    -= ((-gddpsp + gds + FwdSum - gcdsb) + gIdtots)*numberParallel;
  (JMat)[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]
    += -drainConductance*numberParallel;
  (JMat)[li_SourcePrime][ASourcePrimeEquGatePrimeNodeOffset]
    += (gcsgb - Gm + gsspg - gIstotg)*numberParallel;
  (JMat)[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]
    -= (gds - gsspdp + RevSum - gcsdb + gIstotd)*numberParallel;
  (JMat)[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
    += ((sourceConductance + gds + gsspsp + FwdSum + gcssb) + Gmin - gIstots)*numberParallel;
  (JMat)[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]
    -= sourceConductance*numberParallel;
  (JMat)[li_Drain][ADrainEquDrainNodeOffset]
    += drainConductance*numberParallel;
  (JMat)[li_Drain][ADrainEquDrainPrimeNodeOffset]
    -= drainConductance*numberParallel;
  (JMat)[li_Source][ASourceEquSourceNodeOffset]
    += sourceConductance*numberParallel;
  (JMat)[li_Source][ASourceEquSourcePrimeNodeOffset]
    -= sourceConductance*numberParallel;
  if (bodyMod == 1)
  {
    (JMat)[li_Body][ABodyEquExtBodyNodeOffset]
      -= gppp*numberParallel;
    (JMat)[li_ExtBody][AExtBodyEquBodyNodeOffset]
      += gppb*numberParallel;
    (JMat)[li_ExtBody][AExtBodyEquExtBodyNodeOffset]
      += gppp*numberParallel;
  }
  if (selfheat)
  {
    (JMat)[li_DrainPrime][ADrainPrimeEquTemperatureNodeOffset]
      += (GmT + gddpT + gcdT)*numberParallel;
    (JMat)[li_SourcePrime][ASourcePrimeEquTemperatureNodeOffset]
      += (-GmT + gsspT + gcsT)*numberParallel;
    (JMat)[li_Substrate][ASubstrateEquTemperatureNodeOffset]
      += gceT*numberParallel;
    (JMat)[li_GatePrime][AGatePrimeEquTemperatureNodeOffset]
      += (gcgT + gigT_jac)*numberParallel;
    (JMat)[li_Temperature][ATemperatureEquTemperatureNodeOffset]
      += (gTtt  + 1/paramPtr->rth + gcTt)*numberParallel;
    (JMat)[li_Temperature][ATemperatureEquGatePrimeNodeOffset]
      += gTtg*numberParallel;
    (JMat)[li_Temperature][ATemperatureEquDrainPrimeNodeOffset]
      += gTtdp*numberParallel;
    (JMat)[li_Temperature][ATemperatureEquSourcePrimeNodeOffset]
      += gTtsp*numberParallel;
    if (soiMod != 0)
    {
      (JMat)[li_Temperature][ATemperatureEquSubstrateNodeOffset]
        += gTte*numberParallel;
    }
    if (bNode > 0)
    {
      (JMat)[li_Body][ABodyEquTemperatureNodeOffset]
        += (gbbT + gcbT - gigT_jac)*numberParallel;
      (JMat)[li_Temperature][ATemperatureEquBodyNodeOffset]
        += gTtb*numberParallel;
    }
  }

  if( icVDSGiven )
  {
    if( solState.dcopFlag  )
      {
      (JMat)[li_Drain][ADrainEquIdsOffset] += 1.0;
      (JMat)[li_Source][ASourceEquIdsOffset] -= 1.0;
      (JMat)[li_Ids][icVDSEquVdOffset] += 1.0;
      (JMat)[li_Ids][icVDSEquVsOffset] -= 1.0;
    }
    else
    {
      (JMat)[li_Ids][icVDSEquIdsOffset] += 1.0;
    }
  }

  if( icVGSGiven )
  {
    if( solState.dcopFlag  )
    {
      (JMat)[li_Gate][AGateEquIgsOffset] += 1.0;
      (JMat)[li_Source][ASourceEquIgsOffset] -= 1.0;
      (JMat)[li_Igs][icVGSEquVgOffset] += 1.0;
      (JMat)[li_Igs][icVGSEquVsOffset] -= 1.0;
    }
    else
    {
      (JMat)[li_Igs][icVGSEquIgsOffset] += 1.0;
    }
  }

  if( icVBSGiven )
  {
    if( solState.dcopFlag  )
    {
      (JMat)[li_Body][ABodyEquIbsOffset] += 1.0;
      (JMat)[li_Source][ASourceEquIbsOffset] -= 1.0;
      (JMat)[li_Ibs][icVBSEquVbOffset] += 1.0;
      (JMat)[li_Ibs][icVBSEquVsOffset] -= 1.0;
    }
    else
    {
      (JMat)[li_Ibs][icVBSEquIbsOffset] += 1.0;
    }
  }

  if( icVESGiven )
  {
    if( solState.dcopFlag  )
    {
      (JMat)[li_Substrate][ASubstrateEquIesOffset] += 1.0;
      (JMat)[li_Source][ASourceEquIesOffset] -= 1.0;
      (JMat)[li_Ies][icVESEquVeOffset] += 1.0;
      (JMat)[li_Ies][icVESEquVsOffset] -= 1.0;
    }
    else
    {
      (JMat)[li_Ies][icVESEquIesOffset] += 1.0;
    }
  }

  if( icVPSGiven )
  {
    if( solState.dcopFlag  )
    {
      (JMat)[li_ExtBody][AExtBodyEquIpsOffset] += 1.0;
      (JMat)[li_Source][ASourceEquIpsOffset] -= 1.0;
      (JMat)[li_Ips][icVPSEquVpOffset] += 1.0;
      (JMat)[li_Ips][icVPSEquVsOffset] -= 1.0;
    }
    else
    {
      (JMat)[li_Ips][icVPSEquIpsOffset] += 1.0;
    }
  }

  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance:checkModel
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley
// Creation Date : 08/04/04
//-----------------------------------------------------------------------------

bool Instance::checkModel ()
{
  bool bsuccess = true;
  string msg;
  ostringstream err("");

  if (paramPtr->nlx < -paramPtr->leff)
  {
    err << "Fatal: Nlx = " << paramPtr->nlx << " is less than -Leff.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (model_.tox <= 0.0)
  {
    err << "Fatal: Tox = " << model_.tox << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (model_.toxm <= 0.0)
  {
    err << "Fatal: Toxm = " << model_.toxm << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (model_.tox - model_.dtoxcv <= 0.0)
  {
    err << "Fatal: Tox - dtoxcv = " << model_.tox - model_.dtoxcv
        << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (model_.tbox <= 0.0)
  {
    err << "Fatal: Tbox = " << model_.tbox << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->npeak <= 0.0)
  {
    err << "Fatal: Nch = " << paramPtr->npeak << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->ngate < 0.0)
  {
    err << "Fatal: Ngate = " << paramPtr->ngate << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->ngate > 1.e25)
  {
    err << "Fatal: Ngate = %" << paramPtr->ngate << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->dvt1 < 0.0)
  {
    err << "Fatal: Dvt1 = " << paramPtr->dvt1 << " is negative.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->dvt1w < 0.0)
  {
    err << "Fatal: Dvt1w = " << paramPtr->dvt1w << " is negative.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->w0 == -paramPtr->weff)
  {
    msg = "Fatal: (W0 + Weff) = 0 cauing divided-by-zero.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->dsub < 0.0)
  {
    err << "Fatal: Dsub = " << paramPtr->dsub << " is negative.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->b1 == -paramPtr->weff)
  {
    msg = "Fatal: (B1 + Weff) = 0 causing divided-by-zero.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->u0temp <= 0.0)
  {
    err << "Fatal: u0 at current temperature = " << paramPtr->u0temp
        << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->delta < 0.0)
  {
    err << "Fatal: Delta = " << paramPtr->delta << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->vsattemp <= 0.0)
  {
    err << "Fatal: Vsat at current temperature = " << paramPtr->vsattemp
        << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->pclm <= 0.0)
  {
    err << "Fatal: Pclm = " << paramPtr->pclm << " is not positive.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if (paramPtr->drout < 0.0)
  {
    err << "Fatal: Drout = " << paramPtr->drout << " is negative.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }
  if ( model_.unitLengthGateSidewallJctCap > 0.0)
  {
    if (drainPerimeter < paramPtr->weff)
    {
      if (devOptions.verboseLevel > 0)
      {
        err << "Warning: Pd = " << drainPerimeter << " is less than W.\n";
        msg = err.str();
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      }
      drainPerimeter = paramPtr->weff;
    }

    if (sourcePerimeter < paramPtr->weff)
    {
      if (devOptions.verboseLevel > 0)
      {
        err << "Warning: Ps = " << sourcePerimeter << " is less than W.\n";
        msg = err.str();
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      }
      sourcePerimeter =paramPtr->weff;
    }
  }
  if (paramPtr->clc < 0.0)
  {
    err << "Fatal: Clc = " << paramPtr->clc << " is negative.";
    msg = err.str();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    bsuccess = false;
  }

  if (devOptions.verboseLevel > 0)
  {
    if (paramPtr->noff < 0.1)
    {
      err << "Warning: Noff = " << paramPtr->noff << " is too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->noff > 4.0)
    {
      err << "Warning: Noff = " << paramPtr->noff << " is too large.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->moin < 5.0)
    {
      err << "Warning: Moin = " << paramPtr->moin << " is too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->moin > 25.0)
    {
      err << "Warning: Moin = " << paramPtr->moin << " is too large.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (model_.moinFD < 5.0)
    {
      err << "Warning: MoinFD = " << model_.moinFD << " is too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (model_.capMod == 3) {
      if (paramPtr->acde < 0.1)
      {
        err << "Warning: Acde = " << paramPtr->acde << " is too small.";
        msg = err.str();
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      }
      if (paramPtr->acde > 1.6)
      {
        err << "Warning: Acde = " << paramPtr->acde << " is too large.";
        msg = err.str();
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      }
    }
  }
  if (model_.paramChk ==1)
  {
    if (paramPtr->leff <= 5.0e-8)
    {
      err << "Warning: Leff = " << paramPtr->leff << " may be too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->leffCV <= 5.0e-8)
    {
      err << "Warning: Leff for CV = " << paramPtr->leffCV << " may be too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->weff <= 1.0e-7)
    {
      err << "Warning: Weff = " << paramPtr->weff << " may be too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->weffCV <= 1.0e-7)
    {
      err << "Warning: Weff for CV = " << paramPtr->weffCV << " may be too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->nlx < 0.0)
    {
      err << "Warning: Nlx = " << paramPtr->nlx << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (model_.tox < 1.0e-9)
    {
      err << "Warning: Tox = " << model_.tox << " is less than 10A.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->npeak <= 1.0e15)
    {
      err << "Warning: Nch = " << paramPtr->npeak << " may be too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    else if (paramPtr->npeak >= 1.0e21)
    {
      err << "Warning: Nch = " << paramPtr->npeak << " may be too large.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (fabs(paramPtr->nsub) >= 1.0e21)
    {
      err << "Warning: Nsub = " << paramPtr->nsub << " may be too large.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if ((paramPtr->ngate > 0.0) && (paramPtr->ngate <= 1.e18))
    {
      err << "Warning: Ngate = " << paramPtr->ngate << " is less than 1.E18cm^-3.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (paramPtr->dvt0 < 0.0)
    {
      err << "Warning: Dvt0 = " << paramPtr->dvt0 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
    if (fabs(1.0e-6 / (paramPtr->w0 + paramPtr->weff)) > 10.0)
    {
      msg = "Warning: (W0 + Weff) may be too small.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->nfactor < 0.0)
    {
      err << "Warning: Nfactor = " << paramPtr->nfactor << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->cdsc < 0.0)
    {
      err << "Warning: Cdsc = " << paramPtr->cdsc << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->cdscd < 0.0)
    {
      err << "Warning: Cdscd = " << paramPtr->cdscd << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->eta0 < 0.0)
    {
      err << "Warning: Eta0 = " << paramPtr->eta0 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (fabs(1.0e-6 / (paramPtr->b1 + paramPtr->weff)) > 10.0)
    {
      msg = "Warning: (B1 + Weff) may be too small.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->a2 < 0.01)
    {
      err << "Warning: A2 = " << paramPtr->a2 << " is too small. Set to 0.01";
      msg = err.str();
      paramPtr->a2 = 0.01;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    else if (paramPtr->a2 > 1.0)
    {
      err << "Warning: A2 = " << paramPtr->a2 << " is larger than 1. A2 is set to 1 and A1 is set to 0.";
      msg = err.str();
      paramPtr->a2 = 1.0;
      paramPtr->a1 = 0.0;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->rdsw < 0.0)
    {
      err << "Warning: Rdsw = " << paramPtr->rdsw << " is negative. Set to zero.";
      msg = err.str();
      paramPtr->rdsw = 0.0;
      paramPtr->rds0 = 0.0;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    else if ((paramPtr->rds0 > 0.0) && (paramPtr->rds0 < 0.001))
    {
      err << "Warning: Rds at current temperature = " << paramPtr->rds0 << " is less than 0.001 ohm. Set to zero.";
      msg = err.str();
      paramPtr->rds0 = 0.0;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->vsattemp < 1.0e3)
    {
      err << "Warning: Vsat at current temperature = " << paramPtr->vsattemp <<
             " may be too small.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->pdibl1 < 0.0)
    {
      err << "Warning: Pdibl1 = " << paramPtr->pdibl1 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->pdibl2 < 0.0)
    {
      err << "Warning: Pdibl2 = " << paramPtr->pdibl2 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.cgdo < 0.0)
    {
      err << "Warning: cgdo = " << model_.cgdo << " is negative. Set to zero.";
      msg = err.str();
      model_.cgdo = 0.0;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.cgso < 0.0)
    {
      err << "Warning: cgso = " << model_.cgso << " is negative. Set to zero.";
      msg = err.str();
      model_.cgso = 0.0;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.cgeo < 0.0)
    {
      err << "Warning: cgeo = " << model_.cgeo << " is negative. Set to zero.";
      msg = err.str();
      model_.cgeo = 0.0;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.ntun < 0.0)
    {
      err << "Warning: Ntun = " << model_.ntun << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.ndiode < 0.0)
    {
      err << "Warning: Ndiode = " << model_.ndiode << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.isbjt < 0.0)
    {
      err << "Warning: Isbjt = " << model_.isbjt << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.isdif < 0.0)
    {
      err << "Warning: Isdif = " << model_.isdif << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.isrec < 0.0)
    {
      err << "Warning: Isrec = " << model_.isrec << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.istun < 0.0)
    {
      err << "Warning: Istun = " << model_.istun << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.tt < 0.0)
    {
      err << "Warning: Tt = " << model_.tt << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.csdmin < 0.0)
    {
      err << "Warning: Csdmin = " << model_.csdmin << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.csdesw < 0.0)
    {
      err << "Warning: Csdesw = " << model_.csdesw << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.asd < 0.0)
    {
      err << "Warning: Asd = " << model_.asd << " should be within (0, 1).";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.rth0 < 0.0)
    {
      err << "Warning: Rth0 = " << model_.rth0 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.cth0 < 0.0)
    {
      err << "Warning: Cth0 = " << model_.cth0 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.rbody < 0.0)
    {
      err << "Warning: Rbody = " << model_.rbody << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.rbsh < 0.0)
    {
      err << "Warning: Rbsh = " << model_.rbsh << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->nigc <= 0.0)
    {
      err << "Fatal: nigc = " << paramPtr->nigc << " is non-positive.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      bsuccess = false;
    }

    if (paramPtr->poxedge <= 0.0)
    {
      err << "Fatal: poxedge = " << paramPtr->poxedge << " is non-positive.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      bsuccess = false;
    }

    if (paramPtr->pigcd <= 0.0)
    {
      err << "Fatal: pigcd = " << paramPtr->pigcd << " is non-positive.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      bsuccess = false;
    }

    if (model_.wth0 < 0.0)
    {
      err << "Warning:  Wth0 = " << model_.wth0 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.rhalo < 0.0)
    {
      err << "Warning:  Rhalo = " << model_.rhalo << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.ntox < 0.0)
    {
      err << "Warning:  Ntox = " << model_.ntox << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.toxref < 0.0)
    {
      err << "Warning:  Toxref = " << model_.toxref << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      bsuccess = false;
    }

    if (model_.ebg < 0.0)
    {
      err << "Warning:  Ebg = " << model_.ebg << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.vevb < 0.0)
    {
      err << "Warning:  Vevb = " << model_.vevb << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->alphaGB1 < 0.0)
    {
      err << "Warning:  AlphaGB1 = " << paramPtr->alphaGB1 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->betaGB1 < 0.0)
    {
      err << "Warning:  BetaGB1 = " << paramPtr->betaGB1 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.vgb1 < 0.0)
    {
      err << "Warning:  Vgb1 = " << model_.vgb1 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.vecb < 0.0)
    {
      err << "Warning:  Vecb = " << model_.vecb << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->alphaGB2 < 0.0)
    {
      err << "Warning:  AlphaGB2 = " << paramPtr->alphaGB2 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->betaGB2 < 0.0)
    {
      err << "Warning:  BetaGB2 = " << paramPtr->betaGB2 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.vgb2 < 0.0)
    {
      err << "Warning:  Vgb2 = " << model_.vgb2 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.toxqm <= 0.0)
    {
      err << "Fatal: Toxqm = " << model_.toxqm << " is not positive.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
      bsuccess = false;
    }

    if (model_.voxh < 0.0)
    {
      err << "Warning:  Voxh = " << model_.voxh << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.deltavox <= 0.0)
    {
      err << "Fatal: Deltavox = " << model_.deltavox << " is not positive.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.k1w1 < 0.0)
    {
      err << "Warning:  K1w1 = " << model_.k1w1 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.k1w2 < 0.0)
    {
      err << "Warning:  K1w2 = " << model_.k1w2 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.ketas < 0.0)
    {
      err << "Warning:  Ketas = " << model_.ketas << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.dwbc < 0.0)
    {
      err << "Warning:  Dwbc = " << model_.dwbc << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.beta0 < 0.0)
    {
      err << "Warning:  Beta0 = " << model_.beta0 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.beta1 < 0.0)
    {
      err << "Warning:  Beta1 = " << model_.beta1 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.beta2 < 0.0)
    {
      err << "Warning:  Beta2 = " << model_.beta2 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.tii < 0.0)
    {
      err << "Warning:  Tii = " << model_.tii << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.lii < 0.0)
    {
      err << "Warning:  Lii = " << model_.lii << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.sii1 < 0.0)
    {
      err << "Warning:  Sii1 = " << model_.sii1 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.sii2 < 0.0)
    {
      err << "Warning:  Sii2 = " << model_.sii1 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.siid < 0.0)
    {
      err << "Warning:  Siid = " << model_.siid << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.fbjtii < 0.0)
    {
      err << "Warning:  fbjtii = " << model_.fbjtii << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.vrec0 < 0.0)
    {
      err << "Warning:  Vrec0 = " << model_.vrec0 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.vtun0 < 0.0)
    {
      err << "Warning:  Vtun0 = " << model_.vtun0 << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.nbjt < 0.0)
    {
      err << "Warning:  Nbjt = " << model_.nbjt << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.aely < 0.0)
    {
      err << "Warning:  Aely = " << model_.aely << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.ahli < 0.0)
    {
      err << "Warning:  Ahli = " << model_.ahli << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.rbody < 0.0)
    {
      err << "Warning:  Rbody = " << model_.rbody << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.rbsh < 0.0)
    {
      err << "Warning:  Rbsh = " << model_.rbsh << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->ntrecf < 0.0)
    {
      err << "Warning:  Ntrecf = " << paramPtr->ntrecf << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (paramPtr->ntrecr < 0.0)
    {
      err << "Warning:  Ntrecr = " << paramPtr->ntrecr << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.tcjswg < 0.0)
    {
      err << "Warning:  Tcjswg = " << model_.tcjswg << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.tpbswg < 0.0)
    {
      err << "Warning:  Tpbswg = " << model_.tpbswg << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if ((model_.acde < 0.1) || (model_.acde > 1.6))
    {
      err << "Warning:  Acde = " << model_.acde << " is out of range.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if ((model_.moin < 5.0)||(model_.moin > 25.0))
    {
      err << "Warning:  Moin = " << model_.moin << " is out of range.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.dlbg < 0.0)
    {
      err << "Warning:  dlbg = " << model_.dlbg << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }


    if (model_.agidl < 0.0)
    {
      err << "Warning:  Agidl = " << model_.agidl << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.bgidl < 0.0)
    {
      err << "Warning:  Bgidl = " << model_.bgidl << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.ngidl < 0.0)
    {
      err << "Warning:  Ngidl = " << model_.ngidl << " is negative.";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.esatii < 0.0)
    {
      err << "Warning: Esatii = " << model_.esatii << " should be within (0, 1).";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }


    if (paramPtr->xj > model_.tsi)
    {
      err << "Warning: Xj = " << paramPtr->xj << " is thicker than Tsi = " << model_.tsi << ".";
      msg = err.str();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }

    if (model_.capMod < 2)
    {
      msg = "Warning: Warning: capMod < 2 is not supported by BSIM3SOI.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
    }
  }

  return  bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/02
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;

  if( icVDSGiven )
  {
    (extData.currStoVectorRawPtr)[li_store_vds] = icVDS;
    (extData.nextStoVectorRawPtr)[li_store_vds] = icVDS;
  }

  if( icVGSGiven )
  {
    (extData.currStoVectorRawPtr)[li_store_vgs] = icVGS;
    (extData.nextStoVectorRawPtr)[li_store_vgs] = icVGS;
  }

  if( icVBSGiven )
  {
    (extData.currStoVectorRawPtr)[li_store_vbs] = icVBS;
    (extData.nextStoVectorRawPtr)[li_store_vbs] = icVBS;
  }

  if( icVESGiven )
  {
    (extData.currStoVectorRawPtr)[li_store_ves] = icVES;
    (extData.nextStoVectorRawPtr)[li_store_ves] = icVES;
  }

  if( icVPSGiven )
  {
    (extData.currStoVectorRawPtr)[li_store_vps] = icVPS;
    (extData.nextStoVectorRawPtr)[li_store_vps] = icVPS;
  }

  return bsuccess;
}

// Class Model
//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{

  if (!given("DSUB"))
    dsub = drout;
  if (!given("XJ"))
    xj = tsi;
  if (!given("XDIF"))
    xdif = xbjt;
  if (!given("DWC"))
    dwc = Wint;
  if (!given("DLC"))
    dlc = Lint;
  if (!given("DLCIG"))
    dlcig = Lint;

  vcrit   = CONSTvt0 * log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
  factor1 = sqrt(CONSTEPSSI / CONSTEPSOX * tox);

  Vtm0 = CONSTKoverQ * tnom;
  Eg0  = CONSTEg0 - CONSTalphaEg * tnom * tnom / (tnom + CONSTbetaEg);
  ni   = CONSTNi0 * (tnom / CONSTREFTEMP) * sqrt(tnom / CONSTREFTEMP)
    * exp(21.5565981 - Eg0 / (2.0 * Vtm0));
  cbox = 3.453133e-11 / tbox;
  csi = 1.03594e-10 / tsi;

  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams(string param)
{
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                                  SolverState & ss1,
                                                  DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
    version                            ("3.2"),
    dtype                              (CONSTNMOS),
    cbox                               (0.0),
    csi                                (0.0),
    mobMod                             (0),
    capMod                             (0),
    binUnit                            (0),
    paramChk                           (0),
    model_l                            (0.0),
    model_w                            (0.0),
    Lmax                               (0.0),
    Lmin                               (0.0),
    Wmax                               (0.0),
    Wmin                               (0.0),
    dtoxcv                             (0.0),
    npeak                              (0.0),
    pdibl1                             (0.0),
    pdibl2                             (0.0),
    pdiblb                             (0.0),
    shMod                              (0.0),
    tbox                               (0.0),
    tsi                                (0.0),
    rth0                               (0.0),
    cth0                               (0.0),
    ngidl                              (0.0),
    agidl                              (0.0),
    bgidl                              (0.0),
    ndiode                             (0.0),
    xbjt                               (0.0),
    xdif                               (0.0),
    xrec                               (0.0),
    xtun                               (0.0),
    bodyJctGateSideGradingCoeff        (0.0),
    fnoiMod                            (0.0),
    tnoiMod                            (0.0),
    tnoia                              (0.0),
    tnoib                              (0.0),
    rnoia                              (0.0),
    rnoib                              (0.0),
    ntnoi                              (0.0),
    noif                               (0.0),
    k1w1                               (0.0),
    k1w2                               (0.0),
    ketas                              (0.0),
    dwbc                               (0.0),
    beta1                              (0.0),
    beta2                              (0.0),
    vdsatii0                           (0.0),
    tii                                (0.0),
    lii                                (0.0),
    sii0                               (0.0),
    sii1                               (0.0),
    sii2                               (0.0),
    siid                               (0.0),
    fbjtii                             (0.0),
    esatii                             (0.0),
    ntun                               (0.0),
    nrecf0                             (0.0),
    nrecr0                             (0.0),
    isbjt                              (0.0),
    isdif                              (0.0),
    isrec                              (0.0),
    istun                              (0.0),
    ln                                 (0.0),
    vrec0                              (0.0),
    vtun0                              (0.0),
    nbjt                               (0.0),
    lbjt0                              (0.0),
    ldif0                              (0.0),
    vabjt                              (0.0),
    aely                               (0.0),
    ahli                               (0.0),
    rbody                              (0.0),
    rbsh                               (0.0),
    cgeo                               (0.0),
    tt                                 (0.0),
    ndif                               (0.0),
    vsdfb                              (0.0),
    vsdth                              (0.0),
    csdmin                             (0.0),
    asd                                (0.0),
    csdesw                             (0.0),
    ntrecf                             (0.0),
    ntrecr                             (0.0),
    dlcb                               (0.0),
    fbody                              (0.0),
    delvt                              (0.0),
    kb1                                (0.0),
    dlbg                               (0.0),
    igbMod                             (0.0),
    igcMod                             (0.0),
    toxqm                              (0.0),
    wth0                               (0.0),
    rhalo                              (0.0),
    ntox                               (0.0),
    toxref                             (0.0),
    ebg                                (0.0),
    vevb                               (0.0),
    alphaGB1                           (0.0),
    betaGB1                            (0.0),
    vgb1                               (0.0),
    vecb                               (0.0),
    alphaGB2                           (0.0),
    betaGB2                            (0.0),
    vgb2                               (0.0),
    voxh                               (0.0),
    deltavox                           (0.0),
    aigc                               (0.0),
    bigc                               (0.0),
    cigc                               (0.0),
    aigsd                              (0.0),
    bigsd                              (0.0),
    cigsd                              (0.0),
    nigc                               (0.0),
    pigcd                              (0.0),
    poxedge                            (0.0),
    dlcig                              (0.0),
    soiMod                             (0),
    vbs0pd                             (0.0),
    vbs0fd                             (0.0),
    vbsa                               (0.0),
    nofffd                             (0.0),
    vofffd                             (0.0),
    k1b                                (0.0),
    k2b                                (0.0),
    dk2b                               (0.0),
    dvbd0                              (0.0),
    dvbd1                              (0.0),
    moinFD                             (0.0),
    rgateMod                           (0),
    bug1830fix                         (0),
    xrcrg1                             (0.0),
    xrcrg2                             (0.0),
    rshg                               (0.0),
    ngcon                              (0.0),
    xgw                                (0.0),
    xgl                                (0.0),
    lalphaGB1                          (0.0),
    lbetaGB1                           (0.0),
    lalphaGB2                          (0.0),
    lbetaGB2                           (0.0),
    lndif                              (0.0),
    lntrecf                            (0.0),
    lntrecr                            (0.0),
    lxbjt                              (0.0),
    lxdif                              (0.0),
    lxrec                              (0.0),
    lxtun                              (0.0),
    laigc                              (0.0),
    lbigc                              (0.0),
    lcigc                              (0.0),
    laigsd                             (0.0),
    lbigsd                             (0.0),
    lcigsd                             (0.0),
    lnigc                              (0.0),
    lpigcd                             (0.0),
    lpoxedge                           (0.0),
    lnpeak                             (0.0),
    lk1w1                              (0.0),
    lk1w2                              (0.0),
    lketas                             (0.0),
    lpdibl1                            (0.0),
    lpdibl2                            (0.0),
    lpdiblb                            (0.0),
    lfbjtii                            (0.0),
    lbeta1                             (0.0),
    lbeta2                             (0.0),
    lvdsatii0                          (0.0),
    llii                               (0.0),
    lesatii                            (0.0),
    lsii0                              (0.0),
    lsii1                              (0.0),
    lsii2                              (0.0),
    lsiid                              (0.0),
    lkb1                               (0.0),
    lagidl                             (0.0),
    lbgidl                             (0.0),
    lngidl                             (0.0),
    lntun                              (0.0),
    lndiode                            (0.0),
    lnrecf0                            (0.0),
    lnrecr0                            (0.0),
    lisbjt                             (0.0),
    lisdif                             (0.0),
    lisrec                             (0.0),
    listun                             (0.0),
    lvrec0                             (0.0),
    lvtun0                             (0.0),
    lnbjt                              (0.0),
    llbjt0                             (0.0),
    lvabjt                             (0.0),
    laely                              (0.0),
    lahli                              (0.0),
    lvsdfb                             (0.0),
    lvsdth                             (0.0),
    ldelvt                             (0.0),
    lxrcrg1                            (0.0),
    lxrcrg2                            (0.0),
    walphaGB1                          (0.0),
    wbetaGB1                           (0.0),
    walphaGB2                          (0.0),
    wbetaGB2                           (0.0),
    wndif                              (0.0),
    wntrecf                            (0.0),
    wntrecr                            (0.0),
    wxbjt                              (0.0),
    wxdif                              (0.0),
    wxrec                              (0.0),
    wxtun                              (0.0),
    waigc                              (0.0),
    wbigc                              (0.0),
    wcigc                              (0.0),
    waigsd                             (0.0),
    wbigsd                             (0.0),
    wcigsd                             (0.0),
    wnigc                              (0.0),
    wpigcd                             (0.0),
    wpoxedge                           (0.0),
    wnpeak                             (0.0),
    wk1w1                              (0.0),
    wk1w2                              (0.0),
    wkb1                               (0.0),
    wketas                             (0.0),
    wpdibl1                            (0.0),
    wpdibl2                            (0.0),
    wpdiblb                            (0.0),
    wfbjtii                            (0.0),
    wbeta1                             (0.0),
    wbeta2                             (0.0),
    wvdsatii0                          (0.0),
    wlii                               (0.0),
    wesatii                            (0.0),
    wsii0                              (0.0),
    wsii1                              (0.0),
    wsii2                              (0.0),
    wsiid                              (0.0),
    wagidl                             (0.0),
    wbgidl                             (0.0),
    wngidl                             (0.0),
    wntun                              (0.0),
    wndiode                            (0.0),
    wnrecf0                            (0.0),
    wnrecr0                            (0.0),
    wisbjt                             (0.0),
    wisdif                             (0.0),
    wisrec                             (0.0),
    wistun                             (0.0),
    wvrec0                             (0.0),
    wvtun0                             (0.0),
    wnbjt                              (0.0),
    wlbjt0                             (0.0),
    wvabjt                             (0.0),
    waely                              (0.0),
    wahli                              (0.0),
    wvsdfb                             (0.0),
    wvsdth                             (0.0),
    wdelvt                             (0.0),
    wxrcrg1                            (0.0),
    wxrcrg2                            (0.0),
    palphaGB1                          (0.0),
    pbetaGB1                           (0.0),
    palphaGB2                          (0.0),
    pbetaGB2                           (0.0),
    pndif                              (0.0),
    pntrecf                            (0.0),
    pntrecr                            (0.0),
    pxbjt                              (0.0),
    pxdif                              (0.0),
    pxrec                              (0.0),
    pxtun                              (0.0),
    paigc                              (0.0),
    pbigc                              (0.0),
    pcigc                              (0.0),
    paigsd                             (0.0),
    pbigsd                             (0.0),
    pcigsd                             (0.0),
    pnigc                              (0.0),
    ppigcd                             (0.0),
    ppoxedge                           (0.0),
    pnpeak                             (0.0),
    pk1w1                              (0.0),
    pk1w2                              (0.0),
    pkb1                               (0.0),
    pketas                             (0.0),
    ppdibl1                            (0.0),
    ppdibl2                            (0.0),
    ppdiblb                            (0.0),
    pfbjtii                            (0.0),
    pbeta1                             (0.0),
    pbeta2                             (0.0),
    pvdsatii0                          (0.0),
    plii                               (0.0),
    pesatii                            (0.0),
    psii0                              (0.0),
    psii1                              (0.0),
    psii2                              (0.0),
    psiid                              (0.0),
    pagidl                             (0.0),
    pbgidl                             (0.0),
    pngidl                             (0.0),
    pntun                              (0.0),
    pndiode                            (0.0),
    pnrecf0                            (0.0),
    pnrecr0                            (0.0),
    pisbjt                             (0.0),
    pisdif                             (0.0),
    pisrec                             (0.0),
    pistun                             (0.0),
    pvrec0                             (0.0),
    pvtun0                             (0.0),
    pnbjt                              (0.0),
    plbjt0                             (0.0),
    pvabjt                             (0.0),
    paely                              (0.0),
    pahli                              (0.0),
    pvsdfb                             (0.0),
    pvsdth                             (0.0),
    pdelvt                             (0.0),
    pxrcrg1                            (0.0),
    pxrcrg2                            (0.0),
    npeakGiven                         (false),
    csdminGiven                        (false),
    vsdthGiven                         (false),
    vsdfbGiven                         (false),
    gamma1Given                        (false),
    gamma2Given                        (false),
    vbxGiven                           (false),
    vbmGiven                           (false),
    xtGiven                            (false),
    k1Given                            (false),
    k2Given                            (false),
    vcrit                              (0.0),
    vtm                                (0.0),
    tox                                (0.0),
    toxm                               (0.0),
    cdsc                               (0.0),
    cdscb                              (0.0),
    cdscd                              (0.0),
    cit                                (0.0),
    nfactor                            (0.0),
    xj                                 (0.0),
    vsat                               (0.0),
    at                                 (0.0),
    a0                                 (0.0),
    ags                                (0.0),
    a1                                 (0.0),
    a2                                 (0.0),
    keta                               (0.0),
    nsub                               (0.0),
    ngate                              (0.0),
    gamma1                             (0.0),
    gamma2                             (0.0),
    vbx                                (0.0),
    vbm                                (0.0),
    xt                                 (0.0),
    k1                                 (0.0),
    kt1                                (0.0),
    kt1l                               (0.0),
    kt2                                (0.0),
    k2                                 (0.0),
    k3                                 (0.0),
    k3b                                (0.0),
    w0                                 (0.0),
    nlx                                (0.0),
    dvt0                               (0.0),
    dvt1                               (0.0),
    dvt2                               (0.0),
    dvt0w                              (0.0),
    dvt1w                              (0.0),
    dvt2w                              (0.0),
    drout                              (0.0),
    dsub                               (0.0),
    vth0                               (0.0),
    ua                                 (0.0),
    ua1                                (0.0),
    ub                                 (0.0),
    ub1                                (0.0),
    uc                                 (0.0),
    uc1                                (0.0),
    u0                                 (0.0),
    ute                                (0.0),
    voff                               (0.0),
    delta                              (0.0),
    rdsw                               (0.0),
    prwg                               (0.0),
    prwb                               (0.0),
    prt                                (0.0),
    eta0                               (0.0),
    etab                               (0.0),
    pclm                               (0.0),
    pvag                               (0.0),
    wr                                 (0.0),
    dwg                                (0.0),
    dwb                                (0.0),
    b0                                 (0.0),
    b1                                 (0.0),
    alpha0                             (0.0),
    beta0                              (0.0),
    cgsl                               (0.0),
    cgdl                               (0.0),
    ckappa                             (0.0),
    cf                                 (0.0),
    clc                                (0.0),
    cle                                (0.0),
    dwc                                (0.0),
    dlc                                (0.0),
    noff                               (0.0),
    acde                               (0.0),
    moin                               (0.0),
    tcjswg                             (0.0),
    tpbswg                             (0.0),
    lcdsc                              (0.0),
    lcdscb                             (0.0),
    lcdscd                             (0.0),
    lcit                               (0.0),
    lnfactor                           (0.0),
    lxj                                (0.0),
    lvsat                              (0.0),
    lat                                (0.0),
    la0                                (0.0),
    lags                               (0.0),
    la1                                (0.0),
    la2                                (0.0),
    lketa                              (0.0),
    lnsub                              (0.0),
    lngate                             (0.0),
    lk1                                (0.0),
    lkt1                               (0.0),
    lkt1l                              (0.0),
    lkt2                               (0.0),
    lk2                                (0.0),
    lk3                                (0.0),
    lk3b                               (0.0),
    lw0                                (0.0),
    lnlx                               (0.0),
    ldvt0                              (0.0),
    ldvt1                              (0.0),
    ldvt2                              (0.0),
    ldvt0w                             (0.0),
    ldvt1w                             (0.0),
    ldvt2w                             (0.0),
    ldrout                             (0.0),
    ldsub                              (0.0),
    lvth0                              (0.0),
    lua                                (0.0),
    lua1                               (0.0),
    lub                                (0.0),
    lub1                               (0.0),
    luc                                (0.0),
    luc1                               (0.0),
    lu0                                (0.0),
    lute                               (0.0),
    lvoff                              (0.0),
    ldelta                             (0.0),
    lrdsw                              (0.0),
    lprwg                              (0.0),
    lprwb                              (0.0),
    lprt                               (0.0),
    leta0                              (0.0),
    letab                              (0.0),
    lpclm                              (0.0),
    lpvag                              (0.0),
    lwr                                (0.0),
    ldwg                               (0.0),
    ldwb                               (0.0),
    lb0                                (0.0),
    lb1                                (0.0),
    lalpha0                            (0.0),
    lbeta0                             (0.0),
    lcgsl                              (0.0),
    lcgdl                              (0.0),
    lckappa                            (0.0),
    lnoff                              (0.0),
    lacde                              (0.0),
    lmoin                              (0.0),
    wcdsc                              (0.0),
    wcdscb                             (0.0),
    wcdscd                             (0.0),
    wcit                               (0.0),
    wnfactor                           (0.0),
    wxj                                (0.0),
    wvsat                              (0.0),
    wat                                (0.0),
    wa0                                (0.0),
    wags                               (0.0),
    wa1                                (0.0),
    wa2                                (0.0),
    wketa                              (0.0),
    wnsub                              (0.0),
    wngate                             (0.0),
    wk1                                (0.0),
    wkt1                               (0.0),
    wkt1l                              (0.0),
    wkt2                               (0.0),
    wk2                                (0.0),
    wk3                                (0.0),
    wk3b                               (0.0),
    ww0                                (0.0),
    wnlx                               (0.0),
    wdvt0                              (0.0),
    wdvt1                              (0.0),
    wdvt2                              (0.0),
    wdvt0w                             (0.0),
    wdvt1w                             (0.0),
    wdvt2w                             (0.0),
    wdrout                             (0.0),
    wdsub                              (0.0),
    wvth0                              (0.0),
    wua                                (0.0),
    wua1                               (0.0),
    wub                                (0.0),
    wub1                               (0.0),
    wuc                                (0.0),
    wuc1                               (0.0),
    wu0                                (0.0),
    wute                               (0.0),
    wvoff                              (0.0),
    wdelta                             (0.0),
    wrdsw                              (0.0),
    wprwg                              (0.0),
    wprwb                              (0.0),
    wprt                               (0.0),
    weta0                              (0.0),
    wetab                              (0.0),
    wpclm                              (0.0),
    wpvag                              (0.0),
    wwr                                (0.0),
    wdwg                               (0.0),
    wdwb                               (0.0),
    wb0                                (0.0),
    wb1                                (0.0),
    walpha0                            (0.0),
    wbeta0                             (0.0),
    wcgsl                              (0.0),
    wcgdl                              (0.0),
    wckappa                            (0.0),
    wnoff                              (0.0),
    wacde                              (0.0),
    wmoin                              (0.0),
    pcdsc                              (0.0),
    pcdscb                             (0.0),
    pcdscd                             (0.0),
    pcit                               (0.0),
    pnfactor                           (0.0),
    pxj                                (0.0),
    pvsat                              (0.0),
    pat                                (0.0),
    pa0                                (0.0),
    pags                               (0.0),
    pa1                                (0.0),
    pa2                                (0.0),
    pketa                              (0.0),
    pnsub                              (0.0),
    pngate                             (0.0),
    pk1                                (0.0),
    pkt1                               (0.0),
    pkt1l                              (0.0),
    pkt2                               (0.0),
    pk2                                (0.0),
    pk3                                (0.0),
    pk3b                               (0.0),
    pw0                                (0.0),
    pnlx                               (0.0),
    pdvt0                              (0.0),
    pdvt1                              (0.0),
    pdvt2                              (0.0),
    pdvt0w                             (0.0),
    pdvt1w                             (0.0),
    pdvt2w                             (0.0),
    pdrout                             (0.0),
    pdsub                              (0.0),
    pvth0                              (0.0),
    pua                                (0.0),
    pua1                               (0.0),
    pub                                (0.0),
    pub1                               (0.0),
    puc                                (0.0),
    puc1                               (0.0),
    pu0                                (0.0),
    pute                               (0.0),
    pvoff                              (0.0),
    pdelta                             (0.0),
    prdsw                              (0.0),
    pprwg                              (0.0),
    pprwb                              (0.0),
    pprt                               (0.0),
    peta0                              (0.0),
    petab                              (0.0),
    ppclm                              (0.0),
    ppvag                              (0.0),
    pwr                                (0.0),
    pdwg                               (0.0),
    pdwb                               (0.0),
    pb0                                (0.0),
    pb1                                (0.0),
    palpha0                            (0.0),
    pbeta0                             (0.0),
    pcgsl                              (0.0),
    pcgdl                              (0.0),
    pckappa                            (0.0),
    pnoff                              (0.0),
    pacde                              (0.0),
    pmoin                              (0.0),
    tnom                               (0.0),
    cgso                               (0.0),
    cgdo                               (0.0),
    xpart                              (0.0),
    sheetResistance                    (0.0),
    GatesidewallJctPotential           (0.0),
    unitLengthGateSidewallJctCap       (0.0),
    Lint                               (0.0),
    Ll                                 (0.0),
    Llc                                (0.0),
    Lln                                (0.0),
    Lw                                 (0.0),
    Lwc                                (0.0),
    Lwn                                (0.0),
    Lwl                                (0.0),
    Lwlc                               (0.0),
    Wint                               (0.0),
    Wl                                 (0.0),
    Wlc                                (0.0),
    Wln                                (0.0),
    Ww                                 (0.0),
    Wwc                                (0.0),
    Wwn                                (0.0),
    Wwl                                (0.0),
    Wwlc                               (0.0),
    cox                                (0.0),
    factor1                            (0.0),
    oxideTrapDensityA                  (0.0),
    oxideTrapDensityB                  (0.0),
    oxideTrapDensityC                  (0.0),
    em                                 (0.0),
    ef                                 (0.0),
    af                                 (0.0),
    kf                                 (0.0),
    vth0Given                          (false),
    igbModGiven                        (false),
    Vtm0                               (0.0),
    Eg0                                (0.0),
    ni                                 (0.0)
{
  if (getType() != "")
  {
    if (getType() == "NMOS") {
      dtype = CONSTNMOS;
    }
    else if (getType() == "PMOS") {
      dtype = CONSTPMOS;
    }
    else
    {
      string msg = "Could not recognize the type for model ";
      msg += getName();
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

#ifdef Xyce_B3SOI_USE_DEFL
  if (!model_lGiven)
    model_l=devOptions.defl;
  if (!model_wGiven)
    model_w=devOptions.defw;
#endif

  if (!given("TNOM"))
    tnom = devOptions.tnom;
  if (!given("TOXM"))
    toxm = tox;      /* v3.2 */
  if (!given("TOXQM"))
    toxqm = tox;
  if (!given("CF"))
    cf = 2.0 * CONSTEPSOX / M_PI * log(1.0 + 0.4e-6 / tox);
  if (!vth0Given)
    vth0 = (dtype == CONSTNMOS) ? 0.7 : -0.7;
  if (!given("UC"))
    uc = (mobMod == 3) ? -0.0465 : -0.0465e-9;
  if (!given("UC1"))
    uc1 = (mobMod == 3) ? -0.056 : -0.056e-9;
  if (!given("U0"))
    u0 = (dtype == CONSTNMOS) ? 0.067 : 0.025;
  if (!given("VOXH"))
    voxh = 5.0;
  if (!given("DELTAVOX"))
    deltavox = 0.005;
  if (!given("AIGC"))
    aigc = (dtype == CONSTNMOS) ? 0.43 : 0.31;
  if (!given("BIGC"))
    bigc = (dtype == CONSTNMOS) ? 0.054 : 0.024;
  if (!given("CIGC"))
    cigc = (dtype == CONSTNMOS) ? 0.075 : 0.03;
  if (!given("AIGSD"))
    aigsd = (dtype == CONSTNMOS) ? 0.43 : 0.31;
  if (!given("BIGSD"))
    bigsd = (dtype == CONSTNMOS) ? 0.054 : 0.024;
  if (!given("CIGSD"))
    cigsd = (dtype == CONSTNMOS) ? 0.075 : 0.03;
  if (!given("NOIA"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityA = 6.25e41;
    else
      oxideTrapDensityA=6.188e40;
  }
  if (!given("NOIB"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityB = 3.125e26;
    else
      oxideTrapDensityB = 1.5e25;
  }

  // PMC  It is assumed tox does not change in any temperature
  //      interpolation done in Instance::updateTemperature()
  cox = 3.453133e-11 / tox;
  if (!given("CGDO"))
  {
    if (given("DLC") && (dlc > 0.0))
      cgdo = dlc * cox - cgdl ;
    else
      cgdo = 0.6 * xj * cox;
  }
  if (!given("CGSO"))
  {
    if (given("DLC") && (dlc > 0.0))
      cgso = dlc * cox - cgsl ;
    else
      cgso = 0.6 * xj * cox;
  }

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  if ((soiMod != 0) && (soiMod != 1) && (soiMod != 2) && (soiMod != 3))
  {
    soiMod = 0;
    string msg = "soiMod has been set to its default value: 0";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
  }
  else if ((rgateMod != 0) && (rgateMod != 1) && (rgateMod != 2) && (rgateMod != 3))
  {
    rgateMod = 0;
    string msg = "rgateMod has been set to its default value: 0";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
  }
  if ((fnoiMod != 0) && (fnoiMod != 1))
  {
    fnoiMod = 1;
    string msg = "fnoiMod has been set to default value: 1";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
  }
  if ((tnoiMod != 0) && (tnoiMod != 1))
  {
    tnoiMod = 0;
    string msg = "tnoiMod has been set to default value: 0";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING,msg);
  }
  if (GatesidewallJctPotential < 0.1)
  {
    GatesidewallJctPotential = 0.1;
    string msg = "Given pbswg is less than 0.1. Pbswg is set to 0.1.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
  }
  sizeDependParamList.clear();

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
Model::~Model ()
{
  list<SizeDependParam*>::iterator it_dpL =
   sizeDependParamList.begin();
  list<SizeDependParam*>::iterator end_dpL =
   sizeDependParamList.end();
  for( ; it_dpL != end_dpL; ++it_dpL )
    delete (*it_dpL);

  sizeDependParamList.clear ();

  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}


//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << endl;
  os << "    name     modelName  Parameters" << endl;

  for (i=0, iter=first; iter!=last; ++iter,++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;

  return os;
}


//----------------------------------------------------------------------------
// Function       : Model::clearTemperatureData
//
// Purpose        : This is mainly here to delete rid of the size
//                  dependent parameters, which are also temperature dependent.
//
// Special Notes  : This is called right before the circuit temperature is
//                  changed.
//
// Scope          : public
// Creator        : Eric R. Keiter, 9233, computation sciences
// Creation Date  : 03/08/2005
//----------------------------------------------------------------------------
bool Model::clearTemperatureData ()
{
  list<SizeDependParam*>::iterator it_dpL =
   sizeDependParamList.begin();

  list<SizeDependParam*>::iterator end_dpL =
   sizeDependParamList.end();

  for( ; it_dpL != end_dpL; ++it_dpL )
    delete (*it_dpL);

  sizeDependParamList.clear ();

  return true;
}

//-----------------------------------------------------------------------------
// MOSFET B3SOI Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/06/09
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;

  // vector<N_DEV_MOSFET_B3SOIModel*>::iterator iterM;
  // vector<N_DEV_MOSFET_B3SOIModel*>::iterator firstM = modelContainer_.begin();
  // vector<N_DEV_MOSFET_B3SOIModel*>::iterator lastM  = modelContainer_.end();

  // vector<N_DEV_MOSFET_B3SOIInstance*>::iterator iterI;
  // vector<N_DEV_MOSFET_B3SOIInstance*>::iterator firstI;
  // vector<N_DEV_MOSFET_B3SOIInstance*>::iterator lastI;


  // // first loop over the models:
  // for (iterM=firstM; iterM!=lastM; ++iterM)
  // {
  //   // loop over the instances for this model.
  //   firstI = (*iterM)->instanceContainer.begin();
  //   lastI   = (*iterM)->instanceContainer.end();

  //   for (iterI=firstI; iterI!=lastI; ++iterI)
  //   {
  //     N_DEV_MOSFET_B3SOIInstance & mi = *(*iterI);

  for (ModelMap::const_iterator it_model = getModelMap().begin(); it_model != getModelMap().end(); ++it_model)
  {
    for (InstanceVector::const_iterator it = (*it_model).second->getInstanceVector().begin(); it != (*it_model).second->getInstanceVector().end(); ++it)
    {
      Instance & mi = *(*it);

      double * stoVec = mi.extData.nextStoVectorRawPtr;

      // save voltage drops
      bool btmp = mi.updateIntermediateVars ();
      bsuccess = bsuccess && btmp;

      // voltage drops:
      stoVec[mi.li_store_vbd] = mi.vbd;
      stoVec[mi.li_store_vbs] = mi.vbs;
      stoVec[mi.li_store_vgs] = mi.vgs;
      stoVec[mi.li_store_vds] = mi.vds;
      stoVec[mi.li_store_ves] = mi.ves;
      stoVec[mi.li_store_vps] = mi.vps;

      stoVec[mi.li_store_vgp] = mi.Vgp;
      stoVec[mi.li_store_vd] = mi.Vdp;
      stoVec[mi.li_store_vs] = mi.Vsp;
      stoVec[mi.li_store_vp] = mi.Vp;
      stoVec[mi.li_store_ve] = mi.Ve;
      stoVec[mi.li_store_vg] = mi.Vg;
      stoVec[mi.li_store_vgm] = mi.Vgm;
      stoVec[mi.li_store_deltemp] = mi.delTemp;

      stoVec[mi.li_store_vges] = mi.vges;
      stoVec[mi.li_store_vgms] = mi.vgms;

      // intrinsic capacitors:
      staVec[mi.li_state_qb] = mi.qb;
      staVec[mi.li_state_qg] = mi.qg;
      staVec[mi.li_state_qd] = mi.qd;
      staVec[mi.li_state_qe] = mi.qe;
      staVec[mi.li_state_qgmid] = mi.qgmid;
      staVec[mi.li_state_qth] = mi.qth;

      // if this is the first newton step of the first time step
      // of the transient simulation, we need to enforce that the
      // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
      // compatibility.  ERK.

      // Note:  I think this kind of thing is enforced (or should be enforced,
      // anyway) at the time integration level.  So I'm not sure this step is
      // really needed, at least for new-DAE.  Derivatives out of the DCOP
      // are supposed to be zero at the first newton step.

      if (!(getSolverState().dcopFlag) && getSolverState().initTranFlag && getSolverState().newtonIter==0)
      {
        double * currStaVec = mi.extData.currStaVectorRawPtr;

        // intrinsic capacitors:
        currStaVec[mi.li_state_qb] = mi.qb;
        currStaVec[mi.li_state_qg] = mi.qg;
        currStaVec[mi.li_state_qd] = mi.qd;
        currStaVec[mi.li_state_qe] = mi.qe;
        currStaVec[mi.li_state_qgmid] = mi.qgmid;
        currStaVec[mi.li_state_qth] = mi.qth;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/06/09
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & mi = *(*it);

    double * dFdxdVp = mi.extData.dFdxdVpVectorRawPtr;
    double * dQdxdVp = mi.extData.dQdxdVpVectorRawPtr;

    // load F:

    double Coef_f_body=0.0;
    double Coef_f_extBody=0.0;
    double Coef_f_gate=0.0;
    double Coef_f_gatePrime=0.0;
    double Coef_f_gateMid=0.0;
    double Coef_f_drain=0.0;
    double Coef_f_drainPrime=0.0;
    double Coef_f_source=0.0;
    double Coef_f_sourcePrime=0.0;
    double Coef_f_substrate=0.0;
    double Coef_f_temp=0.0;

    mi.Gmin = getDeviceOptions().gmin * 1e-6;
    mi.geltd = mi.grgeltd;

    // I have elected to add in the Gmin-based Jdxp terms right here.
    // In adding them, I've made use of the fact that whether we need voltage
    // limiting or not, we have to add Gmin*somedrop.  If we need voltage
    // limiting, we have to add in Gmin*(somedrop-somedrop_orig) with the
    // opposite sign.  I've simplified that down to the single term that
    // remains, which is correct even if voltage limiting is not needed
    // --- TVR

    double iGmin_bs = mi.model_.dtype*(mi.Gmin*mi.vbs_orig);
    double iGmin_gd = mi.model_.dtype*(mi.Gmin*mi.vgd_orig);

    if (mi.soiMod != 2)
    {
      Coef_f_body  -= (mi.model_.dtype*(mi.ceqbody)
                    + iGmin_bs)* mi.numberParallel;
    }

    Coef_f_gatePrime  -= (mi.model_.dtype*(mi.ceqgate - mi.ceqgcrg)
                        + mi.Igtoteq + iGmin_gd)* mi.numberParallel;

    Coef_f_drainPrime  += (mi.model_.dtype*(mi.ceqbd)  - mi.cdreq + mi.Idtoteq
                        + mi.Idrain + iGmin_gd) * mi.numberParallel;

    Coef_f_sourcePrime  += (mi.cdreq
                          + mi.model_.dtype*(mi.ceqbs)
                          + mi.Istoteq + mi.Isource + iGmin_bs) * mi.numberParallel;

    if (mi.rgateMod == 2)
    {
      Coef_f_gate -= mi.model_.dtype*mi.ceqgcrg * mi.numberParallel;
    }
    else if (mi.rgateMod == 3)
    {
      Coef_f_gateMid -= mi.model_.dtype*(mi.ceqgcrg) * mi.numberParallel;
    }

    if (mi.bodyMod == 1)
    {
      Coef_f_extBody += mi.model_.dtype*mi.ceqbodcon * mi.numberParallel;
    }

    if (mi.selfheat)
    {
      Coef_f_temp -= (mi.ceqth + mi.Ith) * mi.numberParallel;
    }

    ////////////////////////////////////////////////////////////////////////
    // This next section deals with linear resistor currents (mostly) that
    // would not be part of the spice3f5 RHS load.
    //
    // These include the source and drain load resistors, the body tie
    // resistor and the various options for gate resistors.
    //
    // Isource and Idrain correspond to the linear load resistors.

    if (mi.sNodePrime == 1) Coef_f_source  -= mi.Isource * mi.numberParallel;
    if (mi.dNodePrime == 1) Coef_f_drain -= mi.Idrain * mi.numberParallel;

    // Note: extra body nodes were already handled, above.

    // Handle extra gate nodes, if they exist.
    // mi.rgateMod==0   no gate resistor.
    // mi.rgateMod==1   linear gate resistor
    // mi.rgateMod==2   nonlinear gate resistor
    // mi.rgateMod==3   2 gate resistors, in series.
    //
    if (mi.rgateMod > 0)
    {
      Coef_f_gate -= mi.Igate * mi.numberParallel;
      if (mi.rgateMod == 3)
      {
        Coef_f_gateMid -= (mi.IgateMid - mi.Igate) * mi.numberParallel;
        Coef_f_gatePrime += mi.IgateMid * mi.numberParallel;
      }
      else
      {
        Coef_f_gatePrime += mi.Igate * mi.numberParallel;
      }
    }

    // Now place all the calculated values in the RHS vector.
    if(mi.li_Body != -1)
    {
      fVec[mi.li_Body       ] -= Coef_f_body;
    }
    if(mi.li_ExtBody != -1)
    {
      fVec[mi.li_ExtBody    ] -= Coef_f_extBody;
    }
    fVec[mi.li_Gate       ] -= Coef_f_gate;
    fVec[mi.li_GatePrime  ] -= Coef_f_gatePrime;
    fVec[mi.li_GateMid    ] -= Coef_f_gateMid;
    fVec[mi.li_Drain      ] -= Coef_f_drain;
    fVec[mi.li_DrainPrime ] -= Coef_f_drainPrime;
    fVec[mi.li_Source     ] -= Coef_f_source;
    fVec[mi.li_SourcePrime] -= Coef_f_sourcePrime;
    fVec[mi.li_Substrate  ] -= Coef_f_substrate;

    if(mi.li_Temperature != -1)
    {
      fVec[mi.li_Temperature] -= Coef_f_temp;
    }

    if( mi.loadLeadCurrent )
    {

      if( mi.dNodePrime != 1 )
      {
        // drain prime is the same as drain so add in its contribution
        storeLeadF[mi.li_store_dev_id] = -Coef_f_drain - Coef_f_drainPrime;
      }
      else
      {
        storeLeadF[mi.li_store_dev_id] = -Coef_f_drain;
      }
      if( mi.gNodePrime != 1 )
      {
        // gate prime is the same as gate so add in its contribution
        storeLeadF[mi.li_store_dev_ig] = -Coef_f_gate - Coef_f_gatePrime;
      }
      else
      {
        storeLeadF[mi.li_store_dev_ig] = -Coef_f_gate;
      }
      if( mi.sNodePrime != 1 )
      {
        // source prime is the same as source so add in its contribution
        storeLeadF[mi.li_store_dev_is] = -Coef_f_source - Coef_f_sourcePrime;
      }
      else
      {
        storeLeadF[mi.li_store_dev_is] = -Coef_f_source;
      }

      storeLeadF[mi.li_store_dev_ie] = -Coef_f_substrate;
      if(mi.li_Body != -1)
      {
        storeLeadF[mi.li_store_dev_ib] = -Coef_f_body;
      }
      else
      {
        storeLeadF[mi.li_store_dev_ib] = 0.0;
      }
    }

    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      if ( mi.icVDSGiven )
      {
        double coef = (mi.extData.nextSolVectorRawPtr)[mi.li_Ids];
        fVec[mi.li_Drain] += coef;
        fVec[mi.li_Source] += -coef;
        if( mi.loadLeadCurrent )
        {
          storeLeadF[mi.li_store_dev_id] = coef;
          storeLeadF[mi.li_store_dev_is] = -coef;
        }
        double cVs = (mi.extData.nextSolVectorRawPtr)[mi.li_Source];
        double cVd = (mi.extData.nextSolVectorRawPtr)[mi.li_Drain];
        fVec[mi.li_Ids] += (cVd - cVs - mi.icVDS);
      }
      if ( mi.icVGSGiven )
      {
        double coef = (mi.extData.nextSolVectorRawPtr)[mi.li_Igs];
        fVec[mi.li_Gate] += coef;
        fVec[mi.li_Source] += -coef;
        if( mi.loadLeadCurrent )
        {
          storeLeadF[mi.li_store_dev_ig] = coef;
          storeLeadF[mi.li_store_dev_is] = -coef;
        }
        double cVs = (mi.extData.nextSolVectorRawPtr)[mi.li_Source];
        double cVg = (mi.extData.nextSolVectorRawPtr)[mi.li_Gate];
        fVec[mi.li_Igs] += (cVg - cVs - mi.icVGS);
      }
      if ( mi.icVBSGiven )
      {
        double coef = (mi.extData.nextSolVectorRawPtr)[mi.li_Ibs];
        fVec[mi.li_Body] += coef;
        fVec[mi.li_Source] += -coef;
        if( mi.loadLeadCurrent )
        {
          storeLeadF[mi.li_store_dev_ib] = coef;
          storeLeadF[mi.li_store_dev_is] = -coef;
        }
        double cVs = (mi.extData.nextSolVectorRawPtr)[mi.li_Source];
        double cVb = (mi.extData.nextSolVectorRawPtr)[mi.li_Body];
        fVec[mi.li_Ibs] += (cVb - cVs - mi.icVBS);
      }
      if ( mi.icVESGiven )
      {
        double coef = (mi.extData.nextSolVectorRawPtr)[mi.li_Ies];
        fVec[mi.li_Substrate] += coef;
        fVec[mi.li_Source] += -coef;
        if( mi.loadLeadCurrent )
        {
          storeLeadF[mi.li_store_dev_ie] = coef;
          storeLeadF[mi.li_store_dev_is] = -coef;
        }
        double cVs = (mi.extData.nextSolVectorRawPtr)[mi.li_Source];
        double cVe = (mi.extData.nextSolVectorRawPtr)[mi.li_Substrate];
        fVec[mi.li_Ies] += (cVe - cVs - mi.icVES);
      }
      if ( mi.icVPSGiven )
      {
        double coef = (mi.extData.nextSolVectorRawPtr)[mi.li_Ips];
        fVec[mi.li_ExtBody] += coef;
        fVec[mi.li_Source] += -coef;
        if( mi.loadLeadCurrent )
        {
          storeLeadF[mi.li_store_dev_ib] = coef;
          storeLeadF[mi.li_store_dev_is] = -coef;
        }
        double cVs = (mi.extData.nextSolVectorRawPtr)[mi.li_Source];
        double cVp = (mi.extData.nextSolVectorRawPtr)[mi.li_ExtBody];
        fVec[mi.li_Ies] += (cVp - cVs - mi.icVPS);
      }
    }

    // Set up the Jdxp vector:
    double Coef_f_body_Jdxp        = 0.0;
    double Coef_f_extBody_Jdxp     = 0.0;
    double Coef_f_gate_Jdxp        = 0.0;
    double Coef_f_gatePrime_Jdxp   = 0.0;
    double Coef_f_gateMid_Jdxp     = 0.0;
    double Coef_f_drain_Jdxp       = 0.0;
    double Coef_f_drainPrime_Jdxp  = 0.0;
    double Coef_f_source_Jdxp      = 0.0;
    double Coef_f_sourcePrime_Jdxp = 0.0;
    double Coef_f_substrate_Jdxp   = 0.0;
    double Coef_f_temp_Jdxp        = 0.0;

    double iGmin_bs_Jdxp = mi.model_.dtype*mi.Gmin*(mi.vbs-mi.vbs_orig);
    double iGmin_gd_Jdxp = mi.model_.dtype*mi.Gmin*(mi.vgd-mi.vgd_orig);

    if (getDeviceOptions().voltageLimiterFlag)
    {
      if (mi.soiMod != 2)
      {
        Coef_f_body_Jdxp  -=
        (mi.model_.dtype*(mi.ceqbody_Jdxp )
          - iGmin_bs_Jdxp)* mi.numberParallel;
      }

      Coef_f_gatePrime_Jdxp  -=
          (mi.model_.dtype*(mi.ceqgate_Jdxp - mi.ceqgcrg_Jdxp)
          + mi.Igtoteq_Jdxp - iGmin_gd_Jdxp)* mi.numberParallel;

      Coef_f_drainPrime_Jdxp  +=
          (mi.model_.dtype*(mi.ceqbd_Jdxp)
          - mi.cdreq_Jdxp + mi.Idtoteq_Jdxp
          + mi.Idrain_Jdxp - iGmin_gd_Jdxp) * mi.numberParallel;

      Coef_f_sourcePrime_Jdxp  += (mi.cdreq_Jdxp
        + mi.model_.dtype*(mi.ceqbs_Jdxp)
                                + mi.Istoteq_Jdxp
                                + mi.Isource_Jdxp - iGmin_bs_Jdxp
                            ) * mi.numberParallel;

      Coef_f_substrate_Jdxp  -= mi.model_.dtype*mi.ceqqe_Jdxp * mi.numberParallel;

      if (mi.rgateMod == 2)
      {
        Coef_f_gate_Jdxp -=
          mi.model_.dtype*mi.ceqgcrg_Jdxp * mi.numberParallel;
      }
      else if (mi.rgateMod == 3)
      {
        Coef_f_gateMid_Jdxp -=
        mi.model_.dtype*(mi.ceqgcrg_Jdxp) * mi.numberParallel;
      }

      if (mi.bodyMod == 1)
      {
        Coef_f_extBody_Jdxp += mi.model_.dtype*mi.ceqbodcon_Jdxp * mi.numberParallel;
      }

      if (mi.selfheat)
      {
        Coef_f_temp_Jdxp -= ( mi.ceqth_Jdxp ) * mi.numberParallel;
      }

      if (mi.sNodePrime == 1) Coef_f_source_Jdxp  -= mi.Isource_Jdxp * mi.numberParallel;
      if (mi.dNodePrime == 1) Coef_f_drain_Jdxp -= mi.Idrain_Jdxp * mi.numberParallel;

      if (mi.rgateMod > 0)
      {
        Coef_f_gate_Jdxp -= mi.Igate_Jdxp * mi.numberParallel;
        if (mi.rgateMod == 3)
        {
          Coef_f_gateMid_Jdxp -= (mi.IgateMid_Jdxp - mi.Igate_Jdxp) * mi.numberParallel;
          Coef_f_gatePrime_Jdxp += mi.IgateMid_Jdxp * mi.numberParallel;
        }
        else
        {
          Coef_f_gatePrime_Jdxp += mi.Igate_Jdxp * mi.numberParallel;
        }
      }

      if(mi.li_Body != -1)
      {
        dFdxdVp[mi.li_Body       ] += Coef_f_body_Jdxp;
      }
      if(mi.li_ExtBody != -1)
      {
        dFdxdVp[mi.li_ExtBody    ] += Coef_f_extBody_Jdxp;
      }
      dFdxdVp[mi.li_Gate       ] += Coef_f_gate_Jdxp;
      dFdxdVp[mi.li_GatePrime  ] += Coef_f_gatePrime_Jdxp;
      dFdxdVp[mi.li_GateMid    ] += Coef_f_gateMid_Jdxp;
      dFdxdVp[mi.li_Drain      ] += Coef_f_drain_Jdxp;
      dFdxdVp[mi.li_DrainPrime ] += Coef_f_drainPrime_Jdxp;
      dFdxdVp[mi.li_Source     ] += Coef_f_source_Jdxp;
      dFdxdVp[mi.li_SourcePrime] += Coef_f_sourcePrime_Jdxp;
      dFdxdVp[mi.li_Substrate  ] += Coef_f_substrate_Jdxp;
      if(mi.li_Temperature != -1)
      {
        dFdxdVp[mi.li_Temperature] += Coef_f_temp_Jdxp;
      }
    }

    // load Q:

    mi.auxChargeCalculations ();

    double Coef_q_body=0.0;
    double Coef_q_extBody=0.0;
    double Coef_q_gate=0.0;
    double Coef_q_gatePrime=0.0;
    double Coef_q_gateMid=0.0;
    double Coef_q_drainPrime=0.0;
    double Coef_q_sourcePrime=0.0;
    double Coef_q_substrate=0.0;
    double Coef_q_temp=0.0;

    if (mi.soiMod != 2)
    {
      Coef_q_body  -= (mi.model_.dtype*(mi.Qeqqb))* mi.numberParallel;
    }

    Coef_q_gatePrime  -= (mi.model_.dtype*(mi.Qeqqg))* mi.numberParallel;

    Coef_q_drainPrime  += (mi.model_.dtype*(-mi.Qeqqd)) * mi.numberParallel;

    Coef_q_sourcePrime  += ( + mi.model_.dtype*(+ mi.Qeqqg + mi.Qeqqb
                      + mi.Qeqqd + mi.Qeqqe + mi.Qeqqgmid)) * mi.numberParallel;

    Coef_q_substrate  -= mi.model_.dtype*mi.Qeqqe * mi.numberParallel;

    if (mi.rgateMod == 3)
    {
      Coef_q_gateMid -= mi.model_.dtype*(mi.Qeqqgmid) * mi.numberParallel;
    }

    if (mi.selfheat)
    {
      Coef_q_temp -= (mi.Qeqqth) * mi.numberParallel;
    }

    ////////////////////////////////////////////////////////////////////////
    // Now place all the calculated values in the Q vector.

    if(mi.li_Body != -1)
    {
      qVec[mi.li_Body       ] -= Coef_q_body;
    }
    if(mi.li_ExtBody != -1)
    {
      qVec[mi.li_ExtBody    ] -= Coef_q_extBody;
    }
    qVec[mi.li_Gate       ] -= Coef_q_gate;
    qVec[mi.li_GatePrime  ] -= Coef_q_gatePrime;
    qVec[mi.li_GateMid    ] -= Coef_q_gateMid;
    qVec[mi.li_DrainPrime ] -= Coef_q_drainPrime;
    qVec[mi.li_SourcePrime] -= Coef_q_sourcePrime;
    qVec[mi.li_Substrate  ] -= Coef_q_substrate;
    if(mi.li_Temperature != -1)
    {
      qVec[mi.li_Temperature] -= Coef_q_temp;
    }

    if( mi.loadLeadCurrent )
    {

      if( mi.dNodePrime != 1 )
      {
        // drain prime is the same as drain so add in its contribution
        storeLeadQ[mi.li_store_dev_id] = -Coef_q_drainPrime;
      }
      else
      {
        storeLeadQ[mi.li_store_dev_id] = 0;
      }
      if( mi.gNodePrime != 1 )
      {
        // gate prime is the same as gate so add in its contribution
        storeLeadQ[mi.li_store_dev_ig] = -Coef_q_gate - Coef_q_gatePrime;
      }
      else
      {
        storeLeadQ[mi.li_store_dev_ig] = -Coef_q_gate;
      }
      if( mi.sNodePrime != 1 )
      {
        // source prime is the same as source so add in its contribution
        storeLeadQ[mi.li_store_dev_is] = -Coef_q_sourcePrime;
      }
      else
      {
        storeLeadQ[mi.li_store_dev_is] = 0.0;
      }
      storeLeadQ[mi.li_store_dev_ie] = -Coef_q_substrate;
      if(mi.li_Body != -1)
      {
        storeLeadQ[mi.li_store_dev_ib] = -Coef_q_body;
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // Voltage limiting section:
    double Coef_q_body_Jdxp        = 0.0;
    double Coef_q_extBody_Jdxp     = 0.0;
    double Coef_q_gate_Jdxp        = 0.0;
    double Coef_q_gatePrime_Jdxp   = 0.0;
    double Coef_q_gateMid_Jdxp     = 0.0;
    double Coef_q_drainPrime_Jdxp  = 0.0;
    double Coef_q_sourcePrime_Jdxp = 0.0;
    double Coef_q_substrate_Jdxp   = 0.0;
    double Coef_q_temp_Jdxp        = 0.0;

    if (getDeviceOptions().voltageLimiterFlag)
    {
      if (mi.soiMod != 2)
      {
        Coef_q_body_Jdxp  -=
        (mi.model_.dtype*(mi.Qeqqb_Jdxp))* mi.numberParallel;
      }

      Coef_q_gatePrime_Jdxp  -=
        (mi.model_.dtype*(mi.Qeqqg_Jdxp))* mi.numberParallel;

      Coef_q_drainPrime_Jdxp  +=
        (mi.model_.dtype*(- mi.Qeqqd_Jdxp)) * mi.numberParallel;

      Coef_q_sourcePrime_Jdxp  += (+mi.model_.dtype*(mi.Qeqqg_Jdxp + mi.Qeqqb_Jdxp
                          + mi.Qeqqd_Jdxp + mi.Qeqqe_Jdxp + mi.Qeqqgmid_Jdxp)
                          ) * mi.numberParallel;

      Coef_q_substrate_Jdxp  -= mi.model_.dtype*mi.Qeqqe_Jdxp * mi.numberParallel;

      if (mi.rgateMod == 3)
      {
        Coef_q_gateMid_Jdxp -=
        mi.model_.dtype*(mi.Qeqqgmid_Jdxp) * mi.numberParallel;
      }

      if (mi.selfheat)
      {
        Coef_q_temp_Jdxp -= (mi.Qeqqth_Jdxp ) * mi.numberParallel;
      }

      if(mi.li_Body != -1)
      {
        dQdxdVp[mi.li_Body       ] += Coef_q_body_Jdxp;
      }
      if(mi.li_ExtBody != -1)
      {
        dQdxdVp[mi.li_ExtBody    ] += Coef_q_extBody_Jdxp;
      }
      dQdxdVp[mi.li_Gate       ] += Coef_q_gate_Jdxp;
      dQdxdVp[mi.li_GatePrime  ] += Coef_q_gatePrime_Jdxp;
      dQdxdVp[mi.li_GateMid    ] += Coef_q_gateMid_Jdxp;
      dQdxdVp[mi.li_DrainPrime ] += Coef_q_drainPrime_Jdxp;
      dQdxdVp[mi.li_SourcePrime] += Coef_q_sourcePrime_Jdxp;
      dQdxdVp[mi.li_Substrate  ] += Coef_q_substrate_Jdxp;
      if(mi.li_Temperature != -1)
      {
        dQdxdVp[mi.li_Temperature] += Coef_q_temp_Jdxp;
      }
    }

  }

  return true;
}

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/06/09
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
#ifdef _OMP
#pragma omp parallel for
#endif
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & mi = *(*it);

    // F-matrix:
    if (mi.rgateMod == 1)
    {

      *mi.f_GateEquGateNodePtr
        += mi.geltd*mi.numberParallel;

      *mi.f_GatePrimeEquGateNodePtr
        -= mi.geltd*mi.numberParallel;

      *mi.f_GateEquGatePrimeNodePtr
        -= mi.geltd*mi.numberParallel;
      // It seems that this should be here, but it shouldn't.
      // The mi.geltd term is added into the G'-G' element later down in this
      // routine
      //
      //    *mi.f_GatePrime][mi.AGatePrimeEquGatePrimeNodePtr
      //      += mi.geltd*mi.numberParallel;
    }
    else if (mi.rgateMod == 2)
    {

      *mi.f_GateEquGateNodePtr
        += mi.gcrg_jac*mi.numberParallel;

      *mi.f_GateEquGatePrimeNodePtr
        += mi.gcrgg_jac*mi.numberParallel;

      *mi.f_GateEquDrainPrimeNodePtr
        += mi.gcrgd_jac*mi.numberParallel;

      *mi.f_GateEquSourcePrimeNodePtr
        += mi.gcrgs_jac*mi.numberParallel;

      *mi.f_GatePrimeEquGateNodePtr
        -= mi.gcrg_jac*mi.numberParallel;
      if (mi.soiMod != 2)
      {
        *mi.f_GateEquBodyNodePtr
        += mi.gcrgb_jac*mi.numberParallel;
      }
    }
    else if (mi.rgateMod == 3)
    {

      *mi.f_GateEquGateNodePtr
        += mi.geltd*mi.numberParallel;

      *mi.f_GateEquGateMidNodePtr
        += -mi.geltd*mi.numberParallel;

      *mi.f_GateMidEquGateNodePtr
        += -mi.geltd*mi.numberParallel;

      *mi.f_GateMidEquGateMidNodePtr
        += (mi.geltd + mi.gcrg_jac)*mi.numberParallel;

      *mi.f_GateMidEquDrainPrimeNodePtr
        += (mi.gcrgd_jac)*mi.numberParallel;

      *mi.f_GateMidEquGatePrimeNodePtr
        += mi.gcrgg_jac*mi.numberParallel;

      *mi.f_GateMidEquSourcePrimeNodePtr
        += (mi.gcrgs_jac)*mi.numberParallel;

      if (mi.soiMod != 2)
      {
        *mi.f_GateMidEquBodyNodePtr
        += mi.gcrgb_jac*mi.numberParallel;
      }


      *mi.f_GatePrimeEquGateMidNodePtr
        -= mi.gcrg_jac*mi.numberParallel;
    }
    if (mi.soiMod != 0)
    {

      *mi.f_DrainPrimeEquSubstrateNodePtr
        += (mi.Gme + mi.gddpe)*mi.numberParallel;

      *mi.f_SourcePrimeEquSubstrateNodePtr
        += (mi.gsspe - mi.Gme)*mi.numberParallel;
      if (mi.soiMod != 2)
      {
        *mi.f_GatePrimeEquSubstrateNodePtr
        += mi.gige_jac*mi.numberParallel;
        *mi.f_BodyEquSubstrateNodePtr
        -= mi.gige_jac*mi.numberParallel;
      }
    }

    if (mi.soiMod != 2)
    {
      if (mi.rgateMod == 0 || mi.rgateMod == 1)
      {
        *mi.f_GatePrimeEquBodyNodePtr
        -= (-mi.gigb_jac - mi.gIgtotb)*mi.numberParallel;
      }
      else
      {
        *mi.f_GatePrimeEquBodyNodePtr
        += (mi.gigb_jac +mi.gIgtotb - mi.gcrgb_jac)*mi.numberParallel;
      }

      *mi.f_DrainPrimeEquBodyNodePtr
        -= ((-mi.gddpb - mi.Gmbs) + mi.gIdtotb)*mi.numberParallel;

      *mi.f_SourcePrimeEquBodyNodePtr
        -= ((-mi.gsspb + mi.Gmbs) + mi.Gmin + mi.gIstotb)*mi.numberParallel;

      *mi.f_BodyEquSubstrateNodePtr
        += (mi.gbbe)*mi.numberParallel;

      *mi.f_BodyEquGatePrimeNodePtr
        += (-mi.gigg + mi.gbbg)*mi.numberParallel;

      *mi.f_BodyEquDrainPrimeNodePtr
        += (-mi.gigd_jac + mi.gbbdp)*mi.numberParallel;

      *mi.f_BodyEquSourcePrimeNodePtr
        += (mi.gbbsp - mi.Gmin - mi.gigs)*mi.numberParallel;

      *mi.f_BodyEquBodyNodePtr
        += (-mi.gigb_jac + mi.gbbb + mi.Gmin)*mi.numberParallel;
    }
    if (mi.rgateMod == 0)
    {

      *mi.f_GatePrimeEquGatePrimeNodePtr
        += (mi.gigg + mi.Gmin + mi.gIgtotg)*mi.numberParallel;

      *mi.f_GatePrimeEquDrainPrimeNodePtr
        += (mi.gigd_jac - mi.Gmin + mi.gIgtotd)*mi.numberParallel;

      *mi.f_GatePrimeEquSourcePrimeNodePtr
        += (mi.gigs + mi.gIgtots)*mi.numberParallel;
    }
    else if (mi.rgateMod == 1)
    {

      *mi.f_GatePrimeEquGatePrimeNodePtr
        += (mi.gigg_jac + mi.Gmin + mi.gIgtotg + mi.geltd)*mi.numberParallel;

      *mi.f_GatePrimeEquDrainPrimeNodePtr
        += (mi.gigd - mi.Gmin + mi.gIgtotd)*mi.numberParallel;

      *mi.f_GatePrimeEquSourcePrimeNodePtr
        += (mi.gigs_jac + mi.gIgtots)*mi.numberParallel;
    }
    else
    {

      *mi.f_GatePrimeEquGatePrimeNodePtr
        += (mi.gigg + mi.Gmin + mi.gIgtotg - mi.gcrgg_jac)*mi.numberParallel;

      *mi.f_GatePrimeEquDrainPrimeNodePtr
        += (mi.gigd_jac - mi.Gmin + mi.gIgtotd - mi.gcrgd_jac)*mi.numberParallel;

      *mi.f_GatePrimeEquSourcePrimeNodePtr
        += (mi.gigs + mi.gIgtots - mi.gcrgs_jac)*mi.numberParallel;
    }

    *mi.f_DrainPrimeEquGatePrimeNodePtr
      += ((mi.Gm) + mi.gddpg - mi.Gmin - mi.gIdtotg)*mi.numberParallel;

    *mi.f_DrainPrimeEquDrainPrimeNodePtr
      += ((mi.drainConductance + mi.gds + mi.gddpdp + mi.RevSum) + mi.Gmin - mi.gIdtotd)*mi.numberParallel;

    *mi.f_DrainPrimeEquSourcePrimeNodePtr
      -= ((-mi.gddpsp + mi.gds + mi.FwdSum) + mi.gIdtots)*mi.numberParallel;

    *mi.f_DrainPrimeEquDrainNodePtr
      += -mi.drainConductance*mi.numberParallel;

    *mi.f_SourcePrimeEquGatePrimeNodePtr
      += (-mi.Gm + mi.gsspg - mi.gIstotg)*mi.numberParallel;

    *mi.f_SourcePrimeEquDrainPrimeNodePtr
      -= (mi.gds - mi.gsspdp + mi.RevSum + mi.gIstotd)*mi.numberParallel;

    *mi.f_SourcePrimeEquSourcePrimeNodePtr
      += ((mi.sourceConductance + mi.gds + mi.gsspsp + mi.FwdSum) + mi.Gmin - mi.gIstots)*mi.numberParallel;

    *mi.f_SourcePrimeEquSourceNodePtr
      -= mi.sourceConductance*mi.numberParallel;

    *mi.f_DrainEquDrainNodePtr
      += mi.drainConductance*mi.numberParallel;

    *mi.f_DrainEquDrainPrimeNodePtr
      -= mi.drainConductance*mi.numberParallel;

    *mi.f_SourceEquSourceNodePtr
      += mi.sourceConductance*mi.numberParallel;

    *mi.f_SourceEquSourcePrimeNodePtr
      -= mi.sourceConductance*mi.numberParallel;
    if (mi.bodyMod == 1)
    {

      *mi.f_BodyEquExtBodyNodePtr
        -= mi.gppp*mi.numberParallel;

      *mi.f_ExtBodyEquBodyNodePtr
        += mi.gppb*mi.numberParallel;

      *mi.f_ExtBodyEquExtBodyNodePtr
        += mi.gppp*mi.numberParallel;
    }
    if (mi.selfheat)
    {

      *mi.f_DrainPrimeEquTemperatureNodePtr
        += (mi.GmT + mi.gddpT)*mi.numberParallel;

      *mi.f_SourcePrimeEquTemperatureNodePtr
        += (-mi.GmT + mi.gsspT)*mi.numberParallel;


      *mi.f_GatePrimeEquTemperatureNodePtr
        += (mi.gigT_jac)*mi.numberParallel;

      *mi.f_TemperatureEquTemperatureNodePtr
        += (mi.gTtt  + 1/mi.paramPtr->rth)*mi.numberParallel;

      *mi.f_TemperatureEquGatePrimeNodePtr
        += mi.gTtg*mi.numberParallel;

      *mi.f_TemperatureEquDrainPrimeNodePtr
        += mi.gTtdp*mi.numberParallel;

      *mi.f_TemperatureEquSourcePrimeNodePtr
        += mi.gTtsp*mi.numberParallel;
      if (mi.soiMod != 0)
      {
        *mi.f_TemperatureEquSubstrateNodePtr
          += mi.gTte*mi.numberParallel;
      }
      if (mi.bNode > 0)
      {
        *mi.f_BodyEquTemperatureNodePtr
          += (mi.gbbT - mi.gigT_jac)*mi.numberParallel;
        *mi.f_TemperatureEquBodyNodePtr
          += mi.gTtb*mi.numberParallel;
      }
    }
    if( mi.icVDSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_DrainEquIdsPtr += 1.0;
        *mi.f_SourceEquIdsPtr -= 1.0;
        *mi.f_icVDSEquVdPtr += 1.0;
        *mi.f_icVDSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVDSEquIdsPtr += 1.0;
      }
    }

    if( mi.icVGSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_GateEquIgsPtr += 1.0;
        *mi.f_SourceEquIgsPtr -= 1.0;
        *mi.f_icVGSEquVgPtr += 1.0;
        *mi.f_icVGSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVGSEquIgsPtr += 1.0;
      }
    }

    if( mi.icVBSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_BodyEquIbsPtr += 1.0;
        *mi.f_SourceEquIbsPtr -= 1.0;
        *mi.f_icVBSEquVbPtr += 1.0;
        *mi.f_icVBSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVBSEquIbsPtr += 1.0;
      }
    }

    if( mi.icVESGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_SubstrateEquIesPtr += 1.0;
        *mi.f_SourceEquIesPtr -= 1.0;
        *mi.f_icVESEquVePtr += 1.0;
        *mi.f_icVESEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVESEquIesPtr += 1.0;
      }
    }

    if( mi.icVPSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_ExtBodyEquIpsPtr += 1.0;
        *mi.f_SourceEquIpsPtr -= 1.0;
        *mi.f_icVPSEquVpPtr += 1.0;
        *mi.f_icVPSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVPSEquIpsPtr += 1.0;
      }
    }


    // Q-matrix:

    if (mi.rgateMod == 3)
    {

      *mi.q_GateMidEquGateMidNodePtr
        += (mi.CAPcgmgmb)*mi.numberParallel;

      *mi.q_GateMidEquDrainPrimeNodePtr
        += (mi.CAPcgmdb)*mi.numberParallel;


      *mi.q_GateMidEquSourcePrimeNodePtr
        += (mi.CAPcgmsb)*mi.numberParallel;

      *mi.q_GateMidEquSubstrateNodePtr
        += mi.CAPcgmeb*mi.numberParallel;


      *mi.q_DrainPrimeEquGateMidNodePtr
        += mi.CAPcdgmb*mi.numberParallel;


      *mi.q_SourcePrimeEquGateMidNodePtr
        += mi.CAPcsgmb*mi.numberParallel;

      *mi.q_SubstrateEquGateMidNodePtr
        += mi.CAPcegmb*mi.numberParallel;
    }


    *mi.q_SubstrateEquDrainPrimeNodePtr
        += mi.CAPcedb*mi.numberParallel;


    *mi.q_DrainPrimeEquSubstrateNodePtr
        += mi.CAPcdeb*mi.numberParallel;

    *mi.q_SourcePrimeEquSubstrateNodePtr
        += mi.CAPcseb*mi.numberParallel;

    *mi.q_SubstrateEquGatePrimeNodePtr
        += mi.CAPcegb*mi.numberParallel;

    *mi.q_GatePrimeEquSubstrateNodePtr
        += mi.CAPcgeb*mi.numberParallel;
    if (mi.soiMod != 2)
    {

      *mi.q_SubstrateEquBodyNodePtr
        -= (mi.CAPcegb + mi.CAPcedb + mi.CAPcesb + mi.CAPceeb + mi.CAPcegmb)*mi.numberParallel;
      if (mi.rgateMod == 0 || mi.rgateMod == 1)
      {
        *mi.q_GatePrimeEquBodyNodePtr
        -= (+mi.CAPcggb + mi.CAPcgdb + mi.CAPcgsb + mi.CAPcgeb)*mi.numberParallel;
      }
      else
      {
        *mi.q_GatePrimeEquBodyNodePtr
        += (+ mi.CAPcgbb)*mi.numberParallel;
      }

      *mi.q_DrainPrimeEquBodyNodePtr
        -= ((+ mi.CAPcdgb + mi.CAPcddb + mi.CAPcdeb + mi.CAPcdsb) + mi.CAPcdgmb)*mi.numberParallel;

      *mi.q_SourcePrimeEquBodyNodePtr
        -= ((+ mi.CAPcsgb + mi.CAPcsdb + mi.CAPcseb + mi.CAPcssb) + mi.CAPcsgmb)*mi.numberParallel;

      *mi.q_BodyEquSubstrateNodePtr
        += (+ mi.CAPcbeb)*mi.numberParallel;

      *mi.q_BodyEquGatePrimeNodePtr
        += (+ mi.CAPcbgb)*mi.numberParallel;

      *mi.q_BodyEquDrainPrimeNodePtr
        += (+ mi.CAPcbdb)*mi.numberParallel;

      *mi.q_BodyEquSourcePrimeNodePtr
        += (mi.CAPcbsb)*mi.numberParallel;

      *mi.q_BodyEquBodyNodePtr
        += (-mi.CAPcbgb - mi.CAPcbdb - mi.CAPcbsb - mi.CAPcbeb)*mi.numberParallel;
    }

    *mi.q_SubstrateEquSubstrateNodePtr
      += mi.CAPceeb*mi.numberParallel;
    if (mi.rgateMod == 0)
    {

      *mi.q_GatePrimeEquGatePrimeNodePtr
        += (mi.CAPcggb)*mi.numberParallel;

      *mi.q_GatePrimeEquDrainPrimeNodePtr
        += (mi.CAPcgdb)*mi.numberParallel;

      *mi.q_GatePrimeEquSourcePrimeNodePtr
        += (mi.CAPcgsb)*mi.numberParallel;
    }
    else if (mi.rgateMod == 1)
    {

      *mi.q_GatePrimeEquGatePrimeNodePtr
        += (mi.CAPcggb)*mi.numberParallel;

      *mi.q_GatePrimeEquDrainPrimeNodePtr
        += (mi.CAPcgdb)*mi.numberParallel;

      *mi.q_GatePrimeEquSourcePrimeNodePtr
        += (mi.CAPcgsb)*mi.numberParallel;
    }
    else
    {

      *mi.q_GatePrimeEquGatePrimeNodePtr
        += (mi.CAPcggb)*mi.numberParallel;

      *mi.q_GatePrimeEquDrainPrimeNodePtr
        += (mi.CAPcgdb)*mi.numberParallel;

      *mi.q_GatePrimeEquSourcePrimeNodePtr
        += (mi.CAPcgsb)*mi.numberParallel;
    }

    *mi.q_DrainPrimeEquGatePrimeNodePtr
      += ((mi.CAPcdgb))*mi.numberParallel;

    *mi.q_DrainPrimeEquDrainPrimeNodePtr
      += ((mi.CAPcddb) )*mi.numberParallel;

    *mi.q_DrainPrimeEquSourcePrimeNodePtr
      -= ((-mi.CAPcdsb) )*mi.numberParallel;


    *mi.q_SourcePrimeEquGatePrimeNodePtr
      += (mi.CAPcsgb)*mi.numberParallel;

    *mi.q_SourcePrimeEquDrainPrimeNodePtr
      -= (-mi.CAPcsdb )*mi.numberParallel;

    *mi.q_SourcePrimeEquSourcePrimeNodePtr
      += ((mi.CAPcssb) )*mi.numberParallel;

    if (mi.selfheat)
    {

      *mi.q_DrainPrimeEquTemperatureNodePtr
        += (mi.CAPcdT)*mi.numberParallel;

      *mi.q_SourcePrimeEquTemperatureNodePtr
        += (mi.CAPcsT)*mi.numberParallel;

      *mi.q_SubstrateEquTemperatureNodePtr
        += mi.CAPceT*mi.numberParallel;

      *mi.q_GatePrimeEquTemperatureNodePtr
        += (mi.CAPcgT)*mi.numberParallel;

      *mi.q_TemperatureEquTemperatureNodePtr
        += (mi.CAPcTt)*mi.numberParallel;

      if (mi.bNode > 0)
      {
        *mi.q_BodyEquTemperatureNodePtr
        += (mi.CAPcbT)*mi.numberParallel;
      }
    }
  }
  return true;
}
#else
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/06/09
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  int sizeInstances = instanceContainer_.size();

#ifdef _OMP
#pragma omp parallel for
#endif
  for (int i=0; i<sizeInstances; ++i)
  {
    Instance & mi = *(instanceContainer_.at(i));

    // F-matrix:
    if (mi.rgateMod == 1)
    {

      dFdx[mi.li_Gate][mi.AGateEquGateNodeOffset]
        += mi.geltd*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquGateNodeOffset]
        -= mi.geltd*mi.numberParallel;

      dFdx[mi.li_Gate][mi.AGateEquGatePrimeNodeOffset]
        -= mi.geltd*mi.numberParallel;
      // It seems that this should be here, but it shouldn't.
      // The mi.geltd term is added into the G'-G' element later down in this
      // routine
      //
      //    dFdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset]
      //      += mi.geltd*mi.numberParallel;
    }
    else if (mi.rgateMod == 2)
    {

      dFdx[mi.li_Gate][mi.AGateEquGateNodeOffset]
        += mi.gcrg_jac*mi.numberParallel;

      dFdx[mi.li_Gate][mi.AGateEquGatePrimeNodeOffset]
        += mi.gcrgg_jac*mi.numberParallel;

      dFdx[mi.li_Gate][mi.AGateEquDrainPrimeNodeOffset]
        += mi.gcrgd_jac*mi.numberParallel;

      dFdx[mi.li_Gate][mi.AGateEquSourcePrimeNodeOffset]
        += mi.gcrgs_jac*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquGateNodeOffset]
        -= mi.gcrg_jac*mi.numberParallel;
      if (mi.soiMod != 2)
      {
        dFdx[mi.li_Gate][mi.AGateEquBodyNodeOffset]
          += mi.gcrgb_jac*mi.numberParallel;
      }
    }
    else if (mi.rgateMod == 3)
    {

      dFdx[mi.li_Gate][mi.AGateEquGateNodeOffset]
        += mi.geltd*mi.numberParallel;

      dFdx[mi.li_Gate][mi.AGateEquGateMidNodeOffset]
        += -mi.geltd*mi.numberParallel;

      dFdx[mi.li_GateMid][mi.AGateMidEquGateNodeOffset]
        += -mi.geltd*mi.numberParallel;

      dFdx[mi.li_GateMid][mi.AGateMidEquGateMidNodeOffset]
        += (mi.geltd + mi.gcrg_jac)*mi.numberParallel;

      dFdx[mi.li_GateMid][mi.AGateMidEquDrainPrimeNodeOffset]
        += (mi.gcrgd_jac)*mi.numberParallel;

      dFdx[mi.li_GateMid][mi.AGateMidEquGatePrimeNodeOffset]
        += mi.gcrgg_jac*mi.numberParallel;

      dFdx[mi.li_GateMid][mi.AGateMidEquSourcePrimeNodeOffset]
        += (mi.gcrgs_jac)*mi.numberParallel;

      if (mi.soiMod != 2)
      {
        dFdx[mi.li_GateMid][mi.AGateMidEquBodyNodeOffset]
          += mi.gcrgb_jac*mi.numberParallel;
      }


      dFdx[mi.li_GatePrime][mi.AGatePrimeEquGateMidNodeOffset]
        -= mi.gcrg_jac*mi.numberParallel;
    }
    if (mi.soiMod != 0)
    {

      dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquSubstrateNodeOffset]
        += (mi.Gme + mi.gddpe)*mi.numberParallel;

      dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSubstrateNodeOffset]
        += (mi.gsspe - mi.Gme)*mi.numberParallel;
      if (mi.soiMod != 2)
      {
        dFdx[mi.li_GatePrime][mi.AGatePrimeEquSubstrateNodeOffset]
          += mi.gige_jac*mi.numberParallel;
        dFdx[mi.li_Body][mi.ABodyEquSubstrateNodeOffset]
          -= mi.gige_jac*mi.numberParallel;
      }
    }

    if (mi.soiMod != 2)
    {
      if (mi.rgateMod == 0 || mi.rgateMod == 1)
      {
        dFdx[mi.li_GatePrime][mi.AGatePrimeEquBodyNodeOffset]
          -= (-mi.gigb_jac - mi.gIgtotb)*mi.numberParallel;
      }
      else
      {
        dFdx[mi.li_GatePrime][mi.AGatePrimeEquBodyNodeOffset]
          += (mi.gigb_jac +mi.gIgtotb - mi.gcrgb_jac)*mi.numberParallel;
      }

      dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquBodyNodeOffset]
        -= ((-mi.gddpb - mi.Gmbs) + mi.gIdtotb)*mi.numberParallel;

      dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquBodyNodeOffset]
        -= ((-mi.gsspb + mi.Gmbs) + mi.Gmin + mi.gIstotb)*mi.numberParallel;

      dFdx[mi.li_Body][mi.ABodyEquSubstrateNodeOffset]
        += (mi.gbbe)*mi.numberParallel;

      dFdx[mi.li_Body][mi.ABodyEquGatePrimeNodeOffset]
        += (-mi.gigg + mi.gbbg)*mi.numberParallel;

      dFdx[mi.li_Body][mi.ABodyEquDrainPrimeNodeOffset]
        += (-mi.gigd_jac + mi.gbbdp)*mi.numberParallel;

      dFdx[mi.li_Body][mi.ABodyEquSourcePrimeNodeOffset]
        += (mi.gbbsp - mi.Gmin - mi.gigs)*mi.numberParallel;

      dFdx[mi.li_Body][mi.ABodyEquBodyNodeOffset]
        += (-mi.gigb_jac + mi.gbbb + mi.Gmin)*mi.numberParallel;
    }
    if (mi.rgateMod == 0)
    {

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset]
        += (mi.gigg + mi.Gmin + mi.gIgtotg)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquDrainPrimeNodeOffset]
        += (mi.gigd_jac - mi.Gmin + mi.gIgtotd)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquSourcePrimeNodeOffset]
        += (mi.gigs + mi.gIgtots)*mi.numberParallel;
    }
    else if (mi.rgateMod == 1)
    {

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset]
        += (mi.gigg_jac + mi.Gmin + mi.gIgtotg + mi.geltd)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquDrainPrimeNodeOffset]
        += (mi.gigd - mi.Gmin + mi.gIgtotd)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquSourcePrimeNodeOffset]
        += (mi.gigs_jac + mi.gIgtots)*mi.numberParallel;
    }
    else
    {

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset]
        += (mi.gigg + mi.Gmin + mi.gIgtotg - mi.gcrgg_jac)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquDrainPrimeNodeOffset]
        += (mi.gigd_jac - mi.Gmin + mi.gIgtotd - mi.gcrgd_jac)*mi.numberParallel;

      dFdx[mi.li_GatePrime][mi.AGatePrimeEquSourcePrimeNodeOffset]
        += (mi.gigs + mi.gIgtots - mi.gcrgs_jac)*mi.numberParallel;
    }

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquGatePrimeNodeOffset]
      += ((mi.Gm) + mi.gddpg - mi.Gmin - mi.gIdtotg)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset]
      += ((mi.drainConductance + mi.gds + mi.gddpdp + mi.RevSum) + mi.Gmin - mi.gIdtotd)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset]
      -= ((-mi.gddpsp + mi.gds + mi.FwdSum) + mi.gIdtots)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainNodeOffset]
      += -mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquGatePrimeNodeOffset]
      += (-mi.Gm + mi.gsspg - mi.gIstotg)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset]
      -= (mi.gds - mi.gsspdp + mi.RevSum + mi.gIstotd)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset]
      += ((mi.sourceConductance + mi.gds + mi.gsspsp + mi.FwdSum) + mi.Gmin - mi.gIstots)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourceNodeOffset]
      -= mi.sourceConductance*mi.numberParallel;

    dFdx[mi.li_Drain][mi.ADrainEquDrainNodeOffset]
      += mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_Drain][mi.ADrainEquDrainPrimeNodeOffset]
      -= mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_Source][mi.ASourceEquSourceNodeOffset]
      += mi.sourceConductance*mi.numberParallel;

    dFdx[mi.li_Source][mi.ASourceEquSourcePrimeNodeOffset]
      -= mi.sourceConductance*mi.numberParallel;
    if (mi.bodyMod == 1)
    {

      dFdx[mi.li_Body][mi.ABodyEquExtBodyNodeOffset]
        -= mi.gppp*mi.numberParallel;

      dFdx[mi.li_ExtBody][mi.AExtBodyEquBodyNodeOffset]
        += mi.gppb*mi.numberParallel;

      dFdx[mi.li_ExtBody][mi.AExtBodyEquExtBodyNodeOffset]
        += mi.gppp*mi.numberParallel;
    }
    if (mi.selfheat)
    {

      dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquTemperatureNodeOffset]
        += (mi.GmT + mi.gddpT)*mi.numberParallel;

      dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquTemperatureNodeOffset]
        += (-mi.GmT + mi.gsspT)*mi.numberParallel;


      dFdx[mi.li_GatePrime][mi.AGatePrimeEquTemperatureNodeOffset]
        += (mi.gigT_jac)*mi.numberParallel;

      dFdx[mi.li_Temperature][mi.ATemperatureEquTemperatureNodeOffset]
        += (mi.gTtt  + 1/mi.paramPtr->rth)*mi.numberParallel;

      dFdx[mi.li_Temperature][mi.ATemperatureEquGatePrimeNodeOffset]
        += mi.gTtg*mi.numberParallel;

      dFdx[mi.li_Temperature][mi.ATemperatureEquDrainPrimeNodeOffset]
        += mi.gTtdp*mi.numberParallel;

      dFdx[mi.li_Temperature][mi.ATemperatureEquSourcePrimeNodeOffset]
        += mi.gTtsp*mi.numberParallel;
      if (mi.soiMod != 0)
      {
        dFdx[mi.li_Temperature][mi.ATemperatureEquSubstrateNodeOffset]
          += mi.gTte*mi.numberParallel;
      }
      if (mi.bNode > 0)
      {
        dFdx[mi.li_Body][mi.ABodyEquTemperatureNodeOffset]
          += (mi.gbbT - mi.gigT_jac)*mi.numberParallel;
        dFdx[mi.li_Temperature][mi.ATemperatureEquBodyNodeOffset]
          += mi.gTtb*mi.numberParallel;
      }
    }
    if( mi.icVDSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Drain][mi.ADrainEquIdsOffset] += 1.0;
        dFdx[mi.li_Source][mi.ASourceEquIdsOffset] -= 1.0;
        dFdx[mi.li_Ids][mi.icVDSEquVdOffset] += 1.0;
        dFdx[mi.li_Ids][mi.icVDSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Ids][mi.icVDSEquIdsOffset] += 1.0;
      }
    }

    if( mi.icVGSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Gate][mi.AGateEquIgsOffset] += 1.0;
        dFdx[mi.li_Source][mi.ASourceEquIgsOffset] -= 1.0;
        dFdx[mi.li_Igs][mi.icVGSEquVgOffset] += 1.0;
        dFdx[mi.li_Igs][mi.icVGSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Igs][mi.icVGSEquIgsOffset] += 1.0;
      }
    }

    if( mi.icVBSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Body][mi.ABodyEquIbsOffset] += 1.0;
        dFdx[mi.li_Source][mi.ASourceEquIbsOffset] -= 1.0;
        dFdx[mi.li_Ibs][mi.icVBSEquVbOffset] += 1.0;
        dFdx[mi.li_Ibs][mi.icVBSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Ibs][mi.icVBSEquIbsOffset] += 1.0;
      }
    }

    if( mi.icVESGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Substrate][mi.ASubstrateEquIesOffset] += 1.0;
        dFdx[mi.li_Source][mi.ASourceEquIesOffset] -= 1.0;
        dFdx[mi.li_Ies][mi.icVESEquVeOffset] += 1.0;
        dFdx[mi.li_Ies][mi.icVESEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Ies][mi.icVESEquIesOffset] += 1.0;
      }
    }

    if( mi.icVPSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_ExtBody][mi.AExtBodyEquIpsOffset] += 1.0;
        dFdx[mi.li_Source][mi.ASourceEquIpsOffset] -= 1.0;
        dFdx[mi.li_Ips][mi.icVPSEquVpOffset] += 1.0;
        dFdx[mi.li_Ips][mi.icVPSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Ips][mi.icVPSEquIpsOffset] += 1.0;
      }
    }


    // Q-matrix:

    if (mi.rgateMod == 3)
    {

      dQdx[mi.li_GateMid][mi.AGateMidEquGateMidNodeOffset]
        += (mi.CAPcgmgmb)*mi.numberParallel;

      dQdx[mi.li_GateMid][mi.AGateMidEquDrainPrimeNodeOffset]
        += (mi.CAPcgmdb)*mi.numberParallel;


      dQdx[mi.li_GateMid][mi.AGateMidEquSourcePrimeNodeOffset]
        += (mi.CAPcgmsb)*mi.numberParallel;

      dQdx[mi.li_GateMid][mi.AGateMidEquSubstrateNodeOffset]
        += mi.CAPcgmeb*mi.numberParallel;


      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquGateMidNodeOffset]
        += mi.CAPcdgmb*mi.numberParallel;


      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquGateMidNodeOffset]
        += mi.CAPcsgmb*mi.numberParallel;

      dQdx[mi.li_Substrate][mi.ASubstrateEquGateMidNodeOffset]
        += mi.CAPcegmb*mi.numberParallel;
    }


    dQdx[mi.li_Substrate][mi.ASubstrateEquDrainPrimeNodeOffset]
        += mi.CAPcedb*mi.numberParallel;


    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquSubstrateNodeOffset]
        += mi.CAPcdeb*mi.numberParallel;

    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquSubstrateNodeOffset]
        += mi.CAPcseb*mi.numberParallel;

    dQdx[mi.li_Substrate][mi.ASubstrateEquGatePrimeNodeOffset]
        += mi.CAPcegb*mi.numberParallel;

    dQdx[mi.li_GatePrime][mi.AGatePrimeEquSubstrateNodeOffset]
        += mi.CAPcgeb*mi.numberParallel;
    if (mi.soiMod != 2)
    {

      dQdx[mi.li_Substrate][mi.ASubstrateEquBodyNodeOffset]
        -= (mi.CAPcegb + mi.CAPcedb + mi.CAPcesb + mi.CAPceeb + mi.CAPcegmb)*mi.numberParallel;
      if (mi.rgateMod == 0 || mi.rgateMod == 1)
      {
        dQdx[mi.li_GatePrime][mi.AGatePrimeEquBodyNodeOffset]
          -= (+mi.CAPcggb + mi.CAPcgdb + mi.CAPcgsb + mi.CAPcgeb)*mi.numberParallel;
      }
      else
      {
        dQdx[mi.li_GatePrime][mi.AGatePrimeEquBodyNodeOffset]
          += (+ mi.CAPcgbb)*mi.numberParallel;
      }

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquBodyNodeOffset]
        -= ((+ mi.CAPcdgb + mi.CAPcddb + mi.CAPcdeb + mi.CAPcdsb) + mi.CAPcdgmb)*mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquBodyNodeOffset]
        -= ((+ mi.CAPcsgb + mi.CAPcsdb + mi.CAPcseb + mi.CAPcssb) + mi.CAPcsgmb)*mi.numberParallel;

      dQdx[mi.li_Body][mi.ABodyEquSubstrateNodeOffset]
        += (+ mi.CAPcbeb)*mi.numberParallel;

      dQdx[mi.li_Body][mi.ABodyEquGatePrimeNodeOffset]
        += (+ mi.CAPcbgb)*mi.numberParallel;

      dQdx[mi.li_Body][mi.ABodyEquDrainPrimeNodeOffset]
        += (+ mi.CAPcbdb)*mi.numberParallel;

      dQdx[mi.li_Body][mi.ABodyEquSourcePrimeNodeOffset]
        += (mi.CAPcbsb)*mi.numberParallel;

      dQdx[mi.li_Body][mi.ABodyEquBodyNodeOffset]
        += (-mi.CAPcbgb - mi.CAPcbdb - mi.CAPcbsb - mi.CAPcbeb)*mi.numberParallel;
    }

    dQdx[mi.li_Substrate][mi.ASubstrateEquSubstrateNodeOffset]
      += mi.CAPceeb*mi.numberParallel;
    if (mi.rgateMod == 0)
    {

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset]
        += (mi.CAPcggb)*mi.numberParallel;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquDrainPrimeNodeOffset]
        += (mi.CAPcgdb)*mi.numberParallel;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquSourcePrimeNodeOffset]
        += (mi.CAPcgsb)*mi.numberParallel;
    }
    else if (mi.rgateMod == 1)
    {

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset]
        += (mi.CAPcggb)*mi.numberParallel;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquDrainPrimeNodeOffset]
        += (mi.CAPcgdb)*mi.numberParallel;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquSourcePrimeNodeOffset]
        += (mi.CAPcgsb)*mi.numberParallel;
    }
    else
    {

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquGatePrimeNodeOffset]
        += (mi.CAPcggb)*mi.numberParallel;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquDrainPrimeNodeOffset]
        += (mi.CAPcgdb)*mi.numberParallel;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquSourcePrimeNodeOffset]
        += (mi.CAPcgsb)*mi.numberParallel;
    }

    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquGatePrimeNodeOffset]
      += ((mi.CAPcdgb))*mi.numberParallel;

    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset]
      += ((mi.CAPcddb) )*mi.numberParallel;

    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset]
      -= ((-mi.CAPcdsb) )*mi.numberParallel;


    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquGatePrimeNodeOffset]
      += (mi.CAPcsgb)*mi.numberParallel;

    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset]
      -= (-mi.CAPcsdb )*mi.numberParallel;

    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset]
      += ((mi.CAPcssb) )*mi.numberParallel;

    if (mi.selfheat)
    {

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquTemperatureNodeOffset]
        += (mi.CAPcdT)*mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquTemperatureNodeOffset]
        += (mi.CAPcsT)*mi.numberParallel;

      dQdx[mi.li_Substrate][mi.ASubstrateEquTemperatureNodeOffset]
        += mi.CAPceT*mi.numberParallel;

      dQdx[mi.li_GatePrime][mi.AGatePrimeEquTemperatureNodeOffset]
        += (mi.CAPcgT)*mi.numberParallel;

      dQdx[mi.li_Temperature][mi.ATemperatureEquTemperatureNodeOffset]
        += (mi.CAPcTt)*mi.numberParallel;

      if (mi.bNode > 0)
      {
        dQdx[mi.li_Body][mi.ABodyEquTemperatureNodeOffset]
          += (mi.CAPcbT)*mi.numberParallel;
      }
    }
  }
  return true;
}
#endif

} // namespace Resistor
} // namespace Device
} // namespace Xyce
