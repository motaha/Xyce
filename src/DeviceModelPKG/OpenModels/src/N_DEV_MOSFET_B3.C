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
// Filename       : $RCSfile: N_DEV_MOSFET_B3.C,v $
//
// Purpose        : This file implements the BSIM3 MOSFET model.  It
//                  is intended to be compatible with the Berkeley SPICE
//                  (3f5) version, BSIM3 version 3.2.2.
//
// Special Notes  :
//
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/12/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.346.2.3 $
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

// ----------   Xyce Includes   ----------

#ifndef Xyce_USE_BSIM3_CONST
#include <N_DEV_Const.h>
#endif

#include <N_DEV_MOSFET_B3.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

// ---------- BSIM3 constants ---------------
// Normally, these should be obtained from N_DEV_Const.h
// I have them here, for purposes of comparison to spice.
#ifdef Xyce_USE_BSIM3_CONST

#define CONSTMAX_EXP 5.834617425e+14
#define CONSTMIN_EXP 1.713908431e-15

#define CONSTEXP_THRESHOLD 34.0
#define CONSTEPSOX 3.453133e-11
#define CONSTEPSSI 1.03594e-10

#define CONSTQ 1.60219e-19

#define CONSTDELTA_1 0.02
#define CONSTDELTA_2 0.02
#define CONSTDELTA_3 0.02
#define CONSTDELTA_4 0.02

#define CONSTroot2    sqrt(2.0)
#define CONSTCtoK     (273.15)
#define CONSTREFTEMP  (300.15)

#define CONSTboltz    1.3806226e-23
#define CONSTKoverQ   8.617087e-5  // Kb / q  where q = 1.60219e-19

#define CONSTvt0     (CONSTboltz * (27.0 +CONSTCtoK)/CONSTQ)

#define CONSTEg300   (1.1150877) // band gap for Si at T=300.15K (room temp)
#define CONSTEg0     (1.16)      // band gap for Si at T=0K. (eV)
#define CONSTalphaEg (7.02e-4)   // (eV/K)
#define CONSTbetaEg  (1108.0)    // (K)

#define CONSTNi0     (1.45e10)   // carrier concentration at room temp.

#define CONSTNMOS 1
#define CONSTPMOS -1

#endif


namespace Xyce {
namespace Device {

template<>
ParametricData<MOSFET_B3::Instance>::ParametricData()
{
    // Set up configuration constants:
    setNumNodes(4);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(1);
    addModelType("NMOS");
    addModelType("PMOS");

    // Set up double precision variables:
    addPar ("TEMP", 0.0, false, TIME_DEP,
      &MOSFET_B3::Instance::temp,
      NULL,  STANDARD, CAT_NONE, "");

    addPar ("L", 0.0, true, NO_DEP,
      &MOSFET_B3::Instance::l,
      NULL, U_METER,  CAT_GEOMETRY, "Channel length");

    addPar ("W", 0.0, true, NO_DEP,
      &MOSFET_B3::Instance::w,
      NULL, U_METER,  CAT_GEOMETRY, "Channel width");

    addPar ("AD", 0.0, false, NO_DEP,
      &MOSFET_B3::Instance::drainArea,
      NULL, U_METER2,  CAT_GEOMETRY, "Drain diffusion area");

    addPar ("AS", 0.0, false, NO_DEP,
      &MOSFET_B3::Instance::sourceArea,
      NULL, U_METER2,  CAT_GEOMETRY, "Source diffusion area");

    addPar ("NRD", 1.0, false, NO_DEP,
      &MOSFET_B3::Instance::drainSquares,
      NULL, U_SQUARES, CAT_GEOMETRY, "Multiplier for RSH to yield parasitic resistance of drain");

    addPar ("NRS", 1.0, false, NO_DEP,
      &MOSFET_B3::Instance::sourceSquares,
      NULL, U_SQUARES, CAT_GEOMETRY, "Multiplier for RSH to yield parasitic resistance of source");

    addPar ("PD", 0.0, false, NO_DEP,
      &MOSFET_B3::Instance::drainPerimeter,
      NULL, U_METER, CAT_GEOMETRY, "Drain diffusion perimeter");

    addPar ("PS", 0.0, false, NO_DEP,
      &MOSFET_B3::Instance::sourcePerimeter,
      NULL, U_METER, CAT_GEOMETRY, "Source diffusion perimeter");

    addPar ("M", 1.0, false, NO_DEP,
      &MOSFET_B3::Instance::numberParallel,
      NULL, U_NONE, CAT_CONTROL,  "Multiplier for M devices connected in parallel");

    addPar ("IC1", 0.0, false, NO_DEP,
      &MOSFET_B3::Instance::icVDS,
      &MOSFET_B3::Instance::icVDSGiven,
      U_VOLT, CAT_VOLT, "Initial condition on Vds");

    addPar ("IC2", 0.0, false, NO_DEP,
      &MOSFET_B3::Instance::icVGS,
      &MOSFET_B3::Instance::icVGSGiven,
      U_VOLT, CAT_VOLT, "Initial condition on Vgs");

    addPar ("IC3", 0.0, false, NO_DEP,
      &MOSFET_B3::Instance::icVBS,
      &MOSFET_B3::Instance::icVBSGiven,
      U_VOLT, CAT_VOLT, "Initial condition on Vbs");

    // Set up non-double precision variables:
    addPar ("NQSMOD", 0, false, NO_DEP,
            &MOSFET_B3::Instance::nqsMod, NULL,
            U_NONE, CAT_CONTROL, "Flag for NQS model");
    addPar ("OFF",false,false, NO_DEP,
            &MOSFET_B3::Instance::OFF,
            NULL, U_LOGIC, CAT_VOLT,
            "Initial condition of no voltage drops accross device");

    // This tells the parser that IC1, IC2, and IC3 are to be input as a vector of "IC"
    makeVector ("IC", 3);
}

template<>
ParametricData<MOSFET_B3::Model>::ParametricData()
{
    // Set up map for normal (double) param variables:
    addPar ("TOX", 150.e-10, true, NO_DEP,
      &MOSFET_B3::Model::tox,
      NULL, 
       U_METER, CAT_GEOMETRY, "Gate oxide thickness");

    addPar ("TOXM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::toxm,
      NULL, U_METER, CAT_PROCESS, "Gate oxide thickness used in extraction");

    addPar ("CDSC", 2.4e-4, false, NO_DEP,
      &MOSFET_B3::Model::cdsc,
      NULL, U_FARADMM2, CAT_DC, "Drain/source to channel coupling capacitance");

    addPar ("CDSCB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cdscb,
      NULL, U_FVM1MM2, CAT_DC, "Body-bias sensitivity of CDSC");

    addPar ("CDSCD", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cdscd,
      NULL, U_FVM1MM2, CAT_DC, "Drain-bias sensitivity of CDSC");

    addPar ("CIT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cit,
      NULL, U_FARADMM2, CAT_DC, "Interface trap capacitance");

    addPar ("NFACTOR", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::nfactor,
      NULL, U_NONE, CAT_DC, "Subthreshold swing factor");

    addPar ("XJ", 0.15e-6, false, NO_DEP,
      &MOSFET_B3::Model::xj,
      NULL, U_METER, CAT_GEOMETRY, "Junction depth");

    addPar ("VSAT", 8.0e4, false, NO_DEP,
      &MOSFET_B3::Model::vsat,
      NULL, U_MSM1, CAT_DC, "Saturation velocity at temp = TNOM");

    addPar ("AT", 3.3e4, false, NO_DEP,
      &MOSFET_B3::Model::at,
      NULL, U_MSM1, CAT_TEMP, "Temperature coefficient for saturation velocity");

    addPar ("A0", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::a0,
      NULL, U_NONE, CAT_DC, "Bulk charge effect coefficient for channel length");

    addPar ("AGS", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ags,
      NULL, U_VOLTM1, CAT_DC, "Gate-bias coefficient of abulk");

    addPar ("A1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::a1,
      NULL, U_VOLTM1, CAT_DC, "First non-saturation effect parameter");

    addPar ("A2", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::a2,
      NULL, U_NONE, CAT_DC, "Second non-saturation factor");

    addPar ("KETA", -0.047, false, NO_DEP,
      &MOSFET_B3::Model::keta,
      NULL, U_VOLTM1, CAT_DC, "Body-bias coefficient of bulk charge effect");

    addPar ("NSUB", 6.0e16, false, NO_DEP,
      &MOSFET_B3::Model::nsub,
      &MOSFET_B3::Model::nsubGiven,
       
       U_CMM3, CAT_DOPING, "Substrate doping density");

    addPar ("NCH", 1.7e17, false, NO_DEP,
      &MOSFET_B3::Model::npeak,
      &MOSFET_B3::Model::npeakGiven,
       U_CMM3, CAT_PROCESS, "Channel doping concentration");

    addPar ("NGATE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ngate,
      NULL, U_CMM3, CAT_DC, "Poly gate doping concentration");

    addPar ("GAMMA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::gamma1,
      &MOSFET_B3::Model::gamma1Given,
       U_VOLTH, CAT_PROCESS, "Body effect coefficient near the surface");

    addPar ("GAMMA2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::gamma2,
      &MOSFET_B3::Model::gamma2Given,
       U_VOLTH, CAT_PROCESS, "Body effect coefficient in the bulk");

    addPar ("VBX", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::vbx,
      &MOSFET_B3::Model::vbxGiven,
       U_VOLT, CAT_PROCESS, "Vbs at which the depetion region = XT");

    addPar ("VBM", -3.0, false, NO_DEP,
      &MOSFET_B3::Model::vbm,
      &MOSFET_B3::Model::vbmGiven,
       U_VOLT, CAT_DC, "Maximum applied body-bias in threshold voltage calculation");

    addPar ("XT", 1.55e-7, false, NO_DEP,
      &MOSFET_B3::Model::xt,
      &MOSFET_B3::Model::xtGiven,
       U_METER, CAT_PROCESS, "Doping depth");

    addPar ("K1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::k1,
      &MOSFET_B3::Model::k1Given,
       U_VOLTH, CAT_DC, "First-order body effect coefficient");

    addPar ("KT1", -0.11, false, NO_DEP,
      &MOSFET_B3::Model::kt1,
      NULL, U_VOLT, CAT_TEMP, "Themperature coefficient for threshold voltage");

    addPar ("KT1L", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::kt1l,
      NULL, U_VM, CAT_TEMP, "Channel length dependence of the temerature coefficient for the threshold voltage");

    addPar ("KT2", 0.022, false, NO_DEP,
      &MOSFET_B3::Model::kt2,
      NULL, U_NONE, CAT_TEMP, "Body-bias coefficient fo the threshold voltage temperature effect");

    addPar ("K2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::k2,
      &MOSFET_B3::Model::k2Given,
       U_NONE, CAT_DC, "second-order body effect coefficient");

    addPar ("K3", 80.0, false, NO_DEP,
      &MOSFET_B3::Model::k3,
      NULL, U_NONE, CAT_DC, "Narrow width coefficient");

    addPar ("K3B", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::k3b,
      NULL, U_VOLTM1, CAT_DC, "Body effect coefficient of K3");

    addPar ("W0", 2.5e-6, false, NO_DEP,
      &MOSFET_B3::Model::w0,
      NULL, U_METER, CAT_DC, "Narrow-width paameter");

    addPar ("NLX", 1.74e-7, false, NO_DEP,
      &MOSFET_B3::Model::nlx,
      NULL, U_METER, CAT_DC, "Lateral non-uniform doping parameter");

    addPar ("DVT0", 2.2, false, NO_DEP,
      &MOSFET_B3::Model::dvt0,
      NULL, U_NONE, CAT_DC, "First coefficient of short-channel effect effect on threshold voltage");

    addPar ("DVT1", 0.53, false, NO_DEP,
      &MOSFET_B3::Model::dvt1,
      NULL, 	 U_NONE, CAT_DC, "Second coefficient of short-channel effect effect on threshold voltage");

    addPar ("DVT2", -0.032, false, NO_DEP,
      &MOSFET_B3::Model::dvt2,
      NULL, 	 U_VOLTM1, CAT_DC, "Body-bias coefficient of short-channel effect effect on threshold voltage");

    addPar ("DVT0W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::dvt0w,
      NULL, 	 U_METERM1, CAT_DC, "First coefficient of narrow-width effect effect on threshold voltage for small channel length");

    addPar ("DVT1W", 5.3e6, false, NO_DEP,
      &MOSFET_B3::Model::dvt1w,
      NULL, 	 U_METERM1, CAT_DC, "Second coefficient of narrow-width effect effect on threshold voltage for small channel length");

    addPar ("DVT2W", -0.032, false, NO_DEP,
      &MOSFET_B3::Model::dvt2w,
      NULL, 	 U_VOLTM1, CAT_DC, "Body-bias coefficient of narrow-width effect effect on threshold voltage for small channel length");

    addPar ("DROUT", 0.56, false, NO_DEP,
      &MOSFET_B3::Model::drout,
      NULL, U_NONE, CAT_DC, "L-depedance Coefficient of the DIBL correction parameter in Rout");

    addPar ("DSUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::dsub,
      NULL, U_NONE, CAT_DC, "DIBL coefficient exponent in subthreshhold region");

    addPar ("VTH0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::vth0,
      &MOSFET_B3::Model::vth0Given,
       U_VOLT, CAT_DC, "Threshold voltage at Vbs = 0 for large L");

    addPar ("UA", 2.25e-9, false, NO_DEP,
      &MOSFET_B3::Model::ua,
      NULL, U_MVM1, CAT_DC, "First-order mobility degradation coefficient");

    addPar ("UA1", 4.31e-9, false, NO_DEP,
      &MOSFET_B3::Model::ua1,
      NULL, U_MVM1, CAT_TEMP, "Temperature coefficient for UA");

    addPar ("UB", 5.87e-19, false, NO_DEP,
      &MOSFET_B3::Model::ub,
      NULL, U_M2VM2, CAT_DC, "First-order mobility degradation coefficient");

    addPar ("UB1", -7.61e-18, false, NO_DEP,
      &MOSFET_B3::Model::ub1,
      NULL, U_M2VM2, CAT_TEMP, "Temperature coefficient for UB");

    addPar ("UC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::uc,
      NULL, U_MVM2, CAT_DC, "Body effect of mobility degridation coefficient");

    addPar ("UC1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::uc1,
      NULL, U_MVM2DEGCM1, CAT_TEMP, "Temperature coefficient for UC");

    addPar ("U0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::u0,
      NULL, 
       U_CMM2VM1SM1, CAT_PROCESS, "Surface mobility");

    addPar ("UTE", -1.5, false, NO_DEP,
      &MOSFET_B3::Model::ute,
      NULL, U_NONE, CAT_TEMP, "Mobility temerature exponent");

    addPar ("VOFF", -0.08, false, NO_DEP,
      &MOSFET_B3::Model::voff,
      NULL, U_VOLT, CAT_DC, "Offset voltage in the subthreshold region at large W and L");

    addPar ("RDSW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::rdsw,
      NULL, U_OHMMICRON, CAT_DC, "Parasitic resistance per unit width");

    addPar ("PRWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::prwg,
      NULL, U_VOLTM1, CAT_DC, "Gate-bias effect coefficient of RDSW");

    addPar ("PRWB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::prwb,
      NULL, U_VOLTMH, CAT_DC, "Body effect coefficient of RDSW");

    addPar ("PRT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::prt,
      NULL, U_OHMMICRON, CAT_TEMP, "Temerature coefficient for RDSW");

    addPar ("ETA0", 0.08, false, NO_DEP,
      &MOSFET_B3::Model::eta0,
      NULL, U_NONE, CAT_DC, "DIBL coefficient in subthreshold region");

    addPar ("ETAB", -0.07, false, NO_DEP,
      &MOSFET_B3::Model::etab,
      NULL, U_VOLTM1, CAT_DC, "Body-bias coefficient for the subthreshold DIBL effect");

    addPar ("PCLM", 1.3, false, NO_DEP,
      &MOSFET_B3::Model::pclm,
      NULL, U_NONE, CAT_DC, "Channel length modulation parameter");

    addPar ("PDIBLC1", 0.39, false, NO_DEP,
      &MOSFET_B3::Model::pdibl1,
      NULL, U_NONE, CAT_DC, "First output resistance DIBL effect correction parameter");

    addPar ("PDIBLC2", 0.0086, false, NO_DEP,
      &MOSFET_B3::Model::pdibl2,
      NULL, U_NONE, CAT_DC, "Second output resistance DIBL effect correction parameter");

    addPar ("PDIBLCB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdiblb,
      NULL, U_VOLTM1, CAT_DC, "Body effect coefficient of DIBL correction parameter");

    addPar ("PSCBE1", 4.24e8, false, NO_DEP,
      &MOSFET_B3::Model::pscbe1,
      NULL, U_VMM1, CAT_DC, "First substrate current body effect parameter");

    addPar ("PSCBE2", 1.0e-5, false, NO_DEP,
      &MOSFET_B3::Model::pscbe2,
      NULL, U_VMM1, CAT_DC, "second substrate current body effect parameter");

    addPar ("PVAG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvag,
      NULL, U_NONE, CAT_DC, "Gate dependence of early voltage");

    addPar ("DELTA", 0.01, false, NO_DEP,
      &MOSFET_B3::Model::delta,
      NULL, U_VOLT, CAT_DC, "Effective Vds parameter");

    addPar ("WR", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::wr,
      NULL, U_NONE, CAT_DC, "Width offset from Weff for Rds Calculation");

    addPar ("DWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::dwg,
      NULL, U_MVMH, CAT_DC, "Coefficient of gate depedence of Weff");

    addPar ("DWB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::dwb,
      NULL, U_MVMH, CAT_DC, "Coefficient of substrate body bias dependence of Weff");

    addPar ("B0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::b0,
      NULL, U_METER, CAT_DC, "Bulk charge effect coefficient for channel width");

    addPar ("B1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::b1,
      NULL, U_METER, CAT_DC, "Bulk charge effect offset");

    addPar ("ALPHA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::alpha0,
      NULL, U_MVM1, CAT_DC, "First parameter of impact-ionization current");

    addPar ("ALPHA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::alpha1,
      NULL, U_VOLTM1, CAT_DC, "Isub parameter for length scaling");

    addPar ("BETA0", 30.0, false, NO_DEP,
      &MOSFET_B3::Model::beta0,
      NULL, U_VOLT, CAT_DC, "Second parameter of impact-ionization current");

    addPar ("IJTH", 0.1, false, NO_DEP,
      &MOSFET_B3::Model::ijth,
      NULL, U_AMP, CAT_DC, "Diode limiting current");

    addPar ("VFB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::vfb,
      &MOSFET_B3::Model::vfbGiven,
       U_VOLT, CAT_DC, "Flat-band voltage");

    addPar ("ELM", 5.0, false, NO_DEP,
      &MOSFET_B3::Model::elm,
      NULL, U_NONE, CAT_NQS, "Elmore constant of the channel");

    addPar ("CGSL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cgsl,
      NULL, U_FARADMM1, CAT_CAP, "Light-doped source-gate region overlap capacitance");

    addPar ("CGDL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cgdl,
      NULL, U_FARADMM1, CAT_CAP, "Light-doped drain-gate region overlap capacitance");

    addPar ("CKAPPA", 0.6, false, NO_DEP,
      &MOSFET_B3::Model::ckappa,
      NULL, U_FARADMM1, CAT_CAP, "Coefficient for lightly doped region overlap capacitance fireing field capacitance");

    addPar ("CF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cf,
      NULL, U_FARADMM1, CAT_CAP, "Firing field capacitance");

    addPar ("VFBCV", -1.0, false, NO_DEP,
      &MOSFET_B3::Model::vfbcv,
      NULL, U_VOLT, CAT_CAP, "Flat-band voltage parameter (for CAPMOD = 0 only)");

    addPar ("CLC", 0.1e-6, false, NO_DEP,
      &MOSFET_B3::Model::clc,
      NULL, U_METER, CAT_CAP, "Constant term for short-channel model");

    addPar ("CLE", 0.6, false, NO_DEP,
      &MOSFET_B3::Model::cle,
      NULL, U_NONE, CAT_CAP, "Exponetial term for the short-channel model");

    addPar ("DWC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::dwc,
      NULL, U_METER, CAT_CAP, "Width offset fitting parameter from C-V");

    addPar ("DLC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::dlc,
      NULL, U_METER, CAT_CAP, "Length offset fitting parameter from C-V");

    addPar ("NOFF", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::noff,
      NULL, U_NONE, CAT_CAP, "CV parameter in Vgsteff, CV for weak to strong inversion");

    addPar ("VOFFCV", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::voffcv,
      NULL, U_VOLT, CAT_CAP, "CV parameter in Vgsteff, CV for weak to strong inversion");

    addPar ("ACDE", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::acde,
      NULL, U_MVM1, CAT_CAP, "Exponetial coefficient for charge thickness in capmod = 3 for accumulation and depletion regions");

    addPar ("MOIN", 15.0, false, NO_DEP,
      &MOSFET_B3::Model::moin,
      NULL, U_NONE, CAT_CAP, "Coefficient for the gate-bias dependent surface potential");

    addPar ("TCJ", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::tcj,
      NULL, U_KM1, CAT_TEMP, "Temperature coefficient of Cj");

    addPar ("TCJSW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::tcjsw,
      NULL, U_KM1, CAT_TEMP, "Temperature coefficient of Cswj");

    addPar ("TCJSWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::tcjswg,
      NULL, U_KM1, CAT_TEMP, "Temperature coefficient of Cjswg");

    addPar ("TPB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::tpb,
      NULL, U_VKM1, CAT_TEMP, "Temperature coefficient of Pb");

    addPar ("TPBSW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::tpbsw,
      NULL, U_VKM1, CAT_TEMP, "Temperature coefficient of Pbsw");

    addPar ("TPBSWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::tpbswg,
      NULL, U_VKM1, CAT_TEMP, "Temperature coefficient of Pbswg");

    addPar ("LCDSC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lcdsc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of cdsc");

    addPar ("LCDSCB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lcdscb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of cdscb");

    addPar ("LCDSCD", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lcdscd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of cdscd");

    addPar ("LCIT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lcit,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of cit");

    addPar ("LNFACTOR", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lnfactor,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of nfactor");

    addPar ("LXJ", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lxj,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of xj");

    addPar ("LVSAT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lvsat,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of vsat");

    addPar ("LAT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lat,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of at");

    addPar ("LA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::la0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of a0");

    addPar ("LAGS", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lags,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ags");

    addPar ("LA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::la1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of a1");

    addPar ("LA2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::la2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of a2");

    addPar ("LKETA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lketa,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of keta");

    addPar ("LNSUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lnsub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of nsub");

    addPar ("LNCH", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lnpeak,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of nch");

    addPar ("LNGATE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lngate,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ngate");

    addPar ("LGAMMA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lgamma1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of gamma1");

    addPar ("LGAMMA2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lgamma2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of gamma2");

    addPar ("LVBX", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lvbx,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VBX");

    addPar ("LVBM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lvbm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VBM");

    addPar ("LXT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lxt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of XT");

    addPar ("LK1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lk1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of K1");

    addPar ("LKT1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lkt1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of KT1");

    addPar ("LKT1L", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lkt1l,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of KT1L");

    addPar ("LKT2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lkt2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of KT2");

    addPar ("LK2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lk2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of K2");

    addPar ("LK3", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lk3,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of K3");

    addPar ("LK3B", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lk3b,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of K3B");

    addPar ("LW0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lw0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of W0");

    addPar ("LNLX", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lnlx,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NLX");

    addPar ("LDVT0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldvt0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT0");

    addPar ("LDVT1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldvt1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT1");

    addPar ("LDVT2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldvt2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT2");

    addPar ("LDVT0W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldvt0w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT0W");

    addPar ("LDVT1W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldvt1w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT1W");

    addPar ("LDVT2W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldvt2w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DVT2W");

    addPar ("LDROUT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldrout,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DROUT");

    addPar ("LDSUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldsub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of LDSUB");

    addPar ("LVTH0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lvth0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VT0");

    addPar ("LUA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lua,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of UA");

    addPar ("LUA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lua1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of UA1");

    addPar ("LUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of UB");

    addPar ("LUB1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lub1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of UB1");

    addPar ("LUC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::luc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of UC");

    addPar ("LUC1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::luc1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of UC1");

    addPar ("LU0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lu0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of U0");

    addPar ("LUTE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lute,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of UTE");

    addPar ("LVOFF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lvoff,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VOFF");

    addPar ("LRDSW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lrdsw,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of RDSW");

    addPar ("LPRWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lprwg,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PRWG");

    addPar ("LPRWB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lprwb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PRWB");

    addPar ("LPRT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lprt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PRT");

    addPar ("LETA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::leta0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ETA0");

    addPar ("LETAB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::letab,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ETAB");

    addPar ("LPCLM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lpclm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PCLM");

    addPar ("LPDIBLC1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lpdibl1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PDIBLC1");

    addPar ("LPDIBLC2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lpdibl2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PDIBLC2");

    addPar ("LPDIBLCB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lpdiblb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PDIBLCB");

    addPar ("LPSCBE1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lpscbe1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PSCBE1");

    addPar ("LPSCBE2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lpscbe2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PSCBE2");

    addPar ("LPVAG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lpvag,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of PVAG");

    addPar ("LDELTA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldelta,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DELTA");

    addPar ("LWR", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lwr,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of WR");

    addPar ("LDWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldwg,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DWG");

    addPar ("LDWB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ldwb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of DWB");

    addPar ("LB0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lb0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of B0");

    addPar ("LB1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lb1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of B1");

    addPar ("LALPHA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lalpha0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ALPHA0");

    addPar ("LALPHA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lalpha1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ALPHA1");

    addPar ("LBETA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lbeta0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of BETA0");

    addPar ("LVFB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lvfb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VFB");

    addPar ("LELM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lelm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ELM");

    addPar ("LCGSL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lcgsl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CGSL");

    addPar ("LCGDL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lcgdl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CGDL");

    addPar ("LCKAPPA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lckappa,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CKAPPA");

    addPar ("LCF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lcf,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CF");

    addPar ("LCLC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lclc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CLC");

    addPar ("LCLE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lcle,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of CLE");

    addPar ("LVFBCV", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lvfbcv,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VFBCV");

    addPar ("LNOFF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lnoff,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of NOFF");

    addPar ("LVOFFCV", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lvoffcv,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of VOFFCV");

    addPar ("LACDE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lacde,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of ACDE");

    addPar ("LMOIN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::lmoin,
      NULL, U_INVALID, CAT_DEPENDENCY, "Length dependence of MOIN");

    addPar ("WCDSC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wcdsc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CDSC");

    addPar ("WCDSCB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wcdscb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CDSCB");

    addPar ("WCDSCD", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wcdscd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CDSCD");

    addPar ("WCIT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wcit,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CIT");

    addPar ("WNFACTOR", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wnfactor,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NFACTOR");

    addPar ("WXJ", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wxj,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of XJ");

    addPar ("WVSAT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wvsat,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VSAT");

    addPar ("WAT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wat,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of AT");

    addPar ("WA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wa0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of A0");

    addPar ("WAGS", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wags,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of AGS");

    addPar ("WA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wa1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of A1");

    addPar ("WA2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wa2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of A2");

    addPar ("WKETA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wketa,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of KETA");

    addPar ("WNSUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wnsub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NSUB");

    addPar ("WNCH", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wnpeak,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NCH");

    addPar ("WNGATE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wngate,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NGATE");

    addPar ("WGAMMA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wgamma1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of GAMMA1");

    addPar ("WGAMMA2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wgamma2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of GAMMA2");

    addPar ("WVBX", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wvbx,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VBX");

    addPar ("WVBM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wvbm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VBM");

    addPar ("WXT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wxt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of XT");

    addPar ("WK1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wk1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of K1");

    addPar ("WKT1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wkt1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of KT1");

    addPar ("WKT1L", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wkt1l,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of KT1L");

    addPar ("WKT2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wkt2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of KT2");

    addPar ("WK2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wk2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of K2");

    addPar ("WK3", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wk3,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of K3");

    addPar ("WK3B", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wk3b,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of K3B");

    addPar ("WW0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ww0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of W0");

    addPar ("WNLX", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wnlx,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NLX");

    addPar ("WDVT0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdvt0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT0");

    addPar ("WDVT1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdvt1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT1");

    addPar ("WDVT2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdvt2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT2");

    addPar ("WDVT0W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdvt0w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT0W");

    addPar ("WDVT1W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdvt1w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT1W");

    addPar ("WDVT2W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdvt2w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DVT2W");

    addPar ("WDROUT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdrout,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DROUT");

    addPar ("WDSUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdsub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DSUB");

    addPar ("WVTH0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wvth0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VTO");

    addPar ("WUA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wua,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of UA");

    addPar ("WUA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wua1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of UA1");

    addPar ("WUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of UB");

    addPar ("WUB1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wub1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of UB1");

    addPar ("WUC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wuc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of UC");

    addPar ("WUC1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wuc1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of UC1");

    addPar ("WU0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wu0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of U0");

    addPar ("WUTE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wute,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of UTE");

    addPar ("WVOFF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wvoff,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VOFF");

    addPar ("WRDSW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wrdsw,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of RDSW");

    addPar ("WPRWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wprwg,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PRWG");

    addPar ("WPRWB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wprwb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PRWB");

    addPar ("WPRT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wprt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PRT");

    addPar ("WETA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::weta0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ETA0");

    addPar ("WETAB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wetab,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ETAB");

    addPar ("WPCLM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wpclm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PCLM");

    addPar ("WPDIBLC1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wpdibl1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PDIBLC1");

    addPar ("WPDIBLC2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wpdibl2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PDIBLC2");

    addPar ("WPDIBLCB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wpdiblb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PDIBLCB");

    addPar ("WPSCBE1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wpscbe1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PSCBE1");

    addPar ("WPSCBE2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wpscbe2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PSCBE2");

    addPar ("WPVAG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wpvag,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of PVAG");

    addPar ("WDELTA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdelta,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DELTA");

    addPar ("WWR", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wwr,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of WR");

    addPar ("WDWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdwg,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of WG");

    addPar ("WDWB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wdwb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of DWB");

    addPar ("WB0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wb0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of B0");

    addPar ("WB1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wb1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of B1");

    addPar ("WALPHA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::walpha0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ALPHA0");

    addPar ("WALPHA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::walpha1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ALPHA1");

    addPar ("WBETA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wbeta0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of BETA0");

    addPar ("WVFB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wvfb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VFB");

    addPar ("WELM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::welm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ELM");

    addPar ("WCGSL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wcgsl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CGSL");

    addPar ("WCGDL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wcgdl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CGDL");

    addPar ("WCKAPPA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wckappa,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CKAPPA");

    addPar ("WCF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wcf,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CF");

    addPar ("WCLC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wclc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CLC");

    addPar ("WCLE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wcle,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of CLE");

    addPar ("WVFBCV", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wvfbcv,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VFBCV");

    addPar ("WNOFF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wnoff,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of NOFF");

    addPar ("WVOFFCV", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wvoffcv,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of VOFFCV");

    addPar ("WACDE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wacde,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of ACDE");

    addPar ("WMOIN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::wmoin,
      NULL, U_INVALID, CAT_DEPENDENCY, "Width dependence of MOIN");

    addPar ("PCDSC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pcdsc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CDSC");

    addPar ("PCDSCB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pcdscb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CDSCB");

    addPar ("PCDSCD", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pcdscd,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CDSCD");

    addPar ("PCIT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pcit,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CIT");

    addPar ("PNFACTOR", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pnfactor,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NFACTOR");

    addPar ("PXJ", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pxj,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XJ");

    addPar ("PVSAT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvsat,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VSAT");

    addPar ("PAT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pat,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AT");

    addPar ("PA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pa0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of A0");

    addPar ("PAGS", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pags,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of AGS");

    addPar ("PA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pa1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of A1");

    addPar ("PA2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pa2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of A2");

    addPar ("PKETA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pketa,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KETA");

    addPar ("PNSUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pnsub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NSUB");

    addPar ("PNCH", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pnpeak,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NCH");

    addPar ("PNGATE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pngate,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NGATE");

    addPar ("PGAMMA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pgamma1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of GAMMA1");

    addPar ("PGAMMA2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pgamma2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of GAMMA2");

    addPar ("PVBX", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvbx,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VBX");

    addPar ("PVBM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvbm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VBM");

    addPar ("PXT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pxt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of XT");

    addPar ("PK1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pk1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K1");

    addPar ("PKT1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pkt1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KT1");

    addPar ("PKT1L", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pkt1l,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KT1L");

    addPar ("PKT2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pkt2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of KT2");

    addPar ("PK2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pk2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K2");

    addPar ("PK3", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pk3,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K3");

    addPar ("PK3B", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pk3b,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of K3B");

    addPar ("PW0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pw0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of W0");

    addPar ("PNLX", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pnlx,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NLX");

    addPar ("PDVT0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdvt0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT0");

    addPar ("PDVT1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdvt1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT1");

    addPar ("PDVT2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdvt2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT2");

    addPar ("PDVT0W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdvt0w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT0W");

    addPar ("PDVT1W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdvt1w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT1W");

    addPar ("PDVT2W", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdvt2w,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DVT2W");

    addPar ("PDROUT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdrout,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DROUT");

    addPar ("PDSUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdsub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DSUB");

    addPar ("PVTH0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvth0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VT0");

    addPar ("PUA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pua,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UA");

    addPar ("PUA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pua1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UA1");

    addPar ("PUB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pub,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UB");

    addPar ("PUB1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pub1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UB1");

    addPar ("PUC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::puc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UC");

    addPar ("PUC1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::puc1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UC1");

    addPar ("PU0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pu0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of U0");

    addPar ("PUTE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pute,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of UTE");

    addPar ("PVOFF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvoff,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VOFF");

    addPar ("PRDSW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::prdsw,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of RDSW");

    addPar ("PPRWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pprwg,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PRWG");

    addPar ("PPRWB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pprwb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PRWB");

    addPar ("PPRT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pprt,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PRT");

    addPar ("PETA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::peta0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ETA0");

    addPar ("PETAB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::petab,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ETAB");

    addPar ("PPCLM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ppclm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PCLM");

    addPar ("PPDIBLC1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ppdibl1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PDIBLC1");

    addPar ("PPDIBLC2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ppdibl2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PDIBLC2");

    addPar ("PPDIBLCB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ppdiblb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PDIBLCB");

    addPar ("PPSCBE1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ppscbe1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PSCBE1");

    addPar ("PPSCBE2", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ppscbe2,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PSCBE2");

    addPar ("PPVAG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::ppvag,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of PVAG");

    addPar ("PDELTA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdelta,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DELTA");

    addPar ("PWR", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pwr,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of WR");

    addPar ("PDWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdwg,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DWG");

    addPar ("PDWB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pdwb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of DWB");

    addPar ("PB0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pb0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of B0");

    addPar ("PB1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pb1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of B1");

    addPar ("PALPHA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::palpha0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ALPHA0");

    addPar ("PALPHA1", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::palpha1,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ALPHA1");

    addPar ("PBETA0", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pbeta0,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of BETA0");

    addPar ("PVFB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvfb,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VFB");

    addPar ("PELM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pelm,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ELM");

    addPar ("PCGSL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pcgsl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CGSL");

    addPar ("PCGDL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pcgdl,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CGDL");

    addPar ("PCKAPPA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pckappa,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CKAPPA");

    addPar ("PCF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pcf,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CF");

    addPar ("PCLC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pclc,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CLC");

    addPar ("PCLE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pcle,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of CLE");

    addPar ("PVFBCV", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvfbcv,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VFBCV");

    addPar ("PNOFF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pnoff,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of NOFF");

    addPar ("PVOFFCV", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pvoffcv,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of VOFFCV");

    addPar ("PACDE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pacde,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of ACDE");

    addPar ("PMOIN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::pmoin,
      NULL, U_INVALID, CAT_DEPENDENCY, "Cross-term dependence of MOIN");

#ifdef FRINGE_DONE
    addPar ("USEFRINGE", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::useFring,
      NULL, U_UNKNOWN, CAT_INVALID, "NOT in BSIM3");

#endif
    addPar ("TNOM", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::tnom,
      NULL, STANDARD, CAT_NONE, "");

    addPar ("CGSO", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cgso,
      NULL, U_FARADMM1, CAT_CAP, "Non-LLD region source-gate overlap capacitance per unit channel length");

    addPar ("CGDO", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cgdo,
      NULL, U_FARADMM1, CAT_CAP, "Non-LLD region drain-gate overlap capacitance per unit channel length");

    addPar ("CGBO", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::cgbo,
      NULL, U_FARADMM1, CAT_CAP, "Gate-bulk overlap capacitance per unit channel length");

    addPar ("XPART", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::xpart,
      NULL, U_NONE, CAT_CAP, "Charge partitioning rate flag");

    addPar ("RSH", 0.0, false, MIN_RES,
      &MOSFET_B3::Model::sheetResistance,
      NULL, U_OHM, CAT_RES, "Drain, source diffusion sheet resistance");

    addPar ("JS", 1.e-4, false, NO_DEP,
      &MOSFET_B3::Model::jctSatCurDensity,
      NULL, U_AMPMM2, CAT_PROCESS, "Bulk p-n saturation current density");

    addPar ("JSW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::jctSidewallSatCurDensity,
      NULL, U_AMPMM1, CAT_DC, "Sidewall saturation current per unit length");

    addPar ("PB", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::bulkJctPotential,
      NULL, U_VOLT, CAT_VOLT, "Bulk p-n bottom potential");

    addPar ("MJ", 0.5, false, NO_DEP,
      &MOSFET_B3::Model::bulkJctBotGradingCoeff,
      NULL, U_NONE, CAT_DOPING, "Bulk p-n bottom grading coefficient");

    addPar ("PBSW", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::sidewallJctPotential,
      NULL, U_VOLT, CAT_CAP, "Source/drain side junction built-in potential");

    addPar ("PBSWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::GatesidewallJctPotential,
      NULL, U_VOLT, CAT_CAP, "Source/drain gate sidewall junction built-in potential");

    addPar ("MJSW", 0.33, false, NO_DEP,
      &MOSFET_B3::Model::bulkJctSideGradingCoeff,
      NULL, U_NONE, CAT_DOPING, "Bulk p-n sidewall grading coefficient");

    addPar ("CJ", 5.e-4, false, NO_DEP,
      &MOSFET_B3::Model::unitAreaJctCap,
      NULL, U_FARADMM2, CAT_CAP, "Bulk p-n zero-bias bottom capacitance/area");

    addPar ("CJSW", 5.e-10, false, NO_DEP,
      &MOSFET_B3::Model::unitLengthSidewallJctCap,
      NULL, U_FARADMM2, CAT_CAP, "Bulk p-n zero-bias sidewall capacitance/area");

    addPar ("MJSWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::bulkJctGateSideGradingCoeff,
      NULL, U_NONE, CAT_CAP, "Source/grain gate sidewall junction capacitance grading coeficient");

    addPar ("CJSWG", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::unitLengthGateSidewallJctCap,
      NULL, U_FARADMM1, CAT_CAP, "Source/grain gate sidewall junction capacitance per unit width");

    addPar ("NJ", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::jctEmissionCoeff,
      NULL, U_NONE, CAT_TEMP, "Emission coefficient of junction");

    addPar ("XTI", 3.0, false, NO_DEP,
      &MOSFET_B3::Model::jctTempExponent,
      NULL, U_NONE, CAT_TEMP, "Junction current temperature exponent coefficient");

    addPar ("NOIA", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::oxideTrapDensityA,
      NULL, U_NONE, CAT_FLICKER, "Noise parameter a");

    addPar ("NOIB", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::oxideTrapDensityB,
      NULL, U_NONE, CAT_FLICKER, "Noise parameter b");

    addPar ("NOIC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::oxideTrapDensityC,
      NULL, U_NONE, CAT_FLICKER, "Noise parameter c");

    addPar ("EM", 4.1e7, false, NO_DEP,
      &MOSFET_B3::Model::em,
      NULL, U_VMM1, CAT_FLICKER, "Saturation field");

    addPar ("EF", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::ef,
      NULL, U_NONE, CAT_FLICKER, "Flicker exponent");

    addPar ("AF", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::af,
      NULL, U_NONE, CAT_FLICKER, "Flicker noise exponent");

    addPar ("KF", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::kf,
      NULL, U_NONE, CAT_FLICKER, "Flicker noise coefficient");

    addPar ("LINT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Lint,
      NULL, U_METER, CAT_DC, "Length of offset fiting parameter from I-V without bias");

    addPar ("LL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Ll,
      NULL, U_MEXPLL, CAT_GEOMETRY, "Coefficient of length dependence for length offset");

    addPar ("LLC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Llc,
      NULL, U_MEXPLL, CAT_GEOMETRY, "Coefficient of length dependence for CV channel length offset");

    addPar ("LLN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Lln,
      NULL, U_NONE, CAT_GEOMETRY, "Power of length dependence for length offset");

    addPar ("LW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Lw,
      NULL, U_MEXPLW, CAT_GEOMETRY, "Coefficient of width dependence for length offset");

    addPar ("LWC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Lwc,
      NULL, U_MEXPLW, CAT_GEOMETRY, "Coefficient of width dependence for channel length offset");

    addPar ("LWN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Lwn,
      NULL, U_NONE, CAT_GEOMETRY, "Power of width dependence for length offset");

    addPar ("LWL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Lwl,
      NULL, U_MEXPLLLW, CAT_GEOMETRY, "Coefficient of length and width cross term for length offset");

    addPar ("LWLC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Lwlc,
      NULL, U_MEXPLLLW, CAT_GEOMETRY, "Coefficient of length and width dependence for CV channel length offset");

    addPar ("WINT", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wint,
      NULL, U_METER, CAT_DC, "Width-offset fitting parameter from I-V without bias");

    addPar ("WL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wl,
      NULL, U_MEXPWL, CAT_GEOMETRY, "Coefficient of length dependence for width offset");

    addPar ("WLC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wlc,
      NULL, U_MEXPWL, CAT_GEOMETRY, "Coefficient of length dependence for CV channel width offset");

    addPar ("WLN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wln,
      NULL, U_NONE, CAT_GEOMETRY, "Power of length dependece of width offset");

    addPar ("WW", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Ww,
      NULL, U_MEXPWW, CAT_GEOMETRY, "Coefficient of width dependence for width offset");

    addPar ("WWC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wwc,
      NULL, U_MEXPWW, CAT_GEOMETRY, "Coefficient of width dependence for CV channel width offset");

    addPar ("WWN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wwn,
      NULL, U_NONE, CAT_GEOMETRY, "Power of width dependence of width offset");

    addPar ("WWL", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wwl,
      NULL, U_MEXPWLWW, CAT_GEOMETRY, "Coefficient of length and width cross term for width offset");

    addPar ("WWLC", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wwlc,
      NULL, U_MEXPWLWW, CAT_GEOMETRY, "Coefficient of length and width dependence for CV channel width offset");

    addPar ("L", 5.e-6, false, NO_DEP,
      &MOSFET_B3::Model::model_l,
      NULL, U_METER,  CAT_GEOMETRY, "Channel length");

    addPar ("W", 5.e-6, false, NO_DEP,
      &MOSFET_B3::Model::model_w,
      NULL, U_METER,  CAT_GEOMETRY, "Channel width");

    addPar ("LMAX", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::Lmax,
      NULL, U_METER, CAT_BIN, "Maximum channel length");

    addPar ("LMIN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Lmin,
      NULL, U_METER, CAT_BIN, "Minimum channel length");

    addPar ("WMAX", 1.0, false, NO_DEP,
      &MOSFET_B3::Model::Wmax,
      NULL, U_METER, CAT_BIN, "Maximum channel width");

    addPar ("WMIN", 0.0, false, NO_DEP,
      &MOSFET_B3::Model::Wmin,
      NULL, U_METER, CAT_BIN, "Minimum channel width");

    // Set up exceptions (ie variables that are not doubles):
    addPar ("MOBMOD", 1, false, NO_DEP,
            &MOSFET_B3::Model::mobMod, NULL,
            U_NONE, CAT_CONTROL, "Mobility model selector");
    addPar ("BINUNIT", 1, false, NO_DEP,
            &MOSFET_B3::Model::binUnit, NULL,
            U_NONE, CAT_CONTROL, "Binning unit selector");
    addPar ("CAPMOD", 3, false, NO_DEP,
            &MOSFET_B3::Model::capMod, NULL,
            U_NONE, CAT_CONTROL, "Flag for capacitance models");
    addPar ("PARAMCHK", 0, false, NO_DEP,
            &MOSFET_B3::Model::paramChk, NULL,
            U_NONE, CAT_CONTROL, "Parameter value check");
    addPar ("NOIMOD", 1, false, NO_DEP,
            &MOSFET_B3::Model::noiMod, NULL,
            U_NONE, CAT_CONTROL, "Flag for noise models");
    addPar ("VERSION", string("3.2.2"), false, NO_DEP,
            &MOSFET_B3::Model::version, NULL,
            U_NONE, CAT_CONTROL, "Version number");
}

namespace MOSFET_B3 {

vector< vector<int> > Instance::jacStamp_DC_SC;
vector< vector<int> > Instance::jacStamp_DC;
vector< vector<int> > Instance::jacStamp_SC;
vector< vector<int> > Instance::jacStamp;

vector<int> Instance::jacMap_DC_SC;
vector<int> Instance::jacMap_DC;
vector<int> Instance::jacMap_SC;
vector<int> Instance::jacMap;

vector< vector<int> > Instance::jacMap2_DC_SC;
vector< vector<int> > Instance::jacMap2_DC;
vector< vector<int> > Instance::jacMap2_SC;
vector< vector<int> > Instance::jacMap2;

ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{

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

  numIntVars = 0;
  if ( sourceConductance!=0.0 ) ++numIntVars;
  if ( drainConductance!=0.0 ) ++numIntVars;
  if ( nqsMod ) ++numIntVars;

  if (icVDSGiven) ++numIntVars;
  if (icVGSGiven) ++numIntVars;
  if (icVBSGiven) ++numIntVars;


  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/00
//-----------------------------------------------------------------------------
Instance::Instance(
  InstanceBlock & IB,
  Model &model,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)
  : DeviceInstance              (IB,mlData1,ss1,ed1,do1),
    model_                          (model),
    dNode                             (0),
    gNode                             (0),
    sNode                             (0),
    bNode                             (0),
    dNodePrime                        (0),
    sNodePrime                        (0),
    qNode                             (0),
    ueff                              (0.0),
    thetavth                          (0.0),
    von                               (0.0),
    vdsat                             (0.0),
    cgdo                              (0.0),
    cgso                              (0.0),
    vjsm                              (0.0),
    IsEvjsm                           (0.0),
    vjdm                              (0.0),
    IsEvjdm                           (0.0),
    l                                 (devOptions.defl),
    w                                 (devOptions.defw),
    drainArea                         (devOptions.defad),
    sourceArea                        (devOptions.defas),
    drainSquares                      (1.0),
    sourceSquares                     (1.0),
    drainPerimeter                    (0.0),
    sourcePerimeter                   (0.0),
    sourceConductance                 (0.0),
    drainConductance                  (0.0),
    icVBS                             (0.0),
    icVDS                             (0.0),
    icVGS                             (0.0),
    OFF                               (false),
    mode                              (0),
    nqsMod                            (0),
    numberParallel                    (1.0),
    qinv                              (0.0),
    cd                                (0.0),
    cbs                               (0.0),
    cbd                               (0.0),
    csub                              (0.0),
    cdrain                            (0.0),
    gm                                (0.0),
    gds                               (0.0),
    gmbs                              (0.0),
    gbd                               (0.0),
    gbs                               (0.0),
    gbbs                              (0.0),
    gbgs                              (0.0),
    gbds                              (0.0),
    cggb                              (0.0),
    cgdb                              (0.0),
    cgsb                              (0.0),
    cbgb                              (0.0),
    cbdb                              (0.0),
    cbsb                              (0.0),
    cdgb                              (0.0),
    cddb                              (0.0),
    cdsb                              (0.0),
    capbd                             (0.0),
    capbs                             (0.0),
    cqgb                              (0.0),
    cqdb                              (0.0),
    cqsb                              (0.0),
    cqbb                              (0.0),
    qgate                             (0.0),
    qbulk                             (0.0),
    qdrn                              (0.0),
    gtau                              (0.0),
    gtg                               (0.0),
    gtd                               (0.0),
    gts                               (0.0),
    gtb                               (0.0),
    limitedFlag                       (false),
    paramPtr                          (NULL),
    icVBSGiven                        (0),
    icVDSGiven                        (0),
    icVGSGiven                        (0),
    temp                              (devOptions.temp.dVal()),
    ChargeComputationNeeded           (true),
    gcbdb                             (0.0),
    gcbgb                             (0.0),
    gcbsb                             (0.0),
    gcddb                             (0.0),
    gcdgb                             (0.0),
    gcdsb                             (0.0),
    gcgdb                             (0.0),
    gcggb                             (0.0),
    gcgsb                             (0.0),
    gcsdb                             (0.0),
    gcsgb                             (0.0),
    gcssb                             (0.0),
    qgd                               (0.0),
    qgs                               (0.0),
    qgb                               (0.0),
    qgdo                              (0.0),
    qgso                              (0.0),
    qsrc                              (0.0),
    CoxWL                             (0.0),
    Cgg                               (0.0),
    Cgd                               (0.0),
    Cgb                               (0.0),
    Cdg                               (0.0),
    Cdd                               (0.0),
    Cds                               (0.0),
    Csg                               (0.0),
    Csd                               (0.0),
    Css                               (0.0),
    Csb                               (0.0),
    Cbg                               (0.0),
    Cbd                               (0.0),
    Cbb                               (0.0),
    CAPcggb                           (0.0),
    CAPcgdb                           (0.0),
    CAPcgsb                           (0.0),
    CAPcbgb                           (0.0),
    CAPcbdb                           (0.0),
    CAPcbsb                           (0.0),
    CAPcdgb                           (0.0),
    CAPcddb                           (0.0),
    CAPcdsb                           (0.0),
    CAPcsgb                           (0.0),
    CAPcsdb                           (0.0),
    CAPcssb                           (0.0),
    Qeqqd_Jdxp                        (0.0),
    Qeqqb_Jdxp                        (0.0),
    Qeqqg_Jdxp                        (0.0),
    dxpart                            (0.0),
    sxpart                            (0.0),
    ggtg                              (0.0),
    ggtd                              (0.0),
    ggts                              (0.0),
    ggtb                              (0.0),
    ddxpart_dVd                       (0.0),
    ddxpart_dVg                       (0.0),
    ddxpart_dVb                       (0.0),
    ddxpart_dVs                       (0.0),
    dsxpart_dVd                       (0.0),
    dsxpart_dVg                       (0.0),
    dsxpart_dVb                       (0.0),
    dsxpart_dVs                       (0.0),
    gbspsp                            (0.0),
    gbbdp                             (0.0),
    gbbsp                             (0.0),
    gbspg                             (0.0),
    gbspb                             (0.0),
    gbspdp                            (0.0),
    gbdpdp                            (0.0),
    gbdpg                             (0.0),
    gbdpb                             (0.0),
    gbdpsp                            (0.0),
    cdreq                             (0.0),
    ceqbd                             (0.0),
    ceqbs                             (0.0),
    Gm                                (0.0),
    Gmbs                              (0.0),
    FwdSum                            (0.0),
    RevSum                            (0.0),
    T1global                          (0.0),
    dVgst_dVg                         (0.0),
    dVgst_dVb                         (0.0),
    dVgs_eff_dVg                      (0.0),
    dDeltaPhi_dVg                     (0.0),
    dDeltaPhi_dVd                     (0.0),
    dDeltaPhi_dVb                     (0.0),
    vtm                               (0.0),
    jctTempSatCurDensity              (0.0),
    jctSidewallTempSatCurDensity      (0.0),
    unitAreaJctCapTemp                (0.0),
    unitLengthSidewallJctCapTemp      (0.0),
    unitLengthGateSidewallJctCapTemp  (0.0),
    PhiBTemp                          (0.0),
    PhiBSWTemp                        (0.0),
    PhiBSWGTemp                       (0.0),
    Vd                                (0.0),
    Vs                                (0.0),
    Vg                                (0.0),
    Vb                                (0.0),
    Vsp                               (0.0),
    Vdp                               (0.0),
    Qtotal                            (0.0),
    Vddp                              (0.0),
    Vssp                              (0.0),
    Vbsp                              (0.0),
    Vbdp                              (0.0),
    Vgsp                              (0.0),
    Vgdp                              (0.0),
    Vgb                               (0.0),
    Vdpsp                             (0.0),
    Idrain                            (0.0),
    Isource                           (0.0),
    df1dVdp                           (0.0),
    df2dVdp                           (0.0),
    df1dVsp                           (0.0),
    df2dVsp                           (0.0),
    df1dVg                            (0.0),
    df2dVg                            (0.0),
    df1dVb                            (0.0),
    df2dVb                            (0.0),
    vbd_old                           (0.0),
    vbs_old                           (0.0),
    vgs_old                           (0.0),
    vds_old                           (0.0),
    vgs_orig                          (0.0),
    vds_orig                          (0.0),
    vbs_orig                          (0.0),
    vbd_orig                          (0.0),
    vgd_orig                          (0.0),
    newtonIterOld                     (0),
    cqdef             (0.0),
    vgb               (0.0),
    vgd               (0.0),
    vbd               (0.0),
    vbs               (0.0),
    vgs               (0.0),
    vds               (0.0),
    qb                (0.0),
    qg                (0.0),
    qd                (0.0),
    qbs               (0.0),
    qbd               (0.0),
    qcheq             (0.0),
    qcdump            (0.0),
    qdef              (0.0),
    gqdef             (0.0),
    gcqdb             (0.0),
    gcqsb             (0.0),
    gcqgb             (0.0),
    gcqbb             (0.0),
    ScalingFactor     (0.0),
    cqgate            (0.0),
    cqbulk            (0.0),
    cqdrn             (0.0),
    cdreq_Jdxp        (0.0),
    ceqbd_Jdxp        (0.0),
    ceqbs_Jdxp        (0.0),
    cqdef_Jdxp        (0.0),
    ceqqd_Jdxp        (0.0),
    ceqqb_Jdxp        (0.0),
    ceqqg_Jdxp        (0.0),
// matrix and vectors indices:
// state vector: (local indices)
    li_store_vbd             (-1),
    li_store_vbs             (-1),
    li_store_vgs             (-1),
    li_store_vds             (-1),
    li_store_von             (-1),
    li_store_dev_id          (-1),
    li_store_dev_is          (-1),
    li_store_dev_ig          (-1),
    li_store_dev_ib          (-1),
    li_state_qb              (-1),
    li_state_qg              (-1),
    li_state_qd              (-1),
    li_state_qbs             (-1),
    li_state_qbd             (-1),
    li_state_qcheq           (-1),
    li_state_qcdump          (-1),
    li_state_qdef            (-1),
// solution vector: (local indices)
    li_Drain                 (-1),
    li_Gate                  (-1),
    li_Source                (-1),
    li_Bulk                  (-1),
    li_DrainPrime            (-1),
    li_SourcePrime           (-1),
    li_Charge                (-1),
    li_Ibs                   (-1),
    li_Ids                   (-1),
    li_Igs                   (-1),
// matrix offsets:
// jacobian:
//  drain row
    ADrainEquDrainNodeOffset             (-1),
    ADrainEquDrainPrimeNodeOffset        (-1),
    ADrainEquIdsOffset                   (-1),
//  gate row
    AGateEquGateNodeOffset               (-1),
    AGateEquBulkNodeOffset               (-1),
    AGateEquDrainPrimeNodeOffset         (-1),
    AGateEquSourcePrimeNodeOffset        (-1),
    AGateEquChargeVarOffset              (-1),
    AGateEquIgsOffset                    (-1),
//  source row
    ASourceEquSourceNodeOffset           (-1),
    ASourceEquSourcePrimeNodeOffset      (-1),
    ASourceEquIbsOffset                  (-1),
    ASourceEquIdsOffset                  (-1),
    ASourceEquIgsOffset                  (-1),
//  bulk row
    ABulkEquGateNodeOffset               (-1),
    ABulkEquBulkNodeOffset               (-1),
    ABulkEquDrainPrimeNodeOffset         (-1),
    ABulkEquSourcePrimeNodeOffset        (-1),
    ABulkEquChargeVarOffset              (-1),
    ABulkEquIbsOffset                    (-1),
// drain' row
    ADrainPrimeEquDrainNodeOffset        (-1),
    ADrainPrimeEquGateNodeOffset         (-1),
    ADrainPrimeEquBulkNodeOffset         (-1),
    ADrainPrimeEquDrainPrimeNodeOffset   (-1),
    ADrainPrimeEquSourcePrimeNodeOffset  (-1),
    ADrainPrimeEquChargeVarOffset        (-1),
// source' row
    ASourcePrimeEquGateNodeOffset        (-1),
    ASourcePrimeEquSourceNodeOffset      (-1),
    ASourcePrimeEquBulkNodeOffset        (-1),
    ASourcePrimeEquDrainPrimeNodeOffset  (-1),
    ASourcePrimeEquSourcePrimeNodeOffset (-1),
    ASourcePrimeEquChargeVarOffset       (-1),
// Charge row
    AChargeEquChargeVarOffset            (-1),
    AChargeEquDrainPrimeNodeOffset       (-1),
    AChargeEquGateNodeOffset             (-1),
    AChargeEquSourcePrimeNodeOffset      (-1),
    AChargeEquBulkNodeOffset             (-1),
    icVBSEquVsOffset                     (-1),
    icVBSEquVbOffset                     (-1),
    icVBSEquIbsOffset                    (-1),
    icVDSEquVdOffset                     (-1),
    icVDSEquVsOffset                     (-1),
    icVDSEquIdsOffset                    (-1),
    icVGSEquVgOffset                     (-1),
    icVGSEquVsOffset                     (-1),
    icVGSEquIgsOffset                    (-1),
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
//
  // V_d Row:
    f_DrainEquDrainNodePtr(0),             // a
    f_DrainEquDrainPrimeNodePtr(0),        // b
    f_DrainEquIdsPtr(0),                   // i1

  // V_g Row:
    f_GateEquGateNodePtr(0),               // c
    f_GateEquBulkNodePtr(0),               // d
    f_GateEquDrainPrimeNodePtr(0),         // e
    f_GateEquSourcePrimeNodePtr(0),        // f
    f_GateEquChargeVarPtr(0),              // 1
    f_GateEquIgsPtr(0),                    // i2

  // V_s Row:
    f_SourceEquSourceNodePtr(0),           // g
    f_SourceEquSourcePrimeNodePtr(0),      // h
    f_SourceEquIbsPtr(0),                  // i3
    f_SourceEquIdsPtr(0),                  // i4
    f_SourceEquIgsPtr(0),                  // i5

  // V_b Row:
    f_BulkEquGateNodePtr(0),               // i
    f_BulkEquBulkNodePtr(0),               // j
    f_BulkEquDrainPrimeNodePtr(0),         // k
    f_BulkEquSourcePrimeNodePtr(0),        // l
    f_BulkEquChargeVarPtr(0),              // 2
    f_BulkEquIbsPtr(0),                    // i6

  // V_d' Row:
    f_DrainPrimeEquDrainNodePtr(0),        // m
    f_DrainPrimeEquGateNodePtr(0),         // n
    f_DrainPrimeEquBulkNodePtr(0),         // o
    f_DrainPrimeEquDrainPrimeNodePtr(0),   // p
    f_DrainPrimeEquSourcePrimeNodePtr(0),  // q
    f_DrainPrimeEquChargeVarPtr(0),        // 3

  // V_s' Row:
    f_SourcePrimeEquGateNodePtr(0),        // r
    f_SourcePrimeEquSourceNodePtr(0),      // s
    f_SourcePrimeEquBulkNodePtr(0),        // t
    f_SourcePrimeEquDrainPrimeNodePtr(0),  // u
    f_SourcePrimeEquSourcePrimeNodePtr(0), // v
    f_SourcePrimeEquChargeVarPtr(0),       // 4

  // MOSFET charge (Q) Row:
    f_ChargeEquChargeVarPtr(0),            // 5
    f_ChargeEquDrainPrimeNodePtr(0),       // 6
    f_ChargeEquGateNodePtr(0),             // 7
    f_ChargeEquSourcePrimeNodePtr(0),      // 8
    f_ChargeEquBulkNodePtr(0),             // 9

    // icVBS
    f_icVBSEquVsPtr(0),                   // i7
    f_icVBSEquVbPtr(0),                   // i8
    f_icVBSEquIbsPtr(0),                  // i9

    // icVDS
    f_icVDSEquVdPtr(0),                   // i10
    f_icVDSEquVsPtr(0),                   // i11
    f_icVDSEquIdsPtr(0),                  // i12

    // icVGS
    f_icVGSEquVgPtr(0),                   // i13
    f_icVGSEquVsPtr(0),                   // i14
    f_icVGSEquIgsPtr(0),                  // i15

  // V_d Row:
    q_DrainEquDrainNodePtr(0),             // a
    q_DrainEquDrainPrimeNodePtr(0),        // b
    q_DrainEquIdsPtr(0),                   // i1

  // V_g Row:
    q_GateEquGateNodePtr(0),               // c
    q_GateEquBulkNodePtr(0),               // d
    q_GateEquDrainPrimeNodePtr(0),         // e
    q_GateEquSourcePrimeNodePtr(0),        // f
    q_GateEquChargeVarPtr(0),              // 1
    q_GateEquIgsPtr(0),                    // i2

  // V_s Row:
    q_SourceEquSourceNodePtr(0),           // g
    q_SourceEquSourcePrimeNodePtr(0),      // h
    q_SourceEquIbsPtr(0),                  // i3
    q_SourceEquIdsPtr(0),                  // i4
    q_SourceEquIgsPtr(0),                  // i5

  // V_b Row:
    q_BulkEquGateNodePtr(0),               // i
    q_BulkEquBulkNodePtr(0),               // j
    q_BulkEquDrainPrimeNodePtr(0),         // k
    q_BulkEquSourcePrimeNodePtr(0),        // l
    q_BulkEquChargeVarPtr(0),              // 2
    q_BulkEquIbsPtr(0),                    // i6

  // V_d' Row:
    q_DrainPrimeEquDrainNodePtr(0),        // m
    q_DrainPrimeEquGateNodePtr(0),         // n
    q_DrainPrimeEquBulkNodePtr(0),         // o
    q_DrainPrimeEquDrainPrimeNodePtr(0),   // p
    q_DrainPrimeEquSourcePrimeNodePtr(0),  // q
    q_DrainPrimeEquChargeVarPtr(0),        // 3

  // V_s' Row:
    q_SourcePrimeEquGateNodePtr(0),        // r
    q_SourcePrimeEquSourceNodePtr(0),      // s
    q_SourcePrimeEquBulkNodePtr(0),        // t
    q_SourcePrimeEquDrainPrimeNodePtr(0),  // u
    q_SourcePrimeEquSourcePrimeNodePtr(0), // v
    q_SourcePrimeEquChargeVarPtr(0),       // 4

  // MOSFET charge (Q) Row:
    q_ChargeEquChargeVarPtr(0),            // 5
    q_ChargeEquDrainPrimeNodePtr(0),       // 6
    q_ChargeEquGateNodePtr(0),             // 7
    q_ChargeEquSourcePrimeNodePtr(0),      // 8
    q_ChargeEquBulkNodePtr(0),             // 9

    // icVBS
    q_icVBSEquVsPtr(0),                   // i7
    q_icVBSEquVbPtr(0),                   // i8
    q_icVBSEquIbsPtr(0),                  // i9

    // icVDS
    q_icVDSEquVdPtr(0),                   // i10
    q_icVDSEquVsPtr(0),                   // i11
    q_icVDSEquIdsPtr(0),                  // i12

    // icVGS
    q_icVGSEquVgPtr(0),                   // i13
    q_icVGSEquVsPtr(0),                   // i14
    q_icVGSEquIgsPtr(0),                  // i15
//
#endif
    blockHomotopyID                      (0),
    randomPerturb                        (0.0),
    // one last thing
    updateTemperatureCalled_ (false)
{
  numIntVars   = 3;
  numExtVars   = 4;
  numStateVars = 12;
  setNumStoreVars(5);
  numLeadCurrentStoreVars = 4; // drain, gate, source & base lead currents

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 2;
  devConMap[2] = 1;
  devConMap[3] = 3;

  setName(IB.getName());
  setModelName(model_.getName());

  blockHomotopyID =
    devSupport.getGainScaleBlockID(devOptions.numGainScaleBlocks);
  randomPerturb =
    devSupport.getRandomPerturbation();


  if( jacStamp.empty() )
  {
    jacStamp_DC_SC.resize(6);
    jacStamp_DC_SC[0].resize(2);
    jacStamp_DC_SC[0][0]=0;
    jacStamp_DC_SC[0][1]=4;
    jacStamp_DC_SC[1].resize(4);
    jacStamp_DC_SC[1][0]=1;
    jacStamp_DC_SC[1][1]=3;
    jacStamp_DC_SC[1][2]=4;
    jacStamp_DC_SC[1][3]=5;
    jacStamp_DC_SC[2].resize(2);
    jacStamp_DC_SC[2][0]=2;
    jacStamp_DC_SC[2][1]=5;
    jacStamp_DC_SC[3].resize(4);
    jacStamp_DC_SC[3][0]=1;
    jacStamp_DC_SC[3][1]=3;
    jacStamp_DC_SC[3][2]=4;
    jacStamp_DC_SC[3][3]=5;
    jacStamp_DC_SC[4].resize(5);
    jacStamp_DC_SC[4][0]=0;
    jacStamp_DC_SC[4][1]=1;
    jacStamp_DC_SC[4][2]=3;
    jacStamp_DC_SC[4][3]=4;
    jacStamp_DC_SC[4][4]=5;
    jacStamp_DC_SC[5].resize(5);
    jacStamp_DC_SC[5][0]=1;
    jacStamp_DC_SC[5][1]=2;
    jacStamp_DC_SC[5][2]=3;
    jacStamp_DC_SC[5][3]=4;
    jacStamp_DC_SC[5][4]=5;

    jacMap_DC_SC.clear();
    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_DC,    jacMap_DC, jacMap2_DC, 5, 2, 6);

    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_SC,    jacMap_SC, jacMap2_SC, 4, 0, 6);

    jacStampMap(jacStamp_DC, jacMap_DC, jacMap2_DC,
                jacStamp,    jacMap, jacMap2, 4, 0, 6);

#ifdef  Xyce_DEBUG_DEVICE
    if (devOptions.debugLevel > 0)
    {
      cout << "Instance::Instance jacStampMap_DS_SC" << endl;
      for (int k = 0; k<jacMap_DC_SC.size(); ++k )
      {
        cout << "jacStamp_DC_SC[ " << jacMap_DC_SC[k] << " ] = { ";
        for (int q = 0; q < jacMap2_DC_SC[k].size(); ++q )
        {
          cout << jacStamp_DC_SC[jacMap_DC_SC[k]][jacMap2_DC_SC[k][q]] << "  ";
        }
        cout << "}" << endl;
      }

      cout << "Instance::Instance jacStampMap_DS" << endl;
      for (int k = 0; k<jacMap_DC.size(); ++k )
      {
        cout << "jacStamp_DC[ " << jacMap_DC[k] << " ] = { ";
        for (int q = 0; q < jacMap2_DC[k].size(); ++q )
        {
          cout << jacStamp_DC[jacMap_DC[k]][jacMap2_DC[k][q]] << "  ";
        }
        cout << "}" << endl;
      }

      cout << "Instance::Instance jacStampMap_SC" << endl;
      for (int k = 0; k<jacMap_SC.size(); ++k )
      {
        cout << "jacStamp_SC[ " << jacMap_SC[k] << " ] = { ";
        for (int q = 0; q < jacMap2_SC[k].size(); ++q )
        {
          cout << jacStamp_SC[jacMap_SC[k]][jacMap2_SC[k][q]] << "  ";
        }
        cout << "}" << endl;
      }

      cout << "Instance::Instance jacStampMap" << endl;
      for (int k = 0; k<jacMap.size(); ++k )
      {
        cout << "jacStamp[ " << jacMap[k] << " ] = { ";
        for (int q = 0; q < jacMap2[k].size(); ++q )
        {
          cout << jacStamp[jacMap[k]][jacMap2[k][q]] << "  ";
        }
        cout << "}" << endl;
      }
    }
#endif
  }

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
    w = model_.model_w;
  if (!given("AD"))
    drainArea = devOptions.defad;
  if (!given("AS"))
    sourceArea = devOptions.defas;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  if (given("NQSMOD"))
  {
    string msg = "Instance::Instance";
    msg += "  nsqMod = 1.  Not allowed yet.  Setting to 0.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
  }

  if (devOptions.verboseLevel > 0 && (l > model_.Lmax || l < model_.Lmin))
  {
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n"
        << "Channel length out of range for the model " << getName();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, oss.str());
  }

  if (devOptions.verboseLevel > 0 && (w > model_.Wmax || w < model_.Wmin))
  {
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n"
        << "Channel width out of range for the model " << getName();
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, oss.str());
  }

  // if we need to, build the jacobian stamp and operating point
  // jacobian stamp specificly needed for this device's initial
  // conditions

  // There are 4 optional internal variables that may be included in
  // this device: charge Q, Vbs, Vds and Vgs.  We need to assign these
  // variables to columns in our jacobian depending on what's there.
  // For example if they're all required than Q is the 7th variable
  // and goes in column 6 (numbering from zero), Vbs is 8th in column
  // 7, Vds is 9th in column, Vgs is 10th in column 9.  If any
  // variable is not needed, those in higher columns shift
  // down. I.e. if Q isn't needed, then Vbs is in column 6, Vds in
  // column 7 and Vgs in column 8.

  int currentCol = 6;
  int numExtraCol = 0;
  int qCol = -1, icVBSCol = -1, icVDSCol = -1, icVGSCol = -1;
  if( nqsMod )
  {
    qCol = currentCol;
    ++currentCol;
    ++numExtraCol;
  }
  if( icVBSGiven )
  {
    icVBSCol = currentCol;
    ++currentCol;
    ++numExtraCol;
  }
  if( icVDSGiven )
  {
    icVDSCol = currentCol;
    ++currentCol;
    ++numExtraCol;
  }
  if( icVGSGiven )
  {
    icVGSCol = currentCol;
    ++currentCol;
    ++numExtraCol;
  }

  // ok now build this instance's jacobian stamp
  if( nqsMod || icVBSGiven || icVDSGiven || icVGSGiven )
  {
    // we need a special version of the jacStamp for the
    jacStampSpecial.resize(6 + numExtraCol);
    // row KCL d
    int numNonZeros = 2 + (icVDSGiven ? 1:0);
    jacStampSpecial[0].resize( numNonZeros );
    jacStampSpecial[0][0]=0;
    jacStampSpecial[0][1]=4;
    currentCol = 2;
    if( icVDSGiven )
    {
      jacStampSpecial[0][currentCol]=icVDSCol;
      ++currentCol;
    }

    // row KCL g
    numNonZeros = 4 + (nqsMod ? 1:0 ) + (icVGSGiven ? 1:0);
    jacStampSpecial[1].resize( numNonZeros );
    jacStampSpecial[1][0]=1;
    jacStampSpecial[1][1]=3;
    jacStampSpecial[1][2]=4;
    jacStampSpecial[1][3]=5;
    currentCol = 4;
    if( nqsMod )
    {
      jacStampSpecial[1][currentCol]=qCol;
      ++currentCol;
    }
    if( icVGSGiven )
    {
      jacStampSpecial[1][currentCol]=icVGSCol;
    }

    // row KCL s
    numNonZeros = 2 + (icVBSGiven ? 1:0) + (icVDSGiven ? 1:0) + (icVGSGiven ? 1:0);
    jacStampSpecial[2].resize( numNonZeros );
    jacStampSpecial[2][0]=2;
    jacStampSpecial[2][1]=5;
    currentCol = 2;
    if ( icVBSGiven )
    {
      jacStampSpecial[2][currentCol] = icVBSCol;
      ++currentCol;
    }
    if( icVDSGiven )
    {
      jacStampSpecial[2][currentCol] = icVDSCol;
      ++currentCol;
    }
    if( icVGSGiven )
    {
      jacStampSpecial[2][currentCol] = icVGSCol;
      ++currentCol;
    }

    // row KCL b
    numNonZeros = 4 + (nqsMod ? 1:0 ) + (icVBSGiven ? 1:0);
    jacStampSpecial[3].resize( numNonZeros );
    jacStampSpecial[3][0]=1;
    jacStampSpecial[3][1]=3;
    jacStampSpecial[3][2]=4;
    jacStampSpecial[3][3]=5;
    currentCol = 4;
    if( nqsMod )
    {
      jacStampSpecial[3][currentCol] = qCol;
      ++currentCol;
    }
    if( icVBSGiven )
    {
      jacStampSpecial[3][currentCol] = icVBSCol;
      ++currentCol;
    }

    // row KCL d'
    numNonZeros = 5 + (nqsMod ? 1:0 );
    jacStampSpecial[4].resize( numNonZeros );
    jacStampSpecial[4][0]=0;
    jacStampSpecial[4][1]=1;
    jacStampSpecial[4][2]=3;
    jacStampSpecial[4][3]=4;
    jacStampSpecial[4][4]=5;
    currentCol = 5;
    if( nqsMod )
    {
      jacStampSpecial[4][currentCol] = qCol;
      ++currentCol;
    }

    // row KCL s'
    numNonZeros = 5 + (nqsMod ? 1:0 );
    jacStampSpecial[5].resize( numNonZeros );
    jacStampSpecial[5][0]=1;
    jacStampSpecial[5][1]=2;
    jacStampSpecial[5][2]=3;
    jacStampSpecial[5][3]=4;
    jacStampSpecial[5][4]=5;
    currentCol = 5;
    if( nqsMod )
    {
      jacStampSpecial[5][currentCol] = qCol;
      ++currentCol;
    }

    int currentRow = 6;

    // Q row if we need it
    if( nqsMod )
    {
      // add in charge row
      jacStampSpecial[currentRow].resize( 5 );
      jacStampSpecial[currentRow][0] = 1;
      jacStampSpecial[currentRow][1] = 3;
      jacStampSpecial[currentRow][2] = 4;
      jacStampSpecial[currentRow][3] = 5;
      jacStampSpecial[currentRow][4] = qCol;
      ++currentRow;
    }

    // icVBS row if we need it
    if( icVBSGiven )
    {
      jacStampSpecial[currentRow].resize( 3 );
      jacStampSpecial[currentRow][0] = 2;
      jacStampSpecial[currentRow][1] = 3;
      jacStampSpecial[currentRow][2] = icVBSCol;
      ++currentRow;
    }

    // icVDS row if we need it
    if( icVDSGiven )
    {
      jacStampSpecial[currentRow].resize( 3 );
      jacStampSpecial[currentRow][0] = 0;
      jacStampSpecial[currentRow][1] = 2;
      jacStampSpecial[currentRow][2] = icVDSCol;
      ++currentRow;
    }

    // icVGS row if we need it
    if( icVGSGiven )
    {
      jacStampSpecial[currentRow].resize( 3 );
      jacStampSpecial[currentRow][0] = 1;
      jacStampSpecial[currentRow][1] = 2;
      jacStampSpecial[currentRow][2] = icVGSCol;
      ++currentRow;
    }

    // now we just need to merge Vd' and or Vs' nodes if that's called
    // for in this device
    if ( (drainConductance == 0.0) && (sourceConductance == 0.0) )
    {
      // temporary to hold intermediate results
      vector< vector<int> > jacStampSpecialMergedTemp;
      vector<int>           jacSpecialMapTemp;
      vector< vector<int> > jacSpecialMapTemp2;

      jacStampMap( jacStampSpecial,           jacSpecialMap,     jacSpecialMap2,
               jacStampSpecialMergedTemp, jacSpecialMapTemp, jacSpecialMapTemp2,
               5, 2, jacStampSpecial.size() );

      jacStampMap( jacStampSpecialMergedTemp, jacSpecialMapTemp,   jacSpecialMapTemp2,
               jacStampSpecialMerged,     jacSpecialMergedMap, jacSpecialMergedMap2,
               4, 0, jacStampSpecial.size() );

    }
    else if (drainConductance == 0.0)
    {
      jacStampMap( jacStampSpecial,          jacSpecialMap,          jacSpecialMap2,
               jacStampSpecialMerged,    jacSpecialMergedMap,    jacSpecialMergedMap2,
               4, 0, jacStampSpecial.size() );

    }
    else if (sourceConductance == 0.0)
    {
      jacStampMap( jacStampSpecial,       jacSpecialMap,       jacSpecialMap2,
               jacStampSpecialMerged, jacSpecialMergedMap, jacSpecialMergedMap2,
               5, 2, jacStampSpecial.size() );
    }
    else
    {
      // no rows or columns were merged, but we need to initialize
      // jacSpecialMap and jacSpecialMap2 as these will be used to
      // index into the jacobian in registerJacLIDs()
      // copied from DeviceInstance::jacStampMap initialization

      if (jacSpecialMap.size() == 0)
      {
        jacSpecialMap.resize(jacStampSpecial.size());
        jacSpecialMap2.resize(jacStampSpecial.size());
        for (int i=0 ; i<jacStampSpecial.size() ; ++i)
        {
          jacSpecialMap[i] = i;
          jacSpecialMap2[i].resize(jacStampSpecial[i].size());
          for (int j=0 ; j<jacStampSpecial[i].size() ; ++j)
          {
            jacSpecialMap2[i][j] = j;
          }
        }
      }
    }
  }
#ifdef  Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "Instance::Instance jacStampSpecial" << endl;
    for (int k = 0; k<jacSpecialMap.size(); ++k )
    {
      cout << "jacSpecialMap[" << jacSpecialMap[k] << " ] = { ";
      for (int q = 0; q < jacSpecialMap2[k].size(); ++q )
      {
        cout << jacStampSpecial[jacSpecialMap[k]][jacSpecialMap2[k][q]] <<"  ";
      }
      cout << "}" << endl;
    }

    cout << "Instance::Instance jacStampSpecialMerged" << endl;
    for (int k = 0; k<jacSpecialMergedMap.size(); ++k )
    {
      cout << "jacSpecialMap[" << jacSpecialMergedMap[k] << " ] = { ";
      for (int q = 0; q < jacSpecialMergedMap2[k].size(); ++q )
      {
        cout << jacStampSpecialMerged[jacSpecialMergedMap[k]][jacSpecialMergedMap2[k][q]] <<"  ";
      }
      cout << "}" << endl;
    }
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
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
  }
#endif

  numIntVars = 0;
  if ( sourceConductance!=0.0 ) ++numIntVars;
  if ( drainConductance!=0.0 ) ++numIntVars;
  if ( nqsMod ) ++numIntVars;
  if (icVBSGiven) ++numIntVars;
  if (icVDSGiven) ++numIntVars;
  if (icVGSGiven) ++numIntVars;

  if ( numIntVars !=  numInt )
  {
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

  li_Drain = extLIDVec[0];
  li_Gate = extLIDVec[1];
  li_Source = extLIDVec[2];
  li_Bulk = extLIDVec[3];

  int intLoc = 0;

  if( drainConductance != 0.0 )
    li_DrainPrime = intLIDVec[intLoc++];
  else
    li_DrainPrime = li_Drain;

  if( sourceConductance != 0.0 )
    li_SourcePrime = intLIDVec[intLoc++];
  else
    li_SourcePrime = li_Source;

  if( nqsMod )
    li_Charge = intLIDVec[intLoc++];

  if( icVBSGiven )
  {
    if( li_Bulk == li_Source )
    {
       msg = "Instance::registerLIDs:";
       msg += "Tried to specify an initial condition on V_Bulk_Source ";
       msg += "when Bulk and Source nodes are the same node.";
       N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    li_Ibs = intLIDVec[intLoc++];
  }

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

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "\n  local variable indices:\n";
    cout << "  li_Drain         = " << li_Drain << endl;
    cout << "  li_Gate          = " << li_Gate << endl;
    cout << "  li_Source        = " << li_Source << endl;
    cout << "  li_Bulk          = " << li_Bulk << endl;

    if (drainConductance)
      cout << "  li_DrainPrime  = " << li_DrainPrime << endl;
    if (sourceConductance)
      cout << "  li_SourcePrime = " << li_SourcePrime << endl;

    if (nqsMod)
      cout << "  li_Charge      = " << li_Charge << endl;

    if (icVBSGiven)
      cout << "  li_Ibs         = " << li_Ibs << endl;

    if (icVDSGiven)
      cout << "  li_Ids         = " << li_Ids << endl;

    if (icVGSGiven)
      cout << "  li_Igs         = " << li_Igs << endl;
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
    if ( drainConductance!=0.0 )
    {
      tmpstr = getName()+"_drainprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_DrainPrime ] = tmpstr;
    }

    if ( sourceConductance != 0.0 )
    {
      tmpstr = getName()+"_sourceprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_SourcePrime ] = tmpstr;
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
  }

  return intNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 4/3/2013
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
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
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

  int lid=0;
  // Intrinsic capacitors:
  li_state_qb        = staLIDVec[lid++];
  li_state_qg        = staLIDVec[lid++];
  li_state_qd        = staLIDVec[lid++];

  // Parasitic capacitors:
  li_state_qbs       = staLIDVec[lid++];
  li_state_qbd       = staLIDVec[lid++];

  // state variables, cheq
  li_state_qcheq     = staLIDVec[lid++];

  // state variables, cdump
  li_state_qcdump    = staLIDVec[lid++];

  // state variable, qdef
  li_state_qdef      = staLIDVec[lid++];


#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "  Local State indices:" << endl;
    cout << endl;
    cout << "  li_state_qb            = " << li_state_qb  << endl;
    cout << "  li_state_qg            = " << li_state_qg << endl;
    cout << "  li_state_qd            = " << li_state_qd << endl;
    cout << "  li_state_qbs           = " << li_state_qbs << endl;
    cout << "  li_state_qbd           = " << li_state_qbd << endl;
    cout << "  li_state_qcheq         = " << li_state_qcheq << endl;
    cout << "  li_state_qcdump        = " << li_state_qcdump << endl;
    cout << "  li_state_qdef          = " << li_state_qdef << endl;
    cout << dashedline << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/9/11
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const vector<int> & stoLIDVecRef )
{

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
    string msg = "Instance::registerStoreLIDs:";
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

  int lid=0;
  // Voltage drops:
  li_store_vbd       = stoLIDVec[lid++];
  li_store_vbs       = stoLIDVec[lid++];
  li_store_vgs       = stoLIDVec[lid++];
  li_store_vds       = stoLIDVec[lid++];
  li_store_von       = stoLIDVec[lid++];

  if( loadLeadCurrent )
  {
    li_store_dev_id = stoLIDVec[lid++];
    li_store_dev_ig = stoLIDVec[lid++];
    li_store_dev_is = stoLIDVec[lid++];
    li_store_dev_ib = stoLIDVec[lid++];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    cout << "  Local Store indices:" << endl;
    cout << endl;
    cout << "  li_store_vbd           = " << li_store_vbd << endl;
    cout << "  li_store_vbs           = " << li_store_vbs << endl;
    cout << "  li_store_vgs           = " << li_store_vgs << endl;
    cout << "  li_store_vds           = " << li_store_vds << endl;
    cout << "  li_store_von           = " << li_store_von ;
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
// Creation Date : 9/4/02
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  if (icVBSGiven || icVDSGiven || icVGSGiven || nqsMod )
  {
    if( !drainConductance  || !sourceConductance )
    {
      return jacStampSpecialMerged;
    }
    else
    {
      return jacStampSpecial;
    }
  }
  else
  {
    if( drainConductance && sourceConductance && !nqsMod )
      return jacStamp_DC_SC;
    else if( drainConductance && !sourceConductance && !nqsMod )
      return jacStamp_DC;
    else if( !drainConductance && sourceConductance && !nqsMod )
      return jacStamp_SC;
    else if( !drainConductance && !sourceConductance && !nqsMod )
      return jacStamp;
    else
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
              ": NQSMOD not supported for DIRECT MATRIX ACCESS\n" );
  }

  N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
          "Instance::jacobianStamp should not get here!\n" );

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/4/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  vector<int> map;
  vector< vector<int> > map2;


  if (icVBSGiven || icVDSGiven || icVGSGiven || nqsMod )
  {
    if( drainConductance && sourceConductance )
    {
      map = jacSpecialMap;
      map2 = jacSpecialMap2;
    }
    else
    {
      map = jacSpecialMergedMap;
      map2 = jacSpecialMergedMap2;
    }
  }
  else
  {
    if (drainConductance)
    {
      if (sourceConductance)
      {
        map = jacMap_DC_SC;
        map2 = jacMap2_DC_SC;
      }
      else
      {
        map = jacMap_DC;
        map2 = jacMap2_DC;
      }
    }
    else
    {
      if (sourceConductance)
      {
        map = jacMap_SC;
        map2 = jacMap2_SC;
      }
      else
      {
        map = jacMap;
        map2 = jacMap2;
      }
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
    {
      cout << "Instance::registerJacLIDs map selected" << endl;
      for (int k = 0; k<map.size(); ++k )
        {
          cout << "map[ " << k << "] = " << map[k] << "  map2[] = { ";
          for (int q = 0; q < map2[k].size(); ++q )
            {
              cout << map2[k][q] <<"  ";
            }
          cout << "}" << endl;
        }

      for(int k = 0; k<jacLIDVec.size(); ++k )
        {
          cout << "jacLIDVec[ " << k << "] = { ";
          for (int q = 0; q < jacLIDVec[k].size(); ++q )
            {
              cout << jacLIDVec[k][q] <<"  ";
            }
          cout << "}" << endl;
        }
    }
#endif

  int nextColumn = 0;

  // V_d row
  ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
  ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];
  if( icVDSGiven )
  {
    ADrainEquIdsOffset           = jacLIDVec[map[0]][map2[0][2]];
  }

  // V_g row
  AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
  AGateEquBulkNodeOffset               = jacLIDVec[map[1]][map2[1][1]];
  AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][2]];
  AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][3]];
  nextColumn = 4;
  if( nqsMod )
  {
    AGateEquChargeVarOffset          = jacLIDVec[map[1]][map2[1][nextColumn]];
    ++nextColumn;
  }
  if( icVGSGiven )
  {
    AGateEquIgsOffset            = jacLIDVec[map[1]][map2[1][nextColumn]];
    ++nextColumn;
  }

  // V_s row
  ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
  ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];
  nextColumn = 2;
  if( icVBSGiven )
  {
    ASourceEquIbsOffset          = jacLIDVec[map[2]][map2[2][nextColumn]];
    ++nextColumn;
  }
  if( icVDSGiven )
  {
    ASourceEquIdsOffset          = jacLIDVec[map[2]][map2[2][nextColumn]];
    ++nextColumn;
  }
  if( icVGSGiven )
  {
    ASourceEquIgsOffset          = jacLIDVec[map[2]][map2[2][nextColumn]];
    ++nextColumn;
  }

  // V_b row
  ABulkEquGateNodeOffset               = jacLIDVec[map[3]][map2[3][0]];
  ABulkEquBulkNodeOffset               = jacLIDVec[map[3]][map2[3][1]];
  ABulkEquDrainPrimeNodeOffset         = jacLIDVec[map[3]][map2[3][2]];
  ABulkEquSourcePrimeNodeOffset        = jacLIDVec[map[3]][map2[3][3]];
  nextColumn = 4;
  if( nqsMod )
  {
    ABulkEquChargeVarOffset          = jacLIDVec[map[3]][map2[3][nextColumn]];
    ++nextColumn;
  }
  if( icVBSGiven )
  {
    ABulkEquIbsOffset            = jacLIDVec[map[3]][map2[3][nextColumn]];
    ++nextColumn;
  }

  // V_d'
  ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[4]][map2[4][0]];
  ADrainPrimeEquGateNodeOffset         = jacLIDVec[map[4]][map2[4][1]];
  ADrainPrimeEquBulkNodeOffset         = jacLIDVec[map[4]][map2[4][2]];
  ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[4]][map2[4][3]];
  ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[4]][map2[4][4]];
  if( nqsMod )
  {
    ADrainPrimeEquChargeVarOffset    = jacLIDVec[map[4]][map2[4][5]];
  }


  // V_s'
  ASourcePrimeEquGateNodeOffset        = jacLIDVec[map[5]][map2[5][0]];
  ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[5]][map2[5][1]];
  ASourcePrimeEquBulkNodeOffset        = jacLIDVec[map[5]][map2[5][2]];
  ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[5]][map2[5][3]];
  ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[5]][map2[5][4]];
  if( nqsMod )
  {
    ASourcePrimeEquChargeVarOffset   = jacLIDVec[map[5]][map2[5][5]];
  }

  int nextRow = 6;
  if( nqsMod )
  {
    AChargeEquChargeVarOffset         = jacLIDVec[map[nextRow]][map2[nextRow][0]];
    AChargeEquDrainPrimeNodeOffset    = jacLIDVec[map[nextRow]][map2[nextRow][1]];
    AChargeEquGateNodeOffset          = jacLIDVec[map[nextRow]][map2[nextRow][2]];
    AChargeEquSourcePrimeNodeOffset   = jacLIDVec[map[nextRow]][map2[nextRow][3]];
    AChargeEquBulkNodeOffset          = jacLIDVec[map[nextRow]][map2[nextRow][4]];
    ++nextRow;
  }


  if( icVBSGiven )
  {
    icVBSEquVbOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][0]];
    icVBSEquVsOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][1]];
    icVBSEquIbsOffset                 = jacLIDVec[map[nextRow]][map2[nextRow][2]];
    ++nextRow;
  }

  if( icVDSGiven )
  {
    icVDSEquVdOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][0]];
    icVDSEquVsOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][1]];
    icVDSEquIdsOffset                 = jacLIDVec[map[nextRow]][map2[nextRow][2]];
    ++nextRow;
  }

  if( icVGSGiven )
  {
    icVGSEquVgOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][0]];
    icVGSEquVsOffset                  = jacLIDVec[map[nextRow]][map2[nextRow][1]];
    icVGSEquIgsOffset                 = jacLIDVec[map[nextRow]][map2[nextRow][2]];
    ++nextRow;
  }

  if (nqsMod)
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
     ": NQSMOD not supported.\n" );
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  // V_d row
  f_DrainEquDrainNodePtr             = 	  &(dFdx[li_Drain][ADrainEquDrainNodeOffset]);
  f_DrainEquDrainPrimeNodePtr        = 	  &(dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);
  if( icVDSGiven )
  {
    f_DrainEquIdsPtr           = 	    &(dFdx[li_Drain][ADrainEquIdsOffset]);
  }

  // V_g row
  f_GateEquGateNodePtr               = 	  &(dFdx[li_Gate][AGateEquGateNodeOffset]);
  f_GateEquBulkNodePtr               = 	  &(dFdx[li_Gate][AGateEquBulkNodeOffset]);
  f_GateEquDrainPrimeNodePtr         = 	  &(dFdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  f_GateEquSourcePrimeNodePtr        = 	  &(dFdx[li_Gate][AGateEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    f_GateEquChargeVarPtr          = 	    &(dFdx[li_Gate][AGateEquChargeVarOffset]);

  }
  if( icVGSGiven )
  {
    f_GateEquIgsPtr            = 	    &(dFdx[li_Gate][AGateEquIgsOffset]);

  }

  // V_s row	  // V_s row
  f_SourceEquSourceNodePtr           = 	  &(dFdx[li_Source][ASourceEquSourceNodeOffset]);
  f_SourceEquSourcePrimeNodePtr      = 	  &(dFdx[li_Source][ASourceEquSourcePrimeNodeOffset]);
  if( icVBSGiven )
  {
    f_SourceEquIbsPtr          = 	    &(dFdx[li_Source][ASourceEquIbsOffset]);

  }
  if( icVDSGiven )
  {
    f_SourceEquIdsPtr          = 	    &(dFdx[li_Source][ASourceEquIdsOffset]);

  }
  if( icVGSGiven )
  {
    f_SourceEquIgsPtr          = 	    &(dFdx[li_Source][ASourceEquIgsOffset]);

  }

  // V_b row
  f_BulkEquGateNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquGateNodeOffset]);
  f_BulkEquBulkNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquBulkNodeOffset]);
  f_BulkEquDrainPrimeNodePtr         = 	  &(dFdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]);
  f_BulkEquSourcePrimeNodePtr        = 	  &(dFdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    f_BulkEquChargeVarPtr          = 	    &(dFdx[li_Bulk][ABulkEquChargeVarOffset]);

  }
  if( icVBSGiven )
  {
    f_BulkEquIbsPtr            = 	    &(dFdx[li_Bulk][ABulkEquIbsOffset]);

  }

  // V_d'
  f_DrainPrimeEquDrainNodePtr        = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  f_DrainPrimeEquGateNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  f_DrainPrimeEquBulkNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]);
  f_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  f_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    f_DrainPrimeEquChargeVarPtr    = 	    &(dFdx[li_DrainPrime][ADrainPrimeEquChargeVarOffset]);
  }


  // V_s'
  f_SourcePrimeEquGateNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  f_SourcePrimeEquSourceNodePtr      = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  f_SourcePrimeEquBulkNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]);
  f_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  f_SourcePrimeEquSourcePrimeNodePtr = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    f_SourcePrimeEquChargeVarPtr   = 	    &(dFdx[li_SourcePrime][ASourcePrimeEquChargeVarOffset]);
  }

  if( nqsMod )
  {
    f_ChargeEquChargeVarPtr         = 	    &(dFdx[li_SourcePrime][AChargeEquChargeVarOffset]);
    f_ChargeEquDrainPrimeNodePtr    = 	    &(dFdx[li_SourcePrime][AChargeEquDrainPrimeNodeOffset]);
    f_ChargeEquGateNodePtr          = 	    &(dFdx[li_SourcePrime][AChargeEquGateNodeOffset]);
    f_ChargeEquSourcePrimeNodePtr   = 	    &(dFdx[li_SourcePrime][AChargeEquSourcePrimeNodeOffset]);
    f_ChargeEquBulkNodePtr          = 	    &(dFdx[li_SourcePrime][AChargeEquBulkNodeOffset]);
  }


  if( icVBSGiven )
  {
    f_icVBSEquVbPtr                  = 	    &(dFdx[li_Ibs][icVBSEquVbOffset]);
    f_icVBSEquVsPtr                  = 	    &(dFdx[li_Ibs][icVBSEquVsOffset]);
    f_icVBSEquIbsPtr                 = 	    &(dFdx[li_Ibs][icVBSEquIbsOffset]);
  }

  if( icVDSGiven )
  {
    f_icVDSEquVdPtr                  = 	    &(dFdx[li_Ids][icVDSEquVdOffset]);
    f_icVDSEquVsPtr                  = 	    &(dFdx[li_Ids][icVDSEquVsOffset]);
    f_icVDSEquIdsPtr                 = 	    &(dFdx[li_Ids][icVDSEquIdsOffset]);
  }

  if( icVGSGiven )
  {
    f_icVGSEquVgPtr                  = 	    &(dFdx[li_Igs][icVGSEquVgOffset]);
    f_icVGSEquVsPtr                  = 	    &(dFdx[li_Igs][icVGSEquVsOffset]);
    f_icVGSEquIgsPtr                 = 	    &(dFdx[li_Igs][icVGSEquIgsOffset]);
  }



  // V_d row
  q_DrainEquDrainNodePtr             = 	  &(dQdx[li_Drain][ADrainEquDrainNodeOffset]);
  q_DrainEquDrainPrimeNodePtr        = 	  &(dQdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);
  if( icVDSGiven )
  {
    q_DrainEquIdsPtr           = 	    &(dQdx[li_Drain][ADrainEquIdsOffset]);
  }

  // V_g row
  q_GateEquGateNodePtr               = 	  &(dQdx[li_Gate][AGateEquGateNodeOffset]);
  q_GateEquBulkNodePtr               = 	  &(dQdx[li_Gate][AGateEquBulkNodeOffset]);
  q_GateEquDrainPrimeNodePtr         = 	  &(dQdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  q_GateEquSourcePrimeNodePtr        = 	  &(dQdx[li_Gate][AGateEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    q_GateEquChargeVarPtr          = 	    &(dQdx[li_Gate][AGateEquChargeVarOffset]);

  }
  if( icVGSGiven )
  {
    q_GateEquIgsPtr            = 	    &(dQdx[li_Gate][AGateEquIgsOffset]);

  }

  // V_s row	  // V_s row
  q_SourceEquSourceNodePtr           = 	  &(dQdx[li_Source][ASourceEquSourceNodeOffset]);
  q_SourceEquSourcePrimeNodePtr      = 	  &(dQdx[li_Source][ASourceEquSourcePrimeNodeOffset]);
  if( icVBSGiven )
  {
    q_SourceEquIbsPtr          = 	    &(dQdx[li_Source][ASourceEquIbsOffset]);

  }
  if( icVDSGiven )
  {
    q_SourceEquIdsPtr          = 	    &(dQdx[li_Source][ASourceEquIdsOffset]);

  }
  if( icVGSGiven )
  {
    q_SourceEquIgsPtr          = 	    &(dQdx[li_Source][ASourceEquIgsOffset]);

  }

  // V_b row
  q_BulkEquGateNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquGateNodeOffset]);
  q_BulkEquBulkNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquBulkNodeOffset]);
  q_BulkEquDrainPrimeNodePtr         = 	  &(dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]);
  q_BulkEquSourcePrimeNodePtr        = 	  &(dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    q_BulkEquChargeVarPtr          = 	    &(dQdx[li_Bulk][ABulkEquChargeVarOffset]);

  }
  if( icVBSGiven )
  {
    q_BulkEquIbsPtr            = 	    &(dQdx[li_Bulk][ABulkEquIbsOffset]);

  }

  // V_d'
  q_DrainPrimeEquDrainNodePtr        = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  q_DrainPrimeEquGateNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  q_DrainPrimeEquBulkNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]);
  q_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  q_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    q_DrainPrimeEquChargeVarPtr    = 	    &(dQdx[li_DrainPrime][ADrainPrimeEquChargeVarOffset]);
  }


  // V_s'
  q_SourcePrimeEquGateNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  q_SourcePrimeEquSourceNodePtr      = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  q_SourcePrimeEquBulkNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]);
  q_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  q_SourcePrimeEquSourcePrimeNodePtr = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);
  if( nqsMod )
  {
    q_SourcePrimeEquChargeVarPtr   = 	    &(dQdx[li_SourcePrime][ASourcePrimeEquChargeVarOffset]);
  }

  if( nqsMod )
  {
    q_ChargeEquChargeVarPtr         = 	    &(dQdx[li_SourcePrime][AChargeEquChargeVarOffset]);
    q_ChargeEquDrainPrimeNodePtr    = 	    &(dQdx[li_SourcePrime][AChargeEquDrainPrimeNodeOffset]);
    q_ChargeEquGateNodePtr          = 	    &(dQdx[li_SourcePrime][AChargeEquGateNodeOffset]);
    q_ChargeEquSourcePrimeNodePtr   = 	    &(dQdx[li_SourcePrime][AChargeEquSourcePrimeNodeOffset]);
    q_ChargeEquBulkNodePtr          = 	    &(dQdx[li_SourcePrime][AChargeEquBulkNodeOffset]);
  }


  if( icVBSGiven )
  {
    q_icVBSEquVbPtr                  = 	    &(dQdx[li_Ibs][icVBSEquVbOffset]);
    q_icVBSEquVsPtr                  = 	    &(dQdx[li_Ibs][icVBSEquVsOffset]);
    q_icVBSEquIbsPtr                 = 	    &(dQdx[li_Ibs][icVBSEquIbsOffset]);
  }

  if( icVDSGiven )
  {
    q_icVDSEquVdPtr                  = 	    &(dQdx[li_Ids][icVDSEquVdOffset]);
    q_icVDSEquVsPtr                  = 	    &(dQdx[li_Ids][icVDSEquVsOffset]);
    q_icVDSEquIdsPtr                 = 	    &(dQdx[li_Ids][icVDSEquIdsOffset]);
  }

  if( icVGSGiven )
  {
    q_icVGSEquVgPtr                  = 	    &(dQdx[li_Igs][icVGSEquVgOffset]);
    q_icVGSEquVsPtr                  = 	    &(dQdx[li_Igs][icVGSEquVsOffset]);
    q_icVGSEquIgsPtr                 = 	    &(dQdx[li_Igs][icVGSEquIgsOffset]);
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/22/00
//-----------------------------------------------------------------------------
bool Instance::updateTemperature (const double & temp_tmp)
{
  char msg[128];

  double tmp, tmp1, tmp2, tmp3, Eg;
  double T0, T1, T2, T3, T4, T5, Ldrn, Wdrn;
  double delTemp, TRatio, Inv_L, Inv_W, Inv_LW;
  //double Dw, Dl;
  double Tnom;
  double Nvtm, SourceSatCurrent, DrainSatCurrent;

  // stuff from model paramters:
  double Eg0 = model_.Eg0;
  double ni = model_.ni;
  double Vtm0 = model_.Vtm0;

  bool bsuccess = true;

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

  Tnom = model_.tnom;
  TRatio = temp/Tnom;

  vtm = CONSTKoverQ * temp;
  Eg = CONSTEg0 - CONSTalphaEg * temp * temp / (temp + CONSTbetaEg);

  if (temp != Tnom)
  {
    T0 = Eg0 / Vtm0 - Eg/vtm + model_.jctTempExponent*log(temp/Tnom);
    T1 = exp(T0 / model_.jctEmissionCoeff);
    jctTempSatCurDensity         = model_.jctSatCurDensity* T1;
    jctSidewallTempSatCurDensity = model_.jctSidewallSatCurDensity * T1;
  }
  else
  {
    jctTempSatCurDensity         = model_.jctSatCurDensity;
    jctSidewallTempSatCurDensity = model_.jctSidewallSatCurDensity;
  }

  if (jctTempSatCurDensity < 0.0)         jctTempSatCurDensity = 0.0;
  if (jctSidewallTempSatCurDensity < 0.0) jctSidewallTempSatCurDensity = 0.0;


  // Temperature dependence of D/B and S/B diode capacitance begins
  delTemp = temp - Tnom;
  T0 = model_.tcj * delTemp;

  if (T0 >= -1.0)
  {
    unitAreaJctCapTemp = model_.unitAreaJctCap *(1.0 + T0);
  }
  else if (unitAreaJctCapTemp > 0.0)
  {
    unitAreaJctCapTemp = 0.0;

    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, string(
   "Temperature effect has caused cj to be negative. Cj clamped to zero.\n"));
  }

  T0 = model_.tcjsw * delTemp;

  if (T0 >= -1.0)
  {
    unitLengthSidewallJctCapTemp =
      model_.unitLengthSidewallJctCap *(1.0 + T0);
  }
  else if (unitLengthSidewallJctCapTemp > 0.0)
  {
    unitLengthSidewallJctCapTemp = 0.0;
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, string(
"Temperature effect has caused cjsw to be negative. Cjsw clamped to zero.\n"
    ));
  }

  T0 = model_.tcjswg * delTemp;

  if (T0 >= -1.0)
  {
    unitLengthGateSidewallJctCapTemp =
      model_.unitLengthGateSidewallJctCap *(1.0 + T0);
  }
  else if (unitLengthGateSidewallJctCapTemp > 0.0)
  {
    unitLengthGateSidewallJctCapTemp = 0.0;
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, string(
"Temperature effect has caused cjswg to be negative. Cjswg clamped to zero.\n"
   ));
  }

  PhiBTemp = model_.bulkJctPotential - model_.tpb * delTemp;

  if (PhiBTemp < 0.01)
  {
    PhiBTemp = 0.01;
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, string(
   "Temperature effect has caused pb to be < 0.01. Pb clamped to 0.01.\n"
    ));
  }

  PhiBSWTemp = model_.sidewallJctPotential - model_.tpbsw * delTemp;

  if (PhiBSWTemp <= 0.01)
  {
    PhiBSWTemp = 0.01;
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, string(
    "Temperature effect has caused pbsw to be < 0.01. Pbsw clamped to 0.01.\n"
    ));
  }

  PhiBSWGTemp = model_.GatesidewallJctPotential - model_.tpbswg * delTemp;

  if (PhiBSWGTemp <= 0.01)
  {
    PhiBSWGTemp = 0.01;
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, string(
   "Temperature effect has caused pbswg to be < 0.01. Pbswg clamped to 0.01.\n"
    ));
  }
  // End of junction capacitance


  // This next block determines whether  or not to use a previously allocated
  // set of size dependent parameters. These are stored in a list that is
  // owned by the model.  If the values for length and width match those of
  // a previously allocated set, then use the old set.  If not, allocate a new set.

  list<SizeDependParam*>::iterator it_dpL =
    model_.sizeDependParamList.begin();
  list<SizeDependParam*>::iterator end_dpL =
    model_.sizeDependParamList.end();

  paramPtr = NULL;

  for( ; it_dpL != end_dpL; ++it_dpL )
    if( (*it_dpL)->Length == l && (*it_dpL)->Width == w )
      paramPtr = (*it_dpL);

  if ( paramPtr != NULL )
  {
  }
  else
  {
    paramPtr = new SizeDependParam ();

    model_.sizeDependParamList.push_back( paramPtr );
    paramPtr->referenceTemperature = temp_tmp;

    Ldrn = l;
    Wdrn = w;
    paramPtr->Length = Ldrn;
    paramPtr->Width = Wdrn;

    T0 = pow(Ldrn, model_.Lln);
    T1 = pow(Wdrn, model_.Lwn);

    tmp1 = model_.Ll / T0 + model_.Lw / T1
      + model_.Lwl / (T0 * T1);

    paramPtr->dl = model_.Lint + tmp1;

    tmp2 = model_.Llc / T0 + model_.Lwc / T1
      + model_.Lwlc / (T0 * T1);

    paramPtr->dlc = model_.dlc + tmp2;

    T2 = pow(Ldrn, model_.Wln);
    T3 = pow(Wdrn, model_.Wwn);

    tmp1 = model_.Wl / T2 + model_.Ww / T3
      + model_.Wwl / (T2 * T3);

    paramPtr->dw = model_.Wint + tmp1;
    tmp2 = model_.Wlc / T2 + model_.Wwc / T3
      + model_.Wwlc / (T2 * T3);

    paramPtr->dwc = model_.dwc + tmp2;

    paramPtr->leff = l - 2.0 * paramPtr->dl;
    if (paramPtr->leff <= 0.0)
    {
      sprintf(msg,
              ": mosfet %s, model %s: Effective channel length <= 0",
              getName().c_str(), model_.getName().c_str());
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
    }

    paramPtr->weff = w - 2.0 * paramPtr->dw;
    if (paramPtr->weff <= 0.0)
    {
      sprintf(msg,
              ": mosfet %s, model %s: Effective channel width <= 0",
              getName().c_str(), model_.getName().c_str());
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
    }

    paramPtr->leffCV = l - 2.0 * paramPtr->dlc;
    if (paramPtr->leffCV <= 0.0)
    {
      sprintf(msg,
              ": mosfet %s, model %s: Effective channel length for C-V <= 0",
              getName().c_str(), model_.getName().c_str());
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
    }

    paramPtr->weffCV = w - 2.0 * paramPtr->dwc;
    if (paramPtr->weffCV <= 0.0)
    {
      sprintf(msg,
              ": mosfet %s, model %s: Effective channel width for C-V <= 0",
              getName().c_str(), model_.getName().c_str());
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
    }


    if (model_.binUnit == 1)
    {
      Inv_L = 1.0e-6 / paramPtr->leff;
      Inv_W = 1.0e-6 / paramPtr->weff;
      Inv_LW = 1.0e-12 / (paramPtr->leff * paramPtr->weff);
    }
    else
    {
      Inv_L = 1.0 / paramPtr->leff;
      Inv_W = 1.0 / paramPtr->weff;
      Inv_LW = 1.0 / (paramPtr->leff * paramPtr->weff);
    }

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

    paramPtr->cit = model_.cit
      + model_.lcit * Inv_L
      + model_.wcit * Inv_W
      + model_.pcit * Inv_LW;

    paramPtr->nfactor = model_.nfactor
      + model_.lnfactor * Inv_L
      + model_.wnfactor * Inv_W
      + model_.pnfactor * Inv_LW;

    paramPtr->xj = model_.xj
      + model_.lxj * Inv_L
      + model_.wxj * Inv_W
      + model_.pxj * Inv_LW;

    paramPtr->vsat = model_.vsat
      + model_.lvsat * Inv_L
      + model_.wvsat * Inv_W
      + model_.pvsat * Inv_LW;

    paramPtr->at = model_.at
      + model_.lat * Inv_L
      + model_.wat * Inv_W
      + model_.pat * Inv_LW;

    paramPtr->a0 = model_.a0
      + model_.la0 * Inv_L
      + model_.wa0 * Inv_W
      + model_.pa0 * Inv_LW;

    paramPtr->ags = model_.ags
      + model_.lags * Inv_L
      + model_.wags * Inv_W
      + model_.pags * Inv_LW;

    paramPtr->a1 = model_.a1
      + model_.la1 * Inv_L
      + model_.wa1 * Inv_W
      + model_.pa1 * Inv_LW;

    paramPtr->a2 = model_.a2
      + model_.la2 * Inv_L
      + model_.wa2 * Inv_W
      + model_.pa2 * Inv_LW;

    paramPtr->keta = model_.keta
      + model_.lketa * Inv_L
      + model_.wketa * Inv_W
      + model_.pketa * Inv_LW;

    paramPtr->nsub = model_.nsub
      + model_.lnsub * Inv_L
      + model_.wnsub * Inv_W
      + model_.pnsub * Inv_LW;

    paramPtr->npeak = model_.npeak
      + model_.lnpeak * Inv_L
      + model_.wnpeak * Inv_W
      + model_.pnpeak * Inv_LW;

    paramPtr->ngate = model_.ngate
      + model_.lngate * Inv_L
      + model_.wngate * Inv_W
      + model_.pngate * Inv_LW;

    paramPtr->gamma1 = model_.gamma1
      + model_.lgamma1 * Inv_L
      + model_.wgamma1 * Inv_W
      + model_.pgamma1 * Inv_LW;

    paramPtr->gamma2 = model_.gamma2
      + model_.lgamma2 * Inv_L
      + model_.wgamma2 * Inv_W
      + model_.pgamma2 * Inv_LW;

    paramPtr->vbx = model_.vbx
      + model_.lvbx * Inv_L
      + model_.wvbx * Inv_W
      + model_.pvbx * Inv_LW;

    paramPtr->vbm = model_.vbm
      + model_.lvbm * Inv_L
      + model_.wvbm * Inv_W
      + model_.pvbm * Inv_LW;

    paramPtr->xt = model_.xt
      + model_.lxt * Inv_L
      + model_.wxt * Inv_W
      + model_.pxt * Inv_LW;

    paramPtr->vfb = model_.vfb
      + model_.lvfb * Inv_L
      + model_.wvfb * Inv_W
      + model_.pvfb * Inv_LW;

    paramPtr->k1 = model_.k1
      + model_.lk1 * Inv_L
      + model_.wk1 * Inv_W
      + model_.pk1 * Inv_LW;

    paramPtr->kt1 = model_.kt1
      + model_.lkt1 * Inv_L
      + model_.wkt1 * Inv_W
      + model_.pkt1 * Inv_LW;

    paramPtr->kt1l = model_.kt1l
      + model_.lkt1l * Inv_L
      + model_.wkt1l * Inv_W
      + model_.pkt1l * Inv_LW;

    paramPtr->k2 = model_.k2
      + model_.lk2 * Inv_L
      + model_.wk2 * Inv_W
      + model_.pk2 * Inv_LW;

    paramPtr->kt2 = model_.kt2
      + model_.lkt2 * Inv_L
      + model_.wkt2 * Inv_W
      + model_.pkt2 * Inv_LW;

    paramPtr->k3 = model_.k3
      + model_.lk3 * Inv_L
      + model_.wk3 * Inv_W
      + model_.pk3 * Inv_LW;

    paramPtr->k3b = model_.k3b
      + model_.lk3b * Inv_L
      + model_.wk3b * Inv_W
      + model_.pk3b * Inv_LW;

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

    paramPtr->drout = model_.drout
      + model_.ldrout * Inv_L
      + model_.wdrout * Inv_W
      + model_.pdrout * Inv_LW;

    paramPtr->dsub = model_.dsub
      + model_.ldsub * Inv_L
      + model_.wdsub * Inv_W
      + model_.pdsub * Inv_LW;

    paramPtr->vth0 = model_.vth0
      + model_.lvth0 * Inv_L
      + model_.wvth0 * Inv_W
      + model_.pvth0 * Inv_LW;

    paramPtr->ua = model_.ua
      + model_.lua * Inv_L
      + model_.wua * Inv_W
      + model_.pua * Inv_LW;

    paramPtr->ua1 = model_.ua1
      + model_.lua1 * Inv_L
      + model_.wua1 * Inv_W
      + model_.pua1 * Inv_LW;

    paramPtr->ub = model_.ub
      + model_.lub * Inv_L
      + model_.wub * Inv_W
      + model_.pub * Inv_LW;

    paramPtr->ub1 = model_.ub1
      + model_.lub1 * Inv_L
      + model_.wub1 * Inv_W
      + model_.pub1 * Inv_LW;

    paramPtr->uc = model_.uc
      + model_.luc * Inv_L
      + model_.wuc * Inv_W
      + model_.puc * Inv_LW;

    paramPtr->uc1 = model_.uc1
      + model_.luc1 * Inv_L
      + model_.wuc1 * Inv_W
      + model_.puc1 * Inv_LW;

    paramPtr->u0 = model_.u0
      + model_.lu0 * Inv_L
      + model_.wu0 * Inv_W
      + model_.pu0 * Inv_LW;

    paramPtr->ute = model_.ute
      + model_.lute * Inv_L
      + model_.wute * Inv_W
      + model_.pute * Inv_LW;

    paramPtr->voff = model_.voff
      + model_.lvoff * Inv_L
      + model_.wvoff * Inv_W
      + model_.pvoff * Inv_LW;

    paramPtr->delta = model_.delta
      + model_.ldelta * Inv_L
      + model_.wdelta * Inv_W
      + model_.pdelta * Inv_LW;

    paramPtr->rdsw = model_.rdsw
      + model_.lrdsw * Inv_L
      + model_.wrdsw * Inv_W
      + model_.prdsw * Inv_LW;

    paramPtr->prwg = model_.prwg
      + model_.lprwg * Inv_L
      + model_.wprwg * Inv_W
      + model_.pprwg * Inv_LW;

    paramPtr->prwb = model_.prwb
      + model_.lprwb * Inv_L
      + model_.wprwb * Inv_W
      + model_.pprwb * Inv_LW;

    paramPtr->prt = model_.prt
      + model_.lprt * Inv_L
      + model_.wprt * Inv_W
      + model_.pprt * Inv_LW;

    paramPtr->eta0 = model_.eta0
      + model_.leta0 * Inv_L
      + model_.weta0 * Inv_W
      + model_.peta0 * Inv_LW;

    paramPtr->etab = model_.etab
      + model_.letab * Inv_L
      + model_.wetab * Inv_W
      + model_.petab * Inv_LW;

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

    paramPtr->pscbe1 = model_.pscbe1
      + model_.lpscbe1 * Inv_L
      + model_.wpscbe1 * Inv_W
      + model_.ppscbe1 * Inv_LW;

    paramPtr->pscbe2 = model_.pscbe2
      + model_.lpscbe2 * Inv_L
      + model_.wpscbe2 * Inv_W
      + model_.ppscbe2 * Inv_LW;

    paramPtr->pvag = model_.pvag
      + model_.lpvag * Inv_L
      + model_.wpvag * Inv_W
      + model_.ppvag * Inv_LW;

    paramPtr->wr = model_.wr
      + model_.lwr * Inv_L
      + model_.wwr * Inv_W
      + model_.pwr * Inv_LW;

    paramPtr->dwg = model_.dwg
      + model_.ldwg * Inv_L
      + model_.wdwg * Inv_W
      + model_.pdwg * Inv_LW;

    paramPtr->dwb = model_.dwb
      + model_.ldwb * Inv_L
      + model_.wdwb * Inv_W
      + model_.pdwb * Inv_LW;

    paramPtr->b0 = model_.b0
      + model_.lb0 * Inv_L
      + model_.wb0 * Inv_W
      + model_.pb0 * Inv_LW;

    paramPtr->b1 = model_.b1
      + model_.lb1 * Inv_L
      + model_.wb1 * Inv_W
      + model_.pb1 * Inv_LW;

    paramPtr->alpha0 = model_.alpha0
      + model_.lalpha0 * Inv_L
      + model_.walpha0 * Inv_W
      + model_.palpha0 * Inv_LW;

    paramPtr->alpha1 = model_.alpha1
      + model_.lalpha1 * Inv_L
      + model_.walpha1 * Inv_W
      + model_.palpha1 * Inv_LW;

    paramPtr->beta0 = model_.beta0
      + model_.lbeta0 * Inv_L
      + model_.wbeta0 * Inv_W
      + model_.pbeta0 * Inv_LW;

    // CV model
    paramPtr->elm = model_.elm
      + model_.lelm * Inv_L
      + model_.welm * Inv_W
      + model_.pelm * Inv_LW;

    paramPtr->cgsl = model_.cgsl
      + model_.lcgsl * Inv_L
      + model_.wcgsl * Inv_W
      + model_.pcgsl * Inv_LW;

    paramPtr->cgdl = model_.cgdl
      + model_.lcgdl * Inv_L
      + model_.wcgdl * Inv_W
      + model_.pcgdl * Inv_LW;

    paramPtr->ckappa = model_.ckappa
      + model_.lckappa * Inv_L
      + model_.wckappa * Inv_W
      + model_.pckappa * Inv_LW;

    paramPtr->cf = model_.cf
      + model_.lcf * Inv_L
      + model_.wcf * Inv_W
      + model_.pcf * Inv_LW;

    paramPtr->clc = model_.clc
      + model_.lclc * Inv_L
      + model_.wclc * Inv_W
      + model_.pclc * Inv_LW;

    paramPtr->cle = model_.cle
      + model_.lcle * Inv_L
      + model_.wcle * Inv_W
      + model_.pcle * Inv_LW;

    paramPtr->vfbcv = model_.vfbcv
      + model_.lvfbcv * Inv_L
      + model_.wvfbcv * Inv_W
      + model_.pvfbcv * Inv_LW;

    paramPtr->acde = model_.acde
      + model_.lacde * Inv_L
      + model_.wacde * Inv_W
      + model_.pacde * Inv_LW;

    paramPtr->moin = model_.moin
      + model_.lmoin * Inv_L
      + model_.wmoin * Inv_W
      + model_.pmoin * Inv_LW;

    paramPtr->noff = model_.noff
      + model_.lnoff * Inv_L
      + model_.wnoff * Inv_W
      + model_.pnoff * Inv_LW;

    paramPtr->voffcv = model_.voffcv
      + model_.lvoffcv * Inv_L
      + model_.wvoffcv * Inv_W
      + model_.pvoffcv * Inv_LW;

    paramPtr->abulkCVfactor = 1.0
      + pow((paramPtr->clc / paramPtr->leffCV), paramPtr->cle);

    T0 = (TRatio - 1.0);

    paramPtr->ua = paramPtr->ua + paramPtr->ua1 * T0;
    paramPtr->ub = paramPtr->ub + paramPtr->ub1 * T0;
    paramPtr->uc = paramPtr->uc + paramPtr->uc1 * T0;

    if (paramPtr->u0 > 1.0) paramPtr->u0 = paramPtr->u0 / 1.0e4;

    paramPtr->u0temp = paramPtr->u0 * pow(TRatio, paramPtr->ute);

    paramPtr->vsattemp = paramPtr->vsat - paramPtr->at * T0;

    paramPtr->rds0 = (paramPtr->rdsw + paramPtr->prt * T0)
      / pow(paramPtr->weff * 1E6, paramPtr->wr);

#ifdef CHECK_MODEL_DONE
    if (checkModel((*M_iter), iterI, ckt))
    {
      sprintf(msg,
  "Fatal error(s) detected during V3.2 parameter checking for %s in model %s",
        name.c_str(),
       model_.name.c_str());
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_, lineNumber_);
    }
#endif

    paramPtr->cgdo = (model_.cgdo + paramPtr->cf) * paramPtr->weffCV;
    paramPtr->cgso = (model_.cgso + paramPtr->cf) * paramPtr->weffCV;
    paramPtr->cgbo = model_.cgbo * paramPtr->leffCV;

    T0 = paramPtr->leffCV * paramPtr->leffCV;

    paramPtr->tconst = paramPtr->u0temp * paramPtr->elm / (model_.cox
                                                           * paramPtr->weffCV * paramPtr->leffCV * T0);

    if (!model_.npeakGiven && model_.gamma1Given)
    {
      T0 = paramPtr->gamma1 * model_.cox;
      paramPtr->npeak = 3.021E22 * T0 * T0;
    }

    paramPtr->phi     = 2.0 * Vtm0 * log(paramPtr->npeak / ni);
    paramPtr->sqrtPhi = sqrt(paramPtr->phi);
    paramPtr->phis3   = paramPtr->sqrtPhi * paramPtr->phi;

    paramPtr->Xdep0 = sqrt(2.0 * CONSTEPSSI / (CONSTQ * paramPtr->npeak * 1.0e6))
      * paramPtr->sqrtPhi;

    paramPtr->sqrtXdep0 = sqrt(paramPtr->Xdep0);
    paramPtr->litl = sqrt(3.0 * paramPtr->xj * model_.tox);

    paramPtr->vbi = Vtm0 * log(1.0e20 * paramPtr->npeak / (ni * ni));

    paramPtr->cdep0 = sqrt(CONSTQ * CONSTEPSSI * paramPtr->npeak * 1.0e6 / 2.0
                           / paramPtr->phi);

    paramPtr->ldeb = sqrt(CONSTEPSSI * Vtm0 / (CONSTQ
                                               * paramPtr->npeak * 1.0e6)) / 3.0;

    paramPtr->acde *= pow((paramPtr->npeak / 2.0e16), -0.25);


    if (model_.k1Given || model_.k2Given)
    {
      if (!model_.k1Given)
      {
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,string("Warning: k1 should be specified with k2.\n"));
        paramPtr->k1 = 0.53;
      }

      if (!model_.k2Given)
      {
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,string("Warning: k2 should be specified with k1.\n"));
        paramPtr->k2 = -0.0186;
      }

      if (model_.nsubGiven)
      {
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,string("Warning: nsub is ignored because k1 or k2 is given.\n"));
      }

      if (model_.xtGiven)
      {
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,string("Warning: xt is ignored because k1 or k2 is given.\n"));
      }

      if (model_.vbxGiven)
      {
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,string("Warning: vbx is ignored because k1 or k2 is given.\n"));
      }

      if (model_.gamma1Given)
      {
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,string("Warning: gamma1 is ignored because k1 or k2 is given.\n"));
      }

      if (model_.gamma2Given)
      {
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,string("Warning: gamma2 is ignored because k1 or k2 is given.\n"));
      }
    }
    else
    {
      if (!model_.vbxGiven)
        paramPtr->vbx = paramPtr->phi - 7.7348e-4 * paramPtr->npeak
          * paramPtr->xt * paramPtr->xt;

      if (paramPtr->vbx > 0.0)
        paramPtr->vbx = -paramPtr->vbx;

      if (paramPtr->vbm > 0.0)
        paramPtr->vbm = -paramPtr->vbm;

      if (!model_.gamma1Given)
        paramPtr->gamma1 = 5.753e-12 * sqrt(paramPtr->npeak) / model_.cox;

      if (!model_.gamma2Given)
        paramPtr->gamma2 = 5.753e-12 * sqrt(paramPtr->nsub) / model_.cox;

      T0 = paramPtr->gamma1 - paramPtr->gamma2;
      T1 = sqrt(paramPtr->phi - paramPtr->vbx) - paramPtr->sqrtPhi;
      T2 = sqrt(paramPtr->phi * (paramPtr->phi - paramPtr->vbm)) - paramPtr->phi;

      paramPtr->k2 = T0 * T1 / (2.0 * T2 + paramPtr->vbm);
      paramPtr->k1 = paramPtr->gamma2 - 2.0 * paramPtr->k2 * sqrt(paramPtr->phi
                                                                  - paramPtr->vbm);
    }

    if (paramPtr->k2 < 0.0)
    {
      T0 = 0.5 * paramPtr->k1 / paramPtr->k2;
      paramPtr->vbsc = 0.9 * (paramPtr->phi - T0 * T0);

      if (paramPtr->vbsc > -3.0) paramPtr->vbsc = -3.0;
      else if (paramPtr->vbsc < -30.0) paramPtr->vbsc = -30.0;
    }
    else
    {
      paramPtr->vbsc = -30.0;
    }

    if (paramPtr->vbsc > paramPtr->vbm) paramPtr->vbsc = paramPtr->vbm;

    if (!model_.vfbGiven)
    {
      if (model_.vth0Given)
      {
        paramPtr->vfb = model_.dtype * paramPtr->vth0
          - paramPtr->phi - paramPtr->k1 * paramPtr->sqrtPhi;
      }
      else
      {   paramPtr->vfb = -1.0;
      }
    }

    if (!model_.vth0Given)
    {
      paramPtr->vth0 = model_.dtype
        * (paramPtr->vfb + paramPtr->phi + paramPtr->k1
           * paramPtr->sqrtPhi);
    }

    paramPtr->k1ox = paramPtr->k1 * model_.tox / model_.toxm;
    paramPtr->k2ox = paramPtr->k2 * model_.tox / model_.toxm;

    T1 = sqrt(CONSTEPSSI / CONSTEPSOX * model_.tox * paramPtr->Xdep0);
    T0 = exp(-0.5 * paramPtr->dsub * paramPtr->leff / T1);

    paramPtr->theta0vb0 = (T0 + 2.0 * T0 * T0);

    T0 = exp(-0.5 * paramPtr->drout * paramPtr->leff / T1);
    T2 = (T0 + 2.0 * T0 * T0);

    paramPtr->thetaRout = paramPtr->pdibl1 * T2 + paramPtr->pdibl2;

    tmp = sqrt(paramPtr->Xdep0);
    tmp1 = paramPtr->vbi - paramPtr->phi;
    tmp2 = model_.factor1 * tmp;

    T0 = -0.5 * paramPtr->dvt1w * paramPtr->weff * paramPtr->leff / tmp2;

    if (T0 > -CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 * (1.0 + 2.0 * T1);
    }
    else
    {
      T1 = CONSTMIN_EXP;
      T2 = T1 * (1.0 + 2.0 * T1);
    }
    T0 = paramPtr->dvt0w * T2;
    T2 = T0 * tmp1;

    T0 = -0.5 * paramPtr->dvt1 * paramPtr->leff / tmp2;

    if (T0 > -CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T3 = T1 * (1.0 + 2.0 * T1);
    }
    else
    {
      T1 = CONSTMIN_EXP;
      T3 = T1 * (1.0 + 2.0 * T1);
    }

    T3 = paramPtr->dvt0 * T3 * tmp1;

    T4 = model_.tox * paramPtr->phi / (paramPtr->weff + paramPtr->w0);

    T0 = sqrt(1.0 + paramPtr->nlx / paramPtr->leff);
    T5 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi
      + (paramPtr->kt1 + paramPtr->kt1l / paramPtr->leff) * (TRatio - 1.0);

    tmp3 = model_.dtype * paramPtr->vth0 - T2 - T3 + paramPtr->k3 * T4 + T5;

    paramPtr->vfbzb = tmp3 - paramPtr->phi - paramPtr->k1 * paramPtr->sqrtPhi;

  } // End of vfbzb

  cgso = paramPtr->cgso;
  cgdo = paramPtr->cgdo;

  Nvtm = vtm * model_.jctEmissionCoeff;

  if ((sourceArea <= 0.0) &&
      (sourcePerimeter <= 0.0))
  {
    SourceSatCurrent = 1.0e-14;
  }
  else
  {
    SourceSatCurrent = sourceArea * jctTempSatCurDensity
      + sourcePerimeter
      * jctSidewallTempSatCurDensity;
  }

  if ((SourceSatCurrent > 0.0) && (model_.ijth > 0.0))
  {
    vjsm = Nvtm * log(model_.ijth / SourceSatCurrent + 1.0);
    IsEvjsm = SourceSatCurrent * exp(vjsm / Nvtm);
  }

  if ((drainArea <= 0.0) &&
      (drainPerimeter <= 0.0))
  {
    DrainSatCurrent = 1.0e-14;
  }
  else
  {
    DrainSatCurrent = drainArea * jctTempSatCurDensity
      + drainPerimeter
      * jctSidewallTempSatCurDensity;
  }

  if ((DrainSatCurrent > 0.0) && (model_.ijth > 0.0))
  {
    vjdm = Nvtm * log(model_.ijth / DrainSatCurrent + 1.0);
    IsEvjdm = DrainSatCurrent * exp(vjdm / Nvtm);
  }

  updateTemperatureCalled_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/09/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

// begin the b3ld.c parameters:
  double SourceSatCurrent(0.0), DrainSatCurrent(0.0);
  double vgdo(0.0);

  double VgstNVt(0.0), ExpVgst(0.0);

  double czbd(0.0), czbdsw(0.0), czbdswg(0.0), czbs(0.0), czbssw(0.0), czbsswg(0.0);
  double evbd(0.0), evbs(0.0), arg(0.0), sarg(0.0);

  double Vfbeff(0.0), dVfbeff_dVg(0.0), dVfbeff_dVb(0.0), V3(0.0), V4(0.0);

  double MJ(0.0), MJSW(0.0), MJSWG(0.0);

  double qinoi(0.0);
  double Vds(0.0), Vgs(0.0), Vbs(0.0);
  double Vgs_eff(0.0), Vfb(0.0);
  double Phis(0.0), dPhis_dVb(0.0), sqrtPhis(0.0), dsqrtPhis_dVb(0.0);
  double Vth(0.0), dVth_dVb(0.0), dVth_dVd(0.0);
  double Vgst(0.0);

  double Nvtm(0.0);
  double Vtm(0.0);
  double n(0.0), dn_dVb(0.0), dn_dVd(0.0), voffcv(0.0), noff(0.0), dnoff_dVd(0.0), dnoff_dVb(0.0);
  double ExpArg(0.0), V0(0.0), CoxWLcen(0.0), QovCox(0.0), LINK(0.0);
  double DeltaPhi(0.0);

  double Cox(0.0), Tox(0.0), Tcen(0.0), dTcen_dVg(0.0), dTcen_dVd(0.0), dTcen_dVb(0.0);
  double Ccen(0.0), Coxeff(0.0), dCoxeff_dVg(0.0), dCoxeff_dVd(0.0), dCoxeff_dVb(0.0);
  double Denomi(0.0), dDenomi_dVg(0.0), dDenomi_dVd(0.0), dDenomi_dVb(0.0);

  double dueff_dVg(0.0), dueff_dVd(0.0), dueff_dVb(0.0);
  double Esat(0.0);

  double Vdsat(0.0);

  double EsatL(0.0), dEsatL_dVg(0.0), dEsatL_dVd(0.0), dEsatL_dVb(0.0);

  double dVdsat_dVg(0.0), dVdsat_dVb(0.0), dVdsat_dVd(0.0), Vasat(0.0), dAlphaz_dVg(0.0), dAlphaz_dVb(0.0);
  double dVasat_dVg(0.0), dVasat_dVb(0.0), dVasat_dVd(0.0), Va(0.0);

  double dVa_dVd(0.0), dVa_dVg(0.0), dVa_dVb(0.0);
  double Vbseff(0.0), dVbseff_dVb(0.0), VbseffCV(0.0), dVbseffCV_dVb(0.0);
  double Arg1(0.0);

  double One_Third_CoxWL(0.0), Two_Third_CoxWL(0.0), Alphaz(0.0);

  double T0(0.0), dT0_dVg(0.0), dT0_dVd(0.0), dT0_dVb(0.0);
  double T1(0.0), dT1_dVg(0.0), dT1_dVd(0.0), dT1_dVb(0.0);
  double T2(0.0), dT2_dVg(0.0), dT2_dVd(0.0), dT2_dVb(0.0);
  double T3(0.0), dT3_dVg(0.0), dT3_dVd(0.0), dT3_dVb(0.0);
  double T4(0.0);

  double T5(0.0);
  double T6(0.0);
  double T7(0.0);
  double T8(0.0);
  double T9(0.0);
  double T10(0.0);
  double T11(0.0), T12(0.0);

  double tmp(0.0), Abulk(0.0), dAbulk_dVb(0.0), Abulk0(0.0), dAbulk0_dVb(0.0);

  double VACLM(0.0), dVACLM_dVg(0.0), dVACLM_dVd(0.0), dVACLM_dVb(0.0);
  double VADIBL(0.0), dVADIBL_dVg(0.0), dVADIBL_dVd(0.0), dVADIBL_dVb(0.0);

  double Xdep(0.0), dXdep_dVb(0.0), lt1(0.0), dlt1_dVb(0.0), ltw(0.0), dltw_dVb(0.0);
  double Delt_vth(0.0), dDelt_vth_dVb(0.0);

  double Theta0(0.0), dTheta0_dVb(0.0);

  double TempRatio(0.0), tmp1(0.0), tmp2(0.0), tmp3(0.0), tmp4(0.0);

  double DIBL_Sft(0.0), dDIBL_Sft_dVd(0.0);

  double Lambda(0.0), dLambda_dVg(0.0);

  double a1(0.0);

  double Vgsteff(0.0), dVgsteff_dVg(0.0), dVgsteff_dVd(0.0), dVgsteff_dVb(0.0);
  double Vdseff(0.0), dVdseff_dVg(0.0), dVdseff_dVd(0.0), dVdseff_dVb(0.0);
  double VdseffCV(0.0), dVdseffCV_dVg(0.0), dVdseffCV_dVd(0.0), dVdseffCV_dVb(0.0);
  double diffVds(0.0);

  double dAbulk_dVg(0.0);
  double beta(0.0), dbeta_dVg(0.0), dbeta_dVd(0.0), dbeta_dVb(0.0);
  double gche(0.0), dgche_dVg(0.0), dgche_dVd(0.0), dgche_dVb(0.0);
  double fgche1(0.0), dfgche1_dVg(0.0), dfgche1_dVd(0.0), dfgche1_dVb(0.0);
  double fgche2(0.0), dfgche2_dVg(0.0), dfgche2_dVd(0.0), dfgche2_dVb(0.0);
  double Idl(0.0), dIdl_dVg(0.0), dIdl_dVd(0.0), dIdl_dVb(0.0);
  double Idsa(0.0), dIdsa_dVg(0.0), dIdsa_dVd(0.0), dIdsa_dVb(0.0);
  double Ids(0.0);

  double Gds(0.0), Gmb(0.0);
  double Isub(0.0);

  double Gbd(0.0), Gbg(0.0), Gbb(0.0);
  double VASCBE(0.0), dVASCBE_dVg(0.0), dVASCBE_dVd(0.0), dVASCBE_dVb(0.0);
  double CoxWovL(0.0);
  double Rds(0.0), dRds_dVg(0.0), dRds_dVb(0.0), WVCox(0.0), WVCoxRds(0.0);
  double Vgst2Vtm(0.0), VdsatCV(0.0);

  double dVdsatCV_dVg(0.0), dVdsatCV_dVb(0.0);
  double Leff(0.0), Weff(0.0), dWeff_dVg(0.0), dWeff_dVb(0.0);
  double AbulkCV(0.0), dAbulkCV_dVb(0.0);

  double gtau_diff(0.0), gtau_drift(0.0);
  // these shadow member varables and then are unitialized
  // when used in later calculations
  // double qcheq(0.0), cqcheq(0.0), qdef(0.0);

  double Cgg1(0.0), Cgb1(0.0), Cgd1(0.0), Cbg1(0.0), Cbb1(0.0), Cbd1(0.0);

  double Qac0(0.0), Qsub0(0.0);
  double dQac0_dVg(0.0), dQac0_dVb(0.0), dQsub0_dVg(0.0), dQsub0_dVd(0.0), dQsub0_dVb(0.0);
  double von_local(0.0);

  ScalingFactor = 1.0e-9;

  // Don't do charge computations in DC sweeps.
  if (solState.tranopFlag || solState.acopFlag || solState.transientFlag)
  {
    ChargeComputationNeeded = true;
  }
  else
  {
    ChargeComputationNeeded = false;
  }

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Instance::updateIntermediateVars\n";
    cout << "  name = " << getName();
    cout << "  model name = " << model_.getName();
    cout <<"   dtype is " << model_.dtype << endl;
    cout.width(21); cout.precision(13); cout.setf(ios::scientific);
    cout << "  " << endl;
  }
#endif

  int Check = 1;

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
  Vs     = 0.0;
  Vg     = 0.0;
  Vb     = 0.0;
  Vsp    = 0.0;
  Vdp    = 0.0;
  Qtotal = 0.0;

  Vd = (extData.nextSolVectorRawPtr)[li_Drain];
  Vg = (extData.nextSolVectorRawPtr)[li_Gate];
  Vs = (extData.nextSolVectorRawPtr)[li_Source];
  Vb = (extData.nextSolVectorRawPtr)[li_Bulk];
  Vsp = (extData.nextSolVectorRawPtr)[li_SourcePrime];
  Vdp = (extData.nextSolVectorRawPtr)[li_DrainPrime];
  if( nqsMod )
  {
    Qtotal = (extData.nextSolVectorRawPtr)[li_Charge];
  }
  else
  {
    Qtotal = 0.0;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout << "     Vg = " << Vg << endl;
    cout << "     Vb = " << Vb << endl;
    cout << "     Vs = " << Vs << endl;
    cout << "     Vd = " << Vd << endl;
    cout << "    Vsp = " << Vsp << endl;
    cout << "    Vdp = " << Vdp << endl;
    cout << " Qtotal = " << Qtotal << endl;
  }
#endif

  Vddp  = Vd   - Vdp;
  Vssp  = Vs   - Vsp;
  Vbsp  = Vb   - Vsp;
  Vbdp  = Vb   - Vdp;
  Vgsp  = Vg   - Vsp;
  Vgdp  = Vg   - Vdp;
  Vgb   = Vg   - Vb;

  Vdpsp = Vdp  - Vsp;

  // modified from b3ld:  (see lines 221-230)
  vbs  = model_.dtype * Vbsp;
  vgs  = model_.dtype * Vgsp;
  vds  = model_.dtype * Vdpsp;

  qdef = model_.dtype * Qtotal;

  vbd = vbs - vds;
  vgd = vgs - vds;

  origFlag = 1;
  limitedFlag = false;
  vgs_orig = vgs;
  vds_orig = vds;
  vbs_orig = vbs;
  vbd_orig = vbd;
  vgd_orig = vgd;

  // What follows is a block of code designed to impose some  limits,
  //  or initial conditions on the junction voltages.  Initial conditions
  //  should only be imposed on the first Newton step of an operating point.
  //
  // The first possible limit on the  junction voltages has to do with
  // limiting the percent change of junction voltages between  Newton
  // iterations.  The second has to do with avoiding extra floating point
  // operations in the event that the device has in some sense converged
  // (aka BYPASS).  Although the primary point of BYPASS is to reduce
  // neccessary work, it also seems to reduce the number of Newton iterations.
  //
  // NOTE:  We do not support BYPASS.
  //
  // The "old" variables should be the values for the previous
  // Newton iteration, if indeed there was a previous Newton
  // iteration.  If not, just set the  old values equal to
  // the current ones.
  //

  // set an initial condition if appropriate:
  if (solState.initJctFlag && !OFF && devOptions.voltageLimiterFlag)
  {
    if (solState.inputOPFlag)
    {
      N_LAS_Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
      if ((*flagSolVectorPtr)[li_Drain] == 0 || (*flagSolVectorPtr)[li_Gate] == 0 ||
          (*flagSolVectorPtr)[li_Source] == 0 || (*flagSolVectorPtr)[li_SourcePrime] == 0 ||
          (*flagSolVectorPtr)[li_DrainPrime] == 0 || (*flagSolVectorPtr)[li_Bulk] == 0 )
      {
        vbs = 0.0;
        vgs = model_.dtype * paramPtr->vth0 + 0.1;
        vds = 0.1;
        origFlag = 0;
      }
    }
    else
    {
      vbs = 0.0;
      vgs = model_.dtype * paramPtr->vth0 + 0.1;
      vds = 0.1;
      origFlag = 0;
    }
    vbd = vbs - vds;
    vgd = vgs - vds;
    //origFlag = 0;
  }
  else if ((solState.initFixFlag || solState.initJctFlag) && OFF)
  {
    qdef = vbs = vgs = vds = 0;
  }


  if (solState.newtonIter == 0)
  {
    newtonIterOld = 0;

    if (!solState.dcopFlag || (solState.locaEnabledFlag && solState.dcopFlag))
    // ie, first newton step of a transient time step or DCOP continuation step.
    {
      vbs_old = (extData.currStoVectorRawPtr)[li_store_vbs];
      vbd_old = (extData.currStoVectorRawPtr)[li_store_vbd];
      vgs_old = (extData.currStoVectorRawPtr)[li_store_vgs];
      vds_old = (extData.currStoVectorRawPtr)[li_store_vds];
      von_local = (extData.currStoVectorRawPtr)[li_store_von];
    }
    else
    {  // no history
      vbs_old = vbs;
      vbd_old = vbd;
      vgs_old = vgs;
      vds_old = vds;
      von_local = 0.0;
    }
  }
  else
  {
    vbs_old = (extData.nextStoVectorRawPtr)[li_store_vbs];
    vbd_old = (extData.nextStoVectorRawPtr)[li_store_vbd];
    vgs_old = (extData.nextStoVectorRawPtr)[li_store_vgs];
    vds_old = (extData.nextStoVectorRawPtr)[li_store_vds];
    von_local = (extData.nextStoVectorRawPtr)[li_store_von];
  }

  vgdo = vgs_old - vds_old;

  // This next block performs checks on the junction voltages and
  // imposes limits on them if they are too big.
  // Note:  In the level=1 von is multiplied by dtype.  Here it is not.  They
  // are both right.

  if (devOptions.voltageLimiterFlag && !(solState.initFixFlag && OFF))
  {

#ifdef Xyce_DEBUG_DEVICE
    if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
    {
      cout.width(21); cout.precision(13); cout.setf(ios::scientific);
      cout << "  von_local = " << von_local << endl;
      cout << "  CONSTvt0  = " << CONSTvt0 << endl;
      cout << "  vcrit     = " << model_.vcrit << endl;
      cout.width(3);
      cout << solState.newtonIter;
      cout.width(5);cout << getName();
      cout << " old :";
      cout<<" vgs:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vgs_old;
      cout<<" vds:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vds_old;
      cout<<" vbs:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vbs_old;
      cout<<" vbd:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vbd_old << endl;
      cout.width(3);
      cout << solState.newtonIter;
      cout.width(5);cout << getName();
      cout << " Blim:";
      cout<<" vgs:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vgs;
      cout<<" vds:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vds;
      cout<<" vbs:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vbs;
      cout<<" vbd:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vbd << endl;
      cout.width(21); cout.precision(13); cout.setf(ios::scientific);
    }
#endif

    // only do this if we are beyond the first Newton iteration.  On the
    // first newton iteration, the "old" values are from a previous time
    // step.

    if (solState.newtonIter >= 0)
    {
      if (vds_old >= 0.0)
      {
        vgs = devSupport.fetlim( vgs, vgs_old, von_local);
        vds = vgs - vgd;
        vds = devSupport.limvds( vds,  vds_old);
        vgd = vgs - vds;
      }
      else
      {
        vgd = devSupport.fetlim( vgd, vgdo, von_local);
        vds = vgs - vgd;
        vds = -devSupport.limvds( -vds, -vds_old );
        vgs = vgd + vds;
      }

      if (vds >= 0.0)
      {
        vbs = devSupport.pnjlim( vbs, vbs_old, CONSTvt0,
                   model_.vcrit, &Check);
        vbd = vbs - vds;
      }
      else
      {
        vbd = devSupport.pnjlim( vbd, vbd_old, CONSTvt0,
                   model_.vcrit, &Check);
        vbs = vbd + vds;
      }
    }

    // set the origFlag:
#ifdef Xyce_NEW_ORIG_TEST
    double vgs_diff = fabs(vgs - vgs_orig);
    double vds_diff = fabs(vds - vds_orig);
    double vbs_diff = fabs(vbs - vbs_orig);
    double vgd_diff = fabs(vgd - vgd_orig);

    bool noOrigFlag_vgs = 0;
    bool noOrigFlag_vds = 0;
    bool noOrigFlag_vbs = 0;
    bool noOrigFlag_vgd = 0;

    if (vgs_diff != 0.0)
    {
      if (vgs_orig != 0.0)
      {
        if ( fabs(vgs_diff/vgs_orig) > reltol) noOrigFlag_vgs = 1;
      }
      else
      {
        if ( fabs(vgs_diff) > voltTol ) noOrigFlag_vgs = 1;
      }
    }

    if (vds_diff != 0.0)
    {
      if (vds_orig != 0.0)
      {
        if ( fabs(vds_diff/vds_orig) > reltol) noOrigFlag_vds = 1;
      }
      else
      {
        if ( fabs(vds_diff) > voltTol ) noOrigFlag_vds = 1;
      }
    }

    if (vbs_diff != 0.0)
    {
      if (vbs_orig != 0.0)
      {
        if ( fabs(vbs_diff/vbs_orig) > reltol) noOrigFlag_vbs = 1;
      }
      else
      {
        if ( fabs(vbs_diff) > voltTol ) noOrigFlag_vbs = 1;
      }
    }

    if (vgd_diff != 0.0)
    {
      if (vgd_orig != 0.0)
      {
        if ( fabs(vgd_diff/vgd_orig) > reltol) noOrigFlag_vgd = 1;
      }
      else
      {
        if ( fabs(vgd_diff) > voltTol ) noOrigFlag_vgd = 1;
      }
    }

    origFlag = !( noOrigFlag_vgs || noOrigFlag_vds ||
                  noOrigFlag_vbs || noOrigFlag_vgd);

#else
    if (vgs_orig != vgs || vds_orig != vds ||
        vbs_orig != vbs || vbd_orig != vbd || vgd_orig != vgd) origFlag = 0;
#endif

    // for convergence testing:
    if (Check == 1) limitedFlag=true;

#ifdef Xyce_DEBUG_DEVICE
    if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
    {
      cout.width(3);
      cout << solState.newtonIter;
      cout.width(5);cout << getName();
      cout << " Alim:";
      cout<<" vgs:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vgs;
      cout<<" vds:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vds;
      cout<<" vbs:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vbs;
      cout<<" vbd:";cout.width(10);cout.precision(3);cout.setf(ios::scientific);
      cout << vbd;
      if (origFlag) cout << " SAME";
      else          cout << " DIFF";
      cout << endl;
      cout.width(21); cout.precision(13); cout.setf(ios::scientific);
    }
#endif

  } // devOptions.voltageLimiterFlag

  // update the "old" variables:
  if (solState.newtonIter != 0 && solState.newtonIter != newtonIterOld)
  {
    newtonIterOld = solState.newtonIter;
  }

  // Finished with what would have been the series of CKTmode
  // IF statements...

  // file: b3ld.c   line: 345
  // determine DC current and derivatives
  vbd = vbs - vds;
  vgd = vgs - vds;
  vgb = vgs - vbs;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    cout << " mod type = " << model_.modType << endl;
    cout << "dtype = " << model_.dtype << endl;
    cout << "  vbs = " << vbs << endl;
    cout << "  vds = " << vds << endl;
    cout << "  vgs = " << vgs << endl;

    cout << "  vbd = " << vbd << endl;
    cout << "  vgd = " << vgd << endl;
    cout << "  vgb = " << vgs << endl;
    cout << " qdef = " << qdef << endl;
  }
#endif

  // Source/drain junction diode DC model begins
  Nvtm = vtm * model_.jctEmissionCoeff;

  if ((sourceArea <= 0.0) && (sourcePerimeter <= 0.0))
  {
    SourceSatCurrent = 1.0e-14;
  }
  else
  {
    SourceSatCurrent = sourceArea
      * jctTempSatCurDensity
      + sourcePerimeter
      * jctSidewallTempSatCurDensity;
  }

  if (SourceSatCurrent <= 0.0)
  {
    gbs = devOptions.gmin;
    cbs = gbs * vbs;
  }
  else
  {
    if (model_.ijth == 0.0)
    {
      evbs = exp(vbs / Nvtm);
      gbs = SourceSatCurrent * evbs / Nvtm + devOptions.gmin;
      cbs = SourceSatCurrent * (evbs - 1.0) + devOptions.gmin * vbs;
    }
    else
    {
      if (vbs < vjsm)
      {
        evbs = exp(vbs / Nvtm);
        gbs = SourceSatCurrent * evbs / Nvtm + devOptions.gmin;
        cbs = SourceSatCurrent * (evbs - 1.0) + devOptions.gmin * vbs;
      }
      else
      {
        T0 = IsEvjsm / Nvtm;
        gbs = T0 + devOptions.gmin;
        cbs = IsEvjsm - SourceSatCurrent
          + T0 * (vbs - vjsm)
          + devOptions.gmin * vbs;
      }
    }
  }

  if ((drainArea <= 0.0) && (drainPerimeter <= 0.0))
  {
    DrainSatCurrent = 1.0e-14;
  }
  else
  {
    DrainSatCurrent = drainArea
      * jctTempSatCurDensity
      + drainPerimeter
      * jctSidewallTempSatCurDensity;
  }

  if (DrainSatCurrent <= 0.0)
  {
    gbd = devOptions.gmin;
    cbd = gbd * vbd;
  }
  else
  {
    if (model_.ijth == 0.0)
    {
      evbd = exp(vbd / Nvtm);
      gbd = DrainSatCurrent * evbd / Nvtm + devOptions.gmin;
      cbd = DrainSatCurrent * (evbd - 1.0) + devOptions.gmin * vbd;
    }
    else
    {
      if (vbd < vjdm)
      {
        evbd = exp(vbd / Nvtm);
        gbd = DrainSatCurrent * evbd / Nvtm + devOptions.gmin;
        cbd = DrainSatCurrent * (evbd - 1.0) + devOptions.gmin * vbd;
      }
      else
      {
        T0 = IsEvjdm / Nvtm;
        gbd = T0 + devOptions.gmin;
        cbd = IsEvjdm - DrainSatCurrent
          + T0 * (vbd - vjdm)
          + devOptions.gmin * vbd;
      }
    }
  }
  // End of diode DC model

  if (vds >= 0.0)
  {   // normal mode
    mode = 1;
    Vds = vds;
    Vgs = vgs;
    Vbs = vbs;
  }
  else
  {   // inverse mode
    mode = -1;
    Vds = -vds;
    Vgs = vgd;
    Vbs = vbd;
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

    Vds = devSupport.contVds (Vds,solState.nltermScale,devOptions.vdsScaleMin);
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

  T0 = Vbs - paramPtr->vbsc - 0.001;
  T1global = sqrt(T0 * T0 - 0.004 * paramPtr->vbsc);
  Vbseff = paramPtr->vbsc + 0.5 * (T0 + T1global);
  dVbseff_dVb = 0.5 * (1.0 + T0 / T1global);

  if (Vbseff < Vbs) Vbseff = Vbs;

  if (Vbseff > 0.0)
  {
    T0 = paramPtr->phi / (paramPtr->phi + Vbseff);
    Phis = paramPtr->phi * T0;
    dPhis_dVb = -T0 * T0;
    sqrtPhis = paramPtr->phis3 / (paramPtr->phi + 0.5 * Vbseff);
    dsqrtPhis_dVb = -0.5 * sqrtPhis * sqrtPhis / paramPtr->phis3;
  }
  else
  {
    Phis = paramPtr->phi - Vbseff;
    dPhis_dVb = -1.0;
    sqrtPhis = sqrt(Phis);
    dsqrtPhis_dVb = -0.5 / sqrtPhis;
  }

  Xdep = paramPtr->Xdep0 * sqrtPhis / paramPtr->sqrtPhi;
  dXdep_dVb = (paramPtr->Xdep0 / paramPtr->sqrtPhi) * dsqrtPhis_dVb;

  Leff = paramPtr->leff;
  Vtm = vtm;
  // Vth Calculation
  T3 = sqrt(Xdep);
  V0 = paramPtr->vbi - paramPtr->phi;

  T0 = paramPtr->dvt2 * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = paramPtr->dvt2;
  }
  else // Added to avoid any discontinuity problems caused by dvt2
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = paramPtr->dvt2 * T4 * T4;
  }

  lt1 = model_.factor1 * T3 * T1;
  dlt1_dVb = model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = paramPtr->dvt2w * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = paramPtr->dvt2w;
  }
  else // Added to avoid any discontinuity problems caused by dvt2w
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = paramPtr->dvt2w * T4 * T4;
  }

  ltw = model_.factor1 * T3 * T1;
  dltw_dVb = model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = -0.5 * paramPtr->dvt1 * Leff / lt1;
  if (T0 > -CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    Theta0 = T1 * (1.0 + 2.0 * T1);
    dT1_dVb = -T0 / lt1 * T1 * dlt1_dVb;
    dTheta0_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
  }
  else
  {
    T1 = CONSTMIN_EXP;
    Theta0 = T1 * (1.0 + 2.0 * T1);
    dTheta0_dVb = 0.0;
  }

  thetavth = paramPtr->dvt0 * Theta0;
  Delt_vth = thetavth * V0;
  dDelt_vth_dVb = paramPtr->dvt0 * dTheta0_dVb * V0;

  T0 = -0.5 * paramPtr->dvt1w * paramPtr->weff * Leff / ltw;
  if (T0 > -CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 * (1.0 + 2.0 * T1);
    dT1_dVb = -T0 / ltw * T1 * dltw_dVb;
    dT2_dVb = (1.0 + 4.0 * T1) * dT1_dVb;
  }
  else
  {
    T1 = CONSTMIN_EXP;
    T2 = T1 * (1.0 + 2.0 * T1);
    dT2_dVb = 0.0;
  }

  T0 = paramPtr->dvt0w * T2;
  T2 = T0 * V0;
  dT2_dVb = paramPtr->dvt0w * dT2_dVb * V0;

  TempRatio =  temp / model_.tnom - 1.0;
  T0 = sqrt(1.0 + paramPtr->nlx / Leff);
  T1 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi
    + (paramPtr->kt1 + paramPtr->kt1l / Leff
       +  paramPtr->kt2 * Vbseff) * TempRatio;

  tmp2 = model_.tox * paramPtr->phi / (paramPtr->weff + paramPtr->w0);

  T3 = paramPtr->eta0 + paramPtr->etab * Vbseff;
  if (T3 < 1.0e-4) // avoid  discontinuity problems caused by etab
  {
    T9 = 1.0 / (3.0 - 2.0e4 * T3);
    T3 = (2.0e-4 - T3) * T9;
    T4 = T9 * T9;
  }
  else
  {
    T4 = 1.0;
  }

  dDIBL_Sft_dVd = T3 * paramPtr->theta0vb0;
  DIBL_Sft = dDIBL_Sft_dVd * Vds;

  Vth = model_.dtype * paramPtr->vth0 - paramPtr->k1
    * paramPtr->sqrtPhi + paramPtr->k1ox * sqrtPhis
    - paramPtr->k2ox * Vbseff - Delt_vth - T2 + (paramPtr->k3
    + paramPtr->k3b * Vbseff) * tmp2 + T1 - DIBL_Sft;

  von = Vth;

  dVth_dVb = paramPtr->k1ox * dsqrtPhis_dVb - paramPtr->k2ox
    - dDelt_vth_dVb - dT2_dVb + paramPtr->k3b * tmp2
    - paramPtr->etab * Vds * paramPtr->theta0vb0 * T4
    + paramPtr->kt2 * TempRatio;

  dVth_dVd = -dDIBL_Sft_dVd;

  // Calculate n
  tmp2 = paramPtr->nfactor * CONSTEPSSI / Xdep;
  tmp3 = paramPtr->cdsc + paramPtr->cdscb * Vbseff
    + paramPtr->cdscd * Vds;
  tmp4 = (tmp2 + tmp3 * Theta0 + paramPtr->cit) / model_.cox;

  if (tmp4 >= -0.5)
  {
    n = 1.0 + tmp4;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
              + paramPtr->cdscb * Theta0) / model_.cox;
    dn_dVd = paramPtr->cdscd * Theta0 / model_.cox;
  }
  else // avoid  discontinuity problems caused by tmp4
  {
    T0 = 1.0 / (3.0 + 8.0 * tmp4);
    n = (1.0 + 3.0 * tmp4) * T0;
    T0 *= T0;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
              + paramPtr->cdscb * Theta0) / model_.cox * T0;
    dn_dVd = paramPtr->cdscd * Theta0 / model_.cox * T0;
  }

  // Poly Gate Si Depletion Effect
  T0 = paramPtr->vfb + paramPtr->phi;

  // added to avoid the problem caused by ngate
  if ((paramPtr->ngate > 1.e18) && (paramPtr->ngate < 1.e25) && (Vgs > T0))
  {
    T1 = 1.0e6 * CONSTQ * CONSTEPSSI * paramPtr->ngate
      / (model_.cox * model_.cox);
    T4 = sqrt(1.0 + 2.0 * (Vgs - T0) / T1);

    T2 = T1 * (T4 - 1.0);
    T3 = 0.5 * T2 * T2 / T1; // T3 = Vpoly
    T7 = 1.12 - T3 - 0.05;
    T6 = sqrt(T7 * T7 + 0.224);
    T5 = 1.12 - 0.5 * (T7 + T6);
    Vgs_eff = Vgs - T5;
    dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
  }
  else
  {
    Vgs_eff = Vgs;
    dVgs_eff_dVg = 1.0;
  }
  Vgst = Vgs_eff - Vth;

  // Effective Vgst (Vgsteff) Calculation
  T10 = 2.0 * n * Vtm;
  VgstNVt = Vgst / T10;
  ExpArg = (2.0 * paramPtr->voff - Vgst) / T10;

  // MCJ: Very small Vgst
  if (VgstNVt > CONSTEXP_THRESHOLD)
  {
    Vgsteff = Vgst;
    dVgsteff_dVg = dVgs_eff_dVg;
    dVgsteff_dVd = -dVth_dVd;
    dVgsteff_dVb = -dVth_dVb;
  }
  else if (ExpArg > CONSTEXP_THRESHOLD)
  {
    T0 = (Vgst - paramPtr->voff) / (n * Vtm);
    ExpVgst = exp(T0);
    Vgsteff = Vtm * paramPtr->cdep0 / model_.cox * ExpVgst;
    dVgsteff_dVg = Vgsteff / (n * Vtm);
    dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + T0 * Vtm * dn_dVd);
    dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + T0 * Vtm * dn_dVb);
    dVgsteff_dVg *= dVgs_eff_dVg;
  }
  else
  {
    ExpVgst = exp(VgstNVt);
    T1 = T10 * log(1.0 + ExpVgst);
    dT1_dVg = ExpVgst / (1.0 + ExpVgst);
    dT1_dVb = -dT1_dVg * (dVth_dVb + Vgst / n * dn_dVb) + T1 / n * dn_dVb;
    dT1_dVd = -dT1_dVg * (dVth_dVd + Vgst / n * dn_dVd) + T1 / n * dn_dVd;

    dT2_dVg = -model_.cox / (Vtm * paramPtr->cdep0) * exp(ExpArg);
    T2 = 1.0 - T10 * dT2_dVg;

    dT2_dVd = -dT2_dVg * (dVth_dVd - 2.0 * Vtm * ExpArg * dn_dVd)
      + (T2 - 1.0) / n * dn_dVd;

    dT2_dVb = -dT2_dVg * (dVth_dVb - 2.0 * Vtm * ExpArg * dn_dVb)
      + (T2 - 1.0) / n * dn_dVb;

    Vgsteff = T1 / T2;
    T3 = T2 * T2;

    dVgsteff_dVg = (T2 * dT1_dVg - T1 * dT2_dVg) / T3 * dVgs_eff_dVg;
    dVgsteff_dVd = (T2 * dT1_dVd - T1 * dT2_dVd) / T3;
    dVgsteff_dVb = (T2 * dT1_dVb - T1 * dT2_dVb) / T3;
  }

  // Calculate Effective Channel Geometry
  T9 = sqrtPhis - paramPtr->sqrtPhi;
  Weff = paramPtr->weff -2.0 *(paramPtr->dwg * Vgsteff + paramPtr->dwb * T9);
  dWeff_dVg = -2.0 * paramPtr->dwg;
  dWeff_dVb = -2.0 * paramPtr->dwb * dsqrtPhis_dVb;

  if (Weff < 2.0e-8) // to avoid the discontinuity problem due to Weff
  {
    T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
    Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
    T0 *= T0 * 4.0e-16;
    dWeff_dVg *= T0;
    dWeff_dVb *= T0;
  }

  T0 = paramPtr->prwg * Vgsteff + paramPtr->prwb * T9;
  if (T0 >= -0.9)
  {
    Rds = paramPtr->rds0 * (1.0 + T0);
    dRds_dVg = paramPtr->rds0 * paramPtr->prwg;
    dRds_dVb = paramPtr->rds0 * paramPtr->prwb * dsqrtPhis_dVb;
  }
  else // to avoid the discontinuity problem due to prwg and prwb
  {
    T1 = 1.0 / (17.0 + 20.0 * T0);
    Rds = paramPtr->rds0 * (0.8 + T0) * T1;
    T1 *= T1;
    dRds_dVg = paramPtr->rds0 * paramPtr->prwg * T1;
    dRds_dVb = paramPtr->rds0 * paramPtr->prwb * dsqrtPhis_dVb * T1;
  }

  // Calculate Abulk
  T1 = 0.5 * paramPtr->k1ox / sqrtPhis;
  dT1_dVb = -T1 / sqrtPhis * dsqrtPhis_dVb;

  T9 = sqrt(paramPtr->xj * Xdep);
  tmp1 = Leff + 2.0 * T9;
  T5 = Leff / tmp1;
  tmp2 = paramPtr->a0 * T5;
  tmp3 = paramPtr->weff + paramPtr->b1;
  tmp4 = paramPtr->b0 / tmp3;
  T2 = tmp2 + tmp4;
  dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
  T6 = T5 * T5;
  T7 = T5 * T6;

  Abulk0 = 1.0 + T1 * T2;
  dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

  T8 = paramPtr->ags * paramPtr->a0 * T7;
  dAbulk_dVg = -T1 * T8;
  Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
  dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb + 3.0 * T1 * dT2_dVb);

  if (Abulk0 < 0.1) // added to avoid the problems caused by Abulk0
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk0);
    Abulk0 = (0.2 - Abulk0) * T9;
    dAbulk0_dVb *= T9 * T9;
  }

  if (Abulk < 0.1) // added to avoid the problems caused by Abulk
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk);
    Abulk = (0.2 - Abulk) * T9;
    T10 = T9 * T9;
    dAbulk_dVb *= T10;
    dAbulk_dVg *= T10;
  }

  T2 = paramPtr->keta * Vbseff;
  if (T2 >= -0.9)
  {
    T0 = 1.0 / (1.0 + T2);
    dT0_dVb = -paramPtr->keta * T0 * T0;
  }
  else // added to avoid the problems caused by Keta
  {
    T1 = 1.0 / (0.8 + T2);
    T0 = (17.0 + 20.0 * T2) * T1;
    dT0_dVb = -paramPtr->keta * T1 * T1;
  }

  dAbulk_dVg *= T0;
  dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
  dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
  Abulk *= T0;
  Abulk0 *= T0;

  // Mobility calculation
  if (model_.mobMod == 1)
  {
    T0 = Vgsteff + Vth + Vth;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;
    T3 = T0 / model_.tox;
    T5 = T3 * (T2 + paramPtr->ub * T3);
    dDenomi_dVg = (T2 + 2.0 * paramPtr->ub * T3) / model_.tox;
    dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
    dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + paramPtr->uc * T3;
  }
  else if (model_.mobMod == 2)
  {
    T5 = Vgsteff / model_.tox * (paramPtr->ua
         + paramPtr->uc * Vbseff + paramPtr->ub * Vgsteff / model_.tox);

    dDenomi_dVg = (paramPtr->ua + paramPtr->uc * Vbseff
        + 2.0 * paramPtr->ub * Vgsteff / model_.tox) / model_.tox;

    dDenomi_dVd = 0.0;
    dDenomi_dVb = Vgsteff * paramPtr->uc / model_.tox;
  }
  else
  {
    T0 = Vgsteff + Vth + Vth;
    T2 = 1.0 + paramPtr->uc * Vbseff;
    T3 = T0 / model_.tox;
    T4 = T3 * (paramPtr->ua + paramPtr->ub * T3);
    T5 = T4 * T2;

    dDenomi_dVg = (paramPtr->ua + 2.0 * paramPtr->ub * T3) * T2 /model_.tox;
    dDenomi_dVd = dDenomi_dVg * 2.0 * dVth_dVd;
    dDenomi_dVb = dDenomi_dVg * 2.0 * dVth_dVb + paramPtr->uc * T4;
  }

  if (T5 >= -0.8)
  {
    Denomi = 1.0 + T5;
  }
  else // Added to avoid the discontinuity problem caused by ua and ub
  {
    T9 = 1.0 / (7.0 + 10.0 * T5);
    Denomi = (0.6 + T5) * T9;
    T9 *= T9;
    dDenomi_dVg *= T9;
    dDenomi_dVd *= T9;
    dDenomi_dVb *= T9;
  }

  ueff = paramPtr->u0temp / Denomi;
  T9 = -ueff / Denomi;
  dueff_dVg = T9 * dDenomi_dVg;
  dueff_dVd = T9 * dDenomi_dVd;
  dueff_dVb = T9 * dDenomi_dVb;

  // Saturation Drain Voltage  Vdsat
  WVCox = Weff * paramPtr->vsattemp * model_.cox;
  WVCoxRds = WVCox * Rds;

  Esat = 2.0 * paramPtr->vsattemp / ueff;
  EsatL = Esat * Leff;
  T0 = -EsatL /ueff;
  dEsatL_dVg = T0 * dueff_dVg;
  dEsatL_dVd = T0 * dueff_dVd;
  dEsatL_dVb = T0 * dueff_dVb;

  // Sqrt()
  a1 = paramPtr->a1;
  if (a1 == 0.0)
  {
    Lambda = paramPtr->a2;
    dLambda_dVg = 0.0;
  }
  else if (a1 > 0.0) // Added to avoid the discontinuity problem
    // caused by a1 and a2 (Lambda)
  {
    T0 = 1.0 - paramPtr->a2;
    T1 = T0 - paramPtr->a1 * Vgsteff - 0.0001;
    T2 = sqrt(T1 * T1 + 0.0004 * T0);
    Lambda = paramPtr->a2 + T0 - 0.5 * (T1 + T2);
    dLambda_dVg = 0.5 * paramPtr->a1 * (1.0 + T1 / T2);
  }
  else
  {
    T1 = paramPtr->a2 + paramPtr->a1 * Vgsteff - 0.0001;
    T2 = sqrt(T1 * T1 + 0.0004 * paramPtr->a2);
    Lambda = 0.5 * (T1 + T2);
    dLambda_dVg = 0.5 * paramPtr->a1 * (1.0 + T1 / T2);
  }

  Vgst2Vtm = Vgsteff + 2.0 * Vtm;
  if (Rds > 0)
  {
    tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
    tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
  }
  else
  {
    tmp2 = dWeff_dVg / Weff;
    tmp3 = dWeff_dVb / Weff;
  }

  if ((Rds == 0.0) && (Lambda == 1.0))
  {
    T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
    tmp1 = 0.0;
    T1 = T0 * T0;
    T2 = Vgst2Vtm * T0;
    T3 = EsatL * Vgst2Vtm;
    Vdsat = T3 * T0;

    dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
    dT0_dVd = -(Abulk * dEsatL_dVd) * T1;
    dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;

    dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
    dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
    dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
  }
  else
  {
    tmp1 = dLambda_dVg / (Lambda * Lambda);
    T9 = Abulk * WVCoxRds;
    T8 = Abulk * T9;
    T7 = Vgst2Vtm * T9;
    T6 = Vgst2Vtm * WVCoxRds;
    T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
    dT0_dVg = 2.0 * (T8 * tmp2 - Abulk * tmp1
                     + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);

    dT0_dVb = 2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3)
                     + (1.0 / Lambda - 1.0) * dAbulk_dVb);
    dT0_dVd = 0.0;
    T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;

    dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
      + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
                                              + T7 * tmp2 + T6 * dAbulk_dVg);

    dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
      + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);

    dT1_dVd = Abulk * dEsatL_dVd;

    T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
    dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
      + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);

    dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
    dT2_dVd = Vgst2Vtm * dEsatL_dVd;

    T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
    Vdsat = (T1 - T3) / T0;

    dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg)) / T3;
    dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd)) / T3;
    dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb)) / T3;

    dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
                             - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;

    dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
                             - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;

    dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
  }
  vdsat = Vdsat;

  // Effective Vds (Vdseff) Calculation
  T1 = Vdsat - Vds - paramPtr->delta;
  dT1_dVg = dVdsat_dVg;
  dT1_dVd = dVdsat_dVd - 1.0;
  dT1_dVb = dVdsat_dVb;

  T2 = sqrt(T1 * T1 + 4.0 * paramPtr->delta * Vdsat);
  T0 = T1 / T2;
  T3 = 2.0 * paramPtr->delta / T2;
  dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
  dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
  dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

  Vdseff = Vdsat - 0.5 * (T1 + T2);
  dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
  dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
  dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);

  // Added to eliminate non-zero Vdseff at Vds=0.0
  if (Vds == 0.0)
  {
    Vdseff = 0.0;
    dVdseff_dVg = 0.0;
    dVdseff_dVb = 0.0;
  }

  // Calculate VAsat
  tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
  T9 = WVCoxRds * Vgsteff;
  T8 = T9 / Vgst2Vtm;
  T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

  T7 = 2.0 * WVCoxRds * tmp4;
  dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
    - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
            + Vdsat * dAbulk_dVg);

  dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
    - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
  dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

  T9 = WVCoxRds * Abulk;
  T1 = 2.0 / Lambda - 1.0 + T9;
  dT1_dVg = -2.0 * tmp1 +  WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
  dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

  Vasat = T0 / T1;
  dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
  dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
  dVasat_dVd = dT0_dVd / T1;

  if (Vdseff > Vds) Vdseff = Vds;

  diffVds = Vds - Vdseff;

  // Calculate VACLM
  if ((paramPtr->pclm > 0.0) && (diffVds > 1.0e-10))
  {
    T0 = 1.0 / (paramPtr->pclm * Abulk * paramPtr->litl);
    dT0_dVb = -T0 / Abulk * dAbulk_dVb;
    dT0_dVg = -T0 / Abulk * dAbulk_dVg;

    T2 = Vgsteff / EsatL;
    T1 = Leff * (Abulk + T2);
    dT1_dVg = Leff * ((1.0 - T2 * dEsatL_dVg) / EsatL + dAbulk_dVg);
    dT1_dVb = Leff * (dAbulk_dVb - T2 * dEsatL_dVb / EsatL);
    dT1_dVd = -T2 * dEsatL_dVd / Esat;

    T9 = T0 * T1;
    VACLM = T9 * diffVds;
    dVACLM_dVg = T0 * dT1_dVg * diffVds - T9 * dVdseff_dVg
      + T1 * diffVds * dT0_dVg;

    dVACLM_dVb = (dT0_dVb * T1 + T0 * dT1_dVb) * diffVds
      - T9 * dVdseff_dVb;

    dVACLM_dVd = T0 * dT1_dVd * diffVds + T9 * (1.0 - dVdseff_dVd);
  }
  else
  {
    VACLM = CONSTMAX_EXP;
    dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
  }

  // Calculate VADIBL
  if (paramPtr->thetaRout > 0.0)
  {
    T8 = Abulk * Vdsat;
    T0 = Vgst2Vtm * T8;
    dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8
      + Vgst2Vtm * Vdsat * dAbulk_dVg;

    dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
    dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

    T1 = Vgst2Vtm + T8;
    dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
    dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
    dT1_dVd = Abulk * dVdsat_dVd;

    T9 = T1 * T1;
    T2 = paramPtr->thetaRout;

    VADIBL = (Vgst2Vtm - T0 / T1) / T2;
    dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
    dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
    dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

    T7 = paramPtr->pdiblb * Vbseff;
    if (T7 >= -0.9)
    {
      T3 = 1.0 / (1.0 + T7);
      VADIBL *= T3;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = (dVADIBL_dVb - VADIBL * paramPtr->pdiblb) * T3;
      dVADIBL_dVd *= T3;
    }
    else // Added to avoid the discontinuity problem caused by pdiblcb
    {
      T4 = 1.0 / (0.8 + T7);
      T3 = (17.0 + 20.0 * T7) * T4;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = dVADIBL_dVb * T3 - VADIBL * paramPtr->pdiblb * T4 * T4;

      dVADIBL_dVd *= T3;
      VADIBL *= T3;
    }
  }
  else
  {
    VADIBL = CONSTMAX_EXP;
    dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
  }

  // Calculate VA
  T8 = paramPtr->pvag / EsatL;
  T9 = T8 * Vgsteff;
  if (T9 > -0.9)
  {
    T0 = 1.0 + T9;
    dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
    dT0_dVb = -T9 * dEsatL_dVb / EsatL;
    dT0_dVd = -T9 * dEsatL_dVd / EsatL;
  }
  else /* Added to avoid the discontinuity problems caused by pvag */
  {
    T1 = 1.0 / (17.0 + 20.0 * T9);
    T0 = (0.8 + T9) * T1;
    T1 *= T1;
    dT0_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T1;

    T9 *= T1 / EsatL;
    dT0_dVb = -T9 * dEsatL_dVb;
    dT0_dVd = -T9 * dEsatL_dVd;
  }

  tmp1 = VACLM * VACLM;
  tmp2 = VADIBL * VADIBL;
  tmp3 = VACLM + VADIBL;

  T1 = VACLM * VADIBL / tmp3;
  tmp3 *= tmp3;
  dT1_dVg = (tmp1 * dVADIBL_dVg + tmp2 * dVACLM_dVg) / tmp3;
  dT1_dVd = (tmp1 * dVADIBL_dVd + tmp2 * dVACLM_dVd) / tmp3;
  dT1_dVb = (tmp1 * dVADIBL_dVb + tmp2 * dVACLM_dVb) / tmp3;

  Va = Vasat + T0 * T1;
  dVa_dVg = dVasat_dVg + T1 * dT0_dVg + T0 * dT1_dVg;
  dVa_dVd = dVasat_dVd + T1 * dT0_dVd + T0 * dT1_dVd;
  dVa_dVb = dVasat_dVb + T1 * dT0_dVb + T0 * dT1_dVb;

  // Calculate VASCBE
  if (paramPtr->pscbe2 > 0.0)
  {
    if (diffVds > paramPtr->pscbe1 * paramPtr->litl / CONSTEXP_THRESHOLD)
    {
      T0 =  paramPtr->pscbe1 * paramPtr->litl / diffVds;
      VASCBE = Leff * exp(T0) / paramPtr->pscbe2;
      T1 = T0 * VASCBE / diffVds;
      dVASCBE_dVg = T1 * dVdseff_dVg;
      dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
      dVASCBE_dVb = T1 * dVdseff_dVb;
    }
    else
    {
      VASCBE = CONSTMAX_EXP * Leff/paramPtr->pscbe2;
      dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
    }
  }
  else
  {
    VASCBE = CONSTMAX_EXP;
    dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
  }

  // Calculate Ids
  CoxWovL = model_.cox * Weff / Leff;
  beta = ueff * CoxWovL;
  dbeta_dVg = CoxWovL * dueff_dVg + beta * dWeff_dVg / Weff;
  dbeta_dVd = CoxWovL * dueff_dVd;
  dbeta_dVb = CoxWovL * dueff_dVb + beta * dWeff_dVb / Weff;

  T0 = 1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm;
  dT0_dVg = -0.5 * (Abulk * dVdseff_dVg
                    - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg) / Vgst2Vtm;
  dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
  dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff) / Vgst2Vtm;

  fgche1 = Vgsteff * T0;
  dfgche1_dVg = Vgsteff * dT0_dVg + T0;
  dfgche1_dVd = Vgsteff * dT0_dVd;
  dfgche1_dVb = Vgsteff * dT0_dVb;

  T9 = Vdseff / EsatL;
  fgche2 = 1.0 + T9;
  dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
  dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
  dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

  gche = beta * fgche1 / fgche2;
  dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
               - gche * dfgche2_dVg) / fgche2;

  dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
               - gche * dfgche2_dVd) / fgche2;

  dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
               - gche * dfgche2_dVb) / fgche2;

  T0 = 1.0 + gche * Rds;
  T9 = Vdseff / T0;
  Idl = gche * T9;

  dIdl_dVg = (gche * dVdseff_dVg + T9 * dgche_dVg) / T0
    - Idl * gche / T0 * dRds_dVg ;

  dIdl_dVd = (gche * dVdseff_dVd + T9 * dgche_dVd) / T0;
  dIdl_dVb = (gche * dVdseff_dVb + T9 * dgche_dVb
              - Idl * dRds_dVb * gche) / T0;

  T9 =  diffVds / Va;
  T0 =  1.0 + T9;
  Idsa = Idl * T0;
  dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVa_dVg) / Va;
  dIdsa_dVd = T0 * dIdl_dVd + Idl * (1.0 - dVdseff_dVd
                                     - T9 * dVa_dVd) / Va;

  dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVa_dVb) / Va;

  T9 = diffVds / VASCBE;
  T0 = 1.0 + T9;
  Ids = Idsa * T0;

  Gm  = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
  Gds = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd
                                 - T9 * dVASCBE_dVd) / VASCBE;
  Gmb = T0 * dIdsa_dVb - Idsa * (dVdseff_dVb
                                 + T9 * dVASCBE_dVb) / VASCBE;

  Gds += Gm * dVgsteff_dVd;
  Gmb += Gm * dVgsteff_dVb;
  Gm *= dVgsteff_dVg;
  Gmb *= dVbseff_dVb;

  // Substrate current begins
  tmp = paramPtr->alpha0 + paramPtr->alpha1 * Leff;
  if ((tmp <= 0.0) || (paramPtr->beta0 <= 0.0))
  {
    Isub = Gbd = Gbb = Gbg = 0.0;
  }
  else
  {
    T2 = tmp / Leff;
    if (diffVds > paramPtr->beta0 / CONSTEXP_THRESHOLD)
    {
      T0 = -paramPtr->beta0 / diffVds;
      T1 = T2 * diffVds * exp(T0);
      T3 = T1 / diffVds * (T0 - 1.0);
      dT1_dVg = T3 * dVdseff_dVg;
      dT1_dVd = T3 * (dVdseff_dVd - 1.0);
      dT1_dVb = T3 * dVdseff_dVb;
    }
    else
    {
      T3 = T2 * CONSTMIN_EXP;
      T1 = T3 * diffVds;
      dT1_dVg = -T3 * dVdseff_dVg;
      dT1_dVd = T3 * (1.0 - dVdseff_dVd);
      dT1_dVb = -T3 * dVdseff_dVb;
    }
    Isub = T1 * Idsa;
    Gbg = T1 * dIdsa_dVg + Idsa * dT1_dVg;
    Gbd = T1 * dIdsa_dVd + Idsa * dT1_dVd;
    Gbb = T1 * dIdsa_dVb + Idsa * dT1_dVb;

    Gbd += Gbg * dVgsteff_dVd;
    Gbb += Gbg * dVgsteff_dVb;
    Gbg *= dVgsteff_dVg;
    Gbb *= dVbseff_dVb; // bug fixing
  }

  // copy over local drain (channel) current vars to instance vars:
  cdrain = Ids;
  gds = Gds;
  gm = Gm;
  gmbs = Gmb;

  // copy over local substrate current vars to instance vars:
  gbbs = Gbb;
  gbgs = Gbg;
  gbds = Gbd;

  csub = Isub;

  //  thermal noise Qinv calculated from all capMod
  //  * 0, 1, 2 & 3 stored in iterI->qinv 1/1998
  if ((model_.xpart < 0) || (!ChargeComputationNeeded))
  {
    qgate  = qdrn = qsrc = qbulk = 0.0;
    cggb = cgsb = cgdb = 0.0;
    cdgb = cdsb = cddb = 0.0;
    cbgb = cbsb = cbdb = 0.0;
    cqdb = cqsb = cqgb = cqbb = 0.0;

    gtau = 0.0;
    goto finished;
  }
  else if (model_.capMod == 0)
  {
    if (Vbseff < 0.0)
    {
      Vbseff = Vbs;
      dVbseff_dVb = 1.0;
    }
    else
    {
      Vbseff = paramPtr->phi - Phis;
      dVbseff_dVb = -dPhis_dVb;
    }

    Vfb = paramPtr->vfbcv;
    Vth = Vfb + paramPtr->phi + paramPtr->k1ox * sqrtPhis;
    Vgst = Vgs_eff - Vth;
    dVth_dVb = paramPtr->k1ox * dsqrtPhis_dVb;
    dVgst_dVb = -dVth_dVb;
    dVgst_dVg = dVgs_eff_dVg;

    CoxWL = model_.cox * paramPtr->weffCV * paramPtr->leffCV;
    Arg1 = Vgs_eff - Vbseff - Vfb;

    if (Arg1 <= 0.0)
    {
      qgate = CoxWL * Arg1;
      qbulk = -qgate;
      qdrn = 0.0;

      cggb = CoxWL * dVgs_eff_dVg;
      cgdb = 0.0;
      cgsb = CoxWL * (dVbseff_dVb - dVgs_eff_dVg);

      cdgb = 0.0;
      cddb = 0.0;
      cdsb = 0.0;

      cbgb = -CoxWL * dVgs_eff_dVg;
      cbdb = 0.0;
      cbsb = -cgsb;
      qinv = 0.0;
    }
    else if (Vgst <= 0.0)
    {
      T1 = 0.5 * paramPtr->k1ox;
      T2 = sqrt(T1 * T1 + Arg1);
      qgate = CoxWL * paramPtr->k1ox * (T2 - T1);
      qbulk = -qgate;
      qdrn = 0.0;

      T0 = CoxWL * T1 / T2;
      cggb = T0 * dVgs_eff_dVg;
      cgdb = 0.0;
      cgsb = T0 * (dVbseff_dVb - dVgs_eff_dVg);

      cdgb = 0.0;
      cddb = 0.0;
      cdsb = 0.0;

      cbgb = -cggb;
      cbdb = 0.0;
      cbsb = -cgsb;
      qinv = 0.0;
    }
    else
    {
      One_Third_CoxWL = CoxWL / 3.0;
      Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;

      Vdsat = Vgst / AbulkCV;
      dVdsat_dVg = dVgs_eff_dVg / AbulkCV;
      dVdsat_dVb = - (Vdsat * dAbulkCV_dVb + dVth_dVb)/ AbulkCV;

      if (model_.xpart > 0.5)
      { // 0/100 Charge partition model
        if (Vdsat <= Vds)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - T1);

          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.0;

          cggb = One_Third_CoxWL * (3.0 - dVdsat_dVg)* dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;

          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = 0.0;
          cddb = 0.0;
          cdsb = 0.0;

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);

          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          T7 = 2.0 * Vds - T1 - 3.0 * T3;
          T8 = T3 - T1 - 2.0 * Vds;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - 0.5 * (Vds - T3));

          T10 = T4 * T8;
          qdrn = T4 * T7;
          qbulk = -(qgate + qdrn + T10);

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;

          T11 = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + T11 + cgdb);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T7 = T9 * T7;
          T8 = T9 * T8;
          T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
          cdgb = (T7 * dAlphaz_dVg - T9* dVdsat_dVg) * dVgs_eff_dVg;

          T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
          cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
          cdsb = -(cdgb + T12 + cddb);

          T9  = 2.0 * T4 * (1.0 + T5);
          T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg) * dVgs_eff_dVg;
          T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
          T12 = T4 * (2.0 * T2 + T5 - 1.0);
          T0  = -(T10 + T11 + T12);


          cbgb = -(cggb + cdgb + T10);
          cbdb = -(cgdb + cddb + T12);
          cbsb = -(cgsb + cdsb + T0);
          qinv = -(qgate + qbulk);
        }
      }
      else if (model_.xpart < 0.5)
      {   // 40/60 Charge partition model
        if (Vds >= Vdsat)
        { // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - T1);

          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.4 * T2;

          cggb = One_Third_CoxWL* (3.0 - dVdsat_dVg) * dVgs_eff_dVg;

          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          T3 = 0.4 * Two_Third_CoxWL;
          cdgb = -T3 * dVgs_eff_dVg;
          cddb = 0.0;

          T4 = T3 * dVth_dVb;
          cdsb = -(T4 + cdgb);

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - 0.5 * (Vds - T3));

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;
          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds + 1.2 * Vds * Vds;
          T8 = T2 / T1;
          T7 = Vds - T1 - T8 * T6;
          qdrn = T4 * T7;
          T7 *= T9;
          tmp = T8 / T1;
          tmp1 = T4*(2.0 - 4.0 * tmp * T6 + T8 *(16.0 * Vdsat - 6.0 *Vds));

          cdgb = (T7 *dAlphaz_dVg - tmp1 *dVdsat_dVg) *dVgs_eff_dVg;

          T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
          cddb = T4 * (2.0 - (1.0 / (3.0 * T1
                                     * T1) + 2.0 * tmp) * T6 + T8
                       * (6.0 * Vdsat - 2.4 * Vds));

          cdsb = -(cdgb + T10 + cddb);

          T7 = 2.0 * (T1 + T3);
          qbulk = -(qgate - T4 * T7);
          T7 *= T9;
          T0 = 4.0 * T4 * (1.0 - T5);
          T12 = (-T7 * dAlphaz_dVg - cdgb
                 - T0 * dVdsat_dVg) * dVgs_eff_dVg;
          T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
          T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5)
            - cddb;

          tmp = -(T10 + T11 + T12);

          cbgb = -(cggb + cdgb + T12);
          cbdb = -(cgdb + cddb + T11);
          cbsb = -(cgsb + cdsb + tmp);
          qinv = -(qgate + qbulk);
        }
      }
      else
      {   // 50/50 partitioning
        if (Vds >= Vdsat)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.5 * T2;

          cggb = One_Third_CoxWL * (3.0 -dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
          cddb = 0.0;
          T4 = One_Third_CoxWL * dVth_dVb;
          cdsb = -(T4 + cdgb);

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
          qinv = -(qgate + qbulk);
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - 0.5 * (Vds - T3))
            ;

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;

          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T7 = T1 + T3;
          qdrn = -T4 * T7;
          qbulk = - (qgate + qdrn + qdrn);
          T7 *= T9;
          T0 = T4 * (2.0 * T5 - 2.0);

          cdgb = (T0 * dVdsat_dVg - T7 *dAlphaz_dVg) *dVgs_eff_dVg;
          T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
          cddb = T4 * (1.0 - 2.0 * T2 - T5);
          cdsb = -(cdgb + T12 + cddb);

          cbgb = -(cggb + 2.0 * cdgb);
          cbdb = -(cgdb + 2.0 * cddb);
          cbsb = -(cgsb + 2.0 * cdsb);
          qinv = -(qgate + qbulk);
        }
      }
    }
  }
  else
  {
    if (Vbseff < 0.0)
    {
      VbseffCV = Vbseff;
      dVbseffCV_dVb = 1.0;
    }
    else
    {
      VbseffCV = paramPtr->phi - Phis;
      dVbseffCV_dVb = -dPhis_dVb;
    }

    CoxWL = model_.cox * paramPtr->weffCV * paramPtr->leffCV;

    // Seperate VgsteffCV with noff and voffcv
    noff = n * paramPtr->noff;
    dnoff_dVd = paramPtr->noff * dn_dVd;
    dnoff_dVb = paramPtr->noff * dn_dVb;
    T0 = Vtm * noff;
    voffcv = paramPtr->voffcv;
    VgstNVt = (Vgst - voffcv) / T0;

    if (VgstNVt > CONSTEXP_THRESHOLD)
    {
      Vgsteff = Vgst - voffcv;
      dVgsteff_dVg = dVgs_eff_dVg;
      dVgsteff_dVd = -dVth_dVd;
      dVgsteff_dVb = -dVth_dVb;
    }
    else if (VgstNVt < -CONSTEXP_THRESHOLD)
    {
      Vgsteff = T0 * log(1.0 + CONSTMIN_EXP);
      dVgsteff_dVg = 0.0;
      dVgsteff_dVd = Vgsteff / noff;
      dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
      dVgsteff_dVd *= dnoff_dVd;
    }
    else
    {
      ExpVgst = exp(VgstNVt);
      Vgsteff = T0 * log(1.0 + ExpVgst);
      dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
      dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgst - voffcv)
                                      / noff * dnoff_dVd) + Vgsteff / noff * dnoff_dVd;
      dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgst - voffcv)
                                      / noff * dnoff_dVb) + Vgsteff / noff * dnoff_dVb;
      dVgsteff_dVg *= dVgs_eff_dVg;
    } // End of VgsteffCV

    if (model_.capMod == 1)
    {
      Vfb = paramPtr->vfbzb;
      Arg1 = Vgs_eff - VbseffCV - Vfb - Vgsteff;

      if (Arg1 <= 0.0)
      {
        qgate = CoxWL * Arg1;
        Cgg = CoxWL * (dVgs_eff_dVg - dVgsteff_dVg);
        Cgd = -CoxWL * dVgsteff_dVd;
        Cgb = -CoxWL * (dVbseffCV_dVb + dVgsteff_dVb);
      }
      else
      {
        T0 = 0.5 * paramPtr->k1ox;
        T1 = sqrt(T0 * T0 + Arg1);
        T2 = CoxWL * T0 / T1;

        qgate = CoxWL * paramPtr->k1ox * (T1 - T0);

        Cgg = T2 * (dVgs_eff_dVg - dVgsteff_dVg);
        Cgd = -T2 * dVgsteff_dVd;
        Cgb = -T2 * (dVbseffCV_dVb + dVgsteff_dVb);
      }
      qbulk = -qgate;
      Cbg = -Cgg;
      Cbd = -Cgd;
      Cbb = -Cgb;

      One_Third_CoxWL = CoxWL / 3.0;
      Two_Third_CoxWL = 2.0 * One_Third_CoxWL;
      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = Vgsteff / AbulkCV;

      if (VdsatCV < Vds)
      {
        dVdsatCV_dVg = 1.0 / AbulkCV;
        dVdsatCV_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
        T0 = Vgsteff - VdsatCV / 3.0;
        dT0_dVg = 1.0 - dVdsatCV_dVg / 3.0;
        dT0_dVb = -dVdsatCV_dVb / 3.0;
        qgate += CoxWL * T0;
        Cgg1 = CoxWL * dT0_dVg;
        Cgb1 = CoxWL * dT0_dVb + Cgg1 * dVgsteff_dVb;
        Cgd1 = Cgg1 * dVgsteff_dVd;
        Cgg1 *= dVgsteff_dVg;
        Cgg += Cgg1;
        Cgb += Cgb1;
        Cgd += Cgd1;

        T0 = VdsatCV - Vgsteff;
        dT0_dVg = dVdsatCV_dVg - 1.0;
        dT0_dVb = dVdsatCV_dVb;
        qbulk += One_Third_CoxWL * T0;
        Cbg1 = One_Third_CoxWL * dT0_dVg;
        Cbb1 = One_Third_CoxWL * dT0_dVb + Cbg1 * dVgsteff_dVb;
        Cbd1 = Cbg1 * dVgsteff_dVd;
        Cbg1 *= dVgsteff_dVg;
        Cbg += Cbg1;
        Cbb += Cbb1;
        Cbd += Cbd1;

        if (model_.xpart > 0.5)      T0 = -Two_Third_CoxWL;
        else if (model_.xpart < 0.5) T0 = -0.4 * CoxWL;
        else                         T0 = -One_Third_CoxWL;

        qsrc = T0 * Vgsteff;
        Csg = T0 * dVgsteff_dVg;
        Csb = T0 * dVgsteff_dVb;
        Csd = T0 * dVgsteff_dVd;
        Cgb *= dVbseff_dVb;
        Cbb *= dVbseff_dVb;
        Csb *= dVbseff_dVb;
      }
      else
      {
        T0 = AbulkCV * Vds;
        T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.e-20);
        T2 = Vds / T1;
        T3 = T0 * T2;
        dT3_dVg = -12.0 * T2 * T2 * AbulkCV;
        dT3_dVd = 6.0 * T0 * (4.0 * Vgsteff - T0) / T1 / T1 - 0.5;
        dT3_dVb = 12.0 * T2 * T2 * dAbulkCV_dVb * Vgsteff;

        qgate += CoxWL * (Vgsteff - 0.5 * Vds + T3);
        Cgg1 = CoxWL * (1.0 + dT3_dVg);
        Cgb1 = CoxWL * dT3_dVb + Cgg1 * dVgsteff_dVb;
        Cgd1 = CoxWL * dT3_dVd + Cgg1 * dVgsteff_dVd;
        Cgg1 *= dVgsteff_dVg;
        Cgg += Cgg1;
        Cgb += Cgb1;
        Cgd += Cgd1;

        qbulk += CoxWL * (1.0 - AbulkCV) * (0.5 * Vds - T3);
        Cbg1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVg);
        Cbb1 = -CoxWL * ((1.0 - AbulkCV) * dT3_dVb
                         + (0.5 * Vds - T3) * dAbulkCV_dVb)
          + Cbg1 * dVgsteff_dVb;
        Cbd1 = -CoxWL * (1.0 - AbulkCV) * dT3_dVd
          + Cbg1 * dVgsteff_dVd;
        Cbg1 *= dVgsteff_dVg;
        Cbg += Cbg1;
        Cbb += Cbb1;
        Cbd += Cbd1;

        if (model_.xpart > 0.5)
        {   // 0/100 Charge petition model
          T1 = T1 + T1;
          qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
          Csg = -CoxWL * (0.5 + 24.0 * T0 * Vds / T1 / T1 * AbulkCV);
          Csb = -CoxWL * (0.25 * Vds * dAbulkCV_dVb
                          - 12.0 * T0 * Vds / T1 / T1 * (4.0 * Vgsteff - T0)
                          * dAbulkCV_dVb) + Csg * dVgsteff_dVb;
          Csd = -CoxWL * (0.25 * AbulkCV - 12.0 * AbulkCV * T0
                          / T1 / T1 * (4.0 * Vgsteff - T0))
            + Csg * dVgsteff_dVd;
          Csg *= dVgsteff_dVg;
        }
        else if (model_.xpart < 0.5)
        {   // 40/60 Charge petition model
          T1 = T1 / 12.0;
          T2 = 0.5 * CoxWL / (T1 * T1);
          T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
                          * (Vgsteff - 4.0 * T0 / 3.0))
            - 2.0 * T0 * T0 * T0 / 15.0;
          qsrc = -T2 * T3;
          T4 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
            + 0.4 * T0 * T0;
          Csg = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                                                    * Vgsteff - 8.0 * T0 / 3.0)
                                         + 2.0 * T0 * T0 / 3.0);
          Csb = (qsrc / T1 * Vds + T2 * T4 * Vds) * dAbulkCV_dVb
            + Csg * dVgsteff_dVb;
          Csd = (qsrc / T1 + T2 * T4) * AbulkCV
            + Csg * dVgsteff_dVd;
          Csg *= dVgsteff_dVg;
        }
        else
        {   // 50/50 Charge petition model
          qsrc = -0.5 * (qgate + qbulk);
          Csg = -0.5 * (Cgg1 + Cbg1);
          Csb = -0.5 * (Cgb1 + Cbb1);
          Csd = -0.5 * (Cgd1 + Cbd1);
        }
        Cgb *= dVbseff_dVb;
        Cbb *= dVbseff_dVb;
        Csb *= dVbseff_dVb;
      }
      qdrn = -(qgate + qbulk + qsrc);
      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = -(qgate + qbulk);
    }
    else if (model_.capMod == 2)
    {
      Vfb = paramPtr->vfbzb;
      V3 = Vfb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (Vfb <= 0.0)
      {
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * Vfb);
        T2 = -CONSTDELTA_3 / T0;
      }
      else
      {
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * Vfb);
        T2 = CONSTDELTA_3 / T0;
      }

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = Vfb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;
      Qac0 = CoxWL * (Vfbeff - Vfb);
      dQac0_dVg = CoxWL * dVfbeff_dVg;
      dQac0_dVb = CoxWL * dVfbeff_dVb;

      T0 = 0.5 * paramPtr->k1ox;
      T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
      if (paramPtr->k1ox == 0.0)
      {
        T1 = 0.0;
        T2 = 0.0;
      }
      else if (T3 < 0.0)
      {
        T1 = T0 + T3 / paramPtr->k1ox;
        T2 = CoxWL;
      }
      else
      {
        T1 = sqrt(T0 * T0 + T3);
        T2 = CoxWL * T0 / T1;
      }

      Qsub0 = CoxWL * paramPtr->k1ox * (T1 - T0);

      dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
      dQsub0_dVd = -T2 * dVgsteff_dVd;
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb);

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = Vgsteff / AbulkCV;

      V4 = VdsatCV - Vds - CONSTDELTA_4;
      T0 = sqrt(V4 * V4 + 4.0 * CONSTDELTA_4 * VdsatCV);
      VdseffCV = VdsatCV - 0.5 * (V4 + T0);
      T1 = 0.5 * (1.0 + V4 / T0);
      T2 = CONSTDELTA_4 / T0;
      T3 = (1.0 - T1 - T2) / AbulkCV;
      dVdseffCV_dVg = T3;
      dVdseffCV_dVd = T1;
      dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;
      // Added to eliminate non-zero VdseffCV at Vds=0.0
      if (Vds == 0.0)
      {
        VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1e-20);
      T2 = VdseffCV / T1;
      T3 = T0 * T2;

      T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
      T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
      T6 = 12.0 * T2 * T2 * Vgsteff;

      qinoi = -CoxWL * (Vgsteff - 0.5 * T0 + AbulkCV * T3);
      qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
      Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
      Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
        + Cgg1 * dVgsteff_dVb;
      Cgg1 *= dVgsteff_dVg;

      T7 = 1.0 - AbulkCV;
      qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
      T4 = -T7 * (T4 - 1.0);
      T5 = -T7 * T5;
      T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
      Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
      Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
      Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
        + Cbg1 * dVgsteff_dVb;
      Cbg1 *= dVgsteff_dVg;

      if (model_.xpart > 0.5)
      {   // 0/100 Charge petition model
        T1 = T1 + T1;
        qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
        T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
        T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
        T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
        T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
        Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
        Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
        Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
          + Csg * dVgsteff_dVb;
        Csg *= dVgsteff_dVg;
      }
      else if (model_.xpart < 0.5)
      {   // 40/60 Charge petition model
        T1 = T1 / 12.0;
        T2 = 0.5 * CoxWL / (T1 * T1);
        T3 = Vgsteff *(2.0 * T0 *T0/3.0 +Vgsteff *(Vgsteff - 4.0 *T0/ 3.0))
          - 2.0 * T0 * T0 * T0 / 15.0;
        qsrc = -T2 * T3;
        T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
          + 0.4 * T0 * T0;
        T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                                                 * Vgsteff - 8.0 * T0 / 3.0)
                                      + 2.0 * T0 * T0 / 3.0);
        T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
        T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
        Csg = (T4 + T5 * dVdseffCV_dVg);
        Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
        Csb = (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
          + Csg * dVgsteff_dVb;
        Csg *= dVgsteff_dVg;
      }
      else
      {   // 50/50 Charge petition model
        qsrc = -0.5 * (qgate + qbulk);
        Csg = -0.5 * (Cgg1 + Cbg1);
        Csb = -0.5 * (Cgb1 + Cbb1);
        Csd = -0.5 * (Cgd1 + Cbd1);
      }

      qgate += Qac0 + Qsub0;
      qbulk -= (Qac0 + Qsub0);
      qdrn = -(qgate + qbulk + qsrc);

      Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
      Cgd = dQsub0_dVd + Cgd1;
      Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

      Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
      Cbd = Cbd1 - dQsub0_dVd;
      Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

      Cgb *= dVbseff_dVb;
      Cbb *= dVbseff_dVb;
      Csb *= dVbseff_dVb;

      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = qinoi;
    }

    // New Charge-Thickness capMod (CTM) begins
    else if (model_.capMod == 3)
    {
      V3 = paramPtr->vfbzb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (paramPtr->vfbzb <= 0.0)
      {
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * paramPtr->vfbzb);
        T2 = -CONSTDELTA_3 / T0;
      }
      else
      {
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * paramPtr->vfbzb);
        T2 = CONSTDELTA_3 / T0;
      }

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = paramPtr->vfbzb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;

      Cox = model_.cox;
      Tox = 1.0e8 * model_.tox;
      T0 = (Vgs_eff - VbseffCV - paramPtr->vfbzb) / Tox;
      dT0_dVg = dVgs_eff_dVg / Tox;
      dT0_dVb = -dVbseffCV_dVb / Tox;

      tmp = T0 * paramPtr->acde;
      if ((-CONSTEXP_THRESHOLD < tmp) && (tmp < CONSTEXP_THRESHOLD))
      {
        Tcen = paramPtr->ldeb * exp(tmp);
        dTcen_dVg = paramPtr->acde * Tcen;
        dTcen_dVb = dTcen_dVg * dT0_dVb;
        dTcen_dVg *= dT0_dVg;
      }
      else if (tmp <= -CONSTEXP_THRESHOLD)
      {
        Tcen = paramPtr->ldeb * CONSTMIN_EXP;
        dTcen_dVg = dTcen_dVb = 0.0;
      }
      else
      {
        Tcen = paramPtr->ldeb * CONSTMAX_EXP;
        dTcen_dVg = dTcen_dVb = 0.0;
      }

      LINK = 1.0e-3 * model_.tox;
      V3 = paramPtr->ldeb - Tcen - LINK;
      V4 = sqrt(V3 * V3 + 4.0 * LINK * paramPtr->ldeb);
      Tcen = paramPtr->ldeb - 0.5 * (V3 + V4);
      T1 = 0.5 * (1.0 + V3 / V4);
      dTcen_dVg *= T1;
      dTcen_dVb *= T1;

      Ccen = CONSTEPSSI / Tcen;
      T2 = Cox / (Cox + Ccen);
      Coxeff = T2 * Ccen;
      T3 = -Ccen / Tcen;
      dCoxeff_dVg = T2 * T2 * T3;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / Cox;

      Qac0 = CoxWLcen * (Vfbeff - paramPtr->vfbzb);
      QovCox = Qac0 / Coxeff;
      dQac0_dVg = CoxWLcen * dVfbeff_dVg + QovCox * dCoxeff_dVg;
      dQac0_dVb = CoxWLcen * dVfbeff_dVb + QovCox * dCoxeff_dVb;

      T0 = 0.5 * paramPtr->k1ox;
      T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
      if (paramPtr->k1ox == 0.0)
      {
        T1 = 0.0;
        T2 = 0.0;
      }
      else if (T3 < 0.0)
      {
        T1 = T0 + T3 / paramPtr->k1ox;
        T2 = CoxWLcen;
      }
      else
      {
        T1 = sqrt(T0 * T0 + T3);
        T2 = CoxWLcen * T0 / T1;
      }

      Qsub0 = CoxWLcen * paramPtr->k1ox * (T1 - T0);
      QovCox = Qsub0 / Coxeff;
      dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
        + QovCox * dCoxeff_dVg;
      dQsub0_dVd = -T2 * dVgsteff_dVd;
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
        + QovCox * dCoxeff_dVb;

      // Gate-bias dependent delta Phis begins */
      if (paramPtr->k1ox <= 0.0)
      {
        Denomi = 0.25 * paramPtr->moin * Vtm;
        T0 = 0.5 * paramPtr->sqrtPhi;
      }
      else
      {
        Denomi = paramPtr->moin * Vtm * paramPtr->k1ox * paramPtr->k1ox;
        T0 = paramPtr->k1ox * paramPtr->sqrtPhi;
      }
      T1 = 2.0 * T0 + Vgsteff;

      DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff / Denomi);
      dDeltaPhi_dVg = 2.0 * Vtm * (T1 -T0) / (Denomi + T1 * Vgsteff);
      dDeltaPhi_dVd = dDeltaPhi_dVg * dVgsteff_dVd;
      dDeltaPhi_dVb = dDeltaPhi_dVg * dVgsteff_dVb;
      // End of delta Phis

      T3 = 4.0 * (Vth - paramPtr->vfbzb - paramPtr->phi);
      Tox += Tox;
      if (T3 >= 0.0)
      {  T0 = (Vgsteff + T3) / Tox;
      dT0_dVd = (dVgsteff_dVd + 4.0 * dVth_dVd) / Tox;
      dT0_dVb = (dVgsteff_dVb + 4.0 * dVth_dVb) / Tox;
      }
      else
      {  T0 = (Vgsteff + 1.0e-20) / Tox;
      dT0_dVd = dVgsteff_dVd / Tox;
      dT0_dVb = dVgsteff_dVb / Tox;
      }
      tmp = exp(0.7 * log(T0));
      T1 = 1.0 + tmp;
      T2 = 0.7 * tmp / (T0 * Tox);
      Tcen = 1.9e-9 / T1;
      dTcen_dVg = -1.9e-9 * T2 / T1 /T1;
      dTcen_dVd = Tox * dTcen_dVg;
      dTcen_dVb = dTcen_dVd * dT0_dVb;
      dTcen_dVd *= dT0_dVd;
      dTcen_dVg *= dVgsteff_dVg;

      Ccen = CONSTEPSSI / Tcen;
      T0 = Cox / (Cox + Ccen);
      Coxeff = T0 * Ccen;
      T1 = -Ccen / Tcen;
      dCoxeff_dVg = T0 * T0 * T1;
      dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / Cox;

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = (Vgsteff - DeltaPhi) / AbulkCV;
      V4 = VdsatCV - Vds - CONSTDELTA_4;
      T0 = sqrt(V4 * V4 + 4.0 * CONSTDELTA_4 * VdsatCV);
      VdseffCV = VdsatCV - 0.5 * (V4 + T0);
      T1 = 0.5 * (1.0 + V4 / T0);
      T2 = CONSTDELTA_4 / T0;
      T3 = (1.0 - T1 - T2) / AbulkCV;
      T4 = T3 * ( 1.0 - dDeltaPhi_dVg);
      dVdseffCV_dVg = T4;
      dVdseffCV_dVd = T1;
      dVdseffCV_dVb = -T3 * VdsatCV * dAbulkCV_dVb;

      // Added to eliminate non-zero VdseffCV at Vds=0.0
      if (Vds == 0.0)
      {  VdseffCV = 0.0;
      dVdseffCV_dVg = 0.0;
      dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = Vgsteff - DeltaPhi;
      T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
      T3 = T0 / T2;
      T4 = 1.0 - 12.0 * T3 * T3;
      T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
      T6 = T5 * VdseffCV / AbulkCV;

      qgate = qinoi = CoxWLcen * (T1 - T0 * (0.5 - T3));
      QovCox = qgate / Coxeff;
      Cgg1 = CoxWLcen * (T4 * (1.0 - dDeltaPhi_dVg) + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
        * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
        + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


      T7 = 1.0 - AbulkCV;
      T8 = T2 * T2;
      T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
      T10 = T9 * (1.0 - dDeltaPhi_dVg);
      T11 = -T7 * T5 / AbulkCV;
      T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

      qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
      QovCox = qbulk / Coxeff;
      Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
      Cbd1 = CoxWLcen * T11 * dVdseffCV_dVd + Cbg1
        * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
        + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

      if (model_.xpart > 0.5)
      {   // 0/100 partition
        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0 - 0.5 * T0 * T0 / T2);
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

      qgate += Qac0 + Qsub0 - qbulk;
      qbulk -= (Qac0 + Qsub0);
      qdrn = -(qgate + qbulk + qsrc);

      Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
      Cbd = Cbd1 - dQsub0_dVd;
      Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

      Cgg = Cgg1 - Cbg;
      Cgd = Cgd1 - Cbd;
      Cgb = Cgb1 - Cbb;

      Cgb *= dVbseff_dVb;
      Cbb *= dVbseff_dVb;
      Csb *= dVbseff_dVb;

      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
      qinv = -qinoi;
    }  // End of CTM
  }

  finished:
  // Returning Values to Calling Routine
  // COMPUTE EQUIVALENT DRAIN CURRENT SOURCE

  // copy local "cdrain" variable over to instance variable cd.
  cd = cdrain;

  // charge storage elements:
  // bulk-drain and bulk-source depletion capacitances
  // czbd : zero bias drain junction capacitance
  // czbs : zero bias source junction capacitance
  // czbdsw: zero bias drain junction sidewall capacitance
  //         along field oxide
  // czbssw: zero bias source junction sidewall capacitance
  //         along field oxide
  // czbdswg: zero bias drain junction sidewall capacitance
  //          along gate side
  // czbsswg: zero bias source junction sidewall capacitance
  //          along gate side
  if (ChargeComputationNeeded)
  {
    czbd = unitAreaJctCapTemp * drainArea;
    czbs = unitAreaJctCapTemp * sourceArea;
    if (drainPerimeter < paramPtr->weff)
    {
      czbdswg = unitLengthGateSidewallJctCapTemp * drainPerimeter;
      czbdsw = 0.0;
    }
    else
    {
      czbdsw = unitLengthSidewallJctCapTemp
        * (drainPerimeter - paramPtr->weff);

      czbdswg = unitLengthGateSidewallJctCapTemp *  paramPtr->weff;
    }
    if (sourcePerimeter < paramPtr->weff)
    {
      czbssw = 0.0;
      czbsswg = unitLengthGateSidewallJctCapTemp * sourcePerimeter;
    }
    else
    {
      czbssw = unitLengthSidewallJctCapTemp
        * (sourcePerimeter - paramPtr->weff);
      czbsswg = unitLengthGateSidewallJctCapTemp *  paramPtr->weff;
    }

    MJ    = model_.bulkJctBotGradingCoeff;
    MJSW  = model_.bulkJctSideGradingCoeff;
    MJSWG = model_.bulkJctGateSideGradingCoeff;

    // Source Bulk Junction
    if (vbs == 0.0)
    {
      //*(ckt->CKTstate0 + iterI->qbs) = 0.0;
      qbs = 0.0;
      capbs = czbs + czbssw + czbsswg;
    }
    else if (vbs < 0.0)
    {
      if (czbs > 0.0)
      {
        arg = 1.0 - vbs / PhiBTemp;

        if (MJ == 0.5) sarg = 1.0 / sqrt(arg);
        else           sarg = exp(-MJ * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) =
        qbs = PhiBTemp * czbs * (1.0 - arg * sarg) / (1.0 - MJ);

        capbs = czbs * sarg;
      }
      else
      {
        //*(ckt->CKTstate0 + iterI->qbs) = 0.0;
        qbs = 0.0;
        capbs = 0.0;
      }

      if (czbssw > 0.0)
      {
        arg = 1.0 - vbs / PhiBSWTemp;
        if (MJSW == 0.5) sarg = 1.0 / sqrt(arg);
        else             sarg = exp(-MJSW * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) +=
        qbs += PhiBSWTemp * czbssw * (1.0 - arg * sarg) / (1.0 - MJSW);

        capbs += czbssw * sarg;
      }

      if (czbsswg > 0.0)
      {
        arg = 1.0 - vbs / PhiBSWGTemp;
        if (MJSWG == 0.5) sarg = 1.0 / sqrt(arg);
        else              sarg = exp(-MJSWG * log(arg));

        //*(ckt->CKTstate0 + iterI->qbs) +=
        qbs += PhiBSWGTemp * czbsswg * (1.0 - arg * sarg) / (1.0 - MJSWG);

        capbs += czbsswg * sarg;
      }

    }
    else
    {
      T0 = czbs + czbssw + czbsswg;
      T1 = vbs * (czbs * MJ / PhiBTemp + czbssw * MJSW
                  / PhiBSWTemp + czbsswg * MJSWG / PhiBSWGTemp);

      //*(ckt->CKTstate0 + iterI->
      qbs = vbs * (T0 + 0.5 * T1);
      capbs = T0 + T1;
    }

    // Drain Bulk Junction
    if (vbd == 0.0)
    {
      //*(ckt->CKTstate0 + iterI->qbd) = 0.0;
      qbd = 0.0;
      capbd = czbd + czbdsw + czbdswg;
    }
    else if (vbd < 0.0)
    {
      if (czbd > 0.0)
      {
        arg = 1.0 - vbd / PhiBTemp;
        if (MJ == 0.5) sarg = 1.0 / sqrt(arg);
        else           sarg = exp(-MJ * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) =
        qbd = PhiBTemp * czbd * (1.0 - arg * sarg) / (1.0 - MJ);
        capbd = czbd * sarg;
      }
      else
      {
        //*(ckt->CKTstate0 + iterI->qbd) = 0.0;
        qbd = 0.0;
        capbd = 0.0;
      }

      if (czbdsw > 0.0)
      {
        arg = 1.0 - vbd / PhiBSWTemp;
        if (MJSW == 0.5) sarg = 1.0 / sqrt(arg);
        else             sarg = exp(-MJSW * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) +=
        qbd += PhiBSWTemp * czbdsw * (1.0 - arg * sarg) / (1.0 - MJSW);
        capbd += czbdsw * sarg;
      }

      if (czbdswg > 0.0)
      {
        arg = 1.0 - vbd / PhiBSWGTemp;
        if (MJSWG == 0.5) sarg = 1.0 / sqrt(arg);
        else              sarg = exp(-MJSWG * log(arg));

        //*(ckt->CKTstate0 + iterI->qbd) +=
        qbd += PhiBSWGTemp * czbdswg * (1.0 - arg * sarg) / (1.0 - MJSWG);
        capbd += czbdswg * sarg;
      }
    }
    else
    {
      T0 = czbd + czbdsw + czbdswg;
      T1 = vbd * (czbd * MJ / PhiBTemp + czbdsw * MJSW
                  / PhiBSWTemp + czbdswg * MJSWG / PhiBSWGTemp);

      //*(ckt->CKTstate0 + iterI->qbd) = vbd * (T0 + 0.5 * T1);
      qbd = vbd * (T0 + 0.5 * T1);
      capbd = T0 + T1;
    }
  }

  // There is a spice3f5 convergence check that would happen here.
  // (line 2404) skipping...

  // In 3f5, loading a bunch of things into the state vector at this point.
  // (line 2433) skipping...

  // bulk and channel charge plus overlaps

  if (ChargeComputationNeeded)
  {
    // NQS begins
    if (nqsMod)
    {
      qcheq = -(qbulk + qgate);

      cqgb = -(cggb + cbgb);
      cqdb = -(cgdb + cbdb);
      cqsb = -(cgsb + cbsb);
      cqbb = -(cqgb + cqdb + cqsb);

      gtau_drift = fabs(paramPtr->tconst * qcheq) * ScalingFactor;
      T0 = paramPtr->leffCV * paramPtr->leffCV;
      gtau_diff = 16.0 * paramPtr->u0temp * model_.vtm / T0 * ScalingFactor;

      gtau =  gtau_drift + gtau_diff;
    }

    if (model_.capMod == 0)
    {
      if (vgd < 0.0)
      {
        cgdo = paramPtr->cgdo;
        qgdo = paramPtr->cgdo * vgd;
      }
      else
      {
        cgdo = paramPtr->cgdo;
        qgdo =  paramPtr->cgdo * vgd;
      }

      if (vgs < 0.0)
      {
        cgso = paramPtr->cgso;
        qgso = paramPtr->cgso * vgs;
      }
      else
      {
        cgso = paramPtr->cgso;
        qgso =  paramPtr->cgso * vgs;
      }
    }
    else if (model_.capMod == 1)
    {
      if (vgd < 0.0)
      {
        T1 = sqrt(1.0 - 4.0 * vgd / paramPtr->ckappa);
        cgdo = paramPtr->cgdo + paramPtr->weffCV * paramPtr->cgdl / T1;

        qgdo = paramPtr->cgdo * vgd - paramPtr->weffCV * 0.5
          * paramPtr->cgdl * paramPtr->ckappa * (T1 - 1.0);
      }
      else
      {
        cgdo = paramPtr->cgdo + paramPtr->weffCV * paramPtr->cgdl;
        qgdo = (paramPtr->weffCV * paramPtr->cgdl + paramPtr->cgdo) * vgd;
      }

      if (vgs < 0.0)
      {
        T1 = sqrt(1.0 - 4.0 * vgs / paramPtr->ckappa);
        cgso = paramPtr->cgso + paramPtr->weffCV * paramPtr->cgsl / T1;
        qgso = paramPtr->cgso * vgs - paramPtr->weffCV * 0.5
          * paramPtr->cgsl * paramPtr->ckappa * (T1 - 1.0);
      }
      else
      {
        cgso = paramPtr->cgso + paramPtr->weffCV * paramPtr->cgsl;
        qgso = (paramPtr->weffCV * paramPtr->cgsl + paramPtr->cgso) * vgs;
      }
    }
    else
    {
      T0 = vgd + CONSTDELTA_1;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
      T2 = 0.5 * (T0 - T1);

      T3 = paramPtr->weffCV * paramPtr->cgdl;
      T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappa);
      cgdo = paramPtr->cgdo + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);

      qgdo = (paramPtr->cgdo + T3) * vgd - T3 * (T2
                                 + 0.5 * paramPtr->ckappa * (T4 - 1.0));

      T0 = vgs + CONSTDELTA_1;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
      T2 = 0.5 * (T0 - T1);
      T3 = paramPtr->weffCV * paramPtr->cgsl;
      T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappa);
      cgso = paramPtr->cgso + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);
      qgso = (paramPtr->cgso + T3) * vgs - T3 * (T2
                                 + 0.5 * paramPtr->ckappa * (T4 - 1.0));
    }

    //cgdo = cgdo;
    //cgso = cgso;

    setupCapacitors_oldDAE ();
    setupCapacitors_newDAE ();

    //cqdef = cqcheq = 0.0;
    cqdef = 0.0;

    // set some state variables:
    qg = qgate;
    qd = qdrn - qbd;
    qb = qbulk + qbd + qbs;
    if (nqsMod) qcdump = qdef * ScalingFactor;

  } // end of ChargeComputationNeeded if statement.

  // store small signal parameters
  //if (ckt->CKTmode & MODEINITSMSIG) goto line1000;
  // Note: in 3f5, line1000 is at the end of the load, after
  //        the loads to the rhs and the matrix.  So it looks
  //        like this goto essentially means return.


  // Setting up a few currents for the RHS load:
  Idrain = drainConductance * Vddp;
  Isource = sourceConductance * Vssp;

  // Put this kludge in because the matrix load needs T1 but it is used
  // all over the place:
  T1global = T1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
//
// Purpose       : This function sets up the primaray state variables into
//                 the primary state vector.
//
//                 These variables include qb, qg, qd and, in the
//                 event that nqsMod=1, qcdump and qcheq.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/09/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  double * staVec = extData.nextStaVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;

  bsuccess = updateIntermediateVars ();

  // voltage drops:
  stoVec[li_store_vbs] = vbs;
  stoVec[li_store_vgs] = vgs;
  stoVec[li_store_vds] = vds;
  stoVec[li_store_vbd] = vbd;
  stoVec[li_store_von] = von;

  // intrinsic capacitors:
  staVec[li_state_qb] = qb;
  staVec[li_state_qg] = qg;
  staVec[li_state_qd] = qd;

  // parasitic capacitors:
  staVec[li_state_qbs] = qbs;
  staVec[li_state_qbd] = qbd;

  if( nqsMod )
  {
    staVec[li_state_qcheq] = qcheq;
    staVec[li_state_qcdump] = qcdump;
  }

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

    // parasitic capacitors:
    currStaVec[li_state_qbs] = qbs;
    currStaVec[li_state_qbd] = qbd;

    if( nqsMod )
    {
      currStaVec[li_state_qcheq] = qcheq;
      currStaVec[li_state_qcdump] = qcdump;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 bsim3 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * dQdxdVp = extData.dQdxdVpVectorRawPtr;

  auxChargeCalculations ();

  double Qeqqg = 0.0;   // gate charge
  double Qeqqb = 0.0;   // bulk charge
  double Qeqqd = 0.0;   // drain charge
  double Qqdef = 0.0;   // nqs-related charge.
  double Qqcheq = 0.0;  // nqs-related charge.

  // These 3 vars are class variables, and are set up elsewhere.
  //double Qeqqg_Jdxp = 0.0; // limiter, related to gate  cap.
  //double Qeqqb_Jdxp = 0.0; // limiter, related to bulk  cap.
  //double Qeqqd_Jdxp = 0.0; // limiter, related to drain cap.

  if (model_.dtype > 0)
  {
    Qeqqg  = qg;
    Qeqqb  = qb;
    Qeqqd  = qd;
    Qqdef  = qcdump; // this needs to be fixed...
    Qqcheq = qcheq;
  }
  else  // need to convert these to charges.
  {
    Qeqqg  = -qg;
    Qeqqb  = -qb;
    Qeqqd  = -qd;
    Qqdef  = -qcdump;
    Qqcheq = -qcheq;
  }

  qVec[li_Gate] += Qeqqg*numberParallel;
  qVec[li_Bulk] += (Qeqqb)*numberParallel;
  qVec[li_DrainPrime] += (-(-Qeqqd))*numberParallel;
  qVec[li_SourcePrime] += (-(+ Qeqqg + Qeqqb + Qeqqd))*numberParallel;

  if( loadLeadCurrent )
  {
    double * storeLeadQ = extData.storeLeadCurrQCompRawPtr;
    if (drainConductance == 0.0)
    {
      storeLeadQ[li_store_dev_id] = (-(-Qeqqd))*numberParallel;
    }
    if (sourceConductance == 0.0)
    {
      storeLeadQ[li_store_dev_is] = (-(Qeqqg + Qeqqb + Qeqqd))*numberParallel;
    }
    storeLeadQ[li_store_dev_ig] = Qeqqg*numberParallel;
    storeLeadQ[li_store_dev_ib] = (Qeqqb)*numberParallel;
  }

  if (nqsMod)
  {
    // 7 equ. for nqs modification. charge equation.
    qVec[li_Charge] += -(Qqcheq - Qqdef)*numberParallel;
  }

  //////////////////////////////////////////////////
  // limiting section:
  if (devOptions.voltageLimiterFlag)
  {
    // Need the following:
    //  Qeqqg_Jdxp
    //  Qeqqb_Jdxp
    //  Qeqqd_Jdxp
    if (model_.dtype > 0)
    {
      // no-op:
      Qeqqg_Jdxp = Qeqqg_Jdxp;
      Qeqqb_Jdxp = Qeqqb_Jdxp;
      Qeqqd_Jdxp = Qeqqd_Jdxp;
    }
    else
    {
      Qeqqg_Jdxp = -Qeqqg_Jdxp;
      Qeqqb_Jdxp = -Qeqqb_Jdxp;
      Qeqqd_Jdxp = -Qeqqd_Jdxp;
    }

    if (!origFlag)
    {
      dQdxdVp[li_Gate] += -Qeqqg_Jdxp*numberParallel;
      dQdxdVp[li_Bulk] += -(+Qeqqb_Jdxp)*numberParallel;
      dQdxdVp[li_DrainPrime] += (-Qeqqd_Jdxp)*numberParallel;
      dQdxdVp[li_SourcePrime] += (+Qeqqg_Jdxp+Qeqqb_Jdxp+Qeqqd_Jdxp) *numberParallel;
    } // orig flag.
  } // limiter flag

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::auxChargeCalculations
//
// Purpose       : This function does some final "cleanup" calculations
//                 having to do with the capacitors.  
//
// Special Notes : About all this function really does is set up some
//                 voltlim terms, and some unused nqs stuff.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/04
//-----------------------------------------------------------------------------
bool Instance::auxChargeCalculations ()
{
  double T0, T1;

  if (!ChargeComputationNeeded)
  {
    sxpart = (1.0 - (dxpart = (mode > 0) ? 0.4 : 0.6));
    ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
    dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;

    if (nqsMod)
    {
      gtau = 16.0 * paramPtr->u0temp * model_.vtm
        / paramPtr->leffCV / paramPtr->leffCV * ScalingFactor;
    }
    else
    {
      gtau = 0.0;
    }
  }
  else  // ChargeComputation is needed
  {
    double vgb_orig = vgs_orig - vbs_orig;

    Qeqqg_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqg_Jdxp = - CAPcggb * (vgb-vgb_orig)
                   + CAPcgdb * (vbd-vbd_orig)
                   + CAPcgsb * (vbs-vbs_orig);
    }

    Qeqqb_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqb_Jdxp = - CAPcbgb * (vgb-vgb_orig)
                   + CAPcbdb * (vbd-vbd_orig)
                   + CAPcbsb * (vbs-vbs_orig);
    }

    Qeqqd_Jdxp = 0.0;
    if (!origFlag)
    {
      Qeqqd_Jdxp = - CAPcdgb * (vgb-vgb_orig)
                   + CAPcddb * (vbd-vbd_orig)
                   + CAPcdsb * (vbs-vbs_orig);
    }

    // Note:  nqs stuff is not yet finished, especially the voltage
    // limiting aspect.  For limiting, need to re-do T0 and the term
    // added to ceqqd.
    if (nqsMod)
    {
      string msg;
      msg = "Instance::auxChargeCalculations ()";
      msg += " nqsMod=1 is not ready yet.  Re-run with nqsMod=0\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);

      T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
      ceqqg += T0;
      T1 = qdef * gtau;
      ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd
                                   * vbd - ddxpart_dVs * vbs);

      //cqdef = cqcdump - gqdef * qdef;

      //if (!origFlag)
      //{
        //double tmp =   - (gcqgb * (vgb-vgb_orig)
                       //- gcqdb * (vbd-vbd_orig)
                       //- gcqsb * (vbs-vbs_orig)) + T0;
        //cqcheq += tmp;
        //cqcheq_Jdxp = tmp;
      //}
    }
  } // !ChargeComputationNeeded

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupCapacitors_newDAE ()
//
// Purpose       : This takes a lot of the individual capacitive terms and
//                 sums them together for loading into the Jacobian.
//
// Special Notes : This was extracted from updateIntermediateVars.
//                 Different variables are used, and nothing is multiplied
//                 by ag0 = solState.pdt.  The new-dae formulation handles
//                 all the 1/dt - related stuff up in the time integrator.
//
//                 NOTE: nqs not even close to being supported.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/23/04
//-----------------------------------------------------------------------------
bool Instance::setupCapacitors_newDAE ()
{

  if (mode > 0)
  {
    if (nqsMod == 0)
    {
      CAPcggb = (cggb + cgdo + cgso + paramPtr->cgbo );
      CAPcgdb = (cgdb - cgdo);
      CAPcgsb = (cgsb - cgso);

      CAPcdgb = (cdgb - cgdo);
      CAPcddb = (cddb + capbd + cgdo);
      CAPcdsb = cdsb;

      CAPcsgb = -(cggb + cbgb + cdgb + cgso);
      CAPcsdb = -(cgdb + cbdb + cddb);
      CAPcssb = (capbs + cgso - (cgsb + cbsb + cdsb));

      CAPcbgb = (cbgb - paramPtr->cgbo);
      CAPcbdb = (cbdb - capbd);
      CAPcbsb = (cbsb - capbs);
    }
    else  // nqsMode != 0
    {

    } // nqsMod
  }
  else
  {
    if (nqsMod == 0)
    {
      CAPcggb = (cggb + cgdo + cgso + paramPtr->cgbo );
      CAPcgdb = (cgsb - cgdo);
      CAPcgsb = (cgdb - cgso);

      CAPcdgb = -(cggb + cbgb + cdgb + cgdo);
      CAPcddb = (capbd + cgdo - (cgsb + cbsb + cdsb));
      CAPcdsb = -(cgdb + cbdb + cddb);

      CAPcsgb = (cdgb - cgso);
      CAPcsdb = cdsb;
      CAPcssb = (cddb + capbs + cgso);

      CAPcbgb = (cbgb - paramPtr->cgbo);
      CAPcbdb = (cbsb - capbd);
      CAPcbsb = (cbdb - capbs);
    }
    else  // nqsMode != 0
    {

    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupCapacitors_oldDAE ()
// Purpose       : Same as new-DAE version, but including pdt, essentially.
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/17/06
//-----------------------------------------------------------------------------
bool Instance::setupCapacitors_oldDAE ()
{

  double ag0 = solState.pdt;
  double T0 = 0.0;

  // ERK. 12/17/2006.
  // It is necessary to set ag0=0.0, because for the first time step out of
  // the DCOP, all the time derivatives are forced to be zero.  Thus, all
  // their derivatives should also be zero.  If it wasn't for that, then ag0
  // could always be pdt.  (it used to be, before the -jacobian_test capability).
  if (!(solState.dcopFlag) && solState.initTranFlag && solState.newtonIter==0)
  {
    ag0 = 0.0;
  }

  if (mode > 0)
  {
    if (nqsMod == 0)
    {
      gcggb = (cggb + cgdo + cgso + paramPtr->cgbo ) * ag0;
      gcgdb = (cgdb - cgdo) * ag0;
      gcgsb = (cgsb - cgso) * ag0;

      gcdgb = (cdgb - cgdo) * ag0;
      gcddb = (cddb + capbd + cgdo) * ag0;
      gcdsb = cdsb * ag0;

      gcsgb = -(cggb + cbgb + cdgb + cgso) * ag0;
      gcsdb = -(cgdb + cbdb + cddb) * ag0;
      gcssb = (capbs + cgso - (cgsb + cbsb + cdsb)) * ag0;

      gcbgb = (cbgb - paramPtr->cgbo) * ag0;
      gcbdb = (cbdb - capbd) * ag0;
      gcbsb = (cbsb - capbs) * ag0;

      qgd = qgdo;
      qgs = qgso;
      qgb = paramPtr->cgbo * vgb;
      qgate += qgd + qgs + qgb;
      qbulk -= qgb;
      qdrn -= qgd;
      qsrc = -(qgate + qbulk + qdrn);

      ggtg = ggtd = ggtb = ggts = 0.0;
      sxpart = 0.6;
      dxpart = 0.4;
      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    }
    else  // nqsMode != 0
    {
      if (qcheq > 0.0)
      {
        T0 = paramPtr->tconst * qdef * ScalingFactor;
      }
      else
      {
        T0 = -paramPtr->tconst * qdef * ScalingFactor;
      }

      ggtg = gtg = T0 * cqgb;
      ggtd = gtd = T0 * cqdb;
      ggts = gts = T0 * cqsb;
      ggtb = gtb = T0 * cqbb;
      gqdef = ScalingFactor * ag0;

      gcqgb = cqgb * ag0;
      gcqdb = cqdb * ag0;
      gcqsb = cqsb * ag0;
      gcqbb = cqbb * ag0;

      gcggb = (cgdo + cgso + paramPtr->cgbo ) * ag0;
      gcgdb = -cgdo * ag0;
      gcgsb = -cgso * ag0;

      gcdgb = -cgdo * ag0;
      gcddb = (capbd + cgdo) * ag0;
      gcdsb = 0.0;

      gcsgb = -cgso * ag0;
      gcsdb = 0.0;
      gcssb = (capbs + cgso) * ag0;

      gcbgb = -paramPtr->cgbo * ag0;
      gcbdb = -capbd * ag0;
      gcbsb = -capbs * ag0;

      CoxWL = model_.cox * paramPtr->weffCV * paramPtr->leffCV;

      if (fabs(qcheq) <= 1.0e-5 * CoxWL)
      {
        if      (model_.xpart < 0.5) dxpart = 0.4;
        else if (model_.xpart > 0.5) dxpart = 0.0;
        else                               dxpart = 0.5;

        ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      }
      else
      {
        dxpart = qdrn / qcheq;
        Cdd = cddb;
        Csd = -(cgdb + cddb + cbdb);
        ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
        Cdg = cdgb;
        Csg = -(cggb + cdgb + cbgb);
        ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

        Cds = cdsb;
        Css = -(cgsb + cdsb + cbsb);
        ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

        ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
      }
      sxpart = 1.0 - dxpart;
      dsxpart_dVd = -ddxpart_dVd;
      dsxpart_dVg = -ddxpart_dVg;
      dsxpart_dVs = -ddxpart_dVs;
      dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

      qgd = qgdo;
      qgs = qgso;
      qgb = paramPtr->cgbo * vgb;
      qgate = qgd + qgs + qgb;
      qbulk = -qgb;
      qdrn = -qgd;
      qsrc = -(qgate + qbulk + qdrn);
    } // nqsMod
  }
  else
  {
    if (nqsMod == 0)
    {
      gcggb = (cggb + cgdo + cgso + paramPtr->cgbo ) * ag0;
      gcgdb = (cgsb - cgdo) * ag0;
      gcgsb = (cgdb - cgso) * ag0;

      gcdgb = -(cggb + cbgb + cdgb + cgdo) * ag0;
      gcddb = (capbd + cgdo - (cgsb + cbsb + cdsb)) * ag0;
      gcdsb = -(cgdb + cbdb + cddb) * ag0;

      gcsgb = (cdgb - cgso) * ag0;
      gcsdb = cdsb * ag0;
      gcssb = (cddb + capbs + cgso) * ag0;

      gcbgb = (cbgb - paramPtr->cgbo) * ag0;
      gcbdb = (cbsb - capbd) * ag0;
      gcbsb = (cbdb - capbs) * ag0;

      qgd = qgdo;
      qgs = qgso;
      qgb = paramPtr->cgbo * vgb;
      qgate += qgd + qgs + qgb;
      qbulk -= qgb;
      qsrc = qdrn - qgs;
      qdrn = -(qgate + qbulk + qsrc);

      ggtg = ggtd = ggtb = ggts = 0.0;
      sxpart = 0.4;
      dxpart = 0.6;
      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
    }
    else  // nqsMode != 0
    {
      if (qcheq > 0.0)
      {
        T0 = paramPtr->tconst * qdef * ScalingFactor;
      }
      else
      {
        T0 = -paramPtr->tconst * qdef * ScalingFactor;
      }

      ggtg = gtg = T0 * cqgb;
      ggts = gtd = T0 * cqdb;
      ggtd = gts = T0 * cqsb;
      ggtb = gtb = T0 * cqbb;
      gqdef = ScalingFactor * ag0;

      gcqgb = cqgb * ag0;
      gcqdb = cqsb * ag0;
      gcqsb = cqdb * ag0;
      gcqbb = cqbb * ag0;

      gcggb = (cgdo + cgso + paramPtr->cgbo) * ag0;
      gcgdb = -cgdo * ag0;
      gcgsb = -cgso * ag0;

      gcdgb = -cgdo * ag0;
      gcddb = (capbd + cgdo) * ag0;
      gcdsb = 0.0;

      gcsgb = -cgso * ag0;
      gcsdb = 0.0;
      gcssb = (capbs + cgso) * ag0;

      gcbgb = -paramPtr->cgbo * ag0;
      gcbdb = -capbd * ag0;
      gcbsb = -capbs * ag0;

      CoxWL = model_.cox * paramPtr->weffCV * paramPtr->leffCV;

      if (fabs(qcheq) <= 1.0e-5 * CoxWL)
      {
        if      (model_.xpart < 0.5) sxpart = 0.4;
        else if (model_.xpart > 0.5) sxpart = 0.0;
        else                         sxpart = 0.5;

        dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
      }
      else
      {
        sxpart = qdrn / qcheq;
        Css = cddb;
        Cds = -(cgdb + cddb + cbdb);
        dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
        Csg = cdgb;
        Cdg = -(cggb + cdgb + cbgb);
        dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

        Csd = cdsb;
        Cdd = -(cgsb + cdsb + cbsb);
        dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

        dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
      }

      dxpart = 1.0 - sxpart;
      ddxpart_dVd = -dsxpart_dVd;
      ddxpart_dVg = -dsxpart_dVg;
      ddxpart_dVs = -dsxpart_dVs;
      ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

      qgd = qgdo;
      qgs = qgso;
      qgb = paramPtr->cgbo * vgb;
      qgate = qgd + qgs + qgb;
      qbulk = -qgb;
      qsrc = -qgs;
      qdrn = -(qgate + qbulk + qsrc);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 bsim3 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * dFdxdVp = extData.dFdxdVpVectorRawPtr;

  double coef(0.0);

  cdreq_Jdxp = 0.0;
  ceqbd_Jdxp = 0.0;
  ceqbs_Jdxp = 0.0;

  // Do a few auxilliary calculations, derived from  3f5.
  // load current vector
  if (mode >= 0)
  {
    Gm = gm;
    Gmbs = gmbs;
    FwdSum = Gm + Gmbs;
    RevSum = 0.0;

    cdreq =  model_.dtype * (cd);
    ceqbd = -model_.dtype * (csub);

    ceqbs = 0.0;

    gbbdp = -gbds;
    gbbsp = (gbds + gbgs + gbbs);

    gbdpg = gbgs;
    gbdpdp = gbds;
    gbdpb = gbbs;
    gbdpsp = -(gbdpg + gbdpdp + gbdpb);

    gbspg = 0.0;
    gbspdp = 0.0;
    gbspb = 0.0;
    gbspsp = 0.0;
  }
  else
  {
    Gm = -gm;
    Gmbs = -gmbs;
    FwdSum = 0.0;
    RevSum = -(Gm + Gmbs);

    cdreq = -model_.dtype * (cd);
    ceqbs = -model_.dtype * (csub);

    ceqbd = 0.0;

    gbbsp = -gbds;
    gbbdp = (gbds + gbgs + gbbs);

    gbdpg = 0.0;
    gbdpsp = 0.0;
    gbdpb = 0.0;
    gbdpdp = 0.0;

    gbspg = gbgs;
    gbspsp = gbds;
    gbspb = gbbs;
    gbspdp = -(gbspg + gbspsp + gbspb);
  }

  if (model_.dtype > 0)
  {
    ceqbs += (cbs);
    ceqbd += (cbd);
  }
  else
  {
    ceqbs -= (cbs);
    ceqbd -= (cbd);
  }

  if (drainConductance != 0.0)
  {
    fVec[li_Drain] += Idrain*numberParallel;
  }
  if (sourceConductance != 0.0)
  {
    fVec[li_Source] += Isource*numberParallel;
  }
  fVec[li_Bulk] += (ceqbs + ceqbd)*numberParallel;
  fVec[li_DrainPrime] += (-(ceqbd - cdreq)-Idrain)*numberParallel;
  fVec[li_SourcePrime] += (-(cdreq + ceqbs)-Isource)*numberParallel;

  // lead current support
  if( loadLeadCurrent )
  {
    double * storeLeadF = extData.nextStoVectorRawPtr;
    if (drainConductance != 0.0)
    {
      storeLeadF[li_store_dev_id] = Idrain*numberParallel;
    }
    else
    {
      storeLeadF[li_store_dev_id] = (-(ceqbd - cdreq)-Idrain)*numberParallel;
    }
    if (sourceConductance != 0.0)
    {
      storeLeadF[li_store_dev_is] = Isource*numberParallel;
    }
    else
    {
      storeLeadF[li_store_dev_is] = (-(cdreq + ceqbs)-Isource)*numberParallel;
    }
    storeLeadF[li_store_dev_ig] = 0.0;
    storeLeadF[li_store_dev_ib] = (ceqbs + ceqbd)*numberParallel;
  }

  // Initial condition support

  if( solState.dcopFlag && icVDSGiven )
  {
    coef = extData.nextSolVectorRawPtr[li_Ids];
    fVec[li_Drain] += coef;
    fVec[li_Source] += -coef;
    if( loadLeadCurrent )
    {
      double * storeLeadF = extData.nextStoVectorRawPtr;
      storeLeadF[li_store_dev_id]= coef;
      storeLeadF[li_store_dev_is]= -coef;
    }
  }

  if( solState.dcopFlag && icVGSGiven )
  {
    coef = extData.nextSolVectorRawPtr[li_Igs];
    fVec[li_Gate] += coef;
    fVec[li_Source] += -coef;
    if( loadLeadCurrent )
    {
      double * storeLeadF = extData.nextStoVectorRawPtr;
      storeLeadF[li_store_dev_ig]= coef;
      storeLeadF[li_store_dev_is]= -coef;
    }
  }


  if( solState.dcopFlag && icVBSGiven )
  {
    coef = extData.nextSolVectorRawPtr[li_Ibs];
    fVec[li_Bulk] += coef;
    fVec[li_Source] += -coef;
    if( loadLeadCurrent )
    {
      double * storeLeadF = extData.nextStoVectorRawPtr;
      storeLeadF[li_store_dev_ib]= coef;
      storeLeadF[li_store_dev_is]= -coef;
    }
  }



  //////////////////////////////////////////////////
  // limiting section:
  if (devOptions.voltageLimiterFlag)
  {
    if (!origFlag)
    {
      if (mode >= 0)
      {
        // option 1
        double tmp = model_.dtype * (-gds * (vds-vds_orig) -
                         Gm * (vgs-vgs_orig) -
                       Gmbs * (vbs-vbs_orig));

        cdreq_Jdxp += tmp;
        cdreq      += tmp;

        tmp = -model_.dtype * (-gbds * (vds-vds_orig) -
                      gbgs * (vgs-vgs_orig) -
                      gbbs * (vbs-vbs_orig));
        ceqbd_Jdxp += tmp;
        ceqbd      += tmp;
      }
      else
      {
        // option 2
        double tmp = -model_.dtype * (gds * (vds-vds_orig) +
                         Gm * (vgd-vgd_orig) +
                       Gmbs * (vbd-vbd_orig));
        cdreq_Jdxp += tmp;
        cdreq      += tmp;

        tmp = -model_.dtype * (gbds * (vds-vds_orig) -
                     gbgs * (vgd-vgd_orig) -
                     gbbs * (vbd-vbd_orig));
        ceqbd_Jdxp += tmp;
        ceqbd      += tmp;
      }


      if (model_.dtype > 0)
      {
        ceqbs_Jdxp += (-gbs*(vbs-vbs_orig));
        ceqbs      += (-gbs*(vbs-vbs_orig));

        ceqbd_Jdxp += (-gbd*(vbd-vbd_orig));
        ceqbd      += (-gbd*(vbd-vbd_orig));
      }
      else
      {
        ceqbs_Jdxp -= (-gbs*(vbs-vbs_orig) );
        ceqbs      -= (-gbs*(vbs-vbs_orig) );

        ceqbd_Jdxp -= (-gbd*(vbd-vbd_orig) );
        ceqbd      -= (-gbd*(vbd-vbd_orig) );
      }
      dFdxdVp[li_Bulk] += -(ceqbs_Jdxp+ceqbd_Jdxp)*numberParallel;
      dFdxdVp[li_DrainPrime] += (ceqbd_Jdxp-cdreq_Jdxp)*numberParallel;
      dFdxdVp[li_SourcePrime] += (cdreq_Jdxp+ceqbs_Jdxp) *numberParallel;
    } // orig flag.
  } // voltage limiter flag

  // Row associated with icVBS
  if( solState.dcopFlag && icVBSGiven )
  {
    // get the voltage drop from the previous solution
    double cVs = extData.nextSolVectorRawPtr[li_Source];
    double cVb = extData.nextSolVectorRawPtr[li_Bulk];
    fVec[li_Ibs] += (cVb - cVs - icVBS);
  }

  // Row associated with icVDS
  if( solState.dcopFlag && icVDSGiven )
  {
    // get the voltage drop from the previous solution
    double cVd = extData.nextSolVectorRawPtr[li_Drain];
    double cVs = extData.nextSolVectorRawPtr[li_Source];
    fVec[li_Ids] += (cVd - cVs - icVDS);
  }

  // Row associated with icVGS
  if( solState.dcopFlag && icVGSGiven )
  {
    // get the voltage drop from the previous solution
    double cVg = extData.nextSolVectorRawPtr[li_Gate];
    double cVs = extData.nextSolVectorRawPtr[li_Source];

    fVec[li_Igs] += (cVg - cVs - icVGS);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 bsim3 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  if (!(solState.dcopFlag) && solState.initTranFlag && solState.newtonIter==0)
  {
    // do nothing, as for this special case q is always zero.
  }
  else
  {
    N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

    // Row corresponding to the KCL for the drain node: NOTHING

    // Row corresponding to the KCL for the source node: NOTHING

    // Row corresponding to the KCL for the gate node:
    // Check this later.   ERK.  See the comments in the function
    // loadDAEdQdx, regarding ggtg, ggtb, etc.
    //
    // For now I am leaving out the gg terms, as they are zero when
    // nqsMod=0, which is always true.
    //
    dQdx[li_Gate][AGateEquGateNodeOffset]
      += (CAPcggb )*numberParallel;
    dQdx[li_Gate][AGateEquBulkNodeOffset]
      -= (CAPcggb + CAPcgdb + CAPcgsb )*numberParallel;
    dQdx[li_Gate][AGateEquDrainPrimeNodeOffset]
      += (CAPcgdb )*numberParallel;
    dQdx[li_Gate][AGateEquSourcePrimeNodeOffset]
      += (CAPcgsb )*numberParallel;

    // Row corresponding to the KCL for the bulk node:
    dQdx[li_Bulk][ABulkEquGateNodeOffset]
      += (CAPcbgb)*numberParallel;

    dQdx[li_Bulk][ABulkEquBulkNodeOffset]
      += (- CAPcbgb - CAPcbdb - CAPcbsb)*numberParallel;

    dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]
      += (CAPcbdb)*numberParallel;

    dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]
      += (CAPcbsb)*numberParallel;


    // Row corresponding to the KCL for the drain prime node:
    dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]
      -= (+ CAPcdgb + CAPcddb + CAPcdsb )*numberParallel;

    dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]
      += (CAPcdgb) *numberParallel;

    dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]
      += (+ CAPcddb )*numberParallel;

    dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]
      -= (- CAPcdsb) *numberParallel;

    // Row corresponding to the KCL for the source prime node:
    dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]
      += (CAPcsgb) *numberParallel;

    dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]
      -= (+ CAPcsgb + CAPcsdb + CAPcssb) *numberParallel;

    dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]
      -= (- CAPcsdb) *numberParallel;

    dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
      += (+ CAPcssb) *numberParallel;

    // Row associated with the charge equation
    // This is currently not supported.
    if (nqsMod)
    {
      string msg;
      msg = "Instance::loadDAEdQdx";
      msg += " nqsMod=1 is not ready yet.  Re-run with nqsMod=0\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 bsim3 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  // Row corresponding to the KCL for the drain node:
  dFdx [li_Drain][ADrainEquDrainNodeOffset]
    += drainConductance*numberParallel;
  dFdx [li_Drain][ADrainEquDrainPrimeNodeOffset]
    -= drainConductance*numberParallel;

  // Extra term for initial conditions on Vds in operating point
  if( solState.dcopFlag && icVDSGiven )
  {
    dFdx[li_Drain][ADrainEquIdsOffset] += 1.0;
  }

  // Row corresponding to the KCL for the source node:
  dFdx [li_Source][ASourceEquSourceNodeOffset]
    += sourceConductance*numberParallel;
  dFdx [li_Source][ASourceEquSourcePrimeNodeOffset]
    -= sourceConductance*numberParallel;

  // Extra term for initial conditions on Vbs in operating point
  if( solState.dcopFlag && icVBSGiven )
  {
    dFdx[li_Source][ASourceEquIbsOffset] -= 1.0;
  }
  // Extra term for initial conditions on Vds in operating point
  if( solState.dcopFlag && icVDSGiven )
  {
    dFdx[li_Source][ASourceEquIdsOffset] -= 1.0;
  }
  // Extra term for initial conditions on Vgs in operating point
  if( solState.dcopFlag && icVGSGiven )
  {
    dFdx[li_Source][ASourceEquIgsOffset] -= 1.0;
  }

  // Row corresponding to the KCL for the gate node: NOTHING
  // Check this later.   ERK.
  //
  //  The terms beginning with "gc" (gcggb, etc.) definately do NOT
  //  belong here.  I'm not sure aboug the gg terms.  On one hand, the
  //  rhs vector component for the gate node ONLY seems to take
  //  capacitive currents, which implies that all of these are capacitive
  //  conductances.  On the other hand, the gg terms do not appear to
  //  have been created by multiplying by ag0 = pdt = 1/dt.  Generally
  //  capacitive conductances are of the form g = C/dt, and the gg terms
  //  do not have this form.
  //
  //  For now, the gg issue is moot b/c those terms are only nonzero
  //  if nqsMod = 1, which is not a supported option.

  //  However, the gg
  //  terms (ggtg, ggtb, ggtd and ggts)
  //
  //(*JMatPtr)[li_Gate][AGateEquGateNodeOffset]
  //  += (gcggb - ggtg)*numberParallel;
  //(*JMatPtr)[li_Gate][AGateEquBulkNodeOffset]
  //  -= (gcggb + gcgdb + gcgsb + ggtb)*numberParallel;
  //(*JMatPtr)[li_Gate][AGateEquDrainPrimeNodeOffset]
  //  += (gcgdb - ggtd)*numberParallel;
  //(*JMatPtr)[li_Gate][AGateEquSourcePrimeNodeOffset]
  //  += (gcgsb - ggts)*numberParallel;

  // initial conditions on gate node
  // Extra term for initial conditions on Vgs in operating point
  if( solState.dcopFlag && icVGSGiven )
  {
    dFdx[li_Gate][AGateEquIgsOffset] += 1.0;
  }

  // Row corresponding to the KCL for the bulk node:
  dFdx [li_Bulk][ABulkEquGateNodeOffset]
    += (- gbgs)*numberParallel;

  dFdx [li_Bulk][ABulkEquBulkNodeOffset]
    += (gbd + gbs - gbbs)*numberParallel;

  dFdx [li_Bulk][ABulkEquDrainPrimeNodeOffset]
    += (- gbd + gbbdp)*numberParallel;

  dFdx [li_Bulk][ABulkEquSourcePrimeNodeOffset]
    += (- gbs + gbbsp)*numberParallel;

  // Extra term for initial conditions on Vbs in operating point
  if( solState.dcopFlag && icVBSGiven )
  {
    dFdx[li_Bulk][ABulkEquIbsOffset] += 1.0;
  }

  // Row corresponding to the KCL for the drain prime node:
  dFdx [li_DrainPrime][ADrainPrimeEquDrainNodeOffset]
    -= drainConductance*numberParallel;

  dFdx [li_DrainPrime][ADrainPrimeEquBulkNodeOffset]
    -= (gbd - Gmbs - dxpart*ggtb
       - T1global*ddxpart_dVb - gbdpb)*numberParallel;

  dFdx [li_DrainPrime][ADrainPrimeEquGateNodeOffset]
    += (Gm + dxpart*ggtg + T1global*ddxpart_dVg + gbdpg)
    *numberParallel;

  dFdx [li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]
    += (drainConductance + gds + gbd + RevSum + dxpart*ggtd
       + T1global*ddxpart_dVd + gbdpdp)*numberParallel;

  dFdx [li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]
    -= (gds + FwdSum - dxpart*ggts - T1global*ddxpart_dVs - gbdpsp)
    *numberParallel;

  // Row corresponding to the KCL for the source prime node:
  dFdx [li_SourcePrime][ASourcePrimeEquGateNodeOffset]
    += (- Gm + sxpart*ggtg + T1global*dsxpart_dVg + gbspg)
    *numberParallel;

  dFdx [li_SourcePrime][ASourcePrimeEquBulkNodeOffset]
    -= (gbs + Gmbs - sxpart*ggtb
       - T1global*dsxpart_dVb - gbspb)*numberParallel;

  dFdx [li_SourcePrime][ASourcePrimeEquSourceNodeOffset]
    -= sourceConductance*numberParallel;

  dFdx [li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]
    -= (gds + RevSum - sxpart*ggtd - T1global*dsxpart_dVd - gbspdp)
    *numberParallel;

  dFdx [li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
    += (sourceConductance + gds + gbs + FwdSum + sxpart*ggts
       + T1global*dsxpart_dVs + gbspsp)*numberParallel;

  // Row associated with the charge equation
  if (nqsMod)
  {
    string msg;
    msg = "Instance::loadDAEdFdx";
    msg += " nqsMod=1 is not ready yet.  Re-run with nqsMod=0\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // Initial condition rows
  // Row associated with icVBS
  if( icVBSGiven )
  {
    if( solState.dcopFlag  )
    {
      dFdx[li_Ibs][icVBSEquVbOffset] += 1.0;
      dFdx[li_Ibs][icVBSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Ibs][icVBSEquIbsOffset] += 1.0;
    }
  }

  // Row associated with icVDS
  if( icVDSGiven )
  {
    if( solState.dcopFlag  )
    {
      dFdx[li_Ids][icVDSEquVdOffset] += 1.0;
      dFdx[li_Ids][icVDSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Ids][icVDSEquIdsOffset] += 1.0;
    }
  }

  // Row associated with icVGS
  if( icVGSGiven )
  {
    if( solState.dcopFlag  )
    {
      dFdx[li_Igs][icVGSEquVgOffset] += 1.0;
      dFdx[li_Igs][icVGSEquVsOffset] -= 1.0;
    }
    else
    {
      dFdx[li_Igs][icVGSEquIgsOffset] += 1.0;
    }
  }

  return true;
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

  if( icVBSGiven )
  {
    extData.currStoVectorRawPtr[li_store_vbs] = icVBS;
    extData.nextStoVectorRawPtr[li_store_vbs] = icVBS;
  }

  if( icVDSGiven )
  {
    extData.currStoVectorRawPtr[li_store_vds] = icVDS;
    extData.nextStoVectorRawPtr[li_store_vds] = icVDS;
  }

  if( icVGSGiven )
  {
    extData.currStoVectorRawPtr[li_store_vgs] = icVGS;
    extData.nextStoVectorRawPtr[li_store_vgs] = icVGS;
  }

  return bsuccess;
}


// Additional Declarations

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{
  string msg;
  cox = 3.453133e-11 / tox;
  if (!given("TOXM"))    toxm = tox;
  if (!given("DSUB"))  dsub = drout;
  if (!given("LLC"))  Llc = Ll;
  if (!given("LWC"))  Lwc = Lw;
  if (!given("LWLC")) Lwlc = Lwl;
  if (!given("WLC"))  Wlc = Wl;
  if (!given("WWL")) Wwlc = Wwl;
  if (!given("WWLC")) Wwlc = Wwl;
  if (!given("DWC"))  dwc = Wint;
  if (!given("DLC"))  dlc = Lint;

  if (!given("CF"))
  {
    double C1 = 2.0 * CONSTEPSOX;
    double C5 = M_PI;
    double C2 = 1.0 + (0.4e-6 / tox);
    double C3 = log(C2);
    cf = C1*C3/C5;
  }

  if (!given("CGDO"))
  {
    if (given("DLC") && (dlc > 0.0)) cgdo = dlc * cox - cgdl ;
    else                         cgdo = 0.6 * xj * cox;
  }

  if (!given("CGSO"))
  {
    if (given("DLC") && (dlc > 0.0)) cgso = dlc * cox - cgsl ;
    else                         cgso = 0.6 * xj * cox;
  }

  if (!given("CGBO")) cgbo = 2.0 * dwc * cox;

  if (!given("CJSWG"))
    unitLengthGateSidewallJctCap = unitLengthSidewallJctCap ;

  if (!given("PBSWG"))
    GatesidewallJctPotential = sidewallJctPotential;

  if (!given("MJSWG"))
    bulkJctGateSideGradingCoeff = bulkJctSideGradingCoeff;


  // More initializations:  taken from b3temp.c:
  if (bulkJctPotential < 0.1)
  {
    bulkJctPotential = 0.1;
    msg = "Given pb is less than 0.1. Pb is set to 0.1.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
  }

  if (sidewallJctPotential < 0.1)
  {
    sidewallJctPotential = 0.1;
    msg = "Given pbsw is less than 0.1. Pbsw is set to 0.1.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
  }

  if (GatesidewallJctPotential < 0.1)
  {
    GatesidewallJctPotential = 0.1;
    msg = "Given pbswg is less than 0.1. Pbswg is set to 0.1.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
  }

  vcrit   = CONSTvt0 * log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
  factor1 = sqrt(CONSTEPSSI / CONSTEPSOX * tox);

  Vtm0 = CONSTKoverQ * tnom;
  Eg0  = CONSTEg0 - CONSTalphaEg * tnom * tnom / (tnom + CONSTbetaEg);
  ni   = CONSTNi0 * (tnom / CONSTREFTEMP) * sqrt(tnom / CONSTREFTEMP)
    * exp(21.5565981 - Eg0 / (2.0 * Vtm0));

  // If there are any time dependent parameters, set their values at for
  // the current time.

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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/00
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                                  SolverState & ss1,
                                                  DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
    modType       (0),
    dtype         (CONSTNMOS),
    mobMod        (0),
    capMod        (0),
    noiMod        (0),
    binUnit       (0),
    paramChk      (0),
    version       ("3.2.2"),
    tox           (0),
    toxm          (0),
    cdsc          (0),
    cdscb         (0),
    cdscd         (0),
    cit           (0),
    nfactor       (0),
    xj            (0),
    vsat          (0),
    at            (0),
    a0            (0),
    ags           (0),
    a1            (0),
    a2            (0),
    keta          (0),
    nsub          (0),
    npeak         (0),
    ngate         (0),
    gamma1        (0),
    gamma2        (0),
    vbx           (0),
    vbm           (0),
    xt            (0),
    k1            (0),
    kt1           (0),
    kt1l          (0),
    kt2           (0),
    k2            (0),
    k3            (0),
    k3b           (0),
    w0            (0),
    nlx           (0),
    dvt0          (0),
    dvt1          (0),
    dvt2          (0),
    dvt0w         (0),
    dvt1w         (0),
    dvt2w         (0),
    drout         (0),
    dsub          (0),
    vth0          (0),
    ua            (0),
    ua1           (0),
    ub            (0),
    ub1           (0),
    uc            (0),
    uc1           (0),
    u0            (0),
    ute           (0),
    voff          (0),
    delta         (0),
    rdsw          (0),
    prwg          (0),
    prwb          (0),
    prt           (0),
    eta0          (0),
    etab          (0),
    pclm          (0),
    pdibl1        (0),
    pdibl2        (0),
    pdiblb        (0),
    pscbe1        (0),
    pscbe2        (0),
    pvag          (0),
    wr            (0),
    dwg           (0),
    dwb           (0),
    b0            (0),
    b1            (0),
    alpha0        (0),
    alpha1        (0),
    beta0         (0),
    ijth          (0),
    vfb           (0),
    elm           (0),
    cgsl          (0),
    cgdl          (0),
    ckappa        (0),
    cf            (0),
    vfbcv         (0),
    clc           (0),
    cle           (0),
    dwc           (0),
    dlc           (0),
    noff          (0),
    voffcv        (0),
    acde          (0),
    moin          (0),
    tcj           (0),
    tcjsw         (0),
    tcjswg        (0),
    tpb           (0),
    tpbsw         (0),
    tpbswg        (0),
    lcdsc         (0),
    lcdscb        (0),
    lcdscd        (0),
    lcit          (0),
    lnfactor      (0),
    lxj           (0),
    lvsat         (0),
    lat           (0),
    la0           (0),
    lags          (0),
    la1           (0),
    la2           (0),
    lketa         (0),
    lnsub         (0),
    lnpeak        (0),
    lngate        (0),
    lgamma1       (0),
    lgamma2       (0),
    lvbx          (0),
    lvbm          (0),
    lxt           (0),
    lk1           (0),
    lkt1          (0),
    lkt1l         (0),
    lkt2          (0),
    lk2           (0),
    lk3           (0),
    lk3b          (0),
    lw0           (0),
    lnlx          (0),
    ldvt0         (0),
    ldvt1         (0),
    ldvt2         (0),
    ldvt0w        (0),
    ldvt1w        (0),
    ldvt2w        (0),
    ldrout        (0),
    ldsub         (0),
    lvth0         (0),
    lua           (0),
    lua1          (0),
    lub           (0),
    lub1          (0),
    luc           (0),
    luc1          (0),
    lu0           (0),
    lute          (0),
    lvoff         (0),
    ldelta        (0),
    lrdsw         (0),
    lprwg         (0),
    lprwb         (0),
    lprt          (0),
    leta0         (0),
    letab         (0),
    lpclm         (0),
    lpdibl1       (0),
    lpdibl2       (0),
    lpdiblb       (0),
    lpscbe1       (0),
    lpscbe2       (0),
    lpvag         (0),
    lwr           (0),
    ldwg          (0),
    ldwb          (0),
    lb0           (0),
    lb1           (0),
    lalpha0       (0),
    lalpha1       (0),
    lbeta0        (0),
    lvfb          (0),
    lelm          (0),
    lcgsl         (0),
    lcgdl         (0),
    lckappa       (0),
    lcf           (0),
    lclc          (0),
    lcle          (0),
    lvfbcv        (0),
    lnoff         (0),
    lvoffcv       (0),
    lacde         (0),
    lmoin         (0),
    wcdsc         (0),
    wcdscb        (0),
    wcdscd        (0),
    wcit          (0),
    wnfactor      (0),
    wxj           (0),
    wvsat         (0),
    wat           (0),
    wa0           (0),
    wags          (0),
    wa1           (0),
    wa2           (0),
    wketa         (0),
    wnsub         (0),
    wnpeak        (0),
    wngate        (0),
    wgamma1       (0),
    wgamma2       (0),
    wvbx          (0),
    wvbm          (0),
    wxt           (0),
    wk1           (0),
    wkt1          (0),
    wkt1l         (0),
    wkt2          (0),
    wk2           (0),
    wk3           (0),
    wk3b          (0),
    ww0           (0),
    wnlx          (0),
    wdvt0         (0),
    wdvt1         (0),
    wdvt2         (0),
    wdvt0w        (0),
    wdvt1w        (0),
    wdvt2w        (0),
    wdrout        (0),
    wdsub         (0),
    wvth0         (0),
    wua           (0),
    wua1          (0),
    wub           (0),
    wub1          (0),
    wuc           (0),
    wuc1          (0),
    wu0           (0),
    wute          (0),
    wvoff         (0),
    wdelta        (0),
    wrdsw         (0),
    wprwg         (0),
    wprwb         (0),
    wprt          (0),
    weta0         (0),
    wetab         (0),
    wpclm         (0),
    wpdibl1       (0),
    wpdibl2       (0),
    wpdiblb       (0),
    wpscbe1       (0),
    wpscbe2       (0),
    wpvag         (0),
    wwr           (0),
    wdwg          (0),
    wdwb          (0),
    wb0           (0),
    wb1           (0),
    walpha0       (0),
    walpha1       (0),
    wbeta0        (0),
    wvfb          (0),
    welm          (0),
    wcgsl         (0),
    wcgdl         (0),
    wckappa       (0),
    wcf           (0),
    wclc          (0),
    wcle          (0),
    wvfbcv        (0),
    wnoff         (0),
    wvoffcv       (0),
    wacde         (0),
    wmoin         (0),
    pcdsc         (0),
    pcdscb        (0),
    pcdscd        (0),
    pcit          (0),
    pnfactor      (0),
    pxj           (0),
    pvsat         (0),
    pat           (0),
    pa0           (0),
    pags          (0),
    pa1           (0),
    pa2           (0),
    pketa         (0),
    pnsub         (0),
    pnpeak        (0),
    pngate        (0),
    pgamma1       (0),
    pgamma2       (0),
    pvbx          (0),
    pvbm          (0),
    pxt           (0),
    pk1           (0),
    pkt1          (0),
    pkt1l         (0),
    pkt2          (0),
    pk2           (0),
    pk3           (0),
    pk3b          (0),
    pw0           (0),
    pnlx          (0),
    pdvt0         (0),
    pdvt1         (0),
    pdvt2         (0),
    pdvt0w        (0),
    pdvt1w        (0),
    pdvt2w        (0),
    pdrout        (0),
    pdsub         (0),
    pvth0         (0),
    pua           (0),
    pua1          (0),
    pub           (0),
    pub1          (0),
    puc           (0),
    puc1          (0),
    pu0           (0),
    pute          (0),
    pvoff         (0),
    pdelta        (0),
    prdsw         (0),
    pprwg         (0),
    pprwb         (0),
    pprt          (0),
    peta0         (0),
    petab         (0),
    ppclm         (0),
    ppdibl1       (0),
    ppdibl2       (0),
    ppdiblb       (0),
    ppscbe1       (0),
    ppscbe2       (0),
    ppvag         (0),
    pwr           (0),
    pdwg          (0),
    pdwb          (0),
    pb0           (0),
    pb1           (0),
    palpha0       (0),
    palpha1       (0),
    pbeta0        (0),
    pvfb          (0),
    pelm          (0),
    pcgsl         (0),
    pcgdl         (0),
    pckappa       (0),
    pcf           (0),
    pclc          (0),
    pcle          (0),
    pvfbcv        (0),
    pnoff         (0),
    pvoffcv       (0),
    pacde         (0),
    pmoin         (0),
    tnom          (devOptions.tnom),
    cgso          (0),
    cgdo          (0),
    cgbo          (0),
    xpart         (0),
    cFringOut     (0),
    cFringMax     (0),
    sheetResistance              (0),
    jctSatCurDensity             (0),
    jctSidewallSatCurDensity     (0),
    bulkJctPotential             (0),
    bulkJctBotGradingCoeff       (0),
    bulkJctSideGradingCoeff      (0),
    bulkJctGateSideGradingCoeff  (0),
    sidewallJctPotential         (0),
    GatesidewallJctPotential     (0),
    unitAreaJctCap               (0),
    unitLengthSidewallJctCap     (0),
    unitLengthGateSidewallJctCap (0),
    jctEmissionCoeff             (0),
    jctTempExponent              (0),
    Lint     (0),
    Ll       (0),
    Llc      (0),
    Lln      (0),
    Lw       (0),
    Lwc      (0),
    Lwn      (0),
    Lwl      (0),
    Lwlc     (0),
    model_l  (0),
    model_w  (0),
    Lmin     (0),
    Lmax     (0),
    Wmin     (0),
    Wmax     (0),
    Wint     (0),
    Wl       (0),
    Wlc      (0),
    Wln      (0),
    Ww       (0),
    Wwc      (0),
    Wwn      (0),
    Wwl      (0),
    Wwlc     (0),
    vtm      (0),
    cox      (0),
    cof1     (0),
    cof2     (0),
    cof3     (0),
    cof4     (0),
    vcrit    (0),
    factor1  (0),
    PhiB     (0),
    PhiBSW   (0),
    PhiBSWG  (0),
    oxideTrapDensityA (0),
    oxideTrapDensityB (0),
    oxideTrapDensityC (0),
    em       (0),
    ef       (0),
    af       (0),
    kf       (0),
    npeakGiven     (0),
    gamma1Given    (0),
    gamma2Given    (0),
    k1Given        (0),
    k2Given        (0),
    nsubGiven      (0),
    xtGiven        (0),
    vbxGiven       (0),
    vbmGiven       (0),
    vfbGiven       (0),
    vth0Given      (0),
    Vtm0            (0.0),
    Eg0             (0.0),
    ni              (0.0)
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
      std::ostringstream oss;

      oss << "Error in " << netlistLocation() << "\n"
        <<  "Could not recognize the type for model " <<  getName();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
#ifdef Xyce_BSIM3_USE_DEFL
  if (!given("L"))
    model_l=devOptions.defl;
  if (!given("W"))
    model_w=devOptions.defw;
#endif
  if (!given("TNOM"))
    tnom = devOptions.tnom;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  if (!given("VTH0"))
    vth0 = (dtype == CONSTNMOS) ? 0.7 : -0.7;
  if (!given("UC"))
    uc  = (mobMod == 3) ? -0.0465 : -0.0465e-9;
  if (!given("UC1"))
    uc1 = (mobMod == 3) ? -0.056 : -0.056e-9;
  if (!given("U0"))
    u0  = (dtype == CONSTNMOS) ? 0.067 : 0.025;
  if (!given("NOIA"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityA = 1e20;
    else
      oxideTrapDensityA = 9.9e18;
  }
  if (!given("NOIB"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityB = 5e4;
    else
      oxideTrapDensityB = 2.4e3;
  }
  if (!given("NOIC"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityC = -1.4e-12;
    else
      oxideTrapDensityC = 1.4e-12;
  }

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/00
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances (std::ostream &os) const
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
// Creation Date  : 10/26/2004
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
// MOSFET_B3 Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/02/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & mi = *(*it);

    bool btmp = mi.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    // voltage drops:
    double * stoVec = mi.extData.nextStoVectorRawPtr;
    stoVec[mi.li_store_vbs] = mi.vbs;
    stoVec[mi.li_store_vgs] = mi.vgs;
    stoVec[mi.li_store_vds] = mi.vds;
    stoVec[mi.li_store_vbd] = mi.vbd;
    stoVec[mi.li_store_von] = mi.von;

    // intrinsic capacitors:
    staVec[mi.li_state_qb] = mi.qb;
    staVec[mi.li_state_qg] = mi.qg;
    staVec[mi.li_state_qd] = mi.qd;

    // parasitic capacitors:
    staVec[mi.li_state_qbs] = mi.qbs;
    staVec[mi.li_state_qbd] = mi.qbd;

    if( mi.nqsMod )
    {
      staVec[mi.li_state_qcheq] = mi.qcheq;
      staVec[mi.li_state_qcdump] = mi.qcdump;
    }

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
      // re-set the state vector pointer that we are using to the "current"
      // pointer, rather than the "next" pointer.
      double * currStaVec = mi.extData.currStaVectorRawPtr;

      // intrinsic capacitors:
      currStaVec[mi.li_state_qb] = mi.qb;
      currStaVec[mi.li_state_qg] = mi.qg;
      currStaVec[mi.li_state_qd] = mi.qd;

      // parasitic capacitors:
      currStaVec[mi.li_state_qbs] = mi.qbs;
      currStaVec[mi.li_state_qbd] = mi.qbd;

      if( mi.nqsMod )
      {
        currStaVec[mi.li_state_qcheq] = mi.qcheq;
        currStaVec[mi.li_state_qcdump] = mi.qcdump;
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
// Creation Date : 12/02/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & mi = *(*it);

    double * dFdxdVp = mi.extData.dFdxdVpVectorRawPtr;
    double * dQdxdVp = mi.extData.dQdxdVpVectorRawPtr;

    double coef(0.0);
    // F-vector:
    mi.cdreq_Jdxp = 0.0;
    mi.ceqbd_Jdxp = 0.0;
    mi.ceqbs_Jdxp = 0.0;

    // Do a few auxilliary calculations, derived from  3f5.
    // load current vector
    if (mi.mode >= 0)
    {
      mi.Gm = mi.gm;
      mi.Gmbs = mi.gmbs;
      mi.FwdSum = mi.Gm + mi.Gmbs;
      mi.RevSum = 0.0;

      mi.cdreq =  mi.model_.dtype * (mi.cd);
      mi.ceqbd = -mi.model_.dtype * (mi.csub);

      mi.ceqbs = 0.0;

      mi.gbbdp = -mi.gbds;
      mi.gbbsp = (mi.gbds + mi.gbgs + mi.gbbs);

      mi.gbdpg = mi.gbgs;
      mi.gbdpdp = mi.gbds;
      mi.gbdpb = mi.gbbs;
      mi.gbdpsp = -(mi.gbdpg + mi.gbdpdp + mi.gbdpb);

      mi.gbspg = 0.0;
      mi.gbspdp = 0.0;
      mi.gbspb = 0.0;
      mi.gbspsp = 0.0;
    }
    else
    {
      mi.Gm = -mi.gm;
      mi.Gmbs = -mi.gmbs;
      mi.FwdSum = 0.0;
      mi.RevSum = -(mi.Gm + mi.Gmbs);

      mi.cdreq = -mi.model_.dtype * (mi.cd);
      mi.ceqbs = -mi.model_.dtype * (mi.csub);

      mi.ceqbd = 0.0;

      mi.gbbsp = -mi.gbds;
      mi.gbbdp = (mi.gbds + mi.gbgs + mi.gbbs);

      mi.gbdpg = 0.0;
      mi.gbdpsp = 0.0;
      mi.gbdpb = 0.0;
      mi.gbdpdp = 0.0;

      mi.gbspg = mi.gbgs;
      mi.gbspsp = mi.gbds;
      mi.gbspb = mi.gbbs;
      mi.gbspdp = -(mi.gbspg + mi.gbspsp + mi.gbspb);
    }

    if (mi.model_.dtype > 0)
    {
      mi.ceqbs += (mi.cbs);
      mi.ceqbd += (mi.cbd);
    }
    else
    {
      mi.ceqbs -= (mi.cbs);
      mi.ceqbd -= (mi.cbd);
    }

    if (mi.drainConductance != 0.0)
    {
      fVec[mi.li_Drain] += mi.Idrain*mi.numberParallel;
    }
    if (mi.sourceConductance != 0.0)
    {
      fVec[mi.li_Source] += mi.Isource*mi.numberParallel;
    }

    fVec[mi.li_Bulk] += (mi.ceqbs + mi.ceqbd)*mi.numberParallel;
    fVec[mi.li_DrainPrime] += (-(mi.ceqbd - mi.cdreq)-mi.Idrain)*mi.numberParallel;
    fVec[mi.li_SourcePrime] += (-(mi.cdreq + mi.ceqbs)-mi.Isource)*mi.numberParallel;

    if( mi.loadLeadCurrent )
    {
      if (mi.drainConductance != 0.0)
      {
        storeLeadF[mi.li_store_dev_id] = mi.Idrain*mi.numberParallel;
      }
      else
      {
        storeLeadF[mi.li_store_dev_id] = (-(mi.ceqbd - mi.cdreq)-mi.Idrain)*mi.numberParallel;
      }
      if (mi.sourceConductance != 0.0)
      {
        storeLeadF[mi.li_store_dev_is] = mi.Isource*mi.numberParallel;
      }
      else
      {
        storeLeadF[mi.li_store_dev_is] = (-(mi.cdreq + mi.ceqbs)-mi.Isource)*mi.numberParallel;
      }
      storeLeadF[mi.li_store_dev_ig] = 0.0;
      storeLeadF[mi.li_store_dev_ib] = (mi.ceqbs + mi.ceqbd)*mi.numberParallel;
    }

    // Initial condition support

    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Ids];
      fVec[mi.li_Drain] += coef;
      fVec[mi.li_Source] += -coef;
      if( mi.loadLeadCurrent )
      {
        storeLeadF[mi.li_store_dev_id]= coef;
        storeLeadF[mi.li_store_dev_is]= -coef;
      }
    }

    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Igs];
      fVec[mi.li_Gate] += coef;
      fVec[mi.li_Source] += -coef;
      if( mi.loadLeadCurrent )
      {
        storeLeadF[mi.li_store_dev_ig]= coef;
        storeLeadF[mi.li_store_dev_is]= -coef;
      }
    }


    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      coef = mi.extData.nextSolVectorRawPtr[mi.li_Ibs];
      fVec[mi.li_Bulk] += coef;
      fVec[mi.li_Source] += -coef;
      if( mi.loadLeadCurrent )
      {
        storeLeadF[mi.li_store_dev_ib]= coef;
        storeLeadF[mi.li_store_dev_is]= -coef;
      }
    }

    //////////////////////////////////////////////////
    // limiting section:
    if (getDeviceOptions().voltageLimiterFlag)
    {
      if (!mi.origFlag)
      {
        if (mi.mode >= 0)
        {
          // option 1
          double tmp = mi.model_.dtype * (-mi.gds * (mi.vds-mi.vds_orig) -
                          mi.Gm * (mi.vgs-mi.vgs_orig) -
                        mi.Gmbs * (mi.vbs-mi.vbs_orig));

          mi.cdreq_Jdxp += tmp;
          mi.cdreq      += tmp;

          tmp = -mi.model_.dtype * (-mi.gbds * (mi.vds-mi.vds_orig) -
                        mi.gbgs * (mi.vgs-mi.vgs_orig) -
                        mi.gbbs * (mi.vbs-mi.vbs_orig));
          mi.ceqbd_Jdxp += tmp;
          mi.ceqbd      += tmp;
        }
        else
        {
          // option 2
          double tmp = -mi.model_.dtype * (mi.gds * (mi.vds-mi.vds_orig) +
                          mi.Gm * (mi.vgd-mi.vgd_orig) +
                        mi.Gmbs * (mi.vbd-mi.vbd_orig));
          mi.cdreq_Jdxp += tmp;
          mi.cdreq      += tmp;

          tmp = -mi.model_.dtype * (mi.gbds * (mi.vds-mi.vds_orig) -
                      mi.gbgs * (mi.vgd-mi.vgd_orig) -
                      mi.gbbs * (mi.vbd-mi.vbd_orig));
          mi.ceqbd_Jdxp += tmp;
          mi.ceqbd      += tmp;
        }


        if (mi.model_.dtype > 0)
        {
          mi.ceqbs_Jdxp += (-mi.gbs*(mi.vbs-mi.vbs_orig));
          mi.ceqbs      += (-mi.gbs*(mi.vbs-mi.vbs_orig));

          mi.ceqbd_Jdxp += (-mi.gbd*(mi.vbd-mi.vbd_orig));
          mi.ceqbd      += (-mi.gbd*(mi.vbd-mi.vbd_orig));
        }
        else
        {
          mi.ceqbs_Jdxp -= (-mi.gbs*(mi.vbs-mi.vbs_orig) );
          mi.ceqbs      -= (-mi.gbs*(mi.vbs-mi.vbs_orig) );

          mi.ceqbd_Jdxp -= (-mi.gbd*(mi.vbd-mi.vbd_orig) );
          mi.ceqbd      -= (-mi.gbd*(mi.vbd-mi.vbd_orig) );
        }

        dFdxdVp[mi.li_Bulk] += -(mi.ceqbs_Jdxp+mi.ceqbd_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_DrainPrime] += (mi.ceqbd_Jdxp-mi.cdreq_Jdxp)*mi.numberParallel;
        dFdxdVp[mi.li_SourcePrime] += (mi.cdreq_Jdxp+mi.ceqbs_Jdxp)*mi.numberParallel;

      } // orig flag.
    } // voltage limiter flag

    // Row associated with icVBS
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      // get the voltage drop from the previous solution
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];
      double cVb = mi.extData.nextSolVectorRawPtr[mi.li_Bulk];

      fVec[mi.li_Ibs] += (cVb - cVs - mi.icVBS);
    }

    // Row associated with icVDS
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      // get the voltage drop from the previous solution
      double cVd = mi.extData.nextSolVectorRawPtr[mi.li_Drain];
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];

      fVec[mi.li_Ids] += (cVd - cVs - mi.icVDS);
    }

    // Row associated with icVGS
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      // get the voltage drop from the previous solution
      double cVg = mi.extData.nextSolVectorRawPtr[mi.li_Gate];
      double cVs = mi.extData.nextSolVectorRawPtr[mi.li_Source];

      fVec[mi.li_Igs] += (cVg - cVs - mi.icVGS);
    }

    // Q-vector:

    mi.auxChargeCalculations ();

    double Qeqqg = 0.0;   // gate charge
    double Qeqqb = 0.0;   // bulk charge
    double Qeqqd = 0.0;   // drain charge
    double Qqdef = 0.0;   // nqs-related charge.
    double Qqcheq = 0.0;  // nqs-related charge.

    // These 3 vars are class variables, and are set up elsewhere.
    //double Qeqqg_Jdxp = 0.0; // limiter, related to gate  cap.
    //double Qeqqb_Jdxp = 0.0; // limiter, related to bulk  cap.
    //double Qeqqd_Jdxp = 0.0; // limiter, related to drain cap.

    if (mi.model_.dtype > 0)
    {
      Qeqqg  = mi.qg;
      Qeqqb  = mi.qb;
      Qeqqd  = mi.qd;
      Qqdef  = mi.qcdump; // this needs to be fixed...
      Qqcheq = mi.qcheq;
    }
    else  // need to convert these to charges.
    {
      Qeqqg  = -mi.qg;
      Qeqqb  = -mi.qb;
      Qeqqd  = -mi.qd;
      Qqdef  = -mi.qcdump;
      Qqcheq = -mi.qcheq;
    }

    qVec[mi.li_Gate] += Qeqqg*mi.numberParallel;
    qVec[mi.li_Bulk] += (Qeqqb)*mi.numberParallel;
    qVec[mi.li_DrainPrime] += (-(-Qeqqd))*mi.numberParallel;
    qVec[mi.li_SourcePrime] += (-(+ Qeqqg + Qeqqb + Qeqqd))*mi.numberParallel;

    if (mi.nqsMod)
    {
      // 7 equ. for nqs modification. charge equation.
      qVec[mi.li_Charge] += -(Qqcheq - Qqdef)*mi.numberParallel;
    }

    if( mi.loadLeadCurrent )
    {
      if (mi.drainConductance == 0.0)
      {
        storeLeadQ[mi.li_store_dev_id] = (-(-Qeqqd))*mi.numberParallel;
      }
      if (mi.sourceConductance == 0.0)
      {
        storeLeadQ[mi.li_store_dev_is] = (-(Qeqqg + Qeqqb + Qeqqd))*mi.numberParallel;
      }
      storeLeadQ[mi.li_store_dev_ig] = Qeqqg*mi.numberParallel;
      storeLeadQ[mi.li_store_dev_ib] = (Qeqqb)*mi.numberParallel;
    }

    //////////////////////////////////////////////////
    // limiting section:
    if (getDeviceOptions().voltageLimiterFlag)
    {
      // Need the following:
      //  Qeqqg_Jdxp
      //  Qeqqb_Jdxp
      //  Qeqqd_Jdxp
      if (mi.model_.dtype > 0)
      {
#if 0
        // no-op:
        mi.Qeqqg_Jdxp = mi.Qeqqg_Jdxp;
        mi.Qeqqb_Jdxp = mi.Qeqqb_Jdxp;
        mi.Qeqqd_Jdxp = mi.Qeqqd_Jdxp;
#endif
      }
      else
      {
        mi.Qeqqg_Jdxp = -mi.Qeqqg_Jdxp;
        mi.Qeqqb_Jdxp = -mi.Qeqqb_Jdxp;
        mi.Qeqqd_Jdxp = -mi.Qeqqd_Jdxp;
      }

      if (!mi.origFlag)
      {
        dQdxdVp[mi.li_Gate] += -mi.Qeqqg_Jdxp*mi.numberParallel;
        dQdxdVp[mi.li_Bulk] += -(+mi.Qeqqb_Jdxp)*mi.numberParallel;
        dQdxdVp[mi.li_DrainPrime] += (-mi.Qeqqd_Jdxp)*mi.numberParallel;
        dQdxdVp[mi.li_SourcePrime] += (+mi.Qeqqg_Jdxp+mi.Qeqqb_Jdxp+mi.Qeqqd_Jdxp) *mi.numberParallel;
      } // orig flag.
    } // limiter flag
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
// Creation Date : 12/02/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & mi = *(*it);

    // F-matrix:
    // Row corresponding to the KCL for the drain node:

    *mi.f_DrainEquDrainNodePtr
      += mi.drainConductance*mi.numberParallel;

    *mi.f_DrainEquDrainPrimeNodePtr
      -= mi.drainConductance*mi.numberParallel;

    // Extra term for initial conditions on Vds in operating point
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      *mi.f_DrainEquIdsPtr += 1.0;
    }

    // Row corresponding to the KCL for the source node:

    *mi.f_SourceEquSourceNodePtr
      += mi.sourceConductance*mi.numberParallel;

    *mi.f_SourceEquSourcePrimeNodePtr
      -= mi.sourceConductance*mi.numberParallel;

    // Extra term for initial conditions on Vbs in operating point
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      *mi.f_SourceEquIbsPtr -= 1.0;
    }
    // Extra term for initial conditions on Vds in operating point
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      *mi.f_SourceEquIdsPtr -= 1.0;
    }
    // Extra term for initial conditions on Vgs in operating point
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      *mi.f_SourceEquIgsPtr -= 1.0;
    }

    // Row corresponding to the KCL for the gate node: NOTHING
    // Check this later.   ERK.
    //
    //  The terms beginning with "gc" (gcggb, etc.) definately do NOT
    //  belong here.  I'm not sure aboug the gg terms.  On one hand, the
    //  rhs vector component for the gate node ONLY seems to take
    //  capacitive currents, which implies that all of these are capacitive
    //  conductances.  On the other hand, the gg terms do not appear to
    //  have been created by multiplying by ag0 = pdt = 1/dt.  Generally
    //  capacitive conductances are of the form g = C/dt, and the gg terms
    //  do not have this form.
    //
    //  For now, the gg issue is moot b/c those terms are only nonzero
    //  if mi.nqsMod = 1, which is not a supported option.

    //  However, the gg
    //  terms (mi.ggtg, mi.ggtb, ggtd and ggts)
    //
    //(*JMatPtr)[GateEquGateNodePtr
    //  += (gcggb - mi.ggtg)*mi.numberParallel;
    //(*JMatPtr)[GateEquBulkNodePtr
    //  -= (gcggb + gcgdb + gcgsb + mi.ggtb)*mi.numberParallel;
    //(*JMatPtr)[GateEquDrainPrimeNodePtr
    //  += (gcgdb - ggtd)*mi.numberParallel;
    //(*JMatPtr)[GateEquSourcePrimeNodePtr
    //  += (gcgsb - ggts)*mi.numberParallel;

    // initial conditions on gate node
    // Extra term for initial conditions on Vgs in operating point
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      *mi.f_GateEquIgsPtr += 1.0;
    }

    // Row corresponding to the KCL for the bulk node:

    *mi.f_BulkEquGateNodePtr
      += (- mi.gbgs)*mi.numberParallel;


    *mi.f_BulkEquBulkNodePtr
      += (mi.gbd + mi.gbs - mi.gbbs)*mi.numberParallel;


    *mi.f_BulkEquDrainPrimeNodePtr
      += (- mi.gbd + mi.gbbdp)*mi.numberParallel;


    *mi.f_BulkEquSourcePrimeNodePtr
      += (- mi.gbs + mi.gbbsp)*mi.numberParallel;

    // Extra term for initial conditions on Vbs in operating point
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      *mi.f_BulkEquIbsPtr += 1.0;
    }

    // Row corresponding to the KCL for the drain prime node:

    *mi.f_DrainPrimeEquDrainNodePtr
      -= mi.drainConductance*mi.numberParallel;

    *mi.f_DrainPrimeEquBulkNodePtr
      -= (mi.gbd - mi.Gmbs - mi.dxpart*mi.ggtb
        - mi.T1global*mi.ddxpart_dVb - mi.gbdpb)*mi.numberParallel;

    *mi.f_DrainPrimeEquGateNodePtr
      += (mi.Gm + mi.dxpart*mi.ggtg + mi.T1global*mi.ddxpart_dVg + mi.gbdpg)
      *mi.numberParallel;

    *mi.f_DrainPrimeEquDrainPrimeNodePtr
      += (mi.drainConductance + mi.gds + mi.gbd + mi.RevSum + mi.dxpart*mi.ggtd
        + mi.T1global*mi.ddxpart_dVd + mi.gbdpdp)*mi.numberParallel;

    *mi.f_DrainPrimeEquSourcePrimeNodePtr
      -= (mi.gds + mi.FwdSum - mi.dxpart*mi.ggts - mi.T1global*mi.ddxpart_dVs - mi.gbdpsp)
      *mi.numberParallel;

    // Row corresponding to the KCL for the source prime node:
    *mi.f_SourcePrimeEquGateNodePtr
      += (- mi.Gm + mi.sxpart*mi.ggtg + mi.T1global*mi.dsxpart_dVg + mi.gbspg)
      *mi.numberParallel;

    *mi.f_SourcePrimeEquBulkNodePtr
      -= (mi.gbs + mi.Gmbs - mi.sxpart*mi.ggtb
        - mi.T1global*mi.dsxpart_dVb - mi.gbspb)*mi.numberParallel;

    *mi.f_SourcePrimeEquSourceNodePtr
      -= mi.sourceConductance*mi.numberParallel;

    *mi.f_SourcePrimeEquDrainPrimeNodePtr
      -= (mi.gds + mi.RevSum - mi.sxpart*mi.ggtd - mi.T1global*mi.dsxpart_dVd - mi.gbspdp)
      *mi.numberParallel;

    *mi.f_SourcePrimeEquSourcePrimeNodePtr
      += (mi.sourceConductance + mi.gds + mi.gbs + mi.FwdSum + mi.sxpart*mi.ggts
        + mi.T1global*mi.dsxpart_dVs + mi.gbspsp)*mi.numberParallel;

    // Row associated with the charge equation
    if (mi.nqsMod)
    {
      string msg;
      msg = "Instance::loadDAEMatrices";
      msg += " nqsMod=1 is not ready yet.  Re-run with nqsMod=0\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }

    // Initial condition rows
    // Row associated with mi.icVBS
    if( mi.icVBSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_icVBSEquVbPtr += 1.0;
        *mi.f_icVBSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVBSEquIbsPtr += 1.0;
      }
    }

    // Row associated with mi.icVDS
    if( mi.icVDSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_icVDSEquVdPtr += 1.0;
        *mi.f_icVDSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVDSEquIdsPtr += 1.0;
      }
    }

    // Row associated with mi.icVGS
    if( mi.icVGSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        *mi.f_icVGSEquVgPtr += 1.0;
        *mi.f_icVGSEquVsPtr -= 1.0;
      }
      else
      {
        *mi.f_icVGSEquIgsPtr += 1.0;
      }
    }


    // Q-matrix:
    if (!(getSolverState().dcopFlag) && getSolverState().initTranFlag && getSolverState().newtonIter==0)
    {
      // do nothing, as for this special case q is always zero.
    }
    else
    {
      // Row corresponding to the KCL for the drain node: NOTHING

      // Row corresponding to the KCL for the source node: NOTHING

      // Row corresponding to the KCL for the gate node:
      // Check this later.   ERK.  See the comments in the function
      // loadDAE*mi.q_, regarding ggtg, ggtb, etc.
      //
      // For now I am leaving out the gg terms, as they are zero when
      // nqsMod=0, which is always true.
      //
      *mi.q_GateEquGateNodePtr
        += (mi.CAPcggb )*mi.numberParallel;
      *mi.q_GateEquBulkNodePtr
        -= (mi.CAPcggb + mi.CAPcgdb + mi.CAPcgsb )*mi.numberParallel;
      *mi.q_GateEquDrainPrimeNodePtr
        += (mi.CAPcgdb )*mi.numberParallel;
      *mi.q_GateEquSourcePrimeNodePtr
        += (mi.CAPcgsb )*mi.numberParallel;

      // Row corresponding to the KCL for the bulk node:
      *mi.q_BulkEquGateNodePtr
        += (mi.CAPcbgb)*mi.numberParallel;

      *mi.q_BulkEquBulkNodePtr
        += (- mi.CAPcbgb - mi.CAPcbdb - mi.CAPcbsb)*mi.numberParallel;

      *mi.q_BulkEquDrainPrimeNodePtr
        += (mi.CAPcbdb)*mi.numberParallel;

      *mi.q_BulkEquSourcePrimeNodePtr
        += (mi.CAPcbsb)*mi.numberParallel;


      // Row corresponding to the KCL for the drain prime node:
      *mi.q_DrainPrimeEquBulkNodePtr
        -= (+ mi.CAPcdgb + mi.CAPcddb + mi.CAPcdsb )*mi.numberParallel;

      *mi.q_DrainPrimeEquGateNodePtr
        += (mi.CAPcdgb) *mi.numberParallel;

      *mi.q_DrainPrimeEquDrainPrimeNodePtr
        += (+ mi.CAPcddb )*mi.numberParallel;

      *mi.q_DrainPrimeEquSourcePrimeNodePtr
        -= (- mi.CAPcdsb) *mi.numberParallel;

      // Row corresponding to the KCL for the source prime node:
      *mi.q_SourcePrimeEquGateNodePtr
        += (mi.CAPcsgb) *mi.numberParallel;

      *mi.q_SourcePrimeEquBulkNodePtr
        -= (+ mi.CAPcsgb + mi.CAPcsdb + mi.CAPcssb) *mi.numberParallel;

      *mi.q_SourcePrimeEquDrainPrimeNodePtr
        -= (- mi.CAPcsdb) *mi.numberParallel;

      *mi.q_SourcePrimeEquSourcePrimeNodePtr
        += (+ mi.CAPcssb) *mi.numberParallel;

      // Row associated with the charge equation
      // This is currently not supported.
      if (mi.nqsMod)
      {
        string msg;
        msg = "Master::loadDAEMatrices";
        msg += " nqsMod=1 is not ready yet.  Re-run with nqsMod=0\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
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
// Creation Date : 12/02/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  int sizeInstances = instanceContainer_.size();
  for (int i=0; i<sizeInstances; ++i)
  {
    Instance & mi = *(instanceContainer_.at(i));

    int count = 0;
    // F-matrix:
    // Row corresponding to the KCL for the drain node:

    dFdx [mi.li_Drain][mi.ADrainEquDrainNodeOffset]
      += mi.drainConductance*mi.numberParallel;

    dFdx [mi.li_Drain][mi.ADrainEquDrainPrimeNodeOffset]
      -= mi.drainConductance*mi.numberParallel;

    // Extra term for initial conditions on Vds in operating point
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      dFdx[mi.li_Drain][mi.ADrainEquIdsOffset] += 1.0;
    }

    // Row corresponding to the KCL for the source node:

    dFdx [mi.li_Source][mi.ASourceEquSourceNodeOffset]
      += mi.sourceConductance*mi.numberParallel;

    dFdx [mi.li_Source][mi.ASourceEquSourcePrimeNodeOffset]
      -= mi.sourceConductance*mi.numberParallel;

    // Extra term for initial conditions on Vbs in operating point
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      dFdx[mi.li_Source][mi.ASourceEquIbsOffset] -= 1.0;
    }
    // Extra term for initial conditions on Vds in operating point
    if( getSolverState().dcopFlag && mi.icVDSGiven )
    {
      dFdx[mi.li_Source][mi.ASourceEquIdsOffset] -= 1.0;
    }
    // Extra term for initial conditions on Vgs in operating point
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      dFdx[mi.li_Source][mi.ASourceEquIgsOffset] -= 1.0;
    }

    // Row corresponding to the KCL for the gate node: NOTHING
    // Check this later.   ERK.
    //
    //  The terms beginning with "gc" (gcggb, etc.) definately do NOT
    //  belong here.  I'm not sure aboug the gg terms.  On one hand, the
    //  rhs vector component for the gate node ONLY seems to take
    //  capacitive currents, which implies that all of these are capacitive
    //  conductances.  On the other hand, the gg terms do not appear to
    //  have been created by multiplying by ag0 = pdt = 1/dt.  Generally
    //  capacitive conductances are of the form g = C/dt, and the gg terms
    //  do not have this form.
    //
    //  For now, the gg issue is moot b/c those terms are only nonzero
    //  if mi.nqsMod = 1, which is not a supported option.

    //  However, the gg
    //  terms (mi.ggtg, mi.ggtb, ggtd and ggts)
    //
    //(*JMatPtr)[mi.li_Gate][mi.AGateEquGateNodeOffset]
    //  += (gcggb - mi.ggtg)*mi.numberParallel;
    //(*JMatPtr)[mi.li_Gate][mi.AGateEquBulkNodeOffset]
    //  -= (gcggb + gcgdb + gcgsb + mi.ggtb)*mi.numberParallel;
    //(*JMatPtr)[mi.li_Gate][mi.AGateEquDrainPrimeNodeOffset]
    //  += (gcgdb - ggtd)*mi.numberParallel;
    //(*JMatPtr)[mi.li_Gate][mi.AGateEquSourcePrimeNodeOffset]
    //  += (gcgsb - ggts)*mi.numberParallel;

    // initial conditions on gate node
    // Extra term for initial conditions on Vgs in operating point
    if( getSolverState().dcopFlag && mi.icVGSGiven )
    {
      dFdx[mi.li_Gate][mi.AGateEquIgsOffset] += 1.0;
    }

    // Row corresponding to the KCL for the bulk node:

    dFdx [mi.li_Bulk][mi.ABulkEquGateNodeOffset]
      += (- mi.gbgs)*mi.numberParallel;


    dFdx [mi.li_Bulk][mi.ABulkEquBulkNodeOffset]
      += (mi.gbd + mi.gbs - mi.gbbs)*mi.numberParallel;


    dFdx [mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset]
      += (- mi.gbd + mi.gbbdp)*mi.numberParallel;


    dFdx [mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset]
      += (- mi.gbs + mi.gbbsp)*mi.numberParallel;

    // Extra term for initial conditions on Vbs in operating point
    if( getSolverState().dcopFlag && mi.icVBSGiven )
    {
      dFdx[mi.li_Bulk][mi.ABulkEquIbsOffset] += 1.0;
    }

    // Row corresponding to the KCL for the drain prime node:

    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquDrainNodeOffset]
      -= mi.drainConductance*mi.numberParallel;


    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset]
      -= (mi.gbd - mi.Gmbs - mi.dxpart*mi.ggtb
        - mi.T1global*mi.ddxpart_dVb - mi.gbdpb)*mi.numberParallel;


    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquGateNodeOffset]
      += (mi.Gm + mi.dxpart*mi.ggtg + mi.T1global*mi.ddxpart_dVg + mi.gbdpg)
      *mi.numberParallel;


    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset]
      += (mi.drainConductance + mi.gds + mi.gbd + mi.RevSum + mi.dxpart*mi.ggtd
        + mi.T1global*mi.ddxpart_dVd + mi.gbdpdp)*mi.numberParallel;


    dFdx [mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset]
      -= (mi.gds + mi.FwdSum - mi.dxpart*mi.ggts - mi.T1global*mi.ddxpart_dVs - mi.gbdpsp)
      *mi.numberParallel;

    // Row corresponding to the KCL for the source prime node:

    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquGateNodeOffset]
      += (- mi.Gm + mi.sxpart*mi.ggtg + mi.T1global*mi.dsxpart_dVg + mi.gbspg)
      *mi.numberParallel;


    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset]
      -= (mi.gbs + mi.Gmbs - mi.sxpart*mi.ggtb
        - mi.T1global*mi.dsxpart_dVb - mi.gbspb)*mi.numberParallel;


    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquSourceNodeOffset]
      -= mi.sourceConductance*mi.numberParallel;


    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset]
      -= (mi.gds + mi.RevSum - mi.sxpart*mi.ggtd - mi.T1global*mi.dsxpart_dVd - mi.gbspdp)
      *mi.numberParallel;


    dFdx [mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset]
      += (mi.sourceConductance + mi.gds + mi.gbs + mi.FwdSum + mi.sxpart*mi.ggts
        + mi.T1global*mi.dsxpart_dVs + mi.gbspsp)*mi.numberParallel;

    // Row associated with the charge equation
    if (mi.nqsMod)
    {
      string msg;
      msg = "Master::loadDAEMatrices";
      msg += " nqsMod=1 is not ready yet.  Re-run with nqsMod=0\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }

    // Initial condition rows
    // Row associated with mi.icVBS
    if( mi.icVBSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Ibs][mi.icVBSEquVbOffset] += 1.0;
        dFdx[mi.li_Ibs][mi.icVBSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Ibs][mi.icVBSEquIbsOffset] += 1.0;
      }
    }

    // Row associated with mi.icVDS
    if( mi.icVDSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Ids][mi.icVDSEquVdOffset] += 1.0;
        dFdx[mi.li_Ids][mi.icVDSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Ids][mi.icVDSEquIdsOffset] += 1.0;
      }
    }

    // Row associated with mi.icVGS
    if( mi.icVGSGiven )
    {
      if( getSolverState().dcopFlag  )
      {
        dFdx[mi.li_Igs][mi.icVGSEquVgOffset] += 1.0;
        dFdx[mi.li_Igs][mi.icVGSEquVsOffset] -= 1.0;
      }
      else
      {
        dFdx[mi.li_Igs][mi.icVGSEquIgsOffset] += 1.0;
      }
    }


    // Q-matrix:
    if (!(getSolverState().dcopFlag) && getSolverState().initTranFlag && getSolverState().newtonIter==0)
    {
      // do nothing, as for this special case q is always zero.
    }
    else
    {
      // Row corresponding to the KCL for the drain node: NOTHING

      // Row corresponding to the KCL for the source node: NOTHING

      // Row corresponding to the KCL for the gate node:
      // Check this later.   ERK.  See the comments in the function
      // loadDAEdQdx, regarding ggtg, ggtb, etc.
      //
      // For now I am leaving out the gg terms, as they are zero when
      // nqsMod=0, which is always true.
      //
      dQdx[mi.li_Gate][mi.AGateEquGateNodeOffset]
        += (mi.CAPcggb )*mi.numberParallel;
      dQdx[mi.li_Gate][mi.AGateEquBulkNodeOffset]
        -= (mi.CAPcggb + mi.CAPcgdb + mi.CAPcgsb )*mi.numberParallel;
      dQdx[mi.li_Gate][mi.AGateEquDrainPrimeNodeOffset]
        += (mi.CAPcgdb )*mi.numberParallel;
      dQdx[mi.li_Gate][mi.AGateEquSourcePrimeNodeOffset]
        += (mi.CAPcgsb )*mi.numberParallel;

      // Row corresponding to the KCL for the bulk node:
      dQdx[mi.li_Bulk][mi.ABulkEquGateNodeOffset]
        += (mi.CAPcbgb)*mi.numberParallel;

      dQdx[mi.li_Bulk][mi.ABulkEquBulkNodeOffset]
        += (- mi.CAPcbgb - mi.CAPcbdb - mi.CAPcbsb)*mi.numberParallel;

      dQdx[mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset]
        += (mi.CAPcbdb)*mi.numberParallel;

      dQdx[mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset]
        += (mi.CAPcbsb)*mi.numberParallel;


      // Row corresponding to the KCL for the drain prime node:
      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset]
        -= (+ mi.CAPcdgb + mi.CAPcddb + mi.CAPcdsb )*mi.numberParallel;

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquGateNodeOffset]
        += (mi.CAPcdgb) *mi.numberParallel;

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset]
        += (+ mi.CAPcddb )*mi.numberParallel;

      dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset]
        -= (- mi.CAPcdsb) *mi.numberParallel;

      // Row corresponding to the KCL for the source prime node:
      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquGateNodeOffset]
        += (mi.CAPcsgb) *mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset]
        -= (+ mi.CAPcsgb + mi.CAPcsdb + mi.CAPcssb) *mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset]
        -= (- mi.CAPcsdb) *mi.numberParallel;

      dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset]
        += (+ mi.CAPcssb) *mi.numberParallel;

      // Row associated with the charge equation
      // This is currently not supported.
      if (mi.nqsMod)
      {
        string msg;
        msg = "Master::loadDAEMatrices";
        msg += " nqsMod=1 is not ready yet.  Re-run with nqsMod=0\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
      }
    }
  }
  return true;
}
#endif

} // namespace MOSFET_B3
} // namespace Device
} // namespace Xyce

