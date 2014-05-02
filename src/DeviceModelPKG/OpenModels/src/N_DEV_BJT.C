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
// Filename       : $RCSfile: N_DEV_BJT.C,v $
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
// Revision Number: $Revision: 1.381.2.4 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef Xyce_DEBUG_DEVICE
#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif
#endif

// ----------    Xyce Includes  ----------
#include <N_DEV_BJT.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {


namespace BJT {


void Traits::loadInstanceParameters(ParametricData<BJT::Instance> &p)
{
    // Set up configuration constants:
// Set up double precision variables:
    p.addPar ("AREA", 1.0, false, NO_DEP,
      &BJT::Instance::AREA, NULL,
      U_NONE, CAT_GEOMETRY, "Relative device area");

    p.addPar ("IC1",  0.0, false, NO_DEP,
      &BJT::Instance::icVBE,
      &BJT::Instance::IC_GIVEN,
      U_VOLT, CAT_VOLT, "Vector of initial values: Vbe, Vce. Vbe=IC1");

    p.addPar ("IC2",  0.0, false, NO_DEP,
      &BJT::Instance::icVCE,
      &BJT::Instance::IC_GIVEN,
      U_VOLT, CAT_VOLT, "Vector of initial values: Vbe, Vce. Vce=IC2");

    p.addPar ("TEMP", 0.0, false, TIME_DEP,
      &BJT::Instance::TEMP, NULL,
      STANDARD, CAT_UNKNOWN, "Device temperature");

    // Set up non-double precision variables:
    p.addPar ("OFF", false, false, NO_DEP,
            &BJT::Instance::OFF, NULL,
            U_LOGIC, CAT_VOLT, "Initial condition of no voltage drops accross device");

    p.addPar ("LAMBERTW", false, false, NO_DEP,
            &BJT::Instance::lambertWFlag, NULL,
            U_LOGIC, CAT_NONE, "Flag for toggling the use of the lambert-W function instead of exponentials.");

    p.makeVector ("IC", 2);
}

void Traits::loadModelParameters(ParametricData<BJT::Model> &p)
{
    // Set up double precision variables:
    p.addPar ("TNOM", 0.0, false, NO_DEP,
	      &BJT::Model::TNOM, NULL,
       U_DEGC, CAT_UNKNOWN, "Parameter measurement temperature");

    p.addPar ("IS", 1.e-16, false, LOG_T_DEP,
	      &BJT::Model::satCur, NULL,
       U_AMP, CAT_UNKNOWN, "Transport saturation current");

    p.addPar ("BF", 100.0, true, NO_DEP,
	      &BJT::Model::betaF,
	      &BJT::Model::BFgiven,
       U_NONE, CAT_UNKNOWN, "Ideal maximum foward beta");

    p.addPar ("BFM", 100.0, true, NO_DEP,  // hspice
	      &BJT::Model::betaF,
	      &BJT::Model::BFMgiven,
       U_NONE, CAT_UNKNOWN, "Ideal maximum foward beta");

    p.addPar ("NF", 1.0, true, NO_DEP,
	      &BJT::Model::emissionCoeffF, NULL,
       U_NONE, CAT_UNKNOWN, "Foward current emission coefficient");

    // forward early voltage has aliases for different spice flavors.
    p.addPar ("VA", 0.0, false, NO_DEP,  // Hspice/Pspice
	      &BJT::Model::earlyVoltF,
	      &BJT::Model::VAgiven,
       U_VOLT, CAT_UNKNOWN, "Foward early voltage");

    p.addPar ("VBF", 0.0, false, NO_DEP,  // Hspice
	      &BJT::Model::earlyVoltF,
	      &BJT::Model::VBFgiven,
       U_VOLT, CAT_UNKNOWN, "Foward early voltage");

    p.addPar ("VAF", 0.0, false, NO_DEP,  // spice3
	      &BJT::Model::earlyVoltF,
	      &BJT::Model::VAFgiven,
       U_VOLT, CAT_UNKNOWN, "Foward early voltage");

    p.addPar ("IKF", 0.0, false, NO_DEP,
	      &BJT::Model::rollOffF,
	      &BJT::Model::IKFgiven,
       U_AMP, CAT_UNKNOWN, "Corner for foward-beta high-current roll-off");

    p.addPar ("IK", 0.0, false, NO_DEP,
	      &BJT::Model::rollOffF,
	      &BJT::Model::IKgiven,
       U_AMP, CAT_UNKNOWN, "Corner for foward-beta high-current roll-off");

    p.addPar ("JBF", 0.0, false, NO_DEP,
	      &BJT::Model::rollOffF,
	      &BJT::Model::JBFgiven,
       U_AMP, CAT_UNKNOWN, "Corner for foward-beta high-current roll-off");

    // base-emitter leakage sat current has aliases for different spice flavors.
    p.addPar ("ISE", 0.0, false, LOG_T_DEP,
	      &BJT::Model::leakBECurrent,
	      &BJT::Model::leakBECurrentGiven,
       U_AMP, CAT_UNKNOWN, "Base-emitter leakage saturation current");

    p.addPar ("JLE", 0.0, false, LOG_T_DEP,  // Hspice
	      &BJT::Model::leakBECurrent,
	      &BJT::Model::JLEgiven,
       U_AMP, CAT_UNKNOWN, "Base-emitter leakage saturation current");

    p.addPar ("NE", 1.5, false, NO_DEP,
	      &BJT::Model::leakBEEmissionCoeff,
	      &BJT::Model::NEgiven,
       U_NONE, CAT_UNKNOWN, "Base-emitter leakage emission coefficient");

    p.addPar ("NLE", 1.5, false, NO_DEP,
	      &BJT::Model::leakBEEmissionCoeff,
	      &BJT::Model::NLEgiven,
       U_NONE, CAT_UNKNOWN, "Base-emitter leakage emission coefficient");

    // BR/BRM alias:
    p.addPar ("BR", 1.0, false, NO_DEP,
	      &BJT::Model::betaR,
	      &BJT::Model::BRgiven,
       U_NONE, CAT_UNKNOWN, "Ideal maximum reverse beta");

    p.addPar ("BRM", 1.0, false, NO_DEP,
	      &BJT::Model::betaR,
	      &BJT::Model::BRMgiven,
       U_NONE, CAT_UNKNOWN, "Ideal maximum reverse beta");

    p.addPar ("NR", 1.0, true, NO_DEP,
	      &BJT::Model::emissionCoeffR, NULL,
       U_NONE, CAT_UNKNOWN, "Reverse current emission coefficient");

    // reverse early voltage has several aliases:
    // Interestingly, more aliases than for the forward early voltage.
    p.addPar ("VAR", 0.0, false, NO_DEP,  // spice
	      &BJT::Model::earlyVoltR,
	      &BJT::Model::VARgiven,
       U_VOLT, CAT_UNKNOWN, "Reverse early voltage");

    p.addPar ("VB", 0.0, false, NO_DEP,  // pspice, hspice
	      &BJT::Model::earlyVoltR,
	      &BJT::Model::VBgiven,
       U_VOLT, CAT_UNKNOWN, "Reverse early voltage");

    p.addPar ("VRB", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::earlyVoltR,
	      &BJT::Model::VRBgiven,
       U_VOLT, CAT_UNKNOWN, "Reverse early voltage");

    p.addPar ("BV", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::earlyVoltR,
	      &BJT::Model::BVgiven,
       U_VOLT, CAT_UNKNOWN, "Reverse early voltage");

    // Corner for reverse-beta high-current roll-off has an alias:
    p.addPar ("IKR", 0.0, false, NO_DEP,
	      &BJT::Model::rollOffR,
	      &BJT::Model::IKRgiven,
       U_AMP, CAT_UNKNOWN, "Corner for reverse-beta high-current roll-off");

    p.addPar ("JBR", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::rollOffR,
	      &BJT::Model::JBRgiven,
       U_AMP, CAT_UNKNOWN, "Corner for reverse-beta high-current roll-off");

    // Base-collector leakage saturation current has an alias:
    p.addPar ("ISC", 0.0, false, NO_DEP,  //spice/pspice
	      &BJT::Model::leakBCCurrent,
	      &BJT::Model::leakBCCurrentGiven,
       U_AMP, CAT_UNKNOWN, "Base-collector leakage saturation current");

    p.addPar ("JLC", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::leakBCCurrent,
	      &BJT::Model::JLCgiven,
       U_AMP, CAT_UNKNOWN, "Base-collector leakage saturation current");

    p.addPar ("NC", 2.0, false, NO_DEP,
	      &BJT::Model::leakBCEmissionCoeff, NULL,
       U_NONE, CAT_UNKNOWN, "Base-collector leakage emission coefficient");

    p.addPar ("RB", 0.0, false, NO_DEP,
	      &BJT::Model::baseResist, NULL,
       U_OHM, CAT_UNKNOWN, "Zero-bias (maximum) base resistance");

    // IRB (Current at which RB falls off by half)") has an alias:
    p.addPar ("IRB", 0.0, false, NO_DEP,  //spice/pspice
	      &BJT::Model::baseCurrHalfResist,
	      &BJT::Model::IRBgiven,
       U_AMP, CAT_UNKNOWN, "Current at which RB falls off by half");

    p.addPar ("JRB", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::baseCurrHalfResist,
	      &BJT::Model::JRBgiven,
       U_AMP, CAT_UNKNOWN, "Current at which RB falls off by half");

    p.addPar ("IOB", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::baseCurrHalfResist,
	      &BJT::Model::IOBgiven,
       U_AMP, CAT_UNKNOWN, "Current at which RB falls off by half");

    p.addPar ("RBM", 0.0, false, NO_DEP,
	      &BJT::Model::minBaseResist,
	      &BJT::Model::minBaseResistGiven,
       U_OHM, CAT_UNKNOWN, "Maximum base resistance");

    p.addPar ("RE", 0.0, false, MIN_RES,
	      &BJT::Model::emitterResist, NULL,
       U_OHM, CAT_UNKNOWN, "Emitter ohmic resistance");

    p.addPar ("RC", 0.0, false, MIN_RES,
	      &BJT::Model::collectorResist, NULL,
       U_OHM, CAT_UNKNOWN, "Collector ohmic resistance");

    p.addPar ("CJE", 0.0, false, MIN_CAP,
	      &BJT::Model::depCapBE, NULL,
       U_FARAD, CAT_UNKNOWN, "Base-emitter zero-bias p-n capacitance");

    // alias for VJE:
    p.addPar ("VJE", 0.75, false, NO_DEP,  //spice
	      &BJT::Model::potBE,
	      &BJT::Model::VJEgiven,
       U_VOLT, CAT_UNKNOWN, "Base-emitter built-in potential");

    p.addPar ("PE", 0.75, false, NO_DEP,  //pspice/hspice
	      &BJT::Model::potBE,
	      &BJT::Model::PEgiven,
       U_VOLT, CAT_UNKNOWN, "Base-emitter built-in potential");

    // alias for MJE:
    p.addPar ("MJE", 0.33, false, NO_DEP,  //spice
	      &BJT::Model::juncExpBE,
	      &BJT::Model::MJEgiven,
       U_NONE, CAT_UNKNOWN, "Base-emitter p-n grading factor");

    p.addPar ("ME", 0.33, false, NO_DEP,  //pspice/hspice
	      &BJT::Model::juncExpBE,
	      &BJT::Model::MEgiven,
       U_NONE, CAT_UNKNOWN, "Base-emitter p-n grading factor");

    p.addPar ("TF", 0.0, false, NO_DEP,
	      &BJT::Model::transTimeF, NULL,
       U_SECOND, CAT_UNKNOWN, "Ideal foward transit time");

    p.addPar ("XTF", 0.0, false, NO_DEP,
	      &BJT::Model::transTimeBiasCoeffF, NULL,
       U_NONE, CAT_UNKNOWN, "Transit time bias dependence coefficient");

    p.addPar ("VTF", 0.0, false, NO_DEP,
	      &BJT::Model::transTimeFVBC, NULL,
       U_VOLT, CAT_UNKNOWN, "Transit time dependancy on Vbc");

    // alias for ITF:
    p.addPar ("ITF", 0.0, false, NO_DEP,  // spice
	      &BJT::Model::transTimeHighCurrF,
	      &BJT::Model::ITFgiven,
       U_UNKNOWN, CAT_UNKNOWN, "Transit time dependancy on IC");

    p.addPar ("JTF", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::transTimeHighCurrF,
	      &BJT::Model::JTFgiven,
       U_UNKNOWN, CAT_UNKNOWN, "Transit time dependancy on IC");

    p.addPar ("PTF", 0.0, false, NO_DEP,
	      &BJT::Model::excessPhase, NULL,
       U_DEGREE, CAT_UNKNOWN, "Excess Phase at 1/(2pi*TF) Hz");

    p.addPar ("CJC", 0.0, false, MIN_CAP,
	      &BJT::Model::depCapBC, NULL,
       U_FARAD, CAT_UNKNOWN, "Base-collector zero-bias p-n capacitance");

    // alias for VJC:
    p.addPar ("VJC", 0.75, false, NO_DEP,  // spice
	      &BJT::Model::potBC,
	      &BJT::Model::VJCgiven,
       U_VOLT, CAT_UNKNOWN, "Base-collector built-in potential");

    p.addPar ("PC", 0.75, false, NO_DEP,  // pspice, hspice
	      &BJT::Model::potBC,
	      &BJT::Model::PCgiven,
       U_VOLT, CAT_UNKNOWN, "Base-collector built-in potential");

    // alias for MJC:
    p.addPar ("MJC", 0.33, false, NO_DEP, // spice
	      &BJT::Model::juncExpBC,
	      &BJT::Model::MJCgiven,
       U_NONE, CAT_UNKNOWN, "Base-collector p-n grading factor");

    p.addPar ("MC", 0.33, false, NO_DEP, //pspice/hspice
	      &BJT::Model::juncExpBC,
	      &BJT::Model::MCgiven,
       U_NONE, CAT_UNKNOWN, "Base-collector p-n grading factor");

    // alias for XCJC:
    p.addPar ("XCJC", 1.0, false, NO_DEP,
	      &BJT::Model::baseFracBCCap,
	      &BJT::Model::XCJCgiven,
       U_NONE, CAT_UNKNOWN, "Fraction of CJC connected internally to RB");
    p.addPar ("CDIS", 1.0, false, NO_DEP,
	      &BJT::Model::baseFracBCCap,
	      &BJT::Model::CDISgiven,
       U_NONE, CAT_UNKNOWN, "Fraction of CJC connected internally to RB");

    p.addPar ("TR", 0.0, false, NO_DEP,
	      &BJT::Model::transTimeR, NULL,
       U_SECOND, CAT_UNKNOWN, "Ideal reverse transit time");

    // substrate zero bias capacitance has several aliases for different spice flavors.
    p.addPar ("CJS", 0.0, false, NO_DEP,  // spice
	      &BJT::Model::CJS,
	      &BJT::Model::CJSgiven,
       U_FARAD, CAT_UNKNOWN, "Substrate zero-bias p-n capacitance");

    p.addPar ("CCS", 0.0, false, NO_DEP,  // Pspice/Hspice
	      &BJT::Model::CJS,
	      &BJT::Model::CCSgiven,
       U_FARAD, CAT_UNKNOWN, "Substrate zero-bias p-n capacitance");

    p.addPar ("CSUB", 0.0, false, NO_DEP,  // Hspice
	      &BJT::Model::CJS,
	      &BJT::Model::CSUBgiven,
       U_FARAD, CAT_UNKNOWN, "Substrate zero-bias p-n capacitance");

    // aliases for VJS:
    p.addPar ("VJS", 0.75, false, NO_DEP,  // spice
	      &BJT::Model::potSubst,
	      &BJT::Model::VJSgiven,
       U_VOLT, CAT_UNKNOWN, "Substrate built-in potential");

    p.addPar ("PS", 0.75, false, NO_DEP,  // pspice
	      &BJT::Model::potSubst,
	      &BJT::Model::PSgiven,
       U_VOLT, CAT_UNKNOWN, "Substrate built-in potential");

    p.addPar ("PSUB", 0.75, false, NO_DEP,  // hspice
	      &BJT::Model::potSubst,
	      &BJT::Model::PSUBgiven,
       U_VOLT, CAT_UNKNOWN, "Substrate built-in potential");

    // aliases for MJS:
    p.addPar ("MJS", 0.0, false, NO_DEP,
	      &BJT::Model::expSubst,
	      &BJT::Model::MJSgiven,
       U_NONE, CAT_UNKNOWN, "Substrate p-n grading factor");

    p.addPar ("MS", 0.0, false, NO_DEP,
	      &BJT::Model::expSubst,
	      &BJT::Model::MSgiven,
       U_NONE, CAT_UNKNOWN, "Substrate p-n grading factor");

    p.addPar ("ESUB", 0.0, false, NO_DEP,
	      &BJT::Model::expSubst,
	      &BJT::Model::ESUBgiven,
       U_NONE, CAT_UNKNOWN, "Substrate p-n grading factor");

    // alias for XTB:
    p.addPar ("XTB", 0.0, false, NO_DEP,  //spice/pspice
	      &BJT::Model::betaExp,
	      &BJT::Model::XTBgiven,
       U_NONE, CAT_UNKNOWN, "Foward and reverse beta temperature coefficient");

    p.addPar ("TB", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::betaExp,
	      &BJT::Model::TBgiven,
       U_NONE, CAT_UNKNOWN, "Foward and reverse beta temperature coefficient");

    p.addPar ("TCB", 0.0, false, NO_DEP,  // hspice
	      &BJT::Model::betaExp,
	      &BJT::Model::TCBgiven,
       U_NONE, CAT_UNKNOWN, "Foward and reverse beta temperature coefficient");

    p.addPar ("EG", 1.11, false, NO_DEP,
	      &BJT::Model::energyGap, NULL,
       U_EV, CAT_UNKNOWN, "Bandgap voltage (barrier highth)");

    // alias for XTI:
    p.addPar ("XTI", 3.0, false, NO_DEP,  // spice/hspice
	      &BJT::Model::tempExpIS,
	      &BJT::Model::XTIgiven,
       U_UNKNOWN, CAT_UNKNOWN, "Temperature exponent for IS. (synonymous with PT)");

    p.addPar ("PT", 3.0, false, NO_DEP,  // pspice
	      &BJT::Model::tempExpIS,
	      &BJT::Model::PTgiven,
       U_UNKNOWN, CAT_UNKNOWN, "Temperature exponent for IS. (synonymous with XTI)");

    p.addPar ("KF", 0.0, false, NO_DEP,
	      &BJT::Model::fNCoeff, NULL,
       U_NONE, CAT_UNKNOWN, "Flicker noise coefficient");

    p.addPar ("AF", 1.0, false, NO_DEP,
	      &BJT::Model::fNExp, NULL,
       U_UNKNOWN, CAT_UNKNOWN, "Flicker noise exponent");

    p.addPar ("FC", 0.5, false, NO_DEP,
	      &BJT::Model::depCapCoeff, NULL,
       U_NONE, CAT_UNKNOWN, "Foward-bias depletion capacitor coefficient");

    p.addPar ("C2", 0.0, false, NO_DEP,
	      &BJT::Model::c2,
	      &BJT::Model::c2Given,
       U_UNKNOWN, CAT_UNKNOWN, "Coefficient for base-emitter leak current.");

    p.addPar ("C4", 0.0, false, NO_DEP,
	      &BJT::Model::c4,
	      &BJT::Model::c4Given,
       U_UNKNOWN, CAT_UNKNOWN, "Coefficient for base-collector leak current.");

    // alias for NK/NKF:
    p.addPar ("NK", 0.5, false, NO_DEP,
	      &BJT::Model::rollOffExp,
	      &BJT::Model::NKgiven,
       U_NONE, CAT_UNKNOWN, "High current rolloff coefficient");

    p.addPar ("NKF", 0.5, false, NO_DEP,
	      &BJT::Model::rollOffExp,
	      &BJT::Model::NKFgiven,
       U_NONE, CAT_UNKNOWN, "High current rolloff coefficient");

    // Set up non-double precision variables:

    // Thermal model setup:
    DeviceModel::initThermalModel(p);
}



std::vector< std::vector<int> > Instance::jacStamp_RB_RC_RE_;
std::vector< std::vector<int> > Instance::jacStamp_RB_RC_;
std::vector< std::vector<int> > Instance::jacStamp_RB_RE_;
std::vector< std::vector<int> > Instance::jacStamp_RC_RE_;
std::vector< std::vector<int> > Instance::jacStamp_RB_;
std::vector< std::vector<int> > Instance::jacStamp_RC_;
std::vector< std::vector<int> > Instance::jacStamp_RE_;
std::vector< std::vector<int> > Instance::jacStamp_;

std::vector<int> Instance::jacMap_RB_RC_RE_;
std::vector<int> Instance::jacMap_RB_RC_;
std::vector<int> Instance::jacMap_RB_RE_;
std::vector<int> Instance::jacMap_RC_RE_;
std::vector<int> Instance::jacMap_RB_;
std::vector<int> Instance::jacMap_RC_;
std::vector<int> Instance::jacMap_RE_;
std::vector<int> Instance::jacMap_;

std::vector< std::vector<int> > Instance::jacMap2_RB_RC_RE_;
std::vector< std::vector<int> > Instance::jacMap2_RB_RC_;
std::vector< std::vector<int> > Instance::jacMap2_RB_RE_;
std::vector< std::vector<int> > Instance::jacMap2_RC_RE_;
std::vector< std::vector<int> > Instance::jacMap2_RB_;
std::vector< std::vector<int> > Instance::jacMap2_RC_;
std::vector< std::vector<int> > Instance::jacMap2_RE_;
std::vector< std::vector<int> > Instance::jacMap2_;

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  updateTemperature(TEMP);
  return true;
}

//----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra
// Creation Date : 11/30/00
//----------------------------------------------------------------------------

Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       model,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(model),
    AREA(1.0),
    icVBE(0.0),
    icVCE(0.0),
    TEMP(300.0),
    OFF(false),
    lambertWFlag(false),
    IC_GIVEN(false),
    externalNodeMode(false),
    offFlag(false),
    vt          (0.0),
    tSatCur     (0.0),
    tBetaF      (0.0),
    tBetaR      (0.0),
    tBELeakCur  (0.0),
    tBCLeakCur  (0.0),
    tBECap      (0.0),
    tBCCap      (0.0),
    tBEPot      (0.0),
    tBCPot      (0.0),
    tDepCap     (0.0),
    tF1         (0.0),
    tF4         (0.0),
    tF5         (0.0),
    tVCrit      (0.0),
    tleakBEEmissionCoeff(1.5),
    tleakBCEmissionCoeff(2.0),
    tRollOffExp (0.5),
    tInvRollOffF(0.0),
    tInvRollOffR(0.0),
    tInvEarlyVoltF(0.0),
    tInvEarlyVoltR(0.0),
    vEEp        (0.0),
    vBBp        (0.0),
    vCCp        (0.0),
    vBE         (0.0),
    vBC         (0.0),
    vBX         (0.0),
    vCS         (0.0),
    vBE_old     (0.0),
    vBC_old     (0.0),
    vBE_orig    (0.0),
    vBC_orig    (0.0),
    qB          (0.0),
    invqB       (0.0),
    dqBdvEp     (0.0),
    dqBdvBp     (0.0),
    dqBdvCp     (0.0),
    iBE         (0.0),
    iBC         (0.0),
    iBEleak     (0.0),
    iBCleak     (0.0),
    iCE         (0.0),
    iB          (0.0),
    iC          (0.0),
    iE          (0.0),
    iBEhighCurr (0.0),
    gBEhighCurr (0.0),
    gBE         (0.0),
    gBC         (0.0),
    gBEleak     (0.0),
    gBCleak     (0.0),
    gEpr        (0.0),
    gCpr        (0.0),
    gX          (0.0),
    geqCB       (0.0),
    capeqCB     (0.0),
    diBrdvB    (0.0),
    diBrdvEp   (0.0),
    diBrdvCp   (0.0),
    diBrdvBp   (0.0),
    diCEdvEp   (0.0),
    diCEdvCp   (0.0),
    diCEdvBp   (0.0),

    diBEdvEp   (0.0),
    diBEdvCp   (0.0),
    diBEdvBp   (0.0),

    gCapBEdiff (0.0),
    gCapBEdep  (0.0),
    gCapBCdiff (0.0),
    gCapBCdep  (0.0),
    gCapBX     (0.0),
    gCapCS     (0.0),
    gBEtot     (0.0),
    gBCtot     (0.0),
    qBEdiff     (0.0),
    iBEdiff     (0.0),
    capBEdiff   (0.0),
    qBEdep      (0.0),
    iBEdep      (0.0),
    capBEdep    (0.0),
    qCS         (0.0),
    iCS         (0.0),
    capCS       (0.0),
    qBCdiff     (0.0),
    iBCdiff     (0.0),
    capBCdiff   (0.0),
    qBCdep      (0.0),
    iBCdep      (0.0),
    capBCdep    (0.0),
    qBX         (0.0),
    iBX         (0.0),
    capBX       (0.0),


    li_Coll(-1),
    li_CollP(-1),
    li_Base(-1),
    li_BaseP(-1),
    li_Emit(-1),
    li_EmitP(-1),
    li_Subst(-1),
    li_Ifx(-1),
    li_dIfx(-1),
    li_qstateBEdiff(-1),
    li_qstateBEdep(-1),
    li_qstateCS(-1),
    li_qstateBCdiff(-1),
    li_qstateBCdep(-1),
    li_qstateBX(-1),
    li_istateCEXBC(-1),

    li_storevBE(-1),
    li_storevBC(-1),
    li_store_capeqCB(-1),
    li_store_dev_ib(-1),
    li_store_dev_ie(-1),
    li_store_dev_ic(-1),
    li_store_dev_is(-1),

    gcpr     (0.0),
    gepr     (0.0),
    gx       (0.0),
    gm       (0.0),
    go       (0.0),
    gmu      (0.0),
    gpi      (0.0),
    gccs     (0.0),
    geqbx    (0.0),
    geqbc    (0.0),
    nextCexbc(0.0),
    currCexbc(0.0),
    lastCexbc(0.0),
    phaseScalar(0.0),
    dt0        (0.0),
    dt1        (0.0),
    AEmitEquEmitPNodeOffset(-1),
    AEmitPEquEmitNodeOffset(-1),
    ABaseEquBasePNodeOffset(-1),
    ABasePEquBaseNodeOffset(-1),
    ACollEquCollPNodeOffset(-1),
    ACollPEquCollNodeOffset(-1),
    AEmitEquEmitNodeOffset(-1),
    AEmitPEquEmitPNodeOffset(-1),
    ABaseEquBaseNodeOffset(-1),
    ABasePEquBasePNodeOffset(-1),
    ACollEquCollNodeOffset(-1),
    ACollPEquCollPNodeOffset(-1),
    AEmitPEquBasePNodeOffset(-1),
    ABasePEquEmitPNodeOffset(-1),
    AEmitPEquCollPNodeOffset(-1),
    ACollPEquEmitPNodeOffset(-1),
    ABasePEquCollPNodeOffset(-1),
    ACollPEquBasePNodeOffset(-1),
    ABaseEquCollPNodeOffset(-1),
    ACollPEquBaseNodeOffset(-1),
    ASubstEquSubstNodeOffset(-1),
    ASubstEquCollPNodeOffset(-1),
    ACollPEquSubstNodeOffset(-1),
    ABaseEquEmitPNodeOffset(-1) ,
    ACollPEquIfxNodeOffset(-1),
    AEmitPEquIfxNodeOffset(-1),
    AIfxEquCollPNodeOffset(-1),
    AIfxEquBasePNodeOffset(-1),
    AIfxEquEmitPNodeOffset(-1),
    AIfxEquIfxNodeOffset(-1),
    AIfxEqudIfxNodeOffset(-1),
    AdIfxEquCollPNodeOffset(-1),
    AdIfxEquBasePNodeOffset(-1),
    AdIfxEquEmitPNodeOffset(-1),
    AdIfxEquIfxNodeOffset(-1),
    AdIfxEqudIfxNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // f-matrix pointers:
    f_EmitEquEmitPNodePtr(0),
    f_EmitPEquEmitNodePtr(0),
    f_BaseEquBasePNodePtr(0),
    f_BasePEquBaseNodePtr(0),
    f_CollEquCollPNodePtr(0),
    f_CollPEquCollNodePtr(0),
    f_EmitEquEmitNodePtr(0),
    f_EmitPEquEmitPNodePtr(0),
    f_BaseEquBaseNodePtr(0),
    f_BasePEquBasePNodePtr(0),
    f_CollEquCollNodePtr(0),
    f_CollPEquCollPNodePtr(0),
    f_EmitPEquBasePNodePtr(0),
    f_BasePEquEmitPNodePtr(0),
    f_EmitPEquCollPNodePtr(0),
    f_CollPEquEmitPNodePtr(0),
    f_BasePEquCollPNodePtr(0),
    f_CollPEquBasePNodePtr(0),
    f_BaseEquCollPNodePtr(0),
    f_CollPEquBaseNodePtr(0),
    f_SubstEquSubstNodePtr(0),
    f_SubstEquCollPNodePtr(0),
    f_CollPEquSubstNodePtr(0),
    f_BaseEquEmitPNodePtr(0),

    //new offsets for full integration of excess phase term
    f_CollPEquIfxNodePtr(0),
    f_EmitPEquIfxNodePtr(0),

    // ERK.  These 3 are only needed for dcop.
    f_IfxEquCollPNodePtr(0),
    f_IfxEquBasePNodePtr(0),
    f_IfxEquEmitPNodePtr(0),

    f_IfxEquIfxNodePtr(0),
    f_IfxEqudIfxNodePtr(0),

    f_dIfxEquCollPNodePtr(0),
    f_dIfxEquBasePNodePtr(0),
    f_dIfxEquEmitPNodePtr(0),
    f_dIfxEquIfxNodePtr(0),
    f_dIfxEqudIfxNodePtr(0),


    // q-matrix pointers:
    q_EmitEquEmitPNodePtr(0),
    q_EmitPEquEmitNodePtr(0),
    q_BaseEquBasePNodePtr(0),
    q_BasePEquBaseNodePtr(0),
    q_CollEquCollPNodePtr(0),
    q_CollPEquCollNodePtr(0),
    q_EmitEquEmitNodePtr(0),
    q_EmitPEquEmitPNodePtr(0),
    q_BaseEquBaseNodePtr(0),
    q_BasePEquBasePNodePtr(0),
    q_CollEquCollNodePtr(0),
    q_CollPEquCollPNodePtr(0),
    q_EmitPEquBasePNodePtr(0),
    q_BasePEquEmitPNodePtr(0),
    q_EmitPEquCollPNodePtr(0),
    q_CollPEquEmitPNodePtr(0),
    q_BasePEquCollPNodePtr(0),
    q_CollPEquBasePNodePtr(0),
    q_BaseEquCollPNodePtr(0),
    q_CollPEquBaseNodePtr(0),
    q_SubstEquSubstNodePtr(0),
    q_SubstEquCollPNodePtr(0),
    q_CollPEquSubstNodePtr(0),
    q_BaseEquEmitPNodePtr(0),

    //new offsets for full integration of excess phase term
    q_CollPEquIfxNodePtr(0),
    q_EmitPEquIfxNodePtr(0),

    // ERK.  These 3 are only needed for dcop.
    q_IfxEquCollPNodePtr(0),
    q_IfxEquBasePNodePtr(0),
    q_IfxEquEmitPNodePtr(0),

    q_IfxEquIfxNodePtr(0),
    q_IfxEqudIfxNodePtr(0),

    q_dIfxEquCollPNodePtr(0),
    q_dIfxEquBasePNodePtr(0),
    q_dIfxEquEmitPNodePtr(0),
    q_dIfxEquIfxNodePtr(0),
    q_dIfxEqudIfxNodePtr(0),
#endif

    callsOutputPlot (0)

{
  // This basically says, numExtVars = the number of vars specified
  // by the user on the netlist instance line.
  numExtVars   = IB.numExtVars;

  numStateVars = 7;
  setNumStoreVars(3);
  numLeadCurrentStoreVars = 4; // lead current dev_ib, dev_ie, dev_ic and dev_is

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 1;
  devConMap[2] = 1;
  devConMap[3] = 2;

  if (numExtVars > 4)
  {
    devConMap.resize(numExtVars);
  }
  for (int i1=4;i1<numExtVars;++i1)
  {
    devConMap[i1] = 1; // map all the extra ones to the base.
                       // That way I don't have to figure out
                       // which optional nodes have been set
                       // by the user.
  }

  // set up the jacStamp:
  if( jacStamp_.empty() )
  {
    jacStamp_RB_RC_RE_.resize(7);
    jacStamp_RB_RC_RE_[0].resize(2);   // Collector Row
    jacStamp_RB_RC_RE_[0][0] = 0;      // C-C
    jacStamp_RB_RC_RE_[0][1] = 4;      // C-C'

    jacStamp_RB_RC_RE_[1].resize(4);   // Base row
    jacStamp_RB_RC_RE_[1][0] = 1;      // B-B
    jacStamp_RB_RC_RE_[1][1] = 4;      // B-C'
    jacStamp_RB_RC_RE_[1][2] = 5;      // B-B'
    jacStamp_RB_RC_RE_[1][3] = 6;      // B-E'

    jacStamp_RB_RC_RE_[2].resize(2);   // Emitter row
    jacStamp_RB_RC_RE_[2][0] = 2;     // E-E
    jacStamp_RB_RC_RE_[2][1] = 6;     // E-E'

    jacStamp_RB_RC_RE_[3].resize(2);  // Substrate row
    jacStamp_RB_RC_RE_[3][0] = 3;     // S-S
    jacStamp_RB_RC_RE_[3][1] = 4;     // S-C'

    jacStamp_RB_RC_RE_[4].resize(6);  // Collector'
    jacStamp_RB_RC_RE_[4][0] = 0;     // C'-C
    jacStamp_RB_RC_RE_[4][1] = 1;     // C'-B
    jacStamp_RB_RC_RE_[4][2] = 3;     // C'-S
    jacStamp_RB_RC_RE_[4][3] = 4;     // C'-C'
    jacStamp_RB_RC_RE_[4][4] = 5;     // C'-B'
    jacStamp_RB_RC_RE_[4][5] = 6;     // C'-E'

    jacStamp_RB_RC_RE_[5].resize(4);  // Base'
    jacStamp_RB_RC_RE_[5][0] = 1;     // B'-B
    jacStamp_RB_RC_RE_[5][1] = 4;     // B'-C'
    jacStamp_RB_RC_RE_[5][2] = 5;     // B'-B'
    jacStamp_RB_RC_RE_[5][3] = 6;     // B'-E'

    jacStamp_RB_RC_RE_[6].resize(4);  // Emitter'
    jacStamp_RB_RC_RE_[6][0] = 2;     // E'-E
    jacStamp_RB_RC_RE_[6][1] = 4;     // E'-C'
    jacStamp_RB_RC_RE_[6][2] = 5;     // E'-B'
    jacStamp_RB_RC_RE_[6][3] = 6;     // E'-E'

    if (getDeviceOptions().newExcessPhase)
    {
      //Excess Phase Terms
      jacStamp_RB_RC_RE_.resize(9);

      jacStamp_RB_RC_RE_[4].resize(7);   // Collector' row
      jacStamp_RB_RC_RE_[4][6] = 7;      // C'-Ifx

      jacStamp_RB_RC_RE_[6].resize(5);   // Emitter' row
      jacStamp_RB_RC_RE_[6][4] = 7;      // E'-Ifx

      jacStamp_RB_RC_RE_[7].resize(5);   // Ifx row
      jacStamp_RB_RC_RE_[7][0] = 4;      // dIfx-C'
      jacStamp_RB_RC_RE_[7][1] = 5;      // dIfx-B'
      jacStamp_RB_RC_RE_[7][2] = 6;      // dIfx-E'
      jacStamp_RB_RC_RE_[7][3] = 7;      // Ifx-Ifx
      jacStamp_RB_RC_RE_[7][4] = 8;      // Ifx-dIfx

      jacStamp_RB_RC_RE_[8].resize(5);   // dIfx row
      jacStamp_RB_RC_RE_[8][0] = 4;      // dIfx-C'
      jacStamp_RB_RC_RE_[8][1] = 5;      // dIfx-B'
      jacStamp_RB_RC_RE_[8][2] = 6;      // dIfx-E'
      jacStamp_RB_RC_RE_[8][3] = 7;      // dIfx-Ifx
      jacStamp_RB_RC_RE_[8][4] = 8;      // dIfx-dIfx
    }

    jacMap_RB_RC_RE_.clear();

    int StampSize = 0;
    if (getDeviceOptions().newExcessPhase)
    {
      StampSize = 9;
    }
    else
    {
      StampSize = 7;
    }

    jacStampMap(jacStamp_RB_RC_RE_, jacMap_RB_RC_RE_, jacMap2_RB_RC_RE_,
                jacStamp_RB_RE_,    jacMap_RB_RE_, jacMap2_RB_RE_, 4, 0, StampSize );

    jacStampMap(jacStamp_RB_RC_RE_, jacMap_RB_RC_RE_, jacMap2_RB_RC_RE_,
                jacStamp_RC_RE_,    jacMap_RC_RE_, jacMap2_RC_RE_, 5, 1, StampSize );

    jacStampMap(jacStamp_RB_RC_RE_, jacMap_RB_RC_RE_, jacMap2_RB_RC_RE_,
                jacStamp_RB_RC_,    jacMap_RB_RC_, jacMap2_RB_RC_, 6, 2, StampSize );

    jacStampMap(jacStamp_RB_RC_, jacMap_RB_RC_, jacMap2_RB_RC_,
                jacStamp_RB_,    jacMap_RB_, jacMap2_RB_, 4, 0, StampSize );

    jacStampMap(jacStamp_RB_RE_, jacMap_RB_RE_, jacMap2_RB_RE_,
                jacStamp_RE_,    jacMap_RE_, jacMap2_RE_, 4, 1, StampSize );

    jacStampMap(jacStamp_RC_RE_, jacMap_RC_RE_, jacMap2_RC_RE_,
                jacStamp_RC_,    jacMap_RC_, jacMap2_RC_, 5, 2, StampSize );

    jacStampMap(jacStamp_RC_, jacMap_RC_, jacMap2_RC_,
                jacStamp_,    jacMap_, jacMap2_, 4, 0, StampSize );
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    TEMP = getDeviceOptions().temp.getImmutableValue<double>();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // This is tricky.  The prime nodes have 2 issues:
  //   1) do they exist?
  //   2) are they internal or external?
  //
  //   The default BJT behavior is for them to be internal,
  //   if they exist.
  //
  //   Also, node 4 is always "external" even if not specified
  //   by the user!  If the user doesn't specify it, it gets set
  //   to gnd, as an external node.
  //
  //   Nominal order in the netlist, after the original 4 nodes
  //   is collector', base', emitter'.
  //
  //   I've chosen to only allow 2 modes:  all optional variables
  //   are internal, or all optional nodes (that exist) are
  //   external.  No mixing and matching.  This is toggled by the
  //   flag, externalNodeMode.

  int cNode=0;
  int eNode=0;
  int bNode=0;

  cNode = ((model_.collectorResist==0.0)?0:1);
  eNode = ((model_.emitterResist==0.0)?0:1);
  bNode = ((model_.baseResist==0.0)?0:1);

  int numExist = cNode + eNode + bNode;

  if (numExtVars <= 4)
  {
    numExtVars = 4;
    externalNodeMode = false;
    numIntVars = cNode + eNode + bNode;
  }
  else
  {
    if (numExtVars != 4 + numExist)
    {
      std::string msg = ":";
      msg += "  name = " + getName();
      msg += " wrong number of external nodes are set!";
      msg += "  Either set none of them, or set them all.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }

    externalNodeMode = true;
    //numExtVars += cNode + eNode + bNode;  // This is already set in netlist.
  }


  if (getDeviceOptions().newExcessPhase)
  {
    //add in 2 variables for excess phase calc.
    numIntVars += 2;
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra
// Creation Date : 6/13/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                      const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "  BJTInstance::registerLIDs" <<std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
    Xyce::dout() << "  number of internal variables: " << numIntVars << std::endl;
    Xyce::dout() << "  number of external variables: " << numExtVars << std::endl;
    Xyce::dout() << "  numIntVars = " << numIntVars << std::endl;
    Xyce::dout() << "  numExtVars = " << numExtVars << std::endl;
  }
#endif

  // Copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "  Internal LID List" << std::endl;
    for( int i = 0; i < intLIDVec.size(); ++i )
      Xyce::dout() << "  " << intLIDVec[i] << std::endl;
    Xyce::dout() << "  External LID List" << std::endl;
    for( int i = 0; i < extLIDVec.size(); ++i )
      Xyce::dout() << "  " << extLIDVec[i] << std::endl;
  }
#endif

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.  For the matrix indices, first do the
  // rows.

  // First do external variables:
  // it1 at collector, it2 at collectorPrime
  int intIndex = 0;
  int extIndex = 0;

  li_Coll  = extLIDVec[extIndex++];
  li_Base  = extLIDVec[extIndex++];
  li_Emit  = extLIDVec[extIndex++];
  li_Subst = extLIDVec[extIndex++];

  if( model_.collectorResist == 0.0 )
  {
    li_CollP = li_Coll;
  }
  else
  {
    if (externalNodeMode)
      li_CollP = extLIDVec[extIndex++];
    else
      li_CollP = intLIDVec[intIndex++];
  }


  if( model_.baseResist == 0.0 )
  {
    li_BaseP = li_Base;
  }
  else
  {
    if (externalNodeMode)
      li_BaseP = extLIDVec[extIndex++];
    else
      li_BaseP = intLIDVec[intIndex++];
  }


  if( model_.emitterResist == 0.0 )
  {
    li_EmitP = li_Emit;
  }
  else
  {
    if (externalNodeMode)
      li_EmitP = extLIDVec[extIndex++];
    else
      li_EmitP = intLIDVec[intIndex++];
  }


  if (getDeviceOptions().newExcessPhase)
  {
    li_Ifx = intLIDVec[intIndex++];
    li_dIfx = intLIDVec[intIndex++];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "  li_Coll  = " << li_Coll << std::endl;
    Xyce::dout() << "  li_CollP = " << li_CollP << std::endl;
    Xyce::dout() << "  li_Base  = " << li_Base << std::endl;
    Xyce::dout() << "  li_BaseP = " << li_BaseP << std::endl;
    Xyce::dout() << "  li_Emit  = " << li_Emit << std::endl;
    Xyce::dout() << "  li_EmitP = " << li_EmitP << std::endl;
    Xyce::dout() << "  li_Subst = " << li_Subst << std::endl;

    if (getDeviceOptions().newExcessPhase)
    {
      Xyce::dout() << "  li_Ifx   = " << li_Ifx   << std::endl;
      Xyce::dout() << "  li_dIfx  = " << li_dIfx  << std::endl;
    }

    Xyce::dout() << section_divider << std::endl;
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
std::map<int,std::string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  std::string tmpstr;
  if (intNameMap.empty () && !externalNodeMode)
  {
    // Setup the internal node name map:
    if( model_.collectorResist != 0.0 )
    {
      tmpstr = getName()+"_collectorprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_CollP ] = tmpstr;
    }

    if( model_.baseResist != 0.0 )
    {
      tmpstr = getName()+"_baseprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_BaseP ] = tmpstr;
    }

    if( model_.emitterResist != 0.0 )
    {
      tmpstr = getName()+"_emitterprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_EmitP ] = tmpstr;
    }
  }

  if (getDeviceOptions().newExcessPhase)
  {
    tmpstr = getName()+"_ExcessPhase_Ifx";
    spiceInternalName (tmpstr);
    intNameMap [ li_Ifx ] = tmpstr;

    tmpstr = getName()+"_ExcessPhase_dIfx";
    spiceInternalName (tmpstr);
    intNameMap [ li_dIfx ] = tmpstr;
  }

  return intNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical Systems Modeling
// Creation Date : 01/23/13
//-----------------------------------------------------------------------------
std::map<int,std::string> & Instance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    std::string modName(getName());
    spiceInternalName(modName);
    std::string tmpstr;

    tmpstr = modName+":VBE";
    storeNameMap[ li_storevBE ] = tmpstr;
    tmpstr = modName+":VBC";
    storeNameMap[ li_storevBC ] = tmpstr;
    tmpstr = modName+":CAPEQCB";
    storeNameMap[ li_store_capeqCB ] = tmpstr;
    if( loadLeadCurrent )
    {
      tmpstr = modName+":DEV_IB";
      storeNameMap[ li_store_dev_ib ] = tmpstr;
      tmpstr = modName+":DEV_IE";
      storeNameMap[ li_store_dev_ie ] = tmpstr;
      tmpstr = modName+":DEV_IC";
      storeNameMap[ li_store_dev_ic ] = tmpstr;
      tmpstr = modName+":DEV_IS";
      storeNameMap[ li_store_dev_is ] = tmpstr;
    }
  }

  return storeNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra
// Creation Date : 6/13/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int i=0;
  li_qstateBEdiff = staLIDVec[i++];
  li_qstateBEdep = staLIDVec[i++];
  li_qstateCS = staLIDVec[i++];
  li_qstateBCdiff = staLIDVec[i++];
  li_qstateBCdep = staLIDVec[i++];
  li_qstateBX = staLIDVec[i++];
  li_istateCEXBC = staLIDVec[i++];
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/09/11
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int i=0;
  li_storevBE =  stoLIDVec[i++];
  li_storevBC =  stoLIDVec[i++];
  li_store_capeqCB =  stoLIDVec[i++];
  if( loadLeadCurrent )
  {
    li_store_dev_ib = stoLIDVec[i++];
    li_store_dev_ie = stoLIDVec[i++];
    li_store_dev_ic = stoLIDVec[i++];
    li_store_dev_is = stoLIDVec[i++];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra
// Creation Date : 8/21/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (model_.baseResist)
  {
    if (model_.collectorResist)
    {
      if (model_.emitterResist)
        return jacStamp_RB_RC_RE_;
      else
        return jacStamp_RB_RC_;
    }
    else
    {
      if (model_.emitterResist)
        return jacStamp_RB_RE_;
      else
        return jacStamp_RB_;
    }
  }
  else
  {
    if (model_.collectorResist)
    {
      if (model_.emitterResist)
        return jacStamp_RC_RE_;
      else
        return jacStamp_RC_;
    }
    else
    {
      if (model_.emitterResist)
        return jacStamp_RE_;
      else
        return jacStamp_;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra
// Creation Date : 8/27/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  std::vector<int> map;
  std::vector< std::vector<int> > map2;

  if (model_.baseResist)
  {
    if (model_.collectorResist)
    {
      if (model_.emitterResist)
      {
        map = jacMap_RB_RC_RE_;
        map2 = jacMap2_RB_RC_RE_;
      }
      else
      {
        map = jacMap_RB_RC_;
        map2 = jacMap2_RB_RC_;
      }
    }
    else
    {
      if (model_.emitterResist)
      {
        map = jacMap_RB_RE_;
        map2 = jacMap2_RB_RE_;
      }
      else
      {
        map = jacMap_RB_;
        map2 = jacMap2_RB_;
      }
    }
  }
  else
  {
    if (model_.collectorResist)
    {
      if (model_.emitterResist)
      {
        map = jacMap_RC_RE_;
        map2 = jacMap2_RC_RE_;
      }
      else
      {
        map = jacMap_RC_;
        map2 = jacMap2_RC_;
      }
    }
    else
    {
      if (model_.emitterResist)
      {
        map = jacMap_RE_;
        map2 = jacMap2_RE_;
      }
      else
      {
        map = jacMap_;
        map2 = jacMap2_;
      }
    }
  }

  ACollEquCollNodeOffset    = jacLIDVec[map[0]][map2[0][0]];
  ACollEquCollPNodeOffset   = jacLIDVec[map[0]][map2[0][1]];

  ABaseEquBaseNodeOffset    = jacLIDVec[map[1]][map2[1][0]];
  ABaseEquCollPNodeOffset   = jacLIDVec[map[1]][map2[1][1]];
  ABaseEquBasePNodeOffset   = jacLIDVec[map[1]][map2[1][2]];
  ABaseEquEmitPNodeOffset   = jacLIDVec[map[1]][map2[1][3]];

  AEmitEquEmitNodeOffset    = jacLIDVec[map[2]][map2[2][0]];
  AEmitEquEmitPNodeOffset   = jacLIDVec[map[2]][map2[2][1]];

  ASubstEquSubstNodeOffset  = jacLIDVec[map[3]][map2[3][0]];
  ASubstEquCollPNodeOffset  = jacLIDVec[map[3]][map2[3][1]];

  ACollPEquCollNodeOffset   = jacLIDVec[map[4]][map2[4][0]];
  ACollPEquBaseNodeOffset   = jacLIDVec[map[4]][map2[4][1]];
  ACollPEquSubstNodeOffset  = jacLIDVec[map[4]][map2[4][2]];
  ACollPEquCollPNodeOffset  = jacLIDVec[map[4]][map2[4][3]];
  ACollPEquBasePNodeOffset  = jacLIDVec[map[4]][map2[4][4]];
  ACollPEquEmitPNodeOffset  = jacLIDVec[map[4]][map2[4][5]];

  ABasePEquBaseNodeOffset   = jacLIDVec[map[5]][map2[5][0]];
  ABasePEquCollPNodeOffset  = jacLIDVec[map[5]][map2[5][1]];
  ABasePEquBasePNodeOffset  = jacLIDVec[map[5]][map2[5][2]];
  ABasePEquEmitPNodeOffset  = jacLIDVec[map[5]][map2[5][3]];

  AEmitPEquEmitNodeOffset   = jacLIDVec[map[6]][map2[6][0]];
  AEmitPEquCollPNodeOffset  = jacLIDVec[map[6]][map2[6][1]];
  AEmitPEquBasePNodeOffset  = jacLIDVec[map[6]][map2[6][2]];
  AEmitPEquEmitPNodeOffset  = jacLIDVec[map[6]][map2[6][3]];

  if (getDeviceOptions().newExcessPhase)
  {
    //Excess Phase Terms
    ACollPEquIfxNodeOffset    = jacLIDVec[map[4]][map2[4][6]];

    AEmitPEquIfxNodeOffset    = jacLIDVec[map[6]][map2[6][4]];

    AIfxEquCollPNodeOffset    = jacLIDVec[map[7]][map2[7][0]];
    AIfxEquBasePNodeOffset    = jacLIDVec[map[7]][map2[7][1]];
    AIfxEquEmitPNodeOffset    = jacLIDVec[map[7]][map2[7][2]];
    AIfxEquIfxNodeOffset      = jacLIDVec[map[7]][map2[7][3]];
    AIfxEqudIfxNodeOffset     = jacLIDVec[map[7]][map2[7][4]];

    AdIfxEquCollPNodeOffset   = jacLIDVec[map[8]][map2[8][0]];
    AdIfxEquBasePNodeOffset   = jacLIDVec[map[8]][map2[8][1]];
    AdIfxEquEmitPNodeOffset   = jacLIDVec[map[8]][map2[8][2]];
    AdIfxEquIfxNodeOffset     = jacLIDVec[map[8]][map2[8][3]];
    AdIfxEqudIfxNodeOffset    = jacLIDVec[map[8]][map2[8][4]];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/01/08
//-----------------------------------------------------------------------------
void Instance::setupPointers()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  // F-matrix
  f_CollEquCollNodePtr =     &(dFdx[li_Coll][ACollEquCollNodeOffset    ]);
  f_CollEquCollPNodePtr =    &(dFdx[li_Coll][ACollEquCollPNodeOffset   ]);

  f_BaseEquBaseNodePtr =     &(dFdx[li_Base][ABaseEquBaseNodeOffset    ]);
  f_BaseEquCollPNodePtr =    &(dFdx[li_Base][ABaseEquCollPNodeOffset   ]);
  f_BaseEquBasePNodePtr =    &(dFdx[li_Base][ABaseEquBasePNodeOffset   ]);
  f_BaseEquEmitPNodePtr =    &(dFdx[li_Base][ABaseEquEmitPNodeOffset   ]);

  f_EmitEquEmitNodePtr =     &(dFdx[li_Emit][AEmitEquEmitNodeOffset    ]);
  f_EmitEquEmitPNodePtr =    &(dFdx[li_Emit][AEmitEquEmitPNodeOffset   ]);

  f_SubstEquSubstNodePtr =   &(dFdx[li_Subst][ASubstEquSubstNodeOffset  ]);
  f_SubstEquCollPNodePtr =   &(dFdx[li_Subst][ASubstEquCollPNodeOffset  ]);

  f_CollPEquCollNodePtr =    &(dFdx[li_CollP][ACollPEquCollNodeOffset   ]);
  f_CollPEquBaseNodePtr =    &(dFdx[li_CollP][ACollPEquBaseNodeOffset   ]);
  f_CollPEquSubstNodePtr =   &(dFdx[li_CollP][ACollPEquSubstNodeOffset  ]);
  f_CollPEquCollPNodePtr =   &(dFdx[li_CollP][ACollPEquCollPNodeOffset  ]);
  f_CollPEquBasePNodePtr =   &(dFdx[li_CollP][ACollPEquBasePNodeOffset  ]);
  f_CollPEquEmitPNodePtr =   &(dFdx[li_CollP][ACollPEquEmitPNodeOffset  ]);

  f_BasePEquBaseNodePtr =    &(dFdx[li_BaseP][ABasePEquBaseNodeOffset   ]);
  f_BasePEquCollPNodePtr =   &(dFdx[li_BaseP][ABasePEquCollPNodeOffset  ]);
  f_BasePEquBasePNodePtr =   &(dFdx[li_BaseP][ABasePEquBasePNodeOffset  ]);
  f_BasePEquEmitPNodePtr =   &(dFdx[li_BaseP][ABasePEquEmitPNodeOffset  ]);

  f_EmitPEquEmitNodePtr =    &(dFdx[li_EmitP][AEmitPEquEmitNodeOffset   ]);
  f_EmitPEquCollPNodePtr =   &(dFdx[li_EmitP][AEmitPEquCollPNodeOffset  ]);
  f_EmitPEquBasePNodePtr =   &(dFdx[li_EmitP][AEmitPEquBasePNodeOffset  ]);
  f_EmitPEquEmitPNodePtr =   &(dFdx[li_EmitP][AEmitPEquEmitPNodeOffset  ]);


  // Q-matrix:
  q_CollEquCollNodePtr =     &(dQdx[li_Coll][ACollEquCollNodeOffset    ]);
  q_CollEquCollPNodePtr =    &(dQdx[li_Coll][ACollEquCollPNodeOffset   ]);

  q_BaseEquBaseNodePtr =     &(dQdx[li_Base][ABaseEquBaseNodeOffset    ]);
  q_BaseEquCollPNodePtr =    &(dQdx[li_Base][ABaseEquCollPNodeOffset   ]);
  q_BaseEquBasePNodePtr =    &(dQdx[li_Base][ABaseEquBasePNodeOffset   ]);
  q_BaseEquEmitPNodePtr =    &(dQdx[li_Base][ABaseEquEmitPNodeOffset   ]);

  q_EmitEquEmitNodePtr =     &(dQdx[li_Emit][AEmitEquEmitNodeOffset    ]);
  q_EmitEquEmitPNodePtr =    &(dQdx[li_Emit][AEmitEquEmitPNodeOffset   ]);

  q_SubstEquSubstNodePtr =   &(dQdx[li_Subst][ASubstEquSubstNodeOffset  ]);
  q_SubstEquCollPNodePtr =   &(dQdx[li_Subst][ASubstEquCollPNodeOffset  ]);

  q_CollPEquCollNodePtr =    &(dQdx[li_CollP][ACollPEquCollNodeOffset   ]);
  q_CollPEquBaseNodePtr =    &(dQdx[li_CollP][ACollPEquBaseNodeOffset   ]);
  q_CollPEquSubstNodePtr =   &(dQdx[li_CollP][ACollPEquSubstNodeOffset  ]);
  q_CollPEquCollPNodePtr =   &(dQdx[li_CollP][ACollPEquCollPNodeOffset  ]);
  q_CollPEquBasePNodePtr =   &(dQdx[li_CollP][ACollPEquBasePNodeOffset  ]);
  q_CollPEquEmitPNodePtr =   &(dQdx[li_CollP][ACollPEquEmitPNodeOffset  ]);

  q_BasePEquBaseNodePtr =    &(dQdx[li_BaseP][ABasePEquBaseNodeOffset   ]);
  q_BasePEquCollPNodePtr =   &(dQdx[li_BaseP][ABasePEquCollPNodeOffset  ]);
  q_BasePEquBasePNodePtr =   &(dQdx[li_BaseP][ABasePEquBasePNodeOffset  ]);
  q_BasePEquEmitPNodePtr =   &(dQdx[li_BaseP][ABasePEquEmitPNodeOffset  ]);

  q_EmitPEquEmitNodePtr =    &(dQdx[li_EmitP][AEmitPEquEmitNodeOffset   ]);
  q_EmitPEquCollPNodePtr =   &(dQdx[li_EmitP][AEmitPEquCollPNodeOffset  ]);
  q_EmitPEquBasePNodePtr =   &(dQdx[li_EmitP][AEmitPEquBasePNodeOffset  ]);
  q_EmitPEquEmitPNodePtr =   &(dQdx[li_EmitP][AEmitPEquEmitPNodeOffset  ]);

  if (getDeviceOptions().newExcessPhase)
  {
    //F-matrix excess phase Terms
    f_CollPEquIfxNodePtr =       &(dFdx[li_CollP][ACollPEquIfxNodeOffset    ]);

    f_EmitPEquIfxNodePtr =       &(dFdx[li_EmitP][AEmitPEquIfxNodeOffset    ]);

    f_IfxEquCollPNodePtr =       &(dFdx[li_Ifx][AIfxEquCollPNodeOffset    ]);
    f_IfxEquBasePNodePtr =       &(dFdx[li_Ifx][AIfxEquBasePNodeOffset    ]);
    f_IfxEquEmitPNodePtr =       &(dFdx[li_Ifx][AIfxEquEmitPNodeOffset    ]);
    f_IfxEquIfxNodePtr =         &(dFdx[li_Ifx][AIfxEquIfxNodeOffset      ]);
    f_IfxEqudIfxNodePtr =        &(dFdx[li_Ifx][AIfxEqudIfxNodeOffset     ]);

    f_dIfxEquCollPNodePtr =      &(dFdx[li_dIfx][AdIfxEquCollPNodeOffset   ]);
    f_dIfxEquBasePNodePtr =      &(dFdx[li_dIfx][AdIfxEquBasePNodeOffset   ]);
    f_dIfxEquEmitPNodePtr =      &(dFdx[li_dIfx][AdIfxEquEmitPNodeOffset   ]);
    f_dIfxEquIfxNodePtr =        &(dFdx[li_dIfx][AdIfxEquIfxNodeOffset     ]);
    f_dIfxEqudIfxNodePtr =       &(dFdx[li_dIfx][AdIfxEqudIfxNodeOffset    ]);


    // Q-matrix excess phase terms
    q_CollPEquIfxNodePtr =       &(dQdx[li_CollP][ACollPEquIfxNodeOffset    ]);

    q_EmitPEquIfxNodePtr =       &(dQdx[li_EmitP][AEmitPEquIfxNodeOffset    ]);

    q_IfxEquCollPNodePtr =       &(dQdx[li_Ifx][AIfxEquCollPNodeOffset    ]);
    q_IfxEquBasePNodePtr =       &(dQdx[li_Ifx][AIfxEquBasePNodeOffset    ]);
    q_IfxEquEmitPNodePtr =       &(dQdx[li_Ifx][AIfxEquEmitPNodeOffset    ]);
    q_IfxEquIfxNodePtr =         &(dQdx[li_Ifx][AIfxEquIfxNodeOffset      ]);
    q_IfxEqudIfxNodePtr =        &(dQdx[li_Ifx][AIfxEqudIfxNodeOffset     ]);

    q_dIfxEquCollPNodePtr =      &(dQdx[li_dIfx][AdIfxEquCollPNodeOffset   ]);
    q_dIfxEquBasePNodePtr =      &(dQdx[li_dIfx][AdIfxEquBasePNodeOffset   ]);
    q_dIfxEquEmitPNodePtr =      &(dQdx[li_dIfx][AdIfxEquEmitPNodeOffset   ]);
    q_dIfxEquIfxNodePtr =        &(dQdx[li_dIfx][AdIfxEquIfxNodeOffset     ]);
    q_dIfxEqudIfxNodePtr =       &(dQdx[li_dIfx][AdIfxEqudIfxNodeOffset    ]);
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra
// Creation Date : 11/30/00
//-----------------------------------------------------------------------------
bool Instance::updateTemperature( const double & temp )
{

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Start BJTInst::updateTemperature" << std::endl;
    Xyce::dout() << "temp = "<<temp << std::endl;
  }
#endif
  if( temp != -999.0 ) TEMP = temp;
  if (model_.interpolateTNOM(TEMP))
  {
    // make sure interpolation doesn't take any resistance negative
    if(model_.baseResist < 0) model_.baseResist = 0;
    if(model_.emitterResist < 0) model_.emitterResist = 0;
    if(model_.collectorResist < 0) model_.collectorResist = 0;

    // some params may have changed during interpolation
    model_.processParams();
  }

  //Generation of temperature based factors
  vt = TEMP * CONSTKoverQ;
  double TNOM  = model_.TNOM;
  double fact2 = TEMP / CONSTREFTEMP;
  double egfet = CONSTEg0 - ( CONSTalphaEg * TEMP * TEMP ) /
                 ( TEMP + CONSTbetaEg );
  double arg = -egfet / ( 2.0 * CONSTboltz * TEMP ) +
               CONSTEg300 / ( CONSTboltz * ( 2.0 * CONSTREFTEMP ) );
  double pbfact = -2.0 * vt * ( 1.5 * log( fact2 ) + CONSTQ * arg );
  double ratlog = log( TEMP / model_.TNOM );
  double ratio1 = TEMP / model_.TNOM - 1.0;
  double factlog = ratio1 * model_.energyGap / vt
                          + model_.tempExpIS * ratlog;
  double factor  = exp( factlog );
  double bfactor = exp( ratlog * model_.betaExp );

  // Temp. adj. saturation current
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " Factor is " << factor << std::endl;
    Xyce::dout() << " Bfactor is " << bfactor << std::endl;
    Xyce::dout() << " factlog is " << factlog << std::endl;
    Xyce::dout() << " ratio1 is " << ratio1 << std::endl;
    Xyce::dout() << "  TEMP is " << TEMP << std::endl;
    Xyce::dout() << "  TNOM is " << model_.TNOM << std::endl;
    Xyce::dout() << "   Energy gap is " << model_.energyGap << std::endl;
    Xyce::dout() << "   tempExpIS is " << model_.tempExpIS << std::endl;
  }
#endif

  // Temp. adj. zero-bias junction capacitances and built-in potentials
  double fact1 = model_.TNOM / CONSTREFTEMP;
  double pbo = ( model_.potBE - pbfact ) / fact1;
  double gmaold = ( model_.potBE - pbo ) / pbo;

  tBECap = model_.depCapBE / ( 1.0 + model_.juncExpBE *
           ( 4.e-4 * ( model_.TNOM - CONSTREFTEMP ) - gmaold ) );
  tBEPot = fact2 * pbo + pbfact;

  double gmanew = ( tBEPot - pbo ) / pbo;

  tBECap *= 1.0 + model_.juncExpBE *
            ( 4.e-4 * ( TEMP - CONSTREFTEMP ) - gmanew );

  pbo = ( model_.potBC - pbfact ) / fact1;
  gmaold = ( model_.potBC - pbo ) / pbo;

  tBCCap = model_.depCapBC / ( 1.0 + model_.juncExpBC *
           ( 4.e-4 * ( model_.TNOM - CONSTREFTEMP ) - gmaold ) );
  tBCPot = fact2 * pbo + pbfact;

  gmanew = ( tBCPot - pbo ) / pbo;

  tBCCap *= 1.0 + model_.juncExpBC *
            ( 4.e-4 * ( TEMP - CONSTREFTEMP ) - gmanew );

  // Temp. adj. depletion capacitance
  tDepCap = model_.depCapCoeff * tBEPot;

  // Temp. adj. polynomial factors
  double xfc = log( 1.0 - model_.depCapCoeff );
  tF1 = tBEPot * ( 1.0 - exp( ( 1.0 - model_.juncExpBE ) * xfc ) ) /
          ( 1.0 - model_.juncExpBE );

  tF4 = model_.depCapCoeff * tBCPot;
  tF5 = tBCPot * ( 1.0 - exp( ( 1.0 - model_.juncExpBC ) * xfc ) ) /
         ( 1.0 - model_.juncExpBC );

  // Temp. adj. critical voltage
  // tVCrit = vt * log( vt / ( 1.41414 * model_.satCur ) );
  tVCrit = vt * log( vt / ( CONSTroot2 * model_.satCur ) );

  tSatCur = model_.satCur * factor;
   // Temp. adj. beta's
  tBetaF = model_.betaF * bfactor;
  tBetaR = model_.betaR * bfactor;
   // Temp. adj. leakage currents
  if (!model_.leakBECurrentGiven && model_.c2Given)
    model_.leakBECurrent = model_.c2 * model_.satCur;
  if (!model_.leakBCCurrentGiven && model_.c4Given)
    model_.leakBCCurrent = model_.c4 * model_.satCur;
  tBELeakCur = model_.leakBECurrent * exp( factlog / model_.leakBEEmissionCoeff ) / bfactor;
  tBCLeakCur = model_.leakBCCurrent * exp( factlog / model_.leakBCEmissionCoeff ) / bfactor;

  tleakBEEmissionCoeff = model_.leakBEEmissionCoeff;
  tleakBCEmissionCoeff = model_.leakBCEmissionCoeff;
  tRollOffExp      = model_.rollOffExp;
  tBaseResist      = model_.baseResist;
  tCollectorResist = model_.collectorResist;
  tEmitterResist   = model_.emitterResist;

  if (model_.rollOffF != 0)
    tInvRollOffF = 1/model_.rollOffF;
  else
    tInvRollOffF = 0;
  if (model_.rollOffR != 0)
    tInvRollOffR = 1/model_.rollOffR;
  else
    tInvRollOffR = 0;
  if (model_.earlyVoltF != 0)
    tInvEarlyVoltF = 1/model_.earlyVoltF;
  else
    tInvEarlyVoltF = 0;
  if (model_.earlyVoltR != 0)
    tInvEarlyVoltR = 1/model_.earlyVoltR;
  else
    tInvEarlyVoltR = 0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {

    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "In BJTInst::updateTemperature" << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "TEMP       = " << TEMP << std::endl;
    Xyce::dout() << "TNOM       = " << TNOM << std::endl;
    Xyce::dout() << "vt         = " << vt << std::endl;
    Xyce::dout() << "tSatCur    = " << tSatCur << std::endl;
    Xyce::dout() << "tBetaF     = " << tBetaF << std::endl;
    Xyce::dout() << "tBetaR     = " << tBetaR << std::endl;
    Xyce::dout() << "tBELeakCur = " << tBELeakCur << std::endl;
    Xyce::dout() << "tBCLeakCur = " << tBELeakCur << std::endl;
    Xyce::dout() << "tBEPot     = " << tBEPot << std::endl;
    Xyce::dout() << "tBECap     = " << tBECap << std::endl;
    Xyce::dout() << "tBCPot     = " << tBCPot << std::endl;
    Xyce::dout() << "tBCCap     = " << tBCCap << std::endl;
    Xyce::dout() << "tDepCap    = " << tDepCap << std::endl;
    Xyce::dout() << "tF1        = " << tF1 << std::endl;
    Xyce::dout() << "tF4        = " << tF4 << std::endl;
    Xyce::dout() << "tF5        = " << tF5 << std::endl;
    Xyce::dout() << "tleakBEEmissionCoeff = " << tleakBEEmissionCoeff << std::endl;
    Xyce::dout() << "tleakBCEmissionCoeff = " << tleakBCEmissionCoeff << std::endl;
    Xyce::dout() << "tRollOffExp      = " << tRollOffExp << std::endl;
    Xyce::dout() << "tInvRollOffF     = " << tInvRollOffF << std::endl;
    Xyce::dout() << "tInvRollOffR     = " << tInvRollOffR << std::endl;
    Xyce::dout() << "tInvEarlyVoltF   = " << tInvEarlyVoltF << std::endl;
    Xyce::dout() << "tInvEarlyVoltR   = " << tInvEarlyVoltR << std::endl;
    Xyce::dout() << "tBaseResist      = " << tBaseResist << std::endl;
    Xyce::dout() << "tCollectorResist = " << tCollectorResist << std::endl;
    Xyce::dout() << "tEmitterResist   = " << tEmitterResist << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << std::endl;
  }
#endif

  return true;
}


//----------------------------------------------------------------------------
// Function      : Instance::lambertWCurrent
// Purpose       : Determine the current through an individual BJT
// Special Notes :
// Scope         : public
// Creator       : Nick Johnson, Summer Intern
// Creation Date : 7/30/02
//---------------------------------------------------------------------------
bool Instance::lambertWCurrent(double &Id, double &Gd, double Vd, double Vte, double Isat)
{
  double RS = 1.0e-12;
  double arg1 = (Vd + Isat*RS)/Vte;
  arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
  double evd = exp(arg1);
  double lambWArg = Isat*RS*evd/Vte;
  double lambWReturn;
  int ierr;
  double lambWError;
  devSupport.lambertw(lambWArg, lambWReturn, ierr, lambWError);

  Id = -Isat+Vte*(lambWReturn)/RS;
  Gd = lambWReturn / ((1 + lambWReturn)*RS);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDeviceMask
//
// Purpose       : Loads the zero elements of the device mask
//
// Special Notes : elements of the error vector associated with zero
//                 elements of the mask will not be included in weighted
//                 norms by the time integrator.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 01/19/07
//-----------------------------------------------------------------------------
bool Instance::loadDeviceMask ()
{
  bool returnVal=false;
  if (getDeviceOptions().newExcessPhase)
  {
    N_LAS_Vector * maskVectorPtr = extData.deviceMaskVectorPtr;

    (*maskVectorPtr)[li_Ifx] = 0.0;
    (*maskVectorPtr)[li_dIfx] = 0.0;
    returnVal = true;
  }
  return (returnVal);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 BJT instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double td = model_.excessPhaseFac;
  double vbe_diff = vBE - vBE_orig;
  double vbc_diff = vBC - vBC_orig;

  double * solVec = extData.nextSolVectorRawPtr;
  double * qVec = extData.daeQVectorRawPtr;

  qVec[li_Base] -= - model_.TYPE * qBX;
  qVec[li_Subst] -= -model_.TYPE * qCS;
  qVec[li_CollP] -= model_.TYPE * ( qCS + qBX + qBCdep + qBCdiff );
  qVec[li_BaseP] -= -model_.TYPE * ( qBEdep + qBEdiff + qBCdep + qBCdiff );
  qVec[li_EmitP] -= model_.TYPE*( qBEdep + qBEdiff );

  // excess phase ERK-dcop
  if (td != 0 && getDeviceOptions().newExcessPhase)
  {
    qVec[li_Ifx] += solVec[li_Ifx];

    if (!(getSolverState().dcopFlag) )
    {
      qVec[li_dIfx] += solVec[li_dIfx]*td*td;
    }
    else
    {
      qVec[li_dIfx] = 0.0;
    }
  }

  // load up the jdxp terms
  if (getDeviceOptions().voltageLimiterFlag)
  {
    double Cp_Jdxp_q = 0.0;
    double Ep_Jdxp_q = 0.0;
    double Bp_Jdxp_q = 0.0;

    if (!origFlag)
    {
      Cp_Jdxp_q =  -(capBCdep + capBCdiff)*vbc_diff;
      Cp_Jdxp_q *= model_.TYPE;

      Bp_Jdxp_q =  (capBEdep + capBEdiff)*vbe_diff
        + (capBCdiff + capBCdep + capeqCB) *vbc_diff;
      Bp_Jdxp_q *= model_.TYPE;

      Ep_Jdxp_q = - capeqCB * vbc_diff - (capBEdiff + capBEdep)* vbe_diff;
      Ep_Jdxp_q *= model_.TYPE;
    }

    double * dQdxdVp = extData.dQdxdVpVectorRawPtr;
    dQdxdVp[li_CollP] += Cp_Jdxp_q;
    dQdxdVp[li_BaseP] += Bp_Jdxp_q;
    dQdxdVp[li_EmitP] += Ep_Jdxp_q;
  }

  if( loadLeadCurrent )
  {
    double * storeLeadQ = extData.storeLeadCurrQCompRawPtr;

    storeLeadQ[li_store_dev_ic] = -model_.TYPE * ( qCS + qBX + qBCdep + qBCdiff );
    storeLeadQ[li_store_dev_ib] = model_.TYPE * ( qBX + qBEdep + qBEdiff + qBCdep + qBCdiff );
    storeLeadQ[li_store_dev_ie] = -model_.TYPE*( qBEdep + qBEdiff );
    storeLeadQ[li_store_dev_is] = model_.TYPE * qCS;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::auxDAECalculations
//
// Purpose       : This function is supposed to handle a lot of the final
//                 "cleanup" calculations that were done, in the old-DAE
//                 formulation, at the end of updateSecondaryState.
//
//                 For the new-DAE form, they need to be done elsewhere, as
//                 updateSecondaryState isn't called for new-DAE.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences.
// Creation Date : 04/15/04
//-----------------------------------------------------------------------------
void Instance::auxDAECalculations ()
{
  double i_fx;
  double td = model_.excessPhaseFac;
  double * solVec = extData.nextSolVectorRawPtr;

  if ( td != 0 && !(getSolverState().dcopFlag) )
  {
    i_fx = solVec[li_Ifx];
    iCE = i_fx - iBC / qB;
  }
  else
  {
    iCE = (iBE - iBC)/ qB;
  }

  iC = iCE - iBC / tBetaR - iBCleak;
  iB  = iBE / tBetaF + iBEleak + iBC / tBetaR + iBCleak;
  iE = -iC-iB;

  if (td != 0)
  {
    if (!(getSolverState().dcopFlag) )
    {
      diCEdvBp = invqB * (-invqB * iBC  * dqBdvBp - gBC);
      diCEdvCp = invqB * (-invqB * iBC * dqBdvCp + gBC);
      diCEdvEp = -invqB * invqB * iBC * dqBdvEp;
    }
    else // ERK-dcop:  if dcop, use the td=0 case (copied from below)
    {
      diCEdvBp = invqB * (iCE * dqBdvBp + gBE - gBC);
      diCEdvCp = invqB * (iCE * dqBdvCp + gBC);
      diCEdvEp = invqB * (iCE * dqBdvEp - gBE);
    }

    diBEdvBp = invqB *(invqB* iBE * dqBdvBp + gBE);
    diBEdvCp = invqB * invqB * iBE * dqBdvCp;
    diBEdvEp = invqB * (invqB * iBE * dqBdvEp - gBE);
  }
  else
  {
    diCEdvBp = invqB * (iCE * dqBdvBp + gBE - gBC);
    diCEdvCp = invqB * (iCE * dqBdvCp + gBC);
    diCEdvEp = invqB * (iCE * dqBdvEp - gBE);
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 BJT instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double coefC, coefB, coefE, coefCp, coefBp, coefEp;

  double td = model_.excessPhaseFac;
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;

  double vbe_diff = vBE - vBE_orig;
  double vbc_diff = vBC - vBC_orig;
  double vce_diff = vbe_diff - vbc_diff;

  fVec[li_Coll] -= -vCCp * gCpr;
  fVec[li_Base] -= -vBBp * gX;
  fVec[li_Emit] -= -vEEp * gEpr;
  fVec[li_CollP] -= vCCp * gCpr + model_.TYPE * ( - iC );
  fVec[li_BaseP] -= vBBp * gX - model_.TYPE * ( iB );
  fVec[li_EmitP] -= vEEp * gEpr + model_.TYPE * ( - iE );

  // excess phase ERK-dcop.
  double i_fx = 0.0;
  double di_fx = 0.0;

  if (getDeviceOptions().newExcessPhase)
  {
    i_fx = solVec[li_Ifx];
    di_fx = solVec[li_dIfx];

    if (td != 0)
    {
      if (!(getSolverState().dcopFlag) )
      {
        // omega0 = 1/td;
        fVec[li_Ifx] += - di_fx;
        fVec[li_dIfx] += 3 * di_fx*td + 3*i_fx -3 * iBE / qB;
      }
      else
      {
        fVec[li_Ifx] += i_fx -iBE/qB;
        fVec[li_dIfx] = 0.0;
      }
    }
    else
    {
      fVec[li_Ifx] += i_fx;
      fVec[li_dIfx] += di_fx;
    }
  }

  double Cp_Jdxp_f=0.0;
  double Bp_Jdxp_f=0.0;
  double Ep_Jdxp_f=0.0;
  double dIfx_Jdxp_f=0.0 ;
  double Ifx_Jdxp_f=0.0 ;

  // Now for the jdxp terms
  if (getDeviceOptions().voltageLimiterFlag)
  {
    if (!origFlag)
    {
      Cp_Jdxp_f = + diCEdvBp * vbe_diff
                  + diCEdvCp * vce_diff
                  - gBCtot  * vbc_diff;

      Cp_Jdxp_f *= model_.TYPE;

      Bp_Jdxp_f =  gBEtot * vbe_diff + gBCtot  * vbc_diff;
      Bp_Jdxp_f *= model_.TYPE;

      Ep_Jdxp_f = - diCEdvCp * vce_diff - (diCEdvBp + gBEtot) * vbe_diff;
      Ep_Jdxp_f *= model_.TYPE;

      // ERK-dcop.
      if ( td != 0 && getDeviceOptions().newExcessPhase )
      {
        if ( !(getSolverState().dcopFlag) )
        {
          dIfx_Jdxp_f =  -3 *( diBEdvBp*vbe_diff + diBEdvCp * vce_diff);
          dIfx_Jdxp_f *= model_.TYPE;
        }
        else
        {
          Ifx_Jdxp_f =  ( diBEdvBp*vbe_diff + diBEdvCp * vce_diff);
          Ifx_Jdxp_f *= model_.TYPE;
        }
      }
    }

    double * dFdxdVp = extData.dFdxdVpVectorRawPtr;
    dFdxdVp[li_CollP] += Cp_Jdxp_f;
    dFdxdVp[li_BaseP] += Bp_Jdxp_f;
    dFdxdVp[li_EmitP] += Ep_Jdxp_f;

    // ERK-dcop.
    if ( td != 0 && getDeviceOptions().newExcessPhase )
    {
      if ( !(getSolverState().dcopFlag) )
      {
        dFdxdVp[li_dIfx] += dIfx_Jdxp_f;
      }
      else
      {
        dFdxdVp[li_Ifx] += Ifx_Jdxp_f;
      }
    }
  }

  if( loadLeadCurrent )
  {
    double * storeLeadF = extData.nextStoVectorRawPtr;

    storeLeadF[li_store_dev_ic] = model_.TYPE * ( iC );
    storeLeadF[li_store_dev_is] = 0;
    storeLeadF[li_store_dev_ie] = model_.TYPE * ( iE );
    storeLeadF[li_store_dev_ib] = model_.TYPE * ( iB );
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 BJT instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
//                 The "Q" vector contains charges and fluxes.  So, this
//                 matrix will contain only the capacitance terms.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  double td = model_.excessPhaseFac;

  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  dQdx[li_Base][ABaseEquBaseNodeOffset] += capBX;

  dQdx[li_Base][ABaseEquCollPNodeOffset] += -capBX;
  dQdx[li_Subst][ASubstEquSubstNodeOffset] += capCS;
  dQdx[li_Subst][ASubstEquCollPNodeOffset] -= capCS;

  dQdx[li_CollP][ACollPEquBaseNodeOffset] -= capBX;
  dQdx[li_CollP][ACollPEquSubstNodeOffset] -= capCS;
  dQdx[li_CollP][ACollPEquCollPNodeOffset]
    +=  capCS + capBX + capBCdep + capBCdiff;

  dQdx[li_CollP][ACollPEquBasePNodeOffset]
              += -capBCdep - capBCdiff;

  dQdx[li_BaseP][ABasePEquCollPNodeOffset]
    += -capBCdiff - capBCdep - capeqCB;

  dQdx[li_BaseP][ABasePEquBasePNodeOffset]
      += capBEdiff + capBEdep + capBCdiff + capBCdep + capeqCB;

  dQdx[li_BaseP][ABasePEquEmitPNodeOffset]
                += -capBEdiff - capBEdep;
  dQdx[li_EmitP][AEmitPEquCollPNodeOffset] += capeqCB;
  dQdx[li_EmitP][AEmitPEquBasePNodeOffset]
      += -capBEdiff - capBEdep - capeqCB;
  dQdx[li_EmitP][AEmitPEquEmitPNodeOffset]
      += capBEdiff + capBEdep;

  // excess phase terms.  ERK-dcop
  if ( td != 0 && getDeviceOptions().newExcessPhase )
  {
    if (!(getSolverState().dcopFlag) )
    {
      dQdx[li_Ifx][AIfxEquIfxNodeOffset] += 1;
      dQdx[li_dIfx][AdIfxEqudIfxNodeOffset] += 1*td*td;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 BJT instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
//                 This is an algebraic constaint.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  double td = model_.excessPhaseFac;

  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Coll][ACollEquCollNodeOffset] += gCpr;
  dFdx[li_Coll][ACollEquCollPNodeOffset] -= gCpr;
  dFdx[li_Base][ABaseEquBaseNodeOffset] += diBrdvB;
  dFdx[li_Base][ABaseEquCollPNodeOffset] += diBrdvCp;
  dFdx[li_Base][ABaseEquBasePNodeOffset] += diBrdvBp;
  dFdx[li_Base][ABaseEquEmitPNodeOffset] += diBrdvEp;

  dFdx[li_Emit][AEmitEquEmitNodeOffset] += gEpr;

  dFdx[li_Emit][AEmitEquEmitPNodeOffset] -= gEpr;


  dFdx[li_CollP][ACollPEquCollNodeOffset] -= gCpr;
  dFdx[li_CollP][ACollPEquCollPNodeOffset] +=diCEdvCp + gBCtot + gCpr;
  dFdx[li_CollP][ACollPEquBasePNodeOffset]
                += diCEdvBp - gBCtot;
  dFdx[li_CollP][ACollPEquEmitPNodeOffset] += diCEdvEp;

  // excess phase terms.  ERK-dcop
  if ( td != 0 && getDeviceOptions().newExcessPhase && !(getSolverState().dcopFlag) )
    dFdx[li_CollP][ACollPEquIfxNodeOffset] += 1.0 * model_.TYPE;

  dFdx[li_BaseP][ABasePEquBaseNodeOffset] -= diBrdvB;
  dFdx[li_BaseP][ABasePEquCollPNodeOffset]
    += -diBrdvCp - gBCtot;

  dFdx[li_BaseP][ABasePEquBasePNodeOffset]
      += -diBrdvBp + gBEtot + gBCtot;

  dFdx[li_BaseP][ABasePEquEmitPNodeOffset]
                += -diBrdvEp - gBEtot;

  dFdx[li_EmitP][AEmitPEquEmitNodeOffset] -= gEpr;
  dFdx[li_EmitP][AEmitPEquCollPNodeOffset] += -diCEdvCp;
  dFdx[li_EmitP][AEmitPEquBasePNodeOffset]
      += - diCEdvBp - gBEtot;

  // ERK Note:  -diCEdvEp should equal + diCEdvBp + diCEdvCp.  In the
  // original old-DAE form, + diCEdvBp + diCEdvCp is used.  In
  // the more recent new-DAE form, we've used diCEdvEp, but it sometimes
  // results in subtle differences between old- and new-DAE jacobians.
  // This should matter, but .....
  dFdx[li_EmitP][AEmitPEquEmitPNodeOffset]
      //+= gBEtot + gEpr - diCEdvEp;
      += gBEtot + gEpr + diCEdvBp + diCEdvCp; // this is a test.

  // excess phase terms.  ERK-dcop
  if ( td != 0 && getDeviceOptions().newExcessPhase && !(getSolverState().dcopFlag) )
    dFdx[li_EmitP][AEmitPEquIfxNodeOffset] += -1.0 * model_.TYPE;

  // excess phase terms.  ERK-dcop
  if (getDeviceOptions().newExcessPhase)
  {
    if ( td != 0 )
    {
      if ( !(getSolverState().dcopFlag) )
      {
        dFdx[li_Ifx][AIfxEqudIfxNodeOffset]+= -1;
        dFdx[li_dIfx][AdIfxEquCollPNodeOffset] += -3*diBEdvCp * model_.TYPE;
        dFdx[li_dIfx][AdIfxEquBasePNodeOffset] += -3*diBEdvBp * model_.TYPE;
        dFdx[li_dIfx][AdIfxEquEmitPNodeOffset] += -3*diBEdvEp * model_.TYPE;
        dFdx[li_dIfx][AdIfxEqudIfxNodeOffset] += 3* td;
        dFdx[li_dIfx][AdIfxEquIfxNodeOffset] += 3;
      }
      else
      {
        dFdx[li_Ifx][AIfxEquCollPNodeOffset] += -diBEdvCp * model_.TYPE;
        dFdx[li_Ifx][AIfxEquBasePNodeOffset] += -diBEdvBp * model_.TYPE;
        dFdx[li_Ifx][AIfxEquEmitPNodeOffset] += -diBEdvEp * model_.TYPE;
        dFdx[li_Ifx][AIfxEquIfxNodeOffset]+= 1;
        dFdx[li_dIfx][AdIfxEqudIfxNodeOffset] += 1;
      }
    }
    else
    {
      dFdx[li_Ifx][AIfxEquIfxNodeOffset]+= 1;
      dFdx[li_dIfx][AdIfxEqudIfxNodeOffset] += 1;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra
// Creation Date : 1/13/00
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  // Do the bulk of the work in updateIntermediateVars:
  bool bsuccess = updateIntermediateVars();

  double * staVec = extData.nextStaVectorRawPtr;
  double * currStaVec = extData.currStaVectorRawPtr;

  // save voltage drops

  double * stoVec = extData.nextStoVectorRawPtr;
  stoVec[li_storevBE] = vBE;
  stoVec[li_storevBC] = vBC;
  stoVec[li_store_capeqCB] = capeqCB;

  staVec[li_qstateBEdiff] = qBEdiff;
  staVec[li_qstateBEdep] = qBEdep;
  staVec[li_qstateBCdiff] = qBCdiff;
  staVec[li_qstateBCdep] = qBCdep;
  staVec[li_qstateBX] = qBX;
  staVec[li_qstateCS] = qCS;


  // if this is the first newton step of the first time step
  // of the transient simulation, we need to enforce that the
  // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
  // compatibility.  ERK.

  if (!(getSolverState().dcopFlag) && (getSolverState().initTranFlag) && getSolverState().newtonIter==0)
  {
    currStaVec[li_qstateBEdiff] = qBEdiff;
    currStaVec[li_qstateBEdep] = qBEdep;
    currStaVec[li_qstateBCdiff] = qBCdiff;
    currStaVec[li_qstateBCdep] = qBCdep;
    currStaVec[li_qstateBX] = qBX;
    currStaVec[li_qstateCS] = qCS;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra
// Creation Date : 1/13/00
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  double * staDerivVec = extData.nextStaDerivVectorRawPtr;

  // Now that the state vector for time=0 is up-to-date, get the derivative
  // with respect to time of the charge, to obtain the best estimate for
  // the current in the capacitor.

  iBEdiff = staDerivVec[li_qstateBEdiff];
  iBEdep  = staDerivVec[li_qstateBEdep ];
  iCS     = staDerivVec[li_qstateCS    ];
  iBCdiff = staDerivVec[li_qstateBCdiff];
  iBCdep  = staDerivVec[li_qstateBCdep ];
  iBX     = staDerivVec[li_qstateBX    ];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::oldDAEExcessPhaseCalculation1
//
// Purpose       : This handles the parts of the calculation of excess phase
//                 terms that need to be done in the
//                 updatePrimaryState function.
//
// Special Notes : The excess phase calculation is based on a second-order
//                 time derivative.
//
//                 Most of the excess phase work is done in the companion
//                 function, oldDAEExcessPhaseCalculation2.
//
//                 This function mainly enforces a constant history.
//
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences.
// Creation Date : 04/11/04
//-----------------------------------------------------------------------------
void Instance::oldDAEExcessPhaseCalculation1 ()
{
  // do excess phase stuff if we are in transient mode and td is nonzero.
  double td = model_.excessPhaseFac;

  currCexbc = lastCexbc = 0.0;

  if (!(getSolverState().dcopFlag) && td != 0)
  {
    // If there is no cexbc history, create one and use it.
    if (getSolverState().beginIntegrationFlag)
    {
      currCexbc = lastCexbc = iBE/qB;
      (*extData.currStaVectorPtr)[li_istateCEXBC] = currCexbc;
      (*extData.lastStaVectorPtr)[li_istateCEXBC] = lastCexbc;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::oldDAEExcessPhaseCalculation2
//
// Purpose       : This handles the parts of the calculation of excess phase
//                 terms that need to be done in the
//                 updateSecondaryState function.
//
// Special Notes : This function has a companion function,
//                 oldDAEExcessPhaseCalculation1.  Most of the work is done
//                 in this one.
//
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences.
// Creation Date : 04/11/04
//-----------------------------------------------------------------------------
void Instance::oldDAEExcessPhaseCalculation2
   (double & iEX, double & gEX, double & iC_local)
{
  ////////////////////////////////////////////////////////////////
  double td = model_.excessPhaseFac;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "  Excess Phase stuff:" <<std::endl;
    Xyce::dout() << "  name = " << getName() <<std::endl;
    Xyce::dout() << "  td   = " << td <<std::endl;
  }
#endif

  ////////////////////////////////////////////////////////////////
  //  ERK: 07/13/01
  //
  //   weil's approx. for excess phase applied with backward-
  //   euler integration
  //  (copied over from chilespice).
  //
  ////////////////////////////////////////////////////////////////

  iEX=iBE;
  gEX=gBE;

  dt0 = getSolverState().currTimeStep;
  dt1 = getSolverState().lastTimeStep;

  iC_local = 0.0;

  // do excess phase stuff if we are in transient mode and td is nonzero.
  if (!(getSolverState().dcopFlag) && td != 0)
  {
    double arg1, arg2, denom;

    arg1  = dt0 / td;
    arg2  = 3.0 * arg1;
    arg1  = arg2 * arg1;
    denom = 1.0 + arg1 + arg2;
    phaseScalar  = arg1 / denom;

    // Secondary state stuff:
    if (!getSolverState().beginIntegrationFlag)
    {
      currCexbc = (*extData.currStaVectorPtr)[li_istateCEXBC];
      lastCexbc = (*extData.lastStaVectorPtr)[li_istateCEXBC];
    }

    iC_local = ((currCexbc) * (1 + dt0 / dt1 + arg2) - (lastCexbc) * dt0 / dt1) / denom;
    iEX = iBE * phaseScalar;
    gEX = gBE * phaseScalar;

    nextCexbc = iC_local + iEX / qB;
    (*extData.nextStaVectorPtr)[li_istateCEXBC] = nextCexbc;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra
// Creation Date : 1/13/00
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{

  bool bsuccess = true;

  double v_emit, v_emitP;
  double v_base, v_baseP;
  double v_coll, v_collP;
  double v_subst;

  double q1, q2;
  double * solVec = extData.nextSolVectorRawPtr;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel>0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "Instance::updateIntermediateVars " << getName() <<std::endl;
  }
#endif

  // obtain voltages:
  v_emit = solVec[li_Emit];
  v_emitP = solVec[li_EmitP];
  v_base = solVec[li_Base];
  v_baseP = solVec[li_BaseP];
  v_coll = solVec[li_Coll];
  v_collP = solVec[li_CollP];
  v_subst = solVec[li_Subst];

  vEEp = (v_emit - v_emitP);
  vBBp = (v_base - v_baseP);
  vCCp = (v_coll - v_collP);

  vBE = model_.TYPE * (v_baseP - v_emitP);
  vBC = model_.TYPE * (v_baseP - v_collP);
  vBX = model_.TYPE * (v_base - v_collP);
  //vCS = model_.TYPE * (v_collP - v_subst);
  vCS = model_.TYPE * (v_subst - v_collP);

  vBE_orig = vBE;
  vBC_orig = vBC;

  // Reset this to show we're starting from the originals
  origFlag = true;
  offFlag = false;

  if (getSolverState().initJctFlag && !OFF && getDeviceOptions().voltageLimiterFlag)
  {
    if( IC_GIVEN )
    {
      vBE = model_.TYPE * icVBE;
      double vCE = model_.TYPE * icVCE;
      vBC = vBE - vCE;
      vBX = vBC;
      origFlag = false;
    }
    else
    {
      if (getSolverState().inputOPFlag)
      {
        N_LAS_Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
        if ((*flagSolVectorPtr)[li_Emit] == 0 || (*flagSolVectorPtr)[li_EmitP] == 0 ||
            (*flagSolVectorPtr)[li_Base] == 0 || (*flagSolVectorPtr)[li_BaseP] == 0 ||
            (*flagSolVectorPtr)[li_Coll] == 0 || (*flagSolVectorPtr)[li_CollP] == 0 ||
            (*flagSolVectorPtr)[li_Subst] == 0 )
        {
          vBC = 0;
          vBE = tVCrit;
          vBX = vBC;
          origFlag = false;
        }
      }
      else
      {
        vBE=tVCrit;
        vBX = vBC;
        origFlag = false;
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << " Setting device initial condition to Base-Emitter drop=tVCrit (" << tVCrit << ")"<<std::endl;
        }
#endif
      }
    }


#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << "  " <<std::endl;
      Xyce::dout() << "  using UIC.\n";
      Xyce::dout() << "  vBE = " << vBE << std::endl;
      Xyce::dout() << "  vBC = " << vBC << std::endl;
      Xyce::dout() << "  vBX = " << vBX << std::endl;
      Xyce::dout() << "  vCS = " << vCS << std::endl;
    }
#endif
  }
  else if ((getSolverState().initFixFlag || getSolverState().initJctFlag) && OFF)
  {
    vBX = vBC = vBE = 0.0 ;
    // NOTE:  Must not set "origFlag" here, because that flags "non-converged"
    // which causes initFixFlag always to be set, which causes this section
    // always to be executed.  Bad!
    // instead, flag with "offFlag" that will not be used in convergence
    // testing, but WILL be used in the voltage limiting machinery to make
    // this work properly.
    offFlag = true;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << "  " <<std::endl;
      Xyce::dout() << "  BJT explicitly set OFF, zeroing all junction drops.\n";
    }
#endif
  }


  if (getSolverState().newtonIter == 0)
  {
    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      vBE_old = (*extData.currStoVectorPtr)[li_storevBE];
      vBC_old = (*extData.currStoVectorPtr)[li_storevBC];
      capeqCB = (*extData.currStoVectorPtr)[li_store_capeqCB];
    }
    else
    {
      vBE_old = vBE;           // otherwise we have no history
      vBC_old = vBC;
      capeqCB = 0.0;
    }
  }
  else
  {
    vBE_old = (*extData.nextStoVectorPtr)[li_storevBE];
    vBC_old = (*extData.nextStoVectorPtr)[li_storevBC];
    capeqCB = (*extData.nextStoVectorPtr)[li_store_capeqCB];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -2 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "  tVCrit = " << tVCrit << std::endl;
    Xyce::dout().width(3);
    Xyce::dout() << getSolverState().newtonIter;
    Xyce::dout().width(5); Xyce::dout() << getName();
    Xyce::dout() << "  Blim: ";
    Xyce::dout() << "  vBE=";
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << vBE;
    Xyce::dout() << "  vBC=";
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << vBC << std::endl;
  }
#endif

  // if we had bypass implemented, most of it would be here.
  // don't bother limiting voltages if we've already forced OFF.

  if (getDeviceOptions().voltageLimiterFlag && !(getSolverState().initFixFlag&&OFF))
  {
    // "reset" the junction voltages to keep them from changing too much.
    // only do this if this is not the first  newton iteration, and not
    // a by pass step.
    int ichk1, icheck;

    if (getSolverState().newtonIter >= 0)
    {
      ichk1=1;  // don't know what this is for...

      vBE = devSupport.pnjlim(vBE, vBE_old, vt, tVCrit, &icheck);
      vBC = devSupport.pnjlim(vBC, vBC_old, vt, tVCrit, &ichk1);

      if (ichk1 == 1) icheck=1;

      // set the origFlag to zero if things have changed.
      if (icheck==1) origFlag = false;
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > -2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout().width(3);
      Xyce::dout() << getSolverState().newtonIter;
      Xyce::dout().width(5); Xyce::dout() << getName() ;
      Xyce::dout() << "  Alim: ";
      Xyce::dout() << "  vBE=";
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vBE;
      Xyce::dout() << "  vBC=";
      Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << vBC;

      //cout << "  vBE_orig = " << vBE_orig << std::endl;
      //cout << "  vBC_orig = " << vBC_orig << std::endl <<std::endl;
      //cout << "  origFlag = " << origFlag << std::endl;

      if (origFlag) Xyce::dout() << " SAME";
      else          Xyce::dout() << " DIFF";
      Xyce::dout() << std::endl;
    }
#endif
  }

  //Junction Current Calculations

  double csat = tSatCur * AREA;
  double vtF = vt * model_.emissionCoeffF;
  double vtR = vt * model_.emissionCoeffR;
  double vtE = vt * tleakBEEmissionCoeff;
  double vtC = vt * tleakBCEmissionCoeff;
  double iLeakBE = tBELeakCur * AREA;
  double iLeakBC = tBCLeakCur * AREA;

  // determine currents using LambertW
  if (lambertWFlag)
  {
    if( vBE > -5.0 * vtF )
    {
      lambertWCurrent(iBE, gBE, vBE, vtF, csat);

      if( iLeakBE == 0.0 )
        iBEleak = gBEleak = 0.0;
      else
      {
        lambertWCurrent(iBEleak, gBEleak, vBE, vtE, iLeakBE);
      }
    }
    else
    {
      gBE = -csat / vBE + getDeviceOptions().gmin;
      iBE = gBE * vBE;
      gBEleak = -iLeakBE / vBE;
      iBEleak = gBEleak * vBE;
    }

    if( vBC > -5.0 * vtR )
    {
      lambertWCurrent(iBC, gBC, vBC, vtR, csat);

      if( iLeakBC == 0.0 )
        iBCleak = gBCleak = 0.0;
      else
      {
        lambertWCurrent(iBCleak, gBCleak, vBC, vtC, iLeakBC);
      }
    }
    else
    {
      gBC = -csat / vBC + getDeviceOptions().gmin;
      iBC = gBC * vBC;
      gBCleak = -iLeakBC / vBC;
      iBCleak = gBCleak * vBC;
    }
  }

  // determine currents using normal exponential
  else
  {
    if (vBE > -5.0 * vtF)
    {
      double arg1 = vBE / vtF ;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double devBE;
      double evBE;
      if (getSolverState().bjtArtParameterFlag)
      {
        evBE = devSupport.Xexp( arg1 , devBE, getDeviceOptions().exp_order );
      }
      else
      {
         evBE = exp( arg1 );
         devBE=evBE;
      }
      iBE = csat * ( evBE - 1.0 ) + getDeviceOptions().gmin * vBE;
      gBE = csat * devBE / vtF + getDeviceOptions().gmin;

      if (iLeakBE == 0.0)
        iBEleak = gBEleak = 0.0;
      else
      {
        double arg1 = vBE / vtE ;
        arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
        double devBEleak;
        double evBEleak;
        if (getSolverState().bjtArtParameterFlag)
        {
          evBEleak =devSupport.Xexp(arg1,devBEleak,getDeviceOptions().exp_order );
        }
        else
        {
          evBEleak = exp(arg1);
          devBEleak=evBEleak;
        }
        iBEleak = iLeakBE * (evBEleak - 1.0);
        gBEleak = iLeakBE * devBEleak / vtE;
      }
    }
    else
    {
      gBE = -csat / vBE + getDeviceOptions().gmin;
      iBE = gBE * vBE;
      gBEleak = -iLeakBE / vBE;
      iBEleak = gBEleak * vBE;
    }

    if( vBC > -5.0 * vtR )
    {
      double arg1 = vBC / vtR ;
      arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
      double devBC;
      double evBC;
      if (getSolverState().bjtArtParameterFlag)
      {
        evBC = devSupport.Xexp( arg1 , devBC, getDeviceOptions().exp_order );
      }
      else
      {
        evBC = exp( arg1 );
        devBC=evBC;
      }
      iBC = csat * ( evBC - 1.0 ) + getDeviceOptions().gmin * vBC;
      gBC = csat * devBC / vtR + getDeviceOptions().gmin;

      if (iLeakBC == 0.0)
        iBCleak = gBCleak = 0.0;
      else
      {
        double arg1 = vBC / vtC ;
        arg1 = Xycemin(CONSTMAX_EXP_ARG, arg1);
        double devBCleak;
        double evBCleak;
        if (getSolverState().bjtArtParameterFlag)
        {
          evBCleak= devSupport.Xexp(arg1,devBCleak,getDeviceOptions().exp_order);
        }
        else
        {
          evBCleak = exp(arg1);
          devBCleak=evBCleak;
        }
        iBCleak = iLeakBC * (evBCleak - 1.0);
        gBCleak = iLeakBC * devBCleak / vtC;
      }
    }
    else
    {
      gBC = -csat / vBC + getDeviceOptions().gmin;
      iBC = gBC * vBC;
      gBCleak = -iLeakBC / vBC;
      iBCleak = gBCleak * vBC;
    }
  }

  // base charge calculations:
  double rofF = tInvRollOffF / AREA;
  double rofR = tInvRollOffR / AREA;

  q1 = 1.0 / ( 1.0 - tInvEarlyVoltF * vBC - tInvEarlyVoltR * vBE );

  if (rofF == 0.0 && rofR == 0.0)
  {
    qB = q1;
    double q1_qB = q1 * qB;

    dqBdvEp = q1_qB * tInvEarlyVoltR;
    dqBdvCp = q1_qB * tInvEarlyVoltF;
/* replaced with simpler expression after if block
    dqBdvBp = -q1_qB * (model_.invEarlyVoltR + model_.invEarlyVoltF);
*/
  }
  else
  {
    q2 = rofF * iBE  + rofR * iBC;

    // Difference from SPICE: spice puts a max(0,) around the arg definition,
    // which means we use the square root expression all the way down to
    // q2= -1/4.  This is not really good, because it means that the derivative
    // of the square root term goes to negative infinity at the break point
    // between the expressions, and qB will be very discontinuous.
    //  This doesn't have a huge effect on the solutions, but can possibly
    // lead to convergence issues (see bug 1214).  This fix hasn't yet
    // been shown to help anything, but it has not hurt anything either.
    if (q2 >= 0)
    {

      double arg = 1.0 + 4.0 * q2;
      double sqarg = 1.0;

      //    if (fabs(arg) > 0.0) sqarg = sqrt(arg);
      //    double rofF_gBE_invSqarg = rofF * gBE / sqarg;
      //    double rofR_gBC_invSqarg = rofR * gBC / sqarg;
      //  add Pspice compatible rolloff parameter NK = tRollOffExp

      if (fabs(arg) > 0.0) sqarg = pow(arg,tRollOffExp);
      double rofF_gBE_invSqarg = 0.0;
      double rofR_gBC_invSqarg = 0.0;
      if(arg != 0) rofF_gBE_invSqarg = rofF*gBE*2*tRollOffExp*sqarg/arg;
      if(arg != 0) rofR_gBC_invSqarg = rofR*gBC*2*tRollOffExp*sqarg/arg;

      qB = 0.5 * q1 * (1.0 + sqarg);

      dqBdvEp = q1 * (qB * tInvEarlyVoltR + rofF_gBE_invSqarg);
      dqBdvCp = q1 * (qB * tInvEarlyVoltF + rofR_gBC_invSqarg);

      /* replaced with simpler expression after if block
         dqBdvBp = -q1 * (qB * ( model_.invEarlyVoltF + model_.invEarlyVoltR) +
         rofF_gBE_invSqarg + rofR_gBC_invSqarg);
      */
    }
    else
    {
      qB = q1;
      dqBdvEp = q1*qB*tInvEarlyVoltR;
      dqBdvCp = q1*qB*tInvEarlyVoltF;
    }
  }
  dqBdvBp = - (dqBdvEp + dqBdvCp);

  invqB = 1.0 / qB;


  // NOTE: this is a block that (for some reason) was originally in the
  // updatePrimaryState function. It got moved as in general, the "state"
  // functions are only for manipulating the state vector.
  // NOTE: brackets are around this block of code to preserve variable scope.

  double arg, sarg;
  double ctot;
  double czBC;
  double czBX;
  double czBE;
  double czCS;
  double fcpc;
  double fcpe;

  // capacitance calculations and "current" charge value calculations:
  ctot = tBCCap * AREA;
  czBC = ctot * model_.baseFracBCCap;
  czBX = ctot - czBC;
  czBE = tBECap * AREA;
  czCS = model_.CJS * AREA;
  fcpc = model_.depCapCoeff * model_.potBC;
  fcpe = tDepCap;

  // High current perturbation to BE current
  // The high current perturbation to the BE current is only
  // used in determining capacitance values.  It does not enter
  // directly into the RHS KCL equations, or the Jacobian.
  iBEhighCurr = iBE;
  gBEhighCurr = gBE;

  // geqCB is only calculated sometimes (based on the if-statement).
  // If it is not calculated, then the value used will be the one from
  // the state vector, that was pulled out already.
  bool geqCB_recalc(false);
  if (model_.transTimeF != 0.0 && vBE > 0.0 &&
      ((!getSolverState().dcopFlag)||getSolverState().tranopFlag||getSolverState().acopFlag))
  {
    double argtf = 0.0;
    double arg2 = 0.0;
    double arg3 = 0.0;

    if ( model_.transTimeBiasCoeffF != 0.0)
    {
      argtf = model_.transTimeBiasCoeffF;

      if (model_.transTimeVBCFac != 0.0)
        argtf *= exp(vBC * model_.transTimeVBCFac);

      arg2 = argtf;

      if (model_.transTimeHighCurrF != 0.0)
      {
        double tmp = iBEhighCurr /
          (iBEhighCurr + model_.transTimeHighCurrF * AREA);

        argtf *= tmp * tmp;
        arg2   = argtf * (3.0 - 2.0 * tmp);
      }
      arg3 = iBEhighCurr * argtf * model_.transTimeVBCFac;
    }

    iBEhighCurr *= (1.0 + argtf) / qB;
    gBEhighCurr  = (gBEhighCurr * (1.0 + arg2) - iBEhighCurr * dqBdvEp) / qB;

    capeqCB = model_.transTimeF * (arg3 - iBEhighCurr * dqBdvCp) / qB;
    geqCB = capeqCB*getSolverState().pdt;
    geqCB_recalc = true;
  }

  // baseP-emitP depletion capacitance
  if( czBE == 0.0 )
  {
    qBEdep = 0.0;
    capBEdep = 0.0;
  }
  else
  {
    if ( vBE < fcpe )
    {
      arg = 1.0 - vBE / tBEPot;
      sarg = exp( -model_.juncExpBE * log( arg ) );

      qBEdep = tBEPot * czBE * ( 1.0 - arg * sarg ) /
         ( 1.0 - model_.juncExpBE );
      capBEdep = czBE * sarg;
    }
    else
    {
      double czBEf2 = czBE / model_.f2;

      qBEdep = czBE * tF1 + czBEf2 * ( model_.f3 * ( vBE -
                 fcpe ) + ( model_.juncExpBE / ( 2.0 *
                 tBEPot ) ) * ( vBE * vBE - fcpe * fcpe ) );
      capBEdep = czBEf2 * ( model_.f3 + model_.juncExpBE * vBE /
                   tBEPot );
    }
  }

  // baseP-emitP diffusion capacitance
  if (model_.transTimeF == 0.0)
  {
    qBEdiff = capBEdiff = 0.0;
  }
  else
  {
    qBEdiff   = model_.transTimeF * iBEhighCurr;
    capBEdiff = model_.transTimeF * gBEhighCurr;
  }

  // baseP-collP depletion capacitance
  if( czBC == 0.0 )
  {
    qBCdep = 0.0;
    capBCdep = 0.0;
  }
  else
  {
    if ( vBC < fcpc )
    {
      arg = 1.0 - vBC / tBCPot;
      sarg = exp( -model_.juncExpBC * log( arg ) );

      qBCdep = tBCPot * czBC * ( 1.0 - arg * sarg ) /
                 ( 1.0 - model_.juncExpBC );
      capBCdep = czBC * sarg;
    }
    else
    {
      double czBCf2 = czBC / model_.f6;

      qBCdep = czBC * tF5 + czBCf2 * ( model_.f7 * ( vBC -
                 fcpc ) + ( model_.juncExpBC / ( 2.0 *
                 tBCPot ) ) * ( vBC * vBC - fcpc * fcpc ) );
      capBCdep = czBCf2 * ( model_.f7 + model_.juncExpBC * vBC / tBCPot );
    }
  }

  // baseP-collP diffusion capacitance
  if (model_.transTimeR == 0.0)
  {
    qBCdiff = capBCdiff = 0.0;
  }
  else
  {
    qBCdiff   = model_.transTimeR * iBC;
    capBCdiff = model_.transTimeR * gBC;
  }

  // base-collP depletion capacitance
  if( czBX == 0.0 )
  {
    qBX = 0.0;
    capBX = 0.0;
  }
  else
  {
    if (vBX < fcpc)
    {
      arg = 1.0 - vBX / tBCPot;
      sarg = exp( -model_.juncExpBC * log( arg ) );

      qBX = tBCPot * czBX * ( 1.0 - arg * sarg ) /
                 ( 1.0 - model_.juncExpBC );
      capBX = czBX * sarg;
    }
    else
    {
      double czBXf2 = czBX / model_.f6;

      qBX = czBX * tF5 + czBXf2 * ( model_.f7 * ( vBX -
                 fcpc ) + ( model_.juncExpBC / ( 2.0 *
                 tBCPot ) ) * ( vBX * vBX - fcpc * fcpc ) );
      capBX = czBXf2 * ( model_.f7 + model_.juncExpBC * vBX /
                   tBCPot );
    }
  }

  // collP-subst depletion capacitance
  if( czCS == 0.0 )
  {
    qCS = 0.0;
    capCS = 0.0;
  }
  else
  {
    if ( vCS < 0.0 )
    {
      arg = 1.0 - vCS / model_.potSubst;
      sarg = exp( -model_.expSubst * log( arg ) );

      qCS = model_.potSubst * czCS * ( 1.0 - arg * sarg ) /
                 ( 1.0 - model_.expSubst );
      capCS = czCS * sarg;
    }
    else
    {
      qCS = vCS * czCS * ( 1.0 + model_.expSubst * vCS / ( 2.0 *
              model_.potSubst ) );
      capCS = czCS * ( 1.0 + model_.expSubst * vCS / model_.potSubst );
    }
  }

  // Here things get a little confusing.  Generally, this is a logical place
  // to calculate iB, iC, and iE.  iC and iE include the excess phase current.

  // Doing the old excess phase calculation, for both new and old-DAE, out
  // of curiousity to see that they match.  However, if running new-DAE, these
  // calculations will only be used in diagnostics, not the actual loads.
  double iEX_tmp = 0.0; double gEX_tmp = 0.0; double iC_tmp = 0.0;
  oldDAEExcessPhaseCalculation1 ();
  oldDAEExcessPhaseCalculation2 (iEX_tmp,gEX_tmp,iC_tmp);

  if (getDeviceOptions().newExcessPhase)
  {
    auxDAECalculations ();  // iB, iC and iE calculated here.

  }
  else
  {
    if (!(getSolverState().dcopFlag))
    {
      double pdt = getSolverState().pdt;

      // ERK. 12/17/2006.
      // It is necessary to set pdt=0.0, because for the first time step out of
      // the DCOP, all the time derivatives are forced to be zero.  Thus, all
      // their derivatives should also be zero.  If it wasn't for that, then pdt
      // could always be pdt.  (it used to be, before the -jacobian_test capability).
      if (getSolverState().initTranFlag && getSolverState().newtonIter==0)
      {
        pdt = 0.0;
      }

      gCapBEdiff = capBEdiff * pdt;
      gCapBEdep  = capBEdep  * pdt;
      gCapBCdiff = capBCdiff * pdt;
      gCapBCdep  = capBCdep  * pdt;
      gCapBX     = capBX     * pdt;
      gCapCS     = capCS     * pdt;
    }
    else
    {
      gCapBEdiff = 0.0;
      gCapBEdep  = 0.0;
      gCapBCdiff = 0.0;
      gCapBCdep  = 0.0;
      gCapBX     = 0.0;
      gCapCS     = 0.0;
    }

    // Do the iB, iC, and iE calculations for the old-DAE case here.
    double iEX = 0.0;
    double gEX = 0.0;

    iC = iC_tmp;
    iEX = iEX_tmp;
    gEX = gEX_tmp;

    iCE = ( iEX - iBC ) / qB;
    iC += iCE - iBC / tBetaR - iBCleak;

    iB  = iBE / tBetaF + iBEleak + iBC / tBetaR + iBCleak;
    iE = -iC-iB;

    // These 3 derivatives depend upon excess phase terms, old-DAE version, so
    // they have to be set up here.
    diCEdvEp = invqB * (iCE * dqBdvEp - gEX);
    diCEdvCp = invqB * (iCE * dqBdvCp + gBC);
    diCEdvBp = invqB * (iCE * dqBdvBp + gEX - gBC);

  }

  // Some derivatives setup:
  // Emitter and Conductor conductances
  gEpr = model_.emitterConduct * AREA;
  gCpr = model_.collectorConduct * AREA;

  // Generation of gX for b-b' resistance
  double rBpr = model_.minBaseResist / AREA;
  double rBpi = model_.baseResist / AREA - rBpr;
  double xjrB = model_.baseCurrHalfResist * AREA;

  gX = rBpr + rBpi / qB;

  if (fabs(xjrB) > 0.0)
  {
    double arg1 = Xycemax(iB / xjrB, 1.0e-09);
    double arg2 = (-1.0 + sqrt( 1.0 + 14.59025 * arg1)) / 2.4317 / sqrt(arg1);
    arg1 = tan(arg2);
    gX = rBpr + 3.0 * rBpi * (arg1 - arg2) / arg2 / arg1 / arg1;
  }

  if (fabs(gX) > 0.0) gX = 1.0 / gX;

  diBrdvB = gX;
  diBrdvCp = 0.0;
  diBrdvEp = 0.0;
  diBrdvBp = -gX;

  gBEtot = gBE / tBetaF + gBEleak;
  gBCtot = gBC / tBetaR + gBCleak;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/01/06
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles ()
{
  bool bsuccess = true;

  int i;
  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);

  sprintf(filename,"Q_%s.dat",getName().c_str());

  double time = getSolverState().currTime;
  FILE *fp1;

  if (callsOutputPlot <= 0)
  {
    fp1 = fopen(filename,"w");
  }
  else
  {
    fp1 = fopen(filename,"a");
  }

  if (callsOutputPlot <= 0)
  {
    fprintf(fp1,
            " TITLE = \"Debug Excess Phase data: %s \"\n", getName().c_str());
  }

  if (callsOutputPlot <= 0)
  {
    fprintf(fp1,"%s","\tVARIABLES = \"TIME (S)\",\n");

    fprintf(fp1,"%s","\t    \"iBE/qB \",\n");
    fprintf(fp1,"%s","\t    \"currCexbc \",\n");
    fprintf(fp1,"%s","\t    \"lastCexbc \",\n");
    if (getDeviceOptions().newExcessPhase)
    {
      fprintf(fp1,"%s","\t    \"i_fx \",\n");
      fprintf(fp1,"%s","\t    \"di_fx \",\n");
    }

    fprintf(fp1,"%s","\tZONE F=POINT T=\"Excess Phase Data\"\n");
  }

  fprintf(fp1,"  %12.4e",time);
  fprintf(fp1,"  %12.4e",(iBE/qB));
  fprintf(fp1,"  %12.4e",currCexbc);
  fprintf(fp1,"  %12.4e",nextCexbc);

  if (getDeviceOptions().newExcessPhase)
  {
    double * solVec = extData.nextSolVectorRawPtr;
    double i_fx = solVec[li_Ifx];
    double di_fx = solVec[li_dIfx];
    fprintf(fp1,"  %12.4e",i_fx);
    fprintf(fp1,"  %12.4e",di_fx);
  }

  fprintf(fp1,"%s","\n");

  ++callsOutputPlot;
  fclose(fp1);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams ()
{

  if (!leakBECurrentGiven && c2Given)
    leakBECurrent = c2 * satCur;

  if (!leakBCCurrentGiven && c4Given)
    leakBCCurrent = c4 * satCur;

  if (!minBaseResistGiven)
    minBaseResist = baseResist;

  if (VAFgiven && ( earlyVoltF != 0.0 ) )
    invEarlyVoltF = 1.0 / earlyVoltF;
  else
    invEarlyVoltF = 0.0;

  if( IKFgiven && ( rollOffF != 0.0 ) )
    invRollOffF = 1.0 / rollOffF;
  else
    invRollOffF = 0.0;

  if( VARgiven && ( earlyVoltR != 0.0 ) )
    invEarlyVoltR = 1.0 / earlyVoltR;
  else
    invEarlyVoltR = 0.0;

  if( IKRgiven && ( rollOffR != 0.0 ) )
    invRollOffR = 1.0 / rollOffR;
  else
    invRollOffR = 0.0;

  if( collectorResist != 0.0 )
    collectorConduct = 1.0 / collectorResist;
  else
    collectorConduct = 0;

  if( emitterResist != 0.0 )
    emitterConduct = 1.0 / emitterResist;
  else
    emitterConduct = 0;

  if( given("VTF") && ( transTimeFVBC != 0.0 ) )
    transTimeVBCFac = 1.0 / ( transTimeFVBC * 1.44 );
  else
    transTimeVBCFac = 0.0;

  double mpi = M_PI;
  excessPhaseFac = ( excessPhase / ( 180.0 / mpi ) ) * transTimeF;

  if( given("FC") )
  {
    if( depCapCoeff > 0.9999 )
    {
      depCapCoeff = 0.9999;
      Xyce::dout() << "Bad Depletion Capacitance Coefficient" << std::endl;
    }
  }
  else
  {
    depCapCoeff = 0.5;
  }

  double xfc = log( 1.0 - depCapCoeff );

  f2 = exp( ( 1.0 + juncExpBE ) * xfc );
  f3 = 1.0 - depCapCoeff * ( 1.0 + juncExpBE );
  f6 = exp( ( 1.0 + juncExpBC ) * xfc );
  f7 = 1.0 - depCapCoeff * ( 1.0 + juncExpBC );

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
bool Model::processInstanceParams()
{

  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : modelblock constructor
// Special Notes :
// Scope         : public
// Creator       : Laura Boucheron
// Creation Date : 7/12/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    TYPE(1),
    TNOM(300.0),
    satCur(1.0e-16),
    betaF(100.0),
    BFgiven(false),
    BFMgiven(false),
    emissionCoeffF(1.0),
    earlyVoltF(1e99),
    VAgiven(false),
    VAFgiven(false),
    VBFgiven(false),
    rollOffF(1e99),
    IKFgiven(false),
    IKgiven(false),
    JBFgiven(false),
    leakBECurrent(0.0),
    leakBEEmissionCoeff(1.5),
    NEgiven(false),
    NLEgiven(false),
    betaR(1.0),
    BRgiven(false),
    BRMgiven(false),
    emissionCoeffR(1.0),
    earlyVoltR(1e99),
    VARgiven(false),
    VBgiven(false),
    VRBgiven(false),
    BVgiven(false),
    rollOffR(1e99),
    IKRgiven(false),
    JBRgiven(false),

    leakBCCurrent(0.0),
    leakBCEmissionCoeff(2.0),
    baseResist(0.0),
    baseCurrHalfResist(0.0),
    IRBgiven(false),
    JRBgiven(false),
    IOBgiven(false),
    minBaseResist(baseResist),
    emitterResist(0.0),
    collectorResist(0.0),
    depCapBE(0.0),
    potBE(0.75),
    VJEgiven(false),
    PEgiven(false),

    juncExpBE(0.33),
    MJEgiven(false),
    MEgiven(false),

    transTimeF(0.0),
    transTimeBiasCoeffF(0.0),
    transTimeFVBC(1.e99),
    transTimeHighCurrF(0.0),
    ITFgiven(false),
    JTFgiven(false),

    excessPhase(0.0),
    depCapBC(0.0),
    potBC(0.75),
    VJCgiven(false),
    PCgiven(false),

    juncExpBC(0.33),
    MJCgiven(false),
    MCgiven(false),

    baseFracBCCap(1.0),
    XCJCgiven(false),
    CDISgiven(false),
    transTimeR(0.0),
    CJS(0.0),
    CJSgiven(false),
    CCSgiven(false),
    CSUBgiven(false),
    potSubst(0.75),
    VJSgiven(false),
    PSgiven(false),
    PSUBgiven(false),

    expSubst(0.0),
    MJSgiven(false),
    MSgiven(false),
    ESUBgiven(false),

    betaExp(0.0),
    XTBgiven(false),
    TBgiven(false),
    TCBgiven(false),

    energyGap(1.11),
    tempExpIS(3.0),
    XTIgiven(false),
    PTgiven(false),

    depCapCoeff(0.5),
    fNCoeff(0.0),
    fNExp(1.0),
    c2(0.0),
    c4(0.0),
    rollOffExp(0.5),
    NKgiven(false),
    NKFgiven(false),

    leakBCCurrentGiven(false),
    JLCgiven(false),
    leakBECurrentGiven(false),
    JLEgiven(false),
    minBaseResistGiven(false),
    c2Given(false),
    c4Given(false)
{
  // Default value is type=1 (NPN) so only check for pnp
  if (getType() == "pnp" || getType() == "PNP")
  {
    TYPE = -1;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // trap for redunantly specified parameters (some can be specified by several different names)
  // forward beta (BF/BFM)
  if (BFgiven && BFMgiven)
  {
    UserError0(*this) << "Both BF and BFM are set, which is redundant.";
  }
  if (BFgiven || BFMgiven)
  {
    BFgiven = true; // this is the given flag that is actually used later.
  }

  // reverse beta (Re/BRM)
  if (BRgiven && BRMgiven)
  {
    UserError0(*this) << "Both BR and BRM are set, which is redundant.";
  }
  if (BRgiven || BRMgiven)
  {
    BRgiven = true; // this is the given flag that is actually used later.
  }

  // forward early voltage  (VA/VAF/VBF)
  int VAFgivenCount = 0;
  VAFgivenCount += VAgiven?1:0;
  VAFgivenCount += VAFgiven?1:0;
  VAFgivenCount += VBFgiven?1:0;
  if (VAFgivenCount > 1)
  {
    UserError0(*this) << "The forward early voltage is set more than once.  VA, VAF and VBF are aliases.";
  }
  if (VAFgivenCount > 0)
  {
    VAFgiven = true; // this is the given flag that is actually used later.
  }

  // reverse early voltage  (VAR/VB/VRB/BV)
  int VARgivenCount = 0;
  VARgivenCount += VARgiven?1:0;
  VARgivenCount += VRBgiven?1:0;
  VARgivenCount += VBgiven?1:0;
  VARgivenCount += BVgiven?1:0;
  if (VARgivenCount > 1)
  {
    UserError0(*this) << "The reverse early voltage is set more than once.  VAR,VB,VRB and BV are aliases.";
  }
  if (VARgivenCount > 0)
  {
    VARgiven = true; // this is the given flag that is actually used later.
  }

  // high current roll-off (IKF/IK/JBF)
  int IKFgivenCount = 0;
  IKFgivenCount += IKgiven?1:0;
  IKFgivenCount += IKFgiven?1:0;
  IKFgivenCount += JBFgiven?1:0;
  if (IKFgivenCount > 1)
  {
    UserError0(*this) << "High current roll-off is set more than once (IKF, JBF or IK).";
  }
  if (IKFgivenCount > 0)
  {
    IKFgiven = true; // this is the given flag that is actually used later.
  }

  // reverse high current roll-off (IKR/JBR)
  if (IKRgiven && JBRgiven)
  {
    UserError0(*this) << "Both IKR and JBRgiven are set, which is redundant.";
  }
  if (IKRgiven || JBRgiven)
  {
    IKRgiven = true; // this is the flag actually used later.
  }

  // BE leakage saturation current (ISE/JLE).   C2 is sometimes thougt of as equivalent,
  // but it is not.
  if (JLEgiven && leakBECurrentGiven)
  {
    UserError0(*this) << "Both JLE and ISE are set, which is redundant.";
  }
  if (JLEgiven || leakBECurrentGiven)
  {
    leakBECurrentGiven = true; // this is the flag actually used later.
  }

  // BC leakage saturation current (ISC/JLC).   C4 is sometimes thougt of as equivalent,
  // but it is not.
  if (JLCgiven && leakBCCurrentGiven)
  {
    UserError0(*this) << "Both JLC and ISC are set, which is redundant.";
  }
  if (JLCgiven || leakBCCurrentGiven)
  {
    leakBCCurrentGiven = true; // this is the flag actually used later.
  }

  // leakage emission coefficient (NE/NLE)
  if (NLEgiven && NEgiven)
  {
    UserError0(*this) << "Both NLE and NE are set, which is redundant.";
  }
  if (NLEgiven || NEgiven)
  {
    NEgiven = true; // this is the flag actually used later.
  }

  // BE exponential factor (MJE/ME)
  if (MJEgiven && MEgiven)
  {
    UserError0(*this) << " Both MJE and ME are set, which is redundant.";
  }
  if (MJEgiven || MEgiven)
  {
    MJEgiven = true; // this is the flag actually used later.
  }

  // BC exponential factor (MJC/MC)
  if (MJCgiven && MCgiven)
  {
    UserError0(*this) << "Both MJC and MC are set, which is redundant.";
  }
  if (MJCgiven || MCgiven)
  {
    MJCgiven = true; // this is the flag actually used later.
  }

  // zero-bias coll-subst capacitance (CJS/CCS/CSUB)
  int CJSgivenCount = 0;
  CJSgivenCount += CJSgiven?1:0;
  CJSgivenCount += CCSgiven?1:0;
  CJSgivenCount += CSUBgiven?1:0;
  if (CJSgivenCount > 1)
  {
    UserError0(*this) << "The zero-bias collector-substrate capacitance (CJS, CCS or CSUB) is set more than once.";
  }
  if (CJSgivenCount > 0)
  {
    CJSgiven = true; // this is the given flag that is actually used later.
  }

  // current for 1/2 base resist. (IRB/JRB/IOB)
  int IRBgivenCount = 0;
  IRBgivenCount += IRBgiven?1:0;
  IRBgivenCount += JRBgiven?1:0;
  IRBgivenCount += IOBgiven?1:0;
  if (IRBgivenCount > 1)
  {
    UserError0(*this) << "The current for 1/2 base resistance (IRB, JRB or IOB) is set more than once.";
  }
  if (IRBgivenCount > 0)
  {
    IRBgiven = true; // this is the given flag that is actually used later.
  }

  // BE built-in potential (VJE/PE)
  if (VJEgiven && PEgiven)
  {
    UserError0(*this) << "The BE built-in potential (VJE or PE) is set more than once.";
  }
  if (VJEgiven && PEgiven)
  {
    VJEgiven = true; // this is the given flag that is actually used later.
  }

  // BC built-in potential (VJC/PC)
  if (VJCgiven && PCgiven)
  {
    UserError0(*this) << "The BC built-in potential (VJC or PC) is set more than once.";
  }
  if (VJCgiven && PCgiven)
  {
    VJCgiven = true; // this is the given flag that is actually used later.
  }

  // fraction of BC capacitance to int. base node (XCJC/CDIS)
  if (XCJCgiven && CDISgiven)
  {
    UserError0(*this) << "XCJC and CDIS are both set (they are aliases).";
  }
  if (XCJCgiven && CDISgiven)
  {
    XCJCgiven = true; // this is the given flag that is actually used later.
  }

  // substrate junction built-in potential (VJS/PS/PSUB)
  int VJSgivenCount = 0;
  VJSgivenCount += PSgiven?1:0;
  VJSgivenCount += VJSgiven?1:0;
  VJSgivenCount += PSUBgiven?1:0;
  if (VJSgivenCount > 1)
  {
    UserError0(*this) << "The forward early voltage is set more than once.  PS, VJS and PSUB are aliases.";
  }
  if (VJSgivenCount > 0)
  {
    VJSgiven = true; // this is the given flag that is actually used later.
  }

  // substrate junction built-in potential (MJS/MS/ESUB)
  int MJSgivenCount = 0;
  MJSgivenCount += MSgiven?1:0;
  MJSgivenCount += MJSgiven?1:0;
  MJSgivenCount += ESUBgiven?1:0;
  if (MJSgivenCount > 1)
  {
    UserError0(*this) << "The forward early voltage is set more than once.  MS, MJS and ESUB are aliases.";
  }
  if (MJSgivenCount > 0)
  {
    MJSgiven = true; // this is the given flag that is actually used later.
  }

  // (ITF/JTF)
  if (ITFgiven && JTFgiven)
  {
    UserError0(*this) << "ITF and JTF are both set (they are aliases).";
  }
  if (ITFgiven && JTFgiven)
  {
    ITFgiven = true; // this is the given flag that is actually used later.
  }

  // (NK/NKF)
  if (NKgiven && NKFgiven)
  {
    UserError0(*this) << "NK and NKF are both set (they are aliases).";
  }
  if (NKgiven && NKFgiven)
  {
    NKgiven = true; // this is the given flag that is actually used later.
  }

  // beta temperature exponent (XTB/TB/TCB)
  int XTBgivenCount = 0;
  XTBgivenCount += TBgiven?1:0;
  XTBgivenCount += XTBgiven?1:0;
  XTBgivenCount += TCBgiven?1:0;
  if (XTBgivenCount > 1)
  {
    UserError0(*this) << "The forward early voltage is set more than once.  TB, XTB and TCB are aliases.";
  }
  if (XTBgivenCount > 0)
  {
    XTBgiven = true; // this is the given flag that is actually used later.
  }


  // (PT/XTI)
  if (PTgiven && XTIgiven)
  {
    UserError0(*this) << "PT and XTI are both set (they are aliases).";
  }
  if (PTgiven && XTIgiven)
  {
    XTIgiven = true; // this is the given flag that is actually used later.
  }

  // end of redundant parameter traps
  ////////////////////////////////////////////////////////////////////////////////


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    TNOM = getDeviceOptions().tnom;

  if( !NKgiven)
    rollOffExp = 0.5;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::~Model()
{
  std::vector<Instance*>::iterator iterI;
  std::vector<Instance*>::iterator firstI = instanceContainer.begin ();
  std::vector<Instance*>::iterator lastI  = instanceContainer.end ();

  // loop over instances:
  for (iterI = firstI; iterI != lastI; ++iterI)
  {
    delete (*iterI);
  }

}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Laura Boucheron
// Creation Date : 7/11/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << "     name     modelName  Parameters" << std::endl;
  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "        ";
    os << getName();

    os << std::endl;
    os << "  AREA  = " << (*iter)->AREA  << std::endl;
    os << "  icVBE = " << (*iter)->icVBE << std::endl;
    os << "  icVCE = " << (*iter)->icVCE << std::endl;
    os << "  TEMP  = " << (*iter)->TEMP  << std::endl;
    os << "  OFF   = " << (*iter)->OFF   << std::endl;

    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


//-----------------------------------------------------------------------------
// BJT Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);

    double * stoVec = bi.extData.nextStoVectorRawPtr;

    // Do the bulk of the work in updateIntermediateVars:
    bool btmp = bi.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    // save voltage drops
    stoVec[bi.li_storevBE] = bi.vBE;
    stoVec[bi.li_storevBC] = bi.vBC;
    stoVec[bi.li_store_capeqCB] = bi.capeqCB;

    staVec[bi.li_qstateBEdiff] = bi.qBEdiff;
    staVec[bi.li_qstateBEdep] = bi.qBEdep;
    staVec[bi.li_qstateBCdiff] = bi.qBCdiff;
    staVec[bi.li_qstateBCdep] = bi.qBCdep;
    staVec[bi.li_qstateBX] = bi.qBX;
    staVec[bi.li_qstateCS] = bi.qCS;


    // if this is the first newton step of the first time step
    // of the transient simulation, we need to enforce that the
    // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
    // compatibility.  ERK.

    if (!(getSolverState().dcopFlag) && (getSolverState().initTranFlag) && getSolverState().newtonIter==0)
    {
      double * currStaVec = (bi.extData.currStaVectorRawPtr);

      currStaVec[bi.li_qstateBEdiff] = bi.qBEdiff;
      currStaVec[bi.li_qstateBEdep] = bi.qBEdep;
      currStaVec[bi.li_qstateBCdiff] = bi.qBCdiff;
      currStaVec[bi.li_qstateBCdep] = bi.qBCdep;
      currStaVec[bi.li_qstateBX] = bi.qBX;
      currStaVec[bi.li_qstateCS] = bi.qCS;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);

    // Now that the state vector for time=0 is up-to-date, get the derivative
    // with respect to time of the charge, to obtain the best estimate for
    // the current in the capacitor.

    bi.iBEdiff = staDerivVec[bi.li_qstateBEdiff];
    bi.iBEdep  = staDerivVec[bi.li_qstateBEdep ];
    bi.iCS     = staDerivVec[bi.li_qstateCS    ];
    bi.iBCdiff = staDerivVec[bi.li_qstateBCdiff];
    bi.iBCdep  = staDerivVec[bi.li_qstateBCdep ];
    bi.iBX     = staDerivVec[bi.li_qstateBX    ];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);

    double td = bi.model_.excessPhaseFac;
    double vbe_diff = bi.vBE - bi.vBE_orig;
    double vbc_diff = bi.vBC - bi.vBC_orig;
    double vce_diff = vbe_diff - vbc_diff;

    // F-vector:
    fVec[bi.li_Coll] -= -bi.vCCp * bi.gCpr;
    fVec[bi.li_Base] -= -bi.vBBp * bi.gX;
    fVec[bi.li_Emit] -= -bi.vEEp * bi.gEpr;
    fVec[bi.li_CollP] -= bi.vCCp * bi.gCpr + bi.model_.TYPE * ( - bi.iC );
    fVec[bi.li_BaseP] -= bi.vBBp * bi.gX - bi.model_.TYPE * ( bi.iB );
    fVec[bi.li_EmitP] -= bi.vEEp * bi.gEpr + bi.model_.TYPE * ( - bi.iE );

    // excess phase ERK-dcop.
    double i_fx = 0.0;
    double di_fx = 0.0;

    if (getDeviceOptions().newExcessPhase)
    {
      i_fx = solVec[bi.li_Ifx];
      di_fx = solVec[bi.li_dIfx];

      if (td != 0)
      {
        if (!(getSolverState().dcopFlag) )
        {
          // omega0 = 1/td;
          fVec[bi.li_Ifx] += - di_fx;
          fVec[bi.li_dIfx] += 3 * di_fx*td + 3*i_fx -3 * bi.iBE / bi.qB;
        }
        else
        {
          fVec[bi.li_Ifx] += i_fx -bi.iBE/bi.qB;
          fVec[bi.li_dIfx] = 0.0;
        }
      }
      else
      {
        fVec[bi.li_Ifx] += i_fx;
        fVec[bi.li_dIfx] += di_fx;
      }
    }

    // Q-vector:
    qVec[bi.li_Base] -= - bi.model_.TYPE * bi.qBX;
    qVec[bi.li_Subst] -= -bi.model_.TYPE * bi.qCS;
    qVec[bi.li_CollP] -= bi.model_.TYPE * ( bi.qCS + bi.qBX + bi.qBCdep + bi.qBCdiff );
    qVec[bi.li_BaseP] -= -bi.model_.TYPE * ( bi.qBEdep + bi.qBEdiff + bi.qBCdep + bi.qBCdiff );
    qVec[bi.li_EmitP] -= bi.model_.TYPE*( bi.qBEdep + bi.qBEdiff );

    // excess phase ERK-dcop
    if (td != 0 && getDeviceOptions().newExcessPhase)
    {
      qVec[bi.li_Ifx] += solVec[bi.li_Ifx];

      if (!(getSolverState().dcopFlag) )
      {
        qVec[bi.li_dIfx] += solVec[bi.li_dIfx]*td*td;
      }
      else
      {
        qVec[bi.li_dIfx] = 0.0;
      }
    }

    // voltage limiter terms:
    if (getDeviceOptions().voltageLimiterFlag)
    {
      double Cp_Jdxp_f(0.0), Bp_Jdxp_f(0.0), Ep_Jdxp_f(0.0),
             dIfx_Jdxp_f(0.0 ), Ifx_Jdxp_f(0.0 ),
             Cp_Jdxp_q ( 0.0), Ep_Jdxp_q ( 0.0), Bp_Jdxp_q ( 0.0);

      // F-limiters:
      if (!bi.origFlag || bi.offFlag)
      {
        Cp_Jdxp_f = + bi.diCEdvBp * vbe_diff
                    + bi.diCEdvCp * vce_diff
                    - bi.gBCtot  * vbc_diff;

        Cp_Jdxp_f *= bi.model_.TYPE;

        Bp_Jdxp_f =  bi.gBEtot * vbe_diff + bi.gBCtot  * vbc_diff;
        Bp_Jdxp_f *= bi.model_.TYPE;

        Ep_Jdxp_f = - bi.diCEdvCp * vce_diff - (bi.diCEdvBp + bi.gBEtot) * vbe_diff;
        Ep_Jdxp_f *= bi.model_.TYPE;

        // ERK-dcop.
        if ( td != 0 && getDeviceOptions().newExcessPhase )
        {
          if ( !(getSolverState().dcopFlag) )
          {
            dIfx_Jdxp_f =  -3 *( bi.diBEdvBp*vbe_diff + bi.diBEdvCp * vce_diff);
            dIfx_Jdxp_f *= bi.model_.TYPE;
          }
          else
          {
            Ifx_Jdxp_f =  ( bi.diBEdvBp*vbe_diff + bi.diBEdvCp * vce_diff);
            Ifx_Jdxp_f *= bi.model_.TYPE;
          }
        }
      }

      double * dFdxdVp = bi.extData.dFdxdVpVectorRawPtr;
      dFdxdVp[bi.li_CollP] += Cp_Jdxp_f;
      dFdxdVp[bi.li_BaseP] += Bp_Jdxp_f;
      dFdxdVp[bi.li_EmitP] += Ep_Jdxp_f;

      // ERK-dcop.
      if ( td != 0 && getDeviceOptions().newExcessPhase )
      {
        if ( !(getSolverState().dcopFlag) )
        {
          dFdxdVp[bi.li_dIfx] += dIfx_Jdxp_f;
        }
        else
        {
          dFdxdVp[bi.li_Ifx] += Ifx_Jdxp_f;
        }
      }

      // Q-limiters:
      if (!bi.origFlag || bi.offFlag)
      {
        Cp_Jdxp_q =  -(bi.capBCdep + bi.capBCdiff)*vbc_diff;
        Cp_Jdxp_q *= bi.model_.TYPE;

        Bp_Jdxp_q =  (bi.capBEdep + bi.capBEdiff)*vbe_diff
            + (bi.capBCdiff + bi.capBCdep + bi.capeqCB) *vbc_diff;
        Bp_Jdxp_q *= bi.model_.TYPE;

        Ep_Jdxp_q = - bi.capeqCB * vbc_diff - (bi.capBEdiff + bi.capBEdep)* vbe_diff;
        Ep_Jdxp_q *= bi.model_.TYPE;
      }

      double * dQdxdVp = bi.extData.dQdxdVpVectorRawPtr;
      dQdxdVp[bi.li_CollP] += Cp_Jdxp_q;
      dQdxdVp[bi.li_BaseP] += Bp_Jdxp_q;
      dQdxdVp[bi.li_EmitP] += Ep_Jdxp_q;
    }

    if( bi.loadLeadCurrent )
    {
      storeLeadQ[bi.li_store_dev_ic] = -bi.model_.TYPE * ( bi.qCS + bi.qBX + bi.qBCdep + bi.qBCdiff );
      storeLeadQ[bi.li_store_dev_ib] = bi.model_.TYPE * ( bi.qBX + bi.qBEdep + bi.qBEdiff + bi.qBCdep + bi.qBCdiff );
      storeLeadQ[bi.li_store_dev_ie] = -bi.model_.TYPE*( bi.qBEdep + bi.qBEdiff );
      storeLeadQ[bi.li_store_dev_is] = bi.model_.TYPE * bi.qCS;

      storeLeadF[bi.li_store_dev_ic] = bi.model_.TYPE * ( bi.iC );
      storeLeadF[bi.li_store_dev_is] = 0;
      storeLeadF[bi.li_store_dev_ie] = bi.model_.TYPE * ( bi.iE );
      storeLeadF[bi.li_store_dev_ib] = bi.model_.TYPE * ( bi.iB );

    }

  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);

   double td = bi.model_.excessPhaseFac;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
   // F-matrix:

   *bi.f_CollEquCollNodePtr += bi.gCpr;
   *bi.f_CollEquCollPNodePtr -= bi.gCpr;
   *bi.f_BaseEquBaseNodePtr += bi.diBrdvB;
   *bi.f_BaseEquCollPNodePtr += bi.diBrdvCp;
   *bi.f_BaseEquBasePNodePtr += bi.diBrdvBp;
   *bi.f_BaseEquEmitPNodePtr += bi.diBrdvEp;
   *bi.f_EmitEquEmitNodePtr += bi.gEpr;
   *bi.f_EmitEquEmitPNodePtr -= bi.gEpr;
   *bi.f_CollPEquCollNodePtr -= bi.gCpr;
   *bi.f_CollPEquCollPNodePtr +=bi.diCEdvCp + bi.gBCtot + bi.gCpr;
   *bi.f_CollPEquBasePNodePtr += bi.diCEdvBp - bi.gBCtot;
   *bi.f_CollPEquEmitPNodePtr += bi.diCEdvEp;

   // excess phase terms.  ERK-dcop
   if ( td != 0 && getDeviceOptions().newExcessPhase && !(getSolverState().dcopFlag) )
   {
     *bi.f_CollPEquIfxNodePtr += 1.0 * bi.model_.TYPE;
   }
   *bi.f_BasePEquBaseNodePtr -= bi.diBrdvB;
   *bi.f_BasePEquCollPNodePtr += -bi.diBrdvCp - bi.gBCtot;
   *bi.f_BasePEquBasePNodePtr += -bi.diBrdvBp + bi.gBEtot + bi.gBCtot;
   *bi.f_BasePEquEmitPNodePtr += -bi.diBrdvEp - bi.gBEtot;
   *bi.f_EmitPEquEmitNodePtr -= bi.gEpr;
   *bi.f_EmitPEquCollPNodePtr += -bi.diCEdvCp;
   *bi.f_EmitPEquBasePNodePtr += - bi.diCEdvBp - bi.gBEtot;

   // ERK Note:  -bi.diCEdvEp should equal + bi.diCEdvBp + bi.diCEdvCp.  In the
   // original old-DAE form, + bi.diCEdvBp + bi.diCEdvCp is used.  In
   // the more recent new-DAE form, we've used bi.diCEdvEp, but it sometimes
   // results in subtle differences between old- and new-DAE jacobians.
   // This should matter, but .....
   *bi.f_EmitPEquEmitPNodePtr += bi.gBEtot + bi.gEpr + bi.diCEdvBp + bi.diCEdvCp; // this is a test.
       //+= bi.gBEtot + bi.gEpr - bi.diCEdvEp;

   // excess phase terms.  ERK-dcop
   if ( td != 0 && getDeviceOptions().newExcessPhase && !(getSolverState().dcopFlag) )
   {
     *bi.f_EmitPEquIfxNodePtr += -1.0 * bi.model_.TYPE;
   }

   // excess phase terms.  ERK-dcop
   if (getDeviceOptions().newExcessPhase)
   {
     if ( td != 0 )
     {
       if ( !(getSolverState().dcopFlag) )
       {
         *bi.f_IfxEqudIfxNodePtr += -1;

         *bi.f_dIfxEquCollPNodePtr  += -3*bi.diBEdvCp * bi.model_.TYPE;
         *bi.f_dIfxEquBasePNodePtr  += -3*bi.diBEdvBp * bi.model_.TYPE;
         *bi.f_dIfxEquEmitPNodePtr  += -3*bi.diBEdvEp * bi.model_.TYPE;
         *bi.f_dIfxEqudIfxNodePtr  += 3* td;
         *bi.f_dIfxEquIfxNodePtr  += 3;
       }
       else
       {
         *bi.f_IfxEquCollPNodePtr += -bi.diBEdvCp * bi.model_.TYPE;
         *bi.f_IfxEquBasePNodePtr += -bi.diBEdvBp * bi.model_.TYPE;
         *bi.f_IfxEquEmitPNodePtr += -bi.diBEdvEp * bi.model_.TYPE;

         *bi.f_IfxEquIfxNodePtr+= 1;
         *bi.f_dIfxEqudIfxNodePtr += 1;
       }
     }
     else
     {
       *bi.f_IfxEquIfxNodePtr += 1;
       *bi.f_dIfxEqudIfxNodePtr += 1;
     }
   }

   // Q-matrix:

   *bi.q_BaseEquBaseNodePtr += bi.capBX;
   *bi.q_BaseEquCollPNodePtr += -bi.capBX;
   *bi.q_SubstEquSubstNodePtr += bi.capCS;
   *bi.q_SubstEquCollPNodePtr -= bi.capCS;
   *bi.q_CollPEquBaseNodePtr -= bi.capBX;
   *bi.q_CollPEquSubstNodePtr -= bi.capCS;
   *bi.q_CollPEquCollPNodePtr +=  bi.capCS + bi.capBX + bi.capBCdep + bi.capBCdiff;
   *bi.q_CollPEquBasePNodePtr += -bi.capBCdep - bi.capBCdiff;
   *bi.q_BasePEquCollPNodePtr += -bi.capBCdiff - bi.capBCdep - bi.capeqCB;
   *bi.q_BasePEquBasePNodePtr += bi.capBEdiff + bi.capBEdep + bi.capBCdiff + bi.capBCdep + bi.capeqCB;
   *bi.q_BasePEquEmitPNodePtr += -bi.capBEdiff - bi.capBEdep;
   *bi.q_EmitPEquCollPNodePtr += bi.capeqCB;
   *bi.q_EmitPEquBasePNodePtr += -bi.capBEdiff - bi.capBEdep - bi.capeqCB;
   *bi.q_EmitPEquEmitPNodePtr += bi.capBEdiff + bi.capBEdep;

   // excess phase terms.  ERK-dcop
   if ( td != 0 && getDeviceOptions().newExcessPhase )
   {
     if (!(getSolverState().dcopFlag) )
     {
       *bi.q_IfxEquIfxNodePtr += 1;
       *bi.q_dIfxEqudIfxNodePtr += 1*td*td;
     }
   }

#else
   // F-matrix:

   dFdx[bi.li_Coll][bi.ACollEquCollNodeOffset] += bi.gCpr;
   dFdx[bi.li_Coll][bi.ACollEquCollPNodeOffset] -= bi.gCpr;
   dFdx[bi.li_Base][bi.ABaseEquBaseNodeOffset] += bi.diBrdvB;
   dFdx[bi.li_Base][bi.ABaseEquCollPNodeOffset] += bi.diBrdvCp;
   dFdx[bi.li_Base][bi.ABaseEquBasePNodeOffset] += bi.diBrdvBp;
   dFdx[bi.li_Base][bi.ABaseEquEmitPNodeOffset] += bi.diBrdvEp;
   dFdx[bi.li_Emit][bi.AEmitEquEmitNodeOffset] += bi.gEpr;
   dFdx[bi.li_Emit][bi.AEmitEquEmitPNodeOffset] -= bi.gEpr;
   dFdx[bi.li_CollP][bi.ACollPEquCollNodeOffset] -= bi.gCpr;
   dFdx[bi.li_CollP][bi.ACollPEquCollPNodeOffset] +=bi.diCEdvCp + bi.gBCtot + bi.gCpr;
   dFdx[bi.li_CollP][bi.ACollPEquBasePNodeOffset] += bi.diCEdvBp - bi.gBCtot;
   dFdx[bi.li_CollP][bi.ACollPEquEmitPNodeOffset] += bi.diCEdvEp;

   // excess phase terms.  ERK-dcop
   if ( td != 0 && getDeviceOptions().newExcessPhase && !(getSolverState().dcopFlag) )
   {
     dFdx[bi.li_CollP][bi.ACollPEquIfxNodeOffset] += 1.0 * bi.model_.TYPE;
   }

   dFdx[bi.li_BaseP][bi.ABasePEquBaseNodeOffset] -= bi.diBrdvB;
   dFdx[bi.li_BaseP][bi.ABasePEquCollPNodeOffset] += -bi.diBrdvCp - bi.gBCtot;
   dFdx[bi.li_BaseP][bi.ABasePEquBasePNodeOffset] += -bi.diBrdvBp + bi.gBEtot + bi.gBCtot;
   dFdx[bi.li_BaseP][bi.ABasePEquEmitPNodeOffset] += -bi.diBrdvEp - bi.gBEtot;
   dFdx[bi.li_EmitP][bi.AEmitPEquEmitNodeOffset] -= bi.gEpr;
   dFdx[bi.li_EmitP][bi.AEmitPEquCollPNodeOffset] += -bi.diCEdvCp;
   dFdx[bi.li_EmitP][bi.AEmitPEquBasePNodeOffset] += - bi.diCEdvBp - bi.gBEtot;

   // ERK Note:  -bi.diCEdvEp should equal + bi.diCEdvBp + bi.diCEdvCp.  In the
   // original old-DAE form, + bi.diCEdvBp + bi.diCEdvCp is used.  In
   // the more recent new-DAE form, we've used bi.diCEdvEp, but it sometimes
   // results in subtle differences between old- and new-DAE jacobians.
   // This should matter, but .....

   dFdx[bi.li_EmitP][bi.AEmitPEquEmitPNodeOffset] += bi.gBEtot + bi.gEpr + bi.diCEdvBp + bi.diCEdvCp; // this is a test.
       //+= bi.gBEtot + bi.gEpr - bi.diCEdvEp;

   // excess phase terms.  ERK-dcop
   if ( td != 0 && getDeviceOptions().newExcessPhase && !(getSolverState().dcopFlag) )
   {
     dFdx[bi.li_EmitP][bi.AEmitPEquIfxNodeOffset] += -1.0 * bi.model_.TYPE;
   }

   // excess phase terms.  ERK-dcop
   if (getDeviceOptions().newExcessPhase)
   {
     if ( td != 0 )
     {
       if ( !(getSolverState().dcopFlag) )
       {
         dFdx[bi.li_Ifx][bi.AIfxEqudIfxNodeOffset]+= -1;
         dFdx[bi.li_dIfx][bi.AdIfxEquCollPNodeOffset] += -3*bi.diBEdvCp * bi.model_.TYPE;
         dFdx[bi.li_dIfx][bi.AdIfxEquBasePNodeOffset] += -3*bi.diBEdvBp * bi.model_.TYPE;
         dFdx[bi.li_dIfx][bi.AdIfxEquEmitPNodeOffset] += -3*bi.diBEdvEp * bi.model_.TYPE;
         dFdx[bi.li_dIfx][bi.AdIfxEqudIfxNodeOffset] += 3* td;
         dFdx[bi.li_dIfx][bi.AdIfxEquIfxNodeOffset] += 3;
       }
       else
       {
         dFdx[bi.li_Ifx][bi.AIfxEquCollPNodeOffset] += -bi.diBEdvCp * bi.model_.TYPE;
         dFdx[bi.li_Ifx][bi.AIfxEquBasePNodeOffset] += -bi.diBEdvBp * bi.model_.TYPE;
         dFdx[bi.li_Ifx][bi.AIfxEquEmitPNodeOffset] += -bi.diBEdvEp * bi.model_.TYPE;
         dFdx[bi.li_Ifx][bi.AIfxEquIfxNodeOffset]+= 1;
         dFdx[bi.li_dIfx][bi.AdIfxEqudIfxNodeOffset] += 1;
       }
     }
     else
     {
       dFdx[bi.li_Ifx][bi.AIfxEquIfxNodeOffset]+= 1;
       dFdx[bi.li_dIfx][bi.AdIfxEqudIfxNodeOffset] += 1;
     }
   }

   // Q-matrix:

   dQdx[bi.li_Base][bi.ABaseEquBaseNodeOffset] += bi.capBX;
   dQdx[bi.li_Base][bi.ABaseEquCollPNodeOffset] += -bi.capBX;
   dQdx[bi.li_Subst][bi.ASubstEquSubstNodeOffset] += bi.capCS;
   dQdx[bi.li_Subst][bi.ASubstEquCollPNodeOffset] -= bi.capCS;
   dQdx[bi.li_CollP][bi.ACollPEquBaseNodeOffset] -= bi.capBX;
   dQdx[bi.li_CollP][bi.ACollPEquSubstNodeOffset] -= bi.capCS;
   dQdx[bi.li_CollP][bi.ACollPEquCollPNodeOffset] +=  bi.capCS + bi.capBX + bi.capBCdep + bi.capBCdiff;
   dQdx[bi.li_CollP][bi.ACollPEquBasePNodeOffset] += -bi.capBCdep - bi.capBCdiff;
   dQdx[bi.li_BaseP][bi.ABasePEquCollPNodeOffset] += -bi.capBCdiff - bi.capBCdep - bi.capeqCB;
   dQdx[bi.li_BaseP][bi.ABasePEquBasePNodeOffset] += bi.capBEdiff + bi.capBEdep + bi.capBCdiff + bi.capBCdep + bi.capeqCB;
   dQdx[bi.li_BaseP][bi.ABasePEquEmitPNodeOffset] += -bi.capBEdiff - bi.capBEdep;
   dQdx[bi.li_EmitP][bi.AEmitPEquCollPNodeOffset] += bi.capeqCB;
   dQdx[bi.li_EmitP][bi.AEmitPEquBasePNodeOffset] += -bi.capBEdiff - bi.capBEdep - bi.capeqCB;
   dQdx[bi.li_EmitP][bi.AEmitPEquEmitPNodeOffset] += bi.capBEdiff + bi.capBEdep;

   // excess phase terms.  ERK-dcop
   if ( td != 0 && getDeviceOptions().newExcessPhase )
   {
     if (!(getSolverState().dcopFlag) )
     {
       dQdx[bi.li_Ifx][bi.AIfxEquIfxNodeOffset] += 1;
       dQdx[bi.li_dIfx][bi.AdIfxEqudIfxNodeOffset] += 1*td*td;
     }
   }
  #endif
  }
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("q", 1)
    .registerModelType("pnp", 1)
    .registerModelType("npn", 1);
}

} // namespace BJT
} // namespace Device
} // namespace Xyce
