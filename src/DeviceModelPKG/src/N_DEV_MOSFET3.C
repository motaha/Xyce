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
// Filename       : $RCSfile: N_DEV_MOSFET3.C,v $
//
// Purpose        : Implement the MOSFET Level 3 static model
//
// Special Notes  : BADMOS3 option not supported.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.218.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
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

// ---------- Xyce Includes ----------
#include <N_DEV_Const.h>
#include <N_DEV_MOSFET3.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<MOSFET3::Instance>::ParametricData()
{
    // Set up configuration constants:
    setNumNodes(4);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(1);
    addModelType("NMOS");
    addModelType("PMOS");

    // Set up double precision variables:
    addPar ("TEMP", 0.0, false, ParameterType::TIME_DEP,
      &MOSFET3::Instance::temp, NULL,
       STANDARD, CAT_NONE, ""),

    addPar ("L",    0.0, true,    ParameterType::NO_DEP,
      &MOSFET3::Instance::l, NULL,
      U_METER,  CAT_GEOMETRY, "Channel length"),

    addPar ("W",    0.0, true,    ParameterType::NO_DEP,
      &MOSFET3::Instance::w, NULL,
      U_METER,  CAT_GEOMETRY, "Channel width"),

    addPar ("AD",   0.0, false,   ParameterType::NO_DEP,
      &MOSFET3::Instance::drainArea, NULL,
      U_METER2,  CAT_GEOMETRY, "Drain diffusion area"),

    addPar ("AS",   0.0, false,   ParameterType::NO_DEP,
      &MOSFET3::Instance::sourceArea, NULL,
      U_METER2,  CAT_GEOMETRY, "Source diffusion area"),

    addPar ("NRD",  1.0, false,   ParameterType::NO_DEP,
      &MOSFET3::Instance::drainSquares, NULL,
      U_SQUARES, CAT_GEOMETRY, "Multiplier for RSH to yield parasitic resistance of drain"),

    addPar ("NRS",  1.0, false,   ParameterType::NO_DEP,
      &MOSFET3::Instance::sourceSquares, NULL,
      U_SQUARES, CAT_GEOMETRY, "Multiplier for RSH to yield parasitic resistance of source"),

    addPar ("PD",   0.0, false,   ParameterType::NO_DEP,
      &MOSFET3::Instance::drainPerimeter, NULL,
      U_METER, CAT_GEOMETRY, "Drain diffusion perimeter"),

    addPar ("PS",   0.0, false,   ParameterType::NO_DEP,
      &MOSFET3::Instance::sourcePerimeter, NULL,
      U_METER, CAT_GEOMETRY, "Source diffusion perimeter"),

    addPar ("M",    1.0, false,   ParameterType::NO_DEP,
      &MOSFET3::Instance::numberParallel, NULL,
      U_NONE, CAT_CONTROL,  "Multiplier for M devices connected in parallel"),

    // Initial conditions 
    addPar ("IC1", 0.0, false, NO_DEP,
            &MOSFET3::Instance::icVDS,
            &MOSFET3::Instance::IC_GIVEN,
            U_VOLT, CAT_INITIAL, "Initial condition on Drain-Source voltage");

    addPar ("IC2", 0.0, false, NO_DEP,
            &MOSFET3::Instance::icVGS,
            &MOSFET3::Instance::IC_GIVEN,
            U_VOLT, CAT_INITIAL, "Initial condition on Gate-Source voltage");

    addPar ("IC3", 0.0, false, NO_DEP,
            &MOSFET3::Instance::icVBS,
            &MOSFET3::Instance::IC_GIVEN,
            U_VOLT, CAT_INITIAL, "Initial condition on Bulk-Source voltage");

    makeVector ("IC",3);

    // Set up non-double precision variables:
    addPar ("OFF",false,false, ParameterType::NO_DEP,
            &MOSFET3::Instance::OFF,
            NULL, U_LOGIC, CAT_VOLT,
            "Initial condition of no voltage drops across device");
}

template<>
ParametricData<MOSFET3::Model>::ParametricData()
{
    // Set up double precision variables:
    addPar ("L",      1e-4,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::model_l, NULL,
      U_METER, CAT_GEOMETRY, "Default channel length"),

    addPar ("W",      1e-4,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::model_w, NULL,
      U_METER, CAT_GEOMETRY, "Default channel width"),

    addPar ("VTO",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::vt0, NULL,
      U_VOLT, CAT_VOLT, "Zero-bias threshold voltage"),

    addPar ("KP",     2e-5, false, ParameterType::NO_DEP,
      &MOSFET3::Model::transconductance, NULL,
      U_AMPVM2, CAT_PROCESS, "Transconductance coefficient"),

    addPar ("GAMMA",  0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::gamma, NULL,
      U_VOLTH, CAT_PROCESS, "Bulk threshold parameter"),

    addPar ("PHI",    0.6,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::phi, NULL,
      U_VOLT, CAT_PROCESS, "Surface potential"),

    addPar ("RD",     0.0,  false, ParameterType::MIN_RES,
      &MOSFET3::Model::drainResistance, NULL,
      U_OHM, CAT_RES, "Drain ohmic resistance"),

    addPar ("RS",     0.0,  false, ParameterType::MIN_RES,
      &MOSFET3::Model::sourceResistance, NULL,
      U_OHM, CAT_RES, "Source ohmic resistance"),

    addPar ("CBD",    0.0,  false, ParameterType::MIN_CAP,
      &MOSFET3::Model::capBD,
      &MOSFET3::Model::capBDGiven,
      U_FARAD, CAT_CAP, "Zero-bias bulk-drain p-n capacitance"),

    addPar ("CBS",    0.0,  false, ParameterType::MIN_CAP,
      &MOSFET3::Model::capBS,
      &MOSFET3::Model::capBSGiven,
      U_FARAD, CAT_CAP, "Zero-bias bulk-source p-n capacitance"),

    addPar ("IS",     1e-14,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::jctSatCur, NULL,
      U_AMP, CAT_CURRENT, "Bulk p-n saturation current"),

    addPar ("PB",     0.8,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::bulkJctPotential, NULL,
      U_VOLT, CAT_VOLT, "Bulk p-n bottom potential"),

    addPar ("CGSO",   0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::gateSourceOverlapCapFactor, NULL,
      U_FARADMM1, CAT_CAP, "Gate-source overlap capacitance/channel width"),

    addPar ("CGDO",   0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::gateDrainOverlapCapFactor, NULL,
      U_FARADMM1, CAT_CAP, "Gate-drain overlap capacitance/channel width"),

    addPar ("CGBO",   0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::gateBulkOverlapCapFactor, NULL,
      U_FARADMM1, CAT_CAP, "Gate-bulk overlap capacitance/channel length"),

    addPar ("RSH",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::sheetResistance, NULL,
      U_OHM, CAT_RES, "Drain, source diffusion sheet resistance"),

    addPar ("CJ",     0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::bulkCapFactor,
      &MOSFET3::Model::bulkCapFactorGiven,
      U_FARADMM2, CAT_CAP, "Bulk p-n zero-bias bottom capacitance/area"),

    addPar ("MJ",     0.5,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::bulkJctBotGradingCoeff, NULL,
      U_NONE, CAT_DOPING, "Bulk p-n bottom grading coefficient"),

    addPar ("CJSW",   0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::sideWallCapFactor,
      &MOSFET3::Model::sideWallCapFactorGiven,
      U_FARADMM2, CAT_CAP, "Bulk p-n zero-bias sidewall capacitance/area"),

    addPar ("MJSW",   0.33,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::bulkJctSideGradingCoeff, NULL,
      U_NONE, CAT_DOPING, "Bulk p-n sidewall grading coefficient"),

    addPar ("JS",     0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::jctSatCurDensity, NULL,
      U_AMPMM2, CAT_PROCESS, "Bulk p-n saturation current density"),

    addPar ("TOX",    1e-7,  true, ParameterType::NO_DEP,
      &MOSFET3::Model::oxideThickness, NULL,
      U_METER, CAT_GEOMETRY, "Gate oxide thickness"),

    addPar ("LD",     0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::latDiff, NULL,
      U_METER, CAT_DOPING, "Lateral diffusion length"),

    addPar ("UO",     600.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::surfaceMobility, NULL,
      U_CMM2VM1SM1, ParameterCategory(CAT_PROCESS | UNDOCUMENTED), "Surface mobility"),

    addPar ("U0",     600.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::surfaceMobility0, NULL,
      U_CMM2VM1SM1, CAT_PROCESS, "Surface mobility"),

    addPar ("FC",     0.5,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::fwdCapDepCoeff, NULL,
      U_NONE, CAT_CAP, "Bulk p-n forward-bias capacitance coefficient"),

    addPar ("NSUB",   0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::substrateDoping, NULL,
      U_CMM3, CAT_DOPING, "Substrate doping density"),

    addPar ("NSS",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::surfaceStateDensity, NULL,
      U_CMM2, CAT_PROCESS, "Surface state density"),

    addPar ("ETA",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::eta, NULL,
       U_NONE, CAT_PROCESS, "Static feedback"),

    addPar ("DELTA",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::delta, NULL,
       U_NONE, CAT_PROCESS, "Width effect on threshold"),

    addPar ("NFS",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::fastSurfaceStateDensity, NULL,
       U_CMM2, CAT_PROCESS, "Fast surface state density"),

    addPar ("THETA",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::theta, NULL,
       U_VOLTM1, CAT_PROCESS, "Mobility modulation"),

    addPar ("VMAX",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::maxDriftVel, NULL,
       U_MSM1, CAT_PROCESS, "Maximum drift velocity"),

    addPar ("KAPPA",    0.2,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::kappa, NULL,
       U_NONE, CAT_NONE, "Saturation field factor"),

    addPar ("XJ",    0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::junctionDepth, NULL,
       U_METER, CAT_GEOMETRY, "Metallurgical junction depth"),

    addPar ("TNOM",   27.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::tnom, NULL,
      STANDARD, CAT_NONE, ""),

    addPar ("KF",     0.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::fNcoef, NULL,
      U_NONE, CAT_FLICKER, "Flicker noise coefficient"),

    addPar ("AF",     1.0,  false, ParameterType::NO_DEP,
      &MOSFET3::Model::fNexp, NULL,
      U_NONE, CAT_FLICKER, "Flicker noise exponent"),

    // Set up non-double precision variables:
    addPar ("TPG", 1, false, ParameterType::NO_DEP,
      &MOSFET3::Model::gateType, NULL,
      U_NONE, CAT_MATERIAL, "Gate material type (-1 = same as substrate,"
                            " 0 = aluminum, 1 = opposite of substrate)");

    DeviceModel::initThermalModel(*this);
}

namespace MOSFET3 {

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
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{

  // now set the temperature related stuff
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 3/21/01
//-----------------------------------------------------------------------------
Instance::Instance (
  InstanceBlock & IB,
   Model & Miter,
   MatrixLoadData & mlData1,
   SolverState &ss1,
   ExternData  &ed1,
   DeviceOptions &do1)

  : DeviceInstance(IB, mlData1, ss1, ed1,do1),
    model_(Miter),
    dNode(0),
    gNode(0),
    sNode(0),
    bNode(0),
    dNodePrime(0),
    sNodePrime(0),
    l(getDeviceOptions().defl),
    w(getDeviceOptions().defw),
    OFF(false),
    drainArea(getDeviceOptions().defad),
    sourceArea(getDeviceOptions().defas),
    drainSquares(1.0),
    sourceSquares(1.0),
    drainPerimeter(0.0),
    sourcePerimeter(0.0),
    sourceConductance(0.0),
    drainConductance(0.0),
  temp(getDeviceOptions().temp.dVal()),
    numberParallel(1.0),
    tTransconductance(0.0),
    tSurfMob(0.0),
    tPhi(0.0),
    tVto(0.0),
    tSatCur(0.0),
    tSatCurDens(0.0),
    tCbd(0.0),
    tCbs(0.0),
    tCj(0.0),
    tCjsw(0.0),
    tBulkPot(0.0),
    tDepCap(0.0),
    tVbi(0.0),
    icVBS(0.0),
    icVDS(0.0),
    icVGS(0.0),
    von(0.0),
    vdsat(0.0),
    sourceVcrit(0.0),
    drainVcrit(0.0),
    cd(0.0),
    cbs(0.0),
    cbd(0.0),
    gmbs(0.0),
    gm(0.0),
    gds(0.0),
    gbd(0.0),
    gbs(0.0),
    capbd(0.0),
    capbs(0.0),
    Cbd(0.0),
    Cbdsw(0.0),
    Cbs(0.0),
    Cbssw(0.0),
    f2d(0.0),
    f3d(0.0),
    f4d(0.0),
    f2s(0.0),
    f3s(0.0),
    f4s(0.0),
    mode(1),
    mode_low(0.0),
    mode_high(0.0),
    limitedFlag(false),
    IC_GIVEN(false),
    Idrain(0.0),
    Isource(0.0),
    //calculated quantities
    EffectiveLength(0),
    DrainSatCur(0),
    SourceSatCur(0),
    GateSourceOverlapCap(0),
    GateDrainOverlapCap(0),
    GateBulkOverlapCap(0),
    OxideCap(0),
    // local indices
    li_Drain(-1),
    li_DrainPrime(-1),
    li_Source(-1),
    li_SourcePrime(-1),
    li_Gate(-1),
    li_Bulk(-1),
    // matrix and vector pointers:
    // jacobian:
    //  drain row
    ADrainEquDrainNodeOffset(-1),
    ADrainEquDrainPrimeNodeOffset(-1),
    //  gate row
    AGateEquGateNodeOffset(-1),
    AGateEquBulkNodeOffset(-1),
    AGateEquDrainPrimeNodeOffset(-1),
    AGateEquSourcePrimeNodeOffset(-1),
    //  source row
    ASourceEquSourceNodeOffset(-1),
    ASourceEquSourcePrimeNodeOffset(-1),
    //  bulk row
    ABulkEquGateNodeOffset(-1),
    ABulkEquBulkNodeOffset(-1),
    ABulkEquDrainPrimeNodeOffset(-1),
    ABulkEquSourcePrimeNodeOffset(-1),
    // drain' row
    ADrainPrimeEquDrainNodeOffset(-1),
    ADrainPrimeEquGateNodeOffset(-1),
    ADrainPrimeEquBulkNodeOffset(-1),
    ADrainPrimeEquDrainPrimeNodeOffset(-1),
    ADrainPrimeEquSourcePrimeNodeOffset(-1),
    // source' row
    ASourcePrimeEquGateNodeOffset(-1),
    ASourcePrimeEquSourceNodeOffset(-1),
    ASourcePrimeEquBulkNodeOffset(-1),
    ASourcePrimeEquDrainPrimeNodeOffset(-1),
    ASourcePrimeEquSourcePrimeNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // F-vector pointers:
  // V_d Row:
    f_DrainEquDrainNodePtr(0),             // a
    f_DrainEquDrainPrimeNodePtr(0),        // b
  // V_g Row:
    f_GateEquGateNodePtr(0),               // c
    f_GateEquBulkNodePtr(0),               // d
    f_GateEquDrainPrimeNodePtr(0),         // e
    f_GateEquSourcePrimeNodePtr(0),        // f
  // V_s Row:
    f_SourceEquSourceNodePtr(0),           // g
    f_SourceEquSourcePrimeNodePtr(0),      // h
  // V_b Row:
    f_BulkEquGateNodePtr(0),               // i
    f_BulkEquBulkNodePtr(0),               // j
    f_BulkEquDrainPrimeNodePtr(0),         // k
    f_BulkEquSourcePrimeNodePtr(0),        // l
  // V_d' Row:
    f_DrainPrimeEquDrainNodePtr(0),        // m
    f_DrainPrimeEquGateNodePtr(0),         // n
    f_DrainPrimeEquBulkNodePtr(0),         // o
    f_DrainPrimeEquDrainPrimeNodePtr(0),   // p
    f_DrainPrimeEquSourcePrimeNodePtr(0),  // q
  // V_s' Row:
    f_SourcePrimeEquGateNodePtr(0),        // r
    f_SourcePrimeEquSourceNodePtr(0),      // s
    f_SourcePrimeEquBulkNodePtr(0),        // t
    f_SourcePrimeEquDrainPrimeNodePtr(0),  // u
    f_SourcePrimeEquSourcePrimeNodePtr(0), // v
  // Q-vector pointers:
  // V_d Row:
    q_DrainEquDrainNodePtr(0),             // a
    q_DrainEquDrainPrimeNodePtr(0),        // b
  // V_g Row:
    q_GateEquGateNodePtr(0),               // c
    q_GateEquBulkNodePtr(0),               // d
    q_GateEquDrainPrimeNodePtr(0),         // e
    q_GateEquSourcePrimeNodePtr(0),        // f
  // V_s Row:
    q_SourceEquSourceNodePtr(0),           // g
    q_SourceEquSourcePrimeNodePtr(0),      // h
  // V_b Row:
    q_BulkEquGateNodePtr(0),               // i
    q_BulkEquBulkNodePtr(0),               // j
    q_BulkEquDrainPrimeNodePtr(0),         // k
    q_BulkEquSourcePrimeNodePtr(0),        // l
  // V_d' Row:
    q_DrainPrimeEquDrainNodePtr(0),        // m
    q_DrainPrimeEquGateNodePtr(0),         // n
    q_DrainPrimeEquBulkNodePtr(0),         // o
    q_DrainPrimeEquDrainPrimeNodePtr(0),   // p
    q_DrainPrimeEquSourcePrimeNodePtr(0),  // q
  // V_s' Row:
    q_SourcePrimeEquGateNodePtr(0),        // r
    q_SourcePrimeEquSourceNodePtr(0),      // s
    q_SourcePrimeEquBulkNodePtr(0),        // t
    q_SourcePrimeEquDrainPrimeNodePtr(0),  // u
    q_SourcePrimeEquSourcePrimeNodePtr(0), // v
#endif

    vbd(0.0),
    vbs(0.0),
    vgs(0.0),
    vds(0.0),
    vbd_old(0.0),
    vbs_old(0.0),
    vgs_old(0.0),
    vds_old(0.0),
    vbd_orig(0.0),
    vbs_orig(0.0),
    vgs_orig(0.0),
    vds_orig(0.0),
    capgs(0.0),
    qgs(0.0),
    capgd(0.0),
    qgd(0.0),
    capgb(0.0),
    qgb(0.0),
    qbd(0.0),
    qbs(0.0),
    // local indices
    li_store_vbd(-1),
    li_store_vbs(-1),
    li_store_vgs(-1),
    li_store_vds(-1),
    li_store_von(-1),
    li_store_dev_id(-1),
    li_store_dev_is(-1),
    li_store_dev_ig(-1),
    li_store_dev_ib(-1),
    li_state_qgs(-1),
    li_state_qgd(-1),
    li_state_qgb(-1),
    li_state_capgs(-1),
    li_state_capgd(-1),
    li_state_capgb(-1),
    li_state_qbd(-1),
    li_state_qbs(-1),
    blockHomotopyID(0),
    randomPerturb(0.0)
{
  numIntVars   = 2;
  numExtVars   = 4;

  setNumStoreVars(5);
  numLeadCurrentStoreVars = 4; // drain, gate, source & base lead currents
  numStateVars = 8;

  devConMap.resize(4);
  devConMap[0] = 1;
  devConMap[1] = 2;
  devConMap[2] = 1;
  devConMap[3] = 3;

  setName(IB.getName());
  setModelName(model_.getName());

  blockHomotopyID =
    devSupport.getGainScaleBlockID(getDeviceOptions().numGainScaleBlocks);
  randomPerturb =
    devSupport.getRandomPerturbation();

  if( jacStamp.empty() )
  {
    // stamp for RS!=0, RD!=0
    jacStamp_DC_SC.resize(6);
    jacStamp_DC_SC[0].resize(2);  // Drain row
    jacStamp_DC_SC[0][0]=0;       // d-d
    jacStamp_DC_SC[0][1]=4;       // d-d'
    jacStamp_DC_SC[1].resize(4);  // Gate row
    jacStamp_DC_SC[1][0]=1;       // g-g
    jacStamp_DC_SC[1][1]=3;       // g-b
    jacStamp_DC_SC[1][2]=4;       // g-d'
    jacStamp_DC_SC[1][3]=5;       // g-s'
    jacStamp_DC_SC[2].resize(2);  // Source row
    jacStamp_DC_SC[2][0]=2;       // s-s
    jacStamp_DC_SC[2][1]=5;       // s-s'
    jacStamp_DC_SC[3].resize(4);  // Bulk row
    jacStamp_DC_SC[3][0]=1;       // b-g
    jacStamp_DC_SC[3][1]=3;       // b-b
    jacStamp_DC_SC[3][2]=4;       // b-d'
    jacStamp_DC_SC[3][3]=5;       // b-s'
    jacStamp_DC_SC[4].resize(5);  // Drain' row
    jacStamp_DC_SC[4][0]=0;       // d'-d
    jacStamp_DC_SC[4][1]=1;       // d'-g
    jacStamp_DC_SC[4][2]=3;       // d'-b
    jacStamp_DC_SC[4][3]=4;       // d'-d'
    jacStamp_DC_SC[4][4]=5;       // d'-s'
    jacStamp_DC_SC[5].resize(5);  // Source' row
    jacStamp_DC_SC[5][0]=1;       // s'-g
    jacStamp_DC_SC[5][1]=2;       // s'-s
    jacStamp_DC_SC[5][2]=3;       // s'-b
    jacStamp_DC_SC[5][3]=4;       // s'-d'
    jacStamp_DC_SC[5][4]=5;       // s'-s'

    jacMap_DC_SC.clear();
    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_DC,    jacMap_DC, jacMap2_DC, 5, 2, 6);

    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_SC,    jacMap_SC, jacMap2_SC, 4, 0, 6);

    jacStampMap(jacStamp_DC, jacMap_DC, jacMap2_DC,
                jacStamp,    jacMap, jacMap2, 4, 0, 6);
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.dVal();
  if (!given("L"))
    l =model_.model_l;
  if (!given("W"))
    w = model_.model_w;

  if (!given("AD"))
    drainArea = getDeviceOptions().defad;
  if (!given("AS"))
    sourceArea = getDeviceOptions().defas;

  // Note: processParams needs to be before the drain, source resistance calculations.
  processParams ();

  // process source/drain series resistance (from mos1temp)
  if(model_.drainResistance != 0)
  {
    drainConductance = 1/model_.drainResistance;
  }
  else if (model_.given("RSH"))
  {
    if(model_.sheetResistance != 0)
    {
      drainConductance =
        1/(model_.sheetResistance*drainSquares);
    }
    else
    {
      drainConductance = 0;
    }
  }
  else
  {
    drainConductance = 0;
  }
  if(model_.sourceResistance != 0)
  {
    sourceConductance = 1/model_.sourceResistance;
  }
  else if (model_.given("RSH"))
  {
    if(model_.sheetResistance != 0)
    {
      sourceConductance =
        1/(model_.sheetResistance*sourceSquares);
    }
    else
    {
      sourceConductance = 0;
    }
  }
  else
  {
    sourceConductance = 0;
  }

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  if(l - 2 * model_.latDiff <=0)
  {
    string msg = "Effective channel length less than zero.";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  EffectiveLength=l - 2*model_.latDiff;
  GateSourceOverlapCap = model_.gateSourceOverlapCapFactor * w;
  GateDrainOverlapCap = model_.gateDrainOverlapCapFactor * w;
  GateBulkOverlapCap = model_.gateBulkOverlapCapFactor * EffectiveLength;
  OxideCap = model_.oxideCapFactor * EffectiveLength * w;

  numIntVars = (((sourceConductance == 0.0)?0:1)+((drainConductance == 0.0) ? 0:1));
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                                          const vector<int> & extLIDVecRef )
{
#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";
  if (getDeviceOptions().debugLevel > 0)
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
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  number of internal variables: " << numInt << endl;
    cout << "  number of external variables: " << numExt << endl;
  }
#endif

  numIntVars = (((sourceConductance == 0.0)?0:1)+((drainConductance == 0.0) ? 0:1));

  if ( numIntVars != numInt)
  {
    string msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (numExt != numExtVars)
  {
    string msg = "Instance::registerLIDs:";
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

  if( drainConductance )
    li_DrainPrime = intLIDVec[intLoc++];
  else
    li_DrainPrime = li_Drain;

  if( sourceConductance )
    li_SourcePrime = intLIDVec[intLoc];
  else
    li_SourcePrime = li_Source;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "\n variable local indices:\n";
    cout << "  li_Drain       = " << li_Drain << endl;
    cout << "  li_DrainPrime  = " << li_DrainPrime << endl;
    cout << "  li_Source      = " << li_Source << endl;
    cout << "  li_SourcePrime = " << li_SourcePrime << endl;
    cout << "  li_Gate        = " << li_Gate << endl;
    cout << "  li_Bulk        = " << li_Bulk << endl;

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
    // set up the internal name map:
    string tmpstr;
    if ( li_DrainPrime != li_Drain )
    {
      tmpstr = getName()+"_drainprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_DrainPrime ] = tmpstr;
    }

    if ( li_SourcePrime != li_Source )
    {
      tmpstr = getName()+"_sourceprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_SourcePrime ] = tmpstr;
    }
  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_MOSFET3Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 4/3/2013
//-----------------------------------------------------------------------------
map<int,string> & N_DEV_MOSFET3Instance::getStoreNameMap ()
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
// Creator       : Robert Hoekstra, Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
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
    string msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
    cout << "  Number of State LIDs: " << numSta << endl;
#endif

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int lid=0;

  li_state_qgs = staLIDVec[lid++];
  li_state_qgd = staLIDVec[lid++];
  li_state_qgb = staLIDVec[lid++];

  li_state_capgs = staLIDVec[lid++];
  li_state_capgd = staLIDVec[lid++];
  li_state_capgb = staLIDVec[lid++];

  li_state_qbd = staLIDVec[lid++];
  li_state_qbs = staLIDVec[lid++];


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  State local indices:" << endl;
    cout << endl;
    cout << "  li_state_qgs           = " << li_state_qgs ;
    cout << "  li_state_capgs         = " << li_state_capgs;
    cout << "  li_state_capgd         = " << li_state_capgd;
    cout << "  li_state_capgb         = " << li_state_capgb;
    cout << "  li_state_qgd           = " << li_state_qgd;
    cout << "  li_state_qgb           = " << li_state_qgb;
    cout << "  li_state_qbs           = " << li_state_qbs;
    cout << "  li_state_qbd           = " << li_state_qbd;
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
// Creator       : Eric Keiter, SNL
// Creation Date : 12/11/11
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs( const vector<int> & stoLIDVecRef )
{
#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
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
  if (getDeviceOptions().debugLevel > 0)
    cout << "  Number of Store LIDs: " << numSto << endl;
#endif

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int lid=0;
  li_store_vbd = stoLIDVec[lid++];
  li_store_vbs = stoLIDVec[lid++];
  li_store_vgs = stoLIDVec[lid++];
  li_store_vds = stoLIDVec[lid++];
  li_store_von = stoLIDVec[lid++];

  if( loadLeadCurrent )
  {
    li_store_dev_id = stoLIDVec[lid++];
    li_store_dev_ig = stoLIDVec[lid++];
    li_store_dev_is = stoLIDVec[lid++];
    li_store_dev_ib = stoLIDVec[lid++];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  Store local indices:" << endl;
    cout << "  li_store_vbd           = " << li_store_vbd;
    cout << "  li_store_vbs           = " << li_store_vbs;
    cout << "  li_store_vgs           = " << li_store_vgs;
    cout << "  li_store_vds           = " << li_store_vds;
    cout << "  li_store_von           = " << li_store_von;
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
// Creator       : Robert Hoekstra, Computational Sciences
// Creation Date : 9/3/02
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  if( drainConductance != 0.0 && sourceConductance != 0.0 )
    return jacStamp_DC_SC;
  else if( drainConductance != 0.0 && sourceConductance == 0.0 )
    return jacStamp_DC;
  else if( drainConductance == 0.0 && sourceConductance != 0.0 )
    return jacStamp_SC;

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, Computational Sciences
// Creation Date : 9/3/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  vector<int> map;
  vector< vector<int> > map2;

  if (drainConductance != 0.0)
  {
    if (sourceConductance != 0.0)
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
    if (sourceConductance != 0.0)
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
  ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
  ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

  AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
  AGateEquBulkNodeOffset               = jacLIDVec[map[1]][map2[1][1]];
  AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][2]];
  AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][3]];

  ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
  ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

  ABulkEquGateNodeOffset               = jacLIDVec[map[3]][map2[3][0]];
  ABulkEquBulkNodeOffset               = jacLIDVec[map[3]][map2[3][1]];
  ABulkEquDrainPrimeNodeOffset         = jacLIDVec[map[3]][map2[3][2]];
  ABulkEquSourcePrimeNodeOffset        = jacLIDVec[map[3]][map2[3][3]];

  ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[4]][map2[4][0]];
  ADrainPrimeEquGateNodeOffset         = jacLIDVec[map[4]][map2[4][1]];
  ADrainPrimeEquBulkNodeOffset         = jacLIDVec[map[4]][map2[4][2]];
  ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[4]][map2[4][3]];
  ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[4]][map2[4][4]];

  ASourcePrimeEquGateNodeOffset        = jacLIDVec[map[5]][map2[5][0]];
  ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[5]][map2[5][1]];
  ASourcePrimeEquBulkNodeOffset        = jacLIDVec[map[5]][map2[5][2]];
  ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[5]][map2[5][3]];
  ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[5]][map2[5][4]];
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/06/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  // F-pointers:
  f_DrainEquDrainNodePtr             = 	  &(dFdx[li_Drain][ADrainEquDrainNodeOffset]);
  f_DrainEquDrainPrimeNodePtr        = 	  &(dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);

  f_GateEquGateNodePtr               = 	  &(dFdx[li_Gate][AGateEquGateNodeOffset]);
  f_GateEquBulkNodePtr               = 	  &(dFdx[li_Gate][AGateEquBulkNodeOffset]);
  f_GateEquDrainPrimeNodePtr         = 	  &(dFdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  f_GateEquSourcePrimeNodePtr        = 	  &(dFdx[li_Gate][AGateEquSourcePrimeNodeOffset]);

  f_SourceEquSourceNodePtr           = 	  &(dFdx[li_Source][ASourceEquSourceNodeOffset]);
  f_SourceEquSourcePrimeNodePtr      = 	  &(dFdx[li_Source][ASourceEquSourcePrimeNodeOffset]);

  f_BulkEquGateNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquGateNodeOffset]);
  f_BulkEquBulkNodePtr               = 	  &(dFdx[li_Bulk][ABulkEquBulkNodeOffset]);
  f_BulkEquDrainPrimeNodePtr         = 	  &(dFdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]);
  f_BulkEquSourcePrimeNodePtr        = 	  &(dFdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]);

  f_DrainPrimeEquDrainNodePtr        = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  f_DrainPrimeEquGateNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  f_DrainPrimeEquBulkNodePtr         = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]);
  f_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  f_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  f_SourcePrimeEquGateNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  f_SourcePrimeEquSourceNodePtr      = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  f_SourcePrimeEquBulkNodePtr        = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]);
  f_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  f_SourcePrimeEquSourcePrimeNodePtr = 	  &(dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);

  // Q-pointers:
  q_DrainEquDrainNodePtr             = 	  &(dQdx[li_Drain][ADrainEquDrainNodeOffset]);
  q_DrainEquDrainPrimeNodePtr        = 	  &(dQdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);

  q_GateEquGateNodePtr               = 	  &(dQdx[li_Gate][AGateEquGateNodeOffset]);
  q_GateEquBulkNodePtr               = 	  &(dQdx[li_Gate][AGateEquBulkNodeOffset]);
  q_GateEquDrainPrimeNodePtr         = 	  &(dQdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  q_GateEquSourcePrimeNodePtr        = 	  &(dQdx[li_Gate][AGateEquSourcePrimeNodeOffset]);

  q_SourceEquSourceNodePtr           = 	  &(dQdx[li_Source][ASourceEquSourceNodeOffset]);
  q_SourceEquSourcePrimeNodePtr      = 	  &(dQdx[li_Source][ASourceEquSourcePrimeNodeOffset]);

  q_BulkEquGateNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquGateNodeOffset]);
  q_BulkEquBulkNodePtr               = 	  &(dQdx[li_Bulk][ABulkEquBulkNodeOffset]);
  q_BulkEquDrainPrimeNodePtr         = 	  &(dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset]);
  q_BulkEquSourcePrimeNodePtr        = 	  &(dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset]);

  q_DrainPrimeEquDrainNodePtr        = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  q_DrainPrimeEquGateNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  q_DrainPrimeEquBulkNodePtr         = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset]);
  q_DrainPrimeEquDrainPrimeNodePtr   = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  q_DrainPrimeEquSourcePrimeNodePtr  = 	  &(dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  q_SourcePrimeEquGateNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  q_SourcePrimeEquSourceNodePtr      = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  q_SourcePrimeEquBulkNodePtr        = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset]);
  q_SourcePrimeEquDrainPrimeNodePtr  = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  q_SourcePrimeEquSourcePrimeNodePtr = 	  &(dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 diode instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
//  This is similar to the loadRHS function, only this function
//  only loads capacitor charges, and loads them into the daeQ vector.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/06/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double coef(0.0);

  double Qeqbs(0.0),Qeqbd(0.0),Qeqgb(0.0), Qeqgs(0.0), Qeqgd(0.0);
  //double ceqbs,ceqbd,ceqgb, ceqgs, ceqgd; // 3f5 vars
  int Dtype=model_.dtype;

  // do the same Dtype corrections on the charges that
  // are performed on the currents in the loadRHS function.

  // What is cbs and cbd?  They are diode currents (from exponentials),
  // so they are left out of this function.
  Qeqbs = Dtype*(qbs);
  Qeqbd = Dtype*(qbd);
  // These need "Dtype" here because we use them later *without*
  // Dtype, where SPICE uses it *with*
  Qeqgb = Dtype*(qgb);
  Qeqgs = Dtype*(qgs);
  Qeqgd = Dtype*(qgd);

  // 2 KCL for gate node
  coef = (Qeqgs+Qeqgd+Qeqgb);
  qVec[li_Gate] += coef*numberParallel;

  // 4 KCL for bulk node
  coef = Qeqbs + Qeqbd - Qeqgb;
  qVec[li_Bulk] += coef*numberParallel;

  // 5 KCL for drain' node
  coef = -(Qeqbd + Qeqgd);
  qVec[li_DrainPrime] += coef*numberParallel;

  // 6 KCL for source' node
  coef = -(Qeqbs + Qeqgs);
  qVec[li_SourcePrime] += coef*numberParallel;

  if( loadLeadCurrent )
  {
    double * storeLeadQ = extData.storeLeadCurrQCompRawPtr;
    storeLeadQ[li_store_dev_id] = (-(Qeqbd + Qeqgd))*numberParallel;
    storeLeadQ[li_store_dev_is] = (-(Qeqbs + Qeqgs))*numberParallel;
    storeLeadQ[li_store_dev_ig] = (Qeqgs+Qeqgd+Qeqgb)*numberParallel;
    storeLeadQ[li_store_dev_ib] = (Qeqbs + Qeqbd - Qeqgb)*numberParallel;
  }

  // Same as for the loadRHS function, but with capacitive terms:
  //  gcgd = Capgd;
  //  gcgs = Capgs;
  //  gcgb = Capgb;
  //  gcbs = capbs;
  //  gcbd = capbd;
  if(!origFlag)
  {
    // The setup of gcgd, etc. is the same as in the loadDAEdQdxVector
    // function.
    double gcgd, gcgs, gcgb, gcbs, gcbd;
    if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
    {
      gcgd = Capgd;
      gcgs = Capgs;
      gcgb = Capgb;
      // get at the two parasitic caps the same way
      gcbs = capbs;
      gcbd = capbd;
    }
    else
    {
      gcgd = 0.0; gcgs = 0.0; gcgb = 0.0; gcbs = 0.0; gcbd = 0.0;
    }

    // KCL 2
    double coef_Jdxp2 = Dtype*(gcgd*(vgd-vgd_orig)+gcgs*(vgs-vgs_orig)+
                        gcgb*(vgs-vgs_orig-vbs+vbs_orig));

    // 4 KCL for bulk node
    double coef_Jdxp4 = Dtype*(
                       - (gcgb)*(vgs-vgs_orig-vbs+vbs_orig)
                       + (gcgb)*(vbd-vbd_orig)
                       + (gcbs)*(vbs-vbs_orig));

    // 5 KCL for drain' node
    double coef_Jdxp5 = Dtype*(
                        -(gcgd)*(vgd-vgd_orig)
                        -(gcbd)*(vbd-vbd_orig));

    // 6 KCL for source' node
    double coef_Jdxp6 = Dtype*(-gcgs*(vgs-vgs_orig)-(gcbs)*(vbs-vbs_orig));

    double * dQdxdVp = extData.dQdxdVpVectorRawPtr;
    dQdxdVp[li_Gate       ] += coef_Jdxp2*numberParallel;
    dQdxdVp[li_Bulk       ] += coef_Jdxp4*numberParallel;
    dQdxdVp[li_DrainPrime ] += coef_Jdxp5*numberParallel;
    dQdxdVp[li_SourcePrime] += coef_Jdxp6*numberParallel;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 diode instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/06/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double coef(0.0);
  double ceqbs(0.0),ceqbd(0.0),ceqgb(0.0), ceqgs(0.0), ceqgd(0.0); // 3f5 vars
  double gmin1 = getDeviceOptions().gmin;

  int Dtype=model_.dtype;

  // The next few lines are the same as for loadRHS, except
  // the capcitor current terms have been set to zero here.
  ceqbs = Dtype*(cbs);
  ceqbd = Dtype*(cbd);
  // These need "Dtype" here because we use them later *without*
  // Dtype, where SPICE uses it *with*
  ceqgb = 0.0;
  ceqgs = 0.0;
  ceqgd = 0.0;

  // 1 KCL for drain node
  if (drainConductance != 0.0)
  {
    coef = Idrain;
    fVec[li_Drain] += coef*numberParallel;
  }

  // 2 KCL for gate node
  coef = (ceqgs+ceqgd+ceqgb);
  fVec[li_Gate] += coef*numberParallel;

  // 3 KCL for source node
  if (sourceConductance != 0.0)
  {
    coef = Isource;
    fVec[li_Source] += coef*numberParallel;
  }

  // 4 KCL for bulk node
  coef = ceqbs + ceqbd - ceqgb;
  fVec[li_Bulk] += coef*numberParallel;

  // 5 KCL for drain' node
  coef = -Idrain-(ceqbd - cdreq + ceqgd);
  fVec[li_DrainPrime] += coef*numberParallel;

  // 6 KCL for source' node
  coef = -Isource-(ceqbs + cdreq + ceqgs);
  fVec[li_SourcePrime] += coef*numberParallel;

  // Same as for the loadRHS function, but without capacitive terms:
  //  gcgd = Capgd;
  //  gcgs = Capgs;
  //  gcgb = Capgb;
  //  gcbs = capbs;
  //  gcbd = capbd;
  if (!origFlag)
  {
    // 4 KCL for bulk node
    double coef_Jdxp4 = Dtype*(
                       + ((gbd-gmin1))*(vbd-vbd_orig)
                       + ((gbs-gmin1))*(vbs-vbs_orig));

    // 5 KCL for drain' node
    double coef_Jdxp5 = Dtype*(
                        -((gbd-gmin1))*(vbd-vbd_orig)
                        +gds*(vds-vds_orig)
                        +Gm*((mode>0)?(vgs-vgs_orig):(vgd-vgd_orig))
                        +Gmbs*((mode>0)?(vbs-vbs_orig):(vbd-vbd_orig)));

    // 6 KCL for source' node
    double coef_Jdxp6 = Dtype*(
                        -((gbs-gmin1))*(vbs-vbs_orig)
                        -gds*(vds-vds_orig)
                        -Gm*((mode>0)?(vgs-vgs_orig):(vgd-vgd_orig))
                        -Gmbs*((mode>0)?(vbs-vbs_orig):(vbd-vbd_orig)));

    double * dFdxdVp = extData.dFdxdVpVectorRawPtr;
    dFdxdVp[li_Bulk       ] += coef_Jdxp4*numberParallel;
    dFdxdVp[li_DrainPrime ] += coef_Jdxp5*numberParallel;
    dFdxdVp[li_SourcePrime] += coef_Jdxp6*numberParallel;
  }

  if( loadLeadCurrent )
  {
    double * storeLeadF = extData.nextStoVectorRawPtr;
    storeLeadF[li_store_dev_id] = (-(ceqbd - cdreq + ceqgd))*numberParallel;
    storeLeadF[li_store_dev_is] = (-(ceqbs + cdreq + ceqgs))*numberParallel;
    storeLeadF[li_store_dev_ig] = (ceqgs+ceqgd+ceqgb)*numberParallel;
    storeLeadF[li_store_dev_ib] = (ceqbs + ceqbd - ceqgb)*numberParallel;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 diode instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/06/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  double gcgd(0.0);  // d(cqgd)/dVgd
  double gcgs(0.0);  // d(cqgs)/dVgs
  double gcgb(0.0);  // d(cqgb)/dVgb
  double gcbs(0.0);  // d(cqbs)/dVbs
  double gcbd(0.0);  // d(cqbd)/dVbd

  // get at the "conductances" for the gate capacitors with this trick
  //      gcgd = model_.dtype*Capgd;
  //      gcgs = model_.dtype*Capgs;
  //      gcgb = model_.dtype*Capgb;
  //
  //      In the loadRHS function, these would all be multiplied by
  //      getSolverState().pdt.  Here, for dQdx, the pdt term is left out.
  if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
  {
    gcgd = Capgd;
    gcgs = Capgs;
    gcgb = Capgb;
    // get at the two parasitic caps the same way
    gcbs = capbs;
    gcbd = capbd;
  }

  dQdx[li_Gate][AGateEquGateNodeOffset] += (gcgd+gcgs+gcgb)
    *numberParallel;
  dQdx[li_Gate][AGateEquBulkNodeOffset] -= gcgb*numberParallel;
  dQdx[li_Gate][AGateEquDrainPrimeNodeOffset] -= gcgd*numberParallel;
  dQdx[li_Gate][AGateEquSourcePrimeNodeOffset] -= gcgs*numberParallel;

  dQdx[li_Bulk][ABulkEquGateNodeOffset] -= gcgb*numberParallel;
  dQdx[li_Bulk][ABulkEquBulkNodeOffset] += +(gcbs+gcbd+gcgb)
    *numberParallel;
  dQdx[li_Bulk][ABulkEquDrainPrimeNodeOffset] -= +gcbd*numberParallel;
  dQdx[li_Bulk][ABulkEquSourcePrimeNodeOffset] -=
    +gcbs*numberParallel;

  dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset] +=
    -gcgd*numberParallel;
  dQdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset] +=
    -gcbd*numberParallel;
  dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset] +=
    +(gcbd+gcgd)*numberParallel;

  dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset] -=
    gcgs*numberParallel;
  dQdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset] -=
    +gcbs*numberParallel;
  dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]+=
    +(gcbs+gcgs)*numberParallel;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 diode instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/06/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Drain][ADrainEquDrainNodeOffset] +=
    drainConductance*numberParallel;
  dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset] -=
    drainConductance*numberParallel;

  dFdx[li_Source][ASourceEquSourceNodeOffset] +=
    sourceConductance*numberParallel;
  dFdx[li_Source][ASourceEquSourcePrimeNodeOffset] -=
    sourceConductance*numberParallel;

  dFdx[li_Bulk][ABulkEquBulkNodeOffset] +=
    (gbs+gbd)*numberParallel;
  dFdx[li_Bulk][ABulkEquDrainPrimeNodeOffset] -= gbd*numberParallel;
  dFdx[li_Bulk][ABulkEquSourcePrimeNodeOffset] -= gbs*numberParallel;

  dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset] -=
    drainConductance*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset] +=
    Gm*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquBulkNodeOffset] +=
    (-gbd+Gmbs)*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset] +=
    (drainConductance+gds+gbd+revsum)*numberParallel;
  dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset] +=
    (-gds-nrmsum)*numberParallel;

  dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset] -=
    Gm*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset] -=
    sourceConductance*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquBulkNodeOffset] -=
    (gbs+Gmbs)*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset] -=
    (gds+revsum)*numberParallel;
  dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset] +=
    (sourceConductance+gds+gbs+nrmsum)*numberParallel;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 03/01/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;
  // 3f5 likes to use the same variable names in local variables and in
  // structures.  Messes with us!  Define some local versions with capitals
  // instead
  double Von;
  double Vdsat;
  double Beta;
  //
  double evbs;
  double evbd;
  double sarg;
  double sargsw;
  //  double vgs1;
  //  double vgd1;
  //  double vgb1;
  double arg;
  int Check = 1;

  double capgs_old;
  double capgd_old;
  double capgb_old;

  // temporary storage for vgs/vgd/vds so that homotopy doesn't impact
  // voltage limiting "jdxp" terms

  double vgs_save;
  double vgd_save;
  double vds_save;

  // This is one of the vars that are set up at the top of mos3load that
  // should *not* be moved to the instance constructor!  tTransconductance
  // is set by updateTemperature.  Perhaps I should remove it from the
  // instance variables, too, since now it's pretty much a local thing.
  // same goes for the other things that depend on the t* variables!

  if( (tSatCurDens == 0) || (drainArea == 0) || (sourceArea == 0))
  {
    DrainSatCur = tSatCur;
    SourceSatCur = tSatCur;
  }
  else
  {
    DrainSatCur = tSatCurDens * drainArea;
    SourceSatCur = tSatCurDens * sourceArea;
  }
  Beta = tTransconductance * w/EffectiveLength;

  // don't do this anymore, use the one from the manager
  //  newtonIter = nlsMgrPtr->getNonLinearIter ();

  //  we need our solution variables for any of this stuff

  Vd = extData.nextSolVectorRawPtr[li_Drain];
  Vg = extData.nextSolVectorRawPtr[li_Gate];
  Vs = extData.nextSolVectorRawPtr[li_Source];
  Vb = extData.nextSolVectorRawPtr[li_Bulk];
  Vsp = extData.nextSolVectorRawPtr[li_SourcePrime];
  Vdp = extData.nextSolVectorRawPtr[li_DrainPrime];

  // now we need voltage drops
  Vddp = Vd - Vdp;
  Vssp = Vs - Vsp;
  Vbsp = Vb - Vsp;
  Vbdp = Vb - Vdp;
  Vgsp = Vg - Vsp;
  Vgdp = Vg - Vdp;
  Vgb  = Vg - Vb;
  Vdpsp = Vdp - Vsp;

  // Now the things that the 3f5 code really uses (from mos3load's
  // "general iteration" part  at lines 276-295
  vbs = model_.dtype * Vbsp;
  vgs = model_.dtype * Vgsp;
  vds = model_.dtype * Vdpsp;

  vbd = vbs-vds;
  vgd = vgs-vds;

  origFlag = 1;
  limitedFlag=false;
  vgs_orig = vgs;
  vds_orig = vds;
  vbs_orig = vbs;
  vbd_orig = vbd;
  vgd_orig = vgd;

  if (getSolverState().initJctFlag && !OFF && getDeviceOptions().voltageLimiterFlag)
  {
    if (IC_GIVEN)
    {
      vds = model_.dtype*icVDS;
      vgs = model_.dtype*icVGS;
      vbs = model_.dtype*icVBS;
      vbd = vbs - vds;
      vgd = vgs - vds;
      origFlag = false;
    }
    else
    {
      if (getSolverState().inputOPFlag)
      {
        N_LAS_Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
        if ((*flagSolVectorPtr)[li_Drain] == 0 || (*flagSolVectorPtr)[li_Gate] == 0 ||
            (*flagSolVectorPtr)[li_Source] == 0 || (*flagSolVectorPtr)[li_SourcePrime] ||
            (*flagSolVectorPtr)[li_DrainPrime] || (*flagSolVectorPtr)[li_Bulk] )
        {
          vbs = -1;
          vgs = model_.dtype*tVto;
          vds = 0;
          vbd = vbs-vds;
          vgd = vgs-vds;
        }
      }
      else
      {
        vbs = -1;
        vgs = model_.dtype*tVto;
        vds = 0;
        vbd = vbs-vds;
        vgd = vgs-vds;
      }
    }
  }
  else if ((getSolverState().initFixFlag || getSolverState().initJctFlag) && OFF)
  {
    vbs = vgs = vds = 0;
    vbd = vgd = 0;
  }
  
  if (getSolverState().newtonIter == 0)
  {

    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      vbs_old = (*extData.currStoVectorPtr)[li_store_vbs];
      vbd_old = (*extData.currStoVectorPtr)[li_store_vbd];
      vgs_old = (*extData.currStoVectorPtr)[li_store_vgs];
      vds_old = (*extData.currStoVectorPtr)[li_store_vds];
      Von = model_.dtype *
                (*extData.currStoVectorPtr)[li_store_von];
    }
    else
    { // otherwise there is no history
      vbs_old = vbs;
      vbd_old = vbd;
      vgs_old = vgs;
      vds_old = vds;
      Von = 0.0;
    }
    vgd_old = vgs_old-vds_old;
  }
  else
  {
    vbs_old = (*extData.nextStoVectorPtr)[li_store_vbs];
    vbd_old = (*extData.nextStoVectorPtr)[li_store_vbd];
    vgs_old = (*extData.nextStoVectorPtr)[li_store_vgs];
    vds_old = (*extData.nextStoVectorPtr)[li_store_vds];
    Von = model_.dtype *
              (*extData.nextStoVectorPtr)[li_store_von];
    vgd_old = vgs_old-vds_old;
  }

  ////////////////////////////////////////////
  // SPICE-type Voltage Limiting
  ////////////////////////////////////////////
  if (getDeviceOptions().voltageLimiterFlag)
  {
    // Do not do limiting if mode initfix and OFF:
    if (! (getSolverState().initFixFlag && OFF))
    {
      if (vds_old >= 0)
      {
        vgs = devSupport.fetlim( vgs, vgs_old, Von);
        vds = vgs - vgd;
        vds = devSupport.limvds( vds,  vds_old);
        vgd = vgs - vds;
      }
      else
      {
        vgd = devSupport.fetlim( vgd, vgd_old, Von);
        vds = vgs - vgd;
        vds = -devSupport.limvds( -vds, -vds_old );
        vgs = vgd + vds;
      }

      if (vds >= 0.0)
      {
        vbs = devSupport.pnjlim( vbs, vbs_old, vt, sourceVcrit, &Check);
        vbd = vbs - vds;
      }
      else
      {
        vbd = devSupport.pnjlim( vbd, vbd_old, vt, drainVcrit, &Check);
        vbs = vbd + vds;
      }

      // for convergence:
      if (Check == 1) limitedFlag=true;

    }
  }

  ////
  // now all the preliminaries are over - we can start doing the
  //  real work
  ////
  vbd = vbs - vds;
  vgd = vgs - vds;
  Vgb = vgs - vbs;

  // Now set the origFlag
  if (vgs_orig != vgs || vds_orig != vds ||
      vbs_orig != vbs || vbd_orig != vbd || vgd_orig != vgd) origFlag = 0;


  ////
  //  bulk-source and bulk-drain diodes
  //  here we just evaluate the ideal diode current and the
  //   corresponding derivative (conductance).
  ////
  if(vbs <= 0)
  {
    gbs = SourceSatCur/vt;
    gbs += getDeviceOptions().gmin;
    cbs = gbs*vbs;
  }
  else
  {
    evbs = exp(Xycemin(CONSTMAX_EXP_ARG,vbs/vt));
    gbs = (SourceSatCur*evbs/vt + getDeviceOptions().gmin);
    cbs = (SourceSatCur * (evbs-1) + getDeviceOptions().gmin*vbs);
  }
  if(vbd <= 0)
  {
    gbd = DrainSatCur/vt;
    gbd += getDeviceOptions().gmin;
    cbd = gbd *vbd;
  }
  else
  {
    evbd = exp(Xycemin(CONSTMAX_EXP_ARG,vbd/vt));
    gbd = (DrainSatCur*evbd/vt + getDeviceOptions().gmin);
    cbd = (DrainSatCur *(evbd-1) + getDeviceOptions().gmin*vbd);
  }

  // 3f5 does this simple stuff
  if (vds >= 0)
    mode = 1;
  else
    mode = -1;

  {
    // here 3f5 seems to have inserted a c-translation of a fortran routine
    // Sadly, there are GOTO's and labeled statements all over this thing.
    // I'll do what I can to eliminate them all, but I'll start by
    // leaving it VERBATIM and ugly
    //
    // evaluate the drain current, its derivatives and the
    // charges associated with the gate, channel and bulk for mosfets based on
    // semi-empirical equations

    ////
    // * subroutine moseq3(vds,vbs,vgs,gm,gds,gmbs,
    // * qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
    ////

    ////
    //     *     this routine evaluates the drain current, its derivatives and
    //     *     the charges associated with the gate, channel and bulk
    //     *     for mosfets based on semi-empirical equations
    ////

    //
    //      common /mosarg/ vto,beta,gamma,phi,phib,cox,xnsub,xnfs,xd,xj,xld,
    //      1   xlamda,uo,uexp,vbp,utra,vmax,xneff,xl,xw,vbi,von,vdsat,qspof,
    //      2   beta0,beta1,cdrain,xqco,xqc,fnarrw,fshort,lev
    //      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
    //      1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
    //      2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
    //      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
    //      1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
    //      2   pivtol,pivrel
    //

    //    // equivalence (xlamda,alpha),(vbp,theta),(uexp,eta),(utra,xkappa)

    double coeff0 = 0.0631353e0;
    double coeff1 = 0.8013292e0;
    double coeff2 = -0.01110777e0;
    double oneoverxl;  //  1/effective length
    double eta; // eta from model after length factor */
    double phibs;   // phi - vbs */
    double sqphbs;  // square root of phibs */
    double dsqdvb;  //  */
    double sqphis;  // square root of phi */
    double sqphs3;  // square root of phi cubed */
    double wps;
    double oneoverxj;   // 1/junction depth */
    double xjonxl;  // junction depth/effective length */
    double djonxj;
    double wponxj;
    double arga;
    double argb;
    double argc;
    double dwpdvb;
    double dadvb;
    double dbdvb;
    double gammas;
    double fbodys;
    double fbody;
    double onfbdy;
    double qbonco;
    double vbix;
    double wconxj;
    double dfsdvb;
    double dfbdvb;
    double dqbdvb;
    double vth;
    double dvtdvb;
    double csonco;
    double cdonco;
    double dxndvb;
    double dvodvb;
    double dvodvd;
    double vgsx;
    double dvtdvd;
    double onfg;
    double fgate;
    double us;
    double dfgdvg;
    double dfgdvd;
    double dfgdvb;
    double dvsdvg;
    double dvsdvb;
    double dvsdvd;
    double xn;
    double vdsc;
    double onvdsc;
    double dvsdga;
    double vdsx;
    double dcodvb;
    double cdnorm;
    double cdo;
    double cd1;
    double fdrain;
    double fd2;
    double dfddvg;
    double dfddvb;
    double dfddvd;
    double gdsat;
    double cdsat;
    double gdoncd;
    double gdonfd;
    double gdonfg;
    double dgdvg;
    double dgdvd;
    double dgdvb;
    double emax;
    double emongd;
    double demdvg;
    double demdvd;
    double demdvb;
    double delxl;
    double dldvd;
    double dldem;
    double ddldvg;
    double ddldvd;
    double ddldvb;
    double dlonxl;
    double xlfact;
    double diddl;
    double gds0;
    double emoncd;
    double ondvt;
    double onxn;
    double wfact;
    double gms;
    double gmw;
    double fshort;


    // Begin block of mosfet continuation code.
    // This idea is based, loosely, on a paper by Jaijeet
    // Roychowdhury.
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "HOMOTOPY INFO: gainscale   = "
	   << getSolverState().gainScale[blockHomotopyID] << endl;
      cout << "HOMOTOPY INFO: before vds  = " << vds << endl;
      cout << "HOMOTOPY INFO: before vgs  = " << vgs << endl;
    }
#endif

    // Save these before allowing homotopy to tweak them.  It
    // is important to restore them before moving on to
    // calculate RHS, because then the Jdxp terms will attempt to force
    // the external circuit to make these voltage drops the real thing!
    vds_save=vds;
    vgs_save=vgs;
    vgd_save=vgd;

    if (getSolverState().artParameterFlag)
    {
      double alpha = getSolverState().gainScale[blockHomotopyID];
      if (getDeviceOptions().staggerGainScale)
      {
        alpha *= (0.3 * randomPerturb + 1.0);
        if (alpha > 1.0)
        {
          alpha = 1.0;
        }
      }
      double vgstConst = getDeviceOptions().vgstConst;
      if (getDeviceOptions().randomizeVgstConst)
      {
        vgstConst *= randomPerturb;
      }

      vds = devSupport.contVds (vds, getSolverState().nltermScale, getDeviceOptions().vdsScaleMin);

      if (mode==1)
      {
        vgs = devSupport.contVgst (vgs, alpha, vgstConst);
      }
      else
      {
        vgd = devSupport.contVgst (vgd, alpha, vgstConst);
      }
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "HOMOTOPY INFO: after vds   = " << vds << endl;
      cout << "HOMOTOPY INFO: after vgs   = " << vgs << endl;
    }
#endif
    // End of block of mosfet continuation code.

    ////
    // bypasses the computation of charges
    ////

    ////
    // reference cdrain equations to source and
    // charge equations to bulk
    ////
    Vdsat = 0.0;
    oneoverxl = 1.0/EffectiveLength;

    // TVR: Massobrio and Antognetti call this "sigma"
    eta = model_.eta * 8.15e-22/(model_.oxideCapFactor*
                                       EffectiveLength*EffectiveLength*EffectiveLength);
    ////
    //.....square root term
    ////
    if ( (mode==1?vbs:vbd) <=  0.0 )
    {
      if (tPhi > 0)
        sqphis = sqrt(tPhi);
      else
        sqphis = 0;
      sqphs3 = tPhi*sqphis;
      phibs  =  tPhi-(mode==1?vbs:vbd);
      sqphbs  =  sqrt(phibs);
      dsqdvb  =  -0.5/sqphbs;
    }
    else
    {
      sqphis = sqrt(tPhi);
      sqphs3 = tPhi*sqphis;
      sqphbs = sqphis/(1.0+(mode==1?vbs:vbd)/
                       (tPhi+tPhi));
      phibs = sqphbs*sqphbs;
      dsqdvb = -phibs/(sqphs3+sqphs3);
    }
    ////
    //.....short channel effect factor
    ////
    if ( (model_.junctionDepth != 0.0) &&
         (model_.coeffDepLayWidth != 0.0) )
    {
      wps = model_.coeffDepLayWidth*sqphbs;
      oneoverxj = 1.0/model_.junctionDepth;
      xjonxl = model_.junctionDepth*oneoverxl;
      djonxj = model_.latDiff*oneoverxj;
      wponxj = wps*oneoverxj;
      wconxj = coeff0+coeff1*wponxj+coeff2*wponxj*wponxj;
      arga = wconxj+djonxj;
      argc = wponxj/(1.0+wponxj);
      argb = sqrt(1.0-argc*argc);
      fshort = 1.0-xjonxl*(arga*argb-djonxj);
      dwpdvb = model_.coeffDepLayWidth*dsqdvb;
      dadvb = (coeff1+coeff2*(wponxj+wponxj))*dwpdvb*oneoverxj;
      dbdvb = -argc*argc*(1.0-argc)*dwpdvb/(argb*wps);
      dfsdvb = -xjonxl*(dadvb*argb+arga*dbdvb);
    }
    else
    {
      fshort = 1.0;
      dfsdvb = 0.0;
    }
    ////
    //.....body effect
    //
    gammas = model_.gamma*fshort;
    fbodys = 0.5*gammas/(sqphbs+sqphbs);
    fbody = fbodys+model_.narrowFactor/w;

    onfbdy = 1.0/(1.0+fbody);
    dfbdvb = -fbodys*dsqdvb/sqphbs+fbodys*dfsdvb/fshort;
    qbonco =gammas*sqphbs+model_.narrowFactor*phibs/w;
    dqbdvb = gammas*dsqdvb+model_.gamma*dfsdvb*sqphbs-
      model_.narrowFactor/w;
    ////
    //.....static feedback effect
    //
    vbix = tVbi*model_.dtype-eta*(mode*vds);
    ////
    //.....threshold voltage
    //
    vth = vbix+qbonco;
    dvtdvd = -eta;
    dvtdvb = dqbdvb;
    ////
    //.....joint weak inversion and strong inversion
    //
    Von = vth;

    if ( model_.fastSurfaceStateDensity != 0.0 )
    {
      csonco = CONSTQ*model_.fastSurfaceStateDensity *
        1e4  *EffectiveLength*w/OxideCap;  // 1e4 to convert (cm**2/m**2)
      cdonco = qbonco/(phibs+phibs);
      xn = 1.0+csonco+cdonco;
      Von = vth+vt*xn;
      dxndvb = dqbdvb/(phibs+phibs)-qbonco*dsqdvb/(phibs*sqphbs);
      dvodvd = dvtdvd;
      dvodvb = dvtdvb+vt*dxndvb;
    }
    else
    {
      ////
      //.....cutoff region
      //

      if ( (mode==1?vgs:vgd) <= Von ) {
        cdrain = 0.0;
        gm = 0.0;
        gds = 0.0;
        gmbs = 0.0;
        goto innerline1000;
      }
    }
    ////
    //     *.....device is on
    ////

    vgsx = Xycemax((mode==1?vgs:vgd),Von);

    ////
    //.....mobility modulation by gate voltage
    ////
    onfg = 1.0+model_.theta*(vgsx-vth);
    fgate = 1.0/onfg;
    us = tSurfMob * 1e-4  *fgate;  // 1e4 to convert (m**2/cm**2)
    dfgdvg = -model_.theta*fgate*fgate;
    dfgdvd = -dfgdvg*dvtdvd;
    dfgdvb = -dfgdvg*dvtdvb;

    ////
    //     *.....saturation voltage
    ////
    Vdsat = (vgsx-vth)*onfbdy;

    if ( model_.maxDriftVel <= 0.0 )
    {
      dvsdvg = onfbdy;
      dvsdvd = -dvsdvg*dvtdvd;
      dvsdvb = -dvsdvg*dvtdvb-Vdsat*dfbdvb*onfbdy;
    }
    else
    {
      vdsc = EffectiveLength*model_.maxDriftVel/us;
      onvdsc = 1.0/vdsc;
      arga = (vgsx-vth)*onfbdy;
      argb = sqrt(arga*arga+vdsc*vdsc);
      Vdsat = arga+vdsc-argb;
      dvsdga = (1.0-arga/argb)*onfbdy;
      dvsdvg = dvsdga-(1.0-vdsc/argb)*vdsc*dfgdvg*onfg;
      dvsdvd = -dvsdvg*dvtdvd;
      dvsdvb = -dvsdvg*dvtdvb-arga*dvsdga*dfbdvb;
    }
    ////
    //     *.....current factors in linear region
    ////
    vdsx = Xycemin((mode*vds),Vdsat);
    if ( vdsx == 0.0 ) goto line900;

    cdo = vgsx-vth-0.5*(1.0+fbody)*vdsx;
    dcodvb = -dvtdvb-0.5*dfbdvb*vdsx;

    ////
    //     *.....normalized drain current
    ////
    cdnorm = cdo*vdsx;

    gm = vdsx;
    gds = vgsx-vth-(1.0+fbody+dvtdvd)*vdsx;
    gmbs = dcodvb*vdsx;
    ////
    //     *.....drain current without velocity saturation effect
    ////
    cd1 = Beta*cdnorm;
    Beta = Beta*fgate;
    cdrain = Beta*cdnorm;
    gm = Beta*gm+dfgdvg*cd1;
    gds = Beta*gds+dfgdvd*cd1;
    gmbs = Beta*gmbs;

    ////
    //     *.....velocity saturation factor
    ////
    if ( model_.maxDriftVel != 0.0 )
    {
      fdrain = 1.0/(1.0+vdsx*onvdsc);
      fd2 = fdrain*fdrain;
      arga = fd2*vdsx*onvdsc*onfg;
      dfddvg = -dfgdvg*arga;
      dfddvd = -dfgdvd*arga-fd2*onvdsc;
      dfddvb = -dfgdvb*arga;
      ////
      //       *.....drain current
      ////
      gm = fdrain*gm+dfddvg*cdrain;
      gds = fdrain*gds+dfddvd*cdrain;
      gmbs = fdrain*gmbs+dfddvb*cdrain;
      cdrain = fdrain*cdrain;
      Beta = Beta*fdrain;
    }

    ////
    //     *.....channel length modulation
    ////
    if ( (mode*vds) <= Vdsat ) goto line700;
    if ( model_.maxDriftVel <= 0.0 ) goto line510;
    if (model_.alpha == 0.0) goto line700;
    cdsat = cdrain;
    gdsat = cdsat*(1.0-fdrain)*onvdsc;
    gdsat = Xycemax(1.0e-12,gdsat);
    gdoncd = gdsat/cdsat;
    gdonfd = gdsat/(1.0-fdrain);
    gdonfg = gdsat*onfg;
    dgdvg = gdoncd*gm-gdonfd*dfddvg+gdonfg*dfgdvg;
    dgdvd = gdoncd*gds-gdonfd*dfddvd+gdonfg*dfgdvd;
    dgdvb = gdoncd*gmbs-gdonfd*dfddvb+gdonfg*dfgdvb;

    // to have this condition make sense, the option BADMOS3 must be
    // recognized on the options line.  A problem with the KAPPA parameter
    // was detected and fixed in 3f2.  BADMOS3 enables the unfixed deal,
    // just in case parameter fitting is affected.
    // It has never been implemented in Xyce
    //#ifdef BADMOS3_IMPLEMENTED
    //    if (ckt->CKTbadMos3)
    //      emax = cdsat*oneoverxl/gdsat;
    //    else
    //#endif
    emax = model_.kappa * cdsat*oneoverxl/gdsat;
    emoncd = emax/cdsat;
    emongd = emax/gdsat;
    demdvg = emoncd*gm-emongd*dgdvg;
    demdvd = emoncd*gds-emongd*dgdvd;
    demdvb = emoncd*gmbs-emongd*dgdvb;

    arga = 0.5*emax*model_.alpha;
    argc = model_.kappa*model_.alpha;
    argb = sqrt(arga*arga+argc*((mode*vds)-Vdsat));
    delxl = argb-arga;
    dldvd = argc/(argb+argb);
    dldem = 0.5*(arga/argb-1.0)*model_.alpha;
    ddldvg = dldem*demdvg;
    ddldvd = dldem*demdvd-dldvd;
    ddldvb = dldem*demdvb;
    goto line520;
  line510:
    delxl = sqrt(model_.kappa*((mode*vds)-Vdsat)*
                 model_.alpha);
    dldvd = 0.5*delxl/((mode*vds)-Vdsat);
    ddldvg = 0.0;
    ddldvd = -dldvd;
    ddldvb = 0.0;
    ////
    //     *.....punch through approximation
    ////
  line520:
    if ( delxl > (0.5*EffectiveLength) )
    {
      delxl = EffectiveLength-(EffectiveLength*EffectiveLength/
                               (4.0*delxl));
      arga = 4.0*(EffectiveLength-delxl)*(EffectiveLength-delxl)/
        (EffectiveLength*EffectiveLength);
      ddldvg = ddldvg*arga;
      ddldvd = ddldvd*arga;
      ddldvb = ddldvb*arga;
      dldvd =  dldvd*arga;
    }
    ////
    //     *.....saturation region
    ////
    dlonxl = delxl*oneoverxl;
    xlfact = 1.0/(1.0-dlonxl);
    cdrain = cdrain*xlfact;
    diddl = cdrain/(EffectiveLength-delxl);
    gm = gm*xlfact+diddl*ddldvg;
    gds0 = gds*xlfact+diddl*ddldvd;
    gmbs = gmbs*xlfact+diddl*ddldvb;
    gm = gm+gds0*dvsdvg;
    gmbs = gmbs+gds0*dvsdvb;
    gds = gds0*dvsdvd+diddl*dldvd;

    ////
    //     *.....finish strong inversion case
    ////
  line700:
    if ( (mode==1?vgs:vgd) < Von )
    {
      ////
      //       *.....weak inversion
      ////
      onxn = 1.0/xn;
      ondvt = onxn/vt;
      wfact = exp( ((mode==1?vgs:vgd)-Von)*ondvt );
      cdrain = cdrain*wfact;
      gms = gm*wfact;
      gmw = cdrain*ondvt;
      gm = gmw;

      if ((mode*vds) > Vdsat)
      {
        gm = gm+gds0*dvsdvg*wfact;
      }
      gds = gds*wfact+(gms-gmw)*dvodvd;
      gmbs = gmbs*wfact+(gms-gmw)*dvodvb-gmw*
        ((mode==1?vgs:vgd)-Von)*onxn*dxndvb;
    }
    ////
    //     *.....charge computation
    ////
    goto innerline1000;
    ////
    //     *.....special case of vds = 0.0d0
    ////
  line900:

    Beta = Beta*fgate;
    cdrain = 0.0;
    gm = 0.0;
    gds = Beta*(vgsx-vth);
    gmbs = 0.0;
    if ( (model_.fastSurfaceStateDensity != 0.0) &&
         ((mode==1?vgs:vgd) < Von) )
    {
      gds *=exp(((mode==1?vgs:vgd)-Von)/(vt*xn));
    }
  innerline1000:;
    ////
    //     *.....done
    ////
    // end of moseq3
    ////
  }

  // now deal with n vs p polarity

  von = model_.dtype * Von;
  vdsat = model_.dtype * Vdsat;

  ////
  //   *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
  ////

  cd = mode * cdrain - cbd;
  // in 3f5 this is all in a block conditioned on CKTmode, but since
  // it's valid for MODETRAN and MODETRANOP we'll just always do it

  ////
  // * now we do the hard part of the bulk-drain and bulk-source
  // * diode - we evaluate the non-linear capacitance and
  // * charge
  // *
  // * the basic equations are not hard, but the implementation
  // * is somewhat long in an attempt to avoid log/exponential
  // * evaluations
  ////
  ////
  // *  charge storage elements
  // *
  // *.. bulk-drain and bulk-source depletion capacitances
  ////
  // I took out all the CAPBYPASS stuff, and the
  // unnecessary curly braces that wind up there if you do

  // can't bypass the diode capacitance calculations
  if(Cbs != 0 || Cbssw != 0 )
  {
    if (vbs < tDepCap)
    {
      arg=1-vbs/tBulkPot;
      ////
      // * the following block looks somewhat long and messy,
      // * but since most users use the default grading
      // * coefficients of .5, and sqrt is MUCH faster than an
      // * exp(log()) we use this special case code to buy time.
      // * (as much as 10% of total job time!)
      ////
      if(model_.bulkJctBotGradingCoeff==model_.bulkJctSideGradingCoeff)
      {
        if(model_.bulkJctBotGradingCoeff == .5)
        {
          sarg = sargsw = 1/sqrt(arg);
        }
        else
        {
          sarg = sargsw = exp(-model_.bulkJctBotGradingCoeff*log(arg));
        }
      }
      else
      {
        if(model_.bulkJctBotGradingCoeff == .5)
        {
          sarg = 1/sqrt(arg);
        }
        else
        {
          sarg = exp(-model_.bulkJctBotGradingCoeff*log(arg));
        }
        if(model_.bulkJctSideGradingCoeff == .5)
        {
          sargsw = 1/sqrt(arg);
        }
        else
        {
          sargsw =exp(-model_.bulkJctSideGradingCoeff* log(arg));
        }
      }
      qbs = tBulkPot*(Cbs*(1-arg*sarg)/(1-model_.bulkJctBotGradingCoeff)
                      +Cbssw*(1-arg*sargsw)/(1-model_.bulkJctSideGradingCoeff));
      capbs=Cbs*sarg+ Cbssw*sargsw;

    }
    else
    {
      qbs = f4s + vbs*(f2s+vbs*(f3s/2));
      capbs=f2s+f3s*vbs;
    }
  }
  else
  {
    qbs = 0;
    capbs=0;
  }

  //// can't bypass the diode capacitance calculations
  if(Cbd != 0 || Cbdsw != 0 )
  {

    if (vbd < tDepCap)
    {
      arg=1-vbd/tBulkPot;
      ////
      // * the following block looks somewhat long and messy,
      // * but since most users use the default grading
      // * coefficients of .5, and sqrt is MUCH faster than an
      // * exp(log()) we use this special case code to buy time.
      // * (as much as 10% of total job time!)
      ////
      if(model_.bulkJctBotGradingCoeff == .5 &&
         model_.bulkJctSideGradingCoeff == .5)
      {
        sarg = sargsw = 1/sqrt(arg);
      }
      else
      {
        if(model_.bulkJctBotGradingCoeff == .5)
        {
          sarg = 1/sqrt(arg);
        }
        else
        {
          sarg = exp(-model_.bulkJctBotGradingCoeff*log(arg));
        }
        if(model_.bulkJctSideGradingCoeff == .5)
        {
          sargsw = 1/sqrt(arg);
        }
        else
        {
          sargsw =exp(-model_.bulkJctSideGradingCoeff*log(arg));
        }
      }
      qbd =
        tBulkPot*(Cbd*
                  (1-arg*sarg)
                  /(1-model_.bulkJctBotGradingCoeff)
                  +Cbdsw*
                  (1-arg*sargsw)
                  /(1-model_.bulkJctSideGradingCoeff));
      capbd=Cbd*sarg+
        Cbdsw*sargsw;
    }
    else
    {
      qbd = f4d +
        vbd * (f2d + vbd * f3d/2);
      capbd=f2d + vbd * f3d;
    }
  }
  else
  {
    qbd = 0;
    capbd = 0;
  }

  // Now  after a mess of convergence stuff that seems not to apply to us
  // 3f5 saves the vbs, vbd, etc. in the state vector (we don't, it gets
  // saved in updatePrimaryState)

  // Then 3f5 calculates meyer capacitances, capgs, capgb and capgd
  // Careful!  They use the local von and vdsat, which haven't got dtype
  // multiplying them!

  ////
  // *     calculate meyer's capacitors
  ////
  ////
  // * new cmeyer - this just evaluates at the current time,
  // * expects you to remember values from previous time
  // * returns 1/2 of non-constant portion of capacitance
  // * you must add in the other half from previous time
  // * and the constant part
  ////

  if (mode > 0)
  {
    devSupport.qmeyer (vgs,vgd,Vgb,Von,Vdsat,
                       capgs, capgd, capgb, tPhi,OxideCap);
  }
  else
  {
    devSupport.qmeyer (vgd,vgs,Vgb,Von,Vdsat,
                       capgd, capgs, capgb, tPhi,OxideCap);
  }


  //  vgs1 = vgs_old;
  //  vgd1 = vgs1 - vds_old;
  //   vgb1 = vgs1 - vbs_old;

  ////
  // TVR:  Note!  Caution with Capgs vs. capgs.  If I used the instance var
  // capgs on the left hand side, it's OK as long as we never try to implement
  // meyer back averaging, since capgs is going to be recalculated each
  // time through.  But if we do the averaging, the old one will have the
  // constant part and old old part added in, which is not what we wanted.
  // So I'll continue 3f5's way of having a local Capgs that's actually used
  // for the charge computation, and an instance capgs that's saved as state.
  ////

  if((getSolverState().dcopFlag))
//  if(!(getSolverState().tranopFlag))
  {
    Capgs =  2 * capgs + GateSourceOverlapCap ;
    Capgd =  2 * capgd + GateDrainOverlapCap ;
    Capgb =  2 * capgb + GateBulkOverlapCap ;
  }
  else
  {
    capgs_old = (*extData.currStaVectorPtr)[li_state_capgs];
    capgd_old = (*extData.currStaVectorPtr)[li_state_capgd];
    capgb_old = (*extData.currStaVectorPtr)[li_state_capgb];

    Capgs = ( capgs+
              capgs_old +
              GateSourceOverlapCap );
    Capgd = ( capgd+
              capgd_old +
              GateDrainOverlapCap );
    Capgb = ( capgb+
              capgb_old +
              GateBulkOverlapCap );
  }

  // Eric does this kludge in level 1, so I'll do it to make sure capacitances
  // are positive
  Capgs *= (Capgs < 0.0)?-1.0:1.0;
  Capgd *= (Capgd < 0.0)?-1.0:1.0;
  Capgb *= (Capgb < 0.0)?-1.0:1.0;


  // in the sequel, I'll refer to all derivatives of currents w.r.t. voltages
  // as conductances, whether they really are or not.

  // when we get here, cdrain has the channel current
  // cbd and cbs have the currents through the diodes
  // cd has the drain current minus the diode current and seems only to be used
  // for 3f5 convergence stuff
  // qbs and qbd have the charges on the parasitic capacitors
  // capbs and capbd have the capacitances of the parasitics.

  // none of the charges for the gate capacitors have been calculated yet.
  // We've saved the capacitances, so we can get the charges in
  // updatePrimaryState later.

  // Conductances:
  // gbd: the bulk-drain' conductance without the capacitor components
  //      We'll need to get the capacitor contribution in the actual load
  //      using C*dt
  // gbs: bulk-source' without capacitor
  // gm = derivative of channel current w.r.t. gate-source voltage --- gotta
  //     account for mode=normal or mode=reverse when using this!
  // gmbs = derivative of channel current w.r.t bulk-source voltage
  // gds = derivative of channel current w.r.t. drain-source voltage

  // the variables gcgs, gcgb, gcgd should be the conductances for the gate
  // capacitors, but we won't do those here (vide supra), we'll do them
  // in the jacobian load given the capacitances.

  // Now 3f5 doesn't do the resistor currents in the RHS load, because of
  // how they do their numerics.  We do, so let's save those here.

  Idrain = drainConductance * Vddp;
  Isource = sourceConductance * Vssp;

  if (mode >= 0)   // Normal mode
  {
    Gm = gm;      // (xnrm-xrev)*gm  in 3f5
    Gmbs = gmbs;  // (xnrm-xrev)*gmbs in 3f5
    nrmsum = Gm+Gmbs; // xnrm*(gm+gmbs)
    revsum = 0;       // xrev*(gm+gmbs)
    cdreq = model_.dtype*cdrain;
  }
  else
  {
    Gm = -gm;
    Gmbs = -gmbs;
    nrmsum = 0;
    revsum = -(Gm+Gmbs);  // because Gm and Gmbs already have - in them!
    cdreq = -(model_.dtype)*cdrain;
  }

  // It is now essential to restore the vds/vgs/vgd variables that might
  //   have been tweaked by homotopy, lest they have an effect on RHS
  // Jdxp terms.

  vds=vds_save;
  vgs=vgs_save;
  vgd=vgd_save;

  /// CURRENTS to load into RHS:

  // so at this point:

  // current out of drain is
  // Idrain

  // current out of gate:
  // dtype*( (deriv of qgs) + (deriv of qgd) + (deriv of qgb))

  //  the current *out of* the source should be simply
  // Isource.

  // current out of bulk is
  // dtype*(deriv of qbd) + dtype*cbd + dtype*cbs + dtype*(deriv of qbs)
  //  - dtype*(deriv of qgb)

  // current out of drain' is
  // -Idrain - dtype*(deriv of qgd) - (deriv of qbd) - dtype*cbd +
  //  mode*dtype*cdrain

  // the current out of the source' is
  //  -Isource - dtype*(deriv of qgs) - dtype*cbs - (deriv of qbs) -
  //   mode*dtype*cdrain

  //////Conductances to load into Jacobian as they relate to things here:
  /// all of the places where I say Cap/dt I really mean dtype*Cap/dt for
  // the meyer capacitors.  No dtype for the parasitics, though

  // 3f5 handles the mode by doing:
  //        where xnrm=1, xrev=0 if mode>=0, xnrm=0,xrev=1 if mode<0

  // drain-drain = a = drainConductance
  // drain-drain' = b = -drainConductance

  // gate-gate = c = "gcgd+gcgs+gcgb" = Capgd/dt+Capgs/dt+Capgb/dt
  // gate-bulk = d = -gcgb = -Capgb/dt
  // gate-drain' = e = -gcgd = -Capgd/dt
  // gate-source' = f = -gcgs = -Capgs/dt

  // source-source = g = sourceConductance
  // source-source' = h = -sourceConductance

  // bulk-gate = i = -gcgb = -Capgb/dt
  // bulk-bulk = j = gbs+gbd+gcgb+parasitics=gbs+gbd+Capgb/dt+capbs/dt+capbd/dt
  // bulk-drain' = k= -gbd-capbd/dt
  // bulk-source' = l= -gbs-capbs/dt

  // drain'-drain = m = -drainConductance
  // drain'-gate = n = (xnrm-xrev)*gm-Capgd/dt
  // drain'-bulk = o = -gbd-capbd/dt+(xnrm-xrev)*gmbs
  // drain'-drain' = p = drainConductance+gds+gbd+capbd/dt+
  //                   xrev*(gm+gmbs)+ Capgd/dt
  // drain'-source' = q = -gds-xnrm*(gm+gmbs)

  // source'-gate = r = -(xnrm-xrev)*gm-Capgs/dt
  // source'-source = s = -sourceConductance
  // source'-bulk = t = -gbs-capbs/dt-(xnrm-xrev)*gmbs
  // source'-drain' = u= -gds-xrev*(gm+gmbs)
  // source'-source' = v = sourceConductance+gds+gbs+capbs/dt+xnrm*(gm+gmbs)+
  //                       Capgs/dt

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/27/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  // mos3temp vars
  double czbd;    // zero voltage bulk-drain capacitance
  double czbdsw;  // zero voltage bulk-drain sidewall capacitance
  double czbs;    // zero voltage bulk-source capacitance
  double czbssw;  // zero voltage bulk-source sidewall capacitance
  double arg;     // 1 - fc
  double sarg;    // (1-fc) ^^ (-mj)
  double sargsw;  // (1-fc) ^^ (-mjsw)
  double ratio,ratio4;
  double fact2;
  double kt;
  double egfet;
  double pbfact;
  double capfact;
  double phio;
  double pbo;
  double gmanew,gmaold;
  // end of mos3temp stuff

  double tnom;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Instance::Begin of updateTemperature. \n";
    cout <<" name = " << getName() << endl;
    cout << endl;
  }
#endif

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) temp = temp_tmp;

  if (model_.interpolateTNOM(temp))
  {
    model_.processParams();
  }

  tnom = model_.tnom;
  ratio = temp/tnom;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Temperature = "<< temp << endl;
    cout << "tnom = " << tnom << endl;
    cout << "ratio = " << ratio << endl;
  }
#endif

  vt = temp * CONSTKoverQ;
  ratio = temp/tnom;
  fact2 = temp/CONSTREFTEMP;
  kt = temp * CONSTboltz;
  egfet = 1.16-(7.02e-4*temp*temp)/(temp+1108);
  arg = -egfet/(kt+kt)+1.1150877/(CONSTboltz*(CONSTREFTEMP+CONSTREFTEMP));
  pbfact = -2*vt *(1.5*log(fact2)+CONSTQ*arg);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "vt = " << vt << endl;
    cout << "ratio = " << ratio << endl;
    cout << "fact2 = " << fact2 << endl;
    cout << "kt = " << kt << endl;
    cout << "egfet = " << egfet << endl;
    cout << "arg = " << arg << endl;
    cout << "pbfact = " << pbfact << endl;
  }
#endif

  // in mos3temp 3f5 does a bunch of parameter defaulting (lines 155-163)
  //  but we assume that our various parameters have been set already in
  //  the constructors

  // lines 164-203 of mos3temp moved to instance block constructor

  // Here's the entire guts of the mos3temp instance loop, with obvious
  // modifications (here->MOS3 goes away, model->MOS3 turns into model_.)

  ratio4 = ratio * sqrt(ratio);
  tTransconductance = model_.transconductance / ratio4;
  tSurfMob = model_.surfaceMobility/ratio4;
  phio= (model_.phi-model_.pbfact1)/model_.fact1;
  tPhi = fact2 * phio + pbfact;
  tVbi = model_.vt0 - model_.dtype *
    (model_.gamma* sqrt(model_.phi))+.5*(model_.egfet1-egfet)
    + model_.dtype*.5* (tPhi-model_.phi);
  tVto = tVbi + model_.dtype * model_.gamma * sqrt(tPhi);
  tSatCur = model_.jctSatCur* exp(-egfet/vt+model_.egfet1/model_.vtnom);
  tSatCurDens = model_.jctSatCurDensity * exp(-egfet/vt+model_.egfet1/model_.vtnom);
  pbo = (model_.bulkJctPotential - model_.pbfact1)/model_.fact1;
  gmaold = (model_.bulkJctPotential-pbo)/pbo;
  capfact = 1/(1+model_.bulkJctBotGradingCoeff*
               (4e-4*(model_.tnom-CONSTREFTEMP)-gmaold));
  tCbd = model_.capBD * capfact;
  tCbs = model_.capBS * capfact;
  tCj = model_.bulkCapFactor * capfact;
  capfact = 1/(1+model_.bulkJctSideGradingCoeff*
               (4e-4*(model_.tnom-CONSTREFTEMP)-gmaold));
  tCjsw = model_.sideWallCapFactor * capfact;
  tBulkPot = fact2 * pbo+pbfact;
  gmanew = (tBulkPot-pbo)/pbo;
  capfact = (1+model_.bulkJctBotGradingCoeff*(4e-4*(temp-CONSTREFTEMP)-gmanew));
  tCbd *= capfact;
  tCbs *= capfact;
  tCj *= capfact;
  capfact = (1+model_.bulkJctSideGradingCoeff*(4e-4*(temp-CONSTREFTEMP)-gmanew));
  tCjsw *= capfact;
  tDepCap = model_.fwdCapDepCoeff * tBulkPot;

  if( (model_.jctSatCurDensity == 0) || (drainArea == 0) ||
      (sourceArea == 0) )
  {
    sourceVcrit = drainVcrit =
      vt*log(vt/(CONSTroot2*model_.jctSatCur));
  }
  else
  {
    drainVcrit = vt * log( vt / (CONSTroot2 *
                                 model_.jctSatCurDensity * drainArea));
    sourceVcrit = vt * log( vt / (CONSTroot2 *
                                  model_.jctSatCurDensity * sourceArea));
  }
  if(model_.capBDGiven)
  {
    czbd = tCbd;
  }
  else
  {
    if(model_.bulkCapFactorGiven)
    {
      czbd=tCj*drainArea;
    }
    else
    {
      czbd=0;
    }
  }
  if(model_.sideWallCapFactorGiven)
  {
    czbdsw= tCjsw * drainPerimeter;
  }
  else
  {
    czbdsw=0;
  }
  arg = 1-model_.fwdCapDepCoeff;
  sarg = exp( (-model_.bulkJctBotGradingCoeff) * log(arg) );
  sargsw = exp( (-model_.bulkJctSideGradingCoeff) * log(arg) );
  Cbd = czbd;
  Cbdsw = czbdsw;

  // The following lines were taken verbatim from the SPICE level 3 mosfet.
  // There is a mistake in f3d and f4d, which incorrectly use bulkJctPotential
  // instead of tBulkPot.  The result of this mistake is to have an incorrect
  // discontinuity in the charge calculations --- the derivative of charge
  // is right as the voltage drop approaches the discontinuity from either side
  // but when numerical differentiation is done to get currents one obtains a
  // massive contribution to the right hand side.  Since spice never checks
  // the RHS residual norm this seems not to have any effect on SPICE, but
  // Xyce pukes.
  //  f2d = czbd*(1-model_.fwdCapDepCoeff*
  //              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
  //    +  czbdsw*(1-model_.fwdCapDepCoeff*
  //               (1+model_.bulkJctSideGradingCoeff))*
  //    sargsw/arg;
  //  f3d = czbd * model_.bulkJctBotGradingCoeff * sarg/arg/
  //    model_.bulkJctPotential
  //    + czbdsw * model_.bulkJctSideGradingCoeff * sargsw/arg /
  //    model_.bulkJctPotential;
  //  f4d = czbd*model_.bulkJctPotential*(1-arg*sarg)/
  //    (1-model_.bulkJctBotGradingCoeff)
  //    + czbdsw*model_.bulkJctPotential*(1-arg*sargsw)/
  //    (1-model_.bulkJctSideGradingCoeff)
  //    -f3d/2*
  //    (tDepCap*tDepCap)
  //    -tDepCap * f2d;

  // These lines were taken from the equivalent section of the level 1 mosfet
  // and do not have the mistake in the lines above.
  f2d = czbd*(1-model_.fwdCapDepCoeff*
              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
    +  czbdsw*(1-model_.fwdCapDepCoeff*
               (1+model_.bulkJctSideGradingCoeff))*
    sargsw/arg;
  f3d = czbd * model_.bulkJctBotGradingCoeff * sarg/arg/
    tBulkPot
    + czbdsw * model_.bulkJctSideGradingCoeff * sargsw/arg /
    tBulkPot;
  f4d = czbd*tBulkPot*(1-arg*sarg)/
    (1-model_.bulkJctBotGradingCoeff)
    + czbdsw*tBulkPot*(1-arg*sargsw)/
    (1-model_.bulkJctSideGradingCoeff)
    -f3d/2*
    (tDepCap*tDepCap)
    -tDepCap * f2d;
  if(model_.capBSGiven)
  {
    czbs=tCbs;
  }
  else
  {
    if(model_.bulkCapFactorGiven)
    {
      czbs=tCj*sourceArea;
    }
    else
    {
      czbs=0;
    }
  }
  if(model_.sideWallCapFactorGiven)
  {
    czbssw = tCjsw * sourcePerimeter;
  }
  else
  {
    czbssw=0;
  }
  arg = 1-model_.fwdCapDepCoeff;
  sarg = exp( (-model_.bulkJctBotGradingCoeff) * log(arg) );
  sargsw = exp( (-model_.bulkJctSideGradingCoeff) * log(arg) );
  Cbs = czbs;
  Cbssw = czbssw;

  // See the comments above regarding f3d and f4d --- the same error
  // exists in the SPICE mosfet3 code
  //  f2s = czbs*(1-model_.fwdCapDepCoeff*
  //              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
  //    +  czbssw*(1-model_.fwdCapDepCoeff*
  //               (1+model_.bulkJctSideGradingCoeff))*
  //    sargsw/arg;
  //  f3s = czbs * model_.bulkJctBotGradingCoeff * sarg/arg/
  //    model_.bulkJctPotential
  //    + czbssw * model_.bulkJctSideGradingCoeff * sargsw/arg /
  //    model_.bulkJctPotential;
  //  f4s = czbs*model_.bulkJctPotential*(1-arg*sarg)/
  //    (1-model_.bulkJctBotGradingCoeff)
  //    + czbssw*model_.bulkJctPotential*(1-arg*sargsw)/
  //    (1-model_.bulkJctSideGradingCoeff)
  //    -f3s/2*
  //    (tBulkPot*tBulkPot)
  //    -tBulkPot * f2s;

  // so we use these lines from the MOSFET 1 instead.
  f2s = czbs*(1-model_.fwdCapDepCoeff*
              (1+model_.bulkJctBotGradingCoeff))* sarg/arg
    +  czbssw*(1-model_.fwdCapDepCoeff*
               (1+model_.bulkJctSideGradingCoeff))*
    sargsw/arg;
  f3s = czbs * model_.bulkJctBotGradingCoeff * sarg/arg/
    tBulkPot
    + czbssw * model_.bulkJctSideGradingCoeff * sargsw/arg /
    tBulkPot;
  f4s = czbs*tBulkPot*(1-arg*sarg)/
    (1-model_.bulkJctBotGradingCoeff)
    + czbssw*tBulkPot*(1-arg*sargsw)/
    (1-model_.bulkJctSideGradingCoeff)
    -f3s/2*
    (tDepCap*tDepCap)
    -tDepCap * f2s;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/28/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  double * staVector = extData.nextStaVectorRawPtr;
  double * oldstaVector = extData.currStaVectorRawPtr;
  double * stoVector = extData.nextStoVectorRawPtr;
  double * oldstoVector = extData.currStoVectorRawPtr;

  double vgs1, vgd1, vbs1,vgb1, vds1;

  bool bsuccess = updateIntermediateVars ();

  // voltage drops:
  stoVector[li_store_vbd] = vbd;
  stoVector[li_store_vbs] = vbs;
  stoVector[li_store_vgs] = vgs;
  stoVector[li_store_vds] = vds;
  stoVector[li_store_von] = von;

  // now the meyer capacitances
  // we didn't calculate these charges in update IntermediateVars
  // but we did calculate the voltage drops and capacitances.
  // first store the capacitances themselves:
  staVector[li_state_capgs] = capgs;
  staVector[li_state_capgd] = capgd;
  staVector[li_state_capgb] = capgb;

  // now the charges
  // BE CAREFUL!  We can only do Q=CV for DCOP!  Otherwise it's
  // supposed to be *INTEGRATED*:
  // Q = int(t0,t1)C(V)*dV --- and we approximate that by
  // Q(t1)-Q(t0) = CBar*(V(t1)-V(t0)) where CBar is the average.
  // Now with Meyer back averaging, Capxx is the average between the last
  // time step and this one.  So we gotta do the right thing for non-DCOP
  // when backaverage is on.

  if((getSolverState().dcopFlag))
//  if(!(getSolverState().tranopFlag))
  {
    qgs = Capgs*vgs;
    qgd = Capgd*vgd;
    qgb = Capgb*Vgb;
  }
  else
  {
    // get the ones from last time step
    qgs = oldstaVector[li_state_qgs];
    qgd = oldstaVector[li_state_qgd];
    qgb = oldstaVector[li_state_qgb];
    // get the voltage drops, too
    vgs1 = oldstoVector[li_store_vgs];
    vbs1 = oldstoVector[li_store_vbs];
    vds1 = oldstoVector[li_store_vds];

    vgb1 = vgs1-vbs1;
    vgd1 = vgs1-vds1;

    // NOW we can calculate the charge update
    qgs += Capgs*(vgs-vgs1);
    qgd += Capgd*(vgd-vgd1);
    qgb += Capgb*((vgs-vbs)-vgb1);
  }

  staVector[li_state_qgs] = qgs;
  staVector[li_state_qgd] = qgd;
  staVector[li_state_qgb] = qgb;

  // and the diode parasitic capacitors
  // these charges were set in updateIntermediateVars
  staVector[li_state_qbd] = qbd;
  staVector[li_state_qbs] = qbs;

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
  double wkfngs;
  double wkfng;
  double fermig;
  double fermis;
  double vfb;
  double kt1;
  double arg1;

  fact1 = tnom/CONSTREFTEMP;
  vtnom = tnom*CONSTKoverQ;
  kt1 = CONSTboltz * tnom;
  egfet1 = 1.16-(7.02e-4*tnom*tnom)/(tnom+1108);
  arg1 = -egfet1/(kt1+kt1)+1.1150877/(CONSTboltz*(CONSTREFTEMP+CONSTREFTEMP));
  pbfact1 = -2*vtnom *(1.5*log(fact1)+CONSTQ*arg1);

  if(oxideThickness == 0)
  {
    string msg = "Model:: ";
    msg += getName();
    msg +=" has TOX=0";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }
  else
  {
    oxideCapFactor = 3.9 * 8.854214871e-12/oxideThickness;
  }
  if(!given("U0") && !given("UO")) surfaceMobility=600;
  if(!given("KP"))
    transconductance = surfaceMobility * oxideCapFactor * 1e-4;
  if(given("NSUB"))
  {
    if(substrateDoping*1e6 >1.45e16)
    {
      if(!given("PHI"))
      {
        phi = 2*vtnom*
          log(substrateDoping*1e6/1.45e16);
        phi = Xycemax(0.1,phi);
      }
      fermis = dtype * .5 * phi;
      wkfng = 3.2;
      if(!given("TPG")) gateType=1;
      if(gateType != 0)
      {
        fermig = dtype *gateType*.5*egfet1;
        wkfng = 3.25 + .5 * egfet1 - fermig;
      }
      wkfngs = wkfng - (3.25 + .5 * egfet1 +fermis);
      if(!given("GAMMA"))
      {
        gamma = sqrt(2 * 11.70 * CONSTperm0 *
                     CONSTQ * substrateDoping*1e6)/
          oxideCapFactor;
      }
      if(!given("VTO"))
      {
        if(!given("NSS"))
          surfaceStateDensity=0;
        vfb = wkfngs - surfaceStateDensity*1e4*CONSTQ/oxideCapFactor;
        vt0 = vfb + dtype * (gamma * sqrt(phi)+ phi);
      }
      else
      {
        vfb = vt0 - dtype * (gamma*sqrt(phi)+phi);
      }
      alpha = ((11.7*CONSTperm0)+(11.7*CONSTperm0))/
        (CONSTQ*substrateDoping*1e6); //(cm**3/m**3)
      coeffDepLayWidth = sqrt(alpha);
    }
    else
    {
      string msg = "Model:: Nsub < Ni \n";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
  }

  // now model parameter preprocessing
  double mpi = M_PI;
  narrowFactor = delta * 0.5 * mpi * (11.7*CONSTperm0) /oxideCapFactor ;

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
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 2/26/01
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                        SolverState & ss1,
                                        DeviceOptions & do1)
  : DeviceModel(MB, ss1,do1),
    dtype(CONSTNMOS),
    tnom(getDeviceOptions().tnom),
    latDiff(0.0),
    jctSatCurDensity(0.0),
    jctSatCur(0.0),
    drainResistance(0.0),
    sourceResistance(0.0),
    sheetResistance(0.0),
    transconductance(0.0),
    gateSourceOverlapCapFactor(0.0),
    gateDrainOverlapCapFactor(0.0),
    gateBulkOverlapCapFactor(0.0),
    oxideCapFactor(0.0),
    vt0(0.0),
    capBD(0.0),
    capBS(0.0),
    bulkCapFactor(0.0),
    sideWallCapFactor(0.0),
    bulkJctPotential(0.0),
    bulkJctBotGradingCoeff(0.0),
    bulkJctSideGradingCoeff(0.0),
    fwdCapDepCoeff(0.0),
    phi(0.0),
    gamma(0.0),
    substrateDoping(0.0),
    gateType(0),
    surfaceStateDensity(0.0),
    oxideThickness(0.0),
    surfaceMobility(0.0),
    eta(0.0),
    junctionDepth(0.0),
    coeffDepLayWidth(0.0),
    narrowFactor(0.0),
    delta(0.0),
    fastSurfaceStateDensity(0.0),
    theta(0.0),
    maxDriftVel(0.0),
    alpha(0.0),
    kappa(0.0),
    fNcoef(0.0),
    fNexp(0.0),
    capBDGiven(0),
    capBSGiven(0),
    bulkCapFactorGiven(0),
    sideWallCapFactorGiven(0)
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
  if (!given("L"))
    model_l=getDeviceOptions().defl;
  if (!given("W"))
    model_w=getDeviceOptions().defw;
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;
  if (capBD != 0)
    capBDGiven = true;
  if (capBS != 0)
    capBSGiven = true;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  if (given("U0"))
  {
    string msg =  " ******************* \n";
    msg += ": WARNING: You have specified the surface mobility as u0 instead ";
    msg += "of uo.  This is supported, but ill-advised.\n";
    msg += " ***************** \n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
    if (given("UO"))
    {
      string msg = "Not only that, you specified both uo and u0, which is not allowed.";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
    surfaceMobility = surfaceMobility0;
  }

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
Model::~Model ()
{
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
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }
  os << endl;

  return os;
}

//-----------------------------------------------------------------------------
// MOSFET3 Master functions:
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

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & mi = *(*it);

    double * oldstaVectorPtr = mi.extData.currStaVectorRawPtr;
    double * stoVec = mi.extData.nextStoVectorRawPtr;
    double * oldstoVec = mi.extData.currStoVectorRawPtr;

    double vgs1(0.0), vgd1(0.0), vbs1(0.0),vgb1(0.0), vds1(0.0);

    bool btmp = mi.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    // voltage drops:
    stoVec[mi.li_store_vbd] = mi.vbd;
    stoVec[mi.li_store_vbs] = mi.vbs;
    stoVec[mi.li_store_vgs] = mi.vgs;
    stoVec[mi.li_store_vds] = mi.vds;
    stoVec[mi.li_store_von] = mi.von;

    // now the meyer capacitances
    // we didn't calculate these charges in update IntermediateVars
    // but we did calculate the voltage drops and capacitances.
    // first store the capacitances themselves:
    staVec[mi.li_state_capgs] = mi.capgs;
    staVec[mi.li_state_capgd] = mi.capgd;
    staVec[mi.li_state_capgb] = mi.capgb;

    // now the charges
    // BE CAREFUL!  We can only do Q=CV for DCOP!  Otherwise it's
    // supposed to be *INTEGRATED*:
    // Q = int(t0,t1)C(V)*dV --- and we approximate that by
    // Q(t1)-Q(t0) = CBar*(V(t1)-V(t0)) where CBar is the average.
    // Now with Meyer back averaging, Capxx is the average between the last
    // time step and this one.  So we gotta do the right thing for non-DCOP
    // when backaverage is on.

    if((getSolverState().dcopFlag))
//    if(!(getSolverState().tranopFlag))
    {
      mi.qgs = mi.Capgs*mi.vgs;
      mi.qgd = mi.Capgd*mi.vgd;
      mi.qgb = mi.Capgb*mi.Vgb;
    }
    else
    {
      // get the ones from last time step
      mi.qgs = (oldstaVectorPtr)[mi.li_state_qgs];
      mi.qgd = (oldstaVectorPtr)[mi.li_state_qgd];
      mi.qgb = (oldstaVectorPtr)[mi.li_state_qgb];
      // get the voltage drops, too
      vgs1 = oldstoVec[mi.li_store_vgs];
      vbs1 = oldstoVec[mi.li_store_vbs];
      vds1 = oldstoVec[mi.li_store_vds];

      vgb1 = vgs1-vbs1;
      vgd1 = vgs1-vds1;

      // NOW we can calculate the charge update
      mi.qgs += mi.Capgs*(mi.vgs-vgs1);
      mi.qgd += mi.Capgd*(mi.vgd-vgd1);
      mi.qgb += mi.Capgb*((mi.vgs-mi.vbs)-vgb1);
    }

    staVec[mi.li_state_qgs] = mi.qgs;
    staVec[mi.li_state_qgd] = mi.qgd;
    staVec[mi.li_state_qgb] = mi.qgb;

    // and the diode parasitic capacitors
    // these charges were set in updateIntermediateVars
    staVec[mi.li_state_qbd] = mi.qbd;
    staVec[mi.li_state_qbs] = mi.qbs;
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
  double gmin1 = getDeviceOptions().gmin;

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & mi = *(*it);

    int Dtype=mi.getModel().dtype;
    double ceqbs(0.0),ceqbd(0.0),ceqgb(0.0), ceqgs(0.0), ceqgd(0.0);
    double Qeqbs(0.0),Qeqbd(0.0),Qeqgb(0.0), Qeqgs(0.0), Qeqgd(0.0);
    double coef(0.0);

    // F-Vector:
    ceqbs = Dtype*(mi.cbs);
    ceqbd = Dtype*(mi.cbd);
    // These need "Dtype" here because we use them later *without*
    // Dtype, where SPICE uses it *with*
    ceqgb = 0.0;
    ceqgs = 0.0;
    ceqgd = 0.0;

    if (mi.drainConductance != 0.0)
    {
      fVec[mi.li_Drain] += mi.Idrain*mi.numberParallel;
    }

    coef = (ceqgs+ceqgd+ceqgb);

    fVec[mi.li_Gate] += coef*mi.numberParallel;

    if (mi.sourceConductance != 0.0)
    {
      fVec[mi.li_Source] += mi.Isource*mi.numberParallel;
    }

    coef = ceqbs + ceqbd - ceqgb;

    fVec[mi.li_Bulk] += coef*mi.numberParallel;

    coef = -mi.Idrain-(ceqbd - mi.cdreq + ceqgd);

    fVec[mi.li_DrainPrime] += coef*mi.numberParallel;

    coef = -mi.Isource-(ceqbs + mi.cdreq + ceqgs);

    fVec[mi.li_SourcePrime] += coef*mi.numberParallel;

    // Q-Vector:
    Qeqbs = Dtype*(mi.qbs);
    Qeqbd = Dtype*(mi.qbd);
    // These need "Dtype" here because we use them later *without*
    // Dtype, where SPICE uses it *with*
    Qeqgb = Dtype*(mi.qgb);
    Qeqgs = Dtype*(mi.qgs);
    Qeqgd = Dtype*(mi.qgd);

    coef = (Qeqgs+Qeqgd+Qeqgb);

    qVec[mi.li_Gate] += coef*mi.numberParallel;

    coef = Qeqbs + Qeqbd - Qeqgb;

    qVec[mi.li_Bulk] += coef*mi.numberParallel;

    coef = -(Qeqbd + Qeqgd);

    qVec[mi.li_DrainPrime] += coef*mi.numberParallel;

    coef = -(Qeqbs + Qeqgs);

    qVec[mi.li_SourcePrime] += coef*mi.numberParallel;

    // voltage limiters:
    if (!mi.origFlag)
    {
      // F-limiters:
      double coef_Jdxp4 = Dtype*(
            + ((mi.gbd-gmin1))*(mi.vbd-mi.vbd_orig)
            + ((mi.gbs-gmin1))*(mi.vbs-mi.vbs_orig));

      double coef_Jdxp5 = Dtype*(
            -((mi.gbd-gmin1))*(mi.vbd-mi.vbd_orig)
            +mi.gds*(mi.vds-mi.vds_orig)
            +mi.Gm*((mi.mode>0)?(mi.vgs-mi.vgs_orig):(mi.vgd-mi.vgd_orig))
            +mi.Gmbs*((mi.mode>0)?(mi.vbs-mi.vbs_orig):(mi.vbd-mi.vbd_orig)));

      double coef_Jdxp6 = Dtype*(
            -((mi.gbs-gmin1))*(mi.vbs-mi.vbs_orig)
            -mi.gds*(mi.vds-mi.vds_orig)
            -mi.Gm*((mi.mode>0)?(mi.vgs-mi.vgs_orig):(mi.vgd-mi.vgd_orig))
            -mi.Gmbs*((mi.mode>0)?(mi.vbs-mi.vbs_orig):(mi.vbd-mi.vbd_orig)));

      double * dFdxdVp = mi.extData.dFdxdVpVectorRawPtr;
      dFdxdVp[mi.li_Bulk       ] += coef_Jdxp4*mi.numberParallel;
      dFdxdVp[mi.li_DrainPrime ] += coef_Jdxp5*mi.numberParallel;
      dFdxdVp[mi.li_SourcePrime] += coef_Jdxp6*mi.numberParallel;

      // Q-limiters:
      {
        double gcgd(0.0), gcgs(0.0), gcgb(0.0), gcbs(0.0), gcbd(0.0);
        if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
        {
          gcgd = mi.Capgd;
          gcgs = mi.Capgs;
          gcgb = mi.Capgb;
          // get at the two parasitic caps the same way
          gcbs = mi.capbs;
          gcbd = mi.capbd;
        }
        else
        {
          gcgd = 0.0; gcgs = 0.0; gcgb = 0.0; gcbs = 0.0; gcbd = 0.0;
        }

        double coef_Jdxp2 =
          Dtype*(gcgd*(mi.vgd-mi.vgd_orig)+gcgs*(mi.vgs-mi.vgs_orig)+
          gcgb*(mi.vgs-mi.vgs_orig-mi.vbs+mi.vbs_orig));

        double coef_Jdxp4 = Dtype*(
            - (gcgb)*(mi.vgs-mi.vgs_orig-mi.vbs+mi.vbs_orig)
            + (gcgb)*(mi.vbd-mi.vbd_orig)
            + (gcbs)*(mi.vbs-mi.vbs_orig));

        double coef_Jdxp5 = Dtype*(
            -(gcgd)*(mi.vgd-mi.vgd_orig)
            -(gcbd)*(mi.vbd-mi.vbd_orig));

        // 6 KCL for source' node
        double coef_Jdxp6 = Dtype*
          (-gcgs*(mi.vgs-mi.vgs_orig)-(gcbs)*(mi.vbs-mi.vbs_orig));

        double * dQdxdVp = mi.extData.dQdxdVpVectorRawPtr;
        dQdxdVp[mi.li_Gate       ] += coef_Jdxp2*mi.numberParallel;
        dQdxdVp[mi.li_Bulk       ] += coef_Jdxp4*mi.numberParallel;
        dQdxdVp[mi.li_DrainPrime ] += coef_Jdxp5*mi.numberParallel;
        dQdxdVp[mi.li_SourcePrime] += coef_Jdxp6*mi.numberParallel;
      }
    }

    if( mi.loadLeadCurrent )
    {
      storeLeadF[mi.li_store_dev_id] = (-(ceqbd - mi.cdreq + ceqgd))*mi.numberParallel;
      storeLeadQ[mi.li_store_dev_id] = (-(Qeqbd + Qeqgd))*mi.numberParallel;
      storeLeadF[mi.li_store_dev_is] = (-(ceqbs + mi.cdreq + ceqgs))*mi.numberParallel;
      storeLeadQ[mi.li_store_dev_is] = (-(Qeqbs + Qeqgs))*mi.numberParallel;
      storeLeadF[mi.li_store_dev_ig] = (ceqgs+ceqgd+ceqgb)*mi.numberParallel;
      storeLeadQ[mi.li_store_dev_ig] = (Qeqgs+Qeqgd+Qeqgb)*mi.numberParallel;
      storeLeadF[mi.li_store_dev_ib] = (ceqbs + ceqbd - ceqgb)*mi.numberParallel;
      storeLeadQ[mi.li_store_dev_ib] = (Qeqbs + Qeqbd - Qeqgb)*mi.numberParallel;
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
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & mi = *(*it);

    // F-matrix:

    *mi.f_DrainEquDrainNodePtr +=
    mi.drainConductance*mi.numberParallel;

    *mi.f_DrainEquDrainPrimeNodePtr -=
    mi.drainConductance*mi.numberParallel;


    *mi.f_SourceEquSourceNodePtr +=
    mi.sourceConductance*mi.numberParallel;

    *mi.f_SourceEquSourcePrimeNodePtr -=
    mi.sourceConductance*mi.numberParallel;


    *mi.f_BulkEquBulkNodePtr +=
    (mi.gbs+mi.gbd)*mi.numberParallel;

    *mi.f_BulkEquDrainPrimeNodePtr -= mi.gbd*mi.numberParallel;

    *mi.f_BulkEquSourcePrimeNodePtr -= mi.gbs*mi.numberParallel;


    *mi.f_DrainPrimeEquDrainNodePtr -=
    mi.drainConductance*mi.numberParallel;

    *mi.f_DrainPrimeEquGateNodePtr +=
    mi.Gm*mi.numberParallel;

    *mi.f_DrainPrimeEquBulkNodePtr +=
    (-mi.gbd+mi.Gmbs)*mi.numberParallel;

    *mi.f_DrainPrimeEquDrainPrimeNodePtr +=
    (mi.drainConductance+mi.gds+mi.gbd+mi.revsum)*mi.numberParallel;

    *mi.f_DrainPrimeEquSourcePrimeNodePtr +=
    (-mi.gds-mi.nrmsum)*mi.numberParallel;


    *mi.f_SourcePrimeEquGateNodePtr -=
    mi.Gm*mi.numberParallel;

    *mi.f_SourcePrimeEquSourceNodePtr -=
    mi.sourceConductance*mi.numberParallel;

    *mi.f_SourcePrimeEquBulkNodePtr -=
    (mi.gbs+mi.Gmbs)*mi.numberParallel;

    *mi.f_SourcePrimeEquDrainPrimeNodePtr -=
    (mi.gds+mi.revsum)*mi.numberParallel;

    *mi.f_SourcePrimeEquSourcePrimeNodePtr +=
    (mi.sourceConductance+mi.gds+mi.gbs+mi.nrmsum)*mi.numberParallel;

    // Q-matrix:
    double gcgd(0.0);  // d(cqgd)/dVgd
    double gcgs(0.0);  // d(cqgs)/dVgs
    double gcgb(0.0);  // d(cqgb)/dVgb
    double gcbs(0.0);  // d(cqbs)/dVbs
    double gcbd(0.0);  // d(cqbd)/dVbd

    // get at the "conductances" for the gate capacitors with this trick
    //      gcgd = model_.dtype*Capgd;
    //      gcgs = model_.dtype*Capgs;
    //      gcgb = model_.dtype*Capgb;
    //
    //      In the loadRHS function, these would all be multiplied by
    //      getSolverState().pdt.  Here, for *mi.q_, the pdt term is left out.
    if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
    {
      gcgd = mi.Capgd;
      gcgs = mi.Capgs;
      gcgb = mi.Capgb;
      // get at the two parasitic caps the same way
      gcbs = mi.capbs;
      gcbd = mi.capbd;
    }


    *mi.q_GateEquGateNodePtr +=
    (gcgd+gcgs+gcgb)*mi.numberParallel;

    *mi.q_GateEquBulkNodePtr -= gcgb*mi.numberParallel;

    *mi.q_GateEquDrainPrimeNodePtr -= gcgd*mi.numberParallel;

    *mi.q_GateEquSourcePrimeNodePtr -= gcgs*mi.numberParallel;


    *mi.q_BulkEquGateNodePtr -= gcgb*mi.numberParallel;

    *mi.q_BulkEquBulkNodePtr +=
    (+gcbs+gcbd+gcgb)*mi.numberParallel;

    *mi.q_BulkEquDrainPrimeNodePtr -= +gcbd*mi.numberParallel;

    *mi.q_BulkEquSourcePrimeNodePtr -=
    +gcbs*mi.numberParallel;


    *mi.q_DrainPrimeEquGateNodePtr +=
    -gcgd*mi.numberParallel;

    *mi.q_DrainPrimeEquBulkNodePtr +=
    -gcbd*mi.numberParallel;

    *mi.q_DrainPrimeEquDrainPrimeNodePtr +=
    (+gcbd+gcgd)*mi.numberParallel;


    *mi.q_SourcePrimeEquGateNodePtr -=
    gcgs*mi.numberParallel;

    *mi.q_SourcePrimeEquBulkNodePtr -=
    +gcbs*mi.numberParallel;

    *mi.q_SourcePrimeEquSourcePrimeNodePtr+=
    (+gcbs+gcgs)*mi.numberParallel;
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
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  int sizeInstances = instanceContainer_.size();
  for (int i=0; i<sizeInstances; ++i)
  {
    Instance & mi = *(instanceContainer_.at(i));

    // F-matrix:
    dFdx[mi.li_Drain][mi.ADrainEquDrainNodeOffset] +=
    mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_Drain][mi.ADrainEquDrainPrimeNodeOffset] -=
    mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_Source][mi.ASourceEquSourceNodeOffset] +=
    mi.sourceConductance*mi.numberParallel;

    dFdx[mi.li_Source][mi.ASourceEquSourcePrimeNodeOffset] -=
    mi.sourceConductance*mi.numberParallel;

    dFdx[mi.li_Bulk][mi.ABulkEquBulkNodeOffset] +=
    (mi.gbs+mi.gbd)*mi.numberParallel;

    dFdx[mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset] -= mi.gbd*mi.numberParallel;
    dFdx[mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset] -= mi.gbs*mi.numberParallel;


    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainNodeOffset] -=
    mi.drainConductance*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquGateNodeOffset] +=
    mi.Gm*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset] +=
    (-mi.gbd+mi.Gmbs)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset] +=
    (mi.drainConductance+mi.gds+mi.gbd+mi.revsum)*mi.numberParallel;

    dFdx[mi.li_DrainPrime][mi.ADrainPrimeEquSourcePrimeNodeOffset] +=
    (-mi.gds-mi.nrmsum)*mi.numberParallel;


    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquGateNodeOffset] -=
    mi.Gm*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourceNodeOffset] -=
    mi.sourceConductance*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset] -=
    (mi.gbs+mi.Gmbs)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquDrainPrimeNodeOffset] -=
    (mi.gds+mi.revsum)*mi.numberParallel;

    dFdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset] +=
    (mi.sourceConductance+mi.gds+mi.gbs+mi.nrmsum)*mi.numberParallel;

    // Q-matrix:
    double gcgd(0.0);  // d(cqgd)/dVgd
    double gcgs(0.0);  // d(cqgs)/dVgs
    double gcgb(0.0);  // d(cqgb)/dVgb
    double gcbs(0.0);  // d(cqbs)/dVbs
    double gcbd(0.0);  // d(cqbd)/dVbd

    // get at the "conductances" for the gate capacitors with this trick
    //      gcgd = model_.dtype*Capgd;
    //      gcgs = model_.dtype*Capgs;
    //      gcgb = model_.dtype*Capgb;
    //
    //      In the loadRHS function, these would all be multiplied by
    //      getSolverState().pdt.  Here, for dQdx, the pdt term is left out.
    if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
    {
      gcgd = mi.Capgd;
      gcgs = mi.Capgs;
      gcgb = mi.Capgb;
      // get at the two parasitic caps the same way
      gcbs = mi.capbs;
      gcbd = mi.capbd;
    }


    dQdx[mi.li_Gate][mi.AGateEquGateNodeOffset] +=
    (gcgd+gcgs+gcgb)*mi.numberParallel;

    dQdx[mi.li_Gate][mi.AGateEquBulkNodeOffset] -= gcgb*mi.numberParallel;
    dQdx[mi.li_Gate][mi.AGateEquDrainPrimeNodeOffset] -= gcgd*mi.numberParallel;
    dQdx[mi.li_Gate][mi.AGateEquSourcePrimeNodeOffset] -= gcgs*mi.numberParallel;
    dQdx[mi.li_Bulk][mi.ABulkEquGateNodeOffset] -= gcgb*mi.numberParallel;
    dQdx[mi.li_Bulk][mi.ABulkEquBulkNodeOffset] +=
    (+gcbs+gcbd+gcgb)*mi.numberParallel;
    dQdx[mi.li_Bulk][mi.ABulkEquDrainPrimeNodeOffset] -= +gcbd*mi.numberParallel;
    dQdx[mi.li_Bulk][mi.ABulkEquSourcePrimeNodeOffset] -=
    +gcbs*mi.numberParallel;
    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquGateNodeOffset] +=
    -gcgd*mi.numberParallel;
    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquBulkNodeOffset] +=
    -gcbd*mi.numberParallel;
    dQdx[mi.li_DrainPrime][mi.ADrainPrimeEquDrainPrimeNodeOffset] +=
    (+gcbd+gcgd)*mi.numberParallel;
    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquGateNodeOffset] -=
    gcgs*mi.numberParallel;
    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquBulkNodeOffset] -=
    +gcbs*mi.numberParallel;
    dQdx[mi.li_SourcePrime][mi.ASourcePrimeEquSourcePrimeNodeOffset]+=
    (+gcbs+gcgs)*mi.numberParallel;
  }
  return true;
}
#endif

} // namespace MOSFET3
} // namespace Device
} // namespace Xyce
