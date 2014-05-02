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
//
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DiodePDEInstance.C,v $
//
// Purpose        : One dimensional PDE device, instance class
//                  implementation.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/06/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.117.2.2 $
//
// Revision Date  : $Date: 2014/03/17 17:04:42 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <N_UTL_Misc.h>
#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#include <N_DEV_fwd.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DiodePDE.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_RegionData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

namespace Xyce {
namespace Device {
namespace DiodePDE {

// default number of mesh points:
static const int NUM_MESH_POINTS = 11;

// default maximum number of nonzero entries in a matrix row
static const int MAX_COLS_PER_ROW = 40;


void Traits::loadInstanceParameters(ParametricData<DiodePDE::Instance> &p)
{
  p.addPar("ANODE.BC", 0.5, &DiodePDE::Instance::anodebc);
  p.addPar("CATHODE.BC", 0.0, &DiodePDE::Instance::cathodebc);
  p.addPar("EMITTER.BC", 0.5, &DiodePDE::Instance::emitterbc);
  p.addPar("COLLECTOR.BC", 0.0, &DiodePDE::Instance::collectorbc);
  p.addPar("BASE.BC", 0.0, &DiodePDE::Instance::basebc);
  p.addPar("BASE.LOC", 0.5e-3, &DiodePDE::Instance::baseLocation)
    .setGivenMember(&DiodePDE::Instance::baseLocationGiven);
  p.addPar("AREA", 1.0, &DiodePDE::Instance::area);

  // user-specified scaling vars:
  p.addPar("X0", 1.0e-7, &DiodePDE::Instance::x0_user)
    .setExpressionAccess(ParameterType::TIME_DEP);
  p.addPar("C0", 1.0e+15, &DiodePDE::Instance::C0_user)
    .setExpressionAccess(ParameterType::TIME_DEP);
  p.addPar("t0", 1.0e-6, &DiodePDE::Instance::t0_user)
    .setExpressionAccess(ParameterType::TIME_DEP);
  p.addPar("SCALEDENSITYTOMAXDOPING", true, &DiodePDE::Instance::scaleDensityToMaxDoping_)
    .setUnit(U_LOGIC)
    .setDescription("If set the density will be scaled by a fraction of the maximum doping.");
  p.addPar("DENSITYSCALARFRACTION", 1.0e-1, &DiodePDE::Instance::densityScalarFraction_)
    .setUnit(U_LOGIC)
    .setDescription("Fraction of the maximum doping by which density will be scaled.");
  p.addPar("NA", 1.0e+15, &DiodePDE::Instance::Na);
  p.addPar("ND", 1.0e+15, &DiodePDE::Instance::Nd);
  p.addPar("WJ", 1.0e-4, &DiodePDE::Instance::WJ);
  p.addPar("TEMP", 300.15, &DiodePDE::Instance::Temp);
  p.addPar("ANODE.AREA", 0.0, &DiodePDE::Instance::anodeArea);
  p.addPar("EMITTER.AREA", 0.0, &DiodePDE::Instance::emitterArea);
  p.addPar("CATHODE.AREA", 0.0, &DiodePDE::Instance::cathodeArea);
  p.addPar("COLLECTOR.AREA", 0.0, &DiodePDE::Instance::collectorArea);
  p.addPar("BASE.AREA", 0.0, &DiodePDE::Instance::baseArea);
  p.addPar("L", 1.0e-3, &DiodePDE::Instance::length);
  p.addPar("W", 1.0e-3, &DiodePDE::Instance::width);
  p.addPar("MAXVOLTDELTA", 0.025, &DiodePDE::Instance::maxVoltDelta);
  p.addPar("OUTPUTINTERVAL", 0.0, &DiodePDE::Instance::outputInterval)
    .setGivenMember(&DiodePDE::Instance::outputIntervalGiven);
  p.addPar("basex", 0.5e-3, &DiodePDE::Instance::basex);

#ifdef Xyce_DEBUG_DEVICE
  p.addPar("ANODEINDEX", 1, &DiodePDE::Instance::anodeIndex_user)
    .setGivenMember(&DiodePDE::Instance::anodeIndex_userGiven);
  p.addPar("CATHODEINDEX", 0, &DiodePDE::Instance::cathodeIndex_user)
    .setGivenMember(&DiodePDE::Instance::cathodeIndex_userGiven);
#endif

  // Neutron stuff:
  p.addPar("JUNCTIONAREA", 1.0e-4, &DiodePDE::Instance::junctionArea);

  // Set up map for non-double precision variables:
  p.addPar("GRADED", false, &DiodePDE::Instance::gradedJunctionFlag);
  p.addPar("BJTENABLE", false, &DiodePDE::Instance::bjtEnableFlag);
  p.addPar("MOBMODEL", std::string("ARORA"), &DiodePDE::Instance::mobModelName);
  p.addPar("FIELDDEP", false, &DiodePDE::Instance::fieldDependentMobility)
    .setGivenMember(&DiodePDE::Instance::fieldDependentMobilityGiven)
    .setUnit(U_LOGIC)
    .setDescription("If true, use field dependent mobility.");
  p.addPar("BULKMATERIAL", std::string("SI"), &DiodePDE::Instance::bulkMaterial);
  p.addPar("OFFSETOUTPUTVOLTAGE", true, &DiodePDE::Instance::useVoltageOutputOffset_);
  p.addPar("FIRSTELECTRODEOFFSET", false, &DiodePDE::Instance::offsetWithFirstElectrode_);
  p.addPar("DISPLCUR", false, &DiodePDE::Instance::displCurrentFlag);
  p.addPar("OUTPUTNLPOISSON", false, &DiodePDE::Instance::outputNLPoisson);
  p.addPar("OUTPUTREGION", 1, &DiodePDE::Instance::outputRegion);
  p.addPar("MASKVARSTIA", false, &DiodePDE::Instance::maskVarsTIAFlag_)
    .setUnit(U_LOGIC)
    .setDescription("If set to true, then some variables are excluded from the time integration error control calculation.");
  p.addPar("AUGER", true, &DiodePDE::Instance::includeAugerRecomb)
    .setUnit(U_LOGIC);
  p.addPar("SRH", true, &DiodePDE::Instance::includeSRHRecomb)
    .setUnit(U_LOGIC);
  p.addPar("GNUPLOTLEVEL", 1, &DiodePDE::Instance::gnuplotLevel);
  p.addPar("TECPLOTLEVEL", 1, &DiodePDE::Instance::tecplotLevel);
  p.addPar("SGPLOTLEVEL", 0, &DiodePDE::Instance::sgplotLevel);
  p.addPar("VOLTLIM", false, &DiodePDE::Instance::voltLimFlag);
  p.addPar("MESHFILE", std::string("internal.msh"), &DiodePDE::Instance::meshFileName);

// Doping file params:
  p.addPar("DOPING_FILE", std::string("NOFILE"), &DiodePDE::Instance::dopingFileName);
  p.addPar("PDOPE_FILE", std::string("NOFILE"), &DiodePDE::Instance::pdopeFileName);
  p.addPar("NDOPE_FILE", std::string("NOFILE"), &DiodePDE::Instance::ndopeFileName);
  p.addPar("NX", 11, &DiodePDE::Instance::NX);
  p.addPar("DIRICHLETBC", false, &DiodePDE::Instance::dirichletBCFlag)
    .setUnit(U_LOGIC)
    .setDescription("Flag for using Dirichlet boundary conditions for defects.");
  p.addPar("USEOLDNI", false, &DiodePDE::Instance::useOldNi)
    .setGivenMember(&DiodePDE::Instance::useOldNiGiven)
    .setDescription("Flag for using old(inaccurate) intrinsic carrier calculation.");

  p.addComposite("NODE", PDE_1DElectrode::getParametricData(), &DiodePDE::Instance::electrodeMap);
  p.addComposite("DOPINGPROFILES", DopeInfo::getParametricData(), &DiodePDE::Instance::dopeInfoMap);
  p.addComposite("REGION", DopeInfo::getParametricData(), &DiodePDE::Instance::dopeInfoMap);
}


//-----------------------------------------------------------------------------
// Function      : Instance::processParams
//
// Purpose       : This function contains much of the initialization for
//                 the Instance class.  Most of this was
//                 originally in the constructor.
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  updateTemperature(Temp);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       model,
  const FactoryBlock &factory_block)
  : DevicePDEInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(model),
    numMeshPoints(NUM_MESH_POINTS),
    NX(numMeshPoints),
    LX(NX-1),
    maxColsPerRow(MAX_COLS_PER_ROW),
    numElectrodes(2),
    NUMRC(NX*3),
    indicesSetup_(false),
    includeBaseNode_(false),
    useElectrodeSpec_(false),
    maskVarsTIAFlag_(false),
    scaleDensityToMaxDoping_(true),
    densityScalarFraction_(1.0e-1),
    useVoltageOutputOffset_(true),
    offsetWithFirstElectrode_(false),
    VoltageOffset_(0.0),
    Na(1.0e15),
    Nd(1.0e15),
    NnMax(1.0e15),
    NpMax(1.0e15),
    NnMin(1.0e5),  // approx...
    NpMin(1.0e5),
    WJ(1.0e-4),
    XC(0.0),
    XL(0.0),
    XR(0.0),
    Emax(0.0),
    VminExp(0.0),
    VmaxExp(0.0),
    diodeCap(0.0),
    junctionArea(1.0e-5),
    meshFileName(""),
    useOldNi(false),
    useOldNiGiven(false),
    dopingFileName("NOFILE"),
    ndopeFileName("NOFILE"),
    pdopeFileName("NOFILE"),
    width(1.0e-3),
    length(1.0e-3),
    basex(0.5e-3),
    area(1.0),
    anodebc(0.0),
    cathodebc(0.0),

    emitterbc(0.0),
    collectorbc(0.0),
    basebc(0.0),

    anodeArea(0.0),
    cathodeArea(0.0),

    emitterArea(0.0),
    collectorArea(0.0),
    baseArea(0.0),

    baseLocation(0.5e-3),
    baseLocationGiven(false),

    gradedJunctionFlag(false),
    bjtEnableFlag(false),
    displCurrentFlag(false),
    calledBeforeUIVB(false),
    callsOTEC(0),
    callsOSG(0),
    equationSet(0),
    lastOutputTime(-10.0),
    outputInterval(0.0),
    outputIntervalGiven(false),
    outputIndex(0),
    outputNLPoisson(false),
    outputRegion(0),
    tecplotLevel(0),
    sgplotLevel(0),
    voltLimFlag(false),
    includeAugerRecomb(true),
    includeSRHRecomb(true),

#ifdef Xyce_DEBUG_DEVICE
    anodeIndex_user(1),
    anodeIndex_userGiven(false),
    cathodeIndex_user(0),
    cathodeIndex_userGiven(false),
#endif

    maxVoltDelta(0.025), // thermal voltage.
    enableContinuationCalled(false),
    //coupledFlag(true),
    dirichletBCFlag(false),
    columnReorderingFlag(false)
{
  bcVec.clear();

  // these 4 mesh things change later.
  numIntVars   = 3*NUM_MESH_POINTS;
  numExtVars   = 2;
  if (IB.numExtVars != 0)
  {
    numExtVars   = IB.numExtVars;
  }
  numStateVars = 2;

  if (numExtVars < 3)
  {
    includeBaseNode_ = false;
  }
  else if (numExtVars == 3)
  {
    includeBaseNode_ = true;
  }
  else if (numExtVars > 3)
  {
    std::string msg = "DiodePDEInstance:";
    msg += "  name = " + getName();
    msg += " too many external nodes are set!  Set no more than 3.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // check doping files...
  if ( given("PDOPE_FILE") && !given("NDOPE_FILE") )
  {
    std::string msg = "Ndope file specified with no Pdope file.  Exiting.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  if ( !given("PDOPE_FILE") && given("NDOPE_FILE") )
  {
    std::string msg = "Pdope file specified with no Ndope file.  Exiting.";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    Temp = getDeviceOptions().temp.getImmutableValue<double>();

  if (given("MESHFILE"))
  {
    std::string msg = "Instance constructor."
                 "mesh file was specified.  The 1D device doesn't need a mesh file."
                 "  Either add a model statement of level=2, or get rid of the mesh"
                 " file specification.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  if (given("L") && !given("W"))
    width = length;
  if (given("GNUPLOTLEVEL") && !given("TECPLOTLEVEL"))
    tecplotLevel = gnuplotLevel;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  processParams ();

  // calculate dependent (ie computed) params and check for errors:
  ExtendedString tmpName = mobModelName;
  tmpName.toLower();
  mobModelName = tmpName;

  bool bsuccess = true;
  bool bs1 = true;

  bs1 = setupNumVars ();       bsuccess = bsuccess && bs1;
  bs1 = doAllocations ();      bsuccess = bsuccess && bs1;
  bs1 = setupMesh ();          bsuccess = bsuccess && bs1;
  bs1 = setupNodes ();         bsuccess = bsuccess && bs1;
  bs1 = setupDopingProfile (); bsuccess = bsuccess && bs1;
  bs1 = setupMiscConstants (); bsuccess = bsuccess && bs1;
  bs1 = setupScalingVars ();   bsuccess = bsuccess && bs1;

  bs1 = setupJacStamp ();      bsuccess = bsuccess && bs1;
  bs1 = cleanupJacStamp ();    bsuccess = bsuccess && bs1;

  if (!given("AREA")) area = 1.0;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
  for (std::map<std::string, DopeInfo *>::iterator it = dopeInfoMap.begin(); it != dopeInfoMap.end(); ++it)
    delete (*it).second;

  for (std::map<std::string, PDE_1DElectrode *>::iterator it = electrodeMap.begin(); it != electrodeMap.end(); ++it)
    delete (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : Instance::constructComposite
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/14/05
//-----------------------------------------------------------------------------
CompositeParam *
Instance::constructComposite(const std::string & compositeName, const std::string & paramName)
{
  if (compositeName == "DOPINGPROFILES" || compositeName == "REGION")
  {
    DopeInfo *n = new DopeInfo();
    dopeInfoMap[paramName] = n;
    return static_cast<CompositeParam *> (n);
  }
  else if (compositeName == "NODE" || compositeName == "ELECTRODE")
  {
    bcData bc;
    ExtendedString electrodeName = paramName;
    electrodeName.toUpper ();

    bc.eName = electrodeName;
    bc.nName = paramName;
    bc.given = true;
    bc.index = 0;

    if (electrodeName =="ANODE")
    {
      bc.meshIndex = 0;
      bc.neighborNode = 1;
    }
    else
    {
      bc.meshIndex = NUM_MESH_POINTS-1;
      bc.neighborNode = NUM_MESH_POINTS-2;
    }

    if (bc.given) ++numElectrodes;
    if (bc.given) bcVec.push_back(bc);

    PDE_1DElectrode *n = new PDE_1DElectrode();
    electrodeMap[paramName] = n;
    return static_cast<CompositeParam *> (n);
  }
  else
  {
    std::string msg =
      "Instance::constructComposite: unrecognized composite name: ";
    msg += compositeName;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  return NULL;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::doAllocations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/29/10
//-----------------------------------------------------------------------------
bool Instance::doAllocations ()
{
  // Set up a bunch of mesh-based arrays:
  dxVec.resize        (NX,0.0);
  xVec.resize         (NX,0.0);
  CVec.resize         (NX,0.0);
  CdonorVec.resize    (NX,0.0);
  CacceptorVec.resize (NX,0.0);
  VVec.resize         (NX,0.0);
  ExVec.resize        (NX,0.0);
  JnxVec.resize       (NX,0.0);
  JpxVec.resize       (NX,0.0);
  RVec.resize         (NX,0.0);
  SVec.resize         (NX,0.0);
  nnVec.resize        (NX,0.0);
  npVec.resize        (NX,0.0);

  dRdpVec.resize      (NX,0.0);
  dRdnVec.resize      (NX,0.0);

  dJndn1Vec.resize   (NX,0.0);
  dJndn2Vec.resize   (NX,0.0);
  dJndV1Vec.resize   (NX,0.0);
  dJndV2Vec.resize   (NX,0.0);
  dJndp1Vec.resize   (NX,0.0);
  dJndp2Vec.resize   (NX,0.0);

  dJpdn1Vec.resize   (NX,0.0);
  dJpdn2Vec.resize   (NX,0.0);
  dJpdV1Vec.resize   (NX,0.0);
  dJpdV2Vec.resize   (NX,0.0);
  dJpdp1Vec.resize   (NX,0.0);
  dJpdp2Vec.resize   (NX,0.0);

  //unVec.resize(NX,0.0);
  //upVec.resize(NX,0.0);
  tnVec.resize(NX,0.0);
  tpVec.resize(NX,0.0);
  unE_Vec.resize (NX-1,0.0);
  upE_Vec.resize (NX-1,0.0);

  // indexing arrays, local.  jacStamp is resized elsewhere
  li_Vrowarray.resize(NX,0);
  li_Nrowarray.resize(NX,0);
  li_Prowarray.resize(NX,0);

  // displacement current stuff
  stateDispl.resize(NX,0);
  stateDispl_owned.resize(NX,0);
  displCurrent.resize(NX,0.0);
  li_stateDispl.resize(NX,0);

  // set up the boundary stencil:
  boundarySten.resize(NX,0);
  edgeBoundarySten.resize(NX,0);
  internalBoundarySten.resize(NX,0);

  // these will always be set.
  edgeBoundarySten[0]=1;
  edgeBoundarySten[LX]=1;
  boundarySten[0]=1;
  boundarySten[LX]=1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupNodes
//
// Purpose       : This sets up the bcVec container.  bcVec is a vector of
//                 bcData classes, which contain bondary condition
//                 related data.
//
//                 The key issues for a boundary are:
//                    - determine circuit node
//                    - determine mesh boundary location
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/29/10
//-----------------------------------------------------------------------------
bool Instance::setupNodes ()
{
  // If the user did not use the ELECTRODE/NODE vector-composite
  // specification, then set up a default (implicit) set of electrodes.
  //
  // If 2 terminals, assume a diode, with a specification consistent with the
  // SPICE diode noder order:  D cathode anode
  //
  // If 3 terminals, assume a BJT with the node order being the same as
  // the SPICE Gummel-Poon specification:  Q col bas emit
  //
  if ( bcVec.empty() )
  {
    bcVec.clear();
    bcVec.resize(numExtVars);

    useElectrodeSpec_ = false;

    if (includeBaseNode_)
    {
      // collector:
      int collectorIndex=0;
      bcIndexMap["collector"] = collectorIndex;
      bcVec[collectorIndex].eName = "collector";
      bcVec[collectorIndex].Vequ = collectorbc;
      bcVec[collectorIndex].VequGiven = given("COLLECTOR.BC");
      bcVec[collectorIndex].area = collectorArea;
      bcVec[collectorIndex].areaGiven = given("COLLECTOR.AREA");
      bcVec[collectorIndex].meshIndex = LX;
      bcVec[collectorIndex].neighborNode = LX-1;
      if (!given("COLLECTOR.AREA")) bcVec[collectorIndex].area = area;

      // base:
      int baseIndex=1;
      bcIndexMap["base"] = baseIndex;
      bcVec[baseIndex].eName = "base";
      bool found=false;
      int bIndex=0;
      double minDelta = length;
      for (int i=0;i<NX;++i)
      {
        double deltaX=fabs(baseLocation-xVec[i]);
        if (deltaX < minDelta)
        {
          bIndex=i;
          minDelta=deltaX;
        }
      }

      bcVec[baseIndex].Vequ = basebc;
      bcVec[baseIndex].VequGiven = given("BASE.BC");
      bcVec[baseIndex].area = baseArea;
      bcVec[baseIndex].areaGiven = given("BASE.AREA");
      //bcVec[baseIndex].meshIndex = static_cast<int> (LX/2);
      bcVec[baseIndex].meshIndex = bIndex;
      bcVec[baseIndex].neighborNode = bcVec[baseIndex].meshIndex-1;
      if (!given("BASE.AREA")) bcVec[baseIndex].area = area;

      // emitter:
      int emitterIndex=2;
      bcIndexMap["emitter"] = emitterIndex;
      bcVec[emitterIndex].eName = "emitter";
      bcVec[emitterIndex].Vequ = emitterbc;
      bcVec[emitterIndex].VequGiven = given("EMITTER.BC");
      bcVec[emitterIndex].area = emitterArea;
      bcVec[emitterIndex].areaGiven = given("EMITTER.AREA");
      bcVec[emitterIndex].meshIndex = 0;
      bcVec[emitterIndex].neighborNode = 1;
      if (!given("EMITTER.AREA")) bcVec[emitterIndex].area = area;
    }
    else
    {
      // anode:
      int anodeIndex=1;
      bcIndexMap["anode"] = anodeIndex;
      bcVec[anodeIndex].eName = "anode";
      bcVec[anodeIndex].Vequ = anodebc;
      bcVec[anodeIndex].VequGiven = given("ANODE.BC");
      bcVec[anodeIndex].area = anodeArea;
      bcVec[anodeIndex].areaGiven = given("ANODE.AREA");
      bcVec[anodeIndex].meshIndex = 0;
      bcVec[anodeIndex].neighborNode = 1;
      if (!given("ANODE.AREA")) bcVec[anodeIndex].area = area;

      // cathode:
      int cathodeIndex=0;
      bcIndexMap["cathode"] = cathodeIndex;
      bcVec[cathodeIndex].eName = "cathode";
      bcVec[cathodeIndex].Vequ = cathodebc;
      bcVec[cathodeIndex].VequGiven = given("CATHODE.BC");
      bcVec[cathodeIndex].area = cathodeArea;
      bcVec[cathodeIndex].areaGiven = given("CATHODE.AREA");
      bcVec[cathodeIndex].meshIndex = LX;
      bcVec[cathodeIndex].neighborNode = LX-1;
      if (!given("CATHODE.AREA")) bcVec[cathodeIndex].area = area;
    }
  }
  else  // user did use the ELECTRODE/NODE specification.
  {
    useElectrodeSpec_ = true;

    std::vector<int> tmpMeshSten(NX,0);

    for (int iBC=0;iBC<bcVec.size();++iBC)
    {
      PDE_1DElectrode & electrode = *(electrodeMap[bcVec[iBC].nName]);

      bcIndexMap[ bcVec[iBC].eName ] = iBC;

      if (electrode.sideGiven)
      {
        ExtendedString side = electrode.side;
        side.toLower();
        if (side == "left")
        {
          bcVec[iBC].meshIndex = 0;
          bcVec[iBC].neighborNode = 1;
          tmpMeshSten[0] = 1;
        }
        else if (side == "right")
        {
          bcVec[iBC].meshIndex = LX;
          bcVec[iBC].neighborNode = LX-1;
          tmpMeshSten[LX] = 1;
        }
        else if (side == "middle" || side == "mid")
        {

          double location = electrode.location;
          bool found=false;
          int bIndex=0;
          double minDelta = length;
          for (int imesh=0;imesh<NX;++imesh)
          {
            double deltaX=fabs(location-xVec[imesh]);
            if (deltaX < minDelta)
            {
              bIndex=imesh;
              minDelta=deltaX;
            }
          }

          bcVec[iBC].meshIndex = bIndex;
          bcVec[iBC].neighborNode = bIndex-1;
          // assuming current coming from the emitter direciton

          // check to make sure that bIndex isn't already used.
          if (tmpMeshSten[bIndex] == 1)
          {
            std::string msg = "Instance::setupNodes.  Failed to find mesh index for " + bcVec[iBC].eName;
            N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
          }
        }
        else
        {
          std::string msg = "Instance::setupNodes.  unrecognized side specified.";
          N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
        }
      }
      else
      {
        //std::string msg = "Instance::setupNodes.  side NOT specified.";
        //N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
      }

      bcVec[iBC].areaGiven = electrode.areaGiven;
      if (electrode.areaGiven)
      {
        bcVec[iBC].area = electrode.area;
      }
    }
  }

  indicesSetup_=true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << " area               = " << area << std::endl;
    Xyce::dout() << " areaGiven          = " << given("AREA") << std::endl;
    int isize=bcVec.size();
    for (int i=0;i<isize;++i)
    {
      Xyce::dout() << " bcVec["<<i<<"].area      = " << bcVec[i].area << std::endl;
      Xyce::dout() << " bcVec["<<i<<"].areaGiven = " << bcVec[i].areaGiven << std::endl;
      Xyce::dout() << " bcVec["<<i<<"].meshIndex = " << bcVec[i].meshIndex << std::endl;
    }
  }
#endif

  int colmax = maxColsPerRow;
  int bcSize=bcVec.size();
  for (int i=0;i<bcSize;++i)
  {
    bcVec[i].colArray.resize(colmax,-1);
    bcVec[i].dIdXcols.resize(colmax,-1);
    bcVec[i].dIdX.resize(colmax,-1);
    bcVec[i].dFdVckt.resize(colmax,0.0);
  }

  // allocate conductance array:
  numElectrodes = bcVec.size(); // should be n x n,
  // where n=number of terminals.
  condVec.resize(numElectrodes);
  for (int iE=0;iE<numElectrodes;++iE)
  {
    condVec[iE].resize(numElectrodes,0.0);
  }


  // initialize the boundary stencils.
  // Note: two of the points will be at meshIndex=0 and meshIndex=LX.  If there
  // is a 3rd terminal (for the base of a BJT) it will be somewhere in the middle.

  for (int i=0;i<bcSize;++i)
  {
    int meshIndex=bcVec[i].meshIndex;

    if (meshIndex==0 || meshIndex==LX)
    {
      edgeBoundarySten[meshIndex]=1;
    }
    else
    {
      internalBoundarySten[meshIndex]=1;
    }
    boundarySten[meshIndex]=1;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupNumVars
//
// Purpose       : mostly sets up numIntVars.   numExtVars was set earlier
//                 and is easy.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/11/03
//-----------------------------------------------------------------------------
bool Instance::setupNumVars ()
{
  if (given("NX"))
  {
    LX = NX-1;

    if (NX != numMeshPoints)
    {
      numMeshPoints = NX;
      numIntVars    = 3*NX;
    }

    //numExtVars   = IB.numExtVars;
    //numStateVars  = 2 + NX - 1; // the NX-1 is for the displacement current.
    numStateVars  = numExtVars + NX - 1; // the NX-1 is for the displacement current.
    maxColsPerRow = MAX_COLS_PER_ROW;
  }
  else
  {
    std::string msg = "Instance constructor."
                 "  NX parameter was not specified.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupJacStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/11/03
//-----------------------------------------------------------------------------
bool Instance::setupJacStamp ()
{
  int iMeshNode;
  int numVars   = 3;
  int Voffset   = 0;
  int Noffset   = 1;
  int Poffset   = 2;
  int baseIndex = 0;
  int baseIndex_m1 = 0;
  int baseIndex_p1 = 0;

  int Vindex    = 0;
  int Nindex    = 0;
  int Pindex    = 0;

  // "minus one" indices
  int Vindex_m1    = 0;
  int Nindex_m1    = 0;
  int Pindex_m1    = 0;

  // "plus one" indices
  int Vindex_p1    = 0;
  int Nindex_p1    = 0;
  int Pindex_p1    = 0;

  int extVarOffset = numExtVars;
  int iBC;

  int jacSize = numIntVars + extVarOffset;

  meshToLID.clear(); meshToLID.resize(NX,-1);
  jacStamp.clear();   jacStamp.resize(jacSize);

  // set up the meshToLID converter first.
  int bcSize=bcVec.size();
  int lid=0;
  for (int iBC=0;iBC<bcSize;++iBC)
  {
    int meshIndex=bcVec[iBC].meshIndex;
    meshToLID[meshIndex] = lid;
    lid++;
  }

  for (int i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;
    meshToLID[i] = lid;
    lid++;
  }

  // external vars (from circuit) first:
  // coupled mode==1 is handled in the function augJacStampRxnChem, not here.
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    iMeshNode = bcVec[iBC].meshIndex;

    if (edgeBoundarySten[iMeshNode]!=1 &&  internalBoundarySten[iMeshNode]!=1)
    {
      std::string msg = "Instance::setupJacStamp:";
      msg += "Boundary point not in the stencil.";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }

    int iNN   = bcVec[iBC].neighborNode;
    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
    Vindex = baseIndex + Voffset;
    Nindex = baseIndex + Noffset;
    Pindex = baseIndex + Poffset;

    if (iMeshNode > iNN) // i=0, or right-looking
    {
      baseIndex_m1 = numVars*meshToLID[iMeshNode-1] + extVarOffset;

      Vindex_m1 = baseIndex_m1 + Voffset;
      Nindex_m1 = baseIndex_m1 + Noffset;
      Pindex_m1 = baseIndex_m1 + Poffset;

      int col=0;
      jacStamp[iBC].resize(7,-1);
      jacStamp[iBC][col++] = iBC;
      jacStamp[iBC][col++] = Vindex;
      jacStamp[iBC][col++] = Vindex_m1;
      jacStamp[iBC][col++] = Nindex;
      jacStamp[iBC][col++] = Nindex_m1;
      jacStamp[iBC][col++] = Pindex;
      jacStamp[iBC][col++] = Pindex_m1;
    }
    else // i=LX, or left-looking
    {
      baseIndex_p1 = numVars*meshToLID[iMeshNode+1] + extVarOffset;

      Vindex_p1 = baseIndex_p1 + Voffset;
      Nindex_p1 = baseIndex_p1 + Noffset;
      Pindex_p1 = baseIndex_p1 + Poffset;

      int col=0;
      jacStamp[iBC].resize(7,-1);
      jacStamp[iBC][col++] = iBC;
      jacStamp[iBC][col++] = Vindex;
      jacStamp[iBC][col++] = Vindex_p1;
      jacStamp[iBC][col++] = Nindex;
      jacStamp[iBC][col++] = Nindex_p1;
      jacStamp[iBC][col++] = Pindex;
      jacStamp[iBC][col++] = Pindex_p1;
    }
  }

  // Variables associated with the mesh.
  // First do the mesh points that are boundary condition points.
  // Do the anode (BC) mesh point:

  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    iMeshNode = bcVec[iBC].meshIndex;
    int iNN   = bcVec[iBC].neighborNode;
    int NodeIndex = bcIndexMap[bcVec[iBC].eName]; // iBC?

    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
    Vindex = baseIndex + Voffset;
    Nindex = baseIndex + Noffset;
    Pindex = baseIndex + Poffset;

    jacStamp[Vindex].resize(3,-1);
    jacStamp[Nindex].resize(4,-1);
    jacStamp[Pindex].resize(4,-1);

    if (edgeBoundarySten[iMeshNode]==1)
    {
      if (iMeshNode < iNN) // i=0
      {
        baseIndex_p1 = numVars*meshToLID[iMeshNode+1] + extVarOffset;
        Vindex_p1 = baseIndex_p1 + Voffset;
        Nindex_p1 = baseIndex_p1 + Noffset;
        Pindex_p1 = baseIndex_p1 + Poffset;

        int col=0;
        jacStamp[Vindex][col++] = NodeIndex;
        jacStamp[Vindex][col++] = Vindex;
        jacStamp[Vindex][col++] = Vindex_p1;

        col=0;
        jacStamp[Nindex][col++] = Nindex;
        jacStamp[Nindex][col++] = Nindex_p1;
        jacStamp[Nindex][col++] = Pindex;
        jacStamp[Nindex][col++] = Pindex_p1;

        col=0;
        jacStamp[Pindex][col++] = Pindex;
        jacStamp[Pindex][col++] = Pindex_p1;
        jacStamp[Pindex][col++] = Nindex;
        jacStamp[Pindex][col++] = Nindex_p1;
      }
      else // i=LX
      {
        baseIndex_m1 = numVars*meshToLID[iMeshNode-1] + extVarOffset;
        Vindex_m1 = baseIndex_m1 + Voffset;
        Nindex_m1 = baseIndex_m1 + Noffset;
        Pindex_m1 = baseIndex_m1 + Poffset;

        int col=0;
        jacStamp[Vindex][col++] = Vindex_m1;
        jacStamp[Vindex][col++] = Vindex;
        jacStamp[Vindex][col++] = NodeIndex;

        col=0;
        jacStamp[Nindex][col++] = Nindex_m1;
        jacStamp[Nindex][col++] = Nindex;
        jacStamp[Nindex][col++] = Pindex_m1;
        jacStamp[Nindex][col++] = Pindex;

        col=0;
        jacStamp[Pindex][col++] = Pindex_m1;
        jacStamp[Pindex][col++] = Pindex;
        jacStamp[Pindex][col++] = Nindex_m1;
        jacStamp[Pindex][col++] = Nindex;
      }
    }
    else if (internalBoundarySten[iMeshNode]==1) // probably base node
    {
      // base node applies a BC to the potential and majority carrier.
      // The minority carrier does not get a BC, and is treated like an
      // internal point.  So the stamp here is the same as an interior
      // point plus a little extra.  Not all of these will be used.

      baseIndex_m1 = numVars*meshToLID[iMeshNode-1] + extVarOffset;
      baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
      baseIndex_p1 = numVars*meshToLID[iMeshNode+1] + extVarOffset;

      Vindex_m1 = baseIndex_m1 + Voffset;
      Nindex_m1 = baseIndex_m1 + Noffset;
      Pindex_m1 = baseIndex_m1 + Poffset;
      Vindex    = baseIndex    + Voffset;
      Nindex    = baseIndex    + Noffset;
      Pindex    = baseIndex    + Poffset;
      Vindex_p1 = baseIndex_p1 + Voffset;
      Nindex_p1 = baseIndex_p1 + Noffset;
      Pindex_p1 = baseIndex_p1 + Poffset;

      // voltage col arrays:
      int col=0;
      jacStamp[Vindex].resize(6,-1);
      jacStamp[Vindex][col++] = NodeIndex;
      jacStamp[Vindex][col++] = Vindex_m1;
      jacStamp[Vindex][col++] = Vindex;
      jacStamp[Vindex][col++] = Vindex_p1;
      jacStamp[Vindex][col++] = Nindex;
      jacStamp[Vindex][col++] = Pindex;

      // electron col arrays:
      col=0;
      jacStamp[Nindex].resize(9,-1);
      jacStamp[Nindex][col++] = Nindex_m1;
      jacStamp[Nindex][col++] = Nindex;
      jacStamp[Nindex][col++] = Nindex_p1;
      jacStamp[Nindex][col++] = Vindex_m1;
      jacStamp[Nindex][col++] = Vindex;
      jacStamp[Nindex][col++] = Vindex_p1;
      jacStamp[Nindex][col++] = Pindex_m1;
      jacStamp[Nindex][col++] = Pindex;
      jacStamp[Nindex][col++] = Pindex_p1;

      // hole col arrays:
      col=0;
      jacStamp[Pindex].resize(9,-1);
      jacStamp[Pindex][col++] = Pindex_m1;
      jacStamp[Pindex][col++] = Pindex;
      jacStamp[Pindex][col++] = Pindex_p1;
      jacStamp[Pindex][col++] = Vindex_m1;
      jacStamp[Pindex][col++] = Vindex;
      jacStamp[Pindex][col++] = Vindex_p1;
      jacStamp[Pindex][col++] = Nindex_m1;
      jacStamp[Pindex][col++] = Nindex;
      jacStamp[Pindex][col++] = Nindex_p1;

    }
    else// not a boundary.  oops!
    {
      std::string msg = "Instance::setupJacStamp:";
      msg += "Boundary point not in the stencil.";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
  }

  // Now do the non-BC mesh points.
  for (iMeshNode=0;iMeshNode<NX;++iMeshNode)
  {
    if (boundarySten[iMeshNode]==1) continue;

    baseIndex_m1 = numVars*meshToLID[iMeshNode-1] + extVarOffset;
    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
    baseIndex_p1 = numVars*meshToLID[iMeshNode+1] + extVarOffset;

    Vindex_m1 = baseIndex_m1 + Voffset;
    Nindex_m1 = baseIndex_m1 + Noffset;
    Pindex_m1 = baseIndex_m1 + Poffset;
    Vindex    = baseIndex    + Voffset;
    Nindex    = baseIndex    + Noffset;
    Pindex    = baseIndex    + Poffset;
    Vindex_p1 = baseIndex_p1 + Voffset;
    Nindex_p1 = baseIndex_p1 + Noffset;
    Pindex_p1 = baseIndex_p1 + Poffset;

    // voltage col arrays:
    int col=0;
    jacStamp[Vindex].resize(5,-1);
    jacStamp[Vindex][col++] = Vindex_m1;
    jacStamp[Vindex][col++] = Vindex;
    jacStamp[Vindex][col++] = Vindex_p1;
    jacStamp[Vindex][col++] = Nindex;
    jacStamp[Vindex][col++] = Pindex;

    // electron col arrays:
    col=0;
    jacStamp[Nindex].resize(9,-1);
    jacStamp[Nindex][col++] = Nindex_m1;
    jacStamp[Nindex][col++] = Nindex;
    jacStamp[Nindex][col++] = Nindex_p1;
    jacStamp[Nindex][col++] = Vindex_m1;
    jacStamp[Nindex][col++] = Vindex;
    jacStamp[Nindex][col++] = Vindex_p1;
    jacStamp[Nindex][col++] = Pindex_m1;
    jacStamp[Nindex][col++] = Pindex;
    jacStamp[Nindex][col++] = Pindex_p1;

    // hole col arrays:
    col=0;
    jacStamp[Pindex].resize(9,-1);
    jacStamp[Pindex][col++] = Pindex_m1;
    jacStamp[Pindex][col++] = Pindex;
    jacStamp[Pindex][col++] = Pindex_p1;
    jacStamp[Pindex][col++] = Vindex_m1;
    jacStamp[Pindex][col++] = Vindex;
    jacStamp[Pindex][col++] = Vindex_p1;
    jacStamp[Pindex][col++] = Nindex_m1;
    jacStamp[Pindex][col++] = Nindex;
    jacStamp[Pindex][col++] = Nindex_p1;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1)
  {
    // dump the jacStamp to stdout:
    int jacSize = jacStamp.size ();
    Xyce::dout() << "jacStamp size = " << jacSize << std::endl;

    for(int i=0;i<jacSize;++i)
    {
      int colSize = jacStamp[i].size();
      for (int j=0;j<colSize;++j)
      {
        Xyce::dout() << "  jacStamp["<<i<<"]["<<j<<"] = " << jacStamp[i][j] << std::endl;
      }
    }
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::cleanupJacStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/12/10
//-----------------------------------------------------------------------------
bool Instance::cleanupJacStamp ()
{
#if 1
  // set up normal jacMap for when all resistances nonzero
  // If nothing is remapped, this amounts to a null operation when the
  // map is used later.  The maps become important when we start
  // remapping nodes because of zero lead resistances
  jacMap.clear();
  jacMap2.clear();
  jacMap.resize(jacStamp.size());
  jacMap2.resize(jacStamp.size());

  int mapSize = jacMap.size();
  for (int i=0;i<mapSize;++i)
  {
    jacMap[i]=i;
    jacMap2[i].resize(jacStamp[i].size());
    for (int j=0;j<jacStamp[i].size();++j)
    {
      jacMap2[i][j] = j;
    }
  }

  // Now fix the ordering of the columns in the jacStamp.  If the columns in each row
  // are not in ascending order, then the jacStampMap calls below (for removing
  // variables) will not work correctly.
  //
  // NOTE:  This is probably not safe to do, for the PDE devices, so column reordering is off
  // by default.
  if (columnReorderingFlag)
  {
    std::vector< std::vector<int> > tempStamp_eric;
    std::vector< std::vector<int> > tempMap2_eric;
    jacStampMap_fixOrder(jacStamp, jacMap2, tempStamp_eric, tempMap2_eric);
    jacStamp = tempStamp_eric;
    jacMap2 = tempMap2_eric;
  }

#endif // if 1

  return true;
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
  if (intNameMap.empty ())
  {
    int itmp;
    char tmpchar[128];

    for (itmp=0;itmp<128;++itmp) tmpchar[itmp] = 0;

    // set up the name map for this device.
    for (itmp=0;itmp<NX;++itmp)
    {
      if (li_Vrowarray[itmp] != -1)
      {
        sprintf(tmpchar,"%s_V_%d",getName().c_str(), itmp);
        intNameMap[li_Vrowarray[itmp]] = std::string(tmpchar);
      }

      if (li_Nrowarray[itmp] != -1)
      {
        sprintf(tmpchar,"%s_N_%d",getName().c_str(), itmp);
        intNameMap[li_Nrowarray[itmp]] = std::string(tmpchar);
      }

      if (li_Prowarray[itmp] != -1)
      {
        sprintf(tmpchar,"%s_P_%d",getName().c_str(), itmp);
        intNameMap[li_Prowarray[itmp]] = std::string(tmpchar);
      }
    }
  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : DiodeDPEInstance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::registerLIDs:\n";
    Xyce::dout() << "          name = " << getName() << std::endl;
    Xyce::dout() << "        numInt = " << numIntVars << std::endl;
    Xyce::dout() << "        numEXt = " << numExtVars << std::endl;
    Xyce::dout() << "        NX     = " << NX << std::endl;
  }
#endif

  // Copy over the local ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    bcVec[iBC].lid = extLIDVec[iBC];
  }

  // First do the boundary condition mesh points.  There will be imposed
  // boundary conditions at each boundary
  // for potential, hole density and electron density.

  // The electrostatic potential, from the perspective of the device
  // simulation, is not neccessarily the same as the voltage used in the
  // circuit part of the code.  Also, obviously, the densities are not used by
  // the circuit sim., so these boundary conditions are considered internal
  // variables.

  int meshIndex = 0;
  int intIndex = 0;

  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    meshIndex = bcVec[iBC].meshIndex;
    li_Vrowarray[meshIndex] = intLIDVec[intIndex++];
    li_Nrowarray[meshIndex] = intLIDVec[intIndex++];
    li_Prowarray[meshIndex] = intLIDVec[intIndex++];
  }

  // now do the interior points.  These will be blocked (V,N,P) together.
  //meshIndex=1;
  //while (meshIndex < LX)
  for (meshIndex=0;meshIndex<NX;++meshIndex)
  {
    if (boundarySten[meshIndex]==1) continue;

    li_Vrowarray[meshIndex] = intLIDVec[intIndex++];
    li_Nrowarray[meshIndex] = intLIDVec[intIndex++];
    li_Prowarray[meshIndex] = intLIDVec[intIndex++];
  }

#ifdef Xyce_DEBUG_DEVICE

  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << "\n  solution indices:\n";

    for (int i=0;i<NX;++i)
    {
      Xyce::dout() << "     li_Vrowarray["<<i<<"] = " << li_Vrowarray[i];
      Xyce::dout() << "\tli_Nrowarray["<<i<<"] = " << li_Nrowarray[i];
      Xyce::dout() << "\tli_Prowarray["<<i<<"] = " << li_Prowarray[i] << std::endl;
    }
    Xyce::dout() << section_divider << std::endl;
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "  In Instance::registerStateLIDs\n\n";
    Xyce::dout() << "  name             = " << getName() << std::endl;
    Xyce::dout() << "  Number of State LIDs: " << numStateVars << std::endl;
  }
#endif

  // Copy over the local ID lists:
  staLIDVec = staLIDVecRef;

  int i;
  int j;
  for (i=0;i<bcVec.size();++i)
  {
    bcVec[i].li_stateC = staLIDVec[i];
  }

  for (i=0,j=2;i<NX-1;++i,++j)
  {
    li_stateDispl[i] = staLIDVec[j];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "  State indices:" << std::endl;
    Xyce::dout() << std::endl;
    for (i=0;i<bcVec.size();++i)
    {
      Xyce::dout() << "bcVec["<<i<<"].li_stateC = "<<bcVec[i].li_stateC<< std::endl;
    }
    Xyce::dout() << std::endl;

    Xyce::dout() << "  Displacement current state variable local indices:" << std::endl;
    for (i=0;i<NX-1;++i)
    {
      Xyce::dout() << "  li_stateDispl["<<i<<"] = " << li_stateDispl[i] << std::endl;
    }
    Xyce::dout() << section_divider << std::endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233
// Creation Date : 02/11/03
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233
// Creation Date : 02/12/03
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs
( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs ( jacLIDVec );

  int iMeshNode;
  int numVars   = 3;
  int baseIndex;
  int Vindex;
  int Nindex;
  int Pindex;
  int Voffset = 0;
  int Noffset = 1;
  int Poffset = 2;
#ifdef Xyce_DEBUG_DEVICE
  int i,j;
#endif

  int extVarOffset = numExtVars;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::registerJacLIDs" << std::endl;

    int jacLIDSize = jacLIDVec.size();
    Xyce::dout() << "jacLIDSize = " << jacLIDSize << std::endl;
    for (i=0;i<jacLIDSize;++i)
    {
      int jacLIDcolSize = jacLIDVec[i].size();
      Xyce::dout() << std::endl;
      Xyce::dout() << "jacLIDVec["<<i<<"].size = " << jacLIDcolSize << std::endl;
      for (j=0;j<jacLIDcolSize;++j)
      {
        Xyce::dout() << "jacLIDVec["<<i<<"]["<<j<<"] = ";
        Xyce::dout() << jacLIDVec[i][j] << std::endl;
      }
    }
  }
#endif

  li_Vcolarray.resize(NX);
  li_Ncolarray.resize(NX);
  li_Pcolarray.resize(NX);


  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i1=0;
    int numCols = jacLIDVec[iBC].size();
    bcVec[iBC].li_colArray.resize(numCols,-1);

    for (i1=0;i1<numCols;++i1)
    {
      bcVec[iBC].li_colArray[i1] = jacLIDVec[iBC][i1];
    }

    int iMeshNode = bcVec[iBC].meshIndex;

    bcVec[iBC].lidOffset = jacLIDVec[iBC][0];

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << std::endl;
      for(i=0;i<bcVec[iBC].li_colArray.size();++i)
      {
        Xyce::dout() << bcVec[iBC].eName << ": li_colArray["<<i<<"] = "<<bcVec[iBC].li_colArray[i]<< std::endl;
      }
      Xyce::dout() << std::endl;
    }
#endif
  }

  // Do the non-BC mesh points.
  for (iMeshNode=0;iMeshNode<NX;++iMeshNode)
  {
    if (boundarySten[iMeshNode]==1) continue;

    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;

    Vindex = baseIndex + Voffset;
    Nindex = baseIndex + Noffset;
    Pindex = baseIndex + Poffset;

    // voltage col arrays:
    int i1=0;
    int j1=0;
    li_Vcolarray[iMeshNode].resize(5,-1);
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];

    // electron col arrays:
    i1=0; j1=0;
    li_Ncolarray[iMeshNode].resize(9,-1);
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];


    // hole col arrays:
    i1=0; j1=0;
    li_Pcolarray[iMeshNode].resize(9,-1);
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << std::endl;
      Xyce::dout() << "registerJacLIDs: iMeshNode = " << iMeshNode << std::endl;
      Xyce::dout() << "jacLIDVec[Vindex].size = " << jacLIDVec[Vindex].size()<<std::endl;
      Xyce::dout() << "jacLIDVec[Nindex].size = " << jacLIDVec[Nindex].size()<<std::endl;
      Xyce::dout() << "jacLIDVec[Pindex].size = " << jacLIDVec[Pindex].size()<<std::endl;

      for (i=0;i<5;++i)
      {
        Xyce::dout() << " li_Vcolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Vcolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
      for (i=0;i<7;++i)
      {
        Xyce::dout() << " li_Ncolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Ncolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
      for (i=0;i<7;++i)
      {
        Xyce::dout() << " li_Pcolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Pcolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
    }
#endif
  }

  // Do the BC mesh points.  These are actually "first" in the LID array.
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    iMeshNode = bcVec[iBC].meshIndex;
    int iNN   = bcVec[iBC].neighborNode;

    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
    Vindex = baseIndex + Voffset;
    Nindex = baseIndex + Noffset;
    Pindex = baseIndex + Poffset;

    if (edgeBoundarySten[iMeshNode]==1)
    {
      int i1=0;
      int j1=0;
      li_Vcolarray[iMeshNode].resize(3,-1);
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];

      // The stamp for N and P is dependent on the direction.
      if (iMeshNode < iNN)// i=0
      {
        i1=0; j1=0;
        li_Ncolarray[iMeshNode].resize(3,-1);
        li_Ncolarray[iMeshNode][i1++] = -1;
        li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
        li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];

        i1=0; j1=0;
        li_Pcolarray[iMeshNode].resize(3,-1);
        li_Pcolarray[iMeshNode][i1++] = -1;
        li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
        li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      }
      else
      {
        i1=0; j1=0;
        li_Ncolarray[iMeshNode].resize(3,-1);
        li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
        li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];

        i1=0; j1=0;
        li_Pcolarray[iMeshNode].resize(3,-1);
        li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
        li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      }
    }
    else if (internalBoundarySten[iMeshNode]==1) // probably base node
    {
      // voltage col arrays:
      int i1=0;
      int j1=0;
      li_Vcolarray[iMeshNode].resize(6,-1);
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];

      // electron col arrays:
      i1=0; j1=0;
      li_Ncolarray[iMeshNode].resize(9,-1);
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];

      // hole col arrays:
      i1=0; j1=0;
      li_Pcolarray[iMeshNode].resize(9,-1);
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    }
    else// not a boundary.  oops!
    {
      std::string msg = "Instance::registerJacLIDs:";
      msg += "Boundary point not in the stencil.";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << std::endl;
      Xyce::dout() << "registerJacLIDs: ("<<bcVec[iBC].eName<<") iMeshNode = ";
      Xyce::dout() << iMeshNode << std::endl;
      for (i=0;i<3;++i)
      {
        Xyce::dout() << " li_Vcolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Vcolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
      for (i=0;i<3;++i)
      {
        Xyce::dout() << " li_Ncolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Ncolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
      for (i=0;i<3;++i)
      {
        Xyce::dout() << " li_Pcolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Pcolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
    }
#endif
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers()
//
// Purpose       : Sets up raw pointers for optimized matrix loads.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/10
//-----------------------------------------------------------------------------
void Instance::setupPointers()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  fVmatPtr.resize(NX);
  fNmatPtr.resize(NX);
  fPmatPtr.resize(NX);
  qVmatPtr.resize(NX);
  qNmatPtr.resize(NX);
  qPmatPtr.resize(NX);

  for (int i=0;i<NX;++i)
  {
    int Vrow = li_Vrowarray[i];
    int Nrow = li_Nrowarray[i];
    int Prow = li_Prowarray[i];

    int vSize = li_Vcolarray[i].size();
    fVmatPtr[i].resize(vSize);
    qVmatPtr[i].resize(vSize);
    for (int j=0;j<vSize;++j)
    {
      fVmatPtr[i][j] = &(dFdx[Vrow][li_Vcolarray[i][j]]);
      qVmatPtr[i][j] = &(dQdx[Vrow][li_Vcolarray[i][j]]);
    }

    int nSize = li_Ncolarray[i].size();
    fNmatPtr[i].resize(nSize);
    qNmatPtr[i].resize(nSize);
    for (int j=0;j<nSize;++j)
    {
      fNmatPtr[i][j] = &(dFdx[Nrow][li_Ncolarray[i][j]]);
      qNmatPtr[i][j] = &(dQdx[Nrow][li_Ncolarray[i][j]]);
    }

    int pSize = li_Pcolarray[i].size();
    fPmatPtr[i].resize(pSize);
    qPmatPtr[i].resize(pSize);
    for (int j=0;j<pSize;++j)
    {
      fPmatPtr[i][j] = &(dFdx[Prow][li_Pcolarray[i][j]]);
      qPmatPtr[i][j] = &(dQdx[Prow][li_Pcolarray[i][j]]);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermdiateVars
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;
  bool bs1 = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "updateIntermediateVars.  name = " << getName() << std::endl;
  }
#endif

  bs1 = obtainSolution ();        bsuccess = bsuccess && bs1;
  bs1 = calcEfield ();            bsuccess = bsuccess && bs1;
  bs1 = calcMobilities   ();      bsuccess = bsuccess && bs1;
  bs1 = calcRecombination ();     bsuccess = bsuccess && bs1;
  bs1 = calcElectronCurrent ();   bsuccess = bsuccess && bs1;
  bs1 = calcHoleCurrent ();       bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcTerminalCurrents
// Purpose       : Calculates total diode current(s) to be used in the
//                 circuit KCL equations.
//
// Special Notes : Two options:
//
//  1) use the fluxes from the PDE calculation
//  2) use the integrated emission/capture rates from the rxn network.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/04/01
//-----------------------------------------------------------------------------
bool Instance::calcTerminalCurrents ()
{
  bool bsuccess = true;
  double & J0 = scalingVars.J0;
  double & a0 = scalingVars.a0;

  // Calculate the diode current using DD fluxes.
  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    int index = bcVec[iBC].meshIndex;
    int iNN=bcVec[iBC].neighborNode;
    double & area = bcVec[iBC].area;
    double A0=J0*a0*area;

    double sign =  ((iNN > index)?1.0:-1.0);
    int edgeIndex= ((iNN > index)?index:iNN);

    bcVec[iBC].elecCurrent = sign*JnxVec[edgeIndex]*A0;
    bcVec[iBC].holeCurrent = sign*JpxVec[edgeIndex]*A0;

    if (edgeBoundarySten[index]==1)
    {
      bcVec[iBC].currentSum = bcVec[iBC].elecCurrent + bcVec[iBC].holeCurrent;
    }
    else if (internalBoundarySten[index]==1)
    {
      std::string & type = bcVec[iBC].type;

      // only majority carrier goes to the boundary
      if (type=="ntype")
      {
        bcVec[iBC].currentSum = bcVec[iBC].elecCurrent;
      }
      else if (type=="ptype")
      {
        bcVec[iBC].currentSum = bcVec[iBC].holeCurrent;
      }
      else // oops.
      {
        std::string msg = "Instance::calcTerminalCurrents";
        msg += "Unrecognized type on boundary.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }
    else
    {
      std::string msg = "Instance::calcTerminalCurrents";
      msg += "Unrecognized boundary.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }

    if (displCurrentFlag)
    {
      bcVec[iBC].currentSum += bcVec[iBC].displCurrent;
    }
  }

  //FIX THIS:
  LeadCurrent = bcVec[1].currentSum;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::subsection_divider << std::endl
                 << "Calculated currents, etc., coming from the DD calculation:" << std::endl
                 << "  scalingVars.J0 = " << J0<<std::endl
                 << "  scalingVars.a0 = " << a0<<std::endl
                 << Xyce::subsection_divider << std::endl;
    for (int iBC=0;iBC<bcVec.size();++iBC)
    {
      Xyce::dout() << bcVec[iBC];
    }
    Xyce::dout() << Xyce::subsection_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdTerminalCurrents
//
// Purpose       : This function calculates partial derivatives associated
//                 with the Jacobian loads for the KCL equation rows.
//
// Special Notes : Originally, this work was performed in the loadJacDD
//                 function.  However, some of this information is also needed
//                 for the decoupled 2-level Newton, to calculate the
//                 terminal conductances.
//
//                 To calculate the terminal conductances, the following is
//                 needed for each electrode.
//
//                  dIdVckt - derivative of terminal current w.r.t. Vckt.
//                            This is the also the Jacobian contribution
//                            for the (KCL row, KCL col) entry of the matrix.
//
//                  dFdVckt - derivative of the RHS vector w.r.t. Vckt.
//                            This is a vector quantity, and corresponds to
//                            the (*, Vckt) column of the PDE matrix
//                            sub-block.
//
//                  dIdX  - derivative of the terminal current w.r.t. the
//                          vector of PDE solution variables. (ie not
//                          including Vckt, as that is not part of the PDE
//                          domain).  This is a vector quantity.
//                          This corresponds to the (KCL row, *) entry of
//                          the matrix, modulo dIdVckt.
//
//                  With the  "new" boundary conditions, the dFdVckt vector
//                  should only have one nonzero element, which corresponds
//                  to dIdVckt.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/03
//-----------------------------------------------------------------------------
bool Instance::pdTerminalCurrents ()
{
  std::string msg;

  double & J0 = scalingVars.J0;
  double & a0 = scalingVars.a0;

  // first calculate dIdVckt.----------------------------------------
  int iBC;
  int bcSize = bcVec.size();
  for (iBC=0; iBC < bcSize; ++iBC)
  {
    bcVec[iBC].dIdVckt = 0.0;
  } // end if BC loop

  // Now calculate dFdVckt.----------------------------------------
  // dFdVckt is a column of the matrix, which corresponds to the variable
  // Vckt, but only includes rows for internal PDE device variables.
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    int iN  = bcVec[iBC].meshIndex;
    int iNN = bcVec[iBC].neighborNode;
    int dFdVcktIndex = 0;

    // poisson term:
    bcVec[iBC].dFdVckt[dFdVcktIndex] = -scalingVars.rV0;
  }

  // now do dIdX.-------------------------------------------------
  // This is a kludge that I pretty much just copied over from
  // function loadJacKCLDDFormulation.
  // dIdX is pretty much just a KCL row from the matrix, but with only the
  // terms which pertain to internal PDE device variables included.
  // anode:
  int liOffIndex;
  int i;
  int count;

  for (iBC=0;iBC<bcSize;++iBC)
  {
    liOffIndex = 1;
    i=bcVec[iBC].meshIndex; //i is the mesh point

    if (edgeBoundarySten[i]!=1 &&  internalBoundarySten[i]!=1)
    {
      std::string msg = "Instance::pdTerminalCurrents";
      msg += "Unrecognized boundary.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }

    int iNN=bcVec[iBC].neighborNode;
    double & area = bcVec[iBC].area;
    double A0=J0*a0*area;
    std::vector<int> & colA = bcVec[iBC].li_colArray;

    double sign = ((iNN > i)?1.0:-1.0);
    double dJndV = 0.0;
    double dJpdV = 0.0;
    double dJndn = 0.0;
    double dJpdp = 0.0;

    double dJndV_nn = 0.0;
    double dJpdV_nn = 0.0;
    double dJndn_nn = 0.0;
    double dJpdp_nn = 0.0;

    // if looking right, then edge=i, if looking left, then edge=i-1.
    int edgeIndex=i;

    if (i<iNN)
    {
      edgeIndex=i;
      dJndV = dJndV1Vec[edgeIndex];
      dJpdV = dJpdV1Vec[edgeIndex];
      dJndn = dJndn1Vec[edgeIndex];
      dJpdp = dJpdp1Vec[edgeIndex];

      dJndV_nn = dJndV2Vec[edgeIndex];
      dJpdV_nn = dJpdV2Vec[edgeIndex];
      dJndn_nn = dJndn2Vec[edgeIndex];
      dJpdp_nn = dJpdp2Vec[edgeIndex];
    }
    else
    {
      edgeIndex=iNN;
      dJndV = dJndV2Vec[edgeIndex];
      dJpdV = dJpdV2Vec[edgeIndex];
      dJndn = dJndn2Vec[edgeIndex];
      dJpdp = dJpdp2Vec[edgeIndex];

      dJndV_nn = dJndV1Vec[edgeIndex];
      dJpdV_nn = dJpdV1Vec[edgeIndex];
      dJndn_nn = dJndn1Vec[edgeIndex];
      dJpdp_nn = dJpdp1Vec[edgeIndex];
    }

    // center (boundary) point:
    double Vcoef = (sign* dJndV + sign* dJpdV)*A0;
    double Ncoef = (sign* dJndn)*A0;
    double Pcoef = (sign* dJpdp)*A0;

    // neighbor point:
    double Vcoef_nn = (sign* dJndV_nn + sign* dJpdV_nn)*A0;
    double Ncoef_nn = (sign* dJndn_nn)*A0;
    double Pcoef_nn = (sign* dJpdp_nn)*A0;

    if (internalBoundarySten[i]==1)
    {
      std::string & type = bcVec[iBC].type;

      // only majority carrier goes to the boundary
      if (type=="ntype")
      {
        Pcoef=0.0;
        Pcoef_nn=0.0;
      }
      else if (type=="ptype")
      {
        Ncoef=0.0;
        Ncoef_nn=0.0;
      }
      else // oops.
      {
        std::string msg = "Instance::pdTerminalCurrents";
        msg += "Unrecognized type on boundary.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }

    // nn=i+1 for right-looking or nn=i-1 for left-looking

    count = 0;

    liOffIndex=1;
    // derivative w.r.t. V[i]:
    if (bcVec[iBC].li_colArray[liOffIndex] != -1)
    {
      bcVec[iBC].dIdXcols[count] = bcVec[iBC].li_colArray[liOffIndex];
      bcVec[iBC].dIdX[count] = Vcoef;
      ++count;
    }
    ++liOffIndex;

    // derivative w.r.t. V[nn].
    if (bcVec[iBC].li_colArray[liOffIndex] != -1)
    {
      bcVec[iBC].dIdXcols[count] = bcVec[iBC].li_colArray[liOffIndex];
      bcVec[iBC].dIdX[count] = Vcoef_nn;
      ++count;
    }
    ++liOffIndex;

    // derivative w.r.t. n[i]:
    if (bcVec[iBC].li_colArray[liOffIndex] != -1)
    {
      bcVec[iBC].dIdXcols[count] = bcVec[iBC].li_colArray[liOffIndex];
      bcVec[iBC].dIdX[count] = Ncoef;
      ++count;
    }
    ++liOffIndex;

    // derivative w.r.t. n[nn].
    if (bcVec[iBC].li_colArray[liOffIndex] != -1)
    {
      bcVec[iBC].dIdXcols[count] = bcVec[iBC].li_colArray[liOffIndex];
      bcVec[iBC].dIdX[count] = Ncoef_nn;
      ++count;
    }
    ++liOffIndex;

    // derivative w.r.t. p[i]:
    if (bcVec[iBC].li_colArray[liOffIndex] != -1)
    {
      bcVec[iBC].dIdXcols[count] = bcVec[iBC].li_colArray[liOffIndex];
      bcVec[iBC].dIdX[count] = Pcoef;
      ++count;
    }
    ++liOffIndex;

    // derivative w.r.t. p[nn].
    if (bcVec[iBC].li_colArray[liOffIndex] != -1)
    {
      bcVec[iBC].dIdXcols[count] = bcVec[iBC].li_colArray[liOffIndex];
      bcVec[iBC].dIdX[count] = Pcoef_nn;
      ++count;
    }
    ++liOffIndex;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcDXDV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/03
//-----------------------------------------------------------------------------
bool Instance::calcDXDV ()
{

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDFDV
// Purpose       : Load -dfdv into the RHS vector for the specified
//                 electrode.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/03
//-----------------------------------------------------------------------------
bool Instance::loadDFDV (int ielectrode, N_LAS_Vector * dfdvPtr)
{
  bool bsuccess = true;
  bool bs1 = true;
  N_LAS_Vector & dfdv = *(dfdvPtr);

  bcData & bc = bcVec[ielectrode];

  double coef;
  int dFdVindex = 0;

  int inode = bc.neighborNode;

  int Vrow = li_Vrowarray[inode];
  int Nrow = li_Nrowarray[inode];
  int Prow = li_Prowarray[inode];

  // load V term:
  coef = bc.dFdVckt[dFdVindex];
  dfdv[Vrow] = -coef;

  ++dFdVindex;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcConductance
//
// Purpose       : Calculates device conductances for a single electrode of
//                 the PDE device.  This function is for
//                 calculating conductances between extern circuit nodes.
//
//                 The point of calculating these quantities is to provide
//                 a lumped parameter substitute for the full device, when
//                 running 2-level Newton.
//
// Special Notes : This function is (ultimately) invoked from the nonlinear
//                 solver, as that part of the code, when running in
//                 2-level mode, knows when this information is needed.
//
//                 It is assumed that when this function is called, the
//                 "deltaX" vector contains the information:  dXdVckt, where
//                 X is solution vector variables associated with the PDE
//                 device, while Vckt is the voltage on the attached
//                 circuit node.  The reason it is in the "deltaX" vector
//                 is that it was obtained via a linear solver of the
//                 problem:
//
//                 dXdVckt = J^-1 . dFdVckt
//
//                 dFdVckt was calculated previously in pdTerminalCurrents,
//                 and J is the Jacobian.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/03
//-----------------------------------------------------------------------------
bool Instance::calcConductance (int iElectrode, const N_LAS_Vector * dxdvPtr)
{
  bool bsuccess = true;
  const N_LAS_Vector & dxdv = *dxdvPtr;

#ifdef Xyce_DEBUG_DEVICE
  char filename1[256];

  for (int ich = 0; ich < 256; ++ich)
  { filename1[ich] = 0; }

  if (getDeviceOptions().debugLevel > -1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "calcConductances  name = " << getName() << std::endl;
    Xyce::dout() << "electrode = " << bcVec[iElectrode].eName;
    Xyce::dout() << "  dIdVckt = " << bcVec[iElectrode].dIdVckt;
    Xyce::dout() << std::endl;
    Xyce::dout() << std::endl;
  }
#endif

  if (!(bcVec[iElectrode].dxdvAllocated))
  {
    bcVec[iElectrode].dxdvPtr = extData.lasSysPtr->builder().createVector();
    bcVec[iElectrode].dxdvAllocated = true;
  }

  // A linear solve should have just been performed up in the Newton
  // solver.  The result of that solve, dxdv was placed in the RHS vector.
  // dxdv is needed later, so save a copy.
  *(bcVec[iElectrode].dxdvPtr) = *(dxdvPtr);

  // doing the iElectrode Column of the condVec array.
  // This should correspond to the bcVec[iElectrode].gid column of the
  // Jacobian.

  double Gij = 0.0;
  double dIidVj = 0.0;
  double dIidVj_chain = 0.0;  // from the dot product.

  for (int iEqu=0;iEqu< numElectrodes; ++iEqu)
  {
    // conductance Gij .
    //
    // subscript i = variable, which is one of the electrode voltages.
    // subscript j = electrode
    //
    //
    //   Gij = dIi/dVj_chain          + dIi/dVj
    //
    //       = dot( dIi/dX , dX/dVj ) + dIi/dVj
    //
    //  if i != j, then the last term is zero.

    if (iElectrode != iEqu) dIidVj = 0.0;
    else                    dIidVj = bcVec[iEqu].dIdVckt;

    // load dIdX:
    extData.tmpdIdXPtr->putScalar (0.0);
    int DIDXSize = bcVec[iEqu].dIdX.size();
    for (int iDIDX=0;iDIDX<DIDXSize;++iDIDX)
    {
      int index = bcVec[iEqu].dIdXcols[iDIDX];
      double coef = bcVec[iEqu].dIdX[iDIDX];

      if (index < 0) continue;

      (*(extData.tmpdIdXPtr))[index] = coef;
    }

#ifdef Xyce_DEBUG_DEVICE
    sprintf(filename1,"dIdX%02d.txt", iEqu);
    extData.tmpdIdXPtr->writeToFile(filename1);
#endif

    // get dot product:
    dIidVj_chain = dxdv.dotProduct( *(extData.tmpdIdXPtr) );

    Gij = dIidVj_chain + dIidVj;

    condVec[iEqu][iElectrode] = Gij;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > -5 && getSolverState().debugTimeFlag)
    {
      char outstring[128];
      double Itmp = bcVec[iEqu].currentSum;
      double Vtmp = bcVec[iEqu].Vckt - bcVec[iElectrode].Vckt;
      Vtmp *= scalingVars.V0;
      double GV = Gij*Vtmp;
      for(int i=0;i<128;++i) outstring[i] = static_cast<char>(0);
      sprintf(outstring,
              "(%2d,%2d): dotPr=%12.4e G=%12.4e",
              iEqu,iElectrode,dIidVj_chain,Gij);
      Xyce::dout() << std::string(outstring) << std::endl;

      sprintf(outstring,
              "(%2d,%2d): G=%12.4e G*V=%12.4e I=%12.4e V=%12.4e",
              iEqu,iElectrode,Gij,GV,Itmp,Vtmp);
      Xyce::dout() << std::string(outstring) << std::endl;
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/18/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  updateIntermediateVars ();
  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);

  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int li_state=bcVec[iBC].li_stateC;
    staVector[li_state] = bcVec[iBC].currentSum;
  }

  // Now store the dielectric displacement in the state vector for
  // displacement current.
  int i;
  for (i = 0; i< NX-1; ++i)
  {
    double D = eSi * e0 * scalingVars.E0 * ExVec[i];
    staVector[li_stateDispl[i]] = D;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/18/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;
  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
  N_LAS_Vector & staDeriv =  *(extData.nextStaDerivVectorPtr);

  // Now get displacement current.
  for( int i = 0; i< NX-1; ++i)
  {
    displCurrent[i] = staDeriv[li_stateDispl[i]];
  }

  displCurrent[LX] = displCurrent[LX-1];

  calcTerminalCurrents ();

  return bsuccess;
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
//                 For this device, the electrostatic potential is
//                 optionally excluded from time integration norm
//                 calculations.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/07/10
//-----------------------------------------------------------------------------
bool Instance::loadDeviceMask()
{
  if (maskVarsTIAFlag_)
  {
    N_LAS_Vector * maskVectorPtr = extData.deviceMaskVectorPtr;

    for (int i=0;i<NX;++i)
    {
      int Vrow = li_Vrowarray[i];
      (*maskVectorPtr)[Vrow] = 0.0;
      (*maskVectorPtr)[Vrow] = 0.0;
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool Instance::setInitialGuess ()
{
  bool bsuccess = true;
  bool bs1 = true;

  if (variablesScaled)
  {
    bs1 = unScaleVariables ();    bsuccess = bsuccess && bs1;
  }
  bs1 = calcDensityBCs      (); bsuccess = bsuccess && bs1;
  bs1 = calcVequBCs         (); bsuccess = bsuccess && bs1;
  bs1 = calcInitialGuess    (); bsuccess = bsuccess && bs1;
  bs1 = calcMobilities      (); bsuccess = bsuccess && bs1;
  bs1 = calcLifetimes       (); bsuccess = bsuccess && bs1;
  bs1 = scaleVariables      (); bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadVecNLPoisson
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
bool Instance::loadVecNLPoisson (double * rhs)
{
  bool bsuccess = true;
  int i;
  int Vrow, Nrow, Prow;
  double coef, coef2;

  Ut = Vt/scalingVars.V0;

  // KCL equations for the two connecting terminals:
  // For the NL poisson, there is no coupling to the circuit,
  // so nothing to do here.


  // boundary conditions on the mesh:
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    // no coupling to ckt, so effectively Vckt is hardwired to zero.
    //rhs[Vrow] += VVec[i] - (bcVec[iBC].Vckt + bcVec[iBC].Vequ);
    rhs[Vrow] += VVec[i] - (bcVec[iBC].Vequ);
    rhs[Nrow] = 0.0;
    rhs[Prow] = 0.0;
  }

  // mesh points for the PDE problem:
  for (i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;

    ExtendedString semi = bulkMaterial; semi.toLower();

    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    coef  = -(VVec[i+1] - 2*VVec[i] + VVec[i-1]);
    coef /= (dxVec[i-1]*dxVec[i]);
    //coef  *= scalingVars.L0;
    coef  *= scalingVars.L0 * matSupport.getRelPerm(semi);

    double holeDens = getVoltDepHoleDens ( VminExp, VVec[i], Na);
    double elecDens = getVoltDepElecDens ( VmaxExp, VVec[i], Nd);
    coef2 = -(holeDens-elecDens+CVec[i]);

    coef += coef2;
    rhs[Vrow] += coef;

    // Now do electron, hole continuity
    rhs[Nrow] = 0.0;
    rhs[Prow] = 0.0;
  } // row loop...

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadVecDDForm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
bool Instance::loadVecDDForm (double * rhs)
{
  bool bsuccess = true;
  int i;
  int Vrow, Nrow, Prow;
  double coef, coef2;

  // KCL equations for the two connecting terminals:
  // if this is the inner loop of a multilevel Newton solve, don't do the
  // KCL-related loads.
  if ( !(getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM))
  {
    for (int iBC=0;iBC<bcVec.size();++iBC)
    {
      rhs[bcVec[iBC].lid] += bcVec[iBC].currentSum;
    }
  } // end of twoLevelNewtonCouplingMode if statement.

  // mesh points for the PDE problem:

  // boundary conditions:
  // all of these take a Dirchlet BC on voltage.
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    rhs[Vrow] += VVec[i]-bcVec[iBC].Vbc;

    if (edgeBoundarySten[i])
    {
      rhs[Nrow] = nnVec[i]-bcVec[iBC].nnbc;
      rhs[Prow] = npVec[i]-bcVec[iBC].npbc;
    }
    else if (internalBoundarySten[i])
    {
      std::string & type = bcVec[iBC].type;
      double aveDx = 0.5*(dxVec[i-1] + dxVec[i]);

      if (type=="ntype")  // boundary condition on e-, let h+ flow
      {
        rhs[Nrow] = nnVec[i]-bcVec[iBC].nnbc;
        rhs[Prow] = -(JpxVec[i]-JpxVec[i-1])/aveDx - RVec[i];
      }
      else if (type=="ptype")  // boundary condition on h+, let e- flow
      {
        rhs[Nrow] = (JnxVec[i]-JnxVec[i-1])/aveDx - RVec[i];
        rhs[Prow] = npVec[i]-bcVec[iBC].npbc;
      }
      else
      {
        std::string msg = "Instance::loadVecDDForm";
        msg += "Unrecognized type on boundary.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }
    else
    {
      std::string msg = "Instance::loadVecDDForm";
      msg += "Unrecognized stencil on boundary.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }
  }

  // interior mesh points:
  for (i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;

    ExtendedString semi = bulkMaterial; semi.toLower();

    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    coef  = -(VVec[i+1] - 2*VVec[i] + VVec[i-1]);
    coef /= (dxVec[i-1]*dxVec[i]);

    //coef *= scalingVars.L0;
    coef *= scalingVars.L0 * matSupport.getRelPerm(semi);
    coef2 = -(npVec[i]-nnVec[i]+CVec[i]);

    coef += coef2;
    rhs[Vrow] += coef;

    // Now do electron continuity
    // get electron time derivative and scale.
    double aveDx = 0.5*(dxVec[i-1] + dxVec[i]);
    rhs[Nrow] = (JnxVec[i]-JnxVec[i-1])/aveDx - RVec[i];

    // Now do hole continuity
    // get hole time derivative and scale.
    rhs[Prow] = -(JpxVec[i]-JpxVec[i-1])/aveDx - RVec[i];

  } // row loop...

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
//                 It does not bother to check them in any way, or put them
//                 in order.  It only adds them in.
//
// Special Notes : Breakpoints other than those from the photocurrent
//                 in this device are all generated by the reaction network
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints(
  std::vector<N_UTL_BreakPoint> &breakPointTimes)
{

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatNLPoisson
// Purpose       : This function performs an analytic Jacobian matrix load for
//                 the diode-pde class, for the case of solving a nonlinear
//                 poisson equation.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/18/01
//-----------------------------------------------------------------------------
bool Instance::loadMatNLPoisson (N_LAS_Matrix & mat)
{
  bool bsuccess = true;
  int Vrow, Nrow, Prow;
  int i,j;

  Ut = Vt/scalingVars.V0;
  double rUt = 1.0/Ut;

  double elecDens, holeDens;
  double dx1, dx2;
  int iBC;

  // Load the jacobian, row by row.

  // For the NL poisson option, device is not coupled to the ckt,
  // so put 1's on the diagonal.
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    mat[bcVec[iBC].lid][bcVec[iBC].lidOffset] = 1.0;
  }

  // boundary conditions on the mesh:
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    int offset1 = li_Vcolarray[i][0];
    int offset2 = li_Vcolarray[i][1];
    int offset3 = li_Vcolarray[i][2];

    if (i==0)
    {
      mat[Vrow][offset1] = 0.0;
      mat[Vrow][offset2] = 1.0;
      mat[Vrow][offset3] = 0.0;
    }
    else if (i==LX)
    {
      mat[Vrow][offset1] = 0.0;
      mat[Vrow][offset2] = 1.0;
      mat[Vrow][offset3] = 0.0;
    }
    else if (internalBoundarySten[i]==1)
    {
      mat[Vrow][offset1] = 0.0;
      mat[Vrow][offset2] = 0.0;
      mat[Vrow][offset3] = 1.0;
    }

    mat[Nrow][li_Ncolarray[i][1]] = 1.0;
    mat[Prow][li_Pcolarray[i][1]] = 1.0;
  }

  // rows associated with the PDE mesh:
  for (i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;

    ExtendedString semi = bulkMaterial; semi.toLower();

    holeDens = getVoltDepHoleDens ( VminExp   , VVec[i], Na);
    elecDens = getVoltDepElecDens ( VmaxExp   , VVec[i], Nd);

    dx1 = dxVec[i-1];
    dx2 = dxVec[i];
    //double & L0 = scalingVars.L0;
    double L0 = scalingVars.L0 * matSupport.getRelPerm(semi);

    int Vrow = li_Vrowarray[i];
    int Nrow = li_Nrowarray[i];
    int Prow = li_Prowarray[i];

    int offset1 = li_Vcolarray[i][0];
    int offset2 = li_Vcolarray[i][1];
    int offset3 = li_Vcolarray[i][2];

    mat[Vrow][offset1] = -L0/(dx1*dx2);
    mat[Vrow][offset2] = 2.0*L0/(dx1*dx2) + rUt*holeDens + rUt*elecDens;
    mat[Vrow][offset3] = -L0/(dx1*dx2);

    mat[Nrow][li_Ncolarray[i][1]] = 1.0;
    mat[Prow][li_Pcolarray[i][1]] = 1.0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatKCLDDForm
// Purpose       : Loads drift-diffusion-KCL equations into a matrix.
// Special Notes : This function is used for both old and new DAE.
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool Instance::loadMatKCLDDForm (N_LAS_Matrix & mat)
{
  int Vrow;
  int iBC;

  double & J0 = scalingVars.J0;
  double & a0 = scalingVars.a0;

  int i,j;
  int liOffIndex;
  int count  = 0;

  // rows associated with the connecting terminal KCL's:
  int li_row; // circuit node lid

  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    liOffIndex = 1;
    i=bcVec[iBC].meshIndex; //i is the mesh point

    if (edgeBoundarySten[i]!=1 &&  internalBoundarySten[i]!=1)
    {
      std::string msg = "Instance::loadMatKCLForm";
      msg += "Unrecognized boundary.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }

    int iNN=bcVec[iBC].neighborNode;
    li_row = bcVec[iBC].lid;
    double & area = bcVec[iBC].area;
    double A0=J0*a0*area;
    std::vector<int> & colA = bcVec[iBC].li_colArray;

    double sign = ((iNN > i)?1.0:-1.0);
    double dJndV = 0.0;
    double dJpdV = 0.0;
    double dJndn = 0.0;
    double dJndp = 0.0;
    double dJpdp = 0.0;
    double dJpdn = 0.0;

    double dJndV_nn = 0.0;
    double dJpdV_nn = 0.0;
    double dJndn_nn = 0.0;
    double dJndp_nn = 0.0;
    double dJpdp_nn = 0.0;
    double dJpdn_nn = 0.0;

    // if looking right, then edge=i, if looking left, then edge=i-1.
    int edgeIndex=i;

    if (i<iNN)
    {
      edgeIndex=i;
      dJndV = dJndV1Vec[edgeIndex];
      dJpdV = dJpdV1Vec[edgeIndex];
      dJndn = dJndn1Vec[edgeIndex];
      dJndp = dJndp1Vec[edgeIndex];
      dJpdp = dJpdp1Vec[edgeIndex];
      dJpdn = dJpdn1Vec[edgeIndex];

      dJndV_nn = dJndV2Vec[edgeIndex];
      dJpdV_nn = dJpdV2Vec[edgeIndex];
      dJndn_nn = dJndn2Vec[edgeIndex];
      dJndp_nn = dJndp2Vec[edgeIndex];
      dJpdp_nn = dJpdp2Vec[edgeIndex];
      dJpdn_nn = dJpdn2Vec[edgeIndex];
    }
    else
    {
      edgeIndex=iNN;
      dJndV = dJndV2Vec[edgeIndex];
      dJpdV = dJpdV2Vec[edgeIndex];
      dJndn = dJndn2Vec[edgeIndex];
      dJndp = dJndp2Vec[edgeIndex];
      dJpdp = dJpdp2Vec[edgeIndex];
      dJpdn = dJpdn2Vec[edgeIndex];

      dJndV_nn = dJndV1Vec[edgeIndex];
      dJpdV_nn = dJpdV1Vec[edgeIndex];
      dJndn_nn = dJndn1Vec[edgeIndex];
      dJndp_nn = dJndp1Vec[edgeIndex];
      dJpdp_nn = dJpdp1Vec[edgeIndex];
      dJpdn_nn = dJpdn1Vec[edgeIndex];
    }

    // center (boundary) point:
    double Vcoef = sign*(dJndV + dJpdV)*A0;
    double Ncoef = sign*(dJndn + dJpdn)*A0;
    double Pcoef = sign*(dJndp + dJpdp)*A0;

    // neighbor point:
    double Vcoef_nn = sign*(dJndV_nn + dJpdV_nn)*A0;
    double Ncoef_nn = sign*(dJndn_nn + dJpdn_nn)*A0;
    double Pcoef_nn = sign*(dJndp_nn + dJpdp_nn)*A0;

    if (internalBoundarySten[i]==1)
    {
      std::string & type = bcVec[iBC].type;

      // only majority carrier goes to the boundary
      if (type=="ntype")
      {
        Pcoef=0.0;
        Pcoef_nn=0.0;
      }
      else if (type=="ptype")
      {
        Ncoef=0.0;
        Ncoef_nn=0.0;
      }
      else // oops.
      {
        std::string msg = "Instance::loadMatKCLForm";
        msg += "Unrecognized type on boundary.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }

    // nn=i+1 for right-looking or nn=i-1 for left-looking
    // derivative w.r.t. V[i]:
    mat[li_row][colA[liOffIndex++]] += Vcoef;

    // derivative w.r.t. V[nn].
    mat[li_row][colA[liOffIndex++]] += Vcoef_nn;

    // derivative w.r.t. n[i]:
    mat[li_row][colA[liOffIndex++]] += Ncoef;

    // derivative w.r.t. n[nn].
    mat[li_row][colA[liOffIndex++]] += Ncoef_nn;

    // derivative w.r.t. p[i]:
    mat[li_row][colA[liOffIndex++]] += Pcoef;

    // derivative w.r.t. p[nn].
    mat[li_row][colA[liOffIndex++]] += Pcoef_nn;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatCktTrivial
// Purpose       : This function handles rows associated with the connecting
//                 terminal KCL's:
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/29/05
//-----------------------------------------------------------------------------
bool Instance::loadMatCktTrivial (N_LAS_Matrix & mat)
{
  bool bsuccess = true;
  bool bs1 = true;
  int j;
  int iBC;

  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    mat[bcVec[iBC].lid][bcVec[iBC].lidOffset] = 1.0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatDDForm
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
bool Instance::loadMatDDForm (N_LAS_Matrix & mat)
{
  bool bsuccess = true;
  bool bs1 = true;

  int Vrow, Nrow, Prow;

  int i,j;
  double coef;

  // set up some of the partial derivative arrays:
  bs1 = pdRecombination ();    bsuccess = bsuccess && bs1;
  bs1 = pdElectronCurrent ();  bsuccess = bsuccess && bs1;
  bs1 = pdHoleCurrent ();      bsuccess = bsuccess && bs1;
  bs1 = pdTerminalCurrents (); bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -10 && getSolverState().debugTimeFlag)
  {
    for (int i=1;i<LX;++i)
    {
      double aveDx = 0.5*(dxVec[i-1]+dxVec[i]);
      Xyce::dout()<<"\t" <<dJndp1Vec[i]/aveDx<<"\t"<<dJndp2Vec[i]/aveDx
               <<"\t" <<dJpdn1Vec[i]/aveDx<<"\t"<<dJpdn2Vec[i]/aveDx<<std::endl;
    }
  }
#endif

  if ( !(getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM))
  {
    bs1 = loadMatKCLDDForm ( mat );
    bsuccess = bsuccess && bs1;
  }
  else
  {
    bs1 = loadMatCktTrivial ( mat );
    bsuccess = bsuccess && bs1;
  }

  // boundary conditions on the mesh:
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    if (i==0)
    {
      mat[Vrow][li_Vcolarray[i][0]] = -scalingVars.rV0;
      mat[Vrow][li_Vcolarray[i][1]] =  1.0;
      mat[Vrow][li_Vcolarray[i][2]] =  0.0;

      mat[Nrow][li_Ncolarray[i][1]] = +1.0;
      mat[Prow][li_Pcolarray[i][1]] = 1.0;
    }
    else if (i==LX)
    {
      mat[Vrow][li_Vcolarray[i][0]] =  0.0;
      mat[Vrow][li_Vcolarray[i][1]] =  1.0;
      mat[Vrow][li_Vcolarray[i][2]] = -scalingVars.rV0;

      mat[Nrow][li_Ncolarray[i][1]] = +1.0;
      mat[Prow][li_Pcolarray[i][1]] = 1.0;
    }
    else if (internalBoundarySten[i]==1)
    {
      mat[Vrow][li_Vcolarray[i][0]] = -scalingVars.rV0;
      mat[Vrow][li_Vcolarray[i][2]] = 1.0;

      std::string & type = bcVec[iBC].type;

      if (type=="ntype")  // boundary condition on e-, let h+ flow
      {
        mat[Nrow][li_Ncolarray[i][1]] = +1.0;

        double aveDx = 0.5*(dxVec[i-1]+dxVec[i]);

        // derivative w.r.t. npVec[i-1]:
        //mat[Prow][li_Pcolarray[i][0]] = dJpdn1Vec[i-1]/aveDx;
        mat[Prow][li_Pcolarray[i][0]] = dJpdp1Vec[i-1]/aveDx;

        // derivative w.r.t. npVec[i  ]:
        mat[Prow][li_Pcolarray[i][1]] =
          -(dJpdp1Vec[i] - dJpdp2Vec[i-1])/aveDx - dRdpVec[i];

        // derivative w.r.t. npVec[i+1]:
        mat[Prow][li_Pcolarray[i][2]] = -dJpdp2Vec[i]/aveDx;

        // derivative w.r.t.  VVec[i-1]:
        mat[Prow][li_Pcolarray[i][3]] = (dJpdV1Vec[i-1]/aveDx);

        // derivative w.r.t.  VVec[i  ]:
        mat[Prow][li_Pcolarray[i][4]] =
          -(dJpdV1Vec[i] - dJpdV2Vec[i-1])/aveDx;

        // derivative w.r.t.  VVec[i+1]:
        mat[Prow][li_Pcolarray[i][5]] = (-dJpdV2Vec[i]/aveDx);

        // derivative w.r.t. nnVec[i  ]:
        //mat[Prow][li_Pcolarray[i][6]] = -dRdnVec[i];

        // derivative w.r.t. nnVec[i-1]:
        mat[Prow][li_Pcolarray[i][6]] = dJpdn1Vec[i-1]/aveDx;

        // derivative w.r.t. nnVec[i  ]:
        mat[Prow][li_Pcolarray[i][7]] =
          -(dJpdn1Vec[i] - dJpdn2Vec[i-1])/aveDx -dRdnVec[i];

        // derivative w.r.t. nnVec[i+1]:
        mat[Prow][li_Pcolarray[i][8]] = -dJpdn2Vec[i]/aveDx;
      }
      else if (type=="ptype")  // boundary condition on h+, let e- flow
      {
        mat[Prow][li_Pcolarray[i][1]] = 1.0;

        double aveDx = 0.5*(dxVec[i-1]+dxVec[i]);

        // derivative w.r.t. nnVec[i-1]:
        mat[Nrow][li_Ncolarray[i][0]] = -dJndn1Vec[i-1]/aveDx;

        // derivative w.r.t. nnVec[i  ]:
        mat[Nrow][li_Ncolarray[i][1]] =
          (dJndn1Vec[i] - dJndn2Vec[i-1])/aveDx - dRdnVec[i];

        // derivative w.r.t. nnVec[i+1]:
        mat[Nrow][li_Ncolarray[i][2]] = dJndn2Vec[i]/aveDx;

        // derivative w.r.t.  VVec[i-1]:
        mat[Nrow][li_Ncolarray[i][3]] = (-dJndV1Vec[i-1]/aveDx);

        // derivative w.r.t.  VVec[i  ]:
        mat[Nrow][li_Ncolarray[i][4]] =
          (dJndV1Vec[i] - dJndV2Vec[i-1])/aveDx;

        // derivative w.r.t.  VVec[i+1]:
        mat[Nrow][li_Ncolarray[i][5]] = (dJndV2Vec[i]/aveDx);

        // derivative w.r.t. npVec[i-1]:
        mat[Nrow][li_Ncolarray[i][6]] = dJndp1Vec[i-1]/aveDx;

        // derivative w.r.t. npVec[i  ]:
        mat[Nrow][li_Ncolarray[i][7]] =
          -(dJndp1Vec[i] - dJndp2Vec[i-1])/aveDx - dRdpVec[i];

        // derivative w.r.t. npVec[i+1]:
        mat[Nrow][li_Ncolarray[i][8]] = -dJndp2Vec[i]/aveDx;
      }
      else
      {
        std::string msg = "Instance::loadMatDDForm";
        msg += "Unrecognized type on boundary.";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }

    }
  }

  // load the rows associated with the PDE mesh:
  for (i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;

    ExtendedString semi = bulkMaterial; semi.toLower();

    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    // poisson equation row -------------------------------------
    double dx1 = dxVec[i-1];
    double dx2 = dxVec[i];

    //double & L0 = scalingVars.L0;
    double L0 = scalingVars.L0 * matSupport.getRelPerm(semi);

    // del^2 V elements:
    *(fVmatPtr[i][0])=(-L0/(dx1*dx2));
    *(fVmatPtr[i][1])=(2.0*L0/(dx1*dx2));
    *(fVmatPtr[i][2])=(-L0/(dx1*dx2));

    //double dfdn = q/eps;
    // electron density dependency:
    *(fVmatPtr[i][3]) = +1.0;

    // hole density dependency:
    *(fVmatPtr[i][4]) = -1.0;

    // electron continuity row -------------------------------------
    double aveDx = 0.5*(dxVec[i-1]+dxVec[i]);

    // derivative w.r.t. nnVec[i-1]:
    *(fNmatPtr[i][0]) = -dJndn1Vec[i-1]/aveDx;

    // derivative w.r.t. nnVec[i  ]:
    *(fNmatPtr[i][1]) =
      (dJndn1Vec[i] - dJndn2Vec[i-1])/aveDx - dRdnVec[i];


    // derivative w.r.t. nnVec[i+1]:
    *(fNmatPtr[i][2]) = dJndn2Vec[i]/aveDx;

    // derivative w.r.t.  VVec[i-1]:
    *(fNmatPtr[i][3]) = (-dJndV1Vec[i-1]/aveDx);

    // derivative w.r.t.  VVec[i  ]:
    *(fNmatPtr[i][4]) =
      (dJndV1Vec[i] - dJndV2Vec[i-1])/aveDx;

    // derivative w.r.t.  VVec[i+1]:
    *(fNmatPtr[i][5]) = (dJndV2Vec[i]/aveDx);

    // derivative w.r.t. npVec[i-1]:
    *(fNmatPtr[i][6]) = -dJndp1Vec[i-1]/aveDx;

    // derivative w.r.t. npVec[i  ]:
    *(fNmatPtr[i][7]) =
      (dJndp1Vec[i] - dJndp2Vec[i-1])/aveDx -dRdpVec[i];

    // derivative w.r.t. npVec[i+1]:
    *(fNmatPtr[i][8]) = dJndp2Vec[i]/aveDx;

    // hole continuity row -------------------------------------

    // derivative w.r.t. npVec[i-1]:
    *(fPmatPtr[i][0]) = dJpdp1Vec[i-1]/aveDx;

    // derivative w.r.t. npVec[i  ]:
    *(fPmatPtr[i][1]) =
      -(dJpdp1Vec[i] - dJpdp2Vec[i-1])/aveDx - dRdpVec[i];

    // derivative w.r.t. npVec[i+1]:
    *(fPmatPtr[i][2]) = -dJpdp2Vec[i]/aveDx;

    // derivative w.r.t.  VVec[i-1]:
    *(fPmatPtr[i][3]) = (dJpdV1Vec[i-1]/aveDx);

    // derivative w.r.t.  VVec[i  ]:
    *(fPmatPtr[i][4]) =
      -(dJpdV1Vec[i] - dJpdV2Vec[i-1])/aveDx;

    // derivative w.r.t.  VVec[i+1]:
    *(fPmatPtr[i][5]) = (-dJpdV2Vec[i]/aveDx);

    // derivative w.r.t. nnVec[i-1]:
    *(fPmatPtr[i][6]) = -dJpdn1Vec[i-1]/aveDx;

    // derivative w.r.t. nnVec[i  ]:
    *(fPmatPtr[i][7]) =
      (dJpdn1Vec[i] - dJpdn2Vec[i-1])/aveDx -dRdnVec[i];

    // derivative w.r.t. nnVec[i+1]:
    *(fPmatPtr[i][8]) = dJpdn2Vec[i]/aveDx;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcLifetimes
// Purpose       : This function calculates the electron and hole lifetimes
//                 and places them into the tn and tp arrays.
// Special Notes : This function assumes scaling off.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcLifetimes ()
{
  for (int i=0;i<NX;++i)
  {
    tnVec[i] = matSupport.calcLt (false, fabs(CVec[i]));
    tpVec[i] = matSupport.calcLt (true , fabs(CVec[i]));
  }
  return true;
}

#if 0
//-----------------------------------------------------------------------------
// Function      : Instance::calcMobilities
// Purpose       : This function calculates the electron and hole mobilities
//                 and places them into the un and up arrays.
// Special Notes : This function assumes scaling off.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcMobilities ()
{
  bool bsuccess = true;
  int i;

  MobInfo<double> ci;
  ci.mobModelName = mobModelName;
  ci.materialName = bulkMaterial;

  for (i=0;i<NX;++i)
  {
    ci.T = Temp;  // scaling not worked out for temperature...

    ci.N = fabs(CVec[i]);
    ci.N *= ((variablesScaled)?scalingVars.C0:1.0);

    if (ci.N == 0.0) ci.N = 1.0;  // avoid nan's.

    ci.p = npVec[i] *((variablesScaled)?scalingVars.C0:1.0);
    ci.n = nnVec[i] *((variablesScaled)?scalingVars.C0:1.0);

    // electron mobility:
    ci.holeFlag = false;
    unVec[i] = matSupport.calcMob(ci);
    unVec[i] /= ((variablesScaled)?scalingVars.u0:1.0);

    // hole mobility:
    ci.holeFlag = true;
    upVec[i] = matSupport.calcMob(ci);
    upVec[i] /= ((variablesScaled)?scalingVars.u0:1.0);
  }

  // Now do edge mobilities:
  for (i=0;i<LX;++i)
  {
    // option 1
    if (ci.mobModelName != "carr")
    {
      unE_Vec[i] = (unVec[i+1]+unVec[i])/2.0;
      upE_Vec[i] = (upVec[i+1]+upVec[i])/2.0;
    }
    // option 2
    else if(ci.mobModelName == "carr")
    {
      ci.N = (fabs(CVec[i+1])+fabs(CVec[i]))*0.5;
      ci.N *= ((variablesScaled)?scalingVars.C0:1.0);

      if (ci.N == 0.0) ci.N = 1.0;

      // for the carrier densities, do a "product average".
      ci.n = pow((fabs(nnVec[i+1])*fabs(nnVec[i])),0.5)*
             ((variablesScaled)?scalingVars.C0:1.0);

      ci.p = pow((fabs(npVec[i+1])*fabs(npVec[i])),0.5)*
             ((variablesScaled)?scalingVars.C0:1.0);

      //electron mobility
      ci.holeFlag = false;
      unE_Vec[i] = matSupport.calcMob(ci);
      unE_Vec[i] /= ((variablesScaled)?scalingVars.u0:1.0);

      // hole mobility
      ci.holeFlag = true;
      upE_Vec[i] = matSupport.calcMob(ci);
      upE_Vec[i] /= ((variablesScaled)?scalingVars.u0:1.0);
    }
  }

  return bsuccess;
}
#else
//-----------------------------------------------------------------------------
// Function      : Instance::calcMobilities
// Purpose       : This function calculates the electron and hole mobilities
//                 and places them into the un and up arrays.
//
// Special Notes : The mobility functions assume that scaling is off, so
//                 this function unscales and rescales things as needed.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcMobilities ()
{
  bool bsuccess = true;
  int i;

  MobInfo<pdeFadType> mi;
  mi.mobModelName = mobModelName;
  mi.materialName = bulkMaterial;
  mi.fieldDependent = fieldDependentMobility;

  // Now do edge mobilities:
  for (i=0;i<LX;++i)
  {
    // possibly these doping related quantities should be determined
    // using splines instead.
    mi.N = (fabs(CVec[i+1])+fabs(CVec[i]))*0.5;
    mi.N *= ((variablesScaled)?scalingVars.C0:1.0);

    mi.Na = (fabs(CacceptorVec[i])+fabs(CacceptorVec[i+1]))*0.5
            *((variablesScaled)?scalingVars.C0:1.0);

    mi.Nd = (fabs(CdonorVec[i])+fabs(CdonorVec[i+1]))*0.5
            *((variablesScaled)?scalingVars.C0:1.0);

    if (mi.N == 0.0) mi.N = 1.0;

    pdeFadType v1=VVec[i];
    pdeFadType v2=VVec[i+1];
    pdeFadType n1=nnVec[i];
    pdeFadType n2=nnVec[i+1];
    pdeFadType p1=npVec[i];
    pdeFadType p2=npVec[i+1];

    v1.diff(0,6);
    v2.diff(1,6);
    n1.diff(2,6);
    n2.diff(3,6);
    p1.diff(4,6);
    p2.diff(5,6);

    pdeFadType Efield = (-(v2-v1)/dxVec[i]);
#if 0
    // for the carrier densities, do a "product average".
    mi.n = pow(fabs(n2*n1),0.5)*((variablesScaled)?scalingVars.C0:1.0);
    mi.p = pow(fabs(p2*p1),0.5)*((variablesScaled)?scalingVars.C0:1.0);
#else
    // this is the most consistent, as it relies on the SG approximation for the midpoint density.
    mi.n = fabs(nMidpoint(n1,n2,Efield,dxVec[i],-1))*((variablesScaled)?scalingVars.C0:1.0);
    mi.p = fabs(nMidpoint(p1,p2,Efield,dxVec[i],+1))*((variablesScaled)?scalingVars.C0:1.0);
#endif
    mi.epar = fabs(Efield)*((variablesScaled)?scalingVars.E0:1.0);

    //electron mobility
    mi.holeFlag = false;
    unE_Vec[i] = matSupport.calcMob(mi);
    unE_Vec[i] /= ((variablesScaled)?scalingVars.u0:1.0);

    // hole mobility
    mi.holeFlag = true;
    upE_Vec[i] = matSupport.calcMob(mi);
    upE_Vec[i] /= ((variablesScaled)?scalingVars.u0:1.0);

#if 0
    if( mi.fieldDependent )
    {
      Xyce::dout() << "\tfieldDep=true";
    }
    else
    {
      Xyce::dout() << "\tfieldDep=false";
    }
    Xyce::dout() << "\tefield ="<<mi.epar;
    Xyce::dout() << "\tn ="<<mi.n;
    Xyce::dout() << "\tp ="<<mi.p;
    Xyce::dout() << "\tNd ="<<mi.Nd;
    Xyce::dout() << "\tNa ="<<mi.Na;
#endif

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > -10 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "\tunE["<<i<<"]="<<unE_Vec[i];
      Xyce::dout() << "\tupE["<<i<<"]="<<upE_Vec[i];
      Xyce::dout() << std::endl;
    }
#endif
  }

  return bsuccess;
}
#endif

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
//
// Special Notes : This won't quite work yet, because Ni and the bandgap
//                 functions (up in the material support class) are
//                 not yet really temperature dependent.
//
//     Things that change with temperature:
//
//        thermal voltage (Vt)
//        scaling variable, scalingVars.V0 = Vt
//        intrinsic concentration, Ni.
//        bandgap, Eg.
//        mobilities (updated in real time, so no need to change here).
//        density boundary conditions (depend on Ni)
//        Vequ (may depend on Ni and Eg).
//        defect reactions
//        other??
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature(const double & temp_tmp)
{
  bool bsuccess = true;
  bool bs1 = true;

  if (indicesSetup_) // instance block constructor sets this flag,
    // but default constructor does not.  If it is set,
    // then the device is ready to process this function,
    // but not otherwise.
  {
    Temp = temp_tmp;

    // first un-scale everything, if neccessary:
    if (variablesScaled)
    {
      bs1 = unScaleVariables ();    bsuccess = bsuccess && bs1;
    }

    bs1 = setupMiscConstants ();    bsuccess = bsuccess && bs1;
    bs1 = setupScalingVars ();      bsuccess = bsuccess && bs1;

    bs1 = calcDensityBCs       ();  bsuccess = bsuccess && bs1;
    bs1 = calcVequBCs          ();  bsuccess = bsuccess && bs1;
    bs1 = calcMobilities       ();  bsuccess = bsuccess && bs1;

    if (!variablesScaled)
    {
      bs1 = scaleVariables ();      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcVoltDepDensities
// Purpose       : This function calculates electron and hole densities,
//                 based on the electrostatic potential.  It is only to be
//                 called during the initialization phase.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcVoltDepDensities ()
{
  Ut = Vt/scalingVars.V0;

  for (int i=0;i<NX;++i)
  {
    npVec[i] = getVoltDepHoleDens(VminExp , VVec[i], Na);
    nnVec[i] = getVoltDepElecDens(VmaxExp , VVec[i], Nd);
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupDopingProfile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/28/07
//-----------------------------------------------------------------------------
bool Instance::setupDopingProfile ()
{
  bool bsuccess (false);
  bool fromFile (false);
  int i;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::setupDopingProfile\n";
  }
#endif

  if ( dopingFileName != "NOFILE" )
  {
    DopeInfo::readDopingFile (dopingFileName,
                                    xloc_ndope_vec,
                                    ndope_vec, y2_ndope_vec,
                                    pdope_vec, y2_pdope_vec,
                                    devSupport );

    xloc_pdope_vec.clear();
    xloc_pdope_vec.resize( xloc_ndope_vec.size(), 0.0);
    xloc_pdope_vec = xloc_ndope_vec ;
    bsuccess=true;
    fromFile=true;
  }
  else if ( ( ( ndopeFileName != "NOFILE" ) && ( pdopeFileName != "NOFILE") ) )
  {
    DopeInfo::readDopingFile (ndopeFileName,
                                    xloc_ndope_vec, ndope_vec, y2_ndope_vec, devSupport);
    DopeInfo::readDopingFile (pdopeFileName,
                                    xloc_pdope_vec, pdope_vec, y2_pdope_vec, devSupport);

    bsuccess=true;
    fromFile=true;
  }
  else
  {
    bsuccess = calcDopingProfile ();
  }


  // use the N and P dopants to create the C vector.
  if (fromFile)
  {
    Na = 0.0;
    Nd = 0.0;
    for (i=0;i<NX;++i)
    {
      double xtmp = xVec[i];
      double ndopeDopeValue(0.0), pdopeDopeValue(0.0);
      devSupport.splint(xloc_ndope_vec,
                        ndope_vec, y2_ndope_vec, xtmp, ndopeDopeValue);
      devSupport.splint(xloc_pdope_vec,
                        pdope_vec, y2_pdope_vec, xtmp, pdopeDopeValue);
      CVec[i] = ndopeDopeValue-pdopeDopeValue;

      if (Na > CVec[i]) Na = CVec[i];
      if (Nd < CVec[i]) Nd = CVec[i];
    }
    Na = fabs(Na);
    Nd = fabs(Nd);
  }

  // now that we have the C vector, loop over the boundary
  // conditions and determine if n-type or p-type.

  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    if (CVec[i] > 0.0)
    {
      bcVec[iBC].type = "ntype";
    }
    else
    {
      bcVec[iBC].type = "ptype";
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Na = " << Na << std::endl;
    Xyce::dout() << "Nd = " << Nd << std::endl;
    for (i=0;i<NX;++i)
    {
      Xyce::dout() << "x[" << i << "] = " << xVec[i] << "\t";
      Xyce::dout() << "C[" << i << "] = " << CVec[i] << std::endl;
    }

    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcDopingProfile
//
// Purpose       : This function sets up the initial doping profile.
//
// Special Notes : 03/31/03. This function is being modified to handle a
//                 more general doping specification.  The old way of
//                 specifying doping, which assumes a PN junction, with Na
//                 and Nd may eventually be phased out.  For now, both
//                 methods of specification are supported, with the new
//                 style over-riding the old, in the event that both are
//                 specified.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcDopingProfile ()
{
  bool bsuccess = true;
  double midpoint;
  int i;

  // Check which style of doping specification to use.  If the dopeInfoMap
  // is empty, then assume the old method.  If not, use the dopeInfoMap.

  if (dopeInfoMap.empty ())
  {
    // first setup, check the graded junction parameters:
    if (gradedJunctionFlag)
    {
      // if junction width was not specified, set to 1/10 the diode width.
      if (!given("WJ"))
        WJ = 0.1 * width;

      midpoint = width/2.0;
      XL = midpoint - WJ/2.0;
      XR = midpoint + WJ/2.0;
    }

    for (i=0;i<NX;++i)
    {
      if (gradedJunctionFlag)
      {
        if (xVec[i] <= XL) CVec[i] = +Nd;
        else if (xVec[i]>XL && xVec[i]<XR)
          CVec[i] = Nd-(Na+Nd)/(XR-XL)*(xVec[i]-XL);
        else CVec[i] = -Na;
      }
      else
      {
        if (xVec[i] < xVec[LX]/2.0) CVec[i] = Nd;
        else CVec[i] = -Na;
      }
    }
  }
  else
  {
    // loop over the dope info map, and sum contributions from each
    // doping entity into the total doping array, CVec.
    std::map<std::string, DopeInfo *>::iterator iter;
    std::map<std::string, DopeInfo *>::iterator start = dopeInfoMap.begin();
    std::map<std::string, DopeInfo *>::iterator end   = dopeInfoMap.end();

    for ( iter = start; iter != end; ++iter )
    {
      DopeInfo & di = *(iter->second);
      di.setupInfo(CVec,CdonorVec,CacceptorVec,xVec,devSupport);
    }

    Na = 0.0;
    Nd = 0.0;
    for (i=0;i<NX;++i)
    {
      if (Na > CVec[i]) Na = CVec[i];
      if (Nd < CVec[i]) Nd = CVec[i];
    }
    Na = fabs(Na);
    Nd = fabs(Nd);

  } // if statement

  if (Na == 0.0 || Nd == 0.0)
  {
    // Error in given expression.
    char tmpChar[128];
    sprintf(tmpChar, "Mistake in doping. Na=%12.4e  Nd=%12.4e\n",Na,Nd);
    std::string msg(tmpChar);
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMesh
// Purpose       : This function sets up the mesh.  Should only be called once.
// Special Notes : For now this just does a simple 1D uniform mesh.  Later,
//                 this will perhaps be more complex.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::setupMesh ()
{
  // set up the mesh:
  double dx_tmp = width/(static_cast<double>(LX));

  for (int i=0;i<NX;++i)
  {
    xVec[i] = static_cast<double>(i)*dx_tmp;
  }

  for (int i=0;i<LX;++i)
  {
    dxVec[i] = xVec[i+1]-xVec[i];
  }
  dxVec[LX] = dxVec[LX-1];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    for (int i=0;i<NX;++i)
    {
      Xyce::dout() << "x["<<i<<"] = " << xVec[i];
      Xyce::dout() << "\tdx["<<i<<"] = " << dxVec[i];
      Xyce::dout() << std::endl;
    }
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMiscConstants
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/20/03
//-----------------------------------------------------------------------------
bool Instance::setupMiscConstants ()
{
  if (useOldNi)
  {
    Ni = matSupport.getNi_old (bulkMaterial, Temp); // this is not accurate
  }
  else
  {
    Ni = matSupport.getNi (bulkMaterial, Temp);
  }
  Vt = kb*Temp/charge;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupScalingVars
// Purpose       : This function sets up scaling variables.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/28/01
//-----------------------------------------------------------------------------
bool Instance::setupScalingVars ()
{
  bool bsuccess = true;

  // just to make sure:
  Vt = kb*Temp/charge;
  Vbi = Vt*log(Na*Nd/(Ni*Ni));

  if (given("X0"))
    scalingVars.x0 = x0_user;
  else
    scalingVars.x0  = width;// distance scaling (cm)

  // For the 1D device, cross-sectional area is kind of a weird concept.
  // For the equations within the device, area really doesn't factor into
  // the discretization.  It is only important at the electrodes, for
  // calculating the current coming out of the device.
  scalingVars.a0 = scalingVars.x0*scalingVars.x0;

  scalingVars.T0  = Temp;// temperature scaling (K)  (not really used)

  // electrostatic potential scaling (V)
  scalingVars.V0 = Vt;
  scalingVars.rV0 = 1.0/scalingVars.V0;

  // concentration scaling (cm^-3);
  if (given("C0"))
  {
    scalingVars.C0 = C0_user;
  }
  else if (scaleDensityToMaxDoping_)  // this is the default
  {
    if (Na >= Nd) scalingVars.C0  = Na;
    else          scalingVars.C0  = Nd;

    scalingVars.C0 *= densityScalarFraction_; // 1.0e-2 by default
  }
  else
  {
    scalingVars.C0 = 1.0e+17;
  }

  if (given("t0"))
  {
    scalingVars.t0 = t0_user;
    scalingVars.D0  = (scalingVars.x0*scalingVars.x0)/scalingVars.t0;
  }
  else
  {
    // diffusion coefficient scaling (cm^2/s)
    scalingVars.D0  = 35.0;

    // time scaling (s)
    scalingVars.t0  = (scalingVars.x0*scalingVars.x0)/scalingVars.D0;
  }

  // mobility coefficient scaling (cm^2/V/s)
  scalingVars.u0  = scalingVars.D0/scalingVars.V0;

  // recombination rate scaling (cm^-3/s)
  scalingVars.R0  = scalingVars.D0*scalingVars.C0/(scalingVars.x0*scalingVars.x0);
  scalingVars.rR0 = 1.0/scalingVars.R0;

  // electric field scaling (V/cm)
  scalingVars.E0  = scalingVars.V0/scalingVars.x0;

  // particle flux scaling (cm^-2/s)
  scalingVars.F0  = scalingVars.D0*scalingVars.C0/scalingVars.x0;

  // current density scaling (A/cm^2)
  scalingVars.J0  = charge*scalingVars.D0*scalingVars.C0/scalingVars.x0;

  // Laplacian scaling constant
  //scalingVars.L0  = scalingVars.V0*eps/(charge*scalingVars.x0*scalingVars.x0*scalingVars.C0);
  scalingVars.L0  = scalingVars.V0*e0/(charge*scalingVars.x0*scalingVars.x0*scalingVars.C0);  // Laplacian scaling constant

  // rate constant scaling.   k0 = 1/(C0*t0) = cm^3/sec
  scalingVars.rk0 = scalingVars.C0*scalingVars.t0;
  scalingVars.rt0 = 1.0/scalingVars.t0;
  scalingVars.k0 = 1.0/scalingVars.rk0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -2)
  {
    Xyce::dout() << "scalingVars.x0 = " << scalingVars.x0 << std::endl;
    Xyce::dout() << "scalingVars.a0 = " << scalingVars.a0 << std::endl;
    Xyce::dout() << "scalingVars.T0 = " << scalingVars.T0 << std::endl;
    Xyce::dout() << "scalingVars.V0 = " << scalingVars.V0 << std::endl;
    Xyce::dout() << "scalingVars.C0 = " << scalingVars.C0 << std::endl;
    Xyce::dout() << "scalingVars.D0 = " << scalingVars.D0 << std::endl;
    Xyce::dout() << "scalingVars.u0 = " << scalingVars.u0 << std::endl;
    Xyce::dout() << "scalingVars.R0 = " << scalingVars.R0 << std::endl;
    Xyce::dout() << "scalingVars.t0 = " << scalingVars.t0 << std::endl;
    Xyce::dout() << "scalingVars.E0 = " << scalingVars.E0 << std::endl;
    Xyce::dout() << "scalingVars.F0 = " << scalingVars.F0 << std::endl;
    Xyce::dout() << "scalingVars.J0 = " << scalingVars.J0 << std::endl;
    Xyce::dout() << "scalingVars.L0 = " << scalingVars.L0 << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::scaleVariables
//
// Purpose       : This function performs scaling on all the relevant variables.
//
// Special Notes : It should only be called at the end of the initial setup.
//                 Calculations done during the course of the calculation are
//                 performed with the assumption of scaling.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/28/01
//-----------------------------------------------------------------------------
bool Instance::scaleVariables ()
{
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  Na /= scalingVars.C0;
  Nd /= scalingVars.C0;
  Ni /= scalingVars.C0;

  int i;
  for (i=0;i<bcVec.size();++i)
  {
    bcVec[i].Vbc  /= scalingVars.V0;
    bcVec[i].Vequ /= scalingVars.V0;
    bcVec[i].nnbc /= scalingVars.C0;
    bcVec[i].npbc /= scalingVars.C0;

    bcVec[i].area /= scalingVars.a0;
  }
  area /= scalingVars.a0;

  VminExp    /= scalingVars.V0;
  VmaxExp    /= scalingVars.V0;
  Vbi        /= scalingVars.V0;

  maxVoltDelta /= scalingVars.V0;

  for (i=0;i<NX;++i)
  {
    nnVec[i] /= scalingVars.C0;
    npVec[i] /= scalingVars.C0;
    CVec[i]  /= scalingVars.C0;
    CdonorVec[i] /= scalingVars.C0;
    CacceptorVec[i] /= scalingVars.C0;
    VVec[i]  /= scalingVars.V0;
    //unVec[i] /= scalingVars.u0;
    //upVec[i] /= scalingVars.u0;
    tnVec[i] /= scalingVars.t0;
    tpVec[i] /= scalingVars.t0;
    xVec[i]  /= scalingVars.x0;
    dxVec[i] /= scalingVars.x0;

    (*solVectorPtr)[li_Vrowarray[i]] = VVec[i];
    (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
    (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
  }

  variablesScaled = true;

  return true;

}

//-----------------------------------------------------------------------------
// Function      : Instance::unScaleVariables
// Purpose       : This function is the inverse of scaleVariables.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/20/03
//-----------------------------------------------------------------------------
bool Instance::unScaleVariables ()
{
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  Na *= scalingVars.C0;
  Nd *= scalingVars.C0;
  Ni *= scalingVars.C0;

  int i;
  for (i=0;i<bcVec.size();++i)
  {
    bcVec[i].Vbc  *= scalingVars.V0;
    bcVec[i].Vequ *= scalingVars.V0;
    bcVec[i].nnbc *= scalingVars.C0;
    bcVec[i].npbc *= scalingVars.C0;

    bcVec[i].area *= scalingVars.a0;
  }
  area *= scalingVars.a0;

  VminExp    *= scalingVars.V0;
  VmaxExp    *= scalingVars.V0;
  Vbi        *= scalingVars.V0;

  maxVoltDelta *= scalingVars.V0;

  for (i=0;i<NX;++i)
  {
    nnVec[i] *= scalingVars.C0;
    npVec[i] *= scalingVars.C0;
    CVec[i]  *= scalingVars.C0;
    CdonorVec[i] *= scalingVars.C0;
    CacceptorVec[i] *= scalingVars.C0;
    VVec[i]  *= scalingVars.V0;
    //unVec[i] *= scalingVars.u0;
    //upVec[i] *= scalingVars.u0;
    tnVec[i] *= scalingVars.t0;
    tpVec[i] *= scalingVars.t0;
    xVec[i]  *= scalingVars.x0;
    dxVec[i] *= scalingVars.x0;

    (*solVectorPtr)[li_Vrowarray[i]] = VVec[i];
    (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
    (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
  }

  variablesScaled = false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcInitialGuess
// Purpose       : This function calculates the initial e-, h+ densties
//                 and the intial voltage.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcInitialGuess ()
{
  bool bsuccess = true;
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  // set up an initial guess for nn and np, based on the doping profile,
  // and the equilibrium density expressions.  Place these in the
  // solution vector.

  double tmp;
  double Ci;
  double Cisq, Nisq;
  for (int i=0;i<NX;++i)
  {
    Ci = CVec[i];
    Cisq = Ci*Ci;
    Nisq = Ni*Ni;  // Ni is the intrinsic concentration

    // equilibrium electron concentration:
    tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
    nnVec[i] = ((Ci>=0)?(tmp):(0.0)) + ((Ci<0)?(Nisq/tmp):(0.0));

    // equilibrium hole concentration:
    tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
    npVec[i] = ((Ci<=0)?(tmp):(0.0)) + ((Ci>0)?(Nisq/tmp):(0.0));
  }

  // set up initial guess for V, place in solution vector
  double Vmax = -1.0e99;
  double Vmin = +1.0e99;
  for (int i=0;i<NX;++i)
  {
    // the doping is n-type.
    if (nnVec[i]>=npVec[i])
    {
      VVec[i] = + Vt * log(nnVec[i]/Ni);
    }
    // the doping is p-type.
    else
    {
      VVec[i] = - Vt * log(npVec[i]/Ni);
    }

    if (Vmax < VVec[i]) Vmax = VVec[i];
    if (Vmin > VVec[i]) Vmin = VVec[i];
  }

  // get the maximum and minimum potentials.
  VmaxExp = -1.0e99;
  VminExp = +1.0e99;

  for (int i=0;i<NX;++i)
  {
    if (VmaxExp < VVec[i]) VmaxExp = VVec[i];
    if (VminExp > VVec[i]) VminExp = VVec[i];
  }

  for (int i=0;i<NX;++i)
  {
    (*solVectorPtr)[li_Vrowarray[i]] = VVec[i];
    (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
    (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcVequBCs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/17/03
//-----------------------------------------------------------------------------
bool Instance::calcVequBCs ()
{
  bool bsuccess = true;

  Vt = kb*Temp/charge;
  Vbi = Vt*log(Na*Nd/(Ni*Ni));

  double VminBC =+1.0e+99;
  double VmaxBC =-1.0e+99;

  int bcSize=bcVec.size();
  for (int i=0;i<bcSize;++i)
  {
    int mIndex = bcVec[i].meshIndex;
    double Ci = CVec[mIndex];
    double Cisq = Ci*Ci;
    double Nisq = Ni*Ni;  // Ni is the intrinsic concentration
    double tmp, nnTmp, npTmp;

    // equilibrium electron concentration:
    tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
    nnTmp = ((Ci>=0)?(tmp):(0.0)) + ((Ci<0)?(Nisq/tmp):(0.0));

    // equilibrium hole concentration:
    tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
    npTmp = ((Ci<=0)?(tmp):(0.0)) + ((Ci>0)?(Nisq/tmp):(0.0));

    ExtendedString semi = bulkMaterial;
    semi.toLower();
    ExtendedString mater = bcVec[i].material;
    mater.toLower();

    if (bcVec[i].VequGiven != 1)
    {
      if (mater=="neutral")
      {
        // the doping is n-type.
        if (Ci>0)
        {
          bcVec[i].Vequ = + Vt * log(nnTmp/Ni);
        }
        else        // the doping is p-type.
        {
          bcVec[i].Vequ = - Vt * log(npTmp/Ni);
        }
      }
      else // this electrode is a schottky barrier.
      {
        // the doping is n-type.
        if (Ci>0)
        {
          bcVec[i].Vequ = + matSupport.workfunc(mater)
                          - matSupport.affin(semi)
                          - 0.5 * matSupport.bandgap(semi, Temp)
                          + 2.0 * Vt * log(nnTmp/Ni);
        }
        else        // the doping is p-type.
        {
          bcVec[i].Vequ = + matSupport.workfunc(mater)
                          - matSupport.affin(semi)
                          - 0.5 * matSupport.bandgap(semi, Temp)
                          - 2.0 * Vt * log(npTmp/Ni);
        }
      }
    }

    if (VminBC > bcVec[i].Vequ) VminBC = bcVec[i].Vequ;
    if (VmaxBC < bcVec[i].Vequ) VmaxBC = bcVec[i].Vequ;
  }

  VoltageOffset_ = -VminBC;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcDensityBCs
// Purpose       : This function sets up the boundary condition variables
//                 for each electrode.
//
// Special Notes : This function is similar to calcBoundaryConditions, but
//                 this one only calculates BC's on N and P.  Since these
//                 never change, they only need to be calculated once.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/03
//-----------------------------------------------------------------------------
bool Instance::calcDensityBCs ()
{
  bool bsuccess = true;

  NnMax = -1.0e+99;
  NpMax = -1.0e+99;

  NnMin = +1.0e+99;
  NpMin = +1.0e+99;

  // This density boundary condition is from Selberherr,
  // enforcing thermal equilibrium and
  // vanishing space charge at ohmic contacts
  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    int i1 = bcVec[iBC].meshIndex;
    bcVec[iBC].nnbc = 0.5*(sqrt(CVec[i1]*CVec[i1]+4*Ni*Ni)+CVec[i1]);
    bcVec[iBC].npbc = 0.5*(sqrt(CVec[i1]*CVec[i1]+4*Ni*Ni)-CVec[i1]);

    if (NnMax < bcVec[iBC].nnbc) NnMax = bcVec[iBC].nnbc;
    if (NpMax < bcVec[iBC].npbc) NpMax = bcVec[iBC].npbc;
  }

  if (NnMin <= 0) NnMin = 1.56269e-9;  // just a guess.
  if (NpMin <= 0) NpMin = 1.56269e-9;  // just a guess.

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcBoundaryConditions
// Purpose       : This function sets up the boundary condition variables
//                 for each electrode.
//
// Special Notes : If a continuation method is being used, a good parameter
//                 to vary is the voltage boundary condition  on the
//                 electrodes.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool Instance::calcBoundaryConditions ()
{
  bool bsuccess = true;
  int iBC;

  int bcSize=bcVec.size();
  if (getSolverState().PDEcontinuationFlag)
  {
    for (iBC=0;iBC<bcSize;++iBC)
    {
      bcVec[iBC].Vbc = bcVec[iBC].Vckt_ramp + bcVec[iBC].Vequ;
    }
  }
  else
  {
    for (iBC=0;iBC<bcSize;++iBC)
    {
      bcVec[iBC].Vbc = (bcVec[iBC].Vckt + bcVec[iBC].Vequ);
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::obtainNodeVoltages.
//
// Purpose       : This function obtains the nodal voltages from the
//                 solution vector, which are applied as boundary
//                 conditions on the electrodes.
//
// Special Notes : This was originally part of function obtainSolution, but
//                 is needed also by function enableContinuation.  So I've
//                 put it in one place.
//
//                 If voltage limiting is turned on, this is the function
//                 in which to apply it.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/13/02
//-----------------------------------------------------------------------------
bool Instance::obtainNodeVoltages ()
{
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    bcVec[iBC].Vckt = (*solVectorPtr)[bcVec[iBC].lid];
    bcVec[iBC].Vckt /= scalingVars.V0;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::applyVoltageLimiting
//
// Purpose       : if voltage limiting is turned on, this function
//                 applies it to the Vckt values.
//
// Special Notes : This is only really set up to work when the 2-level
//                 Newton is being used.
//
//                 For now, this is just a test capability.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/15/02
//-----------------------------------------------------------------------------
bool Instance::applyVoltageLimiting ()
{
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    double v1     = bcVec[iBC].Vckt * scalingVars.V0;
    double v1_old = bcVec[iBC].Vckt_old * scalingVars.V0;
    double delV1 = v1 - v1_old;

    if ( delV1 > 1.25 )  v1 = v1_old + 1.25;

    if ( delV1 < -0.75) v1 = v1_old - 0.75;

    bcVec[iBC].Vckt       = v1/scalingVars.V0;
    bcVec[iBC].Vckt_final = v1/scalingVars.V0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::obtainSolution
// Purpose       : This function extracts V, nn, and np from the solution
//                 vector and copies them into local arrays.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::obtainSolution ()
{
  bool bsuccess = true;
  bool bs1 = true;
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  // First get the two circuit node voltages:
  bsuccess = obtainNodeVoltages ();

  // set up the V solution array:
  int i;
  for (i=0;i<NX;++i)
  {
    VVec[i] = (*solVectorPtr)[li_Vrowarray[i]];
  }

  // If the previous solution is from the nonlinear Poisson solution,
  // then calculate what the electron and hole densities must be, and
  // place them into the solution vector.

  // If we are past the nonlinear Poisson phase, then simply obtain
  // nn and np from the solution vector and move on.

  if (getSolverState().dcopFlag && getSolverState().doubleDCOPStep==0)
  {
    calcVoltDepDensities ();

    for (i=0;i<NX;++i)
    {
      (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
      (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
    }
  }
  else
  {
    for (i=0;i<NX;++i)
    {
      nnVec[i] = (*solVectorPtr)[li_Nrowarray[i]];

#ifdef Xyce_PDE_DENSITY_CONSTRAINT
      if (nnVec[i] < 0.0) nnVec[i] = 0.0;
#endif

      npVec[i] = (*solVectorPtr)[li_Prowarray[i]];

#ifdef Xyce_PDE_DENSITY_CONSTRAINT
      if (npVec[i] < 0.0) npVec[i] = 0.0;
#endif
    }

    // now set boundary conditions:
    // if the circuit is coupled to the PDE device, then bc's
    // must be updated everytime.
    //
    // If the circuit and PDE device are not coupled, then the
    // circuit node voltages can be considered constant, and the
    // BC's only need updating at the first Newton step.
    if ( !(getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM))
    {
      bs1 = calcBoundaryConditions (); bsuccess = bsuccess && bs1;
    }
    else
    {
      if (getSolverState().newtonIter == 0)
      {
        bs1 = calcBoundaryConditions (); bsuccess = bsuccess && bs1;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/22/03
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles ()
{
  bool bsuccess = true;
  bool bs1 = true;
  bool skipOutput = false;

  // usually, don't bother outputting nonlinear Poisson result.
  if (equationSet == 0 && !(outputNLPoisson))  return bsuccess;

  // If using output interval, check if enough time has passed to do
  // another output.  (only applies for transient - not DCOP).
  if ( !(getSolverState().dcopFlag) &&
       !(getSolverState().forceFinalOutput) &&
       outputIntervalGiven)
  {
    double outMult = static_cast<double> (outputIndex);
    double nextOutputTime = outMult * outputInterval;

    if (nextOutputTime > getSolverState().currTime)
    {
      skipOutput = true;
    }
  }

  // If this is a "forced" final output, make sure that it didn't already output.
  // This can happen if the output interval is an exact multiple of the
  // total simulation time.
  if (getSolverState().forceFinalOutput &&
      getSolverState().currTime==lastOutputTime) skipOutput=true;

  if (skipOutput) return bsuccess;
  ++outputIndex;
  lastOutputTime = getSolverState().currTime;

#ifdef Xyce_DEBUG_DEVICE
  Xyce::dout() << std::endl << "Doing an output at time = " << getSolverState().currTime << std::endl;
#endif

  if (tecplotLevel > 0) {bs1 = outputTecplot (); bsuccess = bsuccess && bs1;}
  if (sgplotLevel  > 0) {bs1 = outputSgplot  (); bsuccess = bsuccess && bs1;}

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputTecplot
// Purpose       : This function outputs a file which is easily plottable
//                 by tecplot.  Simply run tecplot "filename.dat" <return>
//
//
// Special Notes : This file can also be plotted using gnuplot.  If the
//                 name of the file is "Z1DIODE_000.dat", plot inside of
//                 gnuplot using:
//
//                 plot "Z1DIODE_000.dat" using 1:3 w l
//
//                 or, if you want a log plot, for the doping:
//
//                 plot "Z1DIODE_000.dat" using 1:(log($6)) w l
//
//                 The "$" and both pairs of parens are needed for some
//                 reason.
//
// Special Notes : If tecplot level is set to 1, then output each dataset
//                 in a separate file.  If not, then append to a single file.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::outputTecplot ()
{
  int i;
  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);

  if (tecplotLevel == 1)
  {
    sprintf(filename,"%s_%03d.dat",outputName.c_str(),callsOTEC);
  }
  else
  {
    sprintf(filename,"%s.dat",outputName.c_str());
  }

  double time = getSolverState().currTime;
  FILE *fp1;

  if (tecplotLevel == 1)
  {
    fp1 = fopen(filename,"w");
  }
  else
  {
    if (callsOTEC <= 0)
    {
      fp1 = fopen(filename,"w");
    }
    else
    {
      fp1 = fopen(filename,"a");
    }
  }

  if (tecplotLevel == 1)
  {
    if (equationSet == 0)
    {
      fprintf(fp1,
              " TITLE = \"Spatially Dependent data for PDE diode: %s  time = %20.12e seconds. equation set = nonlinear Poisson\",\n",
              outputName.c_str(),time);
    }
    else
    {
      fprintf(fp1,
              " TITLE = \"Spatially Dependent data for PDE diode: %s  time = %20.12e seconds. equation set = drift diffusion\",\n",
              outputName.c_str(),time);
    }
  }
  else
  {
    if (callsOTEC <= 0)
    {
      fprintf(fp1,
              " TITLE = \"Spatially Dependent data for PDE diode: %s  time = %20.12e seconds.\",\n",
              outputName.c_str(),time);
    }
  }

  int rSize=0;
  int cSize=0;

  if (callsOTEC <= 0 || tecplotLevel == 1)
  {
    fprintf(fp1,"%s","\tVARIABLES = \"X \",\n");

    fprintf(fp1,"%s","\t    \"V \",\n");
    fprintf(fp1,"%s","\t    \"nn (electron dens.) \",\n");
    fprintf(fp1,"%s","\t    \"np (hole dens.) \",\n");
    fprintf(fp1,"%s","\t    \"nn*np (carrier product) \",\n");
    fprintf(fp1,"%s","\t    \"Dopant density \",\n");
    fprintf(fp1,"%s","\t    \"fabs(Dopant density)\",\n");
    fprintf(fp1,"%s","\t    \"electron lifetime \",\n");
    fprintf(fp1,"%s","\t    \"hole lifetime \",\n");
    //fprintf(fp1,"%s","\t    \"electron mobility \",\n");
    //fprintf(fp1,"%s","\t    \"hole mobility \",\n");
    fprintf(fp1,"%s","\t    \"Jn \",\n");
    fprintf(fp1,"%s","\t    \"Jp \",\n");
    fprintf(fp1,"%s","\t    \"R  \",\n");
    fprintf(fp1,"%s","\t    \"Ex \",\n");
    fprintf(fp1,"%s","\t    \"Idispl \", \n");

  }

  fprintf(fp1,"\tZONE F=POINT,I=%d", NX);

  if (getSolverState().dcopFlag)
  {
    fprintf(fp1,"  T = \"DCOP step = %d\" \n", callsOTEC);
  }
  else
  {
    fprintf(fp1,"  T = \"time step = %d time = %20.12e\" AUXDATA time = \"%20.12e seconds\" \n", callsOTEC , time, time);
  }

  double vcorrection = 0.0;
  if (useVoltageOutputOffset_)
  {
    if (offsetWithFirstElectrode_) // not the default.  This is here to match Wampler's 1D code.
    {
      vcorrection = -VVec[0]*scalingVars.V0;
    }
    else
    {
      vcorrection = VoltageOffset_;
    }
  }

  for (i=0;i<NX;++i)
  {
    fprintf(fp1,"  %20.12e",xVec[i]*scalingVars.x0);
    fprintf(fp1,"  %20.12e", (VVec[i]*scalingVars.V0 + vcorrection) );
    fprintf(fp1,"  %20.12e",nnVec[i]*scalingVars.C0);
    fprintf(fp1,"%s","\n");
    fprintf(fp1,"  %20.12e",npVec[i]*scalingVars.C0);
    fprintf(fp1,"  %20.12e",nnVec[i]*scalingVars.C0*npVec[i]*scalingVars.C0);
    fprintf(fp1,"  %20.12e",CVec[i]*scalingVars.C0);
    fprintf(fp1,"  %20.12e",fabs(CVec[i]*scalingVars.C0));
    fprintf(fp1,"%s","\n");
    fprintf(fp1,"  %20.12e",tnVec[i]*scalingVars.t0);
    fprintf(fp1,"  %20.12e",tpVec[i]*scalingVars.t0);
    //fprintf(fp1,"  %20.12e",unVec[i]*scalingVars.u0);
    fprintf(fp1,"%s","\n");
    //fprintf(fp1,"  %20.12e",upVec[i]*scalingVars.u0);
    fprintf(fp1,"  %20.12e",JnxVec[i]*scalingVars.J0);
    fprintf(fp1,"  %20.12e",JpxVec[i]*scalingVars.J0);
    fprintf(fp1,"%s","\n");
    fprintf(fp1,"  %20.12e",RVec[i]*scalingVars.R0);
    fprintf(fp1,"  %20.12e",ExVec[i]*scalingVars.E0);
    fprintf(fp1,"  %20.12e",displCurrent[i]*scalingVars.J0);
    fprintf(fp1,"%s","\n");

    fprintf(fp1,"%s","\n");
  }

  ++callsOTEC;
  fclose(fp1);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputSgplot
// Purpose       : Outputs data in a format readable by the simgen plotting
//                 tools.
//
// Special Notes : Despite the name, the output file should be plotted using
//                 the program oneplot(a 1D plotter), not sgplot
//                 (a 2d plotter).
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::outputSgplot ()
{
  int i;
  char fileName[32];

  static const int LEN_IDENT2 = 31;

  for (i = 0 ; i < 32; ++i)
    fileName[i] = static_cast<char>(0);

  sprintf(fileName,"%s_%03d.res",outputName.c_str(),callsOSG);
  ++callsOSG;

  FILE * handle1 = fopen(fileName, "w");

  UINT numArrays  = 3;
  double timeVar = 0.0;

  UINT inx = NX;

  char title[64];
  sprintf(title, "%s", "Xyce diodePDE 1D output");

  fwrite(&inx      , sizeof(  UINT), 1, handle1);  // array size.
  fwrite(&numArrays, sizeof(  UINT), 1, handle1);  // number of arrays, besides x.
  fwrite( title    , sizeof(  char),64, handle1);  // title
  fwrite(&timeVar  , sizeof(double), 1, handle1);  // time.

  char names[3][LEN_IDENT2];
  sprintf(names[0], "%s", "V");
  sprintf(names[1], "%s", "Ne");
  sprintf(names[2], "%s", "Np");

  // output the variable names, other than x:
  for(i=0;i<numArrays;++i)
  {
    fwrite(names[i], sizeof(char),(LEN_IDENT2), handle1);
  }

  double vcorrection = 0.0;
  if (useVoltageOutputOffset_)
  {
    if (offsetWithFirstElectrode_) // not the default.  This is here to match Wampler's 1D code.
    {
      vcorrection = -VVec[0]*scalingVars.V0;
    }
    else
    {
      vcorrection = VoltageOffset_;
    }
  }

  for (i=0;i<inx;++i)
  {
    xVec[i] *= scalingVars.x0;
    VVec[i] *= scalingVars.V0 + vcorrection;
    nnVec[i] *= scalingVars.C0;
    npVec[i] *= scalingVars.C0;
  }

  // output x-axis:
  fwrite( &xVec[0], sizeof(double),inx , handle1 );

  // output V
  fwrite( &VVec[0], sizeof(double),inx , handle1 );

  // output nn
  fwrite( &nnVec[0], sizeof(double),inx , handle1 );

  // output np
  fwrite( &npVec[0], sizeof(double),inx , handle1 );

  for (i=0;i<inx;++i)
  {
    xVec[i] /= scalingVars.x0;
    VVec[i] /= scalingVars.V0;
    nnVec[i] /= scalingVars.C0;
    npVec[i] /= scalingVars.C0;
  }

  fclose(handle1);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcRecombination
// Purpose       :
// Special Notes : This function assumes scaling is turned on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
bool Instance::calcRecombination ()
{
  if (!includeAugerRecomb && !includeSRHRecomb) return true;

  for (int i=0;i<NX;++i)
  {
    double Rsrh=0.0;
    double Raug=0.0;

    double n  = nnVec[i];
    double p  = npVec[i];
    double tn = tnVec[i];
    double tp = tpVec[i];

    if (includeSRHRecomb)
    {
      Rsrh = matSupport.calcRsrh (bulkMaterial, Ni,n,p,tn,tp);
    }

    if (includeAugerRecomb)
    {
      Raug = matSupport.calcRaug (bulkMaterial, Ni,n,p);
    }

    RVec[i] = (Rsrh + Raug);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdRecombination
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with the recombination term.
// Special Notes : This function assumes scaling is turned on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/01
//-----------------------------------------------------------------------------
bool Instance::pdRecombination ()
{
  if (!includeAugerRecomb && !includeSRHRecomb) return true;

  int i;
  double A, B, C;
  double dAdn, dAdp;
  double dBdn, dBdp;

  for (i=0;i<NX;++i)
  {
    double dRsrhdn=0.0;
    double dRsrhdp=0.0;
    double dRaugdn=0.0;
    double dRaugdp=0.0;

    // Rsrh derivatives: checklater.
    // (Rsrch = A*B)

    double n  = nnVec[i];
    double p  = npVec[i];
    double tn = tnVec[i];
    double tp = tpVec[i];

    if (includeSRHRecomb)
    {
      dRsrhdn = matSupport.pdRsrhN(bulkMaterial,Ni,n,p,tn,tp);
      dRsrhdp = matSupport.pdRsrhP(bulkMaterial,Ni,n,p,tn,tp);
    }

    if (includeAugerRecomb)
    {
      dRaugdn = matSupport.pdRaugN(bulkMaterial,Ni,n,p);
      dRaugdp = matSupport.pdRaugP(bulkMaterial,Ni,n,p);
    }

    dRdnVec[i] = dRsrhdn + dRaugdn;
    dRdpVec[i] = dRsrhdp + dRaugdp;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcElectronCurrent
// Purpose       :
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
bool Instance::calcElectronCurrent ()
{
  Ut = Vt/scalingVars.V0;
  for (int i=0;i<LX;++i)
  {
    JnxVec[i] =
      -J_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);
  }
  JnxVec[LX] = JnxVec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdElectronCurrent
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with electron current.
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/01
//-----------------------------------------------------------------------------
bool Instance::pdElectronCurrent ()
{
  Ut = Vt/scalingVars.V0;

  for (int i=0;i<LX;++i)
  {
    dJndn1Vec[i] =
      -dJdn1_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndn2Vec[i] =
      -dJdn2_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndV1Vec[i] =
      -dJdV1_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndV2Vec[i] =
      -dJdV2_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndp1Vec[i] =
      -dJdp1_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndp2Vec[i] =
      -dJdp2_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);
  }

  dJndn1Vec[LX] = dJndn1Vec[LX-1];
  dJndn2Vec[LX] = dJndn2Vec[LX-1];
  dJndV1Vec[LX] = dJndV1Vec[LX-1];
  dJndV2Vec[LX] = dJndV2Vec[LX-1];
  dJndp1Vec[LX] = dJndp1Vec[LX-1];
  dJndp2Vec[LX] = dJndp2Vec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcHoleCurrent
// Purpose       : This function assumes scaling is on.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
bool Instance::calcHoleCurrent ()
{
  Ut = Vt/scalingVars.V0;

  for (int i=0;i<LX;++i)
  {
    JpxVec[i] =
      J_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);
  }

  JpxVec[LX] = JpxVec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdHoleCurrent
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with the hole current.
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/01
//-----------------------------------------------------------------------------
bool Instance::pdHoleCurrent ()
{
  Ut = Vt/scalingVars.V0;

  for (int i=0;i<LX;++i)
  {
    dJpdp1Vec[i] =
      dJdn1_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdp2Vec[i] =
      dJdn2_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdV1Vec[i] =
      dJdV1_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdV2Vec[i] =
      dJdV2_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdn1Vec[i] =
      -dJdp1_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdn2Vec[i] =
      -dJdp2_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);
  }

  dJpdn1Vec[LX] = dJpdn1Vec[LX-1];
  dJpdn2Vec[LX] = dJpdn2Vec[LX-1];
  dJpdV1Vec[LX] = dJpdV1Vec[LX-1];
  dJpdV2Vec[LX] = dJpdV2Vec[LX-1];
  dJpdn1Vec[LX] = dJpdn1Vec[LX-1];
  dJpdn2Vec[LX] = dJpdn2Vec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcEfield
// Purpose       : This function works with or without scaling.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
bool Instance::calcEfield ()
{
  double absEx;
  Emax = 0.0;

  for (int i=0;i<LX;++i)
  {
    ExVec[i] = -(VVec[i+1] - VVec[i])/dxVec[i];

    absEx = fabs(ExVec[i]);
    if (absEx > Emax) Emax = absEx;
  }
  Emax *= scalingVars.E0;

  ExVec[LX] = ExVec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::enablePDEContinuation
// Purpose       : Sets up the various parameters neccessary for a continuation
//                 calculation.  Mainly, it sets up the voltage step size
//                 for all the voltage BC's.
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool Instance::enablePDEContinuation ()
{
  bool bnoChange = true;
  int iBC;
  int bcSize=bcVec.size();

  continuationAlpha = 1.0;

  if (!enableContinuationCalled)
  {
    for (iBC=0;iBC<bcSize;++iBC)
    {
      bcVec[iBC].Vckt_old = bcVec[iBC].Vckt;
    }
  }

  obtainNodeVoltages ();

  for (iBC=0;iBC<bcSize;++iBC)
  {
    bcVec[iBC].Vckt_final = bcVec[iBC].Vckt;
    bcVec[iBC].Vckt_orig  = bcVec[iBC].Vckt;
  }

  // This (voltlim) is a very new thing.  Use carefully...
  if (getDeviceOptions().voltageLimiterFlag && voltLimFlag)
  {
    applyVoltageLimiting ();
  }

  for (iBC=0;iBC<bcSize;++iBC)
  {
    double dV,tmp1V, tmp2V;
    tmp1V = bcVec[iBC].Vckt_final;
    tmp2V = bcVec[iBC].Vckt_old;
    dV    = tmp1V - tmp2V;

    bcVec[iBC].Vckt_delta = dV;

    bcVec[iBC].Vckt_deltaC = dV/
                             (static_cast<double>(getSolverState().maxPDEContinuationSteps));

    // if this deltaC is too big, then we need to change the
    // number of continuation steps.
    double maxDelta = maxVoltDelta;

    if (fabs(bcVec[iBC].Vckt_deltaC) > maxDelta)
    {
      int tmp_steps = static_cast<int>(fabs(dV)/maxDelta) + 1;
      getSolverState().maxPDEContinuationSteps = tmp_steps;

      bcVec[iBC].Vckt_deltaC = dV/
                               static_cast<double>(getSolverState().maxPDEContinuationSteps);
    }

    if (fabs(dV) > 1.0e-3) bnoChange = false;

    bcVec[iBC].Vckt_ramp     = bcVec[iBC].Vckt_old;
    bcVec[iBC].Vckt_ramp_old = bcVec[iBC].Vckt_old;
  }

#if 0
  // now set up any continuation params associated with photogen.
  bool bnoChangePhotogen;

  bnoChangePhotogen = enablePhotogenContinuation ();

  bnoChange = bnoChange && bnoChangePhotogen;
#endif

  if (!enableContinuationCalled) enableContinuationCalled = true;

  // if none of the boundary conditions  have changed, then
  // return a false.
  return (!bnoChange);
}

//-----------------------------------------------------------------------------
// Function      : Instance::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool Instance::disablePDEContinuation ()
{
  int iBC;

  int bcSize=bcVec.size();
  for (iBC=0;iBC<bcSize;++iBC)
  {
    bcVec[iBC].Vckt_old   = bcVec[iBC].Vckt_final;
  }

#if 0
  photoA1_old = photoA1_final;
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setPDEContinuationAlpha
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/18/03
//-----------------------------------------------------------------------------
void Instance::setPDEContinuationAlpha (double alpha)
{
#ifdef Xyce_DEBUG_DEVICE
  Xyce::dout() << section_divider << std::endl;
  Xyce::dout() << "Instance::setPDEContinuationAlpha" << std::endl;
#endif


  // First do the voltage boundary conditions:
  int bcSize=bcVec.size();
  for (int iBC=0;iBC<bcSize;++iBC)
  {
    bcVec[iBC].Vckt_ramp = bcVec[iBC].Vckt_old + (bcVec[iBC].Vckt_delta)*alpha;

    // make sure we haven't gone too far:
    if ((bcVec[iBC].Vckt_delta >  0 && bcVec[iBC].Vckt_ramp >  bcVec[iBC].Vckt_final) ||
        (bcVec[iBC].Vckt_delta <= 0 && bcVec[iBC].Vckt_ramp <= bcVec[iBC].Vckt_final) )
    {
      bcVec[iBC].Vckt_ramp = bcVec[iBC].Vckt_final;
    }

#ifdef Xyce_DEBUG_DEVICE
    Xyce::dout() << "  " << bcVec[iBC].eName << "  Vckt_ramp = " << bcVec[iBC].Vckt_ramp << std::endl;
#endif
  }

#if 0
  // now do the photogeneration term, if neccessary:
  photoA1_ramp =  photoA1_old + photoA1_Delta * alpha;

  // make sure we haven't gone too far:
  if ((photoA1_Delta >  0 && photoA1_ramp >  photoA1_final) ||
      (photoA1_Delta <= 0 && photoA1_ramp <= photoA1_final) )
  {
    photoA1_ramp = photoA1_final;
  }

#ifdef Xyce_DEBUG_DEVICE
  Xyce::dout() << "  photoA1_ramp = " << photoA1_ramp << std::endl;
  Xyce::dout() << section_divider << std::endl;
#endif
#endif
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  Device *device = new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);

  return device;
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("pde", 1)
    .registerModelType("pde", 1)
    .registerModelType("zod", 1);
}

} // namespace DiodePDE
} // namespace Device
} // namespace Xyce
