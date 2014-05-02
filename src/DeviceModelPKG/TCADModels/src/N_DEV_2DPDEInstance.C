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
// Filename       : $RCSfile: N_DEV_2DPDEInstance.C,v $
//
// Purpose        : This file contains a lot of the
//                  implementation of the instance class for the two
//                  dimensional PDE based semiconductor device.
//
//                  Functions pertaining to the initial setup are in other
//                  files, as are functions relating to mesh handling and
//                  parameter handling.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/05/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.152 $
//
// Revision Date  : $Date: 2014/02/24 23:49:18 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ----------  Standard Includes ----------
#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_2DPDE.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_PDE_2DMesh.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_DEV_SourceData.h>

#include <N_DEV_DiodePDE.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {
namespace TwoDPDE {

// note that this macros is only a default - it will probably not be used.
// default maximum number of nonzero entries in a matrix row
static const int MAX_COLS_PER_ROW = 10;


void Traits::loadInstanceParameters(ParametricData<TwoDPDE::Instance> &p)
{
  p.addPar("AREA", 1.0, &TwoDPDE::Instance::area)
    .setExpressionAccess(TIME_DEP);
  p.addPar("NA", 1.0e+15, &TwoDPDE::Instance::Na)
    .setExpressionAccess(TIME_DEP);
  p.addPar("ND", 1.0e+15, &TwoDPDE::Instance::Nd)
    .setExpressionAccess(TIME_DEP);
  p.addPar("WJ", 1.0e-4, &TwoDPDE::Instance::WJ)
    .setExpressionAccess(TIME_DEP);
  p.addPar("TEMP", 0.0, &TwoDPDE::Instance::Temp)
    .setExpressionAccess(TIME_DEP);
  p.addPar("X0", 1.0e-4, &TwoDPDE::Instance::x0_user)
    .setExpressionAccess(TIME_DEP);
  p.addPar("XSTART", 0.0, &TwoDPDE::Instance::xstart)
    .setExpressionAccess(TIME_DEP);
  p.addPar("YSTART", 0.0, &TwoDPDE::Instance::ystart)
    .setExpressionAccess(TIME_DEP);
  p.addPar("XEND", 0.0, &TwoDPDE::Instance::xend)
    .setExpressionAccess(TIME_DEP);
  p.addPar("YEND", 0.0, &TwoDPDE::Instance::yend)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PH.A1", 0.0, &TwoDPDE::Instance::photoA1)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PH.TSTART", 0.0, &TwoDPDE::Instance::photoTstart)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PH.TSTOP", 0.0, &TwoDPDE::Instance::photoTstop)
    .setExpressionAccess(TIME_DEP);
  p.addPar("INT", 0.0, &TwoDPDE::Instance::intensity)
    .setExpressionAccess(TIME_DEP);
  p.addPar("L", 1.0e-3, &TwoDPDE::Instance::deviceLength)
    .setExpressionAccess(TIME_DEP);
  p.addPar("W", 1.0e-3, &TwoDPDE::Instance::deviceWidth)
    .setExpressionAccess(TIME_DEP);
  p.addPar("MAXVOLTDELTA", 0.025, &TwoDPDE::Instance::maxVoltDelta)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PH.TD", 0.0, &TwoDPDE::Instance::photoTd)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PH.TR", 0.0, &TwoDPDE::Instance::photoTr)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PH.TF", 0.0, &TwoDPDE::Instance::photoTf)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PH.PW", 0.0, &TwoDPDE::Instance::photoPw)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PH.PER", 0.0, &TwoDPDE::Instance::photoPer)
    .setExpressionAccess(TIME_DEP);
  p.addPar("OUTPUTINTERVAL", 0.0, &TwoDPDE::Instance::outputInterval)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PENALTYPREFAC", 1.0, &TwoDPDE::Instance::penaltyPrefac)
    .setGivenMember(&TwoDPDE::Instance::penaltyPrefacGiven)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PENALTYPOW", 2.0, &TwoDPDE::Instance::penaltyPow)
    .setGivenMember(&TwoDPDE::Instance::penaltyPowGiven)
    .setExpressionAccess(TIME_DEP);
  p.addPar("PULSEDATA", 0.0, &TwoDPDE::Instance::PulseData)
    .setExpressionAccess(TIME_DEP);

  // Set up map for non-double precision variables:
  p.addPar("MESHFILE", std::string("internal.msh"), &TwoDPDE::Instance::meshFileName);
  p.addPar("GRADED", false, &TwoDPDE::Instance::gradedJunctionFlag);
  p.addPar("MOBMODEL", std::string("ARORA"), &TwoDPDE::Instance::mobModelName);
  p.addPar("BULKMATERIAL", std::string("SI"), &TwoDPDE::Instance::bulkMaterial);
#ifdef Xyce_OXIDE_ENABLED
  p.addPar("ALLOXIDE", false, false, Pars::NO_DEP, &TwoDPDE::Instance::allOxideFlag);
#endif
  p.addPar("DISPLCUR", false, &TwoDPDE::Instance::displCurrentFlag);
  p.addPar("PHOTOGEN", false, &TwoDPDE::Instance::photogenOnFlag);
  p.addPar("PH.TYPE", std::string("UNIFORM"), &TwoDPDE::Instance::photoString);
  p.addPar("TECPLOTLEVEL", 1, &TwoDPDE::Instance::tecplotLevel);
  p.addPar("INTERPGRIDSIZE", 20, &TwoDPDE::Instance::interpGridSize);
  p.addPar("SGPLOTLEVEL", 0, &TwoDPDE::Instance::sgplotLevel);
  p.addPar("GNUPLOTLEVEL", 0, &TwoDPDE::Instance::gnuplotLevel);
  p.addPar("TXTDATALEVEL", 1, &TwoDPDE::Instance::txtDataLevel);
  p.addPar("VOLTLIM", false, &TwoDPDE::Instance::voltLimFlag);
  p.addPar("USEMATRIXGID", false, &TwoDPDE::Instance::useMatrixGIDFlag);
  p.addPar("USEVECTORGID", false, &TwoDPDE::Instance::useVectorGIDFlag);
  p.addPar("CONSTBOUNDARY", false, &TwoDPDE::Instance::constBoundaryFlag);
  p.addPar("NX", 11, &TwoDPDE::Instance::numMeshPointsX);
  p.addPar("NY", 11, &TwoDPDE::Instance::numMeshPointsY);
  p.addPar("CYL", false, &TwoDPDE::Instance::cylGeomFlag);
  p.addPar("OUTPUTNLPOISSON", false, &TwoDPDE::Instance::outputNLPoisson);
  p.addPar("PENALTY", false, &TwoDPDE::Instance::penaltyFlag);
  p.addPar("TYPE", std::string("PNP"), &TwoDPDE::Instance::deviceType);
  p.addPar("USEOLDNI", false, &TwoDPDE::Instance::useOldNi)
    .setGivenMember(&TwoDPDE::Instance::useOldNiGiven)
    .setUnit(U_LOGIC)
    .setDescription("Flag for using old (inaccurate) intrinsic carrier calculation.");

  p.addComposite ("NODE", PDE_2DElectrode::getParametricData(), &TwoDPDE::Instance::electrodeMap);
  p.addComposite ("DOPINGPROFILES", DopeInfo::getParametricData(), &TwoDPDE::Instance::dopeInfoMap);
  p.addComposite ("REGION", DopeInfo::getParametricData(), &TwoDPDE::Instance::dopeInfoMap);
}

// Class Instance

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &               IB,
  Model &                       Miter,
  const FactoryBlock &factory_block)
  : DevicePDEInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    Is(1.0e-14),
    Id(0.0),
    Emax(0.0),
    VminExp(0.0),
    VmaxExp(0.0),
    dIVec(),
    penaltyPrefacGiven(false),
    penaltyPowGiven(false),
    LeadCurrent1(0.0),
    LeadCurrent2(0.0),
    LeadCurrent3(0.0),
    LeadCurrent4(0.0),
    LeadCurrent5(0.0),
    LeadCurrent6(0.0),
    LeadCurrent7(0.0),
    LeadCurrent8(0.0),
    Na(1.0e15),
    Nd(1.0e15),
    Vbi(0.0),
    //Vbi(Vt*log(Na*Nd/(Ni*Ni))),
    WJ(1.0e-4),
    XC(0.0),
    XL(0.0),
    XR(0.0),
    NnMax(1.0e15),
    NpMax(1.0e15),
    NnMin(1.0e5),  // approx...
    NpMin(1.0e5),
    useOldNi(false),
    useOldNiGiven(false),
    meshFileName(""),
    deviceType("PNP"),
    usingInternalMesh(false),

    deviceInitialized(false),
    meshPerturbed (false),
    dopingPerturbed (false),
    photogenPerturbed (false),

    numMeshPointsX(11),
    numMeshPointsY(11),
    deviceLength(1.0e-3),
    deviceWidth(1.0e-3),
    cylGeomFlag(false),
    penaltyFlag(false),
    penaltyPrefac(1.0),
    penaltyPow(2.0),
    area(1.0),
#ifdef Xyce_OXIDE_ENABLED
    allOxideFlag(false),
#endif

    gradedJunctionFlag(false),
    calledBeforeSIGB(false),
    callsOSG (0),
    callsOTEC(0),
    callsOTECvec(0),
    callsOGNU(0),
    callsOTXT(0),

    displCurrentFlag(false),
    constBoundaryFlag(false),
    calcConductanceFlag(false),
    equationSet(0),
    outputInterval(0.0),
    outputIndex(0),
    outputNLPoisson(false),
    lastOutputTime(-10.0),
    tecplotLevel(0),
    sgplotLevel(0),
    gnuplotLevel(0),
    txtDataLevel(1),
    interpGridSize(20),
    voltLimFlag(false),
    useMatrixGIDFlag(false),
    useVectorGIDFlag(false),
    meshContainerPtr(0),
    meshCopyContainerPtr(0),

    xVec(),
    yVec(),
    CVec(),
    minDXVec(),
    areaVec(),
    VVec(),
    nnVec(),
    npVec(),
    totSrcVec(),
    RVec(),
    SVec(),
    elecPenalty(),
    holePenalty(),
    pdElecPenalty(),
    pdHolePenalty(),
    unVec(),
    upVec(),
    unE_Vec(),
    upE_Vec(),
    tnVec(),
    tpVec(),
    EfieldVec(),
    JnVec(),
    JpVec(),
    displPotential(),
    displCurrent(),
    outputVec(),

    dRdpVec(),
    dRdnVec(),
    dJndn1Vec(),
    dJndn2Vec(),
    dJndV1Vec(),
    dJndV2Vec(),
    dJpdn1Vec(),
    dJpdn2Vec(),
    dJpdV1Vec(),
    dJpdV2Vec(),
    boundarySten(),
    boundaryStenV(),
    boundaryStenN(),
    boundaryStenP(),
    boundaryTest(),
    boundaryNeighborSten(),
    Vrowarray(),
    Vcolarray(),
    Nrowarray(),
    Ncolarray(),
    Prowarray(),
    Pcolarray(),
    vOwnVec(), // on processor tag, for voltage variables.
    nnOwnVec(), // on processor tag, for elec. density variables.
    npOwnVec(), // on processor tag, for hole density variables.
    li_Vrowarray(),
    li_Nrowarray(),
    li_Prowarray(),
    li_VoffsetArray(),
    li_NoffsetArray(),
    li_PoffsetArray(),
    MESHtoLID_V(),
    MESHtoLID_N(),
    MESHtoLID_P(),
    aiEdge(),
    aiEdge_nf(),
    iNumPlotEdges(0),
    iNumPlotEdges_nf(0),
    tmpBCmap(),
    labelIndex(),
    labelNameVector(),
    labelDIMap(),
    meshNeighborMultiMap(),
    electrodeMap(),
    stateDispl(),
    stateDispl_owned(),
    li_stateDispl(),

    numMeshPoints(0),
    numInterfaceMeshPoints(0),
    numMeshEdges (0),
    numMeshCells (0),
    numMeshLabels (0),
    maxColsPerRow(MAX_COLS_PER_ROW),
    numElectrodes (0),
    condVec(),
    capVec(),
    pdTermsAllocated(false),
    meshToLID(),
    jacStamp()
{

#if 0
  vector<Param>::const_iterator iter_t;
  vector<Param>::const_iterator begin_t = parameterVec.begin();
  vector<Param>::const_iterator end_t   = parameterVec.end();

  for (iter_t = begin_t; iter_t !=  end_t;  ++iter_t)
  {
    Xyce::dout() << "Tag: " << iter_t->tag();
    Xyce::dout() << "  Value = " << iter_t->sVal() << "  Given: " << iter_t->given() << endl;
  }
#endif


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    Temp = getDeviceOptions().temp.getImmutableValue<double>();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  ExtendedString tmpName = meshFileName;
  tmpName.toLower();
  meshFileName = tmpName;
  tmpName = mobModelName;
  tmpName.toLower();
  mobModelName = tmpName;
  tmpName = bulkMaterial;
  tmpName.toLower();
  bulkMaterial = tmpName;
  if (given("PH.TYPE"))
  {
    if (photoString == "UNIFORM")
      photoType = _DC_DATA;
    else if (photoString == "PULSE")
      photoType = _PULSE_DATA;
    else
    {
      std::string msg = "::: ph.type not recognized.\n";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
    }
  }

  if ( given("TYPE") )
  {
    ExtendedString tmpTypeName = deviceType;
    tmpTypeName.toUpper ();
    deviceType = tmpTypeName;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
    Xyce::dout() << "Doing standard initialization of ."<< std::endl;
#endif
  bool bsuccess = setupMesh ();
  bool bs1 = true;

  bs1 = doAllocations ();           bsuccess = bsuccess && bs1;
  bs1 = setupDINodes ();            bsuccess = bsuccess && bs1;
  bs1 = setupBCEdgeAreas ();        bsuccess = bsuccess && bs1;
  bs1 = checkForElectrodeOverlap ();bsuccess = bsuccess && bs1;
  bs1 = setupBoundaryStencil ();    bsuccess = bsuccess && bs1;
  bs1 = setupNumVars ();            bsuccess = bsuccess && bs1;
  bs1 = setupLabelIndex ();         bsuccess = bsuccess && bs1;
  bs1 = setupMinDXVector ();        bsuccess = bsuccess && bs1;
  bs1 = setupJacStamp ();           bsuccess = bsuccess && bs1;
  bs1 = setupMiscConstants ();      bsuccess = bsuccess && bs1;
  bs1 = setupScalingVars ();        bsuccess = bsuccess && bs1;
  bs1 = calcDopingProfile ();       bsuccess = bsuccess && bs1;
  bs1 = setupPhotogen ();           bsuccess = bsuccess && bs1;

  deviceInitialized = true;
  processParams ();

  devConMap.resize(numExtVars);
  for (int i=0 ; i<numExtVars ; ++i)
    devConMap[i] = i;
}


//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
Instance::~Instance()
{

  if (meshCopyContainerPtr != 0) delete meshCopyContainerPtr;
  if (meshContainerPtr != 0)     delete meshContainerPtr;

  if (DataSaved_ptr != 0) delete DataSaved_ptr;

  // Loop over the dopeInfoMap (if it is not empty) and delete its contents.
  if (!(dopeInfoMap.empty()))
  {
    std::map<std::string,DopeInfo *>::iterator iter;
    std::map<std::string,DopeInfo *>::iterator begin = dopeInfoMap.begin();
    std::map<std::string,DopeInfo *>::iterator end   = dopeInfoMap.end  ();

    for(iter=begin;iter!=end;++iter)
    {
      if (iter->second != 0) delete iter->second;
    }
  }

  // Loop over the electrodeMap (if it is not empty) and delete its contents.
  if (!(electrodeMap.empty()))
  {
    std::map<std::string, PDE_2DElectrode * >::iterator iterE;
    std::map<std::string, PDE_2DElectrode * >::iterator beginE = electrodeMap.begin();
    std::map<std::string, PDE_2DElectrode * >::iterator endE   = electrodeMap.end  ();

    for(iterE=beginE;iterE!=endE;++iterE)
    {
      if (iterE->second != 0) delete iterE->second;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::constructComposite
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/09/05
//-----------------------------------------------------------------------------
CompositeParam *Instance::constructComposite(const std::string & compositeName, const std::string & paramName)
{
  if (compositeName == "DOPINGPROFILES" || compositeName == "REGION")
  {
    DopeInfo *n = new DopeInfo();
    dopeInfoMap[paramName] = n;
    return (static_cast<CompositeParam *> (n));
  }
  if (compositeName == "NODE")
  {
    DeviceInterfaceNode dINode;
    ExtendedString dIName = paramName;
    dIName.toUpper ();

    dINode.eName = dIName;
    dINode.nName = paramName;
    dINode.given = true;
    dINode.index = 0;

    if (dINode.given) ++numElectrodes;
    if (dINode.given) dIVec.push_back(dINode);

    PDE_2DElectrode *n = new PDE_2DElectrode();
    electrodeMap[paramName] = n;
    return (static_cast<CompositeParam *> (n));
  }
  std::string msg =
    "Instance::constructComposite: unrecognized composite name: ";
  msg += compositeName;
  N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);

  return NULL;
}

// Additional Declarations
//
//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermdiateVars
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;
  bool bs1 = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "updateIntermediateVars.  name = " << getName() << std::endl;
  }
#endif

  // This first call, to setInitialGuess, probably won't do
  // anything, as this function only does
  // anything the first time that it is called, and it is usually called
  // from the time integrator at the very beginning of the run.
  bs1 = setInitialGuess (); bsuccess = bsuccess && bs1;
  bs1 = obtainSolution ();       bsuccess = bsuccess && bs1;
  bs1 = calcEfield ();           bsuccess = bsuccess && bs1;
  bs1 = calcMobilities ();       bsuccess = bsuccess && bs1;
  bs1 = calcRecombination ();    bsuccess = bsuccess && bs1;
  bs1 = calcElectronCurrent ();  bsuccess = bsuccess && bs1;
  bs1 = calcHoleCurrent ();      bsuccess = bsuccess && bs1;
  bs1 = calcTerminalCurrents (); bsuccess = bsuccess && bs1;
  bs1 = calcTerminalCharges ();  bsuccess = bsuccess && bs1;

  if (photogenOnFlag)
  {
    bs1 = calcPhotogen ();
    bsuccess = bsuccess && bs1;
  }

  if (penaltyFlag)
  {
    bs1 = calcPenalty ();
    bsuccess = bsuccess && bs1;
  }

  sumSources ();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bs1;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcTerminalCurrents
// Purpose       : Calculates total device current(s) to be used in the
//                 circuit KCL equations.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcTerminalCurrents ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "calcTerminalCurrents.  name = " << getName() << std::endl;
  }
#endif

  // loop over the device interface nodes, sum the currents going into each one.
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end   ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  int ind=1;
  for (iterDI=firstDI; iterDI!=lastDI; ++iterDI, ++ind)
  {
    // loop over the nodes of this device interface node:

    if ( !( meshContainerPtr->labelEdgeType (iterDI->eName) ) ) continue;

    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

    std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
    std::vector<int>::iterator iterI;

    iterDI->currentSum = 0.0;

    int nodeIndex;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {

      // loop over neighbor nodes/edges to get current sum for this node.
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);

      std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
      std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end  ();
      std::vector<EDGEINFO>::iterator iterEI;

      double sum       = 0.0; // total current for this node

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << " --------------- " << std::endl;
        Xyce::dout() << "name = " << iterDI->eName;
        Xyce::dout() << "  node      = " << *iterI << std::endl;
      }
#endif

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        int iedge = iterEI->iedge;
        int neighbor = iterEI->inode;
        mEdge * edgePtr = meshContainerPtr->getEdge(iedge);
        double ilen = edgePtr->ilen;

        double sign = (*iterI < neighbor) ? 1.0 : -1.0;

        sum += (sign*JnVec[iedge] + sign*JpVec[iedge])*ilen;
      }

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << " sum*scalingVars.a0 = "<< sum*scalingVars.a0 << std::endl;
        Xyce::dout() << " sum    = " << sum << std::endl;
        Xyce::dout() << " --------------- " << std::endl;
      }
#endif
      // total scaled current mult. by total scaled area for this node.
      sum *= scalingVars.a0;

      iterDI->currentSum += sum;
    } // iterI loop.

    iterDI->currentSum *= scalingVars.J0;
    if (ind == 1)
      LeadCurrent1 = iterDI->currentSum;
    else if (ind == 2)
      LeadCurrent2 = iterDI->currentSum;
    else if (ind == 3)
      LeadCurrent3 = iterDI->currentSum;
    else if (ind == 4)
      LeadCurrent4 = iterDI->currentSum;
    else if (ind == 5)
      LeadCurrent5 = iterDI->currentSum;
    else if (ind == 6)
      LeadCurrent6 = iterDI->currentSum;
    else if (ind == 7)
      LeadCurrent7 = iterDI->currentSum;
    else if (ind == 8)
      LeadCurrent8 = iterDI->currentSum;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << "  " << iterDI->eName;
      Xyce::dout() << " Ickt = " << iterDI->currentSum << std::endl;
    }
#endif
  }  // iterDI loop.


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
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
//                 needed for each DeviceInterfaceNode:
//
//                  dIdVckt - derivative of terminal current w.r.t. Vckt.
//                            This is the also the Jacobian contribution
//                            for the (KCL row, KCL col) entry of the matrix.
//
//                  dFdVckt - (col. vector) derivative of the RHS vector
//                            w.r.t. Vckt.  This is a vector quantity, and
//                            corresponds to the (*, Vckt) column of the
//                            PDE matrix sub-block.
//
//                  dIdX  - (row vector) derivative of the terminal current
//                          w.r.t. the vector of PDE solution variables. (ie not
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
// Creation Date : 12/01/02
//-----------------------------------------------------------------------------
bool Instance::pdTerminalCurrents ()
{
  bool bsuccess = true;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  std::string msg;
  std::string semi(bulkMaterial);

  int nodeIndex;
  int iedge  = 0;
  int inodeB = 0;
  double sign = 1.0;
  double dJndV = 0.0;
  double dJpdV = 0.0;
  double dJndn = 0.0;
  double dJpdp = 0.0;
  double coef;
  double nodeArea,ilen,elen;

  if (!pdTermsAllocated)
  {
    allocatePDTerms ();
    pdTermsAllocated = true;
  }

  // first calculate dIdVckt.----------------------------------------
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    iterDI->dIdVckt = 0.0;

    // if using the old BC, dIdVckt is just zero.
#ifdef Xyce_NEW_BC

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << Xyce::subsection_divider << std::endl;
    }
#endif

    if (iterDI->gid ==-1) continue;

    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

    // obtain the node indices for the current label, loop over them.
    std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
    std::vector<int>::iterator iterI;

    std::vector<EDGEINFO>::iterator firstEI;
    std::vector<EDGEINFO>::iterator lastEI;
    std::vector<EDGEINFO>::iterator iterEI;

    // for the "new" boundary conditions, the currents coming into
    // the electrode have a direct dependency on the voltage at
    // the circuit node, so they all contribute to the dependence
    // of that KCL equation on that node's voltage.  The contributions
    // here are all the equivalent of those that disappear due to
    // the removal of the boundary mesh nodes from the system
    // of equations.
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // do the center point first.
      coef = 0.0;
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
      {
        iedge = iterEI->iedge;
        inodeB = iterEI->inode;
        mEdge * edgePtr = meshContainerPtr->getEdge(iedge);
        ilen = edgePtr->ilen;

        // if this is a boundary edge, skip it, because the neighbor point
        // is not part of the solution vector.  Only use
        // those which extend into the interior.
        if (boundarySten[inodeB]) continue;

        sign = (*iterI < inodeB) ? 1.0 : -1.0;

        if (*iterI<inodeB)
        {
          dJndV = dJndV1Vec[iedge];
          dJpdV = dJpdV1Vec[iedge];
        }
        else
        {
          dJndV = dJndV2Vec[iedge];
          dJpdV = dJpdV2Vec[iedge];
        }

        coef += (sign* dJndV + sign* dJpdV)*ilen;
      } // end of nn edge loop

      double tmpsum = (scalingVars.rV0)*coef*scalingVars.J0*scalingVars.a0;

      iterDI->dIdVckt += tmpsum;

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        Xyce::dout().setf(std::ios::left);
        Xyce::dout() << iterDI->eName<< ":";
        Xyce::dout() << " KCL pdTerminalCurrent row = " << iterDI->gid;
        Xyce::dout() << " contrib = " << tmpsum;
        Xyce::dout() << " dIdVckt = " << iterDI->dIdVckt << std::endl;
      }
#endif

    } // end of node loop

#endif // Xyce_NEW_BC

  }  // end of DI loop
#ifdef Xyce_NEW_BC
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::subsection_divider << std::endl;
  }
#endif
#endif // Xyce_NEW_BC

  // Now calculate dFdVckt.----------------------------------------
#ifdef Xyce_NEW_BC
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {

    // zero out the dFdVckt vector.
    int itmp=0;
    int size = iterDI->dFdVckt.size();
    for (itmp=0;itmp<size;++itmp)
    {
      iterDI->dFdVckt[itmp] = 0.0;
    }

    int numNeighbor = iterDI->neighborNodes.size();

    int iNeighbor;
    int dFdVindex = 0;
    for (iNeighbor=0;iNeighbor<numNeighbor;++iNeighbor)
    {
      int inode = iterDI->neighborNodes[iNeighbor];

      mNode * nodePtr = meshContainerPtr->getNode(inode);
      nodeArea = nodePtr->area;

      double coef = 0.0;
      int iNN=0;
      for (iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        int inodeB = nodePtr->edgeInfoVector[iNN].inode;

        // if nodeB is not a boundary node, never mind.
        if (boundarySten[inodeB]!=1) continue;

        // if it is a boundary node, but not part of the
        // current boundary, also never mind.
        if (labelNameVector[inodeB]!= iterDI->eName) continue;

        ilen   = nodePtr->edgeInfoVector[iNN].ilen;
        elen   = nodePtr->edgeInfoVector[iNN].elen;
        iedge  = nodePtr->edgeInfoVector[iNN].iedge;

        // poisson equation contribution:
        coef   =  ilen/elen;
        coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;
        coef  *= scalingVars.rV0;

        iterDI->dFdVckt[dFdVindex] = coef;
        ++dFdVindex;

        // electron equation contribution:
        double dJdV = 0.0;
        if (inode>inodeB) { dJdV = dJndV1Vec[iedge]; }
        else              { dJdV = dJndV2Vec[iedge]; }

        coef = ((inode<inodeB)?1.0:-1.0) * dJdV * ilen/nodeArea;
        coef *= scalingVars.rV0;

        iterDI->dFdVckt[dFdVindex] = coef;
        ++dFdVindex;

        // hole equation contribution:
        dJdV = 0.0;
        if (inode>inodeB) { dJdV = dJpdV1Vec[iedge]; }
        else              { dJdV = dJpdV2Vec[iedge]; }

        coef = -((inode<inodeB)?1.0:-1.0) * dJdV * ilen/nodeArea;
        coef *= scalingVars.rV0;

        iterDI->dFdVckt[dFdVindex] = coef;
        ++dFdVindex;
      }
    }
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "-----" << std::endl;
      Xyce::dout() << "neighbor nodes for boundary: " << iterDI->eName << std::endl;
      Xyce::dout() << "-----" << std::endl;
      for (iNeighbor=0;iNeighbor<numNeighbor;++iNeighbor)
      {
        int inode = iterDI->neighborNodes[iNeighbor];
        Xyce::dout() << "\t"<<iNeighbor<<"  "<< inode << std::endl;
      }

      Xyce::dout() << "-----" << std::endl;
      Xyce::dout() << "dFdVckt vector for boundary: " << iterDI->eName << std::endl;
      Xyce::dout() << "-----" << std::endl;

      iNeighbor=0;
      int idf = 0;
      for (iNeighbor=0;iNeighbor<numNeighbor;++iNeighbor)
      {
        int inode = iterDI->neighborNodes[iNeighbor];
        int Vrow = Vrowarray[inode];
        int Nrow = Nrowarray[inode];
        int Prow = Prowarray[inode];
        mNode * nodePtr = meshContainerPtr->getNode(inode);

        double dfdV = 0.0;

        int iNN=0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          int inodeB = nodePtr->edgeInfoVector[iNN].inode;

          // if nodeB is not a boundary node, never mind.
          if (boundarySten[inodeB]!=1) continue;

          // if it is a boundary node, but not part of the
          // current boundary, also never mind.
          if (labelNameVector[inodeB]!= iterDI->eName) continue;

          dfdV = iterDI->dFdVckt[idf];
          Xyce::dout() << "\t"<<idf;
          Xyce::dout() << "  \tv_"<<inode<<"\t"<<Vrow<<"\t"<<dfdV<< std::endl;
          ++idf;

          dfdV = iterDI->dFdVckt[idf];
          Xyce::dout() << "\t"<<idf;
          Xyce::dout() << "  \tn_"<<inode<<"\t"<<Nrow<<"\t"<<dfdV<< std::endl;
          ++idf;

          dfdV = iterDI->dFdVckt[idf];
          Xyce::dout() << "\t"<<idf;
          Xyce::dout() << "  \tp_"<<inode<<"\t"<<Prow<<"\t"<<dfdV<< std::endl;
          ++idf;
        }// end of iNN loop
      } // end of iNeighbor loop
      Xyce::dout() << "-----" << std::endl;
    }
#endif
  }
#else // old BC not set up yet.

#endif // Xyce_NEW_BC

  // now do dIdX.-------------------------------------------------
  double Vcoef;
  double Ncoef;
  double Pcoef;
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

    // zero out the dIdX vector.
    int numdIdX = iterDI->dIdX.size ();
    for (int j=0;j<numdIdX;++j)
    {
      iterDI->dIdX[j] = 0.0;
    }

    // obtain the node indices for the current label, loop over them.
    // for each edge node, add an extra column entry to the colarray.

    std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
    std::vector<int>::iterator iterI;

    std::vector<EDGEINFO>::iterator firstEI;
    std::vector<EDGEINFO>::iterator lastEI;
    std::vector<EDGEINFO>::iterator iterEI;

    int iVcol = 0;
    int iNcol = 0;
    int iPcol = 0;
    int col1;
    int cnt2;
    bool bmatch;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // do the center point first.
      Vcoef = 0.0;
      Ncoef = 0.0;
      Pcoef = 0.0;
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
      {
        iedge = iterEI->iedge;
        inodeB = iterEI->inode;
        sign = (*iterI < inodeB) ? 1.0 : -1.0;
        mEdge * edgePtr = meshContainerPtr->getEdge(iedge);
        ilen = edgePtr->ilen;

        if (*iterI<inodeB)
        {
          dJndV = dJndV1Vec[iedge];
          dJpdV = dJpdV1Vec[iedge];
          dJndn = dJndn1Vec[iedge];
          dJpdp = dJpdn1Vec[iedge];
        }
        else
        {
          dJndV = dJndV2Vec[iedge];
          dJpdV = dJpdV2Vec[iedge];
          dJndn = dJndn2Vec[iedge];
          dJpdp = dJpdn2Vec[iedge];
        }

        Vcoef += (sign* dJndV + sign* dJpdV)*ilen;
        Ncoef += (sign* dJndn)*ilen;
        Pcoef += (sign* dJpdp)*ilen;

      } // end of nn edge loop

      col1 = iterDI->Vcol[iVcol];
      if (col1 != -1)
      {
        // find this column:
        bmatch = false;
        int size = iterDI->dIdXcols.size();
        for (cnt2=0;cnt2<size;++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          msg = "pdTerminalCurrents:  Could not find a column match in dIdXcols";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
        }

        iterDI->dIdX[cnt2] += Vcoef*scalingVars.J0*scalingVars.a0;

        Xyce::dout() << iterDI->eName;
      }

      col1 = iterDI->Ncol[iNcol];
      if (col1 != -1)
      {
        // find this column:
        bmatch = false;
        int size = iterDI->dIdXcols.size();
        for (cnt2=0;cnt2<size;++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          msg = "pdTerminalCurrents:  Could not find a column match in dIdXcols";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
        }

        iterDI->dIdX[cnt2] += Ncoef*scalingVars.J0*scalingVars.a0;
      }

      col1 = iterDI->Pcol[iPcol];
      if (col1 != -1)
      {
        // find this column:
        bmatch = false;
        int size = iterDI->dIdXcols.size();
        for (cnt2=0;cnt2<size;++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          msg = "pdTerminalCurrents:  Could not find a column match in dIdXcols";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
        }

        iterDI->dIdX[cnt2] += Pcoef*scalingVars.J0*scalingVars.a0;
      }
      ++iVcol;
      ++iNcol;
      ++iPcol;

      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iVcol,++iNcol,++iPcol)
      {
        iedge = iterEI->iedge;
        inodeB = iterEI->inode;
        sign = (*iterI < inodeB) ? 1.0 : -1.0;
        mEdge * edgePtr = meshContainerPtr->getEdge(iedge);
        ilen = edgePtr->ilen;

        if (*iterI>inodeB)
        {
          dJndV = dJndV1Vec[iedge];
          dJpdV = dJpdV1Vec[iedge];
          dJndn = dJndn1Vec[iedge];
          dJpdp = dJpdn1Vec[iedge];
        }
        else
        {
          dJndV = dJndV2Vec[iedge];
          dJpdV = dJpdV2Vec[iedge];
          dJndn = dJndn2Vec[iedge];
          dJpdp = dJpdn2Vec[iedge];
        }

        Vcoef = (sign* dJndV + sign* dJpdV)*ilen;
        Ncoef = (sign* dJndn)*ilen;
        Pcoef = (sign* dJpdp)*ilen;

        col1 = iterDI->Vcol[iVcol];
        if (col1 != -1)
        {
          // find this column:
          bmatch = false;
          int size = iterDI->dIdXcols.size();
          for (cnt2=0;cnt2<size;++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            msg = "pdTerminalCurrents:  Could not find a column match in dIdXcols";
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
          }

          iterDI->dIdX[cnt2] += Vcoef*scalingVars.J0*scalingVars.a0;
        }

        col1 = iterDI->Ncol[iNcol];
        if (col1 != -1)
        {
          // find this column:
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            msg = "pdTerminalCurrents:  Could not find a column match in dIdXcols";
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
          }

          iterDI->dIdX[cnt2] += Ncoef*scalingVars.J0*scalingVars.a0;
        }

        col1 = iterDI->Pcol[iPcol];
        if (col1 != -1)
        {
          // find this column:
          bmatch = false;
          int size = iterDI->dIdXcols.size();
          for (cnt2=0;cnt2<size;++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            msg = "pdTerminalCurrents:  Could not find a column match in dIdXcols";
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
          }

          iterDI->dIdX[cnt2] += Pcoef*scalingVars.J0*scalingVars.a0;
        }
      } // end of nn edge loop
    } // end of node loop
  }  // end of DI loop

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
    {
      int size = iterDI->dIdXcols.size();
      int size2= iterDI->dIdX.size ();
      Xyce::dout() << "dIdX for electrode: " << iterDI->eName << std::endl;
      for (int ididx=0;ididx<size;++ididx)
      {
        Xyce::dout() << "\t"<< iterDI->dIdXcols[ididx];
        Xyce::dout() << "\t"<< iterDI->dIdX[ididx] << std::endl;
      }
    }
    Xyce::dout() << "Done with Instance::pdTerminalCurrents" << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcTerminalCharges
// Purpose       : Calculates total terminal charge for each electrode.
//
// Special Notes : This is used in capacitance extraction calculations.
//
//                 The charge will be equal to surface integral of the
//                 normal electric field.
//
//                 This function is still not quite finished.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/05/03
//-----------------------------------------------------------------------------
bool Instance::calcTerminalCharges ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "calcTerminalCharges.  name = " << getName() << std::endl;
  }
#endif

  // loop over the device interface nodes, sum the currents going into each one.
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end   ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI=firstDI; iterDI!=lastDI; ++iterDI)
  {
    // loop over the nodes of this device interface node:

    if ( !( meshContainerPtr->labelEdgeType (iterDI->eName) ) ) continue;

    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

    std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
    std::vector<int>::iterator iterI;

    iterDI->chargeSum = 0.0;

    int nodeIndex;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {

      // loop over neighbor nodes/edges to get charge sum for this node.
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);

      std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
      std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end  ();
      std::vector<EDGEINFO>::iterator iterEI;

      double sum       = 0.0; // total charge for this node

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << " --------------- " << std::endl;
        Xyce::dout() << "name = " << iterDI->eName;
        Xyce::dout() << "  node      = " << *iterI << std::endl;
      }
#endif

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        int iedge = iterEI->iedge;
        int neighbor = iterEI->inode;

        mEdge * edgePtr = meshContainerPtr->getEdge(iedge);
        double ilen = edgePtr->ilen;

        double sign = (*iterI < neighbor) ? +1.0 : -1.0;

        double contrib = sign * eSi * e0 * scalingVars.E0 * EfieldVec[iedge] * ilen;
        sum += contrib;

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "neighbor = "<< neighbor;
          Xyce::dout() << " Efield = " << EfieldVec[iedge];
          Xyce::dout() << " contrib = " << contrib << std::endl;
        }
#endif
      }

      double tmp = scalingVars.a0;
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << " sum*scalingVars.a0 = "<< sum*tmp << std::endl;
        Xyce::dout() << " sum    = " << sum << std::endl;
        Xyce::dout() << " --------------- " << std::endl;
      }
#endif
      // total scaled charge mult. by total scaled area for this node.
      sum *= tmp;

      iterDI->chargeSum += sum;

    } // iterI loop.

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << "  " << iterDI->eName;
      Xyce::dout() << " terminal charge = " << iterDI->chargeSum << std::endl;
    }
#endif
  }  // iterDI loop.


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdTerminalCharges
//
// Purpose       : Sets up derivatives used in calculating capacitance.
//
//                 These quantities are not used for the 2-level Newton,
//                 but can be used to obtain interesting information.  Just
//                 like the 2-level Newton is used to obtain lumped
//                 parameter conductances, this stuff is analogously used
//                 to obtain lumped parameter capacitances.
//
//          dQdVckt - derivative of terminal charge w.r.t. Vckt.
//                  This is analogous to dIdVckt, from pdTerminalCurrents.
//
//          dQdX  - (row vector) derivative of the terminal charge
//                  w.r.t. the vector of PDE solution variables. (ie not
//                  including Vckt, as that is not part of the PDE
//                  domain).  This is a vector quantity.
//
//                  This is analogous to the quantity dIdX, which is
//                  calculated in the function pdTerminalCurrents.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/05/03
//-----------------------------------------------------------------------------
bool Instance::pdTerminalCharges ()
{
  bool bsuccess = true;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  std::string msg;

  int nodeIndex;
  int iedge  = 0;
  int inodeB = 0;
  double sign = 1.0;
  double dEdV = 0.0;
  double coef;

  if (!pdTermsAllocated)
  {
    allocatePDTerms ();
    pdTermsAllocated = true;
  }

  // first calculate dQdVckt.----------------------------------------
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    iterDI->dQdVckt = 0.0;

    // if using the old BC, dQdVckt is just zero.
#ifdef Xyce_NEW_BC

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << Xyce::subsection_divider << std::endl;
    }
#endif

    if (iterDI->gid ==-1) continue;

    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

    // obtain the node indices for the current label, loop over them.
    std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
    std::vector<int>::iterator iterI;

    std::vector<EDGEINFO>::iterator firstEI;
    std::vector<EDGEINFO>::iterator lastEI;
    std::vector<EDGEINFO>::iterator iterEI;

    // for the "new" boundary conditions, the currents coming into
    // the electrode have a direct dependency on the voltage at
    // the circuit node, so they all contribute to the dependence
    // of that KCL equation on that node's voltage.  The contributions
    // here are all the equivalent of those that disappear due to
    // the removal of the boundary mesh nodes from the system
    // of equations.
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // do the center point first.
      coef = 0.0;
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
      {
        iedge = iterEI->iedge;
        inodeB = iterEI->inode;

        mEdge * edgePtr = meshContainerPtr->getEdge(iedge);
        double elen = edgePtr->elen;
        double ilen = edgePtr->ilen;

        // if this is a boundary edge, skip it, because the neighbor point
        // is not part of the solution vector.  Only use
        // those which extend into the interior.
        if (boundarySten[inodeB]) continue;

        //sign = (*iterI < inodeB) ? 1.0 : -1.0;
        sign = 1.0;

        dEdV = (1.0/elen);

        coef += sign* dEdV * ilen;
      } // end of nn edge loop

      double tmpsum = (scalingVars.rV0)*eSi*e0*coef*scalingVars.E0*scalingVars.a0;

      iterDI->dQdVckt += tmpsum;

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        Xyce::dout().setf(std::ios::left);
        Xyce::dout() << iterDI->eName<<":";
        Xyce::dout() << " KCL pdTerminalCharges row = " << iterDI->gid;
        Xyce::dout() << " contrib = " << tmpsum;
        Xyce::dout() << " dQdVckt = " << iterDI->dQdVckt << std::endl;
      }
#endif

    } // end of node loop

#endif // Xyce_NEW_BC

  }  // end of DI loop

#ifdef Xyce_NEW_BC
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::subsection_divider << std::endl;
  }
#endif
#endif // Xyce_NEW_BC

  // now do dQdX.-------------------------------------------------
  double Vcoef;
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

    int numdQdX = iterDI->dQdX.size ();
    for (int j=0;j<numdQdX;++j)
    {
      iterDI->dQdX[j] = 0.0;
    }

    // obtain the node indices for the current label, loop over them.
    // for each edge node, add an extra column entry to the colarray.

    std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
    std::vector<int>::iterator iterI;

    std::vector<EDGEINFO>::iterator firstEI;
    std::vector<EDGEINFO>::iterator lastEI;
    std::vector<EDGEINFO>::iterator iterEI;

    int iVcol = 0;
    int iNcol = 0;
    int iPcol = 0;
    int col1;
    int cnt2;
    bool bmatch;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // do the center point first.
      Vcoef = 0.0;
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
      {
        iedge = iterEI->iedge;
        inodeB = iterEI->inode;
        //sign = (*iterI < inodeB) ? 1.0 : -1.0;
        sign = 1.0;

        mEdge * edgePtr = meshContainerPtr->getEdge(iedge);
        double elen = edgePtr->elen;
        double ilen = edgePtr->ilen;

        dEdV = (1.0/elen);
        Vcoef += sign* dEdV * ilen;

      } // end of nn edge loop

      col1 = iterDI->Vcol[iVcol];
      if (col1 != -1)
      {
        // find this column:
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          msg = "pdTerminalCharges:  Could not find a column match in dIdXcols";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
        }
        iterDI->dQdX[cnt2] += Vcoef*eSi*e0*scalingVars.E0*scalingVars.a0;
      }

      ++iVcol;

      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iVcol,++iNcol,++iPcol)
      {
        iedge = iterEI->iedge;
        inodeB = iterEI->inode;
        //sign = (*iterI < inodeB) ? 1.0 : -1.0;
        sign = -1.0;

        mEdge * edgePtr = meshContainerPtr->getEdge(iedge);
        double elen = edgePtr->elen;
        double ilen = edgePtr->ilen;

        dEdV = (1.0/elen);
        Vcoef = sign* dEdV * ilen;

        col1 = iterDI->Vcol[iVcol];
        if (col1 != -1)
        {
          // find this column:
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            msg = "pdTerminalCharges:  Could not find a column match in dIdXcols";
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
          }

          iterDI->dQdX[cnt2] += Vcoef*eSi*e0*scalingVars.E0*scalingVars.a0;
        }
      } // end of nn edge loop
    } // end of node loop
  }  // end of DI loop

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
    {
      int size = iterDI->dIdXcols.size();
      int size2= iterDI->dQdX.size ();
      Xyce::dout() << "dQdX for electrode: " << iterDI->eName << std::endl;
      for (int ididx=0;ididx<size;++ididx)
      {
        Xyce::dout() << "\t"<< iterDI->dIdXcols[ididx];
        Xyce::dout() << "\t"<< iterDI->dQdX[ididx] << std::endl;
      }
    }
    Xyce::dout() << "Done with Instance::pdTerminalCharges" << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcDXDV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/02/02
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
// Creation Date : 12/02/02
//-----------------------------------------------------------------------------
bool Instance::loadDFDV (int ielectrode, N_LAS_Vector * dfdvPtr)
{
  bool bsuccess = true;
  bool bs1 = true;
  DeviceInterfaceNode & dINode = dIVec[ielectrode];
  N_LAS_Vector & dfdv = *dfdvPtr;

  int numNeighbor = dINode.neighborNodes.size();

  double coef;
  int iNeighbor;
  int dFdVindex = 0;
  for (iNeighbor=0;iNeighbor<numNeighbor;++iNeighbor)
  {
    int inode = dINode.neighborNodes[iNeighbor];
    int Vrow = Vrowarray[inode];
    int Nrow = Nrowarray[inode];
    int Prow = Prowarray[inode];

    mNode * nodePtr = meshContainerPtr->getNode(inode);

    int iNN=0;
    for (iNN=0;iNN<nodePtr->cnode;++iNN)
    {
      int inodeB = nodePtr->edgeInfoVector[iNN].inode;

      // if nodeB is not a boundary node, never mind.
      if (boundarySten[inodeB]!=1) continue;

      // if it is a boundary node, but not part of the
      // current boundary, also never mind.
      if (labelNameVector[inodeB]!= dINode.eName) continue;

      // load V term:
      coef = dINode.dFdVckt[dFdVindex];
      if( useVectorGIDFlag )
      {
        if (Vrow != -1)
        {
          bs1 = dfdv.setElementByGlobalIndex(Vrow, -coef, 0);
          bsuccess = bsuccess && bs1;
        }
      }
      else
      {
        dfdv[Vrow] = -coef;
      }

      ++dFdVindex;

      // load N term:
      coef = dINode.dFdVckt[dFdVindex];
      if( useVectorGIDFlag )
      {
        if (Nrow != -1)
        {
          bs1 = dfdv.setElementByGlobalIndex(Nrow, -coef, 0);
          bsuccess = bsuccess && bs1;
        }
      }
      else
      {
        dfdv[Nrow] = -coef;
      }
      ++dFdVindex;


      // load P term:
      coef = dINode.dFdVckt[dFdVindex];
      if( useVectorGIDFlag )
      {
        if (Prow != -1)
        {
          bs1 = dfdv.setElementByGlobalIndex(Prow, -coef, 0);
          bsuccess = bsuccess && bs1;
        }
      }
      else
      {
        dfdv[Prow] = -coef;
      }
      ++dFdVindex;
    }
  }

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
//   05/06/03:     Adapting this function to also calculate terminal
//                 capacitances.  The form of this calculation is exactly
//                 the same, but instead of using terminal currents to get
//                 conductances, I'm using terminal charges to get
//                 capacitances.  (dQ/dV instead of dI/dV).
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/02
//-----------------------------------------------------------------------------
bool Instance::calcConductance (int iElectrode, const N_LAS_Vector * dxdvPtr)
{
  bool bsuccess = true;
  const N_LAS_Vector & dxdv = *dxdvPtr;

  calcConductanceFlag = true;

#ifdef Xyce_DEBUG_DEVICE
  char filename1[256];

  for (int ich = 0; ich < 256; ++ich)
  { filename1[ich] = 0; }

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "calcConductances  name = " << getName() << std::endl;
    Xyce::dout() << "electrode = " << dIVec[iElectrode].eName;
    Xyce::dout() << "  dIdVckt = " << dIVec[iElectrode].dIdVckt;
    Xyce::dout() << std::endl;
    Xyce::dout() << std::endl;
  }
#endif

  if (!(dIVec[iElectrode].dxdvAllocated))
  {
    dIVec[iElectrode].dxdvPtr = extData.lasSysPtr->builder().createVector();
    dIVec[iElectrode].dxdvAllocated = true;
  }

  // A linear solve should have just been performed up in the Newton
  // solver.  The result of that solve, dxdv was placed in the RHS vector.
  // dxdv is needed later, so save a copy.
  *(dIVec[iElectrode].dxdvPtr) = *(dxdvPtr);

  // doing the iElectrode Column of the condVec array.
  // This should correspond to the dIVec[iElectrode].gid column of the
  // Jacobian.

  double Gij = 0.0;
  double dIidVj = 0.0;
  double dIidVj_chain = 0.0;  // from the dot product.

  double Cij = 0.0;
  double dQidVj = 0.0;
  double dQidVj_chain = 0.0;  // from the dot product.

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

    if (iElectrode != iEqu)
    {
      dIidVj = 0.0;
      dQidVj = 0.0;
    }
    else
    {
      dIidVj = dIVec[iEqu].dIdVckt;
      dQidVj = dIVec[iEqu].dQdVckt;
    }

    // load dIdX, dQdX:
    extData.tmpdIdXPtr->putScalar (0.0);
    extData.tmpdQdXPtr->putScalar (0.0);
    int DIDXSize = dIVec[iEqu].dIdX.size();
    for (int iDIDX=0;iDIDX<DIDXSize;++iDIDX)
    {
      int index = dIVec[iEqu].dIdXcols[iDIDX];
      double coefI = dIVec[iEqu].dIdX[iDIDX];
      double coefQ = dIVec[iEqu].dQdX[iDIDX];

      if( useVectorGIDFlag )
      {
        if (index != -1)
        {
          extData.tmpdIdXPtr->setElementByGlobalIndex (index, coefI, 0);
          extData.tmpdQdXPtr->setElementByGlobalIndex (index, coefQ, 0);
        }
      }
      else
      {
        (*(extData.tmpdIdXPtr))[index] = coefI;
        (*(extData.tmpdQdXPtr))[index] = coefQ;  // CHECK THIS!
      }                                         // this is probably wrong.
                                                // I think this DVA
      // implementation might break in parallel.
    }

#ifdef Xyce_DEBUG_DEVICE
    sprintf(filename1,"dIdX%02d.txt", iEqu);
    extData.tmpdIdXPtr->writeToFile(filename1);

    sprintf(filename1,"dQdX%02d.txt", iEqu);
    extData.tmpdQdXPtr->writeToFile(filename1);
#endif

    // get dot product of dXdv and dIdX:
    dIidVj_chain = dxdv.dotProduct( *(extData.tmpdIdXPtr) );

    // total conductance:
    Gij = dIidVj_chain + dIidVj;
    condVec[iEqu][iElectrode] = Gij;


    // get dot product of dXdv and dQdX:
    dQidVj_chain = dxdv.dotProduct( *(extData.tmpdQdXPtr) );

    // total capacitance:
    Cij = dQidVj_chain + dQidVj;
    capVec[iEqu][iElectrode] = Cij;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      char outstring[128];
      double Itmp = dIVec[iEqu].currentSum;
      double Vtmp = dIVec[iEqu].Vckt - dIVec[iElectrode].Vckt;
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
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
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
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  updateIntermediateVars ();
  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  for (; iterDI!=lastDI;++iterDI)
  {
    if( useVectorGIDFlag )
    {
      if (iterDI->stateC_owned)
      {
        // place this value for the current sum in the state vector.
        (staVectorPtr)->setElementByGlobalIndex( iterDI->stateC,
                                                 iterDI->currentSum, 0);
      }
    }
    else
    {
      (*staVectorPtr)[iterDI->li_stateC] = iterDI->currentSum;
    }
  }

  // Now store the dielectric displacement in the state vector
  // in order to calculate displacement current.
  int i;
#ifdef Xyce_OLD_DISPLACEMENT

  for (i = 0; i< numMeshEdges; ++i)
  {
    double D = eSi * e0 * scalingVars.E0 * EfieldVec[i];

    if( useVectorGIDFlag )
    {
      if( stateDispl_owned[i] )
      {
        // place this value for the charge in the state vector.
        (staVectorPtr)->setElementByGlobalIndex( stateDispl[i], D, 0);
      }
    }
    else
    {
      (*staVectorPtr)[li_stateDispl[i]] = D;
    }
  }

#else

  // ERK NOTE: as VVEC is already synchronized with the solution vector,
  // it is probably unneccessary to also store it in state.  Fix later.
  for (i = 0; i< numMeshPoints; ++i)
  {
    if( useVectorGIDFlag )
    {
      if( stateDispl_owned[i] )
      {
        // place this value for the charge in the state vector.
        (staVectorPtr)->setElementByGlobalIndex( stateDispl[i], (scalingVars.V0*VVec[i]), 0);
      }
    }
    else
    {
      (*staVectorPtr)[li_stateDispl[i]] = scalingVars.V0 * VVec[i];
    }
  }

#endif


  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;

  // If this is a "ckt-only , don't need to recalculate anything.
  if (getSolverState().twoLevelNewtonCouplingMode==OUTER_PROBLEM) return bsuccess;

  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  for (; iterDI!=lastDI;++iterDI)
  {
    // obtain this value from the state vector.
    if( useVectorGIDFlag )
    {
      iterDI->currentSum =
        (staVectorPtr)->getElementByGlobalIndex( iterDI->stateC, 0);
    }
    else
    {
      iterDI->currentSum = (*staVectorPtr)[iterDI->li_stateC];
    }
  }

  // Now get the displacement current:
  int i;
  double dcmax = 0.0;

#ifdef Xyce_OLD_DISPLACEMENT
  for (i = 0; i< numMeshEdges; ++i)
  {
    N_LAS_Vector * staDerivPtr = extData.nextStaDerivVectorPtr;
    if( useVectorGIDFlag )
    {
      displCurrent[i] = staDerivPtr->getElementByGlobalIndex(stateDispl[i], 0);
    }
    else
    {
      displCurrent[i] = (*staDerivPtr)[li_stateDispl[i]];
    }

    if (fabs(displCurrent[i]) > dcmax) dcmax = fabs(displCurrent[i]);
  }

#else

  // first get the "displacement potential"
  for (i = 0; i< numMeshPoints; ++i)
  {
    N_LAS_Vector * staDerivPtr = extData.nextStaDerivVectorPtr;
    if( useVectorGIDFlag )
    {
      displPotential[i] = staDerivPtr->getElementByGlobalIndex(stateDispl[i], 0);
    }
    else
    {
      displPotential[i] = (*staDerivPtr)[li_stateDispl[i]];
    }

    if (fabs(displCurrent[i]) > dcmax) dcmax = fabs(displCurrent[i]);
  }

  // now use it to get the dielectric displacement current.  This should look pretty
  // similar to the calcEfield function from this point on.
  for (i=0;i<numMeshEdges;++i)
  {
    mEdge * edgePtr = meshContainerPtr->getEdge(i);

    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;
    double elen   = edgePtr->elen;

    displCurrent[i] = -(displPotential[inodeB] - displPotential[inodeA])/elen;
    displCurrent[i] *= matSupport.getRelPerm(bulkMaterial) * e0;

    if (fabs(displCurrent[i]) > dcmax) dcmax = fabs(displCurrent[i]);
  }

#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    Xyce::dout() << "  Maximum displacement current:  " << dcmax << std::endl;
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setInitialGuess
//
// Purpose       : This function sets up a bunch of initial information,
//                 that only has to be set up once at the beginning of the
//                 simulation.  That includes the
//                 simple, analytic initial guess, the density boundary
//                 conditions, carrier lifetimes, scaling variables, etc.
//
// Special Notes : Mobilities are set up here, because originally this code
//                 only had mobilities that were doping-dependent ONLY.
//                 Now the code has mobilities that are also dependent on
//                 carrier density, which can change throughout the
//                 simulation.  So, mobilities are also calculated
//                 elsewhere, but it didn't hurt to leave them here.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool Instance::setInitialGuess ()
{
  bool bsuccess = true;
  bool bs1 = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::setInitialGuess\n";
  }
#endif

  // The stuff inside of the "called before" if statement needs to be here,
  // rather than in the instance constructor, because some of these
  // functions need for the GID's to have been setup.
  if (!calledBeforeSIGB)
  {
    bs1 = calcDensityBCs   ();  bsuccess = bsuccess && bs1;
    bs1 = calcVequBCs      ();  bsuccess = bsuccess && bs1;
    bs1 = calcInitialGuess ();  bsuccess = bsuccess && bs1;
    bs1 = calcMobilities   ();  bsuccess = bsuccess && bs1;
    bs1 = calcLifetimes    ();  bsuccess = bsuccess && bs1;
    bs1 = scaleVariables   ();  bsuccess = bsuccess && bs1;
    calledBeforeSIGB = true;

#ifdef Xyce_DEBUG_DEVICE
    // if running with debug turned on, output the initial guess:
    if (getDeviceOptions().debugLevel > 3 && getSolverState().debugTimeFlag && tecplotLevel > 0)
    {
      outputTecplot        ();
      if (tecplotLevel > 2) outputTecplotVectors ();
    }

    if (getDeviceOptions().debugLevel > 3 && getSolverState().debugTimeFlag && sgplotLevel > 0)
    {
      outputSgplot        ();
    }

    if (getDeviceOptions().debugLevel > 3 && getSolverState().debugTimeFlag && gnuplotLevel > 0)
    {
      outputGnuplot        ();
    }
#endif
  }  // calledBeforeSIGB

  return (bsuccess);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadVecNLPoisson
// Purpose       :
//
// Special Notes : This function returns the nonlinear poisson rhs load,
//                 multiplied by a scalar.  The scalar generally will
//                 be 1.0 or -1.0.  -1.0 is used for the new-DAE formualtion,
//                 while 1.0 is the normal thing to do for old-DAE.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/24/05
//-----------------------------------------------------------------------------
bool Instance::loadVecNLPoisson (double scalar, N_LAS_Vector * vecPtr)
{
  bool bsuccess = true;
  bool bs1 = true;
  std::string semi(bulkMaterial);
  int i;
  int Vrow, Nrow, Prow;
  double coef, coef2;
#ifndef Xyce_NEW_BC
  double vtmp;
#endif
  double ilen, elen, nodeArea;
  double holeDens;
  double elecDens;

  Ut = Vt/scalingVars.V0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::loadVecNLPoisson\n";
    Xyce::dout() << "       name = " << getName() <<"\n";

    Xyce::dout() << "      Vt    = " << Vt << "\n";
    Xyce::dout() << "      Ut    = " << Ut << "\n";
    Xyce::dout() << "      scalingVars.V0    = " << scalingVars.V0 << "\n";

  }
#endif

  // mesh points for the PDE problem:
  for (i=0;i<numMeshPoints;++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "--------" << std::endl;
      Xyce::dout() << "Mesh Point i =  " << i;
      Xyce::dout() << "  x = " << xVec[i] *scalingVars.x0;
      Xyce::dout() << "  y = " << yVec[i] *scalingVars.x0 << std::endl;
    }
#endif

    if( useVectorGIDFlag )
    {
      Vrow = Vrowarray[i];
      Nrow = Nrowarray[i];
      Prow = Prowarray[i];
    }
    else
    {
      Vrow = li_Vrowarray[i];
      Nrow = li_Nrowarray[i];
      Prow = li_Prowarray[i];
    }

    // Is this node a node with an explicit boundary condition?
    // If so, apply the BC.  Otherwise, just load the region operators.

#ifdef Xyce_NEW_BC
    // If we are using the "new" boundary conditions, the are simply
    // imposed, rather than being solved, so don't do anything here.
    if (boundarySten[i]) continue;

#else // "old" boundary condition:
    bool doneFlag = false;

    if (boundarySten[i]) // do BC load
    {
      int DIindex = labelDIMap[labelNameVector[i]];
      if (dIVec[DIindex].given)
      {
        vtmp = 0.0;
        if (vOwnVec[i] == 1) vtmp = VVec[i];

        int i1 = dIVec[DIindex].meshGlobalToLocal[i];
        coef = vtmp - dIVec[DIindex].VbcVec[i1];

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "BC  load:  Vrow = " << Vrow << " coef = " << coef << std::endl;
          Xyce::dout() << "    vtmp = " << vtmp << std::endl;
          Xyce::dout() << "    Vckt = " << dIVec[DIindex].Vckt << std::endl;
        }
#endif
        if( useVectorGIDFlag )
        {
          bs1 = vecPtr->sumElementByGlobalIndex(Vrow, -scalar*coef, 0);
          bsuccess = bsuccess && bs1;
          bs1 = vecPtr->sumElementByGlobalIndex(Nrow,  0.0 , 0);
          bsuccess = bsuccess && bs1;
          bs1 = vecPtr->sumElementByGlobalIndex(Prow,  0.0 , 0);
          bsuccess = bsuccess && bs1;
        }
        else
        {
          (*vecPtr)[Vrow] += -scalar*coef;
          (*vecPtr)[Nrow] += 0.0;
          (*vecPtr)[Prow] += 0.0;
        }

        doneFlag = true;

      } // inner if
    } // outer if

    if (doneFlag) continue;
#endif // Xyce_NEW_BC

    // if load is not done yet, then do an interior point load:
    mNode * nodePtr = meshContainerPtr->getNode(i);
    nodeArea = nodePtr->area;
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "--------" << std::endl;
      Xyce::dout() << "Interior:  mesh point: " << i << "  Vrow = " << Vrow;
      Xyce::dout() << std::endl;
    }
#endif

    coef = 0.0;
    for (int iNN=0;iNN<nodePtr->cnode;++iNN)
    {
      int iedge  = nodePtr->edgeInfoVector[iNN].iedge;
      int inodeB = nodePtr->edgeInfoVector[iNN].inode;
      ilen = nodePtr->edgeInfoVector[iNN].ilen;
      elen = nodePtr->edgeInfoVector[iNN].elen;

      double efield_loc = (VVec[i]-VVec[inodeB])/elen;
      coef +=  -efield_loc * ilen ;
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "----" << std::endl;
        Xyce::dout() << "  iedge  = " << iedge << std::endl;
        Xyce::dout() << "  Vlocal = " << VVec[i] << std::endl;
        Xyce::dout() << "  Vneigh = " << VVec[inodeB] << std::endl;
        Xyce::dout() << "  efield = " << efield_loc << std::endl;
        Xyce::dout() << "  elen   = " << elen << std::endl;
        Xyce::dout() << "  ilen   = " << ilen << std::endl;
        Xyce::dout() << "  inodeB = " << inodeB;
        Xyce::dout() << "  x[b] = " << xVec[inodeB]*scalingVars.x0;
        Xyce::dout() << "  y[b] = " << yVec[inodeB]*scalingVars.x0 << std::endl;
        Xyce::dout() << std::endl;
      }
#endif
    }
    coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;

    holeDens = getVoltDepHoleDens ( VminExp, VVec[i], NpMax);
    elecDens = getVoltDepElecDens ( VmaxExp, VVec[i], NnMax);
    coef2    = -(holeDens-elecDens+CVec[i]);

#ifdef Xyce_OXIDE_ENABLED
    if (!allOxideFlag)
    {
      coef += coef2;
    }
#else
    coef += coef2;
#endif

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "--------" << std::endl;
      Xyce::dout() << "  holeD = " << holeDens << std::endl;
      Xyce::dout() << "  elecD = " << elecDens << std::endl;
      Xyce::dout() << "  dopeD = " << CVec[i] << std::endl;
      Xyce::dout() << "  coef2 = " << coef2 << std::endl;
      Xyce::dout() << "  scalingVars.L0    = " << scalingVars.L0 * matSupport.getRelPerm(semi)   << std::endl;
      Xyce::dout() << "  nodeArea  = " << nodeArea  << std::endl;
      Xyce::dout() << "  coef  = " << coef  << std::endl;
      Xyce::dout() << "--------" << std::endl;
    }
#endif

    if( useVectorGIDFlag )
    {
      bs1 = vecPtr->sumElementByGlobalIndex(Vrow, -scalar*coef, 0);
      bsuccess = bsuccess && bs1;
      bs1 = vecPtr->sumElementByGlobalIndex(Nrow,  0.0 , 0);
      bsuccess = bsuccess && bs1;
      bs1 = vecPtr->sumElementByGlobalIndex(Prow,  0.0 , 0);
      bsuccess = bsuccess && bs1;
    }
    else
    {
      (*vecPtr)[Vrow] += -scalar*coef;
      (*vecPtr)[Nrow] += 0.0;
      (*vecPtr)[Prow] += 0.0;
    }

  } // mesh point loop

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    Xyce::dout() << section_divider << std::endl;
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadVecDDForm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/24/05
//-----------------------------------------------------------------------------
bool Instance::loadVecDDForm
(double scalar, double dndtScalar, N_LAS_Vector * vecPtr)
{
  bool bsuccess = true;
  bool bs1 = true;
  std::string semi(bulkMaterial);
  int i;
  int iNN;
  int Vrow, Nrow, Prow;
  double coef, coef2;
  double ilen, elen, nodeArea;
  double holeDens, elecDens;
#ifndef Xyce_NEW_BC
  double vtmp;
  double ntmp;
  double ptmp;
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\n"<<section_divider << std::endl;
    Xyce::dout() << "Instance::loadVecDDForm\n";
    Xyce::dout() << "         name = " << getName()  << "\n";
  }
#endif

  // KCL equations for the various connecting terminals:
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  // if this is the inner loop of a multilevel Newton solve, don't do the
  // KCL-related loads.
  if ( !(getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM))
  {
    for (;iterDI!=lastDI;++iterDI)
    {
      coef = iterDI->currentSum;

      if( useVectorGIDFlag )
      {
        if (iterDI->gid != -1)
        {
          bs1 = vecPtr->sumElementByGlobalIndex(iterDI->gid, -scalar*coef, 0);
          bsuccess = bsuccess && bs1;
        }
      }
      else
      {
        (*vecPtr)[iterDI->lid] += -scalar*coef;
      }

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "KCL for "<< iterDI->eName << ":\n";
        Xyce::dout() << " row = " << iterDI->gid << "\n";
        Xyce::dout() << "coef = " << coef << "\n";
      }
#endif

    }
  } // end of twoLevelNewtonCouplingMode if statement.

  // mesh points for the PDE problem:
  for (i=0;i<numMeshPoints;++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "--------" << std::endl;
      Xyce::dout() << "Mesh Point i =  " << i;
      Xyce::dout() << "  x = " << xVec[i]*scalingVars.x0;
      Xyce::dout() << "  y = " << yVec[i]*scalingVars.x0 << std::endl;
    }
#endif

    if( useVectorGIDFlag )
    {
      Vrow = Vrowarray[i];
      Nrow = Nrowarray[i];
      Prow = Prowarray[i];
    }
    else
    {
      Vrow = li_Vrowarray[i];
      Nrow = li_Nrowarray[i];
      Prow = li_Prowarray[i];
    }

    bool doneFlagV = false;
    bool doneFlagN = false;
    bool doneFlagP = false;
#ifdef Xyce_NEW_BC
    // if doing the "new" boundary conditions, then just skip the load.
    // The boundary conditions are just imposed.
    // If using this type of BC, Neumann conditions are not an option.
    // Only Dirichlet.  (not mixed either)
    if (boundarySten[i]) continue;
#else
    bool doneFlag  = false;

    // Is this node a node with an explicit DIRICHLET boundary condition?
    // If so, apply the BC.  Otherwise, just load the region operators.
    // If the BC has been specified to be neumann, then don't do anything
    // unusual - the region operator will enforce a neumann condition.

    if (boundarySten[i]) // do BC load
    {
      int DIindex = labelDIMap[labelNameVector[i]];

      if (dIVec[DIindex].given)
      {
        int ilocal = dIVec[DIindex].meshGlobalToLocal[i];

        if (dIVec[DIindex].neumannBCFlagV==false)
        {
          vtmp = 0.0;
          if (vOwnVec[i] == 1) vtmp = VVec[i];

          coef = vtmp - dIVec[DIindex].VbcVec[ilocal];

          if( useVectorGIDFlag )
          {
            bs1 = vecPtr->sumElementByGlobalIndex(Vrow, -scalar*coef, 0);
            bsuccess = bsuccess && bs1;
          }
          else
          {
            (*vecPtr)[Vrow] += -scalar*coef;
          }

#ifdef Xyce_DEBUG_DEVICE
          if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
          {
            Xyce::dout() << "BC  load:  Vrow = " << Vrow;
            Xyce::dout() << " coef = " << coef << std::endl;
            Xyce::dout() << "    vtmp = " << vtmp << std::endl;
            Xyce::dout() << "    Vckt = " << dIVec[DIindex].Vckt << std::endl;
            Xyce::dout() << "    Vequ = " << dIVec[DIindex].VequVec[ilocal] << std::endl;
          }
#endif
          doneFlagV = true;
        }  // end of neumann if statement.

        if (dIVec[DIindex].neumannBCFlagN==false)
        {
          ntmp = 0.0;
          if (nnOwnVec[i] == 1) ntmp = nnVec[i];
          coef = ntmp - (dIVec[DIindex].nnbcVec[ilocal]);

          if( useVectorGIDFlag )
          {
            bs1 = vecPtr->sumElementByGlobalIndex(Nrow,  -scalar*coef , 0);
            bsuccess = bsuccess && bs1;
          }
          else
          {
            (*vecPtr)[Nrow] += -scalar*coef;
          }

#ifdef Xyce_DEBUG_DEVICE
          if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
          {
            Xyce::dout() << "BC  load:  Nrow = ";
            Xyce::dout() << Nrow << " coef = " << coef << std::endl;
            Xyce::dout() << "    ntmp = " << ntmp << std::endl;
            Xyce::dout() << "    nnbc = " << dIVec[DIindex].nnbcVec[ilocal] << std::endl;
          }
#endif
          doneFlagN = true;
        } // end of neumann if


        if (dIVec[DIindex].neumannBCFlagP==false)
        {
          ptmp = 0.0;
          if (npOwnVec[i] == 1) ptmp = npVec[i];
          coef = ptmp - (dIVec[DIindex].npbcVec[ilocal]);

          if( useVectorGIDFlag )
          {
            bs1 = vecPtr->sumElementByGlobalIndex(Prow,  -scalar*coef , 0);
            bsuccess = bsuccess && bs1;
          }
          else
          {
            (*vecPtr)[Prow] += -scalar*coef;
          }

#ifdef Xyce_DEBUG_DEVICE
          if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
          {
            Xyce::dout() << "BC  load:  Prow = ";
            Xyce::dout() << Prow << " coef = " << coef << std::endl;
            Xyce::dout() << "    ptmp = " << ptmp << std::endl;
            Xyce::dout() << "    npbc = " << dIVec[DIindex].npbcVec[ilocal] << std::endl;
          }
#endif
          doneFlagP = true;
        } // end of neumann if

        doneFlag = (doneFlagV && doneFlagN && doneFlagP);

      } // inner if (was this DI given by the user in the netlist)
    } // outer if  (finding the label in the label name vector)

    if (doneFlag) continue;

#endif // Xyce_NEW_BC

    // Do the interior mesh points, if we've gotten this far:

    mNode * nodePtr = meshContainerPtr->getNode(i);
    nodeArea = nodePtr->area;
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "--------" << std::endl;
      Xyce::dout() << "Interior:  mesh point: " << i << "  Vrow = " << Vrow;
      Xyce::dout() << std::endl;
    }
#endif

    // do poisson's equation first:
    if (!doneFlagV)
    {
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "--------" << std::endl;
        Xyce::dout() << "  Poisson equ:";
      }
#endif

      coef = 0.0;
      for (iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        int iedge  = nodePtr->edgeInfoVector[iNN].iedge;
        int inodeB = nodePtr->edgeInfoVector[iNN].inode;
        ilen = nodePtr->edgeInfoVector[iNN].ilen;
        elen = nodePtr->edgeInfoVector[iNN].elen;

        //coef +=  ((i<inodeB)?1.0:-1.0) *  EfieldVec[iedge] * ilen ;
        double efield_loc = (VVec[i]-VVec[inodeB])/elen;
        coef +=  -efield_loc * ilen ;
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "----" << std::endl;
          Xyce::dout() << "  iedge  = " << iedge << std::endl;
          Xyce::dout() << "  Vlocal = " << VVec[i] << std::endl;
          Xyce::dout() << "  Vneigh = " << VVec[inodeB] << std::endl;
          Xyce::dout() << "  efield = " << efield_loc << std::endl;
          Xyce::dout() << "  elen   = " << elen << std::endl;
          Xyce::dout() << "  ilen   = " << ilen << std::endl;
          Xyce::dout() << "  inodeB = " << inodeB;
          Xyce::dout() << "  x[b] = " << xVec[inodeB]*scalingVars.x0;
          Xyce::dout() << "  y[b] = " << yVec[inodeB]*scalingVars.x0 << std::endl;
          Xyce::dout() << std::endl;
        }
#endif
      }
      coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;

      holeDens = npVec[i];
      elecDens = nnVec[i];
      coef2    = -(holeDens-elecDens+CVec[i]);

#ifdef Xyce_OXIDE_ENABLED
      if (!allOxideFlag)
      {
        coef += coef2;
      }
#else
      coef += coef2;
#endif

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "--------" << std::endl;
        Xyce::dout() << "  holeD = " << holeDens << std::endl;
        Xyce::dout() << "  elecD = " << elecDens << std::endl;
        Xyce::dout() << "  dopeD = " << CVec[i] << std::endl;
        Xyce::dout() << "  coef2 = " << coef2 << std::endl;
        Xyce::dout() << "  scalingVars.L0    = " << scalingVars.L0 * matSupport.getRelPerm(semi)   << std::endl;
        Xyce::dout() << "  nodeArea  = " << nodeArea  << std::endl;
        Xyce::dout() << "  coef  = " << coef  << std::endl;
        Xyce::dout() << "--------" << std::endl;
      }
#endif
      if ( getSolverState().chargeHomotopy )
      {
        coef *= getSolverState().chargeAlpha;
      }

      if( useVectorGIDFlag )
      {
        bs1 = vecPtr->sumElementByGlobalIndex(Vrow, -scalar*coef, 0);
        bsuccess = bsuccess && bs1;
      }
      else
      {
        (*vecPtr)[Vrow] += -scalar*coef;
      }

      doneFlagV = true;

    } // doneFlagV if statement

    // Now do electron continuity
    if (!doneFlagN)
    {
      // get electron time derivative and scale.

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "--------" << std::endl;
        Xyce::dout() << "  Electron equ:";
      }
#endif

      double dndt = 0.0;

#if 0
      if( useVectorGIDFlag )
      {
        if ( !((getSolverState().dcopFlag) || Nrowarray[i] == -1) )
        {
          dndt = (extData.nextSolDerivVectorPtr)->getElementByGlobalIndex
                 (Nrowarray[i],0);
        }
      }
      else
      {
        if ( !(getSolverState().dcopFlag) )
          dndt = (*extData.nextSolDerivVectorPtr)[li_Nrowarray[i]];
      }

      dndt *= scalingVars.t0*dndtScalar; // if new-DAE, dndtScalar is 0.0
#endif

      coef = 0.0;
      for (iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        int iedge  = nodePtr->edgeInfoVector[iNN].iedge;
        int inodeB = nodePtr->edgeInfoVector[iNN].inode;
        ilen = nodePtr->edgeInfoVector[iNN].ilen;
        elen = nodePtr->edgeInfoVector[iNN].elen;

        coef +=  ((i<inodeB)?1.0:-1.0) *  JnVec[iedge] * ilen ;

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "----" << std::endl;
          Xyce::dout() << "  iedge  = " << iedge << std::endl;
          Xyce::dout() << "  Jn     = " << JnVec[iedge] << std::endl;
          Xyce::dout() << "  nlocal = " << nnVec[i] << std::endl;
          Xyce::dout() << "  nneigh = " << nnVec[inodeB] << std::endl;
          Xyce::dout() << "  elen   = " << elen << std::endl;
          Xyce::dout() << "  ilen   = " << ilen << std::endl;
          Xyce::dout() << "  inodeB = " << inodeB;
          Xyce::dout() << "  x[b] = " << xVec[inodeB]*scalingVars.x0;
          Xyce::dout() << "  y[b] = " << yVec[inodeB]*scalingVars.x0 << std::endl;
          Xyce::dout() << std::endl;
        }
#endif
      }
      coef /= nodeArea;
#if 1
      coef += - totSrcVec[i] - dndt;
#else
      coef += - RVec[i] - dndt;
#endif

      coef += elecPenalty[i];

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "--------" << std::endl;
        Xyce::dout().setf(std::ios::left);
        Xyce::dout() << "  row  = " << Nrow;
        Xyce::dout().setf(std::ios::scientific);
        Xyce::dout() << "  coef=" << coef;
        Xyce::dout() << " nodeArea ="<<nodeArea;
        Xyce::dout() << " R[i]="<<RVec[i];
        Xyce::dout() << " dndt="<<dndt;
        Xyce::dout() << "\n";
        Xyce::dout() << "--------" << std::endl;
      }
#endif

      if( useVectorGIDFlag )
      {
        bs1 = vecPtr->setElementByGlobalIndex(Nrow,-scalar*coef,0);
        bsuccess = bsuccess && bs1;
      }
      else
      {
        (*vecPtr)[Nrow] += -scalar*coef;
      }

      doneFlagN = true;

    } // end of doneFlagN if statement.

    // Now do hole continuity
    if (!doneFlagP)
    {
      // get hole time derivative and scale.
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "--------" << std::endl;
        Xyce::dout() << "  Hole equ:";
      }
#endif

      double dpdt = 0.0;
#if 0
      if( useVectorGIDFlag )
      {
        if ( !((getSolverState().dcopFlag) || Prowarray[i] == -1) )
          dpdt = (extData.nextSolDerivVectorPtr)->getElementByGlobalIndex
                 (Prowarray[i],0);
      }
      else
      {
        if ( !(getSolverState().dcopFlag) )
          dpdt = (*extData.nextSolDerivVectorPtr)[li_Prowarray[i]];
      }

      dpdt *= scalingVars.t0*dndtScalar; // if new-DAE, dndtScalar is 0.0.
#endif

      coef = 0.0;
      for (iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        int iedge  = nodePtr->edgeInfoVector[iNN].iedge;
        int inodeB = nodePtr->edgeInfoVector[iNN].inode;
        ilen = nodePtr->edgeInfoVector[iNN].ilen;
        elen = nodePtr->edgeInfoVector[iNN].elen;

        coef +=  ((i<inodeB)?1.0:-1.0) *  JpVec[iedge] * ilen ;

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "----" << std::endl;
          Xyce::dout() << "  iedge  = " << iedge << std::endl;
          Xyce::dout() << "  Jp     = " << JpVec[iedge] << std::endl;
          Xyce::dout() << "  plocal = " << npVec[i] << std::endl;
          Xyce::dout() << "  pneigh = " << npVec[inodeB] << std::endl;
          Xyce::dout() << "  elen   = " << elen << std::endl;
          Xyce::dout() << "  ilen   = " << ilen << std::endl;
          Xyce::dout() << "  inodeB = " << inodeB;
          Xyce::dout() << "  x[b] = " << xVec[inodeB]*scalingVars.x0;
          Xyce::dout() << "  y[b] = " << yVec[inodeB]*scalingVars.x0 << std::endl;
          Xyce::dout() << std::endl;
        }
#endif
      }
      coef /= -nodeArea;

#if 1
      coef += - totSrcVec[i] - dpdt;
#else
      coef += - RVec[i] - dpdt;
#endif
      coef += holePenalty[i];

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "--------" << std::endl;
        Xyce::dout().setf(std::ios::left);
        Xyce::dout() << "  row  = " << Prow;
        Xyce::dout().setf(std::ios::scientific);
        Xyce::dout() << "  coef=" << coef;
        Xyce::dout() << " nodeArea ="<<nodeArea;
        Xyce::dout() << " RVec[i]="<<RVec[i];
        Xyce::dout() << " dpdt="<<dpdt;
        Xyce::dout() << "\n";
        Xyce::dout() << "--------" << std::endl;
      }
#endif

      if( useVectorGIDFlag )
      {
        bs1 = vecPtr->setElementByGlobalIndex(Prow,-scalar*coef,0);
        bsuccess = bsuccess && bs1;
      }
      else
      {
        (*vecPtr)[Prow] += -scalar*coef;
      }

      doneFlagP = true;

    } // end of doneFlagP

  } // ip_iter, row loop...

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    Xyce::dout() << section_divider << std::endl;
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatNLPoisson
// Purpose       :
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/24/05
//-----------------------------------------------------------------------------
bool Instance::loadMatNLPoisson (N_LAS_Matrix * matPtr)
{
  bool bsuccess = true;
  bool bs1 = true;

  int Vrow, Nrow, Prow;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::loadJacNonlinPoisson" << "\n";
    Xyce::dout() << "  name = " << getName() <<"\n";
    Xyce::dout() << "\n";
  }
#endif

  int i,j;
  int count  = 0;
  std::string semi(bulkMaterial);
  Ut = Vt/scalingVars.V0;
  double pre = 1.0/Ut;

  double elecDens;
  double holeDens;

  int iNN;

  double coef;
  double ilen, elen, nodeArea;

#ifdef Xyce_DEBUG_DEVICE
  double q   = charge;
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "pre = " << pre << "\n";
    Xyce::dout() << "Vt  = " << Vt  << "\n";
    Xyce::dout() << "eps = " << eps << "\n";
    Xyce::dout() << "q   = " << q   << "\n";
    Xyce::dout() << "Na  = " << Na  << "\n";
    Xyce::dout() << "Nd  = " << Nd  << "\n";
    Xyce::dout() << "NpMax  = " << NpMax  << "\n";
    Xyce::dout() << "NnMax  = " << NnMax  << "\n";
  }
#endif

  // note: fix the vals and cols arrays later!
  int numCol = cols.size();
  if (vals.size () < cols.size()) numCol = vals.size();

  // set up some of the partial derivative arrays:
  bs1 = pdRecombination ();     bsuccess = bsuccess && bs1;
  bs1 = pdElectronCurrent ();   bsuccess = bsuccess && bs1;
  bs1 = pdHoleCurrent ();       bsuccess = bsuccess && bs1;

  // Load the jacobian, row by row.
  // rows associated with the connecting terminal KCL's:
  for (j=0;j<numCol;++j)
  {
    cols[j] = -1;
    vals[j] = 0.0;
  }

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  // If this PDE device is not coupled (always true for the nonlinear
  // poisson) to the ckt, put some 1's on the
  // diagonal to be safe.  The method by which the matrix is set up
  // with 1's may or may not have already loaded these.
  vals[0] = 1.0;
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    if (iterDI->gid !=-1)
    {
      cols[0] = iterDI->gid;
      bs1 = matPtr->putRow (iterDI->gid, 1, &vals[0], &cols[0]);
      bsuccess = bsuccess && bs1;
    }
  }

  // rows associated with the PDE mesh:
  if( useMatrixGIDFlag )
  {
    for (i=0;i<numMeshPoints;++i)
    {
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "\nmesh point i = " << i << std::endl;
      }
#endif
      for (j=0;j<numCol;++j)
      {
        cols[j] = -1;
        vals[j] = 0.0;
      }

      Vrow = Vrowarray[i];
      Nrow = Nrowarray[i];
      Prow = Prowarray[i];

#ifdef Xyce_NEW_BC
      if (boundarySten[i]) continue;
#else
      bool doneFlag = false;

      // Is this node a node with an explicit boundary condition?
      // If so, apply the BC.  Otherwise, just load the region operators.

      if (boundarySten[i]) // do BC load
      {
        int DIindex = labelDIMap[labelNameVector[i]];
        if (dIVec[DIindex].given)
        {
          count = 0;
          if (Vrow != -1)
          {
            ++count;
            vals[0] = 1.0;
            cols[0] = Vrow;

            if (dIVec[DIindex].gid != -1)
            {
              ++count;
              vals[1] = -scalingVars.rV0;
              cols[1] = dIVec[DIindex].gid;
            }
          }

          if (count > 0)
          {
            bs1 = matPtr->putRow (Vrow, count, &vals[0], &cols[0]);
            bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
            if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
            {
              Xyce::dout() << "BC  load:  Vrow = " << Vrow << std::endl;
              for (int eric=0;eric<count;++eric)
              {
                Xyce::dout() << "  cols["<<eric<<"] = " << cols[eric];
                Xyce::dout() << "  vals["<<eric<<"] = " << vals[eric];
                Xyce::dout() << std::endl;
              }
            }
#endif
          }

          if (Nrow != -1)
          {
            count = 1;
            vals[0] = 1.0;
            cols[0] = Nrow;
            bs1 = matPtr->putRow (Nrow, 1, &vals[0], &cols[0]);
            bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
            if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
            {
              Xyce::dout() << "BC  load:  Nrow = " << Nrow;
              Xyce::dout() << "  coef = " << vals[0] << std::endl;
            }
#endif
          }

          if (Prow != -1)
          {
            count = 1;
            vals[0] = 1.0;
            cols[0] = Prow;

            bs1 = matPtr->putRow (Prow, 1, &vals[0], &cols[0]);
            bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
            if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
            {
              Xyce::dout() << "BC  load:  Prow = " << Vrow;
              Xyce::dout() << "  coef = " << vals[0] << std::endl;
            }
#endif
          }

          doneFlag = true;

        } // inner if  (given if)
      } // outer if (name if)

      if (doneFlag) continue;
#endif // Xyce_NEW_BC

      // if interior point:

      // if load is not done yet, then do an interior point load:
#ifdef Xyce_OXIDE_ENABLED
      if (!allOxideFlag)
      {
        holeDens = getVoltDepHoleDens ( VminExp, VVec[i], NpMax);
        elecDens = getVoltDepElecDens ( VmaxExp, VVec[i], NnMax);
      }
#else
      holeDens = getVoltDepHoleDens ( VminExp, VVec[i], NpMax);
      elecDens = getVoltDepElecDens ( VmaxExp, VVec[i], NnMax);
#endif

      mNode * nodePtr = meshContainerPtr->getNode(i);
      nodeArea = nodePtr->area;

      // center point:
      coef = 0.0;
      count = 0;
      for (iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        int inodeB = nodePtr->edgeInfoVector[iNN].inode;
        ilen = nodePtr->edgeInfoVector[iNN].ilen;
        elen = nodePtr->edgeInfoVector[iNN].elen;

        coef +=  -ilen/elen;
      }
      coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;

      // now add the terms associated with the electron and hole densities:
#ifdef Xyce_OXIDE_ENABLED
      if (!allOxideFlag)
      {
        coef +=  pre*holeDens + pre*elecDens;
      }
#else
      coef += pre*holeDens + pre*elecDens;
#endif

      if (Vcolarray[i][0] != -1)
      {
        cols[count] = Vcolarray[i][0];
        vals[count] = coef;
        ++count;
      }
      else
      {
        Xyce::dout() << "    center point: i="<<i;
        Xyce::dout() << "    xVec[i] = " << xVec[i];
        Xyce::dout() << "    yVec[i] = " << yVec[i];
        Xyce::dout() << std::endl;
        Xyce::dout() << "    label = " << labelNameVector[i] << std::endl;
        Xyce::dout() << "    boundarySten = " << boundarySten[i] << std::endl;
        Xyce::dout() << "    Vcolarray[i][0] = ";
        Xyce::dout() << Vcolarray[i][0] << std::endl;
      }

      // neighbor points:
      coef = 0.0;
      for (iNN=0;iNN<nodePtr->cnode;++iNN)
      {
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "non-center point: i="<<i;
          Xyce::dout() << "  Vcolarray[i]["<<iNN+1<<"] = ";
          Xyce::dout() << Vcolarray[i][iNN+1] << std::endl;
        }
#endif
        if (Vcolarray[i][iNN+1] == -1) continue;

        int inodeB = nodePtr->edgeInfoVector[iNN].inode;
        ilen = nodePtr->edgeInfoVector[iNN].ilen;
        elen = nodePtr->edgeInfoVector[iNN].elen;

        coef   =  ilen/elen;
        coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;

        cols[count] = Vcolarray[i][iNN+1];
        vals[count] = coef;
        ++count;
      }

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "\n";
        Xyce::dout().setf(std::ios::left);
        Xyce::dout() << " POISSON LOAD: row = " << Vrow;
        Xyce::dout() << " count = " << count << "\n";
        Xyce::dout().setf(std::ios::scientific);
        for (j=0;j<count;++j)
        {
          Xyce::dout() << "  cols["<<j<<"] = " << cols[j];
          Xyce::dout() << "  vals["<<j<<"] = " << vals[j];
          Xyce::dout() << std::endl;
        }
      }
#endif

      if (Vrow != -1 && count > 0)
      {
        bs1 = matPtr->putRow (Vrow, count, &vals[0], &cols[0]);
        bsuccess = bsuccess && bs1;
      }
      else
      {
        Xyce::dout() << "OOOPS!" << std::endl;
        exit(0);
      }


      if (Nrow != -1)
      {
        vals[0] = 1.0;
        cols[0] = Nrow;
        bs1 = matPtr->putRow (Nrow, 1, &vals[0], &cols[0]);
        bsuccess = bsuccess && bs1;
      }

      if (Prow != -1)
      {
        vals[0] = 1.0;
        cols[0] = Prow;
        bs1 = matPtr->putRow (Prow, 1, &vals[0], &cols[0]);
        bsuccess = bsuccess && bs1;
      }

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << subsection_divider << "\n";
      }
#endif
    } // mesh point for loop
  }
  else // direct Matrix Access
  {
    for (i=0;i<numMeshPoints;++i)
    {
      Vrow = li_Vrowarray[i];
      Nrow = li_Nrowarray[i];
      Prow = li_Prowarray[i];

      std::vector<int> & Voff = li_VoffsetArray[i];
      std::vector<int> & Noff = li_NoffsetArray[i];
      std::vector<int> & Poff = li_PoffsetArray[i];

#ifdef Xyce_NEW_BC
      if (boundarySten[i]) continue;
#else
      bool doneFlag = false;

      // Is this node a node with an explicit boundary condition?
      // If so, apply the BC.  Otherwise, just load the region operators.

      if (boundarySten[i]) // do BC load
      {
        int DIindex = labelDIMap[labelNameVector[i]];
        if (dIVec[DIindex].given)
        {
          // poisson BC
          // Assuming the second offset is the (Vrow, ckt_gid) entry.
          // CHECK this!
          (*matPtr)[Vrow][Voff[0]] =  1.0;
          (*matPtr)[Vrow][Voff[1]] = -scalingVars.rV0;

          // electron equation is "off"
          (*matPtr)[Nrow][Noff[0]] =  1.0;

          // hole equation is "off"
          (*matPtr)[Prow][Poff[0]] =  1.0;

          doneFlag = true;
        } // inner if  (given if)
      } // outer if (name if)

      if (doneFlag) continue;
#endif // Xyce_NEW_BC

      // if interior point:

      // if load is not done yet, then do an interior point load:
#ifdef Xyce_OXIDE_ENABLED
      if (!allOxideFlag)
      {
        holeDens = getVoltDepHoleDens ( VminExp, VVec[i], NpMax);
        elecDens = getVoltDepElecDens ( VmaxExp, VVec[i], NnMax);
      }
#else
      holeDens = getVoltDepHoleDens ( VminExp, VVec[i], NpMax);
      elecDens = getVoltDepElecDens ( VmaxExp, VVec[i], NnMax);
#endif

      mNode * nodePtr = meshContainerPtr->getNode(i);
      nodeArea = nodePtr->area;

      // center point:
      coef = 0.0;
      count = 0;
      for (iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        int inodeB = nodePtr->edgeInfoVector[iNN].inode;
        ilen = nodePtr->edgeInfoVector[iNN].ilen;
        elen = nodePtr->edgeInfoVector[iNN].elen;

        coef +=  -ilen/elen;
      }
      coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;

      // now add the terms associated with the electron and hole densities:
#ifdef Xyce_OXIDE_ENABLED
      if (!allOxideFlag)
      {
        coef +=  pre*holeDens + pre*elecDens;
      }
#else
      coef +=  pre*holeDens + pre*elecDens;
#endif

      (*matPtr)[Vrow][Voff[0]] += coef;

      // neighbor points:
      coef = 0.0;
      for (iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        // If the neighbor is a boundary, skip.
        // Do I need to do this if statement? CHECK!
        if (Vcolarray[i][iNN+1] == -1) continue;

        int inodeB = nodePtr->edgeInfoVector[iNN].inode;
        ilen = nodePtr->edgeInfoVector[iNN].ilen;
        elen = nodePtr->edgeInfoVector[iNN].elen;

        coef   =  ilen/elen;
        coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;

        (*matPtr)[Vrow][Voff[iNN+1]] += coef;
      }

      // electron equation is "off"
      (*matPtr)[Nrow][Noff[0]] =  1.0;

      // hole equation is "off"
      (*matPtr)[Prow][Poff[0]] =  1.0;

    } // mesh point loop
  } // direct access

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatKCLDDForm
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/24/05
//-----------------------------------------------------------------------------
bool Instance::loadMatKCLDDForm (N_LAS_Matrix * matPtr)
{
  bool bsuccess = true;
  bool bs1 = true;
  int j;
  int count;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Starting Instance::loadJacKCLDDFormulation" << std::endl;
  }
#endif

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  if( useMatrixGIDFlag )
  {
    int numCol = cols.size();
    if (vals.size () < cols.size()) numCol = vals.size();

    // Initialize the cols and vals arrays:
    for (j=0;j<numCol;++j)
    {
      cols[j] = -1;
      vals[j] = 0.0;
    }

    for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
    {
      if (iterDI->gid ==-1) continue;
      count = 0;
      // dependence of electrode current on ckt node voltage.
      cols[count] = iterDI->gid;
      vals[count] = 0.0; // Setting this to zero is correct for
      // the "old" boundary conditions.  For the
      // the "new" boundary conditions, this will
      // eventually be a nonzero number.
      count=1;

#ifdef Xyce_NEW_BC
      vals[0]  = iterDI->dIdVckt;  // calculated in pdTerminalCurrents

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "\n";
        Xyce::dout().setf(std::ios::left);
        Xyce::dout() << "  KCL new BC LOAD: row = " << iterDI->gid;
        Xyce::dout() << " count = " << count;
        Xyce::dout() << " name = " << iterDI->eName << "\n";

        Xyce::dout().setf(std::ios::scientific);
        Xyce::dout() << "  cols[0] = " << cols[0];
        Xyce::dout() << "  vals[0] = " << vals[0];
        Xyce::dout() << std::endl;
      }
#endif
      bs1 = matPtr->sumIntoRow (iterDI->gid, count, &vals[0], &cols[0]);
      bsuccess = bsuccess && bs1;
#endif // Xyce_NEW_BC

      // Now do the dependence of the electrode current on
      // all the PDE vars.  This was calculated in pdTerminalCurrents.
      count  = iterDI->dIdXcols.size();

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "\n";
        Xyce::dout().setf(std::ios::left);
        Xyce::dout() << "  KCL LOAD: row = " << iterDI->gid;
        Xyce::dout() << " count = " << count;
        Xyce::dout() << " name = " << iterDI->eName << "\n";

        Xyce::dout().setf(std::ios::scientific);
        for (j=0;j<count;++j)
        {
          Xyce::dout() << "  cols["<<j<<"] = " << iterDI->dIdXcols[j];
          Xyce::dout() << "  vals["<<j<<"] = " << iterDI->dIdX[j];
          Xyce::dout() << std::endl;
        }
      }
#endif

      bs1 = matPtr->sumIntoRow
            (iterDI->gid, count, &(iterDI->dIdX[0]), &(iterDI->dIdXcols[0]));
      bsuccess = bsuccess && bs1;

    }  // end of DI loop
  }
  else // use direct matrix access
  {
    for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
    {
      // First, dependence of electrode current on ckt node voltage.
      double coef = 0.0;  // "old" BC value.
      int offset  = iterDI->lidOffset;

#ifdef Xyce_NEW_BC
      coef     = iterDI->dIdVckt;
#endif // Xyce_NEW_BC

      (*matPtr)[iterDI->lid][offset] += coef;

      // Now do the dependence of the electrode current on
      // all the PDE vars.  This was calculated in pdTerminalCurrents.
      int size  = iterDI->dIdXcols.size();
      for (int i=0;i<size;++i)
      {
        coef   = iterDI->dIdX[i];
        offset = iterDI->dIdXoffset[i];
        (*matPtr)[iterDI->lid][offset] += coef;
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "dIdX["<<i<<"] = " << iterDI->dIdX[i];
          Xyce::dout() << "  dIdXcols["<<i<<"] = " << iterDI->dIdXcols[i];
          Xyce::dout() << "  dIdXoffset["<<i<<"] = " << iterDI->dIdXoffset[i];
          Xyce::dout() << std::endl;
        }
#endif
      }
    }  // end of DI loop
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Done with Instance::loadJacKCLDDFormulation" << std::endl;
  }
#endif
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatCktTrivial
// Purpose       :
// Special Notes : if this PDE device is not coupled to the ckt, put
//                 some 1's on the diagonal just to be safe.  The method by
//                 which the matrix is set up with 1's may or may not have
//                 already loaded these.  This may be obsolete later, with
//                 new petra functionality.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/24/05
//-----------------------------------------------------------------------------
bool Instance::loadMatCktTrivial (N_LAS_Matrix * matPtr)
{
  bool bsuccess = true;
  bool bs1 = true;

  // note: fix the vals and cols arrays later!
  int numCol = cols.size();
  if (vals.size () < cols.size()) numCol = vals.size();

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  if( useMatrixGIDFlag )
  {
    vals[0] = 1.0;
    for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
    {
      if (iterDI->gid !=-1)
      {
        cols[0] = iterDI->gid;
        bs1 = matPtr->putRow (iterDI->gid, 1, &vals[0], &cols[0]);
        bsuccess = bsuccess && bs1;
      }
    }
  }
  else // using direct matrix access
  {
    for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
    {
      (*matPtr)[iterDI->lid][iterDI->lidOffset] = 1.0;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatDDForm
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/24/05
//-----------------------------------------------------------------------------
bool Instance::loadMatDDForm
(double dndtScalar, N_LAS_Matrix * matPtr)
{
  bool bsuccess = true;
  bool bs1 = true;

  int Vrow, Nrow, Prow;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\n"<<section_divider << "\n";
    Xyce::dout() << "Instance::loadMatDDForm" << "\n";
    Xyce::dout() << "  name = " << getName() <<"\n";
    Xyce::dout() << "\n";
  }
#endif

  int i,j;
  int iNN;
  int count  = 0;
  int cnt2   = 0;
  int iedge  = 0;
  int inodeB = 0;
  bool bmatch;

  double coef;
  double ilen, elen, nodeArea;

  // note: fix the vals and cols arrays later!
  int numCol = cols.size();
  if (vals.size () < cols.size()) numCol = vals.size();

#ifdef Xyce_DEBUG_DEVICE
#endif

  // obtain partial time derivatives, scale by scalingVars.t0:
  // dndtScalar is either 1.0 (old DAE) or 0.0 (new DAE).
  double dDNDTdn;
  double dDPDTdp;

  if (!(getSolverState().dcopFlag))
  {
    dDNDTdn =  getSolverState().pdt*scalingVars.t0*dndtScalar;
    dDPDTdp =  getSolverState().pdt*scalingVars.t0*dndtScalar;
  }
  else
  {
    dDNDTdn = 0.0;
    dDPDTdp = 0.0;
  }

  // Initialize the cols and vals arrays:
  for (j=0;j<numCol;++j) { cols[j] = -1; vals[j] = 0.0; }

  // load the rows associated with the PDE mesh:
  if( useMatrixGIDFlag )
  {
    for (i=0;i<numMeshPoints;++i)
    {
      Vrow = Vrowarray[i];
      Nrow = Nrowarray[i];
      Prow = Prowarray[i];

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << subsection_divider << "\n";
        Xyce::dout() << "mesh point i = " << i;
        Xyce::dout() << "  Vrow = " << Vrow;
        Xyce::dout() << "  Nrow = " << Nrow;
        Xyce::dout() << "  Prow = " << Prow;
        Xyce::dout() << "\n";
      }
#endif

      bool doneFlagV = false;
      bool doneFlagN = false;
      bool doneFlagP = false;
      bool doneFlag  = false;
#ifdef Xyce_NEW_BC
      if (boundarySten[i]) continue;
#else
      // Is this node a node with an explicit boundary condition?
      // If so, apply the BC.  Otherwise, just load the region operators.

      if (boundarySten[i]) //do BC load
      {
        int DIindex = labelDIMap[labelNameVector[i]];

        if (dIVec[DIindex].given)
        {
          if (dIVec[DIindex].neumannBCFlagV==false)
          {
#ifdef Xyce_DEBUG_DEVICE
            if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
            {
              Xyce::dout() << "BC  load:  Vrow = " << Vrow << std::endl;
            }
#endif
            count = 0;
            if (Vrow != -1)
            {
              ++count;
              vals[0] = 1.0;
              cols[0] = Vrow;

              if (dIVec[DIindex].gid != -1)
              {
                ++count;
                vals[1] = -scalingVars.rV0;
                cols[1] = dIVec[DIindex].gid;
              }
            }

            if (count > 0)
            {
              bs1 = matPtr->putRow (Vrow, count, &vals[0], &cols[0]);
              bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
              if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
              {
                Xyce::dout() << "BC  load:  Vrow = " << Vrow << std::endl;
                for (int eric=0;eric<count;++eric)
                {
                  Xyce::dout() << "  cols["<<eric<<"] = " << cols[eric];
                  Xyce::dout() << "  vals["<<eric<<"] = " << vals[eric];
                  Xyce::dout() << std::endl;
                }
              }
#endif

            } // Vrow if statement.
            doneFlagV = true;
          } // neumann flag if statement.

          if (dIVec[DIindex].neumannBCFlagN==false)
          {
            if (Nrow != -1)
            {
              count = 1;
              vals[0] = 1.0;
              cols[0] = Nrow;
              bs1 = matPtr->putRow (Nrow, 1, &vals[0], &cols[0]);
              bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
              if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
              {
                Xyce::dout() << "BC  load:  Nrow = " << Nrow << std::endl;
                for (int eric=0;eric<count;++eric)
                {
                  Xyce::dout() << "  cols["<<eric<<"] = " << cols[eric];
                  Xyce::dout() << "  vals["<<eric<<"] = " << vals[eric];
                  Xyce::dout() << std::endl;
                }
              }
#endif
            } // Nrow if statement.
            doneFlagN = true;
          } // neumann flag if statement.

          if (dIVec[DIindex].neumannBCFlagP==false)
          {
            if (Prow != -1)
            {
              count = 1;
              vals[0] = 1.0;
              cols[0] = Prow;
              bs1 = matPtr->putRow (Prow, 1, &vals[0], &cols[0]);
              bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
              if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
              {
                Xyce::dout() << "BC  load:  Prow = " << Prow << std::endl;
                for (int eric=0;eric<count;++eric)
                {
                  Xyce::dout() << "  cols["<<eric<<"] = " << cols[eric];
                  Xyce::dout() << "  vals["<<eric<<"] = " << vals[eric];
                  Xyce::dout() << std::endl;
                }
              }
#endif
            } // Prow if statement.
            doneFlagP = true;
          } // neumann flag if statement.

          doneFlag = (doneFlagV && doneFlagN && doneFlagP);
        } // inner if
      } // outer if

      if (doneFlag) continue;

#endif // Xyce_NEW_BC

      // if load is not done yet, then do an interior point load:
      std::string semi(bulkMaterial);
      // Poisson's equation:
      mNode * nodePtr = meshContainerPtr->getNode(i);
      nodeArea = nodePtr->area;

      // center point, with respect to V:
      coef = 0.0;
      count = 0;
      for (j=0;j<numCol;++j) { cols[j] = -1; vals[j] = 0.0; }

      iNN=0;

      if (!doneFlagV)
      {
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          int inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen = nodePtr->edgeInfoVector[iNN].ilen;
          elen = nodePtr->edgeInfoVector[iNN].elen;

          coef +=  -ilen/elen;
        }
        coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;

        if (Vcolarray[i][0] != -1)
        {
          cols[count] = Vcolarray[i][0];
          vals[count] = coef;
          ++count;
        }
#ifdef Xyce_DEBUG_DEVICE
        else
        {
          Xyce::dout() << "    center point: i="<<i;
          Xyce::dout() << "  Vcolarray[i][0] = " << Vcolarray[i][0] << std::endl;
        }
#endif

        // neighbor points, with respect to V:
        coef = 0.0;
        iNN=0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
#ifdef Xyce_DEBUG_DEVICE
          if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
          {
            Xyce::dout() << "non-center point: i="<<i;
            Xyce::dout() << "  Vcolarray[i]["<<iNN+1<<"] = ";
            Xyce::dout() << Vcolarray[i][iNN+1] << std::endl;
          }
#endif
          int tmpCol = Vcolarray[i][iNN+1];

          if (tmpCol == -1) continue;

          int inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen = nodePtr->edgeInfoVector[iNN].ilen;
          elen = nodePtr->edgeInfoVector[iNN].elen;

          coef   =  ilen/elen;
          coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;
#ifdef Xyce_NEW_BC
          coef  *= ((boundarySten[inodeB]==1)?(scalingVars.rV0):1.0);

          // check if node == any previous nodes in the cols array
          if (boundarySten[inodeB]==1)
          {
            bmatch = false;
            for (cnt2=0;cnt2<count;++cnt2)
            {
              if (cols[cnt2] == tmpCol)
              { bmatch = true; break; }
            }
            if (!bmatch) { cnt2 = count;  ++count;}
          }
          else
          {
            cnt2 = count;
            ++count;
          }
#else
          cnt2 = count;
          ++count;
#endif // Xyce_NEW_BC
          cols[cnt2] = tmpCol;
          vals[cnt2] += coef;
        }

        // center point, with repsect to electron density:
        if (Ncolarray[i][0] != -1)
        {
          cols[count] = Ncolarray[i][0];
          vals[count] = +1.0;
          ++count;
        }

        // center point, with respect to hole density:
        if (Pcolarray[i][0] != -1)
        {
          cols[count] = Pcolarray[i][0];
          vals[count] = -1.0;
          ++count;
        }

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "\n";
          Xyce::dout().setf(std::ios::left);
          Xyce::dout() << " POISSON LOAD: row = " << Vrow;
          Xyce::dout() << "  count = " << count << "\n";
          Xyce::dout().setf(std::ios::scientific);
          for (j=0;j<count;++j)
          {
            Xyce::dout() << "  cols["<<j<<"] = " << cols[j];
            Xyce::dout() << "  vals["<<j<<"] = " << vals[j];
            Xyce::dout() << std::endl;
          }
        }
#endif

        if (Vrow != -1 && count > 0)
        {
          bs1 = matPtr->sumIntoRow (Vrow, count, &vals[0], &cols[0]);
          bsuccess = bsuccess && bs1;
        }
        else
        {
          Xyce::dout() << "OOOPS 1!" << std::endl;
          exit(0);
        }
        doneFlagV = true;
      } // doneFlagV if statement.


      // electron continuity row:
      if (!doneFlagN)
      {
        for (j=0;j<numCol;++j) { cols[j] = -1; vals[j] = 0.0; }

        count = 0;
        coef = 0.0;

        // center point first: (with respect to nn)
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdn = 0.0;
          if (i<inodeB)
          {
            dJdn = dJndn1Vec[iedge];
          }
          else
          {
            dJdn = dJndn2Vec[iedge];
          }

          coef += ((i<inodeB)?1.0:-1.0) * dJdn * ilen;
        }

        coef /= nodeArea;
        coef += - dRdnVec[i] - dDNDTdn;

        coef += pdElecPenalty[i];

        vals[count] = coef;
        cols[count] = Ncolarray[i][0];
        ++count;

        // neighbor points, with respect to nn:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          if (Ncolarray[i][iNN+1] == -1) continue;

          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdn = 0.0;
          if (i>inodeB)
          {
            dJdn = dJndn1Vec[iedge];
          }
          else
          {
            dJdn = dJndn2Vec[iedge];
          }

          coef = ((i<inodeB)?1.0:-1.0) * dJdn * ilen/nodeArea;

          vals[count] = coef;
          cols[count] = Ncolarray[i][iNN+1];

          ++count;
        }

        // center point, with respect to V:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdV = 0.0;
          if (i<inodeB)
          {
            dJdV = dJndV1Vec[iedge];
          }
          else
          {
            dJdV = dJndV2Vec[iedge];
          }

          coef += ((i<inodeB)?1.0:-1.0) * dJdV * ilen;
        }
        coef /= nodeArea;

        if (Vcolarray[i][0] != -1)
        {
          vals[count] = coef;
          cols[count] = Vcolarray[i][0];
          ++count;
        }

        // neighbor points, with respect to V:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          int tmpCol = Vcolarray[i][iNN+1];
          if (tmpCol == -1) continue;

          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdV = 0.0;

          if (i>inodeB)
          {
            dJdV = dJndV1Vec[iedge];
          }
          else
          {
            dJdV = dJndV2Vec[iedge];
          }

          coef = ((i<inodeB)?1.0:-1.0) * dJdV * ilen/nodeArea;
#ifdef Xyce_NEW_BC
          coef  *= ((boundarySten[inodeB]==1)?(scalingVars.rV0):1.0);

          // check if node == any previous nodes in the cols array
          if (boundarySten[inodeB]==1)
          {
            bmatch = false;
            for (cnt2=0;cnt2<count;++cnt2)
            {
              if (cols[cnt2] == tmpCol)
              { bmatch = true; break; }
            }
            if (!bmatch) { cnt2 = count;  ++count;}
          }
          else
          {
            cnt2 = count;
            ++count;
          }
#else
          cnt2 = count;
          ++count;
#endif // Xyce_NEW_BC
          vals[cnt2] += coef;
          cols[cnt2] = tmpCol;
        }

        // center point, with respect to np:
        if (Pcolarray[i][0] != -1)
        {
          vals[count] = -dRdpVec[i];
          cols[count] = Pcolarray[i][0];
          ++count;
        }

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "\n";
          Xyce::dout().setf(std::ios::left);
          Xyce::dout() << "ELECTRON LOAD:  row = " << Nrow;
          Xyce::dout() << "  count = " << count<< "\n";
          Xyce::dout().setf(std::ios::scientific);
          for (j=0;j<count;++j)
          {
            Xyce::dout() << "  cols["<<j<<"] = " << cols[j];
            Xyce::dout() << "  vals["<<j<<"] = " << vals[j];
            Xyce::dout() << std::endl;
          }
        }
#endif

        if (Nrow != -1 && count > 0)
        {
          bs1 = matPtr->sumIntoRow (Nrow,count,&vals[0],&cols[0]);
          bsuccess = bsuccess && bs1;
        }
        else
        {
          Xyce::dout() << "OOOPS 2!" << std::endl;
          exit(0);
        }

        doneFlagN = true;

      } // doneFlagN if statement.

      // hole continuity row:
      if (!doneFlagP)
      {
        for (j=0;j<numCol;++j) { cols[j] = -1; vals[j] = 0.0; }

        count = 0;
        coef  = 0.0;

        // center point first: (with respect to np)
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdp = 0.0;
          if (i<inodeB)
          {
            dJdp = dJpdn1Vec[iedge];
          }
          else
          {
            dJdp = dJpdn2Vec[iedge];
          }

          coef += - ((i<inodeB)?1.0:-1.0) * dJdp * ilen;
        }

        coef /= nodeArea;
        coef += - dRdpVec[i] - dDPDTdp;

        coef += pdHolePenalty[i];

        vals[count] = coef;
        cols[count] = Pcolarray[i][0];
        ++count;

        // neighbor points, with respect to np:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          if (Pcolarray[i][iNN+1] == -1) continue;

          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdp = 0.0;
          if (i>inodeB)
          {
            dJdp = dJpdn1Vec[iedge];
          }
          else
          {
            dJdp = dJpdn2Vec[iedge];
          }

          coef =-((i<inodeB)?1.0:-1.0) * dJdp * ilen/nodeArea;

          vals[count] = coef;
          cols[count] = Pcolarray[i][iNN+1];
          ++count;
        }

        // center point, with respect to V:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdV = 0.0;
          if (i<inodeB)
          {
            dJdV = dJpdV1Vec[iedge];
          }
          else
          {
            dJdV = dJpdV2Vec[iedge];
          }

          coef += -((i<inodeB)?1.0:-1.0) * dJdV * ilen;
        }
        coef /= nodeArea;

        if (Vcolarray[i][0] != -1)
        {
          vals[count] = coef;
          cols[count] = Vcolarray[i][0];
          ++count;
        }

        // neighbor points, with respect to V:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          int tmpCol = Vcolarray[i][iNN+1];
          if (tmpCol == -1) continue;

          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdV = 0.0;

          if (i>inodeB)
          {
            dJdV = dJpdV1Vec[iedge];
          }
          else
          {
            dJdV = dJpdV2Vec[iedge];
          }

          coef =-((i<inodeB)?1.0:-1.0) * dJdV * ilen/nodeArea;
#ifdef Xyce_NEW_BC
          coef  *= ((boundarySten[inodeB]==1)?(scalingVars.rV0):1.0);

          // check if node == any previous nodes in the cols array
          if (boundarySten[inodeB]==1)
          {
            bmatch = false;
            for (cnt2=0;cnt2<count;++cnt2)
            {
              if (cols[cnt2] == tmpCol)
              { bmatch = true; break; }
            }
            if (!bmatch) { cnt2 = count;  ++count;}
          }
          else
          {
            cnt2 = count;
            ++count;
          }
#else
          cnt2 = count;
          ++count;
#endif // Xyce_NEW_BC
          vals[cnt2] += coef;
          cols[cnt2] = tmpCol;
        }

        // center point, with respect to nn:
        if (Ncolarray[i][0] != -1)
        {
          vals[count] = -dRdnVec[i];
          cols[count] = Ncolarray[i][0];
          ++count;
        }

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "\n";
          Xyce::dout().setf(std::ios::left);
          Xyce::dout() << "    HOLE LOAD:  row = " << Prow;
          Xyce::dout() << "  count = " << count<< "\n";
          Xyce::dout().setf(std::ios::scientific);
          for (j=0;j<count;++j)
          {
            Xyce::dout() << "  cols["<<j<<"] = " << cols[j];
            Xyce::dout() << "  vals["<<j<<"] = " << vals[j];
            Xyce::dout() << std::endl;
          }
        }
#endif

        if (Prow != -1 && count > 0)
        {
          bs1 = matPtr->sumIntoRow (Prow,count,&vals[0],&cols[0]);
          bsuccess = bsuccess && bs1;
        }
        else
        {
          Xyce::dout() << "OOOPS 3!" << std::endl;
          exit(0);
        }

        doneFlagP = true;

      } //  doneFlagP if statement.

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << subsection_divider << "\n";
      }
#endif
    } // mesh loop.
  }
  else // direct matrix access:
  {
    for (i=0;i<numMeshPoints;++i)
    {
      int ioffset;

      Vrow = li_Vrowarray[i];
      Nrow = li_Nrowarray[i];
      Prow = li_Prowarray[i];

      std::vector<int> & Voff = li_VoffsetArray[i];
      std::vector<int> & Noff = li_NoffsetArray[i];
      std::vector<int> & Poff = li_PoffsetArray[i];

      bool doneFlagV = false;
      bool doneFlagN = false;
      bool doneFlagP = false;
      bool doneFlag  = false;

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "mesh point i = " << i << std::endl;
      }
#endif

#ifdef Xyce_NEW_BC
      if (boundarySten[i]) continue;
#else
      // Is this node a node with an explicit boundary condition?
      // If so, apply the BC.  Otherwise, just load the region operators.

      if (boundarySten[i]) //do BC load
      {
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "  Doing BC load" << std::endl;
        }
#endif
        int DIindex = labelDIMap[labelNameVector[i]];

        if (dIVec[DIindex].given)
        {
          if (dIVec[DIindex].neumannBCFlagV==false)
          {
            // This assumes the first offset is the (Vrow, Vrow) entry,
            // and the second is the (Vrow, gid) entry. CHECK!
            // This is probably an incorrect assumption...
            (*matPtr)[Vrow][Voff[0]] =  1.0;
            (*matPtr)[Vrow][Voff[1]] = -scalingVars.rV0;

            doneFlagV = true;
          } // neumann flag if statement.

          if (dIVec[DIindex].neumannBCFlagN==false)
          {
            // This assumes the first offset is the (Nrow, Nrow) entry.
            // CHECK!
            (*matPtr)[Nrow][Noff[0]] =  1.0;

            doneFlagN = true;
          } // neumann flag if statement.

          if (dIVec[DIindex].neumannBCFlagP==false)
          {
            // This assumes the first offset is the (Prow, Prow) entry.
            // CHECK!
            (*matPtr)[Prow][Poff[0]] =  1.0;

            doneFlagP = true;
          } // neumann flag if statement.

          doneFlag = (doneFlagV && doneFlagN && doneFlagP);
        } // inner if
      } // outer if

      if (doneFlag) continue;

#endif // Xyce_NEW_BC

      // if load is not done yet, then do an interior point load:

      // Poisson's equation:
      mNode * nodePtr = meshContainerPtr->getNode(i);
      nodeArea = nodePtr->area;

      // center point, with respect to V:
      coef = 0.0;
      iNN=0;

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "  Doing Poisson load" << std::endl;
      }
#endif
      std::string semi(bulkMaterial);
      if (!doneFlagV)
      {
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          int inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen = nodePtr->edgeInfoVector[iNN].ilen;
          elen = nodePtr->edgeInfoVector[iNN].elen;

          coef +=  -ilen/elen;
        }
        coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;

        (*matPtr)[Vrow][Voff[0]] += coef;

        // neighbor points, with respect to V:
        coef = 0.0; iNN=0; ioffset=1;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          if (Vcolarray[i][iNN+1] == -1) continue;

          int inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen = nodePtr->edgeInfoVector[iNN].ilen;
          elen = nodePtr->edgeInfoVector[iNN].elen;

          coef   =  ilen/elen;
          coef  *= -scalingVars.L0 * matSupport.getRelPerm(semi)/nodeArea;
#ifdef Xyce_NEW_BC // again CHECK the offset value.
          coef  *= ((boundarySten[inodeB]==1)?(scalingVars.rV0):1.0);
#endif // Xyce_NEW_BC
          (*matPtr)[Vrow][Voff[ioffset]] += coef;
          ++ioffset;
        }

        // center point, with repsect to electron density:
        (*matPtr)[Vrow][Voff[ioffset]] +=  1.0;
        ++ioffset;

        // center point, with respect to hole density:
        (*matPtr)[Vrow][Voff[ioffset]] += -1.0;
        ++ioffset;

        doneFlagV = true;
      } // doneFlagV if statement.

      // electron continuity row:
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "  Doing electron load" << std::endl;
      }
#endif

      if (!doneFlagN)
      {
        coef = 0.0; ioffset=0;
        // center point first: (with respect to nn)
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdn = 0.0;
          if (i<inodeB)
          {
            dJdn = dJndn1Vec[iedge];
          }
          else
          {
            dJdn = dJndn2Vec[iedge];
          }

          coef += ((i<inodeB)?1.0:-1.0) * dJdn * ilen;
        }

        coef /= nodeArea;
        coef += - dRdnVec[i] - dDNDTdn;
        coef += pdElecPenalty[i];

        (*matPtr)[Nrow][Noff[0]] += coef;

        // neighbor points, with respect to nn:
        coef = 0.0; ioffset=1;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          if (Ncolarray[i][iNN+1] == -1) continue;

          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdn = 0.0;
          if (i>inodeB)
          {
            dJdn = dJndn1Vec[iedge];
          }
          else
          {
            dJdn = dJndn2Vec[iedge];
          }

          coef = ((i<inodeB)?1.0:-1.0) * dJdn * ilen/nodeArea;

          (*matPtr)[Nrow][Noff[ioffset]] += coef;
          ++ioffset;
        }

        // center point, with respect to V:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdV = 0.0;
          if (i<inodeB)
          {
            dJdV = dJndV1Vec[iedge];
          }
          else
          {
            dJdV = dJndV2Vec[iedge];
          }

          coef += ((i<inodeB)?1.0:-1.0) * dJdV * ilen;
        }
        coef /= nodeArea;

        (*matPtr)[Nrow][Noff[ioffset]] += coef;
        ++ioffset;

        // neighbor points, with respect to V:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          if (Vcolarray[i][iNN+1] == -1) continue;

          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdV = 0.0;

          if (i>inodeB)
          {
            dJdV = dJndV1Vec[iedge];
          }
          else
          {
            dJdV = dJndV2Vec[iedge];
          }

          coef = ((i<inodeB)?1.0:-1.0) * dJdV * ilen/nodeArea;
#ifdef Xyce_NEW_BC
          coef  *= ((boundarySten[inodeB]==1)?(scalingVars.rV0):1.0);
#endif // Xyce_NEW_BC
          (*matPtr)[Nrow][Noff[ioffset]] += coef;
          ++ioffset;
        }

        // center point, with respect to np:
        (*matPtr)[Nrow][Noff[ioffset]] += -dRdpVec[i];
        ++ioffset;

        doneFlagN = true;

      } // doneFlagN if statement.

      // hole continuity row:
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << "  Doing hole load" << std::endl;
      }
#endif
      if (!doneFlagP)
      {
        coef  = 0.0; ioffset=0;
        // center point first: (with respect to np)
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdp = 0.0;
          if (i<inodeB)
          {
            dJdp = dJpdn1Vec[iedge];
          }
          else
          {
            dJdp = dJpdn2Vec[iedge];
          }

          coef += - ((i<inodeB)?1.0:-1.0) * dJdp * ilen;
        }

        coef /= nodeArea;
        coef += - dRdpVec[i] - dDPDTdp;
        coef += pdHolePenalty[i];

        (*matPtr)[Prow][Poff[0]] += coef;

        // neighbor points, with respect to np:
        coef = 0.0; ioffset=1;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          if (Pcolarray[i][iNN+1] == -1) continue;

          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdp = 0.0;
          if (i>inodeB)
          {
            dJdp = dJpdn1Vec[iedge];
          }
          else
          {
            dJdp = dJpdn2Vec[iedge];
          }

          coef =-((i<inodeB)?1.0:-1.0) * dJdp * ilen/nodeArea;

          (*matPtr)[Prow][Poff[ioffset]] += coef;
          ++ioffset;
        }

        // center point, with respect to V:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdV = 0.0;
          if (i<inodeB)
          {
            dJdV = dJpdV1Vec[iedge];
          }
          else
          {
            dJdV = dJpdV2Vec[iedge];
          }

          coef += -((i<inodeB)?1.0:-1.0) * dJdV * ilen;
        }
        coef /= nodeArea;

        (*matPtr)[Prow][Poff[ioffset]] += coef;
        ++ioffset;

        // neighbor points, with respect to V:
        coef = 0.0;
        for (iNN=0;iNN<nodePtr->cnode;++iNN)
        {
          if (Vcolarray[i][iNN+1] == -1) continue;

          iedge  = nodePtr->edgeInfoVector[iNN].iedge;
          inodeB = nodePtr->edgeInfoVector[iNN].inode;
          ilen   = nodePtr->edgeInfoVector[iNN].ilen;
          elen   = nodePtr->edgeInfoVector[iNN].elen;

          double dJdV = 0.0;

          if (i>inodeB)
          {
            dJdV = dJpdV1Vec[iedge];
          }
          else
          {
            dJdV = dJpdV2Vec[iedge];
          }

          coef =-((i<inodeB)?1.0:-1.0) * dJdV * ilen/nodeArea;
#ifdef Xyce_NEW_BC // CHECK!
          coef  *= ((boundarySten[inodeB]==1)?(scalingVars.rV0):1.0);
#endif // Xyce_NEW_BC
          (*matPtr)[Prow][Poff[ioffset]] += coef;
          ++ioffset;
        }

        // center point, with respect to nn:
        (*matPtr)[Prow][Poff[ioffset]] += -dRdnVec[i];
        ++ioffset;

        doneFlagP = true;

      } //  doneFlagP if statement.
    } // mesh loop.
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcLifetimes
// Purpose       : This function calculates the electron and hole lifetimes
//                 and places them into the tn and tp arrays.
// Special Notes : This function assumes scaling off.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcLifetimes ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcLifetimes" << "\n";
  }
#endif

  for (i=0;i<numMeshPoints;++i)
  {
    tnVec[i] = matSupport.calcLt (false, fabs(CVec[i]));
    tpVec[i] = matSupport.calcLt (true , fabs(CVec[i]));
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    for (i=0;i<numMeshPoints;++i)
    {
      Xyce::dout() << "  tnVec["<<i<<"] = " <<tnVec[i];
      Xyce::dout() << "  tpVec["<<i<<"] = " <<tpVec[i];
      Xyce::dout() << "\n";
    }
  }

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcMobilities
// Purpose       : This function calculates the electron and hole mobilities
//                 and places them into the un and up arrays.
// Special Notes : This function assumes scaling off.
//
//  ERK: 10/25/2012:  The mobility functions have  been updated to use
//  Sacado, and handle field-dependent mobilities.  However, that usage
//  has not been extended to the 2D devices yet.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcMobilities ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcMobilities" << "\n";
  }
#endif

  MobInfo<double> ci;
  ci.mobModelName = mobModelName;
  ci.materialName = bulkMaterial;
  for (i=0;i<numMeshPoints;++i)
  {
    // ci.materialName = labelNameVector[i];
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

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Nodal Mobilities:" << std::endl;
    for (i=0;i<numMeshPoints;++i)
    {
      Xyce::dout() << "  unVec["<<i<<"] = " << unVec[i];
      Xyce::dout() << "  upVec["<<i<<"] = " << upVec[i];
      Xyce::dout() << "\n";
    }
    Xyce::dout() << std::endl;
  }
#endif

  for (i=0;i< meshContainerPtr->getNumEdges ();++i)
  {
    mEdge * edgePtr = meshContainerPtr->getEdge(i);

    int inodeA = edgePtr->inodeA;
    int inodeB = edgePtr->inodeB;

    // option 1
    if (ci.mobModelName != "carr")
    {
      unE_Vec[i] = (unVec[inodeA]+unVec[inodeB])/2.0;
      upE_Vec[i] = (upVec[inodeA]+upVec[inodeB])/2.0;
    }
    // option 2
    else if(ci.mobModelName == "carr")
    {
      ci.N = (fabs(CVec[inodeA])+fabs(CVec[inodeB]))*0.5;
      ci.N *= ((variablesScaled)?scalingVars.C0:1.0);

      if (ci.N == 0.0) ci.N = 1.0;

      // for the carrier densities, do a "product average".
      ci.n = pow((fabs(nnVec[inodeA])*fabs(nnVec[inodeB])),0.5)*
             ((variablesScaled)?scalingVars.C0:1.0);

      ci.p = pow((fabs(npVec[inodeA])*fabs(npVec[inodeB])),0.5)*
             ((variablesScaled)?scalingVars.C0:1.0);

#ifdef Xyce_DEBUG_DEVICE
      if (ci.n != 0.0 && !(ci.n > 0.0) && !(ci.n < 0.0))
      {
        Xyce::dout() << "ci.n is nan" << std::endl;
        Xyce::dout() << "nnVec[A] = " << nnVec[inodeA];
        Xyce::dout() << "nnVec[B] = " << nnVec[inodeB];
        Xyce::dout() << std::endl;
        exit(0);
      }
#endif
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

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Edge Mobilities" << std::endl;
    for (i=0;i<meshContainerPtr->getNumEdges();++i)
    {
      Xyce::dout() << "  unE_Vec["<<i<<"] = " << unE_Vec[i];
      Xyce::dout() << "  upE_Vec["<<i<<"] = " << upE_Vec[i];
      Xyce::dout() << "\n";
    }
  }

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
//
// Purpose       : Re-sets neccessary things for a change in device
//                 temperature.
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
//        other ???
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature(const double & temp_tmp)
{
  bool bsuccess = true;
  bool bs1 = true;

  Temp = temp_tmp;

  // first un-scale everything, if neccessary:
  if (variablesScaled)
  {
    bs1 = unScaleVariables ();   bsuccess = bsuccess && bs1;
  }

  bs1 = setupMiscConstants ();   bsuccess = bsuccess && bs1;
  bs1 = setupScalingVars ();     bsuccess = bsuccess && bs1;

  bs1 = calcDensityBCs       (); bsuccess = bsuccess && bs1;
  bs1 = calcVequBCs          (); bsuccess = bsuccess && bs1;
  bs1 = calcMobilities       (); bsuccess = bsuccess && bs1;

  if (!variablesScaled)
  {
    bs1 = scaleVariables ();
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
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcVoltDepDensities ()
{
  bool bsuccess = true;

  int i;
  Ut = Vt/scalingVars.V0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcVoltDepDensities\n";
  }
#endif

  for (i=0;i<numMeshPoints;++i)
  {
    npVec[i] = getVoltDepHoleDens(VminExp, VVec[i], NpMax);
    nnVec[i] = getVoltDepElecDens(VmaxExp, VVec[i], NnMax);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcDopingProfile
// Purpose       : This function sets up the initial doping profile.  It should
//                 probably only be called once.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcDopingProfile ()
{
  bool bsuccess = true;
  int i;

  // set up the initial doping array:

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcDopingProfile\n";
  }
#endif

  // zero out CVec, first:
  CVec.resize (numMeshPoints,0.0);
  CdonorVec.resize (numMeshPoints,0.0);
  CacceptorVec.resize (numMeshPoints,0.0);
  //for (int itmp=0;itmp<numMeshPoints;++itmp) CVec[itmp] = 0.0;

#ifdef Xyce_OXIDE_ENABLED
  if (allOxideFlag) return bsuccess;
#endif

  // The doping profile may be contained in the mesh file.
  // or it may need to be specified here.  First, check the
  // mesh container to see if it has this info.

  // use the one owned by the mesh object.
  if (meshContainerPtr->dopingVectorExist ())
  {
    meshContainerPtr->getDopingVector(CVec);

    // Obtain Na and Nd:
    // Nd = number of donors    (maximum + doping)
    // Na = number of acceptors (minimum + doping)
    double NaTmp = 1.0e+99;
    double NdTmp =-1.0e+99;

    for (i=0;i<numMeshPoints;++i)
    {
      if(NaTmp > CVec[i]) NaTmp = CVec[i];
      if(NdTmp < CVec[i]) NdTmp = CVec[i];
    }

    Na = fabs(NaTmp);
    Nd = fabs(NdTmp);
  }
  else  // allocate, set it up here.
  {
    // Check which style of doping specification to use.  If the dopeInfoMap
    // is empty, then assume the old method.  If not, use the dopeInfoMap.
    if (dopeInfoMap.empty ())
    {
      // Two electrodes - do a diode sim.  This is the same code as
      // in the DiodePDE device, mostly, except that the junction is
      // perpendicular to the y-direction, rather than the x-direction.
      if (numElectrodes == 2)
      {
        double midpoint;
        midpoint = deviceWidth/2.0;
        XL       = midpoint - WJ/2.0;
        XR       = midpoint + WJ/2.0;

        // first setup, check the graded junction parameters:
        if (gradedJunctionFlag)
        {
          // if junction width was not specified, set to 1/10 the diode width.
          if (!given("WJ")) WJ = 0.1 * deviceWidth;

          midpoint = deviceWidth/2.0;
          XL = midpoint - WJ/2.0;
          XR = midpoint + WJ/2.0;
        }

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "deviceWidth  = " << deviceWidth  << std::endl;
          Xyce::dout() << "midpoint     = " << midpoint << std::endl;
          Xyce::dout() << "XL           = " << XL << std::endl;
          Xyce::dout() << "XR           = " << XR << std::endl;
          Xyce::dout() << "WJ           = " << WJ << std::endl;
        }
#endif

        for (i=0;i<numMeshPoints;++i)
        {
          double yloc = yVec[i];

          if (gradedJunctionFlag)
          {
            if (yloc <= XL)              CVec[i] = +Nd;
            else if (yloc>XL && yloc<XR) CVec[i] = Nd-(Na+Nd)/(XR-XL)*(yloc-XL);
            else                         CVec[i] = -Na;
          }
          else
          {
            if (yloc < deviceWidth /2.0)  CVec[i] = Nd;
            else                          CVec[i] = -Na;
          }
        }  // i loop
      }
      else if (numElectrodes == 3)
      {
        // Assuming a BJT.
        //
        //   emitter         base              collector
        //   A-----B--------C-----D-------------E-----F---------G
        //   |     |  /     |     |    |        |     |         |
        //   |-----|-       |     |   /         |     |         |
        //   |     |        |     |  /          |     |         |
        //   |-----|--------|-----|-            |     |         |
        //   |     |        |     |             |     |         |
        //   |     |        |     |             |     |         |
        //   |     |        |     |             |     |         |
        //   H-----I--------J-----K-------------L-----M---------N
        //

        // hardwired doping constants for now.  Fix this later.
        // Assuming a PNP, not NPN.

        double Ne;  // emitter doping
        double Nb;  // base doping
        double Nc;  // collector doping
        //int nx = numMeshPointsX;

        // Wad - total width of the emitter, base,
        // and the space in between.
        double Wad = 0.4 * deviceLength;

        // We - width of the emittter, by itself.
        double We  = 0.1 * deviceLength;

        double Re  = 0.5 * We;  // emitter doping radius - about 1/2 of We.
        double Rb  = 0.5 * Wad; // base doping radius - about 1/2  of Wad.

        double Ae;
        double Ab;

        if (given("TYPE") && deviceType=="NPN") // if npn
        {
          Ne  = 1.00e+19;
          Nb  = 1.00e+16;
          Nc  = 1.00e+14;
          Ae  = log(Ne/Nb)/(Re*Re);
          Ab  = log(Nb/Nc)/(Rb*Rb);

          for (i=0;i<numMeshPoints;++i)
          {
            double x = xVec[i];
            double y = yVec[i];
            y -= deviceWidth;

            CVec[i] =
              Nc -                                 // collector
              (Nb+Nc)*DopeInfo::ngdep(x,y,2.0*Wad,Ab,Ab) +     // base
              (Ne+Nb)*DopeInfo::ngdep(x,y,2.0*We ,Ae,Ae);      // emitter
          }
        }
        else  // if pnp
        {
          Ne  = 1.00e+19;
          Nb  = 1.00e+16;
          Nc  = 1.00e+14;
          Ae  = log(Ne/Nb)/(Re*Re);
          Ab  = log(Nb/Nc)/(Rb*Rb);

          for (i=0;i<numMeshPoints;++i)
          {
            double x = xVec[i];
            double y = yVec[i];
            y -= deviceWidth;

            CVec[i] =
              - Nc +                                 // collector
              (Nb+Nc)*DopeInfo::ngdep(x,y,2.0*Wad,Ab,Ab) -     // base
              (Ne+Nb)*DopeInfo::ngdep(x,y,2.0*We ,Ae,Ae);      // emitter
          }
        }

        // For the purposes of the rest of the code getting scaling
        // correct, etc., resetting Na, Nd here.  This is neccessary,
        // since the doping levels are hardwired for this problem.

        // Obtain Na and Nd:
        // Nd = number of donors (maximum + doping)
        // Na = number of acceptors (minimum + doping)
        double NaTmp = 1.0e+99;
        double NdTmp =-1.0e+99;

        for (i=0;i<numMeshPoints;++i)
        {
          if(NaTmp > CVec[i]) NaTmp = CVec[i];
          if(NdTmp < CVec[i]) NdTmp = CVec[i];
        }

        Na = fabs(NaTmp);
        Nd = fabs(NdTmp);
      }
      else if (numElectrodes == 4)
      {

        // Assuming a MOSFET
        //
        //     Source              Gate               Drain
        //
        //    |--5--|   |------------3------------|   |--5--|
        //
        // -  C-----D---E-------------------------F---G-----H  -
        // 6  |         |                         |         |  |
        // -  |         I-------------------------J         |  |
        // n  |                                             |  |Noflux
        // o  |                                             |  2
        // f  |                                             |  |
        // l  |                                             |  |
        // u  |                                             |  |
        // x  K---------------------------------------------L  -
        //
        //    |----------------------1----------------------|
        //                          Bulk

        // hardwired doping constants for now.
        // default is NMOS.

        // fix the geometrical stuff later...
        double WDEV = deviceLength;  // device width

        double WOX  = 0.7 * deviceLength; // oxide width

        double NS    = 1.0e16;       // substrate doping       (cm^-3)
        double NC    = 1.0e19;       // contact doping         (cm^-3)
        double WDIFF = (WDEV-WOX)/2; // diffusion width        (cm)
        double DDIFF = (0.125)*deviceWidth; // diffusion depth (cm)
        double DT    = 1.0e-11;      // diffusion coef. * time (cm^2)

        for (i=0;i<numMeshPoints;++i)
        {
          double x = xVec[i];
          double y = yVec[i];
          y -= deviceWidth;

          CVec[i] =
            -NS +
            (NC+NS)*DopeInfo::nsdep(x,     2*WDIFF,DT)*DopeInfo::nsdep(y,2*DDIFF,DT) +
            (NC+NS)*DopeInfo::nsdep(WDEV-x,2*WDIFF,DT)*DopeInfo::nsdep(y,2*DDIFF,DT) ;
        }

        if (given("TYPE") && deviceType=="PMOS") // if pmos, then invert the default.
        {
          for (i=0;i<numMeshPoints;++i)
          {
            CVec[i] = -CVec[i];
          }
        }

        // For the purposes of the rest of the code getting scaling
        // correct, etc., resetting Na, Nd here.  This is neccessary,
        // since the doping levels are hardwired for this problem.

        // Obtain Na and Nd:
        // Nd = number of donors (maximum + doping)
        // Na = number of acceptors (minimum + doping)
        double NaTmp = 1.0e+99;
        double NdTmp =-1.0e+99;

        for (i=0;i<numMeshPoints;++i)
        {
          if(NaTmp > CVec[i]) NaTmp = CVec[i];
          if(NdTmp < CVec[i]) NdTmp = CVec[i];
        }

        Na = fabs(NaTmp);
        Nd = fabs(NdTmp);
      }
      else
        // more than four electrodes?  not ready yet...
      {
        UserFatal(*this) << "Too many electrodes specified.  numElectrodes = " << numElectrodes;
      }
    }
    else // use the "new" doping specification:
    {
      // loop over the dope info map, and sum contributions from each
      // doping entity into the total doping array, CVec.
      std::map<std::string, DopeInfo *>::iterator iter;
      std::map<std::string, DopeInfo *>::iterator start = dopeInfoMap.begin();
      std::map<std::string, DopeInfo *>::iterator end   = dopeInfoMap.end();

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        Xyce::dout() << "dope info map:" << std::endl;
      }
#endif
      for ( iter = start; iter != end; ++iter )
      {
        DopeInfo & di = *(iter->second);

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0)
        {
          Xyce::dout() << di;
        }
#endif
        di.setupInfo2d(CVec,CdonorVec,CacceptorVec,xVec,yVec,devSupport);

      } // dope info map loop

      Na = 0.0;
      Nd = 0.0;
      for (i=0;i<numMeshPoints;++i)
      {
        if (Na > CVec[i]) Na = CVec[i];
        if (Nd < CVec[i]) Nd = CVec[i];
      }
      Na = fabs(Na);
      Nd = fabs(Nd);
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout().width(20);
    Xyce::dout().precision(12);
    Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << std::endl;
    Xyce::dout() << "Na = " << Na << std::endl;
    Xyce::dout() << "Nd = " << Nd << std::endl<< std::endl;
  }

  if (getDeviceOptions().debugLevel > 0)
  {
    for (int inode=0;inode<numMeshPoints;++inode)
    {
      Xyce::dout() << xVec[inode];
      Xyce::dout() << "  " << yVec[inode];
      Xyce::dout() << "  CVec["<<inode<<"] = " << CVec[inode] << std::endl;
    }
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
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
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::setupScalingVars ()
{
  bool bsuccess = true;

  double maxSize = meshContainerPtr->getMaxSize();

  if (given("X0"))
    scalingVars.x0 = x0_user;
  else
    scalingVars.x0  = maxSize;        // distance scaling (cm)

  // if the geometry is cylindrical, then scale electrode areas by area.
  // if cart., then scale by length.
  scalingVars.a0 = scalingVars.x0;
  if (meshContainerPtr->cylGeom) scalingVars.a0 *= scalingVars.x0;

  scalingVars.T0  = Temp;           // temperature scaling (K)  (not really used)

  // electrostatic potential scaling (V)
  scalingVars.V0 = Vt;
  scalingVars.rV0 = 1.0/scalingVars.V0;

  // concentration scaling (cm^-3);
  if (Na >= Nd) scalingVars.C0  = Na;
  else          scalingVars.C0  = Nd;

  scalingVars.D0  = 35.0;                   // diffusion coefficient scaling (cm^2/s)

  scalingVars.u0  = scalingVars.D0/scalingVars.V0;                  // mobility coefficient scaling (cm^2/V/s)

  scalingVars.R0  = scalingVars.D0*scalingVars.C0/(scalingVars.x0*scalingVars.x0);          // recombination rate scaling (cm^-3/s)
  scalingVars.rR0 = 1.0/scalingVars.R0;

  scalingVars.t0  = (scalingVars.x0*scalingVars.x0)/scalingVars.D0;             // time scaling (s)

  scalingVars.E0  = scalingVars.V0/scalingVars.x0;                  // electric field scaling (V/cm)

  scalingVars.F0  = scalingVars.D0*scalingVars.C0/scalingVars.x0;               // particle flux scaling (cm^-2/s)

  scalingVars.J0  = charge*scalingVars.D0*scalingVars.C0/scalingVars.x0;       // current density scaling (A/cm^2)

  scalingVars.L0  = scalingVars.V0*e0/(charge*scalingVars.x0*scalingVars.x0*scalingVars.C0);  // Laplacian scaling constant

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "  scalingVars.x0 = " << scalingVars.x0 << "\n";
    Xyce::dout() << "  scalingVars.a0 = " << scalingVars.a0 << "\n";
    Xyce::dout() << "  scalingVars.T0 = " << scalingVars.T0 << "\n";
    Xyce::dout() << "  scalingVars.V0 = " << scalingVars.V0 << "\n";
    Xyce::dout() << "  scalingVars.C0 = " << scalingVars.C0 << "\n";
    Xyce::dout() << "  scalingVars.D0 = " << scalingVars.D0 << "\n";
    Xyce::dout() << "  scalingVars.u0 = " << scalingVars.u0 << "\n";
    Xyce::dout() << "  scalingVars.R0 = " << scalingVars.R0 << "\n";
    Xyce::dout() << "  scalingVars.t0 = " << scalingVars.t0 << "\n";
    Xyce::dout() << "  scalingVars.E0 = " << scalingVars.E0 << "\n";
    Xyce::dout() << "  scalingVars.J0 = " << scalingVars.J0 << "\n";
    Xyce::dout() << "  scalingVars.L0 = " << scalingVars.L0 << std::endl;
    Xyce::dout() << std::endl;
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
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::scaleVariables ()
{
  bool bsuccess = true;
  int i;
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  Na /= scalingVars.C0;
  Nd /= scalingVars.C0;
  Ni /= scalingVars.C0;
  NnMax /= scalingVars.C0;
  NpMax /= scalingVars.C0;

  maxVoltDelta /= scalingVars.V0;

  VminExp /=scalingVars.V0;
  VmaxExp /=scalingVars.V0;

  Vbi       /= scalingVars.V0;

  // scale area, elen and ilen!
  bsuccess = meshContainerPtr->scaleMesh(scalingVars.x0);

  // scale boundary conditions:
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  for (;iterDI!=lastDI;++iterDI)
  {
    for (i=0;i<iterDI->numBoundaryPoints;++i)
    {
      iterDI->nnbcVec[i] /= scalingVars.C0;
      iterDI->npbcVec[i] /= scalingVars.C0;
      iterDI->VbcVec [i] /= scalingVars.V0;
      iterDI->VequVec[i] /= scalingVars.V0;
    }

    // If running with cylindrical coordinates, scale by area (scalingVars.x0*scalingVars.x0).
    // If running with cart. coordinates, scale by length (scalingVars.x0).
    iterDI->area /= scalingVars.a0;
    int size = iterDI->areaVector.size();
    for (i = 0; i < size; ++i)
    {
      iterDI->areaVector[i] /= scalingVars.a0;
    }
  }

  for (i=0;i<numMeshPoints;++i)
  {
    nnVec[i] /= scalingVars.C0;
    npVec[i] /= scalingVars.C0;
    CVec[i]  /= scalingVars.C0;
    VVec[i]  /= scalingVars.V0;
    unVec[i] /= scalingVars.u0;
    upVec[i] /= scalingVars.u0;
    tnVec[i] /= scalingVars.t0;
    tpVec[i] /= scalingVars.t0;
    xVec[i]  /= scalingVars.x0;
    yVec[i]  /= scalingVars.x0;

#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

    if( useVectorGIDFlag )
    {
      if (Vrowarray[i] != -1)
        (solVectorPtr)->setElementByGlobalIndex( Vrowarray[i], VVec[i], 0);
      if (Nrowarray[i] != -1)
        (solVectorPtr)->setElementByGlobalIndex( Nrowarray[i], nnVec[i], 0);
      if (Prowarray[i] != -1)
        (solVectorPtr)->setElementByGlobalIndex( Prowarray[i], npVec[i], 0);
    }
    else
    {
      (*solVectorPtr)[li_Vrowarray[i]] = VVec[i];
      (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
      (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
    }
  }

  for (i=0;i<meshContainerPtr->getNumEdges();++i)
  {
    unE_Vec[i] /= scalingVars.u0;
    upE_Vec[i] /= scalingVars.u0;
  }


  variablesScaled = true;

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : Instance::unScaleVariables
//
// Purpose       : Reverses the effect of scaleVariables.
//
// Special Notes : It was neccessary to add this after setting sensitivity
//                 calculations based on the doping.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/17/03
//-----------------------------------------------------------------------------
bool Instance::unScaleVariables ()
{
  bool bsuccess = true;
  int i;
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  Na *= scalingVars.C0;
  Nd *= scalingVars.C0;
  Ni *= scalingVars.C0;
  NnMax *= scalingVars.C0;
  NpMax *= scalingVars.C0;

  maxVoltDelta *= scalingVars.V0;

  VminExp *=scalingVars.V0;
  VmaxExp *=scalingVars.V0;

  Vbi     *= scalingVars.V0;

  // scale area, elen and ilen!
  bsuccess = meshContainerPtr->scaleMesh(1.0/scalingVars.x0);

  // scale boundary conditions:
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  for (;iterDI!=lastDI;++iterDI)
  {
    for (i=0;i<iterDI->numBoundaryPoints;++i)
    {
      iterDI->nnbcVec[i] *= scalingVars.C0;
      iterDI->npbcVec[i] *= scalingVars.C0;
      iterDI->VbcVec [i] *= scalingVars.V0;
      iterDI->VequVec[i] *= scalingVars.V0;
    }

    // If running with cylindrical coordinates, scale by area (scalingVars.x0*scalingVars.x0).
    // If running with cart. coordinates, scale by length (scalingVars.x0).
    iterDI->area *= scalingVars.a0;
    int size = iterDI->areaVector.size();
    for (i = 0; i < size; ++i)
    {
      iterDI->areaVector[i] *= scalingVars.a0;
    }
  }

  for (i=0;i<numMeshPoints;++i)
  {
    nnVec[i] *= scalingVars.C0;
    npVec[i] *= scalingVars.C0;
    CVec[i]  *= scalingVars.C0;
    VVec[i]  *= scalingVars.V0;
    unVec[i] *= scalingVars.u0;
    upVec[i] *= scalingVars.u0;
    tnVec[i] *= scalingVars.t0;
    tpVec[i] *= scalingVars.t0;
    xVec[i]  *= scalingVars.x0;
    yVec[i]  *= scalingVars.x0;

#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

    if( useVectorGIDFlag )
    {
      if (Vrowarray[i] != -1)
        (solVectorPtr)->setElementByGlobalIndex( Vrowarray[i], VVec[i], 0);
      if (Nrowarray[i] != -1)
        (solVectorPtr)->setElementByGlobalIndex( Nrowarray[i], nnVec[i], 0);
      if (Prowarray[i] != -1)
        (solVectorPtr)->setElementByGlobalIndex( Prowarray[i], npVec[i], 0);
    }
    else
    {
      (*solVectorPtr)[li_Vrowarray[i]] = VVec[i];
      (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
      (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
    }
  }

  for (i=0;i<meshContainerPtr->getNumEdges();++i)
  {
    unE_Vec[i] *= scalingVars.u0;
    upE_Vec[i] *= scalingVars.u0;
  }

  variablesScaled = false;

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : Instance::scaleDopeVariables
//
// Purpose       : This function performs scaling, similarly to
//                 function scaleVariables, but it only modifies variables
//                 related to the doping.  This is only used in the context
//                 of a sensitivity calculation, in which the sensitivity
//                 parameters that are being perturbed are doping params.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/28/03
//-----------------------------------------------------------------------------
bool Instance::scaleDopeVariables ()
{
  bool bsuccess = true;
  int i;

  Na /= scalingVars.C0;
  Nd /= scalingVars.C0;
  Ni /= scalingVars.C0;
  NnMax /= scalingVars.C0;
  NpMax /= scalingVars.C0;

  // scale boundary conditions:
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  for (;iterDI!=lastDI;++iterDI)
  {
    for (i=0;i<iterDI->numBoundaryPoints;++i)
    {
      iterDI->nnbcVec[i] /= scalingVars.C0;
      iterDI->npbcVec[i] /= scalingVars.C0;
      iterDI->VequVec[i] /= scalingVars.V0;
    }
  }

  for (i=0;i<numMeshPoints;++i)
  {
    CVec[i]  /= scalingVars.C0;
    unVec[i] /= scalingVars.u0;
    upVec[i] /= scalingVars.u0;
    tnVec[i] /= scalingVars.t0;
    tpVec[i] /= scalingVars.t0;
    nnVec[i] /= scalingVars.C0;
    npVec[i] /= scalingVars.C0;
    xVec[i]  /= scalingVars.x0;
    yVec[i]  /= scalingVars.x0;
  }

  for (i=0;i<meshContainerPtr->getNumEdges();++i)
  {
    unE_Vec[i] /= scalingVars.u0;
    upE_Vec[i] /= scalingVars.u0;
  }


  variablesScaled = true;

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : Instance::unScaleDopeVariables
//
// Purpose       : Reverses the effect of scaleDopeVariables.
//
// Special Notes : Almost everything needed to use a new, perturbed doping
//                 profile is just re-generated from scratch, so there
//                 isn't much to do here.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/28/03
//-----------------------------------------------------------------------------
bool Instance::unScaleDopeVariables ()
{
  bool bsuccess = true;
  int i;

  Na *= scalingVars.C0;
  Nd *= scalingVars.C0;
  Ni *= scalingVars.C0;
  NnMax *= scalingVars.C0;
  NpMax *= scalingVars.C0;

  for (i=0;i<numMeshPoints;++i)
  {
    nnVec[i] *= scalingVars.C0;
    npVec[i] *= scalingVars.C0;
    xVec [i] *= scalingVars.x0;
    yVec [i] *= scalingVars.x0;
  }

  variablesScaled = false;

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : Instance::calcInitialGuess
// Purpose       : This function calculates the initial e-, h+ densties
//                 and the intial voltage.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcInitialGuess ()
{
  bool bsuccess = true;
  int i;
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

#if 0
  N_LAS_Vector * solVectorPtr;
  N_LAS_Vector * currSolVectorPtr;
  N_LAS_Vector * lastSolVectorPtr;

  if (extData.initializeAllFlag == true)
  {
    solVectorPtr = extData.nextSolVectorPtr;
    currSolVectorPtr = extData.currSolVectorPtr;
    lastSolVectorPtr = extData.lastSolVectorPtr;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::calcInitialGuess"<< std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

  // set up an initial guess for nnVec and npVec, based on the doping profile,
  // and the equilibrium density expressions.  Place these in the
  // solution vector.

  double tmp(0.0);
  double Ci(0.0);
  double Cisq(0.0), Nisq(0.0);

  for (i=0;i<numMeshPoints;++i)
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

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << "nnVec[" << i << "] = " << nnVec[i];
      Xyce::dout() << "  npVec[" << i << "] = " << npVec[i];
      Xyce::dout() << std::endl;
    }
#endif
  }

  for (i=0;i<numMeshPoints;++i)
  {
    // the doping is n-type.
    if (Ci>0)
    {
      VVec[i] = + Vt * log(nnVec[i]/Ni);
    }
    // the doping is p-type.
    else
    {
      VVec[i] = - Vt * log(npVec[i]/Ni);
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << "VVec[" << i << "] = " << VVec[i] << std::endl;
    }
#endif
  }

  // Scale and offset this initial guess for V so that it matches the
  // Vequ boundary conditions.
  double VminTmp = +1.0e+99;
  double VmaxTmp = -1.0e+99;

  double VminBC  = +1.0e+99;
  double VmaxBC  = -1.0e+99;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);
    std::vector<int>::iterator iterNode = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastNode = labelPtr->mNodeVector.end ();

    for ( ;iterNode!=lastNode;++iterNode)
    {
      if (VminTmp > VVec[*iterNode]) VminTmp = VVec[*iterNode];
      if (VmaxTmp < VVec[*iterNode]) VmaxTmp = VVec[*iterNode];

      int ilocal = iterDI->meshGlobalToLocal[*iterNode];
      if (VminBC  > iterDI->VequVec[ilocal]) VminBC  = iterDI->VequVec[ilocal];
      if (VmaxBC  < iterDI->VequVec[ilocal]) VmaxBC  = iterDI->VequVec[ilocal];
    }
  }

  double Voffset  = -VminTmp;
  double VtotDiff = VmaxTmp-VminTmp;
  double VBCDiff  = VmaxBC-VminBC;
  double Vscale   = (VtotDiff!=0.0)?(VBCDiff/VtotDiff):1.0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "Voffset  = " << Voffset <<  std::endl;
    Xyce::dout() << "VtotDiff = " << VtotDiff << std::endl;
    Xyce::dout() << "VBCDiff  = " << VBCDiff << std::endl;
    Xyce::dout() << "Vscale   = " << Vscale << std::endl;
  }
#endif

  for (i=0;i<numMeshPoints;++i)
  {
    VVec[i] += Voffset;
  }

  for (i=0;i<numMeshPoints;++i)
  {
    VVec[i] *= Vscale;
  }

  // Place all these initial values into the solution vector.
  for (i=0;i<numMeshPoints;++i)
  {

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << "VVec[" << i << "] = " << VVec[i] << std::endl;
    }
#endif

#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

    if (extData.initializeAllFlag == true)
    {
      if( useVectorGIDFlag )
      {
        if (Vrowarray[i] != -1)
          (solVectorPtr)->setElementByGlobalIndex( Vrowarray[i], VVec[i], 0);
        if (Nrowarray[i] != -1)
          (solVectorPtr)->setElementByGlobalIndex( Nrowarray[i], nnVec[i], 0);
        if (Prowarray[i] != -1)
          (solVectorPtr)->setElementByGlobalIndex( Prowarray[i], npVec[i], 0);

#if 0
        if (Vrowarray[i] != -1)
          (currSolVectorPtr)->setElementByGlobalIndex( Vrowarray[i], VVec[i], 0);
        if (Nrowarray[i] != -1)
          (currSolVectorPtr)->setElementByGlobalIndex( Nrowarray[i], nnVec[i], 0);
        if (Prowarray[i] != -1)
          (currSolVectorPtr)->setElementByGlobalIndex( Prowarray[i], npVec[i], 0);

        if (Vrowarray[i] != -1)
          (lastSolVectorPtr)->setElementByGlobalIndex( Vrowarray[i], VVec[i], 0);
        if (Nrowarray[i] != -1)
          (lastSolVectorPtr)->setElementByGlobalIndex( Nrowarray[i], nnVec[i], 0);
        if (Prowarray[i] != -1)
          (lastSolVectorPtr)->setElementByGlobalIndex( Prowarray[i], npVec[i], 0);
#endif
      }
      else
      {
        (*solVectorPtr)[li_Vrowarray[i]] = VVec[i];
        (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
        (*solVectorPtr)[li_Prowarray[i]] = npVec[i];

#if 0
        (*currSolVectorPtr)[li_Vrowarray[i]] = VVec[i];
        (*currSolVectorPtr)[li_Nrowarray[i]] = nnVec[i];
        (*currSolVectorPtr)[li_Prowarray[i]] = npVec[i];

        (*lastSolVectorPtr)[li_Vrowarray[i]] = VVec[i];
        (*lastSolVectorPtr)[li_Nrowarray[i]] = nnVec[i];
        (*lastSolVectorPtr)[li_Prowarray[i]] = npVec[i];
#endif
      }
    }
  }

  // get the maximum and minimum potentials.
  VmaxExp = -1.0e99;
  VminExp = +1.0e99;

  for (int j=0;j<numMeshPoints;++j)
  {
    if (VmaxExp < VVec[j]) VmaxExp = VVec[j];
    if (VminExp > VVec[j]) VminExp = VVec[j];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl;
    Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "Vmax = " << VmaxExp << std::endl;
    Xyce::dout() << "Vmin = " << VminExp << std::endl;
    Xyce::dout() << std::endl;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcVequBCs
// Purpose       : This function sets up the "Vequ" boundary condition
//                 variables for each electrode.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/28/03
//-----------------------------------------------------------------------------
bool Instance::calcVequBCs ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::calcVequBCs"<< std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
    Xyce::dout() << "  Na = "<< Na << std::endl;
    Xyce::dout() << "  Nd = "<< Nd << std::endl;
  }
#endif

  double VminBC =+1.0e+99;
  double VmaxBC =-1.0e+99;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << "DI name = " << iterDI->eName << std::endl;
      Xyce::dout() << "material = " << iterDI->material << std::endl;
      if (iterDI->materialGiven)
      {
        Xyce::dout() << "material was given" << std::endl;
      }
      else
      {
        Xyce::dout() << "material was NOT given" << std::endl;
      }

      if (iterDI->oxideBndryFlag)
      {
        Xyce::dout() << "This is a boundary WITH an oxide." << std::endl;
      }
      else
      {
        Xyce::dout() << "This is a boundary WITHOUT an oxide." << std::endl;
      }

      Xyce::dout().setf(std::ios::scientific);
    }
#endif

    std::string insul = "sio2";  // oxide layers can only be sio2...

    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);
    std::vector<int>::iterator iterNode = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastNode = labelPtr->mNodeVector.end ();

    for (i=0;i<iterDI->numBoundaryPoints;++i,++iterNode)
    {
      double Ci = CVec[*iterNode];
      double ns = nnVec[*iterNode];  // needed for oxide interface potential
      double Cisq = Ci*Ci;
      double Nisq = Ni*Ni;  // Ni is the intrinsic concentration
      double tmp, nnTmp, npTmp;

      // equilibrium electron concentration:
      tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
      nnTmp = ((Ci>=0)?(tmp):(0.0)) + ((Ci<0)?(Nisq/tmp):(0.0));

      // equilibrium hole concentration:
      tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
      npTmp = ((Ci<=0)?(tmp):(0.0)) + ((Ci>0)?(Nisq/tmp):(0.0));

      ExtendedString mater = iterDI->material;
      mater.toLower();

      if (mater=="neutral")
      {
        // the doping is n-type.
        if (Ci>0)
        {
          iterDI->VequVec[i] = + Vt * log(nnTmp/Ni);
        }
        // the doping is p-type.
        else
        {
          iterDI->VequVec[i] = - Vt * log(npTmp/Ni);
        }
      }
      else // this electrode has a schottky barrier.
      {
        // the doping is n-type and a metal contact.
        if (Ci>0)
        {
          iterDI->VequVec[i] = + matSupport.workfunc(iterDI->material)
                               - matSupport.affin(bulkMaterial)
                               - 0.5 * matSupport.bandgap(bulkMaterial, Temp)
                               + 2.0 * Vt * log(nnTmp/Ni);
        }
        else if (Ci<=0)
        {
          iterDI->VequVec[i] = + matSupport.workfunc(iterDI->material)
                               - matSupport.affin(bulkMaterial)
                               - 0.5 * matSupport.bandgap(bulkMaterial, Temp)
                               - 2.0 * Vt * log(npTmp/Ni);
        }
      }

      // if this is a metal-oxide-semiconductor boundary, add an offset due
      // to the charge stored in the oxide.
      //
      //  V_cap = Qi/Ci
      if (iterDI->oxideBndryFlag)
      {
        iterDI->VequVec[i] += - iterDI->oxcharge * charge *iterDI->oxthick/
                              (matSupport.getRelPerm(insul) * e0);
      }

      iterDI->VbcVec [i] = 0.0;

      if (VminBC > iterDI->VequVec[i]) VminBC = iterDI->VequVec[i];
      if (VmaxBC < iterDI->VequVec[i]) VmaxBC = iterDI->VequVec[i];

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        Xyce::dout() << "Vequ["<<i<<"]=" << iterDI->VequVec[i];
        Xyce::dout() << std::endl;
      }
#endif
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "After offset:" << std::endl;
  }
#endif

  // offset these boundary conditions so that the minimum is zero.
  double Voffset = -VminBC;

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << "DI name = " << iterDI->eName << std::endl;
      Xyce::dout().setf(std::ios::scientific);
    }
#endif

    for (i=0;i<iterDI->numBoundaryPoints;++i)
    {
      iterDI->VequVec[i] += Voffset;
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        Xyce::dout() << "VequVec["<<i<<"] = " << iterDI->VequVec[i] << std::endl;
      }
#endif
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

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
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool Instance::calcDensityBCs ()
{
  bool bsuccess = true;
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider<< std::endl;
    Xyce::dout() << std::endl << "In Instance::calcDensityBCs" << std::endl;
  }
#endif

  NnMax = -1.0e+99;
  NpMax = -1.0e+99;

  NnMin = +1.0e+99;
  NpMin = +1.0e+99;

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << iterDI->eName<< ":" << std::endl;
    }
#endif

    // This expression is from Selberherr, enforcing thermal
    // equilibrium and vanishing space charge at ohmic contacts.
    for (int i=0;i<iterDI->numBoundaryPoints;++i)
    {
      mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);
      int nIndex;
      if (constBoundaryFlag)
        nIndex = labelPtr->mNodeVector[0];
      else
        nIndex = labelPtr->mNodeVector[i];

      iterDI->nnbcVec[i] =
        0.5*(sqrt(CVec[nIndex]*CVec[nIndex]+4*Ni*Ni)+CVec[nIndex]);

      iterDI->npbcVec[i] =
        0.5*(sqrt(CVec[nIndex]*CVec[nIndex]+4*Ni*Ni)-CVec[nIndex]);

      if (NnMax < iterDI->nnbcVec[i]) NnMax = iterDI->nnbcVec[i];
      if (NpMax < iterDI->npbcVec[i]) NpMax = iterDI->npbcVec[i];

      if (NnMin > iterDI->nnbcVec[i]) NnMin = iterDI->nnbcVec[i];
      if (NpMin > iterDI->npbcVec[i]) NpMin = iterDI->npbcVec[i];

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        Xyce::dout() << "\tCVec["<<nIndex<<"] = " << CVec[nIndex]<< std::endl;
        Xyce::dout() << "\tnnbc["<<i<<"] = " << iterDI->nnbcVec[i] << std::endl;
        Xyce::dout() << "\tnpbc["<<i<<"] = " << iterDI->npbcVec[i] << std::endl;
      }
#endif
    }

    // Set the boundaries to reflect nnbc and npbc.  This is neccessary so
    // that the correct Vequ is calculated in function calcVequBCs.
    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);
    std::vector<int>::iterator iterNode = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastNode = labelPtr->mNodeVector.end ();

    for ( ;iterNode!=lastNode;++iterNode)
    {
      int i1 = iterDI->meshGlobalToLocal[*iterNode];
      nnVec[*iterNode] = iterDI->nnbcVec[i1]/scalingVars.C0;
      npVec[*iterNode] = iterDI->npbcVec[i1]/scalingVars.C0;
    }
  }

  if (NnMin <= 0) NnMin = 1.56269e-9;  // just a guess.
  if (NpMin <= 0) NpMin = 1.56269e-9;  // just a guess.

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "NnMax = " << NnMax << std::endl;
    Xyce::dout() << "NpMax = " << NpMax << std::endl;
    Xyce::dout() << "NnMin = " << NnMin << std::endl;
    Xyce::dout() << "NpMin = " << NpMin << std::endl;
    Xyce::dout() << section_divider<< std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcBoundaryConditions
// Purpose       : This function sets up the boundary condition variables
//                 for each electrode.  This function only handles the
//                 voltage BC's.  Density BC's are only set up once, in
//                 function calcDensityBCs.
//
// Special Notes : If a continuation method is being used, the electrode
//                 voltages are being slowly varried.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool Instance::calcBoundaryConditions ()
{
  bool bsuccess = true;
  int i;
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << std::endl << "In Instance::calcBoundaryConditions" << std::endl;
  }
#endif

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    if (getSolverState().PDEcontinuationFlag)
    {
      for (i=0;i<iterDI->numBoundaryPoints;++i)
      {
	iterDI->VbcVec[i] = iterDI->Vckt_ramp + iterDI->VequVec[i];
      }
    }
    else
    {
      for (i=0;i<iterDI->numBoundaryPoints;++i)
      {
	iterDI->VbcVec[i] = iterDI->Vckt + iterDI->VequVec[i];
      }
    }

#ifdef Xyce_NEW_BC
    // if using the "new" boundary conditions, go through the V,n,p
    // arrays and impose the electrode boundary values.
    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);
    std::vector<int>::iterator iterNode = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastNode = labelPtr->mNodeVector.end ();

    for ( ;iterNode!=lastNode;++iterNode)
    {
      int i1 = iterDI->meshGlobalToLocal[*iterNode];
      VVec[*iterNode] = iterDI->VbcVec [i1];
      nnVec[*iterNode] = iterDI->nnbcVec[i1];
      npVec[*iterNode] = iterDI->npbcVec[i1];
    }
#endif // Xyce_NEW_BC

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << iterDI->eName << ":" << std::endl;
      for (i=0;i<iterDI->numBoundaryPoints;++i)
      {
        Xyce::dout()<<"Vbc["<<i<<"] = "<<iterDI->VbcVec[i]<<", "<<iterDI->VbcVec[i]*scalingVars.V0;
        Xyce::dout() << " nnbc["<<i<<"] = "<< iterDI->nnbcVec[i];
        Xyce::dout() << ", "<<iterDI->nnbcVec[i]*scalingVars.C0;
        Xyce::dout() << " npbc["<<i<<"] = "<< iterDI->npbcVec[i];
        Xyce::dout() << ", "<<iterDI->npbcVec[i]*scalingVars.C0;
        Xyce::dout() << std::endl;
      }
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -5 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

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
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << std::endl << "In Instance::obtainNodeVoltages" << std::endl;
  }
#endif

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {

    if( useVectorGIDFlag )
    {
      if (iterDI->gid != -1)
        iterDI->Vckt = (solVectorPtr)->getElementByGlobalIndex(iterDI->gid, 0);
      else
        iterDI->Vckt = 0.0;
    }
    else
    {
      iterDI->Vckt = (*solVectorPtr)[iterDI->lid];
    }

    iterDI->Vckt /= scalingVars.V0;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << iterDI->eName << "  Vckt = " << iterDI->Vckt;
      Xyce::dout() << "  Vckt*scalingVars.V0 = " << (iterDI->Vckt * scalingVars.V0) << std::endl;
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
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
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/15/02
//-----------------------------------------------------------------------------
bool Instance::applyVoltageLimiting ()
{
  bool bsuccess = true;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "device:    " << getName() << std::endl;
  }
#endif

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    double v1      = iterDI->Vckt * scalingVars.V0;
    double v1_old  = iterDI->Vckt_old * scalingVars.V0;
#ifdef Xyce_DEBUG_DEVICE
    double v1_orig = v1;
#endif

    double delV1   = v1 - v1_old;

    if ( delV1 >  1.25 ) v1 = v1_old + 1.25;
    if ( delV1 < -0.75 ) v1 = v1_old - 0.75;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "electrode  = " << iterDI->eName << std::endl;
      Xyce::dout() << "v1         = " << v1 << std::endl;
      Xyce::dout() << "v1_old     = " << v1_old << std::endl;
      Xyce::dout() << "v1_orig    = " << v1_orig << std::endl;
      Xyce::dout() << "v1/scalingVars.V0      = " << v1/scalingVars.V0 << std::endl;
      Xyce::dout() << "v1_old/scalingVars.V0  = " << v1_old/scalingVars.V0 << std::endl;
      Xyce::dout() << "v1_orig/scalingVars.V0 = " << v1_orig/scalingVars.V0 << std::endl;
      Xyce::dout() << std::endl;
    }
#endif

    iterDI->Vckt       = v1/scalingVars.V0;
    iterDI->Vckt_final = v1/scalingVars.V0;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::obtainSolution
// Purpose       : This function extracts V, nn, and np from the solution
//                 vector and copies them into local arrays.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::obtainSolution ()
{
  bool bsuccess = true;
  bool bs1 = true;

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::obtainSolution\n";
    Xyce::dout() << "solVectorPtr = " << solVectorPtr << std::endl;
  }
#endif

  bs1 = obtainNodeVoltages ();
  bsuccess = bsuccess && bs1;

  // set up the solution array:
  int i;
  for (i=0;i<numMeshPoints;++i)
  {

#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

    vOwnVec[i] = (Vrowarray[i] != -1);

    if( useVectorGIDFlag )
    {
      if (Vrowarray[i] != -1)
      {
        VVec[i] = (solVectorPtr)->getElementByGlobalIndex(Vrowarray[i], 0);
      }
    }
    else
    {
      VVec[i] = (*solVectorPtr)[li_Vrowarray[i]];
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    for (i=0;i<numMeshPoints;++i)
    {
      Xyce::dout() << "VVec["<<i<<"]=\t";
      Xyce::dout().width(20);
      Xyce::dout().precision(12);
      Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << VVec[i];
      Xyce::dout() << "  " << VVec[i]*scalingVars.V0;
      Xyce::dout() << "\n";
    }
  }
#endif

  // If the previous solution is from the nonlinear Poisson solution,
  // then calculate what the electron and hole densities must be, and
  // place them into the solution vector.

  // If we are past the nonlinear Poisson phase, then simply obtain
  // nn and np from the solution vector and move on.

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "  About to get the density.\n";
    if (getSolverState().dcopFlag) Xyce::dout() << "DCOP load" << std::endl;
    else                   Xyce::dout() << "Transient load" << std::endl;
    Xyce::dout() << "  doubleDCOPStep = " << getSolverState().doubleDCOPStep << "\n";
  }
#endif

  if ((getSolverState().dcopFlag) && getSolverState().doubleDCOPStep==0)
  {
    calcVoltDepDensities ();

    for (i=0;i<numMeshPoints;++i)
    {
#ifdef Xyce_NEW_BC
      if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

      nnOwnVec[i] = (Nrowarray[i] != -1);
      npOwnVec[i] = (Prowarray[i] != -1);

      if( useVectorGIDFlag )
      {
        if (Nrowarray[i] != -1)
          (solVectorPtr)->setElementByGlobalIndex(Nrowarray[i],nnVec[i]);
        if (Prowarray[i] != -1)
          (solVectorPtr)->setElementByGlobalIndex(Prowarray[i],npVec[i]);
      }
      else
      {
        (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
        (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
      }
    }
  }
  else
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout()<<"  Obtaining densities from solution vector (not voltage dep.)\n";
    }
#endif
    for (i=0;i<numMeshPoints;++i)
    {

#ifdef Xyce_NEW_BC
      if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC
      nnOwnVec[i] = (Nrowarray[i] != -1);
      npOwnVec[i] = (Prowarray[i] != -1);

      if( useVectorGIDFlag )
      {
        if (Nrowarray[i] != -1)
          nnVec[i] = (solVectorPtr)->getElementByGlobalIndex(Nrowarray[i], 0);
      }
      else
      {
        nnVec[i] = (*solVectorPtr)[li_Nrowarray[i]];
      }

#ifdef Xyce_PDE_DENSITY_CONSTRAINT
      // if the density is less than zero, force to be zero.
      if (nnVec[i] < 0.0)
      {
        nnVec[i] = 0.0;
        if( useVectorGIDFlag )
        {
          if (Nrowarray[i] != -1)
          {
            (solVectorPtr)->setElementByGlobalIndex(Nrowarray[i], 0.0, 0);
          }
        }
        else
        {
          (*solVectorPtr)[li_Nrowarray[i]] = 0.0;
        }
      }
#endif

      if( useVectorGIDFlag )
      {
        if (Prowarray[i] != -1)
        {
          npVec[i] = (solVectorPtr)->getElementByGlobalIndex(Prowarray[i], 0);
        }
      }
      else
      {
        npVec[i] = (*solVectorPtr)[li_Prowarray[i]];
      }

#ifdef Xyce_PDE_DENSITY_CONSTRAINT
      // if the density is less than zero, force to be zero.
      if (npVec[i] < 0.0)
      {
	npVec[i] = 0.0;
        if( useVectorGIDFlag )
        {
          if (Prowarray[i] != -1)
            (solVectorPtr)->setElementByGlobalIndex(Prowarray[i], 0.0, 0);
        }
        else
        {
          (*solVectorPtr)[li_Prowarray[i]] = 0;
        }
      }
#endif
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    for (i=0;i<numMeshPoints;++i)
    {
      Xyce::dout() << "nnVec["<<i<<"]=\t";
      Xyce::dout().width(14);
      Xyce::dout().precision(6);
      Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << nnVec[i];
      Xyce::dout() << "  " << nnVec[i]*scalingVars.C0;

      Xyce::dout() << "  npVec["<<i<<"]=\t";
      Xyce::dout().width(14);
      Xyce::dout().precision(6);
      Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << npVec[i];
      Xyce::dout() << "  " << npVec[i]*scalingVars.C0;

      Xyce::dout() << std::endl;
    }
  }
#endif

  // now set boundary conditions:
  // if the circuit is coupled to the PDE device, then bc's
  // must be updated everytime.
  //
  // If the circuit and PDE device are not coupled, then the
  // circuit node voltages can be considered constant, and the
  // BC's only need updating at the first Newton step.

  if (!(getSolverState().twoLevelNewtonCouplingMode==INNER_PROBLEM))
  {
    bs1 = calcBoundaryConditions ();
    bsuccess = bsuccess && bs1;
  }
  else  // ... if NOT coupled
  {
    if (getSolverState().newtonIter == 0)
    {
      bs1 = calcBoundaryConditions ();
      bsuccess = bsuccess && bs1;
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcRecombination
// Purpose       :
// Special Notes : This function assumes scaling is turned on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcRecombination ()
{
  bool bsuccess = true;

  int i;
  double Rsrh, Raug;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcRecombination\n";
    Xyce::dout() << "\n";
  }
#endif

  for (i=0;i<numMeshPoints;++i)
  {
    double n  = nnVec[i];
    double p  = npVec[i];
    double tn = tnVec[i];
    double tp = tpVec[i];

    // assuming Si for now.
    Rsrh = matSupport.calcRsrh (bulkMaterial, Ni,n,p,tn,tp);
    Raug = matSupport.calcRaug (bulkMaterial, Ni,n,p);

    RVec[i] = (Rsrh + Raug);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout().precision(4);
      Xyce::dout() << " nnVec="<<n<<" npVec="<<p;
      Xyce::dout() << " tnVec="<<tn<<" tpVec="<<tp;
      Xyce::dout() << " Rsrh="<<Rsrh;
      Xyce::dout() << " Raug="<<Raug;
      Xyce::dout() << " RVec["<<i<<"]="<<RVec[i];
      Xyce::dout() << "\n";
    }

    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout().precision(4);
      Xyce::dout() << " RVec["<<i<<"]="<<RVec[i];
      Xyce::dout() << std::endl;
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::setupPhotogen
//
// Purpose       : Sets up the photogen source.  For now, this is a simple
//                 uniform source.
//
// Special Notes : Not finished yet.   This is an early implementation.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/21/03
//-----------------------------------------------------------------------------
bool Instance::setupPhotogen ()
{
  bool bsuccess = true;
  Param ndParam;
  Util::Expression expr;

  // for now just supporting "pulse" and "const".
  switch (photoType)
  {
    case _PULSE_DATA:
      ndParam.setTag("PULSEDATA");
      expr.set("spice_pulse(V1, V2, TD, TR, TF, PW, PER)");
      expr.make_constant("V1", 0.0);
      expr.make_constant("V2", photoA1);
      expr.make_constant("TD", photoTd);
      expr.make_constant("TR", photoTr);
      expr.make_constant("TF", photoTf);
      expr.make_constant("PW", photoPw);
      expr.make_constant("PER", photoPer);
      ndParam.setVal(expr);
      setDependentParameter (ndParam, &PulseData, ParameterType::TIME_DEP);
      break;

    case _DC_DATA:
      ndParam.setTag("PULSEDATA");
      expr.set("CONST");
      expr.make_constant("CONST", photoA1);
      ndParam.setVal(expr);
      setDependentParameter (ndParam, &PulseData, ParameterType::TIME_DEP);
      break;

    default:
      break;
  }


  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::calcPhotogen
//
// Purpose       : Calculates a term similar to recombination (in that it
//                 is used the same way, on the RHS of the continuity
//                 equations).  This is a generation term, rather than a
//                 recombination term, so it should have the opposite sign
//                 (+) than that of RVec.
//
//                 My goal in setting this up is to start mimicing the
//                 capability of  medici/davicini.  Most of the user params
//                 are similar to those of medici, except for there being
//                 a "ph" precdeding the original name.
//
//                 For now just doing a uniform source.
//
// Special Notes : Assumes scaling is ON.
//
//                 Not finished yet.   This is an early implementation.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/21/03
//-----------------------------------------------------------------------------
bool Instance::calcPhotogen ()
{
  bool bsuccess = true;
  int i;

  double time = getSolverState().currTime;
  double val;

#ifdef Xyce_DEBUG_DEVICE

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcPhotogen\n";
    Xyce::dout() << "time = " << time;

    if (photogenOnFlag &&
        time >= photoTstart &&
        time <= photoTstop)
    { Xyce::dout() << "  photogen is on right now." << std::endl; }
    else
    { Xyce::dout() << "  photogen is off right now." << std::endl; }

    Xyce::dout() << "\n";
  }
#endif

  for (i=0;i<numMeshPoints;++i)
  {
    SVec[i] = 0.0;
  }

  if (photogenOnFlag &&
      time >= photoTstart &&
      time <= photoTstop)
  {
    if (getSolverState().PDEcontinuationFlag && !(photoContinuationFinished))
    {
      val = photoA1_ramp;
    }
    else
    {
      val = PulseData;
    }

    for (i=0;i<numMeshPoints;++i)
    {
      SVec[i] = (val * scalingVars.rR0);
    }

  } // photogenOnFlag, etc.

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Value of photogen source = " << val << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcPenalty
//
// Purpose       : Calculates extra penalty method term to be added to the
//                 continuity equations.
//
// Special Notes : This should represent a smooth function which is a
//                 function of electron/hole density.  If density is
//                 negative, the penalty function should return a large
//                 number.  If positive, it should return a zero.  It needs
//                 to be differentiable w.r.t. density.
//
//                 For now, penalty function is:
//
//                 pf = prefac*pow(dens,N)  if dens <  0.0
//                    = 0.0                 if dens >= 0.0
//
//                 The default value for prefac is scalingVars.C0, and the default
//                 value for N is 2.0.  The user can set both, however.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/03
//-----------------------------------------------------------------------------
bool Instance::calcPenalty ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  bool nonZeroFlag = false;
#endif

  double prefac;
  double power;

  // prefactor doesn't currently get used, but may later...
  if (penaltyPrefacGiven)
  {
    prefac = penaltyPrefac;
  }
  else
  {
    prefac = 1.0;
  }

  if (penaltyPowGiven)
  {
    power = penaltyPow;
  }
  else
  {
    power = 2.0;
  }

  for (i=0;i<numMeshPoints;++i)
  {
    double elec = nnVec[i]*scalingVars.C0;
    double hole = npVec[i]*scalingVars.C0;

    // electron penalty
    double epen = prefac*pow(elec,power);
    elecPenalty[i] = (nnVec[i]<0.0)?epen:0.0;

    // hole penalty
    double hpen = prefac*pow(hole,power);
    holePenalty[i] = (npVec[i]<0.0)?hpen:0.0;

#ifdef Xyce_DEBUG_DEVICE
    if (elecPenalty[i] != 0.0 ||
        holePenalty[i] != 0.0)
      nonZeroFlag = true;
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && nonZeroFlag)
  {
    for (i=0;i<numMeshPoints;++i)
    {
      Xyce::dout() << "  e- Pen["<<i<<"] = " << elecPenalty[i];
      Xyce::dout() << "  e- den["<<i<<"] = " << nnVec[i];

      Xyce::dout() << "  h+ Pen["<<i<<"] = " << holePenalty[i];
      Xyce::dout() << "  h+ den["<<i<<"] = " << npVec[i];
      Xyce::dout() << std::endl;
    }
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdPenalty
//
// Purpose       : Calculates extra penalty method derivative term to be added to the
//                 continuity equation matrix rows.
//
// Special Notes : This should represent a smooth function which is a
//                 function of electron/hole density.  If density is
//                 negative, the penalty function should return a large
//                 number.  If positive, it should return a zero.  It needs
//                 to be differentiable w.r.t. density.
//
//                 For now, penalty function is:
//
//                 pf = prefac*pow(dens,N)  if dens <  0.0
//                    = 0.0                 if dens >= 0.0
//
//                 The default value for prefac is scalingVars.C0, and the default
//                 value for N is 2.0.  The user can set both, however.
//
//                 d(pf)
//                 -------  = N*prefac*pow(dens,N-1)  if dens < 0.0
//                 d(dens)
//
//                 d(pf)
//                 -------  = 0.0                     if dens >= 0.0
//                 d(dens)
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/03
//-----------------------------------------------------------------------------
bool Instance::pdPenalty ()
{
  bool bsuccess = true;
  int i;

  double prefac;
  double power;

  // prefactor doesn't currently get used, but may later...
  if (penaltyPrefacGiven)
  {
    prefac = penaltyPrefac;
  }
  else
  {
    prefac = 1.0;
  }

  if (penaltyPowGiven)
  {
    power = penaltyPow;
  }
  else
  {
    power = 2.0;
  }

  double power2 = power - 1.0;

  for (i=0;i<numMeshPoints;++i)
  {
    double elec = nnVec[i];
    double hole = npVec[i];
    double Cpow = pow(scalingVars.C0,power);

    // electron penalty
    double epen = power*prefac*Cpow*pow(elec,power2);
    pdElecPenalty[i] = (nnVec[i]<0.0)?epen:0.0;

    // hole penalty
    double hpen = power*prefac*Cpow*pow(hole,power2);
    pdHolePenalty[i] = (npVec[i]<0.0)?hpen:0.0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::sumSources
//
// Purpose       : This function sums in all recombination/generation sources.
//                 This will probably be refactored later - right now it
//                 doesn't check what sources are enabled and what aren't -
//                 it just sums them all in.
//
// Special Notes : Assumes scaling is ON.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/21/03
//-----------------------------------------------------------------------------
bool Instance::sumSources ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE

  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::sumSources\n";
    Xyce::dout() << "\n";
  }
#endif

  for (i=0;i<numMeshPoints;++i)
  {
    // Positive R, negative S, because they are assumed to be recombination
    // terms.  (-totSrcVec) will be the source term.
    totSrcVec[i] = RVec[i] - SVec[i];
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "RVec["<<i<<"] = " << RVec[i];
      Xyce::dout() << "  SVec["<<i<<"] = " << SVec[i];
      Xyce::dout() << "  totSrcVec["<<i<<"] = " << totSrcVec[i] << "\n";
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdRecombination
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with the recombination term.
// Special Notes : This function assumes scaling is turned on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::pdRecombination ()
{
  bool bsuccess = true;

  int i;

  double dRsrhdn;
  double dRsrhdp;
  double dRaugdn;
  double dRaugdp;

#if 0
  double Cn,Cp;
  double pn;
  double A1, B1, C1;
  double dAdn, dAdp;
  double dBdn, dBdp;
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::pdRecombination\n";
    Xyce::dout() << "\n";
  }
#endif

  for (i=0;i<numMeshPoints;++i)
  {
#if 0
    Cn = 2.8e-31;
    Cp = 1.2e-31;
    pn = Ni*Ni;

    // Rsrh derivatives: checklater.
    // (Rsrch = A1*B1)

    double arg = CONSTMAX_EXP_ARG;

    A1 = (nnVec[i]*npVec[i]-pn);
    if (A1 >= exp(arg)) A1 = exp(arg);

    dAdn = (npVec[i]);
    dAdp = (nnVec[i]);

    C1 = (tpVec[i]*(nnVec[i]+Ni)+tnVec[i]*(npVec[i]+Ni));
    if (C1 >= exp(arg)) C1 = exp(arg);

    B1 = 1.0/C1;
    dBdn = -1.0/(C1*C1) * tpVec[i];
    dBdp = -1.0/(C1*C1) * tnVec[i];

    dRsrhdn = dAdn * B1 + dBdn * A1;
    dRsrhdp = dAdp * B1 + dBdp * A1;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  dRsrhdn = " << dRsrhdn << "  dRsrhdp = " << dRsrhdp;
      Xyce::dout() << "\n";
      Xyce::dout() << "  nnVec = " << nnVec[i] << "  npVec = " << npVec[i] << "\n";
      Xyce::dout() << "  tnVec = " << tnVec[i] << "  tpVec = " << tpVec[i] << "\n";
      Xyce::dout() << "  A1 = " << A1 << "  B1 = " << B1 << "  C1 = " << C1 << "\n";
      Xyce::dout() << "  dAdn = " << dAdn << "  dBdn = " << dBdn << "\n";
      Xyce::dout() << "  dAdp = " << dAdp << "  dBdp = " << dBdp << "\n";
    }
#endif

    // Raug derivatives:
    // Raug = A1*B1, where A1 is the same as it was for Rsrch).
    B1 = (Cn*nnVec[i]+Cp*npVec[i]);
    if (B1 >= exp(arg)) B1 = exp(arg);

    dBdn = Cn;
    dBdp = Cp;

    dRaugdn = dAdn*B1 + A1*dBdn;
    dRaugdp = dAdp*B1 + A1*dBdp;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  dRaugdn = " << dRaugdn << "  dRaugdp = " << dRaugdp <<"\n";
      Xyce::dout() << "  A1 = " << A1 << "  B1 = " << B1 << "\n";
      Xyce::dout() << "  dAdn = " << dAdn << "  dBdn = " << dBdn <<"\n";
      Xyce::dout() << "  dAdp = " << dAdp << "  dBdp = " << dBdp <<"\n";
    }
#endif

    dRdnVec[i] = dRsrhdn + dRaugdn;
    dRdpVec[i] = dRsrhdp + dRaugdp;

#else

    double n  = nnVec[i];
    double p  = npVec[i];
    double tn = tnVec[i];
    double tp = tpVec[i];

    dRsrhdn = matSupport.pdRsrhN(bulkMaterial,Ni,n,p,tn,tp);
    dRsrhdp = matSupport.pdRsrhP(bulkMaterial,Ni,n,p,tn,tp);

    dRaugdn = matSupport.pdRaugN(bulkMaterial,Ni,n,p);
    dRaugdp = matSupport.pdRaugP(bulkMaterial,Ni,n,p);

#endif

    dRdnVec[i] = (dRsrhdn + dRaugdn);
    dRdpVec[i] = (dRsrhdp + dRaugdp);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  dRdnVec["<<i<<"] = " << dRdnVec[i];
      Xyce::dout() << "  dRdpVec["<<i<<"] = " << dRdpVec[i];
      Xyce::dout() << "\n";
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcElectronCurrent
// Purpose       :
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcElectronCurrent ()
{
  bool bsuccess = true;

  int i;
  Ut = Vt/scalingVars.V0;

  double jnMax = 0.0;

#ifdef Xyce_DEBUG_DEVICE
  int iMaxIndex = 0;
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcElectronCurrent\n";
    Xyce::dout() << "\n";
  }
#endif


  for (i=0; i<numMeshEdges; ++i)
  {

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "\n";
      Xyce::dout() << "i="<<i<<"\n";
    }
#endif

    mEdge * edgePtr = meshContainerPtr->getEdge(i);

    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;
    double elen   = edgePtr->elen;

    JnVec[i] = Jn(nnVec[inodeA], nnVec[inodeB], EfieldVec[i],
                  unE_Vec[i], elen );

    if (jnMax < fabs(JnVec[i]) )
    {
#ifdef Xyce_DEBUG_DEVICE
      iMaxIndex = i;
#endif
      jnMax = fabs(JnVec[i]);
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  J*scalingVars.J0="<<JnVec[i]*scalingVars.J0; Xyce::dout() << "\n";
    }
#endif

  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    double jScale = 1.0;
    double xScale = 1.0;
    if (variablesScaled) { jScale = scalingVars.J0; xScale = scalingVars.x0; }

    Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "  Max Electron current = " << jnMax;
    Xyce::dout() << "  jScale = " << jScale;
    Xyce::dout() << "\n";

    mEdge * edgePtr = meshContainerPtr->getEdge(iMaxIndex);
    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;

    Xyce::dout() << "  max locations A: (x,y) = (" << xVec[inodeA]*xScale<<", ";
    Xyce::dout() << yVec[inodeA]*xScale<<")\n";
    Xyce::dout() << "  max locations B: (x,y) = (" << xVec[inodeB]*xScale<<", ";
    Xyce::dout() << yVec[inodeB]*xScale<<")\n";
    Xyce::dout() << "  nodes: inodeA = "<<inodeA<<"  inodeB = " << inodeB<< std::endl;

    Xyce::dout() << "  VVec["<<inodeA<<"] = " << VVec[inodeA];
    Xyce::dout() << "  VVec["<<inodeB<<"] = " << VVec[inodeB] << std::endl;

    Xyce::dout() << "  nnVec["<<inodeA<<"] = " << nnVec[inodeA];
    Xyce::dout() << "  nnVec["<<inodeB<<"] = " << nnVec[inodeB] << std::endl;

    Xyce::dout() << "  elen = " << edgePtr->elen << std::endl;
    Xyce::dout() << section_divider << std::endl;

  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdElectronCurrent
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with electron current.
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::pdElectronCurrent ()
{
  bool bsuccess = true;

  int i;
  Ut = Vt/scalingVars.V0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::pdElectronCurrent\n";
    Xyce::dout() << "\n";
  }
#endif

  for (i=0;i<numMeshEdges;++i)
  {

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << subsection_divider <<"\n";
      Xyce::dout() << "i="<<i;
    }
#endif
    mEdge * edgePtr = meshContainerPtr->getEdge(i);

    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;
    double elen   = edgePtr->elen;

    dJndn1Vec[i] = dJndn1(nnVec[inodeA], nnVec[inodeB], EfieldVec[i],
                          unE_Vec[i], elen);

    dJndn2Vec[i] = dJndn2(nnVec[inodeA], nnVec[inodeB], EfieldVec[i],
                          unE_Vec[i], elen);

    dJndV1Vec[i] = dJndV1(nnVec[inodeA], nnVec[inodeB], EfieldVec[i],
                          unE_Vec[i], elen);

    dJndV2Vec[i] = dJndV2(nnVec[inodeA], nnVec[inodeB], EfieldVec[i],
                          unE_Vec[i], elen);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " dJndn1="<<dJndn1Vec[i];
      Xyce::dout() << " dJndn2="<<dJndn2Vec[i];
      Xyce::dout() << " dJndV1="<<dJndV1Vec[i];
      Xyce::dout() << " dJndV2="<<dJndV2Vec[i];
      Xyce::dout() << "\n";
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcHoleCurrent
// Purpose       : This function assumes scaling is on.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcHoleCurrent ()
{
  bool bsuccess = true;
  int i;
  Ut = Vt/scalingVars.V0;

  double jpMax = 0.0;

#ifdef Xyce_DEBUG_DEVICE
  int iMaxIndex = 0;
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcHoleCurrent\n";
    Xyce::dout() << "\n";
  }
#endif

  for (i=0;i<numMeshEdges;++i)
  {

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "\n";
      Xyce::dout() << "i="<<i<<"\n";
    }
#endif

    mEdge * edgePtr = meshContainerPtr->getEdge(i);

    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;
    double elen   = edgePtr->elen;

    JpVec[i] = Jp(npVec[inodeA], npVec[inodeB],
                  EfieldVec[i], upE_Vec[i], elen);

    if (jpMax < fabs(JpVec[i]) )
    {
#ifdef Xyce_DEBUG_DEVICE
      iMaxIndex = i;
#endif
      jpMax = fabs(JpVec[i]);
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  J*scalingVars.J0="<<JpVec[i]*scalingVars.J0; Xyce::dout() <<"\n";
    }
#endif

  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    double jScale = 1.0;
    double xScale = 1.0;
    if (variablesScaled) { jScale = scalingVars.J0; xScale = scalingVars.x0; }

    Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "  Max Hole current = " << jpMax;
    Xyce::dout() << "  jScale = " << jScale;
    Xyce::dout() << "\n";

    mEdge * edgePtr = meshContainerPtr->getEdge(iMaxIndex);
    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;

    Xyce::dout() << "  max locations A: (x,y) = (" << xVec[inodeA]*xScale<<", ";
    Xyce::dout() << yVec[inodeA]*xScale<<")\n";
    Xyce::dout() << "  max locations B: (x,y) = (" << xVec[inodeB]*xScale<<", ";
    Xyce::dout() << yVec[inodeB]*xScale<<")\n";
    Xyce::dout() << "  nodes: inodeA = "<<inodeA<<"  inodeB = " << inodeB<< std::endl;

    Xyce::dout() << "  VVec["<<inodeA<<"] = " << VVec[inodeA];
    Xyce::dout() << "  VVec["<<inodeB<<"] = " << VVec[inodeB] << std::endl;

    Xyce::dout() << "  npVec["<<inodeA<<"] = " << npVec[inodeA];
    Xyce::dout() << "  npVec["<<inodeB<<"] = " << npVec[inodeB] << std::endl;

    Xyce::dout() << "  elen = " << edgePtr->elen << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdHoleCurrent
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with the hole current.
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::pdHoleCurrent ()
{
  bool bsuccess = true;

  int i;
  Ut = Vt/scalingVars.V0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::pdHoleCurrent\n";
    Xyce::dout() << "\n";
  }
#endif


  for (i=0; i<numMeshEdges; ++i)
  {

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << subsection_divider <<"\n";
      Xyce::dout() << "i="<<i;
      Xyce::dout() << "\n";
    }
#endif
    mEdge * edgePtr = meshContainerPtr->getEdge(i);

    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;
    double elen   = edgePtr->elen;

    dJpdn1Vec[i] = dJpdn1(npVec[inodeA], npVec[inodeB],
                          EfieldVec[i], upE_Vec[i], elen );

    dJpdn2Vec[i] = dJpdn2(npVec[inodeA], npVec[inodeB],
                          EfieldVec[i], upE_Vec[i], elen );

    dJpdV1Vec[i] = dJpdV1(npVec[inodeA], npVec[inodeB],
                          EfieldVec[i], upE_Vec[i], elen );

    dJpdV2Vec[i] = dJpdV2(npVec[inodeA], npVec[inodeB],
                          EfieldVec[i], upE_Vec[i], elen );

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " dJpdn1="<<dJpdn1Vec[i];
      Xyce::dout() << " dJpdn2="<<dJpdn2Vec[i];
      Xyce::dout() << " dJpdV1="<<dJpdV1Vec[i];
      Xyce::dout() << " dJpdV2="<<dJpdV2Vec[i];
      Xyce::dout() << "\n" << "\n";
    }
#endif

  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcEfield
// Purpose       : This function works with or without scaling.
//
// Special Notes : For all "edge" defined vector variables, the gradients
//                  are based on (nodeB-nodeA).
//
//                When used in the the box integration algorithm, the
//                "local" node is always node A, and the neighbor node is
//                node B.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
bool Instance::calcEfield ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  int iMaxIndex = 0;

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::calcEfield\n";
    Xyce::dout() << "\n";
  }
#endif

  double absEfield;
  Emax = 0.0;

  for (i=0;i<numMeshEdges;++i)
  {
    mEdge * edgePtr = meshContainerPtr->getEdge(i);

    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;
    double elen   = edgePtr->elen;

    EfieldVec[i] = -(VVec[inodeB] - VVec[inodeA])/elen;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  VV[A] = " << VVec[inodeA];
      Xyce::dout() << "  VV[B] = " << VVec[inodeB];
      Xyce::dout() << "  elen  = " << elen;

      Xyce::dout() << "  EfieldVec["<<i<<"] = " << EfieldVec[i];
      Xyce::dout() << "  E*scalingVars.E0 = " << EfieldVec[i]*scalingVars.E0;
      Xyce::dout() << "\n";
    }
#endif
    if (elen <= 0.00)
    {
      Xyce::dout() << "  edge = " << i;
      Xyce::dout() << "  elen = " << elen;
      Xyce::dout() << "  elen less than zero.  Exiting" << std::endl;
      exit(0);
    }

    absEfield = fabs(EfieldVec[i]);
    if (absEfield > Emax)
    {
      Emax = absEfield;
#ifdef Xyce_DEBUG_DEVICE
      iMaxIndex = i;
#endif
    }
  }

  double eScale;
  if (variablesScaled) { eScale = scalingVars.E0;  }
  else                 { eScale = 1.0; }

  Emax *= eScale;

#ifdef Xyce_DEBUG_DEVICE
  double xScale;
  if (variablesScaled) { xScale = scalingVars.x0;  }
  else                 { xScale = 1.0; }

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "  Max Efield = " << Emax;
    Xyce::dout() << "  eScale     = " << eScale;
    Xyce::dout() << "\n";

    mEdge * edgePtr = meshContainerPtr->getEdge(iMaxIndex);
    int    inodeA = edgePtr->inodeA;
    int    inodeB = edgePtr->inodeB;

    Xyce::dout() << "  max locations A: (x,y) = (" << xVec[inodeA]*xScale<<", ";
    Xyce::dout() << yVec[inodeA]*xScale<<")\n";
    Xyce::dout() << "  max locations B: (x,y) = (" << xVec[inodeB]*xScale<<", ";
    Xyce::dout() << yVec[inodeB]*xScale<<")\n";
    Xyce::dout() << "  nodes: inodeA = "<<inodeA<<"  inodeB = " << inodeB<< std::endl;
    Xyce::dout() << "  VVec["<<inodeA<<"] = " << VVec[inodeA];
    Xyce::dout() << "  VVec["<<inodeB<<"] = " << VVec[inodeB] << std::endl;
    Xyce::dout() << "  elen = " << edgePtr->elen << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::enablePDEContinuation
// Purpose       : Sets up the various parameters neccessary for a continuation
//                 calculation.  Mainly, it sets up the voltage step size
//                 for all the voltage BC's.
//
// Special Notes : This function is called before a Newton loop, or set of
//                 Newton loops is called.  Thus, to have the correct Vckt,
//                 it is neccessary obtain it from the solution vector
//                 directly.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool Instance::enablePDEContinuation ()
{
  bool bnoChange = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -5 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "Instance::enableContinuation.  " << outputName;
    Xyce::dout() << std::endl;
  }
#endif

  continuationAlpha = 1.0;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  // Save the old value, if this function has never been called before.
  // The old value really needs to represent what Vckt was the last time
  // the PDE problem was solved.  The external circuit, which supplies
  // Vckt, may have changed a great deal in between PDE solves, when
  // running 2level Newton.  For that reason, Vckt_old is saved
  // at the end of the previous continuation, if there was one.
  if (!enableContinuationCalled)
  {
    for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
    {
      iterDI->Vckt_old = iterDI->Vckt;
    }
  }

  obtainNodeVoltages ();

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    iterDI->Vckt_final = iterDI->Vckt;
    iterDI->Vckt_orig  = iterDI->Vckt;
  }

  // This (voltlim) is a very new thing.  Use carefully...
  if (getDeviceOptions().voltageLimiterFlag && voltLimFlag)
  {
    applyVoltageLimiting ();
  }

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    double dV,tmp1V, tmp2V;
    tmp1V = iterDI->Vckt_final;
    tmp2V = iterDI->Vckt_old;
    dV    = tmp1V - tmp2V;

    iterDI->Vckt_delta = dV;

    iterDI->Vckt_deltaC = dV/
                          (static_cast<double>(getSolverState().maxPDEContinuationSteps));

    // if this deltaC is too big, then we need to change the
    // number of continuation steps.
    double maxDelta = maxVoltDelta;

    if (fabs(iterDI->Vckt_deltaC) > maxDelta)
    {
      int tmp_steps = static_cast<int>(fabs(dV)/maxDelta) + 1;
      getSolverState().maxPDEContinuationSteps = tmp_steps;

      iterDI->Vckt_deltaC = dV/
                            static_cast<double>(getSolverState().maxPDEContinuationSteps);
    }

    if (fabs(dV) > 1.0e-3) bnoChange = false;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > -5 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << std::endl;
      Xyce::dout() << outputName << " ";
      Xyce::dout().width(10);
      Xyce::dout() << iterDI->eName;
      Xyce::dout().width(10); Xyce::dout().precision(2);
      Xyce::dout() << ":  dV = " << dV;
      Xyce::dout() << "  Vckt_final = " << iterDI->Vckt_final;
      Xyce::dout() << "  Vckt_old   = " << iterDI->Vckt_old << std::endl;
      Xyce::dout() << "  delta      = " << iterDI->Vckt_delta;
      Xyce::dout() << "  deltaC     = " << iterDI->Vckt_deltaC;
      Xyce::dout() << "  steps      = " << getSolverState().maxPDEContinuationSteps;
      Xyce::dout() << std::endl;
    }
#endif

    iterDI->Vckt_ramp = iterDI->Vckt_old;
    iterDI->Vckt_ramp_old = iterDI->Vckt_old;
  }

  // now set up any continuation params associated with photogen.
  bool bnoChangePhotogen = true;

  if (photogenOnFlag)
  {
    std::string msg = "Instance::enablePDEContinuation: PhotogenContinuation not supported";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  bnoChange = bnoChange && bnoChangePhotogen;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -5 && getSolverState().debugTimeFlag)
  {
    if (bnoChange) Xyce::dout() << "bnoChange is TRUE" << std::endl;
    else           Xyce::dout() << "bnoChange is FALSE" << std::endl;
  }
#endif

  if (!enableContinuationCalled) enableContinuationCalled = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -5 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif

  // if none of the boundary conditions  have changed, then
  // return a false.
  return (!bnoChange);
}

//-----------------------------------------------------------------------------
// Function      : Instance::disablePDEContinuation
//
// Purpose       : This function mostly sets up the "old" values of Vckt
//                 and A1, so that they are correct for the next time the
//                 continuation loop is enabled.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool Instance::disablePDEContinuation ()
{
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    iterDI->Vckt_old   = iterDI->Vckt_final;
  }

  photoA1_old = photoA1_final;

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
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    iterDI->Vckt_ramp = iterDI->Vckt_old + (iterDI->Vckt_delta)*alpha;

    // make sure we haven't gone too far:
    if ((iterDI->Vckt_delta >  0 && iterDI->Vckt_ramp >  iterDI->Vckt_final) ||
        (iterDI->Vckt_delta <= 0 && iterDI->Vckt_ramp <= iterDI->Vckt_final) )
    {
      iterDI->Vckt_ramp = iterDI->Vckt_final;

#ifdef Xyce_NEW_PDE_CONTINUATION
      // this line allows us to remove "disablePDEContinuation".
      iterDI->Vckt_old   = iterDI->Vckt_final;
#endif
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      Xyce::dout() << outputName << " ";
      Xyce::dout().width(10);
      Xyce::dout() <<iterDI->eName;
      Xyce::dout() << "\tVckt_ramp = " << iterDI->Vckt_ramp;
      Xyce::dout() << "\tVckt_old = " << iterDI->Vckt_old;
      Xyce::dout() << "\talpha = " << alpha;
      Xyce::dout() << std::endl;
    }
#endif
  }

  // now do the photogeneration term, if neccessary:
  photoA1_ramp =  photoA1_old + photoA1_Delta * alpha;

  // make sure we haven't gone too far:
  if ((photoA1_Delta >  0 && photoA1_ramp >  photoA1_final) ||
      (photoA1_Delta <= 0 && photoA1_ramp <= photoA1_final) )
  {
    photoA1_ramp = photoA1_final;

#ifdef Xyce_NEW_PDE_CONTINUATION
    // this line enables removal of function disablePDEContinuation.
    photoA1_old = photoA1_final;
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  Xyce::dout() << "  photoA1_ramp = " << photoA1_ramp << std::endl;
  Xyce::dout() << section_divider << std::endl;
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::setPDEContinuationBeta
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/18/03
//-----------------------------------------------------------------------------
void Instance::setPDEContinuationBeta (double beta)
{
#ifdef Xyce_DEBUG_DEVICE
  Xyce::dout() << section_divider << std::endl;
  Xyce::dout() << "Instance::setPDEContinuationBeta" << std::endl;
#endif

  // First do the voltage boundary conditions:
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    iterDI->Vckt_ramp = iterDI->Vckt*beta;
#ifdef Xyce_DEBUG_DEVICE
    Xyce::dout() << "  " << iterDI->eName << "  Vckt_ramp = " << iterDI->Vckt_ramp << std::endl;
#endif
  }
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  Device *device = new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);

  return device;
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("pde", 2)
    .registerModelType("pde", 2)
    .registerModelType("zod", 2);
}

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce
