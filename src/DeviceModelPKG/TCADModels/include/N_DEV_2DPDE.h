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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_2DPDE.h,v $
//
// Purpose        : This file contains the classes neccessary for a 2D PDE
//                  based simulation.  MOSFETs, BJTs, Diodes, etc.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 11/14/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.86.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:31 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_2DPDE_h
#define Xyce_N_DEV_2DPDE_h

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DevicePDEInstance.h>
#include <N_DEV_DevicePDEModel.h>

#include <N_DEV_DeviceInterfaceNode.h>
#include <N_DEV_PDE_2DMesh.h>

#include <N_DEV_DiodePDE.h>

#include <N_UTL_BreakPoint.h>


// Have the "new" boundary conditions turned on by default.
#if 1
#define Xyce_NEW_BC 1
#endif

namespace Xyce {
namespace Device {
namespace TwoDPDE {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, DiodePDE::Traits>
{
  static const char *name() {return "2D PDE Device";}
  static const char *deviceTypeName () {return "PDE level 2";}
  static const int numNodes() {return 2;}
  static const int numOptionalNodes() {return 100;}
  static const int numFillNodes() {return 2;}
  static const bool modelRequired() {return true;}
  static const bool isPDEDevice() {return true;}
  static const bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : Instance class for .
//
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
class Instance : public DevicePDEInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;
    
  // functions
public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &             IB,
     Model &                     model,
     const FactoryBlock &factory_block);
  ~Instance();

private:
  Instance(const Instance &right);
  Instance &operator=(const Instance &right);

public:
  void registerGIDs (const std::list<index_pair> & intGIDListRef,
                     const std::list<index_pair> & extGIDListRef );

  void setupIntNameMap  ();
  void setupRowColPairs ();

  void registerStateGIDs (const std::list<index_pair> & staGIDListRef);

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  std::map<int,std::string> & getIntNameMap ();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();

  bool processOneTimeParams( Param & ndParam );

  bool processDopingParams (Param & ndParam, std::string param);

  bool processElectrodeParams (Param & ndParam);

  bool setupJacStamp ();

  bool doSensMeshResize ();
  bool undoSensMeshResize ();

  bool setupMesh ();
  bool doAllocations ();

  bool setupDINodes ();
  bool setupBCEdgeAreas ();

  bool setupBoundaryStencil ();
  bool setupNumVars ();

  bool checkForElectrodeOverlap ();

  bool setupLabelIndex ();
  bool setupMinDXVector ();

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  // functions that are used by both new-DAE and old-DAE:
  bool loadVecNLPoisson (double scalar, N_LAS_Vector * vecPtr);
  bool loadMatNLPoisson (N_LAS_Matrix * matPtr);
  bool loadMatKCLDDForm (N_LAS_Matrix * matPtr);
  bool loadMatDDForm (double dndtScalar, N_LAS_Matrix * matPtr);
  bool loadVecDDForm (double scalar,double dndtScalar,N_LAS_Vector *vecPtr);
  bool loadMatCktTrivial (N_LAS_Matrix * matPtr);
  // end of the "both" DAE functions.

  bool setInitialGuess ();
  bool loadRHSNonlinPoisson ();
  bool loadRHSDDFormulation ();
  bool loadRHSExtractedConductance ();

  bool plotfileFlag () {return true;}

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEQDDFormulation ();
  bool loadDAEQExtractedConductance ();

  bool loadDAEFVector ();
  bool loadDAEFNonlinPoisson ();
  bool loadDAEFDDFormulation ();
  bool loadDAEFExtractedConductance ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();

  bool loadDAEdQdxDDFormulation ();
  bool loadDAEdQdxExtractedConductance ();

  bool loadDAEdFdx ();
  bool loadDAEdFdxNonlinPoisson ();
  bool loadDAEdFdxDDFormulation ();
  bool loadDAEdFdxExtractedConductance ();

  bool calcLifetimes ();
  bool calcMobilities ();
  bool updateTemperature(const double & temp_tmp);
  bool calcVoltDepDensities ();

  bool calcDopingProfile ();
  bool calcInitialGuess ();
  bool obtainSolution ();
  bool obtainNodeVoltages ();
  bool applyVoltageLimiting ();
  bool calcVequBCs ();
  bool calcDensityBCs ();
  bool calcBoundaryConditions ();

  bool setupMiscConstants ();
  bool setupScalingVars ();
  bool scaleVariables ();
  bool unScaleVariables ();
  bool scaleDopeVariables ();
  bool unScaleDopeVariables ();

  bool calcRecombination ();


  bool setupPhotogen ();
  bool calcPhotogen ();
  bool calcPenalty  ();
  bool enablePhotogenContinuation ();
  bool sumSources ();

  bool calcElectronCurrent ();
  bool calcHoleCurrent ();
  bool calcEfield ();

  bool calcTerminalCharges ();
  bool calcTerminalCurrents ();
  bool calcConductance (int iElectrode, const N_LAS_Vector * dxdvPtr);
  bool calcDXDV ();
  bool loadDFDV (int ielectrode, N_LAS_Vector * dfdvPtr);

  bool pdRecombination ();
  bool pdElectronCurrent ();
  bool pdHoleCurrent ();
  bool pdTerminalCurrents ();
  bool pdTerminalCharges ();
  bool allocatePDTerms ();

  bool pdPenalty  ();

  bool outputTecplot        ();
  bool outputTecplotVectors ();
  bool tecplotGeomOutput  (FILE  *fp1);
  bool outputSgplot ();
  bool outputGnuplot ();
  bool outputTxtData ();

  bool enablePDEContinuation ();
  bool disablePDEContinuation ();
  void setPDEContinuationAlpha (double alpha);
  void setPDEContinuationBeta  (double beta);

  bool outputPlotFiles ();

  CompositeParam *constructComposite (const std::string & compositeName, const std::string & paramName);

public:
  // iterator reference to the resistor model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

protected:
private:

  // physical constants:
  double Is;     // saturation current
  double Id;     // diode current, analytic form

  double Emax;   // maximum electric field (V/cm)

  double VminExp;   // maximum potential (V), used in exponential expressions
  double VmaxExp;   // minimum potential (V), used in exponential expressions

  // device interface node vector:
  std::vector<DeviceInterfaceNode> dIVec;

  bool penaltyPrefacGiven;
  bool penaltyPowGiven;

  double LeadCurrent1;
  double LeadCurrent2;
  double LeadCurrent3;
  double LeadCurrent4;
  double LeadCurrent5;
  double LeadCurrent6;
  double LeadCurrent7;
  double LeadCurrent8;

  // doping profile constants:
  double Na;     // acceptor concentration on p-side (cm^-3)
  double Nd;     // donor concentration on n-side    (cm^-3)
  double Vbi;    // built-in potential               (V)
  double WJ;     // linearly graded junction width   (cm)
  double XC;     // center of graded junction        (cm)
  double XL;     // start of graded junction         (cm)
  double XR;     // end of graded junction           (cm)

  // boundary condition variables:
  double NnMax;  // maximum electron concentration
  double NpMax;  // maximum hole concentration.
  double NnMin;  // maximum electron concentration
  double NpMin;  // maximum hole concentration.

  // option to use the old intrinsic calculation, to maintain backward compatibility.
  bool useOldNi;
  bool useOldNiGiven;

  std::string meshFileName;
  std::string deviceType;
  bool usingInternalMesh;
  bool deviceInitialized;
  bool meshPerturbed;
  bool dopingPerturbed;
  bool photogenPerturbed;


  // meshing variables, if using internally generated mesh.
  int numMeshPointsX;
  int numMeshPointsY;
  double deviceLength;
  double deviceWidth;

  bool cylGeomFlag;

  bool penaltyFlag;
  double penaltyPrefac;
  double penaltyPow;
  double PulseData;


  // diode cross sectional area:
  double area;

#ifdef Xyce_OXIDE_ENABLED
  bool allOxideFlag;
#endif

  bool gradedJunctionFlag;
  bool calledBeforeSIGB;
  int  callsOSG;
  int  callsOTEC;
  int  callsOTECvec;
  int  callsOGNU;
  int  callsOTXT;

  bool displCurrentFlag;
  bool constBoundaryFlag;

  bool calcConductanceFlag; // Set when calcConductance is
  // called for the 1st time.

  int equationSet;

  double outputInterval;
  int  outputIndex;
  bool outputNLPoisson;
  double lastOutputTime;

  int tecplotLevel;
  int sgplotLevel;
  int gnuplotLevel;
  int txtDataLevel;

  int interpGridSize;

  bool voltLimFlag;

  bool useMatrixGIDFlag;
  bool useVectorGIDFlag;

  // mesh container pointer:
  PDE_2DMesh * meshContainerPtr;
  PDE_2DMesh * meshCopyContainerPtr;

  // array pointers:
  std::vector<double> xVec;   // x locations
  std::vector<double> yVec;   // y locations
  std::vector<double> CVec;   // doping
  std::vector<double> CdonorVec;   // doping
  std::vector<double> CacceptorVec;   // doping

  std::vector<double> minDXVec; // minimum mesh spacing connected to this node.

  std::vector<double> areaVec;

  std::vector<double> VVec;   // electrostatic potential
  std::vector<double> nnVec;  // electron density
  std::vector<double> npVec;  // hole density

  std::vector<double> totSrcVec; // total source term.
  std::vector<double> RVec;      // recombination.
  std::vector<double> SVec;      // radiation source term.

  std::vector<double> elecPenalty;  // penalty term for negative e- density.
  std::vector<double> holePenalty;  // penalty term for negative h+ density.
  std::vector<double> pdElecPenalty; // penalty deriv. term for negative e- density.
  std::vector<double> pdHolePenalty; // penalty deriv. term for negative h+ density.

  std::vector<double> unVec;  // spatially dependent mobility, electron
  std::vector<double> upVec;  // spatially dependent mobility, hole
  std::vector<double> unE_Vec; // mobility along edge, electron
  std::vector<double> upE_Vec; // mobility along edge, hole
  std::vector<double> tnVec;  // spatially dependent lifetimes, electron
  std::vector<double> tpVec;  // spatially dependent lifetimes, hole

  std::vector<double> EfieldVec; // electric field along an edge.

  std::vector<double> JnVec; // electron current density, along an edge
  std::vector<double> JpVec; // hole current density, along an edge

  std::vector<double> displPotential;  // time derivative of potential at a node,
  // used in calculating displacement current.

  std::vector<double> displCurrent;  // displacement current along an edge.

  std::vector<double> outputVec;

  // derivative arrays:

  // derivatives of recombination terms:
  std::vector<double> dRdpVec;  // derivative of R w.r.t. np  (at a node)
  std::vector<double> dRdnVec;  // derivative of R w.r.t. nn  (at a node)

  // derivatives of current density terms:
  std::vector<double> dJndn1Vec;
  std::vector<double> dJndn2Vec;
  std::vector<double> dJndV1Vec;
  std::vector<double> dJndV2Vec;

  std::vector<double> dJpdn1Vec;
  std::vector<double> dJpdn2Vec;
  std::vector<double> dJpdV1Vec;
  std::vector<double> dJpdV2Vec;

  // matrix index arrays:
  // external variable information is in the dIVec data structure.

  // boundary "neighbor" stencil
  // This array is set to 1 if we are on a mesh node which has a
  // nearest neighbor node which is a boundary condition node.
  // 0 if not.
  std::vector<int> boundarySten;
  std::vector<int> boundaryStenV;
  std::vector<int> boundaryStenN;
  std::vector<int> boundaryStenP;

  std::vector<int> boundaryTest;

  // boundary neighbor stencil
  // This array is set to 1 if we are next to a boundary,  but not on
  // a boundary.  0 if not.  A node directly on the boundary
  // does should return a zero.
  std::vector<int> boundaryNeighborSten;

  // internal variable index arrays:

  // index arrays for poisson's equation:
  std::vector<int>           Vrowarray;
  std::vector< std::vector<int> > Vcolarray;

  // index arrays for electron continuity:
  std::vector<int>           Nrowarray;
  std::vector< std::vector<int> > Ncolarray;

  // index arrays for hole continuity:
  std::vector<int>           Prowarray;
  std::vector< std::vector<int> > Pcolarray;

  // variable ownership arrays.
  std::vector<int>  vOwnVec; // on processor tag, for voltage variables.
  std::vector<int> nnOwnVec; // on processor tag, for elec. density variables.
  std::vector<int> npOwnVec; // on processor tag, for hole density variables.

  //local id's (offsets)
  std::vector<int>     li_Vrowarray;
  std::vector<int>     li_Nrowarray;
  std::vector<int>     li_Prowarray;

  std::vector< std::vector<int> > li_VoffsetArray;
  std::vector< std::vector<int> > li_NoffsetArray;
  std::vector< std::vector<int> > li_PoffsetArray;

  std::vector<int> MESHtoLID_V;
  std::vector<int> MESHtoLID_N;
  std::vector<int> MESHtoLID_P;


  // minor arrays used in tecplot output:
  std::vector<UINT> aiEdge;
  std::vector<UINT> aiEdge_nf;
  //std::vector<UINT> electrodeEdge;
  UINT iNumPlotEdges;
  UINT iNumPlotEdges_nf;

  // this map is mostly used for processing the netlist information
  // regarding boundary conditions.
  std::map<std::string,std::string> tmpBCmap;

  // label index array (of numMeshPoints length)
  std::vector<int> labelIndex;
  std::vector<std::string> labelNameVector;
  std::map<std::string,int> labelDIMap;

  // map between a mesh point index and a list of nearest neighbors
  // for that mesh point.
  std::multimap < int, int* > meshNeighborMultiMap;

  // vector of electrode data:
  std::map<std::string,PDE_2DElectrode*> electrodeMap;

  // displacement current state variable information:
  std::vector<int> stateDispl;
  std::vector<int> stateDispl_owned;  // this is an int array because
  // bool arrays are problematic,
  // esp. with some STL implementations.

  //local id's (offsets)
  std::vector<int> li_stateDispl;

  int numMeshPoints;
  int numInterfaceMeshPoints; // number of mesh points that
  // are along electrode boundaries
  int numMeshEdges;
  int numMeshCells;
  int numMeshLabels;
  int maxColsPerRow;
  int numElectrodes;

  // 2d array of conductances, for 2-level "ckt phase" loads.
  std::vector< std::vector<double> > condVec;

  // 2d array of capacitances,
  std::vector< std::vector<double> > capVec;

  bool  pdTermsAllocated;

  // data related to DMA matrix loads.
  std::vector<int> meshToLID;
  std::vector< std::vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
class Model  : public DevicePDEModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Traits;
    
public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &          MB,
     const FactoryBlock &        factory_block);

  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;
    
  virtual std::ostream &printOutInstances(std::ostream &os) const;

  bool processParams ();
  bool processInstanceParams ();


public:
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

  InstanceVector &getInstanceVector() 
  {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const 
  {
    return instanceContainer;
  }

private:
  std::vector<Instance*> instanceContainer;

private:
};

void registerDevice();

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::TwoDPDE::Instance N_DEV_2DPDEInstance;
typedef Xyce::Device::TwoDPDE::Model N_DEV_2DPDEModel;

#endif
