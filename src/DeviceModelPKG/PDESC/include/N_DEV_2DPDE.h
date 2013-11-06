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
// Revision Number: $Revision: 1.70.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:35 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_2DPDE_h
#define Xyce_N_DEV_2DPDE_h

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DevicePDEInstance.h>
#include <N_DEV_DevicePDEModel.h>

#include <N_DEV_DeviceInterfaceNode.h>
#include <N_DEV_PDE_2DMesh.h>

#include <N_UTL_BreakPoint.h>


// Have the "new" boundary conditions turned on by default.
#if 1
#define Xyce_NEW_BC 1
#endif

namespace Xyce {
namespace Device {
namespace TwoDPDE {

class Model;

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

  // functions
public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Instance(InstanceBlock & IB,
                    Model & Miter,
                    MatrixLoadData & mlData1,
                    SolverState &ss1,
                    ExternData  &ed1,
                    DeviceOptions & do1);
  ~Instance();

private:
  Instance(const Instance &right);
  Instance &operator=(const Instance &right);

public:
  void registerGIDs (const list<index_pair> & intGIDListRef,
                     const list<index_pair> & extGIDListRef );

  void setupIntNameMap  ();
  void setupRowColPairs ();

  void registerStateGIDs (const list<index_pair> & staGIDListRef);

  void registerLIDs( const vector<int> & intLIDVecRef,
                     const vector<int> & extLIDVecRef );
  void registerStateLIDs( const vector<int> & staLIDVecRef );
  map<int,string> & getIntNameMap ();

  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

  bool processParams (string param = "");

  bool processOneTimeParams( Param & ndParam );

  bool processDopingParams (Param & ndParam, string param);

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

  CompositeParam *constructComposite (string & compositeName, string & paramName);

public:
  // iterator reference to the resistor model which owns this instance.
  // Getters and setters
  Model &getModel() {
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
  vector<DeviceInterfaceNode> dIVec;

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

  string meshFileName;
  string deviceType;
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
  vector<double> xVec;   // x locations
  vector<double> yVec;   // y locations
  vector<double> CVec;   // doping
  vector<double> CdonorVec;   // doping
  vector<double> CacceptorVec;   // doping

  vector<double> minDXVec; // minimum mesh spacing connected to this node.

  vector<double> areaVec;

  vector<double> VVec;   // electrostatic potential
  vector<double> nnVec;  // electron density
  vector<double> npVec;  // hole density

  vector<double> totSrcVec; // total source term.
  vector<double> RVec;      // recombination.
  vector<double> SVec;      // radiation source term.

  vector<double> elecPenalty;  // penalty term for negative e- density.
  vector<double> holePenalty;  // penalty term for negative h+ density.
  vector<double> pdElecPenalty; // penalty deriv. term for negative e- density.
  vector<double> pdHolePenalty; // penalty deriv. term for negative h+ density.

  vector<double> unVec;  // spatially dependent mobility, electron
  vector<double> upVec;  // spatially dependent mobility, hole
  vector<double> unE_Vec; // mobility along edge, electron
  vector<double> upE_Vec; // mobility along edge, hole
  vector<double> tnVec;  // spatially dependent lifetimes, electron
  vector<double> tpVec;  // spatially dependent lifetimes, hole

  vector<double> EfieldVec; // electric field along an edge.

  vector<double> JnVec; // electron current density, along an edge
  vector<double> JpVec; // hole current density, along an edge

  vector<double> displPotential;  // time derivative of potential at a node,
  // used in calculating displacement current.

  vector<double> displCurrent;  // displacement current along an edge.

  vector<double> outputVec;

  // derivative arrays:

  // derivatives of recombination terms:
  vector<double> dRdpVec;  // derivative of R w.r.t. np  (at a node)
  vector<double> dRdnVec;  // derivative of R w.r.t. nn  (at a node)

  // derivatives of current density terms:
  vector<double> dJndn1Vec;
  vector<double> dJndn2Vec;
  vector<double> dJndV1Vec;
  vector<double> dJndV2Vec;

  vector<double> dJpdn1Vec;
  vector<double> dJpdn2Vec;
  vector<double> dJpdV1Vec;
  vector<double> dJpdV2Vec;

  // matrix index arrays:
  // external variable information is in the dIVec data structure.

  // boundary "neighbor" stencil
  // This array is set to 1 if we are on a mesh node which has a
  // nearest neighbor node which is a boundary condition node.
  // 0 if not.
  vector<int> boundarySten;
  vector<int> boundaryStenV;
  vector<int> boundaryStenN;
  vector<int> boundaryStenP;

  vector<int> boundaryTest;

  // boundary neighbor stencil
  // This array is set to 1 if we are next to a boundary,  but not on
  // a boundary.  0 if not.  A node directly on the boundary
  // does should return a zero.
  vector<int> boundaryNeighborSten;

  // internal variable index arrays:

  // index arrays for poisson's equation:
  vector<int>           Vrowarray;
  vector< vector<int> > Vcolarray;

  // index arrays for electron continuity:
  vector<int>           Nrowarray;
  vector< vector<int> > Ncolarray;

  // index arrays for hole continuity:
  vector<int>           Prowarray;
  vector< vector<int> > Pcolarray;

  // variable ownership arrays.
  vector<int>  vOwnVec; // on processor tag, for voltage variables.
  vector<int> nnOwnVec; // on processor tag, for elec. density variables.
  vector<int> npOwnVec; // on processor tag, for hole density variables.

  //local id's (offsets)
  vector<int>     li_Vrowarray;
  vector<int>     li_Nrowarray;
  vector<int>     li_Prowarray;

  vector< vector<int> > li_VoffsetArray;
  vector< vector<int> > li_NoffsetArray;
  vector< vector<int> > li_PoffsetArray;

  vector<int> MESHtoLID_V;
  vector<int> MESHtoLID_N;
  vector<int> MESHtoLID_P;


  // minor arrays used in tecplot output:
  vector<UINT> aiEdge;
  vector<UINT> aiEdge_nf;
  //vector<UINT> electrodeEdge;
  UINT iNumPlotEdges;
  UINT iNumPlotEdges_nf;

  // this map is mostly used for processing the netlist information
  // regarding boundary conditions.
  map<string,string> tmpBCmap;

  // label index array (of numMeshPoints length)
  vector<int> labelIndex;
  vector<string> labelNameVector;
  map<string,int> labelDIMap;

  // map between a mesh point index and a list of nearest neighbors
  // for that mesh point.
  multimap < int, int* > meshNeighborMultiMap;

  // vector of electrode data:
  map<string,PDE_2DElectrode*> electrodeMap;

  // displacement current state variable information:
  vector<int> stateDispl;
  vector<int> stateDispl_owned;  // this is an int array because
  // bool arrays are problematic,
  // esp. with some STL implementations.

  //local id's (offsets)
  vector<int> li_stateDispl;

  int numMeshPoints;
  int numInterfaceMeshPoints; // number of mesh points that
  // are along electrode boundaries
  int numMeshEdges;
  int numMeshCells;
  int numMeshLabels;
  int maxColsPerRow;
  int numElectrodes;

  // 2d array of conductances, for 2-level "ckt phase" loads.
  vector< vector<double> > condVec;

  // 2d array of capacitances,
  vector< vector<double> > capVec;

  bool  pdTermsAllocated;

  // data related to DMA matrix loads.
  vector<int> meshToLID;
  vector< vector<int> > jacStamp;
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

public:
  static ParametricData<Model> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Model    (const ModelBlock & MB,
                     SolverState & ss1,
                     DeviceOptions & do1);
  ~Model   ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual std::ostream &printOutInstances(std::ostream &os) const;

  bool processParams (string param = "");
  bool processInstanceParams (string param = "");


public:
  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

private:
  vector<Instance*> instanceContainer;

private:
};

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::TwoDPDE::Instance N_DEV_2DPDEInstance;
typedef Xyce::Device::TwoDPDE::Model N_DEV_2DPDEModel;

#endif
