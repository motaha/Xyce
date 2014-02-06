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
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DiodePDE.h,v $
//
// Purpose        : This file contains the classes neccessary for a PDE
//                  based diode simulation.
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
// Revision Number: $Revision: 1.66.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DiodePDE_h
#define Xyce_N_DEV_DiodePDE_h

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DevicePDEInstance.h>
#include <N_DEV_DevicePDEModel.h>
#include <N_DEV_bcData.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_PDE_Electrode.h>

#include <N_DEV_Param.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {
namespace DiodePDE {

class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : Instance class for DiodePDE.
//
// Special Notes :6
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
class Instance : public DevicePDEInstance
{
  friend class Model;
  friend class ParametricData<Instance>;

public:
    static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Instance(
    InstanceBlock &            IB,
    Model & Miter,
    MatrixLoadData &           mlData1,
    SolverState &              ss1,
    ExternData  &              ed1,
    DeviceOptions &            do1);

  Instance(const Instance &right);
  ~Instance();

  CompositeParam *constructComposite (string & ccompositeName, string & paramName);

  map<int,string> & getIntNameMap ();

  void registerLIDs( const vector<int> & intLIDVecRef,
                     const vector<int> & extLIDVecRef );
  void registerStateLIDs( const vector<int> & staLIDVecRef );

  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

  void setupPointers();

  bool processParams (string param = "");

  bool doAllocations ();
  bool setupNodes ();

  bool setupNumVars ();

  bool setupJacStamp ();
  bool cleanupJacStamp ();

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  bool loadDeviceMask ();

  // functions that are used by both new-DAE and old-DAE:
  bool loadVecNLPoisson (double * rhs);
  bool loadMatNLPoisson (N_LAS_Matrix & mat);
  bool loadMatKCLDDForm (N_LAS_Matrix & mat);
  bool loadMatDDForm (N_LAS_Matrix & mat);
  bool loadVecDDForm (double * rhs);
  bool loadMatCktTrivial (N_LAS_Matrix & mat);
  // end of the "both" DAE functions.

  bool setInitialGuess ();
  bool loadRHSNonlinPoisson ();
  bool loadRHSDDFormulation ();
  bool loadRHSExtractedConductance ();

  bool getInstanceBreakPoints( vector<N_UTL_BreakPoint> &breakPointTimes);

  bool plotfileFlag () {return true;}

  bool loadJacNonlinPoisson ();
  bool loadJacKCLDDFormulation ();
  bool loadJacDDFormulation ();
  bool loadJacExtractedConductance ();

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

  bool setEH_inChemistry ();

  bool setupSourceProfile ();

  bool setupDopingProfile ();
  bool calcDopingProfile ();

  bool setupMesh ();
  bool calcInitialGuess ();
  bool obtainSolution ();
  bool obtainNodeVoltages ();
  bool applyVoltageLimiting ();
  bool calcVequBCs ();
  bool calcDensityBCs   ();
  bool calcBoundaryConditions ();
  bool setupMiscConstants ();
  bool setupScalingVars ();
  bool scaleVariables ();
  bool unScaleVariables ();

  bool calcRecombination ();
  bool calcElectronCurrent ();
  bool calcHoleCurrent ();
  bool calcEfield ();

  bool calcTerminalCurrents ();
  bool calcConductance (int iElectrode, const N_LAS_Vector * dxdvPtr);
  bool calcDXDV ();
  bool loadDFDV (int ielectrode, N_LAS_Vector * dfdvPtr);

  bool pdRecombination ();
  bool pdElectronCurrent ();
  bool pdHoleCurrent ();
  bool pdTerminalCurrents ();

  bool outputTecplot ();
  bool outputSgplot  ();

  bool enablePDEContinuation ();
  bool disablePDEContinuation ();
  void setPDEContinuationAlpha (double alpha);

//    bool continuationStatus();
//    void changeContinuationStepSize(double scale);
//    void updateOldContinuationParam();

  bool outputPlotFiles ();

public:
  // iterator reference to the resistor model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  bool indicesSetup_;
  bool includeBaseNode_;
  bool useElectrodeSpec_;
  bool maskVarsTIAFlag_;
  bool scaleDensityToMaxDoping_;
  double densityScalarFraction_;

  bool useVoltageOutputOffset_;
  bool offsetWithFirstElectrode_;
  double VoltageOffset_;

  // physical constants:

  double Emax;     // maximum electric field (V/cm)

  double VminExp;   // maximum potential (V), used in exponential expressions
  double VmaxExp;   // minimum potential (V), used in exponential expressions

  double diodeCap; // estimated diode capacitance (F)

  double LeadCurrent;

  // inputted doping profiles.
  // pdope is often phosphorus
  vector<double> xloc_pdope_vec;
  vector<double> pdope_vec;
  vector<double> y2_pdope_vec;

  // ndope is often boron.
  vector<double> xloc_ndope_vec;
  vector<double> ndope_vec;
  vector<double> y2_ndope_vec;

  // source file arrays
  vector<double> xloc_source_vec;
  vector<double> source_vec;
  vector<double> y2_source_vec;

  // temp spline arrays:
  vector<double> xlocVec;
  vector<double> specVec;
  vector<double> y2Vec;

  // new doping/initial condition storage:
  map <string, vector<double> >  xlocMap;
  map <string, vector<double> >  specMap;

  // vector of electrode data:
  map<string,PDE_1DElectrode*> electrodeMap;

  // boundary condition array.  probably will be of size 2 or 3.
  vector<bcData>  bcVec;

  map<string,int> bcIndexMap;

  // doping profile constants:
  double Na;     // acceptor concentration on p-side (cm^-3)
  double Nd;     // donor concentration on n-side    (cm^-3)
  double Vbi;    // built-in potential               (V)
  double WJ;     // linearly graded junction width   (cm)
  double XC;     // center of graded junction        (cm)
  double XL;     // start of graded junction         (cm)
  double XR;     // end of graded junction           (cm)

  // boundary condition variables:
  // These should eventually replace Nd and Na.
  double NnMax;  // maximum electron concentration
  double NpMax;  // maximum hole concentration.
  double NnMin;  // maximum electron concentration
  double NpMin;  // maximum hole concentration.

  // mesh constants:
  int    NX;     // number of x mesh points
  int    LX;     // index of last x point.

  // continuation parameters:
  double maxVoltDelta;
  bool enableContinuationCalled;

  // option to use the old intrinsic calculation, to maintain backward compatibility.
  bool useOldNi;
  bool useOldNiGiven;

  // some 2D mesh stuff - mostly to catch netlist mistakes (as this is a 1D device model)
  string meshFileName;

  // doping files.
  string dopingFileName;
  string ndopeFileName;
  string pdopeFileName;

  // diode width:
  double width;
  double length;

  double basex; // location of base electrode, if running bjt mode.

  // diode cross sectional area:
  double area;

  // boundary condition variables.  These are superceded by
  // the data in the bcVec and PDE_1DElectrode structures.
  double anodebc;
  double cathodebc;
  double emitterbc;
  double collectorbc;
  double basebc;

  double anodeArea;
  double cathodeArea;

  double emitterArea;
  double collectorArea;
  double baseArea;

  double baseLocation;
  bool baseLocationGiven;

  bool gradedJunctionFlag;
  bool bjtEnableFlag;
  bool calledBeforeUIVB;
  int  callsOTEC;
  int  callsOSG;

  bool displCurrentFlag;

  int equationSet;

  double outputInterval;
  bool outputIntervalGiven;
  int  outputIndex;
  bool outputNLPoisson;
  double lastOutputTime;

  int outputRegion;
  int tecplotLevel;
  int gnuplotLevel;
  int sgplotLevel;

  bool voltLimFlag;
  bool includeAugerRecomb;
  bool includeSRHRecomb;

#ifdef Xyce_DEBUG_DEVICE
  int anodeIndex_user;
  bool anodeIndex_userGiven;
  int cathodeIndex_user;
  bool cathodeIndex_userGiven;
#endif

  int    NUMRC;  // number of row-column pairs.

  vector<double> displCurrent;

  //*************************************
  // Neutron reaction model stuff:
  double junctionArea;

  vector<int> boundarySten;
  vector<int> edgeBoundarySten;
  vector<int> internalBoundarySten;

  vector<int> regBaseIndexVec;
  vector<int> regNumSpecieVec;
  vector<int> regElectronIndexVec;
  vector<int> regHoleIndexVec;

  vector<double> dxVec;  // mesh spacing.
  vector<double> xVec;   // mesh points.
  vector<double> CVec;   // doping
  vector<double> CdonorVec;    // doping
  vector<double> CacceptorVec; // doping
  vector<double> VVec;   // electrostatic potential
  vector<double> ExVec;  // electric field, x-direction.

  vector<double> JnxVec; // electron current density
  vector<double> JpxVec; // hole current density

  vector<double> RVec;   // recombination.
  vector<double> SVec;   // radiation source term.

  vector<double> nnVec;  // electron density
  vector<double> npVec;  // hole density

  //vector<double> unVec;  // spatially dependent mobility, electron
  //vector<double> upVec;  // spatially dependent mobility, hole
#if 0
  vector<double> unE_Vec; // mobility along edge, electron
  vector<double> upE_Vec; // mobility along edge, hole
#else
  vector<pdeFadType> unE_Vec; // mobility along edge, electron
  vector<pdeFadType> upE_Vec; // mobility along edge, hole
#endif

  vector<double> tnVec;  // spatially dependent lifetimes, electron
  vector<double> tpVec;  // spatially dependent lifetimes, hole

  // derivative arrays:

  vector<double> dRdpVec;
  vector<double> dRdnVec;

  vector<double> dJndn1Vec;
  vector<double> dJndn2Vec;
  vector<double> dJndV1Vec;
  vector<double> dJndV2Vec;
  vector<double> dJndp1Vec;
  vector<double> dJndp2Vec;

  vector<double> dJpdn1Vec;
  vector<double> dJpdn2Vec;
  vector<double> dJpdV1Vec;
  vector<double> dJpdV2Vec;
  vector<double> dJpdp1Vec;
  vector<double> dJpdp2Vec;

  // LID indices.
  vector<int> li_Vrowarray;
  vector< vector<int> > li_Vcolarray;

  vector<int> li_Nrowarray;
  vector< vector<int> > li_Ncolarray;

  vector<int> li_Prowarray;
  vector< vector<int> > li_Pcolarray;

  // columns needed for coupledMode==2
  vector< vector<int> > li_N_rxn_colarray;
  vector< vector<int> > li_P_rxn_colarray;

  vector<int> li_stateDispl;

  // matrix pointers
  vector< vector<double *> > fVmatPtr;
  vector< vector<double *> > fNmatPtr;
  vector< vector<double *> > fPmatPtr;
  vector< vector<double *> > qVmatPtr;
  vector< vector<double *> > qNmatPtr;
  vector< vector<double *> > qPmatPtr;

  // map between a mesh point index and a list of nearest neighbors
  // for that mesh point.
  multimap < int, int* > meshNeighborMultiMap;

  // state variable arrays associated with displacement current.
  vector<int> stateDispl;
  vector<int> stateDispl_owned;
  // this is an int array because I think bool arrays
  // might be problematic for some compilers.

  int numMeshPoints;  // this is probably redundant with NX.
  int maxColsPerRow;
  int numElectrodes;

  // 2d array of conductances, for 2-level "ckt phase" loads.
  vector< vector<double> > condVec;

  // data related to DMA matrix loads.
  vector<int> meshToLID;
  vector< vector<int> > jacStamp;
  vector<int> jacMap;
  vector< vector<int> > jacMap2;

  // flags related to chemistry.
  bool dirichletBCFlag;
  bool columnReorderingFlag;

  ScalingVars unscaled_ScalingVars;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
class Model  : public DevicePDEModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class Instance;
  friend class ParametricData<Model>;

public:
    static ParametricData<Model> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Model(
    const ModelBlock & MB,
    SolverState &      ss1,
    DeviceOptions &    do1);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual std::ostream &printOutInstances(std::ostream &os) const;

  bool processParams (string param = "");
  bool processInstanceParams (string param = "");

private:


public:
  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

private:
  vector<Instance*> instanceContainer;
};

} // namespace DiodePDE
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DiodePDE::Instance N_DEV_DiodePDEInstance;
typedef Xyce::Device::DiodePDE::Model N_DEV_DiodePDEModel;

#endif
