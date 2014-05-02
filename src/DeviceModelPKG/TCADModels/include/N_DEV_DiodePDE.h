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
// Revision Number: $Revision: 1.80.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:31 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DiodePDE_h
#define Xyce_N_DEV_DiodePDE_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
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
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "1D PDE Device";}
  static const char *deviceTypeName() {return "PDE level 3";}
  static const int numNodes() {return 2;}
  static const int numOptionalNodes() {return 100;}
  static const bool modelRequired() {return true;}
  static const bool isPDEDevice() {return true;}
  static const bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

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
  friend class Traits;

public:
  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       IB,
     Model &                     model,
     const FactoryBlock &        factory_block);

  Instance(const Instance &right);
  ~Instance();

  CompositeParam *constructComposite (const std::string & ccompositeName, const std::string & paramName);

  std::map<int,std::string> & getIntNameMap ();

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  void setupPointers();

  bool processParams ();

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

  bool getInstanceBreakPoints( std::vector<N_UTL_BreakPoint> &breakPointTimes);

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
  Model &getModel() 
  {
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
  std::vector<double> xloc_pdope_vec;
  std::vector<double> pdope_vec;
  std::vector<double> y2_pdope_vec;

  // ndope is often boron.
  std::vector<double> xloc_ndope_vec;
  std::vector<double> ndope_vec;
  std::vector<double> y2_ndope_vec;

  // source file arrays
  std::vector<double> xloc_source_vec;
  std::vector<double> source_vec;
  std::vector<double> y2_source_vec;

  // temp spline arrays:
  std::vector<double> xlocVec;
  std::vector<double> specVec;
  std::vector<double> y2Vec;

  // new doping/initial condition storage:
  std::map<std::string, std::vector<double> >  xlocMap;
  std::map<std::string, std::vector<double> >  specMap;

  // vector of electrode data:
  std::map<std::string, PDE_1DElectrode*> electrodeMap;

  // boundary condition array.  probably will be of size 2 or 3.
  std::vector<bcData>  bcVec;

  std::map<std::string,int> bcIndexMap;

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
  std::string meshFileName;

  // doping files.
  std::string dopingFileName;
  std::string ndopeFileName;
  std::string pdopeFileName;

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

  std::vector<double> displCurrent;

  //*************************************
  // Neutron reaction model stuff:
  double junctionArea;

  std::vector<int> boundarySten;
  std::vector<int> edgeBoundarySten;
  std::vector<int> internalBoundarySten;

  std::vector<int> regBaseIndexVec;
  std::vector<int> regNumSpecieVec;
  std::vector<int> regElectronIndexVec;
  std::vector<int> regHoleIndexVec;

  std::vector<double> dxVec;  // mesh spacing.
  std::vector<double> xVec;   // mesh points.
  std::vector<double> CVec;   // doping
  std::vector<double> CdonorVec;    // doping
  std::vector<double> CacceptorVec; // doping
  std::vector<double> VVec;   // electrostatic potential
  std::vector<double> ExVec;  // electric field, x-direction.

  std::vector<double> JnxVec; // electron current density
  std::vector<double> JpxVec; // hole current density

  std::vector<double> RVec;   // recombination.
  std::vector<double> SVec;   // radiation source term.

  std::vector<double> nnVec;  // electron density
  std::vector<double> npVec;  // hole density

  //std::vector<double> unVec;  // spatially dependent mobility, electron
  //std::vector<double> upVec;  // spatially dependent mobility, hole
#if 0
  std::vector<double> unE_Vec; // mobility along edge, electron
  std::vector<double> upE_Vec; // mobility along edge, hole
#else
  std::vector<pdeFadType> unE_Vec; // mobility along edge, electron
  std::vector<pdeFadType> upE_Vec; // mobility along edge, hole
#endif

  std::vector<double> tnVec;  // spatially dependent lifetimes, electron
  std::vector<double> tpVec;  // spatially dependent lifetimes, hole

  // derivative arrays:

  std::vector<double> dRdpVec;
  std::vector<double> dRdnVec;

  std::vector<double> dJndn1Vec;
  std::vector<double> dJndn2Vec;
  std::vector<double> dJndV1Vec;
  std::vector<double> dJndV2Vec;
  std::vector<double> dJndp1Vec;
  std::vector<double> dJndp2Vec;

  std::vector<double> dJpdn1Vec;
  std::vector<double> dJpdn2Vec;
  std::vector<double> dJpdV1Vec;
  std::vector<double> dJpdV2Vec;
  std::vector<double> dJpdp1Vec;
  std::vector<double> dJpdp2Vec;

  // LID indices.
  std::vector<int> li_Vrowarray;
  std::vector< std::vector<int> > li_Vcolarray;

  std::vector<int> li_Nrowarray;
  std::vector< std::vector<int> > li_Ncolarray;

  std::vector<int> li_Prowarray;
  std::vector< std::vector<int> > li_Pcolarray;

  // columns needed for coupledMode==2
  std::vector< std::vector<int> > li_N_rxn_colarray;
  std::vector< std::vector<int> > li_P_rxn_colarray;

  std::vector<int> li_stateDispl;

  // matrix pointers
  std::vector< std::vector<double *> > fVmatPtr;
  std::vector< std::vector<double *> > fNmatPtr;
  std::vector< std::vector<double *> > fPmatPtr;
  std::vector< std::vector<double *> > qVmatPtr;
  std::vector< std::vector<double *> > qNmatPtr;
  std::vector< std::vector<double *> > qPmatPtr;

  // map between a mesh point index and a list of nearest neighbors
  // for that mesh point.
  std::multimap< int, int* > meshNeighborMultiMap;

  // state variable arrays associated with displacement current.
  std::vector<int> stateDispl;
  std::vector<int> stateDispl_owned;
  // this is an int array because I think bool arrays
  // might be problematic for some compilers.

  int numMeshPoints;  // this is probably redundant with NX.
  int maxColsPerRow;
  int numElectrodes;

  // 2d array of conductances, for 2-level "ckt phase" loads.
  std::vector< std::vector<double> > condVec;

  // data related to DMA matrix loads.
  std::vector<int> meshToLID;
  std::vector< std::vector<int> > jacStamp;
  std::vector<int> jacMap;
  std::vector< std::vector<int> > jacMap2;

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

private:


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
};

void registerDevice();

} // namespace DiodePDE
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DiodePDE::Instance N_DEV_DiodePDEInstance;
typedef Xyce::Device::DiodePDE::Model N_DEV_DiodePDEModel;

#endif
