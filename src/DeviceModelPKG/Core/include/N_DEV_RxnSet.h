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
// Filename       : $RCSfile: N_DEV_RxnSet.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Thomas V. Russo, SNL, Component Information and Models
//
// Creation Date  : 08/19/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.22.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_RxnSet_h
#define Xyce_N_DEV_RxnSet_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DevicePDEInstance.h>
#include <N_DEV_DevicePDEModel.h>
#include <N_DEV_CompositeParam.h>

#include <N_DEV_Param.h>
#include <N_UTL_BreakPoint.h>
#include <N_DEV_TransportHelper.h>

namespace Xyce {
namespace Device {

namespace RxnSet {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Rxn Effects Device";}
  static const char *deviceTypeName() {return "YRXN level 1 (Rxn Device)";};
  static const int numNodes() {return 2;}
  static const bool modelRequired() {return true;}
  static const bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This class refers to a single instance of the BJT
//                 device.  It contains indices into the matrix equation.
//                 See the comments for the ResistorInstance class for
//                 more details.
//
//                 The bjt will have 4 external nodes: collector, base,
//                 emitter, and substrate, and 3 internal nodes:
//                 collectorPrime, basePrime, and emitterPrime.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
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
     const InstanceBlock &       IB,
     Model &                     it_MB,
     const FactoryBlock &        factory_block);


  Instance(const Instance &right);

  ~Instance();

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & stateLIDVecRef );
  std::map<int,std::string> & getIntNameMap ();

  const std::vector<std::string> & getDepSolnVars();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature (const double & temp = -999.0 );

  bool getInstanceBreakPoints( std::vector<N_UTL_BreakPoint> &breakPointTimes);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  bool loadDeviceMask ();

  bool plotfileFlag () {return true;}

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  // Debug related load functions for the Jacobian:
  bool loadQMatrix (N_LAS_Matrix & dQdxMat);
  bool loadFMatrix (N_LAS_Matrix & dFdxMat);

  // Debugging Excess Phase function, etc.
  bool outputPlotFiles ();
  bool outputTecplot ();
  bool output2DTecplot ();
  bool outputCarrierDensities ();

  void setupJacStamp ();

  void setupMeshUniform ();

  void allocateRegions ();

  void scaleMesh ();
  void setupFluxVec ();

  void setupScalingVars ();
  void initializeChemistry ();

  void setupPointers();

public:
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

private:

  bool haveAnyReactions; // global haveAnyReactions flag for all regions.
  bool reactionFileCopyFlag;
  bool useScaledVariablesFlag;
  bool useDopingArrayData;

  double outputInterval;
  int  outputIndex;
  double lastOutputTime;

  int outputRegion;
  int tecplotLevel;
  int callsOTEC;
  int callsOTECcarrier;

  //external instance params
  double TEMP;  // instance temperature (TEMP)

  int  newtonIterOld;

  //local indexing of solution and state variables
  int li_Pos;
  int li_Neg;

  // reaction region(s):
  std::vector<Region*> regVec;

  std::vector<int> regLastIndexVec;
  std::vector<int> regFirstReactantIndexVec;
  std::vector<int> regNumSpecieVec;

  // these are relative indices for use in the jacStamp setup:
  std::vector< std::vector<int> > APosEqu_SpeciesOffset;
  std::vector< std::vector<int> > ANegEqu_SpeciesOffset;

  std::vector< std::vector<double *> > APosEqu_SpeciesPtr;
  std::vector< std::vector<double *> > ANegEqu_SpeciesPtr;

  std::vector< std::vector<double *> > APosEqu_ConstPtr;
  std::vector< std::vector<double *> > ANegEqu_ConstPtr;

  // mesh variables:
  std::vector<double> xVec;
  std::vector<double> dxVec;

  std::vector<int> xloStencilVec;
  std::vector<int> xhiStencilVec;

  std::vector<TransportHelper> thVec;

  double outputXscalar;

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;

  int ANegEquNegNodeOffset;
  int ANegEquPosNodeOffset;

  std::vector< std::vector<int> > jacStamp;
  std::vector<int> jacMap;
  std::vector< std::vector<int> > jacMap2;

  bool excludeNoSourceRegionsFlag;
  bool excludeNoSourceRegionsFlagGiven;
  bool transportFlagGiven;
  bool transportFlag;
  bool diffusionFlagGiven;
  bool diffusionFlag;
  bool dirichletBCFlag;
  bool columnReorderingFlag;

  int xloIndex;
  int xhiIndex;

  int callsOutputPlot;
};


//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class Model : public DevicePDEModel
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
  CompositeParam *constructComposite (const std::string & cName, const std::string & pName);

  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;
    
  virtual std::ostream &printOutInstances(std::ostream &os) const;

  bool processParams ();
  bool processInstanceParams ();


public:
  void addInstance(Instance *instance) {
    instanceContainer.push_back(instance);
  }

  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

private:
  std::vector<Instance*> instanceContainer;

private:

  //external model params
  double TNOM;               //nominal temperature

  int  userNumRegions;

  // File name for reaction specification:
  std::string rxnFileName;

  //*************************************
  // Rxn reaction model stuff:
  double xlo;
  double xhi;

  double xlo_source; // source region, low bound
  double xhi_source; // source region, high bound
  bool xlo_sourceGiven;
  bool xhi_sourceGiven;

  std::vector<RegionData*> regionDataVec;
  std::map<std::string,CompositeParam *> regionDataMap;
  std::map<std::string,SpecieSource *> defectSourceMap;

  double masterSource;
};

void registerDevice();

} // namespace RxnSet
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::RxnSet::Instance N_DEV_RxnSetInstance;
typedef Xyce::Device::RxnSet::Model N_DEV_RxnSetModel;

#endif
