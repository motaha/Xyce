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
// Filename       : $RCSfile: N_DEV_Capacitor.h,v $
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
// Revision Number: $Revision: 1.123.2.2 $
//
// Revision Date  : $Date: 2014/03/06 21:33:44 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Capacitor_h
#define Xyce_N_DEV_Capacitor_h

#include <N_DEV_fwd.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Capacitor {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Capacitor";}
  static const char *deviceTypeName() {return "C level 1";}

  static const int numNodes() {return 2;}
  static const char *primaryParameter() {return "C";}
  static const char *instanceDefaultParameter() {return "C";}
  static const bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Capacitor::Instance
// Special Notes : A capacitor  will have two circuit nodes.
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
/// Capacitor instance
///
/// This class refers to a single instance of the capacitor device.  It
/// contains indicies into the matrix equation.  See the comments for the
/// Resistor::Instance class for more details.
///
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;
  friend class Master;

public:
  Instance(
     const Configuration &     configuration,
     const InstanceBlock &     instance_block,
     Model &                   model,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef, const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef );

  std::map<int,std::string> & getIntNameMap ();
  std::map<int,std::string> & getStoreNameMap();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars () { return true; };
  bool updatePrimaryState ();

  bool setIC ();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers();

  void varTypes( std::vector<char> & varTypeVec );

  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  Model &       model_;         //< Owning model

  // Stuff for handling solution-variable-dependent capacitance
  Util::Expression * expPtr;
  int                expNumVars;

  std::vector<double> expVarDerivs;

  // user-specified parameters:
  double C;    // User specified capacitance. (Farads)
  double IC;   // Optional initial value capacitor voltage (V).

  // These are for the semiconductor capacitor
  double length;    // capacitor length
  double width;     // capacitor width
  double temp;      // temperature of this instance

  // Genie 121412. temperature dependence parameters
  // these can override values specified in the model
  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.

  // flags used to tell if the user has specified one of these values
  // on the command line.
  bool tempCoeff1Given;
  bool tempCoeff2Given;

  // These are for the age-aware capacitor
  double age;                 ///< age in hours
  double ageCoef;             ///< degradation coeficient.
  double baseCap;             ///< the baseline capacitance before aging

  bool tempGiven;
  bool ICGiven;
  bool solVarDepC;

  // state variables:
  double q0;                  ///< charge in the capacitor
  // now held in the store vector at li_store_dev_i
  double vcap; // voltage drop across capacitor

  //local id's (offsets)
  int li_Pos;
  int li_Neg;
  int li_Bra;                 ///< for the "voltage source" when IC is specified

  int li_QState;

  std::vector<int> li_dQdXState;
  std::vector<int> li_dCdXState;
  int li_vcapState;
  int li_capState;

  int li_store_dev_i;

  // Offsets for Jacobian
  int APosEquPosNodeOffset;
  int ANegEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int ANegEquNegNodeOffset;

  // offsets for when C is solution-variable dependent
  std::vector<int> APosEquDepVarOffsets;
  std::vector<int> ANegEquDepVarOffsets;

  int ABraEquPosNodeOffset;
  int ABraEquNegNodeOffset;
  int ABraEquBraNodeOffset;
  int APosEquBraNodeOffset;
  int ANegEquBraNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointers for Jacobian
  double * qPosEquPosNodePtr;
  double * qNegEquPosNodePtr;
  double * qPosEquNegNodePtr;
  double * qNegEquNegNodePtr;

  double * fBraEquPosNodePtr;
  double * fBraEquNegNodePtr;
  double * fBraEquBraNodePtr;
  double * fPosEquBraNodePtr;
  double * fNegEquBraNodePtr;

  std::vector<double *> qPosEquDepVarsPtrs;
  std::vector<double *> qNegEquDepVarsPtrs;
#endif

  std::vector< std::vector<int> > jacStamp;
  std::vector< std::vector<int> > jacStamp_IC;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
/// Capacitor Model class
///
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Traits;
  friend class Master;

public:
  Model(
     const Configuration &     configuration,
     const ModelBlock &        model_block,
     const FactoryBlock &      factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  bool processParams ();
  bool processInstanceParams ();
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

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

  // for the semiconductor capacitor
  double cj;     // junction bottom capacitance
  double cjsw;   // junction sidewall capacitance
  double defWidth; // default width
  double narrow;   // narrowing due to side etching
  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.
  double baseCap;
  double tnom;

  bool tnomGiven;
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Capacitor Master class
///
/// The "master" class is the one that contains the updateState, loadDAEVectors
/// and loadDAEMatrices methods that are actually called when it is time to
/// compute and load device contributions.
///
/// The default implementations of these methods in the DeviceMaster
/// template class simply loops over all instances and calls their
/// updatePrimaryState, loadDAEFVector/loadDAEQVector, and
/// loadDAEdFdx/loadDAEdQdx methods, respectively.
///
/// For efficiency, the Capacitor class reimplements these methods to do the
/// work directly, instead of calling instance-level functions.
///
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:
  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       ss1,
     const DeviceOptions &     do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
  virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
};

void registerDevice();

} // namespace Capacitor
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Capacitor::Instance N_DEV_CapacitorInstance;
typedef Xyce::Device::Capacitor::Model N_DEV_CapacitorModel;
typedef Xyce::Device::Capacitor::Master N_DEV_CapacitorMaster;

#endif // Xyce_N_DEV_Capacitor_h

