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
// Filename       : $RCSfile: N_DEV_Resistor.h,v $
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
// Revision Number: $Revision: 1.111.2.5 $
//
// Revision Date  : $Date: 2014/03/06 02:55:06 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Resistor_h
#define Xyce_N_DEV_Resistor_h

#include <N_DEV_fwd.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceMaster.h>

namespace Xyce {
namespace Device {
namespace Resistor {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Resistor";}
  static const char *deviceTypeName() {return "R level 1";}
  static const int numNodes() {return 2;}
  static const char *primaryParameter() {return "R";}
  static const char *instanceDefaultParameter() {return "R";}
  static const bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &p);
  static void loadInstanceParameters(ParametricData<Instance> &p);
};

//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Resistor::Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
///
/// Resistor device instance class.
///
/// An instance is created for each occurance of the device in the netlist.
///
/// It contains "unique" resistor information - ie stuff that will be
/// true of only one resistor in the circuit, such as the nodes to
/// which it is connected.  A resistor is connected to only two
/// circuit nodes.
///
/// This class does not directly contain information about its node
/// indices. It contains indices into the 5 parts (dFdx, dQdx, dx, F,
/// and Q) of the matrix problem A*dx = b, and also the solution
/// vector x.  A is the Jacobian matrix that will be formed from dFdx
/// and d(dQ/dt)dx, dx is the update to x, and b is right hand side
/// function vector that will be formed from F and dQ/dt.  These
/// indices are global, and determined by topology during the
/// initialization stage of execution.
///
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;              ///< Allow ParametricData to changes member values
  friend class Model;
  friend class Traits;
  friend class Master;

public:
  Instance(
     const Configuration &     configuration,
     const InstanceBlock &     instance_block,
     Model &                   model,
     const FactoryBlock &      factory_block);

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Instance::~Instance
  // Purpose       : destructor
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 3/16/00
  //---------------------------------------------------------------------------
  ///
  /// Destroys this instance
  ///
  /// @author Eric Keiter, SNL
  /// @date   3/16/00
  ~Instance() {}

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Instance::getModel
  // Purpose       : destructor
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 3/16/00
  //---------------------------------------------------------------------------
  ///
  /// Gets the resistor model that owns this instance.
  ///
  /// @return reference to the owning Resistor::Model
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  /// @date   Mon Aug 12 08:36:37 2013
  Model &getModel()   { return model_;  }

  virtual void registerLIDs(const std::vector<int> & intLIDVecRef, const std::vector<int> & extLIDVecRef) /* override */;
  virtual void registerStateLIDs(const std::vector<int> & staLIDVecRef) /* override */;
  virtual void registerStoreLIDs(const std::vector<int> & stoLIDVecRef) /* override */;
  virtual void registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec) /* override */;
  virtual std::map<int, std::string> &getStoreNameMap() /* override */;

  virtual bool processParams() /* override */;
  virtual bool updateTemperature(const double & temp_tmp) /* override */;
  virtual bool updateIntermediateVars() /* override */;
  virtual bool updatePrimaryState() /* override */;

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Instance::jacobianStamp
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
  // Creation Date : 08/20/01
  //---------------------------------------------------------------------------
  ///
  /// Return Jacobian stamp that informs topology of the layout of the
  /// resistor jacobian.
  ///
  /// The Jacobian stamp describes the shape of the Jacobian to the
  /// Topology subsystem.  The Topology subsystem, in turn, returns
  /// the offsets into the matrix and solution vectors where this
  /// instance data is located.
  ///
  /// @return const reference to a std::vector of std::vector of
  /// integers describing Jacobian stamp shape
  ///
  /// @author Robert Hoekstra
  /// @date 8/20/2001
  virtual const std::vector< std::vector<int> > &jacobianStamp() const  /* override */ {
    return jacStamp;
  }

  virtual bool loadDAEFVector() /* override */;
  virtual bool loadDAEdFdx() /* override */;

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Instance::loadDAEQVector
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter
  // Creation Date : 1/24/03
  //---------------------------------------------------------------------------
  ///
  /// Load Q vector
  ///
  /// @return true on success
  /// 
  /// Since the Resistor does no charge storage, this is a no-op.
  ///
  /// @author Eric Keiter
  /// @date   1/24/03
  virtual bool loadDAEQVector() /* override */ 
  {
    return true;
  }

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Instance::loadDAEdQdx
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter
  // Creation Date : 1/24/03
  //---------------------------------------------------------------------------
  ///
  /// Load derivative of Q vector with respect to solution vector
  ///
  /// @return true on success
  ///
  /// Since the Resistor does no charge storage, this is a no-op.
  ///
  /// @author Eric Keiter
  /// @date   1/24/03
    virtual bool loadDAEdQdx() /* override */ 
  {
    return true;
  }

  virtual void setupPointers() /* override */;

private:
  static std::vector< std::vector<int> >  jacStamp; ///< All Resistor instances have a common Jacobian Stamp
  static void initializeJacobianStamp();

  Model &     model_;                 ///< Owning model

  // User-specified parameters:
  double      R;                      ///< Resistance (ohms)

  // These are for the semiconductor resistor
  double      length;                 ///< Resistor length.
  double      width;                  ///< Resistor width.
  double      temp;                   ///< Temperature of this instance

  // Temperature dependence parameters, these can override values specified in the model
  double      tempCoeff1;             ///< First order temperature coeff.
  double      tempCoeff2;             ///< Second order temperature coeff.
  double      dtemp;                  ///< Externally specified device temperature.
  ///<   NOT used, only here for compatibility in parsing
  ///<   netlist from simulators that support it

  // Flags used to tell if the user has specified one of these values on the command line.
  bool        tempCoeff1Given;        ///< First order temperation value was given in netlist
  bool        tempCoeff2Given;        ///< Second order temperature coeff was given in netlist
  bool        dtempGiven;             ///< Externally specified device temperature was given in netlist

  // Derived parameters:
  double      G;                      ///< Conductance(1.0/ohms)
  double      i0;                     ///< Current(ohms)

  int         li_Pos;                 ///< Index for Positive Node
  int         li_Neg;                 ///< Index for Negative Node
  int         li_store_dev_i;         ///< Index for Lead Current

  // Offset variables corresponding to the above declared indices.
  int         APosEquPosNodeOffset;   ///< Column index into force matrix of Pos/Pos conductance
  int         APosEquNegNodeOffset;   ///< Column index into force matrix of Pos/Neg conductance
  int         ANegEquPosNodeOffset;   ///< Column index into force matrix of Neg/Pos conductance
  int         ANegEquNegNodeOffset;   ///< Column index into force matrix of Neg/Neg conductance

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointers for Jacobian
  double *    f_PosEquPosNodePtr;
  double *    f_PosEquNegNodePtr;
  double *    f_NegEquPosNodePtr;
  double *    f_NegEquNegNodePtr;
#endif
};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Resistor::Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
///
/// Resistor model class
///
class Model : public DeviceModel
{
  friend class ParametricData<Model>;               ///< Allow ParametricData to changes member values
  friend class Instance;                            ///< Don't force a lot of pointless getters
  friend class Traits;
  friend class Master;                              ///< Don't force a lot of pointless getters

public:
  typedef std::vector<Instance *> InstanceVector;

  Model(
     const Configuration &       configuration,
     const ModelBlock &        model_block,
     const FactoryBlock &      factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Model::addInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David Baur
  // Creation Date : 8/12/2013
  //---------------------------------------------------------------------------
  ///
  /// Add an instance to the list of instances associated with this model
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  /// @date   8/12/2013
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Model::getInstanceVector
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David Baur
  // Creation Date : 8/12/2013
  //---------------------------------------------------------------------------
  ///
  /// Get a non-const reference to the vector for all resistor
  /// instances owned by this model.
  ///
  /// @return reference to InstanceVector containing all resistors owned by this model
  /// 
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  /// @date   8/12/2013
  InstanceVector &getInstanceVector() 
  {
    return instanceContainer;
  }

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Model::getInstanceVector
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David Baur
  // Creation Date : 8/12/2013
  //---------------------------------------------------------------------------
  ///
  /// Get a const reference to the instance vector for all resistors
  /// owned by this model.
  /// 
  /// @return const reference to InstanceVector containing all resistors owned by this model
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  /// @date   Mon Aug 12 09:10:00 2013
  const InstanceVector &getInstanceVector() const 
  {
    return instanceContainer;
  }

  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

  virtual bool processParams() /* override */;
  virtual bool processInstanceParams() /* override */;

private:
  InstanceVector      instanceContainer;            ///< List of owned resistor instances

  // Semiconductor resistor parameters
  double      tempCoeff1;     ///< First order temperature coefficient
  double      tempCoeff2;     ///< Second order temperature coefficient
  double      sheetRes;       ///< Sheet resistance
  double      defWidth;       ///< Default width
  double      narrow;         ///< Narrowing due to side etching
  double      tnom;           ///< Parameter measurement temperature
};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Resistor::Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Resistor master
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
/// For efficiency, the Resistor class reimplements these methods to do the
/// work directly, instead of calling instance-level functions.
///
class Master : public DeviceMaster<Traits>
{
  friend class Instance;                            ///< Don't force a lot of pointless getters
  friend class Model;                               ///< Don't force a lot of pointless getters

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Master::Master
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter
  // Creation Date : 11/26/08
  //---------------------------------------------------------------------------
  ///
  /// Construct a Resistor Device.
  ///
  /// @param configuration
  /// @param factory_block
  /// @param solver_state
  /// @param device_options
  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
    : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options)
  {}

  virtual bool updateState(double * solVec, double * staVec, double * stoVec) /* override */;
  virtual bool loadDAEVectors(double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ) /* override */;
  virtual bool loadDAEMatrices(N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx) /* override */;
};

void registerDevice();

} // namespace Resistor
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Resistor::Instance N_DEV_ResistorInstance;
typedef Xyce::Device::Resistor::Model N_DEV_ResistorModel;
typedef Xyce::Device::Resistor::Master N_DEV_ResistorMaster;

#endif // Xyce_N_DEV_Resistor_h
