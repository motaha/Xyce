//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_NeuronPop1.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 06/10/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_NeuronPop1_h
#define Xyce_N_DEV_NeuronPop1_h

#include <N_DEV_Configuration.h>
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_UTL_fwd.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#ifdef HAVE_MATH_H
#include <math.h>
#endif
#include <fstream>

namespace Xyce {
namespace Device {
namespace NeuronPop1 {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "NeuronPopulation";}
  static const char *deviceTypeName() {return "YNEURONPOP level 1";}
  static const int numNodes() {return 2;}
  static const bool modelRequired() {return true;}
  static const bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 Neuron device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the NeuronInstance
//                 class for a more detailed explanation.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;
    
public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       IB,
     Model &                     Miter,
     const FactoryBlock &        factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  std::map<int,std::string> & getIntNameMap ();
  bool loadDeviceMask ();
  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);

  // this function is used to communicate to the device manager when
  // this device is changing (i.e. updating the population)
  bool getInstanceBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes);

  // updates done during the non-linear solve
  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool setIC ();

  void varTypes( std::vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  // enable the interface to produce plot files.
  bool plotfileFlag () {return true;}
  bool outputPlotFiles ();

  // iterator reference to the Neuron model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  void initializePopulation();
  void updatePopulation();

private:

  Model &       model_;         //< Owning model

  std::vector< std::vector<int> > jacStamp;

  // local indices into solution vector
  int liNodeIn;
  int liNodeOut;

  // local indicies into state vector
  std::vector<int> liNeuronPopState;

  // data needed for neurons
  bool populationInitialized;       // a flag for doing the initialization calculations
  double lastPopulationUpdateTime;  // time when we last did an update to avoid doing multiple updates per time
  double lastNeurogenesisUpdateTime;  // time when we last had a neurogenesis event - to avoid doing multiple updates per time
  int neuronPopSize;                // the current size of the neuron population
  //std::vector<float> neuronXpos;
  //std::vector<float> neuronYpos;

  std::vector<std::string> connectionTargetPopulation;
  bool connectionTargetPopulationGiven;

  // data for output of the population
  // output stream for output of internal state if requested by user
  Teuchos::RefCountPtr< std::ofstream > outputFileStreamPtr;
  bool outputPopulationVarsFlag;      // this flag indicates that the user wants to output the population state
  bool newStateToOutput;              // this flag is used to output only at times when the population has changed

  // flags for updating the population
  int numberOfUpdatesDone;

  // jacobian stamp offsets
  int jsOffsetNodeIn;
  int jsOffsetNodeOut;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Model : public DeviceModel
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

  // model parameters
  int neuronsMax;
  bool neuronsMaxGiven;
  int internalMaxConnections;
  bool internalMaxConnectionsGiven;
  int externalMaxConnections;
  bool externalMaxConnectionsGiven;
  double populationNeurogenesisRate;
  bool populationNeurogenesisRateGiven;

  // time at which population statistics are updated.
  double populationUpdatePeriod;
  bool populationUpdatePeriodGiven;

  int outputPopulationVars;  // flag indicating if user wants neuron populaiton output
  // (position, voltage, connectivity etc.)
};

void registerDevice();

} // namespace NeuronPop1
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::NeuronPop1::Instance N_DEV_NeuronPopInstance1;
typedef Xyce::Device::NeuronPop1::Model N_DEV_NeuronPopModel1;

#endif
