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
// Filename       : $RCSfile: N_DEV_ExternDevice.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/09/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.59 $
//
// Revision Date  : $Date: 2014/02/24 23:49:16 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ExternDevice_h
#define Xyce_N_DEV_ExternDevice_h

#include <N_DEV_Configuration.h>
#include <N_DEV_fwd.h>
#include <N_PDS_fwd.h>

#include <N_TIA_TwoLevelError.h>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>
#include <N_DEV_CompositeParam.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : VoltageNode
// Purpose       : This class contains user specification of a voltage node.
//
// Special Notes : This is needed for vector-composite.
//
// The vsrcName and initVal variables come from composite params.  For example:
//
//  yext y1 vdd vss 1a 2a externcode=xyce netlist=sixNode_in.cir
//  + voltlim=1
//  + node={name=vconnectvdd,vconnectvss,vconnect1a,vconnect2a
//  +       initval = 6, 1, 6, 1}
//
//  The vector-composite is specified by "node".  The above specification should
//  result in 4 of these classes being allocated.
//
//  In the first  class, vsrcName = "vconnectvdd" and initval=6
//  In the second class, vsrcName = "vconnectvss" and initval=1
//  etc.
//
// Creator       : Eric Keiter
// Creation Date : 04/03/06
//-----------------------------------------------------------------------------
class VoltageNode : public CompositeParam
{
    friend class ParametricData<VoltageNode>;

  public:
    static ParametricData<VoltageNode> &getParametricData();

    VoltageNode();

    void processParams ();

  private:

  public:
    std::string vsrcName;   // name of the electrode, set by the user.
    double initVal;    // initial value of the node, set by the user.

    double limValHigh; // voltage limiter high value
    double limValLow;  // voltage limiter low  value
};



namespace ExternDevice {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
    static const char *name() {return "External Device";}
    static const char *deviceTypeName() {return "YEXT level 1 (External Device)";};
    static const int numNodes() {return 2;}
    static const int numOptionalNodes() {return 1000;}
    static const bool isLinearDevice() {return false;}

    static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
    static void loadModelParameters(ParametricData<Model> &model_parameters);
    static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
    friend class Traits;friend class DeviceMaster<Traits>;

public:

  Instance(
    const Configuration &       configuration,
    const InstanceBlock &     instance_block,
    Model &             model,
    const FactoryBlock &factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature ( const double & temp = -999.0 );

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  bool isInnerSolveConverged();
  bool getBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes);

  bool runExternalDevice();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setNonlocallyConnected() { locallyConnected_ = false; }
  void setOwningProc(int proc) { owningProc_ = proc; }
  void setComm(N_PDS_Comm* comm) { comm_ = comm; }

  void homotopyStepSuccess(const std::vector<std::string> & paramNames, const std::vector<double> & paramVals);

  void homotopyStepFailure ();

  void stepSuccess (int analysis);
  void stepFailure (int analysis);
  bool getInitialQnorm (N_TIA_TwoLevelError & tle);
  bool getInnerLoopErrorSum (N_TIA_TwoLevelError & tle);
  bool updateStateArrays();
  bool startTimeStep ();
  bool setInternalParam (std::string & name, double val);
  bool initialize ();

  CompositeParam *constructComposite(const std::string &composite_name, const std::string &);

private:
  bool setupVoltageLimiting_ ();
  bool calcVoltLimFactors_ ();

public:
  std::map<std::string,CompositeParam *> nodeMap;

  Model &getModel() {
    return model_;
  }

private:
  std::string externCode_;
  std::string netlistFileName_;
  std::vector< std::vector<int> > jacStamp;
  std::vector< std::vector<int> > jacLIDs;
  N_TIA_TwoLevelError tlError_;

  // vector of electrode data:
  std::vector<VoltageNode *> voltageNodeVec;

  ExternCodeInterface * extCodePtr_;

  // iterator reference to the resistor model which owns this instance.
  // Getters and setters
private:

  Model &       model_;         //< Owning model

  bool initializeFlag_;

  bool innerSolveStatus_;

  std::map<std::string,double> voltageInputMap_;
  std::vector<double> currentOutputVector_;
  std::vector<std::vector<double> > conductanceJacobian_;

  // voltage vectors used in limiting:
  std::vector<double> voltageOld_;
  std::vector<double> voltageOrig_;
  std::vector<double> voltageLastCall_;
  std::vector<double> voltageDiff_;
  std::vector<double> voltageFactor_;

  std::vector<int> voltageStateID_;
  int newtonIterOld;
  bool voltageLimiterFlag;
  bool initJctGiven_;
  bool nodeGiven_;

  bool locallyConnected_;
  int  owningProc_;
  N_PDS_Comm * comm_;
  const IO::CmdParse &      commandLine_;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/05
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
    friend class Traits;friend class DeviceMaster<Traits>;

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
  virtual bool processInstanceParams() {
    return true;
  }

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

};

//-----------------------------------------------------------------------------
// Function      : Instance::isConverged ()
// Purpose       : Return true if the most recent inner solve was
//                 successful, false if not.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
inline bool Instance::isInnerSolveConverged()
{
  return innerSolveStatus_;
}

//-----------------------------------------------------------------------------
// Function      : VoltageNode::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/03/06
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const VoltageNode & vn)
{
  os << "VoltageNode:\n"
     << "  vsrcName   = " << vn.vsrcName << "\n"
     << "  initVal    = " << vn.initVal << "\n"
     << "  limValHigh = " << vn.limValHigh << "\n"
     << "  limValLow  = " << vn.limValLow << "\n"
     << std::endl;

  return os;
}

void registerDevice();

} // namespace ExternDevice
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ExternDevice::Instance N_DEV_ExternDeviceInstance;
typedef Xyce::Device::ExternDevice::Model N_DEV_ExternDeviceModel;

#endif

