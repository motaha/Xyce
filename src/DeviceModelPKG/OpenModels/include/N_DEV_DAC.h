//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DAC.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Lon Waters
//
// Creation Date  : 07/26/2002
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revsion$
//
// Revsion Date   : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DAC_h
#define Xyce_N_DEV_DAC_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>
#include <N_UTL_BreakPoint.h>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Device {
namespace DAC{

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "DAC";}
  static const char *deviceTypeName() {return "YDAC level 1 (Digital to Analog Interface)";};
  static const int numNodes() {return 2;}
  static const bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//----------------------------------------------------------------------------
// Class          : Instance
// Purpose        : This class refers to a single instance of the DAC device.
//                  It contains indices into the matrix equation. See comments
//                  for the ResistorInstance class for more details.
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 07/26/2002
//----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;friend class Master;

public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   DACiter,
     const FactoryBlock &      factory_block);


  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  // Additional Public Declarations
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  std::map<int,std::string> & getIntNameMap ();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  bool updateTVVEC ( std::vector< std::pair<double, double> > const & newPairs );
  bool getInstanceBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes);

  DeviceState * getInternalState();
  bool setInternalState( const DeviceState & state );

  // iterator reference to the model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

public:

  bool loadDAEQVector () {return true;};
  bool loadDAEFVector ();

  bool loadDAEdQdx () {return true;};
  bool loadDAEdFdx ();

  void varTypes( std::vector<char> & varTypeVec );

private:

  bool updateVoltage(double);

private:
  static std::vector< std::vector<int> > jacStamp;

  Model &       model_;         //< Owning model

  // user-specified parameters:
  std::string file;    // User specified file containing DAC data as time and
  //voltage pairs.

  std::vector< std::pair<double, double> > TVVEC; // vector of (time, voltage) pairs
  // read in from "file".

  // state variables:

  // other variables:
  int numTVpairs_;  // number of (time, voltage) pairs in TVVEC
  double v_pos;
  double v_neg;
  double i_bra;
  double vDrop;
  double voltage_;
  int loc_;

  // Indices into the state vector:



  //local indices (offsets)
  int li_Pos;
  int li_Neg;
  int li_Bra;

  //Locally indexed offsets for jacobian
  int ABraEquPosNodeOffset; // Offset, pos. node voltage contribution,
  // branch current equ.

  int ABraEquNegNodeOffset; // Offset, neg. node voltage contribution,
  // branch current equ.

  int APosEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the pos node

  int ANegEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the neg node
};

//----------------------------------------------------------------------------
// Function       : Model
// Purpose        :
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 07/26/2002
//----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Traits;friend class Master;

public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &        MB,
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

private:

  // Model Parameters
  double riseTime;
  double fallTime;
  double R;
  double L;
  double C;
  bool includeTransitionBP_;

  // Data Members for Associations

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

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;
  friend class Traits;

public:
  Master(
     const Configuration &       configuration,
     const FactoryBlock &      factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);
  virtual bool updateSecondaryState (double * staDeriv, double * stoVec);

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
  virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

  bool getDACDeviceNames ( std::vector<std::string> & dacNames);
};

void registerDevice();

} // namespace DAC
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DAC::Instance N_DEV_DACInstance;
typedef Xyce::Device::DAC::Model N_DEV_DACModel;
typedef Xyce::Device::DAC::Master N_DEV_DACMaster;

#endif
