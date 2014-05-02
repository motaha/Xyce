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
// Filename       : $RCSfile: N_DEV_ISRC.h,v $
//
// Purpose        : Independent current source
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
// Revision Number: $Revision: 1.98.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ISRC_h
#define Xyce_N_DEV_ISRC_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>

#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {
namespace ISRC {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Independent Current Source";}
  static const char *deviceTypeName() {return "I level 1";}
  static const int numNodes() {return 2;}
  static const char *primaryParameter() {return "DCV0";}
  static const char *instanceDefaultParameter() {return "DCV0";}
  static const bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
class Instance : public SourceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;friend class Master;

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
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef);
  void registerStateLIDs( const std::vector<int> & staLIDVecRef);

  void registerStoreLIDs(const std::vector<int> & stoLIDVecRef );

  std::map<int,std::string> & getStoreNameMap ();

  const std::vector< std::vector<int> > & jacobianStamp() const;

  bool processParams ();

  bool updateIntermediateVars () { return true;};
  bool updatePrimaryState ();

  bool loadTrivialMatrixStamp ();
  bool loadTrivialDAE_FMatrixStamp ();

  // load functions, residual:
  bool loadDAEQVector () { return true; }
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () { return true; }
  bool loadDAEdFdx () { return true; }

  bool loadBVectorsforAC (double * bVecReal, double * bVecImag);

protected:
private:

public:
  // iterator reference to the ISRC model that owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;

  Model &       model_;         //< Owning model

  // index variables:
  int li_Pos;
  int li_Neg;

  // Store variables
  int li_store_dev_i;

  // Parameters
  double DCV0;
  double par0;
  double par1;
  double par2;
  double par3;
  double par4;
  double par5;
  double par6;
  double par7;
  double REPEATTIME;
  double T;
  double V;
  int NUM;
  bool REPEAT;
  int TRANSIENTSOURCETYPE;
  bool TRANSIENTSOURCETYPEgiven;
  int ACSOURCETYPE;
  bool ACSOURCETYPEgiven;
  int DCSOURCETYPE;
  bool DCSOURCETYPEgiven;
  bool gotParams;


  double ACMAG;
  double ACPHASE;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class Instance;
  friend class Traits;friend class Master;

public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &        MB,
     const FactoryBlock &      factory_block);
  ~Model ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;
  virtual bool processParams();
  virtual bool processInstanceParams();

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

//-----------------------------------------------------------------------------
// Function      : Instance::loadTrivialMatrixStamp ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/24/03
//-----------------------------------------------------------------------------
inline bool Instance::loadTrivialMatrixStamp ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadTrivialDAE_FMatrixStamp ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/24/03
//-----------------------------------------------------------------------------
inline bool Instance::loadTrivialDAE_FMatrixStamp ()
{
  return true;
}
//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:
  Master(
     const Configuration &       configuration,
     const FactoryBlock &      factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);

  // new DAE stuff:
  // new DAE load functions, residual:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);

  // new DAE load functions, Jacobian:
  virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
};

void registerDevice();

} // namespace ISRC
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ISRC::Instance N_DEV_ISRCInstance;
typedef Xyce::Device::ISRC::Model N_DEV_ISRCModel;
typedef Xyce::Device::ISRC::Master N_DEV_ISRCMaster;

#endif

