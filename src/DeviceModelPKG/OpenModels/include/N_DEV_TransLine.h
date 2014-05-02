//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, 2013, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
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
// Filename       : $RCSfile: N_DEV_TransLine.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 12/11/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.25.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_TransLine_h
#define Xyce_N_DEV_TransLine_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>


#define TRANS_MOD_RLC 1
#define TRANS_MOD_LC  2

namespace Xyce {
namespace Device {
namespace TransLine {

// ---------- Forward Declarations ----------
class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Lossy Transmission Line";}
  static const char *deviceTypeName() {return "YTRANSLINE level 1";}
  static const int numNodes() {return 2;}
  static const bool modelRequired() {return true;}
  static const bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

// note that V1 and V3 are shared by neighbor nodes.
struct lumpData
{
public:
  int indexV1;
  int indexV2;
  int indexI;
  int indexV3;

  int li_V1;
  int li_V2;
  int li_I;
  int li_V3;

  int offset_v1_v2m1; // from previous lump
  int offset_v1_iim1; // from previous lump
  int offset_v1_v1;
  int offset_v1_ii;

  int offset_v2_v2;
  int offset_v2_ii;
  int offset_v2_v3;

  int offset_ii_v1;
  int offset_ii_v2;
  int offset_ii_ii;

  int offset_v3_v2; // only used if second external node
  int offset_v3_v3;

  double f0_ind;
  double i0_ind;
  double coef_ind;
  double i0_res;
  double q0_cap;
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;

public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   Citer,
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
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  bool loadDeviceMask ();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers();

  void varTypes( std::vector<char> & varTypeVec );

public:
  // iterator reference to the resistor model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // Data Members for Class Attributes

  // User-specified parameters:
  int numLumps;
  int numTransLineVars;

  double length;

  bool numLumpsGiven;
  bool lengthGiven;

  double L,C,G,R;

  std::vector<lumpData>lumpVec;

  std::vector< std::vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------
class Model  : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Traits;

public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &      MB,
     const FactoryBlock &    factory_block);
  ~Model   ();

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
  // input parameters:
  int elevNumber;
  double resist;
  double induct;
  double conduct;
  double capac;

  bool elevNumberGiven;
  bool resistGiven;
  bool inductGiven;
  bool conductGiven;
  bool capacGiven;

  int specialCase;     // what kind of model (RC, RLC, RL, ...)

  std::vector<Instance*> instanceContainer;

private:

};

void registerDevice();

} // namespace TransLine
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::TransLine::Instance N_DEV_TransLineInstance;
typedef Xyce::Device::TransLine::Model N_DEV_TransLineModel;

#endif // Xyce__h

