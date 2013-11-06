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
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_NewDevice.h,v $
//
// Purpose        : NewDevice classes.
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
// Revision Number: $Revision: 1.9.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_NewDevice_h
#define Xyce_N_DEV_NewDevice_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace NewDevice {

// ---------- Forward Declarations ----------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 NewDevice device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the NewDeviceInstance
//                 class for a more detailed explanation.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Master;

  public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

    Instance(InstanceBlock & IB,
                           Model & Miter,
                           MatrixLoadData & mlData1,
                           SolverState &ss1,
                           ExternData  &ed1,
                           DeviceOptions & do1);


    ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
    void registerLIDs( const vector<int> & intLIDVecRef,
                       const vector<int> & extLIDVecRef );
    void registerStateLIDs( const vector<int> & staLIDVecRef );

    map<int,string> & getIntNameMap ();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");
    bool updateTemperature(const double & temp_tmp);

    bool updateIntermediateVars () {return true;};
    bool updatePrimaryState ();
    bool updateSecondaryState ();
    bool setIC ();

    void varTypes( vector<char> & varTypeVec );

    // load functions, residual:
    bool loadDAEQVector ();
    bool loadDAEFVector ();

    void auxDAECalculations ();

    // load functions, Jacobian:
    bool loadDAEdQdx ();
    bool loadDAEdFdx ();

  public:
    // iterator reference to the NewDevice model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:
    static vector< vector<int> > jacStamp;

  Model &       model_;         //< Owning model

    // parameter variables

    // state variables

    // local state indices (offsets)

    // local solution indices (offsets)
    int li_Pos;      // local index to positive node on this device
    int li_Neg;      // local index to negative node on this device

    // Matrix equation index variables:

    // Offset variables corresponding to the above declared indices.
    int APosEquPosNodeOffset;
    int APosEquNegNodeOffset;
    int ANegEquPosNodeOffset;
    int ANegEquNegNodeOffset;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Master;

  public:
  static ParametricData<Model> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

    Model(const ModelBlock & MB,
                              SolverState & ss1,
                           DeviceOptions & do1);
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
    vector<Instance*> instanceContainer;

 public:
  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

 private:
};

} // namespace NewDevice
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::NewDevice::Instance N_DEV_NewDeviceInstance;
typedef Xyce::Device::NewDevice::Model N_DEV_NewDeviceModel;

#endif
