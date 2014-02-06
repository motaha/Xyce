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
// Filename       : $RCSfile: N_DEV_ACC.h,v $
//
// Purpose        : ACC classes: provide a device that integrates the
//                  equations of motion of an accelerated body to give its
//                  instantaneous position and velocity
//
// Special Notes  : Intended for use in Coil Gun LDRD netlists.
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 10/24/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ACC_h
#define Xyce_N_DEV_ACC_h
// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace ACC {

// ---------- Forward Declarations -------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
//
//	This  is  the  instance class  for accelerated object devices.  It
//	contains "unique" ACC device  information - ie stuff that
//	will be true of only one ACC in the circuit, such
//	as the nodes to which it is connected.  An ACC is
//	connected to only one circuit node, from which it gets the
//      acceleration.
//
//	This class  does not directly contain information about
//	its node indices. It contains indices into the 3 parts
//	(A, dx, and  b) of the matrix  problem A*dx = b, and
//	also x.  A is the Jacobian  matrix, dx is the update to
//	the solution vector x, and b is the right hand side
//	function vector.  These indices are global, and
//	determined by topology during  the initialization stage
//	of execution.
//
// Special Notes :
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date  : 10/24/07
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;

public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Instance(InstanceBlock & IB,
           Model & Riter,
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

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

public:
  // iterator reference to the ACC model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:
  static vector< vector<int> >  jacStamp;

  Model &       model_;         //< Owning model

  // user-specified paramters:
  double v0; // initial velocity
  double x0; // initial position

  // Local indices:
  // Acceleration (an external node)
  int li_Acc;
  // velocity (internal node)
  int li_Velocity;
  // position (internal node)
  int li_Position;

  // state variable IDs:
  int li_state_vel;
  int li_state_pos;

  // Jacobian offsets:
  int AVelEquAccNodeOffset;
  int AVelEquVelNodeOffset;
  int APosEquVelNodeOffset;
  int APosEquPosNodeOffset;

  // Instance variables --- mostly copies of stuff from solution vector
  // to reduce calls to operator[]
  double position;
  double velocity;
  double acceleration;

  //things we take out of the state vector derivatives
  double xdot;
  double vdot;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
//
//
// Special Notes :
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date  : 10/24/07
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;

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

    virtual bool processParams(std::string param = "")
    {}

    virtual bool processInstanceParams(std::string param = "")
    {}


public:
  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

private:
  vector<Instance*> instanceContainer;

private:
};

} // namespace ACC
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ACC::Instance N_DEV_ACCInstance;
typedef Xyce::Device::ACC::Model N_DEV_ACCModel;

#endif // Xyce_N_DEV_ACC_h
