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
//
//-----------------------------------------------------------------------------
// Filename       : N_DEV_OpAmp.h
//
// Purpose        : Ideal Operational Amplifier Model
//
// Special Notes  : This model assumes infinite gain
//                   and unbounded output voltages
//
// Creator        : Brian Fett, SNL
//
// Creation Date  : 08/02/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.31.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_OpAmp_h
#define Xyce_N_DEV_OpAmp_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace OpAmp {

// ---------- Forward Declarations ----------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Brian Fett
// Creation Date : 08/28/05
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
                        Model & iter,
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

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");

    bool updateIntermediateVars ();
    bool updatePrimaryState ();

    // load functions, residual:
    bool loadDAEQVector () {return true;}
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx () {return true;}
    bool loadDAEdFdx ();

    // iterator reference to the OpAmp model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:
  static vector< vector<int> > jacStamp;

  Model &       model_;         //< Owning model

    // state variables:
    double outCurrent;
    double deltaVoltage;


    // Parameters
    double FAKEPARAM;

    // load variables
    double v_pos, v_neg, v_out, i_bra;

    // indices into state vector:
    int istate_I;  // index for i0;

    // Matrix equation index variables:
    //local indices (offsets)
    int li_Pos;
    int li_Neg;
    int li_Out;
    int li_Bra;

    // Jacobian matrix indices:
    //Locally indexed offsets for jacobian
    int ABraEquPosNodeOffset; // Offset, pos. input voltage contribution
                              // branch current equation.

    int ABraEquNegNodeOffset; // Offset, neg. input voltage contribution
                              // branch current equation.

    int AOutEquBraVarOffset;  // Offset, branch current variable
                              // contribution, KCL of output node
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
// Modified by   : Brian Fett, SNL
// Modify Date   : 08/02/05
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

    Model (const ModelBlock & MB,
                            SolverState & ss1,
			    DeviceOptions & do1);
    ~Model ();

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
  vector <Instance*> instanceContainer;

   private:

    // Data Members for Associations
    double FAKEPARAM;
};

} // namespace OpAmp
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::OpAmp::Instance N_DEV_OpAmpInstance;
typedef Xyce::Device::OpAmp::Model N_DEV_OpAmpModel;

#endif
