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
// Filename       : $RCSfile: N_DEV_Resistor3.h,v $
//
// Purpose        : Resistor with zero resistance.
//
// Special Notes  : This is a special case of the resistor device when it is
//                  specified with zero resistance.  It will then behave like
//                  a voltage source with zero voltage difference so that a
//                  branch current is calculated.
//
// Creator        : Richard Schiek, Electrical and Microsystems Modeling.
//
// Creation Date  : 02/01/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.31.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Resistor3_h
#define Xyce_N_DEV_Resistor3_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Resistor.h>

namespace Xyce {
namespace Device {
namespace Resistor3 {

// ---------- Forward Declarations ----------
class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Resistor::Traits>
{
  static const char *name() {return "Resistor";}
  static const char *deviceTypeName() {return "R level 3";}
  static const int numNodes() {return 2;}
  static const char *primaryParameter() {return "R";}
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
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;
  friend class Master;

public:
  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   Viter,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  std::map<int,std::string> & getIntNameMap ();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();

  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  // load functions, residual:
  bool loadDAEQVector () { return true; }
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () { return true; }
  bool loadDAEdFdx ();

  void setupPointers ();

  double getMaxTimeStepSize ();

  void varTypes( std::vector<char> & varTypeVec );

  void getLIDs(int & lpos, int & lneg,int & lbra)
  {lpos = li_Pos; lneg = li_Neg; lbra = li_Bra;}

public:
  // iterator reference to the vsrc model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;
  static std::vector< std::vector<int> > jacStampPDE;
  static ParametricData<Instance>       parMap_;

  Model &       model_;         //< Owning model

  // state variables:
  double srcCurrent;
  double srcVoltage;
  double srcDrop;
  double srcBC;

  // scale factor
  double scale;
  int nlstep;

  // Parameters
  // user-specified paramters:
  double R;  // resistance  (ohms)
  // these are for the semiconductor resistor
  double length;      // resistor length.
  double width;      // resistor width.
  double temp;   // temperature of this instance
  // temperature dependence parameters
  // these can override values specified in the model
  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.
  double dtemp;        // externally specified device temperature.  This parameter
  // is NOT used and is only here for compatibility in parsing
  // netlist from simulators that do support it.
  // flags used to tell if the user has specified one of these values
  // on the command line.
  bool tempCoeff1Given;
  bool tempCoeff2Given;
  bool dtempGiven;

  // load variables
  double source, v_pos, v_neg, i_bra;

  // indices into state vector:
  int istate_I;  // index for i0;

  // Matrix equation index variables:

  //local indices (offsets)
  int li_Pos;
  int li_Neg;
  int li_Bra;

  // Jacobian matrix indices:
  //Locally indexed offsets for jacobian
  int ABraEquPosNodeOffset; // Offset, pos. node voltage contribution,
  // branch current equ.

  int ABraEquNegNodeOffset; // Offset, neg. node voltage contribution,
  // branch current equ.

  int APosEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the pos node

  int ANegEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the neg node

  //  The following jacobian offsets are only neccessary
  // for 2-level newton.
  int APosEquPosNodeOffset;  // Offset, positive node variable
  // contribution, positive node KCL.

  int ANegEquNegNodeOffset;  // Offset, negative node variable
  // contribution, negative node KCL.

  int ABraEquBraVarOffset;  // Offset, branch current variable
  // contribution, branch current equation.

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Jacobian matrix pointers:
  double * fBraEquPosNodePtr;
  double * fBraEquNegNodePtr;
  double * fPosEquBraVarPtr;
  double * fNegEquBraVarPtr;

  //  The following jacobian pointers are only neccessary for 2-level newton.
  double * fPosEquPosNodePtr;
  double * fNegEquNegNodePtr;
  double * fBraEquBraVarPtr;
#endif
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  friend class ParametricData<Model>;
  friend class Instance;
  friend class Traits;
  friend class Master;

  typedef std::vector<Instance *> InstanceVector;

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
  virtual bool processParams() 
  {
    return true;
  }

  virtual bool processInstanceParams() 
  {
    return true;
  }

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

  // This is the dc and transient analysis value of the source.
  double DC_TRAN;
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
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

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
  virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
};

void registerDevice();

} // namespace Resistor3
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Resistor3::Instance N_DEV_Resistor3Instance;
typedef Xyce::Device::Resistor3::Model N_DEV_Resistor3Model;
typedef Xyce::Device::Resistor3::Master N_DEV_Resistor3Master;

#endif
