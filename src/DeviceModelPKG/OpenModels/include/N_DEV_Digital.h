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
// Filename       : $RCSfile: N_DEV_Digital.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 01/05/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.56.2.3 $
//
// Revision Date  : $Date: 2014/03/12 16:50:27 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Digital_h
#define Xyce_N_DEV_Digital_h

// ----------   Xyce Includes   ----------
//Genie 110812
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Digital {

class Model;
class Instance; //Genie 110812

struct Traits : public DeviceTraits<Model, Instance>
{
    static const char *name() {return "Behavioral Digital";}
    static const char *deviceTypeName() {return "Digital level 1";}
    static const int numNodes() {return 2;}
    static const int numOptionalNodes() {return 20;}
    static const bool modelRequired() {return true;}
    static const bool isLinearDevice() {return true;}

    static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
    static void loadModelParameters(ParametricData<Model> &model_parameters);
    static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This class refers to a single instance of a digital
//                 device.  It contains indicies into the matrix equation.
//                 See the comments for the ResistorInstance class for
//                 more details.
//
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
    friend class ParametricData<Instance>;
    friend class Model;
    friend class Traits;
    friend class Master;

    // NOT is deprecated now, and replaced by INV.
    enum gType {INV, NOT, AND, NAND, OR, NOR, ADD, XOR, NXOR, DFF, DLTCH, BUF};

public:

  Instance(
    const Configuration &       configuration,
    const InstanceBlock &       IB,
    Model &                     Diter,
    const FactoryBlock &        factory_block);


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

  bool updateIntermediateVars () { return true; };
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool getInstanceBreakPoints (std::vector<N_UTL_BreakPoint> &);

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  // Used to support U vs. Y syntax
  std::string getDeviceLetter ();

public:
  // iterator reference to the model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // Data Members for Class Attributes
  // For more explanation of these attributes
  //   see the resistor classes.

  // state variables:

  std::vector<double> qlo;      // charge in the capacitor
  std::vector<double> ilo;      // current throught the capacitor
  std::vector<double> vcaplo;   // voltage drop across capacitor
  std::vector<double> qhi;
  std::vector<double> ihi;
  std::vector<double> vcaphi;
  std::vector<double> qref;
  std::vector<double> iref;
  std::vector<double> vcapref;

  std::vector<double> rilo;
  std::vector<double> rihi;
  std::vector<double> riref;
  std::vector<double> currentOut;
  std::vector<double> currentIn;

  std::vector<double> glo;
  std::vector<double> ghi;

  std::vector<double> qInp;      // charge in the capacitor
  std::vector<double> iInp;      // current throught the capacitor
  std::vector<double> vcapInp;   // voltage drop across capacitor

  std::vector<double> currentInp;

  // input params:

  bool ic1;
  bool ic2;
  bool ic3;

  int numInput;       // Number of input leads
  int numOutput;      // Number of output leads
  enum gType gate;

  //local id's (offsets)
  int li_Lo;
  int li_Hi;
  int li_Ref;
  std::vector<int> li_Inp;
  std::vector<int> li_Out;

  // Input state vars
  std::vector<int> li_currentStateInp;
  std::vector<int> li_transitionTimeInp;
  std::vector<int> li_QinpState;
  std::vector<int> li_IinpState;

  // Output state vars
  std::vector<int> li_currentStateOut;
  std::vector<int> li_transitionTimeOut;
  std::vector<int> li_QloState;
  std::vector<int> li_IloState;
  std::vector<int> li_QhiState;
  std::vector<int> li_IhiState;

  std::vector<bool> inpL;
  std::vector<double> iTime;
  std::vector<bool> outL;
  std::vector<double> oTime;

  double breakTime;

  // Offsets for Jacobian
  int row_Lo;
  int row_Hi;
  int row_Ref;
  std::vector< std::vector<int> > li_jac_Ref;
  std::vector< std::vector<int> > li_jac_Lo;
  std::vector< std::vector<int> > li_jac_Hi;

  std::vector< std::vector<int> > jacStamp;

  //Genie 110812. change state var
  //std::vector<bool> changeState;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
class Model  : public DeviceModel
{
    typedef std::vector<Instance *> InstanceVector;
    
    friend class ParametricData<Model>;
    friend class Instance;
    friend class Traits;friend class Master;

public:
  Model(
    const Configuration &       configuration,
    const ModelBlock &          MB,
    const FactoryBlock &        factory_block);
  ~Model   ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  bool processParams ();
  bool processInstanceParams ();
    virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

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

private:

  // Input Parameters
  double vlo;
  double vhi;
  double vref;
  double clo;
  double chi;
  double cload;
  double rload;
  double s0rlo;
  double s0rhi;
  double s0tsw;
  double s0vlo;
  double s0vhi;
  double s1rlo;
  double s1rhi;
  double s1tsw;
  double s1vlo;
  double s1vhi;
  double delay;

  // Dependent Parameters

  double gload;
};


//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Genie Hsieh, SNL, Parallel Computational Sciences
// Creation Date : 10/23/12
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
    friend class Instance;
    friend class Model;

  public:
    Master(
    std::vector< std::pair<std::string,double> > & parNames,
      const Configuration &       configuration,
    const FactoryBlock &        factory_block,
    const SolverState & ss1,
    const DeviceOptions & do1)
      : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  // Genie 110812. For now override nothing and use the template functions of the following
  // from DeviceMaster.h
  /*virtual bool updateState (double * solVec, double * staVec) ;
    virtual bool updateSecondaryState (double * staDeriv);

    // load functions:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
  */

};

void registerDevice();

} // namespace Digital
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Digital::Instance N_DEV_DigitalInstance;
typedef Xyce::Device::Digital::Model N_DEV_DigitalModel;
typedef Xyce::Device::Digital::Master N_DEV_DigitalMaster;

#endif
