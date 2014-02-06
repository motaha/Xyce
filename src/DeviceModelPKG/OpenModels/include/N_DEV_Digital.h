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
// Revision Number: $Revision: 1.36.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Digital_h
#define Xyce_N_DEV_Digital_h

// ----------   Xyce Includes   ----------
//Genie 110812
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Digital {

// ---------- Forward Declarations ----------
class Model;
class Instance; //Genie 110812

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
  friend class Master;

  enum gType {NOT, AND, NAND, OR, NOR, ADD, XOR, NXOR, DFF};

public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Instance(InstanceBlock & IB,
           Model & Diter,
           MatrixLoadData & mlData1,
           SolverState &ss1,
           ExternData  &ed1,
           DeviceOptions & do1);


  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  // Additional Public Declarations
  void registerLIDs( const vector<int> & intLIDVecRef,
                     const vector<int> & extLIDVecRef );
  void registerStateLIDs( const vector<int> & staLIDVecRef );

  map<int,string> & getIntNameMap ();

  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

  bool processParams (string param = "");

  bool updateIntermediateVars () { return true; };
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool getInstanceBreakPoints (vector<N_UTL_BreakPoint> &);

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

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

  vector<double> qlo;      // charge in the capacitor
  vector<double> ilo;      // current throught the capacitor
  vector<double> vcaplo;   // voltage drop across capacitor
  vector<double> qhi;
  vector<double> ihi;
  vector<double> vcaphi;
  vector<double> qref;
  vector<double> iref;
  vector<double> vcapref;

  vector<double> rilo;
  vector<double> rihi;
  vector<double> riref;
  vector<double> currentOut;
  vector<double> currentIn;

  vector<double> glo;
  vector<double> ghi;

  vector<double> qInp;      // charge in the capacitor
  vector<double> iInp;      // current throught the capacitor
  vector<double> vcapInp;   // voltage drop across capacitor

  vector<double> currentInp;

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
  vector<int> li_Inp;
  vector<int> li_Out;

  // Input state vars
  vector<int> li_currentStateInp;
  vector<int> li_transitionTimeInp;
  vector<int> li_QinpState;
  vector<int> li_IinpState;

  // Output state vars
  vector<int> li_currentStateOut;
  vector<int> li_transitionTimeOut;
  vector<int> li_QloState;
  vector<int> li_IloState;
  vector<int> li_QhiState;
  vector<int> li_IhiState;

  vector<bool> inpL;
  vector<double> iTime;
  vector<bool> outL;
  vector<double> oTime;

  double breakTime;

  // Offsets for Jacobian
  int row_Lo;
  int row_Hi;
  int row_Ref;
  vector< vector<int> > li_jac_Ref;
  vector< vector<int> > li_jac_Lo;
  vector< vector<int> > li_jac_Hi;

  vector< vector<int> > jacStamp;

  //Genie 110812. change state var
  //vector<bool> changeState;
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
  friend class Master;

public:
  static ParametricData<Model> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Model    (const ModelBlock & MB,
            SolverState & ss1,
            DeviceOptions & do1);
  ~Model   ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  bool processParams (string param = "");
  bool processInstanceParams (string param = "");
  virtual std::ostream &printOutInstances(std::ostream &os) const;

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
class Master : public Xyce::Device::DeviceTemplate<Model, Instance>
{
public:
  Master (
    const std::string &dn,
    const std::string &cn,
    const std::string &dmName,
    vector< pair<string,double> > & parNames,
    LinearDevice linearDev,
    SolverState & ss1,
    DeviceOptions & do1)
    : Xyce::Device::DeviceTemplate<Model, Instance>(
      dn, cn, dmName, linearDev, ss1, do1)
  {

  }

  // Genie 110812. For now override nothing and use the template functions of the following
  // from DeviceTemplate.h
  /*virtual bool updateState (double * solVec, double * staVec) ;
    virtual bool updateSecondaryState (double * staDeriv);

    // load functions:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
  */

  friend class Instance;
  friend class Model;
};

} // namespace Digital
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Digital::Instance N_DEV_DigitalInstance;
typedef Xyce::Device::Digital::Model N_DEV_DigitalModel;
typedef Xyce::Device::Digital::Master N_DEV_DigitalMaster;

#endif
