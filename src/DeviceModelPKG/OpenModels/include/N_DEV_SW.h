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
// Filename       : $RCSfile: N_DEV_SW.h,v $
//
// Purpose        : Switch base classes
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
// Revision Number: $Revision: 1.77.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_SW_h
#define Xyce_N_DEV_SW_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

class N_UTL_Expression;

namespace Xyce {
namespace Device {
namespace SW {

// ----------   Fwd Declarations  -------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
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

    Instance(InstanceBlock &IB,
                     Model & SWiter,
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
    void registerStoreLIDs( const vector<int> & stoLIDVecRef );

    map<int,string> & getStoreNameMap();
    
    bool processParams (string param = "");

  const std::vector<std::string> & getDepSolnVars();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool updateIntermediateVars ();
    bool updatePrimaryState ();
    bool updateSecondaryState ();

    // load functions, residual:
    bool loadDAEQVector () {return true;}
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx () {return true;}
    bool loadDAEdFdx ();

    void setupPointers();

  public:
    // iterator reference to the switch model which owns this instance
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

    N_UTL_Expression * Exp_ptr;
    int            expNumVars;
    int            expBaseVar;
    int            expNumDdt;
    list<string>   evnList;

    vector<double> expVarDerivs;
    vector<double> myVarVals;
    vector<double> ddtVals;
    double         expVal;


    // user specified parameters
    double R;     // resistance (ohms)
    double CONTROL;   // Value of control expression
    bool ON, OFF;      // whether switch is on or off initially

    // derived parameters
    double G;     // conductance (1.0/ohms)
    double dGdI,dGdV;

    // places to store node voltages
    double v_pos;
    double v_neg;

    double LeadCurrent;

    // double for current state of switch
    double SW_STATE;

    // and a state variable to save the SW_STATE
    double switch_state;

    vector<int>    li_ddt;

    int li_switch_state;

    // local indices (offsets)
    int li_Pos;
    int li_Neg;
    
    // store vector location for device lead current
    int li_store_dev_i;

    // Offset variables corresponding to the above declared indices.
    int APosEquPosNodeOffset;
    int APosEquNegNodeOffset;
    int ANegEquPosNodeOffset;
    int ANegEquNegNodeOffset;

    // Offsets into the control nodes
    vector<int> APosEquControlNodeOffset;
    vector<int> ANegEquControlNodeOffset;

    // Ptr variables corresponding to the above declared indices.
    double * fPosEquPosNodePtr;
    double * fPosEquNegNodePtr;
    double * fNegEquPosNodePtr;
    double * fNegEquNegNodePtr;

    // Ptrs into the control nodes
    vector<double *> fPosEquControlNodePtr;
    vector<double *> fNegEquControlNodePtr;

    vector< vector<int> > jacStamp;
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

    Model(const ModelBlock &MB,
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

    int dtype;      // device type: 1=SWITCH, 2=ISWITCH, 3=VSWITCH
    double VON;
    double VOFF;
    double ION;
    double IOFF;
    double RON;
    double ROFF;
    double ON;
    double OFF;
    double dInv;    // the inverse of (ON-OFF) or 1e-12, if too small.
    double Lm;      // log mean of resistances
    double Lr;      // log ratio of resistor values
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
class Master : public Xyce::Device::DeviceTemplate<Model, Instance>
{
  public:
    Master (
      const std::string &dn,
      const std::string &cn,
      const std::string &dmName,
           LinearDevice linearDev,
           SolverState & ss1,
           DeviceOptions & do1)
      : Xyce::Device::DeviceTemplate<Model, Instance>(
           dn, cn, dmName, linearDev, ss1, do1)
    {

    }

    virtual bool updateState (double * solVec, double * staVec, double * stoVec);
    virtual bool updateSecondaryState (double * staDeriv, double * stoVec);

    // load functions, residual:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);

    // load functions, Jacobian:
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

    friend class Instance;
    friend class Model;
};

} // namespace SW
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::SW::Instance N_DEV_SWInstance;
typedef Xyce::Device::SW::Model N_DEV_SWModel;
typedef Xyce::Device::SW::Master N_DEV_SWMaster;

#endif
