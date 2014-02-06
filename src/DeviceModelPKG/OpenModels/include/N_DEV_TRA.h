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
// Filename       : $RCSfile: N_DEV_TRA.h,v $
//
// Purpose        : Transmission line.
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
// Revision Number: $Revision: 1.74.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_TRA_h
#define Xyce_N_DEV_TRA_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {
namespace TRA {

// ---------- Forward Declarations ----------
class Model;
class History;

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

  public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }


    Instance(InstanceBlock &IB,
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
    void registerStoreLIDs( const vector<int> & st0LIDVecRef );

    map<int,string> & getIntNameMap ();
    map<int,string> & getStoreNameMap ();

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

    bool getInstanceBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes);
    void acceptStep();

    double getMaxTimeStepSize();

    DeviceState * getInternalState();
    bool setInternalState( const DeviceState & state );

  private:
    void pruneHistory(double t);
    void InterpV1V2FromHistory(double t, double * v1p , double * v2p);

  public:
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:
  static vector< vector<int> > jacStamp;

  Model &       model_;         //< Owning model

    double Z0;   // Characteristic impedence
    double ZO;
    double G0;   // Conductance
    double td;   // Time delay  (= NL/freq if given that way)
    double freq; // frequency
    double NL;   // Normalized length
    // Flags
    // not supporting initial condition just yet... this is An Issue, I think
    // double IC_V12, IC_V34, IC_I1, IC_I2;

    bool DCMODE;
    // Matrix equation index variables
    // solution vector indicies
#if 0
    int xPos1Var_J;   // pos and negative port nodes
    int xNeg1Var_J;   // .....
    int xInt1Var_J;   // Internal variable for port 1
    int xIbr1Var_J;   // Branch current
    int xPos2Var_J;   // Duplicate above for second port
    int xNeg2Var_J;
    int xInt2Var_J;
    int xIbr2Var_J;
#endif
    // KCL equation indices
#if 0
    int bPos1Equ_I;
    int bNeg1Equ_I;
    int bInt1Equ_I;
    int bIbr1Equ_I;
    int bPos2Equ_I;
    int bNeg2Equ_I;
    int bInt2Equ_I;
    int bIbr2Equ_I;
#endif
    // local indices (offsets)
    int li_Pos1;
    int li_Neg1;
    int li_Int1;
    int li_Ibr1;
    int li_Pos2;
    int li_Neg2;
    int li_Int2;
    int li_Ibr2;
    
    // indices into store vec for lead currents if needed
    int li_store_dev_i1;
    int li_store_dev_i2;
    

    // Matrix elements
#if 0
    int APos1EquPos1Node_I;
    int APos1EquPos1Node_J;
    int APos1EquInt1Node_I;
    int APos1EquInt1Node_J;

    int AInt1EquPos1Node_I;
    int AInt1EquPos1Node_J;
    int AInt1EquInt1Node_I;
    int AInt1EquInt1Node_J;
    int AInt1EquIbr1Node_I;
    int AInt1EquIbr1Node_J;

    int ANeg1EquIbr1Node_I;
    int ANeg1EquIbr1Node_J;

    int AIbr1EquInt1Node_I;
    int AIbr1EquInt1Node_J;
    int AIbr1EquNeg1Node_I;
    int AIbr1EquNeg1Node_J;

    int APos2EquPos2Node_I;
    int APos2EquPos2Node_J;
    int APos2EquInt2Node_I;
    int APos2EquInt2Node_J;

    int AInt2EquPos2Node_I;
    int AInt2EquPos2Node_J;
    int AInt2EquInt2Node_I;
    int AInt2EquInt2Node_J;
    int AInt2EquIbr2Node_I;
    int AInt2EquIbr2Node_J;

    int ANeg2EquIbr2Node_I;
    int ANeg2EquIbr2Node_J;

    int AIbr2EquInt2Node_I;
    int AIbr2EquInt2Node_J;
    int AIbr2EquNeg2Node_I;
    int AIbr2EquNeg2Node_J;

    // for DC simulations these 6 pairs get filled because v1 and v2 reduce
    // to V1=Ibr2*Z0 and V2=Ibr1*Z0 and no time delay
    int AIbr1EquPos2Node_I;
    int AIbr1EquPos2Node_J;
    int AIbr1EquNeg2Node_I;
    int AIbr1EquNeg2Node_J;
    int AIbr1EquIbr2Node_I;
    int AIbr1EquIbr2Node_J;

    int AIbr2EquPos1Node_I;
    int AIbr2EquPos1Node_J;
    int AIbr2EquNeg1Node_I;
    int AIbr2EquNeg1Node_J;
    int AIbr2EquIbr1Node_I;
    int AIbr2EquIbr1Node_J;
#endif
    // Matrix equation offset variables
    int APos1EquPos1NodeOffset;
    int APos1EquInt1NodeOffset;
    int AInt1EquPos1NodeOffset;
    int AInt1EquInt1NodeOffset;
    int AInt1EquIbr1NodeOffset;
    int ANeg1EquIbr1NodeOffset;
    int AIbr1EquInt1NodeOffset;
    int AIbr1EquNeg1NodeOffset;
    // for DC simulations these 6 pairs get filled because v1 and v2 reduce
    // to V1=Ibr2*Z0 and V2=Ibr1*Z0 and no time delay
    int APos2EquPos2NodeOffset;
    int APos2EquInt2NodeOffset;
    int AInt2EquPos2NodeOffset;
    int AInt2EquInt2NodeOffset;
    int AInt2EquIbr2NodeOffset;
    int ANeg2EquIbr2NodeOffset;
    int AIbr2EquInt2NodeOffset;
    int AIbr2EquNeg2NodeOffset;
    int AIbr1EquPos2NodeOffset;
    int AIbr1EquNeg2NodeOffset;
    int AIbr1EquIbr2NodeOffset;
    int AIbr2EquPos1NodeOffset;
    int AIbr2EquNeg1NodeOffset;
    int AIbr2EquIbr1NodeOffset;

    double Vpos1,Vpos2,Vneg1,Vneg2,Vint1,Vint2,Ibr1,Ibr2; // solution vars
    double last_t;     // "primary" state variables
    double v1;    //
    double v2;    //
    bool first_BP_call_done;

    vector<History> history;
    int newtonIterOld;
    double timeOld;

    bool newBreakPoint;
    double newBreakPointTime;
};

//-----------------------------------------------------------------------------
// Class         : History
// Purpose       : Provide a structure to save internal state history
// Special Notes :
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 6/14/2001
//-----------------------------------------------------------------------------
class History
{
  friend class Instance;
 public:
  History();
  History(const History &right);
  History(double t, double v1, double v2);
  ~History();
  inline bool operator<(const double &test_t) const;

 private:
  double t;
  double v1;
  double v2;
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

};

//-----------------------------------------------------------------------------
// Function      : History::operator<
// Purpose       : compare used in the lower_bound operation.  Returns
//                 true if the time in the history object is less than
//                 the given time (double)
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 7/27/2005
//-----------------------------------------------------------------------------
inline bool History::operator<(const double &test_t) const
{
  return (t < test_t);
}

} // namespace TRA
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::TRA::Instance N_DEV_TRAInstance;
typedef Xyce::Device::TRA::Model N_DEV_TRAModel;
typedef Xyce::Device::TRA::History N_DEV_TRAHistory;

#endif
