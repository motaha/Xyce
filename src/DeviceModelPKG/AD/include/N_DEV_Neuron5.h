//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_Neuron5.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/02/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron5_h
#define Xyce_N_DEV_Neuron5_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#ifdef HAVE_MATH_H
#include <math.h>
#endif

namespace Xyce {
namespace Device {
namespace Neuron5 {

// ---------- Forward Declarations ----------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 Neuron device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the NeuronInstance
//                 class for a more detailed explanation.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;

public:
  static vector< vector<int> > jacStamp;
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
  bool loadDeviceMask();
  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

  bool processParams (string param = "");
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars ();
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
  // iterator reference to the Neuron model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // derrived quantities computed in updateIntermediateVars
  // and used in the load functions
  double kcl1Fvalue, kcl1Qvalue;
  double kcl2Fvalue, kcl2Qvalue;
  double nEquFvalue, nEquQvalue;
  double mEquFvalue, mEquQvalue;
  double hEquFvalue, hEquQvalue;
  double aEquFvalue, aEquQvalue;
  double bEquFvalue, bEquQvalue;
  double M_EquFvalue, M_EquQvalue;
  double H_EquFvalue, H_EquQvalue;
  double cEquFvalue, cEquQvalue;
  double CaEquFvalue, CaEquQvalue;

  double dkcl1F_dV1, dkcl1F_dV2, dkcl1F_dn, dkcl1F_dm, dkcl1F_dh,
    dkcl1F_da, dkcl1F_db, dkcl1F_dM, dkcl1F_dH, dkcl1F_dc, dkcl1Q_dV1, dkcl1Q_dV2;
  double dkcl2F_dV1, dkcl2F_dV2, dkcl2F_dn, dkcl2F_dm, dkcl2F_dh,
    dkcl2F_da, dkcl2F_db, dkcl2F_dM, dkcl2F_dH, dkcl2F_dc, dkcl2Q_dV1, dkcl2Q_dV2;
  double dnF_dV1, dnF_dn, dnQ_dn;
  double dmF_dV1, dmF_dm, dmQ_dm;
  double dhF_dV1, dhF_dh, dhQ_dh;
  double daF_dV1, daF_da, daQ_da;
  double dbF_dV1, dbF_db, dbQ_db;
  double dMF_dV1, dMF_dM, dMQ_dM;
  double dHF_dV1, dHF_dH, dHQ_dH;
  double dcF_dV1, dcF_dc, dcF_dCa, dcQ_dc;
  double dCaF_dV1, dCaF_dV2, dCaF_dM, dCaF_dH, dCaF_dCa, dCaQ_dCa;

  // state variables
  double potassiumCurrent;
  double sodiumCurrent;

  // local state indices (offsets)
  int li_KCurrentState;
  int li_NaCurrentState;

  // local solution indices (offsets)
  int li_Pos;      // local index to positive node on this device
  int li_Neg;      // local index to negative node on this device
  int li_nPro;     // local index to n promoter value (Na current)
  int li_mPro;     // local index to m promoter value (K current)
  int li_hPro;     // local index to h promoter value (K current)
  int li_aPro;     // local index to
  int li_bPro;     // local index
  int li_M_Pro;    // local index
  int li_H_Pro;    // local index
  int li_cPro;     // local index
  int li_CaPro;    // local index

  // Matrix equation index variables:

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int APosEquNNodeOffset;
  int APosEquMNodeOffset;
  int APosEquHNodeOffset;
  int APosEquANodeOffset;
  int APosEquBNodeOffset;
  int APosEquM_NodeOffset;
  int APosEquH_NodeOffset;
  int APosEquCNodeOffset;

  int ANegEquPosNodeOffset;
  int ANegEquNegNodeOffset;
  int ANegEquNNodeOffset;
  int ANegEquMNodeOffset;
  int ANegEquHNodeOffset;
  int ANegEquANodeOffset;
  int ANegEquBNodeOffset;
  int ANegEquM_NodeOffset;
  int ANegEquH_NodeOffset;
  int ANegEquCNodeOffset;

  int ANEquPosNodeOffset;
  int ANEquNNodeOffset;

  int AMEquPosNodeOffset;
  int AMEquMNodeOffset;

  int AHEquPosNodeOffset;
  int AHEquHNodeOffset;

  int AAEquPosNodeOffset;
  int AAEquANodeOffset;

  int ABEquPosNodeOffset;
  int ABEquBNodeOffset;

  int AM_EquPosNodeOffset;
  int AM_EquM_NodeOffset;

  int AH_EquPosNodeOffset;
  int AH_EquH_NodeOffset;

  int ACEquPosNodeOffset;
  int ACEquCNodeOffset;
  int ACEquCaNodeOffset;

  int ACaEquPosNodeOffset;
  int ACaEquNegNodeOffset;
  int ACaEquM_NodeOffset;
  int ACaEquH_NodeOffset;
  int ACaEquCaNodeOffset;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/02/08
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

  bool processParams (string param = "");
  bool processInstanceParams (string param = "");

private:

  // parameter variables
  double cMem;     // membrane capacitance
  double gMem;     // membrane conductance
  double vRest;    // resting potential
  double eNa;      // sodium rest potential
  double gNa;      // sodium base conductance
  double eK;       // potassium rest potential
  double gK;       // potassium base conductance
  double eA;       // a-current rest potential
  double gA;       // a-current base conductance
  double eCa;      // Calcium rest potential
  double gCa;      // Calcium base conductance
  double eKCa;     // potassium-calcium rest potential
  double gKCa;     // potassium-calcium base conductance
  double CaInit;  // initial intra-cellular calcium concentration
  double CaGamma;  // calcium current to concentration multiplier
  double CaTau;    // calcium removal time constant

  // flags that parameters were given
  bool cMemGiven;
  bool gMemGiven;
  bool vRestGiven;
  bool eNaGiven;
  bool gNaGiven;
  bool eKGiven;
  bool gKGiven;
  bool eAGiven;
  bool gAGiven;
  bool eCaGiven;
  bool gCaGiven;
  bool eKCaGiven;
  bool gKCaGiven;
  bool CaInitGiven;
  bool CaGammaGiven;
  bool CaTauGiven;


public:
  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

private:
  vector<Instance*> instanceContainer;
};

} // namespace Neuron5
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Neuron5::Instance N_DEV_NeuronInstance5;
typedef Xyce::Device::Neuron5::Model N_DEV_NeuronModel5;

#endif
