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
// Filename       : $RCSfile: N_DEV_Neuron6.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 06/10/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.17.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron6_h
#define Xyce_N_DEV_Neuron6_h

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_MembraneModel.h>

#ifdef HAVE_MATH_H
#include <math.h>
#endif

namespace Xyce {
namespace Device {
namespace Neuron6 {

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
// Creation Date : 06/10/09
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

  // model level parameters that can be overridden at the instance level
  double rInt;     // intracellular resistivity
  double radius;   // Segment radius
  double length;   // cable length (segment length = length/nSeg)
  double segArea;  // segment area (derrived from radius, length and nSeg)
  int    nSeg;     // number of segments
  bool rIntGiven;
  bool radiusGiven;
  bool lengthGiven;
  bool nSegGiven;

  // conductance between segments -- calculated from radius, rInt and length and number of segments
  double gSeg;

  int numIntVarsPerSegment;
  int numStateVarsPerSegment;

  // storage for local ID's of internal vars and jacobian offsets
  vector< int > li_internalVars;
  vector< vector< int > > jacobianOffsets;

  // derrived quantities computed in updateIntermediateVars
  // and used in the load functions (no q terms on the external nodes)
  double kcl1Fvalue;
  double kcl2Fvalue;
  // internal segments
  vector<double> segFvalue;
  vector<double> segQvalue;
  vector<double> segNEquFvalue, segNEquQvalue;
  vector<double> segMEquFvalue, segMEquQvalue;
  vector<double> segHEquFvalue, segHEquQvalue;
  vector<double> segAEquFvalue, segAEquQvalue;
  vector<double> segBEquFvalue, segBEquQvalue;
  vector<double> segM_EquFvalue, segM_EquQvalue;
  vector<double> segH_EquFvalue, segH_EquQvalue;
  vector<double> segCEquFvalue, segCEquQvalue;
  vector<double> segCaEquFvalue, segCaEquQvalue;

  // jacobian terms
  double dkcl1F_dVin, dkcl1F_dVs0;
  double dkcl2F_dVout, dkcl2F_dVsn;
  // internal equations
  vector<double> segF_dVp, segF_dV, segF_dVn, segF_dn, segF_dm, segF_dh, segF_da, segF_db, segF_dM, segF_dH, segF_dc;
  vector<double> segQ_dV;
  vector<double> dnF_dV, dnF_dn, dnQ_dn;
  vector<double> dmF_dV, dmF_dm, dmQ_dm;
  vector<double> dhF_dV, dhF_dh, dhQ_dh;
  vector<double> daF_dV, daF_da, daQ_da;
  vector<double> dbF_dV, dbF_db, dbQ_db;
  vector<double> dMF_dV, dMF_dM, dMQ_dM;
  vector<double> dHF_dV, dHF_dH, dHQ_dH;
  vector<double> dcF_dV, dcF_dc, dcF_dCa, dcQ_dc;
  vector<double> dCaF_dV, dCaF_dM, dCaF_dH, dCaF_dCa, dCaQ_dCa;

  // state variables
  vector<double> potassiumCurrent;
  vector<double> sodiumCurrent;

  // local state indices (offsets)
  vector<int> li_KCurrentState;
  vector<int> li_NaCurrentState;

  // local solution indices (offsets)
  int li_Pos;      // local index to positive node on this device
  int li_Neg;      // local index to negative node on this device
  // local solution indices for internal vars (variable number of these)
  vector<int> li_Vol;      // local index to segment voltage
  vector<int> li_nPro;     // local index to n promoter value (Na current)
  vector<int> li_mPro;     // local index to m promoter value (K current)
  vector<int> li_hPro;     // local index to h promoter value (K current)
  vector<int> li_aPro;     // local index to a promoter value
  vector<int> li_bPro;     // local index to a promoter value
  vector<int> li_MPro;     // local index to a promoter value
  vector<int> li_HPro;     // local index to a promoter value
  vector<int> li_cPro;     // local index to a promoter value
  vector<int> li_CaPro;     // local index to a promoter value

  // Matrix equation index variables:

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset, APosEquNextNodeOffset;
  int ANegEquNegNodeOffset, ANegEquLastNodeOffset;
  vector<int> SegVEqnVpreOffset;
  vector<int> SegVEqnVsegOffset;
  vector<int> SegVEqnVnexOffset;
  vector<int> SegVEqnNOffset;
  vector<int> SegVEqnMOffset;
  vector<int> SegVEqnHOffset;
  vector<int> SegVEqnAOffset;
  vector<int> SegVEqnBOffset;
  vector<int> SegVEqnM_Offset;
  vector<int> SegVEqnH_Offset;
  vector<int> SegVEqnCOffset;
  vector<int> NEquVNodeOffset;
  vector<int> NEquNNodeOffset;
  vector<int> MEquVNodeOffset;
  vector<int> MEquMNodeOffset;
  vector<int> HEquVNodeOffset;
  vector<int> HEquHNodeOffset;
  vector<int> AEquVNodeOffset;
  vector<int> AEquANodeOffset;
  vector<int> BEquVNodeOffset;
  vector<int> BEquBNodeOffset;
  vector<int> M_EquVNodeOffset;
  vector<int> M_EquM_NodeOffset;
  vector<int> H_EquVNodeOffset;
  vector<int> H_EquH_NodeOffset;
  vector<int> CEquVNodeOffset;
  vector<int> CEquCNodeOffset;
  vector<int> CEquCaNodeOffset;
  vector<int> CaEquVNodeOffset;
  vector<int> CaEquM_NodeOffset;
  vector<int> CaEquH_NodeOffset;
  vector<int> CaEquCaNodeOffset;

  // maps to track the appropriate jacobian offsets for each segment's previous, current, and next segment
  map <int, int> prevMap;
  map<int, int> segMap;
  map<int, int> nextMap;

  vector< vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  friend class ParametricData<Model>;   typedef std::vector<Instance *> InstanceVector;


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
  double rInt;     // intracellular resistivity
  double radius;   // Segment radius
  double length;   // cable length (segment length = length/nSeg)
  string ionChannelModel; // what model will be used for the ion channels
  int    nSeg;     // number of segments

  // Value of current expression for user-defined membranem model
  double I;

  // these are vectors of strings to allow the user to specify independant vars and
  // equations for a given membrane
  vector<string> membraneCurrentEqus;
  vector<string> membraneIndpVars;
  vector<string> membraneIndpFEqus;
  vector<string> membraneIndpQEqus;
  vector<string> membraneFunctions;
  vector<string> membraneParameters;

  // flags that parameters were given
  bool rIntGiven;
  bool radiusGiven;
  bool lengthGiven;
  bool ionChannelModelGiven;
  bool nSegGiven;

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

  bool membraneCurrentEqusGiven;
  bool membraneIndpVarsGiven;
  bool membraneIndpFEqusGiven;
  bool membraneIndpQEqusGiven;
  bool membraneFunctionsGiven;
  bool membraneParametersGiven;

  // these are meta flags.  If one component of the a sodium current is given
  // then all of the required equaitons will be used.  Otherwise they are off
  // Or, if hodgenHuxleyOn_ is true, then all of the H odgenHuxley equations are loaded.
  // by default all of these are off.
  bool hodgenHuxleyOn_;
  bool ConnorStevensOn_;
  bool sodiumOn_;
  bool potassiumOn_;
  bool aCurrentOn_;
  bool calciumOn_;
  RefCountPtr< MembraneModel > membraneModel_;


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

} // namespace Neuron6
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Neuron6::Instance N_DEV_NeuronInstance6;
typedef Xyce::Device::Neuron6::Model N_DEV_NeuronModel6;

#endif
