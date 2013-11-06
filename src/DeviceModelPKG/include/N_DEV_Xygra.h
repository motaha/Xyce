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
// Filename       : $RCSfile: N_DEV_Xygra.h,v $
//
// Purpose        : Xygra classes.
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
// Revision Number: $Revision: 1.10.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Xygra_h
#define Xyce_N_DEV_Xygra_h

// ----------   Xyce Includes   ----------
#include <Sacado.hpp>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : XygraCoilData
// Purpose       : This is class is a CompositeParameter type for managing
//                 coil vector-composite data
// Special Notes :
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
class XygraCoilData : public CompositeParam
{
  friend class ParametricData<XygraCoilData>;

public:
  static ParametricData<XygraCoilData> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  XygraCoilData();

  void processParams();
#ifdef Xyce_DEBUG_DEVICE
  friend ostream & operator<<(ostream & os, const XygraCoilData & xcd);
#endif

private:

private:
  string name;
  int numWindings;

public:
  string getName() const { return name;};
  int getNumWindings() const { return numWindings;};
};


namespace Xygra {

// ---------- Forward Declarations ----------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 Xygra device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the ResistorInstance
//                 class for a more detailed explanation.
// Special Notes :
// Creator       : Tom Russo
// Creation Date : 8/18/08
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;

  typedef Sacado::Fad::DFad<double> XygraFadType;

public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Instance(
    InstanceBlock &            IB,
    Model & Miter,
    MatrixLoadData &           mlData1,
    SolverState &              ss1,
    ExternData  &              ed1,
    DeviceOptions &            do1);

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

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  bool setIC ();

  bool getVoltages(vector<double> &voltageValues);
  bool setConductances(const vector< vector<double> > &conductanceMatrix);
  bool setK(const vector< vector<double> > &kMatrix, const double t=0);
  bool setSources(const vector<double> &sourceVector,const double t=0);
  int getNumNodes();
  int getNumWindings();
  void getCoilWindings(vector<int> &coilWindings);
  void getCoilNames(std::vector<std::string> &coilNames);

  void varTypes( vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  CompositeParam *constructComposite (string &, string &);

protected:
private:
  void setupJacStamp_();
  void interpolateSandK_();

public:
  // iterator reference to the Xygra model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  map <string,CompositeParam *> coilDataMap;

private:
  // parameter variables

  // state variables
  // This device has no state

  // local state indices (offsets)
  // This device has no state

  // local solution indices (offsets)
  // This device uses an array of li_ values instead of individually named
  // variables.
  vector<int> li_Nodes_;

  // Matrix equation index variables:

  // Offset variables.  Again, this device uses an array instead of
  // discrete variables.
  // A_Equ_NodeOffests[equation][node] is the offset for node in
  // equation
  vector< vector<int> > A_Equ_NodeOffsets_;

  vector< vector<int> > jacStamp_;

  // These guys hold the Alegra input
  vector< vector<double> > theConductanceMatrix_;
  vector< vector<double> > theKMatrix_;
  vector< vector<double> > k0_;
  vector< vector<double> > k1_;
  vector<double> theSourceVector_;
  vector<double> s0_;
  vector<double> s1_;
  // times that (s0,k0) and (s1,k1) apply to.
  double t0_;
  double t1_;

  // For vector composite:
  vector<XygraCoilData*> coilDataVec;
  // total number of coils
  int nCoils;
  // number of windngs in each coil
  vector <int> nWindings;
  // names of each coil
  vector <string> coilNames;
  // sum over coils of number of windings per coil
  int totalNumWindings;
  // offsets into global node array of start of each coil's external vars
  vector <int> coilExtStart;
  // offsets into global node array of start of each coil's Internal vars
  vector <int> coilIntStart;
  // vector of pairs of nodes (pos,neg) for every winding
  vector <pair<int,int> > windingNodes;

  // For computation of RHS/F vector and jacobian/dFdX
  // We copy solution vars here so we can differentiate w.r.t them.
  vector<XygraFadType> solutionVars;
  // This is the vector of winding dv's
  vector<XygraFadType> dV;
  // This is the vector of winding currents
  vector<XygraFadType> windingCurrents;
  // and finally the vector of contributions into F:
  vector<XygraFadType> fContributions;
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

  static int numOrig;
  static int numSer;

  // Additional Implementation Declarations
};

//----------------------------------------------------------------------------
// Function      : Instance::getNumNodes
// Purpose       : Return the number of nodes in a given instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/2008
//----------------------------------------------------------------------------
inline int Instance::getNumNodes()
{
  return numExtVars+numIntVars;
}
//----------------------------------------------------------------------------
// Function      : Instance::getNumWindings()
// Purpose       : Return the number of windings in a given instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/2008
//----------------------------------------------------------------------------
inline int Instance::getNumWindings()
{
  return totalNumWindings;
}

} // namespace Resistor
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Xygra::Instance N_DEV_XygraInstance;
typedef Xyce::Device::Xygra::Model N_DEV_XygraModel;
typedef Xyce::Device::XygraCoilData N_DEV_XygraCoilData;

#endif
