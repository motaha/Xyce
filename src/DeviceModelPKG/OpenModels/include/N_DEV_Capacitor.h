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
// Filename       : $RCSfile: N_DEV_Capacitor.h,v $
//
// Purpose        :
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
// Revision Number: $Revision: 1.104.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Capacitor_h
#define Xyce_N_DEV_Capacitor_h

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Capacitor {

class Model;
class Instance;

//-----------------------------------------------------------------------------
// Class         : Capacitor::Instance
// Purpose       : This class refers to a single instance of the capacitor
//                 device.  It contains indicies into the matrix equation.
//                 See the comments for the ResistorInstance class for
//                 more details.
//
// Special Notes : A capacitor  will have two circuit nodes.
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

/** 
 * Capacitor instance
 *
 * This class refers to a single instance of the capacitor device.  It
 * contains indicies into the matrix equation.  See the comments for the
 * Resistor::Instance class for more details.
 *
 */
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

    Instance(
      InstanceBlock &   IB,
      Model &           Citer,
      MatrixLoadData &  mlData1,
      SolverState &     ss1,
      ExternData &      ed1,
      DeviceOptions &   do1);

    ~Instance();

  private:
    Instance(const Instance &);
    Instance &operator=(const Instance &);

  public:
    void registerLIDs( const vector<int> & intLIDVecRef, const vector<int> & extLIDVecRef );
    void registerStateLIDs( const vector<int> & staLIDVecRef );
    void registerStoreLIDs( const vector<int> & stoLIDVecRef );

    map<int,string> & getIntNameMap ();
    map<int,string> & getStoreNameMap();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");
    bool updateTemperature(const double & temp_tmp);

    bool updateIntermediateVars () { return true; };
    bool updatePrimaryState ();
    bool updateSecondaryState ();

    bool setIC ();

    // load functions, residual:
    bool loadDAEQVector ();
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx ();
    bool loadDAEdFdx ();

    void setupPointers();

    void varTypes( vector<char> & varTypeVec );

    // Getters and setters
    Model &getModel() {
      return model_;
    }

  private:
    Model &       model_;         //< Owning model

    // Stuff for handling solution-variable-dependent capacitance
    N_UTL_Expression * expPtr;
    int                expNumVars;

    vector<double> expVarDerivs;

    // user-specified parameters:
    double C;    // User specified capacitance. (Farads)
    double IC;   // Optional initial value capacitor voltage (V).

    // These are for the semiconductor capacitor
    double length;    // capacitor length
    double width;     // capacitor width
    double temp;      // temperature of this instance

    // Genie 121412. temperature dependence parameters
    // these can override values specified in the model
    double tempCoeff1;   // first order temperature coeff.
    double tempCoeff2;   // second order temperature coeff.

    // flags used to tell if the user has specified one of these values
    // on the command line.
    bool tempCoeff1Given;
    bool tempCoeff2Given;

    // These are for the age-aware capacitor
    double age;                 ///< age in hours
    double ageCoef;             ///< degradation coeficient.
    double baseCap;             ///< the baseline capacitance before aging

    bool tempGiven;
    bool ICGiven;
    bool solVarDepC;

    // state variables:
    double q0;                  ///< charge in the capacitor
    // now held in the store vector at li_store_dev_i
    double vcap; // voltage drop across capacitor

    //local id's (offsets)
    int li_Pos;
    int li_Neg;
    int li_Bra;                 ///< for the "voltage source" when IC is specified

    int li_QState;

    vector<int> li_dQdXState;
    vector<int> li_dCdXState;
    int li_vcapState;
    int li_capState;

    int li_store_dev_i;

    // Offsets for Jacobian
    int APosEquPosNodeOffset;
    int ANegEquPosNodeOffset;
    int APosEquNegNodeOffset;
    int ANegEquNegNodeOffset;

    // offsets for when C is solution-variable dependent
    vector<int> APosEquDepVarOffsets;
    vector<int> ANegEquDepVarOffsets;

    int ABraEquPosNodeOffset;
    int ABraEquNegNodeOffset;
    int ABraEquBraNodeOffset;
    int APosEquBraNodeOffset;
    int ANegEquBraNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // Pointers for Jacobian
    double * qPosEquPosNodePtr;
    double * qNegEquPosNodePtr;
    double * qPosEquNegNodePtr;
    double * qNegEquNegNodePtr;

    double * fBraEquPosNodePtr;
    double * fBraEquNegNodePtr;
    double * fBraEquBraNodePtr;
    double * fPosEquBraNodePtr;
    double * fNegEquBraNodePtr;

    vector<double *> qPosEquDepVarsPtrs;
    vector<double *> qNegEquDepVarsPtrs;
#endif

    vector< vector<int> > jacStamp;
    vector< vector<int> > jacStamp_IC;
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

    Model(const ModelBlock & MB, SolverState & ss1, DeviceOptions & do1);
    ~Model();

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

    // for the semiconductor capacitor
    double cj;     // junction bottom capacitance
    double cjsw;   // junction sidewall capacitance
    double defWidth; // default width
    double narrow;   // narrowing due to side etching
    double tempCoeff1;   // first order temperature coeff.
    double tempCoeff2;   // second order temperature coeff.
    double baseCap;
    double tnom;

    bool tnomGiven;
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
    {}

    virtual bool updateState (double * solVec, double * staVec, double * stoVec);
    virtual bool updateSecondaryState (double * staDeriv, double * stoVec);

    // load functions:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

    friend class Instance;
    friend class Model;
};

} // namespace Capacitor
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Capacitor::Instance N_DEV_CapacitorInstance;
typedef Xyce::Device::Capacitor::Model N_DEV_CapacitorModel;
typedef Xyce::Device::Capacitor::Master N_DEV_CapacitorMaster;

#endif // Xyce_N_DEV_Capacitor_h

