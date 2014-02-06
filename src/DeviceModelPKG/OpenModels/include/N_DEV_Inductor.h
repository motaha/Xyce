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
// Filename       : $RCSfile: N_DEV_Inductor.h,v $
//
// Purpose        : Inductor classes.
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
// Revision Number: $Revision: 1.111.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Inductor_h
#define Xyce_N_DEV_Inductor_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Inductor {

// ---------- Forward Declarations ----------
class Model;
class Instance;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 inductor device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the InductorInstance
//                 class for a more detailed explanation.
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

    Instance(InstanceBlock & IB,
             Model & Iiter,
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

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");
    bool updateTemperature(const double & temp_tmp);

    bool updateIntermediateVars () { return true; };
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

    void setupPointers();

  public:
    // iterator reference to the inductor model which owns this instance.
    // Getters and setters
    Model &getModel() {
      return model_;
    }

  private:
    static vector< vector<int> > jacStamp_BASE;


    Model &       model_;         //< Owning model

    // parameter variables
    double L;  // User specified inductance.
    double IC; // Initial condition: initial, time-zero inductor current(A)
    bool ICGiven;
    double baseL;     // the baseline inductance before temperature effects
    double temp;      // temperature of this instance
    bool tempGiven;

    // Genie 121412. temperature dependence parameters
    // these can override values specified in the model
    double tempCoeff1;   // first order temperature coeff.
    double tempCoeff2;   // second order temperature coeff.

    // flags used to tell if the user has specified one of these values
    // on the command line.
    bool tempCoeff1Given;
    bool tempCoeff2Given;

    // state variables
    double f0; // most recent value for the  flux through the inductor.

    // local state indices (offsets)
    int li_fstate;

    // local solution indices (offsets)
    int li_Pos;
    int li_Neg;
    int li_Bra;

    // Matrix equation index variables:
    vector<int> xLBraVar_J;
    vector<int> li_LBra;

    int ABraEquLBraVar_I; // Row index for the branch current
    // contribution of inductors this instance
    // is coupled to.

    // Offset variables for all of the above index variables.
    int ABraEquPosNodeOffset; // Offset, pos. node voltage contribution,
    // branch current equ.

    int ABraEquNegNodeOffset; // Offset, neg. node voltage contribution,
    // branch current equ.

    int ABraEquBraVarOffset;  // Offset, branch current variable
    // contribution, branch current equation.

    int APosEquBraVarOffset;  // Offset, branch current variable
    // contribution, KCL equation of the pos node

    int ANegEquBraVarOffset;  // Offset, branch current variable
    // contribution, KCL equation of the neg node

    int AEPosEquEBraVarOffset;

    int AENegEquEBraVarOffset;

    int AEBraEquEPosNodeOffset;

    int AEBraEquENegNodeOffset;

    int AEBraEquLNegNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // Pointer variables for the Jacobian matrices
    double * fPosEquBraVarPtr;
    double * fNegEquBraVarPtr;
    double * fBraEquPosNodePtr;
    double * fBraEquNegNodePtr;
    double * fBraEquBraVarPtr;
    double * qBraEquBraVarPtr;
#endif

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
    virtual bool processParams (string param = "");
    virtual bool processInstanceParams (string param = "");


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

    double L;  // User specified inductance.
    double IC; // Initial condition: initial, time-zero inductor current(A)
    double tempCoeff1;     // first order temperature coeff.
    double tempCoeff2;     // second order temperature coeff.
    double baseL;
    double tnom;
    bool tc1Given;
    bool tc2Given;
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
    {

    }

    virtual bool updateState (double * solVec, double * staVec, double * stoVec);
    virtual bool updateSecondaryState (double * staDeriv, double * stoVec);

    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
};

} // namespace Inductor
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Inductor::Instance N_DEV_InductorInstance;
typedef Xyce::Device::Inductor::Model N_DEV_InductorModel;
typedef Xyce::Device::Inductor::Master N_DEV_InductorMaster;

#endif
