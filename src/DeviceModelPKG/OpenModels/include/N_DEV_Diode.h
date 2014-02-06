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
// Filename       : $RCSfile: N_DEV_Diode.h,v $
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
// Revision Number: $Revision: 1.106.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Diode_h
#define Xyce_N_DEV_Diode_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>

namespace Xyce {
namespace Device {
namespace Diode {

class Model;
class Instance;

//-----------------------------------------------------------------------------
// Class         : Instance
//
// Purpose       : This class represents a single instance  of the  diode
//                 device.  It mainly contains indices and pointers into
//                 the matrix equation (see the resistor instance class for
//                 more details).
//
//                 The diode will have 1 internal node, in addition to the
//                 nodes which are connected to it.  This is so that  it
//                 can be placed in series with a resistor that represents
//                 the resistance of intrinsic Si.
//
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 02/28/00
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
    void registerStoreLIDs( const vector<int> & stoLIDVecRef);

    map<int,string> & getIntNameMap ();
    std::map<int, std::string> & getStoreNameMap ();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");
    bool updateTemperature ( const double & temp = -999.0 );
    bool lambertWCurrent (double Isat, double Vte, double RS);
    bool lambertWBreakdownCurrent (double Isat, double Vte, double RS);
    bool lambertWLinearReverseBias (double Isat, double Vte, double RS);

    bool updateIntermediateVars ();
    bool updatePrimaryState ();

    // load functions, residual:
    bool loadDAEQVector ();
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx ();
    bool loadDAEdFdx ();

    void setupPointers();

  public:
    // iterator reference to the diode model which owns this instance.
    // Getters and setters
    Model &getModel() {
      return model_;
    }

  private:
    static vector< vector<int> > jacStamp_RS;
    static vector< vector<int> > jacStamp;

    static vector<int> jacMap_RS;
    static vector<int> jacMap;

    static vector< vector<int> > jacMap2_RS;
    static vector< vector<int> > jacMap2;


    Model &       model_;         //< Owning model

    int  off;
    double Area;
    double InitCond;
    double Temp;
    int lambertWFlag;
    bool InitCondGiven;

    double tJctPot;
    double tJctCap;
    double tDepCap;
    double tSatCur;
    double tVcrit;
    double tF1;
    double tBrkdwnV;
    double tSatCurR;
    double tIKF;
    double tRS;
    double tCOND;
    double tIRF;

    double Id;     //diode current
    double Gd;     //diode conductivity
    double Cd;     //depletion capacitance
    double Gcd;    //dep cap conductivity
    double Qd;     //capacitor charge
    double Icd;    //capacitor current
    double Gspr;
    //double LeadCurrent;

    double Vpp;
    double Vp;
    double Vn;
    double Vc;

    double Vd;
    double Vd_old;
    double Vd_orig;

    int newtonIterOld;

    // end of intermediate variables

    // state variables:
    double q0;  // charge in the capacitor
    double i0;  // current throught the capacitor

    //local indices (offsets)
    // int li_QState;

    // for voltage limiting
    int li_storevd;

    // for lead current
    int li_store_dev_i;

    int li_Pos;
    int li_Neg;
    int li_Pri;

    // Matrix equation local offset variables
    int APosEquPosNodeOffset;
    int APosEquPriNodeOffset;
    int ANegEquNegNodeOffset;
    int ANegEquPriNodeOffset;
    int APriEquPosNodeOffset;
    int APriEquNegNodeOffset;
    int APriEquPriNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // Matrix equation local pointer variables
    double * fPosEquPosNodePtr;
    double * fPosEquPriNodePtr;
    double * fNegEquNegNodePtr;
    double * fNegEquPriNodePtr;
    double * fPriEquPosNodePtr;
    double * fPriEquNegNodePtr;
    double * fPriEquPriNodePtr;

    double * qPosEquPosNodePtr;
    double * qPosEquPriNodePtr;
    double * qNegEquNegNodePtr;
    double * qNegEquPriNodePtr;
    double * qPriEquPosNodePtr;
    double * qPriEquNegNodePtr;
    double * qPriEquPriNodePtr;
#endif


    // Flags
    bool TEMP_GIVEN;
    bool AREA_GIVEN;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/00
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

    double IS;   // saturation current (A)
    double RS;   // ohmic resistance (ohms)
    double COND; // corresponding conductance
    double N;    // emission coefficient
    double ISR;  // recombination saturation current (A)
    double NR;   // emission coefficient for ISR
    double IKF;  // high-injection knee current (A)
    double TT;   // transit time (sec)
    double CJO;  // zero-bias junction capacitance (F)
    double VJ;   // built-in junction potential (V)
    double M;    // grading coefficient
    double EG;   // activation  energy (eV).
    //    For Si, EG = 1.11
    //        Ge, EG = 0.67
    //        Sbd, EG = 0.69
    double XTI;  // isaturation-current temp. exp
    double TIKF; // IKF temperature coeff.
    double TBV1; // BV linear temperature coeff.
    double TBV2; // BV quadratic temperature coeff.
    double TRS1; // RS linear temperature coeff.
    double TRS2; // RS quadratic temperature coeff.
    double FC;   // coefficient for forward-bias depletion capacitance
    // formula
    double BV;   // reverse breakdown voltage
    double IBV;  // current at  breakdown voltage (A)
    double IRF;  // adjustment for linear portion of reverse current
    double NBV;  // reverse breakdown ideality factor
    double IBVL; // low-level current at  breakdown voltage (A)
    double NBVL; // low-level reverse breakdown ideality factor
    double F2;
    double F3;
    double TNOM; // parameter measurement temperature (C)
    double KF;   // flicker noise coefficient
    double AF;   // flicker noise exponent

    bool BVGiven;
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

    // new DAE stuff:
    // new DAE load functions, residual:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);

    // new DAE load functions, Jacobian:
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

    friend class Instance;
    friend class Model;
};

} // namespace Diode
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Diode::Instance N_DEV_DiodeInstance;
typedef Xyce::Device::Diode::Model N_DEV_DiodeModel;
typedef Xyce::Device::Diode::Master N_DEV_DiodeMaster;

#endif

