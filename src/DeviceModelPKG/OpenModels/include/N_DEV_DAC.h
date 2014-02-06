//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DAC.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Lon Waters
//
// Creation Date  : 07/26/2002
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revsion$
//
// Revsion Date   : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DAC_h
#define Xyce_N_DEV_DAC_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>
#include <N_UTL_BreakPoint.h>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Device {
namespace DAC{

class Model;

//----------------------------------------------------------------------------
// Class          : Instance
// Purpose        : This class refers to a single instance of the DAC device.
//                  It contains indices into the matrix equation. See comments
//                  for the ResistorInstance class for more details.
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 07/26/2002
//----------------------------------------------------------------------------
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
             Model & DACiter,
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

    bool updateIntermediateVars ();
    bool updatePrimaryState ();
    bool updateSecondaryState ();

    bool updateTVVEC ( vector< pair<double, double> > const & newPairs );
    bool getInstanceBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes);

    DeviceState * getInternalState();
    bool setInternalState( const DeviceState & state );

    // iterator reference to the model which owns this instance.
    // Getters and setters
    Model &getModel() {
      return model_;
    }

  public:

    bool loadDAEQVector () {return true;};
    bool loadDAEFVector ();

    bool loadDAEdQdx () {return true;};
    bool loadDAEdFdx ();

    void varTypes( vector<char> & varTypeVec );

  private:

    bool updateVoltage(double);

  private:
    static vector< vector<int> > jacStamp;

    Model &       model_;         //< Owning model

    // user-specified parameters:
    string file;    // User specified file containing DAC data as time and
    //voltage pairs.

    vector< pair<double, double> > TVVEC; // vector of (time, voltage) pairs
    // read in from "file".

    // state variables:

    // other variables:
    int numTVpairs_;  // number of (time, voltage) pairs in TVVEC
    double v_pos;
    double v_neg;
    double i_bra;
    double vDrop;
    double voltage_;
    int loc_;

    // Indices into the state vector:



    //local indices (offsets)
    int li_Pos;
    int li_Neg;
    int li_Bra;

    //Locally indexed offsets for jacobian
    int ABraEquPosNodeOffset; // Offset, pos. node voltage contribution,
    // branch current equ.

    int ABraEquNegNodeOffset; // Offset, neg. node voltage contribution,
    // branch current equ.

    int APosEquBraVarOffset;  // Offset, branch current variable
    // contribution, KCL equation of the pos node

    int ANegEquBraVarOffset;  // Offset, branch current variable
    // contribution, KCL equation of the neg node
};

//----------------------------------------------------------------------------
// Function       : Model
// Purpose        :
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 07/26/2002
//----------------------------------------------------------------------------
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
    bool processParams (string param = "");
    bool processInstanceParams (string param = "");
    virtual std::ostream &printOutInstances(std::ostream &os) const;

  private:

    // Model Parameters
    double riseTime;
    double fallTime;
    double R;
    double L;
    double C;
    bool includeTransitionBP_;

    // Data Members for Associations

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

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 02/25/2009
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

    // load functions:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

    bool getDACDeviceNames ( std::vector<std::string> & dacNames);

    friend class Instance;
    friend class Model;
};

} // namespace DAC
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DAC::Instance N_DEV_DACInstance;
typedef Xyce::Device::DAC::Model N_DEV_DACModel;
typedef Xyce::Device::DAC::Master N_DEV_DACMaster;

#endif
