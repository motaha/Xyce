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
// Filename       : $RCSfile: N_DEV_Vcvs.h,v $
//
// Purpose        : Voltage  contolled current source classes.
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
// Revision Number: $Revision: 1.66.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Vcvs_h
#define Xyce_N_DEV_Vcvs_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Vcvs {

// ---------- Forward Declarations ----------
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

    Instance(InstanceBlock & IBref,
			   Model & Viter,
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

    bool updateIntermediateVars () { return true; };
    bool updatePrimaryState ();

    // load functions, residual:
    bool loadDAEQVector () {return true;}
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx () {return true;}
    bool loadDAEdFdx ();

    void setupPointers();

    void varTypes( vector<char> & varTypeVec );

  public:
    // iterator reference to the vcvs model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:
  static vector< vector<int> > jacStamp;

  Model &       model_;         //< Owning model

    double Gain;
    double LeadCurrent;

    InstanceBlock IB;

    // Matrix equation index variables:

    // local indices (offsets)
    int li_Pos;
    int li_Neg;
    int li_Bra;

    int li_ContPos;
    int li_ContNeg;

    // Offset variables for all of the above index variables.
    int ABraEquPosNodeOffset;
    int ABraEquNegNodeOffset;
    int ABraEquContPosNodeOffset;
    int ABraEquContNegNodeOffset;
    int APosEquBraVarOffset;
    int ANegEquBraVarOffset;

    // Ptr variables for all of the above index variables.
    double * f_BraEquPosNodePtr;
    double * f_BraEquNegNodePtr;
    double * f_BraEquContPosNodePtr;
    double * f_BraEquContNegNodePtr;
    double * f_PosEquBraVarPtr;
    double * f_NegEquBraVarPtr;
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

    Model (const ModelBlock & MB,
                           SolverState & ss1,
                       DeviceOptions & do1);
    ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual std::ostream &printOutInstances(std::ostream &os) const;

    virtual bool processParams(std::string param = "")
    {}

    virtual bool processInstanceParams(std::string param = "")
    {}


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

    // load functions:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

    friend class Instance;
    friend class Model;
};


} // namespace Vcvs
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Vcvs::Instance N_DEV_VcvsInstance;
typedef Xyce::Device::Vcvs::Model N_DEV_VcvsModel;
typedef Xyce::Device::Vcvs::Master N_DEV_VcvsMaster;

#endif
