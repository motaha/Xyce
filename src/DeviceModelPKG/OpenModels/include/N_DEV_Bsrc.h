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
// Filename       : $RCSfile: N_DEV_Bsrc.h,v $
//
// Purpose        : General expression dependent source
//
// Special Notes  :
//
// Creator        : Robert J Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/05/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.79.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Bsrc_h
#define Xyce_N_DEV_Bsrc_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_Param.h>
#include <N_UTL_BreakPoint.h>

// ---------- Forward Declarations ----------
class N_UTL_Expression;

namespace Xyce {
namespace Device {
namespace Bsrc {

class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
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
             Model & BMiter,
             MatrixLoadData & mlData1,
             SolverState &ss1,
             ExternData  &ed1,
             DeviceOptions & do1);

    ~Instance();

  private:
    Instance(const Instance &);
    Instance &operator=(const Instance &);

  public:
    void registerLIDs(const vector<int> & intLIDVecRef,
                      const vector<int> & extLIDVecRef );
    void registerStateLIDs(const vector<int> & staLIDVecRef);
    void registerStoreLIDs(const vector<int> & stoLIDVecRef );

    map<int,string> & getIntNameMap ();
    map<int,string> & getStoreNameMap ();

    const vector<std::string> & getDepSolnVars();


    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");

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

    void varTypes( vector<char> & varTypeVec );

  public:
    // iterator reference to the bsrc model which owns this instance.
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
    list<std::string>   evnList;

    vector<double> expVarDerivs;
    vector<double> myVarVals;
    vector<double> ddtVals;
    double         expVal;

    InstanceBlock IB;

    // flag for voltage src, needs current variable
    bool isVSRC;

    // Value of voltage or current expression
    double V;
    double I;

    // scale factor
    double scale;
    int nlstep;

    // indices into state vector:
    vector<int>    li_ddt;

    // solution vector indices:
    // rhs vector indices:
    int li_Pos;
    int li_Neg;
    int li_Bra;
    int li_store_branch;  // branch current stored in store vector
    // if it is not part of the solution vector

    // Local offset variables for all of the above index variables.
    int ABraEquPosNodeOffset;
    int ABraEquNegNodeOffset;
    int APosEquBraVarOffset;
    int ANegEquBraVarOffset;

    vector<int> APosEquExpVarOffsets;
    vector<int> ANegEquExpVarOffsets;
    vector<int> ABraEquExpVarOffsets;

    // Local offset variables for all of the above index variables.
    double * fBraEquPosNodePtr;
    double * fBraEquNegNodePtr;
    double * fPosEquBraVarPtr;
    double * fNegEquBraVarPtr;

    vector<double *> fPosEquExpVarPtrs;
    vector<double *> fNegEquExpVarPtrs;
    vector<double *> fBraEquExpVarPtrs;

    vector< vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
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
    ~Model ();

  private:
    Model();
    Model(const Model &);
    Model &operator=(const Model &);

  public:
    virtual std::ostream &printOutInstances(std::ostream &os) const;
    virtual bool processParams(std::string param = "");
    virtual bool processInstanceParams(std::string param = "");

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


    // This is the dc and transient analysis value of the
    // source.
    double DC_TRAN;

    // This is the AC magnitude
    double ACMAG;

    // This is the AC phase.
    double ACPHASE;

    // This parameter is part of the specification that the
    // source has distortion inputs at a frequency of this
    // magnitude.  It is triggered by the DISTOF1 keyword in
    // the netlist.
    double F1MAG;

    // This parameter is part of the specification that the
    // source has distortion inputs at a frequency of this
    // magnitude.  It is triggered by the DISTOF2 keyword in
    // the netlist.
    double F2MAG;

    // This parameter is associated with the specification that
    // the source has a distortion input.  This is the phase
    // associated with DISTOF1.
    double F1PHASE;

    // This parameter is associated with the specification that
    // the source has a distortion input.  This is the phase
    // associated with DISTOF2.
    double F2PHASE;
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

} // namespace Bsrc
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Bsrc::Instance N_DEV_BsrcInstance;
typedef Xyce::Device::Bsrc::Model N_DEV_BsrcModel;
typedef Xyce::Device::Bsrc::Master N_DEV_BsrcMaster;

#endif

