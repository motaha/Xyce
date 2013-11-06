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
// Filename       : $RCSfile: N_DEV_Vsrc.h,v $
//
// Purpose        : Independent voltage source classes.
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
// Revision Number: $Revision: 1.90.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Vsrc_h
#define Xyce_N_DEV_Vsrc_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Vsrc {

// ---------- Forward Declarations ----------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
class Instance : public SourceInstance
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
			   Model & Viter,
                        MatrixLoadData & mlData1,
                        SolverState &ss1,
                        ExternData  &ed1,
                        DeviceOptions & do1);

    Instance(const Instance & right);
    ~Instance();

    void registerLIDs( const vector<int> & intLIDVecRef,
	               const vector<int> & extLIDVecRef );
    void registerStateLIDs( const vector<int> & staLIDVecRef );
    map<int,string> & getIntNameMap ();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");

    bool updateIntermediateVars ();
    bool updatePrimaryState ();

    // load functions, residual:
    bool loadBVectorsforAC (double * bVecReal, double * bVecImag);
    bool loadDAEQVector () { return true; }
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx () { return true; }
    bool loadDAEdFdx ();

    void setupPointers();

    double getMaxTimeStepSize ();

    void varTypes( vector<char> & varTypeVec );

    void getLIDs(int & lpos, int & lneg,int & lbra)
      {lpos = li_Pos; lneg = li_Neg; lbra = li_Bra;}

  public:
    // iterator reference to the vsrc model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

    // state variables:
    double srcCurrent;
    double srcVoltage;
    double srcDrop;
    double srcBC;

    // scale factor
    double scale;
    int nlstep;

    // Parameters
    double DCV0;
    double par0;
    double par1;
    double par2;
    double par3;
    double par4;
    double par5;
    double par6;
    double par7;
    double REPEATTIME;
    double T;
    double V;
    double ACMAG;
    double ACPHASE;

    int NUM;
    int REPEAT;
    int TRANSIENTSOURCETYPE;
    bool TRANSIENTSOURCETYPEgiven;
    int ACSOURCETYPE;
    bool ACSOURCETYPEgiven;
    int DCSOURCETYPE;
    bool DCSOURCETYPEgiven;
    bool gotParams;

    // load variables
    double source, v_pos, v_neg, i_bra;

    // indices into state vector:
    int istate_I;  // index for i0;

    // Matrix equation index variables:

    //local indices (offsets)
    int li_Pos;
    int li_Neg;
    int li_Bra;

    // Jacobian matrix indices:
    //Locally indexed offsets for jacobian
    int ABraEquPosNodeOffset; // Offset, pos. node voltage contribution,
                              // branch current equ.

    int ABraEquNegNodeOffset; // Offset, neg. node voltage contribution,
                              // branch current equ.

    int APosEquBraVarOffset;  // Offset, branch current variable
                              // contribution, KCL equation of the pos node

    int ANegEquBraVarOffset;  // Offset, branch current variable
                              // contribution, KCL equation of the neg node

    //  The following jacobian offsets are only neccessary
    // for 2-level newton.
    int APosEquPosNodeOffset;  // Offset, positive node variable
                              // contribution, positive node KCL.

    int ANegEquNegNodeOffset;  // Offset, negative node variable
                              // contribution, negative node KCL.

    int ABraEquBraVarOffset;  // Offset, branch current variable
                              // contribution, branch current equation.

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // Jacobian matrix pointers:
    double * fBraEquPosNodePtr;
    double * fBraEquNegNodePtr;
    double * fPosEquBraVarPtr;
    double * fNegEquBraVarPtr;

    //  The following jacobian pointers are only neccessary for 2-level newton.
    double * fPosEquPosNodePtr;
    double * fNegEquNegNodePtr;
    double * fBraEquBraVarPtr;
#endif

    static vector< vector<int> > jacStamp;
    static vector< vector<int> > jacStampPDE;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
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

    // This is the dc and transient analysis value of the source.
    double DC_TRAN;
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

//    bool loadBVectorsforAC(double * bVecReal, double * bVecImag );

    friend class Instance;
    friend class Model;
};

} // namespace Vsrc
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Vsrc::Instance N_DEV_VsrcInstance;
typedef Xyce::Device::Vsrc::Model N_DEV_VsrcModel;
typedef Xyce::Device::Vsrc::Master N_DEV_VsrcMaster;

#endif
