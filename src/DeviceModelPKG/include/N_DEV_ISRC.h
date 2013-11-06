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
// Filename       : $RCSfile: N_DEV_ISRC.h,v $
//
// Purpose        : Independent current source
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
// Revision Number: $Revision: 1.79.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ISRC_h
#define Xyce_N_DEV_ISRC_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>

#include <N_UTL_BreakPoint.h>

// ---------- Forward Declarations ----------
namespace Xyce {
namespace Device {
namespace ISRC {

class Model;
class Instance;

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

    Instance(const Instance &right);
    Instance(InstanceBlock & IB,
			   Model & Iiter,
                        MatrixLoadData & mlData1,
                        SolverState &ss1,
                        ExternData  &ed1,
                        DeviceOptions & do1);

    ~Instance();

    void registerLIDs( const vector<int> & intLIDVecRef,
	               const vector<int> & extLIDVecRef);
    void registerStateLIDs( const vector<int> & staLIDVecRef);

    void registerStoreLIDs(const vector<int> & stoLIDVecRef );

    map<int,string> & getStoreNameMap ();

    const vector< vector<int> > & jacobianStamp() const;

    bool processParams (string param = "");

    bool updateIntermediateVars () { return true;};
    bool updatePrimaryState ();

    bool loadTrivialMatrixStamp ();
    bool loadTrivialDAE_FMatrixStamp ();

    // load functions, residual:
    bool loadDAEQVector () { return true; }
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx () { return true; }
    bool loadDAEdFdx () { return true; }

    bool loadBVectorsforAC (double * bVecReal, double * bVecImag);

  protected:
  private:

  public:
    // iterator reference to the ISRC model that owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:
  static vector< vector<int> > jacStamp;

  Model &       model_;         //< Owning model

    // index variables:
    int li_Pos;
    int li_Neg;

    // Store variables
    int li_store_dev_i;

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
    int NUM;
    int REPEAT;
    int TRANSIENTSOURCETYPE;
    bool TRANSIENTSOURCETYPEgiven;
    int ACSOURCETYPE;
    bool ACSOURCETYPEgiven;
    int DCSOURCETYPE;
    bool DCSOURCETYPEgiven;
    bool gotParams;


    double ACMAG;
    double ACPHASE;
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

  friend class Instance;
  friend class Master;

  public:
  static ParametricData<Model> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

    Model  (const ModelBlock & MB,
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
};

//-----------------------------------------------------------------------------
// Function      : Instance::loadTrivialMatrixStamp ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/24/03
//-----------------------------------------------------------------------------
inline bool Instance::loadTrivialMatrixStamp ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadTrivialDAE_FMatrixStamp ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/24/03
//-----------------------------------------------------------------------------
inline bool Instance::loadTrivialDAE_FMatrixStamp ()
{
  return true;
}
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


} // namespace ISRC
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ISRC::Instance N_DEV_ISRCInstance;
typedef Xyce::Device::ISRC::Model N_DEV_ISRCModel;
typedef Xyce::Device::ISRC::Master N_DEV_ISRCMaster;

#endif

