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
// Filename       : $RCSfile: N_DEV_Synapse2.h,v $
//
// Purpose        : Synapse2 classes
//
// Special Notes  :
//
// Creator        : Christy Warrender, SNL, Cognitive Modeling
//
// Creation Date  : 11/18/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Synapse2_h
#define Xyce_N_DEV_Synapse2_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <Sacado.hpp>

namespace Xyce {
namespace Device {
namespace Synapse2 {

// ---------- Forward Declarations -------
class Model;
class Instance;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
//
//	This  is  the  instance class  for Synapse2s.  It
//	contains "unique" Synapse2  information - ie stuff that
//	will be true of only one  Synapse2 in the circuit, such
//	as the nodes to which it is connected.  A Synapse2 is
//	connected to only two circuit nodes.
//
//	This class  does not directly contain information about
//	its node indices. It contains indices into the 3 parts
//	(A, dx, and  b) of the matrix  problem A*dx = b, and
//	also x.  A is the Jacobian  matrix, dx is the update to
//	the solution vector x, and b is the right hand side
//	function vector.  These indices are global, and
//	determined by topology during  the initialization stage
//	of execution.
//
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
    static vector< vector<int> > jacStamp;
    static ParametricData<Instance> &getParametricData();

    virtual const ParametricData<void> &getMyParametricData() const {
      return getParametricData();
    }

    Instance(InstanceBlock & IB,
             Model & Riter,
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

    bool processParams (string param = "");

    bool updateTemperature(const double & temp_tmp);

    bool updateIntermediateVars ();
    bool updatePrimaryState ();
    bool updateSecondaryState ();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    // load functions, residual:
    bool loadDAEQVector ();
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx ();
    bool loadDAEdFdx ();

    void setupPointers();

    Model &getModel() {
      return model_;
    }

  private:

    Model &       model_;         //< Owning model

    // user-specified paramters:


    //Vector local index for Positive Node
    int li_Prev;
    //Vector local index for Negative Node
    int li_Post;
    int li_rVar;

    // Offset variables corresponding to the above declared indices.
    int APostEquPostNodeOffset;
    int APostEquRNodeOffset;
    int AREquPostNodeOffset;
    int AREquRNodeOffset;

    // Pointers for Jacobian
    double *f_PostEquPostNodePtr;
    double *f_PostEquRNodePtr;
    double *f_REquPostNodePtr;
    double *f_REquRNodePtr;

    // vars used for load calculations
    double ipost;  // post Synapse2 current
    double didVpost;
    double didr;
    double rFval;
    double drFdVpre;
    double drFdr;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
//
//
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

    // Synapse2 parameters
    double gMax;
    double eRev;
    double alpha;
    double beta;
    double vP;
    double kP;
    double tMax;
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
      const string dn,
      const string cn,
      const string dmName,
      LinearDevice        linearDev,
      SolverState & ss1,
      DeviceOptions & do1)
      : Xyce::Device::DeviceTemplate<Model, Instance>(
        dn, cn, dmName, linearDev, ss1, do1)
    {}

    virtual bool updateState (double * solVec, double * staVec, double * stoVec);
    virtual bool updateSecondaryState (double * staDeriv, double * stoVec);

    friend class Instance;
    friend class Model;
};


// These are the base equations for the Synapse2 device.
// They're placed here to get compiled by the auto-diff tool Sacado
// These functions represent the equations that need to be solved
// for this device.

template <typename ScalarT>
static ScalarT PostCurrentEqu( const ScalarT Vpost, const ScalarT r, const ScalarT g, const ScalarT Erev)
{
  ScalarT result =  g * r * (Vpost - Erev);
  return result;
}

template <typename ScalarT>
static ScalarT rEquF( const ScalarT V, const ScalarT r, const ScalarT alpha, const ScalarT beta,
                      const ScalarT Tmax, const ScalarT Vthres)
{
  ScalarT result;
  if (V > Vthres)
  {
    result =  (alpha * Tmax * (1.0 - r) - beta * r);
  }
  else
  {
    result = - beta * r;
  }
  return result;
}

} // namespace Synapse2
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Synapse2::Instance N_DEV_SynapseInstance2;
typedef Xyce::Device::Synapse2::Model N_DEV_SynapseModel2;
typedef Xyce::Device::Synapse2::Master N_DEV_Synapse2Master;

#endif
