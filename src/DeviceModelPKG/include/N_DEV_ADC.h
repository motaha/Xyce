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
// Filename       : $RCSfile: N_DEV_ADC.h,v $
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

#ifndef Xyce_N_DEV_ADC_h
#define Xyce_N_DEV_ADC_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>

#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>
#include <N_UTL_BreakPoint.h>


namespace Xyce {
namespace Device {
namespace ADC {

class Model;

//----------------------------------------------------------------------------
// Class          : Instance
// Purpose        : This class refers to a single instance of the ADC device.
//                  It contains indices into the matrix equation. See comments
//                  for the ResistorInstance class for more details.
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 07/26/2002
//----------------------------------------------------------------------------
class Instance : public N_DEV_DeviceInstance
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
           Model & ADCiter,
           MatrixLoadData & mlData1,
           SolverState &ss1,
           ExternData  &ed1,
           DeviceOptions & do1);

  Instance(const Instance &right);

  ~Instance();

  // Additional Public Declarations
  void registerLIDs( const vector<int> & intLIDVecRef,
                     const vector<int> & extLIDVecRef );
  void registerStateLIDs( const vector<int> & staLIDVecRef );

  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

  bool processParams (string param = "");

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  void getTVVEC(vector< pair<double, double> > & TVVVEC_Out);
  void trimTVVEC(double earliestTime);
  bool getInstanceBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes);
  void acceptStep();

  bool getInstanceParamsMap(map<string,double>& paramsMap);
  int getNumberQuantLevels();
  bool setBitVectorWidth(int width);

  // iterator reference to the model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:
  Model &       model_;         //< Owning model

public:

  bool loadDAEQVector () {return true;};
  bool loadDAEFVector ();

  bool loadDAEdQdx () {return true;};
  bool loadDAEdFdx ();

private:
  static vector< vector<int> > jacStamp;

  // parameters:
  double R;

  // derived parameters:
  double G;  // conductance (1.0/ohms)
  double i0; // current (amps)
  double v_pos, v_neg;

  // state variables:

  // other variables:
  vector< pair<double, double> > TVVEC;
  int outputBitVectorWidth_; // number of bits on digital side
  int nQuantLevels_;         // 2^(outputBitVectorWidth_)
  int lastOutputLevel_;      // save last value, so we know when we're going
  // to cause a change

  // Indices into the state vector:

  //local indices (offsets)
  int li_Pos;
  int li_Neg;

  //Locally indexed offsets for jacobian
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int ANegEquPosNodeOffset;
  int ANegEquNegNodeOffset;
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
    double lowerVoltageLimit_;
    double upperVoltageLimit_;
    double settlingTime_;

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
    friend class Instance;
    friend class Model;

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

    bool getBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes);
    bool getADCMap(map<string,map<string,double> >& ADCMap);
    bool getTimeVoltagePairs( map<string,vector<pair<double,double> > >& tvvmap);
};

//-----------------------------------------------------------------------------
// Function      : Instance::getNumberQuantLevels
// Purpose       : return number of possible quantization levels of an ADC
// Special Notes : Assumes that width has already been set!
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/07/2004
//-----------------------------------------------------------------------------

inline int Instance::getNumberQuantLevels()
{
  return nQuantLevels_;
}

} // namespace ADC
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ADC::Instance N_DEV_ADCInstance;
typedef Xyce::Device::ADC::Model N_DEV_ADCModel;
typedef Xyce::Device::ADC::Master N_DEV_ADCMaster;

#endif
