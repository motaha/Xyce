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
// Filename       : $RCSfile: N_DEV_Region.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 07/19/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Region_h
#define Xyce_N_DEV_Region_h

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>
#include <N_DEV_fwd.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_ReactionNetwork.h>
#include <N_DEV_Specie.h>

// ---------- Forward Declarations -------
class N_LAS_Vector;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : Region
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 04/14/11
//-----------------------------------------------------------------------------
class Region
{
  // functions:
public:
  Region ( RegionData & rd,
           DeviceOptions & devOp,
           SolverState & solst, bool sourceOn=true);

  // constructor that allows a non-default reaction network:
  Region ( RegionData & rd,
           DeviceOptions & devOp,
           SolverState & solst,
           N_DEV_ReactionNetwork & reactionNet);

  virtual ~Region ();

  bool outputTecplot ();

///
  void initializeReactionNetwork(ScalingVars & sv);
  void setInitialCondition (const std::string & name, const double val);
  void setRateConstants(double T);
  void setupScalingVars (ScalingVars & sv);
  void scaleVariables ();
  void unscaleVariables ();
  void scaleRateConstants ();
  void unscaleRateConstants ();
  void addSource(string speciesName, N_UTL_Expression *expr);
  void addMasterSource(string speciesName);
  inline void setMasterSourceValue(double msv) {theReactions.setMasterSourceValue(msv);};
  inline void setSimTime(double time) {  theReactions.setSimTime(time);};

  bool reactantExist(string reactantname)
  {
    return theReactions.reactantExist(reactantname);
  };

  bool constantExist(string constantname)
  {
    return theReactions.constantExist(constantname);
  };

  bool getDoNothingFlag ();
  inline string getName () { return name; }

  double getBreakTime();

  void setupJacStamp ( vector< vector<int> > & jacStamp, vector<int> & colDep, int & firstReactant, int & lastIndex );

  void registerLIDs( const vector<int> & intLIDVec, const vector<int> & extLIDVec, int & intIndex);

  void augmentNameMap ( map<int,string> & intNameMap, DeviceInstance & di);

  void registerStateLIDs (const vector<int> & staLIDVecRef, int & i);

  void registerJacLIDs ( const vector< vector<int> > & jacLIDVec,
                         const vector<int> &map,
                         const vector< vector<int> > &map2 );

  void setupPointers (N_LAS_Matrix & dfdx, N_LAS_Matrix & dqdx);

  void updateIntermediateVars ( double * solVector, double * oldSolVector, double time);

  bool loadDAEQVector (double * qVec);
  bool loadDAEFVector (double * fVec);
  bool loadDAEdFdxdV (double * dfdxdv,double vdiff);

  bool loadDAEdQdx (N_LAS_Matrix & dqdx);
  bool loadDAEdFdx (N_LAS_Matrix & dfdx);

  bool loadDeviceMask (N_LAS_Vector & mask);

  bool updateSecondaryState (double * staDeriv);

  // These three simple accessors are here so we can avoid having *any*
  // public data.
  double getStateConcentration(int i) {
    return tempConcentrations[i];
  }

  bool haveAnyReactions();

  int getStateConcentrationLID(int i) {
    return li_state_Concentrations[i];
  }

  void setConstantConcentration(const std::string &constName, double value);

  // This is here to allow us to copy a full ReactionNetwork object from
  // one region to another without having to re-parse every time.
  // We return a const reference so nobody can tinker with our internal
  // data.
  const N_DEV_ReactionNetwork & getReactionNetwork() {
    return theReactions;
  }

  // These query methods here so we can communicate how many equations we're
  // adding to the base device
  int getNumIntVars();

  int getNumSpecies() {
    return theReactions.getNumSpecies();
  }

  int getNumConstants() {
    return theReactions.getNumConstants();
  }

  int getSpeciesLID (const std::string &name);

  double getDiffusionCoefficient (const std::string & name, double temp);
  double getDiffusionCoefficient (int specie, double temp);
  double getConcentrationScaling();
  double getLengthScaling();
  int getSpeciesNum(const std::string & name);
  const std::string & getSpeciesName(int i);
  const std::string & getConstantsName(int i);
  const double getSpeciesVal(int i);
  const double getConstantsVal(int i);

///

private:
  void createDefaultReactionNetwork(const std::string & reactionSpecFile);

public:
  RegionData & regData;

protected:
  // data:
  string name;
  string outputName;

  bool explicitCarrierFlag;
  bool useScaledVariablesFlag;
  bool variablesScaledFlag;
  bool rateConstantsScaledFlag;
  int callsOTEC;

  // reactions
  N_DEV_ReactionNetwork theReactions;
  // vector of constant concentrations (species held fixed)
  vector<double> theConstantConcentrations;
  // working storage for communicating between updateIntermediateVars
  // and updatePrimaryState
  vector<double> tempConcentrations;
  // actual time derivatives of concentrations
  vector<double> tempConcentrationDerivs;
  // initial conditions
  vector<double> initialConcentrations;
  vector<double> ddt;
  vector< vector<double> > tempJac;
  vector< vector<double> > tempAuxJac;

  int baseReactionIndex;

  vector< vector<double *> > dfdxConcEquConcVarPtrs;
  vector< vector<double *> > dqdxConcEquConcVarPtrs;

  vector< vector<double *> > dfdxConcEquAuxVarPtrs;
  vector< vector<double *> > dqdxConcEquAuxVarPtrs;

  vector< vector<int> > AConcentrationEquConcentrationNodeOffsets;
  vector< vector<int> > AConcentrationEquAuxNodeOffsets;

  // reaction species indices:
  vector<int> li_Concentrations;

  // reaction state vars... these are redundant storage, because we also
  // need concentration derivatives
  vector<int> li_state_Concentrations;

  // Rxn set scaling variables.
  double x0;  // distance scaling (cm)
  double a0;  // area scaling (cm^2)
  double C0;  // concentration scaling (cm^-3);
  double D0;  // diffusion coefficient scaling (cm^2/s)
  double u0;  // mobility coefficient scaling (cm^2/V/s)
  double R0;  // recombination rate scaling (cm^-3/s)
  double rR0; // reciprocal of R0
  double t0;  // time scaling (s)
  double k0;  // rate constant scaling (cm^3/s)
  double rt0; // reciprocal
  double rk0; // reciprocal

  bool outputBefore1;
  bool outputBefore2;

  DeviceOptions & devOptions;
  SolverState & solState;
};

//-----------------------------------------------------------------------------
// Function      : Region::setConstantConcentration
// Purpose       : Set value of concentration for a constant species.
// Special Notes : The species had better be a constant already defined.
//                 No error checking is done.  So sue me.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/16/06
//-----------------------------------------------------------------------------

inline void Region::setConstantConcentration(const std::string &constName,
                                             double value)
{
  theConstantConcentrations[theReactions.getConstantNum(constName)]=
    (value)*((variablesScaledFlag)?(1/C0):1.0);
}

//-----------------------------------------------------------------------------
// Function      : Region::getSpeciesLID
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/30/06
//-----------------------------------------------------------------------------
inline int Region::getSpeciesLID (const std::string &namearg)
{
  return li_Concentrations [theReactions.getSpeciesNum(namearg)];
}

//-----------------------------------------------------------------------------
// Function      : Region::getDiffusionCoefficient
// Purpose       :
// Special Notes : If scaling is on, this function returns the scaled value.
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/30/06
//-----------------------------------------------------------------------------
inline double Region::getDiffusionCoefficient
(const std::string & namearg, const double temp)
{
  return
    ((variablesScaledFlag)?(D0):1.0)*
    theReactions.getDiffusionCoefficient(namearg,temp);
}

//-----------------------------------------------------------------------------
// Function      : Region::getDiffusionCoefficient
// Purpose       :
// Special Notes : If scaling is on, this function returns the scaled value.
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/31/06
//-----------------------------------------------------------------------------
inline double Region::getDiffusionCoefficient
(int specie, double temp)
{
  return
    ((variablesScaledFlag)?(D0):1.0)*
    theReactions.getDiffusionCoefficient(specie,temp);
}

//-----------------------------------------------------------------------------
// Function      : Region::getConcentrationScaling
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/30/06
//-----------------------------------------------------------------------------
inline double Region::getConcentrationScaling()
{
  return C0;
}

//-----------------------------------------------------------------------------
// Function      : Region::getLengthScaling
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/30/06
//-----------------------------------------------------------------------------
inline double Region::getLengthScaling()
{
  return x0;
}

//-----------------------------------------------------------------------------
// Function      : Region::getSpeciesNum
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/30/06
//-----------------------------------------------------------------------------
inline int Region::getSpeciesNum(const std::string & namearg)
{
  return theReactions.getSpeciesNum(namearg);
}

//-----------------------------------------------------------------------------
// Function      : Region::getSpeciesameN
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/30/06
//-----------------------------------------------------------------------------
inline const std::string & Region::getSpeciesName(int i)
{
  return theReactions.getSpeciesName(i);
}

//-----------------------------------------------------------------------------
// Function      : Region::getConstantsName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/30/06
//-----------------------------------------------------------------------------
inline const std::string & Region::getConstantsName(int i)
{
  return theReactions.getConstantsName(i);
}

//-----------------------------------------------------------------------------
// Function      : Region::getSpeciesVal
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 04/15/07
//-----------------------------------------------------------------------------
inline const double Region::getSpeciesVal(int i)
{
  return (tempConcentrations[i]*((variablesScaledFlag)?(C0):(1.0)));
}

//-----------------------------------------------------------------------------
// Function      : Region::getConstantsVal
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 04/15/07
//-----------------------------------------------------------------------------
inline const double Region::getConstantsVal(int i)
{
  return (theConstantConcentrations[i] *((variablesScaledFlag)?(C0):(1.0)));
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Region N_DEV_Region;

#endif

