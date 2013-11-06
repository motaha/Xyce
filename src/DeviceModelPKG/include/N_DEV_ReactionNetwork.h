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
// Filename       : $RCSfile: N_DEV_ReactionNetwork.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 03/20/2006
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#ifndef N_DEV_ReactionNetwork_H
#define N_DEV_ReactionNetwork_H
#include <N_UTL_Misc.h>
#include <iosfwd>
#include <vector>
#include <map>
#include <string>
using namespace std;

#include <N_DEV_Reaction.h>
#include <N_DEV_Specie.h>

#include <N_ERH_ErrorMgr.h>

// Note: sstream and strstream are only needed here because of all the
// inlined functions.
#include <sstream>

class N_UTL_Expression;

//-----------------------------------------------------------------------------
// Class         : N_DEV_ReactionNetwork
// Purpose       :
// Special Notes :
// Creator       : Tom Russo, SNL
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
class N_DEV_ReactionNetwork
{
public:

  N_DEV_ReactionNetwork(string name="NoName");
  N_DEV_ReactionNetwork(const N_DEV_ReactionNetwork & right);
  virtual ~N_DEV_ReactionNetwork();

  void setReactionNetworkFromFile(const std::string &fileName);

  void addReaction(string name);
  void addReaction(string name, N_DEV_Reaction &reaction);
  void addReactant(string name, string reactant ,double stoich);
  void addProduct(string name, string reactant ,double stoich);

  // set rate constant calculator for each type of reaction
  void setSimpleCalc(string name, double k);
  void setCaptureCalc(string name, double sigma, double v);
  void setEmissionCalc(string name, double sigma, double v, double N,
                       double E);
  void setElectronCaptureCalc(string name, double sigma);
  void setElectronEmissionCalc(string name, double sigma, double E);
  void setHoleCaptureCalc(string name, double sigma);
  void setHoleEmissionCalc(string name, double sigma, double E);
  void setComplexCalc(string name);
  void setDecomplexCalc(string name, double bindingEnergy,
                        double gammaAB, double gammaA, double gammaB,
                        double concSi);
  //

  void setRateConstant(string name, double k);
  void scaleRateConstant(string name, double kscale);

  void setScaleParams(double c, double t, double x);

  // Use rate constant calculators
  void scaleRateConstantFromCalculator(string name);
  void unscaleRateConstantFromCalculator(string name);
  void setRateConstantFromCalculator(string name,double T);
  void setRateConstantsFromCalc(double T);
  void scaleRateConstantsFromCalc();
  void unscaleRateConstantsFromCalc();

  void setSpecies(vector<N_DEV::Specie> &theSpeciesVect);
  void addSpecie(const N_DEV::Specie &aSpecie);
  void setConstants(vector<N_DEV::Specie> &theConstantsVect);
  void addConstant(const N_DEV::Specie &aConstant);

  int getReactionNum(string name);

  void addSourceTerm(string speciesName,string expressionStr);
  void addSourceTerm(string speciesName,N_UTL_Expression *expression);
  void addMasterSourceTerm(string speciesName);

  void addInitialCondition(string speciesName,double value);
  pair<string,double> getInitialCondition(int i);
  int getNumInitialConditions();

  void setSimTime(double time);
  double getBreakpointTime();
  inline void setSourceScaleFac(double scf) {sourceScaleFac=scf;};
  inline void setMasterSourceValue(double msv) {masterSourceValue=msv;};

  void getDdt(vector<double> &concs,vector<double> &constants,
              vector<double> &ddt);
  void getJac(vector<double> &concs, vector<double> &constants,
              vector<vector<double> >&jac);
  void getDFdConst(const std::string &constantName,
                   vector<double> &concs, vector<double> &constants,
                   vector<double> &dFdConst);

  double getRate(vector<double> &concs,vector<double> &constants,
                 vector<int> &captureVect, vector<int> &emissionVect);
  void getDRateDC(vector<double> &concs,vector<double> &constants,
                    vector<int> &captureVect, vector<int> &emissionVect,
                    vector<double>&dratedc);

  void getDRateDConst(vector<double> &concs,vector<double> &constants,
                      vector<int> &captureVect, vector<int> &emissionVect,
                      vector<double>&dratedc);

  double getCaptureLifetime(vector<double> &concs,vector<double> &constants,
                 vector<int> &captureVect,double &concentration);
  void getCaptureLifetimes(vector<double> &concs,vector<double> &constants,
                           vector<int> &captureVect,double &concentration,
                           vector<double> &lifetimes);
  double getELifetime(vector<double>&concs,vector<double>&constants);
  double getHLifetime(vector<double>&concs,vector<double>&constants);
  void getELifetimes(vector<double>&concs,vector<double>&constants,
                     vector<double> &lifetimes);
  void getHLifetimes(vector<double>&concs,vector<double>&constants,
                     vector<double> &lifetimes);

  double getERate(vector<double> &concs,vector<double> &constants);
  double getHRate(vector<double> &concs,vector<double> &constants);

  void getDERateDC(vector<double> &concs,vector<double> &constants,
                   vector<double>&dratedc);
  void getDHRateDC(vector<double> &concs,vector<double> &constants,
                   vector<double>&dratedc);

  void getDERateDConst(vector<double> &concs,vector<double> &constants,
                       vector<double>&dratedConst);
  void getDHRateDConst(vector<double> &concs,vector<double> &constants,
                       vector<double>&dratedConst);

  int getNumSpecies();
  int getNumConstants();
  const std::string & getSpeciesName(int i);
  const std::string & getConstantsName(int i);
  int getSpeciesNum(string name);
  int getConstantNum(string name);
  int getReactantNum(string name);
  bool reactantExist(string name);
  bool constantExist(string name);

  void setName(string name);
  void clear();
  void output(ostream & os) const;

  double getDiffusionCoefficient (const std::string & name, const double temp);
  double getDiffusionCoefficient (int specie, const double temp);

  int getChargeState (const std::string & name);
  int getChargeState (int specie);

  void setApplySources(bool flag);

private:
  N_DEV_Reaction &getReaction(string name);
  N_DEV_Reaction &getReaction(int i);

  map <string,int> speciesMap;
  vector<N_DEV::Specie> species;
  map <string,int> constantsMap;
  vector<N_DEV::Specie> constants;
  vector<pair<string,double> > initialConditions;
  vector<N_DEV_Reaction> theReactions;
  map <string,int> reactionNamesMap;
  std::vector<std::string> reactionNames;
  string myName;
  vector<int> electronCaptureReactions;
  vector<int> holeCaptureReactions;
  vector<int> electronEmissionReactions;
  vector<int> holeEmissionReactions;
  vector< pair<int,N_UTL_Expression *> > theSourceTerms;
  vector<int> masterSourceSpecies;
  double masterSourceValue;
  double sourceScaleFac;
  double C0;
  double t0;
  double x0;

  bool applySources;
};

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::clear
// Purpose       : make the reaction network an empty one
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::clear()
{
  speciesMap.clear();
  species.clear();
  constantsMap.clear();
  constants.clear();
  theReactions.clear();
  reactionNamesMap.clear();
  setName("NoName");
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setName
// Purpose       : Accessor function to set network name
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setName(string name)
{
  myName=name;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getNumSpecies
// Purpose       : Accessor function to returning number of species recorded
// Special Notes : Only solution species returned, does not include constants
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int N_DEV_ReactionNetwork::getNumSpecies()
{
  return (species.size());
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getNumConstants
// Purpose       : Accessor function to returning number of constant species
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
inline int N_DEV_ReactionNetwork::getNumConstants()
{
  return (constants.size());
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getSpeciesName
// Purpose       : Accessor function to returning name of indicated species
// Special Notes : Only returns solution species names, not constants
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline const std::string & N_DEV_ReactionNetwork::getSpeciesName(int i)
{
  return species[i].getName();
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getConstantsName
// Purpose       : Accessor function to returning name of indicated constant
//                 species
// Special Notes : Only returns solution constants names, not variable species
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline const std::string & N_DEV_ReactionNetwork::getConstantsName(int i)
{
  return constants[i].getName();
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getSpeciesNum
// Purpose       : Accessor function to return number of named variable specie
// Special Notes :  returns -1 if the specified specie is not a variable
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int N_DEV_ReactionNetwork::getSpeciesNum(string name)
{
  map<string,int>::iterator n_i;
  n_i = speciesMap.find(name);

  if (n_i == speciesMap.end())
  {
    return -1;
  }
  else
  {
    return n_i->second;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::reactantExist
// Purpose       : Returns false if the named specie doesn't exist, or if
//                 it is a constant.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 4/24/07
//-----------------------------------------------------------------------------
inline bool N_DEV_ReactionNetwork::reactantExist(string name)
{
  bool retFlag(true);
  int i=getSpeciesNum(name);
  if (i == -1)
  {
    retFlag = false;  // this reactant does not exist.
  }

  return retFlag;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::constantExist
// Purpose       : Returns false if the named specie doesn't exist, or if
//                 it is a constant.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 4/24/07
//-----------------------------------------------------------------------------
inline bool N_DEV_ReactionNetwork::constantExist(string name)
{
  bool retFlag(true);
  int i=getConstantNum(name);
  if (i == -1)
  {
    retFlag = false;  // this constant does not exist.
  }

  return retFlag;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getReactantNum
// Purpose       : Accessor function to return number of named specie
// Special Notes :  Returns negative numbers for constants, positive for
//                  variables
//                 The species in question better exist, coz there's no way
//                 to return an error condition other than bombing.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int N_DEV_ReactionNetwork::getReactantNum(string name)
{
  int i=getSpeciesNum(name);
  if (i == -1)
  {
    i=getConstantNum(name);
    if (i == -1)
    {
      string msg="N_DEV_ReactionNetwork::getReactantNum: invalid species name specified: ";
      msg += name;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }
    i = -(i+1);
  }
  return i;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getConstantNum
// Purpose       : Accessor function to return number of named constant specie
// Special Notes :  returns -1 if the specified specie is not a constant
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int N_DEV_ReactionNetwork::getConstantNum(string name)
{
  map<string,int>::iterator n_i;
  n_i = constantsMap.find(name);

  if (n_i == constantsMap.end())
  {
    return -1;
  }
  else
  {
    return n_i->second;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::addInitialCondition
// Purpose       : add an initial condition to the list
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::addInitialCondition(string speciesName,
                                                     double value)
{
  initialConditions.push_back(pair<string,double>(speciesName,value));
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getNumInitialConditions
// Purpose       : get an initial condition from the list
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline int
N_DEV_ReactionNetwork::getNumInitialConditions()
{
  return(initialConditions.size());
}
//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getInitialCondition
// Purpose       : get an initial condition from the list
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline pair<string,double>
N_DEV_ReactionNetwork::getInitialCondition(int i)
{
  return(initialConditions[i]);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getReactionNum
// Purpose       : Accessor function to return number of named reaction
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int N_DEV_ReactionNetwork::getReactionNum(string name)
{
  map<string,int>::iterator n_i;
  n_i = reactionNamesMap.find(name);

  if (n_i == reactionNamesMap.end())
  {
    return -1;
  }
  else
  {
    return n_i->second;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getReaction
// Purpose       : Accessor function to returning reference to indexed reaction
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline N_DEV_Reaction & N_DEV_ReactionNetwork::getReaction(int i)
{
  return theReactions[i];
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getReaction
// Purpose       : Accessor function to returning reference to named reaction
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline  N_DEV_Reaction
&N_DEV_ReactionNetwork::getReaction(string name)
{
  int ni;
  ni=getReactionNum(name);
  if (ni == -1)
  {
    ostringstream ost;
    ost << " Attempt to access non-existant reaction " << name << endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, ost.str());
  }

  return theReactions[ni];
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setSimpleCalc
// Purpose       : set the named reaction's rate calculator to type Simple
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setSimpleCalc(string name, double k)
{
  getReaction(name).setSimpleRateCalculator(k,C0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setCaptureCalc
// Purpose       : set the named reaction's rate calculator to type Capture
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setCaptureCalc(string name,
                                                      double sigma, double v)
{
  getReaction(name).setCaptureRateCalculator(sigma,v,C0,t0,x0);
}
//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setElectronCaptureCalc
// Purpose       : set the named reaction's rate calculator to type Capture
// Special Notes : Specifically uses v_n=2.3e7 cm/s
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setElectronCaptureCalc(string name,
                                                      double sigma)
{
  getReaction(name).setCaptureRateCalculator(sigma,2.3e7,C0,t0,x0);
}
//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setHoleCaptureCalc
// Purpose       : set the named reaction's rate calculator to type Capture
// Special Notes : uses v_p=1.9e7 cm/s
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setHoleCaptureCalc(string name,
                                                      double sigma)
{
  getReaction(name).setCaptureRateCalculator(sigma,1.9e7,C0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setEmissionCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setEmissionCalc(string name,
                                                        double sigma, double v,
                                                        double N, double E)
{
  getReaction(name).setEmissionRateCalculator(sigma,v,N,E,C0,t0,x0);
}
//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setElectronEmissionCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :Uses v_n=2.3e7 cm/s and N_c=2.86e19 cm^{-3}
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setElectronEmissionCalc(string name,
                                                        double sigma, double E)
{
  getReaction(name).setEmissionRateCalculator(sigma,2.3e7,2.86e19,E,C0,t0,x0);
}
//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setHoleEmissionCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :  uses v_p=1.9e7 cm/s and N_v=2.66e19 cm^{-3}
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setHoleEmissionCalc(string name,
                                                        double sigma, double E)
{
  getReaction(name).setEmissionRateCalculator(sigma,1.9e7,2.66e19,E,C0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setComplexCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setComplexCalc(string name)
{
  getReaction(name).setComplexRateCalculator(species,constants,C0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setDecomplexCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setDecomplexCalc(string name,
                                                         double bindingEnergy,
                                                         double gammaAB,
                                                         double gammaA,
                                                         double gammaB,
                                                         double concSi)
{

  getReaction(name).setDecomplexRateCalculator(species,constants,
                                               bindingEnergy,gammaAB,gammaA,
                                               gammaB,concSi,C0,t0,x0);
}



//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 05/27/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::output(ostream & os) const
{
  int i;

  for (i=0;i< species.size();++i)
  {
    os << "species["<<i<<"] = " << species[i].getName() << endl;
  }
  os << endl;


  for (i=0;i<theReactions.size();++i)
  {
    os << reactionNames[i];
    theReactions[i].output( species, os);
  }

  if (electronCaptureReactions.size() != 0)
  {
    os << "Electron Capture Reactions: " << endl;
    for ( i=0; i<electronCaptureReactions.size(); ++i)
    {
      os << "  Reaction number " << electronCaptureReactions[i] << "("
         << reactionNames[electronCaptureReactions[i]] << ")" << endl;
    }
  }
  if (holeCaptureReactions.size() != 0)
  {
    os << "Hole Capture Reactions: " << endl;
    for ( i=0; i<holeCaptureReactions.size(); ++i)
    {
      os << "  Reaction number " << holeCaptureReactions[i] << "("
         << reactionNames[holeCaptureReactions[i]] << ")" << endl;
    }
  }
  if (electronEmissionReactions.size() != 0)
  {
    os << "Electron Emission Reactions: " << endl;
    for ( i=0; i<electronEmissionReactions.size(); ++i)
    {
      os << "  Reaction number " << electronEmissionReactions[i] << "("
         << reactionNames[electronEmissionReactions[i]] << ")" << endl;
    }
  }
  if (holeEmissionReactions.size() != 0)
  {
    os << "Hole Emission Reactions: " << endl;
    for ( i=0; i<holeEmissionReactions.size(); ++i)
    {
      os << "  Reaction number " << holeEmissionReactions[i] << "("
         << reactionNames[holeEmissionReactions[i]] << ")" << endl;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 05/27/06
//-----------------------------------------------------------------------------
inline ostream & operator<<(ostream & os, const N_DEV_ReactionNetwork & rn)
{
  os << "Reaction Network: " << endl;
  rn.output(os);

  return os;
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getERate
// Purpose       : Compute the total rate at which electrons are "consumed" or
//                 "produced" by all the capture and emission reactions,
//                 if there are any.  This can be used even if the electron
//                 concentration is held fixed (it'll just be the sum of all
//                 reaction rates involving electrons).
// Special Notes : Assumes that all emission and capture reactions are
//                 of the form B=>A+E or A+E=>B and will be incorrect if
//                 any reaction not of this form is given a name with
//                 "_ELECTRON_EMISSION" or "_ELECTRON_CAPTURE" in it.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline double N_DEV_ReactionNetwork::getERate(vector<double> &concs,
                                            vector<double> &constant_vec)
{
  return getRate(concs,constant_vec,electronCaptureReactions,
                 electronEmissionReactions);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getELifetime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
inline double N_DEV_ReactionNetwork::getELifetime(vector<double> &concs,
                                            vector<double> &constant_vec)
{

  double eConc;
  int concNum=getReactantNum("E");
  // Bleah
  if (concNum < 0)
    eConc = constant_vec[-(concNum+1)];
  else
    eConc = concs[concNum];


  return getCaptureLifetime(concs,constant_vec,electronCaptureReactions,
                 eConc);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getELifetimes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::getELifetimes(vector<double> &concs,
                                            vector<double> &constant_vec,
                                            vector<double> &lifetimes)
{

  double eConc;
  int concNum=getReactantNum("E");
  // Bleah
  if (concNum < 0)
    eConc = constant_vec[-(concNum+1)];
  else
    eConc = concs[concNum];


  getCaptureLifetimes(concs,constant_vec,electronCaptureReactions,
                      eConc,lifetimes);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getDERateDC
// Purpose       : compute vector of derivatives of net electron emission
//                 rate (e.g. result of getERate) with respect to concentration
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::getDERateDC(vector<double> &concs,
                                            vector<double> &constant_vec,
                                            vector<double> &dRatedC)
{
  getDRateDC(concs,constant_vec,electronCaptureReactions,
                 electronEmissionReactions,dRatedC);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getDERateDConst
// Purpose       : compute vector of derivatives of net electron emission
//                 rate (e.g. result of getERate) with respect to concentration
//
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 11/15/08
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::getDERateDConst(vector<double> &concs,
                                                        vector<double> &constant_vec,
                                                        vector<double> &dRatedConst)
{
  getDRateDConst(concs,constant_vec,electronCaptureReactions,
                  electronEmissionReactions,dRatedConst);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getHRate
// Purpose       : Compute the total rate at which holes are "consumed" or
//                 "produced" by all the capture and emission reactions,
//                 if there are any.  This can be used even if the electron
//                 concentration is held fixed (it'll just be the sum of all
//                 reaction rates involving electrons).
// Special Notes : Assumes that all emission and capture reactions are
//                 of the form B=>A+H or A+H=>B and will be incorrect if
//                 any reaction not of this form is given a name with
//                 "_HOLE_EMISSION" or "_HOLE_CAPTURE" in it.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline double N_DEV_ReactionNetwork::getHRate(vector<double> &concs,
                                            vector<double> &constant_vec)
{
  return getRate(concs,constant_vec,holeCaptureReactions,
                 holeEmissionReactions);;
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getHLifetime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
inline double N_DEV_ReactionNetwork::getHLifetime(vector<double> &concs,
                                            vector<double> &constant_vec)
{

  double hConc;
  int concNum=getReactantNum("H");
  // Bleah
  if (concNum < 0)
    hConc = constant_vec[-(concNum+1)];
  else
    hConc = concs[concNum];

  return getCaptureLifetime(concs,constant_vec,holeCaptureReactions,
                 hConc);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getHLifetimes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::getHLifetimes(vector<double> &concs,
                                            vector<double> &constant_vec,
                                            vector<double> &lifetimes)
{

  double hConc;
  int concNum=getReactantNum("H");
  // Bleah
  if (concNum < 0)
    hConc = constant_vec[-(concNum+1)];
  else
    hConc = concs[concNum];


  getCaptureLifetimes(concs,constant_vec,holeCaptureReactions,
                      hConc,lifetimes);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getDHRateDC
// Purpose       : compute vector of derivatives of net hole emission
//                 rate (e.g. result of getHRate) with respect to concentration
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::getDHRateDC(vector<double> &concs,
                                            vector<double> &constant_vec,
                                            vector<double> &dRatedC)
{
  getDRateDC(concs,constant_vec,holeCaptureReactions,
                 holeEmissionReactions,dRatedC);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getDHRateDConst
// Purpose       : compute vector of derivatives of net hole emission
//                 rate (e.g. result of getHRate) with respect to concentration
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::getDHRateDConst(vector<double> &concs,
                                                        vector<double> &constant_vec,
                                                        vector<double> &dRatedConst)
{
  getDRateDConst(concs,constant_vec,holeCaptureReactions,
                 holeEmissionReactions,dRatedConst);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getDiffusionCoefficient
// Purpose       :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline double N_DEV_ReactionNetwork::getDiffusionCoefficient
(const std::string & name, const double temp)
{
  int num = getSpeciesNum(name);
  double D = 0.0;

  if (num < 0)
    D = 0.0;
  else
    D = species[num].getDiffusionCoefficient(temp);

  return D;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getDiffusionCoefficient
// Purpose       :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/31/06
//-----------------------------------------------------------------------------
inline double N_DEV_ReactionNetwork::getDiffusionCoefficient
  (int specie, const double temp)
{
  double D = 0.0;

  if (specie < 0)
    D = 0.0;
  else
    D = species[specie].getDiffusionCoefficient(temp);

  return D;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getChargeState
// Purpose       :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 02/22/2012
//-----------------------------------------------------------------------------
inline int N_DEV_ReactionNetwork::getChargeState
(const std::string & name)
{
  int num = getSpeciesNum(name);
  int z = 0;

  if (num < 0)
    z = 0;
  else
    z = species[num].getChargeState();

  return z;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::getChargeState
// Purpose       :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 10/31/06
//-----------------------------------------------------------------------------
inline int N_DEV_ReactionNetwork::getChargeState
  (int specie)
{
  int z = 0;

  if (specie < 0)
    z = 0;
  else
    z = species[specie].getChargeState();

  return z;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ReactionNetwork::setApplySources
// Purpose       :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 04/19/09
//-----------------------------------------------------------------------------
inline void N_DEV_ReactionNetwork::setApplySources(bool flag)
{
  applySources = flag;
}

#endif
