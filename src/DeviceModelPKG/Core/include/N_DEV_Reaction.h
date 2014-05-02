//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_Reaction.h,v $
//
// Purpose        : Provide a generic class for reactions using simple
//                  law-of-mass-action kinetics, i.e.:
//                       A+2B+C -> D + 3E + F
//                  implies that the reaction happens at a rate:
//                       reactionRate = k*[A][B]^2[C]
//                  where [X] means "concentration of species X"
//
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 03/20/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_DEV_REACTION_H
#define N_DEV_REACTION_H
#include <iosfwd>
#include <vector>

#include <N_DEV_Specie.h>
#include <N_DEV_RateConstantCalculators.h>

namespace Xyce {
namespace Device {

class Reaction
{
public:
  Reaction();
  Reaction(std::vector< std::pair<int,double> > & ,
           std::vector< std::pair<int,double> > &,
           double);
  Reaction(const Reaction &right);
  ~Reaction();
  void setRateConstant(double);
  void setRateConstantFromCalculator(double T);
  void scaleRateConstant(double);
  void scaleRateConstantFromCalculator();
  void unscaleRateConstantFromCalculator();
  void setReactants(std::vector< std::pair<int,double> > &products );
  void addReactant(int species,double stoich);
  void setProducts(std::vector< std::pair<int,double> > &products );
  void addProduct(int species,double stoich);
  double getRate(std::vector<double> &concentrations,
                 std::vector<double> &constants);
  void getDdt(std::vector<double> &concentrations,
              std::vector<double> &constants,
              std::vector<double> &ddt);
  void getDRateDC(std::vector<double> &concentrations,
                  std::vector<double> &constants,
                  std::vector<double> &dratedc);
  void getDRateDConst(int constNum,
                      std::vector<double> &concentrations,
                      std::vector<double> &constants,
                      double &dratedc);
  void getJac(std::vector<double> &concentrations,
              std::vector<double> &constants,
              std::vector<std::vector<double> > &jac);
  void getDFdConst(int constantNumber,
                   std::vector<double> &concentrations,
                   std::vector<double> &constants,
                   std::vector<double> &dFdConst);

  void output ( const std::vector<Specie> & species,
                std::ostream & os ) const;

  void setSimpleRateCalculator(double k, double C0, double t0, double x0);
  void setCaptureRateCalculator(double sigma, double v, double C0, double t0,
                                double x0);
  void setEmissionRateCalculator(double sigma, double v, double N,
                                 double Energy, double C0, double t0,
                                 double x0);
  void setComplexRateCalculator(std::vector<Specie> &VariableSpecies,
                                std::vector<Specie> &ConstantSpecies,
                                double C0, double t0, double x0);
  void setDecomplexRateCalculator(std::vector<Specie> &VariableSpecies,
                                  std::vector<Specie> &ConstantSpecies,
                                  double bindingEnergy,
                                  double gammaAB, double gammaA, double gammaB,
                                  double concSi,
                                  double C0, double t0, double x0);

  inline void setScaleFactors(double C0, double t0, double x0)
  {
    if (myRateCalc)
      myRateCalc->setScaleFactors(C0,t0,x0);
  };

  Reaction & operator=(const Reaction & right);

private:
  void setDependency(int cSize);
  void setConstDependency(int cSize);
  std::vector< std::pair<int,double> > theReactants;
  std::vector< std::pair<int,double> > theProducts;
  double theRateConstant;
  int numconcs; // size of vector of concentrations
  int numconsts; // size of vector of constants
  std::vector<int> concDependency;
  std::vector<int> constDependency;
  RateCalculator *myRateCalc;


};

//-----------------------------------------------------------------------------
// Function      : Reaction::setRateConstant
// Purpose       : Accessor function to set rate constant
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline void Reaction::setRateConstant(double rateConst)
{
  theRateConstant=rateConst;
}

//-----------------------------------------------------------------------------
// Function      : Reaction::scaleRateConstant
// Purpose       : Accessor function to scale rate constant
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline void Reaction::scaleRateConstant(double scalar)
{
  theRateConstant *= scalar;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Reaction N_DEV_Reaction;

#endif
