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
// Filename       : $RCSfile: N_DEV_RateConstantCalculators.h,v $
//
// Purpose        : strategy pattern for reactions with different temperature
//                  dependent rate constant styles
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 07/27/2006
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

#ifndef N_DEV_RateConstantCalculators_H
#define N_DEV_RateConstantCalculators_H

// this does a "using namespace std;" so we needn't.
#include <N_UTL_Misc.h>

#include <string>
#include <vector>

namespace N_DEV {
class Specie;
}

namespace Xyce {
namespace Device {

// forward declaration:

enum CalcType
{
  NOCALC,
  SIMPLECALC,
  CAPTURECALC,
  EMISSIONCALC,
  COMPLEXCALC,
  DECOMPLEXCALC
};

// Abstract interface class for "rate calculator" strategy pattern
class RateCalculator
{
public:
  virtual double computeRateConstant(double T) = 0;
  virtual double rateConstantScaleFactor() =  0;
  virtual void setScaleFactors(double C0, double t0, double x0)=0;
  virtual CalcType calcType() =  0;
  virtual RateCalculator *Clone()=0;
  virtual ~RateCalculator() {}

};

// Class for trivial, constant rate constant (independent of temperature)
class SimpleRateCalculator: public RateCalculator
{
public:
  SimpleRateCalculator(double k, double C0, double t0, double x0);
  SimpleRateCalculator(SimpleRateCalculator &right);

  virtual SimpleRateCalculator *Clone();
  virtual double computeRateConstant(double T);
  virtual double rateConstantScaleFactor();
  inline virtual void setScaleFactors(double C0, double t0, double x0)
  {rk0=C0*t0;};

  inline virtual CalcType calcType() {return SIMPLECALC;};
private:
  double K;
  double rk0;
};

// Electron or Hole capture reaction
class CaptureRateCalculator: public RateCalculator
{
public:
  CaptureRateCalculator(double sigma, double v, double C0, double t0,
                        double x0);
  CaptureRateCalculator(CaptureRateCalculator &right);
  virtual CaptureRateCalculator *Clone();

  virtual double computeRateConstant(double T);
  virtual double rateConstantScaleFactor();
  inline virtual void setScaleFactors(double C0, double t0, double x0)
  {rk0=C0*t0;};
  inline virtual CalcType calcType() {return CAPTURECALC;};
private:
  double K;
  double rk0;
};

// Electron or Hole emission reaction
class EmissionRateCalculator: public RateCalculator
{
public:
  EmissionRateCalculator(double sigma, double v, double N, double Energy,
                         double C0, double t0, double x0);
  EmissionRateCalculator(EmissionRateCalculator &right);
  virtual EmissionRateCalculator *Clone();

  virtual double computeRateConstant(double T);
  virtual double rateConstantScaleFactor();
  inline virtual void setScaleFactors(double C0, double t0, double x0)
  {T0=t0;};
  inline virtual CalcType calcType() {return EMISSIONCALC;};
private:
  double K_f;
  double E;
  double T0;
};

// for a irreversible two-species complexing reaction
// A+B->AB
class ComplexRateCalculator: public RateCalculator
{
public:
  ComplexRateCalculator(vector<N_DEV::Specie> &VariableSpecies,
                        vector<N_DEV::Specie> &ConstantSpecies,
                        vector< pair<int,double> > &Reactants,
                        double C0, double t0, double x0);
  ComplexRateCalculator(ComplexRateCalculator &right);
  virtual ComplexRateCalculator *Clone();

  virtual double computeRateConstant(double T);
  virtual double rateConstantScaleFactor();
  inline virtual void setScaleFactors(double C0, double t0, double x0)
  {rk0=C0*t0;};
  inline virtual CalcType calcType() {return COMPLEXCALC;};
private:
  N_DEV::Specie *Specie1,*Specie2;
  double reaction_distance_factor;
  bool Tdep;
  double rk0;
};

// for a irreversible two-species decomplexing reaction
// AB->A+B
class DecomplexRateCalculator: public RateCalculator
{
public:
  DecomplexRateCalculator(vector<N_DEV::Specie> &VariableSpecies,
                          vector<N_DEV::Specie> &ConstantSpecies,
                          vector< pair<int,double> > &Reactants,
                          vector< pair<int,double> > &Products,
                          double bindingEnergy,
                          double degenAB,
                          double degenA,
                          double degenB,
                          double siliconConcentration,
                          double C0, double t0, double x0);
  DecomplexRateCalculator(DecomplexRateCalculator &right);
  virtual DecomplexRateCalculator *Clone();

  virtual double computeRateConstant(double T);
  virtual double rateConstantScaleFactor();
  inline virtual void setScaleFactors(double C0, double t0, double x0)
  {
    // The rate constant includes the concentration of silicon, which
    // means it scales differently.  See comments in
    // N_DEV_RxnRegion::scaleRateConstants for details
    rk0=t0;
    c0=C0;
  };
  inline virtual CalcType calcType() {return DECOMPLEXCALC;};

private:
  N_DEV::Specie *Specie1,*Specie2;
  double reaction_distance_factor;
  bool Tdep;
  double deltaE;
  double gammaA,gammaB,gammaAB;
  double concSi;
  double rk0;
  double c0;
};

} // namespace Device
} // namespace Xyce

#endif
