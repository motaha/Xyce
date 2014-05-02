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
// Revision Number: $Revision: 1.8.2.3 $
//
// Revision Date  : $Date: 2014/03/01 18:39:05 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_DEV_RateConstantCalculators_H
#define N_DEV_RateConstantCalculators_H

#include <N_UTL_Misc.h>

#include <string>
#include <vector>

namespace Xyce {
namespace Device {

// forward declaration:
class Specie;

enum CalcType
  {
    NOCALC,
    SIMPLECALC,
    CAPTURECALC,
    EMISSIONCALC,
    COMPLEXCALC,
    DECOMPLEXCALC
  };

/// \brief Abstract interface class for "rate calculator" strategy pattern
///
/// Each reaction in a reaction network has a rate constant calculator
/// associated with it.  The interface is defined by the RateCalculator
/// class.
class RateCalculator
{
public:
  virtual double computeRateConstant(double T) = 0; ///< return rate constant at given temperature
  virtual double rateConstantScaleFactor() =  0; ///< return current scale factor for rate constant
  virtual void setScaleFactors(double C0, double t0, double x0)=0; ///< set concentration, time, and space scale factors (space scale factor is currently unused by any calculator)
  virtual CalcType calcType() =  0; ///< return the type of this calculator
  virtual RateCalculator *Clone()=0; ///< create a copy of this calculator
  virtual ~RateCalculator() {}

};

/// \brief Class for trivial, constant rate constant (independent of temperature)
///
/// This is the most basic rate calculator that returns a constant value 
/// provided in the reaction network input file.  
///
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

/// \brief Rate constant calculator for Electron or Hole capture reaction 
///
/// These reactions are of the form R + E -> P; electron_capture(\f$\sigma\f$)
/// or R+H->P; hole_capture(\f$\sigma\f$). 
/// The reaction rate is then \f$\sigma \times v\f$, where v is 2.3e7 for 
/// electron capture and 1.9e7 for holes.
///

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

/// \brief Rate constant calculator for Electron or Hole emission reaction 
///
/// These reactions are of the form R -> P+E; electron_emission(\f$\sigma\f$,Energy)
/// or R->P+H; hole_emission(\f$\sigma\f$,Energy). 
/// The reaction rate is then \f$K_f*exp(E/K_bT)\f$, with \f$K_f=\sigma*v*N\f$.
/// N is 2.86e19 for electron emission and 2.66e19 for hole emission.
/// v is 2.3e7 for electron emission and 1.9e7 for hole emission.
///
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

/// \brief Rate constant calculator for  a irreversible two-species complexing reaction
///
/// These are reactions of the form A+B->AB.
/// There must be two separate species on the reactants side, or one specie
/// with a stoichiometric coefficient of exactly 2.0.
///
/// The rate constant formula depends on the product of the charge states of 
/// the two species, which we will call \f$ij\f$.
///  - If \f$ij > 0\f$: the reaction rate is zero
///  - If \f$ij == 0\f$: The reaction rate is \f$\frac{4\pi\times 5\times10^{-8}}{D_1D_2}\f$ where \f$D_n\f$ is the diffusion coefficient of species \f$n\f$.
///  - If \f$ij < 0\f$: The reaction rate is \f$\frac{4\pi\times 5\times10^{-8}}{TD_1D_2}\f$ where \f$D_n\f$ is the diffusion coefficient of species \f$n\f$.
///
/// Note that the diffusion coefficients of species are temperature dependent.
///
class ComplexRateCalculator: public RateCalculator
{
public:
  ComplexRateCalculator(std::vector<Specie> &VariableSpecies,
                        std::vector<Specie> &ConstantSpecies,
                        std::vector< std::pair<int,double> > &Reactants,
                        double C0, double t0, double x0);
  ComplexRateCalculator(ComplexRateCalculator &right);
  virtual ComplexRateCalculator *Clone();

  virtual double computeRateConstant(double T);
  virtual double rateConstantScaleFactor();
  inline virtual void setScaleFactors(double C0, double t0, double x0)
  {rk0=C0*t0;};
  inline virtual CalcType calcType() {return COMPLEXCALC;};
private:
  Specie *Specie1,*Specie2;
  double reaction_distance_factor;
  bool Tdep;
  double rk0;
};

// for a irreversible two-species decomplexing reaction
// AB->A+B
class DecomplexRateCalculator: public RateCalculator
{
public:
  DecomplexRateCalculator(std::vector<Specie> &VariableSpecies,
                          std::vector<Specie> &ConstantSpecies,
                          std::vector< std::pair<int,double> > &Reactants,
                          std::vector< std::pair<int,double> > &Products,
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
  Specie *Specie1,*Specie2;
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
