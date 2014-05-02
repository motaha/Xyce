// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_RateConstantCalculators.C,v $
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
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>



#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <string>
#include <vector>

#include <N_DEV_Const.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Specie.h>
#include <N_DEV_RateConstantCalculators.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::SimpleRateCalculator
// Purpose       : constructor for "simple" reaction rate calculator
//                 This is one that has no temperature dependence, and scales
//                 as concentration*time
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  SimpleRateCalculator::SimpleRateCalculator(double k, double C0, double t0,
                                             double x0)
    : K(k)
  {
    setScaleFactors(C0,t0,x0);
  }

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::SimpleRateCalculator
// Purpose       : Copy constructor for "simple" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  SimpleRateCalculator::SimpleRateCalculator(SimpleRateCalculator &right)
    :K(right.K),
     rk0(right.rk0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::Clone
// Purpose       : Copy self operation for "simple" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  SimpleRateCalculator *SimpleRateCalculator::Clone()
  {
    return new SimpleRateCalculator(*this);
  }

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::computeRateConstant
// Purpose       : returns rate constant for simple rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double SimpleRateCalculator::computeRateConstant(double T)
  {
    return(K);
  }

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for simple rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double SimpleRateCalculator::rateConstantScaleFactor()
  {
    return (rk0);
  }


//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::CaptureRateCalculator
// Purpose       : constructor for capture reaction rate calculator
//                 This is one that models electron or hole capture reactions.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  CaptureRateCalculator::CaptureRateCalculator(double sigma, double v,
                                               double C0, double t0,
                                               double x0)
  {
    K=sigma*v;
    setScaleFactors(C0,t0,x0);
  }

//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::CaptureRateCalculator
// Purpose       : Copy constructor for "capture" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  CaptureRateCalculator::CaptureRateCalculator(CaptureRateCalculator &right)
    :K(right.K),
     rk0(right.rk0)
  {
  }
//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::Clone
// Purpose       : Copy self operation for "capture" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  CaptureRateCalculator *CaptureRateCalculator::Clone()
  {
    return new CaptureRateCalculator(*this);
  }

//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::computeRateConstant
// Purpose       : returns rate constant for capture rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double CaptureRateCalculator::computeRateConstant(double T)
  {
    return(K);
  }

//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for capture rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double CaptureRateCalculator::rateConstantScaleFactor()
  {
    return (rk0);
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::EmissionRateCalculator
// Purpose       : constructor for emission reaction rate calculator
//                 This is one that models electron or hole emission reactions.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  EmissionRateCalculator::EmissionRateCalculator(double sigma, double v,
                                               double N, double Energy,
                                               double C0, double t0,
                                               double x0)
    : E(Energy)
  {
    K_f=sigma*v*N;
    setScaleFactors(C0,t0,x0);
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::EmissionRateCalculator
// Purpose       : Copy constructor for "emission" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  EmissionRateCalculator::EmissionRateCalculator(EmissionRateCalculator &right)
    :T0(right.T0),
     E(right.E),
     K_f(right.K_f)
  {
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::Clone
// Purpose       : Copy self operation for "emission" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  EmissionRateCalculator *EmissionRateCalculator::Clone()
  {
    return new EmissionRateCalculator(*this);
  }
//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::computeRateConstant
// Purpose       : returns rate constant for emission rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double EmissionRateCalculator::computeRateConstant(double T)
  {
    double KbT=CONSTboltz*T/CONSTQ;
    return(K_f*exp(-E/KbT));
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for emission rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double EmissionRateCalculator::rateConstantScaleFactor()
  {
    return (T0);
  }


//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::ComplexRateCalculator
// Purpose       : constructor for "complex" reaction rate calculator
//                 This is one that models two discrete species forming a
//                 complex, e.g.:
//                       V0 + VM -> VVM
//
// Special Notes : For this to work, there must either be two species in the
//                 reactants list, or one species with 2.0 as the stochiometric
//                 coefficient.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  ComplexRateCalculator::ComplexRateCalculator(
    std::vector<Specie> &VariableSpecies, std::vector<Specie> &ConstantSpecies,
     std::vector< std::pair<int,double> > &Reactants,
     double C0, double t0, double x0)
  {
    int ij;

    // Check assumptions:
    if ( ! ((Reactants.size() == 1 && Reactants[0].second == 2.0) ||
            (Reactants.size() == 2 && Reactants[0].second == 1.0 &&
             Reactants[1].second == 1.0)))
    {
      std::string msg;
      msg = "ComplexRateCalculator: Invalid attempt to use complex rate method.";
      msg = "  This method is only valid for binary complexing reactions:\n";
      if (Reactants.size() == 1)
      {
        msg += "   Only one reactant specified, but its stoichimetric coefficient is not 2.\n";
      }
      else if (Reactants.size() == 2)
      {
        msg += "   Two reactants specified, but both stoichimetric coefficient are not 1.\n";
      }
      else
      {
        msg += "   More than two reactants specified.\n";
      }
      // exit with error --- we probably need to do something more careful
      // in parallel
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
    }

    if (Reactants[0].first >= 0)
      Specie1 = &(VariableSpecies[Reactants[0].first]);
    else
      Specie1 = &(ConstantSpecies[-(Reactants[0].first+1)]);

    // Handle case where there's only one species with a coefficient of 2.0
    if (Reactants.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
    else
    {
      if (Reactants[1].first >= 0)
        Specie2 = &(VariableSpecies[Reactants[1].first]);
      else
        Specie2 = &(ConstantSpecies[-(Reactants[1].first+1)]);
    }
    ij    =Specie1->getChargeState();
    ij   *=Specie2->getChargeState();



    // Only divide reaction_distance_factor by T in one special case
    Tdep=false;
    if (ij>0)
      reaction_distance_factor=0.0;
    else if (ij == 0)
      reaction_distance_factor = 4*M_PI*5e-8;
    else
    {
      reaction_distance_factor = 4*M_PI*1.4e-4*(-ij);
      Tdep=true;
    }

    setScaleFactors(C0,t0,x0);
  }


//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::ComplexRateCalculator
// Purpose       : Copy constructor for "complex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  ComplexRateCalculator::ComplexRateCalculator(ComplexRateCalculator &right)
    :rk0(right.rk0),
     reaction_distance_factor(right.reaction_distance_factor),
     Tdep(right.Tdep),
     Specie1(right.Specie1),
     Specie2(right.Specie2)
  {
  }

//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::Clone
// Purpose       : Copy self operation for "complex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  ComplexRateCalculator *ComplexRateCalculator::Clone()
  {
    return new ComplexRateCalculator(*this);
  }

//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::computeRateConstant
// Purpose       : returns rate constant for complex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double ComplexRateCalculator::computeRateConstant(double T)
  {
    if (Tdep)
      return(reaction_distance_factor/T*(Specie1->getDiffusionCoefficient(T)
                                        +Specie2->getDiffusionCoefficient(T)));
    else
      return(reaction_distance_factor*(Specie1->getDiffusionCoefficient(T)
                                       +Specie2->getDiffusionCoefficient(T)));

  }

//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for complex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double ComplexRateCalculator::rateConstantScaleFactor()
  {
   return (rk0);
  }





//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::DecomplexRateCalculator
// Purpose       : constructor for "complex" reaction rate calculator
//                 This is one that models two discrete species decomposing from
//                 a complex, e.g.:
//                       VMM->V0 + VM
//
// Special Notes : For this to work, there must either be two species in the
//                 products list, or one species with 2.0 as the stochiometric
//                 coefficient.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 5/04/09
//-----------------------------------------------------------------------------
  DecomplexRateCalculator::DecomplexRateCalculator(
     std::vector<Specie> &VariableSpecies, std::vector<Specie> &ConstantSpecies,
     std::vector< std::pair<int,double> > &Reactants,
     std::vector< std::pair<int,double> > &Products,
     double bindingEnergy, double degenAB, double degenA, double degenB,
     double siliconConcentration,
     double C0, double t0, double x0)
    :concSi(siliconConcentration),
     c0(C0),
     gammaAB(degenAB),
     gammaA(degenA),
     gammaB(degenB),
     deltaE(bindingEnergy)
  {
    int ij;

    // Check assumptions:
    if ( ! ((Products.size() == 1 && Products[0].second == 2.0) ||
            (Products.size() == 2 && Products[0].second == 1.0 &&
             Products[1].second == 1.0)))
    {
      std::string msg;
      msg = "DeomplexRateCalculator: Invalid attempt to use decomplex rate method.";
      msg = "  This method is only valid for decomplexing reactions with two products:\n";
      if (Products.size() == 1)
      {
        msg += "   Only one product specified, but its stoichimetric coefficient is not 2.\n";
      }
      else if (Products.size() == 2)
      {
        msg += "   Two products specified, but both stoichimetric coefficient are not 1.\n";
      }
      else
      {
        msg += "   More than two products specified.\n";
      }
      // exit with error --- we probably need to do something more careful
      // in parallel
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
    }

    if (Products[0].first >= 0)
      Specie1 = &(VariableSpecies[Products[0].first]);
    else
      Specie1 = &(ConstantSpecies[-(Products[0].first+1)]);

    // Handle case where there's only one species with a coefficient of 2.0
    if (Products.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
    else
    {
      if (Products[1].first >= 0)
        Specie2 = &(VariableSpecies[Products[1].first]);
      else
        Specie2 = &(ConstantSpecies[-(Products[1].first+1)]);
    }
    ij    =Specie1->getChargeState();
    ij   *=Specie2->getChargeState();



    // Only divide reaction_distance_factor by T in one special case
    Tdep=false;
    if (ij>0)
      reaction_distance_factor=0.0;
    else if (ij == 0)
      reaction_distance_factor = 4*M_PI*5e-8;
    else
    {
      reaction_distance_factor = 4*M_PI*1.4e-4*(-ij);
      Tdep=true;
    }

    setScaleFactors(C0,t0,x0);
  }


//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::DecomplexRateCalculator
// Purpose       : Copy constructor for "decomplex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 5/04/09
//-----------------------------------------------------------------------------
  DecomplexRateCalculator::DecomplexRateCalculator(DecomplexRateCalculator &right)
    :Specie1(right.Specie1),
     Specie2(right.Specie2),
     reaction_distance_factor(right.reaction_distance_factor),
     Tdep(right.Tdep),
     deltaE(right.deltaE),
     gammaA(right.gammaA),
     gammaB(right.gammaB),
     gammaAB(right.gammaAB),
     concSi(right.concSi),
     rk0(right.rk0),
     c0(right.c0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::Clone
// Purpose       : Copy self operation for "decomplex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  DecomplexRateCalculator *DecomplexRateCalculator::Clone()
  {
    return new DecomplexRateCalculator(*this);
  }


//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::computeRateConstant
// Purpose       : returns rate constant for decomplex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double DecomplexRateCalculator::computeRateConstant(double T)
  {
    double R;
    double KbT=CONSTboltz*T/CONSTQ;
    double k;
    double D1=Specie1->getDiffusionCoefficient(T);
    double D2=Specie2->getDiffusionCoefficient(T);
    if (Tdep)
      R=(reaction_distance_factor/T);
    else
      R=reaction_distance_factor;

    k=(R*(D1 +D2)*(concSi)*((gammaA*gammaB)/gammaAB)*exp(-deltaE/KbT));

    return k;
  }

//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for decomplex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double DecomplexRateCalculator::rateConstantScaleFactor()
  {
   return (rk0);
  }

} // namespace Device
} // namespace Xyce
