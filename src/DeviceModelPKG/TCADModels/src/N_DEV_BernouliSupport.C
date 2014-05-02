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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_BernouliSupport.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/01/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.22 $
//
// Revision Date  : $Date: 2014/02/24 23:49:18 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif


#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_BernouliSupport.h>
#include <N_ERH_ErrorMgr.h>

// These are the old Bernouli breakpoint numbers.  They were generated
// using the SGF program brkpnts.c . They are appropriate for linux
// They should be removed soon, as Xyce can generate these internally now.
#define BP0_BERN    -3.742945958815751e+01
#define BP1_BERN    -1.948848145621305e-02
#define BP2_BERN    1.230611609815494e-02
#define BP3_BERN    3.742945958815751e+01
#define BP4_BERN    7.451332191019408e+02
#define BP0_DBERN   -4.117119704160766e+01
#define BP1_DBERN   -3.742945958815751e+01
#define BP2_DBERN   -1.848271746976161e-02
#define BP3_DBERN   8.806697697210611e-03
#define BP4_DBERN   3.742945958815751e+01
#define BP5_DBERN   7.451332191019408e+02
#define BP0_AUX1    -8.301056680276218e-03
#define BP1_AUX1    8.301056680276218e-03
#define BP0_DAUX1   -4.826242066078996e-03
#define BP1_DAUX1   4.826242066078996e-03
#define BP0_AUX2    -4.436141955583643e+01
#define BP1_AUX2    3.680808162809191e+01
#define BP2_AUX2    7.451332191019419e+02
#define BP0_DAUX2   -7.451332191019419e+02
#define BP1_DAUX2   -4.436141955583643e+01
#define BP2_DAUX2   3.680808162809191e+01
#define BP3_DAUX2   7.451332191019419e+02
#define BP0_MISC    7.097827128183643e+02

// Other defines.  MAXDOUBLE is probably set up elsewhere in Xyce - check this.
#define PRECISION 1.0e-15
#define MAX_ITERATIONS 100

#ifndef MAXDOUBLE
#define MAXDOUBLE   1.797693E+308
#endif

namespace Xyce {
namespace Device {

namespace {

// ERK, 8/1/2004:  I would prefer to just put these functions inside of the
// Bernouli support
// class, but the functions Bisection, Asymtotic and Secant all rely on
// function pointers, and function pointers are a lot easier to use if the
// functions they point to are old-fashioned C functions.
//
// The brkpnts.c program, which is in C, rather than C++ had a very C-oriented
// style.  I haven't had time to move this stuff closer to a C++ style.

double Bbp0a(double x)  { return exp(x) - 1.0; }
double Bbp0b(double x)  { return - 1.0; }

double Bbp1a(double x)  { return x / (exp(x) - 1.0); }
double Bbp1b(double x)
{
  return 1.0 - x/2.0 * (1.0 - x/6.0 * (1.0 - x*x/60.0));
}

double Bbp2a(double x)
{
  return 1.0 - x/2.0 * (1.0 - x/6.0 * (1.0 - x*x/60.0));
}

double Bbp2b(double x)
{
  return x * exp(-x) / (1.0 - exp(-x));
}

double Bbp3a(double x)  { return 1.0 - exp(-x); }
double Bbp3b(double x)  { return 1.0; }

double Bbp4a(double x)  { return x * exp(-x); }
double Bbp4b(double x)  { return 0.0; }

double dBbp0a(double x) { return (1.0 - x) * exp(x) - 1.0; }
double dBbp0b(double x) { return -1.0; }

double dBbp2a(double x)
{
  return ((1.0 - x) * exp(x) - 1.0) / ((exp(x) - 1.0) * (exp(x) - 1.0));
}

double dBbp2b(double x) { return -0.5 + x/6.0 * (1.0 - x*x/30.0); }

double dBbp3a(double x) { return -0.5 + x/6.0 * (1.0 - x*x/30.0); }
double dBbp3b(double x)
{
  return (exp(-x)*(1.0 - x) - exp(-2.0 * x))/((1.0 - exp(-x))*(1.0 - exp(-x)));
}

double dBbp5a(double x) { return exp(-x) * (1.0 - x) - exp(-2.0 * x); }
double dBbp5b(double x) { return 0.0; }

double AUX1bp0a(double x) { return x / sinh(x); }
double AUX1bp0b(double x) { return 1.0 - x*x/6.0 * (1.0 - 7.0*x*x/60.0); }

double dAUX1bp0a(double x)
{
  return (sinh(x) - x*cosh(x)) / (sinh(x) * sinh(x));
}

double dAUX1bp0b(double x) { return -x/3.0 * (1.0 - 7*x*x/30.0); }

double AUX2bp0a(double x) { return 1.0; }
double AUX2bp0b(double x) { return 1.0 + exp(x); }

double AUX2bp1a(double x) { return 1.0 + exp(x); }
double AUX2bp1b(double x) { return exp(x); }

double AUX2bp2a(double x) { return exp(-x); }
double AUX2bp2b(double x) { return 0.0; }

double dAUX2bp0a(double x) { return exp(x); }
double dAUX2bp0b(double x) { return 0.0; }

} // namespace <empty>

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::BernouliSupport
// Purpose       : constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
BernouliSupport::BernouliSupport ():
  // use the numbers generated by brkpnts.c  (for a single platform, linux)
  // bernouli function breakpoints:
  bp0_BERN   ( BP0_BERN ),
  bp1_BERN   ( BP1_BERN ),
  bp2_BERN   ( BP2_BERN ),
  bp3_BERN   ( BP3_BERN ),
  bp4_BERN   ( BP4_BERN ),

  // bernouli derivative function breakpoints:
  bp0_DBERN  ( BP0_DBERN ),
  bp1_DBERN  ( BP1_DBERN ),
  bp2_DBERN  ( BP2_DBERN ),
  bp3_DBERN  ( BP3_DBERN ),
  bp4_DBERN  ( BP4_DBERN ),
  bp5_DBERN  ( BP5_DBERN ),

  // aux1 function breakpoints:
  bp0_AUX1   ( BP0_AUX1 ),
  bp1_AUX1   ( BP1_AUX1 ),

  // aux1 derivative function breakpoints:
  bp0_DAUX1  ( BP0_DAUX1 ),
  bp1_DAUX1  ( BP1_DAUX1 ),

  // aux2 function breakpoints:
  bp0_AUX2   ( BP0_AUX2 ),
  bp1_AUX2   ( BP1_AUX2 ),
  bp2_AUX2   ( BP2_AUX2 ),

  // aux2 derivative function breakpoints:
  bp0_DAUX2  ( BP0_DAUX2 ),
  bp1_DAUX2  ( BP1_DAUX2 ),
  bp2_DAUX2  ( BP2_DAUX2 ),
  bp3_DAUX2  ( BP3_DAUX2 ),

  // This is supposed to be the log of MAXDOUBLE.
  bp0_MISC   ( BP0_MISC )
{
  // using the hardwired numbers b/c the calculated numbers do not work
  // on the alpha.
#define Xyce_NEW_BERN_CALCULATION 1
#ifndef Xyce_NEW_BERN_CALCULATION
  // generate breakpoints internally.  This is better as it will be specific to
  // the current platform.

#if 0
  // Bernouli breakpoints:
  // these are commented out b/c they aren't needed, and they seem to
  // cause problems on OSX.  I don't know why...
  bp0_BERN  = Asymptotic(Bbp0a, Bbp0b,  0.0e+00, -1.0e+02);
  bp1_BERN  = Secant    (Bbp1a, Bbp1b, -1.0e+00);
  bp2_BERN  = Secant    (Bbp2a, Bbp2b, +1.0e+00);
  bp3_BERN  = Asymptotic(Bbp3a, Bbp3b,  0.0e+00, +1.0e+02);
  bp4_BERN  = Asymptotic(Bbp4a, Bbp4b,  1.0e+00, +1.0e+03);

  // Bernouli derivative breakpoints:
  bp0_DBERN = Asymptotic(dBbp0a, dBbp0b,  0.0e+00, -1.0e+02);
  bp1_DBERN = bp0_BERN;
  bp2_DBERN = Secant    (dBbp2a, dBbp2b, -1.0e+00);
//bp2_DBERN = Bisection (dBbp2a, dBbp2b, dB[1], -1.0e-6);
  bp3_DBERN = Secant    (dBbp3a, dBbp3b, +1.0e+00);
  bp4_DBERN = bp3_BERN;
  bp5_DBERN = Asymptotic(dBbp5a, dBbp5b,  1.0e+00, +1.0e+03);
#endif

  // Aux1 breakpoints:
  bp0_AUX1 = Secant(AUX1bp0a, AUX1bp0b, -1.0e+00);
  bp1_AUX1 = Secant(AUX1bp0a, AUX1bp0b,  1.0e+00);

  // Aux1 derivative breakpoints:
  bp0_DAUX1 = Secant(dAUX1bp0a, dAUX1bp0b, -1.0e+00);
  bp1_DAUX1 = Secant(dAUX1bp0a, dAUX1bp0b,  1.0e+00);

  // Aux2 breakpoints:
  bp0_AUX2 = Asymptotic(AUX2bp0a, AUX2bp0b,  0.0e+00, -1.0e+02);
  bp1_AUX2 = Asymptotic(AUX2bp1a, AUX2bp1b,  0.0e+00,  1.0e+02);
  bp2_AUX2 = Asymptotic(AUX2bp2a, AUX2bp2b,  0.0e+00,  1.0e+02);

  // Aux2 derivative breakpoints:
  bp0_DAUX2 = Asymptotic(dAUX2bp0a, dAUX2bp0b,  0.0e+00, -1.0e+02);
  bp1_DAUX2 = bp0_AUX2;
  bp2_DAUX2 = bp1_AUX2;
  bp3_DAUX2 = bp2_AUX2;

  // miscellaneous
  bp0_MISC = log(MAXDOUBLE);
#endif

#ifdef Xyce_DEBUG_DEVICE
#if 0
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
  Xyce::dout() << std::endl;
  Xyce::dout() << "Bernouli function breakpoints: " <<endl;
  Xyce::dout() << "bp0_AUX1  = " << bp0_AUX1  << std::endl;
  Xyce::dout() << "bp1_AUX1  = " << bp1_AUX1  << std::endl;
  Xyce::dout() << "bp0_DAUX1 = " << bp0_DAUX1 << std::endl;
  Xyce::dout() << "bp1_DAUX1 = " << bp1_DAUX1 << std::endl;

  Xyce::dout() << "bp0_AUX2  = " << bp0_AUX2  << std::endl;
  Xyce::dout() << "bp1_AUX2  = " << bp1_AUX2  << std::endl;
  Xyce::dout() << "bp2_AUX2  = " << bp2_AUX2  << std::endl;
  Xyce::dout() << "bp0_DAUX2 = " << bp0_DAUX2 << std::endl;
  Xyce::dout() << "bp1_DAUX2 = " << bp1_DAUX2 << std::endl;
  Xyce::dout() << "bp2_DAUX2 = " << bp2_DAUX2 << std::endl;
  Xyce::dout() << "bp3_DAUX2 = " << bp3_DAUX2 << std::endl;
  Xyce::dout() << "bp0_MISC  = " << bp0_MISC << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif
#endif
}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::~BernouliSupport
// Purpose       : constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
BernouliSupport::~BernouliSupport ()
{

}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::BernouliSupport
// Purpose       : copy constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
BernouliSupport::BernouliSupport
  (const BernouliSupport & right)
{

}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::BernouliSupport
// Purpose       : copy constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
int BernouliSupport::sign(double x)
{
  if      (x < 0.0) return(-1);
  else if (x > 0.0) return(+1);
  else              return(0);
}

//-----------------------------------------------------------------------------
// Function      : double BernouliSupport::Bisection
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
double BernouliSupport::Bisection
  (FUNC func1, FUNC func2, double Xpos, double Xneg)
{
  double Fpos = func1(Xpos) - func2(Xpos);
  double Fneg = func1(Xneg) - func2(Xneg);
  double Xmid, Fmid, Xlast;

  if    (Fpos == 0.0) return(Xpos);
  else if (Fneg == 0.0) return(Xneg);
  else if ((Fpos > 0.0) && (Fneg < 0.0)) ;
  else if ((Fpos < 0.0) && (Fneg > 0.0))
  {
    Xmid = Xpos;
    Xpos = Xneg;
    Xneg = Xmid;
  }
  else
  {
    std::string msg = "BernouliSupport::Bisection ";
    msg += " Initial interval may not contain a root";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  Xlast = 0.0;
  do
  {
    Xmid = 0.5 * (Xpos + Xneg);
    Fmid = func1(Xmid) - func2(Xmid);
    if      (Fmid < 0.0)  Xneg = Xmid;
    else if (Fmid > 0.0)  Xpos = Xmid;
    if (Xlast == Xmid) return(Xmid);
    else Xlast = Xmid;
  } while (Xneg != Xpos);

  return(Xmid);

}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::Secant
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
double BernouliSupport::Secant(FUNC func1, FUNC func2, double x1)
{
  double slope, dx, x3, f3;
  int    s3, iteration;

  double x2 = 0.9 * x1;
  double f1 = func1(x1) - func2(x1);
  double f2 = func1(x2) - func2(x2);
  int    s2 = sign(x2);

  for(;;)
  {
    iteration = 0;
    slope = (f2 - f1) / (x2 - x1);
    dx = f2 / slope;
    x3 = x2 - dx;
    f3 = func1(x3) - func2(x3);
    s3 = sign(x3);

    while ((fabs(f3) >= fabs(f2)) || (s3 != s2))
    {
      dx /= 2.0;
      x3 += dx;
      f3  = func1(x3) - func2(x3);
      s3  = sign(x3);
      if (++iteration > MAX_ITERATIONS)
      {
        if (fabs(f2) <= 100.0 * PRECISION) return(x2);
	std::string msg = "BernouliSupport::Secant ";
	msg += " method not converging.";
	N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
      }
    }

    x1 = x2;
    x2 = x3;
    f1 = f2;
    f2 = f3;

    if ((fabs(dx / x2) <= PRECISION) || (fabs(f2) <= PRECISION)) break;
  }

  return(x3);
}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::Asymptotic
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
double BernouliSupport::Asymptotic(FUNC func1, FUNC func2, double x, double dx)
{
  double test = 1.0;
  while (1)
  {

    if (x==0.0)  test = 1.0;
    else         test = fabs(dx/x);
    if (test <= PRECISION) return(x);

    while (func1(x) != func2(x))
      x += dx;
    dx *= -0.1;

    if (x==0.0)  test = 1.0;
    else         test = fabs(dx/x);
    if (test <= PRECISION) return(x);

    while (func1(x) == func2(x))
      x += dx;
    dx *= -0.1;
  }
}

} // namespace Device
} // namespace Xyce
