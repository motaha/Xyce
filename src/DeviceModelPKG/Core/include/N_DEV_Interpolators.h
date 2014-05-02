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
// Filename       : $RCSfile: N_DEV_Interpolators.h,v $
//
// Purpose        : Interpolator classes
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/31/12
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Interpolators_h
#define Xyce_N_DEV_Interpolators_h

#include <Sacado.hpp>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : interpolator base class
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
class interpolator
{
public:
  interpolator (){};

  virtual void clear (){};

  virtual void init (const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya){};

  virtual void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y){};

  virtual void eval_deriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx){};

  virtual void eval_deriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp){};

  virtual void eval_integ (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result){};

  inline size_t
  binarySearch(
     const std::vector<ScalarT> & xa,
     const ScalarT & x,
     size_t index_lo,
     size_t index_hi);

  inline ScalarT
  integ_eval (
     const ScalarT & ai, const ScalarT & bi, const ScalarT & ci,
     const ScalarT & di, const ScalarT & xi, const ScalarT & a,
     const ScalarT & b);
};

//-----------------------------------------------------------------------------
// Function      : interpolator<ScalarT>::binarySearch
//
// Purpose       : Perform a binary search of an array of values.
//
// Special Notes : adapted from GSL, version 1.15:
//
// The parameters index_lo and index_hi provide an initial bracket,
// and it is assumed that index_lo < index_hi. The resulting index
// is guaranteed to be strictly less than index_hi and greater than
// or equal to index_lo, so that the implicit bracket [index, index+1]
// always corresponds to a region within the implicit value range of
// the value array.
//
// Note that this means the relationship of 'x' to xa[index]
// and xa[index+1] depends on the result region, i.e. the
// behaviour at the boundaries may not correspond to what you
// expect. We have the following complete specification of the
// behaviour.
// Suppose the input is xa[] = { x0, x1, ..., xN }
//    if ( x == x0 )           then  index == 0
//    if ( x > x0 && x <= x1 ) then  index == 0, and sim. for other interior pts
//    if ( x == xN )           then  index == N-1
//    if ( x > xN )            then  index == N-1
//    if ( x < x0 )            then  index == 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
inline size_t
interpolator<ScalarT>::binarySearch(
   const std::vector<ScalarT> & xa,
   const ScalarT & x,
   size_t index_lo,
   size_t index_hi)
{
  size_t ilo = index_lo;
  size_t ihi = index_hi;
  while(ihi > ilo + 1)
  {
    size_t i = (ihi + ilo)/2;
    if(xa[i] > x)
    {
      ihi = i;
    }
    else
    {
      ilo = i;
    }
  }
  return ilo;
}

//-----------------------------------------------------------------------------
// Function      : interpolator<ScalarT>::integ_eval
//
// Purpose       : function for doing the spline integral evaluation
//                 which is common to both the cspline and akima methods
//
// Special Notes : adapted from GSL, version 1.15
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
inline ScalarT
interpolator<ScalarT>::integ_eval (
   const ScalarT & ai,
   const ScalarT & bi,
   const ScalarT & ci,
   const ScalarT & di,
   const ScalarT & xi,
   const ScalarT & a,
   const ScalarT & b)
{
  const ScalarT r1 = a - xi;
  const ScalarT r2 = b - xi;
  const ScalarT r12 = r1 + r2;
  const ScalarT bterm = 0.5 * bi * r12;
  const ScalarT cterm = (1.0 / 3.0) * ci * (r1 * r1 + r2 * r2 + r1 * r2);
  const ScalarT dterm = 0.25 * di * r12 * (r1 * r1 + r2 * r2);
  return (b - a) * (ai + bterm + cterm + dterm);
}

//-----------------------------------------------------------------------------
// Class         : akima spline class
// Purpose       :
// Special Notes : adapted from the GNU Scientific library (GSL) version 1.15.
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
class akima: public interpolator<ScalarT>
{
public:
  akima () {};

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya);

  void clear () { b.clear(); c.clear(); d.clear(); _m.clear(); };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y);

  void eval_deriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx);

  void eval_deriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp);

  void eval_integ (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result);

  void calc (
     const std::vector<ScalarT> & xa,
     std::vector<ScalarT> & b,
     std::vector<ScalarT> & c,
     std::vector<ScalarT> & d,
     std::vector<ScalarT> & m);

public:
  std::vector<ScalarT>  b;
  std::vector<ScalarT>  c;
  std::vector<ScalarT>  d;
  std::vector<ScalarT>  _m;
};

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::calc
// Purpose       :
//
// Special Notes : adapted from GSL, version 1.15
//
// Note that the 'm' indices are offset by 2 on either side of the array,
// compared with the x indices.  ie, m[2] corresponds to x[0], etc.
//
// The original GSL implementation used indicies starting at -2 for the
// m array, by having the local pointer start at +2.  This was to accomodate
// boundary conditions.  I didn't bother with this as I found it confusing,
// and slightly harder to deal with correctly for STL vectors.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::calc (
   const std::vector<ScalarT> & xa,
   std::vector<ScalarT> & b,
   std::vector<ScalarT> & c,
   std::vector<ScalarT> & d,
   std::vector<ScalarT> & m)
{
  size_t i;
  size_t size = xa.size();

  for (i = 0; i < (size - 1); i++)
  {
    const ScalarT NE = fabs (m[i + 3] - m[i+2]) + fabs (m[i + 1] - m[i]);
    if (NE == 0.0)
    {
      b[i] = m[i+2];
      c[i] = 0.0;
      d[i] = 0.0;
    }
    else
    {
      const ScalarT h_i = xa[i + 1] - xa[i];
      const ScalarT NE_next = fabs (m[i + 4] - m[i + 3]) + fabs (m[i+2] - m[i + 1]);
      const ScalarT alpha_i = fabs (m[i + 1] - m[i]) / NE;
      ScalarT alpha_ip1;
      ScalarT tL_ip1;
      if (NE_next == 0.0)
      {
        tL_ip1 = m[i+2];
      }
      else
      {
        alpha_ip1 = fabs (m[i+2] - m[i + 1]) / NE_next;
        tL_ip1 = (1.0 - alpha_ip1) * m[i+2] + alpha_ip1 * m[i + 3];
      }
      b[i] = (1.0 - alpha_i) * m[i + 1] + alpha_i * m[i+2];
      c[i] = (3.0 * m[i+2] - 2.0 * b[i] - tL_ip1) / h_i;
      d[i] = (b[i] + tL_ip1 - 2.0 * m[i+2]) / (h_i * h_i);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::init
// Purpose       :
// Special Notes : adapted from GSL, version 1.15
//
// Note that the 'm' indices are offset by 2 on either side of the array,
// compared with the x indices.  ie, m[2] corresponds to x[0], etc.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::init (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya)
{
  size_t size = xa.size();

  if (b.size() != size) b.resize(size);
  if (c.size() != size) c.resize(size);
  if (d.size() != size) d.resize(size);
  if (_m.size() != size+4) _m.resize(size+4);

  for (int i = 0; i <= size - 2; i++)
  {
    _m[i+2] = (ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]);
  }

  // non-periodic boundary conditions
  _m[0] = 3.0 * _m[2] - 2.0 * _m[3];
  _m[1] = 2.0 * _m[2] - _m[3];
  _m[size + 1] = 2.0 * _m[size] - _m[size-1];
  _m[size + 2] = 3.0 * _m[size] - 2.0 * _m[size-1];

  calc (xa, b, c, d, _m);
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::eval
// Purpose       :
// Special Notes : adapted from GSL, version 1.15
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::eval (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & y)
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  const ScalarT x_lo = xa[index];
  const ScalarT delx = x - x_lo;
  y = ya[index] + delx * (b[index] + delx * (c[index] + d[index] * delx));
  return;
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::eval_deriv
// Purpose       :
// Special Notes : adapted from GSL, version 1.15
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::eval_deriv (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & dydx)
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  ScalarT x_lo = xa[index];
  ScalarT delx = x - x_lo;
  dydx = b[index] + delx * (2.0 * c[index] + 3.0 * d[index] * delx);
  return;
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::eval_deriv2
// Purpose       :
// Special Notes : adapted from GSL, version 1.15
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::eval_deriv2 (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & ypp)
{
  size_t size = xa.size();
  size_t index = this->binarySearch (xa, x, 0, size - 1);

  const ScalarT x_lo = xa[index];
  const ScalarT delx = x - x_lo;
  ypp = 2.0 * c[index] + 6.0 * d[index] * delx;
  return;
}

//-----------------------------------------------------------------------------
// Function      : akima<ScalarT>::eval_integ
// Purpose       :
// Special Notes : adapted from GSL, version 1.15
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void akima<ScalarT>::eval_integ (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & ai,
   const ScalarT & bi,
   ScalarT & result)
{
  size_t size = xa.size();
  size_t index_a = this->binarySearch (xa, ai, 0, size - 1);
  size_t index_b = this->binarySearch (xa, bi, 0, size - 1);
  result = 0.0;

  // interior intervals
  for(size_t i=index_a; i<=index_b; i++)
  {
    const ScalarT x_hi = xa[i + 1];
    const ScalarT x_lo = xa[i];
    const ScalarT y_lo = ya[i];
    const ScalarT dx = x_hi - x_lo;
    if(dx != 0.0)
    {
      if (i == index_a || i == index_b)
      {
        ScalarT x1 = (i == index_a) ? ai : x_lo;
        ScalarT x2 = (i == index_b) ? bi : x_hi;
        result += this->integ_eval (y_lo, b[i], c[i], d[i], x_lo, x1, x2);
      }
      else
      {
        result += dx * (y_lo + dx*(0.5*b[i] + dx*(c[i]/3.0 + 0.25*d[i]*dx)));
      }
    }
    else
    {
      result = 0.0;
      return;
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Class         : cubic spline class
// Purpose       :
// Special Notes : adapted from Numerical Recipies.
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
class cubicSpline: public interpolator<ScalarT>
{
public:
  cubicSpline () {};

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya);

  void clear () { y2.clear(); };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y);

  void eval_deriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx);

  void eval_deriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp);

  // not implemented
  void eval_integ (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result) {};

public:
  std::vector<ScalarT> y2;
};

//-----------------------------------------------------------------------------
// Function      : cubicSpline<ScalarT>::init
// Purpose       :
// Special Notes : adapted from Numerical Recipies.  This roughly correponds
//                 to the "spline" function, which mainly creates the y2 array.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::init
(const std::vector<ScalarT> & xa,
 const std::vector<ScalarT> & ya)
{
  if (y2.size() != xa.size())
  {
    y2.resize(xa.size());
  }

  ScalarT p=0; ScalarT qn=0; ScalarT sig=0; ScalarT un=0;
  int n = y2.size(); std::vector <ScalarT> u(n-1,0.0);

  // Setting the upper and lower boundary conditions to a
  // "natural boundary condition".
  y2[0] = 0.0;
  y2[n-1] = 0.0;

  // This is the decomposition loop of the tridiagonal
  // algorithm.  y2 and u are used for temporary storage
  // of the decomposed factors.
  for (int i=1; i<n-1; i++)
  {
    sig = (xa[i]-xa[i-1])/(xa[i+1]-xa[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (ya[i+1]-ya[i])/(xa[i+1]-xa[i]) -
      (ya[i]-ya[i-1])/(xa[i]-xa[i-1]);
    u[i] = (6.0*u[i]/(xa[i+1]-xa[i-1]) - sig*u[i-1])/p;
  }

  for (int l=n-2; l>=0; l--)
  {
    y2[l] = y2[l]*y2[l+1]+u[l];
  }
};

//-----------------------------------------------------------------------------
// Function      : spline<ScalarT>::eval
// Purpose       :
// Special Notes : adapted from Numerical Recipies.  This roughly correponds
//                 to the "splint" function.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::eval(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & y_spline)
{
  // This method adapted from "Numerical Recipes in C++"
  int n = xa.size();
  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (xa[k] > x_position) khi=k;
    else klo=k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0)
  {
    // if out of range, then use the formula for dy/dx to extrapolate
    // beyond the range.  (formula 3.3.5 from numerical recipies in C)
    if (khi == 0)
    {
      ScalarT h0 = xa[1]-xa[0];
      ScalarT dx = x_position - xa[0];
      ScalarT dydx = (ya[1]-ya[0])/h0 - h0*y2[0]/3.0 - h0*y2[1]/6.0;
      y_spline = ya[0] + dx * dydx;
    }
    else if (klo == n-1)
    {
      ScalarT h1 = xa[n-1]-xa[n-2];
      ScalarT dx = x_position - xa[n-1];
      ScalarT dydx = (ya[n-1]-ya[n-2])/h1 + h1*y2[n-2]/6.0 + h1*y2[n-1]/3.0;
      y_spline = ya[n-1] + dx * dydx;
    }
  }
  else
  {
    a = (xa[khi] - x_position)/h;
    b = (x_position - xa[klo])/h;
    // cubic spline polynomial: (formula 3.3.3 from numerical recipies in C)
    y_spline = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2[klo] + (b*b*b-b)*y2[khi])*(h*h)/6.0;
  }
}

//-----------------------------------------------------------------------------
// Function      : cubicSpline<ScalarT>::eval_deriv
// Purpose       :
// Special Notes : adapted from Numerical Recipies.  This roughly correponds
//                 to the "splint" function.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::eval_deriv(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & dydx_spline)
{
  // This method adapted from "Numerical Recipes in C++"
  int n = xa.size();
  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (xa[k] > x_position) khi=k;
    else klo=k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0)
  {
    // if out of range, then use the formula for dy/dx to extrapolate
    // beyond the range.  (formula 3.3.5 from numerical recipies in C)
    if (khi == 0)
    {
      ScalarT h0 = xa[1]-xa[0];
      dydx_spline = (ya[1]-ya[0])/h0 - h0*y2[0]/3.0 - h0*y2[1]/6.0;
    }
    else if (klo == n-1)
    {
      ScalarT h1 = xa[n-1]-xa[n-2];
      dydx_spline = (ya[n-1]-ya[n-2])/h1 + h1*y2[n-2]/6.0 + h1*y2[n-1]/3.0;
    }
  }
  else
  {
    a = (xa[khi] - x_position)/h;
    b = (x_position - xa[klo])/h;

    // derivative:  (formula 3.3.5 from numerical recipies in C)
    dydx_spline =  (ya[khi]-ya[klo])/h - ((3*a*a-1)*y2[klo] - (3*b*b-1)*y2[khi])*h/6.0;
  }
}

//-----------------------------------------------------------------------------
// Function      : cubicSpline<ScalarT>::eval_deriv2
// Purpose       :
// Special Notes : adapted from Numerical Recipies.  This roughly correponds
//                 to the "splint" function.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void cubicSpline<ScalarT>::eval_deriv2(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x_position,
   ScalarT & ypp)
{
  // This method adapted from "Numerical Recipes in C++"
  int n = xa.size();
  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (xa[k] > x_position) khi=k;
    else klo=k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0)
  {
    // if out of range, assume no curvature.
    if (khi == 0)
    {
      ypp = 0.0;
    }
    else if (klo == n-1)
    {
      ypp = 0.0;
    }
  }
  else
  {
    a = (xa[khi] - x_position)/h;
    b = (x_position - xa[klo])/h;

    // derivative:  (formula 3.3.6 from numerical recipies in C)
    ypp =  a*y2[klo] + b*y2[khi];
  }
}

//-----------------------------------------------------------------------------
// Class         : linear interpolation class
// Purpose       :
// Special Notes : adapted from Numerical Recipies.
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
class linear: public interpolator<ScalarT>
{
public:
  linear () {};

  void init ( const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya) {};

  void clear () { };

  void eval (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & y);

  void eval_deriv (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & dydx);

  void eval_deriv2 (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & x, ScalarT & ypp);

  void eval_integ (
     const std::vector<ScalarT> & xa, const std::vector<ScalarT> & ya,
     const ScalarT & a, const ScalarT & b, ScalarT & result);

public:

};


//-----------------------------------------------------------------------------
// Function      : linear<ScalarT>::eval
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void linear<ScalarT>::eval(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & y)
{
  int n = xa.size();

  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h = xa[khi] - xa[klo];

  if (h == 0.0)
  {
    if (khi == 0)
    {
      y = xa[khi];
    }
    else if (klo == n-1)
    {
      y = xa[klo];
    }
  }
  else
  {
    ScalarT dx = x - xa[klo];
    ScalarT ya0 = ya[khi] - ya[klo];
    y = (dx/h) * ya0 + ya[klo];
  }
}

//-----------------------------------------------------------------------------
// Function      : linear<ScalarT>::eval_deriv
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void linear<ScalarT>::eval_deriv(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & dydx)
{
  int n = xa.size();

  // Find the right place in the table by means of bisection.
  ScalarT h = 0.0; ScalarT a = 0.0; ScalarT b = 0.0;
  int k = 0; int klo = 0; int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h = xa[khi] - xa[klo];

  if (h == 0.0)
  {
    if (khi == 0)
    {
      dydx = 0.0;
    }
    else if (klo == n-1)
    {
      dydx = 0.0;
    }
  }
  else
  {
    ScalarT dx = xa[khi] - xa[klo];
    ScalarT dy  = ya[khi] - ya[klo];
    dydx = dy/dx;
  }
}

//-----------------------------------------------------------------------------
// Function      : linear<ScalarT>::eval_deriv2
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
//-----------------------------------------------------------------------------
template <typename ScalarT>
void linear<ScalarT>::eval_deriv2(
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & x,
   ScalarT & ypp)
{
  ypp = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : linear<ScalarT>::eval_integ
// Purpose       :
// Special Notes : adapted from GSL, version 1.15
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/2012
// ----------------------------------------------------------------------------
template <typename ScalarT>
void linear<ScalarT>::eval_integ (
   const std::vector<ScalarT> & xa,
   const std::vector<ScalarT> & ya,
   const ScalarT & a,
   const ScalarT & b,
   ScalarT & result)
{

  int size = xa.size();
  int index_a = this->binarySearch (xa, a, 0, size - 1);
  int index_b = this->binarySearch (xa, b, 0, size - 1);

  // endpoints span more than one interval
  result = 0.0;

  // interior intervals
  for(int i=index_a; i<=index_b; i++)
  {
    const ScalarT x_hi = xa[i + 1];
    const ScalarT x_lo = xa[i];
    const ScalarT y_lo = ya[i];
    const ScalarT y_hi = ya[i + 1];
    const ScalarT dx = x_hi - x_lo;

    if(dx != 0.0)
    {
      if (i == index_a || i == index_b)
      {
        ScalarT x1 = (i == index_a) ? a : x_lo;
        ScalarT x2 = (i == index_b) ? b : x_hi;
        const ScalarT D = (y_hi-y_lo)/dx;
        result += (x2-x1) * (y_lo + 0.5*D*((x2-x_lo)+(x1-x_lo)));
      }
      else
      {
        result += 0.5 * dx * (y_lo + y_hi);
      }
    }
  }
  return;
}

} // namespace Device
} // namespace Xyce

#endif

