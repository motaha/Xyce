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
// Filename       : $RCSfile: N_DEV_DeviceSupport.C,v $
//
// Purpose        : This file contains similar functions to the spice3f5
//                  file, devsup.c.  It contains support routines for
//                  device models
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/17/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.43 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#if( defined HAVE_FLOAT_H && defined HAVE__ISNAN_AND__FINITE_SUPPORT )
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#else
using std::isnan;
using std::isinf;
#endif


// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_DEV_DeviceSupport.h>
#include <N_DEV_Const.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::lambertw
// Purpose       : provides a lambert-w function for diodes and BJT's.
// Special Notes :
//
// Purpose.  Evaluate principal branch of Lambert W function at x.
//
// w = w(x) is the value of Lambert's function.
// ierr = 0 indicates a safe return.
// ierr = 1 if x is not in the domain.
// ierr = 2 if the computer arithmetic contains a bug.
// xi may be disregarded (it is the error).
//
// Prototype: void lambertw( double, double, int, double);
//
// Reference:
// T.C. Banwell
// Bipolar transistor circuit analysis using the Lambert W-function,
// IEEE Transactions on Circuits and Systems I: Fundamental Theory
// and Applications
//
// vol. 47, pp. 1621-1633, Nov. 2000.
//
// Scope         : public
// Creator       : David Day,  SNL
// Creation Date : 04/16/02
//-----------------------------------------------------------------------------
void DeviceSupport::lambertw(double x, double &w, int &ierr, double &xi)
{
  int i=0, maxit = 10;
  const double turnpt = -exp(-1.), c1 = 1.5, c2 = .75;
  double r, r2, r3, s, mach_eps, relerr = 1., diff;
  mach_eps = 2.e-15;   // float:2e-7
  ierr = 0;

  if( x > c1)
  {
    w = c2*log(x);
    xi = log( x/ w) - w;
  }
  else
  {
    if( x >= 0.0)
    {
      w = x;
      if( x == 0. ) return;
      if( x < (1-c2) ) w = x*(1.-x + c1*x*x);
      xi = - w;
    }
    else
    {
      if( x >= turnpt)
      {
        if( x > -0.2 )
        {
          w = x*(1.0-x + c1*x*x);
          xi = log(1.0-x + c1*x*x) - w;
        }
        else
        {
          diff = x-turnpt;
          if( diff < 0.0 ) diff = -diff;
          w = -1 + sqrt(2.0*exp(1.))*sqrt(x-turnpt);
          if( diff == 0.0 ) return;
          xi = log( x/ w) - w;
        }
      }
      else
      {
        ierr = 1; // x is not in the domain.
        w = -1.0;
        return;
      }
    }
  }

  while( relerr > mach_eps  && i<maxit)
  {
     r = xi/(w+1.0);   //singularity at w=-1
     r2 = r*r;
     r3 = r2*r;
     s  = 6.*(w+1.0)*(w+1.0);
     w = w * (  1.0 + r + r2/(2.0*( w+1.0)) - (2. * w -1.0)*r3/s  );
     if( w * x < 0.0 ) w = -w;
     xi = log( x/ w) - w;

     if( x>1.0 )
     {
       relerr =  xi / w;
     }
     else
     {
       relerr =  xi;
     }
     if(relerr < 0.0 ) relerr = -relerr;
     ++i;
   }
   if( i == maxit ) ierr = 2;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::limvds
// Purpose       : limit the per-iteration change of VDS
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
double DeviceSupport::limvds ( double vnew, double vold)
{

  if(vold >= 3.5)
  {
    if(vnew > vold) vnew = Xycemin(vnew,(3.0 * vold) + 2.0);
    else
    {
      if (vnew < 3.5) vnew = Xycemax(vnew,2.0);
    }
  }
  else
  {
    if(vnew > vold) vnew = Xycemin(vnew, 4.0);
    else            vnew = Xycemax(vnew,-0.5);
  }
  return(vnew);
}


//-----------------------------------------------------------------------------
// Function      : DeviceSupport::pnjlim
// Purpose       : limit the per-iteration change of PN junction voltages
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
double DeviceSupport::pnjlim (
    double vnew,
    double vold,
    double vt,
    double vcrit,
    int *icheck
  )
{
  double arg;

  if((vnew > vcrit) && (fabs(vnew - vold) > (vt + vt)))
  {
    if(vold > 0)
    {
      arg = 1 + (vnew - vold) / vt;

      if(arg > 0)
      {
        vnew = vold + vt * log(arg);
      }
      else
      {
        vnew = vcrit;
      }
    }
    else
    {
      vnew = vt *log(vnew/vt);
    }

    *icheck = 1;
  }
  else
  {
    *icheck = 0;
  }

  return(vnew);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::pnjlim_new
// Purpose       : limit the per-iteration change of PN junction voltages
// Special Notes : Copied from NGSpice, which has the following comment:
//     This code has been fixed by Alan Gillespie adding limiting
//     for negative voltages
// Scope         : public
// Creator       : Tom Russo, SNL 1445, Electrical Systems Modeling
// Creation Date : 11/19/2012
//-----------------------------------------------------------------------------
double DeviceSupport::pnjlim_new (
    double vnew,
    double vold,
    double vt,
    double vcrit,
    int *icheck
  )
{
  double arg;

  if((vnew > vcrit) && (fabs(vnew - vold) > (vt + vt)))
  {
    if(vold > 0)
    {
      arg = (vnew - vold) / vt;

      if(arg > 0)
      {
        vnew = vold + vt * (2+log(arg-2));
      }
      else
      {
        vnew = vold - vt*(2+log(2-arg));
      }
    }
    else
    {
      vnew = vt *log(vnew/vt);
    }

    *icheck = 1;
  }
  else
  {
    if (vnew < 0)
    {
      if (vold > 0)
      {
        arg= -vold -1;
      }
      else
      {
        arg = 2*vold -1;
      }
      if (vnew < arg)
      {
        vnew = arg;
        *icheck=1;
      }
      else
      {
        *icheck = 0;
      }
    }
    else
    {

      *icheck = 0;
    }
  }

  return(vnew);
}


//-----------------------------------------------------------------------------
// Function      : DeviceSupport::fetlim
// Purpose       : limit the per-iteration change of FET voltages
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
double DeviceSupport::fetlim (
    double vnew,
    double vold,
    double vto
  )
{
  double vtsthi;
  double vtstlo;
  double vtox;
  double delv;
  double vtemp;

#ifdef Xyce_FETLIM_NAN_CHECK
  // Make sure that vto is not a Nan.
  if (finiteNumberTest(vto) < 0)
  {
    vto = 1.0; // kludge.
  }
#endif

  vtsthi = fabs(2*(vold-vto))+2;
  vtstlo = vtsthi/2 +2;
  vtox = vto + 3.5;
  delv = vnew-vold;

  if (vold >= vto)
  {
    if(vold >= vtox)
    {
      if(delv <= 0)
      {
        // going off
        if(vnew >= vtox)
        {
          if(-delv >vtstlo) vnew =  vold - vtstlo;
        }
        else
        {
          vnew = Xycemax(vnew,vto+2.0);
        }
      }
      else
      {
        // staying on
        if(delv >= vtsthi) vnew = vold + vtsthi;
      }
    }
    else
    {
      // middle region
      if(delv <= 0)
      {
        // decreasing
        vnew = Xycemax(vnew,vto-0.5);
      }
      else
      {
        // increasing
        vnew = Xycemin(vnew,vto+4.0);
      }
    }
  }
  else
  {
    // off
    if(delv <= 0)
    {
      if(-delv >vtsthi) vnew = vold - vtsthi;
    }
    else
    {
      vtemp = vto + 0.5;
      if(vnew <= vtemp)
      {
        if(delv >vtstlo) vnew = vold + vtstlo;
      }
      else
      {
        vnew = vtemp;
      }
    }
  }
  return(vnew);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::cmeyer
// Purpose       : Compute the MOS overlap capacitances as functions of the
//                 device terminal voltages
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------

void DeviceSupport::cmeyer (
   double vgs0,    // initial voltage gate-source
   double vgd0,    // initial voltage gate-drain
   double vgb0,    // initial voltage gate-bulk
   double von0,
   double vdsat0,
   double vgs1,    // final voltage gate-source
   double vgd1,    // final voltage gate-drain
   double vgb1,    // final voltage gate-bulk
   double covlgs,  // overlap capacitance gate-source
   double covlgd,  // overlap capacitance gate-drain
   double covlgb,  // overlap capacitance gate-bulk
   register double *cgs,
   register double *cgd,
   register double *cgb,
   double phi,
   double cox,
   double von,
   double vdsat
 )
{
    double vdb;
    double vdbsat;
    double vddif;
    double vddif1;
    double vddif2;
    double vgbt;

    *cgs = 0;
    *cgd = 0;
    *cgb = 0;

    vgbt = vgs1-von;
    if (vgbt <= -phi)
    {
      *cgb = cox;
    }
    else if (vgbt <= -phi/2)
    {
      *cgb = -vgbt*cox/phi;
    }
    else if (vgbt <= 0)
    {
      *cgb = -vgbt*cox/phi;
      *cgs = cox/(7.5e-1*phi)*vgbt+cox/1.5;
    }
    else
    {
      vdbsat = vdsat-(vgs1-vgb1);
      vdb = vgb1-vgd1;
      if (vdbsat <= vdb)
      {
        *cgs = cox/1.5;
      }
      else
      {
        vddif = 2.0*vdbsat-vdb;
        vddif1 = vdbsat-vdb-1.0e-12;
        vddif2 = vddif*vddif;
        *cgd = cox*(1.0-vdbsat*vdbsat/vddif2)/1.5;
        *cgs = cox*(1.0-vddif1*vddif1/vddif2)/1.5;
      }
    }

    vgbt = vgs0-von0;
    if (vgbt <= -phi)
    {
      *cgb += cox;
    }
    else if (vgbt <= -phi/2)
    {
      *cgb += -vgbt*cox/phi;
    }
    else if (vgbt <= 0)
    {
      *cgb += -vgbt*cox/phi;
      *cgs += cox/(7.5e-1*phi)*vgbt+cox/1.5;
    }
    else
    {
      vdbsat = vdsat0-(vgs0-vgb0);
      vdb = vgb0-vgd0;
      if (vdbsat <= vdb)
      {
        *cgs += cox/1.5;
      }
      else
      {
        vddif = 2.0*vdbsat-vdb;
        vddif1 = vdbsat-vdb-1.0e-12;
        vddif2 = vddif*vddif;
        *cgd += cox*(1.0-vdbsat*vdbsat/vddif2)/1.5;
        *cgs += cox*(1.0-vddif1*vddif1/vddif2)/1.5;
      }
    }

    *cgs = *cgs *.5 + covlgs;
    *cgd = *cgd *.5 + covlgd;
    *cgb = *cgb *.5 + covlgb;
}


//-----------------------------------------------------------------------------
// Function      : DeviceSupport::qmeyer
// Purpose       : Compute the MOS overlap capacitances as functions of the
//                 device terminal voltages
//
// Special Notes : ARGSUSED  because vgb is no longer used
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
void DeviceSupport::qmeyer (
   double vgs,    // initial voltage gate-source
   double vgd,    // initial voltage gate-drain
   double vgb,    // initial voltage gate-bulk
   double von,
   double vdsat,
   double & capgs,  // non-constant portion of g-s overlap capacitance
   double & capgd,  // non-constant portion of g-d overlap capacitance
   double & capgb,  // non-constant portion of g-b overlap capacitance
   double phi,
   double cox     // oxide capactiance
 )
{
  double vds;
  double vddif;
  double vddif1;
  double vddif2;
  double vgst;

  //double vgdt;
  //double vdenom;
  //double vdenom2;


  vgst = vgs-von;
  if (vgst <= -phi)
  {
    capgb = cox/2;
    capgs = 0;
    capgd = 0;
  }
  else if (vgst <= -phi/2)
  {
    capgb = -vgst*cox/(2*phi);
    capgs = 0;
    capgd = 0;
  }
  else if (vgst <= 0)
  {
    capgb = -vgst*cox/(2*phi);
    capgs = vgst*cox/(1.5*phi)+cox/3;
    capgd = 0;
  }
  else
  {
    vds = vgs-vgd;
    if (vdsat <= vds)
    {
      capgs = cox/3;
      capgd = 0;
      capgb = 0;
    }
    else
    {

      vddif = 2.0*vdsat-vds;
      vddif1 = vdsat-vds/*-1.0e-12*/;
      vddif2 = vddif*vddif;
      capgd = cox*(1.0-vdsat*vdsat/vddif2)/3;
      capgs = cox*(1.0-vddif1*vddif1/vddif2)/3;
      capgb = 0;


      //vgdt = vgd-von;
      //vdenom=vgst + vgdt;
      //vdenom2=vdenom*vdenom;

      //capgd = cox*(1.0-vgst*vgst/vdenom2)/3.0;
      //capgs = cox*(1.0-vgdt*vgdt/vdenom2)/3.0;

    }
  }

}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::qmeyerderivs
// Purpose       : Computes the partial derivatives of the Meyer capacitances
//                 with respect to various voltages.
//
// Special Notes : We need this in order to make the Meyer model MPDE
//                 compatible.
//
// Scope         : public
// Creator       : Keith R. Santarelli, SNL, Electrical & Microsystems Modeling
// Creation Date : 02/08/08
//-----------------------------------------------------------------------------

void DeviceSupport::qmeyerderivs
    (
       double vgs,    // initial voltage gate-source
       double vgd,    // initial voltage gate-drain
       double vgb,    // initial voltage gate-bulk
       double von,
       double vdsat,
       double & dcapgsdvgs, //partial deriv. of capgs with respect to vgs
       double & dcapgsdvgb, //partial deriv. of capgs with respect to vgb
       double & dcapgsdvgd, //partial deriv. of capgs with respect to vgd
       double & dcapgddvgs, //partial deriv. of capgd with respect to vgs
       double & dcapgddvgb, //partial deriv. of capgd with respect to vgb
       double & dcapgddvgd, //partial deriv. of capgd with respect to vgd
       double & dcapgbdvgs, //partial deriv. of capgb with respect to vgs
       double & dcapgbdvgb, //partial deriv. of capgb with respect to vgb
       double & dcapgbdvgd, //partial deriv. of capgb with respect to vgd
       double phi,
       double cox,     // oxide capactiance
       int Dtype //transistor type
    )
{
  double vgst;
  double vds;
  double vdenom, vdenom3;
  double vgdt;

  vgst = vgs-von;
  if (vgst <= -phi)
  {
    dcapgsdvgs=0;
    dcapgsdvgb=0;
    dcapgsdvgd=0;
    dcapgddvgs=0;
    dcapgddvgb=0;
    dcapgddvgd=0;
    dcapgbdvgs=0;
    dcapgbdvgb=0;
    dcapgbdvgd=0;

  }
  else if (vgst <= -phi/2)
  {
    dcapgsdvgs=0;
    dcapgsdvgb=0;
    dcapgsdvgd=0;
    dcapgddvgs=0;
    dcapgddvgb=0;
    dcapgddvgd=0;
    dcapgbdvgs=-1.0*cox/(2.0*phi);
    dcapgbdvgb=0;
    dcapgbdvgd=0;
  }
  else if (vgst <= 0)
  {
    dcapgsdvgs=cox/(1.5*phi);
    dcapgsdvgb=0;
    dcapgsdvgd=0;
    dcapgddvgs=0;
    dcapgddvgb=0;
    dcapgddvgd=0;
    dcapgbdvgs=-1.0*cox/(2.0*phi);
    dcapgbdvgb=0;
    dcapgbdvgd=0;
  }
  else
  {
    vds = vgs-vgd;
    if (vdsat <= vds)
    {
      dcapgsdvgs=0;
      dcapgsdvgb=0;
      dcapgsdvgd=0;
      dcapgddvgs=0;
      dcapgddvgb=0;
      dcapgddvgd=0;
      dcapgbdvgs=0;
      dcapgbdvgb=0;
      dcapgbdvgd=0;
    }
    else
    {
      vgdt = vgd-von;
      vdenom=vgst + vgdt;
      vdenom3=vdenom*vdenom*vdenom;

      dcapgsdvgs=4.0/3.0*cox*vgdt*vgdt/vdenom3;
      dcapgsdvgb=0;
      dcapgsdvgd=-4.0/3.0*cox*vgst*vgdt/vdenom3;
      dcapgddvgs=-4.0/3.0*cox*vgst*vgdt/vdenom3;
      dcapgddvgb=0;
      dcapgddvgd=4.0/3.0*cox*vgst*vgst/vdenom3;
      dcapgbdvgs=0;
      dcapgbdvgb=0;
      dcapgbdvgd=0;

    }
  }

  //Now have to "type-ize" the cap. derivatives:

  //dcapgsdvgs=Dtype*dcapgsdvgs;
  //dcapgsdvgb=Dtype*dcapgsdvgb;
  //dcapgsdvgd=Dtype*dcapgsdvgd;
  //dcapgddvgs=Dtype*dcapgddvgs;
  //dcapgddvgb=Dtype*dcapgddvgb;
  //dcapgddvgd=Dtype*dcapgddvgd;
  //dcapgbdvgs=Dtype*dcapgbdvgs;
  //dcapgbdvgb=Dtype*dcapgbdvgb;
  //dcapgbdvgd=Dtype*dcapgbdvgd;
}



#ifdef notdef

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::cap
// Purpose       : compute equivalent conductances
//                 divide up the channel charge (1-xqc)/xqc to
//                 source and drain
//
// Special Notes : XXX This is no longer used, apparently
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
void DeviceSupport::cap (
  register CKTcircuit *ckt,
  double vgd,
  double vgs,
  double vgb,
  double covlgd,
  double covlgs,
  double covlgb,
  double capbd,
  double capbs,
  double cggb,
  double cgdb,
  double cgsb,
  double cbgb,
  double cbdb,
  double cbsb,
  double *gcggb,
  double *gcgdb,
  double *gcgsb,
  double *gcbgb,
  double *gcbdb,
  double *gcbsb,
  double *gcdgb,
  double *gcddb,
  double *gcdsb,
  double *gcsgb,
  double *gcsdb,
  double *gcssb,
  double qgate,
  double qchan,
  double qbulk,
  double *qdrn,
  double *qsrc,
  double xqc)
{
  double gcd;
  double gcdxd;
  double gcdxs;
  double gcg;
  double gcgxd;
  double gcgxs;
  double gcs;
  double gcsxd;
  double gcsxs;
  double qgb;
  double qgd;
  double qgs;

  gcg = (cggb+cbgb)*ckt->CKTag[1];
  gcd = (cgdb+cbdb)*ckt->CKTag[1];
  gcs = (cgsb+cbsb)*ckt->CKTag[1];
  gcgxd = -xqc*gcg;
  gcgxs = -(1-xqc)*gcg;
  gcdxd = -xqc*gcd;
  gcdxs = -(1-xqc)*gcd;
  gcsxd = -xqc*gcs;
  gcsxs = -(1-xqc)*gcs;
  *gcdgb = gcgxd-covlgd*ckt->CKTag[1];
  *gcddb = gcdxd+(capbd+covlgd)*ckt->CKTag[1];
  *gcdsb = gcsxd;
  *gcsgb = gcgxs-covlgs*ckt->CKTag[1];
  *gcsdb = gcdxs;
  *gcssb = gcsxs+(capbs+covlgs)*ckt->CKTag[1];
  *gcggb = (cggb+covlgd+covlgs+covlgb)*ckt->CKTag[1];
  *gcgdb = (cgdb-covlgd)*ckt->CKTag[1];
  *gcgsb = (cgsb-covlgs)*ckt->CKTag[1];
  *gcbgb = (cbgb-covlgb)*ckt->CKTag[1];
  *gcbdb = (cbdb-capbd)*ckt->CKTag[1];
  *gcbsb = (cbsb-capbs)*ckt->CKTag[1];
  /*
   *     compute total terminal charges
   */
  qgd = covlgd*vgd;
  qgs = covlgs*vgs;
  qgb = covlgb*vgb;
  qgate = qgate+qgd+qgs+qgb;
  qbulk = qbulk-qgb;
  *qdrn = xqc*qchan-qgd;
  *qsrc = (1-xqc)*qchan-qgs;
  /*
   *     finished
   */
}
#endif

#if 0
double DeviceSupport::pred
  (
    CKTcircuit *ckt,
    int loct
  )
{

  /* predict a value for the capacitor at loct by
   * extrapolating from previous values
   */

#ifndef NEWTRUNC
  double xfact;
  xfact = ckt->CKTdelta/ckt->CKTdeltaOld[1];
  return( ( (1+xfact) * *(ckt->CKTstate1+loct) ) -
          (    xfact  * *(ckt->CKTstate2+loct) )  );
#endif /*NEWTRUNC*/

}
#endif

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::contVds
// Purpose       : continuation adjustment for MOSFET drain-source voltage.
//
// Special Notes : min*vds is the value returned for vds when alpha=0.
//
//                 This idea is based, loosely, on a paper by Jaijeet
//                 Rosychowdhury.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/04/03
//-----------------------------------------------------------------------------
double DeviceSupport::contVds (double vds, double alpha, double min)
{
  if (min <= 0.0) min = 0.3;
  return ( vds * (alpha*(1.0-min) + min) );
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::contVgst
// Purpose       : continuation adjustment for MOSFET drain-source voltage.
//
// Special Notes : vgstConst is the value returned for vgst when alpha=0.
//
//                 This idea is based, loosely, on a paper by Jaijeet
//                 Rosychowdhury.
//
//                 The alpha=0 condition, for which vgs, or vgst is
//                 considered constant, essentially makes the device a
//                 single-state device.   No matter what actual voltage is
//                 applied across Vg-Vs, the device will act as though
//                 there is a fixed applied voltage.  If vgstconst is set
//                 to a high voltage (like the default, 3.0) then the
//                 devices are all in an "on" state for alpha=0.  If
//                 vgstconst is set to zero, then the devices are all in an
//                 off state for alpha=0.
//
//                 As the alpha parameter is swept from zero to 1, at some
//                 point the MOSFET will start behaving, functionally like
//                 a MOSFET.  At alpha=0, it doesn't - it kind of acts like
//                 a capacitor.  At what point in the sweep this
//                 functionality kicks in, depends on the point at which
//                 the Ids vs Vgs curve has an intersection with Ids=0, for
//                 reasonable values of Vgs.  (reasonable values probably
//                 being -1 to +5 volts)
//
//                 If Vgstconst is set to zero, this intersection happens
//                 immediately, on the first step where alpha is nonzero.
//                 If Vgstconst is fairly high (say 4.0 volts), this
//                 intersection happens late.  If it is too high, it never
//                 happens.
//
//                 Unfortunately, during the sweep, once it kicks in, it
//                 kicks in for every MOSFET almost all at once, and there
//                 is a sharp transition in the continuation curve.  The
//                 trick to making this homotopy work well, is to make this
//                 transition as gentle as possible, I think.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/04/03
//-----------------------------------------------------------------------------
double DeviceSupport::contVgst
  (double vgst, double alpha, double vgstConst)
{
  return ((alpha)*vgst + (1.0-alpha)*vgstConst);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::getGainScaleBlockID
// Purpose       : determines which block to associate the mosfet device
//                 This is used for a multi-block gainscale continuation.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 01/28/05
//-----------------------------------------------------------------------------
int DeviceSupport::getGainScaleBlockID(int numBlocks)
{
  double val = u.RandomDouble();

  // get rid of negatives
  val = val * val;
  val = sqrt(val);

  double interval = 1.0 / numBlocks;

  for (int i = 0; i < numBlocks; ++i) {
    double lowBound = ((double) i) * interval;
    double highBound = ((double) (i+1)) * interval;
    if ((val >= lowBound) && (val < highBound)) {
      return i;
    }
  }

  return (numBlocks-1);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::getGainScaleBlockID
// Purpose       : computes random perturbation for stretched homotopy.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 01/28/05
//-----------------------------------------------------------------------------
double DeviceSupport::getRandomPerturbation()
{
  double val = u.RandomDouble();
  val = val * val;
  val = sqrt(val);
  return val;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::getGainScaleBlockID
// Purpose       : sets seed for random number generator used in getRandomPerturbation().
// Scope         : public
// Creator       : Richard Schie, Electrical Systems Modeling
// Creation Date : 10/01/12
//-----------------------------------------------------------------------------
int DeviceSupport::SetSeed(unsigned int seedIn)
{
  int retVal = u.SetSeed( seedIn );
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::Xexp
//
// Purpose       : Exponential function with range checking and truncation
//                 for improved convergence when searching for DCOP
//
// Special Notes : This function provides a Taylor series approximation
//                 to an exponential (out to N=order terms),
//                 unless order >=20, in which case it just returns an
//                 exponential.
//
// Scope         : public
//
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/13/06
//-----------------------------------------------------------------------------
double DeviceSupport::Xexp(double X, double &ddX, double order)
{
  if (X > CONSTMAX_EXP_ARG)
  {
    X = CONSTMAX_EXP_ARG;
  }
#if 1
  if (order > 19.999)
  {
    ddX=exp(X);
    return ddX;
  }
  int i, ord;
  double ans, term;

  double base;

  // It is very important to return and use the correct derivative, otherwise
  // conductances are not correct.
  {
    term = 1;
    ans = 1;
    ddX = 1;
    ord = static_cast<int>( order );
    for (i=1 ; i<=ord ; ++i)
    {
      term *= X/i;
      ans += term;
      if (i<ord)
      {
        ddX += term;
      }
    }
    // If we are partway between two integer orders "ord", linear combination
    // of the order=ord and order=ord+1 series.
    if (order > ord)
    {
      ddX += (order-ord)*term;
      term *= X/(ord+1);
      ans += (order-ord)*term;
    }
  }
  return ans;
#else

  // An attempt at a lambertw replacement of the exponential.
  // At the moment it is not as effective a homotopy as the taylor series.
  // For NAND chains longer than 276 it starts to have big trouble around
  // order 9, getting NaN's in the RHS.
  if (order > 19.999)
  {
    ddX=exp(X);
    return ddX;
  }

  double rs,wreturn,xi;
  int i,ierr;
  // rs=1e-{order}
  rs=pow(10.0,-order);
  lambertw(rs*exp(X+rs),wreturn,ierr,xi);
  if (ierr != 0)
  {
    ddX = exp(X);
    return (ddX);
  }
  else
  {
    ddX = wreturn/(1+wreturn)/rs;
    return (wreturn/rs);
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::finiteNumberTest
//
// Purpose       : The motivation for this function is that as of
//                 this writing, it is possible for some of the MOSFET
//                 models to pass in Nan's for the "von" parameter.
//
// Special Notes : This is copied over from a similar
//                 function in the nonlinear solver.
//
// Scope         : public
// Creator       : Eric Keiter, SNL.
// Creation Date : 01/05/07
//-----------------------------------------------------------------------------
int DeviceSupport::finiteNumberTest(const double x)
{
#if( defined HAVE_NAN_INF_SUPPORT || defined HAVE__ISNAN_AND__FINITE_SUPPORT )
  if (isnan(x))
    return -1;

  if (isinf(x))
    return -2;

#else

  // Can pretty much use any number here
  double tol = 1.0e-6;

  // NaN check
  if (!(x <= tol) && !(x > tol))
    return -1;

  // Inf check:
  // Use the fact that: Inf * 0 = NaN
  double z = 0.0 * x;
  if (!(z <= tol) && !(z > tol))
    return -2;

#endif

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::spline
// Purpose       :
// Special Notes : This method adapted from "Numerical Recipes in C++"
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/12/06
//-----------------------------------------------------------------------------
void DeviceSupport::spline(
      std::vector<double> & x,
      std::vector<double> & y,
      std::vector<double> & y2)
{
  double p=0;
  double qn=0;
  double sig=0;
  double un=0;

  int n = y2.size();
  std::vector<double> u(n-1,0.0);

  // Setting the upper and lower boundary conditions to a
  // "natural boundary condition".
  y2[0] = 0.0;
  y2[n-1] = 0.0;

  // This is the decomposition loop of the tridiagonal
  // algorithm.  y2 and u are used for temporary storage
  // of the decomposed factors.
  for (int i=1; i<n-1; i++)
  {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) -
            (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
  }

  for (int l=n-2; l>=0; l--)
  {
    y2[l] = y2[l]*y2[l+1]+u[l];
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::splint
// Purpose       :
// Special Notes : This method adapted from "Numerical Recipes in C++"
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/12/06
//-----------------------------------------------------------------------------
void DeviceSupport::splint(
    std::vector<double> & xa,
    std::vector<double> & ya,
    std::vector<double> & y2a,
    double x_position,
    double & y_spline)
{
  // This method taken from "Numerical Recipes in C++"
  int n = xa.size();
  // Find the right place in the table by means of bisection.
  double h = 0.0;
  double a = 0.0;
  double b = 0.0;
  int k = 0;
  int klo = 0;
  int khi = n-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (xa[k] > x_position) khi=k;
    else klo=k;
  }

  h = xa[khi] - xa[klo];
  if (h == 0.0)
  {
    //std::string err_msg =
    //"Bad input to cubic spline.";
      //throw charon::Exception("This shouldn't happen! Please notify "
                                //"developers!", __FILE__, __LINE__);

    exit (0);
  }
  a = (xa[khi] - x_position)/h;
  b = (x_position - xa[klo])/h;
  // cubic spline polynomial is now evaluated
  y_spline = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi])*(h*h)/6.0;


  // capping the spline value so that it doesn't extend
  // the bounds of the endpoints
  // this prevents wiggles and non-physical values
  if ( (y_spline > ya[klo]) && (y_spline > ya[khi]) )
  {
    if (ya[klo] > ya[khi]) y_spline = ya[klo];
    else y_spline = ya[khi];
  }

  if ( (y_spline < ya[klo]) && (y_spline < ya[khi]) )
  {
    if (ya[klo] < ya[khi]) y_spline = ya[klo];
    else y_spline = ya[khi];
  }

}

} // namespace Device
} // namespace Xyce
