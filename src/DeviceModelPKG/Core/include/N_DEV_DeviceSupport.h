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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DeviceSupport.h,v $
//
// Purpose        : This file contains the device support class.  Most of the
//                 functions of this class are similar to those of devsup.c
//                 in 3f5.
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
// Revision Number: $Revision: 1.21.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Device_Support_h
#define Xyce_N_DEV_Device_Support_h

#include <Epetra_Util.h>

#include <vector>

namespace Xyce {
namespace Device {

class DeviceSupport
{
public:

  void lambertw (double x, double &w, int &ierr, double &xi);

  double limvds ( double vnew, double vold);

  double pnjlim ( double vnew, double vold, double vt, double vcrit,
                  int *icheck);

  double pnjlim_new ( double vnew, double vold, double vt, double vcrit,
                      int *icheck);

  double fetlim ( double vnew, double vold, double vto);

  void cmeyer
  (
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
   );

  void qmeyer
  (
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
   );


  //KRS, 2/8/08:  Adding this new function to compute the partial
  //derivatives of the Meyer capacitances as a function of various
  //voltages.

  void qmeyerderivs
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
   );


#ifdef notdef
  void cap (
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
    double xqc);

#endif

  double contVds (double vds, double alpha, double min = 0.3);

  double contVgst (double vgst, double alpha, double vgstConst = 3.0);

  int getGainScaleBlockID(int numBlocks); // For homotopy

  double getRandomPerturbation();  // For homotopy
  int SetSeed(unsigned int seedIn); // to set the random seed for getRandomPerturbation().

  double Xexp(double, double &,double);

  int finiteNumberTest(const double x);

  void spline( vector<double> & x, vector<double> & y, vector<double> & y2);
  void splint ( vector<double> & xa, vector<double> & ya, vector<double> & y2a,
                double x_position, double & y_spline);

#if 0
  double pred ( CKTcircuit *ckt, int loct);
#endif

protected:
  // For block gainscale homotopy
  Epetra_Util u;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceSupport N_DEV_DeviceSupport;

#endif

