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
// Filename       : $RCSfile: N_DEV_DopeInfo.h,v $
//
// Purpose        : This file contains the PDE device instance base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DopeInfo_h
#define Xyce_N_DEV_DopeInfo_h

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Const.h>
#include <N_DEV_CompositeParam.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_DopeInfo
// Purpose       : Doping region information.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/28/03
//-----------------------------------------------------------------------------
class DopeInfo : public CompositeParam
{
  friend class ParametricData<DopeInfo>;

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

public:
  static ParametricData<DopeInfo> &getParametricData();

  DopeInfo ();
  bool processParam (Param & ndParam, string & param, DevicePDEInstance & di);
  void processParams ();

  void setupInfo (
    vector<double> & CVec,
    vector<double> & CdonorVec,
    vector<double> & CacceptorVec,
    vector<double>  & xVec, DeviceSupport & devSup);

  void setupInfo2d (
    vector<double> & CVec,
    vector<double> & CdonorVec,
    vector<double> & CacceptorVec,
    vector<double> & xVec, vector<double> &yVec, DeviceSupport & devSup);

  static double nsdep(double x, double W, double Dt);
  static double ngdep(double x, double y, double W, double ax, double ay);
  static double ngdep2(double x, double y, double ax, double ay);
  static double erf(double x);

  static void readDopingFile (string & filename, vector<double> & xloc,
                              vector<double> & nvec, vector<double> & y2, DeviceSupport & devSup);

  static void readDopingFile (string & filename, vector<double> & xloc,
                              vector<double> & nvec, vector<double> & y2_n,
                              vector<double> & pvec, vector<double> & y2_p, DeviceSupport & devSup);

public:
  string name;         // this is also the index into the map.
  string type;         // p-type or n-type
  string funcType;     // uniform, step, gaussian.
  string speciesName;  // boron, phosphorus, etc.
  string fileName;     // doping file.
  string exprString;   // if expression based, store here.

  // location params:
  double xmin;
  double xmax;
  bool xminGiven;
  bool xmaxGiven;

  double xloc;
  double xwidth;

  double ymin;
  double ymax;
  bool yminGiven;
  bool ymaxGiven;

  double yloc;
  double ywidth;

  double Nmax;     // maximum magnitude
  double Nmin;     // minimum magnitude

  // sometimes having a really high doping hurts convergence.
  // This lets the user set truncate it.
  double Nmax_chop;
  bool Nmax_chopGiven;

  int flatX;
  int flatY;

  // arrays for spline fitting, if specified from file:
  vector<double> xlocVec;
  vector<double> dopeVec;
  vector<double> y2Vec;
  vector<double> splintDopeVec;
};

// inline functions
//-----------------------------------------------------------------------------
// Function      : DopeInfo::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, SNL, Parallel Computational Sciences
// Creation Date : 03/31/03
//-----------------------------------------------------------------------------
inline ostream & operator<<(ostream & os, const DopeInfo & di)
{
  os << di.name << ":\n";
  os << "  type     = " << di.type << "\n";
  os << "  funcType = " << di.funcType << "\n";

  if (di.funcType == "expression")
  {
    os << "  exp. string = " << di.exprString << "\n";
  }

  os << "  xloc     = " << di.xloc     << "\n";
  os << "  xwidth   = " << di.xwidth   << "\n";
  os << "  yloc     = " << di.yloc     << "\n";
  os << "  ywidth   = " << di.ywidth   << "\n";

  os << "  Nmax     = " << di.Nmax     << "\n";
  os << "  Nmin     = " << di.Nmin     << "\n";

  os << "  flatX    = " << di.flatX << "\n";
  os << "  flatY    = " << di.flatY << "\n";

  os << endl;

  return os;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DopeInfo N_DEV_DopeInfo;

#endif

