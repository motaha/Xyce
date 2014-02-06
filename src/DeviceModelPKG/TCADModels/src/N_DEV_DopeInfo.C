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
// Filename       : $RCSfile: N_DEV_DopeInfo.C,v $
//
// Purpose        : This file contains the details of the dope info class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/04/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.17.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_DEV_DevicePDEInstance.h>
#include <N_UTL_Expression.h>
#include <N_DEV_DeviceSupport.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<DopeInfo>::ParametricData()
{
  addPar("NMAX", 1.0e+15, false, &DopeInfo::Nmax);
  addPar("NMIN", 0.0, false, &DopeInfo::Nmin);
  addPar("NMAXCHOP", 1.0e+20, false, &DopeInfo::Nmax_chop, &DopeInfo::Nmax_chopGiven);
  addPar("XLOC", 0.0, false, &DopeInfo::xloc);
  addPar("XMIN", 0.0, false, &DopeInfo::xmin, &DopeInfo::xminGiven);
  addPar("XMAX", 0.0, false, &DopeInfo::xmax, &DopeInfo::xmaxGiven);
  addPar("XWIDTH", 1.0e-3, false, &DopeInfo::xwidth);
  addPar("YLOC", 0.0, false, &DopeInfo::yloc);
  addPar("YMIN", 0.0, false, &DopeInfo::ymin, &DopeInfo::yminGiven);
  addPar("YMAX", 0.0, false, &DopeInfo::ymax, &DopeInfo::ymaxGiven);
  addPar("YWIDTH", 1.0e-3, false, &DopeInfo::ywidth);

  // Set up map for non-double precision variables:
  addPar("NAME", std::string("none"), false, &DopeInfo::name);
  addPar("FUNCTION", std::string("uniform"), false, &DopeInfo::funcType);
  addPar("TYPE", std::string("ntype"), false, &DopeInfo::type);
  addPar("FLATX", 0, false, &DopeInfo::flatX);
  addPar("FLATY", 0, false, &DopeInfo::flatY);
  addPar("SPECIES", std::string("none"), false, &DopeInfo::speciesName);
  addPar("FILE", std::string("none"), false, &DopeInfo::fileName);
  addPar("EXPRESSION", std::string("none"), false, &DopeInfo::exprString);
}

ParametricData<DopeInfo> &DopeInfo::getParametricData() {
  static ParametricData<DopeInfo> parMap;

  return parMap;
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::DopeInfo
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/07/05
// ----------------------------------------------------------------------------
DopeInfo::DopeInfo ()
  : CompositeParam (),
    name("reg0"),
    type("ntype"),
    funcType("uniform"),
    speciesName("none"),
    fileName("none"),

    xmin(0.0),
    xmax(0.0),
    xloc(0.0),
    xwidth(0.0),

    ymin(0.0),
    ymax(0.0),
    yloc(0.0),
    ywidth(0.0),

    Nmax(1.0e+15),
    Nmin(1.0e+11),
    Nmax_chop(1.0e+99),
    Nmax_chopGiven(false),
    flatX(0),
    flatY(0)
{

  // // Set up mapping from param names to class variables:
  // if (parMap.empty())
  // {
  //   // Set up map for double precision variables:
  //   addPar ("NMAX", 1.0e+15, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::Nmax),
  //           NULL);

  //   addPar ("NMIN", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::Nmin),
  //           NULL);

  //   addPar ("NMAXCHOP", 1.0e+20, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::Nmax_chop),
  //           static_cast <bool CompositeParam:: *> (&DopeInfo::Nmax_chopGiven));

  //   addPar ("XLOC", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::xloc),
  //           NULL);

  //   addPar ("XMIN", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::xmin),
  //           static_cast <bool CompositeParam:: *> (&DopeInfo::xminGiven));

  //   addPar ("XMAX", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::xmax),
  //           static_cast <bool CompositeParam:: *> (&DopeInfo::xmaxGiven));

  //   addPar ("XWIDTH", 1.0e-3, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::xwidth),
  //           NULL);

  //   addPar ("YLOC", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::yloc),
  //           NULL);

  //   addPar ("YMIN", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::ymin),
  //           static_cast <bool CompositeParam:: *> (&DopeInfo::yminGiven));

  //   addPar ("YMAX", 0.0, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::ymax),
  //           static_cast <bool CompositeParam:: *> (&DopeInfo::ymaxGiven));

  //   addPar ("YWIDTH", 1.0e-3, false,
  //           static_cast <double CompositeParam:: *> (&DopeInfo::ywidth),
  //           NULL);

  //   // Set up map for non-double precision variables:
  //   addPar ("NAME", "none", false,
  //           static_cast <string CompositeParam:: *> (&DopeInfo::name),
  //           NULL);

  //   addPar ("FUNCTION", "uniform", false,
  //           static_cast <string CompositeParam:: *> (&DopeInfo::funcType),
  //           NULL);

  //   addPar ("TYPE", "ntype", false,
  //           static_cast <string CompositeParam:: *> (&DopeInfo::type),
  //           NULL);

  //   addPar ("FLATX", 0, false,
  //           static_cast <int CompositeParam:: *> (&DopeInfo::flatX),
  //           NULL);

  //   addPar ("FLATY", 0, false,
  //           static_cast <int CompositeParam:: *> (&DopeInfo::flatY),
  //           NULL);

  //   addPar ("SPECIES", "none", false,
  //           static_cast <string CompositeParam:: *> (&DopeInfo::speciesName),
  //           NULL);

  //   addPar ("FILE", "none", false,
  //           static_cast <string CompositeParam:: *> (&DopeInfo::fileName),
  //           NULL);

  //   addPar ("EXPRESSION", "none", false,
  //           static_cast <string CompositeParam:: *> (&DopeInfo::exprString),
  //           NULL);


  // }
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::processParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/31/03
// ----------------------------------------------------------------------------
bool DopeInfo::processParam
(Param & ndParam, string & param, DevicePDEInstance & di)
{
  bool bsuccess = true;

  return bsuccess;
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/31/03
// ----------------------------------------------------------------------------
void DopeInfo::processParams()
{
  {
    ParameterMap::const_iterator p_i = (*getPMap()).find(string("FUNCTION"));
    const Pars &p = static_cast<const Pars &>(*(*p_i).second);

    ExtendedString tmp = getValue<std::string, DopeInfo>(*this, p);
    setValue<std::string, DopeInfo>(*this, p, static_cast<std::string>(tmp.toLower()));
  }

  {
    ParameterMap::const_iterator p_i = (*getPMap()).find(string("TYPE"));
    const Pars &p = static_cast<const Pars &>(*(*p_i).second);

    ExtendedString tmp = getValue<std::string, DopeInfo>(*this, p);
    setValue<std::string, DopeInfo>(*this, p, static_cast<std::string>(tmp.toLower()));
  }
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::setupInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/20/08
// ----------------------------------------------------------------------------
void DopeInfo:: setupInfo
(vector<double> & CVec,
 vector<double> & CdonorVec,
 vector<double> & CacceptorVec,
 vector<double> & xVec,
 DeviceSupport & devSupport)
{
  int i(0);
  int NX (CVec.size());
  splintDopeVec.resize(NX,0.0);

  double sign = 0.0;
  if (type == "ptype" || type == "acceptor")
  {
    sign = -1.0;
  }
  else if (type == "ntype" || type == "donor")
  {
    sign = 1.0;
  }

  if (funcType == "uniform")
  {
    for (i=0;i<NX;++i)
    {
      if (xmaxGiven && xminGiven)
      {
        if (xVec[i] > xmax || xVec[i] < xmin) // if outside the range, skip
        {
          continue;
        }
      }

      CVec[i] += sign*Nmax;
      splintDopeVec[i] = Nmax;
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax;
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += Nmax;
      }
    }
  }
  else if (funcType == "gaussian")
  {
    double deltaX = fabs(xwidth);

    double ax = log(Nmax/Nmin)/(deltaX*deltaX);
    for (i=0;i<NX;++i)
    {
      double scalarX = 1.0;
      double abs_dx;

      if (given("XLOC") && given("XWIDTH") && deltaX != 0.0)
      {
        abs_dx = fabs(xVec[i]-xloc);

        // if true gaussian, x-section:
        if (flatX==0)
        {
          scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
        }
        else if (flatX>0) // half-gaussian, x-section:
        {
          bool flatReg = (xVec[i] > xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
        else if (flatX<0) // half-gaussian, x-section:
        {
          bool flatReg = (xVec[i] < xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
      }

      CVec[i] += sign*Nmax*scalarX;
      splintDopeVec[i] = Nmax*scalarX;
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax*scalarX;
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += Nmax*scalarX;
      }
    }
  }
  else if (funcType == "step")
  {
    for (i=0;i<NX;++i)
    {
      double x  = xVec[i];
      bool regOn = true;

      if (given("XLOC"))
      {
        if (flatX ==  0) regOn = true;

        if (flatX == -1)
        {
          if (x > xloc) regOn = false;
          else          regOn = true;
        }

        if (flatX == +1)
        {
          if (x < xloc) regOn = false;
          else          regOn = true;
        }
      }

      CVec[i] += (regOn)?(sign*Nmax):(sign*Nmin);
      splintDopeVec[i] = (regOn)?(Nmax):(Nmin);
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
    }
  }
  else if (funcType == "expression")
  {
    if (exprString == "none")
    {
      string msg = "Dope Region : ";
      msg += name;
      msg += " has specified the expression specification, but not provided an expression.";
      msg += "\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    else
    {
#ifdef Xyce_DEBUG_DEVICE
      cout << "DopeInfo::setupInfo: exprString = " << exprString << endl;
#endif
      N_UTL_Expression expr;
      expr.set(exprString);

      for (i=0;i<NX;++i)
      {
        double dopeValue(0.0);

        expr.set_var(string("#X"), xVec[i]);
        expr.evaluateFunction (dopeValue);
        CVec[i] += sign*dopeValue;
        splintDopeVec[i] = dopeValue;
        if (type == "ptype" || type == "acceptor")
        {
          CacceptorVec[i] += dopeValue;
        }
        else if (type == "ntype" || type == "donor")
        {
          CdonorVec[i] += dopeValue;
        }
      }
    }
  }
  else if (funcType == "file")
  {
    if (fileName == "none")
    {
      string msg = "Dope Region : ";
      msg += name;
      msg += " has specified the file specification, but not specified a file name.";
      msg += "\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    else
    {
      readDopingFile (fileName, xlocVec, dopeVec, y2Vec, devSupport);

      // if the user has requested that this be truncated to a max value,
      // do it here:
      if (Nmax_chopGiven)
      {
        int dopeSize=dopeVec.size();
        for (int id=0;id<dopeSize;++id)
        {
          if (dopeVec[id] > Nmax_chop)
          {
            dopeVec[id] = Nmax_chop;
          }
        }
      }

      for (i=0;i<NX;++i)
      {
        double xtmp = xVec[i];
        double dopeValue(0.0);
        devSupport.splint(xlocVec, dopeVec, y2Vec, xtmp, dopeValue);
        CVec[i] += sign*dopeValue;
        splintDopeVec[i] = dopeValue;
        if (type == "ptype" || type == "acceptor")
        {
          CacceptorVec[i] += dopeValue;
        }
        else if (type == "ntype" || type == "donor")
        {
          CdonorVec[i] += dopeValue;
        }
      }
    }
  }
  else
  {
    string msg = "Unrecognized Dope Region function type:  ";
    msg += funcType;
    msg += "  for region: ";
    msg += name;
    msg += "\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::setupInfo2d
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/20/08
// ----------------------------------------------------------------------------
void DopeInfo:: setupInfo2d
(vector<double> & CVec,
 vector<double> & CdonorVec,
 vector<double> & CacceptorVec,
 vector<double> & xVec,
 vector<double> & yVec,
 DeviceSupport & devSupport)
{
  int i(0);
  int numMeshPoints (CVec.size());

  double sign = 1.0;
  if (type == "ptype" || type == "acceptor") sign = -1.0;

  if (funcType == "uniform")
  {
    for (i=0;i<numMeshPoints;++i)
    {

      if (xmaxGiven && xminGiven)
      {
        if (xVec[i] > xmax || xVec[i] < xmin) // if outside the x-range, skip
        {
          continue;
        }
      }

      if (ymaxGiven && yminGiven)
      {
        if (yVec[i] > ymax || yVec[i] < ymin) // if outside the y-range, skip
        {
          continue;
        }
      }

      CVec[i] += sign*Nmax;
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax;
      }
      else
      {
        CdonorVec[i] += Nmax;
      }
    }
  }
  else if (funcType == "gaussian")
  {
    double deltaX = fabs(xwidth);
    double deltaY = fabs(ywidth);

    double ax = 0.0;
    double ay = 0.0;

    if (deltaX!=0.0)
    {
      ax = log(Nmax/Nmin)/(deltaX*deltaX);
    }

    if (deltaY!=0.0)
    {
      ay = log(Nmax/Nmin)/(deltaY*deltaY);
    }

    for (i=0;i<numMeshPoints;++i)
    {
      double scalarX = 1.0;
      double scalarY = 1.0;
      double abs_dx, abs_dy;

      if (given("XLOC") && given("XWIDTH") && deltaX != 0.0)
      {
        abs_dx = fabs(xVec[i]-xloc);

        // if true gaussian, x-section:
        if (flatX==0)
        {
          scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
        }
        else if (flatX>0) // half-guassian, x-section:
        {
          bool flatReg = (xVec[i] > xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
        else if (flatX<0) // half-guassian, x-section:
        {
          bool flatReg = (xVec[i] < xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
      }

      if (given("YLOC") && given("YWIDTH") && deltaY != 0.0)
      {
        abs_dy = fabs(yVec[i]-yloc);
        // if true gaussian, y-section:
        if (flatY==0)
        {
          scalarY *= ngdep2(0.0, abs_dy, 1.0, ay);
        }
        else if (flatY>0) // half-guassian, y-section:
        {
          bool flatReg = (yVec[i] > yloc);

          if (!flatReg)
          {
            scalarY *= ngdep2(0.0, abs_dy, 1.0, ay);
          }
        }
        else if (flatY<0) // half-guassian, y-section:
        {
          bool flatReg = (yVec[i] < yloc);

          if (!flatReg)
          {
            scalarY *= ngdep2(0.0, abs_dy, 1.0, ay);
          }
        }
      }

      CVec[i] += sign*Nmax*scalarX*scalarY;

      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax*scalarX*scalarY;
      }
      else
      {
        CdonorVec[i] += Nmax*scalarX*scalarY;
      }
    }
  }
  else if (funcType == "step")
  {
    for (i=0;i<numMeshPoints;++i)
    {
      double x = xVec[i];
      double y = yVec[i];
      bool regOnX = true;
      bool regOnY = true;

      if (given("YLOC"))
      {
        if (flatY ==  0) regOnY = true;

        if (flatY == -1)
        {
          if(y > yloc) regOnY = false;
          else            regOnY = true;
        }

        if (flatY == +1)
        {
          if (y < yloc) regOnY = false;
          else             regOnY = true;
        }
      }

      if (given("XLOC"))
      {
        if (flatX ==  0) regOnX = true;

        if (flatX == -1)
        {
          if(x > xloc) regOnX = false;
          else            regOnX = true;
        }

        if (flatX == +1)
        {
          if (x < xloc) regOnX = false;
          else             regOnX = true;
        }
      }
      bool regOn = (regOnX && regOnY);

      CVec[i] += (regOn)?(sign*Nmax):(sign*Nmin);
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
      else
      {
        CdonorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
    }
  }
  else
  {
    string msg = "Unrecognized Dope Region function type:  ";
    msg += funcType;
    msg += "  for region: ";
    msg += name;
    msg += "\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::nsdep
// Purpose       : This function returns an approximate deposition profile
//                 of a step implant driven in an inert environment.
// Special Notes :
//                        1       W/2 + x          W/2 - x
//        nsdep(x,W,Dt) = - (erf(---------) + erf(----------))
//                        2      2*sqrt(Dt)       2*sqrt(Dt)
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DopeInfo::nsdep (double x, double W, double Dt)
{
  double D  = 2.0 * sqrt(Dt);
  double Wh = W / 2.0;
  return 0.5 * (erf((Wh + x)/D) + erf((Wh - x)/D));
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::ngdep
//
// Purpose       : This function returns an approximate Gaussian deposition.
//
//                 Adapted from a similar function in SGF.
//
// Special Notes : This function is a little flakey, in that it assumes
//                 that we're using a cylindrical geomtry
//                 (x = radius, y = height.) It also assumes that the (0,0)
//                 origin is in the upper-left-hand corner of the mesh.
//
//                 So, y=0.0 is the upper surface of the device, and the
//                 implant is coming from the top (above y=0.0).  Hence,
//                 the (y>0) conditional.
//
//                 Also, the width parameter (W), is set up to be a
//                 diameter about x=0, which is the reason for the 0.5*W -
//                 half of this diameter will impact this radius.
//
//                 The implant is completely flat and constant in the
//                 x-direction, as long as fabs(x) is less than W/2.0.
//                 Beyond W/2.0, the gaussian profile kicks in.
//
//                 The parameters ax and ay are scaling parameters, and
//                 correspond to how much you want the doping to vary with
//                 space.  A typical value for either can be set up as:
//
//                 Ax = ln (Nhi / Nlo)/(Rx*Rx)
//
//                 where:
//
//                   Nhi = the max. level of doping
//                   Nlo = the min. level of doping
//                   Rx  = distance over which this doping should vary.
//
//                   Nhi/Nlo = 10^N, where N = # of orders of magnitude
//                   that should vary between x=0 and x=Rx.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DopeInfo::ngdep
(double x, double y, double W, double ax, double ay)
{
  double xprime = fabs(x) - (0.5 * W);
  return ((xprime <= 0.0) ? 1.0 : exp(-ax*xprime*xprime))*
    ((y      >  0.0) ? 0.0 : exp(-ay*y*y));
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::ngdep2
//
// Purpose       : This function returns an approximate Gaussian deposition.
//
// Special Notes : This function is a modification of the original ngdep
//                 (see above), and I designed it to address some of the
//                 peculiarities of ngdep.
//
//                 (1) I've gotten rid of the width, W.  I'm just
//                 going to assume that whoever is calling this function
//                 can set that(the constant region) up on their own.
//
//                 (2) I've removed the conditionals cause things to
//                 be set to zero, or one, or whatever, if you are on one
//                 side or another of the suface.  I'm assuming that
//                 whoever calls this function can do that themselves, if
//                 they need to.
//
//                 (3) I've removed the stuff that sets the retVal to zero
//                 for y>0.  Again, this is the user's problem.
//
//                 ax and ay mean the same as they did for the original
//                 ngdep.  (see above).
//
//                 It is possible to use this for the 1D case, pretty
//                 easily.  Set the xflag to false, hold y fixed at zero,
//                 and have x correspond to the 1D mesh locations. (or, set
//                 xflag to true, hold x fixed at zero, and let y
//                 correspond to 1D mesh locations - either way).
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/28/03
// ----------------------------------------------------------------------------
double DopeInfo::ngdep2
(double x, double y, double ax, double ay)
{
  double retVal = exp(-ax*x*x)* exp(-ay*y*y);
  return retVal;
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::erf
// Purpose       : This function returns the error function.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DopeInfo::erf(double x)
{
  double t1 = 1.0 / (1.0 + 0.3275911 * fabs(x));
  double t2 = t1 * t1;
  double t3 = t2 * t1;
  double t4 = t3 * t1;
  double t5 = t4 * t1;
  double result = 1.0 - (0.254829592*t1 - 0.284496736*t2 + 1.421413741*t3 -
                         1.453152027*t4 + 1.061405429*t5) * exp(-x*x);
  return (x < 0.0) ? -result : result;
}

//-----------------------------------------------------------------------------
// Function      : DopeInfo::readDopingFile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/07
//-----------------------------------------------------------------------------
void DopeInfo::readDopingFile
(string & filename, vector<double> & xloc, vector<double> & nvec,
 vector<double> & y2, DeviceSupport & devSupport)
{
  ifstream input;
  double x_loc(0.0);
  double value(0.0);
  xloc.clear();
  nvec.clear();
  y2.clear();

  input.open( filename.c_str(), ios::in );
  if ( input.good() )
  {
    bool endOfFile = input.eof();
    while (!endOfFile)
    {
      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> x_loc;
      }
      else
      {
        break;
      }

      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> value;
      }
      else
      {
        break;
      }
      xloc.push_back(x_loc);
      nvec.push_back(value);
    }
    input.close();
    y2.resize(xloc.size(),0.0);
    devSupport.spline (xloc, nvec, y2);
  }
  else
  {
    string msg = "Error:  Cannot open doping file: " + filename;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : DopeInfo::readDopingFile
// Purpose       :
// Special Notes : This version assumes 2 dopants are in the file, P and N.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/07
//-----------------------------------------------------------------------------
void DopeInfo::readDopingFile (string & filename,
                               vector<double> & xloc, vector<double> & nvec, vector<double> & y2_n,
                               vector<double> & pvec, vector<double> & y2_p,
                               DeviceSupport & devSupport)
{
  ifstream input;
  double x_loc(0.0);
  double value1(0.0);
  double value2(0.0);
  xloc.clear();
  nvec.clear();
  pvec.clear();
  y2_n.clear();
  y2_p.clear();

  input.open( filename.c_str(), ios::in );
  if ( input.good() )
  {
    bool endOfFile = input.eof();
    while (!endOfFile)
    {
      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> x_loc;
      }
      else
      {
        break;
      }

      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> value1;
      }
      else
      {
        break;
      }
      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> value2;
      }
      else
      {
        break;
      }

      xloc.push_back(x_loc);
      nvec.push_back(value1);
      pvec.push_back(value2);

      //cout << "x="<<x_loc<<"  value1="<<value1<<"  value2="<<value2<<endl;
    }
    input.close();
    y2_n.resize(xloc.size(),0.0);
    y2_p.resize(xloc.size(),0.0);
    devSupport.spline (xloc, nvec, y2_n);
    devSupport.spline (xloc, pvec, y2_p);
  }
  else
  {
    string msg = "Error:  Cannot open doping file: " + filename;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }
}

} // namespace Device
} // namespace Xyce
