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
// Filename      : $RCSfile: N_ANP_SweepParam.h,v $
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 9/4/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIME_SWEEP_PARAM_H
#define Xyce_N_TIME_SWEEP_PARAM_H

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>

// ---------- Forward Declarations ----------

//-----------------------------------------------------------------------------
// Class         : N_ANP_SweepParam
//
// Purpose       : This class contains basic parameter data for parameter
//                 sweeps, for a single parameter.  If there are multiple
//                 parameters in the sweep, each one gets a class like
//                 this.
//
// Special Notes : "Step" here refers to steps in a parameter sweep loop,
//                 not time steps or DC sweep steps.
//
//
// Creator       : Eric Keiter, SNL
// Creation Date : 10/31/03
//-----------------------------------------------------------------------------
class N_ANP_SweepParam
{
public:
  // Default constructor.
  N_ANP_SweepParam () : 
   name(""),
   type("LIN"),
   startVal(0.0),
   stopVal(0.0),
   stepVal(0.0),
   stepMult(0.0),
   currentVal(0.0),
   numSteps(0),
   count(-1),
   maxStep(0),
   interval(1),
   outerStepNumber(0),
   sweepResetFlag_(false),
   lastLocalStepNumber_(-1)
   {};

  // Destructor
  ~N_ANP_SweepParam () {};

  bool updateCurrentVal (int stepNumber);
  bool getSweepResetFlag() {return sweepResetFlag_;};

  string name;
  string type;

  double startVal;
  double stopVal;
  double stepVal;
  double stepMult;

  double currentVal;

  int numSteps;

  int count;
  int maxStep;
  int interval;

  int outerStepNumber;

  vector <double> valList;

 private:
  bool sweepResetFlag_;
  int lastLocalStepNumber_;

};

//-----------------------------------------------------------------------------
// Function      : N_ANP_SweepParam::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/02/03
//-----------------------------------------------------------------------------
inline ostream & operator<<(ostream & os, const N_ANP_SweepParam & sp)
{
  os << "\tname            = " << sp.name;
  os << "\tcurrentVal      = " << sp.currentVal;
  os << endl;
  return os;
}

#endif // Xyce_N_TIME_SWEEP_PARAM_H

