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
// Filename      : $RCSfile: N_ANP_SweepParam.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 9/4/04
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.6.4.2 $
// Revision Date  : $Date: 2013/10/03 17:23:31 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <iostream>

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif


// ----------   Xyce Includes   ----------
#include <N_ANP_SweepParam.h>
#include <N_ERH_ErrorMgr.h>

// ----------- Forward declarations ------------

//-----------------------------------------------------------------------------
// Function      : N_ANP_SweepParam::updateCurrentVal
//
// Purpose       : Updates the values of the parameters used in a sweep.
//
// Special Notes : This is very similar to the "update" function in the
//                 class N_DEV_SweepData.  (which no longer exists).
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 10/31/2003
//-----------------------------------------------------------------------------
bool N_ANP_SweepParam::updateCurrentVal (int stepNumberArg)
{
  outerStepNumber  = stepNumberArg/interval;
  int inum            = outerStepNumber/maxStep;
  int localStepNumber = outerStepNumber - inum*maxStep;

  // We must keep track of whether we're at the first step of our sweep,
  // because some device manager features need to know that.
  // It is important that we only set this when localStepNumber first becoms
  // zero, not every time localStepNumber *is* zero, because an outer loop
  // might remain at 0 for quite some time.

  if (localStepNumber == 0 && localStepNumber != lastLocalStepNumber_)
  {
    sweepResetFlag_=true;
  } 
  else
  {
    sweepResetFlag_=false;
  }
  lastLocalStepNumber_=localStepNumber;

  if (type == "LIN")
  {
    currentVal = startVal + static_cast<double>(localStepNumber)*stepVal;
    ++count;
  }
  else if (type == "DEC" || type == "OCT")
  {
    currentVal = startVal*pow(stepMult, static_cast<double>(localStepNumber) );
    ++count;
  }
  else if (type == "LIST")
  {
    int size=  valList.size();
    int index = (localStepNumber < size)?localStepNumber:(size-1);
    currentVal = valList[index];
    ++count;
  }
  else
  {
    string msg = "N_ANP_SweepParam::updateCurrentVal: ";
    msg += " Unsupported type specified.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

#ifdef Xyce_DEBUG_TIME
  cout << endl;
  cout << "----------------------" << endl;
  cout << "updateCurrentVal" << endl;
  cout << "  name             = " << name << endl;
  cout << "  stepNumberArg    = " << stepNumberArg<< endl;
  cout << "  interval         = " << interval  << endl;
  cout << "  outerStepNumber  = " << outerStepNumber << endl;
  cout << "  localStepNumber  = " << localStepNumber << endl;
  cout << "  inum             = " << inum      << endl;
  cout << "  sweepResetFlag   = " << sweepResetFlag_ << endl;
  cout << "  currentVal       = " << currentVal << endl;
  cout << "----------------------" << endl;
#endif

  return true;
}

