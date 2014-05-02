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
// Filename       : $RCSfile: N_IO_Objective.h,v $
//
// Purpose        : Objective class for Dakota runs
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 03/05/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.24.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Objective_h
#define Xyce_N_IO_Objective_h

#include <list>
#include <string>
#include <vector>
#include <set>
#include <time.h>

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_UTL_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_fwd.h>
#include <N_ANP_fwd.h>
#include <N_UTL_fwd.h>

#ifdef Xyce_Dakota
// ---------- Trilinos Includes ----------
#include <Teuchos_RefCountPtrDecl.hpp>
#include <Epetra_SerialDenseVector.h>
#endif

class N_LAS_Vector;

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : Objective
//
// Purpose       : This class manages a single objective (for Dakota).
//
// Special Notes :
//
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/13/06
//-----------------------------------------------------------------------------
class Objective
{
public:
  Objective();
  ~Objective();
  void readData();
  void printData();

  void reset ();
  double save (const N_LAS_Vector *, const N_LAS_Vector *, const N_LAS_Vector *);
  double save (double, double, const N_LAS_Vector *, const N_LAS_Vector *, const N_LAS_Vector *);
  double evaluate ();

  bool parse(const std::string & in, OutputMgr * outputMgrPtr );
  bool initialize(const Util::OptionBlock & , OutputMgr * outputMgrPtr );

  std::string var1;
  std::string var2;
  std::string file;
  std::string function;                      // the function that will be applied to value={}
  Util::ExpressionData *value;
  Util::ExpressionData *weight;

private:
  bool initializeInternal (std::string & fileStr, std::string & functionStr,
                           std::string & valueStr, std::string & weightStr);

  OutputMgr * outputMgrPtr_;

  double lastResult;                    // last value calculated during simulation
  double minResult;                     // minimum value encountered during simulation
  double maxResult;                     // maximum value encountered during simulation
  double magResult;                     // magnitude of value encountered during simulation
  double rmsResult;                     // root-mean-squared value calculated during simulation
  double lastV1;
  double lastWeight;
  int lastInd2;
  int lastIndInterpolate;
  int n1;                               // number of elements in first independent var array
  int n2;                               // number of elements in second independent var array
  int lastN1, lastN2;                   // last time we saved data, this was where it was inserted.
  // this can speed up the search for insertaiton points.
  bool lastResultValid;
  std::vector<double> var1Vals;              // first independent variable
  std::vector<double> var2Vals;              // second independent variable (if needed)
  std::vector<std::vector<double> > dataVal;      // externally supplied data[var2Index][var1Index]
  std::vector<std::vector<double> > simVal;       // simulation calculated data[var2Index][var1Index]
  std::vector<std::vector<bool> > simValValid;    // flags indicating if simulation was done at point [var1Index][var2Index]
  std::vector<std::vector<double> > weightVal;    // weighting filter function
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::Objective N_IO_Objective;

#endif // Xyce_Objective_h
