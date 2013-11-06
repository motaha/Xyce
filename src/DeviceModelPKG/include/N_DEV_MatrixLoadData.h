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
// Filename       : $RCSfile: N_DEV_MatrixLoadData.h,v $
//
// Purpose        :
//
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/30/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_MatrixLoadData_h
#define Xyce_N_DEV_MatrixLoadData_h

#include <vector>

// ---------- Xyce Includes  ----------
#include <N_UTL_Misc.h>

namespace Xyce {
namespace Device {

class colData;
class valData;

//-----------------------------------------------------------------------------
// Class         : MatrixLoadData
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
class MatrixLoadData
{
public:
  MatrixLoadData ();

  MatrixLoadData (const MatrixLoadData & right);

  ~MatrixLoadData ();

  bool initializeAll (int isizeTmp = 100);

  void resizeTestJacSolData(int size);
  void resizeTestJacQData(int size);
  void resizeTestJacStateData(int size);

public:
  int isize;
  int isizeNumJac;

  // temporary jacobian load structures:
  vector<int>    cols;
  vector<double> vals;
  vector<double> Qvals;

  // temporary numerical jacobian load structures:
  vector<valData> val_local;
  vector<valData> Qval_local;
  vector<colData> col_local;
  vector<int>     row_local;
  vector<int>     internalFlag;

  // Structures used by the "testJacobian" function.
  vector <vector <double> > numJac;
  vector < vector<double> > saveJac;
  vector < vector<double> > devJac;
  vector < vector<double> > diffJac;
  vector < vector<double> > relJac;

  vector <vector <double> > numJacQ;
  vector < vector<double> > saveJacQ;
  vector < vector<double> > devJacQ;
  vector < vector<double> > diffJacQ;
  vector < vector<double> > relJacQ;

  vector <vector <int> > status;
  vector <vector <int> > stencil;
  vector <vector <int> > statusQ;

  vector <double> saveRHS;
  vector <double> pertRHS;
  vector <double> origRHS;
  vector <double> saveQ;
  vector <double> pertQ;
  vector <double> origQ;

  vector <double> saveSoln;
  vector <double> pertSoln;
  vector <double> saveCurrSoln;

  vector <double> saveLastState;
  vector <double> saveCurrState;
  vector <double> saveNextState;
  vector <double> saveStateDerivs;
};

//-----------------------------------------------------------------------------
// Class         : colData
// Purpose       : This class contains a vector of column indices.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class colData
{
public:
  colData (int isizeTmp = 100):
    isize(isizeTmp), col()
  { col.reserve(isizeTmp); }

public:
  int isize;
  vector<int> col;
};

//-----------------------------------------------------------------------------
// Class         : valData
// Purpose       : This class contains a vector of value indices.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/30/02
//-----------------------------------------------------------------------------
class valData
{
private:
protected:
public:
  valData (int isizeTmp = 100):
    isize(isizeTmp), val() { val.reserve(isizeTmp); }

private:
protected:
public:
  int isize;
  vector<double> val;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MatrixLoadData N_DEV_MatrixLoadData;

#endif

