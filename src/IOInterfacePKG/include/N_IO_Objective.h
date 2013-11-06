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
// Revision Number: $Revision: 1.16.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Objective_h
#define Xyce_N_IO_Objective_h

// ---------- Standard Includes ----------
#include <list>
#include <string>
#include <vector>
#include <set>
#include <time.h>

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_UTL_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_fwd.h>

#ifdef Xyce_Dakota
// ---------- Trilinos Includes ----------
#include <Teuchos_RefCountPtrDecl.hpp>
#include <Epetra_SerialDenseVector.h>
#endif

// ---------- Forward Declarations ----------
class N_TOP_Topology;
class N_UTL_ExpressionData;
class N_ANP_AnalysisInterface;
class N_LAS_Vector;

//-----------------------------------------------------------------------------
// Class         : N_IO_Objective
//
// Purpose       : This class manages a single objective (for Dakota).
//
// Special Notes :
//
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/13/06
//-----------------------------------------------------------------------------
class N_IO_Objective
{
  public:
    N_IO_Objective();
    ~N_IO_Objective();
    void readData();
    void printData();

    void reset ();
    double save (const N_LAS_Vector *, const N_LAS_Vector *, const N_LAS_Vector *);
    double save (double, double, const N_LAS_Vector *, const N_LAS_Vector *, const N_LAS_Vector *);
    double evaluate ();

    bool parse(const string & in, N_IO_OutputMgr * outputMgrPtr );
    bool initialize(const N_UTL_OptionBlock & , N_IO_OutputMgr * outputMgrPtr );
/*
    bool parse(const string & in, N_DEV_DeviceInterface * devPtr1, N_ANP_AnalysisInterface *,
#ifdef Xyce_PARALLEL_MPI
	  N_PDS_Comm * pdsCommPtr1,
#endif
	  N_TOP_Topology * topPtr1);

    bool initialize(const N_UTL_OptionBlock & , N_DEV_DeviceInterface *, N_ANP_AnalysisInterface *,
#ifdef Xyce_PARALLEL_MPI
          N_PDS_Comm *,
#endif
          N_TOP_Topology *);
*/

    string var1;
    string var2;
    string file;
    string function;                      // the function that will be applied to value={}
//    map<string, string> match;
    N_UTL_ExpressionData *value;
    N_UTL_ExpressionData *weight;

  private:
    bool initializeInternal (string & fileStr, string & functionStr,
		    string & valueStr, string & weightStr);

   // N_DEV_DeviceInterface * devPtr_;
   // N_ANP_AnalysisInterface * anaIntPtr_;
//#ifdef Xyce_PARALLEL_MPI
//    N_PDS_Comm * pdsCommPtr_;
//#endif
//    N_TOP_Topology * topPtr_;
    N_IO_OutputMgr * outputMgrPtr_;

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
    vector<double> var1Vals;              // first independent variable
    vector<double> var2Vals;              // second independent variable (if needed)
    vector<vector<double> > dataVal;      // externally supplied data[var2Index][var1Index]
    vector<vector<double> > simVal;       // simulation calculated data[var2Index][var1Index]
    vector<vector<bool> > simValValid;    // flags indicating if simulation was done at point [var1Index][var2Index]
    vector<vector<double> > weightVal;    // weighting filter function
};

#endif // Xyce_N_IO_Objective_h
