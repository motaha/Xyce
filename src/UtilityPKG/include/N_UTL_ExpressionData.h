//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Class         : N_UTL_ExpressionData
//
// Purpose       : This class manages a single expression.
//
// Special Notes : This class owns a single expression class pointer,
//                 and all the auxilliary data needed to manage its
//                 usage.
//
//                 I felt that this would be a good way to avoid
//                 bloat in functions like outputPRINT.
//
//                 I made all the data public.  Right now, the
//                 evaluate function just returns the value of the
//                 expression.
//
//                 The derivative values are all in
//                 the valDerivs vector after each evaluation.  If
//                 it becomes necessary later to access derivatives
//                 for some reason, a developer only has to access
//                 the public valDerivs vector, right after an
//                 evaluate function call.
//
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------

#ifndef N_UTL_ExpressionData_H
#define N_UTL_ExpressionData_H

#include <string>
#include <list>
#include <vector>

#include <N_IO_fwd.h>

class N_LAS_Vector;
class N_ANP_AnalysisInterface;
class N_UTL_Param;

#ifdef Xyce_PARALLEL_MPI
class N_PDS_Comm;
#endif /* Xyce_PARALLEL_MPI */

class N_TOP_Topology;
class N_UTL_Expression;

#include <N_UTL_Xyce.h>

class N_UTL_ExpressionData
{
  public:
    N_UTL_ExpressionData (const string & expression1,
        N_IO_OutputMgr * outputMgrPtr
	    );

	N_UTL_ExpressionData (const N_UTL_Expression * expressionPointer,
        N_IO_OutputMgr * outputMgrPtrIn
    );

    ~N_UTL_ExpressionData ();

    int setup ();
    double evaluate (const N_LAS_Vector *solVecPtr, const N_LAS_Vector *stateVecPtr, const N_LAS_Vector * stoVecPtr);
    int getNumUnresolvedStrings() const { return numUnresolvedStrings; };

  public:
    N_UTL_Expression * expPtr;
    N_IO_OutputMgr * outputMgrPtr_;

    string         expression;

    int            numVars;
    int            numSpecialVars;       // non solution vars like TIME
    int            numLeads;
    int            numUnresolvedStrings; // number of unresolved strings after setup() is called.
    vector<string> varNames;
    vector<string> specialVarNames;      // non solution vars like TIME
    vector<string> currentParam;
    vector<int>    varGIDs;
    list<N_UTL_Param> expressionVars_;  // all vars from expression for evaluation by OutputMgr.
    vector<double> varVals;
    double         val;

    bool setupCalled;
    bool numUnresolvedStringsChecked;
  private:
    int procID_;
};

#endif  // N_UTL_ExpressionData_H
