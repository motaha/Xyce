//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_MembraneUserDefined.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Christy Warrender, Cognitive Modeling
//
// Creation Date  : 12/14/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MembraneUserDefined_h
#define Xyce_N_DEV_MembraneUserDefined_h

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_MembraneModel.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MembraneUserDefined
// Purpose       : This is class defines a user-defined membrane mechanism.
// Special Notes :
// Creator       : Christy Warrender, Cognitive Modeling
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
class MembraneUserDefined : public MembraneModel
{
  public:
    MembraneUserDefined(SolverState & ss1, double cMem, double gMem, double vRest,
    	vector<string> & currentEqus, vector<string> & indepVars, vector<string> & fEqs,
    	vector<string> & qEqs, vector<string> & extraFunctions, vector<string> & extraParameters);
    ~MembraneUserDefined() {}

    void updateSecondaryState( double * staDerivVec );
    void setJacStamp( int numExtVars, int segmentNumber, int vOffset, vector< vector< int > > & segmentJacStamp );
    void loadDAEQVector( int segmentNumber, vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeQVecPtr, double segArea);
    void loadDAEFVector( int segmentNumber, vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeFVecPtr, double segArea);
    void loadDAEdQdx( int segmentNumber, int vOffset, vector< int > & lidIndexVector, vector< vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dQdxMatPtr, double segArea);
    void loadDAEdFdx( int segmentNumber, int vOffset, vector< int > & lidIndexVector, vector< vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dFdxMatPtr, double segArea);

    // constitutive parameters
    double cMem_;     // membrane capacitance
    double gMem_;     // membrane conductance
    double vRest_;    // membrane rest voltage

    // values similar to what Bsrc uses, for handling a single current equation
    // TODO - I don't think I need these, although I am currently using at least expNumVars;
    // need to remove or properly initialize each
    N_UTL_Expression * Exp_ptr;
    int            expNumVars;
    int            expBaseVar;
    int            expNumDdt;
    list<string>   evnList;
    vector<double> expVarDerivs;
    vector<double> myVarVals;
    vector<double> ddtVals;
    double         expVal;

  private:
    vector<string> currentEqus_;	    // list of equations for contribution to the membrane current
    vector<string> indepVars_;	      // list of independent variables for this membrane model
    vector<string> fEqs_;		          // list of unparsed F equations for this membrane model
    vector<string> qEqs_;		          // list of unparsed Q equations for this membrane model
    vector<string> extraFunctions_;   // list of unparsed extra functions for this membrane model
    vector<string> extraParameters_;  // list of unparsed parameters for this membrane model

    vector<RefCountPtr<N_UTL_Expression> > currentEqusExpRCP_;	    // list of rcp to N_UTL_Expressions for contribution to the membrane current
    vector<RefCountPtr<N_UTL_Expression> > indepVarsExpRCP_;	      // list of rcp to N_UTL_Expressions for independent variables for this membrane model
    vector<RefCountPtr<N_UTL_Expression> > fEqsExpRCP_;		          // list of rcp to N_UTL_Expressions for F equations for this membrane model
    vector<RefCountPtr<N_UTL_Expression> > qEqsExpRCP_;		          // list of rcp to N_UTL_Expressions for Q equations for this membrane model
    vector<RefCountPtr<N_UTL_Expression> > extraFunctionsExpRCP_;   // list of rcp to N_UTL_Expressions for extra functions for this membrane model
    vector<RefCountPtr<N_UTL_Expression> > extraParametersExpRCP_;  // list of rcp to N_UTL_Expressions for parameters for this membrane model

    vector<string> paramNames_;
    vector<double> paramValues_;
    vector<string> funcNames_;
    vector< RefCountPtr<N_UTL_Expression> > funcExpRCP_;
    vector< int > funcNumArgs_;
    vector<string> userDefinedNames_;     // used to get minimal collection of user defined variables/names
    map< string, int > indepVarOffset_;   // a map connecting a vars name to its local offset in solution, F, Q and jacobian (An LID map)
    map< int, string > offsetToIndepVar_; // an inverse map of indepVarOffset_
    vector< map< string, int > > systemJacOffset_;
        // for each row in the jacobian, this gives a map relating contributing var name to offset.
        // i.e. systemJacOffset[ row of segment jacStamp ][ "var name" ] = second index into jacobianOffset (for that given row and variable).

    // These are used to hold the variable names in the user expressions.
    // the names can be used with indepVarOffset_ to assign values to the variables.
    vector< vector<string> > currentEqusVarNames_;
    vector< vector<string> > fEqsEqusVarNames_;
    vector< vector<string> > qEqsEqusVarNames_;

    // These hold the variable values needed to evaluate an expression. The values have to be
    // in the same order as the names in the vector< vector< string > > containers above.
    vector< vector<double> > currentEqusVarValues_;
    vector< vector<double> > fEqsEqusVarValues_;
    vector< vector<double> > qEqsEqusVarValues_;

    void convertStringsToExpression( vector< string > & stringInput, vector<RefCountPtr<N_UTL_Expression> > & expRCPOut );
    void consolidateExpressions();
    void substituteParameters( vector<RefCountPtr<N_UTL_Expression> > & expRCP_ );
    void substituteFunctions( vector<RefCountPtr<N_UTL_Expression> > & expRCP_ );
    void convertSymbolsToVars( vector<RefCountPtr<N_UTL_Expression> > & expRCP_, vector< vector<string> > & expNames, vector< vector<double> > & expValsVec );
    void makeSymbolSet();

};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MembraneUserDefined N_DEV_MembraneUserDefined;

#endif
