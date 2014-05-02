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
// Revision Number: $Revision: 1.14.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
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
  MembraneUserDefined(const SolverState & ss1, double cMem, double gMem, double vRest,
                      std::vector<std::string> & currentEqus, std::vector<std::string> & indepVars, std::vector<std::string> & fEqs,
                      std::vector<std::string> & qEqs, std::vector<std::string> & extraFunctions, std::vector<std::string> & extraParameters);
  ~MembraneUserDefined() {}

  void updateSecondaryState( double * staDerivVec );
  void setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp );
  void loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeQVecPtr, double segArea);
  void loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeFVecPtr, double segArea);
  void loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dQdxMatPtr, double segArea);
  void loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dFdxMatPtr, double segArea);

  // constitutive parameters
  double cMem_;     // membrane capacitance
  double gMem_;     // membrane conductance
  double vRest_;    // membrane rest voltage

  // values similar to what Bsrc uses, for handling a single current equation
  // TODO - I don't think I need these, although I am currently using at least expNumVars;
  // need to remove or properly initialize each
  Util::Expression * Exp_ptr;
  int            expNumVars;
  int            expBaseVar;
  int            expNumDdt;
  std::list<std::string>   evnList;
  std::vector<double> expVarDerivs;
  std::vector<double> myVarVals;
  std::vector<double> ddtVals;
  double         expVal;

private:
  std::vector<std::string> currentEqus_;	    // list of equations for contribution to the membrane current
  std::vector<std::string> indepVars_;	      // list of independent variables for this membrane model
  std::vector<std::string> fEqs_;		          // list of unparsed F equations for this membrane model
  std::vector<std::string> qEqs_;		          // list of unparsed Q equations for this membrane model
  std::vector<std::string> extraFunctions_;   // list of unparsed extra functions for this membrane model
  std::vector<std::string> extraParameters_;  // list of unparsed parameters for this membrane model

  std::vector<RefCountPtr<N_UTL_Expression> > currentEqusExpRCP_;	    // list of rcp to Util::Expressions for contribution to the membrane current
  std::vector<RefCountPtr<N_UTL_Expression> > indepVarsExpRCP_;	      // list of rcp to Util::Expressions for independent variables for this membrane model
  std::vector<RefCountPtr<N_UTL_Expression> > fEqsExpRCP_;		          // list of rcp to Util::Expressions for F equations for this membrane model
  std::vector<RefCountPtr<N_UTL_Expression> > qEqsExpRCP_;		          // list of rcp to Util::Expressions for Q equations for this membrane model
  std::vector<RefCountPtr<N_UTL_Expression> > extraFunctionsExpRCP_;   // list of rcp to Util::Expressions for extra functions for this membrane model
  std::vector<RefCountPtr<N_UTL_Expression> > extraParametersExpRCP_;  // list of rcp to Util::Expressions for parameters for this membrane model

  std::vector<std::string> paramNames_;
  std::vector<double> paramValues_;
  std::vector<std::string> funcNames_;
  std::vector< RefCountPtr<N_UTL_Expression> > funcExpRCP_;
  std::vector< int > funcNumArgs_;
  std::vector<std::string> userDefinedNames_;     // used to get minimal collection of user defined variables/names
  std::map< std::string, int > indepVarOffset_;   // a map connecting a vars name to its local offset in solution, F, Q and jacobian (An LID map)
  std::map< int, std::string > offsetToIndepVar_; // an inverse map of indepVarOffset_
  std::vector< std::map< std::string, int > > systemJacOffset_;
  // for each row in the jacobian, this gives a map relating contributing var name to offset.
  // i.e. systemJacOffset[ row of segment jacStamp ][ "var name" ] = second index into jacobianOffset (for that given row and variable).

  // These are used to hold the variable names in the user expressions.
  // the names can be used with indepVarOffset_ to assign values to the variables.
  std::vector< std::vector<std::string> > currentEqusVarNames_;
  std::vector< std::vector<std::string> > fEqsEqusVarNames_;
  std::vector< std::vector<std::string> > qEqsEqusVarNames_;

  // These hold the variable values needed to evaluate an expression. The values have to be
  // in the same order as the names in the std::vector< std::vector< string > > containers above.
  std::vector< std::vector<double> > currentEqusVarValues_;
  std::vector< std::vector<double> > fEqsEqusVarValues_;
  std::vector< std::vector<double> > qEqsEqusVarValues_;

  void convertStringsToExpression( std::vector< std::string > & stringInput, std::vector<RefCountPtr<N_UTL_Expression> > & expRCPOut );
  void consolidateExpressions();
  void substituteParameters( std::vector<RefCountPtr<N_UTL_Expression> > & expRCP_ );
  void substituteFunctions( std::vector<RefCountPtr<N_UTL_Expression> > & expRCP_ );
  void convertSymbolsToVars( std::vector<RefCountPtr<N_UTL_Expression> > & expRCP_, std::vector< std::vector<std::string> > & expNames, std::vector< std::vector<double> > & expValsVec );
  void makeSymbolSet();

};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MembraneUserDefined N_DEV_MembraneUserDefined;

#endif
