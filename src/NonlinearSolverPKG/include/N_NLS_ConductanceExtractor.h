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
// Filename       : $RCSfile: N_NLS_ConductanceExtractor.h,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/03/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_ConductanceExtractor_h
#define Xyce_N_NLS_ConductanceExtractor_h

// ---------- Standard Includes ----------
#include<vector>
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_NLS_NonLinearSolver.h>

// ---------- Forward Declarations ----------
class N_TOP_Topology;
class N_IO_CmdParse;
class N_PDS_ParMap;
class N_IO_PkgOptionsMgr;

//-----------------------------------------------------------------------------
// Class         : N_NLS_ConductanceExtractor
// Purpose       :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------

class N_NLS_ConductanceExtractor 
{
public:
  N_NLS_ConductanceExtractor ( N_NLS_NonLinearSolver & nls_, 
                               N_TOP_Topology & top_,
                               N_IO_CmdParse & cp);

  ~N_NLS_ConductanceExtractor ();
  
  bool extract ( 
      const map<string,double> & inputMap,
      vector<double> & outputVector,
      vector< vector<double> > & jacobian );
  
  bool extract ( 
      const string & isoName,
      vector< vector<double> > & jacobian );
  
  
  bool setOptions(const N_UTL_OptionBlock& OB);
  
  // Method to register the package options manager
  bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

protected:
private:
  bool setupIDs_( const map<string,double> & inputMap);
  bool setup_dIdX_Vectors_();
  
  bool setupISO2_IDs_(const string & isoName);
  
  void printJacobian_ (
     const map<string,double> & inputMap,
     vector< vector<double> > & jacobian);
  
  void printPetraObjects_ (const string & varName);
  
  struct N_NLS_ConductanceExtractor_OptionsReg : public N_IO_PkgOptionsReg
  {
   N_NLS_ConductanceExtractor_OptionsReg( N_NLS_ConductanceExtractor * ce )
   : CE(ce)
   {}
  
   bool operator()( const N_UTL_OptionBlock & options )
   { return CE->setOptions( options ); }
  
   N_NLS_ConductanceExtractor * CE;
  };

public:
protected:
private:
  int solutionSize_;
  int debugLevel_;

  // temporary stuff, for use with iso devices:
  map<string, double> varMap_;

  // GID variables
  bool gidsSetUpFlag_;
  vector<int> currentGIDs_;
  vector<int> currentLIDs_;
  vector<int> vsrcPosGIDs_;
  vector<int> vsrcPosLIDs_;

  // package references:
  N_NLS_NonLinearSolver & nls_;
  N_TOP_Topology & top_;
  N_IO_CmdParse & commandLine_;
  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;
  
  // linear system data:
  N_LAS_System  * lasSysPtr_;
  N_ANP_AnalysisInterface * anaIntPtr_;
  N_LOA_Loader  * loaderPtr_;
  N_LAS_Vector  * rhsVectorPtr_;
  N_LAS_Vector  * dfdvVectorPtr_;
  N_LAS_Vector  * NewtonVectorPtr_;
  N_LAS_Vector  * dxdvVectorPtr_;
  N_LAS_Solver  * lasSolverPtr_;

  N_LAS_Vector  * matrixDiagonalPtr_;

  vector <N_LAS_Vector*> dIdxPtrVector_;

  N_LAS_Matrix  * jacobianMatrixPtr_;
  N_LAS_Vector  ** nextSolVectorPtrPtr_;
  N_LAS_Vector  ** currSolVectorPtrPtr_;
  N_LAS_Vector  * savedRHSVectorPtr_;
  N_LAS_Vector  * savedNewtonVectorPtr_;
  N_LAS_Vector  * gradVectorPtr_;

  N_LAS_Vector  * columnVectorPtr_;
  N_PDS_ParMap  * columnMapPtr_;
};

#endif

