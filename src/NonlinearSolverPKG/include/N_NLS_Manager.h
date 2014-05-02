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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_Manager.h,v $
//
// Purpose        : Defines N_NLS_Manager class.
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.110 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_Manager_h
#define Xyce_N_NLS_Manager_h

// ---------- Standard Includes ----------
#include <map>
#include <vector>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_IO_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_OptionBlock.h>
#include <N_NLS_ReturnCodes.h>

#include <N_ANP_fwd.h>
#include <N_IO_PkgOptionsMgr.h>

// Nonlinear solver
class N_NLS_NonLinearSolver;

// Loader for RHS and Jacobian
class N_LOA_Loader;

// Linear Algebra Support
class N_LAS_Matrix;
class N_LAS_Vector;
class N_LAS_System;
class N_LAS_PrecondFactory;

class N_NLS_TwoLevelNewton;

class N_NLS_Sensitivity;
class N_NLS_ConductanceExtractor;
class N_NLS_NonLinInfo;

class N_PDS_Manager;

// ---------- Enum Definitions ----------

//-----------------------------------------------------------------------------
// Enum          : AnalysisMode
// Description   : These anlysis modes influence the choice of parameters for
//                 the nonlinear solver. The mode is set by the
//                 setAnalysisMode() function in the N_NLS_Manager class.
//-----------------------------------------------------------------------------
enum AnalysisMode
{
  DC_OP,
  DC_SWEEP,
  TRANSIENT,
  HB_MODE,
  NUM_MODES
};

//-----------------------------------------------------------------------------
// Class         : N_NLS_Manager
// Purpose       : Interface to the nonlinear solvers (NLS). All communication
//                 between Xyce and the nonlinear solver should be via this
//                 interface.
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------

class N_NLS_Manager
{
public:

  N_NLS_Manager(N_IO_CmdParse & cp);
  ~N_NLS_Manager();

  bool setOptions(const N_UTL_OptionBlock& OB );
  bool setTimeOptions(const N_UTL_OptionBlock& OB );
  bool setTranOptions(const N_UTL_OptionBlock& OB );
  bool setHBOptions(const N_UTL_OptionBlock& OB );
  bool getHBOptions(N_UTL_OptionBlock& HBOB);
  bool setTwoLevelOptions (const N_UTL_OptionBlock & OB);
  bool setTwoLevelTranOptions (const N_UTL_OptionBlock & OB);
  bool setSensOptions (const N_UTL_OptionBlock & OB);
  bool setSensitivityOptions (const N_UTL_OptionBlock & OB);
  bool setLinSolOptions(const N_UTL_OptionBlock& OB );
  bool setLocaOptions(const N_UTL_OptionBlock& OB );
  bool setTwoLevelLocaOptions(const N_UTL_OptionBlock& OB );
  bool setDCOPRestartOptions (const N_UTL_OptionBlock& OB );
  bool setICOptions (const N_UTL_OptionBlock& OB );
  bool setNodeSetOptions (const N_UTL_OptionBlock& OB );
  bool registerLoader(N_LOA_Loader* ptr);
  bool registerOutputMgr(N_IO_OutputMgr * outputPtr);
  bool registerRHSVector(N_LAS_Vector* ptr);
  bool registerLinearSystem(N_LAS_System* ptr);
  bool registerPrecondFactory( const RCP<N_LAS_PrecondFactory>& ptr);

  bool registerAnalysisInterface(N_ANP_AnalysisInterface* ptr);
  bool registerTopology(N_TOP_Topology * topPtr);
  bool registerPkgOptionsMgr( N_IO_PkgOptionsMgr *pkgOptPtr );
  bool registerParallelMgr(N_PDS_Manager * pdsMgrPtr);
  void setReturnCodes (const N_NLS_ReturnCodes & retCodeTmp);
  N_NLS_ReturnCodes getReturnCodes() const;

  bool initializeAll();
  int solve();
  bool isFirstContinuationParam();
  bool isFirstSolveComplete();
  int getContinuationStep ();
  bool getLocaFlag ();
  int getNumIterations();
  int getNumResidualLoads();
  int getNumJacobianLoads();
  int getNumLinearSolves();
  int getNumFailedLinearSolves();
  int getNumJacobianFactorizations();
  unsigned int getTotalNumLinearIters();
  double getTotalLinearSolveTime();
  double getTotalResidualLoadTime();
  double getTotalJacobianLoadTime();
  void setAnalysisMode(AnalysisMode mode);

  // This is a more extensive version of setAnalysisMode.
  // It makes it like we are starting over.
  void resetAll(AnalysisMode mode);

  int  getCouplingMode ();

  bool getTwoLevelSolveFlag ();

  void getNonLinInfo (N_NLS_NonLinInfo & nlInfo);

  bool enableSensitivity ();
  bool calcSensitivity (
      std::vector<double> & objectiveVec,
      std::vector<double> & dOdpVec, std::vector<double> & dOdpAdjVec,
      std::vector<double> & scaled_dOdpVec, std::vector<double> & scaled_dOdpAdjVec);

  bool obtainConductances (
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian );

  bool obtainConductances (
        const std::string & isoName,
        std::vector< std::vector<double> > & jacobian );

  void setMatrixFreeFlag(bool matrixFreeFlag);
  void allocateTranSolver();

  // get the norm of F form last solve
  double getMaxNormF() const;

  // get the vector index associated with the max norm.
  int getMaxNormFindex () const;

protected:
private:
  bool allocateSolver_ ();
  void usingNox_ ();
  bool matrixFreeFlag_;

  bool setupSensitivity_ ();

  struct N_NLS_Manager_OptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_OptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_TranOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_TranOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setTranOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_HBOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_HBOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setHBOptions( options ); }

    N_NLS_Manager * Mgr;
  };


  struct N_NLS_Manager_LSOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_LSOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setLinSolOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_LocaOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_LocaOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setLocaOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_SensOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_SensOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setSensOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_SensitivityOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_SensitivityOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setSensitivityOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_TwoLvlOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_TwoLvlOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setTwoLevelOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_TwoLvlTranOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_TwoLvlTranOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setTwoLevelTranOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_TimeOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_TimeOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setTimeOptions( options ); }

    N_NLS_Manager * Mgr;
  };

  struct N_NLS_Manager_DCOPRestartOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_DCOPRestartOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setDCOPRestartOptions (options); }

    N_NLS_Manager * Mgr;
  };


  struct N_NLS_Manager_ICOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_ICOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setICOptions (options); }

    N_NLS_Manager * Mgr;
  };


  struct N_NLS_Manager_NodeSetOptionsReg : public N_IO_PkgOptionsReg
  {
    N_NLS_Manager_NodeSetOptionsReg( N_NLS_Manager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setNodeSetOptions (options); }

    N_NLS_Manager * Mgr;
  };

public:
protected:

private:

  // Pointer to the nonlinear solver that is being used.
  N_NLS_NonLinearSolver *nlsPtr_;

  N_NLS_ConductanceExtractor * conductanceExtractorPtr_;
  N_NLS_Sensitivity          * nlsSensitivityPtr_;
  N_TOP_Topology             * topPtr_;

  N_ANP_AnalysisInterface  * anaIntPtr_;
  N_LOA_Loader   * loaderPtr_;
  N_LAS_System   * lasSysPtr_;
  N_LAS_Vector   * rhsVecPtr_;
  N_IO_OutputMgr * outputPtr_;
  N_PDS_Manager  * pdsMgrPtr_;
  RCP<N_LAS_PrecondFactory> lasPrecPtr_;

  N_IO_CmdParse & commandLine_;
  // package options manager
  N_IO_PkgOptionsMgr * pkgOptMgrPtr_;

  // Flag to determine if we are doing 2-level newton or not.
  bool twoLevelNewtonFlag_;

  // Flag to determine if NOX is the solver in use.
  bool noxFlag_;
  bool noxFlagInner_; //for 2-level newton, option for inner loop to use nox.
  bool noxFlagTransient_;  // Use nox in transient phase of calculation.

  // container to hold netlist option blocks until we know which
  // solver to allocate.
  std::map<std::string,N_UTL_OptionBlock> optionBlockMap_;

  bool setupSensFlag_;

  bool initializeAllFlag_;

  // Return Codes.
  N_NLS_ReturnCodes retCodes_;

  Teuchos::RefCountPtr<N_UTL_Expression> exprPtr;
};

#endif
