//-----------------------------------------------------------------------------
// File          : NLSTest.h
//
// Purpose       : This function is the header file which contains class
//                 definitions for the nonlinear solver package test program.
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------


#ifndef  _NLSTEST_H
#define  _NLSTEST_H

// ---------- Standard Includes ----------
#include <iostream>
#include <vector>
#include <list>
#include <string>

// ----------   Xyce Includes   ----------
#include <N_NLS_Manager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_NLS_Misc.h>
#include <N_UTL_Xyce.h>

// ---------- Forward Declarations ----------
class N_LOA_Loader;
class N_LOA_LoaderMgr;
class N_DEV_DeviceMgr;
class N_LAS_LAFactory;
class N_LAS_MultiVector;
class N_LAS_Matrix;
class N_TIA_TimeIntegrationAlgorithm;

//-----------------------------------------------------------------------------
// Class         : NLSTestor
// Purpose       : This is the top level class for the nonlinear solver testing
//                 program.  The member function, RunTests, is the "main"
//                 function, essentially.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
class NLSTestor
{
  // functions:
public:
  NLSTestor();
  ~NLSTestor();

  bool runTests(int iargs, char *cargs[]);

protected:

private:
  bool doAllocations();
  bool doRegistrations();
  bool doDeAllocations();

  bool doInitialization();
  bool doSolve();


  // attributes
public:

protected:

private:
  N_NLS_Manager         * NLS_Ptr_;
  N_LAS_IterativeSolver * LAS_SolverPtr_;
  N_LAS_Matrix          * LAS_MatrixPtr_;
  N_LAS_MultiVector     * LAS_RHSVecPtr_;
  N_LAS_MultiVector     * LAS_SolVecPtr_;
  N_LOA_LoaderMgr       * LOA_LoaderMgrPtr_;
  N_LOA_Loader          * LOA_LoaderPtr_;
  N_TIA_TimeIntegrationAlgorithm * TIA_Ptr_;
  N_ERH_ErrorMgr        * ERH_Ptr_;
  N_DEV_DeviceMgr       * DEV_Ptr_;

  int iargs;
  char **cargs;
};

#endif


