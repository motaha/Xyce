//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: TIATest.h,v $
//
// Purpose        : This is the test program for the time integration 
//                  pacakge.
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 6/06/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.1.1 $
//
// Revision Date  : $Date: 2000/09/29 20:30:18 $
//
// Current Owner  : $Author: rjhoeks $
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_TIA_TimeIntegrationAlgorithm.h>
#include <N_ERH_ErrorMgr.h>

// ---------- Forward Declarations -----
class N_DEV_DeviceMgr;
class N_LOA_LoaderMgr;
class N_LOA_Loader;
class N_LAS_Solver;
class N_LAS_System;
class N_NLS_Manager;
class N_NLS_NonLinearSolver;
class N_NLS_NLParams;

//-----------------------------------------------------------------------------
// Class         : TIATestor
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/06/00
//-----------------------------------------------------------------------------
class TIATestor
{
  public:
    TIATestor   () {};
    ~TIATestor  () {};

    bool setTiaParams      ();
    bool setNLParams       ();
    bool doAllocations     ();
    bool doRegistrations   ();
    bool doInitializations ();
    bool doMatrixCreation  ();
    bool doDeAllocations   ();

    bool runTests (int iargs, char *cargs[]);

  protected:
  private:

  public:
  protected:
  private:
    N_LAS_System           * lasSysPtr_; 
    N_NLS_Manager          * nlsMgrPtr_; 
    N_NLS_NonLinearSolver  * nlsPtr_;
    N_LOA_LoaderMgr        * loaderMgrPtr_;
    N_LOA_Loader           * loaderPtr_;
    N_DEV_DeviceMgr        * devPtr_;

    N_TIA_TimeIntegrationAlgorithm tia_;
    N_TIA_TIAParams tiaParams_;
    N_NLS_NLParams * nlParamsPtr_;
};



