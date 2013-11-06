//-----------------------------------------------------------------------------
// File          : LOADTest.h
//
// Purpose       : This function is the header file which contains class
//                 definitions for the nonlinear solver package test 
//                 program.
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------


#ifndef  _LOADTEST_H
#define  _LOADTEST_H

// ---------- Standard Includes ----------
#include <iostream>
#include <vector>
#include <list>
#include <string>

// ----------   Xyce Includes   ----------
#include <N_LOA_LoaderMgr.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Misc.h>

using namespace std;

//-----------------------------------------------------------------------------
// Class         : LOADTestor
// Purpose       : This is the top level class for the Loader Services
//                 testing program.  The member function, RunTests, 
//                 is the "main" function, essentially.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
class LOADTestor
{
  // functions:
  public:
    LOADTestor                ();
    ~LOADTestor               ();

    bool  runTests           (int iargs, char *cargs[]);

  protected:

  private:
    bool  doAllocations      ();
    bool  doRegistrations    ();
    bool  doDeAllocations    ();

    bool  doInitialization   ();
    bool  doLoad             ();

  // attributes
  public:

  protected:

  private:
    N_DEV_DeviceMgr   * DEV_DeviceMgrPtr_;
    N_LOA_LoaderMgr   * LOA_LoaderMgrPtr_;
    N_ERH_ErrorMgr    * ERH_Ptr_;

    int iargs;
    char **cargs;
};

#endif


