//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_PDS_Test.C,v $
//
// Purpose        : This function is the test driver for the ParallelDist
//                  package.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/22/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2001/02/13 19:48:36 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

#include <iostream>
#include <vector>
#include <list>
#include <string>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_PDS_Test.h>
#include <N_PDS_Manager.h>
#include <N_ERH_ErrorMgr.h>

// Some of these are temporary includes - for testing purposes only!

// Class N_PDS_ParTest

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParTest::N_PDS_ParTest
// Purpose       : Default constructor.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

N_PDS_ParTest::N_PDS_ParTest()
{

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParTest::~N_PDS_ParTest
// Purpose       : Default destructor.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

N_PDS_ParTest::~N_PDS_ParTest()
{

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParTest::RunTests
// Purpose       : Performs tests on ParallelDist package.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

int N_PDS_ParTest::RunTests(int iargs, char *cargs[], N_PDS_ParTest *ParTest)

{
  static int iSuccess = false;

  // Error string
  static const string msg("N_PDS_Test::RunTests - ");

#ifdef Xyce_PARALLEL_MPI
  static string              lbMethod = ("PARMETIS");
  static list<string_params> lbParams;
#endif

  cout  << endl <<
    "Welcome to the Xyce(TM) ParallelDist testing program." << endl <<
    endl;

#ifdef Xyce_PARALLEL_MPI

  if (iargs > 3)
  {
    if (iargs%2 != 0)
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg + "wrong number "
                             "of arguments.");

    // Setup the arguments for the manager constructor.
    for (int i = 3; i < iargs; i += 2)
    {
      string_params sp(cargs[i], cargs[i+1]);
      lbParams.push_back(sp);
    }
  }

#endif

  N_PDS_Manager *parMgr = new N_PDS_Manager(atoi(cargs[1])
#ifdef Xyce_PARALLEL_MPI
                                            , lbMethod, lbParams
#endif
                                            );

  iSuccess = (parMgr != NULL);

  parMgr->reportLoadBalance();

  if (iSuccess == false)
    cout << "Test of ParallelDist NOT completed successfully." <<
      endl;
  else
    cout << "Test of ParallelDist completed successfully." << endl;

  delete parMgr;

  return iSuccess;

}

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/12/00
//-----------------------------------------------------------------------------

int main(int iargs, char *cargs[])

{
  N_PDS_ParTest *ParTest = new N_PDS_ParTest();

  if (iargs < 2)
  {
    // Error string
    static const string error_msg("main - Wrong number of arguments.");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, error_msg);
  }

  ParTest->RunTests(iargs, cargs, ParTest);

  delete ParTest;
  exit(0);
}
