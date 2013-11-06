//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_PDS_Test.h,v $
//
// Purpose        : Specification file for testing the ParallelDist package.
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
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2001/02/13 19:48:36 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Test_h
#define Xyce_N_PDS_Test_h

// ---------- Standard Includes ----------

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_PDS_Manager.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>

//-----------------------------------------------------------------------------
// Class         : N_PDS_ParTest
// Purpose       : This is the top level class for the parallel testing
//                 program.  The member function, RunTests, is the "main"
//                 function, essentially.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

class N_PDS_ParTest
{

public:

  // Default constructor & destructor.
  N_PDS_ParTest();
  ~N_PDS_ParTest();

  int RunTests(int iargs, char *cargs[], N_PDS_ParTest *ParTest);

protected:

private:
  N_PDS_Manager* parMgrPtr;

};

//-----------------------------------------------------------------------------

#endif
