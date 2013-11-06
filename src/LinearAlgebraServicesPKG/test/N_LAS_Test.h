//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_LAS_Test.h,v $
//
// Purpose        : Specification file for testing the LinearAlgebraServices
//                  package.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/23/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2001/02/13 19:27:16 $
//
// Current Owner  : $Author $
//-------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Test_h
#define Xyce_N_LAS_Test_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_PDS_Manager.h>
#include <N_LAS_LAFactory.h>
#include <N_UTL_Xyce.h>

//-----------------------------------------------------------------------------
// Class         : N_LAS_LATest
// Purpose       : This is the top level class for the linear-algebra testing
//                 program.  The member function, RunTests, is the "main"
//                 function, essentially.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------

class N_LAS_LATest
{

public:

  // Default constructor & destructor.
  N_LAS_LATest();
  ~N_LAS_LATest();

  int RunTests(int iargs, char *cargs[], N_LAS_LATest *LATest);

  int vectorTests(N_PDS_ParMap *parMap, N_LAS_LAFactory *factory,
                  int numVectors);
  int matrixVectorTests(N_PDS_ParMap *parMap, N_LAS_LAFactory *factory,
                        int numVectors);

  int solverTests(N_PDS_ParMap *parMap, N_LAS_LAFactory *factory,
                  int numVectors);

protected:

private:

};

//-----------------------------------------------------------------------------

#endif
