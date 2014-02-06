//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: IO_OutputMgr.C,v $
// Purpose       : This file contains unit tests for the IO_OutputMgr interfaces
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.4 $
// Revision Date  : $Date: 2011/02/18 21:49:15 $
// Current Owner  : $Author: hkthorn $
//-----------------------------------------------------------------------------

#include<Teuchos_UnitTestHarness.hpp>
#include<N_IO_OutputMgr.h>
#include<Teuchos_RefCountPtr.hpp>

TEUCHOS_UNIT_TEST( IO_OutputMgr, stdHeaderMangler ) {
  { // last column gets no post-padding
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo" );
  }
  // Left justification checks:
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 2, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 3, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 5, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo  " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  // Center justification checks:
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 2, N_IO_HeaderData::JUSTIFICATION_CENTER);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 3, N_IO_HeaderData::JUSTIFICATION_CENTER);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 4, N_IO_HeaderData::JUSTIFICATION_CENTER);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 5, N_IO_HeaderData::JUSTIFICATION_CENTER);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], " foo " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 6, N_IO_HeaderData::JUSTIFICATION_CENTER);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], " foo  " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  // Right justification checks:
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 2, N_IO_HeaderData::JUSTIFICATION_RIGHT);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], " foo" );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 3, N_IO_HeaderData::JUSTIFICATION_RIGHT);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], " foo" );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 4, N_IO_HeaderData::JUSTIFICATION_RIGHT);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], " foo" );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 5, N_IO_HeaderData::JUSTIFICATION_RIGHT);
    std::string bar = "bar";
    headerData.addData( bar, 4, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "  foo" );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  // Last column checks:
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 5, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string bar = "bar";
    headerData.addData( bar, 5, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo  " );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 5, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string bar = "bar";
    headerData.addData( bar, 5, N_IO_HeaderData::JUSTIFICATION_CENTER);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo  " );
    TEST_EQUALITY_CONST( mangledStrings[1], " bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 5, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string bar = "bar";
    headerData.addData( bar, 5, N_IO_HeaderData::JUSTIFICATION_RIGHT);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo  " );
    TEST_EQUALITY_CONST( mangledStrings[1], "  bar" );
  }
  // delimiter checking:
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 5, N_IO_HeaderData::JUSTIFICATION_LEFT);
    std::string bar = "bar";
    headerData.addData( bar, 5, N_IO_HeaderData::JUSTIFICATION_RIGHT);
    std::string delimiter = ",";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo," );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 0, N_IO_HeaderData::JUSTIFICATION_INVALID);
    std::string bar = "bar";
    headerData.addData( bar, 0, N_IO_HeaderData::JUSTIFICATION_INVALID);
    std::string delimiter = "--hello--";
    std::vector<std::string> mangledStrings = stdHeaderMangler(headerData, delimiter);
    TEST_EQUALITY_CONST( mangledStrings[0], "foo--hello--" );
    TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
  // Error report:
  { 
    N_IO_HeaderData headerData;
    std::string foo = "foo";
    headerData.addData( foo, 0, N_IO_HeaderData::JUSTIFICATION_INVALID);
    std::string bar = "bar";
    headerData.addData( bar, 0, N_IO_HeaderData::JUSTIFICATION_INVALID);
    std::string delimiter = "";
    std::vector<std::string> mangledStrings;
    // 09/24/08 tscoffe:  we can't call this because the Error handler will segfault
    //TEST_THROW( mangledStrings = stdHeaderMangler(headerData, delimiter), std::logic_error);
    //TEST_EQUALITY_CONST( mangledStrings[0], "foo" );
    //TEST_EQUALITY_CONST( mangledStrings[1], "bar" );
  }
}

TEUCHOS_UNIT_TEST( IO_OutputMgr, getWidthFromStaticIndex ) {
  int W = 0;
  W = getWidthFromStaticIndex( 0, "" );
  TEST_EQUALITY_CONST( W, 8 );
  W = getWidthFromStaticIndex( 1000000, "" );
  TEST_EQUALITY_CONST( W, 8 );
  W = getWidthFromStaticIndex(9999999, "" );
  TEST_EQUALITY_CONST( W, 8 );
  W = getWidthFromStaticIndex(10000000, "" );
  TEST_EQUALITY_CONST( W, 9 );
  W = getWidthFromStaticIndex(100000000, "" );
  TEST_EQUALITY_CONST( W, 10 );
  W = getWidthFromStaticIndex(100000000, "," );
  TEST_EQUALITY_CONST( W, 0 );
  W = getWidthFromStaticIndex(-5, "" );
  TEST_EQUALITY_CONST( W, 8 );
  W = getWidthFromStaticIndex(-100000000, "" );
  TEST_EQUALITY_CONST( W, 8 );
}

