//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: basecmp.C,v $
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 09/09/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2010/04/29 20:38:11 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

// ----------   Standard Includes   ----------
#include <string>

#include <string.h>
#include <math.h>
#ifdef __FreeBSD__
#include <stdio.h>
#endif
#ifdef _LINUX_
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include "basecmp.h"

using namespace std;

#define LIMIT 1.0e-25

//-----------------------------------------------------------------------------
// Function      : baseCompare::baseCompare
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
baseCompare::baseCompare ()
  : rowmax(0),
    endOfFile(false),
#ifdef BOTH_FILES_ARE_XYCE
    XyceFileFlag(true),
#else
    XyceFileFlag(false),
#endif
    i1(0),
    itmp(0),
    itmp2(0),
    i(0),
    j(0),
    row(0), 
    col(0),
    isize(0),
    val(0.0)
{
  for (i=0;i<128;i++) txt[i] = 0;
}

//-----------------------------------------------------------------------------
// Function      : baseCompare::~baseCompare
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
baseCompare::~baseCompare ()
{
}

//-----------------------------------------------------------------------------
// Function      : baseCompare::readMergedMap
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
void baseCompare::readMergedMap()
{
  // obtain the CtoX map.

#ifndef FILES_HAVE_DIFFERENT_ORDER 

  CtoX[ 0] =  0;  // this is not used.
  XtoC[ 0] =  0;
  CtoX[-1] = -1;  // this is not used.
  XtoC[-1] = -1;


  // Test to see if namesMap exists:
  FILE * testFile;
  if ( (testFile=fopen("namesMap.txt", "r")) == 0)
  {
    cout << "Unable to file namesMap.txt file\n";
    exit(0);
  }
  int ierr = fclose( testFile );


  ifstream * inStream = new ifstream("namesMap.txt");

  endOfFile = (*inStream).eof();

  isize = 0;
  cout << "\nReading the mergedMap.txt (namesMap.txt) file... ";
  cout << endl;
  do
  {
    int int1;
    string name;

    (*inStream) >> int1  >> name;
    cout << "\t"<<int1<<"\t"<<name<<endl;

    endOfFile = (*inStream).eof();
    if (endOfFile) break;

    CtoX[int1] = int1;
    XtoC[int1] = int1;
    xyceNameMap[int1] = name;

    if (int1 > isize) isize = int1;

    endOfFile = (*inStream).eof();
  } while (!endOfFile);

  delete inStream;

  cout << " ...done\n";

  cout << "\nisize = " << isize << endl;

#else
  CtoX[ 0] = -1;  // this is not used.
  XtoC[-1] =  0;

  // Test to see if mergedMap file exists:
  FILE * testFile;
  if ( (testFile=fopen("mergedMap.txt", "r")) == 0)
  {
    cout << "Unable to file mergedMap.txt file\n";
    exit(0);
  }
  int ierr = fclose( testFile );

  ifstream * inStream = new ifstream("mergedMap.txt");

  endOfFile = (*inStream).eof();

  isize = 0;
  cout << "\nReading the mergedMap.txt file... ";
  do 
  {
    int int1, int2;
    string name;

    (*inStream) >> int1  >> int2 >> name;
    cout << "\t"<<int1<<"\t"<<int2<<"\t"<<name<<endl;

    endOfFile = (*inStream).eof();
    if (endOfFile) break;

    CtoX[int1] = int2;
    XtoC[int2] = int1;
    xyceNameMap[int2] = name;

    if (int1 > isize) isize = int1;
    if (int2 > isize) isize = int2;

    endOfFile = (*inStream).eof();
  } while (!endOfFile);

  delete inStream;
  cout << " ...done\n";

  cout << "\nisize = " << isize << endl;
#endif


#if 1


  map<int,int>::iterator iter = XtoC.begin();
  map<int,int>::iterator end = XtoC.end();

  for (int i1=0;iter!=end;++iter, ++i1)
  {
    cout << "XtoC first=" << iter->first << " second = " << iter->second ;
    cout << "   CtoX[XtoC] = " <<  CtoX[iter->second] << endl;
  }
  cout << endl;

  iter = CtoX.begin();
  end = CtoX.end();
  for (int i1=0;iter!=end;++iter, ++i1)
  {
    cout << "CtoX first=" << iter->first << " second = " << iter->second ;
    cout << "   XtoC[CtoX] = " <<  XtoC[iter->second] << endl;
  }
  cout << endl;

#endif


}

