//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: basecmp.h,v $
//
// Purpose        : 
//                  
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
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2008/02/21 05:19:06 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#ifndef Xyce_BaseCompare_h
#define Xyce_BaseCompare_h

#include <iostream>
#include <fstream>
#include <map>
#include <vector>

using namespace std;

// uncomment this to have the utilities  read two Xyce files, rather than
// a chilespice file and a Xyce file.
#if 1
#define BOTH_FILES_ARE_XYCE 1
#endif

#if 1
#define FILES_HAVE_DIFFERENT_ORDER 1
#endif

class baseCompare
{
  public:
    baseCompare ();
    virtual ~baseCompare ();

     void readMergedMap ();
  protected:
  private:

  public:
    bool endOfFile;
    bool XyceFileFlag;

    int i1;
    int itmp;
    int itmp2;
    int rowmax;
    int i,j;
    int row, col;
    double val;

    char txt[128];

    map<int,int> CtoX;
    map<int,int> XtoC;
    map<int,string> xyceNameMap;

    int isize;
  protected:
  private:

};

#endif
