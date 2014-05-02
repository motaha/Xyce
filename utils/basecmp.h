//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, 2013, Sandia Corporation, Albuquerque, NM.
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
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2013/09/18 20:27:38 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_BaseCompare_h
#define Xyce_BaseCompare_h

#include <string>
#include <map>

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

    std::map<int,int> CtoX;
    std::map<int,int> XtoC;
    std::map<int,std::string> xyceNameMap;

    int isize;
  protected:
  private:

};

#endif
