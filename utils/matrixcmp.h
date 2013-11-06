//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: matrixcmp.h,v $
//
// Purpose        : This program compares two matrix files, one from
//                  chilespice, one from Xyce.  It uses a merged map
//                  file created by the program "mapMerge.C".  The
//                  map file is assumed to have the name, "mergedMap.txt".
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
// Revision Number: $Revision: 1.4 $
//
// Revision Date  : $Date: 2007/01/02 19:16:17 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#ifndef Xyce_MatrixCompare_h
#define Xyce_MatrixCompare_h

#include <iostream>
#include <fstream>
#include <map>

class baseCompare;

class matrixCompare : public baseCompare
{
  public:
    matrixCompare ();
    virtual ~matrixCompare ();

     void run (int iargs, char *cargs[]);

  protected:
  private:
     void doAllocations ();
     void readXyceMatrixFile 
       (char *carg, 
        vector< vector<double> > &  Xmatrix, 
        vector< vector<int> > & matfill);

     void readSpiceMatrixFile   
       (char *carg, 
        vector< vector<double> > &  Smatrix, 
        vector< vector<int> > & matfill);

     void outputResult ();

  public:
  protected:
  private:
    vector< vector<double> > Smatrix;
    vector< vector<double> > Xmatrix;
    vector< vector<int> > Sfill;
    vector< vector<int> > Xfill;

    vector< vector<int> > signErrors;

    vector< vector<int> > nzFlag;

    friend class rhsCompare;
};

#endif
