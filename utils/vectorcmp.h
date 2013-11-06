//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: vectorcmp.h,v $
//
// Purpose        : This program compares two vector files, one from
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
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2007/01/02 19:16:18 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#ifndef Xyce_VectorCompare_h
#define Xyce_VectorCompare_h

#include <iostream>
#include <fstream>
#include <map>

class baseCompare;

class vectorCompare : public baseCompare
{
  public:
    vectorCompare ();
    virtual ~vectorCompare ();

     void run (int iargs, char *cargs[]);

  protected:
  private:
     void doAllocations ();
     void readXyceVectorFile (char *carg, vector<double> & Xvector);
     void readSpiceVectorFile   (char *carg, vector<double> & Svector);
     void outputResult ();

  public:
  protected:
  private:
    vector<double> Svector;
    vector<double> Xvector;

    friend class rhsCompare;
};

#endif
