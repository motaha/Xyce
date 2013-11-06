//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: rhscmp.h,v $
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 12/26/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2007/01/02 19:16:17 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#ifndef Xyce_rhsCompare_h
#define Xyce_rhsCompare_h

#include <iostream>
#include <fstream>
#include <map>

class matrixCompare;
class vectorCompare;

class rhsCompare : public matrixCompare, public vectorCompare
{
  public:
    rhsCompare ();
    virtual ~rhsCompare ();

     void run (int iargs, char *cargs[]);

     void doAllocations ();
     void outputResult ();

     void matvec ( vector< vector<double> > & mat, 
                   vector< double > & vec,
                   vector< double > & mv );

     void calcJx ();
     void sumRHS ();

  protected:
  private:

  public:
  protected:
  private:
     vector <double> S_sol_vector;
     vector <double> Jx_vector;

};

#endif
