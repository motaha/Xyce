//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: rhscmp.C,v $
//
// Purpose        : This file compares rhs vectors.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
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
#include "vectorcmp.h"
#include "matrixcmp.h"
#include "rhscmp.h"

using namespace std;

#define LIMIT 1.0e-25

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
int main (int iargs, char *cargs[])
{
  rhsCompare RC;

  RC.run(iargs,cargs);

}

//-----------------------------------------------------------------------------
// Function      : rhsCompare::rhsCompare
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
rhsCompare::rhsCompare ()
  : matrixCompare(), vectorCompare()
{
  S_sol_vector.clear();
}

//-----------------------------------------------------------------------------
// Function      : rhsCompare::~rhsCompare
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
rhsCompare::~rhsCompare ()
{
}

//-----------------------------------------------------------------------------
// Function      : rhsCompare::doAllocations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
void rhsCompare::doAllocations ()
{
  int isize = matrixCompare::isize;
  cout << "Allocating arrays... " << endl;
  S_sol_vector.resize(isize+10,0.0);
  Jx_vector.resize(isize+10,0.0);
  cout << "  ...done." << endl;
}

//-----------------------------------------------------------------------------
// Function      : rhsCompare::matvec
// Purpose       :
// Special Notes : This is goofy, as the Smatrix is set up to start its 
//                 row, col indices from 1, but the vectors (spice and xyce)
//                 start from 0.  So, row and col have to be adjusted during
//                 the matvec.  FIX LATER.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
void rhsCompare::matvec ( vector< vector<double> > & mat, 
                          vector< double > & vec,
                          vector< double > & mv )
{
  int i,j;
  int imax=mat.size()-1;

  for (i=0;i<imax;++i)
  {
    int jmax=mat[i].size()-1;
    mv[i] = 0.0;
    for (j=0;j<jmax;++j)
    {
      int row = i+1;
      int col = j+1;
      double val = (mat[row][col]) * (vec[j]);
      mv[i] += val;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : rhsCompare::calcJx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
void rhsCompare::calcJx ()
{
  matvec(Smatrix, S_sol_vector, Jx_vector);
}


//-----------------------------------------------------------------------------
// Function      : rhsCompare::sumRHS
// Purpose       : sets up:  -f + Jx_new - Jx_old
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
void rhsCompare::sumRHS()
{
  int index;
  int imax=Svector.size();

  for (index=0;index<imax;++index)
  {
    double val = Jx_vector[index];
    Svector[index] -= val;
  }
}

//-----------------------------------------------------------------------------
// Function      : rhsCompare::outputResult
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
void rhsCompare::outputResult ()
{
  // open result file
  ofstream * outStream = new ofstream("rhsresult.txt");
  
  double diffmax = 0.0;
  double scaledDiffMax = 0.0;
  int scaledMaxRow = 0; // This is the index of maximum scaled difference.
  int absMaxRow = 0;    // This is the index of maximum absolute difference.
  int index=0;
  map<int,int> & Xyce2Spice = matrixCompare::XtoC;

  char txt[256];  for(index=0;index<256;index++) txt[index] = 0;

  for (index=0;index<=matrixCompare::rowmax;index++)
  {
#ifdef BOTH_FILES_ARE_XYCE
     int Crow = index;
#else
     int Crow = Xyce2Spice[index]-1;
#endif

     double diff = Xvector[index] - Svector[Crow];

     double scaledDiff = 0.0;
     if (diff!=0.0 && Xvector[index] != 0.0)
       scaledDiff = diff/Xvector[index];

     sprintf(txt,"%4d ",index);
     (*outStream) << txt;
     sprintf(txt,"%25.18e ",Xvector[index]);
     (*outStream) << txt;
     sprintf(txt,"%4d ",Crow);
     (*outStream) << txt;
     sprintf(txt,"%25.18e ",Svector[Crow]);
     (*outStream) << txt;
     sprintf(txt,"%25.18e ",diff);
     (*outStream) << txt;
     sprintf(txt,"%25.18e",scaledDiff);
     (*outStream) << txt;
     (*outStream) << "\t" << matrixCompare::xyceNameMap[index];
     (*outStream) << endl;

     if (fabs(diff) > diffmax) 
     {
       diffmax = fabs(diff);
       absMaxRow = index;
     }

     if (fabs(scaledDiff) > scaledDiffMax) 
     {
       scaledDiffMax = fabs(scaledDiff);
       scaledMaxRow = index;
     }

   }

  delete outStream;

  sprintf(txt,"\nThe maximum difference is: %15.7e\n", diffmax);
  cout << txt;
  sprintf(txt,"Xyce vector entry: %4d %15.7e\n",
               absMaxRow,Xvector[absMaxRow]);
  cout << txt;
#ifdef BOTH_FILES_ARE_XYCE
  sprintf(txt,"Xyce vector entry: %4d %15.7e\n\n",
               Xyce2Spice[absMaxRow],
               Svector[absMaxRow]);
#else
  sprintf(txt,"3f5  vector entry: %4d %15.7e\n\n",
               Xyce2Spice[absMaxRow],
               Svector[Xyce2Spice[absMaxRow]-1]);
#endif
  cout << txt;
  sprintf(txt,"name for this entry: %15s\n",matrixCompare::xyceNameMap[absMaxRow].c_str());
  cout << txt;



  sprintf(txt,"\nThe maximum scaled difference is: %15.7e\n", scaledDiffMax);
  cout << txt;
  sprintf(txt,"Xyce vector entry: %4d %15.7e\n", scaledMaxRow,Xvector[scaledMaxRow]);
  cout << txt;
#ifdef BOTH_FILES_ARE_XYCE
  sprintf(txt,"Xyce vector entry: %4d %15.7e\n\n",
               Xyce2Spice[scaledMaxRow],
               Svector[scaledMaxRow]);
#else
  sprintf(txt,"3f5  vector entry: %4d %15.7e\n\n",
               Xyce2Spice[scaledMaxRow],
               Svector[Xyce2Spice[scaledMaxRow]-1]);
#endif
  cout << txt;
  sprintf(txt,"name for this entry: %15s\n",matrixCompare::xyceNameMap[scaledMaxRow].c_str());
  cout << txt << endl;

}

//-----------------------------------------------------------------------------
// Function      : rhsCompare::run 
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 12/26/06
//-----------------------------------------------------------------------------
void rhsCompare::run (int iargs, char *cargs[])
{
  // first check that there are enough arguments:
  if (iargs != 5)
  {
     cout << "Usage:  rhscmp spice_matrix spice_rhs spice_old_sol Xyce_rhs <return>\n";
     exit(0);
  }

  matrixCompare::readMergedMap ();
  matrixCompare::doAllocations ();

  vectorCompare::isize = matrixCompare::isize;
  vectorCompare::xyceNameMap = matrixCompare::xyceNameMap;
  vectorCompare::doAllocations ();

  doAllocations ();


  // read in the various files.
  readSpiceMatrixFile (cargs[1], Smatrix, Sfill);
  readSpiceVectorFile (cargs[2], Svector);
  readSpiceVectorFile (cargs[3], S_sol_vector);
  readXyceVectorFile (cargs[4], Xvector);

  calcJx ();
  sumRHS();

  char txt[256];  for(int index=0;index<256;index++) txt[index] = 0;

  map<int,int> & Spice2Xyce = matrixCompare::CtoX;
  map<int,string> & xyceNameMap = matrixCompare::xyceNameMap;


  map<int,string>::iterator iter = xyceNameMap.begin();
  map<int,string>::iterator end  = xyceNameMap.end ();

  for (;iter!=end;++iter)
  {
    cout << "xyceNameMap["<< iter->first << "] = " << iter->second << endl;
  }
  cout << endl;

  for (int i1=0;i1< matrixCompare::isize; ++i1)
  {
    int Xrow = Spice2Xyce[i1+1];
    //if (xyceNameMap.find(Xrow)!=xyceNameMap.end() )
    if (Xrow >= 0)
    {
      string name = xyceNameMap[Xrow];
      sprintf(txt ,"%15s\t", name.c_str());
    }
    else
    {
      sprintf(txt ,"xrow=%2d    \t",Xrow);
    }
    cout << txt;

    sprintf(txt ,"Jx[%2d] = %20.13e ",i1, Jx_vector[i1]);
    cout << txt;

    sprintf(txt ," S_rhs[%2d] = %20.13e ",i1, Svector[i1]);
    cout << txt;

    sprintf(txt ," X_rhs[%2d] = %20.13e ",Xrow, Xvector[Xrow]);
    cout << txt;

    double diff = fabs(Xvector[Xrow]-Svector[i1]);
    sprintf(txt ," diff = %20.13e ",diff);
    cout << txt << endl;
  }

  if (matrixCompare::rowmax==0) matrixCompare::rowmax=matrixCompare::isize;
  outputResult ();

}


