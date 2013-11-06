//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: matrixcmp.C,v $
//
// Purpose        : This program compares two matrix files, one from
//                  spice, one from Xyce.  It uses a merged map
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
// Revision Number: $Revision: 1.14 $
//
// Revision Date  : $Date: 2007/03/03 18:48:14 $
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
#include "matrixcmp.h"

using namespace std;

#define LIMIT 1.0e-25

#ifndef RHSCMP_BUILD
//-----------------------------------------------------------------------------
// Function      : main
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
int main (int iargs, char *cargs[])
{
  matrixCompare MC;
  MC.run(iargs,cargs);
}
#endif

//-----------------------------------------------------------------------------
// Function      : matrixCompare::matrixCompare
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
matrixCompare::matrixCompare ()
  : baseCompare()
{

}

//-----------------------------------------------------------------------------
// Function      : matrixCompare::~matrixCompare
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
matrixCompare::~matrixCompare ()
{
  Xmatrix.clear();
  Smatrix.clear();
}

//-----------------------------------------------------------------------------
// Function      : matrixCompare::doAllocations
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
void matrixCompare::doAllocations ()
{
  Smatrix.resize(isize+10);
  Xmatrix.resize(isize+10);
  Sfill.resize(isize+10);
  Xfill.resize(isize+10);
  nzFlag.resize(isize+10);
  signErrors.resize(isize+10);

  for (itmp = 0; itmp < isize+10; itmp++)
  {
    Smatrix[itmp].resize(isize+10,0.0);
    Xmatrix[itmp].resize(isize+10,0.0);
    Sfill[itmp].resize(isize+10,0);
    Xfill[itmp].resize(isize+10,0);
    nzFlag [itmp].resize(isize+10,0);
    signErrors[itmp].resize(isize+10,0);
  }
  cout << "Done with allocations\n";
}

//-----------------------------------------------------------------------------
// Function      : matrixCompare::run 
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
void matrixCompare::run (int iargs, char *cargs[])
{
#ifdef BOTH_FILES_ARE_XYCE 
  cout << "Assuming both files are Xyce files.\n";
#else
  cout << "Assuming one first file is from spice, second file from Xyce.\n";
#endif

  // first check that there are 2 arguments:
  if (iargs != 3)
  {
#ifdef BOTH_FILES_ARE_XYCE 
     cout << "Usage:  matcmp Xyce_matrix_file_gold Xyce_matrix_file_test <return>\n";
#else
     cout << "Usage:  matcmp chilspice_matrix_file Xyce_matrix_file <return>\n";
#endif
     exit(0);
  }

  readMergedMap ();
  doAllocations ();

  // read in the first matrix file:    (by default a spice file)
 
  if (XyceFileFlag)
  {
    readXyceMatrixFile (cargs[1], Smatrix, Sfill);
  }
  else
  {
    readSpiceMatrixFile (cargs[1], Smatrix, Sfill);
  }

  // read in the second matrix file:  (by default a Xyce file)
  XyceFileFlag = true;
  if (XyceFileFlag)
  {
    readXyceMatrixFile (cargs[2], Xmatrix, Xfill);
  }
  else
  {
    readSpiceMatrixFile (cargs[2], Xmatrix, Xfill);
  }

  if (rowmax==0) rowmax=isize;
  outputResult ();
}

//-----------------------------------------------------------------------------
// Function      : matrixCompare::readXyceMatrixFile
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
void matrixCompare::readXyceMatrixFile 
  (char *carg, 
   vector< vector<double> > & Xmatrix1, vector< vector<int> > & matfill)
{
  // open and read xyce file.
  cout << "Reading the Xyce matrix file... ";

  ifstream * inStream = new ifstream(carg);

  int imax;
  
  (*inStream) >> imax;

  endOfFile = (*inStream).eof();
  while (!endOfFile)
  {

    (*inStream) >> row;
    (*inStream) >> col;
    (*inStream) >> val;

    if (row > isize)
    {
      cout << "\nERROR: row index = " << row << " is larger than the maximum size from the map file.\n";
      exit(0);
    }

    if (col > isize)
    {
      cout << "\nERROR: col index = " << col << " is larger than the maximum size from the map file.\n";
      exit(0);
    }

    if(fabs(val) > LIMIT) { Xmatrix1[row][col] = val; }
    else                  { Xmatrix1[row][col] = 0.0; }
    matfill[row][col] = 1;

    nzFlag[row][col] = 1;

    //cout << row<<"\t"<<col<<"\t"<<val<<endl;

    endOfFile = (*inStream).eof();
  }
  delete inStream;

  cout << " ...done\n";
}

//-----------------------------------------------------------------------------
// Function      : matrixCompare::readSpiceMatrixFile
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
void matrixCompare::readSpiceMatrixFile 
  (char *carg, vector< vector<double> > & Smatrix1, vector< vector<int> > & matfill)
{
  // open and read in the 3f5 file.  first get past the first two lines.
  // (this assumes the matrices are same size)

  cout << "Reading the spice matrix file... ";

  ifstream * inStream = new ifstream(carg);

  string tmpstr;
  (*inStream) >> tmpstr; //cout << tmpstr <<endl;
  (*inStream) >> rowmax; //cout << rowmax <<"\t";
  (*inStream) >> tmpstr; //cout << tmpstr <<endl;

  int complexflag = 1;

  if(tmpstr=="real") complexflag = 0;

  endOfFile = (*inStream).eof();
  while (!endOfFile)
  {
    (*inStream) >> row;
    (*inStream) >> col;
    (*inStream) >> val;

    //cout << row<<"\t"<<col<<"\t"<<val;

    if ( !(row==0 || col==0) )
    {
      if(fabs(val) > LIMIT) { Smatrix1[row][col] = val; }
      else                  { Smatrix1[row][col] = 0.0; }
      matfill[row][col] = 1;
    }
    else 
    {
      endOfFile = true;
    }

    if (complexflag == 1) 
    {
      (*inStream) >> val;
      //cout << "\t"<<val;
    }
    //cout << endl;

    endOfFile = (*inStream).eof();
  }

  rowmax -= 1;

  delete inStream;
  cout << " ... done\n";
}

//-----------------------------------------------------------------------------
// Function      : matrixCompare::outputResult
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
void matrixCompare::outputResult ()
{
  // open result file
  ofstream * outStream = new ofstream("matresult.txt");
  ofstream * outSt2    = new ofstream("matsignerrors.txt");

  double diffmax = 0.0;
  double scaledDiffMax = 0.0;
  int srow = 0;
  int scol = 0;
  int Crow, Ccol;
  int mC, mX;

  for (i=0;i<=rowmax;i++)
  {
    for (j=0;j<=rowmax;j++)
    {
       Crow = XtoC[i];
       Ccol = XtoC[j];

       if (nzFlag[i][j])
       {
         double diff = Xmatrix[i][j] - Smatrix[Crow][Ccol];

         double scaledDiff = 0.0;
         if (diff!=0.0 && Xmatrix[i][j] != 0.0) scaledDiff = 
           diff/Xmatrix[i][j];

	 if (diff==0.0 && scaledDiff==0.0) continue;

         sprintf(txt ,"%4d %4d  ",i,j);
         (*outStream) << txt;
         sprintf(txt ,"%20.13e  ",Xmatrix[i][j]);
         (*outStream) << txt;
         sprintf(txt ,"%4d %4d  ",Crow,Ccol);
         (*outStream) << txt;
         sprintf(txt ,"%20.13e  ",Smatrix[Crow][Ccol]);
         (*outStream) << txt;
         sprintf(txt ,"%20.13e  ",diff);
         (*outStream) << txt;
         sprintf(txt ,"%20.13e",scaledDiff);
         (*outStream) << txt;

         sprintf(txt ,"\t(%15s,%15s)\n",xyceNameMap[i].c_str(), xyceNameMap[j].c_str());
         (*outStream) << txt;

         if(fabs(diff) > diffmax) 
         {
           diffmax = fabs(diff);
           row = i;
           col = j;
         }

         if(fabs(scaledDiff) > scaledDiffMax)
         {
           scaledDiffMax = fabs(scaledDiff);
           srow = i;
           scol = j;
         }

         if (  (Xmatrix[i][j] > 0 && Smatrix[Crow][Ccol] < 0) ||
               (Xmatrix[i][j] < 0 && Smatrix[Crow][Ccol] > 0) )
         {
           signErrors[i][j] = 1;
           sprintf(txt ,"%4d %4d  ",i,j);
           (*outSt2) << txt;
           sprintf(txt ,"%20.13e  ",Xmatrix[i][j]);
           (*outSt2) << txt;
           sprintf(txt ,"%4d %4d  ",Crow,Ccol);
           (*outSt2) << txt;
           sprintf(txt ,"%20.13e  ",Smatrix[Crow][Ccol]);
           (*outSt2) << txt;
           sprintf(txt ,"%20.13e  ",diff);
           (*outSt2) << txt;
           sprintf(txt ,"%20.13e",scaledDiff);
           (*outSt2) << txt;

           sprintf(txt ,"\t(%15s,%15s)\n",xyceNameMap[i].c_str(), xyceNameMap[j].c_str());
           (*outSt2) << txt;
         }
       }
    }
  }
  mC = 0;
  mX = 0;
  for (i=0;i<=rowmax;i++)
  {
    for (j=0;j<=rowmax;j++)
    {
      Crow = XtoC[i];
      Ccol = XtoC[j];

      if (Sfill[Crow][Ccol] && !Xfill[i][j])
        mC++;
      if (Xfill[i][j] && !Sfill[Crow][Ccol])
        mX++;
    }
  }
  if (mC > 0) {
#ifdef BOTH_FILES_ARE_XYCE 
    (*outStream) << "The following Elements are specified only in the first matrix:" << endl;
#else
    (*outStream) << "The following Elements are specified only in the spice matrix:" << endl;
#endif
    for (i=0;i<=rowmax;i++)
    {
      for (j=0;j<=rowmax;j++)
      {
        Crow = XtoC[i];
        Ccol = XtoC[j];

        if (Sfill[Crow][Ccol] && !Xfill[i][j])
        {
          sprintf(txt ,"%20.13e  ",Smatrix[Crow][Ccol]);
          (*outStream) << "("<<Crow<<","<<Ccol<<")\t"
          << xyceNameMap[i].c_str() << ", " << xyceNameMap[j].c_str() 
          << " = " << txt << endl;
        }
      }
    }
    cout << endl;
  }
  if (mX > 0) {
#ifdef BOTH_FILES_ARE_XYCE 
    (*outStream) << "The following Elements are specified only in the second matrix:" << endl;
#else
    (*outStream) << "The following Elements are specified only in the Xyce matrix:" << endl;
#endif
    for (i=0;i<=rowmax;i++)
    {
      for (j=0;j<=rowmax;j++)
      {
        Crow = XtoC[i];
        Ccol = XtoC[j];

        if (!Sfill[Crow][Ccol] && Xfill[i][j])
        {
          sprintf(txt ,"%20.13e  ",Xmatrix[i][j]);
          (*outStream) << "("<<i<<","<<j<<")\t"
            <<xyceNameMap[i].c_str() << ", " << xyceNameMap[j].c_str() 
            << " = " << txt << endl;
        }
      }
    }
    cout << endl;
  }

   delete outStream;
   delete outSt2;

   sprintf(txt,"\nThe maximum difference is: %15.7e\n", diffmax);
#ifdef BOTH_FILES_ARE_XYCE
   cout << txt;
   sprintf(txt,"Xyce gold matrix entry: %4d %4d %15.7e\n",
                row,col,Xmatrix[row][col]);
   cout << txt;
   sprintf(txt,"Xyce test matrix entry: %4d %4d %15.7e\n\n",
                XtoC[row],
                XtoC[col],
                Smatrix[XtoC[row]][XtoC[col]]);
#else
   cout << txt;
   sprintf(txt,"Xyce matrix entry: %4d %4d %15.7e\n",
                row,col,Xmatrix[row][col]);
   cout << txt;
   sprintf(txt,"3f5  matrix entry: %4d %4d %15.7e\n\n",
                XtoC[row],
                XtoC[col],
                Smatrix[XtoC[row]][XtoC[col]]);
#endif
   cout << txt;
   sprintf(txt ,"The named indices of this entry: (%15s,%15s)\n",xyceNameMap[row].c_str(), xyceNameMap[col].c_str());
   cout << txt;


   sprintf(txt,"\nThe maximum scaled difference is: %15.7e\n", scaledDiffMax);
#ifdef BOTH_FILES_ARE_XYCE
   cout << txt;
   sprintf(txt,"Xyce gold matrix entry: %4d %4d %15.7e\n",
                srow,scol,Xmatrix[srow][scol]);
   cout << txt;
   sprintf(txt,"Xyce test matrix entry: %4d %4d %15.7e\n\n",
                XtoC[srow],
                XtoC[scol],
                Smatrix[XtoC[srow]][XtoC[scol]]);
#else
   cout << txt;
   sprintf(txt,"Xyce matrix entry: %4d %4d %15.7e\n",
                srow,scol,Xmatrix[srow][scol]);
   cout << txt;
   sprintf(txt,"3f5  matrix entry: %4d %4d %15.7e\n\n",
                XtoC[srow],
                XtoC[scol],
                Smatrix[XtoC[srow]][XtoC[scol]]);
#endif
   cout << txt;
   sprintf(txt ,"The named indices of this entry: (%15s,%15s)\n",xyceNameMap[srow].c_str(), xyceNameMap[scol].c_str());
   cout << txt;
}

