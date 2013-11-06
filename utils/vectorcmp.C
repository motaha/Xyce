//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: vectorcmp.C,v $
//
// Purpose        : This program compares two vector files, one from
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
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2010/04/29 20:38:11 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

// ----------   Standard Includes   ----------
#include <iostream>
#include <fstream>
#include <string>
#include <map>
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

using namespace std;

#define LIMIT 1.0e-25


#ifndef RHSCMP_BUILD
//-----------------------------------------------------------------------------
// Function      : main
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/06/02
//-----------------------------------------------------------------------------
int main (int iargs, char *cargs[])
{
  vectorCompare VC;
  VC.run(iargs,cargs);
}
#endif

//-----------------------------------------------------------------------------
// Function      : vectorCompare::vectorCompare
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/06/02
//-----------------------------------------------------------------------------
vectorCompare::vectorCompare ()
  : baseCompare()
{

}

//-----------------------------------------------------------------------------
// Function      : vectorCompare::~vectorCompare
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
vectorCompare::~vectorCompare ()
{
  Xvector.clear();
  Svector.clear();
}

//-----------------------------------------------------------------------------
// Function      : vectorCompare::doAllocations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/06/02
//-----------------------------------------------------------------------------
void vectorCompare::doAllocations ()
{
  cout << "Allocating arrays... " << endl;
  Svector.resize(isize+10,0.0);
  Xvector.resize(isize+10,0.0);
  cout << "  ...done." << endl;
}

//-----------------------------------------------------------------------------
// Function      : vectorCompare::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
void vectorCompare::run (int iargs, char *cargs[])
{

  // first check that there are 2 arguments:
  if (iargs != 3)
  {
#ifdef BOTH_FILES_ARE_XYCE
    cout << "Usage:  veccmp Xyce_vector_gold Xyce_vector_test <return>\n";     exit(0);
#else
    cout << "Usage:  veccmp spice_vector_file Xyce_vector_file <return>\n";     exit(0);
#endif
  }

  readMergedMap ();
  doAllocations ();

  // read in the first vector file:    
  if (XyceFileFlag)
  {
    readXyceVectorFile (cargs[1], Svector);
  }
  else
  {
    readSpiceVectorFile (cargs[1], Svector);
  }


  // read in the second vector file:  (by default a Xyce file)
  XyceFileFlag = true;
  if (XyceFileFlag)
  {
    readXyceVectorFile (cargs[2], Xvector);
  }
  else
  {
    readSpiceVectorFile (cargs[2], Xvector);
  }

  if (rowmax==0) rowmax=isize;
  outputResult ();
}

//-----------------------------------------------------------------------------
// Function      : vectorCompare::readXyceVectorFile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/06/02
//-----------------------------------------------------------------------------
void vectorCompare::readXyceVectorFile (char *carg, vector<double> & Xvector1)
{
  // open and read xyce file.
  cout << "Reading the Xyce vector file: " << std::string(carg) << endl;

  ifstream * inStream = new ifstream(carg);

  int imax=0;

  int size;
  (*inStream) >> size; // this is the size, as specified by the file.

  do
  {
    int i,j;
    double value;
    (*inStream) >> i;
    (*inStream) >> j;
    (*inStream) >> value;

    //cout <<"\t"<<imax<<"\t"<<value<<endl;

    if ( (*inStream).eof() ) break;

    Xvector1[imax] = value;
    imax++;

  } while (!((*inStream).eof()));

  delete inStream;

  cout << " ...done" << endl;
}

//-----------------------------------------------------------------------------
// Function      : vectorCompare::readSpiceVectorFile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/13/02
//-----------------------------------------------------------------------------
void vectorCompare::readSpiceVectorFile (char *carg, vector<double> & Svector1)
{
  cout << "Reading the Spice vector file: " << std::string(carg) << endl;
  //cout << "Reading the Spice vector file...  " << endl;

  ifstream * inStream = new ifstream(carg);

  int rowmax = 0;

  do
  {
    (*inStream) >> val;
    if( (*inStream).eof() ) break;

    Svector1[rowmax] = val;

    rowmax++;

    //cout <<"\t"<<rowmax<<"\t"<<val<<endl;

  } while ( !((*inStream).eof()));

  delete inStream;
  rowmax -= 1;
  cout << " ...done" << endl;
}

//-----------------------------------------------------------------------------
// Function      : vectorCompare::outputResult
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/06/02
//-----------------------------------------------------------------------------
void vectorCompare::outputResult ()
{
  // open result file
  ofstream * outStream = new ofstream("vecresult.txt");
  
  double diffmax = 0.0;
  double scaledDiffMax = 0.0;
  int srow = 0;

  char txt[256];  for(i=0;i<256;i++) txt[i] = 0;

  for (i=0;i<=rowmax;i++)
  {
#ifdef BOTH_FILES_ARE_XYCE
     int Crow = i;
#ifdef FILES_HAVE_DIFFERENT_ORDER 
     Crow = XtoC[i];
#endif
#else
     int Crow = XtoC[i]-1;
#endif

     double diff = Xvector[i] - Svector[Crow];

     double scaledDiff = 0.0;
     if (diff!=0.0 && Xvector[i] != 0.0)
       scaledDiff = diff/Xvector[i];

     sprintf(txt,"%4d ",i);
     (*outStream) << txt;
     sprintf(txt,"%25.18e ",Xvector[i]);
     (*outStream) << txt;
     sprintf(txt,"%4d ",Crow);
     (*outStream) << txt;
     sprintf(txt,"%25.18e ",Svector[Crow]);
     (*outStream) << txt;
     sprintf(txt,"%25.18e ",diff);
     (*outStream) << txt;
     sprintf(txt,"%25.18e",scaledDiff);
     (*outStream) << txt;
     (*outStream) << "\t" << xyceNameMap[i];
     (*outStream) << endl;

     if (fabs(diff) > diffmax) 
     {
       diffmax = fabs(diff);
       row = i;
     }

     if (fabs(scaledDiff) > scaledDiffMax) 
     {
       scaledDiffMax = fabs(scaledDiff);
       srow = i;
     }

   }

  delete outStream;

  sprintf(txt,"\nThe maximum difference is: %15.7e\n", diffmax);
#ifdef BOTH_FILES_ARE_XYCE
  cout << txt;
  sprintf(txt,"Xyce test vector entry: %4d %15.7e\n",
               row,Xvector[row]);
  cout << txt;

  int row1 = XtoC[row];
  double val1 = Svector[row];
#ifdef FILES_HAVE_DIFFERENT_ORDER 
  val1 = Svector[XtoC[row]];
#endif

  sprintf(txt,"Xyce gold vector entry: %4d %15.7e\n\n",row1,val1);

#else
  cout << txt;
  sprintf(txt,"Xyce  vector entry: %4d %15.7e\n",
               row,Xvector[row]);
  cout << txt;
  sprintf(txt,"spice vector entry: %4d %15.7e\n\n",
               XtoC[row],
               Svector[XtoC[row]-1]);
#endif
  cout << txt;
  sprintf(txt,"name for this entry: %15s\n",xyceNameMap[row].c_str());
  cout << txt;



  sprintf(txt,"\nThe maximum scaled difference is: %15.7e\n", scaledDiffMax);
#ifdef BOTH_FILES_ARE_XYCE
  cout << txt;
  sprintf(txt,"Xyce test vector entry: %4d %15.7e\n", srow,Xvector[srow]);
  cout << txt;

  int row2 = XtoC[srow];
  double val2 = Svector[srow];
#ifdef FILES_HAVE_DIFFERENT_ORDER 
  val2 = Svector[XtoC[srow]];
#endif

  sprintf(txt,"Xyce gold vector entry: %4d %15.7e\n\n",row2,val2);

#else
  cout << txt;
  sprintf(txt,"Xyce  vector entry: %4d %15.7e\n", srow,Xvector[srow]);
  cout << txt;
  sprintf(txt,"spice vector entry: %4d %15.7e\n\n",
               XtoC[srow],
               Svector[XtoC[srow]-1]);
#endif
  cout << txt;
  sprintf(txt,"name for this entry: %15s\n",xyceNameMap[srow].c_str());
  cout << txt << endl;

}

