//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: mapMerge.C,v $
//
// Purpose        : This program is intended to take two name maps - 
//                  one from chilespice, one from Xyce, and create 
//                  a single mapping going from chilespice solution 
//                  indices to Xyce solution indices.
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
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2011/03/15 19:52:22 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <vector>
#include <string>

#include "mapMerge.h"

using namespace std;

void readSpiceNamesFile (ifstream * inStream1, map <int, string> & chileMap1, map <string, int> & chileMap2)
{
  int index;
  string name;
  bool endOfFile = (*inStream1).eof();

  while ( !endOfFile )
  {
    (*inStream1) >> index >> name;
    endOfFile = (*inStream1).eof();
    if (endOfFile) continue;

    // fix name to remove "#" characters, and replace with "_"
    string::iterator iter;
    string::iterator iter2;
    string::iterator first = name.begin ();
    string::iterator last  = name.end   ();
    string::iterator test1("#");
    string::iterator test2("_");

    // are there any ":x" characters?  If so convert to ":"
    first = name.begin ();
    last  = name.end ();
    iter  = first;

    // convert all instances of ":x" to ":".
    int icount;
    for (icount=0,iter=first; iter != last; icount++, iter++)
    {
      if (*iter == *test1)
      {
        *iter = *test2;
      }
    }

    ExtendedString tmpStr(name);
    tmpStr.toLower();
    //cout << "tmpStr = " << tmpStr << endl;

    chileMap1[index] = tmpStr;
    chileMap2[tmpStr] = index;

    endOfFile = (*inStream1).eof();
  }
}

void readXyceNamesFile (ifstream * inStream2, map <int, string> & XyceMap1, map <string, int> & XyceMap2)
{
  int index;
  string name;
  bool endOfFile = (*inStream2).eof();

  (*inStream2) >> name; // this is the HEADER

  while ( !endOfFile )
  {
    (*inStream2) >> index >> name;
    endOfFile = (*inStream2).eof();
    if (endOfFile) continue;

    // save the unmodified name in XyceMap1.
    if (XyceMap1.find( index ) == XyceMap1.end())
    {
      XyceMap1[index] = name;
    }
    else
    {
      cout << "Warning!  skipping duplicate index: " << index << "\t" 
        << name << "\t" << XyceMap1[index] << endl;
    }

#ifndef BOTH_FILES_ARE_XYCE 
    // Modify the name to chilespice format
    string::iterator iter;
    string::iterator iter2;
    string::iterator first = name.begin ();
    string::iterator last  = name.end   ();
    string::iterator test1("x");
    string::iterator test2(":");
    string::iterator test3(":x");

    // first look for subcircuit names, and remove
    // the "x" prefix.

    // Is the first charcter an "x"?
    if (*first == *test1) { name.erase(first); }

    // are there any ":x" characters?  If so convert to ":"
    first = name.begin ();
    last  = name.end ();
    iter  = first;

    // convert all instances of ":x" to ":".
    int icount;
    for (icount=0,iter=first; iter != last; icount++, iter++)
    {
      if (icount != 0)
      {
        iter2 = iter;  iter2++;
        if (iter2 != last)
        {
                    //    :                   x
          if (*iter == *test2 && *iter2 == *test1)
          {
            iter++;
            name.erase(iter);
          }
        }
      }
    }
#endif

    // Save the modified name in XyceMap2. 
    XyceMap2[name] = index;
    endOfFile = (*inStream2).eof();
  }

}

int main (int iargs , char *cargs[] )
{
  map <string, int> chileMap2;
  map <int, string> chileMap1;

  map <string, int> XyceMap2;
  map <int, string> XyceMap1;

  map <int, int> mergedMap;

  ifstream * inStream1;
  ifstream * inStream2;

  ofstream * outStream;
  ofstream * outStream2;

  if (iargs != 3)
  {
#ifdef BOTH_FILES_ARE_XYCE 
     cout << "Usage:  mapMerge Xyce_name_file1 Xyce_name_file2 <return>\n";
#else
     cout << "Usage:  mapMerge chilspice_name_file Xyce_name_file <return>\n";
#endif
     exit(0);
  }


#ifdef BOTH_FILES_ARE_XYCE 
  cout << "Reading in the first Xyce name file: " << string(cargs[1]) <<endl;
#else
  cout << "Reading in the chilespice name file: " << string(cargs[1]) <<endl;
#endif

  inStream1 = new ifstream( cargs[1] );
#ifdef BOTH_FILES_ARE_XYCE 
  readXyceNamesFile (inStream1, chileMap1, chileMap2);
#else
  readSpiceNamesFile (inStream1, chileMap1, chileMap2);
#endif
  delete inStream1;

#ifdef BOTH_FILES_ARE_XYCE 
  cout << "Reading in the second Xyce name file: " << string(cargs[2]) <<endl;
#else
  cout << "Reading in the Xyce name file: " << string(cargs[2]) <<endl;
#endif

  inStream2 = new ifstream( cargs[2] );
  readXyceNamesFile (inStream2, XyceMap1, XyceMap2);
#if 0

  cout << "XyceMap2: " << endl;

  map<string, int>::iterator iterX2 = XyceMap2.begin();
  map<string, int>::iterator endX2  = XyceMap2.end ();
  for( ;iterX2!=endX2;++iterX2)
  {
    cout << iterX2->first << "\t" << iterX2->second << endl;
  }

#endif
  delete inStream2;

  int chileSize1 = chileMap1.size();
  int chileSize2 = chileMap2.size();
  int xyceSize1  =  XyceMap1.size();
  int xyceSize2  =  XyceMap2.size();

  cout << " chile size 1 = " << chileSize1 << endl;
  cout << " chile size 2 = " << chileSize2 << endl;
  cout << "  xyce size 1 = " << xyceSize1 << endl;
  cout << "  xyce size 2 = " << xyceSize2 << endl;

  outStream = new ofstream( "mergedMap.txt" );
  outStream2 = new ofstream( "namePerm.txt" );

  map<int, string>::iterator iter = chileMap1.begin();
  map<int, string>::iterator end  = chileMap1.end ();


  for ( ; iter != end; iter++)
  {
    if (XyceMap2.find( iter->second ) != XyceMap2.end())
    {
       mergedMap[iter->first] = XyceMap2[iter->second];
    }
    else
    {
      cout << "Could not find variable name: " << iter->second;
#ifdef BOTH_FILES_ARE_XYCE 
      cout << " in the second Xyce map\n";
#else
      cout << " in the Xyce map\n";
#endif
      mergedMap[iter->first] = -999999;
    }
  }

  map <int,int>::iterator iterMM = mergedMap.begin();
  map <int,int>::iterator endMM  = mergedMap.end ();

  for ( ; iterMM != endMM; iterMM++)
  {
    (*outStream) << "\t" << iterMM->first;
    (*outStream) << "\t" << iterMM->second;
    if (iterMM->second != -999999)
    {
      (*outStream) << "\t" << XyceMap1[iterMM->second];
    }
    else
    {
      (*outStream) << "\t" << chileMap1[iterMM->first];
    }
    (*outStream) << endl;


    (*outStream2) << "\t" << iterMM->first;
    (*outStream2) << "\t" << iterMM->second;
    (*outStream2) << endl;
  }

  delete outStream;
  delete outStream2;

}


