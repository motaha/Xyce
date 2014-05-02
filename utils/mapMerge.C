//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, 2013, Sandia Corporation, Albuquerque, NM.
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
// Revision Number: $Revision: 1.15 $
//
// Revision Date  : $Date: 2013/09/18 22:32:30 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <vector>
#include <string>

#include "mapMerge.h"

using namespace std;

void readSpiceNamesFile (std::ifstream * inStream1, std::map<int, std::string> & chileMap1, std::map<std::string, int> & chileMap2)
{
  int index;
  std::string name;
  bool endOfFile = (*inStream1).eof();

  while ( !endOfFile )
  {
    (*inStream1) >> index >> name;
    endOfFile = (*inStream1).eof();
    if (endOfFile) continue;

    // fix name to remove "#" characters, and replace with "_"
    std::string::iterator iter;
    std::string::iterator iter2;
    std::string::iterator first = name.begin ();
    std::string::iterator last  = name.end   ();
    std::string::iterator test1("#");
    std::string::iterator test2("_");

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

void readXyceNamesFile (std::ifstream * inStream2, std::map<int, std::string> & XyceMap1, std::map<std::string, int> & XyceMap2)
{
  int index;
  std::string name;
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
    std::string::iterator iter;
    std::string::iterator iter2;
    std::string::iterator first = name.begin ();
    std::string::iterator last  = name.end   ();
    std::string::iterator test1("x");
    std::string::iterator test2(":");
    std::string::iterator test3(":x");

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
  std::map<std::string, int> chileMap2;
  std::map<int, std::string> chileMap1;

  std::map<std::string, int> XyceMap2;
  std::map<int, std::string> XyceMap1;

  std::map<int, int> mergedMap;

  std::ifstream * inStream1;
  std::ifstream * inStream2;

  std::ofstream * outStream;
  std::ofstream * outStream2;

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
  cout << "Reading in the first Xyce name file: " << std::string(cargs[1]) <<endl;
#else
  cout << "Reading in the chilespice name file: " << std::string(cargs[1]) <<endl;
#endif

  inStream1 = new std::ifstream( cargs[1] );
#ifdef BOTH_FILES_ARE_XYCE 
  readXyceNamesFile (inStream1, chileMap1, chileMap2);
#else
  readSpiceNamesFile (inStream1, chileMap1, chileMap2);
#endif
  delete inStream1;

#ifdef BOTH_FILES_ARE_XYCE 
  cout << "Reading in the second Xyce name file: " << std::string(cargs[2]) <<endl;
#else
  cout << "Reading in the Xyce name file: " << std::string(cargs[2]) <<endl;
#endif

  inStream2 = new std::ifstream( cargs[2] );
  readXyceNamesFile (inStream2, XyceMap1, XyceMap2);
#if 0

  cout << "XyceMap2: " << endl;

  std::map<std::string, int>::iterator iterX2 = XyceMap2.begin();
  std::map<std::string, int>::iterator endX2  = XyceMap2.end ();
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

  outStream = new std::ofstream( "mergedMap.txt" );
  outStream2 = new std::ofstream( "namePerm.txt" );

  std::map<int, std::string>::iterator iter = chileMap1.begin();
  std::map<int, std::string>::iterator end  = chileMap1.end ();


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

  std::map<int,int>::iterator iterMM = mergedMap.begin();
  std::map<int,int>::iterator endMM  = mergedMap.end ();

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


