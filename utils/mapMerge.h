//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: mapMerge.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 3/26/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4 $
//
// Revision Date  : $Date: 2010/04/26 23:19:13 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------


#ifndef  _MAPMERGE_H
#define  _MAPMERGE_H

#include <string>

#include "basecmp.h"

using namespace std;

//-----------------------------------------------------------------------------
// Class         : ExtendedString
// Purpose       :
//                  This  class  extends the C++ string class by adding
//                  functions toUpper, toLower, and removeWhiteSpace.
//
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class ExtendedString : public string
{
  public:
    ExtendedString        (const char *X) : string(X) {};
    ExtendedString        (const string X) : string(X) {};
    ExtendedString & toUpper();
    ExtendedString & toLower();
    ExtendedString & removeWhiteSpace();
};


//-----------------------------------------------------------------------------
// Function      : ExtendedString::toUpper
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline ExtendedString & ExtendedString::toUpper(void)
{
  for( string::size_type i = 0; i < length(); i++ )
    (*this)[i] = toupper( (*this)[i] );

  return (*this);
}


//-----------------------------------------------------------------------------
// Function      : ExtendedString::toLower
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline ExtendedString & ExtendedString::toLower(void)
{
  for( string::size_type i = 0; i < length(); i++ )
    (*this)[i] = tolower( (*this)[i] );

  return (*this);
}

//-----------------------------------------------------------------------------
// Function      : ExtendedString::removeWhiteSpace
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline ExtendedString & ExtendedString::removeWhiteSpace(void)
{
  string::size_type N = string::npos - 1;

  while (N != string::npos)
  {
    N = find(" ");
    if (N != string::npos) erase (N, 1);
  }

  return (*this);
}

void readSpiceNamesFile (ifstream * inStream1, map <int, string> & chileMap1, map <string, int> & chileMap2);
void readXyceNamesFile  (ifstream * inStream2, map <int, string> & XyceMap1, map <string, int> & XyceMap2);



#endif
