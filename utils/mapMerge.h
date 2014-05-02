//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, 2013, Sandia Corporation, Albuquerque, NM.
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
// Revision Number: $Revision: 1.5 $
//
// Revision Date  : $Date: 2013/09/18 20:27:38 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------


#ifndef  _MAPMERGE_H
#define  _MAPMERGE_H

#include <string>

#include "basecmp.h"

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
class ExtendedString : public std::string
{
  public:
    ExtendedString        (const char *X) : string(X) {};
    ExtendedString        (const std::string &X) : string(X) {};
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
  for( std::string::size_type i = 0; i < length(); i++ )
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
  for( std::string::size_type i = 0; i < length(); i++ )
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
  std::string::size_type N = std::string::npos - 1;

  while (N != std::string::npos)
  {
    N = find(" ");
    if (N != std::string::npos) erase (N, 1);
  }

  return (*this);
}

void readSpiceNamesFile (ifstream * inStream1, std::map<int, std::string> & chileMap1, std::map<std::string, int> & chileMap2);
void readXyceNamesFile  (ifstream * inStream2, std::map<int, std::string> & XyceMap1, std::map<std::string, int> & XyceMap2);

#endif
