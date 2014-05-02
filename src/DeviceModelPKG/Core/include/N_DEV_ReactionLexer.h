//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_ReactionLexer.h,v $
//
// Purpose : Declaration of the ReactionLexer class for the reaction
//           parsing package
//
// Special Notes  : 
//
// Creator        : Thomas V. Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 08/10/2006
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef N_DEV_ReactionLexer_H
#define N_DEV_ReactionLexer_H

#include <iostream>
#include <string>
#include <map>

namespace Xyce {
namespace Device {

class ReactionLexer: public yyFlexLexer 
{
public:
  ReactionLexer(std::istream *input = 0, std::ostream* output = 0)
    : yyFlexLexer(input,output)
  {}
    
  virtual ~ReactionLexer()
  {}
    
  int getToken(XyceDevice::ReactionParser::semantic_type *lvalp, XyceDevice::location *llocp, std::map<std::string,int> &theSpecies);
};

}
}

#endif
