//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Revision Number: $Revision: 1.2.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef N_DEV_ReactionLexer_H
#define N_DEV_ReactionLexer_H
#include <iostream>
#include <string>
#include <map>

class N_DEV_ReactionLexer: public yyFlexLexer 
{

public:
  N_DEV_ReactionLexer(istream *input = 0, ostream* output = 0): yyFlexLexer(input,output) {};

  virtual ~N_DEV_ReactionLexer() {};

  int getToken(N_DEV::N_DEV_ReactionParser::semantic_type *lvalp, N_DEV::location *llocp,map<string,int> &theSpecies);
};
#endif
