//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : N_IO_SpiceSeparatedFieldTool.H
//
// Purpose       : This file defines the N_IO_SpiceSeparatedFieldTool
//                 class. (for use in parsing the netlist)
//
// Special Notes :
//
// Creator       : Alan Lundin, SNL
//
// Creation Date : 08/31/99
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.20.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef SpiceSeparatedFieldTool_H
#define SpiceSeparatedFieldTool_H

// ---------- Standard Includes ----------
#ifdef NEED_BOOL_H
#include <bool.h>
#endif

#include <string>
#include <vector>
#include <iosfwd>

// ----------   Xyce Includes   ----------
#include <N_UTL_Packable.h>
#include <N_UTL_Xyce.h>

class N_IO_SpiceSeparatedFieldTool
{

public:

  struct StringToken : public Packable
  {
     StringToken(): lineNumber_(0), string_(""){};
     StringToken(string const& fileName, size_t lineNum, string const& str):
        lineNumber_(lineNum), string_(str){};
     ~StringToken(){};

     string string_;
     size_t lineNumber_;

     // Packing functionality.
     Packable * instance() const;

     // Counts bytes needed to pack block.
     int packedByteCount() const;

     // Packs OptionBlock into char buffer using MPI_PACK.
     void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;

     // Unpacks OptionBlock from char buffer using MPI_UNPACK.
     void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);
  };

  N_IO_SpiceSeparatedFieldTool(ifstream & input, string const & fileName, 
      const vector< pair< string, string > > & externalParams );
  // Constructor

  ~N_IO_SpiceSeparatedFieldTool() { }
  // Destructor

  int getLine(vector<StringToken> & line);
  // R int
  // R- Returns 1 if end-of-file has not been reached, 0 otherwise.
  // O str
  // O- A string containing the next line of input.
  // Read a line and split it into fields that are stored in the
  // the private attribute "Fields". The line is stored in the
  // private member attribute "Line".
  // KRS, 10/11/07:  this is the original version of getLine, modified to 
  // handle ground synonym replacements.

  int getLine(vector<StringToken> & line, bool replgndvar);
  // R int
  // R- Returns 1 if end-of-file has not been reached, 0 otherwise.
  // O str
  // O- A string containing the next line of input.
  // Read a line and split it into fields that are stored in the
  // the private attribute "Fields". The line is stored in the
  // private member attribute "Line".
  //  KRS, 10/11/07:  this new version of getLine has been added specifically
  // for the preprocess phase so as to ensure that no changes to the netlist
  // file are made during preprocessing (i.e., if '.PREPROCESS REPLACEGND 
  // true' is in the netlist file halfway through the netlist, we don't want
  // occurences of 'GND' to be replaced by '0' halfway through the netlist.
  // Technically, doing so shouldn't cause errors, but this is more of a 
  // precaution for unforseen circumstances.

  int getLineWithComments(vector<StringToken> & line);
  // R int
  // R- Returns 1 if end-of-file has not been reached, 0 otherwise.
  // O str
  // O- A string containing the next line of input.
  // Read a line and split it into fields that are stored in the
  // the private attribute "Fields". The line is stored in the
  // private member attribute "Line".
  // KRS, 12/05/07:  this is the original version of getLine, except that it
  // does NOT skip comments and blank lines (used for netlist copying).
  
  int getLine2(vector<StringToken> & line);
  // R int
  // R- Returns 1 if end-of-file has not been reached, 0 otherwise.
  // O str
  // O- A string containing the next line of input.
  // Read a line and split it into fields that are stored in the
  // the private attribute "Fields". The line is stored in the
  // private member attribute "Line". This function differs from
  // getLine in that it ignores "=" signs on the input line.
  
  void changeCursorLineNumber(int token);
  // Increase or decrease the line number of cursor by the given amount.
  // Increase if token is positive; decrease if token is negative.
  // If the result cursor line number is less than 1, set it to 1.

  // accessors
  void setLineNumber( int loc ) 
  { 
    ( loc <= 0 ) ? cursorLineNum_ = 1 : cursorLineNum_ = loc;
  }
  int getLineNumber() { return cursorLineNum_; }

  // Return the current character position in file.
  streampos getFilePosition() const;

  // Set the location in the ifstream at which the next input operation
  // will begin.
  bool setLocation(streampos const& startLocation);

  const string & getFileName() const {return fileName_;};

protected:

private :

  ifstream & in_;

  string fileName_;
  size_t cursorLineNum_; //The physical line number of the cursor
  vector< pair< string, string > > externalParams_;

  // R bool
  // R-
  bool NextChar_(char & c);
  void substituteExternalParams(vector<StringToken>& line);
  void skipToEndOfLine_();
  void skipCommentsAndBlankLines_();

public:

};

#endif // SpiceSeparatedFieldTool_H 
