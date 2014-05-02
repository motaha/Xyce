//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : N_IO_SpiceSeparatedFieldTool.C
//
// Purpose       : This file defines the functions contained in the
//                 N_IO_SpiceSeparatedFieldTool class.
//
// Special Notes :
//
//
// Creator       : Alan Lundin
//
// Creation Date : 08/31/99
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.71 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <fstream>

// ----------   Xyce Includes   ----------
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_PDS_Comm.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_fwd.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Misc.h>

namespace Xyce {
namespace IO {

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::SpiceSeparatedFieldTool
// Purpose       : constructor
// Special Notes :
// Creator       : Jian Li
// Creation Date : 07/12/2002
//-------------------------------------------------------------------------
SpiceSeparatedFieldTool::SpiceSeparatedFieldTool(std::ifstream & input,
  std::string const & fileStr, const std::vector< std::pair< std::string, std::string > > & externalParams)
: in_(input), cursorLineNum_(1), fileName_(fileStr), externalParams_(externalParams)
{
}

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::getLine
// Purpose       : Read the line (one character at a time via calling
//                 NextChar_),  and split the line into fields (via calling
//                 SplitLine).  KRS 10/11/07:  I have overloaded getLine to
//                 accept a boolean input to do ground synonym checking.
//                 This version is the original version of getLine that does
//                 NOT do gnd synonym checking (no boolean input).
// Special Notes :
// Scope         : public
// Creator       : Alan Lundin
// Creation Date :
//-------------------------------------------------------------------------
int SpiceSeparatedFieldTool::getLine(std::vector<StringToken>& line)
{
  char c=0;
  const std::string nonid(" \t\n\r(){},='");
  line.clear();
  
  skipCommentsAndBlankLines_();

  bool endOfLine = false; //is the line end?

  //If either NextChar_(c) == false or endOfLine == true
  //then the line is end. But checking endOfLine need to be
  //in front of NextChar_ to avoid reading next line.
  while ( endOfLine == false && NextChar_(c)  )
  {
    StringToken field;
    field.lineNumber_ = cursorLineNum_;
    field.string_.reserve(16);

    switch (c)
    {
      //Skip whitespaces
      case ' ' :
      case '\t':
      case '\r':
      case '\n':
         break;

      //case '*' : //should never be here
         //cout << fileName_ <<":" <<cursorLineNum_ <<": "
           //<< "Unexpected char '*', ignore it.\n";
         //break;

      case '=' :
        field.string_ = "=";
        break;

      case ':' :
        field.string_ = ":";
        break;

      case ',' :
        field.string_ = ",";
        break;

      case '(' :
        field.string_ = "(";
        break;

      case ')' :
        field.string_ = ")";
        break;

      case '"' :
      {
        field.string_ = '"';
        while (NextChar_(c))
        {
          field.string_ += c;
          if (c == '"') break;
        }

        break;
      }

      case '\'' :
      {
        // convert single-quote-enclosed expressions to braces-enclosed
        field.string_ = "{";

        bool scanning = true;
        while ( scanning && NextChar_(c) )
        {
          if ( c == '\'')
          {
            field.string_ += '}';
            scanning = false;         
          }
          else 
          {
            field.string_ += c;
          }
        }

        if ( scanning ) //then NextChar_(c) == false
        {
          endOfLine = true;
        }

        break;
      }

      case '{' :
      {
        // test {block} for expected VECTOR* contents 
        int lineSize = line.size();
        // avoid reading past beginning of line std::vector<> 
        if( lineSize > 2)
        {
          std::string earlierField ( ExtendedString( 
          line[ lineSize - 2 ].string_ ).toUpper() );
          // compare the prior field to VECTOR* tags
          if( earlierField == "DOPINGPROFILES" || earlierField == "SOURCELIST" ||
              earlierField == "LAYER" ||
              earlierField == "REGION" || earlierField == "NODE" ||
              earlierField == "RXNREGION" || earlierField == "COIL" || 
              earlierField.find( "MM_" ) == 0 )                          // earlierField starts with MM_ for membrane model parameter
          {
            // break here to process contents normally
            field.string_ = "{";
            break;
          }
          /*
          if( earlierField == "INDVARS" || earlierField == "INDFEQUS" || earlierField == "INDQEQUS" )
          {
            field.string_ = ",";
            break;
          }
          */

          // look for bad scientific notation 0.0e+{0} and 0.0e-{0}
          // not checking poor spacing & 0.0e{0}; sign must be present
          int lastLength = ( line[lineSize - 1].string_ ).size();
          if( lastLength > 2 )
          {
            std::string last2 ( line[lineSize - 1].string_, lastLength - 2, 2 );

            if( last2 == "e+" || last2 == "e-" || last2 == "E+" || last2 == "E-" )
            {
              Report::UserError0().at(fileName_, line[lineSize - 1].lineNumber_) << "Invalid notation encountered.";
            }
          } 
        }

        // otherwise this is not a VECTOR* type so
        // everything from left brace to matching
        // right brace is the next field
        field.string_ = "{";

        int unmatchedLeftBraces = 1;
        while ( unmatchedLeftBraces > 0 &&  NextChar_(c) )
        {
          field.string_ += c;
          if ( c == '}')
          {
            unmatchedLeftBraces--;
          }

          if (c == '{' )
          {
            ++unmatchedLeftBraces;
          }
        }

        if (unmatchedLeftBraces > 0) //then NextChar_(c) == false
        {
          endOfLine = true;
        }

        break;
      }

      // should only get here if VECTOR* {block} was found earlier
      case '}':                      
      {
        field.string_ = "}";  
        break;
      }

      default : 
      {
        //chars other than above
        // store the char until encount the a nonid char
        field.string_ += c; 
        bool metNonIdChar = false; //did it encounter nonid char?

        while (NextChar_(c))
        {
          if (nonid.find(c) == nonid.npos) //c is a id char
          {
            field.string_ += c; 
          }
          else //c is nonid, end of this field
          {
             metNonIdChar = true;
             if (c == '=' || c == ',' || c == '(' || c == ')' || 
                 c == '{' || c == '}' || c == '\'' )
               in_.putback(c);
             if (c == '\n') cursorLineNum_--;
             break;
          }
        }

        if (metNonIdChar == false) //then NextChar_(c) == false makes loop end
        {
           endOfLine = true;
        }
      }
    }
    if (field.string_.length() > 0)
    {
      line.push_back(field);
    }
  }
  // do any needed substitutions
  substituteExternalParams( line );
  
  return !in_.eof();
}


//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::getLine
// Purpose       : Read the line (one character at a time via calling
//                 NextChar_),  and split the line into fields (via calling
//                 SplitLine).  
//
// Special Notes :
//
//                 KRS 10/11/07:  this version of getLine 
//                 actually does ground synonym replacement (based upon the
//                 value of the boolean variable replgndvar).
//
//                 ERK 05/21/09: This function contains code that will
//                 replace single quotes with curly braces when surrounding
//                 expressions (bug 1358 on charleston).   
//
//                 It will also (to a more limited extent) find expressions 
//                 that have neither curly braces nor single quotes around 
//                 them, and add curly braces. (bug 1692 on charleston)
//
//                 For the no-braces, no-quotes expressions, this will
//                 only work under limited circumstances, since it is difficult
//                 to parse this.  Most "." lines don't support it, and 
//                 an "=" sign is required.  Also, the expression 
//                 needs to not contain any whitespace.
//
// Scope         : public
// Creator       : Alan Lundin
// Creation Date :
//-------------------------------------------------------------------------
int SpiceSeparatedFieldTool::getLine(std::vector<StringToken>& line, bool 
replgndvar)
{
  char c(0);
  const std::string nonid(" \t\n\r(){},='");
  line.clear();

  skipCommentsAndBlankLines_();

  bool endOfLine (false); //is the line end?

  //If either NextChar_(c) == false or endOfLine == true
  //then the line is end. But checking endOfLine need to be
  //in front of NextChar_ to avoid reading next line.
  bool noImplicitExpression(true), firstChar(true), isModel(false);

  char lastNonWhite(0);
  char saveC;
  std::string lastTok("");
  while ( endOfLine == false && NextChar_(c)  )
  {
    saveC = c;
    if (firstChar) 
    {
      firstChar = false;
      if (c == '.')
      {
        noImplicitExpression = true;
      }
      else
      {
        noImplicitExpression = false;
      }
    }
    StringToken field;
    std::string isGndSynonym(""); //used to check for "GND", "GND!", "GROUND"
                             //or stuff like "GrOUnD" and replace them with "0"
    field.lineNumber_ = cursorLineNum_;
    field.string_.reserve(16);
    if (c != ' ' && c != '\t' && c != '\r' && c != '\n') 
    {
      if (!noImplicitExpression && lastNonWhite == '=' &&
         ((!isModel && ((lastTok != "IC")&&(lastTok != "TC")) ) || (isModel && lastTok != "VERSION")) &&
          Util::possibleParam(lastTok) && c != '{' && c != '\'' && c != '\"') 
      {
        int numParen (0);
        if (c == '(')
        {
          numParen++;
        }
        std::string extraString(1,c);
        while (!endOfLine) 
        {
          endOfLine = !NextChar_(c);
          if (c == '(')
          {
            numParen++;
          }
          else if (c == ')')
          {
            numParen--;
          }
          if (endOfLine || numParen < 0 || c == ' ' || c == '\t' || c == '\r' || c == '\n')
          {
            break;
          }
          extraString.push_back(c);
        }
        if (numParen < 0)
        {
          in_.putback(c);
        }
        bool simple = Util::isInt(extraString) 
                   || Util::isValue(extraString) 
                   || Util::possibleParam(extraString) 
                   || Util::isBool(extraString) ;
        if (!simple)
        {
          field.string_ = "{";
        }
        else
        {
          field.string_ = "";
        }
        field.string_ += extraString;
        if (!simple)
        {
          field.string_ += "}";
        }
      }
      else 
      {
        switch (c)
        {
          case '=' :
             field.string_ = "=";
             break;

          case ':' :
             field.string_ = ":";
             break;

          case ',' :
              field.string_ = ",";
              break;

          case '(' :
              field.string_ = "(";
              break;

          case ')' :
              field.string_ = ")";
              break;

          case '"' :
           {
              field.string_ = '"';
              while (NextChar_(c))
              {
                field.string_ += c;
                if (c == '"') break;
              }

              break;
           }

          case '\'' :
          {
            // convert single-quote-enclosed expressions to braces-enclosed
            field.string_ = "{";

            bool scanning = true;
            while ( scanning && NextChar_(c) )
            {
              if ( c == '\'')
              {
                field.string_ += '}';
                scanning = false;         
              }
              else 
              {
                field.string_ += c;
              }
            }

            if ( scanning ) //then NextChar_(c) == false
            {
              endOfLine = true;
            }

            break;
          }

          case '{' :
          {
            // test {block} for expected VECTOR* contents 
            int lineSize = line.size();
            // avoid reading past beginning of line std::vector<> 
            if( lineSize > 2)
            {
              std::string earlierField ( line[ lineSize - 2 ].string_ ) ;
              Util::toUpper(earlierField) ;
              // compare the prior field to VECTOR* tags
              if( earlierField == "DOPINGPROFILES" || earlierField == "SOURCELIST" ||
                  earlierField == "LAYER" ||
                  earlierField == "REGION" || earlierField == "NODE" ||
                  earlierField == "RXNREGION" || earlierField == "COIL" || 
                  earlierField.find( "MM_" ) == 0 )                          // earlierField starts with MM_ for membrane model parameter
              {
                // break here to process contents normally
                field.string_ = "{";
                break;
              }

              // look for bad scientific notation 0.0e+{0} and 0.0e-{0}
              // not checking poor spacing & 0.0e{0}; sign must be present
              int lastLength = ( line[lineSize - 1].string_ ).size();
              if( lastLength > 2 )
              {
                std::string last2 ( line[lineSize - 1].string_, lastLength - 2, 2 );

                if( last2 == "e+" || last2 == "e-" || last2 == "E+" || last2 == "E-" )
                {
                  Report::UserError0().at(fileName_, line[lineSize - 1].lineNumber_) << "Invalid notation encountered.";
                }
              } 
            }

            // otherwise this is not a VECTOR* type so
            // everything from left brace to matching
            // right brace is the next field
            field.string_ = "{";

            int unmatchedLeftBraces = 1;
            while ( unmatchedLeftBraces > 0 &&  NextChar_(c) )
            {
              field.string_ += c;
              if ( c == '}')
              {
                unmatchedLeftBraces--;
              }

              if (c == '{' )
              {
                ++unmatchedLeftBraces;
              }
            }

            if (unmatchedLeftBraces > 0) //then NextChar_(c) == false
            {
              endOfLine = true;
            }

            break;
          }

          // should only get here if VECTOR* {block} was found earlier
          case '}':                      
          {
            field.string_ = "}";  
            break;
          }

          default : 
          {
            //chars other than above
            // store the char until encount the a nonid char
            field.string_ += c; 
            isGndSynonym += islower(c) ? toupper(c) : c; //convert to upper for ground synonym check.
            bool metNonIdChar = false; //did it encounter nonid char?

            while (NextChar_(c))
            {
              if (nonid.find(c) == nonid.npos) //c is a id char
              {
                field.string_ += c; 
                isGndSynonym += islower(c) ? toupper(c) : c; //ditto on the GND synonyms
              }
              else //c is nonid, end of this field
              {
                metNonIdChar = true;
                if (c == '=' || c == ',' || c == '(' || c == ')' || c == '{' || c == '}' || c == '\'' )
                {
                  in_.putback(c);
                }
                if (c == '\n') cursorLineNum_--;
                break;
              }
            }

            if (metNonIdChar == false) //then NextChar_(c) == false makes loop end
            {
              endOfLine = true;
            }
          }
        } 
      }
      lastNonWhite = saveC;
    }

    //Check for ground synonyms here.  If replace_ground_ flag is set to
    // true, and a synonym is found, replace field.string_ with
    //"0" and move on.
    if (replgndvar)
    {
      if (isGndSynonym=="GND" || isGndSynonym=="GND!" || isGndSynonym=="GROUND") 
      {
#ifdef Xyce_DEBUG_IO
        Xyce::dout() << "Node " <<field.string_ << " is a synonym for ground.  " <<  
        "Replacing node name with 0" << std::endl;
#endif

        field.string_ = "0";
      }
    }

    if ( !(field.string_.empty()) )
    {
      std::string ucFieldString(field.string_);
      Util::toUpper(ucFieldString);

      if (line.empty()) 
      {
        if (ucFieldString == "YPDE" || ucFieldString == "YEXT" || ucFieldString == "YXYGRA")
        {
          noImplicitExpression = true;
        }
        else if (ucFieldString == ".MODEL") 
        {
          //noImplicitExpression = false;
          noImplicitExpression = true;
          isModel = true;
        }
      }
      else 
      {
        if (ucFieldString == "TABLE")
        {
          noImplicitExpression = true;
        }
      }
      if (lastNonWhite != '=')
      {
        lastTok = ucFieldString;
      }
      line.push_back(field);
    }
  }
  // do any needed substitutions
  substituteExternalParams( line );
  return !in_.eof();
}

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::getLineWithComments
// Purpose       : Read the line (one character at a time via calling
//                 NextChar_),  and split the line into fields (via calling
//                 SplitLine).  This function is the same as getLine, except
//                 it doesn't skip comments and blank lines.
// Special Notes :
// Scope         : public
// Creator       : Keith Santarelli
// Creation Date : 12/05/2007
//-------------------------------------------------------------------------
int SpiceSeparatedFieldTool::getLineWithComments(std::vector<StringToken>& line)
{
  char c=0;
  bool NextCharFlag;
  const std::string nonid(" \t\n\r");
  line.clear();

  bool endOfLine = false; //is the line end?
  NextCharFlag=in_.get(c);

  while ( endOfLine == false && NextCharFlag  )
  {
    StringToken field;
    field.lineNumber_ = cursorLineNum_;
    field.string_.reserve(16);

    if (nonid.find(c) == nonid.npos) //not a whitespace character
    {
      field.string_ += c;
      while ( (NextCharFlag = in_.get(c)) ) 
        { 
        if (nonid.find(c) == nonid.npos)
          field.string_ += c;
        else
          break;
      }
    }
    else
    { 
      if (c == '\n' || c == '\r')
      {
        endOfLine = true;
        field.string_ += c;
      }
      else
      {
        field.string_ += c;
        NextCharFlag = in_.get(c);
      }
    }

    if (field.string_.length() > 0)
    {
      line.push_back(field);
    }
  }
  // do any needed substitutions
  substituteExternalParams( line );
  return !in_.eof();
}

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::getLine2
// Purpose       : Read the line (one character at a time via calling
//                 NextChar_), and split the line into fields (via calling
//                 SplitLine).  This function differs from getLine in that
//                 it calls SplitLine2.
//
// Special Notes : SplitLine2 ignores the "=" sign.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date :
//-------------------------------------------------------------------------
int SpiceSeparatedFieldTool::getLine2(std::vector<StringToken>& line)
{
  char c=0;
  const std::string nonid(" \t\n\r(),");
  line.clear();

  skipCommentsAndBlankLines_();

  bool endOfLine = false; //is the line end?

  //If either NextChar_(c) == false or endOfLine == true
  //then the line is end. But checking endOfLine need to be
  //in front of NextChar_ to avoid reading next line.
  while ( endOfLine == false && NextChar_(c)  )
  {
    StringToken field;
    field.lineNumber_ = cursorLineNum_;
    field.string_.reserve(16);
    switch (c)
    {
      //Skip whitespaces
      case ' ' :
      case '\t':
      case '\r':
      case '\n':
         break;

         //case '*' : //should never be here
         //cout << fileName_ <<":" <<cursorLineNum_ <<": "
         //  << "Unexpected char '*', ignorng it.\n";
         //break;

      case '=' :
         field.string_ = "=";
         break;

      case ':' :
         field.string_ = ":";
         break;

      case ',' :
          field.string_ = "";
          break;

      case '(' :
          field.string_ = "(";
          break;

      case ')' :
        field.string_ = ")";
        break;

      default : 
      {
        //chars other than above
        // store the char until encount the a nonid char
        field.string_ += c;
        bool metNonIdChar = false; //did it encounter nonid char?

        while (NextChar_(c))
        {
          if (nonid.find(c) == nonid.npos) //c is a id char
          {
             field.string_ += c;
          }
          else //c is nonid, end of this field
          {
            metNonIdChar = true;
            in_.putback(c);
            if (c == '\n') cursorLineNum_--;
            break;
          }
        }

        if (metNonIdChar == false) //then NextChar_(c) == false makes loop end
        {
          endOfLine = true;
        }
      }
    }

    if (field.string_.length() > 0)
    {
     line.push_back(field);
    }
  }
  // do any needed substitutions
  substituteExternalParams( line );
  return !in_.eof();
}

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::changeCursorLineNumber
// Purpose       :
//
// Increase or decrease the line number of cursor by the given amount.
// Increase if token is positive; decrease if token is negative.
// If the result cursor line number is less than 1, set it to 1.
//
// Scope         : public
// Special Notes :
// Creator       : Jian Li
// Creation Date : 07/12/2002
//-------------------------------------------------------------------------
void SpiceSeparatedFieldTool::changeCursorLineNumber(int token)
{
  cursorLineNum_ += token;
  if (cursorLineNum_ < 1) cursorLineNum_ = 1;
}

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::getFilePosition
// Purpose       : Return the current character position in file.
// Scope         : public
// Special Notes :
// Creator       : ??
// Creation Date : ??
//-------------------------------------------------------------------------
std::streampos SpiceSeparatedFieldTool::getFilePosition() const 
{ 
  return in_.tellg(); 
}

//----------------------------------------------------------------------------
// Function       : SpiceSeparatedFieldTool::setLocation
// Purpose        : Set the location in the ifstream at which the next
//                  input operation will begin.
// Special Notes  : 
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/21/2003
//----------------------------------------------------------------------------
bool SpiceSeparatedFieldTool::setLocation(std::streampos const& startPosition)
{
  in_.clear();

  if (startPosition == std::streampos(-1))
  {
    in_.seekg(0,std::ios::end);
  }
  else
  {
    in_.seekg(startPosition);
  }
  return true;
}

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::NextChar_
// Purpose       : Read in a character & check for special characters.
// Special Notes :
// Scope         : private
// Creator       : Alan Lundin
// Creation Date :
//-------------------------------------------------------------------------
bool SpiceSeparatedFieldTool::NextChar_(char& c)
{
  // check common case (i.e., not EOL)
  if (in_.eof()) return false;
  in_.get(c);

  // catch eof w/o newline
  if (in_.eof()) return false;

  // Check for inline comment.
  if (c == ';')
  {
    // Found start of inline comment, gobble up the rest of the line.
    while (c != '\r' && c != '\n')
    {
      if (in_.eof()) return false;
      in_.get(c);
    }
  }

  if (c != '\r' && c != '\n')
  {
    return true;
  }

  // Ok it was end-of-(physical)-line, eat up the line terminator and
  // skip past any subsequent empty lines.
  while ((c == '\r') || (c == '\n'))
  {
    ++cursorLineNum_;
    char oc = c;
    if (in_.eof()) return false;
    in_.get(c);
    if (oc == '\r' && c == '\n')
    {
      if (in_.eof()) return false;
      in_.get(c);
    }
  }

  // Skip subsequent comment lines.
  while ( (c == '*') || (c == ' ') )
  {
    while ( c != '\r' && c != '\n' )
    {
      if (in_.eof()) return false;
      in_.get(c);
    }

    // End of (physical) line, eat up the line terminator and
    // skip past any subsequent empty lines.
    while ((c == '\r') || (c == '\n'))
    {
      ++cursorLineNum_;
      char oc = c;
      if (in_.eof()) return false;
      in_.get(c);
      if ( oc == '\r' && c == '\n' )
      {
        if (in_.eof()) return false;
        in_.get(c);
      }
    }
  }

  // if line is not extended, say line is done
  if (c != '+') 
  {
    in_.putback(c);
    return false;
  }

  // otherwise return a blank for the continuation.
  if (in_.eof()) return false;
  c = ' ';
  return true;
}

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::skipToEndOfLine_
// Purpose       : Helper function to skip to the end of a physical line.
// Special Notes :
// Scope         : private
// Creator       : Raikanta Sahu
// Creation Date : 07/11/2002
//-------------------------------------------------------------------------
void SpiceSeparatedFieldTool::skipToEndOfLine_()
{
  char c(0);
  while ( ! in_.eof() )
  {
    in_.get(c);
    if (in_.eof())
      return;

    if ( c == '\n' )
    {
      ++cursorLineNum_;
      return;
    }
    if ( c == '\r')
    {
      if ( ! in_.eof() )
      {
        in_.get(c);
        if ( c == '\n' )
        {
          ++cursorLineNum_;
          return;
        }
        else
        {
          in_.putback(c);
          return;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------
// Function      : SpiceSeparatedFieldTool::skipCommentsAndBlankLines_
// Purpose       : Helper function to skip the comments and white
//                 spaces before reading a logical Line.
// Special Notes :
// Scope         : private
// Creator       : Raikanta Sahu
// Creation Date : 07/11/2002
//-------------------------------------------------------------------------
void SpiceSeparatedFieldTool::skipCommentsAndBlankLines_()
{
  char c(0);
  while ( ! in_.eof() )
  {
    in_.get(c);
    if (in_.eof())
      continue;
    if ( c == '*' || c == ' ' || c == '\t' )
    {
      skipToEndOfLine_();
      continue;
    }
    if ( c == '\r' )
    {
      continue;
    }
    if ( c == '\n' )
    {
      ++cursorLineNum_ ;
      continue;
    }

    in_.putback(c);
    return;
  }
}

//----------------------------------------------------------------------------
// Function       : SpiceSeparatedFieldTool::substituteExternalParams
// Purpose        : During the parsing of a line to tokens, replace any 
//                  supplied external parameters with the new values 
// Special Notes  : 
// Scope          : 
// Creator        : Richard Schiek, Electrical and MEMS modeling
// Creation Date  : 10/09/2008
//----------------------------------------------------------------------------
void SpiceSeparatedFieldTool::substituteExternalParams(std::vector<StringToken>& line)
{
  bool foundMacroStatement = false;
  std::vector<SpiceSeparatedFieldTool::StringToken>::iterator currField = line.begin();
  std::vector<SpiceSeparatedFieldTool::StringToken>::iterator endField = line.end();
  while( currField != endField )
  {
    std::vector< std::pair< std::string, std::string > >::iterator currentSub = externalParams_.begin();
    std::vector< std::pair< std::string, std::string > >::iterator endSub = externalParams_.end();
    if( currField->string_ == "@selecttext" )
      foundMacroStatement=true;
    while( currentSub != endSub )
    {
      // externalParams can contain section names from dakota like
      // variables 2, responses 4, derivatives 4.
      // If we encounter the header "variables" then only use elements in that 
      // section for substitution.  Exit if we hit one of the other keywords 
      // keywords are "variables" "functions" "derivative_variables" and "analysis_components"
      if( currentSub->first == "functions" ||
          currentSub->first == "derivative_variables" ||
          currentSub->first == "analysis_components" )
      {
        break;  // exit substitution loop now.
      }
      else if( currentSub->first == "variables" )
      {
        // just continue the loop;
      }
      else if( currentSub->first == currField->string_ )
      {
        // found a match so do a substitution.
        currField->string_ = currentSub->second;
      }
      currentSub++;
    }
    currField++;
  }
  
  if( foundMacroStatement )
  {
    // After substitution of values we should be able to process the macro 
    currField = line.begin();
    while( currField != endField )
    {
      // find location of macro 
      if( currField->string_ == "@selecttext" )
      {
        // from her on will need to be replaced so save this location 
        std::vector<SpiceSeparatedFieldTool::StringToken>::iterator resultLoc = currField;
        // look ahead to the selection field
        currField++;
        if( currField != endField )
        {
          currField++;
          if( currField != endField )
          {
            // should be at selection key 
            std::string selectionKey = currField->string_;
            int selectionKeyValue = ExtendedString( currField->string_ ).Value();
            currField++;
            while( currField != endField )
            {
              // since dakota works with values and not strings, we need to compare 
              // these as numbers or "0" will not equal "0.00000000E+00" 
              ExtendedString es( currField->string_ );
              if( es.isValue() )
              {
                int currFieldValue = es.Value();
                if( currFieldValue == selectionKeyValue )
                {
                  // found the right key.  The selected text will be 
                  // at currField+2
                  currField++;
                  currField++;
                  if( currField != endField )
                  {
                    // set value of result location;
                    resultLoc->string_ = currField->string_;
                    // erase the rest of the vector 
                    line.erase( ++resultLoc, endField);
                    break;
                  }
                }
              }
              currField++;
            }
          }
          else
          {
            break;
          }
        }
        else
        {
          break;
        }
        // may have hit the end of the line at this point so the 
        // next currField++ would put us past the end.  If so, then exit the loop 
        if( currField == endField )
          break;
      }
    
      currField++;
    }
  }
  
  return;
}

//----------------------------------------------------------------------------
// Function       : SpiceSeparatedFieldTool::StringToken::instance()
// Purpose        : 
// Special Notes  : 
// Scope          : 
// Creator        : Lon Waters
// Creation Date  : 07/11/2003
//----------------------------------------------------------------------------
Packable* SpiceSeparatedFieldTool::StringToken::instance() const
{
  return new SpiceSeparatedFieldTool::StringToken();
}

//----------------------------------------------------------------------------
// Function       : SpiceSeparatedFieldTool::StringToken::packedByteCount
// Purpose        : 
// Special Notes  : 
// Scope          : 
// Creator        : Lon Waters
// Creation Date  : 07/11/2003
//----------------------------------------------------------------------------
int SpiceSeparatedFieldTool::StringToken::packedByteCount() const
{
  int byteCount(0);
  int length;

  // count string_
  length = string_.length();
  byteCount += length + sizeof(int);

  // count lineNumber_
  byteCount += sizeof(int);

  return byteCount;
}

//----------------------------------------------------------------------------
// Function       : SpiceSeparatedFieldTool::StringToken::pack
// Purpose        : 
// Special Notes  : 
// Scope          : 
// Creator        : Lon Waters
// Creation Date  : 07/11/2003
//----------------------------------------------------------------------------
void SpiceSeparatedFieldTool::StringToken::pack(
    char* buf, int bsize, int& pos, N_PDS_Comm* comm) const
{
  int length;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  // pack string_
  length = string_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( string_.c_str(), length, buf, bsize, pos );

  // pack lineNumber_
  int lineNum = lineNumber_;
  comm->pack( &lineNum, 1, buf, bsize, pos );
#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    DevelFatal(*this, "SpiceSeparatedFieldTool::StringToken::pack") << "Predicted pos does not match actual pos";
  }
#endif
}



//----------------------------------------------------------------------------
// Function       : SpiceSeparatedFieldTool::StringToken::unpack
// Purpose        : 
// Special Notes  : 
// Scope          : 
// Creator        : Lon Waters
// Creation Date  : 07/11/2003
//----------------------------------------------------------------------------
void SpiceSeparatedFieldTool::StringToken::unpack(
    char* pB, int bsize, int& pos, N_PDS_Comm* comm)
{
  int length;

  // unpack string_
  comm->unpack( pB, bsize, pos, &length, 1 );
  string_ = std::string( (pB+pos), length);
  pos += length;

  // unpack lineNumber_
  int lineNum;
  comm->unpack( pB, bsize, pos, &lineNum, 1 );
  lineNumber_ = lineNum;
}

} // namespace IO
} // namespace Xyce
