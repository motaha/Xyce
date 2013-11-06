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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_IO_FunctionBlock.C,v $
//
// Purpose        : Define an N_IO_FunctionBlock instantiations of which are
//                  associated with netlist .FUNC lines.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 12/26/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.28.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:43 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#include <iostream>

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif

// ----------   Xyce Includes   ----------
#include <N_IO_CircuitBlock.h>
#include <N_UTL_Misc.h>
#include <N_IO_FunctionBlock.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Expression.h>
#include <N_PDS_Comm.h>

//-----------------------------------------------------------------------------
// Function      : N_IO_FunctionBlock::N_IO_FunctionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/26/2001
//-----------------------------------------------------------------------------
N_IO_FunctionBlock::N_IO_FunctionBlock(
    string const& fileName,
    vector<N_IO_SpiceSeparatedFieldTool::StringToken> const& parsedInputLine)
: netlistFileName_(fileName),
  parsedLine(parsedInputLine)
{
  int len;

  len = parsedLine[parsedLine.size()-1].string_.size();
  if (parsedLine[parsedLine.size()-1].string_.substr(0,1) != "{" ||
      parsedLine[parsedLine.size()-1].string_.substr(len-1,1) != "}") {
    string msg("In .func line for function: ");
    msg += parsedLine[1].string_ + ", expression must be enclosed by curly braces\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }
}

//----------------------------------------------------------------------------
// Function       : N_IO_FunctionBlock::N_IO_FunctionBlock
// Purpose        : copy constructor
// Special Notes  : 
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/26/2003
//----------------------------------------------------------------------------
N_IO_FunctionBlock::N_IO_FunctionBlock( N_IO_FunctionBlock const& rhsFB )
  : netlistFileName_(rhsFB.netlistFileName_),
    parsedLine(rhsFB.parsedLine),
    functionName(rhsFB.functionName),
    functionNameAndArgs(rhsFB.functionNameAndArgs),
    functionArgs(rhsFB.functionArgs),
    functionBody(rhsFB.functionBody)
{
}

//-----------------------------------------------------------------------------
// Function      : N_IO_FunctionBlock::print
// Purpose       : Output the details of a function block to standard out.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/26/2001
//-----------------------------------------------------------------------------
void N_IO_FunctionBlock::print()
{
  cout << endl;
  cout << "Function Information" << endl;
  cout << "--------------------" << endl;
  cout << "  name: " << functionName << endl;
  cout << "  name and args: " << functionNameAndArgs << endl;
  cout << "  body: " << functionBody << endl;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_FunctionBlock::extractData
// Purpose       : Extract function data from parsed line formed from a
//                 netlist .func statement.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/26/2001
//-----------------------------------------------------------------------------
bool N_IO_FunctionBlock::extractData()
{
  ExtendedString ES1("");

  // Set the function name.
  ES1 = parsedLine[1].string_;
  ES1.toUpper();
  functionName = ES1;

  int arg_start = 2; // Start position of the function's argument list.
  int iend = parsedLine.size();
  int arg_end = iend - 2; // End position of the argument list.

  // The argument list must be enclosed by parentheses.
  if ( (parsedLine[arg_start].string_ != "(") || 
       (parsedLine[arg_end].string_ != ")") )
  {
    string msg(".FUNC argument list must be enclosed by parentheses");
    msg += " in function " + functionName + "\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[arg_start].lineNumber_);
  }

  // Collect up the function name and arguments (with enclosing parentheses),
  // since the expression class needs functions in this form.
  ES1 = functionName + "(";

  // Get the list of arguments.
  for ( int i = arg_start+1; i <= arg_end-1; ++i )
  {
    ES1 += parsedLine[i].string_;

    if (parsedLine[i].string_ != ",") 
      functionArgs.push_back(ExtendedString(parsedLine[i].string_).toUpper());
  }

  // Add the closing parenthese.
  ES1 += ")";

  ES1.toUpper();
  functionNameAndArgs = ES1;

  // Get the function body.
  ES1 =  parsedLine[iend - 1].string_;
  ES1.toUpper();
  functionBody = ES1;

  return true; // Only get here on success.
}


//-----------------------------------------------------------------------------
// Function      : N_IO_FunctionBlock::instance
// Purpose       : implement packable
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Packable * N_IO_FunctionBlock::instance() const
{
  return new N_IO_FunctionBlock();
}


//-----------------------------------------------------------------------------
// Function      : N_IO_FunctionBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int N_IO_FunctionBlock::packedByteCount() const
{
  int byteCount = 0;
  int size, j;

  // count parsedLine
  size = parsedLine.size();
  byteCount += sizeof( int );
  for( j = 0; j < size; ++j)
  {
    byteCount += parsedLine[ j ].packedByteCount();
  }

  // count functionName
  byteCount += sizeof( int );
  byteCount += functionName.length();

  // count functionNameAndArgs
  byteCount += sizeof( int );
  byteCount += functionNameAndArgs.length();

  // count functionArgs
  size = functionArgs.size();
  byteCount += sizeof( int );
  for( j = 0; j < size; ++j)
  {
    byteCount += sizeof( int );
    byteCount += functionArgs[ j ].length();
  }

  // count functionBody
  byteCount += sizeof( int );
  byteCount += functionBody.length();  

  // count netlistFileName_
  byteCount += sizeof( int );
  byteCount += netlistFileName_.length();  
  
  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_FunctionBlock::pack
// Purpose       : Packs function block into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_IO_FunctionBlock::pack( char * buf, int bsize, int & pos,
 N_PDS_Comm * comm ) const
{
  int size, length, j;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  // pack parsedLine
  size = parsedLine.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( int k = 0; k < size; ++k)
  {
    parsedLine[ k ].pack( buf, bsize, pos, comm );
  }

  // pack functionName;
  length = functionName.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( functionName.c_str(), length, buf, bsize, pos );

  // pack functionNameAndArgs;
  length = functionNameAndArgs.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( functionNameAndArgs.c_str(), length, buf, bsize, pos );

  // pack functionArgs
  size = functionArgs.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( int k = 0; k < size; ++k)
  {
    length = functionArgs[ k ].length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( functionArgs[ k ].c_str(), length, buf, bsize, pos );
  }

  // pack functionBody;
  length = functionBody.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( functionBody.c_str(), length, buf, bsize, pos );

  // pack netlistFileName_
  length = netlistFileName_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( netlistFileName_.c_str(), length, buf, bsize, pos );  
#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    string msg("Predicted pos does not match actual pos in N_IO_FunctionBlock::pack");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg );
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : N_IO_FunctionBlock::unpack
// Purpose       : Unpacks function block from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_IO_FunctionBlock::unpack( char * pB, int bsize, int & pos,
 N_PDS_Comm * comm)
{
  int length, size, j;

  // unpack parsedLine
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( j = 0; j < size; ++j )
  {
    N_IO_SpiceSeparatedFieldTool::StringToken aStringToken;
    aStringToken.unpack( pB, bsize, pos, comm );
    parsedLine.push_back( aStringToken );
  }

  // unpack function name
  comm->unpack( pB, bsize, pos, &length, 1 );
  functionName = string( ( pB + pos ), length );
  pos += length;

  // unpack function name and args
  comm->unpack( pB, bsize, pos, &length, 1 );
  functionNameAndArgs = string( ( pB + pos ), length );
  pos += length;

  // unpack function args
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( j = 0; j < size; ++j )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    functionArgs.push_back( string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack function body
  comm->unpack( pB, bsize, pos, &length, 1 );
  functionBody = string( ( pB + pos ), length );
  pos += length;

  // unpack netlistFileName
  comm->unpack( pB, bsize, pos, &length, 1 );
  netlistFileName_ = string( ( pB + pos ), length );
  pos += length;
}

