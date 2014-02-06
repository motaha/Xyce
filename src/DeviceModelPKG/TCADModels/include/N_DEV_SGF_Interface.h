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
// Filename       : $RCSfile: N_DEV_SGF_Interface.h,v $
//
// Purpose        : This file contains auxilliary classes, etc. which 
//                  are neccessary for outputting sgplot files, and 
//                  also for reading SGF mesh files.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/08/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_SGF_INTERFACE_h
#define Xyce_N_DEV_SGF_INTERFACE_h

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------

// ---------- Forward Declarations -------

//-----------------------------------------------------------------------------
//  Legacy class definitions needed for sgplot output.
//-----------------------------------------------------------------------------

#define LEN_IDENT           15          // identifier length
#define TYPE_ICONST         4           // integer constant

class RESHEAD
{
  public:
  char szLogo[64];                      // text logo
  char szSign[16];                      // signature
  char szMeshFile[128];                 // name of mesh file
  UINT cConstant;                       // number of constants
  UINT cVariable;                       // number of variables
  UINT cArray;                          // number of arrays
  UINT cElement;                        // number of elements
  UINT c1DArray;                        // number of 1D arrays
  UINT c2DArray;                        // number of 2D arrays
  UINT c3DArray;                        // number of 3D arrays
  UINT cSet;                            // number of data sets
};

class XLATARRAY
{
  public:
  char szName[LEN_IDENT+1];             // array name
  UINT uOffset;                         // array offset
  UINT cDim;                            // number of dimensions
  UINT acElements[3];                   // elements per dimensions
};

class XLATVAR
{
  public:
  char szName[LEN_IDENT+1];             // variable name
  UINT uOffset;                         // variable offset
};



//-----------------------------------------------------------------------------
// Class         : DAXLATARRAY
// Purpose       : Dynamic array of the XLATARRAY class
//
// Special Notes : This is only used for output sgplot style binary files.
//                 Otherwise, ignore.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/29/02
//-----------------------------------------------------------------------------
class DAXLATARRAY
{ private:
    UINT cSize;               // total number of elements
    UINT uInc;                // increment size
    UINT cElements;           // number of elements in array
    XLATARRAY *dynarray;      // dynamic array

    void IncSize ()
    { if (cElements == cSize)
      { cSize += uInc;
        dynarray =
        (XLATARRAY *) realloc(dynarray, cSize*sizeof(XLATARRAY));
      }
    }

  public:
    void set            (const char *,UINT ,UINT,UINT,UINT,UINT);
    UINT GetElements ()       { return cElements;    }
    UINT Add (XLATARRAY t)    { IncSize(); dynarray[cElements] = t;
                                return cElements++; }

    void RemoveLast ()        { --cElements;         }
    void Flush ()             { cElements = 0;       }

    XLATARRAY *GetPointer (UINT i) { return dynarray + i; }
    XLATARRAY operator[]  (UINT i) { return dynarray[i];  }

  public:
     DAXLATARRAY (int numMeshPoints);
    ~DAXLATARRAY ()                { delete [] dynarray; }
};

#endif // Xyce_N_DEV_SGF_INTERFACE_h

