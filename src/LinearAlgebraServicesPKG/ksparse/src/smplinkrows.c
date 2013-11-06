/*
//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002-2011, Sandia Corporation, Albuquerque, NM, USA.
// Under the terms of Contract DE-AC04-94AL85000, there is a
// non-exclusive license for use of this work by or on behalf of the
// U.S. Government.  Export of this program may require a license from
// the United States Government.
// 
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: smplinkrows.c,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2011/08/17 21:20:54 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
*/

#include "spice.h"
#include <stdio.h>
#include "spmatrix.h"
#include "smpdefs.h"
#include "spdefs.h"

void
SMPlinkRows( Matrix )
SMPmatrix *Matrix;
{
    void spcLinkRowsandCreateInternalVectors(char *);

    spcLinkRowsandCreateInternalVectors( (char *)Matrix );
}

void spcLinkRowsandCreateInternalVectors( eMatrix )
char *eMatrix;
{
    MatrixPtr  Matrix = (MatrixPtr)eMatrix;
    
    if (NOT Matrix->RowsLinked)
        spcLinkRows( Matrix );
    if (NOT Matrix->InternalVectorsAllocated)
        spcCreateInternalVectors( Matrix );
}
