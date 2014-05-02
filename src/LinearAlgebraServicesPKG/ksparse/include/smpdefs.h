/*
//-------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
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
// Filename       : $RCSfile: smpdefs.h,v $
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
// Revision Number: $Revision: 1.4 $
//
// Revision Date  : $Date: 2013/10/03 18:43:30 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
*/

#ifndef KSPARSE_SMP
#define KSPARSE_SMP

typedef  char SMPmatrix;
typedef  struct MatrixElement  *SMPelement;

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "complex.h"
#include <stdio.h>

#ifdef __STDC__
int SMPaddElt( SMPmatrix *, int , int , double );
void SMPcClear( SMPmatrix *);
int SMPcLUfac( SMPmatrix *, double );
int SMPcProdDiag( SMPmatrix *, SPcomplex *, int *);
int SMPcReorder( SMPmatrix * , double , double , int *);
int SMPcSolve( SMPmatrix *, double [], double [], double [], double []);
void SMPclear( SMPmatrix *);
void SMPcolSwap( SMPmatrix * , int , int );
void SMPdestroy( SMPmatrix *);
int SMPfillup( SMPmatrix * );
SMPelement * SMPfindElt( SMPmatrix *, int , int , int );
void SMPgetError( SMPmatrix *, int *, int *);
int SMPluFac( SMPmatrix *, double , double );
double * SMPmakeElt( SMPmatrix * , int , int );
int SMPmatSize( SMPmatrix *);
int SMPnewMatrix( SMPmatrix ** );
int SMPnewNode( int , SMPmatrix *);
int SMPpreOrder( SMPmatrix *);
void SMPprint( SMPmatrix * , char *);
int SMPreorder( SMPmatrix * , double , double , double );
void SMProwSwap( SMPmatrix * , int , int );
int SMPsolve( SMPmatrix *, double [], double []);
#else /* stdc */
int SMPaddElt();
void SMPcClear();
int SMPcLUfac();
int SMPcProdDiag();
int SMPcReorder();
int SMPcSolve();
void SMPclear();
void SMPcolSwap();
void SMPdestroy();
int SMPfillup();
SMPelement * SMPfindElt();
void SMPgetError();
int SMPluFac();
double * SMPmakeElt();
int SMPmatSize();
int SMPnewMatrix();
int SMPnewNode();
int SMPpreOrder();
void SMPprint();
int SMPreorder();
void SMProwSwap();
int SMPsolve();
#endif /* stdc */

#endif /*KSPARSE_SMP*/
