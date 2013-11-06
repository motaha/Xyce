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
// Filename       : $RCSfile: driver.c,v $
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
#include "spmatrix.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define ABS_TOL 1.e-10
#define REL_TOL 1.e-6

void check (char *, double *, double *, int);

int main(int argc, char* argv[])
{
   int i, j, size, err, stat;
   FILE *mat, *vec;
   char *m_name, *v_name;
   void *Matrix;
   double *rhs, *solution, *answer, *p, *val_list, **addr_list;
   double value;
   int list_d, n_vals;

   if (argc < 2) {
     fprintf (stderr, "Please type filename on command line (e.g. 'driver ltra3')\n");
     exit(-1);
   }
   m_name = malloc(strlen(argv[1])+5);
   v_name = malloc(strlen(argv[1])+5);
   strcpy (m_name,argv[1]);
   strcat (m_name,".mat");
   strcpy (v_name,argv[1]);
   strcat (v_name,".vec");

   mat = fopen (m_name,"r");
   vec = fopen (v_name,"r");

   fscanf (mat,"%d",&size);

/* note that sparse ignored the zeroth element of these vectors */
   rhs = (double *) malloc((size+1)*sizeof(double));
   solution = (double *) malloc((size+1)*sizeof(double));
   answer = (double *) malloc((size+1)*sizeof(double));

/* Create matrix structure */
   Matrix = (void *) spCreate (size, 0, &err);
   if (err) {
     fprintf (stderr, "Failed to create matrix\n");
   }

   for (i=1 ; i<=size; i++) {
     (void) fscanf (vec, "%lg %lg", &rhs[i], &answer[i]);
   }

   n_vals = 0;
   list_d = 1000;
   val_list = (double *) malloc(list_d*sizeof(double));
   addr_list = (double **) malloc(list_d*sizeof(double *));
   while (fscanf(mat, "%d %d %lg", &i, &j, &value) == 3) {
     p = (double *) spGetElement(Matrix, i, j);
     *p = value;
     if (n_vals == list_d) {
       list_d *= 2;
       val_list = (double *) realloc(val_list, list_d*sizeof(double));
       addr_list = (double **) realloc(addr_list, list_d*sizeof(double *));
     }
     val_list[n_vals] = value;
     addr_list[n_vals] = p;
     n_vals++;
   }
   fclose(mat);
   fclose(vec);

   spOrderAndFactor (Matrix, rhs, 0.001, 1.e-13, 1);
   spSolve (Matrix, rhs, solution, NULL, NULL);
   check ("initial solve", solution, answer, size);

   memcpy (solution, rhs, (size+1)*sizeof(double));
   spClear (Matrix);
   for (i=0 ; i<n_vals ; i++) {
     *(addr_list[i]) = val_list[i];
   }
   spFactorAndSolve (Matrix, solution);
   check ("second solve", solution, answer, size);

   memcpy (solution, rhs, (size+1)*sizeof(double));
   spClear (Matrix);
   for (i=0 ; i<n_vals ; i++) {
     *(addr_list[i]) = val_list[i];
   }
   spFactorAndSolve (Matrix, solution);
   check ("third solve", solution, answer, size);

   return 0;
}

void check (char *msg, double *solution, double *answer, int size)
{
   int err, i;
   double del;

   err = 0;
   for (i=1 ; i<=size; i++) {
     if (fabs(solution[i]-answer[i]) > ABS_TOL) {
       if (fabs(answer[i]) > ABS_TOL) {
         del = (solution[i]-answer[i]);
         if (fabs(del/answer[i]) > REL_TOL) {
           printf ("possible error at element: %d, my value = %.12lg, answer = %.12lg\n",i,solution[i], answer[i]);
           err++;
         }
       }
     }
   }
   if (err) {
     printf ("%d possible errors reported", err);
   }
   else {
     printf ("No errors found");
   }
   printf (" for %s\n",msg);

   return;
}
