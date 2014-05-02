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
// Filename       : $RCSfile: spice.h,v $
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

/*
 *	Portability, global externs, kitchen sink (yuk).
 *	In future releases this will include things from misc.h and util.h,
 *	   which duplicate each other in places
 */
#ifndef KSPARSE_SPICE_H
#define KSPARSE_SPICE_H

#include <Xyce_config.h>

#ifdef HAVE_MATH_H
#include <math.h>
#endif

#ifndef	M_PI
#  define M_PI		3.14159265358979323846
#endif
#ifndef	M_E
#  define M_E  	   2.7182818284590452354
#endif
#ifndef	M_LOG2E
#  define M_LOG2E		1.4426950408889634074
#endif
#ifndef	M_LOG10E
#  define M_LOG10E        0.43429448190325182765
#endif
#define TS_LIMIT 5

#include "hw.h"
#include "config.h"
#include "capabil.h"

#define	NUMELEMS(ARRAY)	(sizeof(ARRAY)/sizeof(*ARRAY))

extern char *Spice_Exec_Dir;
extern char *Spice_Lib_Dir;
extern char *Spice_Help_Dir;
extern char *Spice_Model_Dir;
extern char Spice_OptChar;
extern char *Def_Editor;
extern char *Bug_Addr;
extern int AsciiRawFile;
extern char *Spice_Host;
extern char *Spiced_Log;

extern char Spice_Notice[ ];
extern char Spice_Version[ ];
extern char Spice_Build_Date[ ];

extern char *News_File;
extern char *Default_MFB_Cap;
extern char *Spice_Path;
extern char *Help_Path;
extern char *Lib_Path;
extern int  Patch_Level;

#ifdef MAIN_PROGRAM
    int report_interval, new_raw_head, load_mode, mat_dense;
    int device_error, model_error, *timer_flag, dev_math_error;
    double *simulation_time, *logic_break, *timer_calibration;
    double *accepted_simulation_time;
    int exp_num, exp_d;
    void **exp_list, **exp_values;

    char *xfile;
#else
    extern int report_interval, new_raw_head, load_mode, mat_dense;
    extern int device_error, model_error, *timer_flag, dev_math_error;
    extern double *simulation_time, *logic_break, *timer_calibration;
    extern double *accepted_simulation_time;
    extern char *xfile;
    extern int exp_num, exp_d;
    extern void **exp_list, **exp_values;
#endif

#ifdef SHARED_MEM
/*
void update_SM (double *, double);
*/
#endif /* SHARED_MEM */

#define LOAD_NORMAL 1
#define LOAD_ENERGY 2

#endif
