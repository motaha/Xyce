// Front definitions 

`include "discipline.h"

// Numerical, physical and model constants 
`define TEN_M40 1.0e-40
`define TEN_M07 1.0e-7
`define C2K     273.15
`define KB      1.3806226e-23
`define QQ      1.6021918e-19
`define KBdivQQ 8.61708691805812512584e-5
`define VDLOW   0.05
`define AJE     3.0
`define AJC     2.0
`define AJS     2.0
`define VEXLIM  80.0
`define TMIN    20.0
`define TMAX    600.0
`define PI      3.1415926

// Desriptions and units
//`ifdef __VAMS_COMPACT_MODELING__
  `define OPP(nam,uni,des) (* desc="des", units="uni" *) real nam;
  `define PAR(des,uni) (* desc="des", units="uni" *) parameter real
  `define PAI(des,uni) (* desc="des", units="uni" *) parameter integer
//`else
//  `define OPP(nam,uni,des);
//  `define PAR(des,uni) parameter real
//  `define PAI(des,uni) parameter integer
//`endif

// Smooth limitting functions
`define max_hyp0(result, x, epsilon)\
     eps2 = epsilon * epsilon;\
     x2 = x * x;\
     if (x < 0.0)\
       result = 0.5 * eps2 / (sqrt(x2 + eps2) - x);\
     else\
       result = 0.5 * (sqrt(x2 + eps2) + x);

`define min_logexp(result, x, x0, a)\
     dxa = (x - x0) / (a);\
     if (x < x0)\
          result = x  - a * ln(1.0 + exp(dxa));\
     else\
          result = x0 - a * ln(1.0 + exp(-dxa));

`define max_logexp(result, x, x0, a)\
     dxa = (x - x0) / (a);\
     if (x < x0)\
          result = x0 + a * ln(1.0 + exp(dxa));\
     else\
          result = x  + a * ln(1.0 + exp(-dxa));

`define expLin(result, x)\
     xtmp=x;\
     if (xtmp < `VEXLIM)\
         result = 1.0;\
     else begin\
         result = xtmp - `VEXLIM + 1.0;\
         xtmp = `VEXLIM;\
     end\
     result = result * exp(xtmp);

`define linLog(result, x, vlim)\
      if (x < vlim)\
          result = x;\
      else\
          result = vlim + ln(1.0 + (x - vlim));

