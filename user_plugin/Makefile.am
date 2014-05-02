# -------------------------------------------------------------------------
#  Copyright Notice
#
#    Copyright 2002 Sandia Corporation. Under the terms
#    of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
#    Government retains certain rights in this software.
#
#     Xyce(TM) Parallel Electrical Simulator
#     Copyright (C) 2002-2014 Sandia Corporation
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
#  Filename       : $RCSfile: Makefile.am,v $
#
#  Purpose        :
#  This automake Makefile deomstrate how to compile and link a runtime
#  loaded sharable object from a Verilog module.  The plugin target can
#  be used to run the verilog file through the ADMS compiler, build the
#  device registration bootstrap object and link them into a sharable
#  object. 
#
#  Special Notes  : 
#
#  Creator        : David G Baur
#
#  Creation Date  : 3/11/2014
#
#  Revision Information:
#  ---------------------
#
#  Revision Number: $Revision $
#
#  Revision Date  : $Date: 2014/03/12 16:40:55 $
#
#  Current Owner  : $Author $
# -------------------------------------------------------------------------

AM_CPPFLAGS = @Xyce_INCS@

lib_LTLIBRARIES = libADMSrlc.la

plugin:	libADMSrlc.la

install-plugin: install

# NOTE:
# The Xyce/ADMS back-end creates files that are named after the *module*
# in the verilog file, and this module name may have no relation to the
# name of the *file* in which the verilog is stored.  The .C and .h file
# produced will have the text "N_DEV_ADMS" prepended to the module name
# This prepended text may change in future releases.
#
# It is for this reason that the rule below is not a template rule for
# converting arbitrary .va files into their associated .C file.  The rule
# depends very much on what the module name is, not what the source file
# name is.
#
N_DEV_ADMSrlc.C N_DEV_ADMSrlc.h: $(srcdir)/rlc.va
	admsXml -e $(srcdir)/../utils/ADMS/xyceVersion.xml -e $(srcdir)/../utils/ADMS/xyceBasicTemplates.xml -e $(srcdir)/../utils/ADMS/xyceHeaderFile.xml -e $(srcdir)/../utils/ADMS/xyceImplementationFile.xml $(srcdir)/rlc.va


# While simple devices can have their .C and .h files generated in a single
# run of ADMS (as , there is currently an issue in the ADMS back-end that 
# makes this fail for some devices. The issue is with how node collapse 
# is implemented in the back end.  The PSP model is the only example we
# have so far that exposes this issue.
#
# In that case, it is necessary to run admsXml twice, once with 
# "xyceImplementationFile.xml" to generate the .C file, and then a second 
# time with "xyceHeaderFile.xml" to generate the .h
# This commented-out section shows how to work around this problem if you 
# encounter it:
#N_DEV_ADMSrlc.C: $(srcdir)/rlc.va
#	admsXml -e $(srcdir)/../utils/ADMS/xyceVersion.xml -e $(srcdir)/../utils/ADMS/xyceBasicTemplates.xml -e $(srcdir)/../utils/ADMS/xyceImplementationFile.xml $(srcdir)/rlc.va
#
#N_DEV_ADMSrlc.h: $(srcdir)/rlc.va
#	admsXml -e $(srcdir)/../utils/ADMS/xyceVersion.xml -e $(srcdir)/../utils/ADMS/xyceBasicTemplates.xml -e $(srcdir)/../utils/ADMS/xyceHeaderFile.xml $(srcdir)/rlc.va
#

N_DEV_ADMSrlc_bootstrap.lo: N_DEV_ADMSrlc.h

CLEANFILES=N_DEV_ADMSrlc.C N_DEV_ADMSrlc.h

libADMSrlc_la_SOURCES = N_DEV_ADMSrlc.C N_DEV_ADMSrlc_bootstrap.C
