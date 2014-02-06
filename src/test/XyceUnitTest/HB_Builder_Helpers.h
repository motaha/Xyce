//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: HB_Builder_Helpers.h,v $
// Purpose       : This file contains some helper functions for create N_HB_Builders.
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/10/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.2 $
// Revision Date  : $Date: 2008/09/18 17:04:27 $
// Current Owner  : $Author: tscoffe $
//-----------------------------------------------------------------------------

#ifndef HB_BUILDER_HELPERS_H
#define HB_BUILDER_HELPERS_H

#include <Teuchos_RefCountPtr.hpp>

using Teuchos::RefCountPtr;

class N_HB_Builder;
RefCountPtr<N_HB_Builder> createHBBuilder(int numMPDEBlocks, int numSolutionVars, int numStateVars);

#endif // HB_BUILDER_HELPERS_H

