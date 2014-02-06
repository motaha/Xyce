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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_2DPDESetup.C,v $
//
// Purpose        : This file mostly contains functions that are called
//                  once, during the initial setup phase of the 2D PDE
//                  device.  There are a couple of exceptions - the mesh
//                  resize functions are called during "re-set-up" phases
//                  of a sensitivity calculation.
//
//                  One very important setup function - processParams - is
//                  *not* in this file.  It gets its own file,
//                  N_DEV_2DPDEParam.C.
//
//                  All of the of the functions pertaining to the global
//                  and local ID setup are in this file.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/05/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.42.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_2DPDE.h>
#include <N_DEV_SolverState.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_DEV_PDE_2DMesh.h>
#include <N_DEV_PDE_Electrode.h>

namespace Xyce {
namespace Device {
namespace TwoDPDE {

//-----------------------------------------------------------------------------
// Function      : Instance::doSensMeshResize
// Purpose       :
// Special Notes : Generally, this will be called for a mesh that was
//                 already scaled, so the resized mesh should be considered
//                 scaled as well.
//
//                 As should be obvious from the name, this function is
//                 designed for perturbing the mesh as part of a
//                 sensitivity calculation.  As such, a copy of the
//                 original mesh is saved, to be restored after the
//                 calculation is over.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
bool Instance::doSensMeshResize ()
{
  bool bsuccess = true;
  bool bs1 = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "In Instance::doSensMeshResize." << endl;
  }
#endif

  // make a copy of the mesh.  This will need to be restored later.
  if (meshCopyContainerPtr == NULL)
  {
    meshCopyContainerPtr = new PDE_2DMesh (*meshContainerPtr);
  }
  else
  {
    *meshCopyContainerPtr = *meshContainerPtr;
  }

  // scale the new size:
  if (variablesScaled)
  {
    deviceLength /= scalingVars.x0;
    deviceWidth  /= scalingVars.x0;
  }

  // get the old size:
  double old_length = meshContainerPtr->getXMax () -
                      meshContainerPtr->getXMin ();

  double old_width  = meshContainerPtr->getYMax () -
                      meshContainerPtr->getYMin ();

#ifdef Xyce_DEBUG_DEVICE
  double lengthRatio  = deviceLength/old_length;
  double widthRatio   = deviceWidth/old_width;

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  x0           = " << scalingVars.x0 << endl;
    cout << "  deviceWidth  = " << deviceWidth  << endl;
    cout << "  old_width    = " << old_width    << endl;
    cout << "  widthRatio   = " << widthRatio   << endl;
    cout << "  deviceLength = " << deviceLength << endl;
    cout << "  old_length   = " << old_length   << endl;
    cout << "  lengthRatio  = " << lengthRatio  << endl;
  }
#endif

  // first resize the mesh in the mesh class:
  meshContainerPtr->resizeMesh(deviceLength, deviceWidth);

  // Then update all the mesh sized stuff used in the Instance
  // class.
  meshContainerPtr->getXVector(xVec);
  meshContainerPtr->getYVector(yVec);
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1)
  {
    for (int i=0;i<numMeshPoints;++i)
    {
      cout << " x["<<i<<"] = " << xVec[i];
      cout << " y["<<i<<"] = " << yVec[i];
      cout << endl;
    }
  }
#endif

  // now update all the mesh stuff in the 2DPDE class:
  bs1 = setupBCEdgeAreas ();  bsuccess = bsuccess && bs1;
  bs1 = setupMinDXVector ();  bsuccess = bsuccess && bs1;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Done with Instance::doSensMeshResize." << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::undoSensMeshResize
//
// Purpose       : This un-does the damage done by doSensMeshResize.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
bool Instance::undoSensMeshResize ()
{
  bool bsuccess = true;
  bool bs1 = true;
  // switch the mesh copy back into the official ptr.:
  PDE_2DMesh * tmpPtr;

  tmpPtr               = meshContainerPtr;
  meshContainerPtr     = meshCopyContainerPtr;
  meshCopyContainerPtr = tmpPtr;

  // Restore all the mesh sized stuff used in the Instance
  // class.
  meshContainerPtr->getXVector(xVec);
  meshContainerPtr->getYVector(yVec);

  // now update all the mesh stuff in the 2DPDE class:
  bs1 = setupBCEdgeAreas (); bsuccess = bsuccess && bs1;
  bs1 = setupMinDXVector (); bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMesh
//
// Purpose       : This function should only be called once.  It handles
//                 most of the stuff associated with initializing the
//                 mesh class, and all the stuff in Instance
//                 that depends on the mesh.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
bool Instance::setupMesh ()
{

  bool bsuccess = true;

  ///////////////////////////////////////////////////////////////////
  // First straighten out the electrode map, if it exists.
  // The electrode map is one of the arguments that need to be passed into
  // the "internal" mesh setup, so it has to be corrected first.
  vector<DeviceInterfaceNode>::iterator first = dIVec.begin ();
  vector<DeviceInterfaceNode>::iterator last  = dIVec.end   ();
  vector<DeviceInterfaceNode>::iterator iterV;

  if (!(electrodeMap.empty ()))
  {
    // First make the names in electrodeMap consistent with those in dIVec.
    for (iterV=first;iterV!=last; ++iterV)
    {
      if (!(iterV->given)) continue;

      if ( electrodeMap.find(iterV->nName) != electrodeMap.end () )
      {
        electrodeMap[iterV->nName]->name = iterV->eName;
      }
      else
      {
        string msg = "Instance::doMeshBasedInitializations."
        "can't find " + iterV->nName + " in the electrode Map\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
      }
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
       cout << "list of user-specified electrodes:" << endl;
       map<string, PDE_2DElectrode*>::iterator mapIter;
       map<string, PDE_2DElectrode*>::iterator mapStart = electrodeMap.begin();
       map<string, PDE_2DElectrode*>::iterator mapEnd = electrodeMap.end();

       // for ( mapIter = mapStart; mapIter != mapEnd; ++mapIter )
       // {
       //  cout << *(mapIter->second);
       // }
    }
#endif
  }

  ///////////////////////////////////////////////////////////////////
  // Allocate the mesh container.
  meshContainerPtr = new PDE_2DMesh(devOptions, sgplotLevel);

  if (!given("MESHFILE"))
  {
    string msg = "Instance::doMeshBasedInitializations."
    "no mesh file specified.  Setting meshfile=internal.msh\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_INFO_0,msg);
  }

  ///////////////////////////////////////////////////////////////////
  // Now initialize the mesh either as internal or external.
  if (meshFileName != "internal" && meshFileName != "internal.msh")
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << endl;
      cout << "Reading mesh file..." << endl;
    }
#endif
    usingInternalMesh = false;
    meshContainerPtr->initializeMesh (meshFileName);
    cylGeomFlag = meshContainerPtr->cylGeom;
  }
  else
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << endl;
      cout << "Generating internal mesh..." << endl;
    }
#endif
    usingInternalMesh = true;

    string outputMeshFileName = outputName + ".msh";
    meshContainerPtr->initializeInternalMesh
      (numMeshPointsX,
       numMeshPointsY,
       deviceLength,
       deviceWidth,
       numElectrodes,
       outputMeshFileName,
       electrodeMap,
       cylGeomFlag);
  }

  ///////////////////////////////////////////////////////////////////
  numMeshPoints = meshContainerPtr->getNumNodes ();
  numMeshEdges  = meshContainerPtr->getNumEdges ();
  numMeshCells  = meshContainerPtr->getNumCells ();
  numMeshLabels = meshContainerPtr->getNumLabels ();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "\n";
    cout << "Done setting up the mesh." << endl;
    cout << "  numMeshPoints      = " << numMeshPoints << "\n";
    cout << "  numMeshEdges       = " << numMeshEdges << "\n";
    cout << "  numMeshCells       = " << numMeshCells << "\n";

    //cout << "  Vbi                = " << Vbi << "\n";
    //cout << "  Na                 = " << Na << "\n";
    //cout << "  Nd                 = " << Nd << "\n";
    //cout << "  gradedJunctionFlag = " << gradedJunctionFlag << "\n";
    //cout << "  displCurrentFlag   = " << displCurrentFlag << "\n";
    //cout << "  junction width(WJ) = " << WJ << "\n";
    //cout << "  deviceWidth        = " << deviceWidth << endl;
    //cout << "  deviceLength       = " << deviceLength << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupDINodes
// Purpose       : This function does some miscellaneous setup of the
//                 device interface nodes (boundary condition class).
//                 Mostly this is does some final, misc. cleanup.
//
//                 Should be called after setupMesh.  Should only be called once.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/15/03
//-----------------------------------------------------------------------------
bool Instance::setupDINodes ()
{
  bool bsuccess = true;

  vector<DeviceInterfaceNode>::iterator first = dIVec.begin ();
  vector<DeviceInterfaceNode>::iterator last  = dIVec.end   ();
  vector<DeviceInterfaceNode>::iterator iterV;

  // loop through the boundary condition name vector and
  // check if it matches the boundary conditions given in the mesh file.

  for (iterV=first;iterV!=last; ++iterV)
  {
    ExtendedString tmpName = iterV->eName;
    tmpName.toUpper ();

    bool edgeLabelExist = meshContainerPtr->labelNameExist(tmpName);

    // If the device interface was given by the user(netlist)
    // (as a boundary condition), then if it doesn't exist,
    // then the user has made a mistake in setting up the input file and
    // the code should exit.
    if ((iterV->given))
    {
      if ( !(edgeLabelExist) )
      {
        meshContainerPtr->printLabels ();
        string msg = "Instance::setupDINodes: "
       "The boundary condition label "+tmpName+
       " doesn't exist in the mesh file.\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
      }
    }
    else // If this wasn't given, and doesn't exist, that just means it
         // is part of the default list of boundary condition names,
         // and should just be removed from the device interface vector.
         // This does NOT represent a mistake in the netlist.
    {
      if ( !(edgeLabelExist) )
      {
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0)
          cout << "Erasing DI: " << iterV->eName << endl;
#endif
        dIVec.erase (iterV);
      }
    }
  } // end of iterV loop.

  // Copy over a few of the material related variables from the electrode
  // class to the device interface node class.
  if (!(electrodeMap.empty ()))
  {
    // First make the names in electrodeMap consistent with those in dIVec.
    for (iterV=first;iterV!=last; ++iterV)
    {
      if (!(iterV->given)) continue;

      if ( electrodeMap.find(iterV->nName) != electrodeMap.end () )
      {
        // material stuff:
        iterV->material       = electrodeMap[iterV->nName]->material;
        iterV->materialGiven  = electrodeMap[iterV->nName]->materialGiven;
        iterV->oxideBndryFlag = electrodeMap[iterV->nName]->oxideBndryFlag;
        iterV->oxthick        = electrodeMap[iterV->nName]->oxthick;
        iterV->oxcharge       = electrodeMap[iterV->nName]->oxcharge;
      }
      else
      {
        string msg = "Instance::doMeshBasedInitializations."
        "can't find " + iterV->nName + " in the electrode Map\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
      }
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        cout << endl;
        cout << "name = " << iterV->eName << endl;
        cout << " material = " << iterV->material << endl;
        cout << " mat. given = " << iterV->materialGiven << endl;
        cout << " oxide boundary flag = " << iterV->oxideBndryFlag << endl;
        cout << " oxide thickness = " << iterV->oxthick << endl;
        cout << " oxide charge = " << iterV->oxcharge << endl;
      }
#endif
    }
  }

  // Loop over the boundaries, and check if each one is a
  // neumann, or mixed boundary condition for each variable.
  first = dIVec.begin ();
  last  = dIVec.end   ();
  for (iterV=first;iterV!=last; ++iterV)
  {
    ExtendedString tmpName = iterV->nName;
    tmpName.toLower ();

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
      cout << "Testing the neumann stuff.  Name = " << tmpName << endl;
#endif

    if ( tmpBCmap.find(tmpName) != tmpBCmap.end () )
    {
      if (tmpBCmap[tmpName] == "NEUMANN")
      {
        iterV->neumannBCFlagV = true;
        iterV->neumannBCFlagN = true;
        iterV->neumannBCFlagP = true;
      }

      if (tmpBCmap[tmpName] == "MIXED")
      {
        iterV->neumannBCFlagV = false;
        iterV->neumannBCFlagN = true;
        iterV->neumannBCFlagP = true;
      }

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        cout << "Setting the neumann flags of " << tmpName << ":\n";
        cout << "  Vflag = " << iterV->neumannBCFlagV << endl;
        cout << "  Nflag = " << iterV->neumannBCFlagN << endl;
        cout << "  Pflag = " << iterV->neumannBCFlagP << endl;
      }
#endif
    }
  } // end of iterV loop.

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    first = dIVec.begin ();
    last  = dIVec.end   ();
    cout << "Final DI list: " << endl;
    for (iterV=first;iterV!=last; ++iterV)
    {
      cout << "DI name:" << iterV->eName;
      cout << "  The neumann flags are:" << endl;
      cout << "  Vflag =";
      if (iterV->neumannBCFlagV) cout <<" true." << endl;
      else                       cout <<" false." << endl;
      cout << "  Nflag =";
      if (iterV->neumannBCFlagN) cout <<" true." << endl;
      else                       cout <<" false." << endl;
      cout << "  Pflag =";
      if (iterV->neumannBCFlagP) cout <<" true." << endl;
      else                       cout <<" false." << endl;
    }
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::doAllocations
// Purpose       : A whole bunch of resizes.  Should be called after
//                 setupMesh.  Should only be called once.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/15/03
//-----------------------------------------------------------------------------
bool Instance::doAllocations ()
{
  bool bsuccess = true;

  // allocate conductance, capacitance array:
  condVec.resize(numElectrodes);
  capVec.resize(numElectrodes);
  for (int iE=0;iE<numElectrodes;++iE)
  {
    condVec[iE].resize(numElectrodes,0.0);
    capVec[iE].resize(numElectrodes,0.0);
  }

  // Set up a bunch of mesh-based arrays:
  // Local allocations:
  xVec.resize    (numMeshPoints);  meshContainerPtr->getXVector(xVec);
  yVec.resize    (numMeshPoints);  meshContainerPtr->getYVector(yVec);
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1)
  {
    for (int i=0;i<numMeshPoints;++i)
    {
      cout << " x["<<i<<"] = " << xVec[i];
      cout << " y["<<i<<"] = " << yVec[i];
      cout << endl;
    }
  }
#endif

  minDXVec.resize (numMeshPoints);

  areaVec.resize (numMeshPoints);

  VVec.resize      (numMeshPoints);
  RVec.resize      (numMeshPoints);
  SVec.resize      (numMeshPoints);
  totSrcVec.resize (numMeshPoints);
  nnVec.resize     (numMeshPoints);
  npVec.resize     (numMeshPoints);
  CVec.resize      (numMeshPoints);
  dRdpVec.resize   (numMeshPoints);
  dRdnVec.resize   (numMeshPoints);

  elecPenalty.resize (numMeshPoints,0.0);
  holePenalty.resize (numMeshPoints,0.0);
  pdElecPenalty.resize (numMeshPoints,0.0);
  pdHolePenalty.resize (numMeshPoints,0.0);

  unVec.resize  (numMeshPoints, 0.0);
  upVec.resize  (numMeshPoints, 0.0);
  unE_Vec.resize (numMeshEdges, 0.0);
  upE_Vec.resize (numMeshEdges, 0.0);
  tnVec.resize  (numMeshPoints, 0.0);
  tpVec.resize  (numMeshPoints, 0.0);

  // assume every node is owned until told otherwise.
  // These are only really used for boundary conditions.
  // They may, in fact, not be needed... ERK, 11/09/02
  vOwnVec.resize  (numMeshPoints, 1);
  nnOwnVec.resize (numMeshPoints, 1);
  npOwnVec.resize (numMeshPoints, 1);

  displPotential.resize(numMeshPoints);

  if (sgplotLevel > 0) outputVec.resize(numMeshPoints,0.0);

  stateDispl.resize(numMeshPoints);
  stateDispl_owned.resize(numMeshPoints,-1);
  li_stateDispl.resize(numMeshPoints,0);

  EfieldVec.resize (numMeshEdges);
  JnVec.resize (numMeshEdges);
  JpVec.resize (numMeshEdges);
  displCurrent.resize(numMeshEdges);

  dJndn1Vec.resize (numMeshEdges);
  dJndn2Vec.resize (numMeshEdges);
  dJndV1Vec.resize (numMeshEdges);
  dJndV2Vec.resize (numMeshEdges);

  dJpdn1Vec.resize (numMeshEdges);
  dJpdn2Vec.resize (numMeshEdges);
  dJpdV1Vec.resize (numMeshEdges);
  dJpdV2Vec.resize (numMeshEdges);

  Vrowarray.resize (numMeshPoints,-1);  Vcolarray.resize (numMeshPoints);
  Nrowarray.resize (numMeshPoints,-1);  Ncolarray.resize (numMeshPoints);
  Prowarray.resize (numMeshPoints,-1);  Pcolarray.resize (numMeshPoints);

  boundarySten.resize(numMeshPoints,0);
  boundaryStenV.resize(numMeshPoints,0);
  boundaryStenN.resize(numMeshPoints,0);
  boundaryStenP.resize(numMeshPoints,0);
  boundaryTest.resize(numMeshPoints,0);

  li_Vrowarray.resize (numMeshPoints,0);
  li_Nrowarray.resize (numMeshPoints,0);
  li_Prowarray.resize (numMeshPoints,0);

  labelIndex.resize (numMeshPoints,0);
  labelNameVector.reserve(numMeshPoints);

  // allocate the boundary condition arrays:
  vector<DeviceInterfaceNode>::iterator first = dIVec.begin ();
  vector<DeviceInterfaceNode>::iterator last  = dIVec.end   ();
  vector<DeviceInterfaceNode>::iterator iterV;

  for (iterV=first;iterV!=last; ++iterV)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterV->eName);
    int numPoints = labelPtr->mNodeVector.size();
    iterV->numBoundaryPoints = numPoints;
    iterV->VequVec.resize (numPoints,0.0);
    iterV->VbcVec.resize  (numPoints,0.0);
    iterV->nnbcVec.resize (numPoints,0.0);
    iterV->npbcVec.resize (numPoints,0.0);

    for (int i=0;i<iterV->numBoundaryPoints;++i)
    {
      mLabel * labelPtr = meshContainerPtr->getLabel(iterV->eName);
      int nodeIndex = labelPtr->mNodeVector[i];
      iterV->meshGlobalToLocal[nodeIndex] = i;
    }
  }

  // setup the aiEdge vector for tecplot:
  aiEdge.resize(numMeshEdges,0);
  aiEdge_nf.resize(numMeshEdges,0);
  //electrodeEdge.resize(numMeshEdges,0);

  // determine the label index for noflux.
  int nofluxIndex = 0;
  int numLabels = meshContainerPtr->getNumLabels();
  for (int i1=0;i1<numLabels;++i1)
  {
    mLabel * lPtr = meshContainerPtr->getLabel(i1);
    if (lPtr->name == "NOFLUX")
    {
      nofluxIndex = i1;
    }
  }

  // finish the aiEdge stuff.
  UINT iE;
  iNumPlotEdges = 0;
  iNumPlotEdges_nf = 0;
  for(iE = 0;iE<numMeshEdges;++iE)
  {
    mEdge * edgePtr = meshContainerPtr->getEdge(iE);
    UINT uLabel = edgePtr->uLabel;

    mLabel * labelPtr = meshContainerPtr->getLabel(uLabel);

    if(uLabel != -1u)
    {
      if(uLabel != nofluxIndex)
      {
        if(labelPtr->uType == TYPE_EDGE)
        {
          aiEdge[iE] = 1;
          ++iNumPlotEdges;
        }
      }
      else
      {
        if(labelPtr->uType == TYPE_EDGE)
        {
          aiEdge_nf[iE] = 1;
          ++iNumPlotEdges_nf;
        }
      }
    } // uLabel != -1u
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupLabelIndex
// Purpose       : This function sets up the vector, labelIndex
//
// Special Notes : labelIndex is a key, which maps a node index
//                 to a label index.  This is needed information, so that the
//                 load routines can determine what form rhs and
//                 Jacobian entries they should sum for a given node. - Should
//                 they set it up as a interior point, or a boundary, for
//                 example.
//
//                 By default, all nodes are assumed to belong to region label.
//                 However, some nodes will also belong to
//                 edges, or other regions, or whatever.  Sometimes a node
//                 will be on the node list of more than one region/edge.
//
//                 To resolve such descrepancies, this function assumes that
//                 an edge label will take priority over a region label.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/29/02
//-----------------------------------------------------------------------------
bool Instance::setupLabelIndex ()
{

  // first, loop over the various non-edge(region) labels, and loop over their
  // node lists.  Set the labelIndex for each of these nodes to correspond
  // the current region label.

  int i;
  vector<int>::iterator firstL;
  vector<int>::iterator lastL;
  vector<int>::iterator iterL;

  for (i=0;i<numMeshLabels;++i)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(i);
    if (labelPtr->uType == TYPE_EDGE) continue;

    firstL = labelPtr->mNodeVector.begin();
    lastL  = labelPtr->mNodeVector.end ();

    for (iterL=firstL; iterL!=lastL; ++iterL)
    {
       labelIndex[*iterL] = i;
    }
  }

  // next, loop through the edge labels, and their node lists, and update
  // the various label indices.  If there are conflicts, the edge index will
  // just write over a previously established  region index.

  for (i=0;i<numMeshLabels;++i)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(i);
    if (labelPtr->uType != TYPE_EDGE) continue;

    firstL = labelPtr->mNodeVector.begin();
    lastL  = labelPtr->mNodeVector.end ();

    for (iterL=firstL; iterL!=lastL; ++iterL)
    {
       labelIndex[*iterL] = i;
    }
  }

  for (i=0;i<numMeshPoints;++i)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(labelIndex[i]);
    labelNameVector.push_back(labelPtr->name);
  }

  int size = dIVec.size();
  for (i = 0;i < size; ++i)
  {
    labelDIMap[dIVec[i].eName] = i;
  }

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << dashedline << endl;
    cout << "In the Intance::setupLabelIndex ";
    cout << "  name = "<< getName() <<endl;

#if 0
    vector<string>::iterator firstM = labelNameVector.begin();
    vector<string>::iterator lastM  = labelNameVector.end ();
    vector<string>::iterator iterM;

    for (i=0, iterM=firstM;iterM!=lastM;++i,++iterM)
    {
      cout << "  i = "<<i<<" name = " << *iterM << endl;
    }
#endif
    map<string,int>::iterator firstDIM = labelDIMap.begin ();
    map<string,int>::iterator lastDIM  = labelDIMap.end ();
    map<string,int>::iterator iterDIM;
    for(iterDIM=firstDIM;iterDIM!=lastDIM;++iterDIM)
    {
      cout << " DI index = "<< iterDIM->second;
      cout << "  name = " << iterDIM->first << endl;
    }
  }

  if (getDeviceOptions().debugLevel > 1)
  {
    for (i=0;i<numMeshPoints;++i)
    {
      mLabel * labelPtr = meshContainerPtr->getLabel(labelIndex[i]);
      cout << "  labelIndex["<<i<<"] = " << labelIndex[i];
      cout << "  name = " << labelPtr->name <<endl;
    }
  }

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << dashedline << endl;
  }

#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupBoundaryStencil
// Purpose       : This function sets up the stencil vector for boundary
//                 nodes.  In this case, a "boundary node" is a node from
//                 the mesh which is part of one of the user-specified
//                 boundary conditions.  So, this would include an
//                 electrode, like "COLLECTOR", but would not include
//                 "NOFLUX".
//
//                 If node i is a boundary node, then the value of
//                 boundarySten[i] = 1.  Otherwise boundarySten[i] = 0.
//
//                 7/15/03.  Revised to handle mixed boundary conditions,
//                 and the 3 new boundary stencils (one for each variable:
//                 (V,N,P)).
//
//                 Note: this only has nonzeros in it if using the "NEW_BC"
//                 version of boundary conditions.  At this point 7/17/03,
//                 the boundary stencils are only used for the new boundary
//                 conditions, not the old.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/07/02
//-----------------------------------------------------------------------------
bool Instance::setupBoundaryStencil ()
{
  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  vector<DeviceInterfaceNode>::iterator iterDI;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "In Instance::setupBoundaryStencil." << endl;
  }
#endif

#ifdef Xyce_NEW_BC
  for (iterDI = firstDI; iterDI!=lastDI; ++iterDI)
  {
    // loop over the nodes of this device interface node,
    // If it is an edge label, not a region label, and it is associated
    // with a boundary condition.

     if ( !( meshContainerPtr->labelEdgeType (iterDI->eName) ) ) continue;

     mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

     vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     vector<int>::iterator iterI;

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       int nodeIndex = *iterI;

       if (!(iterDI->neumannBCFlagV)) boundaryStenV[nodeIndex] = 1;
       if (!(iterDI->neumannBCFlagN)) boundaryStenN[nodeIndex] = 1;
       if (!(iterDI->neumannBCFlagP)) boundaryStenP[nodeIndex] = 1;

       // if this BC is dirichlet for all the variables, then set the
       // stencil.
       if (!(iterDI->neumannBCFlagV) &&
           !(iterDI->neumannBCFlagN) && !(iterDI->neumannBCFlagP))
       {
         boundarySten[nodeIndex] = 1;
       }
     }
  }
#endif // NEW_BC

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::checkForElectrodeOverlap
// Purpose       : The purpose of this function is to make sure that there
//                 are not any nodes associated with multiple electrodes.
//
//                 If there were, the boundary conditions wouldn't make sense.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/2004
//-----------------------------------------------------------------------------
bool Instance::checkForElectrodeOverlap ()
{
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "In Instance::checkForElectrodeOverlap." << endl;
  }
#endif

  for (int iDI=0;iDI<dIVec.size();++iDI)
  {
    // loop over the nodes of this device interface node,
    // If it is an edge label, not a region label, and it is associated
    // with a boundary condition.

     if ( !( meshContainerPtr->labelEdgeType (dIVec[iDI].eName) ) ) continue;

     mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[iDI].eName);

     vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     vector<int>::iterator iterI;

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       int nodeIndex = *iterI;

       if (boundaryTest[nodeIndex] != 0)
       {
         string msg = "Electrodes " + dIVec[iDI].eName + " and " +
       	 dIVec[boundaryTest[nodeIndex]-1].eName + " overlap";
         N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
       }

       boundaryTest[nodeIndex] = iDI+1;
     }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupNumVars
// Purpose       : The purpose of this function is to set up a few
//                 integer variables, such as numIntVars and numExtVars.
//                 These numbers are used by the registerGID , etc.
//                 functions.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/15/03
//-----------------------------------------------------------------------------
bool Instance::setupNumVars ()
{
  bool bsuccess = true;

  numIntVars    = 3*numMeshPoints;   // check this also.

  numExtVars    = numElectrodes; // This is the number of external nodes,
                                 // so it could be just about anything.

  numStateVars  = numElectrodes + numMeshPoints;
  //  numMeshPoints is part of numStateVars for displacement current.

  maxColsPerRow = 20;     // check this out later... depends on
                          // the max. NN count. (*3)

  int totalDirchlet = 0;

#ifdef Xyce_NEW_BC
  // For the new boundary conditions, reduce the size of the problem by
  // the number of mesh points along the boundary
  // (If all BC are dirichlet, then *3, for each equation).
  numInterfaceMeshPoints = 0;

  vector<DeviceInterfaceNode>::iterator first = dIVec.begin ();
  vector<DeviceInterfaceNode>::iterator last  = dIVec.end   ();
  vector<DeviceInterfaceNode>::iterator iterV;
  for (iterV=first;iterV!=last; ++iterV)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterV->eName);
    numInterfaceMeshPoints += labelPtr->mNodeVector.size();
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << iterV->eName;
      cout << ":  numInterfaceMeshPoints = ";
      cout << labelPtr->mNodeVector.size ();
      cout << endl;
    }
#endif
  }

  // revising, because of possibility of mixed BC.
  for (iterV=first;iterV!=last; ++iterV)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterV->eName);
    int numPoints = labelPtr->mNodeVector.size();

    int mult = 0;

    if (!(iterV->neumannBCFlagV)) mult += 1;
    if (!(iterV->neumannBCFlagN)) mult += 1;
    if (!(iterV->neumannBCFlagP)) mult += 1;

    totalDirchlet += (mult * numPoints);
  }

  numIntVars -= totalDirchlet;

#endif // Xyce_NEW_BC

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "\n";
    cout << " numInterfaceMeshPoints   = " << numInterfaceMeshPoints<<endl;
    cout << " numMeshPoints            = " << numMeshPoints<<endl;
    cout << " numElectrodes            = " << numElectrodes<<endl;
    cout << " numIntVars               = " << numIntVars<<endl;
    cout << " 3*numMeshPoints          = " << 3*numMeshPoints<<endl;
    cout << " 3*numInterfaceMeshPoints = "<<3*numInterfaceMeshPoints<<endl;
    cout << " totalDirchlet            = " << totalDirchlet<<endl;
    cout << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::allocatePDTerms.
//
// Purpose       : This function sets up and allocates a number of arrays
//                 that are needed by the function pdTerminalCurrents.
//
// Special Notes : Ordinarily, this function would have been called
//                 earlier, as I prefer to get allocations out of the way
//                 as early as possible.  However, it was easier to set
//                 this up using the colarrays, which are set up by
//                 function registerGIDs.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/10/02
//-----------------------------------------------------------------------------
bool Instance::allocatePDTerms ()
{
  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  vector<DeviceInterfaceNode>::iterator iterDI;
  // now do dIdX.
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

    // obtain the node indices for the current label, loop over them.
    // for each edge node, add an extra column entry to the colarray.

    vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
    vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
    vector<int>::iterator iterI;

    vector<EDGEINFO>::iterator firstEI;
    vector<EDGEINFO>::iterator lastEI;
    vector<EDGEINFO>::iterator iterEI;

    int cnt2  = 0;
    int nodeIndex;
    int col1;
    bool bmatch;
    int iVcol = 0;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // voltage variables:
      // do the center point first.
      col1 = iterDI->Vcol[iVcol];
      if (col1 != -1)
      {
        // check if node == any previous nodes in the cols array
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          iterDI->dIdXcols.push_back(col1);
          iterDI->dIdX.push_back(0.0);
          iterDI->dQdX.push_back(0.0);
        }
      }
      ++iVcol;
      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iVcol)
      {
        col1 = iterDI->Vcol[iVcol];
        if (col1 !=-1)
        {
          // check if node == any previous nodes in the cols array
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            iterDI->dIdXcols.push_back(col1);
            iterDI->dIdX.push_back(0.0);
            iterDI->dQdX.push_back(0.0);
          }
        }
      } // end of nn edge loop
    } // end of node loop

    int iNcol = 0;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // electron variables:
      // do the center point first.
      col1 = iterDI->Ncol[iNcol];
      if (col1 != -1)
      {
        // check if node == any previous nodes in the cols array
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          iterDI->dIdXcols.push_back(col1);
          iterDI->dIdX.push_back(0.0);
          iterDI->dQdX.push_back(0.0);
        }
      }
      ++iNcol;
      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iNcol)
      {
        col1 = iterDI->Ncol[iNcol];
        if (col1 !=-1)
        {
          // check if node == any previous nodes in the cols array
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            iterDI->dIdXcols.push_back(col1);
            iterDI->dIdX.push_back(0.0);
            iterDI->dQdX.push_back(0.0);
          }
        }
      } // end of nn edge loop
    } // end of node loop

    int iPcol = 0;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // hole variables:
      // do the center point first.
      col1 = iterDI->Pcol[iPcol];
      if (col1 != -1)
      {
        // check if node == any previous nodes in the cols array
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          iterDI->dIdXcols.push_back(col1);
          iterDI->dIdX.push_back(0.0);
          iterDI->dQdX.push_back(0.0);
        }
      }
      ++iPcol;
      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iPcol)
      {
        col1 = iterDI->Pcol[iPcol];
        if (col1 !=-1)
        {
          // check if node == any previous nodes in the cols array
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            iterDI->dIdXcols.push_back(col1);
            iterDI->dIdX.push_back(0.0);
            iterDI->dQdX.push_back(0.0);
          }
        }
      } // end of nn edge loop
    } // end of node loop

    // Now add to the neighbor node array.  Assuming that any of the 3
    // columns iPcol, iNcol, iVcol will do for a validity check vs. -1,
    // so just checking the pcol.
    iPcol = 0;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // center point first.
      col1 = iterDI->Pcol[iPcol];
      int meshNode = *iterI;

      if (col1 !=-1)
      {
        // check if node == any previous nodes in the cols array
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->neighborNodes.size();++cnt2)
        {
          if (iterDI->neighborNodes[cnt2] == meshNode)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          iterDI->neighborNodes.push_back(meshNode);
        }
      }

      ++iPcol;

      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iPcol)
      {
        col1 = iterDI->Pcol[iPcol];
        int meshNode = iterEI->inode;

        if (col1 !=-1)
        {
          // check if node == any previous nodes in the cols array
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->neighborNodes.size();++cnt2)
          {
            if (iterDI->neighborNodes[cnt2] == meshNode)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            iterDI->neighborNodes.push_back(meshNode);
          }
        }
      } // end of nn edge loop
    } // end of node loop

    int size1 = iterDI->neighborNodes.size();
    int size3 = iterDI->dIdX.size();
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << endl;
      cout << "number of neighbor nodes for " << iterDI->eName;
      cout << " is " << size1 << endl;
      int  i;
      for (i=0;i<size1;++i)
      {
        cout << "neighborNodes["<<i<<"] = " << iterDI->neighborNodes[i] << endl;
      }
      cout << endl;
      cout << "dIdX size for " << iterDI->eName << "  is " << size3 << endl;
      for (i=0;i<size3;++i)
      {
        cout << "dIdX["<<i<<"] = " << iterDI->dIdXcols[i] << endl;
      }
    }
#endif

#ifdef Xyce_NEW_BC
    // to set up dFdVckt, need to take  it through the same set of loops it
    // will be subject to in function pdTerminalCurrents.
    int numNeighbor = iterDI->neighborNodes.size();
    int iNeighbor;
    int dFdVindex = 0;
    for (iNeighbor=0;iNeighbor<numNeighbor;++iNeighbor)
    {
      int inode = iterDI->neighborNodes[iNeighbor];
      mNode * nodePtr = meshContainerPtr->getNode(inode);

      for (int iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        int inodeB = nodePtr->edgeInfoVector[iNN].inode;

        // if nodeB is not a boundary node, never mind.
        if (boundaryStenV[inodeB]!=1) continue;

        // if it is a boundary node, but not part of the
        // current boundary, also never mind.
        if (labelNameVector[inodeB]!= iterDI->eName) continue;

        iterDI->dFdVckt.push_back(0.0); ++dFdVindex;
        iterDI->dFdVckt.push_back(0.0); ++dFdVindex;
        iterDI->dFdVckt.push_back(0.0); ++dFdVindex;
      }
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << iterDI->eName << ": size of dFdVckt = " << dFdVindex << endl;
    }
#endif

#else
    string msg = "NEW_BC must be turned on to use TWO_LEVEL_NEWTON";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_WARNING_0,msg);
#endif // Xyce_NEW_BC

  }  // end of DI loop

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupBCEdgeAreas
//
// Purpose       : This function sets up the areas associated with each
//                 mesh node along a boundary condition edge.
//
// Special Notes : This is potentially tricky.  The mesh may be 2D
//                 cartesian or 2D cylindrical.
//
//                 Fortunately, the function "lengthAdjust", from the mesh
//                 class, is designed for getting areas - either
//                 cylindrical or cartesian.
//
//                 This whole should be moved to the mesh class later.
//
//                 This function is called before scaling is turned on, but
//                 should work either way.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/09/02
//-----------------------------------------------------------------------------
bool Instance::setupBCEdgeAreas ()
{

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline="--------------------------------------------------"
    "---------------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline << "\n";
    cout << "setupBCEdgeAreas.  name = " << getName() << endl;
    cout.setf(ios::scientific);
  }
#endif

  // now set up the local areas for needed to interface device boundary
  // conditions to the circuit

  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI = firstDI; iterDI!=lastDI; ++iterDI)
  {
    // loop over the nodes of this device interface node:

     if ( !( meshContainerPtr->labelEdgeType (iterDI->eName) ) ) continue;

     mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

     vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     vector<int>::iterator iterI;

     iterDI->area       = 0.0;  // total area for the edge.

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       // loop over neighbor nodes/edges to get area sum for this node.
       mNode * nodePtr = meshContainerPtr->getNode(*iterI);

       vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
       vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end  ();
       vector<EDGEINFO>::iterator iterEI;

       double areaLocal = 0.0; // total "area" for the this node
       double areaTmp   = 0.0; // partial area for the this node(from one edge)

#ifdef Xyce_DEBUG_DEVICE
       if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
       {
         cout << " --------------- " << endl;
         cout << "name = " << iterDI->eName;
         cout << "  node      = " << *iterI <<endl;
       }
#endif

       for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
       {
         int neighbor = iterEI->inode;

         // if this edge is actually on the boundary, then sum the
         // edge length into the "area" associated with this
         // boundary node.  Check this by checking the label index
         // of the neighbor.

         int ilabel = labelIndex[neighbor];
         mLabel * labelPtr = meshContainerPtr->getLabel(ilabel);

         areaTmp = 0.0;
         if (labelPtr->name == iterDI->eName)
         // if this is along the bounadry
         {
           if (meshContainerPtr->cylGeom)
           {
             // get the midpoint location, x2:
             double x1 = xVec[*iterI];
             double y1 = yVec[*iterI];

                   double x2 = xVec[neighbor];
                   double y2 = yVec[neighbor];

             double dx = x2-x1;   dx *= 0.5;
             double dy = y2-y1;   dy *= 0.5;

             x2 = x1 + dx;
             y2 = y1 + dy;

             areaTmp = meshContainerPtr->lengthAdjust(x1,y1,x2,y2);
             areaLocal += areaTmp;
           }
           else // cartesian
           {
             areaTmp    = 0.5 * iterEI->elen;
             areaLocal += 0.5 * iterEI->elen;
           }
         }
#ifdef Xyce_DEBUG_DEVICE
         if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
         {
           cout << "  neighbor node   = " << neighbor <<endl;
           cout << "  areaTmp         = " << areaTmp <<endl;
           cout << "  areaLocal       = " << areaLocal << endl;
           cout << "  elen            = " << iterEI->elen << endl;
           cout << "  label name      = " << labelPtr->name << endl;
           cout << "  DI eName        = " << iterDI->eName << endl;
           cout << "---" << endl;
         }
#endif
       }

#ifdef Xyce_DEBUG_DEVICE
       if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
       {
         cout << " --------------- " << endl;
       }
#endif
       iterDI->area += areaLocal;
       iterDI->areaVector.push_back(areaLocal);

     } // iterI loop.


#ifdef Xyce_DEBUG_DEVICE
     if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout.setf(ios::scientific);
      cout << "  Total area for edge: " << iterDI->area << endl;
    }
#endif
  }  // iterDI loop.

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline << "\n";
  }
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMinDXVector
// Purpose       : This finds the minimum edge length.  I'm not sure if
//                 there's any reason to do this anymore.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/26/02
//-----------------------------------------------------------------------------
bool Instance::setupMinDXVector ()
{

  int i;
  for (i=0;i<numMeshPoints;++i)
  {
    // loop over neighbor nodes/edges to get area sum for this node.
    mNode * nodePtr = meshContainerPtr->getNode(i);

    vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
    vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end  ();
    vector<EDGEINFO>::iterator iterEI;

    double tmpDX = +1.0e99;
    for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
    {
      if (tmpDX > iterEI->elen) tmpDX = iterEI->elen;
    }
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2)
    {
      cout << "i = " << i << "   minDX = " << tmpDX; cout << endl;
    }
#endif
    minDXVec[i] = tmpDX;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupJacStamp
//
// Purpose       : This function sets up the jacobian stamp data structure,
//                 which is a symbolic, local-ID (LID) representation of the
//                 jacobian matrix structure of the device.
//
//                 The jacStamp is the structure that gets passed up to
//                 topology.
//
//                 It is similar to some of what has to happen in the
//                 registerGID's function, but it is not dependent upon
//                 knowing the global ID's for the solution variables.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/23/03
//-----------------------------------------------------------------------------
bool Instance::setupJacStamp ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "In Instance::setupJacStamp" << endl;
    cout << "numIntVars = " << numIntVars << endl;
    cout << "numMeshPoints = " << numMeshPoints << endl;
  }
#endif

  // The number of rows in the jacStamp has to include all the *possible*
  // electrodes (external variables). This is a little confusing.
  // The PDE devices are set up so that there are:
  //   2 required nodes
  //   2 not-required, but not-optional nodes (these are fill nodes)
  //   100 optional nodes.
  //   For some reason (that I don't understand) this means that the
  //   jacStamp always has to have at least 4 external nodes.  Even
  //   if the device only has 2 nodes attached to the circuit, the
  //   jacStamp is required to have 4.
  int numPossibleElectrodes = (numExtVars>4)?numExtVars:4;
  jacStamp.resize(numIntVars + numPossibleElectrodes);

  MESHtoLID_V.resize(numMeshPoints,-1);
  MESHtoLID_N.resize(numMeshPoints,-1);
  MESHtoLID_P.resize(numMeshPoints,-1);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "size of jacStamp = " << jacStamp.size() << endl;
  }
#endif

  // Set up the MESHtoLID converter for internal vars.
  int LIDcounter = numPossibleElectrodes;

  for (i=0;i<numMeshPoints;++i)
  {
#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

    if (!(boundaryStenV[i]))
    {
      MESHtoLID_V[i] = LIDcounter; ++LIDcounter;
    }

    if (!(boundaryStenN[i]))
    {
      MESHtoLID_N[i] = LIDcounter; ++LIDcounter;
    }

    if (!(boundaryStenP[i]))
    {
      MESHtoLID_P[i] = LIDcounter; ++LIDcounter;
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // Do the external variables first:
  // Loop over the boundary condition "device interface" labels.
  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  vector<DeviceInterfaceNode>::iterator iterDI;
  int DIsize = dIVec.size ();
  firstDI = dIVec.begin();
  lastDI  = dIVec.end  ();

  // note that the index, in addition to being the index into the
  // array of "DeviceInterfaceNodes" will also be the LID.
  int index;
  for(index=0,iterDI=firstDI;iterDI!=lastDI;++iterDI,++index)
  {
    // check that this label exists, and is an edge label
    // (everything in the dIVec container should pass these tests
    //   by this point - these two if-statements are just a failsafe.)
     if ( !( meshContainerPtr->labelNameExist(dIVec[index].eName) ) ) continue;
     if ( !( meshContainerPtr->labelEdgeType (dIVec[index].eName) ) ) continue;


     // first do the external node's (row=cktnode,col=cktnode) pair:
     jacStamp[index].push_back(index);

     dIVec[index].numCrossTerms = 0;

     // next do the (row,col) pairs for other ckt nodes.  This is needed for
     // 2-level newton only.
     for (int ind2=0;ind2<DIsize;++ind2)
     {
       if (ind2==index) continue;

       jacStamp[index].push_back(ind2);
       ++(dIVec[index].numCrossTerms);
     }

     mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[index].eName);

     // obtain the node indices for the current label, loop over them.
     // for each edge node, add an extra column entry to the colarray,
     // which corresponds to the global ID of the connected ckt node.

     vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     vector<int>::iterator iterI;

     vector<EDGEINFO>::iterator firstEI;
     vector<EDGEINFO>::iterator lastEI;
     vector<EDGEINFO>::iterator iterEI;

     int imesh;

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       // now do the V edge and edge neighbor nodes.
       imesh = *iterI;
       if (!boundaryStenV[imesh])
         jacStamp[index].push_back(MESHtoLID_V[imesh]);

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();

       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         imesh = iterEI->inode;
         if (!boundaryStenV[imesh])
           jacStamp[index].push_back(MESHtoLID_V[imesh]);
       }
     }

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       // now do the N nodes.
       imesh = *iterI;

       if (!boundaryStenN[imesh])
         jacStamp[index].push_back(MESHtoLID_N[imesh]);

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();
       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         imesh = iterEI->inode;
         if (!boundaryStenN[imesh])
           jacStamp[index].push_back(MESHtoLID_N[imesh]);
       }
     }

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       // now do the P nodes.
       imesh = *iterI;
       if (!boundaryStenP[imesh])
         jacStamp[index].push_back(MESHtoLID_P[imesh]);

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();
       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         imesh = iterEI->inode;
         if (!boundaryStenP[imesh])
           jacStamp[index].push_back(MESHtoLID_P[imesh]);
       }
     }// iterI loop
  } // index (dIVec) loop.

  /////////////////////////////////////////////////////////////////////////
  // Now do the internal variables.  (V,N,P on the mesh)
  for(i=0;i<numMeshPoints;++i)
  {
#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << "  mesh point i = " << i << endl;
    }
#endif

    mNode * nodePtr = meshContainerPtr->getNode(i);
    vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin ();
    vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end   ();
    vector<EDGEINFO>::iterator iterEI;

    // get the temporary LID row indices:
    int Vrow = MESHtoLID_V[i];
    int Nrow = MESHtoLID_N[i];
    int Prow = MESHtoLID_P[i];

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << "  Vrow = " << Vrow << endl;
      cout << "  Nrow = " << Nrow << endl;
      cout << "  Prow = " << Prow << endl;
    }
#endif

    // voltage row:
    if (Vrow != -1)
    {
      // First do the (row, row) pair.
      jacStamp[Vrow].push_back(Vrow);

      // loop  over the neighbor nodes of node i.
      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
#ifdef Xyce_NEW_BC
        // if the neighbor node is on a boundary, need to use the
        // GID of the connected ckt node, rather than the (nonexistant)
        // GID of the boundary mesh node.
        if (boundaryStenV[iterEI->inode])
        {
          // get the gid:
          int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
          jacStamp[Vrow].push_back (DIindex);
        }
        else
        {
          int imesh = iterEI->inode;
          int lid = MESHtoLID_V[imesh];
          jacStamp[Vrow].push_back(lid);
        }
#else // Xyce_NEW_BC
        int imesh = iterEI->inode;
        int lid = MESHtoLID_V[imesh];
        jacStamp[Vrow].push_back(lid);
#endif // Xyce_NEW_BC
      }

      jacStamp[Vrow].push_back(Nrow);
      jacStamp[Vrow].push_back(Prow);
    }

    // electron row:
    if (Nrow != -1)
    {
      jacStamp[Nrow].push_back(Nrow);

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        int imesh = iterEI->inode;
        if (!boundaryStenN[imesh])
        {
          int lid = MESHtoLID_N[imesh];
          jacStamp[Nrow].push_back(lid);
        }
      }

      jacStamp[Nrow].push_back(Vrow);

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
#ifdef Xyce_NEW_BC
        // if the neighbor node is on a boundary, need to use the
        // GID of the connected ckt node, rather than the (nonexistant)
        // GID of the boundary mesh node.
        if (boundaryStenV[iterEI->inode])
        {
          // get the gid:
          int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
          jacStamp[Nrow].push_back (DIindex);
        }
        else
        {
          int imesh = iterEI->inode;
          int lid = MESHtoLID_V[imesh];
          jacStamp[Nrow].push_back(lid);
        }
#else // Xyce_NEW_BC
        int imesh = iterEI->inode;
        int lid = MESHtoLID_V[imesh];
        jacStamp[Nrow].push_back(lid);
#endif // Xyce_NEW_BC
      }
      jacStamp[Nrow].push_back(Prow);
    }

    // hole row:
    if (Prow != -1)
    {
      jacStamp[Prow].push_back(Prow);

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        int imesh = iterEI->inode;
        if(!boundaryStenP[imesh])
        {
          int lid = MESHtoLID_P[imesh];
          jacStamp[Prow].push_back(lid);
        }
      }
      jacStamp[Prow].push_back(Vrow);

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
#ifdef Xyce_NEW_BC
        // if the neighbor node is on a boundary, need to use the
        // GID of the connected ckt node, rather than the (nonexistant)
        // GID of the boundary mesh node.
        if (boundaryStenV[iterEI->inode])
        {
          // get the gid:
          int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
          jacStamp[Prow].push_back (DIindex);
        }
        else
        {
          int imesh = iterEI->inode;
          int lid = MESHtoLID_V[imesh];
          jacStamp[Prow].push_back(lid);
        }
#else // Xyce_NEW_BC
        int imesh = iterEI->inode;
        int lid = MESHtoLID_V[imesh];
        jacStamp[Prow].push_back(lid);
#endif // Xyce_NEW_BC
      }
      jacStamp[Prow].push_back(Nrow);
    }
  } // mesh point loop.

#ifndef Xyce_NEW_BC
  // go back to the nodes which are boundary condition nodes - ones
  // that are connected to the external nodes, and add a few
  // (row, col) pairs .  These nodes will need some extra ones, at
  // least to  handle boundary conditions on the voltage.

  // loop over the boundary condition "device interface" labels.
  firstDI = dIVec.begin();
  lastDI  = dIVec.end  ();

  for(index=0,iterDI=firstDI;iterDI!=lastDI;++iterDI,++index)
  {
    // check that this label exists, and is an edge label
    // (everything in the dIVec container should pass these tests
    //   by this point)
     if ( !( meshContainerPtr->labelNameExist(dIVec[index].eName) ) ) continue;
     if ( !( meshContainerPtr->labelEdgeType (dIVec[index].eName) ) ) continue;

     mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[index].eName);

     // obtain the node indices for the current label, loop over them.
     // for each edge node, add an extra column entry to the colarray,
     // which corresponds to the global ID of the connected ckt node.

     vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     vector<int>::iterator iterI;

     // For the "new" boundary conditions, this step isn't neccessary,
     // as the boundary nodes no longer have any dependent variables.
     // However, it is neccessary for the col array associated with
     // the "boundary condition" equation.
     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       int imesh = *iterI;
       int Vrow = MESHtoLID_V[imesh];
       jacStamp[Vrow].push_back(index);
     }
  }  // end of index (dIVec) loop
#endif // Xyce_NEW_BC

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1)
  {
    int irow, icol;
    for(irow=0;irow<jacStamp.size();++irow)
    {
      cout << "irow = " << irow;
      if (irow < dIVec.size())
        cout << "  " << dIVec[irow].eName << "  KCL";
      cout << endl;
      for (icol=0;icol<jacStamp[irow].size();++icol)
      {
        int jsTmp = jacStamp[irow][icol];
        cout << "   jacStamp["<<irow<<"]["<<icol<<"] = "<<jsTmp;
        cout << endl;
      }
    }
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
void Instance::registerGIDs(
  const list<index_pair> & intGIDListRef,
  const list<index_pair> & extGIDListRef)
{
  string msg;
  string tmpstr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << dashedline << endl;
    cout << "In the Intance::registerGIDs function.  ";
    cout << "  name = "<< getName() <<endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intGIDListRef.size();
  int numExt = extGIDListRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "numInt = " << numInt <<"\n";
    cout << "numExt = " << numExt <<"\n";
    cout << "numMeshPoints = " << numMeshPoints << "\n";
  }
#endif

  // number of internal variables equals the number of
  // mesh points *3.
  if (numInt != numIntVars)
  {
    msg = "Instance::registerGIDs:";
    msg += "numInt != 3*numMeshPoints.  ";
    msg += "Check the metadata in the parser\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  // copy over the global ID lists:
  intGIDList.assign(intGIDListRef.begin(), intGIDListRef.end());
  extGIDList.assign(extGIDListRef.begin(), extGIDListRef.end());

  // Set up the rows first, then do the colarrays.

  // First do the external variables:
  // These will all be voltage nodes connected to the devices.
  list<index_pair>::iterator it1 = extGIDList.begin();

  list<index_pair>::iterator first = extGIDList.begin();
  list<index_pair>::iterator last  = extGIDList.end  ();
  list<index_pair>::iterator iter;
  int index;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Setting up the indices for the external circuit nodes:" << endl;
    cout << "External node list:" << endl;
  }
#endif

  int numExtNodes = 0;
  for(iter=first;iter!=last;++iter)
  {
     if ( (iter->col) && (iter->row != -1) )
     {
       ++numExtNodes;
     }

#ifdef Xyce_DEBUG_DEVICE
     if (getDeviceOptions().debugLevel > 0)
    {
      cout << "node = " << iter->row << " col = ";
      cout << iter->col << "  numExtNodes = " << numExtNodes << endl;
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
    cout << " number of owned external nodes = " << numExtNodes << endl;
#endif

  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  vector<DeviceInterfaceNode>::iterator iterDI;
  int numRealDIs = 0;
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    if (iterDI->given) ++numRealDIs;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << " number of user specified boundary conditions = ";
    cout << numRealDIs << endl;
  }
#endif

  if (numRealDIs < numExtNodes)
  {
    msg = "Instance::registerGIDs:";
    msg += "number of boundary conditions < number of external nodes";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  int isizeDI = dIVec.size();

  for(index=0, iter=first; (iter!=last && index < isizeDI); ++iter, ++index)
  {
    dIVec[index].gid = iter->row;
    if ( !(iter->col) ) dIVec[index].gid = -1;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << "   name = "<<dIVec[index].eName<<" gid = ";
      cout << dIVec[index].gid;
      cout << endl;
    }
#endif

    // Obtain the node index for the first node of the current label.
    // Note: This particular attribute may be obsolete...
    mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[index].eName);
    vector<int>::iterator iterTmp = labelPtr->mNodeVector.begin();

    dIVec[index].firstMeshNodeIndex = *iterTmp;
  }


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "doing internal vars, row arrays:"<<endl;
  }
#endif

  // Do the internal variables.  There should be a lot of these.
  list<index_pair>::iterator intEnd = intGIDList.end();
  list<index_pair>::iterator it2    = intGIDList.begin();
  int i=0;

  // The interior points will be blocked (V,N,P) together.
  while (it2 != intEnd && i < numMeshPoints)
  {
    int Vrow = -1;
    int Nrow = -1;
    int Prow = -1;

#ifdef Xyce_NEW_BC
    // don't set up row, col arrays if this is a boundary node.
    if (boundarySten[i]) { ++i; continue; }

    if (!(boundaryStenV[i]))
    {
      Vrow = it2->row;
      if ( !(it2->col) ) Vrow = -1;
      Vrowarray[i] = Vrow;   // should be Vrowarray[i]
      ++it2;
    }

    if (!(boundaryStenN[i]))
    {
      Nrow = it2->row;
      if ( !(it2->col) ) Nrow = -1;
      Nrowarray[i] = Nrow;   // should be Nrowarray[i]
      ++it2;
    }

    if (!(boundaryStenP[i]))
    {
      Prow = it2->row;
      if ( !(it2->col) ) Prow = -1;
      Prowarray[i] = Prow;   // should be Prowarray[i]
      ++it2;
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel >1)
    {
      cout << "doing row arrays for mesh point " << i << endl;
      if (!(boundaryStenV[i]))
	cout << "  Vrow = " << Vrow << endl;

      if (!(boundaryStenN[i]))
	cout << "  Nrow = " << Nrow << endl;

      if (!(boundaryStenP[i]))
	cout << "  Prow = " << Prow << endl;
    }
#endif

#else

    Vrow = it2->row;
    if ( !(it2->col) ) Vrow = -1;
    Vrowarray[i] = Vrow;   // should be Vrowarray[i]
    ++it2;

    Nrow = it2->row;
    if ( !(it2->col) ) Nrow = -1;
    Nrowarray[i] = Nrow;   // should be Nrowarray[i]
    ++it2;

    Prow = it2->row;
    if ( !(it2->col) ) Prow = -1;
    Prowarray[i] = Prow;   // should be Prowarray[i]
    ++it2;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel >1)
    {
      cout << "doing row arrays for mesh point " << i << endl;
      cout << "  Vrow = " << Vrow << endl;
      cout << "  Nrow = " << Nrow << endl;
      cout << "  Prow = " << Prow << endl;
    }
#endif
#endif // Xyce_NEW_BC

    ++i;
  }

  setupRowColPairs ();

  // Make sure the cols and vals arrays are big enough,
  // based on maxColsPerRow. (the variable maxColsPerRow was
  // setup in function setupRowColPairs).
  getMatrixLoadData().initializeAll(maxColsPerRow);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "maxColsPerRow = " << maxColsPerRow <<endl;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  cout << dashedline << endl;
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupIntNameMap
// Purpose       : Sets up the "internal names map".
//
// Special Notes : This is to help with debugging.  It assigns a unique
//                 string name to each internal variable of this device.
//
//                 This function was once part of registerGID's, but I
//                 decided to break that function up to make it more
//                 readable.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/03
//-----------------------------------------------------------------------------
void Instance::setupIntNameMap ()
{
  int i;
  char tmpchar[128]; for (i=0;i<128;++i) tmpchar[i] = 0;

  // set up the intNameMap, if necessary.
  for (i=0;i<numMeshPoints;++i)
  {
#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

    int Vrow, Nrow, Prow;
    Vrow = Vrowarray[i]; Nrow = Nrowarray[i]; Prow = Prowarray[i];

    if (Vrow != -1)
    { sprintf(tmpchar,"%s_V_%d_%s",
              getName().c_str(), i, labelNameVector[i].c_str());
      intNameMap[Vrow] = string(tmpchar);
    }

    if (Nrow != -1)
    { sprintf(tmpchar,"%s_N_%d_%s",
              getName().c_str(), i, labelNameVector[i].c_str());
      intNameMap[Nrow] = string(tmpchar);
    }

    if (Prow != -1)
    { sprintf(tmpchar,"%s_P_%d_%s",
              getName().c_str(), i, labelNameVector[i].c_str());
      intNameMap[Prow] = string(tmpchar);
    }
  }

  return;
}
//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    setupIntNameMap ();
  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupRowColPairs
//
// Purpose       : This function performs part of the "registerGIDs"
//                 functionality, in that it sets up the Jacobian
//                 matrix (row,col) pairs.
//
// Special Notes : This function was once part of registerGIDs, but I
//                 moved it out b/c that function was so big that it needed
//                 to be broken up. This function is now called from
//                 registerGIDs.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/03
//-----------------------------------------------------------------------------
void Instance::setupRowColPairs ()
{

  int i, j;
  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  vector<DeviceInterfaceNode>::iterator iterDI;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Instance::setupRowColPairs ()" << endl;
    cout << "  doing internal vars, col arrays:"<<endl;
  }
#endif

  /////////////////////////////////////////////////////////////////////////
  // now do the colarrays for each variable type:
  for(i=0;i<numMeshPoints;++i)
  {

#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1)
    {
      cout << "doing col arrays for mesh point " << i<<endl;
    }
#endif

    mNode * nodePtr = meshContainerPtr->getNode(i);
    vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin ();
    vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end   ();
    vector<EDGEINFO>::iterator iterEI;

    // voltage col arrays:
    Vcolarray[i].push_back(Vrowarray[i]);

    // loop  over the neighbor nodes of node i.
    for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
    {
#ifdef Xyce_NEW_BC
      // if the neighbor node is on a boundary, need to use the
      // GID of the connected ckt node, rather than the (nonexistant)
      // GID of the boundary mesh node.
      if (boundaryStenV[iterEI->inode])
      {
        // get the gid:
        int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
        Vcolarray[i].push_back (dIVec[DIindex].gid);
      }
      else
      {
        Vcolarray[i].push_back(Vrowarray[iterEI->inode]);
      }
#else
      Vcolarray[i].push_back(Vrowarray[iterEI->inode]);
#endif // Xyce_NEW_BC
    }

    Vcolarray[i].push_back(Nrowarray[i]);
    Vcolarray[i].push_back(Prowarray[i]);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2)
    {
      cout << "size of vcolarray["<<i<<"] = " << Vcolarray[i].size() <<endl;
      for(int eric=0;eric<Vcolarray[i].size();++eric)
        cout << "  col["<<eric<<"] = " << Vcolarray[i][eric] <<endl;
    }
#endif

    // electron col arrays:
    Ncolarray[i].push_back(Nrowarray[i]);

    for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
    {
      Ncolarray[i].push_back(Nrowarray[iterEI->inode]);
    }

    Ncolarray[i].push_back(Vrowarray[i]);
    for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
    {
#ifdef Xyce_NEW_BC
      // if the neighbor node is on a boundary, need to use the
      // GID of the connected ckt node, rather than the (nonexistant)
      // GID of the boundary mesh node.
      if (boundaryStenV[iterEI->inode])
      {
	// get the gid:
        int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
        Ncolarray[i].push_back (dIVec[DIindex].gid);
      }
      else
      {
        Ncolarray[i].push_back(Vrowarray[iterEI->inode]);
      }
#else
      Ncolarray[i].push_back(Vrowarray[iterEI->inode]);
#endif // Xyce_NEW_BC
    }
    Ncolarray[i].push_back(Prowarray[i]);

    // hole col arrays:
    Pcolarray[i].push_back(Prowarray[i]);

    for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
    {
      Pcolarray[i].push_back(Prowarray[iterEI->inode]);
    }
    Pcolarray[i].push_back(Vrowarray[i]);
    for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
    {
#ifdef Xyce_NEW_BC
      // if the neighbor node is on a boundary, need to use the
      // GID of the connected ckt node, rather than the (nonexistant)
      // GID of the boundary mesh node.
      if (boundaryStenV[iterEI->inode])
      {
	// get the gid:
        int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
        Pcolarray[i].push_back (dIVec[DIindex].gid);
      }
      else
      {
        Pcolarray[i].push_back(Vrowarray[iterEI->inode]);
      }
#else
      Pcolarray[i].push_back(Vrowarray[iterEI->inode]);
#endif // Xyce_NEW_BC
    }
    Pcolarray[i].push_back(Nrowarray[i]);

  } // mesh point loop.

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "doing some boundary condition vars:"<<endl;
  }
#endif

  // go back to the nodes which are boundary condition nodes - ones
  // that are connected to the external nodes, and add a few
  // (row, col) pairs .  These nodes will need some extra ones, at
  // least to  handle boundary conditions on the voltage.

  // loop over the boundary condition "device interface" labels.
  firstDI = dIVec.begin();
  lastDI  = dIVec.end  ();

  int index;

  for(index=0,iterDI=firstDI;iterDI!=lastDI;++iterDI,++index)
  {
    // check that this label exists, and is an edge label
    // (everything in the dIVec container should pass these tests
    //   by this point)
     if ( !( meshContainerPtr->labelNameExist(dIVec[index].eName) ) ) continue;
     if ( !( meshContainerPtr->labelEdgeType (dIVec[index].eName) ) ) continue;

#ifdef Xyce_DEBUG_DEVICE
     if (getDeviceOptions().debugLevel > 1)
    {
      cout << "Device Interface: " << dIVec[index].eName << endl;
    }
#endif

     mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[index].eName);

     // obtain the node indices for the current label, loop over them.
     // for each edge node, add an extra column entry to the colarray,
     // which corresponds to the global ID of the connected ckt node.

     vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     vector<int>::iterator iterI;

     // For the "new" boundary conditions, this step isn't neccessary,
     // as the boundary nodes no longer have any dependent variables.
     // However, it is neccessary for the col array associated with
     // the "boundary condition" equation.
#ifndef Xyce_NEW_BC
     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       Vcolarray[*iterI].push_back(dIVec[index].gid);
     }
#endif // Xyce_NEW_BC

     int itmp;
     // Add indices to the colarray of the external circuit nodes.

     // These are for satisfying the KCL of the node - it will depend on
     // all the currents of all the edge nodes.  Unfortunately, this
     // results in a potentially(probably) dense matrix row.
     //
     // Note that if this is a dielectric boundary, then no current is
     // going to come out of device at this boundary (no conduction
     // current, anyway), and that a lot of these (row, col) pairs will
     // wind up being loaded with 0's in that case.  However, it doesn't
     // hurt to have them, so I'm leaving them in.

     // first push back the (gid,gid) entry.  Then do the entries related
     // to the  nearest neighbors of edge nodes.
     //dIVec[index].col.push_back(dIVec[index].gid);

     vector<EDGEINFO>::iterator firstEI;
     vector<EDGEINFO>::iterator lastEI;
     vector<EDGEINFO>::iterator iterEI;

#ifdef Xyce_DEBUG_DEVICE
     if (getDeviceOptions().debugLevel > 2)
     {
       cout << "V edge and edge neighbor gids:" <<endl;
     }
#endif

     // now do the V edge and edge neighbor nodes.
     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       itmp = Vrowarray[*iterI];
       dIVec[index].Vcol.push_back(itmp);

#ifdef Xyce_DEBUG_DEVICE
       if (getDeviceOptions().debugLevel > 2)
       {
         int ind1 = dIVec[index].Vcol.size()-1;
         cout << "  1Vcol["<<ind1<<"] = " << itmp << endl;
       }
#endif

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();

       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         itmp = Vrowarray[iterEI->inode];
         dIVec[index].Vcol.push_back(itmp);

#ifdef Xyce_DEBUG_DEVICE
         if (getDeviceOptions().debugLevel > 2)
         {
           int ind1 = dIVec[index].Vcol.size()-1;
           cout << "  2Vcol["<<ind1<<"] = " << itmp << endl;
         }
#endif
       }
     } // iterI loop

#ifdef Xyce_DEBUG_DEVICE
     if (getDeviceOptions().debugLevel > 2)
     {
       cout << "N edge and edge neighbor gids:" <<endl;
     }
#endif

     // now do the N nodes.
     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       itmp = Nrowarray[*iterI];
       dIVec[index].Ncol.push_back(itmp);

#ifdef Xyce_DEBUG_DEVICE
       if (getDeviceOptions().debugLevel > 2)
       {
         int ind1 = dIVec[index].Ncol.size()-1;
         cout << " 1Ncol["<<ind1<<"] = " << itmp << endl;
       }
#endif

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();

       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         itmp = Nrowarray[iterEI->inode];
         dIVec[index].Ncol.push_back(itmp);

#ifdef Xyce_DEBUG_DEVICE
         if (getDeviceOptions().debugLevel > 2)
         {
           int ind1 = dIVec[index].Ncol.size()-1;
           cout << " 2Ncol["<<ind1<<"] = " << itmp << endl;
         }
#endif
       }
     } // iterI loop

#ifdef Xyce_DEBUG_DEVICE
     if (getDeviceOptions().debugLevel > 2)
     {
       cout << "P edge and edge neighbor gids:" <<endl;
     }
#endif
     // now do the P nodes.
     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       itmp = Prowarray[*iterI];
       dIVec[index].Pcol.push_back(itmp);

#ifdef Xyce_DEBUG_DEVICE
       if (getDeviceOptions().debugLevel > 2)
         {
           int ind1 = dIVec[index].Pcol.size()-1;
           cout << " 1Pcol["<<ind1<<"] = " << itmp << endl;
         }
#endif

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();

       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         itmp = Prowarray[iterEI->inode];
         dIVec[index].Pcol.push_back(itmp);

#ifdef Xyce_DEBUG_DEVICE
         if (getDeviceOptions().debugLevel > 2)
         {
           int ind1 = dIVec[index].Pcol.size()-1;
           cout << " 2Pcol["<<ind1<<"] = " << itmp << endl;
         }
#endif
       }
     }// iterI loop

     if (maxColsPerRow <
                (dIVec[index].Vcol.size() +
                 dIVec[index].Ncol.size() +
                 dIVec[index].Pcol.size() + 10)
        )
     {
         maxColsPerRow =
                (dIVec[index].Vcol.size() +
                 dIVec[index].Ncol.size() +
                 dIVec[index].Pcol.size() + 10); // extra 10 just in case.
     }
  }  // end of index (dIVec) loop

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Setting up the index pairs."<<endl;
  }
#endif

  // Now that all the row and col arrays are set up, do the
  // index pair list.
  index_pair IP(0,0);

  int nn;

  // external KCL points first:
  firstDI = dIVec.begin();
  lastDI  = dIVec.end  ();
  for(index=0,iterDI=firstDI;iterDI!=lastDI;++iterDI,++index)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2)
      cout << "device interface name: " << iterDI->eName << endl;
#endif
    // check that this label exists, and is an edge label
    // (everything in the dIVec container should pass these tests
    //   by this point)
    if ( !( meshContainerPtr->labelNameExist(dIVec[index].eName) ) ) continue;
    if ( !( meshContainerPtr->labelEdgeType (dIVec[index].eName) ) ) continue;
    if ( dIVec[index].gid == -1) continue;

    mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[index].eName);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2)
      cout << "index = " << index<< "  gid = " << dIVec[index].gid <<endl;
#endif

    // Do the (gid, gid) node.
    IP.row = dIVec[index].gid;
    IP.col = dIVec[index].gid;
    indexPairList.push_back(IP);

    // do the other (gid_1, gid_2) elements.  These are needed for 2-level
    // "ckt-only" matrix loads.
    int size2 = dIVec.size();
    for (int ind2=0;ind2<size2;++ind2)
    {
      if (index == ind2) continue;
      IP.row = dIVec[index].gid;
      IP.col = dIVec[ind2].gid;
      indexPairList.push_back(IP);
    }

    // Loop over the mesh nodes that affect this terminal current.
    // Add row, col pairs for each variable.
    nn = dIVec[index].Vcol.size();
    for (j=0;j<nn;++j)
    {
      // electrostatic potential (row, col) pairs.
      if (dIVec[index].Vcol[j] != -1)
      {
        IP.row = dIVec[index].gid;
        IP.col = dIVec[index].Vcol[j];
        indexPairList.push_back(IP);
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2)
        {
          cout << "Vcol["<<j<<"] = " << dIVec[index].Vcol[j] <<endl;
        }
#endif
      }
    }

    for (j=0;j<nn;++j)
    {
      // electron (row, col) pairs.
      if (dIVec[index].Ncol[j] != -1)
      {
        IP.row = dIVec[index].gid;
        IP.col = dIVec[index].Ncol[j];
        indexPairList.push_back(IP);
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2)
        {
          cout << "Ncol["<<j<<"] = " << dIVec[index].Ncol[j] <<endl;
        }
#endif
      }
    }

    for (j=0;j<nn;++j)
    {
      // hole (row, col) pairs.
      if (dIVec[index].Pcol[j] != -1)
      {
        IP.row = dIVec[index].gid;
        IP.col = dIVec[index].Pcol[j];
        indexPairList.push_back(IP);
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2)
        {
          cout << "Pcol["<<j<<"] = " << dIVec[index].Pcol[j] <<endl;
        }
#endif
      }
    } // end of j loop.
  } // end of DI loop

  // interior points next:
  for (i=0;i<numMeshPoints;++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 2)
    {
      cout << "Mesh Index: i = " << i << "  label = ";
      cout << labelNameVector[i] <<endl;
    }
#endif
    nn = Vcolarray[i].size();

    for (j=0;j<nn;++j)
    {
      if (Vrowarray[i] != -1 && Vcolarray[i][j] != -1)
      {
        IP.row = Vrowarray[i];
        IP.col = Vcolarray[i][j];
        indexPairList.push_back(IP);
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2)
        {
          cout << "  V row,col = " << Vrowarray[i] << ", ";
          cout << Vcolarray[i][j] << "\n";
        }
#endif
      }
    }

    nn = Ncolarray[i].size();
    for (j=0;j<nn;++j)
    {
      if (Nrowarray[i] != -1 && Ncolarray[i][j] != -1)
      {
        IP.row = Nrowarray[i];
        IP.col = Ncolarray[i][j];
        indexPairList.push_back(IP);
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2)
        {
          cout << "  N row,col = " << Nrowarray[i] << ", ";
          cout << Ncolarray[i][j] << "\n";
        }
#endif
      }
    }

    nn = Pcolarray[i].size();
    for (j=0;j<nn;++j)
    {
      if (Prowarray[i] != -1 && Pcolarray[i][j] != -1)
      {
        IP.row = Prowarray[i];
        IP.col = Pcolarray[i][j];
        indexPairList.push_back(IP);
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 2)
        {
          cout << "  P row,col = " << Prowarray[i] << ", ";
          cout << Pcolarray[i][j] << "\n";
        }
#endif
      }
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 2)
  {
    cout << endl;
  }
#endif

  // now that the  index pair list has been set up,
  // print it out to the screen.
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 2)
  {
    list<index_pair>::iterator ip_iter;
    for (i=0,
         ip_iter  = indexPairList.begin();
         ip_iter != indexPairList.end();
         ++ip_iter,++i)
    {
       cout << "i="<<i<<":  (";
       cout.width(6);
       cout << ip_iter->row;
       cout << ", ";
       cout.width(6);
       cout <<  ip_iter->col;
       cout << ")\n";
    }
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateGIDs
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/01
//-----------------------------------------------------------------------------
void Instance::registerStateGIDs(
  const list<index_pair> & staGIDListRef)
{

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "\n";
    cout << dashedline << "\n";
    cout << "  In Instance::registerStateGIDs\n\n";
    cout << "  name             = " << getName() << "\n";
  }
#endif

  string msg;

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staGIDListRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateGIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  Number of State GIDs: " << numSta << "\n";
  }
#endif

  // Copy over the global ID lists:
  staGIDList.assign(staGIDListRef.begin(), staGIDListRef.end());

  list<index_pair>::iterator it1 = staGIDList.begin();
  list<index_pair>::iterator last1 = staGIDList.end ();

  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  int i=0;
  for (; (iterDI!=lastDI && it1!=last1);++iterDI,++it1,++i)
  {
    iterDI->stateC       = (*it1).row;
    iterDI->stateC_owned = (it1->col == 1);
  }

  for (i=0;i<numMeshPoints;++i,++it1)
  {
     stateDispl[i] = (*it1).row;
     stateDispl_owned[i] = static_cast<int>((it1->col == 1));
  }


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1)
  {
    cout << "  State indices:" << "\n";
    cout << "\n";
    for (iterDI=firstDI; iterDI!=lastDI;++iterDI)
    {
      cout << "  ";
      cout.width(12);
      cout.setf(ios::right);
      cout << iterDI->eName;
      cout.setf(ios::left);
      cout << "  stateC = " << iterDI->stateC;
      cout << "  stateC_owned = " << iterDI->stateC_owned << "\n";
    }

    cout << dashedline << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : DiodeDPEInstance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                                        const vector<int> & extLIDVecRef)

{
  string msg;
  string tmpstr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << dashedline << endl;
    cout << "In the Intance::registerLIDs function.  ";
    cout << "  name = "<< getName() <<endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "numInt = " << numInt <<"\n";
    cout << "numExt = " << numExt <<"\n";
    cout << "numMeshPoints = " << numMeshPoints << "\n";
  }
#endif

  // number of internal variables equals the number of
  // mesh points *3.
  if (numInt != numIntVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numInt != 3*numMeshPoints.  ";
    msg += "Check the metadata in the parser\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  // copy over the global ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Set up the rows first, then do the colarrays.

  // First do the external variables:
  // These will all be voltage nodes connected to the devices.
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Setting up the indices for the external circuit nodes:" << endl;
    cout << "External node list:" << endl;
  }
#endif

  int isizeDI = dIVec.size();

  int index;
  for(index=0; index < isizeDI; ++index)
  {
     dIVec[index].lid = extLIDVec[index];

#ifdef Xyce_DEBUG_DEVICE
     if (getDeviceOptions().debugLevel > 1)
     {
       cout << "   name = "<<dIVec[index].eName<<" lid = ";
       cout << dIVec[index].lid;
       cout << endl;
     }
#endif
  }


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Doing internal vars, row arrays:"<<endl;
  }
#endif

  // Do the internal variables.  There should be a lot of these.
  int i=0;
  index = 0;

  // The interior points will be blocked (V,N,P) together.
  while (i < numMeshPoints)
  {
#ifdef Xyce_NEW_BC
    if (boundarySten[i]) { ++i; continue; }

    if (!(boundaryStenV[i]))
    {
      li_Vrowarray[i] = intLIDVec[index];   // should be Vrowarray[i]
      ++index;
    }

    if (!(boundaryStenN[i]))
    {
      li_Nrowarray[i] = intLIDVec[index];
      ++index;
    }

    if (!(boundaryStenP[i]))
    {
      li_Prowarray[i] = intLIDVec[index];
      ++index;
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1)
    {
      cout << "doing lid row arrays for mesh point " << i << endl;
      if (!(boundaryStenV[i]))
	cout << "  li_Vrow = " << li_Vrowarray[i] << endl;

      if (!(boundaryStenN[i]))
	cout << "  li_Nrow = " << li_Nrowarray[i] << endl;

      if (!(boundaryStenP[i]))
	cout << "  li_Prow = " << li_Prowarray[i] << endl;
    }
#endif

#else
    li_Vrowarray[i] = intLIDVec[index];   // should be Vrowarray[i]
    ++index;

    li_Nrowarray[i] = intLIDVec[index];
    ++index;

    li_Prowarray[i] = intLIDVec[index];
    ++index;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 1)
    {
      cout << "doing lid row arrays for mesh point " << i << endl;
      cout << "  li_Vrow = " << li_Vrowarray[i] << endl;
      cout << "  li_Nrow = " << li_Nrowarray[i] << endl;
      cout << "  li_Prow = " << li_Prowarray[i] << endl;
    }
#endif
#endif // Xyce_NEW_BC
    ++i;
  }

#ifdef Xyce_DEBUG_DEVICE
  cout << dashedline << endl;
#endif


}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef)
{

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "\n";
    cout << dashedline << "\n";
    cout << "  In Instance::registerStateLIDs\n\n";
    cout << "  name             = " << getName() << "\n";
  }
#endif

  string msg;

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  Number of State LIDs: " << numSta << "\n";
  }
#endif

  // Copy over the local ID lists:
  staLIDVec = staLIDVecRef;

  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  int i=0,j=0;
  for (; iterDI!=lastDI;++iterDI,++i)
  {
    iterDI->li_stateC = staLIDVec[i];
  }

  for (j=0;j<numMeshPoints;++j,++i)
  {
     li_stateDispl[j] = staLIDVec[i];
  }


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1)
  {
    cout << "  State indices:" << "\n";
    cout << "\n";
    for (iterDI=firstDI; iterDI!=lastDI;++iterDI)
    {
      cout << "  ";
      cout.width(12);
      cout.setf(ios::right);
      cout << iterDI->eName;
      cout.setf(ios::left);
      cout << "  li_stateC = " << iterDI->li_stateC;
      cout << endl;
    }

    cout << "  Displacement state indices:\n";
    for (j=0;j<numMeshPoints;++j,++i)
    {
       cout << "  edge: " << j << "  li_stateDispl = " << li_stateDispl[j];
       cout << endl;
    }

    cout << dashedline << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233
// Creation Date : 02/23/03
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
//
// Purpose       : This function sets up the "local-ID" (LID) jacobian
//                 stamp.  This is neccessary for the DMA=on capability,
//                 which is the default.  This represents the stuff that
//                 comes in from topology.
//
// Special Notes : This needs to be consistent with the function,
//                 Instance::setupJacStamp.
//
//
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233
// Creation Date : 02/23/03
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs
    ( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs ( jacLIDVec );

  int i;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "In Instance::registerJacLIDs" << endl;
  }
#endif

  /////////////////////////////////////////////////////////////////////////
  // Do the external variables first:
  // Loop over the boundary condition "device interface" labels.
  vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  vector<DeviceInterfaceNode>::iterator iterDI;
  int DIsize = dIVec.size ();
  firstDI = dIVec.begin();
  lastDI  = dIVec.end  ();

  int index;
  for(index=0,iterDI=firstDI;iterDI!=lastDI;++iterDI,++index)
  {
    int jacRowSize = jacLIDVec[index].size();
    dIVec[index].dIdXoffset.resize(jacRowSize-dIVec[index].numCrossTerms-1);
    dIVec[index].lidOffset = jacLIDVec[index][0];

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << "index = " << index;
      cout << "  jacRowSize = " << jacRowSize;
      cout << "  name = " << dIVec[index].eName << endl;
      cout << " lidOffset = ";
      cout << dIVec[index].lidOffset << endl;
    }
#endif

    int nCT = dIVec[index].numCrossTerms;

    dIVec[index].crossOffsets.resize(nCT);

    for(int itmp=0;itmp<nCT;++itmp)
    {
      dIVec[index].crossOffsets[itmp] = jacLIDVec[index][itmp+1];
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        cout << "  crossOffsets["<<itmp<<"] = ";
        cout << dIVec[index].crossOffsets[itmp] << endl;
      }
#endif
    }

    for (int ioff=nCT+1;ioff<jacRowSize;++ioff)
    {
      int tmpIndex = ioff - (nCT+1);
      dIVec[index].dIdXoffset[tmpIndex] = jacLIDVec[index][ioff];

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        cout << " dIdXoffset["<<tmpIndex<<"] = ";
        cout << dIVec[index].dIdXoffset[tmpIndex] << endl;
      }
#endif
    }
  } // index (dIVec) loop.

  /////////////////////////////////////////////////////////////////////////
  // Now do the internal variables.  (V,N,P on the mesh)
  li_VoffsetArray.resize(numMeshPoints);
  li_NoffsetArray.resize(numMeshPoints);
  li_PoffsetArray.resize(numMeshPoints);
  for(i=0;i<numMeshPoints;++i)
  {
#ifdef Xyce_NEW_BC
    if (boundarySten[i]) continue;
#endif // Xyce_NEW_BC

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << " mesh point i = " << i << endl;
    }
#endif

    mNode * nodePtr = meshContainerPtr->getNode(i);
    vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin ();
    vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end   ();

    // get the temporary LID row indices:
    int Vrow = MESHtoLID_V[i];
    int Nrow = MESHtoLID_N[i];
    int Prow = MESHtoLID_P[i];

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      cout << "   Vrow = " << Vrow << endl;
      cout << "   Nrow = " << Nrow << endl;
      cout << "   Prow = " << Prow << endl;
    }
#endif

    int ioff;
    // voltage row:
    if (Vrow != -1)
    {
      int VrowSize = jacLIDVec[Vrow].size();
      li_VoffsetArray[i].resize(VrowSize);
      for (ioff=0;ioff<VrowSize;++ioff)
      {
        li_VoffsetArray[i][ioff] = jacLIDVec[Vrow][ioff];
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 1)
        {
          cout << "   li_Voffset["<<i<<"]["<<ioff<<"] = ";
          cout << li_VoffsetArray[i][ioff] << endl;
        }
#endif
      }
    }

    // electron row:
    if (Nrow != -1)
    {
      int NrowSize = jacLIDVec[Nrow].size();
      li_NoffsetArray[i].resize(NrowSize);
      for (ioff=0;ioff<NrowSize;++ioff)
      {
        li_NoffsetArray[i][ioff] = jacLIDVec[Nrow][ioff];
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 1)
        {
          cout << "   li_Noffset["<<i<<"]["<<ioff<<"] = ";
          cout << li_NoffsetArray[i][ioff] << endl;
        }
#endif
      }
    }

    // hole row:
    if (Prow != -1)
    {
      int ProwSize = jacLIDVec[Prow].size();
      li_PoffsetArray[i].resize(ProwSize);
      for (ioff=0;ioff<ProwSize;++ioff)
      {
        li_PoffsetArray[i][ioff] = jacLIDVec[Prow][ioff];
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 1)
        {
          cout << "   li_Poffset["<<i<<"]["<<ioff<<"] = ";
          cout << li_PoffsetArray[i][ioff] << endl;
        }
#endif
      }
    }

  } // mesh point loop
}

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce
