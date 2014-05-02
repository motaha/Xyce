//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2003, 2013, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : Xyce.C
//
// Purpose        : front end for standalone Xyce executable
//
// Special Notes  :
//
// Creator        : Eric Rankin
//
// Creation Date  : 01/28/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2013/09/30 20:25:21 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------


#include <N_CIR_Xygra.h>
#include <N_UTL_fwd.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Device {
void registerDevices();
}}


// Function to be called if memory runs out:
void _new_handler (void)
{
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, "OUT OF MEMORY (error in 'new')");
}

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       : front end for standalone Xyce executable
// Special Notes :
// Scope         :
// Creator       : Eric Rankin
// Creation Date : 01/28/2004
//-----------------------------------------------------------------------------
int main( int iargs, char *cargs[] )
{
  // Set divide by zero, and invalid operation handling on linux
  // Set out of memory detection on all systems
  std::set_new_handler (&_new_handler);

  Xyce::Device::registerDevices();

  N_CIR_Xygra * XycePtr = new N_CIR_Xygra();

  bool bsuccess = XycePtr->initialize(iargs, cargs);
  std::vector <std::string> deviceNames;
  std::vector <double> vN;

  if (bsuccess)
    bsuccess = XycePtr->getDeviceNames("Y%XYGRA%DUMMY",deviceNames);

  Xyce::lout() << "Starting Xyce as library test" << std::endl;

  if (bsuccess)
  {
    std::vector<std::vector<int> >coilWindings;
    std::vector<std::vector<std::string> >coilNames;
    coilWindings.resize(deviceNames.size());
    coilNames.resize(deviceNames.size());
    for (int i=0; i < deviceNames.size(); ++i)
    {
      XycePtr->xygraGetCoilWindings(deviceNames[i],coilWindings[i]);
      XycePtr->xygraGetCoilNames(deviceNames[i],coilNames[i]);
      Xyce::dout() << " Xygra device " << deviceNames[i] << " has "
           << coilWindings[i].size() << " coils " << std::endl;
      for (int j=0; j<coilWindings[i].size(); j++)
      {
        Xyce::dout() << "    coil " << j << " is named " << coilNames[i][j] << " and has " << coilWindings[i][j]
             << "windings" << std::endl;
      }
    }

//        bsuccess = XycePtr->runSimulation();
    double completedTime, timeStep;
    completedTime = 0.0;
    timeStep = 1e-3;
    bool opComplete = false;

    while (!(XycePtr->simulationComplete()) && bsuccess)
    {

      // We'll set each winding to a current source to simulate a rough
      // sinusoidal voltage with one period over the entire run
      for (int i = 0; i<deviceNames.size(); ++i)
      {
        int numWindings = XycePtr->xygraGetNumWindings(deviceNames[i]);

        if (completedTime == 0)
        {
          // set the t=0 version first
          std::vector<double> sV(numWindings,0.0);
          XycePtr->xygraSetSources(deviceNames[i],sV,completedTime);
        }

        double currentValue =  .001*sin(2*3.14159265358979*10.0*(completedTime+timeStep));
        std::vector<double> sV(numWindings,currentValue);
        XycePtr->xygraSetSources(deviceNames[i],sV,completedTime+timeStep);
      }

      bsuccess = XycePtr->simulateUntil(completedTime+timeStep,completedTime);
    }

    if (bsuccess)
      bsuccess = XycePtr->finalize();
  }

  delete XycePtr;

  (bsuccess) ? exit(0) : exit(-1);
}

