#include <string>
#include <iostream>
#include <fstream>
#include <list>

#ifdef Xyce_PARALLEL_MPI
  #include <mpi.h>
#endif

#include <N_TOP_CktNode.h>
#include <N_TOP_CktNodeCreator.h>
#include <N_TOP_CktGraphSupport.h>
#include <N_TOP_CktGraphCreator.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_Topology.h>
//#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>
//#include <IO_psn.h>
#include <N_IO_DistribMgr.h>
#include <N_DEV_DeviceMgr.h>
#include <N_UTL_Xyce.h>

int main (int argc, char* argv[]) {

#ifdef Xyce_PARALLEL_MPI
    MPI_Init( &argc, &argv );
#endif

//    N_PDS_Manager PDSMgr;

    N_PDS_Comm PDSComm;

#ifdef Xyce_PARALLEL_MPI
    N_ERH_ErrorMgr::registerComm( &PDSComm );
#endif

    N_DEV_DeviceMgr* DevMgrPtr = N_DEV_DeviceMgr::factory();

    N_TOP_Topology Topo( DevMgrPtr );

//    ParseSpiceNetlist psn ( argv[1], Topo, *DevMgrPtr, true );

//    int isuccess = psn.PSN_Read();

    N_IO_DistribMgr* DistMgrPtr = N_IO_DistribMgr::factory();
    DistMgrPtr->registerTopology( &Topo );
    DistMgrPtr->registerDeviceMgr( DevMgrPtr );

#ifdef Xyce_PARALLEL_MPI
    DistMgrPtr->registerParallelServices( &PDSComm );
#endif

    DistMgrPtr->generateParser( string( argv[1] ) );
    DistMgrPtr->run();

    cout << Topo << endl;

    Topo.OutputBFSGraphLists();

    Topo.OutputDFSGraphLists();

    DevMgrPtr->printOutLists();

#ifdef Xyce_PARALLEL_MPI
    MPI_Finalize();
#endif

    return 0;
}
