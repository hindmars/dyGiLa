// #define USE_ADGRZ
#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"
//#include "plumbing/memalloc.h"

#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::init(glsol &sol) {

    latticeVolumeWithGhost =
        (lattice->mynode.size[0] + 2) * (lattice->mynode.size[1] + 2) * (lattice->mynode.size[2] + 2);  
    latticeVolume =
        (lattice->mynode.size[0]) * (lattice->mynode.size[1]) * (lattice->mynode.size[2]);

    /*********************************/
    /*    container reserve calls    */
    /*********************************/       
    containerReserve_gapAFETemPMarkerU1lVecSq(sol);
    if (
	sol.config.hdf5_mass_current_output == 1
	// || (sol.config.hdf5_massCurr_exaclice_output == 1)
       ){ containerReserve_massCurrent(); }
    if (
	(sol.config.hdf5_spin_current_output == 1)
	|| (sol.config.hdf5_spinCurr_exaslice_output == 1)
       ) { containerReserve_spinCurrent(); }
    if (
	(sol.config.hdf5_A_matrix_output == 1)
	|| (sol.config.hdf5_AMatrix_exaslice_output == 1)
       ) { containerReserve_Amatrix(); }
    if (sol.config.hdf5_lVector_output == 1) { containerReserve_lVec(); }
    if (
	(sol.config.hdf5_GradPhiVector_output == 1)
	|| (sol.config.hdf5_GradPhiVector_exaslice_output == 1)
	|| (sol.config.hdf5_GradPhiVector_clip_output == 1)
	|| (sol.config.hdf5_GradPhiVector_boxclip_output == 1)		
       )
      { containerReserve_GradientPhiVec(); }        

    // ghost handling
    ghostMask(sol);
    
    /*********************************/
    /*    all describeMesh calls     */
    /*********************************/    
    describeMesh(sol);
    describeMesh_gapA_FEDensity(sol);
    if (sol.config.pario_Temperature_pStream == 1) { describeMesh_Temperature(); }
    describeMesh_phaseMarker(sol);
    describeMesh_DB(sol);
    if (sol.config.pario_compute_U1Phase == 1) { describeMesh_U13phi(); }
    if (sol.config.pario_compute_lVector == 1) { describeMesh_lVecSq(); }

    if (
	( sol.config.hdf5_cheerypick == true )
	|| ((!!sol.config.hdf5Ststart == true) && (!!sol.config.hdf5Stend == true))
       )
      {
       if (
	   sol.config.hdf5_mass_current_output == 1
	   // || (sol.config.hdf5_massCurr_exaclice_output == 1)
	  ) {describeMesh_massCurrent();}
       if (
	   (sol.config.hdf5_spin_current_output == 1)
           || (sol.config.hdf5_spinCurr_exaslice_output == 1)
	  ) {describeMesh_spinCurrent();}
       if (
	   (sol.config.hdf5_A_matrix_output == 1)
	   || (sol.config.hdf5_AMatrix_exaslice_output == 1)	   
	  ) {describeMesh_AMatrix();}
       if (sol.config.hdf5_lVector_output == 1) {describeMesh_lVec();}
       if (
	   (sol.config.hdf5_GradPhiVector_output == 1)
	   || (sol.config.hdf5_GradPhiVector_exaslice_output == 1)
	   || (sol.config.hdf5_GradPhiVector_clip_output == 1)
           || (sol.config.hdf5_GradPhiVector_boxclip_output == 1)	   
	  ) { describeMesh_GradientPhiVec(); }	               
      }

    describeMesh_addGhost_verify();
    /*********************************/
    /*  describeMesh calls end here  */
    /*********************************/    
    
    pio_options["mpi_comm"] = MPI_Comm_c2f(lattice->mpi_comm_lat);
    pio_options["runtime/type"] = "ascent";
#if defined CUDA
    pio_options["runtime/vtkm/backend"] = "cuda";
    pio_options["cuda/init"] = "false";
#endif    
    pio_options["timings"] = "false";    
    pio.open(pio_options);
    pio.publish(mesh);

    /*********************************/
    /*   all defineActions calls     */
    /*********************************/        
    defineActions_insitu(sol);
    if (
	( sol.config.hdf5_cheerypick == true )
	|| ((!!sol.config.hdf5Ststart == true) && (!!sol.config.hdf5Stend == true))	
       )
      {
       if (sol.config.hdf5_lVector_output == 1) {defineActions_lVec();}
       if (sol.config.hdf5_GradPhiVector_output == 1) {defineActions_GradientPhiVec();}
       if (sol.config.hdf5_GradPhiVector_exaslice_output == 1) {defineActions_GradientPhiVec_exaslice(sol);}
       if (sol.config.hdf5_GradPhiVector_clip_output == 1) {defineActions_GradientPhiVec_clip(sol);}
       if (sol.config.hdf5_GradPhiVector_boxclip_output == 1) {defineActions_GradientPhiVec_boxclip(sol);}
       if (sol.config.hdf5_mass_current_output == 1) {defineActions_massCurrent(sol);}
       // if (sol.config.hdf5_massCurr_exaclice_output == 1) {defineActions_massCurrent_exaslice(sol);}       
       if (sol.config.hdf5_spin_current_output == 1) {defineActions_spinCurrent(sol);}
       if (sol.config.hdf5_spinCurr_exaslice_output == 1) {defineActions_spinCurrent_exaslice(sol);}       
       if (sol.config.hdf5_A_matrix_output == 1) {defineActions_AMatrix(sol);}
       if (sol.config.hdf5_AMatrix_exaslice_output == 1) {defineActions_AMatrix_exaslice(sol);}       
       if (sol.config.hdf5_pMarker_output == 1) {defineActions_phaseMarker(sol);}
      }
    
    defineActions_printTree();
    
    /*********************************/
    /*   all defineActions calls     */
    /*********************************/
    
    hila::out0 << "------------------------------------------------------------" << std::endl;
    
} // init() end here

