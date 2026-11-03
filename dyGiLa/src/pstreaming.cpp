#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/globals.h" 

#include "glsol.hpp"
#include "matep_namespace_utils.hpp"
#include "dyGiLa.hpp"

#if defined USE_PARIO 
#include "pario.hpp"
#endif

namespace dyGiLa {
  
void pStreaming(glsol &gl, parIO &paraio, unsigned int &stat_counter, const unsigned int &steps)
{
  const unsigned int modPSSR = (stat_counter / steps) % gl.config.PSSRatio;
  
  // phase marking only for streaming when gamma isn't T-dependent
  // TDependnetgamma == true is hanlded in main(), so you have to
  // have it be handled at here for false case.
  if (gl.config.TDependnetgamma == false)
    {
     // AwT refresh for TDependnetgamma == false case
     // you have to do this because true case is handled elsewhere. 
     if (gl.config.useGaussianLP_filter == true)
       { gl.GaussianLPfilter_matrix(/*gl.AwT*/); }
     // do phase-Marking only for pio when gamma is constant
     gl.phaseMarking();
     if (modPSSR == 0)
       {
        hila::out0 << "gl.t is " << gl.t
		   << ", phaseMarking() call is done for const gamma. "
		   << std::endl;

       }
    }
  
  // reducntions output
  gl.write_energies();
  gl.phaseCounting();
  if (modPSSR == 0)
    {
     hila::out0 << "write_energies(), phaseCounting() call is done "
	        << std::endl;
    }

  // parallel streaming
#if defined USE_PARIO
  if (
      ((gl.config.hdf5_A_matrix_output == 1)
      || (gl.config.hdf5_AMatrix_exaslice_output == 1)       
      || (gl.config.hdf5_pMarker_output == 1)
      || (gl.config.hdf5_mass_current_output == 1)
      || (gl.config.hdf5_spin_current_output == 1)
      || (gl.config.hdf5_spinCurr_exaslice_output == 1)              
      || (gl.config.hdf5_lVector_output == 1)
      || (gl.config.hdf5_GradPhiVector_output == 1)
      || (gl.config.hdf5_GradPhiVector_exaslice_output == 1)
      || (gl.config.hdf5_GradPhiVector_clip_output == 1)
      || (gl.config.hdf5_GradPhiVector_boxclip_output == 1))                  
      && (modPSSR == 0)
     )
    {
      if (gl.config.hdf5_cheerypick == true)
	{
	  for (unsigned int &cp_marker : gl.config.hdf5_cheerypick_array)
	    { if (stat_counter == cp_marker) paraio.pstream(gl, stat_counter); }
	}
      else
	{
         if (gl.t >= gl.config.hdf5Ststart && gl.t <= gl.config.hdf5Stend)
	   paraio.pstream(gl, stat_counter);
	}
    }    
  else if (
           // insitu visualization block, no parallel hd5 stream
           (gl.config.hdf5_A_matrix_output != 1)
	   && (gl.config.hdf5_AMatrix_exaslice_output != 1)
           && (gl.config.hdf5_pMarker_output != 1)
           && (gl.config.hdf5_mass_current_output != 1)
           && (gl.config.hdf5_spin_current_output != 1)
           && (gl.config.hdf5_spinCurr_exaslice_output != 1)	   
           && (gl.config.hdf5_lVector_output != 1)
           && (gl.config.hdf5_GradPhiVector_output != 1)
           && (gl.config.hdf5_GradPhiVector_exaslice_output != 1)
	   && (gl.config.hdf5_GradPhiVector_clip_output != 1)
	   && (gl.config.hdf5_GradPhiVector_boxclip_output != 1)	   
	   && (
               (gl.config.do_gapA_clip == 1)
   	       || (gl.config.do_gapA_slice == 1)
	       || (gl.config.do_fed_clip == 1)
	       || (gl.config.do_Temperature_clip == 1)
	       || (gl.config.do_Temperature_slice == 1)
	       || (gl.config.do_Temperature_isosurface == 1)			   
	       || (gl.config.do_gapA_isosurface == 1)
	       || (gl.config.do_phaseMarker_slice1 == 1)
	       || (gl.config.do_phaseMarker_slice2 == 1)	       
	       || (gl.config.do_DB_slice == 1)
	       || (gl.config.do_phaseMarker_isosurface == 1)
	       || (gl.config.do_phaseMarker_fieldclip == 1)
	       || (gl.config.do_phaseMarker_fieldclip_Bphase == 1)
	       || (gl.config.do_phaseMarker_fieldclip_Aphase == 1)
	       || (gl.config.do_U13phi_slice == 1)
	       || (gl.config.do_l_sq_slice == 1)
              )
           && (modPSSR == 0)	  
          )
	paraio.pstream(gl, stat_counter);		

  if (modPSSR == 0) { hila::out0 << "paraio.pstream() call is done " << std::endl; }
#endif	            
} // pstreaming func ends here

} // dyGiLa namespace ends here
