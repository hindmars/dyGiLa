#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>
#include <sys/stat.h>

#include "plumbing/hila.h"
#include "plumbing/globals.h" 

#include "glsol.hpp"
#include "matep_namespace_utils.hpp"
#include "dyGiLa.hpp"

#if defined USE_PARIO 
#include "pario.hpp"
#endif

namespace {

// Write each configured slice of A as raw binary to
// hila-A-slice/hila-A-slice<k>_t-<stat_counter>.bin, and append "stat_counter t" to
// hila-A-slice/hila-A-slice-times.txt. Collective: every rank must call it.
void writeHilaASlices(glsol &gl, const unsigned int stat_counter)
{
  const std::vector<int> &c = gl.config.hila_A_slice_coords;
  static bool first = true;

  if (first && hila::myrank() == 0)
    {
      mkdir("hila-A-slice", 0755);   // fails harmlessly if it exists
      std::ofstream info("hila-A-slice/hila-A-slice-info.txt");
      info << "# Raw binary slices of the order parameter A, written by Field::write_slice_binary.\n"
	   << "# Each site is a 3x3 complex matrix of " << sizeof(real_t) << "-byte floats:\n"
	   << "# A[0][0].re, A[0][0].im, A[0][1].re, ..., A[2][2].im (row-major, 18 numbers).\n"
	   << "# Sites are in lattice order with x fastest, then y, then z.\n"
	   << "# A coordinate of -1 means the whole range in that direction.\n"
	   << "lattice " << gl.config.lx << " " << gl.config.ly << " " << gl.config.lz << "\n"
	   << "dx " << gl.config.dx << "\n";
      for (size_t k = 0; k < c.size() / 3; k++)
	info << "slice " << k << " " << c[3*k] << " " << c[3*k+1] << " " << c[3*k+2] << "\n";
    }
  first = false;

  for (size_t k = 0; k < c.size() / 3; k++)
    {
      const CoordinateVector slice = {c[3*k], c[3*k+1], c[3*k+2]};
      std::ostringstream fname;
      fname << "hila-A-slice/hila-A-slice" << k << "_t-"
	    << std::setfill('0') << std::setw(9) << stat_counter << ".bin";
      gl.A.write_slice_binary(fname.str(), slice);
    }

  if (hila::myrank() == 0)
    {
      std::ofstream times("hila-A-slice/hila-A-slice-times.txt", std::ios::app);
      times << stat_counter << " " << gl.t << "\n";
    }
}

} // anonymous namespace

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

  // raw A slices, at the same outputs as the in-situ images
  if (gl.config.write_hila_A_slice == 1 && modPSSR == 0)
    {
     writeHilaASlices(gl, stat_counter);
     hila::out0 << "writeHilaASlices() call is done " << std::endl;
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
