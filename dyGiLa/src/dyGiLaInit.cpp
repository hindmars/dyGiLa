//#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>
#include <cmath>

#include "plumbing/hila.h"
#include "plumbing/globals.h" 

#include "glsol.hpp"
#include "matep_namespace_utils.hpp"
#include "dyGiLa.hpp"

// #if defined USE_PARIO 
// #include "pario.hpp"
// #endif

namespace dyGiLa {
  
std::tuple<const std::vector<std::string>, const CoordinateVector *const, const unsigned int> dyGiLaInit(glsol &gl, int &argc, char **&argv) {
  // lattice initialization
  hila::initialize(argc, argv);
  std::vector<std::string> name_files = gl.configure("sim_params.txt", argc, argv);
  CoordinateVector box_dimensions = {gl.config.lx, gl.config.ly, gl.config.lz};
  lattice.setup(box_dimensions);
  hila::seed_random(gl.config.seed);

  // host & device memory initialization, gpuMemcpHostToDevice under hood
  matep::init_wrapper_mp();

  // log subsection 
  hila::out0 << "------------------------------------------------------------" << "\n"
             << "--         dyGiLa Simulation Suits Initial State          --" << "\n"
             << "------------------------------------------------------------" << "\n"
             << std::endl;
  
  // initialize Temperature field
  gl.initializeT();

  // initilize static H-field
  gl.initializeH();

  // initialize OP field
  gl.initialize();

  // coordinate origin
  std::initializer_list<int> coordsList {0,0,0};
  // static: a pointer to it is returned, so it must outlive this function
  static const CoordinateVector originpoints(coordsList);

  // number of steps between reduction streaming
  const unsigned int steps = std::round((gl.config.tEnd - gl.config.tStats)
					/(gl.config.dt * gl.config.nOutputs));
  hila::out0 << "steps is " << steps << " in unit of dt. "
	     << "PSSRatio is " << gl.config.PSSRatio  << "." << std::endl;

  // return initialization list
  return std::make_tuple(name_files, &originpoints, steps);
  
} /* dyGiLaInit() end here */

} // orch namespace ends here

