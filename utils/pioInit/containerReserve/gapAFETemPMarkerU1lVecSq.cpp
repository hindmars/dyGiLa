// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

// #include "plumbing/hila.h"
// //#include "plumbing/fft.h"
// #include "plumbing/memalloc.h"

#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::containerReserve_gapAFETemPMarkerU1lVecSq(glsol &sol) {

    if (sol.config.pario_compute_gapA == 1) { gapAContainer.reserve(latticeVolumeWithGhost); }
    
    if (sol.config.pario_compute_feDensity == 1) { feDensityContainer.reserve(latticeVolumeWithGhost); }    
    if (sol.config.pario_Temperature_pStream == 1) { Temperature.reserve(latticeVolumeWithGhost); }
    
    if (sol.config.pario_compute_phaseMarker == 1) { phaseMarker.reserve(latticeVolumeWithGhost); }
    if (sol.config.pario_compute_DB == 1) { DBContainer.reserve(latticeVolumeWithGhost); }

    if (sol.config.pario_compute_U1Phase == 1) { U1_3phiContainer.reserve(latticeVolumeWithGhost); }

    if (sol.config.pario_compute_lVector == 1) { lsqContainer.reserve(latticeVolumeWithGhost); }    
    
} // containerReserve_gapAFETemPMarkerU1() end here

