#define USE_MPI
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

/*
 * Register D_B = sqrt(sum_a (R_a - R_a^B)^2) with the Ascent mesh as a
 * vertex-associated scalar field named "DB".  The values are filled by
 * pstream() from glsol::DB, which phaseMarking() computes.
 */
void parIO::describeMesh_DB(glsol &sol) {

    if (sol.config.pario_compute_DB == 1)
      {
       mesh["fields/DB/association"] = "vertex";
       mesh["fields/DB/topology"] = "topo";
       mesh["fields/DB/values"].set_external(DBContainer.data(), latticeVolumeWithGhost);
      }

} // describeMesh_DB() end here
