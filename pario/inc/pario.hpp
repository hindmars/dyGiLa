#ifndef PARIO_HPP
#define PARIO_HPP

#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "matep.hpp"
#include "glsol.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


// Definition of the field that we will use
using real_t = float;                          // or double ?
// using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time

class parIO{

public:
  parIO() = default;                     // default constructor
  
  matep::Matep matep;
  
  /*----------------------------------------*/
  /*  parallel-IO memember functions        */ 
  /*----------------------------------------*/

  // called in main() before t-while-loop started
  void xdmf(glsol &);
  void xml_Amatrix(glsol &);
  void xml_pMarker(glsol &);
  void xml_massCurrent(glsol &);
  void xml_spinCurrent(glsol &);  

  // called in main()  
  void init(glsol &);
  void pstream(glsol &, unsigned int &);
  void shutdown();

  // xdmf file fstream
  std::fstream xdmf_out;
  std::fstream xml_out;

private:
  
  /* functions called pario engine_init() */
  /*  mesh descriptions and definations   */
  void describeMesh(glsol &);
  void describeMesh_addGhost_verify();
  
  void describeMesh_gapA_FEDensity(glsol &);
  void describeMesh_Temperature();
  void describeMesh_phaseMarker(glsol &);
  void describeMesh_DB(glsol &);
  void describeMesh_U13phi();
  void describeMesh_lVecSq();  
  void describeMesh_lVec();
  void describeMesh_GradientPhiVec();    
  void describeMesh_massCurrent();
  void describeMesh_spinCurrent();  
  void describeMesh_AMatrix();  

  /* actions definations */
  void defineActions_insitu(glsol &);

  void defineActions_lVec();
  void defineActions_GradientPhiVec();
  void defineActions_GradientPhiVec_exaslice(glsol &);
  void defineActions_GradientPhiVec_clip(glsol &);
  void defineActions_GradientPhiVec_boxclip(glsol &);    
  void defineActions_massCurrent(glsol &);
  //void defineActions_massCurrent_exaslice(glsol &);  
  void defineActions_spinCurrent(glsol &);
  void defineActions_spinCurrent_exaslice(glsol &);  
  void defineActions_AMatrix(glsol &);
  void defineActions_AMatrix_exaslice(glsol &);
  void defineActions_phaseMarker(glsol &);
  
  void defineActions_printTree();

  /* utils */
  void containerReserve_gapAFETemPMarkerU1lVecSq(glsol &);
  void containerReserve_lVec();
  void containerReserve_GradientPhiVec();  
  void containerReserve_Amatrix();
  void containerReserve_massCurrent();
  void containerReserve_spinCurrent();
  void ghostMask(glsol &);

  void U1PhaseStreaming(glsol &);
  void lVectorStreaming(glsol &);
  void GradientPhiVectorStreaming(glsol &);    

  /*----- fields declearations -----*/
  
  Field<real_t> gapA;
  Field<real_t> feDensity;
  // Field<real_t> trA_re, trA_im;
  Field<real_t> u11, u12, u13, u21, u22, u23, u31, u32, u33;
  Field<real_t> v11, v12, v13, v21, v22, v23, v31, v32, v33;
  // Field<real_t> eigAv1, eigAv2, eigAv3;

  //Field<Vector<3,real_t>> lVec;
  
  Field<real_t> jm1, jm2, jm3;
  Field<real_t> U1_3phi;
  Field<real_t> lsq, l_1, l_2, l_3;
  Field<real_t> GPhi_1, GPhi_2, GPhi_3;    
  Field<real_t> js11, js21, js31,
                js12, js22, js32,
                js13, js23, js33;
   
  std::vector<real_t> gapAContainer;
  std::vector<real_t> feDensityContainer;
  std::vector<real_t> Temperature;
  std::vector<real_t> phaseMarker;  
  std::vector<real_t> DBContainer;
  std::vector<real_t> U1_3phiContainer;
  
  // std::vector<real_t> trA_reContainer, trA_imContainer;
  std::vector<real_t> u11Container, u12Container, u13Container,
                      u21Container, u22Container, u23Container,
                      u31Container, u32Container, u33Container;
  std::vector<real_t> v11Container, v12Container, v13Container,
                      v21Container, v22Container, v23Container,
                      v31Container, v32Container, v33Container;
  // std::vector<real_t> eigAv1Container, eigAv2Container, eigAv3Container;
  std::vector<real_t> jm1Container, jm2Container, jm3Container;
  std::vector<real_t> lsqContainer, l1Container, l2Container, l3Container;
  std::vector<real_t> GPhi1Container, GPhi2Container, GPhi3Container;  

  std::vector<real_t> js11Container, js21Container, js31Container,
                      js12Container, js22Container, js32Container,
                      js13Container, js23Container, js33Container;
  
  /*--------------------------------*/
  
  long long ghostVolume;
  long long latticeVolumeWithGhost;
  long long latticeVolume;
  
  unsigned char *ghostCellsMask;

  ascent::Ascent pio;
  conduit::Node pio_options;
  conduit::Node actions;
  conduit::Node mesh;
  
  /*----------------------------------------*/
  /*      pario declearations end here      */
  /*----------------------------------------*/
  
};

#endif
