#ifndef GLSOL_HPP
#define GLSOL_HPP

#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "matep.hpp"
#include "dyGiLa_config.hpp"

// Definition of the field that we will use
using real_t = float;                          // or double ?
using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time

// Define convenient enum for addressing the reductions - could use numbers too

/* enumeration counter type  for matrix reduction */
// enum class Matreduc { i_sumA, N_MatREDUCTION };
enum class matreduc { i_sumA, N_MatREDUCTION };

/* enumeration counter type for phase marker reduction */
//enum class pxacc {
enum {
  p0_acc, p1_acc, p2_acc, p3_acc, p4_acc, 
  p5_acc, p6_acc, p7_acc, p8_acc, p9_acc,
  N_PMREDUCTION
};  

/* enumeration counter type for scalar reduction */
//enum class reduc {
enum {
  i_sumgapA,
  i_suma, i_sumb1, i_sumb2, i_sumb3, i_sumb4, i_sumb5,
  i_suma_we, i_sumb1_we, i_sumb2_we, i_sumb3_we, i_sumb4_we, i_sumb5_we,
  i_sumk1, i_sumk2, i_sumk3, i_sumk1_we, i_sumk2_we, i_sumk3_we,
  i_sumkin, i_sumkin_we,
  N_REDUCTION
};  

// D_B value at masked sites (gap below invariants_gapA_min); D_B itself is >= 0
#define DB_MASKED (-1.0)

// Container for simulation parameters and methods
class glsol{

public:
  // glsol() = default;                     // default constructor
  glsol();
  
  // read configration file and initiate scaling_sim.config 
  const std::vector<std::string> configure(const std::string &fname, int argc, char **argv);

  // OP field initialization
  void initialize();

  // T-field initialization  
  void initializeT();

  // p-field initialization
  // void initializep();
  
  // H-field initialization    
  void initializeH();

  void fstreams_open(const std::vector<std::string> &);
  void fstreams_close();  

  void write_energies();
  void write_positions();
  void write_phases();

  // Gaussian Low Pass filters
  void GaussianLPfilter_matrix(/*Field<phi_t> &*/);
  
  void phaseMarking();
  void phaseCounting();
  
  void next();
  void next_bath();
  void next_AdGRz_bath();  
  void next_bath_UniT_quench();
  void next_bath_UniT_quench_Hfield();
  void next_bath_UniT_quench_AdGRz_Hfield();
  void next_bath_UniT_quench_AdGRz_Hfield_confCatch();    
  void next_bath_hotblob_quench_Hfield();
  void next_bath_Quasi2Dhotblob_quench_AdGR_Hfield();  
  void next_bath_hotblob_quench_Hfield_confCatch();        
  
  Field<phi_t> A, AwT; // OP field A, & its Wirsterass transformation
  Field<phi_t> pi;
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  Field<real_t> T;
  Field<real_t> dT;
  // Field<real_t> dT_from_local_TAB;  
  // Field<real_t> p;

  Field<Vector<3,real_t>> H; // H-field, 3-component column vector field
  Field<real_t> phaseMarker;
  Field<real_t> DB; // D_B = sum_a (R_a - R_a^B)^2, set in phaseMarking()
  
  real_t t = 0., pM5VolR = 0., pM9VolR = 0.;
  real_t tc = 0.;
  

  unsigned int extinguish_t = 0.;  
  matep::Matep MP;
  dyGiLaConf config{}; // value-initialised: parameters not set by configure() are 0

private:
  
  /* utils */
  void dPiGLfe();
  void dPiGLfe_AdGRz();  
  void ABOBA();
  void ABOBA_gBranch();
  void dampAndRelax();
  void relax();
  //void AdGRzTreat(phi_t &, phi_t &, phi_t &, phi_t &, const real_t &);

  /* utils conf init */
  void case_0(), case_1(), case_2();
  void case_3(), case_4(), case_5();
  void case_6(), case_7(), case_8();
  void case_9(), case_10();
  
    
};

#endif
