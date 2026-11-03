#ifndef DYGILA_CONFIG_HPP
#define DYGILA_CONFIG_HPP

//#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "plumbing/hila.h"

using real_t = float;                     

// parameters configuraion class of dyGiLa
struct dyGiLaConf {
  
      int lx;
      int ly;
      int lz;
      real_t dx;
      real_t dt;
      real_t dtdxRatio;

      // tThermalizationWaiting, initilial configuration thermalizing time
      // quench time, measured from Tc to 0
      // Ttd_Qend, detinated temperature Ttd_Qend = T_Qend/Tc
      // has1stQStop, has 1st quench stop point
      // Ttd_Q1st, 1st quench stop Temp, Ttd_Q1st = T_Q1st/Tc
      // tQ1Waiting, waiting time of 1st quench stop Temp
      real_t tThermalizationWaiting;
      real_t tauQ1;
      real_t tauQ2;  
      real_t tauQ;
      unsigned int has1stQStop;
      real_t Ttd_Q1st;
      real_t tQ1Waiting;
      unsigned int use_antiQuench;
      real_t Ttd_Qend;
  
      real_t tStart;
      real_t tEnd;
      real_t tdis;

      //real_t gamma;
      unsigned int TDependnetgamma;
      unsigned int shiftTDependent_gamma_td;
      real_t gamma_td_BaseLine = 0.f;
      Complex<real_t> gamma;
      Complex<real_t> gamma1;
      Complex<real_t> gamma2;
      int gammaoffc;

      /* gamma2 appling time, started from 
       * first moment when blob confCatch() block is entred.
       * in unit of tGL. 
       */
      unsigned int extinguish_off_t_count;
      
      int initialCondition;
      real_t variance_sigma;
      
      int seed;
      real_t IniMod;
      real_t Inilc;

      int item;
        
      // real_t T;
      // real_t dT_from_TAB;
      // real_t p;
    
      real_t alpha;
      real_t beta1;
      real_t beta2;
      real_t beta3;
      real_t beta4;
      real_t beta5;

      // T-filed initialing & control parameters
      int initialConditionT;
      real_t IniT;
      real_t ampT; //sine wave profile
      /*-----------------------------------------*/
      /* hot blob Temeprature profile parameters */
      /* Ttdbx mean T in unit of Tc              */
      bool Blob_Tc_cutoff;
      real_t Ttdb1;
      real_t Ttdb0;
      real_t t1;
      /*-----------------------------------------*/
      // Switch for using customer Diffusivity
      unsigned int use_CustomerDctxi;
      // Diffusivity in unit of xi0GL^2. tGL^-1
      real_t Dctxi;     

      // p-filed initialing & control parameters
      int initialConditionp;
      real_t Inip;

      // H-filed initialing & control parameters
      unsigned int withHfield;  
      int initialConditionH;
      Vector<3, real_t> InitH;
  
      real_t tStats;
      real_t nOutputs;

      // pStreaming Squezeing Ratio; PSSRatio, the pIO steps skip ratio
      unsigned int PSSRatio;

      std::fstream stream;
      std::fstream streampc;
      
      std::string xmf2_Amatrix_fname, xmf2_pMarker_fname,
	          xmf2_massCurrent_fname, xmf2_spinCurrent_fname;

      int positions;
      int npositionout;
      
      int boundaryConditions;
      int BCs1;
      int BCs2;
      int BCchangec;
      // int Wn;
      unsigned int use_AdGRz_surfaces;
      real_t bt;
      
      int useTbath;
      real_t Tbath_start;

      int write_phases;
      int write_eigen;

      real_t ptol;

      /*----------------------------------------*/
      /* Approx. Gaussian LP filter parameters  */
      /*----------------------------------------*/
    
      unsigned int useGaussianLP_filter;
      unsigned int numIterGLPfilter;
      real_t GLPfc1, GLPfc2;

      /*----------------------------------------*/
      /* Approx. GLP filter parameters ends     */  
      /*----------------------------------------*/
    
      int evolveT;
      int Tevolvetype;
      real_t startdiffT;
      real_t diffT;
      int bloob_after;
      real_t theat;
    
      /*----------------------------------------*/
      /*     parallel IO control parameters     */
      /*----------------------------------------*/
      unsigned int pario_compute_gapA, pario_compute_feDensity,
	           pario_compute_phaseMarker, pario_compute_U1Phase,
	           pario_compute_lVector, pario_compute_GPhiVector,
	           pario_compute_massCurr, pario_compute_spinCurr;
      // D_B = sum_a (R_a - R_a^B)^2, distance from the B phase in quartic invariants
      unsigned int pario_compute_DB;
      // The quartic invariants are undefined where A ~ 0. Sites with
      // sqrt(tr AA^dagger) below this get phaseMarker 0 and D_B = DB_MASKED.
      real_t invariants_gapA_min;
      unsigned int pario_Temperature_pStream;

      real_t U1_3phi_DetA_ZeroTol, lVec_SqlTol;
  
      unsigned int do_gapA_clip;
      real_t gapA_clip1_point_x, gapA_clip1_point_y, gapA_clip1_point_z,
             gapA_clip1_norm_x, gapA_clip1_norm_y, gapA_clip1_norm_z;
      real_t gapA_clip2_point_x, gapA_clip2_point_y, gapA_clip2_point_z,
             gapA_clip2_norm_x, gapA_clip2_norm_y, gapA_clip2_norm_z;

      unsigned int do_gapA_slice;
      real_t gapA_slice1_point_x, gapA_slice1_point_y, gapA_slice1_point_z,
             gapA_slice1_norm_x, gapA_slice1_norm_y, gapA_slice1_norm_z;
      real_t gapA_slice2_point_x, gapA_slice2_point_y, gapA_slice2_point_z,
             gapA_slice2_norm_x, gapA_slice2_norm_y, gapA_slice2_norm_z;
 
  
      unsigned int do_fed_clip;
      real_t fed_clip_point_x, fed_clip_point_y, fed_clip_point_z,
             fed_clip_norm_x, fed_clip_norm_y, fed_clip_norm_z;

      unsigned int do_Temperature_clip;
      real_t Temperature_clip_point_x, Temperature_clip_point_y, Temperature_clip_point_z,
	     Temperature_clip_norm_x, Temperature_clip_norm_y, Temperature_clip_norm_z,
             Temperature_clamp;

      unsigned int do_Temperature_slice;
      real_t Temperature_slice_point_x, Temperature_slice_point_y, Temperature_slice_point_z,
	     Temperature_slice_norm_x, Temperature_slice_norm_y, Temperature_slice_norm_z;
  
  
      unsigned int do_Temperature_isosurface;
      std::vector<double> Temperature_iso_values_vector;

      unsigned int do_gapA_isosurface;
      std::vector<double> iso_values_vector;

      unsigned int do_phaseMarker_slice1, do_phaseMarker_slice2;
      real_t pMarker_slice1_point_x, pMarker_slice1_point_y, pMarker_slice1_point_z,
	     pMarker_slice1_norm_x, pMarker_slice1_norm_y, pMarker_slice1_norm_z,
	     pMarker_slice2_point_x, pMarker_slice2_point_y, pMarker_slice2_point_z,
	     pMarker_slice2_norm_x, pMarker_slice2_norm_y, pMarker_slice2_norm_z;

  
      unsigned int do_phaseMarker_isosurface;
      std::vector<double> phaseMarker_iso_values_vector;

      // D_B slice: plane, fixed colour range, one rendering per colour table
      unsigned int do_DB_slice;
      real_t DB_slice_point_x, DB_slice_point_y, DB_slice_point_z,
	     DB_slice_norm_x, DB_slice_norm_y, DB_slice_norm_z;
      real_t DB_range_min, DB_range_max;
      std::vector<std::string> DB_color_tables;
      double DB_mask_rgb[3];   // flat colour of masked (small-gap) sites

      unsigned int remove_screen_annotations;
      unsigned int do_phaseMarker_fieldclip,
	           do_phaseMarker_fieldclip_Bphase,
	           do_phaseMarker_fieldclip_Aphase;

      unsigned int do_U13phi_slice;
      real_t U13phi_slice_point_x, U13phi_slice_point_y, U13phi_slice_point_z,
	     U13phi_slice_norm_x, U13phi_slice_norm_y, U13phi_slice_norm_z;

      unsigned int do_l_sq_slice;
      real_t l_sq_slice_point_x, l_sq_slice_point_y, l_sq_slice_point_z,
	     l_sq_slice_norm_x, l_sq_slice_norm_y, l_sq_slice_norm_z;

      unsigned int do_GPhi_slice_extract, do_GPhi_clip_extract, do_GPhi_boxclip_extract;
      real_t GPhi_slice1_point_x, GPhi_slice1_point_y, GPhi_slice1_point_z,
	     GPhi_slice1_norm_x, GPhi_slice1_norm_y, GPhi_slice1_norm_z,
	     GPhi_boxclip1_min_x, GPhi_boxclip1_min_y, GPhi_boxclip1_min_z,
	     GPhi_boxclip1_max_x, GPhi_boxclip1_max_y, GPhi_boxclip1_max_z,	      	
	     GPhi_clip1_point_x, GPhi_clip1_point_y, GPhi_clip1_point_z,
	     GPhi_clip1_norm_x, GPhi_clip1_norm_y, GPhi_clip1_norm_z;

      unsigned int do_AMatrix_slice_extract;
      real_t AM_slice1_point_x, AM_slice1_point_y, AM_slice1_point_z,
	     AM_slice1_norm_x, AM_slice1_norm_y, AM_slice1_norm_z;

      unsigned int do_spinCurr_slice_extract;
      real_t spinCurr_slice1_point_x, spinCurr_slice1_point_y, spinCurr_slice1_point_z,
	     spinCurr_slice1_norm_x, spinCurr_slice1_norm_y, spinCurr_slice1_norm_z;
  
      unsigned int do_massCurr_slice_extract;  
      real_t massCurr_slice1_point_x, massCurr_slice1_point_y, massCurr_slice1_point_z,
	     massCurr_slice1_norm_x, massCurr_slice1_norm_y, massCurr_slice1_norm_z;
  
      unsigned int hdf5_A_matrix_output, hdf5_AMatrix_exaslice_output;
      //              hdf5_trA_output,
      //              hdf5_eigvA_output,
      unsigned int hdf5_mass_current_output, hdf5_massCurr_exaslice_output,
                   hdf5_spin_current_output, hdf5_spinCurr_exaslice_output,
	           hdf5_pMarker_output,
	           hdf5_lVector_output,
	           hdf5_GradPhiVector_output,
	           hdf5_GradPhiVector_exaslice_output,
	           hdf5_GradPhiVector_clip_output,
	           hdf5_GradPhiVector_boxclip_output;
      unsigned int hdf5_cheerypick = 1;
      std::vector<unsigned int> hdf5_cheerypick_array{0};
      real_t hdf5Ststart = 0., hdf5Stend = 0.;	  
	      
      real_t clamp_bias_gapMin, clamp_bias_gapMax;
      real_t clamp_bias_fed_Min, clamp_bias_fed_Max;

      real_t CBxMin, CBxMax, CByMin, CByMax;
  
      unsigned int image_width1, image_height1,
	           image_width2, image_height2;
      real_t camera1_azi, camera1_ele,
             camera2_azi, camera2_ele;
      real_t zoom1, zoom2;
      /*----------------------------------------*/
      /*    parallel IO parameter end           */
      /*----------------------------------------*/
};
  

#endif
