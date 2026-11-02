#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"


const std::vector<std::string> glsol::configure(const std::string &fname, int argc, char **argv) {
  
    hila::input parameters(fname);
    hila::out0 << "------------------------------------------------------------" << "\n"
               << "-- dyGiLa 3D p-Wave TDGL HPC Simulation Suites Parameters --" << "\n"
               << "------------------------------------------------------------" << std::endl;

    const std::vector<unsigned int> lattice_sizes = parameters.get("lattice_sizes_array");
    config.lx = lattice_sizes[0];
    config.ly = lattice_sizes[1];
    config.lz = lattice_sizes[2];    
    config.dx = parameters.get("dx");
    config.dtdxRatio = parameters.get("dtdxRatio");
    config.dt = config.dx * config.dtdxRatio;
    
    /*********************************************/
    /* >>>>>>  Homogenous quench parameters <<<< */    
    /*********************************************/    
    config.tThermalizationWaiting = parameters.get("tThermalizationWaiting");
    config.tauQ1                  = parameters.get("tauQ1");
    config.tauQ2                  = parameters.get("tauQ2");    
    config.has1stQStop            = parameters.get_item("has1stQStop",{"no", "yes"});
    if (config.has1stQStop == 1)
      {
	const std::vector<real_t> temp = parameters.get("TtdQ1st_tQ1Wainting_array");
	config.Ttd_Q1st = temp[0]; config.tQ1Waiting = temp[1];	
	// config.Ttd_Q1st = parameters.get("Ttd_Q1st");
	// config.tQ1Waiting = parameters.get("tQ1Waiting");
      }
    config.use_antiQuench         = parameters.get_item("use_antiQuench",{"no", "yes"});
    config.Ttd_Qend = parameters.get("Ttd_Qend");

    /*********************************************/
    /* > Homogenous quench parameters end here<< */    
    /*********************************************/    
        
    config.tStart = parameters.get("tStart");
    config.tEnd = parameters.get("tEnd");
    config.tdis = parameters.get("tdis");
    
    
    /******************************************/
    /*----   gamma as a complex number   -----*/
    /*   and gamma changes at certain stages  */
    /*   gamma1 and gamma2 are used to hold   */
    /*   different values of gamma.           */
    /*   gamma change is triggered in main()  */
    /******************************************/

    // switch for turning on the T-dependent gamma
    config.TDependnetgamma = parameters.get_item("TDependnetgamma",{"no", "yes"});

    // switch for shifting the T-dependent gamma base line
    // This is usful when T-dep gamma is too small at low Temperature
    config.shiftTDependent_gamma_td = parameters.get_item("shiftTDependent_gamma_td",{"no", "yes"});
    if (config.shiftTDependent_gamma_td == true)
      { config.gamma_td_BaseLine = parameters.get("gamma_td_BaseLine"); }
    
    //config.gamma = parameters.get("gamma");
    std::vector<real_t> tmp1 = parameters.get("gamma1");
    std::vector<real_t> tmp2 = parameters.get("gamma2");    

    if (tmp1.size() < 1 || tmp1.size() >2)
      {
	hila::out0 << "error: gamma1 must be initialized by at least one real or imagnary parts" << std::endl;
	hila::error("\n");
      }
    else if (tmp1.size() == 1)
      {
	hila::out0 << " tmp1.size() = " << tmp1.size() << "\n";
        config.gamma1.real() = tmp1[0];
        config.gamma1.imag() = 0.;
      }
    else  // tmp.size() == 2
      {
        hila::out0 << " tmp1.size() = " << tmp1.size() << "\n";
        config.gamma1.real() = tmp1[0];
        config.gamma1.imag() = tmp1[1];
      }

    if (tmp2.size() < 1 || tmp2.size() >2)
      {
	hila::out0 << "error: gamma must be initialized by at least one real or imagnary parts" << std::endl;
	hila::error("\n");
      }
    else if (tmp2.size() == 1)
      {
	hila::out0 << " tmp2.size() = " << tmp2.size() << "\n";
        config.gamma2.real() = tmp2[0];
        config.gamma2.imag() = 0.;
      }
    else  // tmp.size() == 2
      {
        hila::out0 << " tmp2.size() = " << tmp2.size() << "\n";
        config.gamma2.real() = tmp2[0];
        config.gamma2.imag() = tmp2[1];
      }

    // conunter for triggering gamma value change 
    config.gammaoffc = parameters.get("gammaoffc");
    /*
     * gamma2 appling time, started from first moment when blob confCatch() block is entred.
     * This for blob confCatch function noise removment. After this time, T-dependent gamma 
     * is recovery without thermal noise. (in unit of tGL) 
     */ 
    config.extinguish_off_t_count = parameters.get("extinguish_off_t_count");
    
    /******************************************/
    /*        gamma as complex number         */
    /*  with gamma1(2) reading end here       */
    /******************************************/        

    config.initialCondition = parameters.get_item("initialCondition",{"gaussrand"             //0
								      ,"kgaussrand"           //1
								      ,"normal_phase_real1"   //2
								      ,"normal_phase_real2"   //3
								      ,"normal_phase_complex" //4
								      ,"Bphase"               //5
								      ,"Aphase_partial1"      //6
								      ,"Aphase_full"          //7
	                                                              ,"hotblob"              //8
                                                                      ,"A-n-B"                //9
								      ,"Wiman2016StripeB"});  //10
                                                                       

    hila::out0 << "--config.initialCondition is "
	       << config.initialCondition
	       << std::endl;
    if (// ---------------------------
	// kGaussrand initialization
        // ---------------------------	
	config.initialCondition == 1)
      {
	const std::vector<real_t> temp = parameters.get("case1_params_array");
        config.IniMod = parameters.get("IniMod");
        config.Inilc = parameters.get("Inilc");
      }    

    config.variance_sigma = parameters.get("sigma");    
    config.seed = parameters.get("seed");

    //initialCondition-T
    config.initialConditionT = parameters.get_item("initialConditionT",{"constant","sine","hotblob"});
    if(config.initialConditionT == 0)
      {
	config.IniT = parameters.get("IniT");
      }
    else if (config.initialConditionT == 1)
      {
	const std::vector<real_t> temp = parameters.get("IniT_ampT_array");
	config.IniT = temp[0];
	config.ampT = temp[1];	
	// config.IniT = parameters.get("IniT");
	// config.ampT = parameters.get("ampT");
      }
    else if (config.initialConditionT == 2)
      {
       /*********************************************/
       /* >>>> spherical hot blob parameters <<<<<  */    
       /*********************************************/
	const std::vector<real_t> temp = parameters.get("HotBlob_params_array");	
        config.Ttdb1         = temp[0];
        config.Ttdb0         = temp[1];
        config.t1            = temp[2];
	config.Blob_Tc_cutoff = parameters.get_item("Blob_Tc_cutoff",{"no", "yes"});
        config.use_CustomerDctxi = parameters.get_item(" use_CustomerDctxi",{"no", "yes"});
        if ( config.use_CustomerDctxi == true )
	  {
            config.Dctxi             = parameters.get("Dctxi");		    
	    hila::out0 << "--use Custom Dctxi. " << std::endl;
	  }

       /*********************************************/
       /* >>>>>> hot bloob parameters end here  <<< */    
       /*********************************************/    	
      }
    
    //initialCondition-p
    config.initialConditionp = parameters.get_item("initialConditionp",{"constant"});
    if (config.initialConditionp == 0) { config.Inip = parameters.get("Inip"); }

    //initialCondition-Hfield
    config.withHfield = parameters.get_item("withHfield",{"no", "yes"});
    if (config.withHfield == 1)
      {
       config.initialConditionH = parameters.get_item("initialConditionH",{"constant"});
       std::vector<real_t> temp	= parameters.get("InitH");
       foralldir(al){ config.InitH.e(al) = temp[al]; }
      }
    else
      {
       // H enters dPiGLfe() unconditionally, so set H = 0 explicitly
       config.initialConditionH = 0;
       foralldir(al){ config.InitH.e(al) = 0.0; }
      }
        
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("nOutputs");

    // pStreaming Squezeing Ratio; PSSRatio, the pIO steps skip ratio
    // This control the frequency which pario.pStream() is called;
    // minimum is 1, which means not steps skip, must be unsigned int
    // together with hdf5start/end, one has Precious control on pIO.
    config.PSSRatio = parameters.get("PSSRatio");

    // output_file is the saving path of output file, which offered in congigration file
    // vector for return by names of stream files
    const std::string output_file    = parameters.get("output_file");
    const std::string pVcounter_file = parameters.get("pVcounter_file");

    std::vector<std::string> name_files{output_file, pVcounter_file};

    // xdmf file name, which is provided through config file
    config.xmf2_Amatrix_fname = parameters.get("xmf2_Amatrix_fname");
    config.xmf2_pMarker_fname = parameters.get("xmf2_pMarker_fname");
    config.xmf2_massCurrent_fname = parameters.get("xmf2_massCurrent_fname");
    config.xmf2_spinCurrent_fname = parameters.get("xmf2_spinCurrent_fname");    

    /*----------------------------------------*/
    /* >>>>>>>  boundary conditions  <<<<<<<<<*/
    /*----------------------------------------*/
    config.BCs1 = parameters.get_item("BCs1",{"periodic",       //0
					      "AB",             //1
					      "PairBreaking",   //2
					      "PB_y",           //3
					      "PairB_yz",       //4
					      "BB",             //5
					      "AdGRz",          //6
					      "phaseVortices",  //7
					      "A-n-B"});        //8
    
    config.BCs2 = parameters.get_item("BCs2",{"periodic",       //0
					      "AB",             //1
					      "PairBreaking",   //2
                                              "PB_y",           //3
                                              "PairB_yz",       //4 
                                              "BB",             //5
					      "AdGRz",          //6
                                              "phaseVortices",  //7
                                              "A-n-B"});        //8
    
    config.BCchangec = parameters.get("BCchangec");
    config.use_AdGRz_surfaces = parameters.get_item("use_AdGRz_surfaces",{"no","yes"});
    if ( config.use_AdGRz_surfaces == 1 ) { config.bt = parameters.get("bt_AdGR"); }
    // config.Wn = parameters.get("BoundaryPhaseWindingNO");
                   
    config.evolveT = parameters.get_item("evolveT",{"no","yes"});
    if(config.evolveT == 1
       && config.initialConditionT != 2 )
      {
	config.Tevolvetype = parameters.get_item("Tevolvetype",{"heat","wave","homogeneousQuench"});
	if (config.Tevolvetype == 0 || config.Tevolvetype == 1)
	  {
	   config.startdiffT = parameters.get("startdiffT");
	   config.diffT = parameters.get("diffT");
	  }
      }

    config.useTbath = parameters.get_item("useTbath",{"no","yes"});
    config.Tbath_start = parameters.get("Tbath_start");

    config.ptol = parameters.get("ptol");

    /*----------------------------------------*/    
    /* Approx. Gaussian LP filter parameters  */
    /*----------------------------------------*/
    config.useGaussianLP_filter = parameters.get_item("useGaussianLP_filter",{"no","yes"});
    if (config.useGaussianLP_filter == 1)
      {
       config.numIterGLPfilter = parameters.get("numIterGLPfilter");
       // config.GLPfc1 = parameters.get("GLPfc1");
       config.GLPfc2 = parameters.get("GLPfc2");
       /* 3d finite difference coefficients relation */
       config.GLPfc1 = 1. - 6.*config.GLPfc2; 
      }
    
    /*----------------------------------------*/
    /* Parallel IO Engine control parameters  */
    /*----------------------------------------*/
    config.pario_compute_gapA    = parameters.get_item("pario_compute_gapA",{"no","yes"});    
    config.pario_compute_feDensity     = parameters.get_item("pario_compute_feDensity",{"no","yes"});
    config.pario_Temperature_pStream     = parameters.get_item("pario_Temperature_pStream",{"no","yes"});
    config.pario_compute_phaseMarker = parameters.get_item("pario_compute_phaseMarker",{"no","yes"});
    config.pario_compute_DB          = parameters.get_item("pario_compute_DB",{"no","yes"});
    if (config.pario_compute_phaseMarker == 1 || config.pario_compute_DB == 1)
      { config.invariants_gapA_min = parameters.get("invariants_gapA_min"); }
    config.pario_compute_U1Phase     = parameters.get_item("pario_compute_U1Phase",{"no","yes"});
    config.pario_compute_lVector     = parameters.get_item("pario_compute_lVector",{"no","yes"});
    config.pario_compute_GPhiVector  = parameters.get_item("pario_compute_GPhiVector",{"no","yes"});
    //config.pario_compute_massCurr = parameters.get_item("pario_compute_massCurr",{"no","yes"});
    config.pario_compute_spinCurr = parameters.get_item("pario_compute_spinCurr",{"no","yes"});
    
    
    config.U1_3phi_DetA_ZeroTol = parameters.get("U1_3phi_DetA_ZeroTol");
    config.lVec_SqlTol = parameters.get("lVec_SqlTol");        
    
    config.hdf5_A_matrix_output        = parameters.get_item("hdf5_A_matrix_output",{"no","yes"});
    config.hdf5_AMatrix_exaslice_output        = parameters.get_item("hdf5_AMatrix_exaslice_output",{"no","yes"});        
    config.hdf5_mass_current_output = parameters.get_item("hdf5_mass_current_output",{"no","yes"});
    //config.hdf5_massCurr_exaslice_output = parameters.get_item("hdf5_massCurr_exaslice_output",{"no","yes"});
    
    config.hdf5_spin_current_output = parameters.get_item("hdf5_spin_current_output",{"no","yes"});
    config.hdf5_spinCurr_exaslice_output = parameters.get_item("hdf5_spinCurr_exaslice_output",{"no","yes"});
    
    config.hdf5_pMarker_output         = parameters.get_item("hdf5_pMarker_output",{"no","yes"});
    config.hdf5_lVector_output         = parameters.get_item("hdf5_lVector_output",{"no","yes"});
    config.hdf5_GradPhiVector_output   = parameters.get_item("hdf5_GradPhiVector_output",{"no","yes"});
    config.hdf5_GradPhiVector_exaslice_output   = parameters.get_item("hdf5_GradPhiVector_exaslice_output",{"no","yes"});
    config.hdf5_GradPhiVector_clip_output   = parameters.get_item("hdf5_GradPhiVector_clip_output",{"no","yes"});
    config.hdf5_GradPhiVector_boxclip_output   = parameters.get_item("hdf5_GradPhiVector_boxclip_output",{"no","yes"});

    
    if (
        (config.hdf5_A_matrix_output == 1)
	|| (config.hdf5_AMatrix_exaslice_output == 1)
	|| (config.hdf5_mass_current_output == 1)
	// || (config.hdf5_massCurr_exaslice_output == 1)	
	|| (config.hdf5_spin_current_output == 1)
	|| (config.hdf5_spinCurr_exaslice_output == 1)	
	|| (config.hdf5_pMarker_output == 1)
	|| (config.hdf5_lVector_output == 1)
	|| (config.hdf5_GradPhiVector_output == 1)
	|| (config.hdf5_GradPhiVector_exaslice_output == 1)
	|| (config.hdf5_GradPhiVector_clip_output == 1)
	|| (config.hdf5_GradPhiVector_boxclip_output == 1)			
       )
      {
	config.hdf5_cheerypick = parameters.get_item("hdf5_cheerypick",{"no","yes"});
	if ( config.hdf5_cheerypick == true )
	  {
	    const std::vector<unsigned int> temp = parameters.get("hdf5_cheerypick_array");
	    config.hdf5_cheerypick_array.assign(temp.begin(), temp.end());
   	    hila::out0 << "--do Cheerypick on hdf5. " << std::endl;
	  } // cheerypick params
	else
	  {	    
            config.hdf5Ststart = parameters.get("hdf5Ststart");
            config.hdf5Stend   = parameters.get("hdf5Stend");       
	  } // streming pip params		
      }
    
    config.do_gapA_clip         = parameters.get_item("do_gapA_clip",{"no","yes"});
    if ( config.do_gapA_clip ==1 )
      {
	const std::vector<real_t> temp = parameters.get("gapAClip_params_array");
        config.gapA_clip1_point_x = temp[0];
	config.gapA_clip1_point_y = temp[1];
	config.gapA_clip1_point_z = temp[2];
        config.gapA_clip1_norm_x = temp[3];
	config.gapA_clip1_norm_y = temp[4];
	config.gapA_clip1_norm_z = temp[5];
		
      } // gapA clip control parammeters

    config.do_gapA_slice         = parameters.get_item("do_gapA_slice",{"no","yes"});
    if ( config.do_gapA_slice ==1 )
      {
	const std::vector<real_t> temp = parameters.get("gapASlice_params_array");	
        config.gapA_slice1_point_x = temp[0];
	config.gapA_slice1_point_y = temp[1];
	config.gapA_slice1_point_z = temp[2];
        config.gapA_slice1_norm_x = temp[3];
	config.gapA_slice1_norm_y = temp[4];
	config.gapA_slice1_norm_z = temp[5];	
      } // gapA slice control parammeters

    
    config.do_fed_clip         = parameters.get_item("do_fed_clip",{"no","yes"});
    if ( config.do_fed_clip ==1 )
      {
	const std::vector<real_t> temp = parameters.get("fedClip_params_array");
        config.fed_clip_point_x = temp[0];
	config.fed_clip_point_y = temp[1];
	config.fed_clip_point_z = temp[2];
        config.fed_clip_norm_x = temp[3];
	config.fed_clip_norm_y = temp[4];
	config.fed_clip_norm_z = temp[5];	
      } // gapA clip control parammeters
    
    
    config.do_gapA_isosurface   = parameters.get_item("do_gapA_isosurface",{"no","yes"});
    if ( config.do_gapA_isosurface == 1 )
      {
        config.iso_values_vector = parameters.get("gapAIsoValues_array");
      } // gapA clip control parammeters

    config.do_Temperature_clip = parameters.get_item("do_Temperature_clip",{"no", "yes"});
    if (config.do_Temperature_clip ==1)
      {
	const std::vector<real_t> temp = parameters.get("TempClip_params_array");
        config.Temperature_clip_point_x = temp[0];
	config.Temperature_clip_point_y = temp[1];
	config.Temperature_clip_point_z = temp[2];
        config.Temperature_clip_norm_x = temp[3];
	config.Temperature_clip_norm_y = temp[4];
	config.Temperature_clip_norm_z = temp[5];
	config.Temperature_clamp       = temp[6];
      }

    config.do_Temperature_slice = parameters.get_item("do_Temperature_slice",{"no", "yes"});
    if (config.do_Temperature_slice ==1)
      {
	const std::vector<real_t> temp = parameters.get("TempSlice_params_array");
        config.Temperature_slice_point_x = temp[0];
	config.Temperature_slice_point_y = temp[1];
	config.Temperature_slice_point_z = temp[2];
        config.Temperature_slice_norm_x = temp[3];
	config.Temperature_slice_norm_y = temp[4];
        config.Temperature_slice_norm_z = temp[5];
      }
    
    config.do_Temperature_isosurface = parameters.get_item("do_Temperature_isosurface",{"no","yes"});
    if ( config.do_Temperature_isosurface ==1 )
      {
	std::vector<real_t> tmp3 = parameters.get("TempIsoValues_params_array");
	real_t TcpmK = MP.Tcp_mK(config.Inip);
        for (auto i : tmp3) { config.Temperature_iso_values_vector.push_back(i * TcpmK); }
      } // gapA clip control parammeters

    config.do_phaseMarker_slice1 = parameters.get_item("do_phaseMarker_slice1",{"no","yes"});
    if (config.do_phaseMarker_slice1 == 1)
      {
	const std::vector<real_t> temp = parameters.get("pMarkerSlice1_params_array");
        config.pMarker_slice1_point_x = temp[0];
	config.pMarker_slice1_point_y = temp[1];
	config.pMarker_slice1_point_z = temp[2];
        config.pMarker_slice1_norm_x = temp[3];
	config.pMarker_slice1_norm_y = temp[4];
	config.pMarker_slice1_norm_z = temp[5];
      }

    config.do_phaseMarker_slice2 = parameters.get_item("do_phaseMarker_slice2",{"no","yes"});
    if (config.do_phaseMarker_slice2 == 1)
      {
	const std::vector<real_t> temp = parameters.get("pMarkerSlice2_params_array");
        config.pMarker_slice2_point_x = temp[0];
	config.pMarker_slice2_point_y = temp[1];
	config.pMarker_slice2_point_z = temp[2];
        config.pMarker_slice2_norm_x = temp[3];
	config.pMarker_slice2_norm_y = temp[4];
	config.pMarker_slice2_norm_z = temp[5];	
      }
    
    config.do_phaseMarker_isosurface = parameters.get_item("do_phaseMarker_isosurface",{"no","yes"});
    if (config.do_phaseMarker_isosurface == 1)
      {
	std::vector<real_t> tmp4 = parameters.get("pMarkerIsoValues_params_array");
        for (auto i : tmp4) { config.phaseMarker_iso_values_vector.push_back(i); }
      }

    config.do_DB_slice = parameters.get_item("do_DB_slice",{"no","yes"});
    if (config.do_DB_slice == 1)
      {
	const std::vector<real_t> temp = parameters.get("DBSlice_params_array");
	config.DB_slice_point_x = temp[0];
	config.DB_slice_point_y = temp[1];
	config.DB_slice_point_z = temp[2];
	config.DB_slice_norm_x = temp[3];
	config.DB_slice_norm_y = temp[4];
	config.DB_slice_norm_z = temp[5];
	const std::vector<real_t> range = parameters.get("DB_range_array");
	config.DB_range_min = range[0];
	config.DB_range_max = range[1];
	// quoted, comma-separated Ascent/Viskores colour table names
	const std::vector<std::string> tables = parameters.get("DB_color_tables");
	config.DB_color_tables = tables;
	// RGB in [0,1] for sites masked by invariants_gapA_min
	const std::vector<real_t> rgb = parameters.get("DB_mask_rgb");
	for (int c = 0; c < 3; c++) config.DB_mask_rgb[c] = rgb[c];
      }

    config.remove_screen_annotations = parameters.get_item("remove_screen_annotations",{"no","yes"});
    config.do_phaseMarker_fieldclip = parameters.get_item("do_phaseMarker_fieldclip",{"no","yes"});
    config.do_phaseMarker_fieldclip_Bphase = parameters.get_item("do_phaseMarker_fieldclip_Bphase",{"no","yes"});
    config.do_phaseMarker_fieldclip_Aphase = parameters.get_item("do_phaseMarker_fieldclip_Aphase",{"no","yes"});    

    config.do_U13phi_slice = parameters.get_item("do_U13phi_slice",{"no","yes"});
    if (config.do_U13phi_slice == 1)
      {
	std::vector<real_t> temp = parameters.get("U13phiSlice_params_array");
        config.U13phi_slice_point_x = temp[0];
	config.U13phi_slice_point_y = temp[1];
	config.U13phi_slice_point_z = temp[2];
        config.U13phi_slice_norm_x = temp[3];
	config.U13phi_slice_norm_y = temp[4];
	config.U13phi_slice_norm_z = temp[5];	
      }

    config.do_l_sq_slice = parameters.get_item("do_l_sq_slice",{"no","yes"});
    if (config.do_l_sq_slice == 1)
      {
	std::vector<real_t> temp = parameters.get("lsqSlice_params_array");
        config.l_sq_slice_point_x = temp[0];
	config.l_sq_slice_point_y = temp[1];
	config.l_sq_slice_point_z = temp[2];
        config.l_sq_slice_norm_x = temp[3];
	config.l_sq_slice_norm_y = temp[4];
	config.l_sq_slice_norm_z = temp[5];	
      }

    config.do_GPhi_slice_extract = parameters.get_item("do_GPhi_slice_extract",{"no","yes"});
    if (config.do_GPhi_slice_extract == 1)
      {
        std::vector<real_t> exaslice_params = parameters.get("GPhiSlice_params_array");        	
        config.GPhi_slice1_point_x = exaslice_params[0];
	config.GPhi_slice1_point_y = exaslice_params[1];
	config.GPhi_slice1_point_z = exaslice_params[2];
        config.GPhi_slice1_norm_x = exaslice_params[3];
	config.GPhi_slice1_norm_y = exaslice_params[4];
	config.GPhi_slice1_norm_z = exaslice_params[5];
      }

    config.do_GPhi_clip_extract = parameters.get_item("do_GPhi_clip_extract",{"no","yes"});
    if (config.do_GPhi_clip_extract == 1)
      {
        std::vector<real_t> clip_params = parameters.get("GPhiClip_params_array");        	
        config.GPhi_clip1_point_x = clip_params[0];
	config.GPhi_clip1_point_y = clip_params[1];
	config.GPhi_clip1_point_z = clip_params[2];
        config.GPhi_clip1_norm_x = clip_params[3];
	config.GPhi_clip1_norm_y = clip_params[4];
	config.GPhi_clip1_norm_z = clip_params[5];
      }

    config.do_GPhi_boxclip_extract = parameters.get_item("do_GPhi_boxclip_extract",{"no","yes"});
    if (config.do_GPhi_boxclip_extract == 1)
      {
        std::vector<real_t> clip_params = parameters.get("GPhiBoxClip_params_array");        	
        config.GPhi_boxclip1_min_x = clip_params[0];
	config.GPhi_boxclip1_min_y = clip_params[1];
	config.GPhi_boxclip1_min_z = clip_params[2];
        config.GPhi_boxclip1_max_x = clip_params[3];
	config.GPhi_boxclip1_max_y = clip_params[4];
	config.GPhi_boxclip1_max_z = clip_params[5];
      }

    config.do_AMatrix_slice_extract = parameters.get_item("do_AMatrix_slice_extract",{"no","yes"});
    if (config.do_AMatrix_slice_extract == 1)
      {
        std::vector<real_t> exaslice_params = parameters.get("AMSlice_params_array");        	
        config.AM_slice1_point_x = exaslice_params[0];
	config.AM_slice1_point_y = exaslice_params[1];
	config.AM_slice1_point_z = exaslice_params[2];
        config.AM_slice1_norm_x = exaslice_params[3];
	config.AM_slice1_norm_y = exaslice_params[4];
	config.AM_slice1_norm_z = exaslice_params[5];
      }

    config.do_spinCurr_slice_extract = parameters.get_item("do_spinCurr_slice_extract",{"no","yes"});
    if (config.do_spinCurr_slice_extract == 1)
      {
        std::vector<real_t> exaslice_params = parameters.get("spinCurrSlice_params_array");        	
        config.spinCurr_slice1_point_x = exaslice_params[0];
	config.spinCurr_slice1_point_y = exaslice_params[1];
	config.spinCurr_slice1_point_z = exaslice_params[2];
        config.spinCurr_slice1_norm_x = exaslice_params[3];
	config.spinCurr_slice1_norm_y = exaslice_params[4];
	config.spinCurr_slice1_norm_z = exaslice_params[5];
      }

    
    
    // config.do_gapA_3slice       = parameters.get_item("do_gapA_3slice",{"no","yes"});
    // config.do_fe_slice          = parameters.get_item("do_fe_slice",{"no","yes"});
    // config.do_gapA_slice        = parameters.get_item("do_gapA_slice",{"no","yes"});            

    //const std::vector<real_t> gapAClampBias_params = parameters.get("gapAClampBias_params_array");
    {
     // config.clamp_bias_gapMin = gapAClampBias_params[0];
     // config.clamp_bias_gapMax = gapAClampBias_params[1];
     
     config.clamp_bias_gapMin = parameters.get("clamp_bias_gapMin");
     config.clamp_bias_gapMax = parameters.get("clamp_bias_gapMax"); 
    }

    //const std::vector<real_t> fedClampBias_params = parameters.get("fedClampBias_params_array");
    {
     // config.clamp_bias_fed_Min = fedClampBias_params[0];
     // config.clamp_bias_fed_Max = fedClampBias_params[1];

     config.clamp_bias_fed_Min = parameters.get("clamp_bias_fed_Min");
     config.clamp_bias_fed_Max = parameters.get("clamp_bias_fed_Max");      
    }

    //const std::vector<real_t> ColorBarPosition_params = parameters.get("ColorBarPosition_params_array");
    {
     // config.CBxMin = ColorBarPosition_params[0];
     // config.CBxMax = ColorBarPosition_params[1];
     // config.CByMin = ColorBarPosition_params[2];
     // config.CByMax = ColorBarPosition_params[3];
    
     config.CBxMin = parameters.get("CBxMin");
     config.CBxMax = parameters.get("CBxMax");
     config.CByMin = parameters.get("CByMin");
     config.CByMax = parameters.get("CByMax");    
    }

    // const std::vector<real_t> ImageSizes_params = parameters.get("ImageSizes_params_array");
    {
     // config.image_width1 = ImageSizes_params[0];
     // config.image_height1 = ImageSizes_params[1];
     // config.image_width2 = ImageSizes_params[2];
     // config.image_height2 = ImageSizes_params[3];

     config.image_width1 = parameters.get("image_width1");
     config.image_height1 = parameters.get("image_height1");
     config.image_width2 = parameters.get("image_width2");
     config.image_height2 = parameters.get("image_height2");     
    }

    // const std::vector<real_t> CameraPosition_params = parameters.get("CameraPosition_params_array");
    {
      // config.camera1_azi = CameraPosition_params[0];
      // config.camera1_ele = CameraPosition_params[1];
      // config.camera2_azi = CameraPosition_params[2];
      // config.camera2_ele = CameraPosition_params[3];
      
     config.camera1_azi = parameters.get("camera1_azi");
     config.camera1_ele = parameters.get("camera1_ele");
     config.camera2_azi = parameters.get("camera2_azi");
     config.camera2_ele = parameters.get("camera2_ele");
    }

    // const std::vector<real_t> CameraZoom_params = parameters.get("CameraZoom_params_array");
    {
      // config.zoom1 = CameraZoom_params[0];
      // config.zoom2 = CameraZoom_params[1];

      config.zoom1 = parameters.get("zoom1");
      config.zoom2 = parameters.get("zoom2");;
    }
    
    /*----------------------------------------*/
    /* Parallel IO Engine parameters end      */
    /*----------------------------------------*/

    t = config.tStart;    
    return name_files; //output_file pVfile;
    
} // allocate() function ends here

