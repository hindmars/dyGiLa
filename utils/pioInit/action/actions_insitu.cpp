#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cctype>
#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"
#include "pario.hpp"
#include "hsvCyclicSchemeRGBA.hpp"
#include "discreteColorBarScheme.hpp"

namespace {
// phaseMarker 0 (normal: gap below invariants_gapA_min) is drawn as a separate
// mid-grey plot, so the 1-9 colour bar is unchanged. "<pl>" keeps cells with
// no vertex at 0 (select cells with any 0 vertex, then invert); "<pl>zero"
// is the complement.
void addMarkerZeroSplit(conduit::Node &pipelines, const std::string &pl) {
  pipelines[pl + "/f2/type"] = "threshold";
  pipelines[pl + "/f2/params/field"] = "phaseMarker";
  pipelines[pl + "/f2/params/min_value"] = -1.0;
  pipelines[pl + "/f2/params/max_value"] = 0.5;
  pipelines[pl + "/f2/params/invert"] = "true";

  pipelines[pl + "zero/f1"] = pipelines[pl + "/f1"];
  pipelines[pl + "zero/f2/type"] = "threshold";
  pipelines[pl + "zero/f2/params/field"] = "phaseMarker";
  pipelines[pl + "zero/f2/params/min_value"] = -1.0;
  pipelines[pl + "zero/f2/params/max_value"] = 0.5;
}

void addMarkerZeroPlot(conduit::Node &plot, const std::string &pl) {
  static double grey[2] = {0.5, 0.5}, alpha[2] = {1.0, 1.0}, pos[2] = {0.0, 1.0};
  plot["type"] = "pseudocolor";
  plot["pipeline"] = pl + "zero";
  plot["field"] = "phaseMarker";
  plot["min_value"] = -1.0;   // avoid a zero-width range
  plot["max_value"] = 1.0;
  conduit::Node &ct = plot["color_table"];
  ct["control_points/r"].set_external(grey, 2);
  ct["control_points/g"].set_external(grey, 2);
  ct["control_points/b"].set_external(grey, 2);
  ct["control_points/a"].set_external(alpha, 2);
  ct["control_points/position"].set_external(pos, 2);
  ct["annotation"] = "false";   // colour bar for phaseMarker 1-9 only
}

// Common annotation settings for all slice renders: colour bar across the
// top (CBxMin..CByMax), and colour bar and title switched off together by
// remove_screen_annotations.
void setSliceAnnotations(conduit::Node &render, const double *cbpos, const glsol &sol) {
  render["color_bar_position"].set(cbpos, 4);
  render["screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";
}
} // namespace

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::defineActions_insitu(glsol &sol) {

    conduit::Node &add_act = actions.append();    
    add_act["action"] = "add_scenes";
    conduit::Node &scenes = add_act["scenes"];

    conduit::Node &add_act2 = actions.append();    
    add_act2["action"] = "add_scenes";
    conduit::Node &scenes2 = add_act2["scenes"];

    conduit::Node &add_act3 = actions.append();
    add_act3["action"] = "add_pipelines";
    conduit::Node &pipelines = add_act3["pipelines"];

    conduit::Node &add_act4 = actions.append();
    add_act4["action"] = "add_pipelines";
    conduit::Node &pipelines2 = add_act4["pipelines"];
    
    
    // background color vecctor
    double bg_colvec[3] = {0.8, 0.8, 0.8};
    // foreground color vector
    double fg_colvec[3] = {0., 0., 0.};

    // color bar position for slice
    double colorbar_position[4] = {sol.config.CBxMin, sol.config.CBxMax, sol.config.CByMin, sol.config.CByMax};       
    
    /* >>>>>>>>>>>>>> pipleline gapA  clip <<<<<<<<<<<<< */
    if (sol.config.pario_compute_gapA == 1)
      {
    if (sol.config.do_gapA_clip == 1)
      {
       pipelines["pl1/f1/type"] = "clip";
       conduit::Node &clip_params = pipelines["pl1/f1/params"];

       clip_params["topology"] = "topo";
       clip_params["plane/point/x"] = sol.config.gapA_clip1_point_x;
       clip_params["plane/point/y"] = sol.config.gapA_clip1_point_y;
       clip_params["plane/point/z"] = sol.config.gapA_clip1_point_z;
       clip_params["plane/normal/x"] = sol.config.gapA_clip1_norm_x;
       clip_params["plane/normal/y"] = sol.config.gapA_clip1_norm_y;
       clip_params["plane/normal/z"] = sol.config.gapA_clip1_norm_z;

       scenes["s2/plots/p1/type"] = "pseudocolor";
       scenes["s2/plots/p1/pipeline"] = "pl1";
       scenes["s2/plots/p1/field"] = "gapA";
       scenes["s2/plots/p1/color_table/name"] = "Jet";

       scenes["s2/plots/p1/min_value"] = 0.0;
       scenes["s2/plots/p1/max_value"]
	 = (//-----------------------------------------
	    // if T-conf is hotBlob, use Ttdb0 for gap
	    //-----------------------------------------
	    sol.config.initialConditionT == 2
	   )
	   ? matep.gap_B_td(sol.config.Inip,
			    sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))
	                   * (1. + sol.config.clamp_bias_gapMax)
	   //---------------------------------------------------
	   // T-conf is homogeneous quench, use Ttd_Qend for gap
	   //---------------------------------------------------	 
	   : matep.gap_B_td(sol.config.Inip,
			    sol.config.Ttd_Qend * matep.Tcp_mK(sol.config.Inip))
	                   * (1. + sol.config.clamp_bias_gapMax);

       scenes["s2/renders/r1/color_bar_position"].set(colorbar_position,4);              
       
       scenes["s2/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s2/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s2/renders/r1/image_prefix"] = "insitu/gapA-clip1/gapA-clip1_t-%09d";
       scenes["s2/renders/r1/camera/azimuth"] = sol.config.camera2_azi;
       scenes["s2/renders/r1/camera/elevation"] = sol.config.camera2_ele;
      }
    /* >>>>>>>>>>>>>> pipleline gapA clip 2 <<<<<<<<<<<<< */

    // if (sol.config.do_gapA_clip == 1)
    //   {
    //    pipelines["pl2/f1/type"] = "clip";
    //    conduit::Node &clip_params2 = pipelines["pl2/f1/params"];

    //    clip_params2["topology"] = "topo";
    //    clip_params2["plane/point/x"] = sol.config.gapA_clip2_point_x;
    //    clip_params2["plane/point/y"] = sol.config.gapA_clip2_point_y;
    //    clip_params2["plane/point/z"] = sol.config.gapA_clip2_point_z;
    //    clip_params2["plane/normal/x"] = sol.config.gapA_clip2_norm_x;
    //    clip_params2["plane/normal/y"] = sol.config.gapA_clip2_norm_y;
    //    clip_params2["plane/normal/z"] = sol.config.gapA_clip2_norm_z;

    //    scenes["s3/plots/p1/type"] = "pseudocolor";
    //    scenes["s3/plots/p1/pipeline"] = "pl2";
    //    scenes["s3/plots/p1/field"] = "gapA";

    //    scenes["s3/plots/p1/min_value"]
    // 	 = matep.gap_A_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMin);

    //    scenes["s3/plots/p1/max_value"]
    // 	 = matep.gap_B_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax);

    //    scenes["s3/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
    //    scenes["s3/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);        
    //    scenes["s3/renders/r1/image_prefix"] = "insitu/gapA-clip2/gapA-clip2_t-%09d";
    //    scenes["s3/renders/r1/camera/azimuth"] = sol.config.camera2_azi/*35.0*/;
    //    scenes["s3/renders/r1/camera/elevation"] = sol.config.camera2_ele/*30.0*/;
    //   }

    /* >>>>>>>>>>>>>> pipleline gapA slice1 <<<<<<<<<<<<< */
    
    if (sol.config.do_gapA_slice == 1)
      {
       pipelines["pl3/f1/type"] = "exaslice";
       conduit::Node &slice_params1 = pipelines["pl3/f1/params"];

       slice_params1["point/x"] = sol.config.gapA_slice1_point_x;
       slice_params1["point/y"] = sol.config.gapA_slice1_point_y;
       slice_params1["point/z"] = sol.config.gapA_slice1_point_z;
       slice_params1["normal/x"] = sol.config.gapA_slice1_norm_x;
       slice_params1["normal/y"] = sol.config.gapA_slice1_norm_y;
       slice_params1["normal/z"] = sol.config.gapA_slice1_norm_z;

       scenes["s4/plots/p1/type"] = "pseudocolor";
       scenes["s4/plots/p1/pipeline"] = "pl3";
       scenes["s4/plots/p1/field"] = "gapA";
       scenes["s4/plots/p1/color_table/name"] = "Jet";       

       scenes["s4/plots/p1/min_value"] = 0.0;
       scenes["s4/plots/p1/max_value"]
	 = (//-----------------------------------------
	    // if T-conf is hotBlob, use Ttdb0 for gap
	    //-----------------------------------------
	    sol.config.initialConditionT == 2
	   )
	   ? matep.gap_B_td(sol.config.Inip,
			    sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))
	                   * (1. + sol.config.clamp_bias_gapMax)
	   //---------------------------------------------------
	   // T-conf is homogeneous quench, use Ttd_Qend for gap
	   //---------------------------------------------------	 
	   : matep.gap_B_td(sol.config.Inip,
			    sol.config.Ttd_Qend * matep.Tcp_mK(sol.config.Inip))
	                   * (1. + sol.config.clamp_bias_gapMax);

       setSliceAnnotations(scenes["s4/renders/r1"], colorbar_position, sol);       
       
       scenes["s4/renders/r1/image_width"]  = sol.config.image_width1;
       scenes["s4/renders/r1/image_height"] = sol.config.image_height1;
       scenes["s4/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s4/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);
       scenes["s4/renders/r1/image_prefix"] = "insitu/gapA-slice1/gapA-slice1_t-%09d";
       scenes["s4/renders/r1/camera/azimuth"] = sol.config.camera1_azi;
       scenes["s4/renders/r1/camera/elevation"] = sol.config.camera1_ele;
       scenes["s4/renders/r1/camera/zoom"] = sol.config.zoom1;                     
      }
    /* >>>>>>>>>>>>>> pipleline gapA slice2 <<<<<<<<<<<<< */
    
    // if (sol.config.do_gapA_slice == 1)
    //   {
    //    pipelines2["pl4/f1/type"] = "slice";
    //    conduit::Node &slice_params2 = pipelines2["pl4/f1/params"];

    //    slice_params2["point/x"] = sol.config.gapA_slice2_point_x;
    //    slice_params2["point/y"] = sol.config.gapA_slice2_point_y;
    //    slice_params2["point/z"] = sol.config.gapA_slice2_point_z;
    //    slice_params2["normal/x"] = sol.config.gapA_slice2_norm_x;
    //    slice_params2["normal/y"] = sol.config.gapA_slice2_norm_y;
    //    slice_params2["normal/z"] = sol.config.gapA_slice2_norm_z;

    //    scenes2["s5/plots/p1/type"] = "pseudocolor";
    //    scenes2["s5/plots/p1/pipeline"] = "pl4";
    //    scenes2["s5/plots/p1/field"] = "gapA";
    //    // scenes["s5/image_prefix"] = "gapA-slice2_t-%09d";       

    //    scenes2["s5/plots/p1/min_value"]
    // 	 = 0.0;
    
    //    scenes2["s5/plots/p1/max_value"]
    // 	 = matep.gap_B_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax);

    //    scenes2["s5/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
    //    scenes2["s5/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
    //    scenes2["s5/renders/r1/image_prefix"] = "gapA-slice2_t-%09d";
    //    // scenes["s5/renders/r1/camera/azimuth"] = 45. /*sol.config.camera1_azi*/;
    //    // scenes["s5/renders/r1/camera/elevation"] = 0./*sol.config.camera1_ele*/;
    //   }

    /* >>>>>>>>>>> pipleline isosurfece <<<<<<<<<<<<< */

    if (sol.config.do_gapA_isosurface == 1)
      {
       pipelines["pl5/f1/type"] = "contour";

       conduit::Node &contour_params = pipelines["pl5/f1/params"];
       contour_params["field"] = "gapA";

       //gapA_iso_list_size has to be constexpr in order to fix VLAs warning from Clang
       //defalt set its size to be 2, this should be carefully obeyed in configuration file.       
       const unsigned int iso_list_size = sol.config.iso_values_vector.size();

       //double iso_vals[Temperature_iso_list_size];
       //fix iso_vals array length; fix Clang 19 VLAs warning       
       double iso_vals[2];
       for (unsigned int i = 0; i<iso_list_size; ++i) {iso_vals[i] = sol.config.iso_values_vector[i];}

       contour_params["iso_values"].set(iso_vals, iso_list_size);

       scenes["s6/plots/p1/type"] = "pseudocolor";
       scenes["s6/plots/p1/pipeline"] = "pl5";
       scenes["s6/plots/p1/field"] = "gapA";
       scenes["s6/renders/r1/image_prefix"] = "insitu/gapA-iso/gapA-iso_t-%09d";

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes["s6/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes["s6/renders/r1/render_bg"] = "true";

       scenes["s6/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s6/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes["s6/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes["s6/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }
    
      } // compute gapA block

    
    if (sol.config.pario_compute_feDensity == 1)
      {
    
    /* >>>>>>>>>>>>>> pipleline clip camp <<<<<<<<<<<<< */

    if (sol.config.do_fed_clip == 1)
      {
       pipelines["pl6/f1/type"] = "clip";
       conduit::Node &clip_params3 = pipelines["pl6/f1/params"];

       clip_params3["topology"] = "topo";
       clip_params3["plane/point/x"] = sol.config.fed_clip_point_x;
       clip_params3["plane/point/y"] = sol.config.fed_clip_point_y;
       clip_params3["plane/point/z"] = sol.config.fed_clip_point_z;
       clip_params3["plane/normal/x"] = sol.config.fed_clip_norm_x;
       clip_params3["plane/normal/y"] = sol.config.fed_clip_norm_y;
       clip_params3["plane/normal/z"] = sol.config.fed_clip_norm_z;

       scenes["s7/plots/p1/type"] = "pseudocolor";
       scenes["s7/plots/p1/pipeline"] = "pl6";
       scenes["s7/plots/p1/field"] = "feDensity";

       scenes["s7/plots/p1/min_value"] = matep.f_B_td(sol.config.Inip, sol.config.Ttd_Qend * matep.Tcp_mK(sol.config.Inip)) * (1. + sol.config.clamp_bias_fed_Min);

       scenes["s7/plots/p1/max_value"] = 1.0 * (1. + sol.config.clamp_bias_fed_Max);    

       scenes["s7/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s7/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);        
       scenes["s7/renders/r1/image_prefix"] = "insitu/feDensity-clip_Camp/feDensity-clip_Camp-t-%09d";
       scenes["s7/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes["s7/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
      }

      } // compute feDensity block

    if (sol.config.pario_Temperature_pStream == 1)
      {

    /* >>>>>>>>>>>>>> pipleline Temperaure clip no camp <<<<<<<<<<<<< */

    if (sol.config.do_Temperature_clip == 1)
      {
       pipelines["pl7/f1/type"] = "clip";
       conduit::Node &clip_params4 = pipelines["pl7/f1/params"];

       clip_params4["topology"] = "topo";
       clip_params4["plane/point/x"] = sol.config.Temperature_clip_point_x;
       clip_params4["plane/point/y"] = sol.config.Temperature_clip_point_y;
       clip_params4["plane/point/z"] = sol.config.Temperature_clip_point_z;
       clip_params4["plane/normal/x"] = sol.config.Temperature_clip_norm_x;
       clip_params4["plane/normal/y"] = sol.config.Temperature_clip_norm_y;
       clip_params4["plane/normal/z"] = sol.config.Temperature_clip_norm_z;

       scenes["s8/plots/p1/type"] = "pseudocolor";
       scenes["s8/plots/p1/pipeline"] = "pl7";
       scenes["s8/plots/p1/field"] = "Temperature";

       scenes["s8/plots/p1/min_value"] = sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip);

       // this upper limit needs to be calculated carefully
       scenes["s8/plots/p1/max_value"] = sol.config.Temperature_clamp * matep.Tcp_mK(sol.config.Inip);

       scenes["s8/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s8/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);        
       scenes["s8/renders/r1/image_prefix"] = "insitu/Temeperature-clip/Temeperature-clip-t-%09d";
       scenes["s8/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes["s8/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
      }

    /* >>>>>>>>>>>>>> pipleline Temperature slice <<<<<<<<<<<<< */
    
    if (sol.config.do_Temperature_slice == 1)
      {
       pipelines2["pl4/f1/type"] = "slice";
       conduit::Node &slice_params2 = pipelines2["pl4/f1/params"];

       slice_params2["point/x"] = sol.config.Temperature_slice_point_x;
       slice_params2["point/y"] = sol.config.Temperature_slice_point_y;
       slice_params2["point/z"] = sol.config.Temperature_slice_point_z;
       slice_params2["normal/x"] = sol.config.Temperature_slice_norm_x;
       slice_params2["normal/y"] = sol.config.Temperature_slice_norm_y;
       slice_params2["normal/z"] = sol.config.Temperature_slice_norm_z;

       scenes2["s5/plots/p1/type"] = "pseudocolor";
       scenes2["s5/plots/p1/pipeline"] = "pl4";
       scenes2["s5/plots/p1/field"] = "Temperature";

       scenes2["s5/plots/p1/min_value"]
	 = sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip);
    
       scenes2["s5/plots/p1/max_value"]
	 = matep.Tcp_mK(sol.config.Inip);

       setSliceAnnotations(scenes2["s5/renders/r1"], colorbar_position, sol);
       scenes2["s5/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s5/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes2["s5/renders/r1/image_prefix"] = "insitu/Temperature-slice/Temperature-slice_t-%09d";
      }

    
    /* >>>>>>>>>>> Tmeperature isosurfece <<<<<<<<<<<<< */

    if (sol.config.do_Temperature_isosurface == 1)
      {
       pipelines["pl8/f1/type"] = "contour";

       conduit::Node &contour_params = pipelines["pl8/f1/params"];
       contour_params["field"] = "Temperature";

       //Temperature_iso_list_size has to be constexpr in order to fix VLAs warning from Clang
       //defalt set its size to be 2, this should be carefully obeyed in configuration file.
       const unsigned int Temperature_iso_list_size = sol.config.Temperature_iso_values_vector.size();
       //double iso_vals[Temperature_iso_list_size];
       //fix iso_vals array length; fix Clang 19 VLAs warning
       double iso_vals[2];
       for (unsigned int i = 0; i<Temperature_iso_list_size; ++i) {iso_vals[i] = sol.config.Temperature_iso_values_vector[i];}

       contour_params["iso_values"].set(iso_vals, Temperature_iso_list_size);

       scenes["s9/plots/p1/type"] = "pseudocolor";
       scenes["s9/plots/p1/pipeline"] = "pl8";
       scenes["s9/plots/p1/field"] = "Temperature";
       scenes["s9/renders/r1/image_prefix"] = "insitu/Temperature-iso/Temperature-iso_t-%09d";

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes["s9/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes["s9/renders/r1/render_bg"] = "true";

       scenes["s9/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s9/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes["s9/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes["s9/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }

      } // Temperature pStream

    if (sol.config.pario_compute_phaseMarker == 1)
      {
    /* >>>>>>>>>>>>>> pipleline phaseMarker slice 1 <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_slice1 == 1)
      {
       pipelines2["pl9/f1/type"] = "exaslice";
       conduit::Node &slice_params3 = pipelines2["pl9/f1/params"];

       slice_params3["point/x"] = sol.config.pMarker_slice1_point_x;
       slice_params3["point/y"] = sol.config.pMarker_slice1_point_y;
       slice_params3["point/z"] = sol.config.pMarker_slice1_point_z;
       slice_params3["normal/x"] = sol.config.pMarker_slice1_norm_x;
       slice_params3["normal/y"] = sol.config.pMarker_slice1_norm_y;
       slice_params3["normal/z"] = sol.config.pMarker_slice1_norm_z;
       addMarkerZeroSplit(pipelines2, "pl9");

       scenes["s10/plots/p1/type"] = "pseudocolor";
       scenes["s10/plots/p1/pipeline"] = "pl9";
       scenes["s10/plots/p1/field"] = "phaseMarker";
       //scenes["s10/plots/p1/color_table/name"] = "Jet";
       scenes["s10/plots/p1/color_table/discrete"] = "true";
       scenes["s10/plots/p1/min_value"] = 1.0f;    
       scenes["s10/plots/p1/max_value"] = 9.0f; 

       conduit::Node control_points;
       colorControlPointsGenerator(control_points);       
       scenes["s10/plots/p1/color_table/control_points"] = control_points;       
       addMarkerZeroPlot(scenes["s10/plots/p2"], "pl9");
       //double colorbar_position[4] = {-0.6, 0.6, 0.8, 0.9};       
       setSliceAnnotations(scenes["s10/renders/r1"], colorbar_position, sol);       
       
       scenes["s10/renders/r1/image_width"]  = sol.config.image_width1;
       scenes["s10/renders/r1/image_height"] = sol.config.image_height1;       
       scenes["s10/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s10/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);
       scenes["s10/renders/r1/image_prefix"] = "insitu/pMarker-slice1/pMarker-slice1_t-%09d";
       scenes["s10/renders/r1/camera/azimuth"] = sol.config.camera1_azi;
       scenes["s10/renders/r1/camera/elevation"] = sol.config.camera1_ele;
       scenes["s10/renders/r1/camera/zoom"] = sol.config.zoom1;              
      }

    /* >>>>>>>>>>>>>> pipleline phaseMarker slice 2 <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_slice2 == 1)
      {
       pipelines2["pl10/f1/type"] = "exaslice";
       conduit::Node &slice_params3 = pipelines2["pl10/f1/params"];

       slice_params3["point/x"] = sol.config.pMarker_slice2_point_x;
       slice_params3["point/y"] = sol.config.pMarker_slice2_point_y;
       slice_params3["point/z"] = sol.config.pMarker_slice2_point_z;
       slice_params3["normal/x"] = sol.config.pMarker_slice2_norm_x;
       slice_params3["normal/y"] = sol.config.pMarker_slice2_norm_y;
       slice_params3["normal/z"] = sol.config.pMarker_slice2_norm_z;
       addMarkerZeroSplit(pipelines2, "pl10");

       scenes["s11/plots/p1/type"] = "pseudocolor";
       scenes["s11/plots/p1/pipeline"] = "pl10";
       scenes["s11/plots/p1/field"] = "phaseMarker";
       //scenes["s11/plots/p1/color_table/name"] = "Jet";
       scenes["s11/plots/p1/color_table/discrete"] = "true";
       scenes["s11/plots/p1/min_value"] = 1.0f;
       scenes["s11/plots/p1/max_value"] = 9.0f; 

       conduit::Node control_points;
       colorControlPointsGenerator(control_points);       
       scenes["s11/plots/p1/color_table/control_points"] = control_points;       
       addMarkerZeroPlot(scenes["s11/plots/p2"], "pl10");
       //double colorbar_position[4] = {-0.6, 0.6, 0.8, 0.9};       
       setSliceAnnotations(scenes["s11/renders/r1"], colorbar_position, sol);       
       
       scenes["s11/renders/r1/image_width"]  = sol.config.image_width2;
       scenes["s11/renders/r1/image_height"] = sol.config.image_height2;       
       scenes["s11/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s11/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);
       scenes["s11/renders/r1/image_prefix"] = "insitu/pMarker-slice2/pMarker-slice2_t-%09d";
       scenes["s11/renders/r1/camera/azimuth"] = sol.config.camera2_azi;
       scenes["s11/renders/r1/camera/elevation"] = sol.config.camera2_ele;
       scenes["s11/renders/r1/camera/zoom"] = sol.config.zoom2;             
      }
    
    /* >>>>>>>>>>> phaseMarker isosurfece <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_isosurface == 1)
      {
       pipelines2["pl11/f1/type"] = "contour";

       conduit::Node &contour_params = pipelines2["pl11/f1/params"];
       contour_params["field"] = "phaseMarker";

       //phaseMarker_iso_list_size has to be constexpr in order to fix VLAs warning from Clang
       //defalt set its size to be 2, this should be carefully obeyed in configuration file.       
       const unsigned int phaseMarker_iso_list_size = sol.config.phaseMarker_iso_values_vector.size();
       //double iso_vals[Temperature_iso_list_size];
       //fix iso_vals array length; fix Clang 19 VLAs warning       
       double iso_vals[2];
       for (unsigned int i = 0; i<phaseMarker_iso_list_size; ++i) {iso_vals[i] = sol.config.phaseMarker_iso_values_vector[i];}

       contour_params["iso_values"].set(iso_vals, phaseMarker_iso_list_size);

       scenes2["s12/plots/p1/type"] = "pseudocolor";
       scenes2["s12/plots/p1/pipeline"] = "pl11";
       scenes2["s12/plots/p1/field"] = "phaseMarker";
       scenes2["s12/plots/p1/color_table/name"] = "Jet";       
       scenes2["s12/renders/r1/image_prefix"] = "insitu/pMarker-iso/pMarker-iso_t-%09d";
       scenes2["s12/renders/r1/screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";       

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes2["s12/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes2["s12/renders/r1/render_bg"] = "true";
       scenes2["s12/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s12/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes2["s12/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes2["s12/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }

    /* >>>>>>>>>>> phaseMarker pMarker < 2 fieldclip <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_fieldclip == 1)
      {
       pipelines2["pl12/f1/type"] = "clip_with_field";

       conduit::Node &clip_params = pipelines2["pl12/f1/params"];
       clip_params["field"] = "phaseMarker";
       clip_params["invert"] = "true";
       clip_params["clip_value"] = 2.;

       scenes2["s13/plots/p1/type"] = "pseudocolor";
       scenes2["s13/plots/p1/pipeline"] = "pl12";
       scenes2["s13/plots/p1/field"] = "phaseMarker";
       scenes2["s13/plots/p1/color_table/name"] = "Default";
       //scenes2["s13/plots/p1/color_table/discrete"] = "true";
       scenes2["s13/renders/r1/image_prefix"] = "insitu/pMarker-fieldclip/pMarker-fieldclip_t-%09d";
       scenes2["s13/renders/r1/screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";       

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes2["s13/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes2["s13/renders/r1/render_bg"] = "true";

       scenes2["s13/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s13/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes2["s13/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes2["s13/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }

    /* >>>>>>>>>>> phaseMarker B-isoVolume <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_fieldclip_Bphase == 1)
      {
       pipelines2["pl13/f1/type"] = "isovolume";

       conduit::Node &clip_params = pipelines2["pl13/f1/params"];
       clip_params["field"] = "phaseMarker";
       clip_params["min_value"] = 3.0;
       clip_params["max_value"] = 6.0;

       scenes2["s14/plots/p1/type"] = "pseudocolor";
       scenes2["s14/plots/p1/pipeline"] = "pl13";
       scenes2["s14/plots/p1/field"] = "phaseMarker";
       scenes2["s14/plots/p1/color_table/name"] = "Green";
       //scenes2["s13/plots/p1/color_table/discrete"] = "true";
       scenes2["s14/renders/r1/image_prefix"] = "insitu/pMarker-isoVolume-Bphase/pMarker-isoVolume-Bphase_t-%09d";
       scenes2["s14/renders/r1/screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";       

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes2["s14/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes2["s14/renders/r1/render_bg"] = "true";

       scenes2["s14/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s14/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes2["s14/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes2["s14/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }

    /* >>>>>>>>>>> phaseMarker A-phase fieldclip <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_fieldclip_Aphase == 1)
      {
       pipelines2["pl14/f1/type"] = "clip_with_field";

       conduit::Node &clip_params = pipelines2["pl14/f1/params"];
       clip_params["field"] = "phaseMarker";
       clip_params["invert"] = "false";
       clip_params["clip_value"] = 8.;

       scenes2["s15/plots/p1/type"] = "pseudocolor";
       scenes2["s15/plots/p1/pipeline"] = "pl14";
       scenes2["s15/plots/p1/field"] = "phaseMarker";
       scenes2["s15/plots/p1/color_table/name"] = "Cold and Hot";
       //scenes2["s15/plots/p1/color_table/discrete"] = "true";
       scenes2["s15/renders/r1/image_prefix"] = "insitu/pMarker-fieldclip-Aphase/pMarker-fieldclip-Aphase_t-%09d";
       scenes2["s15/renders/r1/screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes2["s15/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes2["s15/renders/r1/render_bg"] = "true";

       scenes2["s15/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s15/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes2["s15/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes2["s15/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
      }

      } // compute phaseMarker block

    /* >>>>>>>>>>>>>> pipleline D_B slice <<<<<<<<<<<<< */
    // One exaslice pipeline; one scene per colour table in DB_color_tables,
    // all with the same fixed range so the colour maps can be compared.
    // Images: insitu/DB-slice/DB-slice_<table>_t-<cycle>.png
    if ((sol.config.pario_compute_DB == 1) && (sol.config.do_DB_slice == 1))
      {
       pipelines2["plDB/f1/type"] = "exaslice";
       conduit::Node &slice_params = pipelines2["plDB/f1/params"];
       slice_params["point/x"] = sol.config.DB_slice_point_x;
       slice_params["point/y"] = sol.config.DB_slice_point_y;
       slice_params["point/z"] = sol.config.DB_slice_point_z;
       slice_params["normal/x"] = sol.config.DB_slice_norm_x;
       slice_params["normal/y"] = sol.config.DB_slice_norm_y;
       slice_params["normal/z"] = sol.config.DB_slice_norm_z;
       // Split the slice at masked sites (D_B = DB_MASKED < 0).
       // plDB: cells whose vertices are all valid (select cells with any
       //       masked vertex, then invert), so nothing interpolates to DB_MASKED.
       // plDBmask: the complement, drawn in the flat colour DB_mask_rgb.
       pipelines2["plDB/f2/type"] = "threshold";
       pipelines2["plDB/f2/params/field"] = "DB";
       pipelines2["plDB/f2/params/min_value"] = -1.0e30;
       pipelines2["plDB/f2/params/max_value"] = -0.5;
       pipelines2["plDB/f2/params/invert"] = "true";

       pipelines2["plDBmask/f1"] = pipelines2["plDB/f1"];
       pipelines2["plDBmask/f2/type"] = "threshold";
       pipelines2["plDBmask/f2/params/field"] = "DB";
       pipelines2["plDBmask/f2/params/min_value"] = -1.0e30;
       pipelines2["plDBmask/f2/params/max_value"] = -0.5;

       // one-colour table: the same RGB at both ends
       static double mask_r[2], mask_g[2], mask_b[2];
       static double mask_a[2] = {1.0, 1.0}, mask_pos[2] = {0.0, 1.0};
       for (int j = 0; j < 2; j++)
	 {
	  mask_r[j] = sol.config.DB_mask_rgb[0];
	  mask_g[j] = sol.config.DB_mask_rgb[1];
	  mask_b[j] = sol.config.DB_mask_rgb[2];
	 }

       for (size_t i = 0; i < sol.config.DB_color_tables.size(); i++)
	 {
	  const std::string &table = sol.config.DB_color_tables[i];
	  // file-name-safe version of the table name
	  std::string tag;
	  for (const char c : table) tag += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';

	  const std::string s = "sDB" + std::to_string(i);
	  scenes2[s + "/plots/p1/type"] = "pseudocolor";
	  scenes2[s + "/plots/p1/pipeline"] = "plDB";
	  scenes2[s + "/plots/p1/field"] = "DB";
	  scenes2[s + "/plots/p1/color_table/name"] = table;
	  scenes2[s + "/plots/p1/min_value"] = sol.config.DB_range_min;
	  scenes2[s + "/plots/p1/max_value"] = sol.config.DB_range_max;

	  scenes2[s + "/plots/p2/type"] = "pseudocolor";
	  scenes2[s + "/plots/p2/pipeline"] = "plDBmask";
	  scenes2[s + "/plots/p2/field"] = "DB";
	  scenes2[s + "/plots/p2/min_value"] = -2.0;   // avoid a zero-width range
	  scenes2[s + "/plots/p2/max_value"] = 0.0;
	  conduit::Node &ct = scenes2[s + "/plots/p2/color_table"];
	  ct["control_points/r"].set_external(mask_r, 2);
	  ct["control_points/g"].set_external(mask_g, 2);
	  ct["control_points/b"].set_external(mask_b, 2);
	  ct["control_points/a"].set_external(mask_a, 2);
	  ct["control_points/position"].set_external(mask_pos, 2);
	  ct["annotation"] = "false";   // colour bar for D_B only

	  setSliceAnnotations(scenes2[s + "/renders/r1"], colorbar_position, sol);
	  scenes2[s + "/renders/r1/image_width"]  = sol.config.image_width1;
	  scenes2[s + "/renders/r1/image_height"] = sol.config.image_height1;
	  scenes2[s + "/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
	  scenes2[s + "/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);
	  scenes2[s + "/renders/r1/image_prefix"] = "insitu/DB-slice/DB-slice_" + tag + "_t-%09d";
	  scenes2[s + "/renders/r1/camera/azimuth"] = sol.config.camera1_azi;
	  scenes2[s + "/renders/r1/camera/elevation"] = sol.config.camera1_ele;
	  scenes2[s + "/renders/r1/camera/zoom"] = sol.config.zoom1;
	 }
      } // D_B slice block
    
    /* >>>>>>>>>>>>>> pipleline U(1) 3 \phi slice <<<<<<<<<<<<< */

    if (sol.config.pario_compute_U1Phase == 1)
      {
    
    if (sol.config.do_U13phi_slice == 1)
      {
       pipelines2["pl15/f1/type"] = "exaslice";
       conduit::Node &slice_params4 = pipelines2["pl15/f1/params"];

       slice_params4["point/x"] = sol.config.U13phi_slice_point_x;
       slice_params4["point/y"] = sol.config.U13phi_slice_point_y;
       slice_params4["point/z"] = sol.config.U13phi_slice_point_z;
       slice_params4["normal/x"] = sol.config.U13phi_slice_norm_x;
       slice_params4["normal/y"] = sol.config.U13phi_slice_norm_y;
       slice_params4["normal/z"] = sol.config.U13phi_slice_norm_z;

       /* ----------------------------------------------------------
        * custom hsv style Cyclic color map, 
        * RGBA control points are defined in hsvCyclicSchemeRGBA.hpp
        * ----------------------------------------------------------
        */
       conduit::Node hsv_Cyclic;
       hsv_Cyclic["r"].set_external(hsv_r);
       hsv_Cyclic["g"].set_external(hsv_g);
       hsv_Cyclic["b"].set_external(hsv_b);
       hsv_Cyclic["a"].set_external(hsv_a);
       hsv_Cyclic["position"].set_external(hsv_value_point_position);
       
       scenes["s16/plots/p1/type"] = "pseudocolor";
       scenes["s16/plots/p1/pipeline"] = "pl15";
       scenes["s16/plots/p1/field"] = "U1_3phi";
       scenes["s16/plots/p1/min_value"] = -3.141592653589793;
       scenes["s16/plots/p1/max_value"] = 3.141592653589793;
       //scenes["s16/plots/p1/color_table/name"] = "Blue to Orange";
       scenes["s16/plots/p1/color_table/control_points"] = hsv_Cyclic;
       // scenes["s15/plots/p1/color_table/discrete"] = "true";

       setSliceAnnotations(scenes["s16/renders/r1"], colorbar_position, sol);
       scenes["s16/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s16/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s16/renders/r1/image_prefix"] = "insitu/U1_3phi-Slice/U1_3phi-slice_t-%09d";
       scenes["s16/renders/r1/camera/zoom"] = sol.config.zoom1;   // same as the other z-normal slices
      }

      } // compute U1Phase block


    if (sol.config.pario_compute_lVector == 1)
      {
    /* >>>>>>>>>>>>>> pipleline l^2 slice <<<<<<<<<<<<< */
    
    if (sol.config.do_l_sq_slice == 1)
      {
       pipelines2["pl16/f1/type"] = "exaslice";
       conduit::Node &slice_params5 = pipelines2["pl16/f1/params"];

       slice_params5["point/x"] = sol.config.l_sq_slice_point_x;
       slice_params5["point/y"] = sol.config.l_sq_slice_point_y;
       slice_params5["point/z"] = sol.config.l_sq_slice_point_z;
       slice_params5["normal/x"] = sol.config.l_sq_slice_norm_x;
       slice_params5["normal/y"] = sol.config.l_sq_slice_norm_y;
       slice_params5["normal/z"] = sol.config.l_sq_slice_norm_z;

       scenes["s17/plots/p1/type"] = "pseudocolor";
       scenes["s17/plots/p1/pipeline"] = "pl16";
       scenes["s17/plots/p1/field"] = "l_Sq";
       scenes["s17/plots/p1/color_table/name"] = "Cool to Warm Extended";
       // scenes["s16/plots/p1/color_table/discrete"] = "true";

       scenes["s17/plots/p1/min_value"]
	 = 0.f; /*1.f - sol.config.lVec_SqlTol;*/
    
       scenes["s17/plots/p1/max_value"]
	 = 1.f + sol.config.lVec_SqlTol; 

       setSliceAnnotations(scenes["s17/renders/r1"], colorbar_position, sol);
       scenes["s17/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s17/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s17/renders/r1/image_prefix"] = "insitu/lVec_norm2-Slice/lVec_norm2-slice_t-%09d";
       scenes["s17/renders/r1/camera/zoom"] = sol.config.zoom1;   // same as the other z-normal slices
      }

      } // compute lVector block

    
} // defineActions() call end here

