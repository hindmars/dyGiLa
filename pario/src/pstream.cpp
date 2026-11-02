#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
#include <assert.h>
#include <stdlib.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::pstream(glsol &sol, unsigned int &stat_counter) {

    /*-------------------    sqrt(Tr[A.A^+) ----------------------*/
    if (sol.config.pario_compute_gapA == 1)
      {
       //gapA[ALL] = real(sqrt((sol.AwT[X]*(sol.AwT[X].dagger())).trace()));
       gapA[ALL] = sqrt(((sol.AwT[X]*(sol.AwT[X].dagger())).trace()).real());
       gapA.copy_local_data_with_halo(gapAContainer);
      }

    /*--------------------     feDensity      --------------------*/
    if (sol.config.pario_compute_feDensity == 1) {
      feDensity[ALL] = 0;
      onsites(ALL) {
        Complex<double> a(0),b1(0),b2(0),b3(0),b4(0),b5(0);
        Complex<double> k1(0), k2(0), k3(0);
      
        a = sol.config.alpha * (sol.AwT[X]*sol.AwT[X].dagger()).trace();
        b1 = sol.config.beta1 * ((sol.AwT[X]*sol.AwT[X].transpose()).trace()).squarenorm();
        b2 = sol.config.beta2 * ((sol.AwT[X]*sol.AwT[X].dagger()).trace()*(sol.AwT[X]*sol.AwT[X].dagger()).trace());
        b3 = sol.config.beta3 * ((sol.AwT[X]*sol.AwT[X].transpose()*sol.AwT[X].conj()*sol.AwT[X].dagger()).trace());
        b4 = sol.config.beta4 * ((sol.AwT[X]*sol.AwT[X].dagger()*sol.AwT[X]*sol.AwT[X].dagger()).trace());
        b5 = sol.config.beta5 * ((sol.AwT[X]*sol.AwT[X].dagger()*sol.AwT[X].conj()*sol.AwT[X].transpose()).trace());


	/* -----------------------------------------------------------------------------
         * Warining & TODO: the 2h of differetial devrivative only works for dx=0.5, 
         * this should be treated as bug becuase you gonna forget it when dx is changing.
         * -----------------------------------------------------------------------------
         */
        foralldir(j) foralldir (k) foralldir(al){
	  k1 += (sol.AwT[X + k].column(j) - sol.AwT[X - k].column(j)).e(al)
	        * (sol.AwT[X + k].conj().column(j) - sol.AwT[X - k].conj().column(j)).e(al)/(4.0*sol.config.dx*sol.config.dx);
	  k2 += (sol.AwT[X + j].column(j) - sol.AwT[X - j].column(j)).e(al)
	        * (sol.AwT[X + k].conj().column(k) - sol.AwT[X - k].conj().column(k)).e(al)/(4.0*sol.config.dx*sol.config.dx);
	  k3 += (sol.AwT[X + k].column(j) - sol.AwT[X - k].column(j)).e(al)
	        * (sol.AwT[X + j].conj().column(k) - sol.AwT[X - j].conj().column(k)).e(al)/(4.0*sol.config.dx*sol.config.dx);
        }
       feDensity[X] = real(k1 + k2 + k3 + a + b1 + b2 + b3 + b4 + b5);
      } //onsite(All) end here

      feDensity.copy_local_data_with_halo(feDensityContainer);

    } // sol.config.pario_compute_feDensity == 1 block

    /*-------------------- Temperature field --------------------*/
    if (sol.config.pario_Temperature_pStream == 1) { sol.T.copy_local_data_with_halo(Temperature); }
    
    /*-------------------- phaseMarker field --------------------*/
    if (sol.config.pario_compute_phaseMarker == 1) { sol.phaseMarker.copy_local_data_with_halo(phaseMarker); }

    /*-------------------- D_B field --------------------*/
    if (sol.config.pario_compute_DB == 1) { sol.DB.copy_local_data_with_halo(DBContainer); }

    /*------------------- U(1) phase of AwT ---------------------*/
    if (sol.config.pario_compute_U1Phase == 1) {
        U1PhaseStreaming(sol);
        U1_3phi.copy_local_data_with_halo(U1_3phiContainer);            
    } //{ U1PhaseStreaming(sol); }

    /*--------------------- l-vector of AwT ---------------------*/
    if (sol.config.pario_compute_lVector == 1) {
      lVectorStreaming(sol);
      lsq.copy_local_data_with_halo(lsqContainer);
      if (sol.config.hdf5_lVector_output == 1)
        {
          l_1.copy_local_data_with_halo(l1Container);
          l_2.copy_local_data_with_halo(l2Container);
          l_3.copy_local_data_with_halo(l3Container);
        } // only harvest l_i if want to output them      
    }

    /*--------------- Gradient Phi -vector of AwT ----------------*/
    if (sol.config.pario_compute_GPhiVector == 1) {
      GradientPhiVectorStreaming(sol);
      if (
	  sol.config.hdf5_GradPhiVector_output == 1
	  || sol.config.hdf5_GradPhiVector_exaslice_output == 1
	 )
        {
          GPhi_1.copy_local_data_with_halo(GPhi1Container);
          GPhi_2.copy_local_data_with_halo(GPhi2Container);
          GPhi_3.copy_local_data_with_halo(GPhi3Container);
        } // only harvest l_i if want to output them      
    }
    
    /*---------------- mass current components ------------------*/
    if (sol.config.hdf5_mass_current_output == 1){
      Field<Vector<3,double>> jmX;
      // Field<Complex<real_t>>  phaseExp;

      /*onsites(ALL){
        foralldir(i) foralldir(j) foralldir(al){
	  jmX[X].e(i) = ((sol.AwT[X].conj().column(j)).e(al) * (A[X+i].column(j) - A[X-i].column(j)).e(al)/(2.*config.dx)
	                + (sol.AwT[X].conj().column(j)).e(al) * (A[X+j].column(i) - A[X-j].column(i)).e(al)/(2.*config.dx)
			+ (sol.AwT[X].conj().column(i)).e(al) * (A[X+j].column(j) - A[X-j].column(j)).e(al)/(2.*config.dx)).imag();

        } // foralldir() calls end here

	} //onesite(ALL) call ends here*/

      onsites(ALL) {
	jmX[X] = 0;
	foralldir(i) foralldir(j) foralldir(al) {
	  jmX[X].e(i) += (sol.AwT[X].e(al,j).conj() *(sol.AwT[X+i].e(al,j) - sol.AwT[X-i].e(al,j))
			  + sol.AwT[X].e(al,j).conj() * (sol.AwT[X+j].e(al,i) - sol.AwT[X-j].e(al,i))
			  + sol.AwT[X].e(al,i).conj() * (sol.AwT[X+j].e(al,j) - sol.AwT[X-j].e(al,j))).imag();
	} // foralldir end here, outermost foralldir slowest, inner run earier
	jmX[X] /= 2.*sol.config.dx;
      } // onsites(ALL) end here

      jm1[ALL] = jmX[X].e(0);
      jm2[ALL] = jmX[X].e(1);
      jm3[ALL] = jmX[X].e(2);

      /* >>>>>>>> Modulus, phase angle, phaseExp   <<<<<<< */
      // phaseExp[ALL]        = ((sol.AwT[X].transpose()) * sol.AwT[X]).trace();
      // phaseExp2Re[ALL]    = phaseExp[X].real();
      // phaseExp2Im[ALL]    = phaseExp[X].imag();      
      
      // phaseExpModulus[ALL] = phaseExp[X].abs();
      //phaseExpAngle[ALL]   = phaseExp[X].arg()/2.;
      // phaseExpAngle[ALL]   = std::atan2(phaseExp[X].imag(), phaseExp[X].real())/2.;
      /* >>>>>>>> phase angle and modulus end here   <<<<<<*/

      jm1.copy_local_data_with_halo(jm1Container);
      jm2.copy_local_data_with_halo(jm2Container);
      jm3.copy_local_data_with_halo(jm3Container);

      // phaseExpModulus.copy_local_data_with_halo(phaseExpModulusContainer);
      // phaseExpAngle.copy_local_data_with_halo(phaseExpAngleContainer);
      // phaseExp2Re.copy_local_data_with_halo(phaseExp2ReContainer);
      // phaseExp2Im.copy_local_data_with_halo(phaseExp2ImContainer);            
    }

    // /*------------------ spin current components ------------------*/
    if (
	(sol.config.hdf5_spin_current_output == 1)
	|| (sol.config.hdf5_spinCurr_exaslice_output == 1)
       ){
      Field<Matrix<3,3,double>> jsX; // column is alpha for spin, row is i for spatial

      onsites(ALL) {
	matep::Matep MPonsites;
	jsX[X] = 0;
	foralldir(al) foralldir(i) foralldir(be) foralldir(ga) foralldir(j) {
          jsX[X].e(i,al) += -MPonsites.epsilon(al,be,ga)
	                     * (sol.AwT[X].e(be,i).conj() * (sol.AwT[X+j].e(ga,j) - sol.AwT[X-j].e(ga,j))
			        + sol.AwT[X].e(be,j).conj() * (sol.AwT[X+i].e(ga,j) - sol.AwT[X-i].e(ga,j))
			        + sol.AwT[X].e(be,j).conj() * (sol.AwT[X+j].e(ga,i) - sol.AwT[X-j].e(ga,i))).real();
	  
	} // foralldir end here, outermost foralldir slowest, inner run earier
	jsX[X] /= 2.*sol.config.dx;
      } // onsites(ALL) end here

      js11[ALL] = jsX[X].e(0,0);
      js21[ALL] = jsX[X].e(1,0);
      js31[ALL] = jsX[X].e(2,0);
      js12[ALL] = jsX[X].e(0,1);
      js22[ALL] = jsX[X].e(1,1);
      js32[ALL] = jsX[X].e(2,1);
      js13[ALL] = jsX[X].e(0,2);
      js23[ALL] = jsX[X].e(1,2);
      js33[ALL] = jsX[X].e(2,2);
      
      js11.copy_local_data_with_halo(js11Container);
      js21.copy_local_data_with_halo(js21Container);
      js31.copy_local_data_with_halo(js31Container);
      
      js12.copy_local_data_with_halo(js12Container);
      js22.copy_local_data_with_halo(js22Container);
      js32.copy_local_data_with_halo(js32Container);
      
      js13.copy_local_data_with_halo(js13Container);
      js23.copy_local_data_with_halo(js23Container);
      js33.copy_local_data_with_halo(js33Container);

    }

    /*----------------     A matrix elements ---------------------*/
    if (
	(sol.config.hdf5_A_matrix_output == 1)
	|| (sol.config.hdf5_AMatrix_exaslice_output == 1)
       ){
     u11[ALL] = sol.AwT[X].e(0,0).re; v11[ALL] = sol.AwT[X].e(0,0).im;
     u12[ALL] = sol.AwT[X].e(0,1).re; v12[ALL] = sol.AwT[X].e(0,1).im;
     u13[ALL] = sol.AwT[X].e(0,2).re; v13[ALL] = sol.AwT[X].e(0,2).im;
     u21[ALL] = sol.AwT[X].e(1,0).re; v21[ALL] = sol.AwT[X].e(1,0).im;
     u22[ALL] = sol.AwT[X].e(1,1).re; v22[ALL] = sol.AwT[X].e(1,1).im;
     u23[ALL] = sol.AwT[X].e(1,2).re; v23[ALL] = sol.AwT[X].e(1,2).im;
     u31[ALL] = sol.AwT[X].e(2,0).re; v31[ALL] = sol.AwT[X].e(2,0).im;
     u32[ALL] = sol.AwT[X].e(2,1).re; v32[ALL] = sol.AwT[X].e(2,1).im;
     u33[ALL] = sol.AwT[X].e(2,2).re; v33[ALL] = sol.AwT[X].e(2,2).im;

     u11.copy_local_data_with_halo(u11Container); v11.copy_local_data_with_halo(v11Container);
     u12.copy_local_data_with_halo(u12Container); v12.copy_local_data_with_halo(v12Container);
     u13.copy_local_data_with_halo(u13Container); v13.copy_local_data_with_halo(v13Container);
     u21.copy_local_data_with_halo(u21Container); v21.copy_local_data_with_halo(v21Container);
     u22.copy_local_data_with_halo(u22Container); v22.copy_local_data_with_halo(v22Container);
     u23.copy_local_data_with_halo(u23Container); v23.copy_local_data_with_halo(v23Container);
     u31.copy_local_data_with_halo(u31Container); v31.copy_local_data_with_halo(v31Container);
     u32.copy_local_data_with_halo(u32Container); v32.copy_local_data_with_halo(v32Container);
     u33.copy_local_data_with_halo(u33Container); v33.copy_local_data_with_halo(v33Container);        
    }
    
    /*------------------------------------------------------------*/

    // update mesh[state/cycle]
    mesh["state/time"].set_external(&sol.t);
    mesh["state/cycle"] = stat_counter;

    // re-publish mesh with new cycle
    pio.publish(mesh);
    
    pio.execute(actions);

    system(R"(find . -name "*.root" -delete 2>/dev/null)");
} // pstream() end here

