#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::phaseMarking() {

  //Field<phi_t> AwT = A;
  // if (config.useGaussianLP_filter == 1)
  //   { GaussianLPfilter_matrix(AwT);} // shortwave noise removed A filed preperation
  
  onsites (ALL) {

    real_t R1{0},R2{0},R3{0},R4{0},R5{0};
    
    phi_t Ar = (config.useGaussianLP_filter == 1)
               ? AwT[X]/sqrt((AwT[X]*AwT[X].dagger()).trace())
               : A[X]/sqrt((A[X]*A[X].dagger()).trace());

    R1 = ((Ar*Ar.transpose()).trace()).squarenorm();
    R2 = real(((Ar*Ar.dagger()).trace()*(Ar*Ar.dagger()).trace()));
    R3 = real(((Ar*Ar.transpose()*Ar.conj()*Ar.dagger()).trace()));
    R4 = real(((Ar*Ar.dagger()*Ar*Ar.dagger()).trace()));
    R5 = real(((Ar*Ar.dagger()*Ar.conj()*Ar.transpose()).trace()));

    // The R's come from the normalised Ar, so they are undefined at A = 0 and
    // just noise where A is small: sites with a gap below invariants_gapA_min
    // are treated as normal (phaseMarker 0, D_B masked).
    const real_t gapSq = (config.useGaussianLP_filter == 1)
                         ? real((AwT[X]*AwT[X].dagger()).trace())
                         : real((A[X]*A[X].dagger()).trace());
    const bool smallGap = !(gapSq > config.invariants_gapA_min*config.invariants_gapA_min);

    /* >>>>>>>>> D_B: distance from the B phase in invariant space <<<<<<<<  */
    // B phase invariants: R^B = (1, 1, 1/3, 1/3, 1/3).
    if (config.pario_compute_DB == 1)
      {
       const real_t third = 1.0/3.0;
       // Euclidean distance (not its square), so that phases close to B,
       // e.g. planar at 1/(2 sqrt 3) = 0.289, stand out from B in images.
       DB[X] = (!smallGap && R2 > 0.5)
	       ? sqrt((R1-1.0)*(R1-1.0) + (R2-1.0)*(R2-1.0)
	              + (R3-third)*(R3-third) + (R4-third)*(R4-third) + (R5-third)*(R5-third))
	       : DB_MASKED;
      }
    
    /* >>>>>>>>> phase marking logic <<<<<<<<  */
    if (abs(R1-0.0) <= config.ptol)
      {
	// A-phase OP onsite
	if (
	    (abs(R2-1.0) <= config.ptol)
	    && (abs(R3-0.0) <= config.ptol)
	    && (abs(R4-1.0) <= config.ptol)
	    && (abs(R5-1.0) <= config.ptol)
	   )
	  phaseMarker[X] = 9.0f;
	// alpha-state OP onsite
        else if (
                 (abs(R2-1.0) <= config.ptol)
	         && (abs(R3-(1.0/3.0)) <= config.ptol)
	         && (abs(R4-(1.0/3.0)) <= config.ptol)
	         && (abs(R5-(1.0/3.0)) <= config.ptol)
                )
	  phaseMarker[X] = 6.0f;
	// bipolar-state OP onsite
	else if (
                 (abs(R2-1.0) <= config.ptol)
	         && (abs(R3-(1.0/2.0)) <= config.ptol)
	         && (abs(R4-(1.0/2.0)) <= config.ptol)
	         && (abs(R5-(1.0/2.0)) <= config.ptol)
                )
	  phaseMarker[X] = 8.0f;
	// beta-state OP onsite
	else if (
                 (abs(R2-1.0) <= config.ptol)
	         && (abs(R3-1.0) <= config.ptol)
	         && (abs(R4-1.0) <= config.ptol)
	         && (abs(R5-0.0) <= config.ptol)
                )
	  phaseMarker[X] = 2.0f;
	// gamma-state OP onsite
	else if (
                 (abs(R2-1.0) <= config.ptol)
	         && (abs(R3-0.0) <= config.ptol)
	         && (abs(R4-1.0) <= config.ptol)
	         && (abs(R5-0.0) <= config.ptol)
                )
	  phaseMarker[X] = 4.0f;
	// normal vacuum
	else if (
                 (abs(R2-0.0) <= config.ptol)
	         && (abs(R3-0.0) <= config.ptol)
	         && (abs(R4-0.0) <= config.ptol)
	         && (abs(R5-0.0) <= config.ptol)		 
                )
	  phaseMarker[X] = 0.0f;
	else
	  phaseMarker[X] = 1.0f;	  
      }
    else if (abs(R1-1.0) <= config.ptol)
      {
	// B-phase OP onsite
        if (
            (abs(R2-1.0) <= config.ptol)
	    && (abs(R3-(1.0/3.0)) <= config.ptol)
	    && (abs(R4-(1.0/3.0)) <= config.ptol)
	    && (abs(R5-(1.0/3.0)) <= config.ptol)
           )
	  phaseMarker[X] = 5.0f;
	// planar-phase OP onsite
	else if (
                 (abs(R2-1.0) <= config.ptol)
	         && (abs(R3-(1.0/2.0)) <= config.ptol)
	         && (abs(R4-(1.0/2.0)) <= config.ptol)
	         && (abs(R5-(1.0/2.0)) <= config.ptol)
                )
	  phaseMarker[X] = 3.0f;
	// polar-phase OP onsite
	else if (
                 (abs(R2-1.0) <= config.ptol)
	         && (abs(R3-1.0) <= config.ptol)
	         && (abs(R4-1.0) <= config.ptol)
	         && (abs(R5-1.0) <= config.ptol)
                )
	  phaseMarker[X] = 7.0f;
	else
	  phaseMarker[X] = 1.0f;
      }
    else
      phaseMarker[X] = 1.0f;

    // gap below threshold: normal (overrides the invariant-based marking)
    if (smallGap) phaseMarker[X] = 0.0f;
    /* >>>>>>>>> phase marking logic <<<<<<<<  */
    
  } // onsites block ends here

  //hila::out0 << "Phases marking done \n";
  
} // phaseMarking() function ends here

