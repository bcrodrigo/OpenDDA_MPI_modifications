/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Calculates the various efficiencies

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License 


   MODIFICATION:
   Efficiencies are written into a csv file, instead of the standard output

*/
#include "efficiencies.h"

void efficiencies(void){

   /* Geometrical cross section pi*radius^{2} */
   cross_section.geometrical=dpi*radius.value*radius.value;

   /* Q_{ext} Extinction efficiency */
   efficiency.extinction=cross_section.extinction/cross_section.geometrical;

   /* Q_{abs} Absorption efficiency */
   efficiency.absorption=cross_section.absorption/cross_section.geometrical;

   /* Q_{sca} Scattering efficiency */
   efficiency.scattering=cross_section.scattering/cross_section.geometrical;

   fprintf(q_output,"%+.*g,%+.*g,%+.*g,%+.*g,%+.*g,%d,%+.*g,%+.*g,%+.*g\n",DBLP,wavelength.value,DBLP,radius.value,DBLP,euler_theta.value,DBLP,euler_psi.value,DBLP,euler_phi.value,polarisation_state,DBLP,efficiency.extinction,DBLP,efficiency.absorption,DBLP,efficiency.scattering);
      
}
