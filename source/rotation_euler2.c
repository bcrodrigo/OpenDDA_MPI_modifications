/*  

   Function based on rotation.c that only rotates the scattering vectors 
   from the lab frame to the target frame

   Used by 'scattered_field_MPI.c'

*/
#include "my_new_library.h"

void rotation_euler2(double phi,double costheta,double psi){
/* The negative ZYZ vector rotation matrix R(-psi)R(-theta)R(-phi)=
    _                    _
   | Rx_{x} Rx_{y} Rx_{z} |
   | Ry_{x} Ry_{y} Ry_{z} |
   |_Rz_{x} Rz_{y} Rz_{z}_|

   Rx_{x}=cos(phi)cos(theta)cos(psi)-sin(phi)sin(psi)
   Rx_{y}=sin(phi)cos(theta)cos(psi)+cos(phi)sin(psi)
   Rx_{z}=-sin(theta)cos(psi)
   Ry_{x}=-cos(phi)cos(theta)sin(psi)-sin(phi)cos(psi)
   Ry_{y}=-sin(phi)cos(theta)sin(psi)+cos(phi)cos(psi)
   Ry_{z}=sin(theta)sin(psi)
   Rz_{x}=cos(phi)sin(theta)
   Rz_{y}=sin(phi)sin(theta)
   Rz_{z}=cos(theta)

   The negative ZYZ vector rotation is given by:
          _                    _  _  _
         | Rx_{x} Rx_{y} Rx_{z} || xx |
         | Ry_{x} Ry_{y} Ry_{z} || xy |
         |_Rz_{x} Rz_{y} Rz_{z}_||_xz_| */

   double cosphi,sinphi,theta,sintheta,cospsi,sinpsi,rotation[9];

   /* Euler theta is distributed in cos(theta) */
   theta=radians_to_degrees(acos(costheta));
   cosphi=cos(degrees_to_radians(phi)); /* Convert phi to radians and calculate the cosine */
   sinphi=sin(degrees_to_radians(phi)); /* Convert phi to radians and calculate the sine */
   sintheta=sin(degrees_to_radians(theta)); /* Convert theta to radians and calculate the sine */
   cospsi=cos(degrees_to_radians(psi)); /* Convert psi to radians and calculate the cosine */
   sinpsi=sin(degrees_to_radians(psi)); /* Convert psi to radians and calculate the sine */

   /* Perform a negative rotation */
   rotation[0]=cosphi*costheta*cospsi-sinphi*sinpsi; /* Rx_{x} */
   rotation[1]=sinphi*costheta*cospsi+cosphi*sinpsi; /* Rx_{y} */
   rotation[2]=-sintheta*cospsi; /* Rx_{z} */
   rotation[3]=-cosphi*costheta*sinpsi-sinphi*cospsi; /* Ry_{x} */
   rotation[4]=-sinphi*costheta*sinpsi+cosphi*cospsi; /* Ry_{y} */
   rotation[5]=sintheta*sinpsi; /* Ry_{z} */
   rotation[6]=cosphi*sintheta; /* Rz_{x} */
   rotation[7]=sinphi*sintheta; /* Rz_{y} */
   rotation[8]=costheta; /* Rz_{z} */

/* Rotate the scattering vectors LF to TF */
   
   /* Rotate the scattering vector parallel to the scattering plane LF to TF
   xx=scattering_LF.theta_sca[0]
   xy=scattering_LF.theta_sca[1]
   xz=scattering_LF.theta_sca[2] */
   scattering_TF.theta_sca[0]=scattering_LF.theta_sca[0]*rotation[0]+scattering_LF.theta_sca[1]*rotation[1]+
                              scattering_LF.theta_sca[2]*rotation[2]; /* theta_{sca}x in the target frame */
   scattering_TF.theta_sca[1]=scattering_LF.theta_sca[0]*rotation[3]+scattering_LF.theta_sca[1]*rotation[4]+
                              scattering_LF.theta_sca[2]*rotation[5]; /* theta_{sca}y in the target frame */
   scattering_TF.theta_sca[2]=scattering_LF.theta_sca[0]*rotation[6]+scattering_LF.theta_sca[1]*rotation[7]+
                              scattering_LF.theta_sca[2]*rotation[8]; /* theta_{sca}z in the target frame */
   
   /* Rotate the scattering vector perpendicular to the scattering plane LF to TF
   xx=scattering_LF.phi_sca[0]
   xy=scattering_LF.phi_sca[1]
   xz=scattering_LF.phi_sca[2] */
   scattering_TF.phi_sca[0]=scattering_LF.phi_sca[0]*rotation[0]+scattering_LF.phi_sca[1]*rotation[1]+
                              scattering_LF.phi_sca[2]*rotation[2]; /* phi_{sca}x in the target frame */
   scattering_TF.phi_sca[1]=scattering_LF.phi_sca[0]*rotation[3]+scattering_LF.phi_sca[1]*rotation[4]+
                              scattering_LF.phi_sca[2]*rotation[5]; /* phi_{sca}y in the target frame */
   scattering_TF.phi_sca[2]=scattering_LF.phi_sca[0]*rotation[6]+scattering_LF.phi_sca[1]*rotation[7]+
                              scattering_LF.phi_sca[2]*rotation[8]; /* phi_{sca}z in the target frame */
   
   /* Rotate the scattering direction LF to TF
   xx=scattering_LF.n_sca[0]
   xy=scattering_LF.n_sca[1]
   xz=scattering_LF.n_sca[2] */
   scattering_TF.n_sca[0]=scattering_LF.n_sca[0]*rotation[0]+scattering_LF.n_sca[1]*rotation[1]+
                           scattering_LF.n_sca[2]*rotation[2]; /* n_{sca}x in the target frame */
   scattering_TF.n_sca[1]=scattering_LF.n_sca[0]*rotation[3]+scattering_LF.n_sca[1]*rotation[4]+
                           scattering_LF.n_sca[2]*rotation[5]; /* n_{sca}y in the target frame */
   scattering_TF.n_sca[2]=scattering_LF.n_sca[0]*rotation[6]+scattering_LF.n_sca[1]*rotation[7]+
                           scattering_LF.n_sca[2]*rotation[8]; /* n_{sca}z in the target frame */
}