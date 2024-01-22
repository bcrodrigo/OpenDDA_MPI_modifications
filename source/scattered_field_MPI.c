/*

  Function that calculates the scattered field for a set of scattering angles (phi, theta) 
  It outputs the element of a dimensionless f matrix relating incident and scattered electic fields according to

    B. T. Draine. Astrophys. J., 333:848–872, 1988 (https://doi.org/10.1086/166795)

  Note that the elements of the f matrix are related to the Amplitude scattering matrix S as defined in

    Craig F. Bohren and Donald R. Huffman. Absorption and Scattering of Light by Small Particles. Wiley-VCH Verlag GmbH, 2007.
  
  See my thesis [http://hdl.handle.net/1974/14915], section 3.2, equation 3.26 on how to relate S and f.

*/
#include "scattered_field_MPI.h"
#include "rotation_euler2.h"

void scattered_field_MPI(void){
  // variables needed to broadcast the scattering parameters
  long begin;
  int i,count=3;
  int blocklens[3]={3,3,3};
  MPI_Datatype types[3]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  MPI_Aint displacements[3];
  MPI_Datatype broadcast_scatt_type;

  // Define a structure for broadcasting the physical parameters
  struct broadcast_scatt {
      double phi_TF_value[3];
      double theta_TF_value[3];
      double nsca_TF_value[3];
  };
  struct broadcast_scatt scatt;

  long int index0i,index0k,index2; // some long indices
  int index_occupied[parallel.alloc_vector]; // array that will store the index of occupied lattice sites.
  int x,y,z; // coordinates of the occupied lattice site
  int j,k; // just some indices

  dcomplex Ed[3],Ephi,Etheta; // array that accumulates the field produced by the dipoles
  dcomplex dot_inexp;// dot product between the scattering direction and the j-th dipole location.
  dcomplex exp_dipoles[parallel.alloc_vector]; // complex exponential with dot product between dipole locations and scattering direction
  
  double temp_re,temp_im; // variables that store the value of Re{Ed} and Im{Ed} temporarily
  double k3; // k^3, k the wavenumber
  int phi_index,theta_index;//indices for phi,theta
  double ReEd,ImEd; // Real and imaginary part of Ed, just for reduction purposes

  dcomplex temp_scalephi,temp_scaletheta; 
  
  if (parallel.myrank==0) printf("\n Calculating the Scattered field\n");

  // Create the 'scatt' structure, to broadcast the scattering directions to everyone
  MPI_Get_address(&scatt,&displacements[0]);
  MPI_Get_address(&scatt.theta_TF_value,&displacements[1]);
  MPI_Get_address(&scatt.nsca_TF_value,&displacements[2]);
  begin=displacements[0];
  for(i=0;i<count;i++){ displacements[i]-=begin;}

  MPI_Type_create_struct(count,blocklens,displacements,types,&broadcast_scatt_type);
  MPI_Type_commit(&broadcast_scatt_type);

  // Recover the index of all the occupied lattice sites, on each processor
  index2=0;
  for(k=0;k<target.P;k++){
    if(k%parallel.np==parallel.myrank){
         parallel.plane=((int)((double)k/parallel.dnp));
         index0k=k*target.M;
         for(j=0;j<target.M;j++){
            index0i=index0k+j;
            if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ // Lattice site is occupied 
		index_occupied[index2]=index0i; // store the index that corresponds to an occupied lattice site
		index2++;		
          } // end if target.occupied
         } // end for j
      } //end if k%parallel.np
   } // end for k

  k3=wavenumber*wavenumber;	
  k3=k3*wavenumber;

  // We want to calculate the scattered field, for all the specified scattering angles

  phi.value=phi.limits[0]; // Minimum value stored in control.input file
  for (phi_index=0;phi_index<phi.list_length;phi_index++) {
    // If the flag is 1, we read the value from the list
    if ((parallel.myrank==0)&&(phi.flag==1)){phi.value=phi.list[phi_index];} 

    theta.value = theta.limits[0];
    
    for (theta_index=0;theta_index<theta.list_length;theta_index++) {

      if ((parallel.myrank==0)&&(theta.flag==1)){theta.value=theta.list[theta_index];}

      if (parallel.myrank==0) {

        scattering_vectors(phi.value,theta.value); // Defines the scattering vectors in the lab frame
        rotation_euler2(euler_phi.value,euler_theta.value,euler_psi.value);// Change scattering vectors from LF to TF

        // broadcast all the scattering parameters
        for (i=0;i<3;i++) {
          scatt.phi_TF_value[i]=scattering_TF.phi_sca[i];
          scatt.theta_TF_value[i]=scattering_TF.theta_sca[i];
          scatt.nsca_TF_value[i]=scattering_TF.n_sca[i];
        }// end for

      }// end restrict to master

      MPI_Bcast(&scatt,1,broadcast_scatt_type,0,MPI_COMM_WORLD);

      if (parallel.myrank!=0){

        for (i=0;i<3;i++) {
          scattering_TF.phi_sca[i]=scatt.phi_TF_value[i];
          scattering_TF.theta_sca[i]=scatt.theta_TF_value[i];
          scattering_TF.n_sca[i]=scatt.nsca_TF_value[i];
        } // end for 

      } // end restrict to slaves

      /*  
        Here we are doing component by component the calculation of the scattered field.
        More specifically, we are using equation 10 from:
        Draine, J. Opt. Soc. Am. A, Vol.11, No.4 April 1994 (https://doi.org/10.1364/JOSAA.11.001491)
      */
  
      if (parallel.myrank==0) {
        // initialize the components of Ed
        for (j=0;j<3;j++) {
          Ed[j].dat[0] = Ed[j].dat[1] = 0.0;
        }

    	  Ephi.dat[0] = Ephi.dat[1] = 0.0;
    	  Etheta.dat[0] = Etheta.dat[1] = 0.0;

        temp_scaletheta.dat[0] = 0.0;
        temp_scaletheta.dat[1] = 0.0;
        temp_scalephi.dat[0] = 0.0;
        temp_scalephi.dat[1] = 0.0;
      } // end restrict to master

      // evaluate the complex exponential with the dipole locations
    	dot_inexp.dat[0]=0.0;
    	for(index0i=0;index0i<parallel.alloc_vector;index0i++) {
        // first we calculate the complex exponential of Ed
    	  index2=index_occupied[index0i];
    	  z=(int)(index2/target.M);
    	  y=(int)((index2-z*target.M)/target.K);
    	  x=index2-z*target.M-y*target.K;
    	  dot_inexp.dat[1]=(-1)*(2*dpi*(target.dipole_spacing/wavelength.value))*(x*scattering_TF.n_sca[0]+y*scattering_TF.n_sca[1]+z*scattering_TF.n_sca[2]);
    	  exp_dipoles[index0i]=dcomplex_no_r_exp(dot_inexp); //evaluate the exponential
    	} // end for loop

      // calculate Ed component by component

    	for (j=0;j<3;j++){

        temp_re = temp_im = 0.0;
    	  
        for (index0i=0;index0i<parallel.alloc_vector;index0i++) {
    	    temp_re = temp_re+exp_dipoles[index0i].dat[0]*dipole_polarisation[index0i+j*parallel.alloc_vector].dat[0]-exp_dipoles[index0i].dat[1]*dipole_polarisation[index0i+j*parallel.alloc_vector].dat[1];
    	    temp_im = temp_im+exp_dipoles[index0i].dat[0]*dipole_polarisation[index0i+j*parallel.alloc_vector].dat[1]+exp_dipoles[index0i].dat[1]*dipole_polarisation[index0i+j*parallel.alloc_vector].dat[0];
      	} // end loop index0i, parallel region

    	  MPI_Reduce(&temp_re,&ReEd,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    	  MPI_Reduce(&temp_im,&ImEd,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    	  if (parallel.myrank==0) {

          Ed[j].dat[0]=ReEd; Ed[j].dat[1]=ImEd;	
          temp_scalephi=dcomplex_scale(Ed[j],scattering_TF.phi_sca[j]);
          Ephi=dcomplex_add(Ephi,temp_scalephi); // Accumulate dot product
          temp_scaletheta=dcomplex_scale(Ed[j],scattering_TF.theta_sca[j]);
          Etheta=dcomplex_add(Etheta,temp_scaletheta); // Accumulate dot product
    	  } // end restrict to master

    	} // end loop j (components of Ed)

    	if (parallel.myrank==0){

        // Scale complex numbers by k3
        temp_scaletheta=dcomplex_scale(Etheta,k3);
        temp_scalephi=dcomplex_scale(Ephi,k3);
        
        /*

          Print the following to a file

          wavelength
          effective_radius
          euler_theta
          euler_psi
          euler_phi
          polarization_state (integer)
          scattering_phi
          scattering_theta
        */

        fprintf(fmatrix_output,"%+.*g,%+.*g,%+.*g,%+.*g,%+.*g,%d,%+.*g,%+.*g,",
          DBLP,wavelength.value,
          DBLP,radius.value,
          DBLP,euler_theta.value,
          DBLP,euler_psi.value,
          DBLP,euler_phi.value,
          polarisation_state,
          DBLP,phi.value,
          DBLP,theta.value);

        /*

          Print the following to a file

          f_1p
          f_2p
          
          f matrix elements f_{1,p} and f_{2,p}
          where p labels the incident polarization
          
          Note each element is complex-valued, and is written in the file as
          
          Re{f_1p}+Im{f_1p}j, Re{f_2p}+Im{f_2p}j

          where j is the imaginary unit (j = sqrt(-1))
        */

        fprintf(fmatrix_output,"%+.*g%+.*gj,%+.*g%+.*gj\n",
          DBLP,temp_scaletheta.dat[0],
          DBLP,temp_scaletheta.dat[1],
          DBLP,temp_scalephi.dat[0],
          DBLP,temp_scalephi.dat[1]);

    	}// end restrict to master 

      if((parallel.myrank==0)&&(theta.flag==0)){theta.value+=theta.limits[1];} 
      // We readed the theta values from 'control.input' so we increase theta.value

    }// end for theta_index

    if((parallel.myrank==0)&&(phi.flag==0)){phi.value+=phi.limits[1];} // Same as with theta.flag and theta.value 

  }// end for phi_index

  MPI_Type_free(&broadcast_scatt_type);

} // END OF PROGRAM