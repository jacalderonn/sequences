/************************************************************************** 
*                         sequences.c
* 
* Computes constant neutron star sequences with constant baryonic mass
* 
**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "consts.h"
#include "struct.h"

#include "nrutil.h"
#include "equil.h"
#include "equil_util.h"
#include "findmodel.h"
#include "surface.h"
#include "stableorbit.h"
#include "interpol.h"


/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;
      
  int i, ierr;
  double
    e_min, e_max,
   e_center=1e15,                     /* central en. density */
   B,                            /* Quark Bag Constant */
   K=3.0,                        /* Second parameter in "quark" eos */
   spin_freq=100,                  /* Spin Frequency */
    Gamma_P;                      /* Gamma for polytropic EOS */  
                
  int j;

  int nummodels=20;
  float e_c[4], M_0[4];
  float M0, energy_value, temp_energy, ratio_r = 1.0, ej;
  //double Kfreq, Kfreq_j;

  FILE *fpointer;
  fpointer = fopen("E_M_M0_data.txt", "a");

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */



  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'q':
	/* CHOOSE THE EOS TYPE: EITHER "tab" or "poly" or "quark"
	   (default is tab) */
	sscanf(argv[i+1],"%s",eos_type);
	break;  

      case 'b':
	sscanf(argv[i+1],"%lf",&B);
	B *= 1.602e33*KSCALE;
	break;       

      case 'f':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;

      case 'd':
	/* CHOOSE THE NAME OF THE OUTPUT DIRECTORY */
	sscanf(argv[i+1],"%s",data_dir);
	break;

      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_min);
	if(strcmp(eos_type,"poly")!=0)
	  e_min *= C*C*KSCALE;
	break;

      case 'l':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_max);
	if(strcmp(eos_type,"poly")!=0)
	  e_max *= C*C*KSCALE;
	break;

      case 'n':
	sscanf(argv[i+1],"%d",&nummodels);
	break;

     case 's':
	/* CHOOSE THE SPIN FREQUENCY (HZ) */
	sscanf(argv[i+1],"%lf",&spin_freq);
	printf("spin=%g\n",spin_freq);
	break;

      case 'h': 
	fprintf(stderr,"\nQuick help:\n\n");
	fprintf(stderr,"  -q EOS type (tab)\n"); 
	fprintf(stderr,"     tab   : tabulated \n");
        fprintf(stderr,"     quark : simple quark model \n"); 
	fprintf(stderr,"  -b bag constant in MeV/fm^3 for quark models\n");
	fprintf(stderr,"  -f EOS file \n");
	fprintf(stderr,"  -d directory output goes to \n");
	fprintf(stderr,"  -e lowest central energy density to be used, in gr/cm^3\n");
	fprintf(stderr,"  -l largest central energy density \n");
	fprintf(stderr,"  -h this menu\n\n");
	exit(1);
	break;  
      }
    }


  /* PRINT THE HEADER */
  if(strcmp(eos_type,"tab")==0)
    printf("%s,  MDIVxSDIV=%dx%d\n",eos_file,MDIV,SDIV);
  if(strcmp(eos_type,"quark")==0)
    printf("Quark star with B=%f, MDIVxSDIV=%dx%d\n",B/1.602e33/KSCALE,MDIV,SDIV);

  /* SetUpStar loads in the eos and sets up the grid */
  /* Source code for SetUpStar can be found in findmodel.c */

  ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,
		    &eos, &star);

  //printf("The star infrastructure has been set up! \n");

  e_center = e_min;
  temp_energy = e_center;

  printf("e_c \t Mass \t Mass_0\t Radius\tR-ratio\t Spin\t K freq\n");
  printf("e15 \t Msun \t Msun\t km\t --  \t Hz \t Hz \n");

  //while ( e_center <= e_max ){
  
    ierr = MakeSphere(&eos, &star, e_center);
    rns(ratio_r, e_center, &eos, &star); 
    //ierr = SetSpin(&eos, &star, e_center, spin_freq);
    //ierr = Surface(&eos,&star);
  
    printf("%g \t %4.3f \t %4.3f\t %4.2f \t %.2f \t%4.1f \t%4.1f\n",
	      star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, star.R_e*1e-5, ratio_r, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    fprintf(fpointer, "%g %g %g %g %g %g %g\n", 
        star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, star.R_e*1e-5, ratio_r, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    M0 = star.Mass_0/MSUN;
 while(1){   
 printf("-------------------------------------------------------------\n");
  ratio_r = ratio_r - 0.01;

  for(j=0;j<3;j++){
   ej = temp_energy - 0.01*j; 

   ierr = MakeSphere(&eos, &star, ej);
   rns(ratio_r, ej, &eos, &star); 

   printf("%.4g \t %4.3f \t %4.3f\t %4.2f \t %.2f \t%4.1f \n",
        ej, star.Mass/MSUN, star.Mass_0/MSUN, star.R_e*1e-5, ratio_r, star.Omega/(2.0*PI));

   e_c[j] = ej;
   M_0[j] = star.Mass_0/MSUN;
   }
 //printf("-----------------------------------------------\n");
   //printf("M_0 = %g\n",M0);
   //printf("energies:  %g, %g, %g\n",e_c[0], e_c[1], e_c[2]);
   //printf("M_0's: %g, %g, %g\n",M_0[0], M_0[1], M_0[2]);

   energy_value = polyinter(M0, e_c, M_0);

 printf("-------------------------------------------------------------\n");
    temp_energy = energy_value;
    ierr = MakeSphere(&eos, &star, energy_value);
    rns(ratio_r, energy_value, &eos, &star); 
  
    printf("%.4g \t%4.3f \t %4.3f \t %4.2f \t %.2f \t%4.1f \t%4.1f\n",
        energy_value, star.Mass/MSUN, star.Mass_0/MSUN, star.R_e*1e-5, ratio_r, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    if( (star.Omega/(2.0*PI)) > (star.Omega_K/(2.0*PI))) break;

    fprintf(fpointer, "%g %g %g %g %g %g %g\n", 
        energy_value, star.Mass/MSUN, star.Mass_0/MSUN, star.R_e*1e-5, ratio_r, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    //e_center += (e_max-e_min)/(nummodels*1.0);
  //}
  }
  fclose(fpointer);
  return 0;
}









