#include "stdio.h"
#include "math.h"
#include "mersenne_twister.hpp"

#include "multipass3_param.cpp"

int main(void)
{

  // instantiate random number generator
  MTRand rg(1);

  // file writing
  FILE *fpparam, *fpCx, *fpCy, *fpCa, *fpCb;
  char filename[FILENAME_MAX];
 
  // allocate variables
  double t, A, r0, r1, r2, tau, cj, x, K2, s, a, alpha[10], theta, rho;
  double Cx[Nx+1], Cy[Ny+1], Ca[Na+1], Cb[Nb+1];
  double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10;
  int i, j, nx, ny, na, nb, k, l;

  // stoichiometry
  int Deltanx[] = {+1, -1, -2, +2,  0, +1,  0,  0, -1,  0};
  int Deltany[] = { 0,  0, +1, -1,  0,  0,  0,  0,  0,  0};
  int Deltana[] = { 0,  0,  0,  0, +1,  0, -1,  0,  0,  0};
  int Deltanb[] = { 0,  0,  0,  0,  0,  0,  0, +1,  0, -1};

  // record parameters
  sprintf(filename,"%s.param.dat",base);
  fpparam = fopen(filename,"w");
  fprintf(fpparam,"%f\n",thetamin);
  fprintf(fpparam,"%f\n",thetamax);
  fprintf(fpparam,"%f\n",h);
  fprintf(fpparam,"%f\n",nc);
  fprintf(fpparam,"%f\n",lrhomin);
  fprintf(fpparam,"%f\n",lrhomax);
  fprintf(fpparam,"%f\n",T);
  fprintf(fpparam,"%d\n",Ntheta);
  fprintf(fpparam,"%d\n",Nrho);
  fprintf(fpparam,"%d\n",Nx);
  fprintf(fpparam,"%d\n",Ny);
  fprintf(fpparam,"%d\n",Na);
  fprintf(fpparam,"%d\n",Nb);
  fclose(fpparam);

  // loop over rho
  for (l = 0; l < Nrho; l++){
    rho = pow(10,lrhomin + l*(lrhomax-lrhomin)/(Nrho-1));
    printf("rho = %f\n",rho);

    // loop over theta
    for (k = 0; k < Ntheta; k++){
      theta = thetamin + k*(thetamax-thetamin)/(Ntheta-1);
      printf("theta = %f\n",theta);
      
      // initialize count matrices
      for (i = 0; i <= Nx; i++){
	Cx[i] = 0;
      }
      for (i = 0; i <= Ny; i++){
	Cy[i] = 0;
      }
      for (i = 0; i <= Na; i++){
	Ca[i] = 0;
      }
      for (i = 0; i <= Nb; i++){
	Cb[i] = 0;
      }
      
      // convert from critical to Schlogl parameters
      K2 = 3*(theta+1)*nc*nc;
      s =3*nc;
      a = (3*(theta+h)+1)/3/(theta+1)*nc;
      
      // set rates
      k1 = a;
      k2 = 1;
      k3 = rho/nc;
      k4 = rho;
      k5 = rho;
      k6 = s*nc/K2;
      k7 = rho;
      k8 = rho;
      k9 = nc/K2;
      k10 = rho;
      
      // open
      sprintf(filename,"%s_rho%d_theta%d.Cx.dat",base,l,k);
      fpCx = fopen(filename,"w");
      sprintf(filename,"%s_rho%d_theta%d.Cy.dat",base,l,k);
      fpCy = fopen(filename,"w");
      sprintf(filename,"%s_rho%d_theta%d.Ca.dat",base,l,k);
      fpCa = fopen(filename,"w");
      sprintf(filename,"%s_rho%d_theta%d.Cb.dat",base,l,k);
      fpCb = fopen(filename,"w");
      
      // initialize time and numbers
      t = 0;
      nx = round(nc);
      ny = round(nc);
      na = round(nc);
      nb = round(nc);
      
      // simulate
      while (t < T){
	
	// update propensities
	alpha[0] = k1;
	alpha[1] = k2*nx;
	alpha[2] = k3*nx*(nx-1);
	alpha[3] = k4*ny;
	alpha[4] = k5*ny;
	alpha[5] = k6*na;
	alpha[6] = k7*na;
	alpha[7] = k8*ny;
	alpha[8] = k9*nb*nx;
	alpha[9] = k10*nb;
	A = 0;
	for (i = 0; i < 10; i++){
	  A = A + alpha[i];
	}
	
	// choose time tau
	r1 = rg.randDblExc();
	tau = 1/A*log(1/r1);
	
	// choose reaction j
	r2 = rg.randDblExc();
	j = 0;
	cj = alpha[0]/A;
	while (cj < r2){
	  j++;
	  cj = cj + alpha[j]/A;
	}
	
	// update time
	t = t + tau;
	
	// update count vectors
	if (nx <= Nx){
	  Cx[nx] = Cx[nx] + tau;
	}
	if (ny <= Ny){
	  Cy[ny] = Cy[ny] + tau;
	}
	if (na <= Na){
	  Ca[na] = Ca[na] + tau;
	}
	if (nb <= Nb){
	  Cb[nb] = Cb[nb] + tau;
	}
	
	// update molecule numbers
	nx = nx + Deltanx[j];
	ny = ny + Deltany[j];
	na = na + Deltana[j];
	nb = nb + Deltanb[j];
	
      }
      
      // record
      for (i = 0; i <= Nx; i++){
	fprintf(fpCx,"%f\n",Cx[i]);
      }
      for (i = 0; i <= Ny; i++){
	fprintf(fpCy,"%f\n",Cy[i]);
      }
      for (i = 0; i <= Na; i++){
	fprintf(fpCa,"%f\n",Ca[i]);
      }
      for (i = 0; i <= Nb; i++){
	fprintf(fpCb,"%f\n",Cb[i]);
      }
      
      // close
      fclose(fpCx);
      fclose(fpCy);
      fclose(fpCa);
      fclose(fpCb);
    }
  }
  
  return 0;
  
}
