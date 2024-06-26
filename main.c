/*Nomenclature

Kf            - thermal conductivity of base fluid
Kp            - thermal conductivity of Nano-particle
Kef           - Effective thermal conductivity of Nanofluid
phi           - Volume fraction
beta_f        - Volumetric thermal expansion of base fluid
beta_s        - Volumetric thermal expansion of solid
beta_nf       - Volumetric thermal expansion of  Nanofluid
Muf           - Viscosity of non-Newtonian fluid
Mu_nf         - Effective Viscosity of Non-Newtonian nannofluid
rho_f         - Density of base fluid
rho_s         - Density of solid nano-particle
rho_nf        - Density of nanofluid
Cpf           - Heat capacity of base fluid
Cps           - Heat capacity of solid nano-particle
rhoCp_nf      - Volumetric Heat capacity of nanofluid
alpha_nf      - Thermal Diffusivity of nanofluid
N             - Consistency index
nu_nf         - Kinematic viscosity of nanofluid
Pr_nf         - Prandtl number of nanofluid
PrStar        - prandtl number of non-Newtonian nanofluid
Ra_nf         - Rayleigh number of non-Newtonian nanofluid
Nu_nf         - Nusselt number
n             - power law index
L             - Length of Square Geometry
Th            - Temperature of hot wall
Tc            - Temperature of cold wall
g             - gravity
u             - velocity of nanofluid in x direction
v             - velocity of nanofluid in y direction
tau           - stress
psi           - stream function
zeta          - Vorticity
theta         - Temperature

*/
/*
Code for Non-Newtonian Nanofluid

Solution of Natural convection steady state nanofluid flow inside enclosure

Ajay vallabh (16205402)

Task performed
1. Computation of stream function.
2. Computation of temperature field
3. Velocity profiles
4. Local Nusselt Number on hot wall

Files:
input.txt:
          first line  :  Kf, Kp, phi
          Second line : beta_f, beta_s
          Third line  : consistency index(N) , power law index (n)
          fourth line : Density of base fluid(rho_f) , Density of solid nano-particle (rho_s)
          fifth line  : Heat capacity of base fluid(Cpf),Heat capacity of solid nano-particle(Cps)
          sixth line  : Length(L),Height(H),Temperature of hot wall(Th),Temperature of cold wall(Tc),gravity(g)
          seventh line: delta_x, delta_y ,delta_t (time step for unsteady solution)
          eight line  : Precision, Relaxation_parameter, number of time steps to skip while printing the status on screen during execution

For e.g. :

    0.56 200 0.04
    0.87e-4 42e-6
    2.1 1.0
    1000 19300
    4185.5 125.6
    1.0 1.0 45 25 9.81
    0.025 0.025  0.005
    1e-4 1.5 50

Note: Commas are not included

Output files:

     1. output1 - psi
     2. output2 -zeta
     3. output3 -theta
     4. output4 - u
     5  output5 - v

Format:
	x
	y
	u(x,y)
	theta(x,y)
	psi(x,y)
	v(x,y)
	zeta(x,y)
  Algorithm:
     Gauss Seidel with successive over-relaxation, and transient/pseudo-transient method to reach steady state solution

*/

#include"nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NanofluidProperties.h"
int main()
{
    int k, l, skip_itr,iter ;
	// k is the number of nodes along x-axis and l is the same along y-axis
    // skip_iter - number of iterations to skip before priting running status on screen
double  dx, dy,dt, NuAvg, sume, sumo;
double **psi, **psi_old, **tau_xx, **tau_xy, **tau_yy, **zeta, **zeta_old, **theta, **theta_old,**u,**v,**Nu_nf ;
/*     zeta     - vorticity
	   theta 	- temperature
	   u 		- u_x, x-component of velocity
	   v 		- u_y y-component of velocity
	   psi 		- stream function
       delta_x = delta_y = grid size
	   relaxation - relaxation factor in Gauss Seidel algorithm
	   Precision  - precision required in solution i.e. criterion for convergence of solution
	   Ra_nf 		- Rayleigh_number
	   Pr_nf 		- Prandtl_number
	   PrStar       - Prandtl number for non-Newtonian fluid
	   Nu_nf        - Nusselt number
	   tau_xx       - normal mean stress in x direction
	   tau_yy       - normal mean stress in y direction
	   tau_xy       - mean shear stress in y direction with x plane
*/
int i, j  ;  // loop counters
int f ;
/*fluid properties*/
     double Kf, Kp, Kef, phi, beta_f, beta_s, beta_nf ;
     double N, n, rho_f, rho_s, rho_nf, Cpf, Cps,rhoCp_nf ;
     double Muf, Mu_nf , alpha_nf, nu_nf;
	 double relaxation, precision ;
     double L,H,Th, Tc, gravity ;
     double  PrS, Ra_nf ;

double error1,error2,error3,r1,r2,r3,error;
double ae,aw,as,an,ap;
double ae1,as1,aw1,an1,ap1;
double ae2,as2,aw2,an2,ap2,RHS2,STRESS,STRESS1,STRESS2,STRESS3,STRESS4;
double iue,iuw,ivs,ivn;
double q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13;

FILE  *input,*output1, *output2, *output3, *output4, *output5, *output6,*output7,*output8,*output9,*output10;
output1=fopen("psi.txt","w");
output2=fopen("zeta.txt","w");
output3=fopen("theta.txt","w");
output4=fopen("nusselt.txt","w");
output5=fopen("Uvelocity.txt","w");
output6=fopen("Vvelocity.txt","w");
output7=fopen("NuAvg.txt","w");
output8=fopen("tauxx.txt","w");
output9=fopen("tauyy.txt","w");
output10=fopen("tauxy.txt","w");
input=fopen( "input.txt", "r" ) ;
 fscanf(input, "%lf %lf %lf", &Kf, &Kp, &phi ) ;
 fscanf(input, "%lf %lf", &beta_f,&beta_s ) ;
 fscanf(input, "%lf %lf", &N, &n ) ;
 fscanf(input, "%lf %lf", &rho_f, &rho_s ) ;
 fscanf(input, "%lf %lf", &Cpf, &Cps ) ;
 fscanf(input, "%lf %lf %lf %lf %lf", &L, &H,&Th, &Tc, &gravity ) ;
 fscanf(input, "%lf %lf %lf", &dx, &dy, &dt ) ;
 fscanf(input, "%lf %lf %d", &precision, &relaxation, &skip_itr ) ;
k=ceil(1.0/dx+1);
l=ceil(1.0/dy+1);
f=(k-1)*(k-1);
sume=0;
sumo=0;

     Kef = Effective_Thermal_Conductivity(Kf,Kp,phi);
     rho_nf = Nanofluid_Density(rho_f,rho_s,phi);
     rhoCp_nf = Heat_Capacity(rho_f, Cpf, rho_s, Cps, phi);
     beta_nf = Volumetric_Thr_Exp( rho_f, beta_f, rho_s, beta_s, rho_nf, phi);
     alpha_nf = Thermal_Diffusivity( Kef, rhoCp_nf);
     nu_nf = Kinematic_Viscosity( N, rho_nf );
     PrS = PrandtlStar_Number( N, n, L, alpha_nf, rho_nf);
     Ra_nf = Rayleigh_Number( alpha_nf, nu_nf, beta_nf , Th, Tc, gravity, L, n);
     //PrS=7.02;
     //Ra_nf=1000000.0;
     printf("kef=%g\n",Kef);
     printf("Rhonf=%g\n",rho_nf);
     printf("RhoCp_nf=%g\n",rhoCp_nf);
     printf("betanf=%g\n",beta_nf);
     printf("alphanf=%g\n",alpha_nf);
     printf("viscositynf=%g\n",nu_nf);
     printf("Rayleigh=%lf\n",Ra_nf);
     printf("PrandtlS=%lf\n",PrS);
     printf("grid=%d\n",k);

/*define matrix*/
psi      =dmatrix(1,k,1,l);
psi_old  =dmatrix(1,k,1,l);
zeta     =dmatrix(1,k,1,l);
zeta_old =dmatrix(1,k,1,l);
theta    =dmatrix(1,k,1,l);
theta_old=dmatrix(1,k,1,l);
tau_xx   =dmatrix(1,k,1,l);
tau_xy   =dmatrix(1,k,1,l);
tau_yy   =dmatrix(1,k,1,l);
u        =dmatrix(1,k,1,l);
v        =dmatrix(1,k,1,l);
Nu_nf    =dmatrix(1,k,1,l);

/*initialize 2 D matrix*/
for(i=1;i<=k;i++)
{
    for(j=1;j<=l;j++)
    {
    psi[i][j]=0;
    psi_old[i][j]=0;
    theta[i][j]=0;
    theta_old[i][j]=0;
    zeta[i][j]=0;
    zeta_old[i][j]=0;
    u[i][j]=0;
    v[i][j]=0;
    tau_xx[i][j]   = 0;
    tau_xy[i][j]   = 0;
    tau_yy[i][j]   = 0;
    }
}
/*boundary condition */
for(j=2;j<=l-1;j++)
{
 psi_old[1][j] = 0;
 theta_old[1][j] = 1.0;
 zeta_old[1][j] = -2.0*( psi_old[2][j] )*f;
 psi_old[k][j] = 0;
 theta_old[k][j] = 0;
 zeta_old[k][j] = -2.0*( psi_old[k-1][j] )*f;
}
for(i=1;i<=k;i++)
 {
 theta_old[i][l]=theta_old[i][l-1];
 psi_old[i][l]=0;
 zeta_old[i][l]=-2.0*(psi_old[i][l-1])*f;
 psi_old[i][1]=0;
 theta_old[i][1]=theta_old[i][2];
 zeta_old[i][1]=-2.0*(psi_old[i][2])*f;
 }
 /*copy 2d array*/
 for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
        zeta[i][j]=zeta_old[i][j];
        psi[i][j]=psi_old[i][j];
        theta[i][j]=theta_old[i][j];
     }
 }
 error2=1.0;
 error3=1.0;
 error=1.0;
 iter=0;
 while(error>precision)
 {
     iter++;
 error3=0;
 error2=0;
 error=0;
  error1=1.0;
  while(error1>precision)
  {
      error1=0;
   /*copy 2d array*/
  for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
        psi_old[i][j]=psi[i][j];
 }
 //constants
      ae=-dt*f ;
      aw=-dt*f ;
      as=-dt*f ;
      an=-dt*f ;
      ap=(1-(ae+aw+as+an));
 /*inner iteration begin */

  for(i=2;i<=k-1;i++)
 {
     for(j=2;j<=l-1;j++)
     {
      psi[i][j]=(dt*zeta[i][j] + psi_old[i][j] - (ae * psi[i+1][j] + aw * psi[i-1][j] + as * psi[i][j-1] + an * psi[i][j+1] ))/ap;
     r1=fabs(psi[i][j]-psi_old[i][j]);
     if(r1>error1)
        error1=r1;
     }
 }
  }
 /*calculate velocity profile*/
 for(i=2;i<=k-1;i++)
 {
     for(j=2;j<=l-1;j++)
     {
         u[i][j]=(psi[i][j+1] - psi[i][j-1])*(l-1)/(2.0) ;
         v[i][j]=(psi[i-1][j] - psi[i+1][j])*(k-1)/(2.0) ;
     }
 }

 /*calculate shear stress value*/
  for(i=2;i<=k-1;i++)
 {
     for(j=2;j<=l-1;j++)
     {
      q1=(u[i+1][j]-u[i-1][j])*(k-1);
      q2=(u[i][j+1]-u[i][j-1])*(l-1)/(2.0);
      q3=(v[i+1][j]-v[i-1][j])*(k-1)/(2.0);
      q4=(v[i][j+1]-v[i][j-1])*(l-1);
      q5=q2+q3;
     if(fabs(q5)<0.001)
        tau_xy[i][j]=0;
     else
        tau_xy[i][j]=pow(fabs(q5),(n-1))*q5-q5;

    if(fabs(q1)<0.001)
        tau_xx[i][j]=0;
     else
        tau_xx[i][j]=pow(fabs(q1),(n-1))*q1-q1;

    if(fabs(q4)<0.001)
        tau_yy[i][j]=0;
     else
        tau_yy[i][j]=pow(fabs(q4),(n-1))*q4-q4;
     }
 }
 /*boundary condition for shear stress*/
 /*left wall*/
 for(j=2;j<=l-1;j++)
 {
     q6=(-3.0*u[1][j]+4.0*u[2][j]-u[3][j])*(k-1);
     q7=(-3.0*v[1][j]+4.0*v[2][j]-v[3][j])*(k-1)/(2.0);

    if(fabs(q6)<0.001)
        tau_xx[1][j]=0;
     else
        tau_xx[1][j]=pow(fabs(q6),(n-1))*q6-q6;

    if(fabs(q7)<0.001)
        tau_xy[1][j]=0;
     else
        tau_xy[1][j]=pow(fabs(q7),(n-1))*q7-q7;
 }
  /*right wall*/
 for(j=2;j<=l-1;j++)
 {
     q8=(u[k-2][j]-4.0*u[k-1][j]+3.0*u[k][j])*(k-1);
     q9=(v[k-2][j]-4.0*v[k-1][j]+3.0*v[k][j])*(k-1)/(2.0);

    if(fabs(q8)<0.001)
        tau_xx[k][j]=0;
     else
        tau_xx[k][j]=pow(fabs(q8),(n-1))*q8-q8;

    if(fabs(q9)<0.001)
        tau_xy[k][j]=0;
     else
        tau_xy[k][j]=pow(fabs(q9),(n-1))*q9-q9;
 }
 /*lower wall*/
 for(i=1;i<=k;i++)
 {
  q10=(-3.0*v[i][1]+4.0*v[i][2]-v[i][3])*(l-1);
  q11=(-3.0*u[i][1]+4.0*u[i][2]-u[i][3])*(l-1)/(2.0);

     if(fabs(q10)<0.001)
        tau_yy[i][1]=0;
     else
        tau_yy[i][1]=pow(fabs(q10),(n-1))*q10-q10;

    if(fabs(q11)<0.001)
        tau_xy[i][1]=0;
     else
        tau_xy[i][1]=pow(fabs(q11),(n-1))*q11-q11;
 }
 /*upper wall*/
 for(i=1;i<=k;i++)
 {
  q12=(v[i][l-2]-4.0*v[i][l-1]+3.0*v[i][l])*(l-1);
  q13=(u[i][l-2]-4.0*u[i][l-1]+3.0*u[i][l])*(l-1)/(2.0);

     if(fabs(q12)<0.001)
        tau_yy[i][k]=0;
     else
        tau_yy[i][k]=pow(fabs(q12),(n-1))*q12-q12;

    if(fabs(q13)<0.001)
        tau_xy[i][k]=0;
     else
        tau_xy[i][k]=pow(fabs(q13),(n-1))*q13-q13;
 }

/*vorticity equation*/
for(i=2;i<=k-1;i++)
{
  for(j=2;j<=l-1;j++)
  {
   if(u[i][j]>0)
   {
       iuw=1.0;
       iue=0;
   }
   else
   {
       iuw=0;
       iue=1.0;
   }
   if(v[i][j]>0)
   {
       ivs=1.0;
       ivn=0;
   }
   else
   {
       ivs=0;
       ivn=1.0;
   }
    ae1 =-(dt*f+iue*fabs(u[i][j]*dt*(k-1)));
    aw1 =-(dt*f+iuw*fabs(u[i][j]*dt*(k-1)));
    as1 =-(dt*f+ivs*fabs(v[i][j]*dt*(l-1)));
    an1 =-(dt*f+ivn*fabs(v[i][j]*dt*(l-1)));
    ap1 = (1.0-(ae1+aw1+as1+an1));

    theta[i][j]=(theta_old[i][j]-(ae1*theta[i+1][j]+aw1*theta[i-1][j]+as1*theta[i][j-1]+an1*theta[i][j+1]))/ap1;

    ae2 =-((PrS/pow((1-phi),2.5))*dt*f+iue*fabs(u[i][j]*dt*(k-1)));
    aw2 =-((PrS/pow((1-phi),2.5))*dt*f+iuw*fabs(u[i][j]*dt*(k-1)));
    as2 =-((PrS/pow((1-phi),2.5))*dt*f+ivs*fabs(v[i][j]*dt*(l-1)));
    an2 =-((PrS/pow((1-phi),2.5))*dt*f+ivn*fabs(v[i][j]*dt*(l-1)));
    ap2 = (1.0-(ae2+aw2+as2+an2));

    STRESS1 = 0.25 * (tau_yy[i+1][j+1] - tau_yy[i-1][j+1] - tau_yy[i+1][j-1] + tau_yy[i-1][j-1]);
    STRESS2 = 0.25 * (tau_xx[i+1][j+1] - tau_xx[i-1][j+1] - tau_xx[i+1][j-1] + tau_xx[i-1][j-1]);
    STRESS3 = tau_xy[i+1][j] - 2.0 * tau_xy[i][j] + tau_xy[i-1][j];
    STRESS4 = tau_xy[i][j+1] - 2.0 * tau_xy[i][j] + tau_xy[i][j-1];
    STRESS  = STRESS1-STRESS2+STRESS3-STRESS4 ;
    RHS2=(Ra_nf*PrS)*(theta[i+1][j]-theta[i-1][j])*(k-1)/(2.0) ;
    zeta[i][j] = (RHS2*dt+zeta_old[i][j] - (ae2 * zeta[i+1][j] + aw2*zeta[i-1][j] + as2 * zeta[i][j-1] + an2 * zeta[i][j+1])+(dt*f*PrS/pow((1-phi),2.5))*(STRESS))/ap2;

  }
}
 /*boundary condition*/
 for(j=2;j<=l-1;j++)
 {
     zeta[1][j]=-2.0 * f * (psi[2][j]);
     zeta[k][j]=-2.0 * f * (psi[k-1][j]);
 }
 for(i=1;i<=k;i++)
 {
     zeta[i][1]=-2.0 * f * (psi[i][2]);
     zeta[i][l]=-2.0 * f * (psi[i][l-1]);
     theta[i][1]=(2.0*theta[i][4]-9.0*theta[i][3]+18.0*theta[i][2])/11.0;       //third order forward difference
     theta[i][l]=(2.0*theta[i][l-3]-9.0*theta[i][l-2]+18.0*theta[i][l-1])/11.0;
 }
 /*copy array*/
 for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
         zeta_old[i][j] = relaxation * zeta[i][j] + (1 - relaxation) * zeta_old[i][j] ;
         theta_old[i][j]= relaxation * theta[i][j] + (1 - relaxation) * theta_old[i][j];
     }
 }
 /*convergence criteria */
 for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
  r2=fabs(zeta[i][j]-zeta_old[i][j]);
  r3=fabs(theta[i][j]-theta_old[i][j]);
  if(r2>error2)
    error2=r2;
  if(r3>error3)
    error3=r3;
  if(error2>error3)
    error=error2;
  else
    error=error3;
     }
 }
 if(iter%10==0)
 printf("itr=%d \t error=%lf\n",iter,error);
 if(iter>80000)
    break;
 }
 /*calculate nusselt number at the hot wall*/

 for(j=1; j<=l; j++)
 {
     Nu_nf[1][j]=(k-1)*(2.0*theta[4][j]-9.0*theta[3][j]+18.0*theta[2][j]-11.0*theta[1][j])/(6.0);
 }

/*calculate avg nusselt number*/
for(j=2;j<l-1;j=j+2)
{
    sume=sume+Nu_nf[1][j];
}
for(j=3;j<l-1;j=j+2)
{
    sumo=sumo+Nu_nf[1][j];
}
NuAvg=1.0/((k-1)*3.0)*(Nu_nf[1][1]+Nu_nf[1][l]+4.0*sumo+2.0*sume);
 /*print values */
 for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
         fprintf(output1,"%lf,",psi[i][j]);
         fprintf(output2,"%lf,",zeta[i][j]);
         fprintf(output3,"%lf,",theta[i][j]);
         fprintf(output5,"%lf,",u[i][j]);
         fprintf(output6,"%lf,",v[i][j]);
         fprintf(output8,"%lf,",tau_xx[i][j]);
         fprintf(output9,"%lf,",tau_yy[i][j]);
         fprintf(output10,"%lf,",tau_xy[i][j]);
     }
     fprintf(output1,"\n");
     fprintf(output2,"\n");
     fprintf(output3,"\n");
     fprintf(output5,"\n");
     fprintf(output6,"\n");
     fprintf(output8,"\n");
     fprintf(output9,"\n");
     fprintf(output10,"\n");
 }
for(j=1;j<=l;j++)
fprintf(output4,"%lf \n",Nu_nf[1][j]);

fprintf(output7,"Avg nusselt=%lf",NuAvg);
 fclose(output1);
 fclose(output2);
 fclose(output3);
 fclose(output4);
 fclose(output5);
 fclose(output6);
 fclose(output7);
 fclose(output8);
 fclose(output9);
 fclose(output10);
 fclose(input);
    return 0;
}
