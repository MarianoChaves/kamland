/*
 * conversion between natural and MKSA units
 * length 1/eV = 1.97e-7 m
 * time   1/eV = 6.58e-16 sec
 * mass   1/eV = 1.78e-36 Kg */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <sys/types.h>
#include <string>
#include <sstream>
#include <cstdlib>
#include <complex>
#include <vector>
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <globes/globes.h>   /* GLoBES library */

//char MYFILE1[]="../datfiles/location_juno_ih_heavy.dat";
//FILE *out1 = NULL;
/* ############################################################## */
/* ############################################################## */

using namespace std;

inline double sqr(complex<double> x)
	{ return (pow(real(x),2)+pow(imag(x),2));};


double pi = acos(-1.0);

#define GLB_ReSet 6 
#define GLB_ImSet 7
#define GLB_ReSeu 8
#define GLB_ImSeu 9

#define me 0.511
#define Delta 1.29

#define gT 0.987
#define gA 1.27
#define gS 1.02

/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/

double t12;
double t13;
double t23;
double p13;
double D21;
double D31;
double ReSet;
double ImSet;
double ReSeu;
double ImSeu;

/***************************************************************************
 * Store oscillation parameters in internal data structures.               *
 * For more sophisticated probability engines, this would be the right     *
 * place to pre-compute the mixing matrix and parts of the Hamiltonian in  *
 * order to speed up the calls to the actual probability matrix function.  *
 ***************************************************************************/
int my_set_oscillation_parameters(glb_params p, void *user_data)
{
  t12    = glbGetOscParams(p, GLB_THETA_12);
  t13    = glbGetOscParams(p, GLB_THETA_13);
  t23    = glbGetOscParams(p, GLB_THETA_23);
  p13     = glbGetOscParams(p, GLB_DELTA_CP);
  D21     = glbGetOscParams(p, GLB_DM_21);   
  D31     = glbGetOscParams(p, GLB_DM_31);  
  ReSet   = glbGetOscParams(p, GLB_ReSet); 
  ImSet   = glbGetOscParams(p, GLB_ImSet);
  ReSeu   = glbGetOscParams(p, GLB_ReSeu);
  ImSeu   = glbGetOscParams(p, GLB_ImSeu);

  return 0;
}

double fT(double E)
{
	return(1./(2.12-22.81/pow(E,1./3.))+35.36/sqrt(E)-11.72/E);

}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/
int my_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p, t12, GLB_THETA_12);
  glbSetOscParams(p, t13, GLB_THETA_13);
  glbSetOscParams(p, t23, GLB_THETA_23);
  glbSetOscParams(p, p13, GLB_DELTA_CP);
  glbSetOscParams(p, D21, GLB_DM_21);  /* Convert to eV^2 */
  glbSetOscParams(p, D31, GLB_DM_31);  /* Convert to eV^2 */
  glbSetOscParams(p,ReSet, GLB_ReSet);
  glbSetOscParams(p,ImSet, GLB_ImSet);
  glbSetOscParams(p,ReSeu, GLB_ReSeu);
  glbSetOscParams(p,ImSeu, GLB_ImSeu);
  return 0;
}


/***************************************************************************
 * Calculate oscillation probabilities.                                    *
 * Since for our setup, only P_ee is required, all other entries of P are  *
 * set to zero for simplicity. Furthermore, we neglect matter effects and  *
 * the filter feature (parameter filter_sigma).                            *
 * The formula for P_ee is Eq. (36) from hep-ph/0502147.                   *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:            The buffer where the probabilities are to be stored     *
 *   cp_sign:      +1 if probalities for neutrinos are requested, -1 for   *
 *                 anti-neutrinos.                                         *
 *   E:            The neutrino energy in GeV                              *
 *   psteps:       Number of constant density layers in the matter profile *
 *   length:       The lengths of these layers in km                       *
 *   density:      The individual densities of these layers in g/cm^3      *
 *   filter_sigma: Width of low-pass filter as given in the AEDL file      *
 ***************************************************************************/

int my_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

complex<double> minusX = {-1,0};
complex<double> I = sqrt(minusX);
int k,l;

  for (int i=0; i < 3; i++)
    {for (int j=0; j < 3; j++)
      {P[i][j] = 0.0;}}


	complex<double> Ue[3];
	complex<double> Ueconj[3];
	complex<double> Um[3];
	complex<double> Umconj[3];
	complex<double> Ut[3];
	complex<double> Utconj[3];

	complex<double> UeNSI[3];
	complex<double> UeNSIconj[3];

	double Dm2[3][3];
	complex<double> nsi[3][3];
	
	UeNSI[0]=UeNSI[1]=UeNSI[2]=UeNSIconj[0]=UeNSIconj[1]=UeNSIconj[2]=0.0;

	Dm2[1][0]=D21; Dm2[0][1]=-D21; // dm21
	Dm2[2][0]=D31; Dm2[0][2]=-D31; // dm31

	Dm2[2][1]=Dm2[2][0]-Dm2[1][0]; // dm32=dm31-dm21
	Dm2[1][2]=-Dm2[2][1];

	Dm2[0][0]=Dm2[1][1]=Dm2[2][2]=0.0;


	Ue[0] = cos(t12)*cos(t13);
	Ue[1] = sin(t12)*cos(t13);
	Ue[2] = sin(t13);
	Ueconj[0] = cos(t12)*cos(t13);
	Ueconj[1] = sin(t12)*cos(t13);
	Ueconj[2] = sin(t13);

	Um[0] = -sin(t12);
	Um[1] = cos(t12);
	Um[2] = 0;
	Umconj[0] = -sin(t12);
	Umconj[1] = cos(t12);
	Umconj[2] = 0;

	Ut[0] = -cos(t12)*sin(t13);
	Ut[1] = -sin(t12)*sin(t13);
	Ut[2] = cos(t13);
	Utconj[0] = -cos(t12)*sin(t13);
	Utconj[1] = -sin(t12)*sin(t13);
	Utconj[2] = cos(t13);


//	double Nee;

	nsi[0][0]={0.0,0.0};
	nsi[1][1]={0.0,0.0};
	nsi[2][2]={0.0,0.0};

	nsi[0][1]={ReSeu,ImSeu};
	nsi[1][0]={ReSeu,-ImSeu};

	nsi[0][2]={ReSet,ImSet};
	nsi[2][0]={ReSet,-ImSet};

	nsi[1][2]={0.0,0.0};
	nsi[2][1]={0.0,0.0};


	UeNSI[0]=nsi[0][1]*Um[0]+nsi[0][2]*Ut[0];
	UeNSI[1]=nsi[0][1]*Um[1]+nsi[0][2]*Ut[1];
	UeNSI[2]=nsi[0][1]*Um[2]+nsi[0][2]*Ut[2];

	UeNSIconj[0]=conj(UeNSI[0]);
	UeNSIconj[1]=conj(UeNSI[1]);
	UeNSIconj[2]=conj(UeNSI[2]);

	double dSL=-gS/(1+3*gA*gA)*me/(E-Delta);
	double dSS=gS*gS/(1+3*gA*gA);
	double pTL=-gT/gA*me/fT(E);
	double pTT=gT*gT/(gA*gA);
	double dTL=3*gA*gT/(1+3*gA*gA)*me/(E-Delta);
	double dTT=3*gT*gT/(1+3*gA*gA);

//	Nee=1+dSS*(pow(abs(nsi[1][0]),2)+pow(abs(nsi[2][0]),2));	

	double baseline[20]={160,179,191,138,214,146,88,349,345,295,431,401,561,755,830,783,712,986,735,709};
	double power[20]={24.3,13.7,10.2,4.5,10.6,4.9,1.6,14.2,13.2,3.3,6.5,3.8,6.0,10.1,5.3,3.3,11.5,17.4,9.2,8.2};

	complex<double> Prob=0;

	for (int item=0; item<20; item+=1)
	{
	for(k=0;k<3;k++)
	{
	for(l=0;l<3;l++)
	{	
		double L=baseline[item];
		double Po=(1/181.7)*power[item];
		//Prob+= 0.4698*Po/L/L*(exp(-I*Dm2[k][l]*L/E*(2.*1.267))*(Ueconj[k]*Ue[l])*(Ueconj[l]*Ue[k]+dSL*UeNSI[k]*Ueconj[l]+dSL*UeNSIconj[l]*Ue[k]+dSS*UeNSI[k]*UeNSIconj[l])); //scalar
		Prob+=0.4698*Po/L/L*(exp(-I*Dm2[k][l]*L/E*(2.*1.267))*(Ueconj[k]*Ue[l]+pTL*UeNSIconj[k]*Ue[l]+pTL*Ueconj[k]*UeNSI[l]+pTT*UeNSIconj[k]*UeNSI[l])*(Ueconj[l]*Ue[k]+dTL*UeNSI[k]*Ueconj[l]+dTL*UeNSIconj[l]*Ue[k]+dTT*UeNSI[k]*UeNSIconj[l])); //tensor

	};
	};};

	P[0][0]=real(Prob);

  return 0;
}

/* ############################################################## */
/* ############################################################## */

int main(int argc, char *argv[])
{ 	glbInit(argv[0]);

//out1 = fopen(MYFILE1, "w"); if(out1==NULL){printf("Error opening output file.\n");return -1;}
char AEDLFILE1[]="kamland.glb";
ofstream outfile1("chisq.dat");
ofstream outfile2("bfp.dat");


glbSelectMinimizer(GLB_MIN_NESTED_POWELL);

glbRegisterProbabilityEngine(10, &my_probability_matrix, 
&my_set_oscillation_parameters, &my_get_oscillation_parameters, NULL);

glbClearExperimentList();
glbInitExperiment(AEDLFILE1,&glb_experiment_list[0],&glb_num_of_exps); 

glb_params true_values=glbAllocParams();
glb_params test_values=glbAllocParams();
glbSetDensityParams(true_values,1.0,GLB_ALL);
glbSetDensityParams(test_values,1.0,GLB_ALL);
glb_params input_errors = glbAllocParams();
glb_projection my_projection = glbAllocProjection();  

double t12true=asin(sqrt(0.32)); //double t12stddev=0.02;
double t13true=asin(sqrt(0.0222)); double t13stddev=0.0076;
double dm21true=7.55e-5; double dm21stddev=0.20e-5;
double dm31true=2.42e-3; //double dm31stddev=0.04e-3; 

glbDefineParams(true_values,t12true,t13true,pi/4,0,dm21true,dm31true);
glbSetOscParams(true_values,0.0,GLB_ReSet);
glbSetOscParams(true_values,0.0,GLB_ImSet);
glbSetOscParams(true_values,0.0,GLB_ReSeu);
glbSetOscParams(true_values,0.0,GLB_ImSeu);

glbDefineParams(test_values,t12true,t13true,pi/4,0,dm21true,dm31true);
glbSetOscParams(test_values,0.0,GLB_ReSet);
glbSetOscParams(test_values,0.0,GLB_ImSet);
glbSetOscParams(test_values,0.0,GLB_ReSeu);
glbSetOscParams(test_values,0.0,GLB_ImSeu);

glbSetOscillationParameters(true_values);
glbSetRates();
glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_ON);

glbSetCoefficientInRule(0, 0, 0, GLB_SIG, 1.0);
glbSetCoefficientInRule(0, 0, 1, GLB_SIG, 0.0);
glbSetCoefficientInRule(0, 0, 0, GLB_BG, 1.0);

glbDefineParams(input_errors,0,2*t13stddev,0,0,4*dm21stddev,0);
glbSetOscParams(input_errors, 5.0, GLB_ReSet);
glbSetOscParams(input_errors, 5.0, GLB_ImSet);
glbSetOscParams(input_errors, 5.0, GLB_ReSeu);
glbSetOscParams(input_errors, 5.0, GLB_ImSeu);
glbSetDensityParams(input_errors,0.0,GLB_ALL);
glbSetCentralValues(true_values);
glbSetInputErrors(input_errors);

glbDefineProjection(my_projection,GLB_FIXED,GLB_FIXED,GLB_FIXED,GLB_FIXED,GLB_FIXED,GLB_FIXED);
glbSetDensityProjectionFlag(my_projection, GLB_FIXED, GLB_ALL);
glbSetProjectionFlag(my_projection, GLB_FIXED, GLB_ReSet);
glbSetProjectionFlag(my_projection, GLB_FIXED, GLB_ImSet);
glbSetProjectionFlag(my_projection, GLB_FIXED, GLB_ReSeu);
glbSetProjectionFlag(my_projection, GLB_FIXED, GLB_ImSeu);
glbSetProjection(my_projection); 

	
	
	double t12min, t12max, dt12, ddm21, dm21max, dm21min, dnsi, nsimax, nsimin, dnsi2, nsimax2, nsimin2;
	double t13min, t13max, dt13, ddm31, dm31max, dm31min;

	int M, N1, N2, N3, N4, i;
	M=10000000;
	N1=1000;
	N2=1000;
	N3=1000;
	N4=1000;
	
	t12min=log10(0.1);
	dm21min=0.5e-4;

	t12max=log10(10.);
	dm21max=1.1e-4;


//	t13min=0.05;
//	t13max=0.05;
	t13min=0.05;
	t13max=0.15;

//	dm31min=1.0e-3;
//	dm31max=1.0e-3;

	dm31min=1.0e-3;
	dm31max=4.0e-3;

	nsimin=-1.0;
	nsimax=1.0;

	nsimin2=0.0;
	nsimax2=0.0;



	dt12=(t12max-t12min)/N1;
	dt13=(t13max-t13min)/N1;
	ddm21=(dm21max-dm21min)/N2;
	ddm31=(dm31max-dm31min)/N2;
	dnsi=(nsimax-nsimin)/N3;
	dnsi2=(nsimax2-nsimin2)/N4;	

	srand (time(NULL));


for(i=0; i<M; i++)
{


	double t12test=t12min+dt12*(rand() % N1);
	double t13test=t13min+dt13*(rand() % N1);
	double dm21test=dm21min+ddm21*(rand() % N2);
	double dm31test=dm31min+ddm31*(rand() % N2);
	double nsitest=nsimin+dnsi*(rand() % N3);
	double nsitest2=nsimin2+dnsi2*(rand() % N4);

	glbSetOscParams(test_values,dm21test,GLB_DM_21);
	glbSetOscParams(test_values,dm31test,GLB_DM_31);
	glbSetOscParams(test_values,atan(sqrt(pow(10.0,t12test))),GLB_THETA_12);
	glbSetOscParams(test_values,asin(sqrt(t13test))/2.,GLB_THETA_13);
	glbSetOscParams(test_values,nsitest,GLB_ImSet);
	glbSetOscParams(test_values,nsitest2,GLB_ImSet);

	double chi2=glbChiSys(test_values,0,GLB_ALL);


outfile1<<t12test<<"\t"<<t13test<<"\t"<<dm21test<<"\t"<<dm31test<<"\t"<<nsitest<<"\t"<<nsitest2<<"\t"<<chi2<<endl;

};


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


glbFreeParams(true_values);
glbFreeParams(test_values);
glbFreeParams(input_errors);
glbFreeProjection(my_projection);


return 0;

}

// main program ends!

/* ############################################################## */
/* ############################################################## */


