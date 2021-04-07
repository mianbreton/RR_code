/* ============================================================ *
 * Michel-Andrès Breton 2020				        *
 * 							        *
 * Code for analytical random-random pair counts                *
 * ============================================================ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "cuba.h"

// Error on the integrals
#define EPSABS 0 // absolute error

// Cuba general settings
#define NDIM 3
#define NVEC 1
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 500000 // Change that for more accuracy (default 50,000)

// Vegas specific settings
#define NSTART 10000 // nb of integrand evaluation per interval to start with
#define NINCREASE 50000 // Increase of nb of integrand evaluation
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

// Suave specific settings
#define NNEW 10000 // Nb of new integrand evaluation in each subdivision
#define NMIN 2
#define FLATNESS 25.

// Divonne specific settings
#define KEY1 5000 // 47 by default
#define KEY2 1 // 1 by default
#define KEY3 1 // 1 by default
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

// Cuhre specific settings
#define KEY 0 // Integration rule (11 for 3D integration)

// Global variables and functions
double rmin1, rmax1, rmin2, rmax2, smin, smax, W0_1, W0_2, W1W2, n0, epsrel;
int ns, nmu, ncomp_mu, ncomp_smu, ncomp_legendre;
char dndrfile1[BUFSIZ], dndrfile2[BUFSIZ], wformat[BUFSIZ], wfile[BUFSIZ], angle[BUFSIZ], integration[BUFSIZ], eps_rel[BUFSIZ], output_base[BUFSIZ];
const size_t size_cquad = 1000; // workspace size for cquad
const size_t size_qag = 500; // workspace size for qag
gsl_interp_accel *acc[5]; 				
gsl_spline *spline[5]; //Use for interpolation
struct my_f_params2 {double a; double b;}; 
struct my_f_params3 {double a; double b; double c;}; 
struct my_f_params4 {double a; double b; double c; double d;}; 
struct my_f_params6 {double a; double b; double c; double d; double e; double f;}; 

/// INTEGRATION FUNCTIONS ////

// GSL
double int_qag (double func(double, void*), void *p, double xmin, double xmax)
{
    double result, error;
    gsl_function F;
    F.function = func;
    F.params = p;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(size_qag);
    gsl_integration_qag(&F, xmin, xmax, EPSABS, epsrel, size_qag, 6, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}

double int_cquad (double func(double, void*), void *p, double xmin, double xmax)
{
    double result, error;
    gsl_function F;
    F.function = func;
    F.params = p;

    gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(size_cquad);
    gsl_integration_cquad(&F, xmin, xmax, EPSABS, epsrel, w, &result, &error, NULL);
    gsl_integration_cquad_workspace_free(w);
    return result;
}

// CUBA
void int_vegas(int size_comp, cubareal* result, integrand_t func, void *p)
{
    int neval, fail;
    cubareal error[size_comp], prob[size_comp];

    Vegas(NDIM, size_comp, func, p, NVEC,
        epsrel, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        GRIDNO, STATEFILE, SPIN,
        &neval, &fail, result, error, prob);
}

void int_suave(int size_comp, cubareal* result, integrand_t func, void *p)
{
    int nregions, neval, fail;
    cubareal error[size_comp], prob[size_comp];

    Suave(NDIM, size_comp, func, p, NVEC,
        epsrel, EPSABS, VERBOSE | LAST, SEED,
        MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
        STATEFILE, SPIN,
        &nregions, &neval, &fail, result, error, prob);
}

void int_divonne(int size_comp, cubareal* result, integrand_t func, void *p)
{
    int nregions, neval, fail;
    cubareal error[size_comp], prob[size_comp];

    Divonne(NDIM, size_comp, func, p, NVEC,
        epsrel, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
        BORDER, MAXCHISQ, MINDEVIATION,
        NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
        STATEFILE, SPIN,
        &nregions, &neval, &fail, result, error, prob);
}

void int_cuhre(int size_comp, cubareal* result, integrand_t func, void *p)
{
    int nregions, neval, fail;
    cubareal error[size_comp], prob[size_comp];

    Cuhre(NDIM, size_comp, func, p, NVEC,
        epsrel, EPSABS, VERBOSE | LAST,
        MINEVAL, MAXEVAL, KEY,
        STATEFILE, SPIN,
        &nregions, &neval, &fail, result, error, prob);
}

/// INTERPOLATIONS ////

double wtheta_iphi(double phi)
{
    if (phi >= 0 && phi <= M_PI)  return gsl_spline_eval (spline[0], phi, acc[0]);
    else return 0;
}

double n1_ir(double r)
{
    //if (r > rmin && r <= rmax)  return n0; // uniform   
    if (r > rmin1 && r <= rmax1)  return gsl_spline_eval (spline[1], r, acc[1]);
    else return 0;
}

double dndr1_ir(double r)
{
    if (r > rmin1 && r <= rmax1)  return gsl_spline_eval (spline[2], r, acc[2]);
    else return 0;
}

double n2_ir(double r)
{
    if (r > rmin2 && r <= rmax2)  return gsl_spline_eval (spline[3], r, acc[3]);
    else return 0;
}

double dndr2_ir(double r)
{
    if (r > rmin2 && r <= rmax2)  return gsl_spline_eval (spline[4], r, acc[4]);
    else return 0;
}

//// USEFUL FUNCTIONS ////

double r_2(double r1, double s, double mu)
{
    return r1*sqrt(1. + 2*mu*s/r1 + s*s/(r1*r1));
}

double r(double r1, double s, double mu)
{
    return r1*sqrt(1. + mu*s/r1 + s*s/(4*r1*r1));
}

double Phi(double r_1, double r_2, double s)
{
    return acos((r_1*r_1+r_2*r_2-s*s)/(2*r_1*r_2));
}

double mu_MP_to_mu_EP(double r1, double s, double mu){
    if (mu == 1)
	return 1;
    else if (mu == -1)
	return -1;
    else
	return (-s+s*mu*mu+mu*sqrt(s*s*mu*mu-s*s+4*r1*r1))/(2*r1);
}

double mu_EP_to_mu_MP(double r1, double s, double mu){
    if (mu == 1)
	return 1;
    else if (mu == -1)
	return -1;
    else
	return mu*r1/r(r1, s, mu) + 0.5*s/r(r1, s, mu);
}

double legendre1(double x)
{
    return x;
}
double legendre2(double x)
{
    return 0.5*(3.*x*x-1.);
}
double legendre3(double x)
{
    return 0.5*(5.*x*x*x-3.*x);
}
double legendre4(double x)
{
    return 0.125*(35.*x*x*x*x-30.*x*x+3.);
}
double legendre5(double x)
{
    return (63.*x*x*x*x*x-70.*x*x*x+15.*x)*0.125;
}
double legendre6(double x)
{
    return (231.0*x*x*x*x*x*x-315.0*x*x*x*x+105.0*x*x-5.0)*0.0625;
}
double legendre7(double x)
{
    return (429.*x*x*x*x*x*x*x-693.*x*x*x*x*x+315.*x*x*x-35.*x)*0.0625;
}
double legendre8(double x)
{
    return (6435.0*x*x*x*x*x*x*x*x-12012.0*x*x*x*x*x*x+6930.0*x*x*x*x-1260.0*x*x+35.0)/128.0;
}

//// INTEGRANDS ////

// CUBA
static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *p)
{
    // Get args
    double r1 = xx[0];
    double s = xx[1];
    double mu = xx[2];
    
    // Get limits
    struct my_f_params2 * params = (struct my_f_params2 *)p;
    double ismin = (params->a);
    double ismax = (params->b);
    
    // Rescale integrand for all bounds from 0 to 1
    // Rescale r1
    double r1p = rmin1 + (rmax1 - rmin1)*r1;
    double n1 = n1_ir(r1p);
    double sp = ismin + (ismax - ismin)*s;
    
    // Now loop over (mu) bins
    double dmu = 1./nmu;
    for (int j = 0; j < nmu; j++) {
	/// Positive pairs
	double imu_min = j*dmu;
	double imu_max = imu_min + dmu;
	
	// Rescale s
        // Correction for end-point to mid-point bounds if needed
        double mu_min = 0;
        double mu_max = 0;
        if (strcmp(angle, "mid") == 0) {
            mu_min = mu_MP_to_mu_EP(r1p, sp, imu_min);
            mu_max = mu_MP_to_mu_EP(r1p, sp, imu_max);
        } else if(strcmp(angle, "end") == 0){
	    mu_min = imu_min;
	    mu_max = imu_max;
        } else {
	    printf("In parameter file, angle needs: 'mid' or 'end'\n");
        }
	
	// Rescale mu
        double mup = mu_min + (mu_max - mu_min)*mu;
	
	// Compute r2
        double r2 = r_2(r1p, sp, mup);
	
	// Compute integrand
        ff[2*j] = (rmax1-rmin1)*(ismax-ismin)*(mu_max-mu_min)*r1p*r1p*sp*sp*n1*n2_ir(r2)*wtheta_iphi(Phi(r1p, r2, sp));

	/// Negative pairs
	imu_min = -imu_min;
	imu_max = -imu_max;
	
        // Correction for end-point to mid-point bounds if needed
        if (strcmp(angle, "mid") == 0) {
            mu_min = mu_MP_to_mu_EP(r1p, sp, imu_min);
            mu_max = mu_MP_to_mu_EP(r1p, sp, imu_max);
        } else if(strcmp(angle, "end") == 0){
	    mu_min = imu_min;
	    mu_max = imu_max;
        } else {
	    printf("In parameter file, angle needs: 'mid' or 'end'\n");
        }
	
	// Rescale mu
        mup = mu_min + (mu_max - mu_min)*mu;
	
	// Compute r2
        r2 = r_2(r1p, sp, mup);
	
	// Compute integrand
        ff[2*j + 1] = (rmax1-rmin1)*(ismax-ismin)*(mu_max-mu_min)*r1p*r1p*sp*sp*n1*n2_ir(r2)*wtheta_iphi(Phi(r1p, r2, sp));
    }
    return 0;
}

// Legendre
static int Integrand_legendre(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *p)
{
    // Get args
    double r = xx[0];
    double s = xx[1];
    double mu = xx[2];
    
    // Get limits
    // Rescale integrand for all bounds from 0 to 1
    // Rescale r1
    double r1p = rmin1 + (rmax1 - rmin1)*r;
    
    // Now loop over (s,mu) bins
    double ds = (smax-smin)/ns;
    for (int i = 0; i < ns; i++) {
	// Get (s, mu) bins
	double ismin = smin + i*ds;
	double ismax = ismin + ds;
	
	// Rescale s
        double sp = ismin + (ismax - ismin)*s;
	
        // Correction for mid-point to end-point legendre if needed
        double mup = -1 + 2*mu;
        double mup_star = 0;
        if (strcmp(angle, "mid") == 0) {
            mup_star = mu_EP_to_mu_MP(r1p, sp, mup);
        } else if(strcmp(angle, "end") == 0){
	    mup_star = mup;
        } else {
	    printf("In parameter file, angle needs: 'mid' or 'end'\n");
        }
        double r2 = r_2(r1p, sp, mup);
	double factor = (rmax1-rmin1)*(ismax-ismin)*2*r1p*r1p*sp*sp*n1_ir(r1p)*n2_ir(r2)*wtheta_iphi(Phi(r1p, r2, sp));
        ff[9*i+0] = 0.5*factor;
        ff[9*i+1] = 1.5*factor*legendre1(mup_star);
        ff[9*i+2] = 2.5*factor*legendre2(mup_star);
        ff[9*i+3] = 3.5*factor*legendre3(mup_star);
        ff[9*i+4] = 4.5*factor*legendre4(mup_star);
        ff[9*i+5] = 5.5*factor*legendre5(mup_star);
        ff[9*i+6] = 6.5*factor*legendre6(mup_star);
        ff[9*i+7] = 7.5*factor*legendre7(mup_star);
        ff[9*i+8] = 8.5*factor*legendre8(mup_star);
    }
    return 0;
}

// GSL
double fM(double mu, void *p)
{
    struct my_f_params2 * params = (struct my_f_params2 *)p;
    double r1 = (params->a);
    double s = (params->b);
    double r2 = r_2(r1, s, mu);
    return n1_ir(r1)*n2_ir(r2)*wtheta_iphi(Phi(r1, r2, s));
}

double M(double r, double s, double mu_min_tmp, double mu_max_tmp)
{
    struct my_f_params2 p = {r, s}; 
    if (nmu > 0){
        // Correction for end-point to mid-point bounds if needed
        double mu_min = 0;
        double mu_max = 0;
        if (strcmp(angle, "mid") == 0) {
            mu_min = mu_MP_to_mu_EP(r, s, mu_min_tmp);
            mu_max = mu_MP_to_mu_EP(r, s, mu_max_tmp);
        } else if(strcmp(angle, "end") == 0){
	    mu_min = mu_min_tmp;
	    mu_max = mu_max_tmp;
        } else {
	    printf("In parameter file, angle needs: 'mid' or 'end'\n");
        }
	//printf("M\n");
        return int_cquad(fM, &p, mu_min, mu_max);
        //return int_qag(fM, &p, mu_min, mu_max);
    } else {
        return 0; //int_cquad(fM_legendre, &p, -1, 1);
        //return int_qag(fM_legendre, &p, -1, 1);
    }
}

double fL(double s, void *p)
{
    struct my_f_params3 * params = (struct my_f_params3 *)p;
    double r = (params->a);
    double mu_min = (params->b);
    double mu_max = (params->c);
    return s*s*M(r, s, mu_min, mu_max);
}

double L(double r, double smin, double smax, double mu_min, double mu_max)
{
    struct my_f_params3 p = {r, mu_min, mu_max}; 
    //printf("L %f\n", r);
    return int_cquad(fL, &p, smin, smax);
    //return int_qag(fL, &p, smin, smax);
}

double fK(double r, void *p)
{
    struct my_f_params4 * params = (struct my_f_params4 *)p;
    double smin = (params->a);
    double smax = (params->b);
    double mu_min = (params->c);
    double mu_max = (params->d);
    return r*r*L(r, smin, smax, mu_min, mu_max);
}

double K(double smin, double smax, double mu_min, double mu_max)
{
    struct my_f_params4 p = {smin, smax, mu_min, mu_max}; 
    //printf("K\n");
    return int_cquad(fK, &p, rmin1, rmax1);
    //return int_qag(fK, &p, rmin1, rmax1);
}

double RR_gsl(double smin, double smax, double mu_min, double mu_max)
{
    return 8*M_PI*M_PI*K(smin, smax, mu_min, mu_max);
}

void RR_vegas(cubareal* result, double ismin, double ismax)
{
    struct my_f_params2 p = {ismin, ismax};   
    int_vegas(ncomp_mu, result, Integrand, &p);
    for (int i = 0; i < ncomp_mu; i++)
	result[i] *= 8*M_PI*M_PI;
}

void RR_suave(cubareal* result, double ismin, double ismax)
{
    struct my_f_params2 p = {ismin, ismax};   
    int_suave(ncomp_mu, result, Integrand, &p);
    for (int i = 0; i < ncomp_mu; i++)
	result[i] *= 8*M_PI*M_PI;
}

void RR_divonne(cubareal* result, double ismin, double ismax)
{
    struct my_f_params2 p = {ismin, ismax};   
    int_divonne(ncomp_mu, result, Integrand, &p);
    for (int i = 0; i < ncomp_mu; i++)
	result[i] *= 8*M_PI*M_PI;
}

void RR_cuhre(cubareal* result, double ismin, double ismax)
{
    struct my_f_params2 p = {ismin, ismax};   
    int_cuhre(ncomp_mu, result, Integrand, &p);
    for(int i = 0; i < ncomp_mu; i++)
	result[i] *= 8*M_PI*M_PI;
}

// Legendre 
void RR_vegas_legendre(cubareal* result)
{
    int_vegas(ncomp_legendre, result, Integrand_legendre, NULL);
    for(int i = 0; i < ncomp_legendre; i++)
	result[i] *= 8*M_PI*M_PI;
}

void RR_suave_legendre(cubareal* result)
{
    int_suave(ncomp_legendre, result, Integrand_legendre, NULL);
    for(int i = 0; i < ncomp_legendre; i++)
	result[i] *= 8*M_PI*M_PI;
}

void RR_divonne_legendre(cubareal* result)
{
    int_divonne(ncomp_legendre, result, Integrand_legendre, NULL);
    for(int i = 0; i < ncomp_legendre; i++)
	result[i] *= 8*M_PI*M_PI;
}

void RR_cuhre_legendre(cubareal* result)
{
    int_cuhre(ncomp_legendre, result, Integrand_legendre, NULL);
    for(int i = 0; i < ncomp_legendre; i++)
	result[i] *= 8*M_PI*M_PI;
}

// Counts
double fN1(double r, void *p)
{
    return dndr1_ir(r);
    //return n0;
}

double fN2(double r, void *p)
{
    return dndr2_ir(r);
    //return n0;
}

double N1(double rmin, double rmax)
{
    return int_cquad(fN1, NULL, rmin, rmax);
    //return int_qag(fN1, NULL, rmin, rmax);
}

double N2(double rmin, double rmax)
{
    return int_cquad(fN2, NULL, rmin, rmax);
    //return int_qag(fN2, NULL, rmin, rmax);
}

// Get ASCII file length
int get_lines(char* file)
{
    FILE *f;
    int count = 0;
    char c;
    
    // Open the file 
    f = fopen(file, "r"); 
  
    // Check if file exists 
    if (f == NULL) 
    { 
        printf("Could not open file %s", file); 
        return 0; 
    } 
  
    // Extract characters from file and store in character c 
    for (c = getc(f); c != EOF; c = getc(f)) 
        if (c == '\n') // Increment count if this character is newline 
            count = count + 1; 
  
    // Close the file 
    fclose(f); 
    printf("The file %s has %d lines\n", file, count); 
    return count;
}

//// INPUT FILES AND CONFIGURATION READERS ////

void get_w_theta(char* file)
{
    FILE *f;
    int i;
    int val;
	
    f = fopen(file, "r");
    if (strcmp(wformat, "polspice") == 0) {
	printf("Angular correlation function has PolSpice format\n");
        val = get_lines(file)-1;
        double phi[val+1], wtheta[val+1], phir[val+1], wthetar[val+1];
        double d1;
	
        // Skip header
        fscanf(f,"%*[^\n]\n");

        // Read angular correlation function
        for(i=0; i < val; i++) fscanf(f, "%lf %lf %lf\n", &phi[i], &d1, &wtheta[i]);
	
        // Need to put in right order, and also w(0)
        phir[0] = 0.;
        wthetar[0] = W1W2;
        for(i=0; i < val; i++) {
            phir[1+i] = phi[val-1-i];
            wthetar[1+i] = wtheta[val-1-i];
        }
        acc[0] = gsl_interp_accel_alloc ();
        spline[0] = gsl_spline_alloc(gsl_interp_cspline, val);
        gsl_spline_init (spline[0], phir, wthetar, val);
    } else if (strcmp(wformat, "standard") == 0) {
	printf("Angular correlation function has standard format\n");
        val = get_lines(file);
        double phi[val+1], wtheta[val+1];
	
        for(i=0; i < val; i++) fscanf(f, "%lf %lf\n", &phi[i+1], &wtheta[i+1]);
        phi[0] = 0.;
        wtheta[0] = W1W2;

        acc[0] = gsl_interp_accel_alloc ();
        spline[0] = gsl_spline_alloc(gsl_interp_cspline, val);
        gsl_spline_init (spline[0], phi, wtheta, val);
    } else {
	printf("WARNING! In the parameter file, wformat must be 'standard' or 'polspice'\n");
    }
}

void get_dndr(char* file1, char* file2)
{
    FILE *f;
    
    // dndr(r)_1
    int val = get_lines(file1);
    double chi1[val], dndr1[val], n1[val];
    f = fopen(file1, "r");
    for(int i=0; i < val; i++){
	fscanf(f, "%lf %lf\n", &chi1[i], &dndr1[i]);
    }
    rmin1 = chi1[0];
    rmax1 = chi1[val-1];
    
    // Compute n(r)_1
    for(int i=0; i < val; i++){
	n1[i] = dndr1[i]/(4*M_PI*chi1[i]*chi1[i]*W0_1);
    }
    acc[1] = gsl_interp_accel_alloc ();
    spline[1] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[1], chi1, n1, val);

    acc[2] = gsl_interp_accel_alloc ();
    spline[2] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[2], chi1, dndr1, val);

    fclose(f);

    // dndr(r)_2
    val = get_lines(file2);
    double chi2[val], dndr2[val], n2[val];
    f = fopen(file2, "r");
    for(int i=0; i < val; i++){
	fscanf(f, "%lf %lf\n", &chi2[i], &dndr2[i]);
    }

    // Compute n(r)_2
    for(int i=0; i < val; i++){
	n2[i] = dndr2[i]/(4*M_PI*chi2[i]*chi2[i]*W0_2);
    }
    rmin2 = chi2[0];
    rmax2 = chi2[val-1];

    acc[3] = gsl_interp_accel_alloc ();
    spline[3] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[3], chi2, n2, val);

    acc[4] = gsl_interp_accel_alloc ();
    spline[4] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[4], chi2, dndr2, val);
}

void get_param(FILE* par)
{
    // Check parameter file
    if (par!=0) EXIT_SUCCESS;											
    else {
	printf("Wrong parameter file\n"); 
    }

    // Read parameter file
    fscanf(par,"%*s %s\n", dndrfile1);
    fscanf(par,"%*s %s\n", dndrfile2);
    fscanf(par,"%*s %s\n", wformat);
    fscanf(par,"%*s %s\n", wfile);
    fscanf(par,"%*s %lf\n", &W0_1);
    fscanf(par,"%*s %lf\n", &W0_2);
    fscanf(par,"%*s %lf\n", &W1W2);
    fscanf(par,"%*s %lf\n", &smin);
    fscanf(par,"%*s %lf\n", &smax);
    fscanf(par,"%*s %d\n", &ns);
    fscanf(par,"%*s %d\n", &nmu);
    fscanf(par,"%*s %s\n", angle);
    fscanf(par,"%*s %s\n", integration);
    fscanf(par,"%*s %s\n", eps_rel);
    fscanf(par,"%*s %s\n", output_base);

    epsrel = (double)strtod(eps_rel,NULL);
    ncomp_mu = nmu*2; // times 2 because of positive and negative mu
    ncomp_smu = ns*ncomp_mu;
    ncomp_legendre = ns*9; // Compute 9 multipoles

    if (nmu>0 && ncomp_smu>1024) {
      printf("Error! CUBA maximum internal parallelization reached. Reduce the number of bins to be computed.\n");
      exit(1);
    }

    if (nmu<=0 && ncomp_legendre>1024) {
      printf("Error! CUBA maximum internal parallelization reached. Reduce the number of bins to be computed.\n");
      exit(1);
    }
    
    printf("smin, smax = %f - %f with %d bins, and %d bins in mu\n", smin, smax, ns, nmu);
    printf("Angle definition : %s-point\n", angle);
    printf("Integration : %s, with a relative precision of %s\n", integration, eps_rel);
}

int main(int argc, char *argv[])
{
    // Read parameter file
    FILE *par;
    par = fopen(argv[1],"r");  
    get_param(par);
    
    // Get mask correlation
    get_w_theta(wfile);
    
    // Get n(r), the n(r) must already be in the good range !!!!
    get_dndr(dndrfile1, dndrfile2);
    
    // Make Some calculations
    double Nr = 1; // Normalise the RR
    double V0 = 4./3.*M_PI*(rmax1*rmax1*rmax1 - rmin1*rmin1*rmin1);
    n0 = Nr/(V0*W0_1); // mean density
    //printf("Mean density : n0 = %.6le between chi = [%f-%f]\n", n0, rmin1, rmax1);
    double Ntot1 = N1(rmin1, rmax1);
    double Ntot2 = N2(rmin2, rmax2);
    printf("n1(z) gives N = %f, while we want N = %f. Need to apply correction factor\n", Ntot1, Nr);
    printf("n2(z) gives N = %f, while we want N = %f. Need to apply correction factor\n", Ntot2, Nr);
    double f1 = Nr/Ntot1;
    double f2 = Nr/Ntot1;
    double ff = f1*f2;
    double ds = (smax-smin)/ns;
    double dmu = 1./nmu;
    
    // Open file to write
    FILE *f; 
    char output[BUFSIZ];
    if (nmu > 0) {
        sprintf(output,"%s_%s_%s_%s.txt", output_base, angle, integration, eps_rel);
    } else{
        sprintf(output,"%s_%s_%s_%s_legendre.txt", output_base, angle, integration, eps_rel);
    }
    f = fopen(output, "w+");
    
    // Compute RR and write output
    if(nmu > 0){
	double result[ncomp_smu];
	double result_tmp[ncomp_mu];
    	for (int i = 0.; i < ns; i++){
	    double ismin = smin + i*ds;
	    double ismax = ismin + ds;
	    if(strcmp(integration, "gsl") == 0){
		for(int j = 0; j < nmu; j++){
		    int k = i*nmu + j;
	            double imu_min = j*dmu;
	            double imu_max = imu_min + dmu;
	            result[2*k] = RR_gsl(ismin, ismax, imu_min, imu_max);
	            result[2*k+1] = RR_gsl(ismin, ismax, -imu_min, -imu_max);
	            printf("Progress : %.1f%%\r", (double)(k)/(ns*nmu)*100);
	            fflush(stdout);
	        }
            } else if(strcmp(integration, "vegas") == 0){
	        RR_vegas(result_tmp, ismin, ismax);
		for(int j = 0; j < nmu; j++){
		    result[i*2*nmu + 2*j] = result_tmp[2*j];
		    result[i*2*nmu + 2*j + 1] = result_tmp[2*j+1];
		}
	        printf("Progress : %.1f%%\r", (double)(i)/ns*100);
	        fflush(stdout);
            } else if(strcmp(integration, "suave") == 0){
	        RR_suave(result_tmp, ismin, ismax);
		for(int j = 0; j < nmu; j++){
		    result[i*2*nmu + 2*j] = result_tmp[2*j];
		    result[i*2*nmu + 2*j + 1] = result_tmp[2*j+1];
		}
	        printf("Progress : %.1f%%\r", (double)(i)/ns*100);
	        fflush(stdout);
            } else if(strcmp(integration, "divonne") == 0){
	        RR_divonne(result_tmp, ismin, ismax);
		for(int j = 0; j < nmu; j++){
		    result[i*2*nmu + 2*j] = result_tmp[2*j];
		    result[i*2*nmu + 2*j + 1] = result_tmp[2*j+1];
		}
	        printf("Progress : %.1f%%\r", (double)(i)/ns*100);
	        fflush(stdout);
            } else if(strcmp(integration, "cuhre") == 0){
	        RR_cuhre(result_tmp, ismin, ismax);
		for(int j = 0; j < nmu; j++){
		    result[i*2*nmu + 2*j] = result_tmp[2*j];
		    result[i*2*nmu + 2*j + 1] = result_tmp[2*j+1];
		}
	        printf("Progress : %.1f%%\r", (double)(i)/ns*100);
	        fflush(stdout);
	    } else {
	        printf("In parameter file, integration needs: gsl, vegas, suave, divonne or cuhre.\n");
	    }
	}
	// Finalize
    	for (int k = 0.; k < ns*nmu; k++) {
	    int i = k/nmu;
	    int j = k%nmu;
	    double ismin = smin + i*ds;
	    double imu_min = j*dmu;
	    fprintf(f, "%.13le %.13le %.13le %.13le %.13le %.13le\n", ismin, ismin+ds, imu_min, imu_min + dmu, -ff*result[2*k+1], ff*result[2*k]);	    
	}
    } else {
        // Compute multipoles
	double result[ncomp_legendre];
	
        if(strcmp(integration, "gsl") == 0){
	    printf("GSL NOT IMPLEMENTED YET FOR LEGENDRE\n");
        } else if(strcmp(integration, "vegas") == 0){
	    RR_vegas_legendre(result);
        } else if(strcmp(integration, "suave") == 0){
	    RR_suave_legendre(result);
        } else if(strcmp(integration, "divonne") == 0){
	    RR_divonne_legendre(result);
        } else if(strcmp(integration, "cuhre") == 0){
	    RR_cuhre_legendre(result);
	} else {
	    printf("In parameter file, integration needs: gsl, vegas, suave, divonne or cuhre.\n");
	}
	// Finalize
    	for (int i = 0.; i < ns; i++) {
	    double ismin = smin + i*ds;
	    fprintf(f, "%.13le %.13le %.13le %.13le %.13le %.13le %.13le %.13le %.13le %.13le %.13le\n", ismin, ismin+ds, ff*result[9*i], ff*result[9*i+1], ff*result[9*i+2], ff*result[9*i+3], ff*result[9*i+4], ff*result[9*i+5], ff*result[9*i+6], ff*result[9*i+7], ff*result[9*i+8]);
	}
    }
    printf("\n");
    printf("Run completed !\n");
    
    return 0;
}
