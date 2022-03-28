/*==============================================================================
 HEADER: globaldefs.h		[gsm]
 ==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h


#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
//

#include "../libs/stdinc.h"
#include "../libs/numrec.h"
#include "../libs/diffeqs.h"
#include "../libs/quads.h"
#include "../libs/mathfns.h"
#include "../libs/mathutil.h"
#include "../libs/inout.h"
#include "../libs/vectmath.h"
#include "../libs/getparam.h"
#include "../libs/machines.h"
#include "../libs/strings.h"

#if !defined(global)
#  define global extern
#endif

#define IPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&(param);    \
id[nt++]=INT;}

#define RPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&param;    \
id[nt++]=DOUBLE;}

#define BPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&param;    \
id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)    \
{strcpy(tag[nt],paramtext);    \
param=(string) malloc(n);    \
addr[nt]=param;    \
id[nt++]=STRING;}
//

#include "models.h"


#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#define invH0     2997.92458   //This is H_0^{-1} in Mpc/h units
#define FOURPI2   39.4784176043574    //but see 0903.5321
#define PI2     9.8696044010893586188
#define TWOPI2     19.739208802178716
#define SIXPI2  59.21762640653615
#define INVSQRTDTWOPI 0.39894228040143267794

typedef struct {
// Background cosmology:
    real om;
    string olstr;
    real h;
// bias parameters:
    real b1;
    real b2;
    real bs;
    real c1eft;
    real c2eft;
    real s2eft;
// Differential equations evolution parameters:
	real x;
	string dxstr;
	real xstop; 
    int maxnsteps;
	string integration_method;
    real dxmin;
    real eps;

// k table
    string fnamePS;
    real kmin;
    real kmax;
    int Nk;
//
// 2pcf rsd multipoles array:
    real smin;
    real smax;
    int Ns;
    real gsm_width;
    int gsm_sizeyT;
    int gsm_NGL;
// q functions:
	int NqperLogDecade;
	int Nk_qFunctionsQuad;
// CLPT correlation functions output array:
    real rmin;
    real rmax;
    int Nr;
// Post processing parameters:
    bool postprocessing;
    string options;
//
    string paramfile;
// Modified gravity model parameters:
    string mgmodel;
    string suffixModel;
    string model_paramfile;
    int nHS;
    real fR0;
    real omegaBD;
    real screening;
// DGP:
    real eps_DGP;
    real rc_DGP;
// Quadrature parameters:
    string quadratureMethod;
    int nquadSteps;
    int ngausslegpoints;
    real epsquad;
//
} cmdline_data, *cmdline_data_ptr;




typedef struct {
	real cpuinit;
	real dx;
    int method_int;
    int quadmethod_int;

// Modified gravity model parameters:
    real beta2;
//
    char integration_method_comment[100];
    char quadraturemethod_comment[100];

	string headline0;
	string headline1;
	string headline2;
	string headline3;

    char model_comment[100];

	FILE *outlog;
    
    real ol;

    real xnow;
    real xout;
    real xoutinfo;
    real xstop;

	char mode[2];
// 
    char fnamePSPath[100];
    char logfilePath[100];
    char clptDir[100];
    char tmpDir[100];
    char fpfnamekfun[100];
    char fpfnameSPTPowerSpectrum[100];
    char fpfnameqfunctions[100];
    char fpfnameclptfunctions[100];
    char fpfnamersd[100];
    char fpfnamev12[100];
    char fpfnamesigma12_parallel[100];
    char fpfnamesigma12_perp[100];
    //~ char fpfnameTables[100];
    char fpfnameParams[100];


    char fpfnamekfun2[100];
    char fpfnameclptfunctions2[100];    

    real kf;
    real k1;
    real k2;

    real x;
    real k;
    real p;
    
    real f;
    real Dplus;
    real particles_meanPath;
    real sigma8;
} global_data, *global_data_ptr;


global global_data gd;
global cmdline_data cmd;

global real *yout;
#define NEQS3Order    10
#define NEQS2Order    8
#define NEQS1Order      2

typedef struct _pointPSTable {
    real k;
    real ps;
} pointPSTable, *pointPSTableptr;

global int nPSTable;
global pointPSTableptr PSLCDMtab;
global int nPSLogT;
global pointPSTableptr PSLCDMLogtab;

global int nPSLT;
global pointPSTableptr PSLT;
global int nPSLTLog;
global pointPSTableptr PSLTLog;

global real *kPS;
global real *pPS;
global real *pPS2;


#define kPos(x)    (((pointPSTableptr) (x))->k)
#define PS(x)    (((pointPSTableptr) (x))->ps)


typedef struct {
    real eta;
    real y1;
    real y2;
    real y3;
    real y4;
    real y5;
    real y6;
    real y7;
    real y8;
} global_D2v2, *global_D2v2_ptr;

#define etaD2(x)    (((global_D2v2_ptr) (x))->eta)
#define Dpk1D2(x)    (((global_D2v2_ptr) (x))->y1)
#define Dpk2D2(x)    (((global_D2v2_ptr) (x))->y3)
#define DA2D2(x)    (((global_D2v2_ptr) (x))->y5)
#define DB2D2(x)    (((global_D2v2_ptr) (x))->y7)
//

typedef struct {
    real eta;
    real y1;
    real y2;
    real y3;
    real y4;
    real y5;
    real y6;
    real y7;
    real y8;
    real y9;
    real y10;
} global_D3v2, *global_D3v2_ptr;


#define etaD3(x)    (((global_D3v2_ptr) (x))->eta)
#define DpkD3(x)    (((global_D3v2_ptr) (x))->y1)
#define DppD3(x)    (((global_D3v2_ptr) (x))->y3)
#define D2fD3(x)    (((global_D3v2_ptr) (x))->y5)
#define D2mfD3(x)    (((global_D3v2_ptr) (x))->y7)
#define D3symmD3(x)    (((global_D3v2_ptr) (x))->y9)
//

// GL structure
typedef struct {
    int npts;
    real x1;
    real x2;
    real *xgl;
    real *wgl;
} global_GL, *global_GL_ptr;

// STATIC problem: gcc version 11
//global_GL_ptr pGL;
global global_GL_ptr pGL;
//~ global_GL_ptr pGL;

#define nGL(x)    (((global_GL_ptr) (x))->npts)
#define x1GL(x)    (((global_GL_ptr) (x))->x1)
#define x2GL(x)    (((global_GL_ptr) (x))->x2)
#define xGL(x)    (((global_GL_ptr) (x))->xgl)
#define wGL(x)    (((global_GL_ptr) (x))->wgl)

//
// QRs structure
typedef struct {
    int eta;
    real k;
    real Q1;
    real Q2;
    real Q3;
    real Q8;
    real Q5;
    real Qs2;
    //~ real RI;
    real R1;
    real R2;
} global_QRs, *global_QRs_ptr;

#define etaQRs(x)    (((global_QRs_ptr) (x))->eta)
#define kQRs(x)    (((global_QRs_ptr) (x))->k)
#define Q1(x)    (((global_QRs_ptr) (x))->Q1)
#define Q2(x)    (((global_QRs_ptr) (x))->Q2)
#define Q3(x)    (((global_QRs_ptr) (x))->Q3)
#define Q8(x)    (((global_QRs_ptr) (x))->Q8)
#define Q5(x)    (((global_QRs_ptr) (x))->Q5)
#define Qs2(x)    (((global_QRs_ptr) (x))->Qs2)
//~ #define RI(x)    (((global_QRs_ptr) (x))->RI)
#define R1(x)    (((global_QRs_ptr) (x))->R1)
#define R2(x)    (((global_QRs_ptr) (x))->R2)
//

//
// qfunctions structure
typedef struct {
    real q;
    real UL;
    real Uloop;
    real U11;
    real U20;
    real XL;
    real Xloop;
    real YL;
    real Yloop;
    real X10;
    real Y10;
    real preV;
    real T;
    real dotVa;
    real ddotVa;
    real dotVb;
    real V10;
    real iBessel0;
    real iBessel2;
    real iBessel4;
    real nabla2xi; 
} global_qfunctions, *global_qfunctions_ptr;

#define qqfun(x)    (((global_qfunctions_ptr) (x))->q)
#define ULqfun(x)    (((global_qfunctions_ptr) (x))->UL)
#define Uloopqfun(x)    (((global_qfunctions_ptr) (x))->Uloop)
#define U11qfun(x)    (((global_qfunctions_ptr) (x))->U11)
#define U20qfun(x)    (((global_qfunctions_ptr) (x))->U20)
#define XLqfun(x)    (((global_qfunctions_ptr) (x))->XL)
#define Xloopqfun(x)    (((global_qfunctions_ptr) (x))->Xloop)
#define YLqfun(x)    (((global_qfunctions_ptr) (x))->YL)
#define Yloopqfun(x)    (((global_qfunctions_ptr) (x))->Yloop)
#define X10qfun(x)    (((global_qfunctions_ptr) (x))->X10)
#define Y10qfun(x)    (((global_qfunctions_ptr) (x))->Y10)
#define preVqfun(x)    (((global_qfunctions_ptr) (x))->preV)
#define Tqfun(x)    (((global_qfunctions_ptr) (x))->T)
#define dotVaqfun(x)    (((global_qfunctions_ptr) (x))->dotVa)
#define ddotVaqfun(x)    (((global_qfunctions_ptr) (x))->ddotVa)
#define dotVbqfun(x)    (((global_qfunctions_ptr) (x))->dotVb)
#define V10qfun(x)    (((global_qfunctions_ptr) (x))->V10)
#define iBessel0qfun(x)    (((global_qfunctions_ptr) (x))->iBessel0)
#define iBessel2qfun(x)    (((global_qfunctions_ptr) (x))->iBessel2)
#define iBessel4qfun(x)    (((global_qfunctions_ptr) (x))->iBessel4)
#define nabla2xiqfun(x)    (((global_qfunctions_ptr) (x))->nabla2xi)
//

//
// correlation functions structure
typedef struct {
    real q;
    real xi;
    real Lapxi;
} global_corrfunctions, *global_corrfunctions_ptr;

#define qcorrfun(x)    (((global_corrfunctions_ptr) (x))->q)
#define xicorrfun(x)    (((global_corrfunctions_ptr) (x))->xi)
#define Lapxicorrfun(x)    (((global_corrfunctions_ptr) (x))->Lapxi)
//

// BEGIN :: CLPT correlation auxiliary functions and structures
typedef struct {
    real r;
    real xi;
} global_zacorrfunctions, *global_zacorrfunctions_ptr;

#define rzacorrfun(x)    (((global_zacorrfunctions_ptr) (x))->r)
#define xizacorrfun(x)    (((global_zacorrfunctions_ptr) (x))->xi)

typedef struct {
    real r;
    real xiA;
    real xiW;
    real xi10;
    real xi20;
    real xi01;
    real xi02;
    real xi11;
    real nabla2xi;
    real xi001;
    real xi002;
    real xi101;
    real xi011;
    //
    real v12_00;
    real v12_10;
    real v12_20;
    real v12_01;
    real v12_11;
    real v12_001;
    real v12_101;
    real v12_eft;
	// 
	real s12perp_00;
	real s12perp_10;
	real s12perp_20;
	real s12perp_01;
	real s12perp_001;
	real s12par_00;
	real s12par_10;
	real s12par_20;
	real s12par_01;
	real s12par_001;
} global_clptcorrfunctions, *global_clptcorrfunctions_ptr;

#define rclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->r)
#define xiAclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xiA)
#define xiWclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xiW)
#define xi10clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi10)
#define xi20clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi20)
#define xi01clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi01)
#define xi02clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi02)
#define xi11clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi11)
#define nabla2xiclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->nabla2xi)
#define xi001clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi001)
#define xi002clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi002)
#define xi101clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi101)
#define xi011clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi011)

#define r_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->r)
#define v12_00_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->v12_00)
#define v12_10_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->v12_10)
#define v12_20_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->v12_20)
#define v12_01_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->v12_01)
#define v12_11_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->v12_11)
#define v12_001_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->v12_001)
#define v12_101_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->v12_101)
#define v12_eft_vfun(x)    (((global_clptcorrfunctions_ptr) (x))->v12_eft)

#define s12perp_00_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12perp_00)
#define s12perp_10_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12perp_10)
#define s12perp_20_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12perp_20)
#define s12perp_01_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12perp_01)
#define s12perp_001_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12perp_001)
#define s12par_00_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12par_00)
#define s12par_10_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12par_10)
#define s12par_20_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12par_20)
#define s12par_01_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12par_01)
#define s12par_001_sfun(x)    (((global_clptcorrfunctions_ptr) (x))->s12par_001)

typedef struct {
    real s;
    real xi_0;
    real xi_2;
    real xi_4;
} global_rsdmultipoles, *global_rsdmultipoles_ptr;

#define s_rsdmultipoles(x)    (((global_rsdmultipoles_ptr) (x))->s)
#define xi_0_rsdmultipoles(x)    (((global_rsdmultipoles_ptr) (x))->xi_0)
#define xi_2_rsdmultipoles(x)    (((global_rsdmultipoles_ptr) (x))->xi_2)
#define xi_4_rsdmultipoles(x)    (((global_rsdmultipoles_ptr) (x))->xi_4)




// END :: CLPT correlation auxiliary functions and structures



// NEW for v2

typedef struct {
real kmin;
real kmax;
int Nk;
real qmin;
real qmax;
int Nq;		
real rmin;
real rmax;
int Nr;
real smin;
real smax;
int Ns;
} global_output_lists, *global_output_lists_ptr;

global global_output_lists golists;


typedef struct {
    int eta;
    real *kT;
    real *Q1T;
    real *Q2T;
    real *Q3T;
    real *Q8T;
    real *Q5T;
    real *Qs2T;
    //~ real *RIT;
    real *R1T;
    real *R2T;
    real *pklT;
} global_kArrays, *global_kArrays_ptr;

global global_kArrays kArrays;
global global_kArrays kArraysd;




typedef struct {
real *qTab;
real *XLT;
real *YLT;
real *XloopT;
real *YloopT;
real *preVT;
real *TT;
real *X10T;
real *Y10T;
real *ULT;
real *UloopT;
real *U11T;
real *U20T;
real *dotVaT;
real *ddotVaT;
real *dotVbT;
real *V10T;
real *iBessel0T;
real *iBessel2T;
real *iBessel4T;
real *nabla2xiT;
} global_qArrays, *global_qArrays_ptr;


global global_qArrays qArrays;
global global_qArrays qArraysd;



typedef struct {
    real *rTab;
    real *xi_00T;
    real *xi_10T;
    real *xi_20T;
    real *xi_01T;
    real *xi_02T;
    real *xi_11T;
    real *xi_eftT;
    real *xi_001T;
    real *xi_002T;
    real *xi_101T;
    real *xi_011T; 
    real *xi_AT; 
    real *xi_WT; 
    real *v_00T;
    real *v_10T;
    real *v_20T;
    real *v_01T;
    real *v_11T;
    real *v_001T;
    real *v_101T;
    real *v_eftT;         
    real *spar_00T;
    real *spar_10T;
    real *spar_20T;
    real *spar_01T;
    real *spar_001T;            
    real *sperp_00T;
    real *sperp_10T;
    real *sperp_20T;
    real *sperp_01T;
    real *sperp_001T;
    real *xi_zaT;
    real *xi_LT;
} global_rArrays, *global_rArrays_ptr;


global global_rArrays rArrays;
global global_rArrays rArraysd;







global real *inout_xval;
global real *inout_yval;
global real *inout_zval;
global real *inout_wval;

global double dxsav,*xp,**yp;
global int kmax,kount;
global int nrhs;

global long idum;                // seed for random generators




#endif // ! _globaldefs_h

