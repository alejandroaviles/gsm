/*==============================================================================
 NAME: kfunctions.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com)
 ================================================================================ 
*/



#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"


local void quadrature(real ki);
local void loopQsRs(stream outstr, int imin, int imax, real dk);


#define KMIN    1.0e-20

global void compute_kfunctions(void)
{
    stream outstrQsRs, outtables;
    real dk;
    real bTime;
    real kBAOmin=0.005, kBAOmax=1.0, epsquadsave;
    int iBAOmin, iBAOmax;
    global_D2v2_ptr ptmp;
    global_D3v2_ptr ptmpR1;
 //   local lcdmCF
    real fR0save;

    bTime = second();

    
        ptmp = DsSecondOrder_func(KMIN, KMIN, KMIN);
        KA_LCDM = DA2D2(ptmp) / ( (3.0/7.0) * Dpk1D2(ptmp) * Dpk2D2(ptmp) );
        KB_LCDM = KA_LCDM;
//
        ptmpR1 = DsThirdOrder_func(0.0000001, KMIN, KMIN);
        KR1_LCDM = (21.0/5.0)*D3symmD3(ptmpR1)
        /( DpkD3(ptmpR1)*DppD3(ptmpR1)*DppD3(ptmpR1) );
  //      fprintf(stdout,"\nA_LCDM=%g,  KR1_LCDM = %g",KA_LCDM, KR1_LCDM);
        //~ cmd.fR0 = fR0save;
    

    
    fprintf(stdout,"\nk-functions:");
    fprintf(stdout," Nk=%d values from kmin=%g to kmax=%g ",
            cmd.Nk, cmd.kmin, cmd.kmax);

        
    dk = (rlog10(cmd.kmax) - rlog10(cmd.kmin))/((real)(cmd.Nk - 1));
    
    loopQsRs(outstrQsRs, 1, cmd.Nk, dk);
    
    fprintf(stdout,"...time = %g seconds",second()-bTime);    

}
#undef KMIN


local void loopQsRs(stream outstr, int imin, int imax, real dk)
{
    global_QRs qrs;
    real aTime;
    real kval;
    real ki;
    real pkl;
    //~ int counter;
    int i;


    for (i=imin; i<=imax; i++) {
        kval = rlog10(cmd.kmin) + dk*((real)(i - 1));
        ki = rpow(10.0,kval);
        qrs = k_functions_driver(gd.xstop, ki);

        pkl = psInterpolation_nr(ki, kPS, pPS, nPSLT);
 
        kArrays.kT[i-1]   = ki;                
        kArrays.Q1T[i-1]   =qrs.Q1;
        kArrays.Q2T[i-1]   = qrs.Q2;
        kArrays.Q3T[i-1]   = qrs.Q3;
        kArrays.Q5T[i-1]   = qrs.Q5;
        kArrays.Q8T[i-1]   = qrs.Q8;
        kArrays.Qs2T[i-1] = qrs.Qs2;
        kArrays.R1T[i-1]   = qrs.R1;
        kArrays.R2T[i-1]   = qrs.R2;
        kArrays.pklT[i-1]  = pkl;
    }
}





#define QROMBERG     qromo
#define KK  5





// BEGIN Qs and Rs
// kk is the inner integration moment p. 
// kk = ki * r, so usually kk is called p
global_QRs k_functions(real eta, real ki)
{
    int i, j;
    real KR1, fac;
    global_D2v2_ptr ptmp;
//    global_D3v2_ptr ptmpR1;
    real Dpkmin, Dpk;
//
    real *xxGL, *wwGL, *xGL, *wGL;
    real kmin, kmax;
    
    real Q1p, Q2p, Q3p, Q8p;
    real Q1aA, Q2aA, Q3aA, Q8aA;
    real Q1aB, Q2aB, Q3aB, Q8aB;
    real Q5p, Q5aA, Q5aB;    
    real Qs2p, Qs2aA, Qs2aB;   
    real KQ1, KQ2, KQ3, KQ8, KQs2,KQ5;
     
    real PSLA, PSLB;
    real rmin, rmax;
    real rr, deltar;
    real mumin, mumax;
    real xv, w, psl;
    real psl1;
    real Gamma2, r, kmp2, Gamma2_evR;
    int Nx, nGL;
    real ypi, dk;
    
    real R2p, R2aA, R2aB, KR2;
    real RIp, RIaA, RIaB, KRI;
    real R1aA, R1aB, R1p;
    
    real *kk, *dkk;
    //
    pointPSTableptr p;
    //
    global_QRs_ptr QRstmp;
    
    QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));
    
    kmin = kPS[1];
    kmax = kPS[nPSLT];
    if (cmd.nquadSteps==1) {
        dk = 0.;
    } else {
        dk = (rlog10(kmax) - rlog10(kmin))/((real)(cmd.nquadSteps - 1));
    }
    
    kk=dvector(1,cmd.nquadSteps);
    dkk=dvector(1,cmd.nquadSteps);
    kk[1] = rpow(10.0,rlog10(kmin));
    for (i=2; i<cmd.nquadSteps; i++) {
        ypi = rlog10(kmin) + dk*((real)(i - 1));
        kk[i] = rpow(10.0,ypi);
        dkk[i] = (kk[i]-kk[i-1]);
    }
//
// Q functions

    Q1p =0.0; Q1aA = 0.0; Q1aB = 0.0;
    Q2p =0.0; Q2aA = 0.0; Q2aB = 0.0;
    Q3p =0.0; Q3aA = 0.0; Q3aB = 0.0;
    Q5p =0.0; Q5aA = 0.0; Q5aB = 0.0;    
    Q8p =0.0; Q8aA = 0.0; Q8aB = 0.0;   
    Qs2p =0.0; Qs2aA = 0.0; Qs2aB = 0.0;

//
    PSLA = 0.0;
    rmax = kmax/ki;
    rmin = rmin/ki;
    Nx=10;
    xxGL=dvector(1,Nx);
    wwGL=dvector(1,Nx);
    
    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
        PSLB = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        mumin = MAX( -1.0, (1.0 + rsqr(rr) - rsqr(rmax)) / (2.0*rr)  );
        mumax = MIN(   1.0, (1.0  + rsqr(rr) - rsqr(rmin)) / (2.0*rr)  );
        
        if (rr>=0.5)
            mumax = 0.5/rr;
        gauleg(mumin,mumax,xxGL,wwGL,Nx);
        for (j=1; j<=Nx; j++) {
            xv = xxGL[j];
            w = wwGL[j];
            
            kmp2=1.0 + rr*rr - 2.0 * rr * xv;
            psl = psInterpolation_nr(ki * rsqrt(kmp2), kPS, pPS, nPSLT); 
            Gamma2 = KA_LCDM *(1.0-xv*xv) /  kmp2 ;  
            /*
								Gamma2 = (7/3) k . L^(2)(p,k-p)    (evaluated at p1=p, p2=k-p)
                                   (it is  the longitudinal piece of the second order LPT kernel: 
                                    see eq.(A.8) of 1909.05261)
								Gamma_2 =  1 - (p.(k-p))^2/ p^2 |k-p|^2
          */ 

            KQ1 = rsqr(rr)*Gamma2*Gamma2;
            KQ2 = (rr * xv * (1.0 - rr * xv))/kmp2 * Gamma2;
            KQ3 = rsqr(xv)*rsqr(1.0-rr*xv)/rsqr(kmp2);
            KQ5 = rr * (rr + xv - 2.0 * rr * xv*xv) / kmp2 / 2 * Gamma2;
            KQ8 = rsqr(rr)*Gamma2;
            KQs2 = - rr*rr * Gamma2  * (1. - 2.* rr*rr + 4.*rr*xv - 3.* xv*xv)   /  kmp2  ; //This is eq.D18 of 1609.02908
            //~ KQs2 = rr*rr * (xv*xv -1.0)  * (1. - 2.* rr*rr + 4.*rr*xv - 3.* xv*xv)  /  kmp2 /  kmp2 ;  The same as above
            //
            Q1aB +=   w*KQ1*psl;
            Q2aB +=   w*KQ2*psl;
            Q3aB +=   w*KQ3*psl;
            Q5aB +=   w*KQ5*psl;
            Q8aB +=   w*KQ8*psl;
            Qs2aB += w*KQs2*psl;

        }
        Q1p   += dkk[i]*(Q1aA*PSLA + Q1aB*PSLB)/2.0;
        Q2p   += dkk[i]*(Q2aA*PSLA + Q2aB*PSLB)/2.0;
        Q3p   += dkk[i]*(Q3aA*PSLA + Q3aB*PSLB)/2.0;
        Q5p   += dkk[i]*(Q5aA*PSLA + Q5aB*PSLB)/2.0;
        Q8p   += dkk[i]*(Q8aA*PSLA + Q8aB*PSLB)/2.0;
        Qs2p += dkk[i]*(Qs2aA*PSLA + Qs2aB*PSLB)/2.0;
        
        Q1aA =   Q1aB;      Q1aB = 0.0;
        Q2aA =   Q2aB;      Q2aB = 0.0;
        Q3aA =   Q3aB;      Q3aB = 0.0;
        Q5aA =   Q5aB;      Q5aB = 0.0;
        Q8aA =   Q8aB;      Q8aB = 0.0;
        Qs2aA = Qs2aB;   Qs2aB = 0.0;
        
        
        PSLA = PSLB;
    }
    
    Q1p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q2p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q3p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q5p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q8p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Qs2p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    
    free_dvector(wwGL,1,Nx);
    free_dvector(xxGL,1,Nx);

//  R functions
    R2p = 0.0; R2aA = 0.0; R2aB = 0.0;
    R1p = 0.0; R1aA = 0.0; R1aB  = 0.0; 
    //~ RIp = 0.0;  RIaA  = 0.0; RIaB   = 0.0;   

    nGL=16;
    xGL=dvector(1,nGL);
    wGL=dvector(1,nGL);
    gauleg(-1.0,1.0,xGL,wGL,nGL);
    

    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
        psl = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        for (j=1; j<=nGL(pGL); j++) {
            xv = xGL[j];
            w = wGL[j];     
            Gamma2_evR = KA_LCDM * (1.0 - rsqr(xv)); // = (7/3) (k-p).L^(2) (-p,k)
            kmp2=1.0 + rr*rr - 2.0 * rr * xv;
                  
            KR2 = ( rr*xv*(1.0-rr*xv)/kmp2 )*Gamma2_evR; 
            KR1 =  KR1_LCDM * rsqr(rr) * rsqr(1.0 -   xv*xv) / kmp2;
            //~ KRI = ( rsqr(rr)*(1.0-rsqr(xv))/abskmq )*Gamma2_evR;
            //~ KR1 = ( rsqr(rr)*(1.0-rsqr(xv))/kmp2 );  // Kernels EdS
            
                
            R2aB += w*KR2*psl;
            R1aB += w*KR1*psl;
            //~ RIaB += w*KRI*psl;
        }
        
        R2p    += dkk[i]*(R2aA + R2aB) /  (2.0*ki);
        R1p    += dkk[i]*(R1aA + R1aB ) / (2.0*ki);
        
        R2aA = R2aB;  R1aA = R1aB;
        R2aB = 0.0;  R1aB = 0.0;
      
        //~ RIp     += dkk[i]*(RIaA + RIaB  ) /  (2.0*ki);  
       //~ RIaA = RIaB; RIaB = 0.0;
    }
    fac = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    R2p *= (rpow(ki,3.0)/FOURPI2)*fac;
    //~ RIp *= (rpow(ki,3.0)/FOURPI2)*fac;
    R1p *= (rpow(ki,3.0)/FOURPI2)*fac;



    etaQRs(QRstmp) = eta;
    kQRs(QRstmp)    = ki;
    
    Q1(  QRstmp)      = Q1p;
    Q2(  QRstmp)      = Q2p;
    Q3(  QRstmp)      = Q3p;   
    Q5(  QRstmp)      = Q5p;
    Q8(  QRstmp)      = Q8p; 
    R1(  QRstmp)      = R1p;     
    R2(  QRstmp)      = R2p;
    Qs2(QRstmp)      = Qs2p;
    //~ RI(QRstmp)      = RIp;
    
    free_dvector(dkk,1,cmd.nquadSteps);
    free_dvector(kk,1,cmd.nquadSteps);
    
    return *QRstmp;
}



// END Qs and Rs

global_QRs qrs;

global global_QRs k_functions_driver(real eta, real ki)
{
    quadrature(ki);
    return qrs;
}


#define ROMO 1
#define NULLMETHOD 0
#define TRAPEZOID 2
#define TRAPEZOID3 5

local void quadrature(real ki)
{
    switch (gd.quadmethod_int) {
//
        case TRAPEZOID3:
            qrs = k_functions(gd.xstop, ki);
            break;
//
        case NULLMETHOD:
            qrs = k_functions(gd.xstop, ki);
            break;
//
        default:
            qrs = k_functions(gd.xstop, ki);
            break;
    }
}


void quadraturemethod_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"romberg") == 0) {
        *method_int = ROMO;
        strcpy(gd.quadraturemethod_comment, "romberg open quadrature method");
    }
//
    if (strcmp(method_str,"trapezoid") == 0) {
        *method_int = TRAPEZOID;
        strcpy(gd.quadraturemethod_comment, "trapezoid quadrature method");
    }
//
    if (strcmp(method_str,"trapezoid3") == 0) {
        *method_int = TRAPEZOID3;
        strcpy(gd.quadraturemethod_comment, "trapezoid3 quadrature method");
    }
//
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
        strcpy(gd.quadraturemethod_comment,
               "given null quadrature method ... running deafult (trapezoid)");
        fprintf(stdout,"\n\tintegration: default integration method (trapezoid)...\n");
    }
//
    if (*method_int == -1) {
        *method_int = TRAPEZOID;
        strcpy(gd.quadraturemethod_comment,
               "Unknown quadrature method ... running deafult (trapezoid)");
        fprintf(stdout,"\n\tquadrature: Unknown method... %s ",cmd.quadratureMethod);
        fprintf(stdout,
                "\n\trunning default quadrature method (trapezoid)...\n");
    }
}

#undef ROMO
#undef TRAPEZOID
#undef TRAPEZOID3
#undef NULLMETHOD



#undef KK
#undef QROMBERG












