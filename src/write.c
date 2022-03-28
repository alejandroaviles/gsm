/*==============================================================================
 NAME: main.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com), ...
 *  (other people who collaborated: Mario A. Rodriguez-Meza ...)
 ================================================================================ 
 Use: ./gsm -help
 References:  arXiv:1306.1804, arXiv: 1609.02908, arXiv:1909.05261
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

#define FMTQFUNCTIONSDAT    \
"%e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e\n"
//~ local real sigma8(void)
//~ local real sigma8_int(void)
//~ local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local void print_kfunctions(void);
//~ local void print_kfunctions(char[]);
void print_rfunctions(void);


global void write(void){
	
	stream outstr;
	
	    fprintf(stdout,"\nz = %g\n",cmd.xstop);
	    fprintf(stdout,"Dplus = %g\n",gd.Dplus);
	    fprintf(stdout,"f = %g\n",gd.f);    
	    fprintf(stdout,"2*sigma_psi = %g  Mpc/h (Lagrangian particles mean displacement)\n",
							gd.particles_meanPath); 
							
										
	print_kfunctions();
	//~ print_kfunctions(gd.fpfnamekfun2);
	print_rfunctions();

/*
 * ***********************
 * Print qfunctions
 * **********************
*/
   outstr = stropen(gd.fpfnameqfunctions,"w!");
       //~ fprintf(stdout,"\nstropen %s\n",gd.fpfnameqfunctions2);

    fprintf(outstr,"# InputPklFile=%s, redshift z=%g\n",
                cmd.fnamePS, cmd.xstop);                     
    fprintf(outstr,"# Precision: q-functions: NqperLogDecade=%d, Nk_qFunctionsQuad=%d\n",
                cmd.NqperLogDecade,cmd.Nk_qFunctionsQuad); 
    fprintf(outstr,"%1s%5s%12s%11s%13s%11s%9s%11s%12s%11s%11s%13s%8s%11s%11s%13s%13s%12s%20s%12s%12s%12s",
            "#","1.q","2.XL","3.YL","4.Xloop",
            "5.Yloop","6.preV","7.T","8.X10",
            "9.Y10","10.UL","11.Uloop","12.U11",
            "13.U20","14.dotVa","15.ddotVa","16.dotVb","17.V10","18.xiL(xiell0)","19.xiell2",
            "20.xiell4","22.Lapxi\n");

   for (int i=0; i<golists.Nq; i++) {
	   fprintf(outstr,FMTQFUNCTIONSDAT,
	   	qArrays.qTab[i],
		qArrays.XLT[i],
		qArrays.YLT[i],
		qArrays.XloopT[i],
		qArrays.YloopT[i],
		qArrays.preVT[i],
		qArrays.TT[i],
		qArrays.X10T[i],
		qArrays.Y10T[i],
		qArrays.ULT[i],
		qArrays.UloopT[i],
		qArrays.U11T[i],
		qArrays.U20T[i],
		qArrays.dotVaT[i],
		qArrays.ddotVaT[i],
		qArrays.dotVbT[i], 
		qArrays.V10T[i],
		qArrays.iBessel0T[i],
		qArrays.iBessel2T[i],
		qArrays.iBessel4T[i],
		qArrays.nabla2xiT[i]
       );
   }

    fclose(outstr);




	
}
#undef FMTQFUNCTIONSDAT

#define FMTKOUT    \
"%e %e %e %e %e %e %e %e \
%e %e\n"

void print_kfunctions(void){

	stream outstr;

	outstr = stropen(gd.fpfnamekfun,"w!");

    fprintf(outstr,"# InputPklFile=%s, redshift z=%g\n",
                cmd.fnamePS, cmd.xstop);                     
    fprintf(outstr,"# Precision: k-functions: nquadStetps=%d\n",
                cmd.nquadSteps); 
	fprintf(outstr,"%1s%8s%13s%13s%13s%13s%13s%14s%13s%13s%17s",
            "#","1.k","2.Q1","3.Q2","4.Q3",
            "5.Q5","6.Q8","7.Qs2",
            "8.R1","9.R2","10.pklinear\n");


   for (int i=0; i<golists.Nk; i++) {
	    fprintf(outstr,FMTKOUT ,
			kArrays.kT[i],
			kArrays.Q1T[i],
			kArrays.Q2T[i],
			kArrays.Q3T[i],
			kArrays.Q5T[i],
			kArrays.Q8T[i],
			kArrays.Qs2T[i],
			kArrays.R1T[i],
			kArrays.R2T[i],
			kArrays.pklT[i]	
		);	
	};
    fclose(outstr);	
	
};

#define CLPT_OUT_FMT    \
"%e %e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e %e %e %e\n"
 
 
void print_rfunctions(void){
	
   stream outstr;
   
   outstr = stropen(gd.fpfnameclptfunctions,"w!");

   fprintf(outstr,"# InputPklFile=%s, redshift z=%g, OmegaM=%g, h=%g, Linear growth rate  f=%g\n",
                cmd.fnamePS, cmd.xstop, cmd.om, cmd.h, gd.f);       
               
    fprintf(outstr,"# Precision: q-functions: NqperLogDecade=%d, Nk_qFunctionsQuad=%d\n",
                cmd.NqperLogDecade,cmd.Nk_qFunctionsQuad); 


    fprintf(outstr,"%3s%12s%16s%16s%16s%18s%18s%18s%28s%20s%20s%20s%20s%20s%16s%16s%16s \
    %19s%19s%19s%19s%19s%19s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s",
            "#","1.r[Mpc/h]","2.xiCLPT[1]","3.xi_10[b1]","4.xi_20[b1^2]","5.xi_01[b2] ",
            "6.xi_02[b2^2]","7.xi_11[b1*b2]","8.xi_eft[alpha1eft]",
            "9.xi_001[bs]","10.xi_002[bs^2] ","11.xi_101[b1*bs]","12.xi_011[b2*bs]",
            "13.xi_L(linearCF)  ","14.xi_ZA","15.xi_A", "16.xi_W","17.v12_00[1]","18.v12_10[b1]",
            "19.v12_20[b1^2]","20.v12_01[b2]","21.v12_11[b1*b2]","22.v12_001[bs]","23.v12_101[b1*bs]","24.v12_eft[alpha2eft]",
            "25.s12par_00[1]","26.s12par_10[b1]","27.s12par_20[b1^2]","28.s12par_01[b2]","29.s12par_001[bs]",
            "30.s12perp_00[1]","31.s12perp_10[b1]","32.s12perp_20[b1^2]","33.s12perp_01[b2]","34.s12perp_001[bs]\n");
            
   for (int i=0; i<golists.Nr; i++) {
	    fprintf(outstr,CLPT_OUT_FMT ,	
				rArrays.rTab[i],
				rArrays.xi_00T[i],
				rArrays.xi_10T[i],
				rArrays.xi_20T[i],
				rArrays.xi_01T[i],
				rArrays.xi_02T[i],
				rArrays.xi_11T[i],
				rArrays.xi_eftT[i],
				rArrays.xi_001T[i],
				rArrays.xi_002T[i],
				rArrays.xi_101T[i],
				rArrays.xi_011T[i],         
				rArrays.xi_LT[i],
				rArrays.xi_zaT[i],
				rArrays.xi_AT[i], 
				rArrays.xi_WT[i], 
				rArrays.v_00T[i],
				rArrays.v_10T[i],
				rArrays.v_20T[i],
				rArrays.v_01T[i],
				rArrays.v_11T[i],
				rArrays.v_001T[i],
				rArrays.v_101T[i],
				rArrays.v_eftT[i],         
				rArrays.spar_00T[i],
				rArrays.spar_10T[i],
				rArrays.spar_20T[i],
				rArrays.spar_01T[i],
				rArrays.spar_001T[i],            
				rArrays.sperp_00T[i],
				rArrays.sperp_10T[i],
				rArrays.sperp_20T[i],
				rArrays.sperp_01T[i],
				rArrays.sperp_001T[i]
			);
	};
	    fclose(outstr);	
};





void free_variables(void){

    free(qArrays.nabla2xiT);
    free(qArrays.iBessel4T);
    free(qArrays.iBessel2T);
    free(qArrays.iBessel0T); 
    free(qArrays.V10T);
    free(qArrays.dotVbT);     
    free(qArrays.ddotVaT);   
    free(qArrays.dotVaT);   
    free(qArrays.U20T);
    free(qArrays.U11T);
    free(qArrays.UloopT);
    free(qArrays.ULT);
    free(qArrays.Y10T);
    free(qArrays.X10T);
    free(qArrays.TT);
    free(qArrays.preVT);
    free(qArrays.YloopT);
    free(qArrays.XloopT);
    free(qArrays.YLT);
    free(qArrays.XLT);
    free(qArrays.qTab);

    free(qArraysd.nabla2xiT);
    free(qArraysd.iBessel4T);
    free(qArraysd.iBessel2T);
    free(qArraysd.iBessel0T); 
    free(qArraysd.V10T);
    free(qArraysd.dotVbT);     
    free(qArraysd.ddotVaT);   
    free(qArraysd.dotVaT);   
    free(qArraysd.U20T);
    free(qArraysd.U11T);
    free(qArraysd.UloopT);
    free(qArraysd.ULT);
    free(qArraysd.Y10T);
    free(qArraysd.X10T);
    free(qArraysd.TT);
    free(qArraysd.preVT);
    free(qArraysd.YloopT);
    free(qArraysd.XloopT);
    free(qArraysd.YLT);
    free(qArraysd.XLT);
    free(qArraysd.qTab);

    free(kArrays.kT);
    free(kArrays.Q1T);
    free(kArrays.Q2T);
    free(kArrays.Q3T);
    free(kArrays.Q5T);
    free(kArrays.Q8T);
    free(kArrays.Qs2T);
    free(kArrays.R1T);
    free(kArrays.R2T);
    free(kArrays.pklT);

    free(kArraysd.kT);
    free(kArraysd.Q1T);
    free(kArraysd.Q2T);
    free(kArraysd.Q3T);
    free(kArraysd.Q5T);
    free(kArraysd.Q8T);
    free(kArraysd.Qs2T);
    free(kArraysd.R1T);
    free(kArraysd.R2T);
    free(kArraysd.pklT);



    free(rArrays.rTab);
    free(rArrays.xi_00T);
    free(rArrays.xi_10T);
    free(rArrays.xi_20T);
    free(rArrays.xi_01T);
    free(rArrays.xi_02T);
    free(rArrays.xi_11T);
    free(rArrays.xi_eftT);
    free(rArrays.xi_001T);
    free(rArrays.xi_002T);
    free(rArrays.xi_101T);
    free(rArrays.xi_011T); 
    free(rArrays.xi_AT); 
    free(rArrays.xi_WT); 
    free(rArrays.v_00T);
    free(rArrays.v_10T);
    free(rArrays.v_20T);
    free(rArrays.v_01T);
    free(rArrays.v_11T);
    free(rArrays.v_001T);
    free(rArrays.v_101T);
    free(rArrays.v_eftT);         
    free(rArrays.spar_00T);
    free(rArrays.spar_10T);
    free(rArrays.spar_20T);
    free(rArrays.spar_01T);
    free(rArrays.spar_001T);            
    free(rArrays.sperp_00T);
    free(rArrays.sperp_10T);
    free(rArrays.sperp_20T);
    free(rArrays.sperp_01T);
    free(rArrays.sperp_001T);
    free(rArrays.xi_zaT);
    free(rArrays.xi_LT);

    free(rArraysd.rTab);
    free(rArraysd.xi_00T);
    free(rArraysd.xi_10T);
    free(rArraysd.xi_20T);
    free(rArraysd.xi_01T);
    free(rArraysd.xi_02T);
    free(rArraysd.xi_11T);
    free(rArraysd.xi_eftT);
    free(rArraysd.xi_001T);
    free(rArraysd.xi_002T);
    free(rArraysd.xi_101T);
    free(rArraysd.xi_011T); 
    free(rArraysd.xi_AT); 
    free(rArraysd.xi_WT); 
    free(rArraysd.v_00T);
    free(rArraysd.v_10T);
    free(rArraysd.v_20T);
    free(rArraysd.v_01T);
    free(rArraysd.v_11T);
    free(rArraysd.v_001T);
    free(rArraysd.v_101T);
    free(rArraysd.v_eftT);         
    free(rArraysd.spar_00T);
    free(rArraysd.spar_10T);
    free(rArraysd.spar_20T);
    free(rArraysd.spar_01T);
    free(rArraysd.spar_001T);            
    free(rArraysd.sperp_00T);
    free(rArraysd.sperp_10T);
    free(rArraysd.sperp_20T);
    free(rArraysd.sperp_01T);
    free(rArraysd.sperp_001T);
    free(rArraysd.xi_zaT);
    free(rArraysd.xi_LT);









};




/*


local real sigma8_int(real y)
{
    real p;
    real PSL;

    p = rpow(10.0,y);
    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);

    return p*PSL;
}


local real sigma8(void)
{
    real result;
    real kmin, kmax;
    real ymin, ymax;
    real EPSQ = 0.000001;
    int KK = 5;

    kmin = kPS[1];
    kmax = kPS[nPSLT];
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    result= (1.0/SIXPI2)*rlog(10.0)
    *qromo(sigma8_int,ymin,ymax,midpnt,EPSQ,KK);
    
    
    //~ fprintf(stdout,"\n sigma2v= %g \n",result);

    return result;

}



local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;   
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);  
    return (psftmp);
}

*/
