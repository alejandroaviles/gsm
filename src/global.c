/*==============================================================================
 NAME: global.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com), ...
 ================================================================================ 
*/
#include "globaldefs.h"
#include "protodefs.h"

void global_variables(void){


	
	golists.rmin=1.0;
	golists.rmax= cmd.smax + 40.;
	golists.Nr =(int)( golists.rmax - golists.rmin + 1.0);
	
	golists.qmin = 0.001;     
    golists.qmax = golists.rmax + 50.;    
    golists.Nq = (int)floor(log10(golists.qmax/golists.qmin)*(real)cmd.NqperLogDecade); 
    
    golists.kmin = cmd.kmin;
    golists.kmax = cmd.kmax;
    golists.Nk = cmd.Nk;
    
    
    golists.smin = cmd.smin;
    golists.smax = cmd.smax;
    golists.Ns = cmd.Ns;
	
	//~ fprintf(stdout,"\nglists.qmin=%g, glists.qmax=%g, glists.Nq=%d",golists.rmin, golists.rmax, golists.Nr);
	//~ fprintf(stdout,"\nglists.rmin=%g, glists.rmax=%g, glists.Nr=%d",golists.qmin, golists.qmax, golists.Nq);

	
    kArrays.kT = malloc( golists.Nk * sizeof(real));
    kArrays.Q1T   = malloc( golists.Nk * sizeof(real));
    kArrays.Q2T   = malloc( golists.Nk * sizeof(real));
    kArrays.Q3T   = malloc( golists.Nk * sizeof(real));
    kArrays.Q5T   = malloc( golists.Nk * sizeof(real));
    kArrays.Q8T   = malloc( golists.Nk * sizeof(real));
    kArrays.Qs2T = malloc( golists.Nk * sizeof(real));
    kArrays.R1T   = malloc( golists.Nk * sizeof(real));
    kArrays.R2T   = malloc( golists.Nk * sizeof(real));
    kArrays.pklT  = malloc( golists.Nk * sizeof(real));
    	
    kArraysd.kT  = malloc( golists.Nk * sizeof(real));
    kArraysd.Q1T   = malloc( golists.Nk * sizeof(real));
    kArraysd.Q2T   = malloc( golists.Nk * sizeof(real));
    kArraysd.Q3T   = malloc( golists.Nk * sizeof(real));
    kArraysd.Q5T   = malloc( golists.Nk * sizeof(real));
    kArraysd.Q8T   = malloc( golists.Nk * sizeof(real));
    kArraysd.Qs2T = malloc( golists.Nk * sizeof(real));
    kArraysd.R1T   = malloc( golists.Nk * sizeof(real));
    kArraysd.R2T   = malloc( golists.Nk * sizeof(real));
    kArraysd.pklT  = malloc( golists.Nk * sizeof(real));

    qArrays.qTab = malloc( golists.Nq * sizeof(real));
    qArrays.XLT = malloc( golists.Nq * sizeof(real));
    qArrays.YLT = malloc( golists.Nq * sizeof(real));
    qArrays.XloopT = malloc( golists.Nq * sizeof(real));
    qArrays.YloopT = malloc( golists.Nq * sizeof(real));
    qArrays.preVT = malloc( golists.Nq * sizeof(real));
    qArrays.TT = malloc( golists.Nq * sizeof(real));
    qArrays.X10T = malloc( golists.Nq * sizeof(real));
    qArrays.Y10T = malloc( golists.Nq * sizeof(real));
    qArrays.ULT = malloc( golists.Nq * sizeof(real));
    qArrays.UloopT = malloc( golists.Nq * sizeof(real));
    qArrays.U11T = malloc( golists.Nq * sizeof(real));
    qArrays.U20T = malloc( golists.Nq * sizeof(real));
    qArrays.dotVaT = malloc( golists.Nq * sizeof(real));
    qArrays.ddotVaT = malloc( golists.Nq * sizeof(real));
    qArrays.dotVbT = malloc( golists.Nq * sizeof(real)); 
    qArrays.V10T = malloc( golists.Nq * sizeof(real));
    qArrays.iBessel0T = malloc( golists.Nq * sizeof(real));
    qArrays.iBessel2T = malloc( golists.Nq * sizeof(real));
    qArrays.iBessel4T = malloc( golists.Nq * sizeof(real));
    qArrays.nabla2xiT = malloc( golists.Nq * sizeof(real));

    qArraysd.qTab = malloc( golists.Nq * sizeof(real));
    qArraysd.XLT = malloc( golists.Nq * sizeof(real));
    qArraysd.YLT = malloc( golists.Nq * sizeof(real));
    qArraysd.XloopT = malloc( golists.Nq * sizeof(real));
    qArraysd.YloopT = malloc( golists.Nq * sizeof(real));
    qArraysd.preVT = malloc( golists.Nq * sizeof(real));
    qArraysd.TT = malloc( golists.Nq * sizeof(real));
    qArraysd.X10T = malloc( golists.Nq * sizeof(real));
    qArraysd.Y10T = malloc( golists.Nq * sizeof(real));
    qArraysd.ULT = malloc( golists.Nq * sizeof(real));
    qArraysd.UloopT = malloc( golists.Nq * sizeof(real));
    qArraysd.U11T = malloc( golists.Nq * sizeof(real));
    qArraysd.U20T = malloc( golists.Nq * sizeof(real));
    qArraysd.dotVaT = malloc( golists.Nq * sizeof(real));
    qArraysd.ddotVaT = malloc( golists.Nq * sizeof(real));
    qArraysd.dotVbT = malloc( golists.Nq * sizeof(real)); 
    qArraysd.V10T = malloc( golists.Nq * sizeof(real));
    qArraysd.iBessel0T = malloc( golists.Nq * sizeof(real));
    qArraysd.iBessel2T = malloc( golists.Nq * sizeof(real));
    qArraysd.iBessel4T = malloc( golists.Nq * sizeof(real));
    qArraysd.nabla2xiT = malloc( golists.Nq * sizeof(real));
	
	    rArrays.rTab= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_00T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_10T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_20T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_01T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_02T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_11T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_eftT= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_001T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_002T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_101T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_011T= malloc(golists.Nr * sizeof(real) );         
		rArrays.xi_LT= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_zaT= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_AT= malloc(golists.Nr * sizeof(real) ); 
		rArrays.xi_WT= malloc(golists.Nr * sizeof(real) ); 
		rArrays.v_00T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_10T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_20T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_01T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_11T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_001T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_101T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_eftT= malloc(golists.Nr * sizeof(real) );         
		rArrays.spar_00T= malloc(golists.Nr * sizeof(real) );
		rArrays.spar_10T= malloc(golists.Nr * sizeof(real) );
		rArrays.spar_20T= malloc(golists.Nr * sizeof(real) );
		rArrays.spar_01T= malloc(golists.Nr * sizeof(real) );
		rArrays.spar_001T= malloc(golists.Nr * sizeof(real) );            
		rArrays.sperp_00T= malloc(golists.Nr * sizeof(real) );
		rArrays.sperp_10T= malloc(golists.Nr * sizeof(real) );
		rArrays.sperp_20T= malloc(golists.Nr * sizeof(real) );
		rArrays.sperp_01T= malloc(golists.Nr * sizeof(real) );
		rArrays.sperp_001T= malloc(golists.Nr * sizeof(real) );
		
	    rArraysd.rTab= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_00T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_10T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_20T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_01T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_02T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_11T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_eftT= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_001T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_002T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_101T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_011T= malloc(golists.Nr * sizeof(real) );         
		rArraysd.xi_LT= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_zaT= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_AT= malloc(golists.Nr * sizeof(real) ); 
		rArraysd.xi_WT= malloc(golists.Nr * sizeof(real) ); 
		rArraysd.v_00T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_10T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_20T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_01T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_11T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_001T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_101T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_eftT= malloc(golists.Nr * sizeof(real) );         
		rArraysd.spar_00T= malloc(golists.Nr * sizeof(real) );
		rArraysd.spar_10T= malloc(golists.Nr * sizeof(real) );
		rArraysd.spar_20T= malloc(golists.Nr * sizeof(real) );
		rArraysd.spar_01T= malloc(golists.Nr * sizeof(real) );
		rArraysd.spar_001T= malloc(golists.Nr * sizeof(real) );            
		rArraysd.sperp_00T= malloc(golists.Nr * sizeof(real) );
		rArraysd.sperp_10T= malloc(golists.Nr * sizeof(real) );
		rArraysd.sperp_20T= malloc(golists.Nr * sizeof(real) );
		rArraysd.sperp_01T= malloc(golists.Nr * sizeof(real) );
		rArraysd.sperp_001T= malloc(golists.Nr * sizeof(real) );		
		
	
	
	
	
}
