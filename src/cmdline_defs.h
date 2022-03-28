/*==============================================================================
 HEADER: cmdline_defs.h		[gsm]    
 Alejandro Aviles (other collaborators: Mario A. Rodriguez-Meza ...)
 * ==============================================================================
*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	""
#define HEAD2	"GSM code."
#define HEAD3	"..."
   //~ fprintf(outstr,"InputPklFile=%s, redshift z=%g, OmegaM=%g, h=%g, b1=%g, \
                         //~ b2=%g, bs=%g, c1eft=%g, c2eft=%g, s2FoG=%g,\
                         //~ For Precision: q-funtcionts: NqperLogDecade=%g, Nk_qFunctionsQuad=%g, \
                         //~ gsm: gsm_width=%g,gsm_sizeyT=%d, gsm_NGL=%d\n",
                //~ cmd.fnamePS, cmd.xstop, cmd.om, cmd.h, cmd.b1, cmd.b2, cmd.bs, cmd.c1eft, cmd.c2eft,cmd.s2eft,
               //~ cmd.NqperLogDecade,cmd.Nk_qFunctionsQuad, cmd.gsm_width,cmd.gsm_sizeyT,cmd.gsm_NGL);
string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",                   ";Parameter input file. Overwritten by what follows",
//
// Power spectrum table:
    "fnamePS=psLCDM.in",    ";Input filename power spectrum table (k,P(k)). (At redshift zout)",
    "zout=0.0",                         ";Output redshift value",":z",
    "Om=0.281",                     ";Omega matter value (z=0)",":OmegaM",
    "h=0.697",                         ";Hubble parameter",
    "b1=0.0",                           ";Lagrangian linear bias",
    "b2=0.0",                           ";Lagrangian second order bias",
    "bs=0.0",                           ";Lagrangian tidal bias ",
    "c1eft=0.0",                      ";eft counterterm in xi(r), (almost) degenerate with nabla2 bias",
    "c2eft=0.0",                      ";eft counterterm in pairwise velocity infall (not implemented yet)",
    "sFoG=-2.0",                     ";eft counterterm in the pairwise velocity dispersion (similar to sigma_FoG)",":sigma2eft",
// output gsm multipoles
    "smin=1",                     ";Output rsd 2pcf multipoles s minimum",
    "smax=130",                     ";Output rsd 2pcf multipoles s maximum",
    "Ns=100",                     ";number of s in rsd 2pcf multipoles output",
//  k-functions
    "kmin=1e-3",                    ";kmin in k-functions output",
    "kmax=100",                     ";kmax in k-functions output",
    "Nk=250",                       ";Total number of k in  k-functions output",":nk",
    "nquadSteps=200",               ";Number of k´s from the power spectrum table to integrate (trapezoid)",":nquad",
//  q points: Control precision in q functions:
    "NqperLogDecade=100",     ";Total number of q in q-functions output per log decade",":nq",   
    "Nk_qFunctionsQuad=1200",     ";Total number of k in q-functions quadrature",":nkq",   
//  For GSM integration:
    "gsm_width=100",     ";Width of Gaussian in GSM integration",   
    "gsm_sizeyT=160",     ";Number of points in 1-d GSM integration",   
    "gsm_NGL=16",     ";Number of GL points in angular integration for rsd multipoles",      
    "suffixModel=",                 ";Suffix model to add to output filenames", ":suffix",
    "modelParamfile=",              ";If mgmodel=USER, to use the model in models_user.h", ":mpf",
    NULL,
};

#endif // ! _cmdline_defs_h





/*
//~ // CLPT correlation functions table:
    "rmin=1",                      "; NOT USED rmin of the range for CLPT correlation functions NOT USED",
    "rmax=210",                     "; NOT USED rmax of the range for CLPT correlation functions NOT USED",
    "Nr=210",                       "; NOT USED Total number of r´s in the CLPT correlation function",":nr",
//~ // Modified gravity model parameters:
    "mgModel=LCDM",                 ";Modified gravity model to study (HS or DGP), default is LCDM", ":mgm",
    //~ "suffixModel=",                 ";Suffix model to add to output filenames", ":suffix",
    "fR0=1.0e-5",                   ";Hu-Sawicky f_R0",
    "screening=1.0",                ";set to =0 if you want no screenings", ":sc",
//~ // DGP:
    "epsDGP=-1.0",                  ";epsilon DGP parameter, =-1 normal branch, =1 self-accelerating branch",":epsdgp",
    "rcDGP=1.0",                    ";crossover scale parameter in DGP, (in units of 1/H0)",":rcdgp",
//~ //
    //~ "modelParamfile=",              ";If mgmodel=USER, to use the model in models_user.h", ":mpf",
//~ //
//~ // Background cosmology:
    "Om=0.281",                     ";Omega matter value (z=0)",":om",
    "OL= 1 - Om",                   ";Omega Lambda value (z=0). Only works for DGP!",":ol",
    "h=0.697",                      ";Hubble parameter",
//~ //
//~ // Differential equations evolution parameters:
    "etaini=-4.0",                  ";Initial conformal time value :: Log[1/(1 + zini)]",
    "deta=2/5",                     ";Conformal time integration step",
    "detamin=0.",                   ";Min conformal time integration step size",
    "epsSolver=1.0e-4",             ";Differential equations solver tolerance error parameter",":epssolver",
    "maxnsteps=10000",              ";Maximum number of integration steps", ":maxn",
    "solverMethod=rkqs",	        ";Integration method to use", ":solver",
//~ //
//~ // Quadrature parameters:
    "quadratureMethod=trapezoid3",   ";Quadrature method to use", ":quadm",
    "ngausslegpoints=16",           ";Number of Gauss-Legendre of integration points", ":nglpts",
    "epsquad=1.0e-6",               ";Quadrature tolerance error parameter (open Romberg method: romberg)",
//~ //
//~ // Post processing parameters:
    "postprocessing=false",			";Post processing options", ":pp",
    "options=",                     ";Various control options", ":opt",
*/
