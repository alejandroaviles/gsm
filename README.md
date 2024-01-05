# gsm
A Gaussian Streaming Code to compute tracers 2 pt correlation function in redshift space


#
Alejandro Aviles

avilescervantes@gmail.com

#

Other people who contributed to this code:

Sadi Ramirez-Solano

Mario A. Rodriguez-Meza 

Mariana Vargas-Magaña

Sebastien Fromenteau


#

This code computes (very fast) correlation functions for generic tracers in redshift space using the Gaussian Streaming Model and Convolution Lagrangian Perturbation Theory.



## Run


Git clone

```
git clone https://github.com/alejangroaviles/gsm.git
```

or download it from http://www.github.com/alejandroaviles/gsm


Compile and run:

```
/gsm$ make
/gsm$ ./gsm
```

This will compute the correlation function for the input linear power spectrum /gsm/Input/psLCDM.in (in Mpc/h units), with default parameters.


For help:

```
/MGPT$ ./mgpt -help
```


In help you can see how to change parameters, in the form [option]=[value], for example:

```
/MGPT$ ./gsm fnamePS=pkl_z05.dat zout=0.5 om=0.3 b1=0.7 bs=0.1 sFoG=-10 suffix=_run2
```

computes the 2pcf for the linear input /gsm/Input/fnamePS=pkl_z05.dat at redshift z=0.5 with matter abundance Omega_m = 0.3. Lagrangian bias parameters linear: b1=0.7 and tidal bias bs=0.1 and EFT parameter sFoG=-10 in (Mpc/h)^2 units (this is the similar to Fingers of God EFT parameter). The output files will have a suffix _run2. (z and Omega_m are used only to calculate the logarithmic growth rate f.)  


Alternatively you can run the code with a parameters file:

```
/gsm$ ./gsm parameters.in
```


The main output is gsm/rsd_multipoles.dat file with four columns

| column  | function  |
| ------------: |:--------------------| 
| #1            | r   (in Mpc/h)      |
| #2            | \xi_0  monopole     |  
| #3            | \xi_2  quadrupole   |  
| #4            | \xi_4  hexadecapole |  



Other outputs for intermediate calculations are located in the folder gsm/Output/







## References

If you use this code please refer to this repository.

The theory is based mainly on the following papers:

1. Lile Wang, Beth Reid, Martin White, https://arxiv.org/abs/1306.1804
2. Zvonimir Vlah, Emanuele Castorina, Martin White https://arxiv.org/abs/1609.02908
3. Alejandro Aviles https://arxiv.org/abs/1805.05304
4. Georgios Valogiannis, Rachel Bean, Alejandro Aviles, https://arxiv.org/abs/1909.05261






