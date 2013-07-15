/***************************************************************************
 *   Copyright (C) 2008 by Andreas Kuepper                                 *
 *   akuepper@astro.uni-bonn.de                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/***************************************************************************
 *  star array description:                                                *
 *   star[x][0]    -  name                                                 *
 *   star[x][1-3]  -  x, y, z [pc]                                         *
 *   star[x][4-6]  -  vx, vy, vz [km/s]                                    *
 *   star[x][7]    -  stellar type (see Nbody documentation)               *
 *   star[x][8]    -  luminosity [L_sun]                                   *
 *   star[x][9]    -  radius [R_sun]                                       *
 *   star[x][10]   -  mass [M_sun]                                         *
 *   star[x][11]   -  potential energy [M_sun km^2/s^2]                    *
 *   star[x][12]   -  kinetic energy [M_sun km^2/s^2]                      *
 *   star[x][13]   -  status (see subroutine get_status)                   *
 *   star[x][14]   -  escape status (to avoid multiple escapes)            *
 *                                                                         *
 ***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<sys/stat.h>



/***************
 * DEFINITIONS *
 ***************/

#define G 0.0043009211                  //in pc*km^2/s^2*M_sun
#define Pi 3.14159265
#define PI 3.14159265
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )


//Tidal field standard parameters (will get changed by Nbody6 output)

//point mass bulge
double M1 = 9.565439E+10;               //[M_sun]
double b1 = 387.3;                      //pc

//Gamma model bulge / Hernquist bulge (from Law, Majewski & Johnston 2009)
double M1_GAMMA = 3.4e10;               //[M_sun]
double b1_GAMMA = 700.0;                //[pc]

//Myamoto disk
double M2 = 1.0e11;                     //[M_sun]
double a2 = 6500.0;                     //[pc]
double b2 = 260.0;                      //[pc]

//Halo
double q_halo = 1.0;                    //halo flattening along z axis

//Log Halo
double VCIRC = 0.0;                     //[km/s]
double RCIRC = 1.0;                     //[pc]

//NFW halo
double MNFW = 0.0;                      //[M_sun]
double RNFW = 1.0;                      //[pc]


//Solar parameters
double const vxsun = 11.1;              //+0.69/−0.75 - solar motion with respect to the LSR from Schoenrich, Binney & Dehnen (2010) [km/s]
double const vysun = 12.24;             //+0.47−0.47
double const vzsun = 7.25;              //+0.37−0.36

double const alphaGNP = 192.859508;     //Galactic north pole in J2000 coordinates
double const deltaGNP = 27.128336;
double const PAGNP = 122.932;           //Position angle with respect to equatorial pole

double const vLSR = 239.5;              //circular velocity at solar radius
double const rgalsun = 8330.0;          //solar Galactocentric radius [pc] (standard = 8330.0; Gillessen et al. 2009)




/****************
 * SUBROUNTINES *
 ****************/

int cmpmy(float *x1, float *x2);
int cmpmy2(float *x1, float *x2);
int readin0(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt);
int readin1(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, int outputtypet, int Nmaxt);
int readin2(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt);
int readin3(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt);
int readin4(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt);
int readin5(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt);
int readin6(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt);
void shellsort_reverse_1d(double *array, int N);
void shellsort_1d(double *array, int N);	
void shellsort(double **array, int N, int k);
void shellsort_reverse(double **array, int N, int k);
int bound_star_analysis(int Nt, double **start, int **neighborlistt, int neighborst);
int bound_star_analysisfast(int Nt, double **start, int **neighborlistt, int neighborst, double mxt, double myt, double mzt);
int get_densitycenter(int Nt, double **start, int **neighborlistt, double *mxt, double *myt, double *mzt, int densitylevel);
int get_omega(double *RGt, double *VGt, double *omegat, double *Lt);
int get_rtide(double *rtidet, int Nt, double **start, double omegat, double mxt, double myt, double mzt, double dphit, double *RGt);
int get_status(int Nt, double **start, double rtidet, double mxt, double myt, double mzt, int *nt, int *nbint, double mvxt, double mvyt, double mvzt, double rh, double omega, int quick_and_dirty_analysis);
int get_binaries(int Nt, double **start, int **neighborlistt, int neighborst, int *Ntott, int *multiplicityarray, int binoutput, FILE *bint, double timet, int outputtypet);
int get_velocitydispersion(int Nt,double **start,int statust, double *sigmat, double *mvxt, double *mvyt, double *mvzt);
int get_energy(int Nt, double **start, int statust, double *Etemp, double mvxt, double mvyt, double mvzt);
int radii(int Nt, double **start, int statust, double *LRt, double mxt, double myt, double mzt, double *Mt);
int Kradii(int Nt, double **start, int statust, double *LRt, double mxt, double myt, double mzt, double *Mt, double *racc);
int surfacedensityprofile(int Nt, double **start, double **profilet, double mxt, double myt, double mzt, double mvxt, double mvyt, double mvzt, int viewt, int numberradialbinst, double *VGt, double **profile_inside_rtidet, double **profile_bound_starst);
int get_orbit(double *RGt,double *VGt,double **along_the_orbitt,int pcnumberintervalst,int signt);
void get_force(double *x, double *v, double *a);
void do_step(double dt, double *x, double *v);
int rcsuche(int Nt, double **start, int **neighborlistt, double mxt, double myt, double mzt, int densitylevelt, double *Rct);	
int escape(FILE *esc, int Nt, double **start, double rtidet, double mxt, double myt, double mzt, double mvxt, double mvyt, double mvzt, double omegat, double timet, double sigmat, double kTt, double Mt, int outputtypet, double **timebint, int timebinsizet, double thresholdt);
void convert(double *x, double *v, double *dsun, double *vrsun, double *vr, double *l, double *b, double *lcosb, double *RA, double *DEC, double *mu_alpha, double *mu_alphacosdelta, double *mu_delta, double *mutemp, double *PAtemp, int coordtype, int vcoordtype, int radiococo, double vLSRtemp);





/*****************
 * MAIN ROUNTINE *
 *****************/

int main(int argc , char * argv[]){
	
	/**********************
	 * CONTROL PARAMETERS *
	 **********************/
	
	int code = 2;                           //0: McLuster input; 1: NBODY4 (32bit); 2: NBODY6; 3: Holgers NBODY6; 4: Elhams NBODY6; 5: NBODY4new (64bit); 6: NBODY6 (Holger's cluster)
	int Nmax = 300000;                      //size of data array, minimum is twice the initial star number
	int Mend = 0;                           //stop analysis if M<Mend [M_sun]
	int analysistimestep = 1;               //time step for analysis [Myr]; 1: all available time steps
	int quick_and_dirty_analysis = 0;       //calculation of binding energy and neighbours; 0: detailed; 1: fast
	int standard_densitylevel = 5;          //Casertano & Hut density estimator level (j) (standard = 5)
	int standard_neighbornumber	= 20;       //length of neighbour list [20 is a reasonable size]
	int standard_status = 2;                //calculate basic cluster parameters by taking into account: 1: all bound stars within rtide; 2: all stars within rtide; 3: all stars within 2*rtide
	int binoutput = 0;                      //output of binary data; 0: summary; 1: detailed
	int outputtype = 0;                     //0: output of major quantities to STDOUT; 1: debug mode
	int outputstep = 1;                     //output interval for all.txt, all.bin, .sdp, y.bin and y.txt [Myr]; 1: every available time step
	int alloutput = 0;                      //0: output of all stars to all.txt and all.bin off; 1: on
	int tailanalysis = 0;                   //0: analysis and output of tidal-tail stars to y.bin and y.txt off; 1: on
	int potential = 2;                      //0: isolated; 1: near-field approximation; 2: point-mass potential or full galactic potential (parameters set by Nbody6 output)
	int sunviewtime = 12000;                //output interval for snapshots from the position of the sun [Myr]
	int sdp = 0;                            //0: output of surface-density profiles to .sdp off; 1: on
	int view = 1;                           //projection of surface-density profiles; 1: xy-plane; 2: xz-plane
	int binarysearch = 0;                   //search for binaries and replace them by centre-of-mass particles? 0: no; 1: yes
	int trackescapers = 0;                  //search for escapers and determine escape velocities; 0: off; 1: on
	int timebinsize = 100;                  //size of time bins for dynamical temperature-luminosity diagram [Myr] (standard = 100)
	double threshold = 5.1;                 //threshold between evaporaters and ejecters (see Kuepper et al. 2008 eq. 10)
	double limit = 0.0;                     //detection limit for stellar masses in snapshots [Msun]
	double vlimit = 100000.0;//26.0			//detection limit in V magnitudes
	double RG0[3] = {0.0, 8500.0, 0.0};     //position of cluster if not provided by input file
	double VG0[3] = {220.0, 0.0, 0.0};		//velocity of cluster if not provided by input file
	
	
	/******************
	 * INITIALIZATION *
	 ******************/
	
	clock_t t1, t2;                       //start stop watch
	t1 = clock(); 
	
	srand48((unsigned) time(NULL));
	
	int i,j,l,q;
	int Ntot, N, Ninit;
	double time, tbar, rtide, rbar, vbar, mbar, omega, Etot, rcore, Minit, vr;
	double RG[3], VG[3];
	double mx = 0.0, my = 0.0, mz = 0.0;
	double mxold, myold, mzold;
	int densitylevel; 
	double sigmauncor;
	double mvxuncor = 0.0, mvyuncor = 0.0, mvzuncor = 0.0;
	double sigma;
	double mvx = 0.0, mvy = 0.0, mvz = 0.0;
	double dphi;
	int status;
	double xgal, ygal, zgal, Rtemp, ltemp, btemp;
	
	
	
	/********************************
	 * READ-IN AND ANALYSIS OF DATA *
	 ********************************/
	
	printf("Number of files: %i\n\n",argc-1);
	
	for (q=1;q<argc;q++) {
		
		/*********************
		 * MEMORY ALLOCATION *
		 *********************/
		
		//main data array
		int columns = 15;
		double **star;
		star = (double **)calloc(Nmax,sizeof(double *));
		for (j=0;j<Nmax;j++){
			star[j] = (double *)calloc(columns,sizeof(double));
			if (star[j] == NULL) {
				printf("\nNULL; memory allocation failed!\n");
				return 0;
			}
		}
		
		//neighbor array
		int neighbors = standard_neighbornumber;
		int **neighborlist;
		neighborlist = (int **)calloc(Nmax,sizeof(int *));
		for (j=0;j<Nmax;j++){
			neighborlist[j] = (int *)calloc(neighbors,sizeof(int));
			if (neighborlist[j] == NULL) {
				printf("\nNULL; memory allocation failed!\n");
				return 0;
			}
		}
		
		//surface-density arrays
		int numberradialbins = 100;
		double **profile;
		profile = (double **)calloc(numberradialbins,sizeof(double *));
		for (l=0;l<numberradialbins;l++){
			profile[l] = (double *)calloc(8,sizeof(double));
			if (profile[l] == NULL) {
				printf("\nNULL; memory allocation failed!\n");
				return 0;
			}
		}
		double **profile_inside_rtide;
		profile_inside_rtide = (double **)calloc(numberradialbins,sizeof(double *));
		for (l=0;l<numberradialbins;l++){
			profile_inside_rtide[l] = (double *)calloc(8,sizeof(double));
			if (profile_inside_rtide[l] == NULL) {
				printf("\nNULL; memory allocation failed!\n");
				return 0;
			}
		}
		double **profile_bound_stars;
		profile_bound_stars = (double **)calloc(numberradialbins,sizeof(double *));
		for (l=0;l<numberradialbins;l++){
			profile_bound_stars[l] = (double *)calloc(8,sizeof(double));
			if (profile_bound_stars[l] == NULL) {
				printf("\nNULL; memory allocation failed!\n");
				return 0;
			}
		}
		
		//temperature-luminosity arrays
		i = 20000/timebinsize;
		double **timebin;
		timebin = (double **)calloc(i,sizeof(double *));
		for (l=0;l<i;l++) {
			timebin[l] = (double *)calloc(5,sizeof(double));
			if (timebin[l] == NULL) {
				printf("\nNULL; memory allocation failed!\n");
				return 0;
			}
			//timebin[j][0] = number
			//timebin[j][1] = temperature
			//timebin[j][2] = evaporation luminosity
			//timebin[j][3] = ejection luminosity
			//timebin[j][4] = total mass
		}
		
		
		
		/*********************************
		 * DETERMINE NUMBER OF POS-FILES *
		 *********************************/
		
		FILE *fz;
		char filename[200];
		int poscounter=1;
		if (code) {
			j = 2;
			while (j) {
				sprintf(filename, "%s.POS%i", argv[q],j);
				fz = fopen(filename,"rb");
				if (fz==NULL) {
					j = 0;
				} else {
					poscounter = j;
					j++;
				}
			}
			
			sprintf(filename, "%s.POS", argv[q]);
			fz = fopen(filename,"rb");
			if (fz==NULL) {
				printf("\nFile %s not found!\n", filename);
				return 0;
			} 
		} else {
			sprintf(filename, "%s.txt", argv[q]);
			fz = fopen(filename,"r");
			if (fz==NULL) {
				printf("\nFile %s not found!\n", filename);
				return 0;
			} 
			j = 0;
		}
		int poscountertemp = 1;
		
		
		
		/*********************
		 * INITIALIZE OUTPUT *
		 *********************/
		
		char *outputfilename;
		outputfilename = strtok(filename, ".");            
		mkdir(outputfilename,S_IRWXU|S_IRWXG|S_IRWXO );
		
		FILE *new;
		char summary[200];
		sprintf(summary, "%s.new",outputfilename);
		new = fopen(summary,"w");
		
		FILE *gal;
		char orbitfile[200];
		sprintf(orbitfile, "%s.gal",outputfilename);
		gal = fopen(orbitfile,"w");
		
		FILE *bin;	
		char binaryfile[200];
		sprintf(binaryfile, "%s.bin",outputfilename);
		bin = fopen(binaryfile,"w");
		
		FILE *esc;	
		if (trackescapers) {
			char escaperfile[200];
			sprintf(escaperfile, "%s.esc",outputfilename);
			esc = fopen(escaperfile,"w");
		}
		
		
		
		/*****************************
		 * READ-IN AND ANALYSIS LOOP *
		 *****************************/
		
		double Rh = 1.0;
		int step = 0;
		
		do {
			printf("\n##########################################################\n");
			printf("Beginning analysis of file %s (%i of %i)\n", filename, poscountertemp, poscounter);
			
			long lSize, lPos;
			fseek (fz , 0 , SEEK_END);      // obtain file size
			lSize = ftell (fz);
			rewind (fz);
			if (lSize) {
				printf("File size: %li Bytes\n",lSize);
				printf("##########################################################\n");
				RG[0] = RG0[0];
				RG[1] = RG0[1];
				RG[2] = RG0[2];
				VG[0] = VG0[0];
				VG[1] = VG0[1];
				VG[2] = VG0[2];
				
				do{
					
					/***********
					 * READ-IN *
					 ***********/
					
					for (j=0;j<Nmax;j++) star[j][0] = 0.0;
					
					if (code == 1) {
						readin1(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, star, outputtype, Nmax);
						RG[0]=RG[1]=RG[2]=VG[0]=VG[1]=VG[2]=0;
					} else if (code == 2) {
						readin2(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, star, RG, VG, outputtype, Nmax);
					} else if (code == 3) {
						readin3(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, star, RG, VG, outputtype, Nmax);
					} else if (code == 4) {
						readin4(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, star, RG, VG, outputtype, Nmax);
					} else if (code == 5) {
						readin5(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, star, RG, VG, outputtype, Nmax);
					} else if (code == 6) {
						readin6(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, star, RG, VG, outputtype, Nmax);
					} else {
						readin0(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, star, RG, VG, outputtype, Nmax);
					}
					
					int timetemp;
					timetemp = floor(time+0.1);
					
					
					
					/************
					 * ANALYSIS *
					 ************/
					
					if (timetemp % analysistimestep < 1) {		
						
						if (outputtype) printf("\nT: %g Myrs\n",time);
						
						if (step == 0) Ninit = N;
						Ntot = Ninit;
						
						
						//make neighbour lists, determine binding energies in case of quick_and_dirty_analysis = 0
						if (standard_neighbornumber < N-2)
							neighbors = standard_neighbornumber;
						else
							neighbors = N-2;
						
						if (quick_and_dirty_analysis) 
							bound_star_analysisfast(Ntot, star, neighborlist, neighbors, mx, my, mz);
						else
							bound_star_analysis(Ntot, star, neighborlist, neighbors);
						
						
						//determine cluster centre
						if (standard_densitylevel < N-2) 
							densitylevel = standard_densitylevel;
						else
							densitylevel = N-2;
						mxold = mx;
						myold = my;
						mzold = mz;
						get_densitycenter(Ntot, star, neighborlist, &mx, &my, &mz, densitylevel);
						if (outputtype) printf("mx = %g\tmy = %g\tmz = %g\n",mx,my,mz);
						
                        
						//determine omega
						double L;
						if (potential > 0) {
							get_omega(RG, VG, &omega, &L);
							
							//determine second derivative of potential with respect to R (dphi) [for Kuepper et al 2010, eq. 11]
							if (potential > 1) {
                                
                                double rtemp;
                                double xt[3], vt[3], at[3], at2[3]; //force evaluation for dphi
                                
                                rtemp = sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
                                for (i=0;i<3;i++) {
                                    xt[i] = RG[i]/(1.0*rtemp)*(rtemp-20.0);
                                    vt[i] = 0.0;
                                }
                                get_force(xt,vt,at);
                                for (i=0;i<3;i++) xt[i] = RG[i]/(1.0*rtemp)*(rtemp+20.0);
                                get_force(xt,vt,at2);
                                dphi = (sqrt(at[0]*at[0]+at[1]*at[1]+at[2]*at[2])-sqrt(at2[0]*at2[0]+at2[1]*at2[1]+at2[2]*at2[2]))/40.0;
                                if (outputtype) printf("%g\t%g\t%g\n", at[0], at2[0], dphi);

							} else {
								dphi = 2.0*omega*omega;
							}
							
							//determine rtide
							get_rtide(&rtide, Ntot, star, omega, mx, my, mz, dphi, RG);
							if (rtide == 0.0) {
								rtide = pow(G*Ntot*0.4/sqrt(pow((dphi+omega*omega),2)),1.0/3.0);
							}
							
						} else {
							omega = 1.0E-09;
							rtide = 1.0E09;
						}
						
						if (outputtype) printf("omega = %g\t", omega);
						if (outputtype) printf("rtide = %g\n", rtide);
						
						
						//set status: 0 = normal; 1 = potential escaper; 2 = beyond 1*rtide; 3 = beyond 2*rtide; -nrbin = binary
						//determine binding energy in case of quick_and_dirty_analysis = 1
						int nuncor[4]; int nbinuncor[4];
						for (j=0;j<Nmax;j++) star[j][13] = 0.0;
						get_status(Ntot, star, rtide, mx, my, mz, nuncor, nbinuncor,mvx,mvy,mvz,Rh,omega,quick_and_dirty_analysis);
						
						
						//compute sigma und mv, uncorrected for binaries
						status = standard_status;
						get_velocitydispersion(Ntot,star,status,&sigmauncor,&mvxuncor,&mvyuncor,&mvzuncor);
						if (outputtype) printf("mvx = %g\tmvy = %g\tmvz = %g\t --> sigma = %g (uncorrected)\n",mvxuncor,mvyuncor,mvzuncor,sigmauncor);
						
						
						//determine energies
						double Euncor[3];
						status = standard_status;
						get_energy(Ntot, star, status, Euncor,mvxuncor,mvyuncor,mvzuncor);
						if (outputtype) printf("Ebin = %g\tEkin = %g\tEtot = %g (uncorrected)\n", Euncor[0], Euncor[1], Euncor[2]);
						
						
						//detect and mark new escapers
						for(l=0;l<Ntot;l++) {
							if ((star[l][0]) && (!star[l][14]) && (star[l][13] > 2.0)) {
								star[l][14] = 1.0;
							}
						}
						
						
						//find binaries and calculate f_bin
						//build centre-of-mass particle and set binary-member status to -binnr
						int multiplicityarray[neighbors];
						if (code && binarysearch) get_binaries(Ntot, star, neighborlist, neighbors, &Ntot, multiplicityarray, binoutput,bin,time,outputtype);
						if (binoutput) {
							printf("\n");
							for (j=2;j<10;j++) {
								printf ("%i\t",j);
							}
							printf("\n");
							for (j=2;j<10;j++) {
								printf ("%i\t",multiplicityarray[j]);
							}
							printf("\n");
						} else {
							fprintf (bin,"%.0f\t",time);
							for (j=2;j<10;j++) {
								fprintf (bin,"%i\t",multiplicityarray[j]);
							}
							fprintf (bin,"\n");
						}
						
						
						//set status for all stars including centre-of-mass particles
						//determine binding energy in case of quick_and_dirty_analysis = 1
						int n[4]; int nbin[4];
						for (j=0;j<Nmax;j++) if (star[j][13] >= 0.0) star[j][13] = 0.0;
						get_status(Nmax, star, rtide, mx, my, mz, n, nbin, mvx, mvy, mvz, Rh,omega,quick_and_dirty_analysis);
						if (outputtype) printf("Status (singles) \t0: %i\t1: %i\t2: %i\t3: %i\n", n[0],n[1], n[2], n[3]);
						if (outputtype) printf("Status (binaries)\t0: %i\t2: %i\t3: %i\n", nbin[0], nbin[2], nbin[3]);
						double fbin;
						int ntemp = 0;
						for (j=0; j<status;j++) ntemp += n[j];
						int nbintemp = 0;
						for (j=0; j<status;j++) nbintemp += nbin[j];
						fbin = 1.0*nbintemp/ntemp;
						if (outputtype) printf("fbin = %1.5lf\t%i\n", fbin, nbintemp);
						
						
						//determine corrected sigma and mv
						status = standard_status;
						get_velocitydispersion(Nmax,star,status,&sigma,&mvx,&mvy,&mvz);
						if (outputtype) printf("mvx = %g\tmvy = %g\tmvz = %g\t --> sigma = %g (corrected)\n",mvx,mvy,mvz,sigma);
						
						double mv = 0.0;
						int mvn = 0;
						
						for(l=0;l<Ntot;l++) {
							if ((star[l][0]) && (star[l][13]>=0.0) && (star[l][13]<2.0)) {
								mv += sqrt(pow(star[l][4]-mvx,2)+pow(star[l][5]-mvy,2)+pow(star[l][6]-mvz,2));
								mvn++;
							}
						}
						if (mvn) mv /= mvn;
						
						
						//determine corrected energies
						double E[3], kT;
						get_energy(Nmax, star, status, E, mvx, mvy, mvz);
						kT = E[1]/ntemp;
						
						if (outputtype) printf("Ebin = %g\tEkin = %g\tEtot = %g (corrected)\n", E[0], E[1], E[2]);
						
						
						//determine radii, LR2, LR5, LR10, LR20, LR50/R_h, LR90, R_core
						double LR[6];
						double Rc = 0.0;
						double M;
						status = standard_status;
						radii(Ntot, star, status, LR, mx, my, mz, &M);
						Rh = LR[4];
												
						double KR[11];
						double racc = 0.0; //radius at which a < a0
						if (step == 0) Minit = M;
						
						Kradii(Ninit, star, status, KR, mx, my, mz, &Minit, &racc);
						
						rcsuche(Ntot, star, neighborlist, mx, my, mz, densitylevel, &Rc);
						
						if (outputtype) printf("M = %g\tLR2 = %g\tLR5 = %g\tLR10 = %g\tLR20 = %g\tRh = %g\tLR90 = %g\tRc = %g\n", M, LR[0], LR[1], LR[2], LR[3], LR[4], LR[5], Rc);
						
						
						//determine surface-density profile					
						if (sdp) {
							//double rmax = 50.0;//pc
							//double rmin = 0.025;//pc
							double rmax = 500.0;//pc
							double rmin = 0.1;//pc
							double stepsize;
							stepsize = (log10(rmax)-log10(rmin))/(numberradialbins-1);
							
							for (l=0;l<numberradialbins;l++) {
								profile[l][0] = pow(10.0, log10(rmin) + stepsize*l); //radius
								//profile[l][0] = 1.0*exp(log(rmax-rmin)*(l+1)/numberradialbins)-1.0*exp(log(rmax-rmin)/numberradialbins)+rmin; //radius
								//profile[l][0] = 1.0*pow(base,log10(rmax-rmin)/log10(base)*(l+1)/numberradialbins)-1.0*pow(base,log10(rmax-rmin)/log10(base)/numberradialbins)+rmin; //radius
								//profile[l][0] = 1.0*pow(10.0,0.25*l-2.0); //radius
								profile[l][2] = 0; //mass
								profile[l][3] = 0; //mv
								profile[l][4] = 0; //m(v^2)
								profile[l][5] = 0; //number
								profile[l][6] = 0; //error of velocity dispersion
								profile[l][7] = 0; //luminosity
							}
							for (l=0;l<numberradialbins;l++) {
								if (l == 0)
									profile[l][1] =  2.0*PI*pow(profile[l][0],2);
								else
									profile[l][1] = 2.0*PI*(pow(profile[l][0],2) - pow(profile[l-1][0],2)); //area
							}
							
							for (l=0;l<numberradialbins;l++) {
								profile_inside_rtide[l][0] = pow(10.0, log10(rmin) + stepsize*l); //radius
								profile_inside_rtide[l][2] = 0; //mass
								profile_inside_rtide[l][3] = 0; //mv
								profile_inside_rtide[l][4] = 0; //m(v^2)
								profile_inside_rtide[l][5] = 0; //number
								profile_inside_rtide[l][6] = 0; //error of velocity dispersion
								profile_inside_rtide[l][7] = 0; //luminosity
							}
							for (l=0;l<numberradialbins;l++) {
								if (l == 0)
									profile_inside_rtide[l][1] =  2.0*PI*pow(profile_inside_rtide[l][0],2);
								else
									profile_inside_rtide[l][1] = 2.0*PI*(pow(profile_inside_rtide[l][0],2) - pow(profile_inside_rtide[l-1][0],2)); //area
							}
							
							for (l=0;l<numberradialbins;l++) {
								profile_bound_stars[l][0] = pow(10.0, log10(rmin) + stepsize*l); //radius
								//profile_bound_stars[l][0] = 1.0*pow(10.0,log10(rmax)*(l+1)/numberradialbins); //radius
								profile_bound_stars[l][2] = 0; //mass
								profile_bound_stars[l][3] = 0; //mv
								profile_bound_stars[l][4] = 0; //m(v^2)
								profile_bound_stars[l][5] = 0; //number
								profile_bound_stars[l][6] = 0; //error of velocity dispersion
								profile_bound_stars[l][7] = 0; //luminosity
							}
							for (l=0;l<numberradialbins;l++) {
								if (l == 0)
									profile_bound_stars[l][1] =  2.0*PI*pow(profile_bound_stars[l][0],2);
								else
									profile_bound_stars[l][1] = 2.0*PI*(pow(profile_bound_stars[l][0],2) - pow(profile_bound_stars[l-1][0],2)); //area
							}
							
							surfacedensityprofile(Ninit, star, profile, mx, my, mz,mvx,mvy,mvz, view, numberradialbins, VG, profile_inside_rtide, profile_bound_stars);
						}					
						
						
						//determine escape velocities
						if (trackescapers)
							escape(esc, Nmax, star, rtide, mx, my, mz, mvx, mvy, mvz, omega, time, sigma, kT, M, outputtype, timebin, timebinsize, threshold);
						
						
						//stop analysis?
						if (M<Mend) lSize = 0;
						
						
						//output to STDOUT
						if (outputtype == 0) {
							printf("\n---------------------------------------------------------------------------\n");
							printf("T = %g\tM = %g\tkT = %g\tsigma = %g\tfbin = %g\n",time,M,kT,sigma,fbin);
							printf("Rh = %g\tRtide = %g (d2phirdr2 = %g)\tRc = %g\n",Rh,rtide,-1.0*dphi,Rc);
							printf("Ebin = %g\tEkin = %g\tEtot = %g\tL = %g\n",E[0],E[1],E[2],L);
						}
						
						
						
						/************
						 *  OUTPUT  *
						 ************/
						
						FILE *out;
						char output[200];
						double xtemp[3], vtemp[3], dsuntemp, vrsuntemp, vrtemp, lcosbtemp, RAtemp, DECtemp, mu_alphatemp, mu_alphacosdeltatemp, mu_deltatemp, mutemp, PAtemp;
						
						
						//write summary to .new
						fprintf(new,"%g\t%g\t%i\t%i\t%i\t%i\t%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",time,M,N,n[0],n[1],n[2],n[3],nbin[0],fbin,sigma,sigmauncor,E[0],E[1],E[2],Euncor[0],Euncor[1],Euncor[2],kT,LR[0],LR[1],LR[2],LR[3],LR[4],LR[5],Rc,rtide,omega,dphi,L, mv,KR[0],KR[1],KR[2],KR[3],KR[4],KR[5],KR[6],KR[7],KR[8],KR[9],racc);
						
						
						//write surface-density profile to .sdp
						if (sdp) {
							sprintf(output, "%s/%s%07i.sdp",outputfilename,outputfilename,timetemp);
							out = fopen(output,"w");
							for(l=0;l<numberradialbins;l++) {
								fprintf(out,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",profile[l][0],1.0*profile[l][5]/profile[l][1],profile[l][0]/rtide,profile[l][3],profile[l][4],profile[l][5],profile[l][6],profile_inside_rtide[l][0],1.0*profile_inside_rtide[l][5]/profile_inside_rtide[l][1],profile_inside_rtide[l][0]/rtide,profile_inside_rtide[l][3],profile_inside_rtide[l][4],profile_inside_rtide[l][5],profile_inside_rtide[l][6],profile_bound_stars[l][0],1.0*profile_bound_stars[l][5]/profile_bound_stars[l][1],profile_bound_stars[l][0]/rtide,profile_bound_stars[l][3],profile_bound_stars[l][4],profile_bound_stars[l][5],profile_bound_stars[l][6], 1.0*profile[l][7]/profile[l][1], 1.0*profile_inside_rtide[l][7]/profile[l][1], 1.0*profile_bound_stars[l][7]/profile[l][1], 1.0*profile[l][2]/profile[l][1], 1.0*profile_inside_rtide[l][2]/profile[l][1], 1.0*profile_bound_stars[l][2]/profile[l][1]);
							}
							fclose(out);
							printf("Written to: %s\n\n",output);
						}
						
						
						//write star dumps to all.txt and all.bin
						if ((alloutput) && (timetemp % outputstep < 1)) {
							sprintf(output, "%s/%s%07iall.txt",outputfilename,outputfilename,timetemp);
							out = fopen(output,"w");
							/*						
                             fprintf(out,"# rtide = %g\n",rtide);
							 fprintf(out,"# omega = %g\n",omega);
							 fprintf(out,"# M = %g\n",M);
							 fprintf(out,"# time = %g\n",time);
							 fprintf(out,"# mx = %g\n",mx);
							 fprintf(out,"# my = %g\n",my);
							 fprintf(out,"# mz = %g\n",mz);
							 fprintf(out,"# RGx = %g\n",RG[0]);
							 fprintf(out,"# RGy = %g\n",RG[1]);
							 fprintf(out,"# RGz = %g\n",RG[2]);
                             */
							
							double ybinsize = 20.0;   //in pc [20]
							int yrange = 4000, ybins, yt;   //y in pc [4000]
							int xrange = 4000, xbins, xt;   //x in pc
							ybins = 1.0*yrange/ybinsize;
							xbins = 1.0*xrange/ybinsize;
							
							double **cont;
							cont = (double **)calloc(xbins,sizeof(double *));
							for (i=0;i<xbins;i++){
								cont[i] = (double *)calloc(ybins,sizeof(double));
								if (cont[i] == NULL) {
									printf("\nNULL; memory allocation failed!\n");
									return 0;
								}
							}
							
							//writing all stars to all.txt
							double lummean = 0.0;
							int Nlum = 0;
							for(l=0;l<Ntot;l++) {
								if ((star[l][0]) && (star[l][13]>=0.0)) {
									//fprintf(out,"%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f\n",star[l][0],star[l][1],star[l][2],star[l][3],star[l][4],star[l][5],star[l][6],star[l][7],star[l][8],star[l][9],star[l][10],star[l][11],star[l][12],star[l][13]);
									fprintf(out,"%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",star[l][0],star[l][1]-mx,star[l][2]-my,star[l][3]-mz,star[l][4]-mvx,star[l][5]-mvy,star[l][6]-mvz,star[l][7],star[l][8],star[l][9],star[l][10],star[l][11],star[l][12],star[l][13]);
									if (star[l][7]<10) {
										lummean += star[l][8];
										Nlum++;
									}
									//Binning in x-y-plane
									xt = 1.0*(star[l][1]+0.5*xrange)/ybinsize;
									yt = 1.0*(star[l][2]+0.5*yrange)/ybinsize;//use star[l][3] for perpendicular-to-disk orbits
									if ((yt<ybins) && (yt>=0) && ((xt<xbins) && (xt>=0))) {
										cont[xt][yt]++;
									}
								}
							}
							lummean /= Nlum;
							//printf("\nL_mean = %f\n",lummean);
							fclose(out);
							printf("\nWritten to: %s\n",output);
							
							
							//writing bins to all.bin
							sprintf(output, "%s/%s%07iall.bin",outputfilename,outputfilename,timetemp);
							out = fopen(output,"w");
							int ll;
							for(l=0;l<xbins;l++) {
								for(ll=0;ll<ybins;ll++) {
									fprintf(out,"%f\t%f\t%f\n",(l+0.5)*ybinsize-0.5*xrange,(ll+0.5)*ybinsize-0.5*yrange,cont[l][ll]);
								}
								fprintf(out,"\n");
							}
							
							for (i=0;i<xbins;i++) free(cont[i]);
							free(cont);
							fclose(out);
							printf("Written to: %s\n\n",output);
						}	
						
						
						//write tails to y.txt and y.bin	
						if ((tailanalysis) && (timetemp % outputstep < 1)) {
							
							int pcorbitlength = 4000; // +/- pc
							int pcstepsize = 25.0; //pc
							int pcnumberintervals;
							pcnumberintervals = 1*(pcorbitlength/pcstepsize+1);
							double **along_the_orbit;
							along_the_orbit = (double **)calloc(pcnumberintervals,sizeof(double *));
							for (i=0;i<pcnumberintervals;i++){
								along_the_orbit[i] = (double *)calloc(8,sizeof(double));
								if (along_the_orbit[i] == NULL) {
									printf("\nNULL; memory allocation failed!\n");
									return 0;
								}
							}
							
							//calculate orbit forwards
							int sign = 1;
							for (i=0;i<pcnumberintervals;i++)
								along_the_orbit[i][0] = sign*i*pcstepsize;
							
							double RGcorr[3], VGcorr[3];
							RGcorr[0] = RG[0]+mx;
							RGcorr[1] = RG[1]+my;
							RGcorr[2] = RG[2]+mz;
							VGcorr[0] = VG[0]+mvx;
							VGcorr[1] = VG[1]+mvy;
							VGcorr[2] = VG[2]+mvz;
							
							get_orbit(RGcorr,VGcorr,along_the_orbit,pcnumberintervals,sign);
							
							double orbitbin[pcnumberintervals*2-3][11];
							
							for (i=0;i<pcnumberintervals-1;i++) {
								orbitbin[pcnumberintervals+i-2][0] = along_the_orbit[i][0];
								orbitbin[pcnumberintervals+i-2][1] = along_the_orbit[i][1];
								orbitbin[pcnumberintervals+i-2][2] = along_the_orbit[i][2];
								orbitbin[pcnumberintervals+i-2][3] = along_the_orbit[i][3];
								orbitbin[pcnumberintervals+i-2][4] = 0;
								orbitbin[pcnumberintervals+i-2][5] = along_the_orbit[i][5];
								orbitbin[pcnumberintervals+i-2][6] = along_the_orbit[i][6];
								orbitbin[pcnumberintervals+i-2][7] = along_the_orbit[i][7];
								orbitbin[pcnumberintervals+i-2][8] = 0;
								orbitbin[pcnumberintervals+i-2][9] = 0;
								orbitbin[pcnumberintervals+i-2][10] = 0;
							}
							
							//calculate orbit backwards
							sign = -1;
							for (i=0;i<pcnumberintervals;i++)
								along_the_orbit[i][0] = sign*i*pcstepsize;
							
							get_orbit(RGcorr,VGcorr,along_the_orbit,pcnumberintervals,sign);
							
							for (i=0;i<pcnumberintervals-1;i++) {
								orbitbin[i][0] = along_the_orbit[pcnumberintervals-i-2][0];
								orbitbin[i][1] = along_the_orbit[pcnumberintervals-i-2][1];
								orbitbin[i][2] = along_the_orbit[pcnumberintervals-i-2][2];
								orbitbin[i][3] = along_the_orbit[pcnumberintervals-i-2][3];
								orbitbin[i][4] = 0;
								orbitbin[i][5] = along_the_orbit[pcnumberintervals-i-2][5];
								orbitbin[i][6] = along_the_orbit[pcnumberintervals-i-2][6];
								orbitbin[i][7] = along_the_orbit[pcnumberintervals-i-2][7];
								orbitbin[i][8] = 0;
								orbitbin[i][9] = 0;
								orbitbin[i][10] = 0;
							}
							
							
							//determine next orbit point for each star
							double rtemp;
							double rmin[2];
							for(l=0;l<Ntot;l++) {
								if ((star[l][0]) && (star[l][13]>=0.0) && (star[l][10]>=limit)) {
									rmin[0] = 750.0;
									rmin[1] = 0;
									for(i=0;i<pcnumberintervals*2-3;i++) {
										rtemp = sqrt(pow(star[l][1]+RG[0]-orbitbin[i][1],2)+pow(star[l][2]+RG[1]-orbitbin[i][2],2)+pow(star[l][3]+RG[2]-orbitbin[i][3],2));
										if (rtemp < rmin[0]) {
											rmin[0] = rtemp;
											rmin[1] = i;
										}	
									}
									orbitbin[(int) rmin[1]][4] += 1.0;
									
									orbitbin[(int) rmin[1]][8] += sqrt(pow(star[l][4]+VG[0],2)+pow(star[l][5]+VG[1],2)+pow(star[l][6]+VG[2],2));
									
									orbitbin[(int) rmin[1]][9] += pow(star[l][4]+VG[0],2)+pow(star[l][5]+VG[1],2)+pow(star[l][6]+VG[2],2);
									
								}
							}
							for (i=0;i<pcnumberintervals*2-3;i++) {
								if (orbitbin[i][4]) {
									orbitbin[i][8]/=orbitbin[i][4];
									orbitbin[i][10] = orbitbin[i][8];
									
									orbitbin[i][8] -= sqrt(pow(orbitbin[i][5],2)+pow(orbitbin[i][6],2)+pow(orbitbin[i][7],2));//-orbitbin[i][8];
									
									orbitbin[i][9]/=orbitbin[i][4];
								} else {
									orbitbin[i][8] = 0.01;
									orbitbin[i][9] = 0.0;
									orbitbin[i][10] = 0.0;
								}
							}
							
							
							/*						for (i=0;i<pcnumberintervals*2-3;i++)
							 printf("%f\t%f\t%f\t%f\t%f\t%f\n",orbitbin[i][0],orbitbin[i][1],orbitbin[i][2],orbitbin[i][3],sqrt(pow(orbitbin[i][1]-orbitbin[0][1],2)+pow(orbitbin[i][2]-orbitbin[0][2],2)+pow(orbitbin[i][3]-orbitbin[0][3],2)),orbitbin[i][4]);*/
							
							//write bins to y.bin
							sprintf(output, "%s/%s%07iy.bin",outputfilename,outputfilename,timetemp);
							out = fopen(output,"w");
							for(i=0;i<pcnumberintervals*2-3;i++) {
								fprintf(out,"%f\t%f\t%f\t%f\t%f\n",orbitbin[i][0],orbitbin[i][4],orbitbin[i][8],orbitbin[i][9],orbitbin[i][10]);
							}
							
							for (i=0;i<pcnumberintervals;i++) free(along_the_orbit[i]);
							free(along_the_orbit);
							fclose(out);
							
							printf("Written to: %s\n\n",output);
						}
						
						
						//write snapshot from the sun to sun.txt
						if (timetemp % sunviewtime < 1) {
							FILE *out2, *out3;
							char output2[200], output3[200];
							sprintf(output, "%s/%s%07isun.txt",outputfilename,outputfilename,timetemp);
							out = fopen(output,"w");
							
							sprintf(output2, "%s/%s%07isun.bin",outputfilename,outputfilename,timetemp);
							out2 = fopen(output2,"w");
							
							sprintf(output3, "%s/%s%07isunvr.bin",outputfilename,outputfilename,timetemp);
							out3 = fopen(output3,"w");
							
							
							//get coordinates of cluster center
							xgal = RG[0]+rgalsun+mx;
							ygal = RG[1]+my;
							zgal = RG[2]+mz;
							
							xtemp[0] = RG[0]+mx;
							xtemp[1] = RG[1]+my;
							xtemp[2] = RG[2]+mz;
							vtemp[0] = VG[0];
							vtemp[1] = VG[1];
							vtemp[2] = VG[2];
							convert(xtemp, vtemp, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSR);
							
							Rtemp = dsuntemp;
							vr = vrsuntemp;
							
							if (ltemp>180.0) {
								ltemp -= 360.0;
								lcosbtemp = ltemp*cos(btemp);
							}
                            
                            
							//set telescope specificaltions
							//int fov = 64;      //Field of view in deg
							double fov = 4.0;      //Field of view in deg
							double fovh;
							fovh = 0.5*fov;
							//double ybinsize = 1.0*fov/320.0;   //in deg
							double ybinsize = 1.0*fov/80.0;   //in deg
							int ybins, yt;   //b in deg
							int xbins, xt;   //lcosb in deg
							int vrbins, vrt;   //vr in km/s
							double vrbinsize = 1.0;   //in km/s
							ybins = (360.0+fov)/ybinsize;        //full sky
							xbins = (360.0+fov)/ybinsize;        //full sky
							vrbins = 100.0/vrbinsize;
							int xmaxt, xmint, ymaxt, ymint, vrmaxt, vrmint;
							
							xmint = 1.0*(ltemp-fovh+180.0)/ybinsize;
							xmaxt = 1.0*(ltemp+fovh+180.0)/ybinsize;
							ymaxt = 1.0*(btemp+fovh+180.0)/ybinsize;
							ymint = 1.0*(btemp-fovh+180.0)/ybinsize;
							vrmaxt = 100.0/vrbinsize;
							vrmint = 0.0;
														
							double **cont;
							cont = (double **)calloc(xbins,sizeof(double *));
							for (i=0;i<xbins;i++){
								cont[i] = (double *)calloc(ybins,sizeof(double));
								if (cont[i] == NULL) {
									printf("\nNULL; memory allocation failed!\n");
									return 0;
								}
							}
							
							double **contvr;
							contvr = (double **)calloc(xbins,sizeof(double *));
							for (i=0;i<xbins;i++){
								contvr[i] = (double *)calloc(vrbins,sizeof(double));
								if (contvr[i] == NULL) {
									printf("\nNULL; memory allocation failed!\n");
									return 0;
								}
							}
							
							//for Pal5: M(R<4'), M(R<12') and M(R<16') as in Dehnen et al. (2004)
							double M4pal5 = 0.0, M12pal5 = 0.0, M16pal5 = 0.0;
							double lgal, bgal, rtemppal5, drdeg;
							lgal = ltemp;
							bgal = btemp;
							double x1t, x2t, y1t, y2t, z1t, z2t;
							x1t = xgal;
							y1t = ygal;
							z1t = zgal;
							
							double Teff, lTeff, BV, BC, kb,vmag,abvmag;
							double bvc[8], bcc[8];
							double dvmag, dbmag, dBV;
							double BCsun, abvmagsun;
							
							kb = 5.6704E-08*0.5*1.3914E9*0.5*1.3914E9/3.846E26;  //Stefan-Boltzmann constant in Lsun Rsun^-2 K^-4
							
							bvc[0] = -654597.405559323;
							bvc[1] = 1099118.61158915;
							bvc[2] = -789665.995692672;
							bvc[3] = 314714.220932623;
							bvc[4] = -75148.4728506455;
							bvc[5] = 10751.803394526;
							bvc[6] = -853.487897283685;
							bvc[7] = 28.9988730655392;
							
							bcc[0] = -4222907.80590972;
							bcc[1] = 7209333.13326442;
							bcc[2] = -5267167.04593882;
							bcc[3] = 2134724.55938336;
							bcc[4] = -518317.954642773;
							bcc[5] = 75392.2372207101;
							bcc[6] = -6082.7301194776;
							bcc[7] = 209.990478646363;
							
							BCsun = 0.11;  //sun's bolometric correction
							abvmagsun = 4.83; //sun's absolute V magnitude
							
							//write all stellar coordinates to sun.txt
							for(l=0;l<Ntot;l++) {
								if ((star[l][0]) && (star[l][13]>=0.0) && (star[l][10]>=limit)) {
									xtemp[0] = RG[0]+star[l][1];
									xtemp[1] = RG[1]+star[l][2];
									xtemp[2] = RG[2]+star[l][3];
									vtemp[0] = VG[0]+star[l][4];
									vtemp[1] = VG[1]+star[l][5];
									vtemp[2] = VG[2]+star[l][6];
									convert(xtemp, vtemp, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSR);

									if (ltemp>180.0) {
										ltemp -=360.0;
										lcosbtemp = ltemp*cos(btemp);
									}
									
									xgal = star[l][1]+RG[0]+rgalsun;
									ygal = star[l][2]+RG[1];
									zgal = star[l][3]+RG[2];
									vr = vrsuntemp;
									Rtemp = dsuntemp;
									
									
									//for Pal5: M(R<4'), M(R<12') and M(R<16') as in Dehnen et al. (2004)
									rtemppal5 = sqrt(pow(ltemp-lgal,2)+pow(btemp-bgal,2));
									if (rtemppal5 < 4.0*0.00265) M4pal5+=star[l][10];
									if (rtemppal5 < 12.0*0.00265) M12pal5+=star[l][10];
									if (rtemppal5 < 16.0*0.00265) M16pal5+=star[l][10];
									
									
									//get radial distance of each star from cluster center
									//drdeg = (sqrt(pow(ltemp-lgal,2) + pow(btemp-bgal,2)))*360.0/(2.0*PI);
									x2t = xgal;
									y2t = ygal;
									z2t = zgal;
									drdeg = acos((x1t*x2t+y1t*y2t+z1t*z2t)/(sqrt(x1t*x1t+y1t*y1t+z1t*z1t)*sqrt(x2t*x2t+y2t*y2t+z2t*z2t)));
									
									
									//get CMD
									if (star[l][9] && star[l][8]>=0.0) {
										Teff = pow(star[l][8]/(4.0*PI*star[l][9]*star[l][9]*kb),0.25);
										if ((Teff>3000.0) && (Teff<55000.0)) {
											
											lTeff = log10(Teff);
											
											BV = bvc[0] + bvc[1]*lTeff + bvc[2]*pow(lTeff,2) + bvc[3]*pow(lTeff,3) + bvc[4]*pow(lTeff,4) + bvc[5]*pow(lTeff,5) + bvc[6]*pow(lTeff,6) + bvc[7]*pow(lTeff,7);
											
											BC = bcc[0] + bcc[1]*lTeff + bcc[2]*pow(lTeff,2) + bcc[3]*pow(lTeff,3) + bcc[4]*pow(lTeff,4) + bcc[5]*pow(lTeff,5) + bcc[6]*pow(lTeff,6) + bcc[7]*pow(lTeff,7);
											
											if (star[l][8]) abvmag = -2.5*log10(star[l][8])-BC+BCsun+abvmagsun;
											
											vmag = abvmag + 5.0*log10(Rtemp) - 5.0;
											
											double rand1, rand2, prand;
											do {
												rand1 = 2.0*drand48()-1.0;
												rand2 = 2.0*drand48()-1.0;
											} while (rand1*rand1+rand2*rand2 > 1.0);
											prand = sqrt(-2.0*log(rand1*rand1+rand2*rand2)/(rand1*rand1+rand2*rand2));
											dvmag = rand1*prand*sqrt(pow(0.02,2) + pow(0.07*pow(10.0, 0.4*(vmag-25.0)),2));
											dbmag = rand2*prand*sqrt(pow(0.02,2) + pow(0.07*pow(10.0, 0.4*(vmag-25.0)),2));
											dBV = dvmag + dbmag;//sqrt(pow(dvmag,2) + pow(dbmag,2));
											
										} else {
											vmag = 9999.9;
											abvmag = 9999.9;
											BV = 9999.9;
											dvmag = 0.0;
											dBV = 0.0;
										}
									} else {
										Teff = 0.0;
										vmag = 9999.9;
										abvmag = 9999.9;
										BV = 9999.9;
										dvmag = 0.0;
										dBV = 0.0;
									}
																		
									
									if (vmag < vlimit) {
										fprintf(out,"%.0f\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.0lf\t%.16lf\n",star[l][0],Rtemp,ltemp,btemp,vmag,abvmag,BV,vr,star[l][8],dvmag,star[l][10],dBV,Teff,star[l][13],drdeg);
										
										//Binning in observer coordinates
										xt = 1.0*(ltemp+180.0)/ybinsize;
										yt = 1.0*(btemp+180.0)/ybinsize;
										vrt = 1.0*(vr+100.0)/vrbinsize;
										
										//printf ("%i\t%i\t%i\n\n",xt, yt, vrt);
										if ((yt<ymaxt) && (yt>=ymint) && (xt<xmaxt) && (xt>=xmint)) {
											cont[xt][yt] = cont[xt][yt]+1.0;
										}
										if ((vrt<vrmaxt) && (vrt>=vrmint) && (xt<xmaxt) && (xt>=xmint)) {
											contvr[xt][vrt] = contvr[xt][vrt]+1.0;
										}
									}
								}
							}
							
							fclose(out);
							printf("\nWritten to: %s\n",output);
							
                            
							//writing all bins in observer coordinates to sun.bin
							int ll;
							for(l=xmint;l<xmaxt;l++) {
								for(ll=ymint;ll<ymaxt;ll++) {
									fprintf(out2,"%f\t%f\t%lf\n",(l+0.5)*ybinsize-180.0,(ll+0.5)*ybinsize-180.0,cont[l][ll]);
								}
								fprintf(out2,"\n");
							}
							
							fclose(out2);
							for(l=xmint;l<xmaxt;l++) {
								for(ll=vrmint;ll<vrmaxt;ll++) {
									fprintf(out3,"%f\t%f\t%lf\n",(l+0.5)*ybinsize-180.0,(ll+0.5)*vrbinsize-100.0,contvr[l][ll]);
								}
								fprintf(out3,"\n");
							}
							
							fclose(out3);
							for (l=0;l<xbins;l++) free(cont[l]);
							free(cont);
							for (l=0;l<xbins;l++) free(contvr[l]);
							free(contvr);
							
							printf("Written to: %s and %s\n\n",output2, output3);
						}
						
						
						//write galactic orbit to .gal
						xtemp[0] = RG[0];
						xtemp[1] = RG[1];
						xtemp[2] = RG[2];
						vtemp[0] = VG[0];
						vtemp[1] = VG[1];
						vtemp[2] = VG[2];
						
						convert(xtemp, vtemp, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSR);
						
						xgal = RG[0]+rgalsun;
						ygal = RG[1];
						zgal = RG[2];
						Rtemp = dsuntemp;
						vr = vrsuntemp;
						
						fprintf(gal,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",time,RG[0],RG[1],RG[2],VG[0],VG[1],VG[2],Rtemp,ltemp,btemp,vr);
						
					}
					
					step++;
					lPos = ftell(fz);

					if (outputtype) printf("%li bytes of %li\n",lPos,lSize);
				} while(lPos+20<lSize);
				
				
			} else {
				printf("File %s is empty!\n",filename);
				printf("##########################################################\n");
			}
			fclose(fz);
			poscountertemp++;
			sprintf(filename, "%s.POS%i", argv[q],poscountertemp);
			fz = fopen(filename,"rb");
			outputfilename = strtok(filename, ".");
		} while (poscountertemp <= poscounter);
		
		
		//write dynamical temperature-luminosity data to .tld
		if (trackescapers) {
			FILE *tld;	
			char tldfile[200];
			sprintf(tldfile, "%s.tld",outputfilename);
			tld = fopen(tldfile,"w");
			
			i = 20000/timebinsize;
			for (l=0;l<i;l++) {
				if (timebin[l][0]) {
					timebin[l][1] /= timebin[l][0]; //0.5*velocity dispersion
					timebin[l][2] /= timebinsize; //evaporation luminosity
					timebin[l][3] /= timebinsize; //ejection luminosity
					timebin[l][4] /= timebin[l][0]; //total mass
					fprintf(tld,"%f\t%f\t%f\t%f\t%f\t%f\n",(l+0.5)*timebinsize,timebin[l][0],timebin[l][1],timebin[l][2],timebin[l][3],timebin[l][4]);
				}
			}
			fclose(tld);
		}
		
		for (l=0;l<numberradialbins;l++) free (profile[l]);
		free(profile);
		for (j=0;j<Nmax;j++) free (star[j]);
		free(star);
		for (j=0;j<Nmax;j++) free(neighborlist[j]);
		free(neighborlist);	
		
		if (trackescapers) fclose(esc);
		fclose(new);
		fclose(gal);
		if (binoutput) fclose(bin);
		
		printf("\nSummary written to: %s\n\n",summary);
		
	}
	
	t2 = clock();             
	printf("\nElapsed time: %g sec\n",(double)(t2-t1)/CLOCKS_PER_SEC); //write elapsed time to STDOUT
	return 1;
}






/***************
 * SUBROUTINES *
 ***************/

int escape(FILE *esc, int Nt, double **start, double rtidet, double mxt, double myt, double mzt, double mvxt, double mvyt, double mvzt, double omegat, double timet, double sigmat, double kTt, double Mt, int outputtypet, double **timebint, int timebinsizet, double thresholdt){
	int l,p,i;
	double vx, vy, vz, v, x, y, z, r, dr;
	double dt, Fx, Fy, Fz;
	double xt, yt, zt, timett,drt;//,vt,rt,drt;
	double vr, cosa, ratio;
	float dtlist[200][2];
	
	for (l=0;l<Nt;l++) {		
		if ((start[l][0]) && (start[l][14] < 2.0) && (start[l][14])) {
			
			x = start[l][1]-mxt;
			y = start[l][2]-myt;
			z = start[l][3]-mzt;
			r = sqrt(x*x+y*y+z*z);
			dr = r - 2.0*rtidet;
			vx = start[l][4]-mvxt;
			vy = start[l][5]-mvyt;
			vz = start[l][6]-mvzt;
			v = sqrt(vx*vx+vy*vy+vz*vz);
			timett = timet;
			cosa = (x*vx+y*vy+z*vz)/(r*v);
			vr = v*cosa;
			
			if (outputtypet) printf("   ESCAPE: %.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",timet,start[l][0],r,dr,v,vr,start[l][10],sigmat,kTt);
			
			
			//Initialise integration
			dt = -sqrt(dr/vr*dr/vr);   //dt in Myrs
			p = 0;
			drt = dr;
			if (timett > sqrt(dt*dt)) {
				dtlist[p][0] = sqrt(drt*drt);
				dtlist[p][1] = dt;
			} else {
				do {
					dt = -1.0*sqrt(dr/vr*dr/vr)*drand48();
				} while (timett < sqrt(dt*dt));
				dtlist[p][0] = sqrt(drt*drt);
				dtlist[p][1] = dt;
			}
			Fx = start[l][10]*(G*Mt*x/(r*r*r) + 2.0*omegat*vy + 3.0*omegat*omegat*x);
			Fy = start[l][10]*(G*Mt*y/(r*r*r) - 2.0*omegat*vx);
			Fz = start[l][10]*(G*Mt*z/(r*r*r) + omegat*omegat*z);
			
			
			while ((drt > 0.05*2.0*rtidet) || (drt < -0.05*2.0*rtidet)) {
				p++;
				if (p==100) {
					do {
						dt = -1.0*sqrt(dr/vr*dr/vr)*drand48();
					} while (timett < sqrt(dt*dt));
				}
				
				xt = x+vx*dt+0.5*Fx*dt*dt;
				yt = y+vy*dt+0.5*Fy*dt*dt;
				zt = z+vz*dt+0.5*Fz*dt*dt;
				//printf(" %i: %f ",p,dt);
				
				if (p<200) {
					drt = sqrt(xt*xt+yt*yt+zt*zt)-2.0*rtidet;
					if (timett > sqrt(dt*dt)) {
						dtlist[p][0] = sqrt(drt*drt);
						dtlist[p][1] = dt;
					} else {
						dtlist[p][0] = 1000000.0;
						dtlist[p][1] = 0.0;
					}
				} else {
					drt = 0.0;
					if (outputtypet) printf("CORRECTION FAILED (%.0f)\t",start[l][0]);
					qsort(dtlist, 200, sizeof(dtlist[0]), (void *)cmpmy);
				}
				
				if (drt<0) {
					dt *= 0.99;
				} else {
					dt *= 1.01;
				}
			}
			
			if (p){
				
				if (drt<0) {
					dt /= 0.99;
				} else {
					dt /= 1.01;
				}
				
				if (p==200) {
					dt = dtlist[0][1];
					drt = dtlist[0][0];
					if (outputtypet) printf("DTLIST: dr = %f\tdt = %f\n",drt,dt);
				}
				
				vx += Fx*dt;
				vy += Fy*dt;
				vz += Fz*dt;
				
				r = sqrt(xt*xt+yt*yt+zt*zt);
				v = sqrt(vx*vx+vy*vy+vz*vz);
				
				cosa = (xt*vx+yt*vy+zt*vz)/(r*v);
				vr = v*cosa;
				
				timett = timet+dt;
				if (outputtypet) printf("CORRECTED: %.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",timett,start[l][0],r,drt,v,vr,start[l][10],sigmat,kTt);
			}
			
			fprintf(esc,"%.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",timett,start[l][0],r,drt,v,vr,start[l][10],sigmat,kTt);
			
			i = timett/timebinsizet;
			if ((i >= 0) && (i<20000/timebinsizet) ) {
				timebint[i][0] += 1.0;
				timebint[i][1] += 0.5*sigmat;
				timebint[i][4] += Mt;
				
				ratio = v*v/sigmat;
				
				if (ratio < thresholdt) 
					timebint[i][2] += 0.5*start[l][10]*v*v;
				else
					timebint[i][3] += 0.5*start[l][10]*v*v;
			}
			
			start[l][14]=2.0;
			
			
		}
	}	
	
	
	
	return 1;
}

int rcsuche(int Nt, double **start, int **neighborlistt, double mxt, double myt, double mzt, int densitylevelt, double *Rct){ 
	int l, ll;
	double mtemp, r2temp, rho2temp, rho2ges = 0.0, rho2malr2ges = 0.0, rho2malr2temp;
	
	for (l=0;l<Nt;l++) {
		if ((start[l][0]) && (start[l][13]<3)) {
			mtemp = 0.0;
			for (ll=0;ll<densitylevelt;ll++) {
				mtemp += start[neighborlistt[l][ll]][10];
			}
			rho2temp = pow(mtemp/pow(pow(start[neighborlistt[l][densitylevelt]][1]-start[l][1],2)+pow(start[neighborlistt[l][densitylevelt]][2]-start[l][2],2)+pow(start[neighborlistt[l][densitylevelt]][3]-start[l][3],2),1.5),2);
			r2temp = pow(mxt-start[l][1],2)+pow(myt-start[l][2],2)+pow(mzt-start[l][3],2);
			rho2malr2temp = rho2temp*r2temp;
		} else {
			rho2malr2temp = 0.0;
			rho2temp = 0.0;
		}
		
		rho2ges += rho2temp; 
		rho2malr2ges += rho2malr2temp; 
	}
	
	*Rct = sqrt(rho2malr2ges/rho2ges);
	
	return 1;
}

int get_orbit(double *RGt,double *VGt,double **along_the_orbitt,int pcnumberintervalst, int signt){
	double x[3], v[3];
	int k;
	
	x[0] = RGt[0];
	x[1] = RGt[1];
	x[2] = RGt[2];
	v[0] = VGt[0];
	v[1] = VGt[1];
	v[2] = VGt[2];
	
	double t = 0.0;
	double pc = 0.0;
	double pcmax = along_the_orbitt[pcnumberintervalst-1][0];
	double dpcout = along_the_orbitt[1][0];
	double mdiff = 1.E-11;         //Genauigkeit
	double dtout;
	dtout = 0.01*dpcout/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	
	double pcout,diff,dt,tout;
	double *xe1,*xe2,*ve1,*ve2;
	
	
	
	xe1 = malloc(3*sizeof(double));
	xe2 = malloc(3*sizeof(double));
	ve1 = malloc(3*sizeof(double));
	ve2 = malloc(3*sizeof(double));
	
	pcout = 0.0;                   /* Distanzgrenze fuer die naechste Ausgabe */
	tout = 0.0;
	
	dt = signt*0.00001;                    /* Initial Condition fuer die Zeitschritte */
	
	int step = 0;
	
	while (signt * pc <= signt * pcmax) {
		if (signt * pc >= signt * pcout) {              /* Ausgabe der Werte */
			along_the_orbitt[step][1] = x[0];
			along_the_orbitt[step][2] = x[1];
			along_the_orbitt[step][3] = x[2];
			along_the_orbitt[step][4] = 0;
			along_the_orbitt[step][5] = v[0];
			along_the_orbitt[step][6] = v[1];
			along_the_orbitt[step][7] = v[2];
			step++;
			pcout = along_the_orbitt[step][0];              /* Distanzgrenze fuer Ausgabe wird erhoeht */
		}
		
		do {
			for (k=0;k<3;k++) {
				xe1[k]=x[k];
				xe2[k]=x[k];
				ve1[k]=v[k];
				ve2[k]=v[k];
			}
			do_step(dt,xe1,ve1);      /* Ein ganzer Zeitschritt */
			
			do_step(0.5*dt,xe2,ve2);  /* Zwei halbe Zeitschritte */
			do_step(0.5*dt,xe2,ve2);
			
			diff = sqrt(pow(*xe1 - *xe2,2) + pow(*(xe1+1) - *(xe2+1),2) + pow(*(xe1+2) - *(xe2+2),2));
			
			if (diff<mdiff) {         /* Pruefung, ob Differenz der beiden unter Genauigkeitsgrenze liegt */
				t+=dt;
				for (k=0;k<3;k++) {
					x[k]=xe2[k];          /* Pruefung erfolgreich -> Werte weitergeben und Zeitschritt vergroessern */
					v[k]=ve2[k];
				}
				pc+=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*dt;
			}
			dt = dt*pow(mdiff/diff,0.2); /* Verkleinerung des Zeitschritts */
			
			if ((signt*dt < 0.0000000001) ) {
				printf("Abbruch... dt = %lf\n", dt);
				pc = pcmax*2.0;
				diff = mdiff/2;
			}
			
			if (signt * dt > signt * dtout) {
				dt=dtout;
			}
			
		} while (diff>mdiff);          /* Schleife wird erst einmal ausgefuehrt und muss nur bei zu grosser Differenz wiederholt werden */
	}
	
	return 0;
}

void get_force(double *x, double *v, double *a){
	double r1, r2, r3;
	double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;
    
    
    //Bulge first:
    
    //Point mass
    r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2) + b1*b1);
		
    a1x = -G*M1/(r1*r1*r1)**x;
    a1y = -G*M1/(r1*r1*r1)**(x+1);
    a1z = -G*M1/(r1*r1*r1)**(x+2);

    
    //alternatively
    //Gamma model for bulge, now only Hernquist bulge
    //NOT YET FULLY IMPLEMENTED, GAMMA MAY VARY
    r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));
    
    a1x += -G*M1_GAMMA/((r1+b1_GAMMA)*(r1+b1_GAMMA))**x/r1;
    a1y += -G*M1_GAMMA/((r1+b1_GAMMA)*(r1+b1_GAMMA))**(x+1)/r1;
    a1z += -G*M1_GAMMA/((r1+b1_GAMMA)*(r1+b1_GAMMA))**(x+2)/r1;

    //Disk next:
    
    //Miyamato disk
    r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2 + sqrt(*(x+2) * *(x+2) + b2*b2),2));
		
    a2x = -G*M2/(r2*r2*r2) * *x;
    a2y = -G*M2/(r2*r2*r2) * *(x+1);
    a2z = -G*M2/(r2*r2*r2) * (a2 + sqrt(*(x+2) * *(x+2) + b2*b2))/sqrt(*(x+2) * *(x+2) + b2*b2) * *(x+2);
        
    //Halo last:
    
    //Log Halo from Koposov et al. (2010)
    r3 = *x * *x + *(x+1) * *(x+1) + *(x+2)/q_halo * *(x+2)/q_halo + RCIRC*RCIRC; //R^2!
        
    a3x = -VCIRC *  *x/r3;
    a3y = -VCIRC *  *(x+1)/r3;
    a3z = -VCIRC *  *(x+2)/(q_halo*q_halo*r3);
        
    //alternatively
    //NFW halo
    r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/q_halo * *(x+2)/q_halo);
    
    a3x += -G*MNFW/r3 * (log(1.0 + r3/RNFW)/r3 - 1.0/(RNFW+r3)) * *x/r3;
    a3y += -G*MNFW/r3 * (log(1.0 + r3/RNFW)/r3 - 1.0/(RNFW+r3)) * *(x+1)/r3;
    a3z += -G*MNFW/r3 * (log(1.0 + r3/RNFW)/r3 - 1.0/(RNFW+r3)) * *(x+2)/(q_halo*q_halo*r3);
    
    //total tidal force
	*(a+0) = a1x + a2x + a3x;
	*(a+1) = a1y + a2y + a3y;
	*(a+2) = a1z + a2z + a3z;
    
}

void do_step(double dt, double *x, double *v) {
	double hh, *a0,*a1,*a2,*a3,xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	int k;
	
	a0 = malloc(3*sizeof(double));
	a1 = malloc(3*sizeof(double));
	a2 = malloc(3*sizeof(double));
	a3 = malloc(3*sizeof(double));
	
	hh = dt*0.5;
	
	get_force(x,v,a0);                    /* Zwischenschritte */
	for (k=0;k<3;k++) {                /* erster Halbschritt */
		xt1[k] = *(x+k)+hh**(v+k);
		vt1[k] = *(v+k)+hh**(a0+k);
	}  
	
	get_force(&xt1[0], &vt1[0], a1);
	for (k=0;k<3;k++) {                /* zweiter Halbschritt */
		xt2[k] = *(x+k)+hh*vt1[k];
		vt2[k] = *(v+k)+hh**(a1+k);
	}
	
	get_force(&xt2[0], &vt2[0], a2);
	for (k=0;k<3;k++) {                /* dritter Schritt mit Ergebnissen des zweiten */
		xt3[k] = *(x+k)+dt*vt2[k];
		vt3[k] = *(v+k)+dt**(a2+k);
	}
	
	get_force(&xt3[0], &vt3[0], a3);
	for (k=0;k<3;k++) {                /* Runge-Kutta-Formel */
		*(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
		*(v+k) += dt/6.0*(*(a0+k)+2.0*(*(a1+k)+*(a2+k))+*(a3+k));
	}
	free(a0); free(a1); free(a2); free(a3);   /* gibt Speicher wieder frei */
}

int radii(int Nt, double **start, int statust, double *LRt, double mxt, double myt, double mzt, double *Mt){
	int l;
	double **radiiarray;	
	radiiarray = (double **)calloc(Nt,sizeof(double *));
	for (l=0;l<Nt;l++){
		radiiarray[l] = (double *)calloc(2,sizeof(double));
		if (radiiarray[l] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}	
	
	
	//	float radiiarray[Nt][2];
	double rtemp, Mgestemp = 0.0, Mtemp = 0.0;
	LRt[0] = 0.0; LRt[1] = 0.0; LRt[2] = 0.0; LRt[3] = 0.0; LRt[4] = 0.0; LRt[5] = 0.0;
	
	for (l=0;l<Nt;l++){
		if ((start[l][0]) && (start[l][13]<statust) && (start[l][13]>= 0.0)) {
			rtemp = sqrt(pow(start[l][1]-mxt,2) + pow(start[l][2]-myt,2) + pow(start[l][3]-mzt,2));
			radiiarray[l][0] = rtemp;
			Mgestemp += start[l][10];
			radiiarray[l][1] = start[l][10];
		} else {
			radiiarray[l][0] = 1000000;
			radiiarray[l][1] = 0;
		}
	}
	
	*Mt = Mgestemp;
	
	//qsort(radiiarray, sizeof(radiiarray)/sizeof(radiiarray[0]), sizeof(radiiarray[0]), (void *)cmpmy);
	shellsort_reverse(radiiarray, Nt, 2);
	
	for (l=0;l<Nt;l++) {
		Mtemp += radiiarray[l][1];
		if ((LRt[0] == 0) && (Mtemp >= 0.02*Mgestemp))
			LRt[0] = radiiarray[l][0];
		else if ((LRt[1] == 0) && (Mtemp >= 0.05*Mgestemp))
			LRt[1] = radiiarray[l][0];
		else if ((LRt[2] == 0) && (Mtemp >= 0.1*Mgestemp))
			LRt[2] = radiiarray[l][0];
		else if ((LRt[3] == 0) && (Mtemp >= 0.2*Mgestemp))
			LRt[3] = radiiarray[l][0];
		else if ((LRt[4] == 0) && (Mtemp >= 0.5*Mgestemp))
			LRt[4] = radiiarray[l][0];
		else if ((LRt[5] == 0) && (Mtemp >= 0.9*Mgestemp))
			LRt[5] = radiiarray[l][0];
	}
	
	for (l=0;l<Nt;l++) free (radiiarray[l]);
	free(radiiarray);
	
	return 1;
}

int Kradii(int Nt, double **start, int statust, double *LRt, double mxt, double myt, double mzt, double *Mt, double *racc){
	int l;
	double **radiiarray;	
	radiiarray = (double **)calloc(Nt,sizeof(double *));
	for (l=0;l<Nt;l++){
		radiiarray[l] = (double *)calloc(2,sizeof(double));
		if (radiiarray[l] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}	
	//	float radiiarray[Nt][2];
	double rtemp, Mgestemp = 0.0, Mtemp = 0.0;
	LRt[0] = 0.0; LRt[1] = 0.0; LRt[2] = 0.0; LRt[3] = 0.0; LRt[4] = 0.0; LRt[5] = 0.0; LRt[5] = 0.0; LRt[6] = 0.0; LRt[7] = 0.0; LRt[8] = 0.0; LRt[9] = 0.0; LRt[10] = 0.0;
	
	for (l=0;l<Nt;l++){
		if ((start[l][0]) && (start[l][13]>=0.0)) {
			rtemp = sqrt(pow(start[l][1]-mxt,2) + pow(start[l][2]-myt,2) + pow(start[l][3]-mzt,2));
			radiiarray[l][0] = rtemp;
			radiiarray[l][1] = start[l][10];
		} else {
			radiiarray[l][0] = 1000000.0;
			radiiarray[l][1] = 0;
		}
	}
	
	Mgestemp = *Mt;
	
	//	qsort(radiiarray, sizeof(radiiarray)/sizeof(radiiarray[0]), sizeof(radiiarray[0]), (void *)cmpmy);
	shellsort_reverse(radiiarray, Nt, 2);
	
	double acc; //local acceleration
	*racc = 0.0; //radius at which a < a0
	double a0 = 3.6;//pc/Myr^2   //1.2e-10; //a0 in m/s^2
	
	for (l=0;l<Nt;l++) {
		Mtemp += radiiarray[l][1];
		if ((*racc==0.0) && (l>3)) {
			acc = 1.0*G*Mtemp/pow(radiiarray[l][0],2);
			if (acc < a0) *racc = radiiarray[l][0]; 
		}
		if ((LRt[0] == 0) && (Mtemp >= 0.02*Mgestemp))
			LRt[0] = radiiarray[l][0];
		else if ((LRt[1] == 0) && (Mtemp >= 0.05*Mgestemp))
			LRt[1] = radiiarray[l][0];
		else if ((LRt[2] == 0) && (Mtemp >= 0.1*Mgestemp))
			LRt[2] = radiiarray[l][0];
		else if ((LRt[3] == 0) && (Mtemp >= 0.2*Mgestemp))
			LRt[3] = radiiarray[l][0];
		else if ((LRt[4] == 0) && (Mtemp >= 0.3*Mgestemp))
			LRt[4] = radiiarray[l][0];
		else if ((LRt[5] == 0) && (Mtemp >= 0.4*Mgestemp))
			LRt[5] = radiiarray[l][0];
		else if ((LRt[6] == 0) && (Mtemp >= 0.5*Mgestemp))
			LRt[6] = radiiarray[l][0];
		else if ((LRt[7] == 0) && (Mtemp >= 0.6*Mgestemp))
			LRt[7] = radiiarray[l][0];
		else if ((LRt[8] == 0) && (Mtemp >= 0.7*Mgestemp))
			LRt[8] = radiiarray[l][0];
		else if ((LRt[9] == 0) && (Mtemp >= 0.8*Mgestemp))
			LRt[9] = radiiarray[l][0];
		else if ((LRt[5] == 0) && (Mtemp >= 0.9*Mgestemp))
			LRt[5] = radiiarray[l][0];
	}
	
	for (l=0;l<Nt;l++) free (radiiarray[l]);
	free(radiiarray);
	
	return 1;
}

int surfacedensityprofile(int Nt, double **start, double **profilet, double mxt, double myt, double mzt, double mvxt, double mvyt, double mvzt, int viewt, int numberradialbinst, double *VGt, double **profile_inside_rtidet, double **profile_bound_starst){
	int l;
	double **radiiarray;	
	radiiarray = (double **)calloc(Nt,sizeof(double *));
	for (l=0;l<Nt;l++){
		radiiarray[l] = (double *)calloc(5,sizeof(double));
		if (radiiarray[l] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}	
	double rtemp;
	
	for (l=0;l<Nt;l++){
		if ((start[l][0]) && (start[l][13]>=0.0)) {
			if (viewt == 1) {
				rtemp = sqrt(pow(start[l][1]-mxt,2) + pow(start[l][2]-myt,2));//xy-Ebene
			} else if (viewt == 2) {
				rtemp = sqrt(pow(start[l][1]-mxt,2) + pow(start[l][3]-mzt,2));//xz-Ebene
			}
			radiiarray[l][0] = rtemp;
			radiiarray[l][1] = start[l][10];//Masse
			//			radiiarray[l][2] = pow(start[l][4]+VGt[0],2)+pow(start[l][5]+VGt[1],2)+pow(start[l][6]+VGt[2],2);//v^2
			radiiarray[l][2] = pow(start[l][6]+VGt[2],2);//v
			radiiarray[l][3] = start[l][13];//status
			radiiarray[l][4] = start[l][8];//luminosity
		} else {
			radiiarray[l][0] = 1000000.0;
			radiiarray[l][1] = 0;
			radiiarray[l][2] = 0;
			radiiarray[l][3] = -1;
			radiiarray[l][4] = 0;//luminosity
		}
	}
	
	//qsort(radiiarray, sizeof(radiiarray)/sizeof(radiiarray[0]), sizeof(radiiarray[0]), (void *)cmpmy);
	shellsort_reverse(radiiarray, Nt, 5);
	
	//Prozedur fuer profile, profile_inside_rtide und profile_bound_stars machen
	//profile (all stars)
	int ll = 0;
	l = 0;
	while ((ll<numberradialbinst) && (l<Nt)) {
		if ((radiiarray[l][3]>=0) && (radiiarray[l][0]<profilet[ll][0])) {
			profilet[ll][2]+=radiiarray[l][1];
			profilet[ll][3]+=sqrt(radiiarray[l][2]);
			profilet[ll][4]+=radiiarray[l][2];
			profilet[ll][5]++;
			profilet[ll][7]+=radiiarray[l][4];
			l++;
		} else {
			ll++;
		}
	}
	
	ll = 0;
	l = 0;
	while ((ll<numberradialbinst) && (l<Nt)) {
		if ((radiiarray[l][3]>=0) && (radiiarray[l][0]<profilet[ll][0])) {
			profilet[ll][6] += pow(pow(sqrt(radiiarray[l][2])-profilet[ll][3]/profilet[ll][5],2) - (profilet[ll][4]/profilet[ll][5]-pow(profilet[ll][3]/profilet[ll][5],2)),2);
			l++;
		} else {
			ll++;
		}
	}
	
	
	//profile_inside_rtide (all stars inside rtide)
	ll = 0;
	l = 0;
	while ((ll<numberradialbinst) && (l<Nt)) {
		if ((radiiarray[l][3]>=0) && (radiiarray[l][3]<2) && (radiiarray[l][0]<profile_inside_rtidet[ll][0])) {
			profile_inside_rtidet[ll][2]+=radiiarray[l][1];
			profile_inside_rtidet[ll][3]+=sqrt(radiiarray[l][2]);
			profile_inside_rtidet[ll][4]+=radiiarray[l][2];
			profile_inside_rtidet[ll][5]++;
			profile_inside_rtidet[ll][7]+=radiiarray[l][4];
			l++;
		} else if ((radiiarray[l][3]>1) || (radiiarray[l][3]<0)) {
			l++;
		} else {
			ll++;
		}
	}
	
	ll = 0;
	l = 0;
	while ((ll<numberradialbinst) && (l<Nt)) {
		if ((radiiarray[l][3]>=0) && (radiiarray[l][3]<2) && (radiiarray[l][0]<profile_inside_rtidet[ll][0])) {
			profile_inside_rtidet[ll][6] += pow(pow(sqrt(radiiarray[l][2])-profile_inside_rtidet[ll][3]/profile_inside_rtidet[ll][5],2) - (profile_inside_rtidet[ll][4]/profile_inside_rtidet[ll][5]-pow(profile_inside_rtidet[ll][3]/profile_inside_rtidet[ll][5],2)),2);
			l++;
		} else if ((radiiarray[l][3]>1) || (radiiarray[l][3]<0)) {
			l++;
		} else {
			ll++;
		}
	}
	
	
	//profile_bound_stars (all bound stars)
	ll = 0;
	l = 0;
	while ((ll<numberradialbinst) && (l<Nt)) {
		if ((radiiarray[l][3] == 0.0) && (radiiarray[l][0]<profile_bound_starst[ll][0])) {
			profile_bound_starst[ll][2]+=radiiarray[l][1];
			profile_bound_starst[ll][3]+=sqrt(radiiarray[l][2]);
			profile_bound_starst[ll][4]+=radiiarray[l][2];
			profile_bound_starst[ll][5]++;
			profile_bound_starst[ll][7]+=radiiarray[l][4];
			l++;
		} else if ((radiiarray[l][3]>=1.0) || (radiiarray[l][3]<=-1.0)) {
			l++;
		} else {
			ll++;
		}
	}
	
	ll = 0;
	l = 0;
	while ((ll<numberradialbinst) && (l<Nt)) {
		if ((radiiarray[l][3]<1.0) && (radiiarray[l][3]>-1.0) && (radiiarray[l][0]<profile_bound_starst[ll][0])) {
			profile_bound_starst[ll][6] += pow(pow(sqrt(radiiarray[l][2])-profile_bound_starst[ll][3]/profile_bound_starst[ll][5],2) - (profile_bound_starst[ll][4]/profile_bound_starst[ll][5]-pow(profile_bound_starst[ll][3]/profile_bound_starst[ll][5],2)),2);
			l++;
		} else if ((radiiarray[l][3]>=1.0) || (radiiarray[l][3]<-1.0)) {
			l++;
		} else {
			ll++;
		}
	}
	
	for (l=0;l<Nt;l++) free (radiiarray[l]);
	free(radiiarray);
	
	return 1;
	
}

int get_velocitydispersion(int Nt,double **start,int statust, double *sigmat, double *mvxt, double *mvyt, double *mvzt){
	int l, n = 0;
	double sigmax = 0.0, sigmay  = 0.0, sigmaz = 0.0;
	*sigmat = 0.0;
	*mvxt = 0.0; *mvyt = 0.0; *mvzt = 0.0;
	double mv = 0.0, mv2 = 0.0, mtemp = 0.0;
	double vtemp = 0.0, sigmatemp = 0.0;
	double stdsigma = 5000.0;//optional velocity cut, disabled by default
	int nsigma = 100.0;//optional velocity cut, disabled by default 
	
	for (l=0;l<Nt;l++) {
		vtemp = sqrt(start[l][4]*start[l][4]+start[l][5]*start[l][5]+start[l][6]*start[l][6]);
		if (n>10) sigmatemp = sqrt((mv2-1.0/(1.0*n)*mv*mv)/(1.0*n-1.0));
		else sigmatemp = stdsigma; 
		
		if ((start[l][0]) && (start[l][13]<statust) && (start[l][13] >= 0) && (vtemp < nsigma*sigmatemp)) {
			*mvxt += start[l][10]*start[l][4];
			*mvyt += start[l][10]*start[l][5];
			*mvzt += start[l][10]*start[l][6];
			mtemp += start[l][10];
			mv2 += vtemp*vtemp;
			mv += vtemp;
			n++;
		}
	}
	
	if (n) {
		*mvxt /= mtemp;
		*mvyt /= mtemp;
		*mvzt /= mtemp;
	}
	
	sigmatemp = sqrt((mv2-1.0/(1.0*n)*mv*mv)/(1.0*n-1.0));
	//printf("sigmatemp = %f\n",sigmatemp);
	n = 0;
	
	for (l=0;l<Nt;l++) {
		vtemp = sqrt(pow(start[l][4]-*mvxt,2)+pow(start[l][5]-*mvyt,2)+pow(start[l][6]-*mvzt,2));
		//printf("%f\t%f\n",vtemp,sigmatemp);
		if ((start[l][0]) && (start[l][13]<statust) && (start[l][13] >= 0) && (vtemp < nsigma*sigmatemp)) {
			sigmax += pow(start[l][4]-*mvxt,2);
			sigmay += pow(start[l][5]-*mvyt,2);
			sigmay += pow(start[l][6]-*mvzt,2);
			n++;
		}
	}
	if (n) {
		sigmax /= n;
		sigmay /= n;
		sigmaz /= n;
	}
	*sigmat = (sigmax + sigmay + sigmaz)/3.0;
	
	return 1;
}

int get_binaries(int Nt, double **start, int **neighborlistt, int neighborst, int *Ntott, int *multiplicityarrayt, int binoutputt, FILE *bint, double timet, int outputtypet){
	int l, ll, namell;
	double Ebintemp, atemp, etemp, rtemp,rx,ry,rz,vx,vy,vz;
	
	int name1, name2, multiplicity;
	for (l=0;l<neighborst;l++) multiplicityarrayt[l] = 0;
	
	double **binarray;	
	binarray = (double **)calloc(Nt,sizeof(double *));
	for (l=0;l<Nt;l++){
		binarray[l] = (double *)calloc(3,sizeof(double));
		if (binarray[l] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
		binarray[l][0] = 0.0;
	}	
	
	
	int nrbin = 0;
	
	for (l=0;l<Nt;l++) {
		if (start[l][0]) {
			for (ll=0;ll<1;ll++) {
				//for (ll=0;ll<neighborst;ll++) {
				namell = neighborlistt[l][ll];
				rx = start[l][1]-start[namell][1];
				ry = start[l][2]-start[namell][2];
				rz = start[l][3]-start[namell][3];
				rtemp = sqrt(rx*rx+ry*ry+rz*rz);
				
				vx = start[l][4]-start[namell][4];
				vy = start[l][5]-start[namell][5];
				vz = start[l][6]-start[namell][6];
				
				Ebintemp = 0.5*start[l][10]*start[namell][10]/(start[l][10]+start[namell][10])*(vx*vx+vy*vy+vz*vz) - G*start[l][10]*start[namell][10]/rtemp;
				
				atemp = -0.5*G*start[l][10]*start[namell][10]/Ebintemp;
				
				etemp = sqrt(pow(1.0-rtemp/atemp,2) - pow(rx*vx+ry*vy+rz*vz,2)/(G*(start[l][10]+start[namell][10])*atemp));
				
				if (etemp<1.0) {
					if (binoutputt) printf("Binary:\t1 = %i\t2 = %i (%ith)\tr12 = %g\tEges = %g\ta = %g\tecc = %g\t1: %g\t2: %g\n", l, namell, ll, rtemp, Ebintemp, atemp, etemp, start[l][13],start[namell][13]);
					if (binoutputt) fprintf(bint,"%.0f\t%i\t%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",timet , l, namell, ll, rtemp, Ebintemp, atemp, etemp, start[l][10],start[namell][10],sqrt(pow(start[l][1],2)+pow(start[l][2],2)+pow(start[l][3],2)));
					binarray[nrbin][0] = Ebintemp;
					binarray[nrbin][1] = l;
					binarray[nrbin][2] = namell;
					nrbin++;
				} 
			}
		}
	}
	
	shellsort_reverse(binarray, Nt, 3);
	
	ll = 1;
	for (l=0;l<nrbin;l++) {
		if (binarray[l][1]-binarray[l+1][2]) {
			multiplicity = 2;
			name1 = (int) binarray[l][1];
			name2 = (int) binarray[l][2];
			
			while (start[name1][13] < 0) {
				name1 = (int) Nt-start[name1][13];
				multiplicity++;
			}
			while (start[name2][13] < 0) {
				name2 = (int) Nt-start[name2][13];
				multiplicity++;
			}
			
			//Hoehere Multiplizitaeten auf Bindung mit dem Centre of Mass des stablieren Binaries pruefen
			if (multiplicity>2) {
				rx = start[name1][1]-start[name2][1];
				ry = start[name1][2]-start[name2][2];
				rz = start[name1][3]-start[name2][3];
				rtemp = sqrt(rx*rx+ry*ry+rz*rz);
				
				vx = start[name1][4]-start[name2][4];
				vy = start[name1][5]-start[name2][5];
				vz = start[name1][6]-start[name2][6];
				
				Ebintemp = 0.5*start[name1][10]*start[name2][10]/(start[name1][10]+start[name2][10])*(vx*vx+vy*vy+vz*vz) - G*start[name1][10]*start[name2][10]/rtemp;
				
				atemp = -0.5*G*start[name1][10]*start[name2][10]/Ebintemp;
				
				etemp = sqrt(pow(1.0 - rtemp/atemp,2)+ pow(rx*vx+ry*vy+rz*vz,2)/(G*(start[name1][10]+start[name2][10])*atemp));
				
			} else {
				etemp = 0.0;
			}	
			
			if (etemp<1.0) {
				rx = start[name1][1]-start[name2][1];
				ry = start[name1][2]-start[name2][2];
				rz = start[name1][3]-start[name2][3];
				rtemp = sqrt(rx*rx+ry*ry+rz*rz);
				
				vx = start[name1][4]-start[name2][4];
				vy = start[name1][5]-start[name2][5];
				vz = start[name1][6]-start[name2][6];
				
				Ebintemp = 0.5*start[name1][10]*start[name2][10]/(start[name1][10]+start[name2][10])*(vx*vx+vy*vy+vz*vz) - G*start[name1][10]*start[name2][10]/rtemp;
				
				atemp = -0.5*G*start[name1][10]*start[name2][10]/Ebintemp;
				
				if (binoutputt) printf("%g\t%i\t%i\t%i\t-->\t%i (distance from centre: %fpc, semi-major axis: %fpc)\n", binarray[l][0], name1, name2, multiplicity,Nt+ll,sqrt(pow(start[name1][1],2)+pow(start[name1][2],2)+pow(start[name1][3],2)),atemp);
				start[name1][13] = -1*ll;
				start[name2][13] = -1*ll;
				
				
				//CoM erstellen
				start[Nt+ll][0] = -1.0;
				start[Nt+ll][1] = (start[name1][1]*start[name1][10]+start[name2][1]*start[name2][10])/(start[name1][10]+start[name2][10]);
				start[Nt+ll][2] = (start[name1][2]*start[name1][10]+start[name2][2]*start[name2][10])/(start[name1][10]+start[name2][10]);
				start[Nt+ll][3] = (start[name1][3]*start[name1][10]+start[name2][3]*start[name2][10])/(start[name1][10]+start[name2][10]);
				start[Nt+ll][4] = (start[name1][4]*start[name1][10]+start[name2][4]*start[name2][10])/(start[name1][10]+start[name2][10]);
				start[Nt+ll][5] = (start[name1][5]*start[name1][10]+start[name2][5]*start[name2][10])/(start[name1][10]+start[name2][10]);
				start[Nt+ll][6] = (start[name1][6]*start[name1][10]+start[name2][6]*start[name2][10])/(start[name1][10]+start[name2][10]);
				start[Nt+ll][7] = -1.0;
				start[Nt+ll][8] = start[name1][8]+start[name2][8];
				start[Nt+ll][9] = -1.0;
				start[Nt+ll][10] = start[name1][10]+start[name2][10];
				start[Nt+ll][11] = 0.0;
				start[Nt+ll][12] = 0.5*start[Nt+ll][10]*(pow(start[Nt+ll][4],2)+pow(start[Nt+ll][5],2)+pow(start[Nt+ll][6],2));
				start[Nt+ll][13] = 0.0;
				if ((start[name1][14]) && (start[name1][14] < 2.0) && (start[name2][14]) && (start[name2][14] < 2.0)) {
					start[Nt+ll][14] = 1.0;
					start[name1][14] = 2.0;
					start[name2][14] = 2.0;
					if (outputtypet) printf("BINARY ESCAPE:\t%i\t%i\n", name1, name2);
				} else {
					start[Nt+ll][14] = 0.0;
				}
				
				*Ntott = *Ntott + 1;
				ll++;
				multiplicityarrayt[multiplicity]++;
			}
		}
	}
	
	for (l=Nt;l<*Ntott;l++) {
		if ((start[l][0]) && (start[l][13] >= 0)) {
			for (ll=0;ll<*Ntott;ll++) {
				if ((l-ll) && (start[ll][0]) && (start[ll][13] >= 0.0)) {
					rtemp = sqrt(pow(start[l][1]-start[ll][1],2)+pow(start[l][2]-start[ll][2],2)+pow(start[l][3]-start[ll][3],2));
					if (rtemp) start[l][11] += -G*start[l][10]*start[ll][10]/rtemp;
				}
			}
		}
	}
	
	for (l=0;l<Nt;l++) free (binarray[l]);
	free(binarray);
	
	return 1;
}

int get_energy(int Nt, double **start, int statust, double *Etemp, double mvxt, double mvyt, double mvzt){
	int l;
	Etemp[0] = 0.0; //Ebin
	Etemp[1] = 0.0; //Ekin
	Etemp[2] = 0.0; //Eges
	
	for (l=0;l<Nt;l++) {
		if ((start[l][0]) && (start[l][13] < statust) && (start[l][13] >= 0)) {
			Etemp[0] += start[l][11];
			Etemp[1] += 0.5*start[l][10]*(pow(start[l][4]-mvxt,2)+pow(start[l][5]-mvyt,2)+pow(start[l][6]-mvzt,2));
		}
	}
	Etemp[0] /= 2.0;
	Etemp[2] = Etemp[0] + Etemp[1];
	
	return 1;
}

int get_status(int Nt, double **start, double rtidet, double mxt, double myt, double mzt, int *nt, int *nbint, double mvxt, double mvyt, double mvzt, double rh, double omega, int quick_and_dirty_analysis){
	int l, ltemp;
	nt[0] = 0; nt[1] = 0; nt[2] = 0; nt[3] = 0; nbint[0] = 0; nbint[1] = 0; nbint[2] = 0; nbint[3] = 0;
	double rtemp, etemp,mtemp,ERT,rtemp2;
	mtemp = 0.0;
	
	//Create radial mass profile for rough computation of binding energy
	int gridpoints = 100;
	double massprofile[gridpoints][2];
	for (l=0;l<gridpoints;l++) {
		massprofile[l][0] = rtidet*l/gridpoints;
		massprofile[l][1] = 0.0;
	}
	
	for (l=0;l<Nt;l++) {
		if ((start[l][0]) && (start[l][13] >= 0.0)){
			rtemp = sqrt(pow(start[l][1]-mxt,2)+pow(start[l][2]-myt,2)+pow(start[l][3]-mzt,2));
			if (rtemp < rtidet) {
				start[l][13] = 0.0;
				nt[0]++ ;
				mtemp += start[l][10];
				ltemp = floor(1.0*rtemp/rtidet*gridpoints);
				massprofile[ltemp][1] += start[l][10];
			} else if (rtemp < 2.0*rtidet) {
				start[l][13] = 2.0;
				nt[2]++;
			} else {
				start[l][13] = 3.0;
				nt[3]++;
			}
		} else if ((start[l][0]) && (start[l][13] < 0.0)) {
			rtemp = sqrt(pow(start[l][1]-mxt,2)+pow(start[l][2]-myt,2)+pow(start[l][3]-mzt,2));
			if (rtemp < rtidet) {
				nbint[0]++ ;
				mtemp += start[l][10];
			} else if (rtemp < 2*rtidet) {
				nbint[2]++;
			} else {
				nbint[3]++;
			}
		}	
	}
	for (l=0;l<4;l++) nbint[l] = nbint[l]/2.0;
	
	double mtemp2 = 0.0;
	for (l=0;l<gridpoints;l++) {
		mtemp2 += massprofile[l][1];
		massprofile[l][1] = mtemp2;
	}
	
	
	//Potential escapers suchen & Bindunsenergien fuer quick_and_dirty_analysis berechnen
	ERT = -1.5*G*mtemp/rtidet;
	
	for (l=0;l<Nt;l++) {
		if ((start[l][0]) && (start[l][13] == 0.0)){
			if (quick_and_dirty_analysis) {
				rtemp2 = sqrt(pow(start[l][1]-mxt,2)+pow(start[l][2]-myt,2)+pow(start[l][3]-mzt,2));
				ltemp = floor(1.0*rtemp2/rtidet*gridpoints);
				start[l][11] = -G*start[l][10]*mtemp/rtemp2;
				//start[l][11] = - G*mtemp*start[l][10]*(pow(1.305*rtemp2/rh,3)/pow(1.0+pow(1.305*rtemp2/rh,2),1.5)); //assumed Plummer profile for binding energy!
			}
			
			etemp = 0.5*start[l][10]*(pow(start[l][4]-mvxt,2)+pow(start[l][5]-mvyt,2)+pow(start[l][6]-mvzt,2)) + 0.5*start[l][10]*omega*omega*(pow(start[l][3]-mzt,2) - 2.0*pow(start[l][1]-mxt,2)) + start[l][11]; //assumed near-field approximation for binding energy!
			
			
			if (etemp > start[l][10]*ERT) {
				start[l][13] = 1.0;
				nt[0]-- ;
				nt[1]++ ;
			}
		} 
		
	}
	
	return 1;
}

int get_rtide(double *rtidet, int Nt, double **start, double omegat, double mxt, double myt, double mzt, double dphit, double *RGt){
	int l;
	double rtidetemp, Mtemp = 0.0;
	if (*rtidet > 1.0) { 
		for (l=0;l<Nt;l++) {
			if (start[l][0]) Mtemp += start[l][10];
		}
	} else {
		Mtemp = Nt;
	}
	do{
		rtidetemp = pow(G*Mtemp/sqrt(pow((dphit+omegat*omegat),2)),1.0/3.0);
		Mtemp = 0.0;
		for (l=0;l<Nt;l++) {
			if ((start[l][0]) && (sqrt(pow(start[l][1]-mxt,2)+pow(start[l][2]-myt,2)+pow(start[l][3]-mzt,2)) < rtidetemp)) {
				Mtemp += start[l][10];
			}
		}
		*rtidet = pow(G*Mtemp/sqrt(pow((dphit+omegat*omegat),2)),1.0/3.0);
		printf("+");
	} while (sqrt(pow(*rtidet-rtidetemp,2))/ *rtidet > 0.01);
	printf("\n");
	
	return 1;
}

int get_omega(double *RGt, double *VGt, double *omegat, double *Lt) {
	if (RGt[0]) {
		*Lt = sqrt(pow(RGt[1]*VGt[2]-RGt[2]*VGt[1],2)+pow(RGt[2]*VGt[0]-RGt[0]*VGt[2],2)+pow(RGt[0]*VGt[1]-RGt[1]*VGt[0],2));
		*omegat = *Lt/(RGt[0]*RGt[0]+RGt[1]*RGt[1]+RGt[2]*RGt[2]);
	} 
	
	return 1;
}

int get_densitycenter(int Nt, double **start, int **neighborlistt, double *mxt, double *myt, double *mzt, int densitylevelt){ 
	int l, ll;
	double mtemp, rhotemp, rhogestemp = 0.0;
	double softening = 0.001;
	
	//float rhotemplist[Nt];
	
	*mxt = 0; *myt = 0; *mzt = 0;
	
	for (l=0;l<Nt;l++) {
		if ((start[l][0]) && (start[l][13]<3)) {
			mtemp = 0.0;
			for (ll=0;ll<densitylevelt;ll++) {
				mtemp += start[neighborlistt[l][ll]][10];//Mass
				//mtemp += start[neighborlistt[l][ll]][8];//Luminosity
			}
			rhotemp = mtemp/pow(pow(start[neighborlistt[l][densitylevelt]][1]-start[l][1],2)+pow(start[neighborlistt[l][densitylevelt]][2]-start[l][2],2)+pow(start[neighborlistt[l][densitylevelt]][3]-start[l][3],2)+pow(softening,2),1.5);
		} else {
			rhotemp = 0;
		}
		
		//rhotemplist[l]=rhotemp;
		
		*mxt += start[l][1]*rhotemp;
		*myt += start[l][2]*rhotemp;
		*mzt += start[l][3]*rhotemp;
		rhogestemp += rhotemp; 
	}
	//qsort(rhotemplist, Nt, sizeof(float), (void *)cmpmy2);
	//printf("\n%g\t%g\t%g\t%g\n", rhotemplist[0], rhotemplist[1], rhotemplist[2], rhotemplist[3]);
	
	*mxt /= rhogestemp;
	*myt /= rhogestemp;
	*mzt /= rhogestemp;
	
	return 1;
}

int cmpmy(float *x1, float *x2) { //aufsteigend
	if(*x1<*x2) return -1;
	return 1;
}

int cmpmy2(float *x1, float *x2) { //absteigend
	if(*x1>*x2) return -1;
	return 1;
}

void shellsort(double **array, int N, int k) {//largest up
	int i,j,l,n;
	N = N-1;
	double swap[k];
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			for (l=0; l<k; l++) swap[l] = array[i][l];
			for (j = i; ((j >= n) && (array[j-n][0] < swap[0])); j -= n) {
				for (l=0; l<k; l++) array[j][l] = array[j-n][l];
			}
			for (l=0; l<k; l++) array[j][l] = swap[l];
		}
	}
}

void shellsort_reverse(double **array, int N, int k) {//smallest up
	int i,j,l,n;
	N = N-1;
	double swap[k];
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			for (l=0; l<k; l++) swap[l] = array[i][l];
			for (j = i; ((j >= n) && (array[j-n][0] > swap[0])); j -= n) {
				for (l=0; l<k; l++) array[j][l] = array[j-n][l];
			}
			for (l=0; l<k; l++) array[j][l] = swap[l];
		}
	}
}

void shellsort_1d(double *array, int N) {//largest up
	int i,j,n;
	N = N-1;
	double swap;
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			swap = array[i];
			for (j = i; ((j >= n) && (array[j-n] < swap)); j -= n) {
				array[j] = array[j-n];
			}
			array[j] = swap;
		}
	}
}

void shellsort_reverse_1d(double *array, int N) {//smallest up
	int i,j,n;
	N = N-1;
	double swap;
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			swap = array[i];
			for (j = i; ((j >= n) && (array[j-n] > swap)); j -= n) {
				array[j] = array[j-n];
			}
			array[j] = swap;
		}
	}
}

int bound_star_analysisfast(int Nt, double **start, int **neighborlistt, int neighborst, double mxt, double myt, double mzt) {
	int l, ll, i, j, namei, namej, namel;
	double rtemp;
	double ri, rj;
	
	float rtemplist[Nt][2];
	float rtemplist2[neighborst][2];
	
	for (l=0;l<Nt;l++) {
		if (start[l][0]) {
			rtemp = sqrt(pow(start[l][1] - mxt,2)+pow(start[l][2] - myt,2)+pow(start[l][3] - mzt,2));
		} else {
			rtemp = 1000000.0;
		}
		rtemplist[l][0] = rtemp;
		rtemplist[l][1] = l;
	}
	
	qsort(rtemplist, Nt, sizeof(rtemplist[0]), (void *)cmpmy);
	
	for (l=0;l<Nt;l++) {
		i=1; //oben
		j=1; //unten
		namel = (int) rtemplist[l][1];
		if (start[namel][0]) {
			for (ll=1;ll<=neighborst;ll++) {
				//ll ist die Position, in die ins neighborarray geschrieben wird.
				//dann einen weiteren counter, z.B. i, der nach oben und einen anderen, z.B. j, der nach unten zaehlt.
				//dann beide radii ausrechnen und vergleichen.
				//kleineren der beiden in das neighborarray schreiben und ll erhoehen.
				if (l+i<Nt) {
					namei = (int) rtemplist[l+i][1];
					ri = sqrt(pow(start[namel][1]-start[namei][1],2)+pow(start[namel][2]-start[namei][2],2)+pow(start[namel][3]-start[namei][3],2));
				} else {
					ri = 10000000.0;
				}
				
				if (l-j>=0) {
					namej = (int) rtemplist[l-j][1];
					rj = sqrt(pow(start[namel][1]-start[namej][1],2)+pow(start[namel][2]-start[namej][2],2)+pow(start[namel][3]-start[namej][3],2));
				} else {
					rj = 10000000.0;
				}
				
				if (ri<rj) {
					//neighborlistt[namel][ll-1] = namei;
					rtemplist2[ll-1][1] = namei;
					rtemplist2[ll-1][0] = ri;
					i++;
					if (i>j+6) j++;
				} else {
					//neighborlistt[namel][ll-1] = namej;
					rtemplist2[ll-1][1] = namej;
					rtemplist2[ll-1][0] = rj;
					j++;
					if (j>i+6) i++;
				}
			}
			qsort(rtemplist2, neighborst, sizeof(rtemplist[0]), (void *)cmpmy);
			for (ll=0;ll<neighborst;ll++){
				neighborlistt[namel][ll] = rtemplist2[ll][1];
			}
		}
	}
	return 1;
}

int bound_star_analysis(int Nt, double **start, int **neighborlistt, int neighborst){
	int l, ll;
	double Ebintemp, rtemp;
	
	float rtemplist[Nt][2];
	
	for (l=0;l<Nt;l++) {
		if (start[l][0]) {
			//printf("%i\n",l);
			rtemplist[l][0] = 1000000.0;
			rtemplist[l][1] = l;
			
			for (ll=0;ll<Nt;ll++) {
				if ((l-ll) && (start[ll][0])) {
					//radii berechnen
					rtemp = sqrt(pow(start[l][1]-start[ll][1],2)+pow(start[l][2]-start[ll][2],2)+pow(start[l][3]-start[ll][3],2));
					rtemplist[ll][0] = rtemp;
					rtemplist[ll][1] = ll;
				} else {
					rtemplist[ll][0] = 1000000.0;
					rtemplist[ll][1] = ll;
				}	
			}
			qsort(rtemplist, Nt, sizeof(rtemplist[0]), (void *)cmpmy);
			for (ll=0;ll<neighborst;ll++) {
				neighborlistt[l][ll] = rtemplist[ll][1];
			}
			
			//Bindungsenergien aufsummieren bis ca 0.01% Genauigkeit
			start[l][11] = 0.0;
			ll=0;	
			do{
				Ebintemp = -G*start[l][10]*start[(int) rtemplist[ll][1]][10]/rtemplist[ll][0];
				start[l][11]+=Ebintemp;
				ll++;
			} while ((ll<Nt) && (Ebintemp/start[l][11]>0.0000001));
		}
	}
	
	return 1;
}

void convert(double *xtemp, double *vtemp, double *dsuntemp, double *vrsuntemp, double *vrtemp, double *ltemp, double *btemp, double *lcosbtemp, double *RAtemp, double *DECtemp, double *mu_alphatemp, double *mu_alphacosdeltatemp, double *mu_deltatemp, double *mutemp, double *PAtemp, int coordtype, int vcoordtype, int radiococo, double vLSRtemp){
	
	double x,y,z;        //galactic coordinates [kpc]
	double xsun = -rgalsun/1000.0;//galactocentric distance of sun [kpc]
	double dsun;			//heliocentric radial distance [kpc]
	double dxy;             //heliocentric distance in xy-plane [kpc]
	double vx,vy,vz;        //cluster velocity [km/s]
	double vrx,vry,vrz;     //cluster radial velocity in 3d coordinates [km/s]
	double vtx,vty,vtz;     //cluster tangential velocity in 3d coordinates [km/s]
	double t,d,a;
	double T[3][3], A[3][3], B[3][3];
	double TI[3][3];
	double detT;
	double RArad, DECrad;
	double brad, lrad, lcosbrad;
	double mu_alphacosdelta, mu_delta, mu, PArad;
	double vr, vrsun;
	double RAENP = 0.0, DECENP = PI/2.0;  //equatorial coordinates of equatorial north pole
	double xENP, yENP, zENP, dxyENP; //cartesian vector pointing to the equatorial north pole  
	double bENP, lENP; //galactic coordinates of ENP
	double FAK;
	double xdelta, ydelta, zdelta;
	double nx, ny, nz, dvt;
	double vrLSR, vrGSR;
	
	
	//transformation matrix equ. -> gal. from Johnson & Soderblom (1987)
	t = PAGNP/360.0*2.0*PI;
	d = deltaGNP/360.0*2.0*PI;
	a = alphaGNP/360.0*2.0*PI;
	
	T[0][0] = -cos(t)*sin(d)*cos(a)-sin(t)*sin(a);
	T[0][1] = -cos(t)*sin(d)*sin(a)+sin(t)*cos(a);
	T[0][2] = cos(t)*cos(d);
	
	T[1][0] = -sin(t)*sin(d)*cos(a)+cos(t)*sin(a);
	T[1][1] = -sin(t)*sin(d)*sin(a)-cos(t)*cos(a);
	T[1][2] = sin(t)*cos(d);
	
	T[2][0] = cos(d)*cos(a);
	T[2][1] = cos(d)*sin(a);
	T[2][2] = sin(d);
	
	//invert matrix T in the most general way
	detT = T[0][0]*T[1][1]*T[2][2] + T[1][0]*T[2][1]*T[0][2] + T[2][0]*T[0][1]*T[1][2] - T[0][0]*T[1][2]*T[2][1] - T[1][0]*T[2][2]*T[0][1] - T[2][0]*T[0][2]*T[1][1];
	
	TI[0][0] = (T[1][1]*T[2][2]-T[1][2]*T[2][1])/detT;
	TI[1][0] = (T[1][2]*T[2][0]-T[1][0]*T[2][2])/detT;
	TI[2][0] = (T[1][0]*T[2][1]-T[1][1]*T[2][0])/detT;
	
	TI[0][1] = (T[0][2]*T[2][1]-T[0][1]*T[2][2])/detT;
	TI[1][1] = (T[0][0]*T[2][2]-T[2][0]*T[0][2])/detT;
	TI[2][1] = (T[0][1]*T[2][0]-T[0][0]*T[2][1])/detT;
	
	TI[0][2] = (T[0][1]*T[1][2]-T[0][2]*T[1][1])/detT;
	TI[1][2] = (T[0][2]*T[1][0]-T[0][0]*T[1][2])/detT;
	TI[2][2] = (T[0][0]*T[1][1]-T[0][1]*T[1][0])/detT;
	
	//convert to kpc
	x = xtemp[0]/1000.0;
	y = xtemp[1]/1000.0;
	z = xtemp[2]/1000.0;
	
	dsun = *dsuntemp/1000.0;
	
	vx = vtemp[0];
	vy = vtemp[1];
	vz = vtemp[2];
	
	vr = *vrtemp;
	vrsun = *vrsuntemp;
	
	//convert to radians
	DECrad = *DECtemp/360.0*2.0*PI;
	RArad = *RAtemp/360.0*2.0*PI;
	PArad = *PAtemp/360.0*2.0*PI;
	
	
	//get the galactic coordinates first
	if (coordtype == 1) {
		if (radiococo) printf("\nConverting equatorial to galactic coordinates using the transformation matrix:\n"); 
		if (radiococo) printf("%f\t%f\t%f\n",T[0][0],T[0][1],T[0][2]);
		if (radiococo) printf("%f\t%f\t%f\n",T[1][0],T[1][1],T[1][2]);
		if (radiococo) printf("%f\t%f\t%f\n",T[2][0],T[2][1],T[2][2]);
		
		brad = asin(T[2][0]*cos(DECrad)*cos(RArad) + T[2][1]*cos(DECrad)*sin(RArad) + T[2][2]*sin(DECrad));
		if (asin((T[1][0]*cos(DECrad)*cos(RArad) + T[1][1]*cos(DECrad)*sin(RArad) + T[1][2]*sin(DECrad))/cos(brad))>=0.0) {
			lrad = acos((T[0][0]*cos(DECrad)*cos(RArad) + T[0][1]*cos(DECrad)*sin(RArad) + T[0][2]*sin(DECrad))/cos(brad)); 
		} else {
			lrad = 2.0*PI-acos((T[0][0]*cos(DECrad)*cos(RArad) + T[0][1]*cos(DECrad)*sin(RArad) + T[0][2]*sin(DECrad))/cos(brad)); 			
		}
		lcosbrad = lrad*cos(brad);
	} else if (coordtype == 2) {
		brad = *btemp/360.0*2.0*PI;
		if (*ltemp) {
			lrad = *ltemp/360.0*2.0*PI;
			lcosbrad = lrad*cos(brad);
		} else if (*lcosbtemp) {
			lcosbrad = *lcosbtemp/360.0*2.0*PI;
			lrad = lcosbrad/cos(brad);
		}
	} else if (coordtype == 3) {
		if (y >= 0.0)
			lrad = acos((x-xsun)/sqrt(pow(x-xsun,2)+y*y));
		else
			lrad = 2.0*Pi-acos((x-xsun)/sqrt(pow(x-xsun,2)+y*y));
		brad =  atan(z/sqrt(pow(x-xsun,2)+y*y));
		lcosbrad = lrad*cos(brad);
	}
	
	
	//get 3d position of cluster [kpc] from galactic coordinates
	if (coordtype < 3) {
		z = sin(brad)*dsun;
		dxy = sqrt(dsun*dsun-z*z);
		x = cos(lrad)*dxy + xsun;
		y = sin(lrad)*dxy;
	} else if (coordtype == 3) {
		dsun = sqrt(pow(x-xsun,2)+y*y+z*z);		
	}
	
	
	//finally get the equatorial coordinates from galactic coordinates
	if (coordtype > 1) {
		if (radiococo) printf("\nConverting galactic to equatorial coordinates using the transformation matrix:\n"); 
		if (radiococo) printf("%f\t%f\t%f\n",TI[0][0],TI[0][1],TI[0][2]);
		if (radiococo) printf("%f\t%f\t%f\n",TI[1][0],TI[1][1],TI[1][2]);
		if (radiococo) printf("%f\t%f\t%f\n",TI[2][0],TI[2][1],TI[2][2]);
		
		if (radiococo) {
			//unit matrix B = T * TI
			B[0][0] = T[0][0]*TI[0][0] + T[0][1]*TI[1][0] + T[0][2]*TI[2][0];
			B[0][1] = T[0][0]*TI[0][1] + T[0][1]*TI[1][1] + T[0][2]*TI[2][1];
			B[0][2] = T[0][0]*TI[0][2] + T[0][1]*TI[1][2] + T[0][2]*TI[2][2];
			
			B[1][0] = T[1][0]*TI[0][0] + T[1][1]*TI[1][0] + T[1][2]*TI[2][0];
			B[1][1] = T[1][0]*TI[0][1] + T[1][1]*TI[1][1] + T[1][2]*TI[2][1];
			B[1][2] = T[1][0]*TI[0][2] + T[1][1]*TI[1][2] + T[1][2]*TI[2][2];
			
			B[2][0] = T[2][0]*TI[0][0] + T[2][1]*TI[1][0] + T[2][2]*TI[2][0];
			B[2][1] = T[2][0]*TI[0][1] + T[2][1]*TI[1][1] + T[2][2]*TI[2][1];
			B[2][2] = T[2][0]*TI[0][2] + T[2][1]*TI[1][2] + T[2][2]*TI[2][2];
			
			printf("\nCalculating T*T^{-1} = 1 for consistency check:\n");
			printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
			printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
			printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		}
		
		DECrad = asin(TI[2][0]*cos(brad)*cos(lrad)+TI[2][1]*cos(brad)*sin(lrad)+TI[2][2]*sin(brad));
		if (asin((TI[1][0]*cos(brad)*cos(lrad) + TI[1][1]*cos(brad)*sin(lrad) + TI[1][2]*sin(brad))/cos(DECrad))>=0.0) {
			RArad = acos((TI[0][0]*cos(brad)*cos(lrad) + TI[0][1]*cos(brad)*sin(lrad) + TI[0][2]*sin(brad))/cos(DECrad)); 
		} else {
			RArad = 2.0*PI-acos((TI[0][0]*cos(brad)*cos(lrad) + TI[0][1]*cos(brad)*sin(lrad) + TI[0][2]*sin(brad))/cos(DECrad)); 			
		}
	}
	
	
	
	//get tangential velocity in [km/s] from different sets of velocity-measurement types
	
	//get coordinates of equatorial north pole on great circle
	bENP = asin(T[2][0]*cos(DECENP)*cos(RAENP) + T[2][1]*cos(DECENP)*sin(RAENP) + T[2][2]*sin(DECENP));
	if (asin((T[1][0]*cos(DECENP)*cos(RAENP) + T[1][1]*cos(DECENP)*sin(RAENP) + T[1][2]*sin(DECENP))/cos(bENP))>=0.0) {
		lENP = acos((T[0][0]*cos(DECENP)*cos(RAENP) + T[0][1]*cos(DECENP)*sin(RAENP) + T[0][2]*sin(DECENP))/cos(bENP)); 
	} else {
		lENP = 2.0*PI-acos((T[0][0]*cos(DECENP)*cos(RAENP) + T[0][1]*cos(DECENP)*sin(RAENP) + T[0][2]*sin(DECENP))/cos(bENP)); 			
	}
	if (radiococo) printf("\nCoordinates of equatorial north pole:\n");
	if (radiococo) printf("bENP = %f\tlENP = %f\n", bENP, lENP);
	zENP = sin(bENP)*dsun;
	dxyENP = sqrt(dsun*dsun-zENP*zENP);
	xENP = cos(lENP)*dxyENP + xsun;
	yENP = sin(lENP)*dxyENP;
	if (radiococo) printf("xENP = %f\tyENP = %f\tzENP = %f\n", xENP, yENP, zENP);
	
	
	if (vcoordtype == 1) {
		
		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun*vrsun;
		vry = y/dsun*vrsun;
		vrz = z/dsun*vrsun;
		if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));		
		
		//convert to km/s
		mu = *mutemp*dsun*4.74057;
		
		//compute proper motion components
		mu_alphacosdelta = mu*sin(PArad);
		mu_delta = mu*cos(PArad);
		
		A[0][0] = cos(RArad)*cos(DECrad);
		A[0][1] = -sin(RArad);
		A[0][2] = -cos(RArad)*sin(DECrad);
		
		A[1][0] = sin(RArad)*cos(DECrad);
		A[1][1] = cos(RArad);
		A[1][2] = -sin(RArad)*sin(DECrad);
		
		A[2][0] = sin(DECrad);
		A[2][1] = 0.0;
		A[2][2] = cos(DECrad);
		
		//printf("%f\t%f\t%f\n",A[0][0],A[0][1],A[0][2]);
		//printf("%f\t%f\t%f\n",A[1][0],A[1][1],A[1][2]);
		//printf("%f\t%f\t%f\n",A[2][0],A[2][1],A[2][2]);
		
		//B = T * A
		B[0][0] = T[0][0]*A[0][0] + T[0][1]*A[1][0] + T[0][2]*A[2][0];
		B[0][1] = T[0][0]*A[0][1] + T[0][1]*A[1][1] + T[0][2]*A[2][1];
		B[0][2] = T[0][0]*A[0][2] + T[0][1]*A[1][2] + T[0][2]*A[2][2];
		
		B[1][0] = T[1][0]*A[0][0] + T[1][1]*A[1][0] + T[1][2]*A[2][0];
		B[1][1] = T[1][0]*A[0][1] + T[1][1]*A[1][1] + T[1][2]*A[2][1];
		B[1][2] = T[1][0]*A[0][2] + T[1][1]*A[1][2] + T[1][2]*A[2][2];
		
		B[2][0] = T[2][0]*A[0][0] + T[2][1]*A[1][0] + T[2][2]*A[2][0];
		B[2][1] = T[2][0]*A[0][1] + T[2][1]*A[1][1] + T[2][2]*A[2][1];
		B[2][2] = T[2][0]*A[0][2] + T[2][1]*A[1][2] + T[2][2]*A[2][2];
		
		//printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
		//printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
		//printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		
		vx = vrsun*B[0][0] + mu_alphacosdelta*B[0][1] + mu_delta*B[0][2] +vxsun;
		vy = vrsun*B[1][0] + mu_alphacosdelta*B[1][1] + mu_delta*B[1][2] +vysun+vLSRtemp;
		vz = vrsun*B[2][0] + mu_alphacosdelta*B[2][1] + mu_delta*B[2][2] +vzsun;
		
		
		if (radiococo) printf("\nCartesian velocity:\n");
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));
		
	} else if (vcoordtype == 2) {
		
		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun*vrsun;
		vry = y/dsun*vrsun;
		vrz = z/dsun*vrsun;
		if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));		
		
		if (*mu_alphatemp) *mu_alphacosdeltatemp = *mu_alphatemp*cos(DECrad);
		else if (*mu_alphacosdeltatemp) *mu_alphatemp = *mu_alphacosdeltatemp/cos(DECrad);
		
		//convert to km/s
		mu_alphacosdelta = *mu_alphacosdeltatemp*dsun*4.74057;
		mu_delta = *mu_deltatemp*dsun*4.74057;
		mu = sqrt(mu_alphacosdelta*mu_alphacosdelta+mu_delta*mu_delta);
		
		A[0][0] = cos(RArad)*cos(DECrad);
		A[0][1] = -sin(RArad);
		A[0][2] = -cos(RArad)*sin(DECrad);
		
		A[1][0] = sin(RArad)*cos(DECrad);
		A[1][1] = cos(RArad);
		A[1][2] = -sin(RArad)*sin(DECrad);
		
		A[2][0] = sin(DECrad);
		A[2][1] = 0.0;
		A[2][2] = cos(DECrad);
		
		//printf("%f\t%f\t%f\n",A[0][0],A[0][1],A[0][2]);
		//printf("%f\t%f\t%f\n",A[1][0],A[1][1],A[1][2]);
		//printf("%f\t%f\t%f\n",A[2][0],A[2][1],A[2][2]);
		
		//B = T * A
		B[0][0] = T[0][0]*A[0][0] + T[0][1]*A[1][0] + T[0][2]*A[2][0];
		B[0][1] = T[0][0]*A[0][1] + T[0][1]*A[1][1] + T[0][2]*A[2][1];
		B[0][2] = T[0][0]*A[0][2] + T[0][1]*A[1][2] + T[0][2]*A[2][2];
		
		B[1][0] = T[1][0]*A[0][0] + T[1][1]*A[1][0] + T[1][2]*A[2][0];
		B[1][1] = T[1][0]*A[0][1] + T[1][1]*A[1][1] + T[1][2]*A[2][1];
		B[1][2] = T[1][0]*A[0][2] + T[1][1]*A[1][2] + T[1][2]*A[2][2];
		
		B[2][0] = T[2][0]*A[0][0] + T[2][1]*A[1][0] + T[2][2]*A[2][0];
		B[2][1] = T[2][0]*A[0][1] + T[2][1]*A[1][1] + T[2][2]*A[2][1];
		B[2][2] = T[2][0]*A[0][2] + T[2][1]*A[1][2] + T[2][2]*A[2][2];
		
		//printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
		//printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
		//printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		
		vx = vrsun*B[0][0] + mu_alphacosdelta*B[0][1] + mu_delta*B[0][2] +vxsun;
		vy = vrsun*B[1][0] + mu_alphacosdelta*B[1][1] + mu_delta*B[1][2] +vysun+vLSRtemp;
		vz = vrsun*B[2][0] + mu_alphacosdelta*B[2][1] + mu_delta*B[2][2] +vzsun;
		
		if (radiococo) printf("\nCartesian velocity:\n");
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));
		
		//get position angle of proper motion
		
		//heliocentric transverse velocity
		vtx = vx-vxsun-vrx;
		vty = vy-vysun-vLSRtemp-vry;
		vtz = vz-vzsun-vrz;
		if (radiococo) printf("\nTransverse velocity:\n");
		if (radiococo) printf("vtx = %f\tvty = %f\tvtz = %f\tvt = %f [km/s] (heliocentric)\n", vtx, vty, vtz, sqrt(vtx*vtx+vty*vty+vtz*vtz));
		
		//get tangential vector pointing to ENP
		FAK = -((xENP-xsun)*(x-xsun)+yENP*y+zENP*z)/(pow(x-xsun,2)+y*y+z*z);
		xdelta = FAK*(x-xsun)+(xENP-xsun);
		ydelta = FAK*y+yENP;
		zdelta = FAK*z+zENP;
		
		//determine distance (pos or neg) of Xobject + Vt from plane connecting ENP, Xobject and observer for position angle
		nx = y*zENP-z*yENP;
		ny = z*(xENP-xsun)-(x-xsun)*zENP;
		nz = (x-xsun)*yENP-y*(xENP-xsun);
		dvt = nx*(x+vtx)+ny*(y+vty)+nz*(z+vtz)-nx*xsun;
		
		//get position angle of proper motion with respect to tangential vector pointing to ENP
		if (dvt <= 0) 
			PArad = acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		else 
			PArad = 2.0*PI-acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		
		if (radiococo) printf("\nProper motion and position angle:\n");
		if (radiococo) printf("mu = %f\tPA = %f\n", mu, PArad);
		
	} else if (vcoordtype == 3) {
		
		if (radiococo) printf("\nCartesian velocity:\n");
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));		
		
		//heliocentric radial velocity
		vrsun = ((vx-vxsun)*(x-xsun)+(vy-vysun-vLSRtemp)*y+(vz-vzsun)*z)/sqrt(pow(x-xsun,2)+y*y+z*z);
		
		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun*vrsun;
		vry = y/dsun*vrsun;
		vrz = z/dsun*vrsun;
		if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));				
		
		//get position angle of proper motion
		
		//heliocentric transverse velocity
		vtx = vx-vxsun-vrx;
		vty = vy-vysun-vLSRtemp-vry;
		vtz = vz-vzsun-vrz;
		if (radiococo) printf("\nTransverse velocity:\n");
		if (radiococo) printf("vtx = %f\tvty = %f\tvtz = %f\tvt = %f [km/s] (heliocentric)\n", vtx, vty, vtz, sqrt(vtx*vtx+vty*vty+vtz*vtz));
		
		//get tangential vector pointing to ENP
		FAK = -((xENP-xsun)*(x-xsun)+yENP*y+zENP*z)/(pow(x-xsun,2)+y*y+z*z);
		xdelta = FAK*(x-xsun)+(xENP-xsun);
		ydelta = FAK*y+yENP;
		zdelta = FAK*z+zENP;
		
		//determine distance (pos or neg) of Xobject + Vt from plane connecting ENP, Xobject and observer for position angle
		nx = y*zENP-z*yENP;
		ny = z*(xENP-xsun)-(x-xsun)*zENP;
		nz = (x-xsun)*yENP-y*(xENP-xsun);
		dvt = nx*(x+vtx)+ny*(y+vty)+nz*(z+vtz)-nx*xsun;
		
		//get position angle of proper motion with respect to tangential vector pointing to ENP
		if (dvt <= 0) 
			PArad = acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		else 
			PArad = 2.0*PI-acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		
		if (radiococo) printf("\nProper motion and position angle:\n");
		if (radiococo) printf("mu = %f\tPA = %f\n", mu, PArad);
		
		mu = sqrt(vtx*vtx+vty*vty+vtz*vtz);
		mu_delta = mu*cos(PArad);
		mu_alphacosdelta = mu*sin(PArad);
		
	}
	
	
	if (radiococo) printf("\nProper motion:\n");
	if (radiococo) printf("mu_alphacosdelta  = %f\tmu_delta = %f\tmu = %f [km/s]\t PA = %f\n", mu_alphacosdelta, mu_delta, mu, PArad);
	
	vr = (vx*(x-xsun)+vy*y+vz*z)/sqrt(pow(x-xsun,2)+y*y+z*z);
	if (radiococo) printf("\nRadial velocity:\n");
	if (radiococo) printf("vr = %.3f\tvr = %.3f (heliocentric) [km/s]\n", vr, vrsun);
	
	//consistency check with formula for GSR radial velocity from script of Steven Majewski
	if (radiococo) vrLSR = vrsun + (vxsun*cos(brad)*cos(lrad)+vysun*cos(brad)*sin(lrad)+vzsun*sin(brad));
	if (radiococo) vrGSR = vrLSR + vLSRtemp*cos(brad)*sin(lrad);
	if (radiococo) printf("\nConsistency check with formula for Galactic standard of rest (GSR) radial velocity (should be equal to vr):\n");
	if (radiococo) printf("vr_LSR = %f\tvr_GSR = %.3f [km/s]\n", vrLSR, vrGSR);
	
	
	
	//convert back to input units and write to output
	*xtemp = 1000.0*x;
	*(xtemp+1) = 1000.0*y;
	*(xtemp+2) = 1000.0*z;
	
	*dsuntemp = 1000.0*dsun;
	
	*vtemp = vx;
	*(vtemp+1) = vy;
	*(vtemp+2) = vz;
	
	*vrsuntemp = vrsun;
	*vrtemp = vr;
	
	*DECtemp = DECrad*180.0/PI;
	*RAtemp = RArad*180.0/PI;
	
	*btemp = brad*180.0/PI;
	*ltemp = lrad*180.0/PI;
	*lcosbtemp = *ltemp*cos(brad);
	
	*mutemp = mu/(dsun*4.74057);
	*PAtemp = PArad*180.0/PI;
	*mu_deltatemp = mu_delta/(dsun*4.74057);
	*mu_alphacosdeltatemp = mu_alphacosdelta/(dsun*4.74057);
	*mu_alphatemp = *mu_alphacosdeltatemp/cos(DECrad);
	
}

int readin0(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt) {  //McLuster;
	
	int i, k;
	double m, x[3], v[3], m0, epoch, spin, rad, lum; 
	
	i=0;
	while(fscanf(posfile,"%lf %lf %lf %lf %lf %lf %lf %lf %i %lf %lf %lf %lf\n",&m,&x[0],&x[1],&x[2],&v[0],&v[1],&v[2],&m0, &k, &epoch, &spin, &rad, &lum) == 13) {
		start[i][0] = i+1;	
		start[i][1] = x[0];              //in pc
		start[i][2] = x[1];
		start[i][3] = x[2];
		start[i][4] = v[0];             //in km/s
		start[i][5] = v[1];
		start[i][6] = v[2];
		start[i][7] = k;
		start[i][8] = lum;
		start[i][9] = rad;
		start[i][10] = m;                 //in Sonnenmassen
		start[i][12] = 0.5*m*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));  //Ekin
		
		i++;
	}
	
	RGt[0] = 8500.0;
	RGt[1] = 0.0;
	RGt[2] = 0.0;
	VGt[0] = 0.0;
	VGt[1] = 220.0;
	VGt[2] = 0.0;
	*omegat = 220.0/8500.0;
	
	*timet = 0.0;
	*Nt = i;
	*Ntott = i;
	*rtidet = 0.0;
	
	*Etott = 0;
	*rcoret = 0;
	
	return 1;
}

int readin1(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, int outputtypet, int Nmaxt) {  //NBODY4 32bit;
	int l;
	
	int a[6];
	a[0] = 888;
	a[1] = -111;
	a[2] = 999;
	a[3] = -222;
	a[4] = 777;
	a[5] = -333;
	int b;
	do {
		fread(&b,sizeof(int),1,posfile);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (outputtypet) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[4]);
	
	double as[a[2]];
	fread(as,sizeof(double),a[2],posfile);
	
	if (outputtypet) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posfile);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posfile);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posfile);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	fread(phi4,sizeof(double),a[1],posfile);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posfile);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posfile); 
	double *lsev;
	lsev = (double *)malloc(a[1]*sizeof(double));
	fread(lsev,sizeof(float),a[1],posfile);
	double *rsev;
	rsev = (double *)malloc(a[1]*sizeof(double)); 
	fread(rsev,sizeof(float),a[1],posfile);
	
	*tbart = as[4];
	*timet = as[9];
	*Ntott = a[1];
	*Nt = a[4];
	if (as[24])
		*rtidet = as[24]*as[2];
	else
		*rtidet = 0.0;
	
	*rbart = as[2];
	*vbart = as[11];
	*mbart = as[3];
	if (as[25])
		*omegat = as[25]/as[4];
	else 
		*omegat = 0.0;
	
	*Etott = as[18];
	*rcoret = as[12];
	
	/***********************************
	 * Einzelsterne in star schreiben *
	 ***********************************/
	
	for(l=0;l<a[4];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < Nmaxt)) {
			start[name[l]-1][0] = name[l];	
			start[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			start[name[l]-1][2] = xs[3*l+1]*as[2];
			start[name[l]-1][3] = xs[3*l+2]*as[2];
			start[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			start[name[l]-1][5] = vs[3*l+1]*as[11];
			start[name[l]-1][6] = vs[3*l+2]*as[11];
			start[name[l]-1][7] = kstar[l];
			start[name[l]-1][8] = lsev[l];
			start[name[l]-1][9] = rsev[l];
			start[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			start[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
		} else {
			printf("\nFUNNY: name[%i] - 1 = %i\n", l, name[l]-1);
			*Ntott -= 1;
			*Nt -= 1;
		}
	}
	
	
	/*****************************
	 * Speicher wieder freigeben *
	 *****************************/
	
	free(bodys);
	free(xs);
	free(vs);
	free(name);
	free(kstar);
	free(phi4);
	free(lsev);
	free(rsev);
	
	return 1;
}

int readin2(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt) {  //NBODY6;
	int l;	
	int a[6];
	a[0] = 888;
	a[1] = -111;
	a[2] = 999;
	a[3] = -222;
	a[4] = 777;
	a[5] = -333;
	int b;
	do {
		fread(&b,sizeof(int),1,posfile);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (outputtypet) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[4]);
	
	double as[a[2]];
	fread(&as,sizeof(double),a[2],posfile);
	
	if (outputtypet) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nNULL\n");
	fread(bodys,sizeof(double),a[1],posfile);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posfile);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posfile);
	double *rsev;
	rsev = (double *)malloc(a[1]*sizeof(double)); 
	fread(rsev,sizeof(double),a[1],posfile);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posfile);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posfile); 
	double *lsev;
	lsev = (double *)malloc(a[1]*sizeof(double));
	fread(lsev,sizeof(double),a[1],posfile);
	int end[1];
	fread(&end,sizeof(int),1,posfile);
	//fread(&b,sizeof(int),1,posfile);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	//fread(phi4,sizeof(double),a[1],posfile);
	
    /*
    AS(1) = TTOT
    AS(2) = FLOAT(NPAIRS)
    AS(3) = RBAR
    AS(4) = ZMBAR
    AS(5) = RTIDE
    AS(6) = TIDAL(4)
    AS(7) = RDENS(1)
    AS(8) = RDENS(2)
    AS(9) = RDENS(3)
    AS(10) = TSCALE*TTOT
    AS(11) = TSCALE
    AS(12) = VSTAR
    AS(13) = RC
    AS(14) = NC
    AS(15) = VC
    AS(16) = RHOM
    AS(17) = CMAX
    AS(18) = RSCALE
    AS(19) = RSMIN
    AS(20) = DMIN1
    AS(21) = RG(1)
    AS(22) = RG(2)
    AS(23) = RG(3)
    AS(24) = GMG
    AS(25) = ZDUM(1)
    AS(26) = OMEGA
    AS(27) = VG(1)
    AS(28) = VG(2)
    AS(29) = VG(3)
    AS(30) = DISK
    AS(31) = A
    AS(32) = B
    AS(33) = V02
    AS(34) = RL2
    AS(35) = GMB
    AS(36) = AR
    AS(37) = GAM
    AS(38) = ZDUM(2)
    AS(39) = ZDUM(3)
    AS(40) = ZDUM(4)
*/
    
	*rbart = as[2];
	*vbart = as[11];
	*mbart = as[3];
	*tbart = as[10];
	*timet = as[9];
	*Ntott = a[1];
	*Nt = a[3];
    
	RGt[0] = as[20]**rbart;
	RGt[1] = as[21]**rbart;
	RGt[2] = as[22]**rbart;
	VGt[0] = as[26]**vbart;
	VGt[1] = as[27]**vbart;
	VGt[2] = as[28]**vbart;

    M1 = as[23]**mbart;
    b1 = as[24]**rbart;
    if (b1 <= 0.0) b1 = 1.0;
    
    M1_GAMMA = as[34]**mbart;
    b1_GAMMA = as[35]**rbart;
    if (b1_GAMMA <= 0.0) b1_GAMMA = 1.0;
    
    M2 = as[29]**mbart;
    a2 = as[30]**rbart;
    b2 = as[31]**rbart;
    if (a2 <= 0.0) a2 = 1.0;
    if (b2 <= 0.0) b2 = 1.0;

    q_halo = as[39];
    
    VCIRC = sqrt(as[32]**vbart);
    RCIRC = sqrt(as[33]**rbart);
    if (RCIRC <= 0.0) RCIRC = 1.0;

    MNFW = as[37]**mbart;
    RNFW = as[38]**rbart;
    if (RNFW <= 0.0) RNFW = 1.0;

    
	/************************************
	 * write single stars to star array *
	 ************************************/
	
	for(l=0;l<a[3];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < Nmaxt)) {
			start[name[l]-1][0] = name[l];	
			start[name[l]-1][1] = xs[3*l]*as[2];              //convert to pc
			start[name[l]-1][2] = xs[3*l+1]*as[2];
			start[name[l]-1][3] = xs[3*l+2]*as[2];
			start[name[l]-1][4] = vs[3*l]*as[11];             //convert to (km/s)^2
			start[name[l]-1][5] = vs[3*l+1]*as[11];
			start[name[l]-1][6] = vs[3*l+2]*as[11];
			start[name[l]-1][7] = kstar[l];
			start[name[l]-1][8] = lsev[l];
			start[name[l]-1][9] = rsev[l];
			start[name[l]-1][10] = bodys[l]*as[3];                 //convert to solar masses
			start[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
		} else {
			printf("\nFUNNY: name[%i] - 1 = %i\n", l, name[l]-1);
			*Ntott -= 1;
			*Nt -= 1;
		}
	}
	
		
	free(bodys);
	free(xs);
	free(vs);
	free(name);
	free(kstar);
	free(phi4);
	free(lsev);
	free(rsev);
	
	return 1;
}

int readin3(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt) {  //NBODY6;
	int l;	
	
	int a[6];
	a[0] = 888;
	a[1] = -111;
	a[2] = 999;
	a[3] = -222;
	a[4] = 777;
	a[5] = -333;
	int b;
	do {
		fread(&b,sizeof(int),1,posfile);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (outputtypet) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[4]);
	
	a[2] = 30; //30 output parameter
	double as[a[2]];
	fread(&as,sizeof(double),a[2],posfile);
	
	if (outputtypet) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posfile);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posfile);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posfile);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	fread(phi4,sizeof(double),a[1],posfile);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posfile);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posfile); 
	float *lsev;
	lsev = (float *)malloc(a[1]*sizeof(float));
	fread(lsev,sizeof(float),a[1],posfile);
	float *rsev;
	rsev = (float *)malloc(a[1]*sizeof(float)); 
	fread(rsev,sizeof(float),a[1],posfile);
	
	//assume circular orbit at solar radius
	RGt[0] = 8500.0;
	RGt[1] = 0.0;
	RGt[2] = 0.0;
	VGt[0] = 0.0;
	VGt[1] = 220.0;
	VGt[2] = 0.0;
	
	*tbart = as[10];
	*timet = as[9];
	*Ntott = a[1];
	*Nt = a[3];
	
	*rtidet = 0.0;
	
	*rbart = as[2];
	*vbart = as[11];
	*mbart = as[3];
	*omegat = 220.0/8500.0;
	
	*Etott = 0;
	*rcoret = 0;
	
	/***********************************
	 * Einzelsterne in star schreiben *
	 ***********************************/
	
	for(l=0;l<a[3];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < Nmaxt)) {
			start[name[l]-1][0] = name[l];	
			start[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			start[name[l]-1][2] = xs[3*l+1]*as[2];
			start[name[l]-1][3] = xs[3*l+2]*as[2];
			start[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			start[name[l]-1][5] = vs[3*l+1]*as[11];
			start[name[l]-1][6] = vs[3*l+2]*as[11];
			start[name[l]-1][7] = kstar[l];
			start[name[l]-1][8] = lsev[l];
			start[name[l]-1][9] = rsev[l];
			start[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			start[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
		} else {
			printf("\nFUNNY: name[%i] - 1 = %i\n", l, name[l]-1);
			*Ntott -= 1;
			*Nt -= 1;
		}
	}
	
	
	/*****************************
	 * Speicher wieder freigeben *
	 *****************************/
	
	free(bodys);
	free(xs);
	free(vs);
	free(name);
	free(kstar);
	free(phi4);
	free(lsev);
	free(rsev);
	
	return 1;
}

int readin4(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt) {  //NBODY6;
	int l;	
	int a[6], b[2];
	fread(&b,sizeof(int),1,posfile);
	fread(&a,sizeof(int),6,posfile);
	if (outputtypet) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	fread(&b,sizeof(int),2,posfile);
	//fread(&b,sizeof(int),1,posfile);
	double as[a[2]];
	fread(&as,sizeof(double),a[2],posfile);
	
	if (outputtypet) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posfile);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posfile);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posfile);
	double *rsev;
	rsev = (double *)malloc(a[1]*sizeof(double)); 
	fread(rsev,sizeof(double),a[1],posfile);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posfile);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posfile); 
	double *lsev;
	lsev = (double *)malloc(a[1]*sizeof(double));
	fread(lsev,sizeof(double),a[1],posfile);
	int end[1];
	fread(&end,sizeof(int),1,posfile);
	fread(&b,sizeof(int),1,posfile);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	//fread(phi4,sizeof(double),a[1],posfile);
	
	RGt[0] = as[20]*as[2];
	RGt[1] = as[21]*as[2];
	RGt[2] = as[22]*as[2];
	VGt[0] = as[26]*as[11];
	VGt[1] = as[27]*as[11];
	VGt[2] = as[28]*as[11];
	
	*tbart = as[10];
	*timet = as[9];
	*Ntott = a[1];
	*Nt = a[3];
	if (as[24])
		*rtidet = as[24]*as[2];
	else
		*rtidet = 0.0;
	
	*rbart = as[2];
	*vbart = as[11];
	*mbart = as[3];
	if (as[25])
		*omegat = as[25];
	else 
		*omegat = 220.0/8500.0;
	
	*Etott = 0;
	*rcoret = 0;
	
	/***********************************
	 * Einzelsterne in star schreiben *
	 ***********************************/
	
	for(l=0;l<a[3];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < Nmaxt)) {
			start[name[l]-1][0] = name[l];	
			start[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			start[name[l]-1][2] = xs[3*l+1]*as[2];
			start[name[l]-1][3] = xs[3*l+2]*as[2];
			start[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			start[name[l]-1][5] = vs[3*l+1]*as[11];
			start[name[l]-1][6] = vs[3*l+2]*as[11];
			start[name[l]-1][7] = kstar[l];
			start[name[l]-1][8] = lsev[l];
			start[name[l]-1][9] = rsev[l];
			start[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			start[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
		} else {
			printf("\nFUNNY: name[%i] - 1 = %i\n", l, name[l]-1);
			*Ntott -= 1;
			*Nt -= 1;
		}
	}
	
	
	/*****************************
	 * Speicher wieder freigeben *
	 *****************************/
	
	free(bodys);
	free(xs);
	free(vs);
	free(name);
	free(kstar);
	free(phi4);
	free(lsev);
	free(rsev);
	
	return 1;
}

int readin5(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt) {  //NBODY4new 64bit;
	int l;	
	
	int a[6];
	a[0] = 888;
	a[1] = -111;
	a[2] = 999;
	a[3] = -222;
	a[4] = 777;
	a[5] = -333;
	int b;
	do {
		fread(&b,sizeof(int),1,posfile);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (outputtypet) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[4]);
	
	double as[a[2]];
	fread(&as,sizeof(double),a[2],posfile);
	
	if (outputtypet) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posfile);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posfile);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posfile);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	fread(phi4,sizeof(double),a[1],posfile);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posfile);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posfile); 
	double *lsev;
	lsev = (double *)malloc(a[1]*sizeof(double));
	fread(lsev,sizeof(double),a[1],posfile);
	double *rsev;
	rsev = (double *)malloc(a[1]*sizeof(double)); 
	//fread(rsev,sizeof(float),a[1],posfile);
	//int end[1];
	//fread(&end,sizeof(int),1,posfile);
	//fread(&b,sizeof(int),1,posfile);
	
	RGt[0] = as[20]*as[2];
	RGt[1] = as[21]*as[2];
	RGt[2] = as[22]*as[2];
	VGt[0] = as[26]*as[11];
	VGt[1] = as[27]*as[11];
	VGt[2] = as[28]*as[11];
	
	*tbart = as[10];
	*timet = as[9];
	*Ntott = a[1];
	*Nt = a[3];
	if (as[24])
		*rtidet = as[24]*as[2];
	else
		*rtidet = 0.0;
	
	*rbart = as[2];
	*vbart = as[11];
	*mbart = as[3];
	if (as[25])
		*omegat = as[25];
	else 
		*omegat = 0.0;
	
	*Etott = 0;
	*rcoret = 0;
	
	/***********************************
	 * Einzelsterne in star schreiben *
	 ***********************************/
	
	for(l=0;l<a[3];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < Nmaxt)) {
			start[name[l]-1][0] = name[l];	
			start[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			start[name[l]-1][2] = xs[3*l+1]*as[2];
			start[name[l]-1][3] = xs[3*l+2]*as[2];
			start[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			start[name[l]-1][5] = vs[3*l+1]*as[11];
			start[name[l]-1][6] = vs[3*l+2]*as[11];
			start[name[l]-1][7] = kstar[l];
			start[name[l]-1][8] = lsev[l];
			start[name[l]-1][9] = rsev[l];
			start[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			start[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
		} else {
			printf("\nFUNNY: name[%i] - 1 = %i\n", l, name[l]-1);
			*Ntott -= 1;
			*Nt -= 1;
		}
	}
	
	
	/*****************************
	 * Speicher wieder freigeben *
	 *****************************/
	
	free(bodys);
	free(xs);
	free(vs);
	free(name);
	free(kstar);
	free(phi4);
	free(lsev);
	free(rsev);
	
	return 1;
}

int readin6(FILE *posfile, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **start, double *RGt, double *VGt, int outputtypet, int Nmaxt) {  //NBODY6;
	int l;	
	int a[6];
	a[0] = 888;
	a[1] = -111;
	a[2] = 999;
	a[3] = -222;
	a[4] = 777;
	a[5] = -333;
	int b;
	do {
		fread(&b,sizeof(int),1,posfile);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (outputtypet) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[5]);
	fread(&b,sizeof(int),1,posfile);
	fread(&b,sizeof(int),1,posfile);
	fread(&b,sizeof(int),1,posfile);
	
	double as[a[3]];
	fread(&as,sizeof(double),a[3],posfile);
	
	if (outputtypet) for (l=0;l<a[3];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[2]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[2],posfile);
	double *xs;
	xs = (double *)malloc(3*a[2]*sizeof(double));
	fread(xs,sizeof(double),3*a[2],posfile);
	double *vs;
	vs = (double *)malloc(3*a[2]*sizeof(double));
	fread(vs,sizeof(double),3*a[2],posfile);
	double *rsev;
	rsev = (double *)malloc(a[2]*sizeof(double)); 
	fread(rsev,sizeof(double),a[2],posfile);
	int *name;
	name = (int *)malloc(a[2]*sizeof(int));
	fread(name,sizeof(int),a[2],posfile);
	int *kstar;
	kstar = (int *)malloc(a[2]*sizeof(int));
	fread(kstar,sizeof(int),a[2],posfile); 
	double *lsev;
	lsev = (double *)malloc(a[2]*sizeof(double));
	fread(lsev,sizeof(double),a[2],posfile);
	int end[1];
	fread(&end,sizeof(int),1,posfile);
	//fread(&b,sizeof(int),1,posfile);
	double *phi4;
	phi4 = (double *)malloc(a[2]*sizeof(double));
	//fread(phi4,sizeof(double),a[2],posfile);
	
    /*
     AS(1) = TTOT
     AS(2) = FLOAT(NPAIRS)
     AS(3) = RBAR
     AS(4) = ZMBAR
     AS(5) = RTIDE
     AS(6) = TIDAL(4)
     AS(7) = RDENS(1)
     AS(8) = RDENS(2)
     AS(9) = RDENS(3)
     AS(10) = TSCALE*TTOT
     AS(11) = TSCALE
     AS(12) = VSTAR
     AS(13) = RC
     AS(14) = NC
     AS(15) = VC
     AS(16) = RHOM
     AS(17) = CMAX
     AS(18) = RSCALE
     AS(19) = RSMIN
     AS(20) = DMIN1
     AS(21) = RG(1)
     AS(22) = RG(2)
     AS(23) = RG(3)
     AS(24) = GMG
     AS(25) = ZDUM(1)
     AS(26) = OMEGA
     AS(27) = VG(1)
     AS(28) = VG(2)
     AS(29) = VG(3)
     AS(30) = DISK
     AS(31) = A
     AS(32) = B
     AS(33) = V02
     AS(34) = RL2
     AS(35) = GMB
     AS(36) = AR
     AS(37) = GAM
     AS(38) = ZDUM(2)
     AS(39) = ZDUM(3)
     AS(40) = ZDUM(4)
     */
    
	*rbart = as[2];
	*vbart = as[11];
	*mbart = as[3];
	*tbart = as[10];
	*timet = as[9];
	*Ntott = a[1];
	*Nt = a[3];
    
	RGt[0] = as[20]**rbart;
	RGt[1] = as[21]**rbart;
	RGt[2] = as[22]**rbart;
	VGt[0] = as[26]**vbart;
	VGt[1] = as[27]**vbart;
	VGt[2] = as[28]**vbart;
    
    M1 = as[23]**mbart;
    b1 = as[24]**rbart;
    if (b1 <= 0.0) b1 = 1.0;
    
    M1_GAMMA = as[34]**mbart;
    b1_GAMMA = as[35]**rbart;
    if (b1_GAMMA <= 0.0) b1_GAMMA = 1.0;
    
    M2 = as[29]**mbart;
    a2 = as[30]**rbart;
    b2 = as[31]**rbart;
    if (a2 <= 0.0) a2 = 1.0;
    if (b2 <= 0.0) b2 = 1.0;
    
    q_halo = as[39];
    
    VCIRC = sqrt(as[32]**vbart);
    RCIRC = sqrt(as[33]**rbart);
    if (RCIRC <= 0.0) RCIRC = 1.0;
    
    MNFW = as[37]**mbart;
    RNFW = as[38]**rbart;
    if (RNFW <= 0.0) RNFW = 1.0;

	
	/************************************
	 * write single stars to star array *
	 ************************************/
	
	for(l=0;l<a[4];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < Nmaxt)) {
			start[name[l]-1][0] = name[l];	
			start[name[l]-1][1] = xs[3*l]*as[2];              //convert to pc
			start[name[l]-1][2] = xs[3*l+1]*as[2];
			start[name[l]-1][3] = xs[3*l+2]*as[2];
			start[name[l]-1][4] = vs[3*l]*as[11];             //convert to (km/s)^2
			start[name[l]-1][5] = vs[3*l+1]*as[11];
			start[name[l]-1][6] = vs[3*l+2]*as[11];
			start[name[l]-1][7] = kstar[l];
			start[name[l]-1][8] = lsev[l];
			start[name[l]-1][9] = rsev[l];
			start[name[l]-1][10] = bodys[l]*as[3];                 //convert to solar masses
			start[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
		} else {
			printf("\nFUNNY: name[%i] - 1 = %i\n", l, name[l]-1);
			*Ntott -= 1;
			*Nt -= 1;
		}
	}
	
    
	free(bodys);
	free(xs);
	free(vs);
	free(name);
	free(kstar);
	free(phi4);
	free(lsev);
	free(rsev);
	
	return 1;
}
