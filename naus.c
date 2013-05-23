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
 *  Kubus array description:                                               *
 *   kubus[x][0]    -  name                                                *
 *   kubus[x][1-3]  -  x, y, z [pc]                                        *
 *   kubus[x][4-6]  -  vx, vy, vz [km/s]                                   *
 *   kubus[x][7]    -  stellar type (see Nbody documentation)              *
 *   kubus[x][8]    -  luminosity [L_sun]                                  *
 *   kubus[x][9]    -  radius [R_sun]                                      *
 *   kubus[x][10]   -  mass [M_sun]                                        *
 *   kubus[x][11]   -  potential energy [M_sun km^2/s^2]                   *
 *   kubus[x][12]   -  kinetic energy [M_sun km^2/s^2]                     *
 *   kubus[x][13]   -  status (see subroutine statussetzer)                *
 *   kubus[x][14]   -  escape status (to avoid multiple escapes)           *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<sys/stat.h>



/***************
 * DEFINITIONS *
 ***************/

#define G 0.0043009211 //in pc*km^2/s^2*M_sun
#define Pi 3.14159265
#define PI 3.14159265
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

//Allen & Santillan potential - constants:
#define b1allen 387.3            //pc
#define M1allen 606.0*2.32e07   //9.565439E+10    //[M_sun]

#define a2allen 5317.8
#define b2allen 250.0
#define M2allen 3690.0*2.32e07

#define a3allen 12000.0
#define M3allen 4615.0*2.32e07

//point-mass potential - constants:
#define M1pointmass 1.0E+11//9.565439E+10    //[M_sun]

//double const vxsun = 10.0;    //km/s; solar motion with respect to the LSR from Dehnen & Binney (1998)
//double const vysun = 5.3;
//double const vzsun = 7.2;
double const vxsun = 11.1;//+0.69/−0.75 - solar motion with respect to the LSR from Schoenrich, Binney & Dehnen (2010) [km/s]
double const vysun = 12.24;//+0.47−0.47 
double const vzsun = 7.25;//+0.37−0.36

double const alphaGNP = 192.859508; //Galactic north pole in J2000 coordinates
double const deltaGNP = 27.128336;
double const PAGNP = 122.932;       //Position angle with respect to equatorial pole

double const vLSR = 220.0;      //circular velocity at solar radius
double const rgalsun = 8330.0;  //galactocentric radius of the sun [pc]



/****************
 * SUBROUNTINES *
 ****************/

int cmpmy(float *x1, float *x2);
int cmpmy2(float *x1, float *x2);
int readin0(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott);
int readin1(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, int ausgabetypt, int kubussizetott);
int readin2(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott);
int readin3(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott);
int readin4(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott);
int readin5(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott);
int readin6(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott);
void shellsort_reverse_1d(double *array, int N);
void shellsort_1d(double *array, int N);	
void shellsort(double **array, int N, int k);
void shellsort_reverse(double **array, int N, int k);
int bindungsanalyse(int Nt, double **kubust, int **nachbarlistet, int nachbarnt);
int bindungsanalysefast(int Nt, double **kubust, int **nachbarlistet, int nachbarnt, double mxt, double myt, double mzt);
int zentrumsuche(int Nt, double **kubust, int **nachbarlistet, double *mxt, double *myt, double *mzt, int densitylevel);
int omegasuche(double *RGt, double *VGt, double *omegat, double *Lt);
int rtidesuche(double *rtidet, int Nt, double **kubust, double omegat, double mxt, double myt, double mzt, double dphit, double *RGt);
int statussetzer(int Nt, double **kubust, double rtidet, double mxt, double myt, double mzt, int *nt, int *nbint, double mvxt, double mvyt, double mvzt, double rh, double omega, int schnellanalyse);
int doppelsternpruefer(int Nt, double **kubust, int **nachbarlistet, int nachbarnt, int *Ntott, int *multiplicityarray, int binoutput, FILE *bint, double timet, int ausgabetypt);
int geschwindigkeitsdispersionsberechnungsfunktion(int Nt,double **kubust,int statust, double *sigmat, double *mvxt, double *mvyt, double *mvzt);
int energieberechner(int Nt, double **kubust, int statust, double *Etemp, double mvxt, double mvyt, double mvzt);
int radien(int Nt, double **kubust, int statust, double *LRt, double mxt, double myt, double mzt, double *Mt);
int Kradien(int Nt, double **kubust, int statust, double *LRt, double mxt, double myt, double mzt, double *Mt, double *racc);
int surfacedensityprofile(int Nt, double **kubust, double **profilet, double mxt, double myt, double mzt, double mvxt, double mvyt, double mvzt, int viewt, int anzahlradialbinst, double *VGt, double **profile_inside_rtidet, double **profile_bound_starst);
int orbitausrechner(double *RGt,double *VGt,double **entlangdesorbitst,int pcanzahlintervallet,int vorzt, int potenzialt);	
void getforce(double *x, double *v, double *a, int potenzialt);
void do_step(double dt, double *x, double *v, int potenzialt);
int rcsuche(int Nt, double **kubust, int **nachbarlistet, double mxt, double myt, double mzt, int densitylevelt, double *Rct);	
int dphiberechner(double *RGt, double *VGt, double omegat, double *dphit);
int escape(FILE *esc, int Nt, double **kubust, double rtidet, double mxt, double myt, double mzt, double mvxt, double mvyt, double mvzt, double omegat, double timet, double sigmat, double kTt, double Mt, int ausgabetypt, double **timebint, int timebinsizet, double thresholdt);
void convert(double *x, double *v, double *dsun, double *vrsun, double *vr, double *l, double *b, double *lcosb, double *RA, double *DEC, double *mu_alpha, double *mu_alphacosdelta, double *mu_delta, double *mutemp, double *PAtemp, int coordtype, int vcoordtype, int radiococo, double vLSRtemp);





/*****************
 * MAIN ROUNTINE *
 *****************/

int main(int argc , char * argv[]){
	
	/**********************
	 * CONTROL PARAMETERS *
	 **********************/
	
	int code = 2;                         //0: McLuster input; 1: NBODY4 (32bit); 2: NBODY6; 3: Holgers NBODY6; 4: Elhams NBODY6; 5: NBODY4new (64bit); 6: NBODY6 (Holger's cluster)
	int kubussizetot = 300000;            //size of data array, minimum is twice the initial star number
	int Mende = 0;	                      //stop analysis if M<Mende [M_sun]
	int analysistimestep = 1;             //time step for analysis [Myr]; 1: all available time steps
	int schnellanalyse = 0;               //calculation of binding energy and neighbours; 0: detailed; 1: fast
	int standarddensitylevel = 5;	      //Casertano & Hut density estimator level (j) (standard = 5)
	int standardnachbarzahl	= 20;         //length of neighbour list [20 is a reasonable size]
	int standardstatus = 2;	              //calculate basic cluster parameters by taking into account: 1: all bound stars within rtide; 2: all stars within rtide; 3: all stars within 2*rtide
	int binoutput = 0;                    //output of binary data; 0: summary; 1: detailed
	int ausgabetyp = 0;                   //0: output of major quantities to STDOUT; 1: debug mode
	int ausgabeschritt = 1;               //output interval for all.txt, all.bin, .sdp, y.bin and y.txt [Myr]; 1: every available time step
	int allausgabe = 0;                   //0: output of all stars to all.txt and all.bin off; 1: on
	int tailanalyse = 0;                  //0: analysis and output of tidal-tail stars to y.bin and y.txt off; 1: on
	int potenzial = 2;                    //0: near-field approximation; 1: point-mass potential; 2: full galactic potential; -1: isolated	
	int sunviewtime = 12000;              //output interval for snapshots from the position of the sun [Myr]
	int sdp = 0;						  //0: output of surface-density profiles to .sdp off; 1: on
	int view = 1;                         //projection of surface-density profiles; 1: xy-plane; 2: xz-plane
	int binarysearch = 0;				 //search for binaries and replace them by centre-of-mass particles? 0: no; 1: yes
	int trackescapers = 0;                //search for escapers and determine escape velocities; 0: off; 1: on
	int timebinsize = 100;				  //size of time bins for dynamical temperature-luminosity diagram [Myr] (standard = 100)
	double threshold = 5.1;				  //threshold between evaporaters and ejecters (see Kuepper et al. 2008 eq. 10)
	double limit = 0.0;					 //detection limit for stellar masses in snapshots [Msun]
	double vlimit = 100000.0;//26.0			//detection limit in V magnitudes
	double RG0[3] = {0.0, 8500.0, 0.0};	  //position of cluster if not provided by input file
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
	
	printf("Anzahl der angegebenen Dateien: %i\n\n",argc-1);
	
	for (q=1;q<argc;q++) {
		
		/*********************
		 * MEMORY ALLOCATION *
		 *********************/
		
		//main data kubus
		int spalten = 15;
		double **kubus;
		kubus = (double **)calloc(kubussizetot,sizeof(double *));
		for (j=0;j<kubussizetot;j++){
			kubus[j] = (double *)calloc(spalten,sizeof(double));
			if (kubus[j] == NULL) {
				printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
				return 0;
			}
		}
		
		//neighbour array
		int nachbarn = standardnachbarzahl;
		int **nachbarliste;
		nachbarliste = (int **)calloc(kubussizetot,sizeof(int *));
		for (j=0;j<kubussizetot;j++){
			nachbarliste[j] = (int *)calloc(nachbarn,sizeof(int));
			if (nachbarliste[j] == NULL) {
				printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
				return 0;
			}
		}
		
		//surface-density arrays
		int anzahlradialbins = 100;
		double **profile;
		profile = (double **)calloc(anzahlradialbins,sizeof(double *));
		for (l=0;l<anzahlradialbins;l++){
			profile[l] = (double *)calloc(8,sizeof(double));
			if (profile[l] == NULL) {
				printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
				return 0;
			}
		}
		double **profile_inside_rtide;
		profile_inside_rtide = (double **)calloc(anzahlradialbins,sizeof(double *));
		for (l=0;l<anzahlradialbins;l++){
			profile_inside_rtide[l] = (double *)calloc(8,sizeof(double));
			if (profile_inside_rtide[l] == NULL) {
				printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
				return 0;
			}
		}
		double **profile_bound_stars;
		profile_bound_stars = (double **)calloc(anzahlradialbins,sizeof(double *));
		for (l=0;l<anzahlradialbins;l++){
			profile_bound_stars[l] = (double *)calloc(8,sizeof(double));
			if (profile_bound_stars[l] == NULL) {
				printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
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
				printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
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
		char dateiname[200];
		int poszaehler=1;
		if (code) {
			j = 2;
			while (j) {
				sprintf(dateiname, "%s.POS%i", argv[q],j);
				fz = fopen(dateiname,"rb");
				if (fz==NULL) {
					j = 0;
				} else {
					poszaehler = j;
					j++;
				}
			}
			
			sprintf(dateiname, "%s.POS", argv[q]);
			fz = fopen(dateiname,"rb");
			if (fz==NULL) {
				printf("\nDatei %s nicht gefunden!\n", dateiname);
				return 0;
			} 
		} else {
			sprintf(dateiname, "%s.txt", argv[q]);
			fz = fopen(dateiname,"r");
			if (fz==NULL) {
				printf("\nDatei %s nicht gefunden!\n", dateiname);
				return 0;
			} 
			j = 0;
		}
		int poszaehlertemp = 1;
		
		
		
		/*********************
		 * INITIALIZE OUTPUT *
		 *********************/
		
		char *ausgabedateiname;
		ausgabedateiname = strtok(dateiname, ".");            
		mkdir(ausgabedateiname,S_IRWXU|S_IRWXG|S_IRWXO );
		
		FILE *new;
		char zusammenfassung[200];
		sprintf(zusammenfassung, "%s.new",ausgabedateiname);
		new = fopen(zusammenfassung,"w");
		
		FILE *gal;
		char orbitdatei[200];
		sprintf(orbitdatei, "%s.gal",ausgabedateiname);
		gal = fopen(orbitdatei,"w");
		
		FILE *bin;	
		char binarydatei[200];
		sprintf(binarydatei, "%s.bin",ausgabedateiname);
		bin = fopen(binarydatei,"w");
		
		FILE *esc;	
		if (trackescapers) {
			char escaperdatei[200];
			sprintf(escaperdatei, "%s.esc",ausgabedateiname);
			esc = fopen(escaperdatei,"w");
		}
		
		
		
		/*****************************
		 * READ-IN AND ANALYSIS LOOP *
		 *****************************/
		
		double Rh = 1.0;
		int step = 0;
		
		do {
			printf("\n##########################################################\n");
			printf("Beginne Auswertung mit Datei %s (%i von %i)\n", dateiname, poszaehlertemp, poszaehler);
			
			long lSize, lPos;
			fseek (fz , 0 , SEEK_END);      // obtain file size
			lSize = ftell (fz);
			rewind (fz);
			if (lSize) {
				printf("Groesse der Datei: %li Bytes\n",lSize);
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
					
					for (j=0;j<kubussizetot;j++) kubus[j][0] = 0.0;
					
					if (code == 1) {
						readin1(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, kubus, ausgabetyp, kubussizetot);
						RG[0]=RG[1]=RG[2]=VG[0]=VG[1]=VG[2]=0;
					} else if (code == 2) {
						readin2(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, kubus, RG, VG, ausgabetyp, kubussizetot);
					} else if (code == 3) {
						readin3(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, kubus, RG, VG, ausgabetyp, kubussizetot);
					} else if (code == 4) {
						readin4(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, kubus, RG, VG, ausgabetyp, kubussizetot);
					} else if (code == 5) {
						readin5(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, kubus, RG, VG, ausgabetyp, kubussizetot);
					} else if (code == 6) {
						readin6(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, kubus, RG, VG, ausgabetyp, kubussizetot);
					} else {
						readin0(fz, &Ntot, &N, &time, &tbar, &rtide, &rbar, &vbar, &mbar, &omega, &Etot, &rcore, kubus, RG, VG, ausgabetyp, kubussizetot);
					}
					
					int timetemp;
					timetemp = floor(time+0.1);
					
					
					
					/************
					 * ANALYSIS *
					 ************/
					
					if (timetemp % analysistimestep < 1) {		
						
						if (ausgabetyp) printf("\nT: %g Myrs\n",time);
						
						if (step == 0) Ninit = N;
						Ntot = Ninit;
						
						
						//make neighbour lists, determine binding energies in case of schnellanalyse = 0
						if (standardnachbarzahl < N-2)
							nachbarn = standardnachbarzahl;
						else
							nachbarn = N-2;
						
						if (schnellanalyse) 
							bindungsanalysefast(Ntot, kubus, nachbarliste, nachbarn, mx, my, mz);
						else
							bindungsanalyse(Ntot, kubus, nachbarliste, nachbarn);
						
						
						//determine cluster centre
						if (standarddensitylevel < N-2) 
							densitylevel = standarddensitylevel;
						else
							densitylevel = N-2;
						mxold = mx;
						myold = my;
						mzold = mz;
						zentrumsuche(Ntot, kubus, nachbarliste, &mx, &my, &mz, densitylevel);
						if (ausgabetyp) printf("mx = %g\tmy = %g\tmz = %g\n",mx,my,mz);
						
						
						//determine omega
						double L;
						if (potenzial >= 0) {
							omegasuche(RG, VG, &omega, &L);
							//omega = 220.0/sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]); //approximation for circular orbits
							
							
							//determine second derivative of potential with respect to R (dphi) [for Kuepper et al 2009, eq. 11]
							if (potenzial > 1) {
								dphiberechner(RG, VG, omega, &dphi);
							} else {
								dphi = 2.0*omega*omega;
							}
							
							
							//determine rtide
							rtidesuche(&rtide, Ntot, kubus, omega, mx, my, mz, dphi, RG);
							if (rtide == 0.0) {
								rtide = pow(G*Ntot*0.4/sqrt(pow((dphi+omega*omega),2)),1.0/3.0);
							}
							
						} else {
							omega = 1.0E-09;
							rtide = 1.0E09;
						}
						
						if (ausgabetyp) printf("omega = %g\t", omega);
						if (ausgabetyp) printf("rtide = %g\n", rtide);
						
						
						//set status: 0 = normal; 1 = potential escaper; 2 = beyond 1*rtide; 3 = beyond 2*rtide; -nrbin = binary
						//determine binding energy in case of schnellanalyse = 1
						int nuncor[4]; int nbinuncor[4];
						for (j=0;j<kubussizetot;j++) kubus[j][13] = 0.0;
						statussetzer(Ntot, kubus, rtide, mx, my, mz, nuncor, nbinuncor,mvx,mvy,mvz,Rh,omega,schnellanalyse);
						
						
						//Sigma und mv berechnen, unkorrigiert
						status = standardstatus;
						geschwindigkeitsdispersionsberechnungsfunktion(Ntot,kubus,status,&sigmauncor,&mvxuncor,&mvyuncor,&mvzuncor);
						if (ausgabetyp) printf("mvx = %g\tmvy = %g\tmvz = %g\t --> sigma = %g (unkorrigiert)\n",mvxuncor,mvyuncor,mvzuncor,sigmauncor);
						
						
						//determine energies
						double Euncor[3];
						status = standardstatus;
						energieberechner(Ntot, kubus, status, Euncor,mvxuncor,mvyuncor,mvzuncor);
						if (ausgabetyp) printf("Ebin = %g\tEkin = %g\tEges = %g (unkorrigiert)\n", Euncor[0], Euncor[1], Euncor[2]);
						
						
						//detect and mark new escapers
						for(l=0;l<Ntot;l++) {
							if ((kubus[l][0]) && (!kubus[l][14]) && (kubus[l][13] > 2.0)) {
								kubus[l][14] = 1.0;
							}
						}
						
						
						//find binaries and calculate f_bin
						//build centre-of-mass particle and set binary-member status to -binnr
						int multiplicityarray[nachbarn];
						if (code && binarysearch) doppelsternpruefer(Ntot, kubus, nachbarliste, nachbarn, &Ntot, multiplicityarray, binoutput,bin,time,ausgabetyp);
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
						//determine binding energy in case of schnellanalyse = 1
						int n[4]; int nbin[4];
						for (j=0;j<kubussizetot;j++) if (kubus[j][13] >= 0.0) kubus[j][13] = 0.0;
						statussetzer(kubussizetot, kubus, rtide, mx, my, mz, n, nbin, mvx, mvy, mvz, Rh,omega,schnellanalyse);
						if (ausgabetyp) printf("Status (singles) \t0: %i\t1: %i\t2: %i\t3: %i\n", n[0],n[1], n[2], n[3]);
						if (ausgabetyp) printf("Status (binaries)\t0: %i\t2: %i\t3: %i\n", nbin[0], nbin[2], nbin[3]);
						double fbin;
						int ntemp = 0;
						for (j=0; j<status;j++) ntemp += n[j];
						int nbintemp = 0;
						for (j=0; j<status;j++) nbintemp += nbin[j];
						fbin = 1.0*nbintemp/ntemp;
						if (ausgabetyp) printf("fbin = %1.5lf\t%i\n", fbin, nbintemp);
						
						
						//determine corrected sigma and mv
						status = standardstatus;
						geschwindigkeitsdispersionsberechnungsfunktion(kubussizetot,kubus,status,&sigma,&mvx,&mvy,&mvz);
						if (ausgabetyp) printf("mvx = %g\tmvy = %g\tmvz = %g\t --> sigma = %g (corrected)\n",mvx,mvy,mvz,sigma);
						
						double mv = 0.0;
						int mvn = 0;
						
						for(l=0;l<Ntot;l++) {
							if ((kubus[l][0]) && (kubus[l][13]>=0.0) && (kubus[l][13]<2.0)) {
								mv += sqrt(pow(kubus[l][4]-mvx,2)+pow(kubus[l][5]-mvy,2)+pow(kubus[l][6]-mvz,2));
								mvn++;
							}
						}
						if (mvn) mv /= mvn;
						
						
						//determine corrected energies
						double E[3], kT;
						energieberechner(kubussizetot, kubus, status, E, mvx, mvy, mvz);
						kT = E[1]/ntemp;
						
						if (ausgabetyp) printf("Ebin = %g\tEkin = %g\tEges = %g (corrected)\n", E[0], E[1], E[2]);
						
						
						//determine radii, LR2, LR5, LR10, LR20, LR50/R_h, LR90, R_core
						double LR[6];
						double Rc = 0.0;
						double M;
						status = standardstatus;
						radien(Ntot, kubus, status, LR, mx, my, mz, &M);
						Rh = LR[4];
												
						double KR[11];
						double racc = 0.0; //radius at which a < a0
						if (step == 0) Minit = M;
						
						Kradien(Ninit, kubus, status, KR, mx, my, mz, &Minit, &racc);
						
						rcsuche(Ntot, kubus, nachbarliste, mx, my, mz, densitylevel, &Rc);
						
						if (ausgabetyp) printf("M = %g\tLR2 = %g\tLR5 = %g\tLR10 = %g\tLR20 = %g\tRh = %g\tLR90 = %g\tRc = %g\n", M, LR[0], LR[1], LR[2], LR[3], LR[4], LR[5], Rc);
						
						
						//determine surface-density profile					
						if (sdp) {
							//double rmax = 50.0;//pc
							//double rmin = 0.025;//pc
							double rmax = 500.0;//pc
							double rmin = 0.1;//pc
							double stepsize;
							stepsize = (log10(rmax)-log10(rmin))/(anzahlradialbins-1);
							
							for (l=0;l<anzahlradialbins;l++) {
								profile[l][0] = pow(10.0, log10(rmin) + stepsize*l); //radius
								//profile[l][0] = 1.0*exp(log(rmax-rmin)*(l+1)/anzahlradialbins)-1.0*exp(log(rmax-rmin)/anzahlradialbins)+rmin; //radius
								//profile[l][0] = 1.0*pow(base,log10(rmax-rmin)/log10(base)*(l+1)/anzahlradialbins)-1.0*pow(base,log10(rmax-rmin)/log10(base)/anzahlradialbins)+rmin; //radius
								//profile[l][0] = 1.0*pow(10.0,0.25*l-2.0); //radius
								profile[l][2] = 0; //mass
								profile[l][3] = 0; //mv
								profile[l][4] = 0; //m(v^2)
								profile[l][5] = 0; //number
								profile[l][6] = 0; //error of velocity dispersion
								profile[l][7] = 0; //luminosity
							}
							for (l=0;l<anzahlradialbins;l++) {
								if (l == 0)
									profile[l][1] =  2.0*PI*pow(profile[l][0],2);
								else
									profile[l][1] = 2.0*PI*(pow(profile[l][0],2) - pow(profile[l-1][0],2)); //area
							}
							
							for (l=0;l<anzahlradialbins;l++) {
								profile_inside_rtide[l][0] = pow(10.0, log10(rmin) + stepsize*l); //radius
								profile_inside_rtide[l][2] = 0; //mass
								profile_inside_rtide[l][3] = 0; //mv
								profile_inside_rtide[l][4] = 0; //m(v^2)
								profile_inside_rtide[l][5] = 0; //number
								profile_inside_rtide[l][6] = 0; //error of velocity dispersion
								profile_inside_rtide[l][7] = 0; //luminosity
							}
							for (l=0;l<anzahlradialbins;l++) {
								if (l == 0)
									profile_inside_rtide[l][1] =  2.0*PI*pow(profile_inside_rtide[l][0],2);
								else
									profile_inside_rtide[l][1] = 2.0*PI*(pow(profile_inside_rtide[l][0],2) - pow(profile_inside_rtide[l-1][0],2)); //area
							}
							
							for (l=0;l<anzahlradialbins;l++) {
								profile_bound_stars[l][0] = pow(10.0, log10(rmin) + stepsize*l); //radius
								//profile_bound_stars[l][0] = 1.0*pow(10.0,log10(rmax)*(l+1)/anzahlradialbins); //radius
								profile_bound_stars[l][2] = 0; //mass
								profile_bound_stars[l][3] = 0; //mv
								profile_bound_stars[l][4] = 0; //m(v^2)
								profile_bound_stars[l][5] = 0; //number
								profile_bound_stars[l][6] = 0; //error of velocity dispersion
								profile_bound_stars[l][7] = 0; //luminosity
							}
							for (l=0;l<anzahlradialbins;l++) {
								if (l == 0)
									profile_bound_stars[l][1] =  2.0*PI*pow(profile_bound_stars[l][0],2);
								else
									profile_bound_stars[l][1] = 2.0*PI*(pow(profile_bound_stars[l][0],2) - pow(profile_bound_stars[l-1][0],2)); //area
							}
							
							surfacedensityprofile(Ninit, kubus, profile, mx, my, mz,mvx,mvy,mvz, view, anzahlradialbins, VG, profile_inside_rtide, profile_bound_stars);
						}					
						
						
						//determine escape velocities
						if (trackescapers)
							escape(esc, kubussizetot, kubus, rtide, mx, my, mz, mvx, mvy, mvz, omega, time, sigma, kT, M, ausgabetyp, timebin, timebinsize, threshold);
						
						
						//stop analysis?
						if (M<Mende) lSize = 0;
						
						
						//output to STDOUT
						if (ausgabetyp == 0) {
							printf("\n---------------------------------------------------------------------------\n");
							printf("T = %g\tM = %g\tkT = %g\tsigma = %g\tfbin = %g\n",time,M,kT,sigma,fbin);
							printf("Rh = %g\tRtide = %g (d2phirdr2 = %g)\tRc = %g\n",Rh,rtide,-1.0*dphi,Rc);
							printf("Ebin = %g\tEkin = %g\tEtot = %g\tL = %g\n",E[0],E[1],E[2],L);
						}
						
						
						
						/************
						 *  OUTPUT  *
						 ************/
						
						FILE *aus;
						char ausgabe[200];
						double xtemp[3], vtemp[3], dsuntemp, vrsuntemp, vrtemp, lcosbtemp, RAtemp, DECtemp, mu_alphatemp, mu_alphacosdeltatemp, mu_deltatemp, mutemp, PAtemp;
						
						
						//write summary to .new
						fprintf(new,"%g\t%g\t%i\t%i\t%i\t%i\t%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",time,M,N,n[0],n[1],n[2],n[3],nbin[0],fbin,sigma,sigmauncor,E[0],E[1],E[2],Euncor[0],Euncor[1],Euncor[2],kT,LR[0],LR[1],LR[2],LR[3],LR[4],LR[5],Rc,rtide,omega,dphi,L, mv,KR[0],KR[1],KR[2],KR[3],KR[4],KR[5],KR[6],KR[7],KR[8],KR[9],racc);
						
						
						//write surface-density profile to .sdp
						if (sdp) {
							sprintf(ausgabe, "%s/%s%07i.sdp",ausgabedateiname,ausgabedateiname,timetemp);
							aus = fopen(ausgabe,"w");
							for(l=0;l<anzahlradialbins;l++) {
								fprintf(aus,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",profile[l][0],1.0*profile[l][5]/profile[l][1],profile[l][0]/rtide,profile[l][3],profile[l][4],profile[l][5],profile[l][6],profile_inside_rtide[l][0],1.0*profile_inside_rtide[l][5]/profile_inside_rtide[l][1],profile_inside_rtide[l][0]/rtide,profile_inside_rtide[l][3],profile_inside_rtide[l][4],profile_inside_rtide[l][5],profile_inside_rtide[l][6],profile_bound_stars[l][0],1.0*profile_bound_stars[l][5]/profile_bound_stars[l][1],profile_bound_stars[l][0]/rtide,profile_bound_stars[l][3],profile_bound_stars[l][4],profile_bound_stars[l][5],profile_bound_stars[l][6], 1.0*profile[l][7]/profile[l][1], 1.0*profile_inside_rtide[l][7]/profile[l][1], 1.0*profile_bound_stars[l][7]/profile[l][1], 1.0*profile[l][2]/profile[l][1], 1.0*profile_inside_rtide[l][2]/profile[l][1], 1.0*profile_bound_stars[l][2]/profile[l][1]);
							}
							fclose(aus);
							printf("Geschrieben nach: %s\n\n",ausgabe);
						}
						
						
						//write star dumps to all.txt and all.bin
						if ((allausgabe) && (timetemp % ausgabeschritt < 1)) {
							sprintf(ausgabe, "%s/%s%07iall.txt",ausgabedateiname,ausgabedateiname,timetemp);
							aus = fopen(ausgabe,"w");
							/*						fprintf(aus,"# rtide = %g\n",rtide);
							 fprintf(aus,"# omega = %g\n",omega);
							 fprintf(aus,"# M = %g\n",M);
							 fprintf(aus,"# time = %g\n",time);
							 fprintf(aus,"# mx = %g\n",mx);
							 fprintf(aus,"# my = %g\n",my);
							 fprintf(aus,"# mz = %g\n",mz);
							 fprintf(aus,"# RGx = %g\n",RG[0]);
							 fprintf(aus,"# RGy = %g\n",RG[1]);
							 fprintf(aus,"# RGz = %g\n",RG[2]);*/
							
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
									printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
									return 0;
								}
							}
							
							//Rausschreiben aller Sterne nach all.txt
							double lummean = 0.0;
							int Nlum = 0;
							for(l=0;l<Ntot;l++) {
								if ((kubus[l][0]) && (kubus[l][13]>=0.0)) {
									//fprintf(aus,"%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f\n",kubus[l][0],kubus[l][1],kubus[l][2],kubus[l][3],kubus[l][4],kubus[l][5],kubus[l][6],kubus[l][7],kubus[l][8],kubus[l][9],kubus[l][10],kubus[l][11],kubus[l][12],kubus[l][13]);
									fprintf(aus,"%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",kubus[l][0],kubus[l][1]-mx,kubus[l][2]-my,kubus[l][3]-mz,kubus[l][4]-mvx,kubus[l][5]-mvy,kubus[l][6]-mvz,kubus[l][7],kubus[l][8],kubus[l][9],kubus[l][10],kubus[l][11],kubus[l][12],kubus[l][13]);
									if (kubus[l][7]<10) {
										lummean += kubus[l][8];
										Nlum++;
									}
									//Binning in der x-y-Ebene
									xt = 1.0*(kubus[l][1]+0.5*xrange)/ybinsize;
									yt = 1.0*(kubus[l][2]+0.5*yrange)/ybinsize;//kubus[l][3] fuer perpendicular-to-disk models
									if ((yt<ybins) && (yt>=0) && ((xt<xbins) && (xt>=0))) {
										cont[xt][yt]++;
									}
								}
							}
							lummean /= Nlum;
							//printf("\nL_mean = %f\n",lummean);
							fclose(aus);
							printf("\nGeschrieben nach: %s\n",ausgabe);
							
							
							//Rausschreiben der bins nach ...all.bin
							sprintf(ausgabe, "%s/%s%07iall.bin",ausgabedateiname,ausgabedateiname,timetemp);
							aus = fopen(ausgabe,"w");
							int ll;
							for(l=0;l<xbins;l++) {
								for(ll=0;ll<ybins;ll++) {
									fprintf(aus,"%f\t%f\t%f\n",(l+0.5)*ybinsize-0.5*xrange,(ll+0.5)*ybinsize-0.5*yrange,cont[l][ll]);
								}
								fprintf(aus,"\n");
							}
							
							for (i=0;i<xbins;i++) free(cont[i]);
							free(cont);
							fclose(aus);
							printf("Geschrieben nach: %s\n\n",ausgabe);
						}	
						
						
						//write tails to y.txt and y.bin	
						if ((tailanalyse) && (timetemp % ausgabeschritt < 1)) {
							
							int pcorbitreichweite = 4000; // +/- pc
							int pcschrittweite = 25.0; //pc
							int pcanzahlintervalle;
							pcanzahlintervalle = 1*(pcorbitreichweite/pcschrittweite+1);
							double **entlangdesorbits;
							entlangdesorbits = (double **)calloc(pcanzahlintervalle,sizeof(double *));
							for (i=0;i<pcanzahlintervalle;i++){
								entlangdesorbits[i] = (double *)calloc(8,sizeof(double));
								if (entlangdesorbits[i] == NULL) {
									printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
									return 0;
								}
							}
							
							//Orbit vorwaerts berechnen
							int vorz = 1;
							for (i=0;i<pcanzahlintervalle;i++)
								entlangdesorbits[i][0] = vorz*i*pcschrittweite;
							
							double RGcorr[3], VGcorr[3];
							RGcorr[0] = RG[0]+mx;
							RGcorr[1] = RG[1]+my;
							RGcorr[2] = RG[2]+mz;
							VGcorr[0] = VG[0]+mvx;
							VGcorr[1] = VG[1]+mvy;
							VGcorr[2] = VG[2]+mvz;
							
							orbitausrechner(RGcorr,VGcorr,entlangdesorbits,pcanzahlintervalle,vorz,potenzial);
							
							double orbitbin[pcanzahlintervalle*2-3][11];
							
							for (i=0;i<pcanzahlintervalle-1;i++) {
								orbitbin[pcanzahlintervalle+i-2][0] = entlangdesorbits[i][0];
								orbitbin[pcanzahlintervalle+i-2][1] = entlangdesorbits[i][1];
								orbitbin[pcanzahlintervalle+i-2][2] = entlangdesorbits[i][2];
								orbitbin[pcanzahlintervalle+i-2][3] = entlangdesorbits[i][3];
								orbitbin[pcanzahlintervalle+i-2][4] = 0;
								orbitbin[pcanzahlintervalle+i-2][5] = entlangdesorbits[i][5];
								orbitbin[pcanzahlintervalle+i-2][6] = entlangdesorbits[i][6];
								orbitbin[pcanzahlintervalle+i-2][7] = entlangdesorbits[i][7];
								orbitbin[pcanzahlintervalle+i-2][8] = 0;
								orbitbin[pcanzahlintervalle+i-2][9] = 0;
								orbitbin[pcanzahlintervalle+i-2][10] = 0;
							}
							
							//Orbit rueckwaerts berechnen
							vorz = -1;
							for (i=0;i<pcanzahlintervalle;i++)
								entlangdesorbits[i][0] = vorz*i*pcschrittweite;
							
							orbitausrechner(RGcorr,VGcorr,entlangdesorbits,pcanzahlintervalle,vorz,potenzial);
							
							for (i=0;i<pcanzahlintervalle-1;i++) {
								orbitbin[i][0] = entlangdesorbits[pcanzahlintervalle-i-2][0];
								orbitbin[i][1] = entlangdesorbits[pcanzahlintervalle-i-2][1];
								orbitbin[i][2] = entlangdesorbits[pcanzahlintervalle-i-2][2];
								orbitbin[i][3] = entlangdesorbits[pcanzahlintervalle-i-2][3];
								orbitbin[i][4] = 0;
								orbitbin[i][5] = entlangdesorbits[pcanzahlintervalle-i-2][5];
								orbitbin[i][6] = entlangdesorbits[pcanzahlintervalle-i-2][6];
								orbitbin[i][7] = entlangdesorbits[pcanzahlintervalle-i-2][7];
								orbitbin[i][8] = 0;
								orbitbin[i][9] = 0;
								orbitbin[i][10] = 0;
							}
							
							
							//Naechsten Orbitpunkt fuer jeden Stern suchen und binnen
							double rtemp;
							double rmin[2];
							for(l=0;l<Ntot;l++) {
								if ((kubus[l][0]) && (kubus[l][13]>=0.0) && (kubus[l][10]>=limit)) {
									rmin[0] = 750.0;
									rmin[1] = 0;
									for(i=0;i<pcanzahlintervalle*2-3;i++) {
										rtemp = sqrt(pow(kubus[l][1]+RG[0]-orbitbin[i][1],2)+pow(kubus[l][2]+RG[1]-orbitbin[i][2],2)+pow(kubus[l][3]+RG[2]-orbitbin[i][3],2));
										if (rtemp < rmin[0]) {
											rmin[0] = rtemp;
											rmin[1] = i;
										}	
									}
									orbitbin[(int) rmin[1]][4] += 1.0;
									
									orbitbin[(int) rmin[1]][8] += sqrt(pow(kubus[l][4]+VG[0],2)+pow(kubus[l][5]+VG[1],2)+pow(kubus[l][6]+VG[2],2));
									
									orbitbin[(int) rmin[1]][9] += pow(kubus[l][4]+VG[0],2)+pow(kubus[l][5]+VG[1],2)+pow(kubus[l][6]+VG[2],2);
									
								}
							}
							for (i=0;i<pcanzahlintervalle*2-3;i++) {
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
							
							
							/*						for (i=0;i<pcanzahlintervalle*2-3;i++)
							 printf("%f\t%f\t%f\t%f\t%f\t%f\n",orbitbin[i][0],orbitbin[i][1],orbitbin[i][2],orbitbin[i][3],sqrt(pow(orbitbin[i][1]-orbitbin[0][1],2)+pow(orbitbin[i][2]-orbitbin[0][2],2)+pow(orbitbin[i][3]-orbitbin[0][3],2)),orbitbin[i][4]);*/
							
							//Schreiben der Bins entlang der Tails nach y.bin
							sprintf(ausgabe, "%s/%s%07iy.bin",ausgabedateiname,ausgabedateiname,timetemp);
							aus = fopen(ausgabe,"w");
							for(i=0;i<pcanzahlintervalle*2-3;i++) {
								fprintf(aus,"%f\t%f\t%f\t%f\t%f\n",orbitbin[i][0],orbitbin[i][4],orbitbin[i][8],orbitbin[i][9],orbitbin[i][10]);
							}
							
							for (i=0;i<pcanzahlintervalle;i++) free(entlangdesorbits[i]);
							free(entlangdesorbits);
							fclose(aus);
							
							printf("Geschrieben nach: %s\n\n",ausgabe);
						}
						
						
						//write snapshot from the sun to sun.txt
						if (timetemp % sunviewtime < 1) {
							FILE *aus2, *aus3;
							char ausgabe2[200], ausgabe3[200];
							sprintf(ausgabe, "%s/%s%07isun.txt",ausgabedateiname,ausgabedateiname,timetemp);
							aus = fopen(ausgabe,"w");
							
							sprintf(ausgabe2, "%s/%s%07isun.bin",ausgabedateiname,ausgabedateiname,timetemp);
							aus2 = fopen(ausgabe2,"w");
							
							sprintf(ausgabe3, "%s/%s%07isunvr.bin",ausgabedateiname,ausgabedateiname,timetemp);
							aus3 = fopen(ausgabe3,"w");
							
							
							//Beobachterparameter fuer das Haufenzentrum bestimmen
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
							//Teleskop bauen
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
									printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
									return 0;
								}
							}
							
							double **contvr;
							contvr = (double **)calloc(xbins,sizeof(double *));
							for (i=0;i<xbins;i++){
								contvr[i] = (double *)calloc(vrbins,sizeof(double));
								if (contvr[i] == NULL) {
									printf("\nDu NULL; Speicherinitialisierung fehlgeschlagen\n");
									return 0;
								}
							}
							
							//Fuer Pal5 M4, M12 und M16 bestimmen
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
							
							//Rausschreiben aller Sterne in Beobachterkoordinaten nach sun.txt
							for(l=0;l<Ntot;l++) {
								if ((kubus[l][0]) && (kubus[l][13]>=0.0) && (kubus[l][10]>=limit)) {
									xtemp[0] = RG[0]+kubus[l][1];
									xtemp[1] = RG[1]+kubus[l][2];
									xtemp[2] = RG[2]+kubus[l][3];
									vtemp[0] = VG[0]+kubus[l][4];
									vtemp[1] = VG[1]+kubus[l][5];
									vtemp[2] = VG[2]+kubus[l][6];
									convert(xtemp, vtemp, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSR);

									if (ltemp>180.0) {
										ltemp -=360.0;
										lcosbtemp = ltemp*cos(btemp);
									}
									
									xgal = kubus[l][1]+RG[0]+rgalsun;
									ygal = kubus[l][2]+RG[1];
									zgal = kubus[l][3]+RG[2];
									vr = vrsuntemp;
									Rtemp = dsuntemp;
									
									
									//Fuer Pal5 M4, M12 und M16 bestimmen
									rtemppal5 = sqrt(pow(ltemp-lgal,2)+pow(btemp-bgal,2));
									if (rtemppal5 < 4.0*0.00265) M4pal5+=kubus[l][10];
									if (rtemppal5 < 12.0*0.00265) M12pal5+=kubus[l][10];
									if (rtemppal5 < 16.0*0.00265) M16pal5+=kubus[l][10];
									
									
									//den Abstand in deg oder arcmin berechnen und mit ausgeben lassen, das binning kann dann spaeter geschehen
									//drdeg = (sqrt(pow(ltemp-lgal,2) + pow(btemp-bgal,2)))*360.0/(2.0*PI);
									x2t = xgal;
									y2t = ygal;
									z2t = zgal;
									drdeg = acos((x1t*x2t+y1t*y2t+z1t*z2t)/(sqrt(x1t*x1t+y1t*y1t+z1t*z1t)*sqrt(x2t*x2t+y2t*y2t+z2t*z2t)));
									
									
									//CMD erstellen
									if (kubus[l][9] && kubus[l][8]>=0.0) {
										Teff = pow(kubus[l][8]/(4.0*PI*kubus[l][9]*kubus[l][9]*kb),0.25);
										if ((Teff>3000.0) && (Teff<55000.0)) {
											
											lTeff = log10(Teff);
											
											BV = bvc[0] + bvc[1]*lTeff + bvc[2]*pow(lTeff,2) + bvc[3]*pow(lTeff,3) + bvc[4]*pow(lTeff,4) + bvc[5]*pow(lTeff,5) + bvc[6]*pow(lTeff,6) + bvc[7]*pow(lTeff,7);
											
											BC = bcc[0] + bcc[1]*lTeff + bcc[2]*pow(lTeff,2) + bcc[3]*pow(lTeff,3) + bcc[4]*pow(lTeff,4) + bcc[5]*pow(lTeff,5) + bcc[6]*pow(lTeff,6) + bcc[7]*pow(lTeff,7);
											
											if (kubus[l][8]) abvmag = -2.5*log10(kubus[l][8])-BC+BCsun+abvmagsun;
											
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
										fprintf(aus,"%.0f\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.0lf\t%.16lf\n",kubus[l][0],Rtemp,ltemp,btemp,vmag,abvmag,BV,vr,kubus[l][8],dvmag,kubus[l][10],dBV,Teff,kubus[l][13],drdeg);
										
										//Binning in Beobachterkoordinaten
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
							
							fclose(aus);
							printf("\nGeschrieben nach: %s\n",ausgabe);
							
							//Rausschreiben aller Bins in Beobachterkoordinaten nach sun.bin
							int ll;
							for(l=xmint;l<xmaxt;l++) {
								for(ll=ymint;ll<ymaxt;ll++) {
									fprintf(aus2,"%f\t%f\t%lf\n",(l+0.5)*ybinsize-180.0,(ll+0.5)*ybinsize-180.0,cont[l][ll]);
								}
								fprintf(aus2,"\n");
							}
							
							fclose(aus2);
							for(l=xmint;l<xmaxt;l++) {
								for(ll=vrmint;ll<vrmaxt;ll++) {
									fprintf(aus3,"%f\t%f\t%lf\n",(l+0.5)*ybinsize-180.0,(ll+0.5)*vrbinsize-100.0,contvr[l][ll]);
								}
								fprintf(aus3,"\n");
							}
							
							fclose(aus3);
							for (l=0;l<xbins;l++) free(cont[l]);
							free(cont);
							for (l=0;l<xbins;l++) free(contvr[l]);
							free(contvr);
							
							printf("Geschrieben nach: %s and %s\n\n",ausgabe2, ausgabe3);
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
					//if ((poszaehlertemp < 2) && (lPos>613440000)) lPos = lSize+1; //Option fuer 1213new
					//if (time>=1600) lPos = lSize+1; //Option fuer 1219
					if (ausgabetyp) printf("%li bytes of %li\n",lPos,lSize);
				} while(lPos+20<lSize);
				
				
			} else {
				printf("Datei %s ist leer!\n",dateiname);
				printf("##########################################################\n");
			}
			fclose(fz);
			poszaehlertemp++;
			sprintf(dateiname, "%s.POS%i", argv[q],poszaehlertemp);
			fz = fopen(dateiname,"rb");
			ausgabedateiname = strtok(dateiname, ".");
		} while (poszaehlertemp <= poszaehler);
		
		
		//write dynamical temperature-luminosity data to .tld
		if (trackescapers) {
			FILE *tld;	
			char tlddatei[200];
			sprintf(tlddatei, "%s.tld",ausgabedateiname);
			tld = fopen(tlddatei,"w");
			
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
		
		for (l=0;l<anzahlradialbins;l++) free (profile[l]);
		free(profile);
		for (j=0;j<kubussizetot;j++) free (kubus[j]);
		free(kubus);
		for (j=0;j<kubussizetot;j++) free(nachbarliste[j]);
		free(nachbarliste);	
		
		if (trackescapers) fclose(esc);
		fclose(new);
		fclose(gal);
		if (binoutput) fclose(bin);
		
		printf("\nZusammenfassung geschrieben nach: %s\n\n",zusammenfassung);
		
	}
	
	t2 = clock();             
	printf("\nAuswertungszeit: %g sec\n",(double)(t2-t1)/CLOCKS_PER_SEC); //write elapsed time to STDOUT
	return 1;
}






/***************
 * SUBROUTINES *
 ***************/

int escape(FILE *esc, int Nt, double **kubust, double rtidet, double mxt, double myt, double mzt, double mvxt, double mvyt, double mvzt, double omegat, double timet, double sigmat, double kTt, double Mt, int ausgabetypt, double **timebint, int timebinsizet, double thresholdt){
	int l,p,i;
	double vx, vy, vz, v, x, y, z, r, dr;
	double dt, Fx, Fy, Fz;
	double xt, yt, zt, timett,drt;//,vt,rt,drt;
	double vr, cosa, ratio;
	float dtlist[200][2];
	
	for (l=0;l<Nt;l++) {		
		if ((kubust[l][0]) && (kubust[l][14] < 2.0) && (kubust[l][14])) {
			
			x = kubust[l][1]-mxt;
			y = kubust[l][2]-myt;
			z = kubust[l][3]-mzt;
			r = sqrt(x*x+y*y+z*z);
			dr = r - 2.0*rtidet;
			vx = kubust[l][4]-mvxt;
			vy = kubust[l][5]-mvyt;
			vz = kubust[l][6]-mvzt;
			v = sqrt(vx*vx+vy*vy+vz*vz);
			timett = timet;
			cosa = (x*vx+y*vy+z*vz)/(r*v);
			vr = v*cosa;
			
			if (ausgabetypt) printf("   ESCAPE: %.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",timet,kubust[l][0],r,dr,v,vr,kubust[l][10],sigmat,kTt);
			
			
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
			Fx = kubust[l][10]*(G*Mt*x/(r*r*r) + 2.0*omegat*vy + 3.0*omegat*omegat*x);
			Fy = kubust[l][10]*(G*Mt*y/(r*r*r) - 2.0*omegat*vx);
			Fz = kubust[l][10]*(G*Mt*z/(r*r*r) + omegat*omegat*z);
			
			
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
					if (ausgabetypt) printf("CORRECTION FAILED (%.0f)\t",kubust[l][0]);
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
					if (ausgabetypt) printf("DTLIST: dr = %f\tdt = %f\n",drt,dt);
				}
				
				vx += Fx*dt;
				vy += Fy*dt;
				vz += Fz*dt;
				
				r = sqrt(xt*xt+yt*yt+zt*zt);
				v = sqrt(vx*vx+vy*vy+vz*vz);
				
				cosa = (xt*vx+yt*vy+zt*vz)/(r*v);
				vr = v*cosa;
				
				timett = timet+dt;
				if (ausgabetypt) printf("CORRECTED: %.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",timett,kubust[l][0],r,drt,v,vr,kubust[l][10],sigmat,kTt);
			}
			
			fprintf(esc,"%.0f\t%.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",timett,kubust[l][0],r,drt,v,vr,kubust[l][10],sigmat,kTt);
			
			i = timett/timebinsizet;
			if ((i >= 0) && (i<20000/timebinsizet) ) {
				timebint[i][0] += 1.0;
				timebint[i][1] += 0.5*sigmat;
				timebint[i][4] += Mt;
				
				ratio = v*v/sigmat;
				
				if (ratio < thresholdt) 
					timebint[i][2] += 0.5*kubust[l][10]*v*v;
				else
					timebint[i][3] += 0.5*kubust[l][10]*v*v;
			}
			
			kubust[l][14]=2.0;
			
			
		}
	}	
	
	
	
	return 1;
}

int dphiberechner(double *RGt, double *VGt, double omegat, double *dphit){
	double x, y, z, vx, vy, vz;
	double M1, M2, M3, a2, a3, b1, b2;
	double x2, y2, z2, vx2, vy2, vz2, b12, b22, a32;
	double d2phidr2;
	double r2, R2, r, r3, ra3, z2b22;
	double t1,t2,t3,t4,t5,t6;
	double temp1, temp2, temp3, temp4, temp5, temp6;
	
	x = RGt[0];
	y = RGt[1];
	z = RGt[2];
	vx = VGt[0];
	vy = VGt[1];
	vz = VGt[2];
	
	M1 = M1allen;
	M2 = M2allen;
	M3 = M3allen;
	a2 = a2allen;
	b1 = b1allen;
	b2 = b2allen;
	a3 = a3allen;
	
	x2 = x*x;
	y2 = y*y;
	z2 = z*z;
	vx2 = vx*vx;
	vy2 = vy*vy;
	vz2 = vz*vz;
	b12 = b1*b1;
	b22 = b2*b2;
	a32 = a3*a3;
	
	r2 = x2+y2+z2;
	r = sqrt(r2);
	r3 = pow(r,3);
	R2 = x2 + y2;
	ra3 = 1.0+1.0*r/a3;
	z2b22 = sqrt(z2+b22);
	
	t1 = (-3.0*M1*x2/pow(r2+b12,2.5) +M1/pow(r2+b12,1.5) -3.0*M2*x2/pow(R2+pow(a2+z2b22,2),2.5) +M2/pow(R2+pow(a2+z2b22,2),1.5) +M3*x2/(r3*a32*ra3)  +M3*x2/(r2*pow(a3,3)*pow(ra3,2)) -M3/(r*a32*ra3) -2.0*M3*x2/(pow(a3,4)*pow(ra3,3)*r) +M3/(pow(a3,3)*pow(ra3,2)) -M3*(2.0*x2/(pow(ra3,3)*r2*a32)+x2/(pow(ra3,2)*r3*a3)-1.0/(pow(ra3,2)*r*a3)+x2/(r3*a3*ra3)-1.0/(r*a3*ra3)+x2/(r2*a32*pow(ra3,2)))/a3);
	
	t2 = (-3.0*M1*y2/pow(r2+b12,2.5) +M1/pow(r2+b12,1.5) -3.0*M2*y2/pow(R2+pow(a2+z2b22,2),2.5) +M2/pow(R2+pow(a2+z2b22,2),1.5) +M3*y2/(r3*a32*ra3) +M3*y2/(r2*pow(a3,3)*pow(ra3,2)) -M3/(r*a32*ra3) -2.0*M3*y2/(pow(a3,4)*pow(ra3,3)*r) +M3/(pow(a3,3)*pow(ra3,2)) -M3*(2.0*y2/(pow(ra3,3)*r2*a32)+y2/(pow(ra3,2)*r3*a3)-1.0/(pow(ra3,2)*r*a3)+y2/(r3*a3*ra3)-1.0/(r*a3*ra3)+y2/(r2*a32*pow(ra3,2)))/a3);
	
	t3 = (-3.0*M1*z2/pow(r2+b12,2.5) +M1/pow(r2+b12,1.5)-3.0*M2*pow(a2+z2b22,2)*z2/(pow(R2+pow(a2+z2b22,2),2.5)*(z2+b22)) +M2*z2/(pow(R2+pow(a2+z2b22,2),1.5)*(z2+b22)) -M2*(a2+z2b22)*z2/(pow(R2+pow(a2+z2b22,2),1.5)*pow(z2+b22,1.5)) +M2*(a2+z2b22)/(pow(R2+pow(a2+z2b22,2),1.5)*z2b22) +M3*z2/(r3*a32*ra3) +M3*z2/(r2*pow(a3,3)*pow(ra3,2)) -M3/(r*a32*ra3) -2.0*M3*z2/(pow(a3,4)*pow(ra3,3)*r) +M3/(pow(a3,3)*pow(ra3,2)) -M3*(2.0*z2/(pow(ra3,3)*r2*a32)+z2/(pow(ra3,2)*r3*a3)-1.0/(pow(ra3,2)*r*a3)+z2/(r3*a3*ra3)-1.0/(r*a3*ra3)+z2/(r2*a32*pow(ra3,2)))/a3);
	
	t4 = (-3.0*M1*y*x/pow(r2+b12,2.5) -3.0*M2*y*x/pow(R2+pow(a2+z2b22,2),2.5) +M3*y*x/(r3*a32*ra3) +M3*y*x/(r2*pow(a3,3)*pow(ra3,2)) -2.0*M3*y*x/(pow(a3,4)*pow(ra3,3)*r) -M3*(2.0*y*x/(pow(ra3,3)*r2*a32)+y*x/(pow(ra3,2)*r3*a3)+y*x/(r3*a3*ra3)+y*x/(r2*a32*pow(ra3,2)))/a3);
	
	t5 = (-3.0*M1*z*x/pow(r2+b12,2.5) -3.0*M2*(a2+z2b22)*z*x/(pow(R2+pow(a2+z2b22,2),2.5)*z2b22) +M3*z*x/(r3*a32*ra3)+M3*z*x/(r2*pow(a3,3)*pow(ra3,2)) -2.0*M3*z*x/(pow(a3,4)*pow(ra3,3)*r) -M3*(2.0*z*x/(pow(ra3,3)*r2*a32)+z*x/(pow(ra3,2)*r3*a3)+z*x/(r3*a3*ra3)+z*x/(r2*a32*pow(ra3,2)))/a3);
	
	t6 = (-3.0*M1*z*y/pow(r2+b12,2.5) -3.0*M2*(a2+z2b22)*z*y/(pow(R2+pow(a2+z2b22,2),2.5)*z2b22) +M3*z*y/(r3*a32*ra3) +M3*z*y/(r2*pow(a3,3)*pow(ra3,2)) -2.0*M3*z*y/(pow(a3,4)*pow(ra3,3)*r) -M3*(2.0*z*y/(pow(ra3,3)*r2*a32)+z*y/(pow(ra3,2)*r3*a3)+z*y/(r3*a3*ra3)+z*y/(r2*a32*pow(ra3,2)))/a3);
	//printf("%f\t%f\t%f\t%f\t%f\t%f\n", t1,t2,t3,t4,t5,t6);
	//printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x2,y2,z2,x,y,z,r2);
	
	temp1 = x2*t1;
	temp2 = y2*t2;
	temp3 = z2*t3;
	temp4 = 2.0*x*y*t4;
	temp5 = 2.0*x*z*t5;
	temp6 = 2.0*y*z*t6;
	//printf("%f\t%f\t%f\t%f\t%f\t%f\n", temp1,temp2,temp3,temp4,temp5,temp6);
	d2phidr2 = 1.0*G*(temp1 + temp2 + temp3 + temp4 + temp5 + temp6)/r2;
	
	*dphit = -d2phidr2;
	//printf("d2phidr2 = %f\tomega^2 = %f\n", d2phidr2, omegat*omegat);
	
	return 1;
	
} 

int rcsuche(int Nt, double **kubust, int **nachbarlistet, double mxt, double myt, double mzt, int densitylevelt, double *Rct){ 
	int l, ll;
	double mtemp, r2temp, rho2temp, rho2ges = 0.0, rho2malr2ges = 0.0, rho2malr2temp;
	
	for (l=0;l<Nt;l++) {
		if ((kubust[l][0]) && (kubust[l][13]<3)) {
			mtemp = 0.0;
			for (ll=0;ll<densitylevelt;ll++) {
				mtemp += kubust[nachbarlistet[l][ll]][10];
			}
			rho2temp = pow(mtemp/pow(pow(kubust[nachbarlistet[l][densitylevelt]][1]-kubust[l][1],2)+pow(kubust[nachbarlistet[l][densitylevelt]][2]-kubust[l][2],2)+pow(kubust[nachbarlistet[l][densitylevelt]][3]-kubust[l][3],2),1.5),2);
			r2temp = pow(mxt-kubust[l][1],2)+pow(myt-kubust[l][2],2)+pow(mzt-kubust[l][3],2);
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

int orbitausrechner(double *RGt,double *VGt,double **entlangdesorbitst,int pcanzahlintervallet, int vorzt, int potenzialt){
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
	double pcmax = entlangdesorbitst[pcanzahlintervallet-1][0];
	double dpcout = entlangdesorbitst[1][0];
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
	
	dt = vorzt*0.00001;                    /* Initial Condition fuer die Zeitschritte */
	
	int step = 0;
	
	while (vorzt * pc <= vorzt * pcmax) {
		if (vorzt * pc >= vorzt * pcout) {              /* Ausgabe der Werte */
			entlangdesorbitst[step][1] = x[0];
			entlangdesorbitst[step][2] = x[1];
			entlangdesorbitst[step][3] = x[2];
			entlangdesorbitst[step][4] = 0;
			entlangdesorbitst[step][5] = v[0];
			entlangdesorbitst[step][6] = v[1];
			entlangdesorbitst[step][7] = v[2];
			step++;
			pcout = entlangdesorbitst[step][0];              /* Distanzgrenze fuer Ausgabe wird erhoeht */
		}
		
		do {
			for (k=0;k<3;k++) {
				xe1[k]=x[k];
				xe2[k]=x[k];
				ve1[k]=v[k];
				ve2[k]=v[k];
			}
			do_step(dt,xe1,ve1,potenzialt);      /* Ein ganzer Zeitschritt */
			
			do_step(0.5*dt,xe2,ve2,potenzialt);  /* Zwei halbe Zeitschritte */
			do_step(0.5*dt,xe2,ve2,potenzialt);
			
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
			
			if ((vorzt*dt < 0.0000000001) ) {
				printf("Abbruch... dt = %lf\n", dt);
				pc = pcmax*2.0;
				diff = mdiff/2;
			}
			
			if (vorzt * dt > vorzt * dtout) {
				dt=dtout;
			}
			
		} while (diff>mdiff);          /* Schleife wird erst einmal ausgefuehrt und muss nur bei zu grosser Differenz wiederholt werden */
	}
	
	return 0;
}

void getforce(double *x, double *v, double *a, int potenzialt){ 
	double r1, r2, r3;
	double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;
	
	if (potenzialt == 2) {
		//Point mass
		r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2) + b1allen*b1allen);
		
		a1x = -G*M1allen/(r1*r1*r1)**x;
		a1y = -G*M1allen/(r1*r1*r1)**(x+1);
		a1z = -G*M1allen/(r1*r1*r1)**(x+2);
		
		//Miyamato disk
		r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2allen + sqrt(*(x+2) * *(x+2) + b2allen*b2allen),2));
		
		a2x = -G*M2allen/(r2*r2*r2) * *x;
		a2y = -G*M2allen/(r2*r2*r2) * *(x+1);
		a2z = -G*M2allen/(r2*r2*r2) * (a2allen + sqrt(*(x+2) * *(x+2) + b2allen*b2allen))/sqrt(*(x+2) * *(x+2) + b2allen*b2allen) * *(x+2);
		
		//Log Halo
		r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));
		
		a3x = -G*M3allen/(a3allen*a3allen +a3allen*r3) * *x/r3;
		a3y = -G*M3allen/(a3allen*a3allen +a3allen*r3) * *(x+1)/r3;
		a3z = -G*M3allen/(a3allen*a3allen +a3allen*r3) * *(x+2)/r3;
		
		*(a+0) = a1x + a2x + a3x;
		*(a+1) = a1y + a2y + a3y;
		*(a+2) = a1z + a2z + a3z;
		
	} else if (potenzialt == 1) {
		//Point mass
		r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));
		
		*(a+0) = -G*M1pointmass/(r1*r1*r1)**x;
		*(a+1) = -G*M1pointmass/(r1*r1*r1)**(x+1);
		*(a+2) = -G*M1pointmass/(r1*r1*r1)**(x+2);
	}
	
}

void do_step(double dt, double *x, double *v, int potenzialt) {
	double hh, *a0,*a1,*a2,*a3,xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	int k;
	
	a0 = malloc(3*sizeof(double));
	a1 = malloc(3*sizeof(double));
	a2 = malloc(3*sizeof(double));
	a3 = malloc(3*sizeof(double));
	
	hh = dt*0.5;
	
	getforce(x,v,a0,potenzialt);                    /* Zwischenschritte */
	for (k=0;k<3;k++) {                /* erster Halbschritt */
		xt1[k] = *(x+k)+hh**(v+k);
		vt1[k] = *(v+k)+hh**(a0+k);
	}  
	
	getforce(&xt1[0], &vt1[0], a1,potenzialt);
	for (k=0;k<3;k++) {                /* zweiter Halbschritt */
		xt2[k] = *(x+k)+hh*vt1[k];
		vt2[k] = *(v+k)+hh**(a1+k);
	}
	
	getforce(&xt2[0], &vt2[0], a2,potenzialt);
	for (k=0;k<3;k++) {                /* dritter Schritt mit Ergebnissen des zweiten */
		xt3[k] = *(x+k)+dt*vt2[k];
		vt3[k] = *(v+k)+dt**(a2+k);
	}
	
	getforce(&xt3[0], &vt3[0], a3,potenzialt);
	for (k=0;k<3;k++) {                /* Runge-Kutta-Formel */
		*(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
		*(v+k) += dt/6.0*(*(a0+k)+2.0*(*(a1+k)+*(a2+k))+*(a3+k));
	}
	free(a0); free(a1); free(a2); free(a3);   /* gibt Speicher wieder frei */
}

int radien(int Nt, double **kubust, int statust, double *LRt, double mxt, double myt, double mzt, double *Mt){
	int l;
	double **radienarray;	
	radienarray = (double **)calloc(Nt,sizeof(double *));
	for (l=0;l<Nt;l++){
		radienarray[l] = (double *)calloc(2,sizeof(double));
		if (radienarray[l] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}	
	
	
	//	float radienarray[Nt][2];
	double rtemp, Mgestemp = 0.0, Mtemp = 0.0;
	LRt[0] = 0.0; LRt[1] = 0.0; LRt[2] = 0.0; LRt[3] = 0.0; LRt[4] = 0.0; LRt[5] = 0.0;
	
	for (l=0;l<Nt;l++){
		if ((kubust[l][0]) && (kubust[l][13]<statust) && (kubust[l][13]>= 0.0)) {
			rtemp = sqrt(pow(kubust[l][1]-mxt,2) + pow(kubust[l][2]-myt,2) + pow(kubust[l][3]-mzt,2));
			radienarray[l][0] = rtemp;
			Mgestemp += kubust[l][10];
			radienarray[l][1] = kubust[l][10];
		} else {
			radienarray[l][0] = 1000000;
			radienarray[l][1] = 0;
		}
	}
	
	*Mt = Mgestemp;
	
	//qsort(radienarray, sizeof(radienarray)/sizeof(radienarray[0]), sizeof(radienarray[0]), (void *)cmpmy);
	shellsort_reverse(radienarray, Nt, 2);
	
	for (l=0;l<Nt;l++) {
		Mtemp += radienarray[l][1];
		if ((LRt[0] == 0) && (Mtemp >= 0.02*Mgestemp))
			LRt[0] = radienarray[l][0];
		else if ((LRt[1] == 0) && (Mtemp >= 0.05*Mgestemp))
			LRt[1] = radienarray[l][0];
		else if ((LRt[2] == 0) && (Mtemp >= 0.1*Mgestemp))
			LRt[2] = radienarray[l][0];
		else if ((LRt[3] == 0) && (Mtemp >= 0.2*Mgestemp))
			LRt[3] = radienarray[l][0];
		else if ((LRt[4] == 0) && (Mtemp >= 0.5*Mgestemp))
			LRt[4] = radienarray[l][0];
		else if ((LRt[5] == 0) && (Mtemp >= 0.9*Mgestemp))
			LRt[5] = radienarray[l][0];
	}
	
	for (l=0;l<Nt;l++) free (radienarray[l]);
	free(radienarray);
	
	return 1;
}

int Kradien(int Nt, double **kubust, int statust, double *LRt, double mxt, double myt, double mzt, double *Mt, double *racc){
	int l;
	double **radienarray;	
	radienarray = (double **)calloc(Nt,sizeof(double *));
	for (l=0;l<Nt;l++){
		radienarray[l] = (double *)calloc(2,sizeof(double));
		if (radienarray[l] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}	
	//	float radienarray[Nt][2];
	double rtemp, Mgestemp = 0.0, Mtemp = 0.0;
	LRt[0] = 0.0; LRt[1] = 0.0; LRt[2] = 0.0; LRt[3] = 0.0; LRt[4] = 0.0; LRt[5] = 0.0; LRt[5] = 0.0; LRt[6] = 0.0; LRt[7] = 0.0; LRt[8] = 0.0; LRt[9] = 0.0; LRt[10] = 0.0;
	
	for (l=0;l<Nt;l++){
		if ((kubust[l][0]) && (kubust[l][13]>=0.0)) {
			rtemp = sqrt(pow(kubust[l][1]-mxt,2) + pow(kubust[l][2]-myt,2) + pow(kubust[l][3]-mzt,2));
			radienarray[l][0] = rtemp;
			radienarray[l][1] = kubust[l][10];
		} else {
			radienarray[l][0] = 1000000.0;
			radienarray[l][1] = 0;
		}
	}
	
	Mgestemp = *Mt;
	
	//	qsort(radienarray, sizeof(radienarray)/sizeof(radienarray[0]), sizeof(radienarray[0]), (void *)cmpmy);
	shellsort_reverse(radienarray, Nt, 2);
	
	double acc; //local acceleration
	*racc = 0.0; //radius at which a < a0
	double a0 = 3.6;//pc/Myr^2   //1.2e-10; //a0 in m/s^2
	
	for (l=0;l<Nt;l++) {
		Mtemp += radienarray[l][1];
		if ((*racc==0.0) && (l>3)) {
			acc = 1.0*G*Mtemp/pow(radienarray[l][0],2);
			if (acc < a0) *racc = radienarray[l][0]; 
		}
		if ((LRt[0] == 0) && (Mtemp >= 0.02*Mgestemp))
			LRt[0] = radienarray[l][0];
		else if ((LRt[1] == 0) && (Mtemp >= 0.05*Mgestemp))
			LRt[1] = radienarray[l][0];
		else if ((LRt[2] == 0) && (Mtemp >= 0.1*Mgestemp))
			LRt[2] = radienarray[l][0];
		else if ((LRt[3] == 0) && (Mtemp >= 0.2*Mgestemp))
			LRt[3] = radienarray[l][0];
		else if ((LRt[4] == 0) && (Mtemp >= 0.3*Mgestemp))
			LRt[4] = radienarray[l][0];
		else if ((LRt[5] == 0) && (Mtemp >= 0.4*Mgestemp))
			LRt[5] = radienarray[l][0];
		else if ((LRt[6] == 0) && (Mtemp >= 0.5*Mgestemp))
			LRt[6] = radienarray[l][0];
		else if ((LRt[7] == 0) && (Mtemp >= 0.6*Mgestemp))
			LRt[7] = radienarray[l][0];
		else if ((LRt[8] == 0) && (Mtemp >= 0.7*Mgestemp))
			LRt[8] = radienarray[l][0];
		else if ((LRt[9] == 0) && (Mtemp >= 0.8*Mgestemp))
			LRt[9] = radienarray[l][0];
		else if ((LRt[5] == 0) && (Mtemp >= 0.9*Mgestemp))
			LRt[5] = radienarray[l][0];
	}
	
	for (l=0;l<Nt;l++) free (radienarray[l]);
	free(radienarray);
	
	return 1;
}

int surfacedensityprofile(int Nt, double **kubust, double **profilet, double mxt, double myt, double mzt, double mvxt, double mvyt, double mvzt, int viewt, int anzahlradialbinst, double *VGt, double **profile_inside_rtidet, double **profile_bound_starst){
	int l;
	double **radienarray;	
	radienarray = (double **)calloc(Nt,sizeof(double *));
	for (l=0;l<Nt;l++){
		radienarray[l] = (double *)calloc(5,sizeof(double));
		if (radienarray[l] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}	
	double rtemp;
	
	for (l=0;l<Nt;l++){
		if ((kubust[l][0]) && (kubust[l][13]>=0.0)) {
			if (viewt == 1) {
				rtemp = sqrt(pow(kubust[l][1]-mxt,2) + pow(kubust[l][2]-myt,2));//xy-Ebene
			} else if (viewt == 2) {
				rtemp = sqrt(pow(kubust[l][1]-mxt,2) + pow(kubust[l][3]-mzt,2));//xz-Ebene
			}
			radienarray[l][0] = rtemp;
			radienarray[l][1] = kubust[l][10];//Masse
			//			radienarray[l][2] = pow(kubust[l][4]+VGt[0],2)+pow(kubust[l][5]+VGt[1],2)+pow(kubust[l][6]+VGt[2],2);//v^2
			radienarray[l][2] = pow(kubust[l][6]+VGt[2],2);//v
			radienarray[l][3] = kubust[l][13];//status
			radienarray[l][4] = kubust[l][8];//luminosity
		} else {
			radienarray[l][0] = 1000000.0;
			radienarray[l][1] = 0;
			radienarray[l][2] = 0;
			radienarray[l][3] = -1;
			radienarray[l][4] = 0;//luminosity
		}
	}
	
	//qsort(radienarray, sizeof(radienarray)/sizeof(radienarray[0]), sizeof(radienarray[0]), (void *)cmpmy);
	shellsort_reverse(radienarray, Nt, 5);
	
	//Prozedur fuer profile, profile_inside_rtide und profile_bound_stars machen
	//profile (all stars)
	int ll = 0;
	l = 0;
	while ((ll<anzahlradialbinst) && (l<Nt)) {
		if ((radienarray[l][3]>=0) && (radienarray[l][0]<profilet[ll][0])) {
			profilet[ll][2]+=radienarray[l][1];
			profilet[ll][3]+=sqrt(radienarray[l][2]);
			profilet[ll][4]+=radienarray[l][2];
			profilet[ll][5]++;
			profilet[ll][7]+=radienarray[l][4];
			l++;
		} else {
			ll++;
		}
	}
	
	ll = 0;
	l = 0;
	while ((ll<anzahlradialbinst) && (l<Nt)) {
		if ((radienarray[l][3]>=0) && (radienarray[l][0]<profilet[ll][0])) {
			profilet[ll][6] += pow(pow(sqrt(radienarray[l][2])-profilet[ll][3]/profilet[ll][5],2) - (profilet[ll][4]/profilet[ll][5]-pow(profilet[ll][3]/profilet[ll][5],2)),2);
			l++;
		} else {
			ll++;
		}
	}
	
	
	//profile_inside_rtide (all stars inside rtide)
	ll = 0;
	l = 0;
	while ((ll<anzahlradialbinst) && (l<Nt)) {
		if ((radienarray[l][3]>=0) && (radienarray[l][3]<2) && (radienarray[l][0]<profile_inside_rtidet[ll][0])) {
			profile_inside_rtidet[ll][2]+=radienarray[l][1];
			profile_inside_rtidet[ll][3]+=sqrt(radienarray[l][2]);
			profile_inside_rtidet[ll][4]+=radienarray[l][2];
			profile_inside_rtidet[ll][5]++;
			profile_inside_rtidet[ll][7]+=radienarray[l][4];
			l++;
		} else if ((radienarray[l][3]>1) || (radienarray[l][3]<0)) {
			l++;
		} else {
			ll++;
		}
	}
	
	ll = 0;
	l = 0;
	while ((ll<anzahlradialbinst) && (l<Nt)) {
		if ((radienarray[l][3]>=0) && (radienarray[l][3]<2) && (radienarray[l][0]<profile_inside_rtidet[ll][0])) {
			profile_inside_rtidet[ll][6] += pow(pow(sqrt(radienarray[l][2])-profile_inside_rtidet[ll][3]/profile_inside_rtidet[ll][5],2) - (profile_inside_rtidet[ll][4]/profile_inside_rtidet[ll][5]-pow(profile_inside_rtidet[ll][3]/profile_inside_rtidet[ll][5],2)),2);
			l++;
		} else if ((radienarray[l][3]>1) || (radienarray[l][3]<0)) {
			l++;
		} else {
			ll++;
		}
	}
	
	
	//profile_bound_stars (all bound stars)
	ll = 0;
	l = 0;
	while ((ll<anzahlradialbinst) && (l<Nt)) {
		if ((radienarray[l][3] == 0.0) && (radienarray[l][0]<profile_bound_starst[ll][0])) {
			profile_bound_starst[ll][2]+=radienarray[l][1];
			profile_bound_starst[ll][3]+=sqrt(radienarray[l][2]);
			profile_bound_starst[ll][4]+=radienarray[l][2];
			profile_bound_starst[ll][5]++;
			profile_bound_starst[ll][7]+=radienarray[l][4];
			l++;
		} else if ((radienarray[l][3]>=1.0) || (radienarray[l][3]<=-1.0)) {
			l++;
		} else {
			ll++;
		}
	}
	
	ll = 0;
	l = 0;
	while ((ll<anzahlradialbinst) && (l<Nt)) {
		if ((radienarray[l][3]<1.0) && (radienarray[l][3]>-1.0) && (radienarray[l][0]<profile_bound_starst[ll][0])) {
			profile_bound_starst[ll][6] += pow(pow(sqrt(radienarray[l][2])-profile_bound_starst[ll][3]/profile_bound_starst[ll][5],2) - (profile_bound_starst[ll][4]/profile_bound_starst[ll][5]-pow(profile_bound_starst[ll][3]/profile_bound_starst[ll][5],2)),2);
			l++;
		} else if ((radienarray[l][3]>=1.0) || (radienarray[l][3]<-1.0)) {
			l++;
		} else {
			ll++;
		}
	}
	
	for (l=0;l<Nt;l++) free (radienarray[l]);
	free(radienarray);
	
	return 1;
	
}

int geschwindigkeitsdispersionsberechnungsfunktion(int Nt,double **kubust,int statust, double *sigmat, double *mvxt, double *mvyt, double *mvzt){
	int l, n = 0;
	double sigmax = 0.0, sigmay  = 0.0, sigmaz = 0.0;
	*sigmat = 0.0;
	*mvxt = 0.0; *mvyt = 0.0; *mvzt = 0.0;
	double mv = 0.0, mv2 = 0.0, mtemp = 0.0;
	double vtemp = 0.0, sigmatemp = 0.0;
	double stdsigma = 5000.0;//optional velocity cut, disabled by default
	int nsigma = 100.0;//optional velocity cut, disabled by default 
	
	for (l=0;l<Nt;l++) {
		vtemp = sqrt(kubust[l][4]*kubust[l][4]+kubust[l][5]*kubust[l][5]+kubust[l][6]*kubust[l][6]);
		if (n>10) sigmatemp = sqrt((mv2-1.0/(1.0*n)*mv*mv)/(1.0*n-1.0));
		else sigmatemp = stdsigma; 
		
		if ((kubust[l][0]) && (kubust[l][13]<statust) && (kubust[l][13] >= 0) && (vtemp < nsigma*sigmatemp)) {
			*mvxt += kubust[l][10]*kubust[l][4];
			*mvyt += kubust[l][10]*kubust[l][5];
			*mvzt += kubust[l][10]*kubust[l][6];
			mtemp += kubust[l][10];
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
		vtemp = sqrt(pow(kubust[l][4]-*mvxt,2)+pow(kubust[l][5]-*mvyt,2)+pow(kubust[l][6]-*mvzt,2));
		//printf("%f\t%f\n",vtemp,sigmatemp);
		if ((kubust[l][0]) && (kubust[l][13]<statust) && (kubust[l][13] >= 0) && (vtemp < nsigma*sigmatemp)) {
			sigmax += pow(kubust[l][4]-*mvxt,2);
			sigmay += pow(kubust[l][5]-*mvyt,2);
			sigmay += pow(kubust[l][6]-*mvzt,2);
			n++;
		}
	}
	if (n) {
		sigmax /= n;
		sigmay /= n;
		sigmaz /= n;
	}
	*sigmat = sigmax + sigmay + sigmaz;
	
	return 1;
}

int doppelsternpruefer(int Nt, double **kubust, int **nachbarlistet, int nachbarnt, int *Ntott, int *multiplicityarrayt, int binoutputt, FILE *bint, double timet, int ausgabetypt){
	int l, ll, namell;
	double Ebintemp, atemp, etemp, rtemp,rx,ry,rz,vx,vy,vz;
	
	int name1, name2, multiplicity;
	for (l=0;l<nachbarnt;l++) multiplicityarrayt[l] = 0;
	
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
		if (kubust[l][0]) {
			for (ll=0;ll<1;ll++) {
				//for (ll=0;ll<nachbarnt;ll++) {
				namell = nachbarlistet[l][ll];
				rx = kubust[l][1]-kubust[namell][1];
				ry = kubust[l][2]-kubust[namell][2];
				rz = kubust[l][3]-kubust[namell][3];
				rtemp = sqrt(rx*rx+ry*ry+rz*rz);
				
				vx = kubust[l][4]-kubust[namell][4];
				vy = kubust[l][5]-kubust[namell][5];
				vz = kubust[l][6]-kubust[namell][6];
				
				Ebintemp = 0.5*kubust[l][10]*kubust[namell][10]/(kubust[l][10]+kubust[namell][10])*(vx*vx+vy*vy+vz*vz) - G*kubust[l][10]*kubust[namell][10]/rtemp;
				
				atemp = -0.5*G*kubust[l][10]*kubust[namell][10]/Ebintemp;
				
				etemp = sqrt(pow(1.0-rtemp/atemp,2) - pow(rx*vx+ry*vy+rz*vz,2)/(G*(kubust[l][10]+kubust[namell][10])*atemp));
				
				if (etemp<1.0) {
					if (binoutputt) printf("Binary:\t1 = %i\t2 = %i (%ith)\tr12 = %g\tEges = %g\ta = %g\tecc = %g\t1: %g\t2: %g\n", l, namell, ll, rtemp, Ebintemp, atemp, etemp, kubust[l][13],kubust[namell][13]);
					if (binoutputt) fprintf(bint,"%.0f\t%i\t%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",timet , l, namell, ll, rtemp, Ebintemp, atemp, etemp, kubust[l][10],kubust[namell][10],sqrt(pow(kubust[l][1],2)+pow(kubust[l][2],2)+pow(kubust[l][3],2)));
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
			
			while (kubust[name1][13] < 0) {
				name1 = (int) Nt-kubust[name1][13];
				multiplicity++;
			}
			while (kubust[name2][13] < 0) {
				name2 = (int) Nt-kubust[name2][13];
				multiplicity++;
			}
			
			//Hoehere Multiplizitaeten auf Bindung mit dem Centre of Mass des stablieren Binaries pruefen
			if (multiplicity>2) {
				rx = kubust[name1][1]-kubust[name2][1];
				ry = kubust[name1][2]-kubust[name2][2];
				rz = kubust[name1][3]-kubust[name2][3];
				rtemp = sqrt(rx*rx+ry*ry+rz*rz);
				
				vx = kubust[name1][4]-kubust[name2][4];
				vy = kubust[name1][5]-kubust[name2][5];
				vz = kubust[name1][6]-kubust[name2][6];
				
				Ebintemp = 0.5*kubust[name1][10]*kubust[name2][10]/(kubust[name1][10]+kubust[name2][10])*(vx*vx+vy*vy+vz*vz) - G*kubust[name1][10]*kubust[name2][10]/rtemp;
				
				atemp = -0.5*G*kubust[name1][10]*kubust[name2][10]/Ebintemp;
				
				etemp = sqrt(pow(1.0 - rtemp/atemp,2)+ pow(rx*vx+ry*vy+rz*vz,2)/(G*(kubust[name1][10]+kubust[name2][10])*atemp));
				
			} else {
				etemp = 0.0;
			}	
			
			if (etemp<1.0) {
				rx = kubust[name1][1]-kubust[name2][1];
				ry = kubust[name1][2]-kubust[name2][2];
				rz = kubust[name1][3]-kubust[name2][3];
				rtemp = sqrt(rx*rx+ry*ry+rz*rz);
				
				vx = kubust[name1][4]-kubust[name2][4];
				vy = kubust[name1][5]-kubust[name2][5];
				vz = kubust[name1][6]-kubust[name2][6];
				
				Ebintemp = 0.5*kubust[name1][10]*kubust[name2][10]/(kubust[name1][10]+kubust[name2][10])*(vx*vx+vy*vy+vz*vz) - G*kubust[name1][10]*kubust[name2][10]/rtemp;
				
				atemp = -0.5*G*kubust[name1][10]*kubust[name2][10]/Ebintemp;
				
				if (binoutputt) printf("%g\t%i\t%i\t%i\t-->\t%i (distance from centre: %fpc, semi-major axis: %fpc)\n", binarray[l][0], name1, name2, multiplicity,Nt+ll,sqrt(pow(kubust[name1][1],2)+pow(kubust[name1][2],2)+pow(kubust[name1][3],2)),atemp);
				kubust[name1][13] = -1*ll;
				kubust[name2][13] = -1*ll;
				
				
				//CoM erstellen
				kubust[Nt+ll][0] = -1.0;
				kubust[Nt+ll][1] = (kubust[name1][1]*kubust[name1][10]+kubust[name2][1]*kubust[name2][10])/(kubust[name1][10]+kubust[name2][10]);
				kubust[Nt+ll][2] = (kubust[name1][2]*kubust[name1][10]+kubust[name2][2]*kubust[name2][10])/(kubust[name1][10]+kubust[name2][10]);
				kubust[Nt+ll][3] = (kubust[name1][3]*kubust[name1][10]+kubust[name2][3]*kubust[name2][10])/(kubust[name1][10]+kubust[name2][10]);
				kubust[Nt+ll][4] = (kubust[name1][4]*kubust[name1][10]+kubust[name2][4]*kubust[name2][10])/(kubust[name1][10]+kubust[name2][10]);
				kubust[Nt+ll][5] = (kubust[name1][5]*kubust[name1][10]+kubust[name2][5]*kubust[name2][10])/(kubust[name1][10]+kubust[name2][10]);
				kubust[Nt+ll][6] = (kubust[name1][6]*kubust[name1][10]+kubust[name2][6]*kubust[name2][10])/(kubust[name1][10]+kubust[name2][10]);
				kubust[Nt+ll][7] = -1.0;
				kubust[Nt+ll][8] = kubust[name1][8]+kubust[name2][8];
				kubust[Nt+ll][9] = -1.0;
				kubust[Nt+ll][10] = kubust[name1][10]+kubust[name2][10];
				kubust[Nt+ll][11] = 0.0;
				kubust[Nt+ll][12] = 0.5*kubust[Nt+ll][10]*(pow(kubust[Nt+ll][4],2)+pow(kubust[Nt+ll][5],2)+pow(kubust[Nt+ll][6],2));
				kubust[Nt+ll][13] = 0.0;
				if ((kubust[name1][14]) && (kubust[name1][14] < 2.0) && (kubust[name2][14]) && (kubust[name2][14] < 2.0)) {
					kubust[Nt+ll][14] = 1.0;
					kubust[name1][14] = 2.0;
					kubust[name2][14] = 2.0;
					if (ausgabetypt) printf("BINARY ESCAPE:\t%i\t%i\n", name1, name2);
				} else {
					kubust[Nt+ll][14] = 0.0;
				}
				
				*Ntott = *Ntott + 1;
				ll++;
				multiplicityarrayt[multiplicity]++;
			}
		}
	}
	
	for (l=Nt;l<*Ntott;l++) {
		if ((kubust[l][0]) && (kubust[l][13] >= 0)) {
			for (ll=0;ll<*Ntott;ll++) {
				if ((l-ll) && (kubust[ll][0]) && (kubust[ll][13] >= 0.0)) {
					rtemp = sqrt(pow(kubust[l][1]-kubust[ll][1],2)+pow(kubust[l][2]-kubust[ll][2],2)+pow(kubust[l][3]-kubust[ll][3],2));
					if (rtemp) kubust[l][11] += -G*kubust[l][10]*kubust[ll][10]/rtemp;
				}
			}
		}
	}
	
	for (l=0;l<Nt;l++) free (binarray[l]);
	free(binarray);
	
	return 1;
}

int energieberechner(int Nt, double **kubust, int statust, double *Etemp, double mvxt, double mvyt, double mvzt){
	int l;
	Etemp[0] = 0.0; //Ebin
	Etemp[1] = 0.0; //Ekin
	Etemp[2] = 0.0; //Eges
	
	for (l=0;l<Nt;l++) {
		if ((kubust[l][0]) && (kubust[l][13] < statust) && (kubust[l][13] >= 0)) {
			Etemp[0] += kubust[l][11];
			Etemp[1] += 0.5*kubust[l][10]*(pow(kubust[l][4]-mvxt,2)+pow(kubust[l][5]-mvyt,2)+pow(kubust[l][6]-mvzt,2));
		}
	}
	Etemp[0] /= 2.0;
	Etemp[2] = Etemp[0] + Etemp[1];
	
	return 1;
}

int statussetzer(int Nt, double **kubust, double rtidet, double mxt, double myt, double mzt, int *nt, int *nbint, double mvxt, double mvyt, double mvzt, double rh, double omega, int schnellanalyse){
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
		if ((kubust[l][0]) && (kubust[l][13] >= 0.0)){
			rtemp = sqrt(pow(kubust[l][1]-mxt,2)+pow(kubust[l][2]-myt,2)+pow(kubust[l][3]-mzt,2));
			if (rtemp < rtidet) {
				kubust[l][13] = 0.0;
				nt[0]++ ;
				mtemp += kubust[l][10];
				ltemp = floor(1.0*rtemp/rtidet*gridpoints);
				massprofile[ltemp][1] += kubust[l][10];
			} else if (rtemp < 2.0*rtidet) {
				kubust[l][13] = 2.0;
				nt[2]++;
			} else {
				kubust[l][13] = 3.0;
				nt[3]++;
			}
		} else if ((kubust[l][0]) && (kubust[l][13] < 0.0)) {
			rtemp = sqrt(pow(kubust[l][1]-mxt,2)+pow(kubust[l][2]-myt,2)+pow(kubust[l][3]-mzt,2));
			if (rtemp < rtidet) {
				nbint[0]++ ;
				mtemp += kubust[l][10];
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
	
	
	//Potential escapers suchen & Bindunsenergien fuer Schnellanalyse berechnen
	ERT = -1.5*G*mtemp/rtidet;
	
	for (l=0;l<Nt;l++) {
		if ((kubust[l][0]) && (kubust[l][13] == 0.0)){
			if (schnellanalyse) {
				rtemp2 = sqrt(pow(kubust[l][1]-mxt,2)+pow(kubust[l][2]-myt,2)+pow(kubust[l][3]-mzt,2));
				ltemp = floor(1.0*rtemp2/rtidet*gridpoints);
				kubust[l][11] = -G*kubust[l][10]*mtemp/rtemp2;
				//kubust[l][11] = - G*mtemp*kubust[l][10]*(pow(1.305*rtemp2/rh,3)/pow(1.0+pow(1.305*rtemp2/rh,2),1.5)); //assumed Plummer profile for binding energy!
			}
			
			etemp = 0.5*kubust[l][10]*(pow(kubust[l][4]-mvxt,2)+pow(kubust[l][5]-mvyt,2)+pow(kubust[l][6]-mvzt,2)) + 0.5*kubust[l][10]*omega*omega*(pow(kubust[l][3]-mzt,2) - 2.0*pow(kubust[l][1]-mxt,2)) + kubust[l][11]; //assumed near-field approximation for binding energy!
			
			
			if (etemp > kubust[l][10]*ERT) {
				kubust[l][13] = 1.0;
				nt[0]-- ;
				nt[1]++ ;
			}
		} 
		
	}
	
	return 1;
}

int rtidesuche(double *rtidet, int Nt, double **kubust, double omegat, double mxt, double myt, double mzt, double dphit, double *RGt){
	int l;
	double rtidetemp, Mtemp = 0.0;
	if (*rtidet > 1.0) { 
		for (l=0;l<Nt;l++) {
			if (kubust[l][0]) Mtemp += kubust[l][10];
		}
	} else {
		Mtemp = Nt;
	}
	do{
		rtidetemp = pow(G*Mtemp/sqrt(pow((dphit+omegat*omegat),2)),1.0/3.0);
		Mtemp = 0.0;
		for (l=0;l<Nt;l++) {
			if ((kubust[l][0]) && (sqrt(pow(kubust[l][1]-mxt,2)+pow(kubust[l][2]-myt,2)+pow(kubust[l][3]-mzt,2)) < rtidetemp)) {
				Mtemp += kubust[l][10];
			}
		}
		*rtidet = pow(G*Mtemp/sqrt(pow((dphit+omegat*omegat),2)),1.0/3.0);
		printf("+");
	} while (sqrt(pow(*rtidet-rtidetemp,2))/ *rtidet > 0.01);
	printf("\n");
	
	return 1;
}

int omegasuche(double *RGt, double *VGt, double *omegat, double *Lt) {
	if (RGt[0]) {
		*Lt = sqrt(pow(RGt[1]*VGt[2]-RGt[2]*VGt[1],2)+pow(RGt[2]*VGt[0]-RGt[0]*VGt[2],2)+pow(RGt[0]*VGt[1]-RGt[1]*VGt[0],2));
		*omegat = *Lt/(RGt[0]*RGt[0]+RGt[1]*RGt[1]+RGt[2]*RGt[2]);
	} 
	
	return 1;
}

int zentrumsuche(int Nt, double **kubust, int **nachbarlistet, double *mxt, double *myt, double *mzt, int densitylevelt){ 
	int l, ll;
	double mtemp, rhotemp, rhogestemp = 0.0;
	double softening = 0.001;
	
	//float rhotempliste[Nt];
	
	*mxt = 0; *myt = 0; *mzt = 0;
	
	for (l=0;l<Nt;l++) {
		if ((kubust[l][0]) && (kubust[l][13]<3)) {
			mtemp = 0.0;
			for (ll=0;ll<densitylevelt;ll++) {
				mtemp += kubust[nachbarlistet[l][ll]][10];//Mass
				//mtemp += kubust[nachbarlistet[l][ll]][8];//Luminosity
			}
			rhotemp = mtemp/pow(pow(kubust[nachbarlistet[l][densitylevelt]][1]-kubust[l][1],2)+pow(kubust[nachbarlistet[l][densitylevelt]][2]-kubust[l][2],2)+pow(kubust[nachbarlistet[l][densitylevelt]][3]-kubust[l][3],2)+pow(softening,2),1.5);
		} else {
			rhotemp = 0;
		}
		
		//rhotempliste[l]=rhotemp;
		
		*mxt += kubust[l][1]*rhotemp;
		*myt += kubust[l][2]*rhotemp;
		*mzt += kubust[l][3]*rhotemp;
		rhogestemp += rhotemp; 
	}
	//qsort(rhotempliste, Nt, sizeof(float), (void *)cmpmy2);
	//printf("\n%g\t%g\t%g\t%g\n", rhotempliste[0], rhotempliste[1], rhotempliste[2], rhotempliste[3]);
	
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

int bindungsanalysefast(int Nt, double **kubust, int **nachbarlistet, int nachbarnt, double mxt, double myt, double mzt) {
	int l, ll, i, j, namei, namej, namel;
	double rtemp;
	double ri, rj;
	
	float rtempliste[Nt][2];
	float rtempliste2[nachbarnt][2];
	
	for (l=0;l<Nt;l++) {
		if (kubust[l][0]) {
			rtemp = sqrt(pow(kubust[l][1] - mxt,2)+pow(kubust[l][2] - myt,2)+pow(kubust[l][3] - mzt,2));
		} else {
			rtemp = 1000000.0;
		}
		rtempliste[l][0] = rtemp;
		rtempliste[l][1] = l;
	}
	
	qsort(rtempliste, Nt, sizeof(rtempliste[0]), (void *)cmpmy);
	
	for (l=0;l<Nt;l++) {
		i=1; //oben
		j=1; //unten
		namel = (int) rtempliste[l][1];
		if (kubust[namel][0]) {
			for (ll=1;ll<=nachbarnt;ll++) {
				//ll ist die Position, in die ins Nachbararray geschrieben wird.
				//dann einen weiteren zaehler, z.B. i, der nach oben und einen anderen, z.B. j, der nach unten zaehlt.
				//dann beide Radien ausrechnen und vergleichen.
				//kleineren der beiden in das Nachbararray schreiben und ll erhoehen.
				if (l+i<Nt) {
					namei = (int) rtempliste[l+i][1];
					ri = sqrt(pow(kubust[namel][1]-kubust[namei][1],2)+pow(kubust[namel][2]-kubust[namei][2],2)+pow(kubust[namel][3]-kubust[namei][3],2));
				} else {
					ri = 10000000.0;
				}
				
				if (l-j>=0) {
					namej = (int) rtempliste[l-j][1];
					rj = sqrt(pow(kubust[namel][1]-kubust[namej][1],2)+pow(kubust[namel][2]-kubust[namej][2],2)+pow(kubust[namel][3]-kubust[namej][3],2));
				} else {
					rj = 10000000.0;
				}
				
				if (ri<rj) {
					//nachbarlistet[namel][ll-1] = namei;
					rtempliste2[ll-1][1] = namei;
					rtempliste2[ll-1][0] = ri;
					i++;
					if (i>j+6) j++;
				} else {
					//nachbarlistet[namel][ll-1] = namej;
					rtempliste2[ll-1][1] = namej;
					rtempliste2[ll-1][0] = rj;
					j++;
					if (j>i+6) i++;
				}
			}
			qsort(rtempliste2, nachbarnt, sizeof(rtempliste[0]), (void *)cmpmy);
			for (ll=0;ll<nachbarnt;ll++){
				nachbarlistet[namel][ll] = rtempliste2[ll][1];
			}
		}
	}
	return 1;
}

int bindungsanalyse(int Nt, double **kubust, int **nachbarlistet, int nachbarnt){
	int l, ll;
	double Ebintemp, rtemp;
	
	float rtempliste[Nt][2];
	
	for (l=0;l<Nt;l++) {
		if (kubust[l][0]) {
			//printf("%i\n",l);
			rtempliste[l][0] = 1000000.0;
			rtempliste[l][1] = l;
			
			for (ll=0;ll<Nt;ll++) {
				if ((l-ll) && (kubust[ll][0])) {
					//radien berechnen
					rtemp = sqrt(pow(kubust[l][1]-kubust[ll][1],2)+pow(kubust[l][2]-kubust[ll][2],2)+pow(kubust[l][3]-kubust[ll][3],2));
					rtempliste[ll][0] = rtemp;
					rtempliste[ll][1] = ll;
				} else {
					rtempliste[ll][0] = 1000000.0;
					rtempliste[ll][1] = ll;
				}	
			}
			qsort(rtempliste, Nt, sizeof(rtempliste[0]), (void *)cmpmy);
			for (ll=0;ll<nachbarnt;ll++) {
				nachbarlistet[l][ll] = rtempliste[ll][1];
			}
			
			//Bindungsenergien aufsummieren bis ca 0.01% Genauigkeit
			kubust[l][11] = 0.0;
			ll=0;	
			do{
				Ebintemp = -G*kubust[l][10]*kubust[(int) rtempliste[ll][1]][10]/rtempliste[ll][0];
				kubust[l][11]+=Ebintemp;
				ll++;
			} while ((ll<Nt) && (Ebintemp/kubust[l][11]>0.0000001));
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

int readin0(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott) {  //McLuster;
	
	int i, k;
	double m, x[3], v[3], m0, epoch, spin, rad, lum; 
	
	i=0;
	while(fscanf(posdatei,"%lf %lf %lf %lf %lf %lf %lf %lf %i %lf %lf %lf %lf\n",&m,&x[0],&x[1],&x[2],&v[0],&v[1],&v[2],&m0, &k, &epoch, &spin, &rad, &lum) == 13) {
		kubust[i][0] = i+1;	
		kubust[i][1] = x[0];              //in pc
		kubust[i][2] = x[1];
		kubust[i][3] = x[2];
		kubust[i][4] = v[0];             //in km/s
		kubust[i][5] = v[1];
		kubust[i][6] = v[2];
		kubust[i][7] = k;
		kubust[i][8] = lum;
		kubust[i][9] = rad;
		kubust[i][10] = m;                 //in Sonnenmassen
		kubust[i][12] = 0.5*m*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));  //Ekin
		
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

int readin1(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, int ausgabetypt, int kubussizetott) {  //NBODY4 32bit;
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
		fread(&b,sizeof(int),1,posdatei);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (ausgabetypt) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[4]);
	
	double as[a[2]];
	fread(as,sizeof(double),a[2],posdatei);
	
	if (ausgabetypt) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posdatei);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posdatei);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posdatei);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	fread(phi4,sizeof(double),a[1],posdatei);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posdatei);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posdatei); 
	double *lsev;
	lsev = (double *)malloc(a[1]*sizeof(double));
	fread(lsev,sizeof(float),a[1],posdatei);
	double *rsev;
	rsev = (double *)malloc(a[1]*sizeof(double)); 
	fread(rsev,sizeof(float),a[1],posdatei);
	
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
	 * Einzelsterne in Kubus schreiben *
	 ***********************************/
	
	for(l=0;l<a[4];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < kubussizetott)) {
			kubust[name[l]-1][0] = name[l];	
			kubust[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			kubust[name[l]-1][2] = xs[3*l+1]*as[2];
			kubust[name[l]-1][3] = xs[3*l+2]*as[2];
			kubust[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			kubust[name[l]-1][5] = vs[3*l+1]*as[11];
			kubust[name[l]-1][6] = vs[3*l+2]*as[11];
			kubust[name[l]-1][7] = kstar[l];
			kubust[name[l]-1][8] = lsev[l];
			kubust[name[l]-1][9] = rsev[l];
			kubust[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			kubust[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
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

int readin2(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott) {  //NBODY6;
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
		fread(&b,sizeof(int),1,posdatei);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (ausgabetypt) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[4]);
	
	double as[a[2]];
	fread(&as,sizeof(double),a[2],posdatei);
	
	if (ausgabetypt) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posdatei);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posdatei);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posdatei);
	double *rsev;
	rsev = (double *)malloc(a[1]*sizeof(double)); 
	fread(rsev,sizeof(double),a[1],posdatei);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posdatei);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posdatei); 
	double *lsev;
	lsev = (double *)malloc(a[1]*sizeof(double));
	fread(lsev,sizeof(double),a[1],posdatei);
	int end[1];
	fread(&end,sizeof(int),1,posdatei);
	//fread(&b,sizeof(int),1,posdatei);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	//fread(phi4,sizeof(double),a[1],posdatei);
	
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
	 * Einzelsterne in Kubus schreiben *
	 ***********************************/
	
	for(l=0;l<a[3];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < kubussizetott)) {
			kubust[name[l]-1][0] = name[l];	
			kubust[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			kubust[name[l]-1][2] = xs[3*l+1]*as[2];
			kubust[name[l]-1][3] = xs[3*l+2]*as[2];
			kubust[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			kubust[name[l]-1][5] = vs[3*l+1]*as[11];
			kubust[name[l]-1][6] = vs[3*l+2]*as[11];
			kubust[name[l]-1][7] = kstar[l];
			kubust[name[l]-1][8] = lsev[l];
			kubust[name[l]-1][9] = rsev[l];
			kubust[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			kubust[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
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

int readin3(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott) {  //NBODY6;
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
		fread(&b,sizeof(int),1,posdatei);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (ausgabetypt) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[4]);
	
	a[2] = 30; //30 Ausgabeparameter
	double as[a[2]];
	fread(&as,sizeof(double),a[2],posdatei);
	
	if (ausgabetypt) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posdatei);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posdatei);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posdatei);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	fread(phi4,sizeof(double),a[1],posdatei);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posdatei);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posdatei); 
	float *lsev;
	lsev = (float *)malloc(a[1]*sizeof(float));
	fread(lsev,sizeof(float),a[1],posdatei);
	float *rsev;
	rsev = (float *)malloc(a[1]*sizeof(float)); 
	fread(rsev,sizeof(float),a[1],posdatei);
	
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
	 * Einzelsterne in Kubus schreiben *
	 ***********************************/
	
	for(l=0;l<a[3];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < kubussizetott)) {
			kubust[name[l]-1][0] = name[l];	
			kubust[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			kubust[name[l]-1][2] = xs[3*l+1]*as[2];
			kubust[name[l]-1][3] = xs[3*l+2]*as[2];
			kubust[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			kubust[name[l]-1][5] = vs[3*l+1]*as[11];
			kubust[name[l]-1][6] = vs[3*l+2]*as[11];
			kubust[name[l]-1][7] = kstar[l];
			kubust[name[l]-1][8] = lsev[l];
			kubust[name[l]-1][9] = rsev[l];
			kubust[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			kubust[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
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

int readin4(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott) {  //NBODY6;
	int l;	
	int a[6], b[2];
	fread(&b,sizeof(int),1,posdatei);
	fread(&a,sizeof(int),6,posdatei);
	if (ausgabetypt) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	fread(&b,sizeof(int),2,posdatei);
	//fread(&b,sizeof(int),1,posdatei);
	double as[a[2]];
	fread(&as,sizeof(double),a[2],posdatei);
	
	if (ausgabetypt) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posdatei);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posdatei);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posdatei);
	double *rsev;
	rsev = (double *)malloc(a[1]*sizeof(double)); 
	fread(rsev,sizeof(double),a[1],posdatei);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posdatei);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posdatei); 
	double *lsev;
	lsev = (double *)malloc(a[1]*sizeof(double));
	fread(lsev,sizeof(double),a[1],posdatei);
	int end[1];
	fread(&end,sizeof(int),1,posdatei);
	fread(&b,sizeof(int),1,posdatei);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	//fread(phi4,sizeof(double),a[1],posdatei);
	
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
	 * Einzelsterne in Kubus schreiben *
	 ***********************************/
	
	for(l=0;l<a[3];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < kubussizetott)) {
			kubust[name[l]-1][0] = name[l];	
			kubust[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			kubust[name[l]-1][2] = xs[3*l+1]*as[2];
			kubust[name[l]-1][3] = xs[3*l+2]*as[2];
			kubust[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			kubust[name[l]-1][5] = vs[3*l+1]*as[11];
			kubust[name[l]-1][6] = vs[3*l+2]*as[11];
			kubust[name[l]-1][7] = kstar[l];
			kubust[name[l]-1][8] = lsev[l];
			kubust[name[l]-1][9] = rsev[l];
			kubust[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			kubust[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
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

int readin5(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott) {  //NBODY4new 64bit;
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
		fread(&b,sizeof(int),1,posdatei);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (ausgabetypt) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[4]);
	
	double as[a[2]];
	fread(&as,sizeof(double),a[2],posdatei);
	
	if (ausgabetypt) for (l=0;l<a[2];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[1]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[1],posdatei);
	double *xs;
	xs = (double *)malloc(3*a[1]*sizeof(double));
	fread(xs,sizeof(double),3*a[1],posdatei);
	double *vs;
	vs = (double *)malloc(3*a[1]*sizeof(double));
	fread(vs,sizeof(double),3*a[1],posdatei);
	double *phi4;
	phi4 = (double *)malloc(a[1]*sizeof(double));
	fread(phi4,sizeof(double),a[1],posdatei);
	int *name;
	name = (int *)malloc(a[1]*sizeof(int));
	fread(name,sizeof(int),a[1],posdatei);
	int *kstar;
	kstar = (int *)malloc(a[1]*sizeof(int));
	fread(kstar,sizeof(int),a[1],posdatei); 
	double *lsev;
	lsev = (double *)malloc(a[1]*sizeof(double));
	fread(lsev,sizeof(double),a[1],posdatei);
	double *rsev;
	rsev = (double *)malloc(a[1]*sizeof(double)); 
	//fread(rsev,sizeof(float),a[1],posdatei);
	//int end[1];
	//fread(&end,sizeof(int),1,posdatei);
	//fread(&b,sizeof(int),1,posdatei);
	
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
	 * Einzelsterne in Kubus schreiben *
	 ***********************************/
	
	for(l=0;l<a[3];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < kubussizetott)) {
			kubust[name[l]-1][0] = name[l];	
			kubust[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			kubust[name[l]-1][2] = xs[3*l+1]*as[2];
			kubust[name[l]-1][3] = xs[3*l+2]*as[2];
			kubust[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			kubust[name[l]-1][5] = vs[3*l+1]*as[11];
			kubust[name[l]-1][6] = vs[3*l+2]*as[11];
			kubust[name[l]-1][7] = kstar[l];
			kubust[name[l]-1][8] = lsev[l];
			kubust[name[l]-1][9] = rsev[l];
			kubust[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			kubust[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
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

int readin6(FILE *posdatei, int *Ntott, int *Nt, double *timet, double *tbart, double *rtidet, double *rbart, double *vbart, double *mbart, double *omegat, double *Etott, double *rcoret, double **kubust, double *RGt, double *VGt, int ausgabetypt, int kubussizetott) {  //NBODY6;
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
		fread(&b,sizeof(int),1,posdatei);
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = a[4];
		a[4] = a[5];
		a[5] = b;
		if (ausgabetypt) printf("%i %i %i %i %i %i\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	} while (a[0] != a[5]);
	fread(&b,sizeof(int),1,posdatei);
	fread(&b,sizeof(int),1,posdatei);
	fread(&b,sizeof(int),1,posdatei);
	
	double as[a[3]];
	fread(&as,sizeof(double),a[3],posdatei);
	
	if (ausgabetypt) for (l=0;l<a[3];l++) printf("%i:  %lf\n",l+1,as[l]);
	
	double *bodys;
	bodys = (double *)malloc(a[2]*sizeof(double));
	if (NULL == bodys) printf ("\nDu NULL\n");
	fread(bodys,sizeof(double),a[2],posdatei);
	double *xs;
	xs = (double *)malloc(3*a[2]*sizeof(double));
	fread(xs,sizeof(double),3*a[2],posdatei);
	double *vs;
	vs = (double *)malloc(3*a[2]*sizeof(double));
	fread(vs,sizeof(double),3*a[2],posdatei);
	double *rsev;
	rsev = (double *)malloc(a[2]*sizeof(double)); 
	fread(rsev,sizeof(double),a[2],posdatei);
	int *name;
	name = (int *)malloc(a[2]*sizeof(int));
	fread(name,sizeof(int),a[2],posdatei);
	int *kstar;
	kstar = (int *)malloc(a[2]*sizeof(int));
	fread(kstar,sizeof(int),a[2],posdatei); 
	double *lsev;
	lsev = (double *)malloc(a[2]*sizeof(double));
	fread(lsev,sizeof(double),a[2],posdatei);
	int end[1];
	fread(&end,sizeof(int),1,posdatei);
	//fread(&b,sizeof(int),1,posdatei);
	double *phi4;
	phi4 = (double *)malloc(a[2]*sizeof(double));
	//fread(phi4,sizeof(double),a[2],posdatei);
	
	RGt[0] = as[20]*as[2];
	RGt[1] = as[21]*as[2];
	RGt[2] = as[22]*as[2];
	VGt[0] = as[26]*as[11];
	VGt[1] = as[27]*as[11];
	VGt[2] = as[28]*as[11];
	
	*tbart = as[10];
	*timet = as[9];
	*Ntott = a[2];
	*Nt = a[4];
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
	 * Einzelsterne in Kubus schreiben *
	 ***********************************/
	
	for(l=0;l<a[4];l++) {
		if ((name[l]-1 >= 0) && (name[l]-1 < kubussizetott)) {
			kubust[name[l]-1][0] = name[l];	
			kubust[name[l]-1][1] = xs[3*l]*as[2];              //Umrechnung in pc
			kubust[name[l]-1][2] = xs[3*l+1]*as[2];
			kubust[name[l]-1][3] = xs[3*l+2]*as[2];
			kubust[name[l]-1][4] = vs[3*l]*as[11];             //Umrechnung in (km/s)^2
			kubust[name[l]-1][5] = vs[3*l+1]*as[11];
			kubust[name[l]-1][6] = vs[3*l+2]*as[11];
			kubust[name[l]-1][7] = kstar[l];
			kubust[name[l]-1][8] = lsev[l];
			kubust[name[l]-1][9] = rsev[l];
			kubust[name[l]-1][10] = bodys[l]*as[3];                 //Umrechnung in Sonnenmassen
			kubust[name[l]-1][12] = 0.5*bodys[l]*as[3]*(pow(vs[3*l]*as[11],2)+pow(vs[3*l+1]*as[11],2)+pow(vs[3*l+2]*as[11],2));  //Ekin
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
