
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*
  bondinganalyze - code to analyze the bonding environment of all atoms
  Code originally taken from structfactanalyze and modified from it. 

  Neighbour list calculation originally taken from PARCAS analyzeint.c

  NOTE: All arrays here start from 0, not 1 as in parcas !!

*/

typedef int logical;
typedef double real;


struct ngbr {
  int i;
  real r;
  real s;  /* Strength of neighbour */
};



/* Debug defines */

logical debug,maxdebug;

#define DEBUG(var) if (debug) fprintf(stderr, #var " = %g\n",(double) (var));fflush(stdout); fflush(stderr);
#define DEBUGS(str) if (debug) fprintf(stderr,str); fflush(stderr); fflush(stdout)
#define DEBUGSR(str,var) if (debug) fprintf(stderr, "%s " #var " = %g\n",str,(double) (var));fflush(stdout); fflush(stderr);
#define DEBUGSS(str,str2) if (debug) fprintf(stderr, "%s " #str2 " = %s\n",str,str2);fflush(stdout); fflush(stderr);
#define DEBUGSRR(str,var,var2) if (debug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g\n",str,(double) (var),(double) (var2));fflush(stdout); fflush(stderr);
#define DEBUGSRRR(str,var,var2,var3) if (debug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g\n",str,(double) (var),(double) (var2),(double) (var3));fflush(stdout); fflush(stderr);
#define DEBUGSRRRR(str,var,var2,var3,var4) if (debug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g " #var4 " = %g\n",str,(double) (var),(double) (var2),(double) (var3),(double) (var4));fflush(stdout); fflush(stderr);
#define DEBUGSRRRRR(str,var,var2,var3,var4,var5) if (debug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g " #var4 " = %g " #var5 " = %g\n",str,(double) (var),(double) (var2),(double) (var3),(double) (var4),(double) (var5));fflush(stdout); fflush(stderr);


#define MAXDEBUG(var) if (maxdebug) fprintf(stderr, #var " = %g\n",(double) (var));fflush(stdout); fflush(stderr);
#define MAXDEBUGS(str) if (maxdebug) fprintf(stderr,str); fflush(stderr); fflush(stdout)
#define MAXDEBUGSR(str,var) if (maxdebug) fprintf(stderr, "%s " #var " = %g\n",str,(double) (var));fflush(stdout); fflush(stderr);
#define MAXDEBUGSS(str,str2) if (maxdebug) fprintf(stderr, "%s " #str2 " = %s\n",str,str2);fflush(stdout); fflush(stderr);
#define MAXDEBUGSRR(str,var,var2) if (maxdebug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g\n",str,(double) (var),(double) (var2));fflush(stdout); fflush(stderr);
#define MAXDEBUGSRRR(str,var,var2,var3) if (maxdebug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g\n",str,(double) (var),(double) (var2),(double) (var3));fflush(stdout); fflush(stderr);
#define MAXDEBUGSRRRR(str,var,var2,var3,var4) if (maxdebug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g " #var4 " = %g\n",str,(double) (var),(double) (var2),(double) (var3),(double) (var4));fflush(stdout); fflush(stderr);
#define MAXDEBUGSRRRRR(str,var,var2,var3,var4,var5) if (maxdebug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g " #var4 " = %g " #var5 " = %g\n",str,(double) (var),(double) (var2),(double) (var3),(double) (var4),(double) (var5));fflush(stdout); fflush(stderr);


/* Other defines */

int MAXAT;

#define MAXNGBR 50
#define MAXTYPE 6 
int itypemin,itypemax;  /* Actual minimum and maximum */

#define WSCELLMAX 6

#define True 1
#define False 0
#define pi 3.1415926535897

/*
  Function protypes 
*/


void bondinganalyze(real *x0,
		    int *itype,
		    int *wasliq,
		    real *box,
		    real *pbc,
		    real cutR[MAXTYPE][MAXTYPE],
		    real cutS[MAXTYPE][MAXTYPE],
		    real cutmin, // CB
		    real cutmax, // CB
		    int nat,
		    real a,
		    real dr,
		    real time,
		    real printshort,
		    logical printbonds,
		    logical printangles,
		    logical writenbonds,
		    logical onlyliq,
		    logical noneidrstat,
		    int printnot,
		    int printnot2[MAXTYPE],
		    int roughnessZ1,
		    int roughnessZ2,
		    real roughnessz0,
		    real printeffectivea,
		    real effectivearmax,
		    real effectiveafact		    
		    );


void getngbrs(int iat,
	      real *x0,
	      real *box,
	      real *pbc,
	      int N,
	      real Rneicut,
	      struct ngbr *ngbrs,
	      int *nngot);

real getangle(int i,
	      int j,
	      int k,
	      real *x0,
	      real *box,
	      real *pbc);

struct point {
   double x;
   double y;
   double z;
};


logical details=False;  /* Shall we give details of the analysis as well ? */
int binarycut,binaryt1,binaryt2;  /* Shall we do binary type statistics ? */
double SROr1,SROr2; /* Shall we do SRO analysis ? Implies binary analysis as well */
double cutmin, cutmax; // CB

char name[MAXTYPE][8];
int *in,*itype,*wasliq;

int lattice;
char string[80];

#define FCC 1
#define BCC 2
#define DIA 3
#define HCP 4

main(argc,argv)
int argc;
char **argv;
{

   int iarg;
   int i,j,ii,nat;

   FILE *fp;
   char file[50];
   int nline;

   int stringcol,xdatacol,typedatacol;

   logical usetime,newtime,firstdatafound;
   logical iscoordline,printeffectivea=False;
   logical discardout=False,disc,enforcepbc=False,enforce;
   int ntime,natexpectedp,natexpected;
   int ndiscard,ndiscardx,ndiscardy,ndiscardz;
   int nenforcepbc,nenforcepbcx,nenforcepbcy,nenforcepbcz;
   int ndiscardxfar,ndiscardyfar,ndiscardzfar;
   int printnot,printnot2[MAXTYPE];
   int roughnessZ1,roughnessZ2;
   real roughnessz0;
   double printshort,effectiveafact,effectivearmax;
   logical printbonds,printangles,writenbonds,onlyliq,noneidrstat;

   double xx,yy,zz;

   char buf[160],buforig[160],c1;
   char *argbuf,arg[40][40];
   int narg;

   /* structfactanalyze interface stuff */
   double *x0;
   real box[3],pbc[3];

   double casctime,prevtime;
   double forceminx=-1e30,forceminy=-1e30,forceminz=-1e30;

   /* Tersoff cutoff function values */
   double cutR[MAXTYPE][MAXTYPE],cutS[MAXTYPE][MAXTYPE],cutSmax,cutRgen,cutSgen;

   double dr;
   double newbox[3];
   
   /* 
      Equilibrium unit cell size, many other variables are defined from
      this
   */

   double unitcella;
   int unitcelln;

   if ((argc>1 && strcmp(argv[1],"-h")==0.0) || argc < 2) {
      printf("Usage: bondinganalyze [options] file [stringcol [string [xdatacol [xsize [ysize [zsize]]]]]]\n");
      printf("                                 Material always treated as alloy\n");
      printf("\n"); 
      printf("Options: -debug                  debug\n");
      printf("         -maxdebug               debug even more\n");
      printf("         -d                      give detailed information\n");
      printf("         -a a                    Unit cell a\n");
      printf("         -typecol n              Set atom type column to n\n");
      printf("         -cut R S                Set overall Tersoff cutoff function R and S\n");
      printf("                                 OVERRIDES ANY PRECEDING -cutij option! \n");
      printf("         -cutij i j R S          Set pair Tersoff cutoff function R and S\n");
      printf("         -n ncells               Number of unit cells, gives a\n");
      printf("         -writenbonds            Print time development of number of atoms with n bonds \n");
      printf("         -printangles            Print all bonding angles \n");
      printf("         -printshort r           Print all bonds shorter than r\n");
      printf("         -printbonds             Print all bonds shorter than Tersoff cutoff\n");
      printf("         -printa rmax f          Print effective a for all atoms, a = ave r bond (<rmax)*f\n");
      printf("         -printnot n             Print atoms with number of bonds not=n\n");
      printf("         -printnot2 t n          Print atoms with number of bonds not=n for type t\n");
      printf("         -roughness Z1 Z2 z0     Do surface roughness analysis, assume atoms with Z1<Z<Z2 above z0 are surface atoms\n"); 
      printf("         -dr dr                  Set dr for pair correl. function calc.\n");
      printf("         -binary n t1 t2         Print statistics of binary pairs of type t1 t2\n");
      printf("         -SROpair r1 r2 t1 t2    Do SRO analysis between r1-r2 for types t1 t2 only\n"); 
      printf("         -shellcut cutmin cutmax  Do neighbour analysis (in shell between cutmin->cutmax CB\n"); 
      printf("                                 [J.M. Cowley, Phys. Rev. 77(5), 1950, p. 669]\n");
      printf("         -pbc 0/1 0/1 0/1        (Un)Set periodic boundaries: 1 1 1 assumed\n");
      printf("         -enforcepbc             Enforce periodics at periodic boundaries\n");
      printf("         -discardout             Discard atoms outside cell completely\n");
      printf("         -onlyliq                Only analyze 'liquid' atoms (column 7=1 in xyz) for bonding\n");
      printf("         -noneidrstat            Do not output neidrstat files\n");
      printf("\n");
      printf("         -fcc                    Choose parameters apprpriate for FCC [Default]\n");
      printf("         -bcc                    Choose parameters apprpriate for BCC \n");
      printf("         -dia                    Choose parameters apprpriate for DIA \n");
      printf("\n");
      printf("\n");
      printf("Data is taken as x y z beginning from column xdatacol, from lines in which\n");
      printf("stringcol contains string. If stringcol==0, all lines are assumed to have xyz data\n");
      printf("\n");
      printf("Analyses atom coordinates for structure factor\n\n");
      printf("If string is not empty, takes the time step from any line not containing string but\n");
      printf("containing the string fs and analyses separately for each time. If string is empty\n");
      printf("any line with less than three arguments is interpreted to separate times.\n"); 
      printf("Outputs pair correlation function in cols 1 and 3 in neidrstat.*\n");
      printf("\n");
      printf("xsize, ysize and zsize give the box size\n");
      printf("Default: 1 Au 2 114.24 114.24 114.24\n");
      printf("\n");
      return(0);
   }

   lattice=FCC;

   iarg=1;

   /* Set some basic variables which may be modified by command line arguments */

   unitcella=4.08;
   unitcelln=0;
   dr=0.1;

   stringcol=1;
   sprintf(string,"XX"); 
   xdatacol=2;
   typedatacol=xdatacol+3;

   box[0]=114.24;
   box[1]=114.24;
   box[2]=114.24;

   /* Assume periodics */
   pbc[0]=pbc[1]=pbc[2]=1.0;

   /* Assume hard cutoff of 3.0 */
   cutRgen=3.0; cutSgen=3.0; 
   for (i=0;i<MAXTYPE;i++)  for (j=0;j<MAXTYPE;j++) {
     cutR[i][j]=cutRgen;
     cutS[i][j]=cutSgen;
   }
   binarycut=-2;
   binaryt1=1;
   binaryt2=2;
   SROr1=0.0; SROr2=0.0;
   cutmin=0.0; cutmax=0.0; // CB
   newbox[0]=1e30;
   printshort=-1.0;
   printbonds=False;
   printangles=False;
   writenbonds=False;
   printnot=-1;
   roughnessZ1=-1; roughnessZ2=-1; roughnessz0=0.0;
   onlyliq=False;
   noneidrstat=False;
   for (i=0;i<MAXTYPE;i++) printnot2[i]=-1;
      
   while (strncmp(argv[iarg],"-",1)==0) {
     if (strcmp(argv[iarg],"-d")==0) details=True; 
     else if (strcmp(argv[iarg],"-debug")==0) debug=True; 
     else if (strcmp(argv[iarg],"-maxdebug")==0) maxdebug=True; 
     else if (strcmp(argv[iarg],"-discardout")==0) discardout=True; 
     else if (strcmp(argv[iarg],"-enforcepbc")==0) enforcepbc=True; 
     else if (strcmp(argv[iarg],"-onlyliq")==0) onlyliq=True; 
     else if (strcmp(argv[iarg],"-noneidrstat")==0) noneidrstat=True; 
     else if (strcmp(argv[iarg],"-a")==0) sscanf(argv[++iarg],"%lg",&unitcella);
     else if (strcmp(argv[iarg],"-n")==0) sscanf(argv[++iarg],"%d",&unitcelln);
     else if (strcmp(argv[iarg],"-typecol")==0) sscanf(argv[++iarg],"%d",&typedatacol);

     else if (strcmp(argv[iarg],"-cut")==0) {
       sscanf(argv[++iarg],"%lg",&xx);
       sscanf(argv[++iarg],"%lg",&yy);
       cutRgen=xx; cutSgen=yy;
       for(i=0;i<MAXTYPE;i++) for(j=0;j<MAXTYPE;j++) {
	 cutR[i][j]=xx;
         cutS[i][j]=yy;
       }
       printf("Set all cutoffs between: R %g S %g. OVERRIDES ANY PREVIOUS -cutij!\n",xx,yy);
     }
     else if (strcmp(argv[iarg],"-cutij")==0) {
       sscanf(argv[++iarg],"%d",&i);
       sscanf(argv[++iarg],"%d",&j);
       sscanf(argv[++iarg],"%lg",&xx);
       sscanf(argv[++iarg],"%lg",&yy);
       cutR[i][j]=xx;
       cutS[i][j]=yy;
       cutR[j][i]=xx;
       cutS[j][i]=yy;
       printf("Set cutoff for pair %d %d and %d %d between: R %g S %g\n",i,j,j,i,xx,yy);
     }
     else if (strcmp(argv[iarg],"-dr")==0) {
       sscanf(argv[++iarg],"%lg",&dr);
       printf("Set drstat dr to %g\n",dr);
     }
     else if (strcmp(argv[iarg],"-printangles")==0) printangles=True;
     else if (strcmp(argv[iarg],"-writenbonds")==0) writenbonds=True;
     else if (strcmp(argv[iarg],"-printbonds")==0) printbonds=True;
     else if (strcmp(argv[iarg],"-printshort")==0) {
       sscanf(argv[++iarg],"%lg",&printshort);
       printf("Printing out bonds shorter than %g\n",printshort);
     }
     else if (strcmp(argv[iarg],"-printa")==0) {
       printeffectivea=True;
       sscanf(argv[++iarg],"%lg",&effectivearmax);
       sscanf(argv[++iarg],"%lg",&effectiveafact);
       printf("Printing effective a counted for bonds shorter than %g\n",effectivearmax);
       printf("using factor %g to scale average bond length into a\n",effectiveafact);
     }
     else if (strcmp(argv[iarg],"-printnot")==0) {
       sscanf(argv[++iarg],"%d",&printnot);
       printf("Printing out atoms not having %d bonds\n",printnot);
     }
     else if (strcmp(argv[iarg],"-printnot2")==0) {
       sscanf(argv[++iarg],"%d",&i);
       sscanf(argv[++iarg],"%d",&j);
       printnot2[i]=j;
       printf("Printing out atoms of type %d not having %d bonds\n",i,printnot2[i]);
     }
     else if (strcmp(argv[iarg],"-roughness")==0) {
       sscanf(argv[++iarg],"%d",&roughnessZ1);
       sscanf(argv[++iarg],"%d",&roughnessZ2);
       sscanf(argv[++iarg],"%lg",&roughnessz0);
       printf("Analysing surface roughness. Analysis assumes atoms with coordination in range %d %d\n",
	      roughnessZ1,roughnessZ2);
       printf("... which are above %g are surface atoms.\n",roughnessz0);
     }
     else if (strcmp(argv[iarg],"-binary")==0) {
       sscanf(argv[++iarg],"%d",&binarycut);
       sscanf(argv[++iarg],"%d",&binaryt1);
       sscanf(argv[++iarg],"%d",&binaryt2);
       printf("Getting binary neighbourhoood statistics for types %d and %d\n",binaryt1,binaryt2);
     }
     else if (strcmp(argv[iarg],"-SROpair")==0) {
       sscanf(argv[++iarg],"%lg",&SROr1);
       sscanf(argv[++iarg],"%lg",&SROr2);
       binarycut=0;
       sscanf(argv[++iarg],"%d",&binaryt1);
       sscanf(argv[++iarg],"%d",&binaryt2);
       printf("Getting SRO statistics between %g - %g\n",SROr1,SROr2);
       printf(" ... and also getting binary neighbourhoood statistics for types A=%d and B=%d\n",binaryt1,binaryt2);
     }
     else if (strcmp(argv[iarg],"-shellcut")==0) {
       sscanf(argv[++iarg],"%lg",&cutmin);
       sscanf(argv[++iarg],"%lg",&cutmax);
       printf("Getting neighbours in shell between %g - %g\n",cutmin,cutmax);
       printf("Remember to use -cut with values larger than cutmax!!!\n");
     }
     else if (strcmp(argv[iarg],"-pbc")==0) {
       sscanf(argv[++iarg],"%lg",&(pbc[0]));
       sscanf(argv[++iarg],"%lg",&(pbc[1]));
       sscanf(argv[++iarg],"%lg",&(pbc[2]));
       printf("Set periodics: x %g y %g z %g\n",pbc[0],pbc[1],pbc[2]);
     }
     else if (strcmp(argv[iarg],"-fcc")==0) lattice=FCC;
     else if (strcmp(argv[iarg],"-bcc")==0) lattice=BCC;
     else if (strcmp(argv[iarg],"-dia")==0) lattice=DIA;

     else {
       printf("bondinganalyze ERROR: Unknown option %s. Exiting...\n",argv[iarg]);
       exit(0);
     }
     
      iarg++;
   }

   sscanf(argv[iarg++],"%s",file);
   if (details) printf("file %s\n",file);
   if (strlen(file)==0) { printf("No filename given. Exiting\n"); return(0); }

   DEBUGSR("Args",iarg);
   if (iarg < argc) sscanf(argv[iarg++],"%d",&stringcol);
   if (iarg < argc) sscanf(argv[iarg++],"%s",string);
   if (iarg < argc) sscanf(argv[iarg++],"%d",&xdatacol);
   if (iarg < argc) sscanf(argv[iarg++],"%lg",&(box[0]));
   box[1]=box[0]; box[2]=box[0];
   if (iarg < argc) sscanf(argv[iarg++],"%lg",&(box[1]));
   if (iarg < argc) sscanf(argv[iarg++],"%lg",&(box[2]));

   if (unitcelln>0) { 
      unitcella=(box[0]+box[1]+box[2])/(3.0*unitcelln);
      printf("Number of unitcells %d => using a=%g\n",unitcelln,unitcella);
   }

   printf("Program arguments %d %s %d %lg %lg %lg\n",
	  stringcol,string,xdatacol,box[0],box[1],box[2]);
   fflush(NULL);

   usetime=True; 
   if (stringcol==0) {
     usetime=False;
     printf("Stringcol 0, assuming all atoms are of the same time\n");
   }
   for (i=0;i<MAXTYPE;i++) strncpy(name[i],"NotSet",8);
   ntime=0;
   casctime=0.0;

   cutSmax=0.0;
   printf("Using cutoffs of:\n   overall: R %g S %g\n",cutRgen,cutSgen);
   for (i=0;i<MAXTYPE;i++) for (j=0;j<MAXTYPE;j++) {
     if (cutS[i][j] > cutSmax) cutSmax=cutS[i][j];
     if (cutR[i][j] != cutRgen || cutS[i][j] != cutSgen) {
       printf("   pair %d - %d: R %g S %g\n",i,j,cutR[i][j],cutS[i][j]);
     }
   }

   if (printeffectivea && effectivearmax > cutSmax) {
     printf("ERROR: cutS too small for rmax calc. Do something !\n");
     exit(0);
   }
   if (xx<cutmax || yy<cutmin){
     printf("CB debug: Too small values for Tersoff cutoffs!\n");
     printf("CB debug: Set -cut R S to larger values than shellcutoffs.\n\n");
     exit(0);
   }

   MAXAT=200;

   i=sizeof(double)*MAXAT*3+sizeof(int)*MAXAT;
   i+=sizeof(int)*MAXAT;
   printf("Allocating %d bytes for %d atoms ... ",i,MAXAT);
   x0=(double *) malloc(sizeof(double)*3*MAXAT);
   in=(int *)  malloc(sizeof(int)*MAXAT);
   itype=(int *) malloc(sizeof(int)*MAXAT);
   if (onlyliq) wasliq=(int *) malloc(sizeof(int)*MAXAT);
   printf("... done\n");


   if (strcmp(file,"_")==0) fp=stdin;
   else fp=fopen(file,"r");

   if (fp==NULL) {
      printf("File %s open failed, exiting !\n",file);
      return;
   }

   DEBUGS("File opened\n");

   nat=0;
   nline=0;
   itypemin=MAXTYPE+1; itypemax=-1;
   firstdatafound=False;
   newbox[0]=box[0]; newbox[1]=box[1]; newbox[2]=box[2];
   ndiscard=ndiscardx=ndiscardy=ndiscardz=0;
   nenforcepbc=nenforcepbcx=nenforcepbcy=nenforcepbcz=0;
   ndiscardxfar=ndiscardyfar=ndiscardzfar=0;
   while(fgets(buf,159,fp) != NULL) {
      nline++;

      /* Fix atom type from first xyz atom line if not specified above */
      if (nline==3 && strcmp(string,"XX")==0) {
	 strncpy(string,buf,2);
	 printf("Picked atom name string %s from line %d\n",string,nline);
      }

      strcpy(buforig,buf);
      /* Split line into arguments */
      if (strlen(buf)>160) {
	printf("Warning: long line %d : %s\n",nline,buf);
	 continue;
      }
      narg=0; 
      iscoordline=False;
      if ((argbuf=(char *) strtok(buf," \n\t\0")) != NULL) {
	 narg++;
	 strcpy(arg[narg],argbuf);
	
	 if (narg==stringcol) {
	    /* Check whether line is coordinate line */
	    c1=arg[stringcol][0];
	    if (strlen(arg[stringcol])<3&&c1>='A'&&c1<='Z') {
	      iscoordline=True;
	    }
	 }

	 while((argbuf=(char *) strtok(NULL," \n\t\0")) != NULL) { 
	    narg++; strcpy(arg[narg],argbuf); 
	    if (narg==stringcol) {
	       if (strcmp(arg[stringcol],string)==0) {
		  iscoordline=True;
	       }
	       else {
		  c1=arg[stringcol][0];
		  if (c1>='A'&&c1<='Z'&&strlen(arg[stringcol])<3) {
		     printf("Impurity ignored: %s",buforig);
		  }
	       }
	    }
	 }
	 if (stringcol==0) iscoordline=True;
      }
      
      newtime=False;
      if (usetime) {
	 if (! iscoordline) {
	    DEBUGSS("Not coord line",buforig);
	    /* Check whether this is a new time */
	    if (strlen(string)==0 && narg < xdatacol+2) {
	       newtime=True;
	       prevtime=casctime;
	       ntime++; casctime+=1.0;
	       DEBUGSRR("New time found",ntime,casctime);
	    }
	    if (narg==1) {
	       i=sscanf(arg[1],"%d",&ii);
	       if (i==1) { 
		  natexpectedp=natexpected; natexpected=ii; 
		  if (nline==1) natexpectedp=ii;
		  printf("\nExpecting to read in %d atoms for next time %d\n",
			 natexpected,nline);
	       }
	    }
	    if (strlen(string) != 0 && strstr(buforig,"fs") != NULL) {
	       newtime=True;
	       prevtime=casctime;
	       ntime++; casctime+=1.0;
	       for (i=1;i<=narg;i++) {
		  if (strcmp(arg[i],"fs") == 0) sscanf(arg[i-1],"%lg",&casctime);
	       }
	    }

	    for (i=1;i<=narg;i++) {
	      /* Check whether line has box size as well */
	      if (strcmp(arg[i],"boxsize")==0) {
		if (newbox[0]!=1e30) {
		  box[0]=newbox[0];
		  box[1]=newbox[1];
		  box[2]=newbox[2];
		  if (unitcelln>0) {
		    unitcella=box[0]/unitcelln;
		    if (box[1]/unitcelln<unitcella) unitcella=box[1]/unitcelln;
		    if (box[2]/unitcelln<unitcella) unitcella=box[2]/unitcelln;
		  }
		}
		sscanf(arg[i+1],"%lg",&(newbox[0]));
		sscanf(arg[i+2],"%lg",&(newbox[1]));
		sscanf(arg[i+3],"%lg",&(newbox[2]));
		printf("Picked next cell size %lg %lg %lg from time line\n",
		       newbox[0],newbox[1],newbox[2]);
	      }
	    }

	    DEBUGSRRR("New time found",ntime,casctime,nline);
	   
	 }
      }

      if (!newtime) {
	 /* Read in atoms */
	 if (iscoordline && narg>=xdatacol+2) {
	    firstdatafound=True;
	    sscanf(arg[xdatacol+0],"%lg",&(x0[3*nat+0]));
	    sscanf(arg[xdatacol+1],"%lg",&(x0[3*nat+1]));
	    sscanf(arg[xdatacol+2],"%lg",&(x0[3*nat+2]));
	    sscanf(arg[typedatacol],"%d",&(itype[nat]));
	    if (itype[nat]<0) itype[nat]=-itype[nat];

	    in[nat]=nat;
	    if (narg>=xdatacol+4) sscanf(arg[xdatacol+4],"%d",&(in[nat]));
	    if (onlyliq) {
	      if (narg<typedatacol+2) { 
		printf("bondinganalyze onlyliq error: liquid column missing on atom line %d\n",nat);
		exit(0);
	      }
	      sscanf(arg[typedatacol+2],"%d",&(wasliq[nat]));
	    }

	    if (itype[nat]<itypemin) itypemin=itype[nat];
	    if (itype[nat]>itypemax) itypemax=itype[nat];
	    if (itype[nat]>=MAXTYPE) {
	      printf("MAXTYPE %d too small: increase to %d",MAXTYPE,itype[nat]+1);
	      exit(0);
	    }
	    if (xdatacol>0) strncpy(name[itype[nat]],arg[xdatacol-1],8);


	    if (enforcepbc) {
	      enforce=False;
	      if (pbc[0] == 1.0 && x0[3*nat+0] < -newbox[0]/2) { x0[3*nat+0]+=newbox[0]; enforce=True; nenforcepbcx++; }
	      if (pbc[0] == 1.0 && x0[3*nat+0] >= newbox[0]/2) { x0[3*nat+0]-=newbox[0]; enforce=True; nenforcepbcx++; }
	      if (pbc[1] == 1.0 && x0[3*nat+1] < -newbox[1]/2) { x0[3*nat+1]+=newbox[1]; enforce=True; nenforcepbcy++; }
	      if (pbc[1] == 1.0 && x0[3*nat+1] >= newbox[1]/2) { x0[3*nat+1]-=newbox[1]; enforce=True; nenforcepbcy++; }
	      if (pbc[2] == 1.0 && x0[3*nat+2] < -newbox[2]/2) { x0[3*nat+2]+=newbox[2]; enforce=True; nenforcepbcz++; }
	      if (pbc[2] == 1.0 && x0[3*nat+2] >= newbox[2]/2) { x0[3*nat+2]-=newbox[2]; enforce=True; nenforcepbcz++; }
	      if (enforce) nenforcepbc++;
	    }

	    if (discardout) {
	       disc=False;
	       xx=x0[3*nat+0]; yy=x0[3*nat+1]; zz=x0[3*nat+2];
	       if (xx<-newbox[0]/2 || xx>newbox[0]/2) { disc=True; ndiscardx++; }
	       if (yy<-newbox[1]/2 || yy>newbox[1]/2) { disc=True; ndiscardy++; }
	       if (zz<-newbox[2]/2 || zz>newbox[2]/2) { disc=True; ndiscardz++; }
	       if (xx<-newbox[0]/2-2*unitcella || 
		   xx>newbox[0]/2+2*unitcella) { disc=True; ndiscardxfar++; }
	       if (yy<-newbox[1]/2-2*unitcella || 
		   yy>newbox[1]/2+2*unitcella) { disc=True; ndiscardyfar++; }
	       if (zz<-newbox[2]/2-2*unitcella || 
		   zz>newbox[2]/2+2*unitcella) { disc=True; ndiscardzfar++; }
	       if (disc) {
		  nat--; 
		  ndiscard++;
	       }
	    }

	    nat++;
 	    if (nat>=MAXAT) { 
	       MAXAT*=2;
	       i=sizeof(double)*MAXAT*3+sizeof(int)*MAXAT;
	       i+=sizeof(int)*MAXAT;
	       printf("Reallocating %d bytes for %d atoms ... ",i,MAXAT);
	       x0=(double *) realloc(x0,3*sizeof(double)*MAXAT);
	       in=(int *)  realloc(in,sizeof(int)*MAXAT);
	       itype=(int *) realloc(itype,sizeof(int)*MAXAT);
	       if (onlyliq) wasliq=(int *) realloc(wasliq,sizeof(int)*MAXAT);
	       printf("... done\n");
	    }
	 }
      }

      if (firstdatafound && newtime) {
	 printf(" Read in %d atoms for time %d %lg fs\n",nat,ntime-1,prevtime);
	 if (discardout) {
	    printf("Discarded %d atoms outside box; x %d, y %d, z %d\n",
		   ndiscard,ndiscardx,ndiscardy,ndiscardz);
	    printf("Discarded atoms farther than 2 a0 out: x %d, y %d, z %d\n",
		   ndiscardxfar,ndiscardyfar,ndiscardzfar);
	 }
	 ndiscard=ndiscardx=ndiscardy=ndiscardz=0;
	 ndiscardxfar=ndiscardyfar=ndiscardzfar=0;
	 if (enforcepbc) {
	    printf("Enforced periodics on %d atoms outside box; x %d, y %d, z %d\n",
		   nenforcepbc,nenforcepbcx,nenforcepbcy,nenforcepbcz);
	 }
	 nenforcepbc=nenforcepbcx=nenforcepbcy=nenforcepbcz=0;

	 if (nat<3) { 
	    printf("No or too few atoms to handle, looking for next step %d\n",nat);
	    continue;
	 }
	 if (natexpected!=0) {
	    if (nat!=natexpectedp && !discardout) {
	       printf("Warning: Read in nat %d does not match predicted %d\n",
		      nat,natexpectedp);
	       /* continue; */
	    }
	 }

	 DEBUGSRRR("Last atom",x0[3*(nat-1)],x0[3*(nat-1)+1],x0[3*(nat-1)+2]);
	 printf("\n");
	 printf("--------------------------------------------------------------------\n");
	 printf("--------- Going into defect analysis at time %10g ------------\n",prevtime);
	 printf("--------------------------------------------------------------------\n");
	 printf("\n");
	 
	 bondinganalyze(x0,itype,wasliq,box,pbc,cutR,cutS,cutmin,cutmax,nat,unitcella,dr,prevtime,
			printshort,printbonds,printangles,writenbonds,onlyliq,noneidrstat,printnot,printnot2,
			roughnessZ1,roughnessZ2,roughnessz0,printeffectivea,effectivearmax,effectiveafact);

	 printf("\n");
	 printf("--------------------------------------------------------------------\n");
	 printf("---------- Done with defect analysis at time %10g ------------\n",prevtime);
	 printf("--------------------------------------------------------------------\n");
	 printf("\n");
	 
	 newtime=False;
	 firstdatafound=False;
	 printf(" Starting reading in atoms for next time %g\n",casctime);

	 /* Reset necessary variables */
	 nat=0;

      } 
	 
   } /* End of readin loop */

   printf("Read in all atoms\n");

   if (discardout) {
      printf("Discarded %d atoms outside box; x %d, y %d, z %d\n",
	     ndiscard,ndiscardx,ndiscardy,ndiscardz);
      printf("Discarded atoms farther than 2 a0 out: x %d, y %d, z %d\n",
	     ndiscardxfar,ndiscardyfar,ndiscardzfar);
   }
   ndiscard=ndiscardx=ndiscardy=ndiscardz=0;
   if (enforcepbc) {
     printf("Enforced periodics on %d atoms outside box; x %d, y %d, z %d\n",
	    nenforcepbc,nenforcepbcx,nenforcepbcy,nenforcepbcz);
   }
   nenforcepbc=nenforcepbcx=nenforcepbcy=nenforcepbcz=0;

   if (nat<3) { 
      printf("No or too few atoms to handle, exiting program %d\n",nat);
      return;
   }

   if (natexpected!=0) {
      if (nat!=natexpectedp && !discardout) {
	 printf("Warning: Read in nat %d does not match predicted %d\n",
		nat,natexpectedp);
      }
   }
   
   if (!usetime) casctime=1.0; 

   if (newbox[0]!=1e30) {
      box[0]=newbox[0];
      box[1]=newbox[1];
      box[2]=newbox[2];
      if (unitcelln>0) {
	 unitcella=box[0]/unitcelln;
	 if (box[1]/unitcelln<unitcella) unitcella=box[1]/unitcelln;
	 if (box[2]/unitcelln<unitcella) unitcella=box[2]/unitcelln;
      }
   }

   printf("\n");
   printf("--------------------------------------------------------------------\n");
   printf("--------- Going into defect analysis at time %10g ------------\n",casctime);
   printf("--------------------------------------------------------------------\n");
   printf("\n");
   
   bondinganalyze(x0,itype,wasliq,box,pbc,cutR,cutS,cutmin,cutmax,nat,unitcella,dr,casctime,
		  printshort,printbonds,printangles,writenbonds,onlyliq,noneidrstat,printnot,printnot2,
		  roughnessZ1,roughnessZ2,roughnessz0,printeffectivea,effectivearmax,effectiveafact);
   
   printf("\n");
   printf("--------------------------------------------------------------------\n");
   printf("---------- Done with defect analysis at time %10g ------------\n",casctime);
   printf("--------------------------------------------------------------------\n");
   printf("\n");
	 


} /* End of program */


/*

  Actual analysis of data. 
  
*/

#define MAXNGBRS 2000



#define round(x) ((int) (x>0 ? x+0.5 : x-0.5))

void bondinganalyze(
		    real *x0,
		    int *itype,
		    int *wasliq,
		    real *box,
		    real *pbc,
		    real cutR[MAXTYPE][MAXTYPE],
		    real cutS[MAXTYPE][MAXTYPE],
		    real cutmin,
		    real cutmax,
		    int nat,
		    real a,
		    real dr,
		    real time,
		    real printshort,
		    logical printbonds,
		    logical printangles,
		    logical writenbonds,
		    logical onlyliq,
		    logical noneidrstat,
		    int printnot,
		    int printnot2[MAXTYPE],
		    int roughnessZ1,
		    int roughnessZ2,
		    real roughnessz0,
		    real printeffectivea,
		    real effectivearmax,
		    real effectiveafact		    
		    )

/*
   Remember that due to the nature of arrays in C and Fortran,
   any index in C is one less than that in Fortran !
*/

{

   static int analtime=0;   /* Counter of how many time we've been here */

   /* Array which returns neighbour info */
   struct ngbr ngbrs[MAXNGBRS];
   int nngbrstat[MAXNGBRS];
   int ntype[MAXTYPE];
   int typenngbrstat[MAXTYPE][MAXNGBRS];
   real typenngbrave[MAXTYPE];
   real typenngbraver[MAXTYPE];
   int typenngbravern[MAXTYPE];
   int typenngbrstattypestat[MAXTYPE][MAXNGBRS][MAXTYPE];

   int numberoftypetypepairs[MAXTYPE][MAXNGBRS][MAXNGBRS];

   int oneatomngbrtypestat[MAXTYPE];

   int *nngbrarray;
   int roughnessN;
   real roughnesszsum;
   real roughnesszmean;
   real sqsum,abssum;

   char file[160],bdfile[160];

   logical print;
   int jj,kk,nngot,nrave,it,jt,nnsum=0; // CB
   real r,boxmin,rave,r110ave,r1,r2,r3,theta,Vdr,natideal;
   
   /* Neighbour parameters */
   static real rnmax;
   int npairs;         /* Total number of neighbours */
   real nbinaryt1,nbinaryt2;

   int SROnBnnthis,SROnBAnnthis,SROnBnn,SROnBAnn,SROnB;
   double SROnBAsum,SROfBA1,SROfBA2,SROBA1,SROBA2,cA;
   int SROnAnnthis,SROnABnnthis,SROnAnn,SROnABnn,SROnA;
   double SROnABsum,SROfAB1,SROfAB2,SROAB1,SROAB2,cB;

   int drstat[1000000],ndr;

   /* Cutoff function params */
   real nngbr,nngbrmax,fc;

   int i,i3,j,k;

   static logical firsttime=True;

   FILE *fp,*fpbd,*fpbd2,*fpa;

   if (roughnessZ2>0) {
     printf("Allocating neighbour number array for roughness analysis");
     nngbrarray=malloc(sizeof(real)*(nat+1));
     roughnessN=0;
     roughnesszsum=0.0;
   }

   rnmax=cutS[0][0];
   for (i=0;i<MAXTYPE;i++) for (j=0;j<MAXTYPE;j++) {
     if (cutS[i][j] > rnmax) rnmax=cutS[i][j];
   }
   printf("Starting bondinganalyze at time %g for %d atoms. rcut %g\n",time,nat,rnmax);


   if (SROr2 > 0.0) {
     /* Check that Tersoff cutoff is long enough */
     if (SROr2 > cutR[binaryt1][binaryt1] || SROr2 > cutR[binaryt1][binaryt2] || SROr2 > cutR[binaryt2][binaryt2]) {
       printf("Error: SRO analysis requires lower cutoffs > SROR2!\n");
       printf("SROR2 %g  cutR[binaryt1][binaryt1] %g cutR[binaryt1][binaryt2] %g cutR[binaryt2][binaryt2] %g\n",
              SROr2,cutR[binaryt1][binaryt1],cutR[binaryt1][binaryt2],cutR[binaryt2][binaryt2]);
       exit(0);
     }
   }

   if (firsttime) {

   }

   DEBUGSRRRR("bondinganalyze",nat,a,rnmax,time);

   /* Initialize neighbour list calc. */
   i=-1;
   getngbrs(i,x0,box,pbc,nat,rnmax,ngbrs,&nngot);
   npairs=0;

   /* Initialize statistics */
   for(i=0;i<MAXNGBRS;i++) {
     nngbrstat[i]=0;
   }
   for(i=0;i<MAXTYPE;i++) {
     oneatomngbrtypestat[i]=0;
     typenngbrave[i]=0.0;
     typenngbraver[i]=0.0;
     typenngbravern[i]=0;
     ntype[i]=0;
     for(j=0;j<MAXNGBRS;j++) {
       typenngbrstat[i][j]=0;
       for(k=0;k<MAXTYPE;k++) {
	 typenngbrstattypestat[i][j][k]=0;
       }
       for(k=0;k<MAXNGBRS;k++) {
	 numberoftypetypepairs[i][j][k]=0;
       }
     }
   }
   boxmin=box[0]; if (box[1]<boxmin) boxmin=box[1];
   if (box[2]<boxmin) boxmin=box[2];
   ndr=boxmin/dr+1; if (ndr>1000000) { printf("Increase drstat size %d",ndr); exit(0); }
   for (i=0;i<ndr;i++) drstat[i]=0;
   nngbrmax=0.0;

   if (printnot>=0) {
     sprintf(bdfile,"bonddefects.out");
     if (firsttime) fpbd=fopen(bdfile,"w");
     else fpbd=fopen(bdfile,"a");
   }

   print=False;
   for (i=0;i<MAXTYPE;i++) if (printnot2[i]>=0) print=True;
   if (print>=0) {
     sprintf(bdfile,"bonddefects2.out");
     if (firsttime) fpbd2=fopen(bdfile,"w");
     else fpbd2=fopen(bdfile,"a");
   }

   if (printeffectivea) {
     if (firsttime) fpa=fopen("effectivea.xyz","w");
     fprintf(fpa," %d\n",nat);
     fprintf(fpa," Effective a values at time %g fs obtained from bonds < %g scaled by %g\n",
	     time,effectivearmax,effectiveafact);
   }

   SROnBnn=0; SROnBAnn=0; SROnB=0; SROnBAsum=0.0; 
   SROnAnn=0; SROnABnn=0; SROnA=0; SROnABsum=0.0; 

   /* First loop over all atoms */
   for (i=0;i<nat;i++) {
      if (i>0&&i%10000==0) fprintf(stderr,"%1d",(i%100000)/10000); 
      if (onlyliq) { if (!wasliq[i]) continue; }
      getngbrs(i,x0,box,pbc,nat,rnmax,ngbrs,&nngot);
      npairs+=nngot;
      DEBUGSRR("Got ngbrs",i,nngot);



      if (printshort>=0) {
	for (jj=0;jj<nngot;jj++) {
	  r=ngbrs[jj].r;
	  j=ngbrs[jj].i;
	  if (r<=printshort) {
	    printf("Short bond: %2s %10g %10g %10g %d %8d   %.10g\n",
		   name[itype[i]],x0[i*3],x0[i*3+1],x0[i*3+2],itype[i],i,r);
	    printf("            %2s %10g %10g %10g %d %8d\n",
		   name[itype[j]],x0[j*3],x0[j*3+1],x0[j*3+2],itype[j],j);
	  }
	}
      }		
      if (printangles) {
	for (jj=0;jj<nngot;jj++) {
	  for (kk=jj+1;kk<nngot;kk++) {
	    j=ngbrs[jj].i;
	    k=ngbrs[kk].i;
	    theta=getangle(i,j,k,x0,box,pbc);
	    printf("Theta = %10.4f for atoms %6d %6d %6d of types %2d %2d %2d\n",theta,i,j,k,itype[i],itype[j],itype[k]);
	    if (theta<1e-6) {
	      printf("            %2s %10g %10g %10g %d %8d\n",
		   name[itype[i]],x0[i*3],x0[i*3+1],x0[i*3+2],itype[i],i);      
	      printf("            %2s %10g %10g %10g %d %8d\n",
		   name[itype[j]],x0[j*3],x0[j*3+1],x0[j*3+2],itype[j],j);
	      printf("            %2s %10g %10g %10g %d %8d\n",
		   name[itype[k]],x0[k*3],x0[k*3+1],x0[k*3+2],itype[k],k);
	    }
	  }
	}
      }
      if (printeffectivea) {
	rave=0.0; nrave=0;
	r110ave=0.0; 
 	for (jj=0;jj<nngot;jj++) {
	  j=ngbrs[jj].i;
	  if (ngbrs[jj].r < effectivearmax) { 
	    rave+=ngbrs[jj].r; nrave++; 
	    r1=fabs(x0[j*3+0]-x0[i*3+0]);
	    if (r1+r1>box[0] && pbc[0] == 1.0) r1=box[0]-r1;
	    r2=fabs(x0[j*3+1]-x0[i*3+1]);
	    if (r2+r2>box[1] && pbc[1] == 1.0) r2=box[1]-r2;
	    r110ave+=sqrt(r1*r1+r2*r2);
	  }
	}
	if (nrave>0) { rave/=nrave; r110ave/=nrave; }
	fprintf(fpa,"%2s %10g %10g %10g %g %d %d %g %g\n",name[itype[i]],x0[i*3],x0[i*3+1],x0[i*3+2],
		rave*effectiveafact,nrave,itype[i],r110ave/sqrt(2)*sqrt(3)*effectiveafact,(rave>0.0?1.0/pow(rave*effectiveafact,3.0):0.0));
      }

      /* Get number of neighbours using Tersoff cutoff function */
      nngbr=0.0;
      nbinaryt1=0; nbinaryt2=0;
      SROnBnnthis=0;  SROnAnnthis=0;
      SROnBAnnthis=0; SROnABnnthis=0;
      for (jj=0;jj<nngot;jj++) {
	r=ngbrs[jj].r;
	it=itype[i];
	jt=itype[ngbrs[jj].i];
	if (r <= cutR[it][jt] ) fc=1.0; 
	else if (r>=cutS[it][jt]) fc=0.0;
	else fc=0.5+0.5*cos(pi*(r-cutR[it][jt])/(cutS[it][jt]-cutR[it][jt]));
	if (cutmax>0.0) { // CB
	   if (r > cutmin && r <= cutmax) {// CB
	     fc=1.0;// CB
	   }// CB
	   else fc=0.0; // CB
	}// CB
	ngbrs[jj].s=fc;
	nngbr+=fc;
	if (binarycut>-2 && fc>0.0) {
	  if (itype[ngbrs[jj].i]==binaryt1) nbinaryt1+=fc;
	  if (itype[ngbrs[jj].i]==binaryt2) nbinaryt2+=fc;

          if (SROr2>0.0) {
            if (r > SROr1 && r <= SROr2) {
              /* This pair is within the desired SRO range, do analysis. Discount non-A-B atoms */
              if (it==binaryt2) {
                if (jt==binaryt1 || jt==binaryt2) SROnBnnthis++; 
                if (jt==binaryt1) SROnBAnnthis++;
              }
              if (it==binaryt1) {
                if (jt==binaryt2 || jt==binaryt1) SROnAnnthis++; 
                if (jt==binaryt2) SROnABnnthis++;
              }
              DEBUGSRRRR("SRO",it,jt,SROnBnnthis,SROnBAnnthis);
            }
          }
	}

	if (printbonds && fc > 0.0) {
	  j=ngbrs[jj].i;
	  r1=x0[j*3+0]-x0[i*3+0]; r2=x0[j*3+1]-x0[i*3+1]; r3=x0[j*3+2]-x0[i*3+2]; 
	  printf("Tersoff bond %s - %s r %g vec %g %g %g\n",name[itype[i]],name[itype[j]],r,r1,r2,r3);

	}

      } /* End of loop over neighbours to atom i */

      if (nngbrmax<nngbr) nngbrmax=nngbr;

      if (printnot2[itype[i]]>=0 && (int)(nngbr+0.5)!=printnot2[itype[i]]) {
	/* Print out bonddefects */
	fprintf(fpbd2,"%lg %lg %lg %lg %d %g",x0[i*3],x0[i*3+1],x0[i*3+2],
		time,itype[i],nngbr);
	if (binarycut>-2) fprintf(fpbd2," %g %g",nbinaryt1,nbinaryt2);
	fprintf(fpbd2,"\n");
      }

      if (roughnessZ2>0) {
	nngbrarray[i]=nngbr;
	if (nngbr>=roughnessZ1 && nngbr <=roughnessZ2 && x0[i*3+2] > roughnessz0) {
	  /* This is a surface atom by this analysis criteria, sum up quantities */
	  roughnessN++;
	  roughnesszsum+=x0[i*3+2]; 
	} 
      }

      /* Do statistics of number of neighbours, rounding off to nearest integer */
      nngbrstat[(int)(nngbr+0.5)]++;
      ntype[itype[i]]++;
      typenngbrstat[itype[i]][(int)(nngbr+0.5)]++;
      typenngbrave[itype[i]]+=nngbr;

      for (j=itypemin;j<=itypemax;j++) oneatomngbrtypestat[j]=0;
      /* Do various statistics from half bond length*/
      for (jj=0;jj<nngot;jj++) {
	if (ngbrs[jj].s<0.5) continue; 

	oneatomngbrtypestat[itype[ngbrs[jj].i]]++;
	typenngbrstattypestat[itype[i]][(int)(nngbr+0.5)][itype[ngbrs[jj].i]]++;

	if (ngbrs[jj].r <= boxmin) drstat[(int)(ngbrs[jj].r/dr+0.5)]++;
	typenngbraver[itype[i]]+=ngbrs[jj].r;
	typenngbravern[itype[i]]++;
      }

      if (binarycut>-2) {
	if (itype[i]==binaryt1||itype[i]==binaryt2)
	  numberoftypetypepairs[itype[i]][oneatomngbrtypestat[binaryt1]]
	    [oneatomngbrtypestat[binaryt2]]++;
        
        if (SROr2>0.0) {
          if (itype[i] == binaryt2) {
            SROnBnn+=SROnBnnthis;
            SROnBAnn+=SROnBAnnthis;
            if (SROnBnnthis>0.0) SROnBAsum+=1.0*SROnBAnnthis/SROnBnnthis;
            SROnB++;
          }
          if (itype[i] == binaryt1) {
            SROnAnn+=SROnAnnthis;
            SROnABnn+=SROnABnnthis;
            if (SROnAnnthis>0.0) SROnABsum+=1.0*SROnABnnthis/SROnAnnthis;
            SROnA++;
          }

          if (details) printf("Atom %d has SRO params nBA/nB %g nAB/nA %g\n",i,
                              1.0*SROnBAnnthis/(SROnBnnthis>0?SROnBnnthis:1),
                              1.0*SROnABnnthis/(SROnAnnthis>0?SROnAnnthis:1));
        }
      }

      if (details) {
	printf("Atom %s %.4g %.4g %.4g %d %d has %g neighbours: ",
	       name[itype[i]],x0[i*3],x0[i*3+1],x0[i*3+2],itype[i],in[i],nngbr);
	for(k=itypemin;k<=itypemax;k++) {
	  printf("%d %s ",oneatomngbrtypestat[k],name[k]); 
	  if (k<itypemax) printf("and ");  
	} 
	printf("\n");
	for (jj=0;jj<nngot;jj++) {
	  j=ngbrs[jj].i; 
	  printf("   %s %d %g %g %g %d %d r %g fc %g\n",
		 name[itype[j]],j,x0[j*3],x0[j*3+1],x0[j*3+2],itype[j],in[j],ngbrs[jj].r,ngbrs[jj].s);
	}		
      }
      nnsum+=nngbr;
      

   } /* End of loop over atoms */
   printf("CB: A total of %d neighbours was found in shell %f -> %f A.\n",nnsum,cutmin,cutmax);
   printf("CB: These are sum_j=1^v B(i,j), where v is total amount of vacancies \n    and B is the number of vacancies in the shell.\n");
   printf("CB: Divide with v to get <N(i)> [Current et al. Phil. Mag. A 43 (1981) 130]\n");
   printf("\n");

   printf("Number of neighbours stat: ");
   for(jj=MAXNGBRS-2;jj>0;jj--) if (nngbrstat[jj]>0) break;
   for(i=0;i<=jj+1;i++) printf("  %d: %d",i,nngbrstat[i]);
   printf("\n");

   printf("\nTypewise number of neighbours stat: \n");
   for (j=itypemin;j<=itypemax;j++) {
     printf("%s: ",name[j]);
     for(jj=MAXNGBRS-2;jj>0;jj--) if (typenngbrstat[j][jj]>0) break;
     for(i=0;i<=jj+1;i++) printf("  %d: %d",i,typenngbrstat[j][i]);
     printf("\n");
   }
   printf("\nAverage coordination numbers and bond lengths (using Tersoff cutoff function): \n");
   for (j=itypemin;j<=itypemax;j++) {     
     if (ntype[j]>0&&typenngbravern[j]>0) printf("%s: <Z> %.5f <r> %.6f for %d atoms and %d bonds\n",name[j],typenngbrave[j]/ntype[j],typenngbraver[j]/typenngbravern[j],ntype[j],typenngbravern[j]);
   }

   fflush(NULL);

   if (roughnessZ2>0) { 
     /* Do final roughness analysis */
     printf("Roughness analysis detected %d surface atoms\n",roughnessN);
     if (roughnessN>1) {
       roughnesszmean=roughnesszsum/roughnessN;
       sqsum=0.0;
       abssum=0.0;
       for(i=0;i<nat;i++) {
	 nngbr=nngbrarray[i];
	 if (nngbr>=roughnessZ1 && nngbr <= roughnessZ2 && x0[i*3+2] > roughnessz0) {
	   /* This is a surface atom by this analysis criteria, sum up quantities */
	   sqsum+=(x0[i*3+2]-roughnesszmean)*(x0[i*3+2]-roughnesszmean);
	   abssum+=abs(x0[i*3+2]-roughnesszmean);
	 } 
       } /* End of loop over atoms */
       sqsum=sqsum/(roughnessN-1);
       abssum=abssum/roughnessN;

       printf("Surface roughness at time %g average z %g standard deviation %g abs. dev. %g for %d atoms\n",
	      time,roughnesszmean,sqrt(sqsum),abssum,roughnessN);
     }

     /* free(nngbrarray); */
   }


   printf("\nTypewise number of neighbours type stat: \n");
   for (j=itypemin;j<=itypemax;j++) {
     for(jj=MAXNGBRS-2;jj>0;jj--) if (typenngbrstat[j][jj]>0) break;
     for(i=0;i<=jj;i++) {
       if (typenngbrstat[j][i]==0) continue;
       printf("%2d-fold bonded %5d %2s atoms have: ",i,typenngbrstat[j][i],name[j]);
       for(k=itypemin;k<=itypemax;k++) {
	 printf("%5d %2s ",typenngbrstattypestat[j][i][k],name[k]);
	 if (k<itypemax) printf("and ");  
       }
       printf("bonds\n");
     }
     printf("\n");
   }


   if (writenbonds) {
     printf("\nWriting out number of atoms with n bonds to files nbonds* \n");
     for (j=itypemin;j<=itypemax;j++) {
       for(i=0;i<=20;i++) {
	 sprintf(file,"nbonds%s_%d",name[j],i);
	 if (firsttime) fp=fopen(file,"w");
	 else fp=fopen(file,"a");
	 fprintf(fp,"%g %d\n",time, typenngbrstat[j][i]);
	 fclose(fp);
       }
     }
   }

   if (binarycut>-2) {
     printf("Binary neighbourhood type statistics for types %s %s with n>%d:\n",
	    name[binaryt1],name[binaryt2],binarycut);
     for (i=itypemin;i<=itypemax;i++) {
       if (i==binaryt1||i==binaryt2) {
	 for(j=0;j<=nngbrmax;j++) for(k=0;k<=nngbrmax;k++) {
	   if (numberoftypetypepairs[i][j][k]>binarycut) {
	     printf("Number of %s with %s%d %s%d neighbourhoods: %d\n",
		    name[i],name[binaryt1],j,name[binaryt2],k,
		    numberoftypetypepairs[i][j][k]);
	   }
	 }
       }
     }

     if (SROr2>0.0) {
       /* 
          Note that there is a subtle difference in how thge SRO is can be defined
          depending on whether one calculates the average <nBA> and <nB> for all atoms,
          or the average number of neighbours per atom <nBA/nB>. In practise this
          hardly matters in a crystal but to be on the safe side this code calculates 
          both.
       */
       
       cA=1.0*SROnA/(SROnA+SROnB); 
       cB=1.0*SROnB/(SROnA+SROnB); 
       SROBA1=0.0; SROBA2=0.0; if (SROnB > 0 && cA>0.0) { 
         SROfBA1=(1.0*SROnBAnn/SROnB)/(1.0*SROnBnn/SROnB); 
         SROfBA2=(SROnBAsum/SROnB);
         SROBA1 = 1.0-SROfBA1/cA; 
         SROBA2 = 1.0-SROfBA2/cA; 
       }
       SROAB1=0.0; SROAB2=0.0; if (SROnA > 0 && cB>0.0) { 
         SROfAB1=(1.0*SROnABnn/SROnA)/(1.0*SROnAnn/SROnA); 
         SROfAB2=(SROnABsum/SROnA);
         SROAB1 = 1.0-SROfAB1/cB; 
         SROAB2 = 1.0-SROfAB2/cB; 
       }
       DEBUGSRRRR("SRO",SROnBAnn,SROnB,SROnBnn,cA);
       if (details) printf("SRO SROnBAnn %d SROnB %d SROnBnn %d\n",SROnBAnn,SROnB,SROnBnn);
       if (details) printf("SRO SROnABnn %d SROnA %d SROnAnn %d\n",SROnABnn,SROnA,SROnAnn);
       printf(" Short-range order SRO BA %10.6f using cA %9.6f for <nBA>/<nB> %g\n",SROBA1,cA,SROfBA1);
       printf("(short-range order SRO BA %10.6f using cA %9.6f for <nBA/nB>   %g)\n",SROBA2,cA,SROfBA2);
       printf(" Short-range order SRO AB %10.6f using cB %9.6f for <nAB>/<nA> %g\n",SROAB1,cB,SROfAB1);
       printf("(short-range order SRO AB %10.6f using cB %9.6f for <nAB/nA>   %g)\n",SROAB2,cB,SROfAB2);
     }

   }

   if (! noneidrstat) {
     /* Finally output of data files */
     sprintf(file,"neidrstat.%g",time);
     fp=fopen(file,"w");
     printf("Atom distance and pair correlation function written to %s\n",file);
     for (i=0;i<ndr;i++) {
       /* 
	  Calculate normalization volume as exact thickness of shell, not 4 pi r^2 dr 
	  assuming density is N/box volume. See Allen-Tildesley p. 184
       */
       r1=(i-0.5)*dr;
       r2=(i+0.5)*dr;
       Vdr=4*pi/3*(r2*r2*r2-r1*r1*r1);
       natideal=nat/(box[0]*box[1]*box[2])*Vdr;
       if (r2 > rnmax) break;
       fprintf(fp,"%g %d %g %g\n",i*dr,drstat[i],drstat[i]/(natideal*nat),Vdr);
     }
     fflush(fp);
     fclose(fp);
   }

   if (printnot>=0) fclose(fpbd);

   analtime++;

   fflush(stdout); fflush(stderr); 
   firsttime=False;

}



/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*
*
* Linked list subroutine written by Kai Nordlund. 
*
* Written to be as self-contained as possible to make it easy to move
* around if needed. Gives neighbours for one atom only ! When called with
* atom index -1, initializes link cell matrices
*
* Input:
*         iat        Atom index of atom for which neighbours should be got
*         x0(i)       Atom coordinates 
*         box(3)     Cell size in A (or any units)
*         pbc(3)     Periodic boundary conditions (1.0 if true)
*         N          Number of atoms to deal with
*         rnmax      Neighbour cutoff radius in A        
*
* Output: ngbrs[]    Neighbour data
*
*/


void getngbrs(int iat,
	      real *x0,
	      real *box,
	      real *pbc,
	      int N,
	      real Rneicut,
	      struct ngbr *ngbrs,
	      int *nngot)

{

/*
*  Size of one cell is set to exactly Rneicut.
*  Cells are distributed so that cell 0 is between 0 and Rneicut.
*  However, cell size may vary between calls of the subroutine.
*
*/

#define nind(a,b,c) (((a-Nixmin)*ync+(b-Niymin))*znc+(c-Nizmin))

/*
    Number of cells in each direction
*/
   static int xnc,ync,znc;

/*
   Linked list cell arrays, as in Allen-Tildesley
*/

   static int *list,*head; 
   static int Nallocated;
   int icell;

/* Other variables */

   int Nix,Niy,Niz;
   static int Nixmin,Niymin,Nizmin;
   static int Nixmax,Niymax,Nizmax;
   static int Nixmax0,Niymax0,Nizmax0;
   static int curngbrsize;

   real rsq;
   int i,ii,i3,j,jj,j3,nn,nj;
   real r1,r2,r3;
   int ix,iy,iz;
   int dix,diy,diz;
   int inx,iny,inz;

   static real cellsize[3];
   real Rneicutsq;
   
   static logical firsttime=True;

   int ngbrcmp();

   MAXDEBUGSRRR("Starting getngbrs",Rneicut,iat,N); 

   Rneicutsq=Rneicut*Rneicut;

   Nixmax=(int)(box[0]/2.0/Rneicut-0.5);
   Nixmin=-Nixmax;
   Niymax=(int)(box[1]/2.0/Rneicut-0.5);
   Niymin=-Niymax;
   Nizmax=(int)(box[2]/2.0/Rneicut-0.5);
   Nizmin=-Nizmax;
   xnc=2*Nixmax+1; 
   ync=2*Niymax+1; 
   znc=2*Nizmax+1; 

   MAXDEBUGSRRR(" ",xnc,ync,znc);

   cellsize[0]=box[0]/(2*Nixmax+1.0);
   cellsize[1]=box[1]/(2*Niymax+1.0);
   cellsize[2]=box[2]/(2*Nizmax+1.0);

   if (firsttime) {
      Nixmax0=Nixmax;
      Niymax0=Niymax;
      Nizmax0=Nizmax;
      printf("\nFirst time in getngbr: allocating %lu bytes for linkcell arrays\n",
	     sizeof(int)*xnc*ync*znc+sizeof(int)*(N+1));

      printf("Cell ranges %d - %d, %d - %d, %d - %d\n",
	     Nixmin,Nixmax,Niymin,Niymax,Nizmin,Nizmax);
      printf("Cell size %g %g %g\n\n",cellsize[0],cellsize[1],cellsize[2]);
      fflush(stdout);

      head=(int *) malloc(sizeof(int)*xnc*ync*znc);
      list=(int *) malloc(sizeof(int)*(N+1));
      Nallocated=N;
   }

   if (Nixmax>Nixmax0 || Niymax>Niymax0 || Nizmax>Nizmax0) {
      
      printf("\ngetngbr: cell size has changed, increasing arrays\n");
      printf("Cell ranges %d - %d, %d - %d, %d - %d\n",
	     Nixmin,Nixmax,Niymin,Niymax,Nizmin,Nizmax);
      printf("Cell size %g %g %g\n\n",cellsize[0],cellsize[1],cellsize[2]);
      fflush(stdout);
      Nixmax0=Nixmax;
      Niymax0=Niymax;
      Nizmax0=Nizmax;
      head=(int *) realloc(head,sizeof(int)*xnc*ync*znc);
      list=(int *) realloc(list,sizeof(int)*(N+1));
   }

   if (N>Nallocated) {
      
      printf("\ngetngbr: number of atoms has changed, reallocating list array\n");
      list=(int *) realloc(list,sizeof(int)*(N+1));
      Nallocated=N;
   }
   


   if (iat==-1) {
     /*
	Arrange atoms in cells
     */

      /* Zero head array with -1:s since 0 is a valid atom index ! */
      for (ix=Nixmin;ix<=Nixmax;ix++) {
	 for (iy=Niymin;iy<=Niymax;iy++) {
	    for (iz=Nizmin;iz<=Nizmax;iz++) {
	       head[nind(ix,iy,iz)]=-1;
	    }
	 }
      }
      
      for(i=0;i<N;i++) {
         i3=3*i;
         ix=round(x0[i3]/cellsize[0]);
	 if (ix < Nixmin || ix > Nixmax) {
            if (pbc[0] == 1.0) {
               printf("bondinganalyze.c: This is impossible, atom %d x outside cell %d %g xsize %g\n",i,ix,x0[i3],box[0]);
	    }
            else { 
	      if (ix < Nixmin) ix=Nixmin;
	      if (ix > Nixmax) ix=Nixmax;
	    } 
         }

         iy=round(x0[i3+1]/cellsize[1]);
         if (iy < Niymin || iy > Niymax) {
            if (pbc[1] == 1.0) {
               printf("bondinganalyze.c: This is impossible, atom %d y outside cell %d %g ysize %g\n",i,iy,x0[i3+1],box[1]);
	    }
            else { 
	      if (iy < Niymin) iy=Niymin;
	      if (iy > Niymax) iy=Niymax;
            } 
         }

         iz=round(x0[i3+2]/cellsize[2]);
         if (iz < Nizmin || iz > Nizmax) {
            if (pbc[2] == 1.0) {
	      printf("bondinganalyze.c: This is impossible, atom %d z outside cell %d %g zsize %g\n",i,iz,x0[i3+2],box[2]); 
	    }
            else { 
	      if (iz < Nizmin) iz=Nizmin;
	      if (iz > Nizmax) iz=Nizmax;
            } 
         }
         
	 icell=nind(ix,iy,iz);
	 /* printf("i %d ix iy iz %d %d %d icell %d head[icell] %d\n",i,ix,iy,iz,icell,head[icell]); */
	 list[i]=head[icell];
	 head[icell]=i;

      }

      MAXDEBUGSR("Organized atoms into cells",i);

      firsttime=False;
      return;

   } /* End of cell initialization */


   /* Find neighbours of atom iat */

   for(i=0;i<MAXNGBRS;i++) {
      ngbrs[i].r=0;
      ngbrs[i].i=-1;
   }
   i=iat;
   i3=3*i;

   ix=round(x0[i3+0]/cellsize[0]);
   if (ix < Nixmin || ix > Nixmax) {
     if (pbc[0] == 1.0) {
       printf("bondinganalyze.c: This is impossible, atom %d x outside cell %d %g xsize %g\n",i,ix,x0[i3],box[0]);
     }
     else {
       if (ix < Nixmin) ix=Nixmin;
       if (ix > Nixmax) ix=Nixmax;
     }
   }

   iy=round(x0[i3+1]/cellsize[1]);
   if (iy < Niymin || iy > Niymax) {
     if (pbc[1] == 1.0) {
       printf("bondinganalyze.c: This is impossible, atom %d y outside cell %d %g ysize %g\n",i,iy,x0[i3+1],box[1]);
     }
     else {
       if (iy < Niymin) iy=Niymin;
       if (iy > Niymax) iy=Niymax;
     }
   }
   iz=round(x0[i3+2]/cellsize[2]);
   if (iz < Nizmin || iz > Nizmax) {
     if (pbc[2] == 1.0) {
       printf("bondinganalyze.c: This is impossible, atom %d z outside cell %d %g zsize %g\n",i,iz,x0[i3+2],box[2]); 
     }
     else {
       if (iz < Nizmin) iz=Nizmin;
       if (iz > Nizmax) iz=Nizmax;
     }
   }

   MAXDEBUGSRRRRR("Finding neighbours",iat,i3,ix,iy,iz);

   nn=0;
   for(dix=-1;dix<=1;dix++) {
      inx=ix+dix;
      if (inx < Nixmin && pbc[0] == 1.0) inx=Nixmax;
      if (inx > Nixmax && pbc[0] == 1.0) inx=Nixmin;
      if (inx < Nixmin && pbc[0] != 1.0) continue;
      if (inx > Nixmax && pbc[0] != 1.0) continue;
      for(diy=-1;diy<=1;diy++) {
	 iny=iy+diy;
	 if (iny < Niymin && pbc[1] == 1.0) iny=Niymax;
	 if (iny > Niymax && pbc[1] == 1.0) iny=Niymin;
	 if (iny < Niymin && pbc[1] != 1.0) continue;
	 if (iny > Niymax && pbc[1] != 1.0) continue;

	 for(diz=-1;diz<=1;diz++) {
	    inz=iz+diz;
	    if (inz < Nizmin && pbc[2] == 1.0) inz=Nizmax;
	    if (inz > Nizmax && pbc[2] == 1.0) inz=Nizmin;
	    if (inz < Nizmin && pbc[2] != 1.0) continue;
	    if (inz > Nizmax && pbc[2] != 1.0) continue;
                  
	    j=head[nind(inx,iny,inz)];
	    nj=0;
	    while (True) {
	       nj++;
	       if (nj>1) j=list[j];
	       if (j==-1) break;
	       if (j==i) continue; 

	       j3=j*3;
	       r1=fabs(x0[j3+0]-x0[i3+0]);
	       if (r1+r1>box[0] && pbc[0] == 1.0) r1=box[0]-r1;
	       if (r1 > Rneicut) continue;
	       r2=fabs(x0[j3+1]-x0[i3+1]);
	       if (r2+r2>box[1] && pbc[1] == 1.0) r2=box[1]-r2;
	       if (r2 > Rneicut) continue;
	       r3=fabs(x0[j3+2]-x0[i3+2]);
	       if (r3+r3>box[2] && pbc[2] == 1.0) r3=box[2]-r3;
	       if (r3 > Rneicut) continue;

	       rsq=r1*r1+r2*r2+r3*r3;

	       if (rsq > Rneicutsq) continue;

	       /* Accepted neighbour, add to list */
	       if (nn>=MAXNGBRS) {
		  printf("Too many ngbrs found iat %d nn %d",iat,nn);
		  printf("Increase MAXNGBRS from %d\n",MAXNGBRS);
		  exit(0);
	       }

	       ngbrs[nn].i = j;
	       ngbrs[nn].r = sqrt(rsq);
	       nn++;

	    }
	 }
      }
   }

   
   DEBUGSRRR("Found neighbours",iat,nn,nj);
   *nngot=nn;
   
   firsttime=False;
   
   return;
}

sqdistmax(real *x0,
	  int i,
	  int j,
	  real *box,
	  real *pbc,
	  real max,
	  real maxsq,
	  real *rsq) 
{
   real r1,r2,r3;
   int i3,j3;

   i3=i*3;
   j3=j*3;

   r1=fabs(x0[j3+0]-x0[i3+0]);
   if (r1+r1>box[0] && pbc[0]==1.0) r1=box[0]-r1;
   if (r1 > max) { *rsq=-1.0; return; }
   r2=fabs(x0[j3+1]-x0[i3+1]);
   if (r2+r2>box[1] && pbc[1]==1.0) r2=box[1]-r2;
   if (r2 > max) { *rsq=-1.0; return; }
   r3=fabs(x0[j3+2]-x0[i3+2]);
   if (r3+r3>box[2] && pbc[2]==1.0) r3=box[2]-r3;
   if (r3 > max) { *rsq=-1.0; return; }

   *rsq=r1*r1+r2*r2+r3*r3;
   if (*rsq > maxsq) { *rsq=-1.0; return; }

   return;
}


dist(real x1,
     real y1,
     real z1,
     real x2,
     real y2,
     real z2,
     real *box,
     real *pbc,
     real *r) 
{
   real xij,yij,zij;

   xij=(x1-x2); 
   if (pbc[0]==1.0) { if (xij+xij>box[0]) xij-=box[0]; else if (xij+xij<-box[0]) xij+=box[0]; }
   yij=(y1-y2); 
   if (pbc[1]==1.0) { if (yij+yij>box[1]) yij-=box[1]; else if (yij+yij<-box[1]) yij+=box[1]; }
   zij=(z1-z2); 
   if (pbc[2]==1.0) { if (zij+zij>box[2]) zij-=box[2]; else if (zij+zij<-box[2]) zij+=box[2]; }

   *r=sqrt(xij*xij+yij*yij+zij*zij);

}

real getangle(int i,
	      int j,
	      int k,
	      real *x0,
	      real *box,
	      real *pbc)
{
  real rij,xij,yij,zij;
  real rik,xik,yik,zik;
  real halfx,halfy,halfz;
  real costhetaijk,thetaijk,thetaijkdeg;

  halfx=box[0]/2.0;
  halfy=box[1]/2.0;
  halfz=box[2]/2.0;

  xij=x0[j*3+0]-x0[i*3+0]; 
  if (pbc[0]==1.0) { if (xij>halfx) xij-=box[0]; else if (xij<-halfx) xij+=box[0]; }
  yij=x0[j*3+1]-x0[i*3+1]; 
  if (pbc[1]==1.0) { if (yij>halfy) yij-=box[1]; else if (yij<-halfy) yij+=box[1]; }
  zij=x0[j*3+2]-x0[i*3+2]; 
  if (pbc[2]==1.0) { if (zij>halfz) zij-=box[2]; else if (zij<-halfz) zij+=box[2]; }
  rij=sqrt(xij*xij+yij*yij+zij*zij);

  xik=x0[k*3+0]-x0[i*3+0]; 
  if (pbc[0]==1.0) { if (xik>halfx) xik-=box[0]; else if (xik<-halfx) xik+=box[0]; }
  yik=x0[k*3+1]-x0[i*3+1]; 
  if (pbc[1]==1.0) { if (yik>halfy) yik-=box[1]; else if (yik<-halfy) yik+=box[1]; }
  zik=x0[k*3+2]-x0[i*3+2]; 
  if (pbc[2]==1.0) { if (zik>halfz) zik-=box[2]; else if (zik<-halfz) zik+=box[2]; }
  rik=sqrt(xik*xik+yik*yik+zik*zik);

  costhetaijk=(xij*xik+yij*yik+zij*zik)/(rij*rik);
  if (costhetaijk>1.0) costhetaijk=1.0;
  if (costhetaijk<-1.0) costhetaijk=-1.0;

  thetaijk=acos(costhetaijk);
  thetaijkdeg=thetaijk*180.0/pi;

  return(thetaijkdeg);

}
