
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int maxlevel=0;

/*
  To compile: cc -O  -o ringanalyze ringanalyze.o -lm


  ringanalyze -- code to analyze ring structure of an atom system.

  Based on  
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


struct lit {
  int i[6][3];
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

/* ring defs , ehk??my??emmin parametreiksi*/

#define mxlevel maxlevel
#define mxpaths 1000
#define mxlinks 32


/* ring globals */

int *lnks;
int *sring;
int *nodlnkd;
int *lvldist;
int *lvlref;
int *srtpth;
int *querng;
int *ringstat;
int *queue;


int pths;

int pat;
int dat;


/*
  Function protypes 
*/


void dijk(int nodsrc,int lvlreq, int *lnks,int *nodlnkd,int *lvldist,int nat);
void dijk2(int nodsrc,int lvlreq, int *lnks,int *nodlnkd,int *lvldist,int nat);
void srt_rec(int nodcrt,int lvlcrt,int lvlprim);
void pair_search(int *nodcrt,int *nodgoal,struct lit *limit, int *goal_found);


void ringanal(real *x0,
		    int *itype,
		    real *box,
		    real *pbc,
		    real cutR[MAXTYPE][MAXTYPE],
		    real cutS[MAXTYPE][MAXTYPE],
		    int nat,
		    real a,
		    real dr,
		    real time,
		    real printshort,
		    logical printangles,
		    logical writenbonds,
		    int printnot,
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

char name[MAXTYPE][8];
int *in,*itype;

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
   logical discardout=False,disc;
   int ntime,natexpectedp,natexpected;
   int ndiscard,ndiscardx,ndiscardy,ndiscardz;
   int ndiscardxfar,ndiscardyfar,ndiscardzfar;
   int printnot;
   double printshort,effectiveafact,effectivearmax;
   logical printangles,writenbonds;

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
   double cutR[MAXTYPE][MAXTYPE],cutS[MAXTYPE][MAXTYPE],cutSmax;

   double dr;
   double newbox[3];
   
   /* 
      Equilibrium unit cell size, many other variables are defined from
      this
   */

   double unitcella;
   int unitcelln;

   if ((argc>1 && strcmp(argv[1],"-h")==0.0) || argc < 2) {
      printf("Usage: ringa [options] file [stringcol [string [xdatacol [xsize [ysize [zsize]]]]]]\n");
      printf("                          Material always treated as alloy\n");
      printf("\n"); 
      printf("Options: -debug           debug\n");
      printf("         -maxdebug        debug even more\n");
      printf("         -cut R           set cutoff value\n");
      printf("         -cutij i j R     set cutoff balue for pair\n");
      printf("         -max L           maximum ring size\n");
      printf("         -pbc 0/1 0/1 0/1 (Un)Set periodic boundaries: 1 1 1 assumed\n");
    
      printf("\n");
     
      printf("\n");
      printf("\n");
      printf("Data is taken as x y z beginning from column xdatacol, from lines in which\n");
      printf("stringcol contains string. If stringcol==0, all lines are assumed to have xyz data\n");
      printf("\n");
      printf("Analyses atom coordinates for structure factor\n\n");
      printf("If string is not empty, takes the time step from any line not containing string but\n");
      printf("containing the string fs and analyses separately for each time. If string is empty\n");
      printf("any line with less than three arguments is interpreted to separate times.\n"); 
      printf("\n");
      printf("xsize, ysize and zsize give the box size\n");

      printf("\n");
      printf("The method is the 'prime ring' algorithm from Ref.\n");
      printf("Yuan and Cormack, Comput. Mater. Sci 24 (2002) 343.\n");
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

   maxlevel=9;

   /* Assume periodics */
   pbc[0]=pbc[1]=pbc[2]=1.0;

   /* Assume hard cutoff of 3.0 */
   for (i=0;i<MAXTYPE;i++)  for (j=0;j<MAXTYPE;j++) {
     cutR[i][j]=3.0;
     cutS[i][j]=3.0;
   }
   binarycut=-2;
   binaryt1=1;
   binaryt2=2;
   newbox[0]=1e30;
   printshort=-1.0;
   printangles=False;
   writenbonds=False;
   printnot=-1;
      
   while (strncmp(argv[iarg],"-",1)==0) {
     if (strcmp(argv[iarg],"-d")==0) details=True; 
     else if (strcmp(argv[iarg],"-debug")==0) debug=True; 
     else if (strcmp(argv[iarg],"-maxdebug")==0) maxdebug=True; 
     else if (strcmp(argv[iarg],"-discardout")==0) discardout=True; 
     else if (strcmp(argv[iarg],"-a")==0) sscanf(argv[++iarg],"%lg",&unitcella);
     else if (strcmp(argv[iarg],"-n")==0) sscanf(argv[++iarg],"%d",&unitcelln);
     else if (strcmp(argv[iarg],"-typecol")==0) sscanf(argv[++iarg],"%d",&typedatacol);

     else if (strcmp(argv[iarg],"-cut")==0) {
       sscanf(argv[++iarg],"%lg",&xx);
       //       sscanf(argv[++iarg],"%lg",&yy);
       for(i=0;i<MAXTYPE;i++) for(j=0;j<MAXTYPE;j++) {
	 cutR[i][j]=xx;
         cutS[i][j]=xx;
       }
       printf("Set all cutoffs between: R %g S %g. OVERRIDES ANY PREVIOUS -cutij!\n",xx,yy);
     }
     else if (strcmp(argv[iarg],"-max")==0) {
       sscanf(argv[++iarg],"%d",&maxlevel);       
       if((maxlevel%2)==0) maxlevel++;
       printf("Max ring size set to %d\n",maxlevel);             
       maxlevel++;
     }
     else if (strcmp(argv[iarg],"-cutij")==0) {
       sscanf(argv[++iarg],"%d",&i);
       sscanf(argv[++iarg],"%d",&j);
       sscanf(argv[++iarg],"%lg",&xx);
       //       sscanf(argv[++iarg],"%lg",&xx);
       cutR[i][j]=xx;
       cutS[i][j]=xx;
       cutR[j][i]=xx;
       cutS[j][i]=xx;
       printf("Set cutoff for pair %d %d and %d %d between: R %g S %g\n",i,j,j,i,xx,yy);
     }
     else if (strcmp(argv[iarg],"-dr")==0) {
       sscanf(argv[++iarg],"%lg",&dr);
       printf("Set drstat dr to %g\n",dr);
     }
     else if (strcmp(argv[iarg],"-printangles")==0) printangles=True;
     else if (strcmp(argv[iarg],"-writenbonds")==0) writenbonds=True;
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
     else if (strcmp(argv[iarg],"-binary")==0) {
       sscanf(argv[++iarg],"%d",&binarycut);
       sscanf(argv[++iarg],"%d",&binaryt1);
       sscanf(argv[++iarg],"%d",&binaryt2);
       printf("Getting binary neighbourhoood statistics\n");
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
   printf("Using cutoffs of:\n   overall: R %g S %g\n",cutR[5][5],cutS[5][5]);
   for (i=0;i<MAXTYPE;i++) for (j=0;j<MAXTYPE;j++) {
     if (cutS[i][j] > cutSmax) cutSmax=cutS[i][j];
     if (cutR[i][j] != cutR[5][5] || cutS[i][j] != cutS[5][5]) {
       printf("   pair %d - %d: R %g S %g\n",i,j,cutR[i][j],cutS[i][j]);
     }
   }

   if (printeffectivea && effectivearmax > cutSmax) {
     printf("ERROR: cutS too small for rmax calc. Do something !\n");
     exit(0);
   }

   MAXAT=10000;

   i=sizeof(double)*MAXAT*3+sizeof(int)*MAXAT;
   i+=sizeof(int)*MAXAT;
   printf("Allocating %d bytes for %d atoms ... ",i,MAXAT);
   x0=(double *) malloc(sizeof(double)*3*MAXAT);
   in=(int *)  malloc(sizeof(int)*MAXAT);
   itype=(int *) malloc(sizeof(int)*MAXAT);
   printf("... done\n");


   /* ring */
   ringstat=(int *)malloc(sizeof(int)*(mxlevel+2));
   queue=(int *)malloc(sizeof(int)*(MAXAT+1));
   lnks=(int *)malloc(sizeof(int)*(MAXAT+1));
   sring=(int *)malloc(sizeof(int)*(MAXAT+1)*(mxlinks+1));
   nodlnkd=(int *)malloc(sizeof(int)*(MAXAT+1)*(mxlinks+1));
   lvldist=(int *)malloc(sizeof(int)*(MAXAT+1));
   lvlref=(int *)malloc(sizeof(int)*(MAXAT+1)*5);
   srtpth=(int *)malloc(sizeof(int)*(mxpaths+1)*(mxlevel/2+1));
   querng=(int *)malloc(sizeof(int)*(mxpaths+1)*(mxpaths+1));

   
   


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
	 printf("Warning: long line : %s\n",buf);
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
	       /* Check whether time line has box size as well */
	       if (strcmp(arg[6],"boxsize")==0) {
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
		  sscanf(arg[7],"%lg",&(newbox[0]));
		  sscanf(arg[8],"%lg",&(newbox[1]));
		  sscanf(arg[9],"%lg",&(newbox[2]));
		  printf("Picked next cell size %lg %lg %lg from time line\n",
			 newbox[0],newbox[1],newbox[2]);
	       }
	       DEBUGSRRR("New time found",ntime,casctime,nline);
	    }
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

	    if (itype[nat]<itypemin) itypemin=itype[nat];
	    if (itype[nat]>itypemax) itypemax=itype[nat];
	    if (itype[nat]>=MAXTYPE) {
	      printf("MAXTYPE %d too small: increase to %d",MAXTYPE,itype[nat]+1);
	      exit(0);
	    }
	    if (xdatacol>0) strncpy(name[itype[nat]],arg[xdatacol-1],8);


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

               /* ring */
               queue=(int *)realloc(queue,sizeof(int)*(MAXAT+1));
               lnks=(int *)realloc(lnks,sizeof(int)*(MAXAT+1));               
               sring=(int *)realloc(sring,sizeof(int)*(MAXAT+1)*(mxlinks+1));
               nodlnkd=(int *)realloc(nodlnkd,sizeof(int)*(MAXAT+1)*(mxlinks+1));
               lvldist=(int *)realloc(lvldist,sizeof(int)*(MAXAT+1));
               lvlref=(int *)realloc(lvlref,sizeof(int)*(MAXAT+1)*5);               
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
	 printf("--------- Going into ring analysis at time %10g   ------------\n",prevtime);
	 printf("--------------------------------------------------------------------\n");
	 printf("\n");
	 
	
	 
	 
	          ringanal(x0,itype,box,pbc,cutR,cutS,nat,unitcella,dr,prevtime,
			printshort,printangles,writenbonds,printnot,printeffectivea,
			effectivearmax,effectiveafact);
	 
	 printf("\n");
	 printf("--------------------------------------------------------------------\n");
	 printf("---------- Done with ring analysis at time %10g   ------------\n",prevtime);
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
   
   ringanal(x0,itype,box,pbc,cutR,cutS,nat,unitcella,dr,prevtime,
			printshort,printangles,writenbonds,printnot,printeffectivea,
			effectivearmax,effectiveafact);   

  
   printf("\n");
   printf("--------------------------------------------------------------------\n");
   printf("---------- Done with defect analysis at time %10g ------------\n",casctime);
   printf("--------------------------------------------------------------------\n");
   printf("\n");
	 


} /* End of program */


/*

  Actual analysis of data. 
  
*/

#define MAXNGBRS 120



#define round(x) ((int) (x>0 ? x+0.5 : x-0.5))




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
*         x(i)       Atom coordinates scaled between -0.5 and 0.5
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
      printf("\nFirst time in getngbr: allocating %d bytes for linkcell arrays\n",
	     sizeof(int)*xnc*ync*znc+sizeof(int)*(N+1));

      printf("Cell ranges %d - %d, %d - %d, %d - %d\n",
	     Nixmin,Nixmax,Niymin,Niymax,Nizmin,Nizmax);
      printf("Cell size %g %g %g\n\n",cellsize[0],cellsize[1],cellsize[2]);
      fflush(stdout);

      head=(int *) malloc(sizeof(int)*xnc*ync*znc);
      list=(int *) malloc(sizeof(int)*(N+1));
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
               printf("analyzeint.c: This is impossible, atom %d x outside cell %d %g xsize %g\n",i,ix,x0[i3],box[0]);
	    }
            else {
               if (ix < Nixmin) ix=Nixmin;
               if (ix > Nixmax) ix=Nixmax;
            }
         }

         iy=round(x0[i3+1]/cellsize[1]);
         if (iy < Niymin || iy > Niymax) {
            if (pbc[1] == 1.0) {
               printf("analyzeint.c: This is impossible, atom %d y outside cell %d %g ysize %g\n",i,iy,x0[i3+1],box[1]);
	    }
            else {
               if (iy < Niymin) iy=Niymin;
               if (iy > Niymax) iy=Niymax;
            }
         }

         iz=round(x0[i3+2]/cellsize[2]);
         if (iz < Nizmin || iz > Nizmax) {
            if (pbc[2] == 1.0) {
               printf("analyzeint.c: This is impossible, atom %d z outside cell %d %g zsize %g\n",i,iz,x0[i3+2],box[2]); 
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
       printf("analyzeint.c: This is impossible, atom %d x outside cell %d %g xsize %g\n",i,ix,x0[i3],box[0]);
     }
     else {
       if (ix < Nixmin) ix=Nixmin;
       if (ix > Nixmax) ix=Nixmax;
     }
   }

   iy=round(x0[i3+1]/cellsize[1]);
   if (iy < Niymin || iy > Niymax) {
     if (pbc[1] == 1.0) {
       printf("analyzeint.c: This is impossible, atom %d y outside cell %d %g ysize %g\n",i,iy,x0[i3+1],box[1]);
     }
     else {
       if (iy < Niymin) iy=Niymin;
       if (iy > Niymax) iy=Niymax;
     }
   }
   iz=round(x0[i3+2]/cellsize[2]);
   if (iz < Nizmin || iz > Nizmax) {
     if (pbc[2] == 1.0) {
       printf("analyzeint.c: This is impossible, atom %d z outside cell %d %g zsize %g\n",i,iz,x0[i3+2],box[2]); 
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

void dijk(int nodsrc,int lvlreq, int *lnks,int *nodlnkd,int *lvldist,int nat) {
  int t;
  int quebgn,quend,nodcrt;
  int lnkscrt,nodprb,lvlprb;
  
  for(t=0;t<nat;t++) {lvldist[t]=lvlreq+2;}
  lvldist[nodsrc]=0;
  queue[0]=nodsrc;
  dat=0;


  for(quebgn=0,quend=1;quebgn<=quend;quebgn++) {   
    nodcrt=queue[quebgn];
    lvlprb=lvldist[nodcrt]+1;

    for(lnkscrt=0;lnkscrt<lnks[nodcrt];lnkscrt++) {      
      nodprb=nodlnkd[nodcrt*mxlinks+lnkscrt];
      if(lvldist[nodprb]>lvlprb) {
         lvldist[nodprb]=lvlprb;
         if(lvlprb<lvlreq) {
           queue[quend]=nodprb; // MB: increment 'quend' AFTER adding node to
           quend++;             // queue, otherwise 2nd node will always be
           dat++;               // node 0 (cf. Yuan 2002, arrays start from 1)
         }
      }
    }
    
  }
}




int srtpthx[2048];
void srt_rec(int nodcrt,int lvlcrt,int lvlprim){
  //  int srtpthx[mxlevel/2];

  int lnkcrt,nodprb,lvlprb,lvl;

  int t;
 
  
  for(lnkcrt=0;lnkcrt<lnks[nodcrt];lnkcrt++) {    
     nodprb=nodlnkd[nodcrt*mxlinks+lnkcrt];    
     lvlprb=lvldist[nodprb];
     if(lvlprb==0) {            
       for(lvl=1;lvl<=lvlprim+1;lvl++) {         
	 srtpth[((pths)*(mxlevel/2))+lvl]=srtpthx[lvl];

       }
       pths++;
       } else if(lvlprb==(lvlcrt-1)) {
 	      srtpthx[lvlprb]=nodprb;	     
              srt_rec(nodprb,lvlprb,lvlprim);
              }
     
  }
}


void prime_ring(int iodd,int lvlprim,int rngs,int rat) {

  int i,j,k,u;

  struct lit limit;

  int goal_found,irng,pth1,pth2,nodchk,lvlmax,lvlchk;
  int nodmid,irgx,p1x,p2x,t,rlev;

  goal_found=0;

  for(irng=0;irng<rngs;irng++) {   
    pth1=querng[irng];
    pth2=querng[irng+mxpaths*mxpaths/2];     
    if(pth1>=0) {
      for(lvlmax=lvlprim;lvlmax<=(lvlprim+iodd);lvlmax++) {        
	for(lvlchk=1;lvlchk<=(lvlmax-1);lvlchk++) {           // onko herra bugi t??ll???
	  nodchk=srtpth[pth1*(mxlevel/2)+lvlchk];
          nodmid=srtpth[pth2*(mxlevel/2)+lvlmax-lvlchk];                   

	  for(i=1;i<=4;i++) limit.i[i][1]=lvlref[nodmid*4+i]+lvlprim-1;
	  for(i=1;i<=4;i++) limit.i[i][2]=lvlref[nodmid*4+i]-lvlprim+1;
          limit.i[5][1]=lvldist[nodmid]+lvlprim-1;
          limit.i[5][2]=lvldist[nodmid]-lvlprim+1;
          pair_search(&nodchk,&nodmid,&limit,&goal_found);          
      
          if(goal_found==1) {
            goal_found=0;	    
            for(irgx=irng+1;irgx<rngs;irgx++) {
	      	      p1x=querng[irgx];
	                    p2x=querng[irgx+mxpaths*mxpaths/2];
	      if(p1x>=0) {	       
	      	 if((srtpth[p1x*(mxlevel/2)+lvlchk]==nodchk)&&
	           srtpth[p2x*(mxlevel/2)+lvlmax-lvlchk]==nodmid) querng[irgx]=-1;
	         if((srtpth[p2x*(mxlevel/2)+lvlchk]==nodchk)&&
	           srtpth[p1x*(mxlevel/2)+lvlmax-lvlchk]==nodmid) querng[irgx]=-1;	
	    
	      }	      
	      
            }
          goto nring;
          } // goal found
	}
      } // lvlmax loop
      rlev=2*lvlprim+iodd;
      ringstat[rlev]=ringstat[rlev]+1;     
      j=srtpth[pth1*(mxlevel/2)+1];
      k=srtpth[pth2*(mxlevel/2)+1];     
      for(u=0;u<lnks[rat];u++) {       
	if(nodlnkd[rat*mxlinks+u]==j)
	  if(sring[rat*mxlinks+u]>rlev) sring[rat*mxlinks+u]=rlev;
        if(nodlnkd[rat*mxlinks+u]==k)
	  if(sring[rat*mxlinks+u]>rlev) sring[rat*mxlinks+u]=rlev;
      }
    } // pth1 > 0
  nring:;
  }  
}

void pair_search(int *nodcrt,int *nodgoal,struct lit *limit, int *goal_found) {
  int lnkscrt,nodprb,iref;
  struct lit lmtx;


  if(*nodcrt==*nodgoal) {
    *goal_found=1;       
    
  } else
    {
      for(lnkscrt=0;lnkscrt<lnks[*nodcrt];lnkscrt++) {
	nodprb=nodlnkd[*nodcrt*mxlinks+lnkscrt]; 
	for(iref=1;iref<=4;iref++) {
	    if(lvlref[nodprb*4+iref]>=limit->i[iref][1]||
	    lvlref[nodprb*4+iref]<=limit->i[iref][2]) goto nsear;
       	} // iref loop
        if(lvldist[nodprb]>=limit->i[5][1]||
           lvldist[nodprb]<=limit->i[5][2]) goto nsear;
        for(iref=1;iref<=5;iref++) {
	  lmtx.i[iref][1]=limit->i[iref][1]-1;
          lmtx.i[iref][2]=limit->i[iref][2]+1;
	}
	pair_search(&nodprb,nodgoal,&lmtx,goal_found);	
        if(*goal_found==1) return;
      nsear:;
      } // lnkscrt loop
    } // if nodcrt
}





void ringanal(
		    real *x0,
		    int *itype,
		    real *box,
		    real *pbc,
		    real cutR[MAXTYPE][MAXTYPE],
		    real cutS[MAXTYPE][MAXTYPE],
		    int nat,
		    real a,
		    real dr,
		    real time,
		    real printshort,
		    logical printangles,
		    logical writenbonds,
		    int printnot,
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
   int typenngbrstat[MAXTYPE][MAXNGBRS];
   int typenngbrstattypestat[MAXTYPE][MAXNGBRS][MAXTYPE];
   int numberoftypetypepairs[MAXTYPE][MAXNGBRS][MAXNGBRS];

   int oneatomngbrtypestat[MAXTYPE];

   char file[160],bdfile[160];

   int jj,kk,nngot,nrave,it,jt,u,v;
   real r,boxmin,rave,r110ave,r1,r2,r3,theta;
   
   /* Neighbour parameters */
   static real rnmax;
   int npairs;         /* Total number of neighbours */
   real nbinaryt1,nbinaryt2;

   int drstat[1000000],ndr;

   /* Cutoff function params */
   real nngbr,nngbrmax,fc;

   int i,i3,j,k;
   int t1,t2,t3,t4,t5;

   static logical firsttime=True;

   FILE *fp,*fpbd,*fpa;






   // ***************************************************************
   // RING
   // ***************************************************************

   int nodcrt,lvlcrt,lvlprim;





   rnmax=cutS[0][0];
   for (i=0;i<MAXTYPE;i++) for (j=0;j<MAXTYPE;j++) {
     if (cutS[i][j] > rnmax) rnmax=cutS[i][j];
   }
   printf("Starting ring analysis at time %g for %d atoms. rcut %g\n",time,nat,rnmax);


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



   // Create neighbour list

   printf("Creating neighbour list...\n");
   for (i=0;i<nat;i++) {
     //      if (i>0&&i%10000==0) printf("%1d",(i%100000)/10000); 
      getngbrs(i,x0,box,pbc,nat,rnmax,ngbrs,&nngot);


      // ****************************************************************************************

      lnks[i]=nngot; // ring

      //      printf("-%d -%d :",i,lnks[i]);     
      t1=0;
      t3=itype[i];
      for(jj=0;jj<nngot;jj++)  {
	t2=itype[ngbrs[jj].i];
        if(cutS[t3][t2]>=ngbrs[jj].r) {	
         nodlnkd[i*mxlinks+t1]=ngbrs[jj].i;
         t1++;
 	 }                  
      }
      lnks[i]=t1;
      
      

      // ****************************************************************************************

      npairs+=nngot;
   } // End neighbourlist


   for(i=3;i<mxlevel;i++)
     ringstat[i]=0;

   for(i=0;i<(nat*mxlinks);i++) sring[i]=maxlevel+1;

   printf("Creating 4-point guides\n");

   dijk(nat/4*0,nat,lnks,nodlnkd,lvldist,nat);
   for(i=0;i<nat;i++) lvlref[i*4+1]=lvldist[i];

   dijk(nat/4*1,nat,lnks,nodlnkd,lvldist,nat);
   for(i=0;i<nat;i++) lvlref[i*4+2]=lvldist[i];

   dijk(nat/4*2,nat,lnks,nodlnkd,lvldist,nat);
   for(i=0;i<nat;i++) lvlref[i*4+3]=lvldist[i];

   dijk(nat/4*3,nat,lnks,nodlnkd,lvldist,nat);
   for(i=0;i<nat;i++) lvlref[i*4+4]=lvldist[i]; 
 
   // End 4-point guides

   pat=0;   


   // ring statistics of atom 1
   for(i=0;i<nat;i++){
     if((i%1000)==0)  printf("\n%d",i/1000);
     if((i%100)==0) {printf(".");fflush(stdout);} 

   
   // create shortest distance map
   dijk(i,(mxlevel/2)+1,lnks,nodlnkd,lvldist,nat);


   for(t5=0;t5<=dat;t5++) {
     jj=queue[t5];     
     if(lvldist[jj]<(mxlevel/2)) {          
     v=0;
     for(k=0;k<lnks[jj];k++)
       if((lvldist[jj]-1)==lvldist[nodlnkd[jj*mxlinks+k]]) v++;           
     if(v>1) {
     
        nodcrt=jj;
	lvlcrt=lvldist[jj];
	lvlprim=lvlcrt;
	pths=0;              
        srtpthx[lvlprim]=nodcrt;
	srt_rec(nodcrt,lvlcrt,lvlprim);        
        if(pat<pths) pat=pths;

        t4=0;
        for(t1=0;t1<pths;t1++) 
	  for(t2=t1+1;t2<pths;t2++) {
	    t3=0;
	    for(u=1;u<lvlprim;u++)
	       if( srtpth[((t1)*mxlevel/2)+u]==srtpth[((t2)*mxlevel/2)+u]) {t3=1;u=lvlprim+1;}
	  if(t3==0) {	    
	   
	    querng[t4]=t1;
	    querng[t4+mxpaths*mxpaths/2]=t2;
            t4++;
            if(t4>(mxpaths*mxpaths/2)) {printf("mxpaths overflow\n");exit(-1);}
  	  }       
	}
        if(t4>0) {
               prime_ring(0,lvlprim,t4,i);       
	}	
     }
     
     for(k=0;k<lnks[jj];k++)
       if(lvldist[jj]==lvldist[u=nodlnkd[jj*mxlinks+k]]) 
        if(u<jj) {	
	 nodcrt=jj;
	 lvlcrt=lvldist[jj];
	 lvlprim=lvlcrt;
	 pths=0;      
         srtpthx[lvlprim]=nodcrt;
	 srt_rec(nodcrt,lvlcrt,lvlprim);      
         j=pths;
         nodcrt=u;
         lvlcrt=lvldist[jj];
         lvlprim=lvlcrt;
         srtpthx[lvlprim]=nodcrt;
         srt_rec(nodcrt,lvlcrt,lvlprim);
         if(pat<pths) pat=pths;
         t4=0;
         for(t1=0;t1<j;t1++)
	   for(t2=j;t2<pths;t2++) {
	     t3=0;
	      for(u=1;u<lvlprim;u++) 
	       if( srtpth[((t1)*mxlevel/2)+u]==srtpth[((t2)*mxlevel/2)+u]) {t3=1;u=lvlprim+1;}	    
           if(t3==0) {	    
	   
	    querng[t4]=t1;
	    querng[t4+mxpaths*mxpaths/2]=t2;
            t4++;
            if(t4>(mxpaths*mxpaths/2)) {printf("mxpaths overflow\n");exit(-1);}
           }
           } // t2,t1
	 if(t4>0) {
               prime_ring(1,lvlprim,t4,i);
	 }

	   }
         
   } 

   }   
   }
   printf("\n");   analtime++;
   printf("Ring statistics : \n");
   for(i=3;i<mxlevel;i++)
     printf("Size %d rings %d\n",i,ringstat[i]/i);

   for(i=0;i<mxlevel;i++) ringstat[i]=0;
   for(i=0;i<nat;i++)
     for(t1=0;t1<lnks[i];t1++)
       ringstat[sring[i*mxlinks+t1]]++;

   printf("\nSmallest rings for bonds : \n");
   for(i=3;i<mxlevel;i++)
     printf("SSize %d bonds %d\n",i,ringstat[i]/2);
   



   fflush(stdout); fflush(stderr); 
   firsttime=False;

}
