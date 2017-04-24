
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define VERSION "structfactanalyze 25.6.2008"

/*
  structfactanalyze - code to get the structure factor of all atoms
  Code originally taken from fccanalyze and modified from it. 

  Structure factor code taken from PARCAS analyzeint.c

  NOTE: All arrays here start from 0, not 1 as in parcas !!

  Version history:

  24.6.2008      Added possible angles for non-FCC cells
                 Improved on 'emergency restart' handling
		 Added diasecond option

 25.6.2008       Added skinf parameters to reduce emergency restarts
                 Added element-specific output of structure factors.

*/

typedef int logical;
typedef double real;



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
#define MAXDEF 10000
#define MAXVACMOV 20000

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

void structfactanalyze(real *x0,
		       real *box,
		       real *pbc,
		       int nat,
		       real a,
		       real rnn,
		       real skinf,
		       real time,
		       logical printallsf,
		       real pstmax
		       );


void getngbrs(int iat,
	      real *x0,
	      real *box,
	      real *pbc,
	      int N,
	      real Rneicut,
	      int *ingbr,
	      real *rngbr,
	      int nngbr,
	      int *nngot,
	      real skinf);


struct point {
   double x;
   double y;
   double z;
};


logical details=False;  /* Shall we give details of the analysis as well ? */

char name[MAXTYPE][8];
int *in,*itype;

int lattice;
char string[80];
logical surf0=False;   /* Shall we set surface atom Pst == 0 ?? */
logical surfp=False;   /* Shall we set surface missing neighbours to be 'perfect' */
logical surfs=False;   /* Be silent about surface atoms (no 'emergency restart' messages) */
logical surfi=False;   /* Ignore surface effects completely, use too few atoms */
logical alloy;

int PSTSIZE;

#define FCC 1
#define BCC 2
#define DIA 3
#define GRA 4
#define DIASECOND 5
#define WURTZITE 6

main(argc,argv)
int argc;
char **argv;
{

   int iarg;
   int i,ii,nat;

   FILE *fp;
   char file[50];
   int nline;

   int stringcol,xdatacol;

   logical usetime,newtime,firstdatafound;
   logical iscoordline;
   logical printallsf=False,discardout=False,disc;
   int ntime,natexpectedp,natexpected;
   int ndiscard,ndiscardx,ndiscardy,ndiscardz;
   int ndiscardxfar,ndiscardyfar,ndiscardzfar;

   double xx,yy,zz;

   char buf[160],buforig[160],c1;
   char *argbuf,arg[40][40];
   int narg;

   /* structfactanalyze interface stuff */
   double *x0;
   real box[3],pbc[3];

   double casctime,prevtime;
   double forceminx=-1e30,forceminy=-1e30,forceminz=-1e30;
   real pstmax;

   double newbox[3];
   /* 
      Equilibrium unit cell size, many other variables are defined from
      this
   */

   double unitcella,rnn,skinf;
   int unitcelln;

   printf(" ---- %s -----\n",VERSION);

   if ((argc>1 && strcmp(argv[1],"-h")==0.0) || argc < 2) {
      printf("Usage: structfactanalyze [options] file [stringcol [string [xdatacol [xsize [ysize [zsize]]]]]]\n");
      printf("\n"); 
      printf("Options: -debug           debug\n");
      printf("         -maxdebug        debug even more\n");
      printf("         -d               give detailed information\n");
      printf("         -psf             Print out structure factor of all atoms\n");
      printf("         -a a             Unit cell a\n");
      printf("         -skinf f         Set skin factor, if you get emergency restarts increase this > 1.5 or so. Default 1.0\n");
      printf("         -n ncells        Number of unit cells, gives a\n");
      printf("         -pbc 0/1 0/1 0/1 (Un)Set periodic boundaries: 1 1 1 assumed\n");
      printf("         -discardout      Discard atoms outside cell completely\n");
      printf("         -surf0           Set surface atom Pst==0\n");
      printf("         -surfp           Set surface missing neighbours perf.\n");
      printf("         -surfs           Don't complain about surface atoms\n");
      printf("         -surfi           Ignore: allow too few atoms\n");

      printf("         -Pstsize n       Set nmumber of Pst points to print\n");
      printf("         -printngbr pst   Print neighbours of all atoms with Pst<pst\n");
      printf("         -alloy           Treat material as alloy\n");
      printf("\n");
      printf("         -fcc             Choose parameters apprpriate for FCC [Default]\n");
      printf("         -bcc             Choose parameters apprpriate for BCC \n");
      printf("         -dia             Choose parameters apprpriate for DIA \n");
      printf("         -diasecond       Choose parameters apprpriate for DIA including 2nd nearest neighbours \n");
      printf("         -wurtzite        Choose parameters apprpriate for Wurtzite including 2nd nearest neighbours \n");
      printf("         -gra             Choose parameters apprpriate for graphite. p_u set to 1 \n");
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
      printf("Default: 1 Au 2 114.24 114.24 114.24\n");
      printf("\n");
      return(0);
   }

   lattice=FCC;
   alloy=False;

   iarg=1;

   /* Set some basic variables which may be modified by command line arguments */


   unitcella=4.08;
   skinf=1.0;
   rnn=-1.0;        
   unitcelln=0;
   pstmax=-1;

   stringcol=1;
   sprintf(string,"XX"); 
   xdatacol=2;
   
   box[0]=114.24;
   box[1]=114.24;
   box[2]=114.24;
   PSTSIZE=40;

   /* Assume periodics */
   pbc[0]=pbc[1]=pbc[2]=1.0;

   newbox[0]=1e30;
      
   while (strncmp(argv[iarg],"-",1)==0) {
      if (strcmp(argv[iarg],"-d")==0) details=True; 
      else if (strcmp(argv[iarg],"-debug")==0) debug=True; 
      else if (strcmp(argv[iarg],"-maxdebug")==0) maxdebug=True; 

      else if (strcmp(argv[iarg],"-psf")==0) {
	 printallsf=True; 
	 printf("Printing sf of all atoms to file structurefacts.xyz\n");
      }
      else if (strcmp(argv[iarg],"-discardout")==0) discardout=True; 
      else if (strcmp(argv[iarg],"-surf0")==0) surf0=True; 
      else if (strcmp(argv[iarg],"-surfp")==0) surfp=True; 
      else if (strcmp(argv[iarg],"-surfs")==0) surfs=True; 
      else if (strcmp(argv[iarg],"-surfi")==0) surfi=True; 
      else if (strcmp(argv[iarg],"-alloy")==0) alloy=True; 
      else if (strcmp(argv[iarg],"-a")==0) sscanf(argv[++iarg],"%lg",&unitcella);
      else if (strcmp(argv[iarg],"-skinf")==0) sscanf(argv[++iarg],"%lg",&skinf);
      else if (strcmp(argv[iarg],"-rnn")==0) sscanf(argv[++iarg],"%lg",&rnn);
      else if (strcmp(argv[iarg],"-n")==0) sscanf(argv[++iarg],"%d",&unitcelln);

      else if (strcmp(argv[iarg],"-printngbr")==0) 
	 sscanf(argv[++iarg],"%lg",&pstmax);

      else if (strcmp(argv[iarg],"-Pstsize")==0) 
	 sscanf(argv[++iarg],"%d",&PSTSIZE);

      else if (strcmp(argv[iarg],"-pbc")==0) {
	 sscanf(argv[++iarg],"%lg",&(pbc[0]));
	 sscanf(argv[++iarg],"%lg",&(pbc[1]));
	 sscanf(argv[++iarg],"%lg",&(pbc[2]));
	 printf("Set periodics: x %g y %g z %g\n",pbc[0],pbc[1],pbc[2]);
      }
      else if (strcmp(argv[iarg],"-fcc")==0) lattice=FCC;
      else if (strcmp(argv[iarg],"-bcc")==0) lattice=BCC;
      else if (strcmp(argv[iarg],"-dia")==0) lattice=DIA;
      else if (strcmp(argv[iarg],"-diasecond")==0) lattice=DIASECOND;
      else if (strcmp(argv[iarg],"-wurtzite")==0) lattice=WURTZITE;
      else if (strcmp(argv[iarg],"-gra")==0) lattice=GRA;
      else {
	 fprintf(stderr,"Unknown option %s\n",argv[iarg]);
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
   if (alloy) {
      printf("Treating material as an alloy\n");
   }


   ntime=0;
   casctime=0.0;

   MAXAT=10000;

   i=sizeof(double)*MAXAT*3+sizeof(int)*MAXAT;
   if (alloy) i+=sizeof(int)*MAXAT;
   printf("Allocating %d bytes for %d atoms ... ",i,MAXAT);
   x0=(double *) malloc(sizeof(double)*3*MAXAT);
   in=(int *)  malloc(sizeof(int)*MAXAT);
   if (alloy) itype=(int *) malloc(sizeof(int)*MAXAT);
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
   firstdatafound=False;
   newbox[0]=box[0]; newbox[1]=box[1]; newbox[2]=box[2];
   ndiscard=ndiscardx=ndiscardy=ndiscardz=0;
   ndiscardxfar=ndiscardyfar=ndiscardzfar=0;
   while(fgets(buf,159,fp) != NULL) {
      nline++;

      /* Fix atom type from first xyz atom line if not specified above */
      if (nline==3 && strcmp(string,"XX")==0) {
	 strncpy(string,buf,2);
	 printf("Picked default atom name string %s from line %d\n",string,nline);
	 
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
	       if (!alloy && !(strcmp(arg[stringcol],string)==0))  {
		  printf("Impurity ignored: %s %s %s",arg[stringcol],string,buforig);
	       }
	       else {
		  iscoordline=True;
	       }
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
	    in[nat]=nat;
	    if (narg>=xdatacol+4) sscanf(arg[xdatacol+4],"%d",&(in[nat]));
	    if (alloy) {
	      sscanf(arg[xdatacol+3],"%d",&(itype[nat]));
	      if (itype[nat]<0) itype[nat]=-itype[nat];
	      if (xdatacol>0) strncpy(name[itype[nat]],arg[xdatacol-1],8);
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
	       if (alloy) i+=sizeof(int)*MAXAT;
	       printf("Reallocating %d bytes for %d atoms ... ",i,MAXAT);
	       x0=(double *) realloc(x0,3*sizeof(double)*MAXAT);
	       in=(int *)  realloc(in,sizeof(int)*MAXAT);
	       if (alloy) itype=(int *) realloc(itype,sizeof(int)*MAXAT);
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
	 printf("--------- Going into defect analysis at time %10g ------------\n",prevtime);
	 printf("--------------------------------------------------------------------\n");
	 printf("\n");
	 
	 structfactanalyze(x0,box,pbc,nat,unitcella,rnn,skinf,prevtime,printallsf,pstmax);

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
   
   structfactanalyze(x0,box,pbc,nat,unitcella,rnn,skinf,casctime,printallsf,pstmax);
   
   printf("\n");
   printf("--------------------------------------------------------------------\n");
   printf("---------- Done with defect analysis at time %10g ------------\n",casctime);
   printf("--------------------------------------------------------------------\n");
   printf("\n");
	 


} /* End of program */


/*

  Actual analysis of data. 

  Taken from old PARCAS analyzeint, Ekin and movement stuff thrown away.
  
*/


struct ngbr {
   int i;
   real r;
};



#define round(x) ((int) (x>0 ? x+0.5 : x-0.5))


void structfactanalyze(
		       real *x0,
		       real *box,
		       real *pbc,
		       int nat,
		       real a,
		       real rnn,
		       real skinf,
		       real time,
		       logical printallsf,
		       real pstmax
		       )
/*
   Remember that due to the nature of arrays in C and Fortran,
   any index in C is one less than that in Fortran !
*/

{

   static int analtime=0;   /* Counter of how many time we've been here */

   int nint,ii,jj,kk,nngot;
   logical liqflag,recflag;

   int nin,nout,nsame;
   real r0,r1,rin,rout;

   real r;
   
   /* Neighbour parameters */
   static real rnmax;
   static int nngbr;

   static int *atstr; 
   static real *pst,*rst;
   char buf[80];

   real angles[MAXNGBR][MAXNGBR];          /* Angles between neighbours */
   real dists[MAXNGBR];                    /* distances to neighbours */
   real sortedangles[MAXNGBR*MAXNGBR/2];   /* Angles of neighbours sorted by size */
   int ingbr[MAXNGBR];                     /* Atom indices of neighbours */
   real rngbr[MAXNGBR];    
   int npairs;                       /* Total number of neighbours */

   int nliquid,ndefect;
   int nliquidngbr,ndefngbrngbr;

   /* 
      Analysis rules (simplified a bit from original to be more general):
               ~~~~~

      Initially all atoms assumed to be crystalline
      Any atom with Pst between pstlmin and pstlmax => LIQUID
      Any atom with Pst between pstdmin and pstdmax => DEFECT
      Any atom which is LIQUID and has >= nngbr-nngbr/6 LIQUID neighbours => ISLIQUID
      Any atom which is DEFECT but not ISLIQUID => ISDEFECT

      Array atstr[nat] holds above info encoded bitwise:

         bit 1 set: LIQUID
	 bit 2 set: DEFECT
	 bit 4 set: LIQUID, and most neighbours liquid  => ISLIQUID
	 bit 5 set: DEFECT, but not ISLIQUID => ISDEFECT

   */

#define LIQUID 1
#define DEFECT 2
#define ISLIQUID 8
#define ISDEFECT 16

   static real pstlmin,pstlmax,pstdmin,pstdmax;

   int i,i3,j,t;
   int ipst,*Pststat,*Pststattype[MAXTYPE];
   real Pst,Pstmean;
   int irst,*Rststat,*Rststattype[MAXTYPE];
   real Rst,Rstmean;
   int nsurf;

   static logical firsttime=True;

   FILE *fp,*fptype[MAXTYPE];

   Pststat=(int *) malloc(sizeof(int)*(3*PSTSIZE+1));
   Rststat=(int *) malloc(sizeof(int)*(3*PSTSIZE+1));
   for (t=0;t<MAXTYPE;t++) {
     Pststattype[t]=(int *) malloc(sizeof(int)*(3*PSTSIZE+1));
     Rststattype[t]=(int *) malloc(sizeof(int)*(3*PSTSIZE+1));
   }

   if (alloy) {
      /* Find itypemin, itypemax */
      itypemax=-999999999;
      itypemin=+999999999;
      for(i=0;i<nat;i++) {
	 if (itype[i]>itypemax) itypemax=itype[i];
	 if (itype[i]<itypemin) itypemin=itype[i];
      }
      printf("Alloy mode using itypemin %d itypemax %d\n",itypemin,itypemax);
      printf("Atom types:"); for (t=itypemin;t<=itypemax;t++) printf(" %d %s",t,name[t]);
      printf("\n");
   }

   printf("Starting sfanalyze at time %g for %d atoms using a %g\n",time,nat,a);
   

   if (firsttime) {
      /* Selection of lattice structure properties. Also see getngbrangles */
      if (lattice==FCC) {
	 printf("Using FCC parameters, see Nordlund, PRB 56 (1997) 2421\n");
	 /* Neighbour dist criteria = average of first and second neighbour dist */
	 rnmax= a*(1+sqrt(2.0))/(2*sqrt(2.0));  
	 /* Fix number of nearest neighbours to ideal FCC value */
	 nngbr=12;
	 pstlmin=0.35; pstlmax=0.675;
	 pstdmin=0.2;  pstdmax=2.0;
      }
      else if (lattice==BCC) {
	 printf("Using BCC parameters, see Zhu, Phil. Mag. A 71 (1995) 735\n");
	 printf("Not supported yet, exiting !\n");
	 exit(0);
      }
      else if (lattice==DIA) {
	 pstlmin=0.125; pstlmax=2.0;
	 pstdmin=0.125; pstdmax=2.0;
	 nngbr=4; 
	 rnmax=a*(1.0/4.0*sqrt(3.0)+1.0/sqrt(2.0))/2.0; 
	 printf("Using DIA parameters, see Nordlund, PRB 56 (1997) 2421\n");
      } 
      else if (lattice==DIASECOND) {
	 pstlmin=0.3; pstlmax=2.0;
	 pstdmin=0.3; pstdmax=2.0;
	 nngbr=16; 
	 /* Neighbour dist criteria, same as in FCC ! */
	 rnmax= a*(1.0/sqrt(2.0)+sqrt(11.0)/4.0)/2.0;   
	 printf("Using DIASECOND parameters, see notes in ~knordlun/gecasc/sio2/README.amoanalysis \n");
      } 
      else if (lattice==WURTZITE) {
	 pstlmin=0.3; pstlmax=2.0;
	 pstdmin=0.3; pstdmax=2.0;
	 nngbr=17; 
	 /* Neighbour dist criteria, same as in FCC ! */
	 rnmax= a*(1.0/sqrt(2.0)+sqrt(11.0)/4.0)/2.0;  
	 printf("Using WURTZITE parameters \n");
      } 
      else if (lattice==GRA) {
	 pstlmin=0.14; pstlmax=2.0;
	 pstdmin=0.14; pstdmax=2.0;
	 nngbr=3; 
	 rnmax=a*(1.0/sqrt(3.0)+1.0)/2.0; 
	 printf("Using GRA parameters\n");
	 printf("Liquid=defect atoms, with Pst > 0.14\n");
      } 
      else {
	 printf("ERROR: Unknown lattice %d\n",lattice);
	 exit(0);
      }

      if (rnn > 0.0) {
	rnmax=rnn;
	printf("Set neighbour cutoff to %g (regardless of a and n values)\n",rnmax);
      }

      DEBUGSRRR("Lattice params:",rnmax,nngbr,a);
      DEBUGSRRRR("Pst params",pstlmin,pstlmax,pstdmin,pstdmax);

      printf("Allocating sf analysis arrays\n");
      atstr=(int *) malloc(sizeof(int)*(nat+10));
      pst=(real *) malloc(sizeof(real)*(nat+10));
      rst=(real *) malloc(sizeof(real)*(nat+10));
      
      remove("sfdefects.out");
      remove("sfliquidat.out");
      remove("Pststat.out");
      remove("Rststat.out");
   }
   for (i=0;i<nat;i++) atstr[i]=0;

   DEBUGSRRRRR("sfanalyze",nat,a,rnmax,nngbr,time);

   /* Initialize neighbour list calc. */
   i=-1;
   getngbrs(i,x0,box,pbc,nat,rnmax,ingbr,rngbr,nngbr,&nngot,skinf);


   /* First loop over all atoms: get Pst */
   npairs=0;
   nsurf=0;
   Pstmean=0.0;
   Rstmean=0.0;
   for (i=0;i<=3*PSTSIZE;i++) Pststat[i]=0;
   for (i=0;i<=3*PSTSIZE;i++) Rststat[i]=0;
   if (alloy) {
     for (t=0;t<MAXTYPE;t++) {
       for (i=0;i<=3*PSTSIZE;i++) Pststattype[t][i]=0;
       for (i=0;i<=3*PSTSIZE;i++) Rststattype[t][i]=0;
     }
   }

   for (i=0;i<nat;i++) {
      if (i>0&&i%1000==0) printf("%1d",(i%10000)/1000); 
      getngbrs(i,x0,box,pbc,nat,rnmax,ingbr,rngbr,nngbr,&nngot,skinf);
 //printf("numbers of neigbors: %d %d\n",nngbr,nngot);
      npairs+=nngot;
      DEBUGSRRR("Got ngbrs",i,nngbr,nngot);
      if (nngot==nngbr || surfp || surfi) {
	 if ((surfp||surfi) && nngot < nngbr) nsurf++;
	 getngbrangles(i,x0,box,pbc,ingbr,nngbr,nngot,angles,
		       sortedangles,dists,&Pst,&Rst,a); 
	 DEBUGSRRRR("Got for",i,nngbr,Pst,Rst); 
      }
      else {
	 if (surf0) { Pst=0.0; Rst=0.0; nsurf++; }
	 else {
	    printf("Error: nngot %d != nngbr %d\n",nngot,nngbr);
	    exit(0);
	 }
      }
      pst[i]=Pst;
      rst[i]=Rst;

      if (Pst < pstmax) {
	 printf("Atom %g %g %g %g %g %d has neighbours:\n",
		x0[i*3],x0[i*3+1],x0[i*3+2],Pst,Rst,in[i]);
	 for (jj=0;jj<nngot;jj++) {
	    j=ingbr[jj]; 
	    printf("   %g %g %g %d r %g\n",
		   x0[j*3],x0[j*3+1],x0[j*3+2],in[j],rngbr[jj]);
	 }		
      }

      Pstmean+=Pst;
      ipst=(int) (0.5+Pst*PSTSIZE); if (ipst>3*PSTSIZE) ipst=3*PSTSIZE;

      /* printf("debug i %d  ipst %d Pst %g\n",i,ipst,Pst); */

      Pststat[ipst]++;
      Rstmean+=Rst;
      irst=(int) (0.5+4*Rst*PSTSIZE); if (irst>3*PSTSIZE) irst=3*PSTSIZE;
      Rststat[irst]++;
      if (alloy) {
	Pststattype[itype[i]][ipst]++;
	Rststattype[itype[i]][irst]++;
      }

      if (Pst>pstlmin && Pst<pstlmax) atstr[i]|=LIQUID;
      if (Pst>pstdmin && Pst<pstdmax) atstr[i]|=DEFECT;
   }
   printf("\n");

   /* Do and write out Pst and Rst statistics */
   Pstmean/=nat;
   Rstmean/=nat;
   printf("Pstmean %g Rstmean %g\n",Pstmean,Rstmean); 

   if (firsttime) fp=fopen("Pststat.out","w");
   else fp=fopen("Pststat.out","a");
   if (alloy) {
     for(t=itypemin;t<=itypemax;t++) {
       sprintf(buf,"Pststat_%d.out",t);
       if (firsttime) fptype[t]=fopen(buf,"w");
       else fptype[t]=fopen(buf,"a");
     }
   }
   
   for (i=0;i<=3*PSTSIZE;i++) fprintf(fp,"%.4g %g %d\n",time,1.0*i/PSTSIZE,Pststat[i]);
   fflush(fp); fclose(fp); 
   if (alloy) {
     for(t=itypemin;t<=itypemax;t++) {
       for (i=0;i<=3*PSTSIZE;i++) fprintf(fptype[t],"%.4g %g %d\n",time,1.0*i/PSTSIZE,Pststattype[t][i]);
       fflush(fptype[t]); fclose(fptype[t]); 
     }
   }
   
   if (firsttime) fp=fopen("Rststat.out","w");
   else fp=fopen("Rststat.out","a");
   if (alloy) {
     for(t=itypemin;t<=itypemax;t++) {
       sprintf(buf,"Rststat_%d.out",t);
       if (firsttime) fptype[t]=fopen(buf,"w");
       else fptype[t]=fopen(buf,"a");
     }
   }


   for (i=0;i<=PSTSIZE;i++) fprintf(fp,"%.4g %g %d\n",time,1.0*i/4.0/PSTSIZE,Rststat[i]);
   fflush(fp); fclose(fp); 
   if (alloy) {
     for(t=itypemin;t<=itypemax;t++) {
       for (i=0;i<=PSTSIZE;i++) fprintf(fptype[t],"%.4g %g %d\n",time,1.0*i/4.0/PSTSIZE,Rststattype[t][i]);
       fflush(fptype[t]); fclose(fptype[t]); 
     }
   }
   
   /* Second loop over atoms: determine atom type finally */

   DEBUGSR("Second loop over atoms",nat);
   nliquid=ndefect=0;
   for (i=0;i<nat;i++) {
      if (i>0&&i%1000==0) printf("%1d",(i%10000)/1000);
      getngbrs(i,x0,box,pbc,nat,rnmax,ingbr,rngbr,nngbr,&nngot,skinf);
      
      /* Analyze type of neighbours */
      nliquidngbr=0;
      for (jj=0;jj<nngot;jj++) {
	 j=ingbr[jj]; 
	 if (i==j) continue;
	 if (atstr[j]&LIQUID) nliquidngbr++;
      }
      
      if ((atstr[i]&LIQUID) && nliquidngbr>=nngbr-nngbr/6) { 
	 atstr[i]|=ISLIQUID; nliquid++; 
      }
      if (atstr[i]&DEFECT && !(atstr[i]&ISLIQUID)) {
	 atstr[i]|=ISDEFECT; ndefect++;
      }
   }
   printf("\n");

   printf("Found %d liquid, %d defect and %d surface atoms\n",
	  nliquid,ndefect,nsurf);

   if (printallsf) {
      printf("Printing out sf of all atoms to structurefacts.xyz...");
      if (firsttime) fp=fopen("structurefacts.xyz","w");
      else fp=fopen("structurefacts.xyz","a");
      fprintf(fp," %d\n",nat);
      fprintf(fp," Structure factors of all atoms, time %g\n",time);
      for(i=0;i<nat;i++) {
	 i3=i*3;
	 if (alloy) strncpy(buf,name[itype[i]],8);
	 else strncpy(buf,string,8);
	 fprintf(fp,"%s %.12g %.12g %.12g %.8g %.8g %d %d \n",
		buf,x0[i3+0],x0[i3+1],x0[i3+2],pst[i],rst[i],atstr[i],in[i]);
      }
      fflush(fp);
      fclose(fp);
      printf(" ... Done.\n");
   }


   /* General output of stuff */
   if (details) printf("Writing out liquid and defect data\n");
   
   if (firsttime) fp=fopen("sfliquidat.out","w");
   else fp=fopen("sfliquidat.out","a");
   for (i=0;i< nat;i++) {
      if (atstr[i]&ISLIQUID) {
	 fprintf(fp,"%g %g %g %g %d %g\n",
		 x0[i*3],x0[i*3+1],x0[i*3+2],time,i,pst[i]);
      }
   }
   fflush(fp);
   fclose(fp);

   if (firsttime) fp=fopen("sfdefects.out","w");
   else fp=fopen("sfdefects.out","a");
   for (i=0;i< nat;i++) {
      if (atstr[i]&ISDEFECT) {
	 fprintf(fp,"%g %g %g %g %d %g\n",
		 x0[i*3],x0[i*3+1],x0[i*3+2],time,i,pst[i]);
      }
   }
   fflush(fp);
   fclose(fp);

   analtime++;

   fflush(stdout); fflush(stderr); 
   firsttime=False;

}


/*******************************************************************************/

/*******************************************************************************/

/*******************************************************************************/

getngbrangles(
	      int iat,
	      real *x0,
	      real *box,
	      real *pbc,
	      int *ingbr,
	      int nngbr,
	      int nngot,
	      real angles[MAXNGBR][MAXNGBR],
	      real *sortedangles,
	      real dists[MAXNGBR],
	      real *Pst,
	      real *Rst,
	      real a
	      )

/*
   Gets angles between neighbours and calculates the Zhu 
   structural order parameter Pst [Zhu, Phil. Mag. A 71 (1995) 735]

   Also calculate Rst, similar structure parameter for distances

   Note the differences for different lattice structures, selected by
   global variable lattice = FCC/BCC/DIA/GRA 

   External varibale surfp: if True, pretend non-exiting neighbours
   are perfect.
   External variable surfi: if True, don't even use missing angles

*/

{
   int i,j,k,jj,kk;
   real theta;

   real xij,yij,zij,rij;
   real xik,yik,zik,rik;

   real costhetaijk;
   real thetaijk,thetaijkdeg;
   real halfx,halfy,halfz;
   real d,dmin,angsq;

   int nang;

   static int nangles,npossible;
   static real perfectangles[200],perfectdist;
   static real possibleangles[20];
   static real uniformangdist[200];
   static logical firsttime=True;
   static real pu,ru;

   int intcmp();

   if (firsttime) {
      firsttime=False;
      nangles=(nngbr*(nngbr-1))/2;
      if (lattice == FCC) {
	 for (i=0;i<24     ;i++) perfectangles[i]=pi/3;
	 for (;i<24+12     ;i++) perfectangles[i]=pi/2;
	 for (;i<24+12+24  ;i++) perfectangles[i]=2*pi/3;
	 for (;i<24+12+24+6;i++) perfectangles[i]=pi;
	 possibleangles[0]=pi/3;
	 possibleangles[1]=pi/2;
	 possibleangles[2]=2*pi/3;
	 possibleangles[3]=pi;
	 npossible=4;
	 perfectdist=a/sqrt(2.0);
      }
      else if (lattice == BCC) {
	// printf("getngbrangles: BCC not supported yet, exiting !\n");
	// exit(0);
	perfectdist=a/2.0*sqrt(3.0);
	possibleangles[0]=  
      }
      else if (lattice == DIA) {
	 for(i=0;i<nangles;i++) perfectangles[i]=acos(-1.0/3.0);
	 perfectdist=a/4.0*sqrt(3.0);
	 possibleangles[0]=acos(-1.0/3.0);
	 npossible=1;
      }
      else if (lattice == DIASECOND) {
	 /* In degrees angles are:
	    12 35.2644     
	    24 60.0000     
	    36 90.0000     
	    6 109.4712     
	    24 120.0000     
	    12 144.7356     
	    6 180.0000     
	 */
	 for (i=0;i<12     ;i++) perfectangles[i]=35.2644/180*pi;
	 for (;i<12+24     ;i++) perfectangles[i]=pi/3;
	 for (;i<12+24+36     ;i++) perfectangles[i]=pi/2;
	 for (;i<12+24+36+6   ;i++) perfectangles[i]=acos(-1.0/3.0);
	 for (;i<12+24+36+6+24   ;i++) perfectangles[i]=2*pi/3;
	 for (;i<12+24+36+6+24+12   ;i++) perfectangles[i]=144.735/180*pi;
	 for (;i<12+24+36+6+24+12+6   ;i++) perfectangles[i]=pi;
	 perfectdist=a/4.0*sqrt(3.0);
	 possibleangles[0]=35.2644/180*pi;
	 possibleangles[1]=pi/3;
	 possibleangles[2]=pi/2;
	 possibleangles[3]=acos(-1.0/3.0);
	 possibleangles[4]=2*pi/3;
	 possibleangles[5]=144.735/180*pi;
	 possibleangles[6]=pi;
	 npossible=7;
      }

      else if (lattice == WURTZITE) {
	 /* In degrees angles are:
	    15 35.2644
	    24 60.0000
	     3 70.5288
	     3 74.2068
	    36 90.0000     
	     9 109.4712     
	    18 120.0000
	     6 122.9790
	    12 144.7356 
	     6 146.4427
	     4 180.0000     
	 */
	 for (i=0;i<15     ;i++) perfectangles[i]=35.2644/180.0*pi;
	 for (;i<15+24     ;i++) perfectangles[i]=60.0000/180.0*pi;
	 for (;i<15+24+3     ;i++) perfectangles[i]=70.5288/180.0*pi;
	 for (;i<15+24+3+3   ;i++) perfectangles[i]=74.2068/180.0*pi;
	 for (;i<15+24+3+3+36   ;i++) perfectangles[i]=90.0000/180.0*pi;
	 for (;i<15+24+3+3+36+9   ;i++) perfectangles[i]=109.4712/180.0*pi;
	 for (;i<15+24+3+3+36+9+18   ;i++) perfectangles[i]=120.0000/180.0*pi;
	 for (;i<15+24+3+3+36+9+18+6   ;i++) perfectangles[i]=122.9790/180.0*pi;
	 for (;i<15+24+3+3+36+9+18+6+12   ;i++) perfectangles[i]=144.7356/180.0*pi;
	 for (;i<15+24+3+3+36+9+18+6+12+6   ;i++) perfectangles[i]=146.4427/180.0*pi;
	 for (;i<15+24+3+3+36+9+18+6+12+6+4   ;i++) perfectangles[i]=180.0000/180.0*pi;
	 perfectdist=a/4.0*sqrt(3.0);
	 possibleangles[0]=35.2644/180.0*pi;
	 possibleangles[1]=60.0000/180.0*pi;
	 possibleangles[2]=70.5288/180.0*pi;
	 possibleangles[3]=74.2068/180.0*pi;
	 possibleangles[4]=90.0000/180.0*pi;
	 possibleangles[5]=109.4712/180.0*pi;
	 possibleangles[6]=120.0000/180.0*pi;
	 possibleangles[7]=122.9790/180.0*pi;
	 possibleangles[8]=144.7356/180.0*pi;
	 possibleangles[9]=146.4427/180.0*pi;
	 possibleangles[10]=180.0000/180.0*pi;
	 npossible=11;
      }

      else if (lattice == GRA) {
	 for(i=0;i<nangles;i++) perfectangles[i]=2*pi/3;
	 perfectdist=a;
	 possibleangles[0]=pi/3;
	 npossible=1;
      }

      for (i=0;i<nangles;i++) {
	 uniformangdist[i]=pi/(1.0*nangles)*(i+1);
      }
      pu=0.0;
      for(i=0;i<nangles;i++) { 
	 pu+=(uniformangdist[i]-perfectangles[i])*(uniformangdist[i]-perfectangles[i]);
      }
      pu=sqrt(pu);
      if (lattice == GRA) pu=1.0;
      printf("Uniform atom distr. order factor pu = %g\n",pu);
      
   }
   i=iat;

   halfx=box[0]/2.0;
   halfy=box[1]/2.0;
   halfz=box[2]/2.0;

   DEBUGSRRRR("Starting getngbrangles",i,nngbr,ingbr[0],nangles); 

   nang=0;
   for (jj=0;jj<nngot;jj++) {
      j=ingbr[jj];
      MAXDEBUGSRRR("  ",i,j,jj);
      if (i==j) continue;
     
      xij=x0[j*3+0]-x0[i*3+0]; 
      if (pbc[0]==1.0) { if (xij>halfx) xij-=box[0]; else if (xij<-halfx) xij+=box[0]; }
      yij=x0[j*3+1]-x0[i*3+1]; 
      if (pbc[1]==1.0) { if (yij>halfy) yij-=box[1]; else if (yij<-halfy) yij+=box[1]; }
      zij=x0[j*3+2]-x0[i*3+2]; 
      if (pbc[2]==1.0) { if (zij>halfz) zij-=box[2]; else if (zij<-halfz) zij+=box[2]; }

      rij=sqrt(xij*xij+yij*yij+zij*zij);

      dists[jj]=rij;

      for (kk=jj+1;kk<nngot;kk++) {
	 k=ingbr[kk];
	 MAXDEBUGSRRRR("     ",i,j,k,kk);
	 if (j==k || i==k) continue;
	 
	 xik=x0[k*3+0]-x0[i*3+0]; 
	 if (pbc[0]==1.0) { if (xik>halfx) xik-=box[0]; else if (xik<-halfx) xik+=box[0]; }
	 yik=x0[k*3+1]-x0[i*3+1]; 
	 if (pbc[1]==1.0) { if (yik>halfy) yik-=box[1]; else if (yik<-halfy) yik+=box[1]; }
	 zik=x0[k*3+2]-x0[i*3+2]; 
	 if (pbc[2]==1.0) { if (zik>halfz) zik-=box[2]; else if (zik<-halfz) zik+=box[2]; }
	    
	 rik=sqrt(xik*xik+yik*yik+zik*zik);
	 if (rij==0.0 || rik==0.0) { 
	    printf(" i %d j %d k %d rij %g rik %g\n",i,j,k,rij,rik);
	    printf(" theta: i %d j %d k %d %g\n",i,j,k,thetaijkdeg);
	    printf("  xi %g %g %g\n",x0[i*3+0],x0[i*3+1],x0[i*3+2]);
	    printf("  xj %g %g %g\n",x0[j*3+0],x0[j*3+1],x0[j*3+2]);
	    printf("  xk %g %g %g\n",x0[k*3+0],x0[k*3+1],x0[k*3+2]);
	    printf("HORROR ERROR\n"); fflush(NULL); 
	 }
	 costhetaijk=(xij*xik+yij*yik+zij*zik)/(rij*rik);
	 if (costhetaijk>1.0) costhetaijk=1.0;
	 if (costhetaijk<-1.0) costhetaijk=-1.0;

	 thetaijk=acos(costhetaijk);
	 thetaijkdeg=thetaijk*180.0/pi;
	 angles[jj][kk]=thetaijk;
	 angles[kk][jj]=thetaijk;
	 sortedangles[nang]=thetaijk; 
	 nang++;
	
	 /* bondangledistr[(int) ((thetaijkdeg+bondanglestep/2.0)/bondanglestep)]++;*/
	 if (maxdebug) {
	    printf(" i %d j %d k %d rij %g rik %g\n",i,j,k,rij,rik);
	    printf(" theta: i %d j %d k %d %g\n",i,j,k,thetaijkdeg);
	    printf("  xi %g %g %g\n",x0[i*3+0],x0[i*3+1],x0[i*3+2]);
	    printf("  xj %g %g %g\n",x0[j*3+0],x0[j*3+1],x0[j*3+2]);
	    printf("  xk %g %g %g\n",x0[k*3+0],x0[k*3+1],x0[k*3+2]);
	    fflush(stderr); fflush(stdout);
	 }
      } /* End of k loop */
      
   } /* End of j loop */

   DEBUGSRR("Got angles",nang,nangles);

   /* Sort the stuff using the ANSI C sorting routine */
   if (maxdebug) { 
      printf("Angles before sort:");
      for (j=0;j<nang;j++) printf(" %.3f",sortedangles[j]*180.0/pi); 
      printf("\n"); 
   } 
   qsort(sortedangles,nang,sizeof(real),&intcmp);
   if (maxdebug) { 
      printf("\nAngles after sort:");
      for (j=0;j<nang;j++) printf(" %.3f",sortedangles[j]*180.0/pi); 
      printf("\n"); 
   } 

   /* Get the Pst structure factor */
   *Pst=0.0;
   for(j=0;j<nang;j++) { 
      if (surfi || !surfp || (nngot==nngbr)) {
	angsq=(sortedangles[j]-perfectangles[j])*(sortedangles[j]-perfectangles[j]);
	*Pst+=angsq;
	 
	MAXDEBUGSRRRRR("   ",i,*Pst,angsq,sortedangles[j],perfectangles[j]);
      }
      else {
	 /* Find nearest possible perfect angle for surface atom angles */
	 dmin=1e30;
	 for(k=0;k<npossible;k++) {
	    d=(sortedangles[j]-possibleangles[k])*(sortedangles[j]-possibleangles[k]);
	    if (d<dmin) dmin=d;
	 }
	 DEBUGSRRRR("   ",i,dmin,nngot,nngbr);
	 *Pst+=dmin;
      }
   }
   *Pst=sqrt(*Pst)/pu;
   MAXDEBUGSRR("   ",i,*Pst);

   /* Get the Rst structure factor */
   *Rst=0.0;
   for (j=0;j<nngot;j++) *Rst+=(dists[j]-perfectdist)*(dists[j]-perfectdist);
      

   *Rst=sqrt(*Rst)/nngbr;

   DEBUGSRRR("\ngetngbrangles",*Pst,*Rst,nngbr);

}

/*******************************************************************************/
/*******************************************************************************/

/*
*
* Linked list subroutine written by Kai Nordlund. 
*
*  Special version to find exactly nngbr nearest neighbours.
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
*         nngbr      Number of neighbours to be found, nngbr
*
* Output: ingbr[]    Atom indices of neighbours
*
*/


void getngbrs(int iat,
	      real *x0,
	      real *box,
	      real *pbc,
	      int N,
	      real Rneicut,
	      int *ingbr,
	      real *rngbr,
	      int nngbr,
	      int *nngot,
	      real skinf)

{

/*
*  Size of one cell is set to exactly Rneicut.
*  Cells are distributed so that cell 0 is between 0 and Rneicut.
*  However, cell size may vary between calls of the subroutine.
*
*/

#define nind(a,b,c) (((a-Nixmin)*ync+(b-Niymin))*znc+(c-Nizmin))

#define MAXMAXNGBR 200

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

   struct ngbr ngbrs[MAXMAXNGBR];

   real rsq;
   int i,ii,i3,j,jj,j3,nn,nj;
   real r1,r2,r3;
   int ix,iy,iz;
   int dix,diy,diz;
   int inx,iny,inz;

   static real cellsize[3];
   real Rneicutsq;
   
   static logical firsttime=True;
   logical emergencyrestart;

   int ngbrcmp();

   emergencyrestart=False;

emergency:   /* Emergency restart if not enough neighbours found */

   MAXDEBUGSRRR("Starting getngbrs",Rneicut,iat,N); 

   Rneicutsq=Rneicut*Rneicut;

   Nixmax=(int)(box[0]/2.0/(skinf*Rneicut)-0.5);
   Nixmin=-Nixmax;
   Niymax=(int)(box[1]/2.0/(skinf*Rneicut)-0.5);
   Niymin=-Niymax;
   Nizmax=(int)(box[2]/2.0/(skinf*Rneicut)-0.5);
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

   if (Nixmax > Nixmax0 || Niymax > Niymax0 || Nizmax > Nizmax0) {
      
      printf("\ngetngbr: cell size has changed, changing array sizes\n");
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

   if (iat==-1 || emergencyrestart) {
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
               printf("analyzeint.c: This is impossible, atom x outside cell %g %g\n",ix,x0[i3]);
	    }
            else {
               if (ix < Nixmin) ix=Nixmin;
               if (ix > Nixmax) ix=Nixmax;
            }
         }

         iy=round(x0[i3+1]/cellsize[1]);
         if (iy < Niymin || iy > Niymax) {
            if (pbc[1] == 1.0) {
               printf("analyzeint.c: This is impossible, atom y outside cell %g %g\n",iy,x0[i3+1]);
	    }
            else {
               if (iy < Niymin) iy=Niymin;
               if (iy > Niymax) iy=Niymax;
            }
         }

         iz=round(x0[i3+2]/cellsize[2]);
         if (iz < Nizmin || iz > Nizmax) {
            if (pbc[2] == 1.0) {
               printf("analyzeint.c: This is impossible, atom z outside cell %g %g\n",iz,x0[i3+2]); 
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
      if (!emergencyrestart) return;

   } /* End of cell initialization */

   emergencyrestart=False; 

   /* Find neighbours of atom iat */

   for(i=0;i<nngbr;i++) {
      ngbrs[i].r=0;
      ngbrs[i].i=-1;
   }
   i=iat;
   i3=3*i;

   ix=round(x0[i3+0]/cellsize[0]);
   iy=round(x0[i3+1]/cellsize[1]);
   iz=round(x0[i3+2]/cellsize[2]);

   /* MAXDEBUGSRRRRR("Finding neighbours",iat,i3,ix,iy,iz); */

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
	       if (nn>MAXMAXNGBR) {
		  printf("Too many ngbrs found iat %d nn %d",iat,nn);
		  exit(0);
	       }

	       ngbrs[nn].i = j;
	       ngbrs[nn].r = rsq;
	       nn++;

	    }
	 }
      }
   }

   
   DEBUGSRRR("Found neighbours",iat,nn,nj);
   *nngot=(nn>nngbr?nngbr:nn);
   
   if (nn<nngbr) {
      if (!surf0 && !surfp && !surfi) {
	 /*
	    printf("getngbr error: failed to find %d neighbours for atom %d,\n",nngbr,i);
	    printf("Trying emergency restart with rnew %g = 1.1*%g\n",1.1*Rneicut,Rneicut);
	 */
	 Rneicut=Rneicut*1.05;
 	 if (!surfs) printf("Emergency restart new rcut: %g %d %d\n",Rneicut,nn,i);
	 Rneicutsq=Rneicut*Rneicut;
	 emergencyrestart=True;
	 goto emergency;
      }
   }

   /* Sort neighbours according to r taking i with you */
   qsort(ngbrs,nn,sizeof(struct ngbr),&ngbrcmp);
   /* for(j=0;j<nn;j++) printf(" %5d %9.3f",ngbrs[j].i,ngbrs[j].r); printf("\n"); */

   /* Copy the nearest nngot neighbours to ingbr and rngbr */
   if (maxdebug) printf("Neighbours of %d:",i);
   for(j=0;j<*nngot;j++) {
      ingbr[j]=ngbrs[j].i;
      rngbr[j]=sqrt(ngbrs[j].r);
      if (maxdebug) printf("%d %g, ",ingbr[j],rngbr[j]);
   }

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

/*******************************************************************************/
int intcmp(real *x, real *y)
{
   if (*x < *y) return(-1);
   else if (*x > *y) return(1);
   else return(0);

}

/*******************************************************************************/
int ngbrcmp(struct ngbr *x, struct ngbr *y)
{
   if ((*x).r < (*y).r) return(-1);
   else if ((*x).r > (*y).r) return(1);
   else return(0);

}
/*******************************************************************************/

