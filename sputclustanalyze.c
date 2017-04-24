
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*
  sputclustanalyze - code to analyze the cluster size of sputtered atoms.
  Code originally taken from bondinganalyze and modified from it. 
  Cluster analysis taken from clusteranalyze

  Neighbour list calculation originally taken from PARCAS analyzeint.c

  NOTE: All arrays here start from 0, not 1 as in parcas !!

*/

#define MAXNGBRS 80

typedef short logical;
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

#define MAXTYPE 12
int itypemin,itypemax;  /* Actual minimum and maximum */

#define True 1
#define False 0
#define pi 3.1415926535897

/*
  Function protypes 
*/


void clusteranalyze(real *x0,
		    int *itype,
		    real *box,
		    real *ngbrbox,
		    real *pbc,
		    real rcut,
		    real cutij[MAXTYPE][MAXTYPE],
		    int nat,
		    real surfpos,
		    real surfsign,
		    real time,
		    int ntime,  
		    real *v0); /* Krister */


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


char name[MAXTYPE][8];
int *in,*itype;

char string[80];

logical namesout; 
logical printxyz;
logical byevent;
int byeventtype;


main(argc,argv)
int argc;
char **argv;
{

   int iarg;
   int i,j,ii,nat;
   int ndiscard,nsputtriv;

   FILE *fp,*fpYtriv;
   char file[50];
   int nline;

   int stringcol,xdatacol,typedatacol;

   logical usetime,newtime,firstdatafound;
   logical iscoordline;
   int ntime,natexpectedp,natexpected;

   real xx,yy,zz;

   char buf[160],buforig[160],c1;
   char *argbuf,arg[40][40];
   int narg;

   real *x0, *v0;    /* Krister */
   real box[3],pbc[3];

   real casctime,prevtime;

   /* Cutoff for neighbour recognition */
   real rcut;
   real cutij[MAXTYPE][MAXTYPE];

   real newbox[3],ngbrbox[3];
      
   /* Sputtering limit and surface position and sign */
   real sputlim,sputlim0,surfpos,surfpos0,surfsign;

   if ((argc>1 && strcmp(argv[1],"-h")==0.0) || argc < 2) {
      printf("Usage: sputclustanalyze [options] file [stringcol [string [xdatacol [xsize [ysize [zsize]]]]]]\n");
      printf("\n"); 
      printf("Options: -debug           debug\n");
      printf("         -maxdebug        debug even more\n");
      printf("         -sputlim z       Sputtering limit for simple Y calc. zsize+40Å assumed\n");
      printf("         -surfpos z       Surface position\n");
      printf("         -surfsign +1/-1  Surface is positive or negative? Positive assumed\n");
      printf("         -typecol n       Set atom type column to n\n");
      printf("         -cut rcut        Set neighbour recognition cutoff, 5.55 assumed\n");
      printf("         -cutij i j r     Set pair-specific cutoffs\n");
      printf("         -ngbrbox x y z   Set forced box size for neighbour list calculation\n");
      printf("         -pbc 0/1 0/1 0/1 (Un)Set periodic boundaries: 1 1 0 assumed\n");
      printf("         -namesout        Output sputtered atom names in atomclusters.out\n");
      printf("         -printxyz        Output xyz file of atom information including index of cluster they belong to\n");
      printf("         -byevent         Output atomclusters.out separately for each time \n");
      printf("         -bysputonly i    As byevent but only if atoms of type i sputtered\n");
      printf("\n");
      printf("\n");
      printf("Data is taken as x y z beginning from column xdatacol, from lines in which\n");
      printf("stringcol contains string. If stringcol==0, all lines are assumed to have xyz data\n");
      printf("\n");
      printf("Analyses atom coordinates for sputtered atoms and their cluster distribution\n");
      printf("\n");
      printf("Logic is that if no atom in a cluster is closer than rcut from the surface, \n");
      printf("the cluster is sputtered. \n");
      printf("\n");
      printf("If string is not empty, takes the time step from any line not containing string but\n");
      printf("containing the string fs and analyses separately for each time. If string is empty\n");
      printf("any line with less than three arguments is interpreted to separate times.\n"); 
      printf("\n");
      printf("xsize, ysize and zsize give the box size\n");
      printf("Default: 1 Au 2 114.24 114.24 114.24\n");
      printf("\n");
      printf("If option -ngbrbox this box size is used in neighbour list calc\n");
      printf("at nonperiodic boundaries. Otherwise max and min. coordinate in that\n");
      printf("dimension is used.\n");
      return(0);
   }

   iarg=1;

   /* Set some basic variables which may be modified by command line arguments */

   rcut=5.55;    /* Exact value for Foiles EAM Au */
   for (i=0;i<MAXTYPE;i++)  for (j=0;j<MAXTYPE;j++) {
     cutij[i][j]=-1.0;
   }
 
   surfsign=+1.0;
   surfpos0=-1e30;
   sputlim0=-1e30;

   stringcol=1;
   sprintf(string,"XX"); 
   xdatacol=2;
   typedatacol=xdatacol+3;

   box[0]=10.0;
   box[1]=10.0;
   box[2]=10.0;

   /* Assume periodics in xy */
   pbc[0]=pbc[1]=1.0; pbc[2]=0.0;

   newbox[0]=1e30;
   ngbrbox[0]=1e30; ngbrbox[1]=1e30; ngbrbox[2]=1e30;
      
   namesout=False;
   printxyz=False;
   byevent=False;
   byeventtype=-1;

   while (strncmp(argv[iarg],"-",1)==0) {
     if (strcmp(argv[iarg],"-debug")==0) debug=True; 
     else if (strcmp(argv[iarg],"-maxdebug")==0) maxdebug=True; 
     else if (strcmp(argv[iarg],"-sputlim")==0) sscanf(argv[++iarg],"%lg",&sputlim0);
     else if (strcmp(argv[iarg],"-surfpos")==0) sscanf(argv[++iarg],"%lg",&surfpos0);
     else if (strcmp(argv[iarg],"-surfsign")==0) sscanf(argv[++iarg],"%lg",&surfsign);
     else if (strcmp(argv[iarg],"-typecol")==0) sscanf(argv[++iarg],"%d",&typedatacol);

     else if (strcmp(argv[iarg],"-cut")==0) {
       sscanf(argv[++iarg],"%lg",&rcut);
       printf("Set neigbour cutoff to %g\n",rcut);
     }
     else if (strcmp(argv[iarg],"-cutij")==0) {
       sscanf(argv[++iarg],"%d",&i);
       sscanf(argv[++iarg],"%d",&j);
       sscanf(argv[++iarg],"%lg",&xx);
       cutij[i][j]=xx;
       cutij[j][i]=xx;
       printf("Set cutoff for pair %d %d and %d %d to be %g\n",i,j,j,i,xx);
     }
     else if (strcmp(argv[iarg],"-ngbrbox")==0) {
       sscanf(argv[++iarg],"%lg",&(ngbrbox[0]));
       sscanf(argv[++iarg],"%lg",&(ngbrbox[1]));
       sscanf(argv[++iarg],"%lg",&(ngbrbox[2]));
       printf("Set neighbour box size: x %g y %g z %g\n",ngbrbox[0],ngbrbox[1],ngbrbox[2]);
     }
     else if (strcmp(argv[iarg],"-pbc")==0) {
       sscanf(argv[++iarg],"%lg",&(pbc[0]));
       sscanf(argv[++iarg],"%lg",&(pbc[1]));
       sscanf(argv[++iarg],"%lg",&(pbc[2]));
       printf("Set periodics: x %g y %g z %g\n",pbc[0],pbc[1],pbc[2]);
     }
     else if (strcmp(argv[iarg],"-namesout")==0) namesout=True; 
     else if (strcmp(argv[iarg],"-printxyz")==0) printxyz=True; 
     else if (strcmp(argv[iarg],"-byevent")==0) byevent=True; 
     else if (strcmp(argv[iarg],"-bysputonly")==0) {
       byevent=True; 
       sscanf(argv[++iarg],"%d",&byeventtype);
     }
     else {
       printf("sputclustanalyze ERROR: Unknown option %s. Exiting...\n",argv[iarg]);
       exit(0);
     }
     
      iarg++;
   }
   sscanf(argv[iarg++],"%s",file);
   if (strlen(file)==0) { printf("No filename given. Exiting\n"); return(0); }

   DEBUGSR("Args",iarg);
   if (iarg < argc) sscanf(argv[iarg++],"%d",&stringcol);
   if (iarg < argc) sscanf(argv[iarg++],"%s",string);
   if (iarg < argc) sscanf(argv[iarg++],"%d",&xdatacol);
   if (iarg < argc) sscanf(argv[iarg++],"%lg",&(box[0]));
   box[1]=box[0]; box[2]=box[0];
   if (iarg < argc) sscanf(argv[iarg++],"%lg",&(box[1]));
   if (iarg < argc) sscanf(argv[iarg++],"%lg",&(box[2]));


   printf("Program arguments %d %s %d %lg %lg %lg\n",
	  stringcol,string,xdatacol,box[0],box[1],box[2]);
   fflush(NULL);

   surfpos=surfpos0; if (surfpos0==-1e30) surfpos=surfsign*box[2]/2.0;
   sputlim=sputlim0; if (sputlim0==-1e30) sputlim=surfpos+surfsign*40.0;
   printf("Sputtering limit for simple Y calculation %g \n",sputlim);
   printf("Initial surface position %g and sign %g\n",surfpos,surfsign);

   for (i=0;i<MAXTYPE;i++)  for (j=0;j<MAXTYPE;j++) {
     if (cutij[i][j]>rcut) rcut=cutij[i][j];
   }

   printf("Using max cutoff %g\n",rcut);

   usetime=True; 
   if (stringcol==0) {
     usetime=False;
     printf("Stringcol 0, assuming all atoms are of the same time\n");
   }
   for (i=0;i<MAXTYPE;i++) strncpy(name[i],"NotSet",8);
   ntime=0;
   casctime=0.0;


   MAXAT=230000;

   i=sizeof(real)*MAXAT*3+sizeof(int)*MAXAT;
   i+=sizeof(int)*MAXAT;
   printf("Allocating %d bytes for %d atoms ... ",i,MAXAT);
   x0=(real *) malloc(sizeof(real)*3*MAXAT);
   v0=(real *) malloc(sizeof(real)*3*MAXAT);  /* Krister */
   in=(int *)  malloc(sizeof(int)*MAXAT);
   itype=(int *) malloc(sizeof(int)*MAXAT);
   printf("... done\n");


   if (strcmp(file,"_")==0) fp=stdin;
   else fp=fopen(file,"r");

   if (fp==NULL) {
      printf("File %s open failed, exiting !\n",file);
      return;
   }

   fpYtriv=fopen("Ytriv.tdat","w");
   
   DEBUGS("Files opened\n");



   nat=0;
   nline=0;
   ndiscard=0;
   nsputtriv=0;
   itypemin=MAXTYPE+1; itypemax=-1;
   firstdatafound=False;
   newbox[0]=box[0]; newbox[1]=box[1]; newbox[2]=box[2];
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
		  if (strcmp(arg[i],"fs") == 0) {
		    sscanf(arg[i-1],"%lg",&casctime);
		    DEBUGSRR("Found time on line",i,casctime);
		  }
	       }
	       /* Check whether time line has box size as well */
	       if (strcmp(arg[6],"boxsize")==0) {
		  if (newbox[0]!=1e30) {
		     box[0]=newbox[0];
		     box[1]=newbox[1];
		     box[2]=newbox[2];
		  }
		  sscanf(arg[7],"%lg",&(newbox[0]));
		  sscanf(arg[8],"%lg",&(newbox[1]));
		  sscanf(arg[9],"%lg",&(newbox[2]));
		  if (surfpos0==-1e30) surfpos=surfsign*newbox[2]/2.0;
		  sputlim=sputlim0; if (sputlim==-1e30) sputlim=surfpos+surfsign*40.0;
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

	    /* Krister modifications */
	    if (narg>=xdatacol+4) sscanf(arg[xdatacol+4],"%d",&(in[nat])); 

	    if (narg>=xdatacol+5) sscanf(arg[xdatacol+5],"%lg",&(v0[3*nat+0])); else v0[3*nat+0]=0; 
	    if (narg>=xdatacol+6) sscanf(arg[xdatacol+6],"%lg",&(v0[3*nat+1])); else v0[3*nat+1]=0; 
	    if (narg>=xdatacol+7) sscanf(arg[xdatacol+7],"%lg",&(v0[3*nat+2])); else v0[3*nat+2]=0; 

	    if (itype[nat]<itypemin) itypemin=itype[nat];
	    if (itype[nat]>itypemax) itypemax=itype[nat];
	    if (itype[nat]>=MAXTYPE) {
	      printf("MAXTYPE %d too small: increase to %d\n",MAXTYPE,itype[nat]+1);
	      DEBUGSRRRRR(buforig,nline,nat,typedatacol,x0[3*nat+0],in[nat]);
	      exit(0);
	    }
	    if (xdatacol>0) strncpy(name[itype[nat]],arg[xdatacol-1],8);

	    /* Trivial Y calc. */
	    if ( (surfsign<0.0 && x0[3*nat+2] < sputlim) || (surfsign>0.0 && x0[3*nat+2] > sputlim) ) nsputtriv++;

	    /* Always discard atoms inside surface */
	    if ( (surfsign<0.0 && x0[3*nat+2] > surfpos) || (surfsign>0.0 && x0[3*nat+2] < surfpos) ) {
		  nat--; 
		  ndiscard++;
	    }	      
	    
	    nat++;
 	    if (nat>=MAXAT) { 
	       MAXAT*=2;
	       i=sizeof(real)*MAXAT*3+sizeof(int)*MAXAT;
	       i+=sizeof(int)*MAXAT;
	       printf("Reallocating %d bytes for %d atoms ... ",i,MAXAT);
	       x0=(real *) realloc(x0,3*sizeof(real)*MAXAT);
	       v0=(real *) realloc(v0,3*sizeof(real)*MAXAT);  /* Krister */
	       in=(int *)  realloc(in,sizeof(int)*MAXAT);
	       itype=(int *) realloc(itype,sizeof(int)*MAXAT);
	       printf("... done\n");
	    }
	 }
      }

      if (firstdatafound && newtime) {
	 printf(" Read in %d atoms for time %d %lg fs\n",nat,ntime-1,prevtime);
	 printf("Discarded %d atoms inside surface at %g\n",ndiscard,surfpos);
	 printf("--- At %g fs sputtering Y determined for cutoff %g: %d --- \n",prevtime,sputlim,nsputtriv);
	 fprintf(fpYtriv,"%g %d\n",prevtime,nsputtriv); fflush(NULL);
	 ndiscard=0;
	 nsputtriv=0;

	 if (nat<1) { 
	    printf("No atoms to handle, looking for next step %d\n",nat);
	    continue;
	 }
	 DEBUGSRRR("Last atom",x0[3*(nat-1)],x0[3*(nat-1)+1],x0[3*(nat-1)+2]);
	 printf("\n");
	 printf("--------------------------------------------------------------------\n");
	 printf("--------- Going into defect analysis at time %10g ------------\n",prevtime);
	 printf("--------------------------------------------------------------------\n");
	 printf("\n");
	 
	 clusteranalyze(x0,itype,box,ngbrbox,pbc,rcut,cutij,nat,surfpos,surfsign,prevtime,ntime,v0);  

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

   printf("Discarded %d atoms inside surface at %g\n",ndiscard,surfpos);
   printf("--- At %g fs sputtering Y determined for cutoff %g: %d --- \n",casctime,sputlim,nsputtriv);
   fprintf(fpYtriv,"%g %d\n",prevtime,nsputtriv); fflush(NULL); fclose(fpYtriv);
   ndiscard=0;
   nsputtriv=0;

   if (nat<1) { 
      printf("Warning: no atoms to handle\n");
   }
   
   if (!usetime) casctime=1.0; 

   if (newbox[0]!=1e30) {
      box[0]=newbox[0];
      box[1]=newbox[1];
      box[2]=newbox[2];
   }

   printf("\n");
   printf("--------------------------------------------------------------------\n");
   printf("--------- Going into defect analysis at time %10g ------------\n",casctime);
   printf("--------------------------------------------------------------------\n");
   printf("\n");
   
   clusteranalyze(x0,itype,box,ngbrbox,pbc,rcut,cutij,nat,surfpos,surfsign,casctime,ntime,v0);  /* Krister */
   
   printf("\n");
   printf("--------------------------------------------------------------------\n");
   printf("---------- Done with defect analysis at time %10g ------------\n",casctime);
   printf("--------------------------------------------------------------------\n");
   printf("\n");
	 


} /* End of program */


/*

  Actual analysis of data.

  Taken from clusteranalyze.c which was originally obtained from
  clusteranalyze.atom awk script
  
*/

#define MAXCLUS 30000
#define MAXATINCLUS 3000000

void clusteranalyze(real *x0,
		    int *itype,
		    real *box,
		    real *ngbrbox,
		    real *pbc,
		    real rcut,
		    real cutij[MAXTYPE][MAXTYPE],
		    int nat,
		    real surfpos,
		    real surfsign,
		    real time, int ntime, real *v0  /* Krister */
		    )

{

  real rcluster,rclustersq;
  int *iatcluster;
  int iat,N;
  int n,nline,narg;

  int intflag=1,xdatacol=0; 
  int largeatomclustersize=6;
  static FILE *fplog,*fpclu,*fpclustat,*fpY,*fpat;
  static FILE *fpclusat;  /* Krister */


  /* Analysis params */
  int i,j,jj,k,l,ni,nj,n1,n2,i3,j3;
  int type0,type1;
  int ntype0,ntype1; 
  int natint,natreal;
  int *p;
  real xx,yy,zz;
  real xfirst,yfirst,zfirst;

  int natomclusters,ntrueatomclusters,*prev;
  static int ncluster[MAXCLUS];
  /* static int icluster[MAXCLUS][MAXATINCLUS]; */
  int *icluster[MAXCLUS]; /* Array of pointers for dynamic allocation */
  static int iclustersize[MAXCLUS];
  static int natcluster[MAXCLUS];
  static logical sputtered[MAXCLUS];
  static real xcluster[MAXCLUS],ycluster[MAXCLUS],zcluster[MAXCLUS];
  static real xclustermin[MAXCLUS],yclustermin[MAXCLUS],zclustermin[MAXCLUS];
  static real xclustermax[MAXCLUS],yclustermax[MAXCLUS],zclustermax[MAXCLUS];
  static int atomincluststat[MAXATINCLUS],nclustermax,inclustermax;
  int clustmp[MAXATINCLUS];
  int nisolated,ntype1isolated,ntype0isolated;
   
  int atinlargeatomclusters,nlargeatomclusters;

  int ntypeinclus[MAXTYPE];
  int nsputtype[MAXTYPE];

  real rijx,rijy,rijz,r,rsq;
   
  real sizex,sizey,sizez,halfx,halfy,halfz;
  static logical firsttime=True;

  /* Sputtering variables */
  int Y,Yclus,Yclus10,Yclus100;
  int Nad,Nadclus,Nadclus10,Nadclus100;

  /* Neighbour list variables */
  int nngot;
  real xmax,ymax,zmax;
  struct ngbr ngbrs[MAXNGBRS];

  char buf[80];

  printf("Starting clusteranalyze for %d atoms\n",nat);

  if (firsttime) {
    firsttime=False;
    if (!byevent) {
      fpclu=fopen("atomclusters.out","w");
    }
    fpclustat=fopen("sputatomincluststat","w");
    fpY=fopen("Yclus.tdat","w");
    fpclusat=fopen("clusatoms.out", "w");  /* Krister */
    
    if (printxyz) fpat=fopen("sputclustat.xyz","w");


  }


  for (i=0;i<MAXCLUS;i++) {
    iclustersize[i]=100;
    icluster[i]=(int *) malloc(sizeof(int)*iclustersize[i]);
  }

  /* Open log file, this are separate for each time step, only last one remains */
  fplog=fopen("clusteranalyze.log","w");

  /* Since MAXAT may change, these have to be updated every time */
  iatcluster=(int *) malloc(MAXAT*sizeof(int));


  /* Interface renaming of some variables */
  rcluster=rcut; rclustersq=rcluster*rcluster;
  N=nat;
  sizex=box[0]; halfx=box[0]/2;
  sizey=box[1]; halfy=box[1]/2;
  sizez=box[2]; halfz=box[2]/2;

   /* Initialize neighbour list calc. */
  /* Find min and max coordinates to set ngbrbox, if appropriate */
  xmax=ymax=zmax=0.0;
  for (i=0;i<N;i++) {
    i3=i*3;
    if (fabs(x0[i3+0]) > xmax) xmax=fabs(x0[i3+0]);
    if (fabs(x0[i3+1]) > ymax) ymax=fabs(x0[i3+1]);
    if (fabs(x0[i3+2]) > zmax) zmax=fabs(x0[i3+2]);
  }
  if (pbc[0]==1.0) {
    ngbrbox[0]=box[0];
  }
  else {
     if (ngbrbox[0]==1e30) ngbrbox[0]=xmax+3*rcut;
  }
  if (pbc[1]==1.0) {
    ngbrbox[1]=box[1];
  }
  else {
     if (ngbrbox[1]==1e30) ngbrbox[1]=ymax+3*rcut;
  }
  if (pbc[2]==1.0) {
    ngbrbox[2]=box[2];
  }
  else {
     if (ngbrbox[2]==1e30) ngbrbox[2]=zmax+3*rcut;
  }

  printf("Using box size %g %g %g for neighbour list construction\n",ngbrbox[0],ngbrbox[1],ngbrbox[2]);
  getngbrs(-1,x0,ngbrbox,pbc,N,rcut,ngbrs,&nngot);

  /* Make cluster analysis; this is a ''bit'' tricky ! */

  /* Format of cluster arrays: */
  /* ncluster[i]              Number of atoms that belong to cluster */
  /* natcluster[i]            should be identical to ncluster in this implementation */
  /* icluster[i][ncluster[i]]  Indices of those atoms */

  natomclusters=0;
  /* Initialize arrays */
  for (i=0;i<N;i++) {   
    iatcluster[i]=-1; 
  }

  /*********************** START OF CLUSTER CREATION *****************************/

  /* FORM atomclusters by going through all atoms */
  /* and placing adjacent ones to atomclusters */
  for (i=0;i<N;i++) {   /* Outer loop over atoms */
    if (i%1000==0) { fprintf(stderr,"%1d",(i%10000)/1000); fflush(stderr); }

    /* Find neighbours j of i utilizing linked list subroutine */
    /* Note that this does not account for cutij, that is done separately below */
    getngbrs(i,x0,ngbrbox,pbc,N,rcut,ngbrs,&nngot);
    MAXDEBUGSRRR("Got neighbours",i,rcut,nngot);

    for (jj=0;jj<nngot;jj++) {   /* Inner loop over atoms */
      j=ngbrs[jj].i;
      if (i==j) continue;
      if (iatcluster[i]>natomclusters) { 
	fprintf(stderr,"CRAZY i %d %d %d\n",i,iatcluster[i],natomclusters); exit(0); 
      }
      if (iatcluster[j]>natomclusters) { 
	fprintf(stderr,"CRAZY j %d %d %d\n",j,iatcluster[j],natomclusters); exit(0); 
      }

      MAXDEBUGSRRRR("i j r",i,j,rsq,ngbrs[jj].r*ngbrs[jj].r);
      r=ngbrs[jj].r;
      if (r<rcluster) {

	if (cutij[itype[i]][itype[j]] > 0.0) {
	  if (r > cutij[itype[i]][itype[j]]) {
	    /* printf("Bond %d %d r %g ignored because of cutij %g\n",
	       i,j,r,cutij[itype[i]][itype[j]]); */
	    continue;
	  }
	}

	MAXDEBUGSRRR("success i j r",i,j,r);

	if (iatcluster[i] == -1 && iatcluster[j] == -1) {
	  /* Neither is in cluster => create new cluster */
	  MAXDEBUGSRRR("creating cluster",i,j,natomclusters);
	  n=natomclusters;
	  fprintf(fplog,"Creating cluster %d\n",n);
	  if (n>=MAXCLUS) { fprintf(stderr,"Increase MAXCLUS in clusteranalyze\n"); exit(0); }
	  natcluster[n]=2;
	  ncluster[n]=2;
	  icluster[n][0]=i; icluster[n][1]=j;
	  iatcluster[i]=n; iatcluster[j]=n;
	  natomclusters++;
	}
	else if (iatcluster[i] != -1 && iatcluster[j] == -1) {
	  MAXDEBUGSRRR("placing j in i's cluster",i,j,iatcluster[i]);
	  /* i is in cluster, j not => place j in i:s cluster */
	  n=iatcluster[i];
	  /* print "Adding j",j," to cluster",n  >> "clusteranalyze.log"; */
	  if (ncluster[n]>=iclustersize[n]) { 
	    iclustersize[n]*=10;
	    fprintf(fplog,"Increasing icluster element %d size to %d i\n",n,iclustersize[n]); fflush(fplog);
	    DEBUGSRR("icluster realloc",n,iclustersize[n]);
	    icluster[n]=(int *) realloc(icluster[n],(size_t) iclustersize[n]*sizeof(int));
	  }
	  icluster[n][ncluster[n]]=j;
	  natcluster[n]+=1;
	  ncluster[n]++;
	  iatcluster[j]=n;
	}
	else if (iatcluster[i] == -1 && iatcluster[j] != -1) {
	  MAXDEBUGSRRR("placing i in j's cluster",i,j,iatcluster[j]);
	  /* j is in cluster, i not => place i in j:s cluster */
	  n=iatcluster[j];
	  /* print "Adding i",i," to cluster",n  >> "clusteranalyze.log"; */
	  if (ncluster[n]>=iclustersize[n]) { 
	    iclustersize[n]*=10;
	    fprintf(fplog,"Increasing icluster element %d size to %d j\n",n,iclustersize[n]); 
	    DEBUGSRR("icluster realloc",n,iclustersize[n]);
	    icluster[n]=(int *) realloc(icluster[n],(size_t) iclustersize[n]*sizeof(int));
	  }
	  icluster[n][ncluster[n]]=i;
	  natcluster[n]+=1;
	  ncluster[n]++;
	  iatcluster[i]=n;
	}
	else if (iatcluster[i] != -1 && iatcluster[j] != -1) {
	  MAXDEBUGSRRRR("both in clusters",i,j,iatcluster[i],iatcluster[j]);
	  /* Both are in atomclusters */
	  /* If its the same cluster, do nothing */
	  ni=iatcluster[i];
	  nj=iatcluster[j];
	  /* print "Atoms",i,j,"both in atomclusters", */
	  /*      ni,nj >> "clusteranalyze.log"; */
	  if (ni==nj) continue;
               
	  /* Atomclusters are not the same =>  */
	  /* MERGE atomclusters */
	  /* Pick lower one to be the primary */
	  if (ni < nj) { n1=ni; n2=nj; }
	  else { n1=nj; n2=ni; }

	  fprintf(fplog," ... merging %d %d into %d %d\n",n2,ncluster[n2],n1,ncluster[n1]);
	  /* Copy cluster n2 into n1 */
	  for(k=ncluster[n1];k<ncluster[n1]+ncluster[n2];k++) {
	    if (k>=iclustersize[n1]) { 
	      iclustersize[n1]*=10;
	      fprintf(fplog,"Increasing icluster element %d size to %d merge\n",n,iclustersize[n1]); 
	      DEBUGSRR("icluster realloc",n1,iclustersize[n1]);
	      icluster[n1]=(int *) realloc(icluster[n1],(size_t) iclustersize[n1]*sizeof(int));
	    }
	    icluster[n1][k]=icluster[n2][k-ncluster[n1]];
	    
	    /* if (k>=MAXATINCLUS) { fprintf(stderr,"Increase MAXATINCLUS (merge)\n"); exit(0); } */
	  }
	  natcluster[n1]+=natcluster[n2];
	  ncluster[n1]+=ncluster[n2];
               
	  /* Delete cluster n2 by overwriting with upper atomclusters */
	  /* Take care that lower elements are not too small! */
	  natomclusters--;
	  for(k=n2;k<natomclusters;k++) {
	    /* Copy cluster k+1 to k */
	    ncluster[k]=ncluster[k+1];
	    natcluster[k]=natcluster[k+1];
	    MAXDEBUGSRRR("icluster k copy",k,ncluster[k],ncluster[k+1]);
	    /* Note that this has to be a while, not just an if!! */
	    while (ncluster[k]>=iclustersize[k]) { 
	      DEBUGSRRR("icluster k copy",k,ncluster[k],ncluster[k+1]);
	      iclustersize[k]*=10;
	      fprintf(fplog,"Increasing icluster element %d size to %d del\n",n,iclustersize[k]); 
	      DEBUGSRR("icluster realloc",k,iclustersize[k]);
	      icluster[k]=(int *) realloc(icluster[k],(size_t) iclustersize[k]*sizeof(int));
	      DEBUGSRR("icluster realloc done",k,iclustersize[k]);
	    }
	    for (l=0;l<ncluster[k];l++) {   
	      icluster[k][l]=icluster[k+1][l];
	    }
	  }
	  /* Go trough all atoms and rename cluster indices */
	  for (k=0;k<N;k++) {   
	    if (iatcluster[k]==n2) iatcluster[k]=n1;
	    if (iatcluster[k]>n2) iatcluster[k]--;
	  }

	} /* End of placing atom into correct cluster */

      }
    }
  }

  ntrueatomclusters=natomclusters;
  fprintf(stdout,"\nCreated %d atomclusters with multiple atoms\n",natomclusters);

  /* Create individual atomclusters for non-clustered atoms */
  nisolated=0;
  ntype1isolated=0;
  ntype0isolated=0;
  for (i=0;i<N;i++) {
    if (iatcluster[i]==-1) {
      nisolated++;
      n=natomclusters;
      if (n>=MAXCLUS) { fprintf(stderr,"Increase MAXCLUS in clusteranalyze (isol)\n"); exit(0); }
      iatcluster[i]=n;
      ncluster[n]=1;
      natcluster[n]=1;
      icluster[n][0]=i;         
      natomclusters++;
    }
  }
  fprintf(stdout,"Added %d atomclusters with only one atom\n",natomclusters-ntrueatomclusters);

  /*********************** END OF CLUSTER CREATION *****************************/

  /********** Rest is analysis and cleanup stuff **********/

  /* Get center, min and max of each cluster  */

  if (surfsign < 0.0) printf("Interpreting clusters completely below %g as sputtered\n",surfpos-rcut);
  if (surfsign > 0.0) printf("Interpreting clusters completely above %g as sputtered\n",surfpos+rcut);
  atinlargeatomclusters=0;
  nlargeatomclusters=0;
  Y=0; Yclus=Yclus10=Yclus100=0;
  Nad
=0; Nadclus=Nadclus10=Nadclus100=0;
  for (i=0;i<natomclusters;i++) {
    MAXDEBUGSR("Cluster pos anal",i);
    sputtered[i]=False;
    xcluster[i]=0.0; ycluster[i]=0.0; zcluster[i]=0.0; 
    xclustermin[i]=x0[icluster[i][0]*3+0]; yclustermin[i]=x0[icluster[i][0]*3+1]; zclustermin[i]=x0[icluster[i][0]*3+2]; 
    xclustermax[i]=x0[icluster[i][0]*3+0]; yclustermax[i]=x0[icluster[i][0]*3+1]; zclustermax[i]=x0[icluster[i][0]*3+2]; 
    if (ncluster[i]==0) { fprintf(stderr,"Zero cluster",i); continue; }
    for (j=0;j<ncluster[i];j++) {
      /* 
	 26.8 2002: Handling periodics here is a bit tricky. Hence do the following:
	 take the first atom in each cluster as a reference point and ensure the rest of
	 the atoms are on the same side 
      */
      xx=x0[icluster[i][j]*3+0];
      if (j==0) { xfirst=xx; }
      if (pbc[0]==1.0) { rijx=xx-xfirst; if (rijx>halfx) xx-=sizex; else if (rijx<-halfx) xx+=sizex; }
      yy=x0[icluster[i][j]*3+1];
      if (j==0) { yfirst=yy; }
      if (pbc[1]==1.0) { rijy=yy-yfirst; if (rijy>halfy) yy-=sizey; else if (rijy<-halfy) yy+=sizey; }
      zz=x0[icluster[i][j]*3+2];
      if (j==0) { zfirst=zz; }
      if (pbc[2]==1.0) { rijz=zz-zfirst; if (rijz>halfz) zz-=sizez; else if (rijz<-halfz) zz+=sizez; }
      if (xx < xclustermin[i]) xclustermin[i]=xx; if (xx > xclustermax[i]) xclustermax[i]=xx;
      if (yy < yclustermin[i]) yclustermin[i]=yy; if (yy > yclustermax[i]) yclustermax[i]=yy;
      if (zz < zclustermin[i]) zclustermin[i]=zz; if (zz > zclustermax[i]) zclustermax[i]=zz;
      xcluster[i]+=xx;
      ycluster[i]+=yy;
      zcluster[i]+=zz;
    }
    xcluster[i]/=ncluster[i];       
    if (pbc[0]==1.0) { if (xcluster[i]>halfx) xcluster[i]-=sizex; else if (xcluster[i]<-halfx) xcluster[i]+=sizex; }
    ycluster[i]/=ncluster[i];       
    if (pbc[1]==1.0) { if (ycluster[i]>halfy) ycluster[i]-=sizey; else if (ycluster[i]<-halfy) ycluster[i]+=sizey; }
    zcluster[i]/=ncluster[i];       
    if (pbc[2]==1.0) { if (zcluster[i]>halfz) zcluster[i]-=sizez; else if (zcluster[i]<-halfz) zcluster[i]+=sizez; }
    
    if (natcluster[i]>=largeatomclustersize) {
      nlargeatomclusters++;
      atinlargeatomclusters+=natcluster[i];        
    }  

    /*
      Now comes the crucial part for sputtering: if the cluster is within rcut from surface.
      it is not sputtered: otherwise it is! 
    */
    if (surfsign < 0.0 && zclustermax[i] < surfpos-rcut) sputtered[i]=True; 
    if (surfsign > 0.0 && zclustermin[i] > surfpos+rcut) sputtered[i]=True; 
    if (sputtered[i]) {
      Y+=natcluster[i];
      if (natcluster[i]>1) Yclus+=natcluster[i];
      if (natcluster[i]>10) Yclus10+=natcluster[i];
      if (natcluster[i]>100) Yclus100+=natcluster[i];
    }
    else {
      /*
	These guys are not clustered, but are above surface,
	so they should be adatoms!
      */
      Nad+=natcluster[i];
      if (natcluster[i]>1) Nadclus+=natcluster[i];
      if (natcluster[i]>10) Nadclus10+=natcluster[i];
      if (natcluster[i]>100) Nadclus100+=natcluster[i];
    }
  }   


  
  DEBUGS("Writing out atoms and their clustering information to clusatoms.out\n");  /* Krister */

  for (i=0;i<natomclusters;i++) {  /* Krister */
    for (j=0; j<ncluster[i]; j++) {  /* Krister */
      /* icluster[i][j] = atom index for atom number j in cluster number i */   /* Krister */
      fprintf(fpclusat, "%.10f %.10f %.10f %12g %12g %12g %12g %d %d %d %d",    /* Krister */
	      x0[3*icluster[i][j]+0],   /* Krister */
	      x0[3*icluster[i][j]+1],   /* Krister */
	      x0[3*icluster[i][j]+2],   /* Krister */
	      v0[3*icluster[i][j]+0],   /* Krister */
	      v0[3*icluster[i][j]+1],   /* Krister */
	      v0[3*icluster[i][j]+2],   /* Krister */
	      time,   /* Krister */
	      itype[icluster[i][j]],   /* Krister */
	      in[icluster[i][j]],   /* Krister */
	      i,   /* Krister */
	      ncluster[i]);   /* Krister */
      if (sputtered[i]>0) fprintf(fpclusat," sput\n");
      else fprintf(fpclusat," nonsput\n");
    }   /* Krister */
  }   /* Krister */


  if (printxyz) {
    DEBUGS("Writing out atoms and their clustering information to sputclustat.xyz\n");
    fprintf(fpat," %d\n",N);
    fprintf(fpat,"sputclustanalyze output at time %g fs boxsize %g %g %g surfpos %g\n",time,box[0],box[1],box[2],surfpos);
    for (i=0;i<N;i++) {
      fprintf(fpat,"%s %.12g %.12g %.12g %d %d %d\n",name[itype[i]],x0[3*i+0],x0[3*i+1],x0[3*i+2],iatcluster[i],itype[i],in[i]);
    }
  }


  /* Do sputtering statistics by atom type */
  for (j=0;j<MAXTYPE;j++) nsputtype[j]=0;
  for (i=0;i<natomclusters;i++) {
    if (sputtered[i]>0) {
      for (j=0; j<ncluster[i]; j++) {  
	nsputtype[itype[icluster[i][j]]]++;
      }
    }
  }
  printf("Y per atom type: ");
  for (j=0;j<MAXTYPE;j++) {
    if (nsputtype[j]>0) printf(" %s: %d",name[j],nsputtype[j]);
  }
  printf("\n");

  if (byeventtype>=0 && nsputtype[byeventtype]==0)  {
    printf("No sputtered atoms of type %d: no atomclusters output\n",byeventtype);
  }
  else {

    if (byevent) {
      sprintf(buf,"atomclusters.out.%d",ntime);
      fpclu=fopen(buf,"w");
    }
 
    DEBUGS("Writing out cluster data to atomclusters.out\n");
    /* File written in 'vacancies.out/liquidat.out' format */
    for (i=0;i<natomclusters;i++) {
      fprintf(fpclu,"%g %g %g %g %d zrange %g %g # %d",
	      xcluster[i],ycluster[i],zcluster[i],time,
	      natcluster[i],zclustermin[i],zclustermax[i],i);
      if (namesout) {
	for (j=0;j<MAXTYPE;j++) ntypeinclus[j]=0;
	for (j=0; j<ncluster[i]; j++) {  
	  fprintf(fpclu," %2s",name[itype[icluster[i][j]]]);
	  ntypeinclus[itype[icluster[i][j]]]++;
	}      
      }
      if (sputtered[i]>0) fprintf(fpclu," sput");
      else fprintf(fpclu," nonsput\n");
      if (namesout) {
	fprintf(fpclu," ");
	for (j=0;j<MAXTYPE;j++) {
	  if (ntypeinclus[j]>0) {
	    fprintf(fpclu,"%s%d",name[j],ntypeinclus[j]);
	  }
	}
      }
      fprintf(fpclu,"\n");
    }
    if (byevent) fclose(fpclu);

  }


  /* Do statistics of number of atoms/cluster */
  /* for sputtered clusters! */
  if (N+2 > MAXATINCLUS) { fprintf(stderr,"Increase MAXATINCLUS to %d\n",N+2); exit(0); }
  for(i=0;i<N+2;i++) atomincluststat[i]=0;
  nclustermax=0; inclustermax=0;
  for (i=0;i<natomclusters;i++) {
    if (sputtered[i]) {
      if (ncluster[i] > MAXATINCLUS) { fprintf(stderr,"Increase MAXATINCLUS %d\n",ncluster[i]); exit(0); }
      atomincluststat[ncluster[i]]++;
      if (natcluster[i]>nclustermax) {
	nclustermax=natcluster[i];
	inclustermax=i;
      }
    }
  }
  printf("Largest sputtered cluster has %d atoms centered at %g %g %g\n",
	 nclustermax,xcluster[inclustermax],ycluster[inclustermax],zcluster[inclustermax]);
  
  /* printf "Writing out cluster size statistics to atomincluststat\n"; */
  for (i=1;i<nclustermax+2;i++) {
    fprintf(fpclustat,"%12g %d %d %d\n",time,i,atomincluststat[i],i*atomincluststat[i]);
  }
  
  printf("At %g fs Y for different sizes: Tot %d >1 %d >10 %d >100 %d\n",time,Y,Yclus,Yclus10,Yclus100);
  printf("At %g fs Nadatoms for different sizes: Tot %d >1 %d >10 %d >100 %d\n",time,Nad,Nadclus,Nadclus10,Nadclus100);
  printf("\n--- At %g fs number of sputtered atoms from clustanalyze: %d ---\n",time,Y,N-Y);
  
  fprintf(fpY,"%g %d\n",time,Y);

  free(iatcluster);

  for (i=0;i<MAXCLUS;i++) free(icluster[i]);

  getngbrs(-2,x0,ngbrbox,pbc,N,rcut,ngbrs,&nngot);
  
  fclose(fplog);

  /* fclose(fpclustat); */

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
* In this implementation, should be called with -1 at beginning of every time,
* with -2 at end to free memory!
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

#define round(x) ((int) (x>0 ? x+0.5 : x-0.5))

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

   if (iat==-2) {
     free(head);
     free(list);
     return;
   }
       

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

   if (iat==-1) {

     printf("\ngetngbr: allocating %lu bytes for linkcell arrays\n",
	    sizeof(int)*xnc*ync*znc+sizeof(int)*(N+1));
     
     printf("Cell ranges %d - %d, %d - %d, %d - %d\n",
	    Nixmin,Nixmax,Niymin,Niymax,Nizmin,Nizmax);
     printf("Cell size %g %g %g\n\n",cellsize[0],cellsize[1],cellsize[2]);
     fflush(stdout);
     
     head=(int *) malloc(sizeof(int)*xnc*ync*znc);
     list=(int *) malloc(sizeof(int)*(N+1));
     if (head==NULL || list==NULL) {
       fprintf(stderr,"head/list allocation error\n");
       printf("Tried to allocate head: %ld list: %ld bytes\n",sizeof(int)*xnc*ync*znc,sizeof(int)*(N+1));
       exit(0);
     }

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
	 list[i]=head[icell];
	 head[icell]=i;
	 /* printf("i %d ix iy iz %d %d %d icell %d head[icell] %d\n",i,ix,iy,iz,icell,head[icell]);  */

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
	    /* MAXDEBUGSRRRRR("           head",nind(inx,iny,inz),j,inx,iny,inz); */
	    nj=0;
	    while (True) {
	       nj++;
	       if (nj>1) j=list[j];
	       /* MAXDEBUGSRRRRR("           list",iat,j,inx,iny,inz); */
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


   if (iat < 0) {
     DEBUGSRRR("Found neighbours",iat,nn,nj);
   }
   
   MAXDEBUGSRRR("Found neighbours",iat,nn,nj);
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
