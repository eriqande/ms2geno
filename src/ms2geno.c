/* ms2geno.c   A program to take output from Richard Hudson's program 
               makesamples, and then produce microsatellite or SNP
			   genotype data from it.  Currently we are just going to
			   deal with independent loci, but we could add some
			   more interesting stuff for doing linked microsats, etc.,
			   later.
			   
Much of the code for reading in the ms output was copied directly from
Richard Hudson's "sample_stats.c" source file.

*/
#define UN_EXTERN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ECA_MemAlloc.h"
#include "MathStatRand.h"
#include "ranlib.h"
#include "ECA_Opt3.h"




#define MAX_LINE_LENGTH 10000
#define MAX_NUM_POPS 10000
#define MAX_INPUT_FILENAME_LENGTH 5000


typedef enum {
	MSAT,
	SNP
} marker_types;


/* here is a data struct to hold all the options and things */
typedef struct {

	marker_types MarkerType;
	int NumDataSets;
	
	int Ploidy;  /* the ploidy of the simulated organisms */
	int *NumBase; /* number of individuals for each baseline data set */
	int TotBase;  /* total number for baseline */

	int *NumMix; /* number of individuals to put in each mixture */
	int TotMix;
	
	int NumLoc;  /* number of loci to be included in the data set */
	int *NumAlle; /* to report how many alleles were found amongst the microsat samples */
	
	int NumPops;  /* the number of populations that we are dealing with. */
	
	double MultiMutRate; /* the rate of multistep mutations */
	double MultiStepSize; /* Multistep mutations are 2 + Poisson with mean MultiStepSize. */
	
	
	char InputFileName[MAX_INPUT_FILENAME_LENGTH];  /* in case we want to take input from a file, and not from stdin */
	
	
	
	int *BaseStarts; /* to hold the indexes of the individuals that start off the baseline indivs for each pop */
	int *MixStarts; /* to hold the indexes of the individuals that start off the mixture indivs for each pop */
	int ***Basegenos; /* to store the genotypes of individuals in the baseline indexed by pop,indiv,gene-copy */
	int ***Mixgenos;  /* to store the genotypes of individuals in the mixture indexed by pop,indiv,gene-copy */
	char ***IDs;  /* to store the IDs of individuals indexed by pop,indiv,character*/
	int *MinGeno; /* to store the min microsatellite length (which will be negative---this way we can make them all positive) */
	
	
	/* some snp ascertainment stuff */
	int NumAscPops;  /* number of populations with members in the ascertainment set */
	int **AscSets;  /* array of ascertainment populations.  For the i-th population in the ascertainment set (starting from 0)
						AscSets[i][0] is the index (starting from zero) of that population.  AscSets[i][1] is the number of
						individuals from that population in the ascertainment set, and AscSet[i][2] is the minimum number of copies of the 
						minor allele amongst the individuals from this population required to trigger ascertainment of the SNP. */
	int AllPopsAsc;  /* the number of minor allele copies among all the populations in the ascertainment sets that would trigger an ascertainment event */

	int AllPopsGenoAsc;  /* flag.  If 1 then it means to count up the SNP genotypes in all pops and ascertain if you saw at least one of each of
						the three possible genotypes */
	int AllPopsPseudoAFLP; /* flag.  if 1 then we are ascertaining SNPs for pseudo-aflp's */
	int MinAFLP_Band;
	int MinAFLP_NoBand;
	
	
	
	int NoPerm;  /* this should be set nonzero if you DO NOT want to order the ascertainment pops randomly for each SNP */
	
	
	/* some stuff for ms_pop_override */
	int DoOverride;
	int OverRideNumPops;
	int *OverRideNumsInPops;
	
	
} msgeno_opts;



/* function prototypes */
void PrintDataSet(msgeno_opts *m, int DataSetNum);
msgeno_opts *GetMSGenoOpts(int argc, char *argv[]);
void CheckPopNums(msgeno_opts *m, int npops, int *popsizes);
int AddLocus(msgeno_opts *m, char **list, int segsites, int Loc);
int ListToInt(char c);
int SimLoc(msgeno_opts *m, char **list, int segsites, int Loc);
int FindTheSNP(int segsites, char **list, msgeno_opts *m);

double nucdiv(int, int, char **);
double tajd(int, int, double) ;
double hfay(int, int, char **);
double thetah(int, int, char **);

int maxsites = 100 ; /* start with a nice small number */

int main(int argc,char *argv[])
{
	int nsam, j ,nsites, i,  howmany  ;
	int npops=0, *popsizes=NULL, scan_count, check;  /* number of populations sampled from and number sampled from each */
	char junk[MAX_LINE_LENGTH];
	char **list, **cmatrix(), allele,na, line[MAX_LINE_LENGTH + 1], *lptr ;
	FILE *pf, *fopen(), *pfin ;
	double *posit   ;
	int   segsites, count  , nadv, probflag  ;
	double pi , h, th  ,prob ;
	char dum[20], astr[100] ;
	int  nsegsub, segsub( int nsam, int segsites, char **list ) ;
	msgeno_opts *Opts;
	int NumDataSets = 0, Loc=0;
	
	
	Opts = GetMSGenoOpts(argc,argv);

	/* do the seeds for my stuff */
	SeedFromFile("ms2geno_seeds");

	/* check to see if the user wants to take input from a file.  If so, open it. */
	if(strlen(Opts->InputFileName)>0) {
		if( (pfin=fopen(Opts->InputFileName,"r"))==NULL) {
			fprintf(stderr,"Error! Failed to open file \"%s\".  Please check that the path is correct, and maybe ensure there are no spaces in the pathname.  Exiting...\n");
			exit(1);
		}
	}
	else {
		pfin = stdin ;
	}

/* read in first two lines of output  (parameters and seed) */
  fgets( line, MAX_LINE_LENGTH, pfin);
  sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
  
  
  
  /* now, here ECA added some stuff to parse out the number of sizes of the separate subpopulations */
  lptr = line;
  sprintf(lptr,"%s  GET_OUT_NOW",lptr);  /* put a tag at the end of the command line to know when you've gotten to the end */
  
  /* WE HAVE THE OPTION OF OVERRIDING THIS, IF WE USE THE OVERRIDE OPTION.  THIS IS USEFUL
     WHEN YOU HAVE SIMULATED PANMIXIA AND WANT TO FOOL MS2GENO INTO THIKNING IN TERMS OF
	 SEPARATE POPULATIONS */
  if(Opts->DoOverride) { int i; int tot;
		npops = Opts->OverRideNumPops;
		popsizes = Opts->OverRideNumsInPops;
		
		/* now check to make sure we are OK */
		tot=0;
		for(i=0;i<npops;i++) {
			tot+=popsizes[i];
		}
		if(tot>nsam) {
			fprintf(stderr,"Error.  You issued the --ms-pop-override option; however, the sum of the requested number of chromosomes per population is %d which exceeds the total number of chromosomes (nsam) simulated by ms: %d.   Exiting...\n",
						tot,nsam);
			exit(1);	
		}
  }
  else {  /* IF WE ARE NOT OVERRIDING MS POP NUMS */
	  /* iteratively grab strings off the command line till you hit "-I" or the end tag */
	  do {
		check = sscanf(lptr," %s %n",junk,&scan_count);
		lptr += scan_count;
	  } while (strcmp(junk,"-I")!=0 && strcmp(junk,"GET_OUT_NOW")!=0  );
	  
	  /* now, if you hit a -I, get the number of pops and then the sizes of them */
	  if(strcmp(junk,"-I")==0) {
		sscanf(lptr," %d %n",&npops,&scan_count);
		lptr += scan_count;
		popsizes = (int *)calloc(npops,sizeof(int));
		for(i=0;i<npops;i++) {
			sscanf(lptr," %d %n",&(popsizes[i]),&scan_count);
			lptr += scan_count;
		}
	  }
	  else {  /* if you didn't hit a -I, then there is only one pop */
		npops=1;
		popsizes = (int *)calloc(npops,sizeof(int));
		popsizes[0] = nsam;
	  }
  }
  
  /* just some output while testing */
/*  printf("\nnpops %d:  ",npops);
  for(i=0;i<npops;i++) {
	printf("%d ",popsizes[i]);
  }
  printf("\n");
*/  
  
  CheckPopNums(Opts,npops,popsizes); 
  
  
  
  fgets( line, MAX_LINE_LENGTH, pfin);

	if( argc > 1 ) { 
	   nadv = atoi( argv[1] ) ; 
	}

  list = cmatrix(nsam,maxsites+1);
  posit = (double *)malloc( maxsites*sizeof( double ) ) ;

  count=0;
	probflag = 0 ;
while( howmany-count++ ) {

/* read in a sample */
  do {
     fgets( line, MAX_LINE_LENGTH, pfin);
  }while ( line[0] != '/' );
 
  fscanf(pfin,"  segsites: %d", &segsites );
  if( segsites >= maxsites){
	maxsites = segsites + 10 ;
	posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
        biggerlist(nsam,maxsites, list) ;
        }
   if( segsites > 0) {
		fscanf(pfin," %s", astr);
		if( astr[1] == 'r' ){
		   fscanf(pfin," %lf", &prob ) ;
		   probflag = 1;
		   fscanf(pfin," %*s");
		}
		for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
		for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
		
		
		/* now, down here we can stick these into genotype data sets */
		if(AddLocus(Opts,list,segsites,Loc++)==-1)  {
			Loc--;  /* if it was  a monomorphic locus, then re-do it */
		}
		
		/* if we have enough loci for a complete data set, print it and move on to the next */
		if(Loc == Opts->NumLoc) {
			Loc=0;
			NumDataSets++;
			PrintDataSet(Opts, NumDataSets);
			printf("CompletedDataSetNumber: %d   NumAlleles: ",NumDataSets);
			for(i=0;i<Opts->NumLoc;i++)  {
				printf(" %d",Opts->NumAlle[i]);
			}
			printf("\n");
		}
	    
	}
  }
  
  SeedToFile("ms2geno_seeds");
  return(0);
}	
	
	
/* analyse sample ( do stuff with segsites and list)  */
/*	if( argc > 1 ) nsegsub = segsub( nadv, segsites, list) ;
	pi = nucdiv(nsam, segsites, list) ;
	h = hfay(nsam, segsites, list) ;
	th = thetah(nsam, segsites, list) ;
	if( argc > 1 )
	printf("pi: %lf ss: %d  D: %lf H: %lf thetah: %lf segsub: %d \n", pi,segsites, tajd(nsam,segsites,pi) , h , th, nsegsub ) ;
	else if( probflag == 1 ) 
	  printf("pi:\t%lf\tss:\t%d\tD:\t%lf\tthetaH:\t%lf\tH:\t%lf\tprob:\t%g\n", pi,segsites, tajd(nsam,segsites,pi) , th , h, prob  ) ;
	else 
	  printf("pi:\t%lf\tss:\t%d\tD:\t%lf\tthetaH:\t%lf\tH:\t%lf\n", pi,segsites, tajd(nsam,segsites,pi) , th , h  ) ;
	

  }
*/



/* simple function to make sure there are enough chromosomes simulated using ms to 
fulfill the desired data sets.  This also is used to set BaseStarts and MixStarts.

 */
void CheckPopNums(msgeno_opts *m, int npops, int *popsizes) 
{
	int i,temp, mix, base;
	
	if(npops != m->NumPops) {
		fprintf(stderr,"Error. Command line number of pops is %d, but in output from ms it is %d.  They must match. Exiting.\n",m->NumPops,npops);
		exit(1);
	}
	
	m->BaseStarts = (int *)calloc(m->NumPops,sizeof(int));
	m->MixStarts = (int *)calloc(m->NumPops,sizeof(int));
	m->TotMix = 0;
	m->TotBase = 0;

	for(i=0;i<m->NumPops;i++)  {
		temp = mix = base = 0;
		if(m->NumBase) {
			base = m->NumBase[i];
			temp += base;
			m->TotBase += base;
		}
		if(m->NumMix) {
			mix += m->NumMix[i];
			temp += mix;
			m->TotMix += mix;
		}
		if(m->Ploidy * temp > popsizes[i]) {
			fprintf(stderr,"Error! Calling for Ploidy * (%d + %d) = %d chromosomes from population %d, but ms will only produce %d.  Exiting...\n",base,mix,m->Ploidy*temp,i+1,popsizes[i]);
			exit(1);
		}
		
		if(i>0) {
			m->BaseStarts[i] = m->BaseStarts[i-1] + popsizes[i-1];
		}
		m->MixStarts[i] = m->BaseStarts[i] + m->NumBase[i] * m->Ploidy;
	}
	
		
	printf("BaseStarts : ");
	for(i=0;i<npops;i++) printf("%d ",m->BaseStarts[i]);
	printf("\n");
	
	printf("MixStarts : ");
	for(i=0;i<npops;i++) printf("%d ",m->MixStarts[i]);
	printf("\n");
	
	/* now, we also want to make sure we can handle the SNP ascertainment needs */
	if(m->MarkerType==SNP) {
		for(i=0;i<m->NumAscPops;i++) {
			if(m->AscSets[i][0]>=npops) {
				fprintf(stderr,"Error! Population number %d has been designated a snp ascertainment population, but there are only %d populations in the ms simulation. Exiting...\n",
					m->AscSets[i][0]+1,npops);
				exit(1);
			}
			if(m->AscSets[i][1] > m->NumBase[i]) {
				fprintf(stderr,"Error! You requested that population number %d have %d individuals for snp ascertainment population, but there are only %d in the baseline from that population. Exiting...\n",
					m->AscSets[i][0]+1,m->AscSets[i][1],m->NumBase[i]);
				exit(1);
			}
			if(m->AscSets[i][2] > m->AscSets[i][1] * m->Ploidy / 2.0 ) {
				fprintf(stderr,"Error! Snp ascertainment in population %d calls for ascertainment if ",m->AscSets[i][0]+1);
				fprintf(stderr,"there are more than %d minor alleles in %d gene copies, but that is more than half of all the gene copies! Exiting...\n",
					m->AscSets[i][2], m->AscSets[i][1] * m->Ploidy);
				exit(1);
			}
			
		}
	}
}

/* this adds Locus Loc to the Basegenos and Mixgenos fields.  The mutation/ascertainment model used
is determined by the options in m.  Loc should be 0 up to NumLoc-1.  This function returns a 1 if the locus
was polymorphic (i.e. was an ascertained SNP or had more microsatellite alleles than 1.  It returns a -1 if the locus is monomorphic. */
int AddLocus(msgeno_opts *m, char **list, int segsites, int Loc) 
{
	int i,j;
	
	/* do some memory allocation / error catching */
	if(m->TotBase>0 && m->Basegenos==NULL) {
		if(Loc==0) {
			m->Basegenos = (int ***)calloc(m->NumPops,sizeof(int **));
			for(i=0;i<m->NumPops;i++)  {
				m->Basegenos[i] = (int **)calloc(m->NumBase[i],sizeof(int *));
				for(j=0;j<m->NumBase[i];j++)  {
					m->Basegenos[i][j] = (int *)calloc(m->NumLoc * m->Ploidy, sizeof(int));
				}
			}
		}
		else {
			fprintf(stderr,"Error! We are on locus %d and m->TotBase = %d but m->Basegenos is still NULL. Exiting...\n",Loc,m->TotBase);
		}
	}
	if(m->TotMix>0 && m->Mixgenos==NULL) {
		if(Loc==0) {
			m->Mixgenos = (int ***)calloc(m->NumPops,sizeof(int **));
			for(i=0;i<m->NumPops;i++)  {
				m->Mixgenos[i] = (int **)calloc(m->NumMix[i],sizeof(int *));
				for(j=0;j<m->NumMix[i];j++)  {
					m->Mixgenos[i][j] = (int *)calloc(m->NumLoc * m->Ploidy, sizeof(int));
				}
			}
		}
		else {
			fprintf(stderr,"Error! We are on locus %d and m->TotMix = %d but m->Mixgenos is still NULL. Exiting...\n",Loc,m->TotMix);
		}
	}


	
	return(SimLoc(m,list,segsites,Loc));
	
}



int ListToInt(char c)
{
	switch(c) {
		case '0':
			return(0);
		case '1':
			return(1);
		default:
			fprintf(stderr,"Error! Unrecognized character %c in list.  Exiting...\n",c);
			exit(1);
	}
}


/* using the ascertainment info stored in m, return the index of the
mutation between 0 and segsites-1 that is the first polymorphic site
which satisfies the ascertainment criteria.  If a SNP is not ascertained,
then return a -1 */
int FindTheSNP(int segsites, char **list, msgeno_opts *m)
{
	int a,i,j,k,Start;
	int Cnts[2],totCnts[2], temp;
	int Genos[2][2], totGenos[2][2];
	long *Perm;
	int r,d;  /* hold the values of mutations as to whether they produce bands (d) or not (r) for aflp simulations stuff */
	
	
	
	r=UniformRV(0,1);  /* if r==0 then the 0 is the recessive, and if r==1 then the 1 is recessive.  d is the opposite*/
	d = (r+1)%2;
	
	
	/* each time we call this function, we want to cycle over the 
	   ascertainment populations in a different order.  So we permute
	   the array of their subscripts */
	Perm = (long *)calloc(m->NumAscPops,sizeof(long));
	for(a=0;a<m->NumAscPops;a++)  {
		Perm[a] = a;
	}
	if( !(m->NoPerm) ) {
		genprm(Perm,m->NumAscPops);
	}
	
	
	
	for(i=0;i<segsites;i++)  {
		totCnts[0] = totCnts[1] = 0;
		totGenos[0][0]=totGenos[0][1]=totGenos[1][0]=totGenos[1][1]=0;
		for(a=0;a<m->NumAscPops;a++)  {
			j = Perm[a];  /* j is the index of the population being used for ascertainment */
			Cnts[0] = Cnts[1] = 0;
			Start = m->BaseStarts[m->AscSets[j][0]];  /* the index of the first individual in this population */
			
			
			
			/* if we are just ascertaining on number of alleles, just cycle over chromosomes */
			if(m->AscSets[j][2]>=-1) {
				for(k=Start;k<Start + (m->Ploidy * m->AscSets[j][1]); k++)  {  /* cycle over indivs in the ascertainment set */
					Cnts[ListToInt(list[k][i])]++;
				}
			}
			
			
			/* on the other hand, if we are ascertaining on the occurrence of all three genotypes we have to do something different */
			else if(m->AscSets[j][2]==-2 || m->AscSets[j][2]==-3) {
				if(m->Ploidy != 2) {
					fprintf(stderr,"Error! You've requested to ascertain on the basis of occurrence of all three genotypes, but ploidy is not two.  Exiting...\n");
					exit(1);
				}
				Genos[0][0]=Genos[0][1]=Genos[1][0]=Genos[1][1]=0;
				for(k=Start;k<Start + (2 * m->AscSets[j][1]); k+=2)  {  /* cycle over indivs in the ascertainment set */
					Genos[ ListToInt(list[k][i]) ]  [ ListToInt(list[k+1][i]) ]++;
					Cnts[ListToInt(list[k][i])]++;  /* we still count these in case we want to ascertain on occurrence of _alleles_ in all pops */
					Cnts[ListToInt(list[k+1][i])]++;
				}
				
				/* count up all the total genotypes */
				totGenos[0][0] += Genos[0][0];
				totGenos[0][1] += Genos[0][1];
				totGenos[1][0] += Genos[1][0];
				totGenos[1][1] += Genos[1][1];

			}
			
			
			/* add these quantities to the overall counts */
			totCnts[0] += Cnts[0];
			totCnts[1] += Cnts[1];
			
			/* now assess whether this would be ascertained in the current population j. 
			   recall that if AscSets[j][2]==-1 or -3, then we aren't doing population-specific ascertainment with this pop, and if AscSets[j][2]==-2 we
			   are ascertaining by genotypes. */
			temp = ECA_MIN(Cnts[0],Cnts[1]);
			if(m->AscSets[j][2] > -1 &&  temp >= m->AscSets[j][2]) { /* if yes, return the index of the site */
				free(Perm);
				return(i);
			}
			else if(m->AscSets[j][2]==-2) {
				printf("We Have A Pop (%d) Counting Genotypes! Site %d.  [0][0]= %d  Het= %d  [1][1]= %d\n",j,i,Genos[0][0],Genos[0][1]+Genos[1][0],Genos[1][1]);
				if(Genos[0][0]>0   &&   Genos[0][1]+Genos[1][0]>0   &&  Genos[1][1]>0) {  /* if all three genotypes appear */
					free(Perm);
					return(i);
				}
			}
		}
		/* if you have gotten down to here, then there were no sites that satisfied the population-specific
		ascertainment criteria.  Now, we do the all-pops-asc bit. */ 
		temp = ECA_MIN(totCnts[0],totCnts[1]);
		if(m->AllPopsAsc > -1 && temp >= m->AllPopsAsc) {
			free(Perm);
			return(i);
		}
		
		/* and, if we got here, we are looking for three genotypes from amongst any of the populations that we looked at */
		if(m->AllPopsGenoAsc==1) {
			printf("Counting Total Genotypes:  00= %d  Het= %d  11= %d\n",totGenos[0][0],totGenos[0][1]+totGenos[1][0],totGenos[1][1]);
			if(totGenos[0][0]>0   &&   totGenos[0][1]+totGenos[1][0]>0   &&  totGenos[1][1]>0) {  /* if all three genotypes appear */
				free(Perm);
				return(i);
			}
		}
		else if(m->AllPopsPseudoAFLP==1) { 
			/*printf("Recessive=%d  Dom=%d.  Counting Total Genotypes:  00= %d  Het= %d  11= %d\n",r,d,totGenos[0][0],totGenos[0][1]+totGenos[1][0],totGenos[1][1]); */
			if(totGenos[r][r]>=m->MinAFLP_NoBand   &&   totGenos[r][d]+totGenos[d][r] + totGenos[d][d]>=m->MinAFLP_Band) {  /* if all three genotypes appear */
				free(Perm);
				printf("AFLP_CONVERSION: Recessive= %d\n",r+1);
				return(i);
			}
		}
	}
	free(Perm);
	return(-1);
}

int SimLoc(msgeno_opts *m, char **list, int segsites, int Loc) 
{
	int i,j,k,sign,step,ind,gene,cnt;
	int *Steps;
	int *Hash;
	int min=999999999, max=-9999999;
	int newmin, newmax;
	int theSNP;
	
	/* We have some preliminary stuff do to if this is an MSAT locus */
	if(m->MarkerType==MSAT) {
		/* first, determine if each mutation is a step up or down, and how far */
		Steps = (int *)calloc(segsites, sizeof(int));
		Hash = (int *)calloc(10001,sizeof(int));
		for(i=0;i<segsites;i++)  {
			if(ranf()<.5) {
				sign = -1;
			}
			else {
				sign = 1;
			}
			if(ranf()<m->MultiMutRate) {
				step = 2 + ignpoi((float)m->MultiStepSize);
			}
			else {
				step = 1;
			}
			Steps[i] = sign * step;
		}
	}
	
	
	/* if we are doing a SNP, then we need to determine which of the segregating sites is our SNP.
	   We just go through them until we find the first that satisfies the ascertainment criteria */
	if(m->MarkerType==SNP) {
		theSNP =  FindTheSNP(segsites,list,m);
printf("THESNP : %d\n",theSNP);
		if(theSNP==-1) {
			return(-1);  /* if we didn't find a SNP, just return -1 */
		}
	}
	
	
	/* then, cycle over pops, indivs, and gene copies and fill 'em up.  We do different things with them
	   depending on whether we are doing SNPs or microsatellites.  */
	for(i=0;i<m->NumPops;i++)  {
	
		/* do the Base guys */
		ind = gene = cnt = 0;
		for(j=m->BaseStarts[i];j<m->BaseStarts[i]+m->NumBase[i]*m->Ploidy;j++,cnt++)  {
			ind = cnt / m->Ploidy;
			gene = Loc * m->Ploidy + (cnt % m->Ploidy);
			m->Basegenos[i][ind][gene] = 0;
			
			
			if(m->MarkerType==MSAT)  {  /* if MSAT, then each mutation is interpreted as a step up or down */
				for(k=0;k<segsites;k++) {
					m->Basegenos[i][ind][gene] += ListToInt(list[j][k]) * Steps[k];
				}
			}
			if(m->MarkerType==SNP) {
				m->Basegenos[i][ind][gene] = ListToInt(list[j][theSNP]) + 1; /* we code snps as 1's and 2's */
			}
			
			
			/* here we keep track of the min and max allele length */
			if(m->Basegenos[i][ind][gene] < min) {
				min = m->Basegenos[i][ind][gene];
			}
			if(m->Basegenos[i][ind][gene] > max) {
				max = m->Basegenos[i][ind][gene];
			}
		} 
		
		/* do the Mix guys */
		ind = gene = cnt = 0;
		for(j=m->MixStarts[i];j<m->MixStarts[i]+m->NumMix[i]*m->Ploidy;j++,cnt++)  {
			ind = cnt / m->Ploidy;
			gene =  Loc * m->Ploidy + (cnt % m->Ploidy);
			m->Mixgenos[i][ind][gene] = 0;
			
			
			
			if(m->MarkerType==MSAT)  {  /* if MSAT, then each mutation is interpreted as a step up or down */
				for(k=0;k<segsites;k++) {
					m->Mixgenos[i][ind][gene] += ListToInt(list[j][k]) * Steps[k];
				}
			}
			if(m->MarkerType==SNP) {
				m->Mixgenos[i][ind][gene] = ListToInt(list[j][theSNP]) + 1; /* we code snps as 1's and 2's */
			}
			
			
			
			if(m->Mixgenos[i][ind][gene] < min) {
				min = m->Mixgenos[i][ind][gene];
			}
			if(m->Mixgenos[i][ind][gene] > max) {
				max = m->Mixgenos[i][ind][gene];
			}

		} 
	}
	
	
	/* and then finally, if we are  doing MSAT's we will add  enough to all of them so that
	 the min length is 100.  Maybe eventually check to see if any lengths are > 999 */
	if(m->MarkerType==MSAT) {
		newmin = 100;
		newmax = max+100-min;
		m->NumAlle[Loc]=0;
		for(i=0;i<m->NumPops;i++)  {
			ind = gene = cnt = 0;
			for(j=m->BaseStarts[i];j<m->BaseStarts[i]+m->NumBase[i]*m->Ploidy;j++,cnt++)  {
				ind = cnt / m->Ploidy;
				gene = Loc * m->Ploidy + (cnt % m->Ploidy);
				m->Basegenos[i][ind][gene] += 100 - min;
				Hash[m->Basegenos[i][ind][gene]]++;
			} 
			ind = gene = cnt = 0;
			for(j=m->MixStarts[i];j<m->MixStarts[i]+m->NumMix[i]*m->Ploidy;j++,cnt++)  {
				ind = cnt / m->Ploidy;
				gene = Loc * m->Ploidy + (cnt % m->Ploidy);
				m->Mixgenos[i][ind][gene] += 100 - min;
				Hash[m->Mixgenos[i][ind][gene]]++;
			}
		}
	}
	/* now, count alleles */
	if(m->MarkerType==SNP) {
		m->NumAlle[Loc]=2;
	}
	else if(m->MarkerType==MSAT)  {
		for(j=newmin;j<=newmax;j++)  {
			if(Hash[j]>0) {
				m->NumAlle[Loc]++;
			}
		}
		free(Steps);
		free(Hash);
	}
	
	if(m->NumAlle[Loc]>1) {
		return(1);
	}
	else {
		return(-1);
	}
	
}


void PrintDataSet(msgeno_opts *m, int DataSetNum) 
{
	int i,j,k,a;
	FILE *out=stdout;
	char str[1000];
	int **LocHash;
	
	
	/* print the baseline file */
	sprintf(str,"BaseFile_%d.txt",DataSetNum);
	out = fopen(str,"w");
	
	fprintf(out,"%d  %d\n",m->TotBase,m->NumLoc);
	for(i=0;i<m->NumLoc;i++)  {
		fprintf(out,"Locus_%d\n",i+1);
	}
	for(i=0;i<m->NumPops;i++)  {
	
		fprintf(out,"POP PopNum_%d_\n",i+1);
		for(j=0;j<m->NumBase[i];j++) {
			fprintf(out,"Pop_%d_BaseInd_%d",i+1,j+1);
			for(k=0;k<m->NumLoc;k++)  {
				fprintf(out,"   ");
				for(a=0;a<m->Ploidy;a++)  {
					fprintf(out," %d",m->Basegenos[i][j][k*m->Ploidy+a]);
				}
			}
			fprintf(out,"\n");
		}
	}
	fclose(out);
	
	
	/* print the mix file */
	if(m->TotMix>0) {
		sprintf(str,"MixFile_%d.txt",DataSetNum);
		out = fopen(str,"w");

		fprintf(out,"%d  %d\n",m->TotMix,m->NumLoc);
		for(i=0;i<m->NumLoc;i++)  {
			fprintf(out,"Locus_%d\n",i+1);
		}
		fprintf(out,"POP MixtureSample\n");
		for(i=0;i<m->NumPops;i++)  {
			for(j=0;j<m->NumMix[i];j++) {
				fprintf(out,"Pop_%d_MixInd_%d",i+1,j+1);
				for(k=0;k<m->NumLoc;k++)  {
					fprintf(out,"   ");
					for(a=0;a<m->Ploidy;a++)  {
						fprintf(out," %d",m->Mixgenos[i][j][k*m->Ploidy+a]);
					}
				}
				fprintf(out,"\n");
			}
		}	
		fclose(out);
	}
	
}

msgeno_opts *GetMSGenoOpts(int argc, char *argv[])
{
	msgeno_opts *m = (msgeno_opts *)malloc(sizeof(msgeno_opts));
	int PloidyF=0,
		lociF=0,
		baselineF=0,
		mixtureF=0,
		msatF = 0,
		snpF =0,
		allpopsascF=0,
		noperm_f=0,
		overrideF=0,
		ms_fileF=0,
		allpopsgenoF=0,
		allpopsaflpF=0;
	
	
	
	DECLARE_ECA_OPT_VARS;

	SET_OPT_WIDTH(25);
	SET_ARG_WIDTH(14);
	
	SET_PROGRAM_NAME("ms2geno");
	SET_PROGRAM_SHORT_DESCRIPTION("convert ms output to genotypes in files");
	SET_PROGRAM_LONG_DESCRIPTION(
		This program is designed to convert output from the makesamples (also called ms) program\054 written by Richard Hudson\054
		into data sets of microsatellite or SNP genotypes.  ms yields sequences of 0s and 1s along chromosomes from the infinite sites model---a
		1 means the sampled chromosome inherited a certain mutated nucleotide and a 0 means that received the non-mutated form of the nucleotide.
		It is straightforward to use such strings of 0s and 1s to simulate microsatellites via a stepwise mutation process\054 or as
		the background for ascertaining SNPs\054 especially when each locus is assumed  unlinked and simulated on a separate\054 independent\054 coalescent
		tree. 
		
		\n\nThis version of ms2geno assumes that each locus is independently segregating.  A future version could be made which used the recombination
		simulation features of ms to allow for linked markers.  However\054 it does not do that yet!
				
		\n\nThe code for reading ms output into simple data structures was taken from one of the little ms output analysis programs that Dick Hudson
		distributes with makesamples.  Thanks for providing source code Dick!
		
		\n\nHere is how the program works.  First you have to simulate from ms.  You have to make sure that you have twice as many chromosomes (set with
		the nsam option in ms) as you
		want to have diploid genotypes in your ms2geno output data sets.  If you run ms with samples from multiple subpopulations using the -I option\054 
		ms2geno will be able to detect that and partition the individuals appropriately into the subpopulations.   The howmany option in ms sets how many 
		independent coalescent trees are simulated by ms.  Each one corresponds to a separate locus.  If you set howmany=100 in ms\054 say\054 but choose
		to simulate only 5 loci using the -l/--num-loci option with ms2geno\054 then ms2geno will output 100/5=20 separate data sets\054 each with five
		loci.  
		
		
		\n\nThe output goes into files: BaseFile_XX.txt and MixFile_XX.txt where XX is the data set number.  The format is a sort of quasi-GenePop format
		that works with Steven Kalinowski\047s program GMA\054 and with my program gsi_sim.  It might even work with GenePop.  BaseFile and MixFile are just two separate files that 
		
		
		\n\nYou can control the genetic variability by using the -t (theta) option to ms.  With microsatellites\054 increasing the multi-step mutation rate
		(see option -u/---msat) in ms2geno will also increase the variability.  If you are simulating SNPs with ms2geno it is a good idea to simulate
		more than enough trees in ms for all the loci because SNPs might not get ascertained from some of the simulated trees.  You can use the -s
		option to ms to specify how many segregating sites should be available for SNP discovery.  The details of SNP discovery will be covered elsewhere.
		
		\n\nIf the file ms2geno_seeds is  present in the current working directory then its contents are used to seed the random
		number generator.  If not\054 then seeds are generated from the current time.  The current state of the random number generator
		is written to ms2geno_seeds at the end of execution so it may be used over and over again in a sane fashion with 
		respect to the random number generator. 
		
		\n\nThis has been an incomplete description\054 but I just wanted to get something in there for guiLiner.

		
	)
	SET_VERSION("VERSION: 1.1 Beta.  1 August 2007")
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson (eric.anderson@noaa.gov)");
	SET_VERSION_HISTORY(" Program conceived around November 2006. ")


	/* set some defaults */
	m->Ploidy = 2;
	m->NumPops=0;
	m->Basegenos = NULL;
	m->Mixgenos = NULL;
	m->IDs = NULL;
	m->MinGeno = NULL;
	m->NumBase = NULL;
	m->NumMix = NULL;
	m->NumAlle = NULL;
	m->MultiMutRate = 0;
	m->MultiStepSize = 0;		
	m->NumAscPops = 0;
	m->AllPopsAsc = -1;
	m->AscSets = (int **)calloc(MAX_NUM_POPS,sizeof(int **));
	m->NoPerm = 0;
	m->DoOverride=0;
	m->InputFileName[0]='\0';
	m->AllPopsGenoAsc=0;
	m->AllPopsPseudoAFLP=0;
			
	BEGIN_OPT_LOOP 	 
		
		OPEN_SUBSET(Command Line Switches To Input Options to ms2geno,Simulation Options,These options determine how ms2geno uses output from ms to simulated data sets);
		
		if(OPTION(
			Ploidy,
			PloidyF,
			,
			ploidy,
			J,
			J is the ploidy of the organisms to be simulated,
			Each individual is imagined to carry J copies of each gene.  This option is not required.
			The default is 2.
			) ){
			if(ARGS_EQ(1)){
				m->Ploidy = GET_INT;
			}
		} 
		if(REQUIRED_OPTION(
			Number of Loci,
			lociF,
			l,
			num-loci,
			J,
			Number of loci to be included in the data sets,
			This will be the number of loci included in each data set)) {
			if(ARGS_EQ(1)) {
				m->NumLoc = GET_INT;
				m->NumAlle = (int *)calloc(m->NumLoc,sizeof(int));
			}
		}
		if(REQUIRED_OPTION(
			Baseline Sample Sizes,
			baselineF,
			b,
			baselines,
			J1 ... JK,
			Include J1...JK indivs from the K populations,
			Sets the number of individuals to be sampled from the K simulated populations into
			the baseline samples.  K must be equal to the number of populations given using the -I 
			argument to ms.  If this option and -m/--mixtures are used in the same command line they
			must have the same number of arguments K. )) {
			if(ARGS_GEQ(1)) { int temp; int i;
				temp = COUNT_ARGS;
				if(mixtureF) {
					if(temp != m->NumPops) {
						fprintf(stderr,"Error! Number of arguments for options -m/--mixtures and -b/--baselines don't match. Exiting...\n");
						OPT_ERROR;
					}
				}
				m->NumPops = temp;
				m->NumBase = (int *)calloc(temp,sizeof(int));
				for(i=0;i<temp;i++)  {
					m->NumBase[i] = GET_INT;
				}
			}
		}
		if(REQUIRED_OPTION(
			Mixture Sample Sizes,
			mixtureF,
			m,
			mixtures,
			J1 ... JK,
			Include J1...JK indivs from the K populations in the mixture,
			Sets the number of individuals to be sampled from the K simulated populations into
			the  mixture samples.  K must be equal to the number of populations given using the -I 
			argument to ms.  If this option and -b/--baselines are used in the same command line they
			must have the same number of arguments K. )) {
			if(ARGS_GEQ(1)) { int temp; int i;
				temp = COUNT_ARGS;
				if(baselineF) {
					if(temp != m->NumPops) {
						fprintf(stderr,"Error! Number of arguments for options -m/--mixtures and -b/--baselines don't match. Exiting...\n");
						OPT_ERROR;
					}
				}
				m->NumPops = temp;
				m->NumMix = (int *)calloc(temp,sizeof(int));
				for(i=0;i<temp;i++)  {
					m->NumMix[i] = GET_INT;
				}
			}
		}	
		if(CLASHABLE_OPTION(
			Microsatellite Evolution Parameters,
			msatF,
			u,
			msat,
			R1 R2,
			Parameters for stepwise mutation model,
			R1 is the probability of multistep mutations and R2 is the mean of the Poisson distribution that 
			determines the step size of multistep mutation.  The length of multistep mutations is 2 + a draw from
			a Poisson distribution with mean R2.  Applying this option invokes the microsatellite/stepwise mutation model.
			This is the default.  The default value for R1 is 0 in which case R2 is irrelevant.  This options cannot be used in conjunction with
			the -s/--snp option.,
			snpF,
			is issued and the -s/--snp option is also invoked )) {  
				m->MarkerType = MSAT;
				m->MultiMutRate = GET_DUB;
				m->MultiStepSize = GET_DUB;
		}
		if(MULT_USE_OPTION(
			SNP Ascertainment Criteria,
			snpF,
			s,
			snp-asc,
			J1 J2 J3,
			sets the criteria for SNP ascertainment relative to a population J1,
			This says that a SNP will be ascertained if it there are J3 or more copies of the minor allele among J2 individuals from population J1.
			This option may be given any number of times; however it should not be given multiple times for any single population.
			You may wish to count the number of 0s and 1s in a number of different populations and then ascertain on the basis
			of the number of copies of the minor allele across all those populations.  This can be achieved by calling this option
			for each of the populations that you want to have involved; but for each one use J3=-1.  Then give the overall minimum number
			of minor alleles needed using the --all-pops-asc option.  We also have a quick hack\054 hijacking this option\054 that allows
			us to do ascertainment on the basis of finding at least one of every possible diploid genotype.  In other words\054 you will only
			ascertain the SNPs if the 00 and 11 homozygotes appear in the ascertainment sample\054 along with the 01 or 10 heterozygotes.  To use this
			option\054 you let the J3 argument be -2.  Finally\054 if you want to do the pseudo-AFLP thing then you specify the population with the J1 argument
			and you set the J3 argument to -3.  In this case the ascertainment is not done population by population.  With the pseudo-AFLP thang\054 the 
			ascertainment always depends on the numbers of phenotypes amongst all the populations having J3=-3. ,
			MAX_NUM_POPS)) {
				if(ARGS_EQ(3)) { int temp; int pop;
					temp = GET_INT;
					if(temp <= 0) {
						fprintf(stderr,"Error! Argument J1 to -s/--snp-asc is less than 1, but it must be greater than or equal to 1 because it is an index of a population...\n");
						OPT_ERROR;
						temp=GET_INT;
						temp=GET_INT;
					}
					else {
						pop = temp-1;
						/* allocate memory to store the rest of the info */
						m->AscSets[m->NumAscPops] = (int *)calloc(3,sizeof(int));
						m->AscSets[m->NumAscPops][0] = pop;
						
						temp = GET_INT;
						if(temp<=0) {
							fprintf(stderr,"Error! The J2 argument to -s/--snp-asc option is less than 1.  This doesn't make sense...\n");
							OPT_ERROR;
							temp=GET_INT;
						}
						else {
							m->AscSets[m->NumAscPops][1] = temp;
							
							temp = GET_INT;
							if(temp < -3) {
								fprintf(stderr,"Error! The J3 argument to -s/--snp-asc option is less than -2.  This doesn't make sense...Remember use J3=-1 to do no population-specific ascertainment, and J3=-2 to require all three genotypes in diploids, and J3=-3 to do the pseudo-aflp thing.\n");
								OPT_ERROR;
							}
							else {
								m->AscSets[m->NumAscPops][2] = temp;
								if(temp>=0) {
									printf("AscertainmentSetNumber %d  is %d  Individuals from Population %d,  with ascertainment if minor allele count >= %d\n",
									m->NumAscPops+1,m->AscSets[m->NumAscPops][1],m->AscSets[m->NumAscPops][0]+1,m->AscSets[m->NumAscPops][2]);
								}
								else if(temp==-1) {
									printf("AscertainmentSetNumber %d  is %d  Individuals from Population %d,  with ascertainment depending on the all-pops-asc criterion.\n",
									m->NumAscPops+1,m->AscSets[m->NumAscPops][1],m->AscSets[m->NumAscPops][0]+1);
								}
								else if(temp==-2)  {
									printf("AscertainmentSetNumber %d  is %d  Individuals from Population %d,  with ascertainment occurring if all three genotypes are observed in the ascertainment sample.  This is applicable only to diploids.\n",
									m->NumAscPops+1,m->AscSets[m->NumAscPops][1],m->AscSets[m->NumAscPops][0]+1);
								}
								m->NumAscPops++;
								m->MarkerType = SNP;
							}
						}
					
					}
					
				}
		}
		if(OPTION(
			Pooled SNP Ascertainment,
			allpopsascF,
			,
			all-pops-asc,
			J,
			number of minor alleles amongst all ascertainment pops that trigger ascertainment,
			This option can be used for SNP ascertainment.  All the gene copies in ascertainment sets 
			specified by the -s/--snp-asc option go into a pool and a SNP will be ascertained if the least
			frequent site allele occurs in J or more copies in that pool.  Note that this criterion is applied for each site after
			all of the population-specific criteria given in the -s/--snp-asc options.  To ensure that ascertainment is based 
			only on the total number of genes amongst all the ascertainment populations you should make sure that the J3 argument to 
			all the -s/--snp-asc option calls is -1. )) {
			if(ARGS_EQ(1)) {
				m->AllPopsAsc = GET_INT;
			}
		}
		if(OPTION(
			SNPs via Pooled Genotypes,
			allpopsgenoF,
			,
			all-pops-geno-asc,
			,
			base SNP ascertainment on the occurrence of all three genotypes in multiple populations,
			This option can be used when any of the populations have ascertainment samples that are specified with using 
			the -s or --snp-asc option with the J3 argument equal to -2.  Such an option causes a SNP to be ascertained in such 
			a population if all three genotypes are observed in the population.  If you issue this --all-pops-geno-asc option\054 then
			even if all three genotypes are not seen in any of the ascertainment samples\054 a SNP is still ascertained so long as all
			three genotypes appear amongst the inidividuals in the ascertainment samples that had the J3 argument equal to -2.)) {
				if(ARGS_EQ(0)) {
					m->AllPopsGenoAsc=1;
				}
		 } 
		
		
		if(OPTION(
			SNP ascertainment for pseudo-AFLPs,
			allpopsaflpF,
			,
			all-pops-pseudo-aflp,
			J1 J2,
			base SNP ascertainment on finding certain numbers of pseudo-AFLP SNPs,
			This is a bit of a hack.  I am busy doing some simulations and I need this functionality:  This will produce SNPs and it will
			also produce a list of which SNP allele is the recessive allele in a dominant marker system like ALPs.  With that information\054 
			you can then hybridize individuals using their SNP genotypes and process those results into AFLP phenotypes later.  Using this option\054 
			a SNP will be ascertained if at least J1 recessive (non-band-producing) phenotypes and J2 band-producing phenotypes are seen from amongst
			all the individuals in the populations that get the -s or --snp-asc option having the J3 argument equal to -3.
			)) {
				if(ARGS_EQ(2)) { int temp;
					m->AllPopsPseudoAFLP=1;
					temp = GET_INT;
					if(temp<0) {
						fprintf(stderr,"Error! J1 argument to all-pops-pseudo-aflp < 0.  Weird.  Exiting...\n");
						exit(1);
					}
					m->MinAFLP_NoBand=temp;
					
					temp = GET_INT;
					if(temp<0) {
						fprintf(stderr,"Error! J2 argument to all-pops-pseudo-aflp < 0.  Weird.  Exiting...\n");
						exit(1);
					}
					m->MinAFLP_Band=temp;
				}
		 } 
		
		
		if(OPTION(
			Do Not Permute Ascertainment Pops,
			noperm_f,
			,
			no-perm,
			,
			do not permute ascertainment population labels,
			Invoking this option inhibits the default behavior of permuting the ascertainment population labels for each SNP.  Permuting is the
			default behavior because it seems that you could have some biases if you scan through the ascertainment populations in the same order
			every time.)) {
			if(ARGS_EQ(0)) {
				m->NoPerm = 1;
			}
		}
		if(OPTION(
			Override ms Number of Chromosomes,
			overrideF,
			,
			ms-pop-override,
			J1 ... JK,
			Override the ms command line number of chromosomes,
			This odd little option is used to specify to ms2geno the number of chromosomes in each population simulated by ms WHEN YOU WANT THAT
			NUMBER TO BE DIFFERENT FROM WHAT WAS SPECIFIED ON THE MS COMMAND LINE.  Please note that this should be in terms of the number
			of chromosomes and not in terms of the number of diploid individuals.  This is not an option you normally need to use.  The only time
			it should be invoked is when you have simulated a large number of chromosomes in a SINGLE population with ms but want you ms2geno to carve
			up the individuals from that single population into separate populations.  Obviously the number of populations should match what is given for 
			baseline and mixture files.  ms2geno will check to make sure that enough chromosomes were simulated in ms to cover the ones requested here. ) ) {
			if(ARGS_GT(0)) { int i;
				m->DoOverride=1;
				m->OverRideNumPops = COUNT_ARGS;
				m->OverRideNumsInPops = (int *)calloc(m->OverRideNumPops,sizeof(int));
				for(i=0;i<m->OverRideNumPops;i++)  { 
					m->OverRideNumsInPops[i] = GET_INT;
				}
			}
		}
		
		if(OPTION(
			Input File,
			ms_fileF,
			f,
			ms-file,
			F,
			Take input from a file instead of from standard input,
			Should you prefer to have redirected the ms output into a file and saved it\054 you can use that
			saved file as the input to ms2geno by either redirecting it to standard input (with < file.txt) or
			you can name the file on the command line using this option.  This is included because some operating systems
			may not have good file redirection facilities\054 and because this makes it much easier to use ms2geno with
			the guiLiner front end.
			)) {
			if(ARGS_EQ(1))	{
				GET_STR(m->InputFileName);
			}
		}
		
	CLOSE_SUBSET;
		
	END_OPT_LOOP
	
	return(m);

}
	

/* allocates space for gametes (character strings) */
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}

        int
biggerlist(nsam, nmax, list )
        int nsam ;
        unsigned nmax ;
        char ** list ;
{
        int i;

        maxsites = nmax  ;
        for( i=0; i<nsam; i++){
           list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
           if( list[i] == NULL ) perror( "realloc error. bigger");
           }
}                        


	double
nucdiv( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
   	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(1.0 -p1)*nnm1 ;
		}
	return( pi ) ;
}

/*   thetah - pi   */
	double
hfay( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
   	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ;
		}
	return( -pi ) ;
}

/* Fay's theta_H  */
        double
thetah( int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        double pi, p1, nd, nnm1  ;

        pi = 0.0 ;

        nd = nsam;
        nnm1 = nd/(nd-1.0) ;
        for( s = 0; s <segsites; s++){
                p1 = frequency('1', s,nsam,list) ;
                pi += p1*p1 ; 
                }
        return( pi*2.0/( nd*(nd-1.0) )  ) ;
}


        int
frequency( char allele,int site,int nsam,  char **list)
{
        int i, count=0;
        for( i=0; i<nsam; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
        return( count);
}        

	int
segsub( int nsub, int segsites, char **list )
{
	int i, count = 0 , c1 ;
	int frequency( char, int, int, char**) ;

	for(i=0; i < segsites ; i++){
	  c1 = frequency('1',i,nsub, list);
	  if( ( c1 > 0 ) && ( c1 <nsub )  ) count++;
	  }
	return( count ) ;
}
	
