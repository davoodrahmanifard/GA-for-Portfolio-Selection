#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "conio.h"
#include "stdio.h"
#include "dos.h"
#include "time.h"
#include "process.h"
#include "io.h"


#define maxpop 10000
#define maxstring 300


typedef enum _Boolean {
  False=0,
  True
} Boolean;


/*typedef Boolean allele;     /*Allele = bit position */
typedef Boolean *chromosome;


struct individual1
{
Boolean chrom[maxstring];	/* Genotype = bit string */
float fitness;				/* Objective function value */
Boolean tag;
};


struct individual2
{
Boolean chrom[maxstring];	/* Genotype = bit string */
int code[maxstring];		/* Code of asset */
int percent[maxstring];		/* Phenotype */
float fitness;				/* Objective function value */
float sumprct;				/* Sum of phenotypes */
Boolean tag;				/* Validity of chromosome */
};


typedef struct individual1 *population1;
typedef struct individual2 *population2;


/*   variables   */

population1 oldpop1;
population1 newpop1;
population1 poptemp1;
population1 mpop1,popp11;

population2 oldpop2;
population2 newpop2;
population2 poptemp2;
population2 popp12,mpop2;

int popsize1, lchrom1, gen, maxgen1, nmutation, ncross, maxselect,popsize2, maxgen2;
float pcross1, pmutation1, sumfitness, avg, max, min, cmax=0.0,pcross2, pmutation2;
float r[maxstring],s[maxstring],p[maxstring][maxstring],c[maxstring][maxstring],l[maxstring][maxstring],t[maxstring];
float oldrand[55]={0};
int jrand;

time_t start,finish;
FILE *fpt;
FILE *fpr;
FILE *fdata1;
FILE *fdata2;
FILE *fdata3;
FILE *fcorr,*fdata,*fcorr1;


void decode(population2 pop)
{
	int i,j,k;
	int accum,sum;
	int powerof2;
	
	k=0;
	sum=0;
	for(i=0;i<maxselect;i++) {
		accum=0;
		powerof2=1;
		j=0;
		do {
			if(pop->chrom[k]==True)  accum=accum+powerof2;
			powerof2=powerof2*2;
			j++;
			k++;
		}
		while (j<7);
		pop->percent[i]=accum;
		sum=sum+accum;
	}
	pop->sumprct=sum;
}


float objfunc1(population1 pop)
{

	int i,j;
	float y=0.0,t;
	char ch;

	for(i=0;i<lchrom1;i++)
		for(j=0;j<lchrom1;j++) {
			if (pop->chrom[j]==True && i!=j && pop->chrom[i]==True && r[i]!=0.0) {
				if (s[i]==cmax) t=0;
				else t=log10(1-s[i]/cmax);
				if (i>j) 
					y=y+(l[j][i])+log10(r[i])+t;
				else
					y=y+(l[i][j])+log10(r[i])+t; 
			}
		}
	return (y);
}


double objfunc2(population2 pop)
{

	int i,j;
	double y=0,z=0;
	double x=10;
	char ch;

	for(i=0;i<maxselect;i++) {
		y=y+(pop->percent[i])*(pop->percent[i])*(s[pop->code[i]]*s[pop->code[i]])/(pop->sumprct*pop->sumprct);
		if (y<0){
	printf("\n %.6f %8d %.6f %.6f ",y,pop->percent[i],s[pop->code[i]],pop->sumprct);
		ch=getchar();}
	}
	
	for(i=0;i<maxselect;i++) {
		for(j=0;j<maxselect;j++) {
			if (i!=j){ 
				if (i>j)  y=y+p[pop->code[j]][pop->code[i]]*(s[pop->code[i]]*s[pop->code[j]])*(pop->percent[j])*(pop->percent[i])/(pop->sumprct*pop->sumprct);
		      else  z=z+p[pop->code[i]][pop->code[j]]*(s[pop->code[i]]*s[pop->code[j]])*(pop->percent[j])*(pop->percent[i])/(pop->sumprct*pop->sumprct);

  }
				}}

	return (y);
}



statistics1(population1 pop)
{
	int j;
	population1 npop;
	sumfitness=0.0;
	min = pop->fitness;
	max = pop->fitness;
    npop=pop;

	for(j=1;j<=popsize1;j++)
	{
		sumfitness=sumfitness+pop->fitness;
		if(pop->fitness>max){max=pop->fitness;}
		if(pop->fitness<min){min=pop->fitness;}
		pop++;
	}
	
	avg=sumfitness/popsize1;
	pop=npop;

    return 0;
}


statistics2(population2 pop)
{
	int j;
	population2 npop;
	sumfitness=0.0;
	min = pop->fitness;
	max = pop->fitness;
    npop=pop;

	for(j=1;j<=popsize2;j++)
	{
		sumfitness=sumfitness+pop->fitness;
		if(pop->fitness>max){max=pop->fitness;}
		if(pop->fitness<min){min=pop->fitness;}
		pop++;
	}
	
	avg=sumfitness/popsize2;
	pop=npop;

    return 0;
}


void advance_random(void)
{
	int j1;
	float new_random;
	for (j1=0;j1<24;j1++) {
		new_random=oldrand[j1] - oldrand[j1+31];
		if (new_random<0.0) new_random=new_random +1.0;
		oldrand[j1]=new_random;
	}
	for (j1=24;j1<55;j1++) {
		new_random=oldrand[j1] - oldrand[j1-24];
		if (new_random<0.0) new_random=new_random+1.0;
		oldrand[j1]=new_random;
	}
}


void warmup_random(float random_seed)
{
	int j1,ii;
	float new_random, prev_random;

	oldrand[55]=random_seed;
	new_random=1.0e-9;
	prev_random=random_seed;
	for (j1=0;j1<54;j1++) {
		ii=(21*j1) % 55;
		oldrand[ii]=new_random;
		new_random=prev_random - new_random;
		if (new_random<0.0) new_random=new_random+1.0;
		prev_random=oldrand[ii];
	}
	
	advance_random();
	advance_random();
	advance_random();	
	jrand=0;

}


float random()
{
	jrand++;
	if (jrand>55) {
		jrand=1;
		advance_random();
	}
	return oldrand[jrand];
}

int rnd(int low,int high)
{
	int i;

	if (low>=high) i=low;
	else {
		i=floor(random()*(high-low+1)+low);
		if (i>high) i=high;
	}
	return(i);
}


int select1(population1 pop)
{
	float k,rand,partsum=0.0;
	int j;
	population1 npop;


	npop=pop;
	k=random();
	rand=k*sumfitness;
	j=0;
	do
	{
		partsum=partsum+(pop->fitness);
		pop++;
		j++;
	}
	while ((partsum<rand) && (j!=popsize1));
	pop=npop;
	return(j-1);
}


int select2(population2 pop)
{
	float k,rand,partsum=0.0;
	int j;
	population2 npop;


	npop=pop;
	k=random();
	rand=k*sumfitness;
	j=0;
	do
	{
		partsum=partsum+(pop->fitness);
		pop++;
		j++;
	}
	while ((partsum<rand) && (j!=popsize2));
	pop=npop;
	return(j-1);
}




Boolean FlipCoin(float p)
{

	Boolean flip;

	if (p==1.0) flip=True;
	else
		flip= (random()<=p ? True : False);
	
	return flip;
}


void writechrom(FILE *fp,chromosome chrom,int lch)
{
	int j;
	for (j=0;j<lch;j++)
		if (chrom[j]==True) fprintf(fp,"1");
		else fprintf(fp,"0");
}


void initdata(void)
{

	int i,j;
	char re[10],st[10],cov[8],cor[5],corl[5];
	
/* Interactive data inquiry and setup */


	printf("****Genetic Algorithm Data Entry and Initialization *****\n");
	printf("\n");
	printf("Prepared by DAVOOD RAHMANI FARD.\n");
	printf("Enter chromosome length (number of asset in market)(lchrom1)---->"); scanf("%d",&lchrom1);
	printf("Enter number of asset to be chosen(maxselect)---->"); scanf("%d",&maxselect);
    printf("\n");
	do
	{
	printf("Enter population size for first algorithm (popsize1)---->"); scanf("%d",&popsize1);
	}
	while (popsize1%2!=0);
	printf("Enter max. generations for first algorithm (maxgen1)---->"); scanf("%d",&maxgen1);
	printf("Enter crossover probability for first algorithm (pcr1)-->"); scanf("%f",&pcross1);
	printf("Enter mutation probability for first algorithm (pmu1)-->"); scanf("%f",&pmutation1);
    printf("\n");
	do
	{
	printf("Enter population size for second algorithm (popsize2)---->"); scanf("%d",&popsize2);
	}
	while (popsize2%2!=0);
	printf("Enter max. generations for second algorithm (maxgen2)---->"); scanf("%d",&maxgen2);
	printf("Enter crossover probability for second algorithm (pcr2)-->"); scanf("%f",&pcross2);
	printf("Enter mutation probability for second algorithm (pmu2)-->"); scanf("%f",&pmutation2);


	oldpop1=(individual1 *)calloc(popsize1,sizeof(individual1));
    if (!oldpop1) {
		printf("allocation failure");
		exit(1);
	}

	newpop1=(individual1 *)calloc(popsize1,sizeof(individual1));
    if (!newpop1) {
		printf("allocation failure");
		exit(1);
	}

	poptemp1=(individual1 *)calloc(popsize1*3,sizeof(individual1));
    if (!poptemp1) {
		printf("allocation failure");
		exit(1);
	}

	popp11=(individual1 *)calloc(1,sizeof(individual1));

	mpop1=(individual1 *)calloc(1,sizeof(individual1));
	
	oldpop2=(individual2 *)calloc(popsize2,sizeof(individual2));
    if (!oldpop2) {
		printf("allocation failure");
		exit(1);
	}

	newpop2=(individual2 *)calloc(popsize2,sizeof(individual2));
    if (!newpop2) {
		printf("allocation failure");
		exit(1);
	}

	poptemp2=(individual2 *)calloc(popsize2*3,sizeof(individual2));
    if (!poptemp2) {
		printf("allocation failure");
		exit(1);
	}

	popp12=(individual2 *)calloc(1,sizeof(individual2));

	mpop2=(individual2 *)calloc(1,sizeof(individual2));

	cmax=0.0;
	for (i=0;i<lchrom1;i++) {
		fscanf(fdata1,"%6s %s",re,st);
	    r[i]=atof(re);
		s[i]=atof(st);
		if(s[i]>cmax){cmax=s[i];}
	}

	for (i=0;i<lchrom1;i++)
		for (j=i+1;j<lchrom1;j++) {
			fscanf(fdata2,"%s %s ",cov,cor);
			c[i][j]=atof(cov);
			p[i][j]=atof(cor);}
for (i=0;i<lchrom1;i++)
		for (j=i+1;j<lchrom1;j++) {
			fscanf(fdata3,"%s", corl);
			l[i][j]=atof(corl);
		}

}


void initreport1(void)
{


	fprintf(fpt,"\n     First Genetic Algorithm Parameters      ");
	fprintf(fpt,"\n     ----------------------------------      ");
	fprintf(fpt,"\n");
	fprintf(fpt,"\nPopulation size (popsize1) = %d",popsize1);
	fprintf(fpt,"\nChromosome length (lchrom1) = %d",lchrom1);
	fprintf(fpt,"\nNumber of asset to be chosen (maxselect) = %d",maxselect);
	fprintf(fpt,"\nMax. generations for first algorithm(maxgen1) = %d",maxgen1);
	fprintf(fpt,"\nCrossover probability (pcross1) = %.2f",pcross1);
	fprintf(fpt,"\nMutation probability (pmutation1) = %.2f",pmutation1);
	fprintf(fpt,"\n");
	fprintf(fpt,"\n");
}

void initreport11(void)
{
	fprintf(fpt,"\n");
	fprintf(fpt,"\n     Initial Generation Statistics      ");
    fprintf(fpt,"\n     -----------------------------      ");
	fprintf(fpt,"\n");
	fprintf(fpt,"\nInitial population maximum fitness = %.2f",max);
	fprintf(fpt,"\nInitial population average fitness = %.2f",avg);
	fprintf(fpt,"\nInitial population minimum fitness = %.2f",min);
	fprintf(fpt,"\nInitial population sum of fitness  = %.2f",sumfitness);
	fprintf(fpt,"\n");
}


void initreport2(void)
{


	fprintf(fpr,"\n     Second Genetic Algorithm Parameters      ");
	fprintf(fpr,"\n     -----------------------------------      ");
	fprintf(fpr,"\n");
	fprintf(fpr,"\nPopulation size (popsize) = %d",popsize2);
	fprintf(fpr,"\nNumber of asset to be chosen (maxselect) = %d",maxselect);
	fprintf(fpr,"\nMax. generations (maxgen) = %d",maxgen2);
	fprintf(fpr,"\nCrossover probability (pcross) = %.2f",pcross2);
	fprintf(fpr,"\nMutation probability (pmutation) = %.2f",pmutation2);
	fprintf(fpr,"\n");
}


void initreport21(void)
{
	fprintf(fpr,"\n");
	fprintf(fpr,"\n     Initial Generation Statistics      ");
    fprintf(fpr,"\n     -----------------------------      ");
	fprintf(fpr,"\n");
	fprintf(fpr,"\nInitial population maximum fitness = %.2f",max);
	fprintf(fpr,"\nInitial population average fitness = %.2f",avg);
	fprintf(fpr,"\nInitial population minimum fitness = %.2f",min);
	fprintf(fpr,"\nInitial population sum of fitness  = %.2f",sumfitness);
	fprintf(fpr,"\n");
}


Boolean check1(chromosome chrom1, population1 popp, int i)
{
	int j,j1;
	Boolean ttag,ttag1;

	j=0;
	ttag1=False;
	do
	{
		ttag=True;
		j1=0;
		do 
		{
			if (chrom1[j1]!=popp->chrom[j1]) ttag=False;
			j1++;
		}
		while (ttag==True && j1<maxselect*7);
		if (ttag==True)
		{
			ttag1=True;
		}
		popp++;
		j++;
	}
	while (ttag1==False && j<i);

	return(ttag1);
}


Boolean check2(chromosome chrom1, population2 popp, int i)
{
	int j,j1;
	Boolean ttag,ttag1;

	j=0;
	ttag1=False;
	do
	{
		ttag=True;
		j1=0;
		do 
		{
			if (chrom1[j1]!=popp->chrom[j1]) ttag=False;
			j1++;
		}
		while (ttag==True && j1<maxselect*7);
		if (ttag==True) {ttag1=True;}
		popp++;
		j++;
	}
	while (ttag1==False && j<i);

	return(ttag1);
}


void initpop1(void)

/* Initialize a population at random */

{
	int j,j1,mcount,pos;
	population1 npop;

    npop=oldpop1;
		
	fprintf(fpt,"\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
	fprintf(fpt,"\n&                Initial Poulation Generation\n");               
	fprintf(fpt,"\nRow      Chromosome     Fitness Value     Duplication");       
	fprintf(fpt,"\n***************************************************************");
	for (j=0;j<popsize1;j++) {
		j1=0;
		do 
		{
			pos=rnd(0,lchrom1-1);
			if (oldpop1->chrom[pos]==False)
			{
				oldpop1->chrom[pos]=True;
				j1++;
			}
		}
		while (j1<maxselect);
		
		if (check1(oldpop1->chrom,npop,j)==True && j!=0) oldpop1->tag=False;
		else oldpop1->tag=True;

		oldpop1->fitness=objfunc1(oldpop1);		  // Evaluate initial fitness
		fprintf(fpt,"\n");
		fprintf(fpt,"%2d- ",j+1);
		writechrom(fpt,oldpop1->chrom,lchrom1);
		fprintf(fpt,"     %10.5f   ",oldpop1->fitness);
		if (oldpop1->tag==True) fprintf(fpt,"No");
		else fprintf(fpt,"Yes");

		oldpop1++;
	}
	oldpop1=npop;
}


void initpop2(void)

/* Initialize a population at random */

{
	int i,j;
	population2 npop;

    npop=oldpop2;

	fprintf(fpr,"\n***************************************************************");
	fprintf(fpr,"\n                Initial Poulation Generation\n");
	fprintf(fpr,"\nRow      Chromosome     Fitness Value     Duplication");
	fprintf(fpr,"\n***************************************************************");

	for (i=0;i<popsize2;i++)
	{
		for (j=0;j<maxselect*7;j++) oldpop2->chrom[j]=FlipCoin(0.5);
		
		decode(oldpop2);
		if (check2(oldpop2->chrom,npop,i)==True && i!=0) oldpop2->tag=False;
		else oldpop2->tag=True;

		oldpop2->fitness=objfunc2(oldpop2);		  // Evaluate initial fitness
		fprintf(fpr,"\n");
		fprintf(fpr,"%2d- ",i+1);
		writechrom(fpr,oldpop2->chrom,maxselect*7);
		fprintf(fpr,"     %10.5f   ",oldpop2->fitness);
		if (oldpop2->tag==True) fprintf(fpr,"No");
		else fprintf(fpr,"Yes");

		oldpop2++;
		
	}
	oldpop2=npop;
}



void mutation1(population1 pop1, population1 pop2)
{

	int jmute,i,j;

	
	jmute=rnd(0,lchrom1-1);
	fprintf(fpt,"swappoint = %4d    ", jmute);
	nmutation++;
	j=0;
	for (i=jmute;i<lchrom1;i++) 
	{
		pop2->chrom[j]=pop1->chrom[i];
		j++;
	}

	for (i=0;i<jmute;i++)
	{
		pop2->chrom[j]=pop1->chrom[i];
		j++;
	}
}


void mutation2(population2 pop1, population2 pop2)
{
	int i;

	for (i=0;i<maxselect*7;i++)
		if (FlipCoin(pmutation2)==True) 
		{
			nmutation++;
			if (pop1->chrom[i]==True) pop2->chrom[i]=False;
			else pop2->chrom[i]=True;
		} else pop2->chrom[i]=pop1->chrom[i];
}



void crossover1(population1 parent1, population1 parent2, population1 child1, population1 child2, int *jcross1, int *jcross2)
{

	int j,mcount,j1,pos;
	
	for(j=0;j<=*jcross1;j++)
	{
		child1->chrom[j]=parent2->chrom[j];
		child2->chrom[j]=parent1->chrom[j];
	}

	for(j=*jcross1+1;j<=*jcross2;j++)
	{
		child1->chrom[j]=parent1->chrom[j];
		child2->chrom[j]=parent2->chrom[j];
	}

	for(j=*jcross2+1;j<lchrom1;j++)
	{
		child1->chrom[j]=parent2->chrom[j];
		child2->chrom[j]=parent1->chrom[j];
	}

	mcount=0;
	for (j1=0;j1<lchrom1;j1++)	
		if (child1->chrom[j1]==True)	mcount++;

	if (mcount>maxselect)
		while (mcount!=maxselect)
		{
			pos=rnd(0,lchrom1-1);
			if (child1->chrom[pos]==True) {
					child1->chrom[pos]=False;
					mcount--;
				}
		}
	
	if (mcount<maxselect)
		while (mcount!=maxselect)
		{
			pos=rnd(0,lchrom1-1);
			if (child1->chrom[pos]==False) {
					child1->chrom[pos]=True;
					mcount++;
				}
		}

	mcount=0;
	for (j1=0;j1<lchrom1;j1++)	
		if (child2->chrom[j1]==True)	mcount++;

	if (mcount>maxselect)
		while (mcount!=maxselect)
		{
			pos=rnd(0,lchrom1-1);
			if (child2->chrom[pos]==True) {
					child2->chrom[pos]=False;
					mcount--;
				}
		}
	
	if (mcount<maxselect)
		while (mcount!=maxselect)
		{
			pos=rnd(0,lchrom1-1);
			if (child2->chrom[pos]==False) {
					child2->chrom[pos]=True;
					mcount++;
				}
		}
}



void crossover2(population2 parent1, population2 parent2, population2 child1, population2 child2, int *jcross1, int *jcross2)
{

	int j;
		
	for(j=0;j<=*jcross1;j++)
		{
			child1->chrom[j]=parent2->chrom[j];
			child2->chrom[j]=parent1->chrom[j];
		}

		for(j=*jcross1+1;j<=*jcross2;j++) {
			child1->chrom[j]=parent1->chrom[j];
			child2->chrom[j]=parent2->chrom[j];
		}

		for(j=*jcross2+1;j<maxselect*7;j++) {
			child1->chrom[j]=parent2->chrom[j];
			child2->chrom[j]=parent1->chrom[j];
		}

}



void generation1(int gen)
{
	int k,i, j=0,j1, mate1, mate2, jcross1,jcross2;
	population1 npop1, npop2, popp, npop;

	j=0;
	popp=poptemp1;
	npop=oldpop1;
	for (i=0;i<popsize1;i++) {
		for (j1=0;j1<lchrom1;j1++)    	popp->chrom[j1]=npop->chrom[j1];
		popp->fitness=npop->fitness;
		popp->tag=npop->tag;
		popp++;
		npop++;
		j++;
	}

	i=0;
	npop=oldpop1;
	npop1=popp;
	npop2=popp+1;
	fprintf(fpt,"\n\n********************************************************************");
	fprintf(fpt,"\n                         Crossover Operator");
	fprintf(fpt,"\n********************************************************************");
	do {
		if (FlipCoin(pcross1)==True) { 
			jcross1=rnd(0,lchrom1/2-1);
			jcross2=rnd(jcross1+1,lchrom1-1);
			ncross++;
		} else
		{
			jcross1=lchrom1-1;
			jcross2=lchrom1-1;
		}

		if (jcross1!=lchrom1-1 || jcross2!=lchrom1-1)
		{
			mate1=select1(oldpop1);
			mate2=select1(oldpop1);
			crossover1((oldpop1+mate1), (oldpop1+mate2), npop1, npop2, &jcross1, &jcross2);
			npop1->fitness=objfunc1(npop1);
			if (check1(npop1->chrom,poptemp1,j)==True) npop1->tag=False;
			else npop1->tag=True;

			
			npop2->fitness=objfunc1(npop2);
			if (check1(npop2->chrom,poptemp1,j+1)==True) npop2->tag=False;
			else npop2->tag=True;

			fprintf(fpt,"\n\nparent1  ");
			writechrom(fpt,(oldpop1+mate1)->chrom,lchrom1);
			fprintf(fpt,"    Fitness= %10.5f",(oldpop1+mate1)->fitness);
			fprintf(fpt,"\nparent2  ");
			writechrom(fpt,(oldpop1+mate2)->chrom,lchrom1);
			fprintf(fpt,"    Fitness= %10.5f",(oldpop1+mate2)->fitness);
			fprintf(fpt,"\ncross1= %2d   cross2= %2d ",jcross1+1,jcross2+1);
			fprintf(fpt,"\nchild1   ");
			writechrom(fpt,npop1->chrom,lchrom1);
			fprintf(fpt,"    Fitness= %10.5f",npop1->fitness);
			fprintf(fpt,"   Duplication = ");
			if (npop1->tag==True) fprintf(fpt,"No");
			else fprintf(fpt,"Yes");

			fprintf(fpt,"\nchild2   ");
			writechrom(fpt,npop2->chrom,lchrom1);
			fprintf(fpt,"    Fitness= %10.5f",npop2->fitness);
			fprintf(fpt,"   Duplication = ");
			if (npop2->tag==True) fprintf(fpt,"No");
			else fprintf(fpt,"Yes");

			j+=2;
			npop1=npop2+1;
			npop2=npop1+1;	
		}
		i+=2;			
	}
	while (i<popsize1);

	j1=0;
	npop=oldpop1;
	popp=npop1;

	fprintf(fpt,"\n\n*************************************************************************");
	fprintf(fpt,"\n                         Mutation Operator");
	fprintf(fpt,"\n    parent         Fitness       child         Fitness      Duplication");
	fprintf(fpt,"\n***************************************************************************");

	do
	{
		if (FlipCoin(pmutation1)==True) 
		{
			nmutation++;
			for (k=0;k<lchrom1;k++)    	popp->chrom[k]=npop->chrom[k];
			fprintf(fpt,"\n");
			writechrom(fpt,npop->chrom,lchrom1);
			fprintf(fpt,"   %10.5f",npop->fitness);
			fprintf(fpt,"    ");
			mutation1(npop,popp);
			writechrom(fpt,popp->chrom,lchrom1);
			popp->fitness=objfunc1(popp);
			fprintf(fpt,"   %10.5f",popp->fitness);
			fprintf(fpt,"    ");
			if (check1(popp->chrom,poptemp1,j)==True) popp->tag=False;
			else popp->tag=True;

			if (popp->tag==True) fprintf(fpt,"No");
			else fprintf(fpt,"Yes");

			popp++;
			j++;
		}
		npop++;
		j1++;
	}
	while (j1<popsize1);

	npop1=poptemp1;
	for(i=0;i<(popsize1*3)-1;i++)
	{
		npop2=npop1;
		for (j1=i;j1<popsize1*3;j1++)
		{
			if (npop1->fitness<npop2->fitness)
			{
				for (j=0;j<lchrom1;j++)
				{
					popp11->chrom[j]=npop1->chrom[j];
					npop1->chrom[j]=npop2->chrom[j];
					npop2->chrom[j]=popp11->chrom[j];
				}
				
				popp11->fitness=npop1->fitness;
				npop1->fitness=npop2->fitness;
				npop2->fitness=popp11->fitness;

				popp11->tag=npop1->tag;
				npop1->tag=npop2->tag;
				npop2->tag=popp11->tag;
			
			}
			npop2++;
		}
		npop1++;
	}

	npop1=poptemp1;
	npop2=newpop1;
	j1=0;
	i=0;
	do 
	{
		for (j=0;j<lchrom1;j++)    	npop2->chrom[j]=npop1->chrom[j];
		npop2->fitness=npop1->fitness;
		npop2->tag=npop1->tag;
		npop2++;
		j1++;
		for (j=0;j<lchrom1;j++)    	npop1->chrom[j]=False;
		npop1->fitness=0;
		npop1->tag=False;
		npop1++;
		i++;
	}
	while (j1<popsize1);
	for (j=i;j<popsize1*3;j++) 
	{
		for (j1=0;j1<lchrom1;j1++)    	npop1->chrom[j1]=False;
		npop1->fitness=0;
		npop1->tag=False;
		npop1++;
	}
	
	for (j=0;j<lchrom1;j++) 
			mpop1->chrom[j]=newpop1->chrom[j];
		mpop1->fitness=newpop1->fitness;

		
	fprintf(fpt,"\n\n***************************************************************");
	fprintf(fpt,"\n                        Generation  %3d",gen+1);
	fprintf(fpt,"\nRow            Chromosome            Fitness Value  Duplication");
	fprintf(fpt,"\n***************************************************************");

	npop1=newpop1;
	for (i=0;i<popsize1;i++)
	{
		fprintf(fpt,"\n");
		fprintf(fpt,"%2d- ",i+1);
		writechrom(fpt,npop1->chrom,lchrom1);
		fprintf(fpt,"  %10.5f   ",npop1->fitness);
		if (npop1->tag==True) fprintf(fpt,"No");
		else fprintf(fpt,"Yes");
		npop1++;
	}

}



void generation2(int gen)
{
	int i,j, j1, jcross1,jcross2,mcount, pos;
	population2 npop1, npop2, popp, npop;

	j=0;
	popp=poptemp2;
	npop=oldpop2;
	for (i=0;i<popsize2;i++) {
		for (j1=0;j1<maxselect*7;j1++)    	popp->chrom[j1]=npop->chrom[j1];
		for (j1=0;j1<maxselect;j1++)
		{
			popp->percent[j1]=npop->percent[j1];
			popp->code[j1]=npop->code[j1];
		}
		popp->fitness=npop->fitness;
		popp->sumprct=npop->sumprct;
		popp->tag=npop->tag;
		popp++;
		npop++;
		j++;
	}

	i=0;
	npop=oldpop2;
	npop1=popp;
	npop2=popp+1;

	fprintf(fpr,"\n\n************************************************************************");
	fprintf(fpr,"\n                               Crossover Operator");
	fprintf(fpr,"\n************************************************************************");

	do {

		if (FlipCoin(pcross2)==True)
		{ 
			jcross1=rnd(0,(maxselect*7)-1);
			jcross2=rnd(jcross1+1,maxselect*7-1);

			ncross++;
		} 
		else
		{
			jcross1=maxselect*7-1;
			jcross2=maxselect*7-1;


		}

		if (jcross1!=maxselect*7-1&& jcross2!=maxselect*7-1)
		{
			crossover2(npop, npop+1, npop1, npop2, &jcross1,&jcross2);
	
			decode(npop1);
			if (npop1->sumprct==0) {
				pos=rnd(0,maxselect*7-1);
				npop1->chrom[pos]=True;
				decode(npop1);
				}
			npop1->fitness=objfunc2(npop1);
			if (check2(npop1->chrom,poptemp2,j)==True) npop1->tag=False;
			else npop1->tag=True;
			decode(npop2);
			if (npop2->sumprct==0) {
				pos=rnd(0,maxselect*7-1);
				npop2->chrom[pos]=True;
				decode(npop2);
				}
			npop2->fitness=objfunc2(npop2);
			if (check2(npop2->chrom,poptemp2,j+1)==True) npop2->tag=False;
			else npop2->tag=True;

			fprintf(fpr,"\n\nparent1  ");
			writechrom(fpr,npop->chrom,maxselect*7);
			fprintf(fpr,"    Fitness= %10.5f",npop->fitness);
			fprintf(fpr,"\nparent2  ");
			writechrom(fpr,(npop+1)->chrom,maxselect*7);
			fprintf(fpr,"    Fitness= %10.5f",(npop+1)->fitness);
			fprintf(fpr,"\ncross1= %2d   ",jcross1+1);
			fprintf(fpr,"\nchild1   ");
			writechrom(fpr,npop1->chrom,maxselect*7);
			fprintf(fpr,"    Fitness= %10.5f",npop1->fitness);
			fprintf(fpr,"   Duplication = ");
			if (npop1->tag==True) fprintf(fpr,"No");
			else fprintf(fpr,"Yes");

			fprintf(fpr,"\nchild2   ");
			writechrom(fpr,npop2->chrom,maxselect*7);
			fprintf(fpr,"    Fitness= %10.5f",npop2->fitness);
			fprintf(fpr,"   Duplication = ");
			if (npop2->tag==True) fprintf(fpr,"No");
			else fprintf(fpr,"Yes");

			npop1=npop2+1;
			npop2=npop1+1;	
			j+=2;
		}
		i+=2;			
		npop=npop+2;
	}
	while (i<popsize2);

	j1=0;
	npop=oldpop2;
	popp=npop1;

	fprintf(fpr,"\n\n**********************************************************************************************************************************************************************");
	fprintf(fpr,"\n                                      Mutation Operator");
	fprintf(fpr,"\n                              parent                               Fitness                                 child                                Fitness   Duplication");
	fprintf(fpr,"\n*********************************************************************************************************************************************************************");

	do
	{
		fprintf(fpr,"\n");
		writechrom(fpr,npop->chrom,maxselect*7);
		fprintf(fpr,"   %10.5f",npop->fitness);
		fprintf(fpr,"    ");
		mutation2(npop,popp);

		decode(popp);
		if (popp->sumprct==0) {
			pos=rnd(0,maxselect*7-1);
			popp->chrom[pos]=True;
			decode(popp);
		}
		writechrom(fpr,popp->chrom,maxselect*7);
		popp->fitness=objfunc2(popp);

		fprintf(fpr,"   %10.5f",popp->fitness);
		fprintf(fpr,"    ");

		if (check2(popp->chrom,poptemp2,j)==True) popp->tag=False;
		else popp->tag=True;

		if (popp->tag==True) fprintf(fpr,"No");
		else fprintf(fpr,"Yes");
		
		popp++;
		j++;
		npop++;
		j1++;
	}
	while (j1<popsize2);

	npop1=poptemp2;
	for(i=0;i<(popsize2*3)-1;i++)
	{
		
		npop2=npop1;
		for (j1=i;j1<popsize2*3;j1++)
		{
			if (npop1->fitness>npop2->fitness)
			{
				for (j=0;j<maxselect*7;j++)
				{
					popp12->chrom[j]=npop1->chrom[j];
					npop1->chrom[j]=npop2->chrom[j];
					npop2->chrom[j]=popp12->chrom[j];
				}
				for (j=0;j<maxselect;j++)
				{
					popp12->percent[j]=npop1->percent[j];
					npop1->percent[j]=npop2->percent[j];
					npop2->percent[j]=popp12->percent[j];
				}
				popp12->fitness=npop1->fitness;
				npop1->fitness=npop2->fitness;
				npop2->fitness=popp12->fitness;

				popp12->sumprct=npop1->sumprct;
				npop1->sumprct=npop2->sumprct;
				npop2->sumprct=popp12->sumprct;

				popp12->tag=npop1->tag;
				npop1->tag=npop2->tag;
				npop2->tag=popp12->tag;
			}
			npop2++;
		}
		npop1++;
	}
	


	npop1=poptemp2;
	npop2=newpop2;
	j1=0;
	i=0;
	do 
	{
		
		if (npop1->tag==True)
		{
			for (j=0;j<maxselect*7;j++)    	npop2->chrom[j]=npop1->chrom[j];
			for (j=0;j<maxselect;j++)		npop2->percent[j]=npop1->percent[j];
			npop2->fitness=npop1->fitness;
			npop2->sumprct=npop1->sumprct;
			npop2->tag=npop1->tag;
			npop2++;
			j1++;
		}
		for (j=0;j<maxselect*7;j++)    	npop1->chrom[j]=False;
		for (j=0;j<maxselect;j++)		npop1->percent[j]=0;
		npop1->fitness=0;
		npop1->sumprct=0;
		npop1->tag=False;
		npop1++;
		i++;
	}
	while (j1<popsize2 && i<popsize2*3);
	
	for (j=i;j<popsize2*3;j++) 
	{
		for (j1=0;j1<maxselect*7;j1++)    	npop1->chrom[j1]=False;
		for (j1=0;j1<maxselect;j1++)		npop1->percent[j1]=0;
		npop1->fitness=0;
		npop1->sumprct=0;
		npop1->tag=False;
		npop1++;
	}
	
	for (j=0;j<maxselect*7;j++)    	mpop2->chrom[j]=newpop2->chrom[j];
	for (j=0;j<maxselect;j++)		mpop2->percent[j]=newpop2->percent[j];
	mpop2->fitness=newpop2->fitness;
	mpop2->tag=newpop2->tag;
	mpop2->sumprct=newpop2->sumprct;

	fprintf(fpr,"\n\n********************************************************************");
	fprintf(fpr,"\n                          Generation  %3d",gen+1);
	fprintf(fpr,"\nRow               Chromosome              Fitness Value  Duplication");
	fprintf(fpr,"\n********************************************************************");

	npop1=newpop2;
	for (i=0;i<popsize2;i++)
	{
		fprintf(fpr,"\n");
		fprintf(fpr,"%2d- ",i+1);
		writechrom(fpr,npop1->chrom,maxselect*7);
		fprintf(fpr,"  %10.5f   ",npop1->fitness);
		if (npop1->tag==True) fprintf(fpr,"No");
		else fprintf(fpr,"Yes");
		npop1++;
	}
		
}


void corrcalc(void)
{

	double coef[12];
	float p,cov,l,avg1,avg2,std1,std2;
	float retrn[1000],sum1,sum2,sum11,sum22;
	int year[1000],icode[1000],num,ii,i,j,c,k,ii1,z,t;
	char code[8],yr[8],ret[8],ch;


	coef[0]=0.40;
	coef[1]=0.45;
	coef[2]=0.50;
	coef[3]=0.55;
	coef[4]=0.60;
	coef[5]=0.65;
	coef[6]=0.70;
	coef[7]=0.75;
	coef[8]=0.80;
	coef[9]=0.85;
	coef[10]=0.90;
    coef[11]=0.95;
    
	num=0;
	while (!feof(fdata))
	{
		fscanf(fdata,"%s %s %s",code,ret,yr);
		icode[num]=atoi(code);
        retrn[num]=atof(ret);
		year[num]=atoi(yr);
		num++;
	}

	for (c=0;c<num;c++)
	{

		ii=icode[c];
		k=c;
		t=c;
		while (ii==icode[k] && k!=num) k++;

		if (k==num) c=num;

		for (j=k;j<num;j++)
		{
			ii1=icode[j];
			i=0;
			sum1=0;
			sum11=0;
			sum2=0;
			sum22=0;
			c=t;
			z=j;

			while (ii==icode[c] && ii1==icode[j] && c!=num && j!=num)
			{
				while (year[c]>year[j] && ii1==icode[j] && j<num)  j++;
				while (year[c]<year[j] && ii==icode[c] && c<num)   c++;
				if (ii==icode[c] && ii1==icode[j] && year[c]==year[j])
				{

					sum1=sum1+retrn[c];
					sum2=sum2+retrn[j];
					sum11=sum11+retrn[c]*retrn[c];
					sum22=sum22+retrn[j]*retrn[j];
					i++;
					c++;
					j++;
				}
			}

			if (i!=0)
			{
				avg1=sum1/i;
				avg2=sum2/i;
				std1=sqrt(((sum11)/i)-((sum1/i)*(sum1/i)));
				std2=sqrt(((sum22)/i)-((sum2/i)*(sum2/i)));
			}

			c=t;
			j=z;
			cov=0;
			p=0;

			i=0;
			sum1=0;
			while (ii==icode[c] && ii1==icode[j] && c!=num && j!=num)
			{

				while (year[c]>year[j] && ii1==icode[j] && j!=num)  j++;
				while (year[c]<year[j] && ii==icode[c] && c!=num)   c++;
	
				if (ii==icode[c] && ii1==icode[j] && year[c]==year[j])
				{
					sum1=sum1+coef[year[c]-73]*(retrn[c]-avg1)*(retrn[j]-avg2);
					i++;
					c++;
					j++;
				}

			}
					
			if (i!=0)  cov=sum1/i;
			if (std1!=0 && std2!=0)  p=cov/(std1*std2);
            if (std1!=0 && std2!=0)  l=log10(1-(cov/(std1*std2)));
			fprintf(fcorr,"%10.4f %10.4f \n",cov,p);
			fprintf(fcorr1,"%10.4f \n",l);
			while (ii==icode[c] && c<num) c++;
			c--;
			while (ii1==icode[j] && j<num) j++;
			j--;
		}
	}

}

main()
{
	int gen=0,num=0,i,j;
    population1 ppop;
	population2 ppop1,ppop2;
	double elapsed_time;
	float lavg,epcilon;
	char ch;


    fpt = fopen( "fprintf.out", "w" );
    fpr = fopen( "fprintf1.out", "w" );
    fdata1 = fopen( "risk.txt", "r" );
	fcorr = fopen( "corr.txt", "w" );
	fcorr1=fopen("corr1.txt", "w");
    fdata = fopen( "return.txt", "r" );
    
	printf("Please Wait");

    corrcalc();

	fclose(fcorr);
	fclose(fcorr1);
	fclose(fdata);

    fdata2 = fopen( "corr.txt", "r" );
    fdata3=fopen("corr1.txt", "r");
	initdata();
		
	time(&start);
	initreport1();
	warmup_random(0.4);
	initpop1();
	statistics1(oldpop1);
	initreport11();
	
	nmutation=0;
	ncross=0;

	lavg=avg;
	epcilon=0.0001;
	do
	{	
		
		generation1(gen);
		statistics1(newpop1);

		ppop=oldpop1;
		oldpop1=newpop1;
		initreport11();
		newpop1=ppop;
		if (avg-lavg<epcilon && gen!=0) num++;
		else num=0;
		
		lavg=avg;
		gen++;
	}
	while (gen<maxgen1 && num!=15);
	
	time(&finish);

	elapsed_time=difftime(finish,start);
	printf("\n\nTotal elapsed time for first algorithm:  %6.0f seconds in %4dth generation.",elapsed_time,gen-1);

	printf("\nMaximum Fitness function for first algorithm: %.4f \nSelection : ",mpop1->fitness);
	
	j=0;
	for (i=0;i<lchrom1;i++) 
		if (mpop1->chrom[i]==True) 
		{
			printf("%3d  ",i+1);
			t[j]=i;
			j++;
		}

	ppop1=newpop2;
	ppop2=oldpop2;
	for (i=0;i<popsize2;i++)
	{
		for (j=0;j<maxselect;j++)
		{
			ppop1->code[j]=t[j];
			ppop2->code[j]=t[j];
		}
		ppop1++;
		ppop2++;
	}


	ppop1=poptemp2;
	for (i=0;i<popsize2*3;i++)
	{
		for (j=0;j<maxselect;j++)
			ppop1->code[j]=t[j];
		ppop1++;
	}


	for (j=0;j<maxselect;j++)
		{
			popp12->code[j]=t[j];
			mpop2->code[j]=t[j];
		}

	time(&start);
	initreport2();
	warmup_random(0.4);
	initpop2();
	statistics2(oldpop2);
	initreport21();
	
	nmutation=0;
	ncross=0;
	gen=0;
	num=0;
	lavg=avg;
	epcilon=0.0001;
	do
	{	
		
		generation2(gen);
		statistics2(newpop2);

		ppop1=oldpop2;
		oldpop2=newpop2;
		initreport21();
		
		newpop2=ppop1;
		if (fabs(avg-lavg)<epcilon && gen!=0) num++;
		else num=0;
		
		lavg=avg;
		gen++;
	}
	while (gen<maxgen2 && num!=15);
	
	time(&finish);

	elapsed_time=difftime(finish,start);
	printf("\n\nTotal elapsed time for second algorithm:  %6.0f seconds in %4dth generation.",elapsed_time,gen-1);

	printf("\nMaximum fitness function for second algorithm: %.4f \nPercentage : ",mpop2->fitness);
	for (i=0;i<maxselect;i++) 
		printf("%.4f  ",mpop2->percent[i]/mpop2->sumprct);

	printf("\n\nFull results are in files 'fprintf.out' and 'fprintf1.out'.\n");
		
	ch=getchar();
	ch=getchar();

	fclose(fpt);
	fclose(fpr);
	fclose(fdata1);
	fclose(fdata2);
	fclose(fdata3);
	fclose(fcorr);
	fclose(fcorr1);
	free(oldpop1);
	free(newpop1);
	free(oldpop2);
	free(newpop2);
	free(poptemp1);
	free(poptemp2);
	free(popp11);
	free(popp12);
	free(mpop1);
	free(mpop2);

	return 0;
}


