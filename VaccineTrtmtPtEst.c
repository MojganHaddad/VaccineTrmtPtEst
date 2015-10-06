/***********************************************
File: VaccineTreatmntPtEst.c
Date: 9/16/99
Coded by: Mojgan Haddad
Comments: 
************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>

#define Init -1

typedef char *string;

typedef enum{
	U, T
} distnT;

/* paramT is the type of array that will hold all the info from the parameter file that is read in. 
   Each element in the array has a name (ex. mu), a min value, 
   a peak value (irrelevant in a uniform dist'n, but use 0 as a placeholder in the input file), 
   a max value, and a type (either U or T for uniform or triangular) */
   
typedef struct{
	string name;
	double min;
	double peak;
	double max;
	distnT dType;
} paramT;

/* Since we need the parameter values and the outcome variable values ranked for the 
   sensitivity analysis, I made a small structure, such that each value can easily have a 
   corresponding rank. */
   
typedef struct{
	double value;
	int rank;
} matrixDB;

/* All function prototypes below */
static void Error(string msg);
static int RandomInteger(int low, int high);
static void Randomize(void);
static string CopyString(string s);
static void *GetBlock(size_t nbytes);
static distnT ConvertDistnToType(char dType);
static void FillHypercubeMatrix(matrixDB **lhcMatrix, paramT *pArray, int N, int K);
static void FillHypercube(int N, int element, matrixDB **lhcMatrix, paramT parameter);
static void FillUniformHypercube(int N, matrixDB *putArray, double min, double max);
static void FillTriangularHypercube(int N, matrixDB *putArray, double min, double peak, double max);
static void JumbleMatrixRow(matrixDB *row, int N);
static void GiveRankings(matrixDB* array, int N);
static int FindLowest(matrixDB* array, double least, int N);
static void GiveSameRanking(matrixDB* array, int N);
static void PrintHypercubeSamples(matrixDB **lhcMatrix, paramT *pArray, int N, int K);
static void CleanArray(double **timeMatrix, int nTimes, int boxes);
static void PrintOutTimeArray(double **timeMatrix, int nTimes, matrixDB **lhcMatrix, int sample);
static void FillTimeArray(int sample, matrixDB **lhcMatrix, double **timeMatrix, int nTimes, int boxes);
static void FillInitialConditions(int sample, matrixDB **lhcMatrix, double *firstrow);
static void FillRestThroughRungeKutta(int sample, matrixDB **lhcMatrix, double **array, int nTimes, int boxes);
static void Runge_Kutta(double *inArray, double *derivArray, int n, double x, double h, double *outArray, int sample, matrixDB **lhcMatrix);
static void GetDerivatives(double *array, double *derivArray, int sample, matrixDB **lhcMatrix);
static void FillInOutcomeVariable(int sample, double **timeMatrix, matrixDB **lhcMatrix);//, int arrayPos);
static void UncertaintyAnalysis(matrixDB **lhcMatrix, int N, int K);
static double GetVariance(matrixDB *col, int N, double mean);
static double GetMedian(matrixDB *col, int N);
static double GetTotal(matrixDB *col, int N);
static double GetCu(matrixDB **lhcMatrix, int timeBox, int sample);
/* end of function prototypes */


/* global variables- c: # of sex partners/yr
					 pi: influx of sexually active people into the population
   Both c and pi are derived for each of the N hypercube samples- they come from the
   		prevalence and popthous parameter values. */
   		 				 
double c;
double pi;


/* REAL start of program */
main()
{
	FILE *infile;
	int i, j;
	char name[20];
	char dType;
	paramT *pArray;
	matrixDB **lhcMatrix;
	double **timeMatrix;
	int K;
	int N;
	float min, peak, max;
	int nTimes, boxes;
	
	/* number of times Runge-Kutta is called- since the time increment is 0.0005 (see later in the code),
	   nTimes corresponds to 80000*0.0005 yrs -> 40 yrs */
	nTimes = 80000;
	/* the boxes correspond to:
		box0 = X
		box1 = Yv1
		box2 = Yv2
		box3 = Yw
		box4 = Yvw
		box5 = YvT
		box6 = YwT
		box7 = YvwT
		box8 = A1 from untreated Yw
		box9 = A2 from untreated Yv
		box10 = A3 from untreated Yvw
		box11 = A4 from treated Yw
		box12 = A5 from treated Yv
		box13 = A6 from treated Yvw
		box14 = Dth(t)1 from untreated Yw
		box15 = Dth(t)2 from untreated Yv
		box16 = Dth(t)3 from untreated Yvw
		box17 = Dth(t)4 from treated Yw
		box18 = Dth(t)5 from treated Yv
		box19 = Dth(t)6 from treated Yvw
	*/
	boxes = 20;
	
	// Reading in of the input file with the parameter ranges
	infile = fopen("HIVparPtEst.txt", "r");
	if(infile == NULL) Error("Cannot find file named: HIVparPtEst.txt");
	
	fscanf(infile, "LHC Samples:%d\n", &N);
	
	fscanf(infile, "Number of parameters:%d\n", &K);
	
	pArray = (paramT *) GetBlock(K*sizeof(paramT));
	
	for(i=0; i< K; i++){
		fscanf(infile, "%s %f %f %f %c\n", name, &min, &peak, &max, &dType);
		pArray[i].name = CopyString(name);
		pArray[i].min = (double) min;
		pArray[i].peak = (double) peak;
		pArray[i].max = (double) max;
		pArray[i].dType = ConvertDistnToType(dType);
	}
	
	fclose(infile);
	// END of reading in input file
	
	/* N = number of hypercube samples
	   K = number of parameters to be hypercubed */
	/* The way I stored all the hypercubed samples is via a two-dimensional array:
		(memory allocated below):
		lhcMatrix has K+1 columns (each column holds a parameter), and N rows (each row holds a lhc-ed sample)
		For example lhcMatrix[2][3] holds the third hypercubed sample of the parameter muA.
		Actually it holds the value (lhcMatrix[2][3].value) as well as its rank (lhcMatrix[2][3].rank).
		The nice thing about this setup then is that if you look at any row (across the K columns) you get the point
		in the hypercube (ie. all your parameter values for that run).
	*/ 
	
	lhcMatrix = (matrixDB **) GetBlock((K+1)*sizeof(matrixDB *));
	/* extra space in matrix for outcome variable... so after I choose an outcome variable,
	   I can place it alongside the parameter value and rank it as well */
	
	for(i = 0; i< K+1; i++){
		lhcMatrix[i] = (matrixDB *) GetBlock(N*sizeof(matrixDB));
		for(j=0; j<N; j++){
			lhcMatrix[i][j].value = Init;
			lhcMatrix[i][j].rank = Init;
		}	
	}
	
	/* Now I fill the 2-D matrix with hypercubed values and print out the samples */
	FillHypercubeMatrix(lhcMatrix, pArray, N, K);
	PrintHypercubeSamples(lhcMatrix, pArray, N, K);
	
	/* create a 2-D time matrix- nTimes x boxes,
	   such that for each time, I keep the value of boxes */
	   
	timeMatrix = (double **) GetBlock(nTimes*sizeof(double *));

	for(i = 0; i< nTimes; i++){
		timeMatrix[i] = (double *) GetBlock(boxes*sizeof(double));	
	}
	
	/* I use the same timeMatrix for each of the N simulations (too much memory otherwise),
	   so I like to clean the timeMatrix before reusing it*/
	   
	for(i=0; i < N; i++){
		CleanArray(timeMatrix, nTimes, boxes);
		
		FillTimeArray(i, lhcMatrix, timeMatrix, nTimes, boxes);
		
		/* since I overwrite the timematrix in the next turn 
		   of the loop, I choose my outcome variable
		   and place it in the lhcMatrix now (in the ith row) */
		FillInOutcomeVariable(i, timeMatrix, lhcMatrix);//, K);
		
		PrintOutTimeArray(timeMatrix, nTimes, lhcMatrix, i);
	}
	
}

static void PrintOutTimeArray(double **timeMatrix, int nTimes, matrixDB **lhcMatrix, int sample)
{
	FILE* outfile;
	int j,box;
	double deathsYw, deathsYv, deathsYwv, tot;
	// double fromYw, fromYv, fromYwv, Et, Cv, Cu, prev;
	double CumdeathsYw, CumdeathsYv, CumdeathsYwv, CumdeathsYtw;
	double CumdeathsYtvw, CumdeathsYtv, deathsYtw, deathsYtv, deathsYtvw;
	
	// printing out the annual AIDS deaths per 100,000 for WT, V, WT+V
	outfile = fopen("AidsDeathsPer100000", "a");
	if(outfile == NULL) Error("Cannot open 'AidsDeathsPer100000' for writing");
	
	fprintf(outfile, "UW-Deaths UV-Deaths UVW-Deaths TW-Deaths TV-Deaths TVW-Deaths\n");
	for(j=0; j< nTimes; j+=1000){
		deathsYw = timeMatrix[j][8]*lhcMatrix[2][sample].value;
		deathsYv = timeMatrix[j][9]*lhcMatrix[2][sample].value;
		deathsYwv = timeMatrix[j][10]*lhcMatrix[2][sample].value;
		deathsYtw = timeMatrix[j][11]*lhcMatrix[2][sample].value;
		deathsYtv = timeMatrix[j][12]*lhcMatrix[2][sample].value;
		deathsYtvw = timeMatrix[j][13]*lhcMatrix[2][sample].value;
		
		// avoid negative death rate
		if (deathsYw < 0)  deathsYw=0;
		if (deathsYv < 0)  deathsYv=0;
		if (deathsYwv < 0)  deathsYwv=0;
		if (deathsYtw < 0)  deathsYtw=0;
		if (deathsYtv < 0)  deathsYtv=0;
		if (deathsYtvw < 0)  deathsYtvw=0;	
		
		// total population
		tot = 0;
		for (box=0; box<8; box++) tot += timeMatrix[j][box];
		fprintf(outfile, "%g %g %g %g %g %g\n", deathsYw*100000/tot, deathsYv*100000/tot, deathsYwv*100000/tot, deathsYtw*100000/tot, deathsYtv*100000/tot, deathsYtvw*100000/tot);
	}
	fprintf(outfile, "\n");
	fclose(outfile);
	
	outfile = fopen("CumAidsDeaths", "a");
	if(outfile == NULL) Error("Cannot open 'CumAidsDeaths' for writing");
	
	fprintf(outfile, "CumYwDeaths CumYvDths CumYwvDths CumYtwDths CumYtvDths CumYtvwDths\n");
	for(j=0; j< nTimes; j+=1000){
		CumdeathsYw = timeMatrix[j][14];
		CumdeathsYv = timeMatrix[j][15];
		CumdeathsYwv = timeMatrix[j][16];
		CumdeathsYtw = timeMatrix[j][17];
		CumdeathsYtv = timeMatrix[j][18];
		CumdeathsYtvw = timeMatrix[j][19];
		fprintf(outfile, "%g %g %g %g %g %g\n", CumdeathsYw, CumdeathsYv, CumdeathsYwv, CumdeathsYtw, CumdeathsYtv, CumdeathsYtvw);
	}
	fprintf(outfile, "\n");
	fclose(outfile);
	
	/*outfile = fopen("PREV", "a");
	if(outfile == NULL) Error("Cannot open 'PREVyw' for writing");
	
	fprintf(outfile, "PREVywper100000: ");
	for(j=0; j< nTimes; j++){
			if(j%1000 == 0){
				prev = timeMatrix[j][2];
				tot = timeMatrix[j][0] + timeMatrix[j][1] + timeMatrix[j][2] + timeMatrix[j][3] + timeMatrix[j][4] + timeMatrix[j][6] + timeMatrix[j][7];
				fprintf(outfile, "%g ", prev*100000/tot);
			}
	}
	fprintf(outfile, "\n");
	fprintf(outfile, "PREVyvper100000: ");
	for(j=0; j< nTimes; j++){
			if(j%1000 == 0){
				prev = timeMatrix[j][1];
				tot = timeMatrix[j][0] + timeMatrix[j][1] + timeMatrix[j][2] + timeMatrix[j][3] + timeMatrix[j][4]  + timeMatrix[j][6] + timeMatrix[j][7];
				fprintf(outfile, "%g ", prev*100000/tot);
			}
	}
	fprintf(outfile, "\n");
	fprintf(outfile, "PREVywvper100000: ");
	for(j=0; j< nTimes; j++){
			if(j%1000 == 0){
				prev = timeMatrix[j][3];
				tot = timeMatrix[j][0] + timeMatrix[j][1] + timeMatrix[j][2] + timeMatrix[j][3] + timeMatrix[j][4]  + timeMatrix[j][6] + timeMatrix[j][7];
				fprintf(outfile, "%g ", prev*100000/tot);
			}
	}
	fprintf(outfile, "\n");
	fclose(outfile);
	/*
	outfile = fopen("Et", "a");
	if(outfile == NULL) Error("Cannot open 'Et' for writing");
	
	fprintf(outfile, "Et: ");
	for(j=0; j< nTimes; j++){
			if(j%1000 == 0){
				Cv = timeMatrix[j][5];
				Cu = GetCu(lhcMatrix, j, sample);
				Et = 1 - Cv/Cu;
				fprintf(outfile, "%g ", Et);
			}
	}
	fprintf(outfile, "\n");
	fclose(outfile);*/
	/* new graph */
	/*outfile = fopen("IncomingAIDS", "a");
	if(outfile == NULL) Error("Cannot open 'IncomingAIDS' for writing");
	
	fprintf(outfile, "IncomingAIDS: ");
	for(j=0; j< nTimes; j++){
			if(j%1000 == 0){
				fromYw = lhcMatrix[6][sample].value*timeMatrix[j][2];
				fromYv = lhcMatrix[6][sample].value*lhcMatrix[8][sample].value*timeMatrix[j][1];
				fromYwv = lhcMatrix[6][sample].value*lhcMatrix[7][sample].value*timeMatrix[j][3];
				fprintf(outfile, "%g %g %g\n", fromYw, fromYv, fromYwv);
			}
	}
	fprintf(outfile, "\n");
	fclose(outfile); */
}

static void FillTimeArray(int sample, matrixDB **lhcMatrix, double **timeMatrix, int nTimes, int boxes)
{
	/* When I fill in the time array for a specific run,
	 I get the starting equilib values and then run the 
	 scenario through time */
	
	FillInitialConditions(sample, lhcMatrix, timeMatrix[0]);
	
	FillRestThroughRungeKutta(sample, lhcMatrix, timeMatrix, nTimes, boxes);
}

static void FillInitialConditions(int sample, matrixDB **lhcMatrix, double *firstrow)
{
	double p, mu, Ro, betaW, VwU, muA, prev, pop;

	// initial conditions are at equilib without vaccination and treatment	
	p = lhcMatrix[0][sample].value;
	mu = lhcMatrix[1][sample].value;
	muA = lhcMatrix[2][sample].value;
	
	betaW = lhcMatrix[4][sample].value;
	VwU = lhcMatrix[15][sample].value;
	
	prev = lhcMatrix[22][sample].value;
	pop = lhcMatrix[23][sample].value;
	
	Ro = 1/(1-prev);
	c = Ro*(mu + VwU)/betaW;
	pi = pop * (mu + VwU*prev);

	printf("\nRunning Simulation # %d....\n", sample + 1);
	
	// Initial condition: mass vaccination, no treatment
	firstrow[0] = (1-p)*(1-prev)*pop;	// X
	firstrow[1] = p*(1-prev)*pop;	// Yv1
	firstrow[2] = 0; 		// Yv2
	firstrow[3] = pop*prev; // Yw
	firstrow[4] = 0; 		// Yvw
	firstrow[5] = 0; 		// YvT
	firstrow[6] = 0; 		// YwT
	firstrow[7] = 0; 		// YvwT
	firstrow[8] = firstrow[3]*VwU/(mu + muA); // A1 from untreated Yw 
	firstrow[9] = 0; // A2 from untreated Yv
	firstrow[10] = 0; // A3 from untreated Yvw
	firstrow[11] = 0; // A4 from treated Yw
	firstrow[12] = 0; // A5 from treated Yv
	firstrow[13] = 0; // A5 from treated Yvw
	firstrow[14] = 0; // Dth(t)1 from Yw
	firstrow[15] = 0; // Dth(t)2 from Yv
	firstrow[16] = 0; // Dth(t)3 from Yvw
	firstrow[17] = 0; // Dth(t)4 from Ytw
	firstrow[18] = 0; // Dth(t)5 from Ytv
	firstrow[19] = 0; // Dth(t)5 from Ytvw
	
}

static void FillRestThroughRungeKutta(int sample, matrixDB **lhcMatrix, double **array, int nTimes, int boxes)
{
	double *derivArray;
	double x, timeIncrement;
	int i;
	
	/* This function calls the Runge-Kutta function for 
	   each of the nTimes */
	
	derivArray = (double *) GetBlock(boxes*sizeof(double));
	
	GetDerivatives(array[0], derivArray, sample, lhcMatrix);
	
	x = 0.0;
	
	// timeIncrement in days is 0.0005 years
	timeIncrement = 0.0005;
	for(i=0; i<nTimes-1; i++){	
		
		Runge_Kutta(array[i], derivArray, 20, x, timeIncrement, array[i+1], sample, lhcMatrix);
		x += timeIncrement;
		if(i%1000 == 0) printf(".");
	}
	printf("\n");
	free(derivArray);
}

static void Runge_Kutta(double *inArray, double *derivArray, int n, double x, double h, double *outArray, int sample, matrixDB **lhcMatrix)
{
	int i;
	double xh, hh, h6;
	double *dym;
	double *dyt;
	double *yt;
	
	/* Doing the Runge-Kutta- copied from Numerical Recipes... */
	
	dym = (double *) GetBlock((n*sizeof(double)));
	dyt = (double *) GetBlock((n*sizeof(double)));
	yt = (double *) GetBlock((n*sizeof(double)));

	hh = h*0.5;
	h6 = h/6.0;
	xh = x+hh;
	
	for(i = 0; i < n; i++){
		yt[i] = inArray[i]+hh*derivArray[i];
	}
	GetDerivatives(yt, dyt, sample, lhcMatrix);
	
	for(i=0; i<n; i++){
		yt[i] = inArray[i]+hh*dyt[i];
	}
	GetDerivatives(yt, dym, sample, lhcMatrix);
	
	for(i=0; i<n; i++){
		yt[i] = inArray[i]+h*dym[i];
		dym[i] = dyt[i]+dym[i];
	}
	GetDerivatives(yt, dyt, sample, lhcMatrix);
	
	for(i=0; i<n; i++){
		outArray[i] = inArray[i]+h6*(derivArray[i]+dyt[i]+2.0*dym[i]);
	}
	
	free(dym);
	free(dyt);
	free(yt);
}


static void GetDerivatives(double *array, double *derivArray, int sample, matrixDB **lhcMatrix)
{
	double X, Yv1, Yv2, Yw, YwT, Yvw, YvT, YvwT, TOTAL;
	double AYv, AYw, AYvw, AYtw, AYtv, AYtvw, SIGMAv, SIGMAw, SIGMAvw;
	double dxdt, dyv1dt, dyv2dt, dywdt, dywtdt, dyvwdt, daYtvwdt;
	double daYvdt, daYwdt, daYvwdt, daYtwdt, daYtvdt, dyvtdt,dyvwtdt; 
	double dDth1dt, dDth2dt, dDth3dt, dDth4dt, dDth5dt, dDth6dt;
	double p, phi, mu, muA, BwU, BwT, Bv1, Bv2, BvT, BvwU, BvwT, LaV, LaW;
	double Vv, VvU, VvT, VwU, VwT, VvwU, VvwT, Fw, Fv, Fvw;
	
	
	/* get access to the lhc-ed samples stored in the 
	   2-d lhcMatrix for simualtion # 'sample' */
	p = lhcMatrix[0][sample].value;
	mu = lhcMatrix[1][sample].value;
	muA = lhcMatrix[2][sample].value;
	phi = lhcMatrix[3][sample].value;
	
	// Transmission probabilities : Bettas
	BwU = lhcMatrix[4][sample].value;
	BwT = lhcMatrix[5][sample].value;
	Bv1 = lhcMatrix[6][sample].value;
	Bv2 = lhcMatrix[7][sample].value;
	BvT = lhcMatrix[8][sample].value;
	BvwU = lhcMatrix[9][sample].value;
	BvwT = lhcMatrix[10][sample].value;
	
	// Progression rates : Vs
	// VvAIDS = lhcMatrix[11][sample].value;
	Vv = lhcMatrix[12][sample].value;
	VvU = lhcMatrix[13][sample].value;
	VwU = lhcMatrix[15][sample].value;
	VwT = lhcMatrix[16][sample].value * VwU; // VwTAlpha
	VvT = lhcMatrix[14][sample].value * VwT; // VvTAlpha
	VvwU = lhcMatrix[17][sample].value * VwU; // VvwUAlpha
	VvwT = lhcMatrix[18][sample].value * VwT; // VvwTAlpha
	
	// Fraction of treatment : Fs
	SIGMAw = lhcMatrix[19][sample].value;
	Fw = SIGMAw * (mu + VwU) / (1 - SIGMAw);
	SIGMAv = lhcMatrix[20][sample].value;
	Fv = SIGMAv * (mu + VvU) / (1 - SIGMAv);
	SIGMAvw = lhcMatrix[21][sample].value;
	Fvw = SIGMAvw * (mu + VvwU) / (1 - SIGMAvw);
		
	// get the values in each of the boxes
	X = array[0];
	Yv1 = array[1];
	Yv2 = array[2];
	Yw = array[3];
	Yvw = array[4];
	YvT = array[5];
	YwT = array[6];
	YvwT = array[7];
	AYw = array[8]; // untreated Yw
	AYv = array[9]; // untreated Yv
	AYvw = array[10]; // untreated Yvw
	AYtw = array[11];  // treated Yw
	AYtv = array[12];  // treated Yv
	AYtvw = array[13];  // treated Yvw
	/*
	Dth1 = array[13];
	Dth2 = array[14];
	Dth3 = array[15];
	Dth4 = array[16];
	Dth5 = array[17]; */

	 // total sexually active population
	TOTAL = X + Yv1 + Yv2 + Yw + YwT + Yvw + YvT + YvwT;
	
	// Transmissibility-combined : lambda Vaccine & lambda Wild-type
	LaV = (Bv1*Yv1 + Bv2*Yv2 + BvT*YvT) / TOTAL;
	LaW = (BwU*Yw + BwT*YwT + BvwU*Yvw + BvwT*YvwT) / TOTAL;
	
	
	// the model's differential equations
	dxdt = (1-p)*pi - mu*X - c*LaV*X - c*LaW*X;
	dyv1dt = p*pi + c*LaV*X - mu*Yv1 - (1-phi)*c*LaW*Yv1 - Vv*Yv1;
	dyv2dt = Vv*Yv1 - mu*Yv2 - Fv*Yv2 - (1-phi)*c*LaW*Yv2 - VvU*Yv2;
	dywdt = c*LaW*X - mu*Yw - Fw*Yw - VwU*Yw;
	dyvwdt = (1-phi)*c*LaW*Yv1 + (1-phi)*c*LaW*Yv2 - mu*Yvw - Fvw*Yvw - VvwU*Yvw;
	dyvtdt = Fv*Yv2 - mu*YvT - (1-phi)*c*LaW*YvT - VvT*YvT;
	dywtdt = Fw*Yw - mu*YwT - VwT*YwT;
	dyvwtdt = Fvw*Yvw + (1-phi)*c*LaW*YvT - mu*YvwT - VvwT*YvwT;
	// dadt = VwU*Yw + VvU*Yv2 + VvwU*Yvw + VvT*YvT + VwT*YwT + VvwT*YvwT - mu*A - muA*A;
	daYwdt = VwU*Yw - mu*AYw - muA*AYw;
	daYvdt = VvU*Yv2 - mu*AYv - muA*AYv;
	daYvwdt = VvwU*Yvw - mu*AYvw - muA*AYvw;
	daYtwdt = VwT*YwT - mu*AYtw - muA*AYtw;
	daYtvdt = VvT*YvT - mu*AYtv - muA*AYtv;
	daYtvwdt = VvwT*YvwT - mu*AYtvw - muA*AYtvw;
	dDth1dt = muA*AYw;
	dDth2dt = muA*AYv;
	dDth3dt = muA*AYvw;
	dDth4dt = muA*AYtw;
	dDth5dt = muA*AYtv;
	dDth6dt = muA*AYtvw;

	derivArray[0] = dxdt;
	derivArray[1] = dyv1dt;
	derivArray[2] = dyv2dt;
	derivArray[3] = dywdt;
	derivArray[4] = dyvwdt;
	derivArray[5] = dyvtdt;
	derivArray[6] = dywtdt;
	derivArray[7] = dyvwtdt;
	derivArray[8] = daYwdt;
	derivArray[9] = daYvdt;
	derivArray[10] = daYvwdt;
	derivArray[11] = daYtwdt;
	derivArray[12] = daYtvdt;
	derivArray[13] = daYtvwdt;
	derivArray[14] = dDth1dt;
	derivArray[15] = dDth2dt;
	derivArray[16] = dDth3dt;
	derivArray[17] = dDth4dt;
	derivArray[18] = dDth5dt;
	derivArray[19] = dDth6dt;
}  //end of GetDerivatives


static void CleanArray(double **timeMatrix, int nTimes, int boxes)
{
	int i, j;
	
	for(i = 0; i< nTimes; i++){
		for(j=0; j<boxes; j++){
			timeMatrix[i][j] = Init;
		}	
	}
}

static void *GetBlock(size_t nbytes)
{
	void *result;
	
	result = (void *)malloc(nbytes);
	if(result == NULL) Error("No memory available");
	return(result);
}

static void Error(string msg)
{
	va_list args;
	
	va_start(args, msg);
	fprintf(stderr, "Error: ");
	vfprintf(stderr, msg, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(1);
}

static string CopyString(string s)
{
	string newstr;
	
	if(s==NULL) Error("NULL string passed to CopyString");
	newstr = (string) GetBlock(strlen(s) + 1);
	strcpy(newstr, s);
	return(newstr);
}

static void Randomize(void)
{
	srand((int) time(NULL));
}

static int RandomInteger(int low, int high)
{
	int k;
	double d;
	
	d = (double) rand()/((double) RAND_MAX + 1);
	k = (int) (d*(high - low + 1));
	return(low + k);
}

static distnT ConvertDistnToType(char dType)
{
	switch(dType){
		case 'U': return U;
		case 'T': return T;
		default: Error("Check parameter file: distribution type must be 'T' or 'U'\n"); return(U);
	}
}

static void FillHypercubeMatrix(matrixDB **lhcMatrix, paramT *pArray, int N, int K)
{	
	int element;
	
	// cycle through the parameters and hypercube each one of them
	Randomize();
	for(element = 0; element< K; element++){
		printf("Reading in Parameter %d from input file...\n", element+1);
		FillHypercube(N, element, lhcMatrix, pArray[element]);
	}
}

static void FillHypercube(int N, int element, matrixDB **lhcMatrix, paramT parameter)
{
	/* hypercube by pdf, and then jumble up the samples, 
	   since the samples were initially ordered */
	switch(parameter.dType){
		case U: FillUniformHypercube(N, lhcMatrix[element], parameter.min, parameter.max); break;
		case T: FillTriangularHypercube(N, lhcMatrix[element], parameter.min, parameter.peak, parameter.max); break;
		default: Error("Parameter's distn type must be 'T' or 'U'."); break;
	}
	
	
	JumbleMatrixRow(lhcMatrix[element], N);
	// also give the samples ranks
	if(parameter.min != parameter.max){
		GiveRankings(lhcMatrix[element], N);
	} else {
		GiveSameRanking(lhcMatrix[element], N);
	}
}

static void GiveSameRanking(matrixDB* array, int N)
{
	int i;
	
	for(i=0; i< N; i++){
		array[i].rank = 1;
	}
}


static void GiveRankings(matrixDB* array, int N)
{
	double lastLowest;
	int i, index;
	
	lastLowest = Init;
	
	for(i = 0; i< N; i++){
		index = FindLowest(array, lastLowest, N);
		array[index].rank = i+1;
		lastLowest = array[index].value;
	}
}

static int FindLowest(matrixDB* array, double least, int N)
{
	double max;
	int i, index;
	
	max = 10000;
	for(i=0; i<N; i++){
		if(array[i].value < max && array[i].value > least){
			index = i;
			max = array[i].value;
		}
	}
	return(index);
}


static void FillTriangularHypercube(int N, matrixDB *putArray, double min, double peak, double max)
{
	double base, height, volumeSlice, ximin, m1, m2, m, ximax, b1, b2, b, xintoInverse;
	int i;
	
	/* dices up a triangular distn, and goes from min to max, taking the midpoint from each of the N diced-up segments */
	base = max - min;
	height = 2.0/base;
	volumeSlice = 1.0/N;
	 
	ximin = min;
	m1 = height/(peak - min);
	m2 = -height/(max - peak);
	
	b1 = -height*min/(peak-min);
	b2 = height*max/(max-peak);
	
	m = m1;
	b = b1;
	
	for(i = 0; i< N; i++){
		xintoInverse = m*ximin*ximin/2 + b*ximin + volumeSlice;
		ximax = (-b + sqrt(b*b + 2*m*xintoInverse))/m;
		putArray[i].value = (ximin + ximax)/2;
		ximin = ximax;
		if(ximin > peak){
			ximin = peak;
			m = m2;
			b = b2;
			peak = 2*max;
		}
	}
}


static void FillUniformHypercube(int N, matrixDB *putArray, double min, double max)
{
	double incr;
	int i; 
	
	/* much easier than the triangular dicing, this one dices into N equivolume parts and picks the midpoints
	   from each segment */
	if(max == min){
		for(i=0; i< N; i++){
			putArray[i].value = min;
		}
	} else {
		incr = (max-min)/N;
		
		min += incr/2;
	
		for(i = 0; i< N; i++){
			putArray[i].value = min;
			min += incr;
		}
	}
}

static void JumbleMatrixRow(matrixDB *row, int N)
{
	int i, loc;
	double temp;
	
	/* just jumbles up the ordered samples */
	for(i=0; i<N; i++){
		loc = RandomInteger(i, N-1);
		temp = row[i].value;
		row[i].value = row[loc].value;
		row[loc].value = temp;
	}
}

static void PrintHypercubeSamples(matrixDB **lhcMatrix, paramT *pArray, int N, int K)
{
	int i, j, l;
	FILE *outfile;
	
	/* just in case you want to see the samples in an output file */
	outfile = fopen("LHCmatrixValues", "w");
	if(outfile == NULL) Error("Not enough memory available.");
	
	for(l=0; l<K; l++){
		fprintf(outfile, "%s ", pArray[l].name);
	}
	fprintf(outfile, "\n");
	
	for(i=0; i<N; i++){
		for(j=0; j<K; j++){
			fprintf(outfile, "%g ", lhcMatrix[j][i].value);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}


/* Uncertainty Analysis */
static void UncertaintyAnalysis(matrixDB **lhcMatrix, int N, int K)
{
	FILE* outfile;
	double totalOutcome, meanOutcome, medianOutcome, varOutcome;
	
	GiveRankings(lhcMatrix[K], N);

	totalOutcome = GetTotal(lhcMatrix[K], N);
	meanOutcome = totalOutcome/N;
	medianOutcome = GetMedian(lhcMatrix[K], N);
	varOutcome = GetVariance(lhcMatrix[K], N, meanOutcome);

	outfile = fopen("UncertaintyAnalysis", "w");
	fprintf(outfile, "Mean: %g Median: %g Variance: %g\n", meanOutcome, medianOutcome, varOutcome);

	fclose(outfile);
}

static double GetVariance(matrixDB *col, int N, double mean)
{
	double numer, diff, var;
	int i;
	
	numer = 0;
	for(i=0; i< N; i++){
		diff = col[i].value - mean;
		numer += diff*diff;
	}
	var = numer/N;
	return(var);
}

static double GetMedian(matrixDB *col, int N)
{
	int i;
	
	for(i=0; i< N; i++){
		if(col[i].rank == N/2) return(col[i].value);
	}
	return(Init);
}

static double GetTotal(matrixDB *col, int N)
{
	int i;
	double sum;
	
	sum = 0;
	
	for(i=0; i< N; i++){
		sum += col[i].value;
	}
	return(sum);
}


static void FillInOutcomeVariable(int sample, double **timeMatrix, matrixDB **lhcMatrix)//, int arrayPos)
{
	int timePos;
	double deaths, tot;
	
	/* for this run, I chose the outcome variable to be 
	Annual AIDS deaths in the 15th year after vaccine introduction */
	
	timePos = 15/0.0005;
	deaths = timeMatrix[timePos][7]*lhcMatrix[2][sample].value;
	tot = timeMatrix[timePos][0] + timeMatrix[timePos][1] + timeMatrix[timePos][2] + timeMatrix[timePos][3] + timeMatrix[timePos][4] + timeMatrix[timePos][5] + timeMatrix[timePos][6];
	lhcMatrix[20][sample].value = deaths*100000/tot;
}

static double GetCu(matrixDB **lhcMatrix, int timeBox, int sample)
{
	double VwU, mu, muA, yw, a;
	double dCudt, Cut, prev, pop;
	
	mu = lhcMatrix[1][sample].value;
	muA = lhcMatrix[2][sample].value;
	VwU = lhcMatrix[14][sample].value;
	prev = lhcMatrix[23][sample].value;
	pop = lhcMatrix[24][sample].value;
	
	yw = prev * pop;
	a = yw*VwU/(mu + muA);
	
	dCudt = muA*a;

	Cut = dCudt*timeBox*0.0005;
	return(Cut);

}