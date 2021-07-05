/*   
!                          Program Particles 
!  mimics the behaviour of a system of particles affected by mutual forces
!
!  Final application of the Course "Parallel Computing using MPI and OpenMP"
!
!  This program is meant to be used by course participants for demonstrating
!  the abilities acquired in optimizing and parallelising programs.
!
!  The techniques learnt at the course should be extensively used in order to
!  improve the program response times as much as possible, while gaining
!  the same or very closed results, i.e. the series of final produced images
!  and statistical results.
!
!  The code implemented herewith has been written for course exercise only,
!  therefore source code, algorithms and produced results must not be 
!  trusted nor used for anything different from their original purpose. 
!  
!  Description of the program:
!  a squared grid is hit by a field whose result is the distribution of particles
!  with different properties.
!  After having been created the particles move under the effect of mutual
!  forces.
!
!  Would you please send comments to m.cremonesi@cineca.it 
!
!  Program outline:
!  1 - the program starts reading the initial values (InitGrid)
!  2 - the generating field is computed (GeneratingField)
!  3 - the set of created particles is computed (ParticleGeneration)
!  4 - the evolution of the system of particles is computed (SystemEvolution)
!
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define HANDLE_ERROR( err ) ( HandleError( err, __FILE__, __LINE__ ) )
static void HandleError(cudaError_t err, const char *file, int line) {
	if (err != cudaSuccess) {
		printf("%s in %s at line %d\n", cudaGetErrorString(err),file, line);
		exit(EXIT_FAILURE);
	}
}

struct i2dGrid {
	int EX, EY;
	// extensions in X and Y directions
	double Xs, Xe, Ys, Ye;
	// initial and final value for X and Y directions
	int *Values;
	// 2D matrix of values
} GenFieldGrid, ParticleGrid, GenFieldGridDev, ParticleGridDev;
// le ultime due sono le strutture sul device
void print_i2dGrid(struct i2dGrid g) {
	printf("i2dGrid: EX, EY = %d, %d\n", g.EX, g.EY);
	printf("         Xs, Xe = %lf, %lf; Ys, Ye = %lf, %lf\n", g.Xs, g.Xe, g.Ys, g.Ye);
}

struct particle {
	double weight, x, y, vx, vy, fx, fy;
};

void print_particle(struct particle p) {
	printf("particle: weight=%lf, x,y=(%lf,%lf), vx,vy=(%lf,%lf), fx,fy=(%lf,%lf)\n",
	          p.weight, p.x, p.y, p.vx, p.vy, p.fx, p.fy);
}

struct Population {
    int np;
	double *weight, *x, *y, *vx, *vy;
	// particles have a position and few other properties      // WHAT ABOUT FORCE ? WHY POPULATION ISN'T A COLLECTION OF PARTICLES ?
} Particles;

void print_Population(struct Population p) {
	printf("Population: np = %d\n", p.np);
}

void DumpPopulation(struct Population p, int t) {
	/*
    * save population values on file 
   */
	char fname[80];
	FILE *dump;
	sprintf(fname, "Population%4.4d.dmp\0", t);
	dump = fopen(fname, "w");
	if (dump == NULL) {
		fprintf(stderr, "Error write open file %s\n", fname);
		exit(1);
	}
	fwrite(&p.np, sizeof((int)1), 1, dump);
	fwrite(p.weight, sizeof((double)1.0), p.np, dump);
	fwrite(p.x, sizeof((double)1.0), p.np, dump);
	fwrite(p.y, sizeof((double)1.0), p.np, dump);
	fclose(dump);
}

void ParticleStats(struct Population p, int t) {
	/*
     * Stampa il peso totale della popolazione, il baricentro, ed il peso minimo e massimo delle particelle
     * nella popolazione.
	 * write a file with statistics on population 
	*/
	FILE *stats;
	double w, xg, yg, wmin, wmax;
	int i;
	if (t <= 0)
	      stats = fopen("Population.sta", "w"); else
	      stats = fopen("Population.sta", "a");
	// append new data
	if (stats == NULL) {
		fprintf(stderr, "Error append/open file Population.sta\n");
		exit(1);
	}
	w = xg = yg = 0.0;
	wmin = wmax = p.weight[0];
	for (i = 0; i < p.np; i++) {
		if (wmin > p.weight[i])
		         wmin = p.weight[i];
		if (wmax < p.weight[i])
		         wmax = p.weight[i];
		w = w + p.weight[i];
		xg = xg + (p.weight[i] * p.x[i]);
		yg = yg + (p.weight[i] * p.y[i]);
	}
	xg = xg / w;
	yg =yg / w;
	fprintf(stats, "At iteration %d particles: %d; wmin, wmax = %lf, %lf;\n",
	           t, p.np, wmin, wmax);
	fprintf(stats, "   total weight = %lf; CM = (%10.4lf,%10.4lf)\n",
	           w, xg, yg);
	fclose(stats);
}

#define index2D(i, j, LD1) i + ((j)*LD1) // element position in 2-D arrays

// Parameters
int MaxIters, MaxSteps;
double TimeBit;
// Evolution time steps
//  functions  prototypes
int rowlen(char *riga);
int readrow(char *rg, int nc, FILE *daleg);
void InitGrid(char *InputFile);
void ParticleScreen(struct i2dGrid *pgrid, struct Population p, int s);
void IntVal2ppm(int s1, int s2, int *idata, int *vmin, int *vmax, char *name);
double MaxIntVal(int s, int *a);
double MinIntVal(int s, int *a);
double MaxDoubleVal(int s, double *a);
double MinDoubleVal(int s, double *a);
void newparticle(struct particle *p, double weight, double x, double y, double vx, double vy);
//void GeneratingField(struct i2dGrid *grid, int MaxIt);
void ParticleGeneration(struct i2dGrid grid, struct i2dGrid pgrid, struct Population *pp);
void SystemEvolution(struct i2dGrid *pgrid, struct Population *pp, int mxiter);
void ForceCompt(double *f, struct particle p1, struct particle p2);
//struct copyi2dGridInitialization(struct i2dGrid *hostStruct);
// CUDA kernel
// void GeneratingField(struct i2dGrid *grid, int MaxIt);


struct i2dGrid* copyi2dGridInitialization(struct i2dGrid hostStruct) {
  /**
  * Copia i valori di inizializzazione di hostStruct, tralasciando ciò che si
  * trova nella matrice, in una nuova struttura nel device. Restituisce un puntatore
  * alla nuova matrice sul device.
  */
	struct i2dGrid* deviceStruct;
	int* deviceValues;

	HANDLE_ERROR(cudaMalloc((void**) &deviceStruct, sizeof(hostStruct)));
  	HANDLE_ERROR(cudaMemcpy( deviceStruct, &hostStruct, sizeof(struct i2dGrid), cudaMemcpyHostToDevice ));
	HANDLE_ERROR(cudaMalloc((void**) &deviceValues, sizeof(int) * hostStruct.EX * hostStruct.EY ));
	deviceStruct->Values = deviceValues;
	
  return deviceStruct;
}

void copyBacki2dGridToHost(struct i2dGrid* deviceStruct, struct i2dGrid hostStruct) {
  /**
  * Copia tutto quello che c'è in deviceStruct all'interno di hostStruct
  */
  HANDLE_ERROR(cudaMemcpy(&hostStruct, deviceStruct, 
            sizeof(struct i2dGrid) + (sizeof(int) * hostStruct.EX * hostStruct.EY ), 
            cudaMemcpyDeviceToHost ));
}


void newparticle(struct particle *p, double weight, double x, double y, double vx, double vy) {
	/*
	 * define a new object with passed parameters
	*/
	p->weight = weight;
	p->x = x;
	p->y = y;
	p->vx = vx;
	p->vy = vy;
}
void ForceCompt(double *f, struct particle p1, struct particle p2) {
	/*
	 * Compute force acting on p1 by p1-p2 interactions 
	 * 
	*/
	double force, d, d2, dx, dy;
	static double k = 0.001, tiny = (double)1.0 / (double)1000000.0;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	d2 = dx * dx + dy * dy;
	// what if particles get in touch? Simply avoid the case
	if (d2 < tiny)
	      d2 = tiny;
	force = (k * p1.weight * p2.weight) / d2;
	f[0] = force * dx / sqrt(d2);
	f[1] = force * dy / sqrt(d2);
}
void ComptPopulation(struct Population *p, double *forces) {
	/*
	 * compute effects of forces on particles in a interval time
	 * 
	*/
	int i;
	double x0, x1, y0, y1;
	for (i = 0; i < p->np; i++) {
		x0 = p->x[i];
		y0 = p->y[i];
		p->x[i] = p->x[i] + (p->vx[i] * TimeBit) +
		                (0.5 * forces[index2D(0, i, 2)] * TimeBit * TimeBit / p->weight[i]);
		p->vx[i] = p->vx[i] + forces[index2D(0, i, 2)] * TimeBit / p->weight[i];
		p->y[i] = p->y[i] + (p->vy[i] * TimeBit) +
		                (0.5 * forces[index2D(1, i, 2)] * TimeBit * TimeBit / p->weight[i]);
		p->vy[i] = p->vy[i] + forces[index2D(1, i, 2)] * TimeBit / p->weight[i];
	}
}
void InitGrid(char *InputFile) {
	/* Output:
    * GenFieldGrid, ParticleGrid initialization
    * Maxiters, TimeBit
   */
	int nv, iv;
	double dv;
	char filerow[80];
	FILE *inpunit;
	fprintf(stdout, "Initializing grids ...\n");
	inpunit = fopen(InputFile, "r");
	if (!inpunit) {
		fprintf(stderr, "!!!! Error read access to file %s\n", InputFile);
		exit(-1);
	}
	// Now read measured values; they are read in the following order:
	// GenFieldGrid.EX, GenFieldGrid.EY,
	// GenFieldGrid.Xs, GenFieldGrid.Xe, GenFieldGrid.Ys, GenFieldGrid.Ye
	// ParticleGrid.Xs, ParticleGrid.Xe, ParticleGrid.Ys, ParticleGrid.Ye
	nv = 0;
	iv = 0;
	dv = 0.0;
	while (1) {
		if (readrow(filerow, 80, inpunit) < 1) {
			fprintf(stderr, "Error reading input file\n");
			exit(-1);
		}
		if (filerow[0] == '#')
		         continue;
		if (nv <= 0) {
			if (sscanf(filerow, "%d", &iv) < 1) {
				fprintf(stderr, "Error reading EX from string\n");
				exit(-1);
			}
			GenFieldGrid.EX = iv;
			nv = 1;
			continue;
		}
		if (nv == 1) {
			if (sscanf(filerow, "%d", &iv) < 1) {
				fprintf(stderr, "Error reading EY from string\n");
				exit(-1);
			}
			GenFieldGrid.EY = iv;
			nv++;
			continue;
		}
		if (nv == 2) {
			if (sscanf(filerow, "%lf", &dv) < 1) {
				fprintf(stderr, "Error reading GenFieldGrid.Xs from string\n");
				exit(-1);
			}
			GenFieldGrid.Xs = dv;
			nv++;
			continue;
		}
		if (nv == 3) {
			if (sscanf(filerow, "%lf", &dv) < 1) {
				fprintf(stderr, "Error reading GenFieldGrid.Xe from string\n");
				exit(-1);
			}
			GenFieldGrid.Xe = dv;
			nv++;
			continue;
		}
		if (nv == 4) {
			if (sscanf(filerow, "%lf", &dv) < 1) {
				fprintf(stderr, "Error reading GenFieldGrid.Ys from string\n");
				exit(-1);
			}
			GenFieldGrid.Ys = dv;
			nv++;
			continue;
		}
		if (nv == 5) {
			if (sscanf(filerow, "%lf", &dv) < 1) {
				fprintf(stderr, "Error reading GenFieldGrid.Ye from string\n");
				exit(-1);
			}
			GenFieldGrid.Ye = dv;
			nv++;
			continue;
		}
		if (nv <= 6) {
			if (sscanf(filerow, "%d", &iv) < 1) {
				fprintf(stderr, "Error reading ParticleGrid.EX from string\n");
				exit(-1);
			}
			ParticleGrid.EX = iv;
			nv++;
			continue;
		}
		if (nv == 7) {
			if (sscanf(filerow, "%d", &iv) < 1) {
				fprintf(stderr, "Error reading ParticleGrid.EY from string\n");
				exit(-1);
			}
			ParticleGrid.EY = iv;
			nv++;
			continue;
		}
		if (nv == 8) {
			if (sscanf(filerow, "%lf", &dv) < 1) {
				fprintf(stderr, "Error reading ParticleGrid.Xs from string\n");
				exit(-1);
			}
			ParticleGrid.Xs = dv;
			nv++;
			continue;
		}
		if (nv == 9) {
			if (sscanf(filerow, "%lf", &dv) < 1) {
				fprintf(stderr, "Error reading ParticleGrid.Xe from string\n");
				exit(-1);
			}
			ParticleGrid.Xe = dv;
			nv++;
			continue;
		}
		if (nv == 10) {
			if (sscanf(filerow, "%lf", &dv) < 1) {
				fprintf(stderr, "Error reading ParticleGrid.Ys from string\n");
				exit(-1);
			}
			ParticleGrid.Ys = dv;
			nv++;
			continue;
		}
		if (nv == 11) {
			if (sscanf(filerow, "%lf", &dv) < 1) {
				fprintf(stderr, "Error reading ParticleGrid.Ye from string\n");
				exit(-1);
			}
			ParticleGrid.Ye = dv;
			break;
		}
	}
	/*
        Now read MaxIters
      */
	MaxIters = 0;
	while (1) {
		if (readrow(filerow, 80, inpunit) < 1) {
			fprintf(stderr, "Error reading MaxIters from input file\n");
			exit(-1);
		}
		if (filerow[0] == '#' || rowlen(filerow) < 1)
		         continue;
		if (sscanf(filerow, "%d", &MaxIters) < 1) {
			fprintf(stderr, "Error reading MaxIters from string\n");
			exit(-1);
		}
		printf("MaxIters = %d\n", MaxIters);
		break;
	}
	/*
        Now read MaxSteps
      */
	MaxSteps = 0;
	while (1) {
		if (readrow(filerow, 80, inpunit) < 1) {
			fprintf(stderr, "Error reading MaxSteps from input file\n");
			exit(-1);
		}
		if (filerow[0] == '#' || rowlen(filerow) < 1)
		         continue;
		if (sscanf(filerow, "%d", &MaxSteps) < 1) {
			fprintf(stderr, "Error reading MaxSteps from string\n");
			exit(-1);
		}
		printf("MaxSteps = %d\n", MaxSteps);
		break;
	}
	/*
      ! Now read TimeBit
      */
	TimeBit = 0;
	while (1) {
		if (readrow(filerow, 80, inpunit) < 1) {
			fprintf(stderr, "Error reading TimeBit from input file\n");
			exit(-1);
		}
		if (filerow[0] == '#' || rowlen(filerow) < 1)
		         continue;
		if (sscanf(filerow, "%lf", &TimeBit) < 1) {
			fprintf(stderr, "Error reading TimeBit from string\n");
			exit(-1);
		}
		printf("TimeBit = %lf\n", TimeBit);
		break;
	}
	fclose(inpunit);
	// Grid allocations
	// RISERVARE SPAZIO SU GPU (SHARED MEMORY o GLOBAL MEMORY o ALTRO ?) 
	// GenFieldGrid è 2000 x 2400 double = 38,4 MiB. Troppo grande per la shared memory
	iv = GenFieldGrid.EX * GenFieldGrid.EY;
	//cudaMalloc(&(GenFieldGridDev.Values), iv * sizeof(1));
	GenFieldGrid.Values = (int*) malloc(iv * sizeof(1));
	if (GenFieldGrid.Values == NULL) {
		fprintf(stderr, "Error allocating GenFieldGrid.Values \n");
		exit(-1);
	}
	// ParticleGrid è 1000 x 1200 double = 9.6 MiB
	iv = ParticleGrid.EX * ParticleGrid.EY;
	ParticleGrid.Values = (int*) malloc(iv * sizeof(1));
	if (ParticleGrid.Values == NULL) {
		fprintf(stderr, "Error allocating ParticleGrid.Values \n");
		exit(-1);
	}
	fprintf(stdout, "GenFieldGrid ");
	print_i2dGrid(GenFieldGrid);
	fprintf(stdout, "ParticleGrid ");
	print_i2dGrid(ParticleGrid);
	return;
}
// end InitGrid
__global__ void GeneratingField(struct i2dGrid *grid, int MaxIt) {
	/*
   !  Compute "generating" points 
   !  Output:
   !    *grid.Values 
   */
	int ix, iy, iz;
	double ca, cb, za, zb;
	double rad, zan, zbn;
	double Xinc, Yinc, Sr, Si, Ir, Ii;
	//int izmn, izmx;
	int Xdots, Ydots;
	int offset;
	//fprintf(stdout, "Computing generating field ...\n");
	Xdots = grid->EX;
	Ydots = grid->EY;
	// numero di intersezioni lungo X e Y
	Sr = grid->Xe - grid->Xs;
	// distanza tra i punti della griglia lungo X
	Si = grid->Ye - grid->Ys;
	// distanza tra i punti della griglia lungo Y
	Ir = grid->Xs;
	// coord. del punto della griglia più vicino all'origine lungo X
	Ii = grid->Ys;
	// coord. del punto della griglia più vicino all'origine lungo Y
	Xinc = Sr / (double)Xdots;
	Yinc = Si / (double)Ydots;
	//izmn=9999; izmx=-9;
	iy = blockIdx.y * blockDim.y + threadIdx.y;
	ix = blockIdx.x * blockDim.x + threadIdx.x;
	//offset = ix + iy * blockDim.x * gridDim.x;
	offset = blockDim.x * gridDim.x; //  offset totale viene poi calcolato dalla index2D
	if (iy < Ydots && ix < Xdots) {
		ca = Xinc * ix + Ir;
		// distanza dall'origine su X
		cb = Yinc * iy + Ii;
		// distanza dall'origine su Y
		rad = sqrt(ca * ca * ((double)1.0 + (cb / ca) * (cb / ca)));
		// distanza dall'origine
		zan = 0.0;
		zbn = 0.0;
		// da qui in poi, il codice assegna un valore intero 0 < x < 3000 ad ogni punto della griglia
		for (iz = 1; iz <= MaxIt; iz++) {
			if (rad > (double)2.0)
			            break;
			za = zan;
			zb = zbn;
			zan = ca + (za - zb) * (za + zb);
			zbn = 2.0 * (za * zb + cb / 2.0);
			rad = sqrt(zan * zan * ((double)1.0 + (zbn / zan) * (zbn / zan)));
		}
		//if (izmn > iz) izmn=iz;
		//if (izmx < iz) izmx=iz;
		if (iz >= MaxIt)
		            iz = 0;
		grid->Values[index2D(ix, iy, offset)] = iz;
	}
	return;
}
// end GeneratingField
void ParticleGeneration(struct i2dGrid grid, struct i2dGrid pgrid, struct Population *pp) {
	// A system of particles is generated according to the value distribution of
	// grid.Values
	int vmax, vmin, v;
	int Xdots, Ydots;
	int ix, iy, np, n;
	double p;
	Xdots = grid.EX;
	Ydots = grid.EY;
	vmax = MaxIntVal(Xdots * Ydots, grid.Values);
	vmin = MinIntVal(Xdots * Ydots, grid.Values);
	// Just count number of particles to be generated
	vmin = (double)(1 * vmax + 29 * vmin) / 30.0;
	np = 0;
	for (ix = 0; ix < Xdots; ix++) {
		for (iy = 0; iy < Ydots; iy++) {
			v = grid.Values[index2D(ix, iy, Xdots)];
			if (v <= vmax && v >= vmin)
			            np++;
		}
	}
	// allocate memory space for particles
	pp->np = np;
	pp->weight = (double*) malloc(np * sizeof((double)1.0));
	pp->x = (double*) malloc(np * sizeof((double)1.0));
	pp->y = (double*) malloc(np * sizeof((double)1.0));
	pp->vx = (double*) malloc(np * sizeof((double)1.0));
	pp->vy = (double*) malloc(np * sizeof((double)1.0));
	// Population initialization
	n = 0;
	for (ix = 0; ix < grid.EX; ix++) {
		for (iy = 0; iy < grid.EY; iy++) {
			v = grid.Values[index2D(ix, iy, Xdots)];
			if (v <= vmax && v >= vmin) {
				pp->weight[n] = v * 10.0;
				p = (pgrid.Xe - pgrid.Xs) * ix / (grid.EX * 2.0);
				pp->x[n] = pgrid.Xs + ((pgrid.Xe - pgrid.Xs) / 4.0) + p;
				p = (pgrid.Ye - pgrid.Ys) * iy / (grid.EY * 2.0);
				pp->y[n] = pgrid.Ys + ((pgrid.Ye - pgrid.Ys) / 4.0) + p;
				pp->vx[n] = pp->vy[n] = 0.0;
				// at start particles are still
				n++;
				if (n >= np)
				               break;
			}
		}
		if (n >= np)
		         break;
	}
	print_Population(*pp);
}
// end ParticleGeneration
void SystemEvolution(struct i2dGrid *pgrid, struct Population *pp, int mxiter) {
	double *forces;
	double vmin, vmax;
	struct particle p1, p2;
	double f[2];
	int i, j, t;
	// temporary array of forces
	forces = (double*) malloc(2 * pp->np * sizeof((double)1.0));
	if (forces == NULL) {
		fprintf(stderr, "Error mem alloc of forces!\n");
		exit(1);
	}
	// compute forces acting on each particle step by step
	// In ogni intervallo di tempo:
	// 1. Stampa un'immagine delle particelle nel campo
	// 2. aggiorna le forze tra ogni coppia di particelle
	// 3. poi muovi le particelle in base alle forze calcolate
	// Inoltre, scrivi in un file le statistiche del campo di particelle ogni 4 istanti di tempo
	for (t = 0; t < mxiter; t++) {
		fprintf(stdout, "Step %d of %d\n", t, mxiter);
		ParticleScreen(pgrid, *pp, t);
		// DumpPopulation call frequency may be changed
		if (t % 4 == 0)
		         DumpPopulation(*pp, t);
		ParticleStats(*pp, t);
		for (i = 0; i < 2 * pp->np; i++)
		         forces[i] = 0.0;
		for (i = 0; i < pp->np; i++) {
			newparticle(&p1, pp->weight[i], pp->x[i], pp->y[i], pp->vx[i], pp->vy[i]);
			for (j = 0; j < pp->np; j++) {
				if (j != i) {
					newparticle(&p2, pp->weight[j], pp->x[j], pp->y[j], pp->vx[j], pp->vy[j]);
					ForceCompt(f, p1, p2);
					forces[index2D(0, i, 2)] = forces[index2D(0, i, 2)] + f[0];
					forces[index2D(1, i, 2)] = forces[index2D(1, i, 2)] + f[1];
				}
			}
		}
		ComptPopulation(pp, forces);
	}
	free(forces);
}
// end SystemEvolution
void ParticleScreen(struct i2dGrid *pgrid, struct Population pp, int step) {
	// Distribute a particle population in a grid for visualization purposes
	int ix, iy, Xdots, Ydots;
	int np, n, wp;
	double rmin, rmax;
	int static vmin, vmax;
	double Dx, Dy, wint, wv;
	char name[40];
	Xdots = pgrid->EX;
	Ydots = pgrid->EY;
	for (ix = 0; ix < Xdots; ix++) {
		for (iy = 0; iy < Ydots; iy++) {
			pgrid->Values[index2D(ix, iy, Xdots)] = 0;
		}
	}
	rmin = MinDoubleVal(pp.np, pp.weight);
	rmax = MaxDoubleVal(pp.np, pp.weight);
	wint = rmax - rmin;
	Dx = pgrid->Xe - pgrid->Xs;
	Dy = pgrid->Ye - pgrid->Ys;
	for (n = 0; n < pp.np; n++) {
		// keep a tiny border free anyway
		ix = Xdots * pp.x[n] / Dx;
		if (ix >= Xdots - 1 || ix <= 0)
		         continue;
		iy = Ydots * pp.y[n] / Dy;
		if (iy >= Ydots - 1 || iy <= 0)
		         continue;
		wv = pp.weight[n] - rmin;
		wp = 10.0 * wv / wint;
		pgrid->Values[index2D(ix, iy, Xdots)] = wp;
		pgrid->Values[index2D(ix - 1, iy, Xdots)] = wp;
		pgrid->Values[index2D(ix + 1, iy, Xdots)] = wp;
		pgrid->Values[index2D(ix, iy - 1, Xdots)] = wp;
		pgrid->Values[index2D(ix, iy + 1, Xdots)] = wp;
	}
	sprintf(name, "stage%3.3d\0", step);
	if (step <= 0) {
		vmin = vmax = 0;
	}
	IntVal2ppm(pgrid->EX, pgrid->EY, pgrid->Values, &vmin, &vmax, name);
}
// end ParticleScreen
double MinIntVal(int s, int *a) {
	int v;
	int e;
	v = a[0];
	for (e = 0; e < s; e++) {
		if (v > a[e])
		         v = a[e];
	}
	return (v);
}
double MaxIntVal(int s, int *a) {
	int v;
	int e;
	v = a[0];
	for (e = 0; e < s; e++) {
		if (v < a[e])
		         v = a[e];
	}
	return (v);
}
double MinDoubleVal(int s, double *a) {
	double v;
	int e;
	v = a[0];
	for (e = 0; e < s; e++) {
		if (v > a[e])
		         v = a[e];
	}
	return (v);
}
double MaxDoubleVal(int s, double *a) {
	double v;
	int e;
	v = a[0];
	for (e = 0; e < s; e++) {
		if (v < a[e])
		         v = a[e];
	}
	return (v);
}
int rowlen(char *riga) {
	int lungh;
	char c;
	lungh = strlen(riga);
	while (lungh > 0) {
		lungh--;
		c = *(riga + lungh);
		if (c == '\0')
		         continue;
		if (c == '\40')
		         continue;
		/*  space  */
		if (c == '\b')
		         continue;
		if (c == '\f')
		         continue;
		if (c == '\r')
		         continue;
		if (c == '\v')
		         continue;
		if (c == '\n')
		         continue;
		if (c == '\t')
		         continue;
		return (lungh + 1);
	}
	return (0);
}
int readrow(char *rg, int nc, FILE *daleg) {
	//int rowlen(), lrg;
  int lrg;
	if (fgets(rg, nc, daleg) == NULL)
	      return (-1);
	lrg = rowlen(rg);
	if (lrg < nc) {
		rg[lrg] = '\0';
		lrg++;
	}
	return (lrg);
}
void IntVal2ppm(int s1, int s2, int *idata, int *vmin, int *vmax, char *name) {
	/*
      Simple subroutine to dump double data with fixed min & max values 
         in a PPM format
   */
	int i, j;
	int cm[3][256];
	/* R,G,B, Colour Map */
	FILE *ouni, *ColMap;
	int vp, vs;
	int rmin, rmax, value;
	char fname[80], jname[80], command[80];
	/*
      Define color map: 256 colours
   */
	ColMap = fopen("ColorMap.txt", "r");
	if (ColMap == NULL) {
		fprintf(stderr, "Error read opening file ColorMap.txt\n");
		exit(-1);
	}
	for (i = 0; i < 256; i++) {
		if (fscanf(ColMap, " %3d %3d %3d",
		                 &cm[0][i], &cm[1][i], &cm[2][i]) < 3) {
			fprintf(stderr, "Error reading colour map at line %d: r, g, b =", (i + 1));
			fprintf(stderr, " %3.3d %3.3d %3.3d\n", cm[0][i], cm[1][i], cm[2][i]);
			exit(1);
		}
	}
	/*
      Write on unit 700 with  PPM format
   */
	strcpy(fname, name);
	strcat(fname, ".ppm\0");
	ouni = fopen(fname, "w");
	if (!ouni) {
		fprintf(stderr, "!!!! Error write access to file %s\n", fname);
	}
	/*  Magic code */
	fprintf(ouni, "P3\n");
	/*  Dimensions */
	fprintf(ouni, "%d %d\n", s1, s2);
	/*  Maximum value */
	fprintf(ouni, "255\n");
	/*  Values from 0 to 255 */
	rmin = MinIntVal(s1 * s2, idata);
	rmax = MaxIntVal(s1 * s2, idata);
	if ((*vmin == *vmax) && (*vmin == 0)) {
		*vmin = rmin;
		*vmax = rmax;
	} else {
		rmin = *vmin;
		rmax = *vmax;
	}
	vs = 0;
	for (i = 0; i < s1; i++) {
		for (j = 0; j < s2; j++) {
			value = idata[i * s2 + j];
			if (value < rmin)
			            value = rmin;
			if (value > rmax)
			            value = rmax;
			vp = (int)((double)(value - rmin) * (double)255.0 / (double)(rmax - rmin));
			vs++;
			fprintf(ouni, " %3.3d %3.3d %3.3d", cm[0][vp], cm[1][vp], cm[2][vp]);
			if (vs >= 10) {
				fprintf(ouni, " \n");
				vs = 0;
			}
		}
		fprintf(ouni, " ");
		vs = 0;
	}
	fclose(ouni);
	// the following instructions require ImageMagick tool: comment out if not available
	strcpy(jname, name);
	strcat(jname, ".jpg\0");
	sprintf(command, "convert %s %s\0", fname, jname);
	system(command);
	return;
}

// end IntVal2ppm
int main(int argc, char *argv[]) /* FinalApplication */ {
	#include <time.h>
	   time_t t0, t1;
	time(&t0);
    fprintf(stdout, "Starting at: %s", asctime(localtime(&t0)));
	// HANDLE_ERROR(cudaMalloc((void**) &GenFieldGridDev, sizeof(1)))
	InitGrid("Particles.inp");
	
	// 	copia valori di InitGrid su GPU e alloca spazio per le matrici
	printf("Copy to device...\n");
	struct i2dGrid* GenFieldGrid_dev = copyi2dGridInitialization(GenFieldGrid);
	printf("Copy to device 2...\n");
	struct i2dGrid* ParticleGrid_dev = copyi2dGridInitialization(ParticleGrid);


	// GenFieldGrid initialization
	printf("GeneratingField...\n");
  int TILE_WIDTH = 32;
  dim3 dimBlock (TILE_WIDTH, TILE_WIDTH);
  dim3 dimGrid ((GenFieldGrid.EX-1)/TILE_WIDTH+1, (GenFieldGrid.EX-1)/TILE_WIDTH+1);
	//GeneratingField(&GenFieldGrid, MaxIters);
  GeneratingField<<<dimGrid, dimBlock>>>(GenFieldGrid_dev, MaxIters);
  copyBacki2dGridToHost(GenFieldGrid_dev, GenFieldGrid);
	// Particle population initialization
	printf("ParticleGeneration...\n");

  ParticleGeneration(GenFieldGrid, ParticleGrid, &Particles);
  //ParticleGeneration<<<dimGrid, dimBlock>>>(GenFieldGrid, ParticleGrid, &Particles);
	// Compute evolution of the particle population
	printf("SystemEvolution...\n");
	SystemEvolution(&ParticleGrid, &Particles, MaxSteps);
	time(&t1);
	fprintf(stdout, "Ending   at: %s", asctime(localtime(&t1)));
	fprintf(stdout, "Computations ended in %lf seconds\n", difftime(t1, t0));
	fprintf(stdout, "End of program!\n");
	return (0);
}
// end FinalApplication 
