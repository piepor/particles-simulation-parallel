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
!  a squared grid is hit by a field whose result is the distribution of
particles
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
#include <time.h>

//////////////////////////////       DATA            ////////////////////////////
struct i2dGrid {
    int EX, EY;            // extensions in X and Y directions
    double Xs, Xe, Ys, Ye; // initial and final value for X and Y directions
    int *Values;           // 2D matrix of values
} GenFieldGrid, ParticleGrid;

struct particle {
    double weight, x, y, vx, vy, fx, fy;
};

struct Population {
    int np;
    double *weight, *x, *y, *vx,
            *vy; // particles have a position and few other properties
} Particles;

// Parameters
int MaxIters, MaxSteps;
double TimeBit; // Evolution time steps
//////////////////////////      FUNCTION PROTOTYPES     /////////////////////////
// utility functions
int rowlen(char *riga);
int readrow(char *rg, int nc, FILE *daleg);
double MaxIntVal(int s, int *a);
double MinIntVal(int s, int *a);
double MaxDoubleVal(int s, double *a);
double MinDoubleVal(int s, double *a);
// display functions
void IntVal2ppm(int s1, int s2, int *idata, int *vmin, int *vmax, char *name);
void DumpPopulation(struct Population p, int t, char *fileName);
void print_Population(struct Population p);
void DumpForces(double *forces, int t, int np, char *fileName);
void print_i2dGrid(struct i2dGrid g);
// processing
void InitGrid(char *InputFile);
void ParticleScreen(struct i2dGrid *pgrid, struct Population p, int s);
void newparticle(struct particle *p, double weight, double x, double y,
                 double vx, double vy);
void GeneratingField(struct i2dGrid *grid, int MaxIt);
void ParticleGeneration(struct i2dGrid grid, struct i2dGrid pgrid,
                        struct Population *pp);
void SystemEvolution(struct i2dGrid *pgrid, struct Population *pp, int mxiter);
void ForceCompt(double *f, struct particle p1, struct particle p2);
__global__ void ForceCompt_par(double *f, double *x, double *y, double *weight,
                               int np);
void ComptPopulation(struct Population *p, double *forces);
__global__ void ComptPopulation_par(double *x, double *y, double *vx,
                                    double *vy, double *f, double *weight,
                                    int np, double TimeBit);

//////////////////////////      UTILITY FUNCTIONS      //////////////////////////

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
            continue; /*  space  */
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
    // int rowlen(), lrg;
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

//////////////////////////          DEFINES           ///////////////////////////

#define index2D(i, j, LD1) i + ((j)*LD1) // element position in 2-D arrays

#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))
static void HandleError(cudaError_t err, const char *file, int line) {
    if (err != cudaSuccess) {
        printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
        exit(EXIT_FAILURE);
    }
}

//////////////////////////      DISPLAY FUNCTION     ///////////////////////////

void print_i2dGrid(struct i2dGrid g) {
    printf("i2dGrid: EX, EY = %d, %d\n", g.EX, g.EY);
    printf("         Xs, Xe = %lf, %lf; Ys, Ye = %lf, %lf\n", g.Xs, g.Xe, g.Ys,
           g.Ye);
}

void DumpForces(double *forces, int t, int np, char *fileName) {
    /*
   * save population values on file
   */
    char fname[80];
    FILE *dump;

    sprintf(fname, "%s%4.4d.dmp\0", fileName, t);
    dump = fopen(fname, "wt");
    if (dump == NULL) {
        fprintf(stderr, "Error write open file %s\n", fname);
        exit(1);
    }
    fprintf(dump, "ith\tFx\tFy\n");
    for (int i = 0; i < np; i++) {
        fprintf(dump, "%d\t%lf\t%lf\n", i, forces[2 * i + 0], forces[2 * i + 1]);
    }
    fclose(dump);
}

void print_Population(struct Population p) {
    printf("Population: np = %d\n", p.np);
}

void DumpPopulation(struct Population p, int t, char *fileName) {
    /*
   * save population values on file
   */
    char fname[80];
    FILE *dump;

    sprintf(fname, "%s%4.4d.dmp\0", fileName, t);
    dump = fopen(fname, "wt");
    if (dump == NULL) {
        fprintf(stderr, "Error write open file %s\n", fname);
        exit(1);
    }
    // sostituisco fwrite con fprintf che scrive righe di testo
    /*fwrite(&p.np,sizeof((int)1),1,dump);
     fwrite(p.weight,sizeof((double)1.0),p.np,dump);
     fwrite(p.x,sizeof((double)1.0),p.np,dump);
     fwrite(p.y,sizeof((double)1.0),p.np,dump); */
    fprintf(dump, "np: %d\n", p.np);
    fprintf(dump, "ith\tweight\tx\ty\n");
    for (int i = 0; i < p.np; i++) {
        fprintf(dump, "%d\t%lf\t%lf\t%lf\n", i, p.weight[i], p.x[i], p.y[i]);
    }
    fclose(dump);
}

__global__ ParticleScreen_par():

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

void IntVal2ppm(int s1, int s2, int *idata, int *vmin, int *vmax, char *name) {
    /*
     Simple subroutine to dump double data with fixed min & max values
     in a PPM format
   */
    int i, j;
    int cm[3][256]; /* R,G,B, Colour Map */
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
        if (fscanf(ColMap, " %3d %3d %3d", &cm[0][i], &cm[1][i], &cm[2][i]) < 3) {
            fprintf(stderr, "Error reading colour map at line %d: r, g, b =",
                    (i + 1));
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
            vp =
                    (int)((double)(value - rmin) * (double)255.0 / (double)(rmax - rmin));
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
    // the following instructions require ImageMagick tool: comment out if not
    // available
    strcpy(jname, name);
    strcat(jname, ".jpg\0");
    sprintf(command, "convert %s %s\0", fname, jname);
    system(command);

    return;
}

//////////////////////////        PROCESSING         ///////////////////////////

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

    // Now read measured values; they are read in the followin order:
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
    iv = GenFieldGrid.EX * GenFieldGrid.EY;
    GenFieldGrid.Values = (int *)malloc(iv * sizeof(1));
    if (GenFieldGrid.Values == NULL) {
        fprintf(stderr, "Error allocating GenFieldGrid.Values \n");
        exit(-1);
    }
    iv = ParticleGrid.EX * ParticleGrid.EY;
    ParticleGrid.Values = (int *)malloc(iv * sizeof(1));
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

void GeneratingField(struct i2dGrid *grid, int MaxIt) {
    /*
     !  Compute "generating" points
     !  Output:
     !    *grid.Values
   */

    int ix, iy, iz;
    double ca, cb, za, zb;
    double rad, zan, zbn;
    double Xinc, Yinc, Sr, Si, Ir, Ii;
    int izmn, izmx;
    int Xdots, Ydots;
    fprintf(stdout, "Computing generating field ...\n");

    Xdots = grid->EX;
    Ydots = grid->EY;
    Sr = grid->Xe - grid->Xs;
    Si = grid->Ye - grid->Ys;
    Ir = grid->Xs;
    Ii = grid->Ys;
    Xinc = Sr / (double)Xdots;
    Yinc = Si / (double)Ydots;

    izmn = 9999;
    izmx = -9;
    for (iy = 0; iy < Ydots; iy++) {
        for (ix = 0; ix < Xdots; ix++) {
            ca = Xinc * ix + Ir;
            cb = Yinc * iy + Ii;
            rad = sqrt(ca * ca * ((double)1.0 + (cb / ca) * (cb / ca)));
            zan = 0.0;
            zbn = 0.0;
            for (iz = 1; iz <= MaxIt; iz++) {
                if (rad > (double)2.0)
                    break;
                za = zan;
                zb = zbn;
                zan = ca + (za - zb) * (za + zb);
                zbn = 2.0 * (za * zb + cb / 2.0);
                rad = sqrt(zan * zan * ((double)1.0 + (zbn / zan) * (zbn / zan)));
            }
            if (izmn > iz)
                izmn = iz;
            if (izmx < iz)
                izmx = iz;
            if (iz >= MaxIt)
                iz = 0;
            grid->Values[index2D(ix, iy, Xdots)] = iz;
        }
    }
    return;
}

void ParticleGeneration(struct i2dGrid grid, struct i2dGrid pgrid,
                        struct Population *pp) {
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
    pp->weight = (double *)malloc(np * sizeof((double)1.0));
    pp->x = (double *)malloc(np * sizeof((double)1.0));
    pp->y = (double *)malloc(np * sizeof((double)1.0));
    pp->vx = (double *)malloc(np * sizeof((double)1.0));
    pp->vy = (double *)malloc(np * sizeof((double)1.0));

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

                pp->vx[n] = pp->vy[n] = 0.0; // at start particles are still

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

void newparticle(struct particle *p, double weight, double x, double y,
                 double vx, double vy) {
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
    d2 = dx * dx +
         dy * dy; // what if particles get in touch? Simply avoid the case
    if (d2 < tiny)
        d2 = tiny;
    force = (k * p1.weight * p2.weight) / d2;
    f[0] = force * dx / sqrt(d2);
    f[1] = force * dy / sqrt(d2);
}

__global__ void ForceCompt_par(double *f, double *x, double *y, double *weight,
                               int np) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    double force, d, d2, dx, dy, fx, fy;
    static double k = 0.001, tiny = (double)1.0 / (double)1000000.0;

    if (tid < np) {
        for (int i = 0; i < np; i++) {
            if (i != tid) {
                dx = x[i] - x[tid];
                dy = y[i] - y[tid];
                d2 = dx * dx +
                     dy * dy; // what if particles get in touch? Simply avoid the case
                if (d2 < tiny)
                    d2 = tiny;
                force = (k * weight[tid] * weight[i]) / d2;
                fx = force * dx / sqrt(d2);
                fy = force * dy / sqrt(d2);
                f[0 + tid * 2] = f[0 + tid * 2] + fx;
                f[1 + tid * 2] = f[1 + tid * 2] + fy;
            }
        }
    }
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

        p->x[i] =
                p->x[i] + (p->vx[i] * TimeBit) +
                (0.5 * forces[index2D(0, i, 2)] * TimeBit * TimeBit / p->weight[i]);
        p->vx[i] = p->vx[i] + forces[index2D(0, i, 2)] * TimeBit / p->weight[i];

        p->y[i] =
                p->y[i] + (p->vy[i] * TimeBit) +
                (0.5 * forces[index2D(1, i, 2)] * TimeBit * TimeBit / p->weight[i]);
        p->vy[i] = p->vy[i] + forces[index2D(1, i, 2)] * TimeBit / p->weight[i];
    }
}

__global__ void ComptPopulation_par(double *x, double *y, double *vx,
                                    double *vy, double *f, double *weight,
                                    int np, double TimeBit) {
    /*
   * compute effects of forces on particles in a interval time
   *
   */
    int i;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    // double x0, x1, y0, y1;
    if (tid < np) {
        x[tid] = x[tid] + vx[tid] * TimeBit +
                 0.5 * f[0 + tid * 2] * TimeBit * TimeBit / weight[tid];
        vx[tid] = vx[tid] + f[0 + tid * 2] * TimeBit / weight[tid];
        y[tid] = y[tid] + vy[tid] * TimeBit +
                 0.5 * f[1 + tid * 2] * TimeBit * TimeBit / weight[tid];
        vy[tid] = vy[tid] + f[1 + tid * 2] * TimeBit / weight[tid];
    }
}

//////////////////////////           MAIN            ///////////////////////////

int main(int argc, char *argv[]) {
    time_t t0, t1;
    int TILE_WIDTH = 32;
    // dim3 dimBlock (TILE_WIDTH, TILE_WIDTH);
    // dim3 dimGrid ((GenFieldGrid.EX-1)/TILE_WIDTH+1,
    // (GenFieldGrid.EX-1)/TILE_WIDTH+1);
    const int dimBlock = TILE_WIDTH;

    time(&t0);
    fprintf(stdout, "Starting at: %s", asctime(localtime(&t0)));

    // INITIALISE ON THE CPU
    InitGrid("Particles.inp");

    // GenFieldGrid initialization
    printf("GeneratingField...\n");
    GeneratingField(&GenFieldGrid, MaxIters);
    int j = 0;
    for (int i = 0; i < GenFieldGrid.EX * GenFieldGrid.EY; i++) {
        if (GenFieldGrid.Values[i] != 0)
            j++;
    }
    printf("%d particles in the field\n", j);

    // Particle population initialization
    printf("ParticleGeneration...\n");
    ParticleGeneration(GenFieldGrid, ParticleGrid, &Particles);

    // COMPUTATION ON THE GPU
    // Compute evolution of the particle population
    printf("Malloc data on the device...\n");
    // For the part of pure computation only need particles' weight, position and
    // velocity
    double *posX_dev;
    double *posY_dev;
    double *velX_dev;
    double *velY_dev;
    double *weight_dev;
    double *forces_dev;
    double *grid_dev;
    HANDLE_ERROR(cudaMalloc((void **)&posX_dev, sizeof(double) * Particles.np));
    HANDLE_ERROR(cudaMalloc((void **)&posY_dev, sizeof(double) * Particles.np));
    HANDLE_ERROR(cudaMalloc((void **)&velX_dev, sizeof(double) * Particles.np));
    HANDLE_ERROR(cudaMalloc((void **)&velY_dev, sizeof(double) * Particles.np));
    HANDLE_ERROR(cudaMalloc((void **)&weight_dev, sizeof(double) * Particles.np));
    HANDLE_ERROR(cudaMalloc((void **)&grid_dev, sizeof(double) * ParticleGrid.Xe * ParticleGrid.Ye));
    // allocate space for the force array
    HANDLE_ERROR(
            cudaMalloc((void **)&forces_dev, 2 * sizeof(double) * Particles.np));
    // init forces to 0

    printf("Copy data on the device...\n");
    HANDLE_ERROR(cudaMemcpy(posX_dev, Particles.x, sizeof(double) * Particles.np,
                            cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(posY_dev, Particles.y, sizeof(double) * Particles.np,
                            cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(velX_dev, Particles.vx, sizeof(double) * Particles.np,
                            cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(velY_dev, Particles.vy, sizeof(double) * Particles.np,
                            cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(weight_dev, Particles.weight,
                            sizeof(double) * Particles.np,
                            cudaMemcpyHostToDevice));

    printf("SystemEvolution...\n");
    // SystemEvolution(&ParticleGrid, &Particles, MaxSteps);
    double f[2];
    double *forces;
    struct particle p1, p2;
    forces = (double *)malloc(2 * Particles.np * sizeof((double)1.0));
    memset(forces, 0.0, 2 * Particles.np * sizeof(double));
    int N = Particles.np;
    int dimGrid = (N - 1) / TILE_WIDTH + 1;
    for (int k = 0; k < MaxSteps; k++) {
        HANDLE_ERROR(cudaMemset(forces_dev, 0, 2 * sizeof(double) * Particles.np));
        ForceCompt_par<<<dimGrid, dimBlock>>>(forces_dev, posX_dev, posY_dev,
                                              weight_dev, Particles.np);
//        HANDLE_ERROR(cudaMemcpy(forces, forces_dev,
//                                sizeof(double) * 2 * Particles.np,
//                                cudaMemcpyDeviceToHost));
        ComptPopulation_par<<<dimGrid, dimBlock>>>(posX_dev, posY_dev, velX_dev,
                                                   velY_dev, forces_dev, weight_dev,
                                                   Particles.np, TimeBit);
        DumpPopulation(Particles, k, "par_Population\0");
//        DumpForces(forces, k, Particles.np, "par_forces\0");
        //ParticleScreen(&ParticleGrid, Particles, k);
        HANDLE_ERROR(cudaDeviceSynchronize());
        HANDLE_ERROR(cudaMemcpy(Particles.x, posX_dev,
                                sizeof(double) * Particles.np,
                                cudaMemcpyDeviceToHost));
        HANDLE_ERROR(cudaMemcpy(Particles.y, posY_dev,
                                sizeof(double) * Particles.np,
                                cudaMemcpyDeviceToHost));
        HANDLE_ERROR(cudaMemcpy(Particles.vx, velX_dev,
                                sizeof(double) * Particles.np,
                                cudaMemcpyDeviceToHost));
        HANDLE_ERROR(cudaMemcpy(Particles.vy, velY_dev,
                                sizeof(double) * Particles.np,
                                cudaMemcpyDeviceToHost));
        HANDLE_ERROR(cudaMemcpy(Particles.weight, weight_dev,
                                sizeof(double) * Particles.np,
                                cudaMemcpyDeviceToHost));

        //               memset(forces, 0.0, 2 * Particles.np * sizeof(double));
        //               for (int i=0; i < Particles.np; i++ ) {
        //            newparticle(&p1,Particles.weight[i],Particles.x[i],Particles.y[i],Particles.vx[i],Particles.vy[i]);
        //                      for (int j=0; j < Particles.np; j++ ) {
        //                              if ( j != i ) {
        //                                      newparticle(&p2,Particles.weight[j],Particles.x[j],Particles.y[j],Particles.vx[j],Particles.vy[j]);
        //                                      ForceCompt(f,p1,p2);
        //                                      forces[index2D(0,i,2)] =
        //                                      forces[index2D(0,i,2)] + f[0];
        //                                      forces[index2D(1,i,2)] =
        //                                      forces[index2D(1,i,2)] + f[1];
        //                              };
        //                      };
        //               };
        //               ComptPopulation(&Particles, forces);
        fprintf(stdout, "Step %d of %d\n", k + 1, MaxSteps);
    };

    printf("Free device memory...\n");
    HANDLE_ERROR(cudaFree(posX_dev));
    HANDLE_ERROR(cudaFree(posY_dev));
    HANDLE_ERROR(cudaFree(velX_dev));
    HANDLE_ERROR(cudaFree(velY_dev));
    HANDLE_ERROR(cudaFree(weight_dev));
    HANDLE_ERROR(cudaFree(forces_dev));

    printf("Free host memory...\n");
    free(GenFieldGrid.Values);
    free(ParticleGrid.Values);
    free(Particles.x);
    free(Particles.y);
    free(Particles.vx);
    free(Particles.vy);
    free(Particles.weight);

    time(&t1);
    fprintf(stdout, "Ending   at: %s", asctime(localtime(&t1)));
    fprintf(stdout, "Computations ended in %lf seconds\n", difftime(t1, t0));

    fprintf(stdout, "End of program!\n");

    return (0);
}
