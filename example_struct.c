/******************************************************************************

                            Online C Compiler.
                Code, Compile, Run and Debug C program online.
Write your code in this editor and press "Run" button to compile and execute it.

*******************************************************************************/

#include <stdio.h>

struct i2dGrid
{
   int EX, EY;            // extensions in X and Y directions
   double Xs, Xe, Ys, Ye; // initial and final value for X and Y directions
   int *Values;           // 2D matrix of values
} GenFieldGrid, ParticleGrid, GenFieldGridDev, ParticleGridDev; // le ultime due sono le strutture sul device


void f(struct i2dGrid *p) {
    unsigned long s = sizeof(*p);
    
    
    printf("%lu\n", s);
    printf("%d\n", (*p).EX);
}

int main()
{
    printf("%lu\n", sizeof(GenFieldGrid));
    printf("%lu\n", 2*(sizeof(GenFieldGrid.EX)));
    GenFieldGrid.EX = 0;
    printf("%lu\n", 4*(sizeof(GenFieldGrid.Xs)));
    printf("%lu %p\n", sizeof(GenFieldGrid.Values), GenFieldGrid.Values);
    printf("\n\n");
    
    
    f(&GenFieldGrid);
    

    return 0;
}



