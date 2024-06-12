#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "winbgi2.h"
#include "rk4.h"
#define g 9.81
#define l 1.5
#define m 1.0
#define pi 3.14
#define b 0.2
void fun(double t, double* X, double* F);
double enm(double k, double pk);
int main()
{
    FILE* f1 = fopen("kat.txt", "w"); 
    FILE* f2 = fopen("predkosck.txt", "w");
    FILE* f3 = fopen("energia.txt", "w");
    FILE* f4 = fopen("predkoscodk.txt", "w");
    double an = 10.0;
    double k0 = pi * an / 180.0;
    double pk0 = 1;
    double Y[2] = { k0,pk0 };
    double Y1[2];
    double t = 0;
    double tk = 40;
    double t0 = t;
    int N = 512;
    double dt = (tk - t0) / N;

    while (t < tk)
    {
        vrk4(t, Y, dt, 2, fun, Y1);
        double E = enm(Y[0], Y[1]);
        fprintf(f1, "%lf %lf \n", t, Y1[0]);
        fprintf(f2, "%lf %lf\n", t, Y1[1]);
        fprintf(f3, "%1.6lf %1.6lf \n", t, E);
        fprintf(f4, "%lf %lf\n", Y1[0], Y1[1]);
        Y[0] = Y1[0];
        Y[1] = Y1[1];
        t += dt;
    }
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    return 0;
}
double enm(double k, double pk)
{
    return ((m * pow(l, 2)) / 2) * pow(pk, 2) + m * g * l * (1 - cos(k)); //wzor na energie mechaniczna z rownan
}
void fun(double t, double* X, double* F)
{
    F[0] = X[1];
    F[1] = (-1) * (g / l) * sin(X[0])-(b/(m*l))*X[1]; //wzor z wyprowadzonych rownan
}