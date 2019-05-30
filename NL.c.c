#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void progonka (int M, double * u1, double * u2,
               double C, double D,
               double *rP, double *y, double *XX,
               double *l, double *u, double *v);
void LU (int M, double A, double C, double D,
         double *l, double *u, double *v);

int main(void)
{
    FILE *out;
    double h, t, a, x;
    double *u1, *u2, *uSol;
    int M, N, n, m, i;
    double A, C, D;
    double *l, *v, *u;
    double *y, *rP, *XX;

    h = 0.01;
    t = 0.001;
    a = 0.5;

    D = -((a*t)/(2.*h) + (t*t*a*a)/(2.*h*h));
    A = 1. + (t*t*a*a)/(h*h);
    C = (a*t)/(2.*h) - (t*t*a*a)/(2.*h*h);

    M = 2./h;
    N = 1./t;

    //LU-decomposition
    l = (double *)malloc((M-2)*sizeof(double));
    v = (double *)malloc((M-2)*sizeof(double));
    u = (double *)malloc((M-3)*sizeof(double));

    //progonka
    y = (double *)malloc((M-2)*sizeof(double));
    rP = (double *)malloc((M-2)*sizeof(double));
    XX = (double *)malloc((M-2)*sizeof(double));

    u1 = (double *)malloc(M*sizeof(double));
    u2 = (double *)malloc(M*sizeof(double));
    uSol = (double *)malloc(M*sizeof(double));

    out = fopen("outNL.txt","w");

    //Theory
    x = -1.;
    for (m = 0; m < M; m++)
    {
        x = x + h;
        if (x <= 0.25)
        {
            uSol[m] = 1.;
        }
        if ((x > 0.25)&&(x <= 0.5))
        {
            uSol[m] = -4.*(x-0.5);
        }
        if (x > 0.5)
        {
            uSol[m] = 0.;
        }
    }
    //initial conds
    x = -1.;
    for (m = 0; m < M; m++)
    {
        x = x + h;
        if (x <= -0.25)
        {
            u1[m] = 1.;
        }
        if ((x > -0.25)&&(x <= 0))
        {
            u1[m] = -4.*x;
        }
        if (x > 0)
        {
            u1[m] = 0.;
        }
    }

    LU(M, A, C, D, l, u, v);

    for (n = 0; n < N; n++)
    {
        printf("t = %d\n", n);
        //Boundary conds
        u1[0] = 1.;
        u1[M-1] = 0.;
        u2[0] = 1.;
        u2[M-1] = 0.;

        progonka(M, u1, u2, C, D, rP, y, XX, l, u, v);

        for (i = 1; i < M-1; i++)
        {
            u2[i] = XX[i-1];
        }

        for (m = 0; m < M; m++)
        {
            u1[m] = u2[m];
        }
    }

    x = -1.;
    for(m = 0; m < M; m++)
    {
        x = x + h;
        fprintf(out, "%.4f ", x);
        fprintf(out, "%.4f ", u1[m]);
        fprintf(out, "%.4f ", uSol[m]);
        fprintf(out, "%.4f\n", u1[m]-uSol[m]);
    }

    fclose(out);
    free(u1); free(u2); free(l); free(u); free(v);
    free(y); free(rP); free(XX);
    return 0;
}

//LU decomposition
void LU (int M, double A, double C, double D,
         double *l, double *u, double *v)
{
    for (int i = 0; i < M-3; i++)
    {
        v[i] = D;
    }
    v[M-3] = 1.;

    l[0] = A/v[0];
    u[0] = C/l[0];

    for (int i = 1; i < M - 3; i++)
    {
        l[i] = (A - u[i-1])/v[i];
        u[i] = C/l[i];
    }
    l[M-3] = A - u[M-4];
}

void progonka (int M, double * u1, double * u2,
               double C, double D,
               double *rP, double *y, double *XX,
               double *l, double *u, double *v)
{
    //right part
    rP[0] = -D*u2[0] + u1[1];
    rP[M-3] = -C*u2[M-1] + u1[M-2];
    for (int i=1; i<M-3; i++)
    {
        rP[i] = u1[i+1];
    }

    y[0] = rP[0]/l[0];
    for (int i = 1; i < M-2; i++)
    {
        y[i] = (rP[i] - y[i-1])/l[i];
    }

    XX[M-3] = y[M-3]/v[M-3];
    for (int i = M - 4; i > -1; i--)
    {
        XX[i] = (y[i] - u[i]*XX[i+1])/v[i];
    }
}
