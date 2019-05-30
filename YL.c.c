#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    FILE* out;
    double h, t, a, x;
    double *u1, *u2, *uSol;
    int M, N, n, m;

    out = fopen("outLY.txt", "w");

    h = 0.001;
    t = 0.0001;
    a = 0.5;

    M = 2./h;
    N = 1./t;

    u1 = (double *)malloc(M*sizeof(double));
    u2 = (double *)malloc(M*sizeof(double));

    uSol = (double *)malloc(M*sizeof(double));

    x = -1.;
    //Theory
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
    x = -1.;
    //Initial conds
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

    for (n = 0; n < N; n++)
    {
        u1[0] = 1.; // Boundary cond
        u1[M-1] = 0.; // Boundary cond
        for (m = 1; m < M-1; m++)
        {
            u2[m] = t*( (a+1)*u1[m-1]/(2*h) - (a-1)*u1[m+1]/(2*h) + (1/t - 1/h)*u1[m] );
        }
        for (m = 1; m < M-1; m++)
        {
            u1[m] = u2[m];
            //printf(" %.2f ", u1[m]);
            //fprintf(out, "%.2f ", u1[m]);
        }
        //printf("\n");
        //fprintf(out, "\n");
    }

    fprintf(out, "x ");
    fprintf(out, "u ");
    fprintf(out, "delta\n");
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
    free(u1); free(u2); free(uSol);
    return 0;
}
