#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define M_PI 3.14159265358979323846
// Calculation's accuracy
#define EPS 1e-11
// Parameter of equation
#define ALPHA      1.
//define ALPHA      0.
//#define ALPHA      .01
//#define ALPHA      .1
//#define ALPHA      .5

#define J(i,j) J[3*(i)+(j)]
#define A(i,j) A[3*(i)+(j)]

double func(int i, double t, double *x);

void printMatr(double *A);
double maxVec(double *vec);

void K1(double t, double *x, double h, double *k1);
void K2(double t, double *x, double h, double *k1, double *k2);
void K3(double t, double *x, double h, double *k1, double *k2, double *k3);
void K4(double t, double *x, double h, double *k2, double *k3, double *k4);
void K5(double t, double *x, double h, double *k1, double *k2, double *k4, double *k5);
void K6(double t, double *x, double h, double *k1, double *k2, double *k3, double *k4, double *k5, double *k6);
void dX(double *k1, double *k3, double *k4, double *dx);
void R(double *k1, double *k3, double *k4, double *k5, double *k6, double *r);
double locErr(double *k1, double *k3, double *k4, double *k5, double *k6, double *r, double *err);
int RK6(double h, double t_s, double t_e, double* x_s,
        double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
        double *r, double *err, double *dx,
        double *x);

double discrFunc(double h, double t_s, double t_e, double* x_s,
                double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
                double *r, double *err, double *dx,
                double *x,
                int i, double *beta);
int createSoLAE(double h, double t_s, double t_e, double* x_s,
                double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
                double *r, double *err, double *dx,
                double *x,
                double *beta, double *J, double *rP);
double normFedorenko(double h, double t_s, double t_e, double* x_s,
                    double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
                    double *r, double *err, double *dx,
                    double *x,
                    double *beta, double *J, double *rP, double *kappa);
//double normFedorenko(double *J, double *rP, double *kappa);
void Gauss3x3(double *J, double *rP, double* dh);
double detMatr3x3(double *A);
void Kramer3x3(double *J, double *rP, double* dh, double *A);
// with Gauss
/*int Newton(double h, double t_s, double t_e, double* x_s,
            double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
            double *r, double *err, double *dx,
            double *x,
            double *beta, double *J, double *rP, double *kappa, double *dh, double *root, double *discr);*/
// with Kramer
int Newton(double h, double t_s, double t_e, double* x_s,
            double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
            double *r, double *err, double *dx,
            double *x,
            double *beta, double *J, double *rP, double *kappa, double *dh, double *root, double *A, double *discr);

// With Gauss MAIN
/*int main(void)
{
    int i;
    double start, end;
    double *k1, *k2, *k3, *k4, *k5, *k6;
    double *r, *err, *dx;
    double *x_s, *x;
    double *J, *rP, *kappa, *dh, *root;
    double *beta;
    double t_s = 0.0;
    //double t_e = 1.*M_PI;
    double t_e = 1.;
    double *discr;
    double step = 0.01;
    int cnt;

    printf("S T A R T\n");
    start = clock();

    k1 = (double *)malloc(6*sizeof(double));
    k2 = (double *)malloc(6*sizeof(double));
    k3 = (double *)malloc(6*sizeof(double));
    k4 = (double *)malloc(6*sizeof(double));
    k5 = (double *)malloc(6*sizeof(double));
    k6 = (double *)malloc(6*sizeof(double));
    r = (double *)malloc(6*sizeof(double));
    err = (double *)malloc(6*sizeof(double));
    dx = (double *)malloc(6*sizeof(double));
    x_s = (double *)malloc(6*sizeof(double));
    x = (double *)malloc(6*sizeof(double));
    beta = (double *)malloc(3*sizeof(double));
    J = (double *)malloc(9*sizeof(double));
    rP = (double *)malloc(3*sizeof(double));
    kappa = (double *)malloc(3*sizeof(double));
    dh = (double *)malloc(3*sizeof(double));
    root = (double *)malloc(3*sizeof(double));
    discr = (double *)malloc(1*sizeof(double));

    beta[0] = 2.17;
    beta[1] = 15.97;
    beta[2] = -91.67;

    // True root
    beta[0] = 3.01;
    beta[1] = 16.95;
    beta[2] = -93.84;

    x_s[0] = 0.;
    x_s[1] = 0.;
    x_s[2] = 0.;
    x_s[3] = beta[0];
    x_s[4] = beta[1];
    x_s[5] = beta[2];

    cnt =  Newton(step, t_s, t_e, x_s,
                    k1, k2, k3, k4, k5, k6,
                    r, err, dx,
                    x,
                    beta, J, rP, kappa, dh, root, discr);

    end = clock();
    printf("S U C C E S S !\n");
    for(i = 0; i < 3;  i++) printf("root[%d] = %f\n", i, root[i]);
    printf("\nNumber of Newton's steps: %d\n", cnt);
    printf("time = %.6f seconds\n", (end - start)/(CLOCKS_PER_SEC));

    free(k1); free(k2); free(k3); free(k4); free(k5); free(k6);
    free(r); free(err); free(dx);
    free(x_s); free(x);
    free(beta);
    free(J); free(rP); free(kappa); free(dh);
    free(root); free(discr);

    return 0;
}*/
// With Kramer MAIN
int main(void)
{
    FILE*out;
    int i, j;
    double start, end;
    double *k1, *k2, *k3, *k4, *k5, *k6;
    double *r, *err, *dx;
    double *x_s, *x;
    double *J, *rP, *kappa, *dh, *root;
    double *A;
    double *beta;
    double t_s = 0.0;
    double t_e = 1.*M_PI;
    double *discr;
    //double t_e = 1.;
    double step = 0.01;
    int cnt;

    printf("S T A R T\n");
    start = clock();

    k1 = (double *)malloc(6*sizeof(double));
    k2 = (double *)malloc(6*sizeof(double));
    k3 = (double *)malloc(6*sizeof(double));
    k4 = (double *)malloc(6*sizeof(double));
    k5 = (double *)malloc(6*sizeof(double));
    k6 = (double *)malloc(6*sizeof(double));
    r = (double *)malloc(6*sizeof(double));
    err = (double *)malloc(6*sizeof(double));
    dx = (double *)malloc(6*sizeof(double));
    x_s = (double *)malloc(6*sizeof(double));
    x = (double *)malloc(6*sizeof(double));
    beta = (double *)malloc(3*sizeof(double));
    J = (double *)malloc(9*sizeof(double));
    rP = (double *)malloc(3*sizeof(double));
    kappa = (double *)malloc(3*sizeof(double));
    dh = (double *)malloc(3*sizeof(double));
    root = (double *)malloc(3*sizeof(double));
    A = (double *)malloc(9*sizeof(double));
    discr = (double *)malloc(1*sizeof(double));

    out = fopen("roots.txt", "w");
    // Test
    beta[0] = 0.;
    beta[1] = 0.;
    beta[2] = 0.;

    // True root 22 ALPHA 0.1
    /*beta[0] = 3.01;
    beta[1] = 16.95;
    beta[2] = -93.84;*/

    // True root 21 ALPHA 0.0
    /*beta[0] = 0.92;
    beta[1] = 0.46;
    beta[2] = 0.58;*/

    /*for (i = 0; i < 15; i++)
    {
        beta[0] = 0. + i*1.;
        beta[1] = 0. + i*1.;
        beta[2] = 0. + i*1.;

        x_s[0] = 0.;
        x_s[1] = 0.;
        x_s[2] = 0.;
        x_s[3] = beta[0];
        x_s[4] = beta[1];
        x_s[5] = beta[2];

        cnt =  Newton(step, t_s, t_e, x_s,
                        k1, k2, k3, k4, k5, k6,
                        r, err, dx,
                        x,
                        beta, J, rP, kappa, dh, root, A, discr);

        if (cnt != -1)
        {
            printf("S U C C E S S !\n");
            //for(j = 0; j < 3;  j++) printf("root[%d] = %f\n", j, root[j]);
            //printf("\nNumber of Newton's steps: %d\n", cnt);
            for(j = 0; j < 3;  j++) fprintf(out, "root[%d] = %f\n", j, root[j]);
            fprintf(out, "Number of Newton's steps: %d\n", cnt);
            fprintf(out, "Discrepancy is %.15f\n \n", discr[0]);
        }
        else
        {
            printf("There is some problem with RK method!\nTry another beta!\n");
        }
    }*/
    x_s[0] = 0.;
    x_s[1] = 0.;
    x_s[2] = 0.;
    x_s[3] = beta[0];
    x_s[4] = beta[1];
    x_s[5] = beta[2];

    cnt =  Newton(step, t_s, t_e, x_s,
                    k1, k2, k3, k4, k5, k6,
                    r, err, dx,
                    x,
                    beta, J, rP, kappa, dh, root, A, discr);

    if (cnt != -1)
    {
        printf("S U C C E S S !\n");
        for(i = 0; i < 3;  i++) printf("root[%d] = %.11f\n", i, root[i]);
        printf("\nNumber of Newton's steps: %d\n", cnt);
        printf("Discrepancy is %.15f\n \n", discr[0]);
    }
    else
    {
        printf("There is some problem with RK method!\nTry another beta!\n");
    }

    end = clock();
    printf("time = %.6f seconds\n", (end - start)/(CLOCKS_PER_SEC));

    fclose(out);

    free(k1); free(k2); free(k3); free(k4); free(k5); free(k6);
    free(r); free(err); free(dx);
    free(x_s); free(x);
    free(beta);
    free(J); free(rP); free(kappa); free(dh);
    free(root);
    free(A); free(discr);

    return 0;
}
// Verification of Runge-Kutta's method MAIN
/*int main(void)
{
    int i;
    double start, end;
    double *k1, *k2, *k3, *k4, *k5, *k6;
    double *r, *err, *dx;
    double *x_s, *x;
    double t_s = 0.0;
    double t_e = 1.*M_PI;
    //double t_e = 1.;
    double step = 0.01;
    int cnt;

    start = clock();

    k1 = (double *)malloc(6*sizeof(double));
    k2 = (double *)malloc(6*sizeof(double));
    k3 = (double *)malloc(6*sizeof(double));
    k4 = (double *)malloc(6*sizeof(double));
    k5 = (double *)malloc(6*sizeof(double));
    k6 = (double *)malloc(6*sizeof(double));
    r = (double *)malloc(6*sizeof(double));
    err = (double *)malloc(6*sizeof(double));
    dx = (double *)malloc(6*sizeof(double));
    x_s = (double *)malloc(6*sizeof(double));
    x = (double *)malloc(6*sizeof(double));

    x_s[0] = 0.;
    x_s[1] = 1.;
    x_s[2] = 0.;
    x_s[3] = 0.;
    x_s[4] = 0.;
    x_s[5] = 0.;

    printf("S T A R T\n");
    cnt = RK6(step, t_s, t_e, x_s,
              k1, k2, k3, k4, k5, k6, r, err, dx,
              x);

    end = clock();
    printf("S U C C E S S !\n");
    for (i=0; i < 6; i++) printf("x[%d] = %f\n", i, x[i]);
    printf("\n Number of RK steps: %d\n", cnt);
    printf("time = %.6f seconds\n", (end - start)/(CLOCKS_PER_SEC));

    free(k1); free(k2); free(k3); free(k4); free(k5); free(k6);
    free(r); free(err); free(dx);
    free(x_s); free(x);
    return 0;
}*/
// Functions of differential equation
double func(int i, double t, double *x)
{
    switch (i)
    {
        case 0: // f1
            //return x[1];
            // 21
            return x[3]*exp(ALPHA*x[0]);
            // 22
            //return x[3]*exp(ALPHA*t);
            break;
        case 1: // f2
            //return -x[0];
            // 21
            return x[0]*sin(t);
            // 22
            //return x[0];
            break;
        case 2: // f3
            //return 0.;
            // 21
            return (x[0]*cos(t))/(1. + ALPHA*pow(t, 2.));
            // 22
            //return x[0]*x[0]*x[0]*t/(1. + ALPHA*t*t);
            break;
        case 3: // f4
            //return 0.;
            // 21
            return -x[4]*sin(t) - (x[5]*cos(t))/(1. + ALPHA*pow(t, 2.)) - (1./2.)*ALPHA*pow(x[3], 2.)*exp(ALPHA*x[0]);
            // 22
            //return -x[4] -x[5]*3.*x[0]*x[0]*t/(1. + ALPHA*t*t);
            break;
        case 4: // f5
            return 0.;
            break;
        case 5: // f6
            return 0. + 0.*(t + x[0] + x[1] + x[2] + x[3] + x[4] + x[5]);
            break;
        default:
            return -1.;
    }
}
// Matrix printing
void printMatr(double *J)
{
    int i, j;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            printf ("%.3f ", J(i,j));
        }
        printf("\n");
    }
}
// Maximum of vector's components
double maxVec(double *vec)
{
    int i;
    int n = 6;

    double max = vec[0];
    for (i = 0; i < n; i++)
    {
        if (vec[i] > max)
        {
            max = vec[i];
        }
    }
    return max;
}
// Formulas for Runge-Kutta of 6 order
void K1(double t, double *x, double h, double *k1)
{
    int i;

    for (i = 0; i < 6; i++) k1[i] = h*func(i, t, x);
}
void K2(double t, double *x, double h, double *k1, double *k2)
{
    int i;
    double tau;

    tau = t + h*0.5;
    for (i = 0; i < 6; i++) x[i] = x[i] + k1[i]*(1./2.);
    for (i = 0; i < 6; i++) k2[i] = h*func(i, tau, x);
    for (i = 0; i < 6; i++) x[i] = x[i] - k1[i]*(1./2.);
}
void K3(double t, double *x, double h, double *k1, double *k2, double *k3)
{
    int i;
    double tau;

    tau = t + h*0.5;
    for (i = 0; i < 6; i++) x[i] = x[i] + (k1[i] + k2[i])*(1./4.);
    for (i = 0; i < 6; i++) k3[i] = h*func(i, tau, x);
    for (i = 0; i < 6; i++) x[i] = x[i] - (k1[i] + k2[i])*(1./4.);
}
void K4(double t, double *x, double h, double *k2, double *k3, double *k4)
{
    int i;
    double tau;

    tau = t + h;
    for (i = 0; i < 6; i++) x[i] = x[i] + ((-1.)*k2[i] + 2.*k3[i]);
    for (i = 0; i < 6; i++) k4[i] = h*func(i, tau, x);
    for (i = 0; i < 6; i++) x[i] = x[i] - ((-1.)*k2[i] + 2.*k3[i]);
}
void K5(double t, double *x, double h, double *k1, double *k2, double *k4, double *k5)
{
    int i;
    double tau;

    tau = t + 2.*h/3.;
    for (i = 0; i < 6; i++) x[i] = x[i] + (7.*k1[i] + 10.*k2[i] + k4[i])/27.;
    for (i = 0; i < 6; i++) k5[i] = h*func(i, tau, x);
    for (i = 0; i < 6; i++) x[i] = x[i] - (7.*k1[i] + 10.*k2[i] + k4[i])/27.;
}
void K6(double t, double *x, double h, double *k1, double *k2, double *k3, double *k4, double *k5, double *k6)
{
    int i;
    double tau;

    tau = t + h/5.;
    for (i = 0; i < 6; i++) x[i] = x[i] + (28.*k1[i] - 125.*k2[i] + 546.*k3[i] + 54*k4[i] - 378.*k5[i])/625.;
    for (i = 0; i < 6; i++) k6[i] = h*func(i, tau, x);
    for (i = 0; i < 6; i++) x[i] = x[i] - (28.*k1[i] - 125.*k2[i] + 546.*k3[i] + 54*k4[i] - 378.*k5[i])/625.;
}
void dX(double *k1, double *k3, double *k4, double *dx)
{
    int i;

    for (i = 0; i < 6; i++) dx[i] = (k1[i] + 4*k3[i] + k4[i])/6.;
}
void R(double *k1, double *k3, double *k4, double *k5, double *k6, double *r)
{
    int i;

    for (i = 0; i < 6; i++) r[i] = -1.*(42.*k1[i] + 224.*k3[i] + 21.*k4[i] - 162.*k5[i] - 125.*k6[i])/336.;
}
double locErr(double *k1, double *k3, double *k4, double *k5, double *k6, double *r, double *err)
{
    int i;
    double loc_err;

    R(k1, k3, k4, k5, k6, r);
    for(i = 0; i < 6; i++) err[i] = fabs(r[i])/(pow(2., 6));
    loc_err = maxVec(err);

    return loc_err;
}
int RK6(double h, double t_s, double t_e, double* x_s,
        double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
        double *r, double *err, double *dx,
        double *x)
{
    int i;
    double t, sum_t, phi;
    int cnt = 0;

    t = t_s;
    for (i = 0; i < 6; i++) x[i] = x_s[i];
    sum_t = 0.;

    while (sum_t < (t_e - t_s))
    {
        if ((sum_t + h) > (t_e - t_s))
        {
            h = h/9.;
        }
        K1(t, x, h, k1);
        K2(t, x, h, k1, k2);
        K3(t, x, h, k1, k2, k3);
        K4(t, x, h, k2, k3, k4);
        K5(t, x, h, k1, k2, k4, k5);
        K6(t, x, h, k1, k2, k3, k4, k5, k6);

        R(k1, k3, k4, k5, k6, r);
        dX(k1, k3, k4, dx);
        phi = locErr(k1, k3, k4, k5, k6, r, err);

        //
        //printf("phi = %.15f\n", phi);
        //
        if (phi < EPS)
        {
            for (i = 0; i < 6; i++) x[i] = x[i] + dx[i] + r[i];
            t = t + h;
            sum_t = sum_t + h;
            cnt = cnt + 1;

            if (phi < EPS*0.01)
            {
                h = h*2.;
                //printf ("h enlarged to %f\n", h);
            }
        }
        if (phi > EPS)
        {
            if (h > EPS)
            {
                h = h/2;
                //printf ("h reduced to %f\n", h);
            }
            else
            {
                // FATAL ERROR !!!  FATAL ERROR !!!
                // FATAL ERROR !!!  FATAL ERROR !!!
                return -1;
                // FATAL ERROR !!!  FATAL ERROR !!!
                // FATAL ERROR !!!  FATAL ERROR !!!
            }
            //h = h/2;
            //printf ("h reduced to %f\n", h);
        }
        //printf("%.15f\n", t);
    }

    //printf("\n Number of RK steps: %d\n", cnt);
    //printf("\n Values for t = %f:\n", t);
    //for (i=0; i < 6; i++) printf("x[%d] = %f\n", i, x[i]);

    return cnt;
}
// Discrepancy equation functions
double discrFunc(double h, double t_s, double t_e, double* x_s,
                double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
                double *r, double *err, double *dx,
                double *x,
                int i, double *beta)
{
    FILE*out;
    int cnt;

    out = fopen("RK_steps.txt", "w");

    x_s[3] = beta[0];
    x_s[4] = beta[1];
    x_s[5] = beta[2];

    // Calculation of function's values
    // by Runge-Kutta's method of 6th order
    cnt = RK6(h, t_s, t_e, x_s,
              k1, k2, k3, k4, k5, k6, r, err, dx,
              x);

    // FATAL ERROR !!!  FATAL ERROR !!!
    // FATAL ERROR !!!  FATAL ERROR !!!
    if (cnt == -1)
    {
        return 123.;
    }
    // FATAL ERROR !!!  FATAL ERROR !!!
    // FATAL ERROR !!!  FATAL ERROR !!!

    fprintf(out, "%d\n", cnt);
    fclose(out);

    // Functions of discrepancy
    switch (i)
    {
        case 0:
            // 21
            return x[1] - 1.;
            // 22
            //return x[0] - 1.;
            break;
        case 1:
            // 21
            return x[2];
            // 22
            //return x[1];
            break;
        case 2:
            // 21
            return x[3];
            // 22
            //return x[2];
            break;
        default:
            return -1.;
    }
}
// Right part for System of Linear Algebraic Equations in point
// and
// Jacobi matrix of discrepancy function in point
// for System of Linear Algebraic Equations
int createSoLAE(double h, double t_s, double t_e, double* x_s,
                double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
                double *r, double *err, double *dx,
                double *x,
                double *beta, double *J, double *rP)
{
    int i, j, ver;
    double f, feps;

    for (i = 0; i < 3; i++)
    {
        f = discrFunc(h, t_s, t_e, x_s,
                        k1, k2, k3, k4, k5, k6,
                        r, err, dx,
                        x,
                        i, beta);

        // FATAL ERROR !!!  FATAL ERROR !!!
        // FATAL ERROR !!!  FATAL ERROR !!!
        if (fabs(f - 123.) < EPS)
        {
            //printf("There is some problem with RK method\n");
            ver = -1;
            return ver;
        }
        // FATAL ERROR !!!  FATAL ERROR !!!
        // FATAL ERROR !!!  FATAL ERROR !!!

        // Right part of SoLAE
        rP[i] = -f;
        // Jacobi matrix for SoLAE
        for (j = 0; j < 3; j++)
        {
            beta[j] = beta[j] + EPS;

            feps = discrFunc(h, t_s, t_e, x_s,
                            k1, k2, k3, k4, k5, k6,
                            r, err, dx,
                            x,
                            i, beta);
            J[3*i+j] = (feps - f)/EPS;

            beta[j] = beta[j] - EPS;
        }
    }
    ver = 100;
    return ver;
}
// Fedorenko norm
double normFedorenko(double h, double t_s, double t_e, double* x_s,
                    double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
                    double *r, double *err, double *dx,
                    double *x,
                    double *beta, double *J, double *rP, double *kappa)
//double normFedorenko(double *J, double *rP, double *kappa)
{
    int i, j, ver;
    double norm = 0.;

    for (i = 0; i < 3; i++)
    {
        kappa[i] = EPS;
    }

    ver = createSoLAE(h, t_s, t_e, x_s,
                            k1, k2, k3, k4, k5, k6,
                            r, err, dx,
                            x,
                            beta, J, rP);

    // FATAL ERROR !!!  FATAL ERROR !!!
    // FATAL ERROR !!!  FATAL ERROR !!!
    if (ver == -1)
        return -1.;
    // FATAL ERROR !!!  FATAL ERROR !!!
    // FATAL ERROR !!!  FATAL ERROR !!!

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            kappa[i] = kappa[i] + J(i,j)*J(i,j);
        }
        kappa[i] = sqrt(kappa[i]);
        norm = norm + rP[i]*rP[i]/kappa[i];
    }
    norm = sqrt(norm);

    return norm;
}
// Gauss's method for SoLAE
void Gauss3x3(double *J, double *rP, double* dh)
{
    double max, tmp;
    int k, index, i, j;

    k = 0;
    while (k < 3)
    {
        max = fabs(J(k,k));
        index = k;

        for (i = k + 1; i < 3; i++)
        {
            if (fabs(J[(i)*3]) > max)
            {
                max = J[(i)*3];
                index = i;
            }
        }

        /*if (fabs(max) < EPS*1e-3)
        {
            printf("Bad SoLAE\n");
        }*/

        for (i = 0; i < 3; i++)
        {
            tmp = J(k,i);
            J(k,i) = J(index,i);
            J(index,i) = tmp;
        }

        for (i = 0; i < 3; i++)
        {
            tmp = rP[k];
            rP[k] = rP[index];
            rP[index] = tmp;
        }

        for (i = k; i < 3; i++)
        {
            tmp = J(i, k);
            if (fabs(J(i,k))<EPS) continue;

            for (j = 0; j < 3; j++)
            {
                J(i,j) = J(i,j)/tmp;
            }
            rP[i] = rP[i]/tmp;

            if (i == k) continue;
            for (j = 0; j < 3; j++)
            {
                J(i,j) = J(i,j) - J(k,j);
            }
            rP[i] = rP[i] - rP[k];
        }
        k++;
    }

    for (k = 2; k >= 0; k--)
    {
        dh[k] = rP[k];
        for (i = 0;  i < k; i++)
        {
            rP[i] = rP[i] - J(i,k)*dh[k];
        }
    }
}
// Determinant of matrix
double detMatr3x3(double *A)
{
    double a1, a2, a3;
    a1 = A[4]*A[8] - A[7]*A[5];
    a2 = A[3]*A[8] - A[5]*A[6];
    a3 = A[3]*A[7] - A[4]*A[6];

    return A[0]*a1 - A[1]*a2 + A[2]*a3;
}
// Kramer's method for SoLAE
void Kramer3x3(double *J, double *rP, double* dh, double *A)
{
    int i, j;
    double delta, delta1, delta2, delta3;

    for(i = 0; i < 9; i++) A[i] = J[i];
    delta = detMatr3x3(A);
    if(fabs(delta) < EPS) printf("WARNING!!!\n det = %f\n", delta);

    for(i = 0; i < 9; i++) A[i] = J[i];
    for(j = 0; j < 3; j++) A(j,0) = rP[j];
    delta1 = detMatr3x3(A);

    for(i = 0; i < 9; i++) A[i] = J[i];
    for(j = 0; j < 3; j++) A(j,1) = rP[j];
    delta2 = detMatr3x3(A);

    for(i = 0; i < 9; i++) A[i] = J[i];
    for(j = 0; j < 3; j++) A(j,2) = rP[j];
    delta3 = detMatr3x3(A);

    dh[0] = delta1/delta;
    dh[1] = delta2/delta;
    dh[2] = delta3/delta;

    //printf("det = %.15f\n", delta);
    //printf("det1 = %.15f\n", delta1);
    //printf("det2 = %.15f\n", delta2);
    //printf("det3 = %.15f\n", delta3);
}
// Newton's method for discrepancy equation
// with Gauss
/*int Newton(double h, double t_s, double t_e, double* x_s,
            double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
            double *r, double *err, double *dx,
            double *x,
            double *beta, double *J, double *rP, double *kappa, double *dh, double *root, double *discr)*/
// with Kramer
int Newton(double h, double t_s, double t_e, double* x_s,
            double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,
            double *r, double *err, double *dx,
            double *x,
            double *beta, double *J, double *rP, double *kappa, double *dh, double *root, double *A, double *discr)
{
    int i;
    int cnt = 0;
    double p = 2.;
    double norm1, norm2;

    for (i = 0; i < 3; i++)
    {
        root[i] = beta[i];
    }

    norm1 = normFedorenko(h, t_s, t_e, x_s,
                            k1, k2, k3, k4, k5, k6,
                            r, err, dx,
                            x,
                            root, J, rP, kappa);
    norm2 = norm1;

    // FATAL ERROR !!!  FATAL ERROR !!!
    // FATAL ERROR !!!  FATAL ERROR !!!
    if (fabs(norm2 + 1.) < EPS)
    {
        //printf("There is some problem!\n");
        return -1;
    }
    // FATAL ERROR !!!  FATAL ERROR !!!
    // FATAL ERROR !!!  FATAL ERROR !!!

    while(fabs(norm2) > EPS)
    {
        //
        //printMatr(J);
        //
        // SoLAE was created in normFedorenko
        // Gauss START
        /*Gauss3x3(J, rP, dh);
        //
        for (i = 0; i < 3; i++) printf("%f\n", dh[i]);
        printf("Gauss OKEY\n");*/
        // Gauss END

        // Kramer START
        Kramer3x3(J, rP, dh, A);
        //
        //for (i = 0; i < 3; i++) printf("dh[%d] = %.15f\n", i, dh[i]);
        //printf("Kramer OKEY\n");
        // Kramer END

        // Isaev-Sonin modification START
        for(i = 0; i < 3; i++)
        {
            root[i] = beta[i] + dh[i];
        }
        norm2 = normFedorenko(h, t_s, t_e, x_s,
                                k1, k2, k3, k4, k5, k6,
                                r, err, dx,
                                x,
                                root, J, rP, kappa);
        while ((norm2 - norm1) >= EPS)
        {
            for(i = 0; i < 3; i++)
            {
                root[i] = beta[i] + dh[i]/pow(2, p);
            }
            norm2 = normFedorenko(h, t_s, t_e, x_s,
                                    k1, k2, k3, k4, k5, k6,
                                    r, err, dx,
                                    x,
                                    root, J, rP, kappa);

            // NEW SoLAE was created in normFedorenko
            p = p + 1.;
            if (p > 10) continue;
        }
        p = 2.;
        // Isaev-Sonin modification END

        for (i = 0; i < 3; i++)
        {
            beta[i] = root[i];
        }

        cnt = cnt + 1;
        if(cnt > 100)
        {
            printf("Bad start point!\n");
            return cnt;
        }
        //
        discr[0] = norm2;
        //printf("Number of step is %d\nDiscrepancy is %.15f\n", cnt, norm2);
        //printf("\n");

        //
        // FATAL ERROR !!!  FATAL ERROR !!!
        // FATAL ERROR !!!  FATAL ERROR !!!
        if (fabs(norm2 + 1.) < EPS)
        {
            //printf("There is some problem!\n");
            return -1;
        }
        // FATAL ERROR !!!  FATAL ERROR !!!
        // FATAL ERROR !!!  FATAL ERROR !!!
    }
    return cnt;
}
