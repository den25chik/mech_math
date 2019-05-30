#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#define X_MIN -1.
#define X_MAX 1.
#define T_MIN 0
#define T_MAX 0.125
#define eps_n 1e-3

double u0(double x);
double* f(double T, int m, double H, double *u_prev,double *u_curr);
double* progonka(double* f,double* u_prev,double* u_curr, int m, double H, double T);
double* newton(double *u_prev, int m, double H, double T);

double u0(double x)
{
    if (x <= -0.25)
    {
        return 1.;
    }
    if ((x > -0.25)&&(x <= 0))
    {
        return -4.*x;
    }
    if (x > 0)
    {
        return 0.;
    }
}

double* f(double T, int m, double H, double *u_prev,double *u_curr)
{
	int i;
	double* f;
	f = (double*)malloc((m-2)*sizeof(double));
	for(i = 1; i < m-1; i++)
	{
		f[i-1] = u_curr[i] - u_prev[i] + (u_curr[i+1]*u_curr[i+1] - u_curr[i-1]*u_curr[i-1])*T/(4*H) - (u_prev[i]*(u_curr[i+1]*u_curr[i+1]-u_curr[i]*u_curr[i]) - u_prev[i-1]*(u_curr[i]*u_curr[i]-u_curr[i-1]*u_curr[i-1]))*(T*T/(4*H*H));
	}

	return f;
}

double* progonka(double* f, double* u_prev, double* u_curr, int m, double H, double T)
{
	double *x;
	double* a;
	double* b;
	double* c;
	int i;
	x = (double*)malloc((m-2)*sizeof(double));
	a = (double*)malloc((m-2)*sizeof(double));
	b = (double*)malloc((m-2)*sizeof(double));
	c = (double*)malloc((m-2)*sizeof(double));


	for(i = 1;i < m-1;i++)
	{
		b[i-1] = 1 - (T/(2*H*H))*(u_prev[i] - u_prev[i-1])*u_curr[i];

	}

	for(i = 1;i < m-2;i++)
	{
		c[i-1] = ( (T/(2*H)) - (T*T/(2*H*H))*u_prev[i] )*u_curr[i+1];
		a[i-1] = -( (T/(2*H)) + (T*T/(2*H*H))*u_prev[i-1] )*u_curr[i-1];
	}


	for(i = 1;i < m-2;i++)
	{
		b[i] = b[i] - (a[i-1]/b[i-1])*c[i-1];
		f[i] = f[i] - (a[i-1]/b[i-1])*f[i-1];
	}


	x[m-3] = f[m-3]/b[m-3];

	for(i = m-4;i >= 0;i--)
	{
		x[i]=(f[i]-c[i]*x[i+1])/b[i];
	}

	return x;
}


double* newton(double *u_prev, int m, double H, double T)
{
	int k = 0;
	int i = 0;
	double *delta = (double*)malloc((m-2)*sizeof(double));
	double *func = (double*)malloc((m-2)*sizeof(double));
	double *u_curr = (double*)malloc(m*sizeof(double));
	double *u_curr_new = (double*)malloc((m-2)*sizeof(double));


	double max_mod_razn = 0;
	for(k = 0; k < m-2; k++) u_curr_new[k] = u_prev[k+1];
	for(k = 0; k < m; k++) u_curr[k] = u_prev[k];
	do
	{
		func = f(T, m, H, u_prev, u_curr);
		for(k = 0; k < m-2; k++) func[k] = -func[k];
		delta = progonka(func, u_prev, u_curr, m, H, T);
		for(k = 0; k < m-2; k++) u_curr_new[k] = u_curr_new[k] + delta[k];

		max_mod_razn = 0;
		for(k = 0; k < m-2; k++)
		{
			if(fabs(delta[k]) > max_mod_razn)
			{
				max_mod_razn = fabs(delta[k]);
			}
		}

		u_curr[0] = 1;
		u_curr[m-1] = 0;
		for(k = 1; k < m-1; k++) u_curr[k] = u_curr_new[k-1];

		if(i > 1e+5)
		{
			printf("error in newton\n");
			break;
		}
		i++;
	}
	while(max_mod_razn > eps_n);

	return u_curr;
}



int main(void)
{
  int i, j, k;
  double T = 0.01;
  double H = 0.01;
  double pogr1, pogr2, pogr3,pogr4;
  double max_mod=0, max_mod_razn=0, sum_mod=0, sum_mod_razn=0;
  double *u_curr;//Для решения по схеме
  double *u_next;
  double *w_curr;//Сетка пополам
  double *w_next;
  double *temp;
  double *uSol; //Для аналитического решения
  int n; //По времени
  int m; //По координате
  int m_old;
  int d = 1;


  FILE *OUT;
  OUT = fopen("outNN.txt", "w");

  n=(int)((T_MAX-T_MIN)/T) + 1;
  m=(int)((X_MAX-X_MIN)/H) + 1;

  u_curr=(double *)malloc(m * sizeof(double));
  u_next=(double *)malloc(m * sizeof(double));
  uSol=(double *)malloc(m * sizeof(double));


  //Заполняем первый слой (t=0)
	for (i = 0; i < m; i++)
	{
		u_curr[i]= u0(X_MIN+H*i);
	}

  //Заполняем каждый следующий слой по предыдущему
    for (j = 1; j < n; j++)
	{
		u_next = newton(u_curr, m, H, T);

		temp=u_curr;
		u_curr=u_next;
		u_next=temp;
	}


	for (i = 0; i < m; i++) //Аналитическое решение
	{
		if ( X_MIN+i*H <= -0.125 ) uSol[i] = 1;
		if ( (X_MIN+i*H > -0.125) && (X_MIN+i*H <= 0 ) ) uSol[i] = -8*(X_MIN+i*H);
		if ( X_MIN+i*H > 0 ) uSol[i] = 0;
	}

	for (i = 0; i < m; i++) fprintf(OUT, "%f %f %f %f\n", X_MIN+i*H, u_curr[i], uSol[i], uSol[i] - u_curr[i]);

	//Считаем погрешности расчетов
	for (i = 0; i < m; i++)
	{
		if( fabs(u_curr[i]) > max_mod ) max_mod = fabs(u_curr[i]);
		if( fabs(u_curr[i] - uSol[i]) > max_mod_razn ) max_mod_razn = fabs(u_curr[i] - uSol[i]);
		sum_mod = sum_mod + fabs(u_curr[i]);
		sum_mod_razn = sum_mod_razn + fabs(u_curr[i] - uSol[i]);
	}
	pogr1 = max_mod_razn;
	pogr2 = H*sum_mod_razn;
	pogr3 = max_mod_razn / max_mod;
	pogr4 = sum_mod_razn / sum_mod;

	printf("%.5f & %.3f & %.3f & %.3f & %.3f & %.3f\n", T, H, pogr1, pogr2, pogr3, pogr4);
	m_old = m;

	//Сетка пополам
	for (k = 1; k < 5; k++)
	{
		T = T*0.5;
		H = H*0.5;
		d = d*2;
		n=(T_MAX-T_MIN)/T + 1;
		m=(X_MAX-X_MIN)/H + 1;
		w_curr=(double *)malloc(m * sizeof(double));
		w_next=(double *)malloc(m * sizeof(double));

		for (i = 0; i < m; i++)
		{
			w_curr[i]= u0(X_MIN+H*i);
		}


		for (j = 1; j < n; j++)
		{
			w_next = newton(w_curr, m, H, T);

			temp=w_curr;
			w_curr=w_next;
			w_next=temp;
		}


			max_mod_razn=0;
			sum_mod_razn=0;

		for (i = 0; i < m_old; i++)
		{
			if( fabs(u_curr[i] - w_curr[i*d]) > max_mod_razn ) max_mod_razn = fabs(u_curr[i] - w_curr[i*d]);
			sum_mod_razn = sum_mod_razn + fabs(u_curr[i] - w_curr[i*d]);
		}
		pogr1 = max_mod_razn;
		pogr2 = H*sum_mod_razn;
		pogr3 = max_mod_razn / max_mod;
		pogr4 = sum_mod_razn / sum_mod;

		printf("$v^%d$ & %.3f & %.3f & %.3f & %.3f\n", k, pogr1, pogr2, pogr3, pogr4);

		free(w_curr);
		free(w_next);
	}



	fclose(OUT);

  free(u_curr);
  free(u_next);
  free(uSol);
  return 0;
}
