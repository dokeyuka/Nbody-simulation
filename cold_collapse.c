
/******************************************************************************
 template.c: template file for practice 1
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "g5util.h"

#include "cpgplot.h"

#define NMAX 16384
#define SQR(x) ((x) * (x))

#define WSIZE 6.0 /* size of window in inch */
#define VMAX 1.25 /* limit of variables in window */
#define XPLOT 0   /* coordinate for x */
#define YPLOT 1   /* coordinate for y */
#define TX 0.9    /* x coordinate of time */
#define TY -1.1   /* y coordinate of time */
#define CFG 5     /* foreground color */
#define CBG 0     /* background color */

/* the prototype declaration */
void make_spherical_df(int n, double m[], double x[][3], double v[][3], double r_v, double eps2);
void calc_force(int n, double m[], double x[][3], double a[][3], double eps2);
void leap_frog(int n, double m[], double x[][3], double v[][3], double a[][3], double dt, double eps2);
double calc_energy(int n, double m[], double x[][3], double v[][3], double eps2);
void calc_force_on_gpu(int n, double m[], double x[][3], double a[][3], double pot[], double eps2);

double drand48(void);

/* definition of function */
/*for graphics snapshot*/
void open_window(FILE **gp)
{

  if ((*gp = popen("gnuplot -persist", "w")) == NULL)
  {
    printf("gnuplot open error!!\n");
    exit(EXIT_FAILURE);
  }
  fprintf(*gp, "set size square\n");
  // fprintf(*gp, "set terminal gif animate\n");
  // fprintf(*gp, "set output \"plot.gif\"\n");
  fprintf(*gp, "set xrange [-%f:%f]\n", VMAX, VMAX);
  fprintf(*gp, "set yrange [-%f:%f]\n", VMAX, VMAX);

}

void close_window(FILE **gp)
{

  pclose(*gp);
}

void animated_snapshot(int n, double t, double x[][3], FILE **gp)
{

  int i;
  fprintf(*gp, "set key title \"%f\"\n", t);
  fprintf(*gp, "plot '-' with points pointtype 0 notitle\n");
  for (i = 0; i <= n; i++)
  {
    fprintf(*gp, "%f\t%f\n", x[i][0], x[i][1]);
  }
  fprintf(*gp, "e\n");
}

/**********************************************************************
Gaussian with mean = 0.0 and dispersion = 0.0 by Box-Muller method
******************************************************************/
double gaussian(void)
{
  double x, y, r2;
  double z;

  do
  {
    x = 2.0 * drand48() - 1.0;
    y = 2.0 * drand48() - 1.0;
    r2 = x * x + y * y;
  } while (r2 >= 1.0 || r2 == 0.0);
  z = sqrt(-2.0 * log(r2) / r2) * x; /* discard another Gaussian */

  return (z);
}


/**********************************************************************
 make_spherical_df :
  function to make particle distribution
  arg value    : n, m, c, v, r_v, eps2
  return vlaue :
***********************************************************************/
void make_spherical_df(int n, double m[], double x[][3], double v[][3], double r_v, double eps2)
{

  double K = 0;
  double W = 0;

  int i = 0;
  /*x,y,z : x[i][0], x[i][1], x[i][2] (-1 ~ 1 rand) */
  /*vx,vy,vz : v[i][0], v[i][1], v[i][2]*/
  while (i < n)
  {
    for (int j = 0; j < 3; j++)
    {
      x[i][j] = drand48() * 2.0 - 1.0;
    }
    double r = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2]);
    if (r > 1.0)
    {
      continue;
    }
    else
    {
      i++;
    }
  }

  /*W : potential energy*/
  for (int i = 0; i < n - 1; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      double r2 = SQR(x[j][0] - x[i][0]) + SQR(x[j][1] - x[i][1]) + SQR(x[j][2] - x[i][2]);
      W = W - m[i] * m[j] / sqrt(r2 + eps2);
    }
  }

  double sigma = sqrt(2.0 * r_v * fabs(W) / 3.0);
  fprintf(stderr, "sigma = %lf\n", sigma);

  /* velosity follows gaussian(mean:0, vari:sigma) */
  for (int i = 0; i < n; i++)
  {
    for (int k = 0; k < 3; k++)
    {
      v[i][k] = sigma * gaussian();
    }

  }
}

/**********************************************************************
calc_force :
  function to calculate initial acceleration
  arg value    : n, m, x, a, eps2
  return vlaue :
***********************************************************************/
void calc_force(int n, double m[], double x[][3], double a[][3], double eps2){

  /*initialization of accelration*/
  for (int i = 0; i < n; i++){
    for (int k = 0; k < 3; k++){
      a[i][k] = 0.0;
    }
  }

  for(int i = 0; i < n - 1; i++){
    for (int j = i + 1; j < n; j++){
      double r[3];
      r[0] = x[j][0] - x[i][0];  /* x */
      r[1] = x[j][1] - x[i][1];  /* y */
      r[2] = x[j][2] - x[i][2];  /* z */
      double abs_r = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + eps2);
      double r3inv = 1.0 / (abs_r*abs_r*abs_r);

      for (int k = 0; k < 3; k++)
      {
        a[i][k] += m[j] * r[k] * r3inv;
        a[j][k] += - m[j] * r[k] * r3inv;
      }
    }
  }
}

/**********************************************************************
leap_frog :
  function to calculate leap_frog
  arg value    : n, m, x, a, dt,eps2
  return value :
***********************************************************************/
void leap_frog(int n, double m[], double x[][3], double v[][3], double a[][3], double dt, double eps2)
{
  double v_half[n][3];

  for (int i = 0; i < n; i++){
    for (int k = 0; k < 3; k++){
      v_half[i][k] = v[i][k] + 0.5 * a[i][k] * dt;
      x[i][k] += v_half[i][k] * dt;
    }
  }

  calc_force(n, m, x, a, eps2);

  for (int i = 0; i < n; i++){
    for (int k = 0; k < 3; k++){
      v[i][k] = v_half[i][k] + 0.5 * a[i][k] * dt;
    }
  }
}

double calc_kinetic_energy(int n, double m[], double v[][3]){
  /*K : kinetic energy*/
  double K = 0.0;
  for (int i = 0; i < n; i++){
    K += m[i] * (SQR(v[i][0]) + SQR(v[i][1]) + SQR(v[i][2])) / 2.0;
  }
  return K;
}

double calc_potential_energy(int n, double m[], double x[][3], double eps2){
  double W = 0.0;
  /*W : potential energy*/
  for (int i = 0; i < n - 1; i++){
    for (int j = i + 1; j < n; j++){
      double r2 = SQR(x[j][0] - x[i][0]) + SQR(x[j][1] - x[i][1]) + SQR(x[j][2] - x[i][2]);
      W = W - m[i] * m[j] / sqrt(r2 + eps2);
    }
  }

  return W;
}
    /**********************************************************************
    calc_energy :
      function to calculate kinematic and potential energy
      arg value    : n, m, x, v, eps2
      return value :
    ***********************************************************************/
double calc_energy(int n, double m[], double x[][3], double v[][3], double eps2){
  double K;
  double W;
  double E;

  K = calc_kinetic_energy(n, m, v);
  W = calc_potential_energy(n, m, x, eps2);

  E = K + W;
  return E;
}

/**********************************************************************
calc_r_v :
  function to calculate r_v
  arg value    : n, m, x, v, eps2
  return value :
***********************************************************************/
double calc_r_v(int n, double m[], double x[][3], double v[][3], double eps2){
  double K, W, rv;

  K = calc_kinetic_energy(n, m, v);
  W = calc_potential_energy(n, m, x, eps2);

  rv = K/fabs(W);
  return rv;
}



void calc_force_on_gpu(int n, double m[], double x[][3], double a[][3], double pot[], double eps2){

  g5_open();
  g5_st_range(-256.0, 256.0, m[0]);

  g5_set_jp(0, n, m, x);
  g5_set_eps2_to_all(eps2);
  g5_set_n(n);

  g5_calculate_force_on_x(x, a, pot, n);
  g5_close();
}

int main(void)
{
  /**********************************************************************
    definition of variables
  ***********************************************************************/

  /* number of particles */
  int n;
  /* particle data */
  static double m[NMAX], x[NMAX][3], v[NMAX][3], a[NMAX][3];
  /* system energy, Virtual ratio */
  double e, e_ini, r_v, rv;
  /* current time, tmestep, end time, data interval */
  double t, dt, t_end, t_out;
  int k;
  /* softening parameter */
  double eps;
  /* sqared softening parameter */
  double eps2;

  /**********************************************************************
     set necessary variables
  ***********************************************************************/

  fprintf(stderr, "enter values of (eps k(index of dt) t_end t_out N r_v). eps must be about 1/100*phisycal system. dt must be lower than eps/v & 2^k\n");
  scanf("%lf %d %lf %lf %d %lf", &eps, &k, &t_end, &t_out, &n, &r_v);



  /*visualization for check */
  /*eps = 0.03125;
  t_end = 1.0;
  t_out = 0.25;
  n = 1024;
  r_v = 0.1;
  fprintf(stderr, "enter k\n");
  scanf("%d", &k);
  */

  dt = pow(2.0, k);

  for (int i2 = 0; i2 < n; i2++)
  {
    m[i2] = 1.0 / n;
  }
  eps2 = SQR(eps);


  /*for check*/
  fprintf(stderr, "eps2 = %g\n dt = %g\n t_end = %g\n t_out = %g\n n = %d\n r_v= %g\n m = %g\n", eps2, dt, t_end, t_out, n, r_v, m[0]);

  /* initialization of particle data */
  make_spherical_df(n, m, x, v, r_v, eps2);

  /*for check*/
  /*
  for (int l = 0; l < n; l++){
    printf("%f %f %f\n", x[l][0], x[l][1], x[l][2]);
  }
  */

  /* caluculate initial energy */
  e_ini = calc_energy(n, m, x, v, eps2);
  /* calculate initial acceralation */
  g5_open();
  calc_force_on_gpu(n, m, x, a, eps2);
  g5_close();

  FILE *gp;
  open_window(&gp);

  while(t < t_end){
    /* realtime analysis */
    if (fmod(t, t_out) == 0.0){
      rv = calc_r_v(n, m, x, v, eps2);
      printf("%lf %lf\n", t, rv);
    }

    /* time integration */
    leap_frog(n, m, x, v, a, dt, eps2);
    t += dt;

    /* visualization of particle data */
    animated_snapshot(n, t, x, &gp);

    }
  close_window(&gp);

  /* integration error check
  e = calc_energy(n, m, x, v, r_v, eps2);
  printf("%e %e\n", dt, fabs((e - e_ini) / e_ini));
  */
  return 0;
}
