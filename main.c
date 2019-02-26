#include <cerf.h>
#include <complex.h>
#include <cubature.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef _Complex double cdouble;

struct cubature_data
{
  const double q2;
  const double gamma;
  const int d[3];
  double ghat[3];
  double ghatinv[3];
  const double rtol;
} typedef cubature_data;

void usage() {
  printf("Usage: ./luscher_zeta.x <q^2> <gamma> <dx> <dy> <dz> <tolerance>\n");
}

void check_d(const int d[3])
{
  int sum_d = 0;

  for(int i=0; i<3; i++){ sum_d += (d[i] != 0); }

  if(sum_d > 1){
    printf("Error: boost must be along a single direction.\n");
    exit(-1);
  }
}

void set_ghat(cubature_data* data)
{
  if(data->d[0] != 0) {
    data->ghat[0]    = data->gamma;
    data->ghatinv[0] = 1.0 / data->gamma;
  } else if(data->d[1] != 0) {
    data->ghat[1]    = data->gamma;
    data->ghatinv[1] = 1.0 / data->gamma;
  } else if(data->d[2] != 0) {
    data->ghat[2]    = data->gamma;
    data->ghatinv[2] = 1.0 / data->gamma;
  }
}

double Z00_1(const cubature_data* data)
{
  double Z00 = 0.0;

  int n = 0;
  double Z00_0, r2mq2;
  double r = 1.0e+10;
  while(r > data->rtol)
  {
    Z00_0 = Z00;

    for(int nx=-n; nx<=n; nx++){
    for(int ny=-n; ny<=n; ny++){
    for(int nz=-n; nz<=n; nz++){

      if((nx > -n) && (nx < n) && (ny > -n) && (ny < n) && (nz > -n) && (nz < n)){ continue; }

      r2mq2 = pow(data->ghatinv[0]*(nx-0.5*data->d[0]), 2)
          + pow(data->ghatinv[1]*(ny-0.5*data->d[1]), 2)
          + pow(data->ghatinv[2]*(nz-0.5*data->d[2]), 2)
          - data->q2;
      Z00 += exp(-r2mq2) / r2mq2;

    }}}

    r = fabs(Z00-Z00_0) / fabs(0.5*(Z00+Z00_0));
    n++;
  }

  return Z00 / sqrt(4.0*M_PI);
}

double Z00_2(const cubature_data* data)
{
  const cdouble q = csqrt(data->q2);
  return -M_PI * data->gamma * ( exp(data->q2) - sqrt(M_PI) * creal(q*cerfi(q)) );
}

// double Z00_3_integrand(const double* t, const double* q2, const int d[3], const double ghat[3], const double* rtol)
int Z00_3_integrand(unsigned pdim, const double* s, void* params, unsigned fdim, double* f)
{
  cubature_data* data = (cubature_data*) params;

  double II = 0.0;

  int n = 1;
  double II0;
  double r = 1.0e+10;
  while(r > data->rtol)
  {
    II0 = II;

    for(int nx=-n; nx<=n; nx++){
    for(int ny=-n; ny<=n; ny++){
    for(int nz=-n; nz<=n; nz++){

      if((nx > -n) && (nx < n) && (ny > -n) && (ny < n) && (nz > -n) && (nz < n)){ continue; }

      II += pow(-1.0, nx*data->d[0] + ny*data->d[1] + nz*data->d[2]) *
          exp( -M_PI*M_PI/s[0] * ( data->ghat[0]*data->ghat[0]*nx*nx
          + data->ghat[1]*data->ghat[1]*ny*ny
          + data->ghat[2]*data->ghat[2]*nz*nz ) );

    }}}

    r = fabs(II-II0) / fabs(0.5*(II+II0));
    n++;
  }

  *f = pow(M_PI/s[0], 1.5) * exp(s[0]*data->q2) * II;

  return 0;
}

double Z00_3(cubature_data* data)
{
  double Z00, err;
  double s_min = 0.0;
  double s_max = 1.0;

  hcubature(1, Z00_3_integrand, data, 1, &s_min, &s_max, 0, 0.0, data->rtol, ERROR_INDIVIDUAL, &Z00, &err);

  return data->gamma * Z00 / sqrt(4.0*M_PI);
}

int main(int argc, char** argv)
{
  if(argc != 7){
    usage();
    exit(-1);
  }

  double q2    = atof(argv[1]);
  double gamma = atof(argv[2]);
  int d[3]     = { atoi(argv[3]), atoi(argv[4]), atoi(argv[5]) };
  double rtol  = atof(argv[6]);

  cubature_data data = { atof(argv[1]), atof(argv[2]), { atoi(argv[3]), atoi(argv[4]), atoi(argv[5]) },
                         { 1.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }, atof(argv[6]) };

  check_d(d);

  set_ghat(&data);

  double Z00 = Z00_1(&data) + Z00_2(&data) + Z00_3(&data);
  printf("%1.15e\n", Z00);

  return 1;
}