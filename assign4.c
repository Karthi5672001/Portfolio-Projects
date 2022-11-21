/*****************************************************************
 * Prior to compilation execute following lines on nenneke:
 * module purge
 * module load intel impi imkl
 * Then:
 * Compile:  gcc -o assign4 assign4.c -lm -lmkl -liomp5
 * Run: ./assign4
 * ****************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mkl_lapacke.h>

/* Define structure that holds band matrix information (from bandu_utility.c) */
struct band_mat
{
  long ncol;        /* Number of columns in band matrix */
  long nbrows;      /* Number of rows (bands in original matrix) */
  long nbands_up;   /* Number of bands above diagonal */
  long nbands_low;  /* Number of bands below diagonal */
  double *array;    /* Storage for the matrix in banded format */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information */
};
/* Define a new type band_mat (from bandu_utility.c) */
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters. (from bandu_utility.c) */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns)
{
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL)
  {
    return 0;
  }
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++)
  {
    bmat->array[i] = 0.0;
  }
  return 1;
};

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat)
{
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.    (from bandu_utility.c)       */
double *getp(band_mat *bmat, long row, long column)
{
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol )
  {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Return the value of a location in the band matrix, using
   the row and column indexes of the full matrix.    (from bandu_utility.c)       */
double getv(band_mat *bmat, long row, long column)
{
  return *getp(bmat,row,column);
}

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix. (from bandu_utility.c)*/
double setv(band_mat *bmat, long row, long column, double *val)
{
  *getp(bmat,row,column) = *val;
  return *val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays         (from bandu_utility.c)     */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b)
{
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++)
  {
    for (bandno=0;bandno<bmat->nbrows;bandno++)
    {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat)
{
  long i,j;
  for(i=0; i<bmat->ncol;i++)
  {
    for(j=0; j<bmat->nbrows; j++)
    {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}
/*read input (from bandu_utility.c)*/
void read_input(double *L, int *n, double *D, double *v, double *k_plus, double *k_minus)
{
   FILE *input;

    //allocate memory
   input=fopen("input.txt","r");
   fscanf(input,"%lf %d %lf %lf %lf %lf",L,n,D,v,k_plus,k_minus);
   fclose(input);
}

/*  read coefficients (from assign3.c)   */
void read_coefficients(int *n, double *s, double *delta)
{
   //allocate memory
   FILE *coeff_file;
   if(!(coeff_file=fopen("coefficients.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }


   long j;
   for(j=0;j<*n;j++) {
      fscanf(coeff_file,"%lf %lf", &s[j],&delta[j]);
   }
   fclose(coeff_file);
}

/* An example of how to use the band matrix routines to solve a PDE:
   The equation solved is related to the steady state solution of the heat
   diffusion equation.
*/
int main()
{
  band_mat bmat;
  //turn these into pointers
  double *D = malloc(sizeof(double));
  double *v = malloc(sizeof(double));
  double *k_plus = malloc(sizeof(double));
  double *k_minus = malloc(sizeof(double));
  double *L = malloc(sizeof(double));
  int *n = malloc(sizeof(int));
  /*  read coefficients */
  read_input(L, n, D, v, k_plus, k_minus);
  if(*n < 3)
  {
      printf("Error reading parameters from file\n");
      exit(1);
  }
  if(*v < 0)
  {
      *v = - *v;
  }
  //long = (long*)malloc(sizeof(long));
  long ncols = *n*2;
  double *x = malloc(sizeof(double)*ncols);
  double *b = malloc(sizeof(double)*ncols);

  //double *b = malloc(sizeof(double)*ncols);

  /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1
  */
  long nbands_low = *n;
  long nbands_up  = *n;
  init_band_mat(&bmat, nbands_low, nbands_up, ncols);
  //define dx
  double dx = *L/(*n);
  //allocate memory to pointers
  double *s = malloc(*n*sizeof(double));
  double *a = malloc(sizeof(double));
  double *e = malloc(sizeof(double));
  double *c = malloc(sizeof(double));
  double *d = malloc(sizeof(double));
  double *f = malloc(sizeof(double));
  double *h = malloc(sizeof(double));
  double *g = malloc(sizeof(double));
  double *k = malloc(sizeof(double));
  double *delta = malloc(*n*sizeof(double));

  read_coefficients(n, s, delta);
  //Terms in matrix
  *a = -(*D - *v * dx);
  *c = -(*D);
  *d = -(*D - *v * dx);
  *e = -(-2 * *D + *v * dx - dx * dx * *k_plus);
  *f = -(*D);
  *g = -(dx * dx * *k_plus);
  *k = -(dx * dx * *k_minus);
  //checking inputs
  //printf("%lf\n", dx);
  // printf("%lf %d %lf %lf %lf %lf\n",*L,*n,*D,*v,*k_plus ,*k_minus);
  printf("%lf %lf \n",*a,*c);
  printf("%lf %lf \n",*d,*f);
  printf("%lf %lf %lf \n",*g,*k,*e);
  long i;
  for(i=0; i<ncols; i++)
  {
    b[i] = 0.0;
  }
  for(i=0; i<*n; i++)
  {
    b[i] = dx*dx*s[i];
  }
  for(i=0; i<ncols; i++)
  {
    printf("%lf\n", b[i]);
  }
  //printmat(&bmat);
  //printf("\n---------------------------------\n");
  /* Loop over the equation number and set the matrix
     values equal to the coefficients of the grid values.
     Note that boundaries are treated with special cases */
  for(i=0; i<*n; i++)
  {
    if(i>0)    {setv(&bmat,i,i-1,c);};
                    setv(&bmat,i,i,e);
    if(i<*n-1) {setv(&bmat,i,i+1,a);};
  }
  //printf("\n---------------------------------\n");
  //printmat(&bmat);
  for(i = *n; i<ncols; i++)
  {
    *h = -(-2 * *D + *v * dx - dx*dx* delta[i-*n] - dx * dx * *k_minus);// - dx * dx * *k_plus
    //printf("%lf ",delta[i]);
    //printf("%lf \n",*h);
    if(i>*n)       {setv(&bmat,i,i-1,f);};
                    setv(&bmat,i,i,h);
    if(i<ncols-1) {setv(&bmat,i,i+1,d);};
  }
  //k_minus
  for(i=0; i<*n; i++)
  {

    {setv(&bmat,i+*n,i,g);};
  }
  //k_pluis
  for(i=0; i<*n; i++)
  {
    {setv(&bmat,i,i+*n,k);};
  }
  //Periodic boundary conditions
  setv(&bmat,0,*n-1,c);
  setv(&bmat,*n-1,0,a);
  setv(&bmat,*n,ncols-1,f);
  setv(&bmat,ncols-1,*n,d);

  /*  Print coefficient matrix for debugging: */
  printf("\n---------------------------------\n");
  //printmat(&bmat);

  solve_Ax_eq_b(&bmat, x, b);

  for(i=0; i<*n; i++)
  {
    //printf("%lf %lf %lf \n",i*dx,x[i],x[i+*n]);
  }
  //Store results in a file
  FILE *fp = fopen("output.txt", "w+");
  for(i=0; i<*n; i++)
  {
      fprintf( fp ,"%lf %lf %lf \n",i*dx,x[i],x[i+*n]);
  }
  fclose(fp);
/*free memory*/
finalise_band_mat(&bmat);
free(x);
free(s);
free(delta);
free(L);
free(n);
free(D);
free(v);
free(k_minus);
free(k_plus);
free(a);
free(b);
free(c);
free(d);
free(e);
free(f);
free(h);
free(g);
free(k);
return(0);
}

