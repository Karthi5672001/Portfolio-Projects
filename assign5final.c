/*****************************************************************
 * Prior to compilation execute following lines on nenneke:
 * module purge
 * module load intel impi imkl
 * Then:
 * Compile:  gcc -o assign5 assign5.c -lm -lmkl -liomp5
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
void read_input(int *nx,int *ny,double *Lx,double *Ly,double *tD,double *tf,int *s_x0,int *s_x1,int *s_y0,int *s_y1)
{
   FILE *input;

   input=fopen("input.txt","r");
   fscanf(input,"%d",nx);
   fscanf(input,"%d",ny);
   fscanf(input,"%lf",Lx);
   fscanf(input,"%lf",Ly);
   fscanf(input,"%lf",tf);
   fscanf(input,"%lf",tD);
   fscanf(input,"%d",s_x0);
   fscanf(input,"%d",s_x1);
   fscanf(input,"%d",s_y0);
   fscanf(input,"%d",s_y1);
   fclose(input);
}

/*  read coefficients (from assign3.c)   */
void read_coefficients(int *nx, int *ny, double *x_coefficient)
{
   //allocate memory
   FILE *coeff_file;
   if(!(coeff_file=fopen("coefficients.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }


   long j;
   for(j=0;j<*nx * *ny;j++) {
      fscanf(coeff_file,"%lf", &x_coefficient[j]);
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
  int *nx = malloc(sizeof(int));
  int *ny = malloc(sizeof(int));
  double *Lx = malloc(sizeof(double));
  double *Ly = malloc(sizeof(double));
  double *tD = malloc(sizeof(double));
  double *tf = malloc(sizeof(double));
  //boundary conditions
  int *s_x0 = malloc(sizeof(int));
  int *s_x1 = malloc(sizeof(int));
  int *s_y0 = malloc(sizeof(int));
  int *s_y1 = malloc(sizeof(int));
  if (nx==NULL || ny==NULL || Lx==NULL || Ly==NULL || tD==NULL || tf==NULL || s_x0==NULL || s_x1==NULL || s_y0==NULL || s_y1==NULL)
  {
    printf("Memory allocation failed\n");
    return 1;
  }
  /*  read coefficients */
  read_input(nx, ny, Lx, Ly, tD, tf, s_x0, s_x1, s_y0, s_y1);
  if(*nx < 2)
  {
      printf("nx is too small\n");
      exit(1);
  }
  if(*ny < 2)
  {
      printf("ny is too small\n");
      exit(1);
  }
  int ncols = *nx * *ny;
  //define dx,dy
  double dx = *Lx/(*nx);
  double dy = *Ly/(*ny);
  //No. timestep
  double Nt = *tf/(*tD);
  //printf("%lf\n",Nt);
  double *x = malloc(sizeof(double)*ncols);
  double *y = malloc(sizeof(double)*ncols);
  if (x==NULL || y==NULL)
  {
    printf("Memory allocation failed\n");
    return 1;
  }
  /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1
  */
  long nbands_low = *nx*(*ny-1);
  long nbands_up  = *nx*(*ny-1);
  init_band_mat(&bmat, nbands_low, nbands_up, ncols);
  //allocate memory to pointers
  double *x_coefficient = malloc(ncols*sizeof(double));
  if (x_coefficient==NULL)
  {
    printf("Memory allocation failed\n");
    return 1;
  }
  //double *y_coefficient = malloc(*ny*sizeof(double));
  read_coefficients(nx, ny, x_coefficient/*, y_coefficient*/);
  //Coefficients
  double dxsq = dx*dx;
  double dysq = dy*dy;
  double dxdt = 0.25*(dx/(*tD));
  double dydt = 0.25*(dy/(*tD));
  double *a = malloc(sizeof(double));
  double *b = malloc(sizeof(double));
  double *c = malloc(sizeof(double));
  double *d = malloc(sizeof(double));
  double *e = malloc(sizeof(double));
  double *f = malloc(sizeof(double));
  double *c_d = malloc(sizeof(double));
  double *a_n = malloc(sizeof(double));
  double *b_n = malloc(sizeof(double));
  double *d_n = malloc(sizeof(double));
  double *e_n = malloc(sizeof(double));
  double *a_n1 = malloc(sizeof(double));
  double *b_n1 = malloc(sizeof(double));
  double *d_n1 = malloc(sizeof(double));
  double *e_n1 = malloc(sizeof(double));

  if (a==NULL || b==NULL || c==NULL || d==NULL || e==NULL || f==NULL || c_d==NULL || a_n==NULL || b_n==NULL || d_n==NULL || e_n==NULL || a_n1==NULL || b_n1==NULL || d_n1==NULL || e_n1==NULL)
  {
    printf("Memory allocation failed\n");
    return 1;
  }
  //Terms in matrix
  *a = (1.0/(dxsq)+(1/dx)*dxdt);
  *b = (1.0/(dysq)+(1/dy)*dydt);
  *c = -((2.0/(dxsq))+(2.0/(dysq)));
  *d = (1.0/(dxsq)-(1/dx)*dxdt);
  *e = (1.0/(dysq)-(1/dy)*dydt);
  *f = 0.0;
  //dirichlect
  *c_d = +1-((2.0/(dxsq))+(2.0/(dysq)));
  //neumann
  *a_n = (2.0/(dxsq)+(1/dx)*dxdt);
  *b_n = (2.0/(dysq)+(1/dy)*dydt);
  *d_n = (2.0/(dxsq)-(1/dx)*dxdt);
  *e_n = (2.0/(dysq)-(1/dy)*dydt);
  *a_n1 = ((1/dx)*dxdt);
  *b_n1 = ((1/dy)*dydt);
  *d_n1 = (-(1/dx)*dxdt);
  *e_n1 = (-(1/dy)*dydt);
  printf("%lf %lf %lf %lf %lf\n",*a,*b,*c,*d,*e);
  long i,j,k;
  //Store results in a file
  FILE *fp = fopen("output.txt", "w+");
  //u values for t=0
  for(j=0; j<*ny; j++)
  {for(i=0; i<*nx; i++){
      fprintf( fp ,"0.000000 %lf %lf %lf\n",i*dx,j*dy,x_coefficient[i+*nx*j]);
  }}
  //printmat(&bmat);
  //printf("\n---------------------------------\n");
  //printmat(&bmat)
  for(k=1; k<Nt; k++){
  for(i=0; i<(*ny * *nx); i++) {
      //for(j=0; j<*nx * *ny; j++) {
  if(i/(*nx) != (int)(i/(*nx))) {if(i>0){
        setv(&bmat,i-1,i,d);
  }}
  if(i/(*nx) != (int)(i/(*nx))) {if(i<(*nx * *ny)-1){
        setv(&bmat,i+1,i,a);
  }}
  if(i<((*ny * *nx)-*nx)) {
        setv(&bmat,i,i+*nx,e);
  }
  if(i<((*ny * *nx)-*nx)) {
        setv(&bmat,i+*nx,i,b);
  }
  setv(&bmat,i,i,c);
  /* Left hand side of eqaution */
  y[i] = (-1.0+1/(1+(x_coefficient[i]*x_coefficient[i])));
  }
  /*if(*s_y0 == 0){
      for(i=0; i<*nx; i++){
            y[i] = (-1.0+1/(1+0));
      }}
  if(*s_x0 == 0){
      for(i=0; i<*nx; i++){
            y[i**nx] = (-1.0+1/(1+0));
      }}
  if(*s_y1 == 0){
      for(i=0; i<*nx; i++){
            y[i+(*nx*(*ny-1))] = (-1.0+1/(1+0));
      }}
  if(*s_x1 == 0){
      for(i=0; i<*nx; i++){
            y[*nx*(i+1)-1] = (-1.0+1/(1+0));
      }}*/
  //Matrix Corrections
  for(i=0; i<*ny; i++){
  if(*nx*(i+1)-1<*nx* *ny) {
        setv(&bmat,*nx*(i+1)-1,*nx*i,d);
  }
  if(*nx*(i+1)-1<*nx* *ny) {
        setv(&bmat,*nx*i,*nx*(i+1)-1,a);
  }
  }
  for(i=0; i<*nx; i++){
        setv(&bmat,*nx*(*ny-1)+i,i,e);
     }
  for(i=0; i<*nx; i++){
        setv(&bmat,i,*nx*(*ny-1)+i,b);
     }
  //Neumann boundary for y0
  if(*s_y0 == 1){
        for(i=0; i<*nx; i++){
  if(i>0){setv(&bmat,i-1,i,d_n);}
  if(i<(*nx * *ny)-1){setv(&bmat,i+1,i,a_n1);}
  if(i+(*nx*(*ny-1))<(*nx* *ny)){setv(&bmat,i+(*nx*(*ny-1)),i,e_n);}
  if(i<(*nx)){setv(&bmat,i+*nx,i,b_n1);}}}
  //Neumann boundary for x0
  if(*s_x0 == 1){
        for(i=0; i<*ny; i++){
  if(i>0){if(*nx*i<(*nx* *ny)){setv(&bmat,*nx*i-1,*nx*i,d_n);}}
  if(*nx*i+1<(*nx * *ny)){setv(&bmat,*nx*i+1,*nx*i,a_n1);}
  if(*nx*i+*nx<(*ny * *nx)){setv(&bmat,*nx*i,*nx*i+*nx,e_n);}
  if(*nx*(i+1)<(*nx * *ny)){setv(&bmat,*nx*(i+1),*nx*i,b_n1);}}}
  //Neumann boundary for y1
  if(*s_y1 == 1){
        for(i=0; i<*nx; i++){
  if(i>0){if(i+*nx*(*ny-1)+1<(*nx* *ny)){setv(&bmat,i+*nx*(*ny-1),i+*nx*(*ny-1)+1,d_n1);}}
  if(i+*nx*(*ny-1)+1<(*ny * *nx)){setv(&bmat,i+*nx*(*ny-1)+1,i+*nx*(*ny-1),a_n);}
  if(*nx*i+*nx<(*ny * *nx)){setv(&bmat,*nx*(*ny-2)+i,*nx*(*ny-1)+i,e_n1);}
  if(i+*nx*(*ny-1)<(*ny * *nx)){setv(&bmat,i,i+*nx*(*ny-1),b_n);}}}

  //Neumann boundary for x1
  if(*s_x1 == 1){
        for(i=0; i<*ny; i++){
  if(*nx*(i+2)-1<(*ny * *nx)){setv(&bmat,*nx*(i+2)-1,*nx*(i+1)-1,b_n);}//y
  if(*nx*(i+1)-1<(*ny * *nx)){setv(&bmat,*nx*i,*nx*(i+1)-1,a_n);}//x
  if((*nx*(i+2)-1)<(*ny * *nx)){setv(&bmat,*nx*(i+1)-1,*nx*(i+2)-1,e_n1);}//y
  if(*nx*(i+1)-1<(*ny * *nx)){setv(&bmat,*nx*(i+1)-2,*nx*(i+1)-1,d_n1);}}//x
  }

  //dirichlect boundary conditions for x0
  if(*s_y0 == 0){
      for(i=0; i<*nx; i++){
      //if(i>0){setv(&bmat,i,i-1,f);}//A
      //if(i+1<*nx){setv(&bmat,i,i+1,f);}//D
      //if(i+*nx<ncols-1){setv(&bmat,i,i+*nx,f);}//E
      //if(i+*nx*(*ny-1)<ncols-1){setv(&bmat,i,i+*nx*(*ny-1),f);}//B
      setv(&bmat,i,i,c_d);}}//C
      //setv(&bmat,0,*nx-1,f);//A
      //setv(&bmat,*nx-1,0,f);}//D}
  //dirichlect boundary conditions for y0
  if(*s_x0 == 0){
      for(i=0; i<*ny; i++){
      //if(*nx*i+1<ncols-1){setv(&bmat,*nx*i,*nx*i+1,f);}//D
      //if(*nx*(i+1)<ncols-1){setv(&bmat,*nx*(i+1),*nx*i,f);}//B
      //if(*nx*(i+1)-1<ncols-1){setv(&bmat,*nx*i,*nx*(i+1)-1,f);}//A
      //if(*nx*(i+1)<ncols-1){setv(&bmat,*nx*i,*nx*(i+1),f);}//E
      setv(&bmat,*nx*i,*nx*i,c_d);}}//C
      //setv(&bmat,*nx*(*ny-1),0,f);//E
      //setv(&bmat,0,*nx*(*ny-1),f);}//B
  //dirichlect boundary conditions for x1
  if(*s_y1 == 0){
      for(i=0; i<*nx; i++){
      //if(i+*nx*(*ny-1)+1<ncols-1){setv(&bmat,i+*nx*(*ny-1),i+*nx*(*ny-1)+1,f);}//A
      //if(i+*nx*(*ny-1)+1<ncols-1){setv(&bmat,i+*nx*(*ny-1)+1,i+*nx*(*ny-1),f);}//D
      //if(i+*nx*(*ny-1)<ncols-1){setv(&bmat,i+*nx*(*ny-1),i,f);}//E
      //if(i+*nx*(*ny-1)<ncols-1){setv(&bmat,i+*nx*(*ny-1),i+*nx*(*ny-2),f);}//B
      setv(&bmat,i+*nx*(*ny-1),i+*nx*(*ny-1),c_d);}}//C
      //setv(&bmat,*nx*(*ny-1),(*nx* *ny)-1,f);//A
      //setv(&bmat,(*nx* *ny)-1,*nx*(*ny-1),f);}//D
  //dirichlect boundary conditions for y1
  if(*s_x1 == 0){
      for(i=0; i<*ny; i++){
      //if(*nx*(i+1)-1<ncols-1){setv(&bmat,*nx*(i+1)-1,*nx*(i+1)-2,f);}//A
      //if(*nx*(i+1)-1<ncols-1){setv(&bmat,*nx*(i+1),*nx*i,f);}//B
      //if(*nx*(i+1)-1<ncols-1){setv(&bmat,*nx*(i+1)-1,*nx*i,f);}//D
      //if(*nx*(i+2)-1<ncols-1){setv(&bmat,*nx*(i+1)-1,*nx*(i+2)-1,f);}//E
      setv(&bmat,*nx*(i+1)-1,*nx*(i+1)-1,c_d);}}//C
      //setv(&bmat,(*nx* *ny)-1,*nx-1,f);//E
      //setv(&bmat,*nx-1,(*nx* *ny)-1,f);}//B
  /*  Print coefficient matrix for debugging: */
  //printf("---------------------------------\n");
  //printmat(&bmat);
  for(i=0; i<ncols; i++)
  {
      //printf("%lf\n", y[i]);
  }
  solve_Ax_eq_b(&bmat, x, y);
  /*for(j=0; j<*ny; j++)
  {for(i=0; i<*nx; i++)
    printf("%lf %lf %lf %lf\n",*tD*k,i*dx,j*dy,x[i+*nx*j]);
  }*/
  for(j=0; j<*ny; j++)
  {for(i=0; i<*nx; i++){
      fprintf( fp ,"%lf %lf %lf %lf\n",*tD*k,i*dx,j*dy,x[i+*nx*j]);
  }}
  for(i=0; i<ncols; i++){
        x_coefficient[i] = x[i];
  }
  }
  fclose(fp);
/*free memory*/
finalise_band_mat(&bmat);
//free(x);
free(x);free(a);
free(b);free(c);free(d);free(e);free(f);free(a_n);free(b_n);free(c_d);free(d_n);free(e_n);free(a_n1);free(b_n1);free(e_n1);free(d_n1);
free(x_coefficient);free(y);free(nx);free(ny);
free(Lx);free(Ly);free(tD);free(tf);
free(s_x0);free(s_x1);free(s_y0);free(s_y1);
return(0);
}
