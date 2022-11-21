/*********************************************************
* This program solves coupled equations which describs
* the diffusion-damped system with constant coefficients.
* A time splitting technique by Marchuk–Strang is used.
* For two subproblems, A1 and A2 we solve A1 with time length
* tau=(1/2)*dt then solve the other one, A2, for time length
* dt and solve with A1 again for time length tau=(1/2)*dt.
*
* The read_input function should not be modified, but the
* use of the function in main() may or may not be correct.
*
* To compile: gcc -Wall -Werror -std=c99 -o assign3 assign3.c -lm
*
* List of identified errors:
*--------+--------------------------------------------
*  Line  |     Brief description of a fix
* Number |
*--------+-------------------------------------------
*  93    |     Added & infront of variables in read_input
*--------+-------------------------------------------
*  72    |     Added math.h library
*--------+-------------------------------------------
*  78    |     Change double to int
*--------+-------------------------------------------
*  175   |     Changed sfac to -sfac
*--------+-------------------------------------------
*114to121|     Changed malloc(nx) to malloc((nx+1)*sizeof(double));
*--------+-------------------------------------------
*  163   |     Changed uts1[k] = uc[k] + D * deriv * 0.5*dt; to uts1[k] = uc[k] + (( - vc[k] +  D * deriv ) * 0.5*dt);
*--------+-------------------------------------------
*  167   |     Changed vts1[k] = vc[k] + D * deriv * 0.5*dt; to vts1[k] = vc[k] + (( uc[k] + D * deriv ) * 0.5*dt);
*--------+-------------------------------------------
*  184   |     Changed un[k] = uts2[k] + D * deriv * 0.5*dt; to un[k] = uts2[k] + (( - vts2[k] +( D * deriv ))* 0.5*dt);
*--------+-------------------------------------------
*  188   |     Changed vn[k] = vts2[k] + D * deriv * 0.5*dt; to vn[k] = vts2[k] + (( uts2[k] +( D * deriv ))* 0.5*dt);
*--------+-------------------------------------------
*152, 154|     Changed 2*dt to dt
*--------+-------------------------------------------
*  124   |     Changed = to ==
*--------+-------------------------------------------
*  234   |     Changed free(&uc) etc to
* to 241 |     free(uc) etc switched return 0 to exit(0)
*--------+-------------------------------------------
*  231   |     Change ctime += dt; line from 227
*--------+-------------------------------------------
*        |     Changed copying time values to:
*  203   |         double tmp;
*   to   |         tmp = vn;
*        |         // vn = *vc;
*  219   |         vc = tmp;
*        |
*        |         tmp = vts2;
*        |         // vts2 = *vts1;
*        |         vts1 = tmp;
*        |
*        |         tmp = uts2;
*        |         //uts2 = *uts1;
*        |         uts1 = tmp;
*        |
*        |         tmp = un;
*        |         //un = *uc;
*        |         uc = tmp;
*--------+-------------------------------------------
*  221   |      uc[1]=uc[nx-4];
*  222   |      vc[1]=vc[nx-4];
*        |
*--------+-------------------------------------------
********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI  3.141592

void read_input(double *D, double *L, int *nx, double *t_F);

int main(void) {
  /******************************/
  /* Declarations of parameters */
  /******************************/
  /* Number of grid points */
  int nx;
  /* Length of domain */
  double L;
  /* Equation coefficients */
  double D;
  /* Length of time to run simulation. */
  double t_F;

  /* Read in from file; */
  read_input(&D, &L, &nx, &t_F);//checks out
  //printf("%lf, %lf, %d, %lf,\n ",D,L,nx,t_F);

  /* Grid spacing */
  double dx = L/(nx-1);//because L = (nx − 1) × dx
  //printf(" %lf, ",dx);
  double invdx2 = 1.0/(dx*dx);//checks out
  //printf(" %lf, ",invdx2);
  /* Time step */
  double dt = 0.5;

  /************************************************/
  /* Solution Storage at Current / Next time step */
  /************************************************/
  double *uc, *un, *vc, *vn;
  /* Time splitting solutions */
  double *uts1, *uts2, *vts1, *vts2;
  /* Derivative used in finite difference */
  double deriv;
  nx = nx + 4;

  /* Allocate memory according to size of nx*/
  uc = malloc((nx+1)*sizeof(double));
  un = malloc((nx+1)*sizeof(double));
  vc = malloc((nx+1)*sizeof(double));
  vn = malloc((nx+1)*sizeof(double));
  uts1 = malloc((nx+1)*sizeof(double));
  uts2 = malloc((nx+1)*sizeof(double));
  vts1 = malloc((nx+1)*sizeof(double));
  vts2 = malloc((nx+1)*sizeof(double));

  /* Check the allocation pointers */
  if (uc==NULL || un==NULL || vc==NULL || vn==NULL || uts1==NULL || uts2==NULL || vts1==NULL || vts2==NULL)
  {
    printf("Memory allocation failed\n");
    return 1;
  }
  //printf("Memory allocation success\n");

  int k;
  double x;
  /* Current time */
  double ctime = 0.0;

  /* Initialise arrays Explicit first step because lack of derivatives*/
  for(k = 0; k <= nx; k++) {
    x = (k-1)*dx;//checks out
    //printf(" %lf, ", x);
    uc[k]  = 1.0 + sin(2.0*PI*x/L);//checks out
    //printf("%lf\n",uc[k]);
    vc[k]  = 0.0;
    /* Set other arrays to 0 */
    uts1[k] = 0; uts2[k] = 0;
    vts1[k] = 0; vts2[k] = 0;
  }

  /* Loop over timesteps */
  while (ctime < t_F){

    /* Rotation factors for time-splitting scheme. */
    double cfac = cos(dt);//checks out
    //printf(" %lf, ", cfac);
    double sfac = sin(dt);//checks out
    //printf(" %lf, ", sfac);

    /* First substep for diffusion equation, A_1 */
    for (k = 1; k <= nx; k++) {
      x = (k-1)*dx;
      /* Diffusion at half time step. */
      deriv = (uc[k-1] + uc[k+1] - 2*uc[k])*invdx2;
      //printf(" %lf %lf\n ", x, deriv);
      uts1[k] = uc[k] + (( - vc[k] +  D * deriv ) * 0.5*dt);

      deriv = (vc[k-1] + vc[k+1] - 2*vc[k])*invdx2;
      //printf(" %lf %lf\n ", x, deriv);
      vts1[k] = vc[k] + (( uc[k] + D * deriv ) * 0.5*dt);
      //printf("%lf %lf %lf\n ",x, vts1[k], uts1[k]);
    }

    /* Second substep for decay/growth terms, A_2 */
    for (k = 1; k <= nx; k++) {
      x = (k)*dx;
      /* Apply rotation matrix to u and v, */
      uts2[k] = cfac*uts1[k] - sfac*vts1[k];
      vts2[k] = sfac*uts1[k] + cfac*vts1[k];
    }

    /* Thirds substep for diffusion terms, A_1 */
    for (k = 1; k <= nx-4; k++) {
      x = (k)*dx;
      deriv = (uts2[k-1] + uts2[k+1] - 2*uts2[k])*invdx2;
      //printf(" %lf %lf\n ", x, deriv);
      un[k] = uts2[k] + (( - vts2[k] +( D * deriv ))* 0.5*dt);
      //printf(" %lf, ", un[k]);
      deriv = (vts2[k-1] + vts2[k+1] - 2*vts2[k])*invdx2;
      //printf(" %lf %lf\n ", x, deriv);
      vn[k] = vts2[k] + (( uts2[k] +( D * deriv ))* 0.5*dt);
      //printf("%lf %lf %lf\n ",x,vn[k],un[k]);
    }

    /* Copy next values at timestep to u, v arrays. */
    /*double tmp;

    tmp = vn[k];
    vn[k] = vc[k];
    vc[k] = tmp;

    tmp = un[k];
    un[k] = uc[k];
    uc[k] = tmp;*/

    double *tmp;

    tmp = vn;
    //vn = vc;
    vc = tmp;

    tmp = vts2;
    //*vts2 = *vts1;
    vts1 = tmp;

    tmp = uts2;
    //*uts2 = *uts1;
    uts1 = tmp;

    tmp = un;
    //*un = *uc;
    uc = tmp;
    //periodic boundary conditions at x=0, x=L
    uc[1]=uc[nx-4];
    vc[1]=vc[nx-4];

    /* Increment time. */

    for (k = 1; k <= nx-4; k++ ) {
      x = (k-1)*dx;
      printf("%g %g %g %g\n",ctime,x,uc[k],vc[k]);
      //printf("%g\n",vc[k]);
    }
    ctime += dt;
  }
/* Free allocated memory */
  //free(uc);
  //free(un);
  //free(vc);
  //free(vn);
  //free(uts1);
  //free(uts2);
  //free(vts1);
  //free(vts2);

  //return 0;
  exit(0);//memory is freed through exit(0)
}

// The lines below don't contain any bugs! Don't modify them
void read_input(double *D, double *L, int *nx, double *t_F) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(4!=fscanf(infile,"%lf %lf %d %lf",D,L,nx,t_F)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}
