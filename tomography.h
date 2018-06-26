#include "alloc.h"
#include <math.h>

const float PI=3.14159265359;

//const float tdvx=100;                    //tomographic knot x-interval
//const float tdvz=100;                    //tomographic knot z-interval

const float nearsurface=100;            //definition of surface
const int nearoffset=5;                 //define how much is near offset
const int nlocaldip=3;                  //define how much is local for dip scanning
const int nwaveletlen=5;                //the length of the wavelet
const float maxoffegylimit=0.99;        //define the limit for scanning max effective angle
const float weight0=0.1;                //near offset weight
const float c20=-100;                   //c2 scanning start value
const float dc2=1;                      //c2 step
const int nc2=400;                      //number of dc2
const float h1=0.001;                   //large time step of raytracing 
const float h2=0.0002;                  //small time step of raytracing
const float h1to2=10;                   //the boundary of change h1 to h2
const float y2limit=1;                  //define surface for raytracing
const float eps=1e-10;                  //precision
const float tmax=5;                     //maximium recorded time
const float nmaxray=1000;                        //max number of rays for each point
const float nmaxpoint=5000;              //max number of points for each layer
const int winlen=1;                     //define the length of picking window
const float ddip=1.0;                   //dip step
const int ndip=5;                       //number of ddip
const float dtkangle=1.0;                   //take off dangle of raytracing 

struct validpoint						// store the valid reflection points info
{
        float   px;						// the x coordinate of the valid point 
        float   pz;						// the z coordinate of the valid point 
        float   dip;					// local dip of the point
//        float   c2;						// local c2 of the adcig 
        float   maxangle;				// maxmuim valid angle of the local adcig
	int	ilayer;						// the layer of the point 
	float   rz[61];					//delat z
};

/*	read Tomography parameter from par */
void readPar(char *fn1, char *fn2, char *fn3, char *fn5,char *fn6,int *nlayer, float *cig0, float *cig1, int *dpkcig, float *depth, float *tdvx, float *tdvz, float *nmx, float *nmz, float *reg1, float *reg2, int *iter, float *iflag, int *nps, float *flagtopo);

/*	scan angle gather and get the cig info */
void scanAngleGather(char *fn1, int *nz, float *dz, int *noffset, float *doffset, int *ncig, float* xa0, float* dcig);

/*	scan vel model and get the vel info */
void scanVelModel(char* fn2, float* xv0, int* nvx, int* nvz, float* dvx, float* dvz);

/*	pick cig RMO info  */
void createInvariants(char *fn1, char *fn2, char* fn3, char* fn4, int nlayer,float cig0, float cig1, int dpkcig, int ncig, float dcig, float xa0, int noffset, float doffset, int naz, float daz, int nvx, int nvz, float dvx, float dvz, float xv0, float depth, int *nvdpall ,int *nlayers, float iflag, int nps, float flagtopo);

/*	smooth function */
void smooth2d(float **v, int n1, int n2, float r1, float r2, float rw);

/*	ray tracing for time */
void rayTracing(float point_x, float point_z, float** vv, float sita1, int nvx, int nvz, float dvx, float dvz, float *yout );

/*	ray tracing for Tomography equation L */
void rayTracing2(float point_x, float point_z, float** vv, float sita1, int nvx, int nvz, float dvx, float dvz, float* p1x, float* p1z, int *np1, float* yout , float** topo, int layer);

/*	read parameter from TEMP */
void readTempPar(char* fn2, char* fn4, char* fn5, char *fn6,int* nvx, int* nvz, int* nvdpall, float *adcigmaxangle, float* xv0, float* dvx, float* dvz, float* nmx, float* nmz, float* reg1, float* reg2,int* iter, int* nlayer, float* tdvx, float* tdvz);	

/*	rewrite reflection point info for MPI */
void reWriteRfpt(int np, int nvdpall);


/*	structure tensor picking	*/
void cal_parameter(float sigma, float _n0[],float _n1[],float _n2[],float _n3[],float _d1[],float _d2[],float _d3[],float _d4[]);
void applyN(int nd,int m,float x[], float y[],float _n0[],float _n1[],float _n2[],float _n3[],float _d1[],float _d2[],float _d3[],float _d4[] );
void applyNX(int nd,int m2,int m1,float sigma, float **x, float ** y);
void applyXN(int nd,int m2,int m1,float sigma, float **x, float ** y);
void solve_sysmmetric22(float **a,float **v, float *d);
void linklinem(int ii0, int jj0, int step_num, int nps, float value, int num_trace, int num_long, float **x, float **y ,float *px, float *pz, int *ip);
void pick_peak(float **seigrm, float **process_seigrm, int num_trace, int num_long, float iflag);

/*   3drays */
void rayTracing3d(float point_x, float point_y, float point_z, float*** vv, float sita1, float azimuth1, int nvx, int nvy, int nvz, float dvx, float dvy, float dvz, float* p1x, float *p1y, float* p1z, int *np1, float* yout);
int rk4(float y[], float dydx[], int *n, float *x, float *h, float yout[], float ***vv, int nx, int ny, int nz, int dx, int dy, int dz);
int derivs(float x, float y[], float dydx[], float ***vv, int nx, int ny, int nz, int dx, int dy, int dz);
void deriv1(float ***vv, int ix, int iy, int iz, float *vv_x, float* vv_y, float *vv_z, int nx, int ny, int nz, int dx, int dy, int dz);


