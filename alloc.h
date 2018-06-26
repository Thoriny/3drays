#include "stdlib.h"
#include "stdio.h"
#include "segyhdr.h"

void *alloc1 (size_t n1, size_t size);
void **alloc2 (int n1, int n2, int size);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
float *alloc1float(size_t n1);
float **alloc2float(int n1, int n2);
float ***alloc3float(size_t n1, size_t n2, size_t n3);
int *alloc1int(size_t n1);
int **alloc2int(size_t n1, size_t n2);
char **alloc2char(int n1, int n2);

void free1 (void *p);
void free2 (void **p);
void free3 (void ***p);
void free1float(float *p);
void free2float(float **p);
void free3float(float ***p);
void free2int(int **p);
void free2char(char **p);

void zero1int(int *p, size_t n1);
void zero2int(int **p, size_t n1, size_t n2);
void zero1float(float *p, size_t n1);
void zero2float(float **p, size_t n1, size_t n2);
void zero3float(float ***p, size_t n1, size_t n2, size_t n3);


void cp2float(float **in,float **out,int n1,int n2);
void mis2float(float **in1,float **in2,float **out,int n1,int n2 );
void add2float(float **in1,float **in2,float **out,int n1,int n2);
void readfile2float(float **out,char fn[256],int n1,int n2);
void readfilesu2float(float **out,char fn[256],int n1,int n2);
void outputimage(float **data,char fn[256],int n1,int n2);
void output3float(float ***data, char fn[256], int n1, int n2, int n3);
void output1float(float *data,char fn[256],int n1);
void output2floatsu(float ** data,char fn1[256],char fn2[256],int n1,int n2);

void intepolation2d(float **image,float **map,int nix,int niz,float dix,float diz,int nvx,int nvz,float dvx,float dvz);
void intlin_(int *nin, float xin[], float yin[], float *yinl, float *yinr,	int *nout, float xout[], float yout[]);
void intlin (int nin, float xin[], float yin[], float yinl, float yinr, int nout, float xout[], float yout[]);
void xindex (int nx, float ax[], float x, int *index);



void xindex (int nx, float ax[], float x, int *index)
/*****************************************************************************
determine index of x with respect to an array of x values
******************************************************************************
Input:
nx		number of x values in array ax
ax		array[nx] of monotonically increasing or decreasing x values
x		the value for which index is to be determined
index		index determined previously (used to begin search)

Output:
index		for monotonically increasing ax values, the largest index
		for which ax[index]<=x, except index=0 if ax[0]>x;
		for monotonically decreasing ax values, the largest index
		for which ax[index]>=x, except index=0 if ax[0]<x
******************************************************************************
Notes:
This function is designed to be particularly efficient when called
repeatedly for slightly changing x values; in such cases, the index 
returned from one call should be used in the next.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/25/89
*****************************************************************************/
{
	int lower,upper,middle,step;

	/* initialize lower and upper indices and step */
	lower = *index;
	if (lower<0) lower = 0;
	if (lower>=nx) lower = nx-1;
	upper = lower+1;
	step = 1;

	/* if x values increasing */
	if (ax[nx-1]>ax[0]) {

		/* find indices such that ax[lower] <= x < ax[upper] */
		while (lower>0 && ax[lower]>x) {
			upper = lower;
			lower -= step;
			step += step;
		}
		if (lower<0) lower = 0;
		while (upper<nx && ax[upper]<=x) {
			lower = upper;
			upper += step;
			step += step;
		}
		if (upper>nx) upper = nx;

		/* find index via bisection */
		while ((middle=(lower+upper)>>1)!=lower) {
			if (x>=ax[middle])
				lower = middle;
			else
				upper = middle;
		}

	/* else, if not increasing */
	} else {

		/* find indices such that ax[lower] >= x > ax[upper] */
		while (lower>0 && ax[lower]<x) {
			upper = lower;
			lower -= step;
			step += step;
		}
		if (lower<0) lower = 0;
		while (upper<nx && ax[upper]>=x) {
			lower = upper;
			upper += step;
			step += step;
		}
		if (upper>nx) upper = nx;

		/* find index via bisection */
		while ((middle=(lower+upper)>>1)!=lower) {
			if (x<=ax[middle])
				lower = middle;
			else
				upper = middle;
		}
	}

	/* return lower index */
	*index = lower;
}



/*****************************************************************************
INTLIN - evaluate y(x) via linear interpolation of y(x[0]), y(x[1]), ...

intlin		evaluate y(x) via linear interpolation of y(x[0]), y(x[1]), ...

******************************************************************************
Function Prototype:
void intlin (int nin, float xin[], float yin[], float yinl, float yinr,
	int nout, float xout[], float yout[]);

******************************************************************************
Input:
nin		length of xin and yin arrays
xin		array[nin] of monotonically increasing or decreasing x values
yin		array[nin] of input y(x) values
yinl		value used to extraplate y(x) to left of input yin values
yinr		value used to extraplate y(x) to right of input yin values
nout		length of xout and yout arrays
xout		array[nout] of x values at which to evaluate y(x)

Output:
yout		array[nout] of linearly interpolated y(x) values

******************************************************************************
Notes:
xin values must be monotonically increasing or decreasing.

Extrapolation of the function y(x) for xout values outside the range
spanned by the xin values in performed as follows:

	For monotonically increasing xin values,
		yout=yinl if xout<xin[0], and yout=yinr if xout>xin[nin-1].

	For monotonically decreasing xin values, 
		yout=yinl if xout>xin[0], and yout=yinr if xout<xin[nin-1].

If nin==1, then the monotonically increasing case is used.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/

void intlin (int nin, float xin[], float yin[], float yinl, float yinr, 
	int nout, float xout[], float yout[])
/*****************************************************************************
evaluate y(x) via linear interpolation of y(x[0]), y(x[1]), ...
******************************************************************************
Input:
nin		length of xin and yin arrays
xin		array[nin] of monotonically increasing or decreasing x values
yin		array[nin] of input y(x) values
yinl		value used to extraplate y(x) to left of input yin values
yinr		value used to extraplate y(x) to right of input yin values
nout		length of xout and yout arrays
xout		array[nout] of x values at which to evaluate y(x)

Output:
yout		array[nout] of linearly interpolated y(x) values
******************************************************************************
Notes:
xin values must be monotonically increasing or decreasing.

Extrapolation of the function y(x) for xout values outside the range
spanned by the xin values in performed as follows:

	For monotonically increasing xin values,
		yout=yinl if xout<xin[0], and yout=yinr if xout>xin[nin-1].

	For monotonically decreasing xin values, 
		yout=yinl if xout>xin[0], and yout=yinr if xout<xin[nin-1].

If nin==1, then the monotonically increasing case is used.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	static int idx;
	int jout;
	float x;

	/* if input x values are monotonically increasing, then */
	if (xin[0]<=xin[nin-1]) {
		for (jout=0; jout<nout; jout++) {
			x = xout[jout];
			if (x<xin[0])
				yout[jout] = yinl;
			else if (x>xin[nin-1])
				yout[jout] = yinr;
			else if (x==xin[nin-1] || nin==1)
				yout[jout] = yin[nin-1];
			else {
				xindex(nin,xin,x,&idx);
				yout[jout] = yin[idx]+(x-xin[idx])
					*(yin[idx+1]-yin[idx])
					/(xin[idx+1]-xin[idx]);
			}
		}
	
	/* else, if input x values are monotonically decreasing, then */
	} else {
		for (jout=0; jout<nout; jout++) {
			x = xout[jout];
			if (x>xin[0])
				yout[jout] = yinl;
			else if (x<xin[nin-1])
				yout[jout] = yinr;
			else if (x==xin[nin-1] || nin==1)
				yout[jout] = yin[nin-1];
			else {
				xindex(nin,xin,x,&idx);
				yout[jout] = yin[idx]+(x-xin[idx])
					*(yin[idx+1]-yin[idx])
					/(xin[idx+1]-xin[idx]);
			}
		}
	}
}

void intlin_(int *nin, float xin[], float yin[], float *yinl, float *yinr,	int *nout, float xout[], float yout[]){

	intlin(*nin, xin, yin, *yinl, *yinr, *nout, xout, yout);
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 60

void intepolation2d(float **image,float **map,int nix,int niz,float dix,float diz,int nvx,int nvz,float dvx,float dvz)
{

	//image is input & map is output & we intepolate value in two dimension
	
	     // first is in z direction
		 float *xin=(float *)calloc(sizeof(float),niz);
		 float *yin=(float *)calloc(sizeof(float),niz);
		 float *xout=(float *)calloc(sizeof(float),nvz);
		 float *yout=(float *)calloc(sizeof(float),nvz);
         
		 int i,j;

		 for(i=0;i<niz;i++)
			 xin[i]=i*diz;

		 for(i=0;i<nvz;i++)
			 xout[i]=i*dvz;
         
		for(i=0;i<nix;i++)
		{
            for(j=0;j<niz;j++)
				yin[j]=image[i][j];

			intlin(niz,xin,yin,yin[0],0,nvz,xout,yout);

			for(j=0;j<nvz;j++)
				map[i][j]=yout[j];
		 
        }
        //x direction intepolation
        float *xin2=(float*)calloc(sizeof(float),nix);
		float *yin2=(float*)calloc(sizeof(float),nix);
		float *xout2=(float*)calloc(sizeof(float),nvx);
		float *yout2=(float*)calloc(sizeof(float),nvx);

		for(i=0;i<nix;i++)
			xin2[i]=i*dix;

		for(i=0;i<nvx;i++)
			xout2[i]=i*dvx;

	    for(j=0;j<nvz;j++)
		{
			for(i=0;i<nix;i++)
				yin2[i]=map[i][j];

			intlin(nix,xin2,yin2,yin2[0],0,nvx,xout2,yout2);

            for(i=0;i<nvx;i++)
				map[i][j]=yout2[i];


		}

		free(xin);
		free(yin);
		free(xout);
		free(yout);
		free(xin2);
		free(yin2);
		free(xout2);
		free(yout2);
		
}






void cp2float(float **in,float **out,int n1,int n2)
{
   int i,j;
   for(i=0;i<n2;i++)
	   for(j=0;j<n1;j++)
		   out[i][j]=in[i][j];

}

void mis2float(float **in1,float **in2,float **out,int n1,int n2 )
{
	int i,j;
   for(i=0;i<n2;i++)
	   for(j=0;j<n1;j++)
		   out[i][j]=in1[i][j]-in2[i][j];
	
}

void add2float(float **in1,float **in2,float **out,int n1,int n2)
{
	int i,j;
   for(i=0;i<n2;i++)
	   for(j=0;j<n1;j++)
		   out[i][j]=in1[i][j]+in2[i][j];

}

void readfile2float(float **out,char fn[256],int n1,int n2)
{
	int i,j;
	FILE *fpr;
	fpr=fopen(fn,"r");
	for(i=0;i<n2;i++)
		fread(out[i],sizeof(float),n1,fpr);
	fclose(fpr);

}


void readfilesu2float(float **out,char fn[256],int n1,int n2)
{
	int i,j;
	FILE *fpr;
	fpr=fopen(fn,"r");
	for(i=0;i<n2;i++)
	{
		fseek(fpr, 240, SEEK_CUR);
		fread(out[i],sizeof(float),n1,fpr);
	}
	fclose(fpr);

}

void output3float(float ***data, char fn[256], int n1, int n2, int n3)
{
	int i, j;
	FILE *fp;
	fp = fopen(fn, "wb");
		for (i = 0;i < n3;i++)
			for (j = 0;j < n2;j++)
				fwrite(data[i][j], sizeof(float), n1, fp);
	fclose(fp);
}

void outputimage(float **data,char fn[256],int n1,int n2)
{
	int ix;
	FILE *fpout;
 fpout=fopen(fn,"wb");
 for(ix=0;ix<n2;ix++)
 {
	 fwrite(data[ix],sizeof(float),n1,fpout);
 }
 fclose(fpout);
}

void output1float(float *data,char fn[256],int n1)
{
	int i;
	FILE *fp;
	fp=fopen(fn,"wb");
	for(i=0;i<n1;i++)
		fwrite(&data[i],sizeof(float),1,fp);
	fclose(fp);
}

void output2floatsu(float ** data,char fn1[256],char fn2[256],int n1,int n2)
{
	int i,j;
	segyhdr hdr;
	FILE *fp1,*fp2;
	fp1=fopen(fn1,"rb");
	fp2=fopen(fn2,"wb");
	for(i=0;i<n2;i++)
	{
		fread(&hdr,sizeof(segyhdr),1,fp1);
		fseek(fp1,sizeof(float)*hdr.ns,SEEK_CUR);
		fwrite(&hdr,sizeof(segyhdr),1,fp2);
		fwrite(data[i],sizeof(float),n1,fp2);
		

	}
	fclose(fp1);
	fclose(fp2);





}

int *alloc1int(size_t n1)
{
        return (int*)alloc1(n1,sizeof(int));
}

void free2int(int **p)
{
        free2((void**)p);
}

int **alloc2int(size_t n1, size_t n2)
{
        return (int**)alloc2(n1,n2,sizeof(int));
}

float ***alloc3float(size_t n1, size_t n2, size_t n3)
{
        return (float***)alloc3(n1,n2,n3,sizeof(float));
}

void free3float(float ***p)
{
        free3((void***)p);
}

void zero1float(float *p, size_t n1)
{
     int i;
     for(i=0;i<n1;i++) p[i]=0.0;
}

void zero2float(float **p, size_t n1, size_t n2)
{
     int i, j;
     for(i=0;i<n2;i++)
       for(j=0;j<n1;j++)
         p[i][j]=0.0;
}

void zero3float(float ***p, size_t n1, size_t n2, size_t n3)
{
     int i, j, k;
     for(i=0;i<n3;i++)
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0.0;
}


void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
        size_t i3,i2;
        void ***p;

        if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
                return NULL;
        if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
                free(p);
                return NULL;
        }
        if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
                free(p[0]);
                free(p);
                return NULL;
        }
        for (i3=0; i3<n3; i3++) {
                p[i3] = p[0]+n2*i3;
                for (i2=0; i2<n2; i2++)
                        p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
        }
        return p;
}

void free3 (void ***p)
{
        free(p[0][0]);
        free(p[0]);
        free(p);
}

void zero1int(int *p, size_t n1)
{
     int i;
     for(i=0;i<n1;i++) p[i]=0;
}

void zero2int(int **p, size_t n1, size_t n2)
{
     int i, j;
     for(i=0;i<n2;i++)
       for(j=0;j<n1;j++)
         p[i][j]=0;
}

void free1 (void *p)
{
        free(p);
}

void free1float(float *p)
{
        free1(p);
}

float *alloc1float(size_t n1)
{
        return (float*)alloc1(n1,sizeof(float));
}

void *alloc1 (size_t n1, size_t size)
{
        void *p;

        if ((p=malloc(n1*size))==NULL)
                return NULL;
        return p;
}


void **alloc2 (int n1, int n2, int size)
{
        int i;
        void **p;

        p=(void**)calloc(sizeof(void*),n2);
        p[0]=(void*)calloc(size, n1*n2);
        for (i=0; i<n2; i++)
                p[i] = (char*)p[0]+size*n1*i;
        return p;
}

void free2 (void **p)
{
        free(p[0]);
        free(p);
}
char **alloc2char(int n1, int n2)
{
        return (char**)alloc2(n1,n2,sizeof(char));
}

void free2char(char **p)
{
        free2((void**)p);
}

float **alloc2float(int n1, int n2)
{
        return (float**)alloc2(n1,n2,sizeof(float));
}

void free2float(float **p)
{
        free2((void**)p);
}


