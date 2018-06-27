/*=================================================================*
*==================================================================*
*                    3d  ray tracing in isotropic media            *
*                                                                  *
*==================================================================*
*       School of Ocean & Earth Science, Tongji University, China  *
*==================================================================*
*         		Author: Yu yang                                    *
*                             2018-6-26                            *
*==================================================================*/

#include"tomography.h"

int main(int argc, char **argv)
{
	float x, y, z;
	int i, j, k,ang,iazimuth,np1;
	int ix, iy, iz;
	float *p1x, *p1y, *p1z,*yout;
	float theta,azimuth;
	float ***v,***line;

	int nvx = 100;
	int nvy = 100;
	int nvz = 100;

	float dvx=10.0;
	float dvy=10.0;
	float dvz=10.0;


	int pmax = (int)(tmax / h1);	//max size of every ray
	v = alloc3float(nvz, nvy, nvx);
	line = alloc3float(nvz, nvy, nvx);
	p1x = alloc1float(pmax);
	p1y = alloc1float(pmax);
	p1z = alloc1float(pmax);
	yout = alloc1float(8);


	zero3float(line, nvz, nvy, nvx);


	for (i=0;i<nvx;i++)
		for(j=0;j<nvy;j++)
			for (k = 0;k < nvz;k++)
			{
				if (k < nvz / 2)
					v[i][j][k] = 2000.0;
				else
					v[i][j][k] = 2000.0;
			}

	output3float(v, "vel.dat",nvz, nvy, nvx);
	x = 50*dvx;
	y = 50*dvy;
	z = 50*dvz;
   for(iazimuth=0;iazimuth<22;iazimuth++)
	for (ang = 0; ang< 90;ang++)
	{
		theta = ang * 4 * PI / 180;
		azimuth = iazimuth * 4 * PI / 180;
		rayTracing3d(x, y, z, v, theta, azimuth, nvx, nvy, nvz, dvx, dvy, dvz, p1x, p1y, p1z, &np1, yout);

         
		printf("np1=%d\n", np1);

		for (i = 0;i < np1;i++)
		{
			ix = (int)(p1x[i] / dvx);
			iy = (int)(p1y[i] / dvy);
			iz = (int)(p1z[i] / dvz);
			
			line[ix][iy][iz] =20 ;
			printf("ix=%d iy=%d iz=%d line=%f\n", ix, iy, iz, line[ix][iy][iz]);
		}
	}
	output3float(line, "line.dat", nvz, nvy, nvx);





	return 0;
}




void rayTracing3d(float point_x, float point_y, float point_z, float*** vv, float sita1,float azimuth1, int nvx, int nvy,int nvz,float dvx,float dvy, float dvz, float* p1x, float *p1y,float* p1z, int *np1, float* yout )
{
        int     ii, ix, iy,iz, in;
        int     n=5;
        float   x=0.0, h;
        float   *y, *dydx ;
        const int nmax=999999;

        y=(float*)calloc(sizeof(float),8);
        dydx=(float*)calloc(sizeof(float),8);

        ix=(int)(point_x/dvx);
        iy=(int)(point_y/dvy);
        iz=(int)(point_z/dvz);

        y[1]=point_x;
        y[2]=point_y;
        y[3]=point_z;
        y[4]=sin(sita1)*cos(azimuth1)/vv[ix][iy][iz];
        y[5]=sin(sita1)*sin(azimuth1)/vv[ix][iy][iz];
        y[6]=cos(sita1)/vv[ix][iy][iz];
        y[7]=0;
        x=0.0;
        n=7;
        *np1=0;

        for(in=0;in<=nmax;in++)
        {
                if( fabs( y[2]-0.0 )<10.0 )
                        h=h2;
                else
                {
                        p1x[in]=y[1];
                        p1y[in]=y[2];
                        p1z[in]=y[3];
                        *np1=*np1+1;
                        h=h1;
                }

                if(!derivs(x,y,dydx,vv, nvx, nvy, nvz, dvx,dvy, dvz))
                        break;

                if(!rk4(y,dydx,&n,&x,&h,yout, vv, nvx, nvy, nvz, dvx, dvy, dvz))
                        break;

                for(ii=1;ii<=7;ii++)
                        y[ii]=yout[ii];
        }

        p1x[*np1]=y[1];
        p1y[*np1]=y[2];
        p1z[*np1]=y[3];
        *np1=*np1+1;

        free(y);
        free(dydx);
}
/***********************************************************************/
/*****       using Runge-Kutta fourth order integral              ******/
int rk4(float y[], float dydx[], int *n, float *x, float *h, float yout[], float ***vv, int nx, int ny, int nz, int dx, int dy, int dz)
{
	int i;
	float yt[11], dyt[11], dym[11];
	float hh, h6, xh;
	hh = (*h)*0.5;
	h6 = (*h) / 6.0;
	xh = (*x) + hh;

	for (i = 1;i <= *n;i++)
	{
		yt[i] = y[i] + hh*dydx[i];
	}

	if (!derivs(xh, yt, dyt, vv, nx, ny, nz, dx, dy, dz))
		return 0;

	for (i = 1;i <= *n;i++)
	{
		yt[i] = y[i] + hh*dyt[i];
	}

	if (!derivs(xh, yt, dym, vv, nx, ny, nz, dx, dy, dz))
		return 0;

	for (i = 1;i <= *n;i++)
	{
		yt[i] = y[i] + (*h)*dym[i];
		dym[i] = dyt[i] + dym[i];
	}

	if (!derivs((*x) + (*h), yt, dyt, vv, nx, ny, nz, dx, dy, dz))
		return 0;

	for (i = 1;i <= *n;i++)
	{
		yout[i] = y[i] + h6*(dydx[i] + dyt[i] + 2.0*dym[i]);
	}
	for (i = 1;i <= 10;i++)
	{
		dym[i] = 0.0;
		dyt[i] = 0.0;
		yt[i] = 0.0;
	}

	return 1;
}
/****************************************************************************/
/******     calculate derivatives dydx                              *********/
int derivs(float x, float y[], float dydx[], float ***vv, int nx, int ny, int nz, int dx, int dy, int dz)
{
	int ix, iy, iz;
	float vv_x, vv_y, vv_z;
	ix = (int)(y[1] / dx);
	iy = (int)(y[2] / dy);
	iz = (int)(y[3] / dz);
	if (y[1]<0 || y[1] >= dx*nx || y[2]<0 || y[2] >= dy*ny || y[3]<0 || y[3] >= dz*nz)
	{
		return 0;
	}

	deriv1(vv, ix, iy, iz, &vv_x, &vv_y, &vv_z, nx, ny, nz, dx, dy, dz);

	dydx[1] = vv[ix][iy][iz] * vv[ix][iy][iz] * y[4];     //dx(s)/dt

	dydx[2] = vv[ix][iy][iz] * vv[ix][iy][iz] * y[5]; //dy(s)/dt

	dydx[3] = vv[ix][iy][iz] * vv[ix][iy][iz] * y[6]; //dz(s)/dt

	dydx[4] = -(1.0 / vv[ix][iy][iz])*vv_x;     //dpx(s)/dt

	dydx[5] = -(1.0 / vv[ix][iy][iz])*vv_y;     //dpy(s)/dt

	dydx[6] = -(1.0 / vv[ix][iy][iz])*vv_z;     //dpz(s)/dt

	dydx[7] = y[4] * dydx[1] + y[5] * dydx[2] + y[6] * dydx[3];    //v*v*p*p

	return 1;
}
/****************************************************************************/
/******     calculate derivatives dv/dx  dv/dy                        ******/
void deriv1(float ***vv, int ix, int iy, int iz, float *vv_x, float* vv_y, float *vv_z, int nx, int ny, int nz, int dx, int dy, int dz)
{
	if (ix == 0)
		*vv_x = (vv[ix + 1][iy][iz] - vv[ix][iy][iz]) / dx;
	else if (ix == nx - 1)
		*vv_x = (vv[ix][iy][iz] - vv[ix - 1][iy][iz]) / dx;
	else if (ix>0 && ix<(nx - 1))
		*vv_x = (vv[ix + 1][iy][iz] - vv[ix - 1][iy][iz]) / (2.0*dx);
	if (iy == 0)
		*vv_y = (vv[ix][iy + 1][iz] - vv[ix][iy][iz]) / dy;
	else if (iy == ny - 1)
		*vv_y = (vv[ix][iy][iz] - vv[ix][iy - 1][iz]) / dx;
	else if (iy>0 && iy<(ny - 1))
		*vv_x = (vv[ix][iy + 1][iz] - vv[ix][iy - 1][iz]) / (2.0*dy);
	if (iz == 0)
		*vv_z = (vv[ix][iy][iz + 1] - vv[ix][iy][iz]) / dz;
	else if (iz == nz - 1)
		*vv_z = (vv[ix][iy][iz] - vv[ix][iy][iz - 1]) / dz;
	else if (iz>0 && iz<(nz - 1))
		*vv_z = (vv[ix][iy][iz + 1] - vv[ix][iy][iz - 1]) / (2.0*dz);
}
