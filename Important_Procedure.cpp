// Important_Procedure.cpp
// 2017/5/22 : 2D output is added
// 2017/5/25 : In order to spare stack space, change arguments of function from class type to class referrence type, for example, in 'cal_flux' and 'set_dt' and 'step_on'

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
#include "Variables_Definition.h"
#include "Procedure.h"

// Clarification of mesh-grid
extern double X[], Y[], Z[], X_interval[], Y_interval[], Z_interval[];
double sub_var[8][Grid_Num_x][Grid_Num_y][Grid_Num_z];  //subsidiary_variable
double var_x[8][Grid_Num_x], var_x_plushalfdx[8][Grid_Num_x];          // fx, fxi, used in 2-Lax_Wendroff stepon
// For following procedures

void set_mesh()
{
	int nx,ny,nz;
	X[0]=x_min;
	Y[0]=y_min;
	Z[0]=z_min;
	nx=Grid_Num_x-1;
	ny=Grid_Num_x-1;
	nz=Grid_Num_z-1;
	ofstream out("grid.dat");
	ofstream out_2D("grid_2D.dat");

	out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<X[0];
	out_2D<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<X[0];
	int i,j,k;
	if (uniform_x==True)
	{
		for (i=0;i<nx;i++)
		{			
			X_interval[i]=(x_max-x_min)/nx;	
			X[i+1]=X[i]+X_interval[i];
			out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<X[i+1];
			out_2D<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<X[i+1];
		}
	}
	else
	{
		for (i=0;i<nx;i++)
		{
			X_interval[i]=0.01;
			X[i+1]=X[i]+X_interval[i];
			out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<X[i+1];
			out_2D<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<X[i+1];
		}
	}
	X_interval[Grid_Num_x-1]=X_interval[Grid_Num_x-2];
	out<<endl;
	out_2D<<endl;

	out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Y[0];
	if (uniform_y==True)
	{
		for (j=0;j<ny;j++)
		{			
			Y_interval[j]=(y_max-y_min)/ny;	
			Y[j+1]=Y[j]+Y_interval[j];
			out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Y[j+1];
		}
	}
	else
	{
		for (j=0;j<ny;j++)
		{			
			Y_interval[j]=0.01;
			Y[j+1]=Y[j]+Y_interval[j];
			out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Y[j+1];
		}
	}
	Y_interval[Grid_Num_y-1]=Y_interval[Grid_Num_y-2];
	out<<endl;
	
	out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Z[0];
	out_2D<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Z[0];
	if (uniform_z==True)
	{
		for (k=0;k<nz;k++)
		{			
			Z_interval[k]=(z_max-z_min)/nz;	
			Z[k+1]=Z[k]+Z_interval[k];
			out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Z[k+1];
			out_2D<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Z[k+1];
		}
	}
	else
	{
		for (k=0;k<nz;k++)
		{			
			Z_interval[k]=0.01;
			Z[k+1]=Z[k]+Z_interval[k];
			out<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Z[k+1];
			out_2D<<" "<<setprecision(3)<<setiosflags(ios::fixed)<<Z[k+1];
		}
	}
	Z_interval[Grid_Num_z-1]=Z_interval[Grid_Num_z-2];
	out<<endl;
	out_2D<<endl;

	out.close();
	//cout<<"BASIC_VARIABLE::set_mesh invoked!"<<endl;
	// Set mesh-grid				  
}

/* 
   Initialising variables and setting mesh-grid, where
   pointer points to a variable array, the Main Variables,
   and pressure_obj is a referrence to Pressure variable.
   This function will initialize the Main Variables pointed by pointer,
   and Pressure variable referred by pressure_obj.
*/
void initialize(VARIABLE *pointer, BASIC_VARIABLE &pressure_obj)
{	
	double rm,rs,bm,bs,v0,betam, pressuremtotal;
// Following used are global variables listed in "Basic_Parameter.h"
	rm=rho_m_0;	rs=rho_s_0;
	bm=B_m_0; bs=B_s_0; 
	v0=v_0;
	betam=beta_m; 
	pressuremtotal=betam*0.5*pow(bm,2)+0.5*pow(bm,2);      //betam*0.5*pow(bm,2)=pmsp

	int i,j,k,n;
	double x,y,z, dx,dy,dz, r5;
	double rho, Bx, By, Bz, rhoVx, rhoVy, rhoVz;
	double B_Energy, V_Energy, pressure;
	for (i=0;i<Grid_Num_x;i++)
	{
		for (j=0;j<Grid_Num_y;j++)
		{
			for (k=0;k<Grid_Num_z;k++)
			{
				x=X[i]; dx=X_interval[i];
				y=Y[j]; dy=Y_interval[j];
				z=Z[k]; dz=Z_interval[k];
				rho=0.5*(rm+rs);				
				r5=pow(x*x+y*y+z*z,5/2);				  /// attention here 5/2=2
				Bx=-bm*3*x*z/r5;
				By=-bm*3*y*z/r5;
				Bz=-bm*(2*pow(z,2)-pow(x,2)-pow(y,2))/r5;
				pointer[0].value[i][j][k]=rho;
				pointer[4].value[i][j][k]=Bx;
				pointer[5].value[i][j][k]=By;
				pointer[6].value[i][j][k]=Bz;				
				rhoVx=0;
				rhoVy=0.5*v0*By*(1-tanh(x/width_rho))*rho/bs;
				rhoVz=0.5*v0*Bz*(1-tanh(x/width_rho))*rho/bs;
				pointer[1].value[i][j][k]=rhoVx;
				pointer[2].value[i][j][k]=rhoVy;
				pointer[3].value[i][j][k]=rhoVz;				
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				pressure=pressuremtotal-B_Energy;
				pointer[7].value[i][j][k]=B_Energy+V_Energy+pressure/(phy_gamma-1);
				pressure_obj.value[i][j][k]=pressure;
				/*plus half dx: start*/
				rho=0.5*(rm+rs);
				r5=pow(x*x+y*y+z*z,5/2);	      /// attention here 5/2=2			
				Bx=-bm*3*x*z/r5;
				By=-bm*3*y*z/r5;
				Bz=-bm*(2*pow(z,2)-pow(x,2)-pow(y,2))/r5;
				sub_var[0][i][j][k]=rho;
				sub_var[4][i][j][k]=Bx;
				sub_var[5][i][j][k]=By;
				sub_var[6][i][j][k]=Bz;
				rhoVx=0;
				rhoVy=0.5*v0*By*(1-tanh((x+dx/2.)/width_rho))*rho/bs;
				rhoVz=0.5*v0*Bz*(1-tanh((x+dx/2.)/width_rho))*rho/bs;
				sub_var[1][i][j][k]=rhoVx;
				sub_var[2][i][j][k]=rhoVy;
				sub_var[3][i][j][k]=rhoVz;
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				pressure=pressuremtotal-B_Energy;
				sub_var[7][i][j][k]=B_Energy+V_Energy+pressure/(phy_gamma-1);
				/*plus half dx: end*/

				/*vyi0=0, no current correction*/
			}
		}
	}	
	for (n=0;n<8;n++)
	{
		for (i=0;i<Grid_Num_x;i++)
		{
			var_x[n][i]=pointer[n].value[i][0][0];
			var_x_plushalfdx[n][i]=sub_var[n][i][0][0];
		}
	}
	for (n=0;n<8;n++)
		for (i=0;i<Grid_Num_x;i++)
			for (j=0;j<Grid_Num_y;j++)
				for (k=0;k<Grid_Num_z;k++)
					sub_var[n][i][j][k]=0;
	//cout<<"Initialize invoked! But I really don't know the setup written by Teacher Ma! Waiting to be changed to a symmetric Harris Current Sheet!"<<endl;
}


void cal_current(VARIABLE *current, VARIABLE *pointer, Type T)
{
	double Bx_ym, Bx_zm, By_zm, By_xm, Bz_xm, Bz_ym;
	double Bx_yp, Bx_zp, By_zp, By_xp, Bz_xp, Bz_yp;
	double dx, dy, dz;
	int i,j,k, n;
	for (i=1;i<Grid_Num_x-1-T;i++)
	{
		for (j=1;j<Grid_Num_y-1-T;j++)
		{
			for (k=1;k<Grid_Num_z-1-T;k++)
			{
				Bx_ym=pointer[4].value[i][j-1][k]; Bx_yp=pointer[4].value[i][j+1][k];
				Bx_zm=pointer[4].value[i][j][k-1]; Bx_zp=pointer[4].value[i][j][k+1];

				By_zm=pointer[5].value[i][j][k-1]; By_zp=pointer[5].value[i][j][k+1];
				By_xm=pointer[5].value[i-1][j][k]; By_xp=pointer[5].value[i+1][j][k];

				Bz_xm=pointer[6].value[i-1][j][k]; Bz_xp=pointer[6].value[i+1][j][k];
				Bz_ym=pointer[6].value[i][j-1][k]; Bz_yp=pointer[6].value[i][j+1][k];

				dx=X_interval[i-1]+X_interval[i];
				dy=Y_interval[j-1]+Y_interval[j];
				dz=Z_interval[k-1]+Z_interval[k];

				current[0].value[i][j][k]=(Bz_yp-Bz_ym)/dy- \
					(By_zp-By_zm)/dz;
				current[1].value[i][j][k]=(Bx_zp-Bx_zm)/dz- \
					(Bz_xp-Bz_xm)/dx;
				current[2].value[i][j][k]=(By_xp-By_xm)/dx- \
					(Bx_yp-Bx_ym)/dy;
			}
		}
	}
	// (can write a switch statement to include different extrapolation method
	// x-linear extrapolation yz-equivalent
	// boundary at i=0,Grid_Num_x-1-T
	if(T==Complete)                        // linear extrapolation
	{
		// boundary at i=0,Grid_Num_x-1
		if (half_x==True)
		{
			for (j=1;j<Grid_Num_y-1;j++)
			{
				for(k=1;k<Grid_Num_z-1;k++)
				{
					current[0].value[0][j][k]=2*current[0].value[1][j][k]-current[0].value[2][j][k];
					current[0].value[Grid_Num_x-1][j][k]=current[0].value[Grid_Num_x-3][Grid_Num_y-1-j][k];
					current[1].value[0][j][k]=2*current[1].value[1][j][k]-current[1].value[2][j][k];
					current[1].value[Grid_Num_x-1][j][k]=current[1].value[Grid_Num_x-3][Grid_Num_y-1-j][k];
					current[2].value[0][j][k]=2*current[2].value[1][j][k]-current[2].value[2][j][k];
					current[2].value[Grid_Num_x-1][j][k]=-current[2].value[Grid_Num_x-3][Grid_Num_y-1-j][k];
					
				}
			}
		}
		else
		{
			for (n=0;n<3;n++)
			{
				for (j=1;j<Grid_Num_y-1;j++)
				{
					for(k=1;k<Grid_Num_z-1;k++)
					{
						current[n].value[0][j][k]=2*current[n].value[1][j][k]-current[n].value[2][j][k];
						current[n].value[Grid_Num_x-1][j][k]=2*current[n].value[Grid_Num_x-2][j][k]- \
									current[n].value[Grid_Num_x-3][j][k];						
					}
				}
			}			
		}	

		// boundary at j=0,Grid_Num_y-1
		if (period_y==True)
		{
			for (n=0;n<3;n++)
			{
				for (i=0;i<Grid_Num_x;i++)
				{
					for(k=1;k<Grid_Num_z-1;k++)
					{
						current[n].value[i][0][k]=current[n].value[i][Grid_Num_y-2][k];
						current[n].value[i][Grid_Num_y-1][k]=current[n].value[i][1][k];
					}
				}
			}
		}
		else
		{
			for (n=0;n<3;n++)
			{
				for (i=0;i<Grid_Num_x;i++)
				{
					for(k=1;k<Grid_Num_z-1;k++)
					{
						current[n].value[i][0][k]=current[n].value[i][1][k];
						current[n].value[i][Grid_Num_y-1][k]=current[n].value[i][Grid_Num_y-2][k];
					}
				}
			}
		}

		// boundary at k=0,Grid_Num_z-1
		if (half_z==True)
		{
			for (i=0;i<Grid_Num_x;i++)
			{
				for (j=0;j<Grid_Num_z;j++)
				{
					current[0].value[i][j][0]=current[0].value[i][j][1];
					current[0].value[i][j][Grid_Num_z-1]=current[0].value[i][j][Grid_Num_y-3];
					current[1].value[i][j][0]=current[1].value[i][j][1];
					current[1].value[i][j][Grid_Num_z-1]=current[2].value[i][j][Grid_Num_y-3];
					current[2].value[i][j][0]=current[n].value[i][j][1];
					current[2].value[i][j][Grid_Num_z-1]=-current[2].value[i][j][Grid_Num_y-3];
				}
			}
		}
		else
		{
			for (n=0;n<3;n++)
			{
				for (i=0;i<Grid_Num_x;i++)
				{
					for (j=0;j<Grid_Num_z;j++)
					{
						current[n].value[i][j][0]=current[n].value[i][j][1];
						current[n].value[i][j][Grid_Num_z-1]=current[n].value[i][j][Grid_Num_y-2];
					}
				}
			}
		}
	}


	else // T==Incomplete                              //x-direction equvalent extrapolation
	{
		// boundary at i=0,Grid_Num_x-1-T
		if (half_x==True)
		{
			for (j=1;j<Grid_Num_y-1-T;j++)
			{
				for(k=1;k<Grid_Num_z-1-T;k++)
				{
					current[0].value[0][j][k]=current[0].value[1][j][k];
					current[0].value[Grid_Num_x-1-T][j][k]=current[0].value[Grid_Num_x-2-T][Grid_Num_y-2-j][k];
					current[1].value[0][j][k]=current[1].value[1][j][k];
					current[1].value[Grid_Num_x-1-T][j][k]=current[1].value[Grid_Num_x-2-T][Grid_Num_y-2-j][k];
					current[2].value[0][j][k]=current[2].value[1][j][k];
					current[2].value[Grid_Num_x-1-T][j][k]=-current[2].value[Grid_Num_x-2-T][Grid_Num_y-2-j][k];					
				}
			}
		}
		else
		{
			for (n=0;n<3;n++)
			{
				for (j=1;j<Grid_Num_y-1-T;j++)
				{
					for(k=1;k<Grid_Num_z-1-T;k++)
					{
						current[n].value[0][j][k]=current[n].value[1][j][k];
						current[n].value[Grid_Num_x-1-T][j][k]=current[n].value[Grid_Num_x-2-T][j][k];	
					}
				}
			}
		}	

		// boundary at j=0,Grid_Num_y-1-T
		if (period_y==True)
		{
			for (n=0;n<3;n++)
			{
				for (i=0;i<Grid_Num_x-T;i++)
				{
					for(k=1;k<Grid_Num_z-1-T;k++)
					{
						current[n].value[i][0][k]=current[n].value[i][Grid_Num_y-2-T][k];
						current[n].value[i][Grid_Num_y-1-T][k]=current[n].value[i][1][k];
					}
				}
			}
		}
		else
		{
			for (n=0;n<3;n++)
			{
				for (i=0;i<Grid_Num_x-T;i++)
				{
					for(k=1;k<Grid_Num_z-1-T;k++)
					{
						current[n].value[i][0][k]=current[n].value[i][1][k];
						current[n].value[i][Grid_Num_y-1-T][k]=current[n].value[i][Grid_Num_y-2-T][k];
					}
				}
			}
		}

		// boundary at k=0,Grid_Num_z-1-T
		if (half_z==True)
		{
			for (i=0;i<Grid_Num_x-T;i++)
			{
				for (j=0;j<Grid_Num_z-T;j++)
				{
					current[0].value[i][j][0]=current[0].value[i][j][1];
					current[0].value[i][j][Grid_Num_z-1-T]=current[0].value[i][j][Grid_Num_y-2-T];
					current[1].value[i][j][0]=current[1].value[i][j][1];
					current[1].value[i][j][Grid_Num_z-1-T]=current[1].value[i][j][Grid_Num_y-2-T];
					current[2].value[i][j][0]=current[2].value[i][j][1];
					current[2].value[i][j][Grid_Num_z-1-T]=-current[2].value[i][j][Grid_Num_y-2-T];
				}
			}
		}
		else
		{
			for (n=0;n<3;n++)
			{
				for (i=0;i<Grid_Num_x-T;i++)
				{
					for (j=0;j<Grid_Num_z-T;j++)
					{
						current[n].value[i][j][0]=current[n].value[i][j][1];
						current[n].value[i][j][Grid_Num_z-1-T]=current[n].value[i][j][Grid_Num_y-2-T];
					}
				}
			}
		}
	}					
	//cout<<"Current is calculated form Curl B! It's simple, I know!!!"<<endl;	
}

void cal_pressure(BASIC_VARIABLE &pressure_obj, VARIABLE *pointer, Type T)
{	
	double rho,rhoVx, rhoVy, rhoVz, Bx, By, Bz, E;
	double B_Energy, V_Energy;
	int i,j,k;
	for (i=0;i<Grid_Num_x-T;i++)
	{
		for (j=0;j<Grid_Num_y-T;j++)
		{
			for (k=0;k<Grid_Num_z-T;k++)
			{
				rho=pointer[0].value[i][j][k];
				rhoVx=pointer[1].value[i][j][k];
				rhoVy=pointer[2].value[i][j][k];
				rhoVz=pointer[3].value[i][j][k];
				Bx=pointer[4].value[i][j][k];
				By=pointer[5].value[i][j][k];
				Bz=pointer[6].value[i][j][k];
				E=pointer[7].value[i][j][k];
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				pressure_obj.value[i][j][k]=(phy_gamma-1)* \
					(E-B_Energy-V_Energy);
			}
		}
	}
	make_pressure_positive(pressure_obj, 1e-5);
	//cout<<"Pressure is calculated from various kinds of energy, and doesn't explicitly step on according to equation!"<<endl;
}

// Homogeneous eta, waiting to be amended.
void set_eta(BASIC_VARIABLE &eta_obj, VARIABLE *pointer, VARIABLE *current, double time, Type T)
{
	double etam[Grid_Num_y];                                                 // in l-m-n coordinates-systme m-direction dependent value, in another way eta_y
	double etax, etaz, etal=0., etab=0.005, alpha0=2.0;	                     // etab for eta_background
	double ytrig=15., xtrig=0.0, ztrig=0.0;                                  // ytrig for xlen 
	double widthx=width_rho, widthz=2*width_rho;
	double x,y,z, dx,dy,dz;
	int i,j,k;
	if(T==Complete)
	{		
		for (i=0;i<Grid_Num_x-T;i++)
		{
			x=X[i];
			if(abs(x-xtrig)<=0.5)
				etax=1.;
			else
				etax=1.-pow(tanh((x-xtrig)/widthx),2);            // This is ture
			for (j=0;j<Grid_Num_y-T;j++)
			{
				y=Y[j];
				if(abs(y)<=ytrig)                    // ytrig is somewhat length,and it's for ytrig
					etam[j]=1;                      // This is True
				else
					etam[j]=1.-pow(tanh(abs(y)-ytrig),2);
				for (k=0;k<Grid_Num_z-T;k++)
				{
					z=Z[k];
					if(abs(abs(z)-ztrig)<=1.)
						etaz=1.;                                  // This is True
					else
						etaz=1.-pow(tanh((abs(z)-ztrig)/widthz),2);
					eta_obj.value[i][j][k]=etab+etam[j]*etal*etax*etaz;
				}
			}
		}
	}
	else
	{
		for (i=0;i<Grid_Num_x-T;i++)
		{
			x=X[i]; dx=X_interval[i];
			if(abs(x+dx/2.-xtrig)<=0.5)
				etax=1.;
			else
				etax=1.-pow(tanh((x+dx/2.-xtrig)/widthx),2);            // This is True
			for (j=0;j<Grid_Num_y-T;j++)
			{
				y=Y[j]; dy=Y_interval[j];    
				if(abs(y+dy/2.)<=ytrig)                    
					etam[j]=1;                      // This is True
				else
					etam[j]=1.-pow(tanh(abs(y+dy/2.)-ytrig),2);
				for (k=0;k<Grid_Num_z-T;k++)
				{
					z=Z[k]; dz=Z_interval[k];
					if(abs(abs(z+dz/2.)-ztrig)<=1.)
						etaz=1.;                                  // This is True
					else
						etaz=1.-pow(tanh(((z+dz/2.)-ztrig)/widthz),2);     // ????????? bot tanh((abs(zz(jz)+dz/2.)-ztrig)/aw2)?????
					eta_obj.value[i][j][k]=etab+etam[j]*etal*etax*etaz;
				}
			}
		}
	}	
	//cout<<"Set_eta invoked! But for now eta is simply setted to 0.1! And I'm sorry to inform you that you should change it!"<<endl;
}

// Calculating flux from variables
void cal_flux(BASIC_VARIABLE flux[][3], VARIABLE *pointer, VARIABLE *current,\
	BASIC_VARIABLE &pressure_obj, BASIC_VARIABLE &eta_obj, Type T)
{
	double rho, Vx, Vy, Vz;
	double Bx, By, Bz;
	double Energy, B_Energy, V_Energy;  // V_Energy is not used !
	double pressure, eta;
	double current_x, current_y, current_z;
	double V_dot_B, V_cross_B_x, V_cross_B_y, V_cross_B_z;
	// Assignment
	int i,j,k;
	for (i=0;i<Grid_Num_x-T;i++)
	{
		for (j=0;j<Grid_Num_y-T;j++)
		{
			for (k=0;k<Grid_Num_z-T;k++)
			{
				rho=pointer[0].value[i][j][k];
				Vx=pointer[1].value[i][j][k]/rho;
				Vy=pointer[2].value[i][j][k]/rho;
				Vz=pointer[3].value[i][j][k]/rho;
				Bx=pointer[4].value[i][j][k];
				By=pointer[5].value[i][j][k];
				Bz=pointer[6].value[i][j][k];
				Energy=pointer[7].value[i][j][k];
				pressure=pressure_obj.value[i][j][k];
				current_x=current[0].value[i][j][k];
				current_y=current[1].value[i][j][k];
				current_z=current[2].value[i][j][k];
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*rho*(pow(Vx,2)+pow(Vy,2)+pow(Vz,2));
				eta=eta_obj.value[i][j][k];

				// rho_flux assignment
				flux[0][0].value[i][j][k]=rho*Vx;
				flux[0][1].value[i][j][k]=rho*Vy;
				flux[0][2].value[i][j][k]=rho*Vz;
				
				//Bx_flux
				flux[1][0].value[i][j][k]=rho*Vx*Vx+pressure+B_Energy-Bx*Bx;
				flux[1][1].value[i][j][k]=rho*Vx*Vy-By*Bx;
				flux[1][2].value[i][j][k]=rho*Vx*Vz-Bz*Bx;

				flux[2][0].value[i][j][k]=rho*Vy*Vx-Bx*By;
				flux[2][1].value[i][j][k]=rho*Vy*Vy+pressure+B_Energy-By*By;
				flux[2][2].value[i][j][k]=rho*Vy*Vz-Bz*By;

				flux[3][0].value[i][j][k]=rho*Vz*Vx-Bx*Bz;
				flux[3][1].value[i][j][k]=rho*Vz*Vy-By*Bz;
				flux[3][2].value[i][j][k]=rho*Vz*Vz+pressure+B_Energy-Bz*Bz;

				//Magnetic Induction Eq.
				V_cross_B_x=(Vy-di*current_y)*Bz-(Vz-di*current_z)*By;                     // don't take di into acount, that is di=0
				V_cross_B_y=(Vz-di*current_z)*Bx-(Vx-di*current_x)*Bz;
				V_cross_B_z=(Vx-di*current_x)*By-(Vy-di*current_y)*Bx;

				flux[4][0].value[i][j][k]=0;
				flux[4][1].value[i][j][k]=-V_cross_B_z+eta*current_z;
				flux[4][2].value[i][j][k]=V_cross_B_y-eta*current_y;

				flux[5][0].value[i][j][k]=V_cross_B_z-eta*current_z;
				flux[5][1].value[i][j][k]=0;
				flux[5][2].value[i][j][k]=-V_cross_B_x+eta*current_x;

				flux[6][0].value[i][j][k]=-V_cross_B_y+eta*current_y;
				flux[6][1].value[i][j][k]=V_cross_B_x-eta*current_x;
				flux[6][2].value[i][j][k]=0;

				// Energy_fulx
				V_dot_B=Vx*Bx+Vy*By+Vz*Bz;
				Energy=Energy+pressure+B_Energy;       // Energy flux brought by V
				flux[7][0].value[i][j][k]=Vx*Energy-Bx*V_dot_B+eta* \
					(current_y*Bz-current_z*By);
				flux[7][1].value[i][j][k]=Vy*Energy-By*V_dot_B+eta* \
					(current_z*Bx-current_x*Bz);
				flux[7][2].value[i][j][k]=Vz*Energy-Bz*V_dot_B+eta* \
					(current_x*By-current_y*Bx);
			}
		}
	}
	//cout<<"I'm Calculating flux!"<<endl;
}

// Setting time-interval
double set_dt(VARIABLE *pointer, BASIC_VARIABLE &eta_obj, VARIABLE *current, BASIC_VARIABLE &pressure_obj, double time)
{
	//set_eta(eta_obj, pointer, current, time);
	double dt;
	double dx,dy,dz,dxyz;	
	double Temp_dt, dtmin=1000., test_dt=1000.;
	double rho, rhoVx, rhoVy, rhoVz, Bx, By, Bz, pressure;
	int i,j,k;
	for (i=0;i<Grid_Num_x;i++)
	{
		for (j=0;j<Grid_Num_y;j++)
		{
			for (k=0;k<Grid_Num_z;k++)
			{
				dx=X_interval[i];
				dy=Y_interval[j];
				dz=Z_interval[k];
				dxyz=0.5*dx*dy*dz/ \
					sqrt(dx*dx*dy*dy+dy*dy*dz*dz+dz*dz*dx*dx);
				rho=pointer[0].value[i][j][k];
				rhoVx=pointer[1].value[i][j][k];
				rhoVy=pointer[2].value[i][j][k];
				rhoVz=pointer[3].value[i][j][k];
				Bx=pointer[4].value[i][j][k];
				By=pointer[5].value[i][j][k];
				Bz=pointer[6].value[i][j][k];
				pressure=pressure_obj.value[i][j][k];
				Temp_dt=dxyz/( sqrt( pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2) )/rho+ \
					sqrt( ( Bx*Bx+By*By+Bz*Bz+phy_gamma*pressure )/rho ) );
				if (Temp_dt < test_dt)
				{
					test_dt=Temp_dt;
				}
			}
		}
	}
	dt=0.5*min(dtmin,test_dt);
	//cout<<"Set_dt invoked! And dt="<<dt<<endl;
	return dt;
}

// Step on variables
void step_on(VARIABLE *pointer, BASIC_VARIABLE flux[][3], double time, double time_interv)
{
//	VARIABLE var_intmedit[8], current[3];
//	BASIC_VARIABLE pressure, eta, Temp_flux[8][3];            // intermediate variable.    Changed to following declaration! 
/*In order to save space of stack, use new to allocate memory on heap: start */
	VARIABLE *var_intmedit=new VARIABLE[8];
	VARIABLE *current=new VARIABLE[3];
	BASIC_VARIABLE *pre_eta_pointer=new BASIC_VARIABLE[2];
	BASIC_VARIABLE (*Temp_flux)[3]=new BASIC_VARIABLE[8][3];
	BASIC_VARIABLE &pressure=pre_eta_pointer[0];
	BASIC_VARIABLE &eta=pre_eta_pointer[1];
/*In order to save space of stack, use new to allocate memory on heap: end*/
/*Updating variables: start*/
	exclude_soucrce_half_update(var_intmedit, flux, time_interv, First);       
// Last statement do a forward differentiating, so that the value at boundary has not been updated.
//   And current calculating must be done excluding boundary value. The same with pressure, eta, flux.
	cal_current(current, var_intmedit, Incomplete);       
	set_eta(eta, var_intmedit, current, time, Incomplete);
	cal_pressure(pressure, var_intmedit, Incomplete);
	cal_flux(Temp_flux, var_intmedit, current, pressure, eta, Incomplete);
	exclude_soucrce_half_update(pointer, Temp_flux, time_interv, Second);
	//source_update(pointer, time_interv);
/*Updating variables: end*/

/*Set bundary: start*/
	pointer[0].boundary_set(Positive,Positive);
	pointer[1].boundary_set(Negative,Positive);pointer[2].boundary_set(Negative,Positive);pointer[3].boundary_set(Positive,Negative);
	pointer[4].boundary_set(Positive,Negative);pointer[5].boundary_set(Positive,Negative);pointer[6].boundary_set(Negative,Positive);
	pointer[7].boundary_set(Positive,Positive);
/*Set bundary: end*/
/*delete dynamic memory: start */
	delete []var_intmedit;
	delete []current;
	delete []pre_eta_pointer;
	delete []Temp_flux;
/*delete dynamic memory: end*/
	//cout<<"Step_on Step_on!"<<endl;
}

void smooth(VARIABLE *pointer, double time, int nstep)      // how many times do smooth
{
//	VARIABLE tempvar[8];   changed to dynamic array as followings
	VARIABLE *tempvar=new VARIABLE[8];
//	double caf0=1.-0.25*tanh(time/100.);        // caf0 could be changed to caf
	int i,j,k, n;
	for (n=0;n<8;n++)
	{
		for(i=0;i<Grid_Num_x;i++)
		{
			for(j=0;j<Grid_Num_y;j++)
			{
				for(k=0;k<Grid_Num_z;k++)
				{
					tempvar[n].value[i][j][k]=pointer[n].value[i][j][k]-var_x[n][i];
				}
			}
		}
	}
	if (nstep%3==0)
	{
		for (n=0;n<8;n++)
			tempvar[n].smooth_xyz(1);     // smooth_xyz(times)

		tempvar[0].boundary_set(Positive,Positive);
		tempvar[1].boundary_set(Negative,Positive);tempvar[2].boundary_set(Negative,Positive);tempvar[3].boundary_set(Positive,Negative);
		tempvar[4].boundary_set(Positive,Negative);tempvar[5].boundary_set(Positive,Negative);tempvar[6].boundary_set(Negative,Positive);
		tempvar[7].boundary_set(Positive,Positive);

		for (n=0;n<4;n++)
			tempvar[n].average(0.996);                       // avrg1
		tempvar[7].average(0.996);

		tempvar[0].boundary_set(Positive,Positive);
		tempvar[1].boundary_set(Negative,Positive);tempvar[2].boundary_set(Negative,Positive);tempvar[3].boundary_set(Positive,Negative);
		tempvar[4].boundary_set(Positive,Negative);tempvar[5].boundary_set(Positive,Negative);tempvar[6].boundary_set(Negative,Positive);
		tempvar[7].boundary_set(Positive,Positive);
		for (n=4;n<7;n++)
		{
			tempvar[n].average(0.996);                       // avrg2
		}
		tempvar[0].boundary_set(Positive,Positive);
		tempvar[1].boundary_set(Negative,Positive);tempvar[2].boundary_set(Negative,Positive);tempvar[3].boundary_set(Positive,Negative);
		tempvar[4].boundary_set(Positive,Negative);tempvar[5].boundary_set(Positive,Negative);tempvar[6].boundary_set(Negative,Positive);
		tempvar[7].boundary_set(Positive,Positive);
	}
	
	for (n=0;n<8;n++)
	{
		for(i=0;i<Grid_Num_x;i++)
		{
			for(j=0;j<Grid_Num_y;j++)
			{
				for(k=0;k<Grid_Num_z;k++)
				{
					sub_var[n][i][j][k]=tempvar[n].value[i][j][k];
					pointer[n].value[i][j][k]=tempvar[n].value[i][j][k]+var_x[n][i];
				}
			}
		}
	}
// deleting dynamic memory space
	delete []tempvar;
}
