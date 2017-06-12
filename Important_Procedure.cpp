// Important_Procedure.cpp
// 2017/5/22 : 2D output is added
// 2017/5/25 : In order to spare stack space, change arguments of function from class type to class referrence type, for example, in 'cal_flux' and 'set_dt' and 'step_on'

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
#include "Basic_Parameter.h"
#include "Physics_Parameter.h"
#include "Variables_Definition.h"
#include "Procedure.h"

// Clarification of mesh-grid
extern double X[], Y[], Z[], X_interval[], Y_interval[], Z_interval[];
extern int nstep;
double sub_var[8][Grid_Num_x][Grid_Num_y][Grid_Num_z];  //subsidiary_variable
double var_x[8][Grid_Num_x], var_x_plushalfdx[8][Grid_Num_x];          // fx, fxi, used in 2-Lax_Wendroff stepon
// For following procedures
ofstream max_dt_out("max_dt.txt");
ofstream min_dt_out("min_dt.txt");

void set_mesh()
{
	int nx,ny,nz;
	int out_Grid_x, out_Grid_y, out_Grid_z;
	ofstream out("grid.dat");
	X[0]=x_min;	Y[0]=y_min;	Z[0]=z_min;
	nx=Grid_Num_x-1;	ny=Grid_Num_y-1;	nz=Grid_Num_z-1;
	out_Grid_x=(Grid_Num_x-1)/num_out;	out_Grid_y=(Grid_Num_y-1)/num_out;	out_Grid_z=(Grid_Num_z-1)/num_out;

	int i,j,k;
	if (uniform_x==True)
	{
		for (i=0;i<nx;i++)
		{			
			X_interval[i]=(x_max-x_min)/nx;	
			X[i+1]=X[i]+X_interval[i];
		}
	}
	else
	{
		for (i=0;i<nx;i++)
		{
			X_interval[i]=0.01;
			X[i+1]=X[i]+X_interval[i];
		}
	}
	X_interval[Grid_Num_x-1]=X_interval[Grid_Num_x-2];

	if (uniform_y==True)
	{
		for (j=0;j<ny;j++)
		{			
			Y_interval[j]=(y_max-y_min)/ny;	
			Y[j+1]=Y[j]+Y_interval[j];
		}
	}
	else
	{
		for (j=0;j<ny;j++)
		{			
			Y_interval[j]=0.01;
			Y[j+1]=Y[j]+Y_interval[j];
		}
	}
	Y_interval[Grid_Num_y-1]=Y_interval[Grid_Num_y-2];
	
	if (uniform_z==True)
	{
		for (k=0;k<nz;k++)
		{			
			Z_interval[k]=(z_max-z_min)/nz;	
			Z[k+1]=Z[k]+Z_interval[k];
		}
	}
	else
	{
		for (k=0;k<nz;k++)
		{			
			Z_interval[k]=0.01;
			Z[k+1]=Z[k]+Z_interval[k];
		}
	}
	Z_interval[Grid_Num_z-1]=Z_interval[Grid_Num_z-2];
	
	for (i=0;i<=num_out;i++)
		out<<" "<<X[i*out_Grid_x];
	out<<endl;
	for (j=0;j<=num_out;j++)
		out<<" "<<Y[j*out_Grid_y];
	out<<endl;
	for (k=0;k<=num_out;k++)
		out<<" "<<Z[k*out_Grid_z];
	out<<endl;

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
	double x,y,z, dx,dy,dz, r5;
	double rho, Bx, By, Bz, rhoVx, rhoVy, rhoVz;
	double B_Energy, V_Energy, pressure;
// Following used are global variables listed in "Basic_Parameter.h"
	rm=rho_m_0;	rs=rho_s_0;
	bm=B_m_0; bs=B_s_0; 
	v0=v_0;
	betam=beta_m; 
	pressuremtotal=betam*0.5*pow(bm,2)+0.5*pow(bm,2);      //betam*0.5*pow(bm,2)=pmsp
	
	VARIABLE *current=new VARIABLE[3];
	VARIABLE *sub_mag_field=new VARIABLE[3];               //sub_ for subsidary

	int i,j,k,n;
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
				rhoVx=0.;
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
				rhoVx=0.;
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
			}
		}
	}
	cal_current(current, pointer);
	for (i=0;i<Grid_Num_x;i++)
		for (j=0;j<Grid_Num_y;j++)
			for (k=0;k<Grid_Num_z;k++)
			{
				pointer[2].value[i][j][k]=pointer[2].value[i][j][k]+vyi_0*di*current[1].value[i][j][k];	
				sub_mag_field[0].value[i][j][k]=sub_var[4][i][j][k];
				sub_mag_field[1].value[i][j][k]=sub_var[5][i][j][k];
				sub_mag_field[2].value[i][j][k]=sub_var[6][i][j][k];
			}
	
	cal_current(current,sub_mag_field-4);         // -4 is because that sub_mag_field[i]=pointer[i+4];
	for (i=0;i<Grid_Num_x;i++)
		for (j=0;j<Grid_Num_y;j++)
			for (k=0;k<Grid_Num_z;k++)
				sub_var[2][i][j][k]=sub_var[2][i][j][k]+vyi_0*di*current[1].value[i][j][k];	

	for (n=0;n<8;n++)
		for (i=0;i<Grid_Num_x;i++)
		{
			var_x[n][i]=pointer[n].value[i][0][0];
			var_x_plushalfdx[n][i]=sub_var[n][i][0][0];
		}

	for (n=0;n<8;n++)
		for (i=0;i<Grid_Num_x;i++)
			for (j=0;j<Grid_Num_y;j++)
				for (k=0;k<Grid_Num_z;k++)
					sub_var[n][i][j][k]=pointer[n].value[i][j][k]-var_x[n][i];
	// delete dynamic variable-array
	delete []current;
	delete []sub_mag_field;
	//cout<<"Initialize invoked! But I really don't know the setup written by Teacher Ma! Waiting to be changed to a symmetric Harris Current Sheet!"<<endl;
}

void harris_current_initia(VARIABLE *pointer, BASIC_VARIABLE &pressure_obj)
{	
	double rhoinfinity, Bal_coeff, norm_lambda;
	double x,y,z, dx,dy,dz;
	double rho, Bx, By, Bz, rhoVx, rhoVy, rhoVz;
	double B_Energy, V_Energy, pressure;
// Following used are global variables listed in "Basic_Parameter.h"
	rhoinfinity=rho_infinity;
	Bal_coeff=Balance_coefficient;
	norm_lambda=normalized_lambda;
	
	VARIABLE *current=new VARIABLE[3];
	VARIABLE *sub_mag_field=new VARIABLE[3];               //sub_ for subsidary

	int i,j,k,n;
	for (i=0;i<Grid_Num_x;i++)
	{
		for (j=0;j<Grid_Num_y;j++)
		{
			for (k=0;k<Grid_Num_z;k++)
			{
				x=X[i]; dx=X_interval[i];
				y=Y[j]; dy=Y_interval[j];
				z=Z[k]; dz=Z_interval[k];
				rho=pow(1./cosh(x/norm_lambda),2)+rhoinfinity;
				Bz=tanh(x/norm_lambda);
				By=0;
				Bx=0;
				pointer[0].value[i][j][k]=rho;
				pointer[4].value[i][j][k]=Bx;
				pointer[5].value[i][j][k]=By;
				pointer[6].value[i][j][k]=Bz;				
				rhoVx=0.;
				rhoVy=0;       // rho*Bal_coeff*2./norm_lambda should be present in PIC simulation, not in fluid  
							   // when fluid equation does not take collision into consideration
				rhoVz=0;
				pointer[1].value[i][j][k]=rhoVx;
				pointer[2].value[i][j][k]=rhoVy;
				pointer[3].value[i][j][k]=rhoVz;				
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				pressure=Bal_coeff*rho;
				pointer[7].value[i][j][k]=B_Energy+V_Energy+pressure/(phy_gamma-1);
				pressure_obj.value[i][j][k]=pressure;
				/*plus half dx: start*/
				rho=pow(1./cosh((x+dx/2.)/norm_lambda),2)+rhoinfinity;
				Bz=tanh((x+dx/2.)/norm_lambda);
				By=0;
				Bx=0;
				sub_var[0][i][j][k]=rho;
				sub_var[4][i][j][k]=Bx;
				sub_var[5][i][j][k]=By;
				sub_var[6][i][j][k]=Bz;				
				rhoVx=0.;
				rhoVy=0;       // rho*Bal_coeff*2./norm_lambda should be present in PIC simulation, not in fluid  
							   // when fluid equation does not take collision into consideration
				rhoVz=0;
				sub_var[1][i][j][k]=rhoVx;
				sub_var[2][i][j][k]=rhoVy;
				sub_var[3][i][j][k]=rhoVz;				
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				pressure=Bal_coeff*rho;
				sub_var[7][i][j][k]=B_Energy+V_Energy+pressure/(phy_gamma-1);
				/*plus half dx: end*/
			}
		}
	}
	
	for (n=0;n<8;n++)
		for (i=0;i<Grid_Num_x;i++)
		{
			var_x[n][i]=pointer[n].value[i][0][0];
			var_x_plushalfdx[n][i]=sub_var[n][i][0][0];
		}

	for (n=0;n<8;n++)
		for (i=0;i<Grid_Num_x;i++)
			for (j=0;j<Grid_Num_y;j++)
				for (k=0;k<Grid_Num_z;k++)
					sub_var[n][i][j][k]=pointer[n].value[i][j][k]-var_x[n][i];
	// delete dynamic variable-array
	delete []current;
	delete []sub_mag_field;
	//cout<<"Initialize invoked! But I really don't know the setup written by Teacher Ma! Waiting to be changed to a symmetric Harris Current Sheet!"<<endl;
}

// Extracting electric field from flux of magnetic
void ext_from_var( BASIC_VARIABLE *Electric_field, BASIC_VARIABLE *var, VARIABLE *current, BASIC_VARIABLE &eta_obj)
{
	double rho, Vx, Vy, Vz;
	double Bx, By, Bz;
	double eta;
	double current_x, current_y, current_z;
	double V_cross_B_x, V_cross_B_y, V_cross_B_z;
	// Assignment
	int i,j,k;
	for (i=0;i<Grid_Num_x;i++)
	{
		for (j=0;j<Grid_Num_y;j++)
		{
			for (k=0;k<Grid_Num_z;k++)
			{
				rho=var[0].value[i][j][k];
				Vx=var[1].value[i][j][k]/rho;
				Vy=var[2].value[i][j][k]/rho;
				Vz=var[3].value[i][j][k]/rho;
				Bx=var[4].value[i][j][k];
				By=var[5].value[i][j][k];
				Bz=var[6].value[i][j][k];
				current_x=current[0].value[i][j][k];
				current_y=current[1].value[i][j][k];
				current_z=current[2].value[i][j][k];
				eta=eta_obj.value[i][j][k];

				//Magnetic Induction Eq.
				V_cross_B_x=(Vy-di/rho*current_y)*Bz-(Vz-di/rho*current_z)*By;                     // don't take di into acount, that is di=0
				V_cross_B_y=(Vz-di/rho*current_z)*Bx-(Vx-di/rho*current_x)*Bz;
				V_cross_B_z=(Vx-di/rho*current_x)*By-(Vy-di/rho*current_y)*Bx;

				Electric_field[2].value[i][j][k]=-V_cross_B_z+eta*current_z;
				Electric_field[0].value[i][j][k]=-V_cross_B_x+eta*current_x;
				Electric_field[1].value[i][j][k]=-V_cross_B_y+eta*current_y;
			}
		}
	}	
}

// add fluctuation to B variable
void add_fluc(VARIABLE *var)
{
	double norm_lambda;
	double x, z, kx, kz, fluc;
	double Bx, Bz, Eng, fluc_Bx, fluc_Bz;
	int i, j, k;
	fluc=fluctuation;
	kx=fluc_kx; kz=fluc_kz;
	norm_lambda=normalized_lambda;
	// add fluctuation at z=up and down boundary according to <Hurricane, PoP, 1995>, \delta\psi=fluctuation*cos(k_z*z) 
	if (position_fluc==Boundary)
	{
		for (j=0;j<Grid_Num_y;j++)
			for (k=0; k<Grid_Num_z; k++)
			{
				z=Z[k];
				// at i=0
				Bx=var[4].value[0][j][k]; Eng=var[7].value[0][j][k];
				fluc_Bx=-kz*fluc*sin(kz*z);
				Eng=Eng+1./2*(2*Bx*fluc_Bx+fluc_Bx*fluc_Bx);
				Bx=Bx+fluc_Bx;
				var[4].value[0][j][k]=Bx; var[7].value[0][j][k]=Eng;
				// at i=1
				Bx=var[4].value[1][j][k]; Eng=var[7].value[1][j][k];
				fluc_Bx=-kz*fluc*sin(kz*z);
				Eng=Eng+1./2*(2*Bx*fluc_Bx+fluc_Bx*fluc_Bx);
				Bx=Bx+fluc_Bx;
				var[4].value[1][j][k]=Bx; var[7].value[1][j][k]=Eng;
				// at i=Grid_Num_x-2
				Bx=var[4].value[Grid_Num_x-2][j][k]; Eng=var[7].value[Grid_Num_x-2][j][k];
				fluc_Bx=-kz*fluc*sin(kz*z);
				Eng=Eng+1./2*(2*Bx*fluc_Bx+fluc_Bx*fluc_Bx);
				Bx=Bx+fluc_Bx;
				var[4].value[Grid_Num_x-2][j][k]=Bx; var[7].value[Grid_Num_x-2][j][k]=Eng;
				// at i=Grid_Num_x-1
				Bx=var[4].value[Grid_Num_x-1][j][k]; Eng=var[7].value[Grid_Num_x-1][j][k];
				fluc_Bx=-kz*fluc*sin(kz*z);
				Eng=Eng+1./2*(2*Bx*fluc_Bx+fluc_Bx*fluc_Bx);
				Bx=Bx+fluc_Bx;
				var[4].value[Grid_Num_x-1][j][k]=Bx; var[7].value[Grid_Num_x-1][j][k]=Eng;
			    // Does it need to add an fluctuation on sub_var[[][][]?????
			}
	}
	// add fluctuation at neutral-line according to <karimabadi, JGR, 2004>, 
	else if (position_fluc==Neutral_Line)
	{
		for (i=0; i<Grid_Num_x; i++)
			for (j=0;j<Grid_Num_y;j++)
				for (k=0; k<Grid_Num_z; k++)
				{
					x=X[i]; z=Z[k];
					//if (abs(x)<norm_lambda)
					{
						Bx=var[4].value[i][j][k]; Bz=var[6].value[i][j][k]; Eng=var[7].value[i][j][k];
						fluc_Bx=-kz*fluc*cos(kx*x)*sin(kz*z); fluc_Bz=kx*fluc*sin(kx*x)*cos(kz*z);
						Eng=Eng+1./2*(2*Bx*fluc_Bx+fluc_Bx*fluc_Bx)+1./2*(2*Bz*fluc_Bz+fluc_Bz*fluc_Bz);
						Bx=Bx+fluc_Bx; Bz=Bz+fluc_Bz;
						var[4].value[i][j][k]=Bx; var[6].value[i][j][k]=Bz; var[7].value[i][j][k]=Eng;
					}
					// Does it need to add an fluctuation on sub_var[[][][]?????
				}
	}
	
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
				for (j=0;j<Grid_Num_y;j++)
				{
					current[0].value[i][j][0]=current[0].value[i][j][1];
					current[0].value[i][j][Grid_Num_z-1]=current[0].value[i][j][Grid_Num_z-3];
					current[1].value[i][j][0]=current[1].value[i][j][1];
					current[1].value[i][j][Grid_Num_z-1]=current[1].value[i][j][Grid_Num_z-3];
					current[2].value[i][j][0]=current[2].value[i][j][1];
					current[2].value[i][j][Grid_Num_z-1]=-current[2].value[i][j][Grid_Num_z-3];
				}
			}
		}
		else
		{
			for (n=0;n<3;n++)
			{
				for (i=0;i<Grid_Num_x;i++)
				{
					for (j=0;j<Grid_Num_y;j++)
					{
						current[n].value[i][j][0]=current[n].value[i][j][1];
						current[n].value[i][j][Grid_Num_z-1]=current[n].value[i][j][Grid_Num_z-2];
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
				for (j=0;j<Grid_Num_y-T;j++)
				{
					current[0].value[i][j][0]=current[0].value[i][j][1];
					current[0].value[i][j][Grid_Num_z-1-T]=current[0].value[i][j][Grid_Num_z-2-T];
					current[1].value[i][j][0]=current[1].value[i][j][1];
					current[1].value[i][j][Grid_Num_z-1-T]=current[1].value[i][j][Grid_Num_z-2-T];
					current[2].value[i][j][0]=current[2].value[i][j][1];
					current[2].value[i][j][Grid_Num_z-1-T]=-current[2].value[i][j][Grid_Num_z-2-T];
				}
			}
		}
		else
		{
			for (n=0;n<3;n++)
			{
				for (i=0;i<Grid_Num_x-T;i++)
				{
					for (j=0;j<Grid_Num_y-T;j++)
					{
						current[n].value[i][j][0]=current[n].value[i][j][1];
						current[n].value[i][j][Grid_Num_z-1-T]=current[n].value[i][j][Grid_Num_z-2-T];
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
	if (T==Complete)
		make_pressure_positive(pressure_obj, 1e-5);
	//cout<<"Pressure is calculated from various kinds of energy, and doesn't explicitly step on according to equation!"<<endl;
}

// Homogeneous eta, waiting to be amended.
void set_eta(BASIC_VARIABLE &eta_obj, VARIABLE *pointer, VARIABLE *current, double time, Type T, Type TT)   // Type T is for Complete or Incomplete setting; Type TT is for uniform or non-uniform setting
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
					if (TT==Uniform)
						eta_obj.value[i][j][k]=magnetic_Renolds_Number;
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
					if (TT==Uniform)
						eta_obj.value[i][j][k]=magnetic_Renolds_Number;
				}
			}
		}
	}	
	//cout<<"Set_eta invoked! But for now eta is simply setted to 0.1! And I'm sorry to inform you that you should change it!"<<endl;
}

// accumulate energy from V, B and P
void accumulate_eng(double (*Eng_array)[Grid_Num_y][Grid_Num_z], VARIABLE *pointer, BASIC_VARIABLE &pressure_obj, Type T)
{
	double rho, Bx, By, Bz, rhoVx, rhoVy, rhoVz;
	double B_Energy, V_Energy, pressure;

	int i,j,k;
	for (i=0;i<Grid_Num_x-T;i++)
		for (j=0;j<Grid_Num_y-T;j++)
			for (k=0;k<Grid_Num_z-T;k++)
			{
				rho=pointer[0].value[i][j][k];
				rhoVx=pointer[1].value[i][j][k];
				rhoVy=pointer[2].value[i][j][k];
				rhoVz=pointer[3].value[i][j][k];
				Bx=pointer[4].value[i][j][k];
				By=pointer[5].value[i][j][k];
				Bz=pointer[6].value[i][j][k];
				pressure=pressure_obj.value[i][j][k];
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				Eng_array[i][j][k]=B_Energy+V_Energy+pressure/(phy_gamma-1);
			}
}

// Calculating flux from variables
void cal_7flux(BASIC_VARIABLE flux[][3], VARIABLE *pointer, VARIABLE *current,\
	BASIC_VARIABLE &pressure_obj, BASIC_VARIABLE &eta_obj, Type T)
{
	double rho, Vx, Vy, Vz;
	double Bx, By, Bz;
	double B_Energy;  // V_Energy is not used !
	double pressure, eta;
	double current_x, current_y, current_z;
	double V_cross_B_x, V_cross_B_y, V_cross_B_z;
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
				pressure=pressure_obj.value[i][j][k];
				current_x=current[0].value[i][j][k];
				current_y=current[1].value[i][j][k];
				current_z=current[2].value[i][j][k];
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
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
				V_cross_B_x=(Vy-di/rho*current_y)*Bz-(Vz-di/rho*current_z)*By;                     // don't take di into acount, that is di=0
				V_cross_B_y=(Vz-di/rho*current_z)*Bx-(Vx-di/rho*current_x)*Bz;
				V_cross_B_z=(Vx-di/rho*current_x)*By-(Vy-di/rho*current_y)*Bx;

				flux[4][0].value[i][j][k]=0;
				flux[4][1].value[i][j][k]=-V_cross_B_z+eta*current_z;
				flux[4][2].value[i][j][k]=V_cross_B_y-eta*current_y;

				flux[5][0].value[i][j][k]=V_cross_B_z-eta*current_z;
				flux[5][1].value[i][j][k]=0;
				flux[5][2].value[i][j][k]=-V_cross_B_x+eta*current_x;

				flux[6][0].value[i][j][k]=-V_cross_B_y+eta*current_y;
				flux[6][1].value[i][j][k]=V_cross_B_x-eta*current_x;
				flux[6][2].value[i][j][k]=0;
			}
		}
	}
	//cout<<"I'm Calculating flux!"<<endl;
}

// calculate engergy flux
void cal_eng_flux(BASIC_VARIABLE *eng_flux, VARIABLE *pointer, VARIABLE *current,\
	BASIC_VARIABLE &pressure_obj, BASIC_VARIABLE &eta_obj, Type T)
{
	double rho, Vx, Vy, Vz;
	double Bx, By, Bz;
	double Energy, B_Energy, V_Energy;  // V_Energy is not used !	
	double pressure, eta;
	double current_x, current_y, current_z;
	double V_dot_B;

	double (*Eng_array)[Grid_Num_y][Grid_Num_z]=new double[Grid_Num_x][Grid_Num_y][Grid_Num_z];
	accumulate_eng(Eng_array, pointer, pressure_obj, T);    // accumulate energy from new V, new B and old P
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
				Energy=Eng_array[i][j][k];
				pressure=pressure_obj.value[i][j][k];
				current_x=current[0].value[i][j][k];
				current_y=current[1].value[i][j][k];
				current_z=current[2].value[i][j][k];
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*rho*(pow(Vx,2)+pow(Vy,2)+pow(Vz,2));
				eta=eta_obj.value[i][j][k];

				// Energy_fulx
				V_dot_B=Vx*Bx+Vy*By+Vz*Bz;
				Energy=Energy+pressure+B_Energy;       // Energy flux brought by V
				eng_flux[0].value[i][j][k]=Vx*Energy-Bx*V_dot_B+eta* \
					(current_y*Bz-current_z*By);
				eng_flux[1].value[i][j][k]=Vy*Energy-By*V_dot_B+eta* \
					(current_z*Bx-current_x*Bz);
				eng_flux[2].value[i][j][k]=Vz*Energy-Bz*V_dot_B+eta* \
					(current_x*By-current_y*Bx);
			}
		}
	}
	delete []Eng_array;
	//cout<<"I'm Calculating flux!"<<endl;
}

// Setting time-interval
double set_dt(VARIABLE *pointer, BASIC_VARIABLE &eta_obj, VARIABLE *current, BASIC_VARIABLE &pressure_obj, double time, double last_dt)
{
	//set_eta(eta_obj, pointer, current, time);
	double dt;
	double dx,dy,dz,dxyz;	
	double Temp_dt, dt_min=1000., test_dt=1000.;
	double rho, rhoVx, rhoVy, rhoVz, Bx, By, Bz, pressure;
	int max_dt_i, max_dt_j, max_dt_k;                                       // for diagnostic
	double max_rho, max_Vx, max_Vy, max_Vz, max_Bx, max_By, max_Bz, max_P;  // for diagnostic

	int i,j,k, times=0;
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
				if (Temp_dt > dt_min)
				{
					if (times==0)
					{
						max_dt_out<<endl;
						max_dt_out<<" When dt >= "<<int(dt_min)<<" , time step is "<<nstep<<endl\
							      <<"where location is ( i = "<<setw(3)<<i<<", j = "<<setw(3)<<j<<", k = "<<setw(3)<<k<<" )"<<endl;
						max_dt_out<<setiosflags(ios::scientific)<<setprecision(3);
						times+=1;
					}
					else
						max_dt_out<<"And               (     "<<setw(3)<<i<<",     "<<setw(3)<<j<<",     "<<setw(3)<<k<<" )"<<endl;
					max_dt_out<<"     rho          Vx          Vy          Vz          Bx     "\
						      <<"     By          Bz          P"<<endl;
					max_dt_out<<setw(13)<<rho<<setw(12)<<rhoVx/rho<<setw(12)<<rhoVy/rho<<setw(12)<<rhoVz/rho\
						<<setw(12)<<Bx<<setw(12)<<By<<setw(12)<<Bz<<setw(11)<<pressure<<endl;
				}
				if (Temp_dt < test_dt)
				{
					test_dt=Temp_dt;
					max_dt_i=i; max_dt_j=j; max_dt_k=k;
					max_rho=rho; max_Vx=rhoVx/rho; max_Vy=rhoVy/rho; max_Vz=rhoVz/rho;
					max_Bx=Bx; max_By=By; max_Bz=Bz; max_P=pressure;
				}
			}
		}
	}

	dt=0.5*min(dt_min,test_dt);

	if (time>0 && dt<.5*last_dt)
	{
		min_dt_out<<setiosflags(ios::scientific)<<setprecision(3);
		min_dt_out<<"  Noth that dt changes a lot. This dt is calculated on "\
		          <<"location ( i = "<<setw(3)<<max_dt_i<<", j = "<<setw(3)<<max_dt_j<<", k = "<<setw(3)<<max_dt_k<<" )"<<endl;
		min_dt_out<<"And time step is "<<nstep<<" , with variables are:"<<endl;
		min_dt_out<<"     rho          Vx          Vy          Vz          Bx     "\
				  <<"     By          Bz          P"<<endl;
		min_dt_out<<setw(13)<<max_rho<<setw(12)<<max_Vx<<setw(12)<<max_Vy<<setw(12)<<max_Vz\
			<<setw(12)<<max_Bx<<setw(12)<<max_By<<setw(12)<<max_Bz<<setw(11)<<max_P<<endl<<endl;
	}
	//cout<<"Set_dt invoked! And dt="<<dt<<endl;
	return dt;
}

// Step on variables
void step_on(VARIABLE *pointer, VARIABLE *current, BASIC_VARIABLE &pressure, BASIC_VARIABLE &eta, double time, double time_interv)
{
/*In order to save space of stack, use new to allocate memory on heap: start */
	VARIABLE *var_intmedit=new VARIABLE[8];
	VARIABLE *temp_current=new VARIABLE[3];
	BASIC_VARIABLE *temp_pre_eta_pointer=new BASIC_VARIABLE[2];
	BASIC_VARIABLE (*flux7)[3]=new BASIC_VARIABLE[7][3];
	BASIC_VARIABLE *eng_flux=new BASIC_VARIABLE[3];
	BASIC_VARIABLE &temp_pressure=temp_pre_eta_pointer[0];
	BASIC_VARIABLE &temp_eta=temp_pre_eta_pointer[1];
/*In order to save space of stack, use new to allocate memory on heap: end*/
/*Updating variables: start*/
	cal_7flux(flux7, pointer, current, pressure, eta);
	exclude_soucrce_half_7update(var_intmedit, flux7, time_interv, First);

	var_intmedit[0].right_boundary_set(Positive,Positive);
	var_intmedit[1].right_boundary_set(Negative,Positive);var_intmedit[2].right_boundary_set(Negative,Positive);var_intmedit[3].right_boundary_set(Positive,Negative);
	var_intmedit[4].right_boundary_set(Positive,Negative);var_intmedit[5].right_boundary_set(Positive,Negative);var_intmedit[6].right_boundary_set(Negative,Positive);

	cal_current(temp_current, var_intmedit); 
	set_eta(temp_eta, var_intmedit, temp_current, time);

	cal_eng_flux(eng_flux, var_intmedit, temp_current, pressure, temp_eta);     // ues old P and new V, new B, new current and new eta to accumulate energy flux
	exclude_source_hlaf_update_eng(var_intmedit[7], eng_flux, time_interv, First);  // incomplete time-backward-difference concerning about V and B

	var_intmedit[7].right_boundary_set(Positive,Positive);                  // there is no need

// Last statement do a forward differentiating, so that the value at boundary has not been updated.
//   And current calculating must be done excluding boundary value. The same with pressure, eta, flux. 
	cal_current(temp_current, var_intmedit, Incomplete); 
	set_eta(temp_eta, var_intmedit, temp_current, time, Incomplete);
	cal_pressure(temp_pressure, var_intmedit, Incomplete);

	cal_7flux(flux7, var_intmedit, temp_current, temp_pressure, temp_eta, Incomplete);
	exclude_soucrce_half_7update(pointer, flux7, time_interv, Second);

	pointer[0].left_boundary_set();
	pointer[1].left_boundary_set();pointer[2].left_boundary_set();pointer[3].left_boundary_set();
	pointer[4].left_boundary_set();pointer[5].left_boundary_set();pointer[6].left_boundary_set();

	cal_current(temp_current, pointer, Incomplete); 
	set_eta(temp_eta, pointer, temp_current, time, Incomplete); 

	cal_eng_flux(eng_flux, pointer, temp_current, temp_pressure, temp_eta, Incomplete );     // ues old P and new V, new B, new current and new eta to accumulate energy flux
	exclude_source_hlaf_update_eng(pointer[7], eng_flux, time_interv, Second);              // incomplete time-backward-difference concerning about V and B

	pointer[7].left_boundary_set();
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
	delete []temp_current;
	delete []temp_pre_eta_pointer;
	delete []flux7;
	delete []eng_flux;
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
