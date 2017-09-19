// Important_Procedure.cpp
// 2017/5/22 : 2D output is added
// 2017/5/25 : In order to spare stack space, change arguments of function from class type to class referrence type, for example, in 'cal_flux' and 'set_dt' and 'step_on'
// 2017/7/16 : Add iosthermal process when initial energy with pressure and calculating pressure

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
#include "Basic_Parameter.h"
#include "Physics_Parameter.h"
#include "Variables_Definition.h"
#include "Procedure.h"
#ifndef min
#define min(a,b) (a<b ? a : b)
#endif

// Clarification of mesh-grid
extern double X[], Y[], Z[], X_interval[], Y_interval[], Z_interval[];
extern int nstep;
double sub_var[8][Grid_Num_x][Grid_Num_y][Grid_Num_z];  //subsidiary_variable
double var_x[8][Grid_Num_x], var_x_plushalfdx[8][Grid_Num_x];          // fx, fxi, used in 2-Lax_Wendroff stepon
// For following procedures
ofstream max_dt_out("max_dt.txt");
ofstream min_dt_out("min_dt.txt");
ofstream pre_out("pressure_is_negative.txt");

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
	
	out<<setiosflags(ios::fixed)<<setprecision(6);
	for (i=0;i<=num_out;i++)
		out<<setw(10)<<X[i*out_Grid_x]<<" ";
	out<<endl;
	if (out_Grid_y !=0)
		for (j=0;j<=num_out;j++)
			out<<setw(10)<<Y[j*out_Grid_y]<<" ";
	out<<endl;
	for (k=0;k<=num_out;k++)
		out<<setw(10)<<Z[k*out_Grid_z]<<" ";
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
				pointer[7].value[i][j][k]=B_Energy+V_Energy+pressure/(phy_gamma - 1);
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
				sub_var[7][i][j][k]=B_Energy+V_Energy+pressure/(phy_gamma - 1);
				/*plus half dx: end*/
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
				if (phy_gamma - 1 < iso_therm_coeff)
					pointer[7].value[i][j][k] = B_Energy + V_Energy;
				else
					pointer[7].value[i][j][k]=B_Energy+V_Energy+pressure/(phy_gamma - 1);
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
				if (phy_gamma - 1 < iso_therm_coeff)
					sub_var[7][i][j][k] = B_Energy + V_Energy;
				else
					sub_var[7][i][j][k] = B_Energy + V_Energy + pressure / (phy_gamma - 1);
				/*plus half dx: end*/
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

void shear_flow_harris_initia(VARIABLE *pointer, BASIC_VARIABLE &pressure_obj)
{
	double rhoinfinity, Bal_coeff, norm_lambda;
	double x,y,z, dx,dy,dz;
	double rho, Bx, By, Bz, rhoVx, rhoVy, rhoVz;
	double B_Energy, V_Energy, pressure;
	double v_tanh, v_tanh_0, vb, ls, xs;

// Following used are global variables listed in "Physics_Parameter.h"
	rhoinfinity=rho_infinity;
	Bal_coeff=Balance_coefficient;
	norm_lambda=normalized_lambda;

// shear concerned, listed in "Physics_Parameter.h"
	vb = velocity_boundary;
	ls = shear_length;
	xs = shear_location;
	v_tanh_0 = 0.5*vb*tanh(-xs/ls);

	int i,j,k,n;
	for (i=0;i<Grid_Num_x;i++)
	{
		for (j=0;j<Grid_Num_y;j++)
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
				v_tanh = 0.5*vb*tanh( (abs(x)-xs)/ls );
				if (x>=0)
					rhoVz = rho*(v_tanh-v_tanh_0);        // initia shear flow
				else
					rhoVz = -rho*(v_tanh-v_tanh_0);
				pointer[1].value[i][j][k]=rhoVx;
				pointer[2].value[i][j][k]=rhoVy;
				pointer[3].value[i][j][k]=rhoVz;				
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				pressure=Bal_coeff*rho;
				if (phy_gamma - 1 < iso_therm_coeff)
					pointer[7].value[i][j][k] = B_Energy + V_Energy;
				else
					pointer[7].value[i][j][k] = B_Energy + V_Energy + pressure / (phy_gamma - 1);
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
				v_tanh = 0.5*vb*tanh( (abs(x+dx/2.)-xs)/ls );
				if (x>=0)
					rhoVz = rho*(v_tanh-v_tanh_0);        // initia shear flow
				else
					rhoVz = -rho*(v_tanh-v_tanh_0);
				sub_var[1][i][j][k]=rhoVx;
				sub_var[2][i][j][k]=rhoVy;
				sub_var[3][i][j][k]=rhoVz;				
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				pressure=Bal_coeff*rho;
				if (phy_gamma - 1 < iso_therm_coeff)
					sub_var[7][i][j][k] = B_Energy + V_Energy;
				else
					sub_var[7][i][j][k] = B_Energy + V_Energy + pressure / (phy_gamma - 1);
				/*plus half dx: end*/
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
	
	//cout<<"shear_flow_harris_initia invoked!"<<endl;
}

void write_out(VARIABLE *pointer, int nstop, double time, double time_interval)
{
	int i,j,k, n;
	ofstream timeout("temp_step_to_time.dat");
	ofstream var_out[8];
	ofstream sub_var_out("temp_sub_var.dat");                  // This and the following two is very important for calculating 
	ofstream var_x_out("temp_var_x_out.dat");
	ofstream var_x_plushalfdx_out("temp_var_x_plushalfdx.dat");
	open_var_files(var_out);
	timeout<<setw(6)<<nstop<<setw(6)<<" ";
	timeout<<setiosflags(ios::scientific)<<setprecision(16);
	timeout<<setw(25)<<time<<setw(25)<<time_interval<<endl;
	sub_var_out<<setiosflags(ios::scientific)<<setprecision(16);
	var_x_out<<setiosflags(ios::scientific)<<setprecision(16);
	var_x_plushalfdx_out<<setiosflags(ios::scientific)<<setprecision(16);
	for (n=0;n<8;n++)
	{
		pointer[n].record(var_out[n]);
		for (i=0;i<Grid_Num_x;i++)
		{
			var_x_out<<setw(25)<<var_x[n][i]<<" ";
			var_x_plushalfdx_out<<setw(25)<<var_x_plushalfdx[n][i]<<" ";
			for (j=0;j<Grid_Num_y;j++)
			{
				for (k=0;k<Grid_Num_z;k++)
					sub_var_out<<setw(25)<<sub_var[n][i][j][k]<<" ";
				sub_var_out<<endl;
			}
			var_x_out<<endl;
			var_x_plushalfdx_out<<endl;
			sub_var_out<<endl;
		}
		var_x_out<<endl;
		var_x_plushalfdx_out<<endl;
		sub_var_out<<endl;
	}
	timeout.close();
	close_var_files(var_out);
	sub_var_out.close();
	var_x_out.close();
	var_x_plushalfdx_out.close();
}

void read_in(VARIABLE *pointer, int &nstart, double &time, double &time_interval)
{
	int i,j,k, n;
	ifstream timein("temp_step_to_time.dat");
	ifstream var_in[8];
	ifstream sub_var_in("temp_sub_var.dat");                  // This and the following two is very important for calculating 
	ifstream var_x_in("temp_var_x_out.dat");
	ifstream var_x_plushalfdx_in("temp_var_x_plushalfdx.dat");
	open_var_files(var_in);
	timein>>nstart>>time>>time_interval;
	for (n=0; n<8; n++)
	{
		pointer[n].fill(var_in[n]);
		for (i=0;i<Grid_Num_x;i++)
		{
			var_x_in>>var_x[n][i];
			var_x_plushalfdx_in>>var_x_plushalfdx[n][i];
			for (j=0;j<Grid_Num_y;j++)
			{
				for (k=0;k<Grid_Num_z;k++)
					sub_var_in>>sub_var[n][i][j][k];
			}
		}			
	}
	timein.close();
	close_var_files(var_in);
	sub_var_in.close();
	var_x_in.close();
	var_x_plushalfdx_in.close();
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
	// (can write a switch statement to include different extrapolation method
	// x-linear extrapolation yz-equivalent
	// boundary at i=0,Grid_Num_x-1-T
	if(T==Complete)                        // linear extrapolation
	{
		// boundary at i=0,Grid_Num_x-1
		if (half_x==True)
		{
			for (j=1;j<Grid_Num_y-1;j++)
				for(k=1;k<Grid_Num_z-1;k++)
				{
					current[0].value[0][j][k]=2*current[0].value[1][j][k]-current[0].value[2][j][k];
					current[0].value[Grid_Num_x-1][j][k]=current[0].value[Grid_Num_x-3][Grid_Num_y-1-j][k];           // it's meaning?
					current[1].value[0][j][k]=2*current[1].value[1][j][k]-current[1].value[2][j][k];
					current[1].value[Grid_Num_x-1][j][k]=current[1].value[Grid_Num_x-3][Grid_Num_y-1-j][k];
					current[2].value[0][j][k]=2*current[2].value[1][j][k]-current[2].value[2][j][k];
					current[2].value[Grid_Num_x-1][j][k]=-current[2].value[Grid_Num_x-3][Grid_Num_y-1-j][k];					
				}
		}
		else
		{
			for (n=0;n<3;n++)
				for (j=1;j<Grid_Num_y-1;j++)
					for(k=1;k<Grid_Num_z-1;k++)
					{
						current[n].value[0][j][k]=2*current[n].value[1][j][k]-current[n].value[2][j][k];
						current[n].value[Grid_Num_x-1][j][k]=2*current[n].value[Grid_Num_x-2][j][k]- \
									current[n].value[Grid_Num_x-3][j][k];						
					}		
		}	

		// boundary at j=0,Grid_Num_y-1
		if (period_y==True)
		{
			for (n=0;n<3;n++)
				for (i=0;i<Grid_Num_x;i++)
					for(k=1;k<Grid_Num_z-1;k++)
					{
						current[n].value[i][0][k]=current[n].value[i][Grid_Num_y-2][k];
						current[n].value[i][Grid_Num_y-1][k]=current[n].value[i][1][k];
					}
		}
		else
		{
			for (n=0;n<3;n++)
				for (i=0;i<Grid_Num_x;i++)
					for(k=1;k<Grid_Num_z-1;k++)
					{
						current[n].value[i][0][k]=current[n].value[i][1][k];
						current[n].value[i][Grid_Num_y-1][k]=current[n].value[i][Grid_Num_y-2][k];
					}
		}

		// boundary at k=0,Grid_Num_z-1
		if (half_z==True)
		{
			for (i=0;i<Grid_Num_x;i++)
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
		else
		{
			if(period_z == True)
			{
				for (n=0;n<3;n++)
					for (i=0;i<Grid_Num_x;i++)
						for (j=0;j<Grid_Num_y;j++)
						{
							current[n].value[i][j][0]=current[n].value[i][j][Grid_Num_z-2];
							current[n].value[i][j][Grid_Num_z-1]=current[n].value[i][j][1];
						}
			}
			else
			{
				for (n=0;n<3;n++)
					for (i=0;i<Grid_Num_x;i++)
						for (j=0;j<Grid_Num_y;j++)
						{
							current[n].value[i][j][0]=current[n].value[i][j][1];
							current[n].value[i][j][Grid_Num_z-1]=current[n].value[i][j][Grid_Num_z-2];
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
		else
		{
			for (n=0;n<3;n++)
				for (j=1;j<Grid_Num_y-1-T;j++)
					for(k=1;k<Grid_Num_z-1-T;k++)
					{
						current[n].value[0][j][k]=current[n].value[1][j][k];
						current[n].value[Grid_Num_x-1-T][j][k]=current[n].value[Grid_Num_x-2-T][j][k];	
					}
		}	

		// boundary at j=0,Grid_Num_y-1-T
		if (period_y==True)
		{
			for (n=0;n<3;n++)
				for (i=0;i<Grid_Num_x-T;i++)
					for(k=1;k<Grid_Num_z-1-T;k++)
					{
						current[n].value[i][0][k]=current[n].value[i][Grid_Num_y-2-T][k];
						current[n].value[i][Grid_Num_y-1-T][k]=current[n].value[i][1][k];
					}
		}
		else
		{
			for (n=0;n<3;n++)
				for (i=0;i<Grid_Num_x-T;i++)
					for(k=1;k<Grid_Num_z-1-T;k++)
					{
						current[n].value[i][0][k]=current[n].value[i][1][k];
						current[n].value[i][Grid_Num_y-1-T][k]=current[n].value[i][Grid_Num_y-2-T][k];
					}
		}

		// boundary at k=0,Grid_Num_z-1-T
		if (half_z==True)
		{
			for (i=0;i<Grid_Num_x-T;i++)
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
		else
		{
			if(period_z == True)
			{
				for (n=0;n<3;n++)
					for (i=0;i<Grid_Num_x-T;i++)
						for (j=0;j<Grid_Num_y-T;j++)
						{
							current[n].value[i][j][0]=current[n].value[i][j][Grid_Num_z-2-T];
							current[n].value[i][j][Grid_Num_z-1-T]=current[n].value[i][j][1];
						}
			}
			else
			{
				for (n=0;n<3;n++)
					for (i=0;i<Grid_Num_x-T;i++)
						for (j=0;j<Grid_Num_y-T;j++)
						{
							current[n].value[i][j][0]=current[n].value[i][j][1];
							current[n].value[i][j][Grid_Num_z-1-T]=current[n].value[i][j][Grid_Num_z-2-T];
						}
			}
		}
	}					
	//cout<<"Current is calculated form Curl B! It's simple, I know!!!"<<endl;	
}

// Homogeneous eta, waiting to be amended.
void set_eta(BASIC_VARIABLE &eta_obj, VARIABLE *pointer, VARIABLE *current, double time, Type T)   // Type T is for Complete or Incomplete setting; Type TT is for uniform or non-uniform setting
{
	double etam[Grid_Num_y];													// eta_main direction, in another wrod eta along guid field, that is y-directin
	double etax, etaz, etal, etab, alpha0=2.0;									// etal for eta_localize; etab for eta_background
	double ylength=15., xtrig=0.0, ztrig=0.0;									// 'ylength' is 'xlen' in Fortran version, which is length of localized region in y-direction; xtrig and ztrig is localized resistivity region's center
	double widthx=width_rho, widthz=2*width_rho;
	double x,y,z, dx,dy,dz;
	int i,j,k;
	etab = Lundquist_Number;
	etal = localized_Resistivity;
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
				if(abs(y)<=ylength)					// ylength is the length of localized resitivity region along y-directiont
					etam[j]=1;                      // This is True
				else
					etam[j]=1.-pow(tanh(abs(y)-ylength),2);
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
				if(abs(y+dy/2.)<=ylength)                    
					etam[j]=1;                      // This is True
				else
					etam[j]=1.-pow(tanh(abs(y+dy/2.)-ylength),2);
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

Logic cal_pressure(BASIC_VARIABLE &pressure_obj, VARIABLE *pointer, Type T)
{	
	double rho,rhoVx, rhoVy, rhoVz, Bx, By, Bz, Eng, pressure;
	double B_Energy, V_Energy;
	double positive_value=1e-5;
	int i,j,k, times=0;
	for (i=0;i<Grid_Num_x-T;i++)
	{
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
				Eng=pointer[7].value[i][j][k];
				B_Energy=0.5*(pow(Bx,2)+pow(By,2)+pow(Bz,2));
				V_Energy=0.5*(pow(rhoVx,2)+pow(rhoVy,2)+pow(rhoVz,2))/rho;
				if (phy_gamma - 1 < iso_therm_coeff)               // iso_therm_coeff is a small value
					pressure = Balance_coefficient * rho; 
				else
					pressure = (phy_gamma - 1)*(Eng - B_Energy - V_Energy);
				pressure_obj.value[i][j][k]=pressure;
				if(pressure<0.)
				{
					if (times==0)
					{
						pre_out<<endl<<"  Oops, pressure is negative, and program can be stopped!!!"<<endl;
						pre_out<<"When Type T = "<<T<<" , and time step is nt = "<<nstep<<endl;
						pre_out<<"Located in ( xi = "<<setw(3)<<i<<", yj = "<<setw(3)<<j<<", zk = "<<setw(3)<<k<<" ) ";
						pre_out<<endl;
						pre_out<<setiosflags(ios::scientific)<<setprecision(3);
						times+=1;
					}
					else
					{
						pre_out<<"           ( xi = "<<setw(3)<<i<<", yj = "<<setw(3)<<j<<", zk = "<<setw(3)<<k<<" ) ";
						pre_out<<endl;
					}
					pre_out<<"     rho          Vx          Vy          Vz          Bx     "\
						      <<"     By          Bz          Eng          P"<<endl;
					pre_out<<setw(13)<<rho<<setw(12)<<rhoVx/rho<<setw(12)<<rhoVy/rho<<setw(12)<<rhoVz/rho\
						<<setw(12)<<Bx<<setw(12)<<By<<setw(12)<<Bz<<setw(13)<<Eng<<setw(11)<<pressure<<endl;
				}
				if(T==Complete && pressure<positive_value)
					pressure_obj.value[i][j][k]=positive_value;
			}
	}
	if (times==0)
		return True;
	else
		return False;
	//cout<<"Pressure is calculated from various kinds of energy, and doesn't explicitly step on according to equation!"<<endl;
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
		for (j=0;j<Grid_Num_y;j++)
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

// Setting time-interval
double set_dt(VARIABLE *pointer, BASIC_VARIABLE &eta_obj, VARIABLE *current, BASIC_VARIABLE &pressure_obj, double time, double last_dt)
{
	//set_eta(eta_obj, pointer, current, time);
	double dt;
	double dx,dy,dz,dxyz;	
	double Temp_dt, dt_min=1000., test_dt=1000.;
	double rho, rhoVx, rhoVy, rhoVz, Bx, By, Bz, Eng, pressure;
	int max_dt_i, max_dt_j, max_dt_k;                                       // for diagnostic
	double max_rho, max_Vx, max_Vy, max_Vz, max_Bx, max_By, max_Bz, max_Eng, max_P;  // for diagnostic

	int i,j,k, times=0;
	for (i=0;i<Grid_Num_x;i++)
	{
		for (j=0;j<Grid_Num_y;j++)
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
				Eng=pointer[7].value[i][j][k];
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
						      <<"     By          Bz          Eng          P"<<endl;
					max_dt_out<<setw(13)<<rho<<setw(12)<<rhoVx/rho<<setw(12)<<rhoVy/rho<<setw(12)<<rhoVz/rho\
						<<setw(12)<<Bx<<setw(12)<<By<<setw(12)<<Bz<<setw(13)<<Eng<<setw(11)<<pressure<<endl;
				}
				if (Temp_dt < test_dt)
				{
					test_dt=Temp_dt;
					max_dt_i=i; max_dt_j=j; max_dt_k=k;
					max_rho=rho; max_Vx=rhoVx/rho; max_Vy=rhoVy/rho; max_Vz=rhoVz/rho;
					max_Bx=Bx; max_By=By; max_Bz=Bz; max_Eng=Eng, max_P=pressure;
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
				  <<"     By          Bz          Eng          P"<<endl;
		min_dt_out<<setw(13)<<max_rho<<setw(12)<<max_Vx<<setw(12)<<max_Vy<<setw(12)<<max_Vz\
			<<setw(12)<<max_Bx<<setw(12)<<max_By<<setw(12)<<max_Bz<<setw(13)<<max_Eng<<setw(11)<<max_P<<endl<<endl;
	}
	//cout<<"Set_dt invoked! And dt="<<dt<<endl;
	return dt;
}

// add fluctuation to B variable
void add_fluc(VARIABLE *var, double time)
{
	double norm_lambda;
	double x, z;
	double Bkx, Bkz, fluc_B, Vkx, Vkz, fluc_V, omega;
	double Bx, Bz, Eng, fluc_Bx, fluc_Bz;
	double rho, rhoVx, rhoVz, fluc_rhoVx, fluc_rhoVz;
	int i, j, k;
	fluc_B=fluctuation_mag_psi;
	fluc_V=fluctuation_velocity_phi;
	Bkx=fluc_B_kx; Bkz=fluc_B_kz;
	Vkx = fluc_V_kx; Vkz = fluc_V_kz;
	omega = fluc_B_kz*velocity_phase;
	norm_lambda=normalized_lambda;
	// add fluctuation at z=up and down boundary according to <Hurricane, PoP, 1995>, \delta\psi=fluctuation*cos(k_z*z) 
	if (position_fluc == x_Boundary || position_fluc == upperBoundary)
	{
		for (j=0;j<Grid_Num_y;j++)
			for (k=0; k<Grid_Num_z; k++)
			{
				z=Z[k];
				// at i=Grid_Num_x-2
				Bx=var[4].value[Grid_Num_x-2][j][k]; Eng=var[7].value[Grid_Num_x-2][j][k];
				fluc_Bx = Bkz*fluc_B*sin(Bkz*z - omega*time);
				Eng=Eng+1./2*(2*Bx*fluc_Bx+fluc_Bx*fluc_Bx);

				rho=var[0].value[Grid_Num_x-2][j][k]; rhoVx=var[1].value[Grid_Num_x-2][j][k];
				fluc_rhoVx=rho*Vkz*fluc_V*sin(Vkz*z - omega*time);
				Eng=Eng+1./2*(2*rhoVx*fluc_rhoVx+fluc_rhoVx*fluc_rhoVx)/rho;

				rhoVx=rhoVx+fluc_rhoVx; var[1].value[Grid_Num_x-2][j][k]=rhoVx;
				Bx=Bx+fluc_Bx;
				var[4].value[Grid_Num_x-2][j][k]=Bx; var[7].value[Grid_Num_x-2][j][k]=Eng;
				// at i=Grid_Num_x-1
				Bx=var[4].value[Grid_Num_x-1][j][k]; Eng=var[7].value[Grid_Num_x-1][j][k];
				fluc_Bx = Bkz*fluc_B*sin(Bkz*z - omega*time);
				Eng=Eng+1./2*(2*Bx*fluc_Bx+fluc_Bx*fluc_Bx);
				
				rho=var[0].value[Grid_Num_x-1][j][k]; rhoVx=var[1].value[Grid_Num_x-1][j][k];
				fluc_rhoVx= rho*Vkz*fluc_V*sin(Vkz*z - omega*time);
				Eng=Eng+1./2*(2*rhoVx*fluc_rhoVx+fluc_rhoVx*fluc_rhoVx)/rho;

				rhoVx=rhoVx+fluc_rhoVx; var[1].value[Grid_Num_x-1][j][k]=rhoVx;
				Bx=Bx+fluc_Bx;
				var[4].value[Grid_Num_x-1][j][k]=Bx; var[7].value[Grid_Num_x-1][j][k]=Eng;
			    // Does it need to add an fluctuation on sub_var[[][][]?????
			}
	}

	if (position_fluc == x_Boundary || position_fluc == downBoundary)
	{
		for (j = 0; j<Grid_Num_y; j++)
			for (k = 0; k<Grid_Num_z; k++)
			{
				z = Z[k];
				// at i=0
				Bx = var[4].value[0][j][k]; Eng = var[7].value[0][j][k];
				fluc_Bx = (-1) * (Bkz*fluc_B*sin(Bkz*z - omega*time));
				Eng = Eng + 1. / 2 * (2 * Bx*fluc_Bx + fluc_Bx*fluc_Bx);

				rho = var[0].value[0][j][k]; rhoVx = var[1].value[0][j][k];
				fluc_rhoVx = (-1) * rho*(Vkz*fluc_V*sin(Vkz*z - omega*time));
				Eng = Eng + 1. / 2 * (2 * rhoVx*fluc_rhoVx + fluc_rhoVx*fluc_rhoVx) / rho;

				rhoVx = rhoVx + fluc_rhoVx; var[1].value[0][j][k] = rhoVx;
				Bx = Bx + fluc_Bx;
				var[4].value[0][j][k] = Bx; var[7].value[0][j][k] = Eng;
				// at i=1
				Bx = var[4].value[1][j][k]; Eng = var[7].value[1][j][k];
				fluc_Bx = (-1) * (Bkz*fluc_B*sin(Bkz*z - omega*time));
				Eng = Eng + 1. / 2 * (2 * Bx*fluc_Bx + fluc_Bx*fluc_Bx);

				rho = var[0].value[1][j][k]; rhoVx = var[1].value[1][j][k];
				fluc_rhoVx = (-1) * rho*(Vkz*fluc_V*sin(Vkz*z - omega*time));
				Eng = Eng + 1. / 2 * (2 * rhoVx*fluc_rhoVx + fluc_rhoVx*fluc_rhoVx) / rho;

				rhoVx = rhoVx + fluc_rhoVx; var[1].value[1][j][k] = rhoVx;
				Bx = Bx + fluc_Bx;
				var[4].value[1][j][k] = Bx; var[7].value[1][j][k] = Eng;
				// Does it need to add an fluctuation on sub_var[[][][]?????
			}
	}

	// add fluctuation at left or right boundary
	if (position_fluc == leftBoundary)
	{
		for (i = 0; i<Grid_Num_x; i++)
			for (j=0; j<Grid_Num_y; j++)
			{
				x = X[i];
				// at k=0
				Bx = var[4].value[i][j][0]; Eng = var[7].value[i][j][0];
				fluc_Bx = Bkz*fluc_B*tanh(x)*cos(omega*time);
				Eng = Eng + 1. / 2 * (2 * Bx*fluc_Bx + fluc_Bx*fluc_Bx);

				rho = var[0].value[i][j][0]; rhoVx = var[1].value[i][j][0];
				fluc_rhoVx = rho*Vkz*fluc_V*tanh(x)*cos(omega*time);
				Eng = Eng + 1. / 2 * (2 * rhoVx*fluc_rhoVx + fluc_rhoVx*fluc_rhoVx) / rho;

				rhoVx = rhoVx + fluc_rhoVx; var[1].value[i][j][0] = rhoVx;
				Bx = Bx + fluc_Bx;
				var[4].value[i][j][0] = Bx; var[7].value[i][j][0] = Eng;
			}
	}

	// add fluctuation at neutral-line according to <karimabadi, JGR, 2004>, 
	if (position_fluc == Neutral_Line)
	{
		for (i=0; i<Grid_Num_x; i++)
			for (j=0;j<Grid_Num_y;j++)
				for (k=0; k<Grid_Num_z; k++)
				{
					x=X[i]; z=Z[k];
					//if (abs(x)<norm_lambda)
					{
						Bx=var[4].value[i][j][k]; Bz=var[6].value[i][j][k]; 
						rho = var[0].value[i][j][k];
						rhoVx = var[1].value[i][j][k]; rhoVz = var[3].value[i][j][k];
						Eng=var[7].value[i][j][k];
						fluc_Bx = Bkz*fluc_B*cos(Bkx*x)*sin(Bkz*z); fluc_Bz = -Bkx*fluc_B*sin(Bkx*x)*cos(Bkz*z);
						Eng=Eng+1./2*(2*Bx*fluc_Bx+fluc_Bx*fluc_Bx)+1./2*(2*Bz*fluc_Bz+fluc_Bz*fluc_Bz);
						fluc_rhoVx = rho*(Vkz*fluc_V*cos(Vkx*x)*sin(Vkz*z)); fluc_rhoVz = rho*(-Vkx*fluc_V*sin(Vkx*x)*cos(Vkz*z));
						Eng=Eng + 1. / 2 * (2 * rhoVx*fluc_rhoVx + fluc_rhoVx*fluc_rhoVx)/rho + 1. / 2 * (2 * rhoVz*fluc_rhoVz + fluc_rhoVz*fluc_rhoVz)/rho;
						Bx=Bx+fluc_Bx; Bz=Bz+fluc_Bz;
						rhoVx += fluc_rhoVx; rhoVz += fluc_rhoVz;
						var[4].value[i][j][k]=Bx; var[6].value[i][j][k]=Bz; 
						var[1].value[i][j][k] = rhoVx; var[3].value[i][j][k] = rhoVz;
						var[7].value[i][j][k]=Eng;
					}
					// Does it need to add an fluctuation on sub_var[[][][]?????
				}
	}
	
}

// set fluctuation boundary value
void set_bndry_fluc(VARIABLE *var, double time)
{
	double norm_lambda;
	double x, z;
	double Bkz, omega;
	double time_delay = 10.;
	double Bx, Bz, Eng, delta_Bx, fluc_psi;				// for consistence of P, change in B must be represented in change of Eng
	int i, j, k;
	fluc_psi = fluctuation_mag_psi;
	Bkz = fluc_B_kz;
	omega = fluc_B_kz*velocity_phase;
	// set fluctuated boundary value at z=up or down boundary
	if (position_fluc == x_Boundary || position_fluc == upperBoundary)
	{
		for (j = 0; j<Grid_Num_y; j++)
			for (k = 0; k<Grid_Num_z; k++)
			{
				z = Z[k];
				// at i=Grid_Num_x-2
				Bx = var[4].value[Grid_Num_x - 2][j][k]; Eng = var[7].value[Grid_Num_x - 2][j][k];
				delta_Bx = Bkz*fluc_psi*sin(Bkz*z - omega*time)*(1 - exp(-time / time_delay)) - Bx;
				Eng = Eng + 1. / 2 * (2 * Bx*delta_Bx + delta_Bx*delta_Bx);
				var[4].value[Grid_Num_x - 2][j][k] = Bx + delta_Bx; var[7].value[Grid_Num_x - 2][j][k] = Eng;
				// at i=Grid_Num_x-1
				Bx = var[4].value[Grid_Num_x - 1][j][k]; Eng = var[7].value[Grid_Num_x - 1][j][k];
				delta_Bx = Bkz*fluc_psi*sin(Bkz*z - omega*time)*(1 - exp(-time / time_delay)) - Bx;
				Eng = Eng + 1. / 2 * (2 * Bx*delta_Bx + delta_Bx*delta_Bx);
				var[4].value[Grid_Num_x - 1][j][k] = Bx + delta_Bx; var[7].value[Grid_Num_x - 1][j][k] = Eng;
				// Does it need to add an fluctuation on sub_var[[][][]?????
			}
	}
	// set fluctuated boundary value at z=up or down boundary
	if (position_fluc == x_Boundary || position_fluc == downBoundary)
	{
		for (j = 0; j<Grid_Num_y; j++)
			for (k = 0; k<Grid_Num_z; k++)
			{
				z = Z[k];
				// at i=0
				Bx = var[4].value[0][j][k]; Eng = var[7].value[0][j][k];
				delta_Bx = Bkz*fluc_psi*sin(Bkz*z - omega*time)*(1 - exp(-time / time_delay)) - Bx;
				Eng = Eng + 1. / 2 * (2 * Bx*delta_Bx + delta_Bx*delta_Bx);
				var[4].value[0][j][k] = Bx + delta_Bx; var[7].value[0][j][k] = Eng;
				// at i=1
				Bx = var[4].value[1][j][k]; Eng = var[7].value[1][j][k];
				delta_Bx = Bkz*fluc_psi*sin(Bkz*z - omega*time)*(1 - exp(-time / time_delay)) - Bx;
				Eng = Eng + 1. / 2 * (2 * Bx*delta_Bx + delta_Bx*delta_Bx);
				var[4].value[1][j][k] = Bx + delta_Bx; var[7].value[1][j][k] = Eng;
				// Does it need to add an fluctuation on sub_var[[][][]?????
			}
	}
	// add fluctuation at left or right boundary
	if (position_fluc == leftBoundary)
	{
		norm_lambda = normalized_lambda;
		for (i = 0; i<Grid_Num_x; i++)
			for (j = 0; j<Grid_Num_y; j++)
			{
				x = X[i];
				// at k=0
				Bx = var[4].value[i][j][0]; Eng = var[7].value[i][j][0];
				delta_Bx = Bkz*fluc_psi*tanh(x / norm_lambda)*sin(omega*time)*(1 - exp(-time / time_delay));
				Eng = Eng + 1. / 2 * (2 * Bx*delta_Bx + delta_Bx*delta_Bx);
				var[4].value[i][j][0] = Bx + delta_Bx; var[7].value[i][j][0] = Eng;
				// at k=1
				Bx = var[4].value[i][j][1]; Eng = var[7].value[i][j][1];
				delta_Bx = Bkz*fluc_psi*tanh(x / norm_lambda)*sin(omega*time)*(1 - exp(-time / time_delay));
				Eng = Eng + 1. / 2 * (2 * Bx*delta_Bx + delta_Bx*delta_Bx);
				var[4].value[i][j][1] = Bx + delta_Bx; var[7].value[i][j][1] = Eng;
			}
	}

}

// Step on variables
void step_on(VARIABLE *pointer, VARIABLE *current, BASIC_VARIABLE &pressure, BASIC_VARIABLE &eta, double time, double time_interv)
{
/*In order to save space of stack, use new to allocate memory on heap: start */
	VARIABLE *var_intmedit=new VARIABLE[8];
	VARIABLE *temp_current=new VARIABLE[3];
	BASIC_VARIABLE *temp_pre_eta_pointer=new BASIC_VARIABLE[2];
	BASIC_VARIABLE (*flux)[3]=new BASIC_VARIABLE[8][3];
	BASIC_VARIABLE &temp_pressure=temp_pre_eta_pointer[0];
	BASIC_VARIABLE &temp_eta=temp_pre_eta_pointer[1];
/*In order to save space of stack, use new to allocate memory on heap: end*/
/*Updating variables: start*/
	cal_flux(flux, pointer, current, pressure, eta);
	exclude_soucrce_half_update(var_intmedit, flux, time_interv, First);

	cal_current(temp_current, var_intmedit, Incomplete); 
	set_eta(temp_eta, var_intmedit, temp_current, time, Incomplete);
	cal_pressure(temp_pressure, var_intmedit, Incomplete);
// Last statement do a forward differentiating, so that the value at boundary has not been updated.
//   And current calculating must be done excluding boundary value. The same with pressure, eta, flux. 
	cal_flux(flux, var_intmedit, temp_current, temp_pressure, temp_eta, Incomplete);
	exclude_soucrce_half_update(pointer, flux, time_interv, Second);
	
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
	delete []flux;
/*delete dynamic memory: end*/
	//cout<<"Step_on Step_on!"<<endl;
}

// Calculating flux from variables
void cal_flux(BASIC_VARIABLE flux[][3], VARIABLE *pointer, VARIABLE *current,\
	BASIC_VARIABLE &pressure_obj, BASIC_VARIABLE &eta_obj, Type T)
{
	double rho, Vx, Vy, Vz;
	double Bx, By, Bz;
	double Energy, B_Energy;
	double pressure, eta;
	double current_x, current_y, current_z;
	double V_dot_B, V_cross_B_x, V_cross_B_y, V_cross_B_z;

	int i,j,k;
	for (i=0;i<Grid_Num_x-T;i++)
		for (j=0;j<Grid_Num_y-T;j++)
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
				eta=eta_obj.value[i][j][k];

				// rho_flux assignment
				flux[0][0].value[i][j][k]=rho*Vx;
				flux[0][1].value[i][j][k]=rho*Vy;
				flux[0][2].value[i][j][k]=rho*Vz;
				
				//V_flux
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

				// Energy_fulx
				V_dot_B=(Vx-di/rho*current_x)*Bx+(Vy-di/rho*current_y)*By \
					+(Vz-di/rho*current_z)*Bz;
				Energy=Energy-B_Energy+pressure;       // Energy flux brought by V
				flux[7][0].value[i][j][k]=Vx*Energy \
					+(Vx-di/rho*current_x)*2*B_Energy-Bx*V_dot_B+eta* \
					(current_y*Bz-current_z*By);
				flux[7][1].value[i][j][k]=Vy*Energy \
					+(Vy-di/rho*current_y)*2*B_Energy-By*V_dot_B+eta* \
					(current_z*Bx-current_x*Bz);
				flux[7][2].value[i][j][k]=Vz*Energy \
					+(Vz-di/rho*current_z)*2*B_Energy-Bz*V_dot_B+eta* \
					(current_x*By-current_y*Bx);
			}

	//cout<<"I'm Calculating flux!"<<endl;
}

void smooth(VARIABLE *pointer, double time)      // how many times do smooth
{
//	VARIABLE tempvar[8];   changed to dynamic array as followings
	VARIABLE *tempvar=new VARIABLE[8];
//	double caf0=1.-0.25*tanh(time/100.);        // caf0 could be changed to caf
	int i,j,k, n;
	for (n=0;n<8;n++)
		for(i=0;i<Grid_Num_x;i++)
			for(j=0;j<Grid_Num_y;j++)
				for(k=0;k<Grid_Num_z;k++)
				{
					tempvar[n].value[i][j][k]=pointer[n].value[i][j][k]-var_x[n][i];
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
			tempvar[n].average(0.996);                       // avrg2

		tempvar[0].boundary_set(Positive,Positive);
		tempvar[1].boundary_set(Negative,Positive);tempvar[2].boundary_set(Negative,Positive);tempvar[3].boundary_set(Positive,Negative);
		tempvar[4].boundary_set(Positive,Negative);tempvar[5].boundary_set(Positive,Negative);tempvar[6].boundary_set(Negative,Positive);
		tempvar[7].boundary_set(Positive,Positive);
	}
	
	for (n=0;n<8;n++)
		for(i=0;i<Grid_Num_x;i++)
			for(j=0;j<Grid_Num_y;j++)
				for(k=0;k<Grid_Num_z;k++)
				{
					sub_var[n][i][j][k]=tempvar[n].value[i][j][k];
					pointer[n].value[i][j][k]=tempvar[n].value[i][j][k]+var_x[n][i];
				}
// deleting dynamic memory space
	delete []tempvar;
}

void record(ofstream &timeout, int run_num, int record_step, double system_time, VARIABLE *pointer, \
	BASIC_VARIABLE &Pressure_obj, VARIABLE &current_y, BASIC_VARIABLE &Electric_y)
{
	ofstream out;
	char str[12]="time000.dat";
	int out_Grid_x, out_Grid_y, out_Grid_z;
	int i,j,k, n;

	out_Grid_x=(Grid_Num_x-1)/num_out;
	out_Grid_y=(Grid_Num_y-1)/num_out;
	out_Grid_z=(Grid_Num_z-1)/num_out;

	num2str(str+4, run_num);

	out.open(str);
	out<<setiosflags(ios::scientific)<<setprecision(5);
	out<<"title = \"Hall Magnetic Reconnection\""<<endl;
	out<<"variables = \"x\"";
	if (out_Grid_y != 0)
		out<<", \"y\"";
	out<<", \"z\"";
	out<<", \"rho\", \"rhoVx\", \"rhoVy\", \"rhoVz\"";
	out<<", \"Bx\", \"By\", \"Bz\"";
	out<<", \"Pressure\", \"Jy\", \"Ey\""<<endl;
	out<<"zone t = \" nstep ="<<record_step<<" \""<<endl;
	out<<" strandid = 1, "<<"solutiontime = "<<system_time<<endl;
	out<<" i = "<<num_out+1;
	if (out_Grid_y != 0)
		out<<" , j = "<<num_out+1;
	out<<" , k = "<<num_out+1<<endl;

	for (i=0;i<=num_out;i++)
	{
		if (out_Grid_y==0)
		{
			for (k=0;k<=num_out;k++)
			{
				out<<setw(15)<<X[i*out_Grid_x]<<" "<<setw(13)<<Z[k*out_Grid_z]<<" ";
				for (n=0; n<7; n++)
					out<<setw(13)<<pointer[n].value[i*out_Grid_x][0][k*out_Grid_z]<<" ";
				out<<setw(13)<<Pressure_obj.value[i*out_Grid_x][0][k*out_Grid_z]<<" ";
				out<<setw(13)<<current_y.value[i*out_Grid_x][0][k*out_Grid_z]<<" ";
				out<<setw(13)<<Electric_y.value[i*out_Grid_x][0][k*out_Grid_z]<<" ";
				out<<setw(13)<<endl;
			}
		}
		else
		{
			for (j=0;j<=num_out;j++)
			{
				for (k=0;k<=num_out;k++)
				{
					out<<setw(15)<<X[i*out_Grid_x]<<" "\
						<<setw(13)<<Y[j*out_Grid_y]<<" "<<setw(13)<<Z[k*out_Grid_z]<<" ";
					for (n=0; n<7; n++)
						out<<setw(13)<<pointer[i].value[i*out_Grid_x][j*out_Grid_y][k*out_Grid_z]<<" ";
					out<<setw(13)<<Pressure_obj.value[i*out_Grid_x][j*out_Grid_y][k*out_Grid_z]<<" ";
					out<<setw(13)<<current_y.value[i*out_Grid_x][j*out_Grid_y][k*out_Grid_z]<<" ";
					out<<setw(13)<<Electric_y.value[i*out_Grid_x][j*out_Grid_y][k*out_Grid_z]<<" ";
					out<<setw(13)<<endl;
				}
			}
		}
	}
	timeout<<setw(6)<<record_step<<setw(6)<<" "<<setiosflags(ios::scientific)\
			<<setprecision(5)<<setw(15)<<system_time<<endl;
	out.close();
//	cout<<"record() is called;"<<endl;
}
