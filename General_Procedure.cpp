// General_Procedure.cpp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
#include "Basic_Parameter.h"
#include "Variables_Definition.h"
#include "Procedure.h"


// Clarification of mesh-grid
extern double X[], Y[], Z[], X_interval[], Y_interval[], Z_interval[];
extern double var_x[][Grid_Num_x], var_x_plushalfdx[][Grid_Num_x];          // Declaration of external variables
extern double sub_var[8][Grid_Num_x][Grid_Num_y][Grid_Num_z];

// Integarting variables for half dt from Fluxes using 2-order Lax-Wendroff method
// First time (Order o=1) forward difference; Second time (Order o=2) backward difference
void exclude_soucrce_half_update(VARIABLE *update_var, BASIC_VARIABLE flux[][3], \
	double time_interv, Order o)
{	
	double Temp_var;
	double dtflux_x, dtflux_y, dtflux_z;                    // for fdx, gdy, hdz
	double dt;
	double f000, f100, f010, f001, f110, f101, f011, f111;   // eight point of value
	double v000, v100, v010, v001, v110, v101, v011, v111;
	double dx, dy, dz;                                       // spatial interval
	int i,j,k,n;
	switch (o)
	{
	case First:
		{
			dt=0.5*time_interv;
			for (n=0;n<8;n++)
			{				
				for (i=0;i<Grid_Num_x-1;i++)
					for (j=0;j<Grid_Num_y-1;j++)
						for (k=0;k<Grid_Num_z-1;k++)
						{
							dx=X_interval[i];                      
							dy=Y_interval[j];
							dz=Z_interval[k];
// following used is the most simple addressing method, which could be improved
							f000=flux[n][0].value[i][j][k];        // simple addressing method, be able to improve						
							f100=flux[n][0].value[i+1][j][k];      // but times of computing remains the same
							f010=flux[n][0].value[i][j+1][k];
							f001=flux[n][0].value[i][j][k+1];
							f110=flux[n][0].value[i+1][j+1][k];
							f101=flux[n][0].value[i+1][j][k+1];
							f011=flux[n][0].value[i][j+1][k+1];
							f111=flux[n][0].value[i+1][j+1][k+1];

							dtflux_x=-dt*(f111-f011 \
									    +f101-f001 \
										+f110-f010 \
										+f100-f000)/(4.*dx);      // attention to the mius sign  "-" in front of dt

							f000=flux[n][1].value[i][j][k];       						
							f100=flux[n][1].value[i+1][j][k];      
							f010=flux[n][1].value[i][j+1][k];
							f001=flux[n][1].value[i][j][k+1];
							f110=flux[n][1].value[i+1][j+1][k];
							f101=flux[n][1].value[i+1][j][k+1];
							f011=flux[n][1].value[i][j+1][k+1];
							f111=flux[n][1].value[i+1][j+1][k+1];
							dtflux_y=-dt*(f111-f101 \
									    +f110-f100 \
										+f011-f001 \
										+f010-f000)/(4*dy);

							f000=flux[n][2].value[i][j][k];      				
							f100=flux[n][2].value[i+1][j][k];   
							f010=flux[n][2].value[i][j+1][k];
							f001=flux[n][2].value[i][j][k+1];
							f110=flux[n][2].value[i+1][j+1][k];
							f101=flux[n][2].value[i+1][j][k+1];
							f011=flux[n][2].value[i][j+1][k+1];
							f111=flux[n][2].value[i+1][j+1][k+1];
							dtflux_z=-dt*(f111-f110 \
									    +f011-f010 \
										+f101-f100 \
										+f001-f000)/(4*dz);

							v000=sub_var[n][i][j][k];      					
							v100=sub_var[n][i+1][j][k];      
							v010=sub_var[n][i][j+1][k];
							v001=sub_var[n][i][j][k+1];
							v110=sub_var[n][i+1][j+1][k];
							v101=sub_var[n][i+1][j][k+1];
							v011=sub_var[n][i][j+1][k+1];
							v111=sub_var[n][i+1][j+1][k+1];

							Temp_var=(v000+v100+v010+v001+v110+v101+v011+v111)/8.;

							update_var[n].value[i][j][k]=var_x_plushalfdx[n][i]+Temp_var+dtflux_x+dtflux_y+dtflux_z;						
						}
			}			
			break;
		}
	case Second:
		{
			dt=time_interv;
			for (n=0;n<8;n++)
			{				
				for (i=1;i<Grid_Num_x-1;i++)
					for (j=1;j<Grid_Num_y-1;j++)
						for (k=1;k<Grid_Num_z-1;k++)
						{
							dx=X_interval[i-1];
							dy=Y_interval[j-1];
							dz=Z_interval[k-1];
// following used is the most simple addressing method, which could be improved
							f000=flux[n][0].value[i-1][j-1][k-1];       							
							f100=flux[n][0].value[i-1+1][j-1][k-1];      
							f010=flux[n][0].value[i-1][j-1+1][k-1];
							f001=flux[n][0].value[i-1][j-1][k-1+1];
							f110=flux[n][0].value[i-1+1][j-1+1][k-1];
							f101=flux[n][0].value[i-1+1][j-1][k-1+1];
							f011=flux[n][0].value[i-1][j-1+1][k-1+1];
							f111=flux[n][0].value[i-1+1][j-1+1][k-1+1];

							dtflux_x=-dt*(f111-f011 \
									    +f101-f001 \
										+f110-f010 \
										+f100-f000)/(4*dx);

							f000=flux[n][1].value[i-1][j-1][k-1];       						
							f100=flux[n][1].value[i-1+1][j-1][k-1];     
							f010=flux[n][1].value[i-1][j-1+1][k-1];
							f001=flux[n][1].value[i-1][j-1][k-1+1];
							f110=flux[n][1].value[i-1+1][j-1+1][k-1];
							f101=flux[n][1].value[i-1+1][j-1][k-1+1];
							f011=flux[n][1].value[i-1][j-1+1][k-1+1];
							f111=flux[n][1].value[i-1+1][j-1+1][k-1+1];
							dtflux_y=-dt*(f111-f101 \
									    +f110-f100 \
										+f011-f001 \
										+f010-f000)/(4*dy);

							f000=flux[n][2].value[i-1][j-1][k-1];        						
							f100=flux[n][2].value[i-1+1][j-1][k-1];     
							f010=flux[n][2].value[i-1][j-1+1][k-1];
							f001=flux[n][2].value[i-1][j-1][k-1+1];
							f110=flux[n][2].value[i-1+1][j-1+1][k-1];
							f101=flux[n][2].value[i-1+1][j-1][k-1+1];
							f011=flux[n][2].value[i-1][j-1+1][k-1+1];
							f111=flux[n][2].value[i-1+1][j-1+1][k-1+1];
							dtflux_z=-dt*(f111-f110 \
									    +f011-f010 \
										+f101-f100 \
										+f001-f000)/(4*dz);
							update_var[n].value[i][j][k]=update_var[n].value[i][j][k]+dtflux_x+dtflux_y+dtflux_z;						
						}
			}			
			break;
		}
	default:
		{
			;//cout<<"Do you know what is the default method to calculate divergence of the flux in 3D situation? "<<endl;
			//cout<<"I will do nothing if you call default case!"<<endl;
		}
	}
	//cout <<"Exclude_Soucrce_half_Update invoked!"<<endl;
	switch (o)
	{
	case First:
		;
	case Second:
		;
	default:
		;//cout<<"But it's default call! And I will do nothing until you tell me what is the default method to calculate divergence of the flux in 3D situation? "<<endl;	
	}
}

// Update variables from source term
void open_var_files(ofstream *var_out)
{
	var_out[0].open("rho.dat");
	var_out[1].open("rhoVx.dat");var_out[2].open("rhoVy.dat");var_out[3].open("rhoVz.dat"); 
	var_out[4].open("Bx.dat");   var_out[5].open("By.dat");   var_out[6].open("Bz.dat");
	var_out[7].open("Eng.dat");
}

void close_var_files(ofstream *var_out)
{
	int i;
	for(i=0;i<8;i++)
		var_out[i].close();
}

void open_var_files(ifstream *var_in)
{
	var_in[0].open("rho.dat");
	var_in[1].open("rhoVx.dat");var_in[2].open("rhoVy.dat");var_in[3].open("rhoVz.dat"); 
	var_in[4].open("Bx.dat");   var_in[5].open("By.dat");   var_in[6].open("Bz.dat");
	var_in[7].open("Eng.dat");
}

void close_var_files(ifstream *var_in)
{
	int i;
	for(i=0;i<8;i++)
		var_in[i].close();
}

void source_update(VARIABLE *update_var, double time_interv)
{
	//cout<<"Source_Update invoked! But there is no source term!"<<endl;
}

void copy(VARIABLE *update_var, VARIABLE *mother_var)
{
	int i,j,k,n;
	for (n=0;n<8;n++)
		for (i=0;i<Grid_Num_x;i++)
			for (j=0;j<Grid_Num_y;j++)
				for (k=0;k<Grid_Num_z;k++)
					update_var[n].value[i][j][k]=mother_var[n].value[i][j][k];

	//cout<<"Copy invoked!"<<endl;
}

void num2str( char *str, int n )
{
	int hundred_num, deci_num, digi_num;
	
	hundred_num=n/100;
	deci_num=(n-100*hundred_num)/10;
	digi_num=n-100*hundred_num-10*deci_num;
	hundred_num+=48;
	deci_num+=48;
	digi_num+=48;
	
	str[0] = char(hundred_num);
	str[1] = char(deci_num);
	str[2] = char(digi_num);
}