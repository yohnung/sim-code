#include <cmath>
#include "Variables_Definition.h"
#include "Procedure.h"


// Clarification of mesh-grid
extern double X[], Y[], Z[], X_interval[], Y_interval[], Z_interval[];
extern double var_x[][Grid_Num_x], var_x_plushalfdx[][Grid_Num_x];          // 声明声明声明声明

void make_pressure_positive(BASIC_VARIABLE & pressure_obj, double positive_value)
{
	//char * file_name;
	ofstream out("pressure_is_negative.txt");
	int i,j,k;
	for(i=0;i<Grid_Num_x;i++)
	{
		for(j=0;j<Grid_Num_y;j++)
		{
			for(k=0;k<Grid_Num_z;k++)
			{
				if(pressure_obj.value[i][j][k]<0.)
					out<<"Oops, pressure is negative, and program can be stropped!!!"<<endl;
				if(pressure_obj.value[i][j][k]<positive_value)
					pressure_obj.value[i][j][k]=positive_value;
			}
		}
	}
	;
}

// Integarting variables for half dt from Fluxes using 2-order Lax-Wendroff method
// First time (Order o=1) forward difference; Second time (Order o=2) backward difference
void exclude_soucrce_half_update(VARIABLE *update_var, VARIABLE *parameter_var, \
	BASIC_VARIABLE flux[][3], double time_interv, Order o)
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
				{
					for (j=0;j<Grid_Num_y-1;j++)
					{
						for (k=0;k<Grid_Num_z-1;k++)
						{
							dx=X_interval[i];                      
							dy=Y_interval[j];
							dz=Z_interval[k];
							f000=flux[n][0].value[i][j][k];        // 这是最简单寻址方式，可以改进							
							f100=flux[n][0].value[i+1][j][k];      // 但就是改进了寻址方式，用于计算的数量却没有少
							f010=flux[n][0].value[i][j+1][k];
							f001=flux[n][0].value[i][j][k+1];
							f110=flux[n][0].value[i+1][j+1][k];
							f101=flux[n][0].value[i+1][j][k+1];
							f011=flux[n][0].value[i][j+1][k+1];
							f111=flux[n][0].value[i+1][j+1][k+1];

							dtflux_x=dt*(f111-f011 \
									    +f101-f001 \
										+f110-f010 \
										+f100-f000)/(4*dx);

							f000=flux[n][1].value[i][j][k];        // 这是最简单寻址方式，可以改进							
							f100=flux[n][1].value[i+1][j][k];      // 但就是改进了寻址方式，用于计算的数量却没有少
							f010=flux[n][1].value[i][j+1][k];
							f001=flux[n][1].value[i][j][k+1];
							f110=flux[n][1].value[i+1][j+1][k];
							f101=flux[n][1].value[i+1][j][k+1];
							f011=flux[n][1].value[i][j+1][k+1];
							f111=flux[n][1].value[i+1][j+1][k+1];
							dtflux_y=dt*(f111-f101 \
									    +f110-f100 \
										+f011-f001 \
										+f010-f000)/(4*dy);

							f000=flux[n][2].value[i][j][k];        // 这是最简单寻址方式，可以改进							
							f100=flux[n][2].value[i+1][j][k];      // 但就是改进了寻址方式，用于计算的数量却没有少
							f010=flux[n][2].value[i][j+1][k];
							f001=flux[n][2].value[i][j][k+1];
							f110=flux[n][2].value[i+1][j+1][k];
							f101=flux[n][2].value[i+1][j][k+1];
							f011=flux[n][2].value[i][j+1][k+1];
							f111=flux[n][2].value[i+1][j+1][k+1];
							dtflux_z=dt*(f111-f110 \
									    +f011-f010 \
										+f101-f100 \
										+f001-f000)/(4*dz);

							v000=parameter_var[n].value[i][j][k];        // 这是最简单寻址方式，可以改进							
							v100=parameter_var[n].value[i+1][j][k];      // 但就是改进了寻址方式，用于计算的数量却没有少
							v010=parameter_var[n].value[i][j+1][k];
							v001=parameter_var[n].value[i][j][k+1];
							v110=parameter_var[n].value[i+1][j+1][k];
							v101=parameter_var[n].value[i+1][j][k+1];
							v011=parameter_var[n].value[i][j+1][k+1];
							v111=parameter_var[n].value[i+1][j+1][k+1];

							Temp_var=(v000+v100+v010+v001+v110+v101+v011+v111)/8.;

							update_var[n].value[i][j][k]=var_x_plushalfdx[n][i]+Temp_var+dtflux_x+dtflux_y+dtflux_z;						
						}
					}
				}
			}			
			break;
		}
	case Second:
		{
			dt=0.5*time_interv;
			for (n=0;n<8;n++)
			{				
				for (i=1;i<Grid_Num_x-1;i++)
				{
					for (j=1;j<Grid_Num_y-1;j++)
					{
						for (k=1;k<Grid_Num_z-1;k++)
						{
							dx=X_interval[i-1];
							dy=Y_interval[j-1];
							dz=Z_interval[k-1];
							f000=flux[n][0].value[i-1][j-1][k-1];        // 这是最简单寻址方式，可以改进							
							f100=flux[n][0].value[i-1+1][j-1][k-1];      // 但就是改进了寻址方式，用于计算的数量却没有少
							f010=flux[n][0].value[i-1][j-1+1][k-1];
							f001=flux[n][0].value[i-1][j-1][k-1+1];
							f110=flux[n][0].value[i-1+1][j-1+1][k-1];
							f101=flux[n][0].value[i-1+1][j-1][k-1+1];
							f011=flux[n][0].value[i-1][j-1+1][k-1+1];
							f111=flux[n][0].value[i-1+1][j-1+1][k-1+1];

							dtflux_x=dt*(f111-f011 \
									    +f101-f001 \
										+f110-f010 \
										+f100-f000)/(4*dx);

							f000=flux[n][1].value[i-1][j-1][k-1];        // 这是最简单寻址方式，可以改进							
							f100=flux[n][1].value[i-1+1][j-1][k-1];      // 但就是改进了寻址方式，用于计算的数量却没有少
							f010=flux[n][1].value[i-1][j-1+1][k-1];
							f001=flux[n][1].value[i-1][j-1][k-1+1];
							f110=flux[n][1].value[i-1+1][j-1+1][k-1];
							f101=flux[n][1].value[i-1+1][j-1][k-1+1];
							f011=flux[n][1].value[i-1][j-1+1][k-1+1];
							f111=flux[n][1].value[i-1+1][j-1+1][k-1+1];
							dtflux_y=dt*(f111-f101 \
									    +f110-f100 \
										+f011-f001 \
										+f010-f000)/(4*dy);

							f000=flux[n][2].value[i-1][j-1][k-1];        // 这是最简单寻址方式，可以改进							
							f100=flux[n][2].value[i-1+1][j-1][k-1];      // 但就是改进了寻址方式，用于计算的数量却没有少
							f010=flux[n][2].value[i-1][j-1+1][k-1];
							f001=flux[n][2].value[i-1][j-1][k-1+1];
							f110=flux[n][2].value[i-1+1][j-1+1][k-1];
							f101=flux[n][2].value[i-1+1][j-1][k-1+1];
							f011=flux[n][2].value[i-1][j-1+1][k-1+1];
							f111=flux[n][2].value[i-1+1][j-1+1][k-1+1];
							dtflux_z=dt*(f111-f110 \
									    +f011-f010 \
										+f101-f100 \
										+f001-f000)/(4*dz);
							update_var[n].value[i][j][k]=parameter_var[n].value[i][j][k]+dtflux_x+dtflux_y+dtflux_z;						
						}
					}
				}
			}			
			break;
		}
	default:
		{
			//cout<<"Do you know what is the default method to calculate divergence of the flux in 3D situation? "<<endl;
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
		{
			//cout<<"But it's default call! And I will do nothing until you tell me what is the default method to calculate divergence of the flux in 3D situation? "<<endl;
		}	
	}
}

// Update variables from source term
void source_update(VARIABLE *update_var, double time_interv)
{
	//cout<<"Source_Update invoked! But there is no source term!"<<endl;
}

void boudary_set(VARIABLE &variable, Symmetry_Type sign_x, Symmetry_Type sign_z)
{
	int i,j,k;
	// Magnet0sheath and -pause boundary
	// Inflow boundary ,    equivalent extrapolation!!!!????????
	if (half_x==True)
	{
		for (j=1;j<Grid_Num_y-1;j++)
		{
			for (k=1;k<Grid_Num_z-1;k++)
			{
				variable.value[0][j][k]=variable.value[1][j][k];
				variable.value[Grid_Num_x-1][j][k]=sign_x*variable.value[Grid_Num_x-3][Grid_Num_y-1-j][k];    // minus axial symmetry
			}
		}
	}
	else
	{
		for (j=1;j<Grid_Num_y-1;j++)
		{
			for (k=1;k<Grid_Num_z-1;k++)
			{
				variable.value[0][j][k]=variable.value[1][j][k];
				variable.value[Grid_Num_x-1][j][k]=variable.value[Grid_Num_x-2][j][k];
			}
		}
	}		
		
	// outflow boundary,     equiv. extrap.
	if (period_y==True)
	{		
		for(i=0;i<Grid_Num_x;i++)
		{
			for(k=1;k<Grid_Num_z-1;k++)
			{
				variable.value[i][0][k]=variable.value[i][Grid_Num_z-2][k];
				variable.value[i][Grid_Num_z-1][k]=variable.value[i][1][k];
			}
		}
	}
	else
	{
		for(i=0;i<Grid_Num_x;i++)
		{
			for(k=1;k<Grid_Num_z-1;k++)
			{
				variable.value[i][0][k]=variable.value[i][1][k];
				variable.value[i][Grid_Num_z-1][k]=variable.value[i][Grid_Num_z-2][k];
			}
		}
	}

	if (half_z==True)
	{
		for (i=0;i<Grid_Num_x;i++)
		{
			for (j=0;j<Grid_Num_y;j++)
			{
				variable.value[i][j][0]=variable.value[i][j][1];
				variable.value[i][j][Grid_Num_z-1]=sign_z*variable.value[i][j][Grid_Num_z-3];
			}
		}
	}
	else
	{
		for (i=0;i<Grid_Num_x;i++)
		{
			for (j=0;j<Grid_Num_y;j++)
			{
				variable.value[i][j][0]=variable.value[i][j][1];
				variable.value[i][j][Grid_Num_z-1]=variable.value[i][j][Grid_Num_z-2];
			}
		}
	}    
}

void copy(VARIABLE *update_var, VARIABLE *mother_var)
{
	int i,j,k,n;
	for (n=0;n<8;n++)
	{
		for (i=0;i<Grid_Num_x;i++)
		{
			for (j=0;j<Grid_Num_y;j++)
			{
				for (k=0;k<Grid_Num_z;k++)
				{
					update_var[n].value[i][j][k]=mother_var[n].value[i][j][k];
				}
			}
		}
	}
	//cout<<"Copy invoked!"<<endl;
}

/* 已经改写成了类的成员函数
void smooth_xyz(double var[][Grid_Num_x][Grid_Num_y][Grid_Num_z], int times)
{
	double temp_var[Grid_Num_x][Grid_Num_y][Grid_Num_z];
	int i,j,k, m,n;
	double theta;
	for (m=0;m<times;m++)
	{
		for(n=0;n<8;n++)
		{
			for(i=1;i<Grid_Num_x-1;i++)
			{
				for(j=1;j<Grid_Num_y-1;j++)
				{
					for(k=1;k<Grid_Num_z-1;k++)
					{
						temp_var[i][j][k]=var[n][i-1][j][k]+var[n][i+1][j][k]+var[n][i][j-1][k] +\
							var[n][i][j+1][k]+var[n][i][j][k-1]+var[n][i][j][k+1]-6*var[n][i][j][k];
					}
				}
			}

			// Here starts smooth
			// x-direction
			for(i=1;i<Num_Smooth_x;i++)
			{
				for(j=1;j<Grid_Num_y-1;j++)
				{
					for(k=1;k<Grid_Num_z-1;k++)
					{
						theta=3.1415926*(i-1)/(Num_Smooth_x-3);              //...........??????没懂这是什么意思
						var[n][i][j][k]=var[n][i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
							temp_var[i][j][k];                               //...........??????没懂这是什么意思						
					}
				}
			}
			if(half_x==False)
			{
				for(i=Grid_Num_x-Num_Smooth_x;i<Grid_Num_x-1;i++)
				{
					for(j=1;j<Grid_Num_y-1;j++)
					{
						for(k=1;k<Grid_Num_z-1;k++)
						{
							theta=3.1415926*(Grid_Num_x-i-2)/(Num_Smooth_x-3);              //...........??????没懂这是什么意思
							var[n][i][j][k]=var[n][i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
								temp_var[i][j][k];                               //...........??????没懂这是什么意思						
						}
					}
				}
			}
			// y-direction
			if(period_y==False)
			{
				for(i=1;i<Grid_Num_x-1;i++)
				{
					for(j=1;j<Num_Smooth_y;j++)
					{
						for(k=1;k<Grid_Num_z-1;k++)
						{
							theta=3.1415926*(j-1)/(Num_Smooth_x-3);              //...........??????没懂这是什么意思
							var[n][i][j][k]=var[n][i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
								temp_var[i][j][k];                               //...........??????没懂这是什么意思						
						}
					}
				}
				for(i=1;i<Grid_Num_x-1;i++)
				{
					for(j=Grid_Num_y-Num_Smooth_y;j<Grid_Num_y-1;j++)
					{
						for(k=1;k<Grid_Num_z-1;k++)
						{
							theta=3.1415926*(Grid_Num_y-j-2)/(Num_Smooth_x-3);              //...........??????没懂这是什么意思
							var[n][i][j][k]=var[n][i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
								temp_var[i][j][k];                               //...........??????没懂这是什么意思						
						}
					}
				}
			}
			/*
			if (period_y==True)
			{
			;
			}
			
			// z-direction
			for(i=1;i<Grid_Num_x-1;i++)
			{
				for(j=1;j<Grid_Num_y-1;j++)
				{
					for(k=1;k<Grid_Num_z-1;k++)
					{
						// theta=2*3.1415926*Z[k]/Z_min;              //...........??????没懂这是什么意思
						var[n][i][j][k]=var[n][i][j][k]+(1./48.)* \
							temp_var[i][j][k];                        //...........??????没懂这是什么意思
						// (1./48.)* (2.+cos(theta)/3.*temp_var[i][j][k];       
					}
				}
			}
			if(half_z==False)
			{
				for(i=1;i<Grid_Num_x-1;i++)
				{
					for(j=1;j<Grid_Num_y-1;j++)
					{
						for(k=Grid_Num_z-Num_Smooth_z;k<Grid_Num_z-1;k++)
						{
							theta=3.1415926*(Grid_Num_z-k-2)/(Num_Smooth_x-3);         //...........??????没懂这是什么意思
							var[n][i][j][k]=var[n][i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
								temp_var[i][j][k];                                     //...........??????没懂这是什么意思
						}
					}
				}
			}
		}
		pointer[0].boundary_set(Positive,Positive);
		pointer[1].boundary_set(Negative,Positive);pointer[2].boundary_set(Negative,Positive);pointer[3].boundary_set(Positive,Negative);
		pointer[4].boundary_set(Positive,Negative);pointer[5].boundary_set(Positive,Negative);pointer[6].boundary_set(Negative,Positive);
		pointer[7].boundary_set(Positive,Positive);
	}
}

void average(double var[][Grid_Num_x][Grid_Num_y][Grid_Num_z],double aver_coeff)
{
	double temp_var[Grid_Num_x][Grid_Num_y][Grid_Num_z];
	int i,j,k n;
	for (n=0;n<8;n++)
	{
		for(i=1;i<Grid_Num_x-1;i++)
		{
			for(j=1;j<Grid_Num_y-1;j++)
			{
				for(k=1;k<Grid_Num_z-1;k++)
				{
					temp_var[i][j][k]=aver_coeff*var[n][i][j][k]+ \
						(1-aver_coeff)*( var[n][i-1][j][k]+var[n][i+1][j][k]+ \
						var[n][i][j-1][k]+var[n][i][j+1][k]+ \
						var[n][i][j][k-1]+var[n][i][j][k+1] )/6.;
				}
			}
		}
		for(i=1;i<Grid_Num_x-1;i++)
		{
			for(j=1;j<Grid_Num_y-1;j++)
			{
				for(k=1;k<Grid_Num_z-1;k++)
				{
					var[n][i][j1][k]=temp_var[i][j][k];
				}
			}
		}
	}
}

void average_B(double var[][Grid_Num_x][Grid_Num_y][Grid_Num_z],double aver_coeff)
{
	double temp_var[Grid_Num_x][Grid_Num_y][Grid_Num_z];
	int i,j,k n;
	for (n=4;n<7;n++)
	{
		for(i=1;i<Grid_Num_x-1;i++)
		{
			for(j=1;j<Grid_Num_y-1;j++)
			{
				for(k=1;k<Grid_Num_z-1;k++)
				{
					temp_var[i][j][k]=aver_coeff*var[n][i][j][k]+ \
						(1-aver_coeff)*( var[n][i-1][j][k]+var[n][i+1][j][k]+ \
						var[n][i][j-1][k]+var[n][i][j+1][k]+ \
						var[n][i][j][k-1]+var[n][i][j][k+1] )/6.;
				}
			}
		}
		for(i=1;i<Grid_Num_x-1;i++)
		{
			for(j=1;j<Grid_Num_y-1;j++)
			{
				for(k=1;k<Grid_Num_z-1;k++)
				{
					var[n][i][j1][k]=temp_var[i][j][k];
				}
			}
		}
	}
}
已经改写成了类的成员函数*/