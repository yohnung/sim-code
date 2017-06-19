// Variables_Definition.cpp
// 2D output is addedto output file  

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
#include "Basic_Parameter.h"
#include "Variables_Definition.h"

// declarification of mesh-grid
extern double X[], Y[], Z[], X_interval[], Y_interval[], Z_interval[];

// Recording variables' value
void BASIC_VARIABLE::record(ofstream &out_obj)       // recrd1
			   {
				   int out_Grid_x, out_Grid_y, out_Grid_z;
				   int i,j,k;
				   out_Grid_x=(Grid_Num_x-1)/num_out;
				   out_Grid_y=(Grid_Num_y-1)/num_out;
				   out_Grid_z=(Grid_Num_z-1)/num_out;
				   out_obj<<setiosflags(ios::scientific)<<setprecision(15);
				   for (i=0;i<=num_out;i++)
				   {
					   if (out_Grid_y==0)
					   {
						   for (k=0;k<=num_out;k++)
							   out_obj<<" "<<value[i*out_Grid_x][0][k*out_Grid_z];
						   out_obj<<endl;
					   }
					   else
					   {
						   for (j=0;j<=num_out;j++)
						   {
							   for (k=0;k<=num_out;k++)
								   out_obj<<" "<<value[i*out_Grid_x][j*out_Grid_y][k*out_Grid_z];
							   out_obj<<endl;
						   }
					   }
				   }
				   //cout<<"BASIC_VARIABLE::record invoked! But I don't write anything, just to show up!"<<endl;
				   //Recording the value
			   }

// Setting variables' boundary at 0(Start-point) or Grid_Num_?-1-T(End-point)
void VARIABLE::boundary_set(Symmetry_Type sign_x, Symmetry_Type sign_z)		
			{
				int i,j,k;
				// Magnet0sheath and -pause boundary
				// Inflow boundary ,    equivalent extrapolation!!!!????????
				if (half_x==True)
				{
					if (x_fixed_bndry == False)
					{
						for (j=1;j<Grid_Num_y-1;j++)
							for (k=1;k<Grid_Num_z-1;k++)
							{
								value[0][j][k]=value[1][j][k];
								value[Grid_Num_x-1][j][k]=sign_x*value[Grid_Num_x-3][Grid_Num_y-1-j][k];    // minus axial symmetry
							}	
					}
					else
					{
						for (j=1;j<Grid_Num_y-1;j++)
							for (k=1;k<Grid_Num_z-1;k++)
								value[Grid_Num_x-1][j][k]=sign_x*value[Grid_Num_x-3][Grid_Num_y-1-j][k];    // minus axial symmetry
					}
				}
				else
				{
					if (x_fixed_bndry == False)
						for (j=1;j<Grid_Num_y-1;j++)
							for (k=1;k<Grid_Num_z-1;k++)
							{
								value[0][j][k]=value[1][j][k];
								value[Grid_Num_x-1][j][k]=value[Grid_Num_x-2][j][k];
							}					
				}		
		
				// outflow boundary,     equiv. extrap.
				if (period_y==True)
				{		
					for(i=x_fixed_bndry;i<Grid_Num_x-(1-half_x)*x_fixed_bndry;i++)
						for(k=1;k<Grid_Num_z-1;k++)
						{
							value[i][0][k]=value[i][Grid_Num_y-2][k];
							value[i][Grid_Num_y-1][k]=value[i][1][k];
						}
				}
				else
				{
					for(i=x_fixed_bndry;i<Grid_Num_x-(1-half_x)*x_fixed_bndry;i++)
						for(k=1;k<Grid_Num_z-1;k++)
						{
							value[i][0][k]=value[i][1][k];
							value[i][Grid_Num_y-1][k]=value[i][Grid_Num_y-2][k];
						}
				}

				if (half_z==True)
				{
					for (i=x_fixed_bndry;i<Grid_Num_x-(1-half_x)*x_fixed_bndry;i++)
						for (j=0;j<Grid_Num_y;j++)
						{
							value[i][j][0]=value[i][j][1];
							value[i][j][Grid_Num_z-1]=sign_z*value[i][j][Grid_Num_z-3];
						}
				}
				else
				{
					for (i=x_fixed_bndry;i<Grid_Num_x-(1-half_x)*x_fixed_bndry;i++)
						for (j=0;j<Grid_Num_y;j++)
						{
							value[i][j][0]=value[i][j][1];
							value[i][j][Grid_Num_z-1]=value[i][j][Grid_Num_z-2];
						}
				}   
				//cout<<"VARIABLE::basic_boundary_set invoked! Have set the boundary for you!  for you!!!"<<endl;
			}

// set right boundary 

void VARIABLE::smooth_xyz(int times)
{
	int i,j,k, m;
	double theta;
	double (*temp_var)[Grid_Num_y][Grid_Num_z]= new double[Grid_Num_x][Grid_Num_y][Grid_Num_z];
	for (m=0;m<times;m++)
	{		
		for(i=1;i<Grid_Num_x-1;i++)
		{
			for(j=1;j<Grid_Num_y-1;j++)
			{
				for(k=1;k<Grid_Num_z-1;k++)
				{
					temp_var[i][j][k]=value[i-1][j][k]+value[i+1][j][k]+value[i][j-1][k] +\
						value[i][j+1][k]+value[i][j][k-1]+value[i][j][k+1]-6*value[i][j][k];
				}
			}
		}

		/* Here starts smooth*/
		// x-direction
		for(i=1;i<Num_Smooth_x;i++)
		{
			for(j=1;j<Grid_Num_y-1;j++)
			{
				for(k=1;k<Grid_Num_z-1;k++)
				{
					theta=3.1415926*(i-1)/(Num_Smooth_x-3);              //...........?????
					value[i][j][k]=value[i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
						temp_var[i][j][k];                               //...........?????						
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
						theta=3.1415926*(Grid_Num_x-i-2)/(Num_Smooth_x-3);              //...........?????
						value[i][j][k]=value[i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
							temp_var[i][j][k];                               //...........?????						
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
						theta=3.1415926*(j-1)/(Num_Smooth_y-3);              //...........?????
						value[i][j][k]=value[i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
							temp_var[i][j][k];                               //...........?????						
					}
				}
			}
			for(i=1;i<Grid_Num_x-1;i++)
			{
				for(j=Grid_Num_y-Num_Smooth_y;j<Grid_Num_y-1;j++)
				{
					for(k=1;k<Grid_Num_z-1;k++)
					{
						theta=3.1415926*(Grid_Num_y-j-2)/(Num_Smooth_y-3);              //...........?????
						value[i][j][k]=value[i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
							temp_var[i][j][k];                               //...........?????						
					}
				}
			}
		}
		/*
		if (period_y==True)
		{
		;
		}
		*/
		// z-direction
		for(i=1;i<Grid_Num_x-1;i++)
		{
			for(j=1;j<Grid_Num_y-1;j++)
			{
				for(k=1;k<Grid_Num_z-1;k++)
				{
					// theta=2*3.1415926*Z[k]/z_min;              //...........?????
					value[i][j][k]=value[i][j][k]+(1./48.)* \
						temp_var[i][j][k];                        //...........?????
					// (1./48.)* (2.+cos(theta))/3.*temp_var[i][j][k];       
				}
			}
		}
/*
		if(half_z==False)
		{
			for(i=1;i<Grid_Num_x-1;i++)
			{
				for(j=1;j<Grid_Num_y-1;j++)
				{
					for(k=Grid_Num_z-Num_Smooth_z;k<Grid_Num_z-1;k++)
					{
						theta=3.1415926*(Grid_Num_z-k-2)/(Num_Smooth_z-3);         //...........?????
						value[i][j][k]=value[i][j][k]+(1./96.)*(.5*(1+cos(theta)))* \
							temp_var[i][j][k];                                     //...........?????
					}
				}
			}
		}		
*/
	}
	delete []temp_var;
}

void VARIABLE::average(double aver_coeff)
{
	int i,j,k;
	double (*temp_var)[Grid_Num_y][Grid_Num_z]= new double[Grid_Num_x][Grid_Num_y][Grid_Num_z];	
	for(i=1;i<Grid_Num_x-1;i++)
	{
		for(j=1;j<Grid_Num_y-1;j++)
		{
			for(k=1;k<Grid_Num_z-1;k++)
			{
				temp_var[i][j][k]=aver_coeff*value[i][j][k]+ \
					(1-aver_coeff)*( value[i-1][j][k]+value[i+1][j][k]+ \
					value[i][j-1][k]+value[i][j+1][k]+ \
					value[i][j][k-1]+value[i][j][k+1] )/6.;
			}
		}
	}
	for(i=1;i<Grid_Num_x-1;i++)
	{
		for(j=1;j<Grid_Num_y-1;j++)
		{
			for(k=1;k<Grid_Num_z-1;k++)
			{
				value[i][j][k]=temp_var[i][j][k];
			}
		}
	}
	delete []temp_var;
}	
