//main.cpp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
#include "Basic_Parameter.h"
#include "Runtime_Diagnostic_Parameter.h"
#include "Variables_Definition.h"
#include "Procedure.h"


/*global mesh   X[0]_____X_interval[0]_____X[1]: start*/ 
double  X[Grid_Num_x],Y[Grid_Num_y],Z[Grid_Num_z];
double  X_interval[Grid_Num_x],Y_interval[Grid_Num_y],Z_interval[Grid_Num_z];    // X_interval[Gridn_Num_?-1]=X_interval[Gridn_Num_?-2]
/*global mesh   X[0]_____X_interval[0]_____X[1]: end*/


/* global parameter: start*/
	BASIC_VARIABLE Elec_field[3], p, eta; //                         // Define basic variables, eta for conductivity.
	VARIABLE var[8];                               // Define main variables, and var_int for variable_intermediate
	VARIABLE current[3]; // vorticity[3];          // Define more varibles.
	BASIC_VARIABLE flux[8][3];
	double system_time=0;
	double dt=0.1;                                     // Define time.
	int nstep=0;
/* global paremeter: end*/

int main()
{
	ofstream out[11];	
	ofstream timeout("step_to_time.txt");         // file="stepnm"
	out[0].open("rho.dat");
	out[1].open("rhoVx.dat");out[2].open("rhoVy.dat");out[3].open("rhoVz.dat"); 
	out[4].open("Bx.dat");out[5].open("By.dat");out[6].open("Bz.dat");
	out[7].open("E.dat");
	out[8].open("Electric_Field_x.dat");
	out[9].open("Electric_Field_y.dat");
	out[10].open("Electric_Field_z.dat");	

	int i;                                         // cycle variable

	set_mesh();
	//initialize(var, p);	                           // Initializing variables and pressure 
	harris_current_initia(var,p);
	for (i=0;i<8;i++)                              // and out put
		var[i].record(out[i]);

	/* time step on variables:start */
	for (nstep=nstart;nstep<nend;nstep++)
	{
		cal_current(current, var);                          // Calculating current from Magnetic Field.
		set_eta(eta, var, current, system_time);            // Setting space dependent conductivity. (Can make it depend on current)
		cal_pressure(p, var);                               // Calculating pressure from various kinds of energy.	
		dt=set_dt(var, eta, current, p, system_time, dt);       // Settiing appropriate time-interval from main variables, conductivity and pressure and so on. This statement change dt only.
		if (nstep==10)
			add_fluc(var);
		cal_flux(flux, var, current, p, eta);               // Calculating flux from variables, current and pressure.	 
		ext_from_flux(Elec_field, flux);                     // extractig electric field from flux
		step_on(var, flux, system_time, dt);                // Main procedure to time step on variables from Flux explicitly and from Source implicitly.
		smooth(var,system_time, nstep);                     // ?????????? Havn't understand yet ???????????
		system_time=system_time+dt;		
		cout<<setw(8)<<nstep<<setw(15)<<\
			"time = "<<setiosflags(ios::scientific)<<setprecision(15)<<system_time\
			<<setw(14)<<"dt = "<<dt<<endl;
		timeout<<"nstep = "<<nstep<<endl<<" ,"<<setw(20)<<"and time = "<<setiosflags(ios::scientific)\
			<<setprecision(15)<<system_time<<" ,"<<setw(14)<<"dt = "<<dt<<endl;	
		if ((nstep+1)%12==0)
		{
			for (i=0;i<8;i++)
				var[i].record(out[i]);
			for (i=0;i<3;i++)
				Elec_field[i].record(out[i+8]);
		}
	}
	timeout.close();
	for(i=0;i<11;i++)
		out[i].close();
	/* time step on variables:end */
}
