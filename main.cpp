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
	ofstream timeout("step_to_time.dat");
	double system_time=0;
	double dt=0.1;                                     // Define time.
	int nstep=0;
/* global paremeter: end*/

int main()
{
	int run_num=0, record_step;

	set_mesh();
	if (continue_from_files==False)
	{
		//initialize(var, p);	                           // Initializing variables and pressure 
		harris_current_initia(var,p);
	}
	else
		read_in(var, nstart, system_time);

	/* time step on variables:start */
	for (nstep=nstart;nstep<nstart+nend;nstep++)
	{	 
		cal_current(current, var);                          // Calculating current from Magnetic Field.
		set_eta(eta, var, current, system_time);            // Setting space dependent conductivity. (Can make it depend on current)
		cal_pressure(p, var);                               // Calculating pressure from various kinds of energy.	
		ext_from_var(Elec_field, var, current, eta);         // 重新些一个从var计算Elec的函数                 // extractig electric field from flux
		
		if(system_time==0)
		{
			record_step=nstep+1;
			record(timeout, run_num, record_step, system_time, var, p, current[1], Elec_field[1]);
			run_num+=1;
		}

		dt=set_dt(var, eta, current, p, system_time, dt);       // Settiing appropriate time-interval from main variables, conductivity and pressure and so on. This statement change dt only.
		if (nstep+1==10)
			add_fluc(var);
		step_on(var, current, p, eta, system_time, dt);                // Main procedure to time step on variables from Flux explicitly and from Source implicitly.
		smooth(var,system_time);                     // ?????????? Havn't understand yet ???????????
		system_time=system_time+dt;		
		cout<<setw(8)<<nstep<<setw(15)<<\
			"time = "<<setiosflags(ios::scientific)<<setprecision(16)<<system_time\
			<<setw(14)<<"dt = "<<dt<<endl;
		if(abs(system_time-run_num*average_time_interval)<=dt &&\
			nstep+1-record_step>1)
		{
			record_step=nstep+1;
			record(timeout, run_num, record_step, system_time, var, p, current[1], Elec_field[1]);
			run_num+=1;
			write_out(var, record_step, system_time);
		}
	}
	/* time step on variables:end */
}
