// Variable_Definition.h

//#include "Type_and_Macro_Definition.h"
//#include "Basic_Parameter.h"
//#include "Runtime_Diagnostic_Parameter.h"
//#include <iostream>
//#include <fstream>
//#include <iomanip>
//using namespace std;

class BASIC_VARIABLE
{
public:
	double value[Grid_Num_x][Grid_Num_y][Grid_Num_z];
	void record(ofstream &);
};

class VARIABLE: public BASIC_VARIABLE
{
public:	
	void boundary_set(Symmetry_Type, Symmetry_Type);	
	void smooth_xyz(int times);
	void average(double);
};
