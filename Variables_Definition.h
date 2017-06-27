// Variable_Definition.h


class BASIC_VARIABLE
{
public:
	double value[Grid_Num_x][Grid_Num_y][Grid_Num_z];
	void record(ofstream &);
	void fill(ifstream &);
};

class VARIABLE: public BASIC_VARIABLE
{
public:	
	void boundary_set(Symmetry_Type, Symmetry_Type);
	void smooth_xyz(int times);
	void average(double);
};
