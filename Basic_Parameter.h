//Basic_Parameter.h
// change basic set to simulate a harris-current-layer problem

// Grid Number
const int Grid_Num_x=11;     // VS run, maximum value is 13 
const int Grid_Num_y=9;
const int Grid_Num_z=13;     // set Grid_Num_z=4, and make all 
			     // variables independent on y, we can simulate 2D situation

// Controlling and Logical Parameter
//const Logic lstrt=False;           
const Logic uniform_x=True;
const Logic uniform_y=True;       // Unifrom mesh or non-uniform mesh
const Logic uniform_z=True;       // Fortran doesn't use it!
const Logic period_y=True; // False;        // False;       // Periodic Condition in Y-direction
const Logic half_x=True;  //False;
// const Logic half_y=False;         // Symmetric or antisymmetric simulation
const Logic half_z=False;  //False;

// Spatial Range
const double x_min=1.;
const double x_max=11.;
const double y_min=-5.;
const double y_max=5.;
const double z_min=-5;
const double z_max=5.;

// Physical Parameter
const double phy_gamma=1.66667; //5./3. means physics_gamma;
const double beta_m=0.1;        // 0.01
const double rho_m_0=7.33;            // 1.
const double rho_s_0=7.33;            // 1.
const double width_rho=1.;         // 1.
const double B_m_0=13.;             // 1.
const double B_s_0=11.;             // 1.
const double v_0=7.;               // 0.
const double vyi_0=2.30;             // 0.
const double di=0.3;                // 0.  what's its meaning???????? used in current amendation and V_cross_B in flux calculation

const int Num_Smooth_x=2*Grid_Num_x/3;
const int Num_Smooth_y=3*Grid_Num_y/4;
const int Num_Smooth_z=2*Grid_Num_z/3;
