// Grid Number
const int Grid_Num_x=7;
const int Grid_Num_y=7;
const int Grid_Num_z=7;

// Controlling and Logical Parameter
//const Logic lstrt=False;           
const Logic uniform_x=True;
const Logic uniform_y=True;       // Unifrom mesh or non-uniform mesh
const Logic uniform_z=True;
const Logic period_y=False;       // Periodic Condition in Y-direction
const Logic half_x=False;
// const Logic half_y=False;         // Symmetric or antisymmetric simulation
const Logic half_z=False;

// Spatial Range
const double x_min=1.;
const double x_max=11.;
const double y_min=-5.;
const double y_max=5.;
const double z_min=-5;
const double z_max=5.;

// Physical Parameter
const double gamma=5./3.;
const double beta_m=0.01;
const double rho_m_0=1.;
const double rho_s_0=1.;
const double width_rho=1.;
const double B_m_0=1.;
const double B_s_0=1.;
const double v_0=0.;
//const double v_y_i_0=0.;

const int Num_Smooth_x=2*Grid_Num_x/3;
const int Num_Smooth_y=3*Grid_Num_y/4;
const int Num_Smooth_z=2*Grid_Num_z/3;