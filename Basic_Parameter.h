// Basic_Parameter.h
// Usd in all *.cpp files 

/* Macro and Variable-Type definition: start */
// Macro definition
#define Pi 3.141592653589793
// Logic value
#define True 1
#define False 0
// Type value
#define Incomplete 1      // Need to minus 1 from max grid number, meaning every cycle starts at 0, ends at Grid_Num-1.
#define Complete 0        // Need not to minus 1.
// Order value
#define First 1
#define Second 2
// Symetry_Type value
#define Positive 1
#define Negative -1
// Position value
#define upperBoundary 'u'
#define downBoundary 'd'
#define x_Boundary 'x'
#define leftBoundary 'l'
#define rightBoundary 'r'
#define z_Boundary 'z'
#define Neutral_Line 'n'
// variable-type definition
typedef int Logic;        //Logic=True or False
typedef int Type;         // Type=Complete or Incomplete
typedef int Order;        // Order=First or Second
typedef int Symmetry_Type;// Symmetry_Type=Axial or Dot
typedef char Position;
/* Macro and Variable-Type definition: start */


/* Grid Number: start. 
   In 3D simulation Grid Number should be 4n+1.  
   In 2D simulation Grid_Num_y=7 is best, with others 4n+1.   
   4n+1 is for the sake of 'record()', defined in Variables_Definition.cpp
   7 is for the consideration of 'smooth_xyz()' function.                    */
const int Grid_Num_x=41; 
const int Grid_Num_y=7;
const int Grid_Num_z=41; 
const int num_out=20;
/* Grid Number: end */			    

/* Controlling and Logical Parameter:start */
// const Logic Sim_2D=True;        // Not used now! Sim_2D==True, 2D simlulation; otherwise, 3D simulation. 
// const Logic lstrt=False;        // Not used now!
const Logic uniform_x=True;
const Logic uniform_y=True;       // Unifrom mesh or non-uniform mesh
const Logic uniform_z=True;       // Fortran doesn't use it!
const Logic half_x=False;  
// const Logic half_y=False;       // Not used now! Symmetric or antisymmetric simulation. 
const Logic half_z=False;
const Logic period_y=True; // False;     // Periodic Condition in Y-direction
const Logic period_z=(1-half_z)*False;       // if half_z==True, period_z doesn't work
/* Controlling and Logical Parameter: end */

/* Spatial Range: start */
const double x_min=-5.;
const double x_max=5.;
const double y_min=1.;
const double y_max=7.;
const double z_min=-5.;
const double z_max=5.;
/* Spatial Range: end */

/* smooth parameter: start */
const int Num_Smooth_x=2*Grid_Num_x/3;
const int Num_Smooth_y=3*Grid_Num_y/4;
const int Num_Smooth_z=2*Grid_Num_z/3;
/* smooth parameter: end */
