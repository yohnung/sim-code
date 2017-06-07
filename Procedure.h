// Important_Procedure.h
// ???(Shoud describe more detainly about following function)

// Setting mesh grid
void set_mesh();

// Initialising variables and setting mesh-grid
void initialize(VARIABLE *, BASIC_VARIABLE &);

// initialize the problem as a harris-current-shet probelm
void harris_current_initia(VARIABLE *, BASIC_VARIABLE &);

// add fluctuation at z=up and down boundary according to <Hurricane, PoP, 1995>, \delta\psi=fluctuation*cos(k_z*z) 
void fluc_at_bndry(VARIABLE *, double, double);

// add fluctuation at neutral-line according to <karimabadi, JGR, 2004>, 
void fluc_at_neutral_line(VARIABLE *, double, double, double);

// Calculate current
void cal_current(VARIABLE *, VARIABLE *, Type=Complete); 

// Calculate pressure
void cal_pressure(BASIC_VARIABLE &, VARIABLE *, Type=Complete);

// Setting eta
void set_eta(BASIC_VARIABLE &, VARIABLE *, VARIABLE *, double, Type=Complete, Type=Uniform);

// Calculating fluxes
void cal_flux(BASIC_VARIABLE [][3], VARIABLE *, VARIABLE *, BASIC_VARIABLE &, BASIC_VARIABLE &, Type=Complete);

// Extracting electric field from flux
void ext_from_flux(BASIC_VARIABLE *, BASIC_VARIABLE [][3]);

// Setting time-interval
double set_dt(VARIABLE *, BASIC_VARIABLE &, VARIABLE *, BASIC_VARIABLE &, double, double);

// Step on variables
void step_on(VARIABLE *, BASIC_VARIABLE [][3], double, double);

void smooth(VARIABLE *, double, int);



// General_Procedure.h

void make_pressure_positive(BASIC_VARIABLE &, double, Type);

// Update variables using 2-order Lax-Wendroff method
void exclude_soucrce_half_update(VARIABLE *, BASIC_VARIABLE [][3], double, Order);

// Update variables from source term
void source_update(VARIABLE *, double);

void boudary_set(VARIABLE &, Symmetry_Type, Symmetry_Type);

// Copying action
void copy(VARIABLE *, VARIABLE *);

void smooth_xyz(double [][Grid_Num_x][Grid_Num_y][Grid_Num_z], int);

void average(double [][Grid_Num_x][Grid_Num_y][Grid_Num_z],double);

void average_B(double [][Grid_Num_x][Grid_Num_y][Grid_Num_z],double);


