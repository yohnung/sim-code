// Important_Procedure.h
// ???(Shoud describe more detainly about following function)

// Setting mesh grid
void set_mesh();

// Initialising variables and setting mesh-grid
void initialize(VARIABLE *, BASIC_VARIABLE &);

// initialize the problem as a harris-current-shet probelm
void harris_current_initia(VARIABLE *, BASIC_VARIABLE &);

// write to files
void write_out(VARIABLE *, int, double);
// read from files
void read_in(VARIABLE *, int &, double&);
// Extracting electric field from flux
void ext_from_var(BASIC_VARIABLE *, BASIC_VARIABLE *, VARIABLE *, BASIC_VARIABLE &);

// add fluctuation mainly to B variables, specificly refer to fluc_at_bndry/neutral_line
void add_fluc(VARIABLE *);

// Calculate current
void cal_current(VARIABLE *, VARIABLE *, Type=Complete); 

// Calculate pressure
void cal_pressure(BASIC_VARIABLE &, VARIABLE *, Type=Complete);

// Setting eta
void set_eta(BASIC_VARIABLE &, VARIABLE *, VARIABLE *, double, Type=Complete, Type=Uniform);

// accumulate energy form V, B and P
void accumulate_eng(double [][Grid_Num_y][Grid_Num_z], VARIABLE *, BASIC_VARIABLE &, Type=Complete);

// Calculating fluxes
void cal_flux(BASIC_VARIABLE [][3], VARIABLE *, VARIABLE *, BASIC_VARIABLE &, BASIC_VARIABLE &, Type=Complete);

// Setting time-interval
double set_dt(VARIABLE *, BASIC_VARIABLE &, VARIABLE *, BASIC_VARIABLE &, double, double);

// Step on variables
void step_on(VARIABLE *, VARIABLE *, BASIC_VARIABLE &, BASIC_VARIABLE &, double time, double time_interv);

void smooth(VARIABLE *, double);

void record(ofstream &, int, int, double, VARIABLE *, \
	BASIC_VARIABLE &, VARIABLE &, BASIC_VARIABLE &);

// General_Procedure.h
// open and close writting files
void open_var_files(ofstream *);

void close_var_files(ofstream *);

// open and close reading files
void open_var_files(ifstream *);

void close_var_files(ifstream *);

// Update variables using 2-order Lax-Wendroff method
void exclude_soucrce_half_update(VARIABLE *, BASIC_VARIABLE [][3], double, Order);

// Update variables from source term
void source_update(VARIABLE *, double);

void num2str(char *, int);
