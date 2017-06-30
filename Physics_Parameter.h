// Physics_Parameter.h
// Used in Important_Procedure.cpp

/* Physical Parameter: start */
const double phy_gamma=5./3.;      // 5./3.=1.66667 means physics_gamma;
const double di=1.;                  // a switch of Hall-effect
const double magnetic_Renolds_Number=1./200;
// Magneto-pause parameter, mainly used in 'initialize()'
const double beta_m=0.01;        // 0.01
const double rho_m_0=1.;            // 1.
const double rho_s_0=1.;            // 1.
const double width_rho=1.;         // 1.
const double B_m_0=1.;             // 1.
const double B_s_0=1.;             // 1.
const double v_0=1.;               // 0.    initialize rhoVy
const double vyi_0=1.;             // 0.    amend rhoVy in magtopause problem
// Harris current parameter, mainly used in 'harris_initia()'
const double rho_infinity=0.2;
const double Balance_coefficient=.5;
const double normalized_lambda=0.5;
// shear flow
const double max_shear_velocity = 0;
const double shear_length = 1.;
const double shear_location = 0.;
/* Physical Parameter: end */

/* fluctuation parameter: start */
const Position position_fluc=Boundary;
const double fluctuation=0.1*di;
const double fluc_kx=Pi/(x_max-x_min);
const double fluc_kz=2*Pi/(z_max-z_min);
/* fluctuation parameter: end */
