// Physics_Parameter.h
// Used in Important_Procedure.cpp

/* Physical Parameter: start */
const double iso_therm_coeff = 1e-4;   // when gamma-1 is lower than 1e-5, it's isothermal process
const double phy_gamma=5./3.;      // 5./3.=1.66667 means physics_gamma;
const double di=1.;                  // a switch of Hall-effect
const double Lundquist_Number=1./200;
const double localized_Resistivity = 0.005;
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
const double velocity_boundary = -3;
const double shear_length = 2.;
const double shear_location = 5.;
/* Physical Parameter: end */

/* fluctuation parameter: start */
const Position position_fluc=Boundary;
const double fluctuation_mag= -0.1;     // this will generates an island whose width is 1a (current layer's width)
const double fluctuation_velocity= -0.1;
const double fluc_B_kx=Pi/(x_max-x_min);
const double fluc_B_kz=2*Pi/(z_max-z_min);
const double fluc_V_kx= Pi / (x_max - x_min);
const double fluc_V_kz = 2 * Pi / (z_max - z_min);
/* fluctuation parameter: end */
