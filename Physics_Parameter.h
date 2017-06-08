// Physics_Parameter.h
// Used in Important_Procedure.cpp

/* Physical Parameter: start */
const double phy_gamma=1.66667;      // 5./3. means physics_gamma;
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
const double Balance_coefficient=1.;
const double normalized_lambda=0.5;
/* Physical Parameter: end */

/* fluctuation parameter: start */
const Position position_fluc=Boundary;
const double fluctuation=0.2;
const double fluc_kx=4*3.1415926/x_max;
const double fluc_kz=2*3.1415926/z_max;
/* fluctuation parameter: end */
