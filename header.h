#define PI 3.1415927
#define C 1 //299792458                 // [m]
#define EPS0 1 //8.854188e-12           // [C^2 kg^-1 m^-3 s^2]
#define MU0 1 //1.256637e-6             // [kg m s^-2 A^2]

#define DIPOLE 5e-10                // dipole radius [m]
#define OMEGA 1e-14                 // dipole spin period [s]
#define CHARGE 1.602e-19            // elementary charge

#define SIGMA 0.5                   // Charge Gaussian std.dev [dipole radii]
#define GRIDSTEP 0.1              // grid spacing [dipole radii]
#define TIMESTEP GRIDSTEP/30        // timestep [dipole periods]
// #define GRIDSTEP 0.1*DIPOLE
// #define TIMESTEP OMEGA/30
#define GRIDRADIUS 15               // a "grid radius", grid will be (2n+1)^3; recommended 20 or more
#define LOOPCOUNT 50                 // number of loops to run

#define GRIDSIZE (2*GRIDRADIUS+1)   // actual grid size per side
