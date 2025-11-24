#define PI 3.1415927

#define SIGMA 0.5                   // Charge Gaussian std.dev [dipole radii]
#define GRIDSTEP 0.2              // grid spacing [dipole radii]
#define TIMESTEP GRIDSTEP/2        // timestep [dipole periods]
// #define GRIDSTEP 0.1*DIPOLE
// #define TIMESTEP OMEGA/30
#define GRIDRADIUS 8               // a "grid radius", grid will be (2n+1)^3; recommended 20 or more
#define LOOPCOUNT 50                 // number of loops to run

#define GRIDSIZE (2*GRIDRADIUS+1)   // actual grid size per side
