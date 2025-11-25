#define PI 3.1415927

#define SIGMA 0.5                   // Charge Gaussian std.dev [dipole radii]

#define GRIDRADIUS 25               // a "grid radius", grid will be (2n+1)^3; recommended 20 or more
#define GRIDSTEP (2.*PI)/GRIDRADIUS // grid spacing [dipole radii]
#define TIMESTEP GRIDSTEP/2         // timestep [dipole periods]
// #define GRIDSTEP 0.1*DIPOLE
// #define TIMESTEP OMEGA/30
#define LOOPCOUNT 100                 // number of loops to run

#define GRIDSIZE (2*GRIDRADIUS+1)   // actual grid size per side
