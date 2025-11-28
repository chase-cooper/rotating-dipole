#define PI 3.1415927

#define SIGMA 0.5                   // Charge Gaussian std.dev [dipole radii]
#define OMEGA 0.2                   // Dipole angular period scale. Not representative of anything, jsut
                                    //  adjusts the speed of the dipole

#define GRIDRADIUS 40               // a "grid radius", grid will be (2n+1)^3; recommended 20 or more
#define GRIDSTEP (2.*PI)/(OMEGA*GRIDRADIUS) // grid spacing [dipole radii]
#define TIMESTEP GRIDSTEP/2         // timestep [dipole periods]
// #define GRIDSTEP 0.1*DIPOLE
// #define TIMESTEP OMEGA/30
#define LOOPCOUNT 400                 // number of loops to run

#define GRIDSIZE (2*GRIDRADIUS+1)   // actual grid size per side
