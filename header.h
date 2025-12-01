#define PI 3.1415927

#define LOOPCOUNT 500                          // number of loops to run
#define SPEED 1.0                               // Relative dipole spin velocity
#define GRIDRADIUS 20                           // a "grid radius", grid will be (2n+1)^3; recommended 20 or more
#define GRIDSIZE (2*GRIDRADIUS+1)               // actual grid size per side
#define GRIDSTEP 4.*(2.*PI)/(SPEED*GRIDRADIUS)  // grid spacing [dipole radii]. The leading term dictates the desired radius in
                                                //      propogation wavelengths
#define TIMESTEP GRIDSTEP/2                     // timestep [dipole periods]
#define SIGMA 1.0                               // Charge Gaussian std.dev [dipole radii]

