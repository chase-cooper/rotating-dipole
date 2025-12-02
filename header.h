#define PI 3.1415927

#define LOOPCOUNT 300                           // Number of loops to run
#define SPEED 1.0                               // Relative dipole spin velocity
#define GRIDRADIUS 40                           // A "grid radius", grid will be (2n+1)^3; recommended 20 or more
#define SIGMA 1.0                               // Charge Gaussian standard deviation

#define GRIDSIZE (2*GRIDRADIUS+1)               // Actual grid size per side
#define IMAX GRIDSIZE-1                         // The maximum index of a cell along one axis. We define this because
                                                //      it gets used A LOT for multiple functions.
#define GRIDSTEP 8.*(2.*PI)/(SPEED*GRIDRADIUS)  // Grid spacing [dipole radii]. The leading term dictates the desired radius in
                                                //      propogation wavelengths
#define TIMESTEP GRIDSTEP/2                     // Timestep. We arrive at this relation as a simple solution to the
                                                //      stability condition given in equation 31.

