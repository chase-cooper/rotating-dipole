#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
using namespace std;

#include "header.h"
#include "structs.cpp"
#include "writeout.cpp"

// Grid initial function
Grid3D initGrid() {             // Non-dimensionalized
    Cell emptyCell = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    Grid3D grid(GRIDSIZE,vector<vector<Cell>>(GRIDSIZE,vector<Cell>(GRIDSIZE,emptyCell)));

    // Set grid coordinates and checkerboard pattern
    for (int i = 0; i < GRIDSIZE; i++) {
        for (int j = 0; j < GRIDSIZE; j++) {
            for (int k = 0; k < GRIDSIZE; k++) {
                // Set position of each cell center, units of dipole radii
                grid[i][j][k].x = GRIDSTEP*(i-GRIDRADIUS);
                grid[i][j][k].y = GRIDSTEP*(j-GRIDRADIUS);
                grid[i][j][k].z = GRIDSTEP*(k-GRIDRADIUS);
                // Set color
                grid[i][j][k].color = (i+j+k)%2;
            }
        }
    }
    cout << "Grid initialized" << endl;
    return grid;
}

// Dipole-related functions
DipoleUpdate updateDipole(Cell c,double t) {    // (porbably) Non-dimensionalized
    DipoleUpdate res;
    // Locations of dipole charges, units of dipole radii; time in units of dipole periods
    double xpos = cos(SPEED*t);
    double ypos = sin(SPEED*t);
    double zpos = 0.;
    double xneg = -xpos;
    double yneg = -ypos;
    double zneg = 0.;

    // Cell center displacement from positive charge, units of dipole radii
    double dxpos = c.x - xpos;
    double dypos = c.y - ypos;
    double dzpos = c.z - zpos;
    // Cell center displacement from negative charge, units of dipole radii
    double dxneg = c.x - xneg;
    double dyneg = c.y - yneg;
    double dzneg = c.z - zneg;

    // Calculate each charge as a thin Gaussian distribution
    double Gpos = pow(GRIDSTEP/(sqrt(PI)*SIGMA),3)*exp(-(dxpos*dxpos + dypos*dypos + dzpos*dzpos)/(SIGMA*SIGMA));
    double Gneg = pow(GRIDSTEP/(sqrt(PI)*SIGMA),3)*exp(-(dxneg*dxneg + dyneg*dyneg + dzneg*dzneg)/(SIGMA*SIGMA));

    res.p = Gpos - Gneg;                                // Charge density
    res.Jx = -SPEED*sin(SPEED*t)*(Gpos+Gneg);                 // Current density x-component
    res.Jy = SPEED*cos(SPEED*t)*(Gpos+Gneg);                  // Current density y-component
    return res;
}

Grid3D chargeUpdate(Grid3D grid,double currentTime) {
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                // Get charge density at this cell center
                DipoleUpdate du = updateDipole(grid[i][j][k],currentTime);
                grid[i][j][k].p = du.p;     // probably fine
                // grid[i][j][k].Jx = du.Jx;     
                // grid[i][j][k].Jy = du.Jy;     
            }
        }
    }
    return grid;
}

Grid3D currentUpdate(Grid3D grid,double currentTime) {
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                // Get charge density at this cell center
                DipoleUpdate du = updateDipole(grid[i][j][k],currentTime);
                // grid[i][j][k].p = du.p;     // probably fine
                grid[i][j][k].Jx = du.Jx;     
                grid[i][j][k].Jy = du.Jy;     
            }
        }
    }
    return grid;
}

Grid3D electricFieldRelaxation(Grid3D grid) {   // (likely) finished

    for (int t = 0;t<100;t++) {
        // The below loops iteratively approximate the electric potential over the grid. The cells
        //  are split into red and black squares, and they are calculated separately.

        // First, calculate for all the "red" squares (color=1)
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    grid[i][j][k].oldPhi = grid[i][j][k].phi;
                
                    if (grid[i][j][k].color==0) {continue;}

                    double new_phi = pow(GRIDSTEP,2.)*grid[i][j][k].p;
                    // cout << "GRIDSTEP: " << GRIDSTEP << endl;
                    // cout << "Cell charge density: " << grid[i][j][k].p << endl;
                    if (i != 0) {new_phi += grid[i-1][j][k].phi;}
                    if (i != GRIDSIZE-1) {new_phi += grid[i+1][j][k].phi;}
                    if (j != 0) {new_phi += grid[i][j-1][k].phi;}
                    if (j != GRIDSIZE-1) {new_phi += grid[i][j+1][k].phi;}
                    if (k != 0) {new_phi += grid[i][j][k-1].phi;}
                    if (k != GRIDSIZE-1) {new_phi += grid[i][j][k+1].phi;}
                    grid[i][j][k].phi = new_phi/6.;
                    // cout << "New phi: " << new_phi << endl;
                }
            }
        }
        // Now do the "black" squares (color=0)
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    grid[i][j][k].oldPhi = grid[i][j][k].phi;
                
                    if (grid[i][j][k].color==1) {continue;}

                    double new_phi = pow(GRIDSTEP,2.)*grid[i][j][k].p;
                    if (i != 0) {new_phi += grid[i-1][j][k].phi;}
                    if (i != GRIDSIZE-1) {new_phi += grid[i+1][j][k].phi;}
                    if (j != 0) {new_phi += grid[i][j-1][k].phi;}
                    if (j != GRIDSIZE-1) {new_phi += grid[i][j+1][k].phi;}
                    if (k != 0) {new_phi += grid[i][j][k-1].phi;}
                    if (k != GRIDSIZE-1) {new_phi += grid[i][j][k+1].phi;}
                    grid[i][j][k].phi = new_phi/6.;
                }
            }
        }

        // Update the values of phi using relaxation.
        double omega = 1.85;     // Relaxation factor
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    grid[i][j][k].phi = (1.-omega)*grid[i][j][k].oldPhi + omega*grid[i][j][k].phi;
                }
            }
        }

        // // Calculate residuals (NEEDS WORK)
        double residual;
        double sum_residuals;
        double max_residual = 0;            // holds maximum residual in the grid
        double imax = GRIDSIZE-1;           // shortcut for max grid index
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    double laplacian = 0;
                    if (i != 0) {laplacian += grid[i-1][j][k].phi;}
                    if (j != 0) {laplacian += grid[i][j-1][k].phi;}
                    if (k != 0) {laplacian += grid[i][j][k-1].phi;}
                    if (i != imax) {laplacian += grid[imax][j][k].phi;}
                    if (j != imax) {laplacian += grid[i][imax][k].phi;}
                    if (k != imax) {laplacian += grid[i][j][imax].phi;}
                    laplacian = (laplacian - 6.*grid[i][j][k].phi)/(GRIDSTEP*GRIDSTEP);

                    double residual = -pow(GRIDSTEP,3.)*grid[i][j][k].p - laplacian;
                    residual = residual/pow(grid[i][j][k].p,2.);
                    if (abs(residual) > max_residual) {max_residual = abs(residual);}
                    sum_residuals += residual;
                }
            }
        }
        cout << "Residuals: " << to_string(sum_residuals) << endl;
    }
    // Finally, calculate electric field components for each cell
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                double Ex = grid[i][j][k].phi;
                if (i != GRIDSIZE-1) {Ex -= grid[i+1][j][k].phi;}
                grid[i][j][k].Ex = Ex/(GRIDSTEP);

                double Ey = grid[i][j][k].phi;
                if (j != GRIDSIZE-1) {Ey -= grid[i][j+1][k].phi;}
                grid[i][j][k].Ey = Ey/(GRIDSTEP);

                double Ez = grid[i][j][k].phi;
                if (k != GRIDSIZE-1) {Ez -= grid[i][j][k+1].phi;}
                grid[i][j][k].Ez = Ez/(GRIDSTEP);
            }
        }
    }

    cout << "Charge distributed; initial electric field achieved" << endl;
    return grid;
}

Grid3D electricFieldUpdate(Grid3D grid) {
    for (int i=1;i<GRIDSIZE-1;i++) {
        for (int j=1;j<GRIDSIZE-1;j++) {
            for (int k=1;k<GRIDSIZE-1;k++) {
                // Save previous e-field components
                grid[i][j][k].Ex0 = grid[i][j][k].Ex;
                grid[i][j][k].Ey0 = grid[i][j][k].Ey;
                grid[i][j][k].Ez0 = grid[i][j][k].Ez;

                // New e-field components
                double newEx,newEy,newEz;

                // Ex
                newEx = (grid[i][j][k].Hz - grid[i][j][k].Hy);
                if (j>0) {newEx -= grid[i][j-1][k].Hz;}
                if (k>0) {newEx += grid[i][j][k-1].Hy;}
                grid[i][j][k].Ex += 0.5*(newEx - GRIDSTEP*grid[i][j][k].Jx);
                // Ey
                newEy = (grid[i][j][k].Hx - grid[i][j][k].Hz);
                if (k>0) {newEy -= grid[i][j][k-1].Hx;}
                if (i>0) {newEy += grid[i-1][j][k].Hz;}
                grid[i][j][k].Ey += 0.5*(newEy - GRIDSTEP*grid[i][j][k].Jy);
                // Ez
                newEz = (grid[i][j][k].Hy - grid[i][j][k].Hx);
                if (i>0) {newEz -= grid[i-1][j][k].Hy;}
                if (j>0) {newEz += grid[i][j-1][k].Hx;}
                grid[i][j][k].Ez += 0.5*newEz;

                // E-field constraint
                double econ = 0;
                if (i>0) {econ -= grid[i-1][j][k].Ex;}
                if (i<GRIDSIZE-1) {econ += grid[i+1][j][k].Ex;}
                if (j>0) {econ -= grid[i][j-1][k].Ey;}
                if (j<GRIDSIZE-1) {econ += grid[i][j+1][k].Ey;}
                if (k>0) {econ -= grid[i][j][k-1].Ez;}
                if (k<GRIDSIZE-1) {econ += grid[i][j][k+1].Ez;}
                grid[i][j][k].EConstraint = econ/GRIDSTEP - grid[i][j][k].p;
            }
        }
    }
    return grid;
}

Grid3D magneticFieldUpdate(Grid3D grid) {
    for (int i=1;i<GRIDSIZE-1;i++) {
        for (int j=1;j<GRIDSIZE-1;j++) {
            for (int k=1;k<GRIDSIZE-1;k++) {
                // Save previous h-field components
                grid[i][j][k].Hx0 = grid[i][j][k].Hx;
                grid[i][j][k].Hy0 = grid[i][j][k].Hy;
                grid[i][j][k].Hz0 = grid[i][j][k].Hz;

                double newHx,newHy,newHz;
                // Hx
                newHx = (-grid[i][j][k].Ez + grid[i][j][k].Ey);
                if (j<GRIDSIZE-1) {newHx += grid[i][j+1][k].Ez;}
                if (k<GRIDSIZE-1) {newHx -= grid[i][j][k+1].Ey;}
                grid[i][j][k].Hx -= 0.5*newHx;
                // Hy
                newHy = (-grid[i][j][k].Ex+grid[i][j][k].Ez);
                if (k<GRIDSIZE-1) {newHy += grid[i][j][k+1].Ex;}
                if (i<GRIDSIZE-1) {newHy -= grid[i+1][j][k].Ez;}
                grid[i][j][k].Hy -= 0.5*newHy;
                // Hz
                newHz = (-grid[i][j][k].Ey + grid[i][j][k].Ex);
                if (i<GRIDSIZE-1) {newHz += grid[i+1][j][k].Ey;}
                if (j<GRIDSIZE-1) {newHz -= grid[i][j+1][k].Ex;}
                grid[i][j][k].Hz -= 0.5*newHz;

                // H-field constraint
                double hcon = 0;
                if (i>0) {hcon -= grid[i-1][j][k].Hx;}
                if (i<GRIDSIZE-1) {hcon += grid[i+1][j][k].Hx;}
                if (j>0) {hcon -= grid[i][j-1][k].Hy;}
                if (j<GRIDSIZE-1) {hcon += grid[i][j+1][k].Hy;}
                if (k>0) {hcon -= grid[i][j][k-1].Hz;}
                if (k<GRIDSIZE-1) {hcon += grid[i][j][k+1].Hz;}
                grid[i][j][k].HConstraint = hcon/GRIDSTEP;
            }
        }
    }
    return grid;
}

Grid3D applyABC(Grid3D grid) {

    // Need to be updated all at once, not one by one?
    // double coeff = (TIMESTEP - GRIDSTEP)/(TIMESTEP + GRIDSTEP);
    double coeff = -0.5;
    int imax=GRIDSIZE-1;        // max index (i/j/k) shorthand

    Cell emptyCell = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    Grid3D newgrid(GRIDSIZE,vector<vector<Cell>>(GRIDSIZE,vector<Cell>(GRIDSIZE,emptyCell)));
    // double edgefactor = 1.;

    // update x-faces (not edges or corners)
    for (int j=0; j < GRIDSIZE; j++) {
        for (int k=0; k < GRIDSIZE; k++) {
            // double edgefactor = 1.;         // This factor will change to 0.5 if current cell is on an edge, in which case
            //                                 //  it is double-counted, so update values are halved.
            // if ((j==0)||(k==0)||(j==imax)||(k==imax)) {edgefactor = 0.5;}
            // upper boundary
            newgrid[imax][j][k].Ey += grid[imax-1][j][k].Ey0 + coeff*(grid[imax-1][j][k].Ey0 - grid[imax][j][k].Ey0);
            newgrid[imax][j][k].Ez += grid[imax-1][j][k].Ez0 + coeff*(grid[imax-1][j][k].Ez0 - grid[imax][j][k].Ez0);
            newgrid[imax][j][k].Hy += grid[imax-1][j][k].Hy0 + coeff*(grid[imax-1][j][k].Hy0 - grid[imax][j][k].Hy0);
            newgrid[imax][j][k].Hz += grid[imax-1][j][k].Hz0 + coeff*(grid[imax-1][j][k].Hz0 - grid[imax][j][k].Hz0);
            // upper boundary
            newgrid[0][j][k].Ey += grid[1][j][k].Ey0 + coeff*(grid[1][j][k].Ey0 - grid[0][j][k].Ey0);
            newgrid[0][j][k].Ez += grid[1][j][k].Ez0 + coeff*(grid[1][j][k].Ez0 - grid[0][j][k].Ez0);
            newgrid[0][j][k].Hy += grid[1][j][k].Hy0 + coeff*(grid[1][j][k].Hy0 - grid[0][j][k].Hy0);
            newgrid[0][j][k].Hz += grid[1][j][k].Hz0 + coeff*(grid[1][j][k].Hz0 - grid[0][j][k].Hz0);
        }
    }

    // update y-faces
    for (int i=0; i < GRIDSIZE; i++) {
        for (int k=0; k < GRIDSIZE; k++) {
            // double edgefactor = 1.;         // This factor will change to 0.5 if current cell is on an edge, in which case
            //                                 //  it is double-counted, so update values are halved.
            // if ((i==0)||(k==0)||(i==imax)||(k==imax)) {edgefactor = 0.5;}
            // upper boundary
            newgrid[i][imax][k].Ex += grid[i][imax-1][k].Ex0 + coeff*(grid[i][imax-1][k].Ex0 - grid[i][imax][k].Ex0);
            newgrid[i][imax][k].Ez += grid[i][imax-1][k].Ez0 + coeff*(grid[i][imax-1][k].Ez0 - grid[i][imax][k].Ez0);
            newgrid[i][imax][k].Hx += grid[i][imax-1][k].Hx0 + coeff*(grid[i][imax-1][k].Hx0 - grid[i][imax][k].Hx0);
            newgrid[i][imax][k].Hz += grid[i][imax-1][k].Hz0 + coeff*(grid[i][imax-1][k].Hz0 - grid[i][imax][k].Hz0);
            // lower boundary
            newgrid[i][0][k].Ex += grid[i][1][k].Ex0 + coeff*(grid[i][1][k].Ex0 - grid[i][0][k].Ex0);
            newgrid[i][0][k].Ez += grid[i][1][k].Ez0 + coeff*(grid[i][1][k].Ez0 - grid[i][0][k].Ez0);
            newgrid[i][0][k].Hx += grid[i][1][k].Hx0 + coeff*(grid[i][1][k].Hx0 - grid[i][0][k].Hx0);
            newgrid[i][0][k].Hz += grid[i][1][k].Hz0 + coeff*(grid[i][1][k].Hz0 - grid[i][0][k].Hz0);
        }
    }

    // update z-faces
    for (int i=0; i < GRIDSIZE; i++) {
        for (int j=0; j < GRIDSIZE; j++) {
            // double edgefactor = 1.;         // This factor will change to 0.5 if current cell is on an edge, in which case
            //                                 //  it is double-counted, so update values are halved.
            // if ((i==0)||(j==0)||(i==imax)||(j==imax)) {edgefactor = 0.5;}
            
            // upper boundary
            newgrid[i][j][imax].Ex += grid[i][j][imax-1].Ex0 + coeff*(grid[i][j][imax-1].Ex0 - grid[i][j][imax].Ex0);
            newgrid[i][j][imax].Ey += grid[i][j][imax-1].Ey0 + coeff*(grid[i][j][imax-1].Ey0 - grid[i][j][imax].Ey0);
            newgrid[i][j][imax].Hx += grid[i][j][imax-1].Hx0 + coeff*(grid[i][j][imax-1].Hx0 - grid[i][j][imax].Hx0);
            newgrid[i][j][imax].Hy += grid[i][j][imax-1].Hy0 + coeff*(grid[i][j][imax-1].Hy0 - grid[i][j][imax].Hy0);
            // lower boundary
            newgrid[i][j][0].Ex += grid[i][j][1].Ex0 + coeff*(grid[i][j][1].Ex0 - grid[i][j][0].Ex0);
            newgrid[i][j][0].Ey += grid[i][j][1].Ey0 + coeff*(grid[i][j][1].Ey0 - grid[i][j][0].Ey0);
            newgrid[i][j][0].Hx += grid[i][j][1].Hx0 + coeff*(grid[i][j][1].Hx0 - grid[i][j][0].Hx0);
            newgrid[i][j][0].Hy += grid[i][j][1].Hy0 + coeff*(grid[i][j][1].Hy0 - grid[i][j][0].Hy0);
            
        }
    }

    for (int ind1=0; ind1 < GRIDSIZE; ind1++) {
        for (int ind2=0; ind2 < GRIDSIZE; ind2++) {
            double edgefactor = 1.;         // This factor will change to 0.5 if current cell is on an edge, in which case
                                            //  it is double-counted, so update values are halved.
            // if ((ind1==0)||(ind2==0)||(ind1==imax)||(ind2==imax)) {edgefactor = 0.5;}
            // x-faces
            // upper boundary
            grid[imax][ind1][ind2].Ey = edgefactor * (newgrid[imax][ind1][ind2].Ey);
            grid[imax][ind1][ind2].Ez = edgefactor * (newgrid[imax][ind1][ind2].Ez);
            grid[imax][ind1][ind2].Hy = edgefactor * (newgrid[imax][ind1][ind2].Hy);
            grid[imax][ind1][ind2].Hz = edgefactor * (newgrid[imax][ind1][ind2].Hz);
            // lower boundary
            grid[0][ind1][ind2].Ey = edgefactor * (newgrid[0][ind1][ind2].Ey);
            grid[0][ind1][ind2].Ez = edgefactor * (newgrid[0][ind1][ind2].Ez);
            grid[0][ind1][ind2].Hy = edgefactor * (newgrid[0][ind1][ind2].Hy);
            grid[0][ind1][ind2].Hz = edgefactor * (newgrid[0][ind1][ind2].Hz);
            // y-0
            // upper boundary
            grid[ind1][imax][ind2].Ex = edgefactor * (newgrid[ind1][imax][ind2].Ex);
            grid[ind1][imax][ind2].Ez = edgefactor * (newgrid[ind1][imax][ind2].Ez);
            grid[ind1][imax][ind2].Hx = edgefactor * (newgrid[ind1][imax][ind2].Hx);
            grid[ind1][imax][ind2].Hz = edgefactor * (newgrid[ind1][imax][ind2].Hz);
            // lower boundary
            grid[ind1][0][ind2].Ex = edgefactor * (newgrid[ind1][0][ind2].Ex);
            grid[ind1][0][ind2].Ez = edgefactor * (newgrid[ind1][0][ind2].Ez);
            grid[ind1][0][ind2].Hx = edgefactor * (newgrid[ind1][0][ind2].Hx);
            grid[ind1][0][ind2].Hz = edgefactor * (newgrid[ind1][0][ind2].Hz);
            // z-faces
            // upper boundary
            grid[ind1][ind2][imax].Ex = edgefactor * (newgrid[ind1][ind2][imax].Ex);
            grid[ind1][ind2][imax].Ey = edgefactor * (newgrid[ind1][ind2][imax].Ey);
            grid[ind1][ind2][imax].Hx = edgefactor * (newgrid[ind1][ind2][imax].Hx);
            grid[ind1][ind2][imax].Hy = edgefactor * (newgrid[ind1][ind2][imax].Hy);
            // lower boundary
            grid[ind1][ind2][0].Ex = edgefactor * (newgrid[ind1][ind2][0].Ex);
            grid[ind1][ind2][0].Ey = edgefactor * (newgrid[ind1][ind2][0].Ey);
            grid[ind1][ind2][0].Hx = edgefactor * (newgrid[ind1][ind2][0].Hx);
            grid[ind1][ind2][0].Hy = edgefactor * (newgrid[ind1][ind2][0].Hy);
        }
    }
    

    return grid;
}

// Functions for evaluating the system

double GaussLaw(Grid3D grid) {
    int imax = GRIDSIZE-1;
    double totalFlux;
    for (int ind1=0;ind1<GRIDSIZE;ind1++) {
        for (int ind2=0;ind2<GRIDSIZE;ind2++) {
            totalFlux += grid[imax][ind1][ind2].Ex - grid[0][ind1][ind2].Ex;
            totalFlux += grid[ind1][imax][ind2].Ey - grid[ind1][0][ind2].Ey;
            totalFlux += grid[ind1][ind2][imax].Ez - grid[ind1][ind2][0].Ez;
        }
    }
    return totalFlux;
}

double GridEnergy(Grid3D grid) {
    int imax = GRIDSIZE-1;
    double totalEnergy;
    
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                totalEnergy += 0.5*sqrt(pow(grid[i][j][k].Ex,2.) + pow(grid[i][j][k].Ey,2.) + pow(grid[i][j][k].Ey,2.));
                totalEnergy += 0.5*sqrt(pow(grid[i][j][k].Hx,2.) + pow(grid[i][j][k].Hy,2.) + pow(grid[i][j][k].Hz,2.));
            }
        }
    }

    return totalEnergy * pow(GRIDSTEP,3.);
}

vector<double> PoyntingVector(Grid3D grid,int i,int j,int k) {
    double xcomp = (grid[i][j][k].Ey * grid[i][j][k].Hz) - (grid[i][j][k].Ez * grid[i][j][k].Hy);
    double ycomp = (grid[i][j][k].Ez * grid[i][j][k].Hx) - (grid[i][j][k].Ex * grid[i][j][k].Hz);
    double zcomp = (grid[i][j][k].Ex * grid[i][j][k].Hy) - (grid[i][j][k].Ey * grid[i][j][k].Hx);
    // double largest = max(max(abs(xcomp),abs(ycomp)),abs(zcomp));

    vector<double> res = {xcomp,ycomp,zcomp};
    return res;
}

// Main
int main() {
    // // Initialize grid
    double currentTime = 0.;
    Grid3D grid = initGrid();
    grid = chargeUpdate(grid,currentTime);
    grid = currentUpdate(grid,currentTime);
    grid = electricFieldRelaxation(grid);

    cout << setprecision(10);
    cout << "Gauss' Law initial estimate: q = " << GaussLaw(grid) << endl;

    writeOut(grid,0);
    cout << "Initial conditions written to file. Beginning loop..." << endl;

    ofstream summary("outputs/summary.dat");
    summary << "Step \tGauss' Law \t\tEnergy\t\t\tPoynting Vector (Scaled)" << endl;

    for (int count=1;count <= LOOPCOUNT;count++) {       // count is just the loop counter variable
        // Magnetic field half-step
        currentTime += 0.5*TIMESTEP;

        grid = chargeUpdate(grid,currentTime);
        grid = currentUpdate(grid,currentTime);
        grid = magneticFieldUpdate(grid);

        // Electric field half-step
        currentTime += 0.5*TIMESTEP;
        grid = electricFieldUpdate(grid);
        grid = applyABC(grid);

        writeOut(grid,count);

        double g = GaussLaw(grid);
        double e = GridEnergy(grid);
        cout << "Finished loop " << to_string(count) << " of " << to_string(LOOPCOUNT);
        cout << "; Gauss' Law estimate q = " << g;
        cout << "; Total energy E = " << e << endl;
        summary << count << "\t\t" << g << "\t\t" << e << "\t\t\t";
    }
    

    return 0;
    
}
