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
    double xpos = cos(2.*PI*t);
    double ypos = sin(2.*PI*t);
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
    res.Jx = -sin(2.*PI*t)*(Gpos+Gneg);                 // Current density x-component
    res.Jy = cos(2.*PI*t)*(Gpos+Gneg);                  // Current density y-component
    return res;
}

Grid3D electricFieldRelaxation(Grid3D grid) {   // (likely) finished

    for (int t = 0;t<500;t++) {
        // The below loops iteratively approximate the electric potential over the grid. The cells
        //  are split into red and black squares, and they are calculated separately.

        // First, calculate for all the "red" squares (color=1)
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    grid[i][j][k].oldPhi = grid[i][j][k].phi;
                
                    if (grid[i][j][k].color==0) {continue;}

                    double new_phi = pow(GRIDSTEP,2.)*grid[i][j][k].p/EPS0;
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
        // Now do the "black" squares (color=0)
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    grid[i][j][k].oldPhi = grid[i][j][k].phi;
                
                    if (grid[i][j][k].color==1) {continue;}

                    double new_phi = pow(GRIDSTEP,2.)*grid[i][j][k].p/EPS0;
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

                    double residual = -pow(GRIDSTEP,3.)*grid[i][j][k].p/EPS0 - laplacian;
                    residual = residual/pow(grid[i][j][k].p/EPS0,2.);
                    if (abs(residual) > max_residual) {max_residual = abs(residual);}
                    sum_residuals += residual;
                }
            }
        }
        // cout << to_string(sum_residuals) << endl;
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

Grid3D chargeCurrentUpdate(Grid3D grid,double currentTime) {
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                // Get charge density at this cell center
                DipoleUpdate du = updateDipole(grid[i][j][k],currentTime);
                grid[i][j][k].p = du.p;     // probably fine
                grid[i][j][k].Jx = du.Jx;     
                grid[i][j][k].Jy = du.Jy;     
            }
        }
    }
    return grid;
}

Grid3D electricFieldUpdate(Grid3D grid) {
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
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
                grid[i][j][k].Ex += TIMESTEP*(newEx/GRIDSTEP - grid[i][j][k].Jx)/EPS0;
                // Ey
                newEy = (grid[i][j][k].Hx - grid[i][j][k].Hz);
                if (k>0) {newEy -= grid[i][j][k-1].Hx;}
                if (i>0) {newEy += grid[i-1][j][k].Hz;}
                grid[i][j][k].Ey += TIMESTEP*(newEy/GRIDSTEP - grid[i][j][k].Jy)/EPS0;
                // Ez
                newEz = (grid[i][j][k].Hy - grid[i][j][k].Hx);
                if (i>0) {newEz -= grid[i-1][j][k].Hy;}
                if (j>0) {newEz += grid[i][j-1][k].Hx;}
                grid[i][j][k].Ez += TIMESTEP*newEz/EPS0/GRIDSTEP;

                // E-field constraint
                double econ = 0;
                if (i>0) {econ -= grid[i-1][j][k].Ex;}
                if (i<GRIDSIZE-1) {econ += grid[i+1][j][k].Ex;}
                if (j>0) {econ -= grid[i][j-1][k].Ey;}
                if (j<GRIDSIZE-1) {econ += grid[i][j+1][k].Ey;}
                if (k>0) {econ -= grid[i][j][k-1].Ez;}
                if (k<GRIDSIZE-1) {econ += grid[i][j][k+1].Ez;}
                grid[i][j][k].EConstraint = econ/GRIDSTEP - grid[i][j][k].p/EPS0;
            }
        }
    }
    return grid;
}

Grid3D magneticFieldUpdate(Grid3D grid) {
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                // Save previous h-field components
                grid[i][j][k].Hx0 = grid[i][j][k].Hx;
                grid[i][j][k].Hy0 = grid[i][j][k].Hy;
                grid[i][j][k].Hz0 = grid[i][j][k].Hz;

                double newHx,newHy,newHz;
                // Hx
                newHx = (-grid[i][j][k].Ez + grid[i][j][k].Ey);
                if (j<GRIDSIZE-1) {newHx += grid[i][j+1][k].Ez;}
                if (k<GRIDSIZE-1) {newHx -= grid[i][j][k+1].Ey;}
                grid[i][j][k].Hx -= (TIMESTEP)*newHx/MU0/(GRIDSTEP);
                // Hy
                newHy = (-grid[i][j][k].Ex+grid[i][j][k].Ez);
                if (k<GRIDSIZE-1) {newHy += grid[i][j][k+1].Ex;}
                if (i<GRIDSIZE-1) {newHy -= grid[i+1][j][k].Ez;}
                grid[i][j][k].Hy -= (TIMESTEP)*newHy/MU0/(GRIDSTEP);
                // Hz
                newHz = (-grid[i][j][k].Ey + grid[i][j][k].Ex);
                if (i<GRIDSIZE-1) {newHz += grid[i+1][j][k].Ey;}
                if (j<GRIDSIZE-1) {newHz -= grid[i][j+1][k].Ex;}
                grid[i][j][k].Hz -= (TIMESTEP)*newHz/MU0/(GRIDSTEP);
            }
        }
    }
    return grid;
}

Grid3D applyABC(Grid3D grid) {

    // Need to be updated all at once, not one by one?

    double coeff = (C*TIMESTEP - GRIDSTEP)/(C*TIMESTEP + GRIDSTEP);
    int imax=GRIDSIZE-1;        // max index (i/j/k) shorthand

    Cell emptyCell = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    Grid3D newgrid(GRIDSIZE,vector<vector<Cell>>(GRIDSIZE,vector<Cell>(GRIDSIZE,emptyCell)));

    // update x-faces
    for (int j=0; j < GRIDSIZE; j++) {
        for (int k=0; k < GRIDSIZE; k++) {
            // upper boundary
            grid[imax][j][k].Ey = grid[imax-1][j][k].Ey0 + coeff*(grid[imax-1][j][k].Ey - grid[imax][j][k].Ey0);
            grid[imax][j][k].Ez = grid[imax-1][j][k].Ez0 + coeff*(grid[imax-1][j][k].Ez - grid[imax][j][k].Ez0);
            grid[imax][j][k].Hy = grid[imax-1][j][k].Hy0 + coeff*(grid[imax-1][j][k].Hy - grid[imax][j][k].Hy0);
            grid[imax][j][k].Hz = grid[imax-1][j][k].Hz0 + coeff*(grid[imax-1][j][k].Hz - grid[imax][j][k].Hz0);

            // lower boundary
            grid[0][j][k].Ey = grid[1][j][k].Ey0 + coeff*(grid[1][j][k].Ey - grid[0][j][k].Ey0);
            grid[0][j][k].Ez = grid[1][j][k].Ez0 + coeff*(grid[1][j][k].Ez - grid[0][j][k].Ez0);
            grid[0][j][k].Hy = grid[1][j][k].Hy0 + coeff*(grid[1][j][k].Hy - grid[0][j][k].Hy0);
            grid[0][j][k].Hz = grid[1][j][k].Hz0 + coeff*(grid[1][j][k].Hz - grid[0][j][k].Hz0);
        }
    }

    return grid;
}

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

// Main
int main() {
    // // Initialize grid
    double currentTime = 0.;
    Grid3D grid = initGrid();
    grid = chargeCurrentUpdate(grid,currentTime);
    grid = electricFieldRelaxation(grid);
    cout << "Gauss' Law initial estimate: q = " << to_string(GaussLaw(grid)) << endl;

    writeOut(grid,0);
    cout << "Initial conditions written to file. Beginning loop..." << endl;

    for (int count=1;count <= LOOPCOUNT;count++) {       // count is just the loop counter variable
        // Magnetic field half-step
        currentTime += 0.5*TIMESTEP;
        grid = chargeCurrentUpdate(grid,currentTime);
        grid = magneticFieldUpdate(grid);

        // Electric field half-step
        currentTime += 0.5*TIMESTEP;
        grid = chargeCurrentUpdate(grid,currentTime);
        grid = electricFieldUpdate(grid);
        // grid = applyABC(grid);

        writeOut(grid,count);
        cout << "Finished loop " << to_string(count) << " of " << to_string(LOOPCOUNT);
        cout << "; Gauss' Law estimate q = " << to_string(GaussLaw(grid)) << endl;
    }

    return 0;
    
}
