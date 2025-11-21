#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
using namespace std;

#include "header.h"
#include "structs.cpp"
#include "writeout.cpp"

// Dipole-related functions
struct DipoleUpdate {
    double p;
    double Jx,Jy;
};

DipoleUpdate updateDipole(Cell c,double t) {
    DipoleUpdate res;
    // Locations of dipole charges
    double xpos = cos(OMEGA*t);
    double ypos = sin(OMEGA*t);
    double zpos = 0.;
    double xneg = -xpos;
    double yneg = -ypos;
    double zneg = 0.;

    // Cell center displacement from positive charge
    double dxpos = c.x - xpos;
    double dypos = c.y - ypos;
    double dzpos = c.z - zpos;
    // Cell center displacement from negative charge
    double dxneg = c.x - xneg;
    double dyneg = c.y - yneg;
    double dzneg = c.z - zneg;

    // Calculate each charge as a thin Gaussian distribution
    double Gpos = pow(PI*SIGMA*SIGMA,1.5)*exp(-(dxpos*dxpos + dypos*dypos + dzpos*dzpos)/(SIGMA*SIGMA));
    double Gneg = pow(PI*SIGMA*SIGMA,1.5)*exp(-(dxneg*dxneg + dyneg*dyneg + dzneg*dzneg)/(SIGMA*SIGMA));

    res.p = Gpos - Gneg;                        // Charge density
    res.Jx = -OMEGA*sin(OMEGA*t)*(Gpos+Gneg);   // Current density x-component
    res.Jy = OMEGA*cos(OMEGA*t)*(Gpos+Gneg);    // Current density y-component
    return res;
}

// Grid update functions
Grid3D initGrid() {
    Cell emptyCell = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    Grid3D grid(GRIDSIZE,vector<vector<Cell>>(GRIDSIZE,vector<Cell>(GRIDSIZE,emptyCell)));

    // Set grid coordinates and checkerboard pattern
    double offset = (GRIDSIZE-1)/2.;        // Offset due to grid coordinates being [0,GRIDSIZE-1] but physical
                                            //  locations range [-L,L] for some L
    int index = 0;                          // Used to setup checkerboard alternation for electric field relaxation method
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                // Set position of each cell center
                grid[i][j][k].x = (i-offset)*GRIDSTEP;
                grid[i][j][k].y = (j-offset)*GRIDSTEP;
                grid[i][j][k].z = (k-offset)*GRIDSTEP;

                grid[i][j][k].color = index;
                index = (index+1)%2;
            }
            if (GRIDSIZE%2==0) {index = (index+1)%2;}
        }
        if (GRIDSIZE%2==0) {index = (index+1)%2;}
    }

    return grid;
}

Grid3D chargeUpdate(Grid3D grid,double currentTime) {
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                // Get charge density at this cell center
                DipoleUpdate du = updateDipole(grid[i][j][k],currentTime);
                grid[i][j][k].p = du.p;     // probably fine
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
                grid[i][j][k].Jx = du.Jx;     
                grid[i][j][k].Jy = du.Jy;     
            }
        }
    }
    return grid;
}

Grid3D electricFieldRelaxation(Grid3D grid) {
    for (int t = 0;t<100;t++) {
        // The below loops iteratively approximate the electric potential over the grid. The cells
        //  are split into red and black squares, and they are calculated separately.

        // First, calculate for all the "red" squares (color=1)
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    grid[i][j][k].oldPhi = grid[i][j][k].phi;
                
                    if (grid[i][j][k].color==0) {continue;}
                    double new_phi = -grid[i][j][k].p;
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
                    double new_phi = -grid[i][j][k].p;
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

        // Update the values of phi using relaxation. Additionally, we take this opportunity
        //  to calculate the laplacian of phi and get test for convergence.

        double omega = 1.8;     // Relaxation factor
        double epsilon = 1e+2;  // Convergence criterion
        double conv_total;      // Holds sum of laplacians
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    grid[i][j][k].phi = (1.-omega)*grid[i][j][k].oldPhi + omega*grid[i][j][k].phi;

                    double lp = (6./GRIDSTEP/GRIDSTEP)*(grid[i][j][k].oldPhi-grid[i][j][k].phi);
                    conv_total += pow((-grid[i][j][k].p - lp),2)*pow(GRIDSTEP,3);
                }
            }
        }
        conv_total = sqrt(conv_total);
        if (conv_total < epsilon) {continue;}

    }
    // Finally, calculate electric field components for each cell
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                double Ex = grid[i][j][k].phi;
                if (i != GRIDSIZE-1) {Ex -= grid[i+1][j][k].phi;}
                grid[i][j][k].Ex = Ex/GRIDSTEP;

                double Ey = grid[i][j][k].phi;
                if (j != GRIDSIZE-1) {Ey -= grid[i][j+1][k].phi;}
                grid[i][j][k].Ey = Ey/GRIDSTEP;

                double Ez = grid[i][j][k].phi;
                if (k != GRIDSIZE-1) {Ez -= grid[i][j][k+1].phi;}
                grid[i][j][k].Ez = Ez/GRIDSTEP;
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
                grid[i][j][k].Ex += TIMESTEP*(newEx/GRIDSTEP-grid[i][j][k].Jx)/EPS0;
                // Ey
                newEy = (grid[i][j][k].Hx - grid[i][j][k].Hz);
                if (k>0) {newEy -= grid[i][j][k-1].Hx;}
                if (i>0) {newEy += grid[i-1][j][k].Hz;}
                grid[i][j][k].Ey += TIMESTEP*(newEy/GRIDSTEP-grid[i][j][k].Jy)/EPS0;
                // Ez
                newEz = (grid[i][j][k].Hy - grid[i][j][k].Hx);
                if (i>0) {newEz -= grid[i-1][j][k].Hy;}
                if (j>0) {newEz += grid[i][j-1][k].Hx;}
                grid[i][j][k].Ez += TIMESTEP*newEz/EPS0/GRIDSTEP;
            }
        }
    }
    return grid;
}

Grid3D magneticFieldUpdate(Grid3D grid) {
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                double newHx,newHy,newHz;
                // Hx
                newHx = (-grid[i][j][k].Ez + grid[i][j][k].Ey);
                if (j<GRIDSIZE-1) {newHx += grid[i][j+1][k].Ez;}
                if (k<GRIDSIZE-1) {newHx -= grid[i][j][k+1].Ey;}
                grid[i][j][k].Hx -= TIMESTEP*newHx/MU0/GRIDSTEP;
                // Hy
                newHy = (-grid[i][j][k].Ex+grid[i][j][k].Ez);
                if (k<GRIDSIZE-1) {newHy += grid[i][j][k+1].Ex;}
                if (i<GRIDSIZE-1) {newHy -= grid[i+1][j][k].Ez;}
                grid[i][j][k].Hy -= TIMESTEP*newHy/MU0/GRIDSTEP;
                // Hz
                newHz = (-grid[i][j][k].Ey + grid[i][j][k].Ex);
                if (i<GRIDSIZE-1) {newHz += grid[i+1][j][k].Ey;}
                if (j<GRIDSIZE-1) {newHz -= grid[i][j+1][k].Ex;}
                grid[i][j][k].Hz -= TIMESTEP*newHz/MU0/GRIDSTEP;
            }
        }
    }
    return grid;
}

Grid3D applyABC(Grid3D grid) {

    // Need to be updated all at once, not one by one

    double coeff = (C*TIMESTEP - GRIDSTEP)/(C*TIMESTEP + GRIDSTEP);
    int imax=GRIDSIZE-1;        // max index (i/j/k) shorthand

    Cell emptyCell = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    Grid3D newgrid(GRIDSIZE,vector<vector<Cell>>(GRIDSIZE,vector<Cell>(GRIDSIZE,emptyCell)));

    for (int ind1=0;ind1<GRIDSIZE;ind1++) {
        for (int ind2=0;ind2<GRIDSIZE;ind2++) {
            // x-face upper boundary
            newgrid[imax][ind1][ind2].Ey = grid[imax-1][ind1][ind2].Ey0 + coeff*(grid[imax-1][ind1][ind2].Ey - grid[imax][ind1][ind2].Ey0);
            newgrid[imax][ind1][ind2].Ez = grid[imax-1][ind1][ind2].Ez0 + coeff*(grid[imax-1][ind1][ind2].Ez - grid[imax][ind1][ind2].Ez0);
            newgrid[imax][ind1][ind2].Hy = grid[imax-1][ind1][ind2].Hy0 + coeff*(grid[imax-1][ind1][ind2].Hy - grid[imax][ind1][ind2].Hy0);
            newgrid[imax][ind1][ind2].Hz = grid[imax-1][ind1][ind2].Hz0 + coeff*(grid[imax-1][ind1][ind2].Hz - grid[imax][ind1][ind2].Hz0);
            // x-face lower boundary
            newgrid[0][ind1][ind2].Ey = grid[1][ind1][ind2].Ey0 + coeff*(grid[1][ind1][ind2].Ey - grid[0][ind1][ind2].Ey0);
            newgrid[0][ind1][ind2].Ez = grid[1][ind1][ind2].Ez0 + coeff*(grid[1][ind1][ind2].Ez - grid[0][ind1][ind2].Ez0);
            newgrid[0][ind1][ind2].Hy = grid[1][ind1][ind2].Hy0 + coeff*(grid[1][ind1][ind2].Hy - grid[0][ind1][ind2].Hy0);
            newgrid[0][ind1][ind2].Hz = grid[1][ind1][ind2].Hz0 + coeff*(grid[1][ind1][ind2].Hz - grid[0][ind1][ind2].Hz0);
            // y-face upper boundary
            newgrid[ind1][imax][ind2].Ex = grid[ind1][imax-1][ind2].Ex0 + coeff*(grid[ind1][imax-1][ind2].Ex - grid[ind1][imax][ind2].Ex0);
            newgrid[ind1][imax][ind2].Ez = grid[ind1][imax-1][ind2].Ez0 + coeff*(grid[ind1][imax-1][ind2].Ez - grid[ind1][imax][ind2].Ez0);
            newgrid[ind1][imax][ind2].Hx = grid[ind1][imax-1][ind2].Hx0 + coeff*(grid[ind1][imax-1][ind2].Hx - grid[ind1][imax][ind2].Hx0);
            newgrid[ind1][imax][ind2].Hz = grid[ind1][imax-1][ind2].Hz0 + coeff*(grid[ind1][imax-1][ind2].Hz - grid[ind1][imax][ind2].Hz0);
            // y-face lower boundary
            newgrid[ind1][0][ind2].Ex = grid[ind1][1][ind2].Ex0 + coeff*(grid[ind1][1][ind2].Ex - grid[ind1][0][ind2].Ex0);
            newgrid[ind1][0][ind2].Ez = grid[ind1][1][ind2].Ez0 + coeff*(grid[ind1][1][ind2].Ez - grid[ind1][0][ind2].Ez0);
            newgrid[ind1][0][ind2].Hx = grid[ind1][1][ind2].Hx0 + coeff*(grid[ind1][1][ind2].Hx - grid[ind1][0][ind2].Hx0);
            newgrid[ind1][0][ind2].Hz = grid[ind1][1][ind2].Hz0 + coeff*(grid[ind1][1][ind2].Hz - grid[ind1][0][ind2].Hz0);
            // z-face upper boundary
            newgrid[ind1][ind2][imax].Ex = grid[ind1][ind2][imax-1].Ex0 + coeff*(grid[ind1][ind2][imax-1].Ex - grid[ind1][ind2][imax].Ex0);
            newgrid[ind1][ind2][imax].Ey = grid[ind1][ind2][imax-1].Ey0 + coeff*(grid[ind1][ind2][imax-1].Ey - grid[ind1][ind2][imax].Ey0);
            newgrid[ind1][ind2][imax].Hx = grid[ind1][ind2][imax-1].Hx0 + coeff*(grid[ind1][ind2][imax-1].Hx - grid[ind1][ind2][imax].Hx0);
            newgrid[ind1][ind2][imax].Hy = grid[ind1][ind2][imax-1].Hy0 + coeff*(grid[ind1][ind2][imax-1].Hy - grid[ind1][ind2][imax].Hy0);
            // z-face lower boundary
            newgrid[ind1][ind2][0].Ex = grid[ind1][ind2][1].Ex0 + coeff*(grid[ind1][ind2][1].Ex - grid[ind1][ind2][0].Ex0);
            newgrid[ind1][ind2][0].Ey = grid[ind1][ind2][1].Ey0 + coeff*(grid[ind1][ind2][1].Ey - grid[ind1][ind2][0].Ey0);
            newgrid[ind1][ind2][0].Hx = grid[ind1][ind2][1].Hx0 + coeff*(grid[ind1][ind2][1].Hx - grid[ind1][ind2][0].Hx0);
            newgrid[ind1][ind2][0].Hy = grid[ind1][ind2][1].Hy0 + coeff*(grid[ind1][ind2][1].Hy - grid[ind1][ind2][0].Hy0);
        }
    }

    for (int ind1=0;ind1<GRIDSIZE;ind1++) {
        for (int ind2=0;ind2<GRIDSIZE;ind2++) {
            // Update x-faces
            grid[imax][ind1][ind2].Ey = newgrid[imax][ind1][ind2].Ey;
            grid[imax][ind1][ind2].Ez = newgrid[imax][ind1][ind2].Ez;
            grid[imax][ind1][ind2].Hy = newgrid[imax][ind1][ind2].Hy;
            grid[imax][ind1][ind2].Hz = newgrid[imax][ind1][ind2].Hz;
            grid[0][ind1][ind2].Ey = newgrid[0][ind1][ind2].Ey;
            grid[0][ind1][ind2].Ez = newgrid[0][ind1][ind2].Ez;
            grid[0][ind1][ind2].Hy = newgrid[0][ind1][ind2].Hy;
            grid[0][ind1][ind2].Hz = newgrid[0][ind1][ind2].Hz;
            // Update y-faces
            grid[ind1][imax][ind2].Ex = newgrid[ind1][imax][ind2].Ex;
            grid[ind1][imax][ind2].Ez = newgrid[ind1][imax][ind2].Ez;
            grid[ind1][imax][ind2].Hx = newgrid[ind1][imax][ind2].Hx;
            grid[ind1][imax][ind2].Hz = newgrid[ind1][imax][ind2].Hz;
            grid[ind1][0][ind2].Ex = newgrid[ind1][0][ind2].Ex;
            grid[ind1][0][ind2].Ez = newgrid[ind1][0][ind2].Ez;
            grid[ind1][0][ind2].Hx = newgrid[ind1][0][ind2].Hx;
            grid[ind1][0][ind2].Hz = newgrid[ind1][0][ind2].Hz;
            // Update z-faces
            grid[ind1][ind2][imax].Ex = newgrid[ind1][ind2][imax].Ex;
            grid[ind1][ind2][imax].Ey = newgrid[ind1][ind2][imax].Ey;
            grid[ind1][ind2][imax].Hx = newgrid[ind1][ind2][imax].Hx;
            grid[ind1][ind2][imax].Hy = newgrid[ind1][ind2][imax].Hy;
            grid[ind1][ind2][0].Ex = newgrid[ind1][ind2][0].Ex;
            grid[ind1][ind2][0].Ey = newgrid[ind1][ind2][0].Ey;
            grid[ind1][ind2][0].Hx = newgrid[ind1][ind2][0].Hx;
            grid[ind1][ind2][0].Hy = newgrid[ind1][ind2][0].Hy;
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
    // Initialize grid
    double currentTime = 0.;
    Grid3D grid = initGrid();
    grid = chargeUpdate(grid,currentTime);

    // Electric Field
    grid = electricFieldRelaxation(grid);
    grid = writeOut(grid,0);
    
    // MAIN LOOP
    for (int i=1;i<500;i++) {
        currentTime += 0.5*TIMESTEP;

        grid = currentUpdate(grid,currentTime);
        grid = magneticFieldUpdate(grid);       // Check for issues
        
        currentTime += 0.5*TIMESTEP;

        grid = chargeUpdate(grid,currentTime);
        grid = electricFieldUpdate(grid);       // Check for issues
        // grid = applyABC(grid);

        // Writing all grid data
        grid = writeOut(grid,i);
        double energy = GaussLaw(grid);
        cout << to_string(energy) << endl;
    }

    

    
}
