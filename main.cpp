#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
using namespace std;

#include "header.h"
#include "structs.cpp"
#include "writeout.cpp"

//// Grid initialization function
Grid3D initGrid() {
    // Define empty cell and empty grid
    Cell emptyCell = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    Grid3D grid(GRIDSIZE,vector<vector<Cell>>(GRIDSIZE,vector<Cell>(GRIDSIZE,emptyCell)));

    // Set cell coordinates and checkerboard pattern
    for (int i = 0; i < GRIDSIZE; i++) {
        for (int j = 0; j < GRIDSIZE; j++) {
            for (int k = 0; k < GRIDSIZE; k++) {
                // Set position of each cell center, units of dipole radii
                grid[i][j][k].x = GRIDSTEP*(i-GRIDRADIUS);
                grid[i][j][k].y = GRIDSTEP*(j-GRIDRADIUS);
                grid[i][j][k].z = GRIDSTEP*(k-GRIDRADIUS);

                // Set color. Used in the electric field relaxation function later on.
                grid[i][j][k].color = (i+j+k)%2;
            }
        }
    }
    cout << "Grid initialized" << endl;
    return grid;
}

//// Update functions
vector<double> updateDipole(Cell c,double t) {
    // Initialize results vector
    vector<double> res(3);

    // Calculate locations of dipole charges (equation 4, left)
    double xpos = cos(SPEED*t);
    double ypos = sin(SPEED*t);
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

    // Calculate each charge as a thin Gaussian distribution (equation 5)
    double Gpos = pow(GRIDSTEP/(sqrt(PI)*SIGMA),3)*exp(-(dxpos*dxpos + dypos*dypos + dzpos*dzpos)/(SIGMA*SIGMA));
    double Gneg = pow(GRIDSTEP/(sqrt(PI)*SIGMA),3)*exp(-(dxneg*dxneg + dyneg*dyneg + dzneg*dzneg)/(SIGMA*SIGMA));

    // Update results vector
    res[0] = Gpos - Gneg;                       // Charge density (equation 6)
    res[1] = -SPEED*sin(SPEED*t)*(Gpos+Gneg);   // Current density x-component (equations 4, right, and 7)
    res[2] = SPEED*cos(SPEED*t)*(Gpos+Gneg);    // Current density y-component (equations 4, right, and 7)
    return res;
}

Grid3D chargeCurrentUpdate(Grid3D grid,double currentTime) {
    // This function calculates the charge and current density distribution within the grid at a
    //      given time.

    // Loop over all grid cells
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                // Get charge density, current density at this cell center
                vector<double> du = updateDipole(grid[i][j][k],currentTime);
                grid[i][j][k].p = du[0];
                grid[i][j][k].Jx = du[1];     
                grid[i][j][k].Jy = du[2];     
            }
        }
    }
    return grid;
}

Grid3D electricFieldRelaxation(Grid3D grid) {
    // This function carries out the electric field relaxation scheme laid out in the project
    //      document. Relevant equations are cited where used in the code.

    // Iteratively smooth the electric potential ("phi") over the grid. The variable "t" controls
    //      the number of loops. The cells are updated in an alternating pattern using the Gauss-
    //      Seidel method (equations 20 and 22)
    for (int t = 0;t<100;t++) {
        // First, calculate for all the "red" squares (color=1)
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    // Save previous sweep's phi value
                    grid[i][j][k].oldPhi = grid[i][j][k].phi;
                    
                    // Skip "black" squares (color=0)
                    if (grid[i][j][k].color==0) {continue;}

                    // Calculate initial estimate for cell's new phi value (equation 22)
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

        // Now do the "black" squares (color=0)
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    // Save previous sweep's phi value
                    grid[i][j][k].oldPhi = grid[i][j][k].phi;
                
                    // Skip "red" squares (color=1)
                    if (grid[i][j][k].color==1) {continue;}

                    // Calculate initial estimate for cell's new phi value (equation 22)
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

        // Use successive over-relaxation to accelerate convergence (equation 23). We use a 
        //      relaxation parameter of ω = 1.85.
        double omega = 1.85;
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    grid[i][j][k].phi = (1.-omega)*grid[i][j][k].oldPhi + omega*grid[i][j][k].phi;
                }
            }
        }

        // Calculate residuals. In theory, the laplacian of the electric potential, which is the
        //      same as the divergence of the electric field, should be proportional to the charge
        //      density of a cell. We calculate the root mean square of residuals to measure how
        //      accurate the electric potential relaxation scheme has been (equation 24).
        double residual;
        double sum_residuals;
        double sum_charges;
        for (int i=0;i<GRIDSIZE;i++) {
            for (int j=0;j<GRIDSIZE;j++) {
                for (int k=0;k<GRIDSIZE;k++) {
                    // Calculate the laplacian of the electric potential within this cell using
                    //      equation 26.
                    double laplacian = 0;
                    if (i != 0) {laplacian += grid[i-1][j][k].phi;}
                    if (j != 0) {laplacian += grid[i][j-1][k].phi;}
                    if (k != 0) {laplacian += grid[i][j][k-1].phi;}
                    if (i != IMAX) {laplacian += grid[IMAX][j][k].phi;}
                    if (j != IMAX) {laplacian += grid[i][IMAX][k].phi;}
                    if (k != IMAX) {laplacian += grid[i][j][IMAX].phi;}
                    laplacian = (laplacian - 6.*grid[i][j][k].phi)/(GRIDSTEP*GRIDSTEP);

                    double residual = -pow(GRIDSTEP,3.)*grid[i][j][k].p - laplacian;
                    sum_residuals += residual;
                    sum_charges += pow(grid[i][j][k].p,2.);
                }
            }
        }
        // Note: the residual factor is not required to meet a stopping criterion in the current
        //      version of the code.
        cout << "Residual factor: " << sqrt(sum_residuals/sum_charges) << endl;
    }

    // Finally, calculate electric field components for each cell. Electric field components are
    //      calculated using the values of phi in each cell and its neighbours using equations 27.
    //      We use a Dirichlet boundary condition on boundary cells where phi is assumed to be 0
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
    // This function updates the E-field components of each cell in the grid. Updates are made
    //      using equations 13-15. This function also checks that Maxwell's equations are being
    //      satisfied by calculating the error in Gauss' law within each cell

    for (int i=1;i<GRIDSIZE-1;i++) {
        for (int j=1;j<GRIDSIZE-1;j++) {
            for (int k=1;k<GRIDSIZE-1;k++) {
                // Save previous e-field components
                grid[i][j][k].Ex0 = grid[i][j][k].Ex;
                grid[i][j][k].Ey0 = grid[i][j][k].Ey;
                grid[i][j][k].Ez0 = grid[i][j][k].Ez;

                // New e-field components
                double newEx,newEy,newEz;

                // Calculate Ex (equation 13)
                newEx = (grid[i][j][k].Hz - grid[i][j][k].Hy);
                if (j>0) {newEx -= grid[i][j-1][k].Hz;}
                if (k>0) {newEx += grid[i][j][k-1].Hy;}
                grid[i][j][k].Ex += 0.5*(newEx - GRIDSTEP*grid[i][j][k].Jx);

                // Calculate Ey (equation 14)
                newEy = (grid[i][j][k].Hx - grid[i][j][k].Hz);
                if (k>0) {newEy -= grid[i][j][k-1].Hx;}
                if (i>0) {newEy += grid[i-1][j][k].Hz;}
                grid[i][j][k].Ey += 0.5*(newEy - GRIDSTEP*grid[i][j][k].Jy);

                // Calculate Ez (equation 15)
                newEz = (grid[i][j][k].Hy - grid[i][j][k].Hx);
                if (i>0) {newEz -= grid[i-1][j][k].Hy;}
                if (j>0) {newEz += grid[i][j-1][k].Hx;}
                grid[i][j][k].Ez += 0.5*newEz;

                // Calculate E-field constraint. By Gauss' law (equations 3, left, and 17), it 
                //      should be true that the divergence of the E-field in a cell is equal to the
                //      charge density over ε0. The constraint value is charge density minus charge 
                //      density over ε0. The divergence of the E-field within a cell is calculated
                //      using equation 16.
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
    // This function updates the H-field components of each cell within a grid. Updates are made
    //      using equations 10-12. This function also checks that Maxwell's equations are being
    //      satisfied by calculating the error in Gauss' law for the magnetic field within each
    //      cell: that the divergence of the H-field is zero -- the "no magnetic monopole" constraint.

    for (int i=1;i<GRIDSIZE-1;i++) {
        for (int j=1;j<GRIDSIZE-1;j++) {
            for (int k=1;k<GRIDSIZE-1;k++) {
                // Save previous h-field components
                grid[i][j][k].Hx0 = grid[i][j][k].Hx;
                grid[i][j][k].Hy0 = grid[i][j][k].Hy;
                grid[i][j][k].Hz0 = grid[i][j][k].Hz;

                double newHx,newHy,newHz;
                // Calculate Hx (equation 10)
                newHx = (-grid[i][j][k].Ez + grid[i][j][k].Ey);
                if (j<GRIDSIZE-1) {newHx += grid[i][j+1][k].Ez;}
                if (k<GRIDSIZE-1) {newHx -= grid[i][j][k+1].Ey;}
                grid[i][j][k].Hx -= 0.5*newHx;
                // Calculate Hy (equation 11)
                newHy = (-grid[i][j][k].Ex+grid[i][j][k].Ez);
                if (k<GRIDSIZE-1) {newHy += grid[i][j][k+1].Ex;}
                if (i<GRIDSIZE-1) {newHy -= grid[i+1][j][k].Ez;}
                grid[i][j][k].Hy -= 0.5*newHy;
                // Calculate Hz (equation 12)
                newHz = (-grid[i][j][k].Ey + grid[i][j][k].Ex);
                if (i<GRIDSIZE-1) {newHz += grid[i+1][j][k].Ey;}
                if (j<GRIDSIZE-1) {newHz -= grid[i][j+1][k].Ex;}
                grid[i][j][k].Hz -= 0.5*newHz;

                // Calculate H-field constraint. Magnetic mononpoles are not believed to exist,
                //      therefore the divergence of the magnetic field must be 0 everywhere 
                //      (equation 3, right). The divergence of the H-field is calculated using a 
                //      modified version of equation 16.
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
    // This function applies Mur Absorbing Boundary Conditions (ABC's) on the boundary cells of
    //      the grid. An accurate implementation extrapolates how the E- and H-fields on the grid
    //      boundaries would behave if the grid extended further by solving the discrete wave 
    //      equation at each cell. The equations used are equations 28-29.

    double coeff = -0.5;        // This value is related to our fixed definition of Δt/Δx = 1/2,
                                //      which satisfies the stability condition in equation 31.

    // In this implementation, we first store E- and H-field updates in a temporary grid before
    //      updating the values in the main grid.
    Cell emptyCell = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    Grid3D newgrid(GRIDSIZE,vector<vector<Cell>>(GRIDSIZE,vector<Cell>(GRIDSIZE,emptyCell)));\

    // update x-faces
    for (int j=0; j < GRIDSIZE; j++) {
        for (int k=0; k < GRIDSIZE; k++) {
            // upper boundary
            newgrid[IMAX][j][k].Ey += grid[IMAX-1][j][k].Ey0 + coeff*(grid[IMAX-1][j][k].Ey0 - grid[IMAX][j][k].Ey0);
            newgrid[IMAX][j][k].Ez += grid[IMAX-1][j][k].Ez0 + coeff*(grid[IMAX-1][j][k].Ez0 - grid[IMAX][j][k].Ez0);
            newgrid[IMAX][j][k].Hy += grid[IMAX-1][j][k].Hy0 + coeff*(grid[IMAX-1][j][k].Hy0 - grid[IMAX][j][k].Hy0);
            newgrid[IMAX][j][k].Hz += grid[IMAX-1][j][k].Hz0 + coeff*(grid[IMAX-1][j][k].Hz0 - grid[IMAX][j][k].Hz0);
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
            // upper boundary
            newgrid[i][IMAX][k].Ex += grid[i][IMAX-1][k].Ex0 + coeff*(grid[i][IMAX-1][k].Ex0 - grid[i][IMAX][k].Ex0);
            newgrid[i][IMAX][k].Ez += grid[i][IMAX-1][k].Ez0 + coeff*(grid[i][IMAX-1][k].Ez0 - grid[i][IMAX][k].Ez0);
            newgrid[i][IMAX][k].Hx += grid[i][IMAX-1][k].Hx0 + coeff*(grid[i][IMAX-1][k].Hx0 - grid[i][IMAX][k].Hx0);
            newgrid[i][IMAX][k].Hz += grid[i][IMAX-1][k].Hz0 + coeff*(grid[i][IMAX-1][k].Hz0 - grid[i][IMAX][k].Hz0);
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
            // upper boundary
            newgrid[i][j][IMAX].Ex += grid[i][j][IMAX-1].Ex0 + coeff*(grid[i][j][IMAX-1].Ex0 - grid[i][j][IMAX].Ex0);
            newgrid[i][j][IMAX].Ey += grid[i][j][IMAX-1].Ey0 + coeff*(grid[i][j][IMAX-1].Ey0 - grid[i][j][IMAX].Ey0);
            newgrid[i][j][IMAX].Hx += grid[i][j][IMAX-1].Hx0 + coeff*(grid[i][j][IMAX-1].Hx0 - grid[i][j][IMAX].Hx0);
            newgrid[i][j][IMAX].Hy += grid[i][j][IMAX-1].Hy0 + coeff*(grid[i][j][IMAX-1].Hy0 - grid[i][j][IMAX].Hy0);
            // lower boundary
            newgrid[i][j][0].Ex += grid[i][j][1].Ex0 + coeff*(grid[i][j][1].Ex0 - grid[i][j][0].Ex0);
            newgrid[i][j][0].Ey += grid[i][j][1].Ey0 + coeff*(grid[i][j][1].Ey0 - grid[i][j][0].Ey0);
            newgrid[i][j][0].Hx += grid[i][j][1].Hx0 + coeff*(grid[i][j][1].Hx0 - grid[i][j][0].Hx0);
            newgrid[i][j][0].Hy += grid[i][j][1].Hy0 + coeff*(grid[i][j][1].Hy0 - grid[i][j][0].Hy0);
            
        }
    }

    // Make updates to the main grid
    for (int ind1=0; ind1 < GRIDSIZE; ind1++) {
        for (int ind2=0; ind2 < GRIDSIZE; ind2++) {
            double edgefactor = 1.;

            // If a cell lies on a grid edge, its tangential field components will have been
            //      updated twice. By setting "edgefactor" to 0.5, we effectively take the average
            //      of the two updates to be the new tangential field compoent value.
            if ((ind1==0)||(ind2==0)||(ind1==IMAX)||(ind2==IMAX)) {edgefactor = 0.5;}

            // x-faces
            // upper boundary
            grid[IMAX][ind1][ind2].Ey = edgefactor * (newgrid[IMAX][ind1][ind2].Ey);
            grid[IMAX][ind1][ind2].Ez = edgefactor * (newgrid[IMAX][ind1][ind2].Ez);
            grid[IMAX][ind1][ind2].Hy = edgefactor * (newgrid[IMAX][ind1][ind2].Hy);
            grid[IMAX][ind1][ind2].Hz = edgefactor * (newgrid[IMAX][ind1][ind2].Hz);
            // lower boundary
            grid[0][ind1][ind2].Ey = edgefactor * (newgrid[0][ind1][ind2].Ey);
            grid[0][ind1][ind2].Ez = edgefactor * (newgrid[0][ind1][ind2].Ez);
            grid[0][ind1][ind2].Hy = edgefactor * (newgrid[0][ind1][ind2].Hy);
            grid[0][ind1][ind2].Hz = edgefactor * (newgrid[0][ind1][ind2].Hz);

            // y-faces
            // upper boundary
            grid[ind1][IMAX][ind2].Ex = edgefactor * (newgrid[ind1][IMAX][ind2].Ex);
            grid[ind1][IMAX][ind2].Ez = edgefactor * (newgrid[ind1][IMAX][ind2].Ez);
            grid[ind1][IMAX][ind2].Hx = edgefactor * (newgrid[ind1][IMAX][ind2].Hx);
            grid[ind1][IMAX][ind2].Hz = edgefactor * (newgrid[ind1][IMAX][ind2].Hz);
            // lower boundary
            grid[ind1][0][ind2].Ex = edgefactor * (newgrid[ind1][0][ind2].Ex);
            grid[ind1][0][ind2].Ez = edgefactor * (newgrid[ind1][0][ind2].Ez);
            grid[ind1][0][ind2].Hx = edgefactor * (newgrid[ind1][0][ind2].Hx);
            grid[ind1][0][ind2].Hz = edgefactor * (newgrid[ind1][0][ind2].Hz);

            // z-faces
            // upper boundary
            grid[ind1][ind2][IMAX].Ex = edgefactor * (newgrid[ind1][ind2][IMAX].Ex);
            grid[ind1][ind2][IMAX].Ey = edgefactor * (newgrid[ind1][ind2][IMAX].Ey);
            grid[ind1][ind2][IMAX].Hx = edgefactor * (newgrid[ind1][ind2][IMAX].Hx);
            grid[ind1][ind2][IMAX].Hy = edgefactor * (newgrid[ind1][ind2][IMAX].Hy);
            // lower boundary
            grid[ind1][ind2][0].Ex = edgefactor * (newgrid[ind1][ind2][0].Ex);
            grid[ind1][ind2][0].Ey = edgefactor * (newgrid[ind1][ind2][0].Ey);
            grid[ind1][ind2][0].Hx = edgefactor * (newgrid[ind1][ind2][0].Hx);
            grid[ind1][ind2][0].Hy = edgefactor * (newgrid[ind1][ind2][0].Hy);
        }
    }
    return grid;
}

//// Diagnostic functions
double GaussLaw(Grid3D grid) {
    // This function approximates Gauss' Law for the dipole system by summing outgoing E-field
    //      components over the boundary cells. Because the total enclosed charge of the system is
    //      zero, the result of this function should also be zero.

    double totalFlux;
    for (int ind1=0;ind1<GRIDSIZE;ind1++) {
        for (int ind2=0;ind2<GRIDSIZE;ind2++) {
            totalFlux += grid[IMAX][ind1][ind2].Ex - grid[0][ind1][ind2].Ex;
            totalFlux += grid[ind1][IMAX][ind2].Ey - grid[ind1][0][ind2].Ey;
            totalFlux += grid[ind1][ind2][IMAX].Ez - grid[ind1][ind2][0].Ez;
        }
    }
    return totalFlux;
}

double GridEnergy(Grid3D grid) {
    // This function calculates the total energy stored in the E- and H-fields within the grid.
    //      The energy of such a system is given by equation in section XI, question (1).
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
    // This function calculates the Poynting vector of a specific cell within a given grid. The
    //      Poynting vector descriebs the direction in which energy is transferred by an electro-
    //      magnetic wave as it propogates. Sufficiently far from the dipole, the Poynting vector
    //      should oscillate between pointing towards and away from the dipole center.
    double xcomp = (grid[i][j][k].Ey * grid[i][j][k].Hz) - (grid[i][j][k].Ez * grid[i][j][k].Hy);
    double ycomp = (grid[i][j][k].Ez * grid[i][j][k].Hx) - (grid[i][j][k].Ex * grid[i][j][k].Hz);
    double zcomp = (grid[i][j][k].Ex * grid[i][j][k].Hy) - (grid[i][j][k].Ey * grid[i][j][k].Hx);

    vector<double> res = {xcomp,ycomp,zcomp};
    return res;
}

// Main
int main() {
    // Initialize grid and time
    double currentTime = 0.;
    Grid3D grid = initGrid();

    // Make initial updates to charge, current density, and E-field within the grid.
    grid = chargeCurrentUpdate(grid,currentTime);
    grid = electricFieldRelaxation(grid);

    // Output initial estimate for Gauss' law
    cout << setprecision(10);
    cout << "Gauss' Law initial estimate: q = " << GaussLaw(grid) << endl;
    
    // Save grid intial conditions
    writeOut(grid,0);
    cout << "Initial conditions written to file. Beginning loop..." << endl;

    // Open summary file (diagnostics)
    ofstream summary("outputs/summary.dat");
    summary << "Step \tGauss' Law \t\tEnergy\t\t\tPoynting Vector (Scaled)" << endl;

    // Execute main loop
    for (int count=1;count <= LOOPCOUNT;count++) { 
        // Magnetic field half-step
        currentTime += 0.5*TIMESTEP;
        grid = chargeCurrentUpdate(grid,currentTime);
        grid = magneticFieldUpdate(grid);

        // Electric field half-step
        currentTime += 0.5*TIMESTEP;
        grid = electricFieldUpdate(grid);
        grid = applyABC(grid);

        // Save end-of-step grid conditions
        writeOut(grid,count);

        // Output diagnostics to terminal
        double g = GaussLaw(grid);
        double e = GridEnergy(grid);
        cout << "Finished loop " << to_string(count) << " of " << to_string(LOOPCOUNT);
        cout << "; Gauss' Law estimate q = " << g;
        cout << "; Total energy E = " << e << endl;
        summary << count << "\t\t" << g << "\t\t" << e << "\t\t\t";
    }
    return 0;
}
// QED