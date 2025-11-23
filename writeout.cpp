#include <fstream>
#include <iomanip>

#include "header.h"
// #include "structs.cpp"
using namespace std;

void writeOut(Grid3D grid,int id) {
    // Defining files that hold data
    ofstream exfile("outputs/ex/ex"+to_string(id)+".dat");
    ofstream eyfile("outputs/ey/ey"+to_string(id)+".dat");
    ofstream ezfile("outputs/ez/ez"+to_string(id)+".dat");
    ofstream hxfile("outputs/hx/hx"+to_string(id)+".dat");
    ofstream hyfile("outputs/hy/hy"+to_string(id)+".dat");
    ofstream hzfile("outputs/hz/hz"+to_string(id)+".dat");
    ofstream jxfile("outputs/jx/jx"+to_string(id)+".dat");
    ofstream jyfile("outputs/jy/jy"+to_string(id)+".dat");
    ofstream pfile("outputs/p/p"+to_string(id)+".dat");
    ofstream phifile("outputs/phi/phi"+to_string(id)+".dat");
    ofstream econfile("outputs/econ/econ"+to_string(id)+".dat");
    ofstream hconfile("outputs/hcon/hcon"+to_string(id)+".dat");

    // Write grid side length, for reshaping data later
    exfile << GRIDSIZE << endl;
    eyfile << GRIDSIZE << endl;
    ezfile << GRIDSIZE << endl;
    hxfile << GRIDSIZE << endl;
    hyfile << GRIDSIZE << endl;
    hzfile << GRIDSIZE << endl;
    jxfile << GRIDSIZE << endl;
    jyfile << GRIDSIZE << endl;
    pfile << GRIDSIZE << endl;
    phifile << GRIDSIZE << endl;
    econfile << GRIDSIZE << endl;
    hconfile << GRIDSIZE << endl;

    // Set precision t avoid saving unneccessarily precise data
    exfile << setprecision(4);
    eyfile << setprecision(4);
    ezfile << setprecision(4);
    hxfile << setprecision(4);
    hyfile << setprecision(4);
    hzfile << setprecision(4);
    jxfile << setprecision(4);
    jyfile << setprecision(4);
    pfile << setprecision(4);
    phifile << setprecision(4);
    econfile << setprecision(4);
    hconfile << setprecision(4);

    // Write data
    for (int i=0;i<GRIDSIZE;i++) {
        for (int j=0;j<GRIDSIZE;j++) {
            for (int k=0;k<GRIDSIZE;k++) {
                // Cell c = grid[i][j][k];
                exfile << grid[i][j][k].Ex << "\t";
                eyfile << grid[i][j][k].Ey << "\t";
                ezfile << grid[i][j][k].Ez << "\t";
                hxfile << grid[i][j][k].Hx << "\t";
                hyfile << grid[i][j][k].Hy << "\t";
                hzfile << grid[i][j][k].Hz << "\t";
                jxfile << grid[i][j][k].Jx << "\t";
                jyfile << grid[i][j][k].Jy << "\t";
                pfile << grid[i][j][k].p << "\t";
                phifile << grid[i][j][k].phi << "\t";
                econfile << grid[i][j][k].EConstraint << "\t";
                hconfile << grid[i][j][k].HConstraint << "\t";
            }
        }
    }

    // Close files
    exfile.close();
    eyfile.close();
    ezfile.close();
    hxfile.close();
    hyfile.close();
    hzfile.close();
    jxfile.close();
    jyfile.close();
    pfile.close();
    phifile.close();
    econfile.close();
    hconfile.close();
}
