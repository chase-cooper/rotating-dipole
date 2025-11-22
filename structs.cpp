#include <vector>
using namespace std;

// Define Cell struct that holds variables of each "cell" in the grid
struct Cell {
    double x,y,z;               // Position

    double p;                   // Charge density
    double Ex,Ey,Ez;            // Electric Field
    double Ex0,Ey0,Ez0;         // Previous step electric field
    double phi;                 // Electric potential
    double tempPhi;             // placeholder
    double oldPhi;              // Variable to store previous value of electric potential
    double laplacePhi;          // Laplacian of the electric potential
    double EConstraint;         // monitors electric field divergence constraint

    double Jx,Jy;               // Current density
    double Hx,Hy,Hz;            // Magnetic Field
    double Hx0,Hy0,Hz0;         // Previous step magnetic field
    double HConstraint;         // monitors "no magnetic monopoles" rule
    
    int color;                  // "color" used for electric field relaxation
};
typedef vector<vector<vector<Cell>>> Grid3D;

struct DipoleUpdate {
    double p;
    double Jx,Jy;
};