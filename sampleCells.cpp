/*
* sampleCells.cpp - sample cells for DSMC simulation
*
* Takes Nx3 velocity of molecules (inV)
* and a partition of molecules into cells (partitions)
* and the Vcell, Tcell to fill in 
* and outputs the calculated Vcell, Tcell
*
* Usage : from MATLAB
*         >> [Vcell, Tcell] = sampleCells(inV, partitions, Vcell, Tcell)
*
*
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <cstdlib>
#include <cmath> 
#include <ctime> 

class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        matlab::data::TypedArray<double> V = std::move(inputs[0]);
        matlab::data::TypedArray<double> Partitions = std::move(inputs[1]);
        matlab::data::TypedArray<double> VCell = std::move(inputs[2]);
        matlab::data::TypedArray<double> TCell = std::move(inputs[3]);
        sampleCells(V, Partitions, VCell, TCell);
        outputs[0] = std::move(VCell);
        outputs[1] = std::move(TCell);
    }

    void sampleCells(matlab::data::TypedArray<double>& inV, matlab::data::TypedArray<double>& partitions, matlab::data::TypedArray<double>& Vcell, matlab::data::TypedArray<double>& Tcell) {
        double kB = 1.3806504e-23;
        double m = 6.6465e-27;
        for (int i = 0; i < partitions.getDimensions()[0] - 1; i++) {
            int Nstart = (int)(partitions[i][0] + 0.5);
            int Nend = (int)(partitions[(i + 1)][0] + 0.5);
            int Np = Nend - Nstart;
            double Vx = 0.0;
            double Vy = 0.0;
            double Vz = 0.0;
            double Vx2 = 0.0;
            double Vy2 = 0.0;
            double Vz2 = 0.0;
            for (int j = 0; j < Np; j++) {
                double tempx = inV[(Nstart + j)][0];
                double tempy = inV[(Nstart + j)][1];
                double tempz = inV[(Nstart + j)][2];
                Vx = Vx + tempx;
                Vy = Vy + tempy;
                Vz = Vz + tempz;
                Vx2 = Vx2 + tempx * tempx;
                Vy2 = Vy2 + tempy * tempy;
                Vz2 = Vz2 + tempz * tempz;
            }
            Vcell[i][0] = Vx / Np;
            Vcell[i][1] = Vy / Np;
            Vcell[i][2] = Vz / Np;
            Tcell[i][0] = (Vx2 / Np - Vx * Vx / Np / Np) * m / kB;
            Tcell[i][1] = (Vy2 / Np - Vy * Vy / Np / Np) * m / kB;
            Tcell[i][2] = (Vz2 / Np - Vz * Vz / Np / Np) * m / kB;
        }
    }

};