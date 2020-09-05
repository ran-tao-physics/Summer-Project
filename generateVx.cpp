/*
* generateEntry.cpp - generate entry molecules for DSMC simulation
*
* Takes a stagnation temperature (T)
* and a stream velocity in x direction (V)
* and an empty Nx1 matrix where N is the number of molecules generated (in)
* and outputs a Nx1 matrix of the generated velocities in x direction (output)
*
* Usage : from MATLAB
*         >> output = generateVx(T, V, in)
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
        std::srand(std::time(nullptr));
        double T = inputs[0][0]; 
        double V = inputs[1][0];
        matlab::data::TypedArray<double> in = std::move(inputs[2]);
        generateVx(T,V,in);
        outputs[0] = std::move(in);
    }

    void generateVx(double T, double V, matlab::data::TypedArray<double>& in) {
        int count = in.getDimensions()[0];
        double kB = 1.3806504e-23;
        double m = 6.6465e-27;
        double Cm = std::sqrt(2 * kB * T / m);
        double s = V / Cm;
        while (count > 0) {
                double y = 6 * ((double)rand() / (RAND_MAX)) - 3;
                double R2 = (double)rand() / (RAND_MAX);
                double h = std::sqrt(s * s + 2);
                double K = 2 * std::exp(0.5 + s * (s - h) / 2) / (s + h);
                if (R2 <= K*(y+s)*std::exp(-y*y)) {
                    count -= 1;
                    in[count][0] = V + Cm * y;
                }
        }
    }

};