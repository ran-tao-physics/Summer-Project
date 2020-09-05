/*
* moveMolecules.cpp - move molecules for DSMC simulation
*
* Takes the skimmer surface temperature (T)
* and the time step (dt)
* and positions of the molecules (inX)
* and velocities of the molecules (inV)
* and the output of stlread for the geometry (F,V,N) the normals have to be pointing into the flow region
* and outputs the new positions of the molecules (output)
*
* Usage : from MATLAB
*         >> [X, V] = moveMolecules(T, dt, inX, inV, F, V, N)
*
*
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <cstdlib>
#include <cmath> 
#include <ctime> 

using matlab::mex::ArgumentList;
using namespace matlab::data;


class MexFunction : public matlab::mex::Function {
    // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    ArrayFactory factory;

    // Create an output stream
    std::ostringstream stream;
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        std::srand(std::time(nullptr));
        double Tskim = inputs[0][0];
        double dT = inputs[1][0];
        matlab::data::TypedArray<double> X = std::move(inputs[2]);
        matlab::data::TypedArray<double> V = std::move(inputs[3]);
        matlab::data::TypedArray<double> Fskim = std::move(inputs[4]);
        matlab::data::TypedArray<double> Vskim = std::move(inputs[5]);
        matlab::data::TypedArray<double> Nskim = std::move(inputs[6]);
        moveMolecules(Tskim, dT, X, V, Fskim, Vskim, Nskim);
        outputs[0] = std::move(X);
        outputs[1] = std::move(V);
    }

    void moveMolecules(double T, double dt, matlab::data::TypedArray<double>& inX, matlab::data::TypedArray<double>& inV, matlab::data::TypedArray<double>& F, matlab::data::TypedArray<double>& V, matlab::data::TypedArray<double>& N) {
        int numMol = inX.getDimensions()[0];
        int numSurf = F.getDimensions()[0];
        for (int i = 0; i < numMol; i++) {
            double tremain = dt;
            double oldx = inX[i][0];
            double oldy = inX[i][1];
            double oldz = inX[i][2];
            double oldvx = inV[i][0];
            double oldvy = inV[i][1];
            double oldvz = inV[i][2];
            for (int j = 0; j < numSurf; j++) {
                if (tremain == 0) {
                    break;
                }
                else {
                    /*accomodate for 0-based arrays*/
                    int ind1 = (int)(F[j][0] + 0.5) - 1;
                    int ind2 = (int)(F[j][1] + 0.5) - 1;
                    int ind3 = (int)(F[j][2] + 0.5) - 1;
                    double ax = V[ind1][0];
                    double ay = V[ind1][1];
                    double az = V[ind1][2];
                    double bx = V[ind2][0];
                    double by = V[ind2][1];
                    double bz = V[ind2][2];
                    double cx = V[ind3][0];
                    double cy = V[ind3][1];
                    double cz = V[ind3][2];
                    double nx = N[j][0];
                    double ny = N[j][1];
                    double nz = N[j][2];
                    double direction_check = oldvx * nx + oldvy * ny + oldvz * nz;
                    double t = ((ax - oldx) * nx + (ay - oldy) * ny + (az - oldz) * nz) / direction_check;
                    if (direction_check < 0 && t <= tremain && t > 0) { /*possibility of needing to take the minimum of t for all surface patches, but ignored here, assuming geometry is convex*/
                        double contactx = oldx + oldvx * t;
                        double contacty = oldy + oldvy * t;
                        double contactz = oldz + oldvz * t;
                        if (checkInside(ax, ay, az, bx, by, bz, cx, cy, cz, contactx, contacty, contactz)) {
                            oldx = contactx;
                            oldy = contacty;
                            oldz = contactz;
                            double* newV = generateV(T, nx, ny, nz, ax - bx, ay - by, az - bz);
                            oldvx = newV[0];
                            oldvy = newV[1];
                            oldvz = newV[2];
                            tremain = tremain - t;
                            j = -1; /*restart the cycle, check for secondary collisions*/
                            continue;
                        }
                    }
                }
            }
            inX[i][0] = oldx + oldvx * tremain;
            inX[i][1] = oldy + oldvy * tremain;
            inX[i][2] = oldz + oldvz * tremain;
            inV[i][0] = oldvx;
            inV[i][1] = oldvy;
            inV[i][2] = oldvz;
        }
    }

    bool checkInside(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double p3x, double p3y, double p3z, double x, double y, double z) {
        /*use barycentric coordinates to check if point is inside a triangle, test does not depend on counter-clockwise or clockwise numbering*/
        double area = 0.5 * (p1z * (p2x * p3y - p2y * p3x) + p1y * (-p2x * p3z + p3x * p2z) + p1x * (p2y * p3z - p3y * p2z));
        double alpha = 1 / (2 * area) * (z * (p2x * p3y - p2y * p3x) + y * (-p2x * p3z + p3x * p2z) + x * (p2y * p3z - p3y * p2z));
        double beta = 1 / (2 * area) * (z * (p3x * p1y - p3y * p1x) + y * (-p3x * p1z + p1x * p3z) + x * (p3y * p1z - p1y * p3z));
        double gamma = 1 / (2 * area) * (z * (p1x * p2y - p1y * p2x) + y * (-p1x * p2z + p2x * p1z) + x * (p1y * p2z - p2y * p1z));
        bool out = (alpha >= 0) && (beta >= 0) && (gamma >= 0);
        return out;
    }
    
    double norm(double a1, double a2, double a3) {
        return std::sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    }

    double* generateV(double T, double Nx, double Ny, double Nz, double Sx, double Sy, double Sz) {
        double kB = 1.3806504e-23;
        double m = 6.6465e-27;
        double pi = 2 * std::acos(0.0);
        double Cm = std::sqrt(2 * kB * T / m);
        /*construct the face normal coordinate system*/
        /*first the normal vector*/
        double norm1 = norm(Nx, Ny, Nz);
        Nx = Nx / norm1;
        Ny = Ny / norm1;
        Nz = Nz / norm1;
        /*then one of the side vectors*/
        double norm2 = norm(Sx, Sy, Sz);
        Sx = Sx / norm2;
        Sy = Sy / norm2;
        Sz = Sz / norm2;
        /*the remaining vector from cross product*/
        double rx = Ny * Sz - Nz * Sy;
        double ry = Nz * Sx - Nx * Sz;
        double rz = Nx * Sy - Ny * Sx;
        double Cx_face = Cm * std::sqrt(-std::log((double)rand() / (RAND_MAX)));
        double Cy_face = Cm * std::sin(2 * pi * (double)rand() / (RAND_MAX)) * std::sqrt(-std::log((double)rand() / (RAND_MAX)));
        double Cz_face = Cm * std::sin(2 * pi * (double)rand() / (RAND_MAX)) * std::sqrt(-std::log((double)rand() / (RAND_MAX)));
        double* outC = new double[3];
        outC[0] = Nx * Cx_face + Sx * Cy_face + rx * Cz_face;
        outC[1] = Ny * Cx_face + Sy * Cy_face + ry * Cz_face;
        outC[2] = Nz * Cx_face + Sz * Cy_face + rz * Cz_face;
        return outC;
    }

    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }
};