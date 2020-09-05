/*
* collision.cpp - calculate molecule collisions for DSMC simulation
*
* Takes a reference constant for VHS model (sigma_ref = pi*dref^2*gref^(2*omega-1))
* and the particle weight (Wp)
* and the time step (time_step)
* and the volume of a collision cell (cell_volume)
* and particle position and velocity (inX, inV)
* and the partition of particles into cells (partitions)
* and the (sigma*g)max of each cell (sgmax)
* and the cell mean collision time (mct)
* and the cell mean collision separation (mcs)
* and outputs the after-collision velocities of each molecule (outV)
* and the new mct, mcs matrices
*
* Usage : from MATLAB
*         >> [outV,mct,mcs,sgmax] = collision(sigma_ref,Wp,time_step,cell_volume,inX,inV,partitions,sgmax,mct,mcs)
*
*
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <cstdlib>
#include <cmath> 
#include <ctime> 
#include <algorithm>
#include <random>

using matlab::mex::ArgumentList;
using namespace matlab::data;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        double Sigma_ref = inputs[0][0];
        double WP = inputs[1][0];
        double dT = inputs[2][0];
        double VCell = inputs[3][0];
        matlab::data::TypedArray<double> X = std::move(inputs[4]);
        matlab::data::TypedArray<double> V = std::move(inputs[5]);
        matlab::data::TypedArray<double> Partitions = std::move(inputs[6]);
        matlab::data::TypedArray<double> SGMAX = std::move(inputs[7]);
        matlab::data::TypedArray<double> MCT = std::move(inputs[8]);
        matlab::data::TypedArray<double> MCS = std::move(inputs[9]);
        collision(Sigma_ref, WP, dT, VCell, X, V, Partitions, SGMAX, MCT, MCS);
        outputs[0] = std::move(V);
        outputs[1] = std::move(MCT);
        outputs[2] = std::move(MCS);
        outputs[3] = std::move(SGMAX);
    }

    void collision(double sigma_ref, double Wp, double dt, double Vcell, matlab::data::TypedArray<double>& inX, matlab::data::TypedArray<double>& inV, matlab::data::TypedArray<double>& partitions, matlab::data::TypedArray<double>& sgmax, matlab::data::TypedArray<double>& mct, matlab::data::TypedArray<double>& mcs) {
        double pi = 2 * std::acos(0.0);
        double omega = 0.66; /*power law exponent for He*/
        for (int i = 0; i < partitions.getDimensions()[0] - 1; i++) { /*for each cell*/
            std::srand(std::time(nullptr)); /*seed srand for velocity distribution*/
            std::mt19937 rng(std::time(NULL)); /*seed for integer distribution*/
            /*number of molecules in the current cell declared as double because matlab stored it as double*/
            int Nstart = (int)(partitions[i][0] + 0.5);
            int Nend = (int)(partitions[i + 1][0] + 0.5);
            int Np = Nend - Nstart;
            if (Np != 1) {
                double sgmax_old = sgmax[i][0];
                double nCol = 0.5 * Np * (Np - 1) * sgmax_old * Wp * dt / Vcell;
                int nTest = (int)(std::min(std::floor(nCol+0.5),std::floor(Np/2.0))+0.5); 
                double F = nCol / nTest;
                double sgmax_new = 0;
                double mcs_new = 0;
                int nAccepted = 0;
                for (int j = 0; j < nTest; j++) { /*for each molecule pair tested*/
                    /*generate random number engine for the specified range of indices*/
                    std::uniform_int_distribution<> distrib(Nstart, Nend-1);
                    int ind_a = distrib(rng);/*choose index of the first molecule*/
                    int ind_b = distrib(rng);
                    while (ind_b == ind_a) { /*prevent the two indices being the same*/
                        ind_b = distrib(rng);
                    }
                    double g = vlength(inV[ind_a][0] - inV[ind_b][0], inV[ind_a][1] - inV[ind_b][1], inV[ind_a][2] - inV[ind_b][2]); /*relative speed*/
                    double sg = sigma_ref / std::pow(g, 2 - 2 * omega); /*cross section area times relative speed*/
                    double P = F * sg / sgmax_old; /*probability of acceptance*/
                    double Rtest = (double)rand() / (RAND_MAX);
                    if (Rtest < P) { /*accepted*/
                        nAccepted += 1;
                        /*COM velocities*/
                        double w1 = (inV[ind_a][0] + inV[ind_b][0]) / 2;
                        double w2 = (inV[ind_a][1] + inV[ind_b][1]) / 2;
                        double w3 = (inV[ind_a][2] + inV[ind_b][2]) / 2;
                        /*calculate relative velocities for VHS model*/
                        double R1 = (double)rand() / (RAND_MAX);
                        double R2 = (double)rand() / (RAND_MAX);
                        double g1 = g * (2 * R1 - 1);
                        double g2 = g * std::sqrt(1 - (2 * R1 - 1) * (2 * R1 - 1)) * std::cos(2 * pi * R2);
                        double g3 = g * std::sqrt(1 - (2 * R1 - 1) * (2 * R1 - 1)) * std::sin(2 * pi * R2);
                        /*assignment*/
                        inV[ind_a][0] = w1 + g1 / 2;
                        inV[ind_a][1] = w2 + g2 / 2;
                        inV[ind_a][2] = w3 + g3 / 2;
                        inV[ind_b][0] = w1 - g1 / 2;
                        inV[ind_b][1] = w2 - g2 / 2;
                        inV[ind_b][2] = w3 - g3 / 2;
                        sgmax_new = std::max(sgmax_new, sg); /*update sgmax in this iteration*/
                        mcs_new += vlength(inX[ind_a][0] - inX[ind_b][0], inX[ind_a][1] - inX[ind_b][1], inX[ind_a][2] - inX[ind_b][2]); /*collision separation*/
                    }
                }
                sgmax[i][0] = std::max(sgmax_old, sgmax_new);
                mcs[i][0] = mcs_new / nAccepted;
                mct[i][0] = Np * dt / nAccepted;
            }
        }
    }

    double vlength(double a, double b, double c) {
        return std::sqrt(a * a + b * b + c * c);
    }

};