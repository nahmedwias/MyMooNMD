/**
 * @file This example file contains an upward
 * directed (stationary) flow through a wound up pipe.
 * A large portion of the file deals with methods which are used to
 * transform a cylindrical sandwich geometry to a helically coiled tube
 * geometry. These methods have to be calles onto the TCollection
 * before the rest of the problem (FESpaces, FEFUnctions,
 * Matrices, Parallel Communicators etc.) is set up.
 *
 * @date 2016/04/26
 * @author Various, imported to ParMooN by Clemens Bartsch
 */
#ifndef TWISTED_PIPE_FLOW_
#define TWISTED_PIPE_FLOW_

void ExampleFile(bool time_dependent);

void ExactU1(double x, double y,  double z, double *values);

void ExactU2(double x, double y,  double z, double *values);

void ExactU3(double x, double y,  double z, double *values);

void ExactP(double x, double y,  double z, double *values);

/// Initial solution in U1 direction - Hagen-Poiseuille in inflow piece.
void InitialU1(double x, double y, double z, double *values);

/// Initial solution in U2 direction. No flow, nowhere.
void InitialU2(double x, double y, double z, double *values);

/// Initial solution in U3 direction. No flow, nowhere.
void InitialU3(double x, double y, double z, double *values);

/// Initial pressure. No pressure, nowhere.
void InitialP(double x, double y,  double z, double *values);

void BoundCondition(double x, double y, double z, BoundCond &cond);

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value);

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value);

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value);

void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs);

#endif /*TWISTED_PIPE_FLOW_*/
