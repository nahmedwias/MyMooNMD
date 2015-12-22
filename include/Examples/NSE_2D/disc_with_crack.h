void ExampleFile() {
    OutPut("Example: disc_with_crack.h" << endl);
}

void ExactU1(double x, double y, double *values) {
    double r = sqrt(x * x + y * y);
    double theta = atan(y / x);
    values[0] = (3 * sqrt(r) / 2) * (cos(theta / 2) - cos(3 * theta / 2));
    values[1] = 3 * sin(theta / 2) * (-x * y * (sqrt((x * x) / (y * y) + 1) + 2)) / (2 * r * r * r * r * sqrt(r));
    values[2] = 3 * sin(3 * theta / 2) * (y * y + 3 * x * x + x * r) / (2 * r * r * r * r * sqrt(r));
    values[3] = 0;
}

void ExactU2(double x, double y, double *values) {
    double r = sqrt(x * x + y * y);
    double theta = atan(y / x);
    values[0] = (3 * sqrt(r) / 2) * (3 * sin(theta / 2) - sin(3 * theta / 2));
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactP(double x, double y, double *values) {
    double r = sqrt(x * x + y * y);
    double theta = atan(y / x);
    values[0] = -6 / (sqrt(r)) * cos(theta / 2);
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond) {
    cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value) {
    value = 0;
    if (BdComp == 2) {
        // on circle
        value = (cos(Pi * Param) - cos(3. * Pi * Param));
    } else if(BdComp > 2) {
        cerr << "wrong boundary part number: " << BdComp << endl;
    }
}

void U2BoundValue(int BdComp, double Param, double &value) {
    value = 0;
    if (BdComp == 2) {
        // on circle
        value = (3. * sin(Pi * Param) - sin(3. * Pi * Param));
    } else if(BdComp > 2) {
        cerr << "wrong boundary part number: " << BdComp << endl;
    }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs) {
    const double nu = 1. / TDatabase::ParamDB->RE_NR;
    for (int i = 0; i < n_points; i++) {
        coeffs[i][0] = nu;
        coeffs[i][1] = 0; // f1
        coeffs[i][2] = 0; // f2
    }
}