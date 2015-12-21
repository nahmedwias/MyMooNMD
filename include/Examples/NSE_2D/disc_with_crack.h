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
    if(BdComp == 2) {
        // on circle
        double t = 2*Pi*Param;
        value = (3./2.) * (cos(t/2.0) - cos(3.*t/2.));
    } else if(BdComp == 1) {
        double t = sqrt(Param);
        value = (3./2.) * t * (cos(0./2.) - cos(3. * 0./2.));
    } else if(BdComp == 0) {
        // second part of crack
        double t = sqrt(Param);
        double two_Pi_min_eps = 2.*Pi - 0.00001;
        value = (3./2.) * t * (cos(two_Pi_min_eps/2.) - cos(3. * two_Pi_min_eps/2.));
    } else {
        cerr << "wrong boundary part number: " << BdComp << endl;
    }
}

void U2BoundValue(int BdComp, double Param, double &value) {
    value = 0;
    if(BdComp == 2) {
        // on circle
        double t = 2*Pi*Param;
        value = (3./2.) * (3.*sin(t/2.) - sin(3.*t/2.));
    } else if(BdComp == 1) {
        double t = sqrt(Param);
        value = (3. / 2.) * t * (3. * sin(0./2.) - sin(3. * 0./2.));
    } else if(BdComp == 0) {
        double t = sqrt(Param);
        double two_Pi_min_eps = 2.*Pi - 0.00001;
        value = (3. / 2.) * t * (3. * sin(two_Pi_min_eps/2.) - sin(3. * two_Pi_min_eps/2.));
    }
    if (BdComp > 2) cout << "wrong boundary part number: " << BdComp << endl;
}

void PBoundValue(int BdComp, double Param, double &value) {
    value = 0;
    if(BdComp == 2) {
        double t = 2*Pi*Param;
        value = -6. * cos(t/2.);
    } else if(BdComp == 1) {
        double t = sqrt(Param);
        if(t > 0) {
            value = (-6. / t) * cos(0. / 2.);
        } else {
            value = 0;
        }
    } else if(BdComp == 0) {
        double t = sqrt(Param);
        double two_Pi_min_eps = 2.*Pi - 0.00001;
        if(t > 0) {
            value = (-6. / t) * cos(two_Pi_min_eps / 2.);
        } else {
            value = 0;
        }
    }
    if(BdComp > 2) {
        cout << "wrong boundary part number: " << BdComp << endl;
    }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs) {
    const double nu = 1. / TDatabase::ParamDB->RE_NR;
    double val1[4];
    double val2[4];
    double val3[4];
    for (int i = 0; i < n_points; i++) {
        coeffs[i][0] = nu;

        ExactU1(X[i], Y[i], val1);
        ExactU2(X[i], Y[i], val2);
        ExactP(X[i], Y[i], val3);

        coeffs[i][1] = -nu * val1[3] + val3[1]; // f1
        coeffs[i][2] = -nu * val2[3] + val3[2]; // f2

        if (TDatabase::ParamDB->PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
        {
            coeffs[i][1] += val1[0] * val1[1] + val2[0] * val1[2]; // f1
            coeffs[i][2] += val1[0] * val2[1] + val2[0] * val2[2]; // f2
        }
        coeffs[i][3] = val1[1] + val2[2]; // g (divergence)
    }

}