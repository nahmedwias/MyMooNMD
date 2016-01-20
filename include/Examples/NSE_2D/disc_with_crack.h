void ExampleFile() {
    OutPut("Example: disc_with_crack.h" << endl);
}

double atan2_2pi(double y, double x) {
    double theta = atan2(y, x);
    return theta > 0 ? theta : theta + 2. * Pi;
}

void ExactU1(double x, double y, double *values) {
    double r = sqrt(x * x + y * y);
    double theta = atan2_2pi(y, x);
    values[0] = sqrt(r) * (cos(theta / 2.) - cos(3. * theta / 2.));
    values[1] = (x / (2. * r * sqrt(r))) * (cos(theta / 2.) - cos(3. * theta / 2.)) + (y * sqrt(r) / (2. * r * r)) * (sin(theta / 2.) - 3. * sin(3. * theta / 2.));
    values[2] = (y / (2. * r * sqrt(r))) * (cos(theta / 2.) - cos(3. * theta / 2.)) + (x * sqrt(r) / (2. * r * r)) * (3. * sin(3. * theta / 2.) - sin(theta / 2.));
    double u1_xx = (x * y / (2. * r * r * r * sqrt(r))) * (sin(theta / 2.) - 3. * sin(3. * theta / 2.)) + ((1 / (2. * r * sqrt(r))) - 3. * x * x / (4. * r * r * r * sqrt(r))) * (cos(theta / 2.) - cos(3. * theta / 2.))
                   + sqrt(r) * ((y * y * y * sin(theta / 2.) / (x * r * r * r * r)) - (3. * y * y * y * sin(3. * theta / 2.) / (x * r * r * r * r)) - (y * y * cos(theta / 2.) / (4. * r * r * r * r)) +
                                (9 * y * y * cos(3. * theta / 2.) / (4. * r * r * r * r)) - (y * sin(theta / 2.) / (x * r * r)) + (3. * y * sin(3. * theta / 2.) / (x * r * r)));
    double u1_yy = (x * y / (2. * r * r * r * sqrt(r))) * (3. * sin(3. * theta / 2.) - sin(theta / 2.)) + ((1 / (2. * r * sqrt(r))) - 3. * x * x / (4. * r * r * r * sqrt(r))) * (cos(theta / 2.) - cos(3. * theta / 2.))
                   + sqrt(r) *
                     ((-x * x * cos(theta / 2.) / (4. * r * r * r * r)) + (9 * x * x * cos(3. * theta / 2.) / (4. * r * r * r * r)) + (x * y * sin(theta / 2.) / (r * r * r * r)) - (3. * x * y * sin(3. * theta / 2.) / (r * r * r * r)));
    values[3] = u1_xx + u1_yy;
}

void ExactU2(double x, double y, double *values) {
    double r = sqrt(x * x + y * y);
    double theta = atan2_2pi(y, x);
    values[0] = sqrt(r) * (3. * sin(theta / 2.) - sin(3. * theta / 2.));
    values[1] = (x / (2. * r * sqrt(r))) * (3. * sin(theta / 2.) - sin(3. * theta / 2.)) + (3. * y * sqrt(r) / (2. * r * r)) * (cos(3. * theta / 2.) - cos(theta / 2.));
    values[2] = (y / (2. * r * sqrt(r))) * (3. * sin(theta / 2.) - sin(3. * theta / 2.)) + (3. * x * sqrt(r) / (2. * r * r)) * (cos(theta / 2.) - cos(3. * theta / 2.));
    double u2_xx = (1 / (2. * r * sqrt(r)) - 3. * x * x / (4. * r * r * r * sqrt(r))) * (3. * sin(theta / 2.) - sin(3. * theta / 2.)) + (3. * x * y / (2. * r * r * r * sqrt(r))) * (cos(theta / 2.) - cos(3. * theta / 2.))
                   + 3. * sqrt(r) *
                     ((y * y * y * cos(3. * theta / 2.) / (x * r * r * r * r)) + (3. * y * y * sin(3. * theta / 2.) / (4. * r * r * r * r)) - (y * cos(3. * theta / 2.) / (x * r * r)) - (y * y * y * cos(theta / 2.) / (x * r * r * r * r)) -
                      (y * y * sin(theta / 2.) / (4. * r * r * r * r)) + (y * cos(theta / 2.) / (x * r * r)));
    double u2_yy = (1 / (2. * r * sqrt(r)) - 3. * y * y / (4. * r * r * r * sqrt(r))) * (3. * sin(theta / 2.) - sin(3. * theta / 2.)) + (3. * x * y / (2. * r * r * r * sqrt(r))) * (cos(3. * theta / 2.) - cos(theta / 2.))
                   + 3. * sqrt(r) *
                     ((3. * x * x * sin(3. * theta / 2.) / (4. * r * r * r * r)) + (x * y * cos(3. * theta / 2.) / (r * r * r * r)) - (3. * x * x * sin(theta / 2.) / (4. * r * r * r * r)) - (3. * x * y * cos(theta / 2.) / (r * r * r * r)));
    values[3] = u2_xx + u2_yy;
}

void ExactP(double x, double y, double *values) {
    double r = sqrt(x * x + y * y);
    double theta = atan2_2pi(y, x);
    values[0] = (-4. / (sqrt(r))) * cos(theta / 2.);
    values[1] = (2. / (r * r * sqrt(r))) * (x * cos(theta / 2.) - y * sin(theta / 2.));
    values[2] = (2. / (r * r * sqrt(r))) * (x * sin(theta / 2.) - y * cos(theta / 2.));
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
    } else if (BdComp > 2) {
        cerr << "wrong boundary part number: " << BdComp << endl;
    }
}

void U2BoundValue(int BdComp, double Param, double &value) {
    value = 0;
    if (BdComp == 2) {
        // on circle
        value = (3. * sin(Pi * Param) - sin(3. * Pi * Param));
    } else if (BdComp > 2) {
        cerr << "wrong boundary part number: " << BdComp << endl;
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

        // check if derivatives make sense
        {
            ExactU1(X[i], Y[i], val1);
            ExactU2(X[i], Y[i], val2);
            ExactP(X[i], Y[i], val3);

            double f1 = -nu * val1[3] + val3[1];
            if (fabs(f1) > 1e16) {
                cerr << "wrong f1: " << f1 << endl;
            }
            double f2 = -nu * val2[3] + val3[2];
            if (fabs(f2) > 1e16) {
                cerr << "wrong f1: " << f2 << endl;
            }

            double div = val1[1] + val2[2];
            if (fabs(div) > 1e16) {
                cerr << "wrong div = " << div << endl;
            }
        }

        coeffs[i][0] = nu;
        coeffs[i][1] = 0; // f1
        coeffs[i][2] = 0; // f2
    }
}