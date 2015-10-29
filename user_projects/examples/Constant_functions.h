/**
 * @brief Prints information about the used example.
 */
void ExampleFile()
{
	OutPut("Example: Two coupled Convection-Diffusion-Reaction equations\n");
}

/**
 * @brief First coupling term, in the LaTeX file referred to as \f$ R_1(c_1, c_2) \f$ .
 *
 * @param c[0] function @a c1 evaluated at quadrature point @a x
 * @param c[1] function @a c2 evaluated at quadrature point @a x
 * @return sum of c[0] and c[1]
 */
double couplingTerm1(const double* const c)
{
	return c[0] + c[1];
}

/**
 * @brief Second coupling term, in the LaTeX file referred to as \f$ R_2(c_1, c_2) \f$ .
 *
 * @param c[0] function @a c1 evaluated at quadrature point @a x
 * @param c[1] function @a c2 evaluated at quadrature point @a x
 * @return product of c[0] and c[1]
 */
double couplingTerm2(const double* const c)
{
	return c[0] * c[1];
}

/**
 * @brief exact solution \f$ c_1 \f$ and its derivatives at (x, y)
 *
 * The values get sorted as follows:
 *
 * values[0] = \f$ c_1(x, y) \f$
 * values[1] = \f$ \dfrac{\partial c_1(x, y)}{\partial x} \f$
 * values[2] = \f$ \dfrac{\partial c_1(x, y)}{\partial y} \f$
 * values[3] = \f$ \Delta c_1(x, y) \f$
 *
 * @param[in] x abscissa
 * @param[in] y ordinate
 * @param[out] values evaluated function and derivatives
 */
void ExactC1(double x, double y, double *values)
{
	values[0] = 0.25;     	// value
	values[1] = 0.;     	// x-derivative
	values[2] = 0.;    		// y-derivative
	values[3] = 0.;			// Laplace
}

/**
 * @brief exact solution \f$ c_2 \f$ and its derivatives at (x, y)
 *
 * The values get sorted as follows:
 *
 * values[0] = \f$ c_2(x, y) \f$
 * values[1] = \f$ \dfrac{\partial c_2(x, y)}{\partial x} \f$
 * values[2] = \f$ \dfrac{\partial c_2(x, y)}{\partial y} \f$
 * values[3] = \f$ \Delta c_2(x, y) \f$
 *
 * @param[in] x abscissa
 * @param[in] y ordinate
 * @param[out] values evaluated function and derivatives
 */
void ExactC2(double x, double y, double *values)
{
	values[0] = 0.5;     	// value
	values[1] = 0.;     	// x-derivative
	values[2] = 0.;    		// y-derivative
	values[3] = 0.;			// Laplace
}

// ============================================================================
// boundary conditions
// ============================================================================
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
	cond = DIRICHLET;
}

// ============================================================================
// value of boundary condition for c1
// ============================================================================
void BoundValueC1(int BdComp, double Param, double &value)
{
	value = 0.25;
}

// ============================================================================
// value of boundary condition for c2
// ============================================================================
void BoundValueC2(int BdComp, double Param, double &value)
{
	value = 0.5;
}

/**
 * @brief Gives the coefficients of \f$ c_1 \f$ for the bilinear form.
 *
 * i represents the i-th quadrature point
 * coeffs[i][0] = diffusion coefficient
 * coeffs[i][1] = x coordinate for convection coefficient
 * coeffs[i][2] = y coordinate for convection coefficient
 * coeffs[i][3] = reaction coefficient
 * coeffs[i][4] = right hand side
 *
 * @param[in] n_points number of quadrature points
 * @param[in] X abscissas
 * @param[in] Y ordinates
 * @param[in] parameters - not needed here.
 * @param[out] coeffs returned coefficients
 */
void BilinearCoeffsC1(int n_points, double *X, double *Y, double **parameters,
		double **coeffs)
{
	double *valuesC1 = new double[4];
	double *valuesC2 = new double[4];

	for (int i = 0; i < n_points; i++)
	{
		ExactC1(X[i], Y[i], valuesC1);
		ExactC2(X[i], Y[i], valuesC2);

		coeffs[i][0] = 1.;			// diffusion
		coeffs[i][1] = 1.;			// convection x
		coeffs[i][2] = 1.;			// convection y
		coeffs[i][3] = 0.;			// reaction

		//Prepare input for the coupling function, which expects a double array.
		double couplingArray1[2] = {valuesC1[0], valuesC2[0]};

		coeffs[i][4] = coeffs[i][0]*valuesC1[3]
					 + coeffs[i][1]*valuesC1[1] + coeffs[i][2]*valuesC1[2]
					 + coeffs[i][3]
					 + couplingTerm1(couplingArray1); //right hand side

	}

	delete [] valuesC1;
	delete [] valuesC2;
}

/**
 * @brief Gives the coefficients of \f$ c_2 \f$ for the bilinear form.
 *
 * i represents the i-th quadrature point
 * coeffs[i][0] = diffusion coefficient
 * coeffs[i][1] = x coordinate for convection coefficient
 * coeffs[i][2] = y coordinate for convection coefficient
 * coeffs[i][3] = reaction coefficient
 * coeffs[i][4] = right hand side
 *
 * @param[in] n_points number of quadrature points
 * @param[in] X abscissas
 * @param[in] Y ordinates
 * @param[in] parameters - not needed here
 * @param[out] coeffs returned coefficients
 */
void BilinearCoeffsC2(int n_points, double *X, double *Y, double **parameters,
		double **coeffs)
{
	double *valuesC1 = new double[4];
	double *valuesC2 = new double[4];

	for (int i = 0; i < n_points; i++)
	{
		ExactC1(X[i], Y[i], valuesC1);
		ExactC2(X[i], Y[i], valuesC2);

		coeffs[i][0] = 1.;			// diffusion
		coeffs[i][1] = 1.;			// convection x
		coeffs[i][2] = 1.;			// convection y
		coeffs[i][3] = 0.;			// reaction

		//Prepare input for the coupling function, which expects a double array.
		double couplingArray1[2] = {valuesC1[0], valuesC2[0]};

		coeffs[i][4] = coeffs[i][0]*valuesC2[3]
					 + coeffs[i][1]*valuesC2[1] + coeffs[i][2]*valuesC2[2]
					 + coeffs[i][3]
					 + couplingTerm2(couplingArray1);		// right hand side
	}

	delete [] valuesC1;
	delete [] valuesC2;
}


/****************************************************************************************
 * ParameterFunction and AssemblingFunctions used in the "Linearized Decoupled" strategy.
 ****************************************************************************************/

void ParameterFunction(double* in, double* out){
	// Just skip the first two entries - this is where TAuxParam2D->GetParameters() places the x and y value,
	// and pass on as many values as nCoupled_ - for this example that would be 2.
	for (size_t i = 0; i<2 ; ++i){
		out[i] = in [i+2];
	}
}

void AssemblingFunctionC1(
		double Mult, double *coeff, double *param,
		double hK, double **OrigValues, int *N_BaseFuncts,
		double ***LocMatrices, double **LocRhs)
{
	// For convenience: rename the place to write to (Rhs) and the place to read from (Orig,
	// which suposedly contains values of test functions)
	double* Rhs = LocRhs[0];
	double* Orig = OrigValues[0];

	// Calculate the value of the coupling term at the quad point.
	// This relies on the correct interaction with the ParameterFunction and the
	// input order of fe functions to the aux object.

	// coupled term evaluated at the quad point this
	// AssembleFctParam2D is called upon
	double coupledTerm= couplingTerm1(param);

	//Loop over all local base functions.
	for(int i=0;i<N_BaseFuncts[0];i++)
	{
		Rhs[i] -= Mult*coupledTerm*Orig[i]; //(Mult contains quad weigth, and maybe even the determinant due to trafo)
	}
}

void AssemblingFunctionC2(
		double Mult, double *coeff, double *param,
		double hK, double **OrigValues, int *N_BaseFuncts,
		double ***LocMatrices, double **LocRhs)
{
	// For convenience: rename the place to write to (Rhs) and the place to read from (Orig,
	// which suposedly contains values of test functions)
	double* Rhs = LocRhs[0];
	double* Orig = OrigValues[0];

	// Calculate the value of the coupling term at the quad point.
	// This relies on the correct interaction with the ParameterFunction and the
	// input order of fe functions to the aux object.

	// coupled term evaluated at the quad point this
	// AssembleFctParam2D is called upon
	double coupledTerm= couplingTerm2(param);

	//Loop over all local base functions.
	for(int i=0;i<N_BaseFuncts[0];i++)
	{
		Rhs[i] -= Mult*coupledTerm*Orig[i]; //(Mult contains quad weigth, and maybe even the determinant due to trafo)
	}
}
