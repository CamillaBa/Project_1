double * SpecialTriMaEqSo(double * d, double * beta, int n) {
	//Declaring delta and v
	double *delta = new double[n+2];
	double *v = new double[n+2];

	//Algorithm 1 (row reduction)
	//Declares new values for delta & beta
	delta[1] = d[1];
	for (int i = 2; i < n+1; i++) {
		delta[i]=d[i] + float(delta[i-1])/beta[i-1];
	}

	//Algorithm 2 (Solve linear equations)
	v[n] =  delta[n] / beta[n];
	for (int i = n-1; i > 0; i--) {
		v[i]=(delta[i] + v[i+1]) / beta[i];
	}

	//Adding the boundary values to v
	v[0] = 0.0;
	v[n + 1] = 0.0;
	delete[] delta;
	return v;
}