double * TriMaEqSo(double * a, double * b, double * c, double * d,int n) {
		double * beta= new double [n+2];
		double * delta = new double [n+2];
		double epsilon_iminus1;
	
		//Algorithm 1 (row reduction)
		//Decleares new values for beta and delta
		beta[1] = b[1];
		delta[1] = d[1];
		for (int i = 2; i < n+2; i++) {
			epsilon_iminus1 = a[i-1] / beta[i-1];
			beta[i]=b[i] - epsilon_iminus1 * c[i - 1];
			delta[i] = d[i] - epsilon_iminus1 * delta[i - 1];
		}
	
		//Algorithm 2 (Solve linear equations)
		double * v = new double [n+2];
		v[n] = delta[n]/ beta[n];
		for (int i = n-1; i > 0; i--) {
			v[i]=(delta[i] - c[i] * v[i+1]) / beta[i];
		}

		//Adding the boundary values to v
		v[0]=0.0;
		v[n+1]=0.0;

		delete[] beta; delete[] delta;
		return v;
	}