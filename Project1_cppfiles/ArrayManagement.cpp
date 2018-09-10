#include<cmath>

//Function that represents f in the right hand side (RHS) of -u''=f
inline double f(double x) { return 100 * exp(-10 * x); }

void GetArrays(double* &x, double* &a, double* &b, double* &c, double* &d, double* &beta, int n) {
	double h = 1.0 / (n + 1); //step size
	double hh = h * h;

	//In the equation Av=d, the vectors a, b, and c are vectors that determine the tridiagonal matrix A.
	//x is the partition of [0,1]
	x = new double[n + 2]; a = new double[n + 2]; b = new double[n + 2]; c = new double[n + 2]; d = new double[n + 2];
	for (int i = 0; i < n + 2; i++) { a[i] = -1.0;  c[i] = -1.0; b[i] = 2.0; };

	//Creating beta (the updated diagonal) for the special case when A represents u''
	beta = new double[n + 2];
	beta[0] = 0.0;
	for (int i = 1; i < n + 2; i++) { beta[i] = (i + 1.0) / i; }

	//Creating the vector d that represents the right hand side in Av=d,
	//and a parition of the unit interval [0,1] represented by x_i.
	for (int i = 0; i < n + 2; i++) {
		x[i] = i * h;
		d[i] = f(x[i])*hh;
	}
}

void DelArrays(double* &x, double* &a, double* &b, double* &c, double* &d, double* &beta) {
	delete[] x; delete[] d; delete[] b; delete[] a; delete[] c;  delete[] beta;
}