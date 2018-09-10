#include <cmath>
double Error(double * approximation, double *absolute, int n) {
	double enumerator=0.0;
	double denominator=0.0;
	for (int i = 0; i < n; i++) {
		enumerator += (absolute[i]- approximation[i])*(absolute[i] - approximation[i]);
		denominator += absolute[i]* absolute[i];
	}
	return sqrt(enumerator / denominator);
};