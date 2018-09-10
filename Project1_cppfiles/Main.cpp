//This program solves Poisson's equation -u''=f, and compares different methods in terms of speed
//and in terms of relative error. The main program consits of 5 parts, which can be
//commented / uncommented according to needs. The integer n determines the number of grid points.
//
//1) The first part prints the tuples of the form (x_i, u_i, v_i, w_i) to the terminal.Here x
//is the given partition of[0, 1], u is the analytic solution, v is the approximation to u
//using Thomas' method, and w is the special refinement of Thomas' method.
//
//2) Assuming that n = 4, the second part prints pairs(x_i, v_i) to the file "n_equals_4",
//where v is the approximate solution to Poissons equation. The file "n_equals_4" will consist
//of n + 2 lines. This can for example be read by a python script to create a plot.
//
//3) The third part prints the tuples(log10(j), log10(rel error(v)), log10(rel error(w))
//to the file "relerror.txt", where now j is the number of grid points. We increase j 
//from 10 to 100000000. Again u is the analytic solution, v is the approximation using
//Thomas' method, and w is the special method.
//
//4) The fourth part compares CPU times of the special and the general algorithm for the given n.We run
//10 experiments of each method. The result is printed to the terminal, and saved to the file
//"CPU_time_Special_Vs_General_n_equals_4", if say n = 4. The lines in this file consist of tuples
//(special alg time, general alg time) meassured in seconds.
//
//5) This part compares LU decomposition speed to solve - u'' = f with the special algorithm
//for the given n.As with 4), we run 10 experiments and print the result to
//"LU_vs_Special_Algorithm_n_equals_4.txt", if n = 4.
//
//
//As an example, if we set n = 4, the output of this program in the terminal is:
//"x: u: v:  w:
//0 0 0 0
//0.2 0.664674 0.481265 0.481265
//0.4 0.581703 0.421188 0.421188
//0.6 0.397548 0.28785 0.28785
//0.8 0.199701 0.144596 0.144596
//1 5.88857e-18 0 0
//
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//Special algorithm : 0       General algorithm : 0
//
//LU decomposition : 0.001   Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//LU decomposition : 0       Special algorithm : 0
//
//Success!"
//
//Disclamer : part 4) and 5) do not alway produce the same results, since the time needed to
//process these algorithms may vary.

#include <iostream> //Print to command line
#include <fstream> //Pring to file
#include <string> //For using strings to write files
#include <cmath> //For using the exponential function
#include <ctime> //For meassuring time
#include <armadillo>
#include <vector>


using namespace arma;
using namespace std;

//Function that solves the general tridiagonal matrix equation Av=b and returns solution (See "TridiagonalMatrixEqSolve.cpp")
double * TriMaEqSo(double * a, double * b, double * c, double * d, int n);

//Function that solves the special tridiagonal matrix equation Av=b, where A represents a double derivative and returns solution (See "SpecialTriMaEqSo.cpp")
double * SpecialTriMaEqSo(double * d, double * beta, int n);

//Function that returns the absolute error of a vector and an approximation (See "Error.cpp")
double Error(double * approximation, double *absolute, int n);

//Function that returns the analytic solution
inline double uanalytic(double x) {return 1 - (1 - exp(-10))*x - exp(-10 * x);}

//Function that returns arrays x, a, b, c, d and beta given an n,
//and a function that deletes them afterwards (See ArrayManagement.cpp)
void GetArrays(double* &x, double* &a, double* &b, double* &c, double* &d, double* &beta, int n);
void DelArrays(double* &x, double* &a, double* &b, double* &c, double* &d, double* &beta);

//Main program
int main() {
	int n = 4; //determines the number of grid points
	
	double* x; double* a; double* b; double* c; double* d; double* beta; 
	
	double *v; //special tridiagmatrixeq solution 
	double *w; //general tridiagmatrixeq solution
	double *u; //analytic solution

	clock_t time1;
	clock_t time2;
	double difference;
	int laps; //number of experiments when comparing time consumption of different methods

	//=======================================================//
	//print the tuples (x_i, u_i, v_i, w_i) to terminal      //
	//=======================================================//

	GetArrays(x, a, b, c, d, beta, n);
	
	v = TriMaEqSo(a, b, c, d, n);
	w = SpecialTriMaEqSo(d, beta, n);

	cout << "x:" << " " << "u:" << " " << "v: " << " " << "w:"<<endl;
	for (int i = 0; i < n + 2; i++) {
		cout << x[i] << " " << uanalytic(x[i]) <<" " << v[i] <<  " " << w[i] << endl;
	}

	DelArrays(x, a, b, c, d, beta);
	delete[] v; delete[] w;
	cout << endl;

	//=========================================================//
	//print the pairs (x_i, v_i) to file "n_equals_n.txt"	   //
	//=========================================================//
	
	GetArrays(x, a, b, c, d, beta, n);
	v = TriMaEqSo(a,b,c,d, n);
	ofstream myfile;
	string filename("n_equals_");
	filename += to_string(n);
	filename += ".txt";
	myfile.open(filename);
	for (int i = 0; i < n + 2; i++) {
		myfile << x[i] <<" " << uanalytic(x[i]) << " " << v[i] << "\n";
	}
	myfile.close();
	DelArrays(x, a, b, c, d, beta);
	delete[] v;

	//==================================================================//
	//Computing relative error from approximation vs analytical solution//
	//and comparing the special tridiag matrix solver to the general one//
	//The result is printed to "relerror.txt"                           //
	//Each line in "relerror.txt" is of the form:						//
	//																	//
	// log10(h)  log10(abs error special)   log10(abs error general)    //
	//																	//
	//==================================================================//

	double errv; //special tridiag matrix error
	double errw; //general tridiag matrix error

	ofstream myfile2;
	string filename2("relerror.txt");
	myfile2.open(filename2);
	double h;

	for (int j = 10; j <= 10000000; j *=2) { //increase the number 10000000 at own risk!!
		GetArrays(x, a, b, c, d, beta, j);
		u = new double[j + 2];
		for (int i = 0; i < j + 2; i++) { u[i] = uanalytic(x[i]); };
		v = SpecialTriMaEqSo(d, beta, j);
		errv = Error(v, u, j + 2);
		w = TriMaEqSo(a, b, c, d, j);
		errw = Error(w, u, j + 2);
		delete[] u, v, w;
		DelArrays(x, a, b, c, d, beta);
		h = 1.0 / (j + 1);
		myfile2 << log10(h) <<" " <<log10(errv) << " " << log10(errw) << endl;
	}

	//======================================================================================//
	//Meassure time of general tridiagonal matrix eqsolver vs special tridiag matrix solver //
	//The result is printed to the command line window and saved to							//
	// "CPU_time_Special_Vs_General_n_equals_n.txt"											//
	//======================================================================================//

	GetArrays(x, a, b, c, d, beta, n);
	laps=10;

	ofstream myfile4;
	string filename4("CPU_time_Special_Vs_General_n_equals_");
	filename4 += to_string(n);
	filename4 += ".txt";
	myfile4.open(filename4);

	for (int i = 0; i < laps; i++) {
		time1 = clock();
		v = SpecialTriMaEqSo(d,beta,n);
		time2 = clock();
		difference = ((double) (time2 - time1))/((double) CLOCKS_PER_SEC);
		cout << "Special algorithm:	" << difference << "	";
		myfile4 << difference << " ";
		time1 = clock();
		delete[] v;
		v = TriMaEqSo(a,b,c,d, n);
		time2 = clock();
		difference = ((double)(time2 - time1)) / ((double)CLOCKS_PER_SEC);
		cout << "General algorithm:	" << difference << "	"<<endl;
		myfile4 << difference << endl;
		delete[] v;
	}
	cout << endl;
	DelArrays(x, a, b, c, d, beta);

	//===========================================================================//
	//Comparison of LU decomposition (armadillo solve) time vs Thomas' algorithm //
	//===========================================================================//

	ofstream myfile3;
	string filename3("LU_vs_Special_Algorithm_n_equals");
	filename3 += to_string(n);
	filename3 += ".txt";
	myfile3.open(filename3);


	laps = 10;
	for (int j = 1; j <= laps; j ++) {
		GetArrays(x, a, b, c, d, beta, n);
		sp_mat Asparse(n, n); mat L(n, n); mat U(n, n); vec v_arma(n); vec y_arma(n);
		for (int i = 0; i < n - 1; i++) { Asparse(i, i) = 2.0; Asparse(i, i + 1) = -1.0; Asparse(i + 1, i) = -1.0; }
		Asparse(n - 1, n - 1) = 2.0;
		mat A = conv_to<mat>::from(Asparse);
		vec d_arma(n);
		for (int i = 0; i < n; i++) { d_arma(i) = d[i + 1]; };
		
		time1 = clock();
		lu(L, U, A);
		//Solving the equation Ly=d
		y_arma = solve(L, d_arma);
		//Solving the equation Uv=y
		v_arma = solve(U, y_arma);
		time2 = clock();
		difference = ((double)(time2 - time1)) / ((double)CLOCKS_PER_SEC);
		cout << "LU decomposition:	"  << difference << "	";
		myfile3  << difference <<" ";

		//Special tridiagonal matrix equation
		time1 = clock();
		v = SpecialTriMaEqSo(d, beta, j);
		time2 = clock();
		difference = ((double)(time2 - time1)) / ((double)CLOCKS_PER_SEC);
		cout << "Special algorithm:	" <<  difference << "	" << endl;
		myfile3  << difference << endl;

		DelArrays(x, a, b, c, d, beta);
		delete[] v;
	}
	cout << endl;

	//====================================================//
	//Success message if the program runs.				  //
	//====================================================//

	cout << "Success!";
	cin.get();
}