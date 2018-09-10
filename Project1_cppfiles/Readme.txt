This program solves Poisson's equation -u''=f, and compares different methods in terms of speed
and in terms of relative error. The main program consits of 5 parts, which can be 
commented/uncommented according to needs. The integer n determines the number of grid points.

1) The first part prints the tuples of the form (x_i,u_i,v_i,w_i) to the terminal. Here x
is the given partition of [0,1], u is the analytic solution, v is the approximation to u 
using Thomas' method, and w is the special refinement of Thomas' method.

2) Assuming that n=4, the second part prints pairs (x_i, v_i) to the file "n_equals_4",
where v is the approximate solution to Poissons equation. The file "n_equals_4" will consist
of n+2 lines. This can for example be read by a python script to create a plot.

3) The third part prints the tuples (log10(j), log10(rel error (v)), log10(rel error(w)) 
to the file "relerror.txt", where now j is the number of grid points. We increase j from 10 to
100000000. Again u is the analytic solution, v is the approximation using
Thomas' method, and w is the special method.

4) The fourth part compares CPU times of the special and the general algorithm for the given n. We run 
10 experiments of each method. The result is printed to the terminal, and saved to the file
"CPU_time_Special_Vs_General_n_equals_4", if say n=4. The lines in this file consist of tuples
(special alg time, general alg time) meassured in seconds.

5) This part compares LU decomposition speed to solve -u''=f with the special algorithm 
for the given n. As with 4), we run 10 experiments and print the result to 
"LU_vs_Special_Algorithm_n_equals_4.txt", if n=4.


As an example, if we set n=4, the output of this program in the terminal is:
"x: u: v:  w:
0 0 0 0
0.2 0.664674 0.481265 0.481265
0.4 0.581703 0.421188 0.421188
0.6 0.397548 0.28785 0.28785
0.8 0.199701 0.144596 0.144596
1 5.88857e-18 0 0

Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0
Special algorithm:      0       General algorithm:      0

LU decomposition:       0.001   Special algorithm:      0
LU decomposition:       0       Special algorithm:      0
LU decomposition:       0       Special algorithm:      0
LU decomposition:       0       Special algorithm:      0
LU decomposition:       0       Special algorithm:      0
LU decomposition:       0       Special algorithm:      0
LU decomposition:       0       Special algorithm:      0
LU decomposition:       0       Special algorithm:      0
LU decomposition:       0       Special algorithm:      0
LU decomposition:       0       Special algorithm:      0

Success!"

Disclamer: part 4) and 5) do not alway produce the same results, since the time needed to
process these algorithms may vary.