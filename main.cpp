#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <array>
#include <algorithm>
#include <sstream>


const double mesh_step = 0.01;


std::random_device rd;  // Will be used to obtain a seed for the random number engine.
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd().


double trial_function (double & x);

std::vector<double> function_nodes (std::vector<double> & xx, double f(double & x));

std::vector<double> random_shift (std::vector<double> f, const double & max_shift);

double energy (std::vector<double> & x, std::vector<double> & f);


// Technical functions.

std::vector <double> mesh (double left_border, const double & right_border, const double & step);

void data_file_creation (const std::string & name, std::vector<double> & x, std::vector<double> & y);

void plot (const std::string & name, const int & left, const int & right,
           const std::string & title, const std::string & xlabel, const std::string & ylabel);



int main () {
    std::vector<double> x_mesh = std::move(mesh(-10, 10, mesh_step));
    std::vector<double> f_mesh = std::move(function_nodes(x_mesh, trial_function));
    data_file_creation("trial_function", x_mesh, f_mesh);
    plot("trial_function", -10, 10, "Trial function", "x", "psi");
    double E, E_GS = energy(x_mesh, f_mesh);
    do {
        std::vector<double> f_buf = std::move(random_shift(f_mesh, 0.1));
        E = energy(x_mesh, f_buf);
        if (E <= E_GS) {
            E_GS = E;
            f_mesh = std::move(f_buf);
        }
        std::cout << E_GS << std::endl;
    } while (E_GS > 0.5);
    data_file_creation("result_function", x_mesh, f_mesh);
    plot("result_function", -10, 10, "Result function", "x", "psi");
    return 0;
}


// Returns table-function (std::vector<double> f) with random shifted (dis_shift) random point (dis_num).
std::vector<double> random_shift (std::vector<double> f, const double & max_shift) {
    std::uniform_int_distribution<> dis_num (1, f.size()-2);
    std::uniform_real_distribution<double> dis_shift (-max_shift, max_shift);
    f[dis_num(gen)] += dis_shift(gen);
    return f;
}


// Returns the numeric value of 2nd derivative in node i.
double diff2 (std::vector<double> & f, int & i, const double & dx) {
    if (i == 0 || i == f.size()-1)
        return 0;
    else
        return (f[i-1] - 2.0*f[i] + f[i+1]) / std::pow(dx, 2);
}


// Returns harmonic oscillator potential V = 1/2 x^2.
double potential (double & x) {
    return std::pow(x, 2) / 2.0;
}


// Returns <f|H|f>
double mean_hamiltonian (std::vector<double> & x, std::vector<double> & f) {
    double sum = 0;
    for (int i = 0; i < x.size(); ++i)
        sum += f[i] * ((-diff2(f, i, mesh_step) / 2.0) + potential(x[i]) * f[i]); // <f|H|f> = <f|-1/2 d2/dx2 + x2/2|f>
    return sum;
}


// Returns <f|f>
double normalization (std::vector<double> & function) {
    double sum = 0;
    for (double & f : function)
        sum += std::pow(f, 2);
    return sum;
}


// Returns E = <f|H|f> / <f|f>
double energy (std::vector<double> & x, std::vector<double> & f) {
    return mean_hamiltonian(x, f) / normalization(f);
}


// Returns the trial function value at point x.
double trial_function (double & x) {
    return (std::fabs(x) <= 0.5) ? 1 : 0;
}


// Returns std::vector of function values according to the given argument mesh (xx).
std::vector<double> function_nodes (std::vector<double> & xx, double f(double & x)) {
    std::vector<double> ff (xx.size());
    for (int i = 0; i < xx.size(); ++i)
        ff[i] = f(xx[i]);
    return ff;
}


// There are technical functions that are not related to the solution below.

// Creates mesh from left border to right border with given step.
std::vector <double> mesh (double left_border, const double & right_border, const double & step) {
    std::vector <double> xx ((right_border-left_border) / step);
    xx[0] = left_border;
    std::generate(xx.begin()+1, xx.end(), [&] {left_border += step; return left_border;});
    return xx;
}


// Returns std::string from numeric type (std::to_string not safe enough).
template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


// Creates text file (.txt) from with coordinates (x, y) for plotting.
void data_file_creation (const std::string & name, std::vector<double> & x, std::vector<double> & y) {
    std::ofstream fout;
    fout.open(name + ".txt", std::ios::out | std::ios::trunc);
    for (int i = 0; i < x.size(); ++i)
        fout << toString(x[i]) << '\t' << toString(y[i]) << std::endl;
    fout.close();
}


// Creates plot from text file with set parameters.
void plot (const std::string & name, const int & left, const int & right,
           const std::string & title, const std::string & xlabel, const std::string & ylabel) {
    std::string range = "[" + toString(left) + ":" + toString(right) + "]";
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 700, 700",
                                      "set output \'" + name + ".jpg\'",
                                      "set title \'" + title + "\'",
                                      "set grid xtics ytics",
                                      "set xrange " + range,
                                      "set xlabel \'" + xlabel + "\'",
                                      "set ylabel \'" + ylabel + "\'",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "plot \'" + name + ".txt\' using 1:2 w lines",
                                      "set terminal pop",
                                      "set output",
                                      "replot", "q"};
    for (const auto & it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}
