#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <sstream>


const double mesh_step = 0.01;


typedef std::pair<double, double> data;


std::random_device rd;  // Will be used to obtain a seed for the random number engine.
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd().


double trial_function (double & x);

std::vector <double> mesh (double left_border, const double & right_border, const double & step);

std::vector<double> function_nodes (std::vector<double> & xx, double f(double & x));

std::vector<double> random_shift (std::vector<double> f, const double & max_shift);

double Energy (std::vector<double> & x, std::vector<double> & f);


// Technical functions.

std::vector<data> data_collection (std::vector<double> & x, std::vector<double> & f);

void data_file_creation (const std::string & name, const std::vector<data> &exp_data);

void plot (const std::string & name, const int & left, const int & right,
           const std::string & title, const std::string & xlabel, const std::string & ylabel);



int main () {
    std::vector<double> x_mesh = std::move(mesh(-10, 10, mesh_step));
    std::vector<double> f_mesh = std::move(function_nodes(x_mesh, trial_function));
    data_file_creation("trial_function", data_collection(x_mesh, f_mesh));
    plot("trial_function", -10, 10, "Trial function", "x", "psi");
    double E, E_GS = Energy(x_mesh, f_mesh);
    do {
        std::vector<double> f_buf = std::move(random_shift(f_mesh, 0.1));
        E = Energy(x_mesh, f_buf);
        if (E <= E_GS) {
            E_GS = E;
            f_mesh = std::move(f_buf);
        }
        std::cout << E_GS << std::endl;
    } while (E_GS > 0.5);
    data_file_creation("result_function", data_collection(x_mesh, f_mesh));
    plot("result_function", -10, 10, "Result function", "x", "psi");
}


// Returns table-function (std::vector<double> f) with random shifted (shift) random point (num).
std::vector<double> random_shift (std::vector<double> f, const double & max_shift) {
    std::uniform_int_distribution<> dis_num (1, f.size()-2);
    std::uniform_real_distribution<double> dis_shift (-max_shift, max_shift);
    int num = dis_num(gen);
    double shift = dis_shift(gen);
    f[num] += shift;
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
    for (double f : function)
        sum += std::pow(f, 2);
    return sum;
}


// Returns E = <f|H|f> / <f|f>
double Energy (std::vector<double> & x, std::vector<double> & f) {
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


// Creates mesh from left border to right border with given step.
std::vector <double> mesh (double left_border, const double & right_border, const double & step) {
    std::vector <double> xx ((right_border-left_border) / step);
    std::generate(xx.begin(), xx.end(), [&] {left_border += step; return left_border;});
    return xx;
}


// There are technical functions that are not related to the solution below.


std::vector<data> data_collection (std::vector<double> & x, std::vector<double> & f) {
    std::vector<data> collected (x.size());
    for (int i = 0; i < x.size(); ++i)
        collected[i] = std::move(std::make_pair(x[i], f[i]));
    return collected;
}


// Returns std::string from numeric type (std::to_string not safe enough).
template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


// Returns std::string with numeric tuple in it.
template<typename T, size_t... Is>
std::string tuple_to_string_impl (T const& t, std::index_sequence<Is...>) {
    return ((toString(std::get<Is>(t)) + '\t') + ...);
}

template <class Tuple>
std::string tuple_to_string (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return tuple_to_string_impl(t, std::make_index_sequence<size>{});
}


// Creates text file (.txt) with given name from std::vector of numeric tuples.
void data_file_creation (const std::string & name, const std::vector<data> &exp_data) {
    std::ofstream fout;
    fout.open(name + ".txt", std::ios::out | std::ios::trunc);
    for (auto & i : exp_data)
        fout << tuple_to_string(i) << std::endl;
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