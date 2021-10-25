#include <stdio.h>
#include <math.h>

// x' = a*x - b*xy
// y' = c*xy - d*y
// a = d = 10, b = c = 2

double f_func(double x, double y, double a = 10, double b = 2){
    return a*x - b*x*y;
}

double f_func_shx(double x, double y, double a = 10, double b = 2){
    return a - b*y;
}

double f_func_shy(double x, double y, double a = 10, double b = 2){
    return - b*x;
}

double g_func(double x, double y, double c = 2, double d = 10){
    return c*x*y - d*y;
}

double g_func_shx(double x, double y, double c = 2, double d = 10){
    return c*y;
}

double g_func_shy(double x, double y, double c = 2, double d = 10){
    return c*x - d;
}

//x(n) = x(n-1) + h*( 0.25*f(t, x(n-1), y(n-1)) + 0.75*f(t + 2h/3, x(n-1) + 2h/3 * f(t, x(n-1), y(n-1)), y(n-1)) )
//y(n) = y(n-1) + h*( 0.25*g(t, x(n-1), y(n-1)) + 0.75*g(t + 2h/3, x(n-1), y(n-1) + 2h/3 * g(t, x(n-1), y(n-1))))



void runge_kutta_second_solver(double t, int n, // n - iteration number 
    double(*f)(double, double, double, double), double(*g)(double, double, double, double),
    double(*fx)(double, double, double, double), double(*gx)(double, double, double, double),
    double(*fy)(double, double, double, double), double(*gy)(double, double, double, double)
    ){

        double u = 1, v = 1;
        double h = 1 / (double)n ;

        for(int i = 0; i < n; i++){
            u = h*f(u, v, 2, 10) + h*fx(u, v, 2, 10) + h*fy(u, v, 2, 10);
            v = h*g(u, v, 2, 10) + h*gx(u, v, 2, 10) + h*gy(u, v, 2, 10);
        } 

}

int main(){

    return 0;
}