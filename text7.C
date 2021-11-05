#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double standart_grid(int i, int n, double a, double b){
    double x = (double)i / (double)n;
    return (b - a)*x + a;
}

double f_func(double x, double y, double a = 10, double b = 2){
    return a*x - b*x*y;
}

double g_func(double x, double y, double c = 2, double d = 10){
    return c*x*y - d*y;
}

void runge_kutta_second_solver(
    int i, int n, // n - iteration number
    double *x, double *y,
    double(*f)(double, double, double, double),
    double(*g)(double, double, double, double),        
    double(*grid)(int, int, double, double),
        
    double a = 10, double b = 2, double c = 2, double d = 10,
    double x1 = 0, double x2 = 1       
    ){
        
        // x' = a*x - b*xy
        //y' = c*xy - d*y
        // a = d = 10, b = c = 2

        double u = 1.0;
        double v = 1.0;
        //double h = 1 / (double)n ;

        for(int j = 0; j < i; j++){
            
            double xr = grid(i, n, x1, x2);
            double xl = grid(i-1, n, x1, x2);
            
            double u_buff = u;
            double v_buff = v;
            double h = x2 - x1;

            // x(n) = x(n-1) + h*( 0.25*f(t, x(n-1), y(n-1)) + 0.75*f(t + 2h/3, x(n-1) + 2h/3 * f(t, x(n-1), y(n-1)), y(n-1)) )
            // y(n) = y(n-1) + h*( 0.25*g(t, x(n-1), y(n-1)) + 0.75*g(t + 2h/3, x(n-1), y(n-1) + 2h/3 * g(t, x(n-1), y(n-1))))

            u += h*( 0.25*f(u_buff, v_buff, a, b) + 0.75*f(u_buff + 2*h*f(u_buff, v_buff, a, b) / 3.0, v_buff, a, b) );
            v += h*( 0.25*g(u_buff, v_buff, c, d) + 0.75*g(u_buff, v_buff + 2*h*g(u_buff, v_buff, c, d) / 3.0, c, d) );
        } 

        *x = u;
        *y = v;

}

int main(){

    double *xf = (double*) malloc(sizeof(double));
    double *yf = (double*) malloc(sizeof(double));

    int n = 10;
    
    for(int i = 1; i <= n; i++){
    
        runge_kutta_second_solver(i, n, xf, yf, f_func, g_func, standart_grid);
        printf("x = %f, y = %f\n", xf[0], yf[0]);
    
    }
    
    if(xf!=NULL){
        free(xf);
    }

    if(yf!=NULL){
        free(yf);
    }

    return 0;
}