#include <stdio.h>
#include <math.h>

double minus_func(double x){
    return -x;
}

double standart_grid(int i, int n, double a, double b){
    double x = (double)i / (double)n;
    return (b - a)*x + a;
}

double euler_solver( int i, int n, double a, double b, double(*grid)(int, int, double, double), double(*f)(double) ){

    double y = 1.0;
    //x = a+ ( b - a ) * ( (double) i / (double) n );
    for(int j = 0; j < i; j++){
        double x2 = grid(i, n, a, b);
        double x1 = grid(i-1, n, a, b);
        double h = x2 - x1;
        y += h * f( y );
    }

    return y;

}

double runge_kutta_second( int i, int n, double a, double b, double(*grid)(int, int, double, double), double(*f)(double) ){

    double y = 1.0;

    for(int j = 0; j < i; j++){
        double x2 = grid(i, n, a, b);
        double x1 = grid(i-1, n, a, b);
        double h = x2 - x1;
        y += h * ( 0.25 * f( y ) + 0.75 * f( y + 2 * h * f( y ) / 3.0 ));
    }

    return y;
}

double runge_kutta_fourth(int i, int n, double a, double b, double (*grid)(int, int, double, double), double(*f)(double) ){
    double y = 1.0;

    for(int j = 0; j < i; j++){
        double x2 = grid(i, n, a, b);
        double x1 = grid(i-1, n, a, b);
        double h = x2 - x1;
        
        double k1 = f(y);
        double k2 = f(y + 0.5 * k1 * h);
        double k3 = f(y + 0.5 * k2 * h);
        double k4 = f(y + h * k3);

        y += h * ( k1 + 2 * ( k2 + k3 ) + k4 ) / 6.0;
    }

    return y;

}

int main(){

    double a = 0.0;
    double b = 3.0;
    int n = 4097;
    
    FILE *myfile1;
    const char filenam1[12] = {"out6de.txt"};
    myfile1 = fopen (filenam1, "w");

    for(int i = 1; i <= n; i*=2){
        // double x = (b - a) * ((double) i / (double) n);
        fprintf(myfile1, "%d %.15f\n", i, fabs(euler_solver(i*2, i*2, a, b, standart_grid, minus_func) - euler_solver(i, i, a, b, standart_grid, minus_func)));

    }

    fclose(myfile1);

    FILE *myfile2;
    const char filenam2[12] = {"out6ds.txt"};
    myfile2 = fopen (filenam2, "w");

    for(int i = 1; i <= n; i*=2){
        //double x = (b - a) * ((double) i / (double) n);
        fprintf(myfile2, "%d %.15f\n", i, fabs(runge_kutta_second(i*2, i*2, a, b, standart_grid, minus_func) - runge_kutta_second(i, i, a, b, standart_grid, minus_func)));

    }

    fclose(myfile2);

    FILE *myfile3;
    const char filenam3[12] = {"out6df.txt"};
    myfile3 = fopen (filenam3, "w");

    for(int i = 1; i <= n; i*=2){
        //double x = (b - a) * ((double) i / (double) n);
        fprintf(myfile3, "%d %.15f\n", i, fabs(runge_kutta_fourth(i*2, i*2, a, b, standart_grid, minus_func) - runge_kutta_fourth(i, i, a, b, standart_grid, minus_func)));

    }

    fclose(myfile3);

///////////////

    int m = 10;
    FILE *myfile4;
    const char filenam4[12] = {"out6eг.txt"};
    myfile4 = fopen (filenam4, "w");

    for(int i = 1; i <= m; i++){
        double x = (b - a) * ((double) i / (double) m);
        fprintf(myfile4, "%.15f %.15f\n", x, euler_solver(i, m, a, b, standart_grid, minus_func));

    }

    fclose(myfile4);

    FILE *myfile5;
    const char filenam5[12] = {"out6ds.txt"};
    myfile2 = fopen (filenam2, "w");

    for(int i = 1; i <= n; i*=2){
        //double x = (b - a) * ((double) i / (double) n);
        fprintf(myfile2, "%d %.15f\n", i, fabs(runge_kutta_second(i*2, i*2, a, b, standart_grid, minus_func) - runge_kutta_second(i, i, a, b, standart_grid, minus_func)));

    }

    fclose(myfile2);

    FILE *myfile3;
    const char filenam3[12] = {"out6df.txt"};
    myfile3 = fopen (filenam3, "w");

    for(int i = 1; i <= n; i*=2){
        //double x = (b - a) * ((double) i / (double) n);
        fprintf(myfile3, "%d %.15f\n", i, fabs(runge_kutta_fourth(i*2, i*2, a, b, standart_grid, minus_func) - runge_kutta_fourth(i, i, a, b, standart_grid, minus_func)));

    }

    fclose(myfile3);

    return 0;
}