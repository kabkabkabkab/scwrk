#include <stdio.h>
#include <math.h>

double standart_grid(int i, int n, double a, double b){
    double x = (double)i / (double)n;
    return (b - a)*x + a;
}

double chebyshov_grid( int i, int n, double a, double b ){
    return ( b + a ) / 2.0 + cos( M_PI * (2 * i + 1 ) / ( 2.0 * ( n + 1 ) ) ) * ( b - a ) / 2.0 ;
}

double euler_solver( int i, int n, double a, double b, double (*grid)(int, int, double, double) ){
    // x' = -x, x(0) = 1
    if(i == 0){

        return 1;
    }else{

        //y{i} = y{i-1} + (x{i} - x{i-1})*f(x{i-1}, y{i-1})
        double x2 = grid(i, n, a, b);
        double x1 = grid(i-1, n, a, b);
        double h = x2 - x1;
        //euler_solver(i) = euler_solver(i-1) + (x2 - x1)*f(x1, euler_solver(i-1))
        //f() = - euler_solver(i-1)
        return (1 - h)*euler_solver(i-1, n, a, b, grid);
    }

}

double runge_kutta_second(int i, int n, double a, double b, double (*grid)(int, int, double, double)){
    // x' = -x, x(0) = 1
    if(i == 0){

        return 1;
    }else{
        
        double x2 = grid(i, n, a, b);
        double x1 = grid(i-1, n, a, b);
        // h = x2 - x1
        // a = 3/4
        //y{i} = y{i-1} + h*( 0.25*f(x{i-1}, y{i-1}) + 0.75*f(x{i-1} + 2h/3, y{i-1} + 2h/3 * f(x{i-1}, y{i-1})))
        double y0 = runge_kutta_second(i-1, n, a, b, grid);
        
        double h = (x2 - x1);
        return y0*(1 - h - 1/2 * h * h);
    }
}

double runge_kutta_fourth(int i, int n, double a, double b, double (*grid)(int, int, double, double)){

    if(i == 0){

        return 1;
    }else{
        double x2 = grid(i, n, a, b); // x{i}
        double x1 = grid(i-1, n, a, b); // x{i-1}
        double h = (x2 - x1);
        double y0 = runge_kutta_fourth(i-1, n, a, b, grid);
        //y{i} = y{i} + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        
        //k1 = f(x{i-1}, y{i-1})
        double k1 = - y0;
        
        //k2 = f(x{i-1} + h/2, y{i-1} + h/2 * k1)
        //k2 = - y0 + y0 * h / 2
        double k2 = - (y0 + k1*h / 2.0);
        
        //k3 = f(x{i-1} + h/2, y{i-1} + h/2 * k2)
        //k3 = - y0 + ( y0 - y0 * h / 2) * h / 2
        double k3 = - (y0 + k2*h / 2.0);
        
        //k4 = f(x{i-1} + h, y{i-1} + h * k3)
        //k4 = - y0 - h * ( - y0 + ( y0 - y0 * h / 2) * h / 2 )
        double k4 = - (y0 + h*k3);
        
        // k1 + k4 + 2* (k2 + k3) = 
        // - 6 * y0 + 2 * y0 * h - y0 * h * h + y0 * h * h * h / 4  

        //y0 + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        return ( 1 - h + (2 * h - h * h + h * h * h / 4.0 )*h / 6.0 ) * y0 ;
        
    }

}

int main(){
    int N = 200;
    double a = 0.0;
    double b = 3.0;

    FILE *myfile157471;
    const char filenam1[10] = {"out6.txt"};
    myfile157471 = fopen (filenam1, "w");

    for(int i = 0; i < N; i++){
        double x = (b - a)*( (double)i / (double)N ) + a;
        fprintf(myfile157471, "%f %.15f\n", x, runge_kutta_fourth(i, N, a, b, standart_grid) - runge_kutta_second(i, N, a, b, standart_grid));

    }

    fclose(myfile157471);
    return 0;
}