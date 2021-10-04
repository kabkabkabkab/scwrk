#include <stdio.h>
#include <math.h>

double lagr_pol_coeff(int n, int j ,double x){
    double coeff = 1;
    for(int i = 0; i < n; i++){
        if(i != j){
            coeff = coeff*(x - 1 - ((double) i) / ((double) n));
        }
    }

    return coeff;
}

double lagr_pol_appr(double (*f)(double), double (*g)(int, int, double), int n ,double x){

    double pol = 0;

    for(int i = 0; i < n; i++){

        pol +=f( 1 + ((double) i) / ((double) n) ) * g(n, i, x) / g(n, i, ( 1 + ((double) i) / ((double) n) ) ); 

    }

    return pol;

} 

int main(){

    printf("a = %f\n", lagr_pol_appr(log, lagr_pol_coeff, 10, 2) - log(2));
    return 0;    

}