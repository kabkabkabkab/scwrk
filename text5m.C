#include <stdio.h>
#include <string.h>
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

double func_of(double x){
    return 1 / (1 + 25*x*x); 

}

int main(){

    int n = 100;

    FILE *myfile1111;
    const char filenam1111[18] = {"file_of_data.txt"};
    myfile1111 = fopen (filenam1111, "w");

    for(int i = 5; i < 20; i++){
        for(int j = 0; j < 2*n; j++){
            fprintf(myfile1111, "%d %.15f %.15f\n", i, (-1 + ((double)j) / ((double)n)), lagr_pol_appr(func_of, lagr_pol_coeff, i, -1 + ((double)j) / ((double)n) ) - func_of(-1 + ((double)j) / ((double)n)));
        }

        fprintf(myfile1111, "\n", i);
    }

    fclose(myfile1111);

    return 0;    

}