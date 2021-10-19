#include <stdio.h>
#include <math.h>
#define H_VAL 0.00001
#define BESSEL_NUM 100

double bessel_function(int m, double x){
	double step = M_PI / BESSEL_NUM;
	double val = 0;
	double arg;
	double integr = 0;
	for(int i = 0; i < BESSEL_NUM; i++){
		double arg = i*step;
		integr += cos( m*arg - x*sin(arg) ) + 4.0*cos( m*(arg + 0.5*step) - x*sin(arg + 0.5*step) ) + cos( m*(arg + step) - x*sin(arg + step) );

	}

	return integr * step / (6.0 * M_PI);
}

double derivation_func(double (*f)(int, double), int m, double arg){	
	
	if((arg + H_VAL < 2*M_PI) || (arg - H_VAL > 0)){
		
		return (f(m, arg + H_VAL) - f(m, arg - H_VAL)) / (2.0*H_VAL);

	}

	else{

		if(arg + H_VAL > 2*M_PI){

			return (f(m, arg) - f(m, arg - H_VAL)) / H_VAL;

		}

		else{

			printf("error > 2pi!\n");
			return 0;

		}

		if(arg - H_VAL < 0){
			return (f(m, arg + H_VAL) - f(m, arg)) / H_VAL;
		}

		else{

			printf("error < 0!\n");
			return 0;

		}

	}	
	
}

int main(){
	double arg = 1.0;
	int n = 10000;
	printf("\ndJ0(%f)/dx = %.15f\n", arg, derivation_func(bessel_function, 0, arg) );	
	printf("J1(%f) = %.15f\n", arg, bessel_function(1, arg) );
	printf("dJ0(%f)/dx + J1(%f) = %.15f\n", arg, arg, ( derivation_func(bessel_function, 0, arg) + bessel_function(1, arg) ) );

	FILE *myfile1212;
    const char filenam1212[17] = {"file_of_out.txt"};
    myfile1212 = fopen (filenam1212, "w");

    for(int i = 0; i < n; i++){
		arg = M_PI*((double)i) / ((double)n);
        fprintf(myfile1212, "%.15f %.15f\n", arg, ( bessel_function(0, arg + 1 / ((double)n)) - bessel_function(0, arg + 1 / ((double)n)) / (2*(double)n) + bessel_function(1, arg) ) ); 
		//  bessel_function(1, arg)*1.5 - 0.5*bessel_function(1, -arg)
		//derivation_func(bessel_function, 0, arg) + bessel_function(1, arg)
    }

    fclose(myfile1212);


	return 0;
}
