#include <iostream>
#include <cmath>
#include <string>
#include <cfloat>
#include "Gradient.h"

float F(float x){
	return exp(x);
}
float dF(float x){
    return exp(x);
}

int main(){
    FunctionPair<float> Func;
    SimpleGradient<float> Solver;
    Func.deriv = dF;
    Func.funct = F;
    Result Res = Solver.Solve(Func, 200,0.001f,2,0.005);
    printf("\nResult: %s, X: %f, Error_Code: %d\n",Res.outcome.c_str(), Res.result,Res.E_Code);
    return 0;
}