#include <iostream>
#include <cmath>
#include <functional> 
#include <vector>  
#include "Gradient.h"

constexpr int enum_size = 11;
enum class Functions{
    Quad,
    Cubic,
    Power,
    Sin,
    Cos,
    Exp,
    Log,
    Ration,
    Step,
    Gaussian,
    Poly
    
};
enum class Derivatives{
    dQuad,
    dCubic,
    dPower,
    dSin,
    dCos,
    dExp,
    dLog,
    dRation,
    dStep,
    dGaussian,
    dPoly  
};
std::function<float(float)> getFunction(Functions choose){
    int n = 4;
    switch(choose){
        case Functions::Quad:
            return [](float x){return x*x;};
        case Functions::Cubic:
            return [](float x){return pow(x,3);};
        case Functions::Power:
            return [n](float x){return pow(x,n);};
        case Functions::Sin:
            return [](float x){return sin(x);};
        case Functions::Cos:
            return [](float x){return cos(x);};
        case Functions::Exp:
            return [](float x){return exp(x);};
        case Functions::Log:   
            return [](float x){return log(x);};
        case Functions::Ration:
            return [](float x){return (x/(1+pow(x,2)));};
        case Functions::Step:
            return [](float x){return x>=0 ? x: -x;};
        case Functions::Gaussian:   
            return [](float x){return exp(-x*x);};
        case Functions::Poly:   
            return [](float x){return pow(x,3) - 2*x*x + x - 5;};
        default:
            return {}; 
    }
    return NULL; 
}
std::function<float(float)> getDerivative(Derivatives choose){
    int n = 4;
    switch(choose){
        case Derivatives::dQuad:
            return [](float x){return 2*x;};
        case Derivatives::dCubic:
            return [](float x){return 3*pow(x,2);};
        case Derivatives::dPower:
            return [n](float x){return n*pow(x,n-1);};
        case Derivatives::dSin:
            return [](float x){return cos(x);}; 
        case Derivatives::dCos:
            return [](float x){return -sin(x);};
        case Derivatives::dExp:
            return [](float x){return exp(x);}; 
        case Derivatives::dLog:   
            return [](float x){return 1/(x);}; 
        case Derivatives::dRation:
            return [](float x){return (-x*x + 1)/pow((1+x*x),2);};
        case Derivatives::dStep:
            return [](float x){return x>=0 ? 1: -1;};
        case Derivatives::dGaussian:   
            return [](float x){return exp(-x*x)*(-4*x);};
        case Derivatives::dPoly:   
            return [](float x){return 3*pow(x,2) - 4*x + 1;};
        default:
            return {}; 
    }
    return NULL; 
}
struct FunctionConfig{
    Functions func;
    Derivatives derv;
    float x0;
    float step;
    std::string name;
};
std::vector<FunctionConfig> tests = {
       {Functions::Quad, Derivatives::dQuad, 5, 0.5,"Quadatric"},
    {Functions::Cubic, Derivatives::dCubic, 1, 0.3,"Cubic"},
    {Functions::Power, Derivatives::dPower, 2, 0.11,"Power"},
    {Functions::Sin, Derivatives::dSin, 3, 0.1,"Sine"},
    {Functions::Cos, Derivatives::dCos, 3, 0.1,"Cosine"},
    {Functions::Exp, Derivatives::dExp, 0, 0.02,"Exponential"},
    {Functions::Log, Derivatives::dLog, 1, 0.1,"Logarithmic"},
    {Functions::Ration, Derivatives::dRation, 2, 0.2,"Rational"},
    {Functions::Step, Derivatives::dStep, 2, 0.2,"Step"},
    {Functions::Gaussian, Derivatives::dGaussian, 2, 0.1,"Gaussian"},
    {Functions::Poly, Derivatives::dPoly, 4, 0.05,"Polynomial"}
}; 
int main(){
    SimpleGradient<float> Solver;
    Functions F = Functions::Quad; 
    Derivatives D = Derivatives::dQuad;

    Result Res;
    int MaxIter = 100;
    float Error = 0.001f;
    int Guess   = 5;
    float step  = 0.5;
    std::string FuncName = " ";
    std::cout<<"Verbosity Level: "<<endl;
    std::string v = "";
    Log LOG; 
    scanf("%s",v);
    LOG.Verbosity(v);
    for(int i= 0; i < enum_size; i++){
        FuncName = (tests[i].name);
        Res = Solver.Solve(getFunction(tests[i].func),FuncName, getDerivative(tests[i].derv),MaxIter,Error,tests[i].x0,tests[i].step);
        auto ResDeriv = getDerivative(tests[i].derv);
        printf("Result: %s, x: %f, ẋ:%f, Error_Code: %d\n\n",Res.outcome.c_str(), Res.result,ResDeriv(Res.result),Res.E_Code);
    }

    
    
}