#ifndef Gradient
#define Gradient
#include <string.h>
#include <functional>
#include <cfloat>
#include <cmath>


enum Error{
    None,                   //0
    Max_Iterations,         //1
    NaN_Encountered,        //2
    Diverging_Solution,     //3
    Saddle_Point,           //4
    Out_Of_Bounds,          //5
    Growing_Gradient,       //6
    Inflection_Point,       //7
    NaN,                    //8
    Undefined               //9
};
struct Result{
    std::string outcome = "Fail";
    float result;
    Error E_Code;
};
enum VerbLevel{
    INFO  = 0,
    DEBUG = 1,
    ERROR = 2,
    TRACE = 3
};
class Log{
    private:
        bool info,debug,error,trace;
    public:
        void Verbosity(int level){
            info = debug = error = trace = false;
            switch (level)
            {
            case DEBUG:
                debug = true; break;
            case ERROR:
                error = true; break;
            case TRACE:
                trace = true; break;
            default:
                info = true; break;
            }
        }
};
template <typename C> class SimpleGradient{
    public: 
        SimpleGradient(Log& LOG) : logger(LOG){};
        bool D2Check(float e,C xn,float threshold,std::function<C(C)> Derivative){
                const float h = 1e-4;
                auto D2 = [h,Derivative](C X){ return (Derivative(X + h) - Derivative(X - h)) / (2 * h); };

                    if ( abs(D2(xn)) < threshold ){
                        if trace{
                            std::cout<<"D2 = "<<D2(xn)<<"\n";
                        }
                        return false;
                    }
                    else{
                        return true;
                    }
            }
        
    private:
        C Step(C t, C Xn, C Xn1, std::function<C(C)> Deriv, std::function<C(C)> Func){
            float alpha = 1e-4;
            float beta = 0.5;
            auto dfxn = Deriv(Xn);
            auto fxn =  Func(Xn);
            const int max_attempts = 20;
            int attempts = 0;
            while(attempts<max_attempts){
                if(Func(Xn - t * Deriv(Xn)) <= fxn - alpha * t * pow(dfxn,2)){
                    return t;
                }
                attempts++;
                t*=beta;

            }
            return -1;
        }

    Result Solve(std::function<C(C)> Function,std::string FunctionName, std::function<C(C)> Derivative, int max, float e, C initGuess,C step){
        Result x;
        x.outcome = "Fail";
        x.result = initGuess;
        x.E_Code = None;
        float Memo_x = 0;
        int Tolerance = 0.5;
        float threshold = 8*e;
           
        std::cout<<"Function: "<<FunctionName<<"\n";
        for (int i = 0; i < max; ++i){
            const auto fx = Function(x.result);
            const auto dfx = Derivative(x.result);

            if(std::isnan(fx)){
                x.outcome = "Divergent Solution";
                x.E_Code = Diverging_Solution;
                return x;
            }
            else 
            if(abs((x.result)) > FLT_MAX){
                x.outcome = "Out of Bounds";
                x.E_Code  = Out_Of_Bounds;
                return x;
            }
            else{
                C xn = x.result;
                //ITERATE NEW Xn
                x.result = x.result - step*Derivative(x.result); 
                C xn1 = x.result;
                step = Step(step,xn,xn1,Derivative,Function);
                if (step == -1){
                    x.outcome = "Growing Gradient";
                    x.E_Code = Growing_Gradient;
                    return x; 
                }
                C deltaX = e * (1 + abs(xn));
                auto new_dfx = Derivative(x.result);
                auto fxn = Function(xn);
                auto fxn1 = Function(xn1);
                C deltaF = e * (abs(fxn));
                if (abs(dfx)<e){
                    if(abs(xn1 - xn) < deltaX ){
                        if(D2Check(e,xn,threshold,Derivative)){
                        x.outcome = "SUCCESS";
                        x.E_Code = None;
                        return x;
                        }
                    }
                    else{
                        if(abs(fxn1 - fxn) < deltaF){
                            x.outcome = "SUCCESS";
                            x.E_Code = None;
                            return x;
                        }
                        else{
                            //FAIL WHICH FAIL?????
                        }
                    } 
                }      
                else{ 
                    if(abs(xn1 - xn) < deltaX && abs(fxn1 - fxn) < deltaF){
                        //WOULD THIS BE SUCCESS??? 
                            x.outcome = "SUCCESS";
                            x.E_Code = None;
                            return x;
                    }

                } 
            }   
        }
        x.outcome = "Max Iterations";
        x.E_Code = Max_Iterations;
        return x; 
    }	

};
#endif
    