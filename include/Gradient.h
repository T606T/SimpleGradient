#ifndef Gradient
#define Gradient
#include <string.h>
#include <functional>
#include <cfloat>
#include <cmath>
#include <fstream>


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
    Undefined,              //9
    Stalled                 //10
};
struct Result{
    std::string outcome = "Fail";
    float result;
    Error E_Code;
};
enum VerbLevel{
    DEBUG = 0,
    ERROR = 1,
    TRACE = 2,
    INFO  = 3
};
class Log{
    private:
        bool Information,Debugging,Errors,Tracing;
    public:
        void Verbosity(int level){
            Information = Debugging = Errors = Tracing = false;
            switch (level)
            {
            case DEBUG:
                Debugging = true; break;
            case ERROR:
                Errors = true; break;
            case TRACE:
                Tracing = true; break;
            default:
                Information = true; break;
            }
        }
        //methods to return verbosity status.
        bool debug() const {return Debugging;};
        bool error() const {return Errors;};
        bool trace() const {return Tracing;};
        bool info()  const {return Information;};

};
template <typename C> class SimpleGradient{
    private:
        Log* logger = nullptr; 
        C Step(C t, C Xn, std::function<C(C)> Deriv, std::function<C(C)> Func){
            float alpha = 1e-4;
            float beta = 0.5;
            auto dfxn = Deriv(Xn);
            auto fxn =  Func(Xn);
            const int max_attempts = 20;
            int attempts = 0;
            while(attempts<max_attempts){
                if(Func(Xn - t * Deriv(Xn)) <= fxn - alpha * t * pow(dfxn,2)){
                    if(logger && logger->debug() && DebugFile.is_open()){
                        DebugFile << "[Armijo] t: "<< t << " iteration: " << attempts << "\n";
                    }
                    return t;
                }
                attempts++;
                t*=beta;

            }
            return -1;
        }
    public:
        std::fstream DebugFile;
	SimpleGradient() = default;
        SimpleGradient(Log& logger) : logger(&logger){
            if (logger.debug()){
                DebugFile.open("DebugFile.txt", std::ios::out | std::ios::trunc);
                if (!DebugFile.is_open()){
                    throw std::runtime_error("Failed to open File");
                }
            } 
        };
        ~SimpleGradient(){
            if (logger && logger->debug() && DebugFile.is_open()){
                DebugFile.close();
            }
        };
        bool D2Check(C xn,float threshold,std::function<C(C)> Derivative){
                const float h = 1e-4;
                auto D2 = [h,Derivative](C X){ return (Derivative(X + h) - Derivative(X - h)) / (2 * h); };

                    if ( abs(D2(xn)) < threshold ){
                        if (logger && logger->trace()){
                            std::cout<<"D2 = "<<D2(xn)<<"\n";
                        }
                        return false;
                    }
                    else{
                        return true;
                    }
            }

        Result Solve(std::function<C(C)> Function,std::string FunctionName, std::function<C(C)> Derivative, int max, float e, C initGuess,C step){
            Result x;
            x.outcome = "Fail";
            x.result = initGuess;
            x.E_Code = None;
            auto fx = Function(x.result);         //F(Xn)
            auto dfx = Derivative(x.result);      //dF(Xn)
            const float Tolerance = 1e-5f;        //Tolerance for dF(Xn) when it reaches enough closeness to 0. 
            C xn;
            C xn1;
            step = std::min(step, 1.0f / (1.0f + std::fabs(dfx)));
            float threshold = 8*e;
            if(logger && logger->debug() && DebugFile.is_open()){
                DebugFile << "\n---------------------------------------------------------------\n";
                DebugFile << FunctionName<<"\n";
            }
            
/////////////////////////////////////////////////////////////////////MAIN FOR LOOP ////////////////////////////////////////////////////////////////////////////
            for (int i = 0; i < max; ++i){
                fx = Function(x.result);
                dfx = Derivative(x.result);
                //Iteration Counter
                if(logger && logger->debug() && DebugFile.is_open()){
                    DebugFile << "[Algorythm iteration]: "<< i<<"\n";
                }
                //Checking if fx Is not-a-number value
                if(std::isnan(fx) || !std::isfinite(fx)){
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
                    xn = x.result;
                    //Check the learning rate to fulfill Armijo condition.
                    float memo_step = step;
                    step = Step(step,xn,Derivative,Function);
                    //Should Armijo Fail do a rollback.
                    if(step < 0){
                        step = memo_step * 0.5;
                    }
////////////////////ITERATE NEW Xn///////////////////////////////////////////////////////////
                    x.result = x.result - step*Derivative(x.result); 
                    xn1 = x.result;
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
                    if (logger && logger->debug() && DebugFile.is_open()){
                        DebugFile << "[f(xn)]:"<< fxn << "\n";
                        DebugFile << "[f(xn+1)]:"<< fxn << "\n";
                        DebugFile << "[new_dfx]:"<< new_dfx << "\n";
                    }

                    if (abs(dfx)<e){
                        if(abs(xn1 - xn) < deltaX){
                            if(D2Check(xn,threshold,Derivative)){
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
                                if (logger && logger->trace()){
                                    x.outcome = "Undefined error"; //To later determine.
                                    x.E_Code = Undefined;
                                    return x;
                                }

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
                
                step *= 1.05f;
                if(fabs(dfx) < Tolerance){
                    x.outcome = "SUCCESS";
                    x.E_Code = None;
                    return x;
                }
                //TO fix because the algorythm sometimes confuses this outcome with converged late.
                if (std::fabs(dfx) < 1e-3 && std::fabs(xn1 - xn) < e) {
                    x.outcome = "Flat Region (stalled)";
                    x.E_Code = Stalled;
                    return x;
                }  
            }
            if (std::fabs(Derivative(x.result)) < e * 10) {
                x.outcome = "Converged Late";
                x.E_Code = None;
                return x;
            }

            x.outcome = "Max Iterations";
            x.E_Code = Max_Iterations;
            return x; 
        }	

    };
#endif
    
