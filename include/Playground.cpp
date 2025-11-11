#include <iostream>
#include <cmath>
#include <vector>

int main () {
  auto QuadBowl = [](std::vector<float> in){
    return pow(in[0],2.0) + pow(in[1],2.0);
  };
  auto DxQuadBowl = [](float x){return 2*x;};
  auto DyQuadBowl = [](float y){return 2*y;};

  //SUFFICIENT DECREASE CONDITION
  //F(xn + an*pn) <= F(xn) + c1 *an*gradTFn*pn
  //CURVATURE CONDITION 
  //norm(gradF(xn + an*pn)T * pn) >= c2 * norm(gradTFn * pn)
  // 0 < c1 < c2 < 1 
  int iter = 0;
  float c1 = 10e-4;
  float c2 = 0.9;
  std::vector<float> An = {0.1, 0.1};
  std::vector<float> Pn;
  std::vector<float> Xn = {5, 5};
    while( iter < 1000){
      std::vector<float> gradF = {DxQuadBowl(Xn[0]), DyQuadBowl(Xn[1])};
      std::vector<float> gradFT = {{DxQuadBowl(Xn[0])},{DyQuadBowl(Xn[1])}};

      for (int i = 0; i < Xn.size(); i++){
        std::vector<float> var = Xn[i] + An[i]* gradF[i]; 
      }
      if(QuadBowl(var) <= QuadBowl(Xn) + c1*An*gradFT*gradF){
        Xn1 = QuadBowl(var);
      }


    }
  return 0;
}