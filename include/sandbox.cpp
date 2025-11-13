#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

int main(){
	auto QuadBowl = [](std::vector<float> in){
		return pow(in[0],2.0) + pow(in[1], 2.0);
	};
	auto DxQuadBowl = [](float x){return 2*x;};
	auto DyQuadBowl = [](float y){return 2*y;};

	//SUFFICIENT DECREASE CONDITION
	//F(xn + an*pn) <= F(xn) + c1 * an * gradTFn * pn
	//CURVATURE CONDITION
	//norm(gradF(xn + an*pn) >= c2 * norm(gradTFn * pn)
	//0 < c1 < c2 < 1
	int iter = 0; 
	float c1 = 10e-4;
	float c2 = 0.9;
	float An = 1.0f;
	float Fn1 = 0;
	
	std::vector<float> Xn = {100,100}; //First guesses;
	std::vector<float> Pn = {-DxQuadBowl(Xn[0]),-DyQuadBowl(Xn[1])};
	//std::function <std::vector<float>(float)> gradientFunc = gradF; 
	std::vector<float> Xn1 = Xn;
	float grad_dot_p = Pn[0]*DxQuadBowl(Xn[0]) + Pn[1]*DyQuadBowl(Xn[1]);
	float curv_grad_dot_p =  Pn[0]*DxQuadBowl(Xn[0] + An*Pn[0]) + Pn[1]*DyQuadBowl(Xn[1] + An*Pn[1]);



	while(iter < 1000){
		for(int i = 0; i < 2; i++){
			Xn1[i] = Xn[i] + An * Pn[i];
		}
		std::vector<float> Curv = {Xn[0] + An*Pn[0], Xn[1] + An*Pn[1]};
	//SUFFICIENT DECREASE CONDITION
	//F(xn + an*pn) <= F(xn) + c1 * an * gradTFn * pn
		if(!(QuadBowl(Xn1) <= QuadBowl(Xn) + c1*An*grad_dot_p)){
			An *= 0.5f;
		}
		//CURVATURE CONDITION
		//norm(gradF(xn + an*pn) >= c2 * norm(gradTFn * pn)
		else if(!(curv_grad_dot_p>= c2*grad_dot_p)){
			An *= 1.5f;
		}
		else{
			Fn1 = QuadBowl(Xn1);
			Xn = Xn1;
			std::cout << "Xn: {" << Xn[0] << "," << Xn[1] << "} Fn1: " << Fn1 << " iter: "<< iter<< std::endl;
		}
	
		if(std::sqrt(pow(DxQuadBowl(Xn[0]),2.0) + pow(DyQuadBowl(Xn[1]),2.0)) < 1e-4){
			std::cout << "SUCCESS"<< std::endl;
			break;
		}
	}
	return 0;
}

