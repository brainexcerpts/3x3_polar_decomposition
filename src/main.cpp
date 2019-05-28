#include <iostream>

#include "polar_decomposition_3x3.hpp"
#include "svd_3x3.hpp"
#include "eigen_3x3.hpp"

void print(const Mat3& m){
    std::cout << m.a << ", " << m.b << ", " << m.c << std::endl;
    std::cout << m.d << ", " << m.e << ", " << m.f << std::endl;
    std::cout << m.g << ", " << m.h << ", " << m.i << std::endl;
}

int main() {

    Mat3 rot   = Mat3::rotate(Vec3(1.0f,1.0f,1.0f).normalized(), -1.1f);
    Mat3 scale = Mat3::diagonal(5.0f);
    Mat3 m = scale * rot;

    std::cout << "INPUT\n";
    std::cout << "scale: \n";
    print(scale);
    std::cout << "rotation: \n";
    print(rot);
    std::cout << "m: \n";
    print(m);

    Polar_decomposition<false> decomp(m);
    std::cout << "Decompose\n";
    std::cout << "scale: \n";
    print(decomp.matrix_S());
    std::cout << "rotation: \n";
    print(decomp.matrix_R());
    std::cout << "combine: \n";
    print(decomp.matrix_S()*decomp.matrix_R());

    return 0;
}
