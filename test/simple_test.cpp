#include "MultidimArrays.h"
#include <iostream>

int main() {

	Multidim::Array<int, 3> a(3, 3, 3);
	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			for (int k = 0; k < 3; k++) {
				a.at(i,j,k) = 1;
			}
		}

	}

	std::cout << a << std::endl;


	Multidim::Array<int, 2> b(3, 9);


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 9; j++) {
			b.at(i,j) = i*i+j;
		}

	}

	std::cout << b << std::endl;

	Multidim::Array<int, 2> window1 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(0,3));

	std::cout << window1 << std::endl;

	Multidim::Array<int, 2> window2 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(3,6));

	std::cout << window2 << std::endl;

	Multidim::Array<int, 2> window3 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(6,9));

	std::cout << window3 << std::endl;

	Multidim::Array<int, 2> window4 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(0,9,3));

	std::cout << window4 << std::endl;

	Multidim::Array<int, 2> window5 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(2,9,3));

	std::cout << window5 << std::endl;

	Multidim::Array<int, 2> window6 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(-1,0,-3));

	std::cout << window6 << std::endl;

	Multidim::Array<int, 2> window7 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(-7,0,-1));

	std::cout << window7 << std::endl;

	for (int i = 0; i < 3; i++) {

		Multidim::Array<int, 1> window = b.subView(Multidim::DimIndex(i), Multidim::DimSlice());

		std::cout << window << std::endl;

	}

	std::cout << std::endl;
}
