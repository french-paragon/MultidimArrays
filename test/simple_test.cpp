#include "MultidimArrays.h"
#include <iostream>

int main() {

	Multidim::Array<int, 3> a(3, 3, 3);

	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			for (int k = 0; k < 3; k++) {
				a.at(i,j,k) = 1;
				std::cout << a.value(i,j,k) << " ";
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;

	}

	std::cout << std::endl;


	Multidim::Array<int, 2> b(3, 9);


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 9; j++) {
			b.at(i,j) = i*i+j;
			std::cout << b.value(i,j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;

	Multidim::Array<int, 2> window1 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(0,3));


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {
			std::cout << window1.value(i,j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;

	Multidim::Array<int, 2> window2 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(3,6));


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {
			std::cout << window2.value(i,j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;

	Multidim::Array<int, 2> window3 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(6,9));


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {
			std::cout << window3.value(i,j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;

	Multidim::Array<int, 2> window4 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(0,9,3));


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {
			std::cout << window4.value(i,j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;

	Multidim::Array<int, 2> window5 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(2,9,3));


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {
			std::cout << window5.value(i,j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;

	Multidim::Array<int, 2> window6 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(-1,0,-3));


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {
			std::cout << window6.value(i,j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;

	Multidim::Array<int, 2> window7 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(-7,0,-1));


	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {
			std::cout << window7.value(i,j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;

	for (int i = 0; i < 3; i++) {

		Multidim::Array<int, 1> window = b.subView(Multidim::DimIndex(i), Multidim::DimSlice());

		for (int j = 0; j < 9; j++) {
			std::cout << window.value(j) << " ";
		}
		std::cout << std::endl;

	}

	std::cout << std::endl;
}
