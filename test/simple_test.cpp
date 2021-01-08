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

}
