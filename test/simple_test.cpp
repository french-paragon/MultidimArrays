#include "MultidimArrays.h"
#include "MultidimIndexManipulators.h"
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
	using MarrayB = typeof (b);

	static_assert (!MarrayB::isConstView(), "b is not non const view");

	using MarrayBSub = typeof (b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(0,3)));

	static_assert (!MarrayBSub::isConstView(), "sub is not non const view");

	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 9; j++) {
			b.at(i,j) = i*i+j;
		}

	}

	std::cout << b << std::endl;

	Multidim::Array<int, 2> window1 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(0,3));

	std::cout << window1 << ' ' << window1.isDense() << std::endl;

	Multidim::Array<int, 2> window2 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(3,6));

	std::cout << window2 << ' ' << window2.isDense() << std::endl;

	Multidim::Array<int, 2> window3 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(6,9));

	std::cout << window3 << ' ' << window3.isDense() << std::endl;

	Multidim::Array<int, 2> window4 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(0,9,3));

	std::cout << window4 << ' ' << window4.isDense() << std::endl;

	Multidim::Array<int, 2> window5 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(2,9,3));

	std::cout << window5 << ' ' << window5.isDense() << std::endl;

	Multidim::Array<int, 2> window6 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(-1,0,-3));

	std::cout << window6 << ' ' << window6.isDense() << std::endl;

	Multidim::Array<int, 2> window7 = b.subView(Multidim::DimSlice(0,3), Multidim::DimSlice(-7,0,-1));

	std::cout << window7 << ' ' << window7.isDense() << std::endl;

	for (int i = 0; i < 3; i++) {

		Multidim::Array<int, 1> window = b.subView(Multidim::DimIndex(i), Multidim::DimSlice());

		std::cout << window << ' ' << window.isDense() << std::endl;

	}

	std::cout << std::endl;

	constexpr int testIdxConvDims = 3;
	constexpr int testIdxConvDimsSize = 3;

	std::array<int, testIdxConvDims> shape;

	for (int i = 0; i < testIdxConvDims; i++) {
		shape[i] = testIdxConvDimsSize;
	}

	Multidim::IndexConverter<3> testConverter(shape, Multidim::DimsExclusionSet<testIdxConvDims>(testIdxConvDims-1));


	std::cout << std::endl;

	std::cout << testConverter.numberOfPossibleIndices() << std::endl;

	for (int i = 0; i < testConverter.numberOfPossibleIndices(); i++) {

		auto idx = testConverter.getIndexFromPseudoFlatId(i);

		std::cout << "index " << i << ": \t[ ";
		for (int j = 0; j < testIdxConvDims; j++) {
			std::cout << idx[j] << ' ';
		}
		std::cout << "]" << std::endl;
	}

	std::cout << std::endl;
}
