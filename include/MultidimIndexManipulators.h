#ifndef MULTIDIMINDEXCONVERTERS_H
#define MULTIDIMINDEXCONVERTERS_H

/*
Copyright 2022 Paragon<french.paragon@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <array>
#include <type_traits>

namespace Multidim {

template<int nDim>
/*!
 * \brief The DimsExclusionSet class is used to mark some dimensions of a multi-dimensional volume as being excluded from the others.
 */
class DimsExclusionSet {

public:
    template<typename... Ds, std::enable_if_t<std::conjunction_v<std::is_integral<Ds>...>>* = nullptr>
    DimsExclusionSet(Ds ... maskedDims) {

		static_assert(sizeof...(maskedDims) <= nDim,
				"Number of dimension provided is bigger than the number of dimensions in the set !");

		for (int i = 0; i < nDim; i++) {
			_mask[i] = false;
		}

		for (int masked : std::array<int, sizeof...(maskedDims)>({maskedDims...})) {
			_mask[masked] = true;
		}
	}
    template<std::size_t nElems>
	DimsExclusionSet(std::array<int, nElems> const& maskedDims) {
		static_assert(nElems <= nDim,
				"Number of dimension provided is bigger than the number of dimensions in the set !");

		for (int i = 0; i < nDim; i++) {
			_mask[i] = false;
		}

		for (int masked : maskedDims) {
			_mask[masked] = true;
		}
	}

	DimsExclusionSet(DimsExclusionSet<nDim> const& other) :
		_mask(other._mask)
	{

	}

	DimsExclusionSet<nDim> inverse() const {
		IndexMaskBlock mask;

		for (int i = 0; i < nDim; i++) {
			mask[i] = !_mask[i];
		}

		DimsExclusionSet<nDim> r;
		r._mask = mask;

		return r;
	}

	bool indexIsExcluded(int index) const {
		return _mask[index];
	}

protected:
	typedef std::array<bool, nDim> IndexMaskBlock;

	IndexMaskBlock _mask;
};

template<int nDim>
/*!
 * \brief The MultidimIndexConverter class is used to convert between a pseudo-flat index and multidimensional indices for a mutlidimensional volumetric shapes.
 *
 * The class is very low level, and will blindly give answers even if the provided indices are invalid.
 * You are thus at risk of segmentation fault if you blindly use the indices this class generate for you without checking they are valid first.
 *
 * The flat indices provided might not represent the underlying data storage order used by a given tensor or any other class.
 *
 * Some dimensions can be excluded to build a linear iterator over only a subset of dimensions.
 */
class IndexConverter {

public:
	typedef std::array<int, nDim> ShapeBlock;

	IndexConverter(ShapeBlock const& shape,
						   DimsExclusionSet<nDim> const& excludedDims = DimsExclusionSet<nDim>()) :
		_excludedDims(excludedDims),
		_shape(shape)
	{
		int ps = 1;

		for (int i = 0; i < nDim; i++) {
			if (!_excludedDims.indexIsExcluded(i)) {
				_pseudoStride[i] = ps;
				ps *= _shape[i];
			} else {
				_pseudoStride[i] = -1;
			}
		}
	}

	static constexpr ShapeBlock initialIndex() {
		ShapeBlock initial {};

		for (int i = 0; i < nDim; i++) {
			initial[i] = 0;
		}
		return initial;
	}

	static constexpr ShapeBlock invalidIndex() {
		ShapeBlock invalid {};

		for (int i = 0; i < nDim; i++) {
			invalid[i] = -1;
		}
		return invalid;
	}

	int numberOfPossibleIndices() const {

		if (nDim == 0) {
			return 0;
		}

		int n = 1;

		for (int i = 0; i < nDim; i++) {
			if (!_excludedDims.indexIsExcluded(i)) {
				n *= _shape[i];
			}
		}

		return n;
	}

	ShapeBlock getIndexFromPseudoFlatId(int index, ShapeBlock const& valForExcludedIndex = initialIndex()) const {

		ShapeBlock out {};
		int leftOver = index;

		for (int i = nDim-1; i >= 0; i--) {

			bool isExcluded = _excludedDims.indexIsExcluded(i);

			//garanteed optimizable to branchless !
			out[i] = (isExcluded) ? valForExcludedIndex[i] : leftOver/_pseudoStride[i];
			leftOver %= (isExcluded) ? index+1 : _pseudoStride[i];
		}

		return out;
	}

	int getPseudoFlatIdFromIndex(ShapeBlock const& idx) {

		if (nDim == 0) {
			return 0;
		}

		int flat = 0;

		for (int i = 0; i < nDim; i++) {
			if (!_excludedDims.indexIsExcluded(i)) {
				flat += _pseudoStride[i]*idx[i];
			}
		}

		return flat;
	}

protected:
	DimsExclusionSet<nDim> _excludedDims;
	ShapeBlock _shape;
	ShapeBlock _pseudoStride;

};

}

#endif // MULTIDIMINDEXCONVERTERS_H
