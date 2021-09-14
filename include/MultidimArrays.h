#ifndef MULTIDIM_ARRAY_H
#define MULTIDIM_ARRAY_H

/*
Copyright 2021 Paragon<french.paragon@gmail.com>

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
#include <stdexcept>
#include <cstring>
#include <tuple>
#include <assert.h>

namespace Multidim {

typedef int array_size_t;

template<int d>
constexpr bool allDimsValid() {
	return d > 0;
}

template<int d, int... ds>
constexpr bool allDimsValid() {
	return allDimsValid<d>() and allDimsValid<ds...>();
}

template<int... ds>
constexpr bool staticDimsCheck() {
	return allDimsValid<ds...>();
}

template<>
constexpr bool staticDimsCheck() {
	return true;
}

template<typename ... Ds>
bool allDimsValid(array_size_t d, Ds... dims) {
	return allDimsValid(d) and allDimsValid(dims...);
}

template<>
inline bool allDimsValid(array_size_t d) {
	return d > 0;
}

template<size_t S>
bool allDimsValid(std::array<array_size_t, S> const& array) {
	for (size_t i = 0; i < S; i++) {
		if (!allDimsValid(array[i])) {
			return false;
		}
	}
	return true;
}

enum class AccessCheck {
	Check = 1,
	Nocheck = 0
};

struct DimInfo {
	virtual bool stayDim() const = 0;
};

struct DimSlice : public DimInfo {

	DimSlice(array_size_t sIndex = 0,
			 array_size_t eIndex = -1,
			 array_size_t jump = 1) :
		startIndex(sIndex),
		endIndex(eIndex),
		indexJump(jump)
	{

	}

	array_size_t startIndex;
	array_size_t endIndex;
	array_size_t indexJump;

	virtual bool stayDim() const override {
		return true;
	}
};

struct DimIndex : public DimInfo {

	DimIndex(array_size_t id) :
		index(id)
	{

	}

	array_size_t index;

	virtual bool stayDim() const override {
		return false;
	}
};


template<typename s, typename... slices>
constexpr array_size_t slicedDims() {
	return slicedDims<s>() + slicedDims<slices...>();
}

template<>
constexpr array_size_t slicedDims<DimSlice>() {
	return 1;
}

template<>
constexpr array_size_t slicedDims<DimIndex>() {
	return 0;
}

template<typename T, array_size_t nDim>
class Array {

	static_assert(nDim > 0,
                  "No dimensions provided");

public:

	typedef std::array<array_size_t, nDim> ShapeBlock;

	class IndexBlock : public ShapeBlock {
	public:
		IndexBlock& operator+=(array_size_t scalar) {
			for (int i = 0; i < nDim; i++) {
				(*this)[i] += scalar;
			}
			return *this;
		}

		IndexBlock& operator+=(IndexBlock const& other) {
			for (int i = 0; i < nDim; i++) {
				(*this)[i] += other[i];
			}
			return *this;
		}

		IndexBlock operator+(array_size_t scalar) {
			IndexBlock b = *this;
			b += scalar;
			return b;
		}

		IndexBlock operator+(IndexBlock const& other) {
			IndexBlock b = *this;
			b += other;
			return b;
		}

		IndexBlock& operator-=(array_size_t scalar) {
			for (int i = 0; i < nDim; i++) {
				(*this)[i] -= scalar;
			}
			return *this;
		}

		IndexBlock& operator-=(IndexBlock const& other) {
			for (int i = 0; i < nDim; i++) {
				(*this)[i] -= other[i];
			}
			return *this;
		}

		IndexBlock operator-(array_size_t scalar) {
			IndexBlock b = *this;
			b -= scalar;
			return b;
		}

		IndexBlock operator-(IndexBlock const& other) {
			IndexBlock b = *this;
			b -= other;
			return b;
		}

		void clip(ShapeBlock const& constraintShape) {
			for (int j = 0; j < nDim; j++) {
				if ((*this)[j] < 0) (*this)[j] = 0;
				if ((*this)[j] >= constraintShape[j]) (*this)[j] = constraintShape[j] - 1;
			}
		}

		void moveToNextIndex(ShapeBlock const& constraintShape) {
			(*this)[0] += 1;
			for (int j = 0; j < nDim-1 and (*this)[j] >= constraintShape[j]; j++) {
				(*this)[j] = 0;
				(*this)[j+1] += 1;
			}
		}

		bool isInLimit(ShapeBlock const& constraintShape) {
			for (int i = 0; i < nDim; i++) {
				if ((*this)[i] < 0 or (*this)[i] >= constraintShape[i]) {
					return false;
				}
			}
			return true;
		}
	};

	template<array_size_t... dims>
	Array() :
		_shape({dims...})
	{
		static_assert(sizeof...(dims) == nDim or sizeof...(dims) == 0,
				"Wrong number of dimensions provided");

		static_assert(staticDimsCheck<dims...>(),
			"All dimensions must be greather or equal zero.");

		if (sizeof...(dims) == 0) {
			for (int  i = 0; i < nDim; i++) {
				_shape[i] = 0;
			}
		}

		int s = 1;
		for (int i = 0; i < nDim; i++) {
			_strides[i] = s;
			s *= _shape[i];
		}

		if (s > 0) {
			_manageData = true;
			_data = new T[s];
		} else {
			_manageData = false;
			_data = nullptr;
		}
	}


	template<typename... Ds>
	Array(Ds... dims) :
		_shape({dims...})
	{
		static_assert(sizeof...(dims) == nDim or sizeof...(dims) == 0,
				"Wrong number of dimensions provided");

		if (sizeof...(dims) == 0) {
			for (int  i = 0; i < nDim; i++) {
				_shape[i] = 0;
			}
		} else {
			assert(allDimsValid(dims...));
		}

		int s = 1;
		for (int i = 0; i < nDim; i++) {
			_strides[i] = s;
			s *= _shape[i];
		}

		if (s > 0) {
			_manageData = true;
			_data = new T[s];
		} else {
			_manageData = false;
			_data = nullptr;
		}
	}

	Array(ShapeBlock const& shape) :
		_shape(shape)
	{

		assert(allDimsValid(shape));

		int s = 1;
		for (int i = 0; i < nDim; i++) {
			_strides[i] = s;
			s *= _shape[i];
		}

		if (s > 0) {
			_manageData = true;
			_data = new T[s];
		} else {
			_manageData = false;
			_data = nullptr;
		}
	}

	Array(ShapeBlock const& shape, ShapeBlock const& strds) :
		_shape(shape),
		_strides(strds)
	{

		assert(allDimsValid(shape));

		ShapeBlock idx = _shape;
		for (int i = 0; i < nDim; i++) {
			idx[i] = idx[i]-1;
		}
		int s = flatIndex(idx) + 1;
		int v = 1;
		for (int i = 0; i < nDim; i++) {
			v *= _shape[i];
		}

		assert(s == v);

		if (s > 0) {
			_manageData = true;
			_data = new T[s];
		} else {
			_manageData = false;
			_data = nullptr;
		}
	}

	Array(T* data, ShapeBlock const& shape, ShapeBlock const& strides) :
		_shape(shape),
		_strides(strides),
		_data(data),
		_manageData(false)
	{

	}

	Array(Array<T, nDim>& other) :
		_shape(other._shape),
		_strides(other._strides)
	{
		if (other._manageData) {

			int s = flatLenght();

			if (s > 0) {
				_manageData = true;
				_data = new T[s];

				memcpy(_data, other._data, sizeof (T)*s);
			} else {
				_manageData = false;
				_data = nullptr;
			}
		} else {
			_manageData = false;
			_data = other._data;
		}
	}

	Array(Array<T, nDim> const& other) :
		_shape(other._shape),
		_strides(other._strides)
	{
		int s = flatLenght();

		if (s > 0) {
			_manageData = true;
			_data = new T[s];

			memcpy(_data, other._data, sizeof (T)*s);
		} else {
			_manageData = false;
			_data = nullptr;
		}
	}

	Array(Array<T, nDim>&& other) :
		_shape(other._shape),
		_strides(other._strides),
		_manageData(other._manageData),
		_data(other._data)
	{
		other._data = nullptr;
	}

	~Array() {
		if (_manageData and _data != nullptr) {
			delete [] _data;
		}
	}

	ShapeBlock shape() const {
		return _shape;
	}

	int flatLenght() const {

		int s = 1;
		for (int i = 0; i < nDim; i++) {
			s *= _shape[i];
		}

		return s;
	}

	ShapeBlock strides() const {
		return _strides;
	}

	bool isView() const {
		return _manageData;
	}

	bool empty() const {
		for (int i = 0; i < nDim; i++) {
			if (_shape[i] == 0) {
				return true;
			}
		}
		return false;
	}

private:

	array_size_t sliceFirstIndex(DimSlice const& slice) {
		return slice.startIndex;
	}

	array_size_t sliceFirstIndex(DimIndex const& index) {
		return index.index;
	}

	template<typename... slices>
	ShapeBlock firstIndex(slices... s) {
		ShapeBlock blk = {sliceFirstIndex(s)...};
		for (int i = 0; i < nDim; i++) {
			if (blk[i] < 0) {
				blk[i] = _shape[i] + blk[i];
			}
		}
		return blk;
	}

	template<typename... slices>
	int dataOffset(slices... s) {
		return flatIndex(firstIndex(s...));
	}

	template<typename... slices>
	typename Array<T, slicedDims<slices...>()>::ShapeBlock getSubShape(slices... s) {

		using SBlock = typename Array<T, slicedDims<slices...>()>::ShapeBlock;

		std::array<const DimInfo *, sizeof... (slices)> a = {static_cast<const DimInfo *>(&s)...};
		SBlock r;

		int p = 0;
		for (int i = 0; i < sizeof... (slices); i++) {
			if (a[i]->stayDim()) {
				const DimSlice * s = static_cast<const DimSlice *>(a[i]);

				array_size_t eId = (s->endIndex >= 0) ? s->endIndex : (_shape[i] + s->endIndex);
				array_size_t sId = (s->startIndex >= 0) ? s->startIndex : (_shape[i] + s->startIndex);
				array_size_t step = s->indexJump;

				if (step == 0) {
					throw std::out_of_range("Step 0 is invalid");
				}

				if (s->endIndex >= 0 and s->indexJump > 0) {
					eId = eId-1;
				}

				if (sId < 0 or sId >= _shape[i]) {
					throw std::out_of_range("Start index out of range");
				}

				if (eId < 0 or eId >= _shape[i]) {
					throw std::out_of_range("End index out of range");
				}

				if (sId > eId) {

					if (step >= 0) {
						throw std::out_of_range("Step is invalid");
					}

					array_size_t tmp = eId;
					eId = sId;
					sId = tmp;

					step = -step;

				} else if (step <= 0) {
					throw std::out_of_range("Step is invalid");
				}

				array_size_t range = std::abs(eId - sId) + 1;

				r[p] = 1 + (range-1)/step;
				p++;
			}
		}

		return r;
	}

	template<typename... slices>
	typename Array<T, slicedDims<slices...>()>::ShapeBlock getSubStrides(slices... s) {

		using SBlock = typename Array<T, slicedDims<slices...>()>::ShapeBlock;

		std::array<DimInfo const*, sizeof... (slices)> a = {static_cast<DimInfo const*>(&s)...};
		SBlock r;

		int p = 0;
		for (int i = 0; i < sizeof... (slices); i++) {
			if (a[i]->stayDim()) {

				const DimSlice * s = static_cast<const DimSlice *>(a[i]);

				r[p] = _strides[i]*s->indexJump;
				p++;
			}
		}

		return r;
	}

public:

	template<typename... slices>
	Array<T, slicedDims<slices...>()> subView(slices... s) {

		static_assert (sizeof... (slices) == nDim, "Cannot generate an array with zero dimensions");
		static_assert (slicedDims<slices...>() > 0, "Cannot generate an array with zero dimensions");

		using SubArray = Array<T, slicedDims<slices...>()>;
		using SubShapeBlock = typename SubArray::ShapeBlock;

		int offset = dataOffset(s...);
		SubShapeBlock shape = getSubShape(s...);
		SubShapeBlock strides = getSubStrides(s...);

		return SubArray(_data + offset, shape, strides);
	}

	Array<T, std::max(1,nDim-1)> sliceView(int dim, int coord) {

		if (dim < 0 or dim >= nDim) {
			throw std::out_of_range("Dim index out of range");
		}

		if (coord < 0 or coord >= shape()[dim]) {
			throw std::out_of_range("Index out of range");
		}

		using SubArray = Array<T, std::max(1,nDim-1)>;
		using SubShapeBlock = typename SubArray::ShapeBlock;

		int offset = strides()[dim]*coord;
		SubShapeBlock subshape;
		SubShapeBlock substrides;

		if (nDim == 1) {
			subshape = {1};
			substrides = {strides()[0]};

			return SubArray(_data + offset, subshape, substrides);
		}

		for (int i = 0; i < dim; i++) {
			subshape[i] = shape()[i];
			substrides[i] = strides()[i];
		}

		for (int i = dim+1; i < nDim; i++) {
			subshape[i-1] = shape()[i];
			substrides[i-1] = strides()[i];
		}

		return SubArray(_data + offset, subshape, substrides);
	}


	template<AccessCheck c = AccessCheck::Check>
	int flatIndex(ShapeBlock const& idxs) const {

		if (c == AccessCheck::Check) {
			for (int i = 0; i < nDim; i++) {
				if (idxs[i] < 0 or idxs[i] >= _shape[i]) {
					throw std::out_of_range("Index out of range");
				}
			}
		}

		int fIds = 0;

		for (int i = 0; i < nDim; i++) {
			fIds += idxs[i]*_strides[i];
		}

		return fIds;

	}

	template<AccessCheck c = AccessCheck::Check, typename... Ds>
	int flatIndex(int idx0, Ds... idxs) const {

		static_assert(sizeof...(idxs) == nDim -1 or sizeof...(idxs) == 0,
				"Wrong number of indices provided");

		if (sizeof...(idxs) == 0) {

			if (c == AccessCheck::Check) {
				if (idx0 < 0 or idx0 >= flatLenght()) {
					throw std::out_of_range("Index out of range");
				}
			}

			return idx0*_strides[0];

		} else {

			ShapeBlock idxs_array = {idx0, idxs...};
			return flatIndex<c>(idxs_array);

		}

	}

	bool copyData(Array<T,nDim> const& other) {

		if (other.shape() != shape()) {
			return false;
		}

		int n = flatLenght();

		const T* ptr = other._data;
		if (ptr >= _data and ptr < _data + n) {
			return false;
		}

		memcpy(_data, ptr, n*sizeof (T));

		return true;

	}

	template<AccessCheck c = AccessCheck::Check, typename... Ds>
	T& at(Ds... idxs) {
		int fIds = flatIndex<c>(idxs...);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check>
	T& at(ShapeBlock const& idxs) {
		int fIds = flatIndex<c>(idxs);
		return _data[fIds];
	}

	template<typename... Ds>
	T& atUnchecked(Ds... idxs) {
		int fIds = flatIndex<AccessCheck::Nocheck>(idxs...);
		return _data[fIds];
	}

	T& atUnchecked(ShapeBlock const& idxs) {
		int fIds = flatIndex<AccessCheck::Nocheck>(idxs);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check, typename... Ds>
	T value(Ds... idxs) const {
		int fIds = flatIndex<c>(idxs...);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check>
	T value(ShapeBlock const& idxs) const {
		int fIds = flatIndex<c>(idxs);
		return _data[fIds];
	}

	template<typename... Ds>
	T valueUnchecked(Ds... idxs) const {
		int fIds = flatIndex<AccessCheck::Nocheck>(idxs...);
		return _data[fIds];
	}

	T valueUnchecked(ShapeBlock const& idxs) const {
		int fIds = flatIndex<AccessCheck::Nocheck>(idxs);
		return _data[fIds];
	}

	T valueOrAlt(ShapeBlock const& idxs, T const& alt) const {

		int fIds = flatIndex<AccessCheck::Nocheck>(idxs);
		bool inBound = true;
		for (int i = 0; i < nDim; i++) {
			inBound = inBound and idxs[i] >= 0 and idxs[i] < _shape[i];
		}
		fIds = (inBound) ? fIds : 0;
		T const& val = _data[fIds];
		return (inBound) ? val : alt;

	}
	
protected:

	ShapeBlock _shape;
	ShapeBlock _strides;

	T* _data;
	bool _manageData;
	
};

template<class OutStream, typename T, array_size_t nDim>
OutStream& printArray(OutStream& out, Array<T, nDim> const& array, int nPreSpaces = 0) {


	for (int i = 0; i < nPreSpaces; i++) {
		out << ' ';
	}
	out << '[' << '\n';
	for (int i = 0; i < array.shape().back(); i++) {
		printArray(out, const_cast<Array<T, nDim>*>(&array)->sliceView(nDim-1,i), nPreSpaces+1);
	}
	for (int i = 0; i < nPreSpaces; i++) {
		out << ' ';
	}
	out << ']' << '\n';

	return out;

}

template<class OutStream, typename T>
OutStream& printArray(OutStream& out, Array<T, 1> const& array, int nPreSpaces = 0) {

	for (int i = 0; i < nPreSpaces; i++) {
		out << ' ';
	}
	out << '[';
	for (int i = 0; i < array.flatLenght(); i++) {
		if (i != 0) {
			out << ' ';
		}
		out << array.template value<AccessCheck::Nocheck>(i);
	}
	out << "]\n";

	return out;

}

template<class OutStream, typename T, array_size_t nDim>
OutStream& operator<<(OutStream& out, Array<T, nDim> const& array) {

	return printArray(out, array);

}
	
} //namespace multidim 

#endif // MULTIDIM_ARRAY_H
