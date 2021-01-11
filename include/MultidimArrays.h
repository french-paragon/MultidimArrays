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
#include <assert.h>

namespace Multidim {
	
template<int d>
constexpr bool allDimsValid() {
	return d > 0;
}

template<int d, int... ds>
constexpr bool allDimsValid() {
	return allDimsValid<d>() and allDimsValid<ds...>();
}

bool allDimsValid(int d) {
	return d > 0;
}


template<typename ... Ds>
bool allDimsValid(int d, Ds... dims) {
	return allDimsValid(d) and allDimsValid(dims...);
}

enum class AccessCheck {
	Check = 1,
	Nocheck = 0
};

typedef int array_size_t;

template<typename T, array_size_t nDim>
class Array {

	static_assert(nDim > 0,
                  "No dimensions provided");

public:

	typedef std::array<array_size_t, nDim> ShapeBlock;

	template<array_size_t... dims>
	Array() :
		_shape({dims...})
	{
		static_assert(sizeof...(dims) == nDim,
				"Wrong number of dimensions provided");

		static_assert(allDimsValid<dims...>(),
				"All dimensions must be greather or equal zero.");

		int s = 1;
		for (int i = 0; i < nDim; i++) {
			_strides[i] = s;
			s *= _shape[i];
		}

		_manageData = true;
		_data = new T[s];
	}

	template<typename... Ds>
	Array(Ds... dims) :
		_shape({dims...})
	{
		static_assert(sizeof...(dims) == nDim,
				"Wrong number of dimensions provided");

		assert(allDimsValid(dims...));

		int s = 1;
		for (int i = 0; i < nDim; i++) {
			_strides[i] = s;
			s *= _shape[i];
		}

		_manageData = true;
		_data = new T[s];
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

			_manageData = true;
			_data = new T[s];

			memcpy(_data, other._data, sizeof (T)*s);
		} else {
			_manageData = false;
			_data = other._data;
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

	template<AccessCheck c = AccessCheck::Check, typename... Ds>
	int flatIndex(Ds... idxs) const {

		static_assert(sizeof...(idxs) == nDim,
				"Wrong number of indices provided");

		ShapeBlock idxs_array = {idxs...};
		return flatIndex<c>(idxs_array);

	}

	template<AccessCheck c = AccessCheck::Check>
	int flatIndex(ShapeBlock idxs) const {

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
	T& at(Ds... idxs) {
		int fIds = flatIndex<c>(idxs...);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check>
	T& at(ShapeBlock idxs) {
		int fIds = flatIndex<c>(idxs);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check, typename... Ds>
	T value(Ds... idxs) const {
		int fIds = flatIndex<c>(idxs...);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check, typename... Ds>
	T value(ShapeBlock idxs) const {
		int fIds = flatIndex<c>(idxs);
		return _data[fIds];
	}
	
protected:

	ShapeBlock _shape;
	ShapeBlock _strides;

	T* _data;
	bool _manageData;
	
};
	
} //namespace multidim 

#endif // MULTIDIM_ARRAY_H
