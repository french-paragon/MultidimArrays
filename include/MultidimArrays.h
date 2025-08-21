#ifndef MULTIDIM_ARRAY_H
#define MULTIDIM_ARRAY_H

/*
Copyright 2021-2023 Paragon<french.paragon@gmail.com>

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

// at this point if MULTIDIMARRAYS_STRONG_BOUND_CHECK is defined, then all access to arrays are checked.
// this is usefull for debug purposes

#include <array>
#include <stdexcept>
#include <cstring>
#include <tuple>
#include <assert.h>
#include <iostream>
#include <algorithm>

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

enum ArrayDataAccessConstness {
	ConstView = 1,
	NonConstView = 0
};

constexpr ArrayDataAccessConstness oppositeConstness(ArrayDataAccessConstness constness) {
	if (constness == ConstView) {
		return NonConstView;
	}
	return ConstView;
}

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

template<typename... slices>
constexpr std::array<bool, sizeof... (slices)> slicedMask() {
	return {(slicedDims<slices>() == 1)...};
}

template<typename T, array_size_t nDim, ArrayDataAccessConstness viewConstness = NonConstView>
class Array {

	static_assert(nDim > 0,
                  "No dimensions provided");

	template<typename fT, array_size_t fNDim, ArrayDataAccessConstness fViewConstness>
	friend class Array;

public:

	static constexpr bool isConstView() {
		return viewConstness == ConstView;
	}

	typedef std::array<array_size_t, nDim> ShapeBlock;
	typedef T ScalarT;

	class IndexBlock : public ShapeBlock {
	public:
		IndexBlock() :
			ShapeBlock() {

		}
		IndexBlock(ShapeBlock const& other) :
			ShapeBlock(other) {

		}

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

		void setZero() {
			for (int j = 0; j < nDim; j++) {
				(*this)[j] = 0;
			}
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
	/*!
	 * \brief Array constructor with shape known at compile time.
	 */
	Array() :
		_shape({dims...})
	{
		static_assert(sizeof...(dims) == nDim or sizeof...(dims) == 0,
				"Wrong number of dimensions provided");

		static_assert(viewConstness == NonConstView or sizeof...(dims) == 0,
				"Cannot build a non-empty const view array without data !");

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


    template<typename... Ds, std::enable_if_t<std::conjunction_v<std::is_integral<Ds>...>>* = nullptr>
	/*!
	 * \brief Array constructor with set shape
	 * \param dims the sizes of each dimensions.
	 */
	Array(Ds... dims) :
		_shape({dims...})
	{
		static_assert(sizeof...(dims) == nDim or sizeof...(dims) == 0,
				"Wrong number of dimensions provided");

		static_assert(viewConstness == NonConstView or sizeof...(dims) == 0,
				"Cannot build a non-empty const view array without data !");

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

	template<ArrayDataAccessConstness vc = viewConstness>
	/*!
	 * \brief Array constructor with a set shape
	 * \param shape the shape to use.
	 */
	Array(ShapeBlock const& shape, std::enable_if_t<vc == NonConstView and viewConstness == NonConstView, void*> = nullptr) :
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

	template<ArrayDataAccessConstness vc = viewConstness>
	/*!
	 * \brief Array constructor with set shape and strides
	 * \param shape the shape to use
	 * \param strds the strides to use
	 */
	Array(ShapeBlock const& shape, ShapeBlock const& strds, std::enable_if_t<vc == NonConstView and viewConstness == NonConstView, void*> = nullptr) :
		_shape(shape),
		_strides(strds)
	{

		assert(allDimsValid(shape));

		ShapeBlock idx = _shape;
		for (int i = 0; i < nDim; i++) {
			idx[i] = idx[i]-1;
		}
        long s = flatIndex(idx) + 1;
        long v = 1;
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

	/*!
	 * \brief Array constructor from a block of memory
	 * \param data the pointer to the data segment
	 * \param shape the shape of the array
	 * \param strides the strides of the array
	 * \param manage if the array is expected to delete the data segment or not.
	 */
	Array(T* data, ShapeBlock const& shape, ShapeBlock const& strides, bool manage = false) :
		_shape(shape),
		_strides(strides),
		_data(data),
		_manageData(manage)
	{

	}

	/*!
	 * \brief Array copy constructor
	 * \param other the other array to copy.
	 *
	 * Note that the behavior of this constructor is a bit strange for const view arrays.
	 *
	 * If the constructed array is a const view, it will always get constructed as a view of the other array.
	 * This behavior is to make sure that functions taking a const view array as argument will not copy
	 * non const view array passed as argument, which would be really inneficient...
	 *
	 */
	Array(Array<T, nDim, viewConstness>& other) : // TODO : check if we can share ownership of a data segment using smart pointers...
		_shape(other._shape),
		_strides(other._strides)
	{

		if (viewConstness == ConstView) { //construct a view
			_data = other._data;
			_manageData = false;
			return;
		}

		if (other._manageData) {

			if (!other.isDense()) {
				int s = 1;
				for (int i = 0; i < nDim; i++) {
					_strides[i] = s;
					s *= _shape[i];
				}
			}

			int s = flatLenght();

			if (s > 0) {
				_manageData = true;
				_data = new T[s];

				if (other.isDense()) {

					memcpy(_data, other._data, sizeof (T)*s);

				} else {

					IndexBlock idx;
					idx.setZero();

					for (int i = 0; i < s; i++) {
						atUnchecked(idx) = other.valueUnchecked(idx);
						idx.moveToNextIndex(_shape);
					}
				}

			} else {
				_manageData = false;
				_data = nullptr;
			}
		} else {
			_manageData = false;
			_data = other._data;
		}
	}

	/*!
	 * \brief Array constructor with array of opposite constness
	 * \param other the other array to copy.
	 *
	 * Note that the behavior of this constructor is a bit strange for const view arrays.
	 *
	 * If the constructed array is a const view, it will always get constructed as a view of the other array.
	 * This behavior is to make sure that functions taking a const view array as argument will not copy
	 * non const view array passed as argument, which would be really inneficient...
	 *
	 * Also note that we need to explicitly declare this function, if we use a template
	 * for the data view constness, then some function of the qt meta object system get broken.
	 *
	 */
	Array(Array<T, nDim, oppositeConstness(viewConstness)>& other) : // TODO : check if we can share ownership of a data segment using smart pointers...
		_shape(other._shape),
		_strides(other._strides)
	{

		if (viewConstness == ConstView) { //construct a view
			_data = other._data;
			_manageData = false;
			return;
		}

		if (other._manageData) {

			if (!other.isDense()) {
				int s = 1;
				for (int i = 0; i < nDim; i++) {
					_strides[i] = s;
					s *= _shape[i];
				}
			}

			int s = flatLenght();

			if (s > 0) {
				_manageData = true;
				_data = new T[s];

				if (other.isDense()) {

					memcpy(_data, other._data, sizeof (T)*s);

				} else {

					IndexBlock idx;
					idx.setZero();

					for (int i = 0; i < s; i++) {
						_data[flatIndex(idx)] = other.valueUnchecked(idx);
						idx.moveToNextIndex(_shape);
					}
				}

			} else {
				_manageData = false;
				_data = nullptr;
			}
		} else {
			_manageData = false;
			_data = other._data;
		}
	}

	Array(Array<T, nDim, viewConstness> const& other) :
		_shape(other._shape),
		_strides(other._strides)
	{

		if (viewConstness == ConstView) { //construct a view
			_data = other._data;
			_manageData = false;
			return;
		}

		if (!other.isDense()) {
			int s = 1;
			for (int i = 0; i < nDim; i++) {
				_strides[i] = s;
				s *= _shape[i];
			}
		}

		int s = flatLenght();

		if (s > 0) {
			_manageData = true;
			_data = new T[s];

			if (other.isDense()) {

				memcpy(_data, other._data, sizeof (T)*s);

			} else {

				IndexBlock idx;
				idx.setZero();

				for (int i = 0; i < s; i++) {
					_data[flatIndex(idx)] = other.valueUnchecked(idx);
					idx.moveToNextIndex(_shape);
				}
			}

		} else {
			_manageData = false;
			_data = nullptr;
		}
	}

	Array(Array<T, nDim, oppositeConstness(viewConstness)> const& other) :
		_shape(other._shape),
		_strides(other._strides)
	{

		if (viewConstness == ConstView) { //construct a view
			_data = other._data;
			_manageData = false;
			return;
		}

		if (!other.isDense()) {
			int s = 1;
			for (int i = 0; i < nDim; i++) {
				_strides[i] = s;
				s *= _shape[i];
			}
		}

		int s = flatLenght();

		if (s > 0) {
			_manageData = true;
			_data = new T[s];

			if (other.isDense()) {

				memcpy(_data, other._data, sizeof (T)*s);

			} else {

				IndexBlock idx;
				idx.setZero();

				for (int i = 0; i < s; i++) {
					_data[flatIndex(idx)] = other.valueUnchecked(idx);
					idx.moveToNextIndex(_shape);
				}
			}

		} else {
			_manageData = false;
			_data = nullptr;
		}
	}


	template<ArrayDataAccessConstness otherViewConstness,
			 std::enable_if_t<viewConstness == ConstView or otherViewConstness == NonConstView, bool> = true>
	/*!
	 * \brief Array move constructor for const view array
	 * \param other the movable array to copy.
	 *
	 * This constructor is available only for const view array, or when the movable array is non const view,
	 * else it would be possible to build a non const view array from a const view one by simply using std::move...
	 */
	Array(Array<T, nDim, otherViewConstness>&& other) :
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

	Array& operator=(Array<T, nDim>& other) {
		if (_manageData and _data != nullptr) {
			delete [] _data;
		}
		_shape = other.shape();

		if (other.isDense()) {
			_strides = other.strides();
		} else {
			int s = 1;
			for (int i = 0; i < nDim; i++) {
				_strides[i] = s;
				s *= _shape[i];
			}
		}

		if (other._manageData) {

			int s = flatLenght();

			if (s > 0) {
				_manageData = true;
				_data = new T[s];

				if (other.isDense()) {

					memcpy(_data, other._data, sizeof (T)*s);

				} else {

					IndexBlock idx;
					idx.setZero();

					for (int i = 0; i < s; i++) {
						atUnchecked(idx) = other.atUnchecked(idx);
						idx.moveToNextIndex(_shape);
					}
				}
			} else {
				_manageData = false;
				_data = nullptr;
			}
		} else {
			_manageData = false;
			_data = other._data;
		}
		return *this;
	}

	Array& operator=(Array<T, nDim> const& other) {
		if (_manageData and _data != nullptr) {
			delete [] _data;
		}
		_shape = other.shape();

		if (other.isDense()) {
			_strides = other.strides();
		} else {
			int s = 1;
			for (int i = 0; i < nDim; i++) {
				_strides[i] = s;
				s *= _shape[i];
			}
		}

		int s = flatLenght();

		if (s > 0) {
			if (other.isDense()) {
				_manageData = true;
				_data = new T[s];

				memcpy(_data, other._data, sizeof (T)*s);

			} else {

				_manageData = true;
				_data = new T[s];

				IndexBlock idx;
				idx.setZero();

				for (int i = 0; i < s; i++) {
					atUnchecked(idx) = other.valueUnchecked(idx);
					idx.moveToNextIndex(_shape);
				}
			}
		} else {
			_manageData = false;
			_data = nullptr;
		}
		return *this;
	}

	Array& operator=(Array<T, nDim>&& other) {
		if (_manageData and _data != nullptr) {
			delete [] _data;
		}
		_shape = other.shape();
		_strides = other.strides();
		_manageData = other._manageData;
		_data = other._data;

		other._data = nullptr;
		return *this;
	}

	Array& operator=(const Array<T, nDim>&& other) {
		if (_manageData and _data != nullptr) {
			delete [] _data;
		}
		_shape = other.shape();
		_strides = other.strides();
		_manageData = other._manageData;
		_data = other._data;

		other._data = nullptr;
		return *this;
	}

	ShapeBlock shape() const {
		return _shape;
	}

    long flatLenght() const {

        long s = 1;
		for (int i = 0; i < nDim; i++) {
			s *= _shape[i];
		}

		return s;
	}

	ShapeBlock strides() const {
		return _strides;
	}

	bool isView() const {
		return !_manageData;
	}

	bool isDense() const {

		ShapeBlock idx = _shape;

		for (int i = 0; i < nDim; i++) {

			if (_strides[i] < 1) {
				return false;
			}

			idx[i] -= 1;
		}

		return flatIndex(idx) == flatLenght()-1;
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

    inline array_size_t sliceFirstIndex(DimSlice const& slice) const {
		return slice.startIndex;
	}

    inline array_size_t sliceFirstIndex(DimIndex const& index) const {
		return index.index;
	}

	template<typename... slices>
    inline ShapeBlock firstIndex(slices... s) const {
		ShapeBlock blk = {sliceFirstIndex(s)...};
		for (int i = 0; i < nDim; i++) {
			if (blk[i] < 0) {
				blk[i] = _shape[i] + blk[i];
			}
		}
		return blk;
	}

	template<typename... slices>
    inline long dataOffset(slices... s) const {
		return flatIndex(firstIndex(s...));
	}

	template<typename... slices>
    inline typename Array<T, slicedDims<slices...>()>::ShapeBlock getSubShape(slices... s) const {

		using SBlock = typename Array<T, slicedDims<slices...>()>::ShapeBlock;

		std::array<const DimInfo *, sizeof... (slices)> a = {static_cast<const DimInfo *>(&s)...};
		constexpr std::array<bool, sizeof... (slices)> slicedDimsMask = slicedMask<slices...>();
		SBlock r;

		int p = 0;
		for (int i = 0; i < sizeof... (slices); i++) {
			if (slicedDimsMask[i]) {
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
    inline typename Array<T, slicedDims<slices...>()>::ShapeBlock getSubStrides(slices... s) const {

		using SBlock = typename Array<T, slicedDims<slices...>()>::ShapeBlock;

		std::array<DimInfo const*, sizeof... (slices)> a = {static_cast<DimInfo const*>(&s)...};
		constexpr std::array<bool, sizeof... (slices)> slicedDimsMask = slicedMask<slices...>();
		SBlock r;

		int p = 0;
		for (int i = 0; i < sizeof... (slices); i++) {
			if (slicedDimsMask[i]) {

				const DimSlice * s = static_cast<const DimSlice *>(a[i]);

				r[p] = _strides[i]*s->indexJump;
				p++;
			}
		}

		return r;
	}

public:

	template<typename... slices>
	inline Array<T, slicedDims<slices...>(), viewConstness> subView(slices... s) {

		static_assert (sizeof... (slices) == nDim, "Cannot generate an array with zero dimensions");
		static_assert (slicedDims<slices...>() > 0, "Cannot generate an array with zero dimensions");

		using SubArray = Array<T, slicedDims<slices...>(), viewConstness>;
		using SubShapeBlock = typename SubArray::ShapeBlock;

        long offset = dataOffset(s...);
		SubShapeBlock shape = getSubShape(s...);
		SubShapeBlock strides = getSubStrides(s...);

		return SubArray(_data + offset, shape, strides);
	}

	template<typename... slices>
	inline Array<T, slicedDims<slices...>(), ConstView> subView(slices... s) const {

		static_assert (sizeof... (slices) == nDim, "Cannot generate an array with zero dimensions");
		static_assert (slicedDims<slices...>() > 0, "Cannot generate an array with zero dimensions");

		using SubArray = Array<T, slicedDims<slices...>(), ConstView>;
		using SubShapeBlock = typename SubArray::ShapeBlock;

        long offset = dataOffset(s...);
		SubShapeBlock shape = getSubShape(s...);
		SubShapeBlock strides = getSubStrides(s...);

		return SubArray(_data + offset, shape, strides);
	}

	/*!
	 * \brief sliceView return a view of the array were a single dimension is indexed.
	 * \param dim the dim to index
	 * \param coord the index for dim
	 * \return a view of the array with one less dimension.
	 */
	inline Array<T, std::max(1,nDim-1), viewConstness> sliceView(int dim, int coord) {

		if (dim < 0 or dim >= nDim) {
			throw std::out_of_range("Dim index out of range");
		}

		if (coord < 0 or coord >= shape()[dim]) {
			throw std::out_of_range("Index out of range");
		}

		using SubArray = Array<T, std::max(1,nDim-1), viewConstness>;
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

	inline Array<T, std::max(1,nDim-1), ConstView> sliceView(int dim, int coord) const {

		if (dim < 0 or dim >= nDim) {
			throw std::out_of_range("Dim index out of range");
		}

		if (coord < 0 or coord >= shape()[dim]) {
			throw std::out_of_range("Index out of range");
		}

		using SubArray = Array<T, std::max(1,nDim-1), ConstView>;
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
	inline Array<T, 1, viewConstness> indexDimView(int dim, std::array<array_size_t,nDim-1> const& coords) {

		if (dim < 0 or dim >= nDim) {
			throw std::out_of_range("Dim index out of range");
		}

		ShapeBlock startIdx;

		for (int i = 0; i < nDim; i++) {

			if (i < dim) {
				startIdx[i] = coords[i];
			}

			if (i == dim) {
				startIdx[i] = 0;
			}

			if (i > dim) {
				startIdx[i] = coords[i-1];
			}
		}

        long flatidx = flatIndex<c>(startIdx);

		using SubArray = Array<T, 1, viewConstness>;

		return SubArray(&_data[flatidx], {_shape[dim]}, {_strides[dim]});
	}

	template<AccessCheck c = AccessCheck::Check>
	inline Array<T, 1, ConstView> indexDimView(int dim, std::array<array_size_t,nDim-1> const& coords) const {

		if (dim < 0 or dim >= nDim) {
			throw std::out_of_range("Dim index out of range");
		}

		ShapeBlock startIdx;

		for (int i = 0; i < nDim; i++) {

			if (i < dim) {
				startIdx[i] = coords[i];
			}

			if (i == dim) {
				startIdx[i] = 0;
			}

			if (i > dim) {
				startIdx[i] = coords[i-1];
			}
		}

        long flatidx = flatIndex<c>(startIdx);

		using SubArray = Array<T, 1, ConstView>;

		return SubArray(&_data[flatidx], {_shape[dim]}, {_strides[dim]});
	}

	template<AccessCheck c = AccessCheck::Check>
    inline long flatIndex(ShapeBlock const& idxs) const {

        #ifdef MULTIDIMARRAYS_STRONG_BOUND_CHECK
        for (int i = 0; i < nDim; i++) {
            if (idxs[i] < 0 or idxs[i] >= _shape[i]) {
                throw std::out_of_range("Index out of range");
            }
        }
        #endif

		if (c == AccessCheck::Check) {
			for (int i = 0; i < nDim; i++) {
				if (idxs[i] < 0 or idxs[i] >= _shape[i]) {
					throw std::out_of_range("Index out of range");
				}
			}
		}

        long fIds = 0;

		for (int i = 0; i < nDim; i++) {
            fIds += static_cast<long>(idxs[i])*_strides[i];
		}

		return fIds;

	}

	template<AccessCheck c = AccessCheck::Check, typename... Ds>
    inline long flatIndex(int idx0, Ds... idxs) const {

		static_assert(sizeof...(idxs) == nDim -1 or sizeof...(idxs) == 0,
				"Wrong number of indices provided");

		if (sizeof...(idxs) == 0) {

            #ifdef MULTIDIMARRAYS_STRONG_BOUND_CHECK
            if (idx0 < 0 or idx0 >= flatLenght()) {
                throw std::out_of_range("Index out of range");
            }
            #endif

			if (c == AccessCheck::Check) {
				if (idx0 < 0 or idx0 >= flatLenght()) {
					throw std::out_of_range("Index out of range");
				}
			}

            return static_cast<long>(idx0)*_strides[0];

		} else {

            ShapeBlock idxs_array = {idx0, static_cast<int>(idxs)...};
			return flatIndex<c>(idxs_array);

		}

	}

    /*!
     * \brief indexFromFlat convert a flat index back to a multidimensional shape block
     * \param index the flat index
     * \return a ShapeBlock which, once converted with flatIndex, should give back the original flat index.
     */
    ShapeBlock indexFromFlat(long index) const {

		ShapeBlock out;
        long leftOver = index;

        std::array<int, nDim> strides_order;

        for (int i = 0; i < nDim; i++) {
            strides_order[i] = i;
        }

        std::sort(strides_order.begin(), strides_order.end(),
                  [this] (int i1, int i2) { return _strides[i1] < _strides[i2]; });

		for (int i = nDim-1; i >= 0; i--) {

			//garanteed optimizable to branchless !
            out[strides_order[i]] = leftOver/_strides[strides_order[i]];
            leftOver %= _strides[strides_order[i]];
		}

		return out;
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

	template<AccessCheck c = AccessCheck::Check, ArrayDataAccessConstness vc = viewConstness, typename... Ds>
	inline std::enable_if_t<vc == NonConstView, T&> at(Ds... idxs) {
		static_assert (vc == viewConstness, "invalid view constness");
        long fIds = flatIndex<c>(idxs...);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check, ArrayDataAccessConstness vc = viewConstness>
	inline std::enable_if_t<vc == NonConstView, T&> at(ShapeBlock const& idxs) {
		static_assert (vc == viewConstness, "invalid view constness");
        long fIds = flatIndex<c>(idxs);
		return _data[fIds];
	}

	template<ArrayDataAccessConstness vc = viewConstness, typename... Ds>
	inline std::enable_if_t<vc == NonConstView, T&> atUnchecked(Ds... idxs) {
		static_assert (vc == viewConstness, "invalid view constness");
        long fIds = flatIndex<AccessCheck::Nocheck>(idxs...);
		return _data[fIds];
	}

	template<ArrayDataAccessConstness vc = viewConstness>
	inline std::enable_if_t<vc == NonConstView, T&> atUnchecked(ShapeBlock const& idxs) {
		static_assert (vc == viewConstness, "invalid view constness");
        long fIds = flatIndex<AccessCheck::Nocheck>(idxs);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check, typename... Ds>
	inline T value(Ds... idxs) const {
        long fIds = flatIndex<c>(idxs...);
		return _data[fIds];
	}

	template<AccessCheck c = AccessCheck::Check>
	inline T value(ShapeBlock const& idxs) const {
        long fIds = flatIndex<c>(idxs);
		return _data[fIds];
	}

	template<typename... Ds>
	inline T valueUnchecked(Ds... idxs) const {
        long fIds = flatIndex<AccessCheck::Nocheck>(idxs...);
		return _data[fIds];
	}

	inline T valueUnchecked(ShapeBlock const& idxs) const {
        long fIds = flatIndex<AccessCheck::Nocheck>(idxs);
		return _data[fIds];
	}

    inline T valueOrAlt(ShapeBlock const& idxs, T const& alt) const {
		bool inBound = true;
		for (int i = 0; i < nDim; i++) {
			inBound = inBound and idxs[i] >= 0 and idxs[i] < _shape[i];
		}

        #ifdef MULTIDIMARRAYS_STRONG_BOUND_CHECK
        if (!inBound) { //not branchless, but necessary when bound check is forced even with AccessCheck::Nocheck
            return alt;
        }
        #endif

        long fIds = flatIndex<AccessCheck::Nocheck>(idxs);
		fIds = (inBound) ? fIds : 0;
		T const& val = _data[fIds];
		return (inBound) ? val : alt;

	}

	template<ArrayDataAccessConstness vc = viewConstness>
	/*!
	 * \brief takePointer allows for an external operator to take ownership of the managed pointer
	 * \return the data pointer if the multidimensional array do manage it, else nullptr
	 *
	 * This function is meant for systems that needs to transform a multidimensional array into a different class using some kind of buffer protocol
	 * and make the new class responsible for managing the data.
	 *
	 * In most cases it is not advisable to use it, except if you know exactly what you are doing, as
	 * misuse can lead to memory leak or corruption.
	 */
	std::enable_if_t<vc == NonConstView, T*> takePointer() {

		if (_manageData) {
			_manageData = false;
			return _data;
		}
		return nullptr;
	}

	template<typename T_O>
	/*!
	 * \brief cast generate an array of a different type than the original array.
	 * \return the casted array
	 *
	 * An edge case occur when T_O is the same type as T.
	 * If the template type T_O is equal to the array type t, the array is moved to the returned array.
	 * The old array will still keep a reference to the data, but the new array will be in charge of managing the memory.
	 *
	 */
	inline Multidim::Array<T_O, nDim> cast() {

        if (empty()) {
            return Multidim::Array<T_O, nDim>();
        }

		if (std::is_same<T_O, T>::value) {

			Multidim::Array<T_O,nDim> casteted(reinterpret_cast<T_O*>(_data), shape(), strides(), _manageData);
			_manageData = false;
			return casteted;
		}

		Multidim::Array<T_O,nDim> casteted(shape(), strides());

        if (isDense()) {
            #pragma omp parallel for
            for (int i = 0; i < flatLenght(); i++) {
                casteted._data[i] = static_cast<T_O>(_data[i]);
            }

            return casteted;
        }

		#pragma omp parallel for
		for (int i = 0; i < flatLenght(); i++) {
			ShapeBlock idx = indexFromFlat(i);
			casteted.atUnchecked(idx) = static_cast<T_O>(valueUnchecked(idx));
		}

		return casteted;
	}

	template<typename T_O>
	/*!
	 * \brief cast generate an array of a different type than the original array.
	 * \return the casted array
	 *
     * The const variant is guaranteed to copy the data, even is T_O and T are the same.
	 *
	 */
	inline Multidim::Array<T_O, nDim> cast() const {

        if (empty()) {
            return Multidim::Array<T_O, nDim>();
        }

		Multidim::Array<T_O,nDim> casteted(shape(), strides());

        if (isDense()) {
            #pragma omp parallel for
            for (int i = 0; i < flatLenght(); i++) {
                casteted._data[i] = static_cast<T_O>(_data[i]);
            }

            return casteted;
        }

		#pragma omp parallel for
		for (int i = 0; i < flatLenght(); i++) {
			ShapeBlock idx = indexFromFlat(i);
			casteted.atUnchecked(idx) = static_cast<T_O>(valueUnchecked(idx));
		}

		return casteted;
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

template<typename T, array_size_t nDim>
std::ostream& operator<<(std::ostream& out, Array<T, nDim> const& array) {
	return printArray(out, array);
}
	
} //namespace multidim 

#endif // MULTIDIM_ARRAY_H
