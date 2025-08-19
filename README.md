# MultidimArrays
Simple multi-dimensional arrays in c++

The aim of this simple header library is to provide a templated class to manage multi dimensional buffers. 
Basic constructors will assign dense buffers using the new [] operator, making the class behave as a generic tensor class.
You also have access to a constructor to map any buffer, included memory aligned buffers, gpu buffers, ...

The library also implement some utils classes to manipulate arrays' and tensors' indices in arbitrary dimensions.

Move semantics is supported.

## Usage

The library is header only. It can be built and installed using cmake. Alternatively, you can also use it in cmake projects using the FetchContent cmake module.

To use the library in cmake directly using FetchContent write:

```
include(FetchContent) #if you have not done so before
FetchContent_Declare(
  multidimarrays
  GIT_REPOSITORY https://github.com/french-paragon/MultidimArrays.git
  GIT_TAG [select your tag here]
)
FetchContent_MakeAvailable(multidimarrays)
```

in your CMakeLists.txt

You can then link against the Multidim::Arrays target and the corresponding include directories will be set for your targets.

To use the library just include the corresponding headers:

```
#include <MultidimArrays/MultidimArrays.h> //The Array class
#include <MultidimArrays/MultidimIndexManipulators.h> //Manipulators and other tools to play with indices
```
