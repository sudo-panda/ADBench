//--------------------------------------------------------------------*- C++ -*-
// clad - the C++ Clang-based Automatic Differentiator
// version: $Id$
// author:  Vassil Vassilev <vvasilev-at-cern.ch>
//------------------------------------------------------------------------------

#ifndef CLAD_DIFFERENTIATOR
#define CLAD_DIFFERENTIATOR

#include "Array.h"
#include "ArrayRef.h"
#include "BuiltinDerivatives.h"
#include "CladConfig.h"
#include "Tape.h"

#include <assert.h>
#include <stddef.h>

extern "C" {
  int printf(const char* fmt, ...);
  char* strcpy (char* destination, const char* source);
  size_t strlen(const char*);
#if defined(__APPLE__) || defined(_MSC_VER)
  void* malloc(size_t);
  void free(void *ptr);
#else
  void* malloc(size_t) __THROW __attribute_malloc__ __wur;
  void free(void *ptr) __THROW;
#endif
}

namespace clad {
  /// \returns the size of a c-style string
  CUDA_HOST_DEVICE inline unsigned int GetLength(const char* code) {
    unsigned int count;
    const char* code_copy = code;
    #ifdef __CUDACC__
      count = 0;
      while (*code_copy != '\0') {
        count++;
        code_copy++;
      }
    #else
      count = strlen(code_copy);
    #endif
    return count;
  }
  
  /// Tape type used for storing values in reverse-mode AD inside loops.
  template <typename T>
  using tape = tape_impl<T>;

  /// Add value to the end of the tape, return the same value.
  template <typename T>
  CUDA_HOST_DEVICE T push(tape<T>& to, T val) {
    to.emplace_back(val);
    return val;
  }

  /// Add value to the end of the tape, return the same value.
  /// A specialization for clad::array_ref types to use in reverse mode.
  template <typename T, typename U>
  CUDA_HOST_DEVICE clad::array_ref<T> push(tape<clad::array_ref<T>>& to,
                                           U val) {
    to.emplace_back(val);
    return val;
  }

  /// Remove the last value from the tape, return it.
  template <typename T>
  CUDA_HOST_DEVICE T pop(tape<T>& to) {
    T val = to.back();
    to.pop_back();
    return val;
  }

  /// Access return the last value in the tape.
  template <typename T> CUDA_HOST_DEVICE T& back(tape<T>& of) {
    return of.back();
  }
}
#endif // CLAD_DIFFERENTIATOR

// Enable clad after the header was included.
// FIXME: The header inclusion should be made automatic if the pragma is seen.
#pragma clad ON
