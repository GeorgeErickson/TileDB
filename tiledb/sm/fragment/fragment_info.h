/**
 * @file  fragment_info.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2018 TileDB, Inc.
 * @copyright Copyright (c) 2016 MIT and Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file defines struct FragmentInfo.
 */

#ifndef TILEDB_FRAGMENT_INFO_H
#define TILEDB_FRAGMENT_INFO_H

#include "tiledb/sm/misc/uri.h"

namespace tiledb {
namespace sm {

/** Stores basic information about a fragment. */
struct FragmentInfo {
  /** The fragment URI. */
  URI uri_;
  /** True if the fragment is sparse, and false if it is dense. */
  bool sparse_;
  /** The creation timestamp of the fragment. */
  uint64_t timestamp_;
  /** The size of the entire fragment directory. */
  uint64_t fragment_size_;
  /** The fragment's non-empty domain. */
  void* non_empty_domain_;
  /** The size in bytes of non_empty_domain_.  */
  uint64_t non_empty_domain_size_;

  /** Constructor. */
  FragmentInfo() {
    non_empty_domain_ = nullptr;
  }

  /** Constructor. */
  FragmentInfo(
      const URI& uri,
      bool sparse,
      uint64_t timestamp,
      uint64_t fragment_size,
      const void* non_empty_domain,
      uint64_t non_empty_domain_size)
      : uri_(uri)
      , sparse_(sparse)
      , timestamp_(timestamp)
      , fragment_size_(fragment_size)
      , non_empty_domain_size_(non_empty_domain_size) {
    non_empty_domain_ = std::malloc(non_empty_domain_size);
    std::memcpy(non_empty_domain_, non_empty_domain, non_empty_domain_size);
  }

  /** Copy constructor. */
  FragmentInfo(const FragmentInfo& info) {
    uri_ = info.uri_;
    sparse_ = info.sparse_;
    timestamp_ = info.timestamp_;
    fragment_size_ = info.fragment_size_;
    non_empty_domain_size_ = info.non_empty_domain_size_;
    non_empty_domain_ = std::malloc(non_empty_domain_size_);
    std::memcpy(
        non_empty_domain_, info.non_empty_domain_, non_empty_domain_size_);
  }

  /** Move constructor. */
  FragmentInfo(FragmentInfo&& info) {
    *this = info;
    std::free(info.non_empty_domain_);
    info.non_empty_domain_ = nullptr;
  }

  /** Destructor. */
  ~FragmentInfo() {
    std::free(non_empty_domain_);
  }

  /** Copy assignment operator. */
  FragmentInfo& operator=(const FragmentInfo& info) {
    if (&info != this) {
      uri_ = info.uri_;
      sparse_ = info.sparse_;
      timestamp_ = info.timestamp_;
      fragment_size_ = info.fragment_size_;
      non_empty_domain_size_ = info.non_empty_domain_size_;
      std::free(non_empty_domain_);
      non_empty_domain_ = std::malloc(non_empty_domain_size_);
      std::memcpy(
          non_empty_domain_, info.non_empty_domain_, non_empty_domain_size_);
    }
    return *this;
  }

  /** Move assignment operator. */
  FragmentInfo& operator=(FragmentInfo&& info) {
    if (&info != this) {
      *this = info;
      std::free(info.non_empty_domain_);
      info.non_empty_domain_ = nullptr;
    }
    return *this;
  }
};

}  // namespace sm
}  // namespace tiledb

#endif  // TILEDB_FRAGMENT_INFO_H
