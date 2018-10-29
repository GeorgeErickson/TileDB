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
  /** The creation timestamp of the fragment. */
  uint64_t timestamp_;
  /** The size of the entire fragment directory. */
  uint64_t size_;

  /** Constructor. */
  FragmentInfo(const URI& uri, uint64_t timestamp, uint64_t size)
      : uri_(uri)
      , timestamp_(timestamp)
      , size_(size) {
  }
};

}  // namespace sm
}  // namespace tiledb

#endif  // TILEDB_FRAGMENT_INFO_H
