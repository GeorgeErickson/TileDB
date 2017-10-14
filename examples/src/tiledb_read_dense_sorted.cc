/**
 * @file   tiledb_read_dense_sorted.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
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
 * It shows how to read from a dense array, constraining the read
 * to a specific subarray and subset of attributes. The cells are copied to the
 * input buffers sorted in row-major order within the selected subarray.
 */

#include <cstdio>
#include "tiledb.h"

int main() {
  // Initialize context with the default configuration parameters
  tiledb_ctx_t* ctx;
  tiledb_ctx_create(&ctx);

  // Subarray and attributes
  int64_t subarray[] = {3, 4, 2, 4};
  const char* attributes[] = {"a1"};

  // Prepare cell buffers
  int buffer_a1[3];
  void* buffers[] = {buffer_a1};
  uint64_t buffer_sizes[] = {sizeof(buffer_a1)};

  // Create query
  tiledb_query_t* query;
  tiledb_query_create(
      ctx,
      &query,
      "my_dense_array",
      TILEDB_READ,
      TILEDB_ROW_MAJOR,
      subarray,
      attributes,
      1,
      buffers,
      buffer_sizes);

  // Loop until no overflow
  printf(" a1\n----\n");
  tiledb_query_status_t status;
  do {
    // Read from array_schema
    printf("Reading cells...\n");
    tiledb_query_submit(ctx, query);

    // Print cell values
    int64_t result_num = buffer_sizes[0] / sizeof(int);
    for (int i = 0; i < result_num; ++i)
      printf("%3d\n", buffer_a1[i]);

    // Get overflow
    tiledb_query_get_attribute_status(ctx, query, "a1", &status);
  } while (status == TILEDB_INCOMPLETE);

  // Clean up
  tiledb_query_free(ctx, query);
  tiledb_ctx_free(ctx);

  return 0;
}