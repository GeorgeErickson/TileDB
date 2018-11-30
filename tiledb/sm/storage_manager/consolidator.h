/**
 * @file   consolidator.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2018 TileDB, Inc.
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
 * This file defines class Consolidator.
 */

#ifndef TILEDB_CONSOLIDATOR_H
#define TILEDB_CONSOLIDATOR_H

#include "tiledb/sm/array/array.h"
#include "tiledb/sm/misc/status.h"
#include "tiledb/sm/storage_manager/open_array.h"

#include <vector>

namespace tiledb {
namespace sm {

class ArraySchema;
class Query;
class StorageManager;
class URI;

/** Handles array consolidation. */
class Consolidator {
 public:
  /* ********************************* */
  /*           TYPE DEFINITIONS        */
  /* ********************************* */

  /** Consolidation configuration parameters. */
  struct ConsolidationConfig {
    /** Attribute buffer size. */
    uint64_t buffer_size_;
    /**
     * Number of consolidation steps performed in a single
     * consolidation invocation.
     */
    uint32_t steps_;
    /** Minimum number of fragments to consolidate in a single step. */
    uint32_t min_frags_;
    /** Maximum number of fragments to consolidate in a single step. */
    uint32_t max_frags_;
    /**
     * Minimum size ratio for two fragments to be considered for
     * consolidation.
     */
    float size_ratio_;

    /** Constructor. */
    ConsolidationConfig() {
      steps_ = constants::consolidation_steps;
      buffer_size_ = constants::consolidation_buffer_size;
      min_frags_ = constants::consolidation_step_min_frags;
      max_frags_ = constants::consolidation_step_max_frags;
      size_ratio_ = constants::consolidation_step_size_ratio;
    }
  };

  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  /**
   * Constructor.
   *
   * @param storage_manager The storage manager.
   */
  explicit Consolidator(StorageManager* storage_manager);

  /** Destructor. */
  ~Consolidator();

  /* ********************************* */
  /*                API                */
  /* ********************************* */

  /**
   * Consolidates the fragments of the input array.
   *
   * @param array_name URI of array to consolidate.
   * @param encryption_type The encryption type of the array
   * @param encryption_key If the array is encrypted, the private encryption
   *    key. For unencrypted arrays, pass `nullptr`.
   * @param key_length The length in bytes of the encryption key.
   * @param config Configuration parameters for the consolidation
   *     (`nullptr` means default).
   * @return Status
   */
  Status consolidate(
      const char* array_name,
      EncryptionType encryption_type,
      const void* encryption_key,
      uint32_t key_length,
      const Config* config);

 private:
  /* ********************************* */
  /*        PRIVATE ATTRIBUTES         */
  /* ********************************* */

  /** Connsolidation configuration parameters. */
  ConsolidationConfig config_;

  /** The storage manager. */
  StorageManager* storage_manager_;

  /* ********************************* */
  /*          PRIVATE METHODS           */
  /* ********************************* */

  /**
   * Consolidates the input fragments of the input array. This function
   * implements a single consolidation step.
   *
   * @param array_name URI of array to consolidate.
   * @param to_consolidate The fragments to consolidate in this consolidation
   *     step.
   * @param encryption_type The encryption type of the array
   * @param encryption_key If the array is encrypted, the private encryption
   *    key. For unencrypted arrays, pass `nullptr`.
   * @param key_length The length in bytes of the encryption key.
   * @param new_fragment_uri The URI of the fragment created after
   *     consolidating the `to_consolidate` fragments.
   * @return Status
   */
  Status consolidate(
      const URI& array_uri,
      const std::vector<FragmentInfo>& to_consolidate,
      EncryptionType encryption_type,
      const void* encryption_key,
      uint32_t key_length,
      URI* new_fragment_uri);

  /**
   * Copies the array by reading from the fragments to be consolidated
   * (with `query_r`) and writing to the new fragment (with `query_w`).
   *
   * @param query_r The read query.
   * @param query_w The write query.
   * @return Status
   */
  Status copy_array(Query* query_r, Query* query_w);

  /** Cleans up the inputs. */
  void clean_up(
      void* subarray,
      unsigned buffer_num,
      void** buffers,
      uint64_t* buffer_sizes,
      Query* query_r,
      Query* query_w) const;

  /**
   * Creates the buffers that will be used upon reading the input fragments and
   * writing into the new fragment. It also retrieves the number of buffers
   * created.
   *
   * @param array_schema The array schema.
   * @param buffers The buffers to be created.
   * @param buffer_sizes The corresponding buffer sizes.
   * @param buffer_num The number of buffers to be retrieved.
   * @return Status
   */
  Status create_buffers(
      const ArraySchema* array_schema,
      void*** buffers,
      uint64_t** buffer_sizes,
      unsigned int* buffer_num);

  /**
   * Creates the queries needed for consolidation. It also retrieves
   * the number of fragments to be consolidated and the URI of the
   * new fragment to be created.
   *
   * @param array_for_reads The opened array for reading the fragments
   *     to be consolidated.
   * @param array_for_writes The opened array for writing the
   *     consolidated fragments.
   * @param subarray The subarray to read from (the fragments to consolidate)
   *     and write to (the new fragment).
   * @param layout The layout to read from and write to.
   * @param buffers The buffers to be passed in the queries.
   * @param buffer_sizes The corresponding buffer sizes.
   * @param query_r This query reads from the fragments to be consolidated.
   * @param query_w This query writes to the new consolidated fragment.
   * @param new_fragment_uri The URI of the new fragment to be created.
   * @return Status
   */
  Status create_queries(
      Array* array_for_reads,
      Array* array_for_writes,
      void* subarray,
      Layout layout,
      void** buffers,
      uint64_t* buffer_sizes,
      Query** query_r,
      Query** query_w,
      URI* new_fragment_uri);

  /**
   * Computes the subarray issued as a read and write query during
   * consolidation, as well as the layout the query should be issued
   * in. The decision is taken based on whether the input fragments
   * to consolidate are dense or sparse, as well as their non-empty
   * domains.
   */
  Status compute_subarray_and_layout(
      Array* array,
      const std::vector<FragmentInfo>& to_consolidate,
      void** subarray,
      Layout* layout) const;

  /**
   * Computes the non-empty domain of the input fragments to be consolidated.
   * Note that the input contains at least one dense fragment. The computed
   * non-empty domain is the union of all fragments to consolidate.
   *
   * @param array_schema The array schema.
   * @param to_consolidate The fragments to consolidate.
   * @param non_empty_domain The non-empty domain to compute.
   * @param coincides This is set to `true` by the function if the
   *     computed non-empty domain coincides with the array space tiles.
   * @return Status
   */
  Status compute_non_empty_domain(
      ArraySchema* array_schema,
      const std::vector<FragmentInfo>& to_consolidate,
      void* non_empty_domain,
      bool* coincides) const;

  /**
   * Computes the non-empty domain of the input fragments to be consolidated.
   * Note that the input contains at least one dense fragment. The computed
   * non-empty domain is the union of all fragments to consolidate.
   *
   * @tparam T The domain type.
   * @param array_schema The array schema.
   * @param to_consolidate The fragments to consolidate.
   * @param non_empty_domain The non-empty domain to compute.
   * @param coincides This is set to `true` by the function if the
   *     computed non-empty domain coincides with the array space tiles.
   * @return Status
   */
  template <class T>
  Status compute_non_empty_domain(
      ArraySchema* array_schema,
      const std::vector<FragmentInfo>& to_consolidate,
      T* non_empty_domain,
      bool* coincides) const;

  /**
   * Deletes the fragment metadata files of the old fragments that
   * got consolidated. This renders the old fragments "invisible".
   *
   * @param fragments The old fragments.
   * @return Status
   */
  Status delete_old_fragment_metadata(
      const std::vector<FragmentInfo>& fragments);

  /**
   * Deletes the old fragments that got consolidated.
   *
   * @param fragments The old fragments.
   * @return Status
   */
  Status delete_old_fragments(const std::vector<FragmentInfo>& fragments);

  /**
   * Frees the input buffers.
   *
   * @param buffer_num The number of buffers.
   * @param buffers The buffers to be freed.
   * @param buffer_sizes The corresponding buffer sizes.
   */
  void free_buffers(
      unsigned int buffer_num, void** buffers, uint64_t* buffer_sizes) const;

  /**
   * Based on the input fragment info, this algorithm decides the (sorted) list
   * of fragments to be consolidated in the next consolidation step.
   *
   * @param fragments Information about all the fragments.
   * @param to_consolidate The fragments to consolidate in the next step.
   * @return Status
   */
  Status compute_next_to_consolidate(
      const std::vector<FragmentInfo>& fragments,
      std::vector<FragmentInfo>* to_consolidate) const;

  /**
   * Renames the new fragment URI. The new name has the format
   * `__<thread_id>_<timestamp>_<last_fragment_timestamp>`, where
   * `<thread_id>` is the id of the thread that performs the consolidation,
   * `<timestamp>` is the current timestamp in milliseconds, and
   * `<last_fragment_timestamp>` is the timestamp of the last of the
   * consolidated fragments.
   */
  Status rename_new_fragment_uri(URI* uri) const;

  /** Checks and sets the input configuration parameters. */
  Status set_config(const Config* config);

  /**
   * Sets the buffers to the query, using all the attributes in the
   * query schema. There is a 1-1 correspondence between the input `buffers`
   * and the attributes in the schema, considering also the coordinates
   * if the array is sparse in the end.
   */
  Status set_query_buffers(
      Query* query, void** buffers, uint64_t* buffer_sizes) const;

  /**
   * Updates the `fragment_info` by removing `to_consolidate` and
   * replacing those fragment info objects with `new_fragment_info`.
   *
   * @param to_consolidate Fragment info objects to remove from `fragment_info`.
   * @param new_fragment_info The new object that replaces `to_consolidate`
   *     in `fragment_info`.
   * @param fragment_info The fragment info vector to be updated.
   * @return void
   *
   * @note The algorithm assumes that to_consolidate is a sub vector of
   *     `fragment_info` of size >=0.
   */
  void update_fragment_info(
      const std::vector<FragmentInfo>& to_consolidate,
      const FragmentInfo& new_fragment_info,
      std::vector<FragmentInfo>* fragment_info) const;
};

}  // namespace sm
}  // namespace tiledb

#endif  // TILEDB_FRAGMENT_H
