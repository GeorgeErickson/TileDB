/**
 * @file   consolidator.cc
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
 * This file implements the Consolidator class.
 */

#include "tiledb/sm/storage_manager/consolidator.h"
#include "tiledb/sm/fragment/fragment_info.h"
#include "tiledb/sm/misc/logger.h"
#include "tiledb/sm/misc/utils.h"
#include "tiledb/sm/misc/uuid.h"
#include "tiledb/sm/storage_manager/storage_manager.h"

#include <iostream>
#include <sstream>

/* ****************************** */
/*             MACROS             */
/* ****************************** */

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

namespace tiledb {
namespace sm {

/* ****************************** */
/*   CONSTRUCTORS & DESTRUCTORS   */
/* ****************************** */

Consolidator::Consolidator(StorageManager* storage_manager)
    : storage_manager_(storage_manager) {
}

Consolidator::~Consolidator() = default;

/* ****************************** */
/*               API              */
/* ****************************** */

// TODO: Tests
// TODO: Docs: (i) algorithm, (ii) new configs
// TODO: Perf

Status Consolidator::consolidate(
    const char* array_name,
    EncryptionType encryption_type,
    const void* encryption_key,
    uint32_t key_length,
    const Config* config) {
  // Set config parameters
  RETURN_NOT_OK(set_config(config));

  URI array_uri = URI(array_name);
  std::vector<FragmentInfo> to_consolidate;
  auto timestamp = utils::time::timestamp_now_ms();

  EncryptionKey enc_key;
  RETURN_NOT_OK(enc_key.set_key(encryption_type, encryption_key, key_length));

  // Get array schema
  ObjectType object_type;
  RETURN_NOT_OK(storage_manager_->object_type(array_uri, &object_type));
  bool in_cache;
  auto array_schema = (ArraySchema*)nullptr;
  RETURN_NOT_OK(storage_manager_->load_array_schema(
      array_uri, object_type, enc_key, &array_schema, &in_cache));

  // Get fragment info
  std::vector<FragmentInfo> fragment_info;
  RETURN_NOT_OK_ELSE(
      storage_manager_->get_fragment_info(
          array_schema, timestamp, enc_key, &fragment_info),
      delete array_schema);

  uint32_t step = 0;
  do {
    // No need to consolidate if no more than 1 fragment exist
    if (fragment_info.size() <= 1)
      break;

    // Find the next fragments to be consolidated
    RETURN_NOT_OK_ELSE(
        compute_next_to_consolidate(fragment_info, &to_consolidate),
        delete array_schema);

    // Check if there is anything to consolidate
    if (to_consolidate.size() <= 1)
      break;

    // Consolidate the selected fragments
    URI new_fragment_uri;
    RETURN_NOT_OK_ELSE(
        consolidate(
            array_uri,
            to_consolidate,
            encryption_type,
            encryption_key,
            key_length,
            &new_fragment_uri),
        delete array_schema);

    FragmentInfo new_fragment_info;
    RETURN_NOT_OK_ELSE(
        storage_manager_->get_fragment_info(
            array_schema, enc_key, new_fragment_uri, &new_fragment_info),
        delete array_schema);

    // Update fragment info
    update_fragment_info(to_consolidate, new_fragment_info, &fragment_info);

    // Advance number of steps
    ++step;

  } while (step < config_.steps_);

  delete array_schema;

  return Status::Ok();
}

/* ****************************** */
/*        PRIVATE METHODS         */
/* ****************************** */

Status Consolidator::consolidate(
    const URI& array_uri,
    const std::vector<FragmentInfo>& to_consolidate,
    EncryptionType encryption_type,
    const void* encryption_key,
    uint32_t key_length,
    URI* new_fragment_uri) {
  // Open array for reading
  Array array_for_reads(array_uri, storage_manager_);
  RETURN_NOT_OK(array_for_reads.open(
      QueryType::READ,
      to_consolidate,
      encryption_type,
      encryption_key,
      key_length))

  if (array_for_reads.is_empty()) {
    RETURN_NOT_OK(array_for_reads.close());
    return Status::Ok();
  }

  // Open array for writing
  Array array_for_writes(array_uri, storage_manager_);
  RETURN_NOT_OK_ELSE(
      array_for_writes.open(
          QueryType::WRITE, encryption_type, encryption_key, key_length),
      array_for_reads.close());

  // Get schema
  auto array_schema = array_for_reads.array_schema();

  // Create subarray
  void* subarray = nullptr;
  Layout layout;
  bool all_sparse;
  auto st = compute_subarray_and_layout(
      &array_for_reads, to_consolidate, &subarray, &layout, &all_sparse);
  if (!st.ok()) {
    array_for_reads.close();
    array_for_writes.close();
    return st;
  }

  // Prepare buffers
  void** buffers;
  uint64_t* buffer_sizes;
  unsigned int buffer_num;
  st = create_buffers(array_schema, all_sparse, &buffers, &buffer_sizes, &buffer_num);
  if (!st.ok()) {
    array_for_reads.close();
    array_for_writes.close();
    return st;
  }

  // Create queries
  auto query_r = (Query*)nullptr;
  auto query_w = (Query*)nullptr;
  st = create_queries(
      &array_for_reads,
      &array_for_writes,
      all_sparse,
      subarray,
      layout,
      buffers,
      buffer_sizes,
      &query_r,
      &query_w,
      new_fragment_uri);
  if (!st.ok()) {
    storage_manager_->array_close_for_reads(array_uri);
    storage_manager_->array_close_for_writes(array_uri);
    clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);
    return st;
  }

  // Read from one array and write to the other
  st = copy_array(query_r, query_w);
  if (!st.ok()) {
    storage_manager_->array_close_for_reads(array_uri);
    storage_manager_->array_close_for_writes(array_uri);
    clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);
    return st;
  }

  // Close array for reading
  st = storage_manager_->array_close_for_reads(array_uri);
  if (!st.ok()) {
    storage_manager_->array_close_for_writes(array_uri);
    storage_manager_->vfs()->remove_dir(*new_fragment_uri);
    clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);
    return st;
  }

  // Lock the array exclusively
  st = storage_manager_->array_xlock(array_uri);
  if (!st.ok()) {
    storage_manager_->array_close_for_writes(array_uri);
    storage_manager_->vfs()->remove_dir(*new_fragment_uri);
    clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);
    return st;
  }

  // Finalize both queries
  st = query_w->finalize();
  if (!st.ok()) {
    storage_manager_->array_close_for_writes(array_uri);
    clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);
    storage_manager_->array_xunlock(array_uri);
    bool is_dir;
    auto st2 = storage_manager_->vfs()->is_dir(*new_fragment_uri, &is_dir);
    if (is_dir)
      storage_manager_->vfs()->remove_dir(*new_fragment_uri);
    return st;
  }

  // Close array
  storage_manager_->array_close_for_writes(array_uri);
  if (!st.ok()) {
    storage_manager_->array_xunlock(array_uri);
    clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);
    bool is_dir;
    auto st2 = storage_manager_->vfs()->is_dir(*new_fragment_uri, &is_dir);
    if (is_dir)
      storage_manager_->vfs()->remove_dir(*new_fragment_uri);
    return st;
  }

  // Delete old fragment metadata. This makes the old fragments invisible
  st = delete_old_fragment_metadata(to_consolidate);
  if (!st.ok()) {
    delete_old_fragments(to_consolidate);
    storage_manager_->array_xunlock(array_uri);
    clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);
    return st;
  }

  // Unlock the array
  st = storage_manager_->array_xunlock(array_uri);
  if (!st.ok()) {
    delete_old_fragments(to_consolidate);
    clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);
    return st;
  }

  // Delete old fragments. The array does not need to be locked.
  st = delete_old_fragments(to_consolidate);

  // Clean up
  clean_up(subarray, buffer_num, buffers, buffer_sizes, query_r, query_w);

  return st;
}

Status Consolidator::copy_array(Query* query_r, Query* query_w) {
  do {
    RETURN_NOT_OK(query_r->submit());
    RETURN_NOT_OK(query_w->submit());
  } while (query_r->status() == QueryStatus::INCOMPLETE);

  return Status::Ok();
}

void Consolidator::clean_up(
    void* subarray,
    unsigned buffer_num,
    void** buffers,
    uint64_t* buffer_sizes,
    Query* query_r,
    Query* query_w) const {
  std::free(subarray);
  free_buffers(buffer_num, buffers, buffer_sizes);
  delete query_r;
  delete query_w;
}

Status Consolidator::create_buffers(
    const ArraySchema* array_schema,
    bool sparse_mode,
    void*** buffers,
    uint64_t** buffer_sizes,
    unsigned int* buffer_num) {
  // For easy reference
  auto attribute_num = array_schema->attribute_num();
  auto sparse = !array_schema->dense() || sparse_mode;

  // Calculate number of buffers
  *buffer_num = 0;
  for (unsigned int i = 0; i < attribute_num; ++i)
    *buffer_num += (array_schema->attributes()[i]->var_size()) ? 2 : 1;
  *buffer_num += (sparse) ? 1 : 0;

  // Create buffers
  *buffers = (void**)std::malloc(*buffer_num * sizeof(void*));
  if (*buffers == nullptr) {
    return LOG_STATUS(Status::ConsolidationError(
        "Cannot create consolidation buffers; Memory allocation failed"));
  }
  *buffer_sizes = new uint64_t[*buffer_num];
  if (*buffer_sizes == nullptr) {
    return LOG_STATUS(Status::ConsolidationError(
        "Cannot create consolidation buffer sizes; Memory allocation failed"));
  }

  // Allocate space for each buffer
  bool error = false;
  for (unsigned int i = 0; i < *buffer_num; ++i) {
    (*buffers)[i] = std::malloc(config_.buffer_size_);
    if ((*buffers)[i] == nullptr)  // The loop should continue to
      error = true;                // allocate nullptr to each buffer
    (*buffer_sizes)[i] = config_.buffer_size_;
  }

  // Clean up upon error
  if (error) {
    free_buffers(*buffer_num, *buffers, *buffer_sizes);
    *buffers = nullptr;
    *buffer_sizes = nullptr;
    return LOG_STATUS(Status::ConsolidationError(
        "Cannot create consolidation buffers; Memory allocation failed"));
  }

  // Success
  return Status::Ok();
}

Status Consolidator::create_queries(
    Array* array_for_reads,
    Array* array_for_writes,
    bool sparse_mode,
    void* subarray,
    Layout layout,
    void** buffers,
    uint64_t* buffer_sizes,
    Query** query_r,
    Query** query_w,
    URI* new_fragment_uri) {
  // Create read query
  *query_r = new Query(storage_manager_, array_for_reads);
  if (!(*query_r)->array_schema()->is_kv())
    RETURN_NOT_OK((*query_r)->set_layout(layout));
  RETURN_NOT_OK(set_query_buffers(*query_r, sparse_mode, buffers, buffer_sizes));
  RETURN_NOT_OK((*query_r)->set_subarray(subarray));
  if(sparse_mode)
    RETURN_NOT_OK((*query_r)->set_sparse_mode());
  // TODO: in Reader, check read_all_tiles and filter_all_tiles

  // Get last fragment URI, which will be the URI of the consolidated fragment
  *new_fragment_uri = (*query_r)->last_fragment_uri();
  RETURN_NOT_OK(rename_new_fragment_uri(new_fragment_uri));

  // Create write query
  *query_w = new Query(storage_manager_, array_for_writes, *new_fragment_uri);
  if (!(*query_r)->array_schema()->is_kv())
    RETURN_NOT_OK((*query_w)->set_layout(layout));
  RETURN_NOT_OK((*query_w)->set_subarray(subarray));
  RETURN_NOT_OK(set_query_buffers(*query_w, sparse_mode, buffers, buffer_sizes));

  return Status::Ok();
}

Status Consolidator::compute_subarray_and_layout(
    Array* array,
    const std::vector<FragmentInfo>& to_consolidate,
    void** subarray,
    Layout* layout,
    bool* all_sparse) const {
  // Check if all fragments to consolidate are sparse
  *all_sparse = true;
  for(const auto& f : to_consolidate) {
    if (!f.sparse_) {
      *all_sparse = false;
      break;
    }
  }

  auto array_schema = array->array_schema();
  assert(array_schema != nullptr);
  assert(*subarray == nullptr);

  // Create subarray only for the dense case
  if (!all_sparse) {
    *subarray = std::malloc(2 * array_schema->coords_size());
    if (*subarray == nullptr)
      return LOG_STATUS(Status::ConsolidationError(
          "Cannot create subarray; Failed to allocate memory"));

    bool coincides = false;
    compute_non_empty_domain(
        array_schema, to_consolidate, *subarray, &coincides);
    *layout = coincides ? Layout::GLOBAL_ORDER : array_schema->cell_order();
  } else {  // Sparse array
    *subarray = nullptr;
    *layout = Layout::GLOBAL_ORDER;
  }

  return Status::Ok();
}

Status Consolidator::compute_non_empty_domain(
    ArraySchema* array_schema,
    const std::vector<FragmentInfo>& to_consolidate,
    void* non_empty_domain,
    bool* coincides) const {
  // Compute buffer sizes
  switch (array_schema->coords_type()) {
    case Datatype::INT32:
      return compute_non_empty_domain<int>(
          array_schema,
          to_consolidate,
          static_cast<int*>(non_empty_domain),
          coincides);
    case Datatype::INT64:
      return compute_non_empty_domain<int64_t>(
          array_schema,
          to_consolidate,
          static_cast<int64_t*>(non_empty_domain),
          coincides);
    case Datatype::INT8:
      return compute_non_empty_domain<int8_t>(
          array_schema,
          to_consolidate,
          static_cast<int8_t*>(non_empty_domain),
          coincides);
    case Datatype::UINT8:
      return compute_non_empty_domain<uint8_t>(
          array_schema,
          to_consolidate,
          static_cast<uint8_t*>(non_empty_domain),
          coincides);
    case Datatype::INT16:
      return compute_non_empty_domain<int16_t>(
          array_schema,
          to_consolidate,
          static_cast<int16_t*>(non_empty_domain),
          coincides);
    case Datatype::UINT16:
      return compute_non_empty_domain<uint16_t>(
          array_schema,
          to_consolidate,
          static_cast<uint16_t*>(non_empty_domain),
          coincides);
    case Datatype::UINT32:
      return compute_non_empty_domain<uint32_t>(
          array_schema,
          to_consolidate,
          static_cast<uint32_t*>(non_empty_domain),
          coincides);
    case Datatype::UINT64:
      return compute_non_empty_domain<uint64_t>(
          array_schema,
          to_consolidate,
          static_cast<uint64_t*>(non_empty_domain),
          coincides);
    default:
      return LOG_STATUS(
          Status::StorageManagerError("Cannot compute non-empty domain "
                                      "sizes; Invalid domain type"));
  }

  return Status::Ok();
}

template <class T>
Status Consolidator::compute_non_empty_domain(
    ArraySchema* array_schema,
    const std::vector<FragmentInfo>& to_consolidate,
    T* non_empty_domain,
    bool* coincides) const {
  auto dim_num = array_schema->dim_num();
  auto domain_size = 2 * array_schema->coords_size();
  auto domain = (T*)array_schema->domain()->domain();
  auto tile_extents = (T*)array_schema->domain()->tile_extents();
  assert(!to_consolidate.empty());

  // Initialize subarray to the first fragment non-empty domain
  std::memcpy(
      non_empty_domain, to_consolidate[0].non_empty_domain_, domain_size);

  // Initialize subarray
  for (size_t i = 1; i < to_consolidate.size(); ++i) {
    utils::geometry::expand_mbr_with_mbr<T>(
        non_empty_domain,
        (const T*)to_consolidate[i].non_empty_domain_,
        dim_num);
  }

  *coincides = true;
  for (unsigned i = 0; i < dim_num; ++i) {
    assert(
        (non_empty_domain[2 * i + 1] - domain[2 * i]) !=
        std::numeric_limits<T>::max());
    if (((non_empty_domain[2 * i] - domain[2 * i]) % tile_extents[i]) ||
        (non_empty_domain[2 * i + 1] - domain[2 * i] + 1) % tile_extents[i]) {
      *coincides = false;
      break;
    }
  }

  return Status::Ok();
}

Status Consolidator::delete_old_fragment_metadata(
    const std::vector<FragmentInfo>& fragments) {
  for (auto& fragment : fragments) {
    auto meta_uri =
        fragment.uri_.join_path(constants::fragment_metadata_filename);
    RETURN_NOT_OK(storage_manager_->vfs()->remove_file(meta_uri));
  }

  return Status::Ok();
}

Status Consolidator::delete_old_fragments(
    const std::vector<FragmentInfo>& fragments) {
  for (auto& fragment : fragments)
    RETURN_NOT_OK(storage_manager_->vfs()->remove_dir(fragment.uri_));

  return Status::Ok();
}

void Consolidator::free_buffers(
    unsigned int buffer_num, void** buffers, uint64_t* buffer_sizes) const {
  for (unsigned int i = 0; i < buffer_num; ++i) {
    std::free(buffers[i]);
  }
  std::free(buffers);
  delete[] buffer_sizes;
}

Status Consolidator::compute_next_to_consolidate(
    const std::vector<FragmentInfo>& fragments,
    std::vector<FragmentInfo>* to_consolidate) const {
  // Preparation
  to_consolidate->clear();
  auto min = config_.min_frags_;
  min = (uint32_t)((min > fragments.size()) ? fragments.size() : min);
  auto max = config_.max_frags_;
  max = (uint32_t)((max > fragments.size()) ? fragments.size() : max);
  auto size_ratio = config_.size_ratio_;

  // Trivial case - all fragments
  if (min == max && max == fragments.size() && size_ratio == 0.0f) {
    *to_consolidate = fragments;
    return Status::Ok();
  }

  // Trivial case - no fragments
  if (max == 0)
    return Status::Ok();

  // Prepare the dynamic-programming matrix. The rows are from 1 to max
  // and the columns represent the fragments in `fragments`.
  std::vector<std::vector<uint64_t>> m;
  auto col_num = fragments.size();
  auto row_num = max;
  m.resize(row_num);
  for (auto& row : m)
    row.resize(col_num);

  // TODO remove
  for (auto f : fragments) {
    std::cout << f.uri_.to_string() << ": " << f.fragment_size_ << " "
              << ((f.sparse_) ? "sparse" : "dense") << " " << f.timestamp_
              << "\n";
    auto d = (uint64_t*)f.non_empty_domain_;
    std::cout << d[0] << " " << d[1] << "\n";
  }

  // Entry m[i][j] contains the collective size of fragments
  // fragments[j], ..., fragments[j+i]. If the size ratio
  // of any adjacent pair in the above list is smaller than the
  // defined one, then the size sum of that entry is infinity (UINT64_MAX).
  for (size_t i = 0; i < row_num; ++i) {
    for (size_t j = 0; j < col_num; ++j) {
      if (i == 0) {  // In the first row we store the sizes of `fragments`
        m[i][j] = fragments[j].fragment_size_;
      } else if (i + j >= col_num) {  // Non-valid entries
        m[i][j] = UINT64_MAX;
      } else {  // Every other row is computed using the previous row
        auto ratio = (float)fragments[i + (j - 1)].fragment_size_ /
                     fragments[i + j].fragment_size_;
        ratio = (ratio <= 1.0f) ? ratio : 1.0f / ratio;
        if (ratio >= size_ratio)
          m[i][j] = m[i - 1][j] + fragments[i + j].fragment_size_;
        else
          m[i][j] = UINT64_MAX;
      }
      // TODO remove
      std::cout << "m[" << i << "][" << j << "] = " << m[i][j] << "\n";
    }
  }

  // Choose the maximal set of fragments with cardinality in [min, max]
  // with the minimum size
  uint64_t min_size = UINT64_MAX;
  size_t min_col = 0;
  for (int i = row_num - 1; (i >= 0 && i >= (int)min - 1); --i) {
    min_size = UINT64_MAX;
    for (size_t j = 0; j < col_num; ++j) {
      if (m[i][j] < min_size) {
        min_size = m[i][j];
        min_col = j;
      }
    }

    // TODO remove
    std::cout << "min_col: " << min_col << "\n";
    std::cout << "min_size: " << min_size << "\n";

    // Results found
    if (min_size != UINT64_MAX) {
      for (size_t f = min_col; f <= min_col + i; ++f)
        to_consolidate->emplace_back(fragments[f]);
      break;
    }
  }

  // TODO remove
  for (const auto& f : *to_consolidate)
    std::cout << f.uri_.to_string() << " " << f.fragment_size_ << " "
              << f.timestamp_ << "\n";

  return Status::Ok();
}

Status Consolidator::rename_new_fragment_uri(URI* uri) const {
  // Get timestamp
  std::string name = uri->last_path_part();
  auto timestamp_str = name.substr(name.find_last_of('_') + 1);

  // Get current time
  uint64_t ms = utils::time::timestamp_now_ms();

  // Get uuid
  std::string uuid;
  RETURN_NOT_OK(uuid::generate_uuid(&uuid, false));

  std::stringstream ss;
  ss << uri->parent().to_string() << "/__" << uuid << "_" << ms << "_"
     << timestamp_str;

  *uri = URI(ss.str());
  return Status::Ok();
}

Status Consolidator::set_query_buffers(
    Query* query, bool sparse_mode, void** buffers, uint64_t* buffer_sizes) const {
  auto dense = query->array_schema()->dense();
  auto attributes = query->array_schema()->attributes();
  unsigned bid = 0;
  for (const auto& attr : attributes) {
    if (!attr->var_size()) {
      RETURN_NOT_OK(
          query->set_buffer(attr->name(), buffers[bid], &buffer_sizes[bid]));
      ++bid;
    } else {
      RETURN_NOT_OK(query->set_buffer(
          attr->name(),
          (uint64_t*)buffers[bid],
          &buffer_sizes[bid],
          buffers[bid + 1],
          &buffer_sizes[bid + 1]));
      bid += 2;
    }
  }
  if (!dense || sparse_mode)
    RETURN_NOT_OK(
        query->set_buffer(constants::coords, buffers[bid], &buffer_sizes[bid]));

  return Status::Ok();
}

void Consolidator::update_fragment_info(
    const std::vector<FragmentInfo>& to_consolidate,
    const FragmentInfo& new_fragment_info,
    std::vector<FragmentInfo>* fragment_info) const {
  auto to_consolidate_it = to_consolidate.begin();
  auto fragment_it = fragment_info->begin();
  std::vector<FragmentInfo> updated_fragment_info;
  bool new_fragment_added = false;

  while (fragment_it != fragment_info->end()) {
    // No match - add the fragment info and advance `fragment_it`
    if (to_consolidate_it == to_consolidate.end() ||
        fragment_it->uri_.to_string() != to_consolidate_it->uri_.to_string()) {
      updated_fragment_info.emplace_back(*fragment_it);
      ++fragment_it;
    } else {  // Match - add new fragment only once and advance both iterators
      if (!new_fragment_added) {
        updated_fragment_info.emplace_back(new_fragment_info);
        new_fragment_added = true;
      }
      ++fragment_it;
      ++to_consolidate_it;
    }
  }

  assert(
      updated_fragment_info.size() ==
      fragment_info->size() - to_consolidate.size() + 1);

  *fragment_info = std::move(updated_fragment_info);
}

Status Consolidator::set_config(const Config* config) {
  if (config != nullptr) {
    auto params = config->consolidation_params();
    config_.steps_ = params.steps_;
    config_.buffer_size_ = params.buffer_size_;
    config_.size_ratio_ = params.step_size_ratio_;
    config_.min_frags_ = params.step_min_frags_;
    config_.max_frags_ = params.step_max_frags_;
  }

  // Sanity checks
  if (config_.min_frags_ > config_.max_frags_)
    return LOG_STATUS(Status::ConsolidationError(
        "Invalid configuration; Minimum fragments config parameter is larger "
        "than the maximum"));
  if (config_.size_ratio_ > 1.0f || config_.size_ratio_ < 0.0f)
    return LOG_STATUS(Status::ConsolidationError(
        "Invalid configuration; Step size ratio config parameter must be in "
        "[0.0, 1.0]"));

  return Status::Ok();
}

}  // namespace sm
}  // namespace tiledb
