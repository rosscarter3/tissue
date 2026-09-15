//
// Flat CSR-backed row table used for all simulation state (cell, wall and
// vertex data and derivatives). Replaces the legacy vector<vector<double>>
// DataMatrix: one contiguous allocation (cache friendly, autovectorizable
// solver loops over flat()), while still supporting the legacy semantics the
// simulator relies on:
//   - ragged rows (CenterTriangulation reactions resize individual cell rows)
//   - append-row on cell division
//   - swap-with-last row removal on cell/wall/vertex removal
// Row access reads like the old code: data[i][j].
//
#ifndef TISSUE2_CORE_MATRIX_H
#define TISSUE2_CORE_MATRIX_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

namespace tissue {

class Matrix {
public:
  Matrix() { offsets_.push_back(0); }
  Matrix(size_t rows, size_t cols, double value = 0.0) { assign(rows, cols, value); }

  size_t rows() const { return offsets_.size() - 1; }
  bool empty() const { return rows() == 0; }

  size_t rowSize(size_t r) const { return offsets_[r + 1] - offsets_[r]; }
  // Index of row r's first element within flat(); lets callers precompute
  // flat indices for scattered columns.
  size_t rowBegin(size_t r) const { return offsets_[r]; }

  // Row length of row 0; only meaningful for row-uniform tables (vertex and
  // wall data always are; cell data is unless a CT reaction resized rows).
  size_t cols() const { return empty() ? 0 : rowSize(0); }

  std::span<double> operator[](size_t r) {
    assert(r < rows());
    return {data_.data() + offsets_[r], rowSize(r)};
  }
  std::span<const double> operator[](size_t r) const {
    assert(r < rows());
    return {data_.data() + offsets_[r], rowSize(r)};
  }

  // Entire storage as one flat range (solver axpy/norm loops).
  std::span<double> flat() { return {data_.data(), data_.size()}; }
  std::span<const double> flat() const { return {data_.data(), data_.size()}; }
  size_t size() const { return data_.size(); }

  void assign(size_t rows, size_t cols, double value = 0.0) {
    data_.assign(rows * cols, value);
    offsets_.resize(rows + 1);
    for (size_t r = 0; r <= rows; ++r)
      offsets_[r] = r * cols;
  }

  void clear() {
    data_.clear();
    offsets_.assign(1, 0);
  }

  void fill(double value) { std::fill(data_.begin(), data_.end(), value); }

  // Appends a row with the given content; returns its index.
  size_t appendRow(std::span<const double> row) {
    data_.insert(data_.end(), row.begin(), row.end());
    offsets_.push_back(data_.size());
    return rows() - 1;
  }
  // Appends a zero row of length n.
  size_t appendRow(size_t n) {
    data_.insert(data_.end(), n, 0.0);
    offsets_.push_back(data_.size());
    return rows() - 1;
  }
  // Appends a copy of row src (legacy resize(n+1, data[src]) idiom).
  size_t appendRowCopy(size_t src) {
    // Note: spans into data_ would dangle across the reallocation, so copy
    // via index range after reserving.
    size_t b = offsets_[src], e = offsets_[src + 1];
    data_.reserve(data_.size() + (e - b));
    data_.insert(data_.end(), data_.begin() + b, data_.begin() + e);
    offsets_.push_back(data_.size());
    return rows() - 1;
  }

  // Legacy removal idiom: data[r] = data.back(); data.pop_back().
  // Row r takes the last row's content (and length); all other rows keep
  // their order.
  void removeRowSwap(size_t r) {
    size_t last = rows() - 1;
    assert(r <= last);
    if (r != last) {
      std::vector<double> lastRow(data_.begin() + offsets_[last], data_.end());
      // Erase row r, shifting everything after it left.
      data_.erase(data_.begin() + offsets_[r], data_.begin() + offsets_[r + 1]);
      size_t removedLen = offsets_[r + 1] - offsets_[r];
      for (size_t i = r + 1; i < offsets_.size(); ++i)
        offsets_[i] -= removedLen;
      // The old last row now sits at the end minus its own length; drop it
      // and re-insert its content at position r.
      data_.resize(offsets_[last] - 0);
      data_.insert(data_.begin() + offsets_[r], lastRow.begin(), lastRow.end());
      for (size_t i = r + 1; i < last; ++i)
        offsets_[i] += lastRow.size();
      offsets_[last] = data_.size();
      offsets_.pop_back();
      // Recompute trailing offset consistency in debug builds.
      assert(offsets_.back() == data_.size());
    } else {
      data_.resize(offsets_[last]);
      offsets_.pop_back();
    }
  }

  // Resizes row r to newLen (used by CenterTriangulation initiation), zero
  // filling any new elements.
  void resizeRow(size_t r, size_t newLen) {
    size_t oldLen = rowSize(r);
    if (newLen == oldLen)
      return;
    if (newLen > oldLen) {
      data_.insert(data_.begin() + offsets_[r + 1], newLen - oldLen, 0.0);
      for (size_t i = r + 1; i < offsets_.size(); ++i)
        offsets_[i] += newLen - oldLen;
    } else {
      data_.erase(data_.begin() + offsets_[r] + newLen,
                  data_.begin() + offsets_[r + 1]);
      for (size_t i = r + 1; i < offsets_.size(); ++i)
        offsets_[i] -= oldLen - newLen;
    }
  }

  // Copies content from another table with identical shape.
  void copyFrom(const Matrix &other) {
    assert(offsets_ == other.offsets_);
    std::copy(other.data_.begin(), other.data_.end(), data_.begin());
  }

  // Adopts the shape (row structure) of another table, zero-initialized.
  void reshapeLike(const Matrix &other) {
    offsets_ = other.offsets_;
    data_.assign(other.data_.size(), 0.0);
  }

  bool sameShape(const Matrix &other) const { return offsets_ == other.offsets_; }

private:
  std::vector<double> data_;
  std::vector<size_t> offsets_; // rows()+1 entries; offsets_[r]..offsets_[r+1]
};

} // namespace tissue

#endif
