#pragma once

#include <vector>
#include <functional>
#include <cassert>
#include <sstream>
#include <string>
#include <concepts>

template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;

template <typename T>
class AbstractMatrix {
public:
    AbstractMatrix(const std::vector<T>& data, std::size_t rows, std::size_t cols) : data(data), rows(rows), cols(cols) {
        assert(data.size() == rows * cols);
    }

    explicit AbstractMatrix(const std::vector<std::vector<T>>& data) : rows(data.size()), cols(data[0].size()) {
        assert(data.size() > 0);
        this->data = std::vector<T>();
        for (const std::vector<T>& row : data) {
            assert(row.size() == this->cols);
            for (const T& val : row) {
                this->data.push_back(val);
            }
        }
    }

    T get(std::size_t row, std::size_t col) const {
        assert(this->inbounds(row, col));
        return data[row * this->cols + col];
    }

    void set(std::size_t row, std::size_t col, T value) {
        assert(this->inbounds(row, col));
        data[row * this->cols + col] = value;
    }

    bool inbounds(std::size_t row, std::size_t col) const {
        return this->rowInbounds(row) && this->colInbounds(col);
    }

    void mapInplace(std::function<T (T, std::size_t)> fn) {
        for (std::size_t i = 0; i < this->rows * this->cols; i++) {
            this->data[i] = fn(this->data[i], i);
        }
    }

    bool dimsEqual(const AbstractMatrix<T>& other) const {
        return this->rows == other.rows && this->cols == other.cols;
    }

    const std::vector<T>& dataView() const {
        return this->data;
    }

    virtual std::string str() const = 0;

    bool rowInbounds(std::size_t row) const {
        return row < this->rows;
    }
    bool colInbounds(std::size_t col) const {
        return col < this->cols;
    }

    const std::size_t rows;
    const std::size_t cols;

protected:
    std::vector<T> data;

};

template <typename T>
requires Number<T>
class Matrix : public AbstractMatrix<T> {
public:
    using AbstractMatrix<T>::AbstractMatrix;

    static Matrix<T> identity(std::size_t n) {
        std::vector<T> data;
        data.reserve(n * n);
        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 0; j < n; j++) {
               data.push_back(i == j ? 1 : 0);
            }
        }
        return Matrix<T>(data, n, n);
    }

    static Matrix<T> full(T value, std::size_t rows, std::size_t cols) {
        return Matrix<T>(std::vector<T>(rows * cols, value), rows, cols);
    }

    static Matrix<T> fromScalar(T value) {
        return Matrix({value}, 1, 1);
    }

    static Matrix<T> zeros(std::size_t rows, std::size_t cols) {
        return Matrix<T>::full(0, rows, cols);
    }

    virtual std::string str() const override {
        std::stringstream ss;
        ss << "[";
        for (std::size_t i = 0; i < this->rows; i++) {
            for (std::size_t j = 0; j < this->cols; j++) {
                ss << this->get(i, j);
                if (j < this->cols - 1) {
                    ss << ", ";
                }
            }
            if (i < this->rows - 1) {
                ss << "; ";
            }
        }
        ss << "]";
        return ss.str();
    }

    Matrix<T> map(std::function<T (T, int)> fn) const {
        std::vector<T> result(this->rows * this->cols);
        for (std::size_t i = 0; i < this->data.size(); i++) {
            result[i] = fn(this->data[i], i);
        }
        return Matrix<T>(result, this->rows, this->cols);
    }

    // overload arithmetic operators
    Matrix<T> operator+(const Matrix<T>& other) const {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, std::size_t i) { return x + other.data[i]; };
        return this->map(fn);
    }
    void operator+=(const Matrix<T>& other) {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, std::size_t i) { return x + other.data[i]; };
        this->mapInplace(fn);
    }
    Matrix<T> operator-(const Matrix<T>& other) const {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, std::size_t i) { return x - other.data[i]; };
        return this->map(fn);
    }
    void operator-=(const Matrix<T>& other) {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, std::size_t i) { return x - other.data[i]; };
        this->mapInplace(fn);
    }
    Matrix<T> operator*(const Matrix<T>& other) const {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, std::size_t i) { return x * other.data[i]; };
        return this->map(fn);
    }
    void operator*=(const Matrix<T>& other) {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, std::size_t i) { return x * other.data[i]; };
        this->mapInplace(fn);
    }
    Matrix<T> operator/(const Matrix<T>& other) const {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, std::size_t i) { return x / other.data[i]; };
        return this->map(fn);
    }
    void operator/=(const Matrix<T>& other) {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, std::size_t i) { return x / other.data[i]; };
        this->mapInplace(fn);
    }

    Matrix<T> operator+(T scalar) const {
        auto fn = [scalar](T x, std::size_t _i) { return x + scalar; };
        return this->map(fn);
    }
    void operator+=(T scalar) {
        auto fn = [scalar](T x, std::size_t _i) { return x + scalar; };
        this->mapInplace(fn);
    }
    Matrix<T> operator-(T scalar) const {
        auto fn = [scalar](T x, std::size_t _i) { return x - scalar; };
        return this->map(fn);
    }
    void operator-=(T scalar) {
        auto fn = [scalar](T x, std::size_t _i) { return x - scalar; };
        this->mapInplace(fn);
    }
    Matrix<T> operator*(T scalar) const {
        auto fn = [scalar](T x, std::size_t _i) { return x * scalar; };
        return this->map(fn);
    }
    void operator*=(T scalar) {
        auto fn = [scalar](T x, std::size_t _i) { return x * scalar; };
        this->mapInplace(fn);
    }
    Matrix<T> operator/(T scalar) const {
        auto fn = [scalar](T x, std::size_t _i) { return x / scalar; };
        return this->map(fn);
    }
    void operator/=(T scalar) {
        auto fn = [scalar](T x, std::size_t _i) { return x / scalar; };
        this->mapInplace(fn);
    }

    Matrix<T> transpose() const {
        std::vector<T> result;
        result.reserve(this->rows * this->cols);
        for (std::size_t i = 0; i < this->cols; i++) {
            for (std::size_t j = 0; j < this->rows; j++) {
                result.push_back(this->get(j, i));
            }
        }
        return Matrix<T>(result, this->cols, this->rows);
    }

    // methods for row and column operations
    void mapRowInplace(std::size_t row, std::function<T (T, int)> fn) {
        assert(this->rowInbounds(row));
        for (std::size_t i = 0; i < this->cols; i++) {
            this->data[row * this->cols + i] = fn(this->data[row * this->cols + i], i);
        }
    }
    void swapRowsInplace(std::size_t row1, std::size_t row2) {
        assert(this->rowInbounds(row1));
        assert(this->rowInbounds(row2));
        // store row1 in temp vector
        std::vector<T> temp;
        temp.reserve(this->cols);
        for (std::size_t i = 0; i < this->cols; i++) {
            temp.push_back(this->data[row1 * this->cols + i]);
        }
        // move row2 into row1
        for (std::size_t i = 0; i < this->cols; i++) {
            this->data[row1 * this->cols + i] = this->data[row2 * this->cols + i];
        }
        // move temp into row2
        for (std::size_t i = 0; i < this->cols; i++) {
            this->data[row2 * this->cols + i] = temp[i];
        }
    }
    void addRowInplace(std::size_t rowInto, std::size_t rowFrom, T multiple = 1) {
        // add rowFrom * multiple to rowInto
        assert(this->rowInbounds(rowFrom));
        assert(this->rowInbounds(rowInto));
        for (std::size_t i = 0; i < this->cols; i++) {
            this->data[rowInto * this->cols + i] += this->data[rowFrom * this->cols + i] * multiple;
        }
    }
    void scaleRowInplace(std::size_t row, T scalar) {
        assert(this->rowInbounds(row));
        for (std::size_t i = 0; i < this->cols; i++) {
            this->data[row * this->cols + i] *= scalar;
        }
    }

};


template <typename T>
requires Number<T>
class Vector : public Matrix<T> {
public:
    explicit Vector(const std::vector<T>& data) : Matrix<T>(data, data.size(), 1) {}
};

template <typename T>
class BlockMatrix : public AbstractMatrix<Matrix<T>> {
public:
    BlockMatrix(const std::vector<Matrix<T>>& data, std::size_t rows, std::size_t cols) : AbstractMatrix<Matrix<T>>(data, rows, cols) {
        assert(this->dimsOk());
    };

    explicit BlockMatrix(const std::vector<std::vector<Matrix<T>>>& data) : AbstractMatrix<Matrix<T>>(data) {
        assert(this->dimsOk());
    }; 

    virtual std::string str() const override {
        std::stringstream ss;
        ss << "[";
        for (std::size_t i = 0; i < this->rows; i++) {
            for (std::size_t j = 0; j < this->cols; j++) {
                ss << this->get(i, j).str();
                if (j < this->cols - 1) {
                    ss << ", ";
                }
            }
            if (i < this->rows - 1) {
                ss << "; ";
            }
        }
        ss << "]";
        return ss.str();
    }

    Matrix<T> matrix() const {
        std::size_t rows = 0;
        std::size_t cols = 0;
        for (std::size_t i = 0; i < this->rows; i++) {
            rows += this->get(i, 0).rows;
        }
        for (std::size_t i = 0; i < this->cols; i++) {
            cols += this->get(0, i).cols;
        }
        std::vector<T> result;
        result.reserve(rows * cols);
        for (std::size_t i = 0; i < this->rows; i++) {
            for (std::size_t k = 0; k < this->get(i, 0).rows; k++) {
                for (std::size_t j = 0; j < this->cols; j++) {
                    for (std::size_t l = 0; l < this->get(0, j).cols; l++) {
                        result.push_back(this->get(i, j).get(k, l));
                    }
                }
            }
        }
        return Matrix<T>(result, rows, cols);
    }
private:
    bool dimsOk() const {
        for (std::size_t i = 0; i < this->rows; i++) {
            // check that matrices in the same row have the same number of rows
            for (std::size_t j = 0; j < this->cols - 1; j++) {
                if (this->get(i, j).rows != this->get(i, j+1).rows) return false;
            }
        }
        for (std::size_t j = 0; j < this->cols; j++) {
            // check that matrices in the same column have the same number of columns
            for (std::size_t i = 0; i < this->rows - 1; i++) {
                if (this->get(i, j).cols != this->get(i+1, j).cols) return false;
            }
        }
        return true;
    }
};

template <typename T, template<typename> typename M>
class DimView {
public:
    DimView(const M<T>& data) : data(data) {}

    virtual T get(std::size_t idx) const = 0;
    virtual std::size_t size() const = 0;

    Vector<T> vector() const {
        std::vector<T> result;
        result.reserve(this->size());
        for (std::size_t i = 0; i < this->size(); i++) {
            result.push_back(this->get(i));
        }
        return Vector<T>(result);
    }

    bool isZero() const {
        for (std::size_t i = 0; i < this->size(); i++) {
            if (this->get(i) != 0) return false;
        }
        return true;
    }

    bool isElementaryBasis() const {
        // returns true iff only one element is non-zero and it is 1
        bool foundNonZero = false;
        for (std::size_t i = 0; i < this->size(); i++) {
            if (this->get(i) != 0) {
                if (foundNonZero || this->get(i) != 1) return false;
                foundNonZero = true;
            }
        }
    }
protected:
    const M<T>& data;
};

template <typename T, template<typename> typename M = Matrix>
class RowView : public DimView<T, M> {
public:
    RowView(const M<T>& data, std::size_t row) : DimView<T, M>(data), row(row) {
        assert(data.rowInbounds(row));
    }

    virtual T get(std::size_t col) const override {
        return this->data.get(this->row, col);
    }

    virtual std::size_t size() const override {
        return this->data.cols;
    }

    const std::size_t row;
};

template <typename T, template<typename> typename M = Matrix>
class ColView : public DimView<T, M> {
public:
    ColView(const M<T>& data, std::size_t col) : DimView<T, M>(data), col(col) {
        assert(data.colInbounds(col));
    }

    virtual T get(std::size_t row) const override {
        return this->data.get(row, this->col);
    }

    virtual std::size_t size() const override {
        return this->data.rows;
    }

    const std::size_t col;
};