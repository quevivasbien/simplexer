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
    AbstractMatrix(const std::vector<T>& data, int rows, int cols) : data(data), rows(rows), cols(cols) {
        assert(data.size() == rows * cols);
    }

    explicit AbstractMatrix(const std::vector<std::vector<T>>& data) : rows(data.size()), cols(data[0].size()) {
        assert(data.size() > 0);
        this->data = std::vector<T>();
        for (int i = 0; i < this->rows; i++) {
            assert(data[i].size() == this->cols);
            for (int j = 0; j < this->cols; j++) {
                this->data.push_back(data[i][j]);
            }
        }
    }

    T get(int row, int col) const {
        assert(this->inbounds(row, col));
        return data[row * this->cols + col];
    }

    void set(int row, int col, T value) {
        assert(this->inbounds(row, col));
        data[row * this->cols + col] = value;
    }

    bool inbounds(int row, int col) const {
        return this->rowInbounds(row) && this->colInbounds(col);
    }

    void mapInplace(std::function<T (T, int)> fn) {
        for (int i = 0; i < this->rows * this->cols; i++) {
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

    bool rowInbounds(int row) const {
        return row >= 0 && row < this->rows;
    }
    bool colInbounds(int col) const {
        return col >= 0 && col < this->cols;
    }

    const int rows;
    const int cols;

protected:
    std::vector<T> data;

};

template <typename T>
requires Number<T>
class Matrix : public AbstractMatrix<T> {
public:
    using AbstractMatrix<T>::AbstractMatrix;

    static Matrix<T> identity(int n) {
        std::vector<T> data;
        data.reserve(n * n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
               data.push_back(i == j ? 1 : 0);
            }
        }
        return Matrix<T>(data, n, n);
    }

    static Matrix<T> full(T value, int rows, int cols) {
        std::vector<T> data;
        data.reserve(rows * cols);
        for (int i = 0; i < rows * cols; i++) {
            data.push_back(value);
        }
        return Matrix<T>(data, rows, cols);
    }

    static Matrix<T> fromScalar(T value) {
        return Matrix({value}, 1, 1);
    }

    static Matrix<T> zeros(int rows, int cols) {
        return Matrix<T>::full(0, rows, cols);
    }

    virtual std::string str() const override {
        std::stringstream ss;
        ss << "[";
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
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
        std::vector<T> result;
        result.reserve(this->rows * this->cols);
        for (int i = 0; i < this->rows * this->cols; i++) {
            result.push_back(fn(this->data[i], i));
        }
        return Matrix<T>(result, this->rows, this->cols);
    }

    // overload arithmetic operators
    Matrix<T> operator+(const Matrix<T>& other) const {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, int i) { return x + other.data[i]; };
        return this->map(fn);
    }
    void operator+=(const Matrix<T>& other) {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, int i) { return x + other.data[i]; };
        this->mapInplace(fn);
    }
    Matrix<T> operator-(const Matrix<T>& other) const {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, int i) { return x - other.data[i]; };
        return this->map(fn);
    }
    void operator-=(const Matrix<T>& other) {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, int i) { return x - other.data[i]; };
        this->mapInplace(fn);
    }
    Matrix<T> operator*(const Matrix<T>& other) const {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, int i) { return x * other.data[i]; };
        return this->map(fn);
    }
    void operator*=(const Matrix<T>& other) {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, int i) { return x * other.data[i]; };
        this->mapInplace(fn);
    }
    Matrix<T> operator/(const Matrix<T>& other) const {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, int i) { return x / other.data[i]; };
        return this->map(fn);
    }
    void operator/=(const Matrix<T>& other) {
        assert(this->dimsEqual(other));
        auto fn = [other](T x, int i) { return x / other.data[i]; };
        this->mapInplace(fn);
    }

    Matrix<T> operator+(T scalar) const {
        auto fn = [scalar](T x, int _i) { return x + scalar; };
        return this->map(fn);
    }
    void operator+=(T scalar) {
        auto fn = [scalar](T x, int _i) { return x + scalar; };
        this->mapInplace(fn);
    }
    Matrix<T> operator-(T scalar) const {
        auto fn = [scalar](T x, int _i) { return x - scalar; };
        return this->map(fn);
    }
    void operator-=(T scalar) {
        auto fn = [scalar](T x, int _i) { return x - scalar; };
        this->mapInplace(fn);
    }
    Matrix<T> operator*(T scalar) const {
        auto fn = [scalar](T x, int _i) { return x * scalar; };
        return this->map(fn);
    }
    void operator*=(T scalar) {
        auto fn = [scalar](T x, int _i) { return x * scalar; };
        this->mapInplace(fn);
    }
    Matrix<T> operator/(T scalar) const {
        auto fn = [scalar](T x, int _i) { return x / scalar; };
        return this->map(fn);
    }
    void operator/=(T scalar) {
        auto fn = [scalar](T x, int _i) { return x / scalar; };
        this->mapInplace(fn);
    }

    Matrix<T> transpose() const {
        std::vector<T> result;
        result.reserve(this->rows * this->cols);
        for (int i = 0; i < this->cols; i++) {
            for (int j = 0; j < this->rows; j++) {
                result.push_back(this->get(j, i));
            }
        }
        return Matrix<T>(result, this->cols, this->rows);
    }

    // methods for row and column operations
    void mapRowInplace(int row, std::function<T (T, int)> fn) {
        assert(this->rowInbounds(row));
        for (int i = 0; i < this->cols; i++) {
            this->data[row * this->cols + i] = fn(this->data[row * this->cols + i], i);
        }
    }
    void swapRowsInplace(int row1, int row2) {
        assert(this->rowInbounds(row1));
        assert(this->rowInbounds(row2));
        // store row1 in temp vector
        std::vector<T> temp;
        temp.reserve(this->cols);
        for (int i = 0; i < this->cols; i++) {
            temp.push_back(this->data[row1 * this->cols + i]);
        }
        // move row2 into row1
        for (int i = 0; i < this->cols; i++) {
            this->data[row1 * this->cols + i] = this->data[row2 * this->cols + i];
        }
        // move temp into row2
        for (int i = 0; i < this->cols; i++) {
            this->data[row2 * this->cols + i] = temp[i];
        }
    }
    void addRowInplace(int rowInto, int rowFrom, T multiple = 1) {
        // add rowFrom * multiple to rowInto
        assert(this->rowInbounds(rowFrom));
        assert(this->rowInbounds(rowInto));
        for (int i = 0; i < this->cols; i++) {
            this->data[rowInto * this->cols + i] += this->data[rowFrom * this->cols + i] * multiple;
        }
    }
    void scaleRowInplace(int row, T scalar) {
        assert(this->rowInbounds(row));
        for (int i = 0; i < this->cols; i++) {
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
    BlockMatrix(const std::vector<Matrix<T>>& data, int rows, int cols) : AbstractMatrix<Matrix<T>>(data, rows, cols) {
        assert(this->dimsOk());
    };

    explicit BlockMatrix(const std::vector<std::vector<Matrix<T>>>& data) : AbstractMatrix<Matrix<T>>(data) {
        assert(this->dimsOk());
    }; 

    virtual std::string str() const override {
        std::stringstream ss;
        ss << "[";
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
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
        int rows = 0;
        int cols = 0;
        for (int i = 0; i < this->rows; i++) {
            rows += this->get(i, 0).rows;
        }
        for (int i = 0; i < this->cols; i++) {
            cols += this->get(0, i).cols;
        }
        std::vector<T> result;
        result.reserve(rows * cols);
        for (int i = 0; i < this->rows; i++) {
            for (int k = 0; k < this->get(i, 0).rows; k++) {
                for (int j = 0; j < this->cols; j++) {
                    for (int l = 0; l < this->get(0, j).cols; l++) {
                        result.push_back(this->get(i, j).get(k, l));
                    }
                }
            }
        }
        return Matrix<T>(result, rows, cols);
    }
private:
    bool dimsOk() const {
        for (int i = 0; i < this->rows; i++) {
            // check that matrices in the same row have the same number of rows
            for (int j = 0; j < this->cols - 1; j++) {
                if (this->get(i, j).rows != this->get(i, j+1).rows) return false;
            }
        }
        for (int j = 0; j < this->cols; j++) {
            // check that matrices in the same column have the same number of columns
            for (int i = 0; i < this->rows - 1; i++) {
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

    virtual T get(int idx) const = 0;
    virtual int size() const = 0;

    Vector<T> vector() const {
        std::vector<T> result;
        result.reserve(this->size());
        for (int i = 0; i < this->size(); i++) {
            result.push_back(this->get(i));
        }
        return Vector<T>(result);
    }

    bool isZero() const {
        for (int i = 0; i < this->size(); i++) {
            if (this->get(i) != 0) return false;
        }
        return true;
    }

    bool isElementaryBasis() const {
        // returns true iff only one element is non-zero and it is 1
        bool foundNonZero = false;
        for (int i = 0; i < this->size(); i++) {
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
    RowView(const M<T>& data, int row) : DimView<T, M>(data), row(row) {
        assert(data.rowInbounds(row));
    }

    virtual T get(int col) const override {
        return this->data.get(this->row, col);
    }

    virtual int size() const override {
        return this->data.cols;
    }

    const int row;
};

template <typename T, template<typename> typename M = Matrix>
class ColView : public DimView<T, M> {
public:
    ColView(const M<T>& data, int col) : DimView<T, M>(data), col(col) {
        assert(data.colInbounds(col));
    }

    virtual T get(int row) const override {
        return this->data.get(row, this->col);
    }

    virtual int size() const override {
        return this->data.rows;
    }

    const int col;
};