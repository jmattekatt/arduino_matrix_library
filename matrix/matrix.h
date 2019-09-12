/*
    Arduino Matrix class
    @author: Joe Mattekatt
*/

#include "Arduino.h"

#ifndef MATRIX_H
#define MATRIX_H

#define DEFAULT_VALUE 0
#define START_INDEX   0
#define ASSERT(condition, message) if(!condition){Serial.println(message); Serial.flush(); abort();}


template <class T>
class Matrix {

    public:
        Matrix();
        Matrix(const unsigned int init_rows, const unsigned int init_cols);
        Matrix(const unsigned int init_rows, const unsigned int init_cols, const T fill_value);
        Matrix(const unsigned int init_rows, const unsigned int init_cols, const T* init_vals);
        ~Matrix();
        
        Matrix(const Matrix<T> &other);
        Matrix<T>& operator=(Matrix<T> other);
        Matrix(Matrix<T> &&other);

        void print() const;

        inline T& at(const unsigned int row, const unsigned int col);
        inline T at(const unsigned int row, const unsigned int col) const;

        template <typename U>
        decltype(auto) operator+(const U scalar) const;
        template <typename U>
        decltype(auto) operator-(const U scalar) const;
        template <typename U>
        decltype(auto) operator*(const U scalar) const;
        
        template <typename U>
        decltype(auto) operator+(const Matrix<U> &other) const;
        template <typename U>
        decltype(auto) operator-(const Matrix<U> &other) const;
        template <typename U>
        decltype(auto) operator*(const Matrix<U> &other) const;

        bool operator==(const Matrix<T> &other) const;
        Matrix<T> transpose() const;
        void LU_decompose(Matrix<float> &L, Matrix<float> &U) const;
        Matrix<float> solve_for(Matrix<T> &B) const;

        template <typename V>
        friend class Matrix;

    private:
        T* head;

        int rows;
        int cols;
        int size;

        inline T& at(const unsigned int index);
        inline T at(const unsigned int index) const;

        template <typename U>
        inline decltype(auto) apply_matrix_operator(const Matrix<U> &other, auto funct) const;
        template <typename U>
        inline decltype(auto) apply_scalar_operator(const U scalar, auto funct) const;

        inline void swap(Matrix<T> &other);
};

#endif
