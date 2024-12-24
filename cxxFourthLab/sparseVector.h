#pragma once

#include <iostream>
#include <unordered_map>
#include <iterator>

template <typename T>
class SparseVector {
private:
    std::unordered_map<int, T> data;
    int size;

public:
    SparseVector(int size) : size(size) {}

    // Set value at index
    void set(int index, T value) {
        if (value != 0) {
            data[index] = value;
        }
        else {
            data.erase(index);
        }
    }

    // Get value at index
    T get(int index) const {
        auto it = data.find(index);
        return (it != data.end()) ? it->second : T(0);  // Return zero if not found
    }

    // Get size of the vector
    int getSize() const {
        return size;
    }
    T operator[](int index) const {
        auto it = data.find(index);
        return (it != data.end()) ? it->second : T{};
    }
    SparseVector<T> operator+(const SparseVector<T>& other) const {
        SparseVector<T> result(size);
        for (const auto& pair : data) {
            result.set(pair.first, pair.second+other[pair.first]);
        }
   /*     for (const auto& pair : other.data) {
            result.set(pair.first, result[pair.first] + pair.second);
        }*/
        return result;
    }
    SparseVector<T> operator*(const T scalar) const {
        SparseVector<T> result(size);
        for (const auto& pair : data) {
            result.set(pair.first, pair.second * scalar);
        }

        return result;
    }
    SparseVector<T> operator^(const T exponant) const {
        SparseVector<T> result(size);
        for (const auto& pair : data) {
            result.set(pair.first, std::pow(pair.second,exponant));
        }
       
        return result;
    }

    T dot(const SparseVector<T>& other) const {
        T result = 0;
        for (const auto& pair : data) {
            auto it = other.data.find(pair.first);
            if (it != other.data.end()) {
                result += pair.second * it->second;
            }
        }
        return result;
    }

    void print() const {
        for (const auto& pair : data) {
            std::cout << "Index: " << pair.first << ", Value: " << pair.second << std::endl;
        }
    }




    // Iterator class for SparseVector
    class Iterator {
    private:
        typename std::unordered_map<int, T>::iterator it;

    public:
        Iterator(typename std::unordered_map<int, T>::iterator iterator) : it(iterator) {}

        std::pair<int, T> operator*() {
            return *it;
        }

        Iterator& operator++() {
            ++it;
            return *this;
        }

        bool operator!=(const Iterator& other) const {
            return it != other.it;
        }
    };

    // Begin and end methods for iteration
    Iterator begin() {
        return Iterator(data.begin());
    }

    Iterator end() {
        return Iterator(data.end());
    }
};


