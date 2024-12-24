# Реализация разреженных векторов и матриц на C++

## Обзор
Этот репозиторий содержит реализацию разреженных векторов и матриц на C++. Цель заключалась в оптимизации использования памяти для структур данных, которые в основном содержат нулевые значения. Реализация включает в себя различные операции и сравнения производительности с плотными матрицами и векторами.

### Разделение задач

1. **Структуры данных для разреженных векторов и матриц**
   - Реализованы классы `SparseVector` и `SparseMatrix`, использующие `std::unordered_map` для хранения ненулевых элементов, что эффективно экономит память, не храня нули.
   - Фрагмент кода:
     ```cpp
     template <typename T>
     class SparseVector {
     private:
         std::unordered_map<int, T> data;
         int size;
     public:
         SparseVector(int size) : size(size) {}
         void set(int index, T value) {
             if (value != 0) {
                 data[index] = value;
             } else {
                 data.erase(index);
             }
         }
         T get(int index) const {
             auto it = data.find(index);
             return (it != data.end()) ? it->second : T(0);} };
     ```
   -  и для матриц
     ```cpp
    
    struct pair_hash
    {
        int operator()(const std::pair<int, int>& values) const noexcept
        {
            return (std::hash<int>()(values.first)
                ^ std::hash<int>()(values.second));
        }
    };
    template <typename T>
    class SparseMatrix {
    private:
        std::unordered_map<std::pair<int, int>, T, pair_hash> data;
        int rows;
        int cols;
    
    public:
        SparseMatrix(int rows, int cols) : rows(rows), cols(cols) {}
        SparseMatrix() = default;
    
        // Установить значение в (строка, столбец)
        void set(int row, int col, T value) {
            if (value != 0) {
                data[{row, col}] = value;  // Хранить ненулевые значения
            }
            else {
                data.erase({ row, col });  } }
     ```
   -  и я поддержал свои собственные итераторы, используя класс итератора, подобный этому
      ```cpp
        class Iterator {
          private:
              typename std::unordered_map<std::pair<int, int>, T>::iterator it;
      
          public:
              Iterator(typename std::unordered_map<std::pair<int, int>, T>::iterator iterator) : it(iterator) {}
      
              std::pair<std::pair<int, int>, T> operator*() {
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
      
          // Методы begin и end для итерации
          Iterator begin() {
              return Iterator(data.begin());
          }
      
          Iterator end() {
              return Iterator(data.end());
          }
      };
        
      ```

2. **Унарные и бинарные операции для векторов и матриц**
   - Реализованы операции, такие как сложение, умножение на скаляр, скалярное произведение для векторов и сложение, умножение и транспонирование для матриц. 
   - Включена функциональность для инверсии матриц и возведения в степень.
   - Фрагмент кода для сложения:
   - пример операторов векторов
     ```cpp
     SparseVector<T> operator+(const SparseVector<T>& other) const {
         SparseVector<T> result(size);
         for (const auto& pair : data) {
             result.set(pair.first, pair.second + other[pair.first]);
         }
         return result; }
     ```
   - инверсия матрицы с использованием метода Гаусса-Жордана
    ```cpp
    SparseMatrix<T> inverse() const {

    int aug_rows = rows;
    int aug_cols = cols * 2;
    SparseMatrix<T> aug_matrix(aug_rows, aug_cols);
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is singular!!");
    }
    for (auto element : data) {
        aug_matrix.set(element.first.first, element.first.second, element.second) ;
    }
    for (int i = 0; i < aug_rows; i++) {
      

        aug_matrix.set(i, rows + i, 1);
    }
    // Perform Gauss-Jordan elimination
    for (int i = 0; i < aug_rows; i++) {
        // Find pivot
        T pivot = aug_matrix.get(i, i);
        if (pivot == 0) {
            throw std::invalid_argument("Matrix is singular and cannot be inverted!");
        }

        // Normalize the pivot row
        for (int j = 0; j < aug_cols; j++) {
            aug_matrix.set(i, j, aug_matrix.get(i, j) / pivot);
        }

        // Eliminate the current column in all other rows
        for (int k = 0; k < aug_rows; k++) {
            if (k != i) {
                T factor = aug_matrix.get(k, i);
                for (int j = 0; j < aug_cols; j++) {
                    aug_matrix.set(k, j, aug_matrix.get(k, j) - factor * aug_matrix.get(i, j));
                }
            }
        }
    }

    // Extract the right half (the inverted matrix)
    SparseMatrix<T> inverse_matrix(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = cols; j < aug_cols; j++) {
            T value = aug_matrix.get(i, j);
            if (value != 0) { // Only store non-zero values
                inverse_matrix.set(i, j - cols, value);
            }
        }
    }

    return inverse_matrix;}
    ```
   

   - Я использовал ряд для вычисления экспоненты матрицы 
  
    ```cpp
    SparseMatrix<T> exp() const {
     if (rows != cols) {
         throw std::invalid_argument("Матрица должна быть квадратной для возведения в степень");
     }

     SparseMatrix<T> result(rows, cols);
     SparseMatrix<T> term=(*this);
     for (int i = 0; i < rows; ++i) {
         result.set(i, i, T(1));
     }
     result = result + term;
     int n = 2;

     // Вычисление e^(A) = I + A + (A^2)/2! + (A^3)/3! + ...
     while (n<100){ // Ограничить ряд до разумного числа членов
       
             long long int f = factorial(n);
             term = term * (*this);
             term = term / static_cast<T>(n); 
             result = result + term;
         ++n;
     }
     return result;}
    ```

 - Я также предоставил функцию возведения в степень с двумя вариантами для вещественных и целых экспонентов для матриц
     

3. **Элементные операции**
   - Реализованы операции для умножения на скаляр и поэлементного возведения в степень.
   - Фрагмент кода:
   - для векторов
     ```cpp
     SparseVector<T> operator*(const T scalar) const {
         SparseVector<T> result(size);
         for (const auto& pair : data) {
             result.set(pair.first, pair.second * scalar);
         }
         return result;
     }
     ```
   - и для матриц
     ```cpp
      SparseMatrix<T> operator*(const T& scalar) const {
     SparseMatrix<T> result(rows, cols);
     for (const auto& pair : data) {
         result.set(pair.first.first, pair.first.second, pair.second * scalar);
     }
     return result;}
     ```

4. **Сравнение производительности**
   - Созданы функции для измерения времени выполнения операций как для разреженных, так и для плотных матриц и векторов.
   - Проведены обширные тесты для сравнения их производительности при различных операциях.
   - созданы стандартные плотные векторы и матрицы для сравнения
   - измерено время выполнения для различных операций
   - Фрагмент кода для измерения времени для разреженной матрицы:
     ```cpp
     template<typename Func>
     double measureSparseMatrixOperation(const SparseMatrix<int>& sparseMatrixA,
         const SparseMatrix<int>& sparseMatrixB,
         Func operation) {
         auto start = std::chrono::high_resolution_clock::now();
         operation(sparseMatrixA, sparseMatrixB);
         auto end = std::chrono::high_resolution_clock::now();
         return std::chrono::duration<double, std::milli>(end - start).count();
     }
     ```
   - Фрагмент кода для измерения времени для плотной матрицы:
     ```cpp
     template<typename Func>
double measureDenseMatrixOperation(const std::vector<std::vector<int>>& matrixA,
    const std::vector<std::vector<int>>& matrixB,
    Func operation) {
    auto start = std::chrono::high_resolution_clock::now();
    operation(matrixA, matrixB);
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count(); // Возвращает время в миллисекундах
}
     ```
   -  и я сделал то же самое для векторов.
    

5. **Тестирование и результаты**
   - Созданы тестовые случаи для проверки функциональности как разреженных, так и плотных структур данных.
   - Результаты тестов приведены ниже, демонстрируя разницу в производительности между двумя реализациями.

### Выводы тестирования
```plaintext
====== Тестирование матриц ========
Время сложения плотной матрицы: 30.88 мс Время сложения разреженной матрицы: 15.46 мс
Время умножения плотной матрицы: 34272.80 мс Время умножения разреженной матрицы: 10.73 мс
Время умножения плотной матрицы на скаляр: 34272.80 мс Время умножения разреженной матрицы на скаляр: 10.73 мс
====== Тестирование векторов ========
Время сложения плотного вектора: 0.19 мс  Время сложения разреженного вектора: 0.14 мс
Время скалярного произведения плотного вектора: 0.14 мс  Время скалярного произведения разреженного вектора: 0.03 мс
Время умножения плотного вектора на скаляр: 0.13 мс  Время умножения разреженного вектора на скаляр: 0.09 мс
```

### Комментарии к результатам
- Результаты указывают на значительное преимущество производительности разреженных структур данных по сравнению с их плотными аналогами, особенно в операциях с большими матрицами (1000x1000) или векторами (10000 элементов).
- Сложение разреженной матрицы было примерно в два раза быстрее, чем плотное сложение, тогда как операции умножения показали еще большее расхождение (умножение разреженной матрицы заняло всего 10.73 мс по сравнению с более чем 34 секундами для плотного умножения).
- Результаты подтверждают, что использование разреженных представлений может привести к значительным экономиям памяти и времени при работе с большими наборами данных, преимущественно заполненными нулями.

### Заключение
Эта реализация эффективно решает проблему эффективного хранения и манипуляции разреженными векторами и матрицами на C++. Использование хеш-таблиц для хранения только ненулевых элементов позволяет значительно сэкономить память, в то время как общая производительность демонстрируется через различные операции и тесты.
