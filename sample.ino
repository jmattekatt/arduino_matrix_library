#include <matrix.hpp>

void setup() {
  Serial.begin(9600);
}

void loop() {
  float a_init[] = {1, 1, 1, 4, 3, -1, 3, 5, 3};
  Matrix<float> A(3, 3, a_init); // 3x3 Matrix A

  float b_init[] = {1, 6, 4};
  Matrix<float> B(3, 1, b_init); // 3x1 Matrix B

  Matrix<float> L, U;
  A.LU_decompose(L, U);

  Serial.println("A");
  A.print();
  Serial.println("B");
  B.print();

  Serial.println("L: Lower-triangular Matrix");
  L.print();
  Serial.println("U: Upper-triangular Matrix");
  U.print();
  
  Serial.println("Check L*U = A");
  auto check_A = L*U;
  check_A.print();
  ASSERT((A == check_A), "ERROR: A != check_A");

  auto x = A.solve_for(B);
  Serial.println("x: Answer to A*x = B");
  x.print();
  
  Serial.println("Check A*x = B");
  auto check_B = A*x;
  check_B.print();
  ASSERT((B == check_B), "ERROR: B != check_B");
  
  Serial.println();
  Serial.println();
}
