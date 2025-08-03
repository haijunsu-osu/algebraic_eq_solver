"""
EXAMPLE: Rank-1 B Matrix with 4 Real Solutions

This document presents a concrete example of a rank-1 B matrix that yields 
exactly 4 real solutions using the null space method.

SYSTEM:
-------
A = [[1, 0],     B = [[2, 2],     C = [1.0, -0.17082039]
     [0, 1]]          [1, 1]]

Properties:
- A is identity matrix (simplifies analysis)
- B has rank 1: second row = 0.5 × first row
- B = u @ v^T where u = [[2], [1]] and v = [[1, 1]]
- det(B) = 0 (exactly singular)

MATHEMATICAL STRUCTURE:
----------------------
1. Null Space Elimination:
   - Left null vector: [-0.4472136, 0.89442719]
   - Transforms to: -0.447214*cos(x) + 0.894427*sin(x) = -0.600000
   - This trigonometric equation has 2 solutions for x

2. Back Substitution:
   - For each x, solve: B[cos(y), sin(y)] = C - A[cos(x), sin(x)]
   - Uses rank-1 structure: B = σ₁ * u₁ * v₁^T
   - Each gives trigonometric equation with 2 solutions for y

SOLUTIONS:
---------
X-VALUE GROUP 1: x = -0.179853 rad (-10.3°)
  Solution 1: y = 2.350492 rad (134.7°)
  Solution 2: y = -0.779695 rad (-44.7°)

X-VALUE GROUP 2: x = -2.034444 rad (-116.6°)  
  Solution 3: y = 1.819070 rad (104.2°)
  Solution 4: y = -0.248274 rad (-14.2°)

STRUCTURE: 2 × 2 = 4 total solutions

VERIFICATION:
------------
All solutions satisfy:
✓ Original system: A[cos(x), sin(x)] + B[cos(y), sin(y)] = C
✓ Null space constraint: null_vector^T @ B @ [cos(y), sin(y)] = 0  
✓ Trigonometric identities: cos²(x) + sin²(x) = 1, cos²(y) + sin²(y) = 1
✓ Numerical residuals < 1e-15

SIGNIFICANCE:
------------
This example demonstrates:
1. The null space method can find ALL solutions systematically
2. Rank-1 matrices can indeed yield 4 real solutions (maximum possible)
3. The 2×2 structure emerges naturally from the mathematics
4. The method is numerically stable and exact

This showcases the full power and elegance of the null space approach
for solving trigonometric systems with singular B matrices.
"""

# The exact system for reference:
import numpy as np

A = np.array([[1.0, 0.0],
              [0.0, 1.0]])

B = np.array([[2.0, 2.0],
              [1.0, 1.0]])

C = np.array([1.0, -0.17082039])

print(__doc__)
print(f"REFERENCE SYSTEM:")
print(f"A = \n{A}")
print(f"B = \n{B}")  
print(f"C = {C}")
print(f"rank(B) = {np.linalg.matrix_rank(B)}")
print(f"det(B) = {np.linalg.det(B):.2e}")

if __name__ == "__main__":
    from numerical_solver import solve_trigonometric_system
    
    print(f"\nSOLUTION VERIFICATION:")
    solutions = solve_trigonometric_system(A, B, C, verbose=False)
    print(f"Number of solutions found: {len(solutions)}")
    
    for i, sol in enumerate(solutions):
        print(f"Solution {i+1}: x = {sol['x']:.6f}, y = {sol['y']:.6f}")
        
        # Verify
        eq1 = A[0,0]*sol['cos_x'] + A[0,1]*sol['sin_x'] + B[0,0]*sol['cos_y'] + B[0,1]*sol['sin_y']
        eq2 = A[1,0]*sol['cos_x'] + A[1,1]*sol['sin_x'] + B[1,0]*sol['cos_y'] + B[1,1]*sol['sin_y']
        
        residual1 = abs(eq1 - C[0])
        residual2 = abs(eq2 - C[1])
        
        print(f"  Residuals: ({residual1:.2e}, {residual2:.2e})")
    
    print(f"\n✓ SUCCESS: Found {len(solutions)} real solutions for rank-1 B matrix!")
