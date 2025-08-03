#!/usr/bin/env python3
"""
Debug script for rank-1 B matrix case
"""

import numpy as np
import math
from numerical_solver import solve_trigonometric_system, analyze_singular_matrix

def test_rank1_debug():
    """Debug the rank-1 B case step by step."""
    print("=" * 60)
    print("DEBUGGING RANK-1 B MATRIX CASE")
    print("=" * 60)
    
    # Create a proper rank-1 matrix B
    B = np.array([[1.0, 2.0],
                  [0.5, 1.0]])  # rank 1: second row = 0.5 * first row
    
    print(f"B = \n{B}")
    print(f"det(B) = {np.linalg.det(B):.2e}")
    print(f"rank(B) = {np.linalg.matrix_rank(B)}")
    
    # Verify it's actually rank-1
    U, s, Vt = np.linalg.svd(B)
    print(f"Singular values: {s}")
    print(f"Rank-1 confirmed: {np.sum(s > 1e-12) == 1}")
    
    # Analyze the matrix
    analysis = analyze_singular_matrix(B)
    print(f"\nMatrix analysis:")
    print(f"  Type: {analysis['type']}")
    print(f"  Range direction: {analysis['U'][:, 0]}")
    print(f"  Null direction: {analysis['U'][:, 1]}")
    
    # Create a test case where we KNOW a solution should exist
    # Let's work backwards from a known solution
    A = np.array([[1.0, 0.0],
                  [0.0, 1.0]])  # Identity matrix
    
    # Choose specific angles
    x_known = math.pi / 6  # 30 degrees
    y_known = math.pi / 4  # 45 degrees
    
    cos_x = math.cos(x_known)
    sin_x = math.sin(x_known)
    cos_y = math.cos(y_known)
    sin_y = math.sin(y_known)
    
    print(f"\nKnown solution:")
    print(f"  x = {x_known:.6f} rad ({math.degrees(x_known):.1f}°)")
    print(f"  y = {y_known:.6f} rad ({math.degrees(y_known):.1f}°)")
    print(f"  cos(x) = {cos_x:.6f}, sin(x) = {sin_x:.6f}")
    print(f"  cos(y) = {cos_y:.6f}, sin(y) = {sin_y:.6f}")
    
    # Compute C from the known solution
    C = A @ np.array([cos_x, sin_x]) + B @ np.array([cos_y, sin_y])
    
    print(f"\nSystem: A[cos x, sin x] + B[cos y, sin y] = C")
    print(f"A = \n{A}")
    print(f"B = \n{B}")
    print(f"C = {C}")
    
    # Check if the known solution actually satisfies the equations
    eq1 = A[0,0]*cos_x + A[0,1]*sin_x + B[0,0]*cos_y + B[0,1]*sin_y
    eq2 = A[1,0]*cos_x + A[1,1]*sin_x + B[1,0]*cos_y + B[1,1]*sin_y
    
    print(f"\nVerification of known solution:")
    print(f"  Equation 1: {eq1:.6f} = {C[0]:.6f} (diff: {abs(eq1-C[0]):.2e})")
    print(f"  Equation 2: {eq2:.6f} = {C[1]:.6f} (diff: {abs(eq2-C[1]):.2e})")
    
    # Now check the geometric constraint
    print(f"\n" + "="*50)
    print("CHECKING GEOMETRIC CONSTRAINTS")
    print("="*50)
    
    # For the known solution, compute required_vector = C - A[cos(x), sin(x)]
    required_vector = C - A @ np.array([cos_x, sin_x])
    print(f"Required vector for known solution: {required_vector}")
    
    # Check if this vector is in the range of B
    range_direction = analysis['U'][:, 0]
    print(f"Range direction of B: {range_direction}")
    
    # Check if required_vector is parallel to range_direction
    if np.linalg.norm(required_vector) > 1e-12:
        required_unit = required_vector / np.linalg.norm(required_vector)
        dot_product = abs(np.dot(required_unit, range_direction))
        print(f"Dot product (should be ~1 for parallel): {dot_product:.6f}")
        print(f"Vectors are parallel: {dot_product > 1 - 1e-10}")
        
        # Check if B can produce the required vector
        B_pinv = np.linalg.pinv(B)
        trig_y_from_B = B_pinv @ required_vector
        print(f"Pseudo-inverse result: cos(y)={trig_y_from_B[0]:.6f}, sin(y)={trig_y_from_B[1]:.6f}")
        print(f"Identity check: {trig_y_from_B[0]**2 + trig_y_from_B[1]**2:.6f}")
        
        # Verify: B * trig_y_from_B should equal required_vector
        B_times_trig = B @ trig_y_from_B
        print(f"B * pseudo_inverse_result = {B_times_trig}")
        print(f"Should equal required_vector = {required_vector}")
        print(f"Difference: {np.linalg.norm(B_times_trig - required_vector):.2e}")
    
    # Now run the actual solver
    print(f"\n" + "="*50)
    print("RUNNING SOLVER")
    print("="*50)
    
    solutions = solve_trigonometric_system(A, B, C, verbose=True)
    print(f"\nFound {len(solutions)} solutions")
    
    if len(solutions) == 0:
        print("❌ No solutions found - debugging needed!")
    else:
        for i, sol in enumerate(solutions):
            print(f"\nSolution {i+1}:")
            print(f"  x = {sol['x']:.6f} ({math.degrees(sol['x']):.1f}°)")
            print(f"  y = {sol['y']:.6f} ({math.degrees(sol['y']):.1f}°)")


if __name__ == "__main__":
    test_rank1_debug()
