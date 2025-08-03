#!/usr/bin/env python3
"""
Test the new null space approach for rank-1 B matrices.

This script demonstrates the mathematical elegance of the null space method:
1. Left multiply by null vector to eliminate y terms
2. Solve trigonometric equation for x analytically 
3. Substitute back to solve for y
4. Get up to 4 solutions total (2 for x × 2 for y each)
"""

import numpy as np
import math
from numerical_solver import solve_trigonometric_system, solve_trigonometric_equation

def test_null_space_method():
    print("=" * 80)
    print("TESTING NULL SPACE METHOD FOR RANK-1 B MATRICES")
    print("=" * 80)
    
    # Create a rank-1 matrix B where second row = k * first row
    print("\nTest Case 1: B with second row = 0.5 * first row")
    print("-" * 50)
    
    A = np.array([[1.0, 0.0],
                  [0.0, 1.0]])
    
    B = np.array([[2.0, 1.0],
                  [1.0, 0.5]])  # rank 1: second row = 0.5 * first row
    
    print(f"A = \n{A}")
    print(f"B = \n{B}")
    print(f"rank(B) = {np.linalg.matrix_rank(B)}")
    print(f"det(B) = {np.linalg.det(B):.6f}")
    
    # Use a known solution to generate C
    x_ref = math.pi/6  # 30 degrees
    y_ref = math.pi/4  # 45 degrees
    cos_x_ref = math.cos(x_ref)
    sin_x_ref = math.sin(x_ref)
    cos_y_ref = math.cos(y_ref)
    sin_y_ref = math.sin(y_ref)
    
    C = A @ np.array([cos_x_ref, sin_x_ref]) + B @ np.array([cos_y_ref, sin_y_ref])
    
    print(f"Reference solution: x = {x_ref:.6f}, y = {y_ref:.6f}")
    print(f"C = {C}")
    
    # Solve using our method
    solutions = solve_trigonometric_system(A, B, C, verbose=True)
    
    print(f"\nFound {len(solutions)} solutions:")
    for i, sol in enumerate(solutions):
        x_deg = math.degrees(sol['x'])
        y_deg = math.degrees(sol['y'])
        print(f"  Solution {i+1}: x = {sol['x']:.6f} ({x_deg:.1f}°), y = {sol['y']:.6f} ({y_deg:.1f}°)")
        
        # Check if this matches reference
        x_err = min(abs(sol['x'] - x_ref), abs(sol['x'] - x_ref - 2*math.pi), abs(sol['x'] - x_ref + 2*math.pi))
        y_err = min(abs(sol['y'] - y_ref), abs(sol['y'] - y_ref - 2*math.pi), abs(sol['y'] - y_ref + 2*math.pi))
        
        if x_err < 1e-6 and y_err < 1e-6:
            print(f"    ✓ MATCHES REFERENCE SOLUTION")
    
    # Test Case 2: Different rank-1 matrix
    print("\n" + "=" * 80)
    print("Test Case 2: Different rank-1 B matrix")
    print("-" * 50)
    
    A2 = np.array([[1.5, -0.5],
                   [0.3, 2.0]])
    
    # Create B with rank 1: B = u * v^T
    u = np.array([[1.0], [2.0]])  # 2x1 column vector
    v = np.array([[3.0, -1.0]])   # 1x2 row vector
    B2 = u @ v  # 2x2 rank-1 matrix
    
    print(f"A = \n{A2}")
    print(f"B = \n{B2}")
    print(f"rank(B) = {np.linalg.matrix_rank(B2)}")
    print(f"det(B) = {np.linalg.det(B2):.6f}")
    
    # Use different reference solution
    x_ref2 = math.pi/3  # 60 degrees
    y_ref2 = math.pi/6  # 30 degrees
    cos_x_ref2 = math.cos(x_ref2)
    sin_x_ref2 = math.sin(x_ref2)
    cos_y_ref2 = math.cos(y_ref2)
    sin_y_ref2 = math.sin(y_ref2)
    
    C2 = A2 @ np.array([cos_x_ref2, sin_x_ref2]) + B2 @ np.array([cos_y_ref2, sin_y_ref2])
    
    print(f"Reference solution: x = {x_ref2:.6f}, y = {y_ref2:.6f}")
    print(f"C = {C2}")
    
    # Solve using our method
    solutions2 = solve_trigonometric_system(A2, B2, C2, verbose=True)
    
    print(f"\nFound {len(solutions2)} solutions:")
    for i, sol in enumerate(solutions2):
        x_deg = math.degrees(sol['x'])
        y_deg = math.degrees(sol['y'])
        print(f"  Solution {i+1}: x = {sol['x']:.6f} ({x_deg:.1f}°), y = {sol['y']:.6f} ({y_deg:.1f}°)")
        
        # Check if this matches reference
        x_err = min(abs(sol['x'] - x_ref2), abs(sol['x'] - x_ref2 - 2*math.pi), abs(sol['x'] - x_ref2 + 2*math.pi))
        y_err = min(abs(sol['y'] - y_ref2), abs(sol['y'] - y_ref2 - 2*math.pi), abs(sol['y'] - y_ref2 + 2*math.pi))
        
        if x_err < 1e-6 and y_err < 1e-6:
            print(f"    ✓ MATCHES REFERENCE SOLUTION")


def test_trigonometric_equation_solver():
    """Test the helper function for solving a*cos(θ) + b*sin(θ) + c = 0"""
    print("\n" + "=" * 80)
    print("TESTING TRIGONOMETRIC EQUATION SOLVER")
    print("=" * 80)
    
    test_cases = [
        (1, 1, -1, "cos(θ) + sin(θ) = 1"),
        (3, 4, -5, "3*cos(θ) + 4*sin(θ) = 5"),
        (1, 0, -0.5, "cos(θ) = 0.5"),
        (0, 1, -math.sqrt(3)/2, f"sin(θ) = {math.sqrt(3)/2:.6f}"),
        (1, math.sqrt(3), -2, "cos(θ) + √3*sin(θ) = 2"),
    ]
    
    for i, (a, b, c, description) in enumerate(test_cases):
        print(f"\nTest Case {i+1}: {description}")
        print(f"Equation: {a:.6f}*cos(θ) + {b:.6f}*sin(θ) + {c:.6f} = 0")
        
        solutions = solve_trigonometric_equation(a, b, c, verbose=True)
        
        print(f"Solutions: {len(solutions)} found")
        for j, theta in enumerate(solutions):
            theta_deg = math.degrees(theta)
            print(f"  θ_{j+1} = {theta:.6f} rad ({theta_deg:.1f}°)")
            
            # Verify the solution
            residual = a * math.cos(theta) + b * math.sin(theta) + c
            print(f"    Verification: {a:.3f}*cos({theta:.3f}) + {b:.3f}*sin({theta:.3f}) + {c:.3f} = {residual:.2e}")


if __name__ == "__main__":
    test_trigonometric_equation_solver()
    test_null_space_method()
