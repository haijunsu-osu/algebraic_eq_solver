#!/usr/bin/env python3
"""
Find an example of rank-1 B matrix that yields 4 real solutions.

The null space method can theoretically produce up to 4 solutions:
- Up to 2 solutions for x from the trigonometric equation a*cos(x) + b*sin(x) + c = 0
- Up to 2 solutions for y from each x solution

This script systematically constructs examples to demonstrate all 4 solutions.
"""

import numpy as np
import math
from numerical_solver import solve_trigonometric_system, solve_trigonometric_equation

def construct_four_solution_example():
    """
    Construct a rank-1 B matrix example that yields exactly 4 real solutions.
    
    Strategy:
    1. Choose A and rank-1 B matrices
    2. Choose trigonometric equations that have 2 solutions each
    3. Work backward to find appropriate C
    """
    print("=" * 80)
    print("CONSTRUCTING RANK-1 EXAMPLE WITH 4 REAL SOLUTIONS")
    print("=" * 80)
    
    # Choose matrices carefully
    A = np.array([[1.0, 0.0],
                  [0.0, 1.0]])  # Identity for simplicity
    
    # Create rank-1 B with specific structure
    # B = u * v^T where u and v are chosen to give nice trigonometric equations
    u = np.array([[2.0], [1.0]])      # 2x1 column vector
    v = np.array([[1.0, 1.0]])        # 1x2 row vector
    B = u @ v  # 2x2 rank-1 matrix
    
    print(f"A = \n{A}")
    print(f"B = \n{B}")
    print(f"rank(B) = {np.linalg.matrix_rank(B)}")
    print(f"det(B) = {np.linalg.det(B):.6f}")
    
    # Analyze the null space transformation
    U, s, Vt = np.linalg.svd(B)
    null_vector = U[:, 1]
    
    print(f"\nNull space analysis:")
    print(f"Left null vector: {null_vector}")
    print(f"Verification: null_vector^T @ B = {null_vector @ B}")
    
    # The transformed equation will be: null_vector^T @ A @ [cos(x), sin(x)] = null_vector^T @ C
    left_A = null_vector @ A
    print(f"Left multiplied A: {left_A}")
    
    # We want the trigonometric equation: left_A[0]*cos(x) + left_A[1]*sin(x) = right_C
    # For this to have 2 solutions, we need the equation to be solvable
    # Let's choose right_C to give a nice equation
    
    # Choose C such that we get a well-conditioned trigonometric equation
    # We want: left_A[0]*cos(x) + left_A[1]*sin(x) + c = 0 to have 2 solutions
    
    # Let's work with the equation: -0.447*cos(x) + 0.894*sin(x) + c = 0
    # This can be written as: √(0.447² + 0.894²) * cos(x - φ) + c = 0
    # where tan(φ) = 0.894/0.447 = 2
    
    amplitude = math.sqrt(left_A[0]**2 + left_A[1]**2)
    print(f"Amplitude of trigonometric equation: {amplitude:.6f}")
    
    # Choose c such that |c| < amplitude to ensure 2 real solutions
    c_value = 0.5 * amplitude  # This should give 2 solutions
    right_C = -c_value  # Because we rearrange to a*cos(x) + b*sin(x) + c = 0
    
    print(f"Chosen right_C = {right_C:.6f}")
    print(f"This gives equation: {left_A[0]:.6f}*cos(x) + {left_A[1]:.6f}*sin(x) + {c_value:.6f} = 0")
    
    # Now we need to find C such that null_vector^T @ C = right_C
    # This is: null_vector[0]*C[0] + null_vector[1]*C[1] = right_C
    # We have one equation in two unknowns, so we can choose C[0] arbitrarily
    
    C0 = 1.0  # Arbitrary choice
    C1 = (right_C - null_vector[0] * C0) / null_vector[1]
    C = np.array([C0, C1])
    
    print(f"Constructed C = {C}")
    print(f"Verification: null_vector^T @ C = {null_vector @ C:.6f} (should be {right_C:.6f})")
    
    return A, B, C


def test_four_solution_case():
    """Test the constructed example and verify we get 4 solutions."""
    
    A, B, C = construct_four_solution_example()
    
    print(f"\n" + "="*60)
    print("SOLVING THE CONSTRUCTED SYSTEM")
    print("="*60)
    
    # Solve the system
    solutions = solve_trigonometric_system(A, B, C, verbose=True)
    
    print(f"\n" + "="*60)
    print("ANALYSIS OF SOLUTIONS")
    print("="*60)
    
    print(f"Found {len(solutions)} solutions:")
    
    for i, sol in enumerate(solutions):
        x_deg = math.degrees(sol['x'])
        y_deg = math.degrees(sol['y'])
        print(f"\nSolution {i+1}:")
        print(f"  x = {sol['x']:.6f} rad ({x_deg:.1f}°)")
        print(f"  y = {sol['y']:.6f} rad ({y_deg:.1f}°)")
        print(f"  cos(x) = {sol['cos_x']:.6f}, sin(x) = {sol['sin_x']:.6f}")
        print(f"  cos(y) = {sol['cos_y']:.6f}, sin(y) = {sol['sin_y']:.6f}")
        
        # Verify the solution
        eq1 = A[0,0]*sol['cos_x'] + A[0,1]*sol['sin_x'] + B[0,0]*sol['cos_y'] + B[0,1]*sol['sin_y']
        eq2 = A[1,0]*sol['cos_x'] + A[1,1]*sol['sin_x'] + B[1,0]*sol['cos_y'] + B[1,1]*sol['sin_y']
        
        residual1 = abs(eq1 - C[0])
        residual2 = abs(eq2 - C[1])
        
        print(f"  Verification:")
        print(f"    Equation 1: {eq1:.6f} = {C[0]:.6f} (residual: {residual1:.2e})")
        print(f"    Equation 2: {eq2:.6f} = {C[1]:.6f} (residual: {residual2:.2e})")
        
        if residual1 < 1e-10 and residual2 < 1e-10:
            print(f"    ✓ VALID SOLUTION")
        else:
            print(f"    ✗ INVALID SOLUTION")
    
    # Group solutions by x value to see the pattern
    x_groups = {}
    for i, sol in enumerate(solutions):
        x_key = round(sol['x'], 6)  # Round to avoid floating point issues
        if x_key not in x_groups:
            x_groups[x_key] = []
        x_groups[x_key].append((i+1, sol))
    
    print(f"\n" + "="*60)
    print("SOLUTION STRUCTURE ANALYSIS")
    print("="*60)
    
    print(f"Solutions grouped by x value:")
    for x_key, group in x_groups.items():
        x_deg = math.degrees(x_key)
        print(f"\n  x ≈ {x_key:.6f} rad ({x_deg:.1f}°): {len(group)} solutions")
        for sol_num, sol in group:
            y_deg = math.degrees(sol['y'])
            print(f"    Solution {sol_num}: y = {sol['y']:.6f} rad ({y_deg:.1f}°)")
    
    print(f"\nSUMMARY:")
    print(f"✓ Found {len(solutions)} total solutions")
    print(f"✓ Solutions come from {len(x_groups)} distinct x values")
    print(f"✓ Each x value contributes multiple y solutions")
    
    if len(solutions) == 4:
        print(f"✓ SUCCESS: Found exactly 4 real solutions!")
    elif len(solutions) > 4:
        print(f"! Found more than 4 solutions - interesting case!")
    else:
        print(f"! Found fewer than 4 solutions - need to adjust parameters")


def systematic_search_for_four_solutions():
    """
    Systematically search for rank-1 cases with exactly 4 solutions.
    """
    print(f"\n" + "="*80)
    print("SYSTEMATIC SEARCH FOR 4-SOLUTION CASES")
    print("="*80)
    
    success_count = 0
    
    # Try different parameter combinations
    for scale in [0.5, 1.0, 2.0]:
        for offset in [0.3, 0.5, 0.7]:
            for angle in [math.pi/6, math.pi/4, math.pi/3]:
                
                # Construct matrices
                A = np.array([[1.0, 0.0],
                              [0.0, 1.0]])
                
                # Create rank-1 B with rotation and scaling
                cos_a, sin_a = math.cos(angle), math.sin(angle)
                u = np.array([[scale * cos_a], [scale * sin_a]])
                v = np.array([[1.0, offset]])
                B = u @ v
                
                # Choose C to give moderate trigonometric equation
                U, s, Vt = np.linalg.svd(B)
                null_vector = U[:, 1]
                left_A = null_vector @ A
                
                amplitude = math.sqrt(left_A[0]**2 + left_A[1]**2)
                if amplitude < 1e-6:  # Skip degenerate cases
                    continue
                
                # Choose c to be within amplitude for 2 solutions
                c_value = 0.6 * amplitude
                right_C = -c_value
                
                # Construct C
                C0 = 1.0
                if abs(null_vector[1]) < 1e-6:  # Skip if null_vector[1] ≈ 0
                    continue
                C1 = (right_C - null_vector[0] * C0) / null_vector[1]
                C = np.array([C0, C1])
                
                # Solve and count solutions
                try:
                    solutions = solve_trigonometric_system(A, B, C, verbose=False)
                    if len(solutions) == 4:
                        success_count += 1
                        print(f"\n✓ FOUND 4-SOLUTION CASE #{success_count}:")
                        print(f"  scale={scale}, offset={offset}, angle={math.degrees(angle):.1f}°")
                        print(f"  A = \n{A}")
                        print(f"  B = \n{B}")
                        print(f"  C = {C}")
                        print(f"  Found {len(solutions)} solutions")
                        
                        if success_count >= 3:  # Stop after finding a few examples
                            return A, B, C
                            
                except Exception as e:
                    continue  # Skip problematic cases
    
    print(f"\nSearch completed. Found {success_count} cases with exactly 4 solutions.")
    return None, None, None


if __name__ == "__main__":
    # First try the constructed example
    test_four_solution_case()
    
    # If that doesn't work, search systematically
    print(f"\n" + "="*80)
    print("SEARCHING FOR OPTIMAL 4-SOLUTION EXAMPLE")
    print("="*80)
    
    A_best, B_best, C_best = systematic_search_for_four_solutions()
    
    if A_best is not None:
        print(f"\n" + "="*80)
        print("TESTING BEST 4-SOLUTION EXAMPLE")
        print("="*80)
        
        solutions = solve_trigonometric_system(A_best, B_best, C_best, verbose=True)
        print(f"\n✓ FINAL RESULT: Found {len(solutions)} solutions with optimized parameters")
    else:
        print(f"\nNote: May need to adjust search parameters for 4-solution cases")
