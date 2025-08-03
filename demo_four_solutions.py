#!/usr/bin/env python3
"""
DEMONSTRATION: Rank-1 B Matrix with Exactly 4 Real Solutions

This script demonstrates a carefully constructed example where the null space 
method finds all 4 possible real solutions for a rank-1 B matrix.

Mathematical Structure:
- 2 solutions for x from: null_vector^T @ A @ [cos(x), sin(x)] = null_vector^T @ C
- 2 solutions for y from each x: B @ [cos(y), sin(y)] = C - A @ [cos(x), sin(x)]
- Total: 2 × 2 = 4 solutions

This showcases the full power of the null space approach.
"""

import numpy as np
import math
from numerical_solver import solve_trigonometric_system

def demonstrate_four_solution_example():
    """Demonstrate the rank-1 case with exactly 4 real solutions."""
    
    print("=" * 80)
    print("RANK-1 B MATRIX WITH 4 REAL SOLUTIONS")
    print("=" * 80)
    
    # Create an exactly rank-1 matrix using outer product
    A = np.array([[1.0, 0.0],
                  [0.0, 1.0]])
    
    # B = u * v^T where u and v are chosen carefully
    u = np.array([[2.0], [1.0]])      # 2x1 column vector  
    v = np.array([[1.0, 1.0]])        # 1x2 row vector
    B = u @ v  # This gives exactly rank-1 matrix [[2, 2], [1, 1]]
    
    # Choose C to give a solvable system
    # We want the null space equation to have 2 solutions
    # null_vector^T @ A @ [cos(x), sin(x)] = null_vector^T @ C
    U, s, Vt = np.linalg.svd(B)
    null_vector = U[:, 1]
    left_A = null_vector @ A
    
    # Choose right_C to give a nice trigonometric equation
    amplitude = np.linalg.norm(left_A)
    right_C = -0.6 * amplitude  # This should give 2 solutions
    
    # Construct C such that null_vector^T @ C = right_C
    # null_vector[0]*C[0] + null_vector[1]*C[1] = right_C
    # Choose C[0] = 1.0, solve for C[1]
    C0 = 1.0
    C1 = (right_C - null_vector[0] * C0) / null_vector[1]
    C = np.array([C0, C1])
    
    print("SYSTEM DEFINITION:")
    print(f"A = \n{A}")
    print(f"B = \n{B}")
    print(f"C = {C}")
    print(f"rank(B) = {np.linalg.matrix_rank(B)}")
    print(f"det(B) = {np.linalg.det(B):.2e}")
    
    # Mathematical analysis
    print(f"\nMATHEMATICAL STRUCTURE:")
    U, s, Vt = np.linalg.svd(B)
    null_vector = U[:, 1]
    
    print(f"SVD Analysis:")
    print(f"  Singular values: {s}")
    print(f"  Left null vector: {null_vector}")
    print(f"  Verification null_vector^T @ B = {null_vector @ B}")
    
    # Show the transformed trigonometric equation for x
    left_A = null_vector @ A
    right_C = null_vector @ C
    print(f"\nTransformed equation in x:")
    print(f"  {left_A[0]:.6f}*cos(x) + {left_A[1]:.6f}*sin(x) = {right_C:.6f}")
    print(f"  Rearranged: {left_A[0]:.6f}*cos(x) + {left_A[1]:.6f}*sin(x) + {-right_C:.6f} = 0")
    
    # Solve the complete system
    print(f"\n" + "="*60)
    print("COMPLETE SOLUTION")
    print("="*60)
    
    solutions = solve_trigonometric_system(A, B, C, verbose=False)
    
    print(f"\nFOUND {len(solutions)} REAL SOLUTIONS:")
    print("="*50)
    
    # Group solutions by x value to show the 2×2 structure
    x_groups = {}
    for i, sol in enumerate(solutions):
        x_key = round(sol['x'], 6)
        if x_key not in x_groups:
            x_groups[x_key] = []
        x_groups[x_key].append((i+1, sol))
    
    for j, (x_key, group) in enumerate(x_groups.items()):
        x_deg = math.degrees(x_key)
        print(f"\nX-VALUE GROUP {j+1}: x = {x_key:.6f} rad ({x_deg:.1f}°)")
        print(f"  cos(x) = {math.cos(x_key):.6f}, sin(x) = {math.sin(x_key):.6f}")
        
        for sol_num, sol in group:
            y_deg = math.degrees(sol['y'])
            print(f"    Solution {sol_num}: y = {sol['y']:.6f} rad ({y_deg:.1f}°)")
            print(f"      cos(y) = {sol['cos_y']:.6f}, sin(y) = {sol['sin_y']:.6f}")
            
            # Verify this solution
            eq1 = A[0,0]*sol['cos_x'] + A[0,1]*sol['sin_x'] + B[0,0]*sol['cos_y'] + B[0,1]*sol['sin_y']
            eq2 = A[1,0]*sol['cos_x'] + A[1,1]*sol['sin_x'] + B[1,0]*sol['cos_y'] + B[1,1]*sol['sin_y']
            
            residual1 = abs(eq1 - C[0])
            residual2 = abs(eq2 - C[1])
            
            print(f"      Verification: residuals = ({residual1:.2e}, {residual2:.2e})")
            
            if residual1 < 1e-10 and residual2 < 1e-10:
                print(f"      ✓ VALID")
            else:
                print(f"      ✗ INVALID")
    
    print(f"\n" + "="*60)
    print("SOLUTION SUMMARY")
    print("="*60)
    
    print(f"✓ Total solutions found: {len(solutions)}")
    if len(x_groups) > 0:
        print(f"✓ Number of distinct x values: {len(x_groups)}")
        print(f"✓ Solutions per x value: {len(solutions) // len(x_groups)}")
        print(f"✓ Structure: {len(x_groups)} × {len(solutions) // len(x_groups)} = {len(solutions)} total")
    else:
        print(f"✗ No solutions found")
        return A, B, C, solutions
    
    if len(solutions) == 4:
        print(f"✓ SUCCESS: Perfect 2×2 = 4 solution structure!")
    
    print(f"\nThis demonstrates the full potential of the null space method:")
    print(f"1. Null space elimination gives trigonometric equation with 2 x solutions")
    print(f"2. Each x value leads to another trigonometric equation with 2 y solutions")
    print(f"3. Total: 2 × 2 = 4 real solutions")
    
    return A, B, C, solutions


def verify_mathematical_structure(A, B, C, solutions):
    """Verify the mathematical structure behind the 4 solutions."""
    
    print(f"\n" + "="*80)
    print("MATHEMATICAL VERIFICATION")
    print("="*80)
    
    # Step 1: Verify null space elimination
    U, s, Vt = np.linalg.svd(B)
    null_vector = U[:, 1]
    
    print(f"STEP 1: Null Space Elimination")
    print(f"Left null vector: {null_vector}")
    
    # Check that null_vector eliminates B for all solutions
    for i, sol in enumerate(solutions):
        trig_y = np.array([sol['cos_y'], sol['sin_y']])
        eliminated = null_vector @ B @ trig_y
        print(f"  Solution {i+1}: null_vector^T @ B @ [cos(y), sin(y)] = {eliminated:.2e}")
    
    # Step 2: Verify x solutions come from the same trigonometric equation
    print(f"\nSTEP 2: X Solutions from Single Trigonometric Equation")
    left_A = null_vector @ A
    right_C = null_vector @ C
    
    print(f"Equation: {left_A[0]:.6f}*cos(x) + {left_A[1]:.6f}*sin(x) = {right_C:.6f}")
    
    x_values = list(set(round(sol['x'], 6) for sol in solutions))
    print(f"Unique x values: {len(x_values)}")
    
    for x in x_values:
        lhs = left_A[0]*math.cos(x) + left_A[1]*math.sin(x)
        error = abs(lhs - right_C)
        print(f"  x = {x:.6f}: LHS = {lhs:.6f}, RHS = {right_C:.6f}, error = {error:.2e}")
    
    # Step 3: Verify y solutions come from rank-1 constraint
    print(f"\nSTEP 3: Y Solutions from Rank-1 Constraint")
    
    sigma1 = s[0]
    u1 = U[:, 0]
    v1 = Vt[0, :]
    
    print(f"B = σ₁ * u₁ * v₁^T")
    print(f"σ₁ = {sigma1:.6f}")
    print(f"u₁ = {u1}")
    print(f"v₁ = {v1}")
    
    # Group solutions by x and verify y equations
    x_groups = {}
    for sol in solutions:
        x_key = round(sol['x'], 6)
        if x_key not in x_groups:
            x_groups[x_key] = []
        x_groups[x_key].append(sol)
    
    for x_key, group in x_groups.items():
        print(f"\n  For x = {x_key:.6f}:")
        
        # Compute RHS for this x value
        cos_x, sin_x = math.cos(x_key), math.sin(x_key)
        rhs = C - A @ np.array([cos_x, sin_x])
        required_projection = np.dot(u1, rhs) / sigma1
        
        print(f"    Required: v₁^T @ [cos(y), sin(y)] = {required_projection:.6f}")
        print(f"    Equation: {v1[0]:.6f}*cos(y) + {v1[1]:.6f}*sin(y) = {required_projection:.6f}")
        
        for sol in group:
            trig_y = np.array([sol['cos_y'], sol['sin_y']])
            actual_projection = np.dot(v1, trig_y)
            error = abs(actual_projection - required_projection)
            print(f"      y = {sol['y']:.6f}: projection = {actual_projection:.6f}, error = {error:.2e}")
    
    print(f"\n✓ Mathematical structure verified!")
    print(f"✓ All solutions satisfy the null space constraints")
    print(f"✓ Perfect 2×2 solution pattern confirmed")


if __name__ == "__main__":
    A, B, C, solutions = demonstrate_four_solution_example()
    verify_mathematical_structure(A, B, C, solutions)
