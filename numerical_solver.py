#!/usr/bin/env python3
"""
Numerical Solver for Trigonometric System

Solves the system: A[cos x, sin x] + B[cos y, sin y] = C

This script implements the complete algorithm:
1. Express cos(y) and sin(y) in terms of cos(x) and sin(x) 
2. Apply the identity cos²(y) + sin²(y) = 1
3. Use Weierstrass substitution to get a quartic polynomial in t
4. Solve for t using numpy.roots
5. Convert back to x and y

Dependencies: numpy, math only

Author: AI Assistant  
Date: August 3, 2025
"""

import numpy as np
import math


def solve_trigonometric_system(A, B, C, verbose=False):
    """
    Solve the trigonometric system A[cos x, sin x] + B[cos y, sin y] = C.
    
    Args:
        A: 2x2 numpy array
        B: 2x2 numpy array  
        C: 2x1 numpy array
        verbose: Print detailed steps
        
    Returns:
        List of solution dictionaries, each containing 'x', 'y', 'cos_x', 'sin_x', 'cos_y', 'sin_y'
    """
    solutions = []
    
    if verbose:
        print("Solving system: A[cos x, sin x] + B[cos y, sin y] = C")
        print(f"A = \n{A}")
        print(f"B = \n{B}") 
        print(f"C = {C}")
    
    # Extract matrix elements for clarity
    A00, A01 = A[0, 0], A[0, 1]
    A10, A11 = A[1, 0], A[1, 1]
    B00, B01 = B[0, 0], B[0, 1]
    B10, B11 = B[1, 0], B[1, 1]
    C0, C1 = C[0], C[1]
    
    # Check if B is invertible
    det_B = B00*B11 - B01*B10
    if abs(det_B) < 1e-12:
        print("Warning: Matrix B is singular or nearly singular")
        return solutions
    
    if verbose:
        print(f"\nStep 1: Express cos(y) and sin(y) in terms of cos(x) and sin(x)")
        print(f"From B[cos(y), sin(y)]^T = C - A[cos(x), sin(x)]^T")
    
    # Step 2: Build the quartic polynomial coefficients directly
    # This implements the full symbolic derivation numerically
    
    if verbose:
        print(f"\nStep 2: Building quartic polynomial coefficients")
    
    # These coefficients come from the symbolic derivation:
    # After applying Weierstrass substitution and trigonometric identity
    
    # Denominator for all coefficients  
    denom = det_B**2
    
    # Common terms that appear in the coefficients
    A00_sq = A00**2
    A01_sq = A01**2  
    A10_sq = A10**2
    A11_sq = A11**2
    B00_sq = B00**2
    B01_sq = B01**2
    B10_sq = B10**2
    B11_sq = B11**2
    C0_sq = C0**2
    C1_sq = C1**2
    
    # Build the numerators for each coefficient
    # a4 coefficient (t^4 term)
    a4_num = (A00_sq*B10_sq + A00_sq*B11_sq - 2*A00*A10*B00*B10 - 2*A00*A10*B01*B11 
              - 2*A00*B00*B10*C1 - 2*A00*B01*B11*C1 + 2*A00*B10_sq*C0 + 2*A00*B11_sq*C0 
              + A10_sq*B00_sq + A10_sq*B01_sq + 2*A10*B00_sq*C1 - 2*A10*B00*B10*C0 
              + 2*A10*B01_sq*C1 - 2*A10*B01*B11*C0 - B00_sq*B11_sq + B00_sq*C1_sq 
              + 2*B00*B01*B10*B11 - 2*B00*B10*C0*C1 - B01_sq*B10_sq + B01_sq*C1_sq 
              - 2*B01*B11*C0*C1 + B10_sq*C0_sq + B11_sq*C0_sq)
    
    # a3 coefficient (t^3 term)
    a3_num = 4*(-A00*A01*B10_sq - A00*A01*B11_sq + A00*A11*B00*B10 + A00*A11*B01*B11 
                + A01*A10*B00*B10 + A01*A10*B01*B11 + A01*B00*B10*C1 + A01*B01*B11*C1 
                - A01*B10_sq*C0 - A01*B11_sq*C0 - A10*A11*B00_sq - A10*A11*B01_sq 
                - A11*B00_sq*C1 + A11*B00*B10*C0 - A11*B01_sq*C1 + A11*B01*B11*C0)
    
    # a2 coefficient (t^2 term)  
    a2_num = 2*(-A00_sq*B10_sq - A00_sq*B11_sq + 2*A00*A10*B00*B10 + 2*A00*A10*B01*B11 
                + 2*A01_sq*B10_sq + 2*A01_sq*B11_sq - 4*A01*A11*B00*B10 - 4*A01*A11*B01*B11 
                - A10_sq*B00_sq - A10_sq*B01_sq + 2*A11_sq*B00_sq + 2*A11_sq*B01_sq 
                - B00_sq*B11_sq + B00_sq*C1_sq + 2*B00*B01*B10*B11 - 2*B00*B10*C0*C1 
                - B01_sq*B10_sq + B01_sq*C1_sq - 2*B01*B11*C0*C1 + B10_sq*C0_sq + B11_sq*C0_sq)
    
    # a1 coefficient (t^1 term)
    a1_num = 4*(A00*A01*B10_sq + A00*A01*B11_sq - A00*A11*B00*B10 - A00*A11*B01*B11 
                - A01*A10*B00*B10 - A01*A10*B01*B11 + A01*B00*B10*C1 + A01*B01*B11*C1 
                - A01*B10_sq*C0 - A01*B11_sq*C0 + A10*A11*B00_sq + A10*A11*B01_sq 
                - A11*B00_sq*C1 + A11*B00*B10*C0 - A11*B01_sq*C1 + A11*B01*B11*C0)
    
    # a0 coefficient (constant term)
    a0_num = (A00_sq*B10_sq + A00_sq*B11_sq - 2*A00*A10*B00*B10 - 2*A00*A10*B01*B11 
              + 2*A00*B00*B10*C1 + 2*A00*B01*B11*C1 - 2*A00*B10_sq*C0 - 2*A00*B11_sq*C0 
              + A10_sq*B00_sq + A10_sq*B01_sq - 2*A10*B00_sq*C1 + 2*A10*B00*B10*C0 
              - 2*A10*B01_sq*C1 + 2*A10*B01*B11*C0 - B00_sq*B11_sq + B00_sq*C1_sq 
              + 2*B00*B01*B10*B11 - 2*B00*B10*C0*C1 - B01_sq*B10_sq + B01_sq*C1_sq 
              - 2*B01*B11*C0*C1 + B10_sq*C0_sq + B11_sq*C0_sq)
    
    # Final coefficients
    a4 = a4_num / denom
    a3 = a3_num / denom  
    a2 = a2_num / denom
    a1 = a1_num / denom
    a0 = a0_num / denom
    
    if verbose:
        print(f"Quartic polynomial: {a4:.6f}*t^4 + {a3:.6f}*t^3 + {a2:.6f}*t^2 + {a1:.6f}*t + {a0:.6f} = 0")
    
    # Form coefficient array for numpy.roots (highest degree first)
    coeffs = [a4, a3, a2, a1, a0]
    
    # Remove leading zeros for numerical stability
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-12:
        coeffs = coeffs[1:]
    
    if len(coeffs) == 1:
        if abs(coeffs[0]) < 1e-12:
            if verbose:
                print("Warning: All coefficients are zero - infinite solutions")
            return solutions
        else:
            if verbose:
                print("Warning: No solutions - constant polynomial")
            return solutions
    
    # Step 3: Solve quartic polynomial for t
    if verbose:
        print(f"\nStep 3: Solving polynomial of degree {len(coeffs)-1}")
    
    try:
        t_roots = np.roots(coeffs)
    except np.linalg.LinAlgError:
        if verbose:
            print("Warning: Failed to find polynomial roots")
        return solutions
    
    if verbose:
        print(f"Found {len(t_roots)} roots: {t_roots}")
    
    # Step 4: Process each root
    for i, t in enumerate(t_roots):
        if verbose:
            print(f"\nProcessing root {i+1}: t = {t}")
            
        # Skip complex roots with significant imaginary part
        if np.iscomplex(t) and abs(t.imag) > 1e-10:
            if verbose:
                print("  Skipping complex root")
            continue
            
        t = t.real if np.iscomplex(t) else t
        
        # Step 5: Convert t back to cos(x) and sin(x) using inverse Weierstrass substitution
        denominator = 1 + t**2
        cos_x = (1 - t**2) / denominator
        sin_x = 2*t / denominator
        
        if verbose:
            print(f"  cos(x) = {cos_x:.6f}, sin(x) = {sin_x:.6f}")
        
        # Verify trigonometric identity
        if abs(cos_x**2 + sin_x**2 - 1) > 1e-10:
            if verbose:
                print("  Failed trigonometric identity check for x")
            continue
            
        # Calculate x from cos_x and sin_x
        x = math.atan2(sin_x, cos_x)
        
        # Step 6: Solve for cos(y) and sin(y) using the original equations
        # B * [cos_y, sin_y]^T = C - A * [cos_x, sin_x]^T
        rhs = C - A @ np.array([cos_x, sin_x])
        
        try:
            trig_y = np.linalg.solve(B, rhs)
            cos_y, sin_y = trig_y[0], trig_y[1]
        except np.linalg.LinAlgError:
            if verbose:
                print("  Failed to solve for cos(y), sin(y)")
            continue
            
        if verbose:
            print(f"  cos(y) = {cos_y:.6f}, sin(y) = {sin_y:.6f}")
            
        # Verify trigonometric identity for y
        if abs(cos_y**2 + sin_y**2 - 1) > 1e-10:
            if verbose:
                print("  Failed trigonometric identity check for y")
            continue
            
        # Calculate y from cos_y and sin_y
        y = math.atan2(sin_y, cos_y)
        
        # Step 7: Verify the original equations
        eq1_residual = A[0,0]*cos_x + A[0,1]*sin_x + B[0,0]*cos_y + B[0,1]*sin_y - C[0]
        eq2_residual = A[1,0]*cos_x + A[1,1]*sin_x + B[1,0]*cos_y + B[1,1]*sin_y - C[1]
        
        if verbose:
            print(f"  Equation residuals: {eq1_residual:.2e}, {eq2_residual:.2e}")
        
        if abs(eq1_residual) < 1e-10 and abs(eq2_residual) < 1e-10:
            solutions.append({
                'x': x,
                'y': y,
                'cos_x': cos_x,
                'sin_x': sin_x,
                'cos_y': cos_y,
                'sin_y': sin_y,
                't': t
            })
            if verbose:
                print(f"  ✓ Valid solution: x = {x:.6f}, y = {y:.6f}")
        elif verbose:
            print("  ✗ Failed equation verification")
    
    return solutions


def test_solver():
    """Test the solver with a known solution."""
    print("=" * 60)
    print("TESTING THE NUMERICAL SOLVER")
    print("=" * 60)
    
    # Create a test case with known solution
    x_true = math.pi / 6  # 30 degrees
    y_true = math.pi / 4  # 45 degrees
    
    cos_x_true = math.cos(x_true)
    sin_x_true = math.sin(x_true)
    cos_y_true = math.cos(y_true)
    sin_y_true = math.sin(y_true)
    
    print(f"True solution:")
    print(f"  x = {x_true:.6f} rad ({math.degrees(x_true):.1f}°)")
    print(f"  y = {y_true:.6f} rad ({math.degrees(y_true):.1f}°)")
    print(f"  cos(x) = {cos_x_true:.6f}, sin(x) = {sin_x_true:.6f}")
    print(f"  cos(y) = {cos_y_true:.6f}, sin(y) = {sin_y_true:.6f}")
    
    # Define test matrices
    A = np.array([[2.0, 1.0],
                  [1.0, 3.0]])
    
    B = np.array([[1.5, -0.5],
                  [0.5, 2.0]])
    
    # Compute C from the known solution
    C = A @ np.array([cos_x_true, sin_x_true]) + B @ np.array([cos_y_true, sin_y_true])
    
    print(f"\nTest system:")
    print(f"A = \n{A}")
    print(f"B = \n{B}")
    print(f"C = {C}")
    
    # Solve the system
    print(f"\n{'='*40}")
    print("SOLVING...")
    print(f"{'='*40}")
    
    solutions = solve_trigonometric_system(A, B, C, verbose=True)
    
    print(f"\n{'='*40}")
    print("RESULTS")
    print(f"{'='*40}")
    
    print(f"Found {len(solutions)} valid solutions:")
    
    for i, sol in enumerate(solutions):
        print(f"\nSolution {i+1}:")
        print(f"  x = {sol['x']:.6f} rad ({math.degrees(sol['x']):.1f}°)")
        print(f"  y = {sol['y']:.6f} rad ({math.degrees(sol['y']):.1f}°)")
        print(f"  t = {sol['t']:.6f}")
        
        # Check accuracy against known solution
        x_error = abs(sol['x'] - x_true)
        y_error = abs(sol['y'] - y_true)
        
        # Account for 2π periodicity
        x_error = min(x_error, abs(x_error - 2*math.pi), abs(x_error + 2*math.pi))
        y_error = min(y_error, abs(y_error - 2*math.pi), abs(y_error + 2*math.pi))
        
        print(f"  Error in x: {x_error:.2e}")
        print(f"  Error in y: {y_error:.2e}")
        
        # Verify the original equations
        eq1_check = A[0,0]*sol['cos_x'] + A[0,1]*sol['sin_x'] + B[0,0]*sol['cos_y'] + B[0,1]*sol['sin_y']
        eq2_check = A[1,0]*sol['cos_x'] + A[1,1]*sol['sin_x'] + B[1,0]*sol['cos_y'] + B[1,1]*sol['sin_y']
        
        print(f"  Equation 1: {eq1_check:.6f} = {C[0]:.6f} (residual: {abs(eq1_check - C[0]):.2e})")
        print(f"  Equation 2: {eq2_check:.6f} = {C[1]:.6f} (residual: {abs(eq2_check - C[1]):.2e})")
        
        if x_error < 1e-6 and y_error < 1e-6:
            print("  ✓ MATCHES KNOWN SOLUTION")
        else:
            print("  ? Different solution (may be valid due to periodicity)")


def test_multiple_cases():
    """Test with multiple random cases."""
    print(f"\n{'='*60}")
    print("TESTING WITH MULTIPLE RANDOM CASES")
    print(f"{'='*60}")
    
    np.random.seed(42)  # For reproducible results
    
    success_count = 0
    total_tests = 10
    
    for test_num in range(total_tests):
        print(f"\nTest {test_num + 1}/{total_tests}:")
        
        # Generate random test case
        x_true = np.random.uniform(-math.pi, math.pi)
        y_true = np.random.uniform(-math.pi, math.pi)
        
        cos_x_true = math.cos(x_true)
        sin_x_true = math.sin(x_true)
        cos_y_true = math.cos(y_true)
        sin_y_true = math.sin(y_true)
        
        # Generate random matrices (ensure B is well-conditioned)
        A = np.random.uniform(-2, 2, (2, 2))
        B = np.random.uniform(-2, 2, (2, 2))
        
        # Ensure B is non-singular
        while abs(np.linalg.det(B)) < 0.1:
            B = np.random.uniform(-2, 2, (2, 2))
        
        # Compute C from the known solution
        C = A @ np.array([cos_x_true, sin_x_true]) + B @ np.array([cos_y_true, sin_y_true])
        
        print(f"  True: x = {x_true:.3f}, y = {y_true:.3f}")
        
        # Solve the system
        solutions = solve_trigonometric_system(A, B, C, verbose=False)
        
        # Check if any solution matches the true solution
        found_match = False
        for sol in solutions:
            x_error = min(abs(sol['x'] - x_true), 
                         abs(sol['x'] - x_true - 2*math.pi),
                         abs(sol['x'] - x_true + 2*math.pi))
            y_error = min(abs(sol['y'] - y_true),
                         abs(sol['y'] - y_true - 2*math.pi), 
                         abs(sol['y'] - y_true + 2*math.pi))
            
            if x_error < 1e-6 and y_error < 1e-6:
                found_match = True
                break
        
        if found_match:
            success_count += 1
            print(f"  ✓ Found matching solution among {len(solutions)} solutions")
        else:
            print(f"  ✗ No matching solution found among {len(solutions)} solutions")
    
    print(f"\nSuccess rate: {success_count}/{total_tests} ({100*success_count/total_tests:.1f}%)")


if __name__ == "__main__":
    test_solver()
    test_multiple_cases()
