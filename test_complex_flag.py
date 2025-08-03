#!/usr/bin/env python3
"""
Demonstration of the real_solutions_only flag in solve_trigonometric_system()

This script shows how to use the new flag to control whether complex solutions
are included in the results.
"""

import numpy as np
import math
from numerical_solver import solve_trigonometric_system

def test_real_solutions_only_flag():
    """Demonstrate the real_solutions_only flag functionality."""
    
    print("=" * 70)
    print("TESTING real_solutions_only FLAG")
    print("=" * 70)
    
    # Set up a test system
    A = np.array([[2.0, 1.0],
                  [1.0, 3.0]])
    
    B = np.array([[1.5, -0.5],
                  [0.5, 2.0]])
    
    C = np.array([3.0, 4.0])
    
    print(f"Test system:")
    print(f"A = \n{A}")
    print(f"B = \n{B}")
    print(f"C = {C}")
    print()
    
    # Test 1: Default behavior (real solutions only)
    print("1. DEFAULT: real_solutions_only=True")
    print("-" * 40)
    solutions_real = solve_trigonometric_system(A, B, C, verbose=False, real_solutions_only=True)
    print(f"Found {len(solutions_real)} real solutions:")
    
    for i, sol in enumerate(solutions_real):
        is_complex = sol.get('is_complex', False)
        print(f"  Solution {i+1}: x = {sol['x']:.6f} rad, y = {sol['y']:.6f} rad")
        print(f"               is_complex = {is_complex}")
    print()
    
    # Test 2: Include complex solutions
    print("2. EXTENDED: real_solutions_only=False")
    print("-" * 40)
    solutions_all = solve_trigonometric_system(A, B, C, verbose=False, real_solutions_only=False)
    print(f"Found {len(solutions_all)} total solutions:")
    
    for i, sol in enumerate(solutions_all):
        is_complex = sol.get('is_complex', False)
        if is_complex:
            print(f"  Solution {i+1}: COMPLEX SOLUTION")
            print(f"               x = {sol['x']}")
            print(f"               y = {sol['y']}")
        else:
            print(f"  Solution {i+1}: x = {sol['x']:.6f} rad, y = {sol['y']:.6f} rad")
        print(f"               is_complex = {is_complex}")
    print()
    
    # Summary
    additional_solutions = len(solutions_all) - len(solutions_real)
    print("SUMMARY:")
    print(f"  Real solutions: {len(solutions_real)}")
    print(f"  Total solutions: {len(solutions_all)}")
    print(f"  Additional complex solutions: {additional_solutions}")
    
    if additional_solutions > 0:
        print(f"\n✓ The real_solutions_only flag successfully filtered out {additional_solutions} complex solution(s)")
    else:
        print(f"\n• No complex solutions found for this particular system")
        print("  (Try different A, B, C values to see complex solutions)")

if __name__ == "__main__":
    test_real_solutions_only_flag()
