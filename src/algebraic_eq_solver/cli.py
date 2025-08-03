#!/usr/bin/env python3
"""
Command Line Interface for Algebraic Equation Solver
"""

import argparse
import numpy as np
import sys
from .numerical_solver import solve_trig_sys, test_solver, test_multiple_cases, test_singular_cases


def parse_matrix(matrix_str: str) -> np.ndarray:
    """Parse a matrix string in the format '[[1,2],[3,4]]'"""
    try:
        # Remove spaces and split by rows
        matrix_str = matrix_str.strip()
        if not matrix_str.startswith('[[') or not matrix_str.endswith(']]'):
            raise ValueError("Matrix must be in format [[a,b],[c,d]]")
        
        # Parse the matrix
        matrix_str = matrix_str[2:-2]  # Remove outer brackets
        rows = matrix_str.split('],[')
        
        matrix_data = []
        for row in rows:
            row_data = [float(x.strip()) for x in row.split(',')]
            matrix_data.append(row_data)
        
        return np.array(matrix_data)
    except Exception as e:
        raise ValueError(f"Invalid matrix format: {e}")


def parse_vector(vector_str: str) -> np.ndarray:
    """Parse a vector string in the format '[1,2]'"""
    try:
        vector_str = vector_str.strip()
        if not vector_str.startswith('[') or not vector_str.endswith(']'):
            raise ValueError("Vector must be in format [a,b]")
        
        vector_str = vector_str[1:-1]  # Remove brackets
        vector_data = [float(x.strip()) for x in vector_str.split(',')]
        
        return np.array(vector_data)
    except Exception as e:
        raise ValueError(f"Invalid vector format: {e}")


def solve_command(args):
    """Handle the solve command"""
    try:
        A = parse_matrix(args.A)
        B = parse_matrix(args.B)
        C = parse_vector(args.C)
        
        print("Input System:")
        print(f"A = \n{A}")
        print(f"B = \n{B}")
        print(f"C = {C}")
        print()
        
        solutions = solve_trig_sys(A, B, C, verbose=args.verbose, real_solutions_only=not args.complex)
        
        print(f"Found {len(solutions)} solutions:")
        print("-" * 50)
        
        for i, sol in enumerate(solutions):
            print(f"Solution {i+1}:")
            if sol.get('is_complex', False):
                print(f"  θ₁ = {sol['th1']} rad (complex)")
                print(f"  θ₂ = {sol['th2']} rad (complex)")
            else:
                print(f"  θ₁ = {sol['th1']:.6f} rad ({np.degrees(sol['th1']):.2f}°)")
                print(f"  θ₂ = {sol['th2']:.6f} rad ({np.degrees(sol['th2']):.2f}°)")
                
                # Verify solution
                residual1 = A[0,0]*sol['cos_th1'] + A[0,1]*sol['sin_th1'] + B[0,0]*sol['cos_th2'] + B[0,1]*sol['sin_th2'] - C[0]
                residual2 = A[1,0]*sol['cos_th1'] + A[1,1]*sol['sin_th1'] + B[1,0]*sol['cos_th2'] + B[1,1]*sol['sin_th2'] - C[1]
                print(f"  Residuals: {residual1:.2e}, {residual2:.2e}")
            print()
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    return 0


def test_command(args):
    """Handle the test command"""
    if args.type == 'basic':
        test_solver()
    elif args.type == 'random':
        test_multiple_cases()
    elif args.type == 'singular':
        test_singular_cases()
    elif args.type == 'all':
        print("Running all tests...\n")
        test_solver()
        test_multiple_cases()
        test_singular_cases()
    
    return 0


def main():
    """Main entry point for the CLI"""
    parser = argparse.ArgumentParser(
        description="Algebraic Equation Solver - Solve trigonometric systems A[cos θ₁, sin θ₁] + B[cos θ₂, sin θ₂] = C",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Solve a simple system
  trig-solver solve -A "[[2,1],[1,3]]" -B "[[1.5,-0.5],[0.5,2]]" -C "[1,2]"
  
  # Solve with verbose output
  trig-solver solve -A "[[1,0],[0,1]]" -B "[[2,2],[1,1]]" -C "[1,0.5]" --verbose
  
  # Run tests
  trig-solver test --type all
  
  # Include complex solutions
  trig-solver solve -A "[[1,0],[0,1]]" -B "[[1,1],[1,1]]" -C "[1,1]" --complex
        """)
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Solve command
    solve_parser = subparsers.add_parser('solve', help='Solve a trigonometric system')
    solve_parser.add_argument('-A', required=True, help='Matrix A in format [[a,b],[c,d]]')
    solve_parser.add_argument('-B', required=True, help='Matrix B in format [[a,b],[c,d]]')
    solve_parser.add_argument('-C', required=True, help='Vector C in format [a,b]')
    solve_parser.add_argument('--verbose', '-v', action='store_true', help='Show detailed solution steps')
    solve_parser.add_argument('--complex', action='store_true', help='Include complex solutions')
    solve_parser.set_defaults(func=solve_command)
    
    # Test command
    test_parser = subparsers.add_parser('test', help='Run test suites')
    test_parser.add_argument('--type', choices=['basic', 'random', 'singular', 'all'], 
                           default='basic', help='Type of test to run')
    test_parser.set_defaults(func=test_command)
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
