#!/usr/bin/env python3
"""
Symbolic Derivation of Quartic Polynomial Coefficients

This script derives the symbolic expressions for the coefficients of the quartic 
polynomial obtained by solving the system:
A[cos x, sin x] + B[cos y, sin y] = C

Where A, B are 2x2 matrices and C is a 2x1 vector.

The    solutions = solve_trig_sys(A, B, C)
    
    print(f"Found {len(solutions)} solutions:")proach:
1. Express cos(y) and sin(y) in terms of cos(x) and sin(x)
2. Apply the identity cos²(y) + sin²(y) = 1
3. Use Weierstrass substitution: cos(x) = (1-t²)/(1+t²), sin(x) = 2t/(1+t²)
4. Derive the quartic polynomial in t

Author: AI Assistant
Date: August 3, 2025
"""

import sympy as sp
from sympy import symbols, cos, sin, simplify, expand, collect, solve
import numpy as np

def derive_quartic_coefficients():
    """
    Derive symbolic expressions for quartic polynomial coefficients.
    
    Returns:
        Dictionary containing symbolic expressions for coefficients a4, a3, a2, a1, a0
    """
    print("Starting symbolic derivation...")
    
    # Define symbolic variables
    x, y, t = symbols('x y t')
    
    # Define symbolic matrix elements
    # A = [[A00, A01], [A10, A11]]
    A00, A01, A10, A11 = symbols('A00 A01 A10 A11')
    # B = [[B00, B01], [B10, B11]]  
    B00, B01, B10, B11 = symbols('B00 B01 B10 B11')
    # C = [C0, C1]
    C0, C1 = symbols('C0 C1')
    
    print("Step 1: Setting up the original system of equations")
    # Original system: A[cos x, sin x] + B[cos y, sin y] = C
    # Equation 1: A00*cos(x) + A01*sin(x) + B00*cos(y) + B01*sin(y) = C0
    # Equation 2: A10*cos(x) + A11*sin(x) + B10*cos(y) + B11*sin(y) = C1
    
    eq1 = A00*cos(x) + A01*sin(x) + B00*cos(y) + B01*sin(y) - C0
    eq2 = A10*cos(x) + A11*sin(x) + B10*cos(y) + B11*sin(y) - C1
    
    print("Step 2: Expressing cos(y) and sin(y) in terms of cos(x) and sin(x)")
    
    # Solve for cos(y) and sin(y) from the linear system
    # This gives us a 2x2 linear system in cos(y) and sin(y)
    cos_y, sin_y = symbols('cos_y sin_y')
    
    # Rewrite equations in matrix form: B * [cos_y, sin_y]^T = C - A * [cos_x, sin_x]^T
    rhs1 = C0 - A00*cos(x) - A01*sin(x)
    rhs2 = C1 - A10*cos(x) - A11*sin(x)
    
    # Solve the 2x2 system using Cramer's rule
    det_B = B00*B11 - B01*B10
    
    cos_y_expr = (B11*rhs1 - B01*rhs2) / det_B
    sin_y_expr = (-B10*rhs1 + B00*rhs2) / det_B
    
    print("cos(y) =", cos_y_expr)
    print("sin(y) =", sin_y_expr)
    
    print("Step 3: Applying the trigonometric identity cos²(y) + sin²(y) = 1")
    
    # Apply the identity
    identity_eq = cos_y_expr**2 + sin_y_expr**2 - 1
    identity_eq = simplify(identity_eq)
    
    print("Step 4: Applying Weierstrass substitution")
    print("cos(x) = (1-t²)/(1+t²), sin(x) = 2t/(1+t²)")
    
    # Weierstrass substitution
    cos_x_sub = (1 - t**2) / (1 + t**2)
    sin_x_sub = 2*t / (1 + t**2)
    
    # Substitute into the identity equation
    final_eq = identity_eq.subs([(cos(x), cos_x_sub), (sin(x), sin_x_sub)])
    
    print("Step 5: Simplifying to get quartic polynomial")
    
    # Clear denominators by multiplying by (1+t²)²
    final_eq = final_eq * (1 + t**2)**2
    final_eq = expand(final_eq)
    final_eq = simplify(final_eq)
    
    # Collect coefficients of powers of t
    poly_coeffs = sp.Poly(final_eq, t).all_coeffs()
    
    # Ensure we have a quartic (degree 4) polynomial
    while len(poly_coeffs) < 5:
        poly_coeffs.insert(0, 0)
    
    if len(poly_coeffs) > 5:
        print("Warning: Polynomial degree higher than 4!")
    
    # Extract coefficients (highest degree first)
    a4 = poly_coeffs[0] if len(poly_coeffs) > 4 else 0
    a3 = poly_coeffs[1] if len(poly_coeffs) > 3 else 0
    a2 = poly_coeffs[2] if len(poly_coeffs) > 2 else 0
    a1 = poly_coeffs[3] if len(poly_coeffs) > 1 else 0
    a0 = poly_coeffs[4] if len(poly_coeffs) > 0 else 0
    
    print("\nQuartic polynomial: a4*t⁴ + a3*t³ + a2*t² + a1*t + a0 = 0")
    print("\nSymbolic coefficients:")
    print("a4 =", simplify(a4))
    print("a3 =", simplify(a3))
    print("a2 =", simplify(a2))
    print("a1 =", simplify(a1))
    print("a0 =", simplify(a0))
    
    return {
        'a4': simplify(a4),
        'a3': simplify(a3),
        'a2': simplify(a2),
        'a1': simplify(a1),
        'a0': simplify(a0)
    }

def generate_numerical_solver_code(coefficients):
    """
    Generate Python code for the numerical solver.
    
    Args:
        coefficients: Dictionary of symbolic coefficient expressions
    """
    print("\nGenerating numerical solver code...")
    
    # Convert SymPy expressions to Python code strings
    def sympy_to_python(expr):
        """Convert SymPy expression to Python code string."""
        # Convert the expression to a string and replace SymPy functions
        code = str(expr)
        code = code.replace('**', '**')
        return code
    
    a4_code = sympy_to_python(coefficients['a4'])
    a3_code = sympy_to_python(coefficients['a3'])
    a2_code = sympy_to_python(coefficients['a2'])
    a1_code = sympy_to_python(coefficients['a1'])
    a0_code = sympy_to_python(coefficients['a0'])
    
    solver_code = f'''#!/usr/bin/env python3
"""
Numerical Solver for Trigonometric System

Solves the system: A[cos x, sin x] + B[cos y, sin y] = C

This script uses the pre-derived quartic polynomial coefficients to solve
the system numerically without SymPy dependencies.

Dependencies: numpy, math only

Author: AI Assistant  
Date: August 3, 2025
"""

import numpy as np
import math


def solve_trig_sys(A, B, C):
    """
    Solve the trigonometric system A[cos th1, sin th1] + B[cos th2, sin th2] = C.
    
    Args:
        A: 2x2 numpy array
        B: 2x2 numpy array  
        C: 2x1 numpy array
        
    Returns:
        List of solution dictionaries, each containing 'th1', 'th2', 'cos_th1', 'sin_th1', 'cos_th2', 'sin_th2'
    """
    solutions = []
    
    # Extract matrix elements
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
    
    # Compute quartic polynomial coefficients using derived formulas
    a4 = {a4_code}
    a3 = {a3_code}
    a2 = {a2_code}
    a1 = {a1_code}
    a0 = {a0_code}
    
    # Form coefficient array for numpy.roots (highest degree first)
    coeffs = [a4, a3, a2, a1, a0]
    
    # Remove leading zeros
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-12:
        coeffs = coeffs[1:]
    
    if len(coeffs) == 1:
        if abs(coeffs[0]) < 1e-12:
            print("Warning: All coefficients are zero - infinite solutions")
            return solutions
        else:
            print("Warning: No solutions - constant polynomial")
            return solutions
    
    # Solve quartic polynomial for t
    try:
        t_roots = np.roots(coeffs)
    except np.linalg.LinAlgError:
        print("Warning: Failed to find polynomial roots")
        return solutions
    
    # Process each root
    for t in t_roots:
        # Skip complex roots with significant imaginary part
        if np.iscomplex(t) and abs(t.imag) > 1e-10:
            continue
            
        t = t.real if np.iscomplex(t) else t
        
        # Convert t back to cos(x) and sin(x) using inverse Weierstrass substitution
        denominator = 1 + t**2
        cos_x = (1 - t**2) / denominator
        sin_x = 2*t / denominator
        
        # Verify trigonometric identity
        if abs(cos_x**2 + sin_x**2 - 1) > 1e-10:
            continue
            
        # Calculate x from cos_x and sin_x
        x = math.atan2(sin_x, cos_x)
        
        # Solve for cos(y) and sin(y) using the original equations
        # B * [cos_y, sin_y]^T = C - A * [cos_x, sin_x]^T
        rhs = C - A @ np.array([cos_x, sin_x])
        
        try:
            trig_y = np.linalg.solve(B, rhs)
            cos_y, sin_y = trig_y[0], trig_y[1]
        except np.linalg.LinAlgError:
            continue
            
        # Verify trigonometric identity for y
        if abs(cos_y**2 + sin_y**2 - 1) > 1e-10:
            continue
            
        # Calculate y from cos_y and sin_y
        y = math.atan2(sin_y, cos_y)
        
        # Verify the original equations
        eq1_residual = A[0,0]*cos_x + A[0,1]*sin_x + B[0,0]*cos_y + B[0,1]*sin_y - C[0]
        eq2_residual = A[1,0]*cos_x + A[1,1]*sin_x + B[1,0]*cos_y + B[1,1]*sin_y - C[1]
        
        if abs(eq1_residual) < 1e-10 and abs(eq2_residual) < 1e-10:
            solutions.append({{
                'x': x,
                'y': y,
                'cos_x': cos_x,
                'sin_x': sin_x,
                'cos_y': cos_y,
                'sin_y': sin_y,
                't': t
            }})
    
    return solutions


def test_solver():
    """Test the solver with a known solution."""
    print("Testing the numerical solver...")
    
    # Create a test case with known solution
    x_true = math.pi / 6  # 30 degrees
    y_true = math.pi / 4  # 45 degrees
    
    cos_x_true = math.cos(x_true)
    sin_x_true = math.sin(x_true)
    cos_y_true = math.cos(y_true)
    sin_y_true = math.sin(y_true)
    
    print(f"True solution: x = {{x_true:.6f}} rad ({{math.degrees(x_true):.1f}}°)")
    print(f"               y = {{y_true:.6f}} rad ({{math.degrees(y_true):.1f}}°)")
    
    # Define test matrices
    A = np.array([[2.0, 1.0],
                  [1.0, 3.0]])
    
    B = np.array([[1.5, -0.5],
                  [0.5, 2.0]])
    
    # Compute C from the known solution
    C = A @ np.array([cos_x_true, sin_x_true]) + B @ np.array([cos_y_true, sin_y_true])
    
    print(f"\\nTest matrices:")
    print(f"A = \\n{{A}}")
    print(f"B = \\n{{B}}")
    print(f"C = {{C}}")
    
    # Solve the system
    solutions = solve_trigonometric_system(A, B, C)
    
    print(f"\\nFound {{len(solutions)}} solutions:")
    
    for i, sol in enumerate(solutions):
        print(f"\\nSolution {{i+1}}:")
        print(f"  x = {{sol['x']:.6f}} rad ({{math.degrees(sol['x']):.1f}}°)")
        print(f"  y = {{sol['y']:.6f}} rad ({{math.degrees(sol['y']):.1f}}°)")
        print(f"  t = {{sol['t']:.6f}}")
        
        # Check accuracy
        x_error = abs(sol['x'] - x_true)
        y_error = abs(sol['y'] - y_true)
        
        # Account for 2π periodicity
        x_error = min(x_error, abs(x_error - 2*math.pi), abs(x_error + 2*math.pi))
        y_error = min(y_error, abs(y_error - 2*math.pi), abs(y_error + 2*math.pi))
        
        print(f"  Error in x: {{x_error:.2e}}")
        print(f"  Error in y: {{y_error:.2e}}")
        
        # Verify the original equations
        eq1_check = A[0,0]*sol['cos_x'] + A[0,1]*sol['sin_x'] + B[0,0]*sol['cos_y'] + B[0,1]*sol['sin_y']
        eq2_check = A[1,0]*sol['cos_x'] + A[1,1]*sol['sin_x'] + B[1,0]*sol['cos_y'] + B[1,1]*sol['sin_y']
        
        print(f"  Equation 1 residual: {{abs(eq1_check - C[0]):.2e}}")
        print(f"  Equation 2 residual: {{abs(eq2_check - C[1]):.2e}}")


if __name__ == "__main__":
    test_solver()
'''
    
    return solver_code

def main():
    """Main function to run the symbolic derivation."""
    print("=" * 60)
    print("SYMBOLIC DERIVATION OF QUARTIC POLYNOMIAL COEFFICIENTS")
    print("=" * 60)
    
    # Derive the symbolic coefficients
    coefficients = derive_quartic_coefficients()
    
    # Generate the numerical solver
    solver_code = generate_numerical_solver_code(coefficients)
    
    # Write the numerical solver to file
    with open('numerical_solver.py', 'w', encoding='utf-8') as f:
        f.write(solver_code)
    
    print("\n" + "=" * 60)
    print("SYMBOLIC DERIVATION COMPLETE")
    print("=" * 60)
    print("Generated numerical_solver.py")
    print("Run: python numerical_solver.py to test the solver")

if __name__ == "__main__":
    main()
