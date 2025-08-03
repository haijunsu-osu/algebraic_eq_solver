#!/usr/bin/env python3
"""
Symbolic Quartic Polynomial Coefficients

This script outputs the exact symbolic expressions for the coefficients
of the quartic polynomial in variable t, derived from the trigonometric system:

A[cos x, sin x] + B[cos y, sin y] = C

Using Weierstrass substitution: cos(x) = (1-t²)/(1+t²), sin(x) = 2t/(1+t²)

Author: AI Assistant
Date: August 3, 2025
"""

import sympy as sp
from sympy import symbols, simplify, latex


def display_quartic_coefficients():
    """Display the symbolic expressions for quartic polynomial coefficients."""
    
    print("=" * 100)
    print("SYMBOLIC EXPRESSIONS FOR QUARTIC POLYNOMIAL COEFFICIENTS")
    print("=" * 100)
    
    print("\nSystem: A[cos x, sin x] + B[cos y, sin y] = C")
    print("Where:")
    print("  A = [[A₀₀, A₀₁], [A₁₀, A₁₁]]")
    print("  B = [[B₀₀, B₀₁], [B₁₀, B₁₁]]")
    print("  C = [C₀, C₁]")
    
    print("\nAfter applying:")
    print("1. Linear system solution: B[cos y, sin y]ᵀ = C - A[cos x, sin x]ᵀ")
    print("2. Trigonometric identity: cos²(y) + sin²(y) = 1")
    print("3. Weierstrass substitution: cos(x) = (1-t²)/(1+t²), sin(x) = 2t/(1+t²)")
    
    print("\nThe resulting quartic polynomial is: a₄t⁴ + a₃t³ + a₂t² + a₁t + a₀ = 0")
    
    # Define symbolic variables
    A00, A01, A10, A11 = symbols('A00 A01 A10 A11')
    B00, B01, B10, B11 = symbols('B00 B01 B10 B11')
    C0, C1 = symbols('C0 C1')
    
    # Common denominator
    det_B = B00*B11 - B01*B10
    denom = det_B**2
    
    print(f"\nCommon denominator: Δ² where Δ = det(B) = B₀₀B₁₁ - B₀₁B₁₀")
    
    # Define coefficient numerators
    print("\n" + "="*50)
    print("COEFFICIENT a₄ (t⁴ term)")
    print("="*50)
    
    a4_num = (A00**2*B10**2 + A00**2*B11**2 - 2*A00*A10*B00*B10 - 2*A00*A10*B01*B11 
              - 2*A00*B00*B10*C1 - 2*A00*B01*B11*C1 + 2*A00*B10**2*C0 + 2*A00*B11**2*C0 
              + A10**2*B00**2 + A10**2*B01**2 + 2*A10*B00**2*C1 - 2*A10*B00*B10*C0 
              + 2*A10*B01**2*C1 - 2*A10*B01*B11*C0 - B00**2*B11**2 + B00**2*C1**2 
              + 2*B00*B01*B10*B11 - 2*B00*B10*C0*C1 - B01**2*B10**2 + B01**2*C1**2 
              - 2*B01*B11*C0*C1 + B10**2*C0**2 + B11**2*C0**2)
    
    a4 = a4_num / denom
    print(f"a₄ = ({simplify(a4_num)}) / Δ²")
    
    print("\n" + "="*50) 
    print("COEFFICIENT a₃ (t³ term)")
    print("="*50)
    
    a3_num = 4*(-A00*A01*B10**2 - A00*A01*B11**2 + A00*A11*B00*B10 + A00*A11*B01*B11 
                + A01*A10*B00*B10 + A01*A10*B01*B11 + A01*B00*B10*C1 + A01*B01*B11*C1 
                - A01*B10**2*C0 - A01*B11**2*C0 - A10*A11*B00**2 - A10*A11*B01**2 
                - A11*B00**2*C1 + A11*B00*B10*C0 - A11*B01**2*C1 + A11*B01*B11*C0)
    
    a3 = a3_num / denom
    print(f"a₃ = ({simplify(a3_num)}) / Δ²")
    
    print("\n" + "="*50)
    print("COEFFICIENT a₂ (t² term)")
    print("="*50)
    
    a2_num = 2*(-A00**2*B10**2 - A00**2*B11**2 + 2*A00*A10*B00*B10 + 2*A00*A10*B01*B11 
                + 2*A01**2*B10**2 + 2*A01**2*B11**2 - 4*A01*A11*B00*B10 - 4*A01*A11*B01*B11 
                - A10**2*B00**2 - A10**2*B01**2 + 2*A11**2*B00**2 + 2*A11**2*B01**2 
                - B00**2*B11**2 + B00**2*C1**2 + 2*B00*B01*B10*B11 - 2*B00*B10*C0*C1 
                - B01**2*B10**2 + B01**2*C1**2 - 2*B01*B11*C0*C1 + B10**2*C0**2 + B11**2*C0**2)
    
    a2 = a2_num / denom
    print(f"a₂ = ({simplify(a2_num)}) / Δ²")
    
    print("\n" + "="*50)
    print("COEFFICIENT a₁ (t¹ term)")
    print("="*50)
    
    a1_num = 4*(A00*A01*B10**2 + A00*A01*B11**2 - A00*A11*B00*B10 - A00*A11*B01*B11 
                - A01*A10*B00*B10 - A01*A10*B01*B11 + A01*B00*B10*C1 + A01*B01*B11*C1 
                - A01*B10**2*C0 - A01*B11**2*C0 + A10*A11*B00**2 + A10*A11*B01**2 
                - A11*B00**2*C1 + A11*B00*B10*C0 - A11*B01**2*C1 + A11*B01*B11*C0)
    
    a1 = a1_num / denom
    print(f"a₁ = ({simplify(a1_num)}) / Δ²")
    
    print("\n" + "="*50)
    print("COEFFICIENT a₀ (constant term)")
    print("="*50)
    
    a0_num = (A00**2*B10**2 + A00**2*B11**2 - 2*A00*A10*B00*B10 - 2*A00*A10*B01*B11 
              + 2*A00*B00*B10*C1 + 2*A00*B01*B11*C1 - 2*A00*B10**2*C0 - 2*A00*B11**2*C0 
              + A10**2*B00**2 + A10**2*B01**2 - 2*A10*B00**2*C1 + 2*A10*B00*B10*C0 
              - 2*A10*B01**2*C1 + 2*A10*B01*B11*C0 - B00**2*B11**2 + B00**2*C1**2 
              + 2*B00*B01*B10*B11 - 2*B00*B10*C0*C1 - B01**2*B10**2 + B01**2*C1**2 
              - 2*B01*B11*C0*C1 + B10**2*C0**2 + B11**2*C0**2)
    
    a0 = a0_num / denom
    print(f"a₀ = ({simplify(a0_num)}) / Δ²")
    
    print("\n" + "="*100)
    print("PATTERN ANALYSIS")
    print("="*100)
    
    print("\nObservations:")
    print("1. All coefficients have the same denominator: Δ² = (det(B))²")
    print("2. The coefficients are symmetric: a₄ and a₀ have similar structure")
    print("3. The coefficients a₃ and a₁ have similar structure with factor 4")
    print("4. The coefficient a₂ has factor 2 and includes cross-terms")
    
    print("\n" + "="*100)
    print("FACTORED FORMS (Common Terms)")
    print("="*100)
    
    # Group common terms
    print("\nCommon quadratic forms appearing in coefficients:")
    print("• A-terms: A₀₀², A₀₁², A₁₀², A₁₁², A₀₀A₁₀, A₀₁A₁₁")
    print("• B-terms: B₀₀², B₀₁², B₁₀², B₁₁², B₀₀B₁₁, B₀₁B₁₀") 
    print("• C-terms: C₀², C₁², C₀C₁")
    print("• Mixed AB-terms: AᵢⱼBₖₗ products")
    print("• Mixed AC-terms: AᵢⱼCₖ products")
    print("• Mixed BC-terms: BᵢⱼCₖ products")
    
    return {
        'a4': a4,
        'a3': a3,
        'a2': a2,
        'a1': a1,
        'a0': a0,
        'det_B': det_B
    }


def generate_python_implementation():
    """Generate optimized Python implementation."""
    
    print("\n" + "="*100)
    print("OPTIMIZED PYTHON IMPLEMENTATION")
    print("="*100)
    
    implementation = '''
def compute_quartic_coefficients(A, B, C):
    """
    Compute quartic polynomial coefficients for trigonometric system.
    
    Args:
        A: 2x2 matrix
        B: 2x2 matrix 
        C: 2x1 vector
    
    Returns:
        [a4, a3, a2, a1, a0] coefficients of quartic polynomial
    """
    # Extract matrix elements
    A00, A01, A10, A11 = A[0,0], A[0,1], A[1,0], A[1,1]
    B00, B01, B10, B11 = B[0,0], B[0,1], B[1,0], B[1,1]
    C0, C1 = C[0], C[1]
    
    # Pre-compute common terms for efficiency
    A00_sq, A01_sq, A10_sq, A11_sq = A00**2, A01**2, A10**2, A11**2
    B00_sq, B01_sq, B10_sq, B11_sq = B00**2, B01**2, B10**2, B11**2
    C0_sq, C1_sq = C0**2, C1**2
    
    # Cross products
    A00_A10 = A00*A10
    A01_A11 = A01*A11
    B00_B11 = B00*B11
    B01_B10 = B01*B10
    C0_C1 = C0*C1
    
    # Determinant and its square
    det_B = B00_B11 - B01_B10
    det_B_sq = det_B**2
    
    # Coefficient numerators (optimized with pre-computed terms)
    a4_num = (A00_sq*(B10_sq + B11_sq) - 2*A00_A10*(B00*B10 + B01*B11)
              - 2*A00*(B00*B10 + B01*B11)*C1 + 2*A00*(B10_sq + B11_sq)*C0
              + A10_sq*(B00_sq + B01_sq) + 2*A10*(B00_sq + B01_sq)*C1
              - 2*A10*(B00*B10 + B01*B11)*C0 - B00_sq*B11_sq + B00_sq*C1_sq
              + 2*B00_B11*B01_B10 - 2*B00*B10*C0_C1 - B01_sq*B10_sq
              + B01_sq*C1_sq - 2*B01*B11*C0_C1 + (B10_sq + B11_sq)*C0_sq)
    
    a3_num = 4*(-A00*A01*(B10_sq + B11_sq) + A00*A11*(B00*B10 + B01*B11)
                + A01*A10*(B00*B10 + B01*B11) + A01*(B00*B10 + B01*B11)*C1
                - A01*(B10_sq + B11_sq)*C0 - A10*A11*(B00_sq + B01_sq)
                - A11*(B00_sq + B01_sq)*C1 + A11*(B00*B10 + B01*B11)*C0)
    
    a2_num = 2*(-A00_sq*(B10_sq + B11_sq) + 2*A00_A10*(B00*B10 + B01*B11)
                + 2*A01_sq*(B10_sq + B11_sq) - 4*A01_A11*(B00*B10 + B01*B11)
                - A10_sq*(B00_sq + B01_sq) + 2*A11_sq*(B00_sq + B01_sq)
                - B00_sq*B11_sq + B00_sq*C1_sq + 2*B00_B11*B01_B10
                - 2*B00*B10*C0_C1 - B01_sq*B10_sq + B01_sq*C1_sq
                - 2*B01*B11*C0_C1 + (B10_sq + B11_sq)*C0_sq)
    
    a1_num = 4*(A00*A01*(B10_sq + B11_sq) - A00*A11*(B00*B10 + B01*B11)
                - A01*A10*(B00*B10 + B01*B11) + A01*(B00*B10 + B01*B11)*C1
                - A01*(B10_sq + B11_sq)*C0 + A10*A11*(B00_sq + B01_sq)
                - A11*(B00_sq + B01_sq)*C1 + A11*(B00*B10 + B01*B11)*C0)
    
    a0_num = (A00_sq*(B10_sq + B11_sq) - 2*A00_A10*(B00*B10 + B01*B11)
              + 2*A00*(B00*B10 + B01*B11)*C1 - 2*A00*(B10_sq + B11_sq)*C0
              + A10_sq*(B00_sq + B01_sq) - 2*A10*(B00_sq + B01_sq)*C1
              + 2*A10*(B00*B10 + B01*B11)*C0 - B00_sq*B11_sq + B00_sq*C1_sq
              + 2*B00_B11*B01_B10 - 2*B00*B10*C0_C1 - B01_sq*B10_sq
              + B01_sq*C1_sq - 2*B01*B11*C0_C1 + (B10_sq + B11_sq)*C0_sq)
    
    # Return coefficients (highest degree first)
    return [a4_num/det_B_sq, a3_num/det_B_sq, a2_num/det_B_sq, 
            a1_num/det_B_sq, a0_num/det_B_sq]
'''
    
    print(implementation)


def main():
    """Main function to display all symbolic expressions."""
    coefficients = display_quartic_coefficients()
    generate_python_implementation()
    
    print("\n" + "="*100)
    print("VERIFICATION")
    print("="*100)
    print("These symbolic expressions can be verified by:")
    print("1. Substituting specific numerical values for A, B, C")
    print("2. Computing the coefficients using both symbolic and numerical methods")
    print("3. Solving the quartic polynomial and verifying solutions satisfy original equations")
    print("4. Checking that cos²(x) + sin²(x) = 1 and cos²(y) + sin²(y) = 1")


if __name__ == "__main__":
    main()
