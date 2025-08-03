"""
Algebraic Equation Solver

A robust numerical solver for trigonometric systems of the form:
A[cos θ₁, sin θ₁] + B[cos θ₂, sin θ₂] = C

This package provides comprehensive solutions for both regular and singular
matrix cases with advanced numerical stability features.
"""

from .numerical_solver import (
    solve_trig_sys,
    solve_trig_eq,
    analyze_singular_matrix,
    solve_singular_b_system,
    solve_zero_b_case,
    solve_rank_1_b_case,
    solve_both_singular_case,
    test_solver,
    test_multiple_cases,
    test_singular_cases
)

__version__ = "1.0.0"
__author__ = "AI Assistant"
__email__ = ""

__all__ = [
    # Main solving function
    "solve_trig_sys",
    
    # Utility functions
    "solve_trig_eq",
    "analyze_singular_matrix",
    
    # Specialized solvers
    "solve_singular_b_system",
    "solve_zero_b_case", 
    "solve_rank_1_b_case",
    "solve_both_singular_case",
    
    # Testing functions
    "test_solver",
    "test_multiple_cases", 
    "test_singular_cases",
]
