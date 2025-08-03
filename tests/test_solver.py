"""
Test suite for algebraic_eq_solver package
"""

import pytest
import numpy as np
import math
from algebraic_eq_solver import solve_trig_sys, analyze_singular_matrix


class TestBasicSolver:
    """Test basic solver functionality"""
    
    def test_simple_case(self):
        """Test a simple known case"""
        # Create test case with known solution
        th1_true = math.pi / 6  # 30 degrees
        th2_true = math.pi / 4  # 45 degrees
        
        cos_th1_true = math.cos(th1_true)
        sin_th1_true = math.sin(th1_true)
        cos_th2_true = math.cos(th2_true)
        sin_th2_true = math.sin(th2_true)
        
        A = np.array([[2.0, 1.0], [1.0, 3.0]])
        B = np.array([[1.5, -0.5], [0.5, 2.0]])
        C = A @ np.array([cos_th1_true, sin_th1_true]) + B @ np.array([cos_th2_true, sin_th2_true])
        
        solutions = solve_trig_sys(A, B, C, verbose=False)
        
        # Should find at least one solution
        assert len(solutions) > 0
        
        # At least one solution should match the original
        found_match = False
        for sol in solutions:
            th1_error = min(abs(sol['th1'] - th1_true), 
                          abs(sol['th1'] - th1_true - 2*math.pi),
                          abs(sol['th1'] - th1_true + 2*math.pi))
            th2_error = min(abs(sol['th2'] - th2_true),
                          abs(sol['th2'] - th2_true - 2*math.pi), 
                          abs(sol['th2'] - th2_true + 2*math.pi))
            
            if th1_error < 1e-6 and th2_error < 1e-6:
                found_match = True
                break
        
        assert found_match, "Should find the original solution"
    
    def test_input_validation(self):
        """Test input validation"""
        A = np.array([[1, 0], [0, 1]])
        B = np.array([[1, 0], [0, 1]])
        
        # Wrong dimensions should raise ValueError
        with pytest.raises(ValueError):
            solve_trig_sys(A, B, np.array([1, 2, 3]))  # Wrong C dimension
        
        with pytest.raises(ValueError):
            solve_trig_sys(np.array([[1, 2, 3]]), B, np.array([1, 2]))  # Wrong A dimension


class TestSingularCases:
    """Test singular matrix handling"""
    
    def test_zero_matrix(self):
        """Test B = 0 case"""
        A = np.array([[1.0, 0.5], [0.3, 1.2]])
        B = np.array([[0.0, 0.0], [0.0, 0.0]])
        
        # Create a valid C
        cos_th1_true = 0.8
        sin_th1_true = 0.6  # cos²+sin²=1
        C = A @ np.array([cos_th1_true, sin_th1_true])
        
        solutions = solve_trig_sys(A, B, C, verbose=False)
        
        # Should find solutions (θ₂ is free parameter)
        assert len(solutions) > 0
        
        # All solutions should have the same θ₁
        th1_values = [sol['th1'] for sol in solutions]
        th1_expected = math.atan2(sin_th1_true, cos_th1_true)
        
        for th1 in th1_values:
            th1_error = min(abs(th1 - th1_expected),
                          abs(th1 - th1_expected - 2*math.pi),
                          abs(th1 - th1_expected + 2*math.pi))
            assert th1_error < 1e-6, "All solutions should have the same θ₁"
    
    def test_rank_1_matrix(self):
        """Test rank-1 B matrix"""
        A = np.array([[1.0, 0.0], [0.0, 1.0]])
        B = np.array([[2.0, 2.0], [1.0, 1.0]])  # rank 1
        C = np.array([1.0, -0.17082039])  # Chosen to give 4 solutions
        
        solutions = solve_trig_sys(A, B, C, verbose=False)
        
        # This specific case should have multiple solutions
        assert len(solutions) > 0
        
        # Verify all solutions satisfy the original equations
        for sol in solutions:
            residual1 = A[0,0]*sol['cos_th1'] + A[0,1]*sol['sin_th1'] + B[0,0]*sol['cos_th2'] + B[0,1]*sol['sin_th2'] - C[0]
            residual2 = A[1,0]*sol['cos_th1'] + A[1,1]*sol['sin_th1'] + B[1,0]*sol['cos_th2'] + B[1,1]*sol['sin_th2'] - C[1]
            
            assert abs(residual1) < 1e-10, f"Solution should satisfy equation 1, residual: {residual1}"
            assert abs(residual2) < 1e-10, f"Solution should satisfy equation 2, residual: {residual2}"


class TestAnalysis:
    """Test matrix analysis functions"""
    
    def test_singular_matrix_analysis(self):
        """Test singular matrix analysis"""
        # Zero matrix
        B_zero = np.array([[0, 0], [0, 0]])
        analysis = analyze_singular_matrix(B_zero)
        assert analysis['type'] == 'zero_matrix'
        assert analysis['rank'] == 0
        
        # Rank-1 matrix
        B_rank1 = np.array([[2, 4], [1, 2]])
        analysis = analyze_singular_matrix(B_rank1)
        assert analysis['type'] == 'rank_1'
        assert analysis['rank'] == 1
        
    def test_solution_format(self):
        """Test that solutions have the correct format"""
        A = np.array([[1, 0], [0, 1]])
        B = np.array([[1, 0], [0, 1]])
        C = np.array([1, 1])
        
        solutions = solve_trig_sys(A, B, C, verbose=False)
        
        for sol in solutions:
            # Check required keys
            required_keys = ['th1', 'th2', 'cos_th1', 'sin_th1', 'cos_th2', 'sin_th2']
            for key in required_keys:
                assert key in sol, f"Solution missing key: {key}"
            
            # Check trigonometric identities
            assert abs(sol['cos_th1']**2 + sol['sin_th1']**2 - 1) < 1e-10, "cos²+sin²=1 should hold for θ₁"
            assert abs(sol['cos_th2']**2 + sol['sin_th2']**2 - 1) < 1e-10, "cos²+sin²=1 should hold for θ₂"


if __name__ == '__main__':
    pytest.main([__file__])
