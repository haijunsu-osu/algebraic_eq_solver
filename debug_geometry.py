#!/usr/bin/env python3
"""
Detailed debug for the rank-1 case quadratic solving
"""

import numpy as np
import math

def debug_rank1_geometry():
    """Debug the geometric constraint solving in detail"""
    
    # Use the exact same setup as before
    B = np.array([[1.0, 2.0], [0.5, 1.0]])
    A = np.array([[1.0, 0.0], [0.0, 1.0]])
    
    # Known solution
    x_known = math.pi / 6  # 30 degrees
    y_known = math.pi / 4  # 45 degrees
    cos_x = math.cos(x_known)
    sin_x = math.sin(x_known)
    cos_y = math.cos(y_known)
    sin_y = math.sin(y_known)
    
    C = A @ np.array([cos_x, sin_x]) + B @ np.array([cos_y, sin_y])
    
    # SVD analysis
    U, s, Vt = np.linalg.svd(B)
    range_direction = U[:, 0]
    domain_direction = Vt[0, :]
    
    print("=== DEBUGGING RANK-1 GEOMETRY ===")
    print(f"Known solution: x={x_known:.6f}, y={y_known:.6f}")
    print(f"cos(x)={cos_x:.6f}, sin(x)={sin_x:.6f}")
    print(f"cos(y)={cos_y:.6f}, sin(y)={sin_y:.6f}")
    print(f"C = {C}")
    print(f"Range direction: {range_direction}")
    print(f"Domain direction: {domain_direction}")
    print(f"Singular value: {s[0]}")
    
    # Test the known solution point
    x = x_known
    cos_x_test = math.cos(x)
    sin_x_test = math.sin(x)
    
    required_vector = C - A @ np.array([cos_x_test, sin_x_test])
    print(f"\nFor x = {x:.6f} (known solution):")
    print(f"Required vector: {required_vector}")
    
    # Check parallelism
    required_unit = required_vector / np.linalg.norm(required_vector)
    dot_product = abs(np.dot(required_unit, range_direction))
    print(f"Parallel check: {dot_product:.6f} (should be ~1)")
    
    if dot_product > 1 - 1e-8:
        print("✓ Required vector is in range of B")
        
        # Now solve the constraint
        required_magnitude = np.linalg.norm(required_vector)
        required_dot_product = required_magnitude / s[0]
        
        d0, d1 = domain_direction[0], domain_direction[1]
        k = required_dot_product
        
        print(f"\nConstraint equation:")
        print(f"d0*cos(y) + d1*sin(y) = k")
        print(f"{d0:.6f}*cos(y) + {d1:.6f}*sin(y) = {k:.6f}")
        print(f"cos²(y) + sin²(y) = 1")
        
        # Check if known solution satisfies this
        check_value = d0 * cos_y + d1 * sin_y
        print(f"\nKnown solution check: {d0:.6f}*{cos_y:.6f} + {d1:.6f}*{sin_y:.6f} = {check_value:.6f}")
        print(f"Should equal k = {k:.6f}")
        print(f"Difference: {abs(check_value - k):.2e}")
        
        if abs(d0) > 1e-12:
            # Solve quadratic
            a = d1**2 + d0**2
            b = -2 * k * d1
            c = k**2 - d0**2
            
            print(f"\nQuadratic: {a:.6f}*sin²(y) + {b:.6f}*sin(y) + {c:.6f} = 0")
            
            discriminant = b**2 - 4*a*c
            print(f"Discriminant: {discriminant:.6f}")
            
            if discriminant >= 0:
                sqrt_disc = math.sqrt(discriminant)
                sin_y1 = (-b + sqrt_disc) / (2*a)
                sin_y2 = (-b - sqrt_disc) / (2*a)
                
                print(f"sin(y) solutions: {sin_y1:.6f}, {sin_y2:.6f}")
                print(f"Known sin(y): {sin_y:.6f}")
                
                for i, sin_y_sol in enumerate([sin_y1, sin_y2]):
                    if abs(sin_y_sol) <= 1:
                        cos_y_sol = (k - d1 * sin_y_sol) / d0
                        identity_check = cos_y_sol**2 + sin_y_sol**2
                        print(f"Solution {i+1}: cos(y)={cos_y_sol:.6f}, sin(y)={sin_y_sol:.6f}")
                        print(f"  Identity check: {identity_check:.6f}")
                        print(f"  Valid: {abs(identity_check - 1) < 1e-10}")
                        
                        if abs(identity_check - 1) < 1e-10:
                            y_sol = math.atan2(sin_y_sol, cos_y_sol)
                            print(f"  y = {y_sol:.6f} ({math.degrees(y_sol):.1f}°)")
                            
                            # Final verification
                            residual1 = A[0,0]*cos_x_test + A[0,1]*sin_x_test + B[0,0]*cos_y_sol + B[0,1]*sin_y_sol - C[0]
                            residual2 = A[1,0]*cos_x_test + A[1,1]*sin_x_test + B[1,0]*cos_y_sol + B[1,1]*sin_y_sol - C[1]
                            print(f"  Final residuals: {residual1:.2e}, {residual2:.2e}")
            else:
                print("❌ Negative discriminant!")
        else:
            print("❌ d0 ≈ 0, need different approach")
    else:
        print("❌ Required vector not in range of B")

if __name__ == "__main__":
    debug_rank1_geometry()
