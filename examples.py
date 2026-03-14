import math
from root_find import RootFindingProblem

print("=" * 60)
print("       ROOT FINDING TOOLBOX - EXAMPLES")
print("=" * 60)

# ------------------------------------------------------------------
# 1. BISECTION METHOD
#    f(x) = x^3 - x - 2 on [1, 2]
#    True root ≈ 1.5214
# ------------------------------------------------------------------
print("\n[1] Bisection Method")
print("    f(x) = x^3 - x - 2,  interval [1, 2]")

f1 = lambda x: x**3 - x - 2
p1 = RootFindingProblem(f=f1)
root = p1.solve("bisection", a=1, b=2)
print(f"    Root ≈ {root:.10f}")
print(f"    Verify f(root) = {f1(root):.2e}")

# ------------------------------------------------------------------
# 2. FIXED-POINT ITERATION
#    f(x) = x^3 - x - 2 → rearranged as  g(x) = (x + 2)^(1/3)
#    True root ≈ 1.5214
# ------------------------------------------------------------------
print("\n[2] Fixed-Point Iteration")
print("    g(x) = (x + 2)^(1/3),  x0 = 1.0")

g2 = lambda x: (x + 2) ** (1 / 3)
f2 = lambda x: x**3 - x - 2
p2 = RootFindingProblem(f=f2, g=g2)
root = p2.solve("fixed_point", x0=1.0)
print(f"    Root ≈ {root:.10f}")
print(f"    Verify f(root) = {f2(root):.2e}")

# ------------------------------------------------------------------
# 3. NEWTON'S METHOD
#    f(x) = cos(x) - x,  df(x) = -sin(x) - 1
#    True root ≈ 0.7391 (Dottie number)
# ------------------------------------------------------------------
print("\n[3] Newton's Method")
print("    f(x) = cos(x) - x,  x0 = 0.5")

f3  = lambda x:  math.cos(x) - x
df3 = lambda x: -math.sin(x) - 1
p3  = RootFindingProblem(f=f3, df=df3)
root = p3.solve("newton", x0=0.5)
print(f"    Root ≈ {root:.10f}")
print(f"    Verify f(root) = {f3(root):.2e}")

# ------------------------------------------------------------------
# 4. SECANT METHOD
#    f(x) = e^x - 3x,  two roots exist; starting near x=0 and x=1
#    Root near ≈ 0.6190
# ------------------------------------------------------------------
print("\n[4] Secant Method")
print("    f(x) = e^x - 3x,  x0=0, x1=1")

f4 = lambda x: math.exp(x) - 3 * x
p4 = RootFindingProblem(f=f4)
root = p4.solve("secant", x0=0, x1=1)
print(f"    Root ≈ {root:.10f}")
print(f"    Verify f(root) = {f4(root):.2e}")

# ------------------------------------------------------------------
# 5. FALSE POSITION (Regula Falsi)
#    f(x) = x^2 - 4  on [0, 3]
#    True root = 2.0
# ------------------------------------------------------------------
print("\n[5] False Position (Regula Falsi)")
print("    f(x) = x^2 - 4,  interval [0, 3]")

f5 = lambda x: x**2 - 4
p5 = RootFindingProblem(f=f5)
root = p5.solve("false_position", a=0, b=3)
print(f"    Root ≈ {root:.10f}")
print(f"    Verify f(root) = {f5(root):.2e}")

# ------------------------------------------------------------------
# 6. HORNER'S METHOD (polynomial evaluation)
#    p(x) = 2x^4 - 3x^3 + x^2 - 5x + 7 evaluated at x = 2
#    coeffs = [2, -3, 1, -5, 7]
#    Expected: 2(16) - 3(8) + 1(4) - 5(2) + 7 = 32-24+4-10+7 = 9
# ------------------------------------------------------------------
print("\n[6] Horner's Method (polynomial evaluation)")
print("    p(x) = 2x^4 - 3x^3 + x^2 - 5x + 7,  at x = 2")

coeffs = [2, -3, 1, -5, 7]
p6 = RootFindingProblem()
value = p6.solve("horner", coeffs=coeffs, x=2)
print(f"    p(2) = {value}")
print(f"    Expected = 9")

# Bonus: use Horner inside a solver to find a root of the polynomial
# p(x) = x^3 - 6x^2 + 11x - 6  (roots at x=1,2,3)
print("\n    Bonus: find a root of x^3 - 6x^2 + 11x - 6 using bisection + Horner")
poly_coeffs = [1, -6, 11, -6]
p_horner = RootFindingProblem()

def poly_via_horner(x):
    return p_horner.solve("horner", coeffs=poly_coeffs, x=x)

p6b = RootFindingProblem(f=poly_via_horner)
root = p6b.solve("bisection", a=0.5, b=1.5)
print(f"    Root of x^3-6x^2+11x-6 near [0.5,1.5] ≈ {root:.10f}  (expected 1.0)")

# ------------------------------------------------------------------
# 7. MULLER'S METHOD (complex root)
#    f(x) = x^2 + 1  →  roots are ±i (purely imaginary)
# ------------------------------------------------------------------
print("\n[7] Muller's Method (complex root)")
print("    f(x) = x^2 + 1,  x0=0, x1=1, x2=2")
print("    Expected roots: +i or -i")

import cmath

def f7(x):
    # Accept both real and complex input
    return x**2 + 1

p7 = RootFindingProblem(f=f7)
root = p7.solve("muller", x0=0, x1=1, x2=2)
print(f"    Root ≈ {root}")
print(f"    Verify f(root) = {f7(root):.2e}")

# ------------------------------------------------------------------
# Error handling demos
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("  ERROR HANDLING DEMOS")
print("=" * 60)

# Invalid interval
print("\n[!] Bisection with invalid interval:")
try:
    p1.solve("bisection", a=2, b=3)
except ValueError as e:
    print(f"    ValueError caught → {e}")

# Missing derivative for Newton
print("\n[!] Newton without df:")
try:
    bad = RootFindingProblem(f=lambda x: x**2 - 4)
    bad.solve("newton", x0=1.0)
except ValueError as e:
    print(f"    ValueError caught → {e}")

# Missing g for fixed-point
print("\n[!] Fixed-point without g:")
try:
    bad = RootFindingProblem(f=lambda x: x**2 - 4)
    bad.solve("fixed_point", x0=1.0)
except ValueError as e:
    print(f"    ValueError caught → {e}")

print("\n" + "=" * 60)
print("  All examples completed successfully.")
print("=" * 60)