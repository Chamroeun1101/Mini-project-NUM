# Root Finding Toolbox

A Python implementation of classical numerical root-finding methods, solving equations of the form **f(x) = 0**.  
All algorithms are implemented manually using only standard Python and the `cmath` library.

---

## Project Description

This project provides a clean, object-oriented toolbox for numerical root finding through a single class `RootFindingProblem`.  
It supports both real and complex roots, handles edge cases with meaningful error messages, and exposes all methods through one unified interface: `solve()`.

---

## List of Implemented Methods

| # | Method | Callable Name |
|---|--------|---------------|
| 1 | Bisection Method | `"bisection"` |
| 2 | Fixed-Point Iteration | `"fixed_point"` |
| 3 | Newton's Method | `"newton"` |
| 4 | Secant Method | `"secant"` |
| 5 | False Position (Regula Falsi) | `"false_position"` |
| 6 | Steffensen's Method | `"steffensen"` |
| 7 | Horner's Method | `"horner"` |
| 8 | Muller's Method | `"muller"` |

---

## Algorithms and How They Work

### 1. Bisection Method
Given an interval `[a, b]` where `f(a)` and `f(b)` have opposite signs, the method repeatedly halves the interval and keeps the half where the sign change occurs. Converges slowly but is guaranteed to find a root.

### 2. Fixed-Point Iteration
Rewrites `f(x) = 0` as `x = g(x)`, then iterates `x_new = g(x)` until the value stabilizes. Convergence depends on `|g'(x)| < 1` near the root.

### 3. Newton's Method
Uses the tangent line at the current point to find the next approximation:
```
x_new = x - f(x) / f'(x)
```
Converges quadratically but requires the derivative `f'(x)`.

### 4. Secant Method
Similar to Newton's method but approximates the derivative using two previous points instead of an analytical formula:
```
x_new = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
```
No derivative needed, but convergence is slightly slower than Newton.

### 5. False Position (Regula Falsi)
Like bisection, it keeps an interval `[a, b]` where the sign changes. But instead of taking the midpoint, it uses a linear interpolation between `(a, f(a))` and `(b, f(b))` to pick the next point. Faster than bisection in most cases.

### 6. Steffensen's Method
An accelerated fixed-point method using Aitken's delta-squared process. Achieves quadratic convergence without needing the derivative:
```
x_new = x - f(x)^2 / (f(x + f(x)) - f(x))
```

### 7. Horner's Method
An efficient algorithm for evaluating polynomials. Instead of computing each power separately, it uses nested multiplication:
```
p(x) = a_n*x^n + ... + a_1*x + a_0
     = (...((a_n * x + a_{n-1}) * x + a_{n-2}) * x ... + a_0)
```
Reduces n multiplications and n additions, avoiding expensive power operations.

### 8. Muller's Method
Fits a parabola through three points `(x0, f(x0))`, `(x1, f(x1))`, `(x2, f(x2))` and finds the root of that parabola using the quadratic formula. Uses `cmath` to handle complex roots. Useful for finding complex or repeated roots.

---

## File Structure

```
root-finding-project/
│
├── root_finding.py   # RootFindingProblem class with all algorithms
├── examples.py       # Demonstrations of all methods
└── README.md         # Project documentation
```

---

## How to Run the Examples

Make sure you have **Python 3** installed. No external libraries are required.

**Step 1:** Clone the repository
```bash
git clone https://github.com/your-username/root-finding-project.git
cd root-finding-project
```

**Step 2:** Run the examples
```bash
python examples.py
```

You should see the root found by each method printed in the terminal along with verification.

---

## Short Code Example Using `solve()`

```python
from root_finding import RootFindingProblem
import math

# Define the function and its derivative
f  = lambda x: math.cos(x) - x
df = lambda x: -math.sin(x) - 1

# Create the problem instance
p = RootFindingProblem(f=f, df=df)

# Solve using different methods
print(p.solve("bisection",      a=0, b=1))
print(p.solve("newton",         x0=0.5))
print(p.solve("secant",         x0=0, x1=1))
print(p.solve("false_position", a=0, b=1))
print(p.solve("steffensen",     x0=0.5))

# Polynomial evaluation using Horner's method
# p(x) = x^3 - 6x^2 + 11x - 6  at x = 1
p2 = RootFindingProblem()
print(p2.solve("horner", coeffs=[1, -6, 11, -6], x=1))  # → 0

# Complex root using Muller's method
# f(x) = x^2 + 1  →  roots are ±i
p3 = RootFindingProblem(f=lambda x: x**2 + 1)
print(p3.solve("muller", x0=0, x1=1, x2=2))  # → complex number
```

