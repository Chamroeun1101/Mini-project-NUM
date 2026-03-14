import cmath

class RootFindingProblem:
    def __init__(self, f=None, df=None, g=None):
        self.f = f
        self.df = df
        self.g = g
    def solve(self, method, **kwargs):
        method = method.lower().replace("-", "_").replace(" ", "_")
        dispatch = {
            "bisection": self._bisection,
            "fixed_point": self._fixed_point,
            "newton": self._newton,
            "secant": self._secant,
            "false_position": self._false_position,
            "horner": self._horner_solve,
            "muller": self._muller,
        }

        if method not in dispatch:
            raise ValueError(
                f"Unknown method '{method}'. "
                f"Choose from: {list(dispatch.keys())}"
            )
        return dispatch[method](**kwargs)
    
    def _bisection(self, a, b, tol=1e-10, max_iter=1000):
        if self.f is None:
            raise ValueError(f"f(x) must be provided.")
        
        fa = self.f(a)
        fb = self.f(b)

        if fa * fb > 0:
            raise ValueError(
                f"Invalid interval: f(a) and f(b) must have opposite signs."
                f"Got f({a}) = {fa}, f({b}) = {fb}"
            )
        
        for i in range(max_iter):
            c = (a + b) / 2.0
            fc = self.f(c)

            if abs(fc) < tol or (b - a) / 2.0 < tol:
                return c
            if fa * fb < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc
        raise RuntimeError(
            f"Bisection did not converge after {max_iter} iterations."
        )
    
    def _fixed_point(self, x0, tol=1e-10, max_iter=1000):
        if self.g is None:
            raise ValueError("Fixed-point function g(x) must be provided.")
        
        x = x0
        for i in range(max_iter):
            x_new = self.g(x)
            if abs(x_new - x) < tol:
                return x_new
            x = x_new
        raise RuntimeError(f"Fixed-point iteration did not converge after {max_iter} iterations.")

    def _newton(self, x0, tol=1e-10, max_iter=1000):
        if self.f is None:
            raise ValueError("f(x) must be provided.")
        
        if self.df is None:
            raise ValueError("Derivative df(x) must be provided for Newton's method.")
        
        x = x0
        for i in range(max_iter):
            fx = self.f(x)
            dfx = self.df(x)

            if dfx == 0: raise ZeroDivisionError(f"Derivative is zero at x={x}. Newton's method failed.")

            x_new = x - fx / dfx

            if abs(x_new - x) < tol:
                return x_new
            
            x = x_new
        raise RuntimeError(f"Newton's method did not converge after {max_iter} iterations.")
    
    def _secant(self, x0, x1, tol=1e-10, max_iter=1000):
        if self.f is None:
            raise ValueError("f(x) must be provided.")
        
        for i in range(max_iter):
            f0 = self.f(x0)
            f1 = self.f(x1)
            denom = f1 - f0

            if denom == 0:
                raise ZeroDivisionError(f"f(x1) - f(x0) = 0 at iteration {i}. Secant method failed.")
            
            x2 = x1 - f1 * (x1 - x0) / denom

            if abs(x2 - x1) < tol:
                return x2
            
            x0, x1 = x1, x2

        raise RuntimeError(f"Secant method did not converge after {max_iter} iterations.")

    def _false_position(self, a, b, tol=1e-10, max_iter=1000):
        if self.f is None:
            raise ValueError("f(x) must be provided.")
        
        fa = self.f(a)
        fb = self.f(b)

        if fa * fb > 0:
            raise ValueError(
                f"Invalid interval: f(a) and f(b) must have opposite signs."
                f"Got f({a}) = {fa}, f({b}) = {fb}."
            )
        
        for i in range(max_iter):
            #Regula falsi formula
            denom = fb - fa
            if denom == 0:
                raise ZeroDivisionError("Division by zero in false position method.")
            
            c = b - fb * (b - a) / denom
            fc = self.f(c)

            if abs(fc) < tol or abs(b - a) < tol:
                return c
            
            if fa * fc < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc
        raise RuntimeError(f"False position method did not converge after {max_iter} iterations.")

    def _steffensen(self, x0, tol=1e-10, max_iter=1000):
        """
        Steffensen's method: an accelerated fixed-point iteration using
        Aitken's delta-squared process. Achieves quadratic convergence
        without requiring a derivative.

        Requires f(x) to be provided.
        Formula:
            x_{n+1} = x_n - f(x_n)^2 / (f(x_n + f(x_n)) - f(x_n))
        """
        if self.f is None:
            raise ValueError("f(x) must be provided.")

        x = x0
        for i in range(max_iter):
            fx = self.f(x)
            fx2 = self.f(x + fx)

            denom = fx2 - fx
            if denom == 0:
                raise ZeroDivisionError(
                    f"Division by zero in Steffensen's method at iteration {i}."
                )

            x_new = x - fx**2 / denom

            if abs(x_new - x) < tol:
                return x_new

            x = x_new

        raise RuntimeError(
            f"Steffensen's method did not converge after {max_iter} iterations."
        )

    def _horner(self, coeffs, x):
        if not coeffs:
            raise ValueError("Coefficient list must not be empty.")
        
        result = coeffs[0]
        for coeff in coeffs[1:]:
            result = result * x + coeff

        return result
    
    def _horner_solve(self, coeffs, x):
        return self._horner(coeffs, x)
    
    def _muller(self, x0, x1, x2, tol=1e-10, max_iter=1000):
        if self.f is None:
            raise ValueError("f(x) must be provided.")

        x0 = complex(x0)
        x1 = complex(x1)
        x2 = complex(x2)

        for i in range(max_iter):
            f0 = self.f(x0)
            f1 = self.f(x1)
            f2 = self.f(x2)

            h1 = x1 - x0
            h2 = x2 - x1

            if h1 == 0 or h2 == 0:
                raise ZeroDivisionError(
                    "Two starting points are identical in Muller's method."
                )

            delta1 = (f1 - f0) / h1
            delta2 = (f2 - f1) / h2

            denom_a = h1 + h2
            if denom_a == 0:
                raise ZeroDivisionError(
                    "Division by zero when computing coefficient a in Muller's method."
                )

            a = (delta2 - delta1) / denom_a
            b = a * h2 + delta2
            c = f2

            # Discriminant
            discriminant = b * b - 4 * a * c

            # Choose the sign that gives the larger denominator magnitude
            sqrt_disc = cmath.sqrt(discriminant)
            denom_plus  = b + sqrt_disc
            denom_minus = b - sqrt_disc

            if abs(denom_plus) >= abs(denom_minus):
                denominator = denom_plus
            else:
                denominator = denom_minus

            if denominator == 0:
                raise ZeroDivisionError(
                    "Denominator is zero in Muller's method. "
                    "Try different starting points."
                )

            dx = -2 * c / denominator
            x3 = x2 + dx

            if abs(dx) < tol:
                return x3

            x0, x1, x2 = x1, x2, x3

        raise RuntimeError(
            f"Muller's method did not converge after {max_iter} iterations."
        )
