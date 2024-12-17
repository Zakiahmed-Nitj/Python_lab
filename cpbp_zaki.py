# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.  
"""

#Answer c 

import math

class Fraction:
    def __init__(self, numerator, denominator):
        if denominator == 0:
            raise ValueError("Denominator cannot be zero")
        self.numerator = numerator
        self.denominator = denominator
        self.simplify()
    
    def simplify(self):
        """Simplify the fraction by dividing the numerator and denominator by their greatest common divisor (GCD)."""
        gcd = math.gcd(self.numerator, self.denominator)
        self.numerator //= gcd
        self.denominator //= gcd
        # Ensure the denominator is positive
        if self.denominator < 0:
            self.numerator = -self.numerator
            self.denominator = -self.denominator

    def __add__(self, other):
        """Add two fractions."""
        if not isinstance(other, Fraction):
            return NotImplemented
        numerator = self.numerator * other.denominator + self.denominator * other.numerator
        denominator = self.denominator * other.denominator
        result = Fraction(numerator, denominator)
        return result

    def __sub__(self, other):
        """Subtract two fractions."""
        if not isinstance(other, Fraction):
            return NotImplemented
        numerator = self.numerator * other.denominator - self.denominator * other.numerator
        denominator = self.denominator * other.denominator
        result = Fraction(numerator, denominator)
        return result

    def __mul__(self, other):
        """Multiply two fractions."""
        if not isinstance(other, Fraction):
            return NotImplemented
        numerator = self.numerator * other.numerator
        denominator = self.denominator * other.denominator
        result = Fraction(numerator, denominator)
        return result

    def __truediv__(self, other):
        """Divide two fractions."""
        if not isinstance(other, Fraction):
            return NotImplemented
        # Multiply by the reciprocal of the other fraction
        numerator = self.numerator * other.denominator
        denominator = self.denominator * other.numerator
        if denominator == 0:
            raise ZeroDivisionError("Cannot divide by zero.")
        result = Fraction(numerator, denominator)
        return result

    def __eq__(self, other):
        """Check if two fractions are equal."""
        if not isinstance(other, Fraction):
            return NotImplemented
        return self.numerator == other.numerator and self.denominator == other.denominator

    def __str__(self):
        """Return a string representation of the fraction."""
        return f"{self.numerator}/{self.denominator}"

    def __repr__(self):
        return f"Fraction({self.numerator}, {self.denominator})"

    def decimal_value(self):
        """Return the decimal (float) value of the fraction."""
        return self.numerator / self.denominator

# Example usage:
f1 = Fraction(1, 2)
f2 = Fraction(3, 4)
f3 = f1 + f2  # Addition
f4 = f1 - f2  # Subtraction
f5 = f1 * f2  # Multiplication
f6 = f1 / f2  # Division

print(f"f1: {f1}")  # 1/2
print(f"f2: {f2}")  # 3/4
print(f"f1 + f2: {f3}")  # 5/4
print(f"f1 - f2: {f4}")  # -1/4
print(f"f1 * f2: {f5}")  # 3/8
print(f"f1 / f2: {f6}")  # 2/3
print(f"f1 decimal value: {f1.decimal_value()}")  # 0.5



## ansA


def is_match(string, pattern):
    """
    Checks if the string matches the given pattern where:
    - '?' represents any single character.
    - '*' represents zero or more characters.
    """
    m, n = len(string), len(pattern)
    
    # dp[i][j] will be True if string[0..i-1] matches pattern[0..j-1]
    dp = [[False] * (n + 1) for _ in range(m + 1)]
    
    # Base case: empty string matches empty pattern
    dp[0][0] = True
    
    # Base case: pattern can match an empty string if it consists of '*' only
    for j in range(1, n + 1):
        if pattern[j - 1] == '*':
            dp[0][j] = dp[0][j - 1]
    
    # Fill the dp table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if pattern[j - 1] == '*':
                # '*' can match zero or more characters
                dp[i][j] = dp[i][j - 1] or dp[i - 1][j]
            elif pattern[j - 1] == '?' or pattern[j - 1] == string[i - 1]:
                # '?' matches any single character or exact match
                dp[i][j] = dp[i - 1][j - 1]
    
    # The result is whether the entire string matches the entire pattern
    return dp[m][n]

# Example usage:
print(is_match("abcde", "a*de"))  # True
print(is_match("abcde", "a?c*e")) # True
print(is_match("abcde", "a?d*"))  # False
print(is_match("abcd", "a*d?"))   # False






#answer b

   

import sympy 
   print(sympy.__version__) 

import sympy as sp

def solve_equation(equation: str):
    # Define the symbol 'x' (the unknown variable)
    x = sp.symbols('x')
    
    # Replace '=' with '==' to make it compatible with sympy
    equation = equation.replace('=', '==')
    
    # Convert the equation string to a sympy expression using Eq for equality
    try:
        # Using Eq to create an equation: LHS == RHS
        lhs, rhs = equation.split('==')
        eq = sp.Eq(sp.sympify(lhs), sp.sympify(rhs))
        
        # Solve the equation for 'x'
        solution = sp.solve(eq, x)
        
        # Return the solution (usually a list of solutions, so we return the first one)
        if solution:
            return solution[0]  # return the first solution if it's a single solution
        else:
            return "No solution found"
    
    except Exception as e:
        return f"Error: {str(e)}"

# Example usage:
equation_1 = "2*x + 3 = 7"
equation_2 = "x**2 - 5*x + 6 = 0"
equation_3 = "3*x + 4 = 10"

print(f"Solution for '{equation_1}': {solve_equation(equation_1)}")  # Should return 2
print(f"Solution for '{equation_2}': {solve_equation(equation_2)}")  # Should return [2, 3]
print(f"Solution for '{equation_3}': {solve_equation(equation_3)}")  # Should return 2















     
        
        
   
    
   
         
        
        
       
    
    
    
    


            
    
        
        
    
    
    
    
        
           




    
    
    

    



     
   
       


