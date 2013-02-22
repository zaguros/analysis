import numpy as np
from matplotlib import pyplot as plt
import sympy

class Formula:

    def __init__(self):
        self.formula = None
        self.values = {}
        self.uncertainties = {}

    def set_symbol_value(self, symbol, value):
        self.values[symbol] = value

    def symbol_value(self, symbol):
        return self.values[symbol]

    def set_symbol_uncertainty(self, symbol, value):
        self.uncertainties[symbol] = value

    def symbol_uncertainty(self, symbol):
        return self.uncertainties[symbol]

    def value(self):
        return self.formula.evalf(subs=self.values)

    def uncertainty(self):
        usquared = 0.
        for s in self.uncertainties:
            usquared += float((self.formula.diff(s).evalf(subs=self.values)*\
                    self.uncertainties[s])**2)
       
        return np.sqrt(usquared)


if __name__ == '__main__':
    f = Formula()
    a,b,x = sympy.symbols('a,b,x')
    f.formula = a*x + b
    
    f.values[a] = 2
    f.values[b] = 1
    f.uncertainties[a] = 0.2
    
    yvals = []
    yerrs = []
    xvals = np.arange(10)
    for xval in xvals:
        f.values[x] = xval
        yvals.append(f.value())
        yerrs.append(f.uncertainty())
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(xvals, yvals, yerr=yerrs, fmt='o')

