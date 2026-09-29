#! /usr/bin/env python3

from collections import namedtuple
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import scipy.special
from scipy.optimize import fsolve

def sigmoid(x, A, k, x0):
    return A * scipy.special.expit(k * (x - x0))

class SigmoidalSteps:

    def __init__(self, baseline, list_of_steps):

        _Point = namedtuple('Point', ['x', 'y'])
        _Step = namedtuple('Step', ['p1', 'p2', 'new_level'])

        self.baseline = baseline
        self.sigmoids = []

        last_level = baseline
        for entry in list_of_steps:

            step = _Step(_Point(*entry[0]), _Point(*entry[1]), entry[2])

            amp = step.new_level - last_level

            startvals = [2./(step.p2.x - step.p1.x), 0.5*(step.p1.x + step.p2.x)]
            k, x0 = fsolve(lambda x: [last_level + sigmoid(step.p1.x, amp, x[0], x[1]) - step.p1.y, last_level + sigmoid(step.p2.x, amp, x[0], x[1]) - step.p2.y], startvals)
            
            self.sigmoids.append(partial(sigmoid, A=amp, k=k, x0=x0))
            last_level = step.new_level

    def __call__(self, x):
        v = self.baseline
        for f in self.sigmoids:
            v += f(x)
        return v


if __name__ == '__main__':
    baseline = 0.5
    #steps = [Step(Point(2., 1.0), Point(3., 1.5), 2.0),
    #         Step(Point(22., 0.0), Point(25., -1.0), -3.0)]
    steps = [((2., 1.0), (3., 1.5), 2.0),
             ((22., 0.0), (25., -1.0), -3.0)]
    s = SigmoidalSteps(baseline, steps)
    xx = np.linspace(0., 50., 1000)
    ys = s(xx)

    plt.figure()
    plt.axhline(baseline, ls=':', color='grey')
    plt.plot(xx, ys, 'b-')
    for step in steps:
        plt.axhline(step[2], ls=':', color='lightgrey')
        plt.plot([step[0][0], step[1][0]], [step[0][1], step[1][1]], 'ro')

    plt.show()
