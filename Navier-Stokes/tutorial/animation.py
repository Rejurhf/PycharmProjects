import matplotlib.pyplot as plt
import numpy as np
import os

class Animation:


    def __init__(self):
        print(os.getcwd())


    def runAnimation(self):
        x = np.arange(0, 2.2, 0.2)
        y = np.arange(0, 2.2, 0.2)

        X, Y = np.meshgrid(x, y)
        u = np.cos(X) * Y
        v = np.sin(y) * Y

        tmpU = np.amax(u)
        tmpV = np.amax(v)
        print(tmpU, tmpU)

        n = -2
        M = np.hypot(u, v)

        fig, ax = plt.subplots(figsize=(7, 7))
        ax.quiver(X, Y, u, v, M, alpha=0.8, width=0.007, scale=tmpU/0.05)
        ax.scatter(X, Y, s=1)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.axis([-0.2, 2.3, -0.2, 2.3])
        ax.set_aspect('equal')

        plt.show()