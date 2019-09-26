from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math as mth

class Flow2d:
    def __init__(self):
        print("flow 2d")


    def calcUDecrease(self, u, v):
        minXPos = len(u[0])
        maxXPos = -1
        minYPos = len(u)
        maxYPos = -1
        for i in range(len(u)):
            for x in range(len(u[0])):
                if u[i,x] == 0:
                    if x < minXPos:
                        minXPos = x
                    if x > maxXPos:
                        maxXPos = x
                    if i < minYPos:
                        minYPos = i
                    if i > maxYPos:
                        maxYPos = i

        print("X: {} {}, Y: {} {}".format(minXPos, maxXPos, minYPos, maxYPos))

        mul = 0.9
        iter = 0
        yHeight = maxYPos - minYPos + 1
        begYPos = minYPos + yHeight
        if yHeight % 2 != 0:
            yHeight += 1


        # Do work before obstacle
        # while yHeight - 2*iter > 0:
        #     count = 0
        #     for i in range(int(yHeight/2)-1, iter-1, -1):
        #         tmpMul = mul
        #         for j in range(0, count):
        #             tmpMul *= 0.8
        #         # u[minYPos+i,minXPos-iter-1] = u[minYPos+i,minXPos-iter-1] * (1-tmpMul)
        #         # u[maxYPos-i,minXPos-iter-1] = u[maxYPos-i,minXPos-iter-1] * (1-tmpMul)
        #         count += 1
        #
        #         newMul = (begYPos-(minYPos+i))/(yHeight+1)
        #         # u[minYPos+i,minXPos-iter-1] = u[minYPos+i,minXPos-iter-1] * (1-tmpMul)
        #         u[maxYPos-i,minXPos-iter-1] = u[maxYPos-i,minXPos-iter-1] * (1-tmpMul)
        #
        #         print(minXPos-iter-1, i, newMul, tmpMul, minYPos+i, maxYPos-i)
        #         print("-------------")
        #     iter += 1
        #     mul *= 0.5

        mul = 0.9
        iter = 0
        # Do work after obstacle
        while yHeight - 2*iter > 0:
            count = 0
            for i in range(int(yHeight/2)-1, iter-1, -1):
                tmpMul = mul
                for j in range(0, count):
                    tmpMul *= 0.6
                u[minYPos+i, maxXPos+(2*iter)+1] = u[minYPos+i, maxXPos+(2*iter)+1] * (1-tmpMul)
                u[maxYPos-i, maxXPos+(2*iter)+1] = u[maxYPos-i, maxXPos+(2*iter)+1] * (1-tmpMul)
                u[minYPos+i, maxXPos+(2*iter)+2] = u[minYPos+i, maxXPos+(2*iter)+2] * (1-tmpMul*0.7)
                u[maxYPos-i, maxXPos+(2*iter)+2] = u[maxYPos-i, maxXPos+(2*iter)+2] * (1-tmpMul*0.7)
                count += 1
                # print(minXPos+(2*iter)+1, i, tmpMul, minYPos+i, maxYPos-i)
                # print(u[minYPos+i, minXPos+(2*iter)+1], u[minYPos+i, minXPos+(2*iter)+2])
                # print("-------------")
            iter += 1
            mul *= 0.6

        return u, v


    def calcInObstacleFlow(self, u, v, xSize, ySize, obstacle):
        xMul = (len(u[0])-1)/(xSize)
        yMul = (len(u)-1)/(ySize)

        left = mth.ceil(obstacle[0]*xMul)
        right = mth.floor((obstacle[0]+obstacle[1])*xMul)
        bottom = mth.ceil(obstacle[2]*yMul)
        top = mth.floor((obstacle[2]+obstacle[3])*yMul)

        u[bottom:top+1, left:right+1] = 0

        return u, v


    def runSimulation(self):
        print("run")

        nx = 51
        ny = 26
        xSize = 20
        ySize = 10
        x = np.linspace(0,xSize,nx) # last point included, so exactly nx points
        y = np.linspace(0,ySize,ny) # last point included, so exactly ny points
        X,Y = np.meshgrid(x,y)       # makes 2-dimensional mesh grid

        # Draw in xSize, ySize [left, width, bottom, height]
        obstacle = [7,1.5,2.1,5.9]

        v = np.zeros((ny, nx))
        u = np.ones((ny, nx))    # for u-velocity I initialise to 1 everywhere

        u, v = self.calcInObstacleFlow(u, v, xSize, ySize, obstacle)
        u, v = self.calcUDecrease(u, v)

        p = np.zeros((ny, nx)) # np.add(np.absolute(u), np.absolute(v))

        # u[u == 1] = 0
        self.showPlot(X, Y, u, v, p, obstacle, "Title")


    def showPlot(self, X, Y, u, v, p, obstacle, titleText="no text"):
        # Plot the last figure on screen
        fig = plt.figure(figsize=(100, 50), dpi=25)
        plt.contourf(X, Y, p, alpha=0.5)  # alpha - background intensity
        plt.tick_params(axis='both', which='major', labelsize=80)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=80)
        plt.contour(X, Y, p)
        M = np.hypot(u, v)
        plt.quiver(X, Y, u, v, M, scale=1 / 0.02)  ##plotting velocity
        # plt.scatter(X, Y, color='r')
        plt.broken_barh([(obstacle[0], obstacle[1])], (obstacle[2], obstacle[3]), facecolors='grey', alpha=0.8)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title(titleText, fontsize=80)
        plt.show()
