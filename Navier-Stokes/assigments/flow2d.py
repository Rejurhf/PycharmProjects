from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math as mth

class Flow2d:
    def __init__(self):
        print("flow 2d")

    def changeScope(self, min, max, newMin, newMax, value):
        # reassign value
        out = value
        # calc percentage value of (max - min)
        out = (out - min) / (max - min)
        # get value from new scope using percentage value
        out = out * (newMax - newMin)
        # add value to min
        out = newMin + out

        return out


    def calcObstacleBorders(self, u, xSize, ySize, obstacle, rounded=True):
        xMul = (len(u[0]) - 1) / (xSize)
        yMul = (len(u) - 1) / (ySize)

        if rounded:
            left = mth.ceil(obstacle[0] * xMul)
            right = mth.floor((obstacle[0] + obstacle[1]) * xMul)
            bottom = mth.ceil(obstacle[2] * yMul)
            top = mth.floor((obstacle[2] + obstacle[3]) * yMul)
        else:
            left = obstacle[0] * xMul
            right = (obstacle[0] + obstacle[1]) * xMul
            bottom = obstacle[2] * yMul
            top = (obstacle[2] + obstacle[3]) * yMul

        return left, right, bottom, top


    def calcUDecrease(self, u, v, xSize, ySize, obstacle):
        left, right, bottom, top = self.calcObstacleBorders(u, xSize, ySize, obstacle, rounded=False)

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

        # Calc constants
        mul = 0.9
        iter = 0
        yHeight = maxYPos - minYPos + 1
        r = (top - bottom) / 2
        begYPos = bottom + r
        begXPos = left - r/2
        if yHeight % 2 != 0:
            yHeight += 1


        # Do work before obstacle --------------------------------------------------------

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
        toVisitList = [(int(begXPos), int(begYPos))]
        visitedList = []
        print("Start pos: {}, {} {} {}".format(toVisitList[0], begYPos, minXPos, r))

        while len(toVisitList) > 0:
            currElem = toVisitList.pop()
            visitedList.append(currElem)

            distanceToMidd = mth.sqrt((minXPos - currElem[0])**2 + (begYPos - currElem[1])**2)
            distanceToObs = np.absolute(minXPos - currElem[0])
            mul = (distanceToMidd/(1.5*r))

            # Those ifs makes ellipse range around obstacle
            if distanceToMidd + (16*distanceToObs/r) - r/10 <= r:
                u[currElem[1], currElem[0]] = u[currElem[1], currElem[0]] * self.changeScope(0, 1, 0, 0.6, mul)
            elif distanceToMidd + (8*distanceToObs/r) - r/10 <= r:
                u[currElem[1], currElem[0]] = u[currElem[1], currElem[0]] * self.changeScope(0, 1, 0.125, 0.7, mul)
            elif distanceToMidd + (distanceToObs/r) - r/10 <= r:
                u[currElem[1], currElem[0]] = u[currElem[1], currElem[0]] * self.changeScope(0, 1, 0.25, 0.8, mul)
            else:
                u[currElem[1], currElem[0]] = u[currElem[1], currElem[0]] * self.changeScope(0, 1, 0.5, 0.9, mul)

            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    tmpElem = (currElem[0]+x, currElem[1]+y)
                    # Calculate distance from middle of obstacle
                    tmpDistance = mth.sqrt((begXPos - tmpElem[0])**2 + (begYPos - tmpElem[1])**2)
                    # If new elem is in radios and its value is 1 and it is not in (visitedList and toVisitList) then
                    # add to toVisitList
                    if tmpDistance <= r+0.7 and u[tmpElem[1], tmpElem[0]] == 1 and \
                            tmpElem not in visitedList and tmpElem not in toVisitList:
                        toVisitList.append(tmpElem)
                        # print("{} {} {}".format(currElem, tmpElem, tmpDistance))
            # print(currElem, np.absolute(begYPos - currElem[1]), begYPos, r)


        # Do work after obstacle ---------------------------------------------------
        mul = 0.9
        iter = 0
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

        for i in range(5, 19):
            print("({}, {}) {}".format(i, 13, u[13, i]))

        return u, v


    def calcInObstacleFlow(self, u, v, xSize, ySize, obstacle):
        left, right, bottom, top = self.calcObstacleBorders(u, xSize, ySize, obstacle, rounded=True)

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
        u, v = self.calcUDecrease(u, v, xSize, ySize, obstacle)

        p = np.zeros((ny, nx)) # np.add(np.absolute(u), np.absolute(v))
        p = u
        # u[u == 1] = 0
        u[0,0] = 2.3
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
