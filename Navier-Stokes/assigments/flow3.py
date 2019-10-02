import matplotlib.pyplot as plt
import numpy as np
import math as mth


class Flow3:
    def __init__(self):
        print("flow 3")


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


    def calcInObstacleFlow(self, u, v, xSize, ySize, obstacle):
        left, right, bottom, top = self.calcObstacleBorders(u, xSize, ySize, obstacle, rounded=True)

        u[bottom:top+1, left:right+1] = 0

        return u, v


    def isPointInObstacle(self, point, u, xSize, ySize, obstacle):
        left, right, bottom, top = self.calcObstacleBorders(u, xSize, ySize, obstacle)

        if left <= point[0] <= right and bottom <= point[1] <= top:
            return True
        return False


    def getFlowPath(self, u, v, xSize, ySize, obstacle):
        left, right, bottom, top = self.calcObstacleBorders(u, xSize, ySize, obstacle, rounded=True)

        # Number of x and y points in arrays
        nX = len(u[0])
        nY = len(u)

        print("start")
        # start flow from first points
        for y in range(0, nY):
            # visited list if point in visited list then closed flow
            visitedPoints = []

            # x and y start position of point
            posX = 0
            posY = y

            # Flow to the end
            while posX < nX:
                visitedPoints.append((posX, posY))
                if y != posY and (posY < y and not self.isPointInObstacle((posX, posY+1), u, xSize, ySize, obstacle) or \
                        posY > y and not self.isPointInObstacle((posX, posY-1), u, xSize, ySize, obstacle)):
                    if posY < y:
                        posY += 1
                    else:
                        posY -= 1
                    print(posX, posY)
                elif not self.isPointInObstacle((posX + 1, posY), u, xSize, ySize, obstacle):
                    posX += 1
                elif not self.isPointInObstacle((posX, posY + 1), u, xSize, ySize, obstacle):
                    posY += 1
                elif not self.isPointInObstacle((posX, posY - 1), u, xSize, ySize, obstacle):
                    posY -= 1
                else:
                    posX = nX

            if visitedPoints:
                for point in visitedPoints:
                    if u[point[1], point[0]] == 0:
                        u[point[1], point[0]] = 1
                    else:
                        u[point[1], point[0]] += 0.1

        print("end")

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
        u = np.zeros((ny, nx))    # for u-velocity I initialise to 1 everywhere

        u, v = self.getFlowPath(u, v, xSize, ySize, obstacle)
        # u, v = self.calcInObstacleFlow(u, v, xSize, ySize, obstacle)
        # u, v = self.calcUDecrease(u, v, xSize, ySize, obstacle)

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