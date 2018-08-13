import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import sys
import gzip

blockSizeX=1000
blockSizeY=30
LW=2

class Plotter:

    def __init__(self):

        ############################# PARSE ARGUMENTS ###############################

        import argparse, sys

        parser = argparse.ArgumentParser(description='Polygenic selection')

        parser._optionals.title = "help commands"

        parser.add_argument("-f", "--fileName", help="Posterior file",
                            action="store", type=str, required=True)
        parser.add_argument("-i", "--intervalsInfo", help="Intervals info file",
                            action="store", type=str, required=True)
        parser.add_argument("-p", "--posRange", help="Position start/end (in Mb)",
                            action="store", type=float, nargs=2, required=True)
        parser.add_argument("-t", "--timeRange", help="Time start/end (discrete bin indexed from 1 to N)",
                            action="store", type=int, nargs=2, required=True)
        parser.add_argument("-n", "--norm", help="File conaining average coalescent times (one time interval per row)",
                            action="store", type=str, required=False)
        parser.add_argument("-o", "--outPutFileRoot", help="output file root",
                            action="store", type=str, required=False)

        ############################# PROCESS ARGUMENTS ###############################

        np.set_printoptions(linewidth=300)
        args = parser.parse_args()

        self.fileName = args.fileName
        print("fileName: " + self.fileName)

        self.intervalsInfo = args.intervalsInfo
        print("intervalsInfo: " + self.intervalsInfo)

        self.norm = args.norm
        print("norm: " + self.norm)

        posRange = args.posRange
        self.fromPos = posRange[0]
        self.toPos = posRange[1]
        print("posRange: " + str(self.fromPos) + " " + str(self.toPos))

        timeRange = args.timeRange
        self.fromTime = timeRange[0]
        self.toTime = timeRange[1]
        print("timeRange: " + str(self.fromTime) + " " + str(self.toTime))

        self.outPutFileRoot = args.outPutFileRoot
        if self.outPutFileRoot == None:
            self.outPutFileRoot = self.fileName + "." + str(fromPos) + "." + str(toPos) + "." + str(FROM) + "." + str(TO)
        print("outFile: " + self.outPutFileRoot)

    def roundPosX(self, pos):
        return np.round(pos / blockSizeX)

    def roundPosY(self, pos):
        return np.round(pos / blockSizeY)

    def stretchPosterior(self, post, ranges):
        lastDiscreteGen = ranges[-1,2]
        stretchedPosterior = np.zeros(int(lastDiscreteGen))
        for i in range(len(post)):
            val = post[i]
            fromTime = int(ranges[i,0])
            toTime = int(ranges[i,2])
            stretchedPosterior[fromTime:toTime] = val
        return stretchedPosterior

    def plot(self):

        fromPos = self.fromPos*1e6
        toPos = self.toPos*1e6
        FROM = self.fromTime-1
        TO = self.toTime-1

        f=gzip.open(self.fileName,"r")
        lines=f.readlines()
        X=[]
        C=[]
        for x in lines:
            split = x.split()
            # print(split)
            X.append(int(split[3]))
            C.append(split[4:])
        minPos=X[0]
        maxPos=X[-1]
        X = np.array(X)
        C = np.array(C)
        f.close()
        data = C.T.astype(float)

        norm = np.loadtxt(self.norm)
        norm /= sum(norm)
        norm = norm[FROM:TO]
        N = np.loadtxt(self.intervalsInfo)
        ranges = self.roundPosY(N[FROM:TO,0:3])
        yLabels = N[:,2]
        data = data[FROM:TO,:]
        post = (data.T/norm)

        ###### PLOT
        expectedTimes = ranges[:,1]
        # discretize posterior
        firstPos = X[0]
        lastPos = X[-1]
        lastDiscreteGen = ranges[-1,2]
        states = lastDiscreteGen
        discreteFirstPos = self.roundPosX(firstPos)
        discreteLastPos = self.roundPosX(lastPos)
        discreteLen = int(discreteLastPos + 1)
        discretePost = np.zeros((discreteLen, int(states)))
        last=firstPos - 0
        discreteRanges = []
        for i in range(1, len(X)):
            discreteRanges.append([self.roundPosX(last), self.roundPosX((X[i-1] + X[i]) / 2)])
            last = (X[i-1] + X[i]) / 2
        discreteRanges.append([self.roundPosX(last), self.roundPosX(lastPos)+1])
        curr = 0
        for j in range(len(discreteRanges)):
            r = discreteRanges[j]
            stretchedPosterior = self.stretchPosterior(post[j,:], ranges)
            discretePost[int(max(curr, r[0])) : int(r[1]), :] = stretchedPosterior
            curr += int(r[1]) - max(curr, r[0])
        # plot
        plt.figure(figsize=(15,3.5))
        plt.ticklabel_format(useOffset=False, style='plain')
        plt.gcf().subplots_adjust(bottom=0.2)
        pos1 = self.roundPosX(max(fromPos, minPos))
        pos2 = self.roundPosX(min(toPos, maxPos))
        discretePost[0:int(pos1),:]=0
        discretePost[int(pos2):,:]=0
        heatmap = plt.pcolormesh(discretePost.transpose(),cmap="viridis",linewidth=0,rasterized=True)
        heatmap.set_edgecolor('face')
        plt.ylim((0,ranges[-1,0]))
        plt.xlim(pos1, pos2)
        ax = plt.gca()
        from matplotlib.ticker import MaxNLocator
        locatorY=MaxNLocator(nbins=5)
        ax.yaxis.set_major_locator(locatorY)
        locatorX=MaxNLocator(prune='both', nbins=8)
        ax.xaxis.set_major_locator(locatorX)
        ticks = np.round(ax.get_xticks()/1000,2)
        ax.set_xticklabels(ticks)

        fSize=18
        plt.ylabel(r'Years (thousands)', fontsize=fSize)
        plt.xlabel('Position (Mb)', fontsize=fSize)
        plt.tick_params(axis='both', labelsize=fSize)
        cbar = plt.colorbar()
        cbar.ax.set_title('Enrichment', fontsize=fSize)

        cbar.ax.tick_params(labelsize=fSize)
        print("Saving to " + self.outPutFileRoot + '.heat.pdf')
        plt.savefig(self.outPutFileRoot + '.heat.pdf', format='pdf')
        plt.tight_layout()
        plt.close()


if __name__ == '__main__':
    plotter = Plotter()
    plotter.plot()
