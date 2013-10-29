import numpy as np
import matplotlib.pyplot as plt

"""
Classes used to simulate different data sets
"""

class Simulate(object):
    "A class to simulate multiple types of data"
    
    def __init__(self,seed=101):
        """
        Constructor
        """

        self.seed = seed

    def set_rseed(self,seed=None):
        """
        set the random seed
        """
        
        if seed != None:
            self.seed = seed
        
        np.random.seed(seed)
        

    def draw_u(self,M=10,markerIndices=[0,5,9]):
        """
        u - the observed imput variables in the system
        M - the total number of imput variables
        """
        
        allIndices = range(M)
        for i in markerIndices:
            if i not in allIndices:
                raise Exception("Index %s not in range(%s)"%(i,M))

        nonMarkerIndices = set(range(M)).difference(set(markerIndices))


    def linear_draw(self,n,w0=-0.3,w1=0.5,sigma=0.2):
        """
        draw samples with an approximatly linear relationship

        n     - number of samples to 
        w0    - intercept coefficient (truth)
        w1    - regression coefficient (truth)
        sigma - error variance (truth)
        """

        self.set_rseed()
        trueX = np.random.uniform(-1,1,n)
        trueT = w0 + (w1*trueX)
        return trueX, trueT + np.random.normal(0,sigma,n)

    def sine_draw(self,n,sigma=0.3):
        """
        draw samples with a sine wave relationship
        """

        self.set_rseed()
        trueX = np.linspace(0,1.0*np.pi,n) #np.linspace(0,3*np.pi,n)
        trueT = np.sin(trueX) + np.random.normal(0,sigma,n)

        return trueX,trueT


    def random_draw(self,n):
        """
        draw random uniform samples
        """

        self.set_rseed()
        trueX = np.random.uniform(0,1,n)
        trueT = np.random.uniform(0,1,n)

        return trueX,trueT


    def plot_all(self,n,figName=None):
        """
        plot all of the simulated data
        """
        
        fig = plt.figure()
        ax = fig.add_subplot(1,3,1)
        n = 50
        x,y = self.linear_draw(n)
        ax.scatter(x,y)
        ax.set_aspect(1./ax.get_data_ratio())

        ax = fig.add_subplot(1,3,2)
        x,y = self.sine_draw(n)
        ax.scatter(x,y)
        ax.set_aspect(1./ax.get_data_ratio())

        ax = fig.add_subplot(1,3,3)
        x,y = self.random_draw(n)
        ax.scatter(x,y)
        ax.set_aspect(1./ax.get_data_ratio())

        ## save or show the plots
        fig.tight_layout()
        if figName != None:
            fig.savefig(figName,dpi=200)
        else:
            plt.show()

## used to test
if __name__ == "__main__":
    sim = Simulate()
    sim.draw_u()
    #x,y = sim.linear_draw(10)
    #print 'x',x
    #print 'y',y
    #sim.plot_all(10)

