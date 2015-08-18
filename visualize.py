import numpy as np
import matplotlib.pyplot as plt
import csv
from featureEncoding import *
from sklearn import manifold
from mpl_toolkits.mplot3d import Axes3D

def readArray(outputFilename):
    result = []
    f = open(outputFilename, 'r')  
    for row in csv.reader(f):  
        result.append([float(e) for e in row])  
    f.close()
    return result

if __name__=='__main__':
    featureNum = 100
    X_train, Y_train = readFeature('./svm_train', featureNum)
    colors = Y_train
    X_iso = manifold.Isomap(n_neighbors=10, n_components=3).fit_transform(X_train)
    print("Done.")
    '''
    plt.scatter(X_iso[:,0], X_iso[:,1], c=colors, alpha=0.5)
    plt.show()
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_iso[:,0], X_iso[:,1], X_iso[:,2], c=colors, marker='o',alpha=0.2)    
    plt.show()
