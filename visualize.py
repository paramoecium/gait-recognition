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

def draw_manifold():
    featureNum = 100
    X_train, Y_train = readFeature('./svm_train', featureNum)
    colors = np.random.rand(len(Y_train))
    colors = Y_train
    X_iso = manifold.Isomap(n_neighbors=10, n_components=3).fit_transform(X_train)
    '''
    plt.scatter(X_iso[:,0], X_iso[:,1], c=colors)
    plt.show()
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_iso[:,0], X_iso[:,1], X_iso[:,2], c=colors, marker='o')    
    plt.show()


def draw_cross_similarity():
    sensor_data = data = np.genfromtxt('./ICS_slipperData/Alice0105db.csv', dtype=float, delimiter=',', names=True)
    plt.figure()
    from featureEncoding import dynamicTimeWarp
    title = 'Cross Similarity'
    cmap = plt.cm.Blues
    cost = dynamicTimeWarp(sensor_data['Axis1'],sensor_data['Axis1'])
    plt.imshow(cost, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    #plt.xticks(np.arange(len(label_names)), label_names, rotation=45)
    #plt.yticks(np.arange(len(label_names)), label_names)
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.show()

def draw_signal():
    sensor_data = data = np.genfromtxt('./ICS_slipperData/Alice0105db.csv', dtype=float, delimiter=',', names=True)
    plt.figure()
    timestamp = sensor_data['Timestamp']
    plt.xlim(timestamp[0], timestamp[-1])
    plt.plot(timestamp, sensor_data['Axis1'])
    plt.show()

def draw_signal_interpolate():
    sensor_data = data = np.genfromtxt('./ICS_slipperData/Alice0105db.csv', dtype=float, delimiter=',', names=True)
    plt.figure()
    timestamp = sensor_data['Timestamp']
    signal = sensor_data['Axis1']

    sample_frequency = 500
    timestamp_new = np.arange(timestamp[0], timestamp[-1], 1.0/sample_frequency)
    plt.title("sample frequency={}Hz".format(sample_frequency))
    from scipy.interpolate import interp1d
    f = interp1d(signal, timestamp)
    #f = interp1d(signal, timestamp, kind='cubic')
    plt.xlim(timestamp[0], timestamp[-1])
    plt.plot(timestamp, signal, timestamp_new, f(timestamp_new), '-')
    plt.show()

def draw_example():
    plt.figure()
    seqA = [0, 0, 0, 3, 6, 13, 25, 22, 7, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    seqB = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 12, 24, 23, 8 ,3, 1, 0, 0, 0, 0, 0]
    from featureEncoding import dynamicTimeWarp
    cmap = plt.cm.Blues
    cost = dynamicTimeWarp(seqA,seqB)
    plt.imshow(cost, interpolation='nearest', cmap=cmap)
    plt.colorbar()
    plt.tight_layout()
    plt.show()


if __name__=='__main__':
    #draw_cross_similarity()
    #draw_signal()
    #draw_signal_interpolate()
    #draw_example()

    from preprocessing import *
    pattern_mining()
