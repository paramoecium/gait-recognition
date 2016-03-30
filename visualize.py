import numpy as np
import matplotlib
matplotlib.use('Agg') # for linux server without display
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
    from featureEncoding import dynamicTimeWarp
    cost = dynamicTimeWarp(sensor_data['Axis1'],sensor_data['Axis1'])
    plt.figure()
    title = 'Cross Similarity'
    cmap = plt.cm.Blues
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
    output_dir = './slipper_output/Charlie/'
    sensor_data = data = np.genfromtxt('./ICS_slipperData/Charlie0105db.csv', dtype=float, delimiter=',', names=True)
    
    from scipy.interpolate import interp1d
    print 'original...'
    timestamp = sensor_data['Timestamp']
    signal = sensor_data['Axis1']
    plt.figure(figsize=(16, 6), dpi=100)
    plt.xlim(timestamp[0], timestamp[-1])
    plt.plot(timestamp, signal)
    plt.savefig(output_dir+'signal_original.png',dpi=500)
    
    print 'interpolation...'
    sample_frequency = 1000.0
    timestamp_sampled = np.arange(timestamp[0], timestamp[-1], 1/sample_frequency)
    f = interp1d(timestamp, signal, kind='slinear')
    signal_sampled = f(timestamp_sampled)
    plt.figure(figsize=(16, 6), dpi=100)
    plt.xlim(timestamp[0], timestamp[-1])
    plt.title("sample frequency={}Hz".format(sample_frequency))
    plt.plot(timestamp_sampled, signal_sampled)
    plt.savefig(output_dir+'signal_sample.png',dpi=500)

    print 'filtering...'
    from scipy.signal import butter, lfilter
    nyq = 0.5 * sample_frequency
    order = 4
    highcut = 8.0
    b, a = butter(order, highcut/nyq, btype='low')
    draw_frequency_response(b,a)
    #b, a = butter(order, [lowcut/nyq, highcut/nyq], btype='band')
    signal_filtered = lfilter(b, a, signal_sampled)
    plt.figure(figsize=(16, 6), dpi=100)
    plt.xlim(timestamp[0], timestamp[-1])
    plt.title("sample frequency={}Hz, low-pass cutoff={}Hz".format(sample_frequency,highcut))
    plt.plot(timestamp_sampled, signal_filtered)
    plt.savefig(output_dir+'signal_filtered_{}.png'.format(highcut),dpi=500)
    f = interp1d(timestamp_sampled, signal_filtered, kind='slinear')

    print 'cycle detection...'
    from preprocessing import cycle_extraction
    cycles = cycle_extraction(signal_filtered, sample_frequency)
    median_cycle = np.array(sorted(cycles, key = lambda x : x[2]))[len(cycles)/2]
    print '    median cycle:', median_cycle

    print 'warping...'
    target_cycle_length = int(median_cycle[2])
    signals_warped = []
    plt.figure(figsize=(8, 6), dpi=100)
    plt.xlim(0, target_cycle_length-1)
    for cycle in cycles:
        timestamp_start = timestamp_sampled[cycle[0]]
        timestamp_end = timestamp_sampled[cycle[1]]
        timestamp_warped = np.linspace(timestamp_start, timestamp_end, num=target_cycle_length, endpoint=False)
        assert len(timestamp_warped)==target_cycle_length, '{}!={}'.format(len(timestamp_warped),target_cycle_length)
        signals_warped.append(f(timestamp_warped))
        plt.plot(xrange(target_cycle_length), signals_warped[-1])
        print '    plot:', cycle
    plt.title("warped cycles of Alice0105db Axis1")
    plt.savefig(output_dir+'cycle_warped.png',dpi=500)
    np.save(output_dir+'signals_warped.npy', signals_warped)
    
    signals_warped = np.load('./signals_warped.npy')
    from featureEncoding import dtwDistanceMatrix
    #distance_matrix = dtwDistanceMatrix(signals_warped, metric='euclidean', down_sample = False)
    distance_matrix = dtwDistanceMatrix(signals_warped, metric='dtw', down_sample = True)
    np.save(output_dir+'distance_matrix.npy', distance_matrix)
    
    signals_warped = np.load('./signals_warped.npy')
    distance_matrix = np.load('./distance_matrix.npy')
    means = np.mean(distance_matrix, axis=1)
    plt.figure(figsize=(8, 6), dpi=100)
    plt.title("representative cycle of Alice0105db Axis1")
    plt.plot(signals_warped[np.argmin(means)])
    plt.savefig(output_dir+'cycle_best.png',dpi=500)
    mask = means < np.mean(means)-0.5*np.std(means)
    print '    Selected:', sum(mask)
    signals_selected = signals_warped[mask]
    
    plt.figure(figsize=(8, 6), dpi=100)
    #plt.xlim(0, target_cycle_length-1)
    for signal in signals_selected:
        plt.plot(signal)
    plt.title("selected cycles of Alice0105db Axis1")
    plt.savefig(output_dir+'cycle_selected.png',dpi=500)
        
    signals_warped_stat = np.array([np.mean(signals_selected,axis=0),np.std(signals_selected,axis=0)])
    plt.figure(figsize=(8, 6), dpi=100)
    plt.title("envelope of Alice0105db Axis1")
    plt.plot(signals_warped_stat[0])
    plt.plot(np.add(signals_warped_stat[0],signals_warped_stat[1]),'--')
    plt.plot(np.subtract(signals_warped_stat[0],signals_warped_stat[1]),'--')
    plt.savefig(output_dir+'cycle_envelope.png',dpi=500)

def draw_frequency_response(b,a):
    from scipy.signal import freqs
    w, h = freqs(b, a)
    plt.figure(figsize=(8, 6), dpi=100)
    plt.plot(w, 20 * np.log10(abs(h)))
    plt.xscale('log')
    plt.title('Butterworth filter frequency response')
    plt.xlabel('Frequency [radians / second]')
    plt.ylabel('Amplitude [dB]')
    plt.margins(0, 0.1)
    plt.grid(which='both', axis='both')
    plt.axvline(100, color='green') # cutoff frequency
    plt.savefig('frequency_response.png',dpi=500)

def draw_wavelet_transform():
    from scipy.interpolate import interp1d
    import scipy.fftpack
    from scipy import signal
    sensor_data = data = np.genfromtxt('./ICS_slipperData/Alice0105db.csv', dtype=float, delimiter=',', names=True)
    timestamp = sensor_data['Timestamp']
    signal = sensor_data['Axis1']

    sample_frequency = 500
    timestamp_sampled = np.arange(timestamp[0], timestamp[-1], 1.0/sample_frequency)
    f = interp1d(timestamp, signal, kind='slinear')
    signal_sample = f(timestamp_sampled)
    widths = np.arange(1, 200)
    cwtmatr = scipy.signal.cwt(scipy.fftpack.fft(f(timestamp_sampled)), scipy.signal.ricker, widths)
    print len(timestamp_sampled)
    print cwtmatr.shape
    plt.figure(figsize=(16, 6), dpi=100)
    plt.set_size_inches(8.27, 11.69)
    #plt.xticks(timestamp_sampled)
    plt.imshow(cwtmatr)
    #plt.show()
    plt.savefig('signal_cwt.png', dpi=600, bbox_inches='tight')

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
    draw_signal()
    #draw_wavelet_transform()
    #draw_example()
    '''
    from preprocessing import *
    cross_match()
    '''
