import numpy as np
from math import erfc, sqrt
import matplotlib.pyplot as plt

STEP = 1.0 / 20
PERIOD = int(1 / STEP)

LAMBDA = 0 # THRESHOLD

SIGMA_NOISE = 0.1

T = 1
A = 1
E = 1

def draw(t, signal, xLabel, yLabel, title):
    plt.plot(t, signal)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title)
    plt.show()
    plt.savefig(title + '.jpg')

def sample(bitstream, signal, n):
    samples = np.array([signal[PERIOD - 1 + i * PERIOD] for i in range(n)])
    result = (samples > LAMBDA)
    print("Reconstructed Bitstram:", result)
    print("Total number of bits:", n)
    print("Received Wrong:", np.sum(bitstream != result))
    print("BER:", np.sum(bitstream != result) / n)
    return result

def generateSignal(bitstream):
    # Pulse Amplitude Modulation (PAM)
    P = int(T / STEP)
    pulse = np.ones(P) * A
    
    pulseLength = len(pulse)
    inputLength = len(bitstream)
    
    # Generate Signal
    signal = np.zeros(inputLength * pulseLength)
    
    # Apply the pulse shape along the symbol's interval
    for i in range(inputLength):
        signal[(i * pulseLength):((i + 1) * pulseLength)] = pulse if bitstream[i] == 1 else pulse * -1
        
    return signal            

def getFilters():
    # CASE 1
    # h(t) is a Matched Filter with unit energy
    filter1 = np.ones(int(1 / STEP))
    filterLength = len(filter1)
    W = np.zeros(int(1 / STEP) - filterLength)
    filter1 = np.concatenate((filter1, W))
    
    # CASE 2
    # h(t) is non-existent [h(t) = delta(t)]
    filter2 = np.ones(1)
    filterLength = len(filter2)
    W = np.zeros(int(1 / STEP) - filterLength)
    filter2 = np.concatenate((filter2, W))
    
    # CASE 3
    # h(t) has the following impulse response
    filter3 = np.sqrt(3) * np.arange(0, 1, STEP)
    filterLength = len(filter3)
    W = np.zeros(int(1 / STEP) - filterLength)
    filter3 = np.concatenate((filter3, W))
    
    return [filter1, filter2, filter3]

def applyReceiveFilter(signal, filter, num):
    # Apply Convolution (y(t) = r(t) * h(t)) where r(t) = g(t) + w(t)
    result = np.convolve(signal, filter)
    
    # For the first & third filters, we'll need to multiply by the step
    # This won't be needed for the non-existent filter (Delta)
    if (num != 2):
        result = result * STEP
        
    # Return the receive filter result
    return result

# variance is the noise's variance and f is the filter's E_N. (1 or 2 or 3)
def BER(var, g, filter, num, bitstream, n):
    # Generate AWGN Noise
    signalLength = len(g)
    w = np.random.normal(0, var, signalLength)

    # Add the randomly generated noise to the signal
    r = g + w

    # Apply the receive filter
    y = applyReceiveFilter(r, filter, num)

    # Sample the output signal
    samplingOutput = sample(bitstream, y, n)

    # Compute the sampling error in the generated bitstream & return it
    return np.sum(bitstream != samplingOutput) / n

def Q(x):
    """
    This function returns the value of the Q function for a given input x, which is defined as half the complementary error function (erfc) of x divided by the square root of 2.
    """
    return 0.5 * erfc(x / sqrt(2))

def filtersBER(filters):
    # Set the bitstream length
    n = 1000000
    
    # Input Bitstream
    bitstream = np.random.randint(0, 2, n)
    
    # Generate Signal
    signal = generateSignal(bitstream)
    
    # E/N0
    E_N = np.arange(-10, 20, 1)
    N = 1 / (10 ** (E_N / 10))
    
    # Compute the variance
    variance = np.sqrt(N / 2)
    
    # BER for each filter
    numberOfFilters = len(filters)
    
    BERs_th = []
    BERs_th.append([Q(1 / var) for var in variance])
    BERs_th.append([Q(1 / var) for var in variance])
    BERs_th.append([Q((np.sqrt(3)/2) * (1 / var)) for var in variance])
    
    # Plot Theoretical BER
    plt.semilogy(E_N, BERs_th[0], 'c')
    plt.semilogy(E_N, BERs_th[1], 'm--')
    plt.semilogy(E_N, BERs_th[2], 'y')
    
    BERs = []
    for i in range(numberOfFilters):
        BERs.append([BER(var, signal, filters[i], i + 1, bitstream, n) for var in variance])
        
    # Plot the Simulated BER
    plt.semilogy(E_N, BERs[0], 'r')
    plt.semilogy(E_N, BERs[1], 'g.')
    plt.semilogy(E_N, BERs[2], 'b')
    
    # Labels/Legend/YLim
    plt.xlabel('E/N (db)')
    plt.ylabel('BER (log-scale)')
    plt.title(' BER VS. E/N')
    plt.legend(['Theory 1', 'Theory 2', 'Theory 3', 'Matched Filter (1)', 'No Filter (2)', 'Linear Filter (3)'])
    plt.ylim([10 / n, 1])
    plt.savefig('./BitErrorRate.png')
    plt.show()
    
def run():
    # Define n
    n = 10
    
    # Define the range t
    t = np.arange(0, n, STEP)
    
    # Generate a 1D array of random binary values
    bitstream = np.random.randint(0, 2, n)
    
    # Generate the signal from our bitstream
    signal = generateSignal(bitstream)
    
    # Plot the Signal
    draw(t, signal, "Time in seconds", "g(t)", "Signal")
    
    # Generate random Noise (AWGN)
    # Additive White Gaussian Noise
    signalLength = len(signal)
    noise = np.random.normal(0, SIGMA_NOISE, signalLength)
    
    # Plot the Noise
    draw(t, noise, "Time in seconds", "w(t)", "Noise")
    
    # Add the random generated noise to the signal
    noisySignal = signal + noise
    
    # Plot the Noisy Signal
    draw(t, noisySignal, "Time in seconds", "r(t)", "Noisy Signal")
    
    # Generate the receive filters for all 3 cases
    filters = getFilters()
    
    # Receive Filters
    result1 = applyReceiveFilter(noisySignal, filters[0], 1)
    result2 = applyReceiveFilter(noisySignal, filters[1], 2)
    result3 = applyReceiveFilter(noisySignal, filters[2], 3)
    
    # PLOT: CASE 1
    t = np.arange(0, len(result1) * STEP, STEP)
    draw(t, result1, "Time in seconds", "y(t)", "Filter 1 Output")
    sample(bitstream, result1, n)
    
    # PLOT: CASE 2
    t = np.arange(0, len(result2) * STEP, STEP)
    draw(t, result2, "Time in seconds", "y(t)", "Filter 2 Output")
    sample(bitstream, result2, n)
    
    # PLOT: CASE 3
    t = np.arange(0, len(result3) * STEP, STEP)
    draw(t, result3, "Time in seconds", "y(t)", "Filter 3 Output")
    sample(bitstream, result3, n)
    
    # Compute the BER for each filter h(t)
    filtersBER(filters)
    
if __name__ == '__main__':
    run()