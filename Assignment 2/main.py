import numpy as np
from math import erfc, sqrt
import matplotlib.pyplot as plt

STEP = 1.0 / 20
PERIOD = int(1 / STEP)

LAMBDA = 0  # THRESHOLD

SIGMA_NOISE = 0.1

T = 1
A = 1
E = 1


def draw(t, signal, xLabel, yLabel, title):
    """
    This function plots a given signal against the corresponding time axis, with specified labels for the x-axis and y-axis, and a title. The plot is displayed and saved as a JPG file using the provided title.
    """
    plt.plot(t, signal)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title)
    plt.savefig(title + '.jpg')
    plt.show()


def sample(bitstream, signal, n):
    """
    This function samples a given signal at a specified interval and compares the samples with a threshold value. It reconstructs a binary bitstream based on the sample comparison, calculates the Bit Error Rate (BER), and returns the reconstructed bitstream.
    """
    samples = np.array([signal[PERIOD - 1 + i * PERIOD] for i in range(n)])
    result = (samples > LAMBDA)
    print("Reconstructed Bitstram:", result)
    print("Total number of bits:", n)
    print("Received Wrong:", np.sum(bitstream != result))
    print("BER:", np.sum(bitstream != result) / n)
    return result


def generateSignal(bitstream):
    """
    This function generates a signal waveform based on a given binary bitstream using Pulse Amplitude Modulation (PAM) with a specified pulse shape. The function applies the pulse shape along the symbol's interval and returns the resulting signal.
    """
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
    """
    This function retrieves three different filters for the given communication system.
    Case 1: Matched Filter with unit energy.
    Case 2: Non-existent filter (delta function).
    Case 3: Filter with a specific impulse response.
    The function returns a list containing the three filters.
    """
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
    """
    This function applies a receive filter to a given signal using convolution. The resulting filtered signal is returned. If the filter is the first or third filter, the filtered signal is multiplied by the step size.
    """
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
    """
    This function calculates the Bit Error Rate (BER) for a given communication system. It generates additive white Gaussian noise (AWGN) and adds it to the signal. The signal is then passed through a receive filter, sampled, and compared with the original bitstream to compute the BER. The BER value is returned.
    """
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
    """
    This function calculates and plots the Bit Error Rate (BER) for different filters in a communication system. It generates a random bitstream, generates a signal based on the bitstream, and computes the theoretical and simulated BERs for each filter under different E/N0 values. The BER values are plotted on a log-scale graph along with the theoretical BER values.
    """
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
        BERs.append([BER(var, signal, filters[i], i + 1, bitstream, n)
                    for var in variance])

    # Plot the Simulated BER
    plt.semilogy(E_N, BERs[0], 'r')
    plt.semilogy(E_N, BERs[1], 'g.')
    plt.semilogy(E_N, BERs[2], 'b')

    # Labels/Legend/YLim
    plt.xlabel('E/N (db)')
    plt.ylabel('BER (log-scale)')
    plt.title(' BER VS. E/N')
    plt.legend(['Theory 1', 'Theory 2', 'Theory 3',
               'Matched Filter (1)', 'No Filter (2)', 'Linear Filter (3)'])
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

""" 
Q(5)

The Bit Error Rate (BER) decreases as the energy per bit (E) relative to the noise power spectral density (No) increases, which is evident in the plot. This can be explained in various ways:

Using the theoretical expression for BER and noting that the Q function is a decreasing function, it is apparent that Q(a * sqrt(E/No)) for all the cases mentioned in the problem. Therefore, as sqrt is an increasing function, it can be concluded that BER is a decreasing function of E/No.
"""

""" 
Q(6)

The case that uses a matched filter has the lowest BER because it employs a filter that is specifically designed to minimize the probability of error by matching the filter to the pulse. By doing so, the matched filter maximizes the peak pulse SNR at the sampling instant, which helps reduce the probability of error.

It is important to note that in the theoretical case, using a filter or not yields the same expression due to our assumptions on variance and PSD.
"""
