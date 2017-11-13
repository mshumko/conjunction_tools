import numpy as np
import random
import matplotlib.pylab as plt
from statistics import mode
from datetime import datetime, timedelta

# This function will find the locations when time stamps are contonous.
def time_continuity(t):
    """
    NAME:    time_continuity(t)   
    USE:     Calculates where the supplied time array is continous.
    INPUT:   A numpy time array
    RETURNS: A dictionary with 'startTime', 'endTime', and 'duration's for the 
             time array given. The duration array is given in data points,
             not time!
    AUTHOR: Mykhaylo Shumko
    MOD:    2017-06-06
    """
    # Find the dt between subsequent points.
    dt = t[1:] - t[:-1]
    dt = np.array(list(map(lambda x: x.total_seconds(), dt )))
    
    if len(t) == 1: # special case where there is only one value.
        tDict = {'startTime':np.array([0], dtype = int), 
        'endTime':np.array([1], dtype = int), 
        'duration':np.array([1], dtype = int)}
        return tDict

    #print(t, dt)

    # Subtract the min (should be the cadence) to flag the jumps in dt.
    dt -= min(dt)

    # Find where the jumps occur
    jumpInd = np.where(dt != 0)[0]

    tDict = {'startTime':np.nan*np.ones(len(jumpInd)+1, dtype = int), 
        'endTime':np.nan*np.ones(len(jumpInd)+1, dtype = int), 
        'duration':np.nan*np.ones(len(jumpInd)+1, dtype = int)}

    # Inetrate over every instance where the index jumped.
    for i in range(len(jumpInd)):
        tDict['startTime'][i+1] = jumpInd[i]+1
        tDict['endTime'][i] = jumpInd[i]+1
    
    tDict['startTime'][0] = 0
    tDict['endTime'][-1] = len(t)
    tDict['duration'] = np.subtract(tDict['endTime'], tDict['startTime'])

    return tDict


if __name__ == '__main__':
    # Create sythetic time stamps
    startTime = datetime(2017, 6, 6, 12)
    # Create the deltaT array that is completely periodic
    secArr = np.arange(20).astype(np.float_)

    # Randomy bump up array N times
    N = 5
    for i in range(N):
        # Randomly choose an indicie to bump up
        randomIdx = random.randint(0, 20)
        secArr[randomIdx:] += random.randint(10, 20)

    secArr[1:] += random.randint(10, 20)
    secArr[-1:] += random.randint(10, 20)

    tArr = np.array([startTime + timedelta(seconds = 2*i) for i in secArr])

    tDict = time_continuity(tArr)
    print(tArr)
    print(tDict)

    for i in range(len(tDict['startTime'])):
        print(tArr[tDict['startTime'][i].astype(int):tDict['endTime'][i].astype(int)] ,'\n')

##    t = np.arange(0.0, 2.0, 0.01)
##    s = 1 + np.sin(2*np.pi*t)
##    plt.plot(t, s)

##    plt.xlabel('time (s)')
##    plt.ylabel('voltage (mV)')
##    plt.title('About as simple as it gets, folks')
##    plt.grid(True)
##    plt.savefig("test.png")
##    plt.show()
