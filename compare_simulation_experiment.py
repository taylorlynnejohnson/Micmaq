import numpy as np
import matplotlib.pyplot as plt

def compare_sim_exp(fileNameSim='TCF-8-8-8-7day-big', fileNameExp='12-15-2022-twelvehours', lastPhase=168):
    # Compare simulated and experimental temperature data
    fileName = f"{fileNameSim}-{lastPhase}"
    
    # Load simulation data
    data_sim = np.load(fileName + '.npz')  # Assuming data is saved as NumPy compressed file
    
    # Sample times in hours and locations in meters
    sampleTimes = np.arange(0, lastPhase + 1)  # in hours
    sampleLocs = np.array([10, 16.25, 22.5, 28.75, 35.0]) * 0.0254  # inches to meters
    
    # Plot Temperature vs. Radius at Various Times
    figureName = f"Temperature Vs. Radius at Various Times - {fileNameSim} vs. {fileNameExp}"
    plt.figure(figsize=(12, 8), num=figureName)
    plt.xlabel('Radius (m)')
    plt.ylabel('Temperature (C)')
    plt.title(figureName)
    
    # Plot experimental phases
    colorMap = ['k', 'r', 'r', 'r', 'r', 'g', 'g', 'g', 'g', 'b', 'b', 'b', 'b']
    plot_phases_experimental(fileNameExp, sampleTimes, sampleLocs, colorMap)
    plot_phases_simulation(fileNameSim, lastPhase)
    
    # Plot Temperature vs. Time at Five Radial Locations
    plt.figure(figsize=(10, 6), num='Time History at Five Specified Radial Locations')
    plt.xlabel('Time (hours)')
    plt.ylabel('Temperature (C)')
    plt.title('Temperature Vs. Time at Five Radial Locations')
    plt.grid(True)
    
    # Plot simulation time history
    colorMap = ['r', 'k', 'c', 'b', 'g']
    LineStyle = '-'
    plot_time_hist_simulation(fileNameSim, lastPhase, colorMap, LineStyle)
    
    # Plot experimental time history
    LineStyle = ':'
    plot_time_hist_experimental(fileNameExp, colorMap, LineStyle)
    
    plt.show()

def plot_phases_experimental(fileNameExp, sampleTimes, sampleLocs, colorMap):
    # Load experimental data and plot for given sample times and locations
    ExpData = pd.read_csv(fileNameExp + '.txt', delimiter='\t')
    times = ExpData.iloc[1:, 0] / 60  # in hours
    times = times - times[0]
    data1 = ExpData.iloc[1:, 6:11]
    
    for iTime in range(len(sampleTimes)):
        ind = np.searchsorted(times, sampleTimes[iTime])
        if ind < len(times):
            plt.plot(sampleLocs, data1.iloc[ind, :], linestyle='--', linewidth=2, marker='o', color=colorMap[iTime], label=f'{sampleTimes[iTime]} hours')
    
    plt.legend()

def plot_phases_simulation(fileNameSim, lastPhase):
    # Load simulation data and plot for given sample times and locations
    for iPhase in range(1, lastPhase + 1):
        fileName = f"{fileNameSim}-{iPhase}"
        data = np.load(fileName + '.npz')
        Ts = data['Ts']
        r = data['r']
        M_dot = data['M_dot']
        time = data['time']
        color = 'g' if M_dot == 0 else ('r' if M_dot > 0 else 'b')
        plt.plot(r, Ts, '-', linewidth=1, color=color, label=f'{round(time / 3600, 1)} hours')
    plt.legend(loc='eastoutside')

def plot_time_hist_simulation(fileNameSim, lastPhase, colorMap, LineStyle):
    # Load simulation time history data and plot for specified radial locations
    fileName = f"{fileNameSim}-{lastPhase}"
    data = np.load(fileName + '.npz')
    timeHistory = data['timeHistory']
    
    for iHist in range(1, timeHistory.shape[1] // 2):
        plt.plot(timeHistory[:, 0] / (3600 * 24), timeHistory[:, 1 + 2 * iHist], linewidth=2, color=colorMap[iHist - 1], linestyle=LineStyle)
    
    legendLabels = [
        'Inner Wall, Radius = {:.2f} m'.format(data['tracePositions'][0]),
        'Radius = {:.2f} m'.format(data['tracePositions'][1]),
        'Radius = {:.2f} m'.format(data['tracePositions'][2]),
        'Radius = {:.2f} m'.format(data['tracePositions'][3])
    ]
    plt.legend(legendLabels, loc='eastoutside')

def plot_time_hist_experimental(fileNameExp, colorMap, LineStyle):
    # Load experimental time history data and plot for specified radial locations
    ExpData = pd.read_csv(fileNameExp + '.txt', delimiter='\t')
    times = ExpData.iloc[1:, 0] / 60  # in hours
    times = times - times[0]
    data1 = ExpData.iloc[1:, 6:11]
    
    for iHist in range(data1.shape[1]):
        plt.plot(times, data1.iloc[:, iHist], linewidth=3, color=colorMap[iHist], linestyle=LineStyle)
    plt.legend(['A1', 'A2', 'A3', 'A4', 'A5'])


