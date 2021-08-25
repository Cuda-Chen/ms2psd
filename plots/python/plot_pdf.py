import numpy as np
import matplotlib.pyplot as plt
from obspy.imaging.cm import pqlx

def plot_pdf(data_file, periods_file, lownoise_file, highnoise_file):
    data = np.loadtxt(data_file)
    periods = np.loadtxt(periods_file)
    lownoise = np.loadtxt(lownoise_file)
    highnoise = np.loadtxt(highnoise_file)
    y = np.linspace(-50, -200, data.shape[0])

    ax = plt.subplot()

    # Plot low and high noise
    ax.plot(lownoise[:, 0], lownoise[:, 1], '0.4', linewidth=2, zorder=10)
    ax.plot(highnoise[:, 0], highnoise[:, 1], '0.4', linewidth=2, zorder=10)

    # Plot PDF
    meshgrid = np.meshgrid(periods, y)
    ppsd = ax.pcolormesh(meshgrid[0], meshgrid[1], data, cmap=pqlx, zorder =-1)
    plt.draw()

    ax.semilogx()

    cb = plt.colorbar(ppsd, ax=ax)
    ax.set_xlim(2 * 1e-2, 1e2)

    plt.show()

if __name__ == '__main__':
    data_file = "../../pdf_out.txt"
    periods_file = "../../center_periods_out.txt"
    lownoise_file = "../../lownoise.mod"
    highnoise_file = "../../highnoise.mod"

    plot_pdf(data_file, periods_file, lownoise_file, highnoise_file)
