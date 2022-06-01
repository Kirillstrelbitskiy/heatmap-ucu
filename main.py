"""
A module for demonstration of heat equation work.
"""

import numpy as np
import matplotlib.pyplot as plt


class HeatMap:
    """Class for manipulating with Heating Map"""

    def __init__(self, alpha, delta_x, max_iter_time, plate_length):
        """Initializing with default params"""

        self.delta_x = delta_x

        self.delta_t = (delta_x ** 2)/(4 * alpha)
        self.gamma = (alpha * self.delta_t) / (delta_x ** 2)

        self.max_iter_time = max_iter_time
        self.plate_length = plate_length

        self.field = None

    def create_field(self):
        """Filing the map with default value"""

        self.field = np.empty(
            (self.max_iter_time, self.plate_length, self.plate_length))

        u_initial = 0
        self.field.fill(u_initial)

    def set_borders(self, u_initial, borders):
        """Setting values for borders"""

        (u_top, u_left, u_bottom, u_right) = borders

        self.field[0, :, :] = u_initial

        self.field[:, (self.plate_length-1):, :] = u_top
        self.field[:, :, :1] = u_left
        self.field[:, :1, 1:] = u_bottom
        self.field[:, :, (self.plate_length-1):] = u_right

    def calculate(self):
        """Making main calculation using heat equation"""

        for k in range(0, self.max_iter_time-1, 1):
            for i in range(1, self.plate_length-1, self.delta_x):
                for j in range(1, self.plate_length-1, self.delta_x):
                    self.field[k + 1, i, j] = self.gamma * (self.field[k][i+1][j] +
                                                            self.field[k][i-1][j] +
                                                            self.field[k][i][j+1] +
                                                            self.field[k][i][j-1] -
                                                            4*self.field[k][i][j]) 
                    self.field[k + 1, i, j] += self.field[k][i][j]

    def create_plot(self, k):
        """Creating a plot for displaying"""

        plt.clf()

        plt.title(f"Temperature at t = {k*self.delta_t:.3f} unit time")
        plt.xlabel("x")
        plt.ylabel("y")

        plt.pcolormesh(self.field[k], cmap=plt.cm.jet, vmin=0, vmax=100)
        plt.colorbar()

    def animate(self, step=1, delay=0.01):
        """Animation of heating process"""

        plt.ion()

        while True:
            for time_val in range(0, self.max_iter_time, step):
                self.create_plot(time_val)

                plt.draw()
                plt.pause(delay)


def main():
    """Main function for running class methods"""

    hmp = HeatMap(
        plate_length=50,
        max_iter_time=500,
        alpha=2,
        delta_x=1)

    hmp.create_field()
    hmp.set_borders(0, (100, 70, 20, 10))

    hmp.calculate()
    hmp.animate(5, 0.01)


if __name__ == "__main__":
    main()
