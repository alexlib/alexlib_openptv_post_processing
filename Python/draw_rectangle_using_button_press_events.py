import matplotlib.pyplot as plt
import numpy as np

def onpress(event):
    """ called after a button_press_event """

    global prev_point, rectangles

    if event.button == 1: # left mouse button
        if len(prev_point) == 0:
            prev_point.append((event.xdata, event.ydata))
        elif len(prev_point) == 1:
            line, = plt.plot([prev_point[0][0], event.xdata, event.xdata,
                              prev_point[0][0], prev_point[0][0]],
                             [prev_point[0][1], prev_point[0][1], event.ydata,
                              event.ydata, prev_point[0][1]],
                             color='blue')
            rectangles.append(line)
            plt.draw()

            prev_point = [] # delete previously selected point

    elif event.button == 3: # right mouse button
        if len(rectangles) > 0: # if there exist a selected rectangle
            current_ax.lines.remove(rectangles.pop(-1))
            plt.draw()


prev_point = [] # list to hold previously selected point
rectangles = [] # list to hold drawn rectangles

current_ax = plt.subplot(111)
x = np.linspace(0.0, 10.0, 1000)

plt.plot(x, +np.sin(.2*np.pi*x), lw=3.5, c='b', alpha=.7)  # plot something
plt.plot(x, +np.cos(.2*np.pi*x), lw=3.5, c='r', alpha=.5)
plt.plot(x, -np.sin(.2*np.pi*x), lw=3.5, c='g', alpha=.3)

plt.connect('button_press_event', onpress)
plt.show()
