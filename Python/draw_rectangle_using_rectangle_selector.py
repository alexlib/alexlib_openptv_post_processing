

from matplotlib.widgets import RectangleSelector
import numpy as np
import matplotlib.pyplot as plt

def line_select_callback(eclick, erelease):
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print "(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2)
    print " The button you used were: ", eclick.button, erelease.button

    line, = plt.plot([eclick.xdata, erelease.xdata, erelease.xdata,
                      eclick.xdata, eclick.xdata],
                     [eclick.ydata, eclick.ydata, erelease.ydata,
                      erelease.ydata, eclick.ydata],
                     color='blue')
    rectangles.append(line)
    plt.draw()

def key_press_event_fct(event):
    """ """
    global rectangles, current_ax

    print ' Key >%s< pressed.' % (event.key)

    # delete previously selected rectangle (last item in list 'rectangles')
    if event.key in ['D', 'd']:
        if len(rectangles) > 0: # if there exist a selected rectangle
            current_ax.lines.remove(rectangles.pop(-1))
            plt.draw()

    # toggling rectangle-selector activity
    if event.key in ['Q', 'q'] and key_press_event_fct.RS.active:
        print ' RectangleSelector deactivated.'
        key_press_event_fct.RS.set_active(False)
    if event.key in ['A', 'a'] and not key_press_event_fct.RS.active:
        print ' RectangleSelector activated.'
        key_press_event_fct.RS.set_active(True)
        
        

current_ax = plt.subplot(111)                    # make a new plotingrange
N = 100000                                       # If N is large one can see 
x = np.linspace(0.0, 10.0, N)                    # improvement by use blitting!

plt.plot(x, +np.sin(.2*np.pi*x), lw=3.5, c='b', alpha=.7)  # plot something
plt.plot(x, +np.cos(.2*np.pi*x), lw=3.5, c='r', alpha=.5)
plt.plot(x, -np.sin(.2*np.pi*x), lw=3.5, c='g', alpha=.3)

rectangles = [] # empty list to hold rectangles

print "\n      click  -->  release to select a rectangle"
print " key 'd' to delete previously selected rectangle "
print " key 'a' to activate rectangle selector "
print " key 'q' to deactivate rectangle selector "

# drawtype is 'box' or 'line' or 'none'
key_press_event_fct.RS = RectangleSelector(current_ax, line_select_callback,
                                       drawtype='box', useblit=True,
                                       minspanx=5, minspany=5,
                                       spancoords='pixels')
plt.connect('key_press_event', key_press_event_fct)
plt.show()
