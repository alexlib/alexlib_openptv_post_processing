import numpy as np
import Image
import matplotlib.pylab as plt
from im2array import *


im = Image.open('cam1.10000')

imnp = image2array(im)
plt.figure()
plt.imshow(imnp,cmap = plt.cm.gray)

# now press 2 corners: top-left, right-bottom
x = plt.ginput(2)
# parsing
a,b=x[0]; c,d=x[1]
box1 = (a,b,c,d)

# now press again 2 corners: top-left, right-bottom
x = plt.ginput(2)
# parsing
a,b=x[0]; c,d=x[1]
box2 = (a,b,c,d)

# crop using PIL
ball1 = im.crop(box1)
ball2 = im.crop(box2)

# resize to get integer sizes for sure

ball1 = ball1.resize((int(ball1.size[0]),int(ball1.size[1])))
ball2 = ball2.resize((int(ball2.size[0]),int(ball2.size[1])))


ball1np = image2array(ball1)
ball2np = image2array(ball2)

plt.figure()
plt.imshow(ball1np,cmap=plt.cm.gray)
plt.show()
plt.figure()
plt.imshow(ball2np,cmap=plt.cm.gray)
plt.show()
