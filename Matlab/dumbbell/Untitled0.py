# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline
from scipy import ndimage

# <codecell>

a = imread('cam1.10000').astype(np.ubyte)

# <codecell>

imshow(a,cmap=cm.gray)

# <codecell>

b = where(a>50,1,0).astype(np.ubyte)

# <codecell>

struct = array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
struct = ndimage.iterate_structure(struct,2)
c = ndimage.binary_dilation(b,structure=struct,iterations=5)
c = ndimage.binary_fill_holes(c,structure=struct)

# <codecell>

imshow(np.c_[b,c],cmap=cm.gray)

# <codecell>

n = 10
l = 256
img = ndimage.gaussian_filter(a, sigma=l/(4.*n))
# mask = (im > im.mean()).astype(np.float)
# mask += 0.1 * im
# img = mask + 0.2*np.random.randn(*mask.shape)
# hist, bin_edges = np.histogram(img, bins=60)
# bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

from skimage import filter
val = filter.threshold_otsu(img)

binary_img = img > val

import matplotlib.pyplot as plt

plt.subplot(131)
plt.imshow(img)
plt.axis('off')
plt.subplot(132)
plt.plot(bin_centers, hist, lw=2)
plt.axvline(0.5, color='r', ls='--', lw=2)
plt.text(0.57, 0.8, 'histogram', fontsize=20, transform = plt.gca().transAxes)
plt.yticks([])
plt.subplot(133)
plt.imshow(binary_img, cmap=plt.cm.gray, interpolation='nearest')
plt.axis('off')

plt.subplots_adjust(wspace=0.02, hspace=0.3, top=1, bottom=0.1, left=0, right=1)
plt.show()

# <codecell>


# <codecell>

open_img = ndimage.binary_opening(binary_img,struct)
close_img = ndimage.binary_closing(open_img,struct)

# <codecell>

imshow(binary_img)

# <codecell>

imshow(open_img)

# <codecell>

imshow(close_img)

# <codecell>

from skimage import filter
val = filter.threshold_otsu(img)

# <codecell>

val

# <codecell>

labeled_array, num_features = ndimage.measurements.label(c)

# <codecell>

num_features

# <codecell>

labeled_array
imshow(labeled_array)

# <markdowncell>

#     imshow(labeled_array)

# <codecell>

ndimage.measurements.center_of_mass(labeled_array[0])

# <codecell>

labeled_array

# <codecell>

f = ndimage.find_objects(labeled_array)

# <codecell>

f

# <codecell>

centers = ndimage.measurements.center_of_mass(c,labeled_array,[num_features,1])

# <codecell>

centers[0][1]

# <codecell>

y

# <codecell>

imshow(labeled_array)

# <codecell>

centers

# <codecell>

centers

# <codecell>

x = [c[0] for c in centers]
y = [c[1] for c in centers]

# <codecell>

x,y

# <codecell>

x

# <codecell>

y

# <codecell>

a

# <codecell>

where?

# <codecell>

b = where(a>40,a,0)

# <codecell>

a.dtype

# <codecell>

b.dtype

# <codecell>

imshow(a)

# <codecell>

b = where(a>205,a,0)

# <codecell>

b.dtype

# <codecell>


# <codecell>


# <codecell>

imshow(a)

# <codecell>

imshow(b)

# <codecell>

b = ndimage.gaussian_filter(a,20)

# <codecell>

imshow(np.c_[a,b])

# <codecell>

b = where(b>10,b,0)

# <codecell>

imshow(b)

# <codecell>


