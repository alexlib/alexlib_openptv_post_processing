#!/usr/bin/env python
# -*- coding: utf-8 -*-
# tifffile.py

# Copyright (c) 2008-2011, Christoph Gohlke
# Copyright (c) 2008-2011, The Regents of the University of California
# Produced at the Laboratory for Fluorescence Dynamics
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Read and write image data from and to TIFF files.

Image and meta-data can be read from TIFF, BigTIFF, OME-TIFF, STK, LSM, NIH,
and FluoView files. Only a subset of the TIFF specification is supported,
mainly uncompressed and losslessly compressed 2**(0 to 6) bit integer,
16, 32 and 64-bit float, grayscale and RGB(A) images, which are commonly
used in bio-scientific imaging. Specifically, reading JPEG or CCITT
compressed image data is not implemented. Only primary info records are
read for STK, FluoView, and NIH image formats.

TIFF, the Tagged Image File Format, is under the control of Adobe Systems.
BigTIFF allows for files greater than 4 GB. STK, LSM, FluoView, and OME-TIFF
are custom extensions defined by MetaMorph, Carl Zeiss MicroImaging,
Olympus, and the Open Microscopy Environment consortium respectively.

The API is not stable yet and might change between revisions.
Tested on little-endian platforms only.

For command line usage run ``python tifffile.py --help``

:Authors:
  `Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>`__,
  Laboratory for Fluorescence Dynamics, University of California, Irvine

:Version: 2011.12.15

Requirements
------------

* `Python 2.7 or 3.2 <http://www.python.org>`__
* `Numpy 1.6 <http://numpy.scipy.org>`__
* `Matplotlib 1.1 <http://matplotlib.sourceforge.net>`__
  (optional for plotting)
* `tifffile.c 2011.12.15 <http://www.lfd.uci.edu/~gohlke/>`__
  (optional for faster decoding of PackBits and LZW encoded strings)

Acknowledgements
----------------
*  Egor Zindy, University of Manchester, for cz_lsm_scan_info specifics.
*  Wim Lewis, for a bug fix and some read_cz_lsm functions.

References
----------

(1) TIFF 6.0 Specification and Supplements. Adobe Systems Incorporated.
    http://partners.adobe.com/public/developer/tiff/
(2) TIFF File Format FAQ. http://www.awaresystems.be/imaging/tiff/faq.html
(3) MetaMorph Stack (STK) Image File Format.
    http://support.meta.moleculardevices.com/docs/t10243.pdf
(4) File Format Description - LSM 5xx Release 2.0.
    http://ibb.gsf.de/homepage/karsten.rodenacker/IDL/Lsmfile.doc
(5) BioFormats. http://www.loci.wisc.edu/ome/formats.html
(6) The OME-TIFF format.
    http://www.openmicroscopy.org/site/support/file-formats/ome-tiff

Examples
--------

>>> data = numpy.random.rand(301, 219)
>>> imsave('temp.tif', data)
>>> image = imread('temp.tif')
>>> assert numpy.all(image == data)

>>> tif = TIFFfile('test.tif')
>>> images = tif.asarray()
>>> image0 = tif[0].asarray()
>>> for page in tif:
...     for tag in page.tags.values():
...         t = tag.name, tag.value
...     image = page.asarray()
...     if page.is_rgb: pass
...     if page.is_palette:
...         t = page.color_map
...     if page.is_stk:
...         t = page.mm_uic_tags.number_planes
...     if page.is_lsm:
...         t = page.cz_lsm_info
>>> tif.close()

"""

from __future__ import division, print_function

import sys
import os
import math
import zlib
import time
import struct
import warnings
import datetime
import collections
from contextlib import contextmanager
from xml.etree import cElementTree as ElementTree

import numpy

__all__ = ['imsave', 'imread', 'imshow', 'tifffile', 'TIFFfile']


def imsave(filename, data, photometric=None, planarconfig=None,
           resolution=None, bigtiff=False):
    """Write image data to TIFF file.

    Image data are written uncompressed in one stripe per plane.
    Dimensions larger than 2 (grayscale) or 3 (rgb) are flattened and
    saved as separate pages.

    Arguments
    ---------

    filename : str
        Name of file to write.

    data : array like
        Input image. The last dimensions are assumed to be image height,
        width, and color channels (if rgb).

    bigtiff : bool
        If True, the BigTIFF format is used.
        Default is to auto detect based on data size.

    photometric : str {'minisblack', 'miniswhite', 'rgb'}
        Default is to autodetect based on data shape.

    planarconfig : str {'contig, 'planar'}
        Specifies if color samples are stored contiguous or in separate planes.
        Default is to autodetect based on data shape.
        'contig': last dimension contains color channels.
        'planar': third last dimension contains color channels.

    resolution : ((int, int), (int, int))
        X and Y resolution in dots per inch as rational number.

    Example
    -------

    >>> data = numpy.random.rand(3, 301, 219)
    >>> imsave('temp.tif', data)

    """
    data = numpy.asarray(data)
    data_shape = shape = data.shape
    data = numpy.atleast_2d(data)
    pack = struct.pack

    if not bigtiff and data.size * data.dtype.itemsize < 2**31 - 2**16:
        bigtiff = False
        offset_size = 4
        tag_size = 12
        numtag_format = 'H'
        offset_format = 'I'
        val_format = '4s'
    else:
        bigtiff = True
        offset_size = 8
        tag_size = 20
        numtag_format = 'Q'
        offset_format = 'Q'
        val_format = '8s'

    # unify shape of data
    extrasamples = None
    if photometric is None:
        if data.ndim > 2 and (shape[-3] in (3, 4) or shape[-1] in (3, 4)):
            photometric = 'rgb'
        else:
            photometric = 'minisblack'
    if photometric == 'rgb':
        if planarconfig is None:
            if shape[-1] in (3, 4):
                planarconfig = 'contig'
                extrasamples = shape[-1] == 4
            elif shape[-3] in (3, 4):
                planarconfig = 'planar'
                extrasamples = shape[-3] == 4
        else:
            if ((planarconfig == 'contig' and shape[-1] not in (3, 4)) or
                (planarconfig == 'planar' and shape[-3] not in (3, 4))):
                raise ValueError("photometric or planarconfig arguments "
                                 "don't match shape %s" % str(shape))
        if planarconfig == 'contig':
            data = data.reshape((-1, 1) + shape[-3:])
        else:
            data = data.reshape((-1, ) + shape[-3:] + (1, ))
    else:
        planarconfig = None
        data = data.reshape((-1, 1) + shape[-2:] + (1, ))

    shape = data.shape  # (pages, planes, height, width, contig samples)

    bytestr = bytes if sys.version[0] == '2' else lambda x: bytes(x, 'ascii')
    dtypes = {'B': 1, 's': 2, 'H': 3, 'I': 4, '2I': 5, 'b': 6,
              'h': 8, 'i': 9, 'f': 11, 'd': 12, 'Q': 16, 'q': 17}
    tags = []
    tag_data = []

    def tag(name, dtype, number, value, offset=[0]):
        # append tag binary string to tags
        # append (offset, value as binary string) to tag_data
        # increment offset by tag_size
        if dtype == 's':
            value = bytestr(value) + b'\0'
            number = len(value)
            value = (value, )
        t = [pack('HH', TIFF_TAGS_W[name], dtypes[dtype]),
             pack(offset_format, number)]
        if len(dtype) > 1:
            number *= int(dtype[:-1])
            dtype = dtype[-1]
        if number == 1:
            if isinstance(value, (tuple, list)):
                value = value[0]
            t.append(pack(val_format, pack(dtype, value)))
        elif struct.calcsize(dtype) * number <= offset_size:
            t.append(pack(val_format, pack(str(number)+dtype, *value)))
        else:
            t.append(pack(offset_format, 0))
            tag_data.append((offset[0] + offset_size + 4,
                             pack(str(number)+dtype, *value)))
        tags.append(b''.join(t))
        offset[0] += tag_size

    tag('software', 's', 0, "tifffile.py")
    tag('datetime', 's', 0,
        datetime.datetime.now().strftime("%Y:%m:%d %H:%M:%S"))
    tag('compression', 'H', 1, 1)
    tag('orientation', 'H', 1, 1)
    tag('image_width', 'I', 1, shape[-2])
    tag('image_length', 'I', 1, shape[-3])
    tag('new_subfile_type', 'I', 1, 0 if shape[0] == 1 else 2)
    tag('sample_format', 'H', 1,
        {'u': 1, 'i': 2, 'f': 3, 'c': 6}[data.dtype.kind])
    tag('photometric', 'H', 1,
        {'miniswhite': 0, 'minisblack': 1, 'rgb': 2}[photometric])
    if shape != data_shape:
        tag('image_description', 's', 0,
            "shape=(%s)" % (",".join('%i' % i for i in data_shape)))
    if planarconfig:
        tag('planar_configuration', 'H', 1, 1 if planarconfig=='contig' else 2)
        if extrasamples:
            tag('samples_per_pixel', 'H', 1, 4)
            tag('extra_samples', 'H', 1, 1)
            tag('bits_per_sample', 'H', 4, (data.dtype.itemsize * 8, ) * 4)
        else:
            tag('samples_per_pixel', 'H', 1, 3)
            tag('bits_per_sample', 'H', 3, (data.dtype.itemsize * 8, ) * 3)
    else:
        tag('samples_per_pixel', 'H', 1, 1)
        tag('bits_per_sample', 'H', 1, data.dtype.itemsize * 8)
    if resolution:
        tag('x_resolution', '2I', 1, resolution[0])
        tag('y_resolution', '2I', 1, resolution[1])
        tag('resolution_unit', 'H', 1, 2)
    tag('rows_per_strip', 'I', 1, shape[-3])
    # use one strip per plane
    strip_byte_counts = (data[0, 0].size * data.dtype.itemsize, ) * shape[1]
    tag('strip_byte_counts', offset_format, shape[1], strip_byte_counts)
    # strip_offsets must be the last tag; will be updated later
    tag('strip_offsets', offset_format, shape[1], (0, ) * shape[1])
    numtags = len(tags)
    tags = b''.join(tags)

    fhandle = open(filename, 'wb')
    seek = fhandle.seek
    tell = fhandle.tell

    def write(arg, *args):
        fhandle.write(pack(arg, *args) if args else arg)

    write(b'II' if sys.byteorder == 'little' else b'MM')
    if bigtiff:
        write('HHH', 43, 8, 0)
    else:
        write('H', 42)
    ifd_offset = tell()
    write(offset_format, 0)  # first IFD
    for i in range(shape[0]):
        # update pointer at ifd_offset
        pos = tell()
        seek(ifd_offset)
        write(offset_format, pos)
        seek(pos)
        # write tags
        write(numtag_format, numtags)
        tag_offset = tell()
        write(tags)
        ifd_offset = tell()
        write(offset_format, 0)  # offset to next ifd
        # write extra tag data and update pointers
        for off, dat in tag_data:
            pos = tell()
            seek(tag_offset + off)
            write(offset_format, pos)
            seek(pos)
            write(dat)
        # update strip_offsets
        pos = tell()
        if len(strip_byte_counts) == 1:
            seek(ifd_offset - offset_size)
            write(offset_format, pos)
        else:
            seek(pos - offset_size*shape[1])
            strip_offset = pos
            for size in strip_byte_counts:
                write(offset_format, strip_offset)
                strip_offset += size
        seek(pos)
        # write data
        data[i].tofile(fhandle)  # if this fails, try update Python or numpy
        fhandle.flush()
    fhandle.close()


def imread(filename, *args, **kwargs):
    """Return image data from TIFF file as numpy array.

    The first image series is returned if no arguments are provided.

    Arguments
    ---------

    key : int, slice, or sequence of page indices
        Defines which pages to return as array.

    series : int
        Defines which series of pages to return as array.

    Example
    -------

    >>> image = imread('test.tif', 0)

    """
    with tifffile(filename) as tif:
        return tif.asarray(*args, **kwargs)


@contextmanager
def tifffile(filename):
    """Return TIFFfile instance. Support for `with` statement.

    Example
    -------

    >>> with tifffile('test.tif') as tif:
    ...    image = tif.asarray()

    """
    tiff_file = TIFFfile(filename)
    try:
        yield tiff_file
    finally:
        tiff_file.close()


class lazyattr(object):
    """Lazy object attribute whose value is computed on first access."""
    def __init__(self, func):
        self.func = func

    def __get__(self, instance, owner):
        if instance is None:
            return self
        result = self.func(instance)
        if result is NotImplemented:
            return getattr(super(owner, instance), self.func.__name__)
        setattr(instance, self.func.__name__, result)
        return result


class TIFFfile(object):
    """Read image and meta-data from TIFF, STK, LSM, and FluoView files.

    Attributes
    ----------

    pages : list
        All TIFFpages in file.

    series : list of Records(shape, dtype, axes, TIFFpages)
        TIFF pages with compatible shapes and types.

    All attributes are read-only.

    Example
    -------

    >>> tif = TIFFfile('test.tif')
    ... try:
    ...     images = tif.asarray()
    ... except Exception as e:
    ...     print(e)
    ... finally:
    ...     tif.close()

    """
    def __init__(self, filename):
        """Initialize instance from file."""
        filename = os.path.abspath(filename)
        self.fhandle = open(filename, 'rb')
        self.fname = os.path.basename(filename)
        self.fpath = os.path.dirname(filename)
        self.fstat = os.fstat(self.fhandle.fileno())
        self.tiffs = {self.fname: self}  # cache of TIFFfiles
        self.pages = []
        self.offset_size = 4
        try:
            self._fromfile()
        except Exception:
            self.fhandle.close()
            raise

    def close(self):
        """Close open file handle(s)."""
        if not hasattr(self, 'tiffs'):
            return
        for tif in self.tiffs.values():
            if tif.fhandle:
                tif.fhandle.close()
                tif.fhandle = None

    def _fromfile(self):
        """Read TIFF header and all page records from file."""
        self.fhandle.seek(0)
        try:
            self.byte_order = {b'II': '<', b'MM': '>'}[self.fhandle.read(2)]
        except KeyError:
            raise ValueError("not a valid TIFF file")
        version = struct.unpack(self.byte_order+'H', self.fhandle.read(2))[0]
        if version == 43:  # BigTiff
            self.offset_size, zero = struct.unpack(self.byte_order+'HH',
                                                   self.fhandle.read(4))
            if zero or self.offset_size != 8:
                raise ValueError("not a valid BigTIFF file")
        elif version != 42:
            raise ValueError("not a TIFF file")
        self.pages = []
        while True:
            try:
                page = TIFFpage(self)
                self.pages.append(page)
            except StopIteration:
                break
        if not self.pages:
            raise ValueError("empty TIFF file")

    @lazyattr
    def series(self):
        """Return series of TIFFpage with compatible shape and properties."""
        if self.is_ome:
            series = self._omeseries()
        elif self.is_fluoview:
            dims = {b'X': 'X', b'Y': 'Y', b'Z': 'Z', b'T': 'T',
                    b'WAVELENGTH': 'C', b'TIME': 'T', b'XY': 'R',
                    b'EVENT': 'E', b'EXPOSURE': 'O'}
            mmhd = list(reversed(self.pages[0].mm_header.dimensions))
            series = [Record(
                axes=''.join(dims.get(i[0].strip().upper(), '?')
                             for i in mmhd if i[1] > 1),
                shape=tuple(int(i[1]) for i in mmhd if i[1] > 1),
                pages=self.pages, dtype=numpy.dtype(self.pages[0].dtype))]
        elif self.is_lsm:
            lsmi = self.pages[0].cz_lsm_info
            axes = CZ_SCAN_TYPES[lsmi.scan_type]
            if self.pages[0].is_rgb:
                axes = axes.replace('C', '').replace('XY', 'XYC')
            axes = axes[::-1]
            shape = [getattr(lsmi, CZ_DIMENSIONS[i]) for i in axes]
            pages = [p for p in self.pages if not p.is_reduced]
            series = [Record(axes=axes, shape=shape, pages=pages,
                             dtype=numpy.dtype(pages[0].dtype))]
            if len(pages) != len(self.pages):  # reduced RGB pages
                pages = [p for p in self.pages if p.is_reduced]
                cp = 1
                i = 0
                while cp < len(pages) and i < len(shape)-2:
                    cp *= shape[i]
                    i += 1
                shape = shape[:i] + list(pages[0].shape)
                axes = axes[:i] + 'CYX'
                series.append(Record(axes=axes, shape=shape, pages=pages,
                                     dtype=numpy.dtype(pages[0].dtype)))
        elif self.is_nih:
            series = [Record(pages=self.pages,
                             shape=(len(self.pages),) + self.pages[0].shape,
                             axes='I' + self.pages[0].axes,
                             dtype=numpy.dtype(self.pages[0].dtype))]
        elif self.pages[0].is_shaped:
            shape = self.pages[0].tags['image_description'].value[7:-1]
            shape = tuple(int(i) for i in shape.split(b','))
            series = [Record(pages=self.pages, shape=shape,
                             axes='?' * len(shape),
                             dtype=numpy.dtype(self.pages[0].dtype))]
        else:
            shapes = []
            pages = {}
            for page in self.pages:
                shape = page.shape + (page.axes,
                                      page.compression in TIFF_DECOMPESSORS)
                if not shape in pages:
                    shapes.append(shape)
                    pages[shape] = [page]
                else:
                    pages[shape].append(page)
            series = [Record(pages=pages[s],
                             axes=(('I' + s[-2])
                                   if len(pages[s]) > 1 else s[-2]),
                             dtype=numpy.dtype(pages[s][0].dtype),
                             shape=((len(pages[s]), ) + s[:-2]
                                    if len(pages[s]) > 1 else s[:-2]))
                      for s in shapes]
        return series

    def asarray(self, key=None, series=None):
        """Return image data of multiple TIFF pages as numpy array.

        By default the first image series is returned.

        Arguments
        ---------

        key : int, slice, or sequence of page indices
            Defines which pages to return as array.

        series : int
            Defines which series of pages to return as array.

        """
        if key is None and series is None:
            series = 0
        if series is not None:
            pages = self.series[series].pages
        else:
            pages = self.pages

        if key is None:
            pass
        elif isinstance(key, int):
            pages = [pages[key]]
        elif isinstance(key, slice):
            pages = pages[key]
        elif isinstance(key, collections.Iterable):
            pages = [pages[k] for k in key]
        else:
            raise TypeError('key must be an int, slice, or sequence')

        if len(pages) == 1:
            result = pages[0].asarray()
        elif self.is_nih:
            result = numpy.vstack(p.asarray() for p in pages)
            if pages[0].is_palette:
                result = numpy.take(pages[0].color_map, result, axis=1)
                result = numpy.swapaxes(result, 0, 1)
        else:
            if self.is_ome and any(p is None for p in pages):
                firstpage = next(p for p in pages if p)
                nopage = numpy.zeros_like(
                    firstpage.asarray())
            result = numpy.vstack((p.asarray() if p else nopage)
                                  for p in pages)
            if key is None:
                result.shape = self.series[series].shape
        return result

    def _omeseries(self):
        """Return image series in OME-TIFF files."""
        root = ElementTree.XML(self.pages[0].tags['image_description'].value)
        uuid = root.attrib.get('UUID', None)
        self.tiffs = {uuid: self}
        result = []
        for element in root:
            if element.tag.endswith('BinaryOnly'):
                warnings.warn("not an OME-TIFF master file")
                break
            if not element.tag.endswith('Image'):
                continue
            for pixels in element:
                if not pixels.tag.endswith('Pixels'):
                    continue
                atr = pixels.attrib
                axes = "".join(reversed(atr['DimensionOrder']))
                shape = tuple(int(atr['Size'+ax]) for ax in axes)
                size = numpy.prod(shape[:-2])
                ifds = [None] * size
                for data in pixels:
                    if not data.tag.endswith('TiffData'):
                        continue
                    atr = data.attrib
                    ifd = int(atr.get('IFD', 0))
                    num = int(atr.get('NumPlanes', 1 if 'IFD' in atr else 0))
                    num = int(atr.get('PlaneCount', num))
                    idx = [int(atr.get('First'+ax, 0)) for ax in axes[:-2]]
                    idx = numpy.ravel_multi_index(idx, shape[:-2])
                    for uuid in data:
                        if uuid.tag.endswith('UUID'):
                            if uuid.text not in self.tiffs:
                                fn = uuid.attrib['FileName']
                                try:
                                    tf = TIFFfile(os.path.join(self.fpath, fn))
                                except (IOError, ValueError):
                                    warnings.warn("failed to read %s" % fn)
                                    break
                                self.tiffs[uuid.text] = tf
                            pages = self.tiffs[uuid.text].pages
                            for i in range(num if num else len(pages)):
                                ifds[idx + i] = pages[ifd + i]
                            break
                    else:
                        pages = self.pages
                        for i in range(num if num else len(pages)):
                            ifds[idx + i] = pages[ifd + i]
                result.append(Record(axes=axes, shape=shape, pages=ifds,
                                     dtype=numpy.dtype(ifds[0].dtype)))
        return result

    def __del__(self):
        """Close open files."""
        return self.close()

    def __len__(self):
        """Return number of image pages in file."""
        return len(self.pages)

    def __getitem__(self, key):
        """Return specified page."""
        return self.pages[key]

    def __iter__(self):
        """Return iterator over pages."""
        return iter(self.pages)

    def __str__(self):
        """Return string containing information about file."""
        result = [self.fname.capitalize(),
            "%.2f MB" % (self.fstat[6] / 1048576),
            {'<': 'little endian', '>': 'big endian'}[self.byte_order]]
        if self.is_bigtiff:
            result.append("bigtiff")
        if len(self.pages) > 1:
            result.append("%i pages" % len(self.pages))
        if len(self.series) > 1:
            result.append("%i series" % len(self.series))
        if len(self.tiffs) > 1:
            result.append("%i files" % (len(self.tiffs)))
        return ", ".join(result)

    @lazyattr
    def is_bigtiff(self):
        return self.offset_size != 4

    @lazyattr
    def is_rgb(self):
        return all(p.is_rgb for p in self.pages)

    @lazyattr
    def is_palette(self):
        return all(p.is_palette for p in self.pages)

    @lazyattr
    def is_stk(self):
        return all(p.is_stk for p in self.pages)

    @lazyattr
    def is_lsm(self):
        return self.pages[0].is_lsm

    @lazyattr
    def is_nih(self):
        return self.pages[0].is_nih

    @lazyattr
    def is_fluoview(self):
        return self.pages[0].is_fluoview

    @lazyattr
    def is_ome(self):
        return self.pages[0].is_ome


class TIFFpage(object):
    """A TIFF image file directory (IFD).

    Attributes
    ----------

    index : int
        Index of page in file.

    dtype : str {TIFF_SAMPLE_DTYPES}
        Data type of image, colormapped if applicable.

    shape : tuple
        Dimensions of the image array in TIFF page,
        colormapped and with one alpha channel if applicable.

    axes : str
        Axes label codes:
        'X' width, 'Y' height, 'S' sample, 'P' plane,
        'Z' depth, 'I' image series, 'C' color|wavelength|channel,
        'T' time, 'R' region|tile, '?' unknown, '_' missing

    tags : TiffTags
        Dictionary of tags in page.
        Tag values are also directly accessible as attributes.

    color_map : numpy array
        Color look up table if exists.

    mm_uic_tags: Record(dict)
        Consolidated MetaMorph mm_uic# tags, if exists.

    cz_lsm_scan_info: Record(dict)
        LSM scan info attributes, if exists.

    All attributes are read-only.

    """
    def __init__(self, parent):
        """Initialize instance from file."""
        self.parent = parent
        self.index = len(parent.pages)
        self.shape = self._shape = ()
        self.dtype = self._dtype = None
        self.axes = ""
        self.tags = TiffTags()

        self._fromfile()
        self._process_tags()

    def _fromfile(self):
        """Read TIFF IFD structure and its tags from file.

        File cursor must be at storage position of IFD offset and is left at
        offset to next IFD.

        Raises StopIteration if offset (first bytes read) is 0.

        """
        fhandle = self.parent.fhandle
        byte_order = self.parent.byte_order
        offset_size = self.parent.offset_size

        fmt = {4: 'I', 8: 'Q'}[offset_size]
        offset = struct.unpack(byte_order + fmt, fhandle.read(offset_size))[0]
        if not offset:
            raise StopIteration()

        # read standard tags
        tags = self.tags
        fhandle.seek(offset, 0)
        fmt, size = {4: ('H', 2), 8: ('Q', 8)}[offset_size]
        try:
            numtags = struct.unpack(byte_order + fmt, fhandle.read(size))[0]
        except Exception:
            warnings.warn("corrupted page list")
            raise StopIteration()

        for _ in range(numtags):
            tag = TIFFtag(self.parent)
            tags[tag.name] = tag

        # read custom tags
        for name, readtag in CUSTOM_TAGS.values():
            if name in tags and readtag:
                pos = fhandle.tell()
                value = readtag(fhandle, byte_order, tags[name])
                if isinstance(value, dict):  # numpy.core.records.record
                    value = Record(value)
                tags[name].value = value
                fhandle.seek(pos)

        # read LSM info subrecords
        if self.is_lsm:
            pos = fhandle.tell()
            for name, reader in CZ_LSM_INFO_READERS.items():
                try:
                    offset = self.cz_lsm_info["offset_"+name]
                except KeyError:
                    continue
                if not offset:
                    continue
                fhandle.seek(offset)
                try:
                    setattr(self, "cz_lsm_"+name, reader(fhandle, byte_order))
                except ValueError:
                    pass
            fhandle.seek(pos)

    def _process_tags(self):
        """Validate standard tags and initialize attributes.

        Raise ValueError if tag values are not supported.

        """
        tags = self.tags
        for code, (name, default, dtype, count, validate) in TIFF_TAGS.items():
            if not (name in tags or default is None):
                tags[name] = TIFFtag(code, dtype=dtype, count=count,
                                     value=default, name=name)
            if name in tags and validate:
                try:
                    if tags[name].count == 1:
                        setattr(self, name, validate[tags[name].value])
                    else:
                        setattr(self, name, tuple(validate[value]
                                            for value in tags[name].value))
                except KeyError:
                    raise ValueError("%s.value (%s) not supported" %
                                     (name, tags[name].value))

        tag = tags['bits_per_sample']
        if tag.count == 1:
            self.bits_per_sample = tag.value
        else:
            value = tag.value[:self.samples_per_pixel]
            if any((v-value[0] for v in value)):
                self.bits_per_sample = value
            else:
                self.bits_per_sample = value[0]

        tag = tags['sample_format']
        if tag.count == 1:
            self.sample_format = TIFF_SAMPLE_FORMATS[tag.value]
        else:
            value = tag.value[:self.samples_per_pixel]
            if any((v-value[0] for v in value)):
                self.sample_format = [TIFF_SAMPLE_FORMATS[v] for v in value]
            else:
                self.sample_format = TIFF_SAMPLE_FORMATS[value[0]]

        self.strips_per_image = int(math.floor(float(self.image_length +
                            self.rows_per_strip - 1) / self.rows_per_strip))

        key = (self.sample_format, self.bits_per_sample)
        self.dtype = self._dtype = TIFF_SAMPLE_DTYPES.get(key, None)

        if self.is_stk:
            planes = tags['mm_uic2'].count
            # consolidate mm_uci tags
            self.mm_uic_tags = Record(tags['mm_uic2'].value)
            for key in ('mm_uic3', 'mm_uic4', 'mm_uic1'):
                if key in tags:
                    self.mm_uic_tags.update(tags[key].value)
            if self.planar_configuration == 'contig':
                self._shape = (planes, 1, self.image_length,
                               self.image_width, self.samples_per_pixel)
                self.shape = tuple(self._shape[i] for i in (0, 2, 3, 4))
                self.axes = "PYXS"
            else:
                self._shape = (planes, self.samples_per_pixel,
                               self.image_length, self.image_width, 1)
                self.shape = self._shape[:4]
                self.axes = "PSYX"
        elif self.is_palette:
            self.dtype = self.tags['color_map'].dtype[1]
            self.color_map = numpy.array(self.color_map, self.dtype)
            self.color_map.shape = (3, -1)
            self._shape = (1, 1, self.image_length, self.image_width, 1)
            if self.color_map.shape[1] >= 2**self.bits_per_sample:
                self.shape = (3, self.image_length, self.image_width)
                self.axes = "SYX"
            else:
                # LSM and FluoView
                self.shape = (self.image_length, self.image_width)
                self.axes = "YX"
        elif self.is_rgb:
            if self.planar_configuration == 'contig':
                self._shape = (1, 1, self.image_length, self.image_width,
                               self.samples_per_pixel)
                self.shape = (self.image_length, self.image_width, 3)
                self.axes = "YXS"
            else:
                self._shape = (1, self.samples_per_pixel, self.image_length,
                               self.image_width, 1)
                self.shape = self._shape[1:-1]
                self.axes = "SYX"
            if 'extra_samples' in self.tags:
                extra_samples = self.extra_samples
                if self.tags['extra_samples'].count == 1:
                    extra_samples = (extra_samples, )
                for i, exs in enumerate(extra_samples):
                    if exs in ('unassalpha', 'assocalpha'):
                        if self.planar_configuration == 'contig':
                            self.shape = self.shape[:2] + (4,)
                        else:
                            self.shape = (4,) + self.shape[1:]
                        break
        else:
            self._shape = (1, 1, self.image_length, self.image_width, 1)
            self.shape = self._shape[2:4]
            self.axes = "YX"

        if not self.compression and not 'strip_byte_counts' in tags:
            self.strip_byte_counts = numpy.prod(self.shape) * (
                self.bits_per_sample // 8)

    def asarray(self, squeeze=True, colormapped=True, rgbonly=True):
        """Read image data from file and return as numpy array.

        Raise ValueError if format is unsupported.
        If any argument is False, the shape of the returned array might be
        different from the page shape.

        Arguments
        ---------

        squeeze : bool
            If True all length-1 dimensions (except X and Y) are
            squeezed out from result.

        colormapped : bool
            If True color mapping is applied for palette-indexed images.

        rgbonly : bool
            If True return RGB(A) image without additional extra samples.

        """
        fhandle = self.parent.fhandle
        if not fhandle:
            raise IOError("TIFF file is not open")
        if self.dtype is None:
            raise ValueError("data type not supported: %s%i" % (
                self.sample_format, self.bits_per_sample))
        if self.compression not in TIFF_DECOMPESSORS:
            raise ValueError("cannot decompress %s" % self.compression)
        if ('ycbcr_subsampling' in self.tags and
            self.tags['ycbcr_subsampling'].value not in (1, (1, 1))):
            raise ValueError("YCbCr subsampling not supported")
        tag = self.tags['sample_format']
        if tag.count != 1 and any((i-tag.value[0] for i in tag.value)):
            raise ValueError("sample formats don't match %s" % str(tag.value))

        dtype = self._dtype
        shape = self._shape
        image_width = self.image_width
        image_length = self.image_length
        typecode = self.parent.byte_order + dtype
        bits_per_sample = self.bits_per_sample

        if self.is_tiled:
            if 'tile_offsets' in self.tags:
                byte_counts = self.tile_byte_counts
                offsets = self.tile_offsets
            else:
                byte_counts = self.strip_byte_counts
                offsets = self.strip_offsets
            tile_width = self.tile_width
            tile_length = self.tile_length
            tw = (image_width + tile_width - 1) // tile_width
            tl = (image_length + tile_length - 1) // tile_length
            shape = shape[:-3] + (tl*tile_length, tw*tile_width, shape[-1])
            tile_shape = (tile_length, tile_width, shape[-1])
            runlen = tile_width
        else:
            byte_counts = self.strip_byte_counts
            offsets = self.strip_offsets
            runlen = image_width

        try:
            offsets[0]
        except TypeError:
            offsets = (offsets, )
            byte_counts = (byte_counts, )
        if any(o < 2 for o in offsets):
            raise ValueError("corrupted file")

        if (not self.is_tiled and (self.is_stk or (not self.compression
            and bits_per_sample in (8, 16, 32, 64)
            and all(offsets[i] == offsets[i+1] - byte_counts[i]
                    for i in range(len(offsets)-1))))):
            # contiguous data
            fhandle.seek(offsets[0], 0)
            result = numpy.fromfile(fhandle, typecode, numpy.prod(shape))
            result = result.astype('=' + dtype)
        else:
            if self.planar_configuration == 'contig':
                runlen *= self.samples_per_pixel
            if bits_per_sample in (8, 16, 32, 64, 128):
                if (bits_per_sample * runlen) % 8:
                    raise ValueError("data and sample size mismatch")
                unpack = lambda x: numpy.fromstring(x, typecode)
            elif isinstance(bits_per_sample, tuple):
                unpack = lambda x: unpackrgb(x, typecode, bits_per_sample)
            else:
                unpack = lambda x: unpackints(x, typecode, bits_per_sample,
                                              runlen)
            decompress = TIFF_DECOMPESSORS[self.compression]
            if self.is_tiled:
                result = numpy.empty(shape, dtype)
                tw, tl, pl = 0, 0, 0
                for offset, bytecount in zip(offsets, byte_counts):
                    fhandle.seek(offset, 0)
                    tile = unpack(decompress(fhandle.read(bytecount)))
                    tile.shape = tile_shape
                    result[0, pl, tl:tl+tile_length,
                           tw:tw+tile_width, :] = tile
                    del tile
                    tw += tile_width
                    if tw >= shape[-2]:
                        tw, tl = 0, tl + tile_length
                        if tl >= shape[-3]:
                            tl, pl = 0, pl + 1
                result = result[..., :image_length, :image_width, :]
            else:
                result = numpy.empty(shape, dtype).reshape(-1)
                index = 0
                for offset, bytecount in zip(offsets, byte_counts):
                    fhandle.seek(offset, 0)
                    stripe = unpack(decompress(fhandle.read(bytecount)))
                    size = min(result.size, stripe.size)
                    result[index:index+size] = stripe[:size]
                    del stripe
                    index += size

        result.shape = self._shape

        if self.predictor == 'horizontal':
            numpy.cumsum(result, axis=3, dtype=dtype, out=result)

        if colormapped and self.is_palette:
            if self.color_map.shape[1] >= 2**bits_per_sample:
                # FluoView and LSM might fail here
                result = numpy.take(self.color_map, result, axis=1)
        elif rgbonly and self.is_rgb and 'extra_samples' in self.tags:
            # return only RGB and first alpha channel if exists
            extra_samples = self.extra_samples
            if self.tags['extra_samples'].count == 1:
                extra_samples = (extra_samples, )
            for i, exs in enumerate(extra_samples):
                if exs in ('unassalpha', 'assocalpha'):
                    if self.planar_configuration == 'contig':
                        result = result[..., [0, 1, 2, 3+i]]
                    else:
                        result = result[:, [0, 1, 2, 3+i]]
                    break
            else:
                if self.planar_configuration == 'contig':
                    result = result[..., :3]
                else:
                    result = result[:, :3]

        if squeeze:
            try:
                result.shape = self.shape
            except ValueError:
                pass

        return result

    def __str__(self):
        """Return string containing information about page."""
        s = ', '.join(s for s in (
            ' x '.join(str(i) for i in self.shape),
            str(numpy.dtype(self.dtype)),
            '%s bit' % str(self.bits_per_sample),
            self.photometric,
            self.compression if self.compression else 'raw',
            ','.join(t[3:] for t in ('is_stk', 'is_lsm', 'is_nih', 'is_ome',
                                     'is_fluoview', 'is_reduced', 'is_tiled')
                     if getattr(self, t))) if s)
        return "Page %i: %s" % (self.index, s)

    def __getattr__(self, name):
        """Return tag value."""
        if name in self.tags:
            value = self.tags[name].value
            setattr(self, name, value)
            return value
        raise AttributeError(name)

    @lazyattr
    def is_rgb(self):
        """True if page contains a RGB image."""
        return self.tags['photometric'].value == 2

    @lazyattr
    def is_palette(self):
        """True if page contains a palette-colored image."""
        return self.tags['photometric'].value == 3

    @lazyattr
    def is_tiled(self):
        """True if page contains tiled image."""
        return 'tile_width' in self.tags

    @lazyattr
    def is_reduced(self):
        """True if page is a reduced image of another image."""
        return bool(self.tags['new_subfile_type'].value & 1)

    @lazyattr
    def is_stk(self):
        """True if page contains MM_UIC2 tag."""
        return 'mm_uic2' in self.tags

    @lazyattr
    def is_lsm(self):
        """True if page contains LSM CZ_LSM_INFO tag."""
        return 'cz_lsm_info' in self.tags

    @lazyattr
    def is_fluoview(self):
        """True if page contains FluoView MM_STAMP tag."""
        return 'mm_stamp' in self.tags

    @lazyattr
    def is_nih(self):
        """True if page contains NIH image header."""
        return 'nih_image_header' in self.tags

    @lazyattr
    def is_ome(self):
        """True if page contains OME-XML in image_description tag."""
        return ('image_description' in self.tags and self.tags[
            'image_description'].value.startswith(b'<?xml version='))

    @lazyattr
    def is_shaped(self):
        """True if page contains shape in image_description tag."""
        return ('image_description' in self.tags and self.tags[
            'image_description'].value.startswith(b'shape=('))


class TIFFtag(object):
    """A TIFF tag structure.

    Attributes
    ----------

    name : string
        Attribute name of tag.

    code : int
        Decimal code of tag.

    dtype : str
        Datatype of tag data. One of TIFF_DATA_TYPES.

    count : int
        Number of values.

    value : various types
        Tag data. For codes in CUSTOM_TAGS the 4 bytes file content.

    All attributes are read-only.

    """
    __slots__ = ('code', 'name', 'count', 'dtype', 'value', '_offset')

    def __init__(self, arg, **kwargs):
        """Initialize instance from file or arguments."""
        self._offset = None
        if hasattr(arg, 'fhandle'):
            self._fromfile(arg, **kwargs)
        else:
            self._fromdata(arg, **kwargs)

    def _fromdata(self, code, dtype, count, value, name=None):
        """Initialize instance from arguments."""
        self.code = int(code)
        self.name = name if name else str(code)
        self.dtype = TIFF_DATA_TYPES[dtype]
        self.count = int(count)
        self.value = value

    def _fromfile(self, parent):
        """Read tag structure from open file. Advances file cursor."""
        fhandle = parent.fhandle
        byte_order = parent.byte_order

        self._offset = fhandle.tell()

        fmt, size = {4: ('HHI4s', 12), 8: ('HHQ8s', 20)}[parent.offset_size]
        data = fhandle.read(size)
        code, dtype = struct.unpack(byte_order + fmt[:2], data[:4])
        count, value = struct.unpack(byte_order + fmt[2:], data[4:])

        if code in TIFF_TAGS:
            name = TIFF_TAGS[code][0]
        elif code in CUSTOM_TAGS:
            name = CUSTOM_TAGS[code][0]
        else:
            name = str(code)

        try:
            dtype = TIFF_DATA_TYPES[dtype]
        except KeyError:
            raise ValueError("unknown TIFF tag data type %i" % dtype)

        if not code in CUSTOM_TAGS:
            fmt = '%s%i%s' % (byte_order, count*int(dtype[0]), dtype[1])
            size = struct.calcsize(fmt)
            if size <= parent.offset_size:
                value = struct.unpack(fmt, value[:size])
            else:
                pos = fhandle.tell()
                tof = {4: 'I', 8: 'Q'}[parent.offset_size]
                fhandle.seek(struct.unpack(byte_order+tof, value)[0])
                value = struct.unpack(fmt, fhandle.read(size))
                fhandle.seek(pos)
            if len(value) == 1:
                value = value[0]
            if dtype == '1s':
                value = stripnull(value)

        self.code = code
        self.name = name
        self.dtype = dtype
        self.count = count
        self.value = value

    def __str__(self):
        """Return string containing information about tag."""
        return ' '.join(str(getattr(self, s)) for s in self.__slots__)


class Record(dict):
    """Dictionary with attribute access.

    Can also be initialized with numpy.core.records.record.

    """
    __slots__ = ()

    def __init__(self, arg=None, **kwargs):
        if kwargs:
            arg = kwargs
        elif arg is None:
            arg = {}
        try:
            dict.__init__(self, arg)
        except TypeError:
            for i, name in enumerate(arg.dtype.names):
                v = arg[i]
                self[name] = v if v.dtype.char != 'S' else stripnull(v)

    def __getattr__(self, name):
        return self[name]

    def __setattr__(self, name, value):
        self.__setitem__(name, value)

    def __str__(self):
        """Pretty print Record."""
        s = []
        lists = []
        for k in sorted(self):
            if k.startswith('_'):
                continue
            v = self[k]
            if isinstance(v, (list, tuple)) and len(v):
                if isinstance(v[0], Record):
                    lists.append((k, v))
                    continue
                elif isinstance(v[0], TIFFpage):
                    v = [i.index for i in v]
            s.append(("* %s: %s" % (k, str(v))).split("\n",
                                                      1)[0][:PRINT_LINE_LEN])
        for k, v in lists:
            l = []
            for i, w in enumerate(v):
                l.append("* %s[%i]\n  %s" % (k, i,
                                             str(w).replace("\n", "\n  ")))
            s.append('\n'.join(l))
        return '\n'.join(s)


class TiffTags(Record):
    """Dictionary of TIFFtags with attribute access."""
    def __str__(self):
        """Return string with information about all tags."""
        s = []
        #sortbycode = lambda a, b: cmp(a.code, b.code)
        #for tag in sorted(self.values(), sortbycode):
        for tag in sorted(self.values(), key=lambda x: x.code):
            typecode = "%i%s" % (tag.count * int(tag.dtype[0]), tag.dtype[1])
            line = "* %i %s (%s) %s" % (tag.code, tag.name, typecode,
                                        str(tag.value).split('\n', 1)[0])
            s.append(line[:PRINT_LINE_LEN])
        return '\n'.join(s)


def read_nih_image_header(fhandle, byte_order, tag):
    """Read NIH_IMAGE_HEADER tag from file and return as dictionary."""
    fhandle.seek(12 + struct.unpack(byte_order+'I', tag.value)[0])
    return {'version': struct.unpack(byte_order+'H', fhandle.read(2))[0]}


def read_mm_header(fhandle, byte_order, tag):
    """Read MM_HEADER tag from file and return as numpy.rec.array."""
    fhandle.seek(struct.unpack(byte_order+'I', tag.value)[0])
    return numpy.rec.fromfile(fhandle, MM_HEADER, 1, byteorder=byte_order)[0]


def read_mm_stamp(fhandle, byte_order, tag):
    """Read MM_STAMP tag from file and return as numpy.array."""
    fhandle.seek(struct.unpack(byte_order+'I', tag.value)[0])
    return numpy.fromfile(fhandle, byte_order+'8f8', 1)[0]


def read_mm_uic1(fhandle, byte_order, tag):
    """Read MM_UIC1 tag from file and return as dictionary."""
    fhandle.seek(struct.unpack(byte_order+'I', tag.value)[0])
    t = fhandle.read(8*tag.count)
    t = struct.unpack('%s%iI' % (byte_order, 2*tag.count), t)
    return dict((MM_TAG_IDS[k], v) for k, v in zip(t[::2], t[1::2])
                if k in MM_TAG_IDS)


def read_mm_uic2(fhandle, byte_order, tag):
    """Read MM_UIC2 tag from file and return as dictionary."""
    result = {'number_planes': tag.count}
    fhandle.seek(struct.unpack(byte_order+'I', tag.value)[0])
    values = numpy.fromfile(fhandle, byte_order+'I', 6*tag.count)
    result['z_distance'] = values[0::6] // values[1::6]
    #result['date_created'] = tuple(values[2::6])
    #result['time_created'] = tuple(values[3::6])
    #result['date_modified'] = tuple(values[4::6])
    #result['time_modified'] = tuple(values[5::6])
    return result


def read_mm_uic3(fhandle, byte_order, tag):
    """Read MM_UIC3 tag from file and return as dictionary."""
    fhandle.seek(struct.unpack(byte_order+'I', tag.value)[0])
    t = numpy.fromfile(fhandle, '%sI' % byte_order, 2*tag.count)
    return {'wavelengths': t[0::2] // t[1::2]}


def read_mm_uic4(fhandle, byte_order, tag):
    """Read MM_UIC4 tag from file and return as dictionary."""
    fhandle.seek(struct.unpack(byte_order+'I', tag.value)[0])
    t = struct.unpack(byte_order + 'hI'*tag.count, fhandle.read(6*tag.count))
    return dict((MM_TAG_IDS[k], v) for k, v in zip(t[::2], t[1::2])
                if k in MM_TAG_IDS)


def read_cz_lsm_info(fhandle, byte_order, tag):
    """Read CS_LSM_INFO tag from file and return as numpy.rec.array."""
    fhandle.seek(struct.unpack(byte_order+'I', tag.value)[0])
    result = numpy.rec.fromfile(fhandle, CZ_LSM_INFO, 1,
                                byteorder=byte_order)[0]
    {50350412: '1.3', 67127628: '2.0'}[result.magic_number]  # validation
    return result


def read_cz_lsm_time_stamps(fhandle, byte_order):
    """Read LSM time stamps from file and return as list."""
    size, count = struct.unpack(byte_order+'II', fhandle.read(8))
    if size != (8 + 8 * count):
        raise ValueError("lsm_time_stamps block is too short")
    return struct.unpack(('%s%dd' % (byte_order, count)),
                         fhandle.read(8*count))


def read_cz_lsm_event_list(fhandle, byte_order):
    """Read LSM events from file and return as list of (time, type, text)."""
    count = struct.unpack(byte_order+'II', fhandle.read(8))[1]
    events = []
    while count > 0:
        esize, etime, etype = struct.unpack(byte_order+'IdI', fhandle.read(16))
        etext = stripnull(fhandle.read(esize - 16))
        events.append((etime, etype, etext))
        count -= 1
    return events


def read_cz_lsm_scan_info(fhandle, byte_order):
    """Read LSM scan information from file and return as Record."""
    block = Record()
    blocks = [block]
    unpack = struct.unpack

    if 0x10000000 != struct.unpack(byte_order+"I", fhandle.read(4))[0]:
        raise ValueError("not a lsm_scan_info structure")
    fhandle.read(8)

    while True:
        entry, dtype, size = unpack(byte_order+"III", fhandle.read(12))
        if dtype == 2:
            value = stripnull(fhandle.read(size))
        elif dtype == 4:
            value = unpack(byte_order+"i", fhandle.read(4))[0]
        elif dtype == 5:
            value = unpack(byte_order+"d", fhandle.read(8))[0]
        else:
            value = 0

        if entry in CZ_LSM_SCAN_INFO_ARRAYS:
            blocks.append(block)
            name = CZ_LSM_SCAN_INFO_ARRAYS[entry]
            newobj = []
            setattr(block, name, newobj)
            block = newobj
        elif entry in CZ_LSM_SCAN_INFO_STRUCTS:
            blocks.append(block)
            newobj = Record()
            block.append(newobj)
            block = newobj
        elif entry in CZ_LSM_SCAN_INFO_ATTRIBUTES:
            name = CZ_LSM_SCAN_INFO_ATTRIBUTES[entry]
            setattr(block, name, value)
        elif entry == 0xffffffff:
            block = blocks.pop()
        else:
            setattr(block, "unknown_%x" % entry, value)

        if not blocks:
            break
    return block


def _replace_by(module_function, warn=True):
    """Try replace decorated function by module.function."""

    def decorate(func, module_function=module_function, warn=warn):
        sys.path.append(os.path.dirname(__file__))
        try:
            module, function = module_function.split('.')
            func, oldfunc = getattr(__import__(module), function), func
            globals()['__old_' + func.__name__] = oldfunc
        except Exception:
            if warn:
                warnings.warn("failed to import %s" % module_function)
        sys.path.pop()
        return func

    return decorate


@_replace_by('_tifffile.decodepackbits')
def decodepackbits(encoded):
    """Decompress PackBits encoded byte string.

    PackBits is a simple byte-oriented run-length compression scheme.

    """
    func = ord if sys.version[0] == '2' else lambda x: x
    result = []
    i = 0
    try:
        while True:
            n = func(encoded[i]) + 1
            i += 1
            if n < 129:
                result.extend(encoded[i:i+n])
                i += n
            elif n > 129:
                result.extend(encoded[i:i+1] * (258-n))
                i += 1
    except IndexError:
        pass
    return b''.join(result) if sys.version[0] == '2' else bytes(result)


@_replace_by('_tifffile.decodelzw')
def decodelzw(encoded):
    """Decompress LZW (Lempel-Ziv-Welch) encoded TIFF strip (byte string).

    The strip must begin with a CLEAR code and end with an EOI code.

    This is an implementation of the LZW decoding algorithm described in (1).
    It is not compatible with old style LZW compressed files like quad-lzw.tif.

    """
    unpack = struct.unpack

    if sys.version[0] == '2':
        newtable = [chr(i) for i in range(256)]
    else:
        newtable = [bytes([i]) for i in range(256)]
    newtable.extend((0, 0))

    def next_code():
        """Return integer of `bitw` bits at `bitcount` position in encoded."""
        start = bitcount // 8
        s = encoded[start:start+4]
        try:
            code = unpack('>I', s)[0]
        except Exception:
            code = unpack('>I', s + b'\x00'*(4-len(s)))[0]
        code = code << (bitcount % 8)
        code = code & mask
        return code >> shr

    switchbitch = {  # code: bit-width, shr-bits, bit-mask
        255: (9, 23, int(9*'1'+'0'*23, 2)),
        511: (10, 22, int(10*'1'+'0'*22, 2)),
        1023: (11, 21, int(11*'1'+'0'*21, 2)),
        2047: (12, 20, int(12*'1'+'0'*20, 2)), }
    bitw, shr, mask = switchbitch[255]
    bitcount = 0

    if len(encoded) < 4:
        raise ValueError("strip must be at least 4 characters long")

    if next_code() != 256:
        raise ValueError("strip must begin with CLEAR code")

    code = oldcode = 0
    result = []
    while True:
        code = next_code()  # ~5% faster when inlining this function
        bitcount += bitw
        if code == 257:  # EOI
            break
        if code == 256:  # CLEAR
            table = newtable[:]
            lentable = 258
            bitw, shr, mask = switchbitch[255]
            code = next_code()
            bitcount += bitw
            if code == 257:  # EOI
                break
            result.append(table[code])
        else:
            if code < lentable:
                decoded = table[code]
                newcode = table[oldcode] + decoded[:1]
            else:
                newcode = table[oldcode]
                newcode += newcode[:1]
                decoded = newcode
            result.append(decoded)
            table.append(newcode)
            lentable += 1
        oldcode = code
        if lentable in switchbitch:
            bitw, shr, mask = switchbitch[lentable]

    if code != 257:
        raise ValueError("unexpected end of stream (code %i)" % code)

    return b''.join(result)


@_replace_by('_tifffile.unpackints')
def unpackints(data, dtype, itemsize, runlen=0):
    """Decompress byte string to array of integers of any bit size <= 32.

    data : byte str

    dtype : numpy.dtype or str
        A numpy boolean or integer type.

    itemsize : int
        Number of bits per integer.

    runlen : int
        Number of consecutive integers, after which to start at next byte.

    """
    if itemsize == 1:  # bitarray
        data = numpy.fromstring(data, '|B')
        data = numpy.unpackbits(data)
        if runlen % 8:
            data = data.reshape(-1, runlen+(8-runlen%8))
            data = data[:, :runlen].reshape(-1)
        return data.astype(dtype)

    dtype = numpy.dtype(dtype)
    if itemsize in (8, 16, 32, 64):
        return numpy.fromstring(data, dtype)
    if itemsize < 1 or itemsize > 32:
        raise ValueError("itemsize out of range: %i" % itemsize)
    if dtype.kind not in "biu":
        raise ValueError("invalid dtype")

    itembytes = next(i for i in (1, 2, 4, 8) if 8 * i >= itemsize)
    if itembytes != dtype.itemsize:
        raise ValueError("dtype.itemsize too small")
    if runlen == 0:
        runlen = len(data) // itembytes
    skipbits = runlen*itemsize % 8
    if skipbits:
        skipbits = 8 - skipbits
    shrbits = itembytes*8 - itemsize
    bitmask = int(itemsize*'1'+'0'*shrbits, 2)
    dtypestr = '>' + dtype.char  # dtype always big endian?

    unpack = struct.unpack
    l = runlen * (len(data)*8 // (runlen*itemsize + skipbits))
    result = numpy.empty((l, ), dtype)
    bitcount = 0
    for i in range(len(result)):
        start = bitcount // 8
        s = data[start:start+itembytes]
        try:
            code = unpack(dtypestr, s)[0]
        except Exception:
            code = unpack(dtypestr, s + b'\x00'*(itembytes-len(s)))[0]
        code = code << (bitcount % 8)
        code = code & bitmask
        result[i] = code >> shrbits
        bitcount += itemsize
        if (i+1) % runlen == 0:
            bitcount += skipbits
    return result


def unpackrgb(data, dtype, bitspersample):
    """Unpack RGB565 from byte string."""
    if not (len(bitspersample) == 3
            and all(i <= 8 for i in bitspersample)
            and numpy.sum(bitspersample) == 16):
        raise ValueError("sample size not supported %s" % str(bitspersample))
    data = numpy.fromstring(data, dtype[0]+'H')
    r, g, b = bitspersample
    result = numpy.zeros((data.size, 3), 'uint8')
    result[:, 0] = ((data >> (r+g)) & (int('0b'+'1'*r, 2))) << (8 - r)
    result[:, 1] = ((data >> r) & (int('0b'+'1'*g, 2))) << (8 - g)
    result[:, 2] = (data & (int('0b'+'1'*b, 2))) << (8 - b)
    return result.reshape(-1)


def reorient(image, orientation):
    """Return reoriented view of image array.

    image : numpy array
        Non-squeezed output of asarray() functions.
        Axes -3 and -2 must be image length and width respectively.

    orientation : int or str
        One of TIFF_ORIENTATIONS keys or values

    """
    o = TIFF_ORIENTATIONS.get(orientation, orientation)
    if o == 'top_left':
        return image
    elif o == 'top_right':
        return image[..., ::-1, :]
    elif o == 'bottom_left':
        return image[..., ::-1, :, :]
    elif o == 'bottom_right':
        return image[..., ::-1, ::-1, :]
    elif o == 'left_top':
        return numpy.swapaxes(image, -3, -2)
    elif o == 'right_top':
        return numpy.swapaxes(image, -3, -2)[..., ::-1, :]
    elif o == 'left_bottom':
        return numpy.swapaxes(image, -3, -2)[..., ::-1, :, :]
    elif o == 'right_bottom':
        return numpy.swapaxes(image, -3, -2)[..., ::-1, ::-1, :]


def stripnull(string):
    """Return string truncated at first null character."""
    i = string.find(b'\x00')
    return string if (i < 0) else string[:i]


def datetime_from_timestamp(n, epoch=datetime.datetime.fromordinal(693594)):
    """Return datetime object from timestamp in Excel serial format.

    >>> datetime_from_timestamp(40237.029999999795)
    datetime.datetime(2010, 2, 28, 0, 43, 11, 999982)

    """
    return epoch + datetime.timedelta(n)


def test_tifffile(directory='testimages', verbose=True):
    """Read all images in directory. Print error message on failure.

    >>> test_tifffile(verbose=False)

    """
    import glob

    successful = 0
    failed = 0
    start = time.time()
    for f in glob.glob(os.path.join(directory, '*.*')):
        if verbose:
            print("\n%s>\n" % f.lower(), end='')
        t0 = time.time()
        try:
            tif = TIFFfile(f)
        except Exception as e:
            if not verbose:
                print(f, end=' ')
            print("ERROR:", e)
            failed += 1
            continue
        try:
            img = tif.asarray()
        except ValueError:
            try:
                img = tif[0].asarray()
            except Exception as e:
                if not verbose:
                    print(f, end=' ')
                print("ERROR:", e)
        finally:
            tif.close()
        successful += 1
        if verbose:
            print("%s, %s %s, %s, %.0f ms" % (str(tif), str(img.shape),
                img.dtype, tif[0].compression, (time.time()-t0) * 1e3))

    if verbose:
        print("\nSuccessfully read %i of %i files in %.3f s\n" % (
            successful, successful+failed, time.time()-start))


class TIFF_SUBFILE_TYPES(object):
    def __getitem__(self, key):
        result = []
        if key & 1:
            result.append('reduced_image')
        if key & 2:
            result.append('page')
        if key & 4:
            result.append('mask')
        return tuple(result)


TIFF_PHOTOMETRICS = {
    0: 'miniswhite',
    1: 'minisblack',
    2: 'rgb',
    3: 'palette',
    4: 'mask',
    5: 'separated',
    6: 'cielab',
    7: 'icclab',
    8: 'itulab',
    32844: 'logl',
    32845: 'logluv'}

TIFF_COMPESSIONS = {
    1: None,
    2: 'ccittrle',
    3: 'ccittfax3',
    4: 'ccittfax4',
    5: 'lzw',
    6: 'ojpeg',
    7: 'jpeg',
    8: 'adobe_deflate',
    9: 't85',
    10: 't43',
    32766: 'next',
    32771: 'ccittrlew',
    32773: 'packbits',
    32809: 'thunderscan',
    32895: 'it8ctpad',
    32896: 'it8lw',
    32897: 'it8mp',
    32898: 'it8bl',
    32908: 'pixarfilm',
    32909: 'pixarlog',
    32946: 'deflate',
    32947: 'dcs',
    34661: 'jbig',
    34676: 'sgilog',
    34677: 'sgilog24',
    34712: 'jp2000'}

TIFF_DECOMPESSORS = {
    None: lambda x: x,
    'adobe_deflate': zlib.decompress,
    'deflate': zlib.decompress,
    'packbits': decodepackbits,
    'lzw': decodelzw}

TIFF_DATA_TYPES = {
    1: '1B',   # BYTE 8-bit unsigned integer.
    2: '1s',   # ASCII 8-bit byte that contains a 7-bit ASCII code;
               #   the last byte must be NULL (binary zero).
    3: '1H',   # SHORT 16-bit (2-byte) unsigned integer
    4: '1I',   # LONG 32-bit (4-byte) unsigned integer.
    5: '2I',   # RATIONAL Two LONGs: the first represents the numerator of
               #   a fraction; the second, the denominator.
    6: '1b',   # SBYTE An 8-bit signed (twos-complement) integer.
    7: '1B',   # UNDEFINED An 8-bit byte that may contain anything,
               #   depending on the definition of the field.
    8: '1h',   # SSHORT A 16-bit (2-byte) signed (twos-complement) integer.
    9: '1i',   # SLONG A 32-bit (4-byte) signed (twos-complement) integer.
    10: '2i',  # SRATIONAL Two SLONGs: the first represents the numerator
               #   of a fraction, the second the denominator.
    11: '1f',  # FLOAT Single precision (4-byte) IEEE format.
    12: '1d',  # DOUBLE Double precision (8-byte) IEEE format.
    13: '1I',  # IFD unsigned 4 byte IFD offset.
    #14: '',   # UNICODE
    #15: '',   # COMPLEX
    16: '1Q',  # LONG8 unsigned 8 byte integer (BigTiff)
    17: '1q',  # SLONG8 signed 8 byte integer (BigTiff)
    18: '1Q'}  # IFD8 unsigned 8 byte IFD offset (BigTiff)

TIFF_SAMPLE_FORMATS = {
    1: 'uint',
    2: 'int',
    3: 'float',
    #4: 'void',
    #5: 'complex_int',
    6: 'complex'}

TIFF_SAMPLE_DTYPES = {
    ('uint', 1): '?',  # bitmap
    ('uint', 2): 'B',
    ('uint', 3): 'B',
    ('uint', 4): 'B',
    ('uint', 5): 'B',
    ('uint', 6): 'B',
    ('uint', 7): 'B',
    ('uint', 8): 'B',
    ('uint', 9): 'H',
    ('uint', 10): 'H',
    ('uint', 11): 'H',
    ('uint', 12): 'H',
    ('uint', 13): 'H',
    ('uint', 14): 'H',
    ('uint', 15): 'H',
    ('uint', 16): 'H',
    ('uint', 17): 'I',
    ('uint', 18): 'I',
    ('uint', 19): 'I',
    ('uint', 20): 'I',
    ('uint', 21): 'I',
    ('uint', 22): 'I',
    ('uint', 23): 'I',
    ('uint', 24): 'I',
    ('uint', 25): 'I',
    ('uint', 26): 'I',
    ('uint', 27): 'I',
    ('uint', 28): 'I',
    ('uint', 29): 'I',
    ('uint', 30): 'I',
    ('uint', 31): 'I',
    ('uint', 32): 'I',
    ('uint', 64): 'Q',
    ('int', 8): 'b',
    ('int', 16): 'h',
    ('int', 32): 'i',
    ('int', 64): 'q',
    ('float', 16): 'e',
    ('float', 32): 'f',
    ('float', 64): 'd',
    ('complex', 64): 'F',
    ('complex', 128): 'D',
    ('uint', (5, 6, 5)): 'B'}

TIFF_ORIENTATIONS = {
    1: 'top_left',
    2: 'top_right',
    3: 'bottom_right',
    4: 'bottom_left',
    5: 'left_top',
    6: 'right_top',
    7: 'right_bottom',
    8: 'left_bottom'}

# MetaMorph STK tags
MM_TAG_IDS = {
    0: 'auto_scale',
    1: 'min_scale',
    2: 'max_scale',
    3: 'spatial_calibration',
    #4: 'x_calibration',
    #5: 'y_calibration',
    #6: 'calibration_units',
    #7: 'name',
    8: 'thresh_state',
    9: 'thresh_state_red',
    11: 'thresh_state_green',
    12: 'thresh_state_blue',
    13: 'thresh_state_lo',
    14: 'thresh_state_hi',
    15: 'zoom',
    #16: 'create_time',
    #17: 'last_saved_time',
    18: 'current_buffer',
    19: 'gray_fit',
    20: 'gray_point_count',
    #21: 'gray_x',
    #22: 'gray_y',
    #23: 'gray_min',
    #24: 'gray_max',
    #25: 'gray_unit_name',
    26: 'standard_lut',
    27: 'wavelength',
    #28: 'stage_position',
    #29: 'camera_chip_offset',
    #30: 'overlay_mask',
    #31: 'overlay_compress',
    #32: 'overlay',
    #33: 'special_overlay_mask',
    #34: 'special_overlay_compress',
    #35: 'special_overlay',
    36: 'image_property',
    #37: 'stage_label',
    #38: 'autoscale_lo_info',
    #39: 'autoscale_hi_info',
    #40: 'absolute_z',
    #41: 'absolute_z_valid',
    #42: 'gamma',
    #43: 'gamma_red',
    #44: 'gamma_green',
    #45: 'gamma_blue',
    #46: 'camera_bin',
    47: 'new_lut',
    #48: 'image_property_ex',
    49: 'plane_property',
    #50: 'user_lut_table',
    51: 'red_autoscale_info',
    #52: 'red_autoscale_lo_info',
    #53: 'red_autoscale_hi_info',
    54: 'red_minscale_info',
    55: 'red_maxscale_info',
    56: 'green_autoscale_info',
    #57: 'green_autoscale_lo_info',
    #58: 'green_autoscale_hi_info',
    59: 'green_minscale_info',
    60: 'green_maxscale_info',
    61: 'blue_autoscale_info',
    #62: 'blue_autoscale_lo_info',
    #63: 'blue_autoscale_hi_info',
    64: 'blue_min_scale_info',
    65: 'blue_max_scale_info'}
    #66: 'overlay_plane_color',

# Olymus Fluoview
MM_DIMENSION = [
    ('name', 'a16'),
    ('size', 'i4'),
    ('origin', 'f8'),
    ('resolution', 'f8'),
    ('unit', 'a64')]

MM_HEADER = [
    ('header_flag', 'i2'),
    ('image_type', 'u1'),
    ('image_name', 'a257'),
    ('offset_data', 'u4'),
    ('palette_size', 'i4'),
    ('offset_palette0', 'u4'),
    ('offset_palette1', 'u4'),
    ('comment_size', 'i4'),
    ('offset_comment', 'u4'),
    ('dimensions', MM_DIMENSION, 10),
    ('offset_position', 'u4'),
    ('map_type', 'i2'),
    ('map_min', 'f8'),
    ('map_max', 'f8'),
    ('min_value', 'f8'),
    ('max_value', 'f8'),
    ('offset_map', 'u4'),
    ('gamma', 'f8'),
    ('offset', 'f8'),
    ('gray_channel', MM_DIMENSION),
    ('offset_thumbnail', 'u4'),
    ('voice_field', 'i4'),
    ('offset_voice_field', 'u4')]

# Carl Zeiss LSM
CZ_LSM_INFO = [
    ('magic_number', 'i4'),
    ('structure_size', 'i4'),
    ('dimension_x', 'i4'),
    ('dimension_y', 'i4'),
    ('dimension_z', 'i4'),
    ('dimension_channels', 'i4'),
    ('dimension_time', 'i4'),
    ('dimension_data_type', 'i4'),
    ('thumbnail_x', 'i4'),
    ('thumbnail_y', 'i4'),
    ('voxel_size_x', 'f8'),
    ('voxel_size_y', 'f8'),
    ('voxel_size_z', 'f8'),
    ('origin_x', 'f8'),
    ('origin_y', 'f8'),
    ('origin_z', 'f8'),
    ('scan_type', 'u2'),
    ('spectral_scan', 'u2'),
    ('data_type', 'u4'),
    ('offset_vector_overlay', 'u4'),
    ('offset_input_lut', 'u4'),
    ('offset_output_lut', 'u4'),
    ('offset_channel_colors', 'u4'),
    ('time_interval', 'f8'),
    ('offset_channel_data_types', 'u4'),
    ('offset_scan_information', 'u4'),
    ('offset_ks_data', 'u4'),
    ('offset_time_stamps', 'u4'),
    ('offset_event_list', 'u4'),
    ('offset_roi', 'u4'),
    ('offset_bleach_roi', 'u4'),
    ('offset_next_recording', 'u4'),
    ('display_aspect_x', 'f8'),
    ('display_aspect_y', 'f8'),
    ('display_aspect_z', 'f8'),
    ('display_aspect_time', 'f8'),
    ('offset_mean_of_roi_overlay', 'u4'),
    ('offset_topo_isoline_overlay', 'u4'),
    ('offset_topo_profile_overlay', 'u4'),
    ('offset_linescan_overlay', 'u4'),
    ('offset_toolbar_flags', 'u4')]

# Import functions for LSM_INFO subrecords
CZ_LSM_INFO_READERS = {
    'scan_information': read_cz_lsm_scan_info,
    'time_stamps': read_cz_lsm_time_stamps,
    'event_list': read_cz_lsm_event_list}

# Map cz_lsm_info.scan_type to dimension order
CZ_SCAN_TYPES = {
    0: 'XYZCT',  # x-y-z scan
    1: 'XYZCT',  # z scan (x-z plane)
    2: 'XYZCT',  # line scan
    3: 'XYTCZ',  # time series x-y
    4: 'XYZTC',  # time series x-z
    5: 'XYTCZ',  # time series 'Mean of ROIs'
    6: 'XYZTC',  # time series x-y-z
    7: 'XYCTZ',  # spline scan
    8: 'XYCZT',  # spline scan x-z
    9: 'XYTCZ',  # time series spline plane x-z
    10: 'XYZCT'}  # point mode

# Map dimension codes to cz_lsm_info attribute
CZ_DIMENSIONS = {
    'X': 'dimension_x',
    'Y': 'dimension_y',
    'Z': 'dimension_z',
    'C': 'dimension_channels',
    'T': 'dimension_time'}

# Descriptions of cz_lsm_info.data_type
CZ_DATA_TYPES = {
    0: 'varying data types',
    2: '12 bit unsigned integer',
    5: '32 bit float'}

CZ_LSM_SCAN_INFO_ARRAYS = {
    0x20000000: "tracks",
    0x30000000: "lasers",
    0x60000000: "detectionchannels",
    0x80000000: "illuminationchannels",
    0xa0000000: "beamsplitters",
    0xc0000000: "datachannels",
    0x13000000: "markers",
    0x11000000: "timers"}

CZ_LSM_SCAN_INFO_STRUCTS = {
    0x40000000: "tracks",
    0x50000000: "lasers",
    0x70000000: "detectionchannels",
    0x90000000: "illuminationchannels",
    0xb0000000: "beamsplitters",
    0xd0000000: "datachannels",
    0x14000000: "markers",
    0x12000000: "timers"}

CZ_LSM_SCAN_INFO_ATTRIBUTES = {
    0x10000001: "name",
    0x10000002: "description",
    0x10000003: "notes",
    0x10000004: "objective",
    0x10000005: "processing_summary",
    0x10000006: "special_scan_mode",
    0x10000007: "oledb_recording_scan_type",
    0x10000008: "oledb_recording_scan_mode",
    0x10000009: "number_of_stacks",
    0x1000000a: "lines_per_plane",
    0x1000000b: "samples_per_line",
    0x1000000c: "planes_per_volume",
    0x1000000d: "images_width",
    0x1000000e: "images_height",
    0x1000000f: "images_number_planes",
    0x10000010: "images_number_stacks",
    0x10000011: "images_number_channels",
    0x10000012: "linscan_xy_size",
    0x10000013: "scan_direction",
    0x10000014: "time_series",
    0x10000015: "original_scan_data",
    0x10000016: "zoom_x",
    0x10000017: "zoom_y",
    0x10000018: "zoom_z",
    0x10000019: "sample_0x",
    0x1000001a: "sample_0y",
    0x1000001b: "sample_0z",
    0x1000001c: "sample_spacing",
    0x1000001d: "line_spacing",
    0x1000001e: "plane_spacing",
    0x1000001f: "plane_width",
    0x10000020: "plane_height",
    0x10000021: "volume_depth",
    0x10000023: "nutation",
    0x10000034: "rotation",
    0x10000035: "precession",
    0x10000036: "sample_0time",
    0x10000037: "start_scan_trigger_in",
    0x10000038: "start_scan_trigger_out",
    0x10000039: "start_scan_event",
    0x10000040: "start_scan_time",
    0x10000041: "stop_scan_trigger_in",
    0x10000042: "stop_scan_trigger_out",
    0x10000043: "stop_scan_event",
    0x10000044: "stop_scan_time",
    0x10000045: "use_rois",
    0x10000046: "use_reduced_memory_rois",
    0x10000047: "user",
    0x10000048: "use_bccorrection",
    0x10000049: "position_bccorrection1",
    0x10000050: "position_bccorrection2",
    0x10000051: "interpolation_y",
    0x10000052: "camera_binning",
    0x10000053: "camera_supersampling",
    0x10000054: "camera_frame_width",
    0x10000055: "camera_frame_height",
    0x10000056: "camera_offset_x",
    0x10000057: "camera_offset_y",
    # lasers
    0x50000001: "name",
    0x50000002: "acquire",
    0x50000003: "power",
    # tracks
    0x40000001: "multiplex_type",
    0x40000002: "multiplex_order",
    0x40000003: "sampling_mode",
    0x40000004: "sampling_method",
    0x40000005: "sampling_number",
    0x40000006: "acquire",
    0x40000007: "sample_observation_time",
    0x4000000b: "time_between_stacks",
    0x4000000c: "name",
    0x4000000d: "collimator1_name",
    0x4000000e: "collimator1_position",
    0x4000000f: "collimator2_name",
    0x40000010: "collimator2_position",
    0x40000011: "is_bleach_track",
    0x40000012: "is_bleach_after_scan_number",
    0x40000013: "bleach_scan_number",
    0x40000014: "trigger_in",
    0x40000015: "trigger_out",
    0x40000016: "is_ratio_track",
    0x40000017: "bleach_count",
    0x40000018: "spi_center_wavelength",
    0x40000019: "pixel_time",
    0x40000021: "condensor_frontlens",
    0x40000023: "field_stop_value",
    0x40000024: "id_condensor_aperture",
    0x40000025: "condensor_aperture",
    0x40000026: "id_condensor_revolver",
    0x40000027: "condensor_filter",
    0x40000028: "id_transmission_filter1",
    0x40000029: "id_transmission1",
    0x40000030: "id_transmission_filter2",
    0x40000031: "id_transmission2",
    0x40000032: "repeat_bleach",
    0x40000033: "enable_spot_bleach_pos",
    0x40000034: "spot_bleach_posx",
    0x40000035: "spot_bleach_posy",
    0x40000036: "spot_bleach_posz",
    0x40000037: "id_tubelens",
    0x40000038: "id_tubelens_position",
    0x40000039: "transmitted_light",
    0x4000003a: "reflected_light",
    0x4000003b: "simultan_grab_and_bleach",
    0x4000003c: "bleach_pixel_time",
    # detection_channels
    0x70000001: "integration_mode",
    0x70000002: "special_mode",
    0x70000003: "detector_gain_first",
    0x70000004: "detector_gain_last",
    0x70000005: "amplifier_gain_first",
    0x70000006: "amplifier_gain_last",
    0x70000007: "amplifier_offs_first",
    0x70000008: "amplifier_offs_last",
    0x70000009: "pinhole_diameter",
    0x7000000a: "counting_trigger",
    0x7000000b: "acquire",
    0x7000000c: "point_detector_name",
    0x7000000d: "amplifier_name",
    0x7000000e: "pinhole_name",
    0x7000000f: "filter_set_name",
    0x70000010: "filter_name",
    0x70000013: "integrator_name",
    0x70000014: "detection_channel_name",
    0x70000015: "detection_detector_gain_bc1",
    0x70000016: "detection_detector_gain_bc2",
    0x70000017: "detection_amplifier_gain_bc1",
    0x70000018: "detection_amplifier_gain_bc2",
    0x70000019: "detection_amplifier_offset_bc1",
    0x70000020: "detection_amplifier_offset_bc2",
    0x70000021: "detection_spectral_scan_channels",
    0x70000022: "detection_spi_wavelength_start",
    0x70000023: "detection_spi_wavelength_stop",
    0x70000026: "detection_dye_name",
    0x70000027: "detection_dye_folder",
    # illumination_channels
    0x90000001: "name",
    0x90000002: "power",
    0x90000003: "wavelength",
    0x90000004: "aquire",
    0x90000005: "detchannel_name",
    0x90000006: "power_bc1",
    0x90000007: "power_bc2",
    # beam_splitters
    0xb0000001: "filter_set",
    0xb0000002: "filter",
    0xb0000003: "name",
    # data_channels
    0xd0000001: "name",
    0xd0000003: "acquire",
    0xd0000004: "color",
    0xd0000005: "sample_type",
    0xd0000006: "bits_per_sample",
    0xd0000007: "ratio_type",
    0xd0000008: "ratio_track1",
    0xd0000009: "ratio_track2",
    0xd000000a: "ratio_channel1",
    0xd000000b: "ratio_channel2",
    0xd000000c: "ratio_const1",
    0xd000000d: "ratio_const2",
    0xd000000e: "ratio_const3",
    0xd000000f: "ratio_const4",
    0xd0000010: "ratio_const5",
    0xd0000011: "ratio_const6",
    0xd0000012: "ratio_first_images1",
    0xd0000013: "ratio_first_images2",
    0xd0000014: "dye_name",
    0xd0000015: "dye_folder",
    0xd0000016: "spectrum",
    0xd0000017: "acquire",
    # markers
    0x14000001: "name",
    0x14000002: "description",
    0x14000003: "trigger_in",
    0x14000004: "trigger_out",
    # timers
    0x12000001: "name",
    0x12000002: "description",
    0x12000003: "interval",
    0x12000004: "trigger_in",
    0x12000005: "trigger_out",
    0x12000006: "activation_time",
    0x12000007: "activation_number"}

# Map TIFF tag code to attribute name, default value, type, count, validator
TIFF_TAGS = {
    254: ('new_subfile_type', 0, 4, 1, TIFF_SUBFILE_TYPES()),
    255: ('subfile_type', None, 3, 1,
          {0: 'undefined', 1: 'image', 2: 'reduced_image', 3: 'page'}),
    256: ('image_width', None, 4, 1, None),
    257: ('image_length', None, 4, 1, None),
    258: ('bits_per_sample', 1, 3, 1, None),
    259: ('compression', 1, 3, 1, TIFF_COMPESSIONS),
    262: ('photometric', None, 3, 1, TIFF_PHOTOMETRICS),
    266: ('fill_order', 1, 3, 1, {1: 'msb2lsb', 2: 'lsb2msb'}),
    269: ('document_name', None, 2, None, None),
    270: ('image_description', None, 2, None, None),
    271: ('make', None, 2, None, None),
    272: ('model', None, 2, None, None),
    273: ('strip_offsets', None, 4, None, None),
    274: ('orientation', 1, 3, 1, TIFF_ORIENTATIONS),
    277: ('samples_per_pixel', 1, 3, 1, None),
    278: ('rows_per_strip', 2**32-1, 4, 1, None),
    279: ('strip_byte_counts', None, 4, None, None),
    #280: ('min_sample_value', 0, 3, None, None),
    #281: ('max_sample_value', None, 3, None, None),  # 2**bits_per_sample
    282: ('x_resolution', None, 5, 1, None),
    283: ('y_resolution', None, 5, 1, None),
    284: ('planar_configuration', 1, 3, 1, {1: 'contig', 2: 'separate'}),
    285: ('page_name', None, 2, None, None),
    296: ('resolution_unit', 2, 4, 1, {1: 'none', 2: 'inch', 3: 'centimeter'}),
    305: ('software', None, 2, None, None),
    306: ('datetime', None, 2, None, None),
    315: ('artist', None, 2, None, None),
    316: ('host_computer', None, 2, None, None),
    317: ('predictor', 1, 3, 1, {1: None, 2: 'horizontal'}),
    320: ('color_map', None, 3, None, None),
    322: ('tile_width', None, 4, 1, None),
    323: ('tile_length', None, 4, 1, None),
    324: ('tile_offsets', None, 4, None, None),
    325: ('tile_byte_counts', None, 4, None, None),
    338: ('extra_samples', None, 3, None,
          {0: 'unspecified', 1: 'assocalpha', 2: 'unassalpha'}),
    339: ('sample_format', 1, 3, 1, TIFF_SAMPLE_FORMATS),
    530: ('ycbcr_subsampling', 1, 3, 2, None),
    531: ('ycbcr_positioning', 1, 3, 1, None),
    #700: ('xmp', None, 1, None, None),
    33432: ('copyright', None, 1, None, None),
    32997: ('image_depth', None, 4, 1, None),
    32998: ('tile_depth', None, 4, 1, None),
    34665: ('exif_ifd', None, 4, 1, None)}

TIFF_TAGS_W = dict((val[0], key) for (key, val) in TIFF_TAGS.items())

# Map custom TIFF tag codes to attribute names and import functions
CUSTOM_TAGS = {
    33628: ('mm_uic1', read_mm_uic1),
    33629: ('mm_uic2', read_mm_uic2),
    33630: ('mm_uic3', read_mm_uic3),
    33631: ('mm_uic4', read_mm_uic4),
    34361: ('mm_header', read_mm_header),
    34362: ('mm_stamp', read_mm_stamp),
    34386: ('mm_user_block', None),
    34412: ('cz_lsm_info', read_cz_lsm_info),
    43314: ('nih_image_header', read_nih_image_header)}

# Max line length of printed output
PRINT_LINE_LEN = 79


def imshow(data, title=None, isrgb=True, vmin=0, vmax=None,
           cmap=None, photometric='rgb', interpolation='bilinear',
           dpi=96, figure=None, subplot=111, maxdim=4096, **kwargs):
    """Plot n-dimensional images using matplotlib.pyplot.

    Return figure, subplot and plot axis.
    Requires pyplot already imported ``from matplotlib import pyplot``.

    Arguments
    ---------

    isrgb : bool
        If True data will be displayed as color images if possible.

    photometric : str
        'miniswhite', 'minisblack', 'rgb', or 'palette'

    title : str
        Window and subplot title.

    figure : matplotlib.figure.Figure instance (optional).

    subplot : int
        A matplotlib.pyplot.subplot axis.

    maxdim : int
        maximum image size in any dimension.

    kwargs : optional arguments for matplotlib.pyplot.imshow.

    """
    #if photometric not in ('miniswhite', 'minisblack', 'rgb', 'palette'):
    #    raise ValueError("Can't handle %s photometrics" % photometric)
    data = numpy.atleast_2d(data.squeeze())
    data = data[(slice(0, maxdim), ) * len(data.shape)]

    dims = len(data.shape)
    if dims < 2:
        raise ValueError("not an image")
    if dims == 2:
        dims = 0
        isrgb = False
    else:
        if (isrgb and data.shape[-3] in (3, 4)):
            data = numpy.swapaxes(data, -3, -2)
            data = numpy.swapaxes(data, -2, -1)
        elif (not isrgb and data.shape[-1] in (3, 4)):
            data = numpy.swapaxes(data, -3, -1)
            data = numpy.swapaxes(data, -2, -1)
        isrgb = isrgb and data.shape[-1] in (3, 4)
        dims -= 3 if isrgb else 2

    datamax = data.max()
    if data.dtype in (numpy.int8, numpy.int16, numpy.int32,
                      numpy.uint8, numpy.uint16, numpy.uint32):
        datamax = 2**int(math.ceil(math.log(datamax, 2)))
        if isrgb:
            data *= (255.0 / datamax)  # better use digitize()
            data = data.astype('B')
    elif isrgb:
        data /= datamax

    if not isrgb and vmax is None:
        vmax = datamax

    pyplot = sys.modules['matplotlib.pyplot']

    if figure is None:
        pyplot.rc('font', family='sans-serif', weight='normal', size=8)
        figure = pyplot.figure(dpi=dpi, figsize=(10.3, 6.3), frameon=True,
                               facecolor='1.0', edgecolor='w')
        try:
            figure.canvas.manager.window.title(title)
        except Exception:
            pass
        pyplot.subplots_adjust(bottom=0.03*(dims+2), top=0.9,
                               left=0.1, right=0.95, hspace=0.05, wspace=0.0)
    subplot = pyplot.subplot(subplot)

    if title:
        pyplot.title(title, size=11)

    if cmap is None:
        if photometric == 'miniswhite':
            cmap = pyplot.cm.binary
        elif data.dtype in (numpy.uint8, numpy.bool):
            cmap = pyplot.cm.gray
        else:
            cmap = pyplot.cm.gray  # coolwarm

    image = pyplot.imshow(data[(0, ) * dims].squeeze(), vmin=vmin, vmax=vmax,
                          cmap=cmap, interpolation=interpolation, **kwargs)

    if not isrgb:
        pyplot.colorbar()  # panchor=(0.55, 0.5), fraction=0.05

    def format_coord(x, y):
        # callback function to format coordinate display in toolbar
        x = int(x + 0.5)
        y = int(y + 0.5)
        try:
            if dims:
                return "%s @ %s [%4i, %4i]" % (cur_ax_dat[1][y, x],
                                               current, x, y)
            else:
                return "%s @ [%4i, %4i]" % (data[y, x], x, y)
        except IndexError:
            return ""

    pyplot.gca().format_coord = format_coord

    if dims:
        current = list((0, ) * dims)
        cur_ax_dat = [0, data[tuple(current)].squeeze()]
        sliders = [pyplot.Slider(
            pyplot.axes([0.125, 0.03*(axis+1), 0.725, 0.025]),
            'Dimension %i' % axis, 0, data.shape[axis]-1, 0, facecolor='0.5',
            valfmt='%%.0f [%i]' % data.shape[axis]) for axis in range(dims)]
        for slider in sliders:
            slider.drawon = False

        def set_image(current, sliders=sliders, data=data):
            # change image and redraw canvas
            cur_ax_dat[1] = data[tuple(current)].squeeze()
            image.set_data(cur_ax_dat[1])
            for ctrl, index in zip(sliders, current):
                ctrl.eventson = False
                ctrl.set_val(index)
                ctrl.eventson = True
            figure.canvas.draw()

        def on_changed(index, axis, data=data, current=current):
            # callback function for slider change event
            index = int(round(index))
            cur_ax_dat[0] = axis
            if index == current[axis]:
                return
            if index >= data.shape[axis]:
                index = 0
            elif index < 0:
                index = data.shape[axis] - 1
            current[axis] = index
            set_image(current)

        def on_keypressed(event, data=data, current=current):
            # callback function for key press event
            key = event.key
            axis = cur_ax_dat[0]
            if str(key) in '0123456789':
                on_changed(key, axis)
            elif key == 'right':
                on_changed(current[axis] + 1, axis)
            elif key == 'left':
                on_changed(current[axis] - 1, axis)
            elif key == 'up':
                cur_ax_dat[0] = 0 if axis == len(data.shape)-1 else axis + 1
            elif key == 'down':
                cur_ax_dat[0] = len(data.shape)-1 if axis == 0 else axis - 1
            elif key == 'end':
                on_changed(data.shape[axis] - 1, axis)
            elif key == 'home':
                on_changed(0, axis)

        figure.canvas.mpl_connect('key_press_event', on_keypressed)
        for axis, ctrl in enumerate(sliders):
            ctrl.on_changed(lambda k, a=axis: on_changed(k, a))

    return figure, subplot, image


def main(argv=None):
    """Command line usage main function."""
    if float(sys.version[0:3]) < 2.6:
        print("This script requires Python version 2.6 or better.")
        print("This is Python version %s" % sys.version)
        return 0
    if argv is None:
        argv = sys.argv

    import re
    import optparse

    search_doc = lambda r, d: re.search(r, __doc__).group(1) if __doc__ else d
    parser = optparse.OptionParser(
        usage="usage: %prog [options] path",
        description=search_doc("\n\n([^|]*?)\n\n", ''),
        version="%%prog %s" % search_doc(":Version: (.*)", "Unknown"))
    opt = parser.add_option
    opt('-p', '--page', dest='page', type='int', default=-1,
        help="display single page")
    opt('-s', '--series', dest='series', type='int', default=-1,
        help="display series of pages of same shape")
    opt('--noplot', dest='noplot', action='store_true', default=False,
        help="don't display images")
    opt('--interpol', dest='interpol', metavar='INTERPOL', default='bilinear',
        help="image interpolation method")
    opt('--dpi', dest='dpi', type='int', default=96,
        help="set plot resolution")
    opt('--debug', dest='debug', action='store_true', default=False,
        help="raise exception on failures")
    opt('--test', dest='test', action='store_true', default=False,
        help="try read all images in path")
    opt('--doctest', dest='doctest', action='store_true', default=False,
        help="runs the internal tests")
    opt('-v', '--verbose', dest='verbose', action='store_true', default=True)
    opt('-q', '--quiet', dest='verbose', action='store_false')

    settings, path = parser.parse_args()
    path = ' '.join(path)

    if settings.doctest:
        import doctest
        doctest.testmod()
        return 0
    if not path:
        parser.error("No file specified")
    if settings.test:
        test_tifffile(path, settings.verbose)
        return 0

    print("Reading file structure...", end=' ')
    start = time.time()
    try:
        tif = TIFFfile(path)
    except Exception as e:
        if settings.debug:
            raise
        else:
            print("\n", e)
            sys.exit(0)
    print("%.3f ms" % ((time.time()-start) * 1e3))

    if tif.is_ome:
        settings.norgb = True

    images = [(None, tif[0 if settings.page < 0 else settings.page])]
    if not settings.noplot:
        print("Reading image data... ", end=' ')
        notnone = lambda x: next(i for i in x if i is not None)
        start = time.time()
        try:
            if settings.page >= 0:
                images = [(tif.asarray(key=settings.page),
                           tif[settings.page])]
            elif settings.series >= 0:
                images = [(tif.asarray(series=settings.series),
                           notnone(tif.series[settings.series].pages))]
            else:
                images = []
                for i, s in enumerate(tif.series):
                    try:
                        images.append(
                            (tif.asarray(series=i), notnone(s.pages)))
                    except ValueError as e:
                        images.append((None, notnone(s.pages)))
                        if settings.debug:
                            raise
                        else:
                            print("\n* series %i failed: %s... " % (i, e),
                                  end='')
            print("%.3f ms" % ((time.time()-start) * 1e3))
        except Exception as e:
            if settings.debug:
                raise
            else:
                print(e)
    tif.close()

    print("\nTIFF file:", tif)
    print()
    for i, s in enumerate(tif.series):
        print ("Series %i" % i)
        print(s)
        print()
    for i, page in images:
        print(page)
        print(page.tags)
        if page.is_palette:
            print("\nColor Map:", page.color_map.shape, page.color_map.dtype)
        for attr in ('cz_lsm_info', 'cz_lsm_scan_information',
                     'mm_uic_tags', 'mm_header', 'nih_image_header'):
            if hasattr(page, attr):
                print("", attr.upper(), Record(getattr(page, attr)), sep="\n")
        print()

    if images and not settings.noplot:
        try:
            import matplotlib
            matplotlib.use('TkAgg')
            from matplotlib import pyplot
        except ImportError as e:
            warnings.warn("failed to import matplotlib.\n%s" % e)
        else:
            for img, page in images:
                if img is None:
                    continue
                vmin, vmax = None, None
                if page.is_stk:
                    try:
                        vmin = page.mm_uic_tags['min_scale']
                        vmax = page.mm_uic_tags['max_scale']
                    except KeyError:
                        pass
                    else:
                        if vmax <= vmin:
                            vmin, vmax = None, None
                title = "%s\n %s" % (str(tif), str(page))
                imshow(img, title=title, vmin=vmin, vmax=vmax,
                       photometric=page.photometric,
                       interpolation=settings.interpol,
                       dpi=settings.dpi,
                       isrgb=page.tags['photometric'].value > 1)
            pyplot.show()


# Documentation in HTML format can be generated with Epydoc
__docformat__ = "restructuredtext en"

if __name__ == "__main__":
    sys.exit(main())