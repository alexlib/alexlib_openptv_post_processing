from tables import *
import numpy
import os
import ftplib
import sys

# ------------------------------------------------


class TrajPoint(IsDescription):
    traj_id = Int32Col(pos=0)
#    coordinates
    x = Float32Col(pos=1)
    y = Float32Col(pos=2)
    z = Float32Col(pos=3)
#  velocities
    u = Float32Col(pos=4)
    v = Float32Col(pos=5)
    w = Float32Col(pos=6)
# full (lagrangian) accelerations
    ax = Float32Col(pos=7)
    ay = Float32Col(pos=8)
    az = Float32Col(pos=9)
# spatial quantities           
# vorticity
    w1 = Float32Col(pos=10)
    w2 = Float32Col(pos=11)
    w3 = Float32Col(pos=12)
# strain rate 
    s11= Float32Col(pos=13)
    s12 = Float32Col(pos=14)
    s13 = Float32Col(pos=15)
    s22 = Float32Col(pos=16)
    s23 = Float32Col(pos=17)
    s33 = Float32Col(pos=18)
# local acceleration
    alx = Float32Col(pos=19)
    aly = Float32Col(pos=20)
    alz = Float32Col(pos=21)
#   convective acceleration
    daxdx = Float32Col(pos=22)
    daxdy = Float32Col(pos=23)
    daxdz = Float32Col(pos=24)
    daydx = Float32Col(pos=25)
    daydy = Float32Col(pos=26)
    daydz = Float32Col(pos=27)
    dazdx = Float32Col(pos=28)
    dazdy = Float32Col(pos=29)
    dazdz = Float32Col(pos=30)
# annotation: quality    
    quality = UInt8Col(pos=31)
# global time    
    t = Int32Col(pos=32)
    
# ------------------------------------------------
    

def gettext(ftp, filename, outfile=None):
    # fetch a text file
    if outfile is None:
        outfile = sys.stdout
    # use a lambda to add newlines to the lines read from the server
    ftp.retrlines("RETR " + filename, lambda s, w=outfile.write: w(s+"\n"))
    
# ------------------------------------------------
    

def getbinary(ftp, filename, outfile=None):
    # fetch a binary file
    if outfile is None:
        outfile = sys.stdout
    ftp.retrbinary("RETR " + filename, outfile.write)

# ------------------------------------------------

ftp = ftplib.FTP("ifu-gwh-disk.ethz.ch")
ftp.login("gwh","6WH.lab06")
ftp.cwd("/Hydromechanik/Rotating/PTV08/080827/run2")

data = []
ftp.dir(data.append)

h5file = openFile("run2_080827.h5", mode = "w", title = "27-Aug-2008 \omega 1/5Hz (9.9), instationary ")
table = h5file.createTable('/','trajPoint',TrajPoint,'TrajPoint example')
particle = table.row
traj_id = 0


for line in data:
    filename =  line[line.find("g_")::]
    if len(filename) < 11:
        continue
    basename, extension = filename.split('.')
    # print filename
    f = file('junkfile.txt','w') # f = file(filename,'w+')
    getbinary(ftp, filename,f)
    f.close()
    if os.path.getsize('junkfile.txt') == 0:
			continue
    # print tmp.shape
    for row in numpy.loadtxt('junkfile.txt'):
        if row[31] == 0:
            traj_id+=1
        particle['x'] = row[0]
        particle['y'] = row[1]
        particle['z'] = row[2]
        
        particle['u'] = row[3]
        particle['v'] = row[4]
        particle['w'] = row[5]
        
        particle['ax'] = row[6]
        particle['ay'] = row[7]
        particle['az'] = row[8]
        
        particle['w1'] = row[9]
        particle['w2'] = row[10]
        particle['w3'] = row[11]
        
        particle['s11'] = row[12]
        particle['s12'] = row[13]
        particle['s13'] = row[14]
        particle['s22'] = row[15]
        particle['s23'] = row[16]
        particle['s33'] = row[17]
        
        particle['alx'] = row[18]
        particle['aly'] = row[19]
        particle['alz'] = row[20]
        
        particle['daxdx'] = row[21]
        particle['daxdy'] = row[22]
        particle['daxdz'] = row[23]
        particle['daydx'] = row[24]
        particle['daydy'] = row[25]
        particle['daydz'] = row[26]
        particle['dazdx'] = row[27]
        particle['dazdy'] = row[28]
        particle['dazdz'] = row[29]
        particle['quality'] = row[30]
        particle['t'] = row[31]+eval(extension)
        particle['traj_id'] = traj_id
        particle.append()
    
    

table.flush()
h5file.close() 


    