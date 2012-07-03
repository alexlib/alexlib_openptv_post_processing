from tables import *
import numpy
import os


path = "D:\\liberzon\Desktop\\rotating\\080825\\run6"  # insert the path to the directory of interest

h5file = openFile("run6_080825.h5", mode = "w", title = "Run 4, 25-Aug-08 \omega 7")


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
    
    


table = h5file.createTable('/','trajPoint',TrajPoint,'TrajPoint example')
particle = table.row
traj_id = 0

for filename in os.listdir(path):
    basename, extension = filename.split('.')
    if basename == 'g_trajPoint':
    # if i.split(".")[len(i.split("."))-1] == grep:
        # print ''.join([path,'\\',filename])
        if os.path.getsize(''.join([path,'\\',filename])) == 0:
			continue
        tmp = numpy.loadtxt(''.join([path,'\\',filename]))
        # print tmp.shape
        for row in tmp:
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


    