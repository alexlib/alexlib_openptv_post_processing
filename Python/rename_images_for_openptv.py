from glob import glob
import os

filelist = glob('dumbbell_end/*.TIF')

filelist

cam = lambda x: x.split('_')[3].split('-')[0]
cam(filelist[0])

frame = lambda x: '1'+x.split('_')[-1].split('.')[0]
frame(filelist[0])

newname = lambda x: 'db'+cam(x)+'.'+frame(x)
# newname = lambda x: 'cam'+cam(x)+'.'+frame(x)
print newname(filelist[0])


for f in filelist:
    os.rename(f,os.path.join(os.path.dirname(f),newname(f)))

