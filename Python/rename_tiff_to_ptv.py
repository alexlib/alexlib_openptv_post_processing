filelist = os.listdir('.')
for i,f in enumerate(filelist):
    os.rename(f,'cam1.1%05d' % i)
