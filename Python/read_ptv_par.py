def read_ptv_par():

    f = open ('parameters/ptv.par','r');
	
    n_img = n.double(f.readline())
    
    img_name=[]
    img_cal =[]
    for i in n.arange(n_img):
        img_name.append(f.readline().strip())
        img_cal.append(f.readline().strip())


    hp_flag = n.bool(f.readline())
    allCam_flag = n.bool(f.readline())
    tiff_flag = n.bool(f.readline())
    imx = n.int(f.readline())
    imy = n.int(f.readline())
    pix_x = n.double(f.readline())
    pix_y = n.double(f.readline())
    chfield = n.double(f.readline())
    mmp_n1 = n.double(f.readline())
    mmp_n2 = n.double(f.readline())
    mmp_n3 = n.double(f.readline())
    mmp_d =  n.double(f.readline())
    
    f.close()
