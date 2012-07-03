import numpy as n

# ====================================================================

def _write_ptv_par():
        

	""" ptv.par
	ptv.par:	main parameter file
	4	number of cameras 
	cam3.100	image of first camera 
	kal1	calibration data of first camera 
	cam0.100	image of second camera 
	kal3	calibration data of second camera
	cam1.100	image of third camera 
	kal4	calibration data of third camera 
	cam2.100	image of fourth camera 
	kal5	calibration data of fourth camera 
	1	flag for highpass filtering, use (1) or not use (0) 
	1	flag for TIFF header (1) or raw data (0) 
	720	image width in pixel 
	576	image height in pixel 
	0.009	pixel size horizontal [mm] 
	0.0084	pixel size vertical [mm] 
	0	flag for frame, odd or even fields 
	1.0	refractive index air [no unit] 
	1.5	refractive index glass [no unit] 
	1.0	refractive index water [no unit] 
	9.4	thickness of glass [mm]
	
	"""
	f = open ('parameters/ptv.par','w')
	f.write("%d\n" % n_img)
	for i in n.arange(n_img):
		f.write("%s\n" % img_name[i])
		f.write("%s\n" % img_cal[i])
		
	f.write("%d\n" % hp_flag)
	f.write("%d\n" % allCam_flag)
	f.write("%d\n" % tiff_flag)
	f.write("%d\n" % imx)
	f.write("%d\n" % imy)
	f.write("%d\n" % pix_x)
	f.write("%d\n" % pix_y)
	f.write("%d\n" % chfield)
	f.write("%d\n" % mmp_n1)
	f.write("%d\n" % mmp_n2)
	f.write("%d\n" % mmp_n3)
	f.write("%d\n" % mmp_d)
	
	f.close()

# ====================================================================

def _write_cal_ori_par(n_img, fixp_name, img_cal_name, img_ori, tiff_flag, pair_flag, chfield):
	# calibration parameters
	"""
	cal_ori.par:	calibration plate, images, orientation files
	ptv/ssc_cal.c3d	control point file (point number, X, Y, Z in [mm], ASCII
	kal1	calibration 
	kal1.ori	orientation 
	kal3	calibration 
	kal3.ori	orientation 
	kal4	calibration 
	kal4.ori	orientation 
	kal5	calibration 
	kal5.ori	orientation 
	1	flag for TIFF header (1) or raw data (0)
	0   flag for pairs? 
	0	flag for frame (0), odd (1) or even fields (2)
	"""
	
	f = open ('parameters/cal_ori.par','w');
	f.write("%s\n" % fixp_name) 

	for i in n.arange(n_img):
		f.write("%s\n" % img_cal_name[i])
		f.write("%s\n" % img_ori[i])
	
	
	f.write("%d\n" % tiff_flag)
	f.write("%d\n" % pair_flag)
	f.write("%d\n" % chfield)
	f.close()

	
# ====================================================================

def _write_sequence_par():
	"""
	sequence.par: sequence parameters
	cam0. basename for 1.sequence
	cam1. basename for 2. sequence
	cam2. basename for 3. sequence
	cam3. basename for 4. sequence
	100  first image of sequence
	119  last image of sequence
		 
	"""
	f = open('parameters/sequence.par','w');
	
	for i in n.arange(n_img):
		f.write("%s\n" % base_name[i])
	
	
	f.write("%d\n" % first)
	f.write("%d\n" % last)
	f.close()

#============================================

def _write_criteria_par():
	"""
	criteria.par:	object volume and correspondence parameters
	0.0	illuminated layer data, xmin [mm]
	-10.0	illuminated layer data, zmin [mm]
	0.0	illuminated layer data, zmax [mm]
	10.0	illuminated layer data, xmax [mm]
	-10.0	illuminated layer data, zmin [mm]
	0.0	illuminated layer data, zmax [mm]
	0.02	min corr for ratio nx 
	0.02	min corr for ratio ny 
	0.02	min corr for ratio npix 
	0.02	sum of gv    
	33	min for weighted correlation 
	0.02	tolerance to epipolar line [mm]
	
	"""
	f = open ('parameters/criteria.par','w')
	f.write("%f\n" % X_lay[0])
	f.write("%f\n" % Zmin_lay[0])
	f.write("%f\n" % Zmax_lay[0])
	f.write("%f\n" % X_lay[1])
	f.write("%f\n" % Zmin_lay[1])
	f.write("%f\n" % Zmax_lay[1])
	f.write("%f\n" % cnx)
	f.write("%f\n" % cny)
	f.write("%f\n" % cn)
	f.write("%f\n" % csumg)
	f.write("%f\n" % corrmin)
	f.write("%f\n" % eps0)
	f.close()
	

# ========================================
def _write_targ_rec_par():
	"""
	targ_rec.par:	parameters for particle detection
	12	grey value threshold 1. image 
	12	grey value threshold 2. image 
	12	grey value threshold 3. image 
	12	grey value threshold 4. image 
	50	tolerable discontinuity in grey values 
	25	min npix, area covered by particle 
	400	max npix, area covered by particle
	5	min npix in x, dimension in pixel 
	20	max npix in x, dimension in pixel 
	5	min npix in y, dimension in pixel 
	20	max npix in y, dimension in pixel 
	100	sum of grey value
	1	size of crosses
	"""
	
	f = fopen('parameters/targ_rec.par','w')

	for i in n.arange(n_img):
                f.write("%d\n" % gvthres[i])
			
	f.write("%d\n" % disco)
	f.write("%d\n" % nnmin)
	f.write("%d\n" % nnmax)
	f.write("%d\n" % nxmin)
	f.write("%d\n" % nxmax)
	f.write("%d\n" % nymin)
	f.write("%d\n" % nymax)
	f.write("%d\n" % sumg_min)
	f.write("%d\n" % cr_sz)
	f.close()
# ========================================
def _write_man_ori_par():
	"""
	man_ori.par:	point number for manual pre-orientation
	28	image 1 p1 on target plate (reference body) 
	48	image 1 p2 
	42	image 1 p3 
	22	image 1 p4
	28	image 2 p1 
	48	image 2 p2 
	42	image 2 p3 
	23	image 2 p4 
	28	image 3 p1 
	48	image 3 p2 
	42	image 3 p3 
	22	image 3 p4 
	28	image 4 p1 
	48	image 4 p2 
	42	image 4 p3 
	22	image 4 p4
	"""
	
	f = open('parameters/man_ori.par','w')
	for i in n.arange(n_img):
		for j in n.arange(4):
			f.write("%d\n" % nr[i,j])
			
	f.close()

