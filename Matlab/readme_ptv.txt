Readme for the PTV folder

the attached cal 2.zip includes manually created calibration in which I skipped
obviously wrong step of Orientation. I manually changed ORI files and then used
dumbbells in the scenes 7000 - 7209 to get some tuning. The calibration works pretty fine
for scenes 42 and 36. It still doesn't work well for scene32. However, the 
important is that this calibration puts camera 3 and 4 in "better" orientation, due
to a small number (10) dots in the calibration images and the negative z direction,
the orientation procedure fails if automatically launched. 


Scene 32 is processed as follows:

PTV is used to detect the particles and give them correspondences with large epi-
polar distance (5 mm). This works because the number of particles is small. 
The tracking also worked, but I have decided to try also somewhat different approach
in this folder I attach all the Matlab files that I used to track the small 
number of particles in Matlab. 

- multi_stage_tracking_rt_is.m
- read_rt_is_files.m
- plot_long_trajectories (see the output snapshot)
- link_trajectories_rbf
 - haitao_linking_criteria
 
 and the folder /rbf_interp