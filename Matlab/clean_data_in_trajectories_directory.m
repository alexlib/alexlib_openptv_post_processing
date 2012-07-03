% clean the data in the folder
cd ('/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)')
filelist = dir('small*.mat')
for i = 1:length(filelist)
    load(filelist(i).name)
    keep(filelist(i).name(1:end-4),'filelist','i')
    save(filelist(i).name)
end

filelist = dir('large*.mat')
for i = 1:length(filelist)
    load(filelist(i).name)
    keep(filelist(i).name(1:end-4),'filelist','i')
    save(filelist(i).name)
end