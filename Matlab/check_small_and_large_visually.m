% how to check the small and large particles visually

directory = 'c:\PTV\PyPTV\test1\res_189(10680-10819)small\'
small_10680 = ptv_is_to_traj(directory);
directory = 'c:\PTV\PyPTV\test1\res_189(10680-10819)large\'
large_10680 = ptv_is_to_traj(directory);
plot_together_large_small(large_10680,small_10680)
save large_10680 large_10680
save small_10680 small_10680

% after that:
load large_10680
load small_10680
plot_together_large_small(large_10680,small_10680)

