import pickle
import os
import datetime

file_dir = '/home/daniel/Desktop/MRI_Proposal/vlite-observing-progress/progress_script_20210115'
files = os.listdir(file_dir)
files = sorted(files)

for p_file in files:
    p_path = os.path.join(file_dir, p_file)
    with open(p_path, 'rb') as f:
        p_list = pickle.load(f)
    print('File: {}'.format(p_file))
    for pointing in p_list:
        print('    {}'.format(pointing[2]))

