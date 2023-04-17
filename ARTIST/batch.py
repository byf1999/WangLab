import os, time
def batch(function, file_path, file_format_list, *args):
    '''This function take all the files in the directory as input to run the function, args are the parameters of function'''
    file_list = os.listdir(file_path)
    for file in file_list:
        if os.path.splitext(file)[-1] in file_format_list:
            t0 = time.time()
            print(f'Dealing with {file}...')
            function(file, *args)
            print(f'Finished, cost {time.time() - t0:.2f} seconds\n')