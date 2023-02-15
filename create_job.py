import glob
import os


def create_job(N_flavor, Lambda, kir, sigma_max, grid, mu, T, outfile, task_number):
    with open(f'tasks/{task_number}', 'w') as f:
        f.write(str(N_flavor) + "\n")
        f.write(str(Lambda) + "\n")
        f.write(str(kir) + "\n")
        f.write(str(sigma_max) + "\n")
        f.write(str(grid) + "\n")
        f.write(str(mu) + "\n")
        f.write(str(T) + "\n")
        f.write(str(outfile) + "\n")


def get_starting_index():
    files = sorted(map(lambda x: int(x.replace(
        "./tasks/", "")), glob.glob("./tasks/*")))

    # rename all files
    # for i in files:
    #     os.rename(i, i+".bak")

    # # rename again and put in correct order starting with one
    # for i, filename in enumerate(files):
    #     os.rename(filename+".bak", "./tasks/"+str(i+1))
    if len(files) != 0:
        return files[-1] + 1
    else:
        return 1


def print_last_task_index():
    print(len(glob.glob("./tasks/*")))
