import os

def split_file(file_path, num_files):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        num_lines = len(lines)
        lines_per_file = num_lines // num_files
        for i in range(num_files):
            start = i * lines_per_file
            end = start + lines_per_file
            if i == num_files - 1:
                end = num_lines
            with open(f'{file_path}.{i}', 'w') as f:
                f.writelines(lines[start:end])

file_path = input('Enter path to file: ')
num_files = int(input('Enter number of files to split into: '))
split_file(file_path, num_files)
