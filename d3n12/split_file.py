import argparse

def split_file(file_path, file1_path, file2_path):
    with open(file_path) as f:
        lines = f.readlines()
        half = len(lines) // 2
        with open(file1_path, 'w') as f:
            f.writelines(lines[:half])
        with open(file2_path, 'w') as f:
            f.writelines(lines[half:])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split a text file into two files.')
    parser.add_argument('file_path', type=str, help='path to the file to be split')
    parser.add_argument('file1_path', type=str, help='path to the first output file')
    parser.add_argument('file2_path', type=str, help='path to the second output file')
    args = parser.parse_args()

split_file(args.file_path, args.file1_path, args.file2_path)

