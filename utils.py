# Utils for Project 3
import os
from random import sample

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


# Read a biological sequence or RSA sequence from a file:
def read_sequence(file_path, dir=None):
    if file_path is None:
        return None
    sequence = ''
    if dir:
        file_path = os.path.join(dir, file_path)
    with open(file_path, 'r') as f:
        # Ignore the title line
        title = f.readline()
        for line in f:
            sequence += line.strip()

    return sequence.upper()


def read_directory_contents(path, file_extension):
    if not os.path.isdir(path):
        # This is not a valid directory!
        raise Exception('Not a valid directory!')

    # Return a list of files with the file extension
    ls_dir = os.listdir(path)
    return [file_name for file_name in ls_dir if file_name.endswith(file_extension)]


def read_pssm(file_path, dir=None):
    if dir:
        file_path = os.path.join(dir, file_path)
    pssm = []
    with open(file_path, 'r') as f:
        # Ignore the title line
        title = f.readline()
        print('Title :{}:'.format(title))
        if title in ['', '\n']:
            title = f.readline()

        # Get the list of amino acids on the top
        headers = f.readline().strip().split()[:20]

        # Now, read each line of the matrix into a dictionary
        for line in f:
            if line in ['', '\n']:
                break
            line_list = line.strip().split()[:22]
            row = {'this-acid': line_list[1]}
            for acid_num in range(len(headers)):
                row[headers[acid_num]] = int(line_list[acid_num + 2])
            pssm.append(row)

    # Returns a list of dictionaries, where each dict is a row of the matrix
    return pssm

#divide the .fasta files into training and testing sets
def split_files(fasta_list, ss_list):
    test_correct_fasta_files(fasta_list, ss_list)
    fasta_train = sample(fasta_list, int(0.75 * len(fasta_list)))
    fasta_test = [fasta_name for fasta_name in fasta_list if fasta_name not in fasta_train]
    return fasta_train, fasta_test

#make sure each .fasta has a corresponding .ss
def test_correct_fasta_files(self, fasta_list, ss_list):
        for fasta_name in fasta_list:
            if fasta_name.replace('.fasta', '.ss') not in ss_list:
                raise Exception('FASTA files don\'t match up with .ss files: {}'.format(fasta_name))
