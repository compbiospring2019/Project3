# functions used to process training data and make predictions about testing data
import utils

empty_row = {'A': -1, 'C': -1, 'E': -1, 'D': -1, 'G': -1, 'I': -1, 'H': -1, 'K': -1, 'F': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': -1, 'Y': -1}


# returns the mean and standard deviation of all 3 distributions
# distributions['C']['mean'] = mean of distribution representing P(Xi | Y = C)
# distributions['H']['std_dev'] = standard deviation of distribution representing P(Xi | Y = H)
def train(pssm_files, pssm_dir, ss_dir):
    feature_matrix = build_feature_matrix(pssm_files, pssm_dir, ss_dir)
    return None


def build_feature_matrix(pssm_files, pssm_dir, ss_dir):
    """
    Builds a feature matrix based on PSSM and SS files
    """
    feature_matrix = []

    for pssm_file in pssm_files:
        # For each training file, read in the PSSM matrix and the SS file
        pssm = utils.read_pssm(pssm_file, pssm_dir)
        ss = utils.read_sequence(pssm_file.replace('.pssm', '.ss'), ss_dir)
        for row_num in range(len(pssm)):
            # For each amino acid in the PSSM, build a line for the feature matrix
            feature = {'ss': ss[row_num]}
            for row_offset in range(-2, 3):
                if row_num + row_offset < 0:
                    # We're at the top of the PSSM
                    feature[row_offset] = empty_row
                elif row_num + row_offset >= len(pssm):
                    # We're at the bottom of the PSSM
                    feature[row_offset] = empty_row
                else:
                    # We're somewhere in the middle
                    feature[row_offset] = pssm[row_num + row_offset]
            feature_matrix.append(feature)
    return feature_matrix


def main():
    # get filenames
    pssm_list, ss_list, pssm_dir, ss_dir = utils.parse_args()
    # split data into training and testing sets
    pssm_train, pssm_test = utils.split_files(pssm_list, ss_list)
    # train, creating 3 Gaussian distributions, 1 for each possible SS label
    # each distribution will be stored as pairs of mean and standard deviation
    distributions = train(pssm_train, pssm_dir, ss_dir)


if __name__ == '__main__':
    main()
