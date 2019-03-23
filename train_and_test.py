# functions used to process training data and make predictions about testing data
import utils

empty_row = {'A': -1, 'C': -1, 'E': -1, 'D': -1, 'G': -1, 'I': -1, 'H': -1, 'K': -1, 'F': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': -1, 'Y': -1}


#writes distributions to files
def train(pssm_files, pssm_dir, ss_dir):
    feature_matrix = build_feature_matrix(pssm_files, pssm_dir, ss_dir)


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

#reads .pssm, .ss, and .dist files
#expected class labels stored in a list, not written to file
#returns overall q3 accuracy as float, as well as miscategorizations as a dictionary
def test(pssm_files, pssm_dir, ss_dir):
    #for each sequence
    for pssm_file in pssm_files:
        pssm = utils.read_pssm(pssm_file, pssm_dir)
        ss = utils.read_sequence(pssm_file.replace('.pssm', '.ss'), ss_dir)
        #for each acid in the sequence
        for row_num in range(len(pssm)):
            #find feature values
            feature_values = []
            for row_offset in range(-2, 3):
                if row_num + row_offset < 0 or row_num + row_offset >= len(pssm):
                    #out of bounds
                    feature_values.update([-1] * 10)
                else
                    #not out of bounds
                    row = pssm(row_num + row_offset)
                    feature_values.update([row[k] for k in row.keys() if k != 'this-acid'])
            #all feature values recorded

def main():
    #get filenames
    pssm_list, ss_list, pssm_dir, ss_dir = utils.parse_args()
    #split data into training and testing sets
    pssm_train, pssm_test = utils.split_files(pssm, ss)


if __name__ == '__main__':
    main()
