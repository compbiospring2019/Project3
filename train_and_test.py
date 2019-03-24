# functions used to process training data and make predictions about testing data
import utils
from math import sqrt, exp

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
#returns values used to calculate accuracy
def test(pssm_files, pssm_dir, ss_dir):
    #metrics
    correct_c = 0
    correct_e = 0
    correct_h = 0
    total_c = 0
    total_e = 0
    total_h = 0
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
                else:
                    #not out of bounds
                    row = pssm[row_num + row_offset]
                    feature_values.update([row[k] for k in row.keys() if k != 'this-acid'])
            #all feature values recorded
            #now find the maximum probability these features were observed given C, E, and H
            gnb_c = maximum_likelihood(feature_values, "C.dist")
            gnb_e = maximum_likelihood(feature_values, "E.dist")
            gnb_h = maximum_likelihood(feature_values, "H.dist")
            #prediction
            if max([gnb_c, gnb_e, gnb_h]) == gnb_c:
                prediction = 'C'
            elif max([gnb_c, gnb_e, gnb_h]) == gnb_e:
                prediction = 'E'
            else:
                prediction = 'H'
            actual = ss[rownum]
            if actual == 'C':
                total_c += 1
                if prediction == 'C':
                    correct_c += 1
            if actual == 'E':
                total_e += 1
                if prediction == 'E':
                    correct_e += 1
            else:
                total_h += 1
                if prediction == 'H':
                    correct_h += 1
        return total_c, total_e, total_h, correct_c, correct_e, correct_h

#max_prob - maximum probability the given feature values were observed given the specified class label
def maximum_likelihood(feature_values, dist_file, dir="."):
    dist_file = os.path.join(dir, dist_file)
    with open(dist_file, 'r') as f:
        prior = float(f.readline())
        max_prob = 0.0
        for line in f:
            mean, std_dev = [float(x) for x in line.split()]
            probs = [gnb(value, mean, std_dev) for value in feature_values]
            prob = 1.0
            for p in probs:
                prob *= p
            if prob > max_prob:
                max_prob = prob
    return max_prob

def gnb(value, mean, std_dev):
    return 1 / sqrt(2 * 3.14159 * std_dev ** 2) * exp(-1 * (value - mean) ** 2 / (2 * std_dev ** 2))

def main():
    # get filenames
    pssm_list, ss_list, pssm_dir, ss_dir = utils.parse_args()
    # split data into training and testing sets
    pssm_train, pssm_test = utils.split_files(pssm_list, ss_list)
    # Train the model
    train(pssm_train, pssm_dir, ss_dir)


if __name__ == '__main__':
    main()
