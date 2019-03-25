# functions used to process training data and make predictions about testing data
import utils
from math import sqrt, exp
import os

empty_row = {'A': -1, 'C': -1, 'E': -1, 'D': -1, 'G': -1, 'I': -1, 'H': -1, 'K': -1, 'F': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': -1, 'Y': -1}
acids_list = ['A', 'C', 'E', 'D', 'G', 'I', 'H', 'K', 'F', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


# Writes distributions to files
def train(pssm_files, pssm_dir, ss_dir):
    # Generate the feature matrix
    feature_matrix = build_feature_matrix(pssm_files, pssm_dir, ss_dir)
    # Calculate the mu and sigma and prior values
    model = calculate_model(feature_matrix)
    # Write the model to a file
    write_model(model)


def calculate_model(matrix):
    model = {'C': {}, 'E': {}, 'H': {}}

    for class_label in model.keys():
        # For each class label, calculate sigmas, mus, and prior terms
        features = [row for row in matrix if row['ss'] == class_label]
        model[class_label]['sigma'] = {}
        model[class_label]['mu'] = {}
        model[class_label]['prior'] = float(len(features)) / len(matrix)

        # Calculate sigmas and mus for each feature
        for feature in range(100):
            mu = calc_mu(features, feature)
            model[class_label]['sigma'][feature] = calc_sigma(features, feature, mu)
            model[class_label]['mu'][feature] = mu

    return model


def calc_mu(features, feature_num):
    sum = 0.0
    for feature in features:
        sum += feature[feature_num]
    return sum/len(features)


def calc_sigma(features, feature_num, mu):
    sum = 0.0
    for feature in features:
        sum += (feature[feature_num] - mu) ** 2
    return sqrt(sum / len(features))


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
                    values = empty_row
                elif row_num + row_offset >= len(pssm):
                    # We're at the bottom of the PSSM
                    values = empty_row
                else:
                    # We're somewhere in the middle
                    values = pssm[row_num + row_offset]
                for val_num, acid in enumerate(acids_list):
                    feature[((row_offset + 2) * 20) + val_num] = values[acid]
            feature_matrix.append(feature)
    return feature_matrix


def write_model(model):
    for class_label in model.keys():
        file_name = os.path.join(parent_directory, '{}.dist'.format(class_label))
        with open(file_name, 'w') as file:
            file.write('{}\n'.format(model[class_label]['prior']))
            for feature in range(100):
                file.write('{} {}\n'.format(model[class_label]['mu'][feature], model[class_label]['sigma'][feature]))


# reads .pssm, .ss, and .dist files
# expected class labels stored in a list, not written to file
# returns values used to calculate accuracy
def test(pssm_files, pssm_dir, ss_dir):
    # metrics
    correct_c = 0
    correct_e = 0
    correct_h = 0
    total_c = 0
    total_e = 0
    total_h = 0
    # for each sequence
    for pssm_file in pssm_files:
        pssm = utils.read_pssm(pssm_file, pssm_dir)
        ss = utils.read_sequence(pssm_file.replace('.pssm', '.ss'), ss_dir)
        # for each acid in the sequence
        for row_num in range(len(pssm)):
            # find feature values
            feature_values = []
            for row_offset in range(-2, 3):
                if row_num + row_offset < 0 or row_num + row_offset >= len(pssm):
                    # out of bounds
                    feature_values.extend([-1] * 20)
                else:
                    # not out of bounds
                    row = pssm[row_num + row_offset]
                    feature_values.extend([row[k] for k in acids_list])
            # all feature values recorded
            # now find the maximum probability these features were observed given C, E, and H
            gnb_c = maximum_likelihood(feature_values, "C.dist")
            gnb_e = maximum_likelihood(feature_values, "E.dist")
            gnb_h = maximum_likelihood(feature_values, "H.dist")
            # prediction
            if max([gnb_c, gnb_e, gnb_h]) == gnb_c:
                prediction = 'C'
            elif max([gnb_c, gnb_e, gnb_h]) == gnb_e:
                prediction = 'E'
            else:
                prediction = 'H'
            actual = ss[row_num]
            if actual == 'C':
                total_c += 1
                if prediction == 'C':
                    correct_c += 1
            elif actual == 'E':
                total_e += 1
                if prediction == 'E':
                    correct_e += 1
            else:
                total_h += 1
                if prediction == 'H':
                    correct_h += 1
    return [total_c, total_e, total_h, correct_c, correct_e, correct_h]


dists = {}


# max_prob - maximum probability the given feature values were observed given the specified class label
def maximum_likelihood(feature_values, dist_file, dir="."):
    dist_file = os.path.join(dir, dist_file)
    if dist_file not in dists:
        print('Reading in {}'.format(dist_file))
        dists[dist_file] = {'sigma': {}, 'mu': {}}
        with open(dist_file, 'r') as f:
            dists[dist_file]['prior'] = float(f.readline())
            for line_num in range(100):
                mean, std_dev = [float(x) for x in f.readline().split()]
                dists[dist_file]['sigma'][line_num] = std_dev
                dists[dist_file]['mu'][line_num] = mean

    prob = dists[dist_file]['prior']
    for feat_num in range(100):
        prob *= gnb(feature_values[feat_num], dists[dist_file]['mu'][feat_num], dists[dist_file]['mu'][feat_num])
    return prob


def gnb(value, mean, std_dev):
    return exp(-0.5 * (value - mean) ** 2 / (std_dev ** 2)) / sqrt(2 * 3.14159 * (std_dev ** 2))


def accuracy(metrics):
    print("metrics: " + str(metrics))
    print("Q3 Accuracy")
    print("-----------")
    if metrics[0] == 0:
        print("No C's observed in testing set.")
    else:
        print("C: " + str(float(metrics[3]) / metrics[0]))
    if metrics[1] == 0:
        print("No E's observed in testing set.")
    else:
        print("E: " + str(float(metrics[4]) / metrics[1]))
    if metrics[0] == 0:
        print("No H's observed in testing set.")
    else:
        print("H: " + str(float(metrics[5]) / metrics[2]))
    print("Overall: " + str(float(sum(metrics[3:6])) / sum(metrics[0:3])))

def classify(pssm_classify):
    predictions = []
    for row_num in range(len(pssm_classify)):
        # find feature values
        feature_values = []
        for row_offset in range(-2, 3):
            if row_num + row_offset < 0 or row_num + row_offset >= len(pssm_classify):
                # out of bounds
                feature_values.extend([-1] * 20)
            else:
                # not out of bounds
                row = pssm_classify[row_num + row_offset]
                feature_values.extend([row[k] for k in acids_list])
        # all feature values recorded
        # now find the maximum probability these features were observed given C, E, and H
        gnb_c = maximum_likelihood(feature_values, "C.dist")
        gnb_e = maximum_likelihood(feature_values, "E.dist")
        gnb_h = maximum_likelihood(feature_values, "H.dist")
        # prediction
        if max([gnb_c, gnb_e, gnb_h]) == gnb_c:
            prediction = 'C'
        elif max([gnb_c, gnb_e, gnb_h]) == gnb_e:
            prediction = 'E'
        else:
            prediction = 'H'
        predictions.append(prediction)
    print("Prediction")
    print("----------")
    print("".join(predictions))

def main():
    # get filenames
    pssm_list, ss_list, pssm_dir, ss_dir, pssm_classify = utils.parse_args()
    # split data into training and testing sets
    pssm_train, pssm_test = utils.split_files(pssm_list, ss_list)
    # Train the model
    train(pssm_train, pssm_dir, ss_dir)
    print('Trained the model. Now for testing...')
    # test
    metrics = test(pssm_test, pssm_dir, ss_dir)
    # accuracy
    accuracy(metrics)
    #classify
    classify(pssm_classify)

if __name__ == '__main__':
    main()
