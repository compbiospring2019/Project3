# functions used to process training data and make predictions about testing data

import utils


# returns the mean and standard deviation of all 3 distributions
# distributions['C']['mean'] = mean of distribution representing P(Xi | Y = C)
# distributions['H']['std_dev'] = standard deviation of distribution representing P(Xi | Y = H)
def train(fasta_train):
    return None


def main():
    # get filenames
    fasta, ss = utils.parse_args()
    # split data into training and testing sets
    fasta_train, fasta_test = utils.split_files(fasta, ss)
    # train, creating 3 Gaussian distributions, 1 for each possible SS label
    # each distribution will be stored as pairs of mean and standard deviation
    distributions = train(fasta_train)
