# functions used to process training data and make predictions about testing data

import utils

def main():
    # get filenames
    pssm_list, ss_list, pssm_dir, ss_dir = utils.parse_args()
    # split data into training and testing sets
    pssm_train, pssm_test = utils.split_files(pssm, ss)
