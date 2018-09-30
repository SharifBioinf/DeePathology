import csv
import pickle
import numpy as np

results_folder = "./Data/Result/"


def result(result_folder, object_name, file_name):
    file_trials = open(result_folder + object_name, 'rb')
    trials = pickle.load(file_trials)
    best_network = trials.argmin
    best_trials = trials.best_trial
    losses = trials.losses()
    trials_list = trials._trials

    print(best_network, best_trials)
    file_loss = result_folder + file_name + ".txt"
    loss = np.loadtxt(file_loss)
    loss = loss.reshape(-1, 11)
    np.savetxt(result_folder + file_name + ".csv", loss, delimiter=",")

    with open(result_folder + file_name + "-networks" + ".csv", 'w', newline='') as csv_file:
        fieldnames = trials_list[0]["misc"]["vals"].keys()
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        for i in range(len(trials_list)):
            writer.writerow(trials_list[i]["misc"]["vals"])

    np.savetxt(result_folder + file_name + "_losses" + ".csv", losses, delimiter=",")


result(result_folder=results_folder, object_name="trials_scae.obj", file_name="hyperopt-CAE")
result(result_folder=results_folder, object_name="trials_scae_dropout.obj", file_name="hyperopt-Dropout-CAE")
result(result_folder=results_folder, object_name="trials_svae.obj", file_name="hyperopt-VAE")
result(result_folder=results_folder, object_name="trials_svae_dropout.obj", file_name="hyperopt-Dropout-VAE")
