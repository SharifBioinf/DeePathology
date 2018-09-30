from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import LabelEncoder, normalize
import numpy as np
import pickle
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials, space_eval
from sklearn.linear_model import SGDClassifier
import pandas as pd

seed = 0
np.random.seed(seed)

dataset_folder = "./Data/"
df_m_rna_address = dataset_folder + "fpkm.csv"
df_disease_address = dataset_folder + "disease.csv"

df_m_rna = np.loadtxt(df_m_rna_address, delimiter=",")
df_m_rna = normalize(X=df_m_rna, axis=0, norm="max")

df_disease = np.ravel(pd.DataFrame.as_matrix(pd.read_csv(df_disease_address, delimiter=",", header=None)))

label_encoder_tissue = LabelEncoder()
label_encoder_tissue.fit(df_disease)
encoded_disease= label_encoder_tissue.transform(df_disease)

m_rna = df_m_rna


def hyperopt_train_test(params):
    t = params['type']
    del params['type']

    if t == 'knn':
        clf = KNeighborsClassifier(n_jobs=-1, n_neighbors=params["n_neighbors_knn"], weights=params["weights_knn"])
    elif t == "extra_tree_classifier":
        clf = ExtraTreesClassifier(n_jobs=-1, n_estimators=params["n_estimators_etc"],
                                   max_depth=params["max_depth_etc"], criterion=params["criterion_etc"])
    elif t == "random_forest_classifier":
        clf = RandomForestClassifier(n_jobs=-1, n_estimators=params["n_estimators_rfc"],
                                     max_depth=params["max_depth_rfc"], criterion=params["criterion_rfc"])
    elif t == "sgd_classifier":
        clf = SGDClassifier(n_jobs=-1, penalty="elasticnet", warm_start=True, loss=params["loss_sgd"], l1_ratio=params["l1_ratio_sgd"])
    else:
        return 0
    return cross_val_score(clf, m_rna, encoded_disease, cv=5).mean()


space = hp.choice('classifier_type', [

    {
        'type': 'knn',
        'n_neighbors_knn': hp.choice('n_neighbors_knn', range(5, 55, 5)),
        'weights_knn': hp.choice('weights_knn', ["uniform", "distance"]),
    },
    {
        'type': 'extra_tree_classifier',
        'n_estimators_etc': hp.choice('n_estimators_etc', range(20, 301, 20)),
        'criterion_etc': hp.choice('criterion_etc', ["gini", "entropy"]),
        'max_depth_etc': hp.choice('max_depth_etc', range(1, 51, 10)),
    },
    {
        'type': 'random_forest_classifier',
        'n_estimators_rfc': hp.choice('n_estimators_rfc', range(20, 301, 20)),
        'criterion_rfc': hp.choice('criterion_rfc', ["gini", "entropy"]),
        'max_depth_rfc': hp.choice('max_depth_rfc', range(1, 51, 10)),
    },
    {
        'type': 'sgd_classifier',
        'loss_sgd': hp.choice('loss_sgd', ["hinge", "log", "modified_huber", "squared_hinge", "perceptron"]),
        'l1_ratio_sgd': hp.uniform('l1_ratio_sgd', 0.1, 0.9)
    }

])


def f(params):
    acc = hyperopt_train_test(params.copy())
    print(acc)
    return {"loss": -acc, "status": STATUS_OK}


def run_trials():
    trials_step = 1  # how many additional trials to do after loading saved trials. 1 = save after iteration
    num_trials = 1  # initial max_trials. put something small to not have to wait  # result_folder = '/s/chopin/a/grad/asharifi/e/Behrooz/Results/tcga/'
    result_folder = "/data/Result/"

    try:  # try to load an already saved trials object, and increase the max
        file_trials = open(result_folder + 'trials_classification-cl_disease.obj', 'rb')
        trials = pickle.load(file_trials)

        print("Found saved Trials! Loading...")
        num_trials = len(trials.trials) + trials_step
        print("Rerunning from {} trials to {} (+{}) trials".format(len(trials.trials), num_trials, trials_step))

    except:  # create a new trials object and start searching
        trials = Trials()

    best = fmin(fn=f, space=space, algo=tpe.suggest, max_evals=num_trials, trials=trials)
    print(space_eval(space, best))

    print("Best:", best)
    print(trials.best_trial)

    # save the trials object
    file_trials = open(result_folder + 'trials_classification-cl_disease.obj', 'wb')
    pickle.dump(trials, file_trials)


counter = 1
result_folder = "/data/Result/"
while counter < 100:
    run_trials()
    file_trials = open(result_folder + 'trials_classification-cl_disease.obj', 'rb')
    trials = pickle.load(file_trials)
    counter = len(trials.trials)
