from keras.layers import Input, Dense, BatchNormalization, Dropout, GaussianNoise, GaussianDropout, multiply
from keras.models import Model
import keras.backend as backend
import numpy as np
from os.path import isfile
from sklearn.preprocessing import LabelEncoder, normalize
from keras.utils import np_utils
import pandas as pd
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
import pickle
import os
import time

dataset_folder = "./Data/"
results_folder = "./Data/Result/"
bash_folder = "./code/"


def hyperopt_dropout_cae(data_folder, result_folder, seed=2018):
    seed = seed
    np.random.seed(seed)

    df_m_rna_address = data_folder + "fpkm.csv"
    df_mi_rna_address = data_folder + "miRNA.csv"
    df_tissue_address = data_folder + "tissue.csv"
    df_disease_address = data_folder + "disease.csv"

    df_m_rna = np.loadtxt(df_m_rna_address, delimiter=",")
    df_mi_rna = np.loadtxt(df_mi_rna_address, delimiter=",")
    df_tissue = np.ravel(pd.DataFrame.as_matrix(pd.read_csv(df_tissue_address, delimiter=",", header=None)))
    df_disease = np.ravel(pd.DataFrame.as_matrix(pd.read_csv(df_disease_address, delimiter=",", header=None)))

    label_encoder_tissue = LabelEncoder()
    label_encoder_tissue.fit(df_tissue)
    encoded_tissue = label_encoder_tissue.transform(df_tissue)

    label_encoder_disease = LabelEncoder()
    label_encoder_disease.fit(df_disease)
    encoded_disease = label_encoder_disease.transform(df_disease)

    categorical_tissue = np_utils.to_categorical(encoded_tissue)
    categorical_disease = np_utils.to_categorical(encoded_disease)
    m_rna = df_m_rna
    mi_rna = df_mi_rna

    indices = np.arange(m_rna.shape[0])
    indices = indices[0:10750]
    np.random.shuffle(indices)

    m_rna = m_rna[indices]
    m_rna = normalize(X=m_rna, axis=0, norm="max")

    mi_rna = mi_rna[indices]
    mi_rna = normalize(X=mi_rna, axis=0, norm="max")

    categorical_tissue = categorical_tissue[indices]
    categorical_disease = categorical_disease[indices]

    m_rna_train = m_rna[0:9750, ]
    m_rna_test = m_rna[9750:10750, ]

    mi_rna_train = mi_rna[0:9750, ]
    mi_rna_test = mi_rna[9750:10750, ]

    categorical_tissue_train = categorical_tissue[0:9750, ]
    categorical_tissue_test = categorical_tissue[9750:10750, ]

    categorical_disease_train = categorical_disease[0:9750, ]
    categorical_disease_test = categorical_disease[9750: 10750, ]

    print("data loading has just been finished!")
    print(m_rna.shape, mi_rna.shape, categorical_tissue.shape, categorical_disease.shape)

    space = {

        'layer1': hp.choice('layer1', [512, 1024]),
        'layer3': hp.choice('layer3', [128, 256]),
        'layer4': hp.choice('layer4', [8, 12, 16, 20, 28, 24, 32]),
        'layer5': hp.choice('layer5', [256, 512]),

        'gaussian_noise': hp.choice('gaussian_noise', 0.025 * np.array([0, 1, 2])),
        'gaussian_dropout': hp.choice('gaussian_dropout', 0.025 * np.array([0, 1, 2])),

        'dropout1': hp.choice('dropout1', [0, 0.25, 0.5]),
        'dropout2': hp.choice('dropout2', [0, 0.25, 0.5]),
        'dropout3': hp.choice('dropout3', [0, 0.25, 0.5]),
        'dropout4': hp.choice('dropout4', [0, 0.25, 0.5]),

        'activation1': hp.choice('activation1', ['relu', 'linear', 'softplus', 'elu']),
        'activation3': hp.choice('activation3', ['relu', 'linear', 'softplus', 'elu']),
        'activation4': hp.choice('activation4', ['relu', 'linear', 'softplus', 'elu']),
        'activation5': hp.choice('activation5', ['relu', 'linear', 'softplus', 'elu']),
        'activation6': hp.choice('activation6', ['relu', 'linear', 'softplus', 'elu']),

        "lambda": hp.loguniform('lambda', np.log(1e-6), np.log(1e-2)),

        'batch_size': hp.choice('batch_size', [32, 64, 128]),

        'nb_epochs': hp.choice('nb_epochs', [150, 200, 250]),

        'optimizer': 'adam',

    }

    def objective(params):
        inputs = Input(shape=(m_rna.shape[1],), name="inputs")
        inputs_noise = GaussianNoise(stddev=params['gaussian_noise'])(inputs)
        inputs_noise = GaussianDropout(rate=params['gaussian_dropout'] ** 2 / (1 + params['gaussian_dropout'] ** 2))(inputs_noise)
        inputs_0 = BatchNormalization(name="inputs_0")(inputs_noise)
        inputs_0 = Dropout(rate=params['dropout1'], name='dropout_1')(inputs_0)
        inputs_1 = Dense(units=params['layer1'], activation=params['activation1'], name="inputs_1")(inputs_0)
        inputs_2 = BatchNormalization(name="inputs_2")(inputs_1)
        inputs_2 = Dropout(rate=params['dropout2'], name='dropout_2')(inputs_2)
        inputs_3 = Dense(units=params['layer3'], activation=params['activation3'], name="inputs_3")(inputs_2)
        inputs_4 = BatchNormalization(name="inputs_4")(inputs_3)
        inputs_4 = Dropout(rate=params['dropout3'], name='dropout_3')(inputs_4)

        encoded = Dense(units=params['layer4'], activation=params['activation4'], name='encoded')(inputs_4)
        encoded_softmax = Dense(units=params['layer4'], activation="softmax", name='encoded_softmax')(encoded)
        encoded_attention = multiply([encoded, encoded_softmax])

        inputs_5 = Dense(params['layer5'], activation=params['activation5'], name="inputs_5")(encoded_attention)
        inputs_5 = Dropout(rate=params['dropout4'], name='dropout_4')(inputs_5)

        decoded_tcga = Dense(units=m_rna.shape[1], activation=params['activation6'], name="m_rna")(inputs_5)
        decoded_micro_rna = Dense(units=mi_rna.shape[1], activation=params['activation6'], name="mi_rna")(inputs_5)
        decoded_cl_tissue = Dense(units=categorical_tissue.shape[1], activation="softmax", name="cl_tissue")(encoded_attention)
        decoded_cl_disease = Dense(units=categorical_disease.shape[1], activation="softmax", name="cl_disease")(encoded_attention)

        scae_dropout = Model(inputs=inputs, outputs=[decoded_tcga, decoded_micro_rna, decoded_cl_tissue, decoded_cl_disease])

        lambda_value = params['lambda']

        def contractive_loss(y_pred, y_true):
            mse = backend.mean(backend.square(y_true - y_pred), axis=1)

            w = backend.variable(value=scae_dropout.get_layer('encoded').get_weights()[0])  # N inputs N_hidden
            w = backend.transpose(w)  # N_hidden inputs N
            h = scae_dropout.get_layer('encoded').output
            dh = h * (1 - h)  # N_batch inputs N_hidden

            # N_batch inputs N_hidden * N_hidden inputs 1 = N_batch inputs 1
            contractive = lambda_value * backend.sum(dh ** 2 * backend.sum(w ** 2, axis=1), axis=1)

            return mse + contractive

        scae_dropout.compile(optimizer=params['optimizer'], loss=[contractive_loss, "mse", "cosine_proximity", "cosine_proximity"],
                             loss_weights=[1e-3, 1e-3, 5e-1, 5e-1],
                             metrics={"m_rna": ["mae", "mse"], "mi_rna": ["mae", "mse"], "cl_tissue": "acc", "cl_disease": "acc"})

        scae_dropout.fit(m_rna_train, [m_rna_train, mi_rna_train, categorical_tissue_train, categorical_disease_train],
                         batch_size=params['batch_size'], epochs=params['nb_epochs'], verbose=2)

        score = scae_dropout.evaluate(m_rna_test, [m_rna_test, mi_rna_test, categorical_tissue_test, categorical_disease_test], verbose=0,
                                      batch_size=params['batch_size'])

        print(score)

        with open(result_folder + 'hyperopt-Dropout-CAE.txt', 'ab') as file:
            np.savetxt(file, score, delimiter=",")

        return {'loss': np.mean([score[1], score[2], score[5], score[6], score[7], score[8], -score[9], -score[10]]), 'status': STATUS_OK}

    def run_trials(max_trials=100):
        trials_step = 1  # how many additional trials to do after loading saved trials. 1 = save after iteration
        num_trials = 1  # initial max_trials.

        counter = 0
        while counter < max_trials:
            if isfile(result_folder + 'trials_scae_dropout.obj'):
                file_trials = open(result_folder + 'trials_scae_dropout.obj', 'rb')
                trials = pickle.load(file_trials)

                print("Found saved Trials! Loading...")
                num_trials = len(trials.trials) + trials_step

                if len(trials.trials) < max_trials:
                    print("Rerunning from {} trials to {} (+{}) trials".format(len(trials.trials), num_trials, trials_step))
                    best = fmin(fn=objective, space=space, algo=tpe.suggest, max_evals=num_trials, trials=trials)

                    print("Best:", best)
                    print(trials.best_trial)
                    counter = len(trials.trials)

                    # save the trials object
                    file_trials = open(result_folder + 'trials_scae_dropout.obj', 'wb')
                    pickle.dump(trials, file_trials)
                else:
                    break

            else:  # create a new trials object and start searching
                trials = Trials()
                best = fmin(fn=objective, space=space, algo=tpe.suggest, max_evals=num_trials, trials=trials)

                print("Best:", best)
                print(trials.best_trial)
                counter = len(trials.trials)

                # save the trials object
                file_trials = open(result_folder + 'trials_scae_dropout.obj', 'wb')
                pickle.dump(trials, file_trials)

    run_trials(200)

    print("run has just been finished")


def job():
    free_gpu_status = os.system(bash_folder + "test_gpu.sh")

    if free_gpu_status == 0:
        print("GPU is free, start running now")
        hyperopt_dropout_cae(data_folder=dataset_folder, result_folder=results_folder)
        time.sleep(6)


job()
