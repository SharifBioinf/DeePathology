from keras.layers import Input, Dense, BatchNormalization, Dropout, GaussianNoise, GaussianDropout
from keras.models import Model
import keras.backend as backend
import numpy as np
from sklearn.preprocessing import LabelEncoder, normalize
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_squared_error as mse
from keras.utils import np_utils
from keras.callbacks import CSVLogger, History
import pandas as pd

local_dataset_folder = './Data/'
local_results_folder = './Data/Result/'
dropout_model_name = "-dropout.csv"
gaussian_model_name = "-gaussian.csv"
model_spec = "optimal_" #or hyperopt_


def contractive_dropout_autoencoder(machine_name, local_data_folder, local_result_folder, model_specific, seed=2018):
    seed = seed
    np.random.seed(seed=seed)
    dataset_folder = local_data_folder
    df_m_rna_address = dataset_folder + "fpkm.csv"
    df_mi_rna_address = dataset_folder + "miRNA.csv"
    df_tissue_address = dataset_folder + "tissue.csv"
    df_disease_address = dataset_folder + "disease.csv"

    df_m_rna = np.loadtxt(df_m_rna_address, delimiter=",")
    df_mi_rna = np.loadtxt(df_mi_rna_address, delimiter=",")
    df_tissue = np.ravel(pd.DataFrame.as_matrix(pd.read_csv(df_tissue_address, delimiter=",", header=None)))
    df_disease = np.ravel(pd.DataFrame.as_matrix(pd.read_csv(df_disease_address, delimiter=",", header=None)))

    df_m_rna = normalize(X=df_m_rna, axis=0, norm="max")
    df_mi_rna = normalize(X=df_mi_rna, axis=0, norm="max")

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
    mi_rna = mi_rna[indices]

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

    print("data loading has just been finished")
    print(m_rna.shape, mi_rna.shape, categorical_tissue.shape, categorical_disease.shape)

    batch_size = 64
    nb_epochs = 200

    def create_model():
        inputs = Input(shape=(m_rna.shape[1],), name="inputs")
        inputs_noise = GaussianNoise(stddev=0.025)(inputs)
        inputs_noise = GaussianDropout(rate=0.025 ** 2 / (1 + 0.025 ** 2))(inputs_noise)
        inputs_0 = BatchNormalization(name="inputs_0")(inputs_noise)
        inputs_0 = Dropout(rate=0.0, name='dropout_1')(inputs_0)
        inputs_1 = Dense(1024, activation="softplus", name="inputs_1")(inputs_0)
        inputs_2 = BatchNormalization(name="inputs_2")(inputs_1)
        inputs_2 = Dropout(rate=0.0, name='dropout_2')(inputs_2)
        inputs_3 = Dense(256, activation="softplus", name="inputs_3")(inputs_2)
        inputs_4 = BatchNormalization(name="inputs_4")(inputs_3)
        inputs_4 = Dropout(rate=0.25, name='dropout_3')(inputs_4)

        encoded = Dense(units=12, activation='relu', name='encoded')(inputs_4)

        inputs_5 = Dense(512, activation="linear", name="inputs_5")(encoded)
        inputs_5 = Dropout(rate=0.25, name='dropout_4')(inputs_5)

        decoded_tcga = Dense(units=m_rna.shape[1], activation='relu', name="m_rna")(inputs_5)
        decoded_micro_rna = Dense(units=mi_rna.shape[1], activation='relu', name="mi_rna")(inputs_5)
        cl_0 = Dense(units=categorical_tissue.shape[1], activation="softmax", name="cl_tissue")(encoded)
        cl_2 = Dense(units=categorical_disease.shape[1], activation="softmax", name="cl_disease")(encoded)

        scae = Model(inputs=inputs, outputs=[decoded_tcga, decoded_micro_rna, cl_0, cl_2])

        lambda_value = 9.5581e-3

        def contractive_loss(y_pred, y_true):
            mse = backend.mean(backend.square(y_true - y_pred), axis=1)

            w = backend.variable(value=scae.get_layer('encoded').get_weights()[0])  # N inputs N_hidden
            w = backend.transpose(w)  # N_hidden inputs N
            h = scae.get_layer('encoded').output
            dh = h * (1 - h)  # N_batch inputs N_hidden

            # N_batch inputs N_hidden * N_hidden inputs 1 = N_batch inputs 1
            contractive = lambda_value * backend.sum(dh ** 2 * backend.sum(w ** 2, axis=1), axis=1)

            return mse + contractive

        scae.compile(optimizer='nadam',
                     loss=[contractive_loss, "mse", "cosine_proximity", "cosine_proximity"],
                     loss_weights=[0.001, 0.001, 0.5, 0.5],
                     metrics={"m_rna": ["mae", "mse"], "mi_rna": ["mae", "mse"], "cl_tissue": "acc", "cl_disease": "acc"})

        return scae

    model = create_model()

    result_folder = local_result_folder + model_specific

    file_name = "best-scae-dropout.log"

    csv_logger = CSVLogger(result_folder + file_name)
    history = History()
    model.fit(m_rna_train, [m_rna_train, mi_rna_train, categorical_tissue_train, categorical_disease_train], batch_size=batch_size, epochs=nb_epochs,
              callbacks=[csv_logger, history],
              validation_data=(m_rna_test, [m_rna_test, mi_rna_test, categorical_tissue_test, categorical_disease_test]), verbose=2)

    print(history.history.keys())
    print("fitting has just been finished")
    # save the model and encoded-layer output
    model.save(filepath=result_folder + "scae-dropout.h5")

    layer_name = "encoded"
    encoded_layer_model = Model(inputs=model.input, outputs=model.get_layer(layer_name).output)
    encoded_output = encoded_layer_model.predict(df_m_rna)
    np.savetxt(X=encoded_output, fname=result_folder + "encoded_scae_dropout.csv", delimiter=",")

    # save the result and prediction value

    data_pred = model.predict(m_rna_test, batch_size=batch_size, verbose=2)
    np.savetxt(X=m_rna_test, fname=result_folder + "tcga_genes_scae_dropout.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=mi_rna_test, fname=result_folder + "micro_rna_scae_dropout.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=categorical_tissue_test, fname=result_folder + "categorical_tissue_scae_dropout.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=categorical_disease_test, fname=result_folder + "categorical_disease_scae_dropout.csv", delimiter=",", fmt='%1.3f')

    np.savetxt(X=data_pred[0], fname=result_folder + "tcga_genes_scae_dropout_pred.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=data_pred[1], fname=result_folder + "micro_rna_scae_dropout_pred.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=data_pred[2], fname=result_folder + "categorical_tissue_scae_dropout_pred.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=data_pred[3], fname=result_folder + "categorical_disease_scae_dropout_pred.csv", delimiter=",", fmt='%1.3f')


    singlecell_score = np.c_[mae(df_singlecell[0:500, ], singlecell_pred[0]), mse(df_singlecell[0:500, ], singlecell_pred[0])]
    print(singlecell_score)
    np.savetxt(X=singlecell_score, fname=result_folder + "sc_scae_dropout.csv", delimiter=",", fmt="%s")

    print("prediction process has just been finished")

    for i in range(1, 51):
        print(i)
        dropout_factor = 1 - np.divide(i, 100)
        dropout_matrix = np.random.binomial(n=1, p=dropout_factor, size=m_rna.shape)
        m_rna_dropout = np.multiply(m_rna, dropout_matrix)
        m_rna_temp_test = m_rna_dropout[9750:10750, ]
        score = model.evaluate(m_rna_temp_test, [m_rna_temp_test, mi_rna_test, categorical_tissue_test, categorical_disease_test], verbose=0, batch_size=batch_size)
        print(score)

        with open(result_folder + 'dropout-Dropout-CAE.txt', 'ab') as file:
            np.savetxt(file, score, delimiter=",")

    print("dropout has just been finished")

    for i in range(1, 26):
        print(i)
        noise_factor = np.divide(i, 100)
        noise_matrix = noise_factor * np.random.normal(loc=0.0, scale=1.0, size=m_rna.shape)
        m_rna_noisy = m_rna + noise_matrix
        m_rna_temp_test = m_rna_noisy[9750:10750, ]
        score = model.evaluate(m_rna_temp_test, [m_rna_temp_test, mi_rna_test, categorical_tissue_test, categorical_disease_test], verbose=0, batch_size=batch_size)
        print(score)

        with open(result_folder + 'gaussian-Dropout-CAE.txt', 'ab') as file:
            np.savetxt(file, score, delimiter=",")

    print("gaussian has just been finished")


contractive_dropout_autoencoder(local_data_folder=local_dataset_folder, local_result_folder=local_results_folder, model_specific=model_spec)

print("run has just been finished")
