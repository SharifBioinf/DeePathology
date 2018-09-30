from keras.layers import Input, Dense, BatchNormalization, Dropout, GaussianNoise, GaussianDropout
from keras.models import Model
import keras.backend as backend
import numpy as np
from sklearn.preprocessing import LabelEncoder
from keras.utils import np_utils
from keras.callbacks import CSVLogger, History
import pandas as pd

machine = 'porsche'
local_dataset_folder = './Data/
local_results_folder = './Data/Result/'
dropout_model_name = "-dropout.csv"
gaussian_model_name = "-gaussian.csv"


def cae_one_versus_all(local_data_folder, local_result_folder, index):
    seed = 0
    np.random.seed(seed=seed)
    dataset_folder = local_data_folder
    df_m_rna_in_address = dataset_folder + "m_rna_in_" + str(index) + ".csv"
    df_mi_rna_in_address = dataset_folder + "mi_rna_in" + str(index) + ".csv"
    df_healthy_in_address = dataset_folder + "healthy_in" + str(index) + ".csv"

    df_m_rna_in = np.loadtxt(df_m_rna_in_address, delimiter=",")
    df_mi_rna_in = np.loadtxt(df_mi_rna_in_address, delimiter=",")
    df_healthy_in = np.ravel(pd.DataFrame.as_matrix(pd.read_csv(df_healthy_in_address, delimiter=",", header=None)))

    df_m_rna_out_address = dataset_folder + "m_rna_out_" + str(index) + ".csv"
    df_mi_rna_out_address = dataset_folder + "mi_rna_out_" + str(index) + ".csv"
    df_healthy_out_address = dataset_folder + "healthy_out_" + str(index) + ".csv"

    df_m_rna_out = np.loadtxt(df_m_rna_out_address, delimiter=",")
    df_mi_rna_out = np.loadtxt(df_mi_rna_out_address, delimiter=",")
    df_healthy_out = np.ravel(pd.DataFrame.as_matrix(pd.read_csv(df_healthy_out_address, delimiter=",", header=None)))

    m_rna_l1 = np.loadtxt(dataset_folder + "m_rna_l1.csv", delimiter=",")
    mi_rna_l1 = np.loadtxt(dataset_folder + "mi_rna_l1.csv", delimiter=",")

    m_rna_inv_l1 = 1. / m_rna_l1
    m_rna_inv_l1[np.where(m_rna_inv_l1 == np.inf)] = 0

    df_m_rna_out = np.array([np.kron(m_rna_inv_l1[i], df_m_rna_out[:, i]) for i in df_m_rna_out.shape[1]]).transpose()
    df_mi_rna_out = np.array([np.kron(m_rna_inv_l1[i], df_mi_rna_out[:, i]) for i in df_mi_rna_out.shape[1]]).transpose()

    label_encoder_disease = LabelEncoder()
    label_encoder_disease.fit(df_healthy_out)
    encoded_healthy_out = label_encoder_disease.transform(df_healthy_out)

    categorical_healthy_out = np_utils.to_categorical(encoded_healthy_out)
    m_rna_out = df_m_rna_out
    mi_rna_out = df_mi_rna_out

    print(df_m_rna_out.shape, df_mi_rna_out.shape, categorical_healthy_out.shape)

    batch_size = 64
    nb_epochs = 200

    def create_model():
        inputs = Input(shape=(m_rna_out.shape[1],), name="inputs")
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

        decoded_tcga = Dense(units=m_rna_out.shape[1], activation='relu', name="m_rna_out")(inputs_5)
        decoded_micro_rna = Dense(units=mi_rna_out.shape[1], activation='relu', name="mi_rna_out")(inputs_5)
        healthy = Dense(units=categorical_healthy_out.shape[1], activation="softmax", name="cl_disease")(encoded)

        scae = Model(inputs=inputs, outputs=[decoded_tcga, decoded_micro_rna, healthy])

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
                     metrics={"m_rna_out": ["mae", "mse"], "mi_rna_out": ["mae", "mse"], "cl_tissue": "acc", "cl_disease": "acc"})

        return scae

    model = create_model()

    result_folder = local_result_folder

    file_name = "best-" + str(index) + ".log"

    csv_logger = CSVLogger(result_folder + file_name)
    history = History()
    model.fit(m_rna_out, [m_rna_out, mi_rna_out, categorical_healthy_out], batch_size=batch_size,
              epochs=nb_epochs, allbacks=[csv_logger, history], validation_split=0.2, verbose=2)

    print(history.history.keys())
    print("fitting has just been finished")
    # save the model and encoded-layer output


    # save the result and prediction value

    pred_result = model.evaluate(df_m_rna_in, [df_m_rna_in, df_mi_rna_in, df_healthy_in])
    np.savetxt(X=pred_result, fname=result_folder + "evaluation_result" + str(index) + ".csv", delimiter=",")
    data_pred = model.predict(df_m_rna_in, batch_size=batch_size, verbose=2)

    df_m_rna_in_pred = np.array([np.kron(m_rna_l1[i], data_pred[0][:, i]) for i in range(data_pred[0].shape[1])]).transpose()
    df_mi_rna_in_pred = np.array([np.kron(mi_rna_l1[i], data_pred[1][:, i]) for i in range(data_pred[1].shape[1])]).transpose()

    np.savetxt(X=df_m_rna_in_pred, fname=result_folder + "m_rna_pred" + str(index) + ".csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=df_mi_rna_in_pred, fname=result_folder + "mi_rna" + str(index) + ".csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=data_pred[2], fname=result_folder + "healthy_pred" + str(index) + ".csv", delimiter=",", fmt='%1.3f')


cae_one_versus_all(machine_name=machine, local_data_folder=local_dataset_folder, local_result_folder=local_results_folder, index=0)

print("run has just been finished")
