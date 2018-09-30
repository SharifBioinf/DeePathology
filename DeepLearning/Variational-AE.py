from keras.layers import Input, Dense, BatchNormalization, Lambda, GaussianDropout, GaussianNoise
from keras.models import Model
import keras.backend as backend
import numpy as np
from sklearn.preprocessing import LabelEncoder, normalize
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_squared_error as mse
from keras.utils import np_utils
from keras.callbacks import CSVLogger, History
import pandas as pd

machine = 'porsche'
local_dataset_folder = './Data/'
local_results_folder = './Data/Result/'
dropout_model_name = "-dropout.csv"
gaussian_model_name = "-gaussian.csv"
model_spec = "optimal_"  #or hyperopt_


def variational_autoencoder(local_data_folder, local_result_folder, model_specific, seed=2018):
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

    batch_size = 125
    nb_epochs = 200

    def create_model():
        n_z = 24

        # Q(z|x) -- encoder
        inputs = Input(shape=(m_rna.shape[1],), name="inputs")
        inputs_noise = GaussianNoise(stddev=0.0)(inputs)
        inputs_noise = GaussianDropout(rate=0.0 ** 2 / (1 + 0.0 ** 2))(inputs_noise)
        inputs_0 = BatchNormalization(name="inputs_0")(inputs_noise)
        inputs_1 = Dense(1024, activation='relu', name="inputs_1")(inputs_0)
        inputs_2 = BatchNormalization(name="inputs_2")(inputs_1)
        inputs_3 = Dense(128, activation='relu', name="inputs_3")(inputs_2)
        inputs_4 = BatchNormalization(name="inputs_4")(inputs_3)

        h_q = Dense(n_z, activation="relu", name="h_q")(inputs_4)
        mu = Dense(n_z, activation=None, name="mu")(h_q)
        log_sigma = Dense(n_z, activation=None, name="log_sigma")(h_q)

        def sample_z(args):
            mu_samples, log_sigma_sample = args
            eps = backend.random_normal(shape=(batch_size, n_z), mean=0., stddev=1.)
            return mu_samples + backend.exp(log_sigma_sample / 2) * eps

        # Sample z ~ Q(z|x)
        z = Lambda(sample_z, name="lambda")([mu, log_sigma])

        # P(x|z) -- decoder
        decoder_hidden = Dense(256, activation='linear', name="decoder_hidden")

        decoded_tcga = Dense(m_rna_train.shape[1], activation='relu', name="m_rna")
        decoded_micro_rna = Dense(mi_rna.shape[1], activation='relu', name="mi_rna")
        decoded_cl_0 = Dense(categorical_tissue.shape[1], activation="softmax", name="cl_tissue")
        decoded_cl_2 = Dense(categorical_disease.shape[1], activation="softmax", name="cl_disease")

        h_p = decoder_hidden(z)
        outputs_tcga = decoded_tcga(h_p)
        outputs_micro_rna = decoded_micro_rna(h_p)
        outputs_cl_0 = decoded_cl_0(h_p)
        outputs_cl_2 = decoded_cl_2(h_p)

        lambda_value = 2.5596e-6

        def vae_loss(y_true, y_pred):
            """ Calculate loss = reconstruction loss + KL loss for each data_pred in minibatch """
            # E[log P(x|z)]
            recon = backend.mean(backend.square(y_true - y_pred), axis=1)
            # D_KL(Q(z|x) || P(z|x)); calculate in closed form as both dist. are Gaussian
            kl = 0.5 * backend.sum(backend.exp(log_sigma) + backend.square(mu) - 1. - log_sigma, axis=1)

            return recon + lambda_value * kl

        svae = Model(inputs, [outputs_tcga, outputs_micro_rna, outputs_cl_0, outputs_cl_2])

        svae.compile(optimizer='nadam',
                     loss=[vae_loss, "mse", "cosine_proximity", "cosine_proximity"],
                     loss_weights=[1e-3, 1e-3, 5e-1, 5e-1],
                     metrics={"m_rna": ["mae", "mse"], "mi_rna": ["mae", "mse"], "cl_tissue": "acc", "cl_disease": "acc"})

        # Generator model, generate new data given latent variable z
        d_in = Input(shape=(n_z,))
        d_h = decoder_hidden(d_in)
        d_out = decoded_tcga(d_h)
        generator = Model(d_in, d_out)

        return svae, generator

    model, generator_model = create_model()

    result_folder = local_result_folder + model_specific

    file_name = "best-svae.log"

    csv_logger = CSVLogger(result_folder + file_name)
    history = History()

    model.fit(m_rna_train, [m_rna_train, mi_rna_train, categorical_tissue_train, categorical_disease_train], batch_size=batch_size, epochs=nb_epochs,
              callbacks=[csv_logger, history],
              validation_data=(m_rna_test, [m_rna_test, mi_rna_test, categorical_tissue_test, categorical_disease_test]), verbose=2)

    print(history.history.keys())
    print("fitting has just been finished")

    # save the model and encoded-layer output
    model.save(filepath=result_folder + "svae.h5")

    layer_name = "h_q"
    encoded_layer_model = Model(inputs=model.input, outputs=model.get_layer(layer_name).output)
    encoded_output = encoded_layer_model.predict(df_m_rna)
    np.savetxt(X=encoded_output, fname=result_folder + "encoded_svae.csv", delimiter=",")

    # generate data from normal distribution
    latent_dim = 24  # latent_dim = n_z
    n_samples = 10000
    z_normal = np.random.normal(loc=0.0, scale=1.0, size=(n_samples, latent_dim))
    gen_data = generator_model.predict(z_normal)
    np.savetxt(X=gen_data, fname=result_folder + "gen_data_vae.csv", delimiter=",", fmt='%1.3f')

    # save the result and prediction value

    data_pred = model.predict(m_rna_test, batch_size=batch_size, verbose=2)
    np.savetxt(X=m_rna_test, fname=result_folder + "tcga_genes_svae.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=mi_rna_test, fname=result_folder + "micro_rna_svae.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=categorical_tissue_test, fname=result_folder + "categorical_tissue_svae.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=categorical_disease_test, fname=result_folder + "categorical_disease_svae.csv", delimiter=",", fmt='%1.3f')

    np.savetxt(X=data_pred[0], fname=result_folder + "tcga_genes_svae_pred.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=data_pred[1], fname=result_folder + "micro_rna_svae_pred.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=data_pred[2], fname=result_folder + "categorical_tissue_svae_pred.csv", delimiter=",", fmt='%1.3f')
    np.savetxt(X=data_pred[3], fname=result_folder + "categorical_disease_svae_pred.csv", delimiter=",", fmt='%1.3f')

    print("prediction process has just been finished")

    for i in range(1, 51):
        print(i)
        dropout_factor = 1-np.divide(i, 100)
        dropout_matrix = np.random.binomial(n=1, p=dropout_factor, size=m_rna.shape)
        m_rna_dropout = np.multiply(m_rna, dropout_matrix)
        m_rna_temp_test = m_rna_dropout[9750:10750, ]
        score = model.evaluate(m_rna_temp_test, [m_rna_temp_test, mi_rna_test, categorical_tissue_test, categorical_disease_test], verbose=0, batch_size=batch_size)
        print(score)

        with open(result_folder + 'dropout-VAE.txt', 'ab') as file:
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

        with open(result_folder + 'gaussian-VAE.txt', 'ab') as file:
            np.savetxt(file, score, delimiter=",")

    print("gaussian has just been finished")


variational_autoencoder(local_data_folder=local_dataset_folder, local_result_folder=local_results_folder, model_specific=model_spec)

print("run has just been finished")
