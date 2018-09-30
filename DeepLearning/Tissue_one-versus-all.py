import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder


local_dataset_folder = './Data/'


def tissue_one_versus_all(local_data_folder):
    dataset_folder = local_data_folder
    df_m_rna_address = dataset_folder + "fpkm.csv"
    df_mi_rna_address = dataset_folder + "miRNA.csv"
    df_tissue_address = dataset_folder + "tissue.csv"
    df_disease_address = dataset_folder + "disease.csv"

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

    normal_label_index = np.where(np.unique(df_disease) == "Normal")
    disease_label = np.zeros(shape=(df_m_rna.shape[0],))
    disease_label[np.where(encoded_disease == normal_label_index[0])] = 1
    disease_label.astype(np.int64)

    np.savetxt(X=disease_label, fname=dataset_folder + "healthy_full_data" + ".csv", delimiter=",", fmt='%1.3f')

    for i in range(len(np.unique(encoded_tissue))):
        print(i)
        temp_index = np.where(encoded_tissue == i)

        temp_m_rna_out = np.squeeze(df_m_rna[np.where(encoded_tissue != i), :], axis=0)
        temp_m_rna_in = np.squeeze(df_m_rna[np.where(encoded_tissue == i), :], axis=0)

        temp_mi_rna_out = np.squeeze(df_mi_rna[np.where(encoded_tissue != i), :], axis=0)
        temp_mi_rna_in = np.squeeze(df_mi_rna[np.where(encoded_tissue == i), :], axis=0)

        temp_disease_out = disease_label[np.where(encoded_tissue != i)]
        temp_disease_in = disease_label[np.where(encoded_tissue == i)]

        np.savetxt(X=temp_m_rna_in, fname=dataset_folder + "m_rna_in_" + str(i) + ".csv", delimiter=",", fmt='%1.3f')
        np.savetxt(X=temp_m_rna_out, fname=dataset_folder + "m_rna_out_" + str(i) + ".csv", delimiter=",", fmt='%1.3f')
        np.savetxt(X=temp_mi_rna_in, fname=dataset_folder + "mi_rna_in_" + str(i) + ".csv", delimiter=",", fmt='%1.3f')
        np.savetxt(X=temp_mi_rna_out, fname=dataset_folder + "mi_rna_out_" + str(i) + ".csv", delimiter=",", fmt='%1.3f')
        np.savetxt(X=temp_disease_in, fname=dataset_folder + "healthy_in_" + str(i) + ".csv", delimiter=",", fmt='%1.3f')
        np.savetxt(X=temp_disease_out, fname=dataset_folder + "healthy_out_" + str(i) + ".csv", delimiter=",", fmt='%1.3f')
        np.savetxt(X=temp_index[0], fname=dataset_folder + "index_" + str(i) + ".csv", delimiter=",", fmt='%1.3f')


tissue_one_versus_all(local_data_folder=local_dataset_folder)
