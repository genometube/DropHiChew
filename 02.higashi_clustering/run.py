from higashi.Higashi_wrapper import *
config = "/research/xieyeming1/proj_2023/hicVelocity_20231027/higashi/tanchen_scHiC10x293t_20240505_iced_1M/config.JSON"
higashi_model = Higashi(config)

higashi_model.process_data()

higashi_model.prep_model()
# Stage 1 training
higashi_model.train_for_embeddings()

higashi_model.train_for_imputation_nbr_0()
higashi_model.impute_no_nbr()

higashi_model.train_for_imputation_with_nbr()
higashi_model.impute_with_nbr()

# Visualize embedding results
cell_embeddings = higashi_model.fetch_cell_embeddings()
print (cell_embeddings.shape)


import numpy as np
np.savetxt('cell_embeddings.txt', cell_embeddings, delimiter='\t')
np.save('cell_embeddings.npy', cell_embeddings)

cell_type = higashi_model.label_info['cell_type']
np.savetxt('cell_type.txt', cell_type, delimiter='\t',fmt='%s')
np.save('cell_type.npy', cell_type)

from umap import UMAP
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt

cell_type = higashi_model.label_info['cell_type']
    # fig = plt.figure(figsize=(14, 5))
    # ax = plt.subplot(1, 2, 1)
    # vec = PCA(n_components=2).fit_transform(cell_embeddings)
    # sns.scatterplot(x=vec[:, 0], y=vec[:, 1], hue=cell_type, ax=ax, s=6, linewidth=0)
    # handles, labels = ax.get_legend_handles_labels()
    # labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    # ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
    # ax = plt.subplot(1, 2, 2)
    #
    # vec = UMAP(n_components=5).fit_transform(cell_embeddings)
    # sns.scatterplot(x=vec[:, 0], y=vec[:, 1], hue=cell_type, ax=ax, s=6, linewidth=0)
    # handles, labels = ax.get_legend_handles_labels()
    # labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    # ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
    # plt.tight_layout()
    # plt.savefig("umap_vec_comp5_default.pdf")
    # plt.show()
    # plt.close()
    # np.savetxt('umap_vec_comp5_default.txt', vec, delimiter='\\t')
vec = UMAP(n_components=5).fit_transform(cell_embeddings)
sns.scatterplot(x=vec[:, 0], y=vec[:, 1], hue=cell_type, s=6, linewidth=0)
plt.savefig("umap_vec_comp5_default.pdf")
plt.show()
plt.close()
np.savetxt('umap_vec_comp5_default.txt', vec, delimiter='\\t')

vec1 = UMAP(n_components=5,n_neighbors=5,min_dist=0.01,metric='correlation').fit_transform(cell_embeddings)
sns.scatterplot(x=vec1[:, 0], y=vec1[:, 1], hue=cell_type, s=6, linewidth=0)
plt.savefig("umap_vec_comp5_n5_d0.01_cor.pdf")
plt.show()
plt.close()
np.savetxt('umap_vec_comp5_n5_d0.01_cor.txt', vec1, delimiter='\\t')

vec2 = UMAP(n_components=15,n_neighbors=15,min_dist=0.01,metric='correlation').fit_transform(cell_embeddings)
sns.scatterplot(x=vec2[:, 0], y=vec2[:, 1], hue=cell_type, s=6, linewidth=0)
plt.savefig("umap_vec_comp15_n15_d0.01_cor.pdf")
plt.show()
plt.close()
np.savetxt('umap_vec_comp15_n15_d0.01_cor.txt', vec1, delimiter='\\t')

count = 0
fig = plt.figure(figsize=(6, 2 * 5))
for id_ in np.random.randint(0, 620, 5):
    ori, nbr0, nbr5 = higashi_model.fetch_map("chr3", id_)
    count += 1
    ax = plt.subplot(5, 3, count * 3 - 2)
    ax.imshow(ori.toarray(), cmap='Reds', vmin=0.0, vmax=np.quantile(ori.data, 0.6))
    ax.set_xticks([], [])
    ax.set_yticks([], [])
    if count == 1:
        ax.set_title("raw")

    ax = plt.subplot(5, 3, count * 3 - 1)
    ax.imshow(nbr0.toarray(), cmap='Reds', vmin=0.0, vmax=np.quantile(nbr0.data, 0.95))
    ax.set_xticks([], [])
    ax.set_yticks([], [])
    if count == 1:
        ax.set_title("higashi, k=0")

    ax = plt.subplot(5, 3, count * 3)
    ax.imshow(nbr5.toarray(), cmap='Reds', vmin=0.0, vmax=np.quantile(nbr5.data, 0.95))
    ax.set_xticks([], [])
    ax.set_yticks([], [])
    if count == 1:
        ax.set_title("higashi, k=5")

plt.tight_layout()
plt.savefig("imputation.png")
plt.close()
