from umap import UMAP
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
#%matplotlib inline
import numpy as np
cell_embeddings = np.load('cell_embeddings.npy')
# cell_type = np.load('cell_type.npy',allow_pickle='T')
cell_type_raw = np.loadtxt('higashi_input/cell_type.txt')
cell_type = cell_type_raw * 100
cell_type[cell_type > 13]= 13
cell_type[cell_type < 0]= 0


s_para=6
vec0 = UMAP(n_components=3,n_neighbors=3,min_dist=0.01,metric='correlation').fit_transform(cell_embeddings)
sns.scatterplot(x=vec0[:, 0], y=vec0[:, 1], hue=cell_type, s=s_para, linewidth=0,palette='magma')
plt.tight_layout()
plt.savefig("umap_vec_comp3_n3_d0.01_cor.pdf")
plt.show()
plt.close()
np.savetxt('umap_vec_comp3_n3_d0.01_cor.txt', vec0, delimiter='\t')

vec2 = UMAP(n_components=3,n_neighbors=3,min_dist=0.01).fit_transform(cell_embeddings)
sns.scatterplot(x=vec2[:, 0], y=vec2[:, 1], hue=cell_type, s=s_para, linewidth=0,palette='magma')
plt.tight_layout()
plt.savefig("umap_vec_comp3_n3_d0.01.png")
plt.savefig("umap_vec_comp3_n3_d0.01.pdf")
plt.show()
plt.close()
np.savetxt('umap_vec_comp3_n3_d0.01.txt', vec2, delimiter='\t')

vec1 = UMAP(n_components=5,n_neighbors=5,min_dist=0.01,metric='correlation').fit_transform(cell_embeddings)
sns.scatterplot(x=vec1[:, 0], y=vec1[:, 1], hue=cell_type, s=s_para, linewidth=0,palette='magma')
plt.tight_layout()
plt.savefig("umap_vec_comp5_n5_d0.01_cor.pdf")
plt.savefig("umap_vec_comp5_n5_d0.01_cor.png")
plt.show()
plt.close()
np.savetxt('umap_vec_comp5_n5_d0.01_cor.txt', vec1, delimiter='\t')

vec3 = UMAP(n_components=15,n_neighbors=15,min_dist=0.01,metric='correlation').fit_transform(cell_embeddings)
sns.scatterplot(x=vec3[:, 0], y=vec3[:, 1], hue=cell_type, s=s_para, linewidth=0,palette='magma')
plt.tight_layout()
plt.savefig("umap_vec_comp15_n15_d0.01_cor.pdf")
plt.show()
plt.close()
np.savetxt('umap_vec_comp15_n15_d0.01_cor.txt', vec3, delimiter='\t')

