##This was an attempt I made following tutorials on using heat diffusion kernels on graph topology to embed structurally similar nodes of a graph into R^D. 
## The motivation was to utilize the embedded nodes for further ML processing. However, the embeddings were of poor quality.
## Not because of the algorithm, but the graph itself was noisy.

import networkx as nx
from stellargraph.mapper import GraphWaveGenerator
from stellargraph import StellarGraph
from sklearn.decomposition import PCA
import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse.linalg import eigs
import tensorflow as tf
from tensorflow.keras import backend as K
from sklearn.manifold import TSNE

G = nx.read_edgelist("/home/spirpinias/Desktop/MEGENAgraph")
G = StellarGraph.from_networkx(G)
sample_points = np.linspace(0, 100, 50).astype(np.float32)
#degree20 and scales 5,10

degree = 10
scales = [5, 10]

generator = GraphWaveGenerator(G, scales=scales, degree=degree)

embeddings_dataset = generator.flow(
    node_ids=G.nodes(), sample_points=sample_points, batch_size=10, repeat=False
)

embeddings = [x.numpy() for x in embeddings_dataset]

trans_emb = PCA(n_components=3).fit_transform(np.vstack(embeddings))

plt.scatter(
    trans_emb[:, 0], trans_emb[:, 1], cmap="jet", alpha=0.7,
)
plt.title("graphWave 2 PC's")
plt.show()

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(trans_emb[:, 0], trans_emb[:, 1], trans_emb[:, 2],
           linewidths=1, alpha=.7,
           edgecolor='k',
           s = 200)
plt.show()

tsne = TSNE(n_components=2, random_state=7, perplexity=15)
embeddings_2d = tsne.fit_transform(np.vstack(embeddings))

figure = plt.figure(figsize=(11, 9))

ax = figure.add_subplot(111)

ax.scatter(embeddings_2d[:, 0], embeddings_2d[:, 1])

plt.show()
