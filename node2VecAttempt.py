import networkx as nx
from node2vec import Node2Vec
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

#Import Graph from MEGENA.
graph = nx.read_edgelist("/home/spirpinias/Desktop/MEGENAgraph")

# Precompute probabilities and generate walks
node2vec = Node2Vec(graph, dimensions=20, walk_length=16, num_walks=100, workers=4) 

# Embed nodes
model = node2vec.fit(window=10, min_count=1, batch_words=4) 

#Embeddings
nodes = [x for x in model.wv.vocab]
embeddings = np.array([model.wv[x] for x in nodes])

tsne = TSNE(n_components=2, random_state=7, perplexity=15)
embeddings_2d = tsne.fit_transform(embeddings)

figure = plt.figure(figsize=(11, 9))

ax = figure.add_subplot(111)

ax.scatter(embeddings_2d[:, 0], embeddings_2d[:, 1])

plt.show()

pca=PCA(n_components=10)
embeddings_PCA=pca.fit_transform(StandardScaler().fit_transform(embeddings))

fig = plt.figure(figsize=(8,5))
sing_vals = np.arange(0,10,1)
plt.plot(sing_vals, pca.explained_variance_ratio_, 'ro-', linewidth=2)
plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance')
