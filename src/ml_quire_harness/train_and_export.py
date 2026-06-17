import numpy as np
from sklearn.datasets import load_digits
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

np.random.seed(0)

# Real dataset: 8x8 handwritten digit images, 10 classes
digits = load_digits()
X, y = digits.data, digits.target  # X: (1797, 64)

scaler = StandardScaler()
X = scaler.fit_transform(X)  # normalize -> bounded, ~N(0,1) per-feature, matches "ML tensor" assumption

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)

# Train a small, real MLP (one hidden layer, 32 units) -- actually trained, not synthetic
clf = MLPClassifier(hidden_layer_sizes=(32,), activation='relu', max_iter=3000,
                     random_state=0, alpha=1e-4)
clf.fit(X_train, y_train)
print("train acc:", clf.score(X_train, y_train))
print("test acc :", clf.score(X_test, y_test))

W1, b1 = clf.coefs_[0], clf.intercepts_[0]   # (64,32), (32,)
W2, b2 = clf.coefs_[1], clf.intercepts_[1]   # (32,10), (10,)

print("W1 shape", W1.shape, "W2 shape", W2.shape)
print("|W1| range:", np.abs(W1).min(), np.abs(W1).max(), "std:", W1.std())
print("|W2| range:", np.abs(W2).min(), np.abs(W2).max(), "std:", W2.std())
print("|X_test| range:", np.abs(X_test).min(), np.abs(X_test).max(), "std:", X_test.std())

# Compute real hidden-layer activations (post-ReLU) on real test images -> these feed layer 2
Z1 = X_test @ W1 + b1
A1 = np.maximum(Z1, 0)  # ReLU activations: real, actually-occurring values
print("|A1| (hidden activations) range:", A1.min(), A1.max(), "std:", A1.std(), "frac_zero:", (A1==0).mean())

# Build the set of dot products to test:
# Layer 1: z = x . W1[:,j]   for x in a sample of test images, j over all 32 hidden units
# Layer 2: z = a . W2[:,k]   for a = real hidden activations, k over all 10 output units
rng = np.random.default_rng(0)
n_samples_l1 = 40
idx1 = rng.choice(X_test.shape[0], n_samples_l1, replace=False)

pairs = []  # (layer_tag, x_vector, w_vector)
for i in idx1:
    x = X_test[i]
    for j in range(W1.shape[1]):
        pairs.append(("L1", x, W1[:, j]))

n_samples_l2 = 40
idx2 = rng.choice(A1.shape[0], n_samples_l2, replace=False)
for i in idx2:
    a = A1[i]
    for k in range(W2.shape[1]):
        pairs.append(("L2", a, W2[:, k]))

print("total dot products:", len(pairs))

# Export to a flat text file for the C++ harness: layer_tag, len, x..., w...
with open("dotpairs.txt", "w") as f:
    f.write(f"{len(pairs)}\n")
    for tag, x, w in pairs:
        f.write(f"{tag} {len(x)}\n")
        f.write(" ".join(f"{v:.17g}" for v in x) + "\n")
        f.write(" ".join(f"{v:.17g}" for v in w) + "\n")

print("wrote dotpairs.txt")
