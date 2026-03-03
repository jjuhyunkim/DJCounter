# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: plot_python
#     language: python
#     name: plot_python
# ---

# %%
import statsmodels.api as sm
import os
import pickle
import pandas as pd
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import argparse

# %%
script_path = os.path.abspath(__file__)
print(script_path)

script_dir = os.path.dirname(script_path)
print(script_dir)

# %%
parser = argparse.ArgumentParser()

parser.add_argument("--dj_count", type=str, required=True,
                    help="DJ count for GRCh38")

parser.add_argument("--plt_prefix", type=str,
                    default="dj.distribution",
                    help="Plot prefix for output files")

parser.add_argument("--model", type=str,
                    default=f"{script_dir}/GRCh38_fast_accurate_to_CHM13.pkl",
                    help="Model file for prediction")

args = parser.parse_args()
dj_count = float(args.dj_count)
plt_prefix = args.plt_prefix
model = args.model

print("DJ count:", dj_count)
print("Plot prefix:", plt_prefix)

# %%
with open(model, "rb") as f:
    model = pickle.load(f)

# %%
df = pd.DataFrame({'DJ_GRCh38': [dj_count]})
X_pred = sm.add_constant(df[['DJ_GRCh38']], has_constant='add')
df['DJ_predicted'] = model.predict(X_pred)
print(f"Predicted DJ for GRCh38 with DJ count {dj_count}: {df['DJ_predicted'].iloc[0]:.2f}")


# %%
def calculate_DJ_prediction(df, fixed_means, pis, sigmas):
    df = df.copy()

    means_arr = np.asarray(fixed_means, dtype=float)
    K = len(means_arr)

    x = pd.to_numeric(df["DJ_predicted"], errors="coerce").to_numpy()
    un = np.vstack([
        pis[k] * norm.pdf(x, loc=means_arr[k], scale=sigmas[k])
        for k in range(K)
    ]).T
    probs = un / (un.sum(axis=1, keepdims=True) + 1e-300)

    group_idx = np.argmax(probs, axis=1).astype(int)

    df.loc[:, "group"] = group_idx
    df.loc[:, "group_mean"] = means_arr[group_idx]
    df.loc[:, "group_prob"] = probs[np.arange(len(df)), group_idx]

    group_map = dict(zip(range(len(means_arr)), np.round(means_arr)))
    df.loc[:, "group"] = df["group"].map(group_map)

    return df


# %%
fixed_means = [8,9.1,10.1,10.7,11.1,11.9,13]
pis = [0.00205486, 0.03653676, 0.89998205, 0.03648907, 0.02024702, 0.00344312, 0.00124713]
sigmas = [0.14599045, 0.18270757, 0.16448536, 0.12851647, 0.19678077, 0.19770631, 0.16025023]


# %%
df = calculate_DJ_prediction(df, fixed_means, pis, sigmas)
df.head()

# %%
x = np.linspace(7, 14, 1000)

# Compute each Gaussian component
components = []
for mu, sigma, pi in zip(fixed_means, sigmas, pis):
    component = pi * norm.pdf(x, mu, sigma)
    components.append(component)

# Total mixture
gmm = sum(components)

# Plot
plt.figure()
for component in components:
    plt.plot(x, component)

plt.axvline(df['DJ_predicted'].iloc[0], color='red', linestyle='--', label='Predicted DJ')
plt.legend()
plt.plot(x, gmm)
plt.xlabel("x")
plt.ylabel("Density")
plt.title("Gaussian Mixture Model")
plt.savefig(f"{plt_prefix}.pdf", dpi=300)
plt.show()

# %%
x = np.linspace(7, 14, 1000)

# Compute each Gaussian component
components = []
for mu, sigma, pi in zip(fixed_means, sigmas, pis):
    component = pi * norm.pdf(x, mu, sigma)
    components.append(component)

# Total mixture
gmm = sum(components)

# Plot
plt.figure()
for component in components:
    plt.plot(x, component)
plt.ylim(0, 0.2)
plt.axvline(df['DJ_predicted'].iloc[0], color='red', linestyle='--', label='Predicted DJ')
plt.legend()
plt.plot(x, gmm)
plt.xlabel("x")
plt.ylabel("Density")
plt.title("Gaussian Mixture Model")
plt.savefig(f"{plt_prefix}.zoomin.pdf", dpi=300)
plt.show()

# %%
