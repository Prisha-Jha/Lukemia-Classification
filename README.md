# Leukemia Subtype Classification via Singular Value Decomposition (SVD)

An end-to-end computational linear algebra pipeline that classifies leukemia subtypes—**Acute Lymphoblastic Leukemia (ALL)** versus **Acute Myeloid Leukemia (AML)**—using Singular Value Decomposition (SVD), orthogonal subspace projection, and geometric nearest-centroid classification on high-dimensional microarray gene expression data (7,129 genes).

---

## 📌 Project Overview

Gene expression datasets are characterized by the **curse of dimensionality** ($p \gg n$), containing thousands of gene measurements for a relatively small cohort of patients. Analyzing these datasets directly introduces significant biological noise and collinearity.

This project implements an unsupervised linear algebra pipeline to:
1. Reduce dimensionality from **7,129 gene features down to a 2D diagnostic subspace** ($\mathbb{R}^2$) without forming the explicit $7129 \times 7129$ covariance matrix.
2. Uncover the primary **biomarker genes** responsible for differentiating ALL from AML.
3. Build an **interactive nearest-centroid diagnostic engine** using Euclidean distance ($L_2$ norm) to classify unseen patient samples.

---

## 🧬 Dataset

* **Source:** Golub et al. (1999) Leukemia Benchmark Dataset (accessed via Kaggle / OpenML).
* **Dimensionality:** 7,129 gene expression levels across clinical patient tissue samples.
* **Target Classes:** 
  * `ALL`: Acute Lymphoblastic Leukemia
  * `AML`: Acute Myeloid Leukemia

---

## 📐 Mathematical Pipeline
Raw Data Matrix X (N × 7129)
│
▼
Mean-Centering (X_c = X - μ)
│
▼
Singular Value Decomposition: X_c = U · Σ · V^T
│
▼
Variance Analysis & Truncation (k = 2) via Scree Spectrum
│
▼
Orthogonal Subspace Projection: Z = X_c · V_k  (N × 2)
│
├──► Biomarker Discovery (Ranking |v_1| weights)
│
└──► Nearest-Centroid Classifier (L2 Distance in 2D Subspace)

### 1. Mean-Centering
Translates the origin to the dataset's center of mass:
$$\vec{\mu} = \frac{1}{N} \sum_{i=1}^N \vec{x}_i, \quad X_c = X - \vec{\mu}$$
This isolates disease-specific variations from shared baseline cellular metabolism.

### 2. Singular Value Decomposition (SVD)
Factors the centered matrix directly:
$$X_c = U \Sigma V^T$$
* **$U$ ($N \times N$):** Left singular vectors representing patient coordinates along orthonormal axes.
* **$\Sigma$ ($N$ singular values):** Quantifies variance/spread along each component ($\sigma_1 \ge \sigma_2 \ge \dots$).
* **$V^T$ ($N \times 7129$):** Right singular vectors representing principal gene loading axes (metagenes).

### 3. Dimensionality Reduction & Subspace Projection
Extracts the top $k = 2$ principal basis vectors ($V_k \in \mathbb{R}^{7129 \times 2}$) and projects the centered patient profiles:
$$Z = X_c V_k \in \mathbb{R}^{N \times 2}$$

### 4. Biomarker Discovery
The first right singular vector ($\vec{v}_1^T$, Row 0 of $V^T$) serves as the primary disease-separating axis. Evaluating absolute magnitudes:
$$\text{Impact}_j = \vert{}v_{1, j}\vert{}$$
ranks and isolates the highest-contributing biomarker genes.

### 5. Nearest-Centroid Classifier ($L_2$ Norm)
Computes cluster centroids in the 2D subspace:
$$\vec{\mu}_{\text{ALL}} = \frac{1}{N_{\text{ALL}}} \sum_{i \in \text{ALL}} \vec{z}_i, \quad \vec{\mu}_{\text{AML}} = \frac{1}{N_{\text{AML}}} \sum_{j \in \text{AML}} \vec{z}_j$$

For any patient sample $\vec{x}_{\text{new}}$:
$$\vec{z}_{\text{new}} = (\vec{x}_{\text{new}} - \vec{\mu}) V_k$$
$$\hat{y} = \arg\min_{c \in \{\text{ALL}, \text{AML}\}} \Vert{}\vec{z}_{\text{new}} - \vec{\mu}_c\Vert{}_2$$

---
