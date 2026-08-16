# AI-Guided Phase Stability Prediction of Refractory Multi-Principal Element Alloys

This repository implements the complete machine learning pipeline used for phase stability prediction of refractory multi-principal element alloys (RMPEAs) using **Mixture Density Networks (MDNs)** and **uncertainty-guided active learning**.

The framework automates data generation, Thermo-Calc labeling, feature engineering, probabilistic model training, and active learning to efficiently explore large alloy design spaces.

---

## Features

- Mixture Density Networks (MDNs) for probabilistic phase prediction
- Bayesian optimization for hyperparameter tuning
- Thermo-Calc automated phase labeling
- Feature engineering and feature reduction
- Custom oversampling strategy for imbalanced phases
- Uncertainty-guided active learning
- Reproducible pipeline for manuscript figures

---

# Repository Structure

```
SampleGeneration/
PhaseLabeling/
FeatureLabeling/
FeatureEngineering/
Oversampling/
TrainingProcess/
PerformanceEvaluation/
ActiveLearning/
checkpoint/
```

---

# Workflow

```
Sample Generation
        │
        ▼
Phase Labeling (Thermo-Calc)
        │
        ▼
Feature Labeling
        │
        ▼
Feature Engineering
        │
        ▼
Oversampling
        │
        ▼
Bayesian Hyperparameter Optimization
        │
        ▼
MDN Training
        │
        ▼
Performance Evaluation
        │
        ▼
Active Learning
```

---

# Pipeline Description

## 1. Sample Generation

**Notebook**

`SampleGeneration.ipynb`

### Purpose

Generate the alloy compositions used throughout the framework.

### Output

- `Samples.txt`
- `Samples_T.txt`

The notebook

- generates 70,000 RMPEA compositions at **850 K**
- expands each composition to six temperature conditions
- produces approximately **420,000 training samples**

---

## 2. Phase Labeling

**Notebook**

`PhaseLabeling.ipynb`

### Purpose

Generate equilibrium phase labels using **Thermo-Calc**.

### Input

- Samples.txt
- Samples_T.txt

### Output

- Phases.txt
- Phases_T.txt

Generated labels include

- FCC
- BCC
- HCP
- B2
- Laves
- Sigma
- Heusler
- Liquid

---

## 3. Feature Labeling

**Notebook**

`FeatureLabeling.ipynb`

### Purpose

Generate the 51 material descriptors used for machine learning.

The notebook

- computes primary descriptors
- computes secondary descriptors
- produces

```
labels.txt
labels_extra.txt
```

> **Note**
>
> `Phi.txt` must first be generated using the publicly available ASAP descriptor package.

---

## 4. Feature Engineering

**Notebook**

`FeatureEngineering.ipynb`

### Stage 1 – Data Cleaning

- Initial screening
- Distribution balancing
- Dataset compilation

### Stage 2 – Feature Reduction

- Correlation analysis
- Removal of redundant descriptors
- Generation of final train/test datasets

Outputs include

```
xlo.npy
ylo.npy
xlo_test.npy
ylo_test.npy
```

---

## 5. Oversampling

**Notebook**

`Oversampling.ipynb`

### Purpose

Address class imbalance using the custom oversampling strategy developed in this work.

Produces

- oversampled training data
- validation datasets

---

## 6. Model Training

### Hyperparameter Optimization

Run

```
MDN_Finetuning.ipynb
```

Bayesian optimization searches for the optimal

- dropout
- network depth
- neurons
- batch size
- learning rate
- number of Gaussian mixtures

The best architecture is selected based on the minimum validation loss.

---

### Final Training

Run

```
MDN.ipynb
```

using

- oversampled training data
- validation data
- test data

The best checkpoint corresponds to the minimum validation loss.

---

## 7. Performance Evaluation

**Notebook**

`PerformanceEvaluation.ipynb`

Generates

- Figure 1
- Figures 3–7

The repository also includes trained checkpoints for

- six MDN models
- feature-ablation studies
- ensemble studies

---

## 8. Active Learning

**Notebook**

`ActiveLearning.ipynb`

Implements the uncertainty-guided active learning framework presented in the manuscript.

The notebook includes

- training dataset preparation
- active learning cycles
- uncertainty acquisition
- Figures 8–10

Example workflows for different sample sizes and active learning cycles are also included.

---

# Requirements

- Python
- PyTorch
- NumPy
- Scikit-learn
- Thermo-Calc
- ASAP Descriptor Package

---

# Citation

If you use this repository, please cite:

Shargh, A. K., Stiles, C. D., & El-Awady, J. A. (2026). Uncertainty-aware phase fraction prediction and active-learning-guided out-of-domain discovery of refractory multi-principal element alloys. arXiv preprint arXiv:2604.18322.

---

# Contact

Email: ashargh1@jhu.edu
