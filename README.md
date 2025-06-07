# baez-duarte-kalman
Code for Kalman-filtered Báez–Duarte approximants toward Riemann Hypothesis 
# Báez–Duarte Kalman-Filtered Approximation

This repository contains Python code implementing a Kalman-filtered Báez–Duarte approximation to the constant function 1, in the context of the Nyman–Beurling reformulation of the Riemann Hypothesis.

## 🔍 Overview

The code constructs a linear combination of modular functions r_k(n) = n \bmod k to approximate the constant sequence y_n = 1, applies a smoothing low-pass filter, and finally applies a scalar Kalman filter to enforce convergence.

The project supports the analytic work described in our accompanying arXiv paper:

> “Analytic Validation of the Kalman-Filtered Báez–Duarte Approximation in L^2(0,1)”  
> [[Link to arXiv will go here]]

## 🧠 Mathematical Framework

This work is based on the Báez–Duarte theorem, which states that:

> The Riemann Hypothesis is true if and only if the constant function 1 lies in the closure of the span of  
> r_k(x) := \left\{ \frac{x}{k} \right\}, \quad k \in \mathbb{N}  
> in the Hilbert space L^2(0,1).

We implement this approximation numerically using weighted least squares, moving average smoothing, and Kalman filtering.

## 📂 Files

- `Baez_Duarte_Kalman_Filtered.py` — Main Python script with all computation and plots
- `README.md` — This file
- `.gitignore` — Standard exclusions for Python projects

## 📈 Output

The script prints:
- Raw L² error before filtering
- Error after moving average smoothing
- Kalman-filtered error (typically ≈ 0)

It also generates:
- A plot of the approximants vs. the constant function
- A log-scale convergence plot of L² error vs. basis size M

## 🔗 Citation

If you use this code, please cite the arXiv paper or this GitHub repository directly.

## 📜 License

This project is open-source under the MIT License.
