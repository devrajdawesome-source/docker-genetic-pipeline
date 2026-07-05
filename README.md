# Dockerized Genetic Analysis Pipeline (`ldetect`)

A production-ready containerization blueprint engineered to isolate, compile, and deploy the legacy `ldetect` bioinformatics framework. This configuration eliminates cross-platform environment drift and mitigates OS-level library conflicts across research deployments.

## ⚙️ Bundled Toolchains & Runtimes
* **Languages:** Python (Analytical Kernels), R-Base (Statistical Run-time Layers)
* **Bioinformatics Packages:** `fast1` (Sequence processing), `plink` (Whole-genome association analysis), `bcftools` (Variant calling data manipulation)

## 📦 Deployment Instructions

### 1. Build the Isolated Image
To compile the multi-stage environment blueprint locally, execute:
```bash
docker build -t genetic-pipeline .
