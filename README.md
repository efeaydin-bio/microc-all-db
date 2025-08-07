# Chromatin Interaction Database of BCP ALL

This web application presents a database of chromatin interactions identified through Micro-C experiments on **35 primary B-cell precursor acute lymphoblastic leukemia (BCP ALL)** samples.

Users can explore chromatin loops:
- Around specific **genomic coordinates**
- For the cis-regulatory elements of selected **genes of interest**
- Across **BCP ALL subtypes**
- In a **high-resolution merged map (1 kb resolution)**

---

## 📘 Reference Manual

A detailed guide to the app features, usage examples, and query types is available in the reference manual:

👉 [View the Reference Manual](https://lu.app.box.com/s/55ootvakcpdt0gzr6fytxl302acm44b0)

---

## 📄 Original Study


👉 https://www.biorxiv.org/content/10.1101/2025.07.22.666073v1

---

## 🛠️ Pre-requisites
### ✅ No Installation (Recommended)
The app can be accessed directly online via Streamlit Cloud — no installation or system configuration required.

### 💻 Run Locally (Advanced)

If you prefer to run the app on your own machine, it has been tested on **macOS** and requires the following dependencies:

#### 🔹 System Dependencies

These must be installed **before** installing Python packages:

- `bedtools`
- `zlib1g-dev`
- `build-essential`

🔹 Python packages

The following Python packages (with specific versions) are required:

matplotlib == 3.6.2  
pygenometracks == 3.9  
streamlit == 1.44.0  
pandas == 1.5.3  
numpy == 1.26.4  
pybedtools == 0.9.0  

## 🚀 Usage

### 🌐 Web App (Online)

You can access the app online (no installation required) via Streamlit Cloud:

🔗 [https://microc-db-all.streamlit.app](https://microc-db-all.streamlit.app)

---

### 🐳 Docker Image (Offline)

For offline use and faster local performance, a fully self-contained Docker Image will be uploaded here shortly.

---

### 💻 Manual Run (Offline)

Alternatively, you can clone the repository and run the app locally with:

```streamlit run app.py```
---




This app was developed to support chromatin interaction research in BCP ALL by providing an intuitive interface to query and visualize looping structures across different subtypes.

For questions, feedback, or collaboration inquiries, feel free to reach out.
