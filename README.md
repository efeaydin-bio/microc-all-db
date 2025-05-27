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

👉 [View the Reference Manual](https://data.cyverse.org/dav-anon/iplant/home/efeaydin/user_guide.pdf)

---

## 📄 Citation

If you use this database in your research, please cite the following publication:

👉 [Insert publication link here]

---

## 🚀 Usage

### 🌐 Web App (Online)

You can access the app online (no installation required) via Streamlit Cloud:

🔗 [https://microc-db-all.streamlit.app](https://microc-db-all.streamlit.app)

---

### 🐳 Docker Image (Offline)

For offline use and faster local performance, you can run the app via Docker. This option is fully self-contained and does **not** require internet once downloaded.

#### 🔧 Prerequisites

- Install Docker: [https://docs.docker.com/get-started/get-docker/](https://docs.docker.com/get-started/get-docker/)

#### 📦 Step-by-Step Instructions

1. **Download the appropriate Docker image:**

   - **Mac (Apple Silicon)**: [Insert link for Mac `.tar` file]
   - **Intel/AMD**: [Insert link for Intel `.tar` file]

2. **Load the image into Docker**:

   Open a terminal and run:

   docker load -i myFile.tar
   
4. ** Start the app
   docker run -p 8501:8501 microc-app

5. Access the app in your browser:
   Open http://localhost:8501

This app was developed to support chromatin interaction research in BCP ALL by providing an intuitive interface to query and visualize looping structures across different subtypes.

For questions, feedback, or collaboration inquiries, feel free to reach out.
