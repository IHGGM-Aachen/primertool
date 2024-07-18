# Primertool

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow)](https://opensource.org/licenses/MIT)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.36.0-ff4b4b?logo=streamlit)](https://docs.streamlit.io/)
[![Python](https://img.shields.io/badge/Python-3.8.18-3776AB?logo=python)](https://docs.python.org/release/3.8.18/)
[![Docker](https://img.shields.io/badge/Docker-27.0.3-2496ED?logo=docker)](https://docs.docker.com/engine/install/)

The premise of this package is to generate primers for PCR/Sanger sequencing for either:

- a specific variant ([hgvs](https://varnomen.hgvs.org/) variant nomenclature), either the whole exon or if not in an
  exon for the genomic position
- an exon (transcript number and exon )
- all exons of a transcript (transcript number)
- around a genomic position (chromosome and start/stop position)

This tool allows for primers based in hg19 or hg38.

---

## Table of Contents

- [Setup](#setup)
  - [Docker Installation (Deployment)](#docker-installation-deployment)
  - [Manual Installation (Development)](#manual-installation-development)
- [Documentation](#documentation)
- [License](#license)

## Setup
Clone the repository
    
    git clone https://github.com/IHGGM-Aachen/primertool.git
    
and navigate to the cloned repository.

### Docker Installation (Deployment)
Requires [Docker](https://docs.docker.com/engine/install/) installed.

1. **Build the Docker image**
    ```bash
    docker build -t primertool .
    ```
   
2. **Run the Docker container**
    ```bash
    docker run -p 8501:8501 primertool
    ```

### Manual Installation (Development)
1. **Install requirements**
    ```bash
    pip install requirements.txt
    ```
   
2. **Run the application**
    ```bash
    streamlit run /path/to/primertool/streamlit_main.py
    ```
    The application will now be running on `http://localhost:8501` and can be accessed through a web browser.    

## Documentation
The full documentation can be viewed here: [Primertool Documentation](https://IHGGM-Aachen.github.io/primertool/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
