# Fading-Fidelity Model Predictive Control  
_Exploiting Sensitivity Decay for Computational Efficiency_

This repository contains code, data and experiments accompanying my master thesis.  
See the **`notebooks/`** folder to reproduce figures and interact with several examples of slow-fast systems.

## Setup 

### 0. Prerequisites – Git LFS (needed to pull the data files)

```bash
# macOS
brew install git-lfs
# Debian/Ubuntu
sudo apt-get install git-lfs
# Windows
choco install git-lfs

# Initialise the Git LFS hooks once:
git lfs install
```

### 1. Clone the repo

```bash
git clone git@github.com:l-m-schroth/fading_fidelity_MPC.git
cd fading_fidelity_MPC
```

### 2. (Optional) create a fresh environment

```bash
python -m venv venv          # or: python3 -m venv venv
source venv/bin/activate     # Windows: venv\Scripts\activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Run the setup script

```bash
pip install -e .
```

### 5. acados

Our code makes use of **acados**. If you already have an existing acados installation on your device, try running it.  
We encountered small issues with different multi-phase solver libraries wrongfully influencing each other. If that happens, install acados from our fork—these issues were addressed there: <https://github.com/l-m-schroth/acados>. In any case you need to install the python interface in the newly created venv as highlighted in the acados installation guide:

```bash
pip install -e <acados_root>/interfaces/acados_template
```