# 🧬 Protein Folding using Genetic Algorithm, Simulated Annealing & Particle Swarm Optimization (HP Model)
This project explores protein folding using the **Hydrophobic-Polar (HP) model** on both **2D and 3D lattice structures**, with optimization via **Genetic Algorithm (GA)**, **Simulated Annealing (SA)** and **Particle Swarm Optimization**.  
Developed as part of my research internship at **Manipal University, Mahe, India**.

---

## 🚀 Features

- Lattice Models:
  - 2D Square lattice
  - 2D Triangular lattice
  - 3D Cubic lattice (Genetic Algorithm, experimental)
- Optimization Methods:
  - Genetic Algorithm (GA) with selection, crossover, mutation, elitism
  - Simulated Annealing (SA) with temperature-based acceptance
  - Particle Swarm Optimization (PSO) with 
- New Additions:

✅ FASTA-to-HP sequence conversion \
✅ Benchmark comparison with optimal 2D energies \
✅ Streamlit-based interactive interface for folding visualization and metrics \
✅ Sequence energy & runtime logging \
✅ Support for benchmark sequences dataset

---

## 🧱 Lattice Structures

| Lattice         | Description |
|-----------------|-------------|
| Square (2D)     | Basic 4-directional model |
| Triangular (2D) | 6-directional with 60° turns |
| Cubic (3D)      | Future enhancement using 3D GA methods |

---

## 🖥️ Technologies Used

- Python 3
- Streamlit
- Matplotlib
- NumPy (for future 3D implementation)

---

## 📂 Project Structure

```bash
Protein_Folding_Model/
│── methods/
     └── GA_Method
     └── SA_Method
     └── PSO_Method          
│── data/                   # Benchmark HP sequences & known optimal energies
     └── utils/             # FASTA parser, visualization, helper functions
│── results/                # Output energies, runtimes, visualizations
│── docs/                   # Documentation
│    └── research_paper.pdf # Research paper
│── requirements.txt        # Python dependencies
│── README.md               # Project documentation
```

---

## 📌 Usage

### **Option 1 – Run Locally**

1. Clone the repository:

```bash
git clone https://github.com/yassinalamelden/Protein_Folding_Model.git
cd Protein_Folding_Model
```

2. Install dependencies

```bash
pip install -r requirements.txt
```

3.Run Streamlit interface

```bash
streamlit run streamlit_app/genatic_algorithm_method_streamlit.py
```

4.Select your options
- Input HP sequence or upload FASTA file
- Choose lattice type (2D square, 2D triangular, or 3D cubic)
- Select optimization method (GA or SA)
- Run simulation and view:
The app will:
  - Run the GA folding algorithm
  - Show minimum energy
  - Compare with optimal (if available)
  - Plot the folding structure
  - Runtime statistics

### **Option 2 – Run on the Cloud**

[**Launch Streamlit App**](https://proteinfoldingmodel.streamlit.app/)

(No setup required — just paste your HP sequence or upload a FASTA file to get folding results.)

---

## 📉 Example Output

For Genatic Algoritm (GA):

```bash
# Run the Streamlit app
streamlit run streamlit_app/genatic_algorithm_method_streamlit.py
```

Input
- Sequence: HPHPPHHPHPPHPHHPPHPH (length = 20)
- Lattice: 2D Square
- Method: Genetic Algorithm (GA)

Output

✅ Minimum Energy Found: -8 \
📊 Known Optimal Energy: -9 \
⏱ Runtime: ~2.4 seconds \
🖼 Visualization: Fold plotted on 2D lattice with hydrophobic (H) and polar (P) residues \

For Simulated Annealing (SA):

```bash
# Run SA implementation
python methods/SA_Method/simulated_annealing_method.py
```

Output Example:

```bash
Sequence length: 24
Lattice: 2D Triangular
Minimum Energy Found: -9
Known Optimal Energy: -9
Deviation: 0
Runtime: 3.1s
```

---

## 🔎 Optimization Techniques

We have implemented and compared **three metaheuristic algorithms** for protein folding on the HP model:

- 🧬 **Genetic Algorithm (GA)**  
  - Population-based search with crossover, mutation, elitism  
  - Strong scalability for both 2D and 3D lattices  
  - Consistently achieves near-optimal energies  

- 🔥 **Simulated Annealing (SA)**  
  - Single-solution probabilistic search with temperature-based acceptance  
  - Useful for escaping local minima  
  - Performance depends heavily on cooling schedule  

- 🐦 **Particle Swarm Optimization (PSO)**  
  - Swarm-based search guided by personal best (`pBest`) and global best (`gBest`)  
  - Encodes folding as particle positions/velocities in lattice space  
  - Achieves stable folds but slower convergence than GA on longer sequences  

---

### ✅ Summary
- **GA**: Best trade-off between runtime, accuracy, and scalability  
- **SA**: Flexible, avoids local minima but sensitive to parameters  
- **PSO**: Effective but slower on discrete HP lattice compared to GA  

---

## 📈 Optimal Energies Reference (2D HP Sequences)

| Sequence Length | Known Optimal Energy |
| --------------- | -------------------- |
| 20              | -9                   |
| 36              | -14                  |
| 48              | -23                  |
| 60              | -36                  |
| 85              | -53                  |

---

## 📄 Research Paper
This repository is supported by the following research work:

Optimization-Based Prediction of Protein Structures in the 2D & 3D HP Model: \
A Comparative Study of Metaheuristic Approaches \
Authors: Yassin M. Alam Elden, Omar M. Hegab, Mariam E. Elshamy

---

## 📜 License
This Project is an open source and available under the MIT License.

---

## 📬 Contact
Author: Yassin Mahmoud \
Intern at Manipal University, Mahe, India \
[Gmail](yassin.alamelden@gmail.com) \
[LinkedIn](https://www.linkedin.com/in/yassin-mahmoud-6130b5228)

