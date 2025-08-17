# 🧬 Protein Folding using Genetic Algorithm & Simulated Annealing (HP Model)
This project explores protein folding using the **Hydrophobic-Polar (HP) model** on both **2D and 3D lattice structures**, with optimization via **Genetic Algorithm (GA)** and **Simulated Annealing (SA)**.  
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
│── methods/                # Implementation of GA, SA, PSO (in progress)
│── benchmarks/             # Reference sequences & optimal energies
│── utils/                  # Helper functions (FASTA parser, visualization, etc.)
│── streamlit_app/          # Streamlit-based interactive UI
│── results/                # Output energies, runtimes, visualizations
│── requirements.txt        # Python dependencies
│── README.md               # Project documentation
```

---

## 📌 Usage

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

---

## 📉 Example Output

For Genatic Algoritm (GA):

```bash
# Run the Streamlit app
streamlit run streamlit_app/genatic_algorithm_method_streamlit.py
```

input
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

## 📊 Optimization Techniques
✅ Genetic Algorithm

- Random population initialization
- Tournament selection
- Multi-point crossover
- Mutation with decaying rate
- Compactness penalty to favor compact folds

🔄 Simulated Annealing

- Random fold mutation
- Temperature-controlled probabilistic acceptance
- Iterative convergence to energy minimum

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

## 📄 License
This Project is an open source and available under the MIT License.

---

## 📬 Contact
Author: Yassin Mahmoud \
Intern at Manipal University, Mahe, India \
[Gmail](yassin.alamelden@gmail.com) \
[LinkedIn](https://www.linkedin.com/in/yassin-mahmoud-6130b5228)

