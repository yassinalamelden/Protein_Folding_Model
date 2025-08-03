# 🧬 Protein Folding using Genetic Algorithm & Simulated Annealing (HP Model)
This project explores protein folding using the **Hydrophobic-Polar (HP) model** on both **2D and 3D lattice structures**, with optimization via **Genetic Algorithm (GA)** and **Simulated Annealing (SA)**.  
Developed as part of my research internship at **Manipal University, Mahe**.

---

## 🚀 Features

- ✅ Folding simulation on:
  - 2D **Square** lattice
  - 2D **Triangular** lattice
  - (Planned) 3D Square lattice
- ✅ HP energy evaluation
- ✅ GA with selection, crossover, mutation, and elitism
- ✅ Visualization of folded protein using Streamlit
- ✅ Optimal energy comparison for known sequences
- 🧪 Simulated Annealing integration

---

## 🧱 Lattice Structures

| Lattice       | Description |
|---------------|-------------|
| Square (2D)   | Basic 4-directional model |
| Triangular (2D) | 6-directional with 60° turns |
| Square (3D)   | Future enhancement using 3D GA methods |

---

## 🖥️ Technologies Used

- Python 3
- Streamlit
- Matplotlib
- NumPy (for future 3D implementation)

---

## 📂 Project Structure

```bash
.
├── GA.py               # Main Streamlit app with GA logic
├── utils/              # (planned) Utility functions
├── results/            # Output visualizations or logs
├── sequences/          # Sample protein sequences
└── README.md
```

---

## 📌 Sample Use (Streamlit Interface)

1. Run the app:

```bash
streamlit run GA.py
```

2. Input your HP sequence (e.g., HPHPPHHPHPPHPHHPPHPH)
3. Select lattice type: Square or Triangular
4. Click "Run Simulation"

The app will:

- Run the GA folding algorithm
- Show minimum energy
- Compare with optimal (if available)
- Plot the folding structure

---

## 📉 Example Output

```bash
Minimum HP Energy: -9
Ideal Optimal Energy: -9
Deviation from Optimal: 0
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

## 🙋‍♂️ Author
Yassin Mahmoud \
Intern at Manipal University, Mahe \
[LinkedIn](https://www.linkedin.com/in/yassin-mahmoud-6130b5228)

