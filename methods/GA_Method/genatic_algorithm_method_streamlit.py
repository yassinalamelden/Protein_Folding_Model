# protein_folding_streamlit.py

import streamlit as st
import random
import heapq
import matplotlib.pyplot as plt
import math

# Directions for 2D lattice: Up, Right, Down, Left
square_directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]

# Directions for 2D Triangular lattice (60° angles)
triangular_directions = [
    (1.0, 0.0),
    (0.5, math.sqrt(3)/2),
    (-0.5, math.sqrt(3)/2),
    (-1.0, 0.0),
    (-0.5, -math.sqrt(3)/2),
    (0.5, -math.sqrt(3)/2)
]

optimal_2d = {
    20: -9, 24: -9, 25: -8, 36: -14,
    48: -23, 50: -21, 60: -36, 64: -42, 85: -53
}

def generate_random_fold(length):
    return [random.randint(0, 2) for _ in range(length - 1)]

def square_fold_protein(sequence, moves):
    direction = 0
    pos = (0, 0)
    coords = [pos]
    occupied = {pos}
    for move in moves:
        direction = (direction - 1) % 4 if move == 1 else (direction + 1) % 4 if move == 2 else direction
        dx, dy = square_directions[direction]
        pos = (pos[0] + dx, pos[1] + dy)
        if pos in occupied:
            return None
        coords.append(pos)
        occupied.add(pos)
    return coords

def triangular_fold_protein(sequence, moves):
    direction = 0
    pos = (0.0, 0.0)
    coords = [pos]
    occupied = {pos}
    for move in moves:
        direction = (direction + 1) % 6 if move == 1 else (direction - 1) % 6 if move == 2 else direction
        dx, dy = triangular_directions[direction]
        new_pos = (round(pos[0] + dx, 5), round(pos[1] + dy, 5))
        if new_pos in occupied:
            return None
        coords.append(new_pos)
        occupied.add(new_pos)
        pos = new_pos
    return coords

def true_energy(sequence, coords):
    if not coords or len(coords) != len(sequence): return 0
    contacts = set()
    for i, s in enumerate(sequence):
        if s != 'H': continue
        x, y = coords[i]
        for dx, dy in square_directions:
            neighbor = (x + dx, y + dy)
            if neighbor in coords:
                j = coords.index(neighbor)
                if j < len(sequence) and sequence[j] == 'H' and abs(i - j) > 2:
                    contacts.add(tuple(sorted((i, j))))
    return -len(contacts)

def compactness_penalty(coords):
    if not coords or len(coords) < 2: return 0
    xs, ys = zip(*coords)
    area = (max(xs) - min(xs) + 1) * (max(ys) - min(ys) + 1)
    return 1 / area if area > 0 else 0

def multi_point_crossover(p1, p2):
    if len(p1) < 3: return p1
    cut1, cut2 = sorted(random.sample(range(1, len(p1)), 2))
    return p1[:cut1] + p2[cut1:cut2] + p1[cut2:]

def mutate(moves, rate):
    return [random.randint(0, 2) if random.random() < rate else m for m in moves]

def tournament_selection(scored, size=3, winners=100):
    selected = []
    actual_size = min(size, len(scored))
    if actual_size < 1: return selected
    while len(selected) < winners and len(scored) >= actual_size:
        group = random.sample(scored, actual_size)
        best = min(group, key=lambda x: x[0])[1]
        selected.append(best)
    return selected

def remove_twins(population):
    return list({tuple(m): m for m in population}.values())

def genetic_algorithm(sequence, lattice, generations=5000, pop_size=200):
    best_energy = float('inf')
    best_coords = None
    population = [generate_random_fold(len(sequence)) for _ in range(pop_size)]
    fold_func = square_fold_protein if lattice == "Square" else triangular_fold_protein
    for gen in range(generations):
        rate = max(0.02, 0.2 * (1 - gen / generations))
        scored = []
        for moves in population:
            coords = fold_func(sequence, moves)
            if coords and len(coords) == len(sequence):
                energy = true_energy(sequence, coords)
                penalty = compactness_penalty(coords)
                score = energy - penalty
                if energy < best_energy:
                    best_energy = energy
                    best_coords = coords
                scored.append((score, moves))
        if len(scored) < 2:
            population = [generate_random_fold(len(sequence)) for _ in range(pop_size)]
            continue
        top_folds = heapq.nsmallest(100, scored)
        elite = [top_folds[0][1]]
        parents = tournament_selection(top_folds)
        new_population = elite.copy()
        while len(new_population) < pop_size and len(parents) >= 2:
            p1, p2 = random.sample(parents, 2)
            child = mutate(multi_point_crossover(p1, p2), rate)
            new_population.append(child)
        population = remove_twins(new_population)
    return best_energy, best_coords

def visualize(sequence, coords):
    fig, ax = plt.subplots(figsize=(6, 6))
    xs, ys = zip(*coords)
    colors = ['red' if c == 'H' else 'blue' for c in sequence]
    for i in range(1, len(coords)):
        x1, y1 = coords[i-1]
        x2, y2 = coords[i]
        ax.plot([x1, x2], [y1, y2], 'k-')
    for (x, y), color in zip(coords, colors):
        ax.plot(x, y, 'o', markersize=10, color=color)
    ax.set_title("Protein Fold Visualization (Red=H, Blue=P)")
    ax.set_aspect('equal')
    ax.grid(True)
    st.pyplot(fig)

# 🧬 Streamlit UI
st.title("Protein Folding Simulator")
st.markdown("Enter your HP sequence and choose lattice type to simulate folding using a genetic algorithm.")

sequence_input = st.text_input("Protein sequence (only H and P):", max_chars=200).upper().replace(" ", "")
lattice_choice = st.radio("Lattice type:", ["Square", "Triangular"])
run_button = st.button("Run Simulation")

if run_button and sequence_input:
    if not all(c in "HP" for c in sequence_input):
        st.error("❌ Sequence must contain only 'H' and 'P' characters.")
    else:
        with st.spinner("Running simulation..."):
            energy, coords = genetic_algorithm(sequence_input, lattice_choice)
            optimal = optimal_2d.get(len(sequence_input), None)
            deviation = energy - optimal if optimal is not None else "N/A"
            st.subheader("🧪 Simulation Results")
            st.text(f"Minimum HP Energy: {energy}")
            st.text(f"Ideal Optimal Energy: {optimal if optimal else 'Unknown'}")
            st.text(f"Deviation from Optimal: {deviation}")
            if coords:
                visualize(sequence_input, coords)
