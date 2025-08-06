# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 03:45:17 2025

@author: Yassin, Omar, Mariam
"""

import random
import heapq
import matplotlib.pyplot as plt
import math
import streamlit as st
import re
import sys
from mpl_toolkits.mplot3d import Axes3D  # for 3D plotting
import time

# Directions for 2D lattice: Up, Right, Down, Left
square_directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]

## 3D Cubic Directions
cubic_directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

# Directions for 2D Trianguler Lattice:
trianguler_directions = [
    (1.0, 0.0),  # 0: Right
    (0.5, math.sqrt(3) / 2),  # 1: Up-Right
    (-0.5, math.sqrt(3) / 2),  # 2: Up-Left
    (-1.0, 0.0),  # 3: Left
    (-0.5, -math.sqrt(3) / 2),  # 4: Down-Left
    (0.5, -math.sqrt(3) / 2)  # 5: Down-Right
]

### Optimal Energies (Benchmarks dataset)
optimal_2d = {
    20: -9,
    24: -9,
    25: -8,
    36: -14,
    48: -23,
    50: -21,
    60: -36,
    64: -42,
    85: -53
}

optimal_3d = {
    20: -11,
    24: -13,
    25: -9,
    36: -18,
    48: -31,
    50: -34,
    60: -55,
    64: -59
}


### Converting FASTA to HP Sequence
def fasta_to_hp(fasta_sequence):
    """
    Converts a FASTA amino acid sequence to HP model using:
    H = Hydrophobic {A, G, I, L, M, F, P, W, V}
    P = Polar       {R, N, D, C, E, Q, H, K, S, T, Y}
    """

    h_residues = {'A', 'G', 'I', 'L', 'M', 'F', 'P', 'W', 'V'}
    p_residues = {'R', 'N', 'D', 'C', 'E', 'Q', 'H', 'K', 'S', 'T', 'Y'}

    sequence = fasta_sequence.upper()

    hp_sequence = ''

    # aa : amino acid
    for aa in sequence:
        if aa in h_residues:
            hp_sequence += 'H'
        elif aa in p_residues:
            hp_sequence += 'P'
        else:
            raise ValueError(f"Unkown amino acid '{aa}' in sequence. ")
    return hp_sequence


# Randomizing the fold
def generate_random_fold(length):
    return [random.randint(0, 2) for _ in range(length - 1)]


# Protein Folding: 2D Square
def square_fold_protein(sequence, moves):
    square_direction = 0
    square_position = (0, 0)
    square_coords = [square_position]
    square_occupied = {square_position}
    for move in moves:
        square_direction = (square_direction - 1) % 4 if move == 1 else (square_direction + 1) % 4 if move == 2 else square_direction
        dx, dy = square_directions[square_direction]
        square_position = (square_position[0] + dx, square_position[1] + dy)
        if square_position in square_occupied:
            return None
        square_coords.append(square_position)
        square_occupied.add(square_position)
    return square_coords


## Protein Folding: 3D Cubic Lattice
def rotate_yaw_left(d): x, y, z = d; return (-z, y, x)


def rotate_yaw_right(d): x, y, z = d; return (z, y, -x)


def rotate_pitch_up(d): x, y, z = d; return (x, -z, y)


def rotate_pitch_down(d): x, y, z = d; return (x, z, -y)


def reverse_direction(d): return tuple(-v for v in d)


def fold_protein_3d_full(seq, moves):
    pos, coords, occupied = (0, 0, 0), [(0, 0, 0)], {(0, 0, 0)}
    direction = (1, 0, 0)
    for move in moves:
        if move == 1:
            direction = rotate_yaw_left(direction)
        elif move == 2:
            direction = rotate_yaw_right(direction)
        elif move == 3:
            direction = rotate_pitch_up(direction)
        elif move == 4:
            direction = rotate_pitch_down(direction)
        elif move == 5:
            direction = reverse_direction(direction)
        dx, dy, dz = direction
        pos = (pos[0] + dx, pos[1] + dy, pos[2] + dz)
        if pos in occupied: return None
        coords.append(pos)
        occupied.add(pos)
    return coords


# Protein Folding: 2D Trianguler
def trianguler_fold_protein(sequence, moves):
    trianguler_direction = 0  # Start facing right (index 0)
    trianguler_position = (0.0, 0.0)
    trianguler_coords = [trianguler_position]
    trianguler_occupied = {trianguler_position}

    for move in moves:
        if move == 1:
            trianguler_direction = (trianguler_direction + 1) % 6  # Left turn (60°)
        elif move == 2:
            trianguler_direction = (trianguler_direction - 1) % 6  # Right turn (-60°)

        dx, dy = trianguler_directions[trianguler_direction]
        new_position = (round(trianguler_position[0] + dx, 5), round(trianguler_position[1] + dy, 5))

        if new_position in trianguler_occupied:
            return None

        trianguler_coords.append(new_position)
        trianguler_occupied.add(new_position)
        trianguler_position = new_position

    return trianguler_coords


### Energy Functions
# Actual Energy for 2D Square Lattice
def true_energy(sequence, coords):
    if not coords or len(coords) != len(sequence):
        return 0
    contacts = set()
    for i in range(len(sequence)):
        if sequence[i] != 'H':
            continue
        x, y = coords[i]
        for dx, dy in square_directions:
            neighbor = (x + dx, y + dy)
            if neighbor in coords:
                j = coords.index(neighbor)
                if j < len(sequence) and sequence[j] == 'H' and abs(i - j) > 2:
                    contacts.add(tuple(sorted((i, j))))
    return -len(contacts)


## Actual Energy for 3D Cubic Lattice
def true_3d_energy(seq, coords):
    if not coords or len(coords) != len(seq): return 0
    contacts = set()
    for i, aa in enumerate(seq):
        if aa != 'H': continue
        x, y, z = coords[i]
        for dx, dy, dz in cubic_directions:
            n = (x + dx, y + dy, z + dz)
            if n in coords:
                j = coords.index(n)
                if seq[j] == 'H' and abs(i - j) > 2:
                    contacts.add(tuple(sorted((i, j))))
    return -len(contacts)


### Compactness Penalties
def compactness_penalty(coords):
    if not coords or len(coords) < 2:
        return 0
    xs, ys = zip(*coords)
    area = (max(xs) - min(xs) + 1) * (max(ys) - min(ys) + 1)
    return 1 / area if area > 0 else 0


def compactness_penalty_3d(coords):
    if not coords or len(coords) < 2: return 0
    xs, ys, zs = zip(*coords)
    volume = (max(xs) - min(xs) + 1) * (max(ys) - min(ys) + 1) * (max(zs) - min(zs) + 1)
    return 1 / volume if volume > 0 else 0


### Genetic Operators

### The multi Point Crossover for both 2D and 3D
def multi_point_crossover(p1, p2):
    if len(p1) < 3:
        return p1
    cut1, cut2 = sorted(random.sample(range(1, len(p1)), 2))
    return p1[:cut1] + p2[cut1:cut2] + p1[cut2:]


# Mutation
def mutate(moves, rate):
    return [random.randint(0, 2) if random.random() < rate else m for m in moves]


## Mutation
def mutate_3d(moves, rate):
    return [random.randint(0, 5) if random.random() < rate else m for m in moves]


def tournament_selection(scored, size=3, winners=100):
    selected = []
    actual_size = min(size, len(scored))
    if actual_size < 1:
        return selected
    while len(selected) < winners and len(scored) >= actual_size:
        group = random.sample(scored, actual_size)
        best = min(group, key=lambda x: x[0])[1]
        selected.append(best)
    return selected


def remove_twins(population):
    return list({tuple(m): m for m in population}.values())


## New: Generate a valid folding only
def generate_valid_fold(length, max_attempts=1000):
    for _ in range(max_attempts):
        moves = [random.randint(0, 5) for _ in range(length - 1)]
        coords = fold_protein_3d_full("H" * length, moves)
        if coords and len(coords) == length:
            return moves
    return None


### Genetic Algorithm Implementations

# GA Main Block for Square 2D
def square_genetic_algorithm(sequence, generations=5000, pop_size=200):
    square_best_energy = float('inf')
    square_best_coords = None
    population = [generate_random_fold(len(sequence)) for _ in range(pop_size)]

    for gen in range(generations):
        mutation_rate = max(0.02, 0.2 * (1 - gen / generations))
        square_scored = []

        for moves in population:
            square_coords = square_fold_protein(sequence, moves)
            if square_coords and len(square_coords) == len(sequence):
                square_energy = true_energy(sequence, square_coords)
                penalty = compactness_penalty(square_coords)
                square_score = square_energy - penalty
                if square_energy < square_best_energy:
                    square_best_energy = square_energy
                    square_best_coords = square_coords
                square_scored.append((square_score, moves))

        if len(square_scored) < 2:
            population = [generate_random_fold(len(sequence)) for _ in range(pop_size)]
            continue

        top_folds = heapq.nsmallest(100, square_scored)
        elite = [top_folds[0][1]]
        parents = tournament_selection(top_folds)

        new_population = elite.copy()
        while len(new_population) < pop_size and len(parents) >= 2:
            p1, p2 = random.sample(parents, 2)
            child = mutate(multi_point_crossover(p1, p2), mutation_rate)
            new_population.append(child)

        population = remove_twins(new_population)

    return square_best_energy, square_best_coords


# GA Main Block for Trianguler 2D lattice
def trianguler_genetic_algorithm(sequence, generations=5000, pop_size=500):
    best_energy = float('inf')
    best_coords = None
    population = [generate_random_fold(len(sequence)) for _ in range(pop_size)]

    for gen in range(generations):
        mutation_rate = max(0.02, 0.2 * (1 - gen / generations))
        scored = []

        for moves in population:
            trianguler_coords = trianguler_fold_protein(sequence, moves)
            if trianguler_coords and len(trianguler_coords) == len(sequence):
                energy = true_energy(sequence, trianguler_coords)
                penalty = compactness_penalty(trianguler_coords)
                score = energy - penalty
                if energy < best_energy:
                    best_energy = energy
                    best_coords = trianguler_coords
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
            child = mutate(multi_point_crossover(p1, p2), mutation_rate)
            new_population.append(child)

        population = remove_twins(new_population)

    return best_energy, best_coords


## GA Main Block for 3D Cubic Lattice
def genetic_algorithm_3d(sequence, generations=5000, pop_size=200):
    best_energy = float('inf')
    best_coords = None
    population = []
    while len(population) < pop_size:
        move = generate_valid_fold(len(sequence))
        if move: population.append(move)
    stale_count = 0

    for gen in range(generations):
        mutation_rate = max(0.02, 0.2 * (1 - gen / generations))
        scored = []
        valid_count = 0

        for moves in population:
            coords = fold_protein_3d_full(sequence, moves)
            if coords and len(coords) == len(sequence):
                valid_count += 1
                energy = true_3d_energy(sequence, coords)
                penalty = 0.3 * compactness_penalty_3d(coords)
                score = energy - penalty
                if energy < best_energy:
                    best_energy = energy
                    best_coords = coords
                    stale_count = 0
                scored.append((score, moves))

        if len(scored) < 2:
            population = []
            while len(population) < pop_size:
                move = generate_valid_fold(len(sequence))
                if move: population.append(move)
            continue

        if stale_count >= 200:
            population = []
            while len(population) < pop_size:
                move = generate_valid_fold(len(sequence))
                if move: population.append(move)
            stale_count = 0
            continue

        top = heapq.nsmallest(100, scored)
        elite = [x[1] for x in top[:5]]
        parents = tournament_selection(top, winners=pop_size // 2)
        new_pop = elite[:]
        while len(new_pop) < pop_size and len(parents) >= 2:
            p1, p2 = random.sample(parents, 2)
            child = mutate_3d(multi_point_crossover(p1, p2), mutation_rate)
            coords = fold_protein_3d_full(sequence, child)
            if coords and len(coords) == len(sequence):
                new_pop.append(child)
        population = remove_twins(new_pop)
        stale_count += 1

    return best_energy, best_coords


### Visualization of the Fold (2D and 3D)
def visualize_fold(sequence, coords, min_energy):
    if not coords or len(coords) != len(sequence):
        st.error("Cannot visualize: invalid folding")
        return
    dim = len(coords[0])
    colors = ['red' if c == 'H' else 'blue' for c in sequence]

    if dim == 2:
        xs, ys = zip(*coords)
        fig, ax = plt.subplots(figsize=(8, 8))
        for i in range(1, len(coords)):
            x1, y1 = coords[i - 1]
            x2, y2 = coords[i]
            ax.plot([x1, x2], [y1, y2], 'k-')
        for (x, y), color in zip(coords, colors):
            ax.plot(x, y, 'o', markersize=12, color=color)
        ax.set_title(f"2D Protein Fold Visualization (Red=H, Blue=P)\nMinimum Energy = {min_energy}", fontsize=18)
        ax.set_aspect('equal')
        ax.grid(True)
        st.pyplot(fig)

    elif dim == 3:
        from mpl_toolkits.mplot3d import Axes3D  # Just in case

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        xs, ys, zs = zip(*coords)
        for i in range(1, len(coords)):
            x1, y1, z1 = coords[i - 1]
            x2, y2, z2 = coords[i]
            ax.plot([x1, x2], [y1, y2], [z1, z2], 'k-')
        for (x, y, z), color in zip(coords, colors):
            ax.scatter(x, y, z, color=color, s=100)
        ax.set_title(f"3D Protein Fold Visualization (Red=H, Blue=P)\nMinimum Energy = {min_energy}", fontsize=16)
        ax.set_xlabel("X", fontweight='bold')
        ax.set_ylabel("Y", fontweight='bold')
        ax.set_zlabel("Z", fontweight='bold')
        ax.grid(True)
        st.pyplot(fig)


####---- Main Program ----####


def clean_sequence(seq):
    # Keeps only H and P characters, discards everything else
    return ''.join(re.findall(r'[HP]', seq.upper()))

def main():
    sequence = None
    st.title("Protein Folding Genetic Algorithm")
    input_type = st.radio("Select Input Type:", ["FASTA Sequence", "HP Sequence"])

    if input_type == "FASTA Sequence":
        fasta_sequence = st.text_input("Enter the FASTA sequence:")
        if fasta_sequence:
            try:
                sequence = fasta_to_hp(fasta_sequence)
                sequence = sequence.replace(" ", "")  # Clean any spaces
                st.success(f"Converted HP Sequence: {sequence}")
            except ValueError as e:
                st.error(f"Error: {e}")
                return
    else:
        import re

        raw_sequence = st.text_input("Enter HP sequence (only H and P):")
        if raw_sequence:
            cleaned_sequence = clean_sequence(raw_sequence)
            if not cleaned_sequence:
                st.error("Sequence is empty or invalid. Please enter only 'H' and 'P'.")
                return
            sequence = cleaned_sequence

    if not sequence:
        return

    structure_choice = st.selectbox("Choose Structure Type:", ["2D", "3D"])
    length = len(sequence)
    st.write(f"Sequence Length: {length}")

    if structure_choice == "2D":
        lattice_choice = st.radio("Select 2D Lattice Type:", ["Square Lattice", "Triangular Lattice"])
        if st.button("Run 2D GA"):
            if lattice_choice == "Square Lattice":
                st.write("Running Square 2D Lattice...")
                start_time = time.time()
                min_energy, best_coords = square_genetic_algorithm(sequence)
                end_time = time.time()
                elapsed_time = end_time - start_time
                st.write(f"⏱️ Running Time: {elapsed_time:.2f} seconds")
            else:
                st.write("Running Triangular 2D Lattice...")
                start_time = time.time()
                min_energy, best_coords = trianguler_genetic_algorithm(sequence)
                end_time = time.time()
                elapsed_time = end_time - start_time
                st.write(f"⏱️ Running Time: {elapsed_time:.2f} seconds")
            optimal_energy = optimal_2d.get(length, "N/A")
            st.success(f"Minimum HP Energy: {min_energy}")
            st.info(f"Optimal Energy: {optimal_energy}")

            visualize_fold(sequence, best_coords, min_energy)

    elif structure_choice == "3D":
        if st.button("Run 3D GA"):
            st.write("Running 3D Genetic Algorithm...")
            start_time = time.time()
            min_energy, best_coords = genetic_algorithm_3d(sequence)
            end_time = time.time()
            elapsed_time = end_time - start_time
            st.write(f"⏱️ Running Time: {elapsed_time:.2f} seconds")
            optimal_energy = optimal_3d.get(length, "N/A")
            st.success(f"Minimum HP Energy: {min_energy}")
            st.info(f" Optimal Energy: {optimal_energy}")
            visualize_fold(sequence, best_coords, min_energy)

if __name__ == "__main__":
    main()
