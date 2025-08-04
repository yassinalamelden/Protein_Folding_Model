# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 12:07:23 2025

@author: Yassin
"""

import random
import heapq
import matplotlib.pyplot as plt
import math
import sys


# Directions for 2D lattice: Up, Right, Down, Left
square_directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]

# Directions for 2D Trianguler Lattice:
trianguler_directions = [
    (1.0, 0.0),                                     # 0: Right
    (0.5, math.sqrt(3)/2),                          # 1: Up-Right
    (-0.5, math.sqrt(3)/2),                         # 2: Up-Left
    (-1.0, 0.0),                                    # 3: Left
    (-0.5, -math.sqrt(3)/2),                        # 4: Down-Left
    (0.5, -math.sqrt(3)/2)                          # 5: Down-Right
]


# Known optimal energies for reference
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


# Converting Fasta to HP Sequence
def fasta_to_hp(fasta_sequence):
    """
    Converts a FASTA amino acid sequence to HP model using:
    H = Hydrophobic {A, G, I, L, M, F, P, W, V}
    P = Polar       {R, N, D, C, E, Q, H, K, S, T, Y}
    """
    
    h_residues = {'A', 'G','I', 'L', 'M', 'F', 'P', 'W', 'V'}
    p_residues = {'R', 'N', 'D', 'C', 'E', 'Q', 'H','K', 'S', 'T', 'Y'}
    
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


# Protein Fold Square
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


# Protein Fold Trianguler
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


def compactness_penalty(coords):
    if not coords or len(coords) < 2:
        return 0
    xs, ys = zip(*coords)
    area = (max(xs) - min(xs) + 1) * (max(ys) - min(ys) + 1)
    return 1 / area if area > 0 else 0


def multi_point_crossover(p1, p2):
    if len(p1) < 3: 
        return p1
    cut1, cut2 = sorted(random.sample(range(1, len(p1)), 2))
    return p1[:cut1] + p2[cut1:cut2] + p1[cut2:]


def mutate(moves, rate):
    return [random.randint(0, 2) if random.random() < rate else m for m in moves]


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


# Visualization of the Fold

def visualize_fold(sequence, coords, min_energy):
    if not coords or len(coords) != len(sequence):
        print("Cannot visualize: invalid folding")
        return

    xs, ys = zip(*coords)
    colors = ['red' if c == 'H' else 'blue' for c in sequence]

    plt.figure(figsize=(8, 8))
    for i in range(1, len(coords)):
        x1, y1 = coords[i-1]
        x2, y2 = coords[i]
        plt.plot([x1, x2], [y1, y2], 'k-')  # black line between amino acids

    for (x, y), color in zip(coords, colors):
        plt.plot(x, y, 'o', markersize=12, color=color)

    plt.title(f"Protein Fold Visualization (Red=H, Blue=P)\nMinimum Energy = {min_energy}", fontsize = 18)
    plt.axis('equal')
    plt.xticks(fontsize = 12, fontweight = 'bold')
    plt.yticks(fontsize = 12, fontweight = 'bold')
    plt.grid(True)
    plt.show()

# Run Simulation
if __name__ == "__main__":
    while True:
        sequence_type = input("Input Type:\n1. Fasta Sequence (amino acids)\n2. HP sequence\nEnter your choice: ").strip()
        if sequence_type == '1':
            try:
               fasta_sequence = "".join(input("Enter the fasta sequence: ").upper().split())
               sequence = fasta_to_hp(fasta_sequence)
               print(f"Converted HP Sequence: {sequence}")
            except ValueError as e:
                print("Error", e)
                continue
        elif sequence_type == '2':
            sequence = "".join(input("Enter the sequence (only H and P): ").upper().split())
            if not sequence or any(c not in "HP" for c in sequence):
                print("Invalid sequence. Please enter a string containing only 'H' and 'P'.")
                continue
        else:
            print("Invalid choice. Please enter 1 or 2.")
            continue
        
        length = len(sequence)
        print("Sequence Length: ", length)
        lattice_type = input("Choose from menu: \n1. Square 2D Lattice \n2. Tranguler 2D Lattice \n3. Exit the program \nEnter your choice: ").strip().lower()
        if lattice_type == '1':
            print("Square 2D Lattice: Loading...\n")
            min_energy, best_coords = square_genetic_algorithm(sequence, generations=5000, pop_size=200)
        elif lattice_type == '2':
            print("Trianguler 2D Lattice: Loading...\n")
            min_energy, best_coords = trianguler_genetic_algorithm(sequence, generations=5000, pop_size=200)
        else:
            print("Thank you. Goodbye!")
            sys.exit()
    
        # Final Output
        print("Minimum HP Energy:", min_energy)
        if length in optimal_2d:
            optimal = optimal_2d[length]
            deviation = min_energy - optimal
            print("Ideal Optimal Energy:", optimal)
            print("Deviation from Optimal:", deviation)
        else:
            optimal_2d = min_energy
            print("optimal energy for the sequence: ", optimal_2d)
        
        # Visualization
        visualize_fold(sequence, best_coords, min_energy)
    
        again = input("\nDo you want to run another sequence again? (yes/no): ").lower()
        if again != 'yes':
            print("Exiting Program. Goodbye!")
            sys.exit()
    
    
