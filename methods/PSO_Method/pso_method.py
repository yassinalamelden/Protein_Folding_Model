# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:32:20 2025

@author: Yassin
"""

import random
import math
import time
import sys
import matplotlib.pyplot as plt

square_directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]


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


### Converting FASTA to HP Sequence
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


# Compactness penalty (encourages dense folds)
def compactness_penalty(coords):
    if not coords or len(coords) < 2:
        return 0
    xs, ys = zip(*coords)
    area = (max(xs) - min(xs) + 1) * (max(ys) - min(ys) + 1)
    penalty = 1 / area if area > 0 else 0
    return penalty * 0.2 if len(coords) <= 20 else penalty


### Particle representation
class Particle:
    def __init__(self, sequence):
        self.sequence = sequence
        self.position = generate_random_fold(len(sequence))  # list of moves
        self.coords = square_fold_protein(sequence, self.position)
        self.energy = true_energy(sequence, self.coords) if self.coords else float('inf')
        self.best_position = self.position[:]
        self.best_energy = self.energy
        self.velocity = []


### Exploratory operators


def discrete_velocity_update(position, velocity, pbest, gbest, w, c1, c2):          # pbest : Personal best position , gbest : Global best position , w : inertia weight , c1 : cognitive coefficient , c2 : social coefficient
    new_velocity = []
    for i in range(len(position)):
        if random.random() < w and velocity:
            new_velocity.append(random.choice(velocity))
        if random.random() < c1:
            new_velocity.append((i, pbest[i]))
        if random.random() < c2:
            new_velocity.append((i, gbest[i]))
    return new_velocity


def apply_velocity(position, velocity):
    new_position = position[:]
    for (idx, move) in velocity:
        new_position[idx] = move
    return new_position


def compact_seed_moves(n):
    if n <= 1:
        return []
    moves = []
    straight_budget = max(1, n // 6)
    for i in range(n - 1):
        if i < straight_budget:
            moves.append(0)
        else:
            moves.append(1 if (i % 2 == 0) else 2)

    for _ in range(200):
        if square_fold_protein("H"*n, moves) is not None:
            return moves

        j = random.randrange(0, n - 1)
        moves[j] = random.randint(0, 2)

    return generate_random_fold(n)


### Exploratory operators

# Copies a block of moves from another fold
def segment_pull(base, guide, seg_len):
    if not base or not guide or len(base) != len(guide):
        return base[:]
    L = len(base)
    if L == 0:
        return base[:]
    seg = max(1, min(seg_len, L))
    start = random.randint(0, L - seg)
    new_pos = base[:]
    new_pos[start:start+seg] = guide[start:start+seg]
    return new_pos


# Copies blocks from both personal best and global best
def multi_segment_pull(base, guide1, guide2, seg_lens=(4, 8)):
    pos = segment_pull(base, guide1, seg_lens[0])
    pos = segment_pull(pos,  guide2, seg_lens[1])
    return pos


# Pivot rotation: flips turns beyond a pivot
def pivot_move(moves, pivots=2):
    if not moves:
        return moves[:]
    pos = moves[:]
    for _ in range(pivots):
        i = random.randrange(0, len(pos))
        for j in range(i, len(pos)):
            if pos[j] == 1:
                pos[j] = 2
            elif pos[j] == 2:
                pos[j] = 1
    return pos


def maybe_reset_velocity(p, iteration, period=50, fraction=0.10):
    p.velocity = []


### Particle re-initialization (random restart)
def hybrid_restart_particle(sequence, gbest_position=None, keep_fraction=0.0):
    n = len(sequence)
    if gbest_position and keep_fraction > 0.0:
        base = gbest_position[:]
        L = len(base)
        k = max(1, int(L * (1.0 - keep_fraction)))
        idxs = random.sample(range(L), k)
        for i in idxs:
            base[i] = random.randint(0, 2)
        pos = base
    else:
        pos = compact_seed_moves(n) if random.random() < 0.5 else generate_random_fold(n)

    for _ in range(200):
        coords = square_fold_protein(sequence, pos)
        if coords is not None:
            E = true_energy(sequence, coords)
            return pos, coords, E
        j = random.randrange(0, len(pos))
        pos[j] = random.randint(0, 2)

    pos = generate_random_fold(n)
    coords = square_fold_protein(sequence, pos)
    E = true_energy(sequence, coords) if coords else float('inf')
    return pos, coords, E


# PSO Main Block
def square_pso(sequence, swarm_size=150, max_iter=5000,
               w_start=0.95, w_end=0.30,
               c1_start=1.6, c1_end=2.2,
               c2_start=1.0, c2_end=2.6,
               restart_interval=200, restart_fraction=0.25,
               velocity_reset_period=50, velocity_reset_fraction=0.10):
    
    swarm = []
    n = len(sequence)
    half = swarm_size // 2

    for _ in range(half):
        p = Particle(sequence)
        p.position = compact_seed_moves(n)
        p.coords   = square_fold_protein(sequence, p.position)
        p.energy   = true_energy(sequence, p.coords) if p.coords else float('inf')
        p.best_position = p.position[:]
        p.best_energy   = p.energy
        swarm.append(p)

    for _ in range(swarm_size - half):
        swarm.append(Particle(sequence))

    gbest_particle = min(swarm, key=lambda q: q.energy)
    gbest_position = gbest_particle.position[:]
    gbest_energy   = gbest_particle.energy

    for iteration in range(max_iter):
        # parameters update
        t = iteration / max_iter
        w  = w_start  - (w_start  - w_end)  * t
        c1 = c1_start + (c1_end   - c1_start) * t
        c2 = c2_start + (c2_end   - c2_start) * t

        worst_count = max(1, int(swarm_size * velocity_reset_fraction)) if (iteration % velocity_reset_period == 0 and iteration > 0) else 0
        worst_idxs = []
        if worst_count > 0:
            worst_idxs = sorted(range(swarm_size), key=lambda i: swarm[i].energy, reverse=True)[:worst_count]

        for idx, p in enumerate(swarm):
 
            if idx in worst_idxs:
                maybe_reset_velocity(p, iteration, period=velocity_reset_period, fraction=velocity_reset_fraction)

            candidate = p.position[:]

            if p.velocity and random.random() < w:
                for (j, mv) in p.velocity:
                    candidate[j] = mv

            seg_small = max(2, int(0.10 * len(candidate)))
            seg_large = max(3, int(0.20 * len(candidate)))
            if random.random() < min(1.0, c1 / 2.5):
                candidate = segment_pull(candidate, p.best_position, seg_large)
            if random.random() < min(1.0, c2 / 2.5):
                candidate = segment_pull(candidate, gbest_position, seg_small)

            # Multi-segment pull
            if random.random() < 0.25:
                candidate = multi_segment_pull(candidate, p.best_position, gbest_position, seg_lens=(seg_small, seg_large))
            
            # Pivot mutation
            if random.random() < 0.30:
                candidate = pivot_move(candidate, pivots=1)
            
            # Validate fold
            coords = square_fold_protein(sequence, candidate)
            if coords is None:
                tmp = segment_pull(p.position, p.best_position, seg_small)
                coords = square_fold_protein(sequence, tmp)
                if coords is None:
                    tmp = segment_pull(p.position, gbest_position, seg_small)
                    coords = square_fold_protein(sequence, tmp)
                    if coords is None:
                        tmp = p.position[:]
                        coords = p.coords
                candidate = tmp

            energy = true_energy(sequence, coords) - compactness_penalty(coords)
            
            # Update particle
            p.position = candidate
            p.coords   = coords
            p.energy   = energy

            # Update velocity
            new_vel = []
            for j, (mj_old, mj_new) in enumerate(zip(p.best_position, candidate)):
                if mj_old != mj_new:
                    new_vel.append((j, mj_new))

            if len(new_vel) > max(5, len(candidate)//4):
                new_vel = random.sample(new_vel, max(5, len(candidate)//4))
            p.velocity = new_vel

            # Update personal and global bests
            if energy < p.best_energy:
                p.best_energy   = energy
                p.best_position = candidate[:]

            if energy < gbest_energy:
                gbest_energy   = energy
                gbest_position = candidate[:]

        # Periodic restart of worst particles
        if restart_interval and (iteration + 1) % restart_interval == 0:
            swarm.sort(key=lambda q: q.energy, reverse=True)
            rcount = max(1, int(swarm_size * restart_fraction))
            for i in range(rcount):
                keep_frac = 0.25 if random.random() < 0.5 else 0.0
                pos, coords, E = hybrid_restart_particle(sequence, gbest_position, keep_fraction=keep_frac)
                swarm[i].position = pos
                swarm[i].coords   = coords
                swarm[i].energy   = E
                swarm[i].best_position = pos[:]
                swarm[i].best_energy   = E
                swarm[i].velocity = []
                if E < gbest_energy:
                    gbest_energy   = E
                    gbest_position = pos[:]

    best_coords = square_fold_protein(sequence, gbest_position)
    return int(round(gbest_energy)), best_coords


### Visualization
def visualize_fold(sequence, coords, min_energy):
    if not coords or len(coords) != len(sequence):
        print("Cannot visualize: invalid folding")
        return
    dim = len(coords[0])
    colors = ['red' if c == 'H' else 'blue' for c in sequence]

    if dim == 2:
        xs, ys = zip(*coords)
        plt.figure(figsize=(8, 8))
        for i in range(1, len(coords)):
            x1, y1 = coords[i - 1]
            x2, y2 = coords[i]
            plt.plot([x1, x2], [y1, y2], 'k-')
        for (x, y), color in zip(coords, colors):
            plt.plot(x, y, 'o', markersize=12, color=color)
        plt.title(f"2D Protein Fold Visualization (Red=H, Blue=P)\nMinimum Energy = {min_energy}", fontsize=18)
        plt.xticks(fontsize = 12, fontweight = 'bold')
        plt.yticks(fontsize = 12, fontweight = 'bold')
        plt.axis('equal')
        plt.grid(True)
        plt.show()


####---- Main Program ----####
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
        print("Sequence Length:", length)
        
        structure_type = input("\nChoose from menu: \n1. 2D Structure \n2. Exit the program \nEnter your choice: ").strip().lower()
        if structure_type == '1':
            lattice_type = input("\nChoose from menu: \n1. Square 2D Lattice \n2. Exit the program \nEnter your choice: ").strip().lower()
            if lattice_type == '1':
                print("Running Square 2D Lattice...\n")
                
                start_time = time.time()
                min_energy, best_coords = square_pso(
                    sequence,
                    swarm_size=150,
                    max_iter=5000,
                    w_start=0.95, w_end=0.3,
                    c1_start=1.6, c1_end=2.2,
                    c2_start=1.0, c2_end=2.6,
                    restart_interval=200,
                    restart_fraction=0.25,
                    velocity_reset_period=50,
                    velocity_reset_fraction=0.10
                )

                end_time = time.time()
                total_time = (end_time - start_time) / 60 
                
                print("Minimum HP Energy:", min_energy)
                if length in optimal_2d:
                    optimal = optimal_2d[length]
                    deviation = min_energy - optimal
                    print("Ideal Optimal Energy:", optimal)
                    print("Deviation from Optimal:", deviation)
                else:
                    optimal_2d = min_energy
                    print("optimal energy for the sequence: ", optimal_2d)
            else:
                print("Thank you. Goodbye!")
                sys.exit()

        else:
                print("Thank you. Goodbye!")
                sys.exit()

        print(f"\nTotal Runtime: {total_time:.2f} minutes")
        
        visualize_fold(sequence, best_coords, min_energy)
    
        again = input("\nDo you want to run another sequence again? (yes/no): ").lower()
        if again != 'yes':
            print("Exiting Program. Goodbye!")
            sys.exit()
            
            
