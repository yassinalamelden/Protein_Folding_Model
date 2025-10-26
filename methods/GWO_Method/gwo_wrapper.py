import numpy as np

def continuous_to_moves(x, move_set=("F", "L", "R")):
    x = np.array(x).flatten()
    k = len(move_set)
    if x.size == 0:
        return tuple()

    xmin, xmax = float(x.min()), float(x.max())
    if xmax - xmin < 1e-12:
        return tuple(move_set[0] for _ in x)

    x01 = (x - xmin) / (xmax - xmin + 1e-12)
    idx = np.floor(x01 * k).astype(int)
    idx = np.clip(idx, 0, k - 1)
    return tuple(move_set[i] for i in idx)


def _strip_nonletters(s: str) -> str:

    return "".join(ch for ch in s.upper() if "A" <= ch <= "Z")


def _parse_fasta_text(fasta_text: str) -> str:
    lines = fasta_text.strip().splitlines()
    if not lines:
        return ""
    if lines[0].startswith(">"):
        seq_lines = [ln.strip() for ln in lines[1:] if ln and not ln.startswith(">")]
        return _strip_nonletters("".join(seq_lines))

    return _strip_nonletters("".join(lines))


def _aa_to_hp(seq_aa: str) -> str:
    HYDRO = set("AVILMFYWC" + "G")
    return "".join("H" if aa in HYDRO else "P" for aa in seq_aa)


def _auto_to_hp(seq_or_fasta: str) -> str:
    s = seq_or_fasta.strip()
    if not s:
        return ""

    up = s.upper()
    only_hp = set(_strip_nonletters(up)) <= {"H", "P"} and len(_strip_nonletters(up)) > 0
    if only_hp:
        return _strip_nonletters(up)

    aa_seq = _parse_fasta_text(s)
    return _aa_to_hp(aa_seq)

def hp_or_fasta_energy_factory(sequence_or_fasta, fold_from_moves, energy_hp):
    hp_seq = _auto_to_hp(sequence_or_fasta)
    if not hp_seq:
        raise ValueError("Empty sequence provided to hp_or_fasta_energy_factory().")
    dim = len(hp_seq) - 1

    def hp_energy(x):
        x = np.array(x).flatten()
        if x.size < dim:
            x = np.pad(x, (0, dim - x.size))
        else:
            x = x[:dim]

        moves = continuous_to_moves(x)
        try:
            coords = fold_from_moves(moves)
            e = energy_hp(coords, hp_seq)
            return float(e)
        except Exception:
            return 1e6 

    return hp_energy, dim