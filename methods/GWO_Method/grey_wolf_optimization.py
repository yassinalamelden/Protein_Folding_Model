import numpy as np
from dataclasses import dataclass
import time

@dataclass
class GWOConfiguration:
    pop_size: int = 30
    max_iter: int = 300
    lb: float = -1.0
    ub: float = 1.0
    seed: int | None = None

class GWO:
    def __init__(self, dim, fitness_fn, cfg: GWOConfiguration = GWOConfiguration()):
        self.dim = dim
        self.fitness_fn = fitness_fn
        self.cfg = cfg
        self.rng = np.random.default_rng(cfg.seed)

        self.lb = float(cfg.lb)
        self.ub = float(cfg.ub)

        self.X = None
        self.fitness = None

        self.alpha = None
        self.beta  = None
        self.delta = None
        self.f_alpha = None
        self.f_beta  = None
        self.f_delta = None

        self.history = {"best_f": []}

    def _init_population(self):
        self.X = self.rng.uniform(self.lb, self.ub, size=(self.cfg.pop_size, self.dim))

    def _clip(self, X):
        return np.clip(X, self.lb, self.ub)
    
    def _evaluate(self):
        self.fitness = np.apply_along_axis(self.fitness_fn, 1, self.X)
        order = np.argsort(self.fitness)
        self.alpha = self.X[order[0]].copy()
        self.beta  = self.X[order[1]].copy()
        self.delta = self.X[order[2]].copy()
        self.f_alpha = float(self.fitness[order[0]])
        self.f_beta  = float(self.fitness[order[1]])
        self.f_delta = float(self.fitness[order[2]])

    def minimize(self):
        start = time.time()
        self._init_population()
        self.X = self._clip(self.X)
        self._evaluate()
        self.history["best_f"] = [self.f_alpha]

        for t in range(self.cfg.max_iter):
            a = 2.0 - 2.0 * (t / (self.cfg.max_iter - 1 if self.cfg.max_iter > 1 else 1))

            for i in range(self.cfg.pop_size):
                r1 = self.rng.random(self.dim); r2 = self.rng.random(self.dim)
                A1 = 2*a*r1 - a
                C1 = 2*r2
                D_alpha = np.abs(C1*self.alpha - self.X[i])
                X1 = self.alpha - A1*D_alpha

                r1 = self.rng.random(self.dim); r2 = self.rng.random(self.dim)
                A2 = 2*a*r1 - a
                C2 = 2*r2
                D_beta = np.abs(C2*self.beta - self.X[i])
                X2 = self.beta - A2*D_beta

                r1 = self.rng.random(self.dim); r2 = self.rng.random(self.dim)
                A3 = 2*a*r1 - a
                C3 = 2*r2
                D_delta = np.abs(C3*self.delta - self.X[i])
                X3 = self.delta - A3*D_delta

                self.X[i] = (X1 + X2 + X3) / 3.0

            self.X = self._clip(self.X)

            self._evaluate()
            self.history["best_f"].append(self.f_alpha)

            if (t % 10 == 0) or (t == self.cfg.max_iter - 1):
                print(f"Iter {t+1}/{self.cfg.max_iter}  best_f = {self.f_alpha}")

        elapsed = time.time() - start
        return self.f_alpha, self.alpha.copy(), self.history
    

if __name__ == "__main__":
    def sphere(x):
        return float(np.sum(x**2))
    
    cfg = GWOConfiguration(pop_size=40, max_iter=200, lb=-5.0, ub=5.0, seed=42)
    opt = GWO(dim=10, fitness_fn=sphere, cfg=cfg)
    best_f, best_x, hist = opt.minimize()
    print("Best fitness:", best_f)
    print("Best x (first 5 dims):", best_x[:5])
    print("Iterations logged:", len(hist["best_f"]))
