import sys, time, warnings
sys.path.insert(0, r'C:\Users\Florian Travail\Documents\FlowSom\Perplexity\flowsom_pipeline_pro')
sys.path.insert(0, r'C:\Users\Florian Travail\Documents\FlowSom\Perplexity')
warnings.filterwarnings('ignore')
import logging; logging.disable(logging.CRITICAL)
import numpy as np
from src.core.clustering import FlowSOMClusterer

# Simule un dataset typique : ~80k cellules, 8 marqueurs (ratio=3, ~40k NBM + ~40k patient)
X = np.random.randn(80000, 8).astype('float32')

times = []
for dim in [13, 14, 15]:
    t0 = time.perf_counter()
    cl = FlowSOMClusterer(xdim=dim, ydim=dim, n_metaclusters=8, seed=42, use_gpu=True)
    cl.fit(X)
    dt = time.perf_counter() - t0
    times.append(dt)
    print(f"dim={dim}x{dim} k={dim*dim}: {dt:.2f}s")

avg = sum(times) / len(times)
print(f"\nMoyenne : {avg:.2f}s / run")
print(f"25200 runs : {25200*avg/3600:.1f}h ({25200*avg/60:.0f} min)")
print(f"5040 runs (sans axe MRD) : {5040*avg/3600:.1f}h ({5040*avg/60:.0f} min)")
