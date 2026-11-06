# NACA 0012 — Mach 0.2, AOA 0° (FASE 15)

Mesh originale `QuickStart/mesh_NACA0012_inv.su2`, immutata. Le due configurazioni
usano gamma=1.4, R=1, rho_inf=p_inf=T_inf=1, V_inf=.2 sqrt(1.4), REF_AREA=1.
Non sono configurazioni transoniche e non eseguono il caso AOA=2°.

Dal repository:

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_naca15.py \
  --output-dir TestCases/primitive_euler_examples/naca0012_aoa0/results \
  --primitive-steps 3000 --native-steps 600 --seconds-per-run 180
```

Lo script richiede Python, NumPy e Matplotlib; usa solo l'eseguibile SU2_CFD per
l'evoluzione e il writer ParaView nativo. La directory di destinazione deve
essere nuova: nessuna sovrascrittura di una campagna precedente. Limiti
configurabili per iterazioni e durata di ogni run, smoke test e controlli
negativi prima della campagna. I cfg effettivi e i log rimangono in ogni
sottodirectory. Per rigenerare soltanto analisi e grafici:

```sh
python3 UnitTests/run_primitive_euler_naca15.py \
  --output-dir TestCases/primitive_euler_examples/naca0012_aoa0/results \
  --postprocess-only
```

`primitive.cfg`: Forward Euler globale fisico, CFL .05, WLS/MUSCL affine,
limiter NONE, farfield locale caratteristico. Il TE debole viene attivato
**esplicitamente** con `PRIMITIVE_WEAK_TE_GLOBAL_ID=199`; nessun ID è hardcoded
nel solver. Senza opzione il TE originale resta rifiutato. Opzione disponibile
solo per PRIMITIVE_EULER seriale; la validazione topologica della FASE 13C
rimane invariata. Altri nodi di parete mantengono i vincoli forti esistenti.

`native.cfg`: vero EULER conservativo SU2, Roe/MUSCL affine, WLS, limiter NONE,
Euler implicito, MG=2, CFL=10, farfield nativo. È un confronto di stati
stazionari, non di traiettorie temporali equivalenti.

Output principali sotto `results/`:

- `primitive/iterations.pvd`, `native/iterations.pvd`: aprire direttamente in
  ParaView. Le ascisse PVD sono **indici di iterazione**, non tempi fisici; il
  tempo FE primitivo rimane nella history e nei campi. Le sequenze collegano
  direttamente i VTU volumetrici di SU2, non convertiti dallo script. I VTM
  nativi steady possono contenere anche un blocco non indicizzato: vengono
  lasciati intatti, ma esclusi dalla sequenza per evitare snapshot mescolati.
  Per la superficie aprire i file `*/zone_0/airfoil.vtu` delle stesse iterazioni.
- Campi `Density`, `Velocity` (vettore x/y), `Pressure`, `Mach` per entrambi;
  `Speed` diretto per il primitivo. Per il nativo usare Calculator
  `mag(Velocity)` per il modulo. Il writer usa Float32, i restart Float64.
- `primitive/history.csv`: `RHS_CURRENT`, `RHS_INTERIOR`, `RHS_TE` sono norme
  lineari del residuo **ricalcolato sullo stato corrente**, dopo la reazione di
  parete. `MAX_MACH` è nodale, `MAX_TRACE_MACH` comprende le tracce interne e
  fisiche ricostruite. Le righe includono CL/CD e minimi di rho/p.
- `primitive/surface_primitive_pressure_*.csv`: campioni di pressione e Cp alle
  mezze facce. `surface_primitive_te_*.csv`: stato/RHS al TE, le due fluttuazioni
  riflesse D-minus e le due velocità normali delle tracce.
- `native_initial/`: un'iterazione esplicita nativa con CFL=0, senza avanzamento,
  per esportare lo stato iniziale elaborato dalle BC native. Non è uno snapshot
  analitico sintetico e non è una nuova formulazione numerica.
- `native_recheck/`: stesso controllo a CFL=0 sul restart finale, per il
  residuo spaziale finale. Il postprocessing verifica rho/energia e interno
  invariati, e misura l'eventuale riproiezione nativa del momento a parete
  (richiesta inferiore a 1e-9, non un aggiornamento temporale).
- `comparison/`: quattro grafici PNG, `summary.csv`, controlli degli snapshot e
  confronto Cp per lato. L'interpolazione lineare dei Cp nodali nativi è
  valutata alle ascisse primitive nell'intervallo comune; RMS non pesata sui
  punti confrontabili e scarto massimo sono riportati separatamente.

Il controllo VTU legge header, byte counts e array raw, verifica punti/celle,
connettività e campi finiti, positività e Mach nodale<1. Non sostituisce una
verifica grafica interattiva. Per il nativo gli estremi fisici riguardano gli
snapshot salvati; per il primitivo anche ogni riga della history.

La riduzione del residuo è normalizzata separatamente: RHS primitivo in norma
L2 sui volumi duali e RMS nativa della densità conservativa hanno significati
diversi. Una terminazione regolare non certifica convergenza: leggere il
rapporto della FASE 15 per esito e limiti effettivamente misurati.
