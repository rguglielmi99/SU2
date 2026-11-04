# Esempi delle BC fisiche primitive

Ogni cartella `regular`/`jitter` contiene una mesh triangolare condivisa
N=16 e otto configurazioni indipendenti. Gli esempi sono uniformi: servono
per avviare il driver e confrontare i dati fisici. La validazione non
uniforme e acustica è descritta nel rapporto FASE 9.

Dalla radice del repository:

```sh
./ninja -C build SU2_CFD/src/SU2_CFD -j4
cd TestCases/primitive_euler_boundary_examples/regular
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD wall.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD symmetry.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD farfield_subsonic.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD farfield_supersonic.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD inlet_outlet_subsonic.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD inlet_outlet_supersonic.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD walls_farfield.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD walls_inlet_outlet.cfg
```

Gli stessi comandi funzionano dalla cartella `jitter`. Tutti i 16 casi
sono stati eseguiti con l'eseguibile ordinario. Aprire in ParaView le
serie `<nome>_*.vtm`, Apply, blocco Internal. Gli output hanno prefissi
distinti e non sovrascrivono gli altri casi della stessa cartella.

Per il confronto nativo sostituire `SOLVER=PRIMITIVE_EULER` con `EULER`:
i marker e i dati fisici rimangono gli stessi. Il runner
`UnitTests/check_primitive_euler_boundary_cli.py` esegue entrambi e controlla
i checkpoint Float64. Le opzioni PRIMITIVE_* sono specifiche del metodo;
NODE/AGGREGATED, MUSCL affine e FE rimangono obbligatori. Un rank, un thread,
SYSTEM_MEASUREMENTS=SI, REF_DIMENSIONALIZATION=DIMENSIONAL, gas ideale.

`PRIMITIVE_INITIALIZATION=FREESTREAM` usa pressione, densità oppure
temperatura, Mach e AOA nativi. Nei casi subsonici rho=p=R=1 e u=.4;
nei supersonici u=2. I dati totali dell'inlet sono coerenti con questi
stati. MASS_FLOW, verificato dal runner e dal benchmark, significa rho
e modulo della velocità come nel parser SU2, non portata integrata.

MARKER_FAR richiede `PRIMITIVE_BOUNDARY_MODE=LOCAL_CHARACTERISTIC`: è un
reservoir caratteristico locale con dati freestream nativi, **non ancora
equivalente alla chiusura nativa a ampiezza finita**. Per il vecchio stato
esterno prescritto o il vortice analitico
continuare a usare MARKER_CUSTOM e gli esempi della FASE 8. Un vecchio
FAR con modalità UNIFORM/ISENTROPIC_VORTEX o BOUNDARY_STATE esplicita viene
rifiutato: non cambia silenziosamente significato. Per ora il cfg con FAR
attivo non accetta questi dati globali neppure per un altro marker CUSTOM.
I casi misti analitici sono verificati via API.

Pareti e simmetrie richiedono velocità iniziale nodale tangente. I tratti
devono essere localmente rettilinei; spigoli reali vanno separati per
marker. Suddividere una curva in tanti marker non realizza una parete curva.
Non sono supportati backflow dell'outlet di sola pressione, periodicità,
pareti curve, RIEMANN/GILES/engine/actuator disk o interfacce multizona.
Un dato incompatibile viene rifiutato, senza clipping o cambio di solver.

Formule, risultati e limiti sono nel
[rapporto FASE 9](../../Docs/primitive_euler_phase9_boundary_conditions.md).
