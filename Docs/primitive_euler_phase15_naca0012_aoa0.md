# FASE 15 — NACA 0012 subsonico, AOA=0°

## Esito: C — problema identificato, caso non validato

Il solver nativo converge; il primitivo resta in un transitorio regolare ma
molto breve, con resistenza e Cp fortemente differenti dal nativo. Non è stato
ottenuto il confronto richiesto fra **due soluzioni stazionarie**. CL prossimo
a zero e assenza di errori non bastano a validare il caso. Questa conclusione
non dimostra che il metodo primitivo non ammetta uno stato stazionario, né
attribuisce gli scarti a una causa spaziale già isolata.

## A. Configurazioni e modifiche

Branch `feature/primitive-path-solver`, HEAD iniziale `e69ba734b0`, working tree
inizialmente pulito. Mesh originale: 5233 punti, 10216 triangoli, 200 segmenti
airfoil e 400 mezze facce. SHA256 invariato:
`9094b51c2628d3bb4c865d774b2308dbb59a213ad19c007cfcd7d57e7aaeebeb`.
Nessuna modifica al mesh/cfg QuickStart, nessun caso NACA a incidenza diversa.

Entrambi: M_inf=.2, AOA=0°, gamma=1.4, R=1, rho_inf=p_inf=T_inf=1,
V_inf=.23664319132398465, q_inf=.028, scale DIMENSIONAL, REF_AREA=REF_LENGTH=1,
MARKER_MONITORING=(airfoil), Roe/MUSCL affine, WLS, limiter NONE.

| Solver | Avanzamento | CFL | Multigrid | Farfield |
|---|---|---:|---:|---|
| Primitivo | Forward Euler globale fisico | .05 | 0 | LOCAL_CHARACTERISTIC, reservoir locale linearizzato |
| Nativo | Euler implicito, FGMRES/ILU | 10 | 2, V-cycle | BC farfield nativa |

Traiettorie, BC e rappresentazioni sono differenti: il nativo è un riferimento
numerico, non una soluzione esatta. I dati fisici freestream coincidono;
l'inizializzazione elaborata dalle rispettive pareti può differire.

Produzione modificata soltanto in:

- `Common/include/CConfig.hpp`, `Common/src/CConfig.cpp`: opzione
  `PRIMITIVE_WEAK_TE_GLOBAL_ID`, default disabilitato (-1); valori espliciti
  negativi e uso con altri solver rifiutati.
- `SU2_CFD/include/iteration/primitive_euler_weak_te_config.hpp`: conversione
  GlobalIndex→ID locale univoco/domain, seriale, senza ID NACA hardcoded.
- `SU2_CFD/src/solvers/CSolverFactory.cpp`: passa l'ID risolto al costruttore
  esistente, che esegue la validazione topologica della FASE 13C.
- `SU2_CFD/src/output/CPrimitiveEulerOutput.cpp`: con opt-in TE, diagnostiche
  current-state e CSV delle due mezze facce alle scadenze di output.

Altri file: test boundary, nuovo script `UnitTests/run_primitive_euler_naca15.py`,
due cfg e README in `TestCases/primitive_euler_examples/naca0012_aoa0`, questo
rapporto e risultati generati sotto `results/`. Kernel, geometria, MUSCL,
riflessione, fluttuazioni, FE e solver conservativo nativo **non modificati**.
Nessun nuovo pseudo-tempo locale, Kutta o schema di convergenza.

## B. Convergenza misurata

| Solver | Iterazioni | R iniziale → finale | R_final/R_init | Riduzione | CL finale | CD finale | Max Mach |
|---|---:|---|---:|---:|---:|---:|---:|
| PRIMITIVE_EULER | 3000 | .0470270 → .0339300 | .72150 | .142 decadi | -1.22165e-7 | .697670 | .208370 |
| EULER | 600 | 6.39814e-5 → 6.98320e-13 | 1.09144e-8 | 7.962 decadi | -6.74244e-5 | .000409850 | .233443 |

**Norme diverse**, normalizzate separatamente nei grafici: primitivo L2
volumetrica delle quattro componenti del RHS dopo reazione di parete; nativo
RMS della componente conservativa densità, ricavata da rms[Rho] log10. Non si
confrontano le ampiezze assolute delle due norme.

Il RHS primitivo è ricalcolato a ogni riga della history, anche iniziale e
finale; non si usa dt*RHS. Quello interno passa da .0260361 a .0341839,
**aumentando**. Minimi su tutte le iterazioni primitive: rho=.9632637160,
p=.9489443939. Max Mach ricostruito=.2708199245, includendo entrambi i lati
delle facce interne e le tracce interne sulle facce fisiche. Tutti subsonici.

Tempo fisico finale primitivo=.00794954033, dt_finale=2.6622143e-6;
97.18 s wall-clock. Un tempo convettivo di corda è circa 4.23: il run ne copre
solo lo 0.188%. Al costo osservato, un tempo convettivo richiederebbe circa
1.6 milioni di passi/14 ore: **scala indicativa**, non una previsione di
convergenza. Diagnostiche/output e compilazione concorrente incidono sul costo.
Pseudo-tempo locale valutato come seguito possibile, non aggiunto: richiede
un percorso distinto con dt_i, accettazione transazionale e semantica temporale
verificati; non si è alterato FE per forzare un risultato nel budget.

Nativo: 31.43 s wall-clock; minimi sugli snapshot rho=.982158959,
p=.988416910. Il massimo Mach nativo riguarda tutti i punti degli snapshot
salvati, non ogni iterazione né le sue tracce ricostruite.

Residuo nativo finale ricalcolato caricando il restart e facendo una valutazione
esplicita a CFL=0, MG=0, stessi flussi/gradienti sul livello fine. Audit
Float64: coordinate, rho, energia e momento interno **identici**; piccola
riproiezione del momento solo a parete, max 1.08143e-10. Quindi nessun
avanzamento temporale, ma non invarianza bitwise universale. Questo è un
comportamento nativo misurato, non una correzione applicata al suo schema.

## C. Aerodinamica e trailing edge

Negli ultimi 100 passi primitivi, intervallo CL=6.27e-10 ma intervallo
CD=.0039894. Il piccolo CL segnala simmetria, **non stazionarietà**. Negli
ultimi 60 passi nativi: intervallo CL=7.57e-7, CD=4.92e-12.

Cp superiore/inferiore separati per segno di y; interpolazione lineare dei Cp
nodali nativi alle ascisse delle mezze facce primitive, nell'intervallo comune.
400 confronti: RMS non pesata=**2.95207**, massimo assoluto=**8.53976**.
I Cp usano q=.028 in entrambi; i campi Float32 del writer nativo limitano la
precisione del postprocessing. Differenze troppo grandi per rivendicare
accordo, ma il primitivo non è ancora stazionario: non costituiscono una
misura isolata dell'errore asintotico del suo operatore.

TE primitivo globale 199: inizialmente conserva V_inf, senza proiezione forte;
gli altri 199 nodi wall sono proiettati dalla politica FREESTREAM esistente.

| Passo | rho_TE | u_TE | p_TE | norma RHS_TE |
|---:|---:|---:|---:|---:|
| 0 | 1 | .236643 | 1 | 140.245 |
| 600 | .986520 | .206390 | .981186 | 5.45097 |
| 1800 | .987786 | .195040 | .982961 | 2.30326 |
| 3000 | .988518 | .189776 | .984003 | 1.32344 |

v_TE finale=9.67e-10; le due tracce hanno un≈-.02680456 e le due fluttuazioni
D-minus hanno norma≈6.99509e-6. Stato, quattro RHS e quattro componenti di
ciascuna fluttuazione sono nei CSV TE; riassunto in `comparison/te_summary.csv`.
Nessuna oscillazione ampia emerge ai sei checkpoint; non è prova di stabilità
a lungo termine. Residuo TE ancora persistente, non nascosto dalla norma globale.
Il nativo ha u_TE=0, v_TE≈-1.47e-6: differenza verificata nei vincoli nodali,
non prova che uno dei due trattamenti sia fisicamente esatto o un bug universale.

## D. ParaView e riproducibilità

Aprire queste sequenze, mantenendo intatta la struttura `results/`:

- [Primitivo, 6 snapshot](../TestCases/primitive_euler_examples/naca0012_aoa0/results/primitive/iterations.pvd): passi 0,600,1200,1800,2400,3000.
- [Nativo, 6 snapshot](../TestCases/primitive_euler_examples/naca0012_aoa0/results/native/iterations.pvd): iniziale elaborato dalle BC a CFL=0, poi iterazioni accettate 121,241,361,481,600.
- [Grafici e CSV comparativi](../TestCases/primitive_euler_examples/naca0012_aoa0/results/comparison/summary.csv): `residuals.png`, `lift.png`, `drag.png`, `cp.png`, Cp e controlli snapshot.

Sequenze PVD indicizzate per iterazione, non tempo fisico comune. VTU prodotti
dal writer SU2: Density, Velocity vettoriale, Pressure, Mach; Speed diretto nel
primitivo, Calculator `mag(Velocity)` nel nativo. Superfici disponibili nei
VTU airfoil delle stesse directory; Cp primitivo di mezza faccia nei CSV.

Verificati byte counts, array finiti, connettività, punti/celle attesi, positività,
Mach e riferimenti PVD. Verificata anche coerenza delle primitive native con U
entro precisione Float32. I VTM steady nativi contengono blocchi sia indicizzati
sia non indicizzati: i PVD puntano **direttamente ai VTU indicizzati**, evitando
snapshot mescolati senza modificare writer o file originali. Nessun nuovo
writer dei campi, nessuna apertura GUI ParaView eseguita.

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_naca15.py --output-dir TestCases/primitive_euler_examples/naca0012_aoa0/results --primitive-steps 3000 --native-steps 600 --seconds-per-run 180
python3 UnitTests/run_primitive_euler_naca15.py --output-dir TestCases/primitive_euler_examples/naca0012_aoa0/results --postprocess-only
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase15-tests
/private/tmp/primitive-phase15-tests/test_driver '[primitive-native-te-audit],[primitive-weak-te],[primitive-trailing-edge],[primitive-aerodynamics],[primitive-geometry],[primitive-boundary],[primitive-wall-gradients],Forward Euler*,Rejected time steps*'
python3 UnitTests/check_primitive_euler_polygon_cli.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/primitive-phase15-polygon
git diff --check
```

Build riuscita; **31 casi/332678 asserzioni** passati; CLI poligonale passata
(quattro snapshot/tabelle, history CL/CD verificata). Smoke NACA entrambi
passati; rifiuti espliciti senza opt-in, ID assente, nodo non ammissibile,
opzione sul nativo. Percorso MPI esplicitamente rifiutato dal codice;
nessun run multirank effettuato. Nessuna campagna MMS/vortice completa.
Il postprocessing è stato corretto e rieseguito senza nuove simulazioni dopo
aver misurato riproiezione nativa e manifest misti. Nessun commit/push o altra
operazione Git di scrittura. FASE 16 non eseguita.
