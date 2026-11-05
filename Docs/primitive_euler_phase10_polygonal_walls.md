# FASE 10 — Normale nodale Euler e inizializzazione freestream

## Perimetro

Tre file produttivi modificati: `CPrimitiveEulerSolver.hpp/.cpp` e
`CPrimitiveEulerIteration.cpp`. Nessuna modifica a kernel, fluttuazioni,
MUSCL, Forward Euler, farfield, inlet/outlet o solver nativi. Nessuna nuova
opzione. I dati e i rapporti FASE 9 restano la fotografia della fase precedente.

## Geometria e vincolo

`BuildSlipConstraints` usa le mezze facce già presenti in `CAggregatedGeometry`:

    N_B = sum_f normal_f,       n_B = N_B / ||N_B||.

`normal_f` è già integrata sulla lunghezza della mezza faccia e uscente.
La geometria primitiva verifica già che la sua somma per nodo/marker sia
opposta alla normale nativa `geometry.vertex[marker][v]->GetNormal()`:
quest'ultima segue la convenzione entrante. Non occorre una nuova geometria.
La verifica numerica è ripetuta nel test poligonale, anche invertendo la
connettività di bordo prima del preprocessing SU2.

Per un singolo marker Euler si applica una sola reazione nodale

    rhs_velocity = (I - n_B n_B^T) rhs_velocity_raw.

La reazione resta contabilizzata in `SpatialResult::boundary`, come in FASE 9.
Ogni mezza faccia conserva invece la propria normale: lo stato riflesso e
`ComputeFluctuations(...).minus` sono invariati. **Slip nodale non significa
velocità normale nulla rispetto a ciascuna singola mezza faccia.**

Un marker Euler non rettilineo viene interpretato come approssimazione
poligonale di una superficie liscia: non si inferisce la regolarità da una
soglia angolare. Uno spigolo fisico distinto deve essere identificato con
marker distinti; restano supportati gli incontri di pareti localmente
rettilinee della FASE 9, con i rispettivi vincoli indipendenti. Non si usa
questa possibilità per suddividere artificialmente una curva in marker.

Sono rifiutati: più di due mezze facce dello stesso marker a un nodo,
risultanti nulle/non rappresentabili o cancellate entro il roundoff, e
intersezioni fra marker di slip in cui uno è localmente non rettilineo.
`MARKER_SYM` mantiene il requisito di rettilineità locale. Le tolleranze
256 epsilon servono al riconoscimento del roundoff, non a classificare uno
spigolo come liscio. Un trailing edge acuto non è validato da questa fase.

## Inizializzazione

Il driver, soltanto con `PRIMITIVE_INITIALIZATION=FREESTREAM`, assegna prima
lo stato uniforme e chiama `ProjectFreestreamWallVelocity` al tempo zero:

    velocity_B = velocity_inf - (velocity_inf dot n_B) n_B.

Sono modificati solo i nodi incidenti a `MARKER_EULER`; densità, pressione e
nodi interni restano invariati. Agli spigoli rettilinei esplicitamente
separati si applicano i vincoli già supportati. Il log annuncia il numero
di nodi proiettati prima del primo output/assemblaggio. Non è una correzione
post-step: Forward Euler usa il RHS vincolato già esistente. Nessun clipping
o fallback. UNIFORM, inizializzazioni analitiche e utenti delle API che non
richiedono questa politica devono fornire stati compatibili, come prima.
La proiezione non è disponibile dopo l'avanzamento del tempo.

## Verifica riproducibile

Comandi dalla radice del repository (artefatti soltanto in `/private/tmp`):

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase10-tests
/private/tmp/primitive-phase10-tests/test_driver '[primitive-boundary]'
/private/tmp/primitive-phase10-tests/test_driver '[primitive-residual]~[primitive-spatial-mms],[primitive-time],[primitive-geometry],[primitive-euler]'
/private/tmp/primitive-phase10-tests/test_driver '[primitive-solver]~[primitive-order-audit]~[primitive-spatial-mms]~[primitive-time]~[primitive-boundary]'
python3 UnitTests/check_primitive_euler_polygon_cli.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/primitive-phase10-cli
git diff --check
```

Build riuscita; suite BC: **10 casi, 1422 asserzioni, tutti passati**.
I test coprono risultante analitica (.15,-1), confronto con
le normali native, componente tangenziale non nulla, mancata proiezione
automatica nell'assemblaggio, conservazione iniziale di rho/p e dei nodi
non-wall, cinque passi FE positivi e finiti, e regressioni rettilinee/ruotate.
Una giunzione Euler non-manifold viene rifiutata già dal controllo geometrico
esistente (`Primitive boundary normal inconsistent with SU2`), prima
dell'inizializzazione; non si aggira questo controllo per raggiungere quello
successivo sulle incidenze. È verificato anche il rifiuto della simmetria
non rettilinea.

Smoke integrato: mesh di 6 nodi e 5 triangoli, due segmenti Euler non
collineari e bordo CUSTOM; 5 passi dt=1e-4, 4 restart controllati, slip nodale
entro 2e-13, densità/pressione positive e output finiti. Il checkpoint iniziale
verifica anche la proiezione del driver e i nodi non-wall invariati.

Regressioni selezionate: 42 casi/145234 asserzioni e 36 casi/20877743
asserzioni passati (selezioni sovrapposte, non sommare i conteggi); la seconda
include anche diagnostiche esistenti di consistenza e limiter. Nessuna nuova
campagna di convergenza, nessun confronto con un solver nativo usato come oracolo.

## Limiti

Questa è una verifica funzionale e geometrica, non una validazione
aerodinamica. Restano da studiare accuratezza di bordo e geometria su profili
raffinati, pressioni/forze, gestione del trailing edge e inizializzazione di
problemi fisici più impegnativi. Non si certificano conservazione, ordine
spaziale generale o impermeabilità puntuale di ogni mezza faccia. Il farfield
locale conserva integralmente le limitazioni documentate in FASE 9.
