# FASE 13B — Audit del trailing edge nativo SU2

## Esito e perimetro

**La mesh NACA originale rimane rifiutata da PRIMITIVE_EULER.** Non è stato
implementato un prototipo debole né modificato alcun file produttivo in
questa fase. Le modifiche produttive FASE 13 già presenti restano inalterate.
Modificato solo `primitive_euler_boundary_tests.cpp` con un test
`[primitive-native-te-audit]`; nuovo artefatto: questo rapporto.

L'audit smentisce l'ipotesi che la parete nativa sia soltanto un flusso di
pressione debole: utilizza una proiezione forte della quantità di moto.
Il fatto che il caso NACA sia eseguibile non certifica la correttezza
aerodinamica del singolo valore al trailing edge.

## Percorso nativo verificato nel codice

Riferimenti ai file sotto SU2_CFD, salvo indicazione diversa:

| Passaggio | File, funzione, righe verificate |
|---|---|
| Dispatch forte, dopo gli altri contributi di bordo | `src/integration/CIntegration.cpp`, `Space_Integration`, 201–205 |
| Euler-wall delega alla simmetria | `include/solvers/CFVMFlowSolverBase.hpp`, `BC_Euler_Wall`, 1130–1134 |
| Normale grezza o corretta per intersezioni | `include/solvers/CFVMFlowSolverBase.inl`, `BC_Sym_Plane`, 1193–1207 |
| Proiezione del buffer conservativo precedente | stesso metodo, 1259–1283 |
| Proiezione del residuo di quantità di moto | stesso metodo, 1285–1292 |
| Premoltiplicazione dei blocchi Jacobiani diagonale/vicini per il proiettore | stesso metodo, 1294–1319 |
| Proiezione della truncation error multigrid | stesso metodo, 1322–1329 |
| FE usa l'aggiornamento scalare | `include/solvers/CFVMFlowSolverBase.hpp`, `Explicit_Iteration_impl`, 803–887, in particolare 843–844 |
| Semantica effettiva di quell'aggiornamento | `include/variables/CVariable.hpp`, `AddSolution(iPoint,iVar,...)`, 339–340 |
| Aggiornamento implicito con la stessa overload scalare | `include/solvers/CFVMFlowSolverBase.hpp`, `CompleteImplicitIteration_impl`, 984–1000 |

Per parete statica, BC_Sym_Plane modifica **Solution_Old**, non direttamente
Solution: m_old <- (I-nnᵀ)m_old. Proietta allo stesso modo il residuo di
momento. L'overload scalare AddSolution esegue **Solution=Solution_Old+incremento**,
non Solution+=incremento: la proiezione diventa quindi una modifica dello
stato al successivo aggiornamento. Non c'è una proiezione separata generica
dell'incremento nella BC ordinaria; il Jacobiano implicito è coerentemente
modificato. Densità ed energia non vengono proiettate da questo blocco.

Nel ramo statico non viene aggiunto esplicitamente (0,p*Sx,p*Sy,0).
Non è quindi un trattamento esclusivamente a flusso di pressione. A livello
algebrico, un eventuale p_i*S_i sarebbe parallelo alla normale aggregata e
verrebbe eliminato dalla proiezione del momento: non significa che il codice
calcoli separatamente le pressioni sulle due mezze facce. Il ramo **mobile**
(inl 1211–1257) calcola contributi upwind energetici, irrilevanti qui.

### Normali e inizializzazione

`Common/src/geometry/CPhysicalGeometry.cpp`, `SetBoundControlVolume`,
6836–6922, accumula i contributi alle normali dei vertici. In 2D
`Common/src/geometry/dual_grid/CVertex.cpp`, `SetNodes_Coord`, 48–50,
somma il vettore ruotato del segmento nodo/punto medio. La convenzione
nativa è entrante, opposta alle normali primitive uscenti.

`Common/src/geometry/CGeometry.cpp`, `ComputeModifiedSymmetryNormals`,
2857–2970 e seguito, unisce normali di marker curvi condivisi e gestisce
intersezioni. Su **un solo marker airfoil** non crea un'eccezione al TE.
Il test verifica che `symmetryNormals[wall]` non contenga una sostituzione
per il nodo 199. Nessuna logica Kutta/TE compare nel percorso verificato.

`CEulerSolver.cpp:142,296` costruisce le variabili dal freestream;
`src/variables/CEulerVariable.cpp:53–75` assegna U_inf a Solution e lo copia
in Solution_Old, senza proiezione di parete. La preparazione ordinaria copia
nuovamente lo stato precedente prima dell'integrazione spaziale
(`CMultiGridIntegration.cpp:583–599`); la BC interviene successivamente.

### Multigrid distinto dal percorso ordinario

`CFVMFlowSolverBase.inl`, `MultigridProjectEulerWall`, 3076–3133, proietta
esplicitamente il buffer indicato da use_solution_old. I call site sono
`CMultiGridIntegration.cpp:788` per la correzione coarse in Solution_Old e
`:1074` per la soluzione ristretta. Non sono la spiegazione esclusiva del
vincolo: la BC ordinaria lo applica anche con MGLEVEL=0.
I cfg NACA locali usano MGLEVEL=3 ed EULER_IMPLICIT
(`QuickStart/inv_NACA0012.cfg:127,157`; cfg TestCases:63,75).
La diagnostica seguente isola invece il percorso FE senza multigrid.

## Nodo 199: misura riproducibile

Il test legge direttamente `QuickStart/mesh_NACA0012_inv.su2`, senza
riscriverla, e trova il nodo tramite GlobalIndex=199. Coordinate (1,0),
normale nativa integrata **(+3.632896519e-5,0)**; normale uscente opposta.
Le normali locali sono circa (-.143816,±.989604): il proiettore aggregato
elimina x, non applica separatamente due normali.

Stato diagnostico dimensionale rho=p=1, u=.4, v=0, gamma=1.4:

| Fase nativa misurata | u corrente | u precedente | E conservativa | p ricavata |
|---|---:|---:|---:|---:|
| Costruzione freestream | .4 | .4 | 2.58 | 1 |
| Dopo BC_Euler_Wall | .4 | 0 | 2.58 | 1 nello stato corrente |
| Dopo FE con residuo azzerato | 0 | 0 | 2.58 | 1.032 |

Si immette anche un residuo **sintetico**, (.1,.2,.3,.4), per verificare
il proiettore: dopo BC diventa (.1,0,.3,.4). Non è il residuo fisico NACA.
L'aggiornamento a residuo zero isola l'effetto del buffer precedente; non
è una simulazione aerodinamica né una misura della pressione a convergenza.
Con v iniziale non nullo, la sola proiezione al TE conserverebbe v:
non impone necessariamente velocità totale nulla per ogni incidenza.

Nel primitivo, come dedotto dal percorso produttivo verificato, il driver assegna prima W_inf
(`CPrimitiveEulerIteration.cpp:93–98`), poi ProjectFreestreamWallVelocity
costruisce i vincoli e rifiuta il nodo 199 prima di applicare la proiezione.
Non esiste quindi uno stato evoluto primitivo del TE in questa fase.
La vecchia media FASE 10 avrebbe eliminato u conservando rho/p, mentre il
nativo conserva rho/E: le due operazioni non sono termodinamicamente identiche.

## Ipotesi debole: proposta, non implementazione

Eliminare soltanto al TE il vincolo nodale lascerebbe attive riflessione e
D-minus su ciascuna mezza faccia (`CPrimitiveEulerSolver.cpp:330–340`), con
il nodo determinato dall'intero RHS. Un vertice singolare non richiede
necessariamente una tangente nodale unica in una formulazione debole, ma
ciò **non dimostra** consistenza/stabilità di questa specifica chiusura
path-conservative. Slip nodale e impermeabilità debole delle facce sarebbero
proprietà diverse da misurare. Non sarebbe una copia del trattamento nativo,
che mantiene il vincolo forte; nessun flusso conservativo può essere inserito
direttamente nel residuo primitivo.

La mesh identifica soltanto airfoil/farfield, senza metadati di vertice TE.
Angolo e settore reentrante non distinguono in generale TE, spigolo di altro
tipo o difetto di mesh; massimo x e ID 199 sarebbero euristiche/whitelist.
Prima di un prototipo servono un'identificazione TE esplicita e validata
(vertice, due facce e settore fluido), poi prove dedicate della chiusura
debole, dei bilanci e del raffinamento. Non è stata introdotta tale politica.

## Test ed esito

Comandi dalla radice del repository:

```sh
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase13b-build
/private/tmp/primitive-phase13b-build/test_driver '[primitive-native-te-audit]'
/private/tmp/primitive-phase13b-build/test_driver '[primitive-trailing-edge],[primitive-boundary]'
git diff --check
```

Diagnostica: **1 caso/23 asserzioni**; regressioni **14 casi/1996 asserzioni**,
tutti passati. SHA256 della mesh invariato rispetto alla FASE 13.
Resta attivo il rifiuto di TE acuti e angoli non supportati; nessun prototipo,
simulazione stazionaria completa, commit o push. La NACA originale è
**ancora rifiutata**, non supportata sperimentalmente. Kutta e circolazione
restano fuori dal perimetro.
