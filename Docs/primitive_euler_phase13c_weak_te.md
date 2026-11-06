# FASE 13C — Prototipo debole esplicito del trailing edge

## Opt-in e perimetro

Modificati `CPrimitiveEulerSolver.hpp/.cpp` e
`primitive_euler_boundary_tests.cpp`; nuovi file: questo rapporto e
`primitive_euler_phase13c_weak_te.csv` (sei campioni).
Il working tree conserva anche il lavoro precedente delle FASI 13/13B.

Il costruttore accetta un quarto argomento opzionale, una lista di **ID locali
della geometria originale**:

    CPrimitiveEulerSolver solver(geometry, config, 0, {localTip});

Default: lista vuota. Factory e driver ordinario non passano autorizzazioni;
nessun parametro cfg, nome di marker, coordinata o ID NACA è codificato nel
solver. Il rifiuto produttivo ordinario della NACA resta attivo.

L'autorizzazione è immutabile per istanza, validata nel costruttore prima
dei vincoli/inizializzazione. Questa versione ristretta richiede un nodo
di dominio con esattamente due mezze facce, di **un solo marker Euler**, due
vicini distinti, normali non collineari e non coorientate, risultante non
degenere e settore fluido reentrante. Sono rifiutati duplicati, ID fuori
dominio, nodi interni/farfield/lisci e intersezioni. Non si tratta di un
classificatore automatico: è il chiamante che dichiara esplicitamente il TE.

Al solo nodo validato BuildSlipConstraints restituisce nessun vincolo forte;
ProjectFreestreamWallVelocity lo esclude e non lo conta come proiettato.
Il RHS resta quello completo. Il calcolo del settore fluido FASE 13 è
raccolto in un helper condiviso, senza duplicare operatori. Su tutti gli
altri nodi restano i controlli e la politica precedenti, compresi quelli
non-manifold e di degenerazione.

**Invariati:** riflessione, D-minus, normali, kernel, MUSCL/gradienti/limiter,
positivity scaling, FE, farfield e solver nativo. Nessuna condizione di Kutta,
esenzione implicita o correzione conservativa.

## Geometria e stato iniziale verificati

Il test legge direttamente `QuickStart/mesh_NACA0012_inv.su2`, senza
riscriverla. Prima dell'opt-in trova GlobalIndex=199, verifica (1,0), marker
airfoil, rami verso GlobalIndex 0/198, due punti di quadratura distinti,
normali uscenti (-.143815653,±.989604496) e settore fluido 343.462595°.
Un'incoerenza causa fallimento, non selezione di un altro nodo.

Configurazione: rho=p=1, u=.4, v=0, gamma=1.4, gas constant=1, dimensionale;
NODE/AGGREGATED/NODE_POINT, WLS, MUSCL affine, limiter NONE, nessuna sorgente.
La callback farfield usa il vero adattatore LOCAL_CHARACTERISTIC con dati
coerenti col freestream, senza fallback prescritto.

Senza opt-in, inizializzazione e residuo sono rifiutati. Con opt-in, il TE
resta esattamente (.4,0) dopo l'inizializzazione; gli altri **199** nodi
airfoil ricevono la consueta proiezione, verificata contro la normale nativa.
La vecchia proiezione aggregata al TE avrebbe dato (0,0).

## Evoluzione breve

Cinque passi FE costanti, dt=9.7266253584e-7 ricavato dall'EstimateTimeStep
iniziale con CFL=.02. T finale=4.8633126792e-6: intervallo fisico molto breve.

| Passo | rho TE | u TE | v TE | p TE | norma RHS TE | max norma RHS vicini |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | .400000000 | 0 | 1 | 237.082 | 207.706 |
| 1 | .999867334 | .399967141 | -1.039e-12 | .999814267 | 236.188 | 207.374 |
| 2 | .999735185 | .399934144 | -1.901e-12 | .999629270 | 235.297 | 207.041 |
| 3 | .999603552 | .399901009 | -2.587e-12 | .999445005 | 234.409 | 206.707 |
| 4 | .999472433 | .399867737 | -3.100e-12 | .999261470 | 233.524 | 206.371 |
| 5 | .999341826 | .399834329 | -3.441e-12 | .999078663 | 232.642 | 206.035 |

Le norme sommano i quadrati delle quattro componenti primitive, nelle unità
scelte. I vicini sono i nodi adiacenti dell'arco primale. I residui sono
elevati, ma finiti e decrescenti nel campione; non si impone che siano nulli.
La variazione di u è circa -.0414%, di rho -.0658%, di p -.0921%.
Non appaiono grandi oscillazioni nei sei campioni; non si certifica stabilità
a tempi lunghi o ad altri CFL.

Minimi globali finali rho=.999341826, p=.999078663; tutti gli stati e RHS
sono finiti, nessun passo rifiutato. Massimo slip sui **soli nodi forti**
<=6.25e-15. Il TE non è incluso in questa misura: il suo prodotto con la
normale aggregata uscente (-1,0) vale circa -.4, intenzionalmente non nullo.

### Mezze facce: controllo indipendente della chiusura

Ad ogni campione il test ricostruisce separatamente le due tracce, le riflette
con le rispettive normali e ricalcola entrambe le D-minus. La loro somma
coincide con SpatialResult::boundary al TE; il RHS coincide con
-(faces+volume+somma_Dminus)/V, senza reazione aggiunta. Non si è copiato
un flusso conservativo nel residuo primitivo.

| Quantità per le due mezze facce | Inizio | Fine |
|---|---:|---:|
| u_L dot n_1 | -.05718944773 | -.05716304975 |
| u_L dot n_2 | -.05718944765 | -.05716304967 |
| norma D-minus, faccia 1 | 1.5082593e-5 | 1.5065752e-5 |
| norma D-minus, faccia 2 | 1.5082593e-5 | 1.5065752e-5 |

Entrambe le facce contribuiscono e mantengono geometria distinta. Questi
prodotti delle tracce non sono imposti nulli: chiusura debole non significa
impermeabilità forte del nodo o di ciascuna traccia interna.

Le variazioni globali finali, rispetto allo stato discreto **dopo** le
proiezioni iniziali, sono Delta M=-6.0072e-10 e Delta E=-9.3755e-9.
Sono solo diagnostiche: dominio aperto con farfield, nessuna sottrazione
dei flussi esterni, somme globali soggette a cancellazione. Non dimostrano
conservazione o consistenza della chiusura.

## Confronto nativo e regressioni

Rieseguito l'audit 13B: la BC nativa proietta Solution_Old e residuo del
momento lungo la tangente alla normale aggregata. Nel suo aggiornamento
controllato a residuo zero u diventa zero, mentre rho/E restano fissate e
la pressione ricavata diventa 1.032. Il primitivo debole evolve invece tutte
le primitive tramite il RHS completo. Sono condizioni nodali diverse;
l'audit nativo a residuo zero **non** è un confronto di avanzamenti fisici
a parità di operatore, né un oracolo di correttezza. Non è stato aggiunto
un run NACA nativo a residuo fisico assemblato.

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase13c-build
/private/tmp/primitive-phase13c-build/test_driver '[primitive-weak-te]'
/private/tmp/primitive-phase13c-build/test_driver '[primitive-trailing-edge],[primitive-native-te-audit],[primitive-boundary],[primitive-wall-gradients],Constant primitive fields*,Affine NODE*,Forward Euler applies*,Rejected time steps*'
git diff --check
```

Comandi dalla radice. Build riuscita; prototipo **1 caso/314869 asserzioni**;
regressioni **21 casi/14087 asserzioni**, tutti passati. Nessuna campagna
completa MMS/vortice/bilanci. CSV: righe WEAK_TE senza prefisso.
Mesh originale invariata; nessun commit/push.

**A. Prototipo promettente**, limitatamente all'evoluzione iniziale osservata:
nessuna stagnazione imposta, positività e vincoli forti rimanenti preservati.
La NACA è supportata **solo sperimentalmente via API esplicita**; resta
rifiutata dal driver ordinario. Non sono validati Kutta, circolazione,
impermeabilità debole consistente, accuratezza o convergenza stazionaria.
