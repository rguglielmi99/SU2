# FASE 11 — Accuratezza su parete curva ed evoluzione stazionaria

## Perimetro e riferimento analitico

Modificato soltanto `UnitTests/SU2_CFD/primitive_euler_boundary_tests.cpp`:
tre test `[primitive-curved-wall]`, fixture SU2 e API produttive esistenti.
Nuovi artefatti: questo rapporto e `primitive_euler_phase11_curved_wall.csv`
(72 righe). **Nessun file produttivo modificato**, nessuna accelerazione,
nuova BC, sorgente o confronto con un solver conservativo usato come oracolo.

Con rho=1, u=-Omega*y, v=Omega*x, p=1+Omega²*r²/2 e Omega=.2:
div(u)=0, u·grad(rho)=0, u·grad(p)=0;

    (u·grad)u = -Omega² (x,y),    grad(p)/rho = Omega² (x,y).

Le quattro equazioni primitive hanno quindi prodotto volumetrico nullo.
Il test verifica questa identità con derivate analitiche in ogni nodo,
senza differenze numeriche né forcing.

## Mesh e norme

Fan iniziale di sei triangoli nel disco unitario; raffinamenti uniformi 1→4,
con proiezione sul cerchio dei soli nuovi vertici di bordo. Un unico
MARKER_EULER, nessun altro marker. Triangoli conformi, orientati e non degeneri.

| Segmenti | Nodi | Triangoli | h = massimo lato triangolare |
|---:|---:|---:|---:|
| 12 | 19 | 24 | .619657 |
| 24 | 61 | 96 | .337063 |
| 48 | 217 | 384 | .174919 |

Configurazione: NODE, AGGREGATED, NODE_POINT, MUSCL affine, WLS, limiter NONE,
gamma=1.4. Nessuna proiezione dell'inizializzazione analitica: è già tangente
alla normale nodale risultante, radiale su questi poligoni regolari.

Si misura **SpatialResult::rhs dopo la reazione nodale**, non il residuo
grezzo: E_k=sqrt(sum V_i rhs_ik² / sum V_i). Ogni regione usa il proprio
volume totale; `interior` include tutti i nodi non-wall, non solo quelli profondi.
E=sqrt(sum_k E_k²). Gli ordini sono log(E_old/E_new)/log(h_old/h_new),
non log2: h non dimezza esattamente.

## A — Difetto spaziale

Per simmetria E_u=E_v entro roundoff; la tabella li raggruppa. CSV con tutte
le componenti, errori, regioni, slip, tracce e diagnostiche dei gradienti.

| Regione | Segmenti | E_rho | E_u=E_v | E_p | E | ordine E |
|---|---:|---:|---:|---:|---:|---:|
| globale | 12 | .0175451 | .0494089 | .0263939 | .0767265 | — |
| globale | 24 | .0143253 | .0251674 | .0216864 | .0440717 | .911 |
| globale | 48 | .0112341 | .0147619 | .0169639 | .0291515 | .630 |
| parete | 12 | .0204470 | .0746686 | .0308632 | .111899 | — |
| parete | 24 | .0243205 | .0490939 | .0367158 | .0822189 | .506 |
| parete | 48 | .0261821 | .0374679 | .0394603 | .0710655 | .222 |
| interno | 12 | .0150131 | .00988569 | .0224786 | .0304324 | — |
| interno | 24 | .00908317 | .00830298 | .0138367 | .0202938 | .665 |
| interno | 48 | .00663508 | .00660751 | .0100630 | .0152514 | .435 |

Ordini componenti (rho,u/v,p), 12→24 / 24→48:
globale (.333,1.108,.323)/(.371,.813,.374);
parete (-.285,.689,-.285)/(-.112,.412,-.110);
interno (.825,.287,.797)/(.479,.348,.485).
**Il difetto globale scende, ma densità e pressione al bordo peggiorano
leggermente. Non si dimostra convergenza puntuale al bordo né ordine due.**

Slip nodale iniziale massimo 1.43e-16; velocità tangenziale .2 su tutte le mesh.
Massimo |u_L·n_f| ricostruito: .0517638, .0261052, .0130806. È non nullo e
scende circa linearmente con la lunghezza del segmento. Lo slip nodale non
è impermeabilità di ogni mezza faccia.
Per la traccia analitica nel punto della stessa mezza faccia il massimo è
Omega*L/4=.1*sin(pi/segmenti): .0258819, .0130526, .00654031, metà del valore
ricostruito. La geometria introduce quindi già un difetto non nullo; la
ricostruzione effettiva lo amplifica in questo caso.

### Incompatibilità verificata nel percorso dei gradienti

`CPrimitiveEulerSolver::PrepareReconstruction` richiama il LS nativo con
indice velocità 1 (`CPrimitiveEulerSolver.cpp:99–112`). Questo invoca
`correctGradientsSymmetry` (`computeGradientsLeastSquares.hpp:318–320`),
che tratta anche EULER_WALL, elimina il gradiente normale degli scalari e
gli accoppiamenti normale/tangente delle velocità
(`correctGradientsSymmetry.hpp:59–89,93–115`).

Nel riferimento circolare p_n=.04 e la matrice della velocità è antisimmetrica.
La correzione locale da parete piana annulla il primo e rimuove la seconda:
il test misura, su tutte le mesh, massimo errore Frobenius della matrice
velocità .2828427125 e massimo errore di p_n .04 ai nodi di parete.
Questi errori **non diminuiscono**. È una incompatibilità matematica
verificata con questo equilibrio curvo, non una prova di bug universale
nei solver nativi. Non quantifica da sola tutti i contributi al RHS.
Non è stata aggirata o corretta: serve una revisione mirata del contratto
dei gradienti su pareti curve prima di rivendicare accuratezza generale.

## B — Evoluzione breve e sensibilità temporale

Tempo finale comune T=.1. Passi costanti ricavati dal CFL iniziale .05 e
aggiustati per raggiungere esattamente T; secondo run con dt/2.

| Segmenti | dt | E_sol(T), dt/2 | RHS(T), dt/2 | L2(W_dt/2-W_dt) / E_sol |
|---:|---:|---:|---:|---:|
| 12 | .00270270 | .00696120 | .0637656 | .146% |
| 24 | .00138889 | .00350001 | .0290727 | .161% |

Il confronto usa la differenza effettiva fra i due campi, non solo fra norme.
Il drift diminuisce col raffinamento e la sensibilità temporale è piccola
rispetto al drift; due passi non costituiscono una prova completa del limite
dt→0. La soluzione campionata non è un punto fisso discreto.

## C — Tentativo di equilibrio discreto

Mesh 12 segmenti, 250 passi globali dt=.00276177169 (CFL iniziale .05),
T=.690443. È lo stesso operatore FE autonomo, interpretato come rilassamento:
nessun pseudo-tempo locale o metodo temporale aggiuntivo.

| Iterazione | 0 | 50 | 100 | 150 | 200 | 250 |
|---|---:|---:|---:|---:|---:|---:|
| RHS L2 | .0767265 | .0603108 | .0513722 | .0456710 | .0417018 | .0386076 |

Rapporto finale/iniziale **.503185**. Diminuisce ai checkpoint osservati,
ma non è raggiunto un equilibrio; errore dalla soluzione iniziale .0343570
e velocità tangenziale minima .145641 contro .2 iniziale. Non emerge ancora
stabilizzazione delle variabili. FE è ammissibile nel breve intervallo
campionato, ma non dimostrato efficiente per convergere allo stazionario.
Il decadimento può accompagnare diffusione del moto, non certifica accuratezza.

In tutti i checkpoint A–C: slip <=1.47e-16, rho>=.996028, p>=.999914;
nessun passo rifiutato. Nessuna proprietà di stabilità incondizionata dedotta.

## Riproducibilità, regressioni e seguito

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase11-build
/private/tmp/primitive-phase11-build/test_driver '[primitive-curved-wall]'
/private/tmp/primitive-phase11-build/test_driver '[primitive-boundary],[primitive-geometry],Constant primitive fields*,Affine NODE*,Forward Euler applies*,Rejected time steps*'
git diff --check
```

Build riuscita; nuovi test: **3 casi, 12599 asserzioni**; regressioni ridotte:
**23 casi, 7782 asserzioni**, tutti passati. Nessuna campagna completa
MMS/vortice/FASE 9. Il CSV riproduce le righe `CURVED_WALL,` del primo run,
privandole del prefisso; `dt_difference` contiene nelle colonne `error_*`
la differenza fra campi temporali, le altre righe l'errore sul riferimento
iniziale. Le diagnostiche dei gradienti sono relative a quest'ultimo, anche
dopo l'evoluzione. `rhs_total` è sempre globale e ripetuto nelle tre regioni.

La fase diagnostica è completata, **la validazione di accuratezza della
parete curva no**. Passo successivo raccomandato: isolare la correzione dei
gradienti e distinguere simmetria piana da slip curvo, conservando kernel e
normali geometriche. Richiede una decisione separata; nessuna modifica
produttiva o accelerazione stazionaria applicata. Non è validato un profilo
alare, né un limite asintotico con queste sole tre mesh. Nessun commit o push.
