# FASE 11B — Gradienti Euler senza correzione di simmetria piana

## Politica e file

Quattro file produttivi modificati:

- `correctGradientsSymmetry.hpp`: selettore opzionale `correctEulerWalls=true`;
  SYMMETRY_PLANE è sempre selezionato, EULER_WALL solo se il selettore è true.
- `computeGradientsLeastSquares.hpp`: propagazione del parametro, con default
  true, attraverso wrapper 2D/3D e algoritmo LS/WLS.
- `computeGradientsGreenGauss.hpp`: identica propagazione per Green–Gauss.
- `CPrimitiveEulerSolver.cpp`, `PrepareReconstruction`: passa esplicitamente
  false, mantenendo idxVel=1 e tutti gli altri argomenti.

Le chiamate native non sono modificate e conservano il default true. Non si
disabilita la correzione sulle vere simmetrie, non si riclassificano marker e
non si introducono opzioni di configurazione. Le formule di correzione sono
invariate; cambia soltanto la selezione del marker nel percorso primitivo.
Kernel, fluttuazioni, MUSCL affine, limiter, positivity scaling, FE, normali,
stati riflessi e proiettore nodale non cambiano.

Esteso `primitive_euler_boundary_tests.cpp` con due test
`[primitive-wall-gradients]`. Nuovi artefatti: questo rapporto e
`primitive_euler_phase11b_curved_wall.csv` (51 righe).
Le modifiche e i dati FASE 11 già presenti nel working tree sono conservati;
il suo CSV resta la baseline storica, senza sovrascrittura.

## Gradienti: evidenza analitica e default nativo

Sul disco rotante, LS/WLS recuperano ai nodi di parete la matrice affine
G=((0,-.2),(.2,0)): massimo errore Frobenius <=2.07e-16.
Prima era .2828427125 indipendentemente da h.

Per p=1+.02r², p_n=.04 sul cerchio; il gradiente quadratico non viene
assunto esatto. Il massimo errore normale misurato è:

| Segmenti | Prima | Dopo LS | Dopo WLS |
|---:|---:|---:|---:|
| 12 | .04 | .0149434 | .0155043 |
| 24 | .04 | .00743183 | .00745831 |
| 48 | .04 | .00367819 | .00361676 |

La derivata normale ricostruita è non nulla in ogni nodo wall e si avvicina
al riferimento; non si impone esattezza a precisione macchina.

Su una parete rettilinea si usa u=.4+.3y, v=0, rho=1+.1y, p=1+.2y:
lo stato soddisfa slip su y=0, ma du/dn è non nullo. LS/WLS primitivi
mantengono du/dy=.3. Green–Gauss mantiene una derivata non nulla, senza
imporre esattezza affine del proprio stencil di bordo.

Per ciascuno dei tre metodi, il test confronta il gradiente del solver con
la chiamata generica avente false; ripete con MARKER_SYM, verificando che
la correzione resti attiva. Infine omette il nuovo argomento, come nei call
site nativi invariati: du/dy e dp/dy vengono annullati su entrambi i marker,
in accordo con il comportamento precedente. È un controllo mirato delle
API condivise, non una campagna completa di regressione dei solver nativi.

## Difetto spaziale prima/dopo

Identiche mesh FASE 11 (12/24/48 segmenti, 19/61/217 nodi), stato iniziale,
NODE/AGGREGATED/NODE_POINT, WLS, limiter NONE, nessuna sorgente.
Si usa il RHS effettivo **dopo** la reazione nodale. Norme volumetriche
normalizzate nel volume di ciascuna regione; E=sqrt(sum_k E_k²).
Ordini calcolati con h massimo lato triangolare: .619657/.337063/.174919.

| Regione | Segmenti | E prima | E dopo | ordine prima | ordine dopo |
|---|---:|---:|---:|---:|---:|
| globale | 12 | .0767265 | .0283553 | — | — |
| globale | 24 | .0440717 | .00953475 | .911 | 1.790 |
| globale | 48 | .0291515 | .00327148 | .630 | 1.631 |
| parete | 12 | .111899 | .0415228 | — | — |
| parete | 24 | .0822189 | .0185692 | .506 | 1.322 |
| parete | 48 | .0710655 | .00874409 | .222 | 1.148 |
| interno | 12 | .0304324 | .0107710 | — | — |
| interno | 24 | .0202938 | .00320169 | .665 | 1.992 |
| interno | 48 | .0152514 | .00102010 | .435 | 1.744 |

Massimo slip nodale iniziale prima/dopo <=1.47e-16; velocità tangenziale .2,
quindi non artificialmente annullata. Il massimo |u_L·n_f| diventa:

| Segmenti | Prima | Dopo |
|---:|---:|---:|
| 12 | .0517638 | .0258819 |
| 24 | .0261052 | .0130526 |
| 48 | .0130806 | .00654031 |

La nuova traccia coincide entro roundoff con il campo affine esatto nel punto
della mezza faccia: Omega*L/4=.1*sin(pi/segmenti). Il suo valore non nullo
è un difetto geometrico rispetto alla normale del segmento, distinto dallo
slip nodale. Non è stato eliminato o compensato nel kernel.

## Forward Euler e regressioni

Ripetuta soltanto l'evoluzione breve FASE 11, T=.1, con dt e dt/2;
**non eseguite le 250 iterazioni stazionarie**.

| Segmenti | dt | Errore soluzione a dt/2 | RHS finale a dt/2 | differenza campi / errore |
|---:|---:|---:|---:|---:|
| 12 | .00270270 | .00275716 | .0268861 | .0802% |
| 24 | .00138889 | .000913434 | .00884219 | .0708% |

Nessun passo rifiutato; stati e residui finiti. Ai checkpoint rho>=.999298,
p>=.999958 e slip <=1.47e-16. Piccola sensibilità a dt nel campione,
non una prova di stabilità incondizionata o di convergenza al limite temporale.

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase11b-build
/private/tmp/primitive-phase11b-build/test_driver '[primitive-wall-gradients],Rotating disk measures*,Rotating disk separates*'
/private/tmp/primitive-phase11b-build/test_driver '[primitive-boundary],[primitive-geometry],Constant primitive fields*,Affine NODE*,Forward Euler applies*,Rejected time steps*'
git diff --check
```

Build riuscita. Verifiche mirate: **4 casi, 20808 asserzioni**; regressioni:
**23 casi, 7782 asserzioni**, tutti passati. Nessuna campagna completa MMS,
vortice o FASE 9. Il CSV contiene le righe CURVED_WALL senza prefisso;
stessa interpretazione delle colonne della FASE 11, inclusa dt_difference.
I sei risultati LS/WLS sono inoltre stampati come WALL_GRADIENT dal test.

## Valutazione e limiti

1. **Ricostruzione:** recupero affine LS/WLS verificato; gradiente quadratico
   della pressione convergente. Correzione sulle vere simmetrie preservata.
2. **Residuo sulla parete poligonale:** migliorato e decrescente, ma le tracce
   esatte continuano a non soddisfare lo slip su ciascun segmento. Eventuali
   effetti delle fluttuazioni locali richiedono un esperimento separato:
   questa fase non ne dimostra un bug e non li corregge.
3. **Accuratezza globale:** ordini osservati 1.79 e 1.63; ordine due globale
   non dimostrato. Tre mesh non identificano il regime asintotico, né
   validano profili alari o comportamento stazionario a lungo termine.

Nessun commit, push o altra operazione Git di modifica dello stato.
