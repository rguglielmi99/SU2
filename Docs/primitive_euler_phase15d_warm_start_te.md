# FASE 15D — Warm start nativo, politica TE e ultimo tentativo stazionario

## Decisione

**B — Risultato parziale, non una validazione stazionaria del NACA 0012.**
La combinazione più promettente è il campo convergente EULER nativo come
inizializzazione, TE PROJECTED e pseudo-tempo locale Forward Euler.
A CFL=.05, dopo 9000 aggiornamenti autonomi, il RHS globale passa da
.001060758 a 2.771767e-5; quello interno da .000550626 a 2.848876e-5.
Il TE conserva RHS=.0794207 e i coefficienti variano ancora. Non è stato
generato un confronto quantitativo Cp come fra due soluzioni stazionarie.

Si raccomanda di **chiudere la FASE 15 con questa limitazione documentata**,
senza sviluppare altre accelerazioni per inseguire lo steady su questo caso.
Questo risultato non invalida automaticamente le verifiche instazionarie
precedenti, né certifica correttezza fisica della politica PROJECTED.

## Perimetro e implementazione

Mesh originale, Mach=.2, AOA=0°, gamma=1.4, rho_inf=p_inf=R=1,
u_inf=.2*sqrt(1.4), v_inf=0, REF_AREA=1. NODE, AGGREGATED, NODE_POINT,
MUSCL affine, WLS, limiter NONE e farfield LOCAL_CHARACTERISTIC invariati.
Nessun Kutta, correzione di shock, modifica al kernel, alle fluttuazioni,
alle facce interne, alla ricostruzione o alla geometria. Nessun AOA=2°.
Il solver EULER nativo non è stato modificato o rieseguito.

Modifiche produttive specifiche 15D:

- `Common/include/CConfig.hpp`, `Common/src/CConfig.cpp`: opzioni
  `PRIMITIVE_TE_POLICY=WEAK|PROJECTED` (default WEAK) e
  `PRIMITIVE_INITIAL_FILE`, riservate al primitivo con validazione esplicita.
- `SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp` e relativo `.cpp`:
  inizializzazione FILE tramite lettore binario SU2, proiezione iniziale
  esplicita e vincolo aggregato al solo TE autorizzato per PROJECTED.
- `SU2_CFD/src/iteration/CPrimitiveEulerIteration.cpp`: inizializzazione FILE
  una sola volta, proiezione dichiarata e diagnostica iniziale prima del
  primo passo. Il file non viene riletto durante l'evoluzione.

`AdvanceExplicitEuler` e `AdvanceLocalPseudoTime` non sono stati modificati
in questa fase. Nessun nuovo integratore, CFL adattivo, MG o schema implicito.
Il guard di inizializzazione impedisce di riproiettare il campo dopo passi
accettati. La modalità WEAK rimane disponibile e predefinita.
`PRIMITIVE_WEAK_TE_GLOBAL_ID` resta il nome storico dell'autorizzazione
topologica, anche quando la politica esplicita è PROJECTED.
Nessun ID NACA è hardcoded nel codice produttivo.

## Importazione Float64 e identità dei nodi

Nuovo convertitore `UnitTests/convert_primitive_euler_initial.py`, applicato
a `results/native/native_state_000599.dat` della FASE 15A.
Si usano rho, momentum_x/y ed Energy **Float64**, non i VTU Float32:

    u = momentum_x/rho, v = momentum_y/rho
    p = (gamma-1)*(Energy - (momentum_x²+momentum_y²)/(2*rho))

Nel formato binario nativo l'ID non è una colonna: l'ordine delle righe è
l'ordine GlobalIndex del writer SU2. Il convertitore verifica tale identità
su tutti i 5233 nodi contro le coordinate della mesh e i suoi ID espliciti,
senza ricerca del vicino più prossimo o interpolazione. Nel file convertito
aggiunge GlobalIndex esplicito e il caricatore lo ricontrolla nodo per nodo.

- Copertura completa/univoca degli ID e coordinate: scarto massimo **0.0**.
- Inversa W→U: massimo scarto **8.881784197e-16**.
- Prima della proiezione: rho_min=.98215893, p_min=.98841868, finiti/positivi.
- Gamma, freestream e REF_DIMENSIONALIZATION=DIMENSIONAL verificati dalla
  configurazione archiviata del nativo; il driver primitivo usa riferimenti
  unitari. Il restart originale rimane intatto.
- SHA256 restart originale:
  `cbb014b7a30800f338b690daa8388bf603c4d551ac9a862d554f0cefe3603857`.
- SHA256 mesh:
  `9094b51c2628d3bb4c865d774b2308dbb59a213ad19c007cfcd7d57e7aaeebeb`.

Il file convertito ha schema x,y,Density,Velocity_x,Velocity_y,Pressure,
GlobalIndex e formato binario SU2. `InitializePrimitiveFile` riutilizza
`CSolver::Read_SU2_Restart_Binary`: controlla header, dimensione e payload
**prima** che il lettore possa interpolare una mesh diversa, poi schema,
ID, coordinate e ammissibilità di tutti gli stati prima di scriverne uno.
Il lettore nativo impone un nome inferiore a 100 caratteri: si usa un percorso
relativo corto e si rifiutano percorsi troppo lunghi prima della chiamata.

Questa è un'inizializzazione di campo, **non una validazione del restart di
prosecuzione**. RESTART_SOL resta disabilitato; tempo/counter iniziano da zero.
Ogni prova parte dal freestream o dallo stesso restart nativo originale.

### Proiezioni misurate

Il driver proietta esplicitamente le sole velocità dei nodi Euler forti
all'inizializzazione, mantenendo rho/p. Sul campo nativo importato, per
entrambe le politiche:

| Regione | Nodi | Max correzione velocità | Max correzione rho/p |
|---|---:|---:|---:|
| Parete ordinaria | 199 | 1.062277351e-10 | 0 |
| TE | 1 | 0 | 0 |
| Altri | 5033 | 0 | 0 |

Le differenze sono misurate sui restart Float64 dopo l'inizializzazione,
non dedotte da VTU a precisione ridotta. Nel freestream PROJECTED il TE
ha invece una correzione esplicita della velocità, assente in WEAK.

## Significato della politica TE

La validazione topologica preesistente richiede un TE a due rami Euler
coerenti. Si usa la somma delle due normali integrate reali, senza sostituire
normali alle facce:

    N = S_1 + S_2, n = N/|N|, P = I - n*n^T
    velocità iniziale PROJECTED = P*velocità importata
    RHS_vel PROJECTED = P*RHS_vel grezzo
    boundary_vel += V*n*(n·RHS_vel grezzo)

I due D-minus restano presenti e invariati; l'ultima riga contabilizza la
reazione di vincolo nel residuo integrato. rho/p e loro RHS non sono
proiettati. Il medesimo meccanismo di reazione preesistente resta invariato
sui nodi ordinari. In WEAK il TE non riceve né proiezione iniziale né reazione
nodale; restano le riflessioni delle due mezze facce e il suo RHS libero.

PROJECTED **non equivale** alla BC conservativa nativa: proiettare velocità
a rho,p fissi cambia l'energia conservativa di
`-rho*(velocità·n)^2/2`. Proiettare invece il momento a energia totale fissa
può aumentare p. Sul TE simmetrico di questo caso la tangente aggregata
ammette solo la componente trasversale: ciò può favorire una condizione
vicina alla stagnazione, senza imporre Kutta né dimostrarne l'ottimalità.

Il calo immediato del RHS TE passando da WEAK a PROJECTED è in parte
**costruito dal diverso vincolo**, non è una prova di convergenza. Si misura
quindi anche il decadimento successivo del proprio RHS e del residuo interno.

## Matrice diagnostica prima del primo aggiornamento

Le history all'indice 0 precedono qualsiasi passo. Norme come in 15C:
`E_G² = sum_G V_i*sum_k RHS_ik² / sum_G V_i`, con il RHS spaziale completo
corrente dopo la reazione. Regioni disgiunte: interno 4983, parete forte 199,
TE 1, farfield 50. Sono norme diagnostiche di primitive adimensionali,
non una norma di energia. Farfield iniziale nullo nel freestream: nessuna
divisione per zero per normalizzarlo.

| Inizializzazione/TE | Eglob | Eint | Ewall | ETE | Efar | CL | CD |
|---|---:|---:|---:|---:|---:|---:|---:|
| Freestream/WEAK, riferimento 15C | .0470270 | .0260361 | 13.1941 | 140.245 | 0 | 0 | 0 |
| Freestream/PROJECTED | .0482564 | .0277444 | 13.3367 | 144.273 | 0 | 0 | 0 |
| Nativo/WEAK | .00276515 | .000550626 | .300462 | 290.260 | 8.66769e-7 | -6.73956e-5 | .000305180 |
| Nativo/PROJECTED | .00106076 | .000550626 | .300462 | 16.9857 | 8.66769e-7 | -6.73956e-5 | .000305180 |

Sul campo nativo iniziale, Mmax nodale=.23344177 e tracce=.23380512.
Le due inizializzazioni native coincidono anche al TE. Cambia quindi la
politica di reazione, non l'importazione del campo. In WEAK il singolo TE
pesa 85.58% nel quadrato del RHS globale; in PROJECTED pesa 1.99%, la parete
ordinaria 73.89%, l'interno 24.12%. Il campo nativo con residuo conservativo
molto piccolo **non è un punto fisso del primitivo**. Non si deduce un bug
del nativo: i due operatori e le BC nodali sono differenti a mesh finita.

## Prove brevi, CFL=.05, 300 passi

Freestream/WEAK riutilizza i risultati 15C, senza un'altra simulazione.

| Inizializzazione/TE | Eglob finale | Eint finale | Ewall finale | ETE finale | CL finale | CD finale | Costo s |
|---|---:|---:|---:|---:|---:|---:|---:|
| Freestream/WEAK | .0275565 | .0283484 | 2.08630 | 10.9326 | 3.05303e-6 | .683667 | riusato |
| Freestream/PROJECTED | .0276592 | .0284110 | 2.14322 | 60.0774 | 3.05495e-6 | .682793 | 11.67 |
| Nativo/WEAK | .000833452 | .000745532 | .114661 | 31.3011 | -6.72431e-5 | -.00319828 | 11.24 |
| Nativo/PROJECTED | .000337929 | .000277126 | .0696564 | 3.15552 | -6.72437e-5 | -.00348258 | 10.22 |

Nessun rifiuto e tutti gli stati subsonici/positivi. Minimi rho/p osservati:
freestream/PROJECTED .955014/.937532, nativo/WEAK .981987/.988103,
nativo/PROJECTED .980985/.988103. Massime tracce rispettivamente
.323090, .233805, .233805. Il cambio TE da solo non migliora il transitorio
freestream: l'interno cresce in entrambi i casi. Nel warm start, PROJECTED
dimezza invece il residuo interno, mentre WEAK lo aumenta di circa il 35%.
Questo motiva la selezione di nativo/PROJECTED, non quattro campagne lunghe.

## Prosecuzione selezionata e limiti

Una prova lunga nativo/PROJECTED, CFL=.05, 9000 passi indipendenti dal campo
nativo originale. Costo 312.71 s, circa 34.75 ms/passo, nessun rifiuto,
nessuna riduzione del CFL. Residuo ricalcolato sul campo finale, non norma
dell'aggiornamento. Riferimenti e min/max coprono l'intera history.

| Quantità | Iniziale | Finale |
|---|---:|---:|
| Eglob | .001060758 | 2.771767e-5 |
| Eint | .000550626 | 2.848876e-5 |
| Ewall | .300462 | .00211948 |
| ETE | 16.9857 | .0794207 |
| Efar | 8.66769e-7 | 1.53970e-7 |
| CL | -6.73956e-5 | -3.51653e-5 |
| CD | .000305180 | -.003213625 |

Decadi globali rispetto al proprio iniziale **1.583**, rispetto al freestream
WEAK **3.230**. Decadi interne rispetto al proprio iniziale circa 1.286,
rispetto al freestream circa 2.961. Range ultimi 200 passi: CL=6.52615e-5,
CD=5.14405e-5. Il TE ha ridotto il suo RHS di circa 2.33 decadi, ma conserva
un difetto locale e oscillazioni. Il CD negativo resta un limite fisico
serio del risultato; il residuo piccolo non ne certifica l'accuratezza.

rho_min=.98098544, p_min=.98806223, Mmax nodale=.23356904, tracce=.23380512.
Passi locali finali: minimo/TE=2.838873e-6, massimo=.02415325.
Il tempo fisico rimane zero: non si confrontano questi campi a un tempo
fisico comune. Nessuna soluzione nativa viene imposta dopo l'inizializzazione.

Checkpoint della prova .05 (tutti i componenti e le regioni sono nel CSV):

| Iterazione | Eglob | Eint | ETE | CL | CD |
|---:|---:|---:|---:|---:|---:|
| 0 | .001060758 | .000550626 | 16.9857 | -6.73956e-5 | .000305180 |
| 300 | .000337929 | .000277126 | 3.15552 | -6.72437e-5 | -.00348258 |
| 1000 | .000170767 | .000171241 | .546349 | -6.94749e-5 | -.00312939 |
| 2000 | 8.58331e-5 | 8.66919e-5 | .366533 | -6.50819e-5 | -.00310070 |
| 3000 | 5.69908e-5 | 5.81978e-5 | .289711 | -6.39735e-5 | -.00311243 |
| 5000 | 3.91126e-5 | 4.04931e-5 | .157664 | -4.52523e-5 | -.00318224 |
| 6000 | 3.46165e-5 | 3.53076e-5 | .151449 | -1.95287e-5 | -.00313146 |
| 9000 | 2.77177e-5 | 2.84888e-5 | .0794207 | -3.51653e-5 | -.00321363 |

TE finale: rho=1.00095226, u=0, v=-5.06132e-6, p=1.02039454.
RHS TE: (.0327508, 0, .00888436, -.0718060). La componente normale della
velocità è vincolata; gli altri residui restano misurati e non sono azzerati.
Il piccolo passo al TE non viene scambiato per convergenza.

Differenza finale dal campo nativo importato, stati Float64 e pesi dei
volumi VTU (Float32):

| Componente | Max differenza assoluta | RMS pesata sul volume globale |
|---|---:|---:|
| rho | .00721139 | 9.56131e-7 |
| u | .0262564 | 4.05901e-6 |
| v | .00431533 | 2.00554e-6 |
| p | .00794167 | 1.28532e-6 |

La norma globale piccola è influenzata dall'ampio dominio esterno; i massimi
mostrano cambiamenti locali non trascurabili. Il campo evolve autonomamente,
non è mantenuto artificialmente uguale a quello nativo.

## Controllo CFL finale e budget

Dopo la prova lunga regolare, un solo controllo aggiuntivo nativo/PROJECTED
a CFL=.1, sempre dallo stesso file originale, per 3000 passi. Non è una
prosecuzione e non è stato usato un CFL infinitesimo per ridurre gli update.

| CFL / passi | Eglob finale | Eint finale | ETE finale | Costo s | ms/passo | Decadi proprie/s |
|---|---:|---:|---:|---:|---:|---:|
| .05 / 9000 | 2.77177e-5 | 2.84888e-5 | .0794207 | 312.71 | 34.75 | .00506 |
| .1 / 3000 | 3.56023e-5 | 3.61683e-5 | .150539 | 97.43 | 32.48 | .01513 |

A .1: CL=-1.45871e-5, CD=-.00313286; range ultimi 200 passi
CL=.000240314 e CD=2.58013e-5. Minimi rho/p=.980974/.988076, massimo Mach
nodale=.233464 e ricostruito=.233805, zero rifiuti. Passi finali
5.67462e-6… .0483065. Nessun arresto anticipato, tutte le prove raggiungono
il proprio limite di iterazioni.

Per confronto, .05 a 6000 passi ha Eglob=3.46165e-5, Eint=3.53076e-5 e
ETE=.151449, vicini a .1/3000. Il raddoppio del CFL accelera il rilassamento
iniziale del warm start, mentre in 15C il freestream mostrava risalite.
`CFL*numero_passi` è solo un confronto indicativo: i passi sono locali,
dipendono dallo stato e non definiscono un tempo fisico comune. Le decadi/s
su intervalli diversi non costituiscono un fattore di accelerazione
asintotico o una misura di accuratezza. Il controllo .1 non risolve la
variazione dei coefficienti o il difetto locale TE.

Budget configurabile **480 s**, spesi **443.28 s** complessivi nei tre pilot
e nelle due prove successive, incluso output e sorveglianza. Baseline 15C
riutilizzata, compilazione/test/postprocessing separati. Nessuna necessità
documentata di ridurre ulteriormente il CFL: i warm start sono positivi e
in decadimento globale. Nessun adattamento automatico o SSPRK3 introdotto.

Il filtro prudente di postprocessing richiede Eglob ed Eint sotto 1e-3 dei
rispettivi iniziali freestream WEAK, ETE<1e-3, range CL/CD<1e-5 sugli ultimi
200 passi e snapshot finale presente. Non è una soglia universale di
accuratezza. Nessuna prova lo supera: non si dichiara A soltanto perché il
warm start beneficia di un residuo iniziale minore. Anche la prova migliore
rimane **B**: forte miglioramento e campo ammissibile, ma residuo locale,
oscillazioni dei coefficienti e CD negativo impediscono una validazione.

## Localizzazione del residuo

L'analisi dei VTU 15C è soltanto postprocessing. Si separano i nodi interni
entro distanza euclidea .1 dal TE (1,0), entro .1 dal LE (0,0), e il resto.
Queste sono regioni diagnostiche geometriche, non cambiamenti di stencil.
Quote del **quadrato** della norma globale:

| Campo | Interno vicino TE | Interno vicino LE | Altro interno | Parete incluso TE | Farfield |
|---|---:|---:|---:|---:|---:|
| 15C, .05/3000 | 13.53% | 9.09% | 75.10% | 2.28% | .00016% |
| 15C, .1/3000 | 60.18% | 9.31% | 26.26% | 4.07% | .1846% |
| 15C, .2/6000 | 28.12% | 22.56% | 39.68% | 9.64% | .0068% |
| 15D, nativo/PROJECTED, .05/9000 | 67.16% | 17.20% | 10.19% | 5.45% | .00032% |

Nell'ultimo campo l'interno totale pesa 94.55%; il **singolo** nodo TE solo
.0638%. Il residuo si concentra quindi soprattutto nell'interno vicino al
TE, ma esiste anche una componente vicino al LE e nel resto del dominio.
La correlazione con il TE non prova che la sua BC sia l'unica causa. Non si
può diagnosticare definitivamente, con questi run, un limite temporale FE,
un difetto spaziale o la dinamica accoppiata delle pareti/TE/interno.

## Verifiche, artefatti e riproducibilità

Nuovi test `[primitive-te-policy]` e `[primitive-initial-file]`:
**2 casi, 299 asserzioni**, passati. Verificano default WEAK, proiezione
PROJECTED e scalari invariati, facce/volume e nodi ordinari invariati sullo
stesso campo, contabilizzazione della reazione, aggiornamento sincrono
locale, rollback di candidato negativo, inizializzazione completa e rifiuto
transazionale di ID duplicati, coordinate errate, pressione negativa e
dimensione errata. Regresse mirate: **11 casi, 317783 asserzioni**, passate
(WEAK/TE, pseudo-tempo, Forward Euler e passi rifiutati). CLI poligonale:
4 snapshot e tabelle pressione, CL/CD coerenti con le forze. Build riuscita.
Un primo link dei test, avviato prima della fine della ricompilazione delle
librerie, è stato ripetuto dopo la build; nessun errore di test residuo.

Radice separata `TestCases/primitive_euler_examples/naca0012_aoa0/results15d/`:
`conversion.json`, `native_initial.dat`, `summary.csv`, `checkpoints.csv`,
`baseline_reference.csv`, `budget.json`, history/log/config e restart di ogni
run. `initial_projection.csv` misura le correzioni; `field_difference_native.csv`
misura la deriva in Float64 rispetto al campo nativo (non errore rispetto a
una soluzione esatta), usando volumi VTU per le norme pesate.

In ogni cartella `run_*/`, `iterations.pvd` punta a VTU validati, con
Pseudo_Iteration, Local_Pseudo_Step, rho, velocità, p, Mach e quattro RHS
nodali. Snapshot iniziale, intermedi e finale; tempo PVD = indice di
pseudo-iterazione, nessun campo Physical_Time. XML, riferimenti e payload
VTU verificati, così come ID, conteggi, positività e norme ricostruite da
RHS/volumi (rtol=3e-6 per Float32). Non verificata l'apertura GUI in ParaView.

Grafici `RHS_CURRENT.png`, `RHS_INTERIOR.png`, `RHS_WALL_STRONG.png`,
`RHS_TE.png`, `RHS_normalized.png`, `CL.png`, `CD.png`. I CSV `*_localization.csv`
e `*_largest_residual_nodes.csv` permettono di localizzare i contributi
senza rieseguire simulazioni. Nessun confronto Cp quantitativo quando il
gate di stazionarietà fallisce: `cp_comparison_status.txt`.

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase15d-tests
/private/tmp/primitive-phase15d-tests/test_driver '[primitive-te-policy],[primitive-initial-file]'
/private/tmp/primitive-phase15d-tests/test_driver '[primitive-weak-te],[primitive-pseudo-time],[primitive-trailing-edge],Forward Euler applies*,Rejected time steps*'
python3 UnitTests/run_primitive_euler_steady15d.py --stage pilots
python3 UnitTests/run_primitive_euler_steady15d.py --stage final \
  --label final_native_PROJECTED_cfl005 --initial native --policy PROJECTED \
  --cfl .05 --steps 9000 --budget-seconds 480
python3 UnitTests/run_primitive_euler_steady15d.py --stage final \
  --label cfl_probe_native_PROJECTED_01 --initial native --policy PROJECTED \
  --cfl .1 --steps 3000 --budget-seconds 480
python3 UnitTests/run_primitive_euler_steady15d.py --stage postprocess
python3 UnitTests/check_primitive_euler_polygon_cli.py \
  --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/primitive-phase15d-polygon
git diff --check
```

Lo script rifiuta di sovrascrivere cartelle esistenti per nuovi run; usare
`--output-dir` diverso per una nuova campagna. La modalità postprocess non
esegue simulazioni. Il budget configurabile cumula i secondi dei processi
solver tra invocazioni; compilazione, test e postprocessing sono separati.
Arresti per budget, positività/Mach non ammissibili o quattro medie successive
su 200 passi crescenti tutte oltre 5%; i checkpoint accettati sono conservati.
Tutte le modifiche 15A–15D rimangono nel working tree. Nessun commit o altra
operazione Git di scrittura, nessuna nuova fase avviata.
