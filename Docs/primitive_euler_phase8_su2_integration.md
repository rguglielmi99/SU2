# FASE 8 — Integrazione nel driver standard SU2_CFD

## 8.1 Audit preliminare, registrato prima delle modifiche produttive

Stato iniziale: branch `feature/primitive-path-solver`, working tree pulito.
L'audit seguente deriva dal repository locale, non da un'architettura SU2
assunta. HEAD di partenza e finale: `4058eee641736a11bdc49aa8f24bad1e88abc935`.
I riferimenti di riga della tabella sono relativi a questo commit.

| Sottosistema | Punto di integrazione attuale verificato | Modifica minima prevista | Rischio |
|---|---|---|---|
| Config | `CConfig::SetConfig_Options`, opzioni PRIMITIVE a 1116–1119; postprocessing temporale a 3959–4030 | Modalità iniziale/bordo e controlli specifici prima che TIME_STEPPING azzeri MGLEVEL | Opzioni ignorate se validate troppo tardi |
| Factory solver | `CSolverFactory::CreateSolverContainer` 60–75; `CreateSubSolver` 334–339 crea già il vero primitivo | Conservare factory e metadati NONE per l'integration generica | Non applicare un'integrazione conservativa al layout W |
| Driver | `SU2_CFD.cpp` 90–144 seleziona `CSinglezoneDriver`; `CDriver` 120–125 blocca il primitivo | Sostituire il blocco con validazione fail-fast | Attivare solo seriale, zona unica, mesh statica |
| Geometria | `CDriver::InitializeGeometryFVM` 692–790 legge `.su2`, ordina RCM, crea edge/dual/marker | Nessuna modifica alle formule; verificare equivalenza per coordinate/ID globali | Ordine RCM diverso dal harness |
| Iteration | `CIterationFactory` 116–118 abortisce; `CSinglezoneDriver::Run` chiama `CIteration::Solve` | Aggiungere `CPrimitiveEulerIteration` che chiama una sola volta l'API FE esistente | Evitare doppi assemblaggi e aggiornamenti |
| Integration | `CIntegrationFactory` crea dall'integrationType; primitivo è NONE | Lasciare NONE: FE transazionale già incapsula spazio+tempo | La via generica chiama Source_Residual e BC incompatibili |
| Numerics | `CDriver::InitializeNumerics` 1474–1558 non attiva euler/compressible per PRIMITIVE | Nessun Roe conservativo da creare; contenitori nulli già previsti | Non sostituire il kernel path-conservative |
| Boundary | `EvaluateSpatialResidual` 181–185 ammette CUSTOM/FAR; callback a 221 | Adapter esterno uniforme o vortice analitico sui punti reali, tempo t_n | Forte conservativo non utilizzabile |
| Output | `COutputFactory` usa oggi COutput generico; `CFVMOutput` offre sorter/writer nativi | Aggiungere `CPrimitiveEulerOutput`, layout W esplicito | Etichette, tempo, normalizzazione del residuo |

### Percorso completo verificato

1. `SU2_CFD.cpp` crea `CConfig`, determina driver e chiama `StartSolver`.
2. `CDriver::PreprocessInput` crea la configurazione di zona; oggi segue il
   blocco esplicito del primitivo, prima della lettura della dimensione.
3. `PreprocessOutput` crea output e registra history/volume; poi
   `InitializeGeometryFVM` costruisce la geometria nativa da file.
4. `InitializeSolver` usa la factory e il costruttore primitivo; le variabili
   sono già `CPrimitiveEulerVariable`. `InitializeNumerics` alloca contenitori
   e non attiva un operatore conservativo per il primitivo.
5. `InitializeIntegration` usa i metadati NONE; `PreprocessIteration` richiede
   una physics iteration oggi assente. Il ciclo standard è
   Preprocess → Run/Solve → Postprocess → Update → Monitor → Output.
6. `CSinglezoneDriver::Preprocess` 114–137 imposta tempo e chiama
   SetInitialCondition; il primitivo eredita attualmente il no-op di base.
7. L'API `AdvanceExplicitEuler` 263–287 valuta lo spazio una volta,
   valida tutti i candidati, poi aggiorna W, residuale e tempo. È l'unica
   implementazione numerica da richiamare nell'adapter.
8. `COutput::SetHistoryOutput` 190–214 carica dati comuni e specifici;
   `LoadCommonHistoryData` 2282–2290 ricava un tempo dal contatore. L'output
   primitivo dovrà sovrascriverlo con il tempo realmente accettato.
9. `COutput::SetResultFiles` 829–870 usa sorter FVM e writer selezionati.
   PARAVIEW_XML e PARAVIEW_MULTIBLOCK sono implementati a 525–575;
   quest'ultimo riceve il tempo fisico della history.
10. `CSinglezoneDriver::Monitor` 264 e seguenti arresta per MAX_TIME,
    TIME_ITER o convergenza temporale. L'adapter non deve introdurre
    arresti impliciti basati su norme conservative.

### Piano minimo di implementazione

L'iterazione specializzata appartiene alla factory e al normale driver;
non è un secondo driver. Inizializza una volta, scrive lo stato iniziale
con output nativo, usa un solo FE per passo e lascia il normale monitoraggio
e output al driver. Il bordo resta una callback dell'operatore esistente.
Un output CFVM specifico registra W e norme del RHS primitivo pesate sui
volumi. Le formule del vortice saranno estratte in un helper condiviso,
verificando che il riferimento analitico precedente sia invariato.

Prima campioni uniformi e negativi, poi vortice pilota, equivalenza con
harness e infine N=8,16,32,64. Non si avvia una campagna completa prima
del controllo del primo passo. Nessuna modifica numerica è prevista.

## 8.2 Implementazione e perimetro verificato

La FASE 8 è completata nel perimetro seriale descritto sotto: l'eseguibile
ordinario `SU2_CFD` legge il cfg e la mesh, istanzia il vero solver e produce
output nativo. Non è stato creato un secondo driver. Le modifiche sono
nel working tree, senza nuovo commit.

File produttivi modificati o aggiunti, elenco completo:

| File | Ruolo |
|---|---|
| `Common/include/CConfig.hpp` | Campi/getter delle quattro nuove opzioni e validazione temporale specifica |
| `Common/src/CConfig.cpp` | Parsing e `ConfigurePrimitiveDriverTime`; cattura di MGLEVEL prima del suo azzeramento generico |
| `Common/include/verification/isentropic_vortex.hpp` | Profilo analitico condiviso, estratto senza cambiare le formule FASE 7 |
| `SU2_CFD/src/drivers/CDriver.cpp` | Rimozione del blocco, validazione specifica e diagnostica delle eccezioni |
| `SU2_CFD/src/drivers/CSinglezoneDriver.cpp` | Tempo fisico del solver accettato nel preprocessing primitivo |
| `SU2_CFD/include/iteration/CPrimitiveEulerIteration.hpp` | Iterazione di fisica registrata nel driver ordinario |
| `SU2_CFD/src/iteration/CPrimitiveEulerIteration.cpp` | Inizializzazione, marker, callback esterna, una chiamata FE per passo |
| `SU2_CFD/src/iteration/CIterationFactory.cpp` | Creazione dell'iterazione primitiva |
| `SU2_CFD/include/output/CPrimitiveEulerOutput.hpp` | Output FVM specializzato |
| `SU2_CFD/src/output/CPrimitiveEulerOutput.cpp` | Campi W, history, tempo e audit geometrico opzionale |
| `SU2_CFD/src/output/COutputFactory.cpp` | Creazione dell'output primitivo |
| `SU2_CFD/src/solvers/CSolverFactory.cpp` | Solo commento: integrationType NONE rimane intenzionale |
| `SU2_CFD/src/meson.build` | Registrazione delle due nuove unità di compilazione |

`CPrimitiveEulerSolver`, `CPrimitiveEulerVariable`, kernel, geometria,
ricostruzione, limiter e algoritmo temporale non hanno modifiche. La
physics iteration chiama `AdvanceExplicitEuler` già verificato; non
passa per residui/BC conservative o un `CIntegration` conservativo.
`EstimateTimeStep`, quando richiesto, esegue la stima esistente senza
commettere il residuale; soltanto FE accetta e aggiorna lo stato.

Supportati: singola zona, un rank MPI e un thread OpenMP, 2D triangolare
statico, gas ideale nelle condizioni già ammesse dal solver, NODE,
AGGREGATED, quadrature volumetriche esistenti, FE, MUSCL affine e
gradienti/limiter preesistenti. La campagna vortice usa WLS, NONE e
NODE_POINT. GG, LS, Barth–Jespersen e Venkatakrishnan hanno anche smoke
test uniformi CLI, oltre ai test numerici precedenti.

Rifiutati: implicito/RK/dual-time, multigrid, MPI multiplo/OpenMP multiplo,
3D/elementi non triangolari, moto/deformazione, restart in lettura,
DUAL_MEAN/SEGMENTED e le altre modalità non implementate; differenziazione,
modelli aggiuntivi e sorgenti native, incluso `KIND_VERIFICATION_SOLUTION`.
Il driver richiede `REF_DIMENSIONALIZATION=DIMENSIONAL`, con riferimenti
unitari; non modifica i default globali per gli altri solver.
`CONV_NUM_METHOD_FLOW=ROE` è il selettore già richiesto dal percorso
primitivo: non attiva il flusso Roe conservativo. Le opzioni esplicite
del flusso conservativo come ROE_KAPPA/ENTROPY_FIX_COEFF/LOW_MACH sono
rifiutate per evitare che sembrino operative.

### Configurazione, inizializzazione e bordo

Le quattro opzioni nuove sono:

| Opzione | Significato/default |
|---|---|
| `PRIMITIVE_INITIALIZATION` | UNIFORM oppure ISENTROPIC_VORTEX; default UNIFORM |
| `PRIMITIVE_BOUNDARY_MODE` | UNIFORM oppure ISENTROPIC_VORTEX; default UNIFORM |
| `PRIMITIVE_BOUNDARY_STATE` | Stato esterno rho,u,v,p; default (1,0.5,0.3,1) |
| `PRIMITIVE_GEOMETRY_AUDIT` | Nome CSV opzionale, vuoto per disabilitare |

Lo stato iniziale uniforme usa l'opzione esistente
`PRIMITIVE_INITIAL_STATE`; non è derivato da MACH_NUMBER o dalle opzioni
freestream conservative. Queste ultime restano metadati richiesti dal
parsing SU2, non definiscono W in questo solver. Il vortice sostituisce
lo stato iniziale ai nodi a t=0, con gamma=1.4, beta=2, centro (-.25,-.25),
freestream (1,1,1,1), dominio [-5,5]^2 e profilo esatto FASE 7.
La modalità analitica richiede inizializzazione e bordo entrambi vortice.

Solo `MARKER_CUSTOM` e `MARKER_FAR` sono ammessi. Entrambi qui selezionano
lo stato esterno primitivo prescritto, uniforme oppure analitico, per
la callback debole path-conservative. Non sono una parete né una BC
caratteristica generale di ingresso/uscita. La callback vortice usa
W*(q,t_n) nei punti reali delle mezze facce, senza sostituzione con il
freestream e senza sorgente. Non viene usata `CEulerSolver::BC_Custom`.
L'insieme dei marker mesh deve coincidere con quello configurato;
duplicati, mancanti, extra o tipi non implementati falliscono.

### Tempo fisico e arresto

Richiesti TIME_DOMAIN=YES, TIME_MARCHING=TIME_STEPPING,
TIME_DISCRE_FLOW=EULER_EXPLICIT, MGLEVEL=0.

* **Passo fisso:** UNST_CFL_NUMBER=0, TIME_STEP>0. MAX_TIME esplicito deve
  essere multiplo intero del passo, entro 128 epsilon sul rapporto.
  TIME_ITER assente viene derivato da MAX_TIME; con entrambi si usa il
  minore dei due limiti in passi. Senza MAX_TIME, si pone MAX_TIME=Nt*dt.
  Nessun accorciamento dell'ultimo passo fisso. Per il vortice si impongono
  Nt=2N² e dt=0.5/Nt (sequenza FE_d8 FASE 7).
* **CFL:** UNST_CFL_NUMBER>0, TIME_STEP=0 e TIME_ITER esplicito come cap.
  Si usa la stima esistente, poi dt=min(dt_CFL,MAX_TIME-t_n) se il limite
  temporale interviene. Questo accorciamento finale è esplicito e
  documentato; non è un clipping dello stato né un retry di positività.
  Senza MAX_TIME si arresta al cap dei passi.
* Un passo FE non ammissibile termina con la diagnostica originale;
  nessuna riduzione automatica del passo o modifica delle primitive.

Inizializzazione e output iniziale avvengono una volta. Ogni chiamata
Solve accetta un solo FE, aggiorna CConfig dal tempo reale del solver e
registra il contatore dei passi accettati. Indici output/history 0..Nt:
zero è lo stato iniziale, l'indice n è W^n. Monitor e output restano quelli
di CSinglezoneDriver; non c'è un secondo aggiornamento nell'Update di base.

## 8.3 Geometria e output nativo

Le mesh `.su2` passano da CPhysicalGeometry e dalla costruzione FVM nativa.
Il generatore degli esempi riprende coordinate Float64 e triangoli dagli
archivi FASE 7, con scrittura a 17 cifre. RCM cambia l'ordine locale:
i confronti vengono ricondotti a ID globali e orientazione canonica
dell'edge, invertendo la normale quando si scambiano gli estremi.
L'audit opzionale registra nodi, edge e mezze facce. Per le righe `node`,
le colonne nx/ny contengono **offset del centroide**, non normali;
per edge/boundary sono normali integrate. `second` nel boundary è il marker.

`CPrimitiveEulerOutput` deriva da CFVMOutput e usa sorter e writer SU2.
Gruppi di volume: COORDINATES, SOLUTION (Density, Velocity_x, Velocity_y,
Pressure), PRIMITIVE (Speed, Mach, Temperature, Entropy, Dual_Volume,
Physical_Time, Global_ID). Il writer VTK raggruppa Velocity come vettore.
Le convenzioni sono |v|, Mach=|v|/sqrt(gamma*p/rho),
T=p/(rho*GAS_CONSTANT), Entropy=p/rho^gamma. In particolare T non è theta
del vortice salvo GAS_CONSTANT=1. I riferimenti sono unitari, coerenti con
DIMENSIONAL: nessuna nondimensionalizzazione implicita delle primitive.

Per componente, la history PRIMITIVE_RHS è la norma lineare

    L2_rhs(k) = sqrt(sum_i V_i rhs_i(k)^2 / sum_i V_i),
    rhs_i = -LinSysRes_i/V_i = (W_i^{n+1}-W_i^n)/dt a meno del roundoff.

Include tutti i nodi, anche il bordo, senza logaritmo. La riga al tempo
t_{n+1} riporta il RHS di stadio valutato a t_n, mentre MIN_RHO/MIN_P sono
dello stato aggiornato e MIN_THETA dell'ultima ricostruzione preparata.
La riga iniziale ha RHS zero come placeholder dei work array: **non**
è una valutazione del residuale iniziale. TIME_STEP iniziale è zero.

OUTPUT_WRT_FREQ è una frequenza positiva in passi accettati. Iniziale e
finale sono forzati dal percorso nativo; intermedi agli indici richiesti.
PARAVIEW_MULTIBLOCK produce file `.vtm` e blocchi Internal/bordo, leggibili
come serie in ParaView. Non viene creato un nuovo writer o un PVD speciale.
Il writer XML SU2 usa **Float32**. OUTPUT_FILES=RESTART produce checkpoint
nativi con coordinate e SOLUTION in **Float64**, non i campi derivati;
il token è RESTART, non RESTART_BINARY. Tempo e ID del checkpoint vengono
associati tramite indice file/history e coordinate/CSV geometrico.
Questo output non implica supporto al caricamento restart del primitivo.

## 8.4 Risultati riproducibili A–F

Compilatore: Apple clang 15.0.0 (clang-1500.0.40.1), target
x86_64-apple-darwin22.6.0; build locale esistente, senza nuove dipendenze.
I CSV `primitive_euler_phase8_*` riportano dati completi, tempi wall e
risultati. `primitive_euler_phase8_provenance.json` registra HEAD, compiler
e SHA256 di sorgenti modificate, eseguibile, configurazioni, mesh,
checkpoint, geometrie, history e CSV della campagna.

### A — Configurazione e avvio

I tre esempi inclusi sono stati eseguiti dall'eseguibile ordinario e
terminano Exit Success. L'iterazione verifica con dynamic_cast che il
FLOW_SOL sia realmente CPrimitiveEulerSolver, altrimenti fallisce.
La factory mantiene la creazione originale delle variabili primitive;
nessun CEulerSolver sostitutivo. Il primo passo N=16 regular confrontato
con harness ricompilato differisce al massimo di 2.22044604925e-16.

### B — Stato uniforme

Regular e jitter N=8, 8 passi dt=.01, campioni t=0,.04,.08:
massimo errore di ciascuna delle quattro componenti **zero**; densità e
pressione positive, RHS WLS zero. Ulteriori test CLI GG/LS e limiter
confermano preservazione dello stato: massimo RHS GG=3.68215e-16,
compatibile con il limite di roundoff 64 epsilon max|W|/h=3.41061e-14.
Non è stata richiesta uguaglianza bit a bit di un gradiente GG ottenuto
da somme geometriche; le soglie dei test preesistenti restano invariate.

### C — Driver/harness e identità geometrica

Otto mesh regular/jitter N=8,16,32,64, tre tempi 0,.25,.5 e quattro
componenti: 96 confronti L1_V, L2_V, Linf. Le norme L1/L2 nei CSV sono
normalizzate per il volume totale, come sopra. **Massimo Linf=4.4408920985e-16**.
Soglia dichiarata prima della campagna: 2e-13 sullo stato O(1), stretta
rispetto agli errori di discretizzazione e coerente con i confronti
precedenti. L'accordo osservato è assai più stretto.

Il riferimento è l'archivio Float64 FASE 7, con un controllo aggiuntivo
del vero harness ricompilato contro le librerie correnti: N=16 delle
due famiglie, massimo scarto 4.44089e-16 e 2.22045e-16. Per ciascuna
famiglia, tutti i **1217 record geometrici**, inclusi 128 mezze facce,
hanno differenza massima **zero** dopo la canonicalizzazione. Coordinate
dei checkpoint identiche; volumi verificati anche su tutte le otto mesh
con rtol=2e-14, atol=2e-15. Connettività verificata nel roundtrip VTK.

### D — Driver/analitica

Errori globali normalizzati L2_V a T=.5, colonne rho,u,v,p:

| Famiglia | N | E_rho | E_u | E_v | E_p |
|---|---:|---:|---:|---:|---:|
| regular | 8 | 5.630762e-3 | 2.546016e-2 | 2.575763e-2 | 7.773696e-3 |
| regular | 16 | 2.230226e-3 | 9.121941e-3 | 8.370496e-3 | 3.066722e-3 |
| regular | 32 | 5.693887e-4 | 2.240051e-3 | 1.996820e-3 | 7.794709e-4 |
| regular | 64 | 1.221918e-4 | 4.846532e-4 | 4.379658e-4 | 1.667619e-4 |
| jitter | 8 | 6.140892e-3 | 2.878792e-2 | 2.443754e-2 | 8.525017e-3 |
| jitter | 16 | 2.253626e-3 | 8.536969e-3 | 6.834606e-3 | 3.095602e-3 |
| jitter | 32 | 6.563625e-4 | 2.218890e-3 | 2.141953e-3 | 9.008753e-4 |
| jitter | 64 | 1.627282e-4 | 4.800411e-4 | 4.828401e-4 | 2.233695e-4 |

| Famiglia | Raffinamento | p_rho | p_u | p_v | p_p |
|---|---|---:|---:|---:|---:|
| regular | 8→16 | 1.336140 | 1.480829 | 1.621615 | 1.341903 |
| regular | 16→32 | 1.969704 | 2.025809 | 2.067609 | 1.976130 |
| regular | 32→64 | 2.220266 | 2.208507 | 2.188814 | 2.224705 |
| jitter | 8→16 | 1.446200 | 1.753668 | 1.838169 | 1.461483 |
| jitter | 16→32 | 1.779683 | 1.943886 | 1.673931 | 1.780821 |
| jitter | 32→64 | 2.012029 | 2.208608 | 2.149309 | 2.011895 |

I risultati riproducono FE_d8 FASE 7 entro il roundoff del test C:
massimo scarto delle 32 norme L2 finali rispetto al CSV precedente
1.73472347598e-17.
Si riutilizza quindi il suo controllo temporale d2/d4/d8: non è stata
rilanciata una seconda campagna di quattro dt. Gli ordini sopra sono
osservati sulla soluzione evoluta di questo vortice; non provano ordine
due universale del residuale su ogni mesh irregolare o al bordo.
Tempi del run finale N=64: 57.75 s regular e 62.55 s jitter; il CSV
execution include tutti i tempi (dipendono dal carico della macchina).

### E — Output, history e roundtrip

ParaView 5.13.2/pvpython ha letto **30 snapshot** VTM/VTU, i blocchi,
connettività triangolare e campi. Stato esportato, coordinate, ID e tempo
sono esattamente i cast Float32 dei dati Float64 attesi; tutti i campi
derivati sono finiti, temperatura positiva. Quantizzazione massima dello
stato: 5.96044545e-8. Le norme ricalcolate soddisfano la disuguaglianza
triangolare con il limite di quantizzazione misurato per ogni campo;
massima differenza fra norme Float64 e VTU=2.85388953e-8. Non si applica
al writer Float32 la tolleranza di equivalenza fra solver Float64.

History con Nt+1 righe, indici senza duplicati, campioni e tempi finali
corretti. Controllo indipendente su otto passi N=16 usando differenze
dei checkpoint consecutivi: massimo scarto history/RHS=6.48093e-15,
inferiore al limite di sottrazione Float64 16 epsilon max|W|/dt (~4.76e-12).

### F — Modalità non supportate e precedenze

Tutti i 16 test negative.csv terminano con errore e sottostringa
diagnostica prevista: implicito, multigrid, dual-time, marker sconosciuto,
marker inutilizzato, wall, symmetry, modalità iniziale invalida,
DUAL_MEAN, tempo finale non multiplo, moto, conflitto CFL/passo fisso,
frequenza output nulla, MMS nativa, 3D e quadrilateri. Test aggiuntivo
`mpiexec -n 2` rifiutato esplicitamente. Un FE di dt=100 non ammissibile
viene rifiutato senza checkpoint accettato di indice 1.

Le verifiche temporali native confermano: senza MAX_TIME 8 passi→.08;
senza TIME_ITER ma MAX_TIME=.08 ancora 8; cap finale .04→4 passi;
CFL con cap100 e MAX_TIME=.08→3 passi e tempo .08; CFL con solo cap3
→3 passi e tempo .1020376691454. MARKER_FAR esterno uniforme passa.
EULER e NAVIER_STOKES percorrono il driver ordinario e terminano senza
errore negli smoke test uniformi di otto passi. Questo è un controllo
di compatibilità, non la riesecuzione dell'intera suite SU2 per tutti
i solver e tutte le fisiche.

### Regressioni e integrità dei risultati precedenti

* Suite primitiva ricompilata: **70 casi, 22 253 237 asserzioni**, tutti passati.
* Profilo vortice condiviso: **890 controlli**, inclusi derivati e assenza
  di sorgente; campioni prima/dopo estrazione identici bit a bit.
* MMS analitica esistente: **3300 asserzioni**, passate.
* Test Python output vortice: 4 passati; postprocessing asintotico: 3 passati.
* SHA256 di **98 documenti/artefatti preesistenti** identici al salvataggio
  iniziale, inclusi tutti gli snapshot FASE 7. Non sono stati rigenerati
  o sostituiti i risultati FASI 6–7.

## 8.5 Esempi, riproduzione e limiti

Gli esempi autonomi sono in `TestCases/primitive_euler_examples`:
uniform/uniform.cfg e .su2, vortex/vortex_regular.cfg e .su2,
vortex/vortex_jitter.cfg e .su2, README e .gitignore degli output.
Mesh N=16 dei vortici incluse: nessun harness preliminare necessario.
I nomi di output regular/jitter sono distinti. In ParaView aprire le
serie `vortex_regular_*.vtm` oppure `vortex_jitter_*.vtm`, Apply,
blocco Internal; Physical_Time è disponibile anche come campo.

Comandi effettivamente eseguiti, dalla radice salvo dove indicato:

```sh
./ninja -C build SU2_CFD/src/SU2_CFD -j4
# cwd TestCases/primitive_euler_examples/uniform
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD uniform.cfg
# cwd TestCases/primitive_euler_examples/vortex
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD vortex_regular.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD vortex_jitter.cfg
```

Campagna e diagnostiche, dalla radice; gli script impostano un thread:

```sh
python3 UnitTests/build_primitive_euler_driver_reference.py --output-dir /private/tmp/su2-phase8-reference
python3 UnitTests/run_primitive_euler_driver_e2e.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/su2-phase8-final
python3 UnitTests/run_primitive_euler_driver_e2e.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/su2-phase8-final --reference-only --reference-executable /private/tmp/su2-phase8-reference/reference_vortex
python3 UnitTests/check_primitive_euler_driver_modes.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/su2-phase8-modes-verified --reference-dir /private/tmp/su2-phase8-final
/Applications/ParaView-5.13.2.app/Contents/bin/pvpython UnitTests/validate_primitive_euler_driver_output.py --input-dir /private/tmp/su2-phase8-final --report /private/tmp/su2-phase8-final/vtk.csv
python3 UnitTests/export_primitive_euler_driver_audit.py --input-dir /private/tmp/su2-phase8-final --modes-dir /private/tmp/su2-phase8-modes-verified --output-dir /private/tmp/su2-phase8-package
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-phase8-tests
clang++ -std=c++11 -O2 UnitTests/SU2_CFD/euler_vortex_analytic_tests.cpp -o /private/tmp/su2-phase8-vortex-analytic
/private/tmp/su2-phase8-vortex-analytic
/private/tmp/su2-mms-6b1-final/euler_mms_analytic --reporter compact
python3 UnitTests/test_euler_vortex_output.py
python3 UnitTests/test_euler_asymptotic_postprocessing.py
git diff --check
```

Il generatore `primitive_euler_driver_cases.py` crea le otto risoluzioni
dai riferimenti esistenti; `run_primitive_euler_driver_e2e.py` esegue soltanto
il vero SU2_CFD per il percorso sotto test. Il reference builder aggiunge
export geometrici/primo passo al solo harness di verifica in /private/tmp.
L'exporter raccoglie 11 CSV e un manifest senza rilanciare le simulazioni.
I file di output nativo completi restano nelle directory della campagna
e degli esempi, non vengono aggiunti agli artefatti di sorgente del repository.

Restano fuori perimetro MPI/OpenMP parallelo, 3D, mesh ibride/mobili,
restart in lettura, BC fisiche ulteriori, sorgenti MMS via driver,
driver multizona, adattività e integratori diversi da FE. Non si introduce
alcuna correzione agli shock: velocità di salto e Rankine–Hugoniot non
sono validate da questi test lisci. La baseline resta il metodo
path-conservative originale, con i limiti scientifici già documentati.

### Sicurezza Git e ambiente

Nessun commit, push, fetch, cambio branch o modifica intenzionale della
cronologia/remoti. HEAD resta quello iniziale; tutte le modifiche sono
revisionabili nel working tree. Durante la prima rigenerazione della
build, il preconfigure Meson esistente ha tentato automaticamente una
scrittura Git di configurazione nel submodule Catch2 di CoolProp:
il sandbox l'ha bloccata con errore di lock. Non è stata concessa
escalation Git, né modificato quel meccanismo o i remote per aggirarlo;
la build ha proseguito e le successive invocazioni Ninja sono riuscite.
I programmi MPI, anche a un rank, e pvpython hanno richiesto autorizzazione
di esecuzione per i socket locali; non sono state aggiunte dipendenze.
