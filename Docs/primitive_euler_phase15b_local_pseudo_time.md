# FASE 15B — Iterazione stazionaria in pseudo-tempo locale

## Implementazione e configurazione

Implementato `CPrimitiveEulerSolver::AdvanceLocalPseudoTime`, separato da
`AdvanceExplicitEuler`, il cui corpo è rimasto invariato. Tutto il lavoro della
FASE 15A è conservato nel working tree, inclusi configurazioni e risultati.
Nessun commit o altra operazione Git di scrittura.

Nuova opzione, default NO, disponibile solo per PRIMITIVE_EULER:

```ini
PRIMITIVE_LOCAL_PSEUDO_TIME=YES
CFL_NUMBER=.05
TIME_ITER=300
TIME_STEP=0
UNST_CFL_NUMBER=0
TIME_DOMAIN=YES
TIME_MARCHING=TIME_STEPPING
TIME_DISCRE_FLOW=EULER_EXPLICIT
MGLEVEL=0
```

Non impostare MAX_TIME. TIME_DOMAIN/TIME_STEPPING sono soltanto il contenitore
del driver esistente; TIME_ITER conta iterazioni accettate. Il tempo fisico del
solver e di CConfig resta zero, non viene avanzato. Inizializzazione/BC vortex
dipendenti dal tempo vengono rifiutate dal driver. Le callback C++ dell'API
devono rappresentare dati autonomi: vengono campionate a t=0 e non devono
modificare esternamente lo stato del solver.

Configurazione completa: [primitive_local.cfg](../TestCases/primitive_euler_examples/naca0012_aoa0/primitive_local.cfg).
La validazione temporale comune ammette TIME_STEP=0/UNST_CFL=0 **soltanto** per
questa modalità primitiva; comportamento nativo e default fisico invariati.

## Formula CFL e aggiornamento

Si riutilizza `SpatialResult::spectralRadius`, già calcolato dall'operatore
primitivo per la stima CFL globale. Per una faccia con normale **integrata** S,
`PrimitiveFaceSpeed` prende il massimo fra velocità spettrale del cammino e
stime fisiche dei due stati:

    a_f = max(path_speed, |u_L·S| + c_L |S|, |u_R·S| + c_R |S|)
    sigma_i = sum_incident_faces a_f
    tau_i = CFL * V_i / sigma_i
    W_candidate_i = W_old_i + 2^(-j) tau_i RHS_i(W_old)

Nessun ulteriore fattore geometrico o empirico: |S| è già incluso in a_f.
Ogni faccia interna viene valutata una volta e contribuisce una volta a
ciascun estremo. Ogni mezza faccia fisica contribuisce soltanto al proprio
nodo, usando interno ricostruito ed esterno della BC esistente. Normali
orientate, stati, volumi, somme spettrali e passi devono essere validi,
finiti e positivi dove richiesto. Nessuna somma delle normali medie ai vertex
usata come sostituto delle facce.

Il RHS iniziale completo viene congelato **prima di ogni scrittura nodale**.
Tutti i candidati e i retry usano esclusivamente quello stesso W_old/RHS_old.
Solo dopo avere formato il candidato completo se ne ricalcola il residuo:
questo serve a verificarne slip/ricostruzione/finitudine e a fornire la
diagnostica corrente; non alimenta un aggiornamento asincrono.

Il RHS di velocità sui nodi forti è già tangente: moltiplicarlo per uno scalare
tau_i preserva lo slip entro roundoff. Al TE autorizzato non si introduce
alcuna proiezione; riflessioni, due D-minus e condizioni geometriche restano
quelle della FASE 13C/15A.

## Accettazione e diagnostiche

Tentativi j=0..8 di default: CFL dimezzato globalmente dopo ciascun candidato
non ammissibile, preservando i rapporti dei passi locali. Nessun clipping.
Un errore imprevisto di callback viene rilanciato; non è nascosto come successo.

Se tutti i tentativi falliscono, stato corrente, Solution_Old, LinSysRes,
tempo, gradienti raw/limitati, limiter, Rmatrix e fattori di positività tornano
ai valori d'ingresso. Restano aggiornati soltanto i contatori di rifiuto.
All'accettazione Solution_Old contiene W_old, le cache corrispondono al nuovo
campo e LinSysRes contiene il **residuo corrente ricalcolato**. Questo differisce
volutamente dalla semantica stage del metodo fisico, preservata integralmente.

`GetPseudoDiagnostics()` espone vettore dei passi accettati, minimo/massimo,
rapporto, CFL effettivo, riduzioni dell'ultimo successo e contatori cumulativi.
History: `LOCAL_PSEUDO_TIME` e `PRIMITIVE_STEADY_AUDIT`; RHS globale/interno/TE,
CL/CD, minimi rho/p, Mach nodale e delle tracce. `PSEUDO_REJECTED` conta chiamate
fallite, `PSEUDO_REJECTED_ATTEMPTS` include tentativi scartati prima di successi.
Prima del primo passo i diagnostici dei passi valgono zero, non una stima.

Nei VTU locali `Pseudo_Iteration` e `Local_Pseudo_Step` sostituiscono il campo
Physical_Time; CSV TE etichettati `pseudo_iteration`. Le sequenze PVD hanno
ascisse di iterazione. Il framework può ancora riportare Time_Iter o TimeValue
nei metadati: non sono tempi fisici. Le history del caso locale escludono
Cur_Time/Time_Step. Non si usa mai la norma dell'aggiornamento per certificare
convergenza.

## Test e regressioni

**4 nuovi casi, 388 asserzioni**, passati:

- quiete uniforme invariata con passi non uniformi;
- aggiornamento sincrono e stima locale verificati indipendentemente;
- candidato a passi imposti uguali identico al Forward Euler fisico;
- con CFL=1e-300 l'aggiornamento si arrotonda a zero, ma il residuo resta non
  nullo: nessuna falsa diagnostica di convergenza;
- slip dei nodi forti e TE debole non proiettato su una geometria di prova;
- sink su rho/p: tre dimezzamenti, CFL .1→.0125, rapporti locali conservati;
  esaurimento dei tentativi ed eccezione di callback con rollback delle cache;
- rifiuto di CFL/passi/stati/volumi invalidi e overflow spettrale.

**35 regressioni, 334133 asserzioni**, passate: FE, stima CFL globale, hook
espliciti, costanti/affini, geometria, pareti, TE nativo e debole, aerodinamica.
CLI poligonale passata. Smoke fisico NACA di cinque passi: restart Float64
identici a quelli archiviati in 15A, **max differenza=0**, inclusi iniziale e
snapshot intermedi/finali. Nessuna campagna MMS/vortice completa ripetuta.

## Esperimento NACA breve

Mesh originale immutata, Mach .2, AOA=0°, gamma=1.4, freestream e normalizzazione
di 15A, WLS/MUSCL affine, limiter NONE, LOCAL_CHARACTERISTIC, TE globale 199.
Smoke indipendente di 50 iterazioni seguito da 300 iterazioni da freestream;
budget 120 s per run. Entrambi riusciti: nessun tentativo rifiutato, CFL effettivo
sempre .05. Nessuna campagna lunga o caso AOA=2°.

| Modalità | Iterazioni | RHS globale finale | RHS interno finale | RHS TE finale | CL | CD |
|---|---:|---:|---:|---:|---:|---:|
| Iniziale comune | 0 | .0470270 | .0260361 | 140.245 | 0 | 0 |
| Locale smoke | 50 | .0421143 | .0258864 | 96.0223 | 5.32e-9 | .309271 |
| Globale FE, 15A | 300 | .0456527 | — | 20.7729 | -1.61e-8 | .182687 |
| Locale | 300 | .0275565 | .0283484 | 10.9326 | 3.05e-6 | .683667 |
| Globale FE, 15A | 3000 | .0339300 | .0341839 | 1.32344 | -1.22e-7 | .697670 |

Riduzione globale locale: rapporto finale/iniziale=.58597, solo .232 decadi;
l'interno **aumenta dell'8,9%**. Il TE resta lontano da residuo nullo. Il
miglioramento globale è misurabile, ma non è uniforme per regioni e non
certifica stazionarietà o accuratezza aerodinamica.

Al passo 300: tau_min=tau_TE=2.63891753e-6, tau_max=.0241533021,
rapporto=9152.73. Percentili 10/50/90 del campo salvato:
4.03207e-5 / 2.78073e-4 / .00805206 (Float32). Il TE resta il collo di bottiglia
del passo minimo, non viene artificiosamente accelerato o proiettato.

Run locale: 14.01 s, .0467 s/iterazione inclusi startup, diagnostica e I/O;
FE 15A: 97.18 s/3000, .0324 s/iterazione. Tempi osservati con compilazioni
concorrenti, non benchmark puro. Il locale costa di più per iterazione, ma in
questo campione raggiunge un RHS globale inferiore a quello di 3000 passi FE.

Minimi locali su tutte le iterazioni: rho=.9689123373, p=.9567426016;
max Mach nodale=.2080631643, max tracce=.2708199245. TE finale:
rho=.980873324, u=.212673366, v=-1.70e-9, p=.973346430; tracce normali circa
-.03015698 su entrambe le mezze facce. D-minus e RHS completi nei CSV.
Campi finiti, slip e TE verificati dai test; nessuna stabilità a lungo termine
dedotta da questo esperimento.

## Riproduzione, file e seguito 15C

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase15b-tests
/private/tmp/primitive-phase15b-tests/test_driver '[primitive-pseudo-time]'
/private/tmp/primitive-phase15b-tests/test_driver '[primitive-native-te-audit],[primitive-weak-te],[primitive-trailing-edge],[primitive-aerodynamics],[primitive-geometry],[primitive-boundary],[primitive-wall-gradients],Forward Euler*,Rejected time steps*,CFL estimate*'
/private/tmp/primitive-phase15b-tests/test_driver 'Native explicit hooks*,Constant primitive fields*,Affine NODE*'
python3 UnitTests/run_primitive_euler_local_pseudo15b.py --output-dir TestCases/primitive_euler_examples/naca0012_aoa0/results15b --iterations 300 --cfl .05 --seconds-per-run 120
python3 UnitTests/check_primitive_euler_polygon_cli.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/primitive-phase15b-polygon
git diff --check
```

Lo script richiede una destinazione nuova; `--postprocess-only` rigenera solo
controlli/tabelle. Legge la baseline 15A senza sovrascriverla; se lo smoke mostra
crescita netta del residuo riprova a CFL dimezzato e non estende un caso ancora
in crescita. Il run svolto non ha richiesto tale riduzione.

[Risultati e distribuzioni](../TestCases/primitive_euler_examples/naca0012_aoa0/results15b/summary.csv),
[sequenza locale](../TestCases/primitive_euler_examples/naca0012_aoa0/results15b/diagnostic/iterations.pvd).
Verificati gli array VTU e la loro etichettatura per iterazione, senza nuovo
writer dei campi.

File 15B: CConfig.hpp/.cpp; CPrimitiveEulerSolver.hpp/.cpp;
CPrimitiveEulerIteration.cpp; CPrimitiveEulerOutput.cpp; test time/boundary;
nuovo script locale, primitive_local.cfg, questo rapporto e risultati15b.
Factory/adapter TE e artefatti 15A restano presenti. Nessuna modifica a kernel,
BC, geometria, MUSCL o solver conservativo.

Per 15C: estendere inizialmente a CFL .05 con budget espliciti e checkpoint;
valutare separatamente globale/interno/TE e stabilizzazione CL/CD. Solo dopo
almeno 3–4 decadi di residuo e verifiche di stabilizzazione sarà sensato
confrontare due stati stazionari. Eventuali studi CFL devono mantenere lo stesso
operatore e registrare i rifiuti; un plateau non autorizza correzioni al kernel.
FASE 15C non eseguita e accuratezza aerodinamica ancora aperta.
