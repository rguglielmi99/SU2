# PRIMITIVE_EULER — FASE 5: Forward Euler

## Ambito consegnato

La timeline originale, sezione 16 del PDF, definisce la FASE 5 come
“Tempo: Forward Euler, poi SSP-RK2 con ricalcolo completo e positività”.
Questa consegna implementa **Forward Euler**, nell'ambito semplificato
precedentemente concordato. La ricerca nel codice locale non identifica
un integratore SSP-RK2 direttamente riutilizzabile dal solver primitivo.
SSP-RK2, altri RK e implicito restano rifiutati; non si presenta il RK
multistadio generico SU2 come equivalente SSP-RK2.

Il metodo resta NODE, MUSCL affine, AGGREGATED e primitive W=(rho,u,v,p).
I due volumi preesistenti sono mantenuti; NODE_POINT rimane il riferimento
consigliato. Non si reintroducono fit/Hessiane produttive o SEGMENTED.
Kernel, matrici, segni delle fluttuazioni, limiter e geometria sono invariati.
Gli audit della FASE 4 e i loro limiti rimangono validi.

Branch iniziale `feature/primitive-path-solver`, HEAD `7dcbbdbc18`, working tree
pulito: gli audit precedenti erano già stati committati dall'utente.
Il branch sperimentale non viene modificato. Nessun commit, push o sync
è effettuato da questo lavoro.

## Nuove API e significato del passo

In CPrimitiveEulerSolver:

```cpp
SpatialResult AdvanceExplicitEuler(CGeometry*, const CConfig*,
    double time, double dt, const BoundaryState& boundary = {}, const Source& source = {});
double EstimateTimeStep(CGeometry*, const CConfig*,
    double time, double cfl, const BoundaryState& boundary = {});
double GetPhysicalTime() const;
```

L'avanzamento ricostruisce da zero L(W,t) usando lo stato numerico corrente:
gradienti, limiter, theta, tracce, fluttuazioni, volume, bordo e sorgente.
Le callback vengono valutate al tempo **time**, non time+dt. Poi applica:

```text
W_next_i = W_i + dt * rhs_i
rhs_i = source_i - (faces_i + volume_i + boundary_i)/V_i
```

Non divide nuovamente per il volume e non usa un LinSysRes precedente,
anche se il chiamante lo ha sovrascritto. Il risultato restituito è lo
SpatialResult dello stadio iniziale W,t, non un residuale al tempo finale.
Dopo successo, LinSysRes=-V*rhs dello stesso stadio, Solution_Old=W iniziale,
W=W_next e GetPhysicalTime()=time+dt.

Il parametro time esplicito è autorevole: consente un'inizializzazione a
un tempo scelto dal chiamante. Non implementa il restart SU2, che rimane
rifiutato dalla configurazione. Il clock iniziale del solver è 0.
Passi nulli, negativi, non finiti o non rappresentabili (time+dt<=time,
oppure somma non finita) vengono rifiutati.

Le coordinate, le unità degli stati e del tempo devono essere coerenti:
il solver usa direttamente i valori primitivi e le coordinate della mesh,
come nella FASE 4. Non converte automaticamente un TIME_STEP dimensionale
nelle unità di un'integrazione CFD generale.

## Transazione e positività

ComputeSpatialResidual mantiene il suo contratto: aggiorna il residuale
soltanto dopo una valutazione completa e finita. L'implementazione privata
EvaluateSpatialResidual può inoltre valutare senza commit, per stima CFL
e avanzamento. Una sola guardia RAII impedisce rientro/concorrenza fra
residuale, CFL e avanzamento, inclusi richiami dentro callback.
Le tre API devono essere chiamate fuori da regioni OpenMP.

AdvanceExplicitEuler costruisce tutti i candidati in un vettore temporaneo.
Ogni componente deve essere finita e rho,p devono essere strettamente
positive, usando la stessa validazione termodinamica del kernel.
Un candidato non ammissibile lancia std::domain_error con l'indice del nodo;
non viene modificato **alcun** nodo, neppure quelli validati prima del nodo
che fallisce. Solution_Old, LinSysRes e clock restano quelli precedenti.
Errori di tempo/configurazione/stato/sorgente seguono le validazioni
preesistenti; eccezioni delle callback vengono propagate.

Non si applicano clipping, riduzione nascosta del passo, theta sulle
sorgenti o modifiche delle fluttuazioni. Il chiamante può ridurre dt e
riprovare lo stesso problema; un test verifica rifiuto dopo un passo già
accettato e ripartenza con la stessa sorgente e dt ridotto.
Non si promette che esista un dt accettabile per qualunque sorgente o dato.

La positività di ricostruzione viene applicata prima della valutazione di
ogni stadio; quella nodale viene controllata dopo il candidato temporale.
Al passo seguente la ricostruzione viene ricalcolata sul nuovo W. Come nella
FASE 4, ReconstructAt/GetLimitedGradient descrivono l'ultima preparazione:
dopo un aggiornamento di W la diagnostica non è automaticamente un gradiente
fresco; richiamare ComputeSpatialResidual prima di interpretarla come tale.
Gli array ausiliari possono essere stati ricalcolati anche in un passo
rifiutato. La transazione riguarda W, Solution_Old, residuale e clock;
non promette rollback di ogni array diagnostico. Le callback non devono
mutare direttamente il solver/configurazione durante la valutazione.

## Stima CFL globale

SpatialResult espone anche spectralRadius, un vettore nodale di somme di
velocità caratteristiche integrate (non divise per V). Per una faccia:

```text
sigma_path = abs(ubar dot S) + sqrt(gamma*pbar*beta)*|S|
sigma_L/R  = abs(u_L/R dot S) + sqrt(gamma*p_L/R/rho_L/R)*|S|
sigma_face = max(sigma_path, sigma_L, sigma_R)
Lambda_i = sum_incident_internal sigma_face + sum_physical sigma_boundary
dt = cfl * min_i (V_i / Lambda_i)
```

Beta è l'integrale di 1/rho del percorso originale. La velocità di percorso
è inclusa esplicitamente: non si assume che sia sempre limitata dalle due
velocità del suono fisiche agli estremi. Le normali sono già integrate;
non si moltiplica nuovamente per lunghezza. Una faccia interna contribuisce
alla somma di entrambi i proprietari; le mezze facce fisiche al solo nodo
proprietario. Si usano tracce ricostruite e stato esterno al tempo indicato,
compreso l'effetto delle callback sulle velocità al bordo.

Si richiedono cfl e somme finite e positive. Il passo è **globale**, non
locale per nodo: resta un avanzamento in un unico tempo fisico. In 2D,
sotto omotetia x→h*x, V~h² e Lambda~h, quindi dt~h; i test verificano questa
scala e la linearità in CFL. EstimateTimeStep lascia invariati W, residuale,
Solution_Old e clock; prepara le diagnostiche e valuta il bordo ma non una
sorgente MMS. Il suo costo attuale include la valutazione spaziale completa,
un riuso semplice del percorso validato anziché un secondo assemblatore.

Questa è una **stima caratteristica prudenziale**, non una dimostrazione di
stabilità/positività per l'operatore con volume, limiter e sorgenti. Non limita
la rigidità delle sorgenti, né impedisce oscillazioni finite ma non fisiche.
AdvanceExplicitEuler non tronca un dt scelto esplicitamente alla stima CFL:
il chiamante seleziona il passo e il solver ne verifica l'ammissibilità nodale.
Partire da CFL=0.1; stabilità e dipendenza dal CFL restano da studiare.

## Hook virtuali SU2 e stato del driver generale

SetTime_Step valida la mesh finest, usa il CFL della configurazione corrente
al clock del solver e memorizza **lo stesso dt in tutti i nodi**, aggiornando
Min_Delta_Time e Max_Delta_Time. ExplicitEuler_Iteration richiede un passo
positivo già preparato e uniforme, poi chiama AdvanceExplicitEuler con
ricalcolo completo; non consuma un residuale esterno stantio.
Un test esegue venti passi di stato costante attraverso CSolver*.

Gli hook senza callback possono usare MARKER_FAR con lo stato esterno
preesistente PRIMITIVE_INITIAL_STATE. Per CUSTOM e sorgenti manufactured si
usano le API esplicite con callback. Non sono aggiunte nuove condizioni
fisiche o trasporto accoppiato.

**SU2_CFD config.cfg non è ancora un percorso primitivo evolutivo completo.**
Il driver generale resta bloccato esplicitamente: iteration, integration
factory e output non vengono collegati automaticamente ai percorsi Euler
conservativi, che hanno assunzioni diverse sul layout di Solution.
Il messaggio CLI ora distingue la disponibilità dell'API temporale da questo
limite. L'integration factory del primitivo rimane NONE; gli hook qui
verificati sono chiamabili esplicitamente, non un'integrazione del driver.
Questa consegna implementa il checkpoint temporale della timeline, non il
collegamento generale driver/output. Non dichiara eseguibile una simulazione
primitiva CLI né pronta la validazione scientifica della FASE 6.

## Esempio di uso C++

Con solver, geometry e config già inizializzati sulla mesh statica seriale:

```cpp
const CPrimitiveEulerSolver::BoundaryState boundary =
    [](const PrimitiveEuler2D::Direction&, double,
       const PrimitiveEuler2D::State& inside,
       const PrimitiveEuler2D::Direction&, unsigned short) { return inside; };
double t = 0;
const double end = 0.2;
while (t < end) {
  double dt = std::min(solver.EstimateTimeStep(geometry, config, t, 0.1, boundary), end-t);
  solver.AdvanceExplicitEuler(geometry, config, t, dt, boundary);
  t = solver.GetPhysicalTime();
}
```

La callback qui sopra è estrapolazione di prova, non un sostituto universale
di una condizione fisica. Per MMS usare i valori esterni esatti e passare la
sorgente come ultimo argomento. Per gestire std::domain_error si può dimezzare
dt e riprovare lo stesso t, imponendo un limite al numero di tentativi e
verificando che il passo resti rappresentabile; il solver non decide
silenziosamente una politica adattiva. Le eccezioni diverse dalla
non ammissibilità non vanno nascoste da una riduzione automatica del passo.

## Verifiche temporali

Sei casi nuovi:

1. Stati non uniformi, tempi e coordinate delle callback, confronto esatto
   con W+dt*rhs, segno e singolo fattore V, vecchia soluzione e residuale
   precedente artificialmente sovrascritto; ricalcolo al secondo passo.
2. Somme spettrali su geometria nativa, min globale, scala h/CFL, bordo con
   stato esterno veloce e stima che non modifica stato/residuale.
3. Rifiuti per dt/CFL invalidi, configurazione/geometria, candidati rho/p
   negativi, eccezioni/NaN da callback, guardia di rientro, storia già
   accettata e ripartenza con passo ridotto.
4. Hook virtuali, venti passi costanti, rifiuto di passi locali/non preparati.
5. ODE uniforme manufactured W(t)=W0 exp(t), sorgente analitica W0 exp(t),
   bordo estrapolato: spazio costante, T=0.2, 20/40/80/160 passi. Si misurano
   entrambi i gradienti e tutti e tre i limiter senza errore spaziale.
6. MMS affine accoppiata W(x,t)=exp(t)*(W_base+G*x), con G non banale,
   dati fisici esatti e sorgente W_t+B(W,G). Sullo stato analitico la
   ricostruzione è esatta e il residuale spaziale riproduce il differenziale:
   si isola il tempo pur ricalcolando tutti i contributi su stati evoluti.
   LS/WLS, limiter NONE, T=0.05, 32/64/128/256 passi.

Le norme finali sono L2 sui volumi duali, normalizzate con il volume totale.
I test richiedono ordine temporale fra 0.99 e 1.01 nell'ODE uniforme e fra
0.98 e 1.02 nella MMS accoppiata, dopo il primo livello. Non si pretende
ordine temporale 2 da Forward Euler. Per misurare futuro ordine spaziale 2
con questo integratore, dt=O(h) non basta: serve separare/ridurre l'errore
O(dt), come precisato nella sezione 10 del PDF.

### Dati misurati

Forzante uniforme: LS/WLS e tutti i limiter danno gli stessi valori entro
arrotondamento; qui è riportata WLS+VENKATAKRISHNAN.

| Passi | dt | L2 | Ordine |
| ---: | ---: | ---: | ---: |
| 20 | 0.01 | 1.690583139e-03 | 0.000000 |
| 40 | 0.005 | 8.459971528e-04 | 0.998796 |
| 80 | 0.0025 | 4.231749727e-04 | 0.999399 |
| 160 | 0.00125 | 2.116315854e-04 | 0.999699 |

MMS accoppiata con dati analitici spaziali e temporali:

| G | Passi | dt | L2 | Ordine |
| --- | ---: | ---: | ---: | ---: |
| LEAST_SQUARES | 32 | 0.0015625 | 1.528114613e-04 | 0.000000 |
| LEAST_SQUARES | 64 | 0.00078125 | 7.634049069e-05 | 1.001232 |
| LEAST_SQUARES | 128 | 0.000390625 | 3.815395011e-05 | 1.000616 |
| LEAST_SQUARES | 256 | 0.0001953125 | 1.907290309e-05 | 1.000308 |
| WEIGHTED_LEAST_SQUARES | 32 | 0.0015625 | 1.527971535e-04 | 0.000000 |
| WEIGHTED_LEAST_SQUARES | 64 | 0.00078125 | 7.633299967e-05 | 1.001239 |
| WEIGHTED_LEAST_SQUARES | 128 | 0.000390625 | 3.815011924e-05 | 1.000619 |
| WEIGHTED_LEAST_SQUARES | 256 | 0.0001953125 | 1.907096617e-05 | 1.000310 |

[CSV uniforme](primitive_euler_phase5_temporal.csv), 24 righe dati;
[CSV accoppiata](primitive_euler_phase5_coupled_temporal.csv), 8 righe.
Nessuna verifica di shock o di smooth equivalence conservativa viene
dedotta da questi dati temporali: appartengono alla successiva validazione.

## Build, comandi e stato di consegna

Dalla radice del repository:

```sh
./ninja -C build -j 2 SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-phase5-tests
/private/tmp/su2-phase5-tests/test_driver '[primitive-time]'
/private/tmp/su2-phase5-tests/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-phase5-cli --cli
python3 UnitTests/run_primitive_euler_tests.py --asan --output-dir /private/tmp/su2-phase5-asan
ASAN_OPTIONS=detect_leaks=0 /private/tmp/su2-phase5-asan/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
python3 UnitTests/export_primitive_euler_time.py --log /private/tmp/su2-phase5-results.log
```

I comandi dei test sono stati eseguiti con redirezione dei log sotto
/private/tmp. Dopo gli ultimi controlli di rifiuto/ripartenza è stato
ricompilato il solo primitive_euler_time_tests.cpp e ricollegato ciascun
driver riusando gli argv registrati dal runner (shlex.split/subprocess).
Un'esecuzione ex novo dei comandi sopra include direttamente tutti i
sorgenti correnti. Le dipendenze native erano aggiornate prima degli
ultimi link; la modifica CConfig riguarda soltanto il messaggio informativo.

La build nativa e quella finale incrementale sono terminate con uscita 0.
I quattro smoke CLI passano con le uscite attese 1/0/1/1: guardia primitiva
ancora esplicita, una iterazione del solver Euler conservativo riuscita,
harmonic balance e specie non supportati rifiutati. Non è installato alcun
eseguibile in SU2_bin.

La suite conserva i 64 casi degli audit precedenti e aggiunge sei casi
temporali. I due vecchi test sulle route “evoluzione indisponibile” sono
aggiornati al contratto: Forward Euler e timestep sono ora verificati in
funzione, mentre RK/implicito/centered/source separata rimangono rifiutati.
Nessuna soglia degli audit di accuratezza viene indebolita o riscritta.
Gli artefatti spaziali dei cicli precedenti rimangono invariati.

ASan è parziale: strumenta i test e le classi/factory primitive e basi delle
variabili ricompilate dal runner, non tutti gli archivi SU2/MPI. Le annotazioni
libc++ sono uniformate con _LIBCPP_HAS_NO_ASAN come nei cicli precedenti.
LeakSanitizer non è supportato sul macOS usato ed è disabilitato: non si
certifica assenza di leak. Warning non bloccanti noti: pragma GCC ignorati
su Clang, archivi senza simboli, librerie duplicate, WINDOW_START_ITER dei
test precedenti e bind MPI nel sandbox. Non è eseguita la suite globale SU2,
una simulazione MPI o la validazione scientifica della FASE 6.

Log temporanei principali:
`/private/tmp/su2-phase5-{build,build-final,tests-build,tests-rebuild,time,results,cli,asan-build,asan-rebuild,asan-results}.log`.
I dati CLI dettagliati sono in `/private/tmp/su2-phase5-cli/`.

File interessati:

- CPrimitiveEulerSolver.hpp/.cpp: Forward Euler, stima globale e transazione;
- CDriver.cpp, CIterationFactory.cpp e messaggio in CConfig.cpp: stato corretto
  dell'API temporale e guardia del driver, senza nuove opzioni;
- config_template.cfg: commenti aggiornati, default non modificati;
- primitive_euler_time_tests.cpp e registrazioni UnitTests/Meson/runner;
- adattamento dei due test sulle vecchie route bloccate;
- export_primitive_euler_time.py, questo rapporto, due CSV temporali e
  collegamento nel rapporto della FASE 4.

**Risultato nel perimetro concordato:** il solver primitivo ora avanza nel
tempo tramite Forward Euler, controlla l'ammissibilità dopo ogni passo e
misura ordine temporale 1 su riferimenti analitici. Restano da realizzare
SSP-RK2 e il collegamento generale driver/output, oltre alle verifiche della
soluzione/shock della FASE 6. Nessun limite spaziale già documentato è
presentato come risolto da questa integrazione temporale.

**Esito finale verificato:** 70 casi e **22.253.075 asserzioni superate**, con
uscita 0 sia nel driver normale sia nel driver ASan parziale aggiornato agli
ultimi controlli di rifiuto/ripartenza. Nessuna segnalazione ASan. I quattro
smoke CLI e la build nativa finale terminano come atteso. git diff --check
e la compilazione Python dei runner/esportatore passano; i due CSV hanno
24 e 8 righe dati con errori/ordini finiti. Kernel della FASE 1, geometria,
ricostruzione affine/limiter e variabili primitive non hanno diff rispetto
all'HEAD iniziale. Nessun commit o push.

## Rafforzamento finale dei due test Forward Euler

Il test di rifiuto ora accetta prima una sorgente uniforme non nulla
`(0.2,0.1,-0.05,0.3)` con dt=0.01. Lo snapshot salva tutti i nodi e le
componenti di W, Solution_Old e LinSysRes, oltre al clock: il massimo
residuale assoluto salvato è 0.2 e max|W-Solution_Old| è 0.003.
Dopo il passo non ammissibile si richiede identità esatta di tutto lo
snapshot. La ripartenza a dt=0.001 usa la stessa sorgente del passo rifiutato:
si verificano rhs contro la sorgente nota indipendentemente, W+dt*rhs,
Solution_Old, LinSysRes=-V*rhs e clock. I controlli precedenti su input,
callback e rientro rimangono presenti.

Prima della convergenza temporale, la MMS affine viene valutata sullo stato
analitico a t=0 e t=0.037, con gli stessi dati di bordo e sorgente, confrontando
ogni componente di rhs con W_exact, cioè la derivata temporale analitica.
La tolleranza assoluta è `64*epsilon_double*max(1,max|W_exact|)`:
64 unità di arrotondamento scalate sugli stati O(1) della patch fissa a cinque
nodi, con stencil ben condizionati. Non è una tolleranza di troncamento né
la tolleranza relativa predefinita di Approx.

| Gradiente | max difetto t=0 | max difetto t=0.037 | Tolleranza t=0 / t=0.037 |
| --- | ---: | ---: | ---: |
| LS | 8.881784e-16 | 1.776357e-15 | 5.272227e-14 / 5.470953e-14 |
| WLS | 8.881784e-16 | 1.776357e-15 | 5.272227e-14 / 5.470953e-14 |

I massimi sono compatibili con roundoff; la L2 pesata massima è 1.367015e-15
per LS e 1.247129e-15 per WLS. Il test isola quindi l'errore temporale nella
configurazione verificata. Gli stati iniziali sono ripristinati prima del
ciclo temporale originale; passi, dati analitici e soglie di convergenza
restano invariati. ASan dà massimi <=1.776357e-15, con piccole differenze di
arrotondamento rispetto alla build ottimizzata. Nessun difetto del solver
è emerso e nessun codice produttivo è modificato da questo rafforzamento.

Comandi eseguiti dopo ricompilazione del solo translation unit temporale e
relink dei due driver, riusando gli argv registrati dai runner:

```sh
/private/tmp/su2-phase5-tests/test_driver '[primitive-time]'
/private/tmp/su2-phase5-tests/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
ASAN_OPTIONS=detect_leaks=0 /private/tmp/su2-phase5-asan/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
git diff --check
```

Esiti: 6 casi / 1453 asserzioni per la selezione temporale; **70 casi /
22.253.237 asserzioni** per la suite completa, sia normale sia ASan parziale,
tutti superati con uscita 0. Valgono i limiti ASan già dichiarati sopra.
Log `/private/tmp/su2-phase5-hardening-{tests-build,asan-build,time,results,asan-results}.log`.
Tutti i record diagnostici spaziali della FASE 4 e quelli delle convergenze
temporali coincidono riga per riga con il run precedente. CSV non rigenerati.
Solo questo rapporto e primitive_euler_time_tests.cpp sono aggiornati in
questo intervento; le modifiche preesistenti del working tree sono preservate.
