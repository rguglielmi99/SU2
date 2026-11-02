# Euler 2D primitivo — FASE 4: operatore spaziale

Il [rapporto finale sui limiter](primitive_euler_phase4_limiter_audit.md)
misura direttamente l'operatore affine con NONE, VENKATAKRISHNAN e
BARTH_JESPERSEN, senza modificare produzione o iniziare la FASE 5.

Il nuovo [Passaggio I sull'irregolarità e sul contatto](primitive_euler_phase4_irregular_audit.md)
separa il difetto dissipativo senza MUSCL dalla consistenza debole scalare,
e confronta gradienti/segmenti su jitter e deformazioni lisce. Non modifica
il solver e non dichiara risolti questi limiti in produzione.

Il branch principale mantiene il solver **MUSCL affine NODE / AGGREGATED**,
le correzioni di robustezza e le diagnostiche indipendenti dell'ordine.
L'esperimento quadratico completo è archiviato nel branch locale
`experiment/primitive-quadratic-boundary`, commit `d1732e90e3`.
Fit, Hessiane di produzione, cache e opzione quadratica non appartengono
al solver di questo branch. Vedere il [rapporto di separazione](primitive_euler_branch_archive.md).

L'[audit dell'operatore affine](primitive_euler_phase4_order_audit.md) documenta
ordine globale circa 1.5 sulla famiglia regolare e circa 1 sulle mesh distorte.
Il criterio globale >=2 non è soddisfatto da questo operatore. Restano aperti
il difetto di consistenza MUSCL NO e i limiti di MEAN_COEFF.
La FASE 5 non è implementata; il percorso evolutivo CLI rimane bloccato.

La fase 4 implementa il residuale spaziale per **NODE / AGGREGATED**, su mesh
triangolari 2D statiche e seriali. Le incognite restano esclusivamente
`W=(rho,u,v,p)`. Non cambia il kernel path/Roe della fase 1 né il solver Euler
conservativo. La fase 5 introdurrà l'avanzamento temporale e il collegamento
completo al driver/output: una simulazione primitiva CLI resta esplicitamente
bloccata, anche selezionando Eulero esplicito.

## Operatore e convenzioni

`CPrimitiveEulerSolver::ComputeSpatialResidual` ricostruisce da zero gradienti,
limiter e contributi a ogni chiamata. Restituisce uno `SpatialResult` con
cinque array di stati a quattro componenti:

| Campo | Significato |
| --- | --- |
| `faces` | Somma delle fluttuazioni sulle facce duali interne, già integrate |
| `volume` | Termine interno di volume, già moltiplicato per V |
| `boundary` | Fluttuazione entrante sulle facce fisiche di bordo, già integrata |
| `source` | Sorgente primitiva puntuale S(x_i,t), senza fattore V |
| `rhs` | Operatore semidiscreto `-(faces+volume+boundary)/V + source` |

`LinSysRes` riceve `faces+volume+boundary-V*source`, non `rhs`.
L'assemblaggio in memoria temporanea evita accumuli dalle chiamate precedenti.
Soltanto dopo il controllo di finitezza dell'intero risultato viene aggiornato
LinSysRes; se una validazione o callback lancia, il residuale precedente e W
restano invariati, purché le callback non modifichino esse stesse il solver.
Gli array di lavoro dei gradienti possono essere stati ricalcolati prima
dell'errore: non viene promessa una transazione su tutti gli array ausiliari.

Ogni edge viene valutato una volta con il punto efficace della fase 2,
`Q=(x_i+x_j)/2`, e la normale integrata SU2 orientata i→j. Il kernel riceve
questa normale e applica la misura una sola volta. Si assemblano
`faces[i] += Dminus` e `faces[j] += Dplus`, entrambi con segno positivo.
Non si utilizza il reducer conservativo `+Flux/-Flux`.

L'override virtuale `Upwind_Residual` chiama l'operatore completo senza callback
di sorgente/bordo; può quindi essere usato direttamente con MARKER_FAR.
`Source_Residual` e `Centered_Residual` restano bloccati: il primo non deve
aggiungere una seconda volta sorgenti già incluse e il secondo non definisce
una variante centrata di questo metodo. Timestep e integratori restano bloccati
con `std::logic_error`. La factory conserva `INTEGRATION_TYPE::NONE`.

## Gradienti, limiter e positività

Sono riutilizzati direttamente i template generici SU2:

- `computeGradientsGreenGauss` oppure `computeGradientsLeastSquares`, sulle
  colonne 0–3 di Solution, con indice della prima velocità 1;
- `computeLimiters`, sulle stesse quattro colonne e sul gradiente di
  ricostruzione, passando esplicitamente `GetKind_SlopeLimit_Flow()`.

Il puntatore solver passato ai template è null: non servono comunicazioni in
questo ambito seriale, senza halo e periodicità. Non si passa da CEulerVariable
o dal suo diverso layout delle primitive. Non si usa `GetKind_SlopeLimit()`,
che normalmente viene impostato dall'integrazione tramite SetGlobalParam.

`NUM_METHOD_GRAD_RECON` governa il gradiente usato **sia nelle facce sia nel
volume**. `NUM_METHOD_GRAD` rimane un'opzione SU2 distinta e non viene
sovrascritta. Se la base alloca array distinti, il gradiente generico non è
ricalcolato inutilmente: l'operatore usa esplicitamente quello di ricostruzione.
Sono supportati GREEN_GAUSS, LEAST_SQUARES e WEIGHTED_LEAST_SQUARES. I test
affini verificano LS e WLS; non si attribuisce esattezza affine generale a GG
vicino ai bordi, dove la routine SU2 usa valori nodali.

I limiter supportati sono NONE, BARTH_JESPERSEN e VENKATAKRISHNAN. Le formule,
VENKAT_LIMITER_COEFF e la lunghezza di riferimento rimangono quelle native SU2;
in particolare il ramo denominato BARTH_JESPERSEN usa la propria
specializzazione SU2, senza sostituirla con una formula esterna. Gli altri
limiter vengono rifiutati in modo esplicito, prima di chiamare routine che
richiedono dati non predisposti o limitazione per edge.

La ricostruzione è affine pura: `MUSCL_KAPPA_FLOW` deve essere zero (default).
Valori diversi vengono rifiutati per evitare una diversa proiezione nel limiter
rispetto alla ricostruzione effettiva. Con `MUSCL_FLOW=NO` i gradienti sono nulli
e il volume interno è nullo. Questa variante senza ricostruzione non è
genericamente consistente sulle mesh ammesse, come mostra l'audit finale.

Il gradiente raw resta disponibile nell'array ereditato di ricostruzione.
`GetLimitedGradient(i)` espone separatamente `phi[k]*G_raw[k]*theta_i`.
La contrazione di positività usa un unico theta per tutte le componenti del
nodo, verificando rho e p in tutti i punti effettivamente usati: punti efficaci
interni, punti di bordo e, in MEAN_COEFF, centroide duale. Per ciascuna delle
due componenti si usa la soglia `min(1e-12,0.1*W_i[k])`; quando necessario il
limite su theta viene moltiplicato per 0.99 per restare all'interno positivo.
Non viene alterato né corretto W. Stati nodali invalidi sono rifiutati.
`GetPositivityScale(i)` distingue la contrazione dal limiter nativo.

`ReconstructAt(i,Q)` usa l'ultima ricostruzione predisposta. Dopo aver cambiato
W occorre richiamare ComputeSpatialResidual prima di usarlo come diagnostica.
La positività della ricostruzione non garantisce quella di un futuro
aggiornamento temporale: i controlli per stadio appartengono alla fase 5.

## Volume

In MEAN_COEFF, default, si ricava
`Wmean = W_i + G_limited*(x_centroid-x_i)` con l'offset della fase 2.
Tutti i coefficienti di Ax e Ay dipendono da questo stato: rho, u, v, p e
**1/rho_mean**. In NODE_POINT i coefficienti dipendono da W_i.
In entrambi i casi il prodotto usa lo stesso G_limited delle facce:

```text
B_rho = u*rho_x + v*rho_y + rho*(u_x+v_y)
B_u   = u*u_x   + v*u_y   + p_x/rho
B_v   = u*v_x   + v*v_y   + p_y/rho
B_p   = u*p_x   + v*p_y   + gamma*p*(u_x+v_y)
volume = V*B
```

`VolumeProduct` e la ricostruzione affine sono helper indipendenti dalla
geometria in `primitive_euler_reconstruction.hpp`. Non si usa la beta
logaritmica del path come quadratura del volume. Il reciproco della densità
media resta la chiusura approssimata della specifica, non un integrale esatto
di 1/rho sulla cella. Le due varianti volumetriche producono risultati diversi
anche su campi affini, come previsto; non viene rivendicato ordine 2 misurato.

## Bordi e sorgenti

La geometria aggiunge due facce fisiche per ogni edge di bordo: i segmenti
fra ciascun nodo e il punto medio dell'edge primale. Il punto di quadratura è
il centro di ciascuna metà e la normale integrata è uscente, con misura pari
a metà della lunghezza dell'edge. Il verso è determinato dal triangolo
incidente, senza dipendere dall'ordine dei nodi dell'elemento di bordo.
La somma per nodo/marker è verificata contro l'opposto della normale di
vertice SU2, che nel checkout locale punta verso l'interno. Sono controllate
copertura completa e unicità degli edge esterni. Questo trattamento delle
facce **fisiche di bordo** non introduce la variante SEGMENTED delle facce
duali interne: queste rimangono AGGREGATED.

La prima implementazione supporta:

- MARKER_CUSTOM: callback obbligatoria che restituisce lo stato primitivo
  esterno nel punto di quadratura e al tempo richiesto; permette dati esatti
  prescritti e verifiche successive con manufactured solutions.
- MARKER_FAR: stato esterno uniforme pari a PRIMITIVE_INITIAL_STATE, oppure
  callback fornita esplicitamente. Non è ricavato dal Mach/freestream
  dimensionale del solver conservativo.

La callback riceve `(Q,t,W_inside,n_unit_outward,marker)`; deve restituire rho
e p positive e quattro componenti finite. Il contributo al nodo è soltanto
`Dminus(W_inside,W_outside,S_outward)`. Non viene aggiunto un flusso fisico
conservativo né Dplus per un nodo esterno inesistente. Altri tipi di marker,
incluse pareti e simmetrie, sono rifiutati: non sono riutilizzate implicitamente
le condizioni conservative. La geometria continua a rifiutare periodicità.

La callback sorgente è opzionale, di tipo `State(const Direction& x,double t)`.
Viene valutata nei nodi primali, come richiede NODE, e vale zero se assente.
Si tratta dell'infrastruttura per sorgenti primitive; l'audit finale aggiunge
una MMS solo spaziale. Non è stata introdotta una campagna di evoluzione e
convergenza spazio-temporale della fase 6.

Esempio C++ dopo la costruzione della geometria e del solver:

```cpp
auto boundary = [](const PrimitiveEuler2D::Direction& q, double t,
                   const PrimitiveEuler2D::State& inside,
                   const PrimitiveEuler2D::Direction& normal, unsigned short marker) {
  return PrimitiveEuler2D::State{{1.0, 0.5, 0.3, 1.0}};
};
auto spatial = solver.ComputeSpatialResidual(geometry, config, time, boundary);
// spatial.rhs[i] = dW_i/dt; nessun aggiornamento di W è stato eseguito.
```

Geometria e configurazione devono restare quelle originali, invariate durante
la vita del solver; sono snapshot statici e devono vivere più a lungo del solver.
L'API controlla identità e nullità dei puntatori e finitezza del tempo, senza
pretendere di rilevare ogni modifica in-place di oggetti esterni.

## Ambito seriale e thread safety

La chiamata deve avvenire fuori da regioni OpenMP: una chiamata interna viene
rifiutata prima di modificare gli array. L'assemblaggio degli edge è seriale
e deterministico, quindi non introduce scatter concorrenti sui nodi.
Un atomic_flag impedisce due valutazioni simultanee o rientranti sullo stesso
solver; una guardia RAII lo libera anche su eccezione. Non si promette accesso
concorrente ai setter raw, alle diagnostiche o alle callback che mutano W.
MPI distribuito e accelerazione OpenMP dell'operatore non sono implementati.

## Verifiche iniziali della FASE 4

La selezione conserva i 34 casi delle fasi precedenti e aggiunge 12 casi di
fase 4. Il test storico del blocco di evoluzione controlla ora i sette metodi
ancora indisponibili; Upwind_Residual è verificato separatamente tramite
puntatore CSolver* sul percorso spaziale effettivo.

Copertura nuova:

- chiusura geometrica per ogni duale, perimetro e momenti affini di bordo;
- residuale completo nullo per stato costante con GG, LS e WLS;
- campo affine su mesh irregolare, salti ricostruiti nulli con LS/WLS,
  coefficienti e volume verificati indipendentemente nelle due varianti,
  anche con gamma=1.6 e NUM_METHOD_GRAD distinto;
- tutte e quattro le equazioni espanse del prodotto di volume;
- primo ordine con volume nullo, segni di entrambe le fluttuazioni, singolo
  fattore di misura e somma globale non conservativa;
- limiter NONE/BARTH_JESPERSEN/VENKATAKRISHNAN e un caso positivo ai nodi che
  richiede contrazione dei gradienti per evitare ricostruzioni negative;
- callback di bordo nei veri punti di quadratura, normali uscenti unitarie,
  tempo e solo contributo entrante;
- sorgenti nodali con segno e volume corretti, bilanciamento affine in
  NODE_POINT e differenza attesa in MEAN_COEFF;
- ricalcolo senza accumulo, mancata modifica di W, errori senza commit del
  residuale, rientranza rifiutata e riutilizzo dopo eccezione;
- override Upwind attraverso CSolver*, blocco temporale e livello coarse
  rifiutato; invarianza alla rinumerazione di nodi ed elementi;
- rifiuto esplicito di limiter non predisposti, U-MUSCL non affine e simmetrie.

Comandi dalla radice del repository:

```sh
./ninja -C build -j 2 SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-phase4-tests
/private/tmp/su2-phase4-tests/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-phase4-cli --cli
python3 UnitTests/run_primitive_euler_tests.py --asan --output-dir /private/tmp/su2-phase4-asan
ASAN_OPTIONS=detect_leaks=0 /private/tmp/su2-phase4-asan/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
```

Build CFD e suite normale completate con uscita 0: **46 casi e 144729
asserzioni superate**. I quattro controlli CLI passano: primitivo uscita 1
attesa con il nuovo messaggio di fase 4 e senza Exit Success; EULER conservativo
completa l'iterazione con uscita 0; harmonic balance e specie escono con 1
per il rifiuto specifico di CConfig. Il wrapper CLI termina con uscita 0.

La build ASan usa lo stesso perimetro parziale della fase 3: test, factory,
classi primitive e basi CVariable/CScalarVariable strumentati; librerie native
SU2/MPI non strumentate, annotazioni libc++ uniformate sul Mac tramite
`_LIBCPP_HAS_NO_ASAN`. LeakSanitizer non è supportato su questa piattaforma
e resta disabilitato. Il risultato ASan non può certificare assenza di leak
o sostituire una build completamente sanitizzata.
La build e l'esecuzione ASan sono state completate con uscita 0: **46 casi e
144729 asserzioni superate**, senza segnalazioni AddressSanitizer.

Log in `/private/tmp/`: `su2-phase4-build.log`, `su2-phase4-tests-build.log`,
`su2-phase4-results.log`, `su2-phase4-cli.log`, `su2-phase4-asan-build.log`,
`su2-phase4-asan-results.log`; dettagli CLI in `su2-phase4-cli/`.
Warning non bloccanti: pragma GCC ignorati da Clang, costruttori default
implicitamente eliminati in classi SU2 esistenti, librerie duplicate e archivi
senza simboli. MPI seriale segnala bind TCP negato nel sandbox; le verifiche
terminano comunque correttamente. Le configurazioni temporali dei test
precedenti segnalano l'adeguamento WINDOW_START_ITER/RESTART_ITER.

Non sono stati eseguiti la suite SU2 completa, test MPI distribuiti o misure
d'ordine spazio-temporale. Non sono implementati SEGMENTED, DUAL_MEAN,
SSP-RK2, pareti, periodicità o mesh mobili. Non è stato eseguito ninja install;
usare la build, non presumere che SU2_bin/bin sia aggiornato.

## File interessati

- `primitive_euler_reconstruction.hpp`: helper di ricostruzione e prodotto di volume.
- `primitive_euler_geometry.hpp`: quadratura e validazione dei bordi fisici.
- `CPrimitiveEulerSolver.hpp/.cpp`: pipeline spaziale, diagnostiche, controlli
  di configurazione, positività e commit del residuale.
- `CDriver.cpp`, `CIterationFactory.cpp`, `config_template.cfg`: messaggi e
  descrizioni coerenti con il nuovo checkpoint, mantenendo bloccata la CLI.
- `primitive_euler_residual_tests.cpp`: dodici casi spaziali;
  `primitive_euler_solver_tests.cpp`: aggiornamento del contratto virtuale.
- `primitive_euler_mesh_fixture.hpp`: scelta del tipo di marker nei test.
- `UnitTests/meson.build`, `run_primitive_euler_tests.py`: registrazione della
  nuova suite e messaggio CLI aggiornato.
- `primitive_euler_phase3.md`: rinvio al checkpoint successivo;
  questo documento: decisioni, interfacce, risultati e limiti.

Non sono stati modificati CEulerSolver, CEulerVariable o le formule del kernel
path della fase 1. Nessun commit o push creato.

## Verifiche finali e hardening prima della FASE 5

Questo ciclo aggiunge nove test case e aggiorna il risultato a **55 casi,
148037 asserzioni superate**, uscita 0 nella suite normale. Tutti i 46 casi
precedenti restano presenti. I test che riproducono limiti matematici verificano
e registrano quei limiti: il loro superamento non certifica un metodo consistente.
Non sono stati introdotti integratori, driver completi o nuove condizioni fisiche.

### Gamma effettivo: nessuna incoerenza trovata

L'ispezione di `Common/src/CConfig.cpp`, nel postprocessing `if (standard_air)`,
conferma che STANDARD_AIR riporta gamma a 1.4 e R a 287.058. IDEAL_GAS conserva
il gamma configurato. CPrimitiveEulerVariable memorizza `config->GetGamma()`;
PrepareReconstruction, ReconstructAt e ComputeSpatialResidual leggono lo stesso
getter dalla configurazione originale, che per contratto resta invariata.
PhysicalMatrix, VolumeProduct e ComputeFluctuations ricevono quel valore come
argomento; non contengono una costante 1.4 nascosta.

Il nuovo test usa STANDARD_AIR con GAMMA_VALUE=1.6 (effettivo 1.4) e IDEAL_GAS
con 1.6 e 1.25. Un campo non affine rende non nulli sia i salti ricostruiti
sia il volume. Si controllano suono, energia, entalpia, conversione diagnostica,
matrice fisica e contributi faces/volume/boundary con il gamma effettivo.
Si verifica anche che un gamma diverso darebbe fluttuazioni diverse: il caso
non è insensibile al parametro. Non è stata necessaria una modifica al gamma.

### Geometria con marker multipli

La fixture conserva il comportamento precedente e ora accetta una lista di
marker con edge espliciti. Il nuovo test usa i lati adiacenti sud/est e
nord/ovest di un rettangolo triangolato: i nodi 0 e 2 appartengono a entrambi
i marker. Tutte le otto mezze facce vengono confrontate con valori analitici
di proprietario, tag, punto e normale. Si verificano copertura unitaria,
mezze lunghezze, chiusura di ciascun duale, somme native per nodo/marker e
callback distinte per i due tag, senza raddoppiare contributi agli angoli.

Il controllo viene ripetuto invertendo l'input e invertendo anche gli elementi
di bordo **dopo** la correzione d'orientazione SU2, conservando le metriche native:
la costruzione delle normali rimane corretta senza affidarsi all'ordinamento.
Un altro caso rifiuta edge di bordo inesistenti, duplicati e nodi fuori intervallo.

Sono stati aggiunti controlli espliciti su tutti gli indici FindEdge usati nei
passaggi di bordo prima dell'accesso ai vettori, sui vertici di bordo e sulla
misura fisica finita e positiva. Sono protezioni locali: non è stato riprodotto
un crash precedente su una geometria nativa valida. Le validazioni delle fasi
precedenti già escludevano molti di questi stati incompleti.

### Rinumerazione con MUSCL e positività

LS e WLS sono esercitati con MUSCL attivo, campo non affine su mesh irregolare,
permutazione non banale dei nodi e ordine degli elementi invertito. Un secondo
campo forza theta<1. Dopo la permutazione inversa si confrontano gradienti
limitati, theta, ricostruzioni interne e fisiche di bordo, faces, volume,
boundary, rhs e LinSysRes. Tolleranze relative 2e-11 e assolute 5e-13 ammettono
il riordinamento delle poche somme e le cancellazioni del problema LS;
non si pretende identità bitwise. I confronti passano senza aggiustamenti
successivi delle tolleranze e non mostrano dipendenza numerica dall'ordinamento.

Il test di positività esteso usa rho e p dell'ordine di 1e-12, 1e-100 e 1e-280,
con variazioni di un fattore 100 che richiedono contrazione. Verifica theta
finito in [0,1], positività e finitezza di tutti i punti interni, fisici e
dei centroidi MEAN_COEFF e conservazione esatta degli stati nodali.

È stato identificato un problema nella vecchia espressione: per W=1e308 e
candidate=-1e308, entrambi finiti, `W-candidate` trabocca a infinito e il
fattore diventa zero invece di circa 0.495. `ReconstructionPositivityFactor`
calcola ora lo stesso rapporto con quantità scalate; materializza le divisioni
con variabili volatile per impedirne la rifattorizzazione nella differenza
originale da parte di `-ffast-math`. Il test riproduce l'overflow della vecchia
differenza e verifica 0.495, oltre a un rapporto su valori molto piccoli.

La soglia `min(1e-12,0.1*W)` e il fattore 0.99 rimangono invariati. Se la soglia
non è rappresentabile ed è zero, un candidato non positivo richiede theta=0:
si mantiene il valore nodale positivo con gradiente nullo, evitando una
ricostruzione arrotondata a zero. Il controllo con denorm_min è condizionato
al fatto che la piattaforma lo confronti strettamente positivo; i test a
1e-280 restano normali e vengono sempre eseguiti. Non si promette che ogni
stato fisico estremo sia numericamente utilizzabile: 1/rho, velocità del suono
o prodotti possono ancora eccedere il range double ed essere rifiutati.

La regressione del fattore estremo è un test del calcolo scalare, non una
simulazione completa con rho=1e308. Gli stati completi estremi esercitati
dinamicamente sono quelli piccoli elencati sopra.

### Interfaccia ed eccezioni

Il nuovo test provoca NaN, +Inf e -Inf in ogni componente della callback di
bordo, della sorgente e dello stato nodale scritto attraverso l'accesso raw;
verifica anche rho/p nulle o negative e callback che lanciano runtime_error.
Dopo ogni eccezione controlla LinSysRes invariato e ripete con successo la
valutazione valida, dimostrando il rilascio della guardia evaluating nei percorsi
esercitati. Le sorgenti possono avere segno qualsiasi, ma devono essere finite.
I test precedenti mantengono la verifica della rientranza rifiutata.

ReconstructAt conserva il contratto esistente: prima della prima preparazione
usa il gradiente iniziale nullo e restituisce W_i; dopo una scrittura nodale
combina il W corrente con l'ultimo gradiente disponibile. Non rileva automaticamente
la ricostruzione obsoleta. Il test documenta questo comportamento e verifica che
una nuova chiamata a ComputeSpatialResidual ricalcoli tutto. Nessuna cache di
stato o modifica invasiva dei setter ereditati è stata introdotta.

### Audit algebrico e MMS solo spaziale

**A — Correttezza algebrica.** Il ciclo interno visita ogni edge una volta;
Dminus e Dplus sono sommati con entrambi i segni positivi. Il kernel applica
la misura una volta, il volume viene moltiplicato una volta per V, la sorgente
viene campionata al nodo e sottratta nel residuale integrato, e ogni mezza faccia
di bordo contribuisce soltanto al proprietario. I test iniziali e quelli con
marker multipli esercitano questi fatti. Facce e volume condividono G_limited.
Non sono state cambiate queste formule.

**B — Consistenza spaziale.** La nuova MMS stazionaria campiona:

```text
rho = 2 + 0.15*sin(x)*cos(y)
u   = 0.5 + 0.1*cos(1.3*x)*sin(0.7*y)
v   = 0.3 + 0.08*sin(0.8*x)*cos(1.1*y)
p   = 3 + 0.2*cos(0.9*x)*cos(1.2*y)
```

Si usa S=Ax(W)*W_x+Ay(W)*W_y con derivate analitiche scritte esplicitamente
nel test: il rhs esatto è zero. I bordi ricevono lo stato esatto nei veri punti
di quadratura. Nessuna integrazione temporale e nessuna derivata numerica nel
riferimento. Sul quadrato [0,1]^2 si usano N=8,16,32,64 suddivisioni per lato,
due triangoli per quadrato con diagonale fissa. La seconda famiglia sposta
solo i nodi interni di `0.18*h*sin(2.1*i+1.3*j)` in x e
`0.16*h*cos(1.7*i-0.9*j)` in y; i bordi restano rettilinei.

L'errore è `sqrt(sum_group V*||rhs||^2 / sum_group V)`, distinto in tutti i nodi,
nodi interni e nodi sul bordo fisico. Gli ordini sono log2(E_N/E_2N).
Per la prima mesh il CSV usa 0 come segnaposto dell'ordine non definito.
La MMS è una misura diagnostica: il test controlla errori finiti e stampa i
risultati, senza imporre artificiosamente convergenza o ordine 2.
I dati completi delle 24 valutazioni sono conservati in
[primitive_euler_phase4_spatial_mms.csv](primitive_euler_phase4_spatial_mms.csv).

| Famiglia | Variante | E interno N=64 | p interno 32→64 | E bordo N=64 | p globale 32→64 |
| --- | --- | ---: | ---: | ---: | ---: |
| Regolare | MUSCL NO | 5.978685e-3 | 1.0015 | 2.644312e-1 | 0.4872 |
| Regolare | NODE_POINT | 3.545778e-4 | 1.5234 | 2.635450e-3 | 1.5002 |
| Regolare | MEAN_COEFF | 3.545778e-4 | 1.5234 | 2.630728e-3 | 1.5006 |
| Distorta | MUSCL NO | 5.375299e-2 | -0.0073 | 2.653220e-1 | 0.2351 |
| Distorta | NODE_POINT | 9.791404e-4 | 1.0855 | 2.727544e-3 | 1.1730 |
| Distorta | MEAN_COEFF | 9.804418e-4 | 1.0854 | 2.721559e-3 | 1.1725 |

Con MUSCL NO l'errore interno sulla famiglia distorta ristagna: 0.066436,
0.055122, 0.053482, 0.053753. Sul bordo ristagna anche nella famiglia regolare;
la sua misura totale diminuisce con h, quindi il solo errore globale può
mascherare un difetto locale. Con MUSCL attivo gli errori diminuiscono nelle
famiglie misurate, che costituiscono evidenza sperimentale limitata di
consistenza, non una dimostrazione per mesh arbitrarie.

Un ulteriore caso isola un **difetto di consistenza riproducibile**, senza
dipendere dalla MMS non affine. Nel nodo interno di una patch asimmetrica si
usa `rho=2+0.1*x+0.2*y`, u=0.5, v=0.3, p=3. Il prodotto differenziale esatto
è `(0.11,0,0,0)`. Per questo contatto il kernel si riduce esattamente, sulla
prima componente, a `Dminus=min(u*Sx+v*Sy,0)*(rho_j-rho_i)`.
Le normali uscenti e i salti su tutti i vicini danno R_rho/V=0.074 anziché 0.11.
Con S_rho=0.11, il rhs del nodo è quindi **0.036**, per h=1,1/2,1/4,1/8.
La patch è scalata attorno allo stesso nodo interno, fissato nell'origine;
normali, salti e volumi scalano rispettivamente come h,h,h^2. Il difetto non
si riduce. Il nodo non ha contributi di bordo. Con MUSCL attivo lo stesso
campo affine è riprodotto e il rhs è nullo entro arrotondamento.

Questa è una limitazione della discretizzazione senza ricostruzione su uno
stencil asimmetrico, non un errore corretto scambiando i segni delle fluttuazioni.
Non è stata aggiunta una correzione geometrica o un termine di volume diverso
da quello concordato. La prova locale è più forte della sola mancata convergenza
su quattro mesh: esibisce uno scarto O(1) su una famiglia di stencil simili.

**C — Accuratezza.** I risultati non dimostrano ordine 2: con MUSCL attivo
gli ordini globali finali sono circa 1.50 sulla famiglia regolare e 1.17 sulla
distorta; gli ordini interni includono i nodi adiacenti al bordo. Non si
deduce l'ordine delle soluzioni evolute dall'ordine del residuale: non è stata
eseguita alcuna evoluzione temporale. La campagna di fase 6 resta necessaria.

### Comandi ed esiti dell'audit finale

Comandi effettivamente eseguiti, dalla radice:

```sh
./ninja -C build -j 2 SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-phase4-hardening
/private/tmp/su2-phase4-hardening/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-phase4-hardening-cli --cli
python3 UnitTests/run_primitive_euler_tests.py --asan --output-dir /private/tmp/su2-phase4-hardening-asan
ASAN_OPTIONS=detect_leaks=0 /private/tmp/su2-phase4-hardening-asan/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
```

Durante l'aggiunta dei test è stato ricompilato e ricollegato il solo translation
unit `primitive_euler_residual_tests.cpp` con gli stessi argv del runner,
estratti dai log con shlex.split e passati a subprocess senza shell. La build
CFD termina con uscita 0; la suite normale supera 55 casi/148037 asserzioni,
uscita 0. I quattro controlli CLI passano con i consueti codici 1/0/1/1 e il
wrapper termina con 0. L'eseguibile usato è quello della build nativa.

La strumentazione ASan conserva i limiti già dichiarati: archivi SU2/MPI
parzialmente non strumentati, annotazioni libc++ disabilitate uniformemente,
LeakSanitizer non supportato e detect_leaks=0. Nessun nuovo controllo con
macOS leaks viene dichiarato. La suite SU2 completa e MPI distribuito non
sono stati eseguiti. I warning non bloccanti rimangono quelli descritti sopra.
Anche la build e l'esecuzione ASan terminano con uscita 0 e **55 casi/148037
asserzioni superate**, senza segnalazioni AddressSanitizer. Le 24 righe MMS
e i quattro difetti di contatto coincidono alle cifre stampate con la suite
normale.

Log `/private/tmp/su2-phase4-hardening-{build,tests-build,rebuild,results,cli,asan-build,asan-rebuild,asan-results}.log`;
i file CLI dettagliati sono in `/private/tmp/su2-phase4-hardening-cli/`.
Questi log sono temporanei; il runner, i test e il CSV restano nel repository.

### File modificati in questo ciclo e valutazione

- `primitive_euler_reconstruction.hpp`: rapporto di positività scalato e caso
  limite della soglia non rappresentabile.
- `CPrimitiveEulerSolver.cpp`: usa l'helper per i vincoli di ricostruzione.
- `primitive_euler_geometry.hpp`: controlli locali di indici, vertici e normali.
- `primitive_euler_mesh_fixture.hpp`: marker multipli e inversione degli edge.
- `primitive_euler_residual_tests.cpp`: nove nuovi casi, inclusi MMS spaziale
  e riproduzione del difetto di consistenza del contatto affine.
- Questo documento e il CSV: risultati, distinzione A/B/C e limiti.

**Valutazione storica del baseline prima della correzione quadratica: FASE 4
non pronta per una chiusura matematica e un passaggio incondizionato
alla FASE 5.** L'infrastruttura è utilizzabile e verificata nei percorsi testati,
ma il difetto di consistenza di MUSCL NO deve essere risolto mediante una
scelta numerica esplicita, oppure quella variante deve essere esclusa dall'ambito
accettato del solver evolutivo. Anche limitandosi a MUSCL attivo, il secondo
ordine resta da verificare e non è attestato da questo audit. Non sono state
cambiate le formule per ottenere un ordine prefissato; nessuna fase 5, nessun
commit e nessun push.
