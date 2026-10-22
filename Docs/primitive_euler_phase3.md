# Euler 2D primitivo — FASE 3: classi, configurazione e factory

Questa fase aggiunge lo stato primitivo e la classe solver dedicata, senza
modificare `CEulerSolver` o `CEulerVariable`. Restano confermati AGGREGATED,
NODE ed Eulero esplicito per la prima evoluzione numerica. Residuale, termine
interno di volume, condizioni al bordo, output e avanzamento temporale saranno
implementati nelle fasi successive: questo non è ancora un solver utilizzabile
per simulazioni.

## Configurazione

```ini
SOLVER= PRIMITIVE_EULER
MGLEVEL= 0
TIME_DISCRE_FLOW= EULER_EXPLICIT
PRIMITIVE_REPRESENTATION= NODE
PRIMITIVE_FACE_QUADRATURE= AGGREGATED
PRIMITIVE_VOLUME_QUADRATURE= MEAN_COEFF
PRIMITIVE_INITIAL_STATE= (1.0, 0.5, 0.3, 1.0)
```

Le quattro componenti iniziali sono **rho, u, v, p**, con valori adimensionali
coerenti tra loro e con `GAMMA_VALUE`. Il nuovo stato uniforme è specificato
direttamente; non deriva dalle impostazioni dimensionali di freestream/Mach
usate da Euler conservativo. Default: `(1,0.5,0.3,1)`. Rho e pressione devono
essere positive e tutte le componenti finite.
Per scegliere gamma diverso da 1.4 impostare `FLUID_MODEL=IDEAL_GAS`:
`STANDARD_AIR` mantiene le costanti dell'aria previste da SU2.

Rappresentazione e quadratura delle facce hanno per ora un solo valore
ammissibile, NODE e AGGREGATED. Per il futuro volume sono selezionabili
MEAN_COEFF (default) e NODE_POINT; in questa fase viene solo validata e
memorizzata la scelta, senza calcolare un residuale. MEAN_COEFF utilizzerà
gli offset dei centroidi già disponibili dalla FASE 2.

Le opzioni standard `NUM_METHOD_GRAD`, `NUM_METHOD_GRAD_RECON`,
`SLOPE_LIMITER_FLOW` e i relativi parametri restano quelle SU2. La nuova classe
prepara gli array necessari, ma non calcola ancora gradienti o limiter.
SEGMENTED e DUAL_MEAN vengono rifiutati, non sostituiti automaticamente.

La costruzione richiede gas ideale, griglia triangolare 2D statica, seriale,
non periodica e livello più fine, senza multigrid, restart, Newton-Krylov,
accoppiamenti multizona, turbolenza, specie o assialsimmetria. Sono esclusi
dual-time e harmonic balance. Il supporto AD non è implementato/verificato.
La configurazione viene controllata alla costruzione del solver; errori di
setup generano `std::invalid_argument`. La validazione geometrica rimane in
`CAggregatedGeometry`.

## Stato e proprietà

`CPrimitiveEulerVariable` deriva dalla memoria generica `CScalarVariable`,
usata qui come contenitore di un campo a quattro componenti con gradienti e
limiter. Tale base non introduce conversioni termodinamiche o flussi di uno
schema conservativo. `Solution` e `Primitive` accedono allo **stesso W**:

| Indice | Componente evoluta |
| --- | --- |
| 0 | rho |
| 1 | u |
| 2 | v |
| 3 | p |

`GetPrimitiveState` e `SetPrimitiveState` consentono lettura/scrittura tipizzata;
la scrittura controlla l'ammissibilità prima di modificare i dati. I setter
generici ereditati da SU2 restano accessi di basso livello: la futura evoluzione
dovrà controllare la positività dopo ogni aggiornamento.

Gli accessori densità, velocità e pressione leggono direttamente W. Suono,
energia ed entalpia sono ricavati dallo stato corrente, senza cache da
sincronizzare. `GetEnergy` segue la convenzione SU2 dell'energia totale
**specifica**. `GetConservativeState` restituisce invece
`(rho,rho*u,rho*v,E)` con energia **per volume**. Nessuna conversione cambia
la rappresentazione memorizzata.

Sono inizializzati stato precedente, eventuali stati temporali previsti dalla
base, gradienti nulli, limiter unitari e timestep nullo. Preparare questi array
non equivale a implementare una ricostruzione o un integratore.

## Solver e collegamento a SU2

`CPrimitiveEulerSolver` deriva direttamente da `CSolver` e possiede con RAII
il nuovo stato e la fotografia geometrica `CAggregatedGeometry`. Imposta
quattro incognite, gli array di diagnostica e i buffer dei residuali; registra
correttamente il puntatore ai nodi usato dalla base.

`CSolverFactory::CreateSolverContainer` riconosce PRIMITIVE_EULER e costruisce
solo FLOW_SOL. Il nuovo sottotipo ha metadato di integrazione NONE, coerente
con una classe ancora senza evoluzione. I tipi già esistenti mantengono i
precedenti identificatori del solver principale.

Per evitare un'esecuzione apparentemente riuscita con i metodi vuoti della
base, i metodi di residuale e integrazione della nuova classe generano
`std::logic_error`. Il driver SU2 rifiuta inoltre esplicitamente una simulazione
PRIMITIVE_EULER prima di creare output o numeriche conservative. Questa
protezione andrà rimossa quando le fasi successive avranno un percorso
numerico completo e verificato; non basta eliminare il messaggio per ottenere
un solver valido. La factory può comunque essere esercitata direttamente nei
test di inizializzazione.

## File e verifica

- `Common/include/option_structure.hpp`, `Common/include/CConfig.hpp` e
  `Common/src/CConfig.cpp`: identificatore, parsing e accessori delle opzioni,
  classificazione come problema fluido comprimibile e parametri di flusso.
- `CPrimitiveEulerVariable.hpp/.cpp` e `CPrimitiveEulerSolver.hpp/.cpp`:
  classi dedicate, con proprietà termodinamiche e controlli del contratto.
- `CSolverFactory.hpp/.cpp`, `CDriver.cpp`, `CIterationFactory.cpp`:
  registrazione e protezione delle vie di esecuzione non implementate.
- `SU2_CFD/src/meson.build`, `UnitTests/meson.build`, `config_template.cfg`:
  registrazione delle sorgenti/test e descrizione delle nuove opzioni.
- `UnitTests/SU2_CFD/primitive_euler_mesh_fixture.hpp`: fixture SU2 già usata
  dalla geometria, estratta per condividerla senza duplicare il preprocessing.
- `primitive_euler_solver_tests.cpp`: test di stato, conversione, configurazione,
  factory, metadati, configurazioni rifiutate ed evoluzione non disponibile.
  I test geometrici precedenti conservano la stessa fixture e gli stessi casi.

La selezione completa pertinente è
`[primitive-geometry],[primitive-euler],[primitive-solver]`. Il confronto con
`CEulerVariable` verifica la sola conversione dello stato iniziale, non la
convergenza delle due discretizzazioni verso la medesima soluzione smooth.

Comandi di build e verifica usati dalla radice del repository:

```sh
./ninja -C build -j 2 SU2_CFD/src/SU2_CFD
python3 /private/tmp/su2-build-phase3-tests.py
/private/tmp/su2-phase3-tests '[primitive-geometry],[primitive-euler],[primitive-solver]'
```

Il driver dedicato usa i flag e le macro in `build/compile_commands.json` e le
librerie SU2 aggiornate della build. Non richiede di abilitare l'intera suite
con `enable-tests=true`; con il driver Catch2 standard disponibile si può
usare la stessa selezione di tag.

## Esito della consegna

Compilazione e collegamento di `build/SU2_CFD/src/SU2_CFD` completati, con
codice di uscita 0. Sono stati eseguiti tutti i test della selezione:
**29 test case, 141569 asserzioni superate, codice di uscita 0**. Comprendono
i precedenti 24 test e cinque nuovi test di integrazione del solver. Non è
stata eseguita l'intera suite SU2.

Sono stati inoltre eseguiti due controlli CLI, usando una mesh a triangolo
singolo, con `MGLEVEL=0`, `TIME_DISCRE_FLOW=EULER_EXPLICIT`,
`CONV_NUM_METHOD_FLOW=ROE`, `MACH_NUMBER=0.2`, `CFL_NUMBER=0.1`,
`MUSCL_FLOW=NO` e `ITER=1`:

- PRIMITIVE_EULER: uscita 1 **attesa**, con il messaggio specifico che vieta
  la simulazione incompleta prima di accedere a output o flussi conservativi.
- EULER: una iterazione completata e uscita 0. È un controllo di avvio e
  compatibilità del percorso esistente, non una campagna di regressione CFD.

Il comando per questi controlli è
`python3 /private/tmp/su2-phase3-cli-check.py`; il relativo script verifica
automaticamente il messaggio del prototipo e l'uscita corretta di Euler.
Log della suite: `/private/tmp/su2-phase3-test-results.log`. Log della build
del driver di test: `/private/tmp/su2-phase3-test-build.log`. I driver e gli
script temporanei non fanno parte del repository.

Gli avvisi non bloccanti riguardano pragma GCC non riconosciuti da Clang,
librerie duplicate e oggetti senza simboli negli archivi SU2. Durante la
rigenerazione Meson è stato segnalato anche un accesso negato a un git config
di un sottomodulo nel sandbox; compilazione e link sono comunque riusciti.
L'inizializzazione MPI seriale segnala il bind TCP negato dal sandbox, ma i
test e la simulazione Euler conservativa terminano con codice zero.

Non è stato eseguito `ninja install`: `SU2_bin/bin/SU2_CFD` non viene aggiornato
da questa consegna. Per usare il codice compilato, riferirsi all'eseguibile in
`build/SU2_CFD/src/`. Non sono ancora disponibili simulazioni primitive, MMS,
studi d'ordine o confronti fra velocità di shock; tali risultati richiedono
le fasi numeriche successive.

## Verifiche finali di robustezza prima della FASE 4

Questa verifica integra e aggiorna i risultati precedenti. La selezione conserva
tutti i 29 casi precedenti e aggiunge cinque casi: **34 test case e 142380
asserzioni superate, uscita 0**, sia nel driver normale sia nel driver ASan
parzialmente strumentato descritto sotto. Nessuna funzionalità della fase 4 è
stata introdotta; NODE, AGGREGATED e il default MEAN_COEFF restano invariati.

### Factory: ispezione del codice e correzione circoscritta

`CreateSolverContainer` alloca `new CSolver*[MAX_SOLS]()`; le parentesi
inizializzano tutti gli slot a null. Il ramo PRIMITIVE_EULER costruisce soltanto
FLOW_SOL: non esistono solver precedentemente costruiti in questo ramo.
L'assegnazione allo slot avviene solo al ritorno di `CreateSubSolver`.
Il `catch (...)` libera il contenitore con `delete[] solver` e rilancia con
`throw`, preservando anche eccezioni diverse da `std::invalid_argument`.

Nel costruttore primitivo, `metrics` e `nodes` sono `unique_ptr`; se la
costruzione fallisce, i membri già costruiti e la base vengono distrutti e la
new-expression libera lo spazio dell'oggetto incompleto. L'ispezione di
CSolver/CVariable/CScalarVariable conferma che gli array interessati sono
contenitori proprietari e che i puntatori nativi liberati dalla base partono
da null. La base non possiede una seconda copia dei nodi primitivi: il puntatore
di accesso è un alias. Queste sono conclusioni sull'ownership dal sorgente,
non una misura dinamica di tutte le possibili allocazioni fallite.

È emersa una lacuna concreta nel percorso successivo alla costruzione:
`allocatedSolvers[genericSolver] = metaData` può allocare e lanciare, lasciando
il solver appena costruito senza proprietario. La correzione minima introduce
un `unique_ptr<CSolver> primitiveOwner` in `CreateSubSolver`, usato solo per il
nuovo ramo e rilasciato dopo la registrazione riuscita. Se l'inserimento lancia,
il solver è distrutto prima che il catch esterno liberi il contenitore.
L'API pubblica resta invariata. Non è stato iniettato un `bad_alloc` nella
mappa: questo percorso è verificato mediante analisi dell'ownership, non con
un guasto dinamico forzato.

Il nuovo test `[primitive-factory-unwind]` passa 32 volte dalla factory con
pressione iniziale negativa: la geometria aggregata è già stata allocata quando
la validazione delle variabili lancia. Poi costruisce correttamente un solver
valido e verifica che la geometria della fixture sia ancora utilizzabile.
Il test isolato supera 38 asserzioni; da solo non dimostra assenza di leak.
Gli altri rami della factory, alcuni dei quali costruiscono più solver mediante
puntatori raw, non sono stati rifattorizzati: questa analisi non certifica la
exception safety globale della factory SU2.

### Copertura dinamica degli stati e degli array

- Il setter rifiuta 16 stati: densità/pressione nulle o negative, NaN e
  entrambi gli infiniti in ciascuna delle quattro componenti. Ogni tentativo
  parte da uno stato valido non banale e verifica esattamente tutte le
  componenti invariate. I tre accessori tipizzati di stato/conversione rifiutano
  gli indici nPoint e `ULONG_MAX` con `std::out_of_range`.
- Due stati non banali verificano indipendentemente suono, energia totale
  specifica, entalpia e conversione diagnostica in U con IDEAL_GAS e gamma
  1.25/1.6. STANDARD_AIR con GAMMA_VALUE=1.6 produce gamma effettivo 1.4,
  coerentemente con il postprocessing di CConfig. Il cambio di stato aggiorna
  subito le proprietà; W resta la rappresentazione memorizzata.
- Una mesh triangolare a cinque nodi verifica indipendenza dei nodi, identità
  dei puntatori Solution/Primitive e modifiche bidirezionali. Verifica array
  5x4, gradienti 5x4x2 inizialmente nulli, limiter inizialmente unitario,
  stato precedente e timestep nodale nullo. I buffer temporali hanno dimensione
  zero in steady e 20 in TIME_STEPPING; sono letti solo quando allocati.
  I valori iniziali sono verificati per tutti i nodi e le componenti.
- Le interfacce raw ereditate restano mutabili e possono aggirare il setter:
  un test scrive deliberatamente pressione negativa attraverso Solution e
  verifica il rifiuto della conversione diagnostica. La futura evoluzione
  dovrà controllare la positività dopo ogni aggiornamento numerico.
- Tutti gli otto metodi di evoluzione elencati nel contratto sono chiamati
  attraverso `CSolver*`: ciascuno lancia `std::logic_error` con il messaggio
  specifico atteso, conservando gli stati e un residuale sentinella non nullo.

### Configurazioni e integrazione

La validazione del solver rifiuta dinamicamente implicito, MGLEVEL>0,
restart, dual-time, multizona, gas VW_GAS, turbolenza SA, assialsimmetria,
DUAL_MEAN, SEGMENTED e quadratura di volume sconosciuta. Le opzioni primitive
sono stringhe: in questi casi il parser le accetta e il rifiuto è nel solver.
MGLEVEL e multizona sono impostati anche tramite setter pubblici per verificare
la sola validazione, senza costruire geometrie multigrid o array accoppiati.
Lo stato iniziale non valido è invece rifiutato dalla costruzione delle
variabili. I vincoli geometrici restano coperti dalla suite di fase 2.

Harmonic balance e specie sono esercitati in subprocessi CLI: è CConfig a
terminare con uscita 1, prima di raggiungere la validazione del solver.
Non vengono convertiti artificialmente in eccezioni Catch2. I corrispondenti
controlli difensivi del solver sono stati ispezionati ma non raggiunti
dinamicamente. Analogamente, i controlli per Newton-Krylov, radiazione e
nZone diverso da uno sono stati ispezionati, senza introdurre fixture
accoppiate o modifiche alla configurazione per forzarli.
GREEN_GAUSS, WEIGHTED_LEAST_SQUARES e VENKATAKRISHNAN sono conservati nelle
rispettive opzioni standard e verificati anche costruendo il solver.

I quattro smoke test usano `build/SU2_CFD/src/SU2_CFD` e una mesh minima valida:

| Caso | Uscita | Controllo |
| --- | --- | --- |
| PRIMITIVE_EULER | 1 attesa | Blocco specifico nel driver, nessun Exit Success |
| EULER | 0 | Iterazione 0 effettivamente riportata e Exit Success |
| PRIMITIVE_EULER harmonic balance | 1 attesa | Rifiuto specifico di CConfig |
| PRIMITIVE_EULER specie | 1 attesa | Rifiuto specifico di CConfig |

Sono controlli di avvio, non una validazione CFD. L'eseguibile installato in
SU2_bin non è stato aggiornato.

### Comandi riproducibili e strumenti di memoria

Il nuovo runner persistente sostituisce la necessità degli script temporanei
citati sopra: ricava flag e dipendenze dalla build Meson/Ninja esistente e
stampa i comandi effettivi di compilazione/link. Comandi eseguiti dalla radice:

```sh
./ninja -C build -j 2 SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-final3
/private/tmp/su2-final3/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-final3-cli --cli
python3 UnitTests/run_primitive_euler_tests.py --asan --output-dir /private/tmp/su2-final3-asan
ASAN_OPTIONS=detect_leaks=0 /private/tmp/su2-final3-asan/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
leaks --atExit -- /private/tmp/su2-final3/test_driver '[primitive-factory-unwind]'
```

Build, suite normale, wrapper CLI, build ASan e suite ASan terminano con
uscita 0. La suite normale e quella ASan superano entrambe 34 casi/142380
asserzioni. I primi 29 casi rimangono presenti e superati. Log:
`/private/tmp/su2-final3-build.log`, `su2-final3-results.log`,
`su2-final3-cli.log`, `su2-final3-asan-build.log`,
`su2-final3-asan-results.log`; build CFD in `su2-phase3-final-build.log`.
I quattro log CLI dettagliati sono in `/private/tmp/su2-final3-cli/`.

ASan strumenta test, factory, solver/variabili primitive e basi CVariable e
CScalarVariable; gli archivi SU2 e MPI esistenti restano non strumentati.
Il primo tentativo con annotazioni libc++ miste si è fermato con
`container-overflow` in CConfig::TokenizeString durante confronti di stringhe.
La build nativa disabilita queste annotazioni, quella ASan le abilitava nelle
istanze dei template condivisi. Il runner Darwin uniforma l'impostazione con
`_LIBCPP_HAS_NO_ASAN`: i controlli ASan su heap/stack restano abilitati, quelli
sulla capacità logica dei contenitori libc++ sono esclusi. Il risultato
positivo si riferisce a questa configurazione parziale; una build interamente
strumentata con annotazioni uniformi rimane da eseguire.

LeakSanitizer non è disponibile: anche un eseguibile minimo compilato con
`c++ -fsanitize=address` termina con uscita 134 impostando
`ASAN_OPTIONS=detect_leaks=1`, con messaggio `detect_leaks is not supported on
this platform`. Per questo la suite ASan usa esplicitamente detect_leaks=0.
Il primo tentativo di `leaks` nel sandbox termina con uscita 255 perché non
ottiene il task port. Ripetuto fuori dal sandbox con autorizzazione, il test
Catch2 passa (1 caso/38 asserzioni), ma `leaks` termina con uscita 1 e riporta
**32 leak, 33600 byte**, con stack OpenMPI/PMIx. Il rapporto non mostra stack
del ramo primitivo; non sono stati applicati filtri né introdotti contatori.
Questo non certifica assenza di leak del processo o del nuovo ramo. Rapporto:
`/private/tmp/su2-final3-leaks-escalated.log`.

Warning non bloccanti: pragma GCC ignorati da Clang, librerie duplicate,
archivi con oggetti senza simboli, bind TCP MPI negato nel sandbox e
WINDOW_START_ITER adeguato a RESTART_ITER per la configurazione temporale.
La suite SU2 completa, MPI distribuito e una build completamente sanitizzata
non sono stati eseguiti: la verifica resta focalizzata sulla fase 3 seriale.

### Modifiche e valutazione finale

In questa verifica sono cambiati soltanto:

- `SU2_CFD/src/solvers/CSolverFactory.cpp`: ownership RAII del nuovo solver
  fino alla registrazione dei metadati.
- `UnitTests/SU2_CFD/primitive_euler_solver_tests.cpp`: ampliamento dei test
  di dispatch e cinque nuovi casi di robustezza.
- `UnitTests/run_primitive_euler_tests.py`: runner riproducibile normale,
  ASan parziale e smoke test CLI.
- `Docs/primitive_euler_phase3.md`: questo rapporto.

L'unica correzione alla produzione riguarda il rischio di perdita su
eccezione nella registrazione dei metadati, identificato staticamente;
nessun leak primitivo è stato riprodotto dinamicamente. Le altre modifiche
rafforzano le verifiche. La fase 3 costituisce una base funzionale pronta
per iniziare la fase 4, con i limiti di verifica della memoria esplicitati
sopra. Non risultano fallimenti funzionali residui nella selezione pertinente.
Non è stata iniziata la fase 4 e non sono stati creati commit o eseguiti push.
