# Archiviazione dell’esperimento quadratico e solver principale affine

## Stato iniziale e salvataggio

Branch iniziale `feature/primitive-path-solver`, HEAD
`04e532c513` (fine FASE 3), coincidente con origin al momento dell'operazione.
Nessuna modifica staged; 14 file tracciati modificati e 14 file nuovi.
La FASE 4 non era committata. Nessuno stash è stato usato.

Prima di modificare file sono stati salvati i 28 file in
`/private/tmp/primitive-quadratic-archive/working-tree-before-separation.tar.gz`,
insieme a manifest.json, unstaged.patch e staged.patch. L'archivio è una
precauzione temporanea; il salvataggio persistente è il commit Git seguente.

`experiment/primitive-quadratic-boundary` punta a
`d1732e90e3d2bd18499f4b585bd3c6d6dd4b95d7`:
“Archive phase 4 spatial operator and quadratic boundary experiment”.
Il commit contiene integralmente lo stato iniziale dei 28 file, inclusi
FASE 4, audit, implementazione quadratica, test e artefatti. Si verifica
l'identità di ogni file con l'archivio precedente alla separazione.

Il branch principale riparte dal commit di FASE 3, riceve FASE 4 e diagnostiche
affini in un proprio commit, senza avere il commit sperimentale tra gli antenati.
Non è stato usato cherry-pick del commit quadratico né un reset distruttivo.

## Separazione delle modifiche

FASE 4 e robustezza conservate sul principale:

- CPrimitiveEulerSolver.hpp/.cpp: operatore affine completo, gradienti/limiter,
  positività scalata, callback, sorgente nodale, residuale transazionale,
  guardia contro il rientro e collegamento Upwind_Residual.
- primitive_euler_reconstruction.hpp: Reconstruct affine, VolumeProduct,
  ReconstructionPositivityFactor robusto. Nessuna Hessiana o traccia quadratica.
- primitive_euler_geometry.hpp: mezze facce fisiche e validazioni geometriche.
- CDriver.cpp, CIterationFactory.cpp, CConfig.cpp e config_template.cfg:
  messaggi FASE 4 e blocco esplicito dell'evoluzione; nessuna opzione quadratica.
- Fixture, test del solver e primitive_euler_residual_tests.cpp: tutta la
  copertura precedente, incluso il difetto di consistenza MUSCL NO.
- UnitTests/meson.build e runner: residuale e diagnostiche affini registrati.
- Documenti FASE 3/4 aggiornati e CSV MMS precedente conservato.

Audit numerico mantenuto:

- primitive_euler_accuracy_tests.cpp conserva i tre test indipendenti MMS
  regionale, LS/WLS su campi quadratici e momenti delle facce AGGREGATED.
- Le derivate analitiche sostituite nel riassemblaggio test-only isolano le
  cause d'errore; nessuna derivata analitica viene fornita al solver reale.
  Il tipo Hessian è locale al test, senza API o cache di produzione.
- Sono rimossi da questo branch il fit test-only a due anelli, il suo candidato,
  i test delle API quadratiche e le asserzioni d'ordine della produzione quadratica.
  Quei test appartengono interamente al branch sperimentale.
- Gli altri test non hanno soglie indebolite o risoluzioni ridotte.
- Export CSV e grafico ora mostrano soltanto l'operatore affine; report aggiornati
  senza attribuire al principale l'ordine >=2 ottenuto nell'esperimento.

Esclusivi dell'esperimento:

- SU2_CFD/include/numerics/flow/primitive_euler_boundary_fit.hpp;
- SU2_CFD/src/numerics/flow/primitive_euler_boundary_fit.cpp;
- UnitTests/SU2_CFD/primitive_euler_quadratic_tests.cpp;
- Docs/primitive_euler_phase4_quadratic.md;
- Hessiane/cache/tracce fisiche quadratiche nel solver, registrazione
  PRIMITIVE_BOUNDARY_RECONSTRUCTION e relativo getter CConfig;
- blocchi FitQuadratic/candidato/MMS quadratica in accuracy_tests, registrazioni
  Meson/runner e versioni complete di CSV/SVG e rapporto d'ordine.

Il commit sperimentale contiene anche tutti i file di base e non dipende dal
commit affine successivo. Per elencarlo o recuperare un file senza cambiare branch:

```sh
git show --stat experiment/primitive-quadratic-boundary
git show experiment/primitive-quadratic-boundary:Docs/primitive_euler_phase4_quadratic.md
```

## Limiti numerici conservati

Il principale non soddisfa il criterio globale >=2: MMS affine regolare circa
1.5 e distorta circa 1. Il difetto del contatto MUSCL NO e MEAN_COEFF rimangono
aperti. Questa operazione non introduce correzioni numeriche, FASE 5 o integratori.
Il kernel path-conservative, i segni delle fluttuazioni e il solver conservativo
non vengono modificati dalla separazione.

## Comandi di verifica

```sh
./ninja -C build -j 2 SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-affine-archive-tests
/private/tmp/su2-affine-archive-tests/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-affine-archive-cli --cli
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/su2-order-mpl XDG_CACHE_HOME=/private/tmp/su2-order-font-cache python3 UnitTests/export_primitive_euler_order.py --log /private/tmp/su2-affine-archive-results.log --plot
```

Gli esiti finali sono riportati sotto.
Nessun push, reset --hard, clean -fd o riscrittura di cronologia.

## Inventario dei file principali conservati

- `Common/include/CConfig.hpp`
- `Common/src/CConfig.cpp`
- `Docs/primitive_euler_phase3.md`
- `Docs/primitive_euler_phase4.md`
- `Docs/primitive_euler_phase4_order.svg`
- `Docs/primitive_euler_phase4_order_audit.md`
- `Docs/primitive_euler_phase4_order_moments.csv`
- `Docs/primitive_euler_phase4_order_regions.csv`
- `Docs/primitive_euler_phase4_spatial_mms.csv`
- `SU2_CFD/include/geometry/primitive_euler_geometry.hpp`
- `SU2_CFD/include/numerics/flow/primitive_euler_reconstruction.hpp`
- `SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp`
- `SU2_CFD/src/drivers/CDriver.cpp`
- `SU2_CFD/src/iteration/CIterationFactory.cpp`
- `SU2_CFD/src/meson.build`
- `SU2_CFD/src/solvers/CPrimitiveEulerSolver.cpp`
- `UnitTests/SU2_CFD/primitive_euler_accuracy_tests.cpp`
- `UnitTests/SU2_CFD/primitive_euler_mesh_fixture.hpp`
- `UnitTests/SU2_CFD/primitive_euler_residual_tests.cpp`
- `UnitTests/SU2_CFD/primitive_euler_solver_tests.cpp`
- `UnitTests/export_primitive_euler_order.py`
- `UnitTests/meson.build`
- `UnitTests/run_primitive_euler_tests.py`
- `config_template.cfg`

A questi si aggiunge questo rapporto. Il diff del commit principale rispetto
al padre elenca le modifiche effettive; alcuni file dell’esperimento, come
CConfig.hpp e SU2_CFD/src/meson.build, coincidono di nuovo con la FASE 3.

## Esiti finali

- Build nativa: uscita 0. La modifica CConfig ha richiesto la ricompilazione
  delle dipendenze; warning non bloccanti di Clang/librerie duplicate e il
  precedente avviso Meson sul git config del sottomodulo nel sandbox.
- Suite affine completa: **58 casi, 1.380.527 asserzioni superate**, uscita 0.
  Include i 55 casi precedenti alle diagnostiche d'ordine e i tre test
  diagnostici indipendenti conservati. Nessuna soglia di quei test indebolita.
- Quattro CLI: uscite attese 1/0/1/1, tutte verificate. Primitivo bloccato,
  Euler conservativo riuscito, harmonic balance e specie rifiutati.
- Driver ricollegato agli archivi finali dopo la build; nessuna esecuzione
  valida si basa sugli archivi precedenti alla rimozione del campo CConfig.
- CSV e figura rigenerati: 2280 righe regionali, 20 righe dei momenti.
  Nessun record candidate/production quadratico nei dati del principale.
- git diff --check e compilazione Python dei due script riusciti.
- Ricerca nei sorgenti/configurazione: nessuna API/opzione/cache quadratica
  del primitivo; Hessiane eventualmente presenti in altri modelli SU2 sono
  estranee a questo esperimento e non vengono rimosse.
- Verifica byte per byte: tutti i 28 file del salvataggio iniziale corrispondono
  al commit sperimentale. Kernel, geometria e test residuali non sono cambiati
  durante la separazione.

Log temporanei `/private/tmp/su2-affine-archive-{build,tests-build,tests-relink,results,cli}.log`.
Nessun nuovo run ASan in questa operazione: il rapporto del branch sperimentale
conserva la verifica ASan parziale precedente. Sono stati rieseguiti build,
Catch2 e CLI richiesti; non l'intera suite globale SU2 né esecuzioni MPI.

Il commit principale contiene esclusivamente FASE 4 affine e audit utile;
il suo hash è indicato nel rapporto finale della sessione e si ottiene con
`git rev-parse feature/primitive-path-solver`. HEAD finale è su quel branch,
mentre il riferimento sperimentale resta fermo a d1732e90e3. Entrambi derivano
dalla FASE 3 e il commit sperimentale non è antenato del principale.
L'esperimento è integralmente recuperabile e il working tree finale è pulito.
Nessuna FASE 5 e nessun push.

Il commit di implementazione affine è `feb0b0e49d`. Un controllo staged ha
evidenziato spazi finali emessi da Matplotlib nell’SVG e una riga vuota finale
nel test; sono normalizzati in un commit locale successivo, senza cambiamenti
alle formule o ai test. L’esportatore applica la stessa normalizzazione.
