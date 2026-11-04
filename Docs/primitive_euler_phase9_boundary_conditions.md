# FASE 9 — Condizioni al contorno primitive

## 9.1 Audit preliminare (prima delle modifiche produttive)

Branch `feature/primitive-path-solver`, HEAD `0967645dba` (FASE 8), working
tree inizialmente pulito. Consultato il rapporto FASE 8 e i percorsi locali
di implementazione. I riferimenti seguenti si riferiscono a questo HEAD.

Il dispatch effettivo è `CIntegration::Space_Integration`,
`SU2_CFD/src/integration/CIntegration.cpp:85–208`: interfacce fluide;
BC deboli; valori imposti; CUSTOM/BC termiche forti; accumulo periodico;
infine simmetria e parete. Non tutte le stringhe MARKER del parser sono
BC Euler attive: marker geometrici, monitoraggio, carichi e deformazione
non prescrivono un problema fluidodinamico.

| Marker/modalità | Euler FVM nativo: codice e dati | Fisica/dipendenze | Primitivo iniziale / intervento |
|---|---|---|---|
| MARKER_EULER | CFVMFlowSolverBase.hpp:1130–1134 chiama BC_Sym_Plane | Slip, vincolo nodale e proiezione residuale; non solo ghost | Assente; studiare operatore vincolato |
| MARKER_SYM | CFVMFlowSolverBase.inl:1170–1334 | Stessa routine per Euler; normali anche corrette alle intersezioni, gradienti corretti nel percorso nativo | Assente; condividere formulazione con parete |
| MARKER_FAR | CEulerSolver.cpp:4822–5077, dati freestream, Mach/AoA, gas | Riemann R+/R−, entropia e tangente; stato primitivo poi flusso conservativo | Solo stato esterno FASE 8; distinguere vero farfield |
| MARKER_INLET, TOTAL_CONDITIONS | CEulerSolver.cpp:6954–7242; INLET_TYPE; (tag,Tt,Pt,dx,dy,dz), INLET_USE_NORMAL | Tre dati subsonici e R+ uscente; equazione quadratica, clipping nativo a Mach1 | Assente; costruzione primitiva con diagnostica, senza copiare clipping |
| MARKER_INLET, MASS_FLOW | Stessa routine, 7145–7211; (tag,rho,modulo velocità,direzione) | Il nome non significa portata integrata: densità e velocità prescritte, pressione da R+ | Assente; adattamento locale |
| MARKER_OUTLET | CEulerSolver.cpp:7244–7417, (tag,p statica) | Entropia/tangente/R+ interni; pressione subsonica. Il codice nativo classifica usando Mach totale, non normale | Assente; distinguere classificazione normale e backflow |
| MARKER_SUPERSONIC_INLET | CEulerSolver.cpp:7419–7506; (tag,T,p,u,v,w) | Dati statici completi, rho=p/(RT), flusso upwind | Assente; validare tutte le caratteristiche entranti |
| MARKER_SUPERSONIC_OUTLET | CEulerSolver.cpp:7508–7628, lista tag | Stato interno extrapolato, nessun dato entrante | Assente; validare regime normale |
| MARKER_CUSTOM | CFVMFlowSolverBase.inl:1518–1582 | Con VerificationSolution impone U ai nodi, azzera residuale; altrimenti errore | Completo per callback esterna debole FASE 8; preservare, non dichiarare identità col CUSTOM conservativo |
| MARKER_RIEMANN | CEulerSolver.cpp:5079–5508; tipo, due scalari, vettore | TOTAL_CONDITIONS_PT, STATIC_SUPERSONIC_INFLOW_PT/PD, DENSITY_VELOCITY, STATIC_PRESSURE; proiezione su autovettori conservativi e EOS | Assente; richiede scelta esplicita dell'operatore caratteristico primitivo |
| MARKER_GILES / TurboRiemann | CEulerSolver.cpp:5510–6952 | Medie di sezione, Fourier, rilassamenti, geometria turbo; anche opzioni radiali/mixing | Assente; nuova infrastruttura, rinviato |
| MARKER_ENGINE_INFLOW / EXHAUST | CEulerSolver.cpp:7630–8111 | Feedback Mach/portata/pressione, medie e dati totali; non semplici inlet/outlet statici | Assente; controllo di sistema da rinviare |
| MARKER_ACTDISK | CEulerSolver.cpp:8113–8850 circa | Coppia inlet/outlet, salti, carichi variabili/BEM, donatori | Assente; interfaccia accoppiata, rinviata |
| MARKER_PERIODIC | CConfig.cpp:1598–1605; CFVMFlowSolverBase.inl:1337–1351 | Coppie, centro/angoli/traslazione; comunicazioni PERIODIC_RESIDUAL, volumi e gradienti accoppiati | Assente; CAggregatedGeometry rifiuta periodicità, nuova infrastruttura necessaria |
| MARKER_FLUID_INTERFACE / ZONE_INTERFACE | CFVMFlowSolverBase.inl:1354–1515 e CDriver | Stato donatore interpolato, trasferimenti multizona/sliding | Assente; fuori perimetro seriale singola zona |
| MARKER_NEARFIELD / INTERNAL | CConfig.cpp:1578–1597; nessun ramo BC nel dispatch sopra | Non costituiscono da soli una condizione esterna generica; richiedono verificare accoppiamento/uso geometrico del caso | Nessun supporto fisico; non accettarli come extrapolation |
| ISOTHERMAL / HEATFLUX / HEAT_TRANSFER / CHT / SMOLUCHOWSKI_MAXWELL | Dispatch CIntegration.cpp:161–187 | Fisica viscosa/termica; CConfig.cpp:3701 esclude pareti viscose per Euler | Non applicabili a questo Euler inviscido |
| SPECIES / INLET_TURBULENT / WALL_FUNCTIONS | Parsing aggiuntivo CConfig.cpp:1607–1670 | Dati di modelli accoppiati, non ulteriori BC Euler ideali | Non applicabili |
| PYTHON_CUSTOM / DEFORM / MONITORING / PLOTTING / FLUID_LOAD | Parsing CConfig e metadati marker | Interfacce o selezioni accessorie, non una nuova legge Euler | Fuori perimetro; accettazione sintattica non equivale a supporto |

### Convenzioni e parti riutilizzabili

Il parser nativo è in `Common/src/CConfig.cpp:1572–1798`. Le normali
vertex native sono negate per ottenere l'uscente in CEulerSolver; le
mezze facce primitive hanno già normale uscente integrata. La callback
primitiva riceve la normale unitaria; il kernel riceve quella integrata.
Gli stati native V hanno layout T,u,v,p,rho,H, diverso da W=(rho,u,v,p).
Riutilizzabili: dati fisici/getter del parser, gas ideale, significato
di invarianti, direzioni e unità. Non riutilizzabili direttamente:
flussi, residuali e Jacobiani conservativi, proiezioni in U, codici di
comunicazione periodica che assumono quel layout.

Farfield nativo sceglie le informazioni con il segno del freestream
normale (4822–5005), non con quello locale ricostruito. Il commento storico
sui segni entranti è ambiguo: con normale uscente sono lambda<0 a entrare,
come conferma BC_Riemann:5310–5330. Le differenze di classificazione al
cambio di regime vanno dichiarate, non nascoste dietro lo stesso marker.

### 9.2 Architettura prevista e rischi, prima dell'implementazione

Separare funzioni matematiche di costruzione degli stati da lettura dei
dati CConfig e da adapter del driver; riusare ComputeFluctuations per
D− sulle mezze facce. Nessun secondo operatore interno. CUSTOM resta
compatibile con FASE 8. FAR fisico usa dati freestream nativi; una vecchia
configurazione FAR a stato prescritto deve migrare esplicitamente a CUSTOM,
senza cambio silenzioso di problema fisico.

La sola riflessione non garantisce il vincolo nodale u.n=0: se W_ext=W_int
per uno stato tangente, D−=0, mentre un gradiente normale di pressione
produce comunque RHS_un=−p_n/rho. La routine nativa proietta anche il
residuale totale. Occorre quindi distinguere una parete debole riflessa
da una parete nodale vincolata, senza proclamare equivalenza. Progettare
e verificare separatamente l'eventuale vincolo, inclusi vertici con due
normali, sorgenti e atomicità FE. Non modificare il kernel per nasconderlo.

Rischi ulteriori: dati totali incompatibili con R+; backflow sottoprescritto
con la sola pressione di uscita; numero di caratteristiche che cambia al
sonico; stato intermedio non positivo; riflessione numerica e inaccuratezza
dei gradienti di bordo. Nessun clipping automatico degli stati costruiti.
Periodicità, turbo, actuator disk e interfacce restano problemi distinti.

Ordine operativo: proprietà locali e controesempi; parete/simmetria;
farfield; inlet/outlet; primo passo; benchmark lisci su regular/jitter;
confronto vero CEulerSolver; regressioni FASI 1–8; report con supporto
completo distinto da prototipi e casi rinviati.

## 9.3 Implementazione e stato del supporto

Sono operativi nel driver ordinario parete/simmetria **rettilinee**,
inlet TOTAL_CONDITIONS e MASS_FLOW subsonici, inlet statico supersonico,
outlet di pressione subsonico e outlet supersonico. CUSTOM FASE 8 è
preservato. Il farfield è disponibile solo con scelta esplicita
`PRIMITIVE_BOUNDARY_MODE=LOCAL_CHARACTERISTIC`: è un modello locale
caratteristico verificato nel regime descritto sotto, **non una sostituzione
certificata della chiusura nativa a ampiezza finita**.

**Non si dichiara completata l'intera FASE 9.** Restano pareti curve,
backflow sottoprescritto, transizioni soniche finite, RIEMANN e le BC
accoppiate. Anche gli ordini dei problemi evoluti non sono universalmente
due. Questi limiti non sono compensati modificando il kernel.

File produttivi, elenco completo:

* `SU2_CFD/include/numerics/flow/primitive_euler_boundary.hpp`: funzioni
  pure degli stati, autovalori, validazione di normali/dati/regimi;
* `SU2_CFD/include/iteration/primitive_euler_boundary_config.hpp`: adapter
  dei getter CConfig, marker, unità e dati fisici nativi;
* `SU2_CFD/src/iteration/CPrimitiveEulerIteration.cpp`: dispatch della
  callback, FREESTREAM iniziale, controlli espliciti e opt-in FAR;
* `SU2_CFD/src/solvers/CPrimitiveEulerSolver.cpp`: ammissione dei marker,
  riflessione e **sola reazione nodale di slip** nel contributo di bordo.

Il kernel della FASE 1, matrici e valore assoluto, cammino, segni D±,
geometria, quadrature, variabili NODE, funzioni di ricostruzione/limiter,
VolumeProduct e algoritmo AdvanceExplicitEuler non sono modificati.
Il nuovo vincolo di parete è un operatore di bordo esplicito: non si
presenta come identico al precedente bordo CUSTOM/FAR prescritto.

Test/infrastruttura: `primitive_euler_boundary_tests.cpp`,
`primitive_euler_boundary_audit.cpp`, `primitive_euler_boundary_cases.py`,
`check_primitive_euler_boundary_cli.py`, `run_primitive_euler_boundary_audit.py`,
`export_primitive_euler_boundary_audit.py`; registrazione in UnitTests/Meson
e runner. La fixture ammette marker misti dichiarati esplicitamente e log
non silenziati, mantenendo i default precedenti. Il runner FASE 8 aggiunge
solo `--skip-negative`: gli assert numerici restano identici; i vecchi
test che rifiutavano wall/sym/FAR prescritto sono sostituiti dal contratto
negativo FASE 9, senza riscrivere i CSV delle fasi precedenti.
Il nome/commento del caso in `primitive_euler_residual_tests.cpp` ora
precisa che SYM rifiuta lo stato iniziale non tangente, non il marker:
gli assert e i dati precedenti sono conservati.

### Configurazione e unità

Nessuna nuova opzione CConfig o dipendenza. La stringa esistente
PRIMITIVE_INITIALIZATION accetta anche FREESTREAM: rho e p, oppure p e T,
con EOS rho=p/(RT); velocità da Mach*sqrt(gamma*p/rho) e AOA in gradi.
È lo stesso calcolo Euler inviscido di CEulerSolver.cpp:901–936.
INIT_OPTION=REYNOLDS non ricava la densità da viscosità in Euler inviscido,
come conferma il ramo `if (viscous)` a 936–946. UNIFORM e ISENTROPIC_VORTEX
mantengono il contratto FASE 8; il vortice resta su CUSTOM analitico.

FAR usa esclusivamente i dati freestream nativi e richiede opt-in
LOCAL_CHARACTERISTIC. Un vecchio cfg FAR con modalità UNIFORM/VORTEX o
PRIMITIVE_BOUNDARY_STATE esplicita fallisce: non cambia problema senza
avviso. Il cfg con FAR non può ancora assegnare questi dati globali a un
altro CUSTOM; i casi misti analitici sono supportati dalla callback C++.
Non si dichiara supporto generale ai profili inlet da file.
L'API spaziale a basso livello richiede la callback fisica quando FAR è
LOCAL_CHARACTERISTIC: non ripiega sul vecchio stato esterno prescritto.
Il fallback storico resta solo per le configurazioni API legacy, così da
preservarne il contratto; il driver passa sempre l'adapter configurato.

MARKER_INLET=(tag,Tt,Pt,dx,dy,dz) oppure, con INLET_TYPE=MASS_FLOW,
(tag,rho,modulo_velocità,dx,dy,dz). MASS_FLOW non è una portata integrata.
La direzione viene normalizzata; INLET_USE_NORMAL usa −n locale.
MARKER_SUPERSONIC_INLET=(tag,T,p,u,v,w): velocità dimensionali, non direzione.
La componente z deve essere zero. MARKER_OUTLET=(tag,p_statico).
SYSTEM_MEASUREMENTS=SI, REF_DIMENSIONALIZATION=DIMENSIONAL e riferimenti unitari restano richiesti;
R e gamma sono quelli del gas ideale configurato, senza conversioni
conservative nascoste. Un dato non ammissibile viene rifiutato, non troncato.
Il percorso US è rifiutato: CEulerSolver.cpp:852–861 sostituisce anche R
con 1716.49 per STANDARD_AIR, mentre il postprocessing comune imposta
287.058. Non si riusa quel valore SI per interpretare dati US.

## 9.4 Formulazione matematica e collegamento al residuale

Sia n uscente unitaria; la normale integrata della mezza faccia è S=ell*n.
La callback costruisce B(W,n,dati). Si usa una sola volta il kernel

    r_b = D−(W_int, B, S),       rhs = −(r_faces+r_volume+r_boundary)/V + source.

Non si copia un flusso CEulerSolver. Invertire la connettività locale non
inverte la normale fisica, determinata dalla geometria nativa. Invertire
la normale fisica di un farfield scambia invece ingresso/uscita: non è
un'invarianza da richiedere. Per slip, riflessione e proiettore sono
invarianti sotto n→−n.

### Parete e simmetria: riflessione e reazione

Lo stato riflesso è (rho,u−2(u.n)n,p), positivo in rho/p per costruzione.
Il cammino medio ha velocità normale nulla e, per uno stato tangente,
il salto e D− sono nulli. Questo **non** rende nullo il RHS normale totale:
la pressione può ancora produrre −p_n/rho nel volume.

Perciò, dopo l'assemblaggio, si costruisce una base ortonormale N delle
normali di slip incidenti al nodo e si pone

    rhs_vel_new = (I − N Nᵀ) rhs_vel_raw,
    boundary_vel_new = boundary_vel_raw + V N Nᵀ rhs_vel_raw.

È la reazione del vincolo u.n=0; non è un termine volumetrico artificiale.
Facce interne, volume, sorgente riportata e due equazioni scalari restano
invariati. L'identità rhs=−(faces+volume+boundary)/V+source rimane verificata.
FE usa questo unico RHS, senza proiezioni post-step o clipping. Si richiede
uno stato iniziale già compatibile, entro 256 epsilon scalati sulla velocità;
altrimenti errore esplicito prima della modifica di W/residuale/tempo.
La selezione del rango usa 256 epsilon, non un coefficiente di stabilizzazione.

Due pareti piane indipendenti che si incontrano impongono entrambi i vincoli;
ai vertici di spigolo la velocità deve quindi essere compatibile con entrambi.
Normali non parallele dello **stesso marker** sono rifiutate localmente:
trattarle entrambe come vincoli su una curva imporrebbe erroneamente no-slip.
Separare gli spigoli reali in marker distinti è ammesso; suddividere una
curva in tanti marker non realizza una parete curva. Non sono stati modificati
o rispecchiati i gradienti: la loro accuratezza al bordo resta misurata,
non corretta empiricamente. Euler e SYM condividono esattamente l'operatore.

### Farfield locale: caratteristiche entranti, non Riemann finito nativo

Con c²=gamma*p/rho e un=u.n, gli autovalori sono un−c,un,un,un+c.
Per delta=W_inf−W_int si usano gli autovettori primitivi di A_n(W_int):

    a± = (delta_p ± rho*c*delta_un)/(2*c²),
    a0 = delta_rho − delta_p/c²,       at = delta_ut,
    r± = (1, ±c*nx/rho, ±c*ny/rho, c²),
    r0 = (1,0,0,0),                  rt = (0,−ny,nx,0).

B=W_int+somma(a_j*r_j per lambda_j<0). Sono selezionati 4,3,1,0 dati
rispettivamente in ingresso supersonico, ingresso subsonico, uscita
subsonica, uscita supersonica. Le caratteristiche con velocità nulla
non prescrivono dati; il confronto dei segni usa 64 epsilon sulla scala
di un,c per non classificare roundoff di una rotazione come flusso entrante.
Al tangenziale possono cambiare i segni di un per piccole onde acustiche:
la classificazione viene aggiornata allo stadio successivo.

Questa è una chiusura **linearizzata nello stato interno**, pur usata con
il kernel path-conservative non lineare originale. A ampiezza finita
non conserva esattamente gli invarianti non lineari del native FAR.
Positività di B è obbligatoria; salti che attraversano un regime sonico
non supportato sono rifiutati. Nessun blending o entropy fix è introdotto.

### Inlet e outlet

Per l'outlet subsonico si mantengono K=p/rho^gamma, componente tangenziale
e R+=un+2c/(gamma−1) interni; p_b=p_out,
rho_b=(p_out/K)^(1/gamma), c_b=sqrt(gamma*p_out/rho_b),
un_b=R+−2c_b/(gamma−1). Devono restare 0<un_b<c_b.
Per uscita supersonica non si impongono dati. La classificazione è con
Mach **normale**, mentre il native Outlet usa Mach totale: differenza
esplicita sui casi obliqui a Mach totale supersonico e normale subsonico.
Backflow/tangenza di un outlet di sola pressione sono sottoprescritti
(tre informazioni entranti invece di una) e vengono rifiutati.

Per l'inlet totale si impongono Tt,Pt e direzione d; alfa=d.n<0,
q=|u|, c0²=gamma*R*Tt e si risolve

    R+_int = alfa*q + 2*sqrt(c0²−(gamma−1)*q²/2)/(gamma−1).

La radice positiva dell'equazione quadratica viene verificata nella
relazione non quadrata e nel regime normale subsonico. Si ricavano
T_b=c_b²/(gamma*R), p_b=Pt*(T_b/Tt)^(gamma/(gamma−1)), rho_b=p_b/(R*T_b).
Nessun clipping nativo a Mach1 viene copiato: dati incompatibili falliscono.
MASS_FLOW prescrive rho,q,d; R+_int determina c_b e quindi p_b.
Inlet statico supersonico richiede tutte e quattro le caratteristiche
entranti sia nello stato interno sia in quello prescritto. Outlet
supersonico richiede nessuna caratteristica entrante.

## 9.5 Verifiche locali, geometriche ed end-to-end

Otto test Catch dedicati controllano riflessione involutiva, rho/p e
tangente, normali ruotate/invertite, numero di caratteristiche, proiezioni
acustiche entranti/uscenti, modi entropici/vorticosi, sonico compatibile,
R+ e dati totali, pressione/entropia/tangente all'outlet, dati invalidi.
Le comparazioni Approx nuove fissano esplicitamente epsilon=64 epsilon
double: non si usa la tolleranza relativa predefinita di Catch.

Sulle mesh reali: regular/jitter, rettangolo asimmetrico, rotazione .37 rad,
connettività invertita, mezze facce di diversa lunghezza, angoli misti
wall/CUSTOM e intersezioni di slip. Il canale affine accelera tangenzialmente
di −.1 e mantiene un nullo; la pressione non costante nel test chiuso
dimostra la necessità della reazione. Il centro interno mantiene il
proprio RHS, senza proiezioni. La curva non supportata fallisce esplicitamente.

Il runner CLI finale esegue **41 casi positivi**: 10 configurazioni per
due famiglie e due solver, più uno smoke Navier–Stokes. Comprende
TEMPERATURE_FS e DENSITY_FS, dati totali e MASS_FLOW. Il primitivo preserva
esattamente lo stato uniforme; massimo scarto nativo 2.22045e-16.
Otto negativi verificano FAR legacy, FAR senza opt-in, unità US, pressione negativa,
direzione inlet uscente, outlet supersonico con flusso subsonico,
backflow e periodicità. I 16 esempi distribuiti N=16 sono stati eseguiti
direttamente con SU2_CFD; il generatore ne riproduce esattamente i byte.

## 9.6 Campagna scientifica non uniforme

Si usano **i veri CPrimitiveEulerSolver e CEulerSolver**; per il secondo
si preparano comunicazioni/normalizzazione e profili con SetUniformInlet,
si chiamano Preprocessing, Upwind_Residual, le BC native e
ExplicitEuler_Iteration. CUpwRoe_Flow è soltanto il numerics di bordo nativo,
non un nuovo solver conservativo. Le condizioni iniziali fisiche, mesh,
gamma, unità, dt e tempo finale coincidono. Fallback nativi su nodi o
ricostruzioni non fisiche sono vietati dal test; ogni candidato nativo è
validato prima di un successivo preprocessing che potrebbe ripristinarlo.

Le regioni sono all, boundary(layer0), ring1, ring2, deep3(layer>=3).
La distanza logica è anche distanza dal bordo nel grafo di queste
triangolazioni, invariata dal jitter. Le norme per componente sono

    L1_V = sum(V_i*|e_i|)/sum(V_i),
    L2_V = sqrt(sum(V_i*e_i²)/sum(V_i)),       Linf = max |e_i|.

I CSV regionali compatti sommano L1 delle componenti, usano L2 euclidea
del vettore e Linf massimo componente: non sono una norma di energia fisica.
E_P, E_C confrontano le primitive con il riferimento; E_PC confronta i
due stati. Non si deduce equivalenza a mesh finita.

### Contatto esatto e difetto semidiscreto

Dominio [-5,5]^2, gamma=1.4, rho=1+.02*(1−r²)^4 per r²<1, altrimenti1;
r²=((x−3.5−u0*t)²+y²)/4, u=u0, v=0, p=1. Il profilo è C3, non affine,
soluzione esatta non lineare senza sorgente. T=.5, u0=.4 oppure2 per il
caso supersonico. Un caso ruota geometria, velocità e analitica di .37 rad.
Sette combinazioni BC, due famiglie, N=8,16,32,64: **56 coppie evolute**.
Nt=2N², dt=T/Nt. A t=0 si misura anche tau_P=rhs_P−W_t analitico,
nelle cinque regioni. Nessuna derivata numerica entra nella sorgente
di riferimento e nessuna sorgente è aggiunta al solver.

Per TOTAL/MASS + pressure outlet, il contatto è compatibile con entrambe
le chiusure: a N=64, errori L2 della densità e ordine 32→64:

| Famiglia/caso | E_P | p_P | E_C | p_C | E_PC | p_PC |
|---|---:|---:|---:|---:|---:|---:|
| regular TOTAL/MASS | 7.75196e-6 | 2.30924 | 7.71882e-6 | 2.30876 | 1.71969e-7 | 2.10359 |
| jitter TOTAL/MASS | 1.24856e-5 | 1.99645 | 1.24675e-5 | 1.99481 | 1.83543e-7 | 2.04050 |
| regular supersonic | 3.09644e-5 | 2.20412 | 3.07792e-5 | 2.20471 | 6.44143e-7 | 1.54632 |
| jitter supersonic | 3.14519e-5 | 2.19065 | 3.12856e-5 | 2.18954 | 6.29643e-7 | 1.61611 |

Massimo un sulle pareti nei campioni P/C: 2.96985e-15, inclusa la
rotazione. Euler wall e symmetry danno risultati identici; la rotazione
produce differenze di roundoff. Rho/p restano positivi, anche quando
MUSCL senza limiter produce piccoli undershoot: minimi finali dell'intera
campagna contatto .9993986 e .9997100. Non si certifica una proprietà TVD.

### Risultato discriminante: FAR finito non è ancora compatibile

Per un'uscita n=(1,0), interno (1.02,.4,0,1), freestream (1,.4,0,1),
la proiezione locale lascia esattamente W interno: il salto è entropico
uscente. Il codice nativo combina R+ interno e R− infinito e produce

    rho_b=1.0456275000, un_b=.3708560391, p_b=1.0353508750.

La differenza persiste anche senza mesh: non è un errore del kernel
primitivo e non si elimina cambiando segni o raffinando dt. Nel contatto
con FAR alle uscite, il riferimento analitico non soddisfa dunque
esattamente la legge FAR nativa a ampiezza finita. A N=64 regular,
E_C,rho=1.06570e-4 ed E_PC,rho=1.03762e-4 (ordini .2664 e .06135);
su jitter E_C=1.07879e-4 ed E_PC=1.04331e-4. Il primitivo ha E_P come
la tabella precedente, ma questo **non** prova superiorità su un identico
problema di bordo: le chiusure finite sono diverse. I casi wall/sym con
FAR agli estremi ereditano questo stesso effetto; non va attribuito alla parete.

Per questo FAR resta opt-in e **parziale**, non una compatibilità nativa
certificata. Occorre revisionare separatamente la legge non lineare dei
dati entranti desiderata e confrontare entrambi i solver con un riferimento
che la soddisfi. Riprodurre automaticamente il mixing nativo per annullare
E_PC nasconderebbe questa scelta matematica; non è stato fatto.

### Vorticità convettiva esatta

Tre casi aggiuntivi, N=8,16,32 su regular/jitter (**18 coppie**):
rho=p=1, u=u0, v=.01*F(x−3.5−u0*t), con
F(s)=(1−s²/2.25)^4 per |s|<1.5, altrimenti0. È una soluzione Euler
non lineare esatta, poiché tutti i campi sono indipendenti da y e
v_t+u0*v_x=0. FAR/outlet di pressione/outlet supersonico sono sui lati x;
i lati y usano la traccia analitica CUSTOM, forte nativa e debole primitiva
come in FASE 7. L'oggetto VerificationSolution fornisce solo dati, non una MMS
o una sorgente. Si misurano tau_P e gli errori evoluti per regioni.

A N=32, L2 di v, P/C/PC: regular FAR
(1.17508e-4,1.18897e-4,1.49594e-5), jitter FAR
(1.24257e-4,1.25375e-4,1.59901e-5). Gli ordini P sono 1.6823/1.7217;
per outlet supersonico 1.6644/1.6619. La differenza PC si riduce con ordini
2.12/2.06 per FAR e 2.55/2.57 per supersonico. **Non è evidenza di ordine
due globale della soluzione**: i dati registrano il regime osservato.

### Acustica: riferimento lineare dichiarato, riflessioni misurate

Piccole perturbazioni a=1e-5, stesso F compatto, rho0=p0=1,
c0=sqrt(1.4), u0=.4 (zero per la parete riflettente). Il riferimento
incident/reflected è la soluzione del **sistema linearizzato**, non una
soluzione non lineare esatta. Delta_rho=delta_p/c0²,
delta_u=±delta_p/c0; per MASS_FLOW si include anche l'entropia convettiva
necessaria a mantenere rho e u prescritti. T=6, Nt=8N²; 30 coppie su
N=8,16,32, entrambe le famiglie, più due N=64 regular (parete/outlet).
Le onde sono lontane dalla formazione di shock; la limitazione del
riferimento è controllata dimezzando a, non con sorgenti compensatrici.

Coefficienti lineari R=p_riflessa/p_incidente attesi: FAR0, parete+1,
pressure outlet−1, MASS_FLOW+1, TOTAL−(1−M)/(1+M)=−.49469938.
I pacchetti riflessi usano le rispettive velocità u0±c0, non una specchiatura
con la stessa velocità in flusso medio non nullo. R_fit è la proiezione
L2_V della pressione finale sul pacchetto riflesso atteso, divisa per a.
**Include anche smorzamento e dispersione nel trasporto**, non isola da
solo il coefficiente puntuale della BC. Il CSV riporta anche l'energia
di pressione normalizzata, senza assumere riflessione numerica zero.

| Caso regular | N | R lineare | R_fit P | R_fit C |
|---|---:|---:|---:|---:|
| FAR | 32 | 0 | .000418896 | .000417634 |
| MASS_FLOW | 32 | 1 | .85100146 | .85100073 |
| TOTAL | 32 | −.49469938 | −.44491302 | −.44491332 |
| parete | 32 | 1 | .83070548 | .83069010 |
| parete | 64 | 1 | .95917471 | .95916923 |
| pressure outlet | 32 | −1 | −.62015594 | −.62015593 |
| pressure outlet | 64 | −1 | −.82500357 | −.82500331 |

L'accordo P/C non sostituisce il riferimento indipendente: il pacchetto
dell'outlet resta significativamente distorto a N=64. Non si certifica
un'accuratezza del 5% di R per quell'outlet, né un ordine asintotico da
questi soli dati. Servono raffinamenti ulteriori e separazione fra errore
di trasporto e riflessione prima di specifiche quantitative più forti.

### Controllo temporale e di ampiezza

N=16 regular, terne dt,dt/2,dt/4 per tutti i sette contatti, cinque
acustici e tre vorticosi: **45 coppie evolute**, confronto delle norme
P/C/PC per ogni componente significativa. Criterio: variazione relativa
delle norme <5%, con segnalazione separata dei valori <=1e-12; nessuna
stima d'ordine su valori di roundoff. Tutti i **360 confronti** passano.
Massimi: contatto .1956%, vorticità .3594%, acustica 1.4395% (PC;
P/C <=1.1445%). Variazione assoluta massima di R_fit: .0024895.
Il test quantifica il controllo a N=16, non dimostra indipendenza temporale
uniforme a qualunque risoluzione; dt scala come h² nei raffinamenti.

Cinque ulteriori coppie acustiche a a/2 danno rapporti fra errori P/C
.4999964–.5000514. Il riferimento lineare è quindi adeguato a questa
misura di discretizzazione nel regime campionato; non diventa per questo
una soluzione Euler esatta. Le soglie non sono aumentate per accomodare
i risultati e non si deduce ordine spaziale dal solo dt∝h².

## 9.7 Regressioni FASI 1–8 e integrità

Suite complessiva: 78 casi, 22 254 180 asserzioni; 8 casi/943 asserzioni
nuovi. Analitica vortice: 890 controlli; MMS analitica: 3300 asserzioni;
Python output: 4 test, postprocessing: 3 test. Tutti passati.
Gli assert geometrici, kernel, MUSCL, limiter, MMS e FE precedenti restano.

Regressione driver FASE 8: otto vortici e due uniformi, stessi tempi e
soglia 2e-13, massimo driver/harness 4.4408920985e-16; uniformi invariati.
ParaView 5.13.2 legge i 30 snapshot: campi esattamente pari ai cast Float32
dei checkpoint Float64 e stessa geometria/tempi. Nessuna modifica ai writer.
Il nuovo contratto dei marker ha propri negativi: non si pretende che
la vecchia asserzione «wall non supportato» continui a descrivere FASE 9.
Nessun CSV precedente modificato. Nessun commit, push o altra operazione
Git di scrittura; branch/HEAD rimangono quelli iniziali.

## 9.8 Riproduzione e artefatti

Esempi: `TestCases/primitive_euler_boundary_examples/{regular,jitter}`,
otto cfg per cartella e una sola mesh condivisa per famiglia. README con
comandi reali, stessi marker e dati fisici per EULER/PRIMITIVE_EULER.
Per FAR il flag LOCAL_CHARACTERISTIC è intenzionalmente specifico e
segnala la limitazione matematica sopra, non equivalenza non lineare.

Comandi eseguiti dalla radice; i runner impostano OMP_NUM_THREADS=1:

```sh
./ninja -C build SU2_CFD/src/SU2_CFD -j4
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-phase9-tests-final
OMP_NUM_THREADS=1 /private/tmp/su2-phase9-tests-final/test_driver --reporter compact
python3 UnitTests/check_primitive_euler_boundary_cli.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/su2-phase9-cli-certified
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-contact --build-only
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-contact --skip-build --n 8 16 32 64
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-acoustic --n 8 16 32 --modes acoustic_far acoustic_wall acoustic_outlet acoustic_total acoustic_mass
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-acoustic --skip-build --n 64 --families regular --modes acoustic_outlet acoustic_wall
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-temporal --n 16 --families regular --modes acoustic_far acoustic_wall acoustic_outlet acoustic_total acoustic_mass --divisors 1 2 4
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-contact-temporal --n 16 --families regular --divisors 1 2 4
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-amplitude --n 16 --families regular --modes acoustic_far acoustic_wall acoustic_outlet acoustic_total acoustic_mass --amplitude 0.000005
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-vorticity --n 8 --families regular --modes vorticity_far vorticity_outlet vorticity_supersonic
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-vorticity --skip-build --n 8 16 32 --modes vorticity_far vorticity_outlet vorticity_supersonic
python3 UnitTests/run_primitive_euler_boundary_audit.py --output-dir /private/tmp/su2-phase9-vorticity-temporal --n 16 --families regular --modes vorticity_far vorticity_outlet vorticity_supersonic --divisors 1 2 4
python3 UnitTests/run_primitive_euler_driver_e2e.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/su2-phase9-regression-driver --skip-negative
/Applications/ParaView-5.13.2.app/Contents/bin/pvpython UnitTests/validate_primitive_euler_driver_output.py --input-dir /private/tmp/su2-phase9-regression-driver --report /private/tmp/su2-phase9-regression-driver/vtk.csv
python3 UnitTests/export_primitive_euler_boundary_audit.py --contact-dir /private/tmp/su2-phase9-contact --acoustic-dir /private/tmp/su2-phase9-acoustic --temporal-dir /private/tmp/su2-phase9-temporal --contact-temporal-dir /private/tmp/su2-phase9-contact-temporal --vorticity-dir /private/tmp/su2-phase9-vorticity --vorticity-temporal-dir /private/tmp/su2-phase9-vorticity-temporal --amplitude-dir /private/tmp/su2-phase9-amplitude --cli-dir /private/tmp/su2-phase9-cli-certified --regression-dir /private/tmp/su2-phase9-regression-driver --output-dir /private/tmp/su2-phase9-package
python3 UnitTests/test_euler_vortex_output.py
python3 UnitTests/test_euler_asymptotic_postprocessing.py
/private/tmp/su2-phase8-vortex-analytic
/private/tmp/su2-mms-6b1-final/euler_mms_analytic --reporter compact
git diff --check
```

Il runner costruisce il solo harness diagnostico con le librerie native
esistenti; l'esecuzione degli esempi usa il vero SU2_CFD. L'exporter
`export_primitive_euler_boundary_audit.py` raccoglie 13 CSV e un manifest
SHA256 di sorgenti/config/mesh/output grezzi/eseguibili; nessuna simulazione
è eseguita dal postprocessor. I log e i campi grezzi rimangono sotto
/private/tmp/su2-phase9-*. Il confronto opt-in N=16 verifica che il flag
aggiunto dopo la diagnosi FAR non alteri gli errori numerici della campagna.

I dati archiviati sono `primitive_euler_phase9_{contact,regions,acoustic,
acoustic_regions,reflection,temporal,temporal_reflection,amplitude,
amplitude_reflection,cli,negative,regression,vtk}.csv`. Il manifest
`primitive_euler_phase9_provenance.json` registra gli SHA256; gli archivi
sono stati confrontati con i file grezzi e tutti gli hash verificati.
`contact` include anche i 18 casi vorticosi: 888 righe finali globali;
`regions` contiene 740 misure vettoriali semidiscrete/evolute. Sono riportati
gli ordini osservati per componente; celle vuote non significano ordine zero.

Compilatore Apple clang 15.0.0, build locale esistente; nessuna nuova
dipendenza. MPI seriale e ParaView hanno richiesto socket locali fuori
dal sandbox. La build finale riesce; i warning di pragma GCC ignorati,
archivi senza simboli e librerie MPI duplicate sono quelli già noti.
I tempi per caso sono nei CSV grezzi/colonne seconds: fino a circa 113 s
per un contatto N=64 e 345 s per il pacchetto acustico N=64 su questa macchina.
Non è stata eseguita l'intera suite multiphysics SU2.

## 9.9 Decisioni aperte e prossimi passi

| Stato | BC/attività | Decisione necessaria |
|---|---|---|
| Operativa, verificata nel perimetro | Slip/SYM rettilinee; inlet TOTAL/MASS; inlet/outlet supersonici; pressure outlet senza backflow | Mantenere limiti di dati/regime espliciti e proseguire accuratezza dei pacchetti |
| Parziale, opt-in | FAR locale | Revisionare la chiusura non lineare e il benchmark entropico prima di dichiarare sostituibilità nativa |
| Parziale | Outlet vicino al sonico o con backflow | Definire i dati convettivi entranti mancanti; non extrapolarli arbitrariamente |
| Rinviata | Parete curva | Normale nodale/ricostruzione di bordo coerente con mezze facce, senza falso no-slip |
| Rinviata | MARKER_RIEMANN | Mappatura completa dei tipi/direzioni locali e scelta della proiezione in W: quella in U non è intercambiabile a ampiezza finita |
| Nuova infrastruttura | Periodicità, interfacce, Giles/turbo, actuator disk/engine | Accoppiamento dei DOF, volumi, gradienti, trasformazioni e controlli; non ghost prescritto |
| Fuori perimetro | BC viscose/termiche, MPI, 3D, implicito, mesh mobili | Non abilitate |

La ricerca sulle condizioni di salto/shock rimane separata. Non sono
introdotti aggiustamenti conservativi, nuove ricostruzioni o integratori.
La parte FAR a ampiezza finita viene fermata qui per revisione matematica;
gli altri risultati indipendenti sono conservati e riproducibili.
