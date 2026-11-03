# FASE 6B.1 — MMS comune e difetto semidiscreto

## Risultato e limiti

Sul branch `feature/primitive-path-solver`, base
`5be7a6e72a6c082cce192e73e8dd53eb17f1b13b`, sono stati confrontati i due
operatori produttivi su **36 stati analitici**: tre famiglie di mesh, quattro
risoluzioni e tre tempi. Nessun avanzamento temporale è stato eseguito.

Le sorgenti risultano coerenti nelle due rappresentazioni; geometria e stato
fisico coincidono entro arrotondamento. Le norme trasformate mediante J
permettono un confronto nella stessa rappresentazione conservativa, ma non
rendono identici gli operatori a mesh finita.

Per la densità, sull'ultima transizione 32→64 a t=0.1:

- regular, deep3: p=1.951 primitivo e 1.954 conservativo;
- smooth, deep3: p=1.976 e 1.978;
- jitter, deep3: p=1.060 e 1.060;
- ring1 regular/smooth: circa uno per entrambi;
- tutti i nodi interni regular/smooth: circa 1.59–1.60, non due.

Non si dichiara verificato né l'ordine due globale del residuale né quello
della soluzione evoluta. La FASE 6B.2 resta una proposta da revisionare.
NODE rimane la rappresentazione implementata; la questione matematica generale
NODE/DUAL_MEAN non viene chiusa da questo esperimento.

## Soluzione analitica condivisa

In [euler_common_mms.hpp](../UnitTests/SU2_CFD/euler_common_mms.hpp) si definiscono
gamma=1.4 e, su [0,1]^2,

```
a(t)=0.1 cos(t), b(t)=0.1 sin(t)
rho=1+a sin(2pi x)sin(2pi y)
u=0.5+a cos(2pi x)sin(2pi y)
v=0.3+a sin(2pi x)cos(2pi y)
p=1+b cos(2pi x)cos(2pi y).
```

Evaluate restituisce W, W_t, W_x e W_y analitici. Per ogni tempo reale rho e p
sono >=0.9. In particolare p_t è non nullo già a t=0, mentre il gradiente
spaziale di p è nullo a t=0 ma non ai tempi 0.05 e 0.1.

La conversione e il Jacobiano sono

```
U(W) = (rho, rho*u, rho*v, p/(gamma-1)+rho*(u*u+v*v)/2)
J = [ 1                 0       0       0
      u                 rho     0       0
      v                 0       rho     0
      (u*u+v*v)/2       rho*u   rho*v   1/(gamma-1) ].
```

La sorgente primitiva implementa direttamente le quattro equazioni espanse:

```
S_rho = rho_t + u*rho_x + v*rho_y + rho*(u_x+v_y)
S_u   = u_t + u*u_x + v*u_y + p_x/rho
S_v   = v_t + u*v_x + v*v_y + p_y/rho
S_p   = p_t + u*p_x + v*p_y + gamma*p*(u_x+v_y).
S_U   = J*S_W.
```

La funzione ConservativePDE differenzia indipendentemente i prodotti in
U_t+F_x+G_y; non usa J né PrimitiveSource. Flux implementa i flussi fisici
analitici esclusivamente per verificare queste formule, non un flusso numerico.
Nessuna sorgente è ricavata dal residuale discreto.

### Test analitici

[euler_mms_analytic_tests.cpp](../UnitTests/SU2_CFD/euler_mms_analytic_tests.cpp)
copre 11 punti non speciali e t=0,0.037,0.1,0.283:

- derivate analitiche di W contro differenze finite centrate a cinque punti,
  con h=1e-3,5e-4,2.5e-4 e aritmetica long double;
- J*S_W contro differenziazione esplicita indipendente dei prodotti conservativi;
- sorgente conservativa contro differenze finite di U_t+F_x+G_y;
- S_W contro W_t+A_x W_x+A_y W_y, usando le matrici fisiche del kernel esistente;
- ogni colonna di J contro differenziazione indipendente di U(W).

Risultato: **3.300 assert superati**. Massimo errore sulle derivate di W:
3.14616e-11; massimo disaccordo fra le due formule di S_U: 4.33681e-19;
massimo scarto dell'identità primitiva con le matrici fisiche: 2.22045e-16.
Le tolleranze di arrotondamento della verifica long double tengono conto
dell'epsilon della piattaforma. Le differenze finite non vengono usate come
sorgenti o riferimenti nel confronto semidiscreto.

## Harness e percorsi SU2

[euler_mms_audit.cpp](../UnitTests/SU2_CFD/euler_mms_audit.cpp) istanzia direttamente
CEulerSolver e CPrimitiveEulerSolver, senza sovrascriverne metodi numerici.
Un adattatore **test-only** Manufactured, derivato da CVerificationSolution,
fornisce GetSolution, GetBCState, GetMMSSourceTerm e IsManufacturedSolution=true.
Non è aggiunta alcuna nuova opzione MMS al driver o alla factory produttiva.

Percorso conservativo:

1. CConfig::SetGlobalParam(EULER,RUNTIME_FLOW_SYS), preprocessing P2P/periodico
   della geometria, come nel banco FASE 6A.1;
2. GetInitialCondition ai vertici per t=0; GetSolution per reinserire lo stato
   esatto ai tempi successivi. Set_OldSolution e BC_Custom;
3. CEulerSolver::Preprocessing e Upwind_Residual produttivi;
4. **Source_Residual produttivo**, che richiama l'adattatore MMS e sottrae V*S_U;
5. BC_Custom per applicare i vincoli Dirichlet e verificarli separatamente.

Riferimenti al repository:

- `CEulerSolver.cpp:1678–1708`, preprocessing MUSCL;
- `CEulerSolver.cpp:1803–1816`, selezione effettiva del Roe SIMD;
- `CFVMFlowSolverBase.inl:1582` e seguenti, EdgeFluxResidual;
- `numerics_simd/flow/convection/common.hpp:201–271`, ricostruzione di T,u,v,p
  e ricalcolo di densità/entalpia;
- `numerics_simd/flow/convection/upwind.hpp:203–297`, Roe e dissipazione;
- `CEulerSolver.cpp:2078–2327`, Source_Residual, con MMS alle righe 2290–2318;
- `CFVMFlowSolverBase.inl:1518–1577`, BC_Custom forte;
- `CPrimitiveEulerSolver.cpp:180–238`, operatore primitivo completo e sorgenti.

L'array numerics passato a Source_Residual contiene puntatori nulli: tutte le
sorgenti fisiche che li userebbero sono disabilitate e controllate prima della
chiamata; il percorso MMS non richiede oggetti numerici sorgente.
Non viene implementato un secondo solver conservativo.

Il primitivo riceve lo stesso campo W esatto e callback condivise per stato
esterno e S_W tramite ComputeSpatialResidual. Vengono usati NODE, AGGREGATED,
NODE_POINT, MUSCL affine, WLS, limiter NONE e configurazione Forward Euler.
La configurazione temporale non viene avanzata: il tempo del primitivo è
l'argomento esplicito della chiamata; quello conservativo è SetPhysicalTime.

### Fisica, unità e ricostruzioni

Entrambi usano IDEAL_GAS, gamma=1.4, mesh statica, MGLEVEL=0, Roe,
MUSCL_KAPPA_FLOW=0 e WLS. Nel conservativo restano ROE_KAPPA=0.5 e
ENTROPY_FIX_COEFF=0.001, senza correzione low-Mach. Non si modificano formule,
gradienti, limiter o dissipazione per far coincidere gli operatori.

REF_DIMENSIONALIZATION=DIMENSIONAL con rho_ref=p_ref=velocity_ref=time_ref=1.
Il costruttore CEulerSolver compila questi riferimenti. L'API primitiva usa
direttamente W,x,t e non esegue quella non dimensionalizzazione: il solo harness
dichiara esplicitamente i suoi metadati di riferimento uguali a uno, senza
riscalare né stati né sorgenti. Lo smoke test iniziale ha evidenziato questa
distinzione, che è stata corretta prima della campagna accettata.

MACH_NUMBER=0.4 serve al costruttore conservativo; il campo freestream iniziale
è poi interamente sostituito dalla soluzione analitica, anche sui nodi di bordo.
Non si confonde quel Mach con il Mach locale della MMS. Non si usano far field.

La differenza fra ricostruzioni resta effettiva: il conservativo ricostruisce
T,u,v,p e deriva rho; il primitivo ricostruisce rho,u,v,p. Inoltre cambiano la
formulazione dell'operatore e la chiusura del bordo. Il confronto non isola
uno solo di questi meccanismi e non attribuisce tutte le differenze a T contro rho.

## Geometria, inizializzazione e controlli preliminari

Si riutilizza MeshFixture e si congela, in
[euler_mms_test_support.hpp](../UnitTests/SU2_CFD/euler_mms_test_support.hpp),
una copia dei generatori e della quadratura indipendente già verificati in 6A.1.
I file 6A.1 non vengono modificati. Le famiglie sono:

- regular: due triangoli per quadrato con diagonale fissa;
- smooth: deformazione `(0.08 sin(pi*x)sin(pi*y),0.06 sin(pi*x)sin(pi*y))`;
- jitter: perturbazioni interne deterministiche
  `(0.18h sin(2.1i+1.3j),0.16h cos(1.7i-0.9j))`.

N=8,16,32,64, h=1/N. Nessun triangolo invertito o degenere.
Le due geometrie native sono costruite separatamente dagli stessi dati; vengono
confrontate tutte le coordinate e tutti i volumi. La corrispondenza con gli
indici del generatore è verificata per ogni nodo. Le metriche primitive sono
confrontate a loro volta con i volumi del conservativo.

Massimi sull'intera campagna accettata:

| Verifica | Massimo scarto |
|---|---:|
| Coordinate delle due geometrie | 0 |
| Volumi delle due geometrie e metriche primitive | 0 |
| W iniziale/reinserito primitivo | 0 |
| U iniziale/reinserito conservativo | 8.88178e-16 |
| Primitive ricavate dal vero stato conservativo | 6.66134e-16 |
| Inserimento sorgente conservativa integrata | 2.77556e-17 |
| Sorgente conservativa recuperata dividendo per V | 3.36711e-14 |
| Sorgente primitiva restituita dal solver | 6.66134e-16 |
| Identità numerica J*(L_W-W_t)=J*L_W-U_t | 2.63678e-16 |
| Errore dei tempi di sorgenti/bordi | 0 |
| Stato Dirichlet imposto rispetto al riferimento | 4.44089e-16 |

Si verifica esplicitamente `R_after-R_flux=-V*S_U` per ogni nodo/componente,
prima di applicare il vincolo di bordo; questo controlla segno, un solo fattore
volume, coordinate e tempo. Per il primitivo si controlla result.source rispetto
alla sorgente analitica nodale. Sono controllati il numero di callback di bordo
e di sorgente e i rispettivi tempi, non soltanto i valori risultanti.

Non sono intervenuti fallback non fisici del conservativo; la contrazione di
positività primitiva è esattamente uno in tutti i casi. Tutti gli stati analitici
sono ammissibili. Non sono introdotti clipping o correzioni di residuale.

### Medie duali diagnostiche

Ogni contributo al duale è suddiviso in sottotriangoli vertice–punto medio–
baricentro, controllando i volumi contro SU2. La quadratura Gauss-Duffy è
verificata sui sei monomi 1,x,y,x²,xy,y² con ordini 3,5,8 e momenti baricentrici
esatti. Si calcolano indipendentemente media(U(W)) e media(W), con ordine 8 e
confronto 5-vs-8 per entrambi i campi.

Errore polinomiale massimo: 3.33067e-15. Massima differenza delle medie 5-vs-8:
5.01377e-13. Le medie restano diagnostiche: non inizializzano gli stati e non
entrano nelle sorgenti. I CSV riportano media-punto per U e W e il massimo di
`media(U(W))-U(media(W))`, che raggiunge 7.71067e-4 sulla campagna.
Con u,v,p variabili questa differenza non è più degenere come nel contatto 6A.1.
Non si possono scambiare le due medie o dedurre DUAL_MEAN dalla sola geometria.

## Convenzioni del difetto e maschere di confronto

```
tau_W = result.rhs - W_t
tau_U = -LinSysRes/V - U_t
U_t   = J(W_exact)*W_t
difference_U = J(W_exact)*tau_W - tau_U.
```

Gli operatori includono le rispettive sorgenti. Tau_U è campionato solo sui
DOF liberi. Non si sottrae U_t da un residuale azzerato per Dirichlet dichiarando
il risultato errore spaziale sul bordo.

Il primitivo usa stato esterno analitico sui punti delle mezze facce fisiche e
fluttuazione entrante. Il conservativo impone lo stato puntuale ai vertici di
bordo e azzera il residuale. La riga Dirichlet_constraint misura U_numerico-U_esatto
sui soli nodi vincolati: sono al più dell'ordine dell'arrotondamento
(Linf massimo 4.44089e-16), senza interpretarle come accuratezza di una chiusura
spaziale conservativa debole. Il residuale vincolato è invece esattamente zero.

Per componente si esportano:

```
RMS   = sqrt(sum(e_i²)/count)
L1_V  = sum(V_i*abs(e_i))/sum(V_i)
L2_V  = sqrt(sum(V_i*e_i²)/sum(V_i))
Linf  = max(abs(e_i))
p     = log(E_N/E_2N)/log(2).
```

Regioni: all, boundary, ring1, ring2, deep3 (almeno tre strati esclusi),
fixed ([0.25,0.75]^2 sulle coordinate effettive), interior (tutti i nodi non
vincolati). Per tau_U, J_tau_W e difference_U, `all_free` è il globale sui DOF
liberi ed equivale a interior; non equivale al globale `all` del primitivo.
Non viene inventata una norma globale conservativa che includa difetti di bordo
non definiti. Count e volume riportano esplicitamente il dominio di ogni norma.
Gli ordini di norme identicamente nulle e dei vincoli Dirichlet imposti entro
arrotondamento restano vuoti; non si impone ordine due.

## Tabelle dei risultati

[CSV completo delle norme e degli ordini](primitive_euler_phase6b1_mms_errors.csv):
5.760 righe per componenti, tempi, regioni e tipi di difetto/diagnostica.
[CSV dei controlli](primitive_euler_phase6b1_mms_checks.csv): 744 righe.
Nei controlli la misura scalare è nella colonna rms; le altre colonne non
costituiscono norme/ordini. Nei risultati le componenti W sono rho,u,v,p;
quelle U sono rho,rho*u,rho*v,rho*E. Non si confronta la componente p direttamente
con rho*E: occorre J. Le tabelle seguenti usano L2 volumetrica.

### Densità, t=0.1, N=64 e ordini 32→64

| Famiglia | Regione | tau_W L2 | p | tau_U L2 | p | differenza_U L2 | p |
|---|---|---:|---:|---:|---:|---:|---:|
|regular|all|3.74348e-03|1.596|—|—|—|—|
|regular|boundary|1.64952e-02|1.120|—|—|—|—|
|regular|ring1|8.44038e-03|0.976|8.59172e-03|0.978|5.97896e-04|1.049|
|regular|ring2|8.43385e-04|2.092|8.54592e-04|2.101|6.79581e-05|2.137|
|regular|deep3|1.19918e-03|1.951|1.21285e-03|1.954|7.59210e-05|2.106|
|regular|fixed|1.15842e-03|1.921|1.16679e-03|1.915|7.28909e-05|2.090|
|regular|interior|2.39909e-03|1.591|2.43859e-03|1.593|1.66347e-04|1.676|
|smooth|all|4.41141e-03|1.579|—|—|—|—|
|smooth|boundary|1.95056e-02|1.088|—|—|—|—|
|smooth|ring1|1.00527e-02|0.997|1.02491e-02|0.996|7.55568e-04|1.037|
|smooth|ring2|9.43556e-04|1.987|9.53567e-04|1.992|7.89769e-05|2.130|
|smooth|deep3|1.32542e-03|1.976|1.32626e-03|1.978|7.59202e-05|2.176|
|smooth|fixed|1.36126e-03|1.967|1.34848e-03|1.958|6.91524e-05|2.226|
|smooth|interior|2.81229e-03|1.598|2.85680e-03|1.596|2.02740e-04|1.647|
|jitter|all|7.88195e-03|1.215|—|—|—|—|
|jitter|boundary|1.82501e-02|1.117|—|—|—|—|
|jitter|ring1|1.08205e-02|0.973|1.09446e-02|0.975|6.20801e-04|1.068|
|jitter|ring2|6.84182e-03|1.337|6.85212e-03|1.339|2.04706e-04|1.606|
|jitter|deep3|7.02710e-03|1.060|7.06847e-03|1.060|4.14164e-04|1.041|
|jitter|fixed|6.96958e-03|1.030|7.01206e-03|1.032|3.85224e-04|0.995|
|jitter|interior|7.31120e-03|1.124|7.35815e-03|1.125|4.20775e-04|1.095|

Il globale libero all_free del conservativo coincide con la riga interior; i trattini indicano un difetto non definito sui DOF vincolati, non un errore zero.

### Confronto nello spazio U, interno fisso a t=0.1, N=64

| Famiglia | componente U | J*tau_W L2 | p | tau_U L2 | p | differenza_U L2 | p |
|---|---|---:|---:|---:|---:|---:|---:|
|regular|rho|1.15842e-03|1.921|1.16679e-03|1.915|7.28909e-05|2.090|
|regular|rho*u|1.02833e-03|2.077|1.03440e-03|2.070|3.28366e-05|2.119|
|regular|rho*v|7.86252e-04|2.127|7.88983e-04|2.119|1.90158e-05|2.284|
|regular|rho*E|4.52453e-03|1.994|4.52368e-03|1.991|2.56815e-05|2.969|
|smooth|rho|1.36126e-03|1.967|1.34848e-03|1.958|6.91524e-05|2.226|
|smooth|rho*u|1.26279e-03|2.173|1.24075e-03|2.172|4.77520e-05|2.437|
|smooth|rho*v|9.85878e-04|2.234|9.63036e-04|2.235|3.85212e-05|2.459|
|smooth|rho*E|5.30515e-03|2.037|5.27522e-03|2.034|6.04238e-05|2.619|
|jitter|rho|6.96958e-03|1.030|7.01206e-03|1.032|3.85224e-04|0.995|
|jitter|rho*u|6.45314e-03|1.042|6.51954e-03|1.044|8.10759e-04|1.032|
|jitter|rho*v|4.70422e-03|1.075|4.69369e-03|1.074|4.45964e-04|0.987|
|jitter|rho*E|2.53106e-02|1.033|2.52200e-02|1.033|6.47119e-04|1.051|

### Dipendenza dal tempo: densità, deep3, N=64, ordini 32→64

| Famiglia | t | p(tau_W) | p(tau_U) | p(differenza_U) |
|---|---:|---:|---:|---:|
|regular|0|1.941|1.944|2.107|
|regular|0.05|1.946|1.949|2.106|
|regular|0.1|1.951|1.954|2.106|
|smooth|0|1.969|1.971|2.183|
|smooth|0.05|1.973|1.974|2.180|
|smooth|0.1|1.976|1.978|2.176|
|jitter|0|1.061|1.061|1.032|
|jitter|0.05|1.061|1.061|1.036|
|jitter|0.1|1.060|1.060|1.041|

### Raffinamento completo: densità nell’interno fisso, t=0.1

| Famiglia | N | tau_W L2 | p | tau_U L2 | p | differenza_U L2 | p |
|---|---:|---:|---:|---:|---:|---:|---:|
|regular|8|4.64690e-02|—|4.72236e-02|—|5.01066e-03|—|
|regular|16|1.47848e-02|1.652|1.48662e-02|1.667|1.42109e-03|1.818|
|regular|32|4.38749e-03|1.753|4.39924e-03|1.757|3.10241e-04|2.196|
|regular|64|1.15842e-03|1.921|1.16679e-03|1.915|7.28909e-05|2.090|
|smooth|8|8.72663e-02|—|8.69070e-02|—|5.12480e-03|—|
|smooth|16|1.98206e-02|2.138|1.95988e-02|2.149|1.60085e-03|1.679|
|smooth|32|5.32220e-03|1.897|5.23795e-03|1.904|3.23489e-04|2.307|
|smooth|64|1.36126e-03|1.967|1.34848e-03|1.958|6.91524e-05|2.226|
|jitter|8|8.16750e-02|—|8.18124e-02|—|4.53103e-03|—|
|jitter|16|3.40144e-02|1.264|3.40207e-02|1.266|1.88891e-03|1.262|
|jitter|32|1.42316e-02|1.257|1.43393e-02|1.246|7.67932e-04|1.299|
|jitter|64|6.96958e-03|1.030|7.01206e-03|1.032|3.85224e-04|0.995|

## Interpretazione scientifica

Le verifiche analitiche, fisiche e di inserimento delle sorgenti sono superate.
La differenza fra i due difetti non è una differenza di soluzione manufactured,
normalizzazione o convenzione dei segni: entrambe le rappresentazioni hanno
lo stesso riferimento continuo e sorgenti trasformate coerentemente.

Su regular/smooth il primo anello è circa O(h), ring2 e deep3 sono vicini a
O(h²) per la densità. Una fascia di errore O(h), con misura O(h), contribuisce
O(h^(3/2)) alla L2 globale. Gli ordini circa 1.59–1.60 dei nodi liberi sono
compatibili con questa contaminazione e con contributi interni O(h²) ancora
presenti alle risoluzioni finite. Anche il conservativo ha un primo anello
limitante, pur imponendo esattamente lo stato al bordo. Per il primitivo si
aggiunge la sua chiusura debole, misurata separatamente.

Su jitter il difetto è circa O(h) anche nell'interno profondo per entrambi:
non è un problema confinato al bordo. La deformazione liscia non equivale a
perturbazioni O(h) con frequenza legata agli indici; le famiglie non sono fuse.
Le osservazioni sono coerenti con gli audit dei gradienti precedenti, ma questo
esperimento non effettua sostituzioni causali di gradienti o quadrature.

J_tau_W-tau_U decresce nelle regioni libere. Sulle mesh regolari/deformate
nell'interno la differenza è circa O(h²) per la densità; su jitter è circa O(h).
Questo è un confronto coerente tra difetti, non un'identità discreta o una prova
di equivalenza delle soluzioni evolute. Gli ordini maggiori di due osservati in
alcune differenze/componenti sono cancellazioni/transitori da riportare, non
una nuova dichiarazione d'ordine dello schema.

Per esempio, su regular N=64, t=0.1, fixed, la norma del difetto primitivo di u
è 4.58582e-4. Non va confrontata direttamente con quella di rho*u conservativa
(1.03440e-3). Dopo J, la componente confrontabile vale 1.02833e-3 e la norma
della differenza dei campi di difetto è 3.28366e-5. Quest'ultima non è la
semplice differenza delle due norme: è la norma del loro scarto nodo per nodo.

La verifica ufficiale SU2 Navier–Stokes riguarda errori della soluzione
stazionaria. Il presente difetto semidiscreto Euler non è quella quantità.
Non si trasferisce automaticamente il suo secondo ordine a questo esperimento,
né si deduce l'ordine della soluzione evoluta dal solo ordine di tau.

## Proposta operativa FASE 6B.2 — da revisionare

1. Riutilizzare questi componenti analitici e le stesse mesh; inizializzare
   puntualmente W e U a t=0, mantenendo gamma, unità, schemi e limiter invariati.
2. Verificare a mesh fissata che il ciclo temporale valuti sorgenti e condizioni
   esatte a t_n; per il conservativo gestire esplicitamente anche il vincolo
   Dirichlet a t_(n+1), distinguendolo dal residuale dei DOF liberi.
3. Usare Forward Euler nativo nei due solver, con stesso passo fisico uniforme
   dt=T/Nt e T iniziale 0.1. Non identificare uguaglianza di CFL con uguaglianza
   del passo; non introdurre clipping, riduzioni nascoste o integratori nuovi.
4. Eseguire dt,dt/2,dt/4 a mesh fissata. Richiedere differenza delle soluzioni
   alle ultime due risoluzioni temporali <5% dell'errore spaziale risolvibile,
   per le regioni/componenti rilevanti, e variazione degli ordini <0.1;
   usare una soglia assoluta dichiarata soltanto quando domina l'arrotondamento.
5. Misurare separatamente W_prim-W_exact, W(U_cons)-W_exact,
   U_prim-U_exact, U_cons-U_exact e la differenza fra solver nella stessa
   rappresentazione, con tutte le regioni e norme già definite. Non pretendere
   coincidenza a mesh finita né ordine due globale a priori.
6. Registrare minimi rho/p, fallback/contrazioni, passi effettivi e bilanci
   conservativi diagnostici, contabilizzando il bordo forte. Confrontare i due
   errori con l'analitica prima di interpretare la loro differenza reciproca.

Se il secondo ordine globale della soluzione non emerge, documentare il limite
e fermarsi alla revisione scientifica: questa consegna non autorizza correzioni
al solver o una campagna shock. Nessuna FASE 6B.2 è stata eseguita ora.

## Comandi, regressioni e stato Git

Comandi eseguiti dalla radice del repository:

```sh
python3 UnitTests/run_euler_mms_audit.py --output-dir /private/tmp/su2-mms-6b1-final
/private/tmp/su2-mms-6b1-final/euler_mms_analytic --reporter compact
python3 UnitTests/run_euler_mms_audit.py --output-dir /private/tmp/su2-mms-6b1-complete --target harness
python3 UnitTests/export_euler_mms_audit.py --executable /private/tmp/su2-mms-6b1-complete/euler_mms_audit --output-dir /private/tmp/su2-mms-6b1-complete/campaign --jobs 2
/private/tmp/su2-phase5-tests/test_driver --reporter compact > /private/tmp/su2-mms-6b1-final/regression.log 2>&1
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-mms-6b1-debug/cli --cli
python3 -m py_compile UnitTests/run_euler_mms_audit.py UnitTests/export_euler_mms_audit.py
git diff --check
```

Il build script riutilizza compilatore, flag e dipendenze native della build
release esistente, collegandosi alle vere librerie SU2. Compila separatamente
il test analitico Catch2 e l'harness MPI. `--target harness` evita di ricompilare
il test analitico quando cambiano solo diagnostiche dell'harness. Per riprodurre
tutto si può costruire il target predefinito both in una sola directory.
Gli ordini sono calcolati dall'exporter, senza soglie d'ordine atteso.
I log dei casi e CSV grezzi restano nella directory temporanea della campagna;
`--collect-only` consente di rigenerare il postprocessing senza eseguire SU2.

Risultati delle verifiche:

- test analitico nuovo: 1 caso, **3.300 assert** superati;
- suite precedente non modificata: **70 casi, 22.253.237 assert** superati sul
  driver già costruito (nessun sorgente produttivo dipendente è cambiato);
- quattro controlli CLI: primitivo bloccato come previsto, Euler operativo,
  harmonic balance e species non supportati respinti come previsto;
- campagna finale: tutti i 12 harness completati, tre tempi ciascuno;
- compilazione Python e controlli di formato CSV/diff verificati.

L'unico file preesistente modificato è UnitTests/meson.build, con l'aggiunta
del nuovo test analitico all'elenco dei test. Non si cambia alcun test o soglia
precedente. Tutti i sorgenti produttivi, driver e kernel restano invariati.
Gli undici file non tracciati della FASE 6A.1 sono preservati integralmente.
Il confronto SHA256 prima/dopo conferma che tutti e undici sono invariati.
I nuovi file 6B.1 sono: due header, test analitico, harness, due script, questo
rapporto e due CSV. Tutte le modifiche nel repository sono applicate con
apply_patch. Nessun commit, push, pull, sync o cambio branch; staging vuoto.
