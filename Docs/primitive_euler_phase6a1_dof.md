# FASE 6A.1 — Significato operativo dei DOF del vero CEulerSolver

## Ambito e stato

Esperimento sul branch `feature/primitive-path-solver`, base
`5be7a6e72a6c082cce192e73e8dd53eb17f1b13b` (FASE 5), working tree inizialmente
pulito. Nessuna modifica ai solver produttivi, al kernel primitivo, agli schemi
conservativi o ai test/soglie delle fasi precedenti. Nessun commit/push/pull/sync.
FASE 6B non avviata. La decisione matematica NODE/DUAL_MEAN resta aperta.

## Benchmark e riferimenti esatti

Dominio [0,1]^2, gamma=1.4, u=0.5, v=0.3, p=1, tempo finale T=0.05.
Il benchmark principale è il contatto liscio

```
rho(x,y,t) = 1 + 0.2 sin(2 pi (x-0.5t)) cos(2 pi (y-0.3t)).
```

È una soluzione esatta delle quattro equazioni Euler senza sorgente: rho è
trasportata a velocità costante e pressione/velocità restano uniformi.
Il campo è non affine e dipende dal tempo. Non occorre implementare una nuova
MMS né compensare discretamente il residuale.

Un controllo indipendente, aggiunto dopo la prima campagna e conservato
separatamente, usa `rho=1+0.2(x-0.5t)+0.15(y-0.3t)` sulla famiglia perturbata.
Qui la differenza media-punto della densità è esattamente
`0.2 delta_x + 0.15 delta_y`: la curvatura non può essere la causa di D.
Questo controllo non sostituisce il benchmark sinusoidale.

Per ogni nodo sono calcolati separatamente U_point, media(U(W)) e media(W),
all'inizio e/o al tempo finale come richiesto dall'inizializzazione e dalle
metriche. U=(rho,rho*u,rho*v,p/(gamma-1)+rho*(u*u+v*v)/2).
Non si sostituisce mai U(media(W)) a media(U(W)).

Limite intenzionale: per questi contatti, con u,v,p costanti, U dipende
affinemente da rho. Quindi media(U(W)) e U(media(W)) coincidono analiticamente.
Il loro confronto numerico controlla la quadratura, ma NON discrimina le due
medie primitive/conservative in un flusso generale con velocità variabile.
Le quattro componenti conservative non costituiscono quattro prove indipendenti:
gli errori esatti sono proporzionali a quelli di densità con fattori
1, 0.5, 0.3 e 0.17; il codice esporta comunque tutte le componenti.

## Percorso produttivo effettivamente eseguito

`UnitTests/SU2_CFD/euler_dof_audit.cpp` istanzia direttamente **CEulerSolver**,
senza sovrascriverne metodi numerici. Usa CConfig, CPhysicalGeometry e il lettore
mesh nativo tramite il fixture già esistente. Il preprocessing P2P e periodico
viene eseguito anche con un unico rank, come richiesto dall'infrastruttura SU2.
SetGlobalParam(EULER,RUNTIME_FLOW_SYS) attiva il contesto del driver.

Configurazione: IDEAL_GAS, gamma=1.4, DIMENSIONAL, rho_inf=p_inf=1;
MGLEVEL=0, TIME_STEPPING, EULER_EXPLICIT, ROE, MUSCL=YES, kappa MUSCL=0,
WEIGHTED_LEAST_SQUARES, limiter NONE. Sono mantenuti i default produttivi
ROE_KAPPA=0.5, ENTROPY_FIX_COEFF=0.001, LOW_MACH_CORR=NO.
Nessuna modifica alla dissipazione, alla ricostruzione di T,u,v,p o al fallback.
Riferimenti:

- `CEulerSolver.cpp:1678–1708`: preprocessing e gradienti MUSCL.
- `CEulerSolver.cpp:1803–1816`: Roe gas ideale passa al percorso SIMD.
- `CFVMFlowSolverBase.inl:1582` e seguenti: EdgeFluxResidual produttivo.
- `numerics_simd/flow/convection/common.hpp:201–271`: ricostruzione delle
  primitive SU2 e controlli di ammissibilità.
- `numerics_simd/flow/convection/upwind.hpp:203–297`: Roe e dissipazione.
- `CFVMFlowSolverBase.hpp:803–845`: Forward Euler conservativo R/V.
- `CPhysicalGeometry.cpp:6683–6797`: costruzione di normali e volumi duali.

Il solo adattatore analitico è una sottoclasse **test-only** di
CVerificationSolution assegnata al puntatore pubblico VerificationSolution;
fornisce GetSolution e GetBCState. Nessuna classe MMS produttiva viene modificata.

Ogni passo esegue: Set_OldSolution, BC_Custom al tempo corrente,
Preprocessing, Upwind_Residual, BC_Custom per azzerare i residuali vincolati,
ExplicitEuler_Iteration e BC_Custom al nuovo tempo.
Il passo uniforme è assegnato ai nodi; il metodo di aggiornamento resta quello
produttivo. Source_Residual non è necessario: tutte le sorgenti fisiche sono nulle.
Si verificano gamma, riferimenti rho/p/tempo uguali a 1, tempo finale effettivo,
assenza di fallback nodali/di ricostruzione, finitezza e positività.
Il cambiamento della soluzione è controllato esplicitamente: non è un confronto
della sola inizializzazione.

### Bordo e inizializzazione

Tutti i lati usano il vero BC_Custom manufactured forte:
`CFVMFlowSolverBase.inl:1518–1577`. I nodi di bordo sono imposti al valore
puntuale analitico e il loro residuale è azzerato. L'imposizione al nuovo tempo
è una gestione esplicita del vincolo del banco di prova, non un cambiamento
dell'algoritmo Euler nell'interno. Evita un ritardo di un passo nei dati esatti.

Ogni mesh è ripetuta partendo da U_point oppure dalla media conservativa esatta.
Anche nel secondo caso i nodi Dirichlet vengono subito ricondotti al punto:
le due inizializzazioni differiscono quindi sui DOF liberi, non sul bordo.
Il confronto delle medie sul bordo NON può essere usato come prova autonoma
dell'interpretazione del DOF.

La regione fisica fissa è [0.25,0.75]^2, selezionata sulle coordinate deformate.
T=0.05 limita la distanza caratteristica fisica a meno di 0.1; questa regione è
ben separata dal bordo nel problema continuo. Ciò non costituisce una prova
di perfetto isolamento dal bordo nel problema discreto.

## Geometria e quadratura indipendente

N=8,16,32,64 su tre famiglie; N indica intervalli per lato, h=1/N:

- regular: due triangoli con diagonale costante per ogni quadrato cartesiano;
- smooth: x'=x+0.08 sin(pi*x)sin(pi*y), y'=y+0.06 sin(pi*x)sin(pi*y);
- jitter: soli nodi interni spostati di
  `(0.18h sin(2.1i+1.3j), 0.16h cos(1.7i-0.9j))`.

Ogni triangolo ha area orientata positiva >0.1h^2. Per ogni vertice del triangolo
il contributo al duale è suddiviso nei due sottotriangoli vertice–punto medio–
baricentro. La loro somma riproduce il volume nativo SU2 con tolleranza relativa
2e-12. Questo controllo non usa CAggregatedGeometry né un operatore primitivo.

Si calcolano i momenti attorno al vertice:

```
delta = (1/V) integral (x-x_i) dV
M_ab  = (1/V) integral (x_a-x_i,a)(x_b-x_i,b) dV.
```

Volume e momenti vengono integrati analiticamente attraverso le identità
baricentriche dei triangoli. Per una sottocella di area A e vertici relativi z_k,
`integral z = A sum(z_k)/3` e
`integral z_a z_b = A[(sum z_a)(sum z_b)+sum(z_k,a z_k,b)]/12`.

Le medie dei campi usano invece quadratura tensoriale Gauss-Legendre dopo
trasformazione Duffy sul triangolo, con Jacobiano `2A(1-r)`. Gli ordini 3,5,8
sono verificati separatamente su **1,x,y,x^2,xy,y^2** rispetto ai momenti
analitici. Le medie conservative sinusoidali a t=0,T sono confrontate con
ordini 5 e 8; le medie finali usano ordine 8.

Errore massimo sui monomi, campagna N<=64: 3.33e-15; incluso N=128: 3.45e-15.
Differenza massima delle medie U fra ordine 5 e 8: 3.47e-14.
Per il controllo affine: rispettivamente 2.89e-15 e 1.16e-14.
La soglia di asimmetria significativa è dichiarata prima della classificazione:
`|delta|/h > 0.01`. Il CSV riporta la frazione nodale, non volumetrica.
M2 indica la norma di Frobenius del tensore; il CSV geometrico contiene anche
le norme delle tre componenti xx,xy,yy.

## Norme, regioni e controllo temporale

Per ogni componente: RMS=sqrt(sum e^2/N_g),
L2=sqrt(sum V_i e_i^2/sum V_i), Linf=max|e_i|.
Sono esportati Epoint, Emean e D, per all, boundary, ring1, ring2, deep3
(almeno tre strati esclusi), fixed e interior (tutti i DOF non sul bordo).
Ordini: log(E_N/E_2N)/log(2), senza imporre un valore atteso.

Si usa `steps=ceil(T*N*N/0.25)*divisor`, dt=T/steps, divisor=1,2,4.
I campi finali completi vengono confrontati fra dt e dt/2 e fra dt/2 e dt/4;
le differenze sono esportate per componente, regione e tutte e tre le norme.
Si verifica che l'incertezza temporale sia piccola rispetto sia a Epoint sia a D.
Sul bordo forte la differenza temporale è zero per imposizione, non per
accuratezza dell'operatore; questi rapporti non sono usati per inferenze.

N<=64, massimo delta24/Epoint fuori dal bordo: 1.62% (sinusoidale),
1.35% (affine); massimo delta24/D: 0.78% e 0.14%.
Il passo minimo ha dunque un effetto molto inferiore alla separazione punto/media.
Non si usa Richardson per alterare la soluzione numerica.
Su jitter N=128: massimo delta24/Epoint fuori bordo 0.641%,
massimo delta24/D 0.030%. In fixed, rho, inizializzazione point:
delta12=3.41844e-7, delta24=1.70922e-7, ordine temporale 1.000005.
Con inizializzazione mean: 3.46064e-7, 1.73028e-7, ordine 1.000035.
Entrambe le dimezzature confermano il comportamento Forward Euler.

I 102 casi evolutivi analizzati sono 72 principali, sei aggiuntivi a N=128 e
24 controlli affini; i 12 controlli geometrici finali non avanzano la soluzione.
Nel sinusoidale il minimo rho finale fra i casi registrati è 0.79991281,
il minimo p è 0.9999999999999923. Non sono intervenuti fallback o clipping.
Questo verifica ammissibilità per questi dati, non una proprietà TVD/positivity
preserving generale. L'assenza di limiter consente piccoli overshoot.

## Risultati numerici

Le tabelle complete per componenti e regioni sono negli artefatti CSV allegati
a questo rapporto. Le tabelle riassuntive sotto usano la densità e L2 volumetrica.

- [Risultati sinusoidali e ordini, incluso N=128](primitive_euler_phase6a1_dof_errors.csv).
- [Controllo temporale sinusoidale](primitive_euler_phase6a1_dof_temporal.csv).
- [Effetto dell'inizializzazione](primitive_euler_phase6a1_dof_initialization.csv).
- [Risultati affini](primitive_euler_phase6a1_dof_affine.csv),
  [tempo](primitive_euler_phase6a1_dof_affine_temporal.csv) e
  [inizializzazione](primitive_euler_phase6a1_dof_affine_initialization.csv).
- [Momenti geometrici completi](primitive_euler_phase6a1_dof_geometry.csv).

Nel CSV gli indici delle componenti conservative sono 0=rho,1=rho*u,2=rho*v,
3=rho*E; per moment_component sono 0=xx,1=xy,2=yy. Le colonne p_* sono gli
ordini rispetto alla mesh precedente della medesima famiglia/configurazione.
Per offset_over_h, extra è la frazione di nodi con offset/h>0.01;
per validation, rms è il massimo errore polinomiale e l2 il massimo 5-vs-8;
per evolution, rms=dt, l2=massimo cambiamento di rho, linf=min rho, extra=min p.
Queste due righe di metadati non sono norme/ordini di errore della soluzione.

### Momenti geometrici e separazione analitica (interno)

| Famiglia | N | RMS delta | max delta | RMS delta/h | max delta/h | % asimmetrici | RMS M2 | D L2 | p(D) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|regular|8|2.87078e-17|3.90148e-17|2.29662e-16|3.12119e-16|0.00|2.28753e-03|6.15149e-03|—|
|regular|16|2.93635e-17|4.02125e-17|4.69815e-16|6.43400e-16|0.00|5.71882e-04|1.58162e-03|1.960|
|regular|32|2.96229e-17|3.87327e-17|9.47934e-16|1.23944e-15|0.00|1.42970e-04|3.98134e-04|1.990|
|regular|64|2.93556e-17|3.88839e-17|1.87876e-15|2.48857e-15|0.00|3.57426e-05|9.97093e-05|1.997|
|smooth|8|4.32140e-03|6.68786e-03|3.45712e-02|5.35028e-02|85.71|2.44670e-03|7.18316e-03|—|
|smooth|16|1.08181e-03|1.72318e-03|1.73089e-02|2.75709e-02|71.56|6.13338e-04|1.86807e-03|1.943|
|smooth|32|2.69248e-04|4.34054e-04|8.61593e-03|1.38897e-02|34.55|1.53455e-04|4.71454e-04|1.986|
|smooth|64|6.70672e-05|1.08737e-04|4.29230e-03|6.95916e-03|0.00|3.83716e-05|1.18114e-04|1.997|
|jitter|8|1.45067e-02|2.13830e-02|1.16054e-01|1.71064e-01|97.96|2.46620e-03|1.00725e-02|—|
|jitter|16|7.48318e-03|1.13944e-02|1.19731e-01|1.82310e-01|99.11|6.22356e-04|4.56957e-03|1.140|
|jitter|32|3.78533e-03|5.69720e-03|1.21131e-01|1.82310e-01|99.58|1.56011e-04|2.30670e-03|0.986|
|jitter|64|1.90363e-03|2.85026e-03|1.21832e-01|1.82416e-01|99.52|3.90512e-05|1.17404e-03|0.974|
|jitter|128|9.54563e-04|1.42513e-03|1.22184e-01|1.82416e-01|99.57|9.76897e-06|5.93958e-04|0.983|

### Momenti al bordo (N=64)

| Famiglia | RMS delta | max delta | RMS delta/h | max delta/h | RMS M2 |
|---|---:|---:|---:|---:|---:|
|regular|4.54541e-03|6.44498e-03|2.90906e-01|4.12479e-01|3.57738e-05|
|smooth|4.63355e-03|6.49133e-03|2.96547e-01|4.15445e-01|3.84036e-05|
|jitter|4.58522e-03|6.58796e-03|2.93454e-01|4.21629e-01|3.66476e-05|

### Contatto sinusoidale: regione fisica fissa

| Famiglia | N | init | Epoint | p | Emean | p | D | p |
|---|---:|---|---:|---:|---:|---:|---:|---:|
|regular|8|point|3.54449e-03|—|2.95039e-03|—|6.09313e-03|—|
|regular|8|mean|9.20935e-03|—|3.21619e-03|—|6.09313e-03|—|
|regular|16|point|8.39426e-04|2.078|1.21510e-03|1.280|1.57537e-03|1.951|
|regular|16|mean|2.19732e-03|2.067|8.19148e-04|1.973|1.57537e-03|1.951|
|regular|32|point|1.94796e-04|2.107|3.87277e-04|1.650|3.97606e-04|1.986|
|regular|32|mean|4.91370e-04|2.161|1.93606e-04|2.081|3.97606e-04|1.986|
|regular|64|point|4.70860e-05|2.049|1.06726e-04|1.859|9.96663e-05|1.996|
|regular|64|mean|1.13596e-04|2.113|4.70151e-05|2.042|9.96663e-05|1.996|
|smooth|8|point|5.34190e-03|—|3.48742e-03|—|7.81035e-03|—|
|smooth|8|mean|1.21831e-02|—|4.69035e-03|—|7.81035e-03|—|
|smooth|16|point|1.22675e-03|2.123|9.75495e-04|1.838|1.96689e-03|1.989|
|smooth|16|mean|3.11262e-03|1.969|1.18955e-03|1.979|1.96689e-03|1.989|
|smooth|32|point|2.69836e-04|2.185|3.75358e-04|1.378|5.08657e-04|1.951|
|smooth|32|mean|7.25858e-04|2.100|2.61790e-04|2.184|5.08657e-04|1.951|
|smooth|64|point|6.10892e-05|2.143|1.11143e-04|1.756|1.27991e-04|1.991|
|smooth|64|mean|1.68008e-04|2.111|5.85280e-05|2.161|1.27991e-04|1.991|
|jitter|8|point|3.42615e-03|—|9.81029e-03|—|1.18907e-02|—|
|jitter|8|mean|1.38307e-02|—|4.95698e-03|—|1.18907e-02|—|
|jitter|16|point|1.23689e-03|1.470|4.72215e-03|1.055|5.06343e-03|1.232|
|jitter|16|mean|3.35862e-03|2.042|2.75693e-03|0.846|5.06343e-03|1.232|
|jitter|32|point|3.58844e-04|1.785|2.35483e-03|1.004|2.42394e-03|1.063|
|jitter|32|mean|8.72487e-04|1.945|1.99207e-03|0.469|2.42394e-03|1.063|
|jitter|64|point|1.05451e-04|1.767|1.19374e-03|0.980|1.20807e-03|1.005|
|jitter|64|mean|2.29845e-04|1.924|1.20709e-03|0.723|1.20807e-03|1.005|
|jitter|128|point|2.66705e-05|1.983|5.99818e-04|0.993|6.02217e-04|1.004|
|jitter|128|mean|6.02288e-05|1.932|6.13926e-04|0.975|6.02217e-04|1.004|

### Ultima transizione su jitter: globale e regioni interne (64→128)

| init | Regione | Epoint | p | Emean | p | D | p |
|---|---|---:|---:|---:|---:|---:|---:|
|point|all|3.38002e-05|1.860|6.15174e-04|1.027|6.16207e-04|1.035|
|point|boundary|2.89782e-17|-0.054|1.44316e-03|0.999|1.44316e-03|0.999|
|point|ring1|8.24384e-05|1.610|6.87846e-04|0.992|6.85741e-04|0.993|
|point|ring2|8.01681e-05|1.163|7.99402e-04|0.943|7.97774e-04|0.954|
|point|deep3|2.82875e-05|1.890|5.81393e-04|0.943|5.82732e-04|0.954|
|point|interior|3.40664e-05|1.872|5.92870e-04|0.973|5.93958e-04|0.983|
|point|fixed|2.66705e-05|1.983|5.99818e-04|0.993|6.02217e-04|1.004|
|mean|all|6.77924e-05|1.906|6.29328e-04|1.009|6.16207e-04|1.035|
|mean|boundary|2.89782e-17|-0.054|1.44316e-03|0.999|1.44316e-03|0.999|
|mean|ring1|9.54763e-05|2.077|6.89162e-04|0.918|6.85741e-04|0.993|
|mean|ring2|1.07596e-04|1.431|8.25227e-04|0.967|7.97774e-04|0.954|
|mean|deep3|6.55279e-05|1.885|5.96371e-04|0.926|5.82732e-04|0.954|
|mean|interior|6.83262e-05|1.917|6.07772e-04|0.956|5.93958e-04|0.983|
|mean|fixed|6.02288e-05|1.932|6.13926e-04|0.975|6.02217e-04|1.004|

### Controllo affine su jitter: regione fisica fissa

| N | init | Epoint | p | Emean | p | D | p |
|---:|---|---:|---:|---:|---:|---:|---:|
|8|point|9.48570e-06|—|3.37860e-03|—|3.38363e-03|—|
|8|mean|2.35975e-03|—|1.13880e-03|—|3.38363e-03|—|
|16|point|4.22489e-06|1.167|1.67581e-03|1.012|1.67818e-03|1.012|
|16|mean|8.11852e-04|1.539|9.63101e-04|0.242|1.67818e-03|1.012|
|32|point|1.72173e-06|1.295|8.23847e-04|1.024|8.24876e-04|1.025|
|32|mean|2.02336e-04|2.004|6.88729e-04|0.484|8.24876e-04|1.025|
|64|point|6.10704e-07|1.495|4.12500e-04|0.998|4.12850e-04|0.999|
|64|mean|3.06157e-05|2.724|4.09994e-04|0.748|4.12850e-04|0.999|

### Memoria dell'inizializzazione: differenza dei campi finali (fixed, rho)

| Famiglia | N | L2 U_point-init meno U_mean-init |
|---|---:|---:|
|regular|8|5.79644e-03|
|regular|16|1.56451e-03|
|regular|32|3.97380e-04|
|regular|64|9.96673e-05|
|smooth|8|7.06049e-03|
|smooth|16|1.96637e-03|
|smooth|32|5.18001e-04|
|smooth|64|1.30669e-04|
|jitter|8|1.13015e-02|
|jitter|16|2.93019e-03|
|jitter|32|7.96618e-04|
|jitter|64|2.08486e-04|
|jitter|128|5.46595e-05|

## Interpretazione e limiti

La famiglia regular e la deformazione smooth non discriminano punto/media
nell'interno: D è già O(h^2). Il bordo rende D globale circa O(h^1.5), ma è
vincolato al punto e non può risolvere l'ambiguità matematica.
La famiglia jitter ha invece una frazione persistente e ampia di duali
asimmetrici e D circa O(h) nell'interno: non serve un'altra famiglia per ottenere
la separazione geometrica richiesta.

**Esito del criterio discriminante nell'interno fisso:** sul sinusoidale jitter,
64→128, D ha ordine 1.00435; Epoint ha ordine 1.98326 (point-init) e 1.93214
(mean-init), mentre Emean ha ordine 0.99289 e 0.97539. È evidenza numerica
favorevole al riferimento puntuale, compatibile con l'emergere di ordine due
per questo campo/regione. La separazione D=6.02217e-4 è circa 3500 volte la
differenza temporale delta24 e molto maggiore dell'incertezza di quadratura.
Non si basa soltanto su costanti: i tassi punto e media si separano.

Il vettore finale dipende ancora dall'inizializzazione: in fixed la differenza
L2 delle densità è 2.08486e-4 a N=64 e 5.46595e-5 a N=128 (ordine circa 1.93).
È confrontabile con Epoint, ma decresce più rapidamente di D=O(h). Pertanto
l'inizializzazione cambia un contributo del livello dell'errore puntuale senza
preservare nel campo finale l'intera separazione analitica punto/media.

**Limite globale:** non si dichiara ordine due globale. A N=128 su jitter gli
ordini globali Epoint sono 1.86048 e 1.90559 per le due inizializzazioni;
su regular/smooth, point-init a N=64 dà circa 1.74 globale pur superando 2 in
fixed. Il bordo imposto e gli anelli richiedono analisi distinta. Il risultato
ufficiale Navier–Stokes stazionario non viene esteso a questo Euler instazionario.

I risultati devono essere interpretati con il criterio del prompt, non con il
solo rapporto delle costanti. Il controllo affine mostra in modo particolarmente
chiaro che una distanza molto minore dal punto non prova da sola ordine due.
Anche il decadimento della memoria dell'inizializzazione non autorizza a
identificare automaticamente il DOF con una media o con un punto.
Nel controllo affine point-init, Epoint fixed ha ordine 1.49532 nell'ultima
transizione: il criterio O(h^2) non è stabilito. Mean-init mostra una transizione
rapida (ordine 2.72441) che non va scambiata per ordine asintotico dimostrato.
Questa seconda prova resta formalmente non conclusiva sul tasso puntuale.

Questi test riguardano un campo di contatto, WLS, nessun limiter, Roe produttivo
e bordo Dirichlet forte. Non dimostrano una proprietà universale del solver,
né ordine due per onde acustiche o per ogni triangolazione irregolare.
Non trasferiscono un significato di media conservativa a una media primitiva.
NODE resta implementato; DUAL_MEAN non viene introdotto.
La conclusione operativa è quindi circoscritta: il CEulerSolver testato mostra
una convergenza verso il riferimento puntuale più rapida di quella verso la
media sul contatto sinusoidale in un interno realmente asimmetrico, anche
partendo dalle medie esatte. Non è dimostrato che ogni suo DOF abbia sempre un
significato puntuale puro. La scelta matematica generale resta aperta; prima di
trasferirla al primitivo servono campi con u,v,p variabili e verifica degli altri
regimi, oltre alla campagna MMS comparativa ancora non avviata.

## Riproduzione, verifiche e file

Comandi eseguiti dalla radice del repository (esecuzione seriale MPI; il launcher
imposta OMP_NUM_THREADS=1). Le directory v2/final identificano le build accettate:

```sh
python3 UnitTests/run_euler_dof_audit.py --output-dir /private/tmp/su2-dof-audit-v2
python3 UnitTests/export_euler_dof_audit.py --executable /private/tmp/su2-dof-audit-v2/euler_dof_audit --output-dir /private/tmp/su2-dof-audit-v2/campaign --jobs 2
python3 UnitTests/export_euler_dof_audit.py --executable /private/tmp/su2-dof-audit-v2/euler_dof_audit --output-dir /private/tmp/su2-dof-audit-v2/refinement128 --jobs 2 --sizes 128 --families jitter
python3 UnitTests/export_euler_dof_audit.py --executable /private/tmp/su2-dof-audit-v2/euler_dof_audit --output-dir /private/tmp/su2-dof-audit-v2/refinement128 --collect-only --sizes 128 --families jitter --extra-records /private/tmp/su2-dof-audit-v2/campaign/all_records.csv
python3 UnitTests/run_euler_dof_audit.py --output-dir /private/tmp/su2-dof-audit-final
python3 UnitTests/export_euler_dof_audit.py --executable /private/tmp/su2-dof-audit-final/euler_dof_audit --output-dir /private/tmp/su2-dof-audit-final/affine --jobs 2 --families jitter --benchmark affine
python3 UnitTests/run_euler_dof_audit.py --output-dir /private/tmp/su2-dof-audit-complete
python3 UnitTests/export_euler_dof_audit.py --executable /private/tmp/su2-dof-audit-complete/euler_dof_audit --output-dir /private/tmp/su2-dof-audit-complete/geometry --geometry-only --jobs 2
python3 -m py_compile UnitTests/run_euler_dof_audit.py UnitTests/export_euler_dof_audit.py
/private/tmp/su2-phase5-tests/test_driver --reporter compact
```

Per rigenerare tutto dalla sorgente finale si può usare l'eseguibile complete
anche per le prime tre campagne: aggiunge righe moment_component, senza cambiare
i risultati evolutivi. I due smoke test finali sinusoidale/affine a N=8 sono
stati ripetuti con questa build per confrontare i campi con le build precedenti.

Il build script riusa compilatore, flag e dipendenze native della build `build`;
compila esclusivamente il banco di prova e lo collega a libSU2core e Common.
Non ricompila né sostituisce metodi del solver conservativo.
La campagna salva CSV, log e snapshot `.field` nella directory temporanea;
gli script permettono di rigenerarli. I CSV di sintesi sono conservati in Docs.
`--collect-only` evita il rilancio; `--extra-records` permette di aggregare
all_records.csv di un'altra campagna dello stesso benchmark.

La versione finale aggiunge componenti geometriche e controlli d'ingresso;
i calcoli evolutivi riportati sono stati prodotti dalle successive build isolate
v2/final, con gli stessi metodi e parametri. Non è stata cambiata alcuna formula
del solver tra le campagne. La geometria finale è stata riverificata separatamente.

Suite precedente: **70 casi, 22,253,237 assert superati** sul driver già costruito
della FASE 5, appropriato perché nessun suo sorgente/dipendenza è stato modificato.
Build del nuovo eseguibile riuscita; rimangono solo warning upstream già presenti
su opzioni warning GCC ignorate da clang e librerie duplicate al link.
I primi smoke test hanno rivelato preprocessing P2P e contesto globale mancanti
nel banco di prova; sono stati corretti prima di tutte le campagne accettate.

File nuovi: UnitTests/SU2_CFD/euler_dof_audit.cpp,
UnitTests/run_euler_dof_audit.py, UnitTests/export_euler_dof_audit.py,
questo rapporto e i sette CSV collegati sopra. Tutti inseriti tramite apply_patch.
Stato finale: undici file non tracciati intenzionali, nessun file preesistente
modificato e staging vuoto sul branch originale.
Nessun test delle fasi precedenti è stato cambiato. Nessun commit o sync.
Revisione dei risultati richiesta prima di qualsiasi FASE 6B.
