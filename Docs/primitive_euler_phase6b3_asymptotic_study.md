# FASE 6B.3 — Consolidamento asintotico e difetti firmati

## Risposta scientifica e perimetro effettivamente verificato

La convergenza evolutiva osservata nella 6B.2 è robusta rispetto
all'eliminazione **diagnostica stimata** del principale errore Forward Euler,
alle risoluzioni N≤64: Richardson cambia gli ordini L2 al massimo di
0.00823. Non è dimostrata la persistenza oltre N=64. L'evoluzione a N=128
non è stata avviata dopo la stima preventiva del costo; il difetto
semidiscreto, molto meno costoso, è stato invece misurato fino a N=128
su tutte e tre le famiglie. Non si dichiara ordine globale almeno due.

Su jitter, nell'interno fisico fisso, il difetto primitivo a t=0.1 ha ordini
L2 64→128 di (1.025,1.012,1.021,1.019). Le cinque proiezioni firmate su
funzioni lisce hanno ordini fra 1.882 e 2.206. Il rapporto fra modulo del
momento e integrale del modulo dell'integrando è 0.028–0.116 a N=128.
Questo è un riscontro quantitativo di cancellazione spaziale del difetto,
compatibile con un errore evoluto che decada più rapidamente della sua
norma forte. **Non prova il meccanismo causale dell'evoluzione** e non
dimostra superconvergenza asintotica su mesh jitter generiche.

La campagna comprende:

- riproduzione di tutte le 36 coppie della 6B.2 a T=0.1, N=8,16,32,64;
- terne dt,dt/2,dt/4, stime Richardson e confronto fra estrapolazioni successive;
- 15 valutazioni di mesh semidiscrete, ciascuna a t=0,0.05,0.1, fino a N=128;
- controllo a T=0.2 su regular/jitter, N=8,16,32, tre passi ciascuno:
  ulteriori 18 coppie evolutive;
- contabilizzazione disgiunta di bordo, ring1, ring2 e deep3;
- regressioni e controlli analitici precedenti, senza adattarli ai risultati.

Non sono stati modificati sorgenti produttivi né file delle FASI 6A.1,
6B.1, 6B.2. Nessun nuovo schema, DUAL_MEAN, Hessiana, correzione conservativa,
limiter, modifica ai bordi o al driver generale. Nessuna operazione Git di
scrittura, cambio branch o sincronizzazione.

## Costo preventivo e raffinamento non eseguito

La 6B.2 misurava circa 72 s per una coppia N=64,d=4; il numero di nodi
cresce circa di quattro e quello dei passi di quattro passando a N=128.
Con crescita O(N^4), la stima preventiva era circa 19 min per famiglia
a d=4 e 38 min a d=8. Una terna d=2,4,8, necessaria per controllare anche
Richardson sul nuovo livello, richiederebbe circa 67 min CPU per famiglia.
Regular+jitter richiederebbero circa 134 min CPU, tutte e tre circa 202 min,
prima di eventuali ulteriori dimezzamenti. L'estrapolazione del costo è
una stima, non un tempo misurato a N=128. I run baseline ripetuti sotto
concorrenza hanno impiegato 69–98 s a N=64,d=4: la stima non è ottimistica
rispetto al nuovo carico.

A N=128,T=0.1, d=4 servirebbero 26216 passi e dt=3.81446445e-6;
a d=8, 52432 passi e dt dimezzato. I 16641 nodi hanno un costo di memoria
modesto: il run semidiscreto jitter N=128 ha misurato RSS massimo 88.6 MB,
footprint 60.6 MB; il baseline evolutivo regular N=64,d=2 circa 32.1 MB RSS.
Per l'evoluzione N=128 si stimano prudenzialmente 0.1–0.2 GB per processo,
oltre al runtime MPI. La memoria non è l'ostacolo principale, è il tempo.

La matrice evolutiva estesa è stata quindi fermata **prima** dell'avvio
oneroso, come richiesto. I suoi ordini 64→128 sono **non disponibili**;
non si sostituiscono con quelli semidiscreti. Alternativa minima per una
successiva sessione di calcolo dedicata: regular e jitter, N=128,d=2,4,8,
poi controllo delle frazioni temporali e delle due estrapolazioni; solo dopo
smooth. Mantenere i passi già verificati dei livelli inferiori e ripetere
il controllo se serve d=16. Non si propone di allargare le tolleranze o di
cambiare arbitrariamente T per dichiarare raggiunto il raffinamento.

## Riproduzione e implementazione dei soli diagnostici

`run_euler_asymptotic_study.py` legge i due harness esistenti e genera una
translation unit in scratch. Ogni inserimento richiede che l'ancora testuale
sia presente esattamente una volta; un cambiamento del harness richiede
revisione. Il sorgente temporaneo completo, gli hash degli ingressi e
dell'eseguibile vengono conservati nel manifest di build.

Nel harness evolutivo si aggiunge soltanto l'esportazione dei riferimenti
W*,U* ottenuti dall'implementazione analitica C++ condivisa. Per T=0.2 si
cambia il tempo finale e si raddoppia esattamente il numero di passi della
stessa mesh a T=0.1: il passo fisico resta identico, non il numero di passi.
Nel harness semidiscreto si esportano i vettori già calcolati e i contributi
restituiti da `ComputeSpatialResidual`. Nessun contributo viene sostituito.

La ricompilazione usa i flag effettivi di `CPhysicalGeometry.cpp` e gli
archivi del vero SU2, come le fasi precedenti. Si mantengono:
gamma=1.4, riferimenti unitari, la stessa MMS non affine dipendente dal
tempo, generatori regular/smooth/jitter invariati, sorgenti S_W e J(W)S_W,
NODE/AGGREGATED/NODE_POINT, MUSCL affine WLS/NONE, Roe e FE nativi.
Il primitivo conserva il bordo debole analitico; il conservativo conserva
`BC_Custom` forte a t_n e t_(n+1). Non si impone un bordo forte al primitivo.

Confrontati 36 campi finali completi e tutte le quattro norme dei relativi
record con i run originali: **massima differenza zero**, sia nei campi
che nelle norme. Tolleranza dichiarata per la riproduzione: 2e-13 assoluta
nelle unità unitarie. La trasformazione U↔W per il solo postprocessing è
quella del gas ideale, senza una seconda implementazione della MMS.

## Separazione spazio-tempo

Per ogni mesh a T=0.1 si ripete la regola originale
n_t=d ceil(T N²/0.25), d=1,2,4. Le soluzioni finali sono confrontate
**come campi**, non sottraendo norme degli errori. Con q_d corrispondente
al divisore d:

\[
\delta_1=\|q_1-q_2\|,\quad \delta_2=\|q_2-q_4\|,
\qquad p_t=\log_2(\delta_1/\delta_2).
\]

Per il primitivo si costruiscono W0_f=2W4−W2 e W0_c=2W2−W1, poi si
calcolano U(W0_f),U(W0_c). Per il conservativo si costruiscono invece
U0_f=2U4−U2, U0_c=2U2−U1 e **solo dopo** W(U0_f),W(U0_c).
Non si scambia questa operazione con l'estrapolazione diretta di W(U).
Le stime sono admissibili nei dati analizzati; restano esclusivamente
diagnostiche, non vengono avanzate dal solver.

Si registrano separatamente errore FE, errore estrapolato, correzione
\|q0_f−q4\| e differenza \|q0_f−q0_c\|. Quest'ultima è un controllo
di coerenza delle estrapolazioni, non un bound rigoroso dell'errore residuo.
Per le componenti primarie risolvibili, nelle sette regioni, si ottiene:

| Quantità L2 massima | T=0.1 | T=0.2 |
|---|---:|---:|
| Correzione Richardson / errore FE | 0.011683 | 0.012938 |
| Differenza fra estrapolazioni / errore estrapolato | 0.001005 | 0.001091 |
| Variazione assoluta dell'ordine FE→Richardson | 0.008221 | 0.010328 |

A T=0.1 i p_t L2 risolvibili sono 0.9971–1.0822, mediana 1.0124;
a T=0.2 sono 1.0042–1.1243. Non tutte le terne grossolane sono già
asintotiche in senso stretto, ma il termine temporale corretto è piccolo
e le estrapolazioni successive sono molto più vicine dell'errore spaziale
osservato. È evidenza quantitativa a favore della separazione, non una
prova rigorosa del limite dt→0. Si mantiene il criterio precedente del 5%
e la soglia di roundoff 1e-12, usata soltanto per evitare ordini spuri;
nessun valore è scartato per il suo ordine spaziale.

Gli ordini FE sotto raffinamento sono ordini **completamente discreti lungo
dt∼h²**. Gli ordini Richardson sono stime dell'errore della soluzione
semidiscreta evoluta a T, dopo il controllo temporale. Nessuno dei due è
automaticamente l'ordine del difetto di troncamento istantaneo.

## Errori e successione degli ordini evolutivi

Le tuple W sono sempre (rho,u,v,p); le tuple U nei CSV sono
(rho,rho*u,rho*v,rho*E). `difference_W` confronta W_prim e W(U_cons);
`difference_U` confronta U(W_prim) e U_cons. Non si mescolano componenti
di rappresentazioni diverse. Tutti i valori seguenti sono globali a T=0.1,
con d=4 per FE. N finale 32 significa 16→32; 64 significa 32→64.

| Famiglia | Campo | N finale | L2 FE (rho,u,v,p) | p FE | p Richardson |
|---|---|---:|---|---|---|
| regular | W_primitive | 32 | 4.310e-04, 3.115e-04, 3.225e-04, 4.411e-04 | 2.228, 2.139, 2.163, 2.099 | 2.227, 2.137, 2.161, 2.097 |
| regular | W_primitive | 64 | 1.036e-04, 6.763e-05, 6.831e-05, 1.147e-04 | 2.057, 2.203, 2.239, 1.943 | 2.057, 2.203, 2.239, 1.943 |
| regular | W_conservative | 32 | 4.200e-04, 4.568e-04, 4.656e-04, 4.463e-04 | 1.821, 1.931, 1.933, 1.558 | 1.821, 1.930, 1.932, 1.556 |
| regular | W_conservative | 64 | 1.207e-04, 1.158e-04, 1.182e-04, 1.335e-04 | 1.799, 1.980, 1.977, 1.741 | 1.798, 1.980, 1.977, 1.741 |
| regular | difference_W | 32 | 4.479e-04, 3.192e-04, 3.284e-04, 5.186e-04 | 2.125, 1.772, 1.748, 2.040 | 2.124, 1.771, 1.747, 2.039 |
| regular | difference_W | 64 | 1.103e-04, 8.643e-05, 8.874e-05, 1.313e-04 | 2.022, 1.885, 1.888, 1.982 | 2.022, 1.885, 1.887, 1.982 |
| smooth | W_primitive | 32 | 5.263e-04, 3.703e-04, 3.972e-04, 4.929e-04 | 2.215, 2.127, 2.122, 2.174 | 2.214, 2.125, 2.120, 2.172 |
| smooth | W_primitive | 64 | 1.228e-04, 8.173e-05, 8.480e-05, 1.271e-04 | 2.099, 2.180, 2.228, 1.956 | 2.099, 2.179, 2.227, 1.955 |
| smooth | W_conservative | 32 | 5.478e-04, 5.366e-04, 5.394e-04, 5.437e-04 | 1.887, 1.902, 1.908, 1.696 | 1.886, 1.901, 1.906, 1.694 |
| smooth | W_conservative | 64 | 1.537e-04, 1.366e-04, 1.384e-04, 1.573e-04 | 1.834, 1.974, 1.963, 1.789 | 1.833, 1.974, 1.962, 1.788 |
| smooth | difference_W | 32 | 5.442e-04, 4.053e-04, 4.128e-04, 5.995e-04 | 2.026, 1.749, 1.711, 2.018 | 2.025, 1.748, 1.710, 2.017 |
| smooth | difference_W | 64 | 1.392e-04, 1.103e-04, 1.133e-04, 1.524e-04 | 1.966, 1.877, 1.866, 1.976 | 1.966, 1.877, 1.865, 1.975 |
| jitter | W_primitive | 32 | 5.030e-04, 3.251e-04, 3.536e-04, 4.854e-04 | 2.166, 2.156, 2.114, 2.109 | 2.165, 2.155, 2.112, 2.108 |
| jitter | W_primitive | 64 | 1.203e-04, 7.158e-05, 7.707e-05, 1.240e-04 | 2.064, 2.183, 2.198, 1.968 | 2.064, 2.183, 2.198, 1.968 |
| jitter | W_conservative | 32 | 4.773e-04, 4.619e-04, 4.825e-04, 4.640e-04 | 1.845, 1.965, 1.934, 1.686 | 1.845, 1.964, 1.932, 1.685 |
| jitter | W_conservative | 64 | 1.328e-04, 1.171e-04, 1.220e-04, 1.371e-04 | 1.845, 1.979, 1.984, 1.759 | 1.845, 1.979, 1.984, 1.758 |
| jitter | difference_W | 32 | 4.579e-04, 3.223e-04, 3.362e-04, 5.233e-04 | 2.118, 1.803, 1.738, 2.057 | 2.117, 1.802, 1.737, 2.056 |
| jitter | difference_W | 64 | 1.115e-04, 8.759e-05, 8.973e-05, 1.317e-04 | 2.038, 1.879, 1.906, 1.990 | 2.038, 1.879, 1.905, 1.990 |

Gli errori assoluti estrapolati, tutte le quattro norme, le componenti
conservative, N=8,16 e tutte le regioni sono nel CSV `evolution`.
Per 64→128 non esiste un dato evolutivo. La pressione primitiva resta
sotto due nell'ultima transizione e il conservativo non mostra due per
tutte le primitive. Togliere il principale errore temporale non risolve
questi limiti. Nemmeno ordini maggiori di due su alcuni intervalli sono
una dimostrazione di ordine asintotico superiore.

## Difetto semidiscreto e proiezioni firmate

Per il primitivo si misura tau_W=L_h(I_h W*,t)−W*_t; per il conservativo
tau_U=L_h(I_h U*,t)−J(W*)W*_t **soltanto sui nodi liberi**.
Si esportano anche J(W*)tau_W e la differenza con tau_U. I nodi Dirichlet
conservativi non ricevono un difetto artificiale zero: nel file firmato
la colonna conservativa è vuota su quei nodi. Il difetto debole primitivo
di bordo resta misurato.

Le norme riportano L1_V=ΣV|e|/ΣV, L2_V=(ΣVe²/ΣV)^(1/2), L∞ e RMS
nodale non pesata, con le maschere precedenti. `all_free` indica che
il bordo forte è escluso; non va confuso con `all` primitivo.

Difetto primitivo nella regione fisica fissa [0.25,0.75]² a t=0.1:

| Famiglia | N | L2 tau_W, fixed | p L2 |
|---|---:|---|---|
| regular | 32 | 4.387e-03, 2.180e-03, 2.184e-03, 6.347e-03 | 1.753, 2.515, 2.540, 1.870 |
| regular | 64 | 1.158e-03, 4.586e-04, 4.561e-04, 1.601e-03 | 1.921, 2.249, 2.260, 1.987 |
| regular | 128 | 2.951e-04, 1.089e-04, 1.084e-04, 3.986e-04 | 1.973, 2.074, 2.073, 2.006 |
| smooth | 32 | 5.322e-03, 3.099e-03, 3.120e-03, 7.505e-03 | 1.897, 2.529, 2.531, 1.914 |
| smooth | 64 | 1.361e-03, 6.040e-04, 6.042e-04, 1.859e-03 | 1.967, 2.359, 2.369, 2.013 |
| smooth | 128 | 3.439e-04, 1.348e-04, 1.347e-04, 4.605e-04 | 1.985, 2.164, 2.166, 2.013 |
| jitter | 32 | 1.423e-02, 7.696e-03, 7.171e-03, 1.821e-02 | 1.257, 1.414, 1.481, 1.263 |
| jitter | 64 | 6.970e-03, 3.670e-03, 3.362e-03, 8.911e-03 | 1.030, 1.069, 1.093, 1.031 |
| jitter | 128 | 3.425e-03, 1.819e-03, 1.656e-03, 4.396e-03 | 1.025, 1.012, 1.021, 1.019 |

Per i momenti la normalizzazione è **sempre l'area fisica 0.25**, non la
somma variabile dei volumi dei nodi selezionati:

\[
M_\varphi=\frac1{0.25}\sum_{i\in fixed}V_i\varphi_i\tau_i,
\qquad A_\varphi=\frac1{0.25}\sum_{i\in fixed}V_i|\varphi_i\tau_i|.
\]

Il CSV conserva M firmato, |M|, A e |M|/A, oltre alla misura discreta
campionata. La maschera fisica resta invariata; la sua approssimazione
nodale, in particolare sui margini della finestra, varia con la mesh.
Questo è un limite della diagnostica di momento, non un motivo per
riaggiustare la normalizzazione o filtrare i risultati.

Jitter, t=0.1, N=128, difetto primitivo:

| phi | M_phi (rho,u,v,p) | p abs(M), 64→128 | abs(M)/integrale assoluto |
|---|---|---|---|
| 1 | -7.816e-05, -4.262e-05, -3.748e-05, -1.239e-04 | 1.938, 1.953, 2.038, 2.023 | 0.029, 0.029, 0.028, 0.036 |
| x | -4.746e-05, -2.143e-05, -2.092e-05, -6.177e-05 | 2.013, 2.068, 2.206, 2.084 | 0.034, 0.029, 0.031, 0.036 |
| y | -4.666e-05, -2.364e-05, -1.883e-05, -6.291e-05 | 2.042, 2.166, 2.176, 2.108 | 0.035, 0.031, 0.028, 0.036 |
| sin_sin | 1.025e-04, 3.506e-05, 3.070e-05, 1.549e-04 | 1.923, 2.012, 1.976, 2.006 | 0.094, 0.062, 0.060, 0.116 |
| cos_cos | -5.519e-05, -3.050e-05, -3.461e-05, -8.086e-05 | 1.882, 1.993, 1.990, 1.998 | 0.054, 0.054, 0.068, 0.061 |

Anche nelle componenti conservative coerenti gli ordini dei momenti
64→128 sono 1.882–2.138 per Jtau_W e 1.880–2.122 per tau_U.
Il CSV contiene tutti i segni, tempi, livelli e funzioni test, comprese
eventuali inversioni di segno o ordini negativi. Gli ordini di |M| non
vengono stimati quando il segnale è sotto la soglia dichiarata.

Il file compresso `signed_defects.csv.gz` contiene coordinate, volume,
strato topologico, flag libero, famiglia/N/tempo/componente e difetti
firmati. È stato costruito dai campi nodali, **non dalle norme**.
Gli hash e il percorso scratch sono nel manifest `provenance`; un momento
è stato verificato anche con una somma scalare indipendente del campo
grezzo (jitter N=128, t=0.1, phi=x, rho:
−4.746406498750841e-5).

### Contributi produttivi separati

Per il primitivo si esportano faces=−F/V, volume=−B/V,
boundary=−G/V, source=S e rhs con i loro segni effettivi. La somma
ricostruisce il RHS produttivo, controllata a ogni nodo/componente con
tolleranza 5e-13; massimo misurato 2.23e-16, registrato nei `checks`. Le norme dei singoli
contributi non sono errori di troncamento: i termini possono essere
O(1) e cancellarsi nella somma.

Il conservativo non espone la medesima decomposizione in fluttuazioni
e volume. Si conserva quindi la misura completa del difetto libero e
i controlli nativi di sorgente della 6B.1, senza costruire un nuovo
assemblaggio o assegnare un flusso fisico al residuale Dirichlet azzerato.

## Bordo, anelli e peso globale

Per ciascuna regione si registrano norma locale e
C_R=Σ_(i in R) V_i e_i². La partizione disgiunta è esclusivamente
**boundary + ring1 + ring2 + deep3**. `all`, `interior`, `fixed` sono
diagnostiche sovrapposte, mai sommate nella partizione. Si verifica
C_boundary+C_ring1+C_ring2+C_deep3=C_all entro 5e-13 relativo;
lo stesso vale per i difetti sui nodi liberi, senza il bordo escluso.
Il CSV pubblica misura della regione, C_R, C_all e frazione per componente.

A T=0.1,N=64, ordini locali 32→64 e pesi della soluzione primitiva
estrapolata (le tuple C hanno unità di errore al quadrato per area):

| Famiglia | Regione | p difetto tau_W | p soluzione Richardson | C_R (rho,u,v,p) | C_R/C_all |
|---|---|---|---|---|---|
| regular | boundary | 1.120, 1.015, 1.020, 1.051 | 2.073, 1.845, 1.962, 2.115 | 2.513e-09, 6.836e-10, 6.785e-10, 1.978e-09 | 0.234, 0.149, 0.145, 0.150 |
| regular | ring1 | 0.976, 1.094, 1.113, 0.965 | 1.957, 2.216, 1.865, 2.065 | 5.694e-10, 9.187e-11, 1.161e-10, 5.645e-10 | 0.053, 0.020, 0.025, 0.043 |
| regular | ring2 | 2.092, 2.318, 2.228, 2.178 | 1.758, 2.553, 2.650, 1.691 | 3.802e-10, 1.675e-10, 1.300e-10, 3.094e-10 | 0.035, 0.037, 0.028, 0.024 |
| regular | deep3 | 1.951, 2.249, 2.253, 2.032 | 1.789, 2.207, 2.251, 1.735 | 7.257e-09, 3.633e-09, 3.743e-09, 1.030e-08 | 0.677, 0.794, 0.802, 0.783 |
| jitter | boundary | 1.117, 1.020, 1.019, 1.052 | 2.079, 1.818, 1.974, 2.115 | 2.717e-09, 7.245e-10, 7.336e-10, 2.180e-09 | 0.188, 0.141, 0.123, 0.142 |
| jitter | ring1 | 0.973, 1.084, 1.057, 0.950 | 2.021, 2.148, 1.940, 2.048 | 8.108e-10, 1.491e-10, 1.912e-10, 7.846e-10 | 0.056, 0.029, 0.032, 0.051 |
| jitter | ring2 | 1.337, 0.984, 1.006, 1.245 | 1.931, 2.421, 2.487, 1.893 | 6.321e-10, 2.077e-10, 2.069e-10, 4.946e-10 | 0.044, 0.041, 0.035, 0.032 |
| jitter | deep3 | 1.060, 1.068, 1.092, 1.063 | 1.859, 2.195, 2.199, 1.767 | 1.030e-08, 4.044e-09, 4.809e-09, 1.192e-08 | 0.712, 0.789, 0.810, 0.775 |

Il primo anello ha difetto circa O(h), ma errore evoluto con ordini
1.865–2.216 su regular e 1.940–2.148 su jitter e pesa solo circa
2–6% di C_all. Non è sostenuta l'ipotesi che, da solo, imponga l'ordine
globale evoluto. Il bordo debole pesa circa 12–23%, l'interno deep3
68–81% in queste due famiglie. Sono pesi osservati a N=64, non un
argomento asintotico sufficiente né una prova di indipendenza dal bordo.

Per il conservativo il bordo imposto è esatto per costruzione: nessun
ordine di bordo viene dichiarato. Ring1 pesa circa 5–13%, deep3
74–87%; la possibile influenza del bordo sulle caratteristiche entranti
non si misura dal solo errore nullo dei DOF prescritti. Il primitivo
non viene modificato per uguagliare questo trattamento.

## Controllo mirato a T=0.2

Regular e jitter, N=8,16,32, d=1,2,4. I passi sono esattamente quelli
della stessa mesh a T=0.1, con il doppio degli avanzamenti. Le sorgenti
e i bordi sono valutati al tempo effettivo lungo l'intera traiettoria.
La tabella include 8→16 e 16→32; non sostituisce il raffinamento N=128
richiesto sul tempo originale.

| Famiglia | Campo | N finale | L2 FE (rho,u,v,p) | p FE | p Richardson |
|---|---|---:|---|---|---|
| regular | W_primitive | 16 | 2.904e-03, 1.619e-03, 1.713e-03, 2.234e-03 | 2.329, 1.698, 1.727, 2.416 | 2.327, 1.699, 1.726, 2.414 |
| regular | W_primitive | 32 | 5.604e-04, 4.398e-04, 4.526e-04, 4.523e-04 | 2.374, 1.880, 1.920, 2.304 | 2.373, 1.879, 1.919, 2.303 |
| regular | W_conservative | 16 | 2.956e-03, 2.291e-03, 2.361e-03, 2.385e-03 | 2.257, 1.474, 1.477, 2.357 | 2.254, 1.472, 1.475, 2.352 |
| regular | W_conservative | 32 | 6.051e-04, 6.601e-04, 6.755e-04, 4.667e-04 | 2.288, 1.795, 1.805, 2.354 | 2.287, 1.794, 1.804, 2.352 |
| regular | difference_W | 16 | 2.764e-03, 1.768e-03, 1.796e-03, 3.062e-03 | 2.069, 1.686, 1.651, 1.980 | 2.066, 1.684, 1.649, 1.977 |
| regular | difference_W | 32 | 6.414e-04, 4.913e-04, 4.950e-04, 7.539e-04 | 2.108, 1.848, 1.859, 2.022 | 2.107, 1.847, 1.858, 2.021 |
| jitter | W_primitive | 16 | 3.151e-03, 1.745e-03, 1.832e-03, 2.412e-03 | 2.292, 1.844, 1.744, 2.347 | 2.289, 1.844, 1.743, 2.344 |
| jitter | W_primitive | 32 | 6.235e-04, 4.521e-04, 4.746e-04, 5.045e-04 | 2.337, 1.948, 1.949, 2.257 | 2.337, 1.947, 1.948, 2.257 |
| jitter | W_conservative | 16 | 3.137e-03, 2.385e-03, 2.416e-03, 2.453e-03 | 2.202, 1.497, 1.567, 2.259 | 2.199, 1.495, 1.565, 2.255 |
| jitter | W_conservative | 32 | 6.493e-04, 6.575e-04, 6.816e-04, 4.795e-04 | 2.272, 1.859, 1.825, 2.355 | 2.271, 1.858, 1.824, 2.354 |
| jitter | difference_W | 16 | 2.799e-03, 1.816e-03, 1.813e-03, 3.084e-03 | 2.005, 1.559, 1.611, 1.946 | 2.002, 1.557, 1.609, 1.943 |
| jitter | difference_W | 32 | 6.537e-04, 4.935e-04, 5.015e-04, 7.547e-04 | 2.098, 1.880, 1.854, 2.031 | 2.098, 1.879, 1.853, 2.030 |

La differenza reciproca diminuisce e il controllo temporale resta piccolo,
ma gli ordini dipendono da componente e intervallo. Non si ricava una
garanzia generale dall'aumento del tempo finale su questo sottoinsieme.

## Ammissibilità, bilanci e test

A T=0.1 i minimi su tutte le traiettorie sono quelli riprodotti della
6B.2: rho_prim≥0.8994332, p_prim≥0.9658616,
rho_cons≥0.8994762, p_cons≥0.9826444.
A T=0.2: rho_prim≥0.8994332, p_prim≥0.9607912,
rho_cons≥0.8994762, p_cons≥0.9705509.
In tutti i run theta=1, nessuna contrazione e nessun fallback nativo;
il candidato FE grezzo viene controllato prima del preprocessing successivo.

I tempi di stadio differiscono da quelli richiesti al massimo di 2.78e-17.
I bilanci discreti diagnostici chiudono entro 4.89e-15 a T=0.1 e
2.23e-15 a T=0.2. Per il conservativo includono esplicitamente
l'iniezione del vincolo forte; per il primitivo separano incremento
J(W)rhs e resto nonlineare della trasformazione U(W).
Sono registri dell'aggiornamento, **non** dimostrazioni di conservazione
fisica, corrette velocità di shock o condizioni Rankine–Hugoniot.

Verifiche eseguite:

- regressioni esistenti: 70 casi, 22 253 237 asserzioni superate;
- analitica MMS esistente: 1 caso, 3300 asserzioni superate;
- 3 nuovi test di postprocessing: cancellazione del termine temporale
  dominante con Richardson, non commutatività conversione/estrapolazione,
  additività delle energie regionali su volumi non uniformi;
- 36 riproduzioni di campi/norme, controlli di stadio/sorgenti/ammissibilità,
  ricostruzione dei contributi, partizioni senza sovrapposizioni;
- verifica Python, whitespace, identità dei CSV archiviati e hash dei
  precedenti 28 file delle fasi 6 preservati.

## Conclusione e minimo esperimento discriminante successivo

**Osservato:** l'errore temporale dominante è piccolo; Richardson non
cambia la conclusione evolutiva a N≤64. Il difetto forte jitter resta
di primo ordine fino a N=128, mentre i suoi momenti lisci decadono
quasi come h². I primi anelli non dominano l'energia dell'errore finale.

**Ipotesi compatibile:** una parte importante del difetto jitter è
oscillante e ha piccole proiezioni su modi lisci; l'evoluzione numerica
può attenuarla o trasportarla in modo da produrre un errore finale
più piccolo della stima basata sulla sola norma istantanea del difetto.

**Non dimostrato:** che questa sia la causa sufficiente, che WLS o
quadratura siano il responsabile isolato, che l'ordine evoluto persista
a N=128 o su altre famiglie irregolari. Momenti piccoli non controllano
da soli il propagatore discreto, la dipendenza temporale del difetto o
gli effetti delle condizioni al contorno. L'analisi non include una
stima di stabilità uniforme del semigruppo discreto.

Il minimo seguito per la robustezza asintotica è la terna evolutiva
N=128 descritta nel piano di costo. Per discriminare la spiegazione
causale occorre poi un **esperimento tangente diagnostico**, da approvare
separatamente: lungo l'analitica campionata stimare azioni di L'_h tramite
differenze direzionali dell'operatore produttivo e ricostruire la risposta
linearizzata alla storia firmata tau(t), confrontandola con l'errore reale.
Registrare separatamente risposta alle proiezioni lisce e al resto, usando
le stesse regioni e condizioni originali. Questo distinguerebbe semplice
cancellazione dei momenti da attenuazione dinamica o cancellazione nel
tempo. Non basta mescolare casualmente il difetto preservandone la norma,
perché cambierebbero anche le proiezioni e l'interazione col bordo.
Nessuno di questi nuovi esperimenti è stato implementato qui.

## File, comandi e stato Git

Nuovi strumenti:

- `UnitTests/run_euler_asymptotic_study.py`: compilazione e sole inserzioni diagnostiche;
- `UnitTests/run_euler_asymptotic_campaign.py`: matrice esplicita limitata nei costi;
- `UnitTests/export_euler_asymptotic_study.py`: Richardson, momenti firmati e partizioni;
- `UnitTests/test_euler_asymptotic_postprocessing.py`: verifiche indipendenti.

I nove artefatti `Docs/primitive_euler_phase6b3_*` oltre al rapporto sono:
`evolution.csv`, `temporal.csv`, `defects.csv`, `moments.csv`,
`contributions.csv`, `rings.csv`, `checks.csv`, `reproduction.csv`,
`provenance.json`. I CSV sono pubblicati per tutte le componenti,
regioni, famiglie e livelli effettivamente eseguiti, non solo per i casi
che mostrano un ordine desiderato. I campi firmati compressi sono circa
28 MB: restano nello scratch indicato nel manifest, insieme ai campi
evolutivi, ai sorgenti generati e ai log. Lo scratch non è permanente;
si rigenera con i comandi sotto. Non si presenta il solo manifest come
un archivio permanente dei campi nodali.

```sh
python3 UnitTests/run_euler_asymptotic_study.py defect --output-dir /private/tmp/su2-6b3/defect-build
python3 UnitTests/run_euler_asymptotic_study.py evolution --output-dir /private/tmp/su2-6b3/evolution-build
python3 UnitTests/run_euler_asymptotic_study.py evolution --time 0.2 --output-dir /private/tmp/su2-6b3/long-build
python3 UnitTests/run_euler_asymptotic_campaign.py defect --executable /private/tmp/su2-6b3/defect-build/euler_asymptotic_defect --output-dir /private/tmp/su2-6b3/defects --jobs 2
python3 UnitTests/run_euler_asymptotic_campaign.py baseline --executable /private/tmp/su2-6b3/evolution-build/euler_asymptotic_evolution --output-dir /private/tmp/su2-6b3/baseline --jobs 2
python3 UnitTests/run_euler_asymptotic_campaign.py long --executable /private/tmp/su2-6b3/long-build/euler_asymptotic_evolution --output-dir /private/tmp/su2-6b3/long --jobs 1
python3 UnitTests/export_euler_asymptotic_study.py --root /private/tmp/su2-6b3 --output-dir /private/tmp/su2-6b3/analysis --reference-baseline /private/tmp/su2-mms-6b2-complete/campaign --preservation-manifest /private/tmp/su2-6b3-preserved.json
python3 UnitTests/test_euler_asymptotic_postprocessing.py
/private/tmp/su2-phase5-tests/test_driver --reporter compact
/private/tmp/su2-mms-6b1-final/euler_mms_analytic --reporter compact
python3 -m py_compile UnitTests/run_euler_asymptotic_study.py UnitTests/run_euler_asymptotic_campaign.py UnitTests/export_euler_asymptotic_study.py UnitTests/test_euler_asymptotic_postprocessing.py
git diff --check
git status --short
```

Se la baseline scratch non esiste più, rigenerarla prima con i comandi
6B.2; il confronto completo dei campi non viene silenziosamente saltato.
Le valutazioni dei difetti e delle evoluzioni sono state effettuate con
i comandi sopra; l'esportazione è stata eseguita inizialmente per i
difetti e poi con `--skip-defects` per l'evoluzione, evitando di ricomprimere
gli stessi campi. Il comando completo senza tale flag rigenera tutto.
Lo script impacchetta anche i CSV aggregati del rapporto: concatena i due tempi e riporta
nei checks gli scalari semidiscreti sotto la chiave stage=defect.
Il manifest di preservazione è una mappa percorso→SHA256 acquisita prima
dei run; i suoi valori sono conservati nel `provenance.json` del rapporto.
Per rigenerazioni future può essere ricostruito da quella chiave, oppure
omesso se si vogliono rigenerare soltanto i diagnostici senza ripetere
la verifica storica di preservazione.

Stato Git: branch `feature/primitive-path-solver`, HEAD invariato
`5be7a6e72a6c082cce192e73e8dd53eb17f1b13b`.
L'unica modifica tracciata resta la riga preesistente in `UnitTests/meson.build`;
i file 6A.1/6B.1/6B.2 restano non tracciati e preservati. Questa fase
aggiunge solo i quattro strumenti, il rapporto e i nove artefatti elencati.
Non sono stati eseguiti commit, push, pull, fetch, sync, merge o cambi branch.
Ci si ferma a questo rapporto, con il limite evolutivo N=128 esplicito.
