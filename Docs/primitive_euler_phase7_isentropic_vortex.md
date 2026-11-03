# FASE 7 — Trasporto del vortice isentropico 2D

## Esito e limiti

Eseguiti il vero `CPrimitiveEulerSolver` e il vero `CEulerSolver` SU2,
senza sorgenti, su regular/jitter N=8,16,32,64. Sono state completate
32 coppie di traiettorie con quattro passi temporali per mesh, dopo un
pilota regular N=16. Non è stata eseguita smooth né alcuna mesh N≥128.

Entrambi i solver convergono all'analitica; la loro differenza converge
anch'essa a zero. Nell'ultima transizione gli errori rispetto all'analitica
sono compatibili con il secondo ordine dopo controllo temporale. La
differenza fra solver su jitter **non** mostra due in tutte le componenti:
gli ordini L2 estrapolati sono (1.881,1.715,1.796,2.250).
Non si interpreta questo dato come un errore tecnico da correggere.

Sono disponibili i tre confronti indipendenti, tutte le norme e regioni,
grafici, diagnostiche fisiche e 24 VTU in otto sequenze PVD. ParaView 5.13.2
ha letto effettivamente tutti i file con i suoi lettori VTK/PVD. Non sono
stati eseguiti un rendering o un'apertura grafica interattiva.

Nessun file preesistente, operatore produttivo o test precedente è stato
modificato. Branch `feature/primitive-path-solver`, HEAD invariato rispetto
all'inizio della fase: `5c60ffa9897f10360ee6f35907630860928f44d6`.
Nessun commit, push, pull, fetch, sync, merge, rebase o cambio branch.

## Soluzione analitica omogenea

Gamma=1.4, beta=2, rho_inf=p_inf=u_inf=v_inf=1; centro iniziale
(−0.25,−0.25). Si usa esattamente il profilo richiesto:

\[
X=x+0.25-t,\quad Y=y+0.25-t,\quad
f=e^{(1-X^2-Y^2)/2},\quad c=\frac{\beta}{2\pi},
\]
\[
u=1-cYf,\quad v=1+cXf,\quad
\theta=1-\frac{(\gamma-1)\beta^2}{8\gamma\pi^2}f^2,
\quad \rho=\theta^{5/2},\quad p=\theta^{7/2}.
\]

Lo stato conservativo è (rho,rho*u,rho*v,p/(gamma−1)+rho*(u²+v²)/2).
Le derivate sono differenziate esplicitamente nel nuovo header analitico:
f_x=−Xf, f_y=−Yf, theta_x=2aXf², theta_y=2aYf²,
u_x=cXYf, u_y=−c(1−Y²)f, v_x=c(1−X²)f, v_y=−cXYf;
le derivate termodinamiche seguono dalle potenze di theta e W_t=−W_x−W_y.
Non si ottengono le derivate dalla soluzione numerica.

Il test standalone verifica l'espansione primitiva, il kernel matriciale
produttivo, l'espansione conservativa indipendente del prodotto e un controllo
indipendente mediante differenze centrate di W e dei flussi conservativi.
Riutilizza le trasformazioni fisiche della FASE 6; non ne duplica l'operatore.

| Verifica | Massimo osservato |
|---|---:|
| W_t+A_x W_x+A_y W_y | 6.94e-17 |
| U_t+F_x+G_y, derivate analitiche | 1.95e-16 |
| Differenze centrate delle derivate W, eps=1e-5 | 5.55e-11 |
| PDE conservativa con flussi differenziati numericamente | 2.23e-10 |
| p/rho^gamma−1 | 2.23e-16 |
| Conversione W→U→W | 2.23e-16 |

Passano 890 controlli, inclusi temperatura p/rho>0, rho/p>0,
traslazione del profilo, velocità (1,1) nel centro e minimo di pressione
nel centro continuo. Le soglie analitiche sono dichiarate nel test prima
del confronto con i risultati dei solver. Nessuna sorgente MMS viene usata.

## Dominio, risoluzione e trattamento del bordo

Dominio [−5,5]², tempo finale 0.5. Il centro si sposta da
(−0.25,−0.25) a (0.25,0.25); distanza minima dal bordo 4.75.
Poiché r exp((1−r²)/2) decresce per r≥1, si ricava un envelope valido
per **tutto** l'intervallo temporale:

| Perturbazione analitica al bordo | Limite superiore |
|---|---:|
| Norma della perturbazione di velocità | 3.14273e-5 |
| 1−rho | 1.56340e-11 |
| 1−p | 2.18875e-11 |

Il bordo usa comunque lo stato analitico esatto dipendente dal tempo,
non lo stato uniforme. Nel conservativo si mantiene `BC_Custom` forte
a t_n, dopo l'assemblaggio e a t_(n+1). Nel primitivo resta il trattamento
debole path-conservative sui punti delle mezze facce. Non sono condizioni
numericamente equivalenti a mesh finita.

La geometria è la `Grid` originale delle FASI 6, scalata con x→10x−5;
connettività, identificazione degli strati e costruzione duale sono invariati.
h=10/N; triangoli 2N² e nodi (N+1)². Per jitter, solo gli interni sono
perturbati di 0.18h sin(2.1i+1.3j) in x e 0.16h cos(1.7i−0.9j) in y.
**Non esiste un seed casuale**: le fasi delle funzioni trigonometriche
sono fissate, esattamente come nella FASE 6. Non si introduce un generatore
casuale nuovo. Smooth resta disponibile nel supporto precedente invariato,
ma il nuovo harness ne rifiuta l'esecuzione.

Il diametro del nucleo definito da r≤1 attraversa circa 1.6,3.2,6.4,12.8
spaziature, rispettivamente: N=8 è deliberatamente conservato anche se
sottorisolve il nucleo. Questo limita l'interpretazione asintotica delle
prime transizioni; non si elimina la mesh perché non mostra ordine due.

## Configurazioni native e pilota

I due solver ricevono la stessa mesh, gamma, riferimenti unitari e stato
fisico nodale analitico iniziale. Rimangono Forward Euler nativo, Roe,
MUSCL affine, WLS, limiter NONE; per il primitivo NODE/AGGREGATED/NODE_POINT.
Non si modificano ricostruzione, quadratura, kernel o fluttuazioni.

`run_euler_vortex.py` riutilizza il harness accoppiato della 6B.2 mediante
inserzioni controllate in una translation unit temporanea, conservata con
i suoi hash. Il nuovo adapter restituisce `IsManufacturedSolution=false`;
l'eventuale chiamata a sorgenti MMS genera un errore. `AdvanceExplicitEuler`
primitivo viene chiamato senza sorgente. Il percorso nativo `Source_Residual`
non ha sorgenti attive; entrambi i conteggi sono esattamente zero in tutti
i passi. Nessun solver conservativo sostitutivo è stato implementato.

Il pilota regular N=16 usa 64 passi, dt=0.0078125. Sono verificati
inizializzazione entro 2e-14 assoluto, tempi, candidati FE, positività,
assenza di fallback, bordo e writer. Un controllo preliminare che richiedeva
identità bit per bit dell'inizializzazione ha rilevato solo 2.22e-16 di
roundoff; si è usata la tolleranza esplicita 2e-14 già adottata dal harness
per confrontare rappresentazioni fisiche, senza modificare alcun test
precedente o risultato numerico. Il pilota automatico supera tutti i controlli
e salva un manifest; la campagna rifiuta di partire senza un pilota superato
con lo stesso eseguibile.

## Controllo temporale quantitativo

Si parte da dt≈0.02h², con numero di passi
n_t=d ceil(0.5 N²/(0.02*10²))=d N²/4. I primi divisori sono 1,2,4;
quattro controlli L2 del ring1 a N=8,t=0.25 superano di poco il criterio
del 5% (massimo 5.1155%). Questi dati restano integralmente nel CSV
`temporal_initial`, non vengono cancellati o adattati.

Si aggiunge d=8 a **entrambi** i solver e a tutte le mesh. I risultati
FE principali usano d=8; la terna finale di controllo è d=2,4,8.
Nessun passo locale o diverso fra solver. A N=64: 8192 passi,
dt=6.103515625e-5; a N=8: 128 passi, dt=0.00390625.

Le differenze sono calcolate sui campi finali agli stessi nodi:
delta1=||q2−q4||, delta2=||q4−q8||, p_t=log2(delta1/delta2).
Richardson estrapola W nel primitivo e **U nel conservativo**, prima
della conversione a W. Si confrontano anche le estrapolazioni dalle
coppie (2,4) e (4,8). Sono esclusivamente diagnostici esterni ai solver.

Passano tutti gli 832 controlli L2 risolvibili a t=0.25 e t=0.5,
su componenti e regioni dei due solver: massima correzione Richardson
rispetto all'errore FE=2.5153%. Gli ordini temporali L2 risolvibili sono
0.9785–1.0884, compatibili con FE al primo ordine. Il bordo imposto del
conservativo e i segnali sotto 1e-12 non ricevono un ordine artificiale.
I record non risolvibili restano nei CSV con ordine vuoto.

La massima variazione degli ordini L2 FE→Richardson è 0.01858;
fra gli ultimi due passi FE è 0.01850, per gli errori primari risolvibili
nelle regioni considerate. La correzione è quindi piccola entro il criterio
dichiarato; non è una dimostrazione del limite esatto dt→0. I CSV contengono
anche L1, L∞ e RMS, comprese eventuali sensibilità maggiori in quelle norme.

Gli ordini FE sono **completamente discreti lungo dt∼h²**. Gli ordini
Richardson stimano il comportamento spaziale dell'errore evoluto dopo
controllo temporale; non si chiamano automaticamente «ordine spaziale puro».

## Tre confronti indipendenti: norme e ordini

Si misurano E_P=||W_P−W*||, E_C=||W(U_C)−W*|| e
E_PC=||W_P−W(U_C)|| sulle stesse coordinate e maschere. Le norme sono:
L1_V=ΣV|e|/ΣV, L2_V=(ΣVe²/ΣV)^(1/2), L∞=max|e|,
RMS=(Σe²/n)^(1/2) non pesata. Le tuple sono sempre (rho,u,v,p).
Gli ordini sono log2(E_N/E_2N). La conversione primitiva nativa esportata
è usata direttamente per le norme FE, compreso il suo roundoff: writer
e norme non utilizzano due campi quasi uguali ma distinti.

Risultati globali a T=0.5 (FE d=8 e stima Richardson):

| Mesh | Errore | Transizione | L2 FE finale (rho,u,v,p) | p FE | p Richardson |
|---|---|---|---|---|---|
| regular | W_primitive | 8→16 | 2.230e-03, 9.122e-03, 8.370e-03, 3.067e-03 | 1.336, 1.481, 1.622, 1.342 | 1.332, 1.479, 1.619, 1.338 |
| regular | W_primitive | 16→32 | 5.694e-04, 2.240e-03, 1.997e-03, 7.795e-04 | 1.970, 2.026, 2.068, 1.976 | 1.968, 2.025, 2.066, 1.975 |
| regular | W_primitive | 32→64 | 1.222e-04, 4.847e-04, 4.380e-04, 1.668e-04 | 2.220, 2.209, 2.189, 2.225 | 2.220, 2.210, 2.190, 2.224 |
| regular | W_conservative | 8→16 | 2.242e-03, 9.022e-03, 8.333e-03, 3.352e-03 | 1.326, 1.493, 1.610, 1.305 | 1.322, 1.491, 1.607, 1.301 |
| regular | W_conservative | 16→32 | 5.628e-04, 2.219e-03, 1.981e-03, 8.309e-04 | 1.994, 2.024, 2.073, 2.012 | 1.993, 2.023, 2.071, 2.011 |
| regular | W_conservative | 32→64 | 1.204e-04, 4.824e-04, 4.360e-04, 1.734e-04 | 2.225, 2.202, 2.184, 2.260 | 2.224, 2.203, 2.185, 2.260 |
| regular | difference_W | 8→16 | 1.097e-04, 2.176e-04, 2.547e-04, 6.610e-04 | 1.723, 2.253, 2.293, 1.213 | 1.720, 2.240, 2.281, 1.208 |
| regular | difference_W | 16→32 | 4.438e-05, 4.385e-05, 7.770e-05, 1.376e-04 | 1.305, 2.311, 1.713, 2.265 | 1.301, 2.308, 1.708, 2.259 |
| regular | difference_W | 32→64 | 9.383e-06, 7.083e-06, 1.453e-05, 1.961e-05 | 2.242, 2.630, 2.419, 2.811 | 2.237, 2.624, 2.411, 2.801 |
| jitter | W_primitive | 8→16 | 2.254e-03, 8.537e-03, 6.835e-03, 3.096e-03 | 1.446, 1.754, 1.838, 1.461 | 1.442, 1.751, 1.834, 1.458 |
| jitter | W_primitive | 16→32 | 6.564e-04, 2.219e-03, 2.142e-03, 9.009e-04 | 1.780, 1.944, 1.674, 1.781 | 1.779, 1.943, 1.675, 1.780 |
| jitter | W_primitive | 32→64 | 1.627e-04, 4.800e-04, 4.828e-04, 2.234e-04 | 2.012, 2.209, 2.149, 2.012 | 2.012, 2.210, 2.150, 2.012 |
| jitter | W_conservative | 8→16 | 2.259e-03, 8.497e-03, 6.809e-03, 3.476e-03 | 1.440, 1.759, 1.827, 1.411 | 1.436, 1.756, 1.823, 1.407 |
| jitter | W_conservative | 16→32 | 6.572e-04, 2.204e-03, 2.131e-03, 9.587e-04 | 1.781, 1.947, 1.676, 1.858 | 1.781, 1.946, 1.677, 1.858 |
| jitter | W_conservative | 32→64 | 1.631e-04, 4.797e-04, 4.814e-04, 2.328e-04 | 2.010, 2.200, 2.146, 2.042 | 2.011, 2.201, 2.147, 2.043 |
| jitter | difference_W | 8→16 | 1.822e-04, 3.118e-04, 3.408e-04, 7.748e-04 | 1.218, 1.998, 1.855, 1.064 | 1.218, 1.988, 1.845, 1.061 |
| jitter | difference_W | 16→32 | 6.399e-05, 1.091e-04, 9.716e-05, 1.637e-04 | 1.510, 1.515, 1.810, 2.243 | 1.507, 1.515, 1.808, 2.240 |
| jitter | difference_W | 32→64 | 1.737e-05, 3.322e-05, 2.801e-05, 3.444e-05 | 1.881, 1.715, 1.794, 2.249 | 1.881, 1.715, 1.796, 2.250 |

Confronti nella regione mobile r≤2 centrata sull'analitica allo stesso tempo:

| Mesh | Errore | Transizione | L2 FE finale (rho,u,v,p) | p FE | p Richardson |
|---|---|---|---|---|---|
| regular | W_primitive | 8→16 | 6.082e-03, 2.396e-02, 2.185e-02, 8.353e-03 | 1.265, 1.506, 1.658, 1.270 | 1.262, 1.504, 1.656, 1.266 |
| regular | W_primitive | 16→32 | 1.566e-03, 5.980e-03, 5.243e-03, 2.141e-03 | 1.958, 2.002, 2.059, 1.964 | 1.956, 2.001, 2.058, 1.962 |
| regular | W_primitive | 32→64 | 3.370e-04, 1.322e-03, 1.177e-03, 4.595e-04 | 2.216, 2.178, 2.156, 2.220 | 2.215, 2.179, 2.157, 2.220 |
| regular | W_conservative | 8→16 | 6.117e-03, 2.369e-02, 2.170e-02, 9.195e-03 | 1.253, 1.521, 1.650, 1.257 | 1.249, 1.519, 1.648, 1.254 |
| regular | W_conservative | 16→32 | 1.545e-03, 5.920e-03, 5.183e-03, 2.285e-03 | 1.985, 2.000, 2.065, 2.009 | 1.984, 2.000, 2.064, 2.007 |
| regular | W_conservative | 32→64 | 3.315e-04, 1.315e-03, 1.170e-03, 4.777e-04 | 2.220, 2.170, 2.148, 2.258 | 2.220, 2.172, 2.149, 2.257 |
| regular | difference_W | 8→16 | 2.965e-04, 5.399e-04, 6.579e-04, 1.819e-03 | 0.467, 0.986, 1.856, 1.081 | 0.467, 0.977, 1.849, 1.076 |
| regular | difference_W | 16→32 | 1.238e-04, 1.191e-04, 2.146e-04, 3.813e-04 | 1.260, 2.181, 1.616, 2.254 | 1.255, 2.176, 1.612, 2.249 |
| regular | difference_W | 32→64 | 2.628e-05, 1.946e-05, 4.041e-05, 5.453e-05 | 2.237, 2.613, 2.409, 2.806 | 2.232, 2.607, 2.401, 2.797 |
| jitter | W_primitive | 8→16 | 6.001e-03, 2.242e-02, 1.766e-02, 8.223e-03 | 1.452, 1.799, 1.813, 1.469 | 1.447, 1.797, 1.809, 1.465 |
| jitter | W_primitive | 16→32 | 1.769e-03, 5.915e-03, 5.667e-03, 2.423e-03 | 1.763, 1.923, 1.640, 1.763 | 1.762, 1.922, 1.641, 1.762 |
| jitter | W_primitive | 32→64 | 4.418e-04, 1.304e-03, 1.292e-03, 6.056e-04 | 2.001, 2.182, 2.134, 2.001 | 2.001, 2.183, 2.135, 2.001 |
| jitter | W_conservative | 8→16 | 6.003e-03, 2.232e-02, 1.754e-02, 9.363e-03 | 1.450, 1.806, 1.803, 1.416 | 1.445, 1.803, 1.799, 1.412 |
| jitter | W_conservative | 16→32 | 1.770e-03, 5.873e-03, 5.623e-03, 2.587e-03 | 1.762, 1.926, 1.642, 1.856 | 1.762, 1.926, 1.643, 1.856 |
| jitter | W_conservative | 32→64 | 4.426e-04, 1.303e-03, 1.286e-03, 6.323e-04 | 2.000, 2.172, 2.128, 2.032 | 2.000, 2.174, 2.129, 2.033 |
| jitter | difference_W | 8→16 | 4.929e-04, 8.189e-04, 9.197e-04, 2.158e-03 | -0.387, 1.460, 1.327, 0.918 | -0.390, 1.453, 1.323, 0.915 |
| jitter | difference_W | 16→32 | 1.776e-04, 3.038e-04, 2.687e-04, 4.546e-04 | 1.473, 1.431, 1.775, 2.247 | 1.471, 1.430, 1.773, 2.244 |
| jitter | difference_W | 32→64 | 4.833e-05, 9.288e-05, 7.828e-05, 9.626e-05 | 1.877, 1.710, 1.779, 2.240 | 1.877, 1.709, 1.781, 2.241 |

Tutti i valori assoluti estrapolati e tutte le altre norme/regioni/tempi,
incluse le traiettorie d=1,2,4, sono negli `errors` e `orders`. Il bordo
conservativo è prescritto: non si interpreta il suo errore nullo come
ordine di accuratezza misurato. Il bordo primitivo è invece misurato.

Il controllo triangolare E_PC≤E_P+E_C passa per ogni norma, componente,
regione, tempo e variante. È un controllo del postprocessing, non una
validazione autonoma dei solver. A T=0.5, globalmente in L2, il rapporto
E_PC/(E_P+E_C) è 0.0073–0.1179 e nessun errore reciproco cade vicino
al roundoff. Segnali prossimi a 1e-12 o inferiori a 1e-8 del denominatore
sono comunque marcati nella tabella di controllo, senza attribuire loro
un ordine affidabile. Gli ordini E_PC possono restare sensibili a
cancellazioni dei coefficienti degli errori dei due solver.

Su jitter l'ultimo p_PC di rho/u/v resta sotto due mentre entrambi gli
errori rispetto all'analitica mostrano p vicino a due. Non c'è contraddizione:
si tratta di p su intervalli finiti, non di bounds asintotici uniformi.
Se entrambi gli errori avessero bounds O(h²) dimostrati, la disuguaglianza
triangolare fornirebbe un bound O(h²) anche per la differenza; questo
non fissa il rapporto di errori su queste quattro mesh.

## Positività, stabilità, bordi e costo

Su tutta la campagna, inclusi i passi iniziali più grossolani:
rho_P≥0.9043832, p_P≥0.8687403,
rho_C≥0.9043031, p_C≥0.8685682. Zero fallback e contrazioni;
theta minimo=1. Il massimo dt*lambda_i/V_i primitivo osservato è
0.509773. Questi dati non sono una garanzia di stabilità per altri vortici.

Tempi di stadio esattamente concordi, errore del vincolo forte zero,
aggiornamento nativo FE verificato entro 2.23e-16. I registri discreti
chiudono entro 3.98e-13 (primitivo) e 3.42e-13 (conservativo), su area
totale 100. Come nelle FASI 6 includono l'iniezione del bordo forte e
il resto nonlineare della trasformazione primitiva: non sono prove di
conservazione fisica o di correttezza degli shock.

A N=64,T=0.5 il contributo del bordo a C_R=ΣV e² primitivo è circa
2.7e-7–4.7e-7 del totale; ring1 pesa circa 1.0e-7–1.3e-6.
La regione mobile r≤2 contiene circa il 90–96% dell'energia dell'errore
dei due solver. Per la differenza reciproca ne contiene circa 95–99%.
Questi valori riducono l'evidenza di una contaminazione dominante diretta
dal bordo; non isolano causalmente il suo effetto propagato nell'interno.

N=64,d=4 ha richiesto circa 40 s per coppia; d=8 circa 91 s sotto
concorrenza, ancora sostenibile. La somma dei tempi interni dei 32 run
accoppiati è circa 343.5 s; non è il wall time complessivo della sessione
con due processi concorrenti, build, writer e test. Non sono stati
avviati raffinamenti superiori, come richiesto.

## Diagnostiche fisiche: fase, ampiezza e forma

Il minimo continuo esatto è p=0.8689305, rho=0.9045194;
la massima perturbazione di velocità è beta/(2pi)=0.3183099 a r=1.
I minimi esatti **nodali** dipendono dalla posizione del centro e dalla mesh;
la tabella usa per l'attenuazione il deficit di pressione esatto campionato
sugli stessi nodi, distinto da quello continuo.

A N=64,T=0.5:

| Mesh | Solver | p_min | rho_min | Ampiezza pressione / esatta nodale | Errore centro nodale | L2 entropia | Larghezza | Anisotropia |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| regular | primitive | 0.870978 | 0.90606 | 0.991714 | 0.0883883 | 3.87363e-06 | 0.971469 | 1.00423 |
| regular | conservative | 0.8711 | 0.905952 | 0.990781 | 0.0883883 | 2.53834e-05 | 0.971259 | 1.00405 |
| regular | exact | 0.8699 | 0.90524 | 1 | 0.0883883 | 8.66147e-17 | 0.968264 | 1.00101 |
| jitter | primitive | 0.871609 | 0.906525 | 0.987833 | 0.0940404 | 4.15824e-06 | 0.97132 | 1.00385 |
| jitter | conservative | 0.872059 | 0.906579 | 0.984372 | 0.0940404 | 3.00055e-05 | 0.971172 | 1.00369 |
| jitter | exact | 0.870028 | 0.905335 | 1 | 0.0940404 | 8.44061e-17 | 0.968577 | 1.0013 |

Il centro stimato con argmin p coincide, in entrambi i solver e in tutti
i casi finali riportati, con il nodo del minimo analitico campionato.
L'errore rispetto al centro continuo è circa 0.08839 su regular e
0.09404 su jitter a N=64; è già presente nell'analitica campionata.
Non si attribuisce precisione subgrid né si dichiara nullo l'errore di fase.
Nei CSV sono presenti i nodi, le coordinate e la distanza fra i due minimi.

Si registra inoltre il baricentro pesato da max(1−p,0) nella regione mobile,
come descrittore distinto dal minimo. La matrice di covarianza di questi
pesi fornisce larghezza sqrt(traccia) e anisotropia sqrt(lambda_max/lambda_min).
Anche questi descrittori vanno confrontati con l'analitica campionata, che
ha anisotropia non esattamente uno su mesh finite, specialmente jitter.
Non sono una localizzazione del centro con accuratezza certificata.

Il conservativo attenua maggiormente il deficit di pressione a N=64:
rapporto 0.99078 contro 0.99171 su regular e 0.98437 contro 0.98783
su jitter. L'errore di entropia è maggiore nel conservativo. La larghezza
e l'anisotropia sono però molto simili e leggermente minori nel conservativo;
il massimo di velocità non assegna lo stesso vincitore su ogni famiglia.
Non si conclude che un solver sia universalmente «più diffusivo» sulla
base di un solo indicatore. L'entropia da sola non misura l'errore globale.

## Grafici e visualizzazione ParaView

I due grafici log-log mostrano tutte le componenti in L2 globale,
tre errori distinti e due varianti (FE e Richardson). La linea O(h²)
è un riferimento ancorato al dato fine primitivo, **non un fit imposto**:

- [Regular](primitive_euler_phase7_convergence_regular.svg)
- [Jitter](primitive_euler_phase7_convergence_jitter.svg)

Sequenze principali, disponibili anche per N=8,16,32 nella stessa cartella:

- [Regular N=64](primitive_euler_phase7_vortex/regular_64.pvd)
- [Jitter N=64](primitive_euler_phase7_vortex/jitter_64.pvd)

Ogni sequenza contiene t=0,0.25,0.5 e la medesima mesh effettiva per
primitivo, conservativo e analitica. I campi includono rho/u/v/p dei tre
stati, errori firmati di entrambi, differenza firmata, modulo di velocità,
entropia, ID nodale, volume duale e tempo fisico. Nessun modulo sostituisce
silenziosamente il campo firmato.

Controlli automatici: numero e ordine dei punti, numero/tipo/connettività
dei triangoli, coordinate e jitter effettivi, valori finiti, tempi e indici,
roundtrip esatto dei campi, corrispondenza con i campi delle norme. Il
lettore XML è stato seguito dal vero `vtkXMLUnstructuredGridReader` e da
`PVDReader` di ParaView 5.13.2, con verifica dei valori contro gli snapshot.
Le norme globali sono state inoltre ricalcolate indipendentemente dai
campi errore letti da VTK e confrontate con i CSV (tolleranza relativa
2e-14 e assoluta 2e-16), anche sui file archiviati nel repository:
massima differenza osservata 6.94e-18.
La sandbox inizialmente impediva il socket locale del runtime MPI di
ParaView; la stessa lettura è riuscita dopo autorizzazione, senza cambiare
i file. Nessuna apertura grafica o immagine renderizzata viene dichiarata.

Per la verifica visiva manuale:

1. Aprire il `.pvd` e premere Apply; scegliere Surface.
2. Usare `p_primitive`, poi `p_conservative`, poi `p_exact`, mantenendo
   **la stessa scala** (per esempio 0.86–1.00) e gli stessi livelli di Contour.
3. Per l'errore selezionare `error_p_primitive` o `error_p_conservative`;
   per il confronto reciproco `difference_p_primitive_conservative`.
   Usare una mappa divergente con estremi simmetrici rispetto a zero.
4. Ripetere per rho/u/v; mantenere scala e camera uguali fra i tre stati.
5. Usare Play per i tre snapshot. È una sequenza temporale sparsa,
   non una descrizione risolta di ogni dettaglio del trasporto.

Contorni e animazione aiutano l'interpretazione ma non sostituiscono
norme, controlli temporali e diagnosi quantitative.

## Risposte alle domande scientifiche

**A — Primitivo contro analitica.** Converge su entrambe le famiglie.
Gli ultimi p L2 Richardson sono (2.220,2.210,2.190,2.224) su regular e
(2.012,2.210,2.150,2.012) su jitter. Le transizioni grossolane mostrano
ordini inferiori e sono riportate integralmente.

**B — Conservativo contro analitica.** Converge alla stessa soluzione:
(2.224,2.203,2.185,2.260) su regular e (2.011,2.201,2.147,2.043)
su jitter nell'ultima transizione. Si tratta del solver SU2 effettivo.

**C — Differenza reciproca.** Diminuisce: (2.237,2.624,2.411,2.801)
su regular; (1.881,1.715,1.796,2.250) su jitter. Non è una conferma
del secondo ordine reciproco in tutte le componenti jitter.

**D — Compatibilità col secondo ordine.** I confronti indipendenti con
l'analitica lo supportano nell'ultimo intervallo per entrambi i solver;
la sensibilità temporale misurata è piccola. Quattro mesh, nucleo
grossolanamente risolto ai primi livelli e un solo T non costituiscono
una prova generale dell'ordine asintotico.

**E — Diffusione, fase e deformazione.** Il conservativo attenua di più
la pressione e ha errore entropico maggiore nei dati fini; non emerge
un vincitore universale per fase o deformazione. Il minimo nodale non
risolve differenze subgrid; le covarianze di forma sono vicine.

**F — Attribuzione causale.** Il tempo non spiega le differenze principali
entro i controlli effettuati. Il contributo diretto del bordo è piccolo,
ma i trattamenti forte/debole restano differenti. La maggiore differenza
su jitter è compatibile con una sensibilità alla geometria; non isola
WLS, quadratura o fluttuazioni come causa unica. Servirebbero esperimenti
discriminanti ulteriori, non implementati qui.

Rispetto alle MMS delle FASI 6, questo caso non ha sorgenti né pressioni
manufactured indipendenti e ha il nucleo lontano dal bordo. Conferma
con dati propri la convergenza liscia a questi livelli, senza cancellare
i difetti semidiscreti jitter osservati prima o dedurne automaticamente
il comportamento del vortice. Non offre conclusioni su shock o
condizioni Rankine–Hugoniot.

## Test, file e riproducibilità

Esiti: analitica vortice 890 controlli; writer 4 test, compresi casi
negativi di connettività, tempo e valori non finiti; pilota e 32 coppie
di traiettorie validate; controllo triangolare completo; lettura VTK/PVD
di 24/8 file. Regressioni precedenti: 70 casi, 22 253 237 asserzioni;
MMS analitica precedente: 3300 asserzioni; postprocessing Richardson
precedente: 3 test. Nessuna soglia precedente modificata.

Nuovi file sorgente/strumenti:

- `UnitTests/SU2_CFD/euler_isentropic_vortex.hpp` e `euler_vortex_analytic_tests.cpp`;
- `UnitTests/run_euler_vortex.py` e `run_euler_vortex_campaign.py`;
- `UnitTests/export_euler_vortex.py` e `euler_vortex_output.py`;
- `UnitTests/validate_euler_vortex_paraview.py` e `test_euler_vortex_output.py`.

CSV `primitive_euler_phase7_*`: `physics`, `meshes`, `checks` (tempi e passi
effettivi inclusi), `errors`, `orders`, `temporal`, `temporal_initial`,
`diagnostics`, `triangle`, `vtu_checks`, `paraview_checks`.
Il `provenance.json` conserva hash di build, pilota e dei 41 file precedenti
preservati; il report, i due SVG e la cartella VTU/PVD completano la consegna.
I campi grezzi e il sorgente temporaneo restano nello scratch; i VTU/PVD
finali sono anche archiviati nel repository e non dipendono dalla sua durata.

```sh
clang++ -std=c++11 -O2 UnitTests/SU2_CFD/euler_vortex_analytic_tests.cpp -o /private/tmp/su2-vortex-analytic
/private/tmp/su2-vortex-analytic
python3 UnitTests/run_euler_vortex.py --output-dir /private/tmp/su2-phase7/build
python3 UnitTests/run_euler_vortex_campaign.py --executable /private/tmp/su2-phase7/build/euler_vortex --output-dir /private/tmp/su2-phase7/campaign --pilot
python3 UnitTests/run_euler_vortex_campaign.py --executable /private/tmp/su2-phase7/build/euler_vortex --output-dir /private/tmp/su2-phase7/campaign --run --jobs 2 --resume
python3 UnitTests/export_euler_vortex.py --input-dir /private/tmp/su2-phase7/campaign --output-dir /private/tmp/su2-phase7/results
python3 UnitTests/run_euler_vortex_campaign.py --executable /private/tmp/su2-phase7/build/euler_vortex --output-dir /private/tmp/su2-phase7/campaign --run --jobs 2 --resume --divisors 1 2 4 8
python3 UnitTests/export_euler_vortex.py --input-dir /private/tmp/su2-phase7/campaign --output-dir /private/tmp/su2-phase7/final --divisors 1 2 4 8
/Applications/ParaView-5.13.2.app/Contents/bin/pvpython UnitTests/validate_euler_vortex_paraview.py --vtu-dir /private/tmp/su2-phase7/final/vtu --input-dir /private/tmp/su2-phase7/campaign --divisor 8 --errors /private/tmp/su2-phase7/final/errors.csv --report /private/tmp/su2-phase7/final/paraview_checks.csv
python3 UnitTests/test_euler_vortex_output.py
python3 UnitTests/test_euler_asymptotic_postprocessing.py
/private/tmp/su2-phase5-tests/test_driver --reporter compact
/private/tmp/su2-mms-6b1-final/euler_mms_analytic --reporter compact
git diff --check
git status --short
```

La posizione di ParaView è quella installata in questo ambiente, non una
dipendenza del solver. Build e regression driver preesistenti sono gli
stessi usati nelle FASI 6; gli script ricavano flag e archivi dalla build
locale SU2. La compilazione ha soltanto gli avvisi già noti sulle warning
GCC ignorate da Clang e librerie duplicate. Il primo rendering Matplotlib
ha segnalato cache font non scrivibili; i grafici sono stati generati,
e la configurazione finale usa cache in scratch senza cambiare i dati.

Stato finale: esclusivamente nuovi file della FASE 7, nessuna modifica
tracciata ai file preesistenti e HEAD invariato. Nessuna sincronizzazione
Git e nessun intervento sul precedente problema HTTPS. Ci si ferma qui;
nessuna fase o correzione numerica ulteriore è avviata automaticamente.
