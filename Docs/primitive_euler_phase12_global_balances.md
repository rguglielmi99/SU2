# FASE 12 — Bilanci globali in dominio chiuso

## Perimetro e diagnostiche

Modificato solo `UnitTests/SU2_CFD/primitive_euler_boundary_tests.cpp`: tre
test `[primitive-global-balance]`, fixture e API esistenti. Nuovi file: questo
rapporto e `primitive_euler_phase12_global_balances.csv` (17 righe).
**Nessuna modifica produttiva**, correzione conservativa, sorgente, nuovo
output, confronto nativo o variazione di BC/gradienti/FE.

Si integrano coi volumi duali le densità diagnostiche

    Q(W) = (rho, rho*u, rho*v, p/(gamma-1)+rho*(u²+v²)/2,
            rho*(ln(p)-gamma*ln(rho))), gamma=1.4.
    Q_h = sum_i V_i Q(W_i),    dot(Q_h) = sum_i V_i DQ(W_i) rhs_i.

Il RHS è quello effettivo **dopo** la reazione nodale. La regola della catena
è implementata per tutte e cinque le grandezze, senza differenze temporali.
Un test separato la verifica mediante derivate direzionali centrate su uno
stato arbitrario; queste differenze verificano la diagnostica, non sostituiscono
i tassi. Le grandezze sono integrali discreti da campioni NODE: il loro errore
iniziale di quadratura è distinto dal drift, riferito sempre a Q_h(0).

Per la quantità di moto si calcola F_h=sum_b p_rec(q_b) S_b, con tracce
ricostruite sulle mezze facce reali. Il difetto semidiscreto è dot(P_h)+F_h,
non dot(P_h). Durante FE, I_F=sum_j dt*F_h(W_j) usa gli stessi stati iniziali
di ciascun passo. Il suo integrale temporale è approssimato al primo ordine.

## Casi e mesh

Quadrato [0,1]², quattro marker Euler distinti; triangoli regolari lungo la
stessa diagonale, N=8/16/32, (N+1)² nodi e 2N² triangoli. h=1/N (il massimo
lato è sqrt(2)/N, con gli stessi rapporti). NODE/AGGREGATED/NODE_POINT,
MUSCL affine, WLS, limiter NONE, nessuna sorgente.

**A: quiete**, rho=p=1, u=v=0, N=8. RHS, tassi, forza e somma delle normali
sono nulli nel calcolo; dopo cinque passi dt=.001, M=1, P=0, E=2.5, H=0
restano invariati. Positività e slip verificati. È un equilibrio particolare,
non una prova di conservazione per dati generali.

**B: dato iniziale liscio non stazionario** richiesto:

    rho=1+.03 cos(pi*x) cos(pi*y)
    u=.1 sin(pi*x) cos(pi*y), v=-.1 cos(pi*x) sin(pi*y)
    p=1+.04 cos(pi*x)+.03 cos(pi*y).

Soddisfa slip sui lati e ai vertici; le pressioni sulle pareti danno
analiticamente F=(-.08,-.06), quindi dot(P)=(.08,.06) all'istante iniziale.
Non si assume che questa forza resti costante durante l'evoluzione.

## Tassi semidiscreti iniziali

| N | abs(dot M) | abs(dot E) | ordine E | abs(dot H) | ordine H |
|---:|---:|---:|---:|---:|---:|
| 8 | 2.894e-18 | 1.86547e-4 | — | 1.13953e-5 | — |
| 16 | 6.685e-18 | 2.55769e-5 | 2.867 | 1.69258e-6 | 2.751 |
| 32 | 4.474e-18 | 3.31434e-6 | 2.948 | 2.27955e-7 | 2.892 |

dot(E) è negativo, dot(H) positivo; segni e valori completi nel CSV.
Gli ordini di questi bilanci integrali possono beneficiare di cancellazioni:
non misurano l'ordine puntuale del residuale o l'accuratezza della soluzione.
Non si assegna ordine alla massa iniziale: è al livello del roundoff e
non costituisce un risultato generale, come mostra l'evoluzione successiva.

| N | dot Px | dot Py | dot Px+Fx | dot Py+Fy |
|---:|---:|---:|---:|---:|
| 8 | .07825095 | .05884505 | -1.74905e-3 | -1.15495e-3 |
| 16 | .07958777 | .05972056 | -4.12235e-4 | -2.79444e-4 |
| 32 | .07989932 | .05993022 | -1.00676e-4 | -6.97773e-5 |

Ordini assoluti dei difetti: Px 2.085/2.034, Py 2.047/2.002.
Sia F ricostruita sia la quadratura con p analitica nei medesimi punti
coincidono con (-.08,-.06) entro 3.1e-16: in questo caso il difetto iniziale
di quantità di moto non è spiegato dall'errore della forza ricostruita.
L'esattezza osservata dipende anche dalla forma separabile del dato.

## Evoluzione breve: drift e contributo temporale

T=.05, CFL iniziale .1, dt costante aggiustato per raggiungere T:
51/102/203 passi per N=8/16/32; ripetizione N=16 con 204 passi.

| N | dt | Delta M/M0 | Delta E/E0 | Delta H | Delta Px+Ix | Delta Py+Iy |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 9.80392e-4 | -8.2548e-9 | -3.90718e-6 | 2.36311e-7 | -1.24865e-4 | -8.49869e-5 |
| 16 | 4.90196e-4 | -4.3344e-9 | -5.31947e-7 | 1.82038e-8 | -3.82670e-5 | -2.66237e-5 |
| 16, dt/2 | 2.45098e-4 | -4.3978e-9 | -5.63059e-7 | 1.87505e-8 | -3.83255e-5 | -2.66662e-5 |
| 32 | 2.46305e-4 | -8.2760e-10 | -5.54357e-8 | -3.65276e-9 | -1.25240e-5 | -8.87331e-6 |

Si registra anche T_Q=Q_h(T)-Q_h(0)-sum_j dt*dot(Q_h)(W_j): è il resto
del cambiamento di variabili sotto FE, non il solo errore rispetto a una
soluzione semidiscreta esatta. La traiettoria e la quadratura dei tassi
dipendono a loro volta da dt.

Per N=16, l'integrale del tasso energetico è -1.48417e-6, mentre T_E è
+1.52977e-7; dimezzando dt diventano -1.48553e-6 e +7.64745e-8.
Il contributo temporale positivo compensa parzialmente la perdita spaziale:
il drift energetico totale può quindi peggiorare lievemente riducendo dt.
T_H=-1.80129e-9 passa a -9.03219e-10. Il resto della massa, lineare nelle
primitive, è roundoff; i resti non lineari di P/E sono verificati contro
le formule polinomiali esatte del singolo passo, senza soglie conservative.

Al tempo finale dot(M) è -1.139e-6/-3.514e-7/-5.499e-8 per N=8/16/32:
la cancellazione iniziale non persiste. Delta H cambia segno al raffinamento;
non si deduce una proprietà di entropia discreta. Tutti i passi sono
ammissibili; ai checkpoint min(rho)=.97, min(p)=.929632 circa,
massimo slip=1.23e-17. Nessun rifiuto o fallback.

## Riproducibilità e valutazione

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase12-build
/private/tmp/primitive-phase12-build/test_driver '[primitive-global-balance]'
/private/tmp/primitive-phase12-build/test_driver '[primitive-boundary],[primitive-geometry],[primitive-wall-gradients],Constant primitive fields*,Affine NODE*,Forward Euler applies*,Rejected time steps*'
git diff --check
```

Build riuscita. Nuovi test: 3 casi/3807988 asserzioni; regressioni ridotte:
25 casi/17428 asserzioni, tutti passati. Nessuna campagna estesa. Il CSV
riproduce le righe GLOBAL_BALANCE senza prefisso; `int_*` sono somme dt*tasso,
`time_*` i resti T_Q e `*_impulse_defect` i bilanci integrati della quantità
di moto. Nessuna normalizzazione relativa dell'entropia, inizialmente piccola.

**A. Esattezza discreta:** quiete preservata nel campione; nessuna prova di
conservazione esatta per flussi generali.
**B. Conservazione approssimata:** difetti iniziali energetici/entropici e
bilancio di quantità di moto diminuiscono; i drift si riducono sulle tre mesh,
ma con dt variabile non si assegna loro un ordine puramente spaziale.
**C. Limiti:** difetti non nulli a mesh finita, tasso di massa durante il moto
ed entropia senza segno uniforme richiedono ulteriori verifiche per altre
soluzioni o tempi più lunghi. Nessuna evidenza di plateau persistente da questi
tre raffinamenti. Non si estendono i risultati agli shock. Nessun commit/push.
