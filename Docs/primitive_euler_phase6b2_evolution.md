# FASE 6B.2 — Evoluzione della MMS comune con Forward Euler

## Esito e perimetro

Eseguite 36 coppie di traiettorie, fino a T=0.1: tre famiglie triangolari,
N=8,16,32,64 e tre passi temporali per mesh. Entrambi i solver convergono
verso la soluzione analitica e la loro differenza diminuisce. Il controllo
dell'errore temporale è superato. **Non è dimostrato ordine globale almeno
due per tutte le componenti**: sul raffinamento 32→64, per esempio, la
pressione primitiva ha ordine 1.943 sulla famiglia regolare; il conservativo
ha ordini globali fra 1.741 e 1.980 nelle primitive.

Questi risultati non autorizzano correzioni numeriche né una campagna shock.
Si conclude la campagna prevista e si rimanda alla revisione scientifica
l'interpretazione del regime asintotico. Non si deduce l'ordine della soluzione
evoluta dall'ordine del difetto semidiscreto misurato nella FASE 6B.1.

Branch `feature/primitive-path-solver`, HEAD di riferimento
`5be7a6e72a6c082cce192e73e8dd53eb17f1b13b`. Nessuna modifica ai file di
produzione, nessun commit, push o sync. NODE resta la rappresentazione
implementata; la decisione matematica definitiva NODE/DUAL_MEAN resta aperta.

## Problema e ciclo di avanzamento

Si riutilizzano analitica, Jacobiano, sorgenti e geometrie di
[FASE 6B.1](primitive_euler_phase6b1_mms.md), senza duplicare le derivate:

\[
a(t)=0.1\cos t,\quad b(t)=0.1\sin t,\qquad
W^*=\begin{pmatrix}
1+a\sin(2\pi x)\sin(2\pi y)\\
0.5+a\cos(2\pi x)\sin(2\pi y)\\
0.3+a\sin(2\pi x)\cos(2\pi y)\\
1+b\cos(2\pi x)\cos(2\pi y)
\end{pmatrix}.
\]

Gamma=1.4; coordinate, volumi e unità fisiche coincidono. Si inizializzano
W e U(W) puntualmente a t=0. La sorgente primitiva è
S_W=W_t+A_x W_x+A_y W_y; quella conservativa è S_U=J(W)S_W.
Il riferimento è valutato a T, per ogni componente, nel vertice fisico.
Non si confrontano qui incognite con medie duali analitiche.

Si usano Roe nativo, MUSCL affine, WLS, limiter NONE; per il primitivo
NODE/AGGREGATED/NODE_POINT. Il banco collega realmente `CEulerSolver` e
`CPrimitiveEulerSolver` alle librerie SU2 compilate: non contiene un secondo
operatore conservativo. Non usa il driver generale SU2.

Per il conservativo ogni passo esegue `Set_OldSolution`, `BC_Custom(t_n)`,
`Preprocessing`, `Upwind_Residual`, `Source_Residual`, `BC_Custom(t_n)` e
`ExplicitEuler_Iteration`. Il secondo richiamo al bordo azzera le righe
vincolate dopo l'assemblaggio. Dopo FE si impone `BC_Custom(t_(n+1))`:
il vincolo forte deve seguire il tempo fisico. Si controlla il candidato
grezzo prima di qualsiasi successivo preprocessing che potrebbe ripristinare
stati precedenti. I tempi del verification adapter sono verificati.

Il primitivo chiama `AdvanceExplicitEuler(t_n,dt,exterior,source)`; la
condizione esterna esatta è valutata nei punti delle mezze facce di bordo,
con trattamento debole originale. Nessun vincolo forte viene aggiunto.
Le condizioni sono analiticamente coerenti ma numericamente diverse:
il confronto globale include questa differenza, non dimostra equivalenza
a mesh finita. Il campo conservativo sul bordo forte ha errore nullo salvo
roundoff; il suo ordine di bordo non viene stimato.

`InnerIter=0` disabilita soltanto la stampa nativa degli errori di verifica
durante il loop; gli errori finali sono calcolati indipendentemente a T.
Ogni passo nativo è controllato contro U_new=U_old-dt*R/V. Nessuna riduzione
nascosta di dt, clipping, rollback, ricostruzione nuova o correzione conservativa.

## Risoluzione temporale e criteri

Per ogni mesh e d=1,2,4:

\[
n_t=d\left\lceil T N^2/0.25\right\rceil,\qquad dt=T/n_t.
\]

Il passo fisico uniforme è uguale nei due solver. Il parametro CFL non viene
usato per inferirne l'uguaglianza. I risultati spaziali pubblicati usano d=4:
n_t=104,412,1640,6556 rispettivamente; l'intera campagna copre dt fra
1.5253203173e-5 e 3.8461538462e-3. A mesh fissata si confrontano i campi,
non la semplice differenza dei loro errori scalari:

\[
\delta_c=\|q_{dt}-q_{dt/2}\|,\quad
\delta_f=\|q_{dt/2}-q_{dt/4}\|,\quad
p_t=\log_2(\delta_c/\delta_f).
\]

Per ciascuna componente dei quattro campi primari, regione e mesh, si richiede
delta_f/E_fine ≤ 0.05 nella L2 pesata. Per ogni coppia di mesh consecutive,
si richiede |p_(dt/2)-p_(dt/4)| ≤ 0.1. Il bordo conservativo imposto è escluso;
1e-12 è soltanto una soglia dichiarata per classificare errori di roundoff,
non una modifica al metodo o al riferimento. La differenza reciproca fra
solver è registrata ma non usata come denominatore di accettazione.

Risultati: **1248 verifiche L2 temporali e 936 verifiche di stabilità
dell'ordine, nessun fallimento**. Massimo delta_f/E=0.01168258; massima
variazione d'ordine=0.00815810. Anche L∞ supera gli stessi limiti:
massimi 0.01672143 e 0.01475798. Gli ordini temporali L2 risolvibili sono
0.9971–1.0822, mediana 1.0124: coerenti con Forward Euler.

Questa è una verifica quantitativa di indipendenza entro la tolleranza
dichiarata, non un'eliminazione esatta dell'errore temporale. Poiché dt
decresce come h², la sola curva sotto raffinamento spaziale non sarebbe
sufficiente a separare i due errori: sono essenziali i dimezzamenti a mesh fissa.

## Norme, regioni e risultati spaziali

Per una componente e_i: L1=ΣV_i|e_i|/ΣV_i,
L2=(ΣV_i e_i²/ΣV_i)^(1/2), L∞=max|e_i|. RMS è invece
(Σe_i²/n)^(1/2), non pesata. Si conserva il nome delle colonne precedente.
Gli ordini sono log(E_N/E_2N)/log2.

Le sette regioni sono: tutti i nodi, bordo fisico, primo anello, secondo
anello, interno da tre strati, regione fisica fissa [0.25,0.75]² e tutti
i nodi interni. Gli anelli e deep3 cambiano dominio fisico sotto raffinamento;
la regione fissa aiuta a controllare questa ambiguità. La famiglia `smooth`
deforma dolcemente i vertici, `jitter` introduce perturbazioni O(h) degli
interni; sono le stesse definizioni deterministiche del banco 6B.1.

L2 globale vettoriale delle primitive, somma in quadratura delle quattro
componenti (unità di riferimento unitarie; i CSV restano per componente):

| Famiglia | N | W primitivo | W(U conservativo) | Differenza W |
|---|---:|---:|---:|---:|
| regular | 8 | 1.52670846e-2 | 1.22065809e-2 | 1.28121188e-2 |
| regular | 16 | 3.40784411e-3 | 3.18211953e-3 | 3.28231953e-3 |
| regular | 32 | 7.62468792e-4 | 8.95028357e-4 | 8.24203757e-4 |
| regular | 64 | 1.82010791e-4 | 2.44492977e-4 | 2.11512393e-4 |
| smooth | 8 | 1.71935921e-2 | 1.43951298e-2 | 1.40666260e-2 |
| smooth | 16 | 4.06518422e-3 | 3.91446911e-3 | 3.80717802e-3 |
| smooth | 32 | 9.02675586e-4 | 1.08375244e-3 | 9.95099173e-4 |
| smooth | 64 | 2.12379274e-4 | 2.93535673e-4 | 2.60034534e-4 |
| jitter | 8 | 1.55817703e-2 | 1.21134912e-2 | 1.26127628e-2 |
| jitter | 16 | 3.73078641e-3 | 3.43816886e-3 | 3.34984731e-3 |
| jitter | 32 | 8.48126356e-4 | 9.43037740e-4 | 8.36894935e-4 |
| jitter | 64 | 2.02273097e-4 | 2.55042646e-4 | 2.13316044e-4 |

Ordini L2 globali 32→64, componenti (rho,u,v,p):

| Famiglia | Primitivo | Conservativo | Differenza W |
|---|---|---|---|
| regular | 2.057, 2.203, 2.239, 1.943 | 1.799, 1.980, 1.977, 1.741 | 2.022, 1.885, 1.888, 1.982 |
| smooth | 2.099, 2.180, 2.228, 1.956 | 1.834, 1.974, 1.963, 1.789 | 1.966, 1.877, 1.866, 1.976 |
| jitter | 2.064, 2.183, 2.198, 1.968 | 1.845, 1.979, 1.984, 1.759 | 2.038, 1.879, 1.906, 1.990 |

Sulla regione fisica fissa regolare gli ordini sono rispettivamente
(1.820,2.218,2.232,1.797) e (1.842,2.227,2.244,1.832).
La differenza W nella stessa regione scende a N=64 a
(5.519e-6,3.388e-6,2.581e-6,2.711e-6).
Quindi una differenza piccola nell'interno non implica da sola secondo ordine
di entrambi rispetto all'analitica. Non si attribuisce causalmente l'ordine
inferiore a due a un singolo termine senza ulteriori esperimenti.

I CSV includono anche U(W_prim), U_cons e la loro differenza, oltre alle
regioni di bordo e agli anelli. Nessuna equivalenza nonlineare fra medie
primitive e conservative viene assunta. I risultati evoluti non annullano
i difetti semidiscreti precedentemente osservati su mesh irregolari.

## Positività, tempi e bilanci diagnostici

In tutti i passi: min rho primitivo=0.8994332, min p primitivo=0.9658617;
min rho conservativo=0.8994763, min p conservativo=0.9826445.
Theta minimo=1, zero contrazioni, zero fallback nativi su nodi o facce.
Il massimo dt*lambda_i/V_i primitivo sull'intera campagna è 0.4820583;
non è un'affermazione di stabilità per altri dati o mesh.

Errore massimo del tempo di stadio=1.39e-17; coordinate e volumi fra solver
coincidono esattamente; errore Dirichlet finale=0; errore dell'aggiornamento
nativo FE ≤2.23e-16. Si verificano i conteggi esatti delle chiamate a sorgenti
e bordi, il tempo finale e l'assenza di passi locali diversi dal dt prescritto.

Posto Q=ΣV_i U_i, il bilancio conservativo registrato è
Q^m-Q^0+Σ_n dt R^n-Σ_n I_boundary^n, dove R contiene le righe realmente
usate da FE e I_boundary è l'incremento imposto dal vincolo a t_(n+1).
Il massimo difetto è 4.00e-15. L'iniezione del vincolo è una contabilità
del serbatoio Dirichlet, **non** il flusso fisico integrato sul bordo.

Per il primitivo si registra separatamente ΣV dt J(W^n) rhs e il resto
nonlineare ΣV[U(W^(n+1))-U(W^n)-dt J(W^n)rhs]. La loro somma riproduce
Q^m-Q^0 entro 4.89e-15. È un'identità diagnostica dell'aggiornamento,
non una prova di conservazione fisica e non un termine correttivo applicato.
Gli scarti da ΣV_i U_exact(x_i,T) sono pubblicati per componente, mantenendo
distinta la quadratura del riferimento da un integrale continuo esatto.
Nessuna conclusione Rankine–Hugoniot deriva da questa MMS liscia.

## Riproducibilità, file e verifiche

File nuovi di questa fase:

- `UnitTests/SU2_CFD/euler_mms_evolution.cpp`: banco accoppiato produttivo;
- `UnitTests/run_euler_mms_evolution.py`: compilazione con flag e archivi SU2;
- `UnitTests/export_euler_mms_evolution.py`: campagna, dimezzamenti, ordini e gate;
- questo rapporto e quattro CSV `primitive_euler_phase6b2_{errors,temporal,order_stability,checks}.csv`.

`errors` contiene i 2016 record della risoluzione temporale più fine;
`temporal` contiene differenze dei campi e criteri per ogni norma;
`order_stability` confronta gli ordini alle ultime due risoluzioni temporali;
`checks` contiene i controlli di tutte le 36 traiettorie accoppiate.
I campi grezzi, log e risultati per ogni dt sono nella directory temporanea
della campagna; si rigenerano con i comandi seguenti (la directory temporanea
non è un archivio permanente). Il manifest impedisce di mescolare eseguibili
diversi durante una ripresa. È stato creato alla raccolta finale, dopo i run;
SHA256 dell'eseguibile usato:
`2264c0d7449cacf0a6ef8d7afed2332f5fbba6a44a99cdea275b5f89843b0efd`.

```sh
python3 UnitTests/run_euler_mms_evolution.py --output-dir /private/tmp/su2-mms-6b2-complete
python3 UnitTests/export_euler_mms_evolution.py --executable /private/tmp/su2-mms-6b2-complete/euler_mms_evolution --output-dir /private/tmp/su2-mms-6b2-complete/campaign --jobs 2
python3 UnitTests/export_euler_mms_evolution.py --executable /private/tmp/su2-mms-6b2-complete/euler_mms_evolution --output-dir /private/tmp/su2-mms-6b2-complete/campaign --collect-only
/private/tmp/su2-phase5-tests/test_driver --reporter compact
python3 -m py_compile UnitTests/run_euler_mms_evolution.py UnitTests/export_euler_mms_evolution.py
git diff --check
```

Build del banco riuscita contro gli archivi nativi. Regressioni: 70 casi,
22 253 237 asserzioni superate. I file e risultati 6A.1/6B.1 sono conservati;
la modifica preesistente a `UnitTests/meson.build` non è alterata da questa fase.
Non è stato realizzato un collegamento al driver generale, né output SU2
di produzione, né un caso shock/vortice/contatto. Restano da revisionare
il regime asintotico, l'effetto del bordo forte/debole e il significato dei DOF;
non vengono implementate ulteriori correzioni o fasi in questa consegna.
