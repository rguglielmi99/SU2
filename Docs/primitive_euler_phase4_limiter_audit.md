# FASE 4 — Audit dei limiter sull’operatore affine reale

## Ambito

Branch `feature/primitive-path-solver`, HEAD `3df2a51577`; il working tree
contieneva già l'audit precedente non committato. Le sue modifiche sono
preservate. L'esperimento separato rimane a `d1732e90e3`.
Questo ciclo non modifica solver, kernel, geometria o configurazione.
Non introduce correzioni dei due problemi precedenti, nuove ricostruzioni,
Hessiane produttive o quadratura segmentata. Nessuna FASE 5, commit/push/sync.

La configurazione è NODE / AGGREGATED / NODE_POINT / MUSCL YES, con gradienti
LEAST_SQUARES o WEIGHTED_LEAST_SQUARES e SLOPE_LIMITER_FLOW uguale a NONE,
VENKATAKRISHNAN o BARTH_JESPERSEN. Sono i nomi effettivamente supportati
nel wrapper computeLimiters e nella validazione del solver primitivo.
Tutte le misure usano direttamente ComputeSpatialResidual.

## Implementazione SU2 effettiva

Riferimenti locali: SU2_CFD/include/limiters/computeLimiters.hpp,
computeLimiters_impl.hpp e CLimiterDetails.hpp; PrepareReconstruction in
CPrimitiveEulerSolver.cpp; CConfig.cpp per i valori di riferimento.

Il wrapper NONE assegna 1 a ogni componente. Per gli altri due algoritmi,
il nodo i raccoglie min/max di W nei vicini diretti, includendo W_i, e i
massimi/minimi delle proiezioni del gradiente sul punto medio primale:

```text
p_ij,k = (1/2) G_i,k · (x_j-x_i),  k= rho,u,v,p
p+ = max(0,p_ij), p- = min(0,p_ij)
delta+ = max(W_i,W_j)-W_i, delta- = min(W_i,W_j)-W_i
f(p,delta,e) = (delta² + 2 delta p + e)/(delta² + delta p + 2p² + e)
phi_i,k = min(f(p+,delta+,e), f(p-,delta-,e))
```

Il kappa passato dal solver è zero: umusclProjection non cambia questa
proiezione. Il fattore geometrico è 1 per entrambi i limiter. Non si usano
stati ghost analitici della MMS per determinare min/max o gradienti.

- BARTH_JESPERSEN usa **questa funzione razionale**, con
  e=epsilon(passivedouble), normalmente 2.220446049250313e-16.
  Non è la consueta funzione min(1,delta/p) del Barth–Jespersen classico.
  È la specializzazione presente nel repository, esplicitamente commentata
  come funzione di Venkatakrishnan con epsilon numerico.
- VENKATAKRISHNAN usa la medesima funzione con
  e=max(|K L|³,epsilon), K=VENKAT_LIMITER_COEFF, default 0.05.
  L=GetRefElemLength è impostato a 1 nel sistema metrico e 1/0.3048 in US;
  **non dipende dal volume locale né da h**. In questi test e=0.000125.
  Non è il parametro REF_LENGTH dei momenti aerodinamici e non esiste qui
  un automatismo che sostituisce L con la dimensione della cella.

Questi valori non vengono ritarati fra N=8 e 128 o fra componenti primitive.
Una dipendenza dall'ampiezza/dimensionalizzazione delle variabili è quindi
parte della scelta attuale di regolarizzazione. Non si importa una teoria
basata su e~h³ per giustificare questa implementazione a e costante.

La funzione non ha un clamp esplicito a 1: per r=delta/p>2 ed e piccolo,
f=(r²+2r)/(r²+r+2)>1. Anche il minimo delle due direzioni può superare 1.
L'audit misura separatamente phi<1 e phi>1; non fallisce attribuendo al codice
un contratto [0,1] che non possiede. I valori devono essere finiti e >=0.
Non emerge da questo fatto, da solo, un errore implementativo da correggere.

Il limiter è uno scalare **per primitiva e nodo** e moltiplica entrambe le
direzioni del gradiente di quella primitiva. Non è un singolo coefficiente
comune alle quattro primitive, né una limitazione indipendente della singola
traccia. In seguito il solver calcola theta_i per rho e p sulle ricostruzioni
usate e contrae tutte le componenti insieme:

```text
G_limited_i,k = theta_i * phi_i,k * G_raw_i,k
```

Prima si applica il limiter di pendenza, poi la positività. Le fluttuazioni,
il volume NODE_POINT e le mezze facce fisiche usano lo stesso G_limited.
W_i resta invariato; il limiter non applica una seconda trasformazione al
risultato del kernel. Il commit transazionale del residuale resta intatto.

## Esperimenti e diagnostiche

La MMS originale W, G analitico, sorgente B(W,G) e dati fisici esatti sono
immutati. Famiglie regolare e deformazione liscia fissa già usate nell'audit:
x'=x+0.035 sin(2pi x)sin(2pi y),
y'=y+0.03 sin(2pi x)sin(pi y).
N=8,16,32,64,128, con sei gruppi originali:
all, boundary, ring1, ring2, deep3 (strati >=3) e corners.
La norma è sqrt(sum V_i ||rhs_i||² / sum V_i), senza escludere componenti.
Gli ordini sono log2(E_N/E_2N); si esportano tutte le transizioni.

Per ogni gruppo e ciascuna delle quattro primitive si registrano phi minimo,
percentili 10/50/90 (quantile empirico con indice floor(p*(n-1))), media |1-phi|,
frazione phi<1-1e-8, phi<0.99, <0.9, <0.5 e phi>1+1e-8.
Le frazioni sono contate per nodo, non pesate sui volumi. Il quinto canale
(component=4) usa min_k phi_i,k e misura se almeno una componente è limitata;
non è una quinta incognita. Si registrano anche frazione theta<1-1e-12 e minimo
theta. Le norme degli errori di G_limited contro G analitico aiutano a separare
effetto del gradiente da ordine del residuale completo.

Il test complementare degli estremi usa un campo distinto dalla MMS originale:
W_k=b_k+0.1(k+1)[(x-0.5)²+(y-0.5)²], b=(2,0.5,0.3,3).
È fisicamente ammissibile, con minimo interno comune e gradiente analitico.
La stessa sorgente B(W,G) e il bordo esatto rendono il residuale esatto nullo.
Una norma locale aggiuntiva copre i 25 nodi topologici a distanza <=2 dal
centro (5x5); non cambia nessuna definizione regionale precedente.
Il CSV locale riporta anche phi al centro e il minimo nelle sue vicinanze.
Sul regolare NONE riproduce il residuo quadratico nell'interno fino a roundoff:
gli ordini calcolati su errori ~1e-14 non vanno interpretati come convergenza.

Il test affine usa gradienti costanti non banali e sorgente analitica su
regolare e jitter, LS/WLS e tutti i limiter. Verifica phi=1 entro 2e-11,
G esatto entro 2e-12, theta=1 e rhs nullo entro 2e-11. Per un campo affine,
delta±=2p± e f=1: l'assenza di limitazione è matematicamente prevista.
Il caso p=delta=0 restituisce 1 grazie a e>0.
Lo stato costante ha residuale nullo con tutte le scelte.

Un terzo test confronta patch con dati rho/p ripidi, rinumerazione non banale,
ordine degli elementi invertito e connettività fisica invertita dopo il
preprocessing SU2. Confronta phi, theta, gradienti e rhs; verifica stati
nodali immutati e rho/p positive su tutte le facce interne/fisiche usate.
Non si pretende monotonicità, TVD o positività di un passo temporale inesistente.

## Analisi dell’errore e interpretazione

Pongo epsilon_i,k=G_limited_i,k-G_esatto_i,k. Nella MMS theta=1; quindi

```text
epsilon_i,k = (G_raw_i,k-G_esatto_i,k) + (phi_i,k-1) G_raw_i,k.
```

Nel volume NODE_POINT la perturbazione è esattamente
A_x(W_i)epsilon_x+A_y(W_i)epsilon_y, poiché W_i è il dato analitico esatto.
Al punto medio m_ij delle facce interne il salto è

```text
jump_ij(m) = O(h³) - 0.5 (epsilon_i+epsilon_j) (x_j-x_i).
```

Con epsilon=O(h²), il salto è O(h³), la fluttuazione integrata O(h⁴) e il
residuale normalizzato O(h²). Questa è una condizione sufficiente semplice
nell'interno profondo, purché stati, stencil e velocità caratteristiche
rimangano uniformemente regolari. Non è una condizione necessaria: errori
coerenti di volume/facce possono cancellarsi in famiglie simmetriche.
Non si inferisce l'ordine del residuale dalla sola norma del gradiente.

Se phi resta discosto da 1 con intensità O(1) dove G è O(1), in generale
non si conserva questa accuratezza. Se phi→1 ma 1-phi=O(h), il gradiente
può restare soltanto O(h) e le cancellazioni dipendono dalla regolarità
spaziale di phi e dalla geometria. Una frazione attiva vicina al 100% non
implica da sola forte limitazione: servono intensità e tassi di decadimento.
Vicino a un estremo, G può essere O(h) o zero: limitazione O(1) non significa
necessariamente errore O(1) del gradiente, ma può dare errore locale O(h).
Una zona che si restringe può poi pesare diversamente nella norma globale.

Su una faccia fisica q, l'errore della traccia affine include O(h²) naturale
più epsilon_i(q-x_i). La normale è O(h) e V~h²: resta un possibile errore
O(h) anche con gradiente esatto. La limitazione può aggiungervi un difetto
maggiore. I dati del bordo sono quindi confrontati con NONE prima di
attribuire una perdita d'ordine al limiter.

Con la funzione SU2, per un campo liscio lontano dagli estremi,
delta=2p+O(h²), p=O(h). Con epsilon macchina trascurabile, lo sviluppo di f
attorno a delta/p=2 dà phi-1=O(h), non automaticamente O(h²).
Con regolarizzazione fissa e>0, invece:

```text
1-f = p (2p-delta)/(delta²+delta p+2p²+e).
```

Quando e domina, il numeratore liscio è O(h³), perciò 1-phi può decadere
come h³. Vicino al minimo quadratico p,delta=O(h²), il numeratore è O(h⁴).
Questo spiega perché VENKATAKRISHNAN tende a disattivarsi nei test raffinati;
non garantisce la stessa transizione per un'ampiezza o nondimensionalizzazione
differente. Non sono modificati K o L per ottenere i risultati.

### Particolarità del bordo con BARTH_JESPERSEN

In una direzione in cui il nodo fisico è estremo sui soli vicini, delta-=0
(o delta+=0). Un piccolo errore monolaterale del gradiente può tuttavia dare
p-!=0. Il coefficiente di quella direzione diventa e/(2p-²+e). Con epsilon
macchina esso può essere quasi zero: il minimo delle direzioni spegne
**tutto** il gradiente della componente, compresa la sua derivata tangenziale
non piccola. Nell'intervallo N<=128 la soglia macchina non regolarizza
sufficientemente questa situazione. Se h tende a zero molto oltre queste
risoluzioni, e costante può infine dominare: non si attribuisce ai tassi
misurati una prova di limite asintotico universale.

Le forti riduzioni al bordo della MMS BARTH_JESPERSEN sono quindi
un'interazione fra stencil monolaterale, estremi sui vicini e funzione
regolarizzata quasi senza epsilon, non un segno errato nel kernel. VENKAT
con e=1.25e-4 attenua questo meccanismo. Un vero trattamento alternativo
richiederebbe una scelta numerica esplicita sui dati di ricostruzione al
bordo; questo audit non introduce ghost states, clamp o nuovi parametri.

La funzione razionale e la possibile amplificazione phi>1 sono implementate
esplicitamente nel codice SU2, non un difetto dimostrato del solver primitivo.
Non vi è quindi una correzione di produzione obbligata da applicare in questo
ciclo. Il comportamento osservato viene documentato anziché sostituito con
la teoria di un diverso limiter chiamato con lo stesso nome.

## Risultati quantitativi

Le tabelle seguenti sono estratte dai CSV finali senza selezionare livelli.
I tassi su errori di roundoff del campo quadratico NONE non sono significativi.

### MMS originale: norma globale

| Famiglia | G | Limiter | L2 N=8 | p 8→16 | p 16→32 | p 32→64 | p 64→128 | L2 N=128 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| regular | LS | NONE | 1.455626e-02 | 1.49558 | 1.49752 | 1.49857 | 1.49922 | 2.288827e-04 |
| regular | LS | VENKATAKRISHNAN | 1.404336e-02 | 1.46806 | 1.48395 | 1.49143 | 1.49666 | 2.287289e-04 |
| regular | LS | BARTH_JESPERSEN | 3.750020e-02 | 0.74467 | 0.62576 | 0.56974 | 0.56664 | 6.597932e-03 |
| regular | WLS | NONE | 1.317217e-02 | 1.50248 | 1.50090 | 1.50024 | 1.50005 | 2.052919e-04 |
| regular | WLS | VENKATAKRISHNAN | 1.275548e-02 | 1.47822 | 1.48888 | 1.49355 | 1.49759 | 2.051575e-04 |
| regular | WLS | BARTH_JESPERSEN | 3.693799e-02 | 0.73419 | 0.61997 | 0.57196 | 0.57363 | 6.530941e-03 |
| smooth_distorted | LS | NONE | 1.527705e-02 | 1.50411 | 1.50047 | 1.49889 | 1.49906 | 2.382838e-04 |
| smooth_distorted | LS | VENKATAKRISHNAN | 1.455903e-02 | 1.46854 | 1.48033 | 1.48960 | 1.49583 | 2.380832e-04 |
| smooth_distorted | LS | BARTH_JESPERSEN | 3.836373e-02 | 0.74446 | 0.62156 | 0.56670 | 0.56360 | 6.799071e-03 |
| smooth_distorted | WLS | NONE | 1.346314e-02 | 1.49822 | 1.49993 | 1.49998 | 1.49996 | 2.106399e-04 |
| smooth_distorted | WLS | VENKATAKRISHNAN | 1.291992e-02 | 1.46824 | 1.48314 | 1.49155 | 1.49692 | 2.104697e-04 |
| smooth_distorted | WLS | BARTH_JESPERSEN | 3.760945e-02 | 0.72941 | 0.61453 | 0.56849 | 0.57588 | 6.702582e-03 |

### MMS originale: interno profondo

| Famiglia | G | Limiter | L2 N=8 | p 8→16 | p 16→32 | p 32→64 | p 64→128 | L2 N=128 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| regular | LS | NONE | 7.438037e-04 | 2.03030 | 1.93747 | 1.95815 | 1.97765 | 3.106319e-06 |
| regular | LS | VENKATAKRISHNAN | 6.874457e-04 | 1.31368 | 2.19276 | 2.21074 | 2.05157 | 3.152178e-06 |
| regular | LS | BARTH_JESPERSEN | 7.954172e-04 | 1.70971 | 1.72133 | 1.74806 | 1.72480 | 6.642218e-06 |
| regular | WLS | NONE | 7.102142e-04 | 1.98779 | 1.92241 | 1.95294 | 1.97556 | 3.102454e-06 |
| regular | WLS | VENKATAKRISHNAN | 6.886171e-04 | 1.32182 | 2.19288 | 2.20854 | 2.04976 | 3.148297e-06 |
| regular | WLS | BARTH_JESPERSEN | 8.049064e-04 | 1.72513 | 1.72457 | 1.74809 | 1.72430 | 6.637266e-06 |
| smooth_distorted | LS | NONE | 3.173341e-03 | 2.79747 | 2.50978 | 2.30452 | 2.13266 | 3.699572e-06 |
| smooth_distorted | LS | VENKATAKRISHNAN | 3.146837e-03 | 2.80136 | 2.35387 | 2.37591 | 2.18221 | 3.748562e-06 |
| smooth_distorted | LS | BARTH_JESPERSEN | 3.312041e-03 | 2.81916 | 2.19567 | 1.74469 | 1.66056 | 9.669404e-06 |
| smooth_distorted | WLS | NONE | 8.899947e-04 | 2.18149 | 1.87180 | 1.95719 | 1.98493 | 3.487615e-06 |
| smooth_distorted | WLS | VENKATAKRISHNAN | 1.104722e-03 | 2.05459 | 1.99238 | 2.18302 | 2.05656 | 3.537974e-06 |
| smooth_distorted | WLS | BARTH_JESPERSEN | 1.337692e-03 | 2.07994 | 1.76655 | 1.64148 | 1.63432 | 9.601349e-06 |

### MMS originale: regioni fisiche N=128, errore / p 64→128

| Famiglia | G | Limiter | Bordo | Ring1 | Ring2 |
| --- | --- | --- | ---: | ---: | ---: |
| regular | LS | NONE | 1.4707e-03 / 1.002 | 7.8000e-04 / 1.010 | 4.0630e-06 / 1.982 |
| regular | LS | VENKATAKRISHNAN | 1.4687e-03 / 0.997 | 7.8044e-04 / 1.012 | 4.0698e-06 / 2.006 |
| regular | LS | BARTH_JESPERSEN | 4.3193e-02 / 0.079 | 2.1708e-02 / 0.059 | 9.1415e-06 / 1.614 |
| regular | WLS | NONE | 1.3152e-03 / 1.003 | 7.0321e-04 / 1.011 | 4.0585e-06 / 1.980 |
| regular | WLS | VENKATAKRISHNAN | 1.3134e-03 / 0.997 | 7.0362e-04 / 1.014 | 4.0653e-06 / 2.004 |
| regular | WLS | BARTH_JESPERSEN | 4.2753e-02 / 0.085 | 2.1489e-02 / 0.068 | 9.1374e-06 / 1.613 |
| smooth_distorted | LS | NONE | 1.5306e-03 / 1.000 | 8.1247e-04 / 1.013 | 5.1907e-06 / 2.160 |
| smooth_distorted | LS | VENKATAKRISHNAN | 1.5281e-03 / 0.994 | 8.1288e-04 / 1.015 | 5.2362e-06 / 2.205 |
| smooth_distorted | LS | BARTH_JESPERSEN | 4.4507e-02 / 0.075 | 2.2372e-02 / 0.057 | 1.4214e-05 / 1.555 |
| smooth_distorted | WLS | NONE | 1.3495e-03 / 1.002 | 7.2149e-04 / 1.013 | 4.8584e-06 / 2.035 |
| smooth_distorted | WLS | VENKATAKRISHNAN | 1.3474e-03 / 0.995 | 7.2190e-04 / 1.016 | 4.9004e-06 / 2.087 |
| smooth_distorted | WLS | BARTH_JESPERSEN | 4.3874e-02 / 0.087 | 2.2057e-02 / 0.071 | 1.4066e-05 / 1.522 |

### Attivazione MMS WLS, N=128, canale min delle quattro primitive

Frazione attiva usa phi<1-1e-8; la colonna <0.99 distingue interventi più intensi.

| Famiglia | Limiter | Regione | Frazione attiva | Min phi | p10 / p50 / p90 | Frazione <0.99 | Media abs(1-phi) |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: |
| regular | NONE | boundary | 0.00000 | 1.0000000 | 1.0000000 / 1.0000000 / 1.0000000 | 0.00000 | 0.00000e+00 |
| regular | NONE | deep3 | 0.00000 | 1.0000000 | 1.0000000 / 1.0000000 / 1.0000000 | 0.00000 | 0.00000e+00 |
| regular | VENKATAKRISHNAN | boundary | 0.99414 | 0.9999633 | 0.9999675 / 0.9999940 / 0.9999975 | 0.00000 | 1.17299e-05 |
| regular | VENKATAKRISHNAN | deep3 | 1.00000 | 0.9999500 | 0.9999546 / 0.9999657 / 0.9999908 | 0.00000 | 3.06571e-05 |
| regular | BARTH_JESPERSEN | boundary | 0.99609 | 0.0003684 | 0.0007705 / 0.7823490 / 0.9991746 | 0.54102 | 4.88640e-01 |
| regular | BARTH_JESPERSEN | deep3 | 1.00000 | 0.9483196 | 0.9959352 / 0.9986587 / 0.9992939 | 0.01553 | 2.09124e-03 |
| smooth_distorted | NONE | boundary | 0.00000 | 1.0000000 | 1.0000000 / 1.0000000 / 1.0000000 | 0.00000 | 0.00000e+00 |
| smooth_distorted | NONE | deep3 | 0.00000 | 1.0000000 | 1.0000000 / 1.0000000 / 1.0000000 | 0.00000 | 0.00000e+00 |
| smooth_distorted | VENKATAKRISHNAN | boundary | 0.99609 | 0.9999613 | 0.9999665 / 0.9999930 / 0.9999973 | 0.00000 | 1.21292e-05 |
| smooth_distorted | VENKATAKRISHNAN | deep3 | 1.00000 | 0.9999301 | 0.9999507 / 0.9999715 / 0.9999918 | 0.00000 | 2.97766e-05 |
| smooth_distorted | BARTH_JESPERSEN | boundary | 0.99609 | 0.0003743 | 0.0007164 / 0.7821248 / 0.9991620 | 0.53906 | 4.88247e-01 |
| smooth_distorted | BARTH_JESPERSEN | deep3 | 1.00000 | 0.9473318 | 0.9960004 / 0.9985704 / 0.9992048 | 0.01785 | 2.12041e-03 |

### Differenze fra primitive, MMS regolare WLS N=128

| Limiter | Primitiva | Regione | Min phi | Frazione attiva | Media abs(1-phi) |
| --- | --- | --- | ---: | ---: | ---: |
| VENKATAKRISHNAN | rho | boundary | 0.9999909 | 0.64062 | 3.39599e-06 |
| VENKATAKRISHNAN | u | boundary | 0.9999970 | 0.74414 | 8.35144e-07 |
| VENKATAKRISHNAN | v | boundary | 0.9999986 | 0.72070 | 5.83999e-07 |
| VENKATAKRISHNAN | p | boundary | 0.9999633 | 0.64648 | 1.14376e-05 |
| VENKATAKRISHNAN | rho | deep3 | 0.9999894 | 1.00000 | 5.65475e-06 |
| VENKATAKRISHNAN | u | deep3 | 0.9999958 | 0.99663 | 1.00155e-06 |
| VENKATAKRISHNAN | v | deep3 | 0.9999980 | 1.00000 | 7.91435e-07 |
| VENKATAKRISHNAN | p | deep3 | 0.9999500 | 1.00000 | 3.05084e-05 |
| BARTH_JESPERSEN | rho | boundary | 0.0003684 | 0.89062 | 2.45551e-01 |
| BARTH_JESPERSEN | u | boundary | 0.0007613 | 0.99609 | 2.43804e-01 |
| BARTH_JESPERSEN | v | boundary | 0.0014882 | 0.97656 | 2.39478e-01 |
| BARTH_JESPERSEN | p | boundary | 0.7196365 | 0.65039 | 3.83857e-03 |
| BARTH_JESPERSEN | rho | deep3 | 0.9985442 | 1.00000 | 5.74624e-04 |
| BARTH_JESPERSEN | u | deep3 | 0.9980796 | 1.00000 | 3.57866e-04 |
| BARTH_JESPERSEN | v | deep3 | 0.9985332 | 1.00000 | 3.77602e-04 |
| BARTH_JESPERSEN | p | deep3 | 0.9483196 | 1.00000 | 1.89335e-03 |

### Minimo liscio: campo quadratico, interno profondo

| Famiglia | G | Limiter | L2 N=128 | p 16→32 | p 32→64 | p 64→128 |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| regular | LS | NONE | 7.048758e-14 | roundoff | roundoff | roundoff |
| regular | LS | VENKATAKRISHNAN | 9.805314e-06 | 2.69361 | 2.85316 | 3.24868 |
| regular | LS | BARTH_JESPERSEN | 1.510433e-04 | 2.11089 | 1.93180 | 1.82146 |
| regular | WLS | NONE | 7.452751e-14 | roundoff | roundoff | roundoff |
| regular | WLS | VENKATAKRISHNAN | 9.805314e-06 | 2.69361 | 2.85316 | 3.24868 |
| regular | WLS | BARTH_JESPERSEN | 1.510433e-04 | 2.11089 | 1.93180 | 1.82146 |
| smooth_distorted | LS | NONE | 6.917529e-06 | 2.94438 | 2.97595 | 2.99164 |
| smooth_distorted | LS | VENKATAKRISHNAN | 1.683199e-05 | 2.56164 | 2.74624 | 3.12683 |
| smooth_distorted | LS | BARTH_JESPERSEN | 1.454057e-04 | 2.17286 | 1.95037 | 1.86363 |
| smooth_distorted | WLS | NONE | 2.242676e-06 | 2.96816 | 2.98362 | 2.99424 |
| smooth_distorted | WLS | VENKATAKRISHNAN | 1.428593e-05 | 2.53456 | 2.72690 | 3.17474 |
| smooth_distorted | WLS | BARTH_JESPERSEN | 1.449181e-04 | 2.14601 | 1.92810 | 1.85657 |

### Vicinanze del minimo: 25 nodi, WLS N=128

| Famiglia | Limiter | L2 locale | p 64→128 | Min phi locale | Phi al centro (rho/u/v/p) |
| --- | --- | ---: | ---: | ---: | --- |
| regular | NONE | 8.843515e-14 | roundoff | 1.00000000 | 1.000000 / 1.000000 / 1.000000 / 1.000000 |
| regular | VENKATAKRISHNAN | 4.870603e-07 | 4.98762 | 0.99996187 | 1.000000 / 1.000000 / 1.000000 / 1.000000 |
| regular | BARTH_JESPERSEN | 2.161911e-03 | 1.00016 | 0.75000001 | 1.000000 / 1.000000 / 1.000000 / 1.000000 |
| smooth_distorted | NONE | 1.112322e-06 | 3.27457 | 1.00000000 | 1.000000 / 1.000000 / 1.000000 / 1.000000 |
| smooth_distorted | VENKATAKRISHNAN | 8.560838e-07 | 2.93105 | 0.99997167 | 1.000000 / 1.000000 / 1.000000 / 1.000000 |
| smooth_distorted | BARTH_JESPERSEN | 1.981258e-03 | 1.00366 | 0.01317977 | 0.176068 / 0.050714 / 0.023193 / 0.013180 |

### Dati ripidi: positività e rinumerazione

| G | Limiter | Frazione theta attivo | Min theta | Min phi |
| --- | --- | ---: | ---: | ---: |
| LEAST_SQUARES | NONE | 0.44000 | 5.5951102e-02 | 1.0000000e+00 |
| LEAST_SQUARES | VENKATAKRISHNAN | 0.00000 | 1.0000000e+00 | 1.9963047e-07 |
| LEAST_SQUARES | BARTH_JESPERSEN | 0.00000 | 1.0000000e+00 | 3.5461502e-19 |
| WEIGHTED_LEAST_SQUARES | NONE | 0.44000 | 3.7034607e-02 | 1.0000000e+00 |
| WEIGHTED_LEAST_SQUARES | VENKATAKRISHNAN | 0.00000 | 1.0000000e+00 | 8.7463140e-08 |
| WEIGHTED_LEAST_SQUARES | BARTH_JESPERSEN | 0.00000 | 1.0000000e+00 | 1.5536576e-19 |

## Conclusioni e configurazione consigliata

NONE riproduce il riferimento: MMS regolare WLS p_deep=1.97556,
p_all=1.50005; sulla mappa liscia WLS p_deep=1.98493. I nuovi valori coincidono
con gli audit precedenti entro roundoff (le diverse unità di compilazione
possono cambiare le ultime cifre). Non si attribuisce la perdita globale del
baseline al limiter, che è identicamente 1.

VENKATAKRISHNAN mantiene un interno compatibile con ordine 2 nelle due
famiglie misurate. L'ultima transizione dà LS/WLS 2.05157/2.04976 sul regolare
e 2.18221/2.05656 sulla mappa liscia. Non ogni transizione è già asintotica:
i tassi completi sopra consentono di vedere il regime di regolarizzazione.
Il globale regolare WLS rimane circa 1.49759 per il bordo, come NONE.
Al N=128 regolare WLS, tutti i nodi deep3 hanno almeno una componente phi
inferiore a 1-1e-8, ma nessuna inferiore a 0.99; il minimo è 0.99995002.
“Attivo su tutti i nodi” qui non significa forte degradazione.

BARTH_JESPERSEN dà invece p_deep finale circa 1.72 regolare e 1.63–1.66
sulla mappa liscia per questa MMS. Con WLS, phi_min deep3=0.94831959 sul
regolare; al bordo scende a 0.00036837 e quasi metà dei nodi ha min_k phi<0.5.
Il residuale di bordo WLS è 0.0427528 con p=0.08494 nell'ultima transizione,
contro circa 0.001315 e p~1 con NONE. Il globale passa a p=0.57363.
Questo è un difetto aggiuntivo del percorso limitato ai livelli misurati,
non il solo errore fisico O(h) già noto. Non si certifica un esponente
asintotico universale oltre N=128 o un limite a epsilon macchina dominato.

Nel campo quadratico con minimo interno, il valore di phi esattamente al
centro non descrive le sue vicinanze. Sul regolare WLS tutti i phi centrali
sono 1, ma nei 25 nodi vicini il minimo BARTH è 0.75000001 e p_locale=1.00016;
VENKAT ha phi_min=0.99996187 e p_locale=4.98762. Quest'ultimo tasso elevato
è specifico della quadratica e della regolarizzazione fissa, non un quinto
ordine generale del metodo. Sulla mappa liscia, al centro dove il gradiente
esatto è zero, piccoli errori del gradiente nativo possono far scendere
BARTH sino a 0.01318 nel p; nonostante questo, l'errore locale resta circa
O(h), perché il gradiente stesso è piccolo. NONE sul regolare ha errori
~1e-14, quindi non si interpreta il suo ordine numerico negativo.

La funzione può amplificare leggermente il gradiente (phi>1): 737 delle
3600 righe aggregate di attività hanno una frazione non nulla in quel canale.
Le righe sono gruppi sovrapposti e componenti: **non** sono una frazione di
nodi globali. I CSV per componente permettono di valutare l'effetto senza
confonderlo con la frazione di limitazione. Theta=1 su tutte le 120
configurazioni MMS/estremi: la positività non causa gli ordini osservati.
Nei dati ripidi NONE richiede theta<1 sul 44% dei nodi; i due limiter riducono
i gradienti abbastanza da lasciare theta=1 nella patch misurata.
Questo non sostituisce la positività indipendente del solver.

**Configurazione candidata più appropriata per il principale:**

```ini
PRIMITIVE_REPRESENTATION=NODE
PRIMITIVE_FACE_QUADRATURE=AGGREGATED
PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT
MUSCL_FLOW=YES
NUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES
SLOPE_LIMITER_FLOW=VENKATAKRISHNAN
VENKAT_LIMITER_COEFF=0.05
```

Non vengono cambiati i default del repository. NONE rimane il riferimento
per accuratezza liscia; BARTH_JESPERSEN non è raccomandato come default di
accuratezza in questa configurazione, considerati bordo ed estremi osservati.
WLS+VENKAT è una scelta motivata dai dati del residuale, non una prova che sia
TVD, monotona o ottimale per tutti gli shock. LS+VENKAT è anch'essa compatibile
con il secondo ordine interno sui casi misurati; WLS resta coerente con la
scelta semplice proposta negli audit precedenti.

Non emerge una correzione minima obbligata per un errore implementativo.
Cambiare la formula BARTH, limitare phi a 1 o rendere L locale sarebbe una
nuova scelta numerica che richiede autorizzazione e test propri: non viene
fatta per forzare gli ordini. La specifica qui chiede di osservare l'operatore
attuale, non di sostituire i limiter SU2.

Per la FASE 5 restano da verificare: CFL e stabilità del passo temporale,
positività dopo aggiornamento, shock/contatti e oscillazioni, dipendenza da
ampiezza/nondimensionalizzazione della regolarizzazione fissa, trattamento
delle condizioni fisiche, convergenza della soluzione e interazione con
l'eventuale gestione iterativa/freeze del limiter. ComputeSpatialResidual
ricalcola il limiter a ogni chiamata; non esiste ancora un driver temporale
che ne gestisca l'attivazione. I limiti già accettati di MUSCL NO e jitter
non sono risolti da questo audit. Nessun secondo ordine globale di bordo
è dichiarato.

## Riproduzione e verifiche

Sono stati riutilizzati gli argv del runner dedicato per ricompilare il solo
translation unit primitive_euler_accuracy_tests.cpp e ricollegare il driver;
log `/private/tmp/su2-limiter-build.log`. Per una ricostruzione autonoma:

```sh
./ninja -C build -j 2 SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-limiter-reproduce
/private/tmp/su2-limiter-reproduce/test_driver '[primitive-limiter-audit]'
/private/tmp/su2-limiter-reproduce/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-limiter-cli --cli
python3 UnitTests/export_primitive_euler_limiters.py --log /private/tmp/su2-limiter-results.log
```

I run effettivi usano `/private/tmp/su2-affine-archive-tests/test_driver` con
le selezioni sopra e redirezione dei log:
`/private/tmp/su2-limiter-diagnostics.log` per il primo run di due test,
`/private/tmp/su2-limiter-results.log` per la suite finale. Dopo l'aggiunta
del test di rinumerazione e della norma locale del minimo è stato ricompilato
il translation unit e rieseguita l'intera suite.

Esito finale: **64 casi, 22.251.834 asserzioni superate**, uscita 0. Include
tutti i 61 casi degli audit precedenti e i tre nuovi casi di limiter.
Non sono indebolite soglie o risoluzioni preesistenti. Build nativa: uscita 0,
“no work to do”, coerentemente con produzione invariata. Quattro controlli
CLI riusciti con uscite 1/0/1/1: primitivo evolutivo bloccato, Euler conservativo
riuscito, harmonic balance e specie rifiutati. Log CLI e build:
`/private/tmp/su2-limiter-cli.log`, `/private/tmp/su2-limiter-native-build.log`.
Warning noti non bloccanti: pragma GCC ignorati da Clang, librerie duplicate,
WINDOW_START_ITER dei vecchi test e bind MPI nel sandbox. Non è stata
rieseguita una build ASan o l'intera suite generale SU2.

File aggiunti/modificati in questo ciclo, oltre a preservare quelli già dirty:

- UnitTests/SU2_CFD/primitive_euler_accuracy_tests.cpp: tre casi aggiunti;
- UnitTests/export_primitive_euler_limiters.py: esportazione autonoma;
- questo rapporto e un collegamento nel documento principale della FASE 4;
- [CSV MMS](primitive_euler_phase4_limiter_mms.csv), 720 righe dati;
- [CSV attività](primitive_euler_phase4_limiter_activity.csv), 3600 righe;
- [CSV estremi locali](primitive_euler_phase4_limiter_extrema.csv), 60 righe;
- [CSV dati ripidi](primitive_euler_phase4_limiter_stress.csv), 6 righe.

Il CSV MMS include anche l'errore del gradiente limitato; quello di attività
include LS e WLS, le quattro primitive e tutti i gruppi a ogni livello,
anche quando una tabella compatta sopra mostra soltanto WLS N=128.
I CSV e le MMS degli audit precedenti non sono sovrascritti.
Gli hash dei file produttivi iniziali sono invariati e il diff di Common,
SU2_CFD e config_template.cfg rimane vuoto. Entrambi i branch mantengono
i propri HEAD iniziali; nessun commit, push, sync o FASE 5.
