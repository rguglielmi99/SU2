# FASE 4 — Passaggio I: consistenza senza MUSCL e accuratezza su mesh irregolari

## Ambito e stato di riferimento

Branch `feature/primitive-path-solver`, HEAD iniziale `3df2a51577`, working tree
inizialmente pulito. Il branch sperimentale rimane a `d1732e90e3`.
Questo ciclo aggiunge esclusivamente test e rapporti: solver, geometria,
kernel, CConfig e registrazioni Meson non sono modificati. Nessun commit,
push, sync, FASE 5 o correzione di produzione.

Il riferimento è stato rieseguito prima degli esperimenti: 58 casi,
1.380.527 asserzioni superate. Sulla MMS NODE_POINT, N=128, interno profondo:

| Famiglia | LS: L2 / p 64→128 | WLS: L2 / p 64→128 |
| --- | ---: | ---: |
| Regolare | 3.106319e-6 / 1.97765 | 3.102454e-6 / 1.97556 |
| Distorta a scala h | 1.077514e-3 / 1.00518 | 4.590496e-4 / 1.00064 |

Il solver prepara LS/WLS (oppure GG nel percorso preesistente), applica il
limiter e un unico theta di positività. MUSCL NO azzera G. Ogni faccia interna
usa Q=(x_i+x_j)/2 e la normale SU2 integrata S_ij; i due scatter sono entrambi
positivi. Il volume è V_i B(W_i,G_i) in NODE_POINT, con
B(W,G)=A_x(W)G_x+A_y(W)G_y. Le mezze facce fisiche usano la traccia affine
del proprietario e l'esterno analitico, aggiungendo soltanto Dminus.
La sorgente è S(x_i), senza volume; rhs=S-(faces+volume+boundary)/V_i.
Il commit di LinSysRes avviene soltanto alla fine. Questa pipeline non cambia.

Tutti i nuovi confronti usano NODE / NODE_POINT, limiter NONE e gamma=1.4.
Il fit diagnostico campiona esclusivamente W numerico e restituisce soltanto G:
non fornisce Hessiane/tracce quadratiche/cache al solver reale.

## A. Operatore completo senza MUSCL

Sia W(x)=W_i+G(x-x_i) e d_ij=x_j-x_i. Orientando tutte le normali localmente
verso l'esterno del duale i, l'operatore esatto implementato è:

```text
R_i^int = (1/2) sum_j [Atilde_ij(S_ij)-|Atilde_ij(S_ij)|] G d_ij
R_i^vol = 0
R_i^bdry = sum_b Dminus(W_i, W_ext(q_b), S_ib)
rhs_i = S(x_i) - (R_i^int+R_i^bdry)/V_i
```

L'identità per il secondo proprietario segue da
Dplus(W_i,W_j,S)=Dminus(W_j,W_i,-S). Non si cambia alcun segno.
Per l'esterno affine esatto, W_ext(q_b)=W_i+G(q_b-x_i).
Si separano C_i=(1/2)sum Atilde*jump e Q_i=-(1/2)sum |Atilde|*jump,
comprendendo le mezze facce fisiche nelle rispettive somme.

Per un nodo interno, congelando le matrici in W_i:

```text
M_i(alpha,beta) = (1/(2 V_i)) sum_j S_ij,alpha d_ij,beta
R_i/V_i = sum_(alpha,beta) A_alpha(W_i) G_beta M_i(alpha,beta)
          - (1/(2 V_i)) sum_j |A(S_ij;W_i)| G d_ij + O(h)
```

A(S) è omogenea di grado uno nella normale integrata. Le variazioni delle
matrici lungo il percorso contribuiscono O(h) dopo divisione per V;
nel contatto descritto sotto la formula scalare è invece esatta.
Una condizione necessaria di consistenza puntuale per tutti i campi affini
è l'identità combinata fra questi momenti e il termine differenziale;
M_i=I e annullamento del momento dissipativo sono condizioni sufficienti
per il termine principale a matrici congelate. Lo stato costante è esatto
anche se tali condizioni affini falliscono: ogni salto è allora zero.

### Il momento centrato è corretto sui duali mediani interni

Sui triangoli mediani il bordo geometrico esatto soddisfa
sum_s S_s (q_s-x_i)^T=V_i I. La differenza tra il momento del segmento e
quello al punto medio primale è T_s=S_s(q_s-m_ij)^T.
In un triangolo incidente con x_i=0 e altri vertici a,b, i due segmenti danno

```text
sum_triangle T_s = J (b b^T - a a^T)/24
```

con J rotazione oraria di 90 gradi e normali uscenti da i. La somma sui
triangoli della stella chiusa telescopa: sum T_s=0. Ne segue
(1/2)sum_j S_ij d_ij^T=V_i I, anche su una stella asimmetrica.
I test verificano tutti e quattro gli elementi di M_i per ogni patch e scala.
La chiusura della normale, da sola, non sarebbe una prova di questa identità.
Al bordo la stella non è chiusa: occorre includere esplicitamente i momenti
fisici e i mezzi operatori entranti; non si estende M_i=I interno al bordo.

### Riduzione esatta al contatto

Per rho=2+g·x, g=(0.1,0.2), a=(u,v)=(0.5,0.3), p=3:
Atilde(S)*jump=(a·S)*delta_rho nel solo componente rho e
|Atilde(S)|*jump=|a·S|*delta_rho. Pertanto

```text
C_i/V_i = a^T M_i g = 0.11
Q_i/V_i = -(1/(2 V_i)) sum_j |a·S_ij| (g·d_ij)
rhs_i = 0.11 - C_i/V_i - Q_i/V_i = -Q_i/V_i
```

La condizione dissipativa necessaria per tutti i gradienti di contatto è
sum_j |a·S_ij| d_ij=0. Coppie opposte con uguale peso la soddisfano;
una stella asimmetrica in generale no. Sotto omotetia S~h, d~h, V~h²,
il termine normalizzato resta O(1). Il problema non è un segno errato,
la sorgente, il volume NODE_POINT o la posizione del punto di quadratura:
con stati costanti per cella quel punto non entra nelle tracce.
È il momento della dissipazione upwind su stati NODE e stella asimmetrica.

| Patch, h=1 | C/V | Q/V | Operatore R/V | rhs con S=0.11 | R_segmenti/V |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fan simmetrico, 4 vicini | 0.110000 | 0 | 0.110000 | ~0 | 0.110000 |
| Fan asimmetrico, 4 vicini | 0.110000 | -0.036000 | 0.074000 | 0.036000 | 0.060500 |
| Stella asimmetrica, 6 vicini | 0.110000 | 0.00938833 | 0.11938833 | -0.00938833 | 0.11938833 |

La prima coppia ha la stessa topologia e gli stessi quattro vertici esterni;
si sposta soltanto il nodo centrale da (1,0.5) a (0.8,0.3). Il terzo caso
usa una stella interna della triangolazione distorta con sei vicini.
Per h=1,1/2,1/4,1/8 i valori normalizzati restano uguali entro arrotondamento.
Nel fan asimmetrico V*rhs vale 0.024, 0.006, 0.0015, 0.000375: è O(h²),
ma **non** o(V). Non è quindi consistenza puntuale nascosta dal volume.
La somma di segmenti reali conserva il momento centrato, ma cambia i pesi
assoluti e non annulla il difetto; nel fan lo aumenta.

Tutti i nodi delle patch, inclusi i bordi, sono confrontati con il residuale
completo: sorgente, volume nullo, centrato e dissipativo fisici/interni.
Il nodo del controesempio è strettamente interno e ha boundary=0.
Inserendo G analitico soltanto nella diagnosi MUSCL affine, le tracce di
contatto coincidono e il volume dà 0.11; rhs torna zero su tutte le patch.

### Consistenza forte, integrata e debole: conclusioni diverse

Il difetto integrato su **una cella che si restringe** tende a zero come h²,
come il volume: non basta per provare convergenza. Su un dominio fisso con
jitter a scala h, L2 profonda tende a circa 0.015 e la L1 integrata non tende
rapidamente a zero: il difetto forte non scompare nemmeno globalmente in L1.

Tuttavia sul solo contatto l'operatore ammette esattamente la riscrittura
scalare FV, usando chiusura delle normali:

```text
F*_ij = 0.5 (a·S_ij)(rho_i+rho_j) - 0.5 |a·S_ij|(rho_j-rho_i)
Dminus_ij = F*_ij - (a·S_ij) rho_i
R_i = sum_internal F*_ij + sum_physical F*_ib
```

F* è antisimmetrico quando si inverte orientamento; non è un cambiamento
del solver primitivo o una proposta di flusso conservativo Euler.
È un'identità per il sottoproblema scalare a velocità e pressione costanti.
Per un campo liscio, l'errore di flusso per faccia rispetto al flusso fisico
esatto è O(h²). Nella somma contro una funzione test liscia con traccia zero,
le facce interne si accoppiano con phi_i-phi_j=O(h). Su una famiglia
shape regular con O(h^-2) facce, la stima globale è O(h); anche l'errore della
quadratura nodale volume/sorgente è al più O(h). Questo giustifica
**consistenza debole** del sottoproblema scalare pur senza consistenza forte
locale. La stabilità e il limite della soluzione sono proprietà ulteriori.
Sul contatto il semidiscreto è monotono per dati scalari e sotto il CFL
scalare il passo Euler può preservare il massimo/minimo; la FASE 5 non esiste
qui e non si trasferisce questa proprietà all'intero sistema Euler.

Il test usa phi=256[x(1-x)y(1-y)]²: è liscia e nulla al bordo, non ha supporto
compatto strettamente interno. Non si stima un ordine universale dai suoi
forti annullamenti specifici. A N=32,64,128:

| Famiglia | L2 profonda | sum V rhs phi |
| --- | --- | --- |
| Jitter a scala h | 0.015015, 0.015017, 0.015022 | 2.6722e-7, 1.8073e-9, 1.2711e-10 |
| Mappa liscia | 1.6252e-3, 7.9010e-4, 3.8955e-4 | 8.1395e-5, 4.1164e-5, 2.0640e-5 |

La somma firmata senza phi è circa 1e-16 su queste famiglie; è un bilancio
particolare del contatto e dei bordi, non una misura dell'errore locale.
Non è dimostrata la convergenza della soluzione Euler primitiva con shock,
né può essere dedotta dalla sola MMS del residuale o dall'identità scalare.

## B. Matrice MMS: gradienti e segmenti

La soluzione liscia W e la sorgente S=B(W,G_esatto) sono quelle del precedente
audit, con derivate analitiche esplicite. Il residuale esatto è zero.
Si usa la norma sqrt(sum_gruppo V_i ||rhs_i||² / sum_gruppo V_i), per
all/boundary/ring1/deep3/ring2/corners. Deep3 esclude i primi tre strati.
Si misurano N=8,16,32,64,128; gli ordini sono log2(E_N/E_2N).
Nessun ordine viene imposto o ottenuto cambiando soglie dei vecchi test.

Tre famiglie, tutte con connettività triangolare a diagonale costante:

1. Regolare cartesiana triangolata.
2. Distorta a scala h: delta_x=0.18 h sin(2.1 i+1.3 j),
   delta_y=0.16 h cos(1.7 i-0.9 j), soltanto nei nodi non fisici.
3. Mappa liscia fissa del quadrato: x'=x+0.035 sin(2pi x)sin(2pi y),
   y'=y+0.03 sin(2pi x)sin(pi y).

Le aree orientate doppie/h² sono verificate >0.2 e le lunghezze degli edge
sono uniformemente O(h): non ci sono celle che collassano. Il jitter rimane
geometricamente regolare ma non diventa localmente simmetrico con h;
la mappa liscia ha invece variazioni dello Jacobiano O(h) fra celle vicine.

Cinque esperimenti cambiano solo gradienti e facce interne:

- native_aggregated: chiamata effettiva a ComputeSpatialResidual.
- exact_aggregated: G analitico, tracce affini al punto medio primale.
- fit_gradient_aggregated: G da LS/WLS quadratico a due anelli, tracce affini.
- exact_segmented: G analitico, integrazione sulle facce duali reali.
- native_segmented: G LS/WLS nativo, stessa integrazione geometrica.

Nel fit diagnostico il valore centrale è fissato; coordinate scalate con
ell=max distanza nello stencil, matrice [dx,dy,.5dx²,dxdy,.5dy²], QR pivotata
rango 5, LS o pesi di riga 1/distanza per WLS. Si restituiscono esclusivamente
i due coefficienti di gradiente, non quelli di curvatura. Non si riusa codice
produttivo del branch sperimentale e non esiste una nuova opzione CConfig.

Per SEGMENTED diagnostico, ogni triangolo fornisce il segmento fra punto
medio primale e baricentro del triangolo. La normale integrata è la rotazione
del segmento orientata i→j. Le normali sommano a quella SU2 entro 2e-14.
Si usano due punti Gauss reali su ciascun segmento e metà della sua normale
per ogni punto. Si ricostruiscono le stesse tracce affini in quei punti.
Volume NODE_POINT, sorgente nodale e mezze facce fisiche restano identici.
Questo è il confronto richiesto a volume invariato, non un nuovo metodo
integrato coerente con una rappresentazione DUAL_MEAN.

Il riassemblaggio con gradienti nativi deve prima riprodurre rhs produttivo
entro 2e-12. Si verifica theta=1 su ogni nodo della MMS, rank/finitezza del fit
e tutte le norme regionali. I test di positività, orientamento e transazione
precedenti restano nella suite completa.

### Ordini nell’interno profondo: tre transizioni consecutive

#### Distorta a scala h

| Metodo | Esperimento | p 16→32 | p 32→64 | p 64→128 | L2 N=128 |
| --- | --- | ---: | ---: | ---: | ---: |
| LS | native_aggregated | 1.02365 | 1.00941 | 1.00518 | 1.0775143e-03 |
| LS | exact_aggregated | 1.82270 | 1.92399 | 1.96373 | 3.6160693e-06 |
| LS | fit_gradient_aggregated | 1.95530 | 1.96792 | 1.98087 | 7.2816126e-06 |
| LS | exact_segmented | 0.99439 | 0.97916 | 0.98506 | 6.7869259e-05 |
| LS | native_segmented | 1.02505 | 1.01013 | 1.00554 | 1.1343429e-03 |
| WLS | native_aggregated | 1.01758 | 1.00023 | 1.00064 | 4.5904961e-04 |
| WLS | exact_aggregated | 1.82270 | 1.92399 | 1.96373 | 3.6160693e-06 |
| WLS | fit_gradient_aggregated | 1.94395 | 1.96653 | 1.98070 | 6.1325053e-06 |
| WLS | exact_segmented | 0.99439 | 0.97916 | 0.98506 | 6.7869259e-05 |
| WLS | native_segmented | 1.01894 | 1.00193 | 1.00134 | 5.1591907e-04 |

#### Mappa liscia fissa

| Metodo | Esperimento | p 16→32 | p 32→64 | p 64→128 | L2 N=128 |
| --- | --- | ---: | ---: | ---: | ---: |
| LS | native_aggregated | 2.50978 | 2.30452 | 2.13266 | 3.6995724e-06 |
| LS | exact_aggregated | 1.68420 | 1.88434 | 1.95184 | 3.4222781e-06 |
| LS | fit_gradient_aggregated | 2.12047 | 2.05072 | 2.01930 | 3.5519832e-06 |
| LS | exact_segmented | 1.71787 | 1.88396 | 1.94849 | 3.4091541e-06 |
| LS | native_segmented | 2.54494 | 2.33860 | 2.15235 | 3.7455100e-06 |
| WLS | native_aggregated | 1.87180 | 1.95719 | 1.98493 | 3.4876145e-06 |
| WLS | exact_aggregated | 1.68420 | 1.88434 | 1.95184 | 3.4222781e-06 |
| WLS | fit_gradient_aggregated | 2.00150 | 1.99573 | 1.99680 | 3.5110134e-06 |
| WLS | exact_segmented | 1.71787 | 1.88396 | 1.94849 | 3.4091541e-06 |
| WLS | native_segmented | 2.02257 | 1.99826 | 1.99539 | 3.4977192e-06 |

#### Controllo regolare

| Metodo | Esperimento | p 16→32 | p 32→64 | p 64→128 | L2 N=128 |
| --- | --- | ---: | ---: | ---: | ---: |
| LS | native_aggregated | 1.93747 | 1.95815 | 1.97765 | 3.1063191e-06 |
| LS | exact_aggregated | 1.83457 | 1.92344 | 1.96348 | 3.0797176e-06 |
| LS | fit_gradient_aggregated | 2.20157 | 2.05986 | 2.01677 | 3.1726919e-06 |
| LS | exact_segmented | 1.83480 | 1.92345 | 1.96347 | 3.0648826e-06 |
| LS | native_segmented | 1.93880 | 1.95856 | 1.97779 | 3.0916216e-06 |
| WLS | native_aggregated | 1.92241 | 1.95294 | 1.97556 | 3.1024541e-06 |
| WLS | exact_aggregated | 1.83457 | 1.92344 | 1.96348 | 3.0797176e-06 |
| WLS | fit_gradient_aggregated | 2.11053 | 2.02063 | 2.00151 | 3.1473141e-06 |
| WLS | exact_segmented | 1.83480 | 1.92345 | 1.96347 | 3.0648826e-06 |
| WLS | native_segmented | 1.92362 | 1.95330 | 1.97568 | 3.0877423e-06 |

### Regioni distorte, N=128: errore / ordine 64→128

| Metodo | Esperimento | Tutti | Bordo | Ring1 | Deep3 |
| --- | --- | ---: | ---: | ---: | ---: |
| LS | native_aggregated | 1.0861e-03 / 1.018 | 1.5745e-03 / 1.006 | 1.1205e-03 / 1.012 | 1.0775e-03 / 1.005 |
| LS | exact_aggregated | 1.2268e-05 / 1.588 | 9.3866e-05 / 1.038 | 4.5537e-06 / 1.982 | 3.6161e-06 / 1.964 |
| LS | fit_gradient_aggregated | 1.4703e-05 / 1.766 | 1.0220e-04 / 1.176 | 9.6575e-06 / 1.973 | 7.2816e-06 / 1.981 |
| LS | exact_segmented | 7.0693e-05 / 1.040 | 1.6974e-04 / 0.984 | 6.5221e-05 / 0.988 | 6.7869e-05 / 0.985 |
| LS | native_segmented | 1.1438e-03 / 1.018 | 1.7112e-03 / 1.006 | 1.1474e-03 / 1.014 | 1.1343e-03 / 1.006 |
| WLS | native_aggregated | 4.9984e-04 / 1.107 | 1.3594e-03 / 1.005 | 8.1690e-04 / 1.017 | 4.5905e-04 / 1.001 |
| WLS | exact_aggregated | 1.2268e-05 / 1.588 | 9.3866e-05 / 1.038 | 4.5537e-06 / 1.982 | 3.6161e-06 / 1.964 |
| WLS | fit_gradient_aggregated | 1.3842e-05 / 1.718 | 9.9364e-05 / 1.135 | 7.8510e-06 / 1.970 | 6.1325e-06 / 1.981 |
| WLS | exact_segmented | 7.0693e-05 / 1.040 | 1.6974e-04 / 0.984 | 6.5221e-05 / 0.988 | 6.7869e-05 / 0.985 |
| WLS | native_segmented | 5.5708e-04 / 1.098 | 1.5030e-03 / 1.004 | 8.4000e-04 / 1.018 | 5.1592e-04 / 1.001 |

### Errore dei gradienti sull’interno distorto, N=128

| Metodo | Gradiente | L2 | p 64→128 |
| --- | --- | ---: | ---: |
| LS | native_gradient | 6.1808539e-04 | 1.00112 |
| LS | fit_gradient | 1.6659902e-05 | 2.00100 |
| WLS | native_gradient | 3.2023026e-04 | 1.00106 |
| WLS | fit_gradient | 1.2335853e-05 | 2.00076 |

## C. Derivazione dei termini dominanti con MUSCL affine

Sia m=(x_i+x_j)/2, d=x_j-x_i, r=q-m ed epsilon_i=G_i-G_esatto(x_i).
Per ciascuna componente, sviluppando un campo C³ attorno a m:

```text
jump(q) = H(m)[d,r] + O(h³)
          + epsilon_j (q-x_j) - epsilon_i (q-x_i)
jump(m) = O(h³) - 0.5 (epsilon_i+epsilon_j) d
```

Le tracce affini con gradienti esatti hanno lo stesso errore quadratico
al **punto medio primale**. Il salto è quindi O(h³) anche su una stella
asimmetrica: non serve che ciascuna traccia sia esatta al secondo grado.
La norma della fluttuazione è O(|S|*|jump|), dato rho,p lontani da zero.
Con |S|~h e V~h², questo salto dà un errore O(h²) nel residuale normalizzato.
Il volume NODE_POINT con G esatto è esattamente il termine differenziale
puntuale, senza errore di quadratura. È il motivo per cui
exact_aggregated recupera ~2 nell'interno distorto senza SEGMENTED.

Con errore nativo epsilon=O(h), le tracce producono invece salti O(h²)
e il volume introduce A_x epsilon_x+A_y epsilon_y=O(h). Non esiste in
generale una cancellazione sulle stelle jitter, perché epsilon varia
alla scala della griglia. Entrambi i contributi fanno parte del medesimo
residuale completo: non si attribuisce tutto l'errore alla sola norma di G.
Sulla griglia simmetrica, o asintoticamente simmetrica della mappa liscia,
i momenti dispari della LS si annullano oppure sono un ordine più piccoli.
In termini della matrice LS, l'errore principale è

```text
epsilon ~ [sum w d d^T]^-1 sum w d * (0.5 H[d,d]).
```

Per jitter geometrico generale il numeratore ha l'ordine che dà O(h);
con coppie quasi opposte la componente dispari si riduce e si ottiene O(h²).
Il fit quadratico diagnostico elimina quel termine di curvatura usando dati
nodali: l'errore dei due coefficienti G è O(h²), con matrice scalata di rango
pieno su queste famiglie. Salti e volume tornano O(h²) nel residuale.
La tabella conferma questa previsione **sull'operatore completo**, non solo
sul gradiente. Non è una prova per qualunque mesh o campo non liscio.

### Perché segmentare non è una correzione automatica

Sul segmento reale r=q-m è generalmente O(h), non zero. Anche con G esatto,
H[d,r] è O(h²). Il termine principale della fluttuazione di una faccia è

```text
Dminus_s = 0.5 [A(S_s;W(m))-|A(S_s;W(m))|] H(m)[d,r_s] + O(h⁴)
```

Ha ordine h³ integrato e h dopo normalizzazione. Con due punti Gauss si
integra correttamente il termine lineare in r lungo il segmento, ma esso
**non deve annullarsi** su una stella distorta. Sul regolare simmetrico la
somma può cancellarlo; sulla mappa liscia lo sbilanciamento è un ordine più
piccolo. La tabella exact_segmented riproduce precisamente questi casi.

La parte centrata A(S_s) è lineare in S_s; la parte |A(S_s)| non lo è nella
direzione. In generale sum_s |A(S_s)| != |A(sum_s S_s)|. Perciò dividere la
normale aggregata senza ricostruire la direzione geometrica non sarebbe
un confronto valido. Il test usa segmenti reali e normali individuali.
Le matrici valutate su tracce che differiscono da W(m) di O(h) cambiano
A(S) di O(h²); moltiplicando il salto O(h²) si ottengono termini h⁴,
cioè O(h²) dopo V. La variazione di A non spiega il termine dominante O(h)
che permane con gradienti esatti segmentati. La funzione spettrale assoluta
è localmente Lipschitz nel regime liscio ammissibile; non si assume una
Taylor analitica attraverso ogni cambio di segno di un autovalore.

I momenti T_s della quadratura, misurati nel precedente audit, sono non nulli
su mesh jitter. Questo resta vero, ma non significa che il residuale NODE
puntuale debba usare l'integrale geometrico al posto di m: l'attuale operatore
accoppia la cancellazione del salto al punto m con un volume puntuale.
Una variante realmente integrata richiederebbe derivare insieme volume,
rappresentazione e ricostruzione. Cambiare le sole facce non è tale derivazione.
Non si propone un termine empirico per cancellare l'errore dei segmenti.

### Sintesi causale

- NO MUSCL: momento dissipativo non bilanciato, difetto forte O(1) su stelle
  asimmetriche; momento centrato corretto. Consistenza debole scalare distinta.
- MUSCL, jitter: dominano gli errori O(h) di LS/WLS, tramite facce e volume.
  AGGREGATED non impedisce l'ordine 2 interno quando si migliora G.
- Segmenti reali con volume invariato: possono aggiungere un termine O(h)
  anche con G esatto, per H[d,r] e per i pesi direzionali della dissipazione.
- MUSCL, mappa liscia: G nativo e operatore interno sono già ~2. Il limite
  non riguarda indistintamente ogni mesh non cartesiana o non simmetrica.
- Bordo e ring1 restano misurati: non si cerca di correggerli in questo ciclo.

## D. Alternative e raccomandazione prima del Passaggio II

### Problema A

1. **Scelta più semplice per un solver orientato agli shock: conservare il
   percorso NO MUSCL e qualificare il criterio di consistenza.** È scorretto
   chiamarlo puntualmente consistente su ogni stella NODE; è però corretto
   misurare anche la consistenza debole, dimostrata qui per il sottoproblema
   scalare del contatto. Per accettarlo come metodo evolutivo occorrono in
   futuro convergenza della soluzione, stabilità/CFL e positività del sistema.
   Non serve correggere il kernel per il solo difetto 0.036.
2. Se la consistenza puntuale affine è un requisito irrinunciabile, gli
   invarianti attuali (stati costanti per cella, kernel e normali fissati,
   volume nullo) determinano già completamente R. Non resta un parametro
   locale libero che possa annullare il momento dissipativo. Cambiare
   NODE_POINT in MEAN_COEFF non cambia nulla quando G=0; SEGMENTED non basta.
   Occorre una modifica strutturale: ricostruzione/redistribuzione che preservi
   i campi lineari, con un gradiente e un termine interno coerenti. La variante
   già disponibile MUSCL affine è il riferimento semplice in regime liscio;
   non è più il significato attuale di “NO MUSCL”.
3. **Non raccomandati:** togliere la dissipazione, correggere il solo rho,
   invertire segni, o sottrarre il difetto del caso test. La parte centrata
   risolve il contatto affine, ma perde il meccanismo upwind; non è una scelta
   robusta per shock né una correzione accettabile del sistema.

Il test debole sostitutivo è quindi matematicamente giustificato sul contatto,
ma non sostituisce la futura verifica della soluzione Euler nonconservativa.
Nessuna impossibilità di convergenza della soluzione viene dedotta dal difetto
forte; nessuna convergenza del sistema completo viene promessa.

### Problema B

1. **Compromesso più semplice e robusto:** mantenere WLS/MUSCL affine come
   default attuale, documentando ~2 sulle deformazioni lisce e ~1 sul jitter
   geometricamente regolare. Non introdurre SEGMENTED per questo obiettivo.
   Per un solver con shock, il limiter ridurrà comunque localmente l'ordine;
   questa scelta non introduce nuove oscillazioni o stencil larghi.
2. **Se si richiede ~2 anche sul jitter:** la modifica più piccola supportata
   dai dati è migliorare soltanto la ricostruzione di G, rendendola esatta
   per quadratiche su uno stencil locale sufficientemente ricco. Mantenere
   tracce affini, AGGREGATED, volume NODE_POINT e kernel. Il fit diagnostico
   dimostra fattibilità con due anelli, non impone quella specifica soluzione
   produttiva. Un solo anello può essere valutato dove ha rango sufficiente;
   non è automaticamente valido su valenze basse. Servono controlli di rango,
   condizionamento, costo e comportamento con limiter; nessuna Hessiana di
   produzione o ricostruzione quadratica di bordo è necessaria.
3. Una formulazione integrata segmentata coerente richiede una revisione
   congiunta di volume e rappresentazione. È più costosa e non motivata dai
   risultati come prima correzione per l'ordine interno. Resta fuori ambito.

Costi del miglioramento di G: cinque coefficienti intermedi per primitiva
in una piccola LS invece di due, stencil più ampio e (su mesh statiche)
possibilità di precomputare soltanto i pesi dei due gradienti. Si può evitare
di conservare Hessiane, ma non il costo/dipendenza dai campioni aggiuntivi.
Con shock, una regressione polinomiale non limitata può amplificare le
oscillazioni e il condizionamento: il theta attuale garantisce solo le tracce
rho/p nei punti usati, non TVD, stabilità temporale o una soluzione positiva
al passo successivo. Non si raccomanda di sostituire silenziosamente LS/WLS
per tutti i regimi. Il comportamento vicino a discontinuità richiede una
scelta esplicita e validazione nel Passaggio II/futura integrazione temporale.

**Raccomandazione:** nessun cambiamento di kernel, segni o quadratura. Per A,
riconoscere il limite di consistenza forte e mantenere la misura debole come
verifica complementare; per B, mantenere il metodo semplice se il target
primario sono gli shock, oppure autorizzare un miglioramento mirato dei soli
gradienti in regime liscio se il secondo ordine su jitter è prioritario.
Il secondo ordine globale di bordo non è un obiettivo di questo ciclo.

## E. Riproduzione, file ed esiti

Comandi eseguiti dalla radice (il primo driver era quello affine verificato
nella sessione precedente, corrispondente ai sorgenti iniziali puliti):

```sh
/private/tmp/su2-affine-archive-tests/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]' > /private/tmp/su2-irregular-baseline.log 2>&1
./ninja -C build -j 2 SU2_CFD/src/SU2_CFD
/private/tmp/su2-affine-archive-tests/test_driver '[primitive-irregular-audit]' > /private/tmp/su2-irregular-diagnostics.log 2>&1
/private/tmp/su2-affine-archive-tests/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]' > /private/tmp/su2-irregular-results.log 2>&1
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-irregular-cli --cli
python3 UnitTests/export_primitive_euler_irregular.py --log /private/tmp/su2-irregular-results.log
```

Tra riferimento e nuovi test si è ricompilato il solo translation unit
primitive_euler_accuracy_tests.cpp e ricollegato il driver, riusando gli argv
c++ del runner registrati in `/private/tmp/su2-affine-archive-tests-build.log`
con shlex.split/subprocess, log `/private/tmp/su2-irregular-tests-build.log`.
Per ricostruire da zero un nuovo driver con gli stessi sorgenti correnti:

```sh
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/su2-irregular-reproduce
/private/tmp/su2-irregular-reproduce/test_driver '[primitive-geometry],[primitive-euler],[primitive-solver]'
```

Il primo run diagnostico aveva tre casi e 2.311.153 asserzioni, tutti passati;
è precedente all'aggiunta delle norme G e del controllo di somma delle normali.
Il risultato valido finale è la suite completa successiva:
**61 casi, 4.483.704 asserzioni superate**, uscita 0. I 58 casi precedenti sono
conservati senza modificare soglie. Build nativa: uscita 0, “no work to do”,
perché nessun sorgente di produzione è cambiato. Quattro CLI superati con
uscite attese 1/0/1/1 (primitivo bloccato, Euler conservativo riuscito,
harmonic balance e specie rifiutati). Nessun test evolutivo è stato eseguito.
Warning non bloccanti: pragma GCC ignorati da Clang, librerie duplicate,
WINDOW_START_ITER nei vecchi test e bind MPI negato dal sandbox.
Non è stata rieseguita una build ASan o la suite globale SU2.

File modificati/aggiunti in questo Passaggio I:

- UnitTests/SU2_CFD/primitive_euler_accuracy_tests.cpp: tre nuovi casi;
- UnitTests/export_primitive_euler_irregular.py: esportatore indipendente;
- questo rapporto e i collegamenti nei due rapporti di FASE 4 esistenti;
- [CSV MMS](primitive_euler_phase4_irregular_mms.csv): 900 righe dati;
- [CSV gradienti](primitive_euler_phase4_irregular_gradients.csv): 360 righe;
- [CSV split contatto](primitive_euler_phase4_irregular_contact.csv): 12 righe;
- [CSV debole](primitive_euler_phase4_irregular_weak.csv): 15 righe.

I CSV precedenti non sono sovrascritti e rimangono confrontabili. I nuovi CSV
contengono tutti i livelli e tutti i gruppi, non soltanto quelli tabulati qui.
Le tabelle sono arrotondate; nessun livello è escluso dai dati grezzi.

Gli hash SHA-256 dei sette file produttivi iniziali sono confrontati alla
fine del ciclo e invariati. Branch principale e riferimento sperimentale
rimangono ai rispettivi HEAD iniziali; nessun commit, push o sync.
**Il lavoro si ferma al Passaggio I. Qualsiasi modifica del solver richiede
la conferma esplicitamente richiesta dal prompt dell'utente.**
