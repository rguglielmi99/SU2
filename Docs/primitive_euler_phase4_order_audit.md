# FASE 4 — Audit mirato dell'ordine spaziale

Questo rapporto descrive l'operatore affine mantenuto nel branch principale.
L'esperimento quadratico completo e i suoi risultati restano nel branch locale
`experiment/primitive-quadratic-boundary`, commit `d1732e90e3`.
Qui non sono presenti fit numerici quadratici né test che richiedano le loro API.
Le sostituzioni analitiche di G/H sono esclusivamente diagnostiche locali ai test:
servono a isolare l'errore della traccia e non entrano nell'operatore di produzione.

## Diagnosi prima di qualsiasi modifica all'operatore

Il solver di produzione, il kernel di fase 1 e la geometria restano invariati
in questo ciclo. Le varianti con derivate sostituite sono riassemblaggi locali
ai test: non sono opzioni del solver e non costituiscono una nuova modalità
SEGMENTED o DUAL_MEAN.

Sulla famiglia regolare con WLS, N=128 e NODE_POINT:

| Gruppo | L2 normalizzata sul gruppo | Ordine 64→128 |
| --- | ---: | ---: |
| Tutti | 2.052919e-4 | 1.50005 |
| Bordo fisico | 1.315244e-3 | 1.00272 |
| Primo anello | 7.032091e-4 | 1.01149 |
| Interno da tre strati | 3.102454e-6 | 1.97556 |

Quindi una fascia di spessore O(h), con errore puntuale O(h), contribuisce
O(h^3) al quadrato della norma L2 globale: la norma è O(h^(3/2)).
L'interno profondo ha già errore O(h^2). La sola norma normalizzata sul gruppo
non basta a confrontare i contributi: i CSV registrano anche la misura duale
V_gruppo e sqrt(sum_gruppo V*||e||^2 / V_totale).

## Separazione causale mediante riassemblaggi solo diagnostici

Il riassemblaggio di riferimento riproduce prima rhs di produzione usando
G_limited originale (errore massimo richiesto <2e-12). Poi cambia un ingrediente
per volta, mantenendo stati, normali integrate, segni, volume e sorgente.
Tutti i dati analitici provengono dal campo MMS già utilizzato.

A N=128, WLS, mesh regolare, NODE_POINT:

| Esperimento | L2 globale | Ordine 64→128 |
| --- | ---: | ---: |
| Produzione | 2.052919e-4 | 1.50005 |
| G analitico soltanto al bordo | 1.225713e-5 | 1.60621 |
| G analitico ovunque | 1.205010e-5 | 1.57522 |
| G analitico e traccia interna fisica esatta | 3.218358e-6 | 2.02650 |
| G analitico e traccia fisica Taylor quadratica | 3.218355e-6 | 2.02650 |

La sostituzione dei soli gradienti di bordo porta il primo anello a ordine
1.98170. La propagazione dell'errore dei gradienti di bordo attraverso le facce
interne è dunque la causa dominante nel baseline. Rimane un errore al bordo
fisico: una traccia affine differisce dalla soluzione esatta di O(h^2), ma
Dminus ha misura O(h) e V=O(h^2), producendo un contributo O(h) al rhs.
Una traccia quadratica ha difetto O(h^3) e contributo O(h^2).

MEAN_COEFF mantiene un altro termine: con G esatto e traccia quadratica,
L2 globale=6.212259e-6, ordine=1.75566, mentre il bordo resta vicino al primo
ordine. Su un duale troncato delta=x_c-x_i=O(h). Posto B(W,G)=Ax(W)Gx+Ay(W)Gy:

B(W_i+G*delta,G)-B(W_i,G) = D_W B(W_i,G)[G*delta] + O(h^2).

In NODE la sorgente e il riferimento differenziale sono puntuali: questa
differenza O(h) non scompare al bordo. Sulla mesh regolare l'offset è nullo
nell'interno, perciò il difetto è confinato ai duali troncati.

## Esperimento separato

La ricostruzione con fit quadratico a due anelli è stata implementata e validata
in un ciclo successivo, poi archiviata integralmente sul branch sperimentale
su richiesta dell'utente. Non è una modalità del solver principale.
I dettagli e l'evidenza >=2 sono nel rapporto `Docs/primitive_euler_phase4_quadratic.md`
di quel branch. Non si attribuisce quel risultato all'operatore affine.

## Configurazione, regioni e risultati completi

Tutte le misure usano NODE, AGGREGATED, MUSCL attivo, limiter NONE e gamma
1.4 effettivo. La soluzione e la sorgente sono quelle del precedente MMS,
con derivate analitiche esplicite. Si misurano entrambe le routine native LS
(non pesata) e WLS (peso inverso del quadrato della distanza), entrambe le
varianti di volume e N=8,16,32,64,128, sulle stesse famiglie regolare/distorta.

La classe Grid nei test usa la distanza topologica
`l=min(i,j,N-i,N-j)` sulla griglia triangolata: bordo l=0, primo anello l=1,
secondo anello l=2 e interno profondo l>=3. I primi tre strati 0/1/2 sono quindi
esclusi dall'interno profondo. Sulla famiglia distorta la connettività è la
stessa e gli spostamenti sono limitati da costanti volte h: i gruppi conservano
il significato geometrico di fasce di spessore O(h). I vertici d'angolo sono
registrati separatamente come sottoinsieme del bordo.

Per un gruppo A:

```text
V_A = sum_(i in A) V_i
E_A = sqrt(sum_(i in A) V_i * ||rhs_i||^2 / V_A)
C_A = sqrt(sum_(i in A) V_i * ||rhs_i||^2 / V_totale)
p_A = log2(E_A(N)/E_A(2N))
```

C_A è il contributo alla norma globale, mentre E_A misura l'errore locale nella
regione. Bordo, ring1, ring2 e deep3 formano una partizione; corners è soltanto
una diagnostica aggiuntiva e non va sommato una seconda volta.
Il CSV riporta ordine 0 sulla prima mesh o quando il rapporto non è definito.

### Residuale nativo, WLS, MEAN_COEFF (default)

| N | E tutti | E bordo | E primo anello | E interno profondo |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 1.3175672e-2 | 2.1864691e-2 | 1.2784818e-2 | 7.1021424e-4 |
| 16 | 4.6471079e-3 | 1.0676648e-2 | 5.9585450e-3 | 1.7906209e-4 |
| 32 | 1.6412153e-3 | 5.2850949e-3 | 2.8816128e-3 | 4.7239024e-5 |
| 64 | 5.8000886e-4 | 2.6307276e-3 | 1.4176655e-3 | 1.2201332e-5 |
| 128 | 2.0502873e-4 | 1.3126045e-3 | 7.0320905e-4 | 3.1024541e-6 |
| p 64→128 | 1.50025 | 1.00303 | 1.01149 | 1.97556 |

Sulla famiglia distorta:

| N | E tutti | E bordo | E primo anello | E interno profondo |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 1.4006077e-2 | 2.1692527e-2 | 1.3231394e-2 | 7.1051602e-3 |
| 16 | 5.8715373e-3 | 1.1104102e-2 | 7.0146644e-3 | 3.7292169e-3 |
| 32 | 2.4289367e-3 | 5.4663897e-3 | 3.3156352e-3 | 1.8410063e-3 |
| 64 | 1.0775725e-3 | 2.7215593e-3 | 1.6535726e-3 | 9.2007032e-4 |
| 128 | 5.0029446e-4 | 1.3560464e-3 | 8.1690243e-4 | 4.5975054e-4 |
| p 64→128 | 1.10694 | 1.00503 | 1.01735 | 1.00089 |

L'interpretazione «solo fascia di bordo» è verificata per la famiglia regolare,
non per quella distorta: nella seconda anche l'interno profondo è circa O(h).
I CSV contengono anche LS e NODE_POINT, senza selezionare solo la combinazione
più favorevole. A N=128, regolare WLS/NODE_POINT, bordo più primo anello
rappresentano il **99.9777% del quadrato dell'errore globale**.

### Contributi separati e cancellazioni

Si decompone il rhs in F=-faces/V, T=-boundary/V ed E_V=S-volume/V.
I valori sotto sono norme E_A di questi vettori, non somme di moduli per nodo.
WLS/NODE_POINT, mesh regolare, N=128:

| Termine | Bordo | Primo anello | Interno profondo |
| --- | ---: | ---: | ---: |
| F | 1.0869391e-3 | 7.0766224e-4 | 8.2943724e-6 |
| T | 9.6250324e-5 | 0 | 0 |
| E_V | 1.7701348e-3 | 7.0055004e-6 | 5.2409460e-6 |
| rhs=F+T+E_V | 1.3152443e-3 | 7.0320905e-4 | 3.1024541e-6 |

Le norme dei contributi non si sommano scalarmente: nell'interno vi sono
cancellazioni. Nel primo anello domina F, pur avendo gradienti locali già
O(h^2): le tracce sulla faccia verso un nodo di bordo dipendono dal suo
G_bordo, accurato soltanto O(h). Sostituendo solo quel gradiente, senza
cambiare fluttuazioni o normali, l'ordine del primo anello torna vicino a 2.
Non si tratta quindi di un errore dei segni o del kernel path-conservative.

## Test dei gradienti indipendente dall'operatore

Il secondo nuovo test costruisce quattro campi quadratici positivi
`W_k = a_k + b_k*x + c_k*y + 0.5*Hxx_k*x^2 + Hxy_k*x*y + 0.5*Hyy_k*y^2`.
Le Hessiane sono costanti, diverse per componente e con termini misti non nulli.
Chiama direttamente `computeGradientsLeastSquares` nativo con solver=nullptr,
confrontando con `(b_k+Hxx_k*x+Hxy_k*y,c_k+Hxy_k*x+Hyy_k*y)`.
Non chiama l'assemblaggio del residuale, non usa limiter e non modifica la
soluzione per migliorare il confronto.

A N=128:

| Famiglia / metodo | E_G tutti | E_G bordo | E_G primo anello | E_G profondo |
| --- | ---: | ---: | ---: | ---: |
| Regolare LS | 2.1217081e-4 | 1.7006914e-3 | 2.1758e-14 | 2.2448e-14 |
| Regolare WLS | 1.9690319e-4 | 1.5783112e-3 | 2.2157e-14 | 2.2676e-14 |
| Distorta LS | 1.0867784e-3 | 1.7635342e-3 | 9.5634481e-4 | 1.0760997e-3 |
| Distorta WLS | 5.8442496e-4 | 1.6100944e-3 | 4.9179528e-4 | 5.5489426e-4 |

Gli ordini di bordo sono circa 1 in tutti i casi. Sulla famiglia regolare i
vicini interni sono accoppiati in direzioni opposte: i termini quadratici
nel sistema LS si cancellano e il gradiente quadratico è riprodotto fino
all'arrotondamento. Gli ordini negativi del CSV per errori ~1e-14 descrivono
il rumore floating-point, non una perdita di consistenza. Sulla famiglia
distorta, invece, l'errore profondo è reale e circa di primo ordine.

Sul campo MMS non quadratico, WLS regolare ha a N=128 errore del gradiente
9.100736e-4 al bordo (p=1.00108), 4.300444e-6 nel primo anello (p=2.00046)
e 4.331865e-6 nell'interno profondo (p=2.00074). Questo separa l'errore locale
della ricostruzione dalla propagazione attraverso le facce nel residuale.

## Momenti della quadratura AGGREGATED

Per una faccia duale effettiva composta dai segmenti m, con centro Q_m e
normale integrata S_m, il test confronta:

```text
I_d[f] = sum_m integral_segment_m f(x) * n_(m,d) dl
Iagg_d[f] = f(M) * S_d,  S=sum_m S_m, M=(x_i+x_j)/2
T_(d,a) = sum_m S_(m,d) * (Q_(m,a)-M_a)
```

Espandendo f attorno a M:

```text
I_d[f]-Iagg_d[f] = sum_a f_,a(M)*T_(d,a)
                  + 0.5*sum_ab f_,ab(M)*U_(d,ab) + O(h^4)
U_(d,ab) = sum_m S_(m,d) * average_segment[(x_a-M_a)*(x_b-M_b)]
```

L'approssimazione aggregata riproduce tutte le funzioni affini pesate con la
normale se e solo se T=0. I momenti vengono integrati esattamente con due
punti di Gauss per segmento (polinomi fino al secondo grado); si confrontano
anche le somme delle normali con le metriche native. Non viene aggiunta al
solver una quadratura segmentata.

Per regolari triangoli con diagonale fissa, sulle facce con due triangoli
incidenti i segmenti sono simmetrici attorno a M: T=0 entro arrotondamento.
Sulle facce con un solo triangolo (edge primale sul bordo fisico) tale simmetria
manca. Sulla mesh distorta manca genericamente anche per due triangoli.

Valori medi a N=128, norma euclidea dei componenti indipendenti dei tensori:

| Famiglia | Triangoli incidenti | ||T|| medio / h^2 | ||U|| medio / h^3 |
| --- | ---: | ---: | ---: |
| Regolare | 1 | 0.0694444 | 0.0158132 |
| Regolare | 2 | 6.44e-15 | 0.0235623 |
| Distorta | 1 | 0.0710715 | 0.0166682 |
| Distorta | 2 | 0.0320295 | 0.0255259 |

La colonna U usa le componenti xx,xy,yy senza raddoppiare xy, come nel test.
Non è un'affermazione di esattezza quadratica: AGGREGATED non la possiede.
I momenti sono una verifica geometrica di f*n, non una prova automatica
sull'integrale delle fluttuazioni, che comprende anche |A_n| non lineare nella
direzione. In particolare non si può dedurre l'ordine dell'intero operatore
dai soli momenti. La presenza di T non nullo identifica una limitazione di
quadratura, ma il prototipo archiviato sul branch sperimentale raggiunge ordine 2 regolare mantenendo
AGGREGATED: introdurre SEGMENTED non è necessario per risolvere il difetto
prioritario, né correggerebbe da solo i gradienti o il volume al bordo.

## Verifiche mantenute e riproduzione

I tre test `[primitive-order-audit]` mantengono senza ridurre le risoluzioni:
MMS regionale e separazione causale, accuratezza LS/WLS su campi quadratici,
momenti geometrici delle facce duali. Le famiglie regolari e distorte sono
misurate fino a N=128. Restano i sei gruppi all/boundary/ring1/deep3/ring2/corners.
Il riassemblaggio test-only deve riprodurre rhs del solver entro 2e-12 prima
di sostituire derivate analitiche. I gradienti quadratici negli stencil profondi
simmetrici devono restare esatti entro 1e-10; le normali dei segmenti devono
sommare a quelle SU2 entro 2e-14. Nessun test d'ordine affine viene trasformato
in una pretesa di secondo ordine.

I test del fit, del candidato e della MMS di produzione quadratica sono nel
branch sperimentale. Il test del contatto affine MUSCL NO rimane integralmente
in primitive_euler_residual_tests.cpp. Non sono introdotte correzioni numeriche.

Comandi ed esiti della separazione: [rapporto Git e verifiche](primitive_euler_branch_archive.md).
CSV regionale e figura sono rigenerati dai soli test conservati; il CSV dei
momenti resta disponibile. I risultati del quadratico sono recuperabili dal
commit sperimentale, non mescolati ai nuovi artefatti del principale.
